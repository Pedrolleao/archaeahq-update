#!/usr/bin/env python3
"""
archaeahq_update.py — find archaeal genomes newly deposited at NCBI, run the ArchaeaHQ curation
pipeline on them and recommend which ones should be added to the database.

    archaeahq_update.py setup                      create/verify the tools and the CheckM2 database
    archaeahq_update.py fetch-db                   download ArchaeaHQ v1.0 from figshare → Archaea_HQ-v1.0/
    archaeahq_update.py check                      list new NCBI assemblies per kingdom (no download)
    archaeahq_update.py run                        full evaluation → two result tables
    archaeahq_update.py sketch-db --out DIR        build a reusable skani sketch of the database
    archaeahq_update.py release                    apply the Add/Replace rows → next database version

Database versions live in folders named Archaea_HQ-v<major>.<minor>/ next to the script (or in
--releases-dir), each with fna/, ArchaeaHQ-Info.tsv and release.json. All commands use the newest
version automatically; `release` builds the next one and removes the superseded version once the
new folder has been verified.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "lib"))

import common as C            # noqa: E402
import ui                     # noqa: E402
import envcheck               # noqa: E402
import ncbi                   # noqa: E402
import quality as Q           # noqa: E402
import redundancy as R        # noqa: E402
import rna                    # noqa: E402
import report                 # noqa: E402
import fetchdb                # noqa: E402
from classify_environments import classify  # noqa: E402

STAGE_TITLES = OrderedDict([
    ("env", "Checking tools and databases"),
    ("list", "Listing NCBI assemblies and finding new genomes"),
    ("download", "Downloading new genomes"),
    ("metadata", "Fetching assembly metadata and taxonomy"),
    ("quality", "Estimating completeness and contamination (CheckM2)"),
    ("redundancy", "Comparing with the database and within the batch (skani)"),
    ("rna", "Counting rRNA and tRNA genes (barrnap)"),
    ("report", "Writing recommendations"),
])

LIST_COLUMNS = ["accession", "kingdom", "asm_name", "paired", "status", "organism", "taxid", "release_date",
                "level", "bioproject", "biosample", "twin_of_db", "twin_type", "download"]


# ═══════════════════════════════════════════════════════════════════════════════
# CLI
# ═══════════════════════════════════════════════════════════════════════════════

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="archaeahq_update.py", description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--version", action="version", version=f"archaeahq-update {C.VERSION}")
    sub = p.add_subparsers(dest="command", metavar="COMMAND")

    def common_opts(sp, need_db=False):
        sp.add_argument("--workdir", default="./archaeahq_update_work",
                        help="work directory for runs, caches and the ledger (default: ./archaeahq_update_work)")
        sp.add_argument("--threads", type=int, default=max(1, (os.cpu_count() or 4) - 2),
                        help="CPU threads for CheckM2 / skani / barrnap")
        sp.add_argument("--api-key", default=None, help="NCBI API key (or set $NCBI_API_KEY)")
        sp.add_argument("--conda-env", default=envcheck.DEFAULT_ENV_NAME,
                        help="name of the conda environment to create when tools are missing")
        sp.add_argument("--checkm2-db", default=None, help="folder (or .dmnd file) of the CheckM2 database")
        sp.add_argument("--no-install", action="store_true", help="never create a conda environment; fail instead")
        sp.add_argument("--plain", action="store_true", help="plain text output (no colours, no progress bars)")
        sp.add_argument("--releases-dir", default=str(HERE),
                        help=f"folder holding the database versions {C.RELEASE_PREFIX}<major>.<minor>/ (default: next to the script)")
        sp.add_argument("--db-table", default=None,
                        help="ArchaeaHQ information table (default: the newest release in --releases-dir, "
                             "else the bundled v1.0 data/ArchaeaHQ-Info.tsv)")
        sp.add_argument("--kingdoms", default="eury,tack,dpann,asgard",
                        help="comma-separated subset of: eury,tack,dpann,asgard")
        if need_db:
            g = sp.add_mutually_exclusive_group(required=False)
            g.add_argument("--db-fna", default=None,
                           help="folder with the database genome FASTA files (default: fna/ of the newest release)")
            g.add_argument("--db-sketch", default=None, help="folder produced by `skani sketch` / `sketch-db`")

    s = sub.add_parser("setup", help="create/verify the conda environment and the CheckM2 database")
    common_opts(s)

    s = sub.add_parser("fetch-db", help="download ArchaeaHQ v1.0 from figshare and set it up as the first database version")
    common_opts(s)
    s.add_argument("--from-zip", default=None, help="use an already downloaded fna.zip instead of downloading it")
    s.add_argument("--extras", default="",
                   help="comma-separated extras to download as well: faa (proteins, 6.3 GB), 16s, tables")
    s.add_argument("--keep-zip", action="store_true", help="keep the zip archive(s) after unpacking")
    s.add_argument("--no-verify", action="store_true", help="skip the MD5 check of the downloaded files")
    s.add_argument("--yes", "-y", action="store_true", help="do not ask for confirmation")

    s = sub.add_parser("check", help="list new NCBI assemblies per kingdom (nothing is downloaded)")
    common_opts(s)
    s.add_argument("--reevaluate", action="store_true", help="ignore the ledger of previously evaluated genomes")
    s.add_argument("--accessions-file", default=None, help="evaluate this accession list instead of the NCBI diff")

    s = sub.add_parser("run", help="full evaluation of the new genomes")
    common_opts(s, need_db=True)
    s.add_argument("--max-genomes", type=int, default=None, help="cap the number of new genomes per kingdom (newest first)")
    s.add_argument("--accessions-file", default=None, help="evaluate this accession list instead of the NCBI diff")
    s.add_argument("--reevaluate", action="store_true", help="ignore the ledger of previously evaluated genomes")
    s.add_argument("--include-twins", action="store_true",
                   help="download and fully process GCA/GCF twins of database genomes (default: report them without downloading)")
    s.add_argument("--from-stage", choices=list(STAGE_TITLES), default=None, help="recompute from this stage on")
    s.add_argument("--run-name", default=None, help="name of the run folder (default: run_<date>)")
    s.add_argument("--no-rna", action="store_true", help="skip rRNA/tRNA counting (columns become NA)")
    s.add_argument("--sankey", action="store_true", help="also render the Sankey figure (needs plotly, kaleido)")
    s.add_argument("--lowmem", action="store_true", help="pass --lowmem to CheckM2")

    s = sub.add_parser("sketch-db", help="build a reusable skani sketch of the database genomes")
    common_opts(s)
    s.add_argument("--db-fna", default=None, help="folder with the database FASTA files (default: fna/ of the newest release)")
    s.add_argument("--out", required=True, help="output folder for the sketch")

    s = sub.add_parser("release", help="apply the Add/Replace recommendations of a run: build the next database version")
    common_opts(s)
    s.add_argument("--db-fna", default=None, help="FASTA folder of the current version (default: fna/ of the newest release)")
    s.add_argument("--out", default=None,
                   help=f"folder of the new version (default: <releases-dir>/{C.RELEASE_PREFIX}<next version>)")
    s.add_argument("--keep-previous", action="store_true",
                   help="keep the superseded version (default: its folder is removed once the new one is verified)")
    s.add_argument("--run", default=None, help="run whose recommendations are applied (default: newest run in the work directory)")
    s.add_argument("--table", default=None, help="full result table to apply instead of the one of --run")
    s.add_argument("--new-fna", default=None, help="folder with the FASTA of the new genomes (default: <run>/02_download/fna)")
    s.add_argument("--accept", default=None,
                   help="file with one accession per line: apply only these Add/Replace rows (default: all)")
    s.add_argument("--version", dest="release_version", default=None,
                   help="label of the new version, e.g. v1.2 (default: current version + 0.1)")
    s.add_argument("--copy", action="store_true", help="copy the FASTA files (default: hard links, copy if not possible)")
    s.add_argument("--dry-run", action="store_true", help="select and compare, but write and delete nothing")
    return p


# ═══════════════════════════════════════════════════════════════════════════════
# Run context
# ═══════════════════════════════════════════════════════════════════════════════

class Context:
    def __init__(self, args):
        self.args = args
        self.workdir = Path(args.workdir).expanduser().resolve()
        self.cache_dir = self.workdir / "cache"
        self.runs_dir = self.workdir / "runs"
        self.ledger = self.workdir / "ledger.tsv"
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.cache_dir.mkdir(exist_ok=True)
        self.api_key = args.api_key or os.environ.get("NCBI_API_KEY")
        self.th = C.load_thresholds()
        self.kingdoms = [k for k in C.load_kingdoms() if k["short"] in {x.strip().lower() for x in args.kingdoms.split(",")}]
        if not self.kingdoms:
            raise SystemExit("--kingdoms selected none of eury,tack,dpann,asgard")
        self.kingdom_labels = [k["label"] for k in self.kingdoms]
        self.taxid_to_label = {str(k["taxid"]): k["label"] for k in C.load_kingdoms()}
        self.releases_dir = Path(args.releases_dir).expanduser().resolve()
        self.current_release: Optional[Path] = C.newest_release(self.releases_dir)
        if args.db_table:
            self.db_table = Path(args.db_table).expanduser().resolve()
        elif self.current_release:
            self.db_table = self.current_release / C.RELEASE_TABLE
        else:
            self.db_table = C.DB_TABLE_PATH
        self.db_version = C.release_version_of_table(self.db_table)
        self.db_rows = C.read_tsv(self.db_table)
        self.db_acc = {r["Assembly ID"] for r in self.db_rows}
        self.db_numeric: Dict[str, str] = {}
        for r in self.db_rows:
            self.db_numeric.setdefault(C.acc_numeric(r["Assembly ID"]), r["Assembly ID"])
        self.db_by_acc = {r["Assembly ID"]: r for r in self.db_rows}
        self.db_index = {r["Name"]: {"acc": r["Assembly ID"], "park": C.to_float(r.get("Park Quality_Score"))}
                         for r in self.db_rows}
        self.tools: Dict[str, envcheck.Tool] = {}
        self.checkm2_db: Optional[Path] = None
        self.run_dir: Optional[Path] = None
        self.log_path: Optional[Path] = None
        self.run_name = ""

    def default_db_fna(self) -> Optional[Path]:
        """fna/ of the release the current table belongs to, when it exists."""
        d = self.db_table.parent / "fna"
        return d if d.is_dir() and C.is_release_folder(self.db_table.parent) else None

    def describe_db(self) -> str:
        return f"{self.db_version}: {len(self.db_rows):,} genomes ({self.db_table})"

    def stage_dir(self, stage: str) -> Path:
        idx = list(STAGE_TITLES).index(stage)
        d = self.run_dir / f"{idx:02d}_{stage}"
        d.mkdir(parents=True, exist_ok=True)
        return d

    def known_accessions(self, reevaluate: bool) -> Dict[str, str]:
        """accession → why it is not 'new'."""
        known = {a: "in database" for a in self.db_acc}
        if not reevaluate:
            if C.EVALUATED_PATH.exists():
                for r in C.read_tsv(C.EVALUATED_PATH):
                    known.setdefault(r["Assembly ID"], "evaluated for v1.0")
            if self.ledger.exists():
                for r in C.read_tsv(self.ledger):
                    known.setdefault(r["Assembly ID"], f"evaluated in {r.get('Run', 'a previous run')}")
        return known


# ═══════════════════════════════════════════════════════════════════════════════
# Stages
# ═══════════════════════════════════════════════════════════════════════════════

def stage_env(ctx: Context, required: List[str], need_checkm2_db: bool) -> None:
    a = ctx.args
    ctx.tools = envcheck.resolve_tools(required, envcheck.OPTIONAL_TOOLS, auto_install=not a.no_install,
                                       env_name=a.conda_env, log_path=ctx.log_path, ui=ui)
    rows = [[t.name, t.version or "?", t.source] for t in ctx.tools.values()]
    ui.table("Tools", ["tool", "version", "where"], rows, justify={"version": "left", "where": "left"})
    if need_checkm2_db and "checkm2" in ctx.tools:
        ctx.checkm2_db = envcheck.ensure_checkm2_db(ctx.tools["checkm2"], ctx.cache_dir / "CheckM2_database",
                                                    a.checkm2_db, ctx.log_path, ui=ui)
        ui.ok(f"CheckM2 database: {ctx.checkm2_db}")
    if ctx.run_dir:
        C.save_json(ctx.run_dir / "environment.json", envcheck.environment_summary(ctx.tools, ctx.checkm2_db))


def stage_list(ctx: Context, out_dir: Path) -> List[dict]:
    a = ctx.args
    out_tsv = out_dir / "new_accessions.tsv"
    rows: List[dict] = []
    if getattr(a, "accessions_file", None):
        accs = [ln.strip() for ln in Path(a.accessions_file).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
        ui.info(f"evaluating {len(accs)} accession(s) from {a.accessions_file}")
        with ui.spinner("querying NCBI for the assembly summaries"):
            summ = ncbi.summarize_accessions(ctx.tools, accs, out_dir / "ncbi_accessions.tsv", ctx.api_key, ctx.log_path)
        by_acc = {r["accession"]: r for r in summ}
        for acc in accs:
            r = by_acc.get(acc, {"accession": acc})
            rows.append(_listing_row(ctx, r, "unknown"))
    else:
        known = ctx.known_accessions(a.reevaluate)
        per_k = OrderedDict()
        for k in ctx.kingdoms:
            with ui.spinner(f"listing {k['label']} assemblies at NCBI (taxid {k['taxid']})"):
                listing = ncbi.list_taxon(ctx.tools, k["taxid"], out_dir / f"ncbi_{k['short']}.tsv", ctx.api_key, ctx.log_path)
            n_db = sum(1 for r in listing if known.get(r["accession"]) == "in database")
            n_eval = sum(1 for r in listing if r["accession"] in known and known[r["accession"]] != "in database")
            fresh = [r for r in listing if r["accession"] not in known]
            fresh.sort(key=lambda r: (r.get("release_date", ""), r["accession"]), reverse=True)
            if getattr(a, "max_genomes", None):
                fresh = fresh[: a.max_genomes]
            for r in fresh:
                rows.append(_listing_row(ctx, r, k["label"]))
            mine = [r for r in rows if r["kingdom"] == k["label"]]
            per_k[k["label"]] = [len(listing), n_db, n_eval, len(fresh),
                                 sum(1 for r in mine if r["twin_type"] == "gca_gcf_twin"),
                                 sum(1 for r in mine if r["twin_type"] == "newer_version")]
        tot = [sum(v[i] for v in per_k.values()) for i in range(6)]
        ui.kingdom_table("New assemblies at NCBI",
                         ["at NCBI", "in DB", "evaluated before", "new", "GCA/GCF twins", "newer versions"],
                         per_k, total_row=tot)
        if tot[5]:
            ui.info(f"{tot[5]} database genome(s) have a newer assembly version at NCBI → evaluated for 'Replace'")
    C.write_tsv(out_tsv, rows, LIST_COLUMNS)
    return rows


def _listing_row(ctx: Context, r: dict, kingdom: str) -> dict:
    acc = r.get("accession", "")
    twin = ctx.db_numeric.get(C.acc_numeric(acc), "") if acc not in ctx.db_acc else ""
    twin_type = twin_kind(acc, twin) if twin else ""
    include_twins = getattr(ctx.args, "include_twins", False)
    if not twin or include_twins or twin_type == "newer_version":
        download = "yes"
    else:
        download = "no (GCA/GCF twin of database genome)"
    if r.get("status") and r["status"] != "current":
        download = f"no (assembly status: {r['status']})"
    return {"accession": acc, "kingdom": kingdom, "asm_name": r.get("asm_name", ""), "paired": r.get("paired", ""),
            "status": r.get("status", ""), "organism": r.get("organism", ""), "taxid": r.get("taxid", ""),
            "release_date": r.get("release_date", ""), "level": r.get("level", ""),
            "bioproject": r.get("bioproject", ""), "biosample": r.get("biosample", ""),
            "twin_of_db": twin, "twin_type": twin_type, "download": download}


def twin_kind(acc: str, db_acc: str) -> str:
    """'newer_version' when acc is a later version of the same GCA/GCF record, else 'gca_gcf_twin'."""
    a, d = C.acc_parts(acc), C.acc_parts(db_acc)
    if a and d and a[0] == d[0] and a[2] > d[2]:
        return "newer_version"
    return "gca_gcf_twin"


def stage_download(ctx: Context, listing: List[dict], out_dir: Path) -> Dict[str, Optional[Path]]:
    todo = [r["accession"] for r in listing if r["download"] == "yes"]
    skipped = len(listing) - len(todo)
    if skipped:
        ui.info(f"{skipped} genome(s) not downloaded (GCA/GCF twins of database genomes or non-current assemblies)")
    result: Dict[str, Optional[Path]] = {}
    if todo:
        th = ctx.th["download"]
        with ui.progress() as prog:
            task = prog.add_task("downloading genomes", total=len(todo))
            result = ncbi.download_genomes(ctx.tools, todo, out_dir / "fna", ctx.api_key, ctx.log_path,
                                           batch_size=th["batch_size"], retries=th["retries"],
                                           progress_cb=lambda n: prog.update(task, advance=n))
    ok_n = sum(1 for p in result.values() if p)
    failed = [a for a, p in result.items() if not p]
    C.write_tsv(out_dir / "download_status.tsv",
                ({"accession": a, "name": (p.stem if p else ""), "path": (str(p) if p else ""),
                  "status": "ok" if p else "failed"} for a, p in result.items()),
                ["accession", "name", "path", "status"])
    if failed:
        ui.warn(f"{len(failed)} genome(s) could not be downloaded (see download_status.tsv)")
    ui.ok(f"{ok_n} genome FASTA files in {out_dir / 'fna'}")
    return result


def stage_metadata(ctx: Context, accessions: List[str], out_dir: Path) -> Dict[str, dict]:
    import pypipeline.runner as runner
    tax_dir = ctx.cache_dir / "taxonkit"
    asm_dir = ctx.cache_dir / "assembly_meta"
    asm_out, tax_out, comb_out = out_dir / "assembly_metadata.tsv", out_dir / "taxonomy_lineage.tsv", out_dir / "combined.tsv"
    with ui.spinner(f"fetching metadata + taxonomy for {len(accessions)} accession(s) "
                    f"(NCBI Datasets API; taxonomy dump downloaded on first use)"):
        rc = runner.run(accessions=accessions, assembly_data_dir=str(asm_dir), taxonomy_data_dir=str(tax_dir),
                        api_key=ctx.api_key, rebuild=False, verbose=False, show_progress=False,
                        assembly_out_path=str(asm_out), taxonomy_out_path=str(tax_out), combined_out_path=str(comb_out))
    if rc != 0:
        raise RuntimeError("metadata pipeline returned a non-zero status")
    meta = {r["assembly_accession"]: r for r in C.read_tsv(comb_out)}
    missing = [a for a in accessions if not meta.get(a, {}).get("organism_name")]
    if missing:
        ui.warn(f"{len(missing)} accession(s) returned no metadata from NCBI")
    ui.ok(f"metadata for {len(meta) - len(missing)} genome(s) → {comb_out.name}")
    return meta


def stage_quality(ctx: Context, fnas: List[Path], fna_dir: Path, out_dir: Path) -> Dict[str, dict]:
    if not fnas:
        ui.warn("no genomes to assess")
        return {}
    ui.info(f"CheckM2 on {len(fnas)} genome(s) with {ctx.args.threads} threads — this is the slowest step "
            f"(roughly 1-3 min per 100 genomes once the models are loaded)")
    with ui.spinner("running checkm2 predict"):
        report_path = Q.run_checkm2(ctx.tools, fna_dir, out_dir / "checkm2", ctx.args.threads, ctx.checkm2_db,
                                    ctx.log_path, lowmem=ctx.args.lowmem)
    qual = Q.parse_quality_report(report_path)
    shutil.copy(report_path, out_dir / "quality_report.tsv")
    n_pass = sum(1 for v in qual.values()
                 if Q.passes_completeness(v["Completeness"], ctx.th) and Q.passes_contamination(v["Contamination"], ctx.th))
    ui.ok(f"{n_pass}/{len(qual)} genome(s) pass completeness ≥{ctx.th['quality']['min_completeness']:g}% "
          f"and contamination ≤{ctx.th['quality']['max_contamination']:g}%")
    return qual


def stage_redundancy(ctx: Context, fnas: List[Path], names_acc: Dict[str, str], qual: Dict[str, dict],
                     out_dir: Path) -> Dict[str, dict]:
    a = ctx.args
    db_fna = Path(a.db_fna).expanduser() if a.db_fna else None
    db_sketch = Path(a.db_sketch).expanduser() if a.db_sketch else None
    with ui.spinner("preparing the skani sketch of the database"):
        sketch_dir, n_db = R.ensure_db_sketch(ctx.tools, db_fna, db_sketch, ctx.cache_dir, a.threads, ctx.log_path, ui=ui)
    ui.ok(f"database sketch: {sketch_dir} ({n_db if n_db >= 0 else '?'} genomes)")
    if not fnas:
        ui.warn("no downloaded genomes to compare")
        return {}
    with ui.spinner(f"skani search: {len(fnas)} new genome(s) against the database"):
        R.search_db(ctx.tools, sketch_dir, fnas, out_dir / "new_vs_db.tsv", a.threads, ctx.log_path)
    with ui.spinner("skani triangle: new genomes against each other"):
        R.triangle(ctx.tools, fnas, out_dir / "new_vs_new.tsv", a.threads, ctx.log_path)
    hits_db, hits_new = R.parse_skani(out_dir / "new_vs_db.tsv"), R.parse_skani(out_dir / "new_vs_new.tsv")
    new = {}
    for f in fnas:
        name = C.strip_fna(f.name)
        q = qual.get(name, {})
        eligible = bool(Q.passes_completeness(q.get("Completeness"), ctx.th)) and bool(Q.passes_contamination(q.get("Contamination"), ctx.th))
        new[name] = {"acc": names_acc.get(name, C.name_to_acc(name)), "park": q.get("_park"), "eligible": eligible}
    res = R.evaluate(new, ctx.db_index, hits_db, hits_new, ctx.th)
    C.write_tsv(out_dir / "redundancy.tsv",
                ({"Name": n, "identical_to_db": r["identical_db"] or "", "identical_evidence": r["identical_db_ani"] or "",
                  "same_species_as_db": r["species_db"] or "", "species_evidence": r["species_db_ani"] or "",
                  "note": r["species_db_note"], "batch_representative": "" if r["batch_rep"] is None else ("yes" if r["batch_rep"] else "no"),
                  "batch_partner": r["batch_partner"] or "", "batch_cluster": "; ".join(r["batch_cluster"])} for n, r in res.items()),
                ["Name", "identical_to_db", "identical_evidence", "same_species_as_db", "species_evidence", "note",
                 "batch_representative", "batch_partner", "batch_cluster"])
    n_id = sum(1 for r in res.values() if r["identical_db"])
    n_sp = sum(1 for r in res.values() if r["species_db"] and not r["identical_db"])
    n_batch = sum(1 for r in res.values() if r["batch_rep"] is False)
    ui.ok(f"identical to a database genome: {n_id} · same species as a database genome: {n_sp} · "
          f"redundant within the batch: {n_batch}")
    return res


def stage_rna(ctx: Context, fnas: List[Path], out_dir: Path) -> Dict[str, tuple]:
    if not fnas:
        return {}
    new_b = envcheck.barrnap_has_trna(ctx.tools["barrnap"])
    if not new_b:
        ui.warn("old barrnap (<1.0) found: rRNA from barrnap, tRNA from aragorn")
    with ui.progress() as prog:
        task = prog.add_task("barrnap (rRNA + tRNA)", total=len(fnas))
        counts = rna.run_all(ctx.tools, fnas, out_dir / "gff", ctx.args.threads, new_b, ctx.log_path,
                             progress_cb=lambda n: prog.update(task, advance=n))
    C.write_tsv(out_dir / "rna_counts.tsv", ({"Name": n, "rRNA": r, "tRNA": t} for n, (r, t) in sorted(counts.items())),
                ["Name", "rRNA", "tRNA"])
    return counts


# ═══════════════════════════════════════════════════════════════════════════════
# Report assembly
# ═══════════════════════════════════════════════════════════════════════════════

def _text(v, default=C.MISSING_TEXT) -> str:
    v = (v or "").strip() if isinstance(v, str) else v
    return str(v) if v not in (None, "") else default


def _kingdom_from_lineage(ctx: Context, meta_row: dict, fallback: str, acc: str = "") -> str:
    ids = set((meta_row or {}).get("lineage_taxids", "").split(";"))
    for taxid, label in ctx.taxid_to_label.items():
        if taxid in ids:
            return label
    if fallback and fallback != "unknown":
        return fallback
    # Last resort (e.g. taxid deleted from NCBI Taxonomy): look the accession up in any kingdom
    # listing cached in the work directory by earlier `check` / `run` invocations.
    return _cached_listing_kingdoms(ctx).get(acc, "Other archaea")


def _cached_listing_kingdoms(ctx: Context) -> Dict[str, str]:
    if getattr(ctx, "_listing_kingdoms", None) is None:
        short_to_label = {k["short"]: k["label"] for k in C.load_kingdoms()}
        found: Dict[str, str] = {}
        for tsv in sorted(ctx.workdir.glob("*/*/ncbi_*.tsv")) + sorted(ctx.workdir.glob("*/*/*/ncbi_*.tsv")):
            short = tsv.stem.replace("ncbi_", "")
            if short not in short_to_label:
                continue
            try:
                with open(tsv, encoding="utf-8", errors="replace") as fh:
                    next(fh)
                    for line in fh:
                        acc = line.split("\t", 1)[0]
                        if acc:
                            found[acc] = short_to_label[short]
            except OSError:
                continue
        ctx._listing_kingdoms = found
    return ctx._listing_kingdoms


def build_rows(ctx: Context, listing: List[dict], downloads: Dict[str, Optional[Path]], meta: Dict[str, dict],
               qual: Dict[str, dict], red: Dict[str, dict], rna_counts: Dict[str, tuple]):
    th = ctx.th
    full_rows, filter_rows = [], []
    for lr in listing:
        acc = lr["accession"]
        m = meta.get(acc, {})
        path = downloads.get(acc)
        name = C.strip_fna(path.name) if path else \
            f"{acc}_{C.safe_asm_name(lr.get('asm_name') or m.get('assembly_name') or 'ASM')}_genomic"
        kingdom = _kingdom_from_lineage(ctx, m, lr["kingdom"], acc)
        twin = lr.get("twin_of_db") or ""
        twin_type = lr.get("twin_type") or (twin_kind(acc, twin) if twin else "")
        newer_version = twin_type == "newer_version"
        notes: List[str] = []
        reasons: List[str] = []

        # --- Downloaded
        if lr["download"] != "yes":
            downloaded = None
            reasons.append(f"Downloaded (skipped: {lr['download'][4:-1] if lr['download'].startswith('no (') else lr['download']})")
        else:
            downloaded = path is not None
            if not downloaded:
                reasons.append("Downloaded (download failed)")

        # --- Quality
        q = qual.get(name)
        if q is None and twin and twin in ctx.db_by_acc:
            dbr = ctx.db_by_acc[twin]
            q = {c: dbr.get(c, C.NOT_AVAILABLE) for c in Q.CHECKM2_COLUMNS.values()}
            q["Park Quality_Score"] = dbr.get("Park Quality_Score", C.NOT_AVAILABLE)
            notes.append(f"quality values copied from database twin {twin}")
        comp_ok = Q.passes_completeness(q["Completeness"], th) if q else None
        cont_ok = Q.passes_contamination(q["Contamination"], th) if q else None
        if q:
            if comp_ok is False:
                reasons.append(f"Completeness_ge70 ({q['Completeness']}%)")
            if cont_ok is False:
                reasons.append(f"Contamination_le10 ({q['Contamination']}%)")
        else:
            comp_ok = cont_ok = False
            if downloaded is not False:
                reasons.append("quality not assessed")

        # --- Redundancy
        r = red.get(name)
        if twin and (not newer_version or not downloaded):
            not_identical, not_species = False, False
            what = "newer version of" if newer_version else "GCA/GCF twin of"
            reasons.append(f"Not_identical_to_DB ({what} {twin})")
            reasons.append(f"Not_same_species_as_DB ({twin})")
        elif r:
            not_identical = r["identical_db"] is None
            not_species = r["species_db"] is None
            if not not_identical:
                reasons.append(f"Not_identical_to_DB ({C.name_to_acc(r['identical_db'])}, {r['identical_db_ani']})")
            if not not_species:
                reasons.append(f"Not_same_species_as_DB ({C.name_to_acc(r['species_db'])}, {r['species_db_ani']})")
            if r["species_db_note"]:
                notes.append(r["species_db_note"])
        else:
            not_identical = not_species = False
            if downloaded:
                reasons.append("redundancy not assessed")
        if r and r["batch_rep"] is not None:
            batch_ok = bool(r["batch_rep"])
            if not batch_ok:
                reasons.append(f"Representative_within_batch (redundant with {C.name_to_acc(r['batch_partner'])})")
        else:
            batch_ok = bool(downloaded) and r is not None
            if downloaded and r is None:
                pass

        filters = {"Downloaded": downloaded, "Completeness_ge70": comp_ok, "Contamination_le10": cont_ok,
                   "Not_identical_to_DB": not_identical, "Not_same_species_as_DB": not_species,
                   "Representative_within_batch": batch_ok}
        add = all(v is True for v in filters.values())

        # --- Newer version of a genome already in the database → Replace (if it still passes quality)
        recommendation = "Add" if add else "Do_not_add"
        if newer_version and twin in ctx.db_by_acc:
            old = ctx.db_by_acc[twin]
            change = ""
            if q:
                change = (f"completeness {old.get('Completeness')}→{q['Completeness']}, "
                          f"contamination {old.get('Contamination')}→{q['Contamination']}, "
                          f"genome size {old.get('Genome_Size')}→{q['Genome_Size']}")
            if downloaded and comp_ok and cont_ok:
                recommendation = "Replace"
                reasons = [f"newer version of database genome {twin} (assembly updated at NCBI "
                           f"{lr.get('release_date') or ''}); {change}".rstrip("; ")]
                notes.append(f"supersedes {twin}: replace the FASTA and table row of the old version")
            else:
                reasons = [f"newer version of database genome {twin} but " +
                           ("download failed" if not downloaded else f"fails the quality gate ({change})") +
                           "; the database keeps the previous version"] + \
                          [x for x in reasons if x.startswith(("Completeness", "Contamination"))]

        # --- Environment
        env_in = {"Metagenome_soource": _text(m.get("metagenome_source"), C.MISSING_METAGENOME),
                  "geo_loc_name": _text(m.get("geo_loc_name")), "Isolation_Source": _text(m.get("isolation_source")),
                  "Env_Broad": _text(m.get("env_broad_scale")), "Env_Local": _text(m.get("env_local_scale")),
                  "Env_medium": _text(m.get("env_medium"))}
        env_cat = classify(env_in)

        rr, tt = rna_counts.get(name, (None, None))
        row = {
            "Recommendation": recommendation,
            "Reason": "; ".join(reasons) if reasons else "passes all filters",
            "Name": name, "Assembly ID": acc, "Archaeal_Kingdom": kingdom,
            "Taxonomy": _text(m.get("full_lineage"), C.NOT_AVAILABLE),
            "Completeness": q["Completeness"] if q else C.NOT_AVAILABLE,
            "Contamination": q["Contamination"] if q else C.NOT_AVAILABLE,
            "Park Quality_Score": q["Park Quality_Score"] if q else C.NOT_AVAILABLE,
            "Bioproject": _text(m.get("bioproject_accession") or lr.get("bioproject")),
            "Biosample": _text(m.get("biosample_accession") or lr.get("biosample")),
            "Environment_Category-ArchaeaHQ": env_cat,
            "Metagenome_source": env_in["Metagenome_soource"], "geo_loc_name": env_in["geo_loc_name"],
            "lat_lon": _text(m.get("lat_lon")), "Isolation_Source": env_in["Isolation_Source"],
            "Env_Broad": env_in["Env_Broad"], "Env_Local": env_in["Env_Local"], "Env_medium": env_in["Env_medium"],
            "Coding_Density": q["Coding_Density"] if q else C.NOT_AVAILABLE,
            "Contig_N50": q["Contig_N50"] if q else C.NOT_AVAILABLE,
            "Average_Gene_Length": q["Average_Gene_Length"] if q else C.NOT_AVAILABLE,
            "Genome_Size": q["Genome_Size"] if q else C.NOT_AVAILABLE,
            "GC_Content": q["GC_Content"] if q else C.NOT_AVAILABLE,
            "Total_Coding_Sequences": q["Total_Coding_Sequences"] if q else C.NOT_AVAILABLE,
            "rRNA": C.NOT_AVAILABLE if rr is None or rr < 0 else rr,
            "tRNA": C.NOT_AVAILABLE if tt is None or tt < 0 else tt,
        }
        full_rows.append(row)
        overall = {"Add": "PASS", "Replace": "REPLACE", "Do_not_add": "FAIL"}[recommendation]
        frow = {"Assembly ID": acc, "Archaeal_Kingdom": kingdom, "Overall": overall}
        frow.update({k: report.pf(v) for k, v in filters.items()})
        frow["Notes"] = "; ".join(notes)
        filter_rows.append(frow)
    order = {k: i for i, k in enumerate(ctx.kingdom_labels)}
    rank = {"Add": 0, "Replace": 1, "Do_not_add": 2}
    full_rows.sort(key=lambda r: (rank[r["Recommendation"]], order.get(r["Archaeal_Kingdom"], 99), r["Assembly ID"]))
    orank = {"PASS": 0, "REPLACE": 1, "FAIL": 2}
    filter_rows.sort(key=lambda r: (orank[r["Overall"]], order.get(r["Archaeal_Kingdom"], 99), r["Assembly ID"]))
    return full_rows, filter_rows


def stage_report(ctx: Context, listing, downloads, meta, qual, red, rna_counts, out_dir: Path) -> None:
    full_rows, filter_rows = build_rows(ctx, listing, downloads, meta, qual, red, rna_counts)
    date = time.strftime("%Y-%m-%d")
    full_path = ctx.run_dir / f"ArchaeaHQ_update_{date}_full_table.tsv"
    filt_path = ctx.run_dir / f"ArchaeaHQ_update_{date}_filters.tsv"
    report.write_full_table(full_path, full_rows)
    report.write_filter_table(filt_path, filter_rows)
    kingdoms_present = ctx.kingdom_labels + sorted({r["Archaeal_Kingdom"] for r in filter_rows} - set(ctx.kingdom_labels))
    counts = report.stage_counts_from_rows(kingdoms_present, filter_rows)
    summary = report.write_summary(ctx.run_dir / "summary.txt", ctx.run_name, kingdoms_present, full_rows, filter_rows, counts)
    sankey_in = ctx.run_dir / "sankey_input.tsv"
    report.write_sankey_input(sankey_in, ctx.run_name, kingdoms_present, counts)
    report.append_ledger(ctx.ledger, ctx.run_name, full_rows)

    stages = list(counts.keys())
    per_k = OrderedDict((k, [counts[s].get(k, 0) for s in stages]) for k in kingdoms_present if counts[stages[0]].get(k, 0))
    ui.kingdom_table("Genomes passing each filter (cumulative)", stages, per_k,
                     total_row=[sum(counts[s].values()) for s in stages])
    n_add = sum(1 for r in full_rows if r["Recommendation"] == "Add")
    n_replace = sum(1 for r in full_rows if r["Recommendation"] == "Replace")
    reasons = Counter(r["Reason"].split(";")[0].split(" (")[0].strip() for r in full_rows
                      if r["Recommendation"] == "Do_not_add")
    outputs = [str(full_path), str(filt_path), str(ctx.run_dir / "summary.txt"), str(sankey_in), f"{ctx.ledger} (appended)"]
    if ctx.args.sankey:
        try:
            import sankey_generic  # noqa: F401
            C.run_cmd([sys.executable, str(HERE / "lib" / "sankey_generic.py"), str(sankey_in),
                       str(ctx.run_dir / "sankey")], log_path=ctx.log_path)
            outputs.append(str(ctx.run_dir / "sankey.html") + " / .png")
        except Exception as e:  # optional feature
            ui.warn(f"Sankey figure not rendered: {e}")
    ui.verdict_panel(n_add, len(full_rows) - n_add - n_replace, dict(reasons), outputs, replace=n_replace)


# ═══════════════════════════════════════════════════════════════════════════════
# Commands
# ═══════════════════════════════════════════════════════════════════════════════

def cmd_setup(ctx: Context) -> int:
    ctx.log_path = ctx.workdir / "setup.log"
    C.setup_file_logging(ctx.workdir / "archaeahq_update.log")
    with ui.stage(1, 1, STAGE_TITLES["env"]):
        stage_env(ctx, envcheck.REQUIRED_TOOLS, need_checkm2_db=True)
        C.save_json(ctx.workdir / "environment.json", envcheck.environment_summary(ctx.tools, ctx.checkm2_db))
    ui.ok("everything needed is in place")
    ui.kv("work directory", ctx.workdir)
    ui.kv("database", ctx.describe_db())
    return 0


def cmd_fetch_db(ctx: Context) -> int:
    a = ctx.args
    ctx.log_path = ctx.workdir / "fetch.log"
    C.setup_file_logging(ctx.workdir / "archaeahq_update.log")
    dest = ctx.releases_dir / C.release_folder_name(fetchdb.FIGSHARE_VERSION)
    extras = [x.strip().lower() for x in a.extras.split(",") if x.strip()]
    keys = {v[0]: k for k, v in fetchdb.FILES.items()}
    bad = [x for x in extras if x not in keys or x == "fna"]
    if bad:
        raise SystemExit(f"unknown --extras {', '.join(bad)}; choose from faa, 16s, tables")
    wanted = ["fna.zip"] + [keys[x] for x in extras]

    ui.kv("dataset", f"ArchaeaHQ {fetchdb.FIGSHARE_VERSION}  ·  {fetchdb.FIGSHARE_PAGE}")
    ui.kv("DOI", fetchdb.FIGSHARE_DOI)
    ui.kv("destination", dest)
    if C.is_release_folder(dest):
        n = sum(1 for _ in (dest / "fna").glob("*.fna")) if (dest / "fna").is_dir() else 0
        ui.ok(f"{dest} is already set up ({n:,} genomes) — nothing to download")
        if not extras:
            _fetch_next_steps(ctx)
            return 0
    if ctx.current_release and ctx.current_release != dest:
        ui.warn(f"a newer database version already exists here: {ctx.current_release}")
        if not a.yes and not ui.confirm("download v1.0 anyway?", default=False):
            return 1

    with ui.spinner("asking figshare for the file list"):
        files = fetchdb.article_files()
    missing = [w for w in wanted if w not in files]
    if missing:
        raise SystemExit(f"figshare article {fetchdb.FIGSHARE_ARTICLE} has no file named {', '.join(missing)}")
    rows, dl_gb, unpacked_gb = [], 0.0, 0.0
    for w in wanted:
        f, (key, what, unpacked) = files[w], fetchdb.FILES[w]
        if not (w == "fna.zip" and a.from_zip):
            dl_gb += f["size"] / 1e9
        rows.append([w, what, f"{f['size'] / 1e9:.2f}", f"{unpacked:g}"])
        unpacked_gb += unpacked
    ui.table("Files", ["file", "content", "download (GB)", "unpacked (GB)"], rows,
             justify={"content": "left"})
    free = fetchdb.free_space_gb(dest)
    need_gb = dl_gb + unpacked_gb
    after_gb = need_gb if a.keep_zip else unpacked_gb
    ui.kv("disk needed", f"~{need_gb:.0f} GB during setup, ~{after_gb:.0f} GB afterwards")
    ui.kv("disk free", f"{free:.0f} GB at {dest.parent}")
    if free < need_gb:
        ui.warn("not enough free space for the download plus the unpacked files")
        if not a.yes and not ui.confirm("continue anyway?", default=False):
            return 1
    if not a.yes and not ui.confirm(f"download and set up ArchaeaHQ {fetchdb.FIGSHARE_VERSION} in {dest}?"):
        ui.info("nothing done")
        return 1

    dest.mkdir(parents=True, exist_ok=True)
    total = 2 + len(wanted) + (0 if a.no_verify else 1)
    step = 0
    zips: Dict[str, Path] = {}
    for w in wanted:
        f = files[w]
        step += 1
        with ui.stage(step, total, f"Downloading {w} ({f['size'] / 1e9:.2f} GB)"):
            if w == "fna.zip" and a.from_zip:
                zips[w] = Path(a.from_zip).expanduser().resolve()
                if not zips[w].is_file():
                    raise SystemExit(f"--from-zip {zips[w]} not found")
                ui.info(f"using {zips[w]} (no download)")
                continue
            target = dest / w
            with ui.progress() as prog:
                task = prog.add_task(f"{w} (MB)", total=max(1, f["size"] // 1_000_000))
                got = [0]

                def cb(n, got=got, task=task, prog=prog):
                    got[0] += n
                    prog.update(task, completed=got[0] // 1_000_000)
                fetchdb.download(f["download_url"], target, f["size"], progress_cb=cb)
            zips[w] = target
            ui.ok(f"{target} ({target.stat().st_size / 1e9:.2f} GB)")

    if not a.no_verify:
        step += 1
        with ui.stage(step, total, "Verifying checksums"):
            for w, path in zips.items():
                if w == "fna.zip" and a.from_zip and path.stat().st_size != files[w]["size"]:
                    ui.warn(f"{path.name}: size differs from figshare's copy; skipping the MD5 check")
                    continue
                with ui.spinner(f"MD5 of {path.name}"):
                    got = fetchdb.md5_of(path)
                if files[w]["md5"] and got != files[w]["md5"]:
                    raise SystemExit(f"{path.name}: MD5 {got} does not match figshare's {files[w]['md5']} — "
                                     f"delete the file and run fetch-db again")
                ui.ok(f"{path.name}: MD5 ok")

    step += 1
    with ui.stage(step, total, "Unpacking the genomes"):
        members = fetchdb.zip_members(zips["fna.zip"], suffixes=(".fna",))
        with ui.progress() as prog:
            task = prog.add_task("genomes", total=len(members))
            fetchdb.extract_flat(zips["fna.zip"], dest / "fna", members, progress_cb=lambda n: prog.update(task, advance=n))
        ui.ok(f"{len(members):,} FASTA files in {dest / 'fna'}")
        for w, path in zips.items():
            if w == "faa.zip":
                mem = fetchdb.zip_members(path, suffixes=(".faa",))
                with ui.progress() as prog:
                    task = prog.add_task("proteins", total=len(mem))
                    fetchdb.extract_flat(path, dest / "faa", mem, progress_cb=lambda n: prog.update(task, advance=n))
                ui.ok(f"{len(mem):,} protein files in {dest / 'faa'}")
        if not a.keep_zip:
            for w, path in zips.items():
                if path.parent == dest and path.suffix == ".zip":
                    path.unlink()
            ui.info("zip archive(s) deleted (--keep-zip keeps them)")

    step += 1
    with ui.stage(step, total, f"Registering {fetchdb.FIGSHARE_VERSION} as the current database"):
        table = C.read_tsv(C.DB_TABLE_PATH)
        have = {C.strip_fna(p.name) for p in (dest / "fna").glob("*.fna")}
        want = {C.strip_fna(r["Name"]) for r in table}
        if have != want:
            ui.warn(f"{len(want - have)} table row(s) without FASTA and {len(have - want)} FASTA file(s) not in the table")
        shutil.copy2(C.DB_TABLE_PATH, dest / C.RELEASE_TABLE)
        C.save_json(dest / C.RELEASE_META, {
            "version": fetchdb.FIGSHARE_VERSION, "date": time.strftime("%Y-%m-%d"), "previous_version": None,
            "source": fetchdb.FIGSHARE_PAGE, "doi": fetchdb.FIGSHARE_DOI,
            "files": {w: files[w]["md5"] for w in wanted}, "genomes": len(table), "tool_version": C.VERSION})
        ui.ok(f"{len(table):,} genomes registered in {dest}")
    _fetch_next_steps(ctx)
    return 0


def _fetch_next_steps(ctx: Context) -> None:
    ui.console().print()
    ui.ok("ArchaeaHQ is ready. To bring it up to date with NCBI:")
    ui.kv("1. what is new", "python3 archaeahq_update.py check")
    ui.kv("2. evaluate it", "python3 archaeahq_update.py run --threads 24")
    ui.kv("3. next version", "python3 archaeahq_update.py release")
    ui.info("the tools are installed with `python3 archaeahq_update.py setup` if you have not done it yet")


def cmd_check(ctx: Context) -> int:
    ctx.run_name = f"check_{time.strftime('%Y-%m-%d')}"
    ctx.run_dir = ctx.workdir / "checks" / ctx.run_name
    ctx.run_dir.mkdir(parents=True, exist_ok=True)
    ctx.log_path = ctx.run_dir / "commands.log"
    C.setup_file_logging(ctx.run_dir / "archaeahq_update.log")
    with ui.stage(1, 2, STAGE_TITLES["env"]):
        stage_env(ctx, ["datasets", "dataformat"], need_checkm2_db=False)
    with ui.stage(2, 2, STAGE_TITLES["list"]):
        out = ctx.run_dir
        rows = stage_list(ctx, out)
    ui.ok(f"{len(rows)} new genome(s) listed in {out / 'new_accessions.tsv'}")
    ui.info("run `archaeahq_update.py run` to evaluate them")
    return 0


def cmd_sketch_db(ctx: Context) -> int:
    ctx.log_path = ctx.workdir / "sketch.log"
    C.setup_file_logging(ctx.workdir / "archaeahq_update.log")
    with ui.stage(1, 2, STAGE_TITLES["env"]):
        stage_env(ctx, ["skani"], need_checkm2_db=False)
    with ui.stage(2, 2, "Sketching the database genomes"):
        db_fna = Path(ctx.args.db_fna).expanduser() if ctx.args.db_fna else ctx.default_db_fna()
        if not db_fna:
            raise SystemExit("no release folder found: pass --db-fna <folder with the database .fna files>")
        files = R.db_fna_files(db_fna)
        out = Path(ctx.args.out).expanduser()
        if out.exists():
            shutil.rmtree(out)
        lst = out.parent / (out.name + ".list.txt")
        out.parent.mkdir(parents=True, exist_ok=True)
        C.write_lines(lst, (str(p) for p in files))
        with ui.spinner(f"skani sketch of {len(files):,} genomes"):
            C.run_cmd(ctx.tools["skani"]("sketch", "-l", lst, "-o", out, "-t", ctx.args.threads), log_path=ctx.log_path)
        (out / ".complete").write_text(f"{len(files)} genomes\n")
    ui.ok(f"sketch written to {out} — pass it to `run --db-sketch {out}`")
    return 0


def cmd_run(ctx: Context) -> int:
    a = ctx.args
    if not a.db_fna and not a.db_sketch:
        d = ctx.default_db_fna()
        if not d:
            raise SystemExit("no release folder found: pass --db-fna <folder with the database .fna files> "
                             "or --db-sketch <skani sketch folder>")
        a.db_fna = str(d)
    ctx.run_name = a.run_name or f"run_{time.strftime('%Y-%m-%d')}"
    ctx.run_dir = ctx.runs_dir / ctx.run_name
    ctx.run_dir.mkdir(parents=True, exist_ok=True)
    ctx.log_path = ctx.run_dir / "commands.log"
    C.setup_file_logging(ctx.run_dir / "archaeahq_update.log")
    ui.kv("run folder", ctx.run_dir)
    ui.kv("database", ctx.describe_db())
    ui.kv("database FASTA", a.db_fna or f"sketch {a.db_sketch}")
    ui.kv("threads", a.threads)

    stages = list(STAGE_TITLES)
    if a.from_stage:
        for s in stages[stages.index(a.from_stage):]:
            C.clear_stage(ctx.stage_dir(s))
    total = len(stages)
    state_path = ctx.run_dir / "state.json"
    state = C.load_json(state_path) if state_path.exists() else {}

    def done(stage):
        return C.stage_is_done(ctx.stage_dir(stage))

    # 0 env
    with ui.stage(1, total, STAGE_TITLES["env"]):
        stage_env(ctx, envcheck.REQUIRED_TOOLS, need_checkm2_db=True)
        C.mark_stage_done(ctx.stage_dir("env"))

    # 1 list
    d = ctx.stage_dir("list")
    with ui.stage(2, total, STAGE_TITLES["list"]):
        if done("list"):
            listing = C.read_tsv(d / "new_accessions.tsv")
            ui.info(f"reusing {len(listing)} accession(s) from a previous attempt")
        else:
            listing = stage_list(ctx, d)
            C.mark_stage_done(d, {"n_new": len(listing)})
    if not listing:
        ui.ok("no new genomes to evaluate — the database is up to date for the selected kingdoms")
        return 0
    accessions = [r["accession"] for r in listing]

    # 2 download
    d = ctx.stage_dir("download")
    with ui.stage(3, total, STAGE_TITLES["download"]):
        if done("download"):
            downloads = {r["accession"]: (Path(r["path"]) if r["status"] == "ok" else None)
                         for r in C.read_tsv(d / "download_status.tsv")}
            ui.info(f"reusing {sum(1 for p in downloads.values() if p)} downloaded genome(s)")
        else:
            downloads = stage_download(ctx, listing, d)
            C.mark_stage_done(d)
    fna_dir = d / "fna"
    fnas = sorted(p for p in downloads.values() if p)
    names_acc = {p.stem: acc for acc, p in downloads.items() if p}

    # 3 metadata
    d = ctx.stage_dir("metadata")
    with ui.stage(4, total, STAGE_TITLES["metadata"]):
        if done("metadata") and (d / "combined.tsv").exists():
            meta = {r["assembly_accession"]: r for r in C.read_tsv(d / "combined.tsv")}
            ui.info(f"reusing metadata for {len(meta)} accession(s)")
        else:
            meta = stage_metadata(ctx, accessions, d)
            C.mark_stage_done(d)

    # 4 quality
    d = ctx.stage_dir("quality")
    with ui.stage(5, total, STAGE_TITLES["quality"]):
        if done("quality") and (d / "quality_report.tsv").exists():
            qual = Q.parse_quality_report(d / "quality_report.tsv")
            ui.info(f"reusing CheckM2 results for {len(qual)} genome(s)")
        elif done("quality"):
            qual = {}
        else:
            qual = stage_quality(ctx, fnas, fna_dir, d)
            C.mark_stage_done(d)

    # 5 redundancy
    d = ctx.stage_dir("redundancy")
    with ui.stage(6, total, STAGE_TITLES["redundancy"]):
        if done("redundancy") and (d / "redundancy.tsv").exists():
            red = _load_redundancy(d / "redundancy.tsv")
            ui.info(f"reusing skani results for {len(red)} genome(s)")
        elif done("redundancy"):
            red = {}
        else:
            red = stage_redundancy(ctx, fnas, names_acc, qual, d)
            C.mark_stage_done(d)

    # 6 rna
    d = ctx.stage_dir("rna")
    with ui.stage(7, total, STAGE_TITLES["rna"]):
        if a.no_rna:
            ui.info("skipped (--no-rna)")
            rna_counts = {}
        elif done("rna") and (d / "rna_counts.tsv").exists():
            rna_counts = {r["Name"]: (int(r["rRNA"]), int(r["tRNA"])) for r in C.read_tsv(d / "rna_counts.tsv")}
            ui.info(f"reusing RNA counts for {len(rna_counts)} genome(s)")
        else:
            rna_counts = stage_rna(ctx, fnas, d)
            C.mark_stage_done(d)

    # 7 report
    with ui.stage(8, total, STAGE_TITLES["report"]):
        stage_report(ctx, listing, downloads, meta, qual, red, rna_counts, ctx.stage_dir("report"))
        C.mark_stage_done(ctx.stage_dir("report"))
    C.save_json(state_path, {"run": ctx.run_name, "finished": time.strftime("%Y-%m-%d %H:%M:%S"), "n_new": len(listing)})
    return 0


def _load_redundancy(path: Path) -> Dict[str, dict]:
    out = {}
    for r in C.read_tsv(path):
        out[r["Name"]] = {
            "identical_db": r["identical_to_db"] or None, "identical_db_ani": r["identical_evidence"] or None,
            "species_db": r["same_species_as_db"] or None, "species_db_ani": r["species_evidence"] or None,
            "species_db_note": r["note"],
            "batch_rep": None if r["batch_representative"] == "" else (r["batch_representative"] == "yes"),
            "batch_partner": r["batch_partner"] or None,
            "batch_cluster": [x for x in r["batch_cluster"].split("; ") if x],
        }
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# release — build the next database version from the recommendations
# ═══════════════════════════════════════════════════════════════════════════════

CHANGE_COLUMNS = ["Assembly ID", "Action", "Archaeal_Kingdom", "Name", "Replaces", "Source_run", "Note"]


def _fna_index(folder: Path) -> Dict[str, Path]:
    """genome name (extension stripped) → FASTA path, for every FASTA file directly in `folder`."""
    idx: Dict[str, Path] = {}
    for p in sorted(folder.iterdir()):
        if p.is_file() and (p.suffix in (".fna", ".fa", ".fasta") or p.name.endswith(".fna.gz")):
            idx.setdefault(C.strip_fna(p.name), p)
    return idx


def _newest_run(ctx: Context) -> Optional[Path]:
    runs = [d for d in ctx.runs_dir.glob("*") if d.is_dir() and list(d.glob("*_full_table.tsv"))] \
        if ctx.runs_dir.exists() else []
    return max(runs, key=lambda d: d.stat().st_mtime) if runs else None


def _place(src: Path, dst: Path, copy: bool) -> str:
    """Hard-link (or copy) src to dst; returns 'link' or 'copy'."""
    if dst.exists() or dst.is_symlink():
        dst.unlink()
    if not copy:
        try:
            os.link(src, dst)
            return "link"
        except OSError:
            pass
    shutil.copy2(src, dst)
    return "copy"


def _mean(vals: List[float]) -> Optional[float]:
    return sum(vals) / len(vals) if vals else None


def _median(vals: List[float]) -> Optional[float]:
    if not vals:
        return None
    s = sorted(vals)
    n = len(s)
    return s[n // 2] if n % 2 else (s[n // 2 - 1] + s[n // 2]) / 2


def _db_stats(rows: List[dict], files: Dict[str, Path]) -> Dict[str, Optional[float]]:
    col = lambda c: [v for v in (C.to_float(r.get(c)) for r in rows) if v is not None]  # noqa: E731
    rrna = [C.to_float(r.get("rRNA")) for r in rows]
    with_rrna = [v for v in rrna if v is not None]
    size = 0
    for r in rows:
        p = files.get(C.strip_fna(r["Name"]))
        if p:
            try:
                size += p.stat().st_size
            except OSError:
                pass
    return {
        "genomes": float(len(rows)),
        "mean completeness (%)": _mean(col("Completeness")),
        "mean contamination (%)": _mean(col("Contamination")),
        "mean Park quality score": _mean(col("Park Quality_Score")),
        "median genome size (Mb)": (lambda m: m / 1e6 if m is not None else None)(_median(col("Genome_Size"))),
        "median contig N50 (kb)": (lambda m: m / 1e3 if m is not None else None)(_median(col("Contig_N50"))),
        "mean GC content": _mean(col("GC_Content")),
        "genomes with ≥1 rRNA gene (%)": (100.0 * sum(1 for v in with_rrna if v > 0) / len(with_rrna)) if with_rrna else None,
        "FASTA size (GB)": size / 1e9,
    }


def _fmt_stat(name: str, v: Optional[float]) -> str:
    if v is None:
        return C.NOT_AVAILABLE
    if name == "genomes":
        return f"{int(v):,}"
    return f"{v:,.2f}"


def _pct(old: int, new: int) -> str:
    return f"{100.0 * (new - old) / old:+.1f}%" if old else ""


def _fmt_delta(name: str, old: Optional[float], new: Optional[float]) -> str:
    if old is None or new is None:
        return ""
    d = new - old
    if name == "genomes":
        pct = f" ({100.0 * d / old:+.1f}%)" if old else ""
        return f"{int(d):+,}{pct}"
    return f"{d:+,.2f}"


def _release_compare(ctx: Context, old_label: str, new_label: str, old_rows: List[dict], new_rows: List[dict],
                     old_files: Dict[str, Path], new_files: Dict[str, Path], changes: List[dict]) -> List[str]:
    """Print the old-vs-new comparison; return the same content as plain text lines."""
    lines: List[str] = [f"ArchaeaHQ {old_label} → {new_label}", ""]

    # 1. genomes per kingdom
    old_k = Counter(r["Archaeal_Kingdom"] for r in old_rows)
    new_k = Counter(r["Archaeal_Kingdom"] for r in new_rows)
    added_k = Counter(c["Archaeal_Kingdom"] for c in changes if c["Action"] == "added")
    repl_k = Counter(c["Archaeal_Kingdom"] for c in changes if c["Action"] == "replaced")
    kingdoms = ctx.kingdom_labels + sorted((set(old_k) | set(new_k)) - set(ctx.kingdom_labels))
    per_k = OrderedDict()
    for k in kingdoms:
        if old_k[k] or new_k[k]:
            per_k[k] = [f"{old_k[k]:,}", f"{added_k[k]:,}", f"{repl_k[k]:,}", f"{new_k[k]:,}", _pct(old_k[k], new_k[k])]
    tot = [f"{len(old_rows):,}", f"{sum(added_k.values()):,}", f"{sum(repl_k.values()):,}", f"{len(new_rows):,}",
           _pct(len(old_rows), len(new_rows))]
    cols = [old_label, "added", "replaced", new_label, "growth"]
    ui.kingdom_table("Genomes per kingdom", cols, per_k, total_row=tot)
    lines += _plain_table(["Kingdom"] + cols, [[k] + v for k, v in per_k.items()] + [["Total"] + tot])

    # 2. summary statistics
    so, sn = _db_stats(old_rows, old_files), _db_stats(new_rows, new_files)
    rows = [[name, _fmt_stat(name, so[name]), _fmt_stat(name, sn[name]), _fmt_delta(name, so[name], sn[name])]
            for name in so]
    ui.table("Database summary", ["metric", old_label, new_label, "change"], rows)
    lines += [""] + _plain_table(["metric", old_label, new_label, "change"], rows)

    # 3. environment categories
    old_e = Counter(r.get("Environment_Category-ArchaeaHQ", "") for r in old_rows)
    new_e = Counter(r.get("Environment_Category-ArchaeaHQ", "") for r in new_rows)
    cats = sorted(set(old_e) | set(new_e), key=lambda c: (-new_e[c], c))
    rows = [[c or "(empty)", f"{old_e[c]:,}", f"{new_e[c]:,}", _fmt_delta("genomes", float(old_e[c]), float(new_e[c]))]
            for c in cats]
    ui.table("Environment categories", ["category", old_label, new_label, "change"], rows)
    lines += [""] + _plain_table(["category", old_label, new_label, "change"], rows)
    return lines


def _plain_table(columns: List[str], rows: List[List[str]]) -> List[str]:
    widths = [max(len(str(x)) for x in col) for col in zip(columns, *rows)] if rows else [len(c) for c in columns]
    fmt = lambda r: "  ".join(str(x).ljust(w) if i == 0 else str(x).rjust(w) for i, (x, w) in enumerate(zip(r, widths)))  # noqa: E731
    return [fmt(columns), "  ".join("-" * w for w in widths)] + [fmt(r) for r in rows]


def cmd_release(ctx: Context) -> int:
    a = ctx.args
    ctx.log_path = ctx.workdir / "release.log"
    C.setup_file_logging(ctx.workdir / "archaeahq_update.log")
    total = 3 if a.dry_run else 5

    # ── inputs ───────────────────────────────────────────────────────────────
    db_fna = Path(a.db_fna).expanduser().resolve() if a.db_fna else ctx.default_db_fna()
    if not db_fna:
        raise SystemExit(f"no release folder found in {ctx.releases_dir}: pass --db-fna <FASTA folder of the current version>")
    if not db_fna.is_dir():
        raise SystemExit(f"--db-fna {db_fna} is not a folder")
    old_table = ctx.db_table
    old_label = ctx.db_version
    if a.release_version:
        v = C.parse_version(a.release_version)
        if not v:
            raise SystemExit(f"--version must look like v1.2 (got {a.release_version!r})")
        new_label = C.version_label(v)
    else:
        try:
            new_label = C.next_version(old_label)
        except ValueError as e:
            raise SystemExit(str(e))
    if C.parse_version(new_label) <= (C.parse_version(old_label) or (0, 0)):
        raise SystemExit(f"new version {new_label} is not newer than the current {old_label}")
    out = Path(a.out).expanduser().resolve() if a.out else ctx.releases_dir / C.release_folder_name(new_label)
    if out == db_fna or out in db_fna.parents or db_fna in out.parents:
        raise SystemExit("the new version folder must be separate from the current FASTA folder")
    if out.exists() and any(out.iterdir()) and not a.dry_run:
        raise SystemExit(f"{out} exists and is not empty — remove it or pass --version/--out")
    prev_release = old_table.parent if C.is_release_folder(old_table.parent) and db_fna == old_table.parent / "fna" else None
    if a.table:
        table = Path(a.table).expanduser().resolve()
        run_dir = table.parent
    else:
        run_dir = ctx.runs_dir / a.run if a.run else _newest_run(ctx)
        if run_dir is None or not run_dir.is_dir():
            raise SystemExit(f"no finished run found in {ctx.runs_dir} — run `archaeahq_update.py run` first or pass --table")
        tables = sorted(run_dir.glob("*_full_table.tsv"))
        if not tables:
            raise SystemExit(f"{run_dir} has no *_full_table.tsv (the run did not reach the report stage)")
        table = tables[-1]
    new_fna = Path(a.new_fna).expanduser().resolve() if a.new_fna else run_dir / "02_download" / "fna"
    if not new_fna.is_dir():
        raise SystemExit(f"folder with the new genomes not found: {new_fna} (pass --new-fna)")
    ui.kv("current version", ctx.describe_db())
    ui.kv("current FASTA", db_fna)
    ui.kv("recommendations", table)
    ui.kv("new genomes FASTA", new_fna)
    ui.kv("new version", f"{new_label} → {out}" + ("  [dry run]" if a.dry_run else ""))

    # ── 1. current database ──────────────────────────────────────────────────
    with ui.stage(1, total, "Reading the current database"):
        with ui.spinner("indexing the FASTA folder"):
            old_files = _fna_index(db_fna)
        missing_old = [r["Assembly ID"] for r in ctx.db_rows if C.strip_fna(r["Name"]) not in old_files]
        orphans = set(old_files) - {C.strip_fna(r["Name"]) for r in ctx.db_rows}
        ui.ok(f"{len(ctx.db_rows):,} table rows, {len(old_files):,} FASTA files")
        if missing_old:
            ui.warn(f"{len(missing_old)} table row(s) have no FASTA file in {db_fna}; their rows are kept, no file is written")
        if orphans:
            ui.warn(f"{len(orphans)} FASTA file(s) are not in the table and are not carried over")

    # ── 2. selection ─────────────────────────────────────────────────────────
    with ui.stage(2, total, "Selecting the recommended genomes"):
        rec_rows = C.read_tsv(table)
        new_files = _fna_index(new_fna)
        accept = None
        if a.accept:
            accept = {ln.strip() for ln in Path(a.accept).read_text().splitlines() if ln.strip() and not ln.startswith("#")}
            ui.info(f"{len(accept)} accession(s) in the --accept list")
        changes: List[dict] = []
        additions: List[tuple] = []          # (row, src_path, replaces_acc)
        seen = set()
        for r in rec_rows:
            rec, acc = r.get("Recommendation", ""), r["Assembly ID"]
            if rec not in ("Add", "Replace") or acc in seen:
                continue
            seen.add(acc)
            ch = {"Assembly ID": acc, "Action": "", "Archaeal_Kingdom": r["Archaeal_Kingdom"], "Name": r["Name"],
                  "Replaces": "", "Source_run": run_dir.name, "Note": ""}
            src = new_files.get(C.strip_fna(r["Name"]))
            if accept is not None and acc not in accept:
                ch.update(Action="skipped", Note="not in the --accept list")
            elif acc in ctx.db_acc:
                ch.update(Action="skipped", Note="already in the database")
            elif src is None:
                ch.update(Action="skipped", Note=f"FASTA not found in {new_fna}")
            else:
                old = ctx.db_numeric.get(C.acc_numeric(acc), "")
                if rec == "Replace" and old in ctx.db_by_acc:
                    ch.update(Action="replaced", Replaces=old)
                    additions.append((r, src, old))
                else:
                    ch.update(Action="added", Note="previous version not in the database; added" if rec == "Replace" else "")
                    additions.append((r, src, ""))
            changes.append(ch)
        n_add = sum(1 for c in changes if c["Action"] == "added")
        n_rep = sum(1 for c in changes if c["Action"] == "replaced")
        n_skip = sum(1 for c in changes if c["Action"] == "skipped")
        ui.ok(f"{n_add:,} genome(s) to add, {n_rep} to replace" + (f", {n_skip} skipped" if n_skip else ""))
        for note, n in Counter(c["Note"] for c in changes if c["Action"] == "skipped").items():
            ui.warn(f"{n} skipped: {note}")
        if not additions:
            ui.warn("nothing to apply — the new version would be identical to the current one")
            if not a.dry_run:
                ui.info("no folder written, nothing removed")
                return 0
        replaced = {old: r for r, _, old in additions if old}
        for old in replaced:
            changes.append({"Assembly ID": old, "Action": "removed", "Archaeal_Kingdom": ctx.db_by_acc[old]["Archaeal_Kingdom"],
                            "Name": ctx.db_by_acc[old]["Name"], "Replaces": "", "Source_run": run_dir.name,
                            "Note": f"superseded by {replaced[old]['Assembly ID']}"})

    # new table: old order kept, a replacement takes the place of the old row, additions appended per kingdom
    new_rows: List[dict] = []
    for r in ctx.db_rows:
        acc = r["Assembly ID"]
        new_rows.append({c: replaced[acc].get(c, "") for c in C.DB_COLUMNS} if acc in replaced else r)
    order = {k: i for i, k in enumerate(ctx.kingdom_labels)}
    for r, _, old in sorted(additions, key=lambda t: (order.get(t[0]["Archaeal_Kingdom"], 99), t[0]["Assembly ID"])):
        if not old:
            new_rows.append({c: r.get(c, "") for c in C.DB_COLUMNS})
    src_by_name = {C.strip_fna(r["Name"]): src for r, src, _ in additions}
    new_files_all = {**{C.strip_fna(r["Name"]): old_files[C.strip_fna(r["Name"])] for r in new_rows
                        if C.strip_fna(r["Name"]) in old_files and C.strip_fna(r["Name"]) not in src_by_name},
                     **src_by_name}

    # ── 3. build the folder ──────────────────────────────────────────────────
    stage_no = 3
    if not a.dry_run:
        with ui.stage(stage_no, total, f"Building {out}"):
            stage_no += 1
            fna_out = out / "fna"
            fna_out.mkdir(parents=True, exist_ok=True)
            n_link = n_copy = 0
            with ui.progress() as prog:
                task = prog.add_task("copying FASTA files" if a.copy else "linking FASTA files", total=len(new_files_all))
                for name, src in new_files_all.items():
                    how = _place(src, fna_out / src.name, a.copy)
                    n_link += how == "link"
                    n_copy += how == "copy"
                    prog.update(task, advance=1)
            C.write_tsv(out / "ArchaeaHQ-Info.tsv", new_rows, C.DB_COLUMNS)
            C.write_tsv(out / "release_changes.tsv", changes, CHANGE_COLUMNS)
            C.save_json(out / "release.json", {
                "version": new_label, "date": time.strftime("%Y-%m-%d"), "previous_version": old_label,
                "previous_table": str(old_table), "previous_fna": str(db_fna),
                "previous_removed": False, "source_run": run_dir.name,
                "source_table": str(table), "genomes": len(new_rows),
                "added": n_add, "replaced": n_rep, "removed": len(replaced), "skipped": n_skip,
                "fasta_hardlinked": n_link, "fasta_copied": n_copy, "tool_version": C.VERSION})
            ui.ok(f"{len(new_files_all):,} FASTA files in {fna_out} ({n_link:,} hard-linked, {n_copy:,} copied)")
            ui.ok(f"{len(new_rows):,} rows in {out / 'ArchaeaHQ-Info.tsv'}")

    # ── 4. comparison ────────────────────────────────────────────────────────
    with ui.stage(stage_no, total, f"Comparing {old_label} with {new_label}"):
        lines = _release_compare(ctx, old_label, new_label, ctx.db_rows, new_rows, old_files, new_files_all, changes)
        if not a.dry_run:
            C.write_lines(out / "comparison.txt", lines)
            ui.kv("comparison", out / "comparison.txt")
            ui.kv("changes", out / "release_changes.tsv")
        else:
            ui.info("dry run: nothing written")

    # ── 5. remove the superseded version ─────────────────────────────────────
    if a.dry_run:
        return 0
    with ui.stage(stage_no + 1, total, f"Removing the superseded version {old_label}"):
        removed = _remove_previous(ctx, out, prev_release, db_fna, old_files, new_files_all, orphans, missing_old,
                                   keep=a.keep_previous)
        meta = C.load_json(out / C.RELEASE_META)
        meta["previous_removed"] = removed
        C.save_json(out / C.RELEASE_META, meta)
    ok_next = out.parent == ctx.releases_dir and out.name == C.release_folder_name(new_label)
    if ok_next:
        ui.ok(f"{new_label} is now the current database: `run`, `check` and `release` use {out} automatically")
    else:
        ui.info(f"next runs: `archaeahq_update.py run --db-table {out / C.RELEASE_TABLE} --db-fna {out / 'fna'}`")
    return 0


def _remove_previous(ctx: Context, out: Path, prev_release: Optional[Path], db_fna: Path, old_files: Dict[str, Path],
                     new_files: Dict[str, Path], orphans: set, missing_old: List[str], keep: bool) -> bool:
    """Delete the previous version once every genome it contributed is verified in the new folder. True if removed."""
    target = prev_release or db_fna
    what = f"release folder {target}" if prev_release else f"FASTA folder {target}"
    protected = {C.BUNDLE_DIR.resolve(), C.DATA_DIR.resolve(), ctx.workdir, ctx.releases_dir, out}
    if keep:
        ui.info(f"kept (--keep-previous): {what}")
        return False
    if target.resolve() in protected or any(p == target.resolve() for p in out.parents) or target == Path("/"):
        ui.warn(f"not removed (protected location): {what}")
        return False
    if orphans:
        ui.warn(f"not removed: {what} holds {len(orphans)} FASTA file(s) that are not in the table; delete it by hand")
        return False
    bad = []
    for name, src in old_files.items():
        if name not in new_files:
            continue                          # replaced genome: intentionally dropped
        dst = out / "fna" / src.name
        try:
            if not dst.is_file() or dst.stat().st_size != src.stat().st_size:
                bad.append(name)
        except OSError:
            bad.append(name)
    if bad:
        ui.warn(f"not removed: {len(bad)} genome(s) of the previous version are missing or differ in the new folder")
        return False
    with ui.spinner(f"deleting {what}"):
        shutil.rmtree(target)
    ui.ok(f"removed {what} ({len(old_files):,} FASTA files, {len(set(old_files) - set(new_files))} of them superseded)")
    return True



# ═══════════════════════════════════════════════════════════════════════════════

def main(argv=None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        return 1
    ui.init(plain=args.plain or bool(os.environ.get("NO_COLOR")))
    ui.banner(args.command, C.VERSION)
    try:
        ctx = Context(args)
        t0 = time.time()
        rc = {"setup": cmd_setup, "fetch-db": cmd_fetch_db, "check": cmd_check, "run": cmd_run,
              "sketch-db": cmd_sketch_db, "release": cmd_release}[args.command](ctx)
        ui.console().print(f"\n[faint]total time {ui.fmt_duration(time.time() - t0)}[/]")
        return rc
    except C.CommandError as e:
        hint = ""
        if getattr(args, "command", "") == "run":
            hint = "Fix the problem and re-run the same command: finished stages are reused automatically."
        ui.error_panel("A command failed", f"exit {e.returncode}: {C.cmd_str(e.cmd)}", e.tail, hint)
        return 2
    except (RuntimeError, SystemExit) as e:
        if isinstance(e, SystemExit) and (e.code in (0, None) or isinstance(e.code, int) and not str(e)):
            raise
        ui.error_panel("Stopped", str(e))
        return 2
    except KeyboardInterrupt:
        ui.warn("interrupted — re-run the same command to resume")
        return 130


if __name__ == "__main__":
    sys.exit(main())
