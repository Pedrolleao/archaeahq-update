"""
redundancy.py — skani-based redundancy checks reproducing the dRep post-processing rules of
ArchaeaHQ v1.0 (analyze_drep_clusters.py):

  identity : ANI >= 99.9 and aligned fraction >= 99   (+ GCA/GCF same-numeric-accession rule)
  species  : ANI >= 95   and aligned fraction >= 95

Three comparisons are made:  new vs database (skani search against a cached sketch of the DB),
new vs new (skani triangle), and the accession-twin rule.  Within the batch of new genomes the
representative of each species cluster is the genome with the highest Park quality score
(the dRep scoring used for the database), ties going to the newest RefSeq version.
"""

from __future__ import annotations

import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .common import (acc_is_refseq, acc_numeric, acc_parts, acc_version, log, name_to_acc, run_cmd, sha1_of,
                    strip_fna, write_lines)


class UnionFind:
    def __init__(self):
        self.parent: Dict[str, str] = {}

    def add(self, x: str) -> None:
        self.parent.setdefault(x, x)

    def find(self, x: str) -> str:
        self.add(x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[rb] = ra

    def components(self) -> Dict[str, set]:
        groups: Dict[str, set] = defaultdict(set)
        for x in self.parent:
            groups[self.find(x)].add(x)
        return dict(groups)


# ── skani wrappers ────────────────────────────────────────────────────────────

def db_fna_files(db_fna_dir: Path) -> List[Path]:
    files = sorted(p for p in db_fna_dir.iterdir()
                   if p.is_file() and p.suffix in (".fna", ".fa", ".fasta") or p.name.endswith(".fna.gz"))
    return files


def ensure_db_sketch(tools, db_fna_dir: Optional[Path], db_sketch: Optional[Path], cache_dir: Path,
                     threads: int, log_path: Path, ui=None) -> Tuple[Path, int]:
    """Return (sketch_dir, n_genomes). Builds and caches the sketch of the DB fna folder when needed."""
    if db_sketch:
        if not (db_sketch / "sketches.db").exists() and not (db_sketch / "markers.bin").exists():
            raise RuntimeError(f"{db_sketch} does not look like a `skani sketch` output folder")
        n = _sketch_count(db_sketch)
        return db_sketch, n
    if not db_fna_dir:
        raise RuntimeError("either --db-fna or --db-sketch is required")
    files = db_fna_files(db_fna_dir)
    if not files:
        raise RuntimeError(f"no genome FASTA files found in {db_fna_dir}")
    key = sha1_of(f"{p.name}\t{p.stat().st_size}" for p in files)[:12]
    sketch_dir = cache_dir / f"db_sketch_{key}"
    if (sketch_dir / "sketches.db").exists() and (sketch_dir / ".complete").exists():
        return sketch_dir, len(files)
    if ui:
        ui.info(f"sketching {len(files):,} database genomes with skani (cached for later runs)")
    cache_dir.mkdir(parents=True, exist_ok=True)
    list_file = cache_dir / f"db_list_{key}.txt"
    write_lines(list_file, (str(p) for p in files))
    if sketch_dir.exists():
        import shutil
        shutil.rmtree(sketch_dir)
    run_cmd(tools["skani"]("sketch", "-l", list_file, "-o", sketch_dir, "-t", threads), log_path=log_path)
    (sketch_dir / ".complete").write_text(f"{len(files)} genomes\n")
    return sketch_dir, len(files)


def _sketch_count(sketch_dir: Path) -> int:
    c = sketch_dir / ".complete"
    if c.exists():
        try:
            return int(c.read_text().split()[0])
        except (ValueError, IndexError):
            pass
    return -1


def search_db(tools, sketch_dir: Path, queries: Sequence[Path], out_tsv: Path, threads: int, log_path: Path) -> None:
    ql = out_tsv.with_suffix(".queries.txt")
    write_lines(ql, (str(q) for q in queries))
    run_cmd(tools["skani"]("search", "-d", sketch_dir, "--ql", ql, "--min-af", 0, "-t", threads, "-o", out_tsv),
            log_path=log_path)


def triangle(tools, queries: Sequence[Path], out_tsv: Path, threads: int, log_path: Path) -> None:
    ql = out_tsv.with_suffix(".list.txt")
    write_lines(ql, (str(q) for q in queries))
    if len(queries) < 2:
        out_tsv.write_text("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\n")
        return
    run_cmd(tools["skani"]("triangle", "-l", ql, "-E", "--min-af", 0, "-t", threads, "-o", out_tsv),
            log_path=log_path)


def parse_skani(path: Path) -> List[Tuple[str, str, float, float, float]]:
    """(ref_name, query_name, ANI, AF_ref, AF_query) with names = file basenames without extension."""
    hits = []
    if not path.exists():
        return hits
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            try:
                ani, afr, afq = float(parts[2]), float(parts[3]), float(parts[4])
            except ValueError:
                continue
            hits.append((strip_fna(os.path.basename(parts[0])), strip_fna(os.path.basename(parts[1])), ani, afr, afq))
    return hits


# ── Evaluation ────────────────────────────────────────────────────────────────

def evaluate(new: Dict[str, dict], db_index: Dict[str, dict], hits_db, hits_new, th: dict) -> Dict[str, dict]:
    """
    new      : {Name: {"acc":..., "park": float|None, "eligible": bool}}  (eligible = passed quality gate)
    db_index : {Name: {"acc":..., "park": float|None}} for the database genomes
    Returns {Name: {"identical_db": partner|None, "species_db": partner|None, "species_db_note": str,
                    "batch_rep": bool|None, "batch_partner": str|None, "batch_cluster": [..]}}
    """
    id_ani, id_af = th["identity"]["min_ani"], th["identity"]["min_aligned_fraction"]
    sp_ani, sp_af = th["species"]["min_ani"], th["species"]["min_aligned_fraction"]

    db_by_numeric: Dict[str, str] = {}
    for name, rec in db_index.items():
        db_by_numeric.setdefault(acc_numeric(rec["acc"]), name)

    res = {name: {"identical_db": None, "identical_db_ani": None, "species_db": None, "species_db_ani": None,
                  "species_db_note": "", "batch_rep": None, "batch_partner": None, "batch_cluster": []}
           for name in new}

    # 1. new vs database
    best_species: Dict[str, Tuple[float, float, str]] = {}
    best_identity: Dict[str, Tuple[float, float, str]] = {}
    for ref, query, ani, afr, afq in hits_db:
        if query not in res:
            continue
        af = min(afr, afq)
        if ani >= id_ani and af >= id_af:
            if query not in best_identity or (ani, af) > best_identity[query][:2]:
                best_identity[query] = (ani, af, ref)
        if ani >= sp_ani and af >= sp_af:
            if query not in best_species or (ani, af) > best_species[query][:2]:
                best_species[query] = (ani, af, ref)
    for name, rec in new.items():
        twin = db_by_numeric.get(acc_numeric(rec["acc"]))
        if twin:
            res[name]["identical_db"] = twin
            db_acc = db_index[twin]["acc"]
            a, d = acc_parts(rec["acc"]), acc_parts(db_acc)
            if db_acc == rec["acc"]:
                res[name]["identical_db_ani"] = "same accession"
            elif a and d and a[0] == d[0] and a[2] > d[2]:
                res[name]["identical_db_ani"] = f"newer version of {db_acc}"
            else:
                res[name]["identical_db_ani"] = "GCA/GCF twin"
        elif name in best_identity:
            ani, af, ref = best_identity[name]
            res[name]["identical_db"] = ref
            res[name]["identical_db_ani"] = f"ANI {ani:.2f} AF {af:.1f}"
        if name in best_species:
            ani, af, ref = best_species[name]
            res[name]["species_db"] = ref
            res[name]["species_db_ani"] = f"ANI {ani:.2f} AF {af:.1f}"
            p_new, p_db = rec.get("park"), db_index.get(ref, {}).get("park")
            if res[name]["identical_db"] is None and p_new is not None and p_db is not None and p_new > p_db:
                res[name]["species_db_note"] = (f"higher quality score than DB genome {name_to_acc(ref)} "
                                                f"({p_new:.2f} vs {p_db:.2f})")
        elif twin and name not in best_species:
            res[name]["species_db"] = twin
            res[name]["species_db_ani"] = res[name]["identical_db_ani"]

    # 2. new vs new: species-level clusters among all compared genomes.  The representative of a
    #    cluster is chosen among the members that passed the quality gate and are not redundant with
    #    the database (highest Park score, then newest RefSeq); every other member fails the
    #    within-batch filter and points to the representative.
    uf = UnionFind()
    for n in new:
        uf.add(n)
    for ref, query, ani, afr, afq in hits_new:
        if ref in new and query in new and ref != query:
            if ani >= sp_ani and min(afr, afq) >= sp_af:
                uf.union(ref, query)
    by_numeric: Dict[str, List[str]] = defaultdict(list)
    for n in new:
        by_numeric[acc_numeric(new[n]["acc"])].append(n)
    for names in by_numeric.values():
        for other in names[1:]:
            uf.union(names[0], other)

    def rank(n: str):
        r = new[n]
        redundant_db = res[n]["identical_db"] is not None or res[n]["species_db"] is not None
        park = r.get("park")
        return (not r.get("eligible"), redundant_db, -(park if park is not None else -1e9),
                not acc_is_refseq(r["acc"]), -acc_version(r["acc"]), n)

    for root, members in uf.components().items():
        rep = min(members, key=rank)
        for n in members:
            res[n]["batch_rep"] = (n == rep)
            res[n]["batch_partner"] = None if n == rep else rep
            res[n]["batch_cluster"] = sorted(members)
    return res
