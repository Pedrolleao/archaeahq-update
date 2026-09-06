"""
ncbi.py — NCBI Datasets CLI wrappers: list assemblies per kingdom, summarise a given
accession list, and download genome FASTA files.
"""

from __future__ import annotations

import csv
import shutil
import tempfile
import zipfile
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

from common import CommandError, log, run_cmd, run_pipe, write_lines

# Fields requested from `dataformat tsv genome` and the keys we store them under.
SUMMARY_FIELDS = [
    ("accession", "accession"),
    ("assminfo-name", "asm_name"),
    ("assminfo-paired-assm-accession", "paired"),
    ("assminfo-paired-assm-status", "paired_status"),
    ("assminfo-status", "status"),
    ("assminfo-suppression-reason", "suppression_reason"),
    ("organism-name", "organism"),
    ("organism-tax-id", "taxid"),
    ("assminfo-release-date", "release_date"),
    ("assminfo-level", "level"),
    ("assminfo-bioproject", "bioproject"),
    ("assminfo-biosample-accession", "biosample"),
    ("assminfo-refseq-category", "refseq_category"),
    ("assmstats-total-sequence-len", "total_len"),
]
SUMMARY_KEYS = [k for _, k in SUMMARY_FIELDS]


def _parse_summary_tsv(path: Path) -> List[dict]:
    rows: List[dict] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            parts += [""] * (len(SUMMARY_KEYS) - len(parts))
            rows.append(dict(zip(SUMMARY_KEYS, parts)))
    return rows


def _summary_cmd(tools, kind: str, target: str, api_key: Optional[str]) -> List[str]:
    cmd = tools["datasets"]("summary", "genome", kind, target, "--as-json-lines")
    if kind == "taxon":
        cmd += ["--assembly-source", "all"]
    if api_key:
        cmd += ["--api-key", api_key]
    return cmd


def _dataformat_cmd(tools) -> List[str]:
    fields = ",".join(f for f, _ in SUMMARY_FIELDS)
    return tools["dataformat"]("tsv", "genome", "--fields", fields, "--elide-header")


def list_taxon(tools, taxid: int, out_tsv: Path, api_key: Optional[str], log_path: Path) -> List[dict]:
    """All assemblies (GenBank + RefSeq) below a taxon → rows (also written to out_tsv)."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    run_pipe(_summary_cmd(tools, "taxon", str(taxid), api_key), _dataformat_cmd(tools), out_tsv, log_path)
    rows = _parse_summary_tsv(out_tsv)
    # Write a headed copy for humans
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_KEYS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


def summarize_accessions(tools, accessions: Sequence[str], out_tsv: Path, api_key: Optional[str],
                         log_path: Path) -> List[dict]:
    """Assembly summaries for an explicit accession list (used with --accessions-file)."""
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    rows: List[dict] = []
    with tempfile.TemporaryDirectory(prefix="ahq_ids_") as td:
        for i in range(0, len(accessions), 500):
            chunk = list(accessions[i:i + 500])
            ids = Path(td) / f"ids_{i}.txt"
            write_lines(ids, chunk)
            raw = Path(td) / f"raw_{i}.tsv"
            cmd = tools["datasets"]("summary", "genome", "accession", "--inputfile", str(ids), "--as-json-lines")
            if api_key:
                cmd += ["--api-key", api_key]
            run_pipe(cmd, _dataformat_cmd(tools), raw, log_path)
            rows.extend(_parse_summary_tsv(raw))
    with open(out_tsv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=SUMMARY_KEYS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return rows


# ── Download ──────────────────────────────────────────────────────────────────

def validate_fasta(path: Path) -> bool:
    try:
        if path.stat().st_size < 100:
            return False
        with open(path, "rb") as fh:
            head = fh.read(2)
        return head.startswith(b">")
    except OSError:
        return False


def download_genomes(tools, accessions: Sequence[str], out_dir: Path, api_key: Optional[str], log_path: Path,
                     batch_size: int = 400, retries: int = 2,
                     progress_cb: Optional[Callable[[int], None]] = None) -> Dict[str, Optional[Path]]:
    """
    Download genomic FASTA for each accession into out_dir as <acc>_<asm>_genomic.fna.
    Returns {accession: path or None}. Already-present valid files are not re-downloaded.
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    result: Dict[str, Optional[Path]] = {}
    existing = {p.name.split("_")[0] + "_" + p.name.split("_")[1]: p for p in out_dir.glob("GC*_genomic.fna")
                if validate_fasta(p)}
    pending: List[str] = []
    for acc in accessions:
        if acc in existing:
            result[acc] = existing[acc]
            if progress_cb:
                progress_cb(1)
        else:
            pending.append(acc)

    work = out_dir / "_download_tmp"
    work.mkdir(exist_ok=True)
    attempt = 0
    while pending and attempt <= retries:
        attempt += 1
        still: List[str] = []
        for i in range(0, len(pending), batch_size):
            batch = pending[i:i + batch_size]
            got = _download_batch(tools, batch, out_dir, work, api_key, log_path, f"a{attempt}_b{i // batch_size}")
            for acc in batch:
                p = got.get(acc)
                if p and validate_fasta(p):
                    result[acc] = p
                    if progress_cb:
                        progress_cb(1)
                else:
                    if p:
                        p.unlink(missing_ok=True)
                    still.append(acc)
        pending = still
        if pending:
            log.warning("download attempt %d: %d accession(s) still missing", attempt, len(pending))
    for acc in pending:
        result[acc] = None
        if progress_cb:
            progress_cb(1)
    shutil.rmtree(work, ignore_errors=True)
    return result


def _download_batch(tools, batch: Sequence[str], out_dir: Path, work: Path, api_key: Optional[str],
                    log_path: Path, tag: str) -> Dict[str, Path]:
    ids = work / f"{tag}.txt"
    write_lines(ids, batch)
    zip_path = work / f"{tag}.zip"
    cmd = tools["datasets"]("download", "genome", "accession", "--inputfile", str(ids), "--include", "genome",
                            "--filename", str(zip_path), "--no-progressbar")
    if api_key:
        cmd += ["--api-key", api_key]
    got: Dict[str, Path] = {}
    try:
        run_cmd(cmd, log_path=log_path)
    except CommandError as e:
        log.warning("datasets download failed for batch %s: %s", tag, e)
        if not zip_path.exists():
            return got
    extract = work / tag
    try:
        with zipfile.ZipFile(zip_path) as zf:
            for member in zf.namelist():
                if member.endswith("_genomic.fna") and "/data/" in member:
                    zf.extract(member, extract)
    except zipfile.BadZipFile:
        log.warning("bad zip for batch %s", tag)
        return got
    for fna in extract.glob("ncbi_dataset/data/*/*_genomic.fna"):
        acc = fna.parent.name
        dest = out_dir / fna.name
        shutil.move(str(fna), dest)
        got[acc] = dest
    shutil.rmtree(extract, ignore_errors=True)
    zip_path.unlink(missing_ok=True)
    return got
