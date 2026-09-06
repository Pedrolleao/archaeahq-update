"""
quality.py — CheckM2 wrapper and the quality columns of the ArchaeaHQ table.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Optional

from .common import NOT_AVAILABLE, fmt_num, park_score, read_tsv, run_cmd, to_float

# CheckM2 quality_report.tsv column → ArchaeaHQ column
CHECKM2_COLUMNS = {
    "Completeness": "Completeness",
    "Contamination": "Contamination",
    "Coding_Density": "Coding_Density",
    "Contig_N50": "Contig_N50",
    "Average_Gene_Length": "Average_Gene_Length",
    "Genome_Size": "Genome_Size",
    "GC_Content": "GC_Content",
    "Total_Coding_Sequences": "Total_Coding_Sequences",
}


def run_checkm2(tools, fna_dir: Path, out_dir: Path, threads: int, db_path: Path, log_path: Path,
                extension: str = "fna", lowmem: bool = False) -> Path:
    """Run `checkm2 predict` on every .fna in fna_dir; return the quality_report.tsv path."""
    report = out_dir / "quality_report.tsv"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    cmd = tools["checkm2"]("predict", "--threads", threads, "--input", fna_dir, "-x", extension,
                           "--output-directory", out_dir, "--database_path", db_path, "--force")
    if lowmem:
        cmd.append("--lowmem")
    run_cmd(cmd, log_path=log_path)
    if not report.exists():
        raise RuntimeError(f"CheckM2 finished but {report} is missing")
    return report


def parse_quality_report(report: Path) -> Dict[str, dict]:
    """{Name: {ArchaeaHQ column: value, ...}} with Park score added."""
    out: Dict[str, dict] = {}
    for row in read_tsv(report):
        name = row.get("Name", "")
        rec = {dst: row.get(src, NOT_AVAILABLE) for src, dst in CHECKM2_COLUMNS.items()}
        rec["Park Quality_Score"] = fmt_num(park_score(rec["Completeness"], rec["Contamination"], rec["Contig_N50"]), 6)
        rec["_park"] = park_score(rec["Completeness"], rec["Contamination"], rec["Contig_N50"])
        out[name] = rec
    return out


def mimag_label(completeness, contamination, th: dict) -> str:
    c, x = to_float(completeness), to_float(contamination)
    if c is None or x is None:
        return NOT_AVAILABLE
    m = th["mimag"]
    if c >= m["high_min_completeness"] and x <= m["high_max_contamination"]:
        return "High (MIMAG)"
    if c >= m["medium_min_completeness"] and x <= m["medium_max_contamination"]:
        return "Medium (MIMAG)"
    return "Low"


def passes_completeness(completeness, th: dict) -> Optional[bool]:
    c = to_float(completeness)
    return None if c is None else c >= th["quality"]["min_completeness"]


def passes_contamination(contamination, th: dict) -> Optional[bool]:
    x = to_float(contamination)
    return None if x is None else x <= th["quality"]["max_contamination"]
