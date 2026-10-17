"""
quality.py — CheckM2 wrapper and the quality columns of the ArchaeaHQ table.
"""

from __future__ import annotations

import math
import re
import shutil
from pathlib import Path
from typing import Callable, Dict, Optional

from .common import NOT_AVAILABLE, fmt_num, park_score, read_tsv, run_cmd, run_cmd_live, to_float

# genomes per DIAMOND batch inside CheckM2 (DefaultValues.DIAMOND_DEFAULT_CHUNK_SIZE)
CHECKM2_DIAMOND_CHUNK = 250

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


class CheckM2Progress:
    """
    Turn the log lines of `checkm2 predict` into progress: (phase, done, total). CheckM2 announces each
    phase and counts genomes through gene calling and metadata; DIAMOND itself is silent, so its progress
    is the number of DIAMOND_RESULTS files that have appeared in diamond_output/ (one per batch of
    CHECKM2_DIAMOND_CHUNK genomes).
    """
    PHASES = [
        (re.compile(r"Calling genes in (\d+) bins"), "calling genes (prodigal)"),
        (re.compile(r"Calculating metadata for (\d+) bins"), "genome metadata"),
        (re.compile(r"Annotating input genomes with DIAMOND"), "DIAMOND annotation"),
        (re.compile(r"Processing DIAMOND output"), "processing DIAMOND output"),
        (re.compile(r"Predicting completeness and contamination"), "predicting completeness/contamination"),
        (re.compile(r"Parsing all results"), "writing the report"),
        (re.compile(r"CheckM2 finished successfully"), "finished"),
    ]
    COUNTER = re.compile(r"Finished processing (\d+) of (\d+)")
    N_PHASES = 6  # phases a user sees before "finished"

    def __init__(self, out_dir: Path, callback: Callable[[str, Optional[int], Optional[int]], None]):
        self.out_dir = out_dir
        self.cb = callback
        self.n_genomes: Optional[int] = None
        self.phase_no = 0
        self.phase = "starting CheckM2"
        self.done: Optional[int] = None
        self.total: Optional[int] = None
        self.cb(self.label(), None, None)

    def label(self) -> str:
        if not self.phase_no:
            return self.phase
        if self.phase_no > self.N_PHASES:
            return f"CheckM2 {self.phase}"
        return f"CheckM2 {self.phase_no}/{self.N_PHASES} {self.phase}"

    def _emit(self):
        self.cb(self.label(), self.done, self.total)

    def on_line(self, line: str) -> None:
        # phase announcements first: CheckM2 logs them on stderr and the counters on stdout, so the two
        # can arrive glued into one line
        if self._phase_line(line):
            return
        m = self.COUNTER.search(line)
        if m and self.phase_no in (1, 2):
            self.done, self.total = int(m.group(1)), int(m.group(2))
            self._emit()

    def _phase_line(self, line: str) -> bool:
        for i, (rx, name) in enumerate(self.PHASES, start=1):
            m = rx.search(line)
            if not m:
                continue
            self.phase_no, self.phase = i, name
            n = int(m.group(1)) if m.groups() else None
            if n is not None:
                self.n_genomes = n
            if i in (1, 2):
                self.done, self.total = 0, n or self.n_genomes
            elif i == 3:
                self.total = max(1, math.ceil((self.n_genomes or 1) / CHECKM2_DIAMOND_CHUNK))
                self.done = 0
                if self.total == 1:
                    self.phase += " (one batch; DIAMOND reports no finer progress)"
            else:
                self.done, self.total = None, None
            self._emit()
            return True
        return False

    def poll(self) -> None:
        """Called periodically while CheckM2 runs: count DIAMOND batches started."""
        if self.phase_no != 3:
            return
        d = self.out_dir / "diamond_output"
        if not d.is_dir():
            return
        started = sum(1 for f in d.iterdir() if f.name.startswith("DIAMOND_RESULTS"))
        done = max(0, min(self.total or 1, started - 1))  # the newest file is the batch still running
        if done != self.done:
            self.done = done
            self._emit()


def run_checkm2(tools, fna_dir: Path, out_dir: Path, threads: int, db_path: Path, log_path: Path,
                extension: str = "fna", lowmem: bool = False,
                progress: Optional[Callable[[str, Optional[int], Optional[int]], None]] = None) -> Path:
    """
    Run `checkm2 predict` on every .fna in fna_dir; return the quality_report.tsv path.
    With `progress(label, done, total)` CheckM2's output is streamed and reported phase by phase.
    """
    report = out_dir / "quality_report.tsv"
    if out_dir.exists():
        shutil.rmtree(out_dir)
    cmd = tools["checkm2"]("predict", "--threads", threads, "--input", fna_dir, "-x", extension,
                           "--output-directory", out_dir, "--database_path", db_path, "--force")
    if lowmem:
        cmd.append("--lowmem")
    if progress:
        tracker = CheckM2Progress(out_dir, progress)
        run_cmd_live(cmd, log_path, tracker.on_line, poll=tracker.poll, poll_interval=2.0)
    else:
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
