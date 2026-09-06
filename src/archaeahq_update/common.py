"""
common.py — shared constants, small helpers and the subprocess wrapper used by every stage.

Nothing here refers to paths outside the bundle: the bundle directory is derived from this file.
"""

from __future__ import annotations

import csv
import hashlib
import json
import logging
import math
import os
import re
import shlex
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, List, Optional, Sequence

from . import __version__ as VERSION  # noqa: E402

PACKAGE_DIR = Path(__file__).resolve().parent
BUNDLE_DIR = PACKAGE_DIR                      # kept for callers that protect the install location
DATA_DIR = PACKAGE_DIR / "data"
DB_TABLE_PATH = DATA_DIR / "ArchaeaHQ-Info.tsv"
EVALUATED_PATH = DATA_DIR / "evaluated_accessions.tsv"
ENVIRONMENT_YML = DATA_DIR / "environment.yml"

log = logging.getLogger("archaeahq")

# Exact column order of the ArchaeaHQ information table (Supplementary Table 1).
DB_COLUMNS: List[str] = [
    "Name", "Assembly ID", "Archaeal_Kingdom", "Taxonomy", "Completeness", "Contamination",
    "Park Quality_Score", "Bioproject", "Biosample", "Environment_Category-ArchaeaHQ",
    "Metagenome_source", "geo_loc_name", "lat_lon", "Isolation_Source", "Env_Broad", "Env_Local",
    "Env_medium", "Coding_Density", "Contig_N50", "Average_Gene_Length", "Genome_Size", "GC_Content",
    "Total_Coding_Sequences", "rRNA", "tRNA",
]

# Value used in the DB table when a metadata field is absent.
MISSING_TEXT = "0"
MISSING_METAGENOME = "N/A"
NOT_AVAILABLE = "NA"

# Filter columns of the simple PASS/FAIL table, in order.
FILTER_COLUMNS: List[str] = [
    "Downloaded", "Completeness_ge70", "Contamination_le10",
    "Not_identical_to_DB", "Not_same_species_as_DB", "Representative_within_batch",
]

STAGES: List[str] = ["env", "list", "download", "metadata", "quality", "redundancy", "rna", "report"]

# Database release folders: <releases dir>/Archaea_HQ-v<major>.<minor>/{fna/, ArchaeaHQ-Info.tsv, release.json}
RELEASE_PREFIX = "Archaea_HQ-v"
RELEASE_TABLE = "ArchaeaHQ-Info.tsv"
RELEASE_META = "release.json"
BUNDLED_VERSION = "v1.0"


# ── Release versions ──────────────────────────────────────────────────────────

def parse_version(label: str):
    """'v1.2' / '1.2' / 'Archaea_HQ-v1.2' → (1, 2); None when not of that form."""
    m = re.search(r"v?(\d+)\.(\d+)$", (label or "").strip())
    return (int(m.group(1)), int(m.group(2))) if m else None


def version_label(v) -> str:
    return f"v{v[0]}.{v[1]}"


def next_version(label: str) -> str:
    """'v1.1' → 'v1.2'."""
    v = parse_version(label)
    if v is None:
        raise ValueError(f"cannot derive the next version from {label!r} — pass --version")
    return version_label((v[0], v[1] + 1))


def release_folder_name(label: str) -> str:
    v = parse_version(label)
    return f"{RELEASE_PREFIX}{v[0]}.{v[1]}" if v else f"{RELEASE_PREFIX}{label.lstrip('v')}"


def is_release_folder(path: Path) -> bool:
    return (path / RELEASE_META).is_file() and (path / RELEASE_TABLE).is_file()


def find_releases(releases_dir: Path) -> List[Path]:
    """Release folders in releases_dir, oldest version first."""
    found = []
    if releases_dir.is_dir():
        for d in releases_dir.iterdir():
            if d.is_dir() and d.name.startswith(RELEASE_PREFIX) and is_release_folder(d):
                v = parse_version(d.name)
                if v:
                    found.append((v, d))
    return [d for _, d in sorted(found)]


def newest_release(releases_dir: Path) -> Optional[Path]:
    rel = find_releases(releases_dir)
    return rel[-1] if rel else None


def release_version_of_table(table_path: Path) -> str:
    """Version label of the database a table belongs to (release.json next to it; the bundled table is v1.0)."""
    rj = table_path.parent / RELEASE_META
    if rj.exists():
        try:
            v = load_json(rj).get("version")
            if v:
                return str(v)
        except (OSError, ValueError):
            pass
    if table_path.resolve() == DB_TABLE_PATH.resolve():
        return BUNDLED_VERSION
    v = parse_version(table_path.parent.name)
    return version_label(v) if v else table_path.parent.name


# ── JSON / TSV helpers ────────────────────────────────────────────────────────

def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def save_json(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2)


def load_kingdoms() -> List[dict]:
    return load_json(DATA_DIR / "kingdoms.json")["kingdoms"]


def kingdom_by(key: str) -> dict:
    """Index kingdoms by 'short', 'taxid' or 'label'."""
    return {str(k[key]): k for k in load_kingdoms()}


def load_thresholds() -> dict:
    return load_json(DATA_DIR / "thresholds.json")


def read_tsv(path: Path, encoding: str = "utf-8") -> List[dict]:
    with open(path, newline="", encoding=encoding, errors="replace") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict], columns: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(columns), delimiter="\t",
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        for ln in lines:
            fh.write(f"{ln}\n")


# ── Accession helpers ─────────────────────────────────────────────────────────

ACC_RE = re.compile(r"^(GC[AF])_(\d+)\.(\d+)$")


def acc_parts(acc: str):
    """'GCA_000123.2' → ('GCA', '000123', 2) or None."""
    m = ACC_RE.match(acc.strip())
    if not m:
        return None
    return m.group(1), m.group(2), int(m.group(3))


def acc_numeric(acc: str) -> str:
    """Version- and prefix-agnostic key: 'GCF_000123.2' → '000123'."""
    p = acc_parts(acc)
    return p[1] if p else acc


def acc_version(acc: str) -> int:
    p = acc_parts(acc)
    return p[2] if p else 0


def acc_is_refseq(acc: str) -> bool:
    return acc.startswith("GCF_")


def name_to_acc(name: str) -> str:
    """'GCA_000007185.1_ASM718v1_genomic' → 'GCA_000007185.1'."""
    parts = name.split("_")
    if len(parts) >= 2 and parts[0] in ("GCA", "GCF"):
        return f"{parts[0]}_{parts[1]}"
    return name


def strip_fna(name: str) -> str:
    for ext in (".fna.gz", ".fna", ".fa", ".fasta"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def safe_asm_name(asm_name: str) -> str:
    """Approximate NCBI's file-name sanitisation of an assembly name."""
    return re.sub(r"[^A-Za-z0-9._-]+", "_", asm_name.strip())


# ── Numeric helpers ───────────────────────────────────────────────────────────

def to_float(x) -> Optional[float]:
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def park_score(completeness, contamination, contig_n50) -> Optional[float]:
    c, x, n = to_float(completeness), to_float(contamination), to_float(contig_n50)
    if c is None or x is None or n is None or n <= 0:
        return None
    return c - 5.0 * x + 0.5 * math.log10(n)


def fmt_num(x, digits: int = 2) -> str:
    if x is None:
        return NOT_AVAILABLE
    if isinstance(x, float):
        if x.is_integer():
            return str(int(x))
        return f"{x:.{digits}f}".rstrip("0").rstrip(".") if digits else str(x)
    return str(x)


def sha1_of(items: Iterable[str]) -> str:
    h = hashlib.sha1()
    for it in items:
        h.update(it.encode("utf-8"))
        h.update(b"\n")
    return h.hexdigest()


# ── Subprocess wrapper ────────────────────────────────────────────────────────

class CommandError(RuntimeError):
    def __init__(self, cmd: Sequence[str], returncode: int, tail: str):
        self.cmd = list(cmd)
        self.returncode = returncode
        self.tail = tail
        super().__init__(f"command failed (exit {returncode}): {cmd_str(cmd)}")


def cmd_str(cmd: Sequence[str]) -> str:
    return " ".join(shlex.quote(str(c)) for c in cmd)


def run_cmd(
    cmd: Sequence[str],
    log_path: Optional[Path] = None,
    cwd: Optional[Path] = None,
    check: bool = True,
    capture: bool = False,
    stdin_path: Optional[Path] = None,
    stdout_path: Optional[Path] = None,
    env: Optional[dict] = None,
    quiet: bool = True,
) -> subprocess.CompletedProcess:
    """
    Run a command. stdout/stderr are appended to `log_path` (if given) unless `capture`
    is requested (then stdout is returned in .stdout) or `stdout_path` is given (stdout
    goes to that file, stderr to the log).
    """
    cmd = [str(c) for c in cmd]
    log.info("RUN  %s", cmd_str(cmd))
    t0 = time.time()
    lf = open(log_path, "a", encoding="utf-8") if log_path else None
    if lf:
        lf.write(f"\n$ {cmd_str(cmd)}\n")
        lf.flush()
    stdin = open(stdin_path, "rb") if stdin_path else None
    stdout_f = open(stdout_path, "wb") if stdout_path else None
    try:
        if capture:
            proc = subprocess.run(cmd, cwd=cwd, env=env, stdin=stdin, stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE, text=True)
            if lf:
                if proc.stderr:
                    lf.write(proc.stderr)
                lf.write(f"[exit {proc.returncode}, {time.time() - t0:.1f}s]\n")
        else:
            out = stdout_f if stdout_f else (lf if lf else (subprocess.DEVNULL if quiet else None))
            err = lf if lf else (subprocess.DEVNULL if quiet else None)
            proc = subprocess.run(cmd, cwd=cwd, env=env, stdin=stdin, stdout=out, stderr=err, text=True)
            if lf:
                lf.write(f"[exit {proc.returncode}, {time.time() - t0:.1f}s]\n")
    finally:
        if stdin:
            stdin.close()
        if stdout_f:
            stdout_f.close()
        if lf:
            lf.close()
    if check and proc.returncode != 0:
        tail = ""
        if capture and proc.stderr:
            tail = proc.stderr[-3000:]
        elif log_path and log_path.exists():
            tail = tail_of_file(log_path, 25)
        raise CommandError(cmd, proc.returncode, tail)
    return proc


def run_pipe(
    cmd1: Sequence[str],
    cmd2: Sequence[str],
    stdout_path: Path,
    log_path: Optional[Path] = None,
) -> None:
    """cmd1 | cmd2 > stdout_path, with stderr of both appended to the log."""
    cmd1 = [str(c) for c in cmd1]
    cmd2 = [str(c) for c in cmd2]
    log.info("RUN  %s | %s > %s", cmd_str(cmd1), cmd_str(cmd2), stdout_path)
    lf = open(log_path, "a", encoding="utf-8") if log_path else None
    if lf:
        lf.write(f"\n$ {cmd_str(cmd1)} | {cmd_str(cmd2)} > {stdout_path}\n")
        lf.flush()
    err = lf if lf else subprocess.DEVNULL
    with open(stdout_path, "wb") as out:
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=err)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=out, stderr=err)
        p1.stdout.close()
        rc2 = p2.wait()
        rc1 = p1.wait()
    if lf:
        lf.write(f"[exit {rc1} | {rc2}]\n")
        lf.close()
    if rc1 != 0:
        raise CommandError(cmd1, rc1, tail_of_file(log_path, 25) if log_path else "")
    if rc2 != 0:
        raise CommandError(cmd2, rc2, tail_of_file(log_path, 25) if log_path else "")


def tail_of_file(path: Path, n: int = 25) -> str:
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            lines = fh.readlines()
        return "".join(lines[-n:])
    except OSError:
        return ""


# ── Stage bookkeeping ─────────────────────────────────────────────────────────

def marker_path(stage_dir: Path) -> Path:
    return stage_dir / ".done"


def stage_is_done(stage_dir: Path) -> bool:
    return marker_path(stage_dir).exists()


def mark_stage_done(stage_dir: Path, info: Optional[dict] = None) -> None:
    stage_dir.mkdir(parents=True, exist_ok=True)
    payload = {"finished": time.strftime("%Y-%m-%d %H:%M:%S")}
    if info:
        payload.update(info)
    save_json(marker_path(stage_dir), payload)


def clear_stage(stage_dir: Path) -> None:
    m = marker_path(stage_dir)
    if m.exists():
        m.unlink()


def setup_file_logging(log_file: Path, level=logging.INFO) -> None:
    log_file.parent.mkdir(parents=True, exist_ok=True)
    log.setLevel(level)
    for h in list(log.handlers):
        log.removeHandler(h)
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    log.addHandler(fh)

