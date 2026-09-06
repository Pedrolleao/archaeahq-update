"""
rna.py — rRNA and tRNA gene counts per genome.

barrnap >= 1.0 (bioconda `barrnap>=1.10`, the same lineage used for ArchaeaHQ v1.0) predicts
rRNA with Infernal/Rfam and tRNA with ARAGORN in one run (`--kingdom arc`).  If only the old
barrnap 0.9 is available, rRNA comes from barrnap and tRNA from a direct ARAGORN call.
The GFF files are parsed with the bundled compile_barrnap.py so the counting is identical to the
original build.

Every barrnap process runs in its own temporary working directory: barrnap 1.x writes fixed-name
scratch files (e.g. barrnap.find_operon.gff) into the current directory, which breaks parallel runs.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from .common import CommandError, log, run_cmd, strip_fna
from . import compile_barrnap as cb

_flags_cache: Dict[str, List[str]] = {}


def barrnap_flags(tool) -> List[str]:
    """Flags for a barrnap >= 1.0 build, chosen from its own --help (1.7 dev and 1.10 differ)."""
    key = " ".join(tool.cmd)
    if key in _flags_cache:
        return _flags_cache[key]
    try:
        out = subprocess.run(tool.cmd + ["--help"], capture_output=True, text=True, timeout=120)
        help_text = (out.stdout or "") + (out.stderr or "")
    except (subprocess.SubprocessError, OSError):
        help_text = ""
    flags = ["--kingdom", "arc", "--threads", "1", "--quiet", "--rrna", "--trna"]
    for opt in ("--no-ncrna", "--no-mrna", "--no-operon"):
        if opt[5:] in help_text:          # "ncrna", "mrna", "operon" mentioned in help
            flags.append(opt)
    _flags_cache[key] = flags
    return flags


def _gff_path(out_dir: Path, fna: Path) -> Path:
    return out_dir / f"{fna.name}-barrnap.gff"


def run_one(tools, fna: Path, out_dir: Path, new_barrnap: bool, log_path: Optional[Path]) -> Tuple[str, int, int]:
    """Return (Name, rRNA count, tRNA count) for one genome; writes <name>.fna-barrnap.gff."""
    name = strip_fna(fna.name)
    gff = _gff_path(out_dir, fna)
    if not gff.exists() or gff.stat().st_size == 0:
        tmp = Path(tempfile.mkdtemp(prefix="barrnap_", dir=str(out_dir)))
        try:
            if new_barrnap:
                cmd = tools["barrnap"](*barrnap_flags(tools["barrnap"]), fna.resolve())
            else:
                cmd = tools["barrnap"]("--kingdom", "arc", "--threads", 1, "--quiet", fna.resolve())
            run_cmd(cmd, stdout_path=gff, log_path=log_path, cwd=tmp)
        finally:
            shutil.rmtree(tmp, ignore_errors=True)
    records = cb.parse_gff_file(str(gff))
    r = sum(1 for rec in records if rec["rna_type"] == "rRNA")
    t = sum(1 for rec in records if rec["rna_type"] == "tRNA")
    if not new_barrnap:
        t = count_trna_aragorn(tools, fna, out_dir, log_path)
    return name, r, t


_ARAGORN_LINE = re.compile(r"^\s*\d+\s+tRNA-")


def count_trna_aragorn(tools, fna: Path, out_dir: Path, log_path: Optional[Path]) -> int:
    if "aragorn" not in tools:
        return -1
    out = out_dir / f"{fna.name}-aragorn.txt"
    if not out.exists():
        run_cmd(tools["aragorn"]("-t", "-gc11", "-w", fna.resolve()), stdout_path=out, log_path=log_path)
    n = 0
    with open(out, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if _ARAGORN_LINE.match(line):
                n += 1
    return n


def run_all(tools, fnas: Sequence[Path], out_dir: Path, threads: int, new_barrnap: bool, log_path: Optional[Path],
            progress_cb: Optional[Callable[[int], None]] = None) -> Dict[str, Tuple[int, int]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    results: Dict[str, Tuple[int, int]] = {}
    workers = max(1, min(threads, len(fnas))) if fnas else 1
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(run_one, tools, f, out_dir, new_barrnap, log_path): f for f in fnas}
        for fut in as_completed(futs):
            f = futs[fut]
            try:
                name, r, t = fut.result()
                results[name] = (r, t)
            except Exception as e:  # keep going; the genome gets NA counts
                log.warning("RNA annotation failed for %s: %s", f.name, e)
                gff = _gff_path(out_dir, f)
                if gff.exists() and gff.stat().st_size == 0:
                    gff.unlink()
                results[strip_fna(f.name)] = (-1, -1)
            if progress_cb:
                progress_cb(1)
    return results
