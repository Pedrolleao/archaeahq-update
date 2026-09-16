"""
proteins.py — predicted proteins of a database version.

ArchaeaHQ v1.0 proteins (`faa.zip` on figshare) were predicted with Prodigal 2.6.3 in metagenome
mode on the selected genomes, so `release` does the same for every genome of the new version that
has no protein file yet. CheckM2's own Prodigal proteins are not reused: CheckM2 runs Prodigal in
single mode (and may pick translation table 4), which gives a slightly different gene set.

Each genome is written to `<name>.faa.part` and renamed when Prodigal finishes, so an interrupted
release resumes where it stopped.
"""

from __future__ import annotations

import gzip
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence

from .common import log, run_cmd, strip_fna

PRODIGAL_MODE = "meta"


def strip_faa(name: str) -> str:
    for ext in (".faa.gz", ".faa"):
        if name.endswith(ext):
            return name[: -len(ext)]
    return name


def faa_index(folder: Path) -> Dict[str, Path]:
    """genome name → protein FASTA, for every .faa file directly in `folder` (empty when it does not exist)."""
    idx: Dict[str, Path] = {}
    if folder.is_dir():
        for p in sorted(folder.iterdir()):
            if p.is_file() and (p.name.endswith(".faa") or p.name.endswith(".faa.gz")):
                idx.setdefault(strip_faa(p.name), p)
    return idx


def predict_one(tools, fna: Path, out_dir: Path, log_path: Optional[Path]) -> Path:
    """Run `prodigal -p meta` on one genome; returns <out_dir>/<name>.faa (reused when it exists)."""
    name = strip_fna(fna.name)
    out = out_dir / f"{name}.faa"
    if out.exists():
        return out
    part = out.with_name(out.name + ".part")
    tmp = Path(tempfile.mkdtemp(prefix="prodigal_", dir=str(out_dir)))
    try:
        src = fna
        if fna.name.endswith(".gz"):                       # prodigal reads plain FASTA only
            src = tmp / f"{name}.fna"
            with gzip.open(fna, "rb") as fi, open(src, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        run_cmd(tools["prodigal"]("-p", PRODIGAL_MODE, "-q", "-i", src.resolve(), "-a", part.resolve(),
                                  "-o", (tmp / "genes.out").resolve()), log_path=log_path)
        part.replace(out)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    return out


def predict_all(tools, fnas: Sequence[Path], out_dir: Path, threads: int, log_path: Optional[Path],
                progress_cb: Optional[Callable[[int], None]] = None) -> Dict[str, Path]:
    """Name → protein FASTA for every genome Prodigal succeeded on; failures are logged and left out."""
    out_dir.mkdir(parents=True, exist_ok=True)
    done: Dict[str, Path] = {}
    workers = max(1, min(threads, len(fnas))) if fnas else 1
    with ThreadPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(predict_one, tools, f, out_dir, log_path): f for f in fnas}
        for fut in as_completed(futs):
            f = futs[fut]
            try:
                done[strip_fna(f.name)] = fut.result()
            except Exception as e:  # keep going; the genome is reported without proteins
                log.warning("Prodigal failed for %s: %s", f.name, e)
            if progress_cb:
                progress_cb(1)
    return done
