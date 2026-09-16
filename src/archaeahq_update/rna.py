"""
rna.py — rRNA and tRNA gene counts per genome.

barrnap >= 1.0 (bioconda `barrnap>=1.10`, the same lineage used for ArchaeaHQ v1.0) predicts
rRNA with Infernal/Rfam and tRNA with ARAGORN in one run (`--kingdom arc`).  If only the old
barrnap 0.9 is available, rRNA comes from barrnap and tRNA from a direct ARAGORN call.
The GFF files are parsed with the bundled compile_barrnap.py so the counting is identical to the
original build.

Every barrnap process runs in its own temporary working directory: barrnap 1.x writes fixed-name
scratch files (e.g. barrnap.find_operon.gff) into the current directory, which breaks parallel runs.

The 16S rRNA sequences of a database version (built by `release`) are cut from the genomes at the
GFF coordinates (reverse-complemented on the minus strand) with the header layout of the v1.0
`Archaea_HQ-16S.fasta`: `>Name#scaffold#start-stop#strand#16S_rRNA`.
"""

from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Collection, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

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


# ─── rRNA sequences ──────────────────────────────────────────────────────────

_COMPLEMENT = str.maketrans("ACGTacgtNnRrYyKkMmSsWwBbDdHhVv",
                            "TGCAtgcaNnYyRrMmKkSsWwVvHhDdBb")


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def _rrna_name(info: str) -> str:
    for part in info.split(";"):
        if part.startswith("Name="):
            return part[5:].strip()
    return ""


def _read_scaffolds(fna: Path, wanted: Collection[str]) -> Dict[str, str]:
    """Sequences of the wanted scaffolds only (ID = first header token); stops once all are read."""
    wanted = set(wanted)
    seqs: Dict[str, str] = {}
    cur: Optional[str] = None
    buf: List[str] = []
    opener = gzip.open if fna.name.endswith(".gz") else open
    with opener(fna, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                    if len(seqs) == len(wanted):
                        return seqs
                tok = line[1:].split()
                cur = tok[0] if tok and tok[0] in wanted else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def iter_rrna(fnas: Sequence[Path], gff_dir: Path, rrna: str = "16S_rRNA") -> Iterator[Tuple[str, str, str]]:
    """
    (genome Name, sequence ID, sequence) for every `rrna` gene (barrnap Name=, e.g. 16S_rRNA) in the
    GFFs of `fnas`, genomes in name order. Genomes without a GFF or without a hit yield nothing.
    """
    for fna in sorted(fnas, key=lambda p: p.name):
        name = strip_fna(fna.name)
        gff = _gff_path(gff_dir, fna)
        if not gff.exists():
            continue
        hits = [r for r in cb.parse_gff_file(str(gff))
                if r["rna_type"] == "rRNA" and _rrna_name(r["info"]) == rrna]
        if not hits:
            continue
        scaffolds = _read_scaffolds(fna, {r["scaffold"] for r in hits})
        for r in hits:
            if r["scaffold"] not in scaffolds:
                log.warning("%s: scaffold %s of a %s hit not found in %s", name, r["scaffold"], rrna, fna.name)
                continue
            start, stop = int(r["start"]), int(r["stop"])
            seq = scaffolds[r["scaffold"]][start - 1:stop]          # GFF: 1-based, inclusive
            if r["strand"] == "-":
                seq = reverse_complement(seq)
            yield name, f"{name}#{r['scaffold']}#{start}-{stop}#{r['strand']}#{rrna}", seq


def read_fasta(path: Path) -> Iterator[Tuple[str, str]]:
    """(header without '>', sequence) for every record of a FASTA file."""
    header: Optional[str] = None
    buf: List[str] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(buf)
                header, buf = line[1:].strip(), []
            elif header is not None:
                buf.append(line.strip())
    if header is not None:
        yield header, "".join(buf)


def genome_of_seq_id(seq_id: str) -> str:
    """Genome Name of a `Name#scaffold#start-stop#strand#16S_rRNA` sequence ID."""
    return strip_fna(seq_id.split("#", 1)[0])


RRNA_TABLE_COLUMNS = ["Genome", "Sequence_ID", "Length_bp"]


def write_rrna(records: Iterable[Tuple[str, str, str]], fasta_path: Path, table_path: Path) -> Dict[str, int]:
    """
    Write (Name, sequence ID, sequence) records to a FASTA and to a table with the layout of
    Supplementary Table 5 (Genome, Sequence_ID, Length_bp). Returns Name → sequences written.
    """
    counts: Dict[str, int] = {}
    fasta_tmp = fasta_path.with_name(fasta_path.name + ".part")
    table_tmp = table_path.with_name(table_path.name + ".part")
    fasta_path.parent.mkdir(parents=True, exist_ok=True)
    with open(fasta_tmp, "w", encoding="utf-8") as fa, open(table_tmp, "w", encoding="utf-8") as tb:
        tb.write("\t".join(RRNA_TABLE_COLUMNS) + "\n")
        for name, seq_id, seq in records:
            fa.write(f">{seq_id}\n{seq}\n")
            tb.write(f"{name}\t{seq_id}\t{len(seq)}\n")
            counts[name] = counts.get(name, 0) + 1
    fasta_tmp.replace(fasta_path)
    table_tmp.replace(table_path)
    return counts
