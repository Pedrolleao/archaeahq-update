"""
fetchdb.py — download ArchaeaHQ v1.0 from figshare and unpack it as the first database version.

Everything is done with the standard library: figshare's public API for the file list and
checksums, resumable HTTP downloads (Range requests, `.part` files), MD5 verification and a
streaming unzip of the genome FASTA files.
"""

from __future__ import annotations

import hashlib
import json
import shutil
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path
from typing import Callable, Dict, List, Optional

FIGSHARE_ARTICLE = 32266599
FIGSHARE_API = f"https://api.figshare.com/v2/articles/{FIGSHARE_ARTICLE}"
FIGSHARE_PAGE = "https://figshare.com/articles/dataset/ArchaeaHQ_v_1_00_2026/32266599"
FIGSHARE_DOI = "10.6084/m9.figshare.32266599"
FIGSHARE_VERSION = "v1.0"

# figshare file name → (key used on the command line, what it is, unpacked size in GB)
FILES = {
    "fna.zip": ("fna", "genome FASTA files (the database)", 36.5),
    "faa.zip": ("faa", "predicted proteins", 15.0),
    "Archaea_HQ-16S.fasta": ("16s", "16S rRNA sequences", 0.02),
    "ArchaeaHQ-Supp_tables.xlsx": ("tables", "supplementary tables of the paper", 0.005),
}

CHUNK = 1 << 20  # 1 MiB
ProgressCb = Optional[Callable[[int], None]]


def article_files(timeout: int = 60) -> Dict[str, dict]:
    """name → {size, download_url, md5} for the files of the figshare article."""
    req = urllib.request.Request(FIGSHARE_API, headers={"User-Agent": "archaeahq-update"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        article = json.load(r)
    out = {}
    for f in article.get("files", []):
        out[f["name"]] = {"size": int(f["size"]), "download_url": f["download_url"],
                          "md5": f.get("computed_md5") or f.get("supplied_md5") or ""}
    return out


def download(url: str, dest: Path, size: int, progress_cb: ProgressCb = None, retries: int = 5,
             timeout: int = 120) -> Path:
    """Download url to dest, resuming a previous `.part` file. progress_cb receives bytes written."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = dest.with_name(dest.name + ".part")
    if dest.exists() and dest.stat().st_size == size:
        if progress_cb:
            progress_cb(size)
        return dest
    attempt = 0
    while True:
        have = part.stat().st_size if part.exists() else 0
        if have > size:                       # corrupt leftover
            part.unlink()
            have = 0
        if have == size:
            break
        if progress_cb and have and attempt == 0:
            progress_cb(have)
        headers = {"User-Agent": "archaeahq-update"}
        if have:
            headers["Range"] = f"bytes={have}-"
        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=timeout) as r:
                if have and r.status != 206:  # server ignored the range: start over
                    have = 0
                    if progress_cb:
                        progress_cb(-part.stat().st_size)
                    part.unlink(missing_ok=True)
                with open(part, "ab" if have else "wb") as fh:
                    while True:
                        chunk = r.read(CHUNK)
                        if not chunk:
                            break
                        fh.write(chunk)
                        if progress_cb:
                            progress_cb(len(chunk))
        except (urllib.error.URLError, ConnectionError, TimeoutError, OSError) as e:
            attempt += 1
            if attempt > retries:
                raise RuntimeError(f"download failed after {retries} retries: {e}") from e
            time.sleep(min(60, 5 * attempt))
            continue
        if part.stat().st_size == size:
            break
        attempt += 1
        if attempt > retries:
            raise RuntimeError(f"download incomplete: {part.stat().st_size:,} of {size:,} bytes")
    part.replace(dest)
    return dest


def md5_of(path: Path, progress_cb: ProgressCb = None) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(CHUNK * 8)
            if not b:
                break
            h.update(b)
            if progress_cb:
                progress_cb(len(b))
    return h.hexdigest()


def zip_members(zip_path: Path, suffixes=(".fna", ".faa", ".fa", ".fasta")) -> List[zipfile.ZipInfo]:
    """Regular files with one of the suffixes, wherever they sit inside the archive."""
    with zipfile.ZipFile(zip_path) as z:
        return [i for i in z.infolist() if not i.is_dir() and i.filename.lower().endswith(suffixes)
                and "@eadir" not in i.filename.lower() and not Path(i.filename).name.startswith(".")]


def extract_flat(zip_path: Path, dest_dir: Path, members: List[zipfile.ZipInfo], progress_cb: ProgressCb = None) -> int:
    """Extract the given members into dest_dir without their folder path; files already complete are skipped."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    n = 0
    with zipfile.ZipFile(zip_path) as z:
        for info in members:
            target = dest_dir / Path(info.filename).name
            if target.exists() and target.stat().st_size == info.file_size:
                if progress_cb:
                    progress_cb(1)
                continue
            with z.open(info) as src, open(target, "wb") as dst:
                shutil.copyfileobj(src, dst, CHUNK)
            n += 1
            if progress_cb:
                progress_cb(1)
    return n


def free_space_gb(path: Path) -> float:
    p = path
    while not p.exists():
        p = p.parent
    return shutil.disk_usage(p).free / 1e9
