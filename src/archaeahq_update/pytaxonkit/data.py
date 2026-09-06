"""Download, extract, and parse NCBI taxonomy dump files."""

import os
import sys
import tarfile
import urllib.request
from pathlib import Path
from typing import Generator, Optional, Tuple

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"
NEEDED_FILES = {"nodes.dmp", "names.dmp", "merged.dmp", "delnodes.dmp"}
DB_FILENAME = "taxonomy.db"


def get_data_dir(override: Optional[str] = None) -> Path:
    """Resolve the data directory from override, env var, or default."""
    if override:
        return Path(override).expanduser().resolve()
    env = os.environ.get("TAXONKIT_DB")
    if env:
        return Path(env).expanduser().resolve()
    return Path.home() / ".taxonkit"


def ensure_data_available(
    data_dir: Path,
    force_download: bool = False,
    verbose: bool = False,
) -> None:
    """
    Ensure taxonomy data is available in data_dir.

    Check order:
    1. taxonomy.db already exists (and not force_download) → done
    2. .dmp files already present → build index only
    3. Neither → download, extract, then build index
    """
    from ..pytaxonkit.db import db_is_valid, build_index, open_db

    data_dir.mkdir(parents=True, exist_ok=True)
    db_path = data_dir / DB_FILENAME

    if not force_download and db_path.exists():
        conn = open_db(data_dir)
        if db_is_valid(conn):
            conn.close()
            if verbose:
                print(f"[data] Using existing database: {db_path}", file=sys.stderr)
            return
        conn.close()
        if verbose:
            print("[data] Database invalid or incomplete, rebuilding.", file=sys.stderr)

    dmp_files_present = all((data_dir / f).exists() for f in NEEDED_FILES)
    if not dmp_files_present or force_download:
        download_taxdump(data_dir, verbose=verbose)
        extract_taxdump(data_dir, verbose=verbose)

    build_index(data_dir, verbose=verbose)


def download_taxdump(data_dir: Path, verbose: bool = False, max_retries: int = 10) -> None:
    """
    Download taxdump.tar.gz from NCBI FTP with a progress bar and resume support.

    Uses HTTP Range requests to resume partial downloads after connection drops.
    Retries up to max_retries times.
    """
    dest = data_dir / "taxdump.tar.gz"

    try:
        from tqdm import tqdm
        _has_tqdm = True
    except ImportError:
        _has_tqdm = False

    CHUNK = 1 << 17  # 128 KB per read

    for attempt in range(1, max_retries + 1):
        existing = dest.stat().st_size if dest.exists() else 0

        # Get total size via HEAD request
        try:
            req_head = urllib.request.Request(TAXDUMP_URL, method="HEAD")
            with urllib.request.urlopen(req_head, timeout=30) as resp:
                total_size = int(resp.headers.get("Content-Length", 0))
        except Exception:
            total_size = 0

        if existing > 0 and existing == total_size:
            if verbose:
                print(f"[data] Already complete: {dest}", file=sys.stderr)
            return

        if verbose or attempt > 1:
            print(
                f"[data] Download attempt {attempt}/{max_retries} "
                f"(resuming from {existing:,} / {total_size:,} bytes): {dest}",
                file=sys.stderr,
            )

        headers = {}
        if existing > 0:
            headers["Range"] = f"bytes={existing}-"

        req = urllib.request.Request(TAXDUMP_URL, headers=headers)

        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                mode = "ab" if existing > 0 else "wb"
                bar_initial = existing

                if _has_tqdm:
                    pbar = tqdm(
                        unit="B",
                        unit_scale=True,
                        unit_divisor=1024,
                        miniters=1,
                        desc="Downloading taxdump",
                        total=total_size if total_size > 0 else None,
                        initial=bar_initial,
                        file=sys.stderr,
                    )

                with open(dest, mode) as fh:
                    while True:
                        chunk = resp.read(CHUNK)
                        if not chunk:
                            break
                        fh.write(chunk)
                        if _has_tqdm:
                            pbar.update(len(chunk))
                        else:
                            downloaded = dest.stat().st_size
                            if total_size > 0:
                                pct = min(100, downloaded * 100 // total_size)
                                print(
                                    f"\r[data] Download: {pct}% ({downloaded:,}/{total_size:,} B)",
                                    end="",
                                    file=sys.stderr,
                                    flush=True,
                                )

                if _has_tqdm:
                    pbar.close()
                else:
                    print(file=sys.stderr)

        except Exception as exc:
            if _has_tqdm:
                try:
                    pbar.close()
                except Exception:
                    pass
            print(f"\n[data] Download error (attempt {attempt}): {exc}", file=sys.stderr)
            if attempt == max_retries:
                raise
            import time
            time.sleep(2 * attempt)
            continue

        # Verify completeness
        final_size = dest.stat().st_size
        if total_size > 0 and final_size < total_size:
            print(
                f"[data] Incomplete download ({final_size:,}/{total_size:,}), retrying...",
                file=sys.stderr,
            )
            if attempt == max_retries:
                raise RuntimeError(
                    f"Download incomplete after {max_retries} attempts: "
                    f"{final_size}/{total_size} bytes"
                )
            import time
            time.sleep(2 * attempt)
            continue

        if verbose:
            print(f"[data] Download complete: {dest}", file=sys.stderr)
        return

    raise RuntimeError(f"Failed to download {TAXDUMP_URL} after {max_retries} attempts")


def extract_taxdump(data_dir: Path, verbose: bool = False) -> None:
    """Extract only the 4 needed .dmp files from taxdump.tar.gz."""
    archive = data_dir / "taxdump.tar.gz"
    if verbose:
        print(f"[data] Extracting from {archive}", file=sys.stderr)

    with tarfile.open(archive, "r:gz") as tf:
        members = [m for m in tf.getmembers() if m.name in NEEDED_FILES]
        # Ensure we only extract into data_dir (security: strip any path components)
        for member in members:
            member.name = Path(member.name).name  # strip directories
            if verbose:
                print(f"[data]   Extracting {member.name}", file=sys.stderr)
        tf.extractall(path=data_dir, members=members)

    if verbose:
        print("[data] Extraction complete.", file=sys.stderr)


# ---------------------------------------------------------------------------
# .dmp parsers
# ---------------------------------------------------------------------------

def parse_nodes_dmp(path: Path) -> Generator[Tuple[int, int, str], None, None]:
    """
    Yield (taxid, parent_taxid, rank) from nodes.dmp.

    nodes.dmp field layout (tab-delimited, separator is '\\t|\\t'):
      0: tax_id
      1: parent tax_id
      2: rank
    """
    sep = "\t|\t"
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\t|\n").split(sep)
            if len(parts) < 3:
                continue
            try:
                taxid = int(parts[0])
                parent = int(parts[1])
            except ValueError:
                continue
            yield taxid, parent, parts[2]


def parse_names_dmp(path: Path) -> Generator[Tuple[int, str], None, None]:
    """
    Yield (taxid, scientific_name) from names.dmp.

    names.dmp field layout (tab-delimited, separator is '\\t|\\t'):
      0: tax_id
      1: name_txt
      2: unique_name
      3: name class
    Keep only rows where name class == 'scientific name'.
    """
    sep = "\t|\t"
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\t|\n").split(sep)
            if len(parts) < 4:
                continue
            if parts[3].strip() != "scientific name":
                continue
            try:
                taxid = int(parts[0])
            except ValueError:
                continue
            yield taxid, parts[1]


def parse_merged_dmp(path: Path) -> Generator[Tuple[int, int], None, None]:
    """
    Yield (old_taxid, new_taxid) from merged.dmp.

    merged.dmp field layout:
      0: old_tax_id
      1: new_tax_id
    """
    sep = "\t|\t"
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\t|\n").split(sep)
            if len(parts) < 2:
                continue
            try:
                old = int(parts[0])
                new = int(parts[1])
            except ValueError:
                continue
            yield old, new


def parse_delnodes_dmp(path: Path) -> Generator[int, None, None]:
    """
    Yield taxid from delnodes.dmp.

    delnodes.dmp field layout:
      0: tax_id
    """
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            parts = line.split("\t")
            if not parts:
                continue
            try:
                yield int(parts[0])
            except ValueError:
                continue
