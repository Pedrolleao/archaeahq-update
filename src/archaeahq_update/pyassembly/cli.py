"""
cli.py — argparse definitions and main dispatch.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List, Set

from .db import open_db, get_cached, store_results
from .fetcher import fetch_metadata, BATCH_SIZE
from .formatter import HEADER, format_record

DEFAULT_DATA_DIR = "~/.assembly_meta"


def _read_accessions(source: str) -> List[str]:
    """Read one accession per line from a file path or stdin ('-')."""
    if source == "-":
        lines = sys.stdin.read().splitlines()
    else:
        try:
            with open(source) as fh:
                lines = fh.read().splitlines()
        except FileNotFoundError:
            sys.exit(f"ERROR: Input file not found: {source}")
    accessions = [ln.strip() for ln in lines if ln.strip()]
    if not accessions:
        sys.exit("ERROR: No accessions provided.")
    return accessions


def _load_written_accessions(output_path: str) -> Set[str]:
    """
    Read an existing (partial) output file and return the set of accessions
    already written to it.  Returns an empty set if the file doesn't exist,
    is empty, or has a different header (i.e. a different run).
    """
    written: Set[str] = set()
    try:
        with open(output_path) as fh:
            first = fh.readline().rstrip("\n")
            if first != HEADER:
                # File is from a different format / unrelated run — don't resume
                return set()
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                acc = line.split("\t")[0]
                if acc:
                    written.add(acc)
    except FileNotFoundError:
        pass
    return written


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="assembly_meta",
        description="Retrieve NCBI assembly metadata for GCA/GCF accessions.",
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # ── fetch subcommand ──────────────────────────────────────────────────────
    fetch = sub.add_parser(
        "fetch",
        help="Fetch metadata for a list of assembly accessions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Read GCA_*/GCF_* accessions (one per line) from FILE or stdin ('-'),\n"
            "look them up in the local SQLite cache, fetch any cache misses from\n"
            "the NCBI Datasets API, and write a TSV to stdout (or --output FILE).\n\n"
            "If --output FILE already exists from a previous interrupted run,\n"
            "the tool resumes automatically: it skips accessions already in the\n"
            "file and appends the remaining rows."
        ),
    )
    fetch.add_argument(
        "input",
        metavar="FILE",
        help="Input file with one accession per line, or '-' for stdin.",
    )
    fetch.add_argument(
        "-o", "--output",
        metavar="FILE",
        default="-",
        help="Output file (default: stdout). Supports automatic resume on re-run.",
    )
    fetch.add_argument(
        "--data-dir",
        metavar="DIR",
        default=DEFAULT_DATA_DIR,
        help=f"Cache directory (default: {DEFAULT_DATA_DIR}).",
    )
    fetch.add_argument(
        "--rebuild",
        action="store_true",
        help="Ignore cache and output file; re-fetch everything from the API.",
    )
    fetch.add_argument(
        "--api-key",
        metavar="KEY",
        default=None,
        help="NCBI API key (overrides $NCBI_API_KEY env var).",
    )
    fetch.add_argument(
        "--no-progress",
        action="store_true",
        help="Suppress tqdm progress bar.",
    )
    fetch.add_argument(
        "--verbose",
        action="store_true",
        help="Print status messages to stderr.",
    )

    return parser


def _cmd_fetch(args: argparse.Namespace) -> int:
    accessions = _read_accessions(args.input)

    api_key = args.api_key or os.environ.get("NCBI_API_KEY")
    show_progress = not args.no_progress
    verbose = args.verbose

    if verbose:
        print(f"  Read {len(accessions)} accessions.", file=sys.stderr)

    # ── resume detection ──────────────────────────────────────────────────────
    # Only possible when writing to a named file (not stdout).
    already_written: Set[str] = set()
    resuming = False

    if args.output != "-" and not args.rebuild:
        already_written = _load_written_accessions(args.output)
        if already_written:
            resuming = True
            if verbose:
                print(
                    f"  Resuming: {len(already_written)} accessions already in"
                    f" '{args.output}', skipping them.",
                    file=sys.stderr,
                )

    # ── open output ───────────────────────────────────────────────────────────
    if args.output == "-":
        out_fh = sys.stdout
        opened = False
    elif resuming:
        out_fh = open(args.output, "a")   # append to existing partial file
        opened = True
    else:
        out_fh = open(args.output, "w")   # fresh file
        opened = True

    conn = open_db(args.data_dir)

    try:
        # Write header only for fresh (non-resume) runs
        if not resuming:
            out_fh.write(HEADER + "\n")
            out_fh.flush()

        # ── chunk accessions; process + stream each batch ─────────────────────
        # Filter out already-written accessions so they're skipped in output,
        # but keep them in the full list for counting purposes.
        pending = [a for a in accessions if a not in already_written]

        if not pending:
            if verbose:
                print("  All accessions already in output. Nothing to do.", file=sys.stderr)
            return 0

        batches = [
            pending[i: i + BATCH_SIZE]
            for i in range(0, len(pending), BATCH_SIZE)
        ]

        total_cached = 0
        total_fetched = 0

        try:
            from tqdm import tqdm as _tqdm
            iterator = (
                _tqdm(
                    batches,
                    desc="Batches",
                    unit="batch",
                    file=sys.stderr,
                    initial=len(already_written) // BATCH_SIZE,   # show resume progress
                    total=len(accessions) // BATCH_SIZE + 1,
                )
                if show_progress
                else batches
            )
        except ImportError:
            iterator = batches

        for batch in iterator:
            # Cache lookup for this batch only (≤200 accessions)
            if args.rebuild:
                batch_cache: dict = {}
            else:
                batch_cache = get_cached(conn, batch)

            total_cached += len(batch_cache)

            missing = [a for a in batch if a not in batch_cache]

            if missing:
                new_records = fetch_metadata(
                    missing,
                    api_key=api_key,
                    show_progress=False,   # tqdm managed here
                    verbose=verbose,
                )
                store_results(conn, new_records)
                for rec in new_records:
                    batch_cache[rec["accession"]] = rec
                total_fetched += len(new_records)

            # Write this batch's rows immediately, preserving input order
            for acc in batch:
                rec = batch_cache.get(acc, {"accession": acc})
                out_fh.write(format_record(rec) + "\n")
            out_fh.flush()

        if verbose:
            skipped = len(already_written)
            print(
                f"  Done. Skipped (already written): {skipped}  |"
                f"  Cache hits: {total_cached}  |"
                f"  API fetched: {total_fetched}  |"
                f"  Total rows: {len(accessions)}",
                file=sys.stderr,
            )

    finally:
        conn.close()
        if opened:
            out_fh.close()

    return 0


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "fetch":
        return _cmd_fetch(args)

    parser.print_help()
    return 1
