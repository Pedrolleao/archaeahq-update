"""argparse CLI definition and command dispatch."""

import argparse
import io
import sys
from pathlib import Path
from typing import List, Optional


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="taxonkit",
        description="Python NCBI Taxonomy Lineage Tool",
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # ------------------------------------------------------------------
    # lineage sub-command
    # ------------------------------------------------------------------
    p_lin = sub.add_parser(
        "lineage",
        help="Query taxonomic lineage of given TaxIDs",
        description=(
            "Read TaxIDs from FILE (or stdin if FILE is '-'), "
            "and write their full lineage to stdout."
        ),
    )
    p_lin.add_argument(
        "file",
        nargs="?",
        default="-",
        metavar="FILE",
        help="Input file with one TaxID per line (default: stdin)",
    )
    p_lin.add_argument(
        "--format",
        choices=["lineage", "canonical", "taxids"],
        default="lineage",
        help=(
            "Output format: "
            "'lineage' = full name string; "
            "'canonical' = 8-rank table; "
            "'taxids' = name string + taxid chain "
            "(default: lineage)"
        ),
    )
    p_lin.add_argument(
        "-d",
        "--delimiter",
        default=";",
        metavar="SEP",
        help="Delimiter for lineage strings (default: ';')",
    )
    p_lin.add_argument(
        "-o",
        "--output",
        default=None,
        metavar="FILE",
        help="Output file (default: stdout)",
    )
    p_lin.add_argument(
        "--data-dir",
        default=None,
        metavar="DIR",
        help="Data directory (default: ~/.taxonkit/ or $TAXONKIT_DB)",
    )
    p_lin.add_argument(
        "--rebuild",
        action="store_true",
        help="Force rebuild of the SQLite index",
    )
    p_lin.add_argument(
        "--no-progress",
        action="store_true",
        help="Suppress tqdm progress bar",
    )
    p_lin.add_argument(
        "--verbose",
        action="store_true",
        help="Print verbose status messages to stderr",
    )

    return parser


def run_lineage(args: argparse.Namespace) -> int:
    """Main handler for the 'lineage' sub-command."""
    from ..pytaxonkit.data import get_data_dir, ensure_data_available
    from ..pytaxonkit.db import open_db, load_all
    from ..pytaxonkit.traversal import TaxonomyTree
    from ..pytaxonkit.formatter import format_result, LINEAGE_HEADER, CANONICAL_HEADER, TAXIDS_HEADER

    # 1. Resolve data dir & ensure DB exists
    data_dir = get_data_dir(args.data_dir)
    ensure_data_available(
        data_dir,
        force_download=args.rebuild,
        verbose=args.verbose,
    )

    # 2. Load tree into memory
    if args.verbose:
        print("[cli] Loading taxonomy into memory...", file=sys.stderr)
    conn = open_db(data_dir)
    parent_of, rank_of, name_of, merged, deleted = load_all(conn)
    conn.close()
    if args.verbose:
        print(
            f"[cli] Loaded {len(parent_of):,} nodes, "
            f"{len(name_of):,} names, "
            f"{len(merged):,} merges, "
            f"{len(deleted):,} deleted.",
            file=sys.stderr,
        )

    tree = TaxonomyTree(parent_of, rank_of, name_of, merged, deleted)

    # 3. Read TaxIDs
    taxids: List[int] = []
    input_source = args.file if args.file else "-"

    try:
        if input_source == "-":
            lines = sys.stdin.readlines()
        else:
            with open(input_source, encoding="utf-8") as fh:
                lines = fh.readlines()
    except OSError as exc:
        print(f"[ERROR] Cannot read input: {exc}", file=sys.stderr)
        return 1

    for line in lines:
        line = line.strip()
        if not line:
            continue
        try:
            taxids.append(int(line))
        except ValueError:
            # Skip header lines or non-integer tokens
            continue

    if args.verbose:
        print(f"[cli] {len(taxids):,} TaxIDs to process.", file=sys.stderr)

    # 4. Query lineages
    results = tree.get_lineages_batch(taxids, show_progress=not args.no_progress)

    # 5. Write output (buffered)
    FLUSH_EVERY = 1000
    out: io.IOBase
    if args.output:
        out = open(args.output, "w", encoding="utf-8", buffering=1 << 20)
    else:
        out = sys.stdout

    try:
        headers = {
            "lineage": LINEAGE_HEADER,
            "canonical": CANONICAL_HEADER,
            "taxids": TAXIDS_HEADER,
        }
        out.write(headers[args.format] + "\n")

        buf: List[str] = []
        for i, result in enumerate(results, 1):
            buf.append(format_result(result, mode=args.format, delimiter=args.delimiter))
            if i % FLUSH_EVERY == 0:
                out.write("\n".join(buf) + "\n")
                buf.clear()
        if buf:
            out.write("\n".join(buf) + "\n")
    finally:
        if args.output and out is not sys.stdout:
            out.close()

    return 0


def main(argv: Optional[List[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "lineage":
        return run_lineage(args)

    parser.print_help()
    return 0
