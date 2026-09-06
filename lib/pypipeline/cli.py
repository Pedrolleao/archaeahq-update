"""
cli.py — argparse definitions and main dispatch for the unified pipeline.
"""

import argparse
import os
import sys
from pathlib import Path
from typing import List


DEFAULT_OUTPUT_DIR   = "."
DEFAULT_ASSEMBLY_OUT = "assembly_metadata.tsv"
DEFAULT_TAXONOMY_OUT = "taxonomy_lineage.tsv"
DEFAULT_COMBINED_OUT = "combined.tsv"

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_ASSEMBLY_DATA_DIR = os.path.normpath(os.path.join(_HERE, "..", "assembly_data"))
DEFAULT_TAXONOMY_DATA_DIR = os.path.normpath(os.path.join(_HERE, "..", "taxonkit_data"))


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


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="pipeline",
        description=(
            "Unified pipeline: given GCA_*/GCF_* assembly accessions, retrieve\n"
            "assembly metadata + BioSample environmental context + full NCBI\n"
            "taxonomic lineage. Writes three separate TSV files and one combined\n"
            "file containing all data together."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", metavar="COMMAND")

    # ── run subcommand ────────────────────────────────────────────────────────
    run_p = sub.add_parser(
        "run",
        help="Run the full pipeline on a list of assembly accessions.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "Read GCA_*/GCF_* accessions from FILE (or stdin with '-'), fetch\n"
            "assembly metadata and BioSample environmental attributes from NCBI,\n"
            "look up full taxonomic lineage from the local NCBI taxonomy tree,\n"
            "and write three output files:\n\n"
            "  assembly_metadata.tsv   26-column assembly + BioSample TSV\n"
            "  taxonomy_lineage.tsv    11-column per-taxID lineage TSV\n"
            "  combined.tsv            36-column merged TSV (all data)\n\n"
            "If any output file already exists from a previous interrupted run,\n"
            "the pipeline resumes automatically from where it stopped."
        ),
    )

    # Input
    run_p.add_argument(
        "input",
        metavar="FILE",
        help="Input file with one GCA_*/GCF_* accession per line, or '-' for stdin.",
    )

    # Output locations
    out_group = run_p.add_argument_group("output")
    out_group.add_argument(
        "--output-dir", "-d",
        metavar="DIR",
        default=DEFAULT_OUTPUT_DIR,
        help=(
            f"Directory for all output files (default: current directory).\n"
            f"Individual filenames can be overridden with the flags below."
        ),
    )
    out_group.add_argument(
        "--assembly-out",
        metavar="FILE",
        default=None,
        help=f"Assembly metadata output filename (default: {DEFAULT_ASSEMBLY_OUT}).",
    )
    out_group.add_argument(
        "--taxonomy-out",
        metavar="FILE",
        default=None,
        help=f"Taxonomy lineage output filename (default: {DEFAULT_TAXONOMY_OUT}).",
    )
    out_group.add_argument(
        "--combined-out",
        metavar="FILE",
        default=None,
        help=f"Combined output filename (default: {DEFAULT_COMBINED_OUT}).",
    )

    # Data directories
    data_group = run_p.add_argument_group("data")
    data_group.add_argument(
        "--data-dir",
        metavar="DIR",
        default=DEFAULT_ASSEMBLY_DATA_DIR,
        help="Assembly metadata cache directory (default: ./assembly_data/ inside the project).",
    )
    data_group.add_argument(
        "--taxonomy-dir",
        metavar="DIR",
        default=None,
        help="NCBI taxonomy data directory (default: ./taxonkit_data/ inside the project, or $TAXONKIT_DB if set).",
    )

    # Behaviour
    run_p.add_argument(
        "--rebuild",
        action="store_true",
        help="Ignore all caches and output files; re-fetch and rewrite everything.",
    )
    run_p.add_argument(
        "--api-key",
        metavar="KEY",
        default=None,
        help="NCBI API key (overrides $NCBI_API_KEY; doubles assembly fetch rate limit).",
    )
    run_p.add_argument(
        "--no-progress",
        action="store_true",
        help="Suppress tqdm progress bar.",
    )
    run_p.add_argument(
        "--verbose",
        action="store_true",
        help="Print detailed status messages to stderr.",
    )

    return parser


def _resolve_output_path(output_dir: str, filename_override, default_name: str) -> str:
    """Return absolute path for an output file, honouring --output-dir."""
    if filename_override:
        # If the override is already an absolute path or includes a directory,
        # use it as-is; otherwise place it inside output_dir.
        p = Path(filename_override)
        if p.is_absolute() or p.parent != Path("."):
            return str(p)
        return str(Path(output_dir) / p)
    return str(Path(output_dir) / default_name)


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return 1

    if args.command == "run":
        return _cmd_run(args)

    parser.print_help()
    return 1


def _cmd_run(args: argparse.Namespace) -> int:
    from pypipeline.runner import run

    accessions = _read_accessions(args.input)
    api_key    = args.api_key or os.environ.get("NCBI_API_KEY")
    taxonomy_data_dir = (
        args.taxonomy_dir
        or os.environ.get("TAXONKIT_DB")
        or DEFAULT_TAXONOMY_DATA_DIR
    )

    # Resolve output paths
    output_dir = str(Path(args.output_dir).expanduser())
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    assembly_out = _resolve_output_path(output_dir, args.assembly_out, DEFAULT_ASSEMBLY_OUT)
    taxonomy_out = _resolve_output_path(output_dir, args.taxonomy_out, DEFAULT_TAXONOMY_OUT)
    combined_out = _resolve_output_path(output_dir, args.combined_out, DEFAULT_COMBINED_OUT)

    if args.verbose:
        print(f"  Input accessions : {len(accessions)}", file=sys.stderr)
        print(f"  Assembly output  : {assembly_out}",    file=sys.stderr)
        print(f"  Taxonomy output  : {taxonomy_out}",    file=sys.stderr)
        print(f"  Combined output  : {combined_out}",    file=sys.stderr)

    return run(
        accessions       = accessions,
        assembly_data_dir= args.data_dir,
        taxonomy_data_dir= taxonomy_data_dir,
        api_key          = api_key,
        rebuild          = args.rebuild,
        verbose          = args.verbose,
        show_progress    = not args.no_progress,
        assembly_out_path= assembly_out,
        taxonomy_out_path= taxonomy_out,
        combined_out_path= combined_out,
    )
