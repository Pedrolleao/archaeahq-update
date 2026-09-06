"""Output format functions for lineage results."""

from typing import List

from ..pytaxonkit.traversal import LineageResult

# Canonical 8 ranks used by --format canonical
CANONICAL_RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]

# 2025 NCBI restructuring: map newer/alternate rank names → canonical
# "domain" replaces "superkingdom" for Bacteria/Archaea/Eukaryota
# "acellular root" used for Viruses
RANK_ALIASES: dict = {
    "domain": "superkingdom",
    "acellular root": "superkingdom",
}

# Headers for each format mode
LINEAGE_HEADER = "taxid\tlineage"
CANONICAL_HEADER = (
    "taxid\tsuperkingdom\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies"
)
TAXIDS_HEADER = "taxid\tlineage\ttaxids"


def _rank_to_canonical(rank: str) -> str:
    """Normalize a rank string to a canonical rank (or return as-is)."""
    r = rank.lower()
    return RANK_ALIASES.get(r, r)


def format_lineage(result: LineageResult, delimiter: str = ";") -> str:
    """
    Return a tab-delimited line:  input_taxid\\tname1;name2;...

    For deleted/not_found/empty results the lineage column is empty.
    """
    if not result.names:
        return f"{result.input_taxid}\t"
    lineage_str = delimiter.join(result.names)
    return f"{result.input_taxid}\t{lineage_str}"


def format_canonical(result: LineageResult, delimiter: str = ";") -> str:
    """
    Return a tab-delimited line with 9 columns:
      input_taxid  SK  K  P  C  O  F  G  S

    Missing ranks are represented as empty strings.
    Applies RANK_ALIASES so that "domain" maps to superkingdom column.
    """
    # Build rank→name mapping from result
    rank_to_name: dict = {}
    for rank, name in zip(result.ranks, result.names):
        canonical = _rank_to_canonical(rank)
        if canonical in CANONICAL_RANKS and canonical not in rank_to_name:
            rank_to_name[canonical] = name

    cols = [str(result.input_taxid)]
    for rank in CANONICAL_RANKS:
        cols.append(rank_to_name.get(rank, ""))

    return "\t".join(cols)


def format_taxids(result: LineageResult, delimiter: str = ";") -> str:
    """
    Return a tab-delimited line:
      input_taxid\\tname1;name2;...\\ttaxid1;taxid2;...

    For deleted/not_found/empty results the lineage columns are empty.
    """
    if not result.names:
        return f"{result.input_taxid}\t\t"
    name_str = delimiter.join(result.names)
    id_str = delimiter.join(str(t) for t in result.taxids)
    return f"{result.input_taxid}\t{name_str}\t{id_str}"


def format_result(
    result: LineageResult,
    mode: str = "lineage",
    delimiter: str = ";",
) -> str:
    """Dispatch to the appropriate formatter based on mode."""
    if mode == "lineage":
        return format_lineage(result, delimiter)
    elif mode == "canonical":
        return format_canonical(result, delimiter)
    elif mode == "taxids":
        return format_taxids(result, delimiter)
    else:
        raise ValueError(f"Unknown format mode: {mode!r}")
