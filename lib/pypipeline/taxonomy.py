"""
taxonomy.py — Thin wrapper around pytaxonkit for loading and querying lineage.

Adds pytaxonkit to sys.path automatically so this module works regardless of
where the pipeline folder lives relative to the Taxonkit folder.
"""

import os
import sys

# (bundled in archaeahq-update: pyassembly/pytaxonkit live next to this package under lib/)

from pytaxonkit.data import ensure_data_available, get_data_dir
from pytaxonkit.db import open_db, load_all
from pytaxonkit.traversal import TaxonomyTree
from pytaxonkit.formatter import CANONICAL_RANKS, RANK_ALIASES

# Fields exported in the taxonomy output and appended to the combined output.
# tax_id is the first column in the taxonomy-only output.
TAXONOMY_COLS = [
    "full_lineage",
    "lineage_taxids",
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


def load_tree(
    data_dir_override: str = None,
    rebuild: bool = False,
    verbose: bool = False,
) -> TaxonomyTree:
    """
    Ensure taxonomy data is present (downloading if needed), then load the
    entire NCBI taxonomy tree into memory as a TaxonomyTree.

    One-time cost: ~3 s and ~150 MB RAM.  All lineage queries thereafter are
    pure in-memory dict lookups.
    """
    data_dir = get_data_dir(data_dir_override)
    ensure_data_available(data_dir, force_download=rebuild, verbose=verbose)
    conn = open_db(data_dir)
    parent_of, rank_of, name_of, merged, deleted = load_all(conn)
    conn.close()
    return TaxonomyTree(parent_of, rank_of, name_of, merged, deleted)


def lineage_dict(tree: TaxonomyTree, tax_id) -> dict:
    """
    Return a flat dict with full lineage + canonical ranks for a taxID.

    Keys match TAXONOMY_COLS exactly.  Returns a dict of empty strings if
    tax_id is absent, invalid, deleted, or not found.
    """
    empty = {col: "" for col in TAXONOMY_COLS}

    if not tax_id:
        return empty
    try:
        tid = int(tax_id)
    except (ValueError, TypeError):
        return empty

    result = tree.get_lineage(tid)
    if not result.names:
        return empty

    full_lineage  = ";".join(result.names)
    lineage_taxids = ";".join(str(t) for t in result.taxids)

    # Build canonical rank → name mapping (handles RANK_ALIASES)
    rank_to_name: dict = {}
    for rank, name in zip(result.ranks, result.names):
        canonical = RANK_ALIASES.get(rank.lower(), rank.lower())
        if canonical in CANONICAL_RANKS and canonical not in rank_to_name:
            rank_to_name[canonical] = name

    rec = {"full_lineage": full_lineage, "lineage_taxids": lineage_taxids}
    for r in CANONICAL_RANKS:
        rec[r] = rank_to_name.get(r, "")

    return rec
