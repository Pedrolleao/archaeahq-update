"""
merger.py — Headers and row-format functions for the three pipeline outputs.

Imports pyassembly formatters directly so the pipeline never diverges from
what assembly_meta.py produces on its own.
"""

import os
import sys

# (bundled in archaeahq-update: pyassembly/pytaxonkit live next to this package under lib/)

from pyassembly.formatter import HEADER as _ASSEMBLY_HEADER, _FIELDS as _ASSEMBLY_FIELDS
from pypipeline.taxonomy import TAXONOMY_COLS

# ── Headers ───────────────────────────────────────────────────────────────────

ASSEMBLY_HEADER = _ASSEMBLY_HEADER   # 26-column assembly + biosample TSV header

TAXONOMY_HEADER = "\t".join(["tax_id"] + TAXONOMY_COLS)  # 11 columns

# Combined = all 26 assembly cols + 10 taxonomy cols (tax_id already present)
COMBINED_HEADER = ASSEMBLY_HEADER + "\t" + "\t".join(TAXONOMY_COLS)  # 36 columns


# ── Row formatters ─────────────────────────────────────────────────────────────

def format_assembly_row(rec: dict) -> str:
    """26-column assembly + biosample TSV row."""
    values = []
    for key, _ in _ASSEMBLY_FIELDS:
        val = rec.get(key, "")
        values.append("" if val is None else str(val))
    return "\t".join(values)


def format_taxonomy_row(tax_id, tax_rec: dict) -> str:
    """11-column taxonomy TSV row (tax_id + 10 lineage/rank columns)."""
    cols = [str(tax_id)]
    for field in TAXONOMY_COLS:
        val = tax_rec.get(field, "")
        cols.append("" if val is None else str(val))
    return "\t".join(cols)


def format_combined_row(asm_rec: dict, tax_rec: dict) -> str:
    """36-column combined TSV row (26 assembly + 10 taxonomy)."""
    asm_part = format_assembly_row(asm_rec)
    tax_part = "\t".join(
        "" if tax_rec.get(f) is None else str(tax_rec.get(f, ""))
        for f in TAXONOMY_COLS
    )
    return asm_part + "\t" + tax_part
