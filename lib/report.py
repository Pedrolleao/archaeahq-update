"""
report.py — assemble the two output tables, the run summary, the Sankey input and the ledger.
"""

from __future__ import annotations

import csv
import time
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, List, Sequence

from common import DB_COLUMNS, FILTER_COLUMNS, write_tsv

FULL_COLUMNS = ["Recommendation", "Reason"] + DB_COLUMNS
FILTER_TABLE_COLUMNS = ["Assembly ID", "Archaeal_Kingdom", "Overall"] + FILTER_COLUMNS + ["Notes"]


def pf(value) -> str:
    if value is True:
        return "PASS"
    if value is False:
        return "FAIL"
    return "SKIPPED"


def write_full_table(path: Path, rows: Sequence[dict]) -> None:
    write_tsv(path, rows, FULL_COLUMNS)


def write_filter_table(path: Path, rows: Sequence[dict]) -> None:
    write_tsv(path, rows, FILTER_TABLE_COLUMNS)


def write_summary(path: Path, run_name: str, kingdom_order: Sequence[str], full_rows: Sequence[dict],
                  filter_rows: Sequence[dict], stage_counts: Dict[str, Dict[str, int]]) -> str:
    lines = [f"ArchaeaHQ update — run {run_name}", f"generated {time.strftime('%Y-%m-%d %H:%M')}", ""]
    header = f"{'Kingdom':<34}" + "".join(f"{s:>14}" for s in stage_counts)
    lines.append(header)
    for k in kingdom_order:
        lines.append(f"{k:<34}" + "".join(f"{stage_counts[s].get(k, 0):>14}" for s in stage_counts))
    lines.append(f"{'Total':<34}" + "".join(f"{sum(stage_counts[s].values()):>14}" for s in stage_counts))
    lines.append("")
    n_rep = sum(1 for r in full_rows if r["Recommendation"] == "Replace")
    if n_rep:
        lines.append(f"Replace (newer assembly version of a database genome): {n_rep}")
        for r in full_rows:
            if r["Recommendation"] == "Replace":
                lines.append(f"  {r['Assembly ID']:<18} {r['Reason']}")
        lines.append("")
    reasons = Counter(r["Reason"].split(";")[0].strip() for r in full_rows if r["Recommendation"] == "Do_not_add")
    lines.append("Reasons for 'Do_not_add' (first failed filter):")
    for reason, n in reasons.most_common():
        lines.append(f"  {n:>6}  {reason}")
    text = "\n".join(lines) + "\n"
    path.write_text(text, encoding="utf-8")
    return text


def write_sankey_input(path: Path, run_name: str, kingdom_order: Sequence[str],
                       stage_counts: Dict[str, Dict[str, int]]) -> None:
    """Input file for lib/sankey_generic.py (counts must be non-increasing per row)."""
    stages = list(stage_counts.keys())
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# TITLE: ArchaeaHQ update {run_name} – new NCBI assemblies\n")
        fh.write("# START_LABEL: New at NCBI\n")
        fh.write("# END_LABEL: Recommended to add\n")
        fh.write("group\t" + "\t".join(stages) + "\n")
        for k in kingdom_order:
            vals = [stage_counts[s].get(k, 0) for s in stages]
            # enforce monotonic non-increase (defensive; the stages are nested by construction)
            for i in range(1, len(vals)):
                vals[i] = min(vals[i], vals[i - 1])
            if vals[0] == 0:
                continue
            fh.write(k + "\t" + "\t".join(str(v) for v in vals) + "\n")


LEDGER_COLUMNS = ["Assembly ID", "Archaeal_Kingdom", "Evaluated_on", "Decision", "Reason", "Run"]


def append_ledger(path: Path, run_name: str, full_rows: Sequence[dict]) -> None:
    """Record every evaluated accession; rows of an earlier attempt of the same run are replaced."""
    existing: List[List[str]] = []
    if path.exists():
        with open(path, newline="", encoding="utf-8") as fh:
            rd = csv.reader(fh, delimiter="\t")
            next(rd, None)
            existing = [row for row in rd if row and (len(row) < 6 or row[5] != run_name)]
    today = time.strftime("%Y-%m-%d")
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(LEDGER_COLUMNS)
        w.writerows(existing)
        for r in full_rows:
            w.writerow([r["Assembly ID"], r["Archaeal_Kingdom"], today, r["Recommendation"], r["Reason"], run_name])


def stage_counts_from_rows(kingdom_order: Sequence[str], filter_rows: Sequence[dict]) -> Dict[str, Dict[str, int]]:
    """Nested per-kingdom counts along the filter chain (each stage counts genomes passing all previous ones)."""
    chain = [("New", None), ("Downloaded", "Downloaded"), ("Contam ≤10%", "Contamination_le10"),
             ("Compl ≥70%", "Completeness_ge70"), ("Not in DB", "Not_identical_to_DB"),
             ("New species", "Not_same_species_as_DB"), ("Batch rep.", "Representative_within_batch")]
    counts: Dict[str, Dict[str, int]] = OrderedDict()
    for label, col in chain:
        counts[label] = {k: 0 for k in kingdom_order}
    first = chain[0][0]
    for r in filter_rows:
        k = r["Archaeal_Kingdom"]
        if k not in counts[first]:
            for label, _ in chain:
                counts[label][k] = 0
        alive = True
        for label, col in chain:
            if col is not None and r.get(col) != "PASS":
                alive = False
            if alive:
                counts[label][k] += 1
    return counts
