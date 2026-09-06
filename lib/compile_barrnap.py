#!/usr/bin/env python3
"""
compile_barrnap.py
------------------
Compile barrnap GFF outputs from a folder into:
  1. A single table with all RNA annotations + MAG identifier column.
  2. A per-MAG statistics table of RNA type counts.

Optional (--trna flag):
  3. A detailed tRNA-only table with Info column split into individual fields.
  4. A per-MAG amino acid coverage table (% of tRNAs per amino acid + total count).

Optional (--rrna flag):
  5. A detailed rRNA-only table with Info column split into individual fields.
  6. A per-MAG rRNA type coverage table (% of rRNAs per type + total count).

Usage:
    python compile_barrnap.py <input_folder> [--output-prefix PREFIX] [--trna] [--rrna]

Arguments:
    input_folder      Folder containing barrnap .gff files.
    --output-prefix   Prefix for output files (default: barrnap_compiled).
    --trna            Also generate tRNA-specific output files.
    --rrna            Also generate rRNA-specific output files.

Outputs (always):
    <prefix>_table.tsv           - Full annotation table.
    <prefix>_stats.tsv           - Per-MAG RNA type count statistics.

Outputs (with --trna):
    <prefix>_tRNA_table.tsv      - tRNA-only table with split Info columns.
    <prefix>_tRNA_aa_pct.tsv     - Per-MAG amino acid % coverage + total tRNA count.

Outputs (with --rrna):
    <prefix>_rRNA_table.tsv      - rRNA-only table with split Info columns.
    <prefix>_rRNA_types_pct.tsv  - Per-MAG rRNA type % coverage + total rRNA count.
"""

import os
import re
import sys
import argparse
import glob
from collections import defaultdict


# ─── Reference lists ──────────────────────────────────────────────────────────

# The 20 canonical amino acids in alphabetical order by 3-letter code.
# Non-standard tRNA types are appended after so no information is lost.
AA_20 = [
    "Ala", "Arg", "Asn", "Asp", "Cys",
    "Gln", "Glu", "Gly", "His", "Ile",
    "Leu", "Lys", "Met", "Phe", "Pro",
    "Ser", "Thr", "Trp", "Tyr", "Val",
]
AA_EXTRA = ["Pyl", "SeC", "Stop", "???"]

# Known rRNA subunit types in size order (prokaryotic + eukaryotic).
# Any type found in the data that is not listed here will be appended
# alphabetically at the end, so no information is ever lost.
RRNA_KNOWN = [
    "5S_rRNA",
    "5.8S_rRNA",
    "16S_rRNA",
    "18S_rRNA",
    "23S_rRNA",
    "28S_rRNA",
]


# ─── Filename parsing ──────────────────────────────────────────────────────────

def extract_mag_name(filepath):
    """Extract MAG identifier from the barrnap GFF filename."""
    basename = os.path.basename(filepath)
    if ".fna-barrnap" in basename:
        mag = basename.split(".fna-barrnap")[0]
    else:
        mag = os.path.splitext(basename)[0]
    return mag


# ─── GFF parsing ──────────────────────────────────────────────────────────────

def parse_gff_file(filepath):
    """
    Parse a barrnap GFF file.
    Returns a list of dicts with keys:
        mag, scaffold, rna_type, start, stop, strand, info
    Skips comment lines and bare 'gene' stub rows (no coordinates).
    """
    mag = extract_mag_name(filepath)
    records = []

    with open(filepath, "r") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#") or not line.strip():
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            scaffold = parts[0]
            rna_type = parts[2]
            start    = parts[3].strip()
            stop     = parts[4].strip()
            strand   = parts[6]
            info     = parts[8]

            # Drop aragorn companion 'gene' stubs (no coordinates)
            if rna_type == "gene" and (start == "" or stop == ""):
                continue

            records.append({
                "mag":      mag,
                "scaffold": scaffold,
                "rna_type": rna_type,
                "start":    start,
                "stop":     stop,
                "strand":   strand,
                "info":     info,
            })

    return records


def compile_folder(folder):
    """Read all .gff files in folder and return all records."""
    pattern = os.path.join(folder, "*.gff")
    gff_files = sorted(glob.glob(pattern))

    if not gff_files:
        print(f"ERROR: No .gff files found in '{folder}'", file=sys.stderr)
        sys.exit(1)

    print(f"Found {len(gff_files)} GFF files in '{folder}'", file=sys.stderr)

    all_records = []
    for fp in gff_files:
        records = parse_gff_file(fp)
        all_records.extend(records)
        print(f"  {os.path.basename(fp)}: {len(records)} records", file=sys.stderr)

    return all_records


# ─── Info field parsers ────────────────────────────────────────────────────────

def _split_info(info):
    """Split a GFF9 Info string into a key→value dict."""
    fields = {}
    for part in info.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            fields[k.strip()] = v.strip()
    return fields


def parse_trna_info(info):
    """
    Parse the GFF Info field of a tRNA record.

    Example: Name=tRNA-Ala;product=transfer RNA (cgc)

    Returns dict: tRNA_name, amino_acid, anticodon, product.
    """
    fields = _split_info(info)

    name = fields.get("Name", "")
    m = re.match(r"tRNA-(.+)$", name)
    amino_acid = m.group(1) if m else name

    product = fields.get("product", "")
    m2 = re.search(r"\(([^)]+)\)", product)
    anticodon = m2.group(1) if m2 else ""

    return {
        "tRNA_name":  name,
        "amino_acid": amino_acid,
        "anticodon":  anticodon,
        "product":    product,
    }


def parse_rrna_info(info):
    """
    Parse the GFF Info field of an rRNA record.

    Example: Name=16S_rRNA;Alias=SSU_rRNA_archaea;Dbxref=Rfam:RF01959;product=16S ribosomal RNA

    Returns dict: rRNA_name, alias, rfam_id, product.
    """
    fields = _split_info(info)

    name    = fields.get("Name", "")
    alias   = fields.get("Alias", "")
    product = fields.get("product", "")

    # Extract Rfam accession from Dbxref (e.g. "Rfam:RF01959")
    dbxref  = fields.get("Dbxref", "")
    m = re.search(r"Rfam:(RF\d+)", dbxref)
    rfam_id = m.group(1) if m else ""

    return {
        "rRNA_name": name,
        "alias":     alias,
        "rfam_id":   rfam_id,
        "product":   product,
    }


# ─── Writers: core tables ─────────────────────────────────────────────────────

def write_table(records, out_path):
    """Write the compiled annotation table to a TSV file."""
    header = ["MAG", "Scaffold", "RNA_type", "Start", "Stop", "Strand", "Info"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in records:
            row = [
                r["mag"], r["scaffold"], r["rna_type"],
                r["start"], r["stop"], r["strand"], r["info"],
            ]
            fh.write("\t".join(row) + "\n")
    print(f"Annotation table written to:    {out_path}", file=sys.stderr)


def write_stats(records, out_path):
    """
    Write per-MAG RNA type count statistics.
    Columns: MAG, one per unique RNA type (sorted), Total.
    """
    rna_types_all = set()
    mag_counts = defaultdict(lambda: defaultdict(int))

    for r in records:
        mag_counts[r["mag"]][r["rna_type"]] += 1
        rna_types_all.add(r["rna_type"])

    rna_types_sorted = sorted(rna_types_all)
    mags_sorted = sorted(mag_counts.keys())

    header = ["MAG"] + rna_types_sorted + ["Total"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for mag in mags_sorted:
            counts = mag_counts[mag]
            row_counts = [counts.get(rt, 0) for rt in rna_types_sorted]
            total = sum(row_counts)
            row = [mag] + [str(c) for c in row_counts] + [str(total)]
            fh.write("\t".join(row) + "\n")

    print(f"Statistics table written to:    {out_path}", file=sys.stderr)


# ─── Writers: tRNA-specific tables ────────────────────────────────────────────

def write_trna_table(records, out_path):
    """
    Write a tRNA-only table with the Info column expanded into individual fields.

    Columns: MAG, Scaffold, Start, Stop, Strand, tRNA_name, Amino_acid, Anticodon, Product
    """
    trna_records = [r for r in records if r["rna_type"] == "tRNA"]

    header = [
        "MAG", "Scaffold", "Start", "Stop", "Strand",
        "tRNA_name", "Amino_acid", "Anticodon", "Product",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in trna_records:
            parsed = parse_trna_info(r["info"])
            row = [
                r["mag"], r["scaffold"], r["start"], r["stop"], r["strand"],
                parsed["tRNA_name"], parsed["amino_acid"],
                parsed["anticodon"], parsed["product"],
            ]
            fh.write("\t".join(row) + "\n")

    print(f"tRNA detail table written to:   {out_path}", file=sys.stderr)
    return trna_records


def write_trna_aa_pct(trna_records, out_path):
    """
    Write a per-MAG amino acid coverage table.

    Columns: MAG, <AA_20 + AA_EXTRA>, Total_tRNA
    Values: % of tRNAs in that MAG encoding each amino acid (2 d.p.).
    Last column: total tRNA count for that MAG.
    """
    mag_aa_counts = defaultdict(lambda: defaultdict(int))
    mag_total = defaultdict(int)

    for r in trna_records:
        parsed = parse_trna_info(r["info"])
        mag_aa_counts[r["mag"]][parsed["amino_acid"]] += 1
        mag_total[r["mag"]] += 1

    mags_sorted = sorted(mag_total.keys())
    all_aa_cols = AA_20 + AA_EXTRA
    header = ["MAG"] + all_aa_cols + ["Total_tRNA"]

    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for mag in mags_sorted:
            total = mag_total[mag]
            counts = mag_aa_counts[mag]
            pct_cells = [
                f"{counts.get(aa, 0) / total * 100:.2f}" if total else "0.00"
                for aa in all_aa_cols
            ]
            fh.write("\t".join([mag] + pct_cells + [str(total)]) + "\n")

    print(f"tRNA AA pct table written to:   {out_path}", file=sys.stderr)


# ─── Writers: rRNA-specific tables ────────────────────────────────────────────

def write_rrna_table(records, out_path):
    """
    Write an rRNA-only table with the Info column expanded into individual fields.

    Columns: MAG, Scaffold, Start, Stop, Length, Strand, rRNA_name, Alias, Rfam_ID, Product
    """
    rrna_records = [r for r in records if r["rna_type"] == "rRNA"]

    header = [
        "MAG", "Scaffold", "Start", "Stop", "Length", "Strand",
        "rRNA_name", "Alias", "Rfam_ID", "Product",
    ]
    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rrna_records:
            parsed = parse_rrna_info(r["info"])
            try:
                length = str(int(r["stop"]) - int(r["start"]) + 1)
            except ValueError:
                length = ""
            row = [
                r["mag"], r["scaffold"], r["start"], r["stop"], length, r["strand"],
                parsed["rRNA_name"], parsed["alias"],
                parsed["rfam_id"],   parsed["product"],
            ]
            fh.write("\t".join(row) + "\n")

    print(f"rRNA detail table written to:   {out_path}", file=sys.stderr)
    return rrna_records


def write_rrna_types_pct(rrna_records, out_path):
    """
    Write a per-MAG rRNA type coverage table.

    Columns: MAG, <RRNA_KNOWN + any extra types found in data>, Total_rRNA
    Values: % of rRNAs in that MAG belonging to each type (2 d.p.).
    Last column: total rRNA count for that MAG.

    All types in RRNA_KNOWN are always present as columns even if absent
    from the dataset. Any additional type found in the data is appended
    alphabetically at the end.
    """
    mag_type_counts = defaultdict(lambda: defaultdict(int))
    mag_total = defaultdict(int)
    observed_types = set()

    for r in rrna_records:
        parsed = parse_rrna_info(r["info"])
        rna_name = parsed["rRNA_name"]
        mag_type_counts[r["mag"]][rna_name] += 1
        mag_total[r["mag"]] += 1
        observed_types.add(rna_name)

    # Build column list: known types first (always present), then any extras
    extra_types = sorted(observed_types - set(RRNA_KNOWN))
    all_type_cols = RRNA_KNOWN + extra_types

    mags_sorted = sorted(mag_total.keys())
    header = ["MAG"] + all_type_cols + ["Total_rRNA"]

    with open(out_path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for mag in mags_sorted:
            total = mag_total[mag]
            counts = mag_type_counts[mag]
            pct_cells = [
                f"{counts.get(rt, 0) / total * 100:.2f}" if total else "0.00"
                for rt in all_type_cols
            ]
            fh.write("\t".join([mag] + pct_cells + [str(total)]) + "\n")

    print(f"rRNA types pct table written to: {out_path}", file=sys.stderr)


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Compile barrnap GFF outputs into summary tables.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "input_folder",
        help="Folder containing barrnap .gff output files.",
    )
    parser.add_argument(
        "--output-prefix",
        default="barrnap_compiled",
        help="Prefix for output files (default: barrnap_compiled).",
    )
    parser.add_argument(
        "--trna",
        action="store_true",
        help="Generate detailed tRNA tables (tRNA_table and tRNA AA%% coverage).",
    )
    parser.add_argument(
        "--rrna",
        action="store_true",
        help="Generate detailed rRNA tables (rRNA_table and rRNA type%% coverage).",
    )
    args = parser.parse_args()

    folder = args.input_folder.rstrip("/")
    if not os.path.isdir(folder):
        print(f"ERROR: '{folder}' is not a valid directory.", file=sys.stderr)
        sys.exit(1)

    records = compile_folder(folder)
    print(f"\nTotal records parsed: {len(records)}", file=sys.stderr)

    prefix = args.output_prefix
    write_table(records, f"{prefix}_table.tsv")
    write_stats(records, f"{prefix}_stats.tsv")

    if args.trna:
        trna_records = write_trna_table(records, f"{prefix}_tRNA_table.tsv")
        write_trna_aa_pct(trna_records, f"{prefix}_tRNA_aa_pct.tsv")

    if args.rrna:
        rrna_records = write_rrna_table(records, f"{prefix}_rRNA_table.tsv")
        write_rrna_types_pct(rrna_records, f"{prefix}_rRNA_types_pct.tsv")

    # ── Global summary to stdout ──
    rna_type_totals = defaultdict(int)
    for r in records:
        rna_type_totals[r["rna_type"]] += 1

    print("\n--- Global RNA type summary ---")
    for rt in sorted(rna_type_totals):
        print(f"  {rt}: {rna_type_totals[rt]}")
    print(f"  TOTAL: {len(records)}")

    if args.trna:
        trna_total = rna_type_totals.get("tRNA", 0)
        print("\n--- tRNA amino acid breakdown (all MAGs combined) ---")
        aa_totals = defaultdict(int)
        for r in records:
            if r["rna_type"] == "tRNA":
                aa_totals[parse_trna_info(r["info"])["amino_acid"]] += 1
        for aa in sorted(aa_totals):
            pct = aa_totals[aa] / trna_total * 100 if trna_total else 0
            print(f"  {aa}: {aa_totals[aa]} ({pct:.1f}%)")

    if args.rrna:
        rrna_total = rna_type_totals.get("rRNA", 0)
        print("\n--- rRNA type breakdown (all MAGs combined) ---")
        rrna_totals = defaultdict(int)
        for r in records:
            if r["rna_type"] == "rRNA":
                rrna_totals[parse_rrna_info(r["info"])["rRNA_name"]] += 1
        for rt in sorted(rrna_totals):
            pct = rrna_totals[rt] / rrna_total * 100 if rrna_total else 0
            print(f"  {rt}: {rrna_totals[rt]} ({pct:.1f}%)")


if __name__ == "__main__":
    main()
