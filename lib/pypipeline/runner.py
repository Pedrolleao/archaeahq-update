"""
runner.py — Core pipeline orchestration.

Data flow (per batch of 200 assembly accessions):
  1. Assembly metadata  → SQLite cache or NCBI Datasets API (batch of 200)
  2. Taxonomy lineage   → in-memory TaxonomyTree (instant, memoized)
  3. Write              → assembly_out, taxonomy_out, combined_out
  4. Flush all three output files

All three outputs are streamed batch-by-batch so the files are always in a
valid state — if the run crashes, every row that was written is complete.
"""

import os
import sys
from pathlib import Path
from typing import IO, Optional, Set

# (bundled in archaeahq-update: pyassembly/pytaxonkit live next to this package under lib/)

from pyassembly.db import open_db as open_assembly_db, get_cached, store_results
from pyassembly.fetcher import fetch_metadata, BATCH_SIZE

from pypipeline.taxonomy import load_tree, lineage_dict, TAXONOMY_COLS
from pypipeline.merger import (
    ASSEMBLY_HEADER,
    TAXONOMY_HEADER,
    COMBINED_HEADER,
    format_assembly_row,
    format_taxonomy_row,
    format_combined_row,
)


# ── Resume helpers ────────────────────────────────────────────────────────────

def _read_written_accessions(path: str, expected_header: str) -> Set[str]:
    """
    Read accessions already present in an existing output file.
    Returns empty set if the file doesn't exist or has a different header.
    """
    written: Set[str] = set()
    try:
        with open(path) as fh:
            if fh.readline().rstrip("\n") != expected_header:
                return set()
            for line in fh:
                acc = line.split("\t")[0].strip()
                if acc:
                    written.add(acc)
    except FileNotFoundError:
        pass
    return written


def _read_written_taxids(path: str) -> Set[str]:
    """
    Read tax_ids already present in the taxonomy output file.
    Returns empty set if the file doesn't exist or has a wrong header.
    """
    written: Set[str] = set()
    try:
        with open(path) as fh:
            if fh.readline().rstrip("\n") != TAXONOMY_HEADER:
                return set()
            for line in fh:
                tid = line.split("\t")[0].strip()
                if tid:
                    written.add(tid)
    except FileNotFoundError:
        pass
    return written


# ── Main runner ───────────────────────────────────────────────────────────────

def run(
    accessions: list,
    assembly_data_dir: str,
    taxonomy_data_dir: Optional[str],
    api_key: Optional[str],
    rebuild: bool,
    verbose: bool,
    show_progress: bool,
    assembly_out_path: str,
    taxonomy_out_path: str,
    combined_out_path: str,
) -> int:
    """
    Run the full pipeline: assembly metadata → taxonomy lineage → 3 outputs.

    Returns 0 on success.
    """

    # ── Resume detection ──────────────────────────────────────────────────────
    # Combined output is the canonical resume reference.
    already_written: Set[str] = set()
    written_taxids:  Set[str] = set()
    resuming = False

    if not rebuild:
        already_written = _read_written_accessions(combined_out_path, COMBINED_HEADER)
        written_taxids  = _read_written_taxids(taxonomy_out_path)
        if already_written:
            resuming = True
            if verbose:
                print(
                    f"  Resuming: {len(already_written)} accessions already in "
                    f"'{combined_out_path}', skipping them.",
                    file=sys.stderr,
                )

    # ── Open output files ─────────────────────────────────────────────────────
    mode = "a" if resuming else "w"
    asm_fh  = open(assembly_out_path,  mode)
    tax_fh  = open(taxonomy_out_path,  mode)
    comb_fh = open(combined_out_path,  mode)

    if not resuming:
        asm_fh.write(ASSEMBLY_HEADER  + "\n")
        tax_fh.write(TAXONOMY_HEADER  + "\n")
        comb_fh.write(COMBINED_HEADER + "\n")
        asm_fh.flush(); tax_fh.flush(); comb_fh.flush()

    # ── Load taxonomy tree (once) ─────────────────────────────────────────────
    if verbose:
        print("  Loading NCBI taxonomy tree into memory...", file=sys.stderr)
    tree = load_tree(
        data_dir_override=taxonomy_data_dir,
        rebuild=rebuild,
        verbose=verbose,
    )
    if verbose:
        print("  Taxonomy tree loaded.", file=sys.stderr)

    # ── Open assembly cache ───────────────────────────────────────────────────
    conn = open_assembly_db(assembly_data_dir)

    # ── Determine which accessions still need processing ──────────────────────
    pending = [a for a in accessions if a not in already_written]

    if not pending:
        if verbose:
            print("  All accessions already processed. Nothing to do.", file=sys.stderr)
        conn.close()
        asm_fh.close(); tax_fh.close(); comb_fh.close()
        return 0

    batches = [pending[i: i + BATCH_SIZE] for i in range(0, len(pending), BATCH_SIZE)]

    try:
        from tqdm import tqdm as _tqdm
        iterator = (
            _tqdm(
                batches,
                desc="Batches",
                unit="batch",
                file=sys.stderr,
                initial=len(already_written) // BATCH_SIZE,
                total=(len(accessions) + BATCH_SIZE - 1) // BATCH_SIZE,
            )
            if show_progress
            else batches
        )
    except ImportError:
        iterator = batches

    total_asm_cached = 0
    total_asm_fetched = 0
    total_tax_written = 0

    try:
        for batch in iterator:

            # ── Assembly metadata: cache → API ────────────────────────────────
            if rebuild:
                batch_cache: dict = {}
            else:
                batch_cache = get_cached(conn, batch)

            total_asm_cached += len(batch_cache)
            missing = [a for a in batch if a not in batch_cache]

            if missing:
                new_records = fetch_metadata(
                    missing,
                    api_key=api_key,
                    show_progress=False,
                    verbose=verbose,
                )
                store_results(conn, new_records)
                for rec in new_records:
                    batch_cache[rec["accession"]] = rec
                total_asm_fetched += len(new_records)

            # ── Write rows for this batch ─────────────────────────────────────
            for acc in batch:
                asm_rec  = batch_cache.get(acc, {"accession": acc})
                tax_id   = asm_rec.get("tax_id", "")
                tax_rec  = lineage_dict(tree, tax_id)

                # Assembly output (one row per accession)
                asm_fh.write(format_assembly_row(asm_rec) + "\n")

                # Taxonomy output (one row per unique taxID)
                tid_str = str(tax_id) if tax_id else ""
                if tid_str and tid_str not in written_taxids:
                    tax_fh.write(format_taxonomy_row(tax_id, tax_rec) + "\n")
                    written_taxids.add(tid_str)
                    total_tax_written += 1

                # Combined output (one row per accession, all 36 columns)
                comb_fh.write(format_combined_row(asm_rec, tax_rec) + "\n")

            # Flush all three files after each batch (crash-safe)
            asm_fh.flush()
            tax_fh.flush()
            comb_fh.flush()

    finally:
        conn.close()
        asm_fh.close()
        tax_fh.close()
        comb_fh.close()

    if verbose:
        skipped = len(already_written)
        total   = len(accessions)
        print(
            f"  Done.\n"
            f"    Accessions processed : {len(pending)}\n"
            f"    Accessions skipped   : {skipped}\n"
            f"    Assembly cache hits  : {total_asm_cached}\n"
            f"    Assembly API fetched : {total_asm_fetched}\n"
            f"    Unique taxIDs written: {total_tax_written}\n"
            f"    Total rows (combined): {total}",
            file=sys.stderr,
        )

    return 0
