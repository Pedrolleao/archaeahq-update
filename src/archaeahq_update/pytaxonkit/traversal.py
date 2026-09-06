"""In-memory taxonomy tree and lineage walker."""

import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Set


@dataclass
class LineageResult:
    input_taxid: int
    resolved_taxid: int
    names: List[str]
    ranks: List[str]
    taxids: List[int]
    status: str  # "ok" | "deleted" | "merged" | "not_found"


class TaxonomyTree:
    """
    Holds the entire NCBI taxonomy tree in memory as Python dicts.

    Load once at startup (~2–3 s, ~150 MB RAM); all subsequent lineage
    queries are pure dict lookups with result memoization.
    """

    def __init__(
        self,
        parent_of: Dict[int, int],
        rank_of: Dict[int, str],
        name_of: Dict[int, str],
        merged: Dict[int, int],
        deleted: Set[int],
    ) -> None:
        self._parent_of = parent_of
        self._rank_of = rank_of
        self._name_of = name_of
        self._merged = merged
        self._deleted = deleted
        self._cache: Dict[int, LineageResult] = {}

    def get_lineage(self, taxid: int) -> LineageResult:
        """
        Return the complete taxonomic lineage for taxid (top-down order).

        Handles deleted, merged, and not-found taxIDs gracefully,
        emitting warnings to stderr.
        """
        # Fast path: cache hit (common for clustered inputs)
        if taxid in self._cache:
            cached = self._cache[taxid]
            # If the original query used a different input_taxid (e.g. merged),
            # wrap with the correct input_taxid.
            if cached.input_taxid != taxid:
                return LineageResult(
                    input_taxid=taxid,
                    resolved_taxid=cached.resolved_taxid,
                    names=cached.names,
                    ranks=cached.ranks,
                    taxids=cached.taxids,
                    status=cached.status,
                )
            return cached

        # --- deleted? ---
        if taxid in self._deleted:
            print(
                f"[WARNING] taxid {taxid} was deleted", file=sys.stderr
            )
            result = LineageResult(
                input_taxid=taxid,
                resolved_taxid=0,
                names=[],
                ranks=[],
                taxids=[],
                status="deleted",
            )
            self._cache[taxid] = result
            return result

        # --- merged? ---
        resolved = taxid
        if taxid in self._merged:
            new_id = self._merged[taxid]
            print(
                f"[WARNING] taxid {taxid} was merged into {new_id}",
                file=sys.stderr,
            )
            resolved = new_id
            # Check if the merged target is in cache already
            if resolved in self._cache:
                cached = self._cache[resolved]
                result = LineageResult(
                    input_taxid=taxid,
                    resolved_taxid=cached.resolved_taxid,
                    names=cached.names,
                    ranks=cached.ranks,
                    taxids=cached.taxids,
                    status="merged",
                )
                self._cache[taxid] = result
                return result

        # --- walk the tree bottom-up ---
        names_rev: List[str] = []
        ranks_rev: List[str] = []
        taxids_rev: List[int] = []

        current = resolved
        visited: Set[int] = set()  # cycle guard (should never happen, but safety first)

        while True:
            if current in visited:
                # Cycle detected — shouldn't happen with valid NCBI data
                print(
                    f"[WARNING] cycle detected at taxid {current}", file=sys.stderr
                )
                break
            visited.add(current)

            parent = self._parent_of.get(current)
            if parent is None:
                # taxid genuinely not in database
                if current == resolved:
                    # The original (or merged) taxid itself wasn't found
                    print(
                        f"[WARNING] taxid {current} not found", file=sys.stderr
                    )
                    result = LineageResult(
                        input_taxid=taxid,
                        resolved_taxid=0,
                        names=[],
                        ranks=[],
                        taxids=[],
                        status="not_found",
                    )
                    self._cache[taxid] = result
                    if taxid != resolved:
                        self._cache[resolved] = result
                    return result
                # Ancestor missing — stop here
                break

            names_rev.append(self._name_of.get(current, ""))
            ranks_rev.append(self._rank_of.get(current, "no rank"))
            taxids_rev.append(current)

            if parent == 1:
                # Root reached; include root (taxid 1) as well
                names_rev.append(self._name_of.get(1, "root"))
                ranks_rev.append(self._rank_of.get(1, "no rank"))
                taxids_rev.append(1)
                break

            current = parent

        # Reverse to top-down order
        names_rev.reverse()
        ranks_rev.reverse()
        taxids_rev.reverse()

        status = "merged" if taxid in self._merged else "ok"
        result = LineageResult(
            input_taxid=taxid,
            resolved_taxid=resolved,
            names=names_rev,
            ranks=ranks_rev,
            taxids=taxids_rev,
            status=status,
        )

        self._cache[taxid] = result
        if taxid != resolved:
            # Also cache under the resolved taxid for future hits
            resolved_result = LineageResult(
                input_taxid=resolved,
                resolved_taxid=resolved,
                names=names_rev,
                ranks=ranks_rev,
                taxids=taxids_rev,
                status="ok",
            )
            self._cache[resolved] = resolved_result

        return result

    def get_lineages_batch(
        self,
        taxids: List[int],
        show_progress: bool = True,
    ) -> List[LineageResult]:
        """Process a list of taxIDs, optionally showing a tqdm progress bar."""
        try:
            from tqdm import tqdm
            _has_tqdm = True
        except ImportError:
            _has_tqdm = False

        if show_progress and _has_tqdm:
            it = tqdm(taxids, desc="Querying lineages", unit="taxid", file=sys.stderr)
        else:
            it = taxids

        return [self.get_lineage(t) for t in it]
