"""
fetcher.py — NCBI Datasets REST API calls, batching, and rate limiting.

Uses the public NCBI Datasets v2 REST API via `requests`.
BioSample environmental attributes are extracted from the same response
as the assembly metadata — no secondary API calls required.
"""

import sys
import time
from typing import List, Optional

import requests

BATCH_SIZE = 200
RATE_LIMIT_SLEEP = 0.2    # ~5 req/s; halved to 0.1 with API key
MAX_RETRIES = 5

_BASE_URL = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{accessions}/dataset_report"

# BioSample attribute names to try for each output field.
# The Datasets API returns some fields flattened at the biosample root level
# and all fields in the attributes array.  We check both.
_BS_ATTR_CANDIDATES = {
    "collection_date":  ["collection_date"],
    "geo_loc_name":     ["geo_loc_name", "geographic location"],
    "lat_lon":          ["lat_lon", "latitude and longitude"],
    "isolation_source": ["isolation_source", "isolation source"],
    "env_broad_scale":  ["env_broad_scale", "env_biome", "biome", "environment (biome)"],
    "env_local_scale":  ["env_local_scale", "env_feature", "feature", "environment (feature)"],
    "env_medium":       ["env_medium", "env_material", "material", "environment (material)"],
    "depth":            ["depth"],
    "altitude":         ["altitude", "elev", "elevation"],
    "metagenome_source":["metagenome_source"],
}


def _get(obj, *keys, default=""):
    """Safely walk a nested dict, returning default if any key is missing."""
    cur = obj
    for k in keys:
        if not isinstance(cur, dict):
            return default
        cur = cur.get(k)
        if cur is None:
            return default
    return cur if cur is not None else default


def _extract_biosample_fields(biosample: dict) -> dict:
    """
    Extract environmental/sampling attributes from an assembly_info.biosample object.

    The Datasets v2 API includes:
    - Top-level shortcuts for common fields (collection_date, geo_loc_name, etc.)
    - A full `attributes` array with {name, value} pairs for all submitted fields
    """
    # Build a lookup from the attributes array (lowercase keys for robust matching)
    attrs: dict = {}
    for item in (biosample.get("attributes") or []):
        name = (item.get("name") or "").lower().strip()
        val = item.get("value") or ""
        if name and val:
            attrs[name] = val

    result = {
        "biosample_title": _get(biosample, "description", "title"),
    }

    for field, candidates in _BS_ATTR_CANDIDATES.items():
        # 1. Try biosample top-level first (pre-harmonized by NCBI)
        val = biosample.get(field, "")
        # 2. Fall back to scanning attributes array with each candidate name
        if not val:
            for cand in candidates:
                val = attrs.get(cand.lower(), "")
                if val:
                    break
        result[field] = val or ""

    return result


def _extract_record(report: dict) -> dict:
    """Flatten one dataset_report entry into a plain dict."""
    asm      = report.get("assembly_info") or {}
    stats    = report.get("assembly_stats") or {}
    org      = report.get("organism") or {}
    ann      = report.get("annotation_info") or {}
    gene_counts = _get(ann, "stats", "gene_counts") or {}

    # BioSample is nested under assembly_info.biosample (not a flat field)
    biosample = asm.get("biosample") or {}

    record = {
        "accession":           report.get("accession", ""),
        "organism_name":       org.get("organism_name", ""),
        "tax_id":              org.get("tax_id", ""),
        "assembly_level":      asm.get("assembly_level", ""),
        "assembly_name":       asm.get("assembly_name", ""),
        "bioproject":          asm.get("bioproject_accession", ""),
        "biosample":           biosample.get("accession", ""),   # fixed path
        "submission_date":     asm.get("submission_date", ""),
        "genome_size":         stats.get("total_sequence_length", ""),
        "scaffold_n50":        stats.get("scaffold_n50", ""),
        "contig_n50":          stats.get("contig_n50", ""),
        "gc_percent":          stats.get("gc_percent", ""),
        "gene_count_total":    gene_counts.get("total", ""),
        "chromosome_count":    stats.get("number_of_chromosomes", ""),
        "annotation_provider": ann.get("provider", ""),
    }

    # Merge in BioSample environmental fields
    record.update(_extract_biosample_fields(biosample))

    return record


def fetch_metadata(
    accessions: List[str],
    api_key: Optional[str] = None,
    show_progress: bool = True,
    verbose: bool = False,
) -> List[dict]:
    """
    Fetch assembly metadata (including embedded BioSample attributes) from
    the NCBI Datasets v2 REST API.

    Parameters
    ----------
    accessions:     List of GCA_*/GCF_* accession strings.
    api_key:        Optional NCBI API key (doubles rate limit to 10 req/s).
    show_progress:  Show tqdm progress bar if available.
    verbose:        Print batch info to stderr.

    Returns
    -------
    List of flat metadata dicts, one per accession found.
    """
    try:
        from tqdm import tqdm
        _tqdm_available = True
    except ImportError:
        _tqdm_available = False

    headers = {"Accept": "application/json"}
    if api_key:
        headers["api-key"] = api_key
        sleep_between = 0.1
    else:
        sleep_between = RATE_LIMIT_SLEEP

    batches = [accessions[i: i + BATCH_SIZE] for i in range(0, len(accessions), BATCH_SIZE)]

    iterator = (
        tqdm(batches, desc="Fetching", unit="batch", file=sys.stderr)
        if show_progress and _tqdm_available
        else batches
    )

    results: List[dict] = []

    session = requests.Session()
    session.headers.update(headers)

    for batch in iterator:
        accession_str = ",".join(batch)
        url = _BASE_URL.format(accessions=accession_str)

        attempt = 0
        while True:
            resp = session.get(url, params={"page_size": BATCH_SIZE}, timeout=60)

            if resp.status_code == 429 and attempt < MAX_RETRIES:
                wait = 2 ** attempt
                if verbose:
                    print(f"  Rate limited (429). Retrying in {wait}s...", file=sys.stderr)
                time.sleep(wait)
                attempt += 1
                continue

            if resp.status_code != 200:
                print(
                    f"  WARNING: HTTP {resp.status_code} for batch {batch[:3]}...",
                    file=sys.stderr,
                )
                break

            data = resp.json()
            reports = data.get("reports") or []
            for report in reports:
                results.append(_extract_record(report))

            # Handle pagination
            next_token = data.get("next_page_token")
            while next_token:
                resp2 = session.get(
                    url,
                    params={"page_size": BATCH_SIZE, "page_token": next_token},
                    timeout=60,
                )
                if resp2.status_code != 200:
                    break
                data2 = resp2.json()
                for report in (data2.get("reports") or []):
                    results.append(_extract_record(report))
                next_token = data2.get("next_page_token")

            if verbose:
                print(
                    f"  Fetched batch of {len(batch)} → {len(reports)} records",
                    file=sys.stderr,
                )
            break

        time.sleep(sleep_between)

    return results
