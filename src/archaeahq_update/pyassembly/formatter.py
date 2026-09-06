"""
formatter.py — TSV output formatting and header constant.
"""

HEADER = "\t".join([
    # Assembly fields (15)
    "assembly_accession",
    "organism_name",
    "tax_id",
    "assembly_level",
    "assembly_name",
    "bioproject_accession",
    "biosample_accession",
    "submission_date",
    "genome_size",
    "scaffold_n50",
    "contig_n50",
    "gc_percent",
    "gene_count_total",
    "chromosome_count",
    "annotation_provider",
    # BioSample environmental fields (11)
    "biosample_title",
    "collection_date",
    "geo_loc_name",
    "lat_lon",
    "isolation_source",
    "env_broad_scale",
    "env_local_scale",
    "env_medium",
    "depth",
    "altitude",
    "metagenome_source",
])

# Maps internal dict key → output column (in order)
_FIELDS = [
    # Assembly
    ("accession",           "assembly_accession"),
    ("organism_name",       "organism_name"),
    ("tax_id",              "tax_id"),
    ("assembly_level",      "assembly_level"),
    ("assembly_name",       "assembly_name"),
    ("bioproject",          "bioproject_accession"),
    ("biosample",           "biosample_accession"),
    ("submission_date",     "submission_date"),
    ("genome_size",         "genome_size"),
    ("scaffold_n50",        "scaffold_n50"),
    ("contig_n50",          "contig_n50"),
    ("gc_percent",          "gc_percent"),
    ("gene_count_total",    "gene_count_total"),
    ("chromosome_count",    "chromosome_count"),
    ("annotation_provider", "annotation_provider"),
    # BioSample environmental
    ("biosample_title",     "biosample_title"),
    ("collection_date",     "collection_date"),
    ("geo_loc_name",        "geo_loc_name"),
    ("lat_lon",             "lat_lon"),
    ("isolation_source",    "isolation_source"),
    ("env_broad_scale",     "env_broad_scale"),
    ("env_local_scale",     "env_local_scale"),
    ("env_medium",          "env_medium"),
    ("depth",               "depth"),
    ("altitude",            "altitude"),
    ("metagenome_source",   "metagenome_source"),
]


def format_record(record: dict) -> str:
    """Return a single tab-delimited data row; empty string for missing fields."""
    values = []
    for key, _ in _FIELDS:
        val = record.get(key, "")
        if val is None:
            val = ""
        values.append(str(val))
    return "\t".join(values)
