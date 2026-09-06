"""
db.py — SQLite cache schema and read/write helpers.
"""

import sqlite3
import pathlib
from datetime import datetime, timezone
from typing import Dict, List

SCHEMA = """
CREATE TABLE IF NOT EXISTS assemblies (
    accession           TEXT PRIMARY KEY,
    organism_name       TEXT,
    tax_id              INTEGER,
    assembly_level      TEXT,
    assembly_name       TEXT,
    bioproject          TEXT,
    biosample           TEXT,
    submission_date     TEXT,
    genome_size         INTEGER,
    scaffold_n50        INTEGER,
    contig_n50          INTEGER,
    gc_percent          REAL,
    gene_count_total    INTEGER,
    chromosome_count    INTEGER,
    annotation_provider TEXT,
    biosample_title     TEXT,
    collection_date     TEXT,
    geo_loc_name        TEXT,
    lat_lon             TEXT,
    isolation_source    TEXT,
    env_broad_scale     TEXT,
    env_local_scale     TEXT,
    env_medium          TEXT,
    depth               TEXT,
    altitude            TEXT,
    metagenome_source   TEXT,
    fetched_at          TEXT
);
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT
);
"""

# New columns added after the initial release — added to existing DBs via migration
_NEW_COLUMNS = [
    ("biosample_title",   "TEXT"),
    ("collection_date",   "TEXT"),
    ("geo_loc_name",      "TEXT"),
    ("lat_lon",           "TEXT"),
    ("isolation_source",  "TEXT"),
    ("env_broad_scale",   "TEXT"),
    ("env_local_scale",   "TEXT"),
    ("env_medium",        "TEXT"),
    ("depth",             "TEXT"),
    ("altitude",          "TEXT"),
    ("metagenome_source", "TEXT"),
]

_COLUMNS = [
    "accession", "organism_name", "tax_id", "assembly_level", "assembly_name",
    "bioproject", "biosample", "submission_date", "genome_size", "scaffold_n50",
    "contig_n50", "gc_percent", "gene_count_total", "chromosome_count",
    "annotation_provider",
    "biosample_title", "collection_date", "geo_loc_name", "lat_lon",
    "isolation_source", "env_broad_scale", "env_local_scale", "env_medium",
    "depth", "altitude", "metagenome_source",
    "fetched_at",
]

_SQLITE_MAX_VARS = 900   # SQLite limit is 999; stay safely under it


def open_db(data_dir: str) -> sqlite3.Connection:
    """Open (or create) the SQLite cache at data_dir/cache.db."""
    path = pathlib.Path(data_dir).expanduser()
    path.mkdir(parents=True, exist_ok=True)
    db_path = path / "cache.db"
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA)

    # Migrate existing databases: add any new columns that are missing
    for col_name, col_type in _NEW_COLUMNS:
        try:
            conn.execute(
                f"ALTER TABLE assemblies ADD COLUMN {col_name} {col_type} DEFAULT ''"
            )
        except sqlite3.OperationalError:
            pass  # Column already exists

    conn.commit()
    return conn


def get_cached(conn: sqlite3.Connection, accessions: List[str]) -> Dict[str, dict]:
    """
    Return a dict of {accession: record_dict} for accessions already in cache.
    Missing accessions are simply absent from the returned dict.
    Chunked to stay under SQLite's 999-variable IN-clause limit.
    """
    if not accessions:
        return {}
    result: Dict[str, dict] = {}
    for i in range(0, len(accessions), _SQLITE_MAX_VARS):
        chunk = accessions[i: i + _SQLITE_MAX_VARS]
        placeholders = ",".join("?" * len(chunk))
        rows = conn.execute(
            f"SELECT * FROM assemblies WHERE accession IN ({placeholders})",
            chunk,
        ).fetchall()
        result.update({row["accession"]: dict(row) for row in rows})
    return result


def store_results(conn: sqlite3.Connection, records: List[dict]) -> None:
    """Upsert a list of record dicts into the cache."""
    if not records:
        return
    now = datetime.now(timezone.utc).isoformat()
    placeholders = ",".join("?" * len(_COLUMNS))
    sql = (
        f"INSERT OR REPLACE INTO assemblies ({','.join(_COLUMNS)}) "
        f"VALUES ({placeholders})"
    )
    rows = []
    for rec in records:
        row = [rec.get(col, "") for col in _COLUMNS[:-1]]  # skip fetched_at
        row.append(now)
        rows.append(row)
    conn.executemany(sql, rows)
    conn.commit()
