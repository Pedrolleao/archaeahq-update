"""SQLite schema, index builder, and bulk data loader."""

import sqlite3
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Set, Tuple

DB_FILENAME = "taxonomy.db"
BATCH_SIZE = 10_000

SCHEMA = """
CREATE TABLE IF NOT EXISTS nodes (
    taxid  INTEGER PRIMARY KEY,
    parent INTEGER NOT NULL,
    rank   TEXT    NOT NULL
);
CREATE TABLE IF NOT EXISTS names (
    taxid INTEGER PRIMARY KEY,
    name  TEXT    NOT NULL
);
CREATE TABLE IF NOT EXISTS merged (
    old_taxid INTEGER PRIMARY KEY,
    new_taxid INTEGER NOT NULL
);
CREATE TABLE IF NOT EXISTS deleted (
    taxid INTEGER PRIMARY KEY
);
CREATE TABLE IF NOT EXISTS meta (
    key   TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
"""


def open_db(data_dir: Path) -> sqlite3.Connection:
    """Open (or create) the SQLite database with WAL mode."""
    db_path = data_dir / DB_FILENAME
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.execute("PRAGMA cache_size=-65536")   # 64 MB page cache
    conn.execute("PRAGMA temp_store=MEMORY")
    return conn


def db_is_valid(conn: sqlite3.Connection) -> bool:
    """Return True if the database was successfully built (has a build timestamp)."""
    try:
        row = conn.execute(
            "SELECT value FROM meta WHERE key='build_timestamp'"
        ).fetchone()
        return row is not None
    except sqlite3.OperationalError:
        return False


def _execute_schema(conn: sqlite3.Connection) -> None:
    conn.executescript(SCHEMA)
    conn.commit()


def _bulk_insert(conn: sqlite3.Connection, table: str, rows: list, batch_size: int = BATCH_SIZE) -> None:
    """Insert rows into table in batches within a single transaction."""
    if not rows:
        return
    placeholders = ",".join(["?"] * len(rows[0]))
    sql = f"INSERT OR REPLACE INTO {table} VALUES ({placeholders})"
    with conn:
        for i in range(0, len(rows), batch_size):
            conn.executemany(sql, rows[i : i + batch_size])


def build_index(data_dir: Path, verbose: bool = False) -> None:
    """
    Parse all .dmp files and load them into taxonomy.db.

    This is a one-time operation (~30–60 s) that builds the SQLite index.
    """
    from pytaxonkit.data import (
        parse_nodes_dmp,
        parse_names_dmp,
        parse_merged_dmp,
        parse_delnodes_dmp,
    )

    conn = open_db(data_dir)
    _execute_schema(conn)

    # Clear any previous partial data
    with conn:
        for table in ("nodes", "names", "merged", "deleted", "meta"):
            conn.execute(f"DELETE FROM {table}")

    def _load(label: str, table: str, rows: list) -> None:
        if verbose:
            print(f"[db] Loading {len(rows):,} rows into {table}...", file=sys.stderr)
        _bulk_insert(conn, table, rows)
        if verbose:
            print(f"[db] {table} done.", file=sys.stderr)

    # nodes.dmp
    if verbose:
        print(f"[db] Parsing nodes.dmp ...", file=sys.stderr)
    nodes_rows = list(parse_nodes_dmp(data_dir / "nodes.dmp"))
    _load("nodes", "nodes", nodes_rows)
    del nodes_rows

    # names.dmp
    if verbose:
        print(f"[db] Parsing names.dmp ...", file=sys.stderr)
    names_rows = list(parse_names_dmp(data_dir / "names.dmp"))
    _load("names", "names", names_rows)
    del names_rows

    # merged.dmp
    if verbose:
        print(f"[db] Parsing merged.dmp ...", file=sys.stderr)
    merged_rows = list(parse_merged_dmp(data_dir / "merged.dmp"))
    _load("merged", "merged", merged_rows)
    del merged_rows

    # delnodes.dmp
    if verbose:
        print(f"[db] Parsing delnodes.dmp ...", file=sys.stderr)
    del_rows = [(t,) for t in parse_delnodes_dmp(data_dir / "delnodes.dmp")]
    _load("deleted", "deleted", del_rows)
    del del_rows

    # Stamp build timestamp
    ts = datetime.now(tz=timezone.utc).isoformat()
    with conn:
        conn.execute(
            "INSERT OR REPLACE INTO meta VALUES (?, ?)",
            ("build_timestamp", ts),
        )

    conn.close()
    if verbose:
        print(f"[db] Index build complete at {ts}.", file=sys.stderr)


def load_all(conn: sqlite3.Connection) -> Tuple[
    Dict[int, int],   # parent_of:  taxid → parent_taxid
    Dict[int, str],   # rank_of:    taxid → rank
    Dict[int, str],   # name_of:    taxid → scientific name
    Dict[int, int],   # merged:     old_taxid → new_taxid
    Set[int],         # deleted:    set of deleted taxids
]:
    """Load all taxonomy data from SQLite into Python dicts (hot path)."""
    cur = conn.cursor()

    cur.execute("SELECT taxid, parent, rank FROM nodes")
    parent_of: Dict[int, int] = {}
    rank_of: Dict[int, str] = {}
    for taxid, parent, rank in cur.fetchall():
        parent_of[taxid] = parent
        rank_of[taxid] = rank

    cur.execute("SELECT taxid, name FROM names")
    name_of: Dict[int, str] = {row[0]: row[1] for row in cur.fetchall()}

    cur.execute("SELECT old_taxid, new_taxid FROM merged")
    merged: Dict[int, int] = {row[0]: row[1] for row in cur.fetchall()}

    cur.execute("SELECT taxid FROM deleted")
    deleted: Set[int] = {row[0] for row in cur.fetchall()}

    return parent_of, rank_of, name_of, merged, deleted
