"""
SQLite database engine and session management

The database path is configurable via the SAAP_DB_PATH environment variable so a
shared/lab deployment can point every instance at one file (e.g. on a network
share or server volume). WAL mode is enabled for better concurrent read/write.
"""

import os
from pathlib import Path

from sqlalchemy import create_engine, event, inspect, text
from sqlalchemy.orm import DeclarativeBase, sessionmaker

from .cleavage import is_cleavage_position_any

# Default: a file next to the backend package. Override with SAAP_DB_PATH.
_DEFAULT_DB = Path(__file__).resolve().parent.parent / "saap.db"
DB_PATH = Path(os.environ.get("SAAP_DB_PATH", str(_DEFAULT_DB)))
DATABASE_URL = f"sqlite:///{DB_PATH}"

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False, "timeout": 30},
)


def _sql_at_cleavage_site(protein_sequence, position_in_protein, digests_csv, aa_sub):
    """SQLite UDF wrapper: 1/0/NULL so it composes in a WHERE clause."""
    result = is_cleavage_position_any(protein_sequence, position_in_protein, digests_csv, aa_sub)
    return None if result is None else int(result)


@event.listens_for(engine, "connect")
def _set_sqlite_pragmas(dbapi_conn, _record):
    cur = dbapi_conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL;")     # concurrent reads during writes
    cur.execute("PRAGMA busy_timeout=30000;")   # wait rather than error under contention
    cur.execute("PRAGMA foreign_keys=ON;")
    cur.close()
    # Exposes cleavage.is_cleavage_position_any() as saap_at_cleavage_site(...)
    # inside SQL, so it can be used directly in a WHERE clause (see crud.py).
    dbapi_conn.create_function("saap_at_cleavage_site", 4, _sql_at_cleavage_site)


SessionLocal = sessionmaker(bind=engine, autoflush=False, autocommit=False)


class Base(DeclarativeBase):
    pass


def get_db():
    """FastAPI dependency that yields a scoped session and always closes it."""
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()


def _auto_add_columns():
    """Additive migration: add any model columns missing from existing tables.

    SQLite's create_all never alters existing tables, so when the schema gains a
    new nullable column it is added in place (with its index) without dropping
    data. Only handles additions — not renames/drops/type changes.
    """
    from . import models  # noqa: F401  (registers the tables)

    inspector = inspect(engine)
    existing_tables = set(inspector.get_table_names())
    with engine.begin() as conn:
        for table in Base.metadata.sorted_tables:
            if table.name not in existing_tables:
                continue
            have = {c["name"] for c in inspector.get_columns(table.name)}
            for col in table.columns:
                if col.name not in have:
                    col_type = col.type.compile(dialect=engine.dialect)
                    conn.execute(text(
                        f'ALTER TABLE "{table.name}" ADD COLUMN "{col.name}" {col_type}'
                    ))
            for index in table.indexes:
                index.create(conn, checkfirst=True)


def init_db():
    """Create tables if missing, then apply additive column migrations."""
    from . import models  # noqa: F401  (ensure models are registered)

    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    Base.metadata.create_all(bind=engine)
    _auto_add_columns()
