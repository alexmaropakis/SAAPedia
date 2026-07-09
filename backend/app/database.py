"""SQLite database engine and session management.

The database path is configurable via the SAAP_DB_PATH environment variable so a
shared/lab deployment can point every instance at one file (e.g. on a network
share or server volume). WAL mode is enabled for better concurrent read/write.
"""
import os
from pathlib import Path

from sqlalchemy import create_engine, event, inspect, text
from sqlalchemy.orm import DeclarativeBase, sessionmaker

# Default: a file next to the backend package. Override with SAAP_DB_PATH.
_DEFAULT_DB = Path(__file__).resolve().parent.parent / "saap.db"
DB_PATH = Path(os.environ.get("SAAP_DB_PATH", str(_DEFAULT_DB)))
DATABASE_URL = f"sqlite:///{DB_PATH}"

engine = create_engine(
    DATABASE_URL,
    connect_args={"check_same_thread": False, "timeout": 30},
)


@event.listens_for(engine, "connect")
def _set_sqlite_pragmas(dbapi_conn, _record):
    cur = dbapi_conn.cursor()
    cur.execute("PRAGMA journal_mode=WAL;")     # concurrent reads during writes
    cur.execute("PRAGMA busy_timeout=30000;")   # wait rather than error under contention
    cur.execute("PRAGMA foreign_keys=ON;")
    cur.close()


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
    new nullable column (e.g. `digest`), this adds it in place without dropping
    data. Only handles additions — not renames/drops/type changes.
    """
    from . import models  # noqa: F401

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


def init_db():
    """Create tables if missing, then apply additive column migrations."""
    from . import models  # noqa: F401  (ensure models are registered)

    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    Base.metadata.create_all(bind=engine)
    _auto_add_columns()
