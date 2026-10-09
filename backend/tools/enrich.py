"""
gnomAD population check, cross-species recurrence, then curation. Resumable: gnomAD/Ensembl responses
are cached in .cache/ and only unchecked SAAPs are queried, so it can be
stopped (or the machine shut down) and simply re-run.

    cd backend && nice -n 15 .venv/bin/python tools/enrich.py
"""
import logging
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from app import gnomad, orthology  # noqa: E402
from app.curate import prune  # noqa: E402
from app.database import SessionLocal, init_db  # noqa: E402


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s", datefmt="%H:%M:%S")
    for name in ("httpx", "urllib3"):
        logging.getLogger(name).setLevel(logging.WARNING)
    init_db()
    db = SessionLocal()
    for name, fn in (("gnomAD", lambda: gnomad.check_all(db)),
                     ("cross-species recurrence", lambda: orthology.check_all(db)),
                     ("curation", lambda: prune(db))):
        t = time.time()
        logging.info("START %s", name)
        logging.info("DONE %s %s (%ds)", name, fn(), time.time() - t)


if __name__ == "__main__":
    main()
