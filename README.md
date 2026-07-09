# SAAP Database

An interactive local database for **single amino-acid polymorphisms / substituted peptides (SAAP)**
from mass-spectrometry proteomics. Import a **CSV or Excel (.xlsx)** file, and the app de-duplicates
and normalizes it, lets you browse / filter / sort, attach source-paper DOIs to datasets, and export
any selection of SAAP as a **UniProt-style FASTA** of the variant (substituted) peptides.

Everyone runs their own copy on their own machine — download the repo, start it, open a browser.

- **Backend:** FastAPI + SQLite
- **Frontend:** React (build-free — no Node/npm required)
- **Only requirement:** Python 3.9+

---

## Get it running

Download the repo (green **Code ▸ Download ZIP** on GitHub, or `git clone`), then:

### macOS / Linux

```bash
cd SAAP_Database
./run.sh
```

The first run sets up a Python environment and installs dependencies (~15 s); after that it starts
in a second or two. When you see `running at: http://127.0.0.1:8000`, open that URL in your browser.
Stop it with **Ctrl+C**.

### Windows (or if `run.sh` won't run)

```bash
cd SAAP_Database\backend
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
uvicorn app.main:app
```

Then open **http://127.0.0.1:8000**.

---

## Using the app

Three tabs:

- **Browse** — one row per unique SAAP with rollups (datasets, digests, species, acquisition types,
  # observations, # datasets, best PEP, max positional probability, max evidence fragments). Search
  and filter by dataset / digest / species / acquisition / AAS / flags / min positional probability /
  max PEP; sort any column; click a peptide to see every underlying observation. Select rows to
  **export to FASTA** or **delete** — or export everything matching the current filter.
- **Import** — drop a CSV, TSV, or Excel file. Optionally fill in a **dataset → DOI** map so each
  dataset links to its source paper.
- **Datasets** — see each dataset with counts, edit DOIs, open the linked paper, or clear all data.

### What counts as one SAAP

A SAAP is uniquely identified by `(SAAP, BP, AAS)` — the substituted peptide, its base peptide, and
the amino-acid substitution. Every imported row becomes one **observation** linked to its SAAP.
**# datasets** is computed as the number of distinct **(Dataset, Species)** pairs the SAAP is actually
observed in (any `N Datasets` column in the file is ignored). Re-importing the same file is safe —
exact-duplicate rows are detected and not double-counted.

### Import rules

- A peptide is **dropped** if it has no **UniProt** accession.
- Columns are matched tolerantly (case/spacing/punctuation-insensitive), and the import summary
  reports any **unmapped** columns so nothing is silently mis-read.

Recognized columns: `SAAP, BP, AAS, Dataset, TMT/Tissue, Digest, Species, Data acquisition, PEP,
Positional probability, N evidence fragments, UniProt, RefProteins, Genes, missed_cleavage,
AAS_at_peptide_terminus, greater_than_shared, Immunoglobulin, Trypsin`. To add new header spellings,
edit `backend/app/column_map.py`.

### FASTA export

Each selected SAAP emits one entry whose sequence is the **variant peptide**, with a UniProt-style
header. **Species (`OS=`) is taken from the data.** Tip: to build a species-specific search database
(e.g. human only), filter **Species** first, then **Export all filtered**.

---

## Your data

All data lives in a single SQLite file at `backend/saap.db`, created on first run. It stays on your
machine. Delete that file (or use **Datasets ▸ Clear all data**) to start over. The database is
per-machine — importing on your copy does not affect anyone else's.

## Project layout

```
SAAP_Database/
├── run.sh                 # one-command launcher (macOS/Linux)
└── backend/
    ├── requirements.txt
    └── app/
        ├── main.py        # FastAPI app + routes + serves the UI
        ├── database.py    # SQLite engine/session
        ├── models.py      # SAAP, Observation, DatasetInfo
        ├── column_map.py  # tolerant header → field mapping
        ├── ingest.py      # parse, de-dup, clean
        ├── crud.py        # queries + rollups
        ├── fasta.py       # variant-peptide FASTA
        └── static/        # build-free React app (index.html, app.js, styles.css, vendor/)
```

Interactive API docs are at **http://127.0.0.1:8000/docs** while the app is running.
