# SAAPedia: An interactive local database for **substituted amino acid peptide (SAAP)** sequences derived using mass spectrometry proteomics-based methods.

Import a **CSV, TSV, or Excel (.xlsx)** file; SAAPedia de-duplicates and normalizes it, then lets you
browse, filter, and sort the results, attach source-paper DOIs to datasets, and export any selection
as a **UniProt-style FASTA** of the variant peptides.

Everyone runs their own copy locally — download the repo, start it, open a browser. The repo ships
with a populated `saap.db`, so there's data to explore right away.

- **Backend:** FastAPI + SQLite
- **Frontend:** React
- **Requirement:** Python 3.9+

Find a short introduction to SAAPedia [here](https://docs.google.com/presentation/d/1ZOMBIVukzxbd0Y72HJP_PKo2ZN0wGaXiokHkLZKiyy4/edit?usp=sharing).

---

## Starting the app

Download the repo (green **Code ▸ Download ZIP** on GitHub, or `git clone`), then:

### macOS / Linux

```bash
cd SAAPedia
./run.sh
```

The first run creates a Python environment and installs dependencies (~15 s); later runs start in a
second or two. When you see `running at: http://127.0.0.1:8000`, open that URL. Stop with **Ctrl+C**.

### Windows (or if `run.sh` won't run)

```bash
cd SAAPedia\backend
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
uvicorn app.main:app
```

Then open **http://127.0.0.1:8000**.

---

## Using the app

- **Browse** — one row per unique SAAP with rollups (observations, datasets, species, best PEP, max
  positional probability, and more). Search, filter (dataset, digest, species, acquisition, AAS,
  flags, thresholds), and sort any column. Click a peptide to see every underlying observation.
  Select rows to **export to FASTA** or **delete** — or export everything matching the current filter.
- **Import** — drop a CSV, TSV, or Excel file. Optionally add a **dataset → DOI** map so each dataset
  links to its source paper.
- **Datasets** — view each dataset with counts, edit DOIs, open the linked paper, or clear all data.

### What counts as one SAAP?

A SAAP is uniquely identified by `(SAAP, BP, AAS)` — the substituted peptide, its base peptide, and
the amino-acid substitution. Every imported row becomes one **observation** linked to its SAAP.
**# datasets** is the number of distinct **(Dataset, Species)** pairs the SAAP is actually observed
in (any `N Datasets` column in the file is ignored). Re-importing the same file is safe — exact
duplicate rows are detected and not double-counted.

### Import rules

- A peptide is **dropped** if it has no **UniProt** accession.
- Columns are matched tolerantly (case / spacing / punctuation-insensitive). The import summary
  reports any **unmapped** columns, so nothing is silently mis-read.

Recognized columns: `SAAP, BP, AAS, Dataset, TMT/Tissue, Digest, Species, Data acquisition, PEP,
Positional probability, N evidence fragments, UniProt, RefProteins, Genes, missed_cleavage,
AAS_at_peptide_terminus, greater_than_shared, Immunoglobulin, Trypsin`. Add new header spellings in
`backend/app/column_map.py`.

### FASTA export

Each selected SAAP emits one entry — the **variant peptide** sequence with a UniProt-style header,
including species (`OS=`) from the data. To build a species-specific search database (e.g. human
only), filter **Species** first, then **Export all filtered**.

---

## Where is my data?

All data lives in a single SQLite file at `backend/saap.db`. It stays on your machine and is
per-copy — importing on your machine doesn't affect anyone else's. Delete the file (or use
**Datasets ▸ Clear all data**) to start over.

## Project layout

```
SAAPedia/
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

Interactive API docs live at <http://127.0.0.1:8000/docs> while the app is running.
