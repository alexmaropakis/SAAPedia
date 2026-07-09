# SAAP Database

An interactive database for **single amino-acid polymorphisms / substituted peptides (SAAP)**
detected in mass-spectrometry proteomics data. Import a **CSV or Excel (.xlsx)** file, and the app
de-duplicates and normalizes it into a relational store, lets you browse / filter / sort, manage the
datasets and their source papers (DOIs), bulk-delete, and export any selection of SAAP as a
**UniProt-style FASTA** of the variant (substituted) peptides. It runs single-user locally or as a
shared **lab-wide** instance with password-protected editing.

- **Backend:** FastAPI + SQLite (SQLAlchemy)
- **Frontend:** React (build-free — served as vendored static files, no Node/npm required)
- **Runs with only Python installed.**

---

## Run it

```bash
./run.sh
```

Then open **http://127.0.0.1:8000** in your browser. The script creates a virtualenv, installs
dependencies, and starts the server (which serves both the API and the UI). Press **Ctrl+C** to stop.

To change host/port: `HOST=0.0.0.0 PORT=9000 ./run.sh`.

### Manual start (alternative)

```bash
cd backend
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload
```

Try it immediately by importing `sample_data/sample_saap.csv`.

---

## Using the app

Three tabs:

- **Browse** — one row per unique SAAP with rollups (datasets, digests, species, acquisition types,
  # observations, computed # datasets, best PEP, max positional probability, max evidence fragments).
  Search and filter by dataset / digest / species / acquisition / AAS / Ig / trypsin / min positional
  probability / max PEP; sort any column; click a peptide to see every underlying observation. Select
  rows to **export to FASTA** or **delete**. You can also export everything matching the current filter.
- **Import** — drop a **CSV, TSV, or Excel (.xlsx)** file, and optionally fill in a **dataset → DOI**
  map so each dataset links to its source paper. Dataset chips and observation rows then become
  clickable links.
- **Datasets** — see every dataset with its SAAP/observation counts, edit DOIs, open the linked
  paper, and (danger zone) **clear all data**.

### What counts as one SAAP, and one observation

A **SAAP** is uniquely identified by `(SAAP, BP, AAS)` — the substituted peptide, its base/reference
peptide, and the amino-acid substitution. Every line in an imported file becomes one **Observation**
linked to its SAAP. So a SAAP's **# observations** = the number of rows that resolved to it, and
**# datasets** is **computed** as the number of distinct `Dataset` values it appears in (any
`N Datasets` column in the file is ignored). Re-importing the same file is safe: an exact repeat of a
full row is detected by hash and not double-counted.

### Import rules (data cleaning)

On import, a peptide is **dropped** if it has **no UniProt ID** (`UniProt` accession). The import
summary reports how many were removed. (There is no RAAS in this schema.)

### Column mapping (tolerant)

Headers are matched case-, space-, and underscore-insensitively with aliases, so `MTP_seq`,
`mtp seq`, and `MTPSeq` all map to the same field. The import summary reports any **unmapped**
columns so nothing is silently dropped. To teach it new header spellings, edit
`backend/app/column_map.py` (`ALIASES`).

Recognized columns (canonical field ← accepted headers): `SAAP, BP, AAS, Dataset, TMT/Tissue,
Digest, Species, Data acquisition, PEP, Positional probability, N evidence fragments, UniProt,
Genes, RefProteins, missed_cleavage, AAS_at_peptide_terminus, greater_than_shared, Immunoglobulin,
Trypsin`. A `N Datasets` column is intentionally ignored (that value is computed).

### FASTA export

Each selected SAAP emits one entry whose **sequence is the variant peptide** (`MTP_seq`), under a
UniProt-style header. Default template:

```
>saap|{accession}|{entry_name} {protein} SAAP variant ({aa_sub}) OS={species} GN={gene} BP={bp_seq}
```

The template is editable at export time. Fields:
`{accession} {entry_name} {gene} {protein} {aa_sub} {sub_compact} {bp_seq} {mtp_seq} {species} {id}`.
`sub_compact` renders `"V to P"` as `V2P`.

---

## Project layout

```
SAAP_Database/
├── run.sh                     # one-command launcher
├── sample_data/
│   └── sample_saap.csv        # example input
└── backend/
    ├── requirements.txt
    └── app/
        ├── main.py            # FastAPI app + routes + static serving
        ├── database.py        # SQLite engine/session
        ├── models.py          # SAAP, Observation, DatasetInfo tables
        ├── column_map.py      # tolerant header → field mapping
        ├── ingest.py          # CSV parse, de-dup, cleaning, DOI upsert
        ├── crud.py            # list/detail/filter/delete queries + rollups
        ├── fasta.py           # variant-peptide FASTA generation
        ├── util.py            # DOI → URL helper
        ├── schemas.py         # request models
        └── static/            # build-free React SPA (index.html, app.js, styles.css, vendor/)
```

The SQLite database is written to `backend/saap.db`. Delete that file (or use Datasets → Clear all
data) for a clean slate.

---

## API reference

| Method | Path | Purpose |
|---|---|---|
| POST | `/api/upload` | Ingest a CSV (multipart `file`, optional `dataset_doi_map` JSON) |
| GET | `/api/saap` | List SAAP with filters/sort/pagination + rollups |
| GET | `/api/saap/{id}` | One SAAP with all its observations |
| POST | `/api/saap/delete` | Bulk delete: `{ids:[...]}` or `{all:true}` |
| GET | `/api/facets` | Distinct datasets & substitutions (filter options) |
| GET | `/api/datasets` | Datasets with DOI/link and counts |
| POST | `/api/datasets` | Save dataset DOIs: `{map:{dataset:doi}}` |
| GET | `/api/stats` | Overview counts |
| POST | `/api/export/fasta` | FASTA for `ids` or `filters` |

Plus `GET /api/config` → `{write_protected: bool}`. When a write password is configured, the
mutating endpoints (`/api/upload`, `/api/datasets`, `/api/saap/delete`) require an
`X-Write-Password` header; reads and export stay public.

Interactive API docs: **http://127.0.0.1:8000/docs**.

---

## Lab-wide / shared deployment

To run one shared instance for a lab (public browsing, password-protected editing), use
`serve.sh` with environment variables:

| Variable | Purpose |
|---|---|
| `SAAP_DB_PATH` | Path to the shared SQLite file (so all instances use one database) |
| `SAAP_WRITE_PASSWORD` | Password required to import / delete / edit DOIs (reads stay open) |
| `HOST` / `PORT` | Bind address (default `0.0.0.0`) and port (default `8000`) |

```bash
SAAP_DB_PATH=/work/yourlab/saap/saap.db \
SAAP_WRITE_PASSWORD='a-strong-password' \
./serve.sh
```

WAL mode is enabled for concurrent access, and additive schema changes (new columns) migrate
automatically on startup without wiping data. Step-by-step instructions for **Northeastern
Explorer (HPC)** — shared filesystem location, tmux/Slurm, and SSH-tunnel or reverse-proxy
access — are in **[DEPLOY_EXPLORER.md](DEPLOY_EXPLORER.md)**.

---

## Notes

- The frontend uses in-browser Babel (vendored under `static/vendor/`) so it needs no build step.
  If you later install Node, the same components port directly to a Vite project for a production build.
- Species is set at export time (defaults to *Homo sapiens*) since the source CSV has no species
  column; add a `species` column and it will map automatically.
