# SAAPedia

SAAPedia is an interactive local database for substituted amino acid peptide (SAAP) sequences derived from mass spectrometry proteomics data. 

After importing a .csv, .tsv, or .xlsx file containing SAAP-base peptide (BP) pairs and related information, SAAPedia de-duplicates and allows you to browse, filter, and sort the results, and export any selection as UniProt-style FASTA files. 

Each copy is run locally on the terminal after downloading this repository. The repo ships with a populated `saap.db`. 

- **Backend:** FastAPI + SQLite
- **Frontend:** React
- **Requirement:** Python 3.9+

---

## Starting the app

Download the repo (green **Code ▸ Download ZIP** on GitHub, or `git clone`), then:

### macOS / Linux

```bash
cd SAAPedia
./run.sh
```

The first run creates a Python environment and installs dependencies (~15 s); later runs start in a
second or two. When you see `running at: http://127.0.0.1:8000`, open that URL. Stop with **Ctrl+C** on Windows or **Option+C** on Mac.

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

- **Browse** — one row per unique SAAP. Search, filter, sort, choose columns, select rows to export,
  annotate or delete. Immunoglobulin, trypsin and cleavage-site SAAPs are hidden by default (shown as
  removable filter chips). Click a row to open the SAAP.
- **Proteins** — every annotated protein ranked by SAAP count, with a site map of where its
  substitutions fall.
- **Protein / SAAP pages** — an NCBI-style protein map (SAAP sites, UniProt domains, regions,
  topology, PTMs, known variants, AlphaFold pLDDT and AlphaMissense tracks; drag to zoom), a
  GenPept-style sequence view, the AlphaFold 3D model with sites highlighted, and a BP/SAAP alignment
  in protein context.
- **Datasets** — trends and a per-dataset breakdown. Click a dataset to browse it.
- **Import** — drop a CSV, TSV or Excel file, then annotate.

### Substitution scores

| Score | Source | Needs network |
|---|---|---|
| BLOSUM62 | NCBI BLAST matrix | no |
| Grantham distance | Grantham 1974 (AAindex GRAR740104), classes per Li et al. 1984 | no |
| pLDDT at site | AlphaFold DB | yes |
| AlphaMissense | AlphaFold DB (human proteome only) | yes |
| Known variant | UniProt natural variants | yes |
| Genome-encoded | Exact match of the SAAP peptide against the reference proteome FASTAs in `backend/reference/` (I = L) | no |
| gnomAD AF | gnomAD v4 exomes + genomes: germline variants producing the same substitution, via an Ensembl transcript whose protein is identical to the annotated one (human only) | yes |

The gnomAD check and curation run in bulk with one resumable command (it can be stopped and re-run at
any time; fetched data is cached in `backend/.cache/`):

```bash
cd backend && nice -n 15 .venv/bin/python tools/enrich.py
```

Remote tracks are only shown when the remote sequence is identical to the annotated protein
sequence (e.g. non-canonical isoforms have no AlphaFold model), so positions never silently shift.
Responses are cached in `backend/.cache/` (override with `SAAP_CACHE_DIR`).

### Curation

Applied after every import and annotation run (`backend/app/curate.py`); the Import page shows the log.

- **Removed:** immunoglobulins (file flag, IG gene symbol, or protein name); contaminants (peptide in a
  cRAP protein — protease reagents, BSA, keratins — unless it is also the sample species' own sequence);
  genome-encoded SAAPs (substituted peptide verbatim in the reference proteome); SAAPs reproducible by a
  gnomAD variant with AF ≥ 0.0001 (likely polymorphisms); SAAPs matching a UniProt natural variant
  unless gnomAD confirms them rare (AF < 0.0001, or absent).
- **Logged, kept:** positional probability < 0.9 or missing (filter with *Min PosProb*).

### What counts as one SAAP?

A SAAP is uniquely identified by `(SAAP, BP, AAS)` — the substituted peptide, its base peptide, and
the amino-acid substitution. Every imported row becomes one **observation** linked to its SAAP.
**# datasets** is the number of distinct **(Dataset, Species)** pairs the SAAP is actually observed
in (any `N Datasets` column in the file is ignored). Re-importing the same file is safe — exact
duplicate rows are detected and not double-counted.

### Import rules

- A peptide with no **UniProt** accession is kept, not dropped — run **Annotate** afterward to
  resolve it from its gene symbol (or, failing that, its base peptide sequence).
- Columns are matched tolerantly (case / spacing / punctuation-insensitive). The import summary
  reports any **unmapped** columns, so nothing is silently mis-read.

Recognized columns: `SAAP, BP, AAS, Dataset, TMT/Tissue, Digest, Species, Data acquisition, PEP,
Positional probability, N evidence fragments, UniProt, RefProteins, Genes, missed_cleavage,
AAS_at_peptide_terminus, greater_than_shared, Immunoglobulin, Trypsin`. Add new header spellings in
`backend/app/column_map.py`.

## Exports
### FASTA
A major utility of SAAPedia is the ability to export all or a select number of SAAP sequences into a MaxQuant, FragPipe, or other quantitative proteomics engine-compatible UniProt-style FASTA.

Options:
1. Export SAAP sequences (peptide or substitution in whole protein) alone as their own FASTA
2. Export SAAP sequences appended to a reference proteome (uploaded from local desktop)
3. Export SAAP sequences only with reverse decoys
4. Export an entire FASTA containing reference proteome and SAAP sequences with all reverse decoys appended 

### .CSV
There is also the option to export the SAAP database as a .csv file for downstream analysis. 

---

## Where is my data?

All data lives in a single SQLite file at `backend/saap.db`. It stays on your machine and is
per-copy — importing on your machine doesn't affect anyone else's. Delete the file (or use
**Datasets ▸ Clear all data**) to start over. Set `SAAP_DB_PATH` to point at a different file. Local caches (`backend/.cache/`), reference proteomes
(`backend/reference/`) and personal utilities/backups (`utility/`) are not tracked in git.

## Project layout

```
SAAPedia/
├── run.sh                  # one-command launcher (macOS/Linux)
└── backend/
    ├── requirements.txt
    └── app/
        ├── main.py         # FastAPI routes; serves the UI
        ├── database.py     # SQLite engine/session, additive migrations
        ├── models.py       # SAAP, Observation
        ├── column_map.py   # tolerant header → field mapping
        ├── ingest.py       # parse, de-dup, persist
        ├── annotate.py     # UniProt / Ensembl annotation and positioning
        ├── cleavage.py     # digest-aware cleavage-site check
        ├── crud.py         # queries and rollups
        ├── investigate.py  # protein views; maps features and scores onto sites
        ├── external.py     # cached UniProt / AlphaFold clients
        ├── scoring.py      # BLOSUM62, Grantham
        ├── proteome.py     # local reference-proteome match
        ├── gnomad.py       # gnomAD population-variant check
        ├── curate.py       # removal / logging rules
        ├── data/contaminants.fasta  # cRAP (minus UPS standards) + Arg-C, Lys-C
        ├── fasta.py        # FASTA export
        ├── static/         # build-free React app (index.html, js/, styles.css, vendor/)
        └── ../tools/enrich.py     # resumable gnomAD check + curation
```

Interactive API docs live at <http://127.0.0.1:8000/docs> while the app is running.

## References

Chen, S., Francioli, L. C., Goodrich, J. K., Collins, R. L., Kanai, M., Wang, Q., Alföldi, J., Watts, N. A., Vittal, C., Gauthier, L. D., Poterba, T., Wilson, M. W., Tarasova, Y., Phu, W., Grant, R., Yohannes, M. T., Koenig, Z., Farjoun, Y., Banks, E., . . . Karczewski, K. J. (2024). A genomic mutational constraint map using variation in 76,156 human genomes. Nature, 625(7993), 92–100. https://doi.org/10.1038/s41586-023-06045-0

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T., Pritzel, A., Wong, L. H., Zielinski, M., Sargeant, T., Schneider, R. G., Senior, A. W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023). Accurate proteome-wide missense variant effect prediction with AlphaMissense. Science, 381(6664), Article eadg7492. https://doi.org/10.1126/science.adg7492

Dyer, S. C., Austine-Orimoloye, O., Azov, A. G., Barba, M., Barnes, I., Barrera-Enriquez, V. P., Becker, A., Bennett, R., Beracochea, M., Berry, A., Bhai, J., Bhurji, S. K., Boddu, S., Branco Lins, P. R., Brooks, L., Ramaraju, S. B., Campbell, L. I., Martinez, M. C., Charkhchi, M., . . . Yates, A. D. (2025). Ensembl 2025. Nucleic Acids Research, 53(D1), D948–D957. https://doi.org/10.1093/nar/gkae1071

The Global Proteome Machine Organization. (n.d.). cRAP protein sequences [Data set]. Retrieved October 8, 2026, from https://www.thegpm.org/crap/

Grantham, R. (1974). Amino acid difference formula to help explain protein evolution. Science, 185(4154), 862–864. https://doi.org/10.1126/science.185.4154.862

Guez, J., Goodrich, J. K., Moldovan, M. A., Chao, K. R., Kar, P., Panchal, R., Wilson, M. W., Laricchia, K. M., Rohlicek, G., Biba, D., Marten, D., He, Q., Darnowsky, P. W., Grant, R., Weisburd, B., Baxter, S. M., Nadeau, J., Lu, W., Jahl, S., . . . Karczewski, K. J. (2026). Integrating 730,947 exome sequences with clinical literature improves gene discovery [Preprint]. medRxiv. https://doi.org/10.64898/2026.03.23.26349081

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from protein blocks. Proceedings of the National Academy of Sciences, 89(22), 10915–10919. https://doi.org/10.1073/pnas.89.22.10915

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O., Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A., Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B., Nikolov, S., Jain, R., Adler, J., . . . Hassabis, D. (2021). Highly accurate protein structure prediction with AlphaFold. Nature, 596(7873), 583–589. https://doi.org/10.1038/s41586-021-03819-2

Karczewski, K. J., Francioli, L. C., Tiao, G., Cummings, B. B., Alföldi, J., Wang, Q., Collins, R. L., Laricchia, K. M., Ganna, A., Birnbaum, D. P., Gauthier, L. D., Brand, H., Solomonson, M., Watts, N. A., Rhodes, D., Singer-Berk, M., England, E. M., Seaby, E. G., Kosmicki, J. A., . . . MacArthur, D. G. (2020). The mutational constraint spectrum quantified from variation in 141,456 humans. Nature, 581(7809), 434–443. https://doi.org/10.1038/s41586-020-2308-7

Kawashima, S., Pokarowski, P., Pokarowska, M., Kolinski, A., Katayama, T., & Kanehisa, M. (2008). AAindex: Amino acid index database, progress report 2008. Nucleic Acids Research, 36(Database issue), D202–D205. https://doi.org/10.1093/nar/gkm998

Li, W.-H., Wu, C.-I., & Luo, C.-C. (1984). Nonrandomness of point mutation as reflected in nucleotide substitutions in pseudogenes and its evolutionary implications. Journal of Molecular Evolution, 21(1), 58–71. https://doi.org/10.1007/BF02100628

Rego, N., & Koes, D. (2015). 3Dmol.js: Molecular visualization with WebGL. Bioinformatics, 31(8), 1322–1324. https://doi.org/10.1093/bioinformatics/btu829

The UniProt Consortium. (2025). UniProt: The Universal Protein Knowledgebase in 2025. Nucleic Acids Research, 53(D1), D609–D617. https://doi.org/10.1093/nar/gkae1010

Varadi, M., Bertoni, D., Magana, P., Paramval, U., Pidruchna, I., Radhakrishnan, M., Tsenkov, M., Nair, S., Mirdita, M., Yeo, J., Kovalevskiy, O., Tunyasuvunakool, K., Laydon, A., Žídek, A., Tomlinson, H., Hariharan, D., Abrahamson, J., Green, T., Jumper, J., . . . Velankar, S. (2024). AlphaFold Protein Structure Database in 2024: Providing structure coverage for over 214 million protein sequences. Nucleic Acids Research, 52(D1), D368–D375. https://doi.org/10.1093/nar/gkad1011
