"""
FastAPI application: ingestion, querying, investigation and export for SAAPedia.

Serves the build-free React single-page app from ./static as well.
"""

from __future__ import annotations
import csv
import io
import json
from pathlib import Path
from typing import Mapping

from fastapi import Depends, FastAPI, File, Form, HTTPException, Query, Request, UploadFile
from fastapi.responses import PlainTextResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy.orm import Session

from . import annotate as annotate_mod
from . import crud, curate, external, investigate, proteome
from .database import get_db, init_db
from .fasta import BASE_HEADER, DEFAULT_HEADER, PROTEIN_HEADER, generate_fasta
from .ingest import ingest_file
from .models import SAAP, Observation
from .samples import PRIVATE
from .schemas import AnnotateRequest, ExportRequest

# Observation fields that identify a study; only served in private mode.
PRIVATE_FIELDS = {"dataset", "tmt_tissue", "source_file", "row_hash"}


def require_private():
    """Writes (import, annotate, delete) are local-only; the public site is read-only."""
    if not PRIVATE:
        raise HTTPException(403, "Read-only")

STATIC_DIR = Path(__file__).resolve().parent / "static"

app = FastAPI(title="SAAPedia", version="2.0.0")


@app.on_event("startup")
def _startup():
    init_db()


# --- filters -----------------------------------------------------------------
_STR_FILTERS = {"q", "tissue", "digest", "species", "acquisition_type", "aa_sub"} | ({"dataset"} if PRIVATE else set())
_BOOL_FILTERS = {"trypsin", "missed_cleavage", "aas_at_peptide_terminus",
                 "greater_than_shared", "at_cleavage_site", "in_gnomad", "cross_species"}
_NUM_FILTERS = {"min_pos_prob", "max_pep"}


def parse_filters(raw: Mapping | None) -> dict:
    """Whitelist and coerce filter values from query params or a JSON body."""
    out = {}
    for k, v in (raw or {}).items():
        if v is None or v == "":
            continue
        if k in _STR_FILTERS:
            out[k] = str(v)
        elif k in _BOOL_FILTERS:
            out[k] = str(v).lower() in {"true", "1", "yes", "y"}
        elif k in _NUM_FILTERS:
            try:
                out[k] = float(v)
            except (TypeError, ValueError):
                pass
    return out


def _csv_response(columns: list[tuple[str, str]], rows: list[dict], filename: str):
    buf = io.StringIO()
    writer = csv.writer(buf)
    writer.writerow([label for _, label in columns])
    for r in rows:
        writer.writerow(["; ".join(map(str, v)) if isinstance(v, list) else ("" if v is None else v)
                         for v in (r.get(k) for k, _ in columns)])
    return StreamingResponse(
        io.BytesIO(buf.getvalue().encode("utf-8-sig")), media_type="text/csv",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'})


# --- data --------------------------------------------------------------------
@app.post("/api/upload")
async def upload(file: UploadFile = File(...), db: Session = Depends(get_db), _=Depends(require_private)):
    if not file.filename or not file.filename.lower().endswith((".csv", ".tsv", ".txt", ".xlsx")):
        raise HTTPException(400, "Upload a .csv, .tsv, or .xlsx file.")
    try:
        result = ingest_file(db, await file.read(), file.filename)
    except Exception as exc:  # surface parse errors to the UI
        raise HTTPException(400, f"Could not import file: {exc}") from exc
    return result.as_dict()


@app.post("/api/saap/delete")
def delete_saap(payload: dict, db: Session = Depends(get_db), _=Depends(require_private)):
    """Body: {ids:[...]} for specific SAAP, or {all:true} to clear everything."""
    ids, wipe_all = payload.get("ids"), bool(payload.get("all"))
    if not wipe_all and not ids:
        raise HTTPException(400, "Provide ids or all=true.")
    return {"deleted": crud.delete_saap(db, ids=ids, wipe_all=wipe_all)}


@app.get("/api/stats")
def stats(db: Session = Depends(get_db)):
    return crud.stats(db)


@app.get("/api/facets")
def facets(db: Session = Depends(get_db)):
    return crud.distinct_values(db)


@app.get("/api/curation")
def curation(db: Session = Depends(get_db)):
    """Curation log: positional probability and gnomAD status (see curate.py)."""
    return curate.summary(db)


@app.get("/api/datasets")
def datasets(db: Session = Depends(get_db)):
    return crud.overview(db)


# --- SAAP --------------------------------------------------------------------
@app.get("/api/saap")
def list_saap(
    request: Request,
    sort: str = "n_observations",
    order: str = "desc",
    page: int = 1,
    page_size: int = Query(50, ge=1, le=500),
    db: Session = Depends(get_db),
):
    return crud.list_saap(db, sort=sort, order=order, page=page, page_size=page_size,
                          **parse_filters(request.query_params))


@app.get("/api/saap/{saap_id}")
def saap_detail(saap_id: int, db: Session = Depends(get_db)):
    saap, observations = crud.get_saap_detail(db, saap_id) or (None, None)
    if saap is None:
        raise HTTPException(404, "SAAP not found")
    fields = {c.name: getattr(saap, c.name) for c in SAAP.__table__.columns
              if c.name != "protein_sequence"}
    return {
        "saap": {
            **fields,
            **investigate.site_info(saap),
            "proteome_matches": proteome.describe(saap.proteome_hits),
            "at_cleavage_site": crud.saap_cleavage_flag(saap, observations),
            "max_positional_probability": max(
                (v for v in [saap.source_positional_probability,
                             *(o.positional_probability for o in observations)] if v is not None), default=None),
        },
        "observations": [{c.name: getattr(o, c.name) for c in Observation.__table__.columns
                          if PRIVATE or c.name not in PRIVATE_FIELDS} for o in observations],
    }


# --- proteins ----------------------------------------------------------------
@app.get("/api/proteins")
def list_proteins(
    q: str | None = None,
    sort: str = "n_saap",
    order: str = "desc",
    page: int = 1,
    page_size: int = Query(50, ge=1, le=500),
    db: Session = Depends(get_db),
):
    return crud.list_proteins(db, q=q, sort=sort, order=order, page=page, page_size=page_size)


def _check_accession(accession: str):
    if not external.valid_accession(accession):
        raise HTTPException(400, "Invalid UniProt accession")


@app.get("/api/proteins/{accession}")
def protein_detail(accession: str, db: Session = Depends(get_db)):
    _check_accession(accession)
    view = investigate.protein_view(db, accession)
    if view is None:
        raise HTTPException(404, "No annotated SAAP on this protein")
    return view


@app.get("/api/proteins/{accession}/annotations")
def protein_annotations(accession: str, db: Session = Depends(get_db)):
    """UniProt features, AlphaFold pLDDT and AlphaMissense, mapped to the
    cached sequence. Fetched once per protein and cached on disk."""
    _check_accession(accession)
    result = investigate.protein_annotations(db, accession)
    if result is None:
        raise HTTPException(404, "No annotated SAAP on this protein")
    return result


@app.get("/api/proteins/{accession}/function")
def protein_function(accession: str):
    """UniProt function, subcellular location, GO terms and Reactome pathways."""
    _check_accession(accession)
    try:
        return external.uniprot_function(accession)
    except external.ExternalError as exc:
        raise HTTPException(502, f"UniProt unavailable: {exc}") from exc


@app.get("/api/substitutions")
def substitutions(species: str | None = None, tissue: str | None = None, db: Session = Depends(get_db)):
    """Substitution matrix: distinct SAAP per reference -> substituted residue."""
    return investigate.substitution_matrix(db, species=species, tissue=tissue)


@app.get("/api/proteins/{accession}/structure.pdb", response_class=PlainTextResponse)
def protein_structure(accession: str):
    _check_accession(accession)
    try:
        pdb = external.alphafold_pdb(accession)
    except external.ExternalError as exc:
        raise HTTPException(502, f"AlphaFold DB unavailable: {exc}") from exc
    if pdb is None:
        raise HTTPException(404, "No AlphaFold model")
    return PlainTextResponse(pdb, media_type="chemical/x-pdb")


# --- annotation --------------------------------------------------------------
@app.post("/api/annotate")
def annotate(req: AnnotateRequest, db: Session = Depends(get_db), _=Depends(require_private)):
    """Resolve Ensembl IDs, protein sequence and substitution position from
    UniProt. Network failures are reported in the response, not raised."""
    return annotate_mod.annotate_saaps(db, ids=req.ids, overwrite=req.overwrite,
                                       limit=req.limit).as_dict()


@app.get("/api/annotate/status")
def annotate_status(db: Session = Depends(get_db)):
    return crud.annotation_status(db)


# --- export ------------------------------------------------------------------
@app.post("/api/export/fasta")
async def export_fasta(
    payload: str = Form(...),
    reference: UploadFile | None = File(None),
    db: Session = Depends(get_db),
):
    """`payload` is the JSON ExportRequest; `reference` an optional proteome
    FASTA appended after the SAAP entries (and decoyed with them)."""
    try:
        req = ExportRequest(**json.loads(payload))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise HTTPException(400, f"Bad export payload: {exc}") from exc

    filters = parse_filters(req.filters)
    saaps = crud.get_saap_for_export(db, req.ids, filters)
    ref_text = None
    if reference is not None and reference.filename:
        ref_text = (await reference.read()).decode("utf-8", errors="replace")
    if not saaps and not ref_text:
        raise HTTPException(400, "Nothing to export.")

    ids = [s.id for s in saaps]
    # A single filtered species/tissue (or dataset, privately) stamps every
    # header; otherwise each SAAP carries its own.
    species = filters.get("species") or req.species
    dataset = filters.get("dataset") or filters.get("tissue")
    protein_mode = (req.entry_mode or "").lower() == "protein"
    skipped: list[str] = []
    fasta = generate_fasta(
        saaps,
        species_by_id=None if species else crud.species_by_saap(db, ids),
        default_species=species,
        token=dataset or "",
        token_by_id=None if dataset else crud.tokens_by_saap(db, ids),
        include_decoys=req.decoys,
        include_base_peptides=req.base_peptides,
        entry_mode="protein" if protein_mode else "peptide",
        reference_fasta=ref_text,
        line_width=req.line_width or 60,
        header_template=req.header_template or (PROTEIN_HEADER if protein_mode else DEFAULT_HEADER),
        base_header_template=req.base_header_template or BASE_HEADER,
        skipped=skipped,
    )
    if protein_mode and not fasta.strip():
        raise HTTPException(400, "No SAAP could be written as a full-length protein — "
                                 "annotate first." + (f" Example: {skipped[0]}" if skipped else ""))
    filename = (f"saap_{'proteins' if protein_mode else 'peptides'}_{len(saaps)}"
                f"{'_bp' if req.base_peptides else ''}{'_ref' if ref_text else ''}"
                f"{'_decoys' if req.decoys else ''}.fasta")
    return StreamingResponse(
        io.BytesIO(fasta.encode("utf-8")), media_type="text/x-fasta",
        headers={"Content-Disposition": f'attachment; filename="{filename}"',
                 "X-SAAP-Skipped": str(len(skipped)),
                 "Access-Control-Expose-Headers": "X-SAAP-Skipped"})


_ROLLUP_COLUMNS = [
    ("mtp_seq", "SAAP"), ("bp_seq", "BP"), ("aa_sub", "AAS"),
    ("source_gene", "Gene"), ("source_accession", "UniProt"), ("ref_proteins", "RefProteins"),
    ("protein_accession", "Protein accession"), ("positions_all", "Position in protein"),
    ("n_observations", "N Observations"), ("n_tissues", "N Tissues"),
    ("tissues", "Tissues / cell types"), *([("datasets", "Datasets")] if PRIVATE else []),
    ("digests", "Digests"), ("species", "Species"),
    ("acquisition_types", "Data acquisition"),
    ("best_saap_pep", "Best PEP"), ("max_positional_probability", "Max positional probability"),
    ("max_evidence_fragments", "Max evidence fragments"),
    ("trypsin", "Trypsin"),
    ("missed_cleavage", "Missed cleavage"),
    ("aas_at_peptide_terminus", "AAS at peptide terminus"),
    ("greater_than_shared", "Greater than shared"),
    ("at_cleavage_site", "AAS at cleavage site (computed)"),
    ("proteome_hits", "Encoded by (reference proteome)"),
    ("gnomad_status", "gnomAD"), ("gnomad_af", "gnomAD AF (same substitution)"),
    ("known_variant", "UniProt natural variant (same substitution)"),
    ("cross_species", "Recurs in other species"), ("cross_species_detail", "Other-species match"),
]

_PAIR_COLUMNS = [
    ("saap", "SAAP"), ("bp", "BP"), ("substitution", "Substitution"), ("swap", "Swap (BP>SAAP)"),
    ("position_in_protein", "Position in protein"), ("peptide_start", "Peptide start"),
    ("n_positions", "Occurrences"), ("ensembl_gene", "Ensembl gene ID"),
    ("ensembl_transcript", "Ensembl transcript ID"), ("ensembl_protein", "Ensembl protein ID"),
    ("gene", "Gene"), ("protein_accession", "Protein accession"),
    ("protein_description", "Protein description"), ("protein_length", "Protein length"),
    ("annotation_source", "Annotation source"),
]


@app.post("/api/export/csv")
def export_csv(req: ExportRequest, db: Session = Depends(get_db)):
    rows = crud.get_rollups_for_export(db, req.ids, parse_filters(req.filters))
    if not rows:
        raise HTTPException(400, "Nothing to export.")
    return _csv_response(_ROLLUP_COLUMNS, rows, f"saap_export_{len(rows)}.csv")


@app.post("/api/export/pairs")
def export_pairs(req: ExportRequest, db: Session = Depends(get_db)):
    rows = crud.get_pairs_for_export(db, req.ids, parse_filters(req.filters))
    if not rows:
        raise HTTPException(400, "Nothing to export.")
    return _csv_response(_PAIR_COLUMNS, rows, f"saap_bp_pairs_{len(rows)}.csv")


@app.get("/api/health")
def health():
    return {"status": "ok", "version": app.version}


# --- single-page app ---------------------------------------------------------
class _NoCacheStaticFiles(StaticFiles):
    """app.js is transpiled in-browser without a cache-busting query, so force
    revalidation on every load or browsers keep running a stale copy."""

    def is_not_modified(self, response_headers, request_headers) -> bool:
        return False

    async def get_response(self, path: str, scope):
        response = await super().get_response(path, scope)
        response.headers["Cache-Control"] = "no-cache, no-store, must-revalidate"
        return response


if STATIC_DIR.exists():
    app.mount("/", _NoCacheStaticFiles(directory=str(STATIC_DIR), html=True), name="static")
