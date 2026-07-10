"""FastAPI application: ingestion, querying, and FASTA export for the SAAP DB.

Serves the build-free React single-page app from ./static as well.
"""
from __future__ import annotations

import io
import json
from pathlib import Path

from fastapi import Depends, FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import StreamingResponse
from fastapi.staticfiles import StaticFiles
from sqlalchemy.orm import Session

from . import crud
from .database import get_db, init_db
from .fasta import DEFAULT_HEADER, generate_fasta
from .ingest import ingest_file
from .models import SAAP, Observation
from .schemas import ExportRequest

STATIC_DIR = Path(__file__).resolve().parent / "static"

app = FastAPI(title="SAAP Database", version="1.0.0")


@app.on_event("startup")
def _startup():
    init_db()


# ----------------------------- API: ingestion -----------------------------
@app.post("/api/upload")
async def upload(
    file: UploadFile = File(...),
    dataset_doi_map: str | None = Form(None),
    db: Session = Depends(get_db),
):
    if not file.filename or not file.filename.lower().endswith((".csv", ".tsv", ".txt", ".xlsx")):
        raise HTTPException(400, "Please upload a .csv, .tsv, or .xlsx file.")
    doi_map = None
    if dataset_doi_map:
        try:
            parsed = json.loads(dataset_doi_map)
            if isinstance(parsed, dict):
                doi_map = {str(k): str(v) for k, v in parsed.items()}
        except json.JSONDecodeError:
            raise HTTPException(400, "dataset_doi_map must be valid JSON.")
    raw = await file.read()
    try:
        result = ingest_file(db, raw, file.filename, dataset_doi_map=doi_map)
    except Exception as exc:  # surface parse errors to the UI
        raise HTTPException(400, f"Failed to ingest CSV: {exc}") from exc
    return result.as_dict()


@app.get("/api/datasets")
def datasets(db: Session = Depends(get_db)):
    return crud.list_datasets(db)


@app.post("/api/datasets")
def save_dataset_dois(payload: dict, db: Session = Depends(get_db)):
    mapping = payload.get("map") if isinstance(payload, dict) else None
    if not isinstance(mapping, dict):
        raise HTTPException(400, "Expected {\"map\": {dataset: doi}}.")
    saved = crud.upsert_dataset_dois(db, mapping)
    return {"saved": saved}


@app.post("/api/saap/delete")
def delete_saap(payload: dict, db: Session = Depends(get_db)):
    """Bulk delete. Body: {ids:[...]} to delete specific SAAP, or {all:true} to
    wipe every SAAP (and their observations)."""
    ids = payload.get("ids") if isinstance(payload, dict) else None
    wipe_all = bool(payload.get("all")) if isinstance(payload, dict) else False
    if not wipe_all and not ids:
        raise HTTPException(400, "Provide ids or all=true.")
    deleted = crud.delete_saap(db, ids=ids, wipe_all=wipe_all)
    return {"deleted": deleted}


# ------------------------------ API: querying ------------------------------
def _parse_bool(v: str | None):
    if v is None or v == "":
        return None
    return v.lower() in {"true", "1", "yes", "y"}


@app.get("/api/saap")
def list_saap(
    q: str | None = None,
    dataset: str | None = None,
    digest: str | None = None,
    species: str | None = None,
    acquisition_type: str | None = None,
    aa_sub: str | None = None,
    immunoglobulin: str | None = None,
    trypsin: str | None = None,
    missed_cleavage: str | None = None,
    aas_at_peptide_terminus: str | None = None,
    greater_than_shared: str | None = None,
    min_pos_prob: float | None = None,
    max_pep: float | None = None,
    sort: str = "n_observations",
    order: str = "desc",
    page: int = 1,
    page_size: int = Query(50, ge=1, le=500),
    db: Session = Depends(get_db),
):
    return crud.list_saap(
        db, q=q, dataset=dataset, digest=digest, species=species,
        acquisition_type=acquisition_type, aa_sub=aa_sub,
        immunoglobulin=_parse_bool(immunoglobulin), trypsin=_parse_bool(trypsin),
        missed_cleavage=_parse_bool(missed_cleavage),
        aas_at_peptide_terminus=_parse_bool(aas_at_peptide_terminus),
        greater_than_shared=_parse_bool(greater_than_shared),
        min_pos_prob=min_pos_prob, max_pep=max_pep,
        sort=sort, order=order, page=page, page_size=page_size,
    )


@app.get("/api/saap/{saap_id}")
def saap_detail(saap_id: int, db: Session = Depends(get_db)):
    result = crud.get_saap_detail(db, saap_id)
    if result is None:
        raise HTTPException(404, "SAAP not found")
    saap, observations = result
    return {
        "saap": {
            "id": saap.id, "mtp_seq": saap.mtp_seq, "bp_seq": saap.bp_seq,
            "aa_sub": saap.aa_sub, "source_accession": saap.source_accession,
            "source_gene": saap.source_gene, "ref_proteins": saap.ref_proteins,
            "immunoglobulin": saap.immunoglobulin, "trypsin": saap.trypsin,
            "missed_cleavage": saap.missed_cleavage,
            "aas_at_peptide_terminus": saap.aas_at_peptide_terminus,
            "greater_than_shared": saap.greater_than_shared,
        },
        "observations": [_obs_dict(o) for o in observations],
    }


def _obs_dict(o: Observation) -> dict:
    return {c.name: getattr(o, c.name) for c in Observation.__table__.columns}


@app.get("/api/facets")
def facets(db: Session = Depends(get_db)):
    return crud.distinct_values(db)


@app.get("/api/stats")
def stats(db: Session = Depends(get_db)):
    return crud.stats(db)


# ------------------------------ API: export --------------------------------
@app.post("/api/export/fasta")
async def export_fasta(
    payload: str = Form(...),
    reference: UploadFile | None = File(None),
    db: Session = Depends(get_db),
):
    """Build a FASTA. `payload` is the JSON ExportRequest; `reference` is an
    optional proteome FASTA to append (and decoy alongside the SAAP entries)."""
    try:
        req = ExportRequest(**json.loads(payload))
    except (json.JSONDecodeError, TypeError, ValueError) as exc:
        raise HTTPException(400, f"Bad export payload: {exc}") from exc

    filters = _clean_filters(req.filters) if req.filters else None
    saaps = crud.get_saap_for_export(db, req.ids, filters)

    ref_text = None
    if reference is not None and reference.filename:
        ref_text = (await reference.read()).decode("utf-8", errors="replace")

    if not saaps and not ref_text:
        raise HTTPException(400, "Nothing to export — no SAAP selected and no reference file.")

    saap_ids = [s.id for s in saaps]

    # If scoped to one species (via filter or override), stamp every SAAP header
    # with that species; otherwise use each SAAP's own species.
    override_species = (filters or {}).get("species") or (req.species or "")
    species_map = None if override_species else crud.species_by_saap(db, saap_ids)

    # The plex/pool token comes from the Dataset column: one dataset filtered ->
    # that dataset for all; otherwise each SAAP's own dataset(s).
    override_dataset = (filters or {}).get("dataset")
    token_map = None if override_dataset else crud.datasets_by_saap(db, saap_ids)
    default_token = override_dataset or ""

    fasta = generate_fasta(
        saaps,
        species_by_id=species_map,
        default_species=override_species,
        token=default_token,
        token_by_id=token_map,
        include_decoys=bool(req.decoys),
        reference_fasta=ref_text,
        line_width=req.line_width or 60,
        header_template=req.header_template or DEFAULT_HEADER,
    )
    filename = f"saap_export_{len(saaps)}{'_withref' if ref_text else ''}{'_decoys' if req.decoys else ''}.fasta"
    return StreamingResponse(
        io.BytesIO(fasta.encode("utf-8")),
        media_type="text/x-fasta",
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


_ALLOWED_FILTERS = {"q", "dataset", "digest", "species", "acquisition_type", "aa_sub",
                    "immunoglobulin", "trypsin", "missed_cleavage",
                    "aas_at_peptide_terminus", "greater_than_shared",
                    "min_pos_prob", "max_pep"}
_BOOL_FILTERS = {"immunoglobulin", "trypsin", "missed_cleavage",
                 "aas_at_peptide_terminus", "greater_than_shared"}


def _clean_filters(filters: dict) -> dict:
    out = {}
    for k in _ALLOWED_FILTERS:
        if k in filters and filters[k] not in (None, ""):
            out[k] = filters[k]
    for k in _BOOL_FILTERS:
        if k in out:
            out[k] = _parse_bool(str(out[k]))
    for numk in ("min_pos_prob", "max_pep"):
        if numk in out:
            try:
                out[numk] = float(out[numk])
            except (TypeError, ValueError):
                out.pop(numk)
    return out


@app.get("/api/health")
def health():
    return {"status": "ok"}


# ------------------------- static single-page app --------------------------
# Mounted last so /api/* routes take precedence.
if STATIC_DIR.exists():
    app.mount("/", StaticFiles(directory=str(STATIC_DIR), html=True), name="static")
