"""Query helpers: list SAAP with per-SAAP rollups, detail, filters, stats."""
from __future__ import annotations

from sqlalchemy import Select, and_, distinct, exists, func, or_, select
from sqlalchemy.orm import Session

from .models import DatasetInfo, Observation, SAAP
from .util import doi_to_url

# Sort keys the API accepts.
SORTABLE = {
    "mtp_seq", "bp_seq", "aa_sub", "source_gene", "ref_proteins", "source_accession",
    "n_observations", "n_datasets", "best_saap_pep", "max_positional_probability",
    "max_evidence_fragments",
}

# Aggregate columns, in a fixed order, exposed on the subquery.
_AGG_NAMES = ["n_observations", "n_datasets", "datasets", "digests", "species",
              "acquisition_types", "best_saap_pep", "max_positional_probability",
              "max_evidence_fragments"]


def _aggregate_subquery():
    # A "dataset" for counting purposes is the (dataset, species) pair, so the
    # same dataset name observed in two species counts as two. Key is NULL when
    # dataset is absent, so those rows are ignored (matches COUNT DISTINCT).
    dataset_species_key = Observation.dataset.concat("\x1f").concat(
        func.coalesce(Observation.species, "")
    )
    return (
        select(
            Observation.saap_id.label("saap_id"),
            func.count(Observation.id).label("n_observations"),
            func.count(distinct(dataset_species_key)).label("n_datasets"),
            func.group_concat(distinct(Observation.dataset)).label("datasets"),
            func.group_concat(distinct(Observation.digest)).label("digests"),
            func.group_concat(distinct(Observation.species)).label("species"),
            func.group_concat(distinct(Observation.acquisition_type)).label("acquisition_types"),
            func.min(Observation.saap_pep).label("best_saap_pep"),
            func.max(Observation.positional_probability).label("max_positional_probability"),
            func.max(Observation.n_evidence_fragments).label("max_evidence_fragments"),
        )
        .group_by(Observation.saap_id)
        .subquery()
    )


def _apply_filters(stmt: Select, agg, *, q=None, dataset=None, digest=None, species=None,
                   acquisition_type=None, aa_sub=None, immunoglobulin=None, trypsin=None,
                   missed_cleavage=None, aas_at_peptide_terminus=None, greater_than_shared=None,
                   min_pos_prob=None, max_pep=None) -> Select:
    if q:
        like = f"%{q}%"
        stmt = stmt.where(or_(
            SAAP.mtp_seq.ilike(like),
            SAAP.bp_seq.ilike(like),
            SAAP.aa_sub.ilike(like),
            SAAP.source_gene.ilike(like),
            SAAP.ref_proteins.ilike(like),
            SAAP.source_accession.ilike(like),
        ))

    def _obs_exists(col, value):
        return exists().where(and_(Observation.saap_id == SAAP.id, col == value))

    if dataset:
        stmt = stmt.where(_obs_exists(Observation.dataset, dataset))
    if digest:
        stmt = stmt.where(_obs_exists(Observation.digest, digest))
    if species:
        stmt = stmt.where(_obs_exists(Observation.species, species))
    if acquisition_type:
        stmt = stmt.where(_obs_exists(Observation.acquisition_type, acquisition_type))
    if aa_sub:
        stmt = stmt.where(SAAP.aa_sub == aa_sub)
    for flag_val, col in (
        (immunoglobulin, SAAP.immunoglobulin),
        (trypsin, SAAP.trypsin),
        (missed_cleavage, SAAP.missed_cleavage),
        (aas_at_peptide_terminus, SAAP.aas_at_peptide_terminus),
        (greater_than_shared, SAAP.greater_than_shared),
    ):
        if flag_val is not None:
            stmt = stmt.where(col.is_(flag_val))
    if min_pos_prob is not None:
        stmt = stmt.where(agg.c.max_positional_probability >= min_pos_prob)
    if max_pep is not None:
        stmt = stmt.where(agg.c.best_saap_pep <= max_pep)
    return stmt


def _split(concat):
    if not concat:
        return []
    return sorted(v for v in concat.split(",") if v)


def _row_to_dict(row) -> dict:
    saap: SAAP = row[0]
    agg = {name: row[i + 1] for i, name in enumerate(_AGG_NAMES)}
    return {
        "id": saap.id,
        "mtp_seq": saap.mtp_seq,
        "bp_seq": saap.bp_seq,
        "aa_sub": saap.aa_sub,
        "source_accession": saap.source_accession,
        "source_gene": saap.source_gene,
        "ref_proteins": saap.ref_proteins,
        "immunoglobulin": saap.immunoglobulin,
        "trypsin": saap.trypsin,
        "missed_cleavage": saap.missed_cleavage,
        "aas_at_peptide_terminus": saap.aas_at_peptide_terminus,
        "greater_than_shared": saap.greater_than_shared,
        "n_observations": agg["n_observations"],
        "n_datasets": agg["n_datasets"],
        "datasets": _split(agg["datasets"]),
        "digests": _split(agg["digests"]),
        "species": _split(agg["species"]),
        "acquisition_types": _split(agg["acquisition_types"]),
        "best_saap_pep": agg["best_saap_pep"],
        "max_positional_probability": agg["max_positional_probability"],
        "max_evidence_fragments": agg["max_evidence_fragments"],
    }


def list_saap(db: Session, *, sort="n_observations", order="desc", page=1, page_size=50, **filters):
    agg = _aggregate_subquery()
    agg_cols = [agg.c[name] for name in _AGG_NAMES]
    base = select(SAAP, *agg_cols).join(agg, agg.c.saap_id == SAAP.id)
    base = _apply_filters(base, agg, **filters)

    total = db.scalar(select(func.count()).select_from(base.subquery())) or 0

    if sort not in SORTABLE:
        sort = "n_observations"
    if sort in {"mtp_seq", "bp_seq", "aa_sub", "source_gene", "ref_proteins", "source_accession"}:
        sort_col = getattr(SAAP, sort)
    else:
        sort_col = agg.c[sort]
    sort_col = sort_col.desc() if order == "desc" else sort_col.asc()
    stmt = base.order_by(sort_col, SAAP.id.asc())

    page = max(page, 1)
    stmt = stmt.limit(page_size).offset((page - 1) * page_size)

    items = [_row_to_dict(row) for row in db.execute(stmt).all()]
    return {"items": items, "total": total, "page": page, "page_size": page_size}


def get_saap_detail(db: Session, saap_id: int):
    saap = db.get(SAAP, saap_id)
    if saap is None:
        return None
    obs = db.scalars(
        select(Observation).where(Observation.saap_id == saap_id).order_by(Observation.id)
    ).all()
    return saap, obs


def get_saap_for_export(db: Session, ids: list[int] | None, filters: dict | None):
    if ids:
        return db.scalars(select(SAAP).where(SAAP.id.in_(ids)).order_by(SAAP.id)).all()
    agg = _aggregate_subquery()
    stmt = select(SAAP).join(agg, agg.c.saap_id == SAAP.id)
    stmt = _apply_filters(stmt, agg, **(filters or {})).order_by(SAAP.id)
    return db.scalars(stmt).all()


def species_by_saap(db: Session, ids: list[int]) -> dict[int, str]:
    """Map saap_id -> species string (from the data). Multiple species for one
    SAAP are joined with '/'."""
    if not ids:
        return {}
    rows = db.execute(
        select(Observation.saap_id, func.group_concat(distinct(Observation.species)))
        .where(Observation.saap_id.in_(ids), Observation.species.is_not(None))
        .group_by(Observation.saap_id)
    ).all()
    out: dict[int, str] = {}
    for sid, concat in rows:
        vals = sorted({v for v in (concat or "").split(",") if v})
        if vals:
            out[sid] = "/".join(vals)
    return out


def distinct_values(db: Session):
    def col_values(col):
        return list(db.scalars(
            select(distinct(col)).where(col.is_not(None), col != "").order_by(col)
        ).all())
    return {
        "datasets": col_values(Observation.dataset),
        "digests": col_values(Observation.digest),
        "species": col_values(Observation.species),
        "acquisition_types": col_values(Observation.acquisition_type),
        "aa_subs": list(db.scalars(
            select(distinct(SAAP.aa_sub)).where(SAAP.aa_sub != "").order_by(SAAP.aa_sub)
        ).all()),
    }


def delete_saap(db: Session, *, ids=None, wipe_all: bool = False) -> int:
    """Delete SAAP (and their observations). Returns the number of SAAP removed."""
    from sqlalchemy import delete as sa_delete

    if wipe_all:
        n = db.scalar(select(func.count(SAAP.id))) or 0
        db.execute(sa_delete(Observation))
        db.execute(sa_delete(SAAP))
        db.commit()
        return n

    ids = [int(i) for i in (ids or [])]
    if not ids:
        return 0
    n = db.scalar(select(func.count(SAAP.id)).where(SAAP.id.in_(ids))) or 0
    db.execute(sa_delete(Observation).where(Observation.saap_id.in_(ids)))
    db.execute(sa_delete(SAAP).where(SAAP.id.in_(ids)))
    db.commit()
    return n


def list_datasets(db: Session):
    """All datasets known to the DB — present in observations and/or with a saved
    DOI — merged with their paper link and counts."""
    rows = db.execute(
        select(
            Observation.dataset,
            func.count(Observation.id),
            func.count(distinct(Observation.saap_id)),
        )
        .where(Observation.dataset.is_not(None))
        .group_by(Observation.dataset)
    ).all()
    counts = {name: (n_obs, n_saap) for name, n_obs, n_saap in rows}
    doi_by_name = {di.name: di.doi for di in db.scalars(select(DatasetInfo)).all()}

    out = []
    for name in sorted(set(counts) | set(doi_by_name)):
        n_obs, n_saap = counts.get(name, (0, 0))
        doi = doi_by_name.get(name)
        out.append({
            "name": name, "doi": doi, "url": doi_to_url(doi),
            "n_observations": n_obs, "n_saap": n_saap,
        })
    return out


def upsert_dataset_dois(db: Session, mapping: dict) -> int:
    from .ingest import upsert_dataset_dois as _upsert
    return _upsert(db, mapping)


def stats(db: Session):
    return {
        "n_saap": db.scalar(select(func.count(SAAP.id))) or 0,
        "n_observations": db.scalar(select(func.count(Observation.id))) or 0,
        "n_datasets": db.scalar(select(func.count(distinct(Observation.dataset)))) or 0,
        "n_genes": db.scalar(select(func.count(distinct(SAAP.source_gene)))) or 0,
    }
