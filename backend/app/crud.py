"""Query helpers: list SAAP with per-SAAP rollups, detail, filters, stats."""
from __future__ import annotations

from sqlalchemy import Select, and_, case, distinct, exists, func, or_, select
from sqlalchemy.orm import Session

from .annotate import NEEDS_ANNOTATION
from .cleavage import is_cleavage_position_any
from .curate import best_positional_probability
from .models import Observation, SAAP
from .samples import PRIVATE, organ

# Sort keys the API accepts.
SORTABLE = {
    "mtp_seq", "bp_seq", "aa_sub", "source_gene", "ref_proteins", "source_accession",
    "ensembl_gene", "ensembl_transcript", "ensembl_protein", "position_in_protein",
    "n_observations", "n_tissues", "best_saap_pep", "max_positional_probability",
    "max_evidence_fragments", "gnomad_af",
}

# Aggregate columns, in a fixed order, exposed on the subquery.
_AGG_NAMES = ["n_observations", "n_tissues", "tissues", "datasets", "digests", "species",
              "acquisition_types", "best_saap_pep", "max_positional_probability",
              "max_evidence_fragments"]


def _aggregate_subquery():
    # A tissue for counting purposes is the (tissue, species) pair, so mouse
    # and human lung count as two.
    tissue_species_key = Observation.tissue.concat("\x1f").concat(
        func.coalesce(Observation.species, "")
    )
    return (
        select(
            Observation.saap_id.label("saap_id"),
            func.count(Observation.id).label("n_observations"),
            func.count(distinct(tissue_species_key)).label("n_tissues"),
            func.group_concat(distinct(Observation.tissue)).label("tissues"),
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


def _apply_filters(stmt: Select, agg, *, q=None, dataset=None, tissue=None, digest=None, species=None,
                   acquisition_type=None, aa_sub=None, immunoglobulin=None, trypsin=None,
                   missed_cleavage=None, aas_at_peptide_terminus=None, greater_than_shared=None,
                   at_cleavage_site=None, in_gnomad=None, cross_species=None, min_pos_prob=None,
                   max_pep=None) -> Select:
    if q:
        like = f"%{q}%"
        stmt = stmt.where(or_(
            SAAP.mtp_seq.ilike(like),
            SAAP.bp_seq.ilike(like),
            SAAP.aa_sub.ilike(like),
            SAAP.source_gene.ilike(like),
            SAAP.ref_proteins.ilike(like),
            SAAP.source_accession.ilike(like),
            SAAP.ensembl_gene.ilike(like),
            SAAP.ensembl_transcript.ilike(like),
            SAAP.ensembl_protein.ilike(like),
        ))

    def _obs_exists(col, value):
        return exists().where(and_(Observation.saap_id == SAAP.id, col == value))

    if dataset:
        stmt = stmt.where(_obs_exists(Observation.dataset, dataset))
    if tissue:  # an organ ("Brain") matches all of its sub-sites ("Brain (cortex)")
        stmt = stmt.where(exists().where(and_(Observation.saap_id == SAAP.id, or_(
            Observation.tissue == tissue, Observation.tissue.like(tissue.replace("%", "") + " (%")))))
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
        if flag_val is True:
            stmt = stmt.where(col.is_(True))
        elif flag_val is False:
            # NULL means the source file never populated this column for that
            # row (common — most datasets only set some of these flags), not
            # "confirmed false". Treat it as "not flagged" rather than
            # excluding it, so filtering to No doesn't silently drop rows that
            # were simply never evaluated for this flag.
            stmt = stmt.where(or_(col.is_(False), col.is_(None)))
    if at_cleavage_site is not None:
        # Computed live from protein_sequence + position_in_protein + digest
        # (see cleavage.py) rather than a stored column — SQLite UDF registered
        # in database.py. NULL means "not evaluable" (no protein/position, or
        # an unrecognized digest); those rows are kept for both true and false
        # so an unresolved case is never silently dropped.
        expr = func.saap_at_cleavage_site(
            SAAP.protein_sequence, SAAP.position_in_protein, agg.c.digests, SAAP.aa_sub
        )
        stmt = stmt.where(expr == 1) if at_cleavage_site else stmt.where(or_(expr == 0, expr.is_(None)))
    if cross_species is not None:  # True: same substitution recurs in the other species
        stmt = stmt.where(SAAP.cross_species == "same") if cross_species else stmt.where(
            or_(SAAP.cross_species != "same", SAAP.cross_species.is_(None)))
    if in_gnomad is not None:
        stmt = stmt.where(SAAP.gnomad_status == ("present" if in_gnomad else "absent"))
    if min_pos_prob is not None:
        stmt = stmt.where(best_positional_probability(agg.c.max_positional_probability) >= min_pos_prob)
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
        "ensembl_gene": saap.ensembl_gene,
        "ensembl_transcript": saap.ensembl_transcript,
        "ensembl_protein": saap.ensembl_protein,
        "position_in_protein": saap.position_in_protein,
        "positions_all": saap.positions_all,
        "n_positions": saap.n_positions,
        "protein_description": saap.protein_description,
        "protein_accession": saap.protein_accession,
        "immunoglobulin": saap.immunoglobulin,
        "trypsin": saap.trypsin,
        "missed_cleavage": saap.missed_cleavage,
        "aas_at_peptide_terminus": saap.aas_at_peptide_terminus,
        "greater_than_shared": saap.greater_than_shared,
        "at_cleavage_site": is_cleavage_position_any(
            saap.protein_sequence, saap.position_in_protein, agg["digests"], saap.aa_sub
        ),
        "n_observations": agg["n_observations"],
        "n_tissues": agg["n_tissues"],
        "tissues": _split(agg["tissues"]),
        **({"datasets": _split(agg["datasets"])} if PRIVATE else {}),
        "digests": _split(agg["digests"]),
        "species": _split(agg["species"]),
        "acquisition_types": _split(agg["acquisition_types"]),
        "best_saap_pep": agg["best_saap_pep"],
        "max_positional_probability": max(
            (v for v in (agg["max_positional_probability"], saap.source_positional_probability) if v is not None),
            default=None),
        "proteome_hits": saap.proteome_hits,
        "gnomad_status": saap.gnomad_status,
        "cross_species": saap.cross_species,
        "cross_species_detail": saap.cross_species_detail,
        "gnomad_af": saap.gnomad_af,
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
    if sort in {"mtp_seq", "bp_seq", "aa_sub", "source_gene", "ref_proteins", "source_accession",
                "ensembl_gene", "ensembl_transcript", "ensembl_protein", "position_in_protein", "gnomad_af"}:
        sort_col = getattr(SAAP, sort)
    elif sort == "max_positional_probability":
        sort_col = best_positional_probability(agg.c[sort])
    else:
        sort_col = agg.c[sort]
    sort_col = sort_col.desc() if order == "desc" else sort_col.asc()
    stmt = base.order_by(sort_col, SAAP.id.asc())

    page = max(page, 1)
    stmt = stmt.limit(page_size).offset((page - 1) * page_size)

    items = [_row_to_dict(row) for row in db.execute(stmt).all()]
    return {"items": items, "total": total, "page": page, "page_size": page_size}


def saap_cleavage_flag(saap: SAAP, observations: list[Observation]) -> bool | None:
    """Same computed cleavage-site check as the Browse rollup (`at_cleavage_site`
    in `_row_to_dict`), for the single-SAAP detail view — evaluated against the
    digest(s) actually seen in `observations` rather than a pre-aggregated string."""
    digests_csv = ",".join(sorted({o.digest for o in observations if o.digest}))
    return is_cleavage_position_any(saap.protein_sequence, saap.position_in_protein, digests_csv, saap.aa_sub)


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


def get_rollups_for_export(db: Session, ids: list[int] | None, filters: dict | None):
    """Rollup rows (same shape as the Browse table) for CSV export — all matches,
    no pagination. Explicit ids take precedence over filters."""
    agg = _aggregate_subquery()
    agg_cols = [agg.c[name] for name in _AGG_NAMES]
    stmt = select(SAAP, *agg_cols).join(agg, agg.c.saap_id == SAAP.id)
    if ids:
        stmt = stmt.where(SAAP.id.in_(ids))
    else:
        stmt = _apply_filters(stmt, agg, **(filters or {}))
    stmt = stmt.order_by(SAAP.id)
    return [_row_to_dict(row) for row in db.execute(stmt).all()]


def _swap_notation(aa_sub: str | None, bp_seq: str | None, mtp_seq: str | None) -> str:
    """The substitution in BP>SAAP form, e.g. 'V>P'.

    Prefers the AAS column; falls back to comparing the base and substituted
    peptides when the two differ at exactly one residue.
    """
    from .annotate import parse_substitution, substitution_offset

    frm, to = parse_substitution(aa_sub)
    if frm and to:
        return f"{frm}>{to}"
    offset = substitution_offset(bp_seq, mtp_seq)
    if offset is not None:
        return f"{bp_seq[offset]}>{mtp_seq[offset]}"
    return ""


def _pair_row(saap: SAAP) -> dict:
    return {
        "saap": saap.mtp_seq,
        "bp": saap.bp_seq,
        "substitution": saap.aa_sub,
        "swap": _swap_notation(saap.aa_sub, saap.bp_seq, saap.mtp_seq),
        "position_in_protein": saap.positions_all or saap.position_in_protein,
        "peptide_start": saap.peptide_starts_all or saap.peptide_start,
        "n_positions": saap.n_positions,
        "ensembl_gene": saap.ensembl_gene,
        "ensembl_transcript": saap.ensembl_transcript,
        "ensembl_protein": saap.ensembl_protein,
        "gene": saap.source_gene,
        "protein_accession": saap.source_accession,
        "protein_description": saap.protein_description or saap.ref_proteins,
        "protein_length": saap.protein_length,
        "annotation_source": saap.annotation_source,
    }


def get_pairs_for_export(db: Session, ids: list[int] | None, filters: dict | None):
    """SAAP-BP pair rows (one per SAAP) for the pairs CSV export."""
    saaps = get_saap_for_export(db, ids, filters)
    return [_pair_row(s) for s in saaps]


def annotation_status(db: Session) -> dict:
    """Counts of how many SAAP carry Ensembl IDs / positions.

    `n_needs_annotation` uses the exact same condition as an "annotate new"
    run (`NEEDS_ANNOTATION`), so the count shown in the UI always matches how
    many rows that button will actually touch.
    """
    total = db.scalar(select(func.count(SAAP.id))) or 0
    with_position = db.scalar(
        select(func.count(SAAP.id)).where(SAAP.position_in_protein.is_not(None))
    ) or 0
    needs_annotation = db.scalar(
        select(func.count(SAAP.id)).where(NEEDS_ANNOTATION)
    ) or 0
    return {
        "n_saap": total,
        "n_with_position": with_position,
        "n_needs_annotation": needs_annotation,
    }


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


def tokens_by_saap(db: Session, ids: list[int]) -> dict[int, str]:
    """Map saap_id -> FASTA header token: its distinct datasets (private mode)
    or tissues (public) joined with '_'."""
    if not ids:
        return {}
    col = Observation.dataset if PRIVATE else Observation.tissue
    rows = db.execute(
        select(Observation.saap_id, func.group_concat(distinct(col)))
        .where(Observation.saap_id.in_(ids), col.is_not(None))
        .group_by(Observation.saap_id)
    ).all()
    out: dict[int, str] = {}
    for sid, concat in rows:
        vals = sorted({v for v in (concat or "").split(",") if v})
        if vals:
            out[sid] = "_".join(vals)
    return out


def distinct_values(db: Session):
    def col_values(col):
        return list(db.scalars(
            select(distinct(col)).where(col.is_not(None), col != "").order_by(col)
        ).all())
    return {
        **({"datasets": col_values(Observation.dataset)} if PRIVATE else {}),
        "tissues": sorted({t for v in col_values(Observation.tissue) for t in (v, organ(v))}),
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


def _ranked(db: Session, col) -> list[dict]:
    """Distinct-SAAP counts for every value of `col`, most common first."""
    n_saap = func.count(distinct(Observation.saap_id))
    rows = db.execute(
        select(col, n_saap)
        .where(col.is_not(None), col != "")
        .group_by(col)
        .order_by(n_saap.desc())
    ).all()
    return [{"label": label, "n": n} for label, n in rows]


def _group_overview(db: Session, col) -> list[dict]:
    """Per-group (tissue or dataset) x species summary: distinct SAAP, observations,
    annotation coverage (same "has a position" bar as `annotation_status`)."""
    annotated = case((SAAP.position_in_protein.is_not(None), Observation.saap_id))
    rows = db.execute(
        select(col, Observation.species, func.min(Observation.sample_type),
               func.count(Observation.id), func.count(distinct(Observation.saap_id)),
               func.count(distinct(annotated)),
               func.group_concat(distinct(Observation.digest)),
               func.group_concat(distinct(Observation.acquisition_type)))
        .join(SAAP, SAAP.id == Observation.saap_id)
        .where(col.is_not(None))
        .group_by(col, Observation.species)
        .order_by(Observation.species, func.count(distinct(Observation.saap_id)).desc())
    ).all()
    return [{"name": name, "species": species, "sample_type": kind, "n_observations": n_obs, "n_saap": n_saap,
             "n_annotated": n_ann, "digests": _split(dig), "acquisition_types": _split(acq)}
            for name, species, kind, n_obs, n_saap, n_ann, dig, acq in rows]


def overview(db: Session) -> dict:
    """Tissue / cell-type summary plus database-wide trends (dataset breakdown
    only in private mode)."""
    top_substitutions = db.execute(
        select(SAAP.aa_sub, func.count(SAAP.id))
        .where(SAAP.aa_sub.is_not(None), SAAP.aa_sub != "")
        .group_by(SAAP.aa_sub)
        .order_by(func.count(SAAP.id).desc())
        .limit(12)
    ).all()

    return {
        "tissues": _group_overview(db, Observation.tissue),
        **({"datasets": _group_overview(db, Observation.dataset)} if PRIVATE else {}),
        "tissue_distribution": _ranked(db, Observation.tissue),
        "species_distribution": _ranked(db, Observation.species),
        "digest_distribution": _ranked(db, Observation.digest),
        "acquisition_distribution": _ranked(db, Observation.acquisition_type),
        "top_substitutions": [{"label": a, "n": n} for a, n in top_substitutions],
    }


_PROTEIN_SORT = {"gene", "protein_accession", "length", "n_saap", "n_sites", "n_observations"}


def list_proteins(db: Session, *, q=None, sort="n_saap", order="desc", page=1, page_size=50):
    """One row per annotated protein (protein_accession), with its SAAP sites."""
    n_obs = (select(Observation.saap_id, func.count(Observation.id).label("n"))
             .group_by(Observation.saap_id).subquery())
    cols = {
        "protein_accession": SAAP.protein_accession,
        "gene": func.min(SAAP.source_gene),
        "description": func.min(SAAP.protein_description),
        "length": func.max(SAAP.protein_length),
        "n_saap": func.count(SAAP.id),
        "n_sites": func.count(distinct(SAAP.position_in_protein)),
        "n_observations": func.sum(n_obs.c.n),
        "positions": func.group_concat(SAAP.position_in_protein),
    }
    stmt = (select(*[c.label(k) for k, c in cols.items()])
            .join(n_obs, n_obs.c.saap_id == SAAP.id)
            .where(SAAP.protein_accession.is_not(None))
            .group_by(SAAP.protein_accession))
    if q:
        like = f"%{q}%"
        stmt = stmt.where(or_(SAAP.protein_accession.ilike(like), SAAP.source_gene.ilike(like),
                              SAAP.protein_description.ilike(like)))
    total = db.scalar(select(func.count()).select_from(stmt.subquery())) or 0
    key = sort if sort in _PROTEIN_SORT else "n_saap"
    sort_col = cols[key].desc() if order == "desc" else cols[key].asc()
    stmt = stmt.order_by(sort_col, SAAP.protein_accession).limit(page_size).offset((max(page, 1) - 1) * page_size)
    items = []
    for r in db.execute(stmt).mappings():
        item = dict(r)
        item["gene"] = (item["gene"] or "").split(";")[0] or None
        item["positions"] = sorted({int(p) for p in (r["positions"] or "").split(",") if p})
        items.append(item)
    return {"items": items, "total": total, "page": page, "page_size": page_size}


def protein_saaps(db: Session, accession: str) -> list[dict]:
    """Rollup rows for every SAAP mapped to one protein, by position."""
    agg = _aggregate_subquery()
    stmt = (select(SAAP, *[agg.c[n] for n in _AGG_NAMES])
            .join(agg, agg.c.saap_id == SAAP.id)
            .where(SAAP.protein_accession == accession)
            .order_by(SAAP.position_in_protein.is_(None), SAAP.position_in_protein, SAAP.id))
    return [{**_row_to_dict(row), "peptide_start": row[0].peptide_start}
            for row in db.execute(stmt).all()]


def stats(db: Session):
    return {
        "n_saap": db.scalar(select(func.count(SAAP.id))) or 0,
        "n_observations": db.scalar(select(func.count(Observation.id))) or 0,
        "n_tissues": db.scalar(select(func.count(distinct(Observation.tissue)))) or 0,
        "n_cross_species": db.scalar(select(func.count(SAAP.id)).where(SAAP.cross_species == "same")) or 0,
        "private": PRIVATE,
        "n_genes": db.scalar(select(func.count(distinct(SAAP.source_gene)))) or 0,
        "n_proteins": db.scalar(select(func.count(distinct(SAAP.protein_accession)))) or 0,
    }
