# File ingestion (CSV / TSV / XLSX): parse, map columns, de-duplicate, persist

from __future__ import annotations
import csv
import hashlib
import io
from dataclasses import dataclass, field
from sqlalchemy import func, or_, select
from sqlalchemy.orm import Session
from . import column_map
from .curate import prune
from .proteome import check_all
from .models import Observation, SAAP
from .samples import sample_label
from .util import normalize_dataset, normalize_species

# String value -> bool for flag columns.
_TRUE = {"yes", "true", "1", "y", "t"}
_FALSE = {"no", "false", "0", "n", "f"}


@dataclass
class IngestResult:
    filename: str
    rows_read: int = 0
    saap_created: int = 0
    observations_created: int = 0
    duplicate_observations_skipped: int = 0
    rows_skipped_no_identity: int = 0
    columns_mapped: dict[str, str] = field(default_factory=dict)
    columns_unmapped: list[str] = field(default_factory=list)
    saap_pending_uniprot: int = 0
    rows_skipped_immunoglobulin: int = 0
    curation: dict = field(default_factory=dict)

    def as_dict(self) -> dict:
        return self.__dict__.copy()


def _to_float(value: str | None):
    if value is None:
        return None
    value = value.strip()
    if value == "" or value.lower() in {"na", "nan", "null", "none", "-"}:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def _to_int(value: str | None):
    f = _to_float(value)
    return int(f) if f is not None else None


def _to_bool(value: str | None):
    if value is None:
        return None
    v = value.strip().lower()
    if v in _TRUE:
        return True
    if v in _FALSE:
        return False
    return None


def _clean_str(value: str | None):
    if value is None:
        return None
    v = value.strip()
    return v or None


def _split_peptides(value: str) -> list[str]:
    # Function to split a peptide cell that may hold several alternatives 
    # if row has more than one candidate sequence ('SAVSGLWGK;ASVSGLWGK'),

    if not value:
        return []
    return [p.strip().upper() for p in value.replace(",", ";").split(";") if p.strip()]


def _single_sub(mtp: str, bp: str) -> str:
    # Function to define sub as 'X to Y' for one peptide pair
    if not mtp or not bp or len(mtp) != len(bp):
        return ""
    diffs = [i for i, (a, b) in enumerate(zip(bp, mtp)) if a != b]
    if len(diffs) != 1:
        return ""
    i = diffs[0]
    return f"{bp[i].upper()} to {mtp[i].upper()}"


def _derive_aa_sub(mtp: str, bp: str) -> str:
    """
    Function to derive 'X to Y' from the peptide pair, when the file has no AAS column.

    The base and substituted peptides differ at exactly one residue, so the
    substitution is recoverable. Deriving it here (rather than leaving the
    column blank) keeps the identity key stable: a later import of the same
    peptides *with* an AAS column then matches the existing SAAP instead of
    creating a second row for it.

    A cell may hold several candidate peptides ('SAVSGLWGK;ASVSGLWGK') where the
    search could not place the substitution. Every combination is evaluated and
    the result is used only if they all agree on the same residue change — which
    they typically do, differing only in position.

    Returns "" when the pair doesn't resolve to a single change
    """
    mtps = _split_peptides(mtp)
    bps = _split_peptides(bp)
    if not mtps or not bps:
        return ""
    if len(mtps) == 1 and len(bps) == 1:
        return _single_sub(mtps[0], bps[0])

    subs = {s for m in mtps for b in bps if (s := _single_sub(m, b))}
    # Only trust it when every viable pairing gives the same substitution.
    return subs.pop() if len(subs) == 1 else ""


def _saap_identities(record: dict) -> list[tuple[str, str, str]]:
    """
    Function to Identity keys for one row: (mtp_seq, bp_seq, aa_sub) per SAAP peptide

    A row may list several substituted peptides against one base peptide
    ('SAVSGLWGK;ASVSGLWGK'). They are distinct SAAPs — the same residue change
    at different positions — so each becomes its own row with its own AAS.

    An explicit AAS applies to every peptide in the cell; otherwise each one is
    derived from its own pairing with the base peptide.
    """
    mtp_raw = _clean_str(record.get("mtp_seq"))
    if not mtp_raw:
        return []
    bp = _clean_str(record.get("bp_seq")) or ""
    explicit = _clean_str(record.get("aa_sub")) or ""

    peptides = _split_peptides(mtp_raw)
    if not peptides:
        return []

    out: list[tuple[str, str, str]] = []
    for pep in peptides:
        aa_sub = explicit or _derive_aa_sub(pep, bp)
        out.append((pep, bp, aa_sub))
    # De-duplicate in case the cell repeats a peptide.
    return list(dict.fromkeys(out))


def _cell_to_str(v) -> str:
    """
    Function to normalize an XLSX cell value to the string form our parsers expect.
    Integer-valued floats (e.g. intensities read as 6.14e8) drop the '.0'.
    """
    if v is None:
        return ""
    if isinstance(v, float) and v.is_integer():
        return str(int(v))
    return str(v)


def _parse_delimited(raw_bytes: bytes, name: str) -> tuple[list[str], list[dict]]:
    # Function to parse .csv files
    text = raw_bytes.decode("utf-8-sig", errors="replace")
    delimiter = "\t" if name.endswith(".tsv") else ","
    reader = csv.DictReader(io.StringIO(text), delimiter=delimiter)
    headers = reader.fieldnames or []
    return headers, [dict(row) for row in reader]


def _parse_xlsx(raw_bytes: bytes) -> tuple[list[str], list[dict]]:
    # Function to parse .xlsx files
    from openpyxl import load_workbook

    wb = load_workbook(io.BytesIO(raw_bytes), read_only=True, data_only=True)
    try:
        ws = wb.active
        rows_iter = ws.iter_rows(values_only=True)
        try:
            header_row = next(rows_iter)
        except StopIteration:
            return [], []
        headers = [("" if h is None else str(h).strip()) for h in header_row]
        rows: list[dict] = []
        for r in rows_iter:
            if r is None or all(c is None for c in r):
                continue
            rec = {}
            for i, h in enumerate(headers):
                if not h:
                    continue
                rec[h] = _cell_to_str(r[i] if i < len(r) else None)
            rows.append(rec)
        return headers, rows
    finally:
        wb.close()


def _parse_file(raw_bytes: bytes, filename: str) -> tuple[list[str], list[dict]]:
    name = filename.lower()
    if name.endswith(".xlsx"):
        return _parse_xlsx(raw_bytes)
    return _parse_delimited(raw_bytes, name)


def ingest_file(
    db: Session,
    raw_bytes: bytes,
    filename: str,
) -> IngestResult:
    headers, raw_rows = _parse_file(raw_bytes, filename)

    mapping, unmapped = column_map.map_headers(headers)
    result = IngestResult(
        filename=filename,
        columns_mapped=mapping,
        columns_unmapped=unmapped,
    )

    # Cache SAAP identities already seen in this ingest to avoid extra queries
    # Collapses in-file duplicates 
    saap_cache: dict[tuple[str, str, str], SAAP] = {}
    seen_hashes: set[str] = set()

    for raw_row in raw_rows:
        result.rows_read += 1

        # Translate raw headers -> canonical fields
        record: dict[str, str] = {}
        for raw_header, canonical in mapping.items():
            record[canonical] = raw_row.get(raw_header)

        if _to_bool(record.get("immunoglobulin")):
            result.rows_skipped_immunoglobulin += 1
            continue
        identities = _saap_identities(record)
        if not identities:
            result.rows_skipped_no_identity += 1
            continue

        # A row listing several substituted peptides becomes several SAAPs,
        # each getting its own observation from this row
        for identity in identities:
            saap = saap_cache.get(identity)
            if saap is None:
                saap = db.scalar(
                    select(SAAP).where(
                        SAAP.mtp_seq == identity[0],
                        SAAP.bp_seq == identity[1],
                        SAAP.aa_sub == identity[2],
                    )
                )
                if saap is None:
                    saap = SAAP(
                        mtp_seq=identity[0],
                        bp_seq=identity[1],
                        aa_sub=identity[2],
                        source_accession=_clean_str(record.get("source_accession")),
                        source_gene=_clean_str(record.get("source_gene")),
                        ref_proteins=_clean_str(record.get("ref_proteins")),
                        ensembl_gene=_clean_str(record.get("ensembl_gene")),
                        ensembl_transcript=_clean_str(record.get("ensembl_transcript")),
                        ensembl_protein=_clean_str(record.get("ensembl_protein")),
                        protein_description=_clean_str(record.get("protein_description")),
                        protein_length=_to_int(record.get("protein_length")),
                        position_in_protein=_to_int(record.get("position_in_protein")),
                        peptide_start=_to_int(record.get("peptide_start")),
                        annotation_source=("file" if _clean_str(record.get("ensembl_gene"))
                                           or _clean_str(record.get("position_in_protein"))
                                           else None),
                        immunoglobulin=_to_bool(record.get("immunoglobulin")),
                        trypsin=_to_bool(record.get("trypsin")),
                        missed_cleavage=_to_bool(record.get("missed_cleavage")),
                        aas_at_peptide_terminus=_to_bool(record.get("aas_at_peptide_terminus")),
                        greater_than_shared=_to_bool(record.get("greater_than_shared")),
                    )
                    db.add(saap)
                    db.flush()  # assign PK
                    result.saap_created += 1
                else:
                    # Backfill source metadata if the existing row lacks it
                    _backfill_source(saap, record)
                saap_cache[identity] = saap
            else:
                _backfill_source(saap, record)

            row_hash = _hash_row(identity, record)
            if row_hash in seen_hashes or db.scalar(
                select(Observation.id).where(Observation.row_hash == row_hash)
            ) is not None:
                result.duplicate_observations_skipped += 1
                continue
            seen_hashes.add(row_hash)

            db.add(Observation(
                saap_id=saap.id,
                dataset=normalize_dataset(_clean_str(record.get("dataset"))),
              tissue=sample_label(record.get("dataset"), record.get("tmt_tissue"))[0],
              sample_type=sample_label(record.get("dataset"), record.get("tmt_tissue"))[1],
                tmt_tissue=_clean_str(record.get("tmt_tissue")),
                digest=_clean_str(record.get("digest")),
                species=normalize_species(_clean_str(record.get("species"))),
                acquisition_type=_clean_str(record.get("acquisition_type")),
                saap_pep=_to_float(record.get("saap_pep")),
                positional_probability=_to_float(record.get("positional_probability")),
                n_evidence_fragments=_to_int(record.get("n_evidence_fragments")),
                source_file=filename,
                row_hash=row_hash,
            ))
            result.observations_created += 1

    db.commit()
    result.curation = prune(db)
    check_all(db)

    # SAAPs without a UniProt accession are KEPT. They carry a gene name and/or
    # protein description, which the annotation step resolves to an accession
    # later — deleting them here would throw away real observations.
    result.saap_pending_uniprot = count_saap_without_uniprot(db)
    return result


def count_saap_without_uniprot(db: Session) -> int:
    # Function to count how many SAAP still lack a UniProt accession during annotation
    return db.scalar(
        select(func.count(SAAP.id)).where(or_(SAAP.source_accession.is_(None),
                                              SAAP.source_accession == ""))
    ) or 0


def _backfill_source(saap: SAAP, record: dict) -> None:
    # Function to fill in per-peptide attributes on an existing SAAP if previously blank
    for attr in ("source_accession", "source_gene", "ref_proteins",
                 "ensembl_gene", "ensembl_transcript", "ensembl_protein",
                 "protein_description"):
        if getattr(saap, attr) is None:
            val = _clean_str(record.get(attr))
            if val:
                setattr(saap, attr, val)
    for attr in ("protein_length", "position_in_protein", "peptide_start"):
        if getattr(saap, attr) is None:
            val = _to_int(record.get(attr))
            if val is not None:
                setattr(saap, attr, val)
    for attr in ("immunoglobulin", "trypsin", "missed_cleavage",
                 "aas_at_peptide_terminus", "greater_than_shared"):
        if getattr(saap, attr) is None:
            b = _to_bool(record.get(attr))
            if b is not None:
                setattr(saap, attr, b)


def _hash_row(identity: tuple[str, str, str], record: dict) -> str:
    """Stable hash over identity + all observation fields, so re-uploading the
    exact same line is recognized as a duplicate observation."""
    parts = list(identity)
    for f in sorted(column_map.OBSERVATION_FIELDS):
        parts.append(f"{f}={(record.get(f) or '').strip()}")
    return hashlib.sha1("\x1f".join(parts).encode("utf-8")).hexdigest()
