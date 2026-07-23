"""Ensembl / positional annotation for SAAP records.

Two things are resolved here, both keyed off the UniProt accession already
stored on each SAAP:

  1. Cross-references — Ensembl gene (ENSG), transcript (ENST) and protein
     (ENSP) IDs, plus the protein description and length.
  2. Position — the 1-based index of the substituted residue within the full
     protein, found by locating `bp_seq` (the base peptide) in the canonical
     sequence and adding the offset of the substituted residue within it.

Annotation is an *optional enrichment*: it needs network access to UniProt, so
every entry point degrades gracefully. If a lookup fails the SAAP keeps whatever
it already had and is simply reported as unresolved — export and CSV never fail
because annotation could not run.

Values already present from the imported file are never overwritten by a lookup
(file columns win; see `annotate_saaps(overwrite=False)`).
"""
from __future__ import annotations

import re
import time
from dataclasses import dataclass, field
from typing import Iterable, Optional

from sqlalchemy import or_, select
from sqlalchemy.orm import Session

from .models import SAAP

UNIPROT_BASE = "https://rest.uniprot.org/uniprotkb"
# Only the fields we actually consume, to keep responses small.
UNIPROT_FIELDS = "accession,id,protein_name,gene_names,length,sequence,xref_ensembl"

DEFAULT_BATCH_SIZE = 100   # accessions per UniProt search query
DEFAULT_TIMEOUT = 30       # seconds per HTTP request
DEFAULT_PAUSE = 0.2        # polite delay between batches

# "V to P" / "V->P" / "V/P" / "V2P" -> ("V", "P")
_SUB_SPLIT = re.compile(r"\s*(?:to|->|>|/|→|2)\s*", re.IGNORECASE)
_ENS_RE = re.compile(r"^ENS[A-Z]*[GTP]\d+", re.IGNORECASE)


# ----------------------------- data containers -----------------------------
@dataclass
class ProteinRecord:
    """The subset of a UniProt entry we care about."""
    accession: str
    description: Optional[str] = None
    gene: Optional[str] = None
    sequence: Optional[str] = None
    length: Optional[int] = None
    ensembl_gene: Optional[str] = None
    ensembl_transcript: Optional[str] = None
    ensembl_protein: Optional[str] = None


@dataclass
class AnnotationResult:
    requested: int = 0
    resolved: int = 0            # got a UniProt record back
    positioned: int = 0          # base peptide located -> position computed
    unmatched_peptide: int = 0   # record found, but bp_seq not in the sequence
    not_found: int = 0           # accession returned nothing
    failed: int = 0              # network/parse error
    errors: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["errors"] = self.errors[:10]  # cap: this goes to the UI
        return d


# ------------------------------ substitutions ------------------------------
def parse_substitution(aa_sub: str | None) -> tuple[Optional[str], Optional[str]]:
    """'V to P' -> ('V', 'P'). Returns (None, None) if unparseable."""
    if not aa_sub:
        return (None, None)
    parts = [p.strip().upper() for p in _SUB_SPLIT.split(aa_sub.strip()) if p.strip()]
    if len(parts) == 2 and all(len(p) == 1 and p.isalpha() for p in parts):
        return (parts[0], parts[1])
    return (None, None)


def substitution_offset(bp_seq: str | None, mtp_seq: str | None) -> Optional[int]:
    """0-based offset of the substituted residue within the peptide.

    Prefers a direct base-vs-variant comparison (exactly one differing residue).
    Equal-length sequences differing at one position give an unambiguous answer;
    anything else returns None and the caller falls back to the AAS letters.
    """
    if not bp_seq or not mtp_seq or len(bp_seq) != len(mtp_seq):
        return None
    diffs = [i for i, (a, b) in enumerate(zip(bp_seq, mtp_seq)) if a != b]
    return diffs[0] if len(diffs) == 1 else None


def _offset_from_aas(bp_seq: str, aa_sub: str | None) -> Optional[int]:
    """Fallback: locate the substituted residue using the AAS 'from' letter.
    Only trusted when that residue occurs exactly once in the peptide."""
    frm, _ = parse_substitution(aa_sub)
    if not frm:
        return None
    hits = [i for i, ch in enumerate(bp_seq.upper()) if ch == frm]
    return hits[0] if len(hits) == 1 else None


def locate_peptide(protein_seq: str | None, bp_seq: str | None) -> Optional[int]:
    """1-based start of `bp_seq` within the protein; None if absent or ambiguous.

    Leucine/isoleucine are indistinguishable by mass, so a peptide may be written
    with either. If the literal sequence is not found we retry with I and L
    treated as equivalent.
    """
    if not protein_seq or not bp_seq:
        return None
    prot, pep = protein_seq.upper(), bp_seq.upper()
    idx = prot.find(pep)
    if idx >= 0:
        # Ambiguous if the peptide occurs more than once.
        if prot.find(pep, idx + 1) >= 0:
            return None
        return idx + 1
    prot_x, pep_x = prot.replace("I", "L"), pep.replace("I", "L")
    idx = prot_x.find(pep_x)
    if idx >= 0 and prot_x.find(pep_x, idx + 1) < 0:
        return idx + 1
    return None


def compute_position(protein_seq: str | None, saap: SAAP) -> tuple[Optional[int], Optional[int]]:
    """Return (peptide_start, position_in_protein), both 1-based, or (None, None)."""
    start = locate_peptide(protein_seq, saap.bp_seq)
    if start is None:
        return (None, None)
    offset = substitution_offset(saap.bp_seq, saap.mtp_seq)
    if offset is None:
        offset = _offset_from_aas(saap.bp_seq or "", saap.aa_sub)
    if offset is None:
        return (start, None)
    return (start, start + offset)


# ------------------------------ UniProt client ------------------------------
def _first_ensembl(xrefs: list[dict]) -> tuple[Optional[str], Optional[str], Optional[str]]:
    """Pull (gene, transcript, protein) Ensembl IDs out of UniProt xrefs.

    In UniProt's JSON an Ensembl cross-reference carries the transcript as `id`,
    with the protein and gene as nested properties.
    """
    gene = transcript = protein = None
    for xref in xrefs or []:
        if (xref.get("database") or "").lower() != "ensembl":
            continue
        xid = xref.get("id")
        if xid and transcript is None and _ENS_RE.match(xid):
            transcript = xid.split(".")[0]
        for prop in xref.get("properties") or []:
            key = (prop.get("key") or "").lower()
            val = (prop.get("value") or "").split(".")[0]
            if not val or not _ENS_RE.match(val):
                continue
            if "gene" in key and gene is None:
                gene = val
            elif "protein" in key and protein is None:
                protein = val
        if gene and transcript and protein:
            break
    return (gene, transcript, protein)


def _parse_entry(entry: dict) -> ProteinRecord:
    acc = entry.get("primaryAccession") or ""
    desc = None
    pd = entry.get("proteinDescription") or {}
    rec_name = (pd.get("recommendedName") or {}).get("fullName") or {}
    if rec_name.get("value"):
        desc = rec_name["value"]
    else:
        subs = pd.get("submissionNames") or []
        if subs:
            desc = ((subs[0].get("fullName") or {}).get("value")) or None

    genes = entry.get("genes") or []
    gene = None
    if genes:
        gene = ((genes[0].get("geneName") or {}).get("value")) or None

    seq_block = entry.get("sequence") or {}
    sequence = seq_block.get("value")
    length = seq_block.get("length")

    ens_g, ens_t, ens_p = _first_ensembl(entry.get("uniProtKBCrossReferences") or [])
    return ProteinRecord(
        accession=acc, description=desc, gene=gene,
        sequence=sequence, length=length,
        ensembl_gene=ens_g, ensembl_transcript=ens_t, ensembl_protein=ens_p,
    )


def fetch_uniprot(
    accessions: Iterable[str],
    *,
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    pause: float = DEFAULT_PAUSE,
    session=None,
) -> tuple[dict[str, ProteinRecord], list[str]]:
    """Fetch UniProt entries for `accessions`.

    Returns ({accession: ProteinRecord}, [error strings]). Network failures are
    collected rather than raised so a partial result is still usable. `requests`
    is imported lazily so the app runs without it when annotation is unused.
    """
    accs = [a for a in dict.fromkeys(a.strip() for a in accessions if a and a.strip())]
    if not accs:
        return ({}, [])

    try:
        import requests  # noqa: PLC0415  (optional dependency)
    except ImportError:
        return ({}, ["The 'requests' package is required for annotation "
                     "(pip install -r backend/requirements.txt)."])

    sess = session or requests.Session()
    out: dict[str, ProteinRecord] = {}
    errors: list[str] = []

    for i in range(0, len(accs), batch_size):
        chunk = accs[i:i + batch_size]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        try:
            resp = sess.get(
                f"{UNIPROT_BASE}/search",
                params={"query": query, "fields": UNIPROT_FIELDS,
                        "format": "json", "size": str(len(chunk))},
                timeout=timeout,
                headers={"Accept": "application/json"},
            )
            resp.raise_for_status()
            for entry in resp.json().get("results", []):
                rec = _parse_entry(entry)
                if rec.accession:
                    out[rec.accession] = rec
                # Index secondary accessions too, so a query by an old ID resolves.
                for sec in entry.get("secondaryAccessions") or []:
                    out.setdefault(sec, rec)
        except Exception as exc:  # network, HTTP, or JSON problem
            errors.append(f"{type(exc).__name__} for {chunk[0]}…{chunk[-1]}: {exc}")
        if pause and i + batch_size < len(accs):
            time.sleep(pause)

    return (out, errors)


# ------------------------------ orchestration ------------------------------
def _base_accession(value: str | None) -> str:
    """First accession of a possibly ';'-separated list, minus any isoform suffix."""
    if not value:
        return ""
    for part in value.split(";"):
        part = part.strip()
        if part:
            return part.split("-")[0]
    return ""


def apply_record(saap: SAAP, rec: ProteinRecord, *, overwrite: bool = False) -> bool:
    """Copy a ProteinRecord onto a SAAP and compute its position.

    With overwrite=False (the default) only blank fields are filled, so values
    that came from the imported file are preserved. Returns True if anything
    changed.
    """
    changed = False

    def _set(attr: str, value):
        nonlocal changed
        if value in (None, ""):
            return
        if overwrite or getattr(saap, attr) in (None, ""):
            if getattr(saap, attr) != value:
                setattr(saap, attr, value)
                changed = True

    _set("ensembl_gene", rec.ensembl_gene)
    _set("ensembl_transcript", rec.ensembl_transcript)
    _set("ensembl_protein", rec.ensembl_protein)
    _set("protein_description", rec.description)
    _set("protein_length", rec.length)
    # Cached so full-protein FASTA export works offline, without re-querying.
    _set("protein_sequence", rec.sequence)
    _set("ref_proteins", rec.description)
    _set("source_gene", rec.gene)

    start, pos = compute_position(rec.sequence, saap)
    _set("peptide_start", start)
    _set("position_in_protein", pos)
    return changed


def apply_substitution(saap: SAAP) -> tuple[Optional[str], Optional[str]]:
    """Build the full-length protein sequence carrying this SAAP's substitution.

    Returns (variant_sequence, error). Exactly one of the two is set.

    The residue at `position_in_protein` is replaced with the substituted amino
    acid. Before writing, the residue currently at that position is checked
    against what the substitution says should be there — a mismatch means the
    position and the sequence disagree (wrong isoform, stale annotation), so the
    entry is skipped rather than silently emitting a wrong protein.

    The substituted residue is taken from the SAAP/BP comparison where possible
    and from the AAS column otherwise, mirroring how the position was derived.
    """
    seq = saap.protein_sequence
    pos = saap.position_in_protein
    if not seq:
        return (None, "no cached protein sequence — run annotation first")
    if pos is None:
        return (None, "no position in protein")
    if pos < 1 or pos > len(seq):
        return (None, f"position {pos} outside protein (length {len(seq)})")

    frm, to = parse_substitution(saap.aa_sub)
    # Prefer the actual peptide comparison; it reflects the observed data.
    offset = substitution_offset(saap.bp_seq, saap.mtp_seq)
    if offset is not None:
        frm = (saap.bp_seq or "")[offset].upper()
        to = (saap.mtp_seq or "")[offset].upper()
    if not to:
        return (None, f"cannot determine substituted residue from AAS {saap.aa_sub!r}")

    actual = seq[pos - 1].upper()
    if frm and actual != frm:
        # I/L are isobaric, so treat them as interchangeable before rejecting.
        if not (actual in "IL" and frm in "IL"):
            return (None, f"residue at {pos} is {actual!r}, substitution expects {frm!r}")

    return (seq[:pos - 1] + to + seq[pos:], None)


def annotate_saaps(
    db: Session,
    *,
    ids: list[int] | None = None,
    only_missing: bool = True,
    overwrite: bool = False,
    limit: int | None = None,
    batch_size: int = DEFAULT_BATCH_SIZE,
    timeout: int = DEFAULT_TIMEOUT,
    session=None,
) -> AnnotationResult:
    """Annotate SAAP rows from UniProt.

    ids          — restrict to these SAAP ids (default: all).
    only_missing — skip rows that already have an Ensembl gene and a position.
    overwrite    — replace existing values instead of filling blanks only.
    limit        — cap the number of rows processed (useful for a trial run).
    """
    stmt = select(SAAP).where(SAAP.source_accession.is_not(None),
                              SAAP.source_accession != "")
    if ids:
        stmt = stmt.where(SAAP.id.in_(ids))
    if only_missing and not overwrite:
        stmt = stmt.where(or_(SAAP.ensembl_gene.is_(None),
                              SAAP.position_in_protein.is_(None),
                              SAAP.protein_sequence.is_(None)))
    stmt = stmt.order_by(SAAP.id)
    if limit:
        stmt = stmt.limit(limit)

    saaps = list(db.scalars(stmt).all())
    result = AnnotationResult(requested=len(saaps))
    if not saaps:
        return result

    by_acc: dict[str, list[SAAP]] = {}
    for s in saaps:
        acc = _base_accession(s.source_accession)
        if acc:
            by_acc.setdefault(acc, []).append(s)

    records, errors = fetch_uniprot(
        by_acc.keys(), batch_size=batch_size, timeout=timeout, session=session
    )
    result.errors.extend(errors)

    for acc, group in by_acc.items():
        rec = records.get(acc)
        if rec is None:
            if errors:
                result.failed += len(group)
            else:
                result.not_found += len(group)
                for s in group:
                    s.annotation_source = s.annotation_source or "uniprot:not-found"
            continue
        for s in group:
            result.resolved += 1
            apply_record(s, rec, overwrite=overwrite)
            if s.position_in_protein is not None:
                result.positioned += 1
                s.annotation_source = "uniprot"
            else:
                result.unmatched_peptide += 1
                s.annotation_source = "uniprot:unmatched"

    db.commit()
    return result
