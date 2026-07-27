"""
Ensembl / positional annotation for SAAP records.

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
from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Optional
from sqlalchemy import func, or_, select
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


# data containers 
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
    resolved_by_gene: int = 0    # had no accession; found via gene symbol
    no_identifier: int = 0       # neither accession nor gene -> unresolvable
    resolved_by_peptide: int = 0 # identifiers copied from a known base peptide
    resolved_by_sequence: int = 0 # protein found by searching the peptide sequence
    species_corrected: int = 0   # annotation re-resolved into the right species
    aas_filled: int = 0          # AAS string derived from the peptide pair
    merged_duplicates: int = 0   # empty shells folded into an existing SAAP
    unmatched_examples: list[str] = field(default_factory=list)
    errors: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = self.__dict__.copy()
        d["errors"] = self.errors[:10]  # cap: this goes to the UI
        d["unmatched_examples"] = self.unmatched_examples[:10]
        return d


# substitutions 
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


def locate_peptide_all(protein_seq: str | None, bp_seq: str | None) -> list[int]:
    """Every 1-based start of `bp_seq` within the protein, in order.

    Repeat-rich proteins (collagens, low-complexity regions) legitimately
    contain the same tryptic peptide several times. Each occurrence is a valid
    candidate site, so all are returned rather than discarding the peptide.

    Leucine/isoleucine are indistinguishable by mass, so a peptide may be
    recorded with either. The literal sequence is tried first; only if that
    finds nothing are I and L treated as equivalent.
    """
    if not protein_seq or not bp_seq:
        return []
    prot, pep = protein_seq.upper(), bp_seq.upper()

    def _scan(haystack: str, needle: str) -> list[int]:
        out: list[int] = []
        i = haystack.find(needle)
        while i >= 0:
            out.append(i + 1)
            i = haystack.find(needle, i + 1)
        return out

    hits = _scan(prot, pep)
    if hits:
        return hits
    return _scan(prot.replace("I", "L"), pep.replace("I", "L"))


def locate_peptide(protein_seq: str | None, bp_seq: str | None) -> Optional[int]:
    """1-based start of `bp_seq`, or None when absent.

    Returns the first occurrence; use `locate_peptide_all` when every site
    matters. Kept for callers that need a single number.
    """
    hits = locate_peptide_all(protein_seq, bp_seq)
    return hits[0] if hits else None


def contains_peptide(protein_seq: str | None, bp_seq: str | None) -> bool:
    """Whether the protein contains this peptide at all.

    Distinct from `locate_peptide`, which yields a position only when the
    peptide occurs exactly once. A peptide repeated within a protein (common in
    collagens and other repeat-rich sequences) still identifies that protein
    correctly — we just cannot say which copy carries the substitution. Used to
    accept a protein match while leaving the position blank.
    """
    if not protein_seq or not bp_seq:
        return False
    prot, pep = protein_seq.upper(), bp_seq.upper()
    if pep in prot:
        return True
    # Leucine/isoleucine are isobaric and may be recorded interchangeably.
    return pep.replace("I", "L") in prot.replace("I", "L")


def compute_positions(protein_seq: str | None, saap: SAAP) -> tuple[list[int], list[int]]:
    """Return (peptide_starts, substitution_positions), both 1-based lists.

    A peptide occurring several times in the protein yields one candidate site
    per occurrence; all are returned. Empty lists mean the peptide was not found
    (or the substituted residue could not be placed within it).
    """
    starts = locate_peptide_all(protein_seq, saap.bp_seq)
    if not starts:
        return ([], [])
    offset = substitution_offset(saap.bp_seq, saap.mtp_seq)
    if offset is None:
        offset = _offset_from_aas(saap.bp_seq or "", saap.aa_sub)
    if offset is None:
        return (starts, [])
    return (starts, [s + offset for s in starts])


def compute_position(protein_seq: str | None, saap: SAAP) -> tuple[Optional[int], Optional[int]]:
    """First (peptide_start, position_in_protein), or (None, None)."""
    starts, positions = compute_positions(protein_seq, saap)
    return (starts[0] if starts else None, positions[0] if positions else None)


# - UniProt client -
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


# - orchestration -
def _accession_candidates(value: str | None) -> list[str]:
    """All accessions to try for one SAAP, in priority order.

    The source field may hold several accessions ('E5RHP7;P00915') and may carry
    an isoform suffix ('O00429-2'). Both forms are kept: the isoform is tried
    first because the peptide was matched against that specific sequence, with
    the canonical entry as a fallback.
    """
    if not value:
        return []
    out: list[str] = []
    for part in str(value).replace(",", ";").split(";"):
        part = part.strip()
        if not part:
            continue
        out.append(part)             # e.g. 'O00429-2' (isoform-specific)
        base = part.split("-")[0]
        if base != part:
            out.append(base)         # e.g. 'O00429' (canonical fallback)
    # De-duplicate, preserving order.
    return list(dict.fromkeys(out))


def _base_accession(value: str | None) -> str:
    """First accession of a possibly ';'-separated list, minus any isoform suffix."""
    if not value:
        return ""
    for part in value.split(";"):
        part = part.strip()
        if part:
            return part.split("-")[0]
    return ""


def gene_species(gene: str | None) -> str | None:
    """'homo sapiens' / 'mus musculus' implied by a gene symbol's casing.

    Human symbols are all-caps (SPTAN1), mouse are title-case (Sptan1). Returns
    None when the symbol fits neither pattern.
    """
    g = (gene or "").split(";")[0].strip()
    if not g:
        return None
    if re.fullmatch(r"[A-Z0-9][A-Z0-9-]*", g) and not re.fullmatch(r"[A-Z][a-z].*", g):
        return "homo sapiens"
    if re.fullmatch(r"[A-Z][a-z0-9][A-Za-z0-9-]*", g):
        return "mus musculus"
    return None


def annotation_matches_species(saap: SAAP, species: str | None) -> bool:
    """Whether the stored gene symbol belongs to `species`.

    Used to decide if an existing annotation should be replaced: a mouse SAAP
    carrying a human symbol has been annotated against the wrong proteome and
    must be re-resolved rather than preserved.
    """
    if not species:
        return True
    implied = gene_species(saap.source_gene)
    return implied is None or implied == species


def _first_gene(value: str | None) -> str:
    """First usable gene symbol from a ';'-separated list.

    Source files sometimes give several synonyms (';GIG42;ALB;PRO2044'), often
    with empty leading fields. The first non-empty token is used.
    """
    if not value:
        return ""
    for part in value.replace(",", ";").split(";"):
        part = part.strip()
        # Skip obvious internal placeholders that aren't real symbols.
        if part and not part.startswith(("UNQ", "PRO", "GIG")):
            return part
    for part in value.replace(",", ";").split(";"):
        if part.strip():
            return part.strip()
    return ""


def _saap_species_name(db: Session, saap: SAAP) -> str | None:
    """Species the annotation must belong to.

    The rule is fixed, not evidence-weighted: a SAAP observed only in mouse is
    annotated against mouse, only in human against human, and one seen in both
    defaults to human. A single SAAP row carries one gene and accession, so a
    peptide conserved across both species has to pick one proteome — human is
    the reference of record here, which also keeps orthologous peptides from
    landing on mouse or human arbitrarily depending on observation counts.
    """
    from .models import Observation

    names = {
        (r[0] or "").strip().lower()
        for r in db.execute(
            select(Observation.species)
            .where(Observation.saap_id == saap.id, Observation.species.is_not(None))
            .distinct()
        ).all()
    }
    names.discard("")
    if not names:
        return None
    human = {n for n in names if n.startswith("homo")}
    if human:
        return sorted(human)[0]      # human wins outright, including "both"
    return sorted(names)[0]


def _saap_organism(db: Session, saap: SAAP) -> str | None:
    """NCBI taxon id for a SAAP, from its own observations.

    The database mixes species (human and mouse here), so a single global
    organism filter would send mouse genes to human entries. Each SAAP is
    resolved against the species carrying most of its observations; when two
    species tie, no filter is applied rather than picking one arbitrarily.
    """
    name = _saap_species_name(db, saap)
    return TAXA.get(name) if name else None


TAXA = {
    "human": "9606", "homo sapiens": "9606",
    "mouse": "10090", "mus musculus": "10090",
    "rat": "10116", "rattus norvegicus": "10116",
}


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

    starts, positions = compute_positions(rec.sequence, saap)
    _set("peptide_start", starts[0] if starts else None)
    _set("position_in_protein", positions[0] if positions else None)
    # Every candidate site, for peptides that repeat within the protein.
    _set("peptide_starts_all", ",".join(str(x) for x in starts) if starts else None)
    _set("positions_all", ",".join(str(x) for x in positions) if positions else None)
    _set("n_positions", len(positions) if positions else (len(starts) or None))
    return changed


def fetch_uniprot_by_gene(
    genes: Iterable[str],
    *,
    organism_id: str | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    pause: float = DEFAULT_PAUSE,
    session=None,
) -> tuple[dict[str, ProteinRecord], list[str]]:
    """Resolve gene symbols to UniProt records, one query per gene.

    Used for SAAPs that arrived without an accession. Reviewed (Swiss-Prot)
    entries are preferred, and the search is restricted to a single organism
    when given, since the same symbol exists across species. Genes are queried
    individually because a combined OR query gives no way to tell which entry
    answered which symbol.

    Returns ({gene_upper: ProteinRecord}, [errors]).
    """
    syms = [g for g in dict.fromkeys(g.strip() for g in genes if g and g.strip())]
    if not syms:
        return ({}, [])

    try:
        import requests  # noqa: PLC0415
    except ImportError:
        return ({}, ["The 'requests' package is required for annotation "
                     "(pip install -r backend/requirements.txt)."])

    sess = session or requests.Session()
    out: dict[str, ProteinRecord] = {}
    errors: list[str] = []

    for i, sym in enumerate(syms):
        # Swiss-Prot first; many mouse proteins only exist as TrEMBL entries, so
        # retry without the reviewed filter rather than losing them.
        attempts = ["AND reviewed:true", ""]
        for suffix in attempts:
            query = f"gene_exact:{sym} {suffix}".strip()
            if organism_id:
                query += f" AND organism_id:{organism_id}"
            try:
                resp = sess.get(
                    f"{UNIPROT_BASE}/search",
                    params={"query": query, "fields": UNIPROT_FIELDS,
                            "format": "json", "size": "1"},
                    timeout=timeout,
                    headers={"Accept": "application/json"},
                )
                resp.raise_for_status()
                results = resp.json().get("results", [])
                if results:
                    out[sym.upper()] = _parse_entry(results[0])
                    break
            except Exception as exc:
                errors.append(f"{type(exc).__name__} for gene {sym}: {exc}")
                break
        if pause and i + 1 < len(syms):
            time.sleep(pause)

    return (out, errors)


def fetch_isoforms(
    accessions: Iterable[str],
    *,
    timeout: int = DEFAULT_TIMEOUT,
    pause: float = DEFAULT_PAUSE,
    session=None,
) -> tuple[dict[str, ProteinRecord], list[str]]:
    """Fetch specific isoform entries (e.g. 'O00429-2').

    UniProt's search endpoint only returns canonical sequences, so isoforms are
    retrieved one at a time from the entry endpoint. Peptides observed on a
    non-canonical isoform will not be found in the canonical sequence, so this
    is what makes those SAAPs positionable.
    """
    accs = [a for a in dict.fromkeys(a.strip() for a in accessions if a and "-" in a)]
    if not accs:
        return ({}, [])
    try:
        import requests  # noqa: PLC0415
    except ImportError:
        return ({}, ["The 'requests' package is required for annotation."])

    sess = session or requests.Session()
    out: dict[str, ProteinRecord] = {}
    errors: list[str] = []
    for i, acc in enumerate(accs):
        try:
            resp = sess.get(
                f"{UNIPROT_BASE}/{acc}",
                params={"fields": UNIPROT_FIELDS, "format": "json"},
                timeout=timeout,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 404:
                continue
            resp.raise_for_status()
            rec = _parse_entry(resp.json())
            # Key by the requested isoform id, not the primary accession.
            rec.accession = acc
            out[acc] = rec
        except Exception as exc:
            errors.append(f"{type(exc).__name__} for isoform {acc}: {exc}")
        if pause and i + 1 < len(accs):
            time.sleep(pause)
    return (out, errors)


PEPTIDE_SEARCH_BASE = "https://peptidesearch.uniprot.org/asyncrest/"


def fetch_by_peptide(
    peptides: Iterable[str],
    *,
    organism_id: str | None = None,
    timeout: int = DEFAULT_TIMEOUT,
    pause: float = DEFAULT_PAUSE,
    session=None,
) -> tuple[dict[str, ProteinRecord], list[str]]:
    """Find the protein containing each base peptide, via UniProt Peptide Search.

    Last resort for SAAPs imported with no accession and no gene symbol: the
    peptide sequence is the only identifier left.

    This uses the dedicated Peptide Search service
    (https://peptidesearch.uniprot.org/asyncrest/), NOT the UniProtKB `search`
    endpoint — UniProtKB has no query field for "sequence contains this
    stretch", so searching there returns nothing at all. The service is
    asynchronous: POST the peptides, poll the returned job URL, then read the
    matched accessions and fetch their full records.

    Leucine/isoleucine are treated as equivalent (`lEQi`), matching how the
    rest of this module handles isobaric residues. Swiss-Prot-only is tried
    first, then the search is repeated across all of UniProt, since many mouse
    proteins exist only as TrEMBL entries.

    A peptide is accepted only when its matches agree on one gene; anything
    shared across genuinely different genes is left unassigned.

    Returns ({peptide_upper: ProteinRecord}, [errors]).
    """
    peps = [p for p in dict.fromkeys(p.strip().upper() for p in peptides if p and p.strip())]
    if not peps:
        return ({}, [])

    try:
        import requests  # noqa: PLC0415
    except ImportError:
        return ({}, ["The 'requests' package is required for annotation."])

    sess = session or requests.Session()
    out: dict[str, ProteinRecord] = {}
    errors: list[str] = []
    ambiguous: dict[str, list[str]] = {}

    # Peptides are searched one at a time: the service returns a flat list of
    # accessions with no indication of which peptide produced which hit.
    for i, pep in enumerate(peps):
        if len(pep) < 3:
            continue  # the service requires three or more residues
        accessions: list[str] = []
        for sp_only in ("on", ""):
            try:
                accessions = _peptide_search_job(
                    sess, pep, organism_id, sp_only, timeout=timeout
                )
            except Exception as exc:
                errors.append(f"{type(exc).__name__} for peptide {pep}: {exc}")
                accessions = []
                break
            if accessions:
                break

        if not accessions:
            continue

        # Resolve the matched accessions to full records, dropping isoform
        # suffixes (the service returns e.g. 'Q13740-2' alongside 'Q13740').
        bases = list(dict.fromkeys(a.split("-")[0] for a in accessions))[:25]
        records, errs = fetch_uniprot(bases, timeout=timeout, session=sess, pause=0)
        errors.extend(errs)
        hits = [r for r in records.values()
                if r.sequence and contains_peptide(r.sequence, pep)]
        if not hits:
            continue
        genes = {(h.gene or "").upper() for h in hits if h.gene}
        if len(hits) == 1 or len(genes) == 1:
            out[pep] = hits[0]
        else:
            ambiguous[pep] = sorted(genes)[:5]

        if pause and i + 1 < len(peps):
            time.sleep(pause)

    # Report peptides that matched several distinct genes, so a blank cell is
    # explainable rather than mysterious.
    for pep, genes in list(ambiguous.items())[:10]:
        if pep not in out:
            errors.append(f"{pep}: maps to several genes ({', '.join(genes)}) — left unassigned")

    return (out, errors)


def _peptide_search_job(
    sess,
    peptide: str,
    organism_id: str | None,
    sp_only: str,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    max_polls: int = 20,
    poll_wait: float = 3.0,
) -> list[str]:
    """Run one Peptide Search job and return the matched accessions.

    POST -> 202 with a job Location -> poll until 200 (the service answers 303
    with Retry-After while the job is still running). Returns [] if the job
    yields nothing within `max_polls`.
    """
    data = {"peps": peptide, "lEQi": "on"}
    if organism_id:
        data["taxIds"] = organism_id
    if sp_only:
        data["spOnly"] = sp_only

    resp = sess.post(PEPTIDE_SEARCH_BASE, data=data, timeout=timeout,
                     allow_redirects=False)
    if resp.status_code not in (200, 202):
        resp.raise_for_status()
    job_url = resp.headers.get("Location")
    if not job_url:
        # Some deployments answer synchronously with the accession list.
        body = (resp.text or "").strip()
        return [a for a in body.split(",") if a] if body else []

    for _ in range(max_polls):
        job = sess.get(job_url, timeout=timeout, allow_redirects=False)
        if job.status_code == 200:
            body = (job.text or "").strip()
            return [a for a in body.split(",") if a] if body else []
        if job.status_code in (301, 302, 303):
            # Still running; honour Retry-After when the server supplies it.
            wait = job.headers.get("Retry-After")
            time.sleep(min(float(wait), 10.0) if wait and str(wait).isdigit() else poll_wait)
            continue
        job.raise_for_status()
        break
    return []


def stored_positions(saap: SAAP) -> list[int]:
    """Every candidate substitution site recorded for this SAAP.

    Reads the comma-separated `positions_all` when present, falling back to the
    single `position_in_protein` for rows annotated before that column existed.
    """
    raw = (saap.positions_all or "").strip()
    if raw:
        out = []
        for part in raw.split(","):
            part = part.strip()
            if part.isdigit():
                out.append(int(part))
        if out:
            return out
    return [saap.position_in_protein] if saap.position_in_protein else []


def apply_substitutions_all(saap: SAAP) -> tuple[list[tuple[int, str]], Optional[str]]:
    """One variant protein per candidate site.

    Returns ([(position, sequence), ...], error). A peptide repeating within its
    protein yields several entries — each a full-length sequence with the
    substitution applied at one of the sites — since any of them could be the
    real one.
    """
    seq = saap.protein_sequence
    if not seq:
        return ([], "no cached protein sequence — run annotation first")
    positions = stored_positions(saap)
    if not positions:
        from .ingest import _split_peptides

        if len(_split_peptides(saap.mtp_seq or "")) > 1:
            return ([], "ambiguous: alternative peptides give different positions")
        return ([], "no position in protein")

    frm, to = parse_substitution(saap.aa_sub)
    offset = substitution_offset(saap.bp_seq, saap.mtp_seq)
    if offset is not None:
        frm = (saap.bp_seq or "")[offset].upper()
        to = (saap.mtp_seq or "")[offset].upper()
    if not to:
        return ([], f"cannot determine substituted residue from AAS {saap.aa_sub!r}")

    out: list[tuple[int, str]] = []
    skipped: list[str] = []
    for pos in positions:
        if pos < 1 or pos > len(seq):
            skipped.append(f"position {pos} outside protein (length {len(seq)})")
            continue
        actual = seq[pos - 1].upper()
        if frm and actual != frm and not (actual in "IL" and frm in "IL"):
            skipped.append(f"residue at {pos} is {actual!r}, substitution expects {frm!r}")
            continue
        out.append((pos, seq[:pos - 1] + to + seq[pos:]))

    if not out:
        return ([], skipped[0] if skipped else "no usable position")
    return (out, None)


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
        # Several candidate peptides place the substitution at different
        # positions, so there is no single residue to change.
        from .ingest import _split_peptides

        if len(_split_peptides(saap.mtp_seq or "")) > 1:
            return (None, "ambiguous: alternative peptides give different positions")
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


def _merge_duplicate(db: Session, saap: SAAP, aa_sub: str) -> bool:
    """Fold `saap` into an existing SAAP that already has this exact AAS.

    Filling in a blank AAS can make a row identical to one already stored
    (same SAAP/BP/AAS), which the unique constraint forbids. That duplicate is
    the *same* peptide imported without its annotation columns, so its
    observations are moved onto the existing row and the empty shell deleted.

    Returns True if a merge happened.
    """
    from .models import Observation

    twin = db.scalar(
        select(SAAP).where(
            SAAP.mtp_seq == saap.mtp_seq,
            SAAP.bp_seq == saap.bp_seq,
            SAAP.aa_sub == aa_sub,
            SAAP.id != saap.id,
        )
    )
    if twin is None:
        return False

    # Re-point observations, skipping any that would duplicate an existing row.
    existing_hashes = set(
        db.scalars(select(Observation.row_hash).where(Observation.saap_id == twin.id)).all()
    )
    for obs in db.scalars(select(Observation).where(Observation.saap_id == saap.id)).all():
        if obs.row_hash in existing_hashes:
            db.delete(obs)
        else:
            obs.saap_id = twin.id
            existing_hashes.add(obs.row_hash)
    db.delete(saap)
    return True


def resolve_from_database(db: Session, saaps: list[SAAP]) -> int:
    """Fill in identifiers by matching base peptides already known in the DB.

    An import that omits the gene/UniProt columns still carries the base
    peptide, and the same peptide is often already annotated from an earlier
    import. Copying those identifiers across costs no network calls and lets the
    normal UniProt path take over from there.

    Only exact base-peptide matches are used, and only when every annotated SAAP
    sharing that peptide agrees on the accession — a peptide found in several
    different proteins is ambiguous, so it is left alone.

    Returns the number of SAAP updated.
    """
    targets = [s for s in saaps
               if not (s.source_accession or "").strip()
               and not (s.source_gene or "").strip()
               and (s.bp_seq or "").strip()]
    if not targets:
        return 0

    wanted = {s.bp_seq for s in targets}
    donors: dict[str, list[SAAP]] = {}
    # Chunked IN() lookup: SQLite caps variables per statement.
    peps = list(wanted)
    for i in range(0, len(peps), 400):
        chunk = peps[i:i + 400]
        rows = db.scalars(
            select(SAAP).where(
                SAAP.bp_seq.in_(chunk),
                SAAP.source_accession.is_not(None),
                SAAP.source_accession != "",
            )
        ).all()
        for r in rows:
            donors.setdefault(r.bp_seq, []).append(r)

    updated = 0
    for s in targets:
        group = donors.get(s.bp_seq)
        if not group:
            continue

        # A peptide conserved across orthologues maps to several accessions
        # (e.g. mouse Hbb-bs and human HBB). Prefer donors observed in the same
        # species as this SAAP before treating the match as ambiguous.
        if len({(d.source_accession or "").strip() for d in group}) > 1:
            species = _saap_species_name(db, s)
            same = [d for d in group if _saap_species_name(db, d) == species] if species else []
            if not same:
                continue  # no species to disambiguate with
            accs = [(d.source_accession or "").strip() for d in same]
            if len(set(accs)) == 1:
                group = same
            else:
                # Real peptides are often shared by paralogues (actin B vs G).
                # Take the accession carrying most of the evidence, but only
                # when it clearly leads — a near-tie stays ambiguous.
                counts = Counter(accs).most_common()
                if len(counts) > 1 and counts[0][1] < counts[1][1] * 2:
                    continue
                winner = counts[0][0]
                group = [d for d in same if (d.source_accession or "").strip() == winner]

        donor = group[0]
        s.source_accession = donor.source_accession
        if not (s.source_gene or "").strip() and donor.source_gene:
            s.source_gene = donor.source_gene
        if not (s.ref_proteins or "").strip() and donor.ref_proteins:
            s.ref_proteins = donor.ref_proteins
        s.annotation_source = "peptide-match"
        updated += 1

    if updated:
        db.commit()
    return updated


def infer_substitution(saap: SAAP) -> str | None:
    """Derive the 'X to Y' AAS string from the peptide pair.

    Imports sometimes omit the AAS column. Since the base and substituted
    peptides differ at exactly one residue, the substitution can be recovered
    directly.

    A cell may list several candidate peptides ('SAVSGLWGK;ASVSGLWGK') when the
    search could not decide which position carries the change; the result is
    used only when every pairing agrees on the same residue substitution.
    Returns None when the pair doesn't resolve to a single change.
    """
    from .ingest import _derive_aa_sub

    derived = _derive_aa_sub(saap.mtp_seq or "", saap.bp_seq or "")
    return derived or None


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

    Annotating the whole database (no `ids`) is a full refresh: every SAAP is
    re-fetched and its gene, protein description, accession and Ensembl fields
    are rewritten from UniProt, so the whole table ends up in one consistent
    format rather than a mix of imported and resolved values.
    """
    if not ids and only_missing and not overwrite:
        only_missing = False
        overwrite = True
    # Every SAAP is eligible: those without an accession are resolved from their
    # gene symbol, so import no longer has to discard them.
    stmt = select(SAAP)
    if ids:
        stmt = stmt.where(SAAP.id.in_(ids))
    if only_missing and not overwrite:
        stmt = stmt.where(or_(SAAP.ensembl_gene.is_(None),
                              SAAP.position_in_protein.is_(None),
                              SAAP.protein_sequence.is_(None),
                              SAAP.source_accession.is_(None),
                              SAAP.source_accession == "",
                              SAAP.aa_sub.is_(None),
                              SAAP.aa_sub == ""))
    stmt = stmt.order_by(SAAP.id)
    if limit:
        stmt = stmt.limit(limit)

    saaps = list(db.scalars(stmt).all())
    result = AnnotationResult(requested=len(saaps))
    if not saaps:
        return result

    # Fill the AAS column where the import omitted it — derivable offline from
    # the peptide pair, so it happens regardless of network access.
    merged_ids: set[int] = set()
    for s_ in saaps:
        if not (s_.aa_sub or "").strip():
            inferred = infer_substitution(s_)
            if not inferred:
                continue
            # The completed row may already exist (same peptide imported twice,
            # once without its annotation columns) — merge instead of colliding.
            if _merge_duplicate(db, s_, inferred):
                merged_ids.add(s_.id)
                result.merged_duplicates += 1
            else:
                s_.aa_sub = inferred
                result.aas_filled += 1
    if result.aas_filled or merged_ids:
        db.commit()
    if merged_ids:
        saaps = [x for x in saaps if x.id not in merged_ids]

    # Rows imported without gene/UniProt columns: copy identifiers from the same
    # base peptide elsewhere in the database before spending any network calls.
    result.resolved_by_peptide = resolve_from_database(db, saaps)

    by_acc: dict[str, list[SAAP]] = {}
    by_gene: dict[str, list[SAAP]] = {}
    by_peptide: dict[tuple[str, str | None], list[SAAP]] = {}
    saap_candidates: dict[int, list[str]] = {}
    for s in saaps:
        target = _saap_species_name(db, s)
        # A stored gene from the wrong species means the whole annotation
        # (accession included) points at the wrong proteome. Re-resolve from the
        # gene symbol in the correct species rather than trusting the accession.
        if not annotation_matches_species(s, target):
            wrong = _first_gene(s.source_gene)
            if wrong:
                result.species_corrected += 1
                by_gene.setdefault((wrong.upper(), TAXA.get(target or "")), []).append(s)
                continue

        cands = _accession_candidates(s.source_accession)
        if cands:
            saap_candidates[s.id] = cands
            for acc in cands:
                by_acc.setdefault(acc, []).append(s)
            continue
        # No accession: fall back to the gene symbol.
        gene = _first_gene(s.source_gene)
        if gene:
            # Key by (gene, organism): the same symbol exists in several species.
            by_gene.setdefault((gene.upper(), _saap_organism(db, s)), []).append(s)
        elif (s.bp_seq or "").strip():
            # No accession and no gene: the peptide sequence is all we have.
            by_peptide.setdefault(
                ((s.bp_seq or "").strip().upper(), _saap_organism(db, s)), []).append(s)
        else:
            result.no_identifier += 1
            s.annotation_source = s.annotation_source or "no-identifier"

    # Canonical entries come from the batched search endpoint; isoform-specific
    # ids ('O00429-2') must be fetched individually.
    canonical = [a for a in by_acc if "-" not in a]
    isoform_ids = [a for a in by_acc if "-" in a]

    records, errors = fetch_uniprot(
        canonical, batch_size=batch_size, timeout=timeout, session=session
    )
    result.errors.extend(errors)
    iso_records, iso_errors = fetch_isoforms(
        isoform_ids, timeout=timeout, session=session
    )
    records.update(iso_records)
    result.errors.extend(iso_errors)

    for s in saaps:
        cands = saap_candidates.get(s.id)
        if not cands:
            continue
        # Prefer the accession whose sequence actually contains the peptide;
        # fall back to the first entry that resolved at all.
        chosen: ProteinRecord | None = None
        fallback: ProteinRecord | None = None
        for acc in cands:
            rec = records.get(acc)
            if rec is None:
                continue
            if fallback is None:
                fallback = rec
            if rec.sequence and contains_peptide(rec.sequence, s.bp_seq):
                chosen = rec
                break
        rec = chosen or fallback
        if rec is None:
            if errors or iso_errors:
                result.failed += 1
            else:
                result.not_found += 1
                s.annotation_source = s.annotation_source or "uniprot:not-found"
            continue
        result.resolved += 1
        apply_record(s, rec, overwrite=overwrite)
        if s.position_in_protein is not None:
            result.positioned += 1
            s.annotation_source = "uniprot"
        else:
            result.unmatched_peptide += 1
            s.annotation_source = "uniprot:unmatched"
            if len(result.unmatched_examples) < 10:
                result.unmatched_examples.append(
                    f"{s.source_accession or s.source_gene}: {s.bp_seq} not in "
                    f"{rec.accession} ({len(rec.sequence or '')} aa)")

    # Resolve the accession-less SAAPs by gene symbol and back-fill the
    # accession itself, so they behave like any other row from now on.
    if by_gene:
        gene_records: dict[str, ProteinRecord] = {}
        gene_errors: list[str] = []
        # One call per organism so each gene is queried in the right species.
        by_organism: dict[str | None, list[str]] = {}
        for (gene, organism) in by_gene:
            by_organism.setdefault(organism, []).append(gene)
        for organism, genes in by_organism.items():
            recs, errs = fetch_uniprot_by_gene(
                genes, organism_id=organism, timeout=timeout, session=session
            )
            for g, r in recs.items():
                gene_records[(g, organism)] = r
            gene_errors.extend(errs)
        result.errors.extend(gene_errors)
        for key, group in by_gene.items():
            rec = gene_records.get(key)
            if rec is None:
                if gene_errors:
                    result.failed += len(group)
                else:
                    result.not_found += len(group)
                    for s in group:
                        s.annotation_source = s.annotation_source or "uniprot:gene-not-found"
                continue
            for s in group:
                # Fill a missing accession, and replace one that came from the
                # wrong species — the lookup that produced `rec` was constrained
                # to this SAAP's correct organism.
                if rec.accession and (not s.source_accession
                                      or not annotation_matches_species(
                                          s, _saap_species_name(db, s))):
                    s.source_accession = rec.accession
                    s.source_gene = rec.gene or s.source_gene
                result.resolved += 1
                result.resolved_by_gene += 1
                apply_record(s, rec, overwrite=overwrite)
                if s.position_in_protein is not None:
                    result.positioned += 1
                    s.annotation_source = "uniprot:gene"
                else:
                    result.unmatched_peptide += 1
                    s.annotation_source = "uniprot:gene-unmatched"

    # Last resort: search UniProt by the peptide sequence itself.
    if by_peptide:
        seq_records: dict[tuple[str, str | None], ProteinRecord] = {}
        # Track which organism batches actually errored, so a single failed
        # peptide doesn't get every other peptide marked as "failed".
        errored_orgs: set[str | None] = set()
        by_org: dict[str | None, list[str]] = {}
        for (pep, organism) in by_peptide:
            by_org.setdefault(organism, []).append(pep)
        for organism, peps in by_org.items():
            recs, errs = fetch_by_peptide(
                peps, organism_id=organism, timeout=timeout, session=session
            )
            for pep, rec in recs.items():
                seq_records[(pep, organism)] = rec
            if errs:
                errored_orgs.add(organism)
            result.errors.extend(errs)
        for key, group in by_peptide.items():
            rec = seq_records.get(key)
            if rec is None:
                # Distinguish a real "no match" from a service/network failure,
                # and ALWAYS write a marker so the row never stays unattempted
                # (annotation_source is how the UI reports what happened).
                if key[1] in errored_orgs:
                    result.failed += len(group)
                    for s in group:
                        s.annotation_source = "uniprot:peptide-error"
                else:
                    result.not_found += len(group)
                    for s in group:
                        s.annotation_source = "uniprot:peptide-not-found"
                continue
            for s in group:
                if rec.accession and not s.source_accession:
                    s.source_accession = rec.accession
                result.resolved += 1
                result.resolved_by_sequence += 1
                apply_record(s, rec, overwrite=overwrite)
                if s.position_in_protein is not None:
                    result.positioned += 1
                    s.annotation_source = "uniprot:peptide"
                else:
                    result.unmatched_peptide += 1
                    s.annotation_source = "uniprot:peptide-unmatched"

    db.commit()
    return result
