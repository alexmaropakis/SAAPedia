"""Generate FragPipe/Philosopher-safe mock-UniProt FASTA for each SAAP's
substituted peptide.

Header format (matches the lab's build pipeline):

    >sp|{accession}-SAAP{id}-{token}|{gene}-mut {gene} substituted SAAP{id} \
       OS={Species} OX={taxid} GN={gene} PE=1 SV=1

The accession component is made unique per entry by the internal SAAP id (SAAP{id}),
and `token` is the pool/plex label chosen at export time.

Terminology: these peptides are always described as "substituted" — never
"mistranslated". Substitution is the observation; mistranslation is only one of
several possible mechanisms behind it, so the neutral term is used throughout
(headers, templates, and the base-peptide entries below).

Base peptides can optionally be emitted alongside the substituted peptides
(`include_base_peptides`), using BASE_HEADER so the two are easy to tell apart
downstream.
"""
from __future__ import annotations
import re
from .models import SAAP

DEFAULT_HEADER = (
    ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

# Header for the unmodified base peptide (BP) of a SAAP. Kept parallel to
# DEFAULT_HEADER but marked "base" and given a BP{id} identifier so substituted
# and base entries never collide in a search database.
BASE_HEADER = (
    ">sp|{accession}-{bid}-{tok}|{gene}-base {gene} base peptide {bid} "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

# Header for a full-length protein carrying one substitution. Used when
# entry_mode="protein": the sequence is the whole reference protein with the
# SAAP's substitution applied in place, so any protease's peptides are
# searchable — not just the one peptide that was observed.
PROTEIN_HEADER = (
    ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} {sub_compact}@{position} "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

# Header for the unmodified reference protein, emitted once per accession
# alongside the variant proteins in protein mode.
PROTEIN_BASE_HEADER = (
    ">sp|{accession}|{gene}-base {gene} reference protein "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

DEFAULT_LINE_WIDTH = 60

# species (lowercased, first token if combined) -> (OS name, OX taxonomy id)
_SPECIES_INFO = {
    "homo sapiens": ("Homo sapiens", "9606"),
    "human": ("Homo sapiens", "9606"),
    "mus musculus": ("Mus musculus", "10090"),
    "mouse": ("Mus musculus", "10090"),
    "rattus norvegicus": ("Rattus norvegicus", "10116"),
    "rat": ("Rattus norvegicus", "10116"),
    "saccharomyces cerevisiae": ("Saccharomyces cerevisiae", "559292"),
}

_SUB_SPLIT = re.compile(r"\s*(?:to|->|>|/|→)\s*", re.IGNORECASE)


def compact_sub(aa_sub: str | None) -> str:
    """'V to P' -> 'V2P'; falls back to alphanumerics of the raw value."""
    if not aa_sub:
        return "sub"
    parts = [p for p in _SUB_SPLIT.split(aa_sub.strip()) if p]
    if len(parts) == 2:
        return f"{parts[0]}2{parts[1]}"
    return re.sub(r"[^A-Za-z0-9]", "", aa_sub) or "sub"


def sanitize_token(token: str | None) -> str:
    """Lowercase, keep alphanumerics/underscores (matches the pipeline's plex token)."""
    return re.sub(r"[^A-Za-z0-9]+", "_", (token or "").strip().lower()).strip("_")


def first_value(value: str | None) -> str:
    """Gene and protein fields can hold several ';'-separated entries (with leading
    or empty segments). A FASTA header may name only one, so use the first
    non-empty entry. Note: ',' and '/' occur *within* a single name and are kept."""
    if not value:
        return ""
    for part in value.split(";"):
        part = part.strip()
        if part:
            return part
    return ""


def _resolve_species(species_raw: str) -> tuple[str, str]:
    """Return (OS name, OX taxid). Uses the first species if combined; unknown
    species keep their name with a blank taxid."""
    first = (species_raw or "").split("/")[0].strip().lower()
    if first in _SPECIES_INFO:
        return _SPECIES_INFO[first]
    return (first.capitalize() if first else "", "")


def _fields(saap: SAAP, species: str, token: str, seq_no: int) -> dict:
    accession = saap.source_accession or f"SAAP{seq_no}"
    gene = first_value(saap.source_gene) or "-"
    os_name, ox = _resolve_species(species)
    mid = f"SAAP{seq_no}"
    return {
        "id": seq_no,
        "mid": mid,
        "bid": f"BP{seq_no}",          # identifier for the base-peptide entry
        "accession": accession,
        "tok": token,
        "gene": gene,
        "species": os_name,
        "taxid": ox,
        "entry_name": f"{gene}-mut",
        "sub_compact": compact_sub(saap.aa_sub),
        "aa_sub": saap.aa_sub or "?",
        "bp_seq": saap.bp_seq or "-",
        "mtp_seq": saap.mtp_seq,
        "protein": first_value(saap.ref_proteins),
        # Ensembl / positional annotation (blank when not yet annotated), so
        # custom header templates can reference them.
        "ensembl_gene": saap.ensembl_gene or "",
        "ensembl_transcript": saap.ensembl_transcript or "",
        "ensembl_protein": saap.ensembl_protein or "",
        "position": saap.position_in_protein if saap.position_in_protein is not None else "",
        "protein_description": saap.protein_description or "",
    }


def parse_fasta(text: str) -> list[tuple[str, str]]:
    """Parse FASTA text into [(header, sequence)]; header keeps its leading '>'."""
    entries: list[tuple[str, str]] = []
    header, seq = None, []
    for line in text.splitlines():
        line = line.rstrip()
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq)))
            header, seq = line, []
        elif header is not None:
            seq.append(line.strip())
    if header is not None:
        entries.append((header, "".join(seq)))
    return entries


def _wrap(seq: str, width: int) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def generate_fasta(
    saaps: list[SAAP],
    *,
    species_by_id: dict[int, str] | None = None,
    default_species: str = "",
    token: str = "",
    token_by_id: dict[int, str] | None = None,
    include_decoys: bool = False,
    include_base_peptides: bool = False,
    entry_mode: str = "peptide",
    reference_fasta: str | None = None,
    line_width: int = DEFAULT_LINE_WIDTH,
    header_template: str = DEFAULT_HEADER,
    base_header_template: str = BASE_HEADER,
    skipped: list[str] | None = None,
) -> str:
    """Build the FASTA text.

    entry_mode
      "peptide" (default) — emit the substituted peptide sequence itself. Fine
          when the search uses the same protease the SAAPs were observed with.
      "protein" — emit the full-length reference protein with the substitution
          applied at its position, one entry per SAAP, plus the unmodified
          reference protein once per accession. Use this for multi-digest
          searches: the variant residue is then reachable by whatever peptides
          each protease produces, not only the originally observed peptide.
          Requires annotation (position + cached sequence); SAAPs lacking it are
          skipped and reported via `skipped`.

    include_base_peptides — peptide mode only; in protein mode the unmodified
        reference protein already plays that role.

    skipped — optional list that collects a human-readable reason for every
        SAAP omitted in protein mode.
    """
    species_by_id = species_by_id or {}
    token_by_id = token_by_id or {}
    skipped = skipped if skipped is not None else []

    def _render(template: str, fields: dict, fallback: str) -> str:
        try:
            header = template.format(**fields)
        except (KeyError, IndexError, ValueError):
            # Bad custom template -> fall back to the default so export never fails.
            header = fallback.format(**fields)
        return header if header.startswith(">") else ">" + header

    entries: list[tuple[str, str]] = []

    if entry_mode == "protein":
        from .annotate import apply_substitution

        # Reference proteins are emitted once per accession, even though many
        # SAAPs may map to the same protein.
        seen_reference: set[str] = set()
        for seq_no, saap in enumerate(saaps, start=1):
            species = species_by_id.get(saap.id) or default_species
            tok = sanitize_token(token_by_id.get(saap.id) or token)
            fields = _fields(saap, species, tok, seq_no)

            variant, error = apply_substitution(saap)
            if error:
                skipped.append(f"SAAP {saap.id} ({saap.mtp_seq}): {error}")
                continue

            acc = saap.source_accession or ""
            if acc and acc not in seen_reference and saap.protein_sequence:
                seen_reference.add(acc)
                entries.append((_render(PROTEIN_BASE_HEADER, fields, PROTEIN_BASE_HEADER),
                                saap.protein_sequence))
            entries.append((_render(header_template, fields, PROTEIN_HEADER), variant))
    else:
        # Forward (target) entries: the substituted peptides. The SAAP number is a
        # 1-based export-order index (not the DB id), so it stays contiguous and its
        # max equals the entry count. Note: a given SAAP's number can therefore differ
        # between exports as the selected set changes.
        seen_base: set[str] = set()
        for seq_no, saap in enumerate(saaps, start=1):
            species = species_by_id.get(saap.id) or default_species
            tok = sanitize_token(token_by_id.get(saap.id) or token)
            fields = _fields(saap, species, tok, seq_no)
            entries.append((_render(header_template, fields, DEFAULT_HEADER), saap.mtp_seq))

            if include_base_peptides:
                bp = (saap.bp_seq or "").strip()
                # Skip when absent, unchanged, or already emitted: several SAAPs can
                # share one base peptide, and duplicate FASTA entries break some
                # search engines' protein inference.
                if bp and bp != saap.mtp_seq and bp not in seen_base:
                    seen_base.add(bp)
                    entries.append((_render(base_header_template, fields, BASE_HEADER), bp))

    # ... then the reference proteome (if supplied), passed through unchanged.
    if reference_fasta:
        entries.extend(parse_fasta(reference_fasta))

    out: list[str] = []
    for header, seq in entries:
        out.append(header)
        out.append(_wrap(seq, line_width))
    if include_decoys:
        # '>rev_' + original header, sequence reversed, for EVERY target (SAAP
        # and reference), appended after all forwards (matches the pipeline).
        for header, seq in entries:
            out.append(">rev_" + header[1:])
            out.append(_wrap(seq[::-1], line_width))
    return "\n".join(out) + ("\n" if out else "")
