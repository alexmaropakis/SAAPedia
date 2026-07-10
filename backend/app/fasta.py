"""Generate FragPipe/Philosopher-safe mock-UniProt FASTA for each SAAP's variant
(substituted / mistranslated) peptide.

Header format (matches the lab's build pipeline):

    >sp|{accession}-MTP{id}-{token}|{gene}-mut {gene} mistranslated MTP{id} \
       OS={Species} OX={taxid} GN={gene} PE=1 SV=1

The accession component is made unique per entry by the internal SAAP id (MTP{id}),
and `token` is the pool/plex label chosen at export time.
"""
from __future__ import annotations

import re

from .models import SAAP

DEFAULT_HEADER = (
    ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} mistranslated {mid} "
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


def _resolve_species(species_raw: str) -> tuple[str, str]:
    """Return (OS name, OX taxid). Uses the first species if combined; unknown
    species keep their name with a blank taxid."""
    first = (species_raw or "").split("/")[0].strip().lower()
    if first in _SPECIES_INFO:
        return _SPECIES_INFO[first]
    return (first.capitalize() if first else "", "")


def _fields(saap: SAAP, species: str, token: str, seq_no: int) -> dict:
    accession = saap.source_accession or f"SAAP{seq_no}"
    gene = saap.source_gene or "-"
    os_name, ox = _resolve_species(species)
    mid = f"MTP{seq_no}"
    return {
        "id": seq_no,
        "mid": mid,
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
        "protein": saap.ref_proteins or "",
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
    reference_fasta: str | None = None,
    line_width: int = DEFAULT_LINE_WIDTH,
    header_template: str = DEFAULT_HEADER,
) -> str:
    species_by_id = species_by_id or {}
    token_by_id = token_by_id or {}

    # Forward (target) entries: the SAAP variant peptides. The MTP number is a
    # 1-based export-order index (not the DB id), so it stays contiguous and its
    # max equals the entry count. Note: a given SAAP's number can therefore differ
    # between exports as the selected set changes.
    entries: list[tuple[str, str]] = []
    for seq_no, saap in enumerate(saaps, start=1):
        species = species_by_id.get(saap.id) or default_species
        tok = sanitize_token(token_by_id.get(saap.id) or token)
        fields = _fields(saap, species, tok, seq_no)
        try:
            header = header_template.format(**fields)
        except (KeyError, IndexError, ValueError):
            # Bad custom template -> fall back to the default so export never fails.
            header = DEFAULT_HEADER.format(**fields)
        if not header.startswith(">"):
            header = ">" + header
        entries.append((header, saap.mtp_seq))

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
