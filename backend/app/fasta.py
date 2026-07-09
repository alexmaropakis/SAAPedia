"""Generate UniProt-style FASTA for the variant (substituted) peptide of each SAAP.

Each selected SAAP emits one entry whose sequence is its `mtp_seq`. The header
follows UniProt conventions:

    >db|ACCESSION|ENTRY_NAME description OS=<species> GN=<gene> ...

Because the sequence is the substituted peptide (not the full protein), the
accession/entry-name are derived from the source protein plus the substitution,
and the internal SAAP id guarantees uniqueness.
"""
from __future__ import annotations

import re

from .models import SAAP

DEFAULT_HEADER = (
    ">saap|{accession}|{entry_name} {protein} SAAP variant ({aa_sub}) "
    "OS={species} GN={gene} BP={bp_seq}"
)
DEFAULT_LINE_WIDTH = 60
DEFAULT_SPECIES = "Homo sapiens"

_SUB_SPLIT = re.compile(r"\s*(?:to|->|>|/|→)\s*", re.IGNORECASE)


def compact_sub(aa_sub: str | None) -> str:
    """'V to P' -> 'V2P'; falls back to alphanumerics of the raw value."""
    if not aa_sub:
        return "sub"
    parts = [p for p in _SUB_SPLIT.split(aa_sub.strip()) if p]
    if len(parts) == 2:
        return f"{parts[0]}2{parts[1]}"
    return re.sub(r"[^A-Za-z0-9]", "", aa_sub) or "sub"


def _fields(saap: SAAP, species: str) -> dict:
    accession = saap.source_accession or f"SAAP{saap.id}"
    sub_c = compact_sub(saap.aa_sub)
    return {
        "id": saap.id,
        "accession": accession,
        "entry_name": f"SAAP{saap.id}_{sub_c}",
        "sub_compact": sub_c,
        "gene": saap.source_gene or "-",
        "protein": saap.ref_proteins or "Uncharacterized SAAP",
        "aa_sub": saap.aa_sub or "?",
        "bp_seq": saap.bp_seq or "-",
        "mtp_seq": saap.mtp_seq,
        "species": species,
    }


def _wrap(seq: str, width: int) -> str:
    if width <= 0:
        return seq
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def generate_fasta(
    saaps: list[SAAP],
    *,
    species_by_id: dict[int, str] | None = None,
    default_species: str = "",
    line_width: int = DEFAULT_LINE_WIDTH,
    header_template: str = DEFAULT_HEADER,
) -> str:
    species_by_id = species_by_id or {}
    out: list[str] = []
    for saap in saaps:
        species = species_by_id.get(saap.id) or default_species
        fields = _fields(saap, species)
        try:
            header = header_template.format(**fields)
        except (KeyError, IndexError, ValueError):
            # Bad custom template -> fall back to the default so export never fails.
            header = DEFAULT_HEADER.format(**fields)
        if not header.startswith(">"):
            header = ">" + header
        out.append(header)
        out.append(_wrap(saap.mtp_seq, line_width))
    return "\n".join(out) + ("\n" if out else "")
