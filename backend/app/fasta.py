"""
Generate UniProt-style FASTAs for each SAAP 

Header format:
    >sp|{accession}-SAAP{id}-{tok}|{gene}-mut {gene} substituted SAAP{id} \
        OS={Species} OX={taxid} GN={gene} PE=1 SV=1 

Accession is made unique per entry by internal SAAP ID
Token/tok is pool/plex label chosen at export time 

Base peptides can optionally be emitted alongside SAAPs 
        
"""

from __future__ import annotations
import re
from .models import SAAP

# SAAP header
DEFAULT_HEADER = (
    ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

# Header for the unmodified base peptide (BP) of a SAAP
BASE_HEADER = (
    ">sp|{accession}-{bid}-{tok}|{gene}-base {gene} base peptide {bid} "
    "OS={species} OX={taxid} GN={gene} PE=1 SV=1"
)

# Header for a full-length protein carrying one substitution
# Used when entry_mode="protein"; sequence is whole ref prot with sub applied in place 
# This is useful for multi-protease searches where we want to append the protein, 
# not the tryptic/etc. peptide
# Unmodified ref prot come from reference proteome uploaded by the user
PROTEIN_HEADER = (
    ">sp|{accession}-{mid}-{tok}|{gene}-mut {gene} substituted {mid} {sub_compact}@{position} "
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


def split_peptide_cell(value: str | None) -> list[str]:
    """Alternative sequences held in one peptide cell.

    A row may carry several candidate peptides ('SAVSGLWGK;ASVSGLWGK') when the
    search could not place the substituted residue. A FASTA sequence cannot
    contain a separator, so each alternative becomes its own entry.
    """
    if not value:
        return []
    return [p.strip().upper() for p in str(value).replace(",", ";").split(";") if p.strip()]


def _entry_block_reason(saap: SAAP, species: str) -> str | None:
    """Why this SAAP must not be written to a FASTA, or None if it is fine.

    Two conditions block an entry:
      * no gene — there is no protein identity to assert, so the header would
        carry 'GN=-' and a synthetic accession;
      * the gene symbol belongs to a different species than the entry's OS/OX.

    Both produce records that misidentify the protein, which is worse for a
    search database than the peptide simply being absent.
    """
    if not first_value(saap.source_gene):
        return "no gene — cannot assert a protein identity"
    return _annotation_is_consistent(saap, species)


def _annotation_is_consistent(saap: SAAP, species: str) -> str | None:
    """Reason the annotation does not match `species`, or None if it is fine.

    A SAAP carries one gene/accession, but may be observed in several species.
    Emitting a mouse entry that carries a human gene symbol produces a FASTA
    record whose OS/OX and GN disagree, which corrupts downstream protein
    inference. Symbol casing is the reliable signal: human symbols are
    all-caps (SPTAN1), mouse are title-case (Sptan1).
    """
    from .annotate import gene_species

    gene = first_value(saap.source_gene)
    if not gene:
        return None
    implied = gene_species(gene)
    sp = (species or "").strip().lower()
    if not implied or not sp:
        return None
    if sp.startswith("mus") and implied != "mus musculus":
        return f"gene {gene!r} is a human symbol but the entry is mouse"
    if sp.startswith("homo") and implied != "homo sapiens":
        return f"gene {gene!r} is a mouse symbol but the entry is human"
    return None


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
          applied at its position, one entry per SAAP. Only the substituted
          proteins are written: pair this with your own reference proteome
          (uploaded separately, or concatenated afterwards), which supplies the
          unmodified sequences. Use this for multi-digest searches: the variant
          residue is then reachable by whatever peptides each protease produces,
          not only the originally observed peptide. Requires annotation
          (position + cached sequence); SAAPs lacking it are skipped and
          reported via `skipped`.

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
        from .annotate import apply_substitutions_all

        for seq_no, saap in enumerate(saaps, start=1):
            species = species_by_id.get(saap.id) or default_species
            tok = sanitize_token(token_by_id.get(saap.id) or token)
            reason = _entry_block_reason(saap, species)
            if reason:
                skipped.append(f"SAAP {saap.id} ({saap.mtp_seq}): {reason}")
                continue

            fields = _fields(saap, species, tok, seq_no)

            variants, error = apply_substitutions_all(saap)
            if error:
                skipped.append(f"SAAP {saap.id} ({saap.mtp_seq}): {error}")
                continue

            # Only the substituted protein is emitted. The unmodified reference
            # comes from the reference proteome the user supplies, so writing it
            # here too would duplicate entries and skew protein inference/FDR.
            # A peptide repeating within its protein gives one entry per
            # candidate site, each labelled with its own position.
            for vi, (pos, variant) in enumerate(variants, start=1):
                vfields = dict(fields)
                vfields["position"] = pos
                if len(variants) > 1:
                    vfields["mid"] = f"{fields['mid']}-p{pos}"
                entries.append((_render(header_template, vfields, PROTEIN_HEADER), variant))
    else:
        # Forward (target) entries: the substituted peptides. The SAAP number is a
        # 1-based export-order index (not the DB id), so it stays contiguous and its
        # max equals the entry count. Note: a given SAAP's number can therefore differ
        # between exports as the selected set changes.
        seen_base: set[str] = set()
        for seq_no, saap in enumerate(saaps, start=1):
            species = species_by_id.get(saap.id) or default_species
            tok = sanitize_token(token_by_id.get(saap.id) or token)
            # A FASTA entry asserts a protein identity. Without a gene there is
            # nothing to assert, and the header degrades to 'GN=-' with a fake
            # accession — unusable in a search database, so the entry is skipped.
            reason = _entry_block_reason(saap, species)
            if reason:
                skipped.append(f"SAAP {saap.id} ({saap.mtp_seq}): {reason}")
                continue

            fields = _fields(saap, species, tok, seq_no)
            # A cell may hold several candidate peptides; each becomes its own
            # entry, suffixed -1, -2, ... so the headers stay unique.
            variants = split_peptide_cell(saap.mtp_seq)
            for vi, seq in enumerate(variants, start=1):
                vfields = dict(fields)
                if len(variants) > 1:
                    vfields["mid"] = f"{fields['mid']}-{vi}"
                entries.append((_render(header_template, vfields, DEFAULT_HEADER), seq))

            if include_base_peptides:
                for bi, bp in enumerate(split_peptide_cell(saap.bp_seq), start=1):
                    # Skip when unchanged or already emitted: several SAAPs can
                    # share one base peptide, and duplicate FASTA entries break
                    # some search engines' protein inference.
                    if bp in variants or bp in seen_base:
                        continue
                    seen_base.add(bp)
                    bfields = dict(fields)
                    if len(split_peptide_cell(saap.bp_seq)) > 1:
                        bfields["bid"] = f"{fields['bid']}-{bi}"
                    entries.append((_render(base_header_template, bfields, BASE_HEADER), bp))

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
