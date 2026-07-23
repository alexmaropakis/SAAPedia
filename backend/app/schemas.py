"""Pydantic request/response models."""
from __future__ import annotations

from typing import Optional

from pydantic import BaseModel


class ExportRequest(BaseModel):
    # Provide either explicit SAAP ids, or a filter block to export everything
    # matching the current view. ids take precedence when both are present.
    ids: Optional[list[int]] = None
    filters: Optional[dict] = None
    species: str = ""          # fallback only; species is taken from the data
    token: str = ""            # pool/plex label used in the header accession
    decoys: bool = False       # append rev_ reversed-sequence decoys
    base_peptides: bool = False  # also emit each SAAP's base peptide (BP)
    # "peptide" = substituted peptide sequences (default);
    # "protein" = full-length reference protein with the substitution applied,
    # which is what a multi-digest search needs.
    entry_mode: str = "peptide"
    line_width: int = 60
    header_template: Optional[str] = None
    base_header_template: Optional[str] = None


class AnnotateRequest(BaseModel):
    """Resolve Ensembl IDs and substitution positions from UniProt."""
    ids: Optional[list[int]] = None   # restrict to these SAAP (default: all)
    only_missing: bool = True         # skip rows already annotated
    overwrite: bool = False           # replace existing values
    limit: Optional[int] = None       # cap rows processed (trial runs)
