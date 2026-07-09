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
    token: str = "pooled"      # pool/plex label used in the header accession
    decoys: bool = False       # append rev_ reversed-sequence decoys
    line_width: int = 60
    header_template: Optional[str] = None
