"""Small shared helpers."""
from __future__ import annotations


def doi_to_url(doi: str | None) -> str | None:
    """Turn a raw DOI or URL into a resolvable https link.

    Accepts '10.1038/xyz', 'doi:10.1038/xyz', or a full 'https://doi.org/...'.
    Returns None for blank input.
    """
    if not doi:
        return None
    d = doi.strip()
    if not d:
        return None
    if d.lower().startswith(("http://", "https://")):
        return d
    if d.lower().startswith("doi:"):
        d = d[4:].strip()
    return f"https://doi.org/{d}"
