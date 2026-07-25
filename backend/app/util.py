# Utility functions

from __future__ import annotations

def normalize_dataset(name: str | None) -> str | None:
    # Function to capitalize the first letter of a dataset label 
    # if it is not already capitalized
    if name is None:
        return None
    cleaned = str(name).strip()
    if not cleaned:
        return None
    return cleaned[0].upper() + cleaned[1:]


def doi_to_url(doi: str | None) -> str | None:
    # Function to turn raw DOI/URL into a resolvable https
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
