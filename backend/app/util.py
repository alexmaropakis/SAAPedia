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


def normalize_species(name: str | None) -> str | None:
    # Normalize casing to the standard binomial form, e.g. "homo sapiens" /
    # "HOMO SAPIENS" / "Homo Sapiens" all become "Homo sapiens", so the same
    # species is never split into separate facet/filter values by casing alone.
    if name is None:
        return None
    cleaned = str(name).strip()
    if not cleaned:
        return None
    lowered = cleaned.lower()
    return lowered[0].upper() + lowered[1:]
