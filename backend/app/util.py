# Utility functions

from __future__ import annotations
import re

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


# Immunoglobulin genes: V/D/J segments and constant regions of the heavy,
# kappa and lambda loci (human IGHV3-23, IGKC, IGHG1; mouse Ighv1-72, Igkc,
# Ighg2b), plus the J chain. Anchored so look-alikes such as IGHMBP2 or
# IGLON5 (not immunoglobulins) are not matched.
_IG_GENE = re.compile(
    r"^(IG[HKL][VDJ]\d[\w-]*|IGHG[1-4P]?|IGHG2[ABC]|IGHA[12]?|IGHM|IGHD|IGHE|IGKC|IGLC\d*|IGLL\d|JCHAIN|IGJ)$",
    re.IGNORECASE)
_IG_DESCRIPTION = re.compile(r"^(Immunoglobulin (heavy|kappa|lambda|J)\b|Ig (gamma|kappa|lambda|mu|alpha)\b)", re.IGNORECASE)


def is_immunoglobulin(gene: str | None, description: str | None = None) -> bool:
    genes = [g.strip() for g in (gene or "").replace(",", ";").split(";") if g.strip()]
    return any(_IG_GENE.match(g) for g in genes) or bool(_IG_DESCRIPTION.match((description or "").strip()))
