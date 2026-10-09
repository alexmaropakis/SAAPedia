"""
Digest-aware cleavage-site check for a substitution.

This is deliberately independent of the `missed_cleavage` /
`aas_at_peptide_terminus` columns that ship with imported files — how those
were originally computed is opaque and not necessarily consistent across
source datasets. Instead this derives the answer directly from what we
actually have on hand: the full protein sequence, where the substituted
residue sits in it, and which protease(s) generated the peptide.

A substitution counts as sitting "at a cleavage site" for a given enzyme when
the residue at that position — either the reference (BP) residue or the
substituted (SAAP) residue — is one the enzyme cuts on, and the neighboring
residue doesn't block the cut (e.g. trypsin does not cut before a proline).
Either direction is flagged: a wrong amino-acid call can just as easily
manufacture a new cleavage-defined peptide boundary as erase a real one, and
both are the same artifact mechanism.
"""
from __future__ import annotations
from typing import Optional

# side: "C" = enzyme cuts C-terminal to a `cut` residue (the bond between it
# and the next residue); "N" = cuts N-terminal to a `cut` residue (the bond
# between the previous residue and it). `block` residues on the far side of
# that bond prevent the cut (e.g. trypsin does not cut before a proline).
ENZYME_RULES: dict[str, dict] = {
    "trypsin":      {"side": "C", "cut": set("KR"),   "block": set("P")},
    "lysc":         {"side": "C", "cut": set("K"),    "block": set()},
    "lysn":         {"side": "N", "cut": set("K"),    "block": set()},
    "argc":         {"side": "C", "cut": set("R"),    "block": set()},
    "chymotrypsin": {"side": "C", "cut": set("FYWL"), "block": set("P")},
    "gluc":         {"side": "C", "cut": set("E"),    "block": set()},
    "aspn":         {"side": "N", "cut": set("D"),    "block": set()},
}


def _normalize_digest(digest: str | None) -> Optional[str]:
    if not digest:
        return None
    key = "".join(ch for ch in digest.lower() if ch.isalnum())
    return key or None


def _cuts_here(protein_seq: str, idx0: int, residue: str, rule: dict) -> bool:
    """Would `rule`'s enzyme cut adjacent to `residue` sitting at 0-based
    idx0, given the surrounding sequence?"""
    if residue not in rule["cut"]:
        return False
    if rule["side"] == "C":
        neighbor = protein_seq[idx0 + 1] if idx0 + 1 < len(protein_seq) else None
    else:
        neighbor = protein_seq[idx0 - 1] if idx0 - 1 >= 0 else None
    return not (neighbor and neighbor in rule["block"])


def is_cleavage_position(
    protein_seq: str | None,
    position_in_protein: int | None,
    digest: str | None,
    ref_aa: str | None,
    sub_aa: str | None,
) -> Optional[bool]:
    """True/False if this single enzyme can be evaluated at this position,
    else None — missing protein sequence/position, or an enzyme we have no
    rule for (blank digest, or something non-specific)."""
    rule = ENZYME_RULES.get(_normalize_digest(digest))
    if rule is None or not protein_seq or not position_in_protein:
        return None
    idx0 = position_in_protein - 1
    if not (0 <= idx0 < len(protein_seq)):
        return None
    ref_aa = (ref_aa or protein_seq[idx0]).upper()
    hits_as_ref = _cuts_here(protein_seq, idx0, ref_aa, rule)
    hits_as_sub = bool(sub_aa) and _cuts_here(protein_seq, idx0, sub_aa.upper(), rule)
    return hits_as_ref or hits_as_sub


def is_cleavage_position_any(
    protein_seq: str | None,
    position_in_protein: int | None,
    digests_csv: str | None,
    aa_sub: str | None,
) -> Optional[bool]:
    """Same as `is_cleavage_position`, but for a SAAP observed under several
    digests (a comma-separated list, as rolled up per SAAP). True if it's a
    cleavage-site substitution under ANY of them; None if none of the listed
    digests are ones we can evaluate."""
    from .annotate import parse_substitution  # local import: avoid an import cycle

    if not digests_csv:
        return None
    ref_aa, sub_aa = parse_substitution(aa_sub)
    evaluated = False
    for digest in digests_csv.split(","):
        result = is_cleavage_position(protein_seq, position_in_protein, digest.strip(), ref_aa, sub_aa)
        if result is None:
            continue
        evaluated = True
        if result:
            return True
    return False if evaluated else None
