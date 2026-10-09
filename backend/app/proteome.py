"""
Exact match of substituted peptides against reference proteome FASTAs.

A SAAP should not be encoded anywhere in the genome; if its sequence occurs
verbatim in another protein (or isoform) of the same species, the apparent
substitution is better explained by that protein. Matching is I = L, since the
two are isobaric.

Reference FASTAs (UniProt format, OS= in headers) are read from
SAAP_REFERENCE_DIR (default backend/reference/). Without them the check is
skipped and `proteome_hits` stays NULL.

References
----------
The UniProt Consortium. (2025). UniProt: The Universal Protein Knowledgebase in
    2025. Nucleic Acids Research, 53(D1), D609–D617.
    https://doi.org/10.1093/nar/gkae1010
"""
from __future__ import annotations

import bisect
import os
import re
from functools import lru_cache
from pathlib import Path

from sqlalchemy import select
from sqlalchemy.orm import Session

from .models import Observation, SAAP

REFERENCE_DIR = Path(os.environ.get(
    "SAAP_REFERENCE_DIR", Path(__file__).resolve().parent.parent / "reference"))
MAX_HITS = 50
_OS = re.compile(r"\bOS=(.+?)(?:\s+OX=|\s+GN=|\s+PE=|$)")
_GN = re.compile(r"\bGN=(\S+)")


class Reference:
    """One species' proteome as a single I->L normalized string."""

    def __init__(self):
        self.parts: list[str] = []
        self.starts: list[int] = []
        self.accessions: list[str] = []
        self.meta: dict[str, tuple[str | None, bool]] = {}  # accession -> (gene, reviewed)
        self._len = 0
        self.text = ""

    def add(self, accession: str, seq: str, gene: str | None = None, reviewed: bool = False):
        self.meta[accession] = (gene, reviewed)
        self.starts.append(self._len)
        self.accessions.append(accession)
        self.parts.append(seq.upper().replace("I", "L") + "\n")
        self._len += len(seq) + 1

    def sequence(self, accession: str) -> str | None:
        """I->L normalized sequence of one entry."""
        try:
            i = self.accessions.index(accession)
        except ValueError:
            return None
        end = self.starts[i + 1] - 1 if i + 1 < len(self.starts) else len(self.text) - 1
        return self.text[self.starts[i]:end]

    def find(self, peptide: str) -> list[str]:
        pep, hits, i = peptide.upper().replace("I", "L"), [], -1
        while len(hits) < MAX_HITS:
            i = self.text.find(pep, i + 1)
            if i < 0:
                break
            acc = self.accessions[bisect.bisect_right(self.starts, i) - 1]
            if acc not in hits:
                hits.append(acc)
        return hits


@lru_cache(maxsize=1)
def references() -> dict[str, Reference]:
    refs: dict[str, Reference] = {}
    for path in sorted(REFERENCE_DIR.glob("*.fa*")):
        header, seq = None, []

        def flush():
            if header:
                m, gn = _OS.search(header), _GN.search(header)
                acc = header[1:].split("|")[1] if header.count("|") >= 2 else header[1:].split()[0]
                refs.setdefault((m.group(1) if m else "").strip().lower(), Reference()).add(
                    acc, "".join(seq), gn.group(1) if gn else None, header.startswith(">sp|"))

        with path.open() as fh:
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    flush()
                    header, seq = line, []
                elif line:
                    seq.append(line)
        flush()
    for ref in refs.values():
        ref.text, ref.parts = "".join(ref.parts), []
    return refs


def match(peptide: str, species: set[str]) -> list[str] | None:
    """Accessions containing `peptide` in the given species' proteomes, or
    None when no reference is available for any of them."""
    refs = references()
    usable = [refs[s] for s in species if s in refs]
    if not usable:
        return None
    return [a for ref in usable for a in ref.find(peptide)][:MAX_HITS]


def resolve_unidentified(db: Session, saaps: list[SAAP], species_of) -> int:
    """Give SAAPs with no accession and no gene an identity by locating their
    base peptide in the reference proteome. Accepted only when every hit
    belongs to one gene; Swiss-Prot canonical entries are preferred."""
    if not references():
        return 0
    n = 0
    for s in saaps:
        if (s.source_accession or "").strip() or (s.source_gene or "").strip() or not s.bp_seq:
            continue
        for sp in species_of(s):
            ref = references().get(sp)
            if ref is None:
                continue
            hits = ref.find(s.bp_seq)
            genes = {ref.meta[a][0] for a in hits if ref.meta[a][0]}
            if not hits or len(genes) != 1:
                continue
            best = sorted(hits, key=lambda a: (not ref.meta[a][1], "-" in a))[0]
            s.source_accession, s.source_gene = best, genes.pop()
            s.annotation_source = "reference-match"
            n += 1
            break
    db.commit()
    return n


def describe(hits: str | None) -> list[dict] | None:
    """[{accession, gene, reviewed}] for a stored proteome_hits value."""
    if hits is None:
        return None
    meta = {a: m for ref in references().values() for a, m in ref.meta.items()}
    return [{"accession": a, "gene": (meta.get(a) or (None,))[0], "reviewed": bool((meta.get(a) or (None, False))[1])}
            for a in hits.split(",") if a]


def check_all(db: Session, *, overwrite: bool = False) -> int:
    """Fill `proteome_hits` for every SAAP (unchecked ones only, by default)."""
    if not references():
        return 0
    stmt = select(SAAP)
    if not overwrite:
        stmt = stmt.where(SAAP.proteome_hits.is_(None))
    saaps = db.scalars(stmt).all()
    species: dict[int, set[str]] = {}
    for sid, sp in db.execute(select(Observation.saap_id, Observation.species).distinct()):
        if sp:
            species.setdefault(sid, set()).add(sp.strip().lower())
    n = 0
    for s in saaps:
        hits = match(s.mtp_seq or "", species.get(s.id, set()))
        if hits is not None:
            s.proteome_hits = ",".join(hits)
            n += 1
    db.commit()
    return n
