"""
Cross-species recurrence: is a substitution also observed at the equivalent
site of the orthologous protein in the other species (human <-> mouse)?

Orthologues are paired by gene symbol (case-insensitive). A site is mapped onto
the other protein by its sequence context: the ±WINDOW residues around it are
matched against the other sequence (exactly when the proteins are identical,
otherwise the best ungapped placement with >= MIN_IDENTITY). Per SAAP:

  cross_species = "same"  the same substitution is observed in the other
                          species at the equivalent site (or the identical
                          SAAP is observed in both species)
                  "site"  a different substitution at the equivalent site
                  ""      nothing at that site in the other species
                  NULL    not evaluable (no orthologue with SAAPs, or the
                          site could not be mapped)

cross_species_detail is a JSON list of the matching SAAPs in the other species:
{"species", "gene", "accession", "variant" (numbered in that protein), "saap_id"},
or [{"identical": true}] when the same peptide was itself observed in both.
"""
from __future__ import annotations

import json
from collections import defaultdict

from sqlalchemy import select
from sqlalchemy.orm import Session

from .annotate import stored_positions
from .investigate import residues
from .models import Observation, SAAP

SPECIES = ("Homo sapiens", "Mus musculus")
COMMON = {"Homo sapiens": "Human", "Mus musculus": "Mouse"}
WINDOW = 10
MIN_IDENTITY = 0.7


def _dump(matches: list[dict]) -> str:
    """JSON list of matches, de-duplicated by SAAP."""
    seen, out = set(), []
    for m in matches:
        key = m.get("saap_id", "identical")
        if key not in seen:
            seen.add(key)
            out.append(m)
    return json.dumps(out)


def map_position(seq_a: str, pos: int, seq_b: str) -> int | None:
    """Equivalent 1-based position of seq_a[pos] in seq_b, by local context."""
    if seq_a == seq_b:
        return pos
    lo, hi = max(0, pos - 1 - WINDOW), min(len(seq_a), pos + WINDOW)
    win, offset = seq_a[lo:hi], pos - 1 - lo
    i = seq_b.find(win)
    if i >= 0 and seq_b.find(win, i + 1) < 0:
        return i + offset + 1
    best, best_i, tie = len(win) + 1, None, False
    for j in range(len(seq_b) - len(win) + 1):
        d = sum(1 for a, b in zip(win, seq_b[j:j + len(win)]) if a != b)
        if d < best:
            best, best_i, tie = d, j, False
        elif d == best:
            tie = True
    if best_i is None or tie or 1 - best / len(win) < MIN_IDENTITY:
        return None
    return best_i + offset + 1


def check_all(db: Session) -> dict:
    species_of: dict[int, set[str]] = defaultdict(set)
    for sid, sp in db.execute(select(Observation.saap_id, Observation.species).distinct()):
        if sp in SPECIES:
            species_of[sid].add(sp)
    saaps = db.scalars(select(SAAP).where(SAAP.protein_sequence.is_not(None),
                                          SAAP.position_in_protein.is_not(None))).all()
    # gene -> species -> [(species, saap, sequence, positions, ref, alt)]
    events: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for s in saaps:
        gene = (s.source_gene or "").split(";")[0].strip().upper()
        ref, alt = residues(s)
        for sp in species_of.get(s.id, ()):
            events[gene][sp].append((sp, s, s.protein_sequence, stored_positions(s), ref, alt))

    counts = defaultdict(int)
    for s in saaps:
        sps = species_of.get(s.id, set())
        if set(SPECIES) <= sps:
            s.cross_species, s.cross_species_detail = "same", _dump([{"identical": True}])
            counts["same"] += 1
            continue
        gene = (s.source_gene or "").split(";")[0].strip().upper()
        others = [e for sp in SPECIES if sp not in sps for e in events.get(gene, {}).get(sp, [])]
        if not gene or not others:
            s.cross_species = s.cross_species_detail = None
            continue
        ref, alt = residues(s)
        pos = (stored_positions(s) or [None])[0]
        mapped_any, same, site = False, [], []
        cache: dict[str, int | None] = {}
        for osp, o, seq, positions, oref, oalt in others:
            if seq not in cache:
                cache[seq] = map_position(s.protein_sequence, pos, seq) if pos else None
            q = cache[seq]
            if q is None:
                continue
            mapped_any = True
            if q in positions:
                match = {"species": COMMON[osp], "gene": (o.source_gene or "").split(";")[0].strip() or None,
                         "accession": o.protein_accession, "variant": f"{oref}{q}{oalt}", "saap_id": o.id}
                (same if oalt == alt else site).append(match)
        if same:
            s.cross_species, s.cross_species_detail = "same", _dump(same)
        elif site:
            s.cross_species, s.cross_species_detail = "site", _dump(site)
        elif mapped_any:
            s.cross_species, s.cross_species_detail = "", None
        else:
            s.cross_species = s.cross_species_detail = None
        counts[s.cross_species or ("none" if s.cross_species == "" else "not evaluable")] += 1
    db.commit()
    return dict(counts)
