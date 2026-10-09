"""
Dataset curation, applied after every import and annotation run.

  * Immunoglobulins are removed — flagged by the file's Immunoglobulin column
    or recognized from the gene symbol, the annotated protein name, or the
    file's RefProteins name.
  * Contaminants are removed — the base or substituted peptide occurs in a
    common-contaminant protein (data/contaminants.fasta: protease reagents,
    BSA, casein, skin/hair/wool keratins, ...), or the file flags the row as
    Trypsin.
  * Genome-encoded SAAPs are removed — the substituted peptide occurs
    verbatim in the reference proteome (proteome.py), so no substitution is
    needed to explain it.
  * SAAPs reproducible by a gnomAD germline variant at allele frequency
    >= gnomad.MAX_POPULATION_AF are removed as likely polymorphisms.
  * SAAPs matching a known genetic variant (UniProt natural variant with the
    same residue change) are kept only when gnomAD confirms them rare
    (present below MAX_POPULATION_AF, or absent); unconfirmed ones are removed.
  * Positional probability is only logged, never used to remove anything:
    the counts below MIN_POSITIONAL_PROBABILITY (or with none at all) are
    reported, and the Browse "Min PosProb" filter applies the threshold.

References
----------
The Global Proteome Machine Organization. (n.d.). cRAP protein sequences [Data
    set]. Retrieved October 8, 2026, from https://www.thegpm.org/crap/

Guez, J., Goodrich, J. K., Moldovan, M. A., Chao, K. R., Kar, P., Panchal, R.,
    Wilson, M. W., Laricchia, K. M., Rohlicek, G., Biba, D., Marten, D., He,
    Q., Darnowsky, P. W., Grant, R., Weisburd, B., Baxter, S. M., Nadeau, J.,
    Lu, W., Jahl, S., . . . Karczewski, K. J. (2026). Integrating 730,947 exome
    sequences with clinical literature improves gene discovery [Preprint].
    medRxiv. https://doi.org/10.64898/2026.03.23.26349081
"""
from __future__ import annotations

import logging
from collections import Counter
from functools import lru_cache
from pathlib import Path

from sqlalchemy import delete, func, select
from sqlalchemy.orm import Session

from .models import Observation, SAAP
from .util import is_immunoglobulin

MIN_POSITIONAL_PROBABILITY = 0.9
log = logging.getLogger("saapedia.curate")


def best_positional_probability(agg_max):
    """Best of the observations' max and the source-data value (NULL-safe)."""
    src = SAAP.source_positional_probability
    return func.max(func.coalesce(agg_max, src), func.coalesce(src, agg_max))


CONTAMINANTS_FASTA = Path(__file__).resolve().parent / "data" / "contaminants.fasta"


@lru_cache(maxsize=1)
def contaminant_sequences() -> dict[str, str]:
    """cRAP common contaminants (protease reagents, BSA, casein, keratins from
    skin/hair/wool, ...) without the UPS human standards, plus Arg-C; I->L."""
    out, name = {}, None
    for line in CONTAMINANTS_FASTA.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            out[name] = ""
        elif name and not line.startswith("#"):
            out[name] += line.strip().upper().replace("I", "L")
    return out


_SPECIES_SUFFIX = {"homo sapiens": "_HUMAN", "mus musculus": "_MOUSE"}


def contaminant_hits(db: Session) -> dict[int, str]:
    """SAAP id -> contaminant explaining it, plus rows the file flags as Trypsin.

    A peptide shared with a contaminant counts only when it cannot be the
    sample's own protein: a contaminant from the sample's species (human
    keratins, amylase in human samples) always counts; one from another
    species (bovine, rabbit, sheep, ...) only if the base peptide is absent
    from the sample's reference proteome — conserved peptides such as
    aldolase or cytochrome c are endogenous, not contamination.
    """
    from .proteome import references

    seqs, refs = contaminant_sequences(), references()
    species: dict[int, set[str]] = {}
    for sid, sp in db.execute(select(Observation.saap_id, Observation.species).distinct()):
        if sp:
            species.setdefault(sid, set()).add(sp.strip().lower())
    hits = {sid: "Trypsin (file flag)" for sid, flag in db.execute(select(SAAP.id, SAAP.trypsin)) if flag}
    for sid, bp, mtp in db.execute(select(SAAP.id, SAAP.bp_seq, SAAP.mtp_seq)):
        peps = [p.upper().replace("I", "L") for p in (bp, mtp) if p]
        for name, seq in seqs.items():
            if not any(p in seq for p in peps):
                continue
            sample = species.get(sid, set())
            own = any(name.endswith(_SPECIES_SUFFIX.get(sp, "?")) for sp in sample)
            endogenous = bp and any(sp in refs and refs[sp].find(bp) for sp in sample)
            if own or not endogenous:
                hits[sid] = name
                break
    return hits


def mark_known_variants(db: Session) -> int:
    """Record UniProt natural variants that reproduce each SAAP's substitution."""
    import json

    from .annotate import stored_positions
    from .external import CACHE_DIR
    from .investigate import residues

    cache: dict[str, list | None] = {}
    n = 0
    for s in db.scalars(select(SAAP).where(SAAP.protein_accession.is_not(None))):
        acc = s.protein_accession
        if acc not in cache:
            path = CACHE_DIR / "uniprot" / f"{acc}.json"
            entry = json.loads(path.read_text()) if path.exists() else {}
            cache[acc] = ([f for f in entry.get("features", []) if f["track"] == "variant"]
                          if entry.get("sequence") == s.protein_sequence else None)
        if cache[acc] is None:
            continue
        _, alt = residues(s)
        positions = set(stored_positions(s))
        hits = [f.get("id") or f"{f['ref']}{f['start']}{alt}" for f in cache[acc]
                if f["start"] in positions and alt and alt in (f.get("alt") or "").split(",")]
        s.known_variant = ",".join(hits)
        n += bool(hits)
    db.commit()
    return n


def prune(db: Session) -> dict:
    from .gnomad import MAX_POPULATION_AF

    mark_known_variants(db)

    ig = {sid for sid, flag, gene, desc, ref in db.execute(
        select(SAAP.id, SAAP.immunoglobulin, SAAP.source_gene, SAAP.protein_description, SAAP.ref_proteins))
        if flag or is_immunoglobulin(gene, desc) or is_immunoglobulin(None, ref)}
    common = set(db.scalars(select(SAAP.id).where(SAAP.gnomad_af >= MAX_POPULATION_AF)))
    # Known genetic variants must pass the AF filter explicitly to stay.
    rare_confirmed = (SAAP.gnomad_status == "absent") | (
        (SAAP.gnomad_status == "present") & (SAAP.gnomad_af < MAX_POPULATION_AF))
    known_unconfirmed = set(db.scalars(select(SAAP.id).where(
        SAAP.known_variant != "", SAAP.known_variant.is_not(None), ~rare_confirmed)))
    contam = contaminant_hits(db)
    log.info("contaminants: %s", Counter(contam.values()).most_common())
    contaminants = set(contam) - ig
    encoded = set(db.scalars(select(SAAP.id).where(SAAP.proteome_hits != ""))) - ig - contaminants
    known_unconfirmed -= ig | contaminants | encoded | common
    victims = sorted(ig | contaminants | common | encoded | known_unconfirmed)
    for i in range(0, len(victims), 500):
        chunk = victims[i:i + 500]
        db.execute(delete(Observation).where(Observation.saap_id.in_(chunk)))
        db.execute(delete(SAAP).where(SAAP.id.in_(chunk)))
    db.commit()

    result = {"removed_immunoglobulin": len(ig),
              "removed_contaminant": len(contaminants),
              "removed_genome_encoded": len(encoded),
              "removed_population_variant": len(common - ig - contaminants - encoded),
              "removed_known_variant_unconfirmed": len(known_unconfirmed), **summary(db)}
    log.info("curation: %s", result)
    return result


def summary(db: Session) -> dict:
    """Read-only positional-probability log."""
    obs_max = (select(Observation.saap_id, func.max(Observation.positional_probability).label("pp"))
               .group_by(Observation.saap_id).subquery())
    pp = best_positional_probability(obs_max.c.pp)
    base = select(func.count(SAAP.id)).outerjoin(obs_max, obs_max.c.saap_id == SAAP.id)
    return {
        "n_saap": db.scalar(select(func.count(SAAP.id))) or 0,
        "below_min_positional_probability": db.scalar(base.where(pp < MIN_POSITIONAL_PROBABILITY)) or 0,
        "no_positional_probability": db.scalar(base.where(pp.is_(None))) or 0,
        "min_positional_probability": MIN_POSITIONAL_PROBABILITY,
        "gnomad": dict(db.execute(select(SAAP.gnomad_status, func.count(SAAP.id)).group_by(SAAP.gnomad_status))
                       .all()),
    }
