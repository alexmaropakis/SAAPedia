"""
Protein-level investigation: SAAP sites on their protein, with UniProt
features, AlphaFold confidence and AlphaMissense scores mapped onto them.

Remote tracks are only attached when the remote sequence is identical to the
cached protein sequence — otherwise positions would silently shift.

References
----------
Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T.,
    Pritzel, A., Wong, L. H., Zielinski, M., Sargeant, T., Schneider, R. G.,
    Senior, A. W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023).
    Accurate proteome-wide missense variant effect prediction with
    AlphaMissense. Science, 381(6664), Article eadg7492.
    https://doi.org/10.1126/science.adg7492

Grantham, R. (1974). Amino acid difference formula to help explain protein
    evolution. Science, 185(4154), 862–864.
    https://doi.org/10.1126/science.185.4154.862

Henikoff, S., & Henikoff, J. G. (1992). Amino acid substitution matrices from
    protein blocks. Proceedings of the National Academy of Sciences, 89(22),
    10915–10919. https://doi.org/10.1073/pnas.89.22.10915

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O.,
    Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A.,
    Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B.,
    Nikolov, S., Jain, R., Adler, J., . . . Hassabis, D. (2021). Highly
    accurate protein structure prediction with AlphaFold. Nature, 596(7873),
    583–589. https://doi.org/10.1038/s41586-021-03819-2

The UniProt Consortium. (2025). UniProt: The Universal Protein Knowledgebase in
    2025. Nucleic Acids Research, 53(D1), D609–D617.
    https://doi.org/10.1093/nar/gkae1010

Varadi, M., Bertoni, D., Magana, P., Paramval, U., Pidruchna, I.,
    Radhakrishnan, M., Tsenkov, M., Nair, S., Mirdita, M., Yeo, J.,
    Kovalevskiy, O., Tunyasuvunakool, K., Laydon, A., Žídek, A., Tomlinson, H.,
    Hariharan, D., Abrahamson, J., Green, T., Jumper, J., . . . Velankar, S.
    (2024). AlphaFold Protein Structure Database in 2024: Providing structure
    coverage for over 214 million protein sequences. Nucleic Acids Research,
    52(D1), D368–D375. https://doi.org/10.1093/nar/gkad1011
"""
from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from . import crud, external
from .annotate import parse_substitution, stored_positions, substitution_offset
from .models import SAAP
from .scoring import AA, substitution_scores


def residues(saap: SAAP) -> tuple[str | None, str | None]:
    """(reference, substituted) residue — from the peptide pair when it
    differs at one position, else from the AAS string."""
    offset = substitution_offset(saap.bp_seq, saap.mtp_seq)
    if offset is not None:
        return saap.bp_seq[offset].upper(), saap.mtp_seq[offset].upper()
    return parse_substitution(saap.aa_sub)


def site_info(saap: SAAP) -> dict:
    ref, alt = residues(saap)
    start = saap.peptide_start
    return {
        "ref": ref, "alt": alt,
        "positions": stored_positions(saap),
        "peptide_start": start,
        "peptide_end": start + len(saap.bp_seq) - 1 if start and saap.bp_seq else None,
        **substitution_scores(ref, alt),
    }


def protein_view(db: Session, accession: str) -> dict | None:
    first = db.scalar(select(SAAP).where(SAAP.protein_accession == accession,
                                         SAAP.protein_sequence.is_not(None)).limit(1))
    if first is None:
        return None
    rows = crud.protein_saaps(db, accession)
    by_id = {s.id: s for s in db.scalars(select(SAAP).where(SAAP.protein_accession == accession))}
    sites = [{**r, **site_info(by_id[r["id"]])} for r in rows]
    return {
        "protein": {
            "accession": accession,
            "gene": (first.source_gene or "").split(";")[0] or None,
            "description": first.protein_description,
            "length": len(first.protein_sequence),
            "sequence": first.protein_sequence,
            "ensembl_gene": first.ensembl_gene,
            "ensembl_protein": first.ensembl_protein,
        },
        "sites": sites,
    }


def _remote(fetch, accession: str) -> dict:
    try:
        return fetch(accession)
    except external.ExternalError as exc:
        return {"available": False, "error": str(exc)}


NEAREST = 3


def _nearest(pos: int, functional: list[dict], coords: dict) -> list[dict]:
    """Closest annotated functional residues to `pos` in the AlphaFold model
    (C-alpha distance, Å); the site's own residue is included at 0 Å."""
    if pos not in coords:
        return []
    x0 = coords[pos]
    by_residue: dict[int, dict] = {}
    for f in functional:
        for r in {f["start"], f["end"]}:
            if r not in coords:
                continue
            hit = by_residue.setdefault(r, {
                "residue": r, "annotations": [], "sequence_separation": abs(r - pos),
                "distance": round(sum((a - b) ** 2 for a, b in zip(x0, coords[r])) ** 0.5, 1)})
            label = f"{f['type']}: {f['description']}" if f["description"] else f["type"]
            if label not in hit["annotations"]:
                hit["annotations"].append(label)
    return sorted(by_residue.values(), key=lambda h: h["distance"])[:NEAREST]


def substitution_matrix(db: Session, species: str | None = None, tissue: str | None = None) -> dict:
    """Distinct SAAP per (reference, substituted) residue, optionally within a species / tissue."""
    from sqlalchemy import and_, exists, or_

    from .models import Observation

    stmt = select(SAAP)
    conds = []
    if species:
        conds.append(Observation.species == species)
    if tissue:
        conds.append(or_(Observation.tissue == tissue, Observation.tissue.like(tissue.replace("%", "") + " (%")))
    if conds:
        stmt = stmt.where(exists().where(and_(Observation.saap_id == SAAP.id, *conds)))
    counts: dict[str, dict[str, int]] = {}
    total = 0
    for s in db.scalars(stmt):
        ref, alt = residues(s)
        if ref in AA and alt in AA and ref != alt:
            counts.setdefault(ref, {}).setdefault(alt, 0)
            counts[ref][alt] += 1
            total += 1
    return {"residues": AA, "counts": counts, "total": total}


def protein_annotations(db: Session, accession: str) -> dict | None:
    view = protein_view(db, accession)
    if view is None:
        return None
    seq = view["protein"]["sequence"]

    up = _remote(external.uniprot_entry, accession)
    up_match = up.get("available") and up.get("sequence") == seq
    features = up.get("features", []) if up_match else []

    af = _remote(external.alphafold_entry, accession)
    af_match = af.get("available") and af.get("sequence") == seq and af.get("start", 1) == 1
    plddt = af.get("plddt") if af_match else None
    am = af.get("alphamissense") if af_match else None
    am_mean = None
    if am:
        am_mean = []
        for pos in range(1, len(seq) + 1):
            vals = [c[0] for c in am.get(str(pos), []) if c]
            am_mean.append(round(sum(vals) / len(vals), 3) if vals else None)

    variants = [f for f in features if f["track"] == "variant"]
    # Functional residues for 3D proximity: active/binding/other sites and PTMs.
    functional = [f for f in features if f["track"] in ("site", "ptm")]
    coords = {}
    if af_match and af.get("pdb_url") and functional:
        try:
            coords = external.ca_coordinates(accession)
        except external.ExternalError:
            coords = {}
    sites = {}
    for s in view["sites"]:
        pos, alt = (s["positions"] or [None])[0], s["alt"]
        if not pos:
            continue
        cell = (am.get(str(pos)) or [None] * 20)[AA.index(alt)] if am and alt in AA else None
        sites[s["id"]] = {
            "nearest": _nearest(pos, functional, coords),
            "plddt": plddt[pos - 1] if plddt and pos <= len(plddt) else None,
            "am_pathogenicity": cell[0] if cell else None,
            "am_class": cell[1] if cell else None,
            "known_variants": [v for v in variants if v["start"] <= pos <= v["end"]],
        }

    def status(remote, matched):
        if remote.get("error"):
            return remote["error"]
        if not remote.get("available"):
            return "no entry"
        return None if matched else "sequence differs from annotated protein"

    return {
        "uniprot": {
            "reviewed": up.get("reviewed"), "organism": up.get("organism"),
            "features": features, "issue": status(up, up_match),
        },
        "alphafold": {
            "model_id": af.get("model_id"), "version": af.get("version"),
            "global_plddt": af.get("global_plddt"),
            "plddt": plddt, "am_mean": am_mean,
            "has_structure": bool(af_match and af.get("pdb_url")),
            "issue": status(af, af_match),
        },
        "sites": sites,
    }
