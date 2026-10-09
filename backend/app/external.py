"""
Cached lookups against UniProt and the AlphaFold Protein Structure Database.

Responses are reduced to what the protein views need and cached as JSON under
SAAP_CACHE_DIR (default backend/.cache). Definitive answers — including "no
entry" — are cached; network failures raise ExternalError and are retried on
the next request.

Positional data (features, pLDDT, AlphaMissense) is only meaningful when the
remote sequence is identical to the SAAP's cached protein sequence, so callers
must check `sequence` before mapping positions.

References
----------
Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T.,
    Pritzel, A., Wong, L. H., Zielinski, M., Sargeant, T., Schneider, R. G.,
    Senior, A. W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023).
    Accurate proteome-wide missense variant effect prediction with
    AlphaMissense. Science, 381(6664), Article eadg7492.
    https://doi.org/10.1126/science.adg7492

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

import csv
import io
import json
import os
import re
import threading
from pathlib import Path

from .scoring import AA

CACHE_DIR = Path(os.environ.get(
    "SAAP_CACHE_DIR", Path(__file__).resolve().parent.parent / ".cache"))
UNIPROT = "https://rest.uniprot.org/uniprotkb"
ALPHAFOLD = "https://alphafold.ebi.ac.uk/api/prediction"
TIMEOUT = 30

_ACCESSION = re.compile(r"^[A-Z0-9]{6,10}(-\d+)?$")

# UniProt feature type -> display track. Types not listed are dropped.
FEATURE_TRACKS = {
    "Domain": "domain", "Repeat": "domain", "Zinc finger": "domain",
    "Region": "region", "Motif": "region", "Coiled coil": "region",
    "Compositional bias": "region",
    "Signal": "topology", "Transit peptide": "topology", "Propeptide": "topology",
    "Topological domain": "topology", "Transmembrane": "topology",
    "Intramembrane": "topology",
    "Modified residue": "ptm", "Glycosylation": "ptm", "Lipidation": "ptm",
    "Cross-link": "ptm", "Disulfide bond": "ptm",
    "Active site": "site", "Binding site": "site", "Site": "site",
    "Natural variant": "variant",
}


class ExternalError(Exception):
    """A remote service could not be reached or answered with an error."""


def valid_accession(acc: str) -> bool:
    return bool(_ACCESSION.match(acc or ""))


def _get(url: str, **params):
    import requests

    try:
        resp = requests.get(url, params=params or None, timeout=TIMEOUT)
    except requests.RequestException as exc:
        raise ExternalError(f"{type(exc).__name__}: {exc}") from exc
    if resp.status_code in (400, 404):
        return None
    if not resp.ok:
        raise ExternalError(f"HTTP {resp.status_code} from {url}")
    return resp


def _cached(kind: str, key: str, build):
    path = CACHE_DIR / kind / f"{key}.json"
    if path.exists():
        return json.loads(path.read_text())
    data = build()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")  # concurrent requests
    tmp.write_text(json.dumps(data, separators=(",", ":")))
    tmp.replace(path)
    return data


# --- UniProt -----------------------------------------------------------------
def _features(entry: dict) -> list[dict]:
    out = []
    for f in entry.get("features") or []:
        track = FEATURE_TRACKS.get(f.get("type"))
        loc = f.get("location") or {}
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if not track or start is None or end is None:
            continue
        base = {"type": f["type"], "track": track, "description": f.get("description") or ""}
        if f["type"] == "Disulfide bond":  # two bonded residues, not a span
            out += [{**base, "start": p, "end": p} for p in {start, end}]
            continue
        if track == "variant":
            alt = f.get("alternativeSequence") or {}
            base["ref"] = alt.get("originalSequence") or ""
            base["alt"] = ",".join(alt.get("alternativeSequences") or [])
            base["id"] = f.get("featureId") or ""
        out.append({**base, "start": start, "end": end})
    return out


def uniprot_entry(acc: str) -> dict:
    def build():
        resp = _get(f"{UNIPROT}/{acc}.json")
        if resp is None:
            return {"available": False}
        e = resp.json()
        return {
            "available": True,
            "reviewed": "Swiss-Prot" in (e.get("entryType") or ""),
            "organism": (e.get("organism") or {}).get("scientificName"),
            "sequence": (e.get("sequence") or {}).get("value"),
            "features": _features(e),
        }
    return _cached("uniprot", acc, build)


# --- AlphaFold DB ------------------------------------------------------------
def _alphamissense(url: str) -> dict[str, list]:
    """{position: [score per AA in scoring.AA order, None for the reference]}
    plus {position: [class per AA]} — compact enough to cache per protein."""
    resp = _get(url)
    if resp is None:
        return {}
    table: dict[str, list] = {}
    for row in csv.DictReader(io.StringIO(resp.text)):
        v = row["protein_variant"]
        pos, alt = v[1:-1], v[-1]
        if alt not in AA:
            continue
        scores = table.setdefault(pos, [None] * 20)
        scores[AA.index(alt)] = [round(float(row["am_pathogenicity"]), 4), row["am_class"]]
    return table


def alphafold_entry(acc: str) -> dict:
    def build():
        resp = _get(f"{ALPHAFOLD}/{acc}")
        models = resp.json() if resp is not None else []
        model = next((m for m in models if m.get("uniprotAccession") == acc), None)
        if model is None:
            return {"available": False}
        conf = _get(model["plddtDocUrl"])
        am_url = model.get("amAnnotationsUrl")
        return {
            "available": True,
            "model_id": model.get("modelEntityId"),
            "version": model.get("latestVersion"),
            "sequence": model.get("uniprotSequence") or model.get("sequence"),
            "start": model.get("uniprotStart") or model.get("sequenceStart") or 1,
            "global_plddt": model.get("globalMetricValue"),
            "plddt": conf.json().get("confidenceScore") if conf is not None else None,
            "pdb_url": model.get("pdbUrl"),
            "alphamissense": _alphamissense(am_url) if am_url else None,
        }
    return _cached("alphafold", acc, build)


def alphafold_pdb(acc: str) -> str | None:
    path = CACHE_DIR / "alphafold" / f"{acc}.pdb"
    if path.exists():
        return path.read_text()
    entry = alphafold_entry(acc)
    if not entry.get("available") or not entry.get("pdb_url"):
        return None
    resp = _get(entry["pdb_url"])
    if resp is None:
        return None
    path.write_text(resp.text)
    return resp.text
