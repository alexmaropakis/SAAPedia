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
Ashburner, M., Ball, C. A., Blake, J. A., Botstein, D., Butler, H., Cherry, J.
    M., Davis, A. P., Dolinski, K., Dwight, S. S., Eppig, J. T., Harris, M. A.,
    Hill, D. P., Issel-Tarver, L., Kasarskis, A., Lewis, S., Matese, J. C.,
    Richardson, J. E., Ringwald, M., Rubin, G. M., & Sherlock, G. (2000). Gene
    Ontology: Tool for the unification of biology. Nature Genetics, 25(1),
    25–29. https://doi.org/10.1038/75556

Cheng, J., Novati, G., Pan, J., Bycroft, C., Žemgulytė, A., Applebaum, T.,
    Pritzel, A., Wong, L. H., Zielinski, M., Sargeant, T., Schneider, R. G.,
    Senior, A. W., Jumper, J., Hassabis, D., Kohli, P., & Avsec, Ž. (2023).
    Accurate proteome-wide missense variant effect prediction with
    AlphaMissense. Science, 381(6664), Article eadg7492.
    https://doi.org/10.1126/science.adg7492

The Gene Ontology Consortium, Aleksander, S. A., Balhoff, J., Carbon, S.,
    Cherry, J. M., Drabkin, H. J., Ebert, D., Feuermann, M., Gaudet, P.,
    Harris, N. L., Hill, D. P., Lee, R., Mi, H., Moxon, S., Mungall, C. J.,
    Muruganugan, A., Mushayahama, T., Sternberg, P. W., Thomas, P. D., . . .
    Westerfield, M. (2023). The Gene Ontology knowledgebase in 2023. Genetics,
    224(1), Article iyad031. https://doi.org/10.1093/genetics/iyad031

Jumper, J., Evans, R., Pritzel, A., Green, T., Figurnov, M., Ronneberger, O.,
    Tunyasuvunakool, K., Bates, R., Žídek, A., Potapenko, A., Bridgland, A.,
    Meyer, C., Kohl, S. A. A., Ballard, A. J., Cowie, A., Romera-Paredes, B.,
    Nikolov, S., Jain, R., Adler, J., . . . Hassabis, D. (2021). Highly
    accurate protein structure prediction with AlphaFold. Nature, 596(7873),
    583–589. https://doi.org/10.1038/s41586-021-03819-2

Milacic, M., Beavers, D., Conley, P., Gong, C., Gillespie, M., Griss, J., Haw,
    R., Jassal, B., Matthews, L., May, B., Petryszak, R., Ragueneau, E.,
    Rothfels, K., Sevilla, C., Shamovsky, V., Stephan, R., Tiwari, K., Varusai,
    T., Weiser, J., . . . D’Eustachio, P. (2024). The Reactome Pathway
    Knowledgebase 2024. Nucleic Acids Research, 52(D1), D672–D678.
    https://doi.org/10.1093/nar/gkad1025

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


def _cached(kind: str, key: str, build, version: int | None = None):
    """Cached JSON; an entry written by an older `version` of `build` is rebuilt."""
    path = CACHE_DIR / kind / f"{key}.json"
    if path.exists():
        data = json.loads(path.read_text())
        if version is None or data.get("v") == version:
            return data
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
        ligand = (f.get("ligand") or {}).get("name")
        base = {"type": f["type"], "track": track, "description": f.get("description") or ligand or ""}
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


UNIPROT_ENTRY_VERSION = 2  # 2: binding sites carry their ligand name


def uniprot_entry(acc: str) -> dict:
    def build():
        resp = _get(f"{UNIPROT}/{acc}.json")
        if resp is None:
            return {"v": UNIPROT_ENTRY_VERSION, "available": False}
        e = resp.json()
        return {
            "v": UNIPROT_ENTRY_VERSION,
            "available": True,
            "reviewed": "Swiss-Prot" in (e.get("entryType") or ""),
            "organism": (e.get("organism") or {}).get("scientificName"),
            "sequence": (e.get("sequence") or {}).get("value"),
            "features": _features(e),
        }
    return _cached("uniprot", acc, build, version=UNIPROT_ENTRY_VERSION)


def uniprot_function(acc: str) -> dict:
    """Function summary, subcellular locations, GO terms, Reactome pathways and
    keywords. Entry-level annotations, so isoforms use their canonical entry."""
    base = acc.split("-")[0]

    def build():
        resp = _get(f"{UNIPROT}/{base}", fields="cc_function,cc_subcellular_location,go,cc_pathway,xref_reactome,keyword",
                    format="json")
        if resp is None:
            return {"available": False}
        e = resp.json()
        out = {"available": True, "function": [], "locations": [], "pathways": [], "keywords": [],
               "go": {"P": [], "F": [], "C": []}, "reactome": []}
        for c in e.get("comments", []):
            if c.get("commentType") == "FUNCTION":
                out["function"] += [t["value"] for t in c.get("texts", [])]
            elif c.get("commentType") == "SUBCELLULAR LOCATION":
                out["locations"] += [l["location"]["value"] for l in c.get("subcellularLocations", []) if l.get("location")]
            elif c.get("commentType") == "PATHWAY":
                out["pathways"] += [t["value"] for t in c.get("texts", [])]
        for x in e.get("uniProtKBCrossReferences", []):
            props = {p["key"]: p["value"] for p in x.get("properties", [])}
            if x["database"] == "GO" and ":" in props.get("GoTerm", ""):
                aspect, term = props["GoTerm"].split(":", 1)
                out["go"].setdefault(aspect, []).append({"id": x["id"], "term": term})
            elif x["database"] == "Reactome":
                out["reactome"].append({"id": x["id"], "name": props.get("PathwayName", x["id"])})
        out["keywords"] = [k["name"] for k in e.get("keywords", [])]
        out["locations"] = list(dict.fromkeys(out["locations"]))
        return out
    return _cached("uniprot_function", base, build)


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


def ca_coordinates(acc: str) -> dict[int, tuple[float, float, float]]:
    """C-alpha coordinates by residue number from the AlphaFold model ({} if none)."""
    pdb = alphafold_pdb(acc)
    out = {}
    for line in (pdb or "").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            out[int(line[22:26])] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return out


def alphafold_pdb(acc: str) -> str | None:
    path = CACHE_DIR / "alphafold" / f"{acc}.pdb"
    if path.exists():
        return path.read_text()
    path.parent.mkdir(parents=True, exist_ok=True)
    entry = alphafold_entry(acc)
    if not entry.get("available") or not entry.get("pdb_url"):
        return None
    resp = _get(entry["pdb_url"])
    if resp is None:
        return None
    path.write_text(resp.text)
    return resp.text
