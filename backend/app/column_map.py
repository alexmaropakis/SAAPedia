"""Tolerant file-header -> canonical field mapping.

Headers are normalized (lowercased, non-alphanumerics stripped) before matching,
so "SAAP", "AAS", "TMT/Tissue", "Positional Probability", "UniProt" etc. map
regardless of case/spacing/punctuation. Both the current and older header
spellings are accepted. Add new aliases to ALIASES as source formats vary.
"""
from __future__ import annotations

import re

# canonical field -> set of accepted normalized aliases
ALIASES: dict[str, set[str]] = {
    # --- SAAP identity ---
    "mtp_seq": {"saap", "mtpseq", "mtp", "saapseq", "variantpeptide", "substitutedpeptide", "modseq"},
    "bp_seq": {"bp", "bpseq", "basepeptide", "referencepeptide"},
    "aa_sub": {"aas", "aasub", "aasubstitution", "substitution", "sub", "aachange"},
    # --- per-peptide attributes ---
    "source_accession": {"uniprot", "uniprotaccession", "accession", "sourceaccession", "uniprotid", "proteinaccession"},
    "source_gene": {"genes", "gene", "genename", "sourcegene", "genesymbol"},
    "ref_proteins": {"refproteins", "refprotein", "refprot", "referenceproteins", "referenceprotein", "sourceprotein", "protein"},
    "immunoglobulin": {"immunoglobulin", "ig", "isig"},
    "trypsin": {"trypsin", "istrypsin"},
    "missed_cleavage": {"missedcleavage", "missedcleavages", "misscleavage"},
    "aas_at_peptide_terminus": {"aasatpeptideterminus", "aasatterminus", "aasterminus", "peptideterminus", "aaatterminus"},
    "greater_than_shared": {"greaterthanshared", "gtshared", "greaterthanshareds"},
    # --- per-observation metrics ---
    "dataset": {"dataset", "tissue", "celltype", "tissuecelltype"},
    "tmt_tissue": {"tmttissue", "tmtset", "tmt", "set", "tissuetmt"},
    "digest": {"digest", "protease", "enzyme", "digestion", "digestenzyme"},
    "species": {"species", "organism"},
    "acquisition_type": {"dataacquisitiontype", "dataacquisition", "acquisitiontype", "acquisition", "acqtype", "acquisitionmode"},
    "saap_pep": {"pep", "saappep", "mtppep", "posteriorerrorprobability"},
    "positional_probability": {"positionalprobability", "positionalprob", "posprob", "localizationprobability", "locprob"},
    "n_evidence_fragments": {"nevidencefragments", "evidencefragments", "numevidencefragments", "nfragments", "nevidence"},
}

# Headers we intentionally ignore.
#  - pandas index / row-number columns
#  - "N datasets": recomputed during de-dup, never taken from the file
IGNORE: set[str] = {"rownumber", "index", "unnamed0", "",
                    "ndatasets", "numdatasets", "ndataset", "datasetscount"}

# Fields stored on the SAAP (identity + per-peptide attributes).
SAAP_FIELDS = {
    "mtp_seq", "bp_seq", "aa_sub",
    "source_accession", "source_gene", "ref_proteins",
    "immunoglobulin", "trypsin", "missed_cleavage",
    "aas_at_peptide_terminus", "greater_than_shared",
}
# Everything else canonical is an observation-level field.
OBSERVATION_FIELDS = set(ALIASES) - SAAP_FIELDS

_NORMALIZE_RE = re.compile(r"[^a-z0-9]")


def normalize(header: str) -> str:
    return _NORMALIZE_RE.sub("", (header or "").strip().lower())


# Build reverse lookup: normalized alias -> canonical field.
_REVERSE: dict[str, str] = {}
for _canonical, _aliases in ALIASES.items():
    _REVERSE[normalize(_canonical)] = _canonical
    for _alias in _aliases:
        _REVERSE[normalize(_alias)] = _canonical


def map_headers(headers: list[str]) -> tuple[dict[str, str], list[str]]:
    """Map raw headers to canonical field names.

    Returns (mapping, unmapped) where `mapping` is {raw_header: canonical_field}
    and `unmapped` lists headers we could not place (and are not ignored).
    """
    mapping: dict[str, str] = {}
    unmapped: list[str] = []
    for raw in headers:
        norm = normalize(raw)
        if norm in IGNORE:
            continue
        canonical = _REVERSE.get(norm)
        if canonical is None:
            unmapped.append(raw)
        else:
            mapping[raw] = canonical
    return mapping, unmapped
