"""
Public sample labels: every observation is described only by its tissue or
cell type and its species. Dataset/study names, TMT plexes and source files
stay private (see PRIVATE); datasets of the same tissue are grouped, and
disease cohorts (AD, PD, cancers, ...) are labelled by tissue alone. Labels
read "Organ (sub-site)" where a sub-site is known, e.g. "Brain (frontal cortex)".
"""
from __future__ import annotations

import os
import re

# Expose dataset names, plexes and files (local use only — never on a public site).
PRIVATE = os.environ.get("SAAP_PRIVATE", "").lower() in {"1", "true", "yes"}

# Labels read "Organ (sub-site)" where the data records a sub-site; filtering on
# the organ ("Brain") includes every sub-site (see crud._apply_filters).
CELL_TYPES = {"bcells": "Blood (B cells)", "monocytes": "Blood (monocytes)", "nkcells": "Blood (NK cells)",
              "tcells": "Blood (T cells)", "neutrophils": "Blood (neutrophils)", "hepatocytes": "Liver (hepatocytes)"}

# Dataset (or per-row tissue) name, normalized -> public tissue.
TISSUES = {
    # brain: regions, fractions and brain-bank cohorts
    "anteriorcingulategyrus": "Brain (anterior cingulate gyrus)", "frontalcortex": "Brain (frontal cortex)",
    "cerebralcortex": "Brain (cerebral cortex)", "cerebellum": "Brain (cerebellum)", "cortex": "Brain (cortex)",
    "hippocampus": "Brain (hippocampus)", "striatum": "Brain (striatum)",
    "insolubleproteomematched": "Brain (insoluble)",
    "pooledbrain": "Brain", "brain": "Brain", "ptipd2026": "Brain", "ptiad2026": "Brain",
    # tumours -> tissue of origin
    "brca": "Breast", "ccrcc": "Kidney", "luad": "Lung", "lscc": "Lung", "pdac": "Pancreas", "ucec": "Endometrium",
    # tissues
    "lung": "Lung", "kidney": "Kidney", "pancreas": "Pancreas", "endometrium": "Endometrium",
    "muscle": "Skeletal muscle", "plasma": "Blood (plasma)", "liver": "Liver", "heart": "Heart",
    "fat": "Adipose tissue", "adiposetissue": "Adipose tissue", "aorta": "Aorta", "skin": "Skin",
    "spleen": "Spleen", "placenta": "Placenta", "stomach": "Stomach", "gallbladder": "Gallbladder",
    "thyroidgland": "Thyroid gland", "salivasecretinggland": "Salivary gland",
    "fallopiantube": "Fallopian tube", "esophagus": "Esophagus", "urinarybladder": "Urinary bladder",
    "ovary": "Ovary", "prostategland": "Prostate", "testis": "Testis", "lymphnode": "Lymph node",
    "bonemarrow": "Bone marrow",
    "duodenum": "Small intestine (duodenum)", "smallintestine": "Small intestine",
    "colon": "Large intestine (colon)", "vermiformappendix": "Large intestine (appendix)",
}
UNSPECIFIED = "Unspecified tissue"

_norm = lambda s: re.sub(r"[^a-z0-9]", "", (s or "").lower())


def organ(label: str) -> str:
    return label.split(" (")[0]


def sample_label(dataset: str | None, tmt_tissue: str | None) -> tuple[str, str]:
    """(tissue or cell type, "tissue" | "cell type") for one observation."""
    for name in (tmt_tissue, dataset):  # a per-row tissue (pan-tissue sets) wins over the dataset
        key = _norm(name)
        if key in CELL_TYPES:
            return CELL_TYPES[key], "cell type"
        if key in TISSUES:
            return TISSUES[key], "tissue"
    return UNSPECIFIED, "tissue"


def backfill(db) -> int:
    """Label every observation that has no public tissue yet."""
    from sqlalchemy import select

    from .models import Observation

    n = 0
    for o in db.scalars(select(Observation).where(Observation.tissue.is_(None))):
        o.tissue, o.sample_type = sample_label(o.dataset, o.tmt_tissue)
        n += 1
    db.commit()
    return n
