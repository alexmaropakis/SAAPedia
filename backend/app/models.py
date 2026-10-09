"""
ORM models

A SAAP (single amino-acid polymorphism / substituted peptide) is the de-dup
grain: one row per unique (mtp_seq, bp_seq, aa_sub). Every line in an imported
file becomes an Observation linked to its SAAP, carrying the per-dataset metrics.
"N datasets" is not stored from the file — it is computed as the number of
distinct datasets the SAAP appears in.
"""

from __future__ import annotations
from typing import Optional
from sqlalchemy import Boolean, Float, ForeignKey, Integer, String, Text, UniqueConstraint
from sqlalchemy.orm import Mapped, mapped_column, relationship
from .database import Base


class SAAP(Base):
    __tablename__ = "saap"
    __table_args__ = (
        UniqueConstraint("mtp_seq", "bp_seq", "aa_sub", name="uq_saap_identity"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True)

    # Identity (de-dup key)
    mtp_seq: Mapped[str] = mapped_column(String, index=True)       # SAAP (variant peptide)
    bp_seq: Mapped[str] = mapped_column(String, index=True)        # BP (base peptide)
    aa_sub: Mapped[str] = mapped_column(String, index=True)        # AAS, e.g. "M to D"

    # Per-peptide attributes (constant for the peptide; back-filled on ingest)
    source_accession: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)  # UniProt
    source_gene: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)        # Genes
    ref_proteins: Mapped[Optional[str]] = mapped_column(String, nullable=True)                   # RefProteins

    # --- Ensembl / positional annotation -------------------------------------
    # Populated either from columns in the imported file or by the annotation
    # step (see annotate.py), which resolves them from the UniProt accession.
    ensembl_gene: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)      # ENSG...
    ensembl_transcript: Mapped[Optional[str]] = mapped_column(String, nullable=True)            # ENST...
    ensembl_protein: Mapped[Optional[str]] = mapped_column(String, nullable=True)               # ENSP...
    protein_description: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    protein_length: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    # Canonical protein sequence from UniProt, cached at annotation time. Needed
    # to emit full-length protein entries (with the substitution applied in
    # place) rather than bare peptides — see fasta.py `entry_mode="protein"`.
    protein_sequence: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    # The UniProt accession (isoform-specific when applicable) that
    # protein_sequence was taken from — the key SAAPs are grouped by in the
    # protein views, and what AlphaFold / UniProt features are looked up with.
    protein_accession: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)
    # 1-based index of the substituted residue within the full protein, and the
    # 1-based offset of the peptide's first residue. Both derived by locating
    # bp_seq in the canonical protein sequence. When the peptide occurs more
    # than once these hold the FIRST occurrence, and the *_all columns list
    # every occurrence as a comma-separated string.
    position_in_protein: Mapped[Optional[int]] = mapped_column(Integer, index=True, nullable=True)
    peptide_start: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    positions_all: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    peptide_starts_all: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    # How many times bp_seq occurs in the protein (1 for the normal case).
    n_positions: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)
    # Reference-proteome accessions whose sequence already contains the
    # substituted peptide (I = L), comma-separated; "" = none, NULL = unchecked.
    # A hit means the "substitution" is genome-encoded elsewhere (see proteome.py).
    proteome_hits: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    # gnomAD v4 population check (see gnomad.py): "absent" / "present" /
    # "unmapped" (no Ensembl transcript identical to the protein); NULL =
    # not checked (e.g. non-human). gnomad_af sums every variant producing
    # this exact substitution; gnomad_variants lists their gnomAD IDs.
    # UniProt natural variants producing this exact substitution ("VAR_xxx,..."),
    # "" = none, NULL = not checked (no UniProt entry matching the sequence).
    known_variant: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    gnomad_status: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    gnomad_af: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    gnomad_variants: Mapped[Optional[str]] = mapped_column(Text, nullable=True)
    # Best positional probability for this SAAP in the original source tables
    # (Project_AAS processed data, Tsour et al.), for observations imported
    # without one. Combined with the observations' own values in rollups.
    source_positional_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    # Provenance for the annotation: e.g. "file", "uniprot", or "uniprot:unmatched".
    annotation_source: Mapped[Optional[str]] = mapped_column(String, nullable=True)

    immunoglobulin: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    trypsin: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    missed_cleavage: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    aas_at_peptide_terminus: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)
    greater_than_shared: Mapped[Optional[bool]] = mapped_column(Boolean, nullable=True)

    observations: Mapped[list["Observation"]] = relationship(
        back_populates="saap", cascade="all, delete-orphan"
    )


class Observation(Base):
    __tablename__ = "observation"
    __table_args__ = (
        UniqueConstraint("row_hash", name="uq_observation_row_hash"),
    )

    id: Mapped[int] = mapped_column(Integer, primary_key=True)
    saap_id: Mapped[int] = mapped_column(ForeignKey("saap.id"), index=True)

    dataset: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)
    tmt_tissue: Mapped[Optional[str]] = mapped_column(String, nullable=True)          # TMT/Tissue
    digest: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)
    species: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)
    acquisition_type: Mapped[Optional[str]] = mapped_column(String, index=True, nullable=True)  # Data acquisition

    saap_pep: Mapped[Optional[float]] = mapped_column(Float, nullable=True)            # PEP
    positional_probability: Mapped[Optional[float]] = mapped_column(Float, nullable=True)
    n_evidence_fragments: Mapped[Optional[int]] = mapped_column(Integer, nullable=True)

    source_file: Mapped[Optional[str]] = mapped_column(String, nullable=True)
    # Hash of identity + observation fields -> guards against duplicate rows on re-import.
    row_hash: Mapped[str] = mapped_column(String, index=True)

    saap: Mapped["SAAP"] = relationship(back_populates="observations")
