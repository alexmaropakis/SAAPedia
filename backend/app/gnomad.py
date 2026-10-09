"""
Population-polymorphism check against gnomAD v4 (exomes + genomes, GRCh38).

A SAAP whose substitution is also a germline missense variant may simply be a
genetic polymorphism carried by the sample. For each human SAAP, gnomAD's
variants are read for an Ensembl transcript whose protein is *identical* to the
SAAP's annotated protein sequence (so positions map exactly), and any variant
producing the same amino-acid change at the same position (HGVSp) is recorded.

The stored frequency is the summed allele frequency of all such variants
(AC / AN over exomes + genomes). Curation removes SAAPs at or above
MAX_POPULATION_AF; rarer matches are kept and annotated.

References
----------
Chen, S., Francioli, L. C., Goodrich, J. K., Collins, R. L., Kanai, M., Wang,
    Q., Alföldi, J., Watts, N. A., Vittal, C., Gauthier, L. D., Poterba, T.,
    Wilson, M. W., Tarasova, Y., Phu, W., Grant, R., Yohannes, M. T., Koenig,
    Z., Farjoun, Y., Banks, E., . . . Karczewski, K. J. (2024). A genomic
    mutational constraint map using variation in 76,156 human genomes. Nature,
    625(7993), 92–100. https://doi.org/10.1038/s41586-023-06045-0

Dyer, S. C., Austine-Orimoloye, O., Azov, A. G., Barba, M., Barnes, I.,
    Barrera-Enriquez, V. P., Becker, A., Bennett, R., Beracochea, M., Berry,
    A., Bhai, J., Bhurji, S. K., Boddu, S., Branco Lins, P. R., Brooks, L.,
    Ramaraju, S. B., Campbell, L. I., Martinez, M. C., Charkhchi, M., . . .
    Yates, A. D. (2025). Ensembl 2025. Nucleic Acids Research, 53(D1),
    D948–D957. https://doi.org/10.1093/nar/gkae1071

Guez, J., Goodrich, J. K., Moldovan, M. A., Chao, K. R., Kar, P., Panchal, R.,
    Wilson, M. W., Laricchia, K. M., Rohlicek, G., Biba, D., Marten, D., He,
    Q., Darnowsky, P. W., Grant, R., Weisburd, B., Baxter, S. M., Nadeau, J.,
    Lu, W., Jahl, S., . . . Karczewski, K. J. (2026). Integrating 730,947 exome
    sequences with clinical literature improves gene discovery [Preprint].
    medRxiv. https://doi.org/10.64898/2026.03.23.26349081

Karczewski, K. J., Francioli, L. C., Tiao, G., Cummings, B. B., Alföldi, J.,
    Wang, Q., Collins, R. L., Laricchia, K. M., Ganna, A., Birnbaum, D. P.,
    Gauthier, L. D., Brand, H., Solomonson, M., Watts, N. A., Rhodes, D.,
    Singer-Berk, M., England, E. M., Seaby, E. G., Kosmicki, J. A., . . .
    MacArthur, D. G. (2020). The mutational constraint spectrum quantified from
    variation in 141,456 humans. Nature, 581(7809), 434–443.
    https://doi.org/10.1038/s41586-020-2308-7
"""
from __future__ import annotations

import logging
import re
import time
from concurrent.futures import ThreadPoolExecutor

from sqlalchemy import select
from sqlalchemy.orm import Session

from . import external
from .annotate import stored_positions
from .investigate import residues
from .models import Observation, SAAP

log = logging.getLogger("saapedia.gnomad")

GNOMAD_API = "https://gnomad.broadinstitute.org/api"
ENSEMBL_API = "https://rest.ensembl.org"
DATASET = "gnomad_r4"
MAX_POPULATION_AF = 1e-4
WORKERS = 4

_HGVSP = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")
_AA1 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
        "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
        "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}


def _post(url: str, payload: dict, tries: int = 5):
    import requests

    for attempt in range(tries):
        try:
            resp = requests.post(url, json=payload, timeout=120,
                                 headers={"Content-Type": "application/json", "Accept": "application/json"})
        except requests.RequestException as exc:
            err = exc
        else:
            if resp.status_code not in (429, 500, 502, 503, 504):
                if not resp.ok:
                    raise external.ExternalError(f"HTTP {resp.status_code} from {url}")
                return resp.json()
            err = f"HTTP {resp.status_code}"
        time.sleep(min(60, 5 * 2 ** attempt))
    raise external.ExternalError(f"{url}: {err}")


def _graphql(query: str) -> dict:
    body = _post(GNOMAD_API, {"query": query})
    if body.get("errors") and not body.get("data"):
        raise external.ExternalError(body["errors"][0].get("message", "gnomAD error"))
    return body.get("data") or {}


def gene_transcripts(gene_id: str | None, symbol: str | None) -> list[str]:
    """Transcripts of a gene, canonical / MANE Select first."""
    if not gene_id and not symbol:
        return []
    key = gene_id or f"symbol_{symbol}"
    arg = f'gene_id: "{gene_id}"' if gene_id else f'gene_symbol: "{symbol}"'

    def build():
        data = _graphql(f"{{ gene({arg}, reference_genome: GRCh38) {{ canonical_transcript_id "
                        f"mane_select_transcript {{ ensembl_id }} transcripts {{ transcript_id }} }} }}")
        g = data.get("gene")
        if not g:
            return {"transcripts": []}
        first = [g.get("canonical_transcript_id"), (g.get("mane_select_transcript") or {}).get("ensembl_id")]
        rest = [t["transcript_id"] for t in g.get("transcripts") or []]
        return {"transcripts": list(dict.fromkeys(t.split(".")[0] for t in first + rest if t))}
    return external._cached("gnomad_gene", key, build)["transcripts"]


def protein_sequence(transcript: str) -> str:
    def build():
        try:
            data = _post(f"{ENSEMBL_API}/sequence/id?type=protein", {"ids": [transcript]})
        except external.ExternalError as exc:
            if "HTTP 400" in str(exc) or "HTTP 404" in str(exc):  # non-coding or retired transcript
                return {"seq": ""}
            raise
        return {"seq": (data[0].get("seq") if data else "") or ""}
    return external._cached("ensembl_protein", transcript, build)["seq"]


def prefetch_proteins(transcripts: list[str]) -> None:
    """Batch-cache Ensembl protein sequences (one request per 10 transcripts)."""
    todo = [t for t in dict.fromkeys(transcripts)
            if t and not (external.CACHE_DIR / "ensembl_protein" / f"{t}.json").exists()]
    for i in range(0, len(todo), 10):  # large batches time out when Ensembl is slow
        chunk = todo[i:i + 10]
        try:
            data = _post(f"{ENSEMBL_API}/sequence/id?type=protein", {"ids": chunk})
        except external.ExternalError:
            # One retired ID fails the whole batch: resolve this chunk one by one.
            for t in chunk:
                try:
                    protein_sequence(t)
                except external.ExternalError:
                    pass
            continue
        seqs = {d.get("query"): d.get("seq") or "" for d in data}
        for t in chunk:
            if t in seqs:
                external._cached("ensembl_protein", t, lambda t=t: {"seq": seqs[t]})


_TX_FIELDS = f"variants(dataset: {DATASET}) {{ variant_id hgvsp exome {{ ac an }} genome {{ ac an }} }}"


def _missense_table(variants: list[dict]) -> dict[str, dict[str, list]]:
    table: dict[str, dict[str, list]] = {}
    for v in variants or []:
        m = _HGVSP.match(v.get("hgvsp") or "")
        if not m or m.group(3) not in _AA1 or m.group(1) == m.group(3):
            continue
        ac = sum((v.get(k) or {}).get("ac") or 0 for k in ("exome", "genome"))
        an = sum((v.get(k) or {}).get("an") or 0 for k in ("exome", "genome"))
        if ac and an:
            table.setdefault(m.group(2), {}).setdefault(_AA1[m.group(3)], []).append([v["variant_id"], ac, an])
    return table


def prefetch_missense(transcripts: list[str], batch: int = 2) -> None:
    """Cache missense tables for many transcripts, several per GraphQL request.
    gnomAD's API is rate limited per request and caps query cost at 25; one
    transcript costs 11, so two fit in a request."""
    todo = [t for t in dict.fromkeys(transcripts)
            if t and not (external.CACHE_DIR / "gnomad_transcript" / f"{t}.json").exists()]
    for i in range(0, len(todo), batch):
        if i and i % 200 == 0:
            log.info("gnomAD: variant tables %d/%d", i, len(todo))
        chunk = todo[i:i + batch]
        query = " ".join(f't{j}: transcript(transcript_id: "{t}", reference_genome: GRCh38) {{ {_TX_FIELDS} }}'
                         for j, t in enumerate(chunk))
        try:
            data = _graphql("{ " + query + " }")
        except external.ExternalError:
            continue  # fetched one by one later
        for j, t in enumerate(chunk):
            if f"t{j}" in data:
                table = _missense_table((data[f"t{j}"] or {}).get("variants"))
                external._cached("gnomad_transcript", t, lambda table=table: table)


def transcript_missense(transcript: str) -> dict[str, dict[str, list]]:
    """{position: {alt residue: [[variant_id, ac, an], ...]}} for missense variants."""
    def build():
        data = _graphql(f'{{ transcript(transcript_id: "{transcript}", reference_genome: GRCh38) {{ {_TX_FIELDS} }} }}')
        return _missense_table((data.get("transcript") or {}).get("variants"))
    return external._cached("gnomad_transcript", transcript, build)


def matching_transcript(sequence: str, gene_id: str | None, symbol: str | None, hint: str | None) -> str | None:
    """First transcript whose Ensembl protein is identical to `sequence`."""
    if hint and protein_sequence(hint) == sequence:
        return hint
    # Canonical / MANE Select come first and are nearly always the match.
    others = [t for t in gene_transcripts(gene_id, symbol) if t != hint][:5]
    prefetch_proteins(others)
    return next((t for t in others if protein_sequence(t) == sequence), None)


def _transcript_for(saaps: list[SAAP]) -> str | None:
    first = saaps[0]
    symbol = (first.source_gene or "").split(";")[0].strip() or None
    return matching_transcript(first.protein_sequence, first.ensembl_gene, symbol, first.ensembl_transcript)


def _score(saaps: list[SAAP], table: dict | None) -> dict[int, tuple]:
    out = {}
    for s in saaps:
        if table is None:
            out[s.id] = ("unmapped", None, None)
            continue
        _, alt = residues(s)
        hits = [h for p in stored_positions(s) for h in table.get(str(p), {}).get(alt or "", [])]
        out[s.id] = (("present", sum(ac / an for _, ac, an in hits), ",".join(v for v, _, _ in hits))
                     if hits else ("absent", 0.0, None))
    return out


_FAILED = object()


def _safe(fn, arg):
    try:
        return fn(arg)
    except external.ExternalError:
        return _FAILED


def check_all(db: Session, *, overwrite: bool = False) -> dict:
    """Annotate every positioned human SAAP (unchecked ones only, by default)."""
    human = select(Observation.saap_id).where(Observation.species == "Homo sapiens")
    stmt = select(SAAP).where(SAAP.id.in_(human), SAAP.protein_sequence.is_not(None),
                              SAAP.protein_accession.is_not(None), SAAP.position_in_protein.is_not(None))
    if not overwrite:
        stmt = stmt.where(SAAP.gnomad_status.is_(None))
    by_protein: dict[str, list[SAAP]] = {}
    for s in db.scalars(stmt):
        by_protein.setdefault(s.protein_accession, []).append(s)

    prefetch_proteins([g[0].ensembl_transcript for g in by_protein.values()])
    transcripts: dict[str, str | None] = {}
    with ThreadPoolExecutor(WORKERS) as pool:
        for i, (acc, tx) in enumerate(zip(by_protein, pool.map(lambda g: _safe(_transcript_for, g), by_protein.values())), 1):
            transcripts[acc] = tx
            if i % 200 == 0:
                log.info("gnomAD: transcripts resolved %d/%d", i, len(by_protein))
    prefetch_missense([t for t in transcripts.values() if t and t is not _FAILED])
    results: dict[int, tuple] = {}
    for acc, group in by_protein.items():
        tx = transcripts[acc]
        if tx is _FAILED:
            continue  # network failure: left unchecked, retried next run
        table = _safe(transcript_missense, tx) if tx else None
        if table is not _FAILED:
            results.update(_score(group, table))
    for s in (x for group in by_protein.values() for x in group):
        if s.id in results:
            s.gnomad_status, s.gnomad_af, s.gnomad_variants = results[s.id]
    db.commit()
    counts = {k: sum(1 for r in results.values() if r[0] == k) for k in ("absent", "present", "unmapped")}
    return {"checked": len(results), "unchecked": sum(map(len, by_protein.values())) - len(results), **counts}
