"""
Nomination: turn scored regions into the run bundle eCR_predictor consumes.

The navigator's job ends at "which loci matter, in what order". This module applies the
nomination policy (docs/validation_summary.md) to a scored region set and writes the run
bundle (docs/run_bundle_contract.md):

    <run_id>/
      manifest.json      the verdict (Gate-1 admissibility, Gate-2 PRIMARY) + provenance
      weights.tsv        DENSE, unchanged -> eCR_predictor offtarget.py Tier-2
      nominations.tsv    SPARSE, ranked   -> eCR_predictor target.py

Two things this module exists to enforce:

  * The PRIMARY score is resolved HERE, not downstream. Gate 2 picks GET or signed-Delta
    per transition; the bundle carries one `nomination_score` and names its source only in
    the manifest. A future self-trained model that wins Gate 2 changes nothing downstream.

  * REFUSAL IS AN OUTPUT. A transition that fails Gate 1 yields a valid bundle with ZERO
    nominations and a reason -- not a confident-looking top-1% of noise. On the v1
    benchmark panel that is the correct answer for 3 of 6 transitions; MyoD is the case
    that motivates it (GET *anti*-informative there, AUROC 0.288).

Gate 1 rejects the transition; Gate 2 only chooses WHICH score to nominate from. So
PRIMARY=signed-Delta is a normal, non-refusing outcome (e.g. MEF->mES: admissible, but
signed-Delta already captures the signal).
"""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .model import _normalize
from .weights import RegionWeight, write_region_weights

BUNDLE_VERSION = "1"
DEFAULT_TOP_FRAC = 0.01   # the validated confidence band (GET top ~1%); see below
NOMINATION_HEADER = ["chrom", "start", "end", "rank", "nomination_score", "direction"]


@dataclass
class Nomination:
    chrom: str
    start: int
    end: int
    rank: int                     # 1-based, by nomination_score descending
    nomination_score: float       # [0,1], from the Gate-2 PRIMARY scorer
    direction: float | None       # [-1,1] measured signed-delta; None = unmeasured


def nominate(rows: list[RegionWeight], gate1: dict | None, gate2: dict | None,
             top_frac: float = DEFAULT_TOP_FRAC, score_norm: str = "rank"):
    """Apply the nomination policy. Returns (nominations, manifest_nomination_block).

    `gate1` / `gate2` are the manifest sections (see docs/run_bundle_contract.md), or None
    if that gate did not run. Refuses -- zero nominations, with a reason -- unless BOTH
    gates ran and Gate 1 admitted:

      * Gate 1 missing        -> admissibility unverified, so nothing is trustworthy.
      * Gate 1 admit = False  -> the transition itself is unusable.
      * Gate 2 missing        -> no anchors, so no score was selected. Deliberately no
                                 `--assume-primary` fallback: defaulting to an unmeasured
                                 GET top-1% reintroduces exactly the MyoD failure.

    NOTE this is stricter than scripts/preflight.py, whose Gate 1 is optional (there, Gate 2
    is "the decision"). preflight is a diagnostic; nomination is the production path, and an
    unverified endpoint pair should not produce targets.

    `top_frac` is the validated CONFIDENCE BAND, not a design budget -- GET's top ~1% is
    enriched ~9-10x on a strong clean transition. That is ~3,300 regions for human iN, far
    more than anyone builds eCRs against; the consumer sets its own budget by cutting on
    `rank` (which is why rank is a column).
    """
    refusal = _refusal_reason(gate1, gate2)
    if refusal is not None:
        return [], {"score_source": None, "score_norm": None, "top_frac": None,
                    "n_nominated": 0, "refused": True, "refusal_reason": refusal}

    primary = gate2["primary"]
    scored, norm = _primary_scores(rows, primary, score_norm)

    # stable sort: equal scores keep input order, so ranks are deterministic. Ties are
    # expected and unavoidable under score_norm='rank' -- driver_score is a percentile, so
    # a top_frac band spans [1-top_frac, 1.0] by construction.
    ranked = sorted(scored, key=lambda t: -t[1])
    # ceil, not round: round() would banker's-round a 50-region universe at top_frac=0.01
    # down to ZERO nominations while still reporting refused=False -- a bundle state the
    # contract does not allow (refusal is the only way to nominate nothing).
    k = max(1, math.ceil(len(scored) * top_frac))
    # full precision in memory; rounding is a serialization concern (as in weights.py)
    noms = [Nomination(chrom=r.chrom, start=r.start, end=r.end, rank=i,
                       nomination_score=float(s), direction=r.direction)
            for i, (r, s) in enumerate(ranked[:k], 1)]

    return noms, {"score_source": primary, "score_norm": norm, "top_frac": top_frac,
                  "n_nominated": len(noms), "refused": False, "refusal_reason": None}


def _refusal_reason(gate1: dict | None, gate2: dict | None) -> str | None:
    if gate1 is None:
        return ("Gate-1 did not run (no per-replicate endpoint matrix): admissibility is "
                "unverified, so no score is trustworthy.")
    if not gate1.get("admit"):
        return "Gate-1 REJECT: " + "; ".join(gate1.get("reasons") or ["not admissible"])
    if gate2 is None:
        return ("Gate-2 did not run (no target-cell master-TF anchors): no PRIMARY score "
                "was selected, and there is no fallback.")
    return None


def _primary_scores(rows: list[RegionWeight], primary: str, score_norm: str):
    """Regions paired with their PRIMARY score, plus the normalization actually used.

    PRIMARY=signed-Delta scores by |measured signed-delta| -- which `direction` already
    carries (its normalizations are monotonic in |delta|, so the ranking is preserved).
    Unmeasured regions (direction=None) have no signed-Delta score and are dropped.
    """
    if primary == "signed-Delta":
        measured = [r for r in rows if r.direction is not None]
        if not measured:
            raise ValueError("PRIMARY=signed-Delta but no region carries a direction")
        mag = _normalize(np.abs(np.asarray([r.direction for r in measured], float)), "rank")
        return list(zip(measured, mag)), "rank"
    # any model scorer: driver_score is already normalized by navigate.py
    return [(r, r.driver_score) for r in rows], score_norm


def write_nominations(noms: list[Nomination], path: str | Path) -> None:
    """Write nominations.tsv. A refusal writes the header and nothing else -- that is a
    valid bundle, not an error. An unmeasured direction is an EMPTY field, never 0.0
    (which would mean measured-and-flat), matching the region-weight contract."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(NOMINATION_HEADER)
        for n in noms:
            d = ("" if n.direction is None
                 else round(max(-1.0, min(1.0, float(n.direction))), 4))
            score = max(0.0, min(1.0, float(n.nomination_score)))   # clamp, as weights.py
            w.writerow([n.chrom, int(n.start), int(n.end), n.rank,
                        round(score, 4), d])


def read_nominations(path: str | Path) -> list[Nomination]:
    """Read nominations.tsv back (empty list for a refusal bundle)."""
    out: list[Nomination] = []
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            raw = row.get("direction")
            out.append(Nomination(
                chrom=row["chrom"], start=int(row["start"]), end=int(row["end"]),
                rank=int(row["rank"]), nomination_score=float(row["nomination_score"]),
                direction=float(raw) if raw not in (None, "") else None,
            ))
    return out


def write_bundle(outdir: str | Path, manifest: dict, weights: list[RegionWeight],
                 noms: list[Nomination]) -> None:
    """Write the three-part run bundle. `weights` is the FULL region universe (emitted even
    for a refusal -- Tier-2 off-target weighting only needs relative accessibility
    importance, not a trustworthy driver ranking).

    `bundle_version` is stamped here rather than taken from the caller: the writer is what
    determines the on-disk format, so it is the only thing that can honestly declare it."""
    d = Path(outdir)
    d.mkdir(parents=True, exist_ok=True)
    write_region_weights(weights, d / "weights.tsv")
    write_nominations(noms, d / "nominations.tsv")
    with open(d / "manifest.json", "w") as fh:
        json.dump({**manifest, "bundle_version": BUNDLE_VERSION}, fh, indent=2)
