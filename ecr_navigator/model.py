"""
Zero-shot driver scorer: turn per-region embedding shift into a driver_score in
[0, 1] for the region-weight contract.

Two-point (MEF vs mES) regime — no time-course, so no supervised head. The score
is a normalization of the model's per-region embedding shift between the two cell
states (model-agnostic: ChromBERT, GET, ATACformer, ChromFound, … all feed the same
shift through here). When time-course data arrives, a supervised/hybrid head slots
in here behind the same output contract.
"""
from __future__ import annotations

import numpy as np

from .weights import RegionWeight


def _normalize(shift: np.ndarray, method: str) -> np.ndarray:
    if method == "rank":
        # percentile rank -> uniform [0,1], robust to outlier shifts
        order = shift.argsort().argsort()
        return order / max(len(shift) - 1, 1)
    if method == "minmax":
        lo, hi = float(shift.min()), float(shift.max())
        if hi <= lo:
            return np.zeros_like(shift)
        return (shift - lo) / (hi - lo)
    raise ValueError("unknown normalization %r (use 'rank' or 'minmax')" % method)


def driver_scores(chrom, start, end, shift, method: str = "rank") -> list[RegionWeight]:
    """
    Map embedding shifts to driver_score in [0,1] and return contract rows.

    method:
      'rank'   — percentile rank of the shift (default; robust, uniform spread).
      'minmax' — linear min-max scaling of the raw shift magnitude.
    """
    score = _normalize(np.asarray(shift, dtype=float), method)
    return [
        RegionWeight(chrom=str(c), start=int(s), end=int(e), driver_score=float(v))
        for c, s, e, v in zip(chrom, start, end, score)
    ]
