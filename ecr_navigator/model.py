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


def attach_direction(rows: list[RegionWeight], delta, method: str = "maxabs") -> None:
    """
    Set each row's signed `direction ∈ [-1, 1]` from the per-region accessibility
    change `delta` (= signal_b - signal_a; +opens / -closes), in place. `rows` and
    `delta` must be in the same region order (both come from the same shared-region
    alignment). This is the direction CHANNEL — orthogonal to the magnitude
    driver_score, which is untouched.

    method (magnitude scaling only; the SIGN is always the sign of delta):
      'maxabs' — delta / max(|delta|); preserves relative magnitudes, bounds to
                 [-1,1]. Default; works for any signal scale.
      'raw'    — clamp delta itself to [-1,1]; use when the signal is already a
                 probability (delta in [-1,1]), keeping direction interpretable as
                 ΔP(accessible).
      'signed-rank' — sign(delta) * percentile-rank(|delta|); robust to outliers,
                 mirrors the 'rank' magnitude normalization.
    """
    delta = np.asarray(delta, dtype=float)
    if method == "maxabs":
        scale = np.abs(delta).max()
        direction = delta / scale if scale > 0 else np.zeros_like(delta)
    elif method == "raw":
        direction = np.clip(delta, -1.0, 1.0)
    elif method == "signed-rank":
        mag = np.abs(delta)
        rank = mag.argsort().argsort() / max(len(mag) - 1, 1)
        direction = np.sign(delta) * rank
    else:
        raise ValueError("unknown direction method %r "
                         "(use 'maxabs', 'raw', or 'signed-rank')" % method)
    if len(direction) != len(rows):
        raise ValueError("direction length %d != %d rows" % (len(direction), len(rows)))
    for r, d in zip(rows, direction):
        r.direction = float(d)
