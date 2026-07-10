"""
Region-weight contract I/O — the stable interface eCR_predictor consumes.

This is deliberately model-agnostic: whatever scorer produces driver scores
(supervised, zero-shot LM, hybrid), it emits rows through here, and eCR_predictor
reads them back the same way. See docs/region_weight_contract.md.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

CONTRACT_HEADER = ["chrom", "start", "end", "driver_score"]
DIRECTION_COL = "direction"


@dataclass
class RegionWeight:
    chrom: str
    start: int
    end: int
    driver_score: float             # [0, 1], magnitude — higher = more driver-like
    direction: float | None = None  # [-1, 1], signed: +1 open, -1 close; None if the
    #                                 model resolves importance but not direction


def write_region_weights(rows: list[RegionWeight], path: str | Path) -> None:
    """Write scored regions as the contract TSV.

    Always emits `chrom, start, end, driver_score` (driver_score is the stable
    primary column eCR_predictor consumes). A fifth `direction` column is added
    ONLY when at least one row carries a direction — direction-capable models
    (e.g. EpiAgent's signed predicted-accessibility) fill it; region-only models
    omit it, keeping the 4-column form byte-compatible with existing consumers.

    A row whose `direction is None` in a directioned file (an unmeasured region in an
    otherwise direction-capable run) is written as an EMPTY field, not 0.0 — so
    "unmeasured" stays distinct from a measured "flat" (0.0). read_region_weights maps
    the empty field back to None.
    """
    has_dir = any(r.direction is not None for r in rows)
    header = CONTRACT_HEADER + ([DIRECTION_COL] if has_dir else [])
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        for r in rows:
            score = max(0.0, min(1.0, float(r.driver_score)))  # clamp to [0,1]
            out = [r.chrom, int(r.start), int(r.end), round(score, 4)]
            if has_dir:
                d = ("" if r.direction is None
                     else round(max(-1.0, min(1.0, float(r.direction))), 4))
                out.append(d)
            w.writerow(out)


def read_region_weights(path: str | Path) -> list[RegionWeight]:
    """Read a contract TSV back into RegionWeight rows (direction if present)."""
    out: list[RegionWeight] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            raw_dir = row.get(DIRECTION_COL)
            out.append(RegionWeight(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                driver_score=float(row["driver_score"]),
                direction=float(raw_dir) if raw_dir not in (None, "") else None,
            ))
    return out
