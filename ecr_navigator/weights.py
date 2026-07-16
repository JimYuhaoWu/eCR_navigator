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
UNMEASURED = ""     # how an unmeasured direction is serialized: an EMPTY field
PRECISION = 4       # decimal places for every score/direction the contract emits


# --- contract encoding (the rules every file in a run bundle obeys) -------------------
# These live here, in the contract I/O module, because nominations.tsv (ecr_navigator/
# nominate.py) must encode its columns EXACTLY as weights.tsv does. Duplicating the rules
# per writer let them drift, and the unmeasured-vs-flat rule below is one a consumer
# silently misreads if it drifts.

def fmt_score(value: float) -> float:
    """Serialize a [0,1] score (driver_score, nomination_score). Clamped, not trusted: a
    model emitting an unnormalized score must not leak an out-of-range value into a file
    whose contract declares [0,1]."""
    return round(max(0.0, min(1.0, float(value))), PRECISION)


def fmt_direction(value: float | None):
    """Serialize a [-1,1] direction. `None` (UNMEASURED) becomes an EMPTY field, NEVER 0.0
    -- 0.0 means measured-and-flat (accessible, no change between states), so collapsing
    the two fabricates an open/close call from a missing measurement."""
    if value is None:
        return UNMEASURED
    return round(max(-1.0, min(1.0, float(value))), PRECISION)


def parse_direction(raw: str | None) -> float | None:
    """Read a direction field back. Empty (or column absent) -> None = unmeasured. The
    inverse of fmt_direction; keep the two together so the round-trip cannot drift."""
    return float(raw) if raw not in (None, UNMEASURED) else None


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
            out = [r.chrom, int(r.start), int(r.end), fmt_score(r.driver_score)]
            if has_dir:
                out.append(fmt_direction(r.direction))
            w.writerow(out)


def read_region_weights(path: str | Path) -> list[RegionWeight]:
    """Read a contract TSV back into RegionWeight rows (direction if present)."""
    out: list[RegionWeight] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(RegionWeight(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                driver_score=float(row["driver_score"]),
                direction=parse_direction(row.get(DIRECTION_COL)),
            ))
    return out
