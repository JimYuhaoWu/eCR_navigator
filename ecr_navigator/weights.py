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


@dataclass
class RegionWeight:
    chrom: str
    start: int
    end: int
    driver_score: float   # [0, 1], higher = more driver-like


def write_region_weights(rows: list[RegionWeight], path: str | Path) -> None:
    """Write scored regions as the contract TSV (chrom, start, end, driver_score)."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(CONTRACT_HEADER)
        for r in rows:
            score = max(0.0, min(1.0, float(r.driver_score)))  # clamp to [0,1]
            w.writerow([r.chrom, int(r.start), int(r.end), round(score, 4)])


def read_region_weights(path: str | Path) -> list[RegionWeight]:
    """Read a contract TSV back into RegionWeight rows."""
    out: list[RegionWeight] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            out.append(RegionWeight(
                chrom=row["chrom"],
                start=int(row["start"]),
                end=int(row["end"]),
                driver_score=float(row["driver_score"]),
            ))
    return out
