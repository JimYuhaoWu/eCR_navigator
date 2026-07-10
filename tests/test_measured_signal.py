"""
Tests for scripts/attach_measured_signal.py — the MEASURED-external direction path
(maps an accessibility track onto an artifact's regions by max overlap), and the
review's Finding 1: uncovered regions currently collapse to signal 0, which is
indistinguishable from a measured zero and fabricates a direction on the diff.
"""
from __future__ import annotations

import os
import tempfile

import numpy as np

from _runner import add_repo_paths, run

add_repo_paths()

from attach_measured_signal import map_signal, read_intensity  # noqa: E402  (scripts/)


def _q(starts, ends, chrom="chr1"):
    return (np.array([chrom] * len(starts)),
            np.array(starts, np.int64), np.array(ends, np.int64))


def _by_chrom(starts, ends, vals, chrom="chr1"):
    """A single-chrom intensity index in the shape map_signal expects (sorted)."""
    return {chrom: (np.array(starts, np.int64),
                    np.array(ends, np.int64),
                    np.array(vals, np.float64))}


# --- map_signal: overlap + max selection -----------------------------------

def test_map_signal_max_over_overlapping_intervals():
    by = _by_chrom([10, 100, 200], [50, 160, 260], [1.0, 2.0, 3.0])
    c, s, e = _q([20, 105, 55, 40], [30, 150, 95, 110])
    out = map_signal(c, s, e, by)
    assert out[0] == 1.0                 # [20,30)  -> [10,50)=1
    assert out[1] == 2.0                 # [105,150)-> [100,160)=2
    assert np.isnan(out[2])              # [55,95)  -> no overlap -> unmeasured
    assert out[3] == 2.0                 # [40,110) -> overlaps both -> max(1,2)=2


def test_map_signal_halfopen_boundaries_do_not_overlap():
    by = _by_chrom([10, 100], [50, 160], [1.0, 2.0])
    c, s, e = _q([50, 160], [60, 200])   # touch at the exclusive end -> no overlap
    out = map_signal(c, s, e, by)
    assert np.isnan(out).all()


def test_map_signal_missing_chrom_is_unmeasured():
    by = _by_chrom([10], [50], [1.0], chrom="chr1")
    c, s, e = _q([10], [20], chrom="chr2")
    assert np.isnan(map_signal(c, s, e, by)[0])


# --- read_intensity: header-name and integer-index column selection --------

def test_read_intensity_by_name_and_by_index():
    with tempfile.TemporaryDirectory() as td:
        p = os.path.join(td, "atpm.tsv")
        with open(p, "w") as fh:
            fh.write("chrom\tstart\tend\tatpm_kidney\n")
            fh.write("chr1\t10\t50\t0.7\n")
            fh.write("chr1\t100\t160\t0.2\n")
        by_name = read_intensity(p, "atpm_kidney", has_header=True)
        by_index = read_intensity(p, "3", has_header=True)
    for by in (by_name, by_index):
        s, e, v = by["chr1"]
        assert list(s) == [10, 100] and list(e) == [50, 160]
        assert list(v) == [0.7, 0.2]


# --- Finding 1 (fixed): uncovered region is unmeasured (NaN), not a measured 0 ----

def test_uncovered_region_is_unmeasured_nan():
    """
    Region [1000,1050) has NO overlapping intensity interval, so its signal is
    'unmeasured' (NaN) — not 0.0. Otherwise a real 0.7 in the other state would diff
    to -0.7 and be reported as a confident 'close' that is really just missing data.
    """
    by = _by_chrom([10], [50], [0.7])          # nothing near 1000
    c, s, e = _q([1000], [1050])
    out = map_signal(c, s, e, by)
    assert np.isnan(out[0]), "uncovered region should be NaN, not a measured 0.0"


def test_measured_zero_stays_zero_not_nan():
    """A region that DOES overlap an interval whose value is 0.0 is a measured zero
    (accessible-but-flat), and must remain 0.0 — distinct from NaN/unmeasured."""
    by = _by_chrom([10], [50], [0.0])          # measured, value 0
    c, s, e = _q([20], [30])
    out = map_signal(c, s, e, by)
    assert out[0] == 0.0 and not np.isnan(out[0])


if __name__ == "__main__":
    run(globals())
