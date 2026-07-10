"""
Tests for the `direction` channel added in PR #2 (epiagent-direction):
  features.signed_delta, model.attach_direction, and the row-for-row alignment
  invariant that lets navigate.py join the two decoupled channels
  (magnitude driver_score + signed direction) by position alone.
"""
from __future__ import annotations

import csv
import os
import tempfile

import numpy as np

from _runner import add_repo_paths, run

add_repo_paths()

from ecr_navigator.features import EmbeddingArtifact, embedding_shift, signed_delta
from ecr_navigator.model import attach_direction, driver_scores
from ecr_navigator.weights import read_region_weights, write_region_weights, RegionWeight


# --- helpers ---------------------------------------------------------------

def _art(starts, emb, signal=None, assembly="mm10", state="s"):
    """Build an in-memory artifact on chr1; each region is [start, start+50)."""
    starts = np.asarray(starts, dtype=np.int64)
    return EmbeddingArtifact(
        chrom=np.array(["chr1"] * len(starts)),
        start=starts,
        end=starts + 50,
        embedding=np.asarray(emb, dtype=np.float32),
        meta={"assembly": assembly, "cell_state": state},
        signal=None if signal is None else np.asarray(signal, dtype=np.float32),
    )


def _dir_by_start(rows):
    return {int(r.start): (None if r.direction is None else round(r.direction, 6))
            for r in rows}


# --- attach_direction: the three scalings ----------------------------------

def test_attach_direction_maxabs_sign_and_scale():
    rows = driver_scores(["chr1"] * 3, [100, 200, 300], [150, 250, 350], [1, 1, 1])
    attach_direction(rows, [0.5, 0.0, -0.6], method="maxabs")  # scale = 0.6
    d = _dir_by_start(rows)
    assert d[100] == round(0.5 / 0.6, 6)   # opens, magnitude preserved
    assert d[200] == 0.0                    # flat stays 0
    assert d[300] == -1.0                   # biggest |delta| closes -> hits -1
    assert all(-1.0 <= r.direction <= 1.0 for r in rows)


def test_attach_direction_raw_clamps_to_unit():
    rows = driver_scores(["chr1"] * 4, [0, 1, 2, 3], [1, 2, 3, 4], [1, 1, 1, 1])
    attach_direction(rows, [0.3, -2.0, 1.5, -0.4], method="raw")
    assert [round(r.direction, 4) for r in rows] == [0.3, -1.0, 1.0, -0.4]


def test_attach_direction_signed_rank():
    rows = driver_scores(["chr1"] * 4, [0, 1, 2, 3], [1, 2, 3, 4], [1, 1, 1, 1])
    # mags [0.2,0.4,0.9,0.1] -> percentile ranks /3 = [1/3, 2/3, 1.0, 0.0]
    attach_direction(rows, [-0.2, 0.4, -0.9, 0.1], method="signed-rank")
    got = [round(r.direction, 6) for r in rows]
    assert got == [round(-1 / 3, 6), round(2 / 3, 6), -1.0, 0.0]


def test_attach_direction_nan_delta_left_unset_and_excluded_from_scale():
    """A NaN delta (unmeasured region) -> direction None, and it must NOT enter the
    maxabs scale, so the finite regions still scale by their own max |delta|."""
    rows = driver_scores(["chr1"] * 3, [0, 1, 2], [1, 2, 3], [1, 1, 1])
    attach_direction(rows, [0.3, float("nan"), -0.6], method="maxabs")  # scale = 0.6
    assert abs(rows[0].direction - 0.5) < 1e-9   # 0.3 / 0.6
    assert rows[1].direction is None          # unmeasured, left unset
    assert rows[2].direction == -1.0          # finite max, unaffected by the NaN


def test_attach_direction_all_zero_delta_no_div0():
    rows = driver_scores(["chr1"] * 3, [0, 1, 2], [1, 2, 3], [1, 1, 1])
    attach_direction(rows, [0.0, 0.0, 0.0], method="maxabs")
    assert [r.direction for r in rows] == [0.0, 0.0, 0.0]


def test_attach_direction_length_mismatch_raises():
    rows = driver_scores(["chr1"] * 3, [0, 1, 2], [1, 2, 3], [1, 1, 1])
    try:
        attach_direction(rows, [0.1, 0.2], method="maxabs")  # 2 != 3
    except ValueError:
        return
    assert False, "expected ValueError on length mismatch"


def test_attach_direction_unknown_method_raises():
    rows = driver_scores(["chr1"] * 1, [0], [1], [1])
    try:
        attach_direction(rows, [0.5], method="bogus")
    except ValueError:
        return
    assert False, "expected ValueError on unknown method"


# --- signed_delta ----------------------------------------------------------

def test_signed_delta_none_when_a_signal_missing():
    a = _art([100, 200], [[0, 0], [1, 1]], signal=None)
    b = _art([100, 200], [[0, 0], [1, 1]], signal=[0.2, 0.3])
    assert signed_delta(a, b) is None


def test_signed_delta_propagates_nan_from_unmeasured_state():
    """A region measured in A but NaN (unmeasured) in B diffs to NaN, so it can be
    left out of direction rather than scored as a spurious open/close."""
    a = _art([100, 200], [[0, 0], [1, 1]], signal=[0.7, 0.4])
    b = _art([100, 200], [[0, 0], [1, 1]], signal=[float("nan"), 0.9])
    delta = signed_delta(a, b)
    assert np.isnan(delta[0]) and abs(delta[1] - 0.5) < 1e-6   # float32 signal


# --- the alignment invariant (most important) ------------------------------

def test_alignment_invariant_shuffled_partial_overlap():
    """
    embedding_shift and signed_delta must return the SAME shared regions in the
    SAME order, so attach_direction can pair direction[k] with driver_score[k] by
    position. B is shuffled, drops one A-region (400) and adds an extra (500).
    """
    # A order: 100,200,300,400 ; B order: 300,500,100,200
    a = _art([100, 200, 300, 400],
             [[0, 0], [1, 0], [0, 1], [2, 2]],
             signal=[0.1, 0.5, 0.9, 0.2])
    b = _art([300, 500, 100, 200],
             [[0, 1], [9, 9], [3, 4], [1, 0]],
             signal=[0.3, 0.7, 0.6, 0.5])

    chrom, start, end, shift = embedding_shift(a, b)
    # shared regions, in A's order:
    assert list(start) == [100, 200, 300]

    rows = driver_scores(chrom, start, end, shift, method="rank")
    delta = signed_delta(a, b)
    assert delta is not None
    # delta is signal_b - signal_a for [100,200,300]:
    assert [round(float(x), 6) for x in delta] == [0.5, 0.0, -0.6]

    attach_direction(rows, delta, method="maxabs")  # scale 0.6
    d = _dir_by_start(rows)
    assert d[100] == round(0.5 / 0.6, 6)   # 100 opens
    assert d[200] == 0.0
    assert d[300] == -1.0                    # 300 closes hardest
    # row order still A-order, and each row's coords match its own region
    assert [int(r.start) for r in rows] == [100, 200, 300]


# --- contract emission (ties into weights.py) ------------------------------

def test_contract_emits_five_columns_and_roundtrips():
    rows = driver_scores(["chr1"] * 2, [100, 200], [150, 250], [0.2, 0.8])
    attach_direction(rows, [0.5, -0.5], method="raw")
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "out.tsv")
        write_region_weights(rows, path)
        with open(path) as fh:
            header = next(csv.reader(fh, delimiter="\t"))
        assert header == ["chrom", "start", "end", "driver_score", "direction"]
        back = read_region_weights(path)
    assert [round(r.direction, 4) for r in back] == [0.5, -0.5]


def test_contract_direction_absent_without_signal():
    """driver_scores alone (no attach_direction) -> 4-column, byte-compatible."""
    rows = driver_scores(["chr1"] * 2, [100, 200], [150, 250], [0.2, 0.8])
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "out.tsv")
        write_region_weights(rows, path)
        with open(path) as fh:
            header = next(csv.reader(fh, delimiter="\t"))
    assert header == ["chrom", "start", "end", "driver_score"]


def test_contract_none_direction_is_empty_field_not_zero():
    """In a directioned file, an unmeasured (None) row is an EMPTY field, distinct
    from a measured flat (0.0), and round-trips back to None."""
    rows = [RegionWeight("chr1", 100, 150, 0.5, direction=0.0),   # measured flat
            RegionWeight("chr1", 200, 250, 0.5, direction=None)]  # unmeasured
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "out.tsv")
        write_region_weights(rows, path)
        lines = [ln.rstrip("\n") for ln in open(path)]
        assert lines[0].split("\t") == ["chrom", "start", "end", "driver_score", "direction"]
        assert lines[1].endswith("\t0.0")   # measured flat -> 0.0
        assert lines[2].endswith("\t")      # unmeasured -> empty trailing field
        back = read_region_weights(path)
    assert back[0].direction == 0.0 and back[1].direction is None


def test_contract_write_clamps_out_of_range_direction():
    rows = [RegionWeight("chr1", 100, 150, 0.5, direction=2.0),
            RegionWeight("chr1", 200, 250, 0.5, direction=-3.0)]
    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "out.tsv")
        write_region_weights(rows, path)
        back = read_region_weights(path)
    assert [r.direction for r in back] == [1.0, -1.0]


if __name__ == "__main__":
    run(globals())
