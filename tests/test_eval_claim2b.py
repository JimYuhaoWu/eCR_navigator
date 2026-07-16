"""Self-tests for scripts/eval_direction_claim2b.py (measured-direction trustworthiness).

The properties the analysis rests on:
  * balanced accuracy scores a TRIVIAL everything-opens predictor at 0.50 regardless of the
    opening rate -- the reason the test is two-sided;
  * a genuinely correct sign scores near 1.0 on both sets;
  * accuracy is measured vs the marginal opening rate, not 50%;
  * stratification finds the |direction| where accuracy decays to chance;
  * gene clustering resamples genes, so duplicated anchors of one gene do not inflate n.
"""
import tempfile
from pathlib import Path

from _runner import run, add_repo_paths

add_repo_paths()

import numpy as np
from eval_direction_claim2b import (assign_gene, evaluate, gene_clustered_ci,
                                    sign_accuracy, stratify_by_magnitude)


def _bed(rows):
    d = tempfile.mkdtemp()
    p = Path(d) / "a.bed"
    p.write_text("".join(f"{c}\t{s}\t{e}\t{g}\n" for c, s, e, g in rows))
    return str(p)


# ------------------------------------------------------------------ unit
def test_sign_accuracy_ignores_measured_flat():
    d = np.array([0.5, -0.3, 0.0, 0.2])          # one exact 0.0 -> no call
    acc, n = sign_accuracy(d, +1)
    assert n == 3 and abs(acc - 2 / 3) < 1e-9     # 0.5 and 0.2 open, -0.3 wrong


def test_stratify_locates_the_decay():
    # big |direction| always correct; small |direction| random -> accuracy decays to ~0.5
    rng = np.random.default_rng(0)
    big = np.full(400, 0.4)                        # correct (expected +)
    small = rng.choice([0.01, -0.01], 400)        # coin-flip sign at tiny magnitude
    d = np.concatenate([big, small])
    rows = stratify_by_magnitude(d, +1, np.array([0, 0.05, 1.01]))
    (_, _, _, acc_small), (_, _, _, acc_big) = rows
    assert acc_big == 1.0
    assert 0.4 < acc_small < 0.6                   # chance in the small-magnitude bin


def test_gene_clustering_counts_genes_not_regions():
    # 20 regions but all ONE gene, all wrong-signed: resampling genes can only ever pick that
    # gene, so accuracy is 0 with no spurious confidence from n=20.
    d = np.full(20, -0.5)                          # all close, expected +1 -> all wrong
    genes = np.array(["G"] * 20)
    m, lo, hi = gene_clustered_ci(d, +1, genes, n_boot=200)
    assert m == 0.0 and (lo != lo)                 # single gene -> CI is NaN, not [0,0]


def test_assign_gene_maps_region_to_overlapping_anchor():
    bc, bs, be, bg = (np.array(["chr1", "chr1"]), np.array([100, 500]),
                      np.array([200, 600]), np.array(["A", "B"]))
    g = assign_gene(np.array(["chr1", "chr1", "chr1"]), np.array([150, 550, 900]),
                    np.array([160, 560, 910]), bc, bs, be, bg)
    assert list(g) == ["A", "B", None]


def test_assign_gene_finds_a_longer_earlier_interval():
    """Variable-width anchors (neighborhood = gene +/-50kb): a long interval that STARTS before
    a short one but overlaps a downstream query must be found. A start-sorted sweep that stops
    at the first end<=start would miss it."""
    bc = np.array(["chr1", "chr1"]); bs = np.array([100, 120])
    be = np.array([600, 130]); bg = np.array(["LONG", "short"])   # LONG starts first, spans far
    g = assign_gene(np.array(["chr1"]), np.array([400]), np.array([410]), bc, bs, be, bg)
    assert list(g) == ["LONG"], "must not be missed behind the short interval"


# ------------------------------------------------------------------ integration
def _synth(opening_rate, dest_sign_correct, source_sign_correct, n=6000, seed=1):
    """A contract where background opens at `opening_rate`, plus dest anchors (expected +) and
    source anchors (expected -) whose measured sign is correct with the given probability."""
    rng = np.random.default_rng(seed)
    chrom = np.array(["chr1"] * n)
    start = np.arange(n) * 100
    end = start + 50
    direction = np.where(rng.random(n) < opening_rate, rng.uniform(0, 1, n),
                         -rng.uniform(0, 1, n))
    dest = list(range(0, 60)); src = list(range(60, 120))
    for i in dest:                                 # expected OPEN
        direction[i] = rng.uniform(0.1, 1) * (1 if rng.random() < dest_sign_correct else -1)
    for i in src:                                  # expected CLOSE
        direction[i] = -rng.uniform(0.1, 1) * (1 if rng.random() < source_sign_correct else -1)
    dbed = _bed([("chr1", start[i], end[i], f"D{i}") for i in dest])
    sbed = _bed([("chr1", start[i], end[i], f"S{i}") for i in src])
    return chrom, start, end, direction, dbed, sbed


def test_trivial_everything_opens_scores_balanced_half():
    # opening_rate 0.9, and source anchors ALSO (wrongly) open -> destination looks perfect,
    # source fails; balanced accuracy must expose it near 0.5, NOT the inflated 0.9.
    c, s, e, d, dbed, sbed = _synth(0.9, dest_sign_correct=1.0, source_sign_correct=0.0)
    out = evaluate(c, s, e, d, dbed, sbed)
    assert out["sets"]["destination"]["accuracy"] > 0.95
    assert out["sets"]["source"]["accuracy"] < 0.05
    assert abs(out["balanced_accuracy"] - 0.5) < 0.05


def test_correct_direction_scores_high_on_both_sets():
    c, s, e, d, dbed, sbed = _synth(0.55, dest_sign_correct=0.95, source_sign_correct=0.95)
    out = evaluate(c, s, e, d, dbed, sbed)
    assert out["balanced_accuracy"] > 0.9
    assert out["sets"]["source"]["beats_base"]     # closes beat the 1-opening_rate base


def test_accuracy_is_measured_against_opening_rate_not_half():
    c, s, e, d, dbed, sbed = _synth(0.8, dest_sign_correct=0.85, source_sign_correct=0.85)
    out = evaluate(c, s, e, d, dbed, sbed)
    assert abs(out["opening_rate"] - 0.8) < 0.03
    assert abs(out["sets"]["destination"]["base_rate"] - 0.8) < 0.03
    assert abs(out["sets"]["source"]["base_rate"] - 0.2) < 0.03


if __name__ == "__main__":
    run(globals())
