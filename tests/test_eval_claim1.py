"""Self-tests for the Claim 1 evaluation machinery (scripts/eval_driver_claim1.py).

The point of Claim 1 is to distinguish a driver_score that finds *drivers* from one
that merely finds *big movers* (regions with large |Delta aTPM|). These tests plant a
known confound and a known signal on synthetic data and assert the matched-background
evaluation reports the truth in each regime:

  - signal present  -> matched-background AUROC clearly > 0.5
  - confound only   -> matched-background AUROC ~ 0.5 EVEN THOUGH random-background is
                       inflated (this is the whole reason for matching)
  - label shuffle   -> AUROC ~ 0.5
  - overlap + AUROC + matched sampler behave correctly on hand-checkable inputs
"""
from _runner import run, add_repo_paths

add_repo_paths()

import numpy as np
from eval_driver_claim1 import (
    overlap_labels, matched_negative_indices, auroc, evaluate, topk_fold_enrichment,
)


def _synthetic(n=6000, n_pos=600, alpha=1.0, beta=1.0, seed=1):
    """Return (scores, labels, confound).

    confound |Delta aTPM| ~ U(0,1). Positives are drawn with a bias toward HIGH confound
    (so a naive/random-background test is confounded). score = alpha*is_pos + beta*confound
    + noise. alpha is the true driver-specific signal; beta is the big-mover confound.
    """
    rng = np.random.default_rng(seed)
    confound = rng.uniform(0, 1, n)
    # positives biased toward high confound: sample with weight = confound
    w = confound / confound.sum()
    pos_idx = rng.choice(n, size=n_pos, replace=False, p=w)
    labels = np.zeros(n, dtype=bool)
    labels[pos_idx] = True
    noise = rng.normal(0, 0.25, n)
    scores = alpha * labels.astype(float) + beta * confound + noise
    return scores, labels, confound


# ------------------------------------------------------------------ unit-level
def test_overlap_labels_halfopen_and_chrom():
    chrom = np.array(["chr1", "chr1", "chr1", "chr2"])
    start = np.array([100, 200, 400, 100])
    end = np.array([150, 250, 450, 150])
    # positive covers chr1:140-210 (overlaps region0 [100,150) and region1 [200,250))
    # and chr2:1000-2000 (overlaps nothing here)
    pc = np.array(["chr1", "chr2"]); psx = np.array([140, 1000]); pe = np.array([210, 2000])
    lab = overlap_labels(chrom, start, end, pc, psx, pe)
    assert lab.tolist() == [True, True, False, False], lab.tolist()


def test_overlap_touching_is_not_overlap():
    # half-open: positive [150,200) touches region [100,150) at the boundary -> no overlap
    chrom = np.array(["chr1"]); start = np.array([100]); end = np.array([150])
    lab = overlap_labels(chrom, start, end, np.array(["chr1"]), np.array([150]), np.array([200]))
    assert lab.tolist() == [False], lab.tolist()


def test_auroc_perfect_and_chance():
    s = np.array([0.1, 0.2, 0.3, 0.9, 0.8, 0.7])
    l = np.array([False, False, False, True, True, True])
    assert abs(auroc(s, l) - 1.0) < 1e-9
    # identical scores -> chance
    assert abs(auroc(np.ones(6), l) - 0.5) < 1e-9


def test_matched_sampler_matches_confound_distribution():
    scores, labels, confound = _synthetic()
    neg = matched_negative_indices(labels, confound, n_per_pos=3, n_bins=20, seed=0)
    assert len(neg) > 0
    # matched negatives should have a mean confound close to the positives' (biased high),
    # and clearly higher than the global negative mean.
    pos_mean = confound[labels].mean()
    neg_matched_mean = confound[neg].mean()
    global_neg_mean = confound[~labels].mean()
    assert abs(neg_matched_mean - pos_mean) < 0.08, (neg_matched_mean, pos_mean)
    assert neg_matched_mean > global_neg_mean + 0.05


# ------------------------------------------------------------------ regime tests
def test_signal_present_is_detected_above_confound():
    # true driver signal (alpha=1.0) on top of the confound -> matched AUROC clearly > 0.5
    scores, labels, confound = _synthetic(alpha=1.0, beta=1.0, seed=2)
    res = evaluate(scores, labels, confound, n_boot=400, n_perm=400, seed=3)
    assert res.auroc_matched > 0.6, res
    assert res.ci_matched[0] > 0.5, res.ci_matched          # CI excludes chance
    assert res.perm_p < 0.01, res.perm_p


def test_confound_only_is_rejected_by_matching():
    # NO driver-specific signal (alpha=0): score is PURE big-mover confound.
    # Random background is fooled (AUROC>0.5); matched background must report ~chance.
    scores, labels, confound = _synthetic(alpha=0.0, beta=1.0, seed=4)
    res = evaluate(scores, labels, confound, n_boot=400, n_perm=400, seed=5)
    assert res.auroc_random > 0.6, res.auroc_random          # naive test is confounded
    assert abs(res.auroc_matched - 0.5) < 0.06, res.auroc_matched   # matching removes it
    assert res.confound_gap > 0.1, res.confound_gap          # gap flags the confound


def test_label_shuffle_is_chance():
    scores, labels, confound = _synthetic(alpha=1.0, beta=1.0, seed=6)
    rng = np.random.default_rng(0)
    shuffled = np.zeros_like(labels)
    shuffled[rng.choice(len(labels), size=labels.sum(), replace=False)] = True
    res = evaluate(scores, shuffled, confound, n_boot=200, n_perm=400, seed=7)
    assert abs(res.auroc_matched - 0.5) < 0.06, res.auroc_matched


def _synthetic_opening(n=8000, n_pos=600, alpha=1.0, seed=1):
    """signed delta ~ U(-1,1): ~half open (>0), ~half close (<0). Positives are OPENING
    regions that are also drivers, biased toward large opening. score = alpha*is_pos +
    |delta| + noise (so magnitude is a confound that opening-matching must remove)."""
    rng = np.random.default_rng(seed)
    signed = rng.uniform(-1, 1, n)
    opening = signed > 0
    op = np.where(opening)[0]
    w = signed[op] / signed[op].sum()
    pos_idx = rng.choice(op, size=n_pos, replace=False, p=w)
    labels = np.zeros(n, dtype=bool); labels[pos_idx] = True
    scores = alpha * labels.astype(float) + np.abs(signed) + rng.normal(0, 0.25, n)
    return scores, labels, signed


def test_opening_only_detects_tf_among_equally_opening():
    # drivers score higher among opening regions, at matched opening magnitude
    scores, labels, signed = _synthetic_opening(alpha=1.0, seed=2)
    res = evaluate(scores, labels, np.abs(signed), signed=signed, opening_only=True,
                   n_boot=300, n_perm=300, seed=3)
    assert res.auroc_matched > 0.6, res.auroc_matched
    assert res.ci_matched[0] > 0.5, res.ci_matched


def test_opening_only_confound_only_is_chance():
    # no TF-specific signal (alpha=0): score is pure opening magnitude -> matched ~0.5
    scores, labels, signed = _synthetic_opening(alpha=0.0, seed=4)
    res = evaluate(scores, labels, np.abs(signed), signed=signed, opening_only=True,
                   n_boot=300, n_perm=300, seed=5)
    assert abs(res.auroc_matched - 0.5) < 0.06, res.auroc_matched


def test_opening_only_excludes_closing_positives():
    # positives placed on CLOSING regions must be dropped by the opening filter
    scores, labels, signed = _synthetic_opening(alpha=1.0, seed=6)
    # add closing positives
    closing = np.where(signed < 0)[0][:300]
    labels2 = labels.copy(); labels2[closing] = True
    res = evaluate(scores, labels2, np.abs(signed), signed=signed, opening_only=True,
                   n_boot=200, n_perm=200, seed=7)
    n_opening_pos = int((labels2 & (signed > 0)).sum())
    assert res.n_pos == n_opening_pos, (res.n_pos, n_opening_pos)


def test_topk_fold_enrichment_direction():
    # strong signal -> top-k enriched (>1); pure noise -> ~1
    scores, labels, confound = _synthetic(alpha=3.0, beta=0.0, seed=8)
    assert topk_fold_enrichment(scores, labels, k_frac=0.05) > 3.0
    noise = np.random.default_rng(0).normal(0, 1, len(labels))
    assert 0.5 < topk_fold_enrichment(noise, labels, k_frac=0.05) < 1.8


if __name__ == "__main__":
    run(globals())
