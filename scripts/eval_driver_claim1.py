#!/usr/bin/env python
"""
Claim 1 evaluation: is driver_score informative (drivers vs passengers)?

The core methodological point (see docs/direction.md and the progress notes):
a driver_score is an embedding-shift MAGNITUDE, so "high score overlaps regions
that changed a lot" is nearly tautological. To claim we found *drivers* rather
than merely *big movers*, we test enrichment of a labeled driver set (e.g. master-TF
binding regions) against a background MATCHED on the change magnitude |Delta aTPM|.

This module is pure numpy so it runs anywhere (no scipy/pandas required). It is
exercised end-to-end by tests/test_eval_claim1.py on synthetic data with a planted
signal + a planted confound, so the machinery is trusted before it touches real data.

Primary readout: AUROC of driver_score for classifying P (positives) vs a
change-magnitude-matched negative set, with a bootstrap CI, plus the gap to a
RANDOM negative set (that gap quantifies the big-mover confound). Nulls:
label-shuffle and (upstream) same-state should both collapse to AUROC ~ 0.5.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass

import numpy as np


# ---------------------------------------------------------------- interval overlap
def overlap_labels(chrom, start, end, pos_chrom, pos_start, pos_end):
    """Boolean per query region: does it overlap ANY positive interval (half-open)?

    Per-chrom sweep: sort positives by start, and for each query find positives whose
    start < query_end via searchsorted, then check the max end-so-far reaches query_start.
    O((n+m) log m). All inputs are 1-D array-likes; coordinates are half-open [start,end).
    """
    chrom = np.asarray(chrom)
    start = np.asarray(start, dtype=np.int64)
    end = np.asarray(end, dtype=np.int64)
    out = np.zeros(len(chrom), dtype=bool)
    if len(pos_chrom) == 0:
        return out
    pos_chrom = np.asarray(pos_chrom)
    pos_start = np.asarray(pos_start, dtype=np.int64)
    pos_end = np.asarray(pos_end, dtype=np.int64)

    for c in np.unique(chrom):
        q = np.where(chrom == c)[0]
        p = np.where(pos_chrom == c)[0]
        if len(p) == 0:
            continue
        ps = pos_start[p]
        pe = pos_end[p]
        order = np.argsort(ps, kind="mergesort")
        ps = ps[order]
        pe = pe[order]
        # prefix max of end so we can ask "does any positive starting before q_end reach q_start?"
        pe_cummax = np.maximum.accumulate(pe)
        qs = start[q]
        qe = end[q]
        # positives with ps < qe : indices [0, k)
        k = np.searchsorted(ps, qe, side="left")
        has = k > 0
        # among those k, the largest end is pe_cummax[k-1]; overlap iff that end > qs
        idx = np.clip(k - 1, 0, len(pe_cummax) - 1)
        reach = pe_cummax[idx] > qs
        out[q] = has & reach
    return out


# ------------------------------------------------------------- matched background
def matched_negative_indices(labels, confound, n_per_pos=3, n_bins=20, seed=0):
    """Sample negatives whose `confound` (e.g. |Delta aTPM|) matches the positives'.

    Bin `confound` into quantile bins; within each bin draw n_per_pos * (#pos in bin)
    negatives (or all available if fewer). This makes the negative set's confound
    distribution track the positives' bin-for-bin, so downstream score differences
    cannot be explained by the confound alone.
    """
    labels = np.asarray(labels, dtype=bool)
    confound = np.asarray(confound, dtype=float)
    rng = np.random.default_rng(seed)

    finite = np.isfinite(confound)
    # quantile bin edges from the finite confound values
    qs = np.linspace(0, 1, n_bins + 1)
    edges = np.unique(np.quantile(confound[finite], qs))
    if len(edges) < 2:
        edges = np.array([confound[finite].min(), confound[finite].max() + 1e-9])
    bin_id = np.digitize(confound, edges[1:-1], right=False)  # 0..len(edges)-2
    bin_id[~finite] = -1

    chosen = []
    for b in np.unique(bin_id):
        if b < 0:
            continue
        in_bin = bin_id == b
        pos_pool = np.where(in_bin & labels)[0]
        neg_pool = np.where(in_bin & ~labels)[0]
        if len(pos_pool) == 0 or len(neg_pool) == 0:
            continue
        want = min(len(neg_pool), n_per_pos * len(pos_pool))
        chosen.append(rng.choice(neg_pool, size=want, replace=False))
    if not chosen:
        return np.array([], dtype=int)
    return np.concatenate(chosen)


# --------------------------------------------------------------------------- auroc
def _avg_ranks(x):
    """Average ranks (1-based), ties averaged."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    sx = x[order]
    i = 0
    n = len(x)
    while i < n:
        j = i + 1
        while j < n and sx[j] == sx[i]:
            j += 1
        ranks[order[i:j]] = (i + j - 1) / 2.0 + 1.0  # average of 1-based ranks
        i = j
    return ranks


def auroc(scores, labels):
    """Mann-Whitney AUROC that P (label True) scores higher than N. NaN if degenerate."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    n_pos = int(labels.sum())
    n_neg = int((~labels).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    r = _avg_ranks(scores)
    return (r[labels].sum() - n_pos * (n_pos + 1) / 2.0) / (n_pos * n_neg)


def bootstrap_auroc_ci(scores, labels, n_boot=2000, seed=0, alpha=0.05):
    """Percentile bootstrap CI for AUROC, resampling positives and negatives separately."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    rng = np.random.default_rng(seed)
    pi = np.where(labels)[0]
    ni = np.where(~labels)[0]
    if len(pi) == 0 or len(ni) == 0:
        return (float("nan"), float("nan"))
    vals = np.empty(n_boot)
    for b in range(n_boot):
        ps = rng.choice(pi, size=len(pi), replace=True)
        ns = rng.choice(ni, size=len(ni), replace=True)
        idx = np.concatenate([ps, ns])
        lab = np.concatenate([np.ones(len(ps), bool), np.zeros(len(ns), bool)])
        vals[b] = auroc(scores[idx], lab)
    return (float(np.quantile(vals, alpha / 2)), float(np.quantile(vals, 1 - alpha / 2)))


def shuffle_null_auroc(scores, labels, n_perm=1000, seed=0):
    """Distribution of AUROC under label permutation; should center on 0.5.
    Returns (mean, one-sided permutation p-value that observed >= null)."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    rng = np.random.default_rng(seed)
    obs = auroc(scores, labels)
    n = len(labels)
    k = int(labels.sum())
    null = np.empty(n_perm)
    for i in range(n_perm):
        perm = np.zeros(n, dtype=bool)
        perm[rng.choice(n, size=k, replace=False)] = True
        null[i] = auroc(scores, perm)
    p = (1 + np.sum(null >= obs)) / (n_perm + 1)
    return float(null.mean()), float(p)


def topk_fold_enrichment(scores, labels, k_frac=0.05):
    """Fold-enrichment of positives in the top k_frac of scores vs base rate."""
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    n = len(scores)
    k = max(1, int(round(k_frac * n)))
    top = np.argsort(-scores, kind="mergesort")[:k]
    base = labels.mean()
    if base == 0:
        return float("nan")
    return float(labels[top].mean() / base)


# --------------------------------------------------------------------------- report
@dataclass
class Claim1Result:
    n_pos: int
    n_matched_neg: int
    auroc_matched: float
    ci_matched: tuple
    auroc_random: float
    confound_gap: float          # auroc_random - auroc_matched (bigger = more big-mover confound)
    shuffle_null_mean: float
    perm_p: float
    topk_fold: float


def evaluate(scores, labels, confound, signed=None, opening_only=False, open_thresh=0.0,
             n_per_pos=3, n_bins=20, seed=0, n_boot=2000, n_perm=1000, k_frac=0.05):
    """Enrichment of driver_score for positives vs a confound-matched background.

    Default: positives vs a |confound|-matched negative set over all regions.
    opening_only=True: restrict to regions OPENING toward the end state
    (signed delta > open_thresh) and match negatives on the opening magnitude, so
    positives (opening AND labelled) are compared only to non-labelled regions that
    open by a similar amount — isolating the label (TF binding) from the accessibility
    gain. Requires `signed` (the signed delta; `confound` is its magnitude).
    """
    scores = np.asarray(scores, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    confound = np.asarray(confound, dtype=float)

    keep = np.isfinite(confound)
    if opening_only:
        if signed is None:
            raise ValueError("opening_only requires signed delta")
        keep = keep & (np.asarray(signed, dtype=float) > open_thresh)
    scores, labels, confound = scores[keep], labels[keep], confound[keep]

    pos = np.where(labels)[0]
    neg_matched = matched_negative_indices(labels, confound, n_per_pos, n_bins, seed)
    # random background of the same size, drawn from all non-positives
    rng = np.random.default_rng(seed + 1)
    neg_pool = np.where(~labels)[0]
    neg_random = rng.choice(neg_pool, size=min(len(neg_matched), len(neg_pool)), replace=False)

    def _auroc_against(neg):
        idx = np.concatenate([pos, neg])
        lab = np.concatenate([np.ones(len(pos), bool), np.zeros(len(neg), bool)])
        return idx, lab, auroc(scores[idx], lab)

    idx_m, lab_m, au_m = _auroc_against(neg_matched)
    _, _, au_r = _auroc_against(neg_random)
    ci_m = bootstrap_auroc_ci(scores[idx_m], lab_m, n_boot=n_boot, seed=seed)
    null_mean, perm_p = shuffle_null_auroc(scores[idx_m], lab_m, n_perm=n_perm, seed=seed)
    fold = topk_fold_enrichment(scores, labels, k_frac=k_frac)

    return Claim1Result(
        n_pos=len(pos), n_matched_neg=len(neg_matched),
        auroc_matched=au_m, ci_matched=ci_m, auroc_random=au_r,
        confound_gap=au_r - au_m, shuffle_null_mean=null_mean,
        perm_p=perm_p, topk_fold=fold)


# ------------------------------------------------------------------------ file I/O
def load_contract(path):
    """Read chrom,start,end,driver_score[,direction] TSV (header optional)."""
    chrom, start, end, score, direction = [], [], [], [], []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith("#"):
                continue
            f = ln.split("\t")
            if f[0] in ("chrom", "chr") and not f[1].lstrip("-").isdigit():
                continue  # header
            chrom.append(f[0]); start.append(int(f[1])); end.append(int(f[2]))
            score.append(float(f[3]))
            direction.append(float(f[4]) if len(f) > 4 and f[4] != "" else np.nan)
    return (np.array(chrom), np.array(start), np.array(end),
            np.array(score, float), np.array(direction, float))


def load_bed(path):
    chrom, start, end = [], [], []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.split("\t")
            chrom.append(f[0]); start.append(int(f[1])); end.append(int(f[2]))
    return np.array(chrom), np.array(start), np.array(end)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contract", required=True, help="driver-score contract TSV")
    ap.add_argument("--positives", required=True, help="driver-label BED (e.g. master-TF peaks)")
    ap.add_argument("--confound", default=None,
                    help="TSV chrom,start,end,value giving |Delta aTPM| per region; "
                         "if omitted, |direction| from the contract is used")
    ap.add_argument("--n-per-pos", type=int, default=3)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--k-frac", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--opening-only", action="store_true",
                    help="restrict to regions OPENING toward the end state (signed delta > "
                         "--open-thresh) and match on opening magnitude; tests the label "
                         "(TF binding) among equally-opening regions")
    ap.add_argument("--open-thresh", type=float, default=0.0,
                    help="signed-delta threshold for 'opening' (default 0.0)")
    args = ap.parse_args()

    chrom, start, end, score, direction = load_contract(args.contract)
    pc, ps, pe = load_bed(args.positives)
    labels = overlap_labels(chrom, start, end, pc, ps, pe)

    if args.confound:
        cc, cs, ce, cv, _ = load_contract(args.confound)  # reuse: value in col 4
        # join by exact coordinate key; keep SIGNED delta (sign selects opening regions)
        key = {(a, int(b), int(d)): v for a, b, d, v in zip(cc, cs, ce, cv)}
        signed = np.array([key.get((a, int(b), int(d)), np.nan)
                           for a, b, d in zip(chrom, start, end)])
    else:
        signed = direction
    confound = np.abs(signed)
    if not np.isfinite(confound).any():
        raise SystemExit("no finite confound values; pass --confound with Delta aTPM")

    res = evaluate(score, labels, confound, signed=signed,
                   opening_only=args.opening_only, open_thresh=args.open_thresh,
                   n_per_pos=args.n_per_pos, n_bins=args.n_bins, seed=args.seed,
                   k_frac=args.k_frac)
    if args.opening_only:
        print(f"[opening-only: signed delta > {args.open_thresh}]")
    print(f"positives (regions overlapping label): {res.n_pos}")
    print(f"matched negatives:                     {res.n_matched_neg}")
    print(f"AUROC vs matched background:  {res.auroc_matched:.3f}  "
          f"CI[{res.ci_matched[0]:.3f}, {res.ci_matched[1]:.3f}]")
    print(f"AUROC vs random background:   {res.auroc_random:.3f}")
    print(f"confound gap (random-matched):{res.confound_gap:+.3f}  "
          f"(large => driver_score mostly tracks change magnitude)")
    print(f"label-shuffle null mean:      {res.shuffle_null_mean:.3f}  perm p={res.perm_p:.4g}")
    print(f"top-{args.k_frac:.0%} fold enrichment:      {res.topk_fold:.2f}x")


if __name__ == "__main__":
    main()
