#!/usr/bin/env python
"""
Claim 2A evaluation: does `driver_score` beat / add over a plain signed-Delta
accessibility baseline at recovering drivers?

Claim 1 showed a capable model's driver_score recovers master-TF loci even after
matching the |Delta aTPM| MAGNITUDE — but the surviving signal is "largely
directional" (it mostly says *these regions open*). So the honest next question is
whether driver_score carries information beyond the signed-Delta you already measure
from the ATAC data. (This is the same test that decided AlphaGenome's two-state case:
if signed-Delta from your own data dominates, the model's real contribution is
regulatory-region prioritization, not direction.)

This module promotes signed-Delta from Claim 1's *confound* to a *competing scorer*
and reports, on ONE matched sample (positives + a |signed-Delta|-matched, opening-only
background — the same construction as Claim 1):

  1. AUROC(driver_score)  vs  AUROC(signed_delta)         -- head to head
  2. Delta AUROC = AUROC(driver) - AUROC(signed), paired bootstrap CI (same
     resampled indices applied to both scorers, so the CI is of the DIFFERENCE)
  3. Incremental likelihood-ratio test: does driver_score improve a logistic model
     that already has signed_delta?  LR = 2*(ll_full - ll_reduced), with a
     permutation p-value (shuffle driver_score, refit) and the standardized
     driver_score coefficient.
  4. top-5% fold enrichment for both scorers.

Verdict logic: driver_score's AUROC CI above signed_delta's AND a significant
incremental LR -> the model adds regulatory information beyond "what opens" (Claim 2A
supported). Overlapping CIs / null increment -> driver_score is dominated by
signed-Delta for this endpoint pair (the honest null).

Pure numpy (no scipy/sklearn). The primitives (overlap, matched sampler, AUROC,
top-k) are shared with Claim 1 by importing scripts/eval_driver_claim1.py; the
logistic fit + LR test are here. Exercised end-to-end by tests/test_eval_claim2.py on
synthetic data with a planted directional signal AND a planted signed-Delta confound.
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_driver_claim1 import (          # noqa: E402  (shared primitives)
    auroc, load_bed, load_contract, matched_negative_indices, overlap_labels,
    topk_fold_enrichment,
)


def _standardize(x):
    x = np.asarray(x, dtype=float)
    sd = x.std()
    return (x - x.mean()) / sd if sd > 0 else np.zeros_like(x)


def logistic_fit(X, y, l2=1e-6, iters=100, tol=1e-9):
    """Newton-Raphson logistic regression. X includes an intercept column; y in {0,1}.

    A tiny L2 ridge keeps the Hessian invertible under separation. Returns
    (coefficients, log-likelihood).
    """
    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    n, k = X.shape
    beta = np.zeros(k)
    for _ in range(iters):
        p = 1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30)))
        w = np.clip(p * (1.0 - p), 1e-9, None)
        grad = X.T @ (y - p) - l2 * beta
        hess = -(X.T * w) @ X - l2 * np.eye(k)
        step = np.linalg.solve(hess, grad)
        beta_new = beta - step
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    p = np.clip(1.0 / (1.0 + np.exp(-np.clip(X @ beta, -30, 30))), 1e-12, 1 - 1e-12)
    ll = float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))
    return beta, ll


def incremental_lr_test(driver, signed, labels, n_perm=1000, seed=0):
    """Does driver_score improve a logistic model that already has signed_delta?

    Reduced model: label ~ 1 + signed. Full: label ~ 1 + signed + driver (both
    standardized). LR = 2*(ll_full - ll_reduced); permutation p-value shuffles driver
    among the sample and refits (null distribution of the LR statistic). Returns
    (lr_stat, standardized_driver_coef, perm_p).
    """
    y = np.asarray(labels, dtype=float)
    z_s = _standardize(signed)
    z_d = _standardize(driver)
    ones = np.ones(len(y))
    x_red = np.column_stack([ones, z_s])
    _, ll_red = logistic_fit(x_red, y)
    beta_full, ll_full = logistic_fit(np.column_stack([ones, z_s, z_d]), y)
    lr = 2.0 * (ll_full - ll_red)
    coef_driver = float(beta_full[-1])

    rng = np.random.default_rng(seed)
    null = np.empty(n_perm)
    for i in range(n_perm):
        _, ll_i = logistic_fit(np.column_stack([ones, z_s, rng.permutation(z_d)]), y)
        null[i] = 2.0 * (ll_i - ll_red)
    perm_p = float((1 + np.sum(null >= lr)) / (n_perm + 1))
    return lr, coef_driver, perm_p


def paired_delta_auroc(driver, signed, labels, n_boot=2000, seed=0, alpha=0.05):
    """AUROC of each scorer for the same labels + a paired-bootstrap CI of the
    DIFFERENCE (driver - signed): each bootstrap resamples the positive and negative
    indices ONCE and scores both, so the CI reflects the head-to-head gap."""
    driver = np.asarray(driver, dtype=float)
    signed = np.asarray(signed, dtype=float)
    labels = np.asarray(labels, dtype=bool)
    au_d = auroc(driver, labels)
    au_s = auroc(signed, labels)
    pi = np.where(labels)[0]
    ni = np.where(~labels)[0]
    rng = np.random.default_rng(seed)
    diffs = np.empty(n_boot)
    for b in range(n_boot):
        idx = np.concatenate([rng.choice(pi, len(pi), replace=True),
                              rng.choice(ni, len(ni), replace=True)])
        lab = np.concatenate([np.ones(len(pi), bool), np.zeros(len(ni), bool)])
        diffs[b] = auroc(driver[idx], lab) - auroc(signed[idx], lab)
    ci = (float(np.quantile(diffs, alpha / 2)), float(np.quantile(diffs, 1 - alpha / 2)))
    return au_d, au_s, ci


@dataclass
class Claim2Result:
    n_pos: int
    n_neg: int
    auroc_driver: float
    auroc_signed: float
    delta_auroc: float          # auroc_driver - auroc_signed
    delta_ci: tuple             # paired-bootstrap CI of the difference
    lr_stat: float
    driver_coef: float          # standardized driver coef in label ~ signed + driver
    perm_p: float
    topk_fold_driver: float
    topk_fold_signed: float


def evaluate_claim2(driver, signed, labels, opening_only=True, open_thresh=0.0,
                    n_per_pos=3, n_bins=20, seed=0, n_boot=2000, n_perm=1000, k_frac=0.05):
    """Head-to-head driver_score vs signed_delta on a matched, opening-only sample.

    `signed` is the signed-Delta accessibility per region (its magnitude is the
    matching confound; its value is also the baseline scorer). Restricts to regions
    with finite signed and (if opening_only) signed > open_thresh, matches the negative
    background on |signed|, then compares both scorers on positives + matched negatives.
    """
    driver = np.asarray(driver, dtype=float)
    signed = np.asarray(signed, dtype=float)
    labels = np.asarray(labels, dtype=bool)

    keep = np.isfinite(signed) & np.isfinite(driver)
    if opening_only:
        keep = keep & (signed > open_thresh)
    driver, signed, labels = driver[keep], signed[keep], labels[keep]

    confound = np.abs(signed)
    pos = np.where(labels)[0]
    neg = matched_negative_indices(labels, confound, n_per_pos, n_bins, seed)
    if len(pos) == 0 or len(neg) == 0:
        raise SystemExit("no positives or no matched negatives after filtering")

    sample = np.concatenate([pos, neg])
    lab = np.concatenate([np.ones(len(pos), bool), np.zeros(len(neg), bool)])
    d_s, s_s = driver[sample], signed[sample]

    au_d, au_s, dci = paired_delta_auroc(d_s, s_s, lab, n_boot=n_boot, seed=seed)
    lr, coef, perm_p = incremental_lr_test(d_s, s_s, lab, n_perm=n_perm, seed=seed)
    # top-k fold uses the FULL opening-only ranking (not just the matched sample)
    fold_d = topk_fold_enrichment(driver, labels, k_frac=k_frac)
    fold_s = topk_fold_enrichment(signed, labels, k_frac=k_frac)

    return Claim2Result(
        n_pos=len(pos), n_neg=len(neg), auroc_driver=au_d, auroc_signed=au_s,
        delta_auroc=au_d - au_s, delta_ci=dci, lr_stat=lr, driver_coef=coef,
        perm_p=perm_p, topk_fold_driver=fold_d, topk_fold_signed=fold_s)


def _load_signed_track(path, chrom, start, end):
    """Join a signed-Delta track (chrom,start,end,value TSV) onto the contract regions
    by exact coordinate key; missing -> NaN."""
    cc, cs, ce, cv, _ = load_contract(path)   # reuse loader: value in col 4
    key = {(a, int(b), int(d)): v for a, b, d, v in zip(cc, cs, ce, cv)}
    return np.array([key.get((a, int(b), int(d)), np.nan)
                     for a, b, d in zip(chrom, start, end)])


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--contract", required=True, help="driver-score contract TSV")
    ap.add_argument("--positives", required=True, help="driver-label BED (master-TF loci/binding)")
    ap.add_argument("--signed", default=None,
                    help="signed-Delta track TSV chrom,start,end,value (the baseline "
                         "scorer AND the matching confound); if omitted, the contract's "
                         "`direction` column is used")
    ap.add_argument("--all-regions", action="store_true",
                    help="do NOT restrict to opening regions (default is opening-only, "
                         "matching Claim 1's directional test)")
    ap.add_argument("--open-thresh", type=float, default=0.0)
    ap.add_argument("--n-per-pos", type=int, default=3)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--k-frac", type=float, default=0.05)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    chrom, start, end, score, direction = load_contract(args.contract)
    pc, ps, pe = load_bed(args.positives)
    labels = overlap_labels(chrom, start, end, pc, ps, pe)
    signed = (_load_signed_track(args.signed, chrom, start, end)
              if args.signed else direction)
    if not np.isfinite(signed).any():
        raise SystemExit("no finite signed-Delta values; pass --signed with the track")

    res = evaluate_claim2(score, signed, labels, opening_only=not args.all_regions,
                          open_thresh=args.open_thresh, n_per_pos=args.n_per_pos,
                          n_bins=args.n_bins, seed=args.seed, k_frac=args.k_frac)

    print(f"[{'all-regions' if args.all_regions else 'opening-only'}]  "
          f"positives={res.n_pos}  matched negatives={res.n_neg}")
    print(f"AUROC driver_score : {res.auroc_driver:.3f}")
    print(f"AUROC signed-Delta : {res.auroc_signed:.3f}   (baseline)")
    print(f"Delta AUROC        : {res.delta_auroc:+.3f}   "
          f"CI[{res.delta_ci[0]:+.3f}, {res.delta_ci[1]:+.3f}]  "
          f"({'driver beats baseline' if res.delta_ci[0] > 0 else 'not distinguishable / baseline wins'})")
    print(f"incremental LR     : {res.lr_stat:.2f}  driver_coef={res.driver_coef:+.3f}  "
          f"perm p={res.perm_p:.4g}  "
          f"({'adds beyond signed-Delta' if res.perm_p < 0.05 else 'no incremental value'})")
    print(f"top-{args.k_frac:.0%} fold      : driver {res.topk_fold_driver:.2f}x   "
          f"signed {res.topk_fold_signed:.2f}x")


if __name__ == "__main__":
    main()
