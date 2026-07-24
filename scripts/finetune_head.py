#!/usr/bin/env python
"""
Fine-tuned driver head (Regime 1, head-only probe) — does supervised driver-region
labelling sharpen the score beyond zero-shot driver_score AND the signed-Δ baseline,
and does it generalize to an UNSEEN driver gene within the same transition?

This is the first fine-tuning cut (see docs/finetune_plan.md): a low-capacity head on
FROZEN embeddings, so it reuses the .npz artifacts and cannot overfit a 768-d space on
~50 positives. Endpoint-only by construction — the only features are the two-state
embedding SHIFT and the signed-Δ accessibility, both derivable from the start+end
states alone; no trajectory, no TF identity.

Design (loci target):
  features  = [ |shift| | PCA(emb_end - emb_start, k) direction | signed-Δ ]
              (|shift| ≈ the zero-shot driver_score, so the head starts from it and the
               PCA direction comps can only ADD; --no-magnitude drops it. PCA is
               label-blind, fit once.)
  head      = L2-logistic
  labels    = master-TF loci (positives) vs |signed-Δ|-matched background (negatives)
  validation= LEAVE-ONE-GENE-OUT: a locus of gene g is scored only by a head trained
              without ANY of g's loci; pool out-of-fold scores -> one held-out ranking.
  verdict   = paired ΔAUROC of the head vs zero-shot driver_score and vs signed-Δ, on
              the pooled held-out set. Head must beat BOTH to justify supervision.

opening-only by default (matches the Claim 1/2 directional test): the embedding shift
VECTOR is still full-dimensional, so the head can use embedding direction even though
accessibility direction is fixed.

Pure numpy; self-contained .npz loading (no ecr_navigator import) so it runs on a bare
mirror. Shares primitives with the Claim 1/2 evaluators (upload those two alongside).
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_driver_claim1 import (          # noqa: E402
    auroc, load_contract, matched_negative_indices, overlap_labels,
)
from eval_driver_claim2 import logistic_fit, paired_delta_auroc   # noqa: E402


# --------------------------------------------------------------- artifact + shift
def _load_npz(path):
    z = np.load(path, allow_pickle=False)
    return z["chrom"], z["start"].astype(np.int64), z["end"].astype(np.int64), \
        z["embedding"].astype(np.float64)


def load_shift(npz_a, npz_b):
    """Two state artifacts -> (chrom, start, end, shift[N,D]) over shared regions,
    shift = emb_b - emb_a, aligned by (chrom,start,end) key."""
    ca, sa, ea, xa = _load_npz(npz_a)
    cb, sb, eb, xb = _load_npz(npz_b)
    idx_b = {(c, int(s), int(e)): i for i, (c, s, e) in enumerate(zip(cb, sb, eb))}
    rows_a, rows_b = [], []
    for i, (c, s, e) in enumerate(zip(ca, sa, ea)):
        j = idx_b.get((c, int(s), int(e)))
        if j is not None:
            rows_a.append(i); rows_b.append(j)
    if not rows_a:
        raise SystemExit("no shared regions between the two artifacts")
    ra = np.array(rows_a); rb = np.array(rows_b)
    return ca[ra], sa[ra], ea[ra], xb[rb] - xa[ra]


def join_values(chrom, start, end, vc, vs, ve, vv):
    """Join a (chrom,start,end,value) table onto regions by exact key; missing -> NaN."""
    key = {(c, int(s), int(e)): v for c, s, e, v in zip(vc, vs, ve, vv)}
    return np.array([key.get((c, int(s), int(e)), np.nan)
                     for c, s, e in zip(chrom, start, end)])


def load_named_bed(path):
    """BED with a name in col 4 -> (chrom, start, end, name)."""
    c, s, e, nm = [], [], [], []
    with open(path) as fh:
        for ln in fh:
            ln = ln.rstrip("\n")
            if not ln or ln.startswith(("#", "track", "browser")):
                continue
            f = ln.split("\t")
            c.append(f[0]); s.append(int(f[1])); e.append(int(f[2]))
            nm.append(f[3] if len(f) > 3 else "")
    return np.array(c), np.array(s), np.array(e), np.array(nm)


def gene_of_region(chrom, start, end, gc, gs, ge, gname):
    """For each region, the name of the FIRST overlapping named interval ('' if none).
    Used to tag each positive locus with its driver gene for leave-one-gene-out."""
    out = np.full(len(chrom), "", dtype=object)
    order = {}
    for c in np.unique(gc):
        m = np.where(gc == c)[0]
        order[c] = (gs[m], ge[m], gname[m])
    for i in range(len(chrom)):
        rec = order.get(chrom[i])
        if rec is None:
            continue
        gss, gee, gnn = rec
        hit = np.where((gss < end[i]) & (gee > start[i]))[0]
        if len(hit):
            out[i] = gnn[hit[0]]
    return out


# ------------------------------------------------------------------------- PCA
def pca_transform(shift, k):
    """Label-blind PCA of the shift matrix to k components (centered, via SVD).
    Fit on ALL regions (unsupervised -> no label leakage). Returns Z[N,k]."""
    xc = shift - shift.mean(axis=0, keepdims=True)
    # economy SVD; components are right-singular vectors
    _, _, vt = np.linalg.svd(xc, full_matrices=False)
    comps = vt[:k]
    z = xc @ comps.T
    # standardize columns for a well-conditioned logistic fit
    sd = z.std(axis=0, keepdims=True)
    return z / np.where(sd > 0, sd, 1.0)


# ---------------------------------------------------------- leave-one-gene-out CV
def leave_one_gene_out_scores(feats, labels, genes, confound, n_per_pos=3, n_bins=20,
                              seed=0):
    """Out-of-fold head scores via leave-one-GENE-out CV.

    Positives are grouped by `genes`; negatives (|confound|-matched background) are
    partitioned into as many folds as there are positive genes. Fold g trains a
    logistic head on every OTHER gene's positives + the other negative folds, and
    scores gene g's positives + its held-out negative fold. Returns
    (eval_idx, oof_score, eval_label): the pooled held-out regions, their
    out-of-fold head score, and their label. A region a head trained on is never
    scored by it.
    """
    labels = np.asarray(labels, dtype=bool)
    genes = np.asarray(genes, dtype=object)
    rng = np.random.default_rng(seed)

    pos = np.where(labels)[0]
    neg = matched_negative_indices(labels, confound, n_per_pos, n_bins, seed)
    pos_genes = np.array([genes[i] for i in pos], dtype=object)
    uniq = [g for g in np.unique(pos_genes) if g != ""]
    if len(uniq) < 2:
        raise SystemExit("need >=2 distinct genes for leave-one-gene-out")

    neg_fold = rng.integers(0, len(uniq), size=len(neg))   # random negative partition
    ones = np.ones  # alias

    eval_idx, oof, ev_lab = [], [], []
    for fi, g in enumerate(uniq):
        te_pos = pos[pos_genes == g]
        te_neg = neg[neg_fold == fi]
        tr_pos = pos[pos_genes != g]
        tr_neg = neg[neg_fold != fi]
        if len(te_pos) == 0 or len(tr_pos) == 0 or len(tr_neg) == 0:
            continue
        tr = np.concatenate([tr_pos, tr_neg])
        y = np.concatenate([np.ones(len(tr_pos)), np.zeros(len(tr_neg))])
        X = np.column_stack([ones(len(tr)), feats[tr]])
        beta, _ = logistic_fit(X, y)
        te = np.concatenate([te_pos, te_neg])
        s = np.column_stack([ones(len(te)), feats[te]]) @ beta
        eval_idx.append(te)
        oof.append(s)
        ev_lab.append(np.concatenate([np.ones(len(te_pos), bool),
                                      np.zeros(len(te_neg), bool)]))
    return (np.concatenate(eval_idx), np.concatenate(oof), np.concatenate(ev_lab))


def transfer_scores(feats, labels_train, labels_test, confound, n_per_pos=3, n_bins=20,
                    seed=0):
    """Cross-PANEL transfer: train a head on one driver panel, test on ANOTHER.

    Train on train-panel positives + |confound|-matched negatives; score test-panel
    positives + their matched negatives. Any region used in TRAINING (positive or
    negative) is excluded from the test set, so a site bound by BOTH panels cannot leak
    — the head is only credited for recognizing panel-B drivers it never saw. Returns
    (eval_idx, score, eval_label).
    """
    labels_train = np.asarray(labels_train, dtype=bool)
    labels_test = np.asarray(labels_test, dtype=bool)
    pos_tr = np.where(labels_train)[0]
    neg_tr = matched_negative_indices(labels_train, confound, n_per_pos, n_bins, seed)
    train = np.concatenate([pos_tr, neg_tr])
    y = np.concatenate([np.ones(len(pos_tr)), np.zeros(len(neg_tr))])
    beta, _ = logistic_fit(np.column_stack([np.ones(len(train)), feats[train]]), y)

    used = np.zeros(len(labels_test), dtype=bool)
    used[train] = True                                   # exclude training regions
    pos_te = np.where(labels_test & ~used)[0]
    neg_te = matched_negative_indices(labels_test & ~used, confound, n_per_pos, n_bins,
                                      seed + 1)
    neg_te = neg_te[~used[neg_te]]
    te = np.concatenate([pos_te, neg_te])
    score = np.column_stack([np.ones(len(te)), feats[te]]) @ beta
    ev_lab = np.concatenate([np.ones(len(pos_te), bool), np.zeros(len(neg_te), bool)])
    return te, score, ev_lab


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--emb-a", required=True, help="state-A artifact .npz (e.g. MEF)")
    ap.add_argument("--emb-b", required=True, help="state-B artifact .npz (e.g. mES)")
    ap.add_argument("--contract", required=True, help="zero-shot driver_score contract TSV (baseline 1)")
    ap.add_argument("--signed", required=True, help="signed-Δ track TSV chrom,start,end,value (baseline 2 + confound)")
    ap.add_argument("--loci", required=True,
                    help="positives BED. Leave-one-gene-out needs a gene name in col 4; "
                         "in --test-loci transfer mode this is the TRAIN panel (col 4 ignored)")
    ap.add_argument("--test-loci", default=None,
                    help="if given, TRANSFER mode: train the head on --loci, test on this "
                         "panel (regions shared with training are excluded). No CV.")
    ap.add_argument("--pca-k", type=int, default=15, help="PCA components of the shift")
    ap.add_argument("--no-magnitude", action="store_true",
                    help="exclude the shift L2 norm (the zero-shot signal) from features; "
                         "then the head sees only shift DIRECTION + signed-Δ")
    ap.add_argument("--all-regions", action="store_true", help="do not restrict to opening regions")
    ap.add_argument("--open-thresh", type=float, default=0.0)
    ap.add_argument("--n-per-pos", type=int, default=3)
    ap.add_argument("--n-bins", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    chrom, start, end, shift = load_shift(a.emb_a, a.emb_b)

    cc, cs, ce, cv, _ = load_contract(a.contract)
    driver = join_values(chrom, start, end, cc, cs, ce, cv)
    sc, ss, se, sv, _ = load_contract(a.signed)
    signed = join_values(chrom, start, end, sc, ss, se, sv)

    gc, gs, ge, gname = load_named_bed(a.loci)
    labels = overlap_labels(chrom, start, end, gc, gs, ge)
    genes = (gene_of_region(chrom, start, end, gc, gs, ge, gname)
             if a.test_loci is None else None)     # gene tags only needed for LOGO
    if a.test_loci is not None:
        tc, ts, te_, _ = load_named_bed(a.test_loci)
        labels_test = overlap_labels(chrom, start, end, tc, ts, te_)

    keep = np.isfinite(signed) & np.isfinite(driver)
    if not a.all_regions:
        keep = keep & (signed > a.open_thresh)
    chrom, shift, driver, signed, labels = \
        chrom[keep], shift[keep], driver[keep], signed[keep], labels[keep]
    genes = genes[keep] if genes is not None else None
    if a.test_loci is not None:
        labels_test = labels_test[keep]

    z = pca_transform(shift, a.pca_k)

    def _std(v):
        s = v.std()
        return (v - v.mean()) / (s if s > 0 else 1.0)

    cols = [z, _std(signed)[:, None]]
    if not a.no_magnitude:
        # the shift L2 NORM is (up to normalization) the zero-shot driver_score itself;
        # include it so the head starts from the zero-shot signal and the PCA *direction*
        # comps can only ADD to it. Without this, standardizing the PCA columns strips
        # magnitude and the head can underperform the zero-shot score by construction.
        cols.insert(1, _std(np.linalg.norm(shift, axis=1))[:, None])
    feats = np.column_stack(cols)

    if a.test_loci is not None:
        eval_idx, head_score, ev_lab = transfer_scores(
            feats, labels, labels_test, np.abs(signed), a.n_per_pos, a.n_bins, a.seed)
        mode = f"TRANSFER train={os.path.basename(a.loci)} -> test={os.path.basename(a.test_loci)}"
    else:
        eval_idx, head_score, ev_lab = leave_one_gene_out_scores(
            feats, labels, genes, np.abs(signed), a.n_per_pos, a.n_bins, a.seed)
        n_genes = len(set(g for g in genes[labels] if g != ""))
        mode = f"leave-one-gene-out over {n_genes} genes"
    base_driver = driver[eval_idx]
    base_signed = signed[eval_idx]

    au_head = auroc(head_score, ev_lab)
    _, au_drv, ci_hd = paired_delta_auroc(head_score, base_driver, ev_lab, seed=a.seed)
    _, au_sgn, ci_hs = paired_delta_auroc(head_score, base_signed, ev_lab, seed=a.seed)

    print(f"[{'all' if a.all_regions else 'opening-only'}] pca_k={a.pca_k}  {mode}")
    print(f"held-out positives={int(ev_lab.sum())}  negatives={int((~ev_lab).sum())}")
    print(f"AUROC head (out-of-sample)      : {au_head:.3f}")
    print(f"AUROC zero-shot driver_score    : {au_drv:.3f}")
    print(f"AUROC signed-Δ baseline         : {au_sgn:.3f}")
    def _call(ci):
        if ci[0] > 0:
            return "head wins"
        if ci[1] < 0:
            return "head significantly WORSE"
        return "not distinguishable"
    print(f"Δ head - driver : {au_head - au_drv:+.3f}  CI[{ci_hd[0]:+.3f}, {ci_hd[1]:+.3f}]  ({_call(ci_hd)})")
    print(f"Δ head - signed : {au_head - au_sgn:+.3f}  CI[{ci_hs[0]:+.3f}, {ci_hs[1]:+.3f}]  ({_call(ci_hs)})")
    print("VERDICT: " + ("head beats BOTH baselines — supervision generalizes to unseen genes"
                         if ci_hd[0] > 0 and ci_hs[0] > 0 else
                         "head does NOT clear both baselines — supervision adds nothing transferable"))


if __name__ == "__main__":
    main()
