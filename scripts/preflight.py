#!/usr/bin/env python
"""Nomination preflight: which score to trust for a given transition?

The navigator feeds a small top-k of regions to eCR_predictor. Claim 1/2A showed the
choice of scorer is *transition-dependent*: on a strong, clean transition GET's
`driver_score` beats a plain signed-Delta accessibility baseline and its top ~1% is
sharply enriched for master-TF loci (human iN); on a weak/already-directional one it does
not, and signed-Delta is as good or better (mouse MEF->mES). No endpoint-only indicator we
tried (PC1 variance, driver-vs-magnitude rank divergence) separated those two cases — so we
do not try to *predict* "clean enough"; we *measure* it per transition with known biology.

Two gates, run before nominating:

  GATE 1 -- admissibility (endpoint-only, no model, no drivers). Screens the "nothing works"
    failure mode (weak/partial transition, e.g. the dropped iCM system): requires >=2
    replicates per state, within-state correlation above across-state (coherence margin), and
    the START<->END axis to dominate variance (PC1 fraction). This is a reliable REJECT, not a
    reliable ADMIT -- clearing it does NOT mean GET will beat signed-Delta (mouse cleared it
    with the highest PC1 and still lost). Gate 2 is the actual admit decision.

  GATE 2 -- score selection (the decision). Because eCR design always targets a KNOWN cell
    type, that cell type's canonical master-TF loci are known biology, independent of this
    transition. Drop them into the Claim-2A harness as positives and read the *stable*
    statistics (paired Delta-AUROC CI + incremental LR, NOT a few-hit top-1% fold):
      driver beats signed  (Delta-AUROC CI excludes 0 in driver's favour  OR  incremental-LR
      p < 0.05 with positive driver coef)   ->  PRIMARY = GET driver_score top ~1%
      otherwise                              ->  PRIMARY = signed-Delta top-k (GET supplementary)
    In both cases the measured signed-Delta is attached per region for open/close direction.

Thresholds are first-pass, calibrated on n=2 transitions (iN positive, mouse null) -- tighten
as more transitions accrue. Pure numpy; reuses the Claim-1/2A primitives.
"""
from __future__ import annotations

import argparse
import os
import sys
from dataclasses import dataclass

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from eval_driver_claim1 import load_bed, load_contract, overlap_labels  # noqa: E402
from eval_driver_claim2 import evaluate_claim2                          # noqa: E402

# --- Gate-1 thresholds (endpoint admissibility) ---
MIN_REPS = 2            # replicates required per state
MIN_COHERENCE = 0.10    # within-state minus across-state mean correlation
MIN_PC1 = 0.80          # fraction of variance on the START<->END axis
#   calibration: iN PC1~0.89 (admit), mouse MEF->mES ~0.98 (admit), dropped iCM ~0.68 (reject)


# ============================================================ GATE 1: admissibility
@dataclass
class Admissibility:
    n_a: int
    n_b: int
    coherence_margin: float
    pc1_frac: float
    admit: bool
    reasons: list


def admissibility(matrix, is_state_a, min_reps=MIN_REPS, min_coherence=MIN_COHERENCE,
                  min_pc1=MIN_PC1):
    """Endpoint admissibility from a per-replicate accessibility matrix.

    `matrix`: (n_regions, n_samples) raw counts or signal. `is_state_a`: bool per sample
    (True = start state, False = end state). Normalises to log1p CPM, then:
      - coherence_margin = mean within-state sample corr - mean across-state sample corr
      - pc1_frac         = variance on PC1 of the region-centred matrix (the transition axis)
    Admit iff both states have >=min_reps replicates AND coherence_margin>=min_coherence AND
    pc1_frac>=min_pc1.
    """
    m = np.asarray(matrix, dtype=float)
    a = np.asarray(is_state_a, dtype=bool)
    n_a, n_b = int(a.sum()), int((~a).sum())

    # library-size normalise (CPM) then log1p, so correlations aren't depth-driven
    colsum = m.sum(axis=0, keepdims=True)
    colsum[colsum == 0] = 1.0
    logcpm = np.log1p(m / colsum * 1e6)

    # sample-sample correlation; split within-state vs across-state pairs
    c = np.corrcoef(logcpm.T)
    s = len(a)
    within, across = [], []
    for i in range(s):
        for j in range(i + 1, s):
            (within if a[i] == a[j] else across).append(c[i, j])
    coherence = (float(np.mean(within)) - float(np.mean(across))
                 if within and across else float("nan"))

    # PC1 fraction: variance explained by the first component of the region-centred matrix
    centred = logcpm - logcpm.mean(axis=1, keepdims=True)
    sv = np.linalg.svd(centred, full_matrices=False, compute_uv=False)
    var = sv ** 2
    pc1 = float(var[0] / var.sum()) if var.sum() > 0 else float("nan")

    reasons = []
    if n_a < min_reps or n_b < min_reps:
        reasons.append(f"too few replicates ({n_a}/{n_b} < {min_reps})")
    if not (coherence >= min_coherence):
        reasons.append(f"low replicate coherence ({coherence:.3f} < {min_coherence})")
    if not (pc1 >= min_pc1):
        reasons.append(f"transition axis weak (PC1 {pc1:.3f} < {min_pc1})")
    return Admissibility(n_a, n_b, coherence, pc1, not reasons, reasons)


# ============================================================ GATE 2: score selection
@dataclass
class ScoreSelection:
    primary: str            # "GET" or "signed-Delta"
    reason: str
    auroc_driver: float
    auroc_signed: float
    delta_auroc: float
    delta_ci: tuple
    perm_p: float
    driver_coef: float


def select_score(driver, signed, labels, opening_only=True, **kw):
    """Run the Claim-2A comparison on the target-cell master-TF anchors and apply the rule:
    driver is PRIMARY iff its paired Delta-AUROC CI excludes 0 in its favour, OR the
    incremental LR is significant (p<0.05) with a positive driver coefficient."""
    r = evaluate_claim2(driver, signed, labels, opening_only=opening_only, **kw)
    ci_beats = r.delta_ci[0] > 0
    lr_beats = (r.perm_p < 0.05) and (r.driver_coef > 0)
    if ci_beats or lr_beats:
        primary = "GET"
        why = []
        if ci_beats:
            why.append(f"Delta-AUROC CI[{r.delta_ci[0]:+.3f},{r.delta_ci[1]:+.3f}] excludes 0")
        if lr_beats:
            why.append(f"incremental LR p={r.perm_p:.3g} (driver coef {r.driver_coef:+.2f})")
        reason = "driver beats signed-Delta: " + "; ".join(why)
    else:
        primary = "signed-Delta"
        reason = (f"driver does NOT beat signed-Delta "
                  f"(Delta-AUROC {r.delta_auroc:+.3f} CI[{r.delta_ci[0]:+.3f},"
                  f"{r.delta_ci[1]:+.3f}], LR p={r.perm_p:.3g})")
    return ScoreSelection(primary, reason, r.auroc_driver, r.auroc_signed, r.delta_auroc,
                          r.delta_ci, r.perm_p, r.driver_coef)


# ============================================================ matrix loader (Gate 1 input)
def load_matrix(path, state_a, state_b):
    """Read a per-replicate matrix TSV (header row of sample names; leading non-numeric
    metadata columns like featureCounts' Geneid/Chr/Start/End/Strand/Length are skipped
    automatically). Sample columns are selected by substring: those containing `state_a`
    -> state A, those containing `state_b` -> state B. Returns (matrix, is_state_a)."""
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [ln.rstrip("\n").split("\t") for ln in fh if ln.strip()]
    cols_a = [i for i, h in enumerate(header) if state_a in h]
    cols_b = [i for i, h in enumerate(header) if state_b in h]
    if not cols_a or not cols_b:
        raise SystemExit(f"no columns matched '{state_a}' ({len(cols_a)}) / "
                         f"'{state_b}' ({len(cols_b)}) in header")
    use = cols_a + cols_b
    mat = np.array([[float(r[i]) for i in use] for r in rows], dtype=float)
    is_a = np.array([True] * len(cols_a) + [False] * len(cols_b))
    return mat, is_a


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    # Gate 2 (required — the decision)
    ap.add_argument("--contract", required=True, help="GET driver-score contract TSV")
    ap.add_argument("--signed", required=True,
                    help="signed-Delta track TSV chrom,start,end,value (baseline + confound)")
    ap.add_argument("--anchors", required=True,
                    help="target-cell master-TF loci BED (known biology, e.g. promoters)")
    ap.add_argument("--all-regions", action="store_true",
                    help="do NOT restrict to opening regions (default opening-only)")
    # Gate 1 (optional — admissibility; needs a per-replicate matrix)
    ap.add_argument("--matrix", help="per-replicate accessibility matrix TSV (Gate 1)")
    ap.add_argument("--state-a", help="substring identifying start-state sample columns")
    ap.add_argument("--state-b", help="substring identifying end-state sample columns")
    args = ap.parse_args()

    print("=" * 64)
    if args.matrix:
        if not (args.state_a and args.state_b):
            raise SystemExit("--matrix requires --state-a and --state-b")
        mat, is_a = load_matrix(args.matrix, args.state_a, args.state_b)
        adm = admissibility(mat, is_a)
        verdict = "ADMIT" if adm.admit else "REJECT"
        print(f"GATE 1 admissibility : {verdict}")
        print(f"  replicates         : {adm.n_a} (A) / {adm.n_b} (B)")
        print(f"  coherence margin   : {adm.coherence_margin:+.3f}  (>= {MIN_COHERENCE})")
        print(f"  PC1 fraction       : {adm.pc1_frac:.3f}  (>= {MIN_PC1})")
        if not adm.admit:
            print("  reasons            : " + "; ".join(adm.reasons))
            print("  -> transition not admissible; do not trust any score. Get better endpoints.")
    else:
        print("GATE 1 admissibility : SKIPPED (no --matrix); Gate 2 is the decision")

    chrom, start, end, score, direction = load_contract(args.contract)
    pc, ps, pe = load_bed(args.anchors)
    labels = overlap_labels(chrom, start, end, pc, ps, pe)
    sc_c, sc_s, sc_e, sc_v, _ = load_contract(args.signed)
    key = {(a, int(b), int(d)): v for a, b, d, v in zip(sc_c, sc_s, sc_e, sc_v)}
    signed = np.array([key.get((a, int(b), int(d)), np.nan)
                       for a, b, d in zip(chrom, start, end)])

    sel = select_score(score, signed, labels, opening_only=not args.all_regions)
    print(f"GATE 2 score selection")
    print(f"  anchors in-universe: {int(labels.sum())}  "
          f"(warn: <20 is underpowered)" if labels.sum() < 20 else
          f"  anchors in-universe: {int(labels.sum())}")
    print(f"  AUROC driver/signed: {sel.auroc_driver:.3f} / {sel.auroc_signed:.3f}")
    print(f"  {sel.reason}")
    print("=" * 64)
    print(f"PRIMARY = {sel.primary}"
          + ("  (nominate driver_score top ~1%; attach signed-Delta for direction)"
             if sel.primary == "GET" else
             "  (nominate signed-Delta top-k; GET supplementary)"))


if __name__ == "__main__":
    main()
