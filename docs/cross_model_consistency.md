# Cross-model consistency of `driver_score` (magnitude) — mm10 MEF→mES

**Date:** 2026-07-13 · **Pair:** mm10 MEF→mES · **Channel:** `driver_score` only
(the signed `direction` channel is analysed separately — see below).

## TL;DR

Three independent foundation models produce **mutually uncorrelated** `driver_score`
rankings on this pair (all pairwise Spearman ≈ 0). There is **no cross-model
consensus on driver magnitude.** This does *not* invalidate the scores — it means a
single "the models agree" argument cannot be used as evidence, and each model must be
validated **per-model against ground truth** (the Claim 1 matched-background
enrichment). The `direction` channel is unaffected: it is measured aTPM and is
cross-model-identical by construction.

## Result table

Pairwise agreement of the rank-normalised `driver_score` contracts. Regions are
paired by exact coordinate where the region sets share a coordinate basis (GET and
ChromFound are both on the peak-union), and by **maximum-overlap interval join** where
they differ (ChromBERT lives on ChromBERT's fixed 1 kb grid). Both scores are
re-ranked *within the paired subset* before correlating, so a difference in the
original rank-norm denominators cannot create or hide agreement. Because rank-norm is
monotonic, Spearman is invariant to it — this measures the raw embedding-shift ranking.

| Pair | Paired regions | Spearman | Top-10% overlap (fold vs chance) |
|---|---:|---:|---:|
| ChromFound vs GET | 87,030 | **+0.052** | 1.20× |
| ChromBERT vs GET | 39,950 | **+0.002** | 1.13× |
| ChromBERT vs ChromFound | 20,696 | **−0.021** | 0.94× |

Concordance decays toward the bulk and is weakly positive only at the extreme top for
the one pair that shares a coordinate basis (ChromFound vs GET): top-1% overlap 1.62×,
top-5% 1.34×, top-10% 1.20×. Shuffled-label control Spearman ≈ 0.000 (−0.0004),
confirming the near-zero values are real, not a broken estimator.

Machine-readable copy: [`cross_model_consistency.mm10.tsv`](cross_model_consistency.mm10.tsv).

### Provenance

| Model | Contract | Regions | Coord basis | Emb dim |
|---|---|---:|---|---:|
| GET | `get_driver_scores.mm10.tsv` | 206,313 | peak-union (native mm10) | 768 |
| ChromFound | `chromfound_driver_scores.mm10.tsv` | 87,615 | peak-union (hg38 liftOver) | 128 |
| ChromBERT | `chrombert_driver_scores.mm10.tsv` | 39,950 | 1 kb grid (native mm10) | 768 |
| ATACformer | `atacformer_driver_scores.mm10.tsv` | 1,335 | — | 192 |

- Contracts staged at `/yutiancheng/yuhao/eCR/artifacts/` on the model mirror
  (172.16.78.10). All rank-norm, `--direction off` for this comparison.
- **ChromBERT mm10 was generated for this analysis** (2026-07-13): native mm10
  `chrombert_make_dataset` → `chrombert_get_region_emb` (`--mask mm10_5k_mask_matrix.tsv`)
  on both MEF and mES peak sets (125,467 / 128,790 peaks), then `navigate.py --norm
  rank`. `navigate.py` intersects the two states on the 1 kb grid → 39,950 shared
  regions.
- **ATACformer is excluded** from the correlation: its liftOver overlap with the union
  is only 1,335 regions — too few to estimate rank agreement.

## Interpretation

Two non-exclusive explanations, with different remedies:

1. **Endpoint heterogeneity (leading hypothesis).** A prior PCA of the replicates
   shows the MEF and mES states do **not** cleanly separate — several mES replicates
   sit inside the MEF cloud. If the pseudobulk endpoints are contaminated, the
   embedding-shift is noise-dominated, and noise does not correlate across models.
   *Test:* curate the endpoints to QC-passed, tightly-clustering replicates only
   (Step 0 of the Claim 1 plan), re-run, and check whether the correlations rise.

2. **Intrinsic model differences.** Each encoder captures a different regulatory axis
   (GET: motif + aTPM; ChromFound: scATAC OCR accessibility; ChromBERT: TRN/cistrome
   context from peak presence), so "how far the embedding moves" is not the same
   quantity across models. If (1) is ruled out and correlations stay ≈ 0, there is no
   universal `driver_score` and each model must be scored on its own merits.

**The arbiter for both is the Claim 1 matched-background enrichment against real driver
labels (master-TF binding).** Cross-model agreement was only ever a sanity check: a
model can rank *its own* true drivers at the top without agreeing with another model's
full ranking. So this result reshapes — but does not defeat — Claim 1: the enrichment
must be run and reported **per model**, not as a pooled/consensus score.

## Consequences

- **Claim 1** is evaluated **per model** (`scripts/eval_driver_claim1.py`), one AUROC +
  matched-background result per model, not a single consensus number.
- **Direction is separate and robust.** The signed `direction` channel is differenced
  measured aTPM (tier-1 input-measured for GET/ChromFound); it is cross-model-identical
  and is **not** what this note is about.
- **Slide caveat:** the "GET ≡ ChromFound consistency" highlight in the July progress
  figures (`fig3_direction_tiers.svg`) refers to the **direction** split only. It must
  **not** be read as `driver_score` agreement, which is ≈ 0 across all model pairs.
