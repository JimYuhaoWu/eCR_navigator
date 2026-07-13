# Claim 1 — is `driver_score` informative? (mm10 MEF→mES, first pass)

**Date:** 2026-07-13 · **Pair:** mm10 MEF→mES · **Channel:** `driver_score` (magnitude).
Machine-readable copy: [`claim1_results.mm10.tsv`](claim1_results.mm10.tsv).

## Question

Does a model's `driver_score` rank master-TF-bound regions above other regions, **beyond**
what the change-magnitude `|ΔaTPM|` alone would give? (See `scripts/eval_driver_claim1.py`
and its self-tests `tests/test_eval_claim1.py`.)

## Ground truth

- **Positive set (drivers):** ChIP-Atlas mm10, Pluripotent-stem-cell class, threshold
  Q<1E-20, six pluripotency master TFs — **Pou5f1(Oct4), Sox2, Nanog, Esrrb, Klf4, Myc**.
  Per factor, peaks reproducible across **≥3 independent experiments** were merged to a
  consensus; the six consensus sets were unioned → **132,766 master-TF regions** (46 Mbp).
  A union region is a positive if it overlaps this set. Downloaded to PeiLab2
  `/mnt3/wuyuhao/chip_atlas_mm10/`.
- **Confound (matched):** `|ΔaTPM|` = |mES − MEF| ATAC signal per region, mapped onto each
  contract's regions from `MEF_mESC_3/qc/{MEF,mES}.SignalOnPeaks.txt` via `bedtools map`.
- **Test pair endpoints are the SAME (uncurated) pseudobulks the driver contracts were
  built from** — see caveats.

## Result

| Model | positives | AUROC vs matched bg | 95% CI | AUROC vs random bg | confound gap | top-5% fold |
|---|---:|---:|---:|---:|---:|---:|
| ChromFound | 22,366 | **0.586** | [0.582, 0.591] | 0.583 | −0.003 | **1.49×** |
| GET | 60,418 | 0.507 | [0.504, 0.510] | 0.507 | −0.000 | 1.01× |
| ChromBERT | 20,340 | 0.497 | [0.491, 0.503] | 0.497 | +0.000 | 0.98× |

The label-shuffle null centres at 0.500 for all three. The permutation p-values are tiny
for GET/ChromFound (e.g. GET p≈1e-3) **only because n≈200k makes a 0.507 AUROC
"significant"** — report the effect size (AUROC, fold), not the p-value.

## Read

- **ChromFound is the only model with a real (if modest) signal:** AUROC 0.586, top-5%
  enrichment 1.49×, CI excludes 0.5. Its confound gap ≈ 0, so this is **not** just tracking
  `|ΔaTPM|` — it is driver-specific.
- **GET (0.507) and ChromBERT (0.497) are at chance** on this positive set — their
  `driver_score` does not distinguish master-TF-bound regions from a magnitude-matched
  background.
- Consistent with [`cross_model_consistency.md`](cross_model_consistency.md): if the scores
  do not recover a shared biological ground truth, their mutual disagreement (Spearman ≈ 0)
  is expected.

## Caveats / next steps (do not over-read this first pass)

1. **Uncurated endpoints.** These runs use the existing merged MEF/mES pseudobulks, which
   the replicate PCA shows are heterogeneous (some mES replicates sit in the MEF cloud).
   Endpoint noise degrades the embedding-shift and could wash out driver signal. **Step 0
   (curate endpoints to QC-passed replicates) then re-run** is the primary follow-up; it
   also tests the "endpoint-noise vs intrinsic model-difference" question.
2. **Positive-set breadth.** Positives are *all* mES master-TF binding, including sites that
   are constitutively accessible and do not change in the transition — those are TF-bound
   but not drivers-of-change, and dilute AUROC. A tighter definition (master-TF-bound **and**
   opening toward mES) is worth testing, with care not to re-confound with `|ΔaTPM|`.
3. **Large-n significance.** Ignore the permutation p-values here; use AUROC + CI + fold.
4. Per-model, not pooled — there is no consensus driver_score (see cross-model note).
