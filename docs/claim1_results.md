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

## Endpoint curation (Step 0) — replicate QC + ChromBERT re-run

Joint clustering of all 13 ATAC libraries (genome-bin `SignalAcrossGenome`, common basis
across folders 1 & 2) pinpoints the contamination **entirely in the mES endpoint**:

- **Coherent mESC = 2 libraries only:** `SRR29564546` + `SRR29564555` (r = 0.90).
- `SRR29564556`, `SRR29564557` (labelled mESC) **cluster with MEF** (across-state corr >
  within-state) — MEF-like, dropped.
- `SRR30151579` (lone folder-1 mESC) is low quality (r ≈ 0.24 to the good pair) — dropped.
- The 8 MEF libraries are broadly consistent.

**Controlled re-run (isolate the mES fix), all three models:** MEF endpoint held fixed
(reused each model's original MEF embedding); mES rebuilt from the 2 coherent libraries
only, re-embedded, re-navigated. ChromBERT from `mES.curated.bed` (150,177 peaks); GET
and ChromFound from curated mES aTPM on the union (`atpm_union.curated_mES.tsv`,
99th-pct [0,1]). Confound = curated `|ΔaTPM|` from the coherent-replicate means.

| Model | uncurated AUROC | **curated-mES AUROC** | 95% CI (curated) | top-5% (curated) |
|---|---:|---:|---:|---:|
| GET | 0.507 | **0.490** | [0.487, 0.492] | 0.99× |
| ChromBERT | 0.497 | **0.492** | [0.487, 0.497] | 0.93× |
| ChromFound | 0.586 | **0.440** | [0.435, 0.445] | 1.39× |

**Curation did not sharpen any model — the opposite.** GET and ChromBERT stay at chance.
**ChromFound's only above-chance signal (0.586) did NOT survive curation — it dropped to
0.440, below chance** (though its top-5% enrichment persists at 1.39×, i.e. the
score↔binding relationship is non-monotonic after curation). Two readings, both honest:
(a) the original 0.586 was partly an artifact of the contaminated endpoint, not robust
driver recovery; and/or (b) the curated mES (only 2 libraries, shallower) is itself a
different, lower-coverage endpoint, so ChromFound is simply **sensitive to endpoint
definition**. Either way, **driver_score informativeness is not robustly established on
this pair** — the result flips sign under a reasonable change of endpoints.

> **This does not support the hypothesis that the models give informative driver
> scores.** Recorded in full per the instruction to report all results honestly. The
> direction channel is unaffected (measured aTPM); this concerns magnitude only.

## Reframed hypothesis — opening-only, signed-matched (quick look, uncurated endpoints)

Positives restricted to master-TF-bound regions that **open toward mES** (signed ΔaTPM > 0);
negatives = non-TF opening regions matched on opening magnitude. This isolates TF binding
from the accessibility gain (`--opening-only` in `eval_driver_claim1.py`; self-tested).

| Model | positives | AUROC (opening, matched) | 95% CI | top-5% fold |
|---|---:|---:|---:|---:|
| ChromFound | 12,322 | **0.664** | [0.658, 0.670] | **2.00×** |
| GET | 40,465 | 0.512 | [0.508, 0.516] | 1.04× |
| ChromBERT | 5,166 | 0.501 | [0.490, 0.512] | 1.03× |

**The reframing sharpens ChromFound** (0.586 all-regions → **0.664** opening-only; top-5%
1.49× → 2.00×): among regions opening by the same amount, TF-bound ones score meaningfully
higher. GET and ChromBERT stay at chance even reframed. Caveats: (1) this is on the
**uncurated** endpoints (which still contain the 2 MEF-like mES libraries), so it needs a
clean dataset to confirm; (2) it is **one model** — not a general property of driver_score.
Encouraging enough to justify obtaining a clean MEF→mESC/iPSC dataset and re-testing.

## Clean-endpoint re-run (GSE201577 — verified concordant: 3 MEF + 3 mESC clones)

Full re-embed of both endpoints on the clean dataset (new 86,956-peak mm10 union, per-state
aTPM from the count matrix). GET and ChromFound done; ChromBERT pending its mirror.

| Model | metric | uncurated | curated | **clean** | top-5% (clean) |
|---|---|---:|---:|---:|---:|
| GET | all-regions | 0.507 | 0.490 | **0.581** | 1.27× |
| GET | opening-only | 0.512 | — | 0.536 | 1.09× |
| ChromFound | all-regions | 0.586 | 0.440 | 0.521 | 1.39× |
| ChromFound | opening-only | **0.664** | — | **0.643** | 1.62× |

**Two robust signals emerge on clean data:**
- **GET all-regions rescued to 0.581** (CI [0.577, 0.585]) — GET was at chance on the
  contaminated endpoints; clean endpoints make its driver_score informative. This is the
  clearest evidence that **endpoint quality was a real limiter**.
- **ChromFound opening-only holds at 0.643** (top-5% 1.62×) — robust across all three
  datasets (uncurated 0.664 → clean 0.643). The opening reframing is the durable signal.

Effect sizes are **moderate** (AUROC 0.58–0.64), not strong — driver_score is *weakly-to-
moderately* informative, model- and framing-dependent, not a universal driver detector.
ChromFound's all-regions signal is NOT robust (0.586→0.521); only its opening-only is.
ChromBERT clean run pending (mirror `:35963` down at run time).

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
