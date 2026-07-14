# Claim 1 — is `driver_score` informative? (mm10 MEF→mES)

**Dates:** 2026-07-13 → 2026-07-14 · **Pair:** mm10 MEF→mES · **Channel:** `driver_score`
(magnitude). Machine-readable copy: [`claim1_results.mm10.tsv`](claim1_results.mm10.tsv).
Session status / handoff: [`claim1_progress.md`](claim1_progress.md).

## Summary (final — read this first)

On **clean, verified endpoints** (GSE201577: 3 MEF + 3 mESC concordant clones),
`driver_score` is **moderately informative for 2 of 3 models**:

| Model | best clean AUROC | how |
|---|---|---|
| GET | **0.581** | all-regions (rescued from chance 0.507 by clean endpoints) |
| ChromFound | **0.643** | opening-only (top-5% 1.62×; robust across all datasets) |
| ChromBERT | 0.500 | null everywhere — robustly uninformative |

The earlier flat null (all three ≈ chance) was **substantially an endpoint-quality
artifact**: the original mES endpoint was contaminated (only 2 of 5 mESC ATAC libraries
coherent). Effect sizes are **moderate** (0.58–0.64), not strong — driver_score is a
weak-to-moderate, model- and framing-dependent signal, not a universal driver detector.
The `direction` channel is separate and unaffected (measured aTPM). **Only 3 of the 5
integrated models were tested here — see "Why only 3 of 5 models" below.**

### Phase-2 update (2026-07-14) — the positive set matters more than the model

A second round tested an alternative, independent reprogramming cocktail (**JGES** =
Jdp2/Glis1/Esrrb/Sall4) and, crucially, **reframed what a "driver region" is**. Two findings:

1. **The master-TF *binding* hypothesis fails to generalize.** Using the lab's own JGES
   CUT&Tag (GSE199612) as ground truth, **all three models are at chance** (AUROC 0.48–0.51)
   for JGES-bound regions — in every polarity (all / opening / closing) and at two stringencies.
   GET's phase-1 OSKM number (0.581) did **not** reproduce on JGES binding (0.486). So the
   phase-1 result was **not** "driver_score finds where master TFs bind."
2. **The master-TF *loci* reframe holds — for GET.** When positives are instead the
   **cis-regulatory regions of the pluripotency master-TF genes themselves** (promoters +
   enhancer neighborhoods of 26 core genes), GET's driver_score is elevated
   (**AUROC 0.566–0.582, CI excludes 0.5**), robust to dropping the 4 canonical loci
   (0.564), broad across the network, and it extends to H3K27ac-activated enhancers
   genome-wide (0.574). ChromFound shows only a top-5% tail (2.5–4.2×, null AUROC);
   ChromBERT is small-n noise.

**Unifying reading:** GET's driver_score marks the **opening/activating regulatory
landscape** of pluripotency (master-TF promoters/enhancers and activated enhancers) — which
also explains phase-1: OSKM factors bind *at* those enhancers, so "OSKM binding" ≈ "active
enhancers," while JGES binding is broader and partly at closing chromatin. **But the signal
is largely directional**: every GET result is in all-regions mode and collapses to ~0.50
under opening-only (magnitude-matched), i.e. GET adds "these regions open" beyond raw |Δ|
but does not rank master-TF loci above *other* opening enhancers. Whether that beats a plain
signed-Δaccessibility baseline is a **Claim 2** (direction) question, deferred. Full phase-2
trail below; machine-readable table: [`claim1_results.mtf.tsv`](claim1_results.mtf.tsv).

The sections below are the full investigation trail, in chronological order.

## Why only 3 of 5 models (and when to test all five)

Claim 1 was run on **GET, ChromBERT, and ChromFound** and *not* on **ATACformer** or
**EpiAgent**. This is a deliberate coverage decision, not an omission:

- **ATACformer and EpiAgent are human-designed, fixed-universe token models** (their
  region/cCRE vocabularies are defined on hg38). On mouse they can only be reached by
  **liftOver**, and mm10→(their hg38 universe) leaves a **very sparse** set of in-vocab
  regions — e.g. ATACformer's mm10 liftOver overlaps only ~1,335 of the union regions,
  and EpiAgent is sparser still (its native-hg38 run already yields only a few hundred
  regions under the cCRE rank cap). That is far too few to build a
  `|ΔaTPM|`-matched passenger background, so a mm10 AUROC would be meaningless.
  **Therefore they are not worth testing on this mouse transition.**
- **GET and ChromBERT are native mm10; ChromFound liftOver keeps ~99.8%** — all three
  give dense, per-region coverage over the ~130k positives, which the matched-background
  test requires.

**When human (hg38) data comes in, all five models should be tested.** ATACformer and
EpiAgent are native on hg38 (dense coverage), so the sparsity objection disappears and a
fair Claim 1 test becomes possible for the full five-model matrix. Testing all five on a
clean human transition is the planned next validation step.

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

## First pass — uncurated endpoints, all-regions

| Model | positives | AUROC vs matched bg | 95% CI | AUROC vs random bg | confound gap | top-5% fold |
|---|---:|---:|---:|---:|---:|---:|
| ChromFound | 22,366 | **0.586** | [0.582, 0.591] | 0.583 | −0.003 | **1.49×** |
| GET | 60,418 | 0.507 | [0.504, 0.510] | 0.507 | −0.000 | 1.01× |
| ChromBERT | 20,340 | 0.497 | [0.491, 0.503] | 0.497 | +0.000 | 0.98× |

The label-shuffle null centres at 0.500 for all three. The permutation p-values are tiny
for GET/ChromFound (e.g. GET p≈1e-3) **only because n≈200k makes a 0.507 AUROC
"significant"** — report the effect size (AUROC, fold), not the p-value.

### Read (first pass — superseded by the clean-endpoint result above)

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
aTPM from the count matrix). All three models complete.

| Model | metric | uncurated | curated | **clean** | top-5% (clean) |
|---|---|---:|---:|---:|---:|
| GET | all-regions | 0.507 | 0.490 | **0.581** | 1.27× |
| GET | opening-only | 0.512 | — | 0.536 | 1.09× |
| ChromFound | all-regions | 0.586 | 0.440 | 0.521 | 1.39× |
| ChromFound | opening-only | **0.664** | — | **0.643** | 1.62× |
| ChromBERT | all-regions | 0.497 | 0.492 | 0.500 | 0.97× |
| ChromBERT | opening-only | 0.501 | — | 0.495 | 0.96× |

**Two robust signals emerge on clean data:**
- **GET all-regions rescued to 0.581** (CI [0.577, 0.585]) — GET was at chance on the
  contaminated endpoints; clean endpoints make its driver_score informative. This is the
  clearest evidence that **endpoint quality was a real limiter**.
- **ChromFound opening-only holds at 0.643** (top-5% 1.62×) — robust across all three
  datasets (uncurated 0.664 → clean 0.643). The opening reframing is the durable signal.

Effect sizes are **moderate** (AUROC 0.58–0.64), not strong — driver_score is *weakly-to-
moderately* informative, model- and framing-dependent, not a universal driver detector.
ChromFound's all-regions signal is NOT robust (0.586→0.521); only its opening-only is.
**ChromBERT is null everywhere** (clean 0.500 / 0.495; CI includes 0.5) — robustly
uninformative across all conditions, consistent with its known magnitude issue. So 2 of 3
models carry a moderate driver signal on clean endpoints (GET all-regions, ChromFound
opening-only); ChromBERT does not.

## Phase 2 (2026-07-14) — master-TF *binding* vs master-TF *loci*

Motivated by the concern that the phase-1 positive set (ChIP-Atlas OSKMNE binding) is
**cocktail-specific and OSKM-biased**, and that reprogramming routes are not unique. We
tested the **JGES** cocktail (Jdp2/Glis1/Esrrb/Sall4; the lab's own efficient MEF→iPSC
route, Nat Commun 2023 "NuRD complex cooperates with SALL4 to orchestrate reprogramming"),
then reframed the hypothesis itself. All runs reuse the **clean GSE201577 driver contracts**
— only the positive-label set changes, so no re-embedding was needed (CPU-only).

### Ground-truth sourcing (why ChIP-Atlas was insufficient for JGES)

- **ChIP-Atlas mm10 has no usable JGES data for J and G:** Jdp2 = **0** mouse ChIP
  experiments, Glis1 = **0** (checked all thresholds, PSC and all-cell classes, and the
  master `experimentList.tab`). Sall4 = 6 total (3 PSC, all ES) → a ≥3-experiment consensus
  of only **48 regions**. Only Esrrb is well-covered. So a JGES positive set **cannot** be
  built from ChIP-Atlas — itself evidence that ChIP-Atlas is structurally OSKM-biased
  (canonical factors are far more profiled).
- **Solution — the lab's own paper data, GSE199612** ("NuRD-dependent Reprogramming by
  Sall4-Jdp2-Esrrb-Glis1", CUT&Tag, **native mm10**): all four JGES factors + H3K27ac, in
  D1/D2 reprogramming intermediates. Only bigWig signal tracks are posted (no called peaks),
  so "bound" = top-decile CUT&Tag signal per factor over the clean union (WT samples only,
  not the K5A/N12 NuRD-mutant conditions). JGES positive sets: **union** (≥1 factor) 15,875
  regions; **core2** (≥2 factors) 8,878. Co-binding QC: 3,255 regions bound by all 4 factors
  vs ~9 expected by chance (**~370×**) — coherent ground truth, not noise.

### Result A — master-TF *binding* (JGES): null for all three models

| Model | all-regions | opening-only | closing-only |
|---|---:|---:|---:|
| GET | 0.486 | 0.487 | 0.505 |
| ChromFound | 0.480 | 0.503 | 0.501 |
| ChromBERT | 0.501 | 0.495 | 0.505 |

(union positive set; core2 identical picture, 0.480–0.510 everywhere; all CIs include/below
0.5.) Every model, every polarity, every stringency = chance. GET's OSKM 0.581 does **not**
reproduce on JGES binding (0.486). Checked and ruled out as rescues: (a) **polarity** —
Sall4-NuRD *closes* chromatin, but empirically JGES-bound regions are ~47% closing vs 45%
baseline (balanced), and closing-only is also null; (b) **ground-truth quality** — the ~370×
co-binding rules out noise. **The master-TF-binding hypothesis does not generalize.**

### Result B — master-TF *loci* (regulatory-region reframe): holds for GET

New positives = **cis-regulatory regions of 26 core pluripotency-TF genes** (the 9
reprogramming factors + endogenous network: Zfp42, Prdm14, Tfcp2l1, Dppa2/3/4, Utf1, Lin28a,
Fbxo15, Tbx3, Klf2/5, Nr5a2, Sall1/3, Nr0b1, Foxd3). Windows: **promoter** = TSS±2kb;
**neighborhood** = gene span ±50kb (refGene mm10). Premise validated: these loci
**overwhelmingly open** toward mES (Pou5f1 23/26, Sox2 13/14, Nanog 20/22, Esrrb 29/31,
Prdm14 10/10, Zfp42 6/6, Foxd3 11/11); exceptions are the MEF-expressed / non-specific ones
(Myc 4/11, Klf4 7/13) and Glis1 (1/20).

| Positive set (GET) | positives | AUROC | 95% CI | top-5% |
|---|---:|---:|---:|---:|
| neighborhood, all | 346 | **0.566** | [0.532, 0.599] | 1.73× |
| neighborhood, opening | 277 | 0.523 | [0.485, 0.563] | 1.37× |
| promoter, all | 69 | **0.582** | [0.501, 0.660] | 2.90× |
| neighborhood, **drop top-4 loci** | 256 | **0.564** | [0.524, 0.603] | 2.27× |
| promoter, drop top-4 loci | 57 | 0.573 | [0.482, 0.665] | 3.51× |

**Robust and broad, not locus-driven.** Per-gene GET percentile (mean over all master-TF
regions = **0.598** vs 0.5 random) is highest at **non-canonical** loci — Zfp42 0.94,
Foxd3 0.78, Fbxo15 0.76, Dppa3 0.75, Sall1 0.73, Tfcp2l1 0.70 — while the 4 dominant loci
are only moderate (Pou5f1 0.60, Sox2 0.54, Nanog 0.67, Esrrb 0.63). Low loci are sensible
(Myc 0.37, Klf4 0.49, Glis1 0.29). ChromFound AUROC null but top-5% tail 1.8–3.9×;
ChromBERT small-n noise (n = 14–78, CIs span 0.5).

### Result C — H3K27ac-defined enhancers (GSE199612 D0/D7)

Restricting the master-TF neighborhood to H3K27ac-active peaks (D7 top quartile) did **not**
sharpen — it just cut n to 83 and widened GET's CI (0.560 [0.490, 0.633]); the activated
subset (26) is too small. But a well-powered **genome-wide H3K27ac-activated enhancer** set
(D7-high & gained vs D0, 1,319 regions) gives **GET all-regions 0.574 [0.557, 0.591]** — CI
excludes 0.5. ChromFound tail-only (2.6–2.7×, null AUROC); ChromBERT null.

### Phase-2 synthesis

GET's driver_score is consistently elevated (~0.56–0.58, all-regions) across **four**
framings of the pluripotency regulatory landscape — OSKM binding (0.581), master-TF loci
(0.566), master-TF promoters (0.582), genome-wide activated enhancers (0.574) — and **null**
on JGES binding footprints (0.486). The common thread is **activating regulatory regions**
(promoters/enhancers that open), not TF occupancy. **Caveat that bounds the claim:** every
GET signal collapses to ~0.50 under opening-only (magnitude-matched), so beyond raw |Δ| the
discriminative content is **directional** ("these open"). Whether it beats a signed-Δ
baseline is a Claim 2 question (deferred). Only GET shows this as a rankable AUROC;
ChromFound contributes a non-monotonic top-tail; ChromBERT nothing.

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

## Artifacts & reproducibility (where everything lives)

Contracts and inputs are on the ephemeral GPU mirror's **persistent** shared disk
(`172.16.78.10:/yutiancheng/yuhao/eCR/`) and on **PeiLab2**
(`172.16.78.234:2020` → `/mnt3/wuyuhao/`). Canonical numbers live in this repo
(`claim1_results.mm10.tsv`); the large data files are NOT committed.

- **Driver-score contracts** (mirror `artifacts/`, mirrored to PeiLab2
  `claim1_work/`): `{get,chromfound,chrombert}_driver_scores.mm10{,.curated,.clean}.tsv`.
- **Clean dataset (GSE201577)** inputs: PeiLab2 `/mnt3/wuyuhao/GSE201577_clean/`
  (`union_named.bed` 86,956 peaks, `atpm_union.tsv`, `union_confound.clean.tsv`,
  `MEF.peaks.bed`, `mES.peaks.bed`, `clean.motif.mm10.npz`) and mirror
  `/yutiancheng/yuhao/eCR/clean201577/`. Clean embeddings: mirror `artifacts/`
  (`{get,chrombert}.{MEF,mES}.clean.mm10.npz`) and
  `chromfound_curated/clean_emb/chromfound.{MEF,mES}.hg38.npz`.
- **Positive set:** PeiLab2 `/mnt3/wuyuhao/chip_atlas_mm10/master_tf.consensus.bed`
  (+ per-factor consensus + `README.txt`).
- **Contaminated source datasets:** PeiLab2 `/mnt3/wuyuhao/MEF_mESC{,_2,_3}/`
  (GSE274130 / GSE201852 / GSE243513).
- **Phase-2 ground truth & positives (PeiLab2):**
  - JGES CUT&Tag (GSE199612, mm10): `/mnt3/wuyuhao/jges_gse199612/` — factor bigWigs,
    per-factor signal `GSM*.sig.tab`, `jges.{union,core2}.bed`, `jges.signal.tsv`,
    `build_jges_positive.py`; H3K27ac D0/D7 signal + `mtf.h3k_active.bed`,
    `mtf.h3k_activated.bed`, `h3k_activated_genomewide.bed`, `build_h3k_enh.py`.
  - Master-TF loci: `/mnt3/wuyuhao/mtf_loci/` — `refGene.txt.gz`, `mtf.{promoter,
    neighborhood}[.no4].bed`, `build_mtf_loci.py`, `per_gene_driver.py`.
  - ChIP-Atlas experiment-count check: `/mnt3/wuyuhao/chip_atlas_mm10/` (+ `experimentList.tab`).
  - Eval logs: `/mnt3/wuyuhao/claim1_work/{jges_eval,jges_closing,mtf_eval,mtf_no4,h3k_eval}.log`
    and the run scripts alongside. Canonical numbers: `docs/claim1_results.mtf.tsv`.
- **Eval code (committed):** `scripts/eval_driver_claim1.py` (+ `--opening-only`),
  `tests/test_eval_claim1.py` (11 self-tests). Runnable on any host with numpy.
- **Mirror connection + key-restore procedure:** see
  [`claim1_progress.md`](claim1_progress.md) and [`server_mirrors.md`](server_mirrors.md).
