# Fine-tuning ground-truth corpus — multi-route MEF→mESC labels + augmentation plan

**Purpose.** The fine-tuning bottleneck is *labeled-data size* (a head-only probe overfits at the
current scale — see [`finetune_results.md`](finetune_results.md), [`limited_data_strategy.md`](limited_data_strategy.md)).
This doc records the plan to build a larger, lower-noise **labeled driver corpus** for the
MEF→mESC transition by exploiting **three independent reprogramming routes** plus their
**trajectories (time-course)**, and how we augment the positive set without leaking or
re-teaching the magnitude confound.

Data host: PeiLab2 `/mnt3/wuyuhao/eCR/eCR_predictor/mef_mes_engineer/atac_gt/`.

## 1. The core idea — route-consensus + trajectory labels

Single-transition, prior-knowledge positives are few and noisy. Two orthogonal signals fix this:

- **Cross-route consensus.** A region that changes **concordantly across independent routes**
  (chemical, JGES, OSK) is a high-confidence driver; route disagreement demotes it to passenger.
  Consensus both **de-noises** and **expands** the positive set beyond curated master-TF loci.
- **Trajectory (time-course).** Each route has intermediate timepoints; a driver opens/closes
  *progressively*. `(state_t vs MEF)` for each t is an independent, graded observation of the
  same locus → real Tier-A augmentation (not synthetic noise).

Already computed (`crossroute.log`, `arm_a_results.txt`, `nom_class.tsv`): navigator direction vs
experimental ATAC is **99.8% concordant on OSK (Li 2017)** and **98.6% on JGES** — the consensus
backbone is real and strong.

## 2. Data inventory — three routes (downloaded)

| Route | Series | Samples / timepoints | Processed | Cross-route concordance |
|---|---|---|---|---|
| **Chemical (CiPS)** | GSE110264 (Cao) | CiPS-A/B **D6** (early) + **D40** (late) bigWig; MEF/CiPS endpoints | `cao_chem/` (early_D6, late_D40, to_open/close); `raw/peaks` D6a/b,D40a/b | — |
| **JGES** | GSE199609 (Wang) | J_N12_GES / JGES(K5A) **D1-3**, **D7**; MEF | `wang_jges/` (early_MEF1→late_JGESD7); `raw/peaks` JGESD7 | **98.6%** (352/357) |
| **OSK** | GSE93029 + **Li 2017** (ARM A) | MEF, **mES-OSK-D7**, iPS, ESCs (rep1/2); Li OSKM time-course | `arm_a_results.txt` | **99.8%** (416/417) |
| (endpoints) | GSE201577 | clean MEF / mESC (the navigator bundle source) | `../bundle/` | — |

`GSE199613` downloaded empty (placeholder — recheck or drop).

## 3. Positive-set construction & augmentation (mapped to this data)

Ranked by value (Tier A = real new observations; the only kind that fixes scarcity):

1. **Route-consensus positives (Tier A).** Label a region positive if it changes in the correct
   direction in **≥2 of 3 routes** (tunable: 2/3 vs 3/3 = precision/recall knob). This is the
   primary augmentation — orthogonal evidence, not prior-knowledge re-use.
2. **Trajectory intermediate states (Tier A).** Use every timepoint as its own `(t vs MEF)` shift:
   chemical D6 & D40, JGES D1-3 & D7, OSK D7 + Li time series. Same driver locus, multiple graded
   labels. **Keep all timepoints of one locus in the same CV fold.**
3. **Replicate re-pairing / pseudobulk subsampling (Tier A).** Each state has ≥2 reps; every
   A-rep × B-rep pairing is an independent shift. Highest yield once we add a **single-cell**
   trajectory set (§5) — subsample cells → many pseudobulks per state.
4. **Whole regulatory unit per gene (Tier A).** Include all linked cCREs (promoter + enhancers)
   of each driver gene, not one region — gene-grouped for CV.
5. **Cross-species ortholog transfer (Tier A, later).** Map mouse↔human master-TF loci.
6. **Input perturbations (Tier B):** reverse-complement + ±bp window jitter for GET (recompute
   both endpoints). **Feature-space mixup/jitter (Tier C):** regularization only; jitter σ
   calibrated to measured replicate variance from #3.

## 4. Negative-set — the trap to avoid

Housekeeping-gene negatives are being validated in a separate session. **Warning to carry there:**
housekeeping loci are constitutively open and *don't change* (shift ≈ 0), so a housekeeping-only
negative set collapses the task to "does this region change at all" — which is exactly the
**signed-Δ baseline that already dominates** ([`claim2_results.md`](claim2_results.md)); the head
would win the benchmark while learning nothing beyond magnitude. **We also need passenger
negatives:** regions that change substantially but are *not* route-consensus drivers,
**|Δaccessibility|-matched** to the positives. Housekeeping = easy negatives; passengers = the
hard negatives that force the head to learn driver *identity* at fixed magnitude. The multi-route
data gives passengers for free: changed-in-one-route-only / route-discordant regions.

## 5. Additional data to acquire — trajectory-focused (search 2026-07-29)

Priorities: (a) **dense bulk OSKM time-course** for more intermediate states; (b) **single-cell
(scATAC / multiome) reprogramming trajectories** — these unlock pseudobulk-subsampling
augmentation (#3) and pseudotime intermediates. Accessions from search snippets — **confirm on GEO
before download**:

| Candidate | Type | Why | Accession (verify) |
|---|---|---|---|
| Knaupp 2017 *Cell Stem Cell* "Chromatin Accessibility Dynamics during iPSC Reprogramming" | bulk ATAC, days 0/1/3/5/7 | dense OSKM trajectory | (confirm; maybe GSE-linked to that paper) |
| "Diversification of reprogramming trajectories" *Sci Adv* aba1190 | parallel scRNA + scATAC | single-cell trajectory → pseudobulk aug + pseudotime | (confirm GEO) |
| "TF stoichiometry/motif affinity/syntax … single-cell chromatin dynamics" (PMC10592962) | scATAC, fibroblast→iPSC | recent dense scATAC trajectory | (confirm GEO) |
| MEF reprogramming RNA+ATAC | bulk | extra route/timepoints | GSE213225 (verify) |
| "OSKM factors cooperatively engage chromatin" | OSKM binding/access | early-stage engagement | GSE36570 (verify; may be pre-ATAC assay) |

Not primary (different lineage, keep for cross-transition breadth only): direct cardiac
reprogramming scATAC (PMID 34509499); chemical ESC→totipotent GSE166216.

## 6. Guardrails (non-negotiable)

- **Leakage under leave-one-gene-out:** every augmented copy — timepoint, replicate re-pairing,
  RC/jitter, linked cCRE — stays in the **same gene's fold**. Augment inside folds, never across.
- **Confound matching:** whatever multiplies positives is mirrored on negatives, kept
  **|Δaccessibility|-matched**, or the head just relearns magnitude.
- **Consensus threshold is a hyperparameter** (2/3 vs 3/3), tuned on held-out known genes, not on
  the CV that reports the final number.

## 7. Status / next

- **Done:** 3 routes downloaded + peak-called; cross-route direction concordance computed
  (OSK 99.8%, JGES 98.6%); per-nomination `nom_class.tsv`.
- **Next:** (1) formalize the route-consensus positive set (2/3 & 3/3) + |Δ|-matched passenger
  negatives; (2) fold in trajectory timepoints as graded shifts; (3) confirm + acquire the §5
  single-cell trajectories for replicate/pseudobulk augmentation; (4) re-run the head probe on the
  augmented corpus with strict gene-grouped CV and compare to zero-shot.
