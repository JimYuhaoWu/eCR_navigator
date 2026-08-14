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
| **OSK** | **GSE93029 = Li 2017** (PMID 29220666; = "ARM A") | downloaded: MEF, mES-OSK-D7, iPS, ESCs (rep1/2). **In-series but NOT yet pulled: OSK-D0/D1/D3/D5/D7 (rep1/2) full trajectory + single/pairwise-factor (O/S/K/OK/OS/SK-D1) + GFP-sorted fractions** | `arm_a_results.txt` | **99.8%** (416/417) |
| (endpoints) | GSE201577 | clean MEF / mESC (the navigator bundle source) | `../bundle/` | — |

**Note — no double-count:** "ARM A / Li 2017" and the OSK route are the *same* series
(GSE93029, GEO superseries "Chromatin Open/Close Logic in Cellular Reprogramming"), not two
datasets. Its **OSK-D0→D1→D3→D5→D7 trajectory is already in-series and un-downloaded** — the
cheapest trajectory win we have. `GSE199613` downloaded empty (placeholder — recheck or drop).

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

**First, and free:** pull the **GSE93029 OSK-D1/D3/D5 intermediates** we already own but skipped
(§2) — that alone gives the OSK trajectory. Then acquire new sets (GEO pages checked 2026-07-29):

| Candidate | Organism / assay | Timepoints | Value / caveat | Accession |
|---|---|---|---|---|
| **Knaupp 2017** *Cell Stem Cell* 21:834 (PMID **29220667**) — "Transient & Permanent Reconfiguration…" | **mouse** MEF→iPSC, bulk ATAC (+ matched ChIP/RNA/BS-seq), **FACS-sorted *successful* intermediates** | MEF, D3/D6 SSEA1+, D9/D12 SSEA1+/cKIT+, iPSC GFP+ (2 reps) | **highest-value NEW set** — a *different* OSKM trajectory that isolates reprogramming-*competent* cells (cleaner drivers than bulk, which averages in dead-end cells) | **GSE101905** (confirmed) |
| **GSE100345** "dynamic changes in expression & accessibility during reprogramming" | **HUMAN** BJ→iPSC, **scRNA + scATAC** (single-cell) | D0/D2/D8/D16 | single-cell → unlocks pseudobulk-subsampling aug (#3) + pseudotime. **Human → cross-species/ortholog use only, NOT the mouse route-consensus** | **GSE100345** (confirmed) |
| GSE213225 "somatic reprogramming of transformed tumorigenic cells" | mouse, bulk ATAC | D0/D3/D6/D9/D12/D15 | dense timepoints **but transformed + genetically-modified starting cells** → off-target state; low priority | **GSE213225** (confirmed) |
| PMC10592962 "TF stoichiometry/motif/syntax … single-cell chromatin dynamics" | fibroblast→iPSC scATAC | — | scATAC trajectory | confirm GEO |
| GSE36570 "OSKM cooperatively engage chromatin" | OSKM occupancy | — | likely **pre-ATAC** (nuclease/ChIP, Soufi-era) → verify assay before use | **GSE36570** (verify assay) |

Not primary (different lineage / off-target, keep only for cross-transition breadth): direct
cardiac reprogramming scATAC (PMID 34509499); chemical ESC→totipotent GSE166216.

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
- **Next:** (1) **pull GSE93029 OSK-D1/D3/D5 intermediates** (already-owned trajectory, free win);
  (2) formalize the route-consensus positive set (2/3 & 3/3) + |Δ|-matched passenger negatives;
  (3) fold in trajectory timepoints as graded shifts; (4) acquire Knaupp-2017
  (**GSE101905**, mouse sorted intermediates) + evaluate GSE100345 (human) for
  cross-species/pseudobulk use; (5) re-run the head
  probe on the augmented corpus with strict gene-grouped CV and compare to zero-shot.
