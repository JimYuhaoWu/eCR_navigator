# Claim 1 validation — session progress / handoff

Living status doc so this work can be resumed after a context reset. Updated 2026-07-14
(clean-endpoint re-run complete for all three models).

## Goal

Test whether the navigator's `driver_score` (embedding-shift magnitude) is **informative**
— i.e. ranks true drivers above passengers — on the mm10 MEF→mES transition, using
master-TF ChIP binding as ground truth. (Claim 2 = direction; not started.)

## Headline findings so far (all committed)

- **No cross-model consensus on driver magnitude.** GET/ChromFound/ChromBERT driver_score
  rankings are mutually uncorrelated (Spearman ≈ 0). See `cross_model_consistency.md`.
- **First-pass enrichment (uncurated endpoints, all-regions):** only ChromFound above
  chance (AUROC 0.586); GET 0.507, ChromBERT 0.497. See `claim1_results.md`.
- **Endpoint curation did NOT sharpen** — GET 0.490, ChromBERT 0.492, and ChromFound
  *flipped* to 0.440. The signal was not robust to endpoint definition.
- **Reframed (opening-only, signed-ΔaTPM-matched) sharpens ChromFound to 0.664
  (top-5% 2.0×)**; GET/ChromBERT stay at chance. This is the most promising result.
- **Root cause identified:** the current mES endpoint is contaminated — of 5 "mESC"
  ATAC libraries only 2 are coherent (SRR29564546, SRR29564555); 2 cluster with MEF
  (SRR29564556/557), 1 is low quality (SRR30151579). Plus a strong batch axis. So the
  driver_score results rest on a weak endpoint.

## Scope: why only 3 of 5 models (mm10)

Claim 1 covers **GET, ChromBERT, ChromFound** only. **ATACformer and EpiAgent are
human-designed, fixed-universe (hg38) token models**; on mouse they are reachable only by
liftOver, which leaves a **very sparse** in-vocabulary region set (ATACformer ~1,335 union
regions; EpiAgent sparser) — too few for a `|ΔaTPM|`-matched background, so a mm10 AUROC is
meaningless. Hence they are **not worth testing on the mouse transition**. GET/ChromBERT are
native mm10 and ChromFound liftOver keeps ~99.8%, so all three give dense coverage.
**When human (hg38) data arrives, ALL FIVE models should be tested** — ATACformer/EpiAgent
are native on hg38, so the sparsity problem disappears and a fair five-model test is
possible. That is the planned next validation step.

## STATUS 2026-07-14 (phase 2): master-TF BINDING vs LOCI — COMPLETE

Second round, driven by the OSKM-bias concern. Reuses the clean driver contracts (positive
labels only change → CPU-only, no re-embedding). Full results in `claim1_results.md`
(Phase 2) + `claim1_results.mtf.tsv`.

- **Master-TF *binding* hypothesis FALSIFIED (generalization).** JGES cocktail
  (Jdp2/Glis1/Esrrb/Sall4) ground truth from the lab's own GSE199612 CUT&Tag (native mm10;
  ChIP-Atlas has 0 mouse ChIP for Jdp2/Glis1, thin Sall4). All three models at chance
  (0.48–0.51), every polarity (all/opening/closing), both stringencies. GET's OSKM 0.581 did
  NOT reproduce (0.486).
- **Master-TF *loci* reframe HOLDS for GET.** Positives = cis-regulatory regions (promoter
  ±2kb + ±50kb neighborhood) of 26 core pluripotency-TF genes → GET AUROC 0.566–0.582 (CI
  excludes 0.5), robust to dropping the 4 canonical loci (0.564), broad across the network
  (per-gene mean pct 0.598; highest at Zfp42/Foxd3/Fbxo15/Dppa3, not the big 4). Extends to
  genome-wide H3K27ac-activated enhancers (GET 0.574). ChromFound top-5% tail only (2.5–4.2×,
  null AUROC); ChromBERT small-n noise.
- **Bound claim:** GET marks the OPENING/ACTIVATING regulatory landscape (explains phase-1:
  OSKM binds at those enhancers). But every GET signal collapses to ~0.50 opening-only →
  largely DIRECTIONAL. Whether it beats a signed-Δ baseline = Claim 2 (deferred).
- **Artifacts:** `/mnt3/wuyuhao/{jges_gse199612,mtf_loci,claim1_work}/` (see
  `claim1_results.md` Artifacts). Claim 2 (direction) still not started.

## STATUS 2026-07-14: clean-endpoint re-run COMPLETE (all 3 models)

All three models re-embedded on GSE201577 clean endpoints and enriched (see
`claim1_results.md`): **GET all-regions rescued to 0.581** (was at chance on contaminated
endpoints), **ChromFound opening-only holds at 0.643** (robust), **ChromBERT null
everywhere** (0.500/0.495). Bottom line: on clean endpoints driver_score is *moderately*
informative for 2 of 3 models (GET all-regions; ChromFound opening-only) — the earlier
flat null was substantially an endpoint-quality artifact. Clean inputs staged at
`/yutiancheng/yuhao/eCR/clean201577/` (shared) and `/mnt3/wuyuhao/GSE201577_clean/`;
clean contracts `*_driver_scores.mm10.clean.tsv` in mirror `artifacts/` +
`/mnt3/wuyuhao/claim1_work/`. Claim 1 clean-endpoint test DONE. Claim 2 (direction) not
started (deprioritized).

## Original clean-endpoint plan on GSE201577

`GSE201577` (ATAC subseries of GSE201578, "RNA-seq and ATAC-seq in MEF and mESC") is a
**clean, verified** matched dataset: 3 MEF clones (C1/C2/C4) + 3 mESC clones (C1/C3/C5),
bulk ATAC. Verified concordant from the count matrix: within-state r 0.96–0.98,
across-state −0.65, every sample coherence margin +1.6, PC1 = 98% = the MEF↔mES axis,
**no outliers, no batch-dominated PC**. This is the first fair test of Claim 1.

Plan (re-embed ALL models from scratch — both endpoints are new, cannot reuse old npz):
1. Build clean inputs from the count matrix (already downloaded, small):
   - union peak set (87,027 peaks from `GSE201577_*_MEFS_mESC.counts.txt.gz`),
   - per-state aTPM (MEF = mean C1/C2/C4, mES = mean C1/C3/C5; CPM → 99th-pct [0,1]),
   - signed ΔaTPM confound.
   - CONFIRM assembly is mm10 before embedding.
2. Download per-sample narrowPeaks (~3.5 MB each) → per-state peak sets for ChromBERT.
3. Re-embed + navigate per model → `*_driver_scores.mm10.clean.tsv`:
   - GET (native mm10): build a NEW 282-motif matrix for the new peaks, then embed
     MEF & mES with the new aTPM. Reuse `get_curated/` scaffolding.
   - ChromBERT (native mm10): per-state peak sets → embed MEF & mES.
   - ChromFound (hg38 liftOver): `chromfound_build_input.py` on new union+aTPM (both
     states) → embed → navigate (hg38) → back-lift hg38→mm10. Reuse `chromfound_curated/`.
4. Enrichment (all-regions AND `--opening-only`) vs `master_tf.consensus.bed`.
5. Record results in `claim1_results.md` + `.tsv`; update this doc.

## Infrastructure (recurring — keys wiped on mirror reboot)

- **PeiLab2** (ground-truth + CPU work): `ssh -p 2020 wuyuhao@172.16.78.234`, key
  `~/.ssh/ecr_navigator` (works). Data under `/mnt3/wuyuhao/`. Has bedtools,
  bigWigAverageOverBed, numpy, internet.
- **GPU mirrors**: `172.16.78.10`, key `~/.ssh/ecr_navigator`, host-key pinning off.
  Ports change per reboot; `authorized_keys` wiped each boot → user must open the
  instance web terminal once to restore the key (see `server_mirrors.md`).
  - ChromBERT: `-p 35963` (env `/opt/conda` base; native mm10; mask
    `~/.cache/chrombert/data/config/mm10_5k_mask_matrix.tsv`).
  - GET / models-zoo: `-p 38524` (conda `get` at
    `/yutiancheng/yuhao/miniconda3`; init `source .../miniconda3/etc/profile.d/conda.sh`;
    ckpt `/yutiancheng/yuhao/models/get/pretrain_fetal_adult/checkpoint-799.pth`;
    repo `/yutiancheng/yuhao/get_model`).
  - ChromFound: `-p 38824` (env `/opt/conda` env `chromfound`, needs
    `pip install pyliftover`; model `/yutiancheng/yuhao/models/chromFound`; repo
    `/root/ChromFound`; chains in `/yutiancheng/yuhao/eCR/refs/`).
  - `/yutiancheng/yuhao/eCR/` is SHARED across mirror instances (persistent data disk).

## Key paths / staged assets

- Navigator code on mirror: `/yutiancheng/yuhao/eCR/nav_run/` (navigate.py + ecr_navigator pkg).
- GET scaffold: `/yutiancheng/yuhao/eCR/get_curated/`. ChromFound scaffold:
  `/yutiancheng/yuhao/eCR/chromfound_curated/`.
- Positive set (ChIP-Atlas mm10 master-TF): `/mnt3/wuyuhao/chip_atlas_mm10/master_tf.consensus.bed`
  (132,766 regions; Pou5f1/Sox2/Nanog/Esrrb/Klf4/Myc, ≥3-experiment consensus, TH20).
- Eval + work dir: `/mnt3/wuyuhao/claim1_work/` (eval_driver_claim1.py, confounds, contracts).
- Clean dataset: `/mnt3/wuyuhao/GSE201577_check/` (count matrices + conc.py).
- Contaminated source datasets: `/mnt3/wuyuhao/MEF_mESC{,_2,_3}/` (GSE274130/201852/243513).

## Eval machinery

`scripts/eval_driver_claim1.py` (pure numpy; 11 self-tests in `tests/test_eval_claim1.py`):
matched-background AUROC + bootstrap CI + shuffle null + top-k fold + confound gap.
`--opening-only` restricts to opening regions (signed ΔaTPM > 0) and matches negatives on
opening magnitude — isolates TF binding from accessibility gain. Report AUROC/CI/fold,
NOT the permutation p (n≈200k makes trivial effects "significant").

## Commits (branch `validation-claim1`)

`86254f1` opening-only + reframed quick look · `f6d4913` curated re-run (negative) ·
`26fac19` curation Step 0 · `99e913e` first-pass enrichment · `46906a9` cross-model
consistency · `c1ee6cb` eval machinery.
