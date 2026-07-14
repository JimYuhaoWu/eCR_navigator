# Claim 1 validation — HUMAN (hg38) progress / handoff

Living status doc for the human Claim-1 test (is `driver_score` informative?), the
definitive multi-model round: **all five models are hg38-native**, so this is where
ATACformer and EpiAgent finally get a fair test (mouse liftOver made them too sparse).
Started 2026-07-14. Mirrors the mouse method (`claim1_results.md`): matched-background
AUROC over master-TF **loci** (regulatory regions of the target-cell master TFs) and
**binding** (ChIP/CUT&Tag), all/opening/closing polarities.

## Two direct-reprogramming systems chosen

1. **iN — fibroblast → induced neuron** (Ascl1/E12a). Data: GSE299923 (ATAC, hg19) +
   GSE299920 (in-study Ascl1 ChIP, hg19). **PRIMARY human test bed.**
2. **iCM — fibroblast → induced cardiomyocyte** (GATA4/MEF2C/TBX5/MESP1/MYOCD).
   Data: GSE179011 (ATAC, hg38). **DROPPED early — see below.**

## DECISION 2026-07-14: iCM system dropped (dataset not good enough)

**The iCM dataset (GSE179011) is too weak a transition for a fair Claim-1 test.** QC over
1.06M hg38 cCREs (BJ fibroblast/5F arm; urine/4F arm is a separate cell and anti-correlates,
so cannot be pooled):
- **PC1 = 68% is the source cell type (fibroblast vs urine), NOT the transition.** The
  fibroblast→iCM axis is only **PC2 = 6.4%**.
- **D0 ↔ D14 correlation stays r ≈ 0.65** — the endpoint shift is small (iCM reprogramming
  is partial/incomplete), vs iN's r ≈ 0 and the clean mouse GSE201577's PC1 = 98%.
- **No biological replicates** (n=1 per condition/timepoint).

Consequences in the results: all models null on AUROC; **GET showed only a suggestive
top-5% tail** on cardiac loci (2.3–3.1×, non-monotonic — the same weak pattern ChromFound
gave on mouse), not a rankable signal. So iCM cannot cleanly confirm or refute Claim 1 —
low power by construction. **We stop iCM here and focus compute on iN** (strong transition).
The iCM artifacts/results are kept for the record but not pursued further.

### iCM results captured before dropping (clean GSE179011 fib/5F arm, 233k cCRE universe)

| Model | loci (neighborhood) AUROC | top-5% | binding AUROC | note |
|---|---:|---:|---:|---|
| GET | 0.492 [0.469, 0.515] | 2.65× | 0.485 | null AUROC, suggestive top-tail |
| ATACformer | 0.521 [0.473, 0.567] | 0.63× | 0.477 | null |
| EpiAgent | 0.58 (n=6) | — | 0.37 (n=5) | too sparse (8190-cap) — inconclusive |

## iN system — PRIMARY (in progress)

- **QC: strong, clean.** PC1 = 89% = the fibroblast→iN transition; D0 ↔ D7 r ≈ 0;
  3 tight replicates/state (within r 0.94–0.96); monotonic D0→D2→D4→D7; no outliers.
  Comparable to the clean mouse GSE201577.
- **Endpoints:** start = D0 ×3 (WT/AAA/QA constructs = fibroblast baselines);
  end = D7 WT ×3 (converted iN). hg19 → bridged to hg38 by liftOver of the cCRE universe
  (99.8% mapped); signal measured in hg19, universe/embedding in hg38.
- **Universe:** 348,786 accessible hg38 cCREs (58% opening) + aTPM + signed-Δ confound
  (`/mnt3/wuyuhao/in_clean/`).
- **Ground truth:** 26-gene neural master-TF loci (ASCL1, POU3F2, MYT1L, NEUROD1/2,
  NEUROG1/2, SOX2, PAX6, OLIG1/2, DLX1/2, FOXG1, …) + **Ascl1 binding** (91,328 cCREs,
  in-study WT-D2 ChIP, input-subtracted) — `/mnt3/wuyuhao/neural_gt/`.

### iN results so far

| Model | loci AUROC | top-5% | binding AUROC | note |
|---|---:|---:|---:|---|
| ATACformer | 0.426 [0.375, 0.476] | 0.28× | 0.473 | null / slightly below chance |
| EpiAgent | 0.57 (n=5) | — | — (n=0 prom) | **too sparse (8190-cap) — inconclusive, both systems** |
| **GET** | **0.668 (promoter) / 0.537 (nbhd)** | **2.1× (prom)** | 0.533 (n=80k) | **POSITIVE — see below** |
| ChromFound | 0.498 (null) | 0.45× | **0.572 (n large)** | split: null on loci, **positive on Ascl1 binding** |
| ChromBERT | 0.448 / 0.392 (below) | 1.5× | 0.519 (n=50k, trivial) | null / below-chance |
| ATACformer | 0.426 (below) | 0.28× | 0.473 | null / below-chance |
| EpiAgent | 0.57 (n=5) | — | — | too sparse |

### Final human synthesis (all five models, iN)

**GET is the clear winner and the definitive positive.** On the strong clean iN transition,
only GET recovers drivers at the master-TF regulatory regions — **promoters AUROC 0.668**
(CI [0.606, 0.727]), neighborhood 0.537, binding 0.533 — robust to opening-only. This is
the human confirmation of the mouse phase-2 reframe (GET marks the cis-regulatory regions of
the target-cell master-TF genes).

The other four:
- **ChromFound** — null on loci (0.498) but **positive on Ascl1 binding (0.572)**: it tracks
  the pioneer factor's occupancy, not the master-TF gene loci (opposite emphasis to GET).
- **ChromBERT** — null/below-chance (loci 0.39–0.45; binding 0.519 trivial), same as mouse.
- **ATACformer** — null/below-chance on both loci and binding (its first fair, dense
  hg38-native test — dense coverage ≠ informative driver_score).
- **EpiAgent** — too sparse to test (8,190-cCRE rank cap → too few overlap the loci).

**Bottom line:** the informative driver_score is **model-specific (GET)** and needs a
**clean, strong transition** (iN, not the weak iCM). A capable model on a strong endpoint
pair recovers master-TF-locus drivers on human; weak endpoints (iCM) or weaker models
(ChromBERT/ATACformer/EpiAgent) do not. Direction vs a signed-Δ baseline remains a Claim 2
question. Full numbers: [`claim1_results.human.tsv`](claim1_results.human.tsv).

**GET iN verdict: clear positive (the headline).** On the clean strong iN transition (330k
regions, 23/24 chroms; chr5 excluded — Altius stream kept hanging), GET's driver_score is
**clearly elevated at neural master-TF promoters: AUROC 0.668 [0.606, 0.727] all-regions,
0.664 opening-only** (top-5% 2.1–2.5×); neighborhood 0.537–0.557 (CI excludes 0.5); Ascl1
binding 0.533 (n=80,483, CI excludes 0.5). **Robust to opening-only** — unlike iCM, the
signal is not merely directional; it survives magnitude+direction matching on a strong,
replicated transition. This **confirms the mouse phase-2 reframe on human**: GET marks the
cis-regulatory regions (esp. promoters) of the target-cell master-TF genes. Contrast the two
non-GET models on the same iN data (ATACformer null, EpiAgent too sparse) and the same GET
model on the weak iCM data (null AUROC): the positive needs *both* a capable model (GET) and
a clean strong transition. GPU health confirmed fine (earlier stall was a duplicate-job
deadlock on GPU init, not hardware; single clean relaunch completed normally).

**EpiAgent verdict (both systems): structurally too sparse.** Even hg38-native with dense
input, the 8,190-cCRE rank cap yields ~3.3–3.5k contract regions → only 4–6 overlap the
master-TF loci → no testable AUROC. Confirms the mouse-era concern persists on its own
assembly. **ATACformer verdict (both systems): null** — dense coverage, but driver_score
does not recover drivers (iCM null; iN null/slightly-below-chance).

## Top-5% (top-tail) enrichment — a separate, operationally relevant signal

Across BOTH species and all three ground-truth types, the **top 5% of `driver_score` is
2–4× enriched for master-TF regulatory regions — often even when the confound-matched AUROC
is ~0.5.** This matters because the navigator's job is to **nominate the top handful of driver
regions**, so top-tail enrichment is arguably the metric that matches the use case better than
global AUROC. Compiled `top5_fold` (from `claim1_results.{mm10,mtf,human}.tsv`):

| Ground truth | model / setting | AUROC (matched) | **top-5% fold** | n_pos |
|---|---|---:|---:|---:|
| **promoters** | GET mouse (JGES-panel) | 0.58 | **2.7–2.9×** | 51–69 |
| promoters | GET mouse, drop-top-4-loci | 0.57 | **3.5×** | 40–57 |
| promoters | **GET human iN** | **0.668** | **2.1–2.5×** | 65–85 |
| promoters | ChromFound mouse (opening) | 0.51 (null) | **3.9×** | 31 |
| promoters | GET human iCM (weak system) | 0.47 (null) | **2.3×** | 164 |
| **loci / neighborhood** | GET mouse | 0.57 | 1.7–2.3× | 256–346 |
| loci | **GET human iN** | 0.54 | **2.1–2.4×** | 327–461 |
| loci | **GET human iCM** (weak system) | **0.49 (null)** | **2.3–3.1×** | 372–778 |
| loci | ChromFound mouse (opening) | 0.50 (null) | 2.2× | 109 |
| activated enhancers (gw) | ChromFound mouse | 0.50 (null) | **2.6–2.7×** | 805–880 |
| activated enhancers | ChromFound mouse H3K27ac (opening) | 0.50 (null) | 3.9× | 41 |
| **TF binding** | ChromFound mouse OSKM (opening) | 0.66 | 2.0× | 12k |
| binding | ChromFound human Ascl1 | 0.57 | 1.8× | (large) |
| binding | ChromFound mouse JGES | 0.48 (null) | 1.4–1.65× | 2–10k |
| binding | GET human/mouse | 0.53/0.58 | 1.1–1.3× | large |

**Pattern.** Top-tail enrichment is (a) **strongest for loci/promoters (2–4×)**, weaker for TF
binding (1.2–2×) — the same "loci beat footprints" ordering seen in the AUROC; and (b) present
for **both GET and ChromFound**, and even for the *weak* iCM system and *null-AUROC* settings —
the very top scores concentrate on master-TF regulatory regions across the board.

**Honest caveat — top-5% fold is NOT confound-matched** (unlike the AUROC, which controls for
|Δaccessibility|; `topk_fold_enrichment` in `eval_driver_claim1.py` runs on the raw full
ranking). So when AUROC≈0.5 but top-5% fold is high, the top-scoring regions are enriched for
the ground truth **largely because master-TF loci are big-change regions and the top
driver_scores are big-change regions** — i.e. much of the top-tail signal is a change-magnitude
effect, not confound-controlled evidence of driver-specific recovery. Two regimes:
- **Both AUROC and top-5% high** (GET on loci/promoters, mouse + human iN): real, model-specific
  driver signal (survives |Δ| matching) *and* concentrated at the top — the strong case.
- **Only top-5% high** (GET iCM cardiac loci 2.3–3.1×; ChromFound activated enhancers 2.6×):
  the top nominations are enriched, but it's mostly a |Δ| effect — operationally useful for
  target nomination, but weaker as evidence that driver_score beats a plain |Δaccessibility|
  ranking. (That comparison is exactly the deferred **Claim 2** control.)

**Takeaway for the platform.** For target *nomination* (pick the top-k regions), several
models — not just GET — put ~2–4× more master-TF regulatory regions in their top 5% than
chance. For the scientific claim that driver_score adds information beyond change magnitude,
the matched AUROC remains the rigorous test, and there GET-on-a-clean-transition is the signal.

## Infrastructure notes

- **Motif matrices (GET only)** must be built on **PeiLab2** (`get_regionmotif_matrix.py`,
  remote `tabix -R` over the Vierstra hg38 archetype file): the Model-Zoo mirror's route to
  the Altius server is throttled ~165× (3.5 KB/s vs PeiLab2's 580 KB/s). Long remote streams
  occasionally fail (BGZF error) → just re-run. ~40–60 min per universe.
- **Model Zoo** (`172.16.78.10:38524`): envs `get`, `atacformer`, `EpiAgent` (no ChromFound
  here). ATACformer script `/yutiancheng/yuhao/get_scripts/atac_embed_regions.py`; EpiAgent
  build `epiagent_build_input.py --no-lift` (hg38 native) — needs bedtools on PATH (append
  `.../envs/get/bin`, don't prepend or it shadows the env python).
- **Per-model confound**: `bedtools map` the cCRE signed-Δ onto each model's contract
  regions (contracts differ per model); `run_{human,neural}_eval.sh` do this + all 3 polarities.
- Eval tool: `scripts/eval_driver_claim1.py` (`--opening-only`; negate confound for closing).
