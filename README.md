# eCR_navigator

Part of the **eCR** engineered-chromatin-regulator platform for MEF→iPSC
reprogramming design, alongside
[`eCR_mod_lib`](https://github.com/JimYuhaoWu/eCR_mod_lib) (module library) and
[`eCR_predictor`](https://github.com/JimYuhaoWu/eCR_predictor) (DBD → structure →
fusion → off-target pipeline).

## What this is

A **region-weighting and target-nomination front-end**. During reprogramming,
some genomic regions open/close as **drivers** of the cell-state transition while
others move as **passengers** (downstream consequences). eCR_navigator assigns
each region a **driver-importance score** and can nominate the regions worth
targeting with an engineered chromatin regulator.

## How it connects to the rest of the platform

1. **Off-target severity (Tier 2).** `eCR_predictor`'s off-target module scores a
   DBD's unintended binding across accessible regions as `Σ (site strength ×
   region weight)`. Today it uses a provisional weight from accessibility
   *dynamics* (|Δaccessibility| proxy). eCR_navigator replaces that with a real
   **driver score** — an off-target hit in a driver region is far more dangerous
   than one in a passenger region. It plugs in through a stable region-weight
   contract; no change to eCR_predictor's scoring code.
2. **Target nomination.** Driver regions can be fed as targets into
   `eCR_predictor` Step 1.

## Output contract (what eCR_predictor consumes)

A TSV of scored regions — `chrom  start  end  driver_score` (driver_score in
[0, 1], higher = more driver-like) — on the **same genome assembly as the peaks**
for that run (mm10 for the mouse work, hg38 for human). See
[docs/region_weight_contract.md](docs/region_weight_contract.md). This contract is
the stable interface; the model behind it can evolve freely.

## Status

**Zero-shot driver scoring working across five foundation models (2026-07).** Each
model runs in its own GPU mirror and emits a per-region embedding artifact
([docs/embedding_artifact.md](docs/embedding_artifact.md)); `navigate.py` diffs two
cell states into the driver-score contract — all behind the *same* artifact contract,
so adding a model is a new embed script, not a new pipeline. Full MEF→mES (mm10) runs
are done for all five, and **all five are now validated end-to-end on human (hg38)**
too — kidney vs pancreas (ENCODE): **ChromBERT** 42,305 regions
([docs/chrombert_pipeline.md](docs/chrombert_pipeline.md)), **GET** full 167,488
([docs/get_pipeline.md](docs/get_pipeline.md)), **ChromFound** full 167,488
([docs/chromfound_pipeline.md](docs/chromfound_pipeline.md)), **ATACformer** 40,753
([docs/atacformer_pipeline.md](docs/atacformer_pipeline.md)), and **EpiAgent** 415
(runs clean but sparse — the 8,190-cCRE rank cap,
[docs/epiagent_pipeline.md](docs/epiagent_pipeline.md)). Coverage on human tracks the
model's design: coordinate/motif models (GET, ChromFound) keep every peak; the
fixed-universe token models (ChromBERT grid, ATACformer, EpiAgent) subset to their
vocabulary.

### Model matrix

Integrated (zero-shot embedding-shift → `driver_score ∈ [0,1]`, rank/minmax-normalized):

| Model | Input data | Species (native / bridged) | Zero-shot output | Fine-tuned |
|---|---|---|---|---|
| **ChromBERT** | ATAC peaks → TRN (regulator) context | **hg38 + mm10 both native** (`-g`) | 768-d TRN-embedding shift **+ measured `direction`** (aTPM-Δ, external attach) | ⬜ planned |
| **GET** | 282 Vierstra motif scores **+ aTPM** accessibility | **hg38 + mm10 both native** | 768-d embedding shift **+ input-measured `direction`** (aTPM-Δ, native input channel) | ⬜ planned |
| **ATACformer** | scATAC region **tokens** (fixed 890k-region hg38 universe) | hg38 native · mm10 via liftOver (sparse) | 192-d embedding shift **+ measured `direction`** (aTPM-Δ, external attach) | ⬜ |
| **ChromFound** | scATAC OCRs, coordinate-based **+ continuous accessibility** | hg38 native · mm10 via liftOver (~99.8%) | 128-d embedding shift **+ input-measured `direction`** (accessibility-Δ, native input channel) | ⬜ |
| **EpiAgent** | scATAC **cCRE tokens** (fixed 1.35M-cCRE hg38 universe) | hg38 native · mm10 via liftOver (sparse) | 512-d embedding shift **+ predicted `direction`** (SR head, model-native) | ⬜ |

All five now feed the contract's optional signed `direction` column; the three provenance
tiers (input-measured, predicted-model-native, external-attach) and their trust caveats are
in [`docs/direction.md`](docs/direction.md).

### Validation — is the `driver_score` actually informative? (in progress)

Separate from "does the pipeline run" (above), we test whether `driver_score` recovers
true drivers against a change-magnitude-matched background, on **clean, verified endpoints**
(GSE201577, mm10 MEF→mES). Two rounds:

- **Phase 1 — master-TF ChIP binding (ChIP-Atlas OSKMNE):** GET all-regions AUROC 0.581,
  ChromFound opening-only 0.643, ChromBERT null. (The earlier flat null was substantially an
  endpoint-quality artifact.)
- **Phase 2 — generalization + reframe.** The OSKM binding set is cocktail-biased, so we
  tested an independent route (**JGES** = Jdp2/Glis1/Esrrb/Sall4, from the lab's own
  GSE199612 CUT&Tag) and reframed the target. **(a)** Master-TF *binding* does **not**
  generalize — all three models are at chance on JGES footprints (0.48–0.51); GET's 0.581
  did not reproduce (0.486). **(b)** Master-TF *loci* — the cis-regulatory regions
  (promoter/enhancer) of 26 pluripotency-TF genes — **do** carry elevated GET driver_score
  (0.566–0.582, robust, broad), extending to H3K27ac-activated enhancers genome-wide (0.574).
  The unifying signal is the **opening/activating regulatory landscape**, not TF occupancy —
  but it is **largely directional** (collapses to ~0.50 opening-only), so whether it beats a
  signed-Δaccessibility baseline is a Claim 2 question (deferred).

Only GET shows a rankable AUROC; ChromFound contributes a non-monotonic top-tail; ChromBERT
nothing. **Only 3 of the 5 models are tested on mouse:** ATACformer and EpiAgent are
human-designed fixed-universe models whose mm10 liftOver is too sparse for a matched-background
test — when human (hg38) data comes in, all five should be tested (native assembly, dense
coverage). Full trail, effect sizes, and honest caveats:
[`docs/claim1_results.md`](docs/claim1_results.md) (+ machine-readable
[`docs/claim1_results.mm10.tsv`](docs/claim1_results.mm10.tsv),
[`docs/claim1_results.mtf.tsv`](docs/claim1_results.mtf.tsv)); session handoff / reproduction
paths: [`docs/claim1_progress.md`](docs/claim1_progress.md); cross-model magnitude
consistency: [`docs/cross_model_consistency.md`](docs/cross_model_consistency.md).

- **Phase 3 — HUMAN (hg38), all five models.** The definitive multi-model round on two direct
  reprogramming systems: **fibroblast→iN** (induced neuron, GSE299923/GSE299920; strong clean
  transition) and **fibroblast→iCM** (induced cardiomyocyte, GSE179011; *dropped* — weak
  partial transition). On iN, **only GET is positive** — master-TF **promoters AUROC 0.668**
  (CI excludes 0.5), robust to opening-only — confirming the phase-2 loci reframe on human.
  ChromFound is null on loci but positive on pioneer (Ascl1) binding (0.572); ChromBERT and
  ATACformer are null/below-chance; EpiAgent is too sparse (8,190-cap). The informative
  signal is **model-specific (GET) and needs a clean strong transition**. Full trail:
  [`docs/claim1_human_progress.md`](docs/claim1_human_progress.md) (+ machine-readable
  [`docs/claim1_results.human.tsv`](docs/claim1_results.human.tsv)).

Scoped / candidate (not integrated — see the per-model docs):

| Model | Input data | Species | Output score type | Status |
|---|---|---|---|---|
| **AlphaGenome** | DNA sequence → predicted DNase/ATAC tracks | **hg38 + mm10 both native** | predicted-accessibility Δ **+ signed direction** | scoping — **deprioritized** (its two-state direction is dominated by measured peak-Δ and by EpiAgent's SR head; revisit as an in-silico-mutagenesis / endpoint-generalization track) |
| **scDIFF** | scATAC (cells×peaks) + DNA seq + histone marks | mouse-native (seq generic) | cell-type annotation; repurposable = CACNN accessibility-shift | paused (needs single-cell MEF→mES) |
| **Evo2** | raw DNA sequence **only** | any (hg38 + mm10 native) | sequence-intrinsic constraint/importance **prior** (not a two-state score) | scope-only (needs A100) |

### Two-axis roadmap: zero-shot vs fine-tuned

The **Fine-tuned** column above is the planned second axis. The governing constraint is
**endpoint-only**: to design an eCR for a *novel* transition, the navigator has only the
**start and end** states — never the trajectory between them (creating it is the point).
So *inference is always endpoint-only*, and the two axes differ only in **training-time
supervision** (see [CLAUDE.md](CLAUDE.md) "Endpoint-only principle"):

- **Zero-shot** (today, every score) — endpoint embedding-shift, no labels. The default,
  and the *only* option when no driver supervision exists.
- **Fine-tuned** (not yet built) — the same pretrained backbone with **driver
  supervision** added, from either **completed transitions** (validated drivers of e.g.
  MEF→iPSC, trained across a corpus and deployed on new endpoint pairs) or **wet-lab
  feedback** (perturbation = causal driver labels). This is **not** trajectory
  supervision (retired — it would violate endpoint-only). Gated on assembling a
  driver-labeled corpus.

The output contract keeps `driver_score` as the zero-shot default and carries the
fine-tuned score in an optional second column / `method` tag, so eCR_predictor's
consumption never breaks. This README table is the scoreboard for that comparison as it
fills in. (Fine-tunes that need trajectory / single-cell-at-inference data — ChromFound's
shipped cell-type classifier, scDIFF — are out of scope by the endpoint-only principle.)

Species: hg38 + mm10 are first-class; hg38-only models are bridged by liftOver (see the
per-model docs and CLAUDE.md).
