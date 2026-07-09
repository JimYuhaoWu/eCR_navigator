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
are done for all five; the human (hg38) path is now under test on kidney/pancreas
peak sets.

### Model matrix

Integrated (zero-shot embedding-shift → `driver_score ∈ [0,1]`, rank/minmax-normalized):

| Model | Input data | Species (native / bridged) | Zero-shot output | Fine-tuned |
|---|---|---|---|---|
| **ChromBERT** | ATAC peaks → TRN (regulator) context | **hg38 + mm10 both native** (`-g`) | 768-d TRN-embedding shift | ⬜ planned |
| **GET** | 282 Vierstra motif scores **+ aTPM** accessibility | **hg38 + mm10 both native** | 768-d embedding shift | ⬜ planned |
| **ATACformer** | scATAC region **tokens** (fixed 890k-region hg38 universe) | hg38 native · mm10 via liftOver (sparse) | 192-d embedding shift | ⬜ |
| **ChromFound** | scATAC OCRs, coordinate-based **+ continuous accessibility** | hg38 native · mm10 via liftOver (~99.8%) | 128-d embedding shift | ⬜ |
| **EpiAgent** | scATAC **cCRE tokens** (fixed 1.35M-cCRE hg38 universe) | hg38 native · mm10 via liftOver (sparse) | 512-d embedding shift | ⬜ |

Scoped / candidate (not integrated — see the per-model docs):

| Model | Input data | Species | Output score type | Status |
|---|---|---|---|---|
| **AlphaGenome** | DNA sequence → predicted DNase/ATAC tracks | **hg38 + mm10 both native** | predicted-accessibility Δ **+ signed direction** | scoping — next to build |
| **scDIFF** | scATAC (cells×peaks) + DNA seq + histone marks | mouse-native (seq generic) | cell-type annotation; repurposable = CACNN accessibility-shift | paused (needs single-cell MEF→mES) |
| **Evo2** | raw DNA sequence **only** | any (hg38 + mm10 native) | sequence-intrinsic constraint/importance **prior** (not a two-state score) | scope-only (needs A100) |

### Two-axis roadmap: zero-shot vs fine-tuned

The **Fine-tuned** column above is the planned second axis. Today every score is
**zero-shot** (an embedding shift between two cell states — needs no labels). Once
transdifferentiation **time-course** data is in hand we will **fine-tune** the
native-transdiff-capable models (ChromBERT/GET first) on ATAC log2FC → regulator/region
attribution and report *both* scores per region. The output contract keeps
`driver_score` as the zero-shot default and carries the fine-tuned score in an optional
second column / `method` tag, so eCR_predictor's consumption never breaks (see
[CLAUDE.md](CLAUDE.md) "keep BOTH scores"). This README table is the scoreboard for that
zero-shot-vs-fine-tuned comparison as it fills in.

Species: hg38 + mm10 are first-class; hg38-only models are bridged by liftOver (see the
per-model docs and CLAUDE.md).
