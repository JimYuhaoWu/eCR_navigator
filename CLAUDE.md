# CLAUDE.md — eCR_navigator

## Coding principles

- **No features beyond what was asked.** No speculative abstractions, configurability, or error handling for impossible scenarios.
- **Surgical changes.** Touch only what the task requires. Don't improve adjacent code, comments, or formatting. Match existing style.
- **Surface tradeoffs before coding.** If multiple interpretations exist, present them — don't pick silently. If something is unclear, ask.
- **If you notice unrelated dead code or issues, mention them — don't silently fix them.**

## What this project is

Driver-vs-passenger **genomic region weighting** for cell-fate transitions, and
the target-nomination front-end of the **eCR** platform. It produces per-region
**driver scores** that (a) weight off-target severity in `eCR_predictor` (Tier 2)
and (b) nominate target loci for engineered chromatin regulators.

## Multi-species requirement (FIRST-CLASS — do not hardcode mouse)

eCR must work for **at least two species: Homo sapiens (hg38) and Mus musculus
(mm10)**. MEF→mESC is only the *current* example; human transdifferentiation is
planned. Therefore:

- **Every model integration must support both species.** When wiring a genomic
  foundation model, check which assemblies/checkpoints it ships and treat species
  as a parameter, not a constant. Do not bake `mm10` into code paths.
- **If a model only supports one species**, that is not a dead end — bridge it
  (e.g. **liftOver** the region set to the supported assembly, run the model,
  lift scores back), and document the bridge and its losses.
- The output contract stays per-run single-assembly (`region_weight_contract.md`):
  the assembly is whatever the peaks/genome for *that* run use. Just make sure the
  pipeline can produce either.
- ChromBERT case: **both hg38 and mm10 are native** (mm10 = 5k regulators, hg38 =
  6k), selected by a `-g` flag — no liftOver needed. The mm10 model+data were
  fetched and verified on the mirror. So liftOver is a fallback reserved for models
  that truly lack a species, not ChromBERT. See `docs/server_mirrors.md`.

## Sibling repos (all under github.com/JimYuhaoWu)

- `eCR_mod_lib` — module library (default branch **main**)
- `eCR_predictor` — prediction/design pipeline (default branch **master**)
- `eCR_navigator` — this repo (default branch **main**)

They sit as siblings; dev on Windows/Mac, test on the Linux HPCC.

## Server mirrors (one per model)

Each genomic foundation model runs in its own GPU mirror. To connect to a new one
(this will recur as more models are added), follow the reusable playbook in
[`docs/mirror_onboarding.md`](docs/mirror_onboarding.md); per-model runtime details
and the snapshot policy live in [`docs/server_mirrors.md`](docs/server_mirrors.md).

## Inputs on hand

- **ATAC differential accessibility** — MEF vs mES MACS peaks, 5-col BED, **mm10**
  (e.g. `MEF.e7_peaks.bed`, `mES.e7_peaks.bed`). The union of these already drives
  eCR_predictor's off-target background; eCR_navigator adds the *importance* layer.
- **RNA differential expression** — MEF vs mES.
- **Public reprogramming time-course** data — to be collected (drivers change
  early/causally, passengers late/downstream — temporal ordering is signal).
- `atac_embed.pkl` — an ATAC embedding of unknown provenance; confirm before use.

## Open design decision — DECIDE THIS FIRST (new session)

How to score driver vs passenger:

- **Supervised on the reprogramming trajectory.** Use ATAC/RNA time-course:
  drivers change *early and causally*, passengers *late*. Needs temporal (or
  perturbation) data with enough resolution to order events.
- **Zero-shot / pretrained genomic LM** (Enformer / Borzoi-style). Predict a
  region's regulatory impact from sequence; no training data, but not
  reprogramming-specific.
- **Hybrid** — LM embeddings as features + a light supervised head on the
  time-course.

Pick based on what time-course data is actually in hand. Everything downstream
(`model.py`) sits behind the stable output contract, so this choice is swappable.

## Output contract (STABLE — eCR_predictor depends on it)

TSV: `chrom, start, end, driver_score` with `driver_score ∈ [0, 1]`, same assembly
as the peaks for that run (mm10 for the mouse work, hg38 for human). See
`docs/region_weight_contract.md` and
`ecr_navigator/weights.py`. eCR_predictor's off-target Tier-2 maps its union
regions onto this table via a `weight_fn(region)`.

## Layout

```
ecr_navigator/
  features.py   # load embedding artifacts, align two states, compute shift  ✔
  model.py      # embedding shift -> driver_score [0,1] (rank|minmax)         ✔
  weights.py    # read/write the region-weight contract TSV                   ✔
  io.py         # load ATAC / RNA / peaks                                     (todo)
navigate.py     # entrypoint: two artifacts -> driver-weight contract TSV     ✔
scripts/        # server-mirror workflow — kept in-repo (mirror isn't persistent)
  embedding_artifact.py          # SHARED .npz writer — every model emits through this
  mirror_env.sh / setup_mirror.sh   # SSH/conda/env + idempotent mirror setup
  <model>_embed_regions.py       # ONE per model: state peaks -> per-region .npz
                                 #   (get_, atac_, chromfound_; ChromBERT via
                                 #    run_chrombert_region_emb.sh + hdf5_to_artifact.py)
  get_regionmotif_matrix.py      # GET-specific input prep (motif matrix)
  chromfound_build_input.py      # ChromFound-specific input prep (h5ad)
  compute_atpm.sh                # aTPM (accessibility channel) from bigWigs
docs/
  region_weight_contract.md      # STABLE output contract (navigator -> predictor)
  embedding_artifact.md          # internal .npz contract (mirror -> navigator)
  server_mirrors.md              # mirror access + per-model runtime notes
  {get,atacformer,chromfound}_pipeline.md   # per-model run + species notes
  mirror_onboarding.md           # reusable playbook for the next model
```

**Pattern (model-agnostic).** Each foundation model runs in its own GPU mirror and
emits one `.npz` embedding artifact per cell state through the *shared*
`scripts/embedding_artifact.py` writer; `navigate.py` diffs two states into the
driver-score contract. Adding a model = a new `<model>_embed_regions.py` (+ any input
prep) that calls the shared writer — **not** a new pipeline, and never a copied save
block. Integrated so far: ChromBERT, GET, ATACformer, ChromFound (see the per-model
`docs/*_pipeline.md`). Native-mm10 (ChromBERT, GET) are primary for mouse; hg38-only
scATAC models (ATACformer, ChromFound) are primary for human and a liftOver-bridged
reference for mouse.

**Driver-score readout — keep BOTH scores (decided 2026-07-03):**
- **Zero-shot** — embedding-shift between cell states. Always available; needs no
  labels. Implemented (`features.py` + `model.py`).
- **Fine-tuned** — the native transdifferentiation recipe (fine-tune on ATAC
  log2FC → regulator attribution / region score). Only when fine-tuning data
  exists (not always). NOT yet built.

Rationale: fine-tuning data is available for some transitions but not others, so
the zero-shot score is the universal fallback and the fine-tuned score is the
sharper signal when trainable. Report both per region; do not collapse them
prematurely. This means the output will carry two scores — extend the region-weight
contract with an optional second column (or a `method` tag) rather than
overwriting `driver_score`. Keep `driver_score` as the primary/default so
eCR_predictor's existing consumption never breaks.
