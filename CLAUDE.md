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

## Endpoint-only principle (FIRST-CLASS design constraint)

The navigator exists to design engineered chromatin regulators (eCRs) that *create*
a cell-fate transition. For a **novel** transition, the intermediate trajectory does
not exist — producing it is the entire point. Therefore:

> **At inference, the navigator has only the START state and the END state
> (bulk or reference accessibility ± expression). It never has the trajectory
> between them.**

This is the governing filter for every method choice:

- **Any method that requires trajectory / time-course / single-cell-along-the-
  transition data *at inference* is disqualified.** This retires the old
  "supervised on the reprogramming trajectory (drivers change early, passengers
  late)" idea — temporal ordering is signal we will never have at deploy time. It
  is also why ChromFound's shipped fine-tune (a single-cell cell-type classifier)
  and scDIFF (trajectory HMM) do not fit: they violate endpoint-only.
- **Zero-shot is the correct default, not a fallback.** The embedding-shift between
  two endpoint states is endpoint-only *by construction*; its power is the
  pretrained foundation-model prior (regulatory logic learned from millions of
  cells), which distinguishes upstream/regulatory changes (drivers) from
  downstream/consequential ones (passengers) without the path between them.

### The three regimes (differ only by *training*-time supervision; inference is always endpoint-only)

1. **Known completed transition exists → fine-tune, deploy on a new path**
   (cross-transition transfer). The asset a completed transition (e.g. MEF→iPSC via
   OSKM) provides is the **validated driver labels** (the master TFs and their target
   regions), not the trajectory per se. Train across *many* completed transitions to
   learn general driver principles; transfer is strongest for biologically related
   paths. Needs a driver-labeled corpus (the open data-assembly problem).
2. **No known path → zero-shot** (the common case for novel eCR design). Pretrained
   prior + endpoint shift. Implemented (`features.py` + `model.py`), validated on
   five models, both species.
3. **Wet-lab feedback → fine-tune (design-build-test-learn loop).** Each experiment
   yields new supervision for *this* system — outcome data (partial trajectory) or,
   best, **perturbation data** (target region X, measure effect = direct causal
   driver labels). Over iterations this converts a regime-2 problem into a regime-1
   one.

Unifying rule: **inference is always endpoint-only; fine-tune iff driver-relevant
supervision exists (regime 1 or 3), else zero-shot (regime 2).** Zero-shot and
fine-tuned are not rivals — the fine-tuned model is the same pretrained backbone with
supervision added; regimes 1/3 build on regime 2's prior. Everything downstream
(`model.py`) sits behind the stable output contract, so the regime is swappable.

## Output contract (STABLE — eCR_predictor depends on it)

The unified deliverable is the **run bundle** (`docs/run_bundle_contract.md`,
`ecr_navigator/nominate.py`) — one directory per run:

- `weights.tsv` — the DENSE region-weight contract, **unchanged**: `chrom, start, end,
  driver_score[, direction]`, `driver_score ∈ [0,1]`, same assembly as that run's peaks.
  eCR_predictor's off-target Tier-2 maps its union regions onto this via `weight_fn(region)`.
  See `docs/region_weight_contract.md` and `ecr_navigator/weights.py`.
- `nominations.tsv` — SPARSE, ranked; the instruction to eCR_predictor's `target.py`
  (region → target site). **Order by `rank`, never threshold on `nomination_score`** — it's
  a percentile, so a top-1% band spans 0.99–1.0 by construction.
- `manifest.json` — the Gate-1/Gate-2 verdict + provenance.

Two invariants: **the navigator resolves GET-vs-signed-Δ** (the predictor never learns which
model won, so a future self-trained model changes nothing downstream), and **refusal is a
first-class output** (a Gate-1 reject emits a valid bundle with zero nominations and a
reason — the correct answer for 3 of the 6 v1 benchmark transitions).

## Layout

```
ecr_navigator/
  features.py   # load embedding artifacts, align two states, compute shift  ✔
  model.py      # embedding shift -> driver_score [0,1] (rank|minmax)         ✔
  weights.py    # read/write the region-weight contract TSV                   ✔
  nominate.py   # nomination policy -> the run bundle (docs/run_bundle_contract.md) ✔
  io.py         # load ATAC / RNA / peaks                                     (todo)
navigate.py     # entrypoint: two artifacts -> --out contract TSV, and/or     ✔
                #   --bundle (runs both gates + nomination -> run bundle)
scripts/        # server-mirror workflow — kept in-repo (mirror isn't persistent)
  embedding_artifact.py          # SHARED .npz writer — every model emits through this
  mirror_env.sh / setup_mirror.sh   # SSH/conda/env + idempotent mirror setup
  <model>_embed_regions.py       # ONE per model: state peaks -> per-region .npz
                                 #   (get_, atac_, chromfound_; ChromBERT via
                                 #    run_chrombert_region_emb.sh + hdf5_to_artifact.py)
  get_regionmotif_matrix.py      # GET-specific input prep (motif matrix)
  chromfound_build_input.py      # ChromFound-specific input prep (h5ad)
  compute_atpm.sh                # aTPM (accessibility channel) from bigWigs
  attach_measured_signal.py      # attach MEASURED accessibility (e.g. aTPM) as the
                                 #   artifact `signal` -> direction, for models with no
                                 #   accessibility head (ATACformer/GET/ChromBERT/ChromFound)
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
block. Integrated so far: ChromBERT, GET, ATACformer, ChromFound, EpiAgent (see the per-model
`docs/*_pipeline.md`). Native-mm10 (ChromBERT, GET) are primary for mouse; hg38-only
scATAC models (ATACformer, ChromFound, EpiAgent) are primary for human and a
liftOver-bridged reference for mouse. All five are validated end-to-end on both mm10
(MEF→mES) and hg38 (kidney vs pancreas).

**Validation — is `driver_score` any good? → see [`docs/validation_summary.md`](docs/validation_summary.md)
(canonical) and the `docs/README.md` index.** Bottom lines to keep in mind here:
- **Claim 1 (informative?):** YES for **GET only**, on the master-TF **loci** reframe
  (promoters/enhancers of the target-cell master-TF *genes*, NOT where TFs bind), and only on a
  **strong clean transition** (mouse loci 0.57–0.58; human iN promoters 0.668). Other models
  null/tail-only. Binding footprints don't generalize across cocktails.
- **Claim 2A (beats signed-Δaccessibility?):** system-dependent — **YES on human iN** (GET
  promoters ΔAUROC +0.170, incremental-LR p=0.001), **NO on mouse** (signed-Δ dominates). GET's
  value is regulatory-region *prioritization*, not direction.
- **Nomination policy (`scripts/preflight.py`):** no endpoint-only indicator predicts which
  score wins (PC1/cleanliness fails — mouse is cleanest yet GET loses). So **Gate 1** rejects
  unusable transitions (≥2 reps, coherence ≥0.10, PC1 ≥0.80; a **coarse reject-only** screen)
  and **Gate 2** *measures* it — run the Claim-2A harness on the known target-cell master-TF
  loci; GET PRIMARY iff it beats signed-Δ, else signed-Δ. Trust **GET's top ~1%** on a strong
  clean transition (front-loaded ~9–10× enrichment); no other model earns a top-k.
  **Gate 1 statistics are computed on a FIXED universe (top 50k most accessible regions) — do
  not change this casually.** PC1/coherence rise as low-signal regions are dropped, so raw
  values aren't comparable across transitions; on its full 1.06M-cCRE universe iN scored 0.792
  and *rejected*, despite being the transition GET demonstrably wins. Thresholds calibrated on
  the n=6 benchmark panel at that fixed universe (`docs/benchmark_v1_results.md`).
- **`driver_score` is a magnitude, not a signed call** — always attach the measured signed-Δ
  for open/close direction. **Claim 2B** (is the direction *column* itself correct?) is deferred
  (circular for all current models; only testable on prediction-head models EpiAgent/AlphaGenome).

**Driver-score readout — keep BOTH scores (decided 2026-07-03; reframed 2026-07-09):**
- **Zero-shot** — embedding-shift between endpoint states. Always available; needs no
  labels; endpoint-only. The default (regime 2 above). Implemented (`features.py` +
  `model.py`).
- **Fine-tuned** — the same pretrained backbone with **driver supervision** added,
  where that supervision comes from **regime 1** (validated drivers of *completed*
  transitions, trained across a corpus, deployed on new endpoint pairs) or **regime 3**
  (wet-lab feedback, ideally perturbation = causal driver labels). NOT trajectory
  supervision (retired — see the endpoint-only principle). Inference stays
  endpoint-only. NOT yet built; gated on assembling a driver-labeled corpus.

Rationale: driver supervision exists for some transitions but not others, so the
zero-shot score is the universal fallback and the fine-tuned score is the sharper
signal when trainable. Report both per region; do not collapse them prematurely. The
output carries two scores — extend the region-weight contract with an optional second
column (or a `method` tag) rather than overwriting `driver_score`. Keep `driver_score`
as the primary/default so eCR_predictor's existing consumption never breaks.

**Do NOT pursue** fine-tunes that need trajectory / single-cell-along-the-transition
data at inference (ChromFound's shipped cell-type classifier, scDIFF) — they violate
the endpoint-only principle.
