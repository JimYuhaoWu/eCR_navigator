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

**Claim 1 validation — is `driver_score` INFORMATIVE? (merged, PRs #5/#6; docs
`claim1_results.md` + `claim1_human_progress.md`, machine-readable `claim1_results.*.tsv`).**
Pipeline-runs ≠ informative. Tested by matched-background AUROC (control for
|Δaccessibility|) against master-TF ground truth on clean, strong transitions. Bottom
line: **GET is the one informative model, on the master-TF *loci* reframe** — the
cis-regulatory regions (esp. **promoters**) of the target-cell master-TF genes, NOT where
TFs bind — and only on a clean, strong endpoint pair.

| Model | mouse MEF→mES (clean GSE201577) | human fib→iN (GSE299923) |
|---|---|---|
| **GET** | ✅ master-TF loci **0.57–0.58** (robust, broad) | ✅ master-TF **promoters 0.668**, robust to opening-only |
| ChromFound | loci tail-only (top-5% 2–4×, null AUROC); OSKM opening-only 0.643 | null on loci, **positive on pioneer (Ascl1) binding 0.572** |
| ChromBERT | null | null / below-chance |
| ATACformer | not tested (mm10 liftOver too sparse) | null / below-chance (first fair hg38-native test) |
| EpiAgent | not tested (too sparse) | too sparse (8,190-cCRE rank cap) |

Lessons that shape method choices: (1) master-TF **binding footprints do NOT work** — they
fail to generalize across reprogramming cocktails (mouse OSKM→JGES, human); the informative
target is the master-TF **gene loci** (their promoters/enhancers, which must open). (2) The
signal is **model-specific (GET)** and needs a **strong transition** — the weak/partial human
iCM system (GSE179011) was dropped; all models null there. (3) Whether GET's signal beats a
plain signed-Δaccessibility baseline is **Claim 2A** — now **answered** (see below). Do not
overclaim `driver_score` beyond GET-on-a-clean-transition.

**Claim 2A — does `driver_score` add over signed-Δaccessibility? (branch `validation-claim2`;
docs `claim2_plan.md` + `claim2_results.md`, machine-readable `claim2_results.tsv`).** GET only
(sole rankable Claim-1 signal). signed-Δ promoted from Claim-1 *confound* to a *competing
scorer*; head-to-head paired ΔAUROC + an incremental logistic-LR test. **Result is
system-dependent:** on human **iN** (strong clean transition) driver_score **beats and adds
over** signed-Δ — decisively on master-TF **promoters** (opening-only AUROC 0.664 vs 0.494,
ΔAUROC +0.170 CI[+0.081,+0.264], LR p=0.001; the increment *survives* the tougher all-regions
baseline where the head-to-head AUROC washes out). On mouse **MEF→mES** it does **not** — the
opening-only point estimate favors driver but is underpowered (n=51), and all-regions signed-Δ
wins outright (0.633 vs 0.582). So GET's value is **regulatory-region prioritization at matched
magnitude+direction** on strong/clean pairs, and negligible over signed-Δ on
weaker/already-directional ones (use signed-Δ there — don't pay for the model). This *partly
overturns* the "largely directional" read from Claim 1. **2B** (is the `direction` *column*
correct?) stays deferred — circular for GET/ChromFound/ATACformer/ChromBERT (their direction
IS the measured aTPM-Δ); only a real test for prediction-head models (EpiAgent / AlphaGenome).

**Nomination policy — which score to trust per transition (`scripts/preflight.py`,
`claim2_results.md` §Nomination policy).** The scorer choice is transition-dependent and
**no endpoint-only indicator predicts it** (PC1 cleanliness fails — mouse is cleanest yet GET
loses there; rank-divergence from magnitude fails too). So don't *predict* "clean enough" —
*measure* it per transition with two gates. **Gate 1 (admissibility, endpoint-only):** ≥2
reps/state, replicate-coherence margin ≥0.10, PC1 ≥0.80 — a reliable REJECT of the "nothing
works" mode (dropped iCM), NOT a reliable admit. **Gate 2 (decision):** eCR design always
targets a *known* cell type, so its canonical master-TF loci (known biology) go into the
Claim-2A harness as positives; if GET `driver_score` beats signed-Δ (ΔAUROC CI excludes 0 **or**
incremental-LR p<0.05, driver coef>0) → **PRIMARY = GET top ~1%**, else → **PRIMARY = signed-Δ
top-k**; measured signed-Δ always rides along for direction. Validated: human iN → GET; mouse
MEF→mES → signed-Δ (admissible but driver-not-primary — why Gate 1 alone is insufficient).
Top-k confidence is front-loaded (GET human-iN top-1% ~9–10×, ~1.7× by top-10%) and
GET-specific (all other models null/mid-tail/too-sparse at the top, both species). Thresholds
first-pass (n=2 transitions). **`driver_score` stays a magnitude, not a signed call** — pair
each nominated region with its measured signed-Δ so eCR_predictor knows the intended open/close.

Separate **top-tail signal** (see `claim1_human_progress.md`): the **top 5% of `driver_score`
is 2–4× enriched for master-TF loci/promoters** for several models (GET and ChromFound, both
species), *even when the matched AUROC is ~0.5*. This is the metric that matches the
navigator's job (nominate the top-k regions) — but top-5% fold is NOT |Δ|-matched, so much of
it is a change-magnitude effect. Use it for target-nomination intuition; use the matched AUROC
for the scientific claim.

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
