# docs/ index

Map of the docs — start here to find things fast. **Bold = read first for that topic.**

## Architecture & contracts (the stable interfaces)
| Doc | What it is / read when |
|---|---|
| [`run_bundle_contract.md`](run_bundle_contract.md) | **THE unified output** — `manifest.json` + `weights.tsv` + `nominations.tsv`. The boundary object with eCR_predictor; carries the Gate-1/Gate-2 verdict so the predictor never learns which model won. **Read first for anything about outputs.** Fixtures: `examples/run_bundle/`. |
| [`region_weight_contract.md`](region_weight_contract.md) | The dense `chrom,start,end,driver_score[,direction]` table (= the bundle's `weights.tsv`). Unchanged; consumed by `offtarget.py` Tier-2. |
| [`embedding_artifact.md`](embedding_artifact.md) | Internal `.npz` contract (mirror → navigator). Read when adding a model's embed script. |
| [`direction.md`](direction.md) | The optional signed `direction` column + its three provenance tiers (input-measured / predicted / external-attach) and trust caveats. |

## Validation — is `driver_score` any good? (**read the summary first**)
| Doc | What it is / read when |
|---|---|
| [`validation_summary.md`](validation_summary.md) | **CANONICAL current-state summary of Claim 1 + Claim 2 + the nomination policy.** The one doc to read; everything below is deep provenance. |
| [`claim1_results.md`](claim1_results.md) | Claim 1 mouse full trail (phases 1–2: OSKM binding, JGES, master-TF loci reframe, H3K27ac). Deep dive. |
| [`claim1_human_progress.md`](claim1_human_progress.md) | Claim 1 human (all 5 models, iN + dropped iCM) + top-tail section. Deep dive. |
| [`claim1_progress.md`](claim1_progress.md) | Claim 1 session handoff / reproduction paths (historical). |
| [`cross_model_consistency.md`](cross_model_consistency.md) | Why the 5 models' magnitude rankings don't agree (Spearman≈0). |
| [`claim2_plan.md`](claim2_plan.md) | Claim 2 scope: 2A (done); 2B spec; Claim 3 parked. |
| [`claim2b_results.md`](claim2b_results.md) | **Claim 2B results** — measured `direction` is trustworthy for the ED call at **\|direction\| ≥ 0.05** (consistent on 3 admit transitions; Gate-1-reject control shows no threshold). Two-sided source arm is system-dependent. |
| [`claim2_results.md`](claim2_results.md) | Claim 2A results, top-k sweeps, per-model×species confidence, **nomination policy**. Deep dive. |
| [`benchmark_v1_results.md`](benchmark_v1_results.md) | **v1 benchmark scorecard** — 6 transitions × GET (+ChromBERT on mouse). GET wins on both strong (iN, C/EBPα); fails on all weak. Read with `validation_summary.md`. |
| [`benchmark_spec.md`](benchmark_spec.md) | Benchmark *design*: frozen transition panel so future/self-trained models run the same scorecard. v1 as-built + v2 backlog. |
| [`finetune_plan.md`](finetune_plan.md) · [`finetune_results.md`](finetune_results.md) | **Regime-1 head-only fine-tune probe** (does supervision beat zero-shot?). Result: **no** — head < zero-shot on mouse loci, human iN loci, and OSKMN↔JGES binding transfer. The model is best used **zero-shot**. |
| [`limited_data_strategy.md`](limited_data_strategy.md) | **What to do under scarce driver labels** — ≤3-param combiner, better zero-shot readouts, and the two real unlocks (breadth corpus / Regime-3 perturbation labels). |
| `*.tsv` | Machine-readable results: `claim1_results.{mm10,mtf,human}.tsv`, `claim2_results.tsv`, `finetune_results.{mm10,human,mtfbind}.tsv`, `cross_model_consistency.mm10.tsv`. |

## Per-model pipelines (run + species notes, one per model)
[`get_pipeline.md`](get_pipeline.md) · [`chrombert_pipeline.md`](chrombert_pipeline.md) ·
[`chromfound_pipeline.md`](chromfound_pipeline.md) · [`atacformer_pipeline.md`](atacformer_pipeline.md) ·
[`epiagent_pipeline.md`](epiagent_pipeline.md) · [`alphagenome_pipeline.md`](alphagenome_pipeline.md)
(scoping) · [`evo2_scoping.md`](evo2_scoping.md) (scope-only).
AlphaGenome track lists: `alphagenome_mouse_{atac,dnase}_tracks.tsv`.

## Mirror / server operations
| Doc | What it is / read when |
|---|---|
| [`server_mirrors.md`](server_mirrors.md) | **Mirror access (ports, keys, envs) + per-model runtime notes + artifact locations.** Read when connecting to a GPU mirror or PeiLab2. |
| [`mirror_onboarding.md`](mirror_onboarding.md) | Reusable playbook for wiring up a new model's mirror. |
| [`model_runtime_matrix.md`](model_runtime_matrix.md) | Measured peak GPU mem + torch/CUDA/driver per model. |

> Convention: `*_progress.md` are **session handoff / status** docs (historical once the work
> merges) — the current bottom line always lives in `validation_summary.md`. `*_results.md`
> are the durable deep-dive results.
