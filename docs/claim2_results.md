# Claim 2A results — does `driver_score` beat / add over signed-Δaccessibility?

Test + method: [`claim2_plan.md`](claim2_plan.md). Harness: `scripts/eval_driver_claim2.py`
(promotes Claim 1's signed-Δ *confound* to a competing *scorer*; head-to-head AUROC +
paired-bootstrap ΔAUROC CI + incremental logistic LR test). Model = **GET** only (the
sole Claim-1-positive model). All runs **opening-only, seed 0**, on the same contracts /
positives / signed-Δ tracks Claim 1 used (`/mnt3/wuyuhao/claim1_work/`).

Read the numbers as: does `driver_score` add anything **beyond "rank by how much a region
opens"** (signed-Δ)? `ΔAUROC` CI above 0 **and** incremental LR p<0.05 = yes.

## Mouse MEF→mES (GSE201577 clean) — done 2026-07-15

Contract `get_driver_scores.mm10.clean.tsv`, signed-Δ `get_clean_conf.tsv`,
positives `mtf_loci/mtf.{promoter,neighborhood}.bed`.

| Positives | n_pos | AUROC driver | AUROC signed-Δ | ΔAUROC [95% CI] | incr. LR p | top-5% driver / signed | verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| promoter | 51 | 0.579 | 0.507 | **+0.071 [−0.045, +0.188]** | 0.13 | 2.74× / **5.10×** | not distinguishable (underpowered) |
| neighborhood | 277 | 0.523 | 0.500 | +0.023 [−0.026, +0.073] | 0.21 | 1.37× / **2.38×** | null |

**Mouse read:** GET's `driver_score` does **not** significantly beat or add over the
signed-Δ baseline. On promoters it edges it numerically (+0.071) but the paired CI
straddles 0 and the incremental LR is n.s. (p=0.13, only n=51 positives — underpowered);
neighborhood is flatly null. Notably the signed-Δ **top-5%** is *more* master-TF-promoter
enriched than `driver_score`'s (5.1× vs 2.7×) — i.e. for the top nominations, "rank by
opening" beats the model here.

Caveat: opening-only already conditions on opening, so within opening regions signed-Δ ≈
*opening magnitude*, whose AUROC ≈0.50 means it barely re-ranks promoters — yet
`driver_score` does not clearly do better either. This is the honest, expected outcome
flagged in the plan: on this endpoint pair the model's contribution is **not
distinguishable from a signed-Δ ranking** for master-TF-locus recovery.

## Human iN (fibroblast→induced neuron) — PENDING

The better-powered test (65–85 promoter positives in Claim 1, GET AUROC 0.668). Files:
neural_gt positives + the iN GET contract (`get_driver_scores.in.hg38.tsv`) + the iN
signed-Δ track. To be run and appended here.

## Bottom line so far

On mouse, **`driver_score` does not add over signed-Δ** for master-TF-locus recovery —
consistent with Claim 1's finding that the signal is "largely directional." Whether the
stronger human iN transition (where GET's Claim-1 AUROC was highest, 0.668) shows an
*incremental* edge is the deciding test, pending its inputs.
