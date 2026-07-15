# Claim 2A results — does `driver_score` add over signed-Δaccessibility?

Ran 2026-07-15 (harness `scripts/eval_driver_claim2.py`, committed with the plan). Follows
[`claim2_plan.md`](claim2_plan.md). No new embedding compute — re-uses the on-disk GET
driver contracts + signed-Δ tracks built for Claim 1 (PeiLab2
`/mnt3/wuyuhao/{claim1_work,in_clean}/`). Machine-readable: [`claim2_results.tsv`](claim2_results.tsv).

## Question

Claim 1's GET positive (recovers master-TF loci vs a |Δaccessibility|-**magnitude**-matched
background) was "largely directional" — it mostly said *these regions open*. So 2A promotes
**signed-Δaccessibility** from Claim 1's *confound* to a *competing scorer* and asks whether
`driver_score` carries anything **beyond the signed-Δ you already measure from the ATAC data**.
GET only — the sole model with a rankable Claim-1 signal.

Two readouts on one matched sample (positives + |signed-Δ|-matched negatives):
- **Head-to-head ΔAUROC** with a *paired* bootstrap CI (same resampled indices for both scorers).
- **Incremental LR test** — does `driver_score` improve a logistic model that already has
  signed-Δ? (`LR = 2(ll_full − ll_reduced)`, permutation p, standardized driver coef.)

Two modes: **opening-only** (primary — same construction as Claim 1: background matched on
|signed-Δ| *and* restricted to opening regions, which neutralizes signed-Δ's magnitude by
design) and **all-regions** (a *tougher* baseline — signed-Δ keeps its full magnitude power).

## Bottom line

**2A is SUPPORTED on the strong, clean human iN transition, and NOT on mouse (underpowered /
signed-Δ dominates).** On human iN, GET's `driver_score` adds regulatory-region information
**beyond "what opens"** — decisively on **promoters**. This *partly overturns* the pessimistic
"largely directional" expectation from Claim 1: the foundation-model prior is not merely a
proxy for signed accessibility change, at least on the best system.

| System | Positive | Mode | AUROC driver / signed | ΔAUROC [CI] | incremental LR (perm p) | Verdict |
|---|---|---|---|---|---|---|
| **iN** | **promoter** | opening | **0.664 / 0.494** | **+0.170 [+0.081, +0.264]** | **18.6 (p=0.001)** | ✅ **driver beats + adds** |
| iN | neighborhood | opening | 0.557 / 0.499 | +0.057 [+0.012, +0.103] | 10.3 (p=0.002) | ✅ driver beats + adds |
| iN | promoter | all | 0.668 / 0.599 | +0.069 [−0.009, +0.151] | 21.0 (p=0.001) | ⚠️ AUROC wash, **LR still sig** |
| iN | neighborhood | all | 0.537 / 0.555 | −0.018 [−0.057, +0.021] | 4.6 (p=0.023) | ⚠️ AUROC wash, LR sig |
| mouse | promoter | opening | 0.579 / 0.507 | +0.071 [−0.045, +0.188] | 2.2 (p=0.13) | ➖ null (underpowered, n=51) |
| mouse | neighborhood | opening | 0.523 / 0.500 | +0.023 [−0.026, +0.073] | 1.5 (p=0.21) | ➖ null |
| mouse | promoter | all | 0.582 / **0.633** | −0.051 [−0.137, +0.037] | 0.8 (p=0.37) | ❌ signed-Δ wins |

## Reading it

- **Human iN promoters are the clean win.** Opening-only, `driver_score` beats signed-Δ on
  *both* the paired ΔAUROC (CI excludes 0) *and* the incremental LR (p=0.001, driver coef
  +0.67) — driver contributes over and above signed accessibility, not just re-stating it.
- **The increment survives the tougher all-regions baseline.** When signed-Δ is allowed its
  full magnitude power, the *head-to-head* AUROC gap closes to a wash — but the **incremental
  LR stays significant on human iN** (promoter p=0.001, neighborhood p=0.023). So even where
  signed-Δ alone ranks about as well, `driver_score` still adds *orthogonal* regulatory
  information. This is the honest, defensible form of the claim: GET's value on iN is
  **regulatory-region prioritization at matched magnitude+direction**, not just direction.
- **Mouse MEF→mES does not support 2A.** The opening-only point estimate favors driver
  (+0.07 on promoters) but the CI includes 0 and the LR is null (n=51 — underpowered); in
  all-regions signed-Δ **wins outright** (0.633 vs 0.582). Mouse master-TF promoters
  co-occur with the largest openings, so signed-Δ magnitude captures them — the mouse GET
  signal really is "largely directional," consistent with Claim 1's opening-only collapse.
- **top-5% fold agrees for nomination.** On human, `driver_score`'s top-5% is 2.1–2.5×
  enriched vs signed-Δ's 0.8–1.5× — driver nominates a better top-k. On mouse, signed-Δ's
  top-5% (5×) exceeds driver's (2.7–2.9×) — again the mouse-is-directional story.

## What this means for the platform

The question 2A was built to answer — *is direction worth model compute at all?* — resolves
**system-dependent**:
- On a **strong, clean transition (human iN)**, the GET foundation-model prior earns its
  compute: it prioritizes the right master-TF regulatory regions **beyond** the signed
  accessibility change measured directly from ATAC.
- On a **weaker / already-directional transition (mouse MEF→mES here)**, a plain signed-Δ
  baseline is as good or better — don't pay for the model where the ATAC delta already tells
  the story.

Use `driver_score` for target nomination where the endpoint pair is strong and clean; fall
back to signed-Δ where it isn't. This does **not** settle **2B** (is the `direction` *column*
itself correct?) — that stays deferred and is only non-circular for prediction-head models
(EpiAgent / AlphaGenome); GET's direction *is* the measured aTPM-Δ, so testing it against
signed-Δ is circular. See [`claim2_plan.md`](claim2_plan.md) §2B.

## Provenance / re-run

```bash
# on PeiLab2, /mnt3/wuyuhao/claim1_work
# human iN promoter (the headline), opening-only:
python3 eval_driver_claim2.py --contract /mnt3/wuyuhao/in_clean/get_in.tsv \
    --positives /mnt3/wuyuhao/neural_gt/neural.promoter.bed \
    --signed /mnt3/wuyuhao/in_clean/get_in.conf.tsv
# add --all-regions for the tougher baseline; swap neural.neighborhood.bed for neighborhood.
# mouse: --contract get_driver_scores.mm10.clean.tsv --signed get_clean_conf.tsv
#        --positives /mnt3/wuyuhao/mtf_loci/mtf.{promoter,neighborhood}.bed
```
