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

## Top-nomination sweep (the metric that matches the navigator's job)

The navigator nominates a **limited top-k** of regions for `eCR_predictor`, so the
decision-relevant question is not average AUROC but: *how enriched for master-TF loci are
the very top regions under `driver_score` vs under the directional signed-Δ score?* Fold
enrichment over base rate, computed on the **full** contract ranking (PeiLab2
`claim1_work/topk_sweep.py`; raw ranking, not |Δ|-matched — but the driver-vs-signed
head-to-head is fair because signed-Δ *is* the magnitude scorer):

| System | Positive | | top-0.5% | top-1% | top-2% | top-5% | top-10% |
|---|---|---|---|---|---|---|---|
| **iN** | promoter | driver | 4.7× | **9.4×** | 4.7× | 2.1× | 1.7× |
| | | signed-Δ | 0.0× | 2.4× | 1.2× | 1.2× | 0.7× |
| **iN** | neighborhood | driver | 1.7× | **9.8×** | 4.9× | 2.1× | 1.4× |
| | | signed-Δ | 0.4× | 1.1× | 0.9× | 0.8× | 0.8× |
| mouse | promoter | driver | 5.8× | 2.9× | 4.3× | 2.9× | 2.5× |
| | | signed-Δ | 5.8× | **7.2×** | **9.4×** | **5.2×** | 3.3× |
| mouse | neighborhood | driver | 2.9× | 1.4× | 1.9× | 1.7× | 1.5× |
| | | signed-Δ | 5.2× | **4.9×** | **4.3×** | **2.7×** | 2.5× |

- **Human iN — `driver_score` owns the top.** top-1% is **~9–10×** enriched for master-TF
  loci; the directional score's top regions sit **at or below base rate** (0–2.4×). At the
  cut that matters, the model beats direction by 4–9×. Enrichment is **concentrated in the
  top ~1%** and decays fast (9.4× → 2.1× by top-5% → 1.7× by top-10%) — evidence that a
  *small* nomination set is where the confidence lives.
- **Mouse — reversed.** signed-Δ's top nominations are the more enriched (5–9× at the top)
  than driver's (1.4–2.9×). Mouse master-TF promoters simply *are* the biggest-opening
  regions, so ranking by accessibility change finds them better than the model — the
  concrete mechanism behind mouse being "largely directional."

### Per-model, per-species top-k confidence (all models, driver_score only)

Extends the sweep beyond GET to every model tested, per species — the confidence that a
top-k nomination is a real master-TF locus, at each cut (fold enrichment over that model's
own base rate; `claim1_work/topk_allmodels.py`). `npos` = master-TF positives in that
model's universe; small `npos` ⇒ noisy fold.

**MOUSE mm10 — master-TF promoter**
| model | N | npos | top-0.5% | top-1% | top-2% | top-5% | top-10% |
|---|---|---|---|---|---|---|---|
| **GET** | 86,956 | 69 | 5.8× | **2.9×** | 4.3× | 2.9× | 2.5× |
| ChromFound | 40,365 | 49 | 0× | 0× | 4.1× | 2.5× | 2.0× |
| ChromBERT | 29,343 | 26 | 0× | 0× | 0× | 0× | 0× |

**MOUSE mm10 — master-TF neighborhood**
| model | N | npos | top-0.5% | top-1% | top-2% | top-5% | top-10% |
|---|---|---|---|---|---|---|---|
| **GET** | 86,956 | 346 | 2.9× | 1.4× | 1.9× | 1.7× | 1.5× |
| ChromFound | 40,365 | 154 | 1.3× | 1.3× | 3.3× | 1.8× | 1.2× |
| ChromBERT | 29,343 | 80 | 0× | 0× | 0× | 0× | 0× |

**HUMAN iN hg38 — master-TF promoter**
| model | N | npos | top-0.5% | top-1% | top-2% | top-5% | top-10% |
|---|---|---|---|---|---|---|---|
| **GET** | 329,983 | 85 | 4.7× | **9.4×** | 4.7× | 2.1× | 1.7× |
| ChromBERT | 144,659 | 39 | 0× | 0× | 3.9× | 1.5× | 0.8× |
| ChromFound | 348,786 | 88 | 0× | 0× | 0.6× | 0.5× | 0.6× |
| ATACformer | 112,920 | 22 | 0× | 0× | 0× | 0× | 0× |
| EpiAgent | 3,346 | **0** | — | — | — | — | — |

**HUMAN iN hg38 — master-TF neighborhood**
| model | N | npos | top-0.5% | top-1% | top-2% | top-5% | top-10% |
|---|---|---|---|---|---|---|---|
| **GET** | 329,983 | 461 | 1.7× | **9.8×** | 4.9× | 2.1× | 1.4× |
| ChromBERT | 144,659 | 204 | 0× | 0× | 2.5× | 1.0× | 0.7× |
| ChromFound | 348,786 | 494 | 0× | 0× | 0.5× | 0.5× | 0.7× |
| ATACformer | 112,920 | 143 | 0× | 0× | 0.4× | 0.3× | 0.4× |
| EpiAgent | 3,346 | 5 | 39×¹ | 20×¹ | 10×¹ | 8×¹ | 4×¹ |

¹ EpiAgent's large folds are **1–2 hits out of 5 positives** — noise, not signal; its
8,190-cCRE cap leaves too few master-TF loci in-universe to evaluate (0 promoters overlap).

**Read across models:** **GET is the only model with usable top-k confidence**, and it is
sharpest on human iN (top-1% ~9–10×, front-loaded). ChromFound has a *weak mid-tail* mouse
signal (top-2–5%, ~2–4×) but **nothing at the very top** and is at/below base rate on human.
ChromBERT and ATACformer are **null at every cut** (their occasional top-2% blips are ≤3
hits). EpiAgent is **too sparse to nominate** (0–5 in-universe positives). So for target
nomination, trust **GET's top ~1% on a strong clean transition**; no other model earns a
top-k nomination on this evidence.

## Nomination policy — which score to trust for a new transition (`scripts/preflight.py`)

The scores above are transition-dependent, so the platform needs a *preflight* that decides,
per transition, whether to nominate from GET `driver_score` or from signed-Δ — runnable
**before** nominating, from data available at inference (endpoint accessibility + known
target-cell biology; no drivers of *this* transition, no wet-lab ground truth).

**First, a negative result that shapes the design.** No endpoint-only indicator we tried
separates "GET wins" from "signed-Δ wins":
- **Cleanliness (PC1 / endpoint separation) does NOT predict it.** Recorded PC1: iN ≈0.89,
  mouse MEF→mES ≈0.96–0.98, dropped iCM ≈0.68. Mouse is the *cleanest* yet GET *lost* there.
- **Rank-divergence from magnitude does NOT predict it.** GET is ~equally weakly correlated
  with |signed-Δ| in both systems (Spearman ≈0.21 mouse / ≈0.16 human; top-1% overlap ≈0.02
  both). GET departs from "what opens" in *both* — divergence doesn't say *toward* drivers.

So we do **not predict** "clean enough"; we **measure** it per transition. Two gates:

**Gate 1 — admissibility (endpoint-only; a reliable REJECT, not a reliable ADMIT).** From a
per-replicate accessibility matrix: require ≥2 replicates/state, within-state minus
across-state correlation ≥ 0.10, and PC1 ≥ 0.80. Screens the "nothing works" failure mode
(weak/partial transition, e.g. the dropped iCM). **Clearing it does not mean GET will win** —
mouse clears it (PC1 0.956, coherence +0.845) and GET still loses. Gate 2 is the real decision.

**Gate 2 — score selection (the decision).** eCR design always targets a *known* cell type,
so its canonical master-TF loci are known biology, independent of this transition. Drop them
into the Claim-2A harness as positives and read the **stable** statistics (paired ΔAUROC CI +
incremental LR — not a few-hit top-1% fold):
- driver **beats** signed (ΔAUROC CI excludes 0 in driver's favour **or** incremental-LR
  p<0.05 with positive driver coef) → **PRIMARY = GET `driver_score` top ~1%**.
- otherwise → **PRIMARY = signed-Δ top-k** (GET supplementary).
- In both, the measured signed-Δ rides along per region for open/close **direction**.

**Validated on both transitions (`preflight.py`, run 2026-07-15):**

| Transition | Gate 1 | Gate 2 (85/69 anchors) | PRIMARY |
|---|---|---|---|
| human iN | (per-rep matrix not on disk; QC PC1≈0.89) | driver 0.664 vs signed 0.494, ΔAUROC CI[+0.081,+0.264], LR p=0.001 | **GET top ~1%** |
| mouse MEF→mES | ADMIT (PC1 0.956, coherence +0.845) | driver 0.579 vs signed 0.507, ΔAUROC CI[−0.045,+0.188], LR p=0.13 | **signed-Δ top-k** |

Mouse is the instructive case: **admissible yet driver-not-primary** — exactly why Gate 1
alone is insufficient. Thresholds (PC1≥0.80, coherence≥0.10, "CI>0 or LR p<0.05") are
first-pass, calibrated on **n=2** transitions; tighten as more accrue. Needs ≥~20 in-universe
anchors for a reliable Gate-2 verdict (widen promoter→+neighborhood if sparse).

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
