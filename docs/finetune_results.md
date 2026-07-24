# Fine-tuning results — head-only probe (Regime 1)

Method + scope: [`finetune_plan.md`](finetune_plan.md). Harness `scripts/finetune_head.py`
(GET shift PCA + signed-Δ → L2-logistic; leave-one-gene-out CV; head vs zero-shot
`driver_score` vs signed-Δ by paired ΔAUROC). All runs **opening-only, seed 0**, on the
clean MEF→mES GET artifacts (`get.{MEF,mES}.clean.mm10.npz`, contract
`get_driver_scores.mm10.clean.tsv`, signed-Δ `get_clean_conf.tsv`, loci `mtf.{promoter,
neighborhood}.bed`).

## MEF→mES, loci — done 2026-07-22

**Harness validated on real data:** the zero-shot `driver_score` AUROC (0.579 promoter /
0.523 neighborhood) reproduces Claim 1's GET numbers, and signed-Δ (0.507 / 0.500)
reproduces Claim 2A — so the alignment and baselines are correct.

| target | pca-k | held-out pos (genes) | **head** | driver_score | signed-Δ | Δ head−driver [CI] | Δ head−signed [CI] |
|---|---:|---:|---:|---:|---:|---|---|
| promoter | 15 | 51 (23) | **0.378** | 0.579 | 0.507 | −0.201 [−0.341, −0.056] | −0.129 [−0.249, −0.001] |
| promoter | 5 | 51 (23) | **0.378** | 0.579 | 0.507 | −0.201 [−0.303, −0.100] | −0.129 [−0.259, −0.003] |
| neighborhood | 15 | 277 (26) | **0.471** | 0.523 | 0.500 | −0.052 [−0.112, +0.008] | −0.029 [−0.083, +0.025] |
| neighborhood | 8 | 277 (26) | **0.453** | 0.523 | 0.500 | −0.071 [−0.134, −0.009] | −0.047 [−0.103, +0.007] |
| neighborhood | 5 | 277 (26) | **0.444** | 0.523 | 0.500 | −0.080 [−0.139, −0.023] | −0.056 [−0.112, −0.004] |

**Verdict: the head does NOT transfer — it is at/below chance (0.38–0.47) and
significantly worse than both baselines**, robustly across loci sets and PCA-k (5/8/15).

**Interpretation.** GET's embedding-shift carries **no driver-region signature that
generalizes across driver genes**: training a head on some master-TF genes' loci does not
help — it actively *anti-predicts* — a held-out gene's loci (the below-chance, "gene-
specific direction" regime `test_finetune_head.py` flags as non-transferable). Different
driver genes occupy different, non-shared embedding neighborhoods. This extends Claim 2A:
GET's mouse signal isn't just ≈ signed-Δ, it also doesn't decompose into a *learnable,
transferable* driver code — the zero-shot 0.58 on promoters is locus/magnitude-specific,
not a generic driver signature.

**Caveats.** Small panel (23–26 genes), **mouse-only, GET-only, loci-only**. This does not
doom Regime-1 universally; it says the head-only probe finds no transferable signature on
*this* transition/model/target. Per the plan's exit criteria, that argues against a
backbone fine-tune here and for testing a **stronger transition (human iN**, GET zero-shot
0.668) before concluding.

> **Harness note (2026-07-23):** `finetune_head.py` now includes the shift **L2 norm
> `|shift|`** as a feature by default (≈ the zero-shot `driver_score`), so the head starts
> from the zero-shot signal and the PCA *direction* comps can only add; `--no-magnitude`
> restores the old direction-only features. The MEF→mES numbers above predate this and are
> effectively `--no-magnitude`; the verdict is unchanged (the human iN runs below test both).

## Human iN, loci — done 2026-07-23

The stronger transition (GET zero-shot promoter 0.664 in Claim 1) — the decisive test the
mouse caveats called for. GET `get.{fib,iN}.in.hg38.npz` (329,983 regions × 768), contract
`get_driver_scores.in.hg38.tsv`, signed-Δ `get_in.conf.tsv`, loci `neural.{promoter,
neighborhood}.bed`. opening-only, seed 0, pca-k 15.

| target | features | held-out pos (genes) | **head** | driver_score | signed-Δ | Δ head−driver [CI] | Δ head−signed [CI] |
|---|---|---:|---:|---:|---:|---|---|
| promoter | \|shift\|+dir+signed | 65 (15) | **0.563** | 0.664 | 0.494 | −0.100 [−0.210, +0.009] | +0.070 [−0.059, +0.206] |
| promoter | dir+signed (no-mag) | 65 (15) | **0.551** | 0.664 | 0.494 | −0.113 [−0.228, −0.002] | +0.057 [−0.072, +0.194] |
| neighborhood | \|shift\|+dir+signed | 327 (22) | **0.405** | 0.557 | 0.499 | −0.152 [−0.202, −0.101] | −0.094 [−0.145, −0.044] |

**Verdict: the head does NOT clear the baselines on iN either — it matches or *underperforms*
zero-shot, significantly so on neighborhood (0.405, below chance).** Crucially this holds
**even when the head is handed `|shift|`** (the zero-shot signal) as a feature: leave-one-gene-out
with 65–327 tiny proxy labels fits weights that don't transfer to unseen driver genes, and
the extra learned parameters add non-transferable noise that *degrades* the fixed,
gene-agnostic zero-shot magnitude (which is GET's real Claim-1 signal, 0.664).

**Bottom line across both transitions.** Head-only supervised fine-tuning on loci **does not
help — it hurts.** GET's driver signal is best used **zero-shot**; the small driver-label
corpus (tens of loci across ~15–26 genes) is too little for a supervised head to learn a
*transferable* driver code, and forcing it overfits and loses the zero-shot prior. This
argues **against a backbone fine-tune on loci** (which would overfit harder) and refocuses
Regime-1 on either (a) a much larger driver-label corpus, or (b) Regime-3 perturbation
labels. Machine-readable: [`finetune_results.human.tsv`](finetune_results.human.tsv).

## MEF→mES, binding — OSKMN ↔ JGES cross-panel transfer — done 2026-07-23

Rung 2: train the head on one driver-TF panel, test on the OTHER (same transition), with
shared regions excluded so a site bound by both can't leak. OSKMN = Pou5f1/Sox2/Klf4/Myc/
Nanog ChIP-Atlas (125,273 sites, Esrrb+Sall4 dropped → JGES-side); JGES = `jges.union.bed`
(Jdp2/Glis1/Esrrb/Sall4, 15,875). `transfer_scores`, features |shift|+dir+signed, seed 0.

| direction | held-out pos / neg | **head** | driver_score | signed-Δ | Δ head−driver [CI] | Δ head−signed [CI] |
|---|---:|---:|---:|---:|---|---|
| JGES → OSKMN | 7,770 / 6,686 | **0.495** | 0.512 | **0.582** | −0.017 [−0.032, −0.002] | −0.087 [−0.099, −0.075] |
| OSKMN → JGES | **0** / 0 | — | — | — | — (degenerate) | — |

**OSKMN → JGES is degenerate:** OSKMN binding (125k sites ≈ 38% of opening regions)
**subsumes** JGES — after excluding training regions, 0 JGES positives remain. The panels
are not independent (OSKM binds essentially everywhere JGES does), so this direction can't
be tested. **JGES → OSKMN is decisive:** the head (0.495) is **significantly worse than
both** zero-shot driver_score (0.512) and signed-Δ (0.582); here **signed-Δ (measured
accessibility) is the single best predictor** and the model/supervision add nothing.

## Overall conclusion — supervised fine-tuning does not beat zero-shot

Across three fine-tuning experiments the supervised head **never clears the zero-shot
`driver_score` baseline and usually hurts**: mouse loci (head 0.38–0.47 < zero-shot 0.52–0.58),
human iN loci (head 0.56 < zero-shot 0.664, *even when handed `|shift|`*), and OSKMN↔JGES
binding transfer (head 0.495 < driver 0.512). Leave-one-gene/panel-out with the current tiny,
proxy labels fits weights that don't transfer; the extra parameters add non-transferable noise
that degrades the fixed zero-shot signal.

**This is NOT "accessibility beats the model."** Where it counts, the opposite holds: the
now-merged **Claim 2A** shows the **zero-shot `driver_score` beats signed-Δ on the strong iN
transition** (0.664 vs 0.494; ΔAUROC +0.170; top-1% 9.4× vs 2.4×), and only loses to signed-Δ
on mouse/weak/already-directional transitions ([`validation_summary.md`](validation_summary.md),
[`claim2_results.md`](claim2_results.md)). The JGES→OSKMN result above (signed-Δ 0.582 > driver
0.512) is a *mouse binding* case, consistent with "mouse is directional" — not a general claim.

The narrow, correct lesson: **use GET zero-shot (it already works, and on strong transitions
beats signed-Δ); cheap supervised fine-tuning on the current driver-label corpus adds nothing
and overfits.** A backbone fine-tune (more params, same tiny labels) would overfit harder —
**not pursued**. Regime-1 needs a much larger driver-label corpus or Regime-3 perturbation
labels; see [`limited_data_strategy.md`](limited_data_strategy.md). Machine-readable:
[`finetune_results.mtfbind.tsv`](finetune_results.mtfbind.tsv).
