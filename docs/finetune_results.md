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

## Next (per plan, not yet run)
- Rung 2: OSKMN(−Esrrb) ↔ JGES **binding**, same transition (train one panel, test the
  other).
- Stronger transition: human iN loci (needs the iN GET `.npz` + neural_gt gene-tagged loci).
