# Fine-tuning plan (Regime 1) — head-only probe, loci, leave-one-gene-out

Follows the zero-shot validation (Claim 1: only GET is informative on clean transitions;
[Claim 2A](claim2_results.md): GET's `driver_score` does not beat signed-Δaccessibility).
Fine-tuning is **Regime 1** from [CLAUDE.md](../CLAUDE.md): train on validated driver
labels from *completed* transitions, deploy **endpoint-only** on new pairs. Inference
stays endpoint-only; only *training-time supervision* is added.

**The bar (higher than zero-shot):** a fine-tuned model must beat **both** zero-shot
`driver_score` **and** the signed-Δ baseline — else supervision bought nothing over
measured accessibility.

## Scope (decided with the user)

- **Head-only probe first** — a low-capacity supervised head on FROZEN embeddings (reuses
  the `.npz` artifacts, can't overfit 768-d on ~50 positives). Backbone fine-tune only if
  the probe shows headroom.
- **Loci target first** (master-TF *genes'* regulatory regions — GET's only Claim-1-positive
  signal); binding (OSKMN↔JGES) is a separate later target.
- **One transition at a time** — MEF→mES first. **No cross-cell-type transfer yet.**
- Generalization test within a single transition + panel = **leave-one-driver-gene-out**:
  a locus of gene *g* is scored only by a head trained without any of *g*'s loci.

## Method (`scripts/finetune_head.py`)

```
features  = [ PCA(emb_end - emb_start, k) | signed-Δ ]     # PCA label-blind, fit once
head      = L2-logistic
positives = master-TF loci (gene-tagged BED col 4)
negatives = |signed-Δ|-matched background (Claim-1 sampler)
CV        = leave-one-GENE-out; pool out-of-fold scores -> one held-out ranking
verdict   = paired ΔAUROC of head vs zero-shot driver_score AND vs signed-Δ
```
Opening-only by default (matches Claim 1/2): the shift **vector** is still full-dim, so
the head can exploit embedding *direction* even though accessibility direction is fixed.
Pure numpy, self-contained `.npz` loading (runs on a bare mirror); shares primitives with
`eval_driver_claim1.py` / `eval_driver_claim2.py`. Tests: `tests/test_finetune_head.py`
(a **shared** driver direction transfers; **gene-specific** directions do NOT — the CV
must expose that, which it does).

## Transfer ladder (only the first rung is in scope now)

1. within-transition, leave-one-gene-out — **MEF→mES loci (DONE — see finetune_results.md)**
2. same transition, unseen TF panel — OSKMN↔JGES **binding** (next)
3. cross-cell-type (MEF→mES ↔ fib→iN) — deferred (user: not yet)

## Honest expectation / exit criteria

If the head can't clear the leave-one-gene-out gate on the strongest zero-shot signal
(GET loci), Regime-1 supervision does not generalize *here* and a backbone fine-tune is
not yet warranted — pivot to a stronger transition (human iN) or to Regime 3 (perturbation
labels). Results: [`finetune_results.md`](finetune_results.md).
