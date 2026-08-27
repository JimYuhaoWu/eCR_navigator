# Embedding artifact contract (mirror → navigator)

The **internal** interface between a server model mirror (ChromBERT, ChromFound,
GET, …) and eCR_navigator. Distinct from the *output* contract
(`region_weight_contract.md`, which eCR_predictor consumes). Each model runs in
its own conda env / GPU runtime on the HPCC; they cannot share one Python
process, so they hand off **files**, not objects.

## Why a file, not an import

Each foundation model has incompatible deps (torch versions, CUDA, custom ops)
and lives in a separate server mirror. eCR_navigator's own env stays light
(`numpy`/`pandas`/`sklearn` — no torch). A mirror runs the heavy embedding step,
writes an artifact; eCR_navigator reads it. Any model can be swapped in behind
this format.

## Format: `.npz` (numpy, dependency-light)

One artifact per **(model, cell_state)**. numpy-only so it loads in both the
mirror env and the light navigator env — no pyarrow/parquet dependency.

Filename convention: `{model}.{cell_state}.{assembly}.npz`
e.g. `chrombert.MEF.mm10.npz`, `get.MEF.mm10.npz`, `chromfound.MEF.hg38.npz`.

Arrays inside:

| Key | dtype | shape | Meaning |
|---|---|---|---|
| `chrom` | `<U…` (str) | `(N,)` | chromosome, matching assembly |
| `start` | `int64` | `(N,)` | 0-based BED start |
| `end` | `int64` | `(N,)` | exclusive BED end |
| `embedding` | `float32` | `(N, D)` | per-region model embedding |
| `signal` | `float32` | `(N,)` | **optional** per-region *scalar* accessibility for this state; present only for direction-capable models |
| `meta` | `<U…` (0-d) | `()` | JSON: `{model, cell_state, assembly, dim, source, has_signal, provenance}` |

### The optional `signal` array (direction channel)

Most models emit only `embedding`; the driver score is the L2 shift between two
states' embeddings — a **magnitude** (unsigned by construction; a norm cannot tell
"opens" from "closes"). A model that additionally has a **scalar per-state readout**
(one accessibility value per region) writes it as `signal`. navigate.py differences
the two states (`signal_b − signal_a`) into the signed `direction ∈ [-1, 1]` column of
the output contract — the open/close instruction. Rules:

- `signal[i]` corresponds to `embedding[i]` (same region, same order).
- Omit it (`signal=None`) and the direction column is simply absent — 4-column
  output, byte-compatible with existing consumers. `meta.has_signal` records which.
- A per-region `NaN` means **unmeasured** (the measured-external path could not cover
  that region in this state), which is kept distinct from a measured `0.0`. `navigate.py`
  differences the two states, and any region NaN in either state gets **no** `direction`
  (an empty field in the contract), never a fabricated open/close.
- **Two ways to fill it:**
  - *Model-native* — **EpiAgent**'s Signal-Reconstruction head
    (`sigmoid(signal_decoder(cell_embedding))`, a predicted-accessibility probability in
    `[0,1]`); AlphaGenome's DNase head. The model's own accessibility output.
  - *Measured (external)* — `scripts/attach_measured_signal.py` maps a real
    accessibility track (e.g. GET's aTPM) onto an artifact's regions by overlap. Used
    for models with **no** accessibility head (ATACformer, GET, ChromBERT, ChromFound);
    the model gives magnitude, the track gives direction.
- **Caveat (per model):** direction is principled from a model *designed* to predict
  accessibility (EpiAgent SR, AlphaGenome DNase) or from *measured* data (aTPM Δ — the
  trusted external case). A direction **synthesized from the embedding itself** (e.g.
  projecting the shift) for a non-directional model is a **navigator-side modification,
  not validated** — mark it as such wherever produced (see `region_weight_contract.md`).

**The dtypes are load-critical.** The navigator reads with `allow_pickle=False`
(fast, and it refuses to unpickle arbitrary objects), so `chrom` must be a unicode
## `meta.provenance` — which runtime produced this

Added 2026-08-27, when the models gained a second home (Hygon DCU alongside the
NVIDIA A800). The two backends give slightly but systematically different
numbers, so an artifact has to say how it was made.

```json
"provenance": {
  "host": "instance-...-wuyuhao-77258-0",
  "torch": "2.4.1",
  "torch_build": "2.4.1+das.opt1.dtk25041",
  "device": "cuda", "accelerator": "K100_AI",
  "arch": "gfx928:sramecc+:xnack-",
  "hip": "6.3.25211", "cuda": null,
  "dtype": "fp16", "attn": "flash"
}
```

Filled automatically by `scripts/ecr_runtime.py:runtime_provenance()`;
`dtype`/`attn` are supplied by the caller because only it knows how it ran the
forward pass. **Those two matter as much as the accelerator**: measured on
EpiAgent, like-for-like fp32 agreement between an A800 and a K100_AI is
~1.7e-06, whereas fp16 autocast or a different flash-attention kernel each shift
results by ~1e-01 at the tail.

`torch_build` is recorded separately because `torch.__version__` drops the local
build tag — the DTK build reports a bare `"2.4.1"`, indistinguishable from stock
torch.

Notes:
- **Additive and optional.** Existing readers that pull `model`/`cell_state`/…
  are unaffected; the arrays are unchanged.
- **Artifacts written before this date have no `provenance` key.** All of them
  were produced on the A800 mirrors. Treat a missing key as "A800, unknown
  dtype", not as "unknown machine".
- `write_embedding_artifact(..., provenance=...)` overrides the auto-capture,
  for the case where the artifact is *assembled* on a different machine than the
  one that computed the embedding — e.g. ChromBERT, where the mirror emits hdf5
  and `scripts/hdf5_to_artifact.py` converts it elsewhere. Auto-capture would
  otherwise record the converting host, which is wrong.

string array and `meta` a 0-d unicode-string array — an object array (e.g. a bare
`pandas.Series.to_numpy()`) fails to load. **Do not hand-roll the `np.savez`.** Every
model writes through the one shared helper `scripts/embedding_artifact.py`
(`write_embedding_artifact`), which coerces these dtypes so no model can drift; add a
new model by calling it, not by copying another model's save block.

## Alignment rule

Regions are aligned across cell states by the **`(chrom, start, end)` key**, not
by row order. A driver score is computed only for regions present in *both* cell
states' artifacts; the region set should be (at least) the peak union of the run's
assembly so every off-target region gets a weight.

## Assembly

Per-run, single-assembly — whatever the peaks/genome for *that* run use (e.g. **mm10**
for the MEF→mES mouse work, **hg38** for human, or a lifted assembly). Recorded in
`meta`; the loader enforces that the *two artifacts being diffed share an assembly*
(it does not pin a specific one). See CLAUDE.md "Multi-species requirement".
