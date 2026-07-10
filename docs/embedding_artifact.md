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
| `meta` | `<U…` (0-d) | `()` | JSON: `{model, cell_state, assembly, dim, source, has_signal}` |

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
- **EpiAgent** fills it from its Signal-Reconstruction head
  (`sigmoid(signal_decoder(cell_embedding))`, a model-native predicted-accessibility
  probability in `[0,1]`). See `docs/epiagent_pipeline.md`.
- **Caveat (per model):** for a model *designed* to predict accessibility (EpiAgent's
  SR head, AlphaGenome's DNase head) the direction is principled. For a model with **no
  scalar accessibility readout**, any direction we synthesize (e.g. projecting the
  embedding shift, or borrowing a measured-peak Δ) is a **navigator-side modification,
  not validated** — mark it as such wherever it's produced (see
  `region_weight_contract.md`).

**The dtypes are load-critical.** The navigator reads with `allow_pickle=False`
(fast, and it refuses to unpickle arbitrary objects), so `chrom` must be a unicode
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
