# HANDOFF → EpiAgent implementation (other Mac)

**Delete this file after reading** (`git rm HANDOFF_EPIAGENT.md && git commit`). It's a
one-shot note to fold EpiAgent into eCR_navigator the same way the other models went in.
Everything below is already on `main` as of commit 6882882.

## The pattern you're plugging into (don't reinvent it)

Each foundation model runs in its own GPU mirror and emits **one `.npz` embedding
artifact per cell state**; `navigate.py` diffs two states into the driver-score
contract. So EpiAgent = "produce `epiagent.MEF.<asm>.npz` + `epiagent.mES.<asm>.npz`,
then run navigate.py." No new pipeline.

- **Artifact contract:** `docs/embedding_artifact.md`. Arrays: `chrom` (`<U` str),
  `start`/`end` (int64), `embedding` (N×D float32), `meta` (0-d str JSON).
- **DO NOT hand-roll `np.savez`.** Call the shared writer
  `scripts/embedding_artifact.py::write_embedding_artifact(out, chrom, start, end, emb,
  model="epiagent", cell_state=..., assembly=..., source="epiagent_embed_regions.py")`.
  It coerces the load-critical dtypes — a bare `pandas.to_numpy()` gives object arrays
  that crash the navigator's `allow_pickle=False` load (this bit GET; the writer exists
  precisely to stop that recurring).
- **Readout:** per-region embedding shift MEF vs mES (skip any cell-level pooling; keep
  per-region/per-cCRE vectors), same as `chromfound_embed_regions.py` — copy that script
  as the closest template.
- **navigate.py flags are `--emb-a` / `--emb-b`** (not `--a/--b`).

## What's already staged on the servers (persistent `/yutiancheng`, all instances)

- `/yutiancheng/yuhao/eCR/artifacts/epiagent.{MEF,mES}.h5ad` — the EpiAgent cell×peak
  inputs you built earlier (~80 MB each). Confirm their assembly + var schema.
- `/yutiancheng/yuhao/eCR/epiagent_make_input.py` — your input builder.
- `/yutiancheng/yuhao/eCR/refs/` — `cCRE.bed`, `cCRE_document_frequency.npy` (EpiAgent
  refs), plus `mm10ToHg38.over.chain.gz` / `hg38ToMm10.over.chain.gz`.
- `EpiAgent` conda env on the model-zoo mirror (base `/yutiancheng/yuhao/miniconda3`).
- Reusable per-model scripts live in `/yutiancheng/yuhao/get_scripts/` (incl.
  `embedding_artifact.py`). liftOver binary that works despite glibc: the bioconda one
  at `/yutiancheng/yuhao/miniconda3/envs/scDIFF/bin/liftOver`.

## Species framing (decided, applies to EpiAgent too)

If EpiAgent is hg38/human-trained scATAC (very likely — check its cCRE universe), it is
**primary for human, reference for mouse** (bridge mm10→hg38 by liftOver; expect large
coverage loss on the mouse side). Native-mm10 primary tracks stay **GET + ChromBERT**.
Mirror the framing table in `docs/atacformer_pipeline.md` / `docs/chromfound_pipeline.md`.

## Deliverables to match the others

1. `scripts/epiagent_embed_regions.py` (+ any `epiagent_build_input.py`) — per-region
   embeddings → `.npz` via the shared writer.
2. `docs/epiagent_pipeline.md` — run steps + the mouse/human role note.
3. Stage results to `/yutiancheng/yuhao/eCR/artifacts/epiagent_driver_scores.<asm>.tsv`.

## Merge etiquette

`main` moved a lot since you branched (GET, ATACformer, ChromFound, the shared-writer
refactor). **Rebase your EpiAgent branch on latest `main` and don't force-push over the
model integrations.** The only shared-surface file you should touch is a *new*
`epiagent_*` script + a new doc — nothing else needs editing.

_Full running notes: memory `chromfound-scoping.md`, `ecr-extensions-roadmap.md`;
per-model docs under `docs/*_pipeline.md`._
