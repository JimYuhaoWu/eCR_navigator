# Server model mirrors (HPCC) — access & layout

The genomic foundation models (ChromBERT, ChromFound, GET, EpiAgent, …) each run
in their **own conda env / GPU runtime mirror** on the server. They cannot share
one Python process, so they hand off **files** (see `embedding_artifact.md`), not
objects. This doc records how to reach them and what persists.

## Access — use the SSH key, not a password

A dedicated keypair was installed for automated access:

- Local private key: `~/.ssh/ecr_navigator` (ed25519).
- The **public** key lives in the mirror's `/root/.ssh/authorized_keys`.

Connect:

```bash
ssh -i ~/.ssh/ecr_navigator -p <PORT> root@<IP>
```

### Why key auth (important — restart behavior)

- **The login password rotates on every mirror restart.** Do NOT rely on it.
- **`/root` persists unchanged across restarts** — so `/root/.ssh/authorized_keys`
  (and everything under `/root`: the models, checkpoints, envs) survives. Key auth
  therefore keeps working through password rotations with no re-setup.
- The **only** thing a restart may change that a new session needs from the user
  is the **IP and/or port**. Ask for the current `IP:PORT`; the key handles the
  rest. If key auth ever fails after a restart, re-append the public key once
  (using that session's password) and it persists again.

Last known endpoint: `root@172.16.78.10 -p 35963` (verify — may have rotated).

## ChromBERT mirror (verified 2026-07-03)

- Env: `/opt/conda` (base). torch 2.3, CUDA on an **A100-80GB**. ChromBERT v1.1.1
  at `/root/ChromBERT`; installed package under `/opt/conda/.../chrombert`.
- **Import quirk:** set `PROTOCOL_BUFFERS_PYTHON_IMPLEMENTATION=python` before
  `import chrombert` (stale protobuf/onnx conflict), else the import fails.
  `scripts/mirror_env.sh` sets this for you.
- **Both species are NATIVE — no liftOver for mouse.** ChromBERT ships hg38
  (6k regulators) and **mm10 (5k regulators)**, both at 1 kb. `-g hg38` / `-g mm10`
  is a first-class CLI flag. The mm10 checkpoint + grid + 8.1 GB cistrome hdf5 were
  fetched via `chrombert_prepare_env --genome mm10` and confirmed working. liftOver
  is therefore NOT needed for ChromBERT; keep it only as a fallback for other
  models that genuinely lack a species.
- **Runtime deps beyond the conda env:**
  - `bedtools` — needed by `chrombert_make_dataset`; NOT in the conda env. Installed
    as a **static binary at `/root/bin/bedtools`** (persists under /root; put
    `/root/bin` on PATH). Do NOT `conda install bedtools` — it drags in a python
    rebuild that can perturb the ChromBERT env.
  - HF downloads: use `HF_ENDPOINT=https://hf-mirror.com` (huggingface.co is
    unreliable from this network; direct hits get connection-reset).

### Region-embedding workflow (what `scripts/run_chrombert_region_emb.sh` runs)

For one cell state, on the mirror (GPU):

```
chrombert_make_dataset  <peaks.bed> -g <genome> -o dataset.tsv
    # overlaps peaks onto the 1 kb grid. DO NOT pass --no-filter — that emits the
    # ENTIRE genome grid (~2.1M rows hg38 / ~530k+ mm10, ~11 h). Default filters
    # to only peak-overlapping bins. tsv cols: chrom,start,end,build_region_index,label
chrombert_get_region_emb  dataset.tsv -g <genome> -o emb.hdf5
    # writes region (N,4 int64 = [chrom_id,start,end,build_region_index]) + emb (N,768 f16)
```

`scripts/hdf5_to_artifact.py` then joins `build_region_index` back to the tsv for
the chrom string and writes the `.npz` embedding artifact. Verified end-to-end on
hg38 (5-region smoke test) 2026-07-03.

- Underlying API if you script it directly: `chrombert.ChromBERTEmbedding` →
  `get_region_embedding()` (768-dim mean-pooled TRN embedding per 1 kb locus).
  Cell-state conditioning is via **cell prompts** or **fine-tuning**, not the bare
  pretrained model — relevant when we pick the driver-score readout.
- Directly relevant example: `examples/tutorials/
  tutorial_transdifferentiation_chromatin_accessibility.ipynb` — fibroblast→
  myoblast, structurally the same as MEF→mESC. Recipe: merge both states' peaks →
  align to 1 kb bins → log2 ATAC fold-change label → fine-tune → key-regulator
  attribution + up/unchanged region classification.

### Persistence caveat

The mirror does **not** save changes unless you snapshot it. The static
`/root/bin/bedtools` and the downloaded `~/.cache/chrombert` data (incl. the 8.1 GB
mm10 hdf5) live under paths that persist *only if snapshotted*. If a restart loses
them, `scripts/setup_mirror.sh <genome>` re-provisions idempotently.
