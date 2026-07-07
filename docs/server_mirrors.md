# Server model mirrors (HPCC) — access & layout

> **Onboarding a new mirror?** Follow the reusable, model-agnostic playbook in
> [`mirror_onboarding.md`](mirror_onboarding.md). This file holds the persistence
> details and the per-model specifics (ChromBERT below; add a section per model).

The genomic foundation models (ChromBERT, ChromFound, GET, EpiAgent, …) each run
in their **own conda env / GPU runtime mirror** on the server. They cannot share
one Python process, so they hand off **files** (see `embedding_artifact.md`), not
objects. This doc records how to reach them and what persists.

## Mirrors index

| Mirror | Endpoint | conda base | Models | Notes |
|---|---|---|---|---|
| ChromBERT | `172.16.78.10:35963` | `/opt/conda` | ChromBERT (hg38+mm10) | A100-80GB; verified working |
| Models zoo | `172.16.78.10:35364` | `/yutiancheng/yuhao/miniconda3` | GET, EpiAgent, atacformer, alphagenome, scDIFF, SPG, BindCraft | key onboarded; `source /yutiancheng/yuhao/yuhao.sh` |

Same IP, different ports and host keys per instance → `mirror_env.sh` disables
host-key pinning. Port may change on restart (`MIRROR_PORT=<port>`).

## Access — use the SSH key, not a password

A dedicated keypair was installed for automated access:

- Local private key: `~/.ssh/ecr_navigator` (ed25519).
- The **public** key lives in the mirror's `/root/.ssh/authorized_keys`.

Connect:

```bash
ssh -i ~/.ssh/ecr_navigator -p <PORT> root@<IP>
```

### Restart behavior (important — the platform wipes the SSH key)

This is a university-built managed GPU platform. Empirically:

- **The data disk persists** across restarts: `~/.cache/chrombert` (models),
  `/root/bin`, `/root/.bashrc`, etc. all survive — even the 8 GB mm10 hdf5.
- **BUT `/root/.ssh/authorized_keys` is regenerated on every boot**, discarding
  our key — snapshot or not. (Confirmed: saved a snapshot with the key present,
  restarted, key was gone.) So key auth breaks after every restart unless
  re-injected.
- **Self-heal via `.bashrc`:** since `.bashrc` persists but `authorized_keys`
  doesn't, `setup_mirror.sh` adds a silent, idempotent block to `/root/.bashrc`
  that re-appends our key whenever a shell starts. `.bashrc` can't fix the *first*
  SSH login (chicken-and-egg), but the platform's **web terminal** logs in without
  SSH and sources `.bashrc` — so **after each restart, open the web terminal once**
  and the key is restored; then automated SSH works for the rest of the session.
  (If the platform runs a login shell at boot, this becomes fully automatic.)
- Endpoint: **IP fixed** at `172.16.78.10`; **port may change** — override with
  `MIRROR_PORT=<port>`. Last known `-p 35963`.
- A `.bashrc` backup is kept at `~/.bashrc.ecr_bak` in case the block ever needs
  reverting.

Alternative durable fix: if the platform dashboard ever exposes an SSH-public-key
or startup-script field, register the key there instead — that removes the
open-the-terminal-once step.

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

### Persistence & snapshot policy

The mirror does **not** save changes unless you snapshot it, and snapshotting is
slow — so snapshot **rarely**. The trick is separating the *persistent environment
layer* (snapshot once) from *per-run transients* (pushed from this repo at runtime,
never snapshotted).

**Bake into ONE snapshot (changes ~never):**
- `~/.cache/chrombert` model/data (mm10 + hg38) — the expensive ~8 GB download
- `/root/bin/bedtools` (static binary)
- ChromBERT install (`/opt/conda`, `/root/ChromBERT`)
- `/root/.ssh/authorized_keys` — **critical**: without it, key auth is lost on every
  restart and must be re-added by hand (see Access above)
- `/root/.bashrc` provisioning block (PATH + protobuf/HF env) — added idempotently
  by `scripts/setup_mirror.sh`

**Never needs snapshotting (transient):** the `scripts/` (live in this repo), and
all per-run files (peak BEDs, datasets, hdf5, artifacts) — created in `/tmp`,
fetched back to the local repo.

**Re-snapshot only when the persistent layer grows:** adding a new species/model
checkpoint, standing up a second model mirror (ChromFound/GET with its own env), or
pinning a fine-tuned checkpoint to reuse. Routine embedding runs never require it.
If you produce a fine-tuned checkpoint (~218 MB), prefer scp'ing it back to
versioned storage over re-snapshotting.

If a restart ever loses the environment layer, `scripts/setup_mirror.sh <genome>`
re-provisions it idempotently (bashrc + bedtools + data download).

## GET mirror (models zoo, port 35364) — scoped 2026-07-04

- Env `get` (conda at `/yutiancheng/yuhao/miniconda3`): py3.12, torch 2.6,
  **`get_model` v0.1.0** (GET-Foundation), scanpy/anndata/zarr. Repo at
  `/yutiancheng/yuhao/get_model` (tutorials for prepare/finetune/infer).
- **Multi-species: mm10 IS supported** — `zarr_dataset.py` maps `mm10 -> M36`
  (GENCODE); dataloader takes an `assembly` param (`hg38`/`mm10`). GET operates on
  per-region **motif-score matrices (RegionMotif)** + ATAC, not raw coordinates, so
  cross-species transfer is more natural than for sequence models.
- Framework: **hydra config + Lightning**. Inference via `RegionMotifDataset` /
  `InferenceRegionMotifDataset` + a config yaml; not a one-line CLI like ChromBERT.
- Input is a **zarr** (regions x motif features + ATAC). Demo `pbmc10k_multiome.zarr`
  (254 MB) is already on the instance; tutorial data also on Zenodo
  (astrocyte / pbmc `.zarr.tar`).
- **BLOCKER — pretrained weights:** checkpoints live ONLY in a **Requester-Pays S3
  bucket** `s3://2023-get-xf2217/get_demo/` (prefixes
  `pretrain_human_bingren_shendure_apr2023/fetal_adult/`,
  `Interpretation_all_hg38_allembed_v4_natac/`). No public HTTPS/Zenodo/HF mirror
  exists (checked). Needs AWS credentials (`aws s3 cp --request-payer requester …`)
  or an author-provided copy (Xi Fu, xf2217@cumc.columbia.edu / GET-Foundation
  GitHub). `get_model/utils.py` also accepts an `https://` checkpoint URL if a
  mirror is found.
- **Status:** integration scoped; blocked on obtaining pretrained weights. Once a
  checkpoint lands: inspect per-region outputs (embedding vs TF-importance) → map to
  the `.npz` embedding-artifact contract → same `navigate.py` path as ChromBERT.

### GET checkpoint provisioned (2026-07-07)

- **Pretrained foundation weights on the instance:**
  `/yutiancheng/yuhao/models/get/pretrain_fetal_adult/checkpoint-799.pth`
  (978 MB, `/yutiancheng` = persistent shared storage → no snapshot needed).
- **Sourcing:** checkpoints are in the Requester-Pays S3 bucket
  `s3://2023-get-xf2217/get_demo/checkpoints/regulatory_inference_checkpoint_fetal_adult/pretrain_fetal_adult/checkpoint-799.pth`.
  S3 egress to the HPCC is throttled (~74 KiB/s) so we downloaded to the local PC
  (`aws s3 cp --request-payer requester`, ~178 KiB/s) then sftp-resumed
  (`put -a`) up to the instance. Verified: exact size 1025979643 bytes; loads via
  torch (model/optimizer/epoch/args; 85.5M params).
- **Architecture (from ckpt args):** `model=pretrain_geneformer_base`,
  `input_dim=283` (282 motif clusters + 1 ATAC channel), `num_region_per_sample=200`,
  `use_natac=True`. Layers for interpretation: `region_embed`, `encoder.blocks.0..11`,
  `encoder.norm`.
- **Per-region readouts available (for driver scoring):** (1) per-region embedding
  (region_embed / encoder block) → shift MEF vs mES [parallels ChromBERT, same .npz
  artifact path]; (2) predicted regulatory-activity shift; (3) in-silico region
  perturbation importance (zero a region's input, measure output change — most causal).
- **Next:** run inference on the present `pbmc10k_multiome.zarr` demo to see real
  outputs, then pick the readout. AWS CLI is set up locally (`python -m awscli`,
  creds in `~/.aws`) for any further GET data (motif clustering / interpretation).

### GET RegionMotif prep (mm10 + hg38) — 2026-07-07

GET's input is `region_motif = [282 motif-archetype scores | 1 ATAC (aTPM)] = 283`.
The 282-motif matrix must be built for ANY new peak set — **human needs this exact
step too**; the demo zarr just ships it precomputed. Only the assembly differs.

- **Motif source:** Vierstra archetype motifs, pre-scanned genome-wide, tabix-indexed:
  `https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/{assembly}.archetype_motifs.v1.0.bed.gz`
  (hg38 12 GB, **mm10 9.5 GB** — both exist). We query REMOTELY (`tabix -R`), fetching
  only blocks over the peaks → no full download.
- **`scripts/get_regionmotif_matrix.py`** (runs in the `get` env; needs tabix+bedtools):
  peaks + `--assembly {hg38,mm10}` + canonical 282 `motif_names.txt` → per-peak
  282-motif `.npz`. Validated on tiny mm10 peaks (2×282, sensible top motifs).
  Canonical motif order = the pbmc zarr's `motif_names` (alphabetical, AHR…ZSCAN4).
- **Remaining input gap — aTPM (the 283rd channel):** normalized ATAC per peak per
  cell state (MEF, mES). `add_atpm` reads a BED with `Name` + `aTPM` columns. We have
  MACS peaks (score = -log10 q, NOT aTPM) → need ATAC counts/CPM per peak from the
  fragment/signal data to compute true aTPM. This is the cell-state channel that makes
  MEF vs mES embeddings differ, so it must be real accessibility, not the MACS score.
- **Then:** region_motif = [motif | aTPM] → window 200 regions → GET encoder → per-region
  768-d embedding (proven by the probe) → shift MEF vs mES → `.npz` artifact → navigate.py.
