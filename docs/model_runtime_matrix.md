# Model runtime & GPU matrix (5 integrated models)

Per-model **inference** runtime facts for the eCR_navigator embedding step: peak GPU
memory, the GPU/driver it was measured on, and each conda env's PyTorch/CUDA/Python.
Use it to pick a GPU when standing up a mirror. Measured **2026-07-13**.

## The table

| Model | Peak GPU mem (inference) | Measured on (input) | Recommended GPU | PyTorch | torch CUDA build | Python |
|---|---|---|---|---|---|---|
| **EpiAgent** | **7,706 MiB (~7.5 GB)** | 8,190-token fwd, flash-attn on | **≥16 GB, Ampere+ (sm_80)** for flash-attn (V100 falls back, slower) | 2.0.1+cu117 | 11.7 | 3.11.13 |
| **ChromBERT** | 2,388 MiB (~2.3 GB) | 7,295 regions, hg38-6k | ≥8 GB, any arch | 2.3.0a0+…nv24.04 (NGC) | 12.4 | 3.10.12 |
| **ChromFound** | 1,842 MiB (~1.8 GB) | 167,488 OCRs × 2 states | ≥6 GB, any arch | 2.2.2+cu121 | 12.1 | 3.10.19 |
| **ATACformer** | 1,216 MiB (~1.2 GB) | 66,573 regions | ≥6 GB, any arch | 2.6.0+cu124 | 12.4 | 3.12.11 |
| **GET** | 896 MiB (~0.9 GB) | 167,488 regions | ≥6 GB, any arch | 2.6.0+cu124 | 12.4 | 3.12.12 |

**Host GPU / driver (identical across all mirrors measured):** NVIDIA **A800-SXM4-80GB**,
driver **550.54.14**, host CUDA **12.4** (from `nvidia-smi`). Every peak above was
measured on this same A800-80GB / 550.54.14 / CUDA-12.4 host, so the numbers are directly
comparable.

**Checkpoint sizes (on the shared mount):** EpiAgent 5.5 GB · GET 979 MB · ATACformer
335 MB · ChromFound 1.9 MB · ChromBERT hg38-6k (cached).

## One-line recommendation

A **single 16 GB Ampere-or-newer GPU** (e.g. A4000/A5000/A100/A800) runs **every** model
in the stack for two-state scoring. **EpiAgent is the sizing driver** (~7.5 GB, and it
wants Ampere for flash-attn); the other four all fit under ~2.5 GB and are
architecture-agnostic. The A800-80GB the mirrors currently use is generous headroom, not
a requirement.

## How measured

Each model's real embed script was run on the A800 while `nvidia-smi
--query-gpu=memory.used` was sampled at 5 Hz; peak − idle (idle = 0 MiB on these dedicated
boxes) is reported. ChromBERT/ChromFound on their own mirrors (ports 35963 / 38824); GET,
EpiAgent, ATACformer via their `/yutiancheng/yuhao/miniconda3` envs (a shared mount) run on
the ChromFound A800 box because the models-zoo box (38524) was powered down — same GPU
class, driver, and CUDA host, so the figures are consistent with the other two.

## Worth noting (gaps / caveats)

1. **Inference only, default batch.** These are peak *inference* memory at each script's
   default batch on the input sizes shown. Fine-tuning/training would need substantially
   more. EpiAgent's peak is essentially **fixed** (its 8,190-cCRE rank cap bounds the
   forward regardless of how many cCREs are accessible); the coordinate/region models
   (GET/ChromFound/ATACformer) grew only mildly from tens-of-thousands to 167k regions, so
   their peaks are stable.
2. **EpiAgent env is the odd one out** — torch **2.0.1 / cu117** while the others are
   2.2–2.6 / cu12. It runs on the 12.4 host driver via CUDA backward-compat, but if you
   ever standardize envs, this is the one to re-pin. Its flash-attn path needs **Ampere
   (sm_80+)**; on a V100 (sm_70) it uses the non-flash fallback (same weights, slower,
   similar memory).
3. **GPU allocation varies per session.** Today all three boxes are A800-80GB; earlier
   notes had ChromBERT on an A100. The driver/CUDA are tied to whatever the platform
   allocates, not fixed — re-check `nvidia-smi` after a restart.
4. **Non-GPU costs not in the table:** EpiAgent loads a 5.5 GB checkpoint (needs ≥~8–12 GB
   host RAM); ATACformer's `gtars` tokenization is CPU-bound and slow (minutes for ~66k
   peaks) — a real wall-clock cost with zero GPU use.
5. **AlphaGenome (candidate 6th model) is the real GPU sizing driver** and is **not** in
   this table (not integrated). 450M params but 1 Mb-context activations want a
   large/Hopper GPU — see [`alphagenome_pipeline.md`](alphagenome_pipeline.md).
6. **cuDNN / exact NGC tag** not captured per env; ChromBERT's torch is an NVIDIA NGC
   24.04 container build (`2.3.0a0+…nv24.04`), not a stock wheel — note this if
   reproducing its env exactly.
