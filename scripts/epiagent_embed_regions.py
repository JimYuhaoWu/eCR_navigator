#!/usr/bin/env python
"""
EpiAgent per-cCRE embedding for one cell state -> embedding artifact (.npz).

EpiAgent (xy-chen16/EpiAgent) is a transformer over cCRE *tokens*: a "cell" is its
set of accessible cCREs, drawn from a fixed universe of 1,355,445 hg38 cCREs, ranked
by TF-IDF into a "cell sentence". The shipped inference mean-pools the encoder's CLS
position into a per-CELL vector; for eCR we skip that pool and keep the encoder's
per-cCRE contextual rows — one embedding artifact per state, in the same .npz
contract as GET/ChromBERT/ATACformer/ChromFound, so navigate.py diffs two states
(MEF vs mES) into driver scores unchanged.

We also emit a per-cCRE `signal` array: sigmoid(signal_decoder(CLS)), EpiAgent's
Signal-Reconstruction predicted accessibility (probability in [0,1]) sampled at the
emitted cCREs. navigate.py differences the two states' signals into the contract's
signed `direction` column (open/close) — a model-native readout, not a synthesized
one. Embedding drives magnitude; signal drives direction (decoupled channels).

Token layout (epiagent.dataset.CellDataset): input_ids = [CLS=1] + sentence[:L] +
[SEP=2], where each cCRE token id = (cCRE row index in cCRE.bed) + 4, and the
sentence is TF-IDF-sorted. So transformer_outputs[0, 1:1+L, :] are the per-cCRE
contextual embeddings, and cell sentence token t maps to var row t-4 — whose
coordinate is adata.var_names[t-4] ("chrom:start-end"). Because both states draw
cCREs from the SAME universe, the two artifacts share coordinates and navigate.py
aligns them by chrom:start-end with no extra work.

Assembly note: EpiAgent is **hg38-only**. Mouse (mm10) peaks are lifted to hg38 in
epiagent_build_input.py before this step, so the input h5ad — and thus these
coordinates — are already hg38. `--assembly` is recorded in the artifact meta only.

Coverage limit: EpiAgent's positional (rank) embedding caps a cell sentence at
max_length-2 = 8190 cCREs. A bulk pseudobulk state has far more accessible cCREs
than that, so only the top-8190 by TF-IDF are embedded per state (deterministic;
is_random NOT used). navigate.py then scores the intersection of the two states'
top-8190 sets — see docs/epiagent_pipeline.md.

Usage (EpiAgent env, on the model-zoo mirror):
  python epiagent_embed_regions.py \
      --h5ad /yutiancheng/yuhao/eCR/artifacts/epiagent.MEF.h5ad --state MEF \
      --ckpt /yutiancheng/yuhao/models/EpiAgent/pretrained_EpiAgent.pth \
      --assembly hg38 --out epiagent.MEF.hg38.npz
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch

from embedding_artifact import write_embedding_artifact

CLS, SEP = 1, 2          # EpiAgent special token ids ([CLS], [SEP])
N_SPECIAL = 4            # tokens 0..3 are special; cCRE tokens start at 4
MAX_RANK = 8192         # EpiAgent max_rank_embeddings (positional cap)


def parse_region(v: str):
    """'chr1:9848-10355' -> ('chr1', 9848, 10355)."""
    chrom, se = v.split(":")
    start, end = se.split("-")
    return chrom, int(start), int(end)


def load_model(ckpt: str, use_flash_attn: bool, device: str):
    from epiagent.model import EpiAgent
    model = EpiAgent(
        vocab_size=1355449, num_layers=18, embedding_dim=512,
        num_attention_heads=8, max_rank_embeddings=MAX_RANK,
        use_flash_attn=use_flash_attn,
        pos_weight_for_RLM=torch.tensor(1.), pos_weight_for_CCA=torch.tensor(1.),
    )
    model.load_state_dict(torch.load(ckpt, map_location="cpu"))
    return model.to(device).eval()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="one state's tokenized input (epiagent_build_input.py)")
    ap.add_argument("--state", required=True, help="cell-state name (recorded in meta)")
    ap.add_argument("--ckpt", required=True, help="pretrained_EpiAgent.pth")
    ap.add_argument("--assembly", default="hg38", help="recorded in meta; model is hg38-only")
    ap.add_argument("--use-flash-attn", choices=["auto", "yes", "no"], default="auto",
                    help="flash-attn needs Ampere+ (sm_80). 'auto' picks by GPU capability")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import scanpy as sc

    if args.use_flash_attn == "auto":
        cap = torch.cuda.get_device_capability()
        use_flash = cap[0] >= 8          # Ampere (A100 sm_80) or newer
        print(f">> GPU sm_{cap[0]}{cap[1]}: use_flash_attn={use_flash}")
    else:
        use_flash = args.use_flash_attn == "yes"

    adata = sc.read_h5ad(args.h5ad)
    if adata.n_obs != 1:
        raise SystemExit(f"{args.h5ad} has {adata.n_obs} cells; expected 1 pseudobulk state")
    var_names = adata.var_names.to_numpy()

    sentence = json.loads(adata.obs["cell_sentences"].iloc[0])   # TF-IDF-sorted cCRE token ids
    tokens = sentence[:MAX_RANK - 2]                              # deterministic top-8190
    if len(sentence) > len(tokens):
        print(f">> {args.state}: {len(sentence)} accessible cCREs; keeping top {len(tokens)} "
              f"(rank cap); {len(sentence) - len(tokens)} dropped")

    model = load_model(args.ckpt, use_flash, args.device)
    input_ids = torch.tensor([[CLS] + tokens + [SEP]], dtype=torch.long, device=args.device)
    with torch.no_grad(), torch.cuda.amp.autocast():
        out = model(input_ids, return_transformer_output=True)
        cls = out["transformer_outputs"][:, 0, :]           # (1, 512) cell embedding
        # signal_decoder: cell embedding -> accessibility logit for EVERY cCRE in the
        # universe (output col j == var row j == token j+4). This is EpiAgent's own
        # Signal-Reconstruction head, so it's a model-native accessibility readout,
        # not a navigator-invented direction. sigmoid -> P(accessible) in [0,1].
        signal_all = torch.sigmoid(model.signal_decoder(cls.float()))[0].cpu().numpy()
    # positions 1..L are the cCRE tokens (0 is CLS, L+1 is SEP)
    emb = out["transformer_outputs"][0, 1:1 + len(tokens), :].float().cpu().numpy()

    # map each cCRE token back to its universe coordinate (var row = token - 4)
    rows = [t - N_SPECIAL for t in tokens]
    coords = [parse_region(var_names[r]) for r in rows]
    chrom = [c for c, _, _ in coords]
    start = [s for _, s, _ in coords]
    end = [e for _, _, e in coords]
    signal = signal_all[rows]           # predicted accessibility per emitted cCRE

    n, d = write_embedding_artifact(
        args.out, chrom, start, end, emb,
        model="epiagent", cell_state=args.state, assembly=args.assembly,
        source="epiagent_embed_regions.py", signal=signal)
    print(f"wrote {args.out}: {n} cCREs x {d} dims ({args.state}, {args.assembly}); "
          f"signal range [{signal.min():.3f}, {signal.max():.3f}]")


if __name__ == "__main__":
    main()
