#!/usr/bin/env python
"""
ChromFound per-region (per-OCR) embedding for each cell state -> .npz artifacts.

ChromFound's encoder yields per-OCR contextual representations; the shipped
cell_embedding.py mean-pools them (axis=-1) into a cell vector. For eCR we skip
that pool and keep the per-OCR 128-d vectors — one embedding artifact per state,
in the same .npz contract as GET/ChromBERT/ATACformer, so navigate.py diffs two
states (e.g. MEF vs mES) into driver scores unchanged.

forward(value, chromosome, hg38_start, hg38_end) = backbone(embedding(...))
returns (batch, 128, n_OCR); per-OCR embedding for a cell = that[0].T (n_OCR,128).
Because DataProcessorForPad does no reorder and max_length = n_OCR (no pad/trunc),
per-OCR row i corresponds to var row i.

Run from the ChromFound repo root (needs its src on the path). GPU required.

Usage (chromfound env, on the mirror):
  cd /root/ChromFound && python chromfound_embed_regions.py \
      --h5ad chromfound_input.MEF_mES.h5ad \
      --repo /root/ChromFound \
      --ckpt-dir /yutiancheng/yuhao/models/chromFound \
      --assembly hg38 --out-dir /tmp/chromfound_emb
"""
from __future__ import annotations

import argparse
import os
import sys

import numpy as np
import torch

from embedding_artifact import write_embedding_artifact


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--h5ad", required=True, help="input built by chromfound_build_input.py")
    ap.add_argument("--repo", default="/root/ChromFound")
    ap.add_argument("--ckpt-dir", default="/yutiancheng/yuhao/models/chromFound")
    ap.add_argument("--model-file", default="model.pt")
    ap.add_argument("--config-file", default="chromfd_pretrain.yaml")
    ap.add_argument("--assembly", default="hg38", help="recorded in meta (model is hg38)")
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    # repo imports (dataset_ds uses a bare `from tf_data_process import ...`)
    for p in ["", "src", "src/data", "src/models", "src/utils"]:
        sys.path.insert(0, os.path.join(args.repo, p))
    import yaml
    import scanpy as sc
    from src.cell_embedding import EmbeddingModel
    from src.data.dataset_ds import DatasetMultiPad
    from src.utils.model_utils import ModelUtils

    device = torch.device("cuda")

    with open(os.path.join(args.ckpt_dir, args.config_file)) as f:
        cfg = yaml.safe_load(f)
    m_args, d_args = cfg["model_args"], cfg["data_args"]
    vocab = ModelUtils.get_chromosome_vocab(os.path.join(args.ckpt_dir, "chromosome_vocab.yaml"))
    d_args["chromosome_vocab"] = vocab

    adata = sc.read_h5ad(args.h5ad)
    n_ocr = adata.shape[1]
    states = list(adata.obs["celltype"])
    ct_map = {c: i for i, c in enumerate(sorted(set(states)))}

    d_args.update(cell_type_map=ct_map, cell_type_col="celltype",
                  feature_num=n_ocr, max_length=n_ocr, return_batch_label=False)
    m_args.update(cell_type_num=len(ct_map), feature_num=n_ocr, batch_size=1,
                  max_length=n_ocr, device=device, add_cls=d_args["add_cls"],
                  mask_ratio=0.0)

    model = EmbeddingModel(**m_args)
    state = torch.load(os.path.join(args.ckpt_dir, args.model_file), map_location="cpu")
    model.load_state_dict(state["module"])
    model = model.to(device).eval()

    ds = DatasetMultiPad(adata, **d_args)          # one item per state (cell)
    chrom = adata.var["#Chromosome"].to_numpy().astype(str)
    start = adata.var["hg38_Start"].to_numpy().astype(np.int64)
    end = adata.var["hg38_End"].to_numpy().astype(np.int64)

    os.makedirs(args.out_dir, exist_ok=True)
    with torch.no_grad():
        for i, st in enumerate(states):
            value, chromosome, ps, pe, _ = ds[i]
            out = model(value[None].to(device), chromosome[None].to(device),
                        ps[None].to(device), pe[None].to(device))   # (1, n_OCR, 128)
            e = out[0].float().cpu().numpy()
            if e.shape[0] != n_ocr and e.shape[1] == n_ocr:         # orient to (n_OCR, dim)
                e = e.T
            emb = e[:n_ocr]                                          # (n_OCR, 128)
            outp = os.path.join(args.out_dir, f"chromfound.{st}.{args.assembly}.npz")
            n, d = write_embedding_artifact(
                outp, chrom, start, end, emb,
                model="chromfound", cell_state=st, assembly=args.assembly,
                source="chromfound_embed_regions.py")
            print(f"wrote {outp}: {n} OCRs x {d} dims ({st})")


if __name__ == "__main__":
    main()
