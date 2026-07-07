#!/usr/bin/env python
"""
GET per-region embedding for one cell state -> embedding artifact (.npz).

Genome-parameterized (hg38 OR mm10) — the ONLY species difference is the assembly
passed to the motif builder (which picks the hg38/mm10 Vierstra archetype file).
Runs on the GET model mirror in the `get` conda env.

region_motif is assembled to match GET's training pipeline exactly:
  - 282 motif scores from get_regionmotif_matrix.build_matrix (remote tabix),
    NORMALIZED per column by its max (RegionMotif.normalized_data, motif_scaler=1).
  - aTPM as the LAST (283rd) column, in [0,1] (quantitative_atac=True path).
Regions are tiled into consecutive `--window` (200) blocks; the encoder's per-region
output (dropping the cls token) is the 768-d embedding.

The output .npz matches docs/embedding_artifact.md, so navigate.py diffs two states
(e.g. MEF vs mES) into driver scores exactly as for ChromBERT.

Usage (on the mirror, get env):
  python get_embed_regions.py \
      --peaks union_named.bed --assembly mm10 \
      --atpm-tsv atpm_union.tsv --state MEF \
      --checkpoint /yutiancheng/yuhao/models/get/pretrain_fetal_adult/checkpoint-799.pth \
      --get-repo /yutiancheng/yuhao/get_model \
      --out chrombert_like/get.MEF.mm10.npz
"""
from __future__ import annotations

import argparse
import json
import sys

import numpy as np
import torch

from get_regionmotif_matrix import build_matrix   # same scripts/ dir


def load_get_model(checkpoint: str, get_repo: str):
    """Instantiate GETRegionPretrain and load our checkpoint (exact match)."""
    sys.path.insert(0, get_repo)          # so `import get_model.config` resolves
    from hydra.utils import instantiate
    from get_model.config.config import load_config
    from get_model import utils as U

    cfg = load_config("pretrain_tutorial")
    print("model target:", cfg.model.get("_target_", "?"))
    model = instantiate(cfg.model)        # GETRegionPretrain, input_dim 283
    ck = U.load_checkpoint(checkpoint, model_key=getattr(cfg.finetune, "model_key", "model"))
    ck = U.extract_state_dict(ck)
    try:
        ck = U.rename_state_dict(ck, cfg.finetune.rename_config)   # maps ckpt->model keys
    except Exception as e:
        print("rename skip:", e)
    res = model.load_state_dict(ck, strict=False)
    print(f"load: missing={len(res.missing_keys)} unexpected={len(res.unexpected_keys)}")
    if res.missing_keys:
        print("  missing sample:", res.missing_keys[:5])
    return model.eval()


def embed(model, region_motif: np.ndarray, window: int, device: str) -> np.ndarray:
    """Tile regions into windows, return per-region 768-d encoder embedding (N,768)."""
    store = {}

    def hook(_m, _i, o):
        store["enc"] = o[0] if isinstance(o, tuple) else o   # (1, win+1, 768)

    handle = model.encoder.register_forward_hook(hook)
    model.to(device)
    N = region_motif.shape[0]
    out = np.zeros((N, model.cfg.embed_dim if hasattr(model, "cfg") else 768), dtype=np.float32)
    with torch.no_grad():
        for s in range(0, N, window):
            e = min(s + window, N)
            x = torch.tensor(region_motif[s:e][None], device=device)         # (1,w,283)
            mask = torch.zeros(x.shape[0], x.shape[1], 1, dtype=torch.bool, device=device)
            try:
                model(x, mask)          # head-selection crashes on all-zero mask; hook already ran
            except Exception:
                pass
            enc = store["enc"]          # (1, w+1, 768), index 0 = cls
            out[s:e] = enc[0, 1:(e - s) + 1].float().cpu().numpy()
    handle.remove()
    return out


def load_peaks(bed: str):
    import pandas as pd
    return pd.read_csv(bed, sep="\t", header=None, usecols=[0, 1, 2],
                       names=["chrom", "start", "end"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peaks", required=True, help="region BED (the state-shared union)")
    ap.add_argument("--assembly", required=True, choices=["hg38", "mm10"])
    ap.add_argument("--atpm-tsv", required=True, help="chrom,start,end,atpm_<state>... table")
    ap.add_argument("--state", required=True, help="cell-state name; uses column atpm_<state>")
    ap.add_argument("--motif-names", required=True, help="canonical 282 motif names, one/line")
    ap.add_argument("--motif-file", default=None,
                    help="local tabix-indexed {assembly}.archetype_motifs.v1.0.bed.gz "
                         "(strongly recommended for full runs; default = slow remote URL)")
    ap.add_argument("--checkpoint", required=True)
    ap.add_argument("--get-repo", default="/yutiancheng/yuhao/get_model")
    ap.add_argument("--window", type=int, default=200)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    import pandas as pd
    import tempfile

    motif_names = [l.strip() for l in open(args.motif_names) if l.strip()]
    peaks = load_peaks(args.peaks)

    # 1. motif matrix (raw summed scores), then per-column max-normalize (motif_scaler=1)
    with tempfile.TemporaryDirectory() as wd:
        _, motif = build_matrix(args.peaks, args.assembly, motif_names, wd,
                                motif_file=args.motif_file)
    col_max = motif.max(axis=0)
    col_max[col_max == 0] = 1.0
    motif_norm = (motif / col_max).astype(np.float32)

    # 2. aTPM (already [0,1]) joined by region key, as the last column
    atpm = pd.read_csv(args.atpm_tsv, sep="\t")
    key = peaks["chrom"].astype(str) + ":" + peaks["start"].astype(str) + "-" + peaks["end"].astype(str)
    akey = atpm["chrom"].astype(str) + ":" + atpm["start"].astype(str) + "-" + atpm["end"].astype(str)
    col = f"atpm_{args.state}"
    if col not in atpm.columns:
        raise SystemExit(f"{col} not in {args.atpm_tsv} (have {list(atpm.columns)})")
    amap = dict(zip(akey, atpm[col]))
    atpm_vec = np.array([amap.get(k, 0.0) for k in key], dtype=np.float32).reshape(-1, 1)

    region_motif = np.concatenate([motif_norm, atpm_vec], axis=1)   # (N, 283)
    assert region_motif.shape[1] == len(motif_names) + 1

    # 3. embed
    model = load_get_model(args.checkpoint, args.get_repo)
    emb = embed(model, region_motif, args.window, args.device)

    meta = json.dumps({"model": "get", "cell_state": args.state,
                       "assembly": args.assembly, "dim": int(emb.shape[1]),
                       "source": "get_embed_regions.py"})
    np.savez_compressed(args.out, chrom=peaks["chrom"].to_numpy(),
                        start=peaks["start"].to_numpy(), end=peaks["end"].to_numpy(),
                        embedding=emb, meta=np.array(meta))
    print(f"wrote {args.out}: {emb.shape[0]} regions x {emb.shape[1]} dims "
          f"({args.state}, {args.assembly})")


if __name__ == "__main__":
    main()
