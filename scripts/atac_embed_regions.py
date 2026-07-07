#!/usr/bin/env python
"""
ATACformer per-region embedding for one cell state -> embedding artifact (.npz).

ATACformer (databio, geniml.atacformer) is a transformer over ATAC region *tokens*:
a "cell" is the SET of accessible regions, drawn from a fixed universe of 890,704
regions. The state signal is *which* regions are present, so we embed one state's
accessible peak set as one bag-of-regions and read the encoder's per-region output.

Assembly note: the pretrained model `databio/atacformer-base-hg38` and its universe
are **hg38-only**. Mouse (mm10) peaks must be lifted to hg38 first (see
docs/get_pipeline.md / liftOver); pass the lifted BED here. `--assembly` is recorded
in the artifact meta only.

Region coordinates in the output are the UNIVERSE regions the peaks snap to (not the
raw peaks), so two states embedded this way share coordinates exactly and navigate.py
aligns them by chrom:start-end with no extra work.

Usage (on the atacformer mirror, atacformer env):
  python atac_embed_regions.py \
      --peaks MEF.hg38.bed --state MEF --assembly hg38 \
      --model-dir /yutiancheng/yuhao/models/atacformer \
      --out chrombert_like/atac.MEF.hg38.npz
"""
from __future__ import annotations

import argparse
import json

import numpy as np
import torch


def load(model_dir: str):
    """Return (gtars Tokenizer, AtacformerModel) from a downloaded model dir."""
    from gtars.tokenizers import Tokenizer
    from geniml.atacformer import AtacformerModel

    tok = Tokenizer.from_bed(f"{model_dir}/universe.bed.gz")
    model = AtacformerModel.from_pretrained(model_dir).eval()
    return tok, model


def tokens_for_peaks(tok, peaks_bed: str):
    """Snap peaks to universe regions -> (unique region strings, their token ids).

    A cell is a SET of regions, so duplicate hits are collapsed (first occurrence
    kept); peaks overlapping no universe region become the unk token and are dropped.
    """
    regions = tok.tokenize(peaks_bed)          # ["chr:s-e", ...] snapped to universe
    seen, uniq = set(), []
    for r in regions:
        if r not in seen:
            seen.add(r); uniq.append(r)
    ids = tok.encode(uniq)                      # parallel token ids
    unk = tok.unk_token_id
    keep = [(r, i) for r, i in zip(uniq, ids) if i != unk]
    regs = [r for r, _ in keep]
    tok_ids = [i for _, i in keep]
    return regs, tok_ids


def embed(model, token_ids: list[int], window: int, device: str) -> np.ndarray:
    """Per-region encoder embedding (N, hidden). Chunked into <=window blocks."""
    dim = model.config.hidden_size
    N = len(token_ids)
    out = np.zeros((N, dim), dtype=np.float32)
    model.to(device)
    with torch.no_grad():
        for s in range(0, N, window):
            e = min(s + window, N)
            ids = torch.tensor([token_ids[s:e]], dtype=torch.long, device=device)
            h = model(ids)                     # (1, e-s, hidden)
            out[s:e] = h[0].float().cpu().numpy()
    return out


def parse_region(r: str):
    chrom, se = r.split(":")
    start, end = se.split("-")
    return chrom, int(start), int(end)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--peaks", required=True, help="one state's accessible peaks BED (hg38)")
    ap.add_argument("--state", required=True, help="cell-state name (recorded in meta)")
    ap.add_argument("--assembly", default="hg38", help="recorded in meta; model is hg38-only")
    ap.add_argument("--model-dir", default="/yutiancheng/yuhao/models/atacformer")
    ap.add_argument("--window", type=int, default=2048,
                    help="tokens per encoder chunk (<= max_position_embeddings)")
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    tok, model = load(args.model_dir)
    regs, ids = tokens_for_peaks(tok, args.peaks)
    if not regs:
        raise SystemExit(f"no peaks in {args.peaks} snap to the hg38 universe "
                         f"(is the BED hg38? mm10 must be lifted first)")
    emb = embed(model, ids, args.window, args.device)

    coords = np.array([parse_region(r) for r in regs], dtype=object)
    chrom = coords[:, 0].astype(str)
    start = coords[:, 1].astype(np.int64)
    end = coords[:, 2].astype(np.int64)

    meta = json.dumps({"model": "atacformer", "cell_state": args.state,
                       "assembly": args.assembly, "dim": int(emb.shape[1]),
                       "source": "atac_embed_regions.py"})
    np.savez_compressed(args.out, chrom=chrom, start=start, end=end,
                        embedding=emb, meta=np.array(meta))
    print(f"wrote {args.out}: {emb.shape[0]} regions x {emb.shape[1]} dims "
          f"({args.state}, {args.assembly})")


if __name__ == "__main__":
    main()
