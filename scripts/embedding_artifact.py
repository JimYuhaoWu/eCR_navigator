#!/usr/bin/env python
"""
Shared writer for the embedding-artifact contract (docs/embedding_artifact.md).

Every model mirror (ChromBERT, GET, ATACformer, ChromFound, and whatever comes
next) emits ONE .npz per (model, cell_state) through this single function, so the
on-disk format cannot drift per model. In particular it pins the dtypes the
navigator's loader requires under `allow_pickle=False`:

  chrom     -> unicode string array  (NOT object; a pandas .to_numpy() gives object
                                       arrays, which fail to load — a real past bug)
  start/end -> int64
  embedding -> float32
  meta      -> 0-d unicode-string array holding JSON

numpy-only (no torch/pandas) so it imports in any mirror env. Model scripts import
it as a sibling, exactly like get_embed_regions imports get_regionmotif_matrix;
upload this file alongside the embed script when running on a mirror.
"""
from __future__ import annotations

import json

import numpy as np


def write_embedding_artifact(out, chrom, start, end, embedding, *,
                             model: str, cell_state: str, assembly: str,
                             source: str, signal=None) -> tuple[int, int]:
    """Write one embedding artifact, enforcing the contract dtypes.

    chrom/start/end/embedding may be any array-like (list, pandas Series, ndarray);
    they are coerced to the contract dtypes here so no caller can get them wrong.
    Returns (n_regions, dim).

    signal (optional): a per-region SCALAR accessibility value for this state (one
    float per region, same order as `embedding`). Only models with a scalar
    per-state readout supply it — e.g. EpiAgent's signal-decoder predicted
    accessibility. navigate.py differences the two states' signals into the signed
    `direction` column of the contract; embedding-only models pass signal=None and
    the direction column is simply absent. See docs/embedding_artifact.md.
    """
    chrom = np.asarray(chrom).astype(str)          # -> '<U' unicode, never object
    start = np.asarray(start).astype(np.int64)
    end = np.asarray(end).astype(np.int64)
    embedding = np.asarray(embedding).astype(np.float32)
    n, d = embedding.shape
    if not (len(chrom) == len(start) == len(end) == n):
        raise ValueError(f"length mismatch: chrom {len(chrom)} start {len(start)} "
                         f"end {len(end)} embedding {n}")
    arrays = dict(chrom=chrom, start=start, end=end, embedding=embedding)
    if signal is not None:
        signal = np.asarray(signal).astype(np.float32)
        if len(signal) != n:
            raise ValueError(f"signal length {len(signal)} != {n} regions")
        arrays["signal"] = signal
    meta = json.dumps({"model": model, "cell_state": cell_state,
                       "assembly": assembly, "dim": int(d), "source": source,
                       "has_signal": signal is not None})
    np.savez_compressed(out, meta=np.array(meta), **arrays)
    return n, d
