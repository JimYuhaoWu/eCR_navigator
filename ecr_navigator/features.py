"""
Load embedding artifacts (docs/embedding_artifact.md) and compute the per-region
embedding shift between two cell states — the zero-shot driver signal.

numpy-only: no torch, no model deps. The heavy embedding step already ran in a
server mirror; this consumes its .npz output.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass
class EmbeddingArtifact:
    chrom: np.ndarray            # (N,) str
    start: np.ndarray            # (N,) int64
    end: np.ndarray              # (N,) int64
    embedding: np.ndarray        # (N, D) float32
    meta: dict
    signal: np.ndarray | None = None  # (N,) float32 per-region scalar accessibility,
    #                                   present only for direction-capable models

    @property
    def keys(self) -> np.ndarray:
        """(N,) array of 'chrom:start-end' region keys for alignment."""
        return np.array([f"{c}:{s}-{e}"
                         for c, s, e in zip(self.chrom, self.start, self.end)])


def load_artifact(path: str | Path) -> EmbeddingArtifact:
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    signal = z["signal"].astype(np.float32) if "signal" in z.files else None
    return EmbeddingArtifact(
        chrom=z["chrom"], start=z["start"], end=z["end"],
        embedding=z["embedding"].astype(np.float32), meta=meta, signal=signal,
    )


def _aligned_rows(a: EmbeddingArtifact, b: EmbeddingArtifact):
    """Row indices of the regions shared by both artifacts, aligned by region key."""
    idx_b = {k: i for i, k in enumerate(b.keys)}
    rows_a, rows_b = [], []
    for i, k in enumerate(a.keys):
        j = idx_b.get(k)
        if j is not None:
            rows_a.append(i)
            rows_b.append(j)
    if not rows_a:
        raise ValueError("no shared regions between the two artifacts")
    return np.array(rows_a), np.array(rows_b)


def embedding_shift(a: EmbeddingArtifact, b: EmbeddingArtifact):
    """
    Align two artifacts by region key and return (chrom, start, end, shift) for
    the regions present in both. `shift` = L2 norm of (emb_b - emb_a) per region.

    Drivers reorganize their embedding most between states; passengers move little
    (model-agnostic — whatever foundation model produced the artifacts). Assemblies
    must match.
    """
    if a.meta.get("assembly") != b.meta.get("assembly"):
        raise ValueError("assembly mismatch: %s vs %s"
                         % (a.meta.get("assembly"), b.meta.get("assembly")))
    if a.embedding.shape[1] != b.embedding.shape[1]:
        raise ValueError("embedding dim mismatch: %d vs %d"
                         % (a.embedding.shape[1], b.embedding.shape[1]))

    rows_a, rows_b = _aligned_rows(a, b)
    diff = b.embedding[rows_b] - a.embedding[rows_a]
    shift = np.linalg.norm(diff, axis=1)
    return (a.chrom[rows_a], a.start[rows_a], a.end[rows_a], shift)


def signed_delta(a: EmbeddingArtifact, b: EmbeddingArtifact):
    """
    Per-region signed change in the scalar accessibility signal, state A -> B,
    aligned to the SAME shared regions (same order) as `embedding_shift`.

    Returns `delta = signal_b - signal_a` (>0 opens, <0 closes) or None if either
    artifact lacks a signal. This feeds the contract's `direction` column; the
    magnitude driver_score still comes from `embedding_shift`. The two are decoupled
    channels and are joined by navigate.py on the identical shared-region order.
    """
    if a.signal is None or b.signal is None:
        return None
    rows_a, rows_b = _aligned_rows(a, b)
    return b.signal[rows_b] - a.signal[rows_a]
