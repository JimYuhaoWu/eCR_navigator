"""Self-tests for the fine-tuned driver head (scripts/finetune_head.py).

The point of leave-one-GENE-out is to reward a head that learned a GENERALIZABLE
driver-region signature and to punish one that merely memorized specific loci. These
plant the two regimes in the embedding-shift space and assert the CV reports the truth:

  - SHARED driver direction across genes -> out-of-fold AUROC clearly > 0.5 (transfers)
  - GENE-SPECIFIC directions (no shared component) -> out-of-fold AUROC ~ 0.5 EVEN
    THOUGH an in-fold fit would be high (this is the whole reason for leave-one-gene-out)
  - unit checks for the loaders / PCA / gene tagging
"""
from _runner import run, add_repo_paths

add_repo_paths()

import numpy as np
from finetune_head import (
    pca_transform, leave_one_gene_out_scores, transfer_scores, gene_of_region,
    load_named_bed, join_values,
)
from eval_driver_claim1 import auroc


def _synth(shared, n=6000, n_genes=12, per_gene=25, dim=40, strength=3.0, seed=1):
    """Return (feats, labels, genes, signed) with a planted driver signal in the
    embedding SHIFT. `shared`: one driver direction for all genes (transferable) vs a
    distinct direction per gene (memorizable but not transferable). Positives are also
    biased toward high opening (signed-Δ) so the confound must be matched away."""
    rng = np.random.default_rng(seed)
    # all-opening (signed>0), mirroring finetune_head.main()'s opening-only filter, so
    # the matched signed-Δ feature cannot separate positives from negatives by SIGN —
    # the head must win on the embedding shift, not the accessibility confound.
    signed = rng.uniform(0, 1, n)
    shift = rng.normal(0, 1, (n, dim))
    labels = np.zeros(n, dtype=bool)
    genes = np.full(n, "", dtype=object)

    w = signed / signed.sum()                        # bias positives toward big opening
    chosen = rng.choice(n, size=n_genes * per_gene, replace=False, p=w)
    u_shared = rng.normal(0, 1, dim); u_shared /= np.linalg.norm(u_shared)
    for gi in range(n_genes):
        idx = chosen[gi * per_gene:(gi + 1) * per_gene]
        labels[idx] = True
        genes[idx] = f"gene{gi}"
        u = u_shared if shared else _unit(rng, dim)
        shift[idx] += strength * u
    feats = np.column_stack([pca_transform(shift, 15),
                             (signed - signed.mean()) / signed.std()])
    return feats, labels, genes, signed


def _unit(rng, dim):
    v = rng.normal(0, 1, dim)
    return v / np.linalg.norm(v)


# ------------------------------------------------------------------ unit-level
def test_pca_transform_shape_and_scale():
    rng = np.random.default_rng(0)
    z = pca_transform(rng.normal(0, 1, (2000, 30)), 8)
    assert z.shape == (2000, 8)
    assert np.allclose(z.std(axis=0), 1.0, atol=0.05)      # columns standardized


def test_gene_tagging_and_named_bed(tmp=None):
    chrom = np.array(["chr1", "chr1", "chr2"])
    start = np.array([100, 5000, 100]); end = np.array([200, 5100, 200])
    gc = np.array(["chr1", "chr2"]); gs = np.array([150, 90]); ge = np.array([250, 210])
    gn = np.array(["Sox2", "Pou5f1"])
    tag = gene_of_region(chrom, start, end, gc, gs, ge, gn)
    assert tag.tolist() == ["Sox2", "", "Pou5f1"], tag.tolist()


def test_join_values_exact_key():
    chrom = np.array(["chr1", "chr1"]); start = np.array([10, 20]); end = np.array([15, 25])
    v = join_values(chrom, start, end, np.array(["chr1"]), np.array([20]),
                    np.array([25]), np.array([0.7]))
    assert np.isnan(v[0]) and abs(v[1] - 0.7) < 1e-9, v


# ------------------------------------------------------------------ regime tests
def test_shared_signature_transfers():
    feats, labels, genes, signed = _synth(shared=True, seed=2)
    _, score, lab = leave_one_gene_out_scores(feats, labels, genes, np.abs(signed), seed=3)
    au = auroc(score, lab)
    assert au > 0.65, au                    # generalizes to genes it never trained on


def test_gene_specific_signature_does_not_transfer():
    # each gene's driver direction is unique -> a head trained on OTHER genes cannot
    # score the held-out gene above chance, even though an in-fold fit would ace it.
    feats, labels, genes, signed = _synth(shared=False, seed=4)
    _, score, lab = leave_one_gene_out_scores(feats, labels, genes, np.abs(signed), seed=5)
    au = auroc(score, lab)
    assert au < 0.60, au                    # leave-one-gene-out exposes the non-transfer


def test_leave_one_gene_out_is_truly_out_of_fold():
    # held-out positives count equals the total planted positives (each scored once,
    # by a fold that excluded its gene)
    feats, labels, genes, signed = _synth(shared=True, seed=6)
    idx, score, lab = leave_one_gene_out_scores(feats, labels, genes, np.abs(signed), seed=7)
    assert int(lab.sum()) == int(labels.sum()), (lab.sum(), labels.sum())
    assert len(idx) == len(score) == len(lab)


def test_transfer_shared_direction_generalizes_across_panels():
    # panel A and panel B positives share a driver direction -> a head trained on A
    # recovers B (held out, shared regions excluded)
    rng = np.random.default_rng(20)
    n, dim = 8000, 40
    signed = rng.uniform(0, 1, n)
    shift = rng.normal(0, 1, (n, dim))
    u = _unit(rng, dim)
    a_pos = rng.choice(n, 400, replace=False, p=signed / signed.sum())
    rest = np.setdiff1d(np.arange(n), a_pos)
    b_pos = rng.choice(rest, 400, replace=False, p=signed[rest] / signed[rest].sum())
    lab_a = np.zeros(n, bool); lab_a[a_pos] = True
    lab_b = np.zeros(n, bool); lab_b[b_pos] = True
    shift[a_pos] += 3.0 * u; shift[b_pos] += 3.0 * u        # shared driver direction
    feats = np.column_stack([pca_transform(shift, 15), (signed - signed.mean()) / signed.std()])
    _, score, lab = transfer_scores(feats, lab_a, lab_b, np.abs(signed), seed=21)
    assert auroc(score, lab) > 0.65, auroc(score, lab)


def test_transfer_panel_specific_direction_does_not_generalize():
    # panel A and panel B have DIFFERENT directions -> no cross-panel transfer
    rng = np.random.default_rng(22)
    n, dim = 8000, 40
    signed = rng.uniform(0, 1, n)
    shift = rng.normal(0, 1, (n, dim))
    a_pos = rng.choice(n, 400, replace=False, p=signed / signed.sum())
    rest = np.setdiff1d(np.arange(n), a_pos)
    b_pos = rng.choice(rest, 400, replace=False, p=signed[rest] / signed[rest].sum())
    lab_a = np.zeros(n, bool); lab_a[a_pos] = True
    lab_b = np.zeros(n, bool); lab_b[b_pos] = True
    shift[a_pos] += 3.0 * _unit(rng, dim); shift[b_pos] += 3.0 * _unit(rng, dim)
    feats = np.column_stack([pca_transform(shift, 15), (signed - signed.mean()) / signed.std()])
    _, score, lab = transfer_scores(feats, lab_a, lab_b, np.abs(signed), seed=23)
    assert auroc(score, lab) < 0.60, auroc(score, lab)


if __name__ == "__main__":
    run(globals())
