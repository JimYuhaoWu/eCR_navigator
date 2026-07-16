"""Self-tests for ecr_navigator/nominate.py (the run-bundle producer).

Covers the two invariants of docs/run_bundle_contract.md:
  * the PRIMARY score is resolved here (GET vs signed-Delta) and named only in the manifest;
  * REFUSAL IS AN OUTPUT -- a Gate-1 reject yields a VALID bundle with zero nominations.
Plus the contract details that bit us building the fixtures: rank ordering under a
percentile score (ties), and unmeasured direction never becoming 0.0.
"""
import json
import tempfile
from pathlib import Path

from _runner import run, add_repo_paths

add_repo_paths()

from ecr_navigator.nominate import (nominate, read_nominations, write_bundle,
                                    write_nominations)
from ecr_navigator.weights import RegionWeight

ADMIT = {"admit": True, "pc1_frac": 0.89, "coherence_margin": 0.31, "n_a": 2, "n_b": 2,
         "reasons": []}
REJECT = {"admit": False, "pc1_frac": 0.798, "coherence_margin": 0.80, "n_a": 2, "n_b": 2,
          "reasons": ["transition axis weak (PC1 0.798 < 0.8)"]}
GATE2_GET = {"primary": "GET", "anchor_set": "promoter", "n_anchors": 85,
             "delta_auroc": 0.170, "delta_auroc_ci": [0.081, 0.264], "lr_p": 0.001,
             "lr_coef": 0.67}
GATE2_SIGNED = {"primary": "signed-Delta", "anchor_set": "promoter", "n_anchors": 69,
                "delta_auroc": 0.071, "delta_auroc_ci": [-0.045, 0.188], "lr_p": 0.13,
                "lr_coef": 0.26}


def _rows(n=1000):
    """driver_score ascending, direction alternating sign with |delta| ascending -- so the
    GET and signed-Delta rankings are DIFFERENT and a test can tell which one was used."""
    out = []
    for i in range(n):
        out.append(RegionWeight(chrom="chr1", start=i * 100, end=i * 100 + 50,
                                driver_score=i / (n - 1),
                                direction=(1 if i % 2 else -1) * (n - 1 - i) / (n - 1)))
    return out


# ------------------------------------------------------------------ refusal is an output
def test_refuses_on_gate1_reject():
    noms, block = nominate(_rows(), REJECT, GATE2_GET)
    assert noms == []
    assert block["refused"] is True
    assert block["n_nominated"] == 0
    assert block["score_source"] is None
    assert "Gate-1 REJECT" in block["refusal_reason"]
    assert "PC1" in block["refusal_reason"], "the reason must carry WHY, not just that"


def test_refuses_when_gate1_missing():
    noms, block = nominate(_rows(), None, GATE2_GET)
    assert noms == [] and block["refused"] is True
    assert "Gate-1 did not run" in block["refusal_reason"]


def test_refuses_when_gate2_missing_no_fallback():
    """No anchors -> no PRIMARY -> no nominations. Never silently fall back to GET top-1%."""
    noms, block = nominate(_rows(), ADMIT, None)
    assert noms == [] and block["refused"] is True
    assert "Gate-2 did not run" in block["refusal_reason"]


def test_refusal_bundle_is_valid_and_readable():
    with tempfile.TemporaryDirectory() as d:
        noms, block = nominate(_rows(), REJECT, GATE2_GET)
        write_bundle(d, {"nomination": block}, _rows(), noms)
        p = Path(d)
        assert (p / "manifest.json").exists()
        assert (p / "nominations.tsv").exists()
        assert (p / "weights.tsv").exists(), "weights are emitted even on a refusal"
        assert read_nominations(p / "nominations.tsv") == []
        assert json.load(open(p / "manifest.json"))["nomination"]["refused"] is True


# ------------------------------------------------------------------ PRIMARY selection
def test_get_primary_nominates_top_driver_scores():
    rows = _rows(1000)
    noms, block = nominate(rows, ADMIT, GATE2_GET, top_frac=0.01)
    assert block["score_source"] == "GET"
    assert block["refused"] is False
    assert len(noms) == 10 == block["n_nominated"]
    # driver_score ascends with index -> the last rows win
    assert noms[0].nomination_score == 1.0
    assert [n.start for n in noms] == [(999 - i) * 100 for i in range(10)]


def test_signed_primary_nominates_by_measured_delta_not_driver():
    """PRIMARY=signed-Delta must rank by |direction|, NOT driver_score -- the fixture rows
    order those oppositely, so this fails loudly if the wrong score is used."""
    rows = _rows(1000)
    noms, block = nominate(rows, ADMIT, GATE2_SIGNED, top_frac=0.01)
    assert block["score_source"] == "signed-Delta"
    assert block["score_norm"] == "rank", "|delta| is rank-normalized into [0,1]"
    # |direction| DEscends with index -> the first rows win (opposite of the GET case)
    assert [n.start for n in noms] == [i * 100 for i in range(10)]


def test_signed_primary_drops_unmeasured_regions():
    rows = _rows(100)
    for r in rows[:50]:
        r.direction = None          # unmeasured -> no signed-Delta score exists
    noms, _ = nominate(rows, ADMIT, GATE2_SIGNED, top_frac=0.5)
    assert all(n.direction is not None for n in noms)
    assert len(noms) == 25, "top_frac applies to the 50 MEASURED regions, not all 100"


# ------------------------------------------------------------------ rank / direction detail
def test_rank_is_dense_and_deterministic_under_ties():
    """driver_score is a percentile, so a top_frac band is full of ties (17 at 1.0 in the
    real iN fixture). Ranks must still be 1..k, stable, and tie-broken by input order."""
    rows = [RegionWeight("chr1", i * 10, i * 10 + 5, 1.0, 0.5) for i in range(20)]
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=0.5)
    assert [n.rank for n in noms] == list(range(1, 11))
    assert [n.start for n in noms] == [i * 10 for i in range(10)], "stable: input order"


def test_unmeasured_direction_stays_empty_never_zero():
    """An unmeasured region must not be reported as measured-and-flat (0.0)."""
    rows = _rows(10)
    rows[-1].direction = None
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=1.0)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "nominations.tsv"
        write_nominations(noms, p)
        cells = [ln.split("\t") for ln in p.read_text().splitlines()[1:]]
        assert cells[0][5] == "", "unmeasured -> empty field"
        assert read_nominations(p)[0].direction is None


def test_roundtrip_preserves_rows():
    """Values survive a write/read cycle to the contract's 4dp serialization precision
    (rounding happens on write, as in weights.py -- in-memory rows keep full precision)."""
    rows = _rows(200)
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=0.05)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "nominations.tsv"
        write_nominations(noms, p)
        back = read_nominations(p)
    assert len(back) == len(noms)
    for a, b in zip(back, noms):
        assert (a.chrom, a.start, a.end, a.rank) == (b.chrom, b.start, b.end, b.rank)
        assert abs(a.nomination_score - b.nomination_score) < 1e-4
        assert abs(a.direction - b.direction) < 1e-4


if __name__ == "__main__":
    run(globals())
