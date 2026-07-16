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

from ecr_navigator.nominate import (BUNDLE_VERSION, nominate, read_nominations,
                                    write_bundle, write_nominations)
from ecr_navigator.weights import (RegionWeight, fmt_direction, fmt_score,
                                   parse_direction, read_region_weights,
                                   write_region_weights)

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


# ------------------------------------------------------------------ zero-nomination guard
def test_small_universe_never_yields_a_silent_empty_nomination():
    """A non-refused bundle must never nominate nothing. round() would banker's-round a
    50-region universe at top_frac=0.01 down to k=0 while reporting refused=False -- a
    bundle state the contract does not allow."""
    for n in (1, 30, 50, 99):
        rows = [RegionWeight("chr1", i, i + 1, i / max(n - 1, 1), 0.1) for i in range(n)]
        noms, block = nominate(rows, ADMIT, GATE2_GET, top_frac=0.01)
        assert block["refused"] is False
        assert len(noms) >= 1, f"{n} regions -> {len(noms)} nominations with refused=False"
        assert block["n_nominated"] == len(noms)


def test_top_frac_uses_ceiling_not_rounding():
    """250 * 0.01 = 2.5, where ceil gives 3 but banker's round gives 2 -- a case that
    actually discriminates the two (1.5 would not: both give 2)."""
    rows = [RegionWeight("chr1", i, i + 1, i / 249, 0.1) for i in range(250)]
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=0.01)
    assert len(noms) == 3, len(noms)


# ------------------------------------------------------------------ rank / direction detail
def test_nomination_score_is_clamped_like_driver_score():
    """weights.py clamps driver_score to [0,1]; nominations.tsv must not be laxer -- the
    contract declares nomination_score is [0,1]."""
    rows = [RegionWeight("chr1", 0, 1, 1.7, 0.5), RegionWeight("chr1", 2, 3, -0.4, 0.5)]
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=1.0)
    with tempfile.TemporaryDirectory() as d:
        p = Path(d) / "nominations.tsv"
        write_nominations(noms, p)
        scores = [float(ln.split("\t")[4]) for ln in p.read_text().splitlines()[1:]]
    assert scores == [1.0, 0.0], scores


def test_write_bundle_stamps_bundle_version():
    with tempfile.TemporaryDirectory() as d:
        noms, block = nominate(_rows(100), ADMIT, GATE2_GET)
        write_bundle(d, {"nomination": block}, _rows(100), noms)   # caller omits it
        m = json.load(open(Path(d) / "manifest.json"))
        assert m["bundle_version"] == BUNDLE_VERSION


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


# ------------------------------------------------------- shared contract encoding
# weights.tsv and nominations.tsv must encode identically -- both go to the same consumer,
# and a drift in the unmeasured-vs-flat rule is one it silently misreads. The rules live in
# weights.py; these pin that nominate.py actually uses them.
def test_both_writers_encode_direction_identically():
    vals = [0.7412, 0.0, None, -0.6, 1.9, -2.4]          # incl. out-of-range and unmeasured
    rows = [RegionWeight("chr1", i * 10, i * 10 + 5, 0.5, v) for i, v in enumerate(vals)]
    noms, _ = nominate(rows, ADMIT, GATE2_GET, top_frac=1.0)
    with tempfile.TemporaryDirectory() as d:
        write_region_weights(rows, Path(d) / "w.tsv")
        write_nominations(noms, Path(d) / "n.tsv")
        wcol = [ln.split("\t")[4] for ln in (Path(d) / "w.tsv").read_text().splitlines()[1:]]
        ncol = [ln.split("\t")[5] for ln in (Path(d) / "n.tsv").read_text().splitlines()[1:]]
    assert wcol == ncol, (wcol, ncol)
    assert wcol[2] == "", "unmeasured must be empty in BOTH files"
    assert wcol[1] == "0.0", "measured-flat must stay 0.0, distinct from unmeasured"
    assert wcol[4] == "1.0" and wcol[5] == "-1.0", "both writers must clamp to [-1,1]"


def test_score_encoding_is_shared_and_clamped():
    assert fmt_score(1.7) == 1.0 and fmt_score(-0.4) == 0.0
    assert fmt_score(0.123456) == 0.1235


def test_direction_encoding_round_trips_through_the_shared_rules():
    for v in (0.5, 0.0, -1.0, None):
        assert parse_direction(fmt_direction(v)) == v
    assert parse_direction(fmt_direction(None)) is None, "unmeasured survives a round trip"


# ------------------------------------------------------- shipped boundary fixtures
# examples/run_bundle/ is eCR_predictor's only navigator-produced test material. These pin
# the guarantees its Tier-2 adapter and target.py are written against.
FIXTURES = Path(__file__).resolve().parent.parent / "examples" / "run_bundle"


def test_every_bundle_ships_weights_including_the_refusal():
    """The contract says weights.tsv is emitted for EVERY run, including refusals -- a
    Gate-1 reject is precisely when it is the only useful output."""
    for run_id in ("in_gse299923", "myod_gse186271"):
        w = read_region_weights(FIXTURES / run_id / "weights.tsv")
        assert len(w) > 1000, run_id
    myod = FIXTURES / "myod_gse186271"
    assert read_nominations(myod / "nominations.tsv") == [], "myod is the refusal bundle"
    assert len(read_region_weights(myod / "weights.tsv")) > 1000, \
        "a refusal must still ship usable off-target weights"


def test_fixture_nominations_are_a_subset_of_weights():
    """The mapping a Tier-2 adapter performs: every nominated region must be findable in
    the dense table."""
    d = FIXTURES / "in_gse299923"
    w = {(r.chrom, r.start, r.end) for r in read_region_weights(d / "weights.tsv")}
    noms = read_nominations(d / "nominations.tsv")
    assert noms and all((n.chrom, n.start, n.end) in w for n in noms)


def test_edge_case_fixture_keeps_unmeasured_distinct_from_flat():
    """The contract-critical distinction, pinned as bytes: empty field = unmeasured,
    0.0 = measured-and-flat. A consumer reading our writer must see both."""
    rows = read_region_weights(FIXTURES / "_encoding_edge_cases" / "weights.tsv")
    assert any(r.direction is None for r in rows), "no unmeasured row"
    assert any(r.direction == 0.0 for r in rows), "no measured-flat row"
    raw = (FIXTURES / "_encoding_edge_cases" / "weights.tsv").read_text().splitlines()
    assert any(ln.endswith("\t") for ln in raw[1:]), "unmeasured must serialize as EMPTY"


def test_edge_case_fixture_covers_the_four_column_form():
    """Every ChromBERT contract in the benchmark is 4-column: a missing `direction` column
    is distinct from an empty field, and consumers must handle it."""
    p = FIXTURES / "_encoding_edge_cases" / "weights.no_direction.tsv"
    assert p.read_text().splitlines()[0] == "chrom\tstart\tend\tdriver_score"
    assert all(r.direction is None for r in read_region_weights(p))


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
