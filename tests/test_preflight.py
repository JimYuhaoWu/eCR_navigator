"""Self-tests for scripts/preflight.py (the nomination preflight gates).

GATE 1 (admissibility): a clean, strongly-separated two-state matrix with tight replicates
is ADMITTED; a noise-dominated one (weak transition axis, incoherent replicates) is REJECTED.
GATE 2 (score selection): reuses the Claim-2A harness -> PRIMARY=signed-Delta when driver is
only signed-Delta+noise, PRIMARY=GET when driver carries a real increment.
"""
from _runner import run, add_repo_paths

add_repo_paths()

import numpy as np
from preflight import admissibility, select_score


def _matrix(effect, noise, R=1500, seed=0):
    """(counts, is_state_a) for 3+3 replicates. `effect` scales the per-region A-vs-B
    contrast (the transition axis); `noise` scales per-sample scatter."""
    rng = np.random.default_rng(seed)
    mu = rng.uniform(2.0, 6.0, R)              # per-region baseline (removed by centring)
    d = rng.normal(0.0, 1.0, R)                # per-region transition direction/magnitude
    signs = np.array([+1, +1, +1, -1, -1, -1], dtype=float)
    is_a = np.array([True, True, True, False, False, False])
    L = (mu[:, None] + effect * d[:, None] * signs[None, :]
         + rng.normal(0.0, noise, (R, len(signs))))
    return np.exp(L), is_a                      # exp -> positive "counts"; loader log1p-CPMs


# ------------------------------------------------------------------ GATE 1
def test_admissibility_admits_clean_transition():
    m, is_a = _matrix(effect=1.5, noise=0.05, seed=1)
    adm = admissibility(m, is_a)
    assert adm.admit, adm.reasons
    assert adm.pc1_frac >= 0.80, adm.pc1_frac
    assert adm.coherence_margin >= 0.10, adm.coherence_margin


def test_admissibility_rejects_weak_transition():
    m, is_a = _matrix(effect=0.03, noise=0.9, seed=2)
    adm = admissibility(m, is_a)
    assert not adm.admit, (adm.pc1_frac, adm.coherence_margin)
    assert adm.pc1_frac < 0.80, adm.pc1_frac


def test_admissibility_rejects_too_few_replicates():
    m, is_a = _matrix(effect=1.5, noise=0.05, seed=3)
    is_a = is_a.copy(); is_a[:] = [True, False, False, False, False, False]  # 1 vs 5
    adm = admissibility(m, is_a)
    assert not adm.admit
    assert any("replicate" in r for r in adm.reasons), adm.reasons


# ------------------------------------------------------------------ GATE 2
def _scores(extra, n=9000, n_pos=700, seed=1):
    """driver/signed/labels: positives are large-opening regions (signed-Delta alone ranks
    them high); driver = signed + extra*is_pos + noise. extra=0 -> driver adds nothing."""
    rng = np.random.default_rng(seed)
    signed = rng.uniform(-1, 1, n)
    op = np.where(signed > 0)[0]
    w = signed[op] / signed[op].sum()
    pos = rng.choice(op, size=n_pos, replace=False, p=w)
    labels = np.zeros(n, dtype=bool); labels[pos] = True
    driver = signed + extra * labels.astype(float) + rng.normal(0, 0.25, n)
    return driver, signed, labels


def test_select_signed_when_driver_adds_nothing():
    driver, signed, labels = _scores(extra=0.0, seed=10)
    sel = select_score(driver, signed, labels, n_boot=400, n_perm=300, seed=11)
    assert sel.primary == "signed-Delta", (sel.primary, sel.reason)


def test_select_get_when_driver_adds_signal():
    driver, signed, labels = _scores(extra=1.2, seed=4)
    sel = select_score(driver, signed, labels, n_boot=500, n_perm=400, seed=5)
    assert sel.primary == "GET", (sel.primary, sel.reason)
    assert sel.delta_ci[0] > 0 or sel.perm_p < 0.05, sel


# ------------------------------------------------------- GATE 1 fixed universe (2026-07-16)
def _diluted(effect, noise, signal_regions=800, dead_regions=6000, seed=3):
    """A matrix where only `signal_regions` carry the transition and `dead_regions` are
    low-signal noise -- the shape of a real cCRE universe, where most regions are closed in
    both states. Diluting with those regions is what drags PC1 down."""
    rng = np.random.default_rng(seed)
    is_a = np.array([True, True, True, False, False, False])
    signs = np.array([+1, +1, +1, -1, -1, -1], dtype=float)
    mu = rng.uniform(4.0, 7.0, signal_regions)
    d = rng.normal(0.0, 1.0, signal_regions)
    hot = (mu[:, None] + effect * d[:, None] * signs[None, :]
           + rng.normal(0.0, noise, (signal_regions, 6)))
    cold = rng.uniform(0.0, 0.4, (dead_regions, 6)) + rng.normal(0, noise, (dead_regions, 6))
    return np.expm1(np.clip(np.vstack([hot, cold]), 0, None)), is_a


def test_pc1_is_universe_dependent_without_the_fix():
    """The finding that motivated the fixed universe: the SAME transition scores lower as
    low-signal regions are added, so a raw PC1 is not comparable across region universes."""
    m, is_a = _diluted(effect=1.0, noise=0.35)
    wide = admissibility(m, is_a, universe_n=len(m)).pc1_frac       # everything
    tight = admissibility(m, is_a, universe_n=800).pc1_frac         # signal regions only
    assert tight > wide + 0.05, (tight, wide)


def test_fixed_universe_makes_pc1_stable_across_universe_size():
    """With the universe fixed, padding the matrix with dead regions must NOT move the
    verdict -- that is the whole point of the convention."""
    m, is_a = _diluted(effect=1.0, noise=0.35, dead_regions=2000)
    m2, _ = _diluted(effect=1.0, noise=0.35, dead_regions=20000)
    a = admissibility(m, is_a, universe_n=800)
    b = admissibility(m2, is_a, universe_n=800)
    assert abs(a.pc1_frac - b.pc1_frac) < 0.02, (a.pc1_frac, b.pc1_frac)
    assert a.admit == b.admit


def test_universe_truncation_is_recorded():
    """A universe smaller than N is used whole, but flagged: its value is not comparable."""
    m, is_a = _matrix(effect=1.0, noise=0.2, R=1500)
    a = admissibility(m, is_a, universe_n=50000)
    assert a.universe_n == 1500 and a.universe_truncated is True
    b = admissibility(m, is_a, universe_n=1000)
    assert b.universe_n == 1000 and b.universe_truncated is False


def test_universe_ranking_is_depth_normalized():
    """Which regions count as 'most accessible' must not be decided by a deeply-sequenced
    sample -- rank on CPM, not raw signal. Scaling one column 10x must not move the gate."""
    m, is_a = _diluted(effect=1.0, noise=0.35)
    deep = m.copy()
    deep[:, 0] *= 10.0                      # one sample sequenced 10x deeper
    a = admissibility(m, is_a, universe_n=800)
    b = admissibility(deep, is_a, universe_n=800)
    assert abs(a.pc1_frac - b.pc1_frac) < 0.02, (a.pc1_frac, b.pc1_frac)


if __name__ == "__main__":
    run(globals())
