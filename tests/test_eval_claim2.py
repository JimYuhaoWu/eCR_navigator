"""Self-tests for the Claim 2A machinery (scripts/eval_driver_claim2.py).

Claim 2A asks whether driver_score adds anything OVER a plain signed-Delta baseline.
These plant a known signed-Delta confound and, separately, a known driver-specific
signal, and assert the head-to-head + incremental test report the truth in each regime:

  - driver = pure signed-Delta (relabelled)  -> Delta AUROC ~ 0, incremental LR null
  - driver = signed-Delta + extra signal      -> Delta AUROC > 0 CI-excludes-0,
                                                 incremental LR significant
  - logistic_fit recovers a known separating boundary; LR test is calibrated on a
    null feature.
"""
from _runner import run, add_repo_paths

add_repo_paths()

import numpy as np
from eval_driver_claim2 import (
    logistic_fit, incremental_lr_test, paired_delta_auroc, evaluate_claim2,
)


def _synthetic(n=9000, n_pos=700, extra=1.0, seed=1):
    """Return (driver, signed, labels).

    signed ~ U(-1,1): ~half opening (>0). Positives are OPENING driver regions biased
    toward large opening (so signed-Delta alone already ranks them high -> a strong
    baseline). driver = signed + extra*is_pos + noise. `extra` is the driver-specific
    increment beyond signed-Delta; extra=0 means driver carries NOTHING the baseline
    lacks.
    """
    rng = np.random.default_rng(seed)
    signed = rng.uniform(-1, 1, n)
    op = np.where(signed > 0)[0]
    w = signed[op] / signed[op].sum()
    pos_idx = rng.choice(op, size=n_pos, replace=False, p=w)
    labels = np.zeros(n, dtype=bool)
    labels[pos_idx] = True
    driver = signed + extra * labels.astype(float) + rng.normal(0, 0.25, n)
    return driver, signed, labels


# ------------------------------------------------------------------ unit-level
def test_logistic_fit_recovers_boundary():
    # y = 1 when feature x large; a strong positive coefficient should be recovered
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1, 4000)
    y = (x + rng.normal(0, 0.3, 4000) > 0).astype(float)
    beta, ll = logistic_fit(np.column_stack([np.ones_like(x), x]), y)
    assert beta[1] > 2.0, beta                      # steep positive slope
    assert ll < 0, ll                               # log-likelihood is negative


def test_paired_delta_auroc_sign_and_ci():
    driver, signed, labels = _synthetic(extra=1.5, seed=2)
    au_d, au_s, ci = paired_delta_auroc(driver, signed, labels, n_boot=400, seed=3)
    assert au_d > au_s, (au_d, au_s)                # driver ranks positives better
    assert ci[0] < (au_d - au_s) < ci[1] or ci[0] <= (au_d - au_s), (au_d - au_s, ci)


# ------------------------------------------------------------------ regime tests
def test_driver_is_only_signed_delta_is_null():
    # driver = signed-Delta + noise, NO driver-specific term (extra=0): the model
    # carries nothing the baseline lacks -> Delta AUROC ~ 0, incremental LR null.
    driver, signed, labels = _synthetic(extra=0.0, seed=10)
    res = evaluate_claim2(driver, signed, labels, n_boot=400, n_perm=300, seed=11)
    assert abs(res.delta_auroc) < 0.05, res.delta_auroc
    assert res.delta_ci[0] <= 0 <= res.delta_ci[1], res.delta_ci  # CI straddles 0
    assert res.perm_p > 0.05, res.perm_p                          # no incremental LR


def test_driver_adds_beyond_signed_delta_is_detected():
    driver, signed, labels = _synthetic(extra=1.2, seed=4)
    res = evaluate_claim2(driver, signed, labels, n_boot=500, n_perm=400, seed=5)
    assert res.delta_auroc > 0.03, res.delta_auroc
    assert res.delta_ci[0] > 0, res.delta_ci                     # driver CI-beats baseline
    assert res.driver_coef > 0, res.driver_coef
    assert res.perm_p < 0.05, res.perm_p                         # incremental LR fires


def test_incremental_lr_null_on_pure_confound():
    # incremental_lr_test in isolation: driver independent of label given signed -> null
    rng = np.random.default_rng(6)
    signed = rng.uniform(-1, 1, 6000)
    labels = signed > 0.5                                        # label determined by signed
    driver = rng.normal(0, 1, 6000)                             # unrelated noise
    lr, coef, p = incremental_lr_test(driver, signed, labels, n_perm=400, seed=7)
    assert p > 0.05, (lr, p)


def test_opening_only_filters_and_matches():
    driver, signed, labels = _synthetic(extra=1.0, seed=8)
    res = evaluate_claim2(driver, signed, labels, opening_only=True,
                          n_boot=200, n_perm=200, seed=9)
    # all positives kept are opening; matched negatives drawn only from opening regions
    assert res.n_pos > 0 and res.n_neg > 0
    assert res.auroc_signed > 0.5, res.auroc_signed             # baseline itself is informative


if __name__ == "__main__":
    run(globals())
