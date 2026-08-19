"""Regression tests for AnnDataBuilder._compute_vst_ (roadmap.md Phase 1: "VST hanging/failing on large datasets")."""
import signal
import time

import numpy as np
import pandas as pd
import pytest

from trnagraph.modules.adataBuild import AnnDataBuilder

N_VAR = 1140  # matches the real 76 Sprinzl positions x 15 coverage types column count
ALARM_S = 10  # worst-case bound for the timed test; fixed code finishes in well under 1s
FAST_THRESHOLD_S = 5.0  # generous pass bar, still far below ALARM_S and further below a real hang


def _make_counts(n_obs, n_var, seed=0):
    rng = np.random.default_rng(seed)
    raw = rng.negative_binomial(n=2, p=0.3, size=(n_obs, n_var)).astype(float)
    obs_index = pd.Index([f"sample{i}" for i in range(n_obs)])
    var_index = pd.Index([f"var{j}" for j in range(n_var)])
    group = pd.Series(["A" if i % 2 == 0 else "B" for i in range(n_obs)], index=obs_index)
    sizefactors = rng.uniform(0.7, 1.3, size=n_obs)
    x_norm = raw / sizefactors[:, None]
    return x_norm, raw, sizefactors, group, obs_index, var_index


def test_compute_vst_stays_fast_at_moderate_sample_count():
    """
    PyDESeq2's vst_fit() only skips its own fit_size_factors() recompute
    when BOTH obsm['size_factors'] is set AND self.logmeans is not None.
    Without also setting logmeans, the pre-computed tRNA-control size
    factors _compute_vst_ passed in were silently discarded, and because
    tRNA coverage data is zero-heavy, PyDESeq2 fell back to its "iterative"
    size-factor method -- a scipy.optimize Powell-method search over one
    parameter per sample, whose cost explodes non-linearly with sample
    count. n_obs=100 used to hang (>30s and climbing) with the bug present;
    n_obs=150 here (comfortably past that threshold) should stay under a
    second with the fix in place.

    Isolating this in a subprocess with a hard kill was tried and dropped:
    in this environment, running PyDESeq2/numpy code inside ANY
    multiprocessing child (spawn, fork, or forkserver alike) deadlocks
    regardless of the VST bug, which made that approach actively less
    reliable than not isolating at all. A plain in-process signal alarm
    around the call still bounds worst-case runtime, though: if it fires
    mid-computation, that exception propagates up through PyDESeq2 and gets
    caught by _compute_vst_'s own broad `except Exception`, which falls
    back to a fast log1p path and returns anyway -- so elapsed time stays
    bounded near ALARM_S even on a regression, and the assertion below
    still catches it (just not as instantly as a hard kill would).
    """
    n_obs = 150
    x_norm, x_raw, sizefactors, group, obs_index, var_index = _make_counts(n_obs, N_VAR)
    builder = AnnDataBuilder.__new__(AnnDataBuilder)

    def _handler(signum, frame):
        raise TimeoutError(f"_compute_vst_ exceeded the {ALARM_S}s alarm bound")

    old_handler = signal.signal(signal.SIGALRM, _handler)
    signal.alarm(ALARM_S)
    t0 = time.time()
    try:
        result = builder._compute_vst_(x_norm, x_raw, sizefactors, group, "vst", "parametric", obs_index, var_index)
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, old_handler)
    elapsed = time.time() - t0

    assert elapsed < FAST_THRESHOLD_S, (
        f"_compute_vst_ took {elapsed:.1f}s for n_obs={n_obs} -- should be well under a "
        f"second; this is a regression of the VST hang (see docstring)."
    )
    assert np.asarray(result).shape == (n_obs, N_VAR)
    assert np.isfinite(result).all()


def test_compute_vst_uses_provided_size_factors():
    """
    _compute_vst_ must actually use the tRNA-control size factors it's
    given, not silently recompute its own via PyDESeq2's internal
    fit_size_factors(). Two very different size-factor vectors applied to
    the SAME counts should yield different VST output.
    """
    n_obs = 40
    x_norm, x_raw, _, group, obs_index, var_index = _make_counts(n_obs, N_VAR)
    builder = AnnDataBuilder.__new__(AnnDataBuilder)

    sf_uniform = np.ones(n_obs)
    sf_skewed = np.concatenate([np.full(n_obs // 2, 0.3), np.full(n_obs - n_obs // 2, 3.0)])

    result_uniform = np.asarray(
        builder._compute_vst_(x_norm, x_raw, sf_uniform, group, "vst", "parametric", obs_index, var_index)
    )
    result_skewed = np.asarray(
        builder._compute_vst_(x_norm, x_raw, sf_skewed, group, "vst", "parametric", obs_index, var_index)
    )

    assert not np.allclose(result_uniform, result_skewed), (
        "VST output did not change when very different size factors were supplied -- "
        "the pre-computed size factors are being silently discarded and recomputed internally."
    )
