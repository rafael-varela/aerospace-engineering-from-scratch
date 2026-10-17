"""Public tests for P02.19 Random Processes & Power Spectra. They assert properties, with scipy.signal as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Every random quantity is seeded. Statistical tolerances are derived in a comment next to each one
from the record length and the correlation time of the process; they are set at four or more
standard deviations so that a correct implementation fails by chance less than once in 10^4 runs.
"""
from functools import lru_cache

import numpy as np
import numpy.testing as npt
import pytest
from scipy import signal
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

SW, SU, L, V = 1.5, 1.5, 300.0, 60.0        # sigma_w, sigma_u [m/s], scale length [m], airspeed [m/s]
TAU = 0.35                                  # heave time constant [s] (the lesson aircraft gives 0.3501 s)
DT, DURATION, SEED = 0.02, 2000.0, 2019     # the record used by the statistical tests


@lru_cache(maxsize=None)
def gust_record() -> tuple[np.ndarray, np.ndarray]:
    """The 2000 s synthetic w record (generated once per test session, inside the tests rather than at import)."""
    t, w = impl.generate_gust(DURATION, DT, SW, L, V, "w", seed=SEED)
    return t, w


# ---------------------------------------------------------------- the Dryden model

def test_dryden_psd_integrates_to_the_variance_for_both_components():
    for comp in ("u", "w", "v"):
        var, _ = quad(lambda om: float(impl.dryden_psd(om, SW, L, V, comp)), 0.0, np.inf, limit=500)
        assert var == pytest.approx(SW**2, rel=1e-8)
        assert impl.integrated_variance(lambda om: float(impl.dryden_psd(om, SW, L, V, comp))) == pytest.approx(SW**2, rel=1e-8)
    # Doubling sigma quadruples the spectrum; doubling L doubles the low-frequency level and halves the break frequency.
    om = np.array([0.0, 0.1, 1.0])
    npt.assert_allclose(impl.dryden_psd(om, 2 * SW, L, V), 4.0 * impl.dryden_psd(om, SW, L, V))
    assert impl.dryden_psd(0.0, SW, 2 * L, V) == pytest.approx(2.0 * impl.dryden_psd(0.0, SW, L, V))
    for bad in (dict(sigma=-1.0, L=L, V=V), dict(sigma=SW, L=0.0, V=V), dict(sigma=SW, L=L, V=-5.0)):
        with pytest.raises(ValueError):
            impl.dryden_psd(om, **bad)
    with pytest.raises(ValueError):
        impl.dryden_psd(om, SW, L, V, "x")


def test_dryden_psd_limits_zero_frequency_and_inverse_square_tail():
    b = V / L
    assert impl.dryden_psd(0.0, SW, L, V, "w") == pytest.approx(SW**2 * L / (np.pi * V))
    assert impl.dryden_psd(0.0, SU, L, V, "u") == pytest.approx(2.0 * SU**2 * L / (np.pi * V))
    # Far above the break frequency both spectra fall as omega**-2; w carries 1.5 times the energy of u there.
    om = np.array([100.0, 200.0, 400.0]) * b
    pw, pu = impl.dryden_psd(om, SW, L, V, "w"), impl.dryden_psd(om, SU, L, V, "u")
    npt.assert_allclose(pw[1:] / pw[:-1], 0.25, rtol=2e-4)
    npt.assert_allclose(pw * om**2, 3.0 * SW**2 * V / (np.pi * L), rtol=3e-4)
    npt.assert_allclose(pu * om**2, 2.0 * SU**2 * V / (np.pi * L), rtol=3e-4)
    npt.assert_allclose(pw / pu, 1.5, rtol=2e-4)
    # The w spectrum is not monotone: it rises slightly above its zero-frequency value before falling.
    om_fine = np.linspace(0.0, 3.0 * b, 2000)
    assert impl.dryden_psd(om_fine, SW, L, V, "w").max() > impl.dryden_psd(0.0, SW, L, V, "w")
    assert np.all(np.diff(impl.dryden_psd(om_fine, SU, L, V, "u")) <= 0.0)


def test_psd_in_hertz_is_the_same_spectrum_rescaled():
    f = np.linspace(0.0, 25.0, 50001)
    G = impl.dryden_psd_hz(f, SW, L, V, "w")
    npt.assert_allclose(G, 2.0 * np.pi * impl.dryden_psd(2.0 * np.pi * f, SW, L, V, "w"))
    # integral G df up to 25 Hz misses only the omega**-2 tail above it: 3 sigma^2 V / (pi L omega_N).
    tail = 3.0 * SW**2 * V / (np.pi * L * (2.0 * np.pi * 25.0))
    assert np.trapezoid(G, f) == pytest.approx(SW**2 - tail, rel=2e-3)


def test_wiener_khinchin_autocorrelation_transforms_to_the_psd():
    for comp in ("u", "w"):
        R = lambda s, c=comp: float(impl.dryden_autocorrelation(s, SW, L, V, c))
        assert R(0.0) == pytest.approx(SW**2)
        for om in (0.0, 0.05, 0.2, 1.0, 3.0):
            assert impl.psd_from_autocorrelation(om, R, 40.0 * L / V) == pytest.approx(
                float(impl.dryden_psd(om, SW, L, V, comp)), rel=1e-6, abs=1e-9)
    tau = np.linspace(0.0, 12.0 * L / V, 6001)
    Ru, Rw = impl.dryden_autocorrelation(tau, SU, L, V, "u"), impl.dryden_autocorrelation(tau, SW, L, V, "w")
    assert np.all(Ru > 0.0) and np.all(np.diff(Ru) < 0.0)                       # exponential: positive, decreasing
    assert impl.dryden_autocorrelation(2.0 * L / V, SW, L, V, "w") == pytest.approx(0.0, abs=1e-15)
    assert Rw.min() < 0.0                                                          # the negative lobe of the w correlation
    npt.assert_allclose(impl.dryden_autocorrelation(-tau, SW, L, V, "w"), Rw)     # even function
    # The integral time scale of w is L/(2V), half that of u (truncating at 12 L/V costs 3e-5 relative).
    assert np.trapezoid(Rw, tau) / SW**2 == pytest.approx(L / (2.0 * V), rel=1e-3)
    assert np.trapezoid(Ru, tau) / SU**2 == pytest.approx(L / V, rel=1e-3)


# ---------------------------------------------------------------- shaping filter and simulation

def test_shaping_filter_gain_squared_is_pi_times_the_psd():
    om = np.array([0.0, 0.05, 0.2, 1.0, 5.0, 50.0])
    for comp in ("u", "w"):
        A, B, C, D = impl.shaping_filter(SW, L, V, comp)
        n = A.shape[0]
        H = np.array([(C @ np.linalg.solve(1j * w * np.eye(n) - A, B) + D)[0, 0] for w in om])
        npt.assert_allclose(np.abs(H) ** 2, np.pi * impl.dryden_psd(om, SW, L, V, comp), rtol=1e-10)
        assert np.all(np.real(np.linalg.eigvals(A)) < 0.0)                         # stable filter
    assert impl.shaping_filter(SW, L, V, "u")[0].shape == (1, 1)
    assert impl.shaping_filter(SW, L, V, "w")[0].shape == (2, 2)
    with pytest.raises(ValueError):
        impl.shaping_filter(SW, L, V, "q")


def test_zoh_discretisation_matches_scipy_and_the_scalar_closed_form():
    b, dt = V / L, 0.02
    Ad, Bd = impl.discretize_zoh(np.array([[-b]]), np.array([[b]]), dt)
    assert Ad[0, 0] == pytest.approx(np.exp(-b * dt))
    assert Bd[0, 0] == pytest.approx(1.0 - np.exp(-b * dt))
    A, B, C, D = impl.shaping_filter(SW, L, V, "w")
    Ad, Bd = impl.discretize_zoh(A, B, dt)
    Ad_s, Bd_s, *_ = signal.cont2discrete((A, B, C, D), dt, method="zoh")
    npt.assert_allclose(Ad, Ad_s, rtol=1e-12, atol=1e-15)
    npt.assert_allclose(Bd, Bd_s, rtol=1e-12, atol=1e-15)
    assert Ad.shape == (2, 2) and Bd.shape == (2, 1)


def test_simulate_lti_matches_dlsim_and_the_first_order_step_response():
    dt, n = 0.02, 400
    A, B, C, D = impl.shaping_filter(SW, L, V, "w")
    Ad, Bd = impl.discretize_zoh(A, B, dt)
    u = np.random.default_rng(7).normal(size=n)
    y = impl.simulate_lti(Ad, Bd, C, D, u)
    _, y_s, _ = signal.dlsim((Ad, Bd, C, D, dt), u)
    npt.assert_allclose(y, y_s.ravel(), rtol=1e-10, atol=1e-12)
    # A first-order lag driven by a unit step rises as 1 - exp(-t/tau), exactly at the sample instants.
    Ad1, Bd1 = impl.discretize_zoh(np.array([[-1.0 / TAU]]), np.array([[1.0 / TAU]]), dt)
    y1 = impl.simulate_lti(Ad1, Bd1, np.array([[1.0]]), np.array([[0.0]]), np.ones(n))
    t = np.arange(n) * dt
    npt.assert_allclose(y1, 1.0 - np.exp(-t / TAU), atol=1e-12)
    # D passes the input straight through when the state is zero.
    y2 = impl.simulate_lti(Ad1, Bd1, np.array([[0.0]]), np.array([[2.0]]), u)
    npt.assert_allclose(y2, 2.0 * u)


def test_white_noise_is_seeded_and_has_variance_one_over_dt():
    n, dt = 100_000, 0.02
    e = impl.white_noise(n, dt, np.random.default_rng(3))
    # Sample variance of n independent Gaussians has relative SD sqrt(2/n) = 0.45 %; 3 % is 6.7 sigma.
    assert e.var() == pytest.approx(1.0 / dt, rel=0.03)
    assert abs(e.mean()) < 5.0 / np.sqrt(n * dt)                                   # mean has SD 1/sqrt(n dt)
    npt.assert_array_equal(e, impl.white_noise(n, dt, np.random.default_rng(3)))   # same seed, same noise
    assert not np.allclose(e[:100], impl.white_noise(n, dt, np.random.default_rng(4))[:100])
    with pytest.raises(ValueError):
        impl.white_noise(0, dt, np.random.default_rng(0))


def test_generated_gust_is_reproducible_stationary_and_has_the_model_variance():
    t, w = gust_record()
    assert t.shape == w.shape == (int(DURATION / DT),)
    assert t[0] == 0.0 and t[1] - t[0] == pytest.approx(DT)
    _, w2 = impl.generate_gust(DURATION, DT, SW, L, V, "w", seed=SEED)
    npt.assert_array_equal(w, w2)
    # Variance: relative SD of the sample variance is sqrt(5 L / (4 V T)) = 5.6 % for T = 2000 s, so 25 % is 4.5 sigma.
    assert w.var() == pytest.approx(SW**2, rel=0.25)
    # Mean: Var(mean) = sigma^2 (L/V) / T, SD = 0.075 m/s; 0.3 m/s is 4 sigma.
    assert abs(w.mean()) < 0.3
    # Stationarity: the first and second halves have the same variance within 4 sigma of their own scatter (7.9 % each).
    v1, v2 = w[: w.size // 2].var(), w[w.size // 2:].var()
    assert abs(v1 - v2) / SW**2 < 0.45
    # Gaussianity: the fraction beyond 2 sigma is 4.55 %; with about T V / L = 400 independent stretches the SD is ~1 %.
    assert np.mean(np.abs(w) > 2.0 * SW) == pytest.approx(0.0455, abs=0.04)
    # The u component is generated by the same machinery and has the same variance.
    _, u = impl.generate_gust(400.0, 0.05, SU, L, V, "u", seed=11)
    assert u.var() == pytest.approx(SU**2, rel=0.5)                                # rel SD sqrt(L/(V T)) = 11 % -> 4.5 sigma


# ---------------------------------------------------------------- estimators

def test_welch_and_periodogram_match_scipy_exactly():
    _, w = gust_record()
    x = w[:20000]
    f, P = impl.welch_psd(x, DT, nperseg=1024)
    f_s, P_s = signal.welch(x, fs=1.0 / DT, nperseg=1024)
    npt.assert_allclose(f, f_s)
    npt.assert_allclose(P, P_s, rtol=1e-10)
    f, P = impl.periodogram(x, DT)
    f_s, P_s = signal.periodogram(x, fs=1.0 / DT, window="boxcar", detrend=False)
    npt.assert_allclose(f, f_s)
    npt.assert_allclose(P, P_s, rtol=1e-10)
    # Parseval: the periodogram integrates (sum times df) to the mean square of the record.
    assert np.sum(P) * f[1] == pytest.approx(np.mean(x**2), rel=1e-10)
    with pytest.raises(ValueError):
        impl.welch_psd(x, DT, nperseg=x.size + 1)


def test_welch_estimate_matches_the_dryden_spectrum_band_by_band():
    _, w = gust_record()
    f, P = impl.welch_psd(w, DT, nperseg=8192)
    G = impl.dryden_psd_hz(f, SW, L, V, "w")
    # 23 Hann segments give each bin a relative SD of about 0.2. Averaging the ~15 bins of the lowest
    # band (which are correlated two by two) leaves ~7 %, the 147 bins of the middle band ~3 % plus the
    # segment-to-segment scatter they share (observed 5-10 % over seeds), the 1470 bins of the top band ~1 %.
    assert impl.band_average(f, P, 0.01, 0.1) / impl.band_average(f, G, 0.01, 0.1) == pytest.approx(1.0, abs=0.30)
    assert impl.band_average(f, P, 0.1, 1.0) / impl.band_average(f, G, 0.1, 1.0) == pytest.approx(1.0, abs=0.30)
    assert impl.band_average(f, P, 1.0, 10.0) / impl.band_average(f, G, 1.0, 10.0) == pytest.approx(1.0, abs=0.10)
    # The estimate reproduces the omega**-2 slope between 1 and 10 Hz: a log-log fit gives -2 within 0.15.
    m = (f >= 1.0) & (f < 10.0)
    slope = np.polyfit(np.log(f[m]), np.log(P[m]), 1)[0]
    assert slope == pytest.approx(-2.0, abs=0.15)
    with pytest.raises(ValueError):
        impl.band_average(f, P, 100.0, 200.0)


def test_periodogram_is_inconsistent_and_welch_averaging_reduces_the_scatter():
    _, w = gust_record()
    f, P = impl.periodogram(w, DT)
    band = (f > 0.1) & (f < 1.0)
    rel = P[band] / impl.dryden_psd_hz(f[band], SW, L, V, "w")
    # Each periodogram bin is (asymptotically) an exponential variable: SD/mean = 1 whatever N.
    assert 0.7 < np.std(rel) / np.mean(rel) < 1.4
    assert np.mean(rel) == pytest.approx(1.0, abs=0.15)                           # unbiased on average: 1800 bins
    fw, Pw = impl.welch_psd(w, DT, nperseg=8192)
    bw = (fw > 0.1) & (fw < 1.0)
    relw = Pw[bw] / impl.dryden_psd_hz(fw[bw], SW, L, V, "w")
    # K ~ 23 segments: SD/mean ~ 1/sqrt(K) ~ 0.2; well below the periodogram's 1.
    assert np.std(relw) / np.mean(relw) < 0.4


def test_autocorrelation_estimate_matches_scipy_and_the_model():
    _, w = gust_record()
    x = w[:50000]
    max_lag = 1000
    R = impl.autocorrelation_estimate(x, max_lag)
    xc = x - x.mean()
    R_s = signal.correlate(xc, xc, mode="full")[xc.size - 1: xc.size + max_lag] / xc.size
    npt.assert_allclose(R, R_s, rtol=1e-10, atol=1e-12)
    assert R[0] == pytest.approx(x.var())
    assert np.all(np.abs(R[1:]) <= R[0])
    # Against the model: the estimate at lag k has SD of order sigma^2 sqrt(L/(V T)) ~ 0.08 sigma^2 for T = 1000 s;
    # 0.3 sigma^2 is comfortably outside the scatter and still tight enough to catch a wrong normalisation.
    lags = np.arange(max_lag + 1) * DT
    Rm = impl.dryden_autocorrelation(lags, SW, L, V, "w")
    assert np.max(np.abs(R - Rm)) < 0.3 * SW**2
    # The w correlation goes negative beyond 2L/V = 10 s, and the estimate shows a negative lobe there.
    assert R[lags > 1.5 * L / V].min() < 0.0
    with pytest.raises(ValueError):
        impl.autocorrelation_estimate(x, x.size)


# ---------------------------------------------------------------- aircraft heave response

def test_heave_time_constant_value_and_scaling():
    tau = impl.heave_time_constant(1000.0, 1.19, 60.0, 16.0, 5.0)
    assert tau == pytest.approx(2000.0 / (1.19 * 60.0 * 16.0 * 5.0))
    assert tau == pytest.approx(0.3501, abs=5e-4)
    assert impl.heave_time_constant(2000.0, 1.19, 60.0, 16.0, 5.0) == pytest.approx(2.0 * tau)   # heavier: slower
    assert impl.heave_time_constant(1000.0, 1.19, 120.0, 16.0, 5.0) == pytest.approx(0.5 * tau)  # faster: quicker
    with pytest.raises(ValueError):
        impl.heave_time_constant(1000.0, 0.0, 60.0, 16.0, 5.0)


def test_heave_frequency_response_is_a_first_order_lag_and_its_derivative():
    om = np.array([0.0, 1.0 / TAU, 10.0 / TAU])
    Hv, Ha = impl.heave_response(om, TAU)
    npt.assert_allclose(np.abs(Hv), [1.0, 1.0 / np.sqrt(2.0), 1.0 / np.sqrt(101.0)])
    assert np.angle(Hv[1]) == pytest.approx(-np.pi / 4.0)
    npt.assert_allclose(Ha, 1j * om * Hv)
    assert abs(Ha[2]) == pytest.approx(1.0 / TAU, rel=0.01)                        # high-frequency gain -> 1/tau


def test_response_variances_closed_form_matches_quadrature_and_limits():
    var_v, var_a = impl.heave_response_variances(SW, L, V, TAU)
    assert var_v == pytest.approx(impl.response_variance_numeric(SW, L, V, TAU, "v"), rel=1e-8)
    assert var_a == pytest.approx(impl.response_variance_numeric(SW, L, V, TAU, "a"), rel=1e-8)
    assert 0.0 < var_v < SW**2                                                      # a lag can only remove energy
    # tau -> 0: the aircraft follows the gust, var_v -> sigma_w^2; tau -> inf: it never responds.
    assert impl.heave_response_variances(SW, L, V, 1e-6)[0] == pytest.approx(SW**2, rel=1e-4)
    assert impl.heave_response_variances(SW, L, V, 1e6)[0] < 1e-5 * SW**2
    # var_v falls and the RMS acceleration (times tau) falls as the aircraft gets slower.
    taus = np.array([0.05, 0.2, 0.5, 2.0, 10.0])
    vv = np.array([impl.heave_response_variances(SW, L, V, t)[0] for t in taus])
    va = np.array([impl.heave_response_variances(SW, L, V, t)[1] for t in taus])
    assert np.all(np.diff(vv) < 0.0) and np.all(np.diff(va) < 0.0)
    # Sharp-edged gust bound: the RMS acceleration is below sigma_w / tau (the instant response to a step of one sigma).
    assert np.sqrt(var_a) < SW / TAU


def test_simulated_heave_response_matches_the_closed_form_variances():
    t, w = gust_record()
    v, a = impl.simulate_heave(t, w, TAU)
    var_v, var_a = impl.heave_response_variances(SW, L, V, TAU)
    # sigma_v inherits the gust's 5.6 % variance scatter (2.8 % in SD): 15 % is > 5 sigma.
    assert v.std() == pytest.approx(np.sqrt(var_v), rel=0.15)
    # a decorrelates in ~tau, so its variance scatter is ~1.8 % (0.9 % in SD) and the missing band above
    # 25 Hz removes 0.7 %: 8 % is > 7 sigma and still catches a zero-order-hold bias (+5 %).
    assert a.std() == pytest.approx(np.sqrt(var_a), rel=0.08)
    # a is the time derivative of v: compare with a central difference (relative RMS error well below 1).
    dv = np.gradient(v, t)
    assert np.sqrt(np.mean((a - dv) ** 2)) / a.std() < 0.5
    # Sharp-edged gust: a constant w gives a(0) = w/tau and v -> w.
    tt = np.arange(0, 5.0, DT)
    v1, a1 = impl.simulate_heave(tt, np.full(tt.size, 2.0), TAU)
    assert a1[0] == pytest.approx(2.0 / TAU)
    assert v1[-1] == pytest.approx(2.0, rel=1e-5)
    assert np.all(np.diff(v1) >= 0.0)


def test_rice_crossing_rate_matches_the_counted_crossings():
    t, w = gust_record()
    v, _ = impl.simulate_heave(t, w, TAU)
    var_v, var_a = impl.heave_response_variances(SW, L, V, TAU)
    nu0 = impl.rice_upcrossing_rate(0.0, np.sqrt(var_v), np.sqrt(var_a))
    assert nu0 == pytest.approx(np.sqrt(var_a / var_v) / (2.0 * np.pi))
    assert impl.rice_upcrossing_rate(np.sqrt(var_v), np.sqrt(var_v), np.sqrt(var_a)) == pytest.approx(nu0 * np.exp(-0.5))
    # About 300 zero upcrossings are expected in 2000 s; their count scatters by roughly sqrt(300) = 17 (6 %),
    # and the record's own variance scatter adds a few per cent: 25 % is about 4 sigma.
    assert impl.count_upcrossings(v, 0.0) == pytest.approx(nu0 * DURATION, rel=0.25)
    n1 = impl.count_upcrossings(v, np.sqrt(var_v))
    assert n1 == pytest.approx(nu0 * np.exp(-0.5) * DURATION, rel=0.25)
    # Counting on a known signal: a sine crosses its mean upward once per period.
    s = np.sin(2.0 * np.pi * np.arange(0, 10.0, 0.001) - 0.3)
    assert impl.count_upcrossings(s, 0.0) == 10
    assert impl.count_upcrossings(s, 2.0) == 0
    with pytest.raises(ValueError):
        impl.rice_upcrossing_rate(0.0, 0.0, 1.0)


def test_gaussian_exceedance_fraction_matches_the_simulated_record():
    assert impl.gaussian_exceedance_fraction(1.0, 1.0) == pytest.approx(0.3173, abs=1e-4)
    assert impl.gaussian_exceedance_fraction(2.0, 1.0) == pytest.approx(0.0455, abs=1e-4)
    assert impl.gaussian_exceedance_fraction(0.0, 1.0) == pytest.approx(1.0)
    t, w = gust_record()
    _, a = impl.simulate_heave(t, w, TAU)
    var_a = impl.heave_response_variances(SW, L, V, TAU)[1]
    level = 0.2 * 9.80665
    p = impl.gaussian_exceedance_fraction(level, np.sqrt(var_a))
    # The simulated fraction has ~T/(2 tau) ~ 3000 effective samples (binomial SD ~0.6 %) and inherits the
    # record's sigma_a scatter (a 1 % change in sigma moves p by 4 %): 25 % relative is > 4 sigma.
    assert np.mean(np.abs(a) > level) == pytest.approx(p, rel=0.25)
    with pytest.raises(ValueError):
        impl.gaussian_exceedance_fraction(1.0, 0.0)


def test_variance_estimate_scatter_follows_the_record_length_law():
    T, dt, n_seeds = 400.0, 0.05, 16
    vars_ = np.array([impl.generate_gust(T, dt, SW, L, V, "w", seed=500 + k)[1].var() for k in range(n_seeds)])
    predicted = impl.variance_estimate_relative_sd(T, L, V)
    assert predicted == pytest.approx(np.sqrt(5.0 * L / (4.0 * V * T)))
    measured = np.std(vars_) / SW**2
    # The SD of 16 samples is itself uncertain by about 1/sqrt(2*15) = 18 %, and the law holds only for
    # T >> L/V (here T = 80 L/V); a factor-of-two window is a real test of the 1/sqrt(T) scaling and the 5/4 constant.
    assert 0.5 < measured / predicted < 1.8
    assert impl.record_length_for_variance(predicted, L, V) == pytest.approx(T)
    assert impl.record_length_for_variance(0.01, L, V) == pytest.approx(4.0 * impl.record_length_for_variance(0.02, L, V))
    with pytest.raises(ValueError):
        impl.variance_estimate_relative_sd(0.0, L, V)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "gust_spectrum.png").stat().st_size > 1000
