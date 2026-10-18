"""Public tests for P08.13 Turbulence Fundamentals. They assert statistics, not magic numbers: the
Reynolds decomposition is exact; the moments of known signals are right; the autocorrelation of a
sinusoid and of white noise are what theory says; Welch's spectrum integrates to the variance and
agrees with scipy.signal.welch; and, on a SEEDED SYNTHETIC signal with a von Karman model spectrum,
the analysis recovers the known variance, integral scale and dissipation rate within the scatter
that a finite record allows (tolerances justified beside each assertion; they are about four
standard deviations of the estimator, measured over 40 seeds).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.signal import welch

from conftest import load_impl

impl = load_impl(__file__)

U, SIGMA, L, FS, N = 10.0, 0.3, 0.05, 20000.0, 2 ** 20      # 52 s record of SYNTHETIC grid turbulence


@pytest.fixture(scope="module")
def record():
    t, u = impl.synthetic_signal(N, FS, U, SIGMA, L, seed=7)
    return t, u


# ---------------------------------------------------------------- decomposition and moments

def test_reynolds_decomposition_is_exact():
    u = np.random.default_rng(0).normal(5.0, 2.0, 1000)
    Um, up = impl.reynolds_decomposition(u)
    assert Um == pytest.approx(np.mean(u))
    assert abs(np.mean(up)) < 1e-12                                   # <u'> = 0 by construction
    npt.assert_allclose(Um + up, u, rtol=0, atol=1e-12)


def test_moments_of_known_signals():
    # a square wave +-a about a mean: variance a^2, skewness 0, flatness 1
    u = 3.0 + 0.5 * np.tile([1.0, -1.0], 500)
    m = impl.moments(u)
    assert m["mean"] == pytest.approx(3.0) and m["var"] == pytest.approx(0.25)
    assert m["skew"] == pytest.approx(0.0, abs=1e-12) and m["flat"] == pytest.approx(1.0)
    # a sinusoid: flatness 1.5
    x = np.sin(2 * np.pi * np.arange(10000) / 100)
    assert impl.moments(x)["flat"] == pytest.approx(1.5, rel=1e-6)
    # a one-sided spike train is positively skewed
    assert impl.moments(np.r_[np.zeros(99), 1.0])["skew"] > 5


# ---------------------------------------------------------------- autocorrelation

def test_autocorrelation_of_sinusoid_and_white_noise():
    n = 20000
    x = np.cos(2 * np.pi * np.arange(n) / 50)
    rho = impl.autocorrelation(x - x.mean(), 100)
    assert rho[0] == 1.0
    # biased estimator: rho(k) = (1 - k/n) cos(2 pi k / 50) for a pure cosine (to O(1/n))
    k = np.arange(101)
    npt.assert_allclose(rho, (1 - k / n) * np.cos(2 * np.pi * k / 50), atol=2e-3)
    w = np.random.default_rng(3).standard_normal(n)
    rw = impl.autocorrelation(w - w.mean(), 20)
    assert np.all(np.abs(rw[1:]) < 4 / np.sqrt(n))                    # white noise: rho(k > 0) ~ N(0, 1/n)


def test_autocorrelation_is_linear_not_circular():
    # a ramp's circular correlation would wrap the end onto the start; the linear one does not
    x = np.linspace(-1, 1, 64)
    rho = impl.autocorrelation(x, 63)
    direct = np.array([np.sum(x[: 64 - k] * x[k:]) for k in range(64)]) / np.sum(x * x)
    npt.assert_allclose(rho, direct, atol=1e-12)


def test_integral_time_scale_of_an_exponential():
    dt, T = 1e-4, 5e-3
    tau = np.arange(0, 0.1, dt)
    rho = np.exp(-tau / T)
    assert impl.integral_time_scale(rho, dt) == pytest.approx(T, rel=1e-3)
    rho2 = np.r_[np.exp(-tau[:200] / T), -0.1, 0.5, 0.5]             # tail after the first zero is ignored
    assert impl.integral_time_scale(rho2, dt) < T


# ---------------------------------------------------------------- spectrum

def test_welch_integrates_to_variance_and_matches_scipy(record):
    t, u = record
    f, S = impl.welch_psd(u, FS, 2 ** 14)
    fs_, Ss = welch(u, FS, window="hann", nperseg=2 ** 14, noverlap=2 ** 13, detrend="constant")
    npt.assert_allclose(f, fs_)
    npt.assert_allclose(S, Ss, rtol=1e-8, atol=1e-14)
    # Parseval: the spectrum integrates to the variance, minus what the per-segment mean removal and
    # the window take from the lowest bins (under 1 % here at df = 1.2 Hz; about 5 % at nperseg 2^12)
    assert np.sum(S) * (f[1] - f[0]) == pytest.approx(np.var(u), rel=0.03)


def test_welch_of_a_pure_tone_puts_its_power_at_the_tone():
    fs, n = 1000.0, 2 ** 14
    x = 2.0 * np.sin(2 * np.pi * 125.0 * np.arange(n) / fs)
    f, S = impl.welch_psd(x, fs, 256)
    assert f[np.argmax(S)] == pytest.approx(125.0)
    assert np.sum(S) * (f[1] - f[0]) == pytest.approx(2.0, rel=1e-3)  # variance of 2 sin = 2


# ---------------------------------------------------------------- the synthetic record (labelled synthetic)

def test_synthetic_signal_is_seeded_and_has_the_requested_mean(record):
    t, u = record
    assert t.size == N and t[1] - t[0] == pytest.approx(1 / FS)
    _, u2 = impl.synthetic_signal(N, FS, U, SIGMA, L, seed=7)
    npt.assert_array_equal(u, u2)                                     # same seed, same record
    assert np.mean(u) == pytest.approx(U, abs=1e-9)                   # zero-frequency bin left empty


def test_recovers_the_known_variance(record):
    # relative std of the variance estimate is sqrt(2 int rho^2 / T_rec) = 1.3 % for this record
    # (int rho^2 dtau = 0.866 L/U for the von Karman model); tolerance 6 % is > 4 standard deviations
    m = impl.moments(record[1])
    bv = impl.band_variance(N, FS, U, SIGMA, L)
    assert m["var"] == pytest.approx(bv, rel=0.06)
    assert bv < SIGMA ** 2 and bv > 0.97 * SIGMA ** 2                # the sampled band holds ~99 % of sigma^2
    assert abs(m["skew"]) < 0.05 and m["flat"] == pytest.approx(3.0, abs=0.1)   # Gaussian by construction


def test_recovers_the_integral_scale_two_ways(record):
    # over 40 seeds the autocorrelation estimate of L has mean 51.4 mm and std 2.4 mm, the
    # spectral-plateau estimate 47.6 mm and 2.8 mm: 20 % tolerances are about 4 standard deviations
    t, u = record
    m = impl.moments(u)
    _, up = impl.reynolds_decomposition(u)
    T = impl.integral_time_scale(impl.autocorrelation(up, 4000), 1 / FS)
    assert U * T == pytest.approx(L, rel=0.20)                         # Taylor: L = U T
    f, S = impl.welch_psd(u, FS, 2 ** 14)
    assert impl.integral_scale_from_spectrum(f, S, m["var"], U) == pytest.approx(L, rel=0.20)


def test_recovers_minus_five_thirds_and_the_dissipation_rate(record):
    # over 40 seeds: slope std 0.004, epsilon std 0.35 %; tolerances 0.05 and 3 %
    f, S = impl.welch_psd(record[1], FS, 2 ** 14)
    slope, eps = impl.fit_inertial_range(f, S, U, 300.0, 3000.0)
    assert slope == pytest.approx(-5 / 3, abs=0.05)
    assert eps == pytest.approx(impl.model_epsilon(SIGMA, L), rel=0.03)


def test_kolmogorov_scales_are_consistent():
    nu, eps = 1.5e-5, 0.4
    k = impl.kolmogorov_scales(eps, nu)
    assert k["eta"] * k["v"] / nu == pytest.approx(1.0)                # Reynolds number of the smallest eddies is 1
    assert k["eta"] / k["tau"] == pytest.approx(k["v"])
    assert k["v"] ** 2 / k["tau"] == pytest.approx(eps)               # they dissipate eps
    # 16 times the dissipation halves eta
    assert impl.kolmogorov_scales(16 * eps, nu)["eta"] == pytest.approx(k["eta"] / 2)


def test_record_length_scales_inversely_with_error_squared():
    T1 = impl.record_length_for(5e-3, 0.05)
    assert T1 == pytest.approx(4.0)
    assert impl.record_length_for(5e-3, 0.025) == pytest.approx(4 * T1)


def test_analysis_is_unit_invariant():
    # the same synthetic record in cm/s: variance x 1e4, length scale and slope unchanged
    t, u = impl.synthetic_signal(2 ** 16, FS, U, SIGMA, L, seed=11)
    u_cm = 100 * u
    assert impl.moments(u_cm)["var"] == pytest.approx(1e4 * impl.moments(u)["var"])
    _, up = impl.reynolds_decomposition(u)
    _, upc = impl.reynolds_decomposition(u_cm)
    npt.assert_allclose(impl.autocorrelation(upc, 200), impl.autocorrelation(up, 200), atol=1e-12)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "turbulence_synthetic.png").stat().st_size > 10_000
