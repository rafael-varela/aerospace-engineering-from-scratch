"""Public tests for P06.16 Random Vibration. They assert physics and statistics, not magic numbers: the
exact log-log area of a breakpoint spectrum against closed forms and quadrature; the dB/octave slopes;
the SDOF transmissibility's limits (1 at f = 0, 1 at r = sqrt 2 for every damping, about Q^2 at
resonance, falling at high frequency); Miles' equation as the light-damping limit of the full |H|^2 G
integral and the exact white-noise result (pi/2) fn Q G (1 + 4 zeta^2); the relative-displacement RMS
against its own closed form; the synthesized SYNTHETIC record's variance, mean, Gaussian tails and
spectrum; the exact discretization against sinusoidal steady state and a static input; and the
independent check of the lesson, the analytic response Grms against a seeded time-domain simulation.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

G0 = 9.80665
F_BP = np.array([20.0, 80.0, 500.0, 2000.0])       # illustrative minimum-workmanship-style spectrum
G_BP = np.array([0.01, 0.04, 0.04, 0.01])


# ---------------------------------------------------------------- the input spectrum

def test_area_of_breakpoint_spectrum_is_exact():
    # +3 dB/oct (G ~ f), flat and -3 dB/oct (G ~ 1/f, the logarithmic case) by hand
    exact = 0.01 / 20 * (80 ** 2 - 20 ** 2) / 2 + 0.04 * 420 + 0.04 * 500 * np.log(4.0)
    assert impl.psd_area(F_BP, G_BP) == pytest.approx(exact, rel=1e-12)
    assert impl.grms(F_BP, G_BP) == pytest.approx(np.sqrt(exact), rel=1e-12)
    # an arbitrary spectrum against adaptive quadrature of psd_loglog
    f_bp = np.array([10.0, 45.0, 300.0, 700.0, 1500.0])
    G_bp = np.array([0.002, 0.03, 0.03, 0.2, 0.004])
    num = sum(quad(lambda f: impl.psd_loglog(np.array([f]), f_bp, G_bp)[0], a, b, epsrel=1e-12)[0]
              for a, b in zip(f_bp[:-1], f_bp[1:]))
    assert impl.psd_area(f_bp, G_bp) == pytest.approx(num, rel=1e-9)


def test_loglog_interpolation_and_slopes():
    npt.assert_allclose(impl.psd_loglog(F_BP, F_BP, G_BP), G_BP, rtol=1e-12)
    assert impl.psd_loglog(np.array([10.0, 2500.0]), F_BP, G_BP).tolist() == [0.0, 0.0]   # zero outside
    slopes = impl.segment_slopes_db_per_octave(F_BP, G_BP)
    npt.assert_allclose(slopes, [10 * np.log10(2), 0.0, -10 * np.log10(2)], atol=1e-12)
    # one octave up a +3.01 dB/oct segment doubles G
    assert impl.psd_loglog(np.array([40.0]), F_BP, G_BP)[0] == pytest.approx(0.02, rel=1e-12)


def test_grms_scales_as_square_root_of_level():
    assert impl.grms(F_BP, 4.0 * G_BP) == pytest.approx(2.0 * impl.grms(F_BP, G_BP), rel=1e-12)


# ---------------------------------------------------------------- SDOF response

def test_transmissibility_limits():
    fn = 100.0
    for zeta in (0.01, 0.05, 0.3):
        assert impl.transmissibility_sq(np.array([0.0]), fn, zeta)[0] == pytest.approx(1.0)
        # every curve crosses 1 at r = sqrt 2: the start of isolation
        assert impl.transmissibility_sq(np.array([np.sqrt(2) * fn]), fn, zeta)[0] == pytest.approx(1.0, rel=1e-12)
        # at resonance |H|^2 = (1 + 4 zeta^2) / (4 zeta^2), about Q^2
        assert impl.transmissibility_sq(np.array([fn]), fn, zeta)[0] == pytest.approx(1 + 1 / (4 * zeta ** 2))
    # above resonance it falls; less damping isolates better at high frequency
    f = np.array([300.0, 1000.0, 3000.0])
    T1, T3 = impl.transmissibility_sq(f, fn, 0.01), impl.transmissibility_sq(f, fn, 0.3)
    assert np.all(np.diff(T1) < 0) and np.all(T1 < T3)


def test_miles_is_the_white_noise_limit_of_the_full_integral():
    # flat input from far below to far above fn: exact result (pi/2) fn Q G (1 + 4 zeta^2), minus the
    # tails outside the band (|H|^2 ~ 1 below and ~ (2 zeta fn / f)^2 above; computed here, tiny)
    fn, G = 200.0, 0.05
    f_bp, G_bp = np.array([0.2, 2.0e5]), np.array([G, G])
    for zeta in (0.01, 0.05):
        exact_sq = 0.5 * np.pi * fn / (2 * zeta) * G * (1 + 4 * zeta ** 2)
        tails = G * 0.2 + G * (2 * zeta * fn) ** 2 / 2.0e5
        full = impl.response_grms(f_bp, G_bp, fn, zeta)
        assert full ** 2 == pytest.approx(exact_sq - tails, rel=2e-4)
        assert impl.miles_grms(fn, zeta, G) ** 2 == pytest.approx(exact_sq / (1 + 4 * zeta ** 2), rel=1e-12)


def test_miles_close_to_full_integral_when_spectrum_is_flat_near_fn():
    # fn = 150 Hz sits on the flat plateau: Miles within 1 % for light damping; and the ratio -> 1 as zeta -> 0
    ratios = [impl.miles_grms(150.0, z, 0.04) / impl.response_grms(F_BP, G_BP, 150.0, z) for z in (0.005, 0.02, 0.05)]
    assert all(abs(r - 1) < 0.01 for r in ratios)
    assert abs(ratios[0] - 1) < abs(ratios[2] - 1)


def test_response_grms_converges_with_grid():
    a = impl.response_grms(F_BP, G_BP, 150.0, 0.02, n_per_decade=2000)
    b = impl.response_grms(F_BP, G_BP, 150.0, 0.02, n_per_decade=8000)
    assert a == pytest.approx(b, rel=5e-6)          # trapezoid error falls as 1/n^2


def test_relative_displacement_matches_its_closed_form_for_white_noise():
    # flat G over a wide band: integral of 1/((wn^2 - w^2)^2 + (2 zeta wn w)^2) df = 1 / (64 pi^3 zeta fn^3)
    fn, G, zeta = 120.0, 0.02, 0.04
    z = impl.relative_displacement_rms(np.array([0.05, 1.0e5]), np.array([G, G]), fn, zeta)
    exact = G0 * np.sqrt(G / (64 * np.pi ** 3 * zeta * fn ** 3))
    assert z == pytest.approx(exact, rel=1e-3)
    # displacement falls as fn^-1.5 for a flat input: stiffer mounts move less
    z2 = impl.relative_displacement_rms(np.array([0.05, 1.0e5]), np.array([G, G]), 2 * fn, zeta)
    assert z / z2 == pytest.approx(2 ** 1.5, rel=2e-3)


def test_upcrossing_rate_near_fn_for_a_narrow_band_response():
    nu = impl.upcrossing_rate(F_BP, G_BP, 150.0, 0.02)
    assert 0.95 * 150 < nu < 1.1 * 150


def test_gaussian_exceedance_and_peak_factor():
    assert impl.gaussian_exceedance(1.0) == pytest.approx(0.3173, abs=1e-4)
    assert impl.gaussian_exceedance(3.0) == pytest.approx(0.0027, abs=1e-5)
    p1, p2 = impl.expected_peak_factor(150.0, 10.0), impl.expected_peak_factor(150.0, 100.0)
    assert 3.0 < p1 < p2 < 5.0                       # a long test sees peaks beyond 3 sigma


# ---------------------------------------------------------------- time domain (SYNTHETIC)

def test_synthesized_record_has_the_target_statistics():
    fs, dur = 8192.0, 8.0
    t, a = impl.synthesize_random(F_BP, G_BP, fs, dur, seed=1)
    assert t.size == a.size == int(fs * dur) and t[1] - t[0] == pytest.approx(1 / fs)
    var = impl.psd_area(F_BP, G_BP)
    # sample variance: a sum of ~ 15 000 independent chi^2_2 bins, relative sd about 1 % (seed fixed);
    # 4 % is four standard deviations
    assert np.var(a) == pytest.approx(var, rel=0.04)
    assert abs(np.mean(a)) < 1e-9 * np.sqrt(var)                       # no DC by construction
    # Gaussian tails: fraction beyond 2 sigma is 4.55 %; with ~ 65 000 samples but correlated, allow 0.6 %
    assert np.mean(np.abs(a) > 2 * np.std(a)) == pytest.approx(0.0455, abs=0.006)
    # the same seed gives the same record; a different seed a different one
    npt.assert_array_equal(a, impl.synthesize_random(F_BP, G_BP, fs, dur, seed=1)[1])
    assert not np.allclose(a, impl.synthesize_random(F_BP, G_BP, fs, dur, seed=2)[1])


def test_synthesized_record_reproduces_the_spectrum_band_by_band():
    from scipy.signal import welch
    fs = 8192.0
    _, a = impl.synthesize_random(F_BP, G_BP, fs, 16.0, seed=4)
    f, P = welch(a, fs=fs, nperseg=4096)
    # band averages over ~ 50-1000 Welch bins of ~ 64 segments each: scatter well below 5 %
    for lo, hi in ((100, 400), (600, 1800), (40, 70)):
        sel = (f >= lo) & (f < hi)
        assert np.mean(P[sel]) == pytest.approx(np.mean(impl.psd_loglog(f[sel], F_BP, G_BP)), rel=0.05)


def test_discretization_reproduces_sinusoidal_steady_state_and_statics():
    fs, fn, zeta = 8192.0, 150.0, 0.05
    t = np.arange(int(1.5 * fs)) / fs
    for f in (40.0, 150.0, 400.0):
        y = impl.simulate_sdof_base(np.sin(2 * np.pi * f * t), 1 / fs, fn, zeta)
        tail = y[t > 1.0]                       # transient e^{-zeta wn t} ~ 1e-20 after 1 s
        amp = 0.5 * (tail.max() - tail.min())
        # peak picking on samples under-reads by at most 1 - cos(pi f / fs) (< 0.8 % at 400 Hz)
        assert amp == pytest.approx(np.sqrt(impl.transmissibility_sq(np.array([f]), fn, zeta)[0]), rel=0.008)
    # a steady base acceleration is transmitted unchanged; no input, no output
    y = impl.simulate_sdof_base(np.full(t.size, 2.5), 1 / fs, fn, zeta)
    assert y[-1] == pytest.approx(2.5, rel=1e-9)
    assert np.all(impl.simulate_sdof_base(np.zeros(100), 1 / fs, fn, zeta) == 0.0)


def test_discretization_agrees_with_scipy_lsim_after_startup():
    from scipy.signal import lsim
    fs, fn, zeta = 8192.0, 150.0, 0.05
    wn = 2 * np.pi * fn
    t, a = impl.synthesize_random(F_BP, G_BP, fs, 1.0, seed=7)
    _, yl, _ = lsim(([[0, 1], [-wn ** 2, -2 * zeta * wn]], [[0], [-1]], [[-wn ** 2, -2 * zeta * wn]], [[0]]), a, t)
    y = impl.simulate_sdof_base(a, 1 / fs, fn, zeta)
    # the two start differently (input before t = 0); that transient decays as e^{-zeta wn t}
    assert np.max(np.abs(y[t > 0.5] - yl[t > 0.5])) < 1e-9 * np.std(y)


def test_analytic_grms_matches_time_domain_simulation():
    """The lesson's independent check: integral of |H|^2 G against a seeded SYNTHETIC simulation."""
    fs, dur, fn, zeta = 8192.0, 24.0, 150.0, 0.05
    _, a = impl.synthesize_random(F_BP, G_BP, fs, dur, seed=16)
    y = impl.simulate_sdof_base(a, 1 / fs, fn, zeta)[int(0.5 * fs):]
    rms = np.sqrt(np.mean(y ** 2))
    full = impl.response_grms(F_BP, G_BP, fn, zeta)
    # relative sd of the RMS estimate is half that of the variance, 1/(2 sqrt(zeta wn T)) ~ 1.5 %;
    # allow 4 standard deviations
    u = impl.variance_estimate_relative_sd(fn, zeta, dur - 0.5) / 2
    assert 0.01 < u < 0.03
    assert rms == pytest.approx(full, rel=4 * u)
    # the 3-sigma exceedance fraction of the simulated response is close to the Gaussian 0.27 %
    assert np.mean(np.abs(y) > 3 * rms) == pytest.approx(impl.gaussian_exceedance(3.0), abs=0.001)


def test_main_writes_figure(tmp_path, capsys):
    impl.main(tmp_path)
    assert (tmp_path / "random_vibration.png").stat().st_size > 10_000
    assert "SYNTHETIC" in capsys.readouterr().out
