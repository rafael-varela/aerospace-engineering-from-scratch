"""Public tests for P06.19 Home Lab: Pendulum System Identification. They assert physics and statistics,
not magic numbers: the exact amplitude correction against its series and a solve_ivp period; the
prediction's uncertainty propagation; the cycle-averaged damping laws against direct simulation for
viscous, quadratic and Coulomb damping; energy decay; the SYNTHETIC generator's reproducibility, noise
level and labelled CSV round trip; zero crossings, FFT peak and interpolated peaks on known signals; the
logarithmic decrement on an exact exponential; the lesson's independent check, recovery of the known
SYNTHETIC truth by log decrement and by nonlinear least squares within the stated uncertainty, with the
covariance checked against a seeded bootstrap and against the scatter of repeat runs; the envelope
test that tells Coulomb, viscous and quadratic damping apart; recovery of viscous and quadratic
coefficients together; the amplitude-dependent frequency; and the committed SYNTHETIC CSV.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
from pathlib import Path

import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

L, G = 0.80, 9.81
W0 = np.sqrt(G / L)
ROOT = Path(__file__).resolve().parents[4]


# ---------------------------------------------------------------- prediction and damping laws

def test_amplitude_correction():
    assert impl.amplitude_frequency_ratio(0.0) == pytest.approx(1.0, abs=1e-15)
    A = np.radians([2.0, 5.0, 10.0])
    npt.assert_allclose(impl.amplitude_frequency_ratio(A), 1 - A ** 2 / 16, atol=2e-5)   # series to A^2
    # against the period of the nonlinear pendulum from solve_ivp: the first downward zero crossing is T/4
    th0 = np.radians(30.0)
    ev = lambda t, y: y[0]
    ev.terminal, ev.direction = True, -1
    sol = solve_ivp(lambda t, y: [y[1], -W0 ** 2 * np.sin(y[0])], (0, 10), [th0, 0.0], events=ev,
                    rtol=1e-12, atol=1e-14)
    T = 4 * sol.t_events[0][0]
    assert impl.amplitude_frequency_ratio(th0) == pytest.approx(2 * np.pi / W0 / T, rel=1e-9)


def test_prediction_and_its_uncertainty():
    p = impl.predict_frequency(L, 0.004, G, 0.0, np.radians(20.0))
    assert p["f0"] == pytest.approx(W0 / (2 * np.pi), rel=1e-14)
    assert p["u_f0"] / p["f0"] == pytest.approx(0.5 * 0.004 / L, rel=1e-12)       # f ~ L^(-1/2)
    assert p["f_amp"] < p["f0"]
    q = impl.predict_frequency(L, 0.004, G, 0.02)
    assert q["u_f0"] > p["u_f0"]
    # air drag coefficient scales as d^2 L / m
    b = impl.quadratic_drag_coefficient(1.2, 0.47, 0.05, 0.15, L)
    assert impl.quadratic_drag_coefficient(1.2, 0.47, 0.10, 0.30, 2 * L) == pytest.approx(4 * b, rel=1e-12)


def test_equivalent_damping_signatures():
    A = np.array([0.05, 0.1, 0.2])
    zv = impl.equivalent_zeta(A, W0, b1=0.01)
    zq = impl.equivalent_zeta(A, W0, b2=0.01)
    zc = impl.equivalent_zeta(A, W0, b0=0.01)
    npt.assert_allclose(zv, 0.01 / (2 * W0), rtol=1e-12)          # viscous: independent of amplitude
    npt.assert_allclose(zq / A, zq[0] / A[0], rtol=1e-12)          # quadratic: grows with amplitude
    npt.assert_allclose(zc * A, zc[0] * A[0], rtol=1e-12)          # Coulomb: grows as amplitude falls


@pytest.mark.parametrize("kind", ["viscous", "quadratic", "coulomb"])
def test_cycle_averaged_decay_matches_simulation(kind):
    c = {"viscous": dict(b1=0.02), "quadratic": dict(b2=0.03), "coulomb": dict(b0=0.02)}[kind]
    t = np.linspace(0, 40, 8001)
    th = impl.simulate_swing(t, L, G, np.radians(20.0), **c)
    tp, ap = impl.find_peaks_interp(t, th)
    w = W0 * float(np.mean(impl.amplitude_frequency_ratio(ap)))
    sol = solve_ivp(lambda tt, a: impl.envelope_rate(a, w, **c), (0, tp[-1] - tp[0]), [ap[0]],
                    t_eval=tp - tp[0], rtol=1e-10)
    # averaging over a cycle neglects O(zeta) ripple within each cycle: 3 % is generous for zeta < 0.01
    npt.assert_allclose(sol.y[0], ap, rtol=0.03)


def test_energy_decays_and_is_conserved_without_damping():
    t = np.linspace(0, 20, 16001)
    E = lambda th: 0.5 * np.gradient(th, t) ** 2 + W0 ** 2 * (1 - np.cos(th))
    e0 = E(impl.simulate_swing(t, L, G, np.radians(30.0)))[5:-5]
    assert np.ptp(e0) < 3e-5 * e0.mean()                           # the O(dt^2) error of np.gradient only
    e1 = E(impl.simulate_swing(t, L, G, np.radians(30.0), b1=0.05, b2=0.05))[5:-5]
    assert e1[-1] < 0.5 * e1[0]


# ---------------------------------------------------------------- SYNTHETIC data

def test_synthetic_generator_is_reproducible_and_labelled(tmp_path):
    d1 = impl.make_synthetic_swing(3, duration=20.0)
    d2 = impl.make_synthetic_swing(3, duration=20.0)
    npt.assert_array_equal(d1["theta"], d2["theta"])
    assert not np.allclose(d1["theta"], impl.make_synthetic_swing(4, duration=20.0)["theta"])
    assert np.all(np.diff(d1["t"]) > 0) and np.mean(np.diff(d1["t"])) == pytest.approx(1 / 30, rel=1e-3)
    clean = impl.simulate_swing(d1["t"], L, G, np.radians(20.0), b1=0.003, b2=0.004)
    noise = d1["theta"] - clean - np.radians(0.5)
    # the timestamps are rounded to 1 ms, which adds up to ~0.2 deg of apparent error at peak speed
    assert np.degrees(np.std(noise)) == pytest.approx(0.3, rel=0.1)
    p = impl.write_synthetic_csv(tmp_path / "SYNTHETIC_seed3_swing.csv", 3, duration=20.0)
    text = p.read_text()
    assert text.startswith("# SYNTHETIC") and "NOT a measurement" in text
    t, th = impl.read_swing_csv(p)
    npt.assert_allclose(t, d1["t"], atol=5e-4)
    npt.assert_allclose(np.degrees(th), np.degrees(d1["theta"]), atol=6e-4)


# ---------------------------------------------------------------- frequency and peaks on known signals

def test_crossings_fft_and_peaks_on_known_signals():
    rng = np.random.default_rng(0)
    t = np.sort(np.arange(1800) / 30 + rng.normal(0, 0.002, 1800))
    f, z, wn = 0.55, 0.003, 2 * np.pi * 0.55 / np.sqrt(1 - 0.003 ** 2)
    y = impl.damped_cosine(t, 0.3, z, wn, 0.4, 0.0)
    T, uT = impl.period_from_crossings(impl.zero_crossings(t, y))
    assert 1 / T == pytest.approx(f, rel=2e-5)                     # linear interpolation error only
    assert impl.fft_peak_frequency(t, y) == pytest.approx(f, rel=2e-4)
    tp, ap = impl.find_peaks_interp(t, y)
    # peaks are measured from the record's mean; at a true maximum the damped cosine equals its envelope
    # times wd / sqrt(wd^2 + (zeta wn)^2)
    wd, sig = wn * np.sqrt(1 - z * z), z * wn
    npt.assert_allclose(ap + np.mean(y), 0.3 * np.exp(-sig * tp) * wd / np.hypot(wd, sig), rtol=2e-4)
    assert np.all(np.diff(tp) == pytest.approx(T, rel=1e-3))


def test_log_decrement_on_an_exact_exponential():
    zeta = 0.004
    delta = 2 * np.pi * zeta / np.sqrt(1 - zeta ** 2)
    peaks = 0.3 * np.exp(-delta * np.arange(40))
    d, z = impl.log_decrement(peaks, 25)
    assert d == pytest.approx(delta, rel=1e-12) and z == pytest.approx(zeta, rel=1e-12)
    fit = impl.log_decrement_fit(peaks)
    assert fit["zeta"] == pytest.approx(zeta, rel=1e-10) and fit["u_zeta"] < 1e-12


# ---------------------------------------------------------------- independent check: recover the SYNTHETIC truth

@pytest.fixture(scope="module")
def viscous_record():
    zeta = 0.002
    d = impl.make_synthetic_swing(0, theta0_deg=4.0, duration=60.0, b1=2 * zeta * W0, b2=0.0)
    return d, zeta


def test_nls_and_log_decrement_recover_the_known_truth(viscous_record):
    d, zeta = viscous_record
    t, th = d["t"], d["theta"]
    fit = impl.fit_damped_cosine(t, th)
    A, z, wn, phi, c = fit["p"]
    eA, ez, ew, ephi, ec = fit["perr"]
    # the fitted frequency is the amplitude-averaged one: compare against the exact law at the peaks
    tp, ap = impl.find_peaks_interp(t, th)
    w_true = W0 * float(np.mean(impl.amplitude_frequency_ratio(ap)))
    assert abs(z - zeta) < 3 * ez
    assert abs(wn - w_true) < 3 * ew
    assert abs(c - np.radians(0.5)) < 3 * ec                        # the protractor offset
    assert np.degrees(fit["rms"]) == pytest.approx(0.3, rel=0.1)  # residuals = the noise: model adequate
    ld = impl.log_decrement_fit(ap)
    assert abs(ld["zeta"] - zeta) < 3 * ld["u_zeta"]
    assert abs(ld["zeta"] - z) < 3 * np.hypot(ld["u_zeta"], ez)     # the two methods agree


def test_covariance_agrees_with_bootstrap(viscous_record):
    d, _ = viscous_record
    fit = impl.fit_damped_cosine(d["t"], d["theta"])
    B = impl.bootstrap_fit(d["t"], d["theta"], fit, n_boot=60, seed=3)
    ratio = B.std(axis=0, ddof=1)[1:3] / fit["perr"][1:3]
    # 60 resamples: the sd of an sd estimate is ~ 1/sqrt(2 x 59) = 9 %; allow 0.67-1.5
    assert np.all((ratio > 0.67) & (ratio < 1.5))


def test_covariance_agrees_with_repeat_run_scatter():
    b1 = 2 * 0.002 * W0
    est = []
    for seed in range(12):
        d = impl.make_synthetic_swing(seed, theta0_deg=4.0, duration=60.0, b1=b1, b2=0.0)
        f = impl.fit_damped_cosine(d["t"], d["theta"])
        est.append((f["p"][1], f["perr"][1], f["p"][2], f["perr"][2]))
    est = np.array(est)
    r_z = est[:, 0].std(ddof=1) / est[:, 1].mean()
    r_w = est[:, 2].std(ddof=1) / est[:, 3].mean()
    # 12 runs: s/sigma lies in [0.49, 1.56] with 99 % probability (chi distribution, 11 dof)
    assert 0.45 < r_z < 1.6 and 0.45 < r_w < 1.6


# ---------------------------------------------------------------- diagnosis

@pytest.mark.parametrize("kind,p_expected", [("coulomb", 0.0), ("viscous", 1.0), ("quadratic", 2.0)])
def test_envelope_exponent_tells_the_damping_apart(kind, p_expected):
    c = {"viscous": dict(b1=0.02, b2=0.0), "quadratic": dict(b1=0.0, b2=0.03),
         "coulomb": dict(b1=0.0, b2=0.0, b0=0.02)}[kind]
    d = impl.make_synthetic_swing(5, theta0_deg=30.0, duration=60.0, noise_deg=0.1, offset_deg=0.0, **c)
    tp, ap = impl.find_peaks_interp(d["t"], d["theta"])
    e = impl.fit_envelope(tp, ap)
    assert abs(e["p"] - p_expected) < max(3 * e["u_p"], 0.15)
    assert e["u_p"] < 0.3                                          # decisive: the three are > 3 sigma apart


def test_viscous_and_quadratic_coefficients_recovered_together():
    d = impl.make_synthetic_swing(7, theta0_deg=40.0, duration=150.0, b1=0.01, b2=0.02, noise_deg=0.2,
                                  offset_deg=0.0)
    tp, ap = impl.find_peaks_interp(d["t"], d["theta"])
    vq = impl.fit_envelope_vq(tp, ap, W0)
    assert abs(vq["b1"] - 0.01) < 3 * vq["u_b1"] and abs(vq["b2"] - 0.02) < 3 * vq["u_b2"]
    assert vq["u_b1"] < 0.1 * 0.01 and vq["u_b2"] < 0.1 * 0.02       # a 40 -> 6 deg decay separates them


def test_frequency_falls_with_amplitude_as_the_exact_law():
    d = impl.make_synthetic_swing(9, theta0_deg=40.0, duration=120.0, b1=0.0, b2=0.03, noise_deg=0.1,
                                  offset_deg=0.0)
    A, f = impl.frequency_vs_amplitude(d["t"], d["theta"], cycles=6)
    assert A[0] > A[-1] and f[0] < f[-1]                            # larger swing, lower frequency
    pred = W0 / (2 * np.pi) * impl.amplitude_frequency_ratio(A)
    # window-averaged amplitude against the exact law: agreement to 0.05 %
    npt.assert_allclose(f, pred, rtol=5e-4)


# ---------------------------------------------------------------- the committed SYNTHETIC file and main

def test_committed_synthetic_csv_is_labelled_and_reproducible():
    p = ROOT / "projects" / "p06-pendulum-sysid" / "data" / "SYNTHETIC_seed619_swing.csv"
    assert p.exists() and "SYNTHETIC" in p.name
    assert p.read_text().startswith("# SYNTHETIC")
    t, th = impl.read_swing_csv(p)
    d = impl.make_synthetic_swing(619)
    assert t.size == d["t"].size
    npt.assert_allclose(np.degrees(th), np.degrees(d["theta"]), atol=6e-4)   # 3 decimals in the file


def test_main_writes_figure(tmp_path, capsys):
    impl.main(tmp_path)
    assert (tmp_path / "pendulum_sysid.png").stat().st_size > 10_000
    assert "SYNTHETIC" in capsys.readouterr().out
