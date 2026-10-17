"""Public tests for P06.17 Vibration Isolation. They assert physics, not magic numbers: stiffness and
frequency are inverse operations; the transmissibility is 1 at r = 0 and at r = sqrt 2 for every damping,
below 1 beyond, 1/|1 - r^2| without damping, with the right phase; viscous damping rolls off at
-20 dB/decade and hysteretic damping at -40; the exact resonant peak agrees with a brute-force search and
tends to 1/(2 zeta); the inverse design functions invert; the closed form agrees with a solve_ivp steady
state (the lesson's first independent check); the random response tends to the input for a rigid mount,
to Miles for light damping and to a seeded SYNTHETIC time-domain simulation (the second check); and the
design search returns a feasible optimum and refuses an infeasible specification.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G0 = 9.80665
F_BP = np.array([20.0, 80.0, 500.0, 2000.0])
G_BP = np.array([0.01, 0.04, 0.04, 0.01])


# ---------------------------------------------------------------- isolator basics

def test_stiffness_frequency_and_damping_are_consistent():
    m = 0.35
    k = impl.stiffness_for_frequency(65.0, m)
    assert impl.natural_frequency(k, m) == pytest.approx(65.0, rel=1e-14)
    c = impl.damping_coefficient(0.2, k, m)
    assert c / (2 * np.sqrt(k * m)) == pytest.approx(0.2, rel=1e-14)
    # static sag under n g is n g0 / wn^2 = n g0 m / k
    assert impl.static_deflection(5.0, 65.0) == pytest.approx(5.0 * G0 * m / k, rel=1e-12)


def test_transmissibility_limits_and_crossover():
    r = np.array([1e-6, np.sqrt(2), 2.0, 5.0])
    for zeta in (0.0, 0.05, 0.3, 1.0):
        T = impl.transmissibility(r, zeta)
        assert T[0] == pytest.approx(1.0, abs=1e-9)
        assert T[1] == pytest.approx(1.0, rel=1e-12)              # isolation starts at r = sqrt 2 for every zeta
        assert np.all(T[2:] < 1.0)
    rr = np.array([0.3, 0.8, 1.7, 4.0])
    npt.assert_allclose(impl.transmissibility(rr, 0.0), 1 / np.abs(1 - rr ** 2), rtol=1e-12)


def test_transmissibility_phase():
    zeta = 0.1
    ph = impl.transmissibility_phase(np.array([1e-4, 1.0, 1e4]), zeta)
    assert ph[0] == pytest.approx(0.0, abs=1e-3)
    assert ph[1] == pytest.approx(np.arctan(1 / (2 * zeta)), rel=1e-12)   # H(1) = 1 - i/(2 zeta)
    assert ph[2] == pytest.approx(np.pi / 2, abs=1e-3)                    # viscous: lag tends to 90 deg
    # the hysteretic model has no damper in series with the base, so it rolls off faster
    r = np.array([100.0, 1000.0])
    slope_v = np.diff(np.log10(impl.transmissibility(r, 0.1)))[0]
    slope_h = np.diff(np.log10(impl.transmissibility_hysteretic(r, 0.2)))[0]
    assert slope_v == pytest.approx(-1.0, abs=0.01)                       # -20 dB per decade
    assert slope_h == pytest.approx(-2.0, abs=0.01)                       # -40 dB per decade


def test_more_damping_lowers_the_peak_but_raises_high_frequency_transmission():
    r_hi = np.array([5.0, 10.0])
    T_lo, T_hi = impl.transmissibility(r_hi, 0.05), impl.transmissibility(r_hi, 0.3)
    assert np.all(T_hi > T_lo)
    assert impl.peak_transmissibility(0.3)[1] < impl.peak_transmissibility(0.05)[1]


def test_peak_transmissibility_is_exact():
    for zeta in (0.01, 0.1, 0.3, 0.6):
        rp, Tp = impl.peak_transmissibility(zeta)
        r = np.linspace(0.5 * rp, min(1.5 * rp, 1.4), 400001)
        T = impl.transmissibility(r, zeta)
        assert Tp >= T.max() * (1 - 1e-13)                      # nothing on a fine grid beats it
        assert Tp == pytest.approx(T.max(), rel=1e-7)            # grid step limits the search
        assert rp == pytest.approx(r[np.argmax(T)], abs=1e-5)
        assert rp < 1.0
    # light damping: peak -> Q = 1/(2 zeta) at r -> 1
    rp, Tp = impl.peak_transmissibility(0.005)
    assert Tp == pytest.approx(1 / 0.01, rel=1e-3) and rp == pytest.approx(1.0, abs=1e-4)


def test_inverse_design_functions():
    z = impl.zeta_for_peak(4.0)
    assert impl.peak_transmissibility(z)[1] == pytest.approx(4.0, rel=1e-10)
    fn = impl.fn_for_attenuation(500.0, 0.05, 0.15)
    assert impl.transmissibility(np.array([500.0 / fn]), 0.15)[0] == pytest.approx(0.05, rel=1e-10)
    # undamped: T = 1/(r^2 - 1)  =>  r = sqrt(1 + 1/T)
    assert impl.fn_for_attenuation(500.0, 0.05, 0.0) == pytest.approx(500.0 / np.sqrt(21.0), rel=1e-10)
    # damping costs attenuation, so the same target needs a softer mount
    assert impl.fn_for_attenuation(500.0, 0.05, 0.3) < fn


# ---------------------------------------------------------------- independent check 1: solve_ivp

@pytest.mark.parametrize("r,zeta", [(0.5, 0.1), (1.0, 0.1), (np.sqrt(2), 0.3), (3.0, 0.05), (8.0, 0.2)])
def test_closed_form_matches_solve_ivp_steady_state(r, zeta):
    amp, lag = impl.steady_state_ivp(r, zeta)
    # DOP853 at rtol 1e-10 with the transient decayed to 1e-8: agreement to 1e-6 is generous
    assert amp == pytest.approx(float(impl.transmissibility(np.array([r]), zeta)[0]), rel=1e-6)
    assert lag == pytest.approx(float(impl.transmissibility_phase(np.array([r]), zeta)[0]), abs=1e-6)


# ---------------------------------------------------------------- random response

def test_random_response_limits():
    g_in = impl.psd_grms(F_BP, G_BP)
    stiff, z_stiff = impl.random_response(F_BP, G_BP, 1.0e5, 0.05)        # effectively rigid
    assert stiff == pytest.approx(g_in, rel=1e-3) and z_stiff < 1e-8
    soft, _ = impl.random_response(F_BP, G_BP, 2.0, 0.05)                 # mode far below the band
    assert soft < 0.02 * g_in
    g = [impl.random_response(F_BP, G_BP, fn, 0.15)[0] for fn in (40.0, 80.0, 160.0)]
    z = [impl.random_response(F_BP, G_BP, fn, 0.15)[1] for fn in (40.0, 80.0, 160.0)]
    assert g[0] < g[1] < g[2] and z[0] > z[1] > z[2]                      # the isolation trade


def test_random_response_against_miles_and_white_noise_closed_form():
    # on the flat plateau with light damping, Miles is within 1 %
    g, _ = impl.random_response(F_BP, G_BP, 200.0, 0.02)
    assert impl.miles_grms(200.0, 0.02, 0.04) == pytest.approx(g, rel=0.01)
    # white input: relative displacement RMS = g0 sqrt(G / (64 pi^3 zeta fn^3)) (P06.16)
    fn, G, zeta = 70.0, 0.03, 0.15
    _, zr = impl.random_response(np.array([0.05, 1.0e5]), np.array([G, G]), fn, zeta)
    assert zr == pytest.approx(G0 * np.sqrt(G / (64 * np.pi ** 3 * zeta * fn ** 3)), rel=2e-3)


def test_random_response_matches_synthetic_time_domain():
    """Independent check 2: the integral against a seeded SYNTHETIC simulation (given helpers)."""
    fs, dur, fn, zeta = 8192.0, 24.0, 65.0, 0.2
    _, a = impl.synthesize_random(F_BP, G_BP, fs, dur + 0.5, seed=17)
    acc, z = impl.simulate_sdof_base(a, 1 / fs, fn, zeta)
    sk = int(0.5 * fs)
    g, zr = impl.random_response(F_BP, G_BP, fn, zeta)
    # relative sd of an RMS over T is about 1/(2 sqrt(zeta wn T)) = 1.1 % here; allow 4 sd
    u = 0.5 / np.sqrt(zeta * 2 * np.pi * fn * dur)
    assert np.sqrt(np.mean(acc[sk:] ** 2)) == pytest.approx(g, rel=4 * u)
    assert G0 * np.sqrt(np.mean(z[sk:] ** 2)) == pytest.approx(zr, rel=4 * u)


# ---------------------------------------------------------------- design

def test_design_sweep_trends():
    fns, zetas = np.array([50.0, 100.0, 150.0]), np.array([0.05, 0.2])
    sw = impl.design_sweep(fns, zetas, F_BP, G_BP, n_static=5.0, f_hf=500.0)
    assert sw["grms"].shape == sw["sway3"].shape == sw["tpeak"].shape == sw["thf"].shape == (3, 2)
    assert np.all(np.diff(sw["sway3"], axis=0) < 0)          # stiffer mounts sway less
    assert np.all(np.diff(sw["grms"], axis=0) > 0)           # ... and transmit more
    assert np.all(sw["thf"][:, 1] > sw["thf"][:, 0])         # damping costs high-frequency isolation
    assert np.all(sw["tpeak"][:, 1] < sw["tpeak"][:, 0])     # ... and buys a lower peak
    # with n_static = 0 the sway is only the dynamic part, 3 sigma
    sw0 = impl.design_sweep(fns[:1], zetas[:1], F_BP, G_BP)
    assert sw0["sway3"][0, 0] == pytest.approx(3 * impl.random_response(F_BP, G_BP, 50.0, 0.05)[1], rel=1e-12)


def test_choose_isolator_returns_the_feasible_optimum():
    fns = np.arange(40.0, 161.0, 10.0)
    zetas = np.array([0.05, 0.1, 0.15, 0.2, 0.25])
    kw = dict(fn_min=60.0, sway_max=0.8e-3, tpeak_max=4.0, thf_max=0.06, n_static=5.0, f_hf=500.0)
    best = impl.choose_isolator(fns, zetas, F_BP, G_BP, 0.35, **kw)
    assert best["fn"] >= 60.0 and best["sway3"] <= 0.8e-3 and best["tpeak"] <= 4.0 and best["thf"] <= 0.06
    assert impl.natural_frequency(best["k"], 0.35) == pytest.approx(best["fn"], rel=1e-12)
    assert best["c"] == pytest.approx(impl.damping_coefficient(best["zeta"], best["k"], 0.35))
    sw = impl.design_sweep(fns, zetas, F_BP, G_BP, 5.0, 500.0)
    F = np.repeat(fns[:, None], zetas.size, axis=1)
    ok = (F >= 60.0) & (sw["sway3"] <= 0.8e-3) & (sw["tpeak"] <= 4.0) & (sw["thf"] <= 0.06)
    assert best["grms"] == pytest.approx(sw["grms"][ok].min(), rel=1e-12)
    with pytest.raises(ValueError):
        impl.choose_isolator(fns, zetas, F_BP, G_BP, 0.35, fn_min=60.0, sway_max=0.1e-3, tpeak_max=4.0)


def test_main_writes_figure(tmp_path, capsys):
    impl.main(tmp_path)
    assert (tmp_path / "vibration_isolation.png").stat().st_size > 10_000
    assert "SYNTHETIC" in capsys.readouterr().out
