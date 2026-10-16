"""Public tests for P02.02 Second-Order Linear ODEs. They assert properties, with scipy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks

from conftest import load_impl

impl = load_impl(__file__)

P = impl.landing_gear()
M, K, V0, G = P["m"], P["k"], P["sink_rate"], P["g"]
XS = M * G / K
T = np.linspace(0.0, 3.0, 1501)


def ode_residual(t, m, c, k, x0, v0, dt=1e-5):
    """m x'' + c x' + k x from the closed form, with x'' by central differences of x'."""
    x, v = impl.free_response(t, m, c, k, x0, v0)
    _, vp = impl.free_response(t + dt, m, c, k, x0, v0)
    _, vm = impl.free_response(t - dt, m, c, k, x0, v0)
    return m * (vp - vm) / (2 * dt) + c * v + k * x


# ---------------------------------------------------------------- parameters and roots

def test_parameters_have_the_right_values_scalings_and_units():
    wn = impl.natural_frequency(M, K)
    assert wn == pytest.approx(np.sqrt(K / M))
    assert impl.natural_frequency(4 * M, K) == pytest.approx(wn / 2)                 # omega_n ~ 1/sqrt(m)
    assert impl.natural_frequency(M, 4 * K) == pytest.approx(2 * wn)                 # omega_n ~ sqrt(k)
    c = impl.damper_for_ratio(M, K, 0.7)
    assert impl.damping_ratio(M, c, K) == pytest.approx(0.7)
    assert impl.damped_frequency(M, c, K) == pytest.approx(wn * np.sqrt(1 - 0.49))
    assert impl.damped_frequency(M, impl.damper_for_ratio(M, K, 1.5), K) == 0.0
    with pytest.raises(ValueError):
        impl.natural_frequency(-1.0, K)
    with pytest.raises(ValueError):
        impl.damping_ratio(M, -1.0, K)


def test_characteristic_roots_match_numpy_and_classify_the_regimes():
    for zeta, regime in ((0.0, "undamped"), (0.3, "underdamped"), (1.0, "critical"), (1.8, "overdamped")):
        c = impl.damper_for_ratio(M, K, zeta)
        r = impl.characteristic_roots(M, c, K)
        npt.assert_allclose(np.sort_complex(r), np.sort_complex(np.roots([M, c, K])), rtol=1e-6, atol=1e-6)   # a double root is known only to sqrt(eps)
        assert impl.classify_damping(M, c, K) == regime
        # Vieta: sum = -c/m, product = k/m, whatever the regime.
        assert r.sum() == pytest.approx(-c / M, abs=1e-9)
        assert np.prod(r) == pytest.approx(K / M, rel=1e-10)
    wn = impl.natural_frequency(M, K)
    r = impl.characteristic_roots(M, impl.damper_for_ratio(M, K, 0.3), K)
    npt.assert_allclose(r.real, -0.3 * wn)                                             # decay rate zeta omega_n
    npt.assert_allclose(np.abs(r.imag), wn * np.sqrt(1 - 0.09))                        # damped frequency


# ---------------------------------------------------------------- closed forms

@pytest.mark.parametrize("zeta", [0.0, 0.15, 0.7, 1.0, 1.6])
def test_free_response_satisfies_the_ode_and_the_initial_conditions(zeta):
    c = impl.damper_for_ratio(M, K, zeta)
    x0, v0 = -0.05, 2.5
    x, v = impl.free_response(0.0, M, c, K, x0, v0)
    assert float(x) == pytest.approx(x0) and float(v) == pytest.approx(v0)
    t = np.array([0.05, 0.2, 0.7, 1.5])
    res = ode_residual(t, M, c, K, x0, v0)
    scale = K * abs(x0) + M * abs(v0) * impl.natural_frequency(M, K)
    assert np.max(np.abs(res)) < 1e-6 * scale
    # x' from the closed form equals the derivative of x.
    dt = 1e-6
    xp, _ = impl.free_response(t + dt, M, c, K, x0, v0)
    xm, _ = impl.free_response(t - dt, M, c, K, x0, v0)
    npt.assert_allclose((xp - xm) / (2 * dt), impl.free_response(t, M, c, K, x0, v0)[1], rtol=1e-6, atol=1e-9)


def test_three_regimes_join_continuously_at_critical_damping():
    x_lo, _ = impl.free_response(T, M, impl.damper_for_ratio(M, K, 1 - 1e-6), K, -XS, V0)
    x_cr, _ = impl.free_response(T, M, impl.damper_for_ratio(M, K, 1.0), K, -XS, V0)
    x_hi, _ = impl.free_response(T, M, impl.damper_for_ratio(M, K, 1 + 1e-6), K, -XS, V0)
    npt.assert_allclose(x_lo, x_cr, atol=1e-6)
    npt.assert_allclose(x_hi, x_cr, atol=1e-6)


def test_superposition_and_homogeneity():
    c = impl.damper_for_ratio(M, K, 0.25)
    xa, va = impl.free_response(T, M, c, K, 0.1, 0.0)
    xb, vb = impl.free_response(T, M, c, K, 0.0, 2.0)
    xab, vab = impl.free_response(T, M, c, K, 0.1, 2.0)
    npt.assert_allclose(xa + xb, xab, atol=1e-12)
    npt.assert_allclose(va + vb, vab, atol=1e-12)
    x3, _ = impl.free_response(T, M, c, K, 0.3, 6.0)
    npt.assert_allclose(x3, 3 * xab, atol=1e-12)


def test_undamped_motion_conserves_energy_and_damped_motion_loses_it_monotonically():
    x, v = impl.free_response(T, M, 0.0, K, -XS, V0)
    E = impl.mechanical_energy(x, v, M, K)
    npt.assert_allclose(E, E[0], rtol=1e-10)
    assert E[0] == pytest.approx(0.5 * M * V0 ** 2 + 0.5 * K * XS ** 2)
    c = impl.damper_for_ratio(M, K, 0.3)
    x, v = impl.free_response(T, M, c, K, -XS, V0)
    E = impl.mechanical_energy(x, v, M, K)
    assert np.all(np.diff(E) <= 1e-9 * E[0])
    # dE/dt = -c v^2 exactly: check by differences.
    dEdt = np.gradient(E, T)
    i = slice(20, -20)
    npt.assert_allclose(dEdt[i], -c * v[i] ** 2, atol=0.02 * np.max(c * v ** 2))


def test_wronskian_follows_abels_identity():
    for zeta in (0.0, 0.4, 1.0, 2.0):
        c = impl.damper_for_ratio(M, K, zeta)
        W = impl.wronskian(M, c, K, T)
        npt.assert_allclose(W, np.exp(-c * T / M), rtol=1e-8, atol=1e-12)
        assert np.all(W[np.exp(-c * T / M) > 1e-10] > 0.0)                             # never zero: independent solutions (above cancellation level)


# ---------------------------------------------------------------- numerical integration

def test_rk4_converges_at_fourth_order_to_the_closed_form():
    c = impl.damper_for_ratio(M, K, 0.2)
    errs = []
    for n in (150, 300, 600, 1200):
        ts, ys = impl.rk4_free_response(M, c, K, -XS, V0, 3.0, n)
        xe, ve = impl.free_response(ts, M, c, K, -XS, V0)
        errs.append(np.max(np.abs(ys[:, 0] - xe)))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(np.abs(rates - 4.0) < 0.3), rates
    assert errs[-1] < 1e-8
    with pytest.raises(ValueError):
        impl.rk4_free_response(M, c, K, 0.0, 1.0, 1.0, 0)


def test_first_order_form_matches_solve_ivp_and_its_eigenvalues_are_the_roots():
    c = impl.damper_for_ratio(M, K, 0.5)
    A = impl.to_first_order(M, c, K)
    assert A.shape == (2, 2)
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(A)), np.sort_complex(impl.characteristic_roots(M, c, K)), rtol=1e-10)
    sol = solve_ivp(lambda t, y: A @ y, (0.0, 3.0), [-XS, V0], rtol=1e-10, atol=1e-12, t_eval=T)
    x, v = impl.free_response(T, M, c, K, -XS, V0)
    npt.assert_allclose(sol.y[0], x, atol=1e-8)
    npt.assert_allclose(sol.y[1], v, atol=1e-7)


# ---------------------------------------------------------------- reading a record

def test_peak_finder_agrees_with_scipy_and_the_decrement_recovers_zeta():
    for zeta in (0.05, 0.2, 0.4):
        c = impl.damper_for_ratio(M, K, zeta)
        x, _ = impl.free_response(T, M, c, K, -XS, V0)
        tp, xp = impl.find_peaks_simple(T, x)
        ip, _ = find_peaks(x)
        npt.assert_allclose(tp, T[ip])
        npt.assert_allclose(xp, x[ip])
        delta = impl.logarithmic_decrement(xp)
        assert delta == pytest.approx(2 * np.pi * zeta / np.sqrt(1 - zeta ** 2), rel=2e-2)    # grid resolution of the peaks
        assert impl.zeta_from_decrement(delta) == pytest.approx(zeta, rel=2e-2)
        # Peak spacing is the damped period.
        assert np.mean(np.diff(tp)) == pytest.approx(2 * np.pi / impl.damped_frequency(M, c, K), rel=1e-2)
    assert impl.zeta_from_decrement(2 * np.pi * 0.3 / np.sqrt(1 - 0.09)) == pytest.approx(0.3)   # exact inversion
    with pytest.raises(ValueError):
        impl.logarithmic_decrement(np.array([1.0]))
    with pytest.raises(ValueError):
        impl.logarithmic_decrement(np.array([1.0, -0.5]))


def test_peak_finder_with_separation_rejects_noise_ripples():
    t = np.linspace(0.0, 1.0, 4001)
    clean = np.cos(2 * np.pi * 10 * t)
    rng = np.random.default_rng(7)
    noisy = clean + 0.05 * rng.normal(size=t.size)
    tp_all, _ = impl.find_peaks_simple(t, noisy)
    tp_sep, xp_sep = impl.find_peaks_simple(t, noisy, min_separation=0.07)
    assert tp_all.size > 3 * tp_sep.size                                                # noise maxima were everywhere
    assert 9 <= tp_sep.size <= 11                                                       # one per period, 10 periods
    assert np.min(np.diff(tp_sep)) >= 0.07
    npt.assert_allclose(tp_sep, np.round(tp_sep * 10) / 10, atol=0.01)                 # on the true peaks t = 0.1 k


def test_synthetic_record_estimates_recover_the_true_damping():
    """Seeded synthetic data: at 20 sigma threshold the estimate is within 5 % of zeta = 0.015 for noise up to 2 % of a0.
    The tolerance follows from the peak-height bias (noise adds about +0.4 sigma to each peak) being under 2 % of the
    smallest accepted peak, which is 20 sigma."""
    t, a = impl.synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, 0.5, seed=202)
    assert t.size == a.size == 4000
    # Same seed, same record; different seed, different noise.
    _, a2 = impl.synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, 0.5, seed=202)
    _, a3 = impl.synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, 0.5, seed=203)
    npt.assert_allclose(a, a2)
    assert not np.allclose(a, a3)
    est = impl.estimate_damping_from_record(t, a, 42.0, 20 * 0.5)
    assert est["f_d"] == pytest.approx(42.0 * np.sqrt(1 - 0.015 ** 2), rel=0.02)
    assert est["zeta_fit"] == pytest.approx(0.015, rel=0.05)
    assert est["zeta_decrement"] == pytest.approx(0.015, rel=0.05)
    # Noise-free record: both estimators are exact to the grid.
    t0, a0 = impl.synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, 0.0, seed=1)
    est0 = impl.estimate_damping_from_record(t0, a0, 42.0, 0.5)
    assert est0["zeta_fit"] == pytest.approx(0.015, rel=0.01)
    with pytest.raises(ValueError):
        impl.estimate_damping_from_record(t, a, 42.0, 1e3)


def test_settling_time_definition():
    t = np.linspace(0.0, 10.0, 1001)
    x = np.exp(-t)
    assert impl.settling_time(t, x, 0.02) == pytest.approx(-np.log(0.02), abs=0.02)
    assert impl.settling_time(t, 0.01 * np.ones_like(t), 0.02) == 0.0
    assert impl.settling_time(t, np.ones_like(t), 0.02) == float("inf")


# ---------------------------------------------------------------- the landing gear

def test_touchdown_response_starts_unloaded_and_relaxes_to_the_static_deflection():
    for zeta in (0.3, 0.7, 1.2):
        x, v = impl.touchdown_response(T, M, K, zeta, V0, G)
        assert float(x[0]) == pytest.approx(0.0, abs=1e-12) and float(v[0]) == pytest.approx(V0)
        assert x[-1] == pytest.approx(XS, abs=1e-4)                                    # 3 s is many time constants
        assert XS == pytest.approx(M * G / K)


def test_max_stroke_is_the_true_maximum_and_scales_with_the_energy():
    for zeta in (0.0, 0.3, 0.7, 1.0):
        t_s, x_s = impl.max_stroke(M, K, zeta, V0, G)
        x, _ = impl.touchdown_response(T, M, K, zeta, V0, G)
        assert x_s >= np.max(x) - 1e-9
        assert x_s == pytest.approx(np.max(x), rel=1e-3)
        assert 0.0 < t_s < np.pi / impl.natural_frequency(M, K)                        # before the first half period
    # Undamped: energy gives the stroke exactly, sqrt(x_s^2 + m v0^2 / k) + x_s.
    _, x_un = impl.max_stroke(M, K, 0.0, V0, G)
    assert x_un == pytest.approx(np.sqrt(XS ** 2 + M * V0 ** 2 / K) + XS, rel=1e-8)
    # More damping, less stroke; zero sink rate, no overshoot of the static value beyond the free drop.
    strokes = [impl.max_stroke(M, K, z, V0, G)[1] for z in (0.1, 0.3, 0.6, 1.0)]
    assert np.all(np.diff(strokes) < 0.0)
    t_inf, x_inf = impl.max_stroke(M, K, 2.0, V0, G)
    assert t_inf == float("inf") and x_inf == pytest.approx(XS)


def test_bounce_detection_and_the_strut_force():
    x_light, _ = impl.touchdown_response(T, M, K, 0.1, V0, G)
    x_heavy, _ = impl.touchdown_response(T, M, K, 0.6, V0, G)
    assert impl.loses_contact(x_light) and not impl.loses_contact(x_heavy)
    for zeta in (0.3, 0.7):
        F = impl.strut_force(T, M, K, zeta, V0, G)
        c = impl.damper_for_ratio(M, K, zeta)
        assert F[0] == pytest.approx(c * V0)                                            # damper alone at touchdown
        assert F[-1] == pytest.approx(M * G, rel=1e-3)                                  # supports the weight at rest
        # Newton: m x'' = -(F - m g) along the response, i.e. the force is consistent with the motion.
        x, v = impl.touchdown_response(T, M, K, zeta, V0, G)
        acc = np.gradient(v, T)
        i = slice(5, -5)
        npt.assert_allclose(M * acc[i], M * G - F[i], atol=0.02 * np.max(np.abs(F)))


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "landing_gear_response.png").stat().st_size > 1000
