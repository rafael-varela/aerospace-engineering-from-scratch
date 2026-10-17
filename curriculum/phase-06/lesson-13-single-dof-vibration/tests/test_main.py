"""Public tests for P06.13 Single-DOF Vibration. They assert vibration physics, not magic numbers: the
bracket stiffness and Rayleigh effective mass; the free response in all three regimes satisfies its ODE
and initial conditions and matches solve_ivp; the logarithmic decrement returns the damping that made a
clean record and a noisy SYNTHETIC one within a stated tolerance; the three FRFs have their static,
resonant and high-frequency limits; the half-power method recovers zeta; the transmissibilities cross 1
at r = sqrt(2), peak near 1/(2 zeta) and have the right limits; the closed-form base-sine response from
rest matches solve_ivp and tends to the steady amplitude; the Coulomb decay is linear, frequency-
preserving, sticks, and matches an event-driven solve_ivp; and the equivalent viscous damping of
friction falls with amplitude.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

WN = 2 * np.pi * 78.0


# ---------------------------------------------------------------- model

def test_bracket_stiffness_and_effective_mass():
    E, b, t, L = 68.9e9, 0.04, 0.006, 0.1
    I = b * t ** 3 / 12
    k = impl.cantilever_stiffness(E, I, L)
    assert k == pytest.approx(3 * E * I / L ** 3)
    assert impl.cantilever_stiffness(E, I, 2 * L) == pytest.approx(k / 8)          # L^-3
    assert impl.effective_mass(0.6, 0.0) == 0.6
    assert impl.effective_mass(0.0, 1.0) == pytest.approx(33 / 140)
    assert impl.natural_frequency(4.0, 1.0) == pytest.approx(2.0)
    assert impl.damping_ratio(2 * np.sqrt(3.0 * 5.0), 3.0, 5.0) == pytest.approx(1.0)


@pytest.mark.parametrize("zeta", [0.0, 0.05, 0.7, 1.0, 2.5])
def test_free_response_satisfies_ode_and_matches_solve_ivp(zeta):
    x0, v0 = 1e-3, -0.4
    t = np.linspace(0, 0.1, 2001)
    x = impl.free_response(t, x0, v0, WN, zeta)
    assert x[0] == pytest.approx(x0)
    h = 1e-6
    v = (impl.free_response(t + h, x0, v0, WN, zeta) - impl.free_response(t - h, x0, v0, WN, zeta)) / (2 * h)
    assert v[0] == pytest.approx(v0, rel=1e-5)
    sol = solve_ivp(lambda tt, u: [u[1], -2 * zeta * WN * u[1] - WN ** 2 * u[0]], (0, 0.1), [x0, v0],
                    t_eval=t, method="DOP853", rtol=1e-11, atol=1e-15)
    npt.assert_allclose(x, sol.y[0], atol=1e-10)


def test_log_decrement_recovers_damping():
    t = np.linspace(0, 0.5, 25001)
    for zeta in (0.005, 0.02, 0.1):
        x = impl.free_response(t, 1e-4, 0.0, WN, zeta)
        _, pk = impl.find_peaks_refined(t, x, 78.0)
        assert impl.zeta_from_log_decrement(impl.log_decrement(pk, 5)) == pytest.approx(zeta, rel=1e-4)
    # exact inversion of delta = 2 pi zeta / sqrt(1 - zeta^2)
    z = 0.3
    assert impl.zeta_from_log_decrement(2 * np.pi * z / np.sqrt(1 - z * z)) == pytest.approx(z, rel=1e-12)


def test_log_decrement_on_synthetic_noisy_record():
    t, x = impl.synthetic_tap_test(seed=4, fn=78.0, zeta=0.02)
    _, pk = impl.find_peaks_refined(t, x, 78.0)
    # 1 um noise on peaks that fall from 100 to ~23 um over 10 cycles: about 1 % per peak, 0.1 % in zeta * 2 pi * 10
    assert impl.zeta_from_log_decrement(impl.log_decrement(pk, 10)) == pytest.approx(0.02, abs=0.0015)


# ---------------------------------------------------------------- FRF

def test_frf_limits_and_relations():
    m, k, zeta = 0.6, 1.5e5, 0.02
    c = 2 * zeta * np.sqrt(k * m)
    wn = np.sqrt(k / m)
    assert impl.frf(1e-6 * wn, m, c, k) == pytest.approx(1 / k, rel=1e-9)                  # static
    Hr = impl.frf(wn, m, c, k)
    assert abs(Hr) == pytest.approx(1 / (2 * zeta * k)) and np.angle(Hr) == pytest.approx(-np.pi / 2)
    w = 1e3 * wn
    assert abs(impl.frf(w, m, c, k, "accelerance")) == pytest.approx(1 / m, rel=1e-5)       # mass line
    ws = np.linspace(0.2, 3, 50) * wn
    npt.assert_allclose(impl.frf(ws, m, c, k, "mobility"), 1j * ws * impl.frf(ws, m, c, k))
    # mobility peaks exactly at w_n with value 1/c
    assert abs(impl.frf(wn, m, c, k, "mobility")) == pytest.approx(1 / c)
    with pytest.raises(ValueError):
        impl.frf(ws, m, c, k, "compliance-ish")


@pytest.mark.parametrize("zeta", [0.01, 0.02, 0.05])
def test_half_power_method(zeta):
    m, k = 0.6, 1.5e5
    c = 2 * zeta * np.sqrt(k * m)
    wn = np.sqrt(k / m)
    w = np.linspace(0.5, 1.5, 20001) * wn
    wn_est, z_est = impl.half_power_zeta(w, impl.frf(w, m, c, k))
    assert wn_est == pytest.approx(wn, rel=2 * zeta ** 2)
    assert z_est == pytest.approx(zeta, rel=2 * zeta + 1e-3)       # the method is accurate to O(zeta^2)


# ---------------------------------------------------------------- base excitation

def test_transmissibility_properties():
    r = np.array([1e-4, np.sqrt(2), 1e3])
    for zeta in (0.01, 0.1, 0.5):
        T = impl.transmissibility(r, zeta)
        assert T[0] == pytest.approx(1.0, rel=1e-6) and T[1] == pytest.approx(1.0, rel=1e-12)
        assert T[2] < 1.1 * 2 * zeta / 1e3 + 1e-6 + 1 / 1e6
    rr = np.linspace(0.5, 1.5, 200001)
    for zeta in (0.01, 0.02, 0.05):
        assert np.max(impl.transmissibility(rr, zeta)) == pytest.approx(1 / (2 * zeta), rel=2 * zeta)
    # above sqrt(2) more damping transmits more; below it, less
    assert impl.transmissibility(3.0, 0.3) > impl.transmissibility(3.0, 0.05)
    assert impl.transmissibility(1.0, 0.3) < impl.transmissibility(1.0, 0.05)
    # relative motion: 0 at low frequency, 1 at high frequency, 1/(2 zeta) at resonance
    Rt = impl.relative_transmissibility(np.array([1e-3, 1.0, 1e3]), 0.02)
    assert Rt[0] < 1e-5 and Rt[1] == pytest.approx(25.0) and Rt[2] == pytest.approx(1.0, rel=1e-5)


@pytest.mark.parametrize("ratio", [0.3, 1.0, 2.2])
def test_base_sine_from_rest_matches_solve_ivp(ratio):
    A, zeta = 2 * 9.80665, 0.02
    w = ratio * WN
    t = np.linspace(0, 0.3, 12001)
    z, a = impl.base_sine_response(t, A, w, WN, zeta)
    zn, an = impl.simulate_base_sine(t, A, w, WN, zeta)
    scale = A / WN ** 2 / (2 * zeta)
    npt.assert_allclose(z, zn, atol=1e-8 * scale)
    npt.assert_allclose(a, an, atol=1e-7 * A / (2 * zeta))
    assert z[0] == 0.0 or abs(z[0]) < 1e-15


def test_base_sine_tends_to_steady_amplitude():
    A, zeta = 2 * 9.80665, 0.02
    for ratio in (0.5, 1.0, 1.8):
        w = ratio * WN
        t = np.linspace(0, 12 / (zeta * WN), 60001)               # 12 time constants: transient below e^-12
        z, a = impl.base_sine_response(t, A, w, WN, zeta)
        tail = t > 10 / (zeta * WN)
        Zss = A / WN ** 2 * impl.relative_transmissibility(ratio, zeta) / ratio ** 2
        assert np.max(np.abs(z[tail])) == pytest.approx(Zss, rel=1e-3)
        assert np.max(np.abs(a[tail])) == pytest.approx(A * impl.transmissibility(ratio, zeta), rel=1e-3)


def test_root_stress_from_deflection():
    E, t, L = 68.9e9, 0.006, 0.1
    z = 1e-3
    I = 0.04 * t ** 3 / 12
    P = impl.cantilever_stiffness(E, I, L) * z
    assert impl.bracket_root_stress(z, E, t, L) == pytest.approx(P * L * (t / 2) / I)


# ---------------------------------------------------------------- Coulomb damping

def test_coulomb_decay_is_linear_and_keeps_frequency():
    d, x0 = 2e-6, 1e-4
    T = 2 * np.pi / WN
    t = np.linspace(0, 10 * T, 20001)
    x = impl.coulomb_free_response(t, x0, WN, d)
    # turning points at multiples of half the undamped period, amplitude falling 2d each
    for n in range(1, 8):
        assert abs(impl.coulomb_free_response(np.array([n * T / 2]), x0, WN, d)[0]) == pytest.approx(x0 - 2 * n * d, rel=1e-9)
    xs = impl.simulate_coulomb(t, x0, WN, d)
    npt.assert_allclose(x, xs, atol=1e-11)


def test_coulomb_sticks():
    d, x0 = 3e-6, 1e-5                          # 1e-5 -> 4e-6 after one half cycle... then sticks
    t = np.linspace(0, 20 * 2 * np.pi / WN, 8001)
    x = impl.coulomb_free_response(t, x0, WN, d)
    assert np.ptp(x[-1000:]) == 0.0 and abs(x[-1]) <= d
    xs = impl.simulate_coulomb(t, x0, WN, d)
    npt.assert_allclose(x, xs, atol=1e-11)


def test_equivalent_viscous_damping_of_friction():
    d = 2e-6
    z1, z2 = impl.coulomb_equivalent_zeta(d, 1e-4), impl.coulomb_equivalent_zeta(d, 2e-5)
    assert z2 == pytest.approx(5 * z1)          # inversely proportional to amplitude
    # energy check: 4 F X per cycle equals pi c_eq w_n X^2
    m, k, X = 0.6, 1.5e5, 1e-4
    c_eq = 2 * impl.coulomb_equivalent_zeta(d, X) * np.sqrt(k * m)
    assert np.pi * c_eq * np.sqrt(k / m) * X ** 2 == pytest.approx(4 * (k * d) * X)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "sdof_vibration.png").stat().st_size > 10_000
