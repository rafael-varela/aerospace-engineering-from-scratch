"""Public tests for P02.07 Complex Numbers & Phasors. They assert properties, with NumPy and SciPy as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(207)
CASE = impl.cabin_vibration_case()


# ---------------------------------------------------------------- complex arithmetic and Euler's formula

def test_polar_round_trip_and_angle_convention():
    for _ in range(10):
        z = complex(*RNG.normal(size=2))
        r, th = impl.to_polar(z)
        assert r == pytest.approx(abs(z))
        assert -np.pi < th <= np.pi
        assert impl.from_polar(r, th) == pytest.approx(z, abs=1e-12)
    assert impl.to_polar(-1.0)[1] == pytest.approx(np.pi)          # the branch cut is on the negative real axis
    with pytest.raises(ValueError):
        impl.from_polar(-1.0, 0.0)


def test_exponential_series_matches_numpy_and_gives_eulers_formula():
    for z in (0.0, 1.0, -2.5, 1j, 2 + 3j, -4 - 1j):
        assert impl.exp_series(z) == pytest.approx(np.exp(z), rel=1e-13, abs=1e-14)
    # Euler: e^{i theta} = cos theta + i sin theta, observed from the series, not assumed.
    for th in np.linspace(-3.0, 3.0, 13):
        e = impl.exp_series(1j * th)
        assert e.real == pytest.approx(np.cos(th), abs=1e-14)
        assert e.imag == pytest.approx(np.sin(th), abs=1e-14)
        assert abs(e) == pytest.approx(1.0, abs=1e-14)             # pure imaginary exponent: unit modulus
        assert impl.euler(th) == pytest.approx(e, abs=1e-13)
    assert impl.exp_series(1j * np.pi) == pytest.approx(-1.0, abs=1e-14)   # e^{i pi} + 1 = 0
    # Too few terms must fail visibly: truncation error of the series, not rounding.
    assert abs(impl.exp_series(3.0, n_terms=4) - np.exp(3.0)) > 1.0


def test_exponential_turns_sums_into_products():
    for _ in range(5):
        a, b = complex(*RNG.normal(size=2)), complex(*RNG.normal(size=2))
        assert impl.exp_series(a + b) == pytest.approx(impl.exp_series(a) * impl.exp_series(b), rel=1e-12)


def test_de_moivre_and_nth_roots():
    for th, n in ((0.3, 5), (-1.2, 7), (2.0, 0)):
        assert impl.de_moivre(th, n) == pytest.approx(impl.euler(n * th), abs=1e-12)
    for z, n in ((8.0, 3), (-1.0, 2), (1j, 4), (3 - 4j, 5)):
        roots = impl.nth_roots(z, n)
        assert roots.shape == (n,)
        npt.assert_allclose(roots ** n, z, atol=1e-10)                  # each root raised to n gives z
        npt.assert_allclose(np.abs(roots), abs(z) ** (1.0 / n), rtol=1e-12)
        angles = np.sort(np.angle(roots) % (2 * np.pi))
        npt.assert_allclose(np.diff(angles), 2 * np.pi / n, atol=1e-10)  # equally spaced
        assert abs(np.sum(roots)) < 1e-10 if n > 1 else True           # roots of unity sum to zero
    npt.assert_allclose(np.sort_complex(impl.nth_roots(-1.0, 2)), np.sort_complex(np.array([1j, -1j])), atol=1e-12)
    with pytest.raises(ValueError):
        impl.nth_roots(1.0, 0)
    with pytest.raises(ValueError):
        impl.de_moivre(1.0, -1)


def test_complex_multiplication_is_the_rotation_matrix():
    for _ in range(5):
        v = RNG.normal(size=2)
        th = RNG.uniform(-np.pi, np.pi)
        R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
        npt.assert_allclose(impl.rotate_2d(v, th), R @ v, atol=1e-12)
        assert np.linalg.norm(impl.rotate_2d(v, th)) == pytest.approx(np.linalg.norm(v))   # length preserved
    npt.assert_allclose(impl.rotate_2d(np.array([1.0, 0.0]), np.pi / 2), [0.0, 1.0], atol=1e-12)


# ---------------------------------------------------------------- phasors and sinusoids

def test_phasor_reproduces_the_sinusoid_and_differentiation_is_multiplication_by_i_omega():
    t = np.linspace(0.0, 0.1, 2001)
    omega, A, phi = 2 * np.pi * 40.0, 1.7, -0.6
    Z = impl.phasor(A, phi)
    npt.assert_allclose(impl.from_phasor(Z, omega, t), impl.sinusoid(t, A, omega, phi), atol=1e-12)
    # d/dt [A cos(wt + phi)] = -A w sin(wt + phi) = A w cos(wt + phi + pi/2): phasor i w Z.
    dZ = impl.phasor_derivative(Z, omega)
    assert abs(dZ) == pytest.approx(A * omega)
    assert np.angle(dZ / Z) == pytest.approx(np.pi / 2)
    x = impl.from_phasor(Z, omega, t)
    dx_num = np.gradient(x, t)
    npt.assert_allclose(impl.from_phasor(dZ, omega, t)[5:-5], dx_num[5:-5], rtol=1e-4, atol=1e-6 * A * omega)
    with pytest.raises(ValueError):
        impl.phasor(-1.0, 0.0)


def test_phasor_addition_equals_pointwise_addition_of_the_signals():
    t = np.linspace(0.0, 0.05, 1001)
    omega = 2 * np.pi * 40.0
    comps = [(float(a), float(p)) for a, p in zip(RNG.uniform(0.1, 2.0, 4), RNG.uniform(-np.pi, np.pi, 4))]
    A, phi = impl.add_sinusoids(comps)
    direct = sum(impl.sinusoid(t, a, omega, p) for a, p in comps)
    npt.assert_allclose(impl.sinusoid(t, A, omega, phi), direct, atol=1e-12)
    assert A >= 0.0 and -np.pi < phi <= np.pi


def test_two_source_amplitude_is_the_cosine_rule_with_the_right_extremes():
    a1, a2 = CASE["a1"], CASE["a2"]
    lo, hi = impl.amplitude_range(a1, a2)
    assert lo == pytest.approx(a1 - a2) and hi == pytest.approx(a1 + a2)
    assert impl.combined_amplitude(a1, a2, 0.0) == pytest.approx(hi)
    assert impl.combined_amplitude(a1, a2, np.pi) == pytest.approx(lo)
    assert impl.combined_amplitude(a1, a2, np.pi / 2) == pytest.approx(np.hypot(a1, a2))   # quadrature
    for d in np.linspace(-np.pi, np.pi, 9):
        A, _ = impl.add_sinusoids([(a1, 0.0), (a2, d)])
        assert A == pytest.approx(impl.combined_amplitude(a1, a2, d))
        assert lo - 1e-12 <= A <= hi + 1e-12
        assert impl.combined_amplitude(a1, a2, -d) == pytest.approx(A)                      # symmetric in the sign of the phase


def test_fit_sinusoid_recovers_amplitude_and_phase_from_a_noisy_record():
    t = np.linspace(0.0, 0.5, 5001)
    omega, A, phi = 2 * np.pi * 40.0, 0.8, 2.1
    clean = impl.sinusoid(t, A, omega, phi)
    A_f, phi_f = impl.fit_sinusoid(t, clean, omega)
    assert A_f == pytest.approx(A, rel=1e-12) and phi_f == pytest.approx(phi, abs=1e-12)
    # Synthetic zero-mean noise of rms 0.1 on 5001 samples over 20 cycles: the fit error is about
    # 0.1 * sqrt(2 / 5001) = 0.002 in each of a, b; 10 sigma is a safe bound at a fixed seed.
    noisy = clean + 0.1 * RNG.normal(size=t.size)
    A_n, phi_n = impl.fit_sinusoid(t, noisy, omega)
    assert A_n == pytest.approx(A, abs=0.02) and phi_n == pytest.approx(phi, abs=0.03)


# ---------------------------------------------------------------- complex gain and transmissibility

def test_complex_gain_limits_and_resonance():
    zeta = 0.08
    assert impl.complex_gain(0.0, zeta) == pytest.approx(1.0)                       # static: x = F0/k
    assert abs(impl.complex_gain(50.0, zeta)) < 1e-3                                # mass-controlled: x -> 0
    mag, lag = impl.gain_and_phase_lag(1.0, zeta)
    assert mag == pytest.approx(1.0 / (2 * zeta)) and lag == pytest.approx(np.pi / 2)   # resonance: 1/(2 zeta), 90 deg
    assert impl.gain_and_phase_lag(0.0, zeta)[1] == pytest.approx(0.0, abs=1e-12)
    assert impl.gain_and_phase_lag(50.0, zeta)[1] == pytest.approx(np.pi, abs=0.01)
    # The lag grows monotonically with frequency from 0 to pi.
    lags = np.array([impl.gain_and_phase_lag(r, zeta)[1] for r in np.linspace(0.0, 5.0, 200)])
    assert np.all(np.diff(lags) >= 0.0)
    # Peak magnification occurs just below r = 1 at r = sqrt(1 - 2 zeta^2).
    rs = np.linspace(0.9, 1.1, 20001)
    assert rs[np.argmax(np.abs(impl.complex_gain(rs, zeta)))] == pytest.approx(np.sqrt(1 - 2 * zeta ** 2), abs=2e-5)


def test_phasor_steady_state_matches_the_time_domain_integration():
    m, fn, zeta = CASE["m_tray"], CASE["f_mount"], CASE["zeta"]
    wn = 2 * np.pi * fn
    k, c = m * wn ** 2, 2 * zeta * m * wn
    F0, omega = 5.0, 2 * np.pi * CASE["f_drive"]
    t = np.linspace(0.0, 1.0, 10001)
    sol = solve_ivp(impl.forced_oscillator_rhs, (0.0, 1.0), [0.0, 0.0], t_eval=t, args=(F0, omega, m, c, k),
                    rtol=1e-10, atol=1e-13)
    x_ss = impl.steady_state_displacement(t, F0, omega, m, c, k)
    tail = t > 0.6                                     # transient decays as exp(-zeta wn t): 1e-4 by 0.6 s
    X = F0 / (k - m * omega ** 2 + 1j * c * omega)
    npt.assert_allclose(sol.y[0][tail], x_ss[tail], atol=2e-4 * abs(X))
    A_fit, phi_fit = impl.fit_sinusoid(t[tail], sol.y[0][tail], omega)
    assert A_fit == pytest.approx(abs(X), rel=1e-3)
    assert phi_fit == pytest.approx(np.angle(X), abs=1e-3)
    # The phasor solution satisfies the ODE exactly (second-order finite differences at dt = 0.1 ms).
    dt = t[1] - t[0]
    xdd = (x_ss[2:] - 2 * x_ss[1:-1] + x_ss[:-2]) / dt ** 2
    xd = (x_ss[2:] - x_ss[:-2]) / (2 * dt)
    res = m * xdd + c * xd + k * x_ss[1:-1] - F0 * np.cos(omega * t[1:-1])
    assert np.max(np.abs(res)) < 2e-3 * F0


def test_transmissibility_crosses_one_at_sqrt_two_for_every_damping():
    for zeta in (0.0, 0.05, 0.3, 1.0):
        assert impl.transmissibility(0.0, zeta) == pytest.approx(1.0)
        assert impl.transmissibility(np.sqrt(2.0), zeta) == pytest.approx(1.0, rel=1e-12)
        assert impl.transmissibility(1.2, zeta) > 1.0 and impl.transmissibility(1.6, zeta) < 1.0
    # Above sqrt 2 more damping transmits *more*; at resonance it transmits less.
    assert impl.transmissibility(3.0, 0.3) > impl.transmissibility(3.0, 0.05)
    assert impl.transmissibility(1.0, 0.3) < impl.transmissibility(1.0, 0.05)
    assert impl.transmissibility(1.0, 0.05) == pytest.approx(np.sqrt(1 + 0.01) / 0.1)
    # Undamped: T = 1 / |1 - r^2| and goes as 1/r^2 at high frequency.
    assert impl.transmissibility(4.0, 0.0) == pytest.approx(1.0 / 15.0)


def test_mount_frequency_design_hits_the_target_and_is_monotone():
    f_drive, zeta = CASE["f_drive"], CASE["zeta"]
    fs = [impl.mount_frequency_for_transmissibility(T, f_drive, zeta) for T in (0.5, 0.25, 0.1)]
    for T, f in zip((0.5, 0.25, 0.1), fs):
        assert impl.transmissibility(f_drive / f, zeta) == pytest.approx(T, rel=1e-9)
        assert f < f_drive / np.sqrt(2.0)
    assert fs[0] > fs[1] > fs[2]                       # better isolation needs a softer mount
    with pytest.raises(ValueError):
        impl.mount_frequency_for_transmissibility(1.5, f_drive, zeta)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "phasor_calculator.png").stat().st_size > 1000
