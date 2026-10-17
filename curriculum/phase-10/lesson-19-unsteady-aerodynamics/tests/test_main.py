"""Public tests for P10.19 Unsteady Aerodynamics. They assert properties, not magic numbers:
Theodorsen's function has its limits (C -> 1 as k -> 0, C -> 1/2 as k -> infinity, F and G in
their textbook ranges); the lift tends to the quasi-steady 2 pi alpha_3/4 as k -> 0; the C(k)
implied by Jones's Wagner approximation (through the Laplace transform) agrees with the exact C(k)
to the approximation's known accuracy and with a direct numerical transform; Duhamel's integral
converges at second order to a closed-form convolution; a time-domain simulation of harmonic pitch
reproduces Theodorsen's amplitude and phase after the transient; linearity; pure plunge is always
damped; a pivot ahead of the quarter chord can be fed energy at low k; and main() writes its figure.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- Theodorsen's function

def test_theodorsen_limits_and_known_shape():
    assert impl.theodorsen(0.0) == pytest.approx(1.0)
    C = impl.theodorsen(np.array([1e-4, 1e-3, 1e-2]))
    d = np.abs(1 - C)
    assert d[0] < d[1] < d[2] < 0.06                     # tends to 1 as k -> 0 (slowly: ~ k ln k)
    assert abs(impl.theodorsen(200.0) - 0.5) < 2e-3      # high-frequency limit 1/2
    k = np.linspace(0.01, 5, 200)
    C = impl.theodorsen(k)
    assert np.all(np.diff(C.real) < 0)                   # F falls monotonically from 1 to 1/2
    assert np.all((C.real > 0.5) & (C.real < 1.0))
    assert np.all(C.imag < 0) and C.imag.min() > -0.2    # G negative: the wake always lags
    # magnitude-phase identity C = H1/(H1 + i H0) evaluated independently at one point
    from scipy.special import hankel2
    H1, H0 = hankel2(1, 0.3), hankel2(0, 0.3)
    assert impl.theodorsen(0.3) == pytest.approx(H1 / (H1 + 1j * H0), rel=1e-12)


def test_theodorsen_lift_tends_to_quasi_steady_as_k_to_zero():
    for a in (-0.5, 0.0, 0.3):
        errs = []
        for k in (1e-2, 1e-3, 1e-4):
            for alpha_hat, h_hat in ((1.0, 0.0), (0.0, 1.0 / k), (0.5, 0.3 / k)):  # plunge scaled so hdot/V stays O(1)
                c, nc = impl.theodorsen_lift(k, alpha_hat, h_hat, a)
                qs = impl.quasi_steady_lift(k, alpha_hat, h_hat, a)
                errs.append(abs(c + nc - qs) / abs(qs))
        errs = np.array(errs).reshape(3, 3).max(axis=1)
        assert errs[0] > errs[1] > errs[2] and errs[2] < 2e-3
    # steady pitch: 2 pi alpha, no matter where the pivot is
    c, nc = impl.theodorsen_lift(0.0, 0.1, 0.0, 0.4)
    assert c + nc == pytest.approx(2 * np.pi * 0.1)


def test_lift_and_moment_are_linear_in_the_motion():
    k, a = 0.4, -0.3
    c1, n1 = impl.theodorsen_lift(k, 0.2, 0.0, a)
    c2, n2 = impl.theodorsen_lift(k, 0.0, 0.7, a)
    c3, n3 = impl.theodorsen_lift(k, 0.6, 2.1, a)
    assert c3 + n3 == pytest.approx(3 * (c1 + n1 + c2 + n2))
    m1 = impl.theodorsen_moment(k, 0.2, 0.0, a)
    m2 = impl.theodorsen_moment(k, 0.0, 0.7, a)
    assert impl.theodorsen_moment(k, 0.2, 0.7, a) == pytest.approx(m1 + m2)


def test_moment_about_quarter_chord_has_no_circulatory_part():
    # the circulatory lift acts at the quarter chord, so about a = -1/2 only apparent-mass terms remain
    k = 0.5
    m = impl.theodorsen_moment(k, 1.0, 0.0, -0.5, C=0.0)
    assert impl.theodorsen_moment(k, 1.0, 0.0, -0.5) == pytest.approx(m)
    assert impl.theodorsen_moment(0.0, 1.0, 0.0, -0.5) == pytest.approx(0.0, abs=1e-14)
    # moving the pivot transfers the circulatory lift as a moment: about mid-chord, cm_qs = cl/4
    assert impl.theodorsen_moment(0.0, 0.1, 0.0, 0.0) == pytest.approx(2 * np.pi * 0.1 / 4)


# ---------------------------------------------------------------- Wagner, Kussner, transform

def test_indicial_functions_start_and_end_values():
    assert impl.wagner_jones(0.0) == pytest.approx(0.5)
    assert impl.wagner_jones(1e4) == pytest.approx(1.0)
    assert impl.kussner_sears(0.0) == pytest.approx(0.0)
    assert impl.kussner_sears(1e4) == pytest.approx(1.0)
    s = np.linspace(0, 100, 1001)
    assert np.all(np.diff(impl.wagner_jones(s)) > 0) and np.all(np.diff(impl.kussner_sears(s)) > 0)
    assert impl.rise_distance(0.9) == pytest.approx(12.70, abs=0.05)
    assert impl.wagner_jones(impl.rise_distance(0.75)) == pytest.approx(0.75)


def test_jones_transform_approximates_exact_theodorsen():
    k = np.linspace(1e-4, 5.0, 5000)
    Cj = impl.theodorsen_from_indicial(k)
    C = impl.theodorsen(k)
    assert np.max(np.abs(Cj - C)) < 0.02                 # Jones's fit is good to about 1.5 % of C
    assert impl.theodorsen_from_indicial(0.0) == pytest.approx(1.0)
    assert impl.theodorsen_from_indicial(1e6) == pytest.approx(0.5, abs=1e-6)  # 1 - sum A_i = phi(0)


def test_transform_formula_matches_a_direct_numerical_transform():
    # C(k) = phi(0) + int_0^inf phi'(tau) e^{-ik tau} dtau (Duhamel for e^{iks}, integrated by parts)
    tau = np.linspace(0, 600, 600001)
    dphi = np.gradient(impl.wagner_jones(tau), tau, edge_order=2)
    for k in (0.05, 0.3, 1.0):
        num = impl.wagner_jones(0.0) + np.trapezoid(dphi * np.exp(-1j * k * tau), tau)
        assert abs(num - impl.theodorsen_from_indicial(k)) < 2e-5


# ---------------------------------------------------------------- time domain

def test_duhamel_of_a_step_returns_the_indicial_function():
    s = np.linspace(0, 30, 301)
    x = np.full_like(s, 0.1)
    npt.assert_allclose(impl.duhamel(s, x, impl.wagner_jones), 0.1 * impl.wagner_jones(s), rtol=1e-12)


def test_duhamel_converges_at_second_order():
    k = 0.3
    A, B = (0.165, 0.335), (0.0455, 0.300)

    def exact(s):
        y = np.sin(k * s)
        for a_, b_ in zip(A, B):
            y = y - a_ * k * (b_ * np.cos(k * s) + k * np.sin(k * s) - b_ * np.exp(-b_ * s)) / (b_ ** 2 + k ** 2)
        return y

    errs = []
    for n in (100, 200, 400, 800):
        s = np.linspace(0, 40, n + 1)
        errs.append(np.max(np.abs(impl.duhamel(s, np.sin(k * s), impl.wagner_jones) - exact(s))))
    p = impl.observed_order(errs)
    assert np.all(np.abs(p - 2.0) < 0.1)


def test_time_domain_pitch_reproduces_theodorsen_amplitude_and_phase():
    a = -0.5
    for k in (0.1, 0.3427, 1.0):
        sim = impl.simulate_harmonic_pitch(k, np.radians(2.0), a)
        cj = sum(impl.theodorsen_lift(k, 1.0, 0.0, a, C=impl.theodorsen_from_indicial(k)))
        ce = sum(impl.theodorsen_lift(k, 1.0, 0.0, a))
        td = sim["total"]
        assert abs(td - cj) / abs(cj) < 1e-3              # same indicial function: discretisation only
        assert abs(td - ce) / abs(ce) < 0.03              # exact Theodorsen: within Jones's fit error
        assert abs(np.degrees(np.angle(td) - np.angle(ce))) < 1.5


def test_time_domain_lift_is_linear_and_zero_for_no_motion():
    s = np.linspace(0, 50, 2001)
    alpha = 0.05 * np.sin(0.2 * s) * (1 - np.exp(-s))
    c1, n1 = impl.time_domain_lift(s, alpha, np.zeros_like(s), -0.5)
    c2, n2 = impl.time_domain_lift(s, 3 * alpha, np.zeros_like(s), -0.5)
    npt.assert_allclose(c2 + n2, 3 * (c1 + n1), rtol=1e-10, atol=1e-12)
    c0, n0 = impl.time_domain_lift(s, np.zeros_like(s), np.zeros_like(s), 0.2)
    assert np.max(np.abs(c0 + n0)) == 0.0


def test_constant_plunge_velocity_equals_angle_of_attack():
    # descending at hdot/V = 0.05 is an angle of attack of 0.05 rad: same circulatory lift history
    s = np.linspace(0, 40, 801)
    ca, _ = impl.time_domain_lift(s, np.full_like(s, 0.05), np.zeros_like(s), -0.5)
    ch, _ = impl.time_domain_lift(s, np.zeros_like(s), 0.05 * s, -0.5)
    npt.assert_allclose(ch[5:], ca[5:], rtol=1e-6)


# ---------------------------------------------------------------- loops, energy, gusts

def test_lift_loop_area_equals_the_harmonic_formula():
    k, a = 0.3427, -0.5
    al, cl = impl.lift_loop(k, np.radians(4), np.radians(2), a, n=4001)
    c, nc = impl.theodorsen_lift(k, -1j * np.radians(2), 0.0, a)
    assert impl.loop_area(al, cl) == pytest.approx(impl.aero_work_per_cycle(-1j * np.radians(2), c + nc), rel=1e-5)
    # steady flow: no hysteresis
    al0, cl0 = impl.lift_loop(1e-9, np.radians(4), np.radians(2), a)
    assert abs(impl.loop_area(al0, cl0)) < 1e-9


def test_pure_plunge_is_always_damped():
    for k in np.linspace(0.01, 3, 40):
        c, nc = impl.theodorsen_lift(k, 0.0, 1.0, -0.5)
        # lift is up, h is down: the work done by lift on the downward motion is -oint L dh
        assert impl.aero_work_per_cycle(1.0, -(c + nc)) < 0


def test_single_dof_pitch_flutter_needs_a_forward_pivot_and_low_k():
    assert impl.pitch_flutter_boundary(-0.5) is None       # quarter-chord pivot: damped
    assert impl.pitch_flutter_boundary(0.0) is None        # mid-chord pivot: damped
    k1 = impl.pitch_flutter_boundary(-0.8)
    k2 = impl.pitch_flutter_boundary(-1.0)
    assert k1 is not None and k2 is not None and 0 < k1 < k2 < 0.2
    assert impl.aero_work_per_cycle(1.0, impl.theodorsen_moment(0.5 * k1, 1.0, 0.0, -0.8)) > 0


def test_gust_response_tends_to_quasi_steady_and_lags():
    s = np.linspace(0, 400, 40001)
    wv = np.full_like(s, 0.1)
    cl = impl.gust_response(s, wv)
    assert cl[0] == pytest.approx(0.0, abs=1e-12)
    assert cl[-1] == pytest.approx(2 * np.pi * 0.1, rel=1e-3)
    assert np.all(cl <= 2 * np.pi * 0.1 + 1e-12)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "unsteady_aero.png").stat().st_size > 10_000
