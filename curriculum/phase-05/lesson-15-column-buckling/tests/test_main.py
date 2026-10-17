"""Public tests for P05.15 Column Buckling. They assert properties, not magic numbers: the Euler load
scales as EI/(KL)^2 and the four end conditions stand in the ratios 1/4 : 1 : 2.046 : 4; the
finite-difference eigenvalue of E I w'''' + P w'' = 0 converges at second order to every closed form; the
pinned mode is a half sine and satisfies the equilibrium E I w'' + P w = 0; the Johnson parabola starts at
F_cy, meets Euler at the transition slenderness with equal value and slope, and never exceeds Euler or
F_cy; the tangent-modulus stress solves its own equation, lies below Euler and tends to it for slender
columns.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

E, FCY = 205e9, 435e6
L = 0.5
EI = E * 3.385e-9


# ---------------------------------------------------------------- Euler

def test_euler_load_scaling():
    P = impl.euler_load(E, 1e-8, 1.0)
    assert P == pytest.approx(np.pi ** 2 * E * 1e-8)
    assert impl.euler_load(2 * E, 1e-8, 1.0) == pytest.approx(2 * P)
    assert impl.euler_load(E, 3e-8, 1.0) == pytest.approx(3 * P)
    assert impl.euler_load(E, 1e-8, 2.0) == pytest.approx(P / 4)
    assert impl.euler_load(E, 1e-8, 1.0, K=2.0) == pytest.approx(P / 4)


def test_effective_length_factors_and_load_ratios():
    K = {e: impl.effective_length_factor(e) for e in ("pinned-pinned", "fixed-free", "fixed-pinned", "fixed-fixed")}
    assert K["pinned-pinned"] == 1.0 and K["fixed-free"] == 2.0 and K["fixed-fixed"] == 0.5
    beta = np.pi / K["fixed-pinned"]
    assert np.tan(beta) == pytest.approx(beta, rel=1e-9)            # tan(beta L) = beta L
    ratio = impl.euler_load(E, 1e-8, 1.0, K["fixed-pinned"]) / impl.euler_load(E, 1e-8, 1.0)
    assert ratio == pytest.approx(2.046, abs=1e-3)


def test_euler_stress_is_load_over_area():
    p = impl.tube_properties(0.0222, 0.000889)
    s = 0.7 / p["r"]
    assert impl.euler_stress(E, s) * p["A"] == pytest.approx(impl.euler_load(E, p["I"], 0.7))
    assert p["r"] == pytest.approx(np.sqrt((0.0222 ** 2 + (0.0222 - 2 * 0.000889) ** 2) / 16))


# ---------------------------------------------------------------- finite differences (the independent check)

@pytest.mark.parametrize("ends", ["pinned-pinned", "fixed-free", "fixed-pinned", "fixed-fixed"])
def test_fd_eigenvalue_converges_at_second_order_to_euler(ends):
    exact = impl.euler_load(E, 3.385e-9, L, impl.effective_length_factor(ends))
    errs = [abs(impl.fd_buckling(L, EI, ends, n)[0] - exact) / exact for n in (40, 80, 160)]
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.8) and np.all(orders < 2.3)
    assert errs[-1] < 5e-4


def test_fd_pinned_mode_is_a_half_sine_in_equilibrium():
    P, x, w = impl.fd_buckling(L, EI, "pinned-pinned", 200)
    w = w * np.sign(w[len(w) // 2])
    npt.assert_allclose(w, np.sin(np.pi * x / L), atol=1e-6)
    assert w[0] == pytest.approx(0.0, abs=1e-12) and w[-1] == pytest.approx(0.0, abs=1e-12)
    # at the critical load the bent shape satisfies E I w'' + P w = 0 (moment equilibrium of a cut piece)
    assert impl.euler_mode_residual(x, w, P, EI) < 1e-3
    # at any other load it does not: equilibrium in the bent shape exists only at P_cr
    assert impl.euler_mode_residual(x, w, 1.2 * P, EI) > 0.1


def test_fd_fixed_free_mode_has_zero_slope_at_the_base_and_maximum_at_the_tip():
    P, x, w = impl.fd_buckling(L, EI, "fixed-free", 200)
    w = w * np.sign(w[-1])
    assert w[0] == pytest.approx(0.0, abs=1e-12)
    assert abs(w[1] - w[0]) < 1e-3                                    # w'(0) = 0
    assert np.argmax(w) == len(w) - 1
    npt.assert_allclose(w, 1 - np.cos(np.pi * x / (2 * L)), atol=1e-4)


def test_fd_load_is_proportional_to_EI():
    P1 = impl.fd_buckling(L, EI, "fixed-fixed", 60)[0]
    P2 = impl.fd_buckling(L, 3 * EI, "fixed-fixed", 60)[0]
    assert P2 == pytest.approx(3 * P1, rel=1e-9)


# ---------------------------------------------------------------- Johnson

def test_johnson_parabola_meets_euler_tangentially():
    sc = impl.transition_slenderness(E, FCY)
    assert impl.johnson_stress(E, FCY, sc) == pytest.approx(impl.euler_stress(E, sc), rel=1e-12)
    assert impl.euler_stress(E, sc) == pytest.approx(FCY / 2)
    h = 1e-4 * sc
    dj = (impl.johnson_stress(E, FCY, sc + h) - impl.johnson_stress(E, FCY, sc - h)) / (2 * h)
    de = (impl.euler_stress(E, sc + h) - impl.euler_stress(E, sc - h)) / (2 * h)
    assert dj == pytest.approx(de, rel=1e-6)
    assert impl.johnson_stress(E, FCY, 0.0) == pytest.approx(FCY)


def test_column_curve_is_continuous_monotone_and_bounded():
    s = np.linspace(1.0, 250.0, 2000)
    c = impl.critical_stress(E, FCY, s)
    assert np.all(np.diff(c) < 0)
    assert np.all(c <= FCY) and np.all(c <= impl.euler_stress(E, s) * (1 + 1e-12))
    assert np.max(np.abs(np.diff(c))) < 0.01 * FCY                    # no jump at the transition


def test_strut_margin():
    assert impl.strut_margin(15e3, 10e3) == pytest.approx(0.0)
    assert impl.strut_margin(30e3, -10e3) == pytest.approx(1.0)


# ---------------------------------------------------------------- tangent modulus

def test_ramberg_osgood_and_its_tangent():
    n = 15.0
    assert impl.ramberg_osgood_strain(FCY, E, FCY, n) == pytest.approx(FCY / E + 0.002)
    s = np.array([50e6, 200e6, 400e6])
    h = 1e3
    num = 2 * h / (impl.ramberg_osgood_strain(s + h, E, FCY, n) - impl.ramberg_osgood_strain(s - h, E, FCY, n))
    npt.assert_allclose(impl.tangent_modulus(s, E, FCY, n), num, rtol=1e-5)
    assert impl.tangent_modulus(1e6, E, FCY, n) == pytest.approx(E, rel=1e-9)   # elastic at low stress
    assert np.all(np.diff(impl.tangent_modulus(np.linspace(1e6, 500e6, 50), E, FCY, n)) <= 0)   # non-increasing


def test_tangent_modulus_stress_solves_its_equation_and_tends_to_euler():
    n = 15.0
    for sl in (20.0, 50.0, 80.0, 150.0):
        s = impl.tangent_modulus_stress(sl, E, FCY, n)
        assert s == pytest.approx(np.pi ** 2 * impl.tangent_modulus(s, E, FCY, n) / sl ** 2, rel=1e-9)
        assert s <= impl.euler_stress(E, sl) * (1 + 1e-12)
    assert impl.tangent_modulus_stress(150.0, E, FCY, n) == pytest.approx(impl.euler_stress(E, 150.0), rel=1e-3)
    assert impl.tangent_modulus_stress(10.0, E, FCY, n) < 1.1 * FCY


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "column_buckling.png").stat().st_size > 10_000
