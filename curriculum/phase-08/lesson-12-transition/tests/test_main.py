"""Public tests for P08.12 Transition. They assert physics and numerics, not magic numbers: the
Chebyshev matrix differentiates polynomials exactly; the Orr-Sommerfeld solver converges spectrally
to Orszag's published eigenvalue (alpha = 1, Re = 10 000) and finds the published critical point
Re_c = 5772.22, alpha_c = 1.02056 (labelled published); Poiseuille flow is stable below Re_c and
unstable above it; viscous decay rates scale as 1/Re at low Re; and the e^N envelope built on the
supplied approximate-fit chart behaves as the theory requires (linear in Re_theta beyond the
critical point, later transition for a higher N_crit or a favourable pressure gradient, and
dimensional results unchanged by a change of units).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

RE_CRIT, ALPHA_CRIT = 5772.22, 1.02056            # published (Orszag 1971)
C_ORSZAG = 0.23752649 + 0.00373967j               # published eigenvalue at alpha = 1, Re = 10 000


# ---------------------------------------------------------------- Chebyshev differentiation

def test_cheb_points_and_exact_derivatives_of_polynomials():
    N = 12
    D, x = impl.cheb(N)
    assert D.shape == (N + 1, N + 1) and x[0] == pytest.approx(1.0) and x[-1] == pytest.approx(-1.0)
    npt.assert_allclose(D @ np.ones(N + 1), 0.0, atol=1e-12)
    for k in (1, 3, 7, N):
        npt.assert_allclose(D @ x ** k, k * x ** (k - 1), atol=1e-9)


def test_cheb_converges_spectrally_for_a_smooth_function():
    errs = []
    for N in (8, 12, 16):
        D, x = impl.cheb(N)
        errs.append(np.max(np.abs(D @ np.exp(np.sin(x)) - np.cos(x) * np.exp(np.sin(x)))))
    assert errs[1] < errs[0] / 50 and errs[2] < errs[1] / 50          # far faster than any fixed power


# ---------------------------------------------------------------- Orr-Sommerfeld

def test_os_matrices_shape_and_boundary_conditions_built_in():
    A, B = impl.os_matrices(20, 1.0, 1000.0)
    assert A.shape == B.shape == (19, 19) and np.iscomplexobj(A)
    # B = D2 - alpha^2 on interior points is real; A carries the i/(alpha Re) viscous term
    assert np.allclose(B.imag, 0.0)


def test_orszag_eigenvalue_published():
    c = impl.least_stable(1.0, 1e4, 60)
    assert abs(c - C_ORSZAG) < 1e-7
    assert c.imag > 0                                                # Re = 10 000 is above critical: it grows


def test_spectral_convergence_in_N():
    errs = [abs(impl.least_stable(1.0, 1e4, N) - C_ORSZAG) for N in (20, 30, 40)]
    ratios = impl.observed_decay(errs)
    assert np.all(ratios > 20), ratios                               # each 10 extra points: > 20x smaller


def test_eigenvalues_sorted_least_stable_first_and_phase_speeds_physical():
    c = impl.os_eigenvalues(40, 1.0, 2000.0)
    assert np.all(np.diff(c.imag) <= 1e-12)
    top = c[:10]
    assert np.all((top.real > 0) & (top.real < 1.0))                # waves travel slower than the centreline


def test_stable_below_critical_unstable_above():
    for a in (0.8, 1.0, 1.2):
        assert impl.least_stable(a, 5000.0, 50).imag < 0
    assert impl.least_stable(ALPHA_CRIT, 6500.0, 50).imag > 0
    assert impl.least_stable(ALPHA_CRIT, 5600.0, 50).imag < 0


def test_viscous_decay_scales_as_one_over_re_at_low_re():
    # at Re -> 0 the flow drops out: alpha c Re tends to an eigenvalue of the Stokes operator
    g1 = impl.least_stable(1.0, 1.0, 40).imag * 1.0
    g2 = impl.least_stable(1.0, 0.5, 40).imag * 0.5
    assert g1 < 0 and g1 == pytest.approx(g2, rel=1e-3)


def test_neutral_point_has_zero_growth():
    Re_n = impl.neutral_reynolds(1.0, 50)
    assert abs(impl.least_stable(1.0, Re_n, 50).imag) < 1e-9
    with pytest.raises(ValueError):
        impl.neutral_reynolds(1.0, 50, 1000.0, 3000.0)               # no sign change: no root in the bracket


def test_critical_point_published():
    Re_c, a_c = impl.critical_point(50)
    assert Re_c == pytest.approx(RE_CRIT, rel=1e-5)
    assert a_c == pytest.approx(ALPHA_CRIT, abs=1e-3)
    # the critical point is the minimum of the neutral curve
    assert impl.neutral_reynolds(a_c - 0.05, 50) > Re_c and impl.neutral_reynolds(a_c + 0.05, 50) > Re_c


# ---------------------------------------------------------------- e^N envelope (approximate-fit chart)

def test_envelope_is_zero_before_critical_and_linear_after():
    H = 2.591
    R0 = float(impl.drela_giles_retheta0(H))
    Rt = np.array([0.5 * R0, R0, R0 + 100.0, R0 + 300.0])
    N = impl.en_envelope(Rt, H)
    assert N[0] == 0.0 and N[1] == 0.0
    assert N[3] - N[2] == pytest.approx(200.0 * float(impl.drela_giles_dn_dretheta(H)))


def test_flat_plate_transition_is_consistent_and_unit_invariant():
    nu, U = 1.863e-5, 60.0
    r = impl.flat_plate_transition(U, nu, 9.0)
    k = 0.664115
    assert r["Re_x_tr"] == pytest.approx((r["Re_theta_tr"] / k) ** 2)
    assert float(impl.en_envelope(r["Re_theta_tr"], 2.591)) == pytest.approx(9.0)
    assert r["x_tr"] == pytest.approx(r["Re_x_tr"] * nu / U) and r["x0"] < r["x_tr"]
    assert 1e6 < r["Re_x_tr"] < 5e6                                  # the familiar flat-plate range
    r_mm = impl.flat_plate_transition(U * 1e3, nu * 1e6, 9.0)        # millimetres per second
    assert r_mm["x_tr"] == pytest.approx(r["x_tr"] * 1e3) and r_mm["Re_x_tr"] == pytest.approx(r["Re_x_tr"])


def test_higher_ncrit_and_favourable_gradient_delay_transition():
    nu, U = 1.863e-5, 60.0
    xs = [impl.flat_plate_transition(U, nu, Nc)["x_tr"] for Nc in (5.0, 7.0, 9.0)]
    assert xs[0] < xs[1] < xs[2]
    assert impl.mack_ncrit(0.0007) == pytest.approx(9.0, abs=0.02)
    rex = []
    for m in (-0.03, 0.0, 0.05):
        k, H = impl.thwaites_wedge(m)
        rex.append(impl.transition_reynolds(H, k, 9.0)["Re_x_tr"])
    assert rex[0] < rex[1] < rex[2]                                   # adverse gradient: earlier; favourable: later


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "transition.png").stat().st_size > 10_000
