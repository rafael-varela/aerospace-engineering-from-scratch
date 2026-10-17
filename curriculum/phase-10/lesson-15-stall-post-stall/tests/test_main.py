"""Public tests for P10.15 Stall & Post-Stall. They assert physics, not magic numbers: the Kirchhoff
model reduces to thin-airfoil lift a0 alpha for f = 1 and to a quarter of it for f = 0; the
separation-point sigmoid is 1/2 at alpha1 and symmetric in alpha; the least-squares fit recovers
known parameters from noiseless model data and fits the authentic NACA 0012 lift curve to a few
hundredths; the flat-plate limit has its maximum lift at 45 deg and cd = 2 at 90 deg; the linear
lifting line matches the elliptic closed form; the nonlinear lifting line reproduces the linear one
for a straight-line section and below stall, converges in the number of Fourier terms, and agrees
with the textbook fixed point where that converges; rectangular wings stall at the root and tapered
wings further out.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import pytest

from conftest import load_impl

impl = load_impl(__file__)

SPAN, S_REF = 10.0, 12.0
AR = SPAN ** 2 / S_REF


@pytest.fixture(scope="module")
def fit():
    d = impl.read_ladson_block("I", 2.0e6)
    return impl.fit_kirchhoff(d["alpha"], d["cl"])


# ---------------------------------------------------------------- section model

def test_kirchhoff_limits():
    a = np.radians(np.linspace(-10, 10, 21))
    np.testing.assert_allclose(impl.kirchhoff_cl(a, np.ones_like(a)), 2 * np.pi * a, rtol=1e-14)   # thin airfoil
    np.testing.assert_allclose(impl.kirchhoff_cl(a, np.zeros_like(a), 5.5), 5.5 * a / 4, rtol=1e-14)
    f = np.linspace(0, 1, 11)
    assert np.all(np.diff(impl.kirchhoff_cl(0.2, f)) > 0)                                       # more attached, more lift


def test_separation_point_sigmoid():
    a1, S = np.radians(15.0), np.radians(1.0)
    assert impl.separation_point(a1, a1, S) == pytest.approx(0.5)
    assert impl.separation_point(-0.1, a1, S) == pytest.approx(impl.separation_point(0.1, a1, S))
    assert impl.separation_point(0.0, a1, S) > 0.999 and impl.separation_point(np.radians(30), a1, S) < 1e-3
    a = np.radians(np.linspace(0, 30, 61))
    assert np.all(np.diff(impl.separation_point(a, a1, S)) <= 0)


def test_fit_recovers_known_parameters():
    a_deg = np.linspace(-4, 20, 40)
    a = np.radians(a_deg)
    true = (5.9, np.radians(16.0), np.radians(1.2))
    cl = impl.kirchhoff_cl(a, impl.separation_point(a, true[1], true[2]), true[0])
    ft = impl.fit_kirchhoff(a_deg, cl)
    assert ft["a0"] == pytest.approx(true[0], rel=1e-6)
    assert ft["alpha1"] == pytest.approx(true[1], rel=1e-6)
    assert ft["S"] == pytest.approx(true[2], rel=1e-5)
    assert ft["rms"] < 1e-8


def test_fit_to_authentic_data(fit):
    d = impl.read_ladson_block("I", 2.0e6)
    assert len(d["alpha"]) == 14
    assert fit["rms"] < 0.05
    assert 0.85 * 2 * np.pi < fit["a0"] < 1.02 * 2 * np.pi                      # a little below thin-airfoil 2 pi
    assert abs(fit["cl_max_model"] - np.max(d["cl"])) < 0.1
    a_st, cl_st = impl.section_stall_angle(fit)
    assert cl_st == pytest.approx(fit["cl_max_model"], rel=1e-6)
    assert np.degrees(a_st) < np.degrees(fit["alpha1"])                         # peak before f = 1/2


def test_flat_plate_limit():
    a = np.radians(np.linspace(0, 90, 901))
    cl, cd = impl.flat_plate(a)
    assert np.degrees(a[np.argmax(cl)]) == pytest.approx(45.0, abs=0.1)
    assert cl.max() == pytest.approx(1.0) and cd[-1] == pytest.approx(2.0)
    np.testing.assert_allclose(np.hypot(cl, cd), 2 * np.sin(a), atol=1e-12)     # resultant = normal force


def test_full_range_curve_blends(fit):
    small = np.radians(np.array([2.0, 6.0]))
    np.testing.assert_allclose(impl.full_range_cl(small, fit),
                               impl.kirchhoff_cl(small, impl.separation_point(small, fit["alpha1"], fit["S"]), fit["a0"]),
                               rtol=5e-3)
    big = np.radians(np.array([50.0, 70.0, 90.0]))
    np.testing.assert_allclose(impl.full_range_cl(big, fit), impl.flat_plate(big)[0], atol=1e-3)


# ---------------------------------------------------------------- lifting line

def test_linear_lifting_line_elliptic_closed_form():
    ell = lambda y: 4 * S_REF / (np.pi * SPAN) * np.sqrt(np.clip(1 - (2 * np.asarray(y) / SPAN) ** 2, 0, None))
    zero = lambda y: np.zeros_like(np.asarray(y, dtype=float))
    for a_deg in (2.0, 5.0):
        r = impl.lifting_line_linear(np.radians(a_deg), ell, zero, 2 * np.pi)
        assert r["CL"] == pytest.approx(2 * np.pi * np.radians(a_deg) * AR / (AR + 2), rel=1e-10)
        assert r["e"] == pytest.approx(1.0, abs=1e-10)
        assert r["CDi"] == pytest.approx(r["CL"] ** 2 / (np.pi * AR), rel=1e-10)


def test_nonlinear_reduces_to_linear_for_a_straight_section():
    ch, tw = impl.twisted_wing(3.0)
    lin = impl.lifting_line_linear(np.radians(6.0), ch, tw, 5.8)
    nl = impl.lifting_line_nonlinear(np.radians(6.0), ch, tw, lambda a: 5.8 * a, a0=2 * np.pi)
    assert nl["converged"]
    assert nl["CL"] == pytest.approx(lin["CL"], rel=1e-9)
    np.testing.assert_allclose(nl["cl"], lin["cl"], rtol=1e-8)


def test_nonlinear_matches_linear_below_stall(fit):
    ch, tw = impl.rectangular_wing()
    curve = lambda a: impl.full_range_cl(a, fit)
    nl = impl.lifting_line_nonlinear(np.radians(4.0), ch, tw, curve, a0=fit["a0"])
    lin = impl.lifting_line_linear(np.radians(4.0), ch, tw, fit["a0"])
    assert nl["converged"] and nl["CL"] == pytest.approx(lin["CL"], rel=1e-5)
    np.testing.assert_allclose(nl["cl"], nl["cl_section"], atol=1e-8)          # lifting-line equation satisfied


def test_iteration_converges_and_agrees_with_fixed_point(fit):
    ch, tw = impl.rectangular_wing()
    curve = lambda a: impl.full_range_cl(a, fit)
    # just below stall, iterating with a straight line of the wrong slope (2 pi) on purpose
    nl = impl.lifting_line_nonlinear(np.radians(15.5), ch, tw, curve, a0=2 * np.pi, relax=0.5)
    assert nl["converged"]
    h = nl["history"]
    assert h[-1] < 1e-9 and h[-1] < h[0] * 1e-6
    np.testing.assert_allclose(nl["cl"], nl["cl_section"], atol=1e-8)
    fp = impl.lifting_line_fixed_point(np.radians(15.5), ch, tw, curve, relax=0.05)
    assert fp["converged"] and fp["CL"] == pytest.approx(nl["CL"], rel=1e-8)
    lim = impl.fixed_point_relaxation_limit(ch, fit["a0"])
    assert 0.03 < lim < 0.1
    bad = impl.lifting_line_fixed_point(np.radians(15.5), ch, tw, curve, relax=2.5 * lim, max_iter=3000)
    assert not bad["converged"]                                                 # beyond the limit it fails


def test_fourier_terms_convergence(fit):
    ch, tw = impl.rectangular_wing()
    curve = lambda a: impl.full_range_cl(a, fit)
    CL = [impl.lifting_line_nonlinear(np.radians(12.0), ch, tw, curve, a0=fit["a0"], n_terms=n)["CL"] for n in (4, 8, 16, 32)]
    d = np.abs(np.diff(CL))
    assert np.all(d[1:] < d[:-1] / 4)                                           # changes shrink fast


def test_stall_pattern_rectangular_root_tapered_outboard(fit):
    a_st, _ = impl.section_stall_angle(fit)
    curve = lambda a: impl.full_range_cl(a, fit)
    rect = impl.stall_onset_alpha(*impl.rectangular_wing(), curve, fit["a0"], a_st, tol_deg=0.01)
    tap = impl.stall_onset_alpha(*impl.tapered_wing(0.45), curve, fit["a0"], a_st, tol_deg=0.01)
    wash = impl.stall_onset_alpha(*impl.twisted_wing(4.0), curve, fit["a0"], a_st, tol_deg=0.01)
    assert rect["eta"] < 0.2 and tap["eta"] > rect["eta"] + 0.2
    assert wash["eta"] <= rect["eta"] + 1e-9
    for r in (rect, tap, wash):
        assert r["CL"] < fit["cl_max_model"]                                    # wing CL_max below section cl_max
        so = impl.stall_onset(r["result"], fit["cl_max_model"])
        assert so["ratio"] > 0.95


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "stall.png").stat().st_size > 10_000
