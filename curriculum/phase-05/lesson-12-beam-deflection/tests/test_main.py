"""Public tests for P05.12 Beam Deflection. They assert mechanics, not magic numbers: each closed
form satisfies EI w'''' = q and its boundary conditions; superposition holds; Maxwell's reciprocal
theorem holds; numerical integration of M/EI reproduces the closed forms; an independent
finite-difference solution of (EI w'')'' = q converges to them at second order for smooth loads,
handles variable EI and an indeterminate propped cantilever, and agrees with superposition;
zero load gives zero deflection and doubling EI halves it.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

L, EI = 2.0, 5.0e4


def d4(f, x, h=1e-2):
    """Fourth derivative by central differences (exact for polynomials up to degree 5)."""
    return (f(x - 2 * h) - 4 * f(x - h) + 6 * f(x) - 4 * f(x + h) + f(x + 2 * h)) / h ** 4


def d1(f, x, h=1e-6):
    return (f(x + h) - f(x - h)) / (2 * h)


# ---------------------------------------------------------------- closed forms

def test_closed_forms_satisfy_the_beam_equation_and_boundary_conditions():
    q = 300.0
    x = np.array([0.3, 0.9, 1.6])
    npt.assert_allclose(EI * d4(lambda z: impl.cantilever_udl(z, L, q, EI), x), q, rtol=1e-4)
    npt.assert_allclose(EI * d4(lambda z: impl.simply_supported_udl(z, L, q, EI), x), q, rtol=1e-4)
    # cantilever: w = w' = 0 at the clamp; moment EI w'' and shear EI w''' vanish at the free end
    f = lambda z: impl.cantilever_udl(z, L, q, EI)
    h = 1e-3
    assert f(0.0) == 0.0 and abs(d1(f, 0.0)) < 1e-9
    assert abs(EI * (f(L + h) - 2 * f(L) + f(L - h)) / h ** 2) < 1e-4 * q * L ** 2
    assert abs(EI * (f(L + 2 * h) - 2 * f(L + h) + 2 * f(L - h) - f(L - 2 * h)) / (2 * h ** 3)) < 1e-4 * q * L
    assert f(L) == pytest.approx(q * L ** 4 / (8 * EI))
    g = lambda z: impl.simply_supported_udl(z, L, q, EI)
    assert abs(g(0.0)) < 1e-15 and abs(g(L)) < 1e-15
    assert g(L / 2) == pytest.approx(5 * q * L ** 4 / (384 * EI))


def test_point_load_forms_have_textbook_extremes_and_continuity():
    P = 1000.0
    assert impl.cantilever_point(L, L, L, P, EI) == pytest.approx(P * L ** 3 / (3 * EI))
    assert impl.simply_supported_point(L / 2, L, L / 2, P, EI) == pytest.approx(P * L ** 3 / (48 * EI))
    a = 0.7
    eps = 1e-9
    for f in (lambda z: impl.cantilever_point(z, L, a, P, EI), lambda z: impl.simply_supported_point(z, L, a, P, EI)):
        assert f(a - eps) == pytest.approx(f(a + eps), rel=1e-7)                     # continuous deflection
        assert d1(f, a - 1e-5, 1e-6) == pytest.approx(d1(f, a + 1e-5, 1e-6), rel=1e-3)   # continuous slope
    # beyond a cantilever load the beam is straight
    xs = np.array([1.0, 1.5, 2.0])
    w = impl.cantilever_point(xs, L, a, P, EI)
    assert w[2] - w[1] == pytest.approx(w[1] - w[0])
    assert impl.cantilever_end_moment(L, L, 500.0, EI) == pytest.approx(500.0 * L ** 2 / (2 * EI))


def test_maxwell_reciprocity():
    a, b, P = 0.6, 1.7, 1.0
    assert impl.cantilever_point(b, L, a, P, EI) == pytest.approx(impl.cantilever_point(a, L, b, P, EI))
    assert impl.simply_supported_point(b, L, a, P, EI) == pytest.approx(impl.simply_supported_point(a, L, b, P, EI))


def test_superposition_solves_the_propped_cantilever():
    q = 2000.0
    R = impl.propped_cantilever_reaction(L, q)
    assert R == pytest.approx(3 * q * L / 8)
    x = np.linspace(0, L, 41)
    w = impl.propped_cantilever_udl(x, L, q, EI)
    assert abs(w[0]) < 1e-15 and abs(w[-1]) < 1e-12 * np.max(np.abs(w))             # clamped and propped
    npt.assert_allclose(w, q * x ** 2 * (3 * L ** 2 - 5 * L * x + 2 * x ** 2) / (48 * EI), rtol=1e-10, atol=1e-15)
    # FD (no superposition used) agrees, and its end moment is q L^2 / 8
    xf, wf, mf = impl.fd_beam(L, EI, lambda z: np.full_like(z, q), 200, bc="clamped-pinned")
    npt.assert_allclose(wf, impl.propped_cantilever_udl(xf, L, q, EI), atol=2e-4 * np.max(np.abs(wf)))
    assert mf[0] == pytest.approx(q * L ** 2 / 8, rel=1e-3)


# ---------------------------------------------------------------- moment and integration

def test_bending_moment_from_load_matches_statics():
    x = np.linspace(0, L, 401)
    q = 150.0
    M = impl.bending_moment(x, np.full_like(x, q), support="cantilever")
    npt.assert_allclose(M, q * (L - x) ** 2 / 2, rtol=1e-4, atol=1e-6)
    M = impl.bending_moment(x, np.full_like(x, q), support="simply_supported")
    npt.assert_allclose(M, -q * x * (L - x) / 2, rtol=1e-4, atol=1e-6)        # upward load: hogging
    M = impl.bending_moment(x, np.zeros_like(x), [(0.5, 100.0)], support="cantilever")
    npt.assert_allclose(M, 100.0 * np.clip(0.5 - x, 0, None), atol=1e-9)


def test_integration_reproduces_closed_forms():
    x = np.linspace(0, L, 2001)
    q = 300.0
    M = impl.bending_moment(x, np.full_like(x, q))
    th, w = impl.deflection_by_integration(x, M, EI)
    npt.assert_allclose(w, impl.cantilever_udl(x, L, q, EI), atol=1e-6 * np.max(w))
    assert th[0] == 0.0 and w[0] == 0.0
    M = impl.bending_moment(x, np.zeros_like(x), [(1.2, -800.0)], support="simply_supported")
    th, w = impl.deflection_by_integration(x, M, EI, support="simply_supported")
    npt.assert_allclose(w, impl.simply_supported_point(x, L, 1.2, -800.0, EI), atol=1e-5 * np.max(np.abs(w)))


def test_elliptic_closed_form_matches_integration():
    s, Ls, EIw = 5.0, 13974.5, 6.958e5
    y = np.linspace(0, s, 4001)
    M = impl.bending_moment(y, impl.elliptic_lift(y, Ls, s))
    th, w = impl.deflection_by_integration(y, M, EIw)
    th_ex, w_ex = impl.elliptic_tip(Ls, s, EIw)
    assert w[-1] == pytest.approx(w_ex, rel=1e-4) and th[-1] == pytest.approx(th_ex, rel=1e-4)
    # bracketed by the same lift concentrated at its centroid and spread uniformly
    assert impl.cantilever_point(s, s, 4 * s / (3 * np.pi), Ls, EIw) < w_ex < impl.cantilever_udl(s, s, Ls / s, EIw)


# ---------------------------------------------------------------- independent check: finite differences

@pytest.mark.parametrize("bc,exact", [
    ("clamped-free", lambda x, q: impl.cantilever_udl(x, L, q, EI)),
    ("pinned-pinned", lambda x, q: impl.simply_supported_udl(x, L, q, EI)),
    ("clamped-clamped", lambda x, q: q * x ** 2 * (L - x) ** 2 / (24 * EI)),
    ("clamped-pinned", lambda x, q: impl.propped_cantilever_udl(x, L, q, EI)),
])
def test_finite_differences_converge_at_second_order(bc, exact):
    q = 400.0
    errs = []
    for n in (20, 40, 80, 160):
        x, w, _ = impl.fd_beam(L, EI, lambda z: np.full_like(z, q), n, bc=bc)
        errs.append(np.max(np.abs(w - exact(x, q))))
    orders = impl.observed_order(errs)
    assert np.all(orders > 1.8) and np.all(orders < 2.2)


def test_fd_point_load_and_superposition():
    P, a = -900.0, 1.0
    x, w, m = impl.fd_beam(L, EI, lambda z: np.zeros_like(z), 200, bc="pinned-pinned", point_loads=[(a, P)])
    npt.assert_allclose(w, impl.simply_supported_point(x, L, a, P, EI), atol=1e-3 * np.max(np.abs(w)))
    assert m[100] == pytest.approx(-P * L / 4, rel=1e-10)                     # nodal moment exact (sagging)
    # superposition: response to q1 + q2 is the sum of responses
    q1 = lambda z: 100.0 * np.sin(np.pi * z / L)
    q2 = lambda z: 50.0 + 20.0 * z
    for bc in ("clamped-free", "clamped-pinned"):
        _, w1, _ = impl.fd_beam(L, EI, q1, 100, bc=bc)
        _, w2, _ = impl.fd_beam(L, EI, q2, 100, bc=bc)
        _, w12, _ = impl.fd_beam(L, EI, lambda z: q1(z) + q2(z), 100, bc=bc)
        npt.assert_allclose(w12, w1 + w2, rtol=1e-10, atol=1e-14)


def test_fd_point_load_near_and_at_a_free_end():
    # a load one node from the tip puts a kink in m next to the boundary row; the shear condition
    # must not straddle it (a three-point one-sided difference does, and is wrong by ~50 % here)
    P, n = -100.0, 100
    for a in (L * (n - 1) / n, L):
        x, w, m = impl.fd_beam(L, EI, lambda z: np.zeros_like(z), n, bc="clamped-free", point_loads=[(a, P)])
        assert w[-1] == pytest.approx(float(impl.cantilever_point(L, L, a, P, EI)), rel=1e-3)
        assert m[0] == pytest.approx(P * a, rel=1e-9)                        # root moment from statics


def test_fd_limits_zero_load_and_stiffness_scaling():
    x, w, m = impl.fd_beam(L, EI, lambda z: np.zeros_like(z), 50)
    assert np.allclose(w, 0.0) and np.allclose(m, 0.0)
    _, w1, m1 = impl.fd_beam(L, EI, lambda z: np.full_like(z, 100.0), 50)
    _, w2, m2 = impl.fd_beam(L, 2 * EI, lambda z: np.full_like(z, 100.0), 50)
    npt.assert_allclose(w2, w1 / 2, rtol=1e-10)                               # stiffer beam, half the deflection
    npt.assert_allclose(m2, m1, rtol=1e-10, atol=1e-9)                        # same moments: determinate


def test_fd_variable_stiffness_matches_integration():
    s, Ls, EIr = 5.0, 13974.5, 6.958e5
    EIt = impl.tapered_EI(EIr, 0.7, s)
    q = lambda y: 1.5 * Ls / s * (1 - (y / s) ** 2)                            # smooth load
    y = np.linspace(0, s, 20001)
    M = impl.bending_moment(y, q(y))
    _, w_int = impl.deflection_by_integration(y, M, EIt(y))
    errs = []
    for n in (50, 100, 200):
        x, w, _ = impl.fd_beam(s, EIt, q, n)
        errs.append(abs(w[-1] - w_int[-1]))
    assert errs[-1] < 1e-4 * w_int[-1]
    assert np.all(impl.observed_order(errs) > 1.8)
    # a tapered spar (weaker outboard) deflects more than a uniform one with the root stiffness
    _, w_u = impl.deflection_by_integration(y, M, EIr)
    assert w_int[-1] > w_u[-1]


def test_elliptic_fd_converges_more_slowly_because_of_the_tip():
    s, Ls, EIw = 5.0, 13974.5, 6.958e5
    _, w_ex = impl.elliptic_tip(Ls, s, EIw)
    errs = [abs(impl.fd_beam(s, EIw, lambda y: impl.elliptic_lift(y, Ls, s), n)[1][-1] - w_ex) for n in (100, 200, 400)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 1.2) and np.all(orders < 1.8)                     # sqrt(s - y) load: below 2
    assert errs[-1] < 1e-3 * w_ex


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "wing_deflection.png").stat().st_size > 10_000
