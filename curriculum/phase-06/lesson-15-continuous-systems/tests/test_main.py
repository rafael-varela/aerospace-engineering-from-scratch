"""Public tests for P06.15 Continuous Systems. They assert beam dynamics, not magic numbers: the roots of
the frequency equations, their limits and their asymptotes; the mode shapes satisfy w'''' = beta^4 w and
their boundary conditions and are orthogonal; the element matrices are symmetric, have rigid-body null
spaces and the right total mass; the finite-element eigenproblem converges to the closed form at fourth
order for three supports and with a tip mass; the finite-difference pinned beam converges at second
order and equals its own discrete closed form; the SDOF (Rayleigh) estimate is an upper bound that
becomes exact for a heavy tip mass; and the FE stiffness reproduces the static cantilever of P05.12.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

EI, RA, L = 6.958e5, 9.0, 5.0


# ---------------------------------------------------------------- closed form

def test_roots_solve_the_frequency_equations():
    c = impl.beta_L("clamped-free", 5)
    assert np.all(np.abs(1 + np.cos(c) * np.cosh(c)) < 1e-9 * np.cosh(c))
    assert c[0] == pytest.approx(1.875104, abs=1e-6)
    npt.assert_allclose(impl.beta_L("pinned-pinned", 4), np.pi * np.arange(1, 5), rtol=1e-12)
    cc = impl.beta_L("clamped-clamped", 4)
    assert np.all(np.abs(np.cos(cc) * np.cosh(cc) - 1.0) < 1e-9 * np.cosh(cc))
    # asymptotes: (2n - 1) pi / 2 for the cantilever, (2n + 1) pi / 2 for clamped-clamped
    assert c[4] == pytest.approx(9 * np.pi / 2, rel=1e-6)
    assert cc[3] == pytest.approx(9 * np.pi / 2, rel=1e-6)
    with pytest.raises(ValueError):
        impl.frequency_equation(1.0, "sliding")


def test_frequency_scaling():
    w = impl.natural_frequencies(EI, RA, L, "clamped-free", 2)
    assert w[0] == pytest.approx(1.875104 ** 2 * np.sqrt(EI / (RA * L ** 4)), rel=1e-6)
    w2 = impl.natural_frequencies(4 * EI, RA, 2 * L, "clamped-free", 2)
    npt.assert_allclose(w2, w / 2, rtol=1e-10)                 # sqrt(EI) / L^2


def test_tip_mass_limits_and_sdof_estimate():
    w0 = impl.natural_frequencies(EI, RA, L, "clamped-free", 1)[0]
    assert impl.natural_frequencies(EI, RA, L, "clamped-free", 1, m_tip=0.0)[0] == pytest.approx(w0)
    prev = w0
    for mt in (5.0, 40.0, 400.0):
        wt = impl.natural_frequencies(EI, RA, L, "clamped-free", 1, m_tip=mt)[0]
        est = impl.rayleigh_tip_mass_estimate(EI, RA, L, mt)
        assert wt < prev and est >= wt                           # heavier: lower; Rayleigh: upper bound
        prev = wt
    big = 1e4
    wt = impl.natural_frequencies(EI, RA, L, "clamped-free", 1, m_tip=big)[0]
    assert wt == pytest.approx(np.sqrt(3 * EI / L ** 3 / big), rel=1e-3)
    assert impl.rayleigh_tip_mass_estimate(EI, RA, L, 0.0) == pytest.approx(w0, rel=0.02)


@pytest.mark.parametrize("bc", ["clamped-free", "pinned-pinned", "clamped-clamped"])
def test_mode_shapes_satisfy_ode_bcs_and_orthogonality(bc):
    bl = impl.beta_L(bc, 3)
    x = np.linspace(0, L, 20001)
    dx = x[1] - x[0]
    phis = [impl.mode_shape(x, L, b, bc) for b in bl]
    # Coarser grid for the fourth difference: round-off in cosh - s sinh grows as eps cosh(bL) / h^4, which at
    # h = L/1000 is ~0.8 of the tolerance for the third clamped-clamped mode and platform-dependent; at h = L/400
    # round-off is ~40x smaller and the O(h^2) truncation error stays below 0.07 of the tolerance.
    xc = np.linspace(0, L, 401)
    dc = xc[1] - xc[0]
    for b, ph in zip(bl, phis):
        assert np.max(np.abs(ph)) == pytest.approx(1.0, abs=1e-6)
        pc = impl.mode_shape(xc, L, b, bc)
        d4 = (pc[:-4] - 4 * pc[1:-3] + 6 * pc[2:-2] - 4 * pc[3:-1] + pc[4:]) / dc ** 4
        npt.assert_allclose(d4, (b / L) ** 4 * pc[2:-2], atol=3e-3 * (b / L) ** 4)   # round-off in cosh - cos / h^4
        assert abs(ph[0]) < 1e-10
        if bc != "pinned-pinned":
            assert abs((ph[1] - ph[0]) / dx) < 1e-3 * b / L         # clamped: zero slope
        d2 = lambda i: (ph[i - 1] - 2 * ph[i] + ph[i + 1]) / dx ** 2
        if bc == "clamped-free":
            assert abs(d2(-2)) < 2e-3 * (b / L) ** 2                  # free end: zero moment
        if bc == "pinned-pinned":
            assert abs(ph[-1]) < 1e-10 and abs(d2(1)) < 2e-3 * (b / L) ** 2
    G = np.array([[np.trapezoid(p * q, x) for q in phis] for p in phis])
    off = G - np.diag(np.diag(G))
    assert np.max(np.abs(off)) < 1e-6 * np.min(np.diag(G))


# ---------------------------------------------------------------- finite differences

def test_fd_pinned_matches_its_discrete_closed_form_and_converges_at_second_order():
    for n in (9, 30):
        npt.assert_allclose(impl.fd_pinned_frequencies(EI, RA, L, n, 3), impl.fd_pinned_exact(EI, RA, L, n, 3), rtol=1e-10)
    ex = impl.natural_frequencies(EI, RA, L, "pinned-pinned", 1)[0]
    errs = [abs(impl.fd_pinned_frequencies(EI, RA, L, n, 1)[0] / ex - 1) for n in (9, 19, 39, 79)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 1.9) and np.all(orders < 2.1)
    assert impl.fd_pinned_frequencies(EI, RA, L, 79, 1)[0] < ex     # the central difference underestimates


# ---------------------------------------------------------------- finite elements

def test_element_matrices():
    k, m = impl.beam_element(EI, RA, 0.4)
    npt.assert_allclose(k, k.T)
    npt.assert_allclose(m, m.T)
    npt.assert_allclose(k @ np.array([1.0, 0.0, 1.0, 0.0]), 0.0, atol=1e-6)           # rigid translation
    npt.assert_allclose(k @ np.array([0.0, 1.0, 0.4, 1.0]), 0.0, atol=1e-6)           # rigid rotation
    u = np.array([1.0, 0.0, 1.0, 0.0])
    assert u @ m @ u == pytest.approx(RA * 0.4)                                       # total mass
    assert np.all(np.linalg.eigvalsh(m) > 0)


def test_fem_mass_and_point_mass():
    nodes, K, M = impl.fem_matrices(L, EI, RA, 10, point_masses=[(2.0, 40.0)])
    u = np.zeros(K.shape[0])
    u[0::2] = 1.0
    assert u @ M @ u == pytest.approx(RA * L + 40.0)
    npt.assert_allclose(K, K.T, atol=1e-6)


@pytest.mark.parametrize("bc", ["clamped-free", "pinned-pinned", "clamped-clamped"])
def test_fem_converges_at_fourth_order(bc):
    ex = impl.natural_frequencies(EI, RA, L, bc, 3)
    errs = []
    for n in (8, 16, 32):
        w, _, _ = impl.fem_frequencies(L, EI, RA, n, bc, 3)
        assert np.all(w >= ex * (1 - 1e-12))                       # consistent mass: upper bounds
        errs.append(w[2] / ex[2] - 1)
    orders = impl.observed_order(errs)
    assert np.all(orders > 3.5) and np.all(orders < 4.5)
    w, _, _ = impl.fem_frequencies(L, EI, RA, 32, bc, 3)
    npt.assert_allclose(w, ex, rtol=1e-5)


def test_fem_mode_shapes_match_closed_form():
    w, nodes, W = impl.fem_frequencies(L, EI, RA, 16, "clamped-free", 3)
    bl = impl.beta_L("clamped-free", 3)
    for i in range(3):
        npt.assert_allclose(W[:, i], impl.mode_shape(nodes, L, bl[i], "clamped-free"), atol=2e-3)


def test_fem_with_tip_mass_matches_frequency_equation():
    for mt in (10.0, 40.0):
        ex = impl.natural_frequencies(EI, RA, L, "clamped-free", 2, m_tip=mt)
        w, _, _ = impl.fem_frequencies(L, EI, RA, 20, "clamped-free", 2, point_masses=[(L, mt)])
        npt.assert_allclose(w, ex, rtol=1e-5)


def test_fem_static_cantilever_matches_p05_12():
    n = 8
    nodes, K, M = impl.fem_matrices(L, EI, RA, n)
    f = impl.free_dofs(n, "clamped-free")
    q = 1000.0
    h = L / n
    F = np.zeros(K.shape[0])
    for e in range(n):                                            # consistent nodal loads of a uniform q
        F[2 * e: 2 * e + 4] += q * np.array([h / 2, h * h / 12, h / 2, -h * h / 12])
    u = np.zeros_like(F)
    u[f] = np.linalg.solve(K[np.ix_(f, f)], F[f])
    w_exact = q * nodes ** 2 * (6 * L ** 2 - 4 * L * nodes + nodes ** 2) / (24 * EI)
    npt.assert_allclose(u[0::2], w_exact, rtol=1e-9, atol=1e-15)   # cubic elements: nodally exact


def test_variable_properties_and_tank():
    w0 = impl.natural_frequencies(EI, RA, L, "clamped-free", 2)
    wt, _, _ = impl.fem_frequencies(L, EI, RA, 40, "clamped-free", 2, point_masses=[(2.0, 40.0)])
    assert np.all(wt < w0)                                        # added mass lowers every frequency
    assert wt[1] / w0[1] < wt[0] / w0[0]                           # a tank at 2 m sits near mode 2's antinode region
    # stiffness that tapers outboard with mass that tapers outboard
    wv, _, _ = impl.fem_frequencies(L, lambda x: EI * (1 - 0.7 * x / L), lambda x: RA * 1.5 * (1 - 2 / 3 * x / L),
                                    80, "clamped-free", 1)
    wv2, _, _ = impl.fem_frequencies(L, lambda x: EI * (1 - 0.7 * x / L), lambda x: RA * 1.5 * (1 - 2 / 3 * x / L),
                                     160, "clamped-free", 1)
    assert wv2[0] == pytest.approx(wv[0], rel=1e-4)                # converged
    # same EI everywhere scaled by 2: frequencies by sqrt(2)
    w2, _, _ = impl.fem_frequencies(L, 2 * EI, RA, 16, "clamped-free", 2)
    w1, _, _ = impl.fem_frequencies(L, EI, RA, 16, "clamped-free", 2)
    npt.assert_allclose(w2, np.sqrt(2) * w1, rtol=1e-10)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "beam_modes.png").stat().st_size > 10_000
