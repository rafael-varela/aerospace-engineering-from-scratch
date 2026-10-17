"""Public tests for P02.04 Systems of ODEs & Phase Portraits. They assert properties, with numpy.linalg and scipy as oracles.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from conftest import load_impl

impl = load_impl(__file__)

RNG = np.random.default_rng(204)
A_SP = impl.short_period_matrix()
A_PH = impl.phugoid_matrix()
A_SADDLE = impl.short_period_matrix(m_alpha=2.0)
A_CENTER = np.array([[-0.8, 1.0], [-4.0, 0.8]])
A_STAR = -1.5 * np.eye(2)
A_DEGEN = impl.second_order_matrix(1.0, 4.0, 4.0)
A_NODE = impl.second_order_matrix(1.0, 5.0, 4.0)
CASES = {
    "stable spiral": A_SP, "saddle": A_SADDLE, "center": A_CENTER, "stable star": A_STAR,
    "stable degenerate node": A_DEGEN, "stable node": A_NODE,
    "unstable node": np.array([[2.0, 1.0], [0.0, 3.0]]), "unstable spiral": np.array([[0.5, -2.0], [2.0, 0.5]]),
    "unstable star": 0.7 * np.eye(2), "unstable degenerate node": np.array([[1.0, 1.0], [0.0, 1.0]]),
    "non-isolated": np.array([[1.0, 2.0], [2.0, 4.0]]),
}


# ---------------------------------------------------------------- matrices, trace, determinant, eigenstructure

def test_trace_determinant_and_the_characteristic_polynomial():
    for A in (A_SP, A_PH, A_SADDLE, RNG.normal(size=(2, 2))):
        tau, delta = impl.trace_determinant(A)
        assert tau == pytest.approx(np.trace(A)) and delta == pytest.approx(np.linalg.det(A))
        lam = np.linalg.eigvals(A)
        assert lam.sum() == pytest.approx(tau) and np.prod(lam) == pytest.approx(delta)
        for l in lam:
            assert abs(l * l - tau * l + delta) < 1e-10
    with pytest.raises(ValueError):
        impl.trace_determinant(np.eye(3))


def test_eigen_structure_matches_numpy_and_satisfies_A_v_equals_lambda_v():
    for A in (A_SP, A_SADDLE, A_NODE, A_CENTER, A_PH, RNG.normal(size=(2, 2)), RNG.normal(size=(2, 2))):
        lam, V = impl.eigen_structure(A)
        npt.assert_allclose(np.sort_complex(lam.astype(complex)), np.sort_complex(np.linalg.eigvals(A)), rtol=1e-9, atol=1e-12)
        for k in range(2):
            npt.assert_allclose(A @ V[:, k], lam[k] * V[:, k], atol=1e-10)
            assert np.linalg.norm(V[:, k]) == pytest.approx(1.0)
    # Real eigenvalues come back as a real array; complex ones as complex.
    assert np.isrealobj(impl.eigen_structure(A_SADDLE)[0]) and np.iscomplexobj(impl.eigen_structure(A_SP)[0])
    # Repeated eigenvalue: one eigenvector direction (degenerate) or any direction (star).
    lam, V = impl.eigen_structure(A_DEGEN)
    npt.assert_allclose(lam, [-2.0, -2.0])
    assert abs(abs(V[:, 0] @ V[:, 1]) - 1.0) < 1e-12                                   # same line
    lam, V = impl.eigen_structure(A_STAR)
    npt.assert_allclose(lam, [-1.5, -1.5])
    npt.assert_allclose(A_STAR @ V, -1.5 * V)


# ---------------------------------------------------------------- classification

@pytest.mark.parametrize("name", list(CASES))
def test_classification_by_trace_and_determinant(name):
    assert impl.classify(CASES[name]) == name


def test_classification_is_consistent_with_the_eigenvalues():
    for _ in range(40):
        A = RNG.normal(size=(2, 2)) * 2.0
        lam = np.linalg.eigvals(A)
        kind = impl.classify(A)
        if np.iscomplexobj(lam) and np.any(np.abs(lam.imag) > 1e-9):
            assert "spiral" in kind or kind == "center"
            assert kind.startswith("stable") == (lam.real[0] < 0) or kind == "center"
        elif np.prod(lam.real) < 0:
            assert kind == "saddle"
        else:
            assert "node" in kind or "star" in kind
            assert kind.startswith("unstable") == (lam.real.max() > 0)


def test_modal_parameters_of_the_aircraft_modes():
    sp = impl.modal_parameters(A_SP)
    tau, delta = impl.trace_determinant(A_SP)
    assert sp["omega_n"] == pytest.approx(np.sqrt(delta))
    assert sp["zeta"] == pytest.approx(-tau / (2 * np.sqrt(delta)))
    lam = np.linalg.eigvals(A_SP)
    assert sp["omega_d"] == pytest.approx(abs(lam[0].imag))
    assert sp["t_half"] == pytest.approx(np.log(2) / abs(lam[0].real))
    assert sp["period"] == pytest.approx(2 * np.pi / abs(lam[0].imag))
    # Lanchester's phugoid: omega_n = sqrt(2) g / V0 and zeta = 1 / (sqrt(2) L/D), for any V0 and L/D.
    for v0, ld in ((230.0, 17.0), (60.0, 10.0), (250.0, 20.0)):
        ph = impl.modal_parameters(impl.phugoid_matrix(v0, ld))
        assert ph["omega_n"] == pytest.approx(np.sqrt(2) * 9.81 / v0)
        assert ph["zeta"] == pytest.approx(1 / (np.sqrt(2) * ld))
        assert ph["period"] == pytest.approx(2 * np.pi / (ph["omega_n"] * np.sqrt(1 - ph["zeta"] ** 2)))
    with pytest.raises(ValueError):
        impl.modal_parameters(A_SADDLE)
    with pytest.raises(ValueError):
        impl.phugoid_matrix(-1.0, 10.0)
    # The second-order matrix reproduces the P02.02 roots.
    A2 = impl.second_order_matrix(1200.0, 2 * 0.3 * np.sqrt(150e3 * 1200), 150e3)
    npt.assert_allclose(np.sort_complex(np.linalg.eigvals(A2)), np.sort_complex(np.roots([1200.0, 2 * 0.3 * np.sqrt(150e3 * 1200), 150e3])))


# ---------------------------------------------------------------- solutions

@pytest.mark.parametrize("A", [A_SP, A_SADDLE, A_NODE, A_CENTER, A_STAR, A_DEGEN, A_PH])
def test_closed_form_solution_matches_expm_and_the_ode(A):
    x0 = np.array([0.3, -0.2])
    t_end = 300.0 if A is A_PH else 4.0
    t = np.linspace(0.0, t_end, 61)
    x = impl.solve_linear(A, x0, t)
    assert x.shape == (61, 2)
    npt.assert_allclose(x[0], x0, atol=1e-12)
    for i, tt in enumerate(t):
        npt.assert_allclose(x[i], expm(A * tt) @ x0, rtol=1e-8, atol=1e-10)
    # dx/dt = A x by central differences.
    dt = 1e-5 * t_end
    xp = impl.solve_linear(A, x0, t[5:-5] + dt)
    xm = impl.solve_linear(A, x0, t[5:-5] - dt)
    npt.assert_allclose((xp - xm) / (2 * dt), x[5:-5] @ A.T, rtol=1e-6, atol=1e-8)


def test_rk4_converges_at_fourth_order_and_solve_ivp_agrees():
    x0 = np.array([np.radians(5.0), 0.0])
    errs = []
    for n in (50, 100, 200, 400):
        ts, xs = impl.rk4_trajectory(A_SP, x0, 8.0, n)
        errs.append(np.max(np.abs(xs - impl.solve_linear(A_SP, x0, ts))))
    rates = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(np.abs(rates - 4.0) < 0.3), rates
    sol = solve_ivp(lambda t, y: A_SP @ y, (0.0, 8.0), x0, rtol=1e-10, atol=1e-12, t_eval=ts)
    npt.assert_allclose(sol.y.T, impl.solve_linear(A_SP, x0, ts), atol=1e-8)
    with pytest.raises(ValueError):
        impl.rk4_trajectory(A_SP, x0, 1.0, 0)


def test_saddle_trajectories_follow_the_eigen_directions():
    lam, V = impl.eigen_structure(A_SADDLE)
    assert lam[0] < 0 < lam[1]
    t = np.linspace(0.0, 3.0, 31)
    # A start on the stable eigenvector decays to the origin along that line.
    xs = impl.solve_linear(A_SADDLE, 0.2 * V[:, 0], t)
    npt.assert_allclose(xs, np.outer(0.2 * np.exp(lam[0] * t), V[:, 0]), atol=1e-12)
    assert np.linalg.norm(xs[-1]) < 1e-3
    # A start on the unstable eigenvector grows along that line.
    xu = impl.solve_linear(A_SADDLE, 0.01 * V[:, 1], t)
    npt.assert_allclose(xu, np.outer(0.01 * np.exp(lam[1] * t), V[:, 1]), atol=1e-12)
    # Any other start ends up parallel to the unstable direction.
    xg = impl.solve_linear(A_SADDLE, np.array([0.1, -0.3]), np.array([12.0]))[0]
    cos = abs(xg @ V[:, 1]) / np.linalg.norm(xg)
    assert cos > 1 - 1e-6


def test_center_conserves_the_quadratic_invariant_and_spiral_does_not():
    t = np.linspace(0.0, 20.0, 2001)
    x = impl.solve_linear(A_CENTER, np.array([0.3, 0.1]), t)
    Q = impl.center_invariant(A_CENTER, x)
    npt.assert_allclose(Q, Q[0], rtol=1e-9)
    # For the undamped oscillator the invariant is -(2/m) times the energy.
    m, k = 2.0, 8.0
    A = impl.second_order_matrix(m, 0.0, k)
    x = impl.solve_linear(A, np.array([0.1, 0.5]), t)
    E = 0.5 * m * x[:, 1] ** 2 + 0.5 * k * x[:, 0] ** 2
    npt.assert_allclose(impl.center_invariant(A, x), -2 * E / m, rtol=1e-9)
    # The same quadratic along the damped short period decays monotonically in magnitude (it is not conserved).
    xs = impl.solve_linear(A_SP, np.array([0.3, 0.1]), np.linspace(0.0, 4.0, 401))
    Qs = impl.center_invariant(A_SP, xs)
    assert abs(Qs[-1]) < 0.01 * abs(Qs[0])


def test_spiral_decay_rate_and_period_match_the_eigenvalues():
    x0 = np.array([np.radians(5.0), 0.0])
    t = np.linspace(0.0, 10.0, 100001)
    x = impl.solve_linear(A_SP, x0, t)
    mp = impl.modal_parameters(A_SP)
    # Zero crossings of alpha are half a damped period apart.
    zc = t[np.where(np.diff(np.sign(x[:, 0])) != 0)[0]]
    npt.assert_allclose(np.diff(zc), mp["period"] / 2, rtol=1e-3)
    # The radius in a suitably scaled plane decays as e^{Re(lambda) t}: check |x| at successive periods.
    r = np.linalg.norm(x, axis=1)
    i0, i1 = np.argmin(np.abs(t - 2.0)), np.argmin(np.abs(t - 2.0 - mp["period"]))
    assert r[i1] / r[i0] == pytest.approx(np.exp(-mp["zeta"] * mp["omega_n"] * mp["period"]), rel=1e-3)


def test_vector_field_nullclines_and_portrait_helpers():
    g = np.linspace(-1.0, 1.0, 5)
    X, Y = np.meshgrid(g, g)
    U, V = impl.vector_field(A_SP, X, Y)
    npt.assert_allclose(U, A_SP[0, 0] * X + A_SP[0, 1] * Y)
    npt.assert_allclose(V, A_SP[1, 0] * X + A_SP[1, 1] * Y)
    s1, s2 = impl.nullcline_slopes(A_SP)
    # On the x1-nullcline the field is vertical, on the x2-nullcline horizontal.
    p1 = np.array([1.0, s1])
    assert (A_SP @ p1)[0] == pytest.approx(0.0, abs=1e-12) and (A_SP @ p1)[1] != 0.0
    p2 = np.array([1.0, s2])
    assert (A_SP @ p2)[1] == pytest.approx(0.0, abs=1e-12) and (A_SP @ p2)[0] != 0.0
    assert impl.nullcline_slopes(np.array([[1.0, 0.0], [1.0, 0.0]])) == (float("inf"), float("inf"))
    trajs = impl.phase_portrait(A_SP, np.array([[0.2, 0.0], [0.0, 0.2], [-0.2, 0.1]]), 5.0, 100)
    assert len(trajs) == 3 and all(tr.shape == (101, 2) for tr in trajs)
    for tr, x0 in zip(trajs, [[0.2, 0.0], [0.0, 0.2], [-0.2, 0.1]]):
        npt.assert_allclose(tr[0], x0, atol=1e-15)
        assert np.linalg.norm(tr[-1]) < 0.01 * np.linalg.norm(tr[0])                   # stable: everything goes home


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "phase_portraits.png").stat().st_size > 1000
