"""Public tests for P06.14 Multi-DOF Vibration & Modes. They assert structural dynamics, not magic numbers:
the chain matrices are symmetric with the right static behaviour; the Cholesky eigen-solution matches
scipy.linalg.eigh(K, M), is mass-normalised and M- and K-orthogonal; a two-mass closed form; the Rayleigh
quotient is exact at a mode and an upper bound otherwise; effective masses sum to the total mass; Rayleigh
damping hits its two targets and is diagonalised by the modes while a single damper is not; modal
superposition equals direct integration for the free and forced response and the direct inverse for the
FRF; truncation errors fall as modes are added; the interface apparent mass tends to the total mass and
agrees with the direct inertia force.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.linalg import eigh

from conftest import load_impl

impl = load_impl(__file__)

MS = np.array([150.0, 600.0, 350.0])
KS = np.array([6.0e8, 4.0e8, 2.5e8])


def _stack(zeta=0.02):
    M, K = impl.chain_matrices(MS, KS)
    w, Phi = impl.modes(K, M)
    a, b = impl.rayleigh_coefficients(w[0], w[1], zeta, zeta)
    C = a * M + b * K
    z, _ = impl.modal_damping(Phi, C, w)
    return M, K, C, w, Phi, z


# ---------------------------------------------------------------- matrices and modes

def test_chain_matrices():
    M, K = impl.chain_matrices(MS, KS)
    npt.assert_allclose(M, np.diag(MS))
    npt.assert_allclose(K, K.T)
    assert np.all(np.linalg.eigvalsh(K) > 0)                      # fixed base: no rigid-body mode
    # a unit force on the top mass stretches every spring by 1/k_i: static flexibility
    x = np.linalg.solve(K, np.array([0.0, 0.0, 1.0]))
    npt.assert_allclose(x, np.cumsum(1 / KS), rtol=1e-12)


def test_modes_match_eigh_and_are_orthogonal():
    M, K = impl.chain_matrices(MS, KS)
    w, Phi = impl.modes(K, M)
    lam, _ = eigh(K, M)
    npt.assert_allclose(w, np.sqrt(lam), rtol=1e-12)
    assert np.all(np.diff(w) > 0)
    npt.assert_allclose(Phi.T @ M @ Phi, np.eye(3), atol=1e-12)
    npt.assert_allclose(Phi.T @ K @ Phi, np.diag(w ** 2), atol=1e-9 * w[-1] ** 2)
    npt.assert_allclose(K @ Phi, M @ Phi @ np.diag(w ** 2), rtol=1e-9, atol=1e-6)
    # first mode has no node: all masses move the same way
    assert np.all(Phi[:, 0] > 0)


def test_two_mass_closed_form():
    m1, m2, k1, k2 = 2.0, 1.0, 300.0, 100.0
    M, K = impl.chain_matrices([m1, m2], [k1, k2])
    w, _ = impl.modes(K, M)
    # det(K - w^2 M) = m1 m2 w^4 - (m1 k2 + m2 (k1 + k2)) w^2 + k1 k2 = 0
    a, b, c = m1 * m2, -(m1 * k2 + m2 * (k1 + k2)), k1 * k2
    npt.assert_allclose(w ** 2, np.sort(np.roots([a, b, c])), rtol=1e-12)


def test_rayleigh_quotient():
    M, K = impl.chain_matrices(MS, KS)
    w, Phi = impl.modes(K, M)
    for i in range(3):
        assert impl.rayleigh_quotient(K, M, Phi[:, i]) == pytest.approx(w[i] ** 2, rel=1e-12)
    rng = np.random.default_rng(5)
    for _ in range(20):
        x = rng.normal(size=3)
        R = impl.rayleigh_quotient(K, M, x)
        assert w[0] ** 2 * (1 - 1e-12) <= R <= w[-1] ** 2 * (1 + 1e-12)
    # the static deflection under gravity is a good trial vector (within 1 % in frequency here)
    xs = np.linalg.solve(K, MS)
    assert np.sqrt(impl.rayleigh_quotient(K, M, xs)) == pytest.approx(w[0], rel=1e-2)


def test_effective_masses_sum_to_total():
    M, K = impl.chain_matrices(MS, KS)
    w, Phi = impl.modes(K, M)
    G, meff = impl.participation(Phi, M)
    assert np.sum(meff) == pytest.approx(MS.sum(), rel=1e-12)
    assert np.all(meff >= 0) and meff[0] > 0.5 * MS.sum()
    # a different base direction vector: effective masses sum to r^T M r
    r = np.array([1.0, 2.0, -1.0])
    _, me2 = impl.participation(Phi, M, r)
    assert np.sum(me2) == pytest.approx(r @ M @ r)


# ---------------------------------------------------------------- damping

def test_rayleigh_damping_targets_and_proportionality():
    M, K, C, w, Phi, z = _stack(0.02)
    npt.assert_allclose(z[:2], 0.02, rtol=1e-12)
    a, b = impl.rayleigh_coefficients(w[0], w[1], 0.02, 0.02)
    assert z[2] == pytest.approx(a / (2 * w[2]) + b * w[2] / 2)
    assert z[2] > 0.02                                             # stiffness term dominates at high frequency
    _, coup = impl.modal_damping(Phi, C, w)
    assert coup < 1e-10
    # a single damper between the service module and the payload is not proportional
    Cd = np.zeros((3, 3))
    cd = 2e4
    Cd[1, 1] = Cd[2, 2] = cd
    Cd[1, 2] = Cd[2, 1] = -cd
    _, coup_d = impl.modal_damping(Phi, Cd, w)
    assert coup_d > 0.05


# ---------------------------------------------------------------- superposition vs direct integration

def test_modal_free_response_equals_direct_integration():
    M, K, C, w, Phi, z = _stack()
    t = np.linspace(0, 0.1, 2001)
    x0 = np.linalg.solve(K, MS * 5 * 9.80665)
    v0 = np.array([0.0, 0.01, -0.02])
    xm = impl.modal_free_response(t, x0, v0, w, Phi, M, z)
    xd = impl.direct_integration(M, C, K, lambda tt: np.zeros(3), x0, v0, t)
    npt.assert_allclose(xm[:, 0], x0, atol=1e-15)
    npt.assert_allclose(xm, xd, atol=1e-9 * np.max(np.abs(x0)))


def test_modal_forced_response_and_truncation():
    M, K, C, w, Phi, z = _stack()
    t = np.linspace(0, 0.12, 2401)
    f = lambda tt: np.array([0.0, 0.0, 3000.0 * np.sin(w[0] * tt)])
    xd = impl.direct_integration(M, C, K, f, np.zeros(3), np.zeros(3), t)
    errs = [np.max(np.abs(impl.modal_integration(w, Phi, z, f, t, n) - xd)) for n in (1, 2, 3)]
    assert errs[2] < 1e-8 * np.max(np.abs(xd))                     # all modes: exact
    assert errs[0] > errs[1] > errs[2]                             # truncation error falls with modes
    assert errs[0] < 0.05 * np.max(np.abs(xd))                     # driven at f1: mode 1 dominates


def test_modal_frf_equals_direct_inverse():
    M, K, C, w, Phi, z = _stack()
    for om in (0.3 * w[0], w[0], 0.5 * (w[1] + w[2]), w[2], 2 * w[2]):
        npt.assert_allclose(impl.frf_modal(om, w, Phi, z), impl.frf_direct(om, M, C, K), rtol=1e-9, atol=1e-18)
    # static limit: K^-1
    npt.assert_allclose(impl.frf_modal(0.0, w, Phi, z), np.linalg.inv(K), rtol=1e-10)
    # symmetric (reciprocity)
    H = impl.frf_direct(1.3 * w[0], M, C, K)
    npt.assert_allclose(H, H.T, rtol=1e-12)


def test_apparent_mass_limits_and_direct_check():
    M, K, C, w, Phi, z = _stack()
    _, meff = impl.participation(Phi, M)
    assert abs(impl.apparent_mass(1e-6 * w[0], w, z, meff)[0]) == pytest.approx(MS.sum(), rel=1e-9)
    for om in (0.5 * w[0], w[0], w[1], 1.5 * w[2]):
        _, F = impl.base_sine_accelerations(om, 1.0, M, C, K)
        assert impl.apparent_mass(om, w, z, meff)[0] == pytest.approx(F, rel=1e-9)
    # truncated to mode 1 with the residual mass: within 1 % at f1, exact at zero frequency
    a1 = impl.apparent_mass(np.array([1e-6, w[0]]), w, z, meff, n_modes=1)
    assert abs(a1[0]) == pytest.approx(MS.sum(), rel=1e-9)
    assert abs(a1[1]) == pytest.approx(abs(impl.apparent_mass(w[0], w, z, meff)[0]), rel=1e-2)


def test_base_sine_accelerations_low_frequency_is_rigid():
    M, K, C, w, Phi, z = _stack()
    acc, F = impl.base_sine_accelerations(1e-3, 9.80665, M, C, K)
    npt.assert_allclose(acc, 9.80665, rtol=1e-8)
    assert F.real == pytest.approx(MS.sum() * 9.80665, rel=1e-8)


def test_uniform_fixed_free_chain_closed_form():
    # n equal masses and springs, fixed at one end: w_j = 2 sqrt(k/m) sin((2j - 1) pi / (2 (2n + 1)))
    n, m, k = 8, 3.0, 1.2e4
    M, K = impl.chain_matrices(np.full(n, m), np.full(n, k))
    w, Phi = impl.modes(K, M)
    j = np.arange(1, n + 1)
    npt.assert_allclose(w, 2 * np.sqrt(k / m) * np.sin((2 * j - 1) * np.pi / (2 * (2 * n + 1))), rtol=1e-11)
    # mode j has j - 1 sign changes (nodes)
    for jj in range(n):
        assert np.sum(np.diff(np.sign(Phi[:, jj])) != 0) == jj


def test_undamped_modal_free_response_conserves_energy():
    M, K = impl.chain_matrices(MS, KS)
    w, Phi = impl.modes(K, M)
    t = np.linspace(0, 0.05, 5001)
    x0 = np.array([1e-4, -2e-4, 3e-4])
    x = impl.modal_free_response(t, x0, np.zeros(3), w, Phi, M, np.zeros(3))
    dt = t[1] - t[0]
    v = np.gradient(x, dt, axis=1, edge_order=2)
    E = 0.5 * np.einsum("it,ij,jt->t", v, M, v) + 0.5 * np.einsum("it,ij,jt->t", x, K, x)
    assert np.ptp(E[5:-5]) < 1e-3 * E[0]                         # finite-difference velocities: O(dt^2)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "modal_analysis.png").stat().st_size > 10_000
