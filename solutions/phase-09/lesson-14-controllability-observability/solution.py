"""P09.14 Controllability & Observability: reference solution.

Rank tests from scratch: the controllability matrix [B, AB, ..., A^{n-1}B] and the observability matrix,
the Popov-Belevitch-Hautus (PBH) eigenvalue test that says *which* mode is lost, the controllability and
observability Gramians from Lyapunov equations (solved here by Kronecker vectorisation) and over a finite
horizon (Van Loan's block exponential), and minimum-energy steering from the finite-horizon Gramian.

The aerospace case is a small satellite with reaction wheels, linearised about an inertial attitude hold
(three decoupled double integrators J_i theta_i'' = tau_i): an orthogonal three-wheel set loses an axis when
one wheel fails, a four-wheel pyramid does not; a star tracker that measures only two axes leaves the
boresight angle unobservable even with gyros, while a full star tracker keeps everything observable after
a gyro failure.

Conventions:
    * State x = (theta_x, theta_y, theta_z, omega_x, omega_y, omega_z) [rad, rad/s] for the satellite.
    * Wheel torques u [N m]; B = [[0], [J^-1 W]] with W the 3 x n_wheels matrix of wheel spin axes.
    * Ranks use the SVD with a relative tolerance (tol x largest singular value).

Run:  python solution.py          (prints the checks, writes controllability.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-09/lesson-14-controllability-observability
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# ---------------------------------------------------------------- given

def expm(M: np.ndarray, n_terms: int = 20) -> np.ndarray:
    """Return e^M by scaling and squaring a truncated Taylor series (the method of P02.05)."""
    M = np.asarray(M, dtype=float)
    norm = np.linalg.norm(M, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    X = M / 2 ** s
    term = np.eye(len(M))
    E = np.eye(len(M))
    for k in range(1, n_terms):
        term = term @ X / k
        E = E + term
    for _ in range(s):
        E = E @ E
    return E


def numerical_rank(M: np.ndarray, tol: float = 1e-9) -> int:
    """Return the number of singular values of M above tol x the largest one (0 for a zero matrix)."""
    sv = np.linalg.svd(np.atleast_2d(M), compute_uv=False)
    if sv.size == 0 or sv[0] == 0.0:
        return 0
    return int(np.sum(sv > tol * sv[0]))


def satellite_inertia() -> np.ndarray:
    """Return a representative small-satellite inertia matrix diag(4, 5, 3) kg m^2 (principal axes)."""
    return np.diag([4.0, 5.0, 3.0])


def wheel_axes(config: str = "orthogonal", failed: tuple[int, ...] = ()) -> np.ndarray:
    """Return W (3 x n): unit spin axes of the working wheels, columns of failed wheels removed.

    'orthogonal': three wheels on x, y, z. 'pyramid': four wheels with spin axes at 54.74 deg from +z
    (cos^2 = 1/3) and azimuths 0, 90, 180, 270 deg, the classic redundant arrangement.
    """
    if config == "orthogonal":
        W = np.eye(3)
    elif config == "pyramid":
        b = np.arccos(1.0 / np.sqrt(3.0))
        az = np.radians([0.0, 90.0, 180.0, 270.0])
        W = np.vstack([np.sin(b) * np.cos(az), np.sin(b) * np.sin(az), np.cos(b) * np.ones(4)])
    else:
        raise ValueError(config)
    keep = [i for i in range(W.shape[1]) if i not in failed]
    return W[:, keep]


def satellite_model(J: np.ndarray, W: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (A, B) of the linearised attitude hold: theta' = omega, J omega' = W u (no gyroscopic or orbital terms)."""
    A = np.zeros((6, 6))
    A[:3, 3:] = np.eye(3)
    B = np.zeros((6, W.shape[1]))
    B[3:, :] = np.linalg.solve(J, W)
    return A, B


def short_period_model() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (A, B, C) of the stable short-period pitch model of P09.11/13: x = (alpha, q), input elevator, output q."""
    A = np.array([[-1.2, 1.0], [-6.0, -2.5]])
    return A, np.array([[0.0], [12.0]]), np.array([[0.0, 1.0]])


# ---------------------------------------------------------------- rank tests

def ctrb(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Return the controllability matrix [B, AB, A^2 B, ..., A^{n-1} B] (n x n m)."""
    n = A.shape[0]
    blocks = [B]
    for _ in range(n - 1):
        blocks.append(A @ blocks[-1])
    return np.hstack(blocks)


def obsv(A: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Return the observability matrix [C; CA; ...; C A^{n-1}] (n p x n)."""
    n = A.shape[0]
    blocks = [C]
    for _ in range(n - 1):
        blocks.append(blocks[-1] @ A)
    return np.vstack(blocks)


def is_controllable(A: np.ndarray, B: np.ndarray, tol: float = 1e-9) -> bool:
    """Return True when rank ctrb(A, B) = n."""
    return numerical_rank(ctrb(A, B), tol) == A.shape[0]


def is_observable(A: np.ndarray, C: np.ndarray, tol: float = 1e-9) -> bool:
    """Return True when rank obsv(A, C) = n."""
    return numerical_rank(obsv(A, C), tol) == A.shape[0]


def pbh_uncontrollable_modes(A: np.ndarray, B: np.ndarray, tol: float = 1e-9) -> list[complex]:
    """Return the eigenvalues lambda of A for which rank [lambda I - A, B] < n (PBH test); empty list = controllable.

    Repeated eigenvalues are tested once each (distinct values within 1e-8 relative are merged).
    """
    n = A.shape[0]
    lam = np.linalg.eigvals(A)
    distinct: list[complex] = []
    for li in lam:
        if not any(abs(li - d) <= 1e-8 * max(1.0, abs(li)) for d in distinct):
            distinct.append(complex(li))
    scale = max(np.linalg.norm(A, 2), np.linalg.norm(B, 2), 1e-300)
    bad = []
    for li in distinct:
        M = np.hstack([li * np.eye(n) - A, B]) / scale
        sv = np.linalg.svd(M, compute_uv=False)
        if np.sum(sv > tol) < n:
            bad.append(li)
    return bad


def pbh_unobservable_modes(A: np.ndarray, C: np.ndarray, tol: float = 1e-9) -> list[complex]:
    """Return the eigenvalues for which rank [lambda I - A; C] < n (by duality: PBH controllability of (A^T, C^T))."""
    return pbh_uncontrollable_modes(A.T, C.T, tol)


def unobservable_subspace(A: np.ndarray, C: np.ndarray, tol: float = 1e-9) -> np.ndarray:
    """Return an orthonormal basis (columns) of the null space of the observability matrix: states the outputs never reveal."""
    O = obsv(A, C)
    _, sv, Vt = np.linalg.svd(O)
    r = numerical_rank(O, tol)
    return Vt[r:].T


# ---------------------------------------------------------------- Gramians

def lyap(A: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Return X solving A X + X A^T + Q = 0 by Kronecker vectorisation: (I kron A + A kron I) vec(X) = -vec(Q).

    Column-major vec: vec(A X) = (I kron A) vec X and vec(X A^T) = (A kron I) vec X. O(n^6) work: fine for n <= ~20.
    """
    n = A.shape[0]
    I = np.eye(n)
    K = np.kron(I, A) + np.kron(A, I)
    x = np.linalg.solve(K, -Q.reshape(-1, order="F"))
    X = x.reshape((n, n), order="F")
    return 0.5 * (X + X.T)


def ctrb_gramian(A: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Return W = int_0^inf e^{At} B B^T e^{A^T t} dt from A W + W A^T + B B^T = 0. Raise ValueError unless A is stable."""
    if np.max(np.linalg.eigvals(A).real) >= 0:
        raise ValueError("the infinite-horizon Gramian needs every eigenvalue of A in the open left half plane")
    return lyap(A, B @ B.T)


def obsv_gramian(A: np.ndarray, C: np.ndarray) -> np.ndarray:
    """Return M = int_0^inf e^{A^T t} C^T C e^{At} dt from A^T M + M A + C^T C = 0 (the dual of ``ctrb_gramian``)."""
    return ctrb_gramian(A.T, C.T)


def ctrb_gramian_finite(A: np.ndarray, B: np.ndarray, T: float) -> np.ndarray:
    """Return W(T) = int_0^T e^{At} B B^T e^{A^T t} dt by Van Loan's method (valid for any A, stable or not).

    With H = [[-A, B B^T], [0, A^T]], e^{H T} = [[F11, F12], [0, F22]] and W(T) = F22^T F12.
    """
    n = A.shape[0]
    H = np.zeros((2 * n, 2 * n))
    H[:n, :n] = -A
    H[:n, n:] = B @ B.T
    H[n:, n:] = A.T
    E = expm(H * T)
    W = E[n:, n:].T @ E[:n, n:]
    return 0.5 * (W + W.T)


def min_energy_input(A: np.ndarray, B: np.ndarray, x0: np.ndarray, xf: np.ndarray, T: float):
    """Return (u, energy): the input u(t) of least int_0^T |u|^2 dt that takes x0 to xf at time T, and that energy.

    u(t) = B^T e^{A^T (T - t)} W(T)^{-1} (xf - e^{AT} x0) and energy = d^T W(T)^{-1} d with d = xf - e^{AT} x0.
    ``u`` is a function of a scalar or array t returning shape (m,) or (len(t), m). Raise ValueError if W(T) is singular.
    """
    W = ctrb_gramian_finite(A, B, T)
    if numerical_rank(W, 1e-12) < A.shape[0]:
        raise ValueError("W(T) is singular: the target is not reachable for every x0, xf (uncontrollable)")
    d = np.asarray(xf, dtype=float) - expm(A * T) @ np.asarray(x0, dtype=float)
    lam = np.linalg.solve(W, d)

    def u(t):
        ts = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.array([B.T @ expm(A.T * (T - ti)) @ lam for ti in ts])
        return out[0] if np.ndim(t) == 0 else out

    return u, float(d @ lam)


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Rank tests, PBH and Gramians for the reaction-wheel satellite; save controllability.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    J = satellite_inertia()
    print("representative small satellite, J = diag(4, 5, 3) kg m^2, linearised inertial attitude hold")
    for config, failed in (("orthogonal", ()), ("orthogonal", (2,)), ("pyramid", ()), ("pyramid", (0,)),
                           ("pyramid", (3,)), ("pyramid", (0, 2)), ("pyramid", (0, 1))):
        W = wheel_axes(config, failed)
        A, B = satellite_model(J, W)
        r = numerical_rank(ctrb(A, B))
        bad = pbh_uncontrollable_modes(A, B)
        WT = ctrb_gramian_finite(A, B, 20.0)
        ev = np.linalg.eigvalsh(WT)
        print(f"  {config:10s} failed {str(failed):7s}: rank ctrb = {r}/6, PBH uncontrollable eigenvalues {np.round(bad, 6).tolist()},"
              f" sigma_min(W W^T) = {np.linalg.svd(W @ W.T, compute_uv=False).min():.4f},"
              f" lambda_min W(20 s) = {max(ev.min(), 0.0):.3e}")

    A, B = satellite_model(J, wheel_axes("orthogonal"))
    print("\nsensors (orthogonal wheels):")
    I6 = np.eye(6)
    sensors = {
        "3-axis tracker + 3 gyros": I6,
        "3-axis tracker + gyros y, z (x gyro failed)": I6[[0, 1, 2, 4, 5]],
        "2-axis tracker (boresight z) + 3 gyros": I6[[0, 1, 3, 4, 5]],
        "3 gyros only": I6[3:],
        "2-axis tracker only": I6[:2],
    }
    for name, C in sensors.items():
        N = unobservable_subspace(A, C)
        print(f"  {name:44s}: rank obsv = {numerical_rank(obsv(A, C))}/6, PBH unobservable {np.round(pbh_unobservable_modes(A, C), 6).tolist()},"
              f" unobservable directions {[np.round(N[:, k], 3).tolist() for k in range(N.shape[1])]}")

    Asp, Bsp, Csp = short_period_model()
    Wc = ctrb_gramian(Asp, Bsp)
    tq = np.linspace(0.0, 15.0, 15001)
    step = expm(Asp * (tq[1] - tq[0]))
    Phi = np.eye(2)
    integrand = np.empty((len(tq), 2, 2))
    for k in range(len(tq)):                     # e^{A t_k} = (e^{A h})^k: one exponential, then products
        integrand[k] = Phi @ Bsp @ Bsp.T @ Phi.T
        Phi = step @ Phi
    Wq = np.trapezoid(integrand, tq, axis=0)
    print(f"\nshort-period Gramian (Lyapunov) = {np.array2string(Wc, precision=5).replace(chr(10), '')};"
          f" max |Lyapunov - quadrature to 15 s| = {np.abs(Wc - Wq).max():.2e}")

    Jz = J[2, 2]
    dth = np.radians(10.0)
    tau_max = 0.01
    print(f"\nminimum-energy rest-to-rest slew of {np.degrees(dth):.0f} deg about z (J = {Jz} kg m^2), wheel limit {tau_max * 1e3:.0f} mN m:")
    A1 = np.array([[0.0, 1.0], [0.0, 0.0]])
    B1 = np.array([[0.0], [1.0 / Jz]])
    for T in (10.0, 20.0, 40.0):
        u, E = min_energy_input(A1, B1, np.zeros(2), np.array([dth, 0.0]), T)
        tt = np.linspace(0.0, T, 2001)
        peak = np.abs(u(tt)).max()
        print(f"  T = {T:4.0f} s: peak torque {1e3 * peak:.3f} mN m (6 J dtheta/T^2 = {1e3 * 6 * Jz * dth / T ** 2:.3f}),"
              f" energy {E:.4e} N^2 m^2 s (12 J^2 dtheta^2/T^3 = {12 * Jz ** 2 * dth ** 2 / T ** 3:.4e}),"
              f" bang-bang energy {16 * Jz ** 2 * dth ** 2 / T ** 3:.4e}")
    T_min = np.sqrt(6 * Jz * dth / tau_max)
    print(f"  shortest minimum-energy slew within the wheel limit: T = sqrt(6 J dtheta / tau_max) = {T_min:.2f} s")

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    labels, mins = [], []
    for config, failed in (("orthogonal", ()), ("orthogonal", (2,)), ("pyramid", ()), ("pyramid", (0,)), ("pyramid", (0, 2))):
        Am, Bm = satellite_model(J, wheel_axes(config, failed))
        labels.append(f"{config[:4]} -{list(failed)}" if failed else f"{config[:4]} all")
        mins.append(max(np.linalg.eigvalsh(ctrb_gramian_finite(Am, Bm, 20.0)).min(), 1e-20))
    axes[0].bar(range(len(mins)), mins)
    axes[0].set_yscale("log"); axes[0].set_xticks(range(len(mins))); axes[0].set_xticklabels(labels, fontsize=7, rotation=20)
    axes[0].set_ylabel("lambda_min of W(20 s)"); axes[0].set_title("how controllable: smallest Gramian eigenvalue", fontsize=9)
    tt = np.linspace(0.0, 20.0, 400)
    u, _ = min_energy_input(A1, B1, np.zeros(2), np.array([dth, 0.0]), 20.0)
    axes[1].plot(tt, 1e3 * u(tt)[:, 0], label="minimum energy")
    bang = np.where(tt < 10.0, 1.0, -1.0) * 4 * Jz * dth / 20.0 ** 2
    axes[1].plot(tt, 1e3 * bang, "--", label="bang-bang, same slew")
    axes[1].axhline(1e3 * tau_max, color="k", lw=0.6); axes[1].axhline(-1e3 * tau_max, color="k", lw=0.6)
    axes[1].set_xlabel("t [s]"); axes[1].set_ylabel("wheel torque [mN m]"); axes[1].legend(fontsize=8); axes[1].grid(lw=0.4)
    axes[1].set_title("10 deg slew in 20 s", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / "controllability.png", dpi=150)
    plt.close(fig)
    print(f"saved {out_dir / 'controllability.png'}")


if __name__ == "__main__":
    main()
