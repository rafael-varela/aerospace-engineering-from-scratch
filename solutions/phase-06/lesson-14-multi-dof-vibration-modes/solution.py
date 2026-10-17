"""P06.14 Multi-DOF Vibration & Modes: reference solution.

Modal analysis of a three-mass spacecraft stack on its launch-vehicle interface (axial direction):

* the chain model: three lumped masses on three axial springs (and dampers), fixed at the interface;
* the generalized eigenproblem K phi = w^2 M phi solved by a Cholesky reduction to a symmetric standard
  eigenproblem, mass-normalised modes, M- and K-orthogonality, the Rayleigh quotient;
* modal participation factors and effective masses for base excitation (they sum to the total mass);
* proportional (Rayleigh) damping C = a M + b K fitted to two modal damping ratios, and the test of
  whether a damping matrix is diagonalised by the modes;
* modal superposition of the free response (closed form), of a forced response (each modal SDOF
  integrated separately) and of the frequency response, each checked against direct integration of the
  coupled equations or against the direct inverse (K - w^2 M + i w C)^-1;
* the base-driven steady response: absolute accelerations, and the interface force through the
  "apparent mass" sum of effective masses times SDOF transmissibilities; modal truncation.

The stack's numbers are illustrative (typical of a small spacecraft), not from a particular vehicle.

Conventions: x_i displacement of mass i relative to the inertial frame, z_i relative to the moving
interface; mass 1 is next to the interface; SI units; w in rad/s.

Run:  python solution.py [out_dir]     (prints the tables, writes modal_analysis.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-14-multi-dof-vibration-modes
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import solve_ivp

G0 = 9.80665


def spacecraft_stack() -> dict[str, np.ndarray]:
    """Illustrative axial model, interface upwards: adapter ring 150 kg, service module 600 kg, payload
    module 350 kg; axial stiffnesses 6e8, 4e8 and 2.5e8 N/m (interface-adapter, adapter-service,
    service-payload)."""
    return dict(m=np.array([150.0, 600.0, 350.0]), k=np.array([6.0e8, 4.0e8, 2.5e8]))


# ---------------------------------------------------------------- matrices

# LEARNER
def chain_matrices(m: np.ndarray, k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Mass and stiffness matrices of a chain fixed at one end: spring k[0] joins the base to mass 0,
    spring k[i] joins mass i-1 to mass i. M = diag(m); K is tridiagonal with
    K[i, i] = k[i] + k[i+1] (k[n] = 0) and K[i, i+1] = K[i+1, i] = -k[i+1]."""
    m, k = np.asarray(m, dtype=float), np.asarray(k, dtype=float)
    n = len(m)
    K = np.zeros((n, n))
    for i in range(n):
        K[i, i] += k[i]
        if i > 0:
            K[i - 1, i - 1] += k[i]
            K[i - 1, i] -= k[i]
            K[i, i - 1] -= k[i]
    return np.diag(m), K


# LEARNER
def modes(K: np.ndarray, M: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Natural frequencies (rad/s, ascending) and mass-normalised mode shapes (columns, Phi^T M Phi = I)
    from K phi = w^2 M phi, by Cholesky reduction: M = L L^T, A = L^-1 K L^-T (symmetric), A u = w^2 u
    with numpy.linalg.eigh, phi = L^-T u. Each shape's sign is fixed so its largest entry is positive."""
    L = np.linalg.cholesky(M)
    Li = np.linalg.inv(L)
    A = Li @ K @ Li.T
    lam, U = np.linalg.eigh(0.5 * (A + A.T))
    Phi = Li.T @ U
    for j in range(Phi.shape[1]):
        if Phi[np.argmax(np.abs(Phi[:, j])), j] < 0:
            Phi[:, j] = -Phi[:, j]
    return np.sqrt(lam), Phi


# LEARNER
def rayleigh_quotient(K: np.ndarray, M: np.ndarray, x: np.ndarray) -> float:
    """R(x) = x^T K x / x^T M x: equal to w_i^2 when x is mode i, and never below w_1^2."""
    x = np.asarray(x, dtype=float)
    return float(x @ K @ x / (x @ M @ x))


# LEARNER
def participation(Phi: np.ndarray, M: np.ndarray, r: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Modal participation factors Gamma = Phi^T M r for a rigid-body base motion r (default all ones:
    every mass moves with the interface) and effective masses Gamma^2 (kg, for mass-normalised modes).
    The effective masses sum to r^T M r, the total mass."""
    n = M.shape[0]
    r = np.ones(n) if r is None else np.asarray(r, dtype=float)
    G = Phi.T @ M @ r
    return G, G ** 2


# ---------------------------------------------------------------- damping

# LEARNER
def rayleigh_coefficients(w1: float, w2: float, z1: float, z2: float) -> tuple[float, float]:
    """Coefficients a, b of C = a M + b K that give damping ratios z1 at w1 and z2 at w2, from
    zeta(w) = a/(2 w) + b w/2 (two linear equations)."""
    A = np.array([[1 / (2 * w1), w1 / 2], [1 / (2 * w2), w2 / 2]])
    a, b = np.linalg.solve(A, [z1, z2])
    return float(a), float(b)


# LEARNER
def modal_damping(Phi: np.ndarray, C: np.ndarray, w: np.ndarray) -> tuple[np.ndarray, float]:
    """Modal damping ratios zeta_i = (Phi^T C Phi)_ii / (2 w_i) and a coupling measure, the largest
    off-diagonal |(Phi^T C Phi)_ij| / sqrt(c_ii c_jj); zero when C is proportional (diagonalised by the modes)."""
    Cm = Phi.T @ C @ Phi
    d = np.diag(Cm)
    off = Cm - np.diag(d)
    scale = np.sqrt(np.outer(np.abs(d), np.abs(d)))
    coupling = float(np.max(np.abs(off) / np.where(scale > 0, scale, 1.0))) if len(d) > 1 else 0.0
    return d / (2 * w), coupling


# ---------------------------------------------------------------- superposition and the direct route

# LEARNER
def modal_free_response(t: np.ndarray, x0: np.ndarray, v0: np.ndarray, w: np.ndarray, Phi: np.ndarray,
                        M: np.ndarray, zeta: np.ndarray) -> np.ndarray:
    """Free response by modal superposition, shape (n_dof, len(t)): modal initial conditions
    eta0 = Phi^T M x0, etadot0 = Phi^T M v0; each eta_i is an underdamped SDOF free response at w_i, zeta_i;
    x = Phi eta."""
    t = np.asarray(t, dtype=float)
    e0, ed0 = Phi.T @ M @ x0, Phi.T @ M @ v0
    eta = np.zeros((len(w), len(t)))
    for i in range(len(w)):
        s = zeta[i] * w[i]
        wd = w[i] * np.sqrt(1 - zeta[i] ** 2)
        eta[i] = np.exp(-s * t) * (e0[i] * np.cos(wd * t) + (ed0[i] + s * e0[i]) / wd * np.sin(wd * t))
    return Phi @ eta


# LEARNER
def direct_integration(M: np.ndarray, C: np.ndarray, K: np.ndarray, force: Callable[[float], np.ndarray],
                       x0: np.ndarray, v0: np.ndarray, t: np.ndarray, rtol: float = 1e-10) -> np.ndarray:
    """Integrate M x'' + C x' + K x = f(t) as a first-order system with solve_ivp (DOP853,
    max_step a tenth of the shortest period). Returns x at t, shape (n_dof, len(t))."""
    n = M.shape[0]
    Minv = np.linalg.inv(M)
    wmax = np.sqrt(np.max(np.linalg.eigvals(Minv @ K).real))
    scale = max(np.max(np.abs(x0)), 1e-12)

    def f(tt, u):
        x, v = u[:n], u[n:]
        return np.concatenate([v, Minv @ (force(tt) - C @ v - K @ x)])

    sol = solve_ivp(f, (t[0], t[-1]), np.concatenate([x0, v0]), method="DOP853", t_eval=t, rtol=rtol,
                    atol=1e-12 * scale, max_step=0.1 * 2 * np.pi / wmax)
    return sol.y[:n]


# LEARNER
def modal_integration(w: np.ndarray, Phi: np.ndarray, zeta: np.ndarray, force: Callable[[float], np.ndarray],
                      t: np.ndarray, n_modes: int | None = None, rtol: float = 1e-10) -> np.ndarray:
    """Forced response from rest by modal superposition: each modal coordinate obeys
    eta_i'' + 2 zeta_i w_i eta_i' + w_i^2 eta_i = phi_i^T f(t); integrate the first n_modes of them
    separately with solve_ivp and superpose, x = sum phi_i eta_i."""
    n_modes = len(w) if n_modes is None else n_modes
    x = np.zeros((Phi.shape[0], len(t)))
    for i in range(n_modes):
        ph = Phi[:, i]
        g = lambda tt, u, i=i, ph=ph: [u[1], ph @ force(tt) - 2 * zeta[i] * w[i] * u[1] - w[i] ** 2 * u[0]]
        sol = solve_ivp(g, (t[0], t[-1]), [0.0, 0.0], method="DOP853", t_eval=t, rtol=rtol, atol=1e-16,
                        max_step=0.1 * 2 * np.pi / w[i])
        x += np.outer(ph, sol.y[0])
    return x


# LEARNER
def frf_direct(omega: float, M: np.ndarray, C: np.ndarray, K: np.ndarray) -> np.ndarray:
    """Receptance matrix H(w) = (K - w^2 M + i w C)^-1 (m/N) by a direct solve."""
    return np.linalg.inv(K - omega ** 2 * M + 1j * omega * C)


# LEARNER
def frf_modal(omega: float, w: np.ndarray, Phi: np.ndarray, zeta: np.ndarray, n_modes: int | None = None) -> np.ndarray:
    """Receptance matrix by modal superposition, sum over modes of phi_i phi_i^T / (w_i^2 - w^2 + 2 i zeta_i w_i w);
    exact for proportional damping when all modes are kept."""
    n_modes = len(w) if n_modes is None else n_modes
    H = np.zeros((Phi.shape[0], Phi.shape[0]), dtype=complex)
    for i in range(n_modes):
        H += np.outer(Phi[:, i], Phi[:, i]) / (w[i] ** 2 - omega ** 2 + 2j * zeta[i] * w[i] * omega)
    return H


# LEARNER
def apparent_mass(omega: np.ndarray, w: np.ndarray, zeta: np.ndarray, m_eff: np.ndarray,
                  n_modes: int | None = None) -> np.ndarray:
    """Interface force per unit interface acceleration (kg, complex) for a base-driven structure with
    proportional damping: sum over modes of m_eff,i (1 + 2 i zeta_i r_i) / (1 - r_i^2 + 2 i zeta_i r_i),
    r_i = w / w_i. Equal to the total mass as w -> 0. If modes are truncated, the residual effective mass
    of the omitted modes is added as rigid mass (it moves with the interface)."""
    omega = np.atleast_1d(np.asarray(omega, dtype=float))
    n_modes = len(w) if n_modes is None else n_modes
    out = np.zeros(omega.shape, dtype=complex)
    for i in range(n_modes):
        r = omega / w[i]
        out += m_eff[i] * (1 + 2j * zeta[i] * r) / (1 - r ** 2 + 2j * zeta[i] * r)
    out += np.sum(m_eff[n_modes:])
    return out


# LEARNER
def base_sine_accelerations(omega: float, A: float, M: np.ndarray, C: np.ndarray, K: np.ndarray) -> tuple[np.ndarray, complex]:
    """Steady absolute acceleration amplitudes (complex, m/s^2) of every mass and the interface force
    (complex, N) for an interface acceleration A e^{i w t}, by a direct solve of the relative motion
    (K - w^2 M + i w C) Z = -M r A (r = ones). The absolute accelerations are A - w^2 Z, and the force
    the interface must supply is the total inertia force, sum_i m_i x_i''."""
    n = M.shape[0]
    r = np.ones(n)
    Z = np.linalg.solve(K - omega ** 2 * M + 1j * omega * C, -M @ r * A)
    acc = A - omega ** 2 * Z
    return acc, complex(np.sum(M @ acc))


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save modal_analysis.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    s = spacecraft_stack()
    M, K = chain_matrices(s["m"], s["k"])
    w, Phi = modes(K, M)
    G, meff = participation(Phi, M)
    print(f"stack: masses {s['m']} kg, springs {s['k']} N/m, total {s['m'].sum():.0f} kg")
    print(" mode   f (Hz)    shape (adapter, service, payload), mass-normalised      Gamma    m_eff (kg)   share")
    for i in range(3):
        print(f"  {i + 1}   {w[i] / 2 / np.pi:8.3f}   [{Phi[0, i]:+.5f} {Phi[1, i]:+.5f} {Phi[2, i]:+.5f}]   "
              f"{G[i]:9.4f}   {meff[i]:9.2f}   {meff[i] / s['m'].sum() * 100:5.1f} %")
    print(f" sum of effective masses {meff.sum():.6f} kg; max|Phi^T M Phi - I| {np.max(np.abs(Phi.T @ M @ Phi - np.eye(3))):.1e}; "
          f"max off-diagonal |Phi^T K Phi| {np.max(np.abs(Phi.T @ K @ Phi - np.diag(w ** 2))):.1e}")
    print(f" Rayleigh quotient of the static 1 g deflection: f = {np.sqrt(rayleigh_quotient(K, M, np.linalg.solve(K, s['m']))) / 2 / np.pi:.3f} Hz "
          f"(f1 = {w[0] / 2 / np.pi:.3f} Hz)")

    a, b = rayleigh_coefficients(w[0], w[1], 0.02, 0.02)
    C = a * M + b * K
    zeta, coup = modal_damping(Phi, C, w)
    print(f"\nRayleigh damping fitted to 2 % in modes 1 and 2: a = {a:.4f} 1/s, b = {b:.4e} s; "
          f"modal zeta = {np.round(zeta, 5)}; coupling {coup:.1e}")

    # free response: modal vs direct
    t = np.linspace(0.0, 0.2, 4001)
    x0 = np.linalg.solve(K, s["m"] * 5 * G0)            # released from a 5 g static deflection (thrust cut-off)
    xm = modal_free_response(t, x0, np.zeros(3), w, Phi, M, zeta)
    xd = direct_integration(M, C, K, lambda tt: np.zeros(3), x0, np.zeros(3), t)
    print(f"free response from a 5 g static deflection (max {np.max(np.abs(x0)) * 1e3:.3f} mm): modal vs direct max diff "
          f"{np.max(np.abs(xm - xd)):.1e} m")

    # forced: a 1 g-equivalent sine force on the payload at mode-1 frequency, from rest
    F0 = s["m"][2] * G0
    fvec = lambda tt: np.array([0.0, 0.0, F0 * np.sin(w[0] * tt)])
    xf_d = direct_integration(M, C, K, fvec, np.zeros(3), np.zeros(3), t)
    for nm in (1, 2, 3):
        xf_m = modal_integration(w, Phi, zeta, fvec, t, n_modes=nm)
        print(f"   sine force on the payload at f1, {nm} mode(s): max |modal - direct| {np.max(np.abs(xf_m - xf_d)):.2e} m "
              f"(response peak {np.max(np.abs(xf_d)) * 1e3:.4f} mm)")

    # base sine at mode 1: accelerations and interface force
    A = 1.0 * G0
    acc, Fi = base_sine_accelerations(w[0], A, M, C, K)
    print(f"\n1 g interface sine at f1 = {w[0] / 2 / np.pi:.2f} Hz: |acceleration| (g) adapter {abs(acc[0]) / G0:.2f}, "
          f"service {abs(acc[1]) / G0:.2f}, payload {abs(acc[2]) / G0:.2f}; interface force {abs(Fi) / 1e3:.3f} kN "
          f"(rigid-body {s['m'].sum() * G0 / 1e3:.2f} kN)")
    for nm in (1, 2, 3):
        Fm = apparent_mass(w[0], w, zeta, meff, n_modes=nm)[0] * A
        print(f"   apparent mass with {nm} mode(s) + residual: interface force {abs(Fm) / 1e3:.3f} kN")

    # figure
    f = np.linspace(5, 400, 3000)
    om = 2 * np.pi * f
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    for i in range(3):
        ax[0].plot(np.r_[0, Phi[:, i]], np.arange(4), "o-", label=f"mode {i + 1}: {w[i] / 2 / np.pi:.1f} Hz")
    ax[0].set_yticks(range(4), ["interface", "adapter", "service", "payload"])
    ax[0].set_xlabel("mass-normalised shape")
    ax[0].set_title("axial modes")
    ax[0].legend(fontsize=8)
    Hd = np.array([abs(frf_direct(o, M, C, K)[2, 2]) for o in om])
    Hm = np.array([abs(frf_modal(o, w, Phi, zeta)[2, 2]) for o in om])
    H1 = np.array([abs(frf_modal(o, w, Phi, zeta, 1)[2, 2]) for o in om])
    ax[1].semilogy(f, Hd, label="direct inverse")
    ax[1].semilogy(f, Hm, "--", label="modal, 3 modes")
    ax[1].semilogy(f, H1, ":", label="modal, 1 mode")
    ax[1].set_xlabel("f (Hz)")
    ax[1].set_ylabel("|H_33| (m/N)")
    ax[1].set_title("payload driving-point receptance")
    ax[1].legend(fontsize=8)
    ax[2].semilogy(f, np.abs(apparent_mass(om, w, zeta, meff)), label="3 modes")
    ax[2].semilogy(f, np.abs(apparent_mass(om, w, zeta, meff, 1)), ":", label="1 mode + residual mass")
    ax[2].axhline(s["m"].sum(), color="0.6", lw=0.8)
    ax[2].set_xlabel("f (Hz)")
    ax[2].set_ylabel("|F / a| (kg)")
    ax[2].set_title("interface apparent mass")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "modal_analysis.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'modal_analysis.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
