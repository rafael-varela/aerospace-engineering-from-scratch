"""P02.04 Systems of ODEs & Phase Portraits: starter code.

Planar linear systems dx/dt = A x classified by the trace and determinant of A (saddle, node, spiral,
center, degenerate cases), solved in closed form through eigenvalues and eigenvectors, integrated by an
RK4 stepper written from scratch, and drawn as phase portraits. The examples are the two longitudinal
modes of an aircraft: the short-period mode (angle of attack and pitch rate) and the phugoid (speed and
pitch angle), plus the saddle that appears when the centre of gravity moves too far aft.

Conventions:
    * x = (x1, x2) is the state; A is 2x2 with real entries; t in seconds.
    * Short-period state (alpha [rad], q [rad/s]); phugoid state (du/V0 [-], dtheta [rad]).
    * ``classify`` names the fixed point at the origin from tau = tr A and Delta = det A.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the mode table, writes phase_portraits.png)
Test: pytest curriculum/phase-02/lesson-04-systems-and-phase-portraits
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

G0 = 9.81  # m/s^2


def short_period_matrix(z_alpha_over_v: float = -0.8, m_alpha: float = -4.0, m_q: float = -1.2) -> np.ndarray:
    """Return A for the short-period approximation d/dt (alpha, q) = A (alpha, q):
    A = [[Z_alpha / V0, 1], [M_alpha, M_q]] with the defaults of a jet transport in cruise (units 1/s, 1, 1/s^2, 1/s).

    Z_alpha / V0 < 0: lift from angle of attack rotates the velocity vector; M_alpha < 0: static stability
    (a nose-up disturbance makes a nose-down moment); M_q < 0: pitch damping. The 1 couples q into alpha.
    """
    return np.array([[z_alpha_over_v, 1.0], [m_alpha, m_q]])


def phugoid_matrix(v0: float = 230.0, lift_to_drag: float = 17.0, g: float = G0) -> np.ndarray:
    """Return A for Lanchester's phugoid approximation d/dt (du/V0, dtheta) = A (du/V0, dtheta):
    A = [[-2 g / (V0 L/D), -g / V0], [2 g / V0, 0]].

    Speed above trim raises lift (2 g / V0 term) and pitches the path up; a pitched-up path loses speed to
    gravity (-g / V0); drag at the higher speed removes energy (-2 g / (V0 L/D)). Natural frequency
    sqrt(2) g / V0, damping ratio 1 / (sqrt(2) L/D), independent of everything but speed and L/D.
    """
    if v0 <= 0.0 or lift_to_drag <= 0.0:
        raise ValueError("speed and L/D must be positive")
    return np.array([[-2.0 * g / (v0 * lift_to_drag), -g / v0], [2.0 * g / v0, 0.0]])


def second_order_matrix(m: float, c: float, k: float) -> np.ndarray:
    """Return A = [[0, 1], [-k/m, -c/m]], the first-order form of m x'' + c x' + k x = 0 (P02.02)."""
    return np.array([[0.0, 1.0], [-k / m, -c / m]])


# ---------------------------------------------------------------- classification

def trace_determinant(A: np.ndarray) -> tuple[float, float]:
    """Return (tau, Delta) = (tr A, det A). The characteristic polynomial is lambda^2 - tau lambda + Delta."""
    raise NotImplementedError


def eigen_structure(A: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (eigenvalues, eigenvectors as columns) of A from the trace and determinant, by hand, not ``numpy.linalg.eig``.

    lambda = tau/2 +- sqrt(tau^2/4 - Delta) (complex when the discriminant is negative). For each eigenvalue
    the eigenvector is a non-zero column of adj(A - lambda I): (b, lambda - a) or (lambda - d, c), whichever is
    larger, normalised to unit length. For a repeated eigenvalue with A != lambda I only one eigenvector exists
    and it is returned twice.
    """
    raise NotImplementedError


def classify(A: np.ndarray, tol: float = 1e-12) -> str:
    """Return the type of the fixed point at the origin from (tau, Delta):

    Delta < 0: 'saddle'; Delta = 0: 'non-isolated' (a line of fixed points); tau^2 > 4 Delta > 0: 'stable node' /
    'unstable node'; tau^2 < 4 Delta: 'stable spiral' / 'unstable spiral' / 'center' (tau = 0);
    tau^2 = 4 Delta: 'stable star' / 'unstable star' if A is a multiple of I, else 'stable degenerate node' / 'unstable degenerate node'.
    """
    raise NotImplementedError


def modal_parameters(A: np.ndarray) -> dict[str, float]:
    """Return the engineer's description of an oscillatory mode: omega_n = sqrt(Delta), zeta = -tau / (2 sqrt(Delta)),
    damped frequency, period, time to half (or double) amplitude ln 2 / |Re lambda|. Raise ``ValueError`` unless the eigenvalues are complex."""
    raise NotImplementedError


# ---------------------------------------------------------------- solutions

def solve_linear(A: np.ndarray, x0: np.ndarray, t: np.ndarray | float) -> np.ndarray:
    """Return x(t) solving dx/dt = A x, x(0) = x0, in closed form; shape (len(t), 2).

    Distinct eigenvalues: x = V diag(e^{lambda t}) V^{-1} x0 (complex arithmetic, real part taken; for a complex
    pair this is e^{sigma t}(cos, sin) along the real and imaginary parts of the eigenvector). Repeated eigenvalue:
    x = e^{lambda t} (I + (A - lambda I) t) x0, which covers both the star and the degenerate node, because
    (A - lambda I)^2 = 0 then (Cayley-Hamilton).
    """
    raise NotImplementedError


def rk4_trajectory(A: np.ndarray, x0: np.ndarray, t_end: float, n_steps: int) -> tuple[np.ndarray, np.ndarray]:
    """Integrate dx/dt = A x with classical RK4 on a uniform grid; return (times (n+1,), states (n+1, 2))."""
    raise NotImplementedError


def vector_field(A: np.ndarray, X: np.ndarray, Y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (U, V) = (dx1/dt, dx2/dt) on the mesh (X, Y): U = a X + b Y, V = c X + d Y."""
    raise NotImplementedError


def nullcline_slopes(A: np.ndarray) -> tuple[float, float]:
    """Return (dx2/dx1 along the x1-nullcline, along the x2-nullcline): the lines a x1 + b x2 = 0 and c x1 + d x2 = 0.

    Slopes are -a/b and -c/d; ``inf`` when the line is vertical. Trajectories cross the x1-nullcline vertically
    and the x2-nullcline horizontally; the two lines meet only at the origin when Delta != 0.
    """
    raise NotImplementedError


def center_invariant(A: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Return Q = c x1^2 - 2 a x1 x2 - b x2^2 along states ``x`` (shape (n, 2)) for A = [[a, b], [c, d]].

    When tr A = 0 (a center), dQ/dt = -2 b (a + d) x2^2 = 0: Q is conserved and the orbits are the ellipses
    Q = const. For the undamped oscillator A = [[0, 1], [-k/m, 0]] this is -(k/m) x^2 - v^2, minus 2/m times the energy.
    """
    raise NotImplementedError


def phase_portrait(A: np.ndarray, starts: np.ndarray, t_end: float, n_steps: int = 400) -> list[np.ndarray]:
    """Return a list of closed-form trajectories (each (n_steps + 1, 2)) from each row of ``starts`` over [0, t_end]."""
    raise NotImplementedError


# ---------------------------------------------------------------- the example

def main(out_dir: str | Path = ".") -> None:
    """Print the classification and modal parameters of the aircraft modes and save phase_portraits.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    cases = {
        "short period (jet transport, cruise)": short_period_matrix(),
        "phugoid (V0 = 230 m/s, L/D = 17)": phugoid_matrix(),
        "short period, CG too far aft (M_alpha = +2)": short_period_matrix(m_alpha=2.0),
        "short period with tau forced to 0 (M_q = +0.8): center": np.array([[-0.8, 1.0], [-4.0, 0.8]]),
        "undamped oscillator m x'' + k x = 0 (center)": second_order_matrix(1.0, 0.0, 4.0),
        "critically damped oscillator (degenerate node)": second_order_matrix(1.0, 4.0, 4.0),
        "overdamped oscillator (stable node)": second_order_matrix(1.0, 5.0, 4.0),
    }
    print(f"{'case':<52} {'tau':>7} {'Delta':>7}  {'type':<24} eigenvalues")
    for name, A in cases.items():
        tau, delta = trace_determinant(A)
        lam, _ = eigen_structure(A)
        print(f"{name:<52} {tau:7.3f} {delta:7.3f}  {classify(A):<24} {np.array2string(lam, precision=4)}")
    for name in ("short period (jet transport, cruise)", "phugoid (V0 = 230 m/s, L/D = 17)"):
        mp = modal_parameters(cases[name])
        print(f"{name}: omega_n = {mp['omega_n']:.4f} rad/s, zeta = {mp['zeta']:.4f}, period = {mp['period']:.2f} s,"
              f" time to half amplitude = {mp['t_half']:.2f} s")
    A_sp = short_period_matrix()
    lam, V = eigen_structure(A_sp)
    print("short-period eigenvalues by hand:", lam, " numpy:", np.linalg.eigvals(A_sp))
    # Closed form vs RK4 vs solve_ivp on the short period from a 5 degree alpha disturbance.
    x0 = np.array([np.radians(5.0), 0.0])
    t = np.linspace(0.0, 8.0, 801)
    x_cf = solve_linear(A_sp, x0, t)
    _, x_rk = rk4_trajectory(A_sp, x0, 8.0, 800)
    sol = solve_ivp(lambda tt, y: A_sp @ y, (0.0, 8.0), x0, rtol=1e-10, atol=1e-12, t_eval=t)
    print(f"short period from alpha = 5 deg: max |RK4 - closed form| = {np.max(np.abs(x_rk - x_cf)):.2e},"
          f" max |solve_ivp - closed form| = {np.max(np.abs(sol.y.T - x_cf)):.2e}")
    print(f"peak pitch rate {np.degrees(np.min(x_cf[:, 1])):.2f} deg/s at t = {t[np.argmin(x_cf[:, 1])]:.2f} s;"
          f" alpha within 0.5 deg after t = {t[np.max(np.where(np.abs(np.degrees(x_cf[:, 0])) > 0.5)[0]) + 1]:.2f} s")
    A_sad = short_period_matrix(m_alpha=2.0)
    lam_s, V_s = eigen_structure(A_sad)
    print(f"aft-CG saddle: eigenvalues {lam_s}, unstable direction {V_s[:, 1]} (alpha, q), doubling time {np.log(2) / lam_s[1]:.2f} s")

    # Figure: four portraits.
    fig, axes = plt.subplots(2, 2, figsize=(10, 9))
    panels = [
        (axes[0, 0], "short period: stable spiral", A_sp, 0.3, 8.0),
        (axes[0, 1], "aft CG: saddle", A_sad, 0.3, 2.5),
        (axes[1, 0], "phugoid: slow stable spiral", phugoid_matrix(), 0.3, 300.0),
        (axes[1, 1], "no damping: center", np.array([[-0.8, 1.0], [-4.0, 0.8]]), 0.3, 6.0),
    ]
    for ax, title, A, span, t_end in panels:
        g = np.linspace(-span, span, 17)
        X, Y = np.meshgrid(g, g)
        U, Vf = vector_field(A, X, Y)
        ax.streamplot(X, Y, U, Vf, color="0.75", density=0.8, linewidth=0.6, arrowsize=0.7)
        ang = np.linspace(0.0, 2 * np.pi, 9)[:-1]
        starts = 0.9 * span * np.column_stack([np.cos(ang), np.sin(ang)])
        for traj in phase_portrait(A, starts, t_end):
            ax.plot(traj[:, 0], traj[:, 1], lw=1.0)
        lam, V = eigen_structure(A)
        if np.isrealobj(lam):
            for k in range(2):
                ax.plot([-span * V[0, k], span * V[0, k]], [-span * V[1, k], span * V[1, k]], "k--", lw=0.8)
        ax.set_xlim(-span, span); ax.set_ylim(-span, span)
        ax.set_title(f"{title}: {classify(A)}", fontsize=9)
        ax.set_xlabel("x1"); ax.set_ylabel("x2"); ax.grid(lw=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "phase_portraits.png", dpi=150)
    print(f"saved {out_dir / 'phase_portraits.png'}")


if __name__ == "__main__":
    main()
