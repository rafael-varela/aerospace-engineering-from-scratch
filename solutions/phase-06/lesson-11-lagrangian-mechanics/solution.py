"""P06.11 Lagrangian Mechanics: reference solution.

A double pendulum, derived the analytical way and checked the vectorial way:

* Lagrange's equations d/dt(dL/dq') - dL/dq = 0 for any number of coordinates, built in SymPy;
* the planar double pendulum (two point masses on massless rigid links, a fixed pivot): kinetic and
  potential energy from the Cartesian positions, the equations of motion solved for the two angular
  accelerations, and the same equations written as M(q) q'' = f(q, q');
* an independent Newton-Euler model: free-body diagrams of the two masses with the two link tensions
  as unknowns, solved as a 4 x 4 linear system for (theta1'', theta2'', T1, T2);
* integration with `scipy.integrate.solve_ivp`, conservation of energy, the small-angle normal modes
  (the generalized eigenproblem K phi = w^2 M phi and its closed form), and the sensitivity to initial
  conditions (a finite-time Lyapunov exponent by Benettin renormalisation) at large amplitude.

The aerospace example is a helicopter slung load in hover: a hook block on a long cable with a load on
a short sling below it, the hover point treated as fixed. Numbers are illustrative, not from a
particular aircraft.

Conventions: angles theta1, theta2 from the downward vertical, positive counter-clockwise; the pivot
at the origin, y up; state y = [theta1, theta2, omega1, omega2]; SI units, radians.

Run:  python solution.py [out_dir]     (prints the tables, writes double_pendulum.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-11-lagrangian-mechanics
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

G0 = 9.80665


@dataclass(frozen=True)
class DoublePendulum:
    """Two point masses m1, m2 on massless rigid links of lengths l1 (pivot to m1) and l2 (m1 to m2)."""
    m1: float
    m2: float
    l1: float
    l2: float
    g: float = G0


def slung_load() -> DoublePendulum:
    """Illustrative helicopter slung load in hover: an 80 kg hook block on a 12 m cable carrying a
    1500 kg container whose centre of mass hangs 4 m below the hook on its sling."""
    return DoublePendulum(m1=80.0, m2=1500.0, l1=12.0, l2=4.0)


def lab_pendulum() -> DoublePendulum:
    """Equal masses and links (1 kg, 1 m): the textbook double pendulum used for the chaos study."""
    return DoublePendulum(m1=1.0, m2=1.0, l1=1.0, l2=1.0)


# ---------------------------------------------------------------- Lagrange, symbolically

# LEARNER
def euler_lagrange(L: sp.Expr, qs: list[sp.Function], t: sp.Symbol) -> list[sp.Expr]:
    """Lagrange's equations for the coordinates ``qs`` (SymPy functions of ``t``): returns the list of
    expressions d/dt(dL/dq_i') - dL/dq_i, each of which is zero along a motion."""
    out = []
    for q in qs:
        qd = sp.diff(q, t)
        out.append(sp.diff(sp.diff(L, qd), t) - sp.diff(L, q))
    return out


# LEARNER
def double_pendulum_lagrangian() -> tuple[sp.Expr, list[sp.Function], sp.Symbol, tuple[sp.Symbol, ...]]:
    """The Lagrangian T - V of the planar double pendulum, built from the Cartesian positions
    x1 = l1 sin th1, y1 = -l1 cos th1, x2 = x1 + l2 sin th2, y2 = y1 - l2 cos th2.
    Returns (L, [th1(t), th2(t)], t, (m1, m2, l1, l2, g))."""
    t = sp.Symbol("t")
    m1, m2, l1, l2, g = sp.symbols("m1 m2 l1 l2 g", positive=True)
    th1, th2 = sp.Function("th1")(t), sp.Function("th2")(t)
    x1, y1 = l1 * sp.sin(th1), -l1 * sp.cos(th1)
    x2, y2 = x1 + l2 * sp.sin(th2), y1 - l2 * sp.cos(th2)
    T = sp.Rational(1, 2) * m1 * (sp.diff(x1, t) ** 2 + sp.diff(y1, t) ** 2) \
        + sp.Rational(1, 2) * m2 * (sp.diff(x2, t) ** 2 + sp.diff(y2, t) ** 2)
    V = m1 * g * y1 + m2 * g * y2
    return T - V, [th1, th2], t, (m1, m2, l1, l2, g)


# LEARNER
@lru_cache(maxsize=None)
def derive_equations() -> tuple[sp.Matrix, sp.Matrix, tuple[sp.Symbol, ...]]:
    """Lagrange's equations of the double pendulum written as M(q) q'' = f(q, q').

    Build the Lagrangian, apply ``euler_lagrange``, replace th_i'' by symbols d1, d2, th_i' by w1, w2
    and th_i by a1, a2, expand, split the two equations (which are linear in d1, d2) with
    ``sympy.linear_eq_to_matrix``, and tidy each entry with ``trigsimp``.
    Returns (M, f, (a1, a2, w1, w2, m1, m2, l1, l2, g)). Cached: SymPy takes about a second."""
    L, (th1, th2), t, (m1, m2, l1, l2, g) = double_pendulum_lagrangian()
    eqs = euler_lagrange(L, [th1, th2], t)
    a1, a2, w1, w2, d1, d2 = sp.symbols("a1 a2 w1 w2 d1 d2")
    eqs = [e.subs({sp.diff(th1, t, 2): d1, sp.diff(th2, t, 2): d2}) for e in eqs]
    eqs = [e.subs({sp.diff(th1, t): w1, sp.diff(th2, t): w2}) for e in eqs]
    eqs = [sp.expand(e.subs({th1: a1, th2: a2})) for e in eqs]
    M, f = sp.linear_eq_to_matrix(eqs, [d1, d2])
    tidy = lambda e: sp.factor(sp.trigsimp(e))           # entry by entry: fast, and readable output
    return M.applyfunc(tidy), f.applyfunc(tidy), (a1, a2, w1, w2, m1, m2, l1, l2, g)


# LEARNER
@lru_cache(maxsize=None)
def derive_accelerations() -> tuple[tuple[sp.Expr, sp.Expr], tuple[sp.Symbol, ...]]:
    """Solve M(q) q'' = f from ``derive_equations`` for the two angular accelerations.
    Returns ((th1_dd, th2_dd), (a1, a2, w1, w2, m1, m2, l1, l2, g)) in plain symbols."""
    M, f, syms = derive_equations()
    sol = M.LUsolve(f)
    return (sol[0], sol[1]), syms


@lru_cache(maxsize=None)
def _lagrange_function() -> Callable:
    (e1, e2), syms = derive_accelerations()
    return sp.lambdify(syms, [e1, e2], "numpy")


def lagrange_accelerations(y: np.ndarray, p: DoublePendulum) -> np.ndarray:
    """Angular accelerations [th1'', th2''] from the SymPy derivation, evaluated numerically at the
    state y = [th1, th2, w1, w2] (one state or an array of shape (4, n))."""
    f = _lagrange_function()
    return np.array(f(y[0], y[1], y[2], y[3], p.m1, p.m2, p.l1, p.l2, p.g), dtype=float)


# ---------------------------------------------------------------- the same equations by hand

# LEARNER
def mass_matrix(theta: np.ndarray, p: DoublePendulum) -> np.ndarray:
    """Configuration-dependent mass matrix M(q) of the double pendulum, with D = th1 - th2:
    [[(m1 + m2) l1^2, m2 l1 l2 cos D], [m2 l1 l2 cos D, m2 l2^2]]."""
    c = np.cos(theta[0] - theta[1])
    return np.array([[(p.m1 + p.m2) * p.l1 ** 2, p.m2 * p.l1 * p.l2 * c],
                     [p.m2 * p.l1 * p.l2 * c, p.m2 * p.l2 ** 2]])


# LEARNER
def generalized_forces(y: np.ndarray, p: DoublePendulum) -> np.ndarray:
    """Right-hand side f(q, q') of M(q) q'' = f: the velocity-squared (centripetal) terms and gravity,
    f1 = -m2 l1 l2 sin D w2^2 - (m1 + m2) g l1 sin th1,  f2 = m2 l1 l2 sin D w1^2 - m2 g l2 sin th2."""
    s = np.sin(y[0] - y[1])
    return np.array([-p.m2 * p.l1 * p.l2 * s * y[3] ** 2 - (p.m1 + p.m2) * p.g * p.l1 * np.sin(y[0]),
                     p.m2 * p.l1 * p.l2 * s * y[2] ** 2 - p.m2 * p.g * p.l2 * np.sin(y[1])])


# LEARNER
def rhs(t: float, y: np.ndarray, p: DoublePendulum) -> np.ndarray:
    """State derivative [w1, w2, th1'', th2''] from M(q) q'' = f(q, q') (2 x 2 solve by Cramer's rule)."""
    M = mass_matrix(y[:2], p)
    f = generalized_forces(y, p)
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    d1 = (f[0] * M[1, 1] - M[0, 1] * f[1]) / det
    d2 = (M[0, 0] * f[1] - M[1, 0] * f[0]) / det
    return np.array([y[2], y[3], d1, d2])


# ---------------------------------------------------------------- independent route: Newton-Euler

# LEARNER
def newton_euler(y: np.ndarray, p: DoublePendulum) -> tuple[np.ndarray, np.ndarray]:
    """Accelerations and link tensions from free-body diagrams, with no energy or Lagrangian used.

    Unit vectors along each link (pointing from its upper end to its lower end) e_i = (sin th_i, -cos th_i)
    and perpendicular to it n_i = (cos th_i, sin th_i). Kinematics of rigid links:
    a1 = l1 (th1'' n1 - w1^2 e1),  a2 = a1 + l2 (th2'' n2 - w2^2 e2).
    Kinetics, with link tensions T1 (pivot-m1) and T2 (m1-m2) pulling each mass towards the other end:
    m2 a2 = -T2 e2 + m2 g_vec,  m1 a1 = -T1 e1 + T2 e2 + m1 g_vec,  g_vec = (0, -g).
    Four scalar equations in (th1'', th2'', T1, T2); returns ([th1'', th2''], [T1, T2])."""
    th1, th2, w1, w2 = y
    e1, n1 = np.array([np.sin(th1), -np.cos(th1)]), np.array([np.cos(th1), np.sin(th1)])
    e2, n2 = np.array([np.sin(th2), -np.cos(th2)]), np.array([np.cos(th2), np.sin(th2)])
    gv = np.array([0.0, -p.g])
    # unknowns u = [d1, d2, T1, T2]; a1 = l1 d1 n1 - l1 w1^2 e1 ; a2 = a1 + l2 d2 n2 - l2 w2^2 e2
    A = np.zeros((4, 4))
    b = np.zeros(4)
    # mass 2: m2 (l1 d1 n1 + l2 d2 n2) + T2 e2 = m2 g + m2 (l1 w1^2 e1 + l2 w2^2 e2)
    A[0:2, 0] = p.m2 * p.l1 * n1
    A[0:2, 1] = p.m2 * p.l2 * n2
    A[0:2, 3] = e2
    b[0:2] = p.m2 * gv + p.m2 * (p.l1 * w1 ** 2 * e1 + p.l2 * w2 ** 2 * e2)
    # mass 1: m1 l1 d1 n1 + T1 e1 - T2 e2 = m1 g + m1 l1 w1^2 e1
    A[2:4, 0] = p.m1 * p.l1 * n1
    A[2:4, 2] = e1
    A[2:4, 3] = -e2
    b[2:4] = p.m1 * gv + p.m1 * p.l1 * w1 ** 2 * e1
    u = np.linalg.solve(A, b)
    return u[:2], u[2:]


# ---------------------------------------------------------------- integration and energy

# LEARNER
def simulate(y0: np.ndarray, t_eval: np.ndarray, p: DoublePendulum,
             rtol: float = 1e-10, atol: float = 1e-10) -> np.ndarray:
    """Integrate the double pendulum with DOP853 and return the states at ``t_eval``, shape (4, n)."""
    sol = solve_ivp(rhs, (t_eval[0], t_eval[-1]), np.asarray(y0, dtype=float), method="DOP853",
                    t_eval=t_eval, rtol=rtol, atol=atol, args=(p,))
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y


# LEARNER
def energy(y: np.ndarray, p: DoublePendulum) -> np.ndarray:
    """Total mechanical energy T + V (J), with V = 0 at the pivot height; y has shape (4,) or (4, n)."""
    th1, th2, w1, w2 = y
    T = 0.5 * (p.m1 + p.m2) * p.l1 ** 2 * w1 ** 2 + 0.5 * p.m2 * p.l2 ** 2 * w2 ** 2 \
        + p.m2 * p.l1 * p.l2 * w1 * w2 * np.cos(th1 - th2)
    V = -(p.m1 + p.m2) * p.g * p.l1 * np.cos(th1) - p.m2 * p.g * p.l2 * np.cos(th2)
    return T + V


def positions(y: np.ndarray, p: DoublePendulum) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Cartesian positions (x1, y1, x2, y2) of the two masses (given)."""
    x1, y1 = p.l1 * np.sin(y[0]), -p.l1 * np.cos(y[0])
    return x1, y1, x1 + p.l2 * np.sin(y[1]), y1 - p.l2 * np.cos(y[1])


# ---------------------------------------------------------------- small-angle normal modes

# LEARNER
def linearized_matrices(p: DoublePendulum) -> tuple[np.ndarray, np.ndarray]:
    """Mass and stiffness matrices of the motion linearised about the hanging equilibrium:
    M0 = M(0) and K0 = diag((m1 + m2) g l1, m2 g l2), so that M0 q'' + K0 q = 0."""
    M0 = mass_matrix(np.zeros(2), p)
    K0 = np.diag([(p.m1 + p.m2) * p.g * p.l1, p.m2 * p.g * p.l2])
    return M0, K0


# LEARNER
def normal_modes(p: DoublePendulum) -> tuple[np.ndarray, np.ndarray]:
    """Natural frequencies (rad/s, ascending) and mode shapes (columns, scaled so theta1 = 1) of the
    small-angle motion, from the generalized eigenproblem K0 phi = w^2 M0 phi."""
    M0, K0 = linearized_matrices(p)
    lam, phi = eigh(K0, M0)
    return np.sqrt(lam), phi / phi[0]


# LEARNER
def normal_mode_frequencies_closed_form(p: DoublePendulum) -> np.ndarray:
    """Roots of det(K0 - w^2 M0) = 0, which reduces to
    m1 l1 l2 w^4 - (m1 + m2) g (l1 + l2) w^2 + (m1 + m2) g^2 = 0."""
    a = p.m1 * p.l1 * p.l2
    b = -(p.m1 + p.m2) * p.g * (p.l1 + p.l2)
    c = (p.m1 + p.m2) * p.g ** 2
    disc = np.sqrt(b * b - 4 * a * c)
    w2 = np.array([2 * c / (-b + disc), (-b + disc) / (2 * a)])   # stable form for the small root
    return np.sqrt(w2)


def dominant_frequency(t: np.ndarray, x: np.ndarray) -> float:
    """Angular frequency (rad/s) of the largest FFT peak of a uniformly sampled signal, refined by a
    parabola through the peak bin and its neighbours (given)."""
    x = np.asarray(x, dtype=float) - np.mean(x)
    n = len(x)
    X = np.abs(np.fft.rfft(x * np.hanning(n)))
    k = int(np.argmax(X[1:]) + 1)
    if 1 <= k < len(X) - 1:
        a, b, c = np.log(X[k - 1] + 1e-300), np.log(X[k] + 1e-300), np.log(X[k + 1] + 1e-300)
        k = k + 0.5 * (a - c) / (a - 2 * b + c)
    return 2 * np.pi * k / (n * (t[1] - t[0]))


# ---------------------------------------------------------------- sensitivity to initial conditions

# LEARNER
def lyapunov_estimate(y0: np.ndarray, p: DoublePendulum, t_end: float, d0: float = 1e-8,
                      tau: float = 0.5) -> tuple[float, np.ndarray, np.ndarray]:
    """Finite-time largest Lyapunov exponent (1/s) by Benettin's method: integrate a reference and a
    neighbour d0 away for a time tau, record ln(d / d0), pull the neighbour back to distance d0 along
    the current separation, repeat. The separation is measured in the state space with angles in rad
    and rates scaled by sqrt(l1/g) so both have the same units. Returns (lambda, t_k, running mean)."""
    scale = np.array([1.0, 1.0, np.sqrt(p.l1 / p.g), np.sqrt(p.l1 / p.g)])
    y = np.asarray(y0, dtype=float)
    z = y + np.array([d0, 0.0, 0.0, 0.0])

    def f2(t, s):
        return np.concatenate([rhs(t, s[:4], p), rhs(t, s[4:], p)])

    n = int(round(t_end / tau))
    logs, times, running = 0.0, [], []
    for k in range(n):
        sol = solve_ivp(f2, (0.0, tau), np.concatenate([y, z]), method="DOP853", rtol=1e-10, atol=1e-12)
        y, z = sol.y[:4, -1], sol.y[4:, -1]
        dv = (z - y) * scale
        d = np.linalg.norm(dv)
        logs += np.log(d / d0)
        z = y + dv / d * d0 / scale
        times.append((k + 1) * tau)
        running.append(logs / ((k + 1) * tau))
    return running[-1], np.array(times), np.array(running)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save double_pendulum.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = slung_load()
    Ms, fs, _ = derive_equations()
    print("Lagrange (SymPy): M(q) =", Ms.tolist())
    print("                  f(q, q') =", fs.T.tolist()[0])
    w, phi = normal_modes(p)
    wc = normal_mode_frequencies_closed_form(p)
    print(f"\nslung load: hook {p.m1:.0f} kg on {p.l1:.0f} m, load {p.m2:.0f} kg on {p.l2:.0f} m")
    for i in range(2):
        print(f"  mode {i + 1}: w = {w[i]:.5f} rad/s (closed form {wc[i]:.5f}), period {2 * np.pi / w[i]:.3f} s, "
              f"shape th2/th1 = {phi[1, i]:+.4f}")
    print(f"  single pendulum of length l1 + l2: period {2 * np.pi * np.sqrt((p.l1 + p.l2) / p.g):.3f} s")

    # Lagrange vs Newton-Euler at random states
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(200):
        y = np.concatenate([rng.uniform(-np.pi, np.pi, 2), rng.uniform(-3, 3, 2)])
        aL = lagrange_accelerations(y, p)
        aN, _ = newton_euler(y, p)
        worst = max(worst, np.max(np.abs(aL - aN)) / max(1.0, np.max(np.abs(aN))))
    print(f"\nLagrange vs Newton-Euler at 200 random states: max relative difference {worst:.1e}")

    # a 10 degree swing of the load in mode 1: tensions and energy
    t = np.linspace(0, 60, 6001)
    y0 = np.array([np.radians(10.0), np.radians(10.0) * phi[1, 0], 0, 0])
    Y = simulate(y0, t, p)
    E = energy(Y, p)
    E_swing = E[0] - energy(np.zeros(4), p)          # energy of the swing above the hanging equilibrium
    T1 = np.array([newton_euler(Y[:, k], p)[1][0] for k in range(0, len(t), 10)])
    W = (p.m1 + p.m2) * p.g
    print(f"10 deg mode-1 swing over 60 s: energy drift {np.max(np.abs(E - E[0])) / E_swing:.1e} of the swing energy; "
          f"cable tension {T1.min() / W:.4f}-{T1.max() / W:.4f} x static weight")
    print(f"  simulated frequency {dominant_frequency(t, Y[0]):.5f} rad/s vs linear {w[0]:.5f}")

    # chaos: lab pendulum
    q = lab_pendulum()
    for deg in (5.0, 30.0, 90.0, 120.0):
        lam, _, _ = lyapunov_estimate(np.radians([deg, deg, 0, 0]) * np.array([1, 1, 0, 0]), q, 30.0)
        print(f"lab pendulum released at th1 = th2 = {deg:5.1f} deg: finite-time Lyapunov exponent {lam:.3f} 1/s")

    # figure
    t2 = np.linspace(0, 20, 4001)
    ya = simulate(np.radians([120.0, 120.0, 0, 0]) * np.array([1, 1, 0, 0]), t2, q)
    yb = simulate(np.radians([120.0, 120.0, 0, 0]) * np.array([1, 1, 0, 0]) + np.array([1e-8, 0, 0, 0]), t2, q)
    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    x1, y1, x2, y2 = positions(ya, q)
    ax[0].plot(x2, y2, lw=0.5)
    ax[0].set_aspect("equal")
    ax[0].set_title("lab double pendulum, 120 deg release: path of m2")
    ax[0].set_xlabel("x (m)")
    ax[0].set_ylabel("y (m)")
    ax[1].semilogy(t2, np.abs(ya[1] - yb[1]) + 1e-16)
    ax[1].set_xlabel("t (s)")
    ax[1].set_ylabel("|delta theta2| (rad)")
    ax[1].set_title("two releases 1e-8 rad apart")
    ax[2].plot(t, np.degrees(Y[0]), label="theta1 (cable)")
    ax[2].plot(t, np.degrees(Y[1]), label="theta2 (sling)")
    ax[2].set_xlabel("t (s)")
    ax[2].set_ylabel("angle (deg)")
    ax[2].set_title("slung load, 10 deg swing in mode 1")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "double_pendulum.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'double_pendulum.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
