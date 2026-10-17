"""P03.10 Nonlinear Systems: reference solution.

Newton-Raphson for systems F(x) = 0, applied to the longitudinal trim of an aircraft (three unknowns:
angle of attack, elevator deflection, thrust) and to a steady heat-conduction problem whose
conductivity depends on temperature (256 unknowns):

* Newton with an analytic Jacobian, a finite-difference Jacobian and a complex-step Jacobian;
* damped Newton with an Armijo backtracking line search on the merit function ||F||^2 / 2;
* Broyden's quasi-Newton update (one residual evaluation per iteration, no Jacobian);
* continuation (warm starts) for a sweep of trim conditions over speed;
* Picard (successive substitution) against Newton on the nonlinear conduction problem.

Oracles: ``scipy.optimize.root`` (hybr / lm) and ``scipy.optimize.newton_krylov``.

Unknown scaling for trim: x = (alpha [rad], delta_e [rad], T/W), so every unknown is O(0.1-1) and the
residuals are divided by the weight (and the chord for the moment), so every equation is O(1) too.

Run:  python solution.py          (prints the verification tables, writes nonlinear_systems.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-10-nonlinear-systems
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665


# ---------------------------------------------------------------- the aircraft and its trim equations

@dataclass
class Aircraft:
    """Longitudinal linear-aerodynamics model of a regional jet (coefficients per radian)."""
    mass: float = 40_000.0      # kg
    S: float = 80.0             # m^2 reference area
    c: float = 3.0              # m mean aerodynamic chord
    CL0: float = 0.25
    CLa: float = 5.5
    CLde: float = 0.40
    CD0: float = 0.025
    k: float = 0.045            # induced-drag factor: CD = CD0 + k CL^2
    Cm0: float = 0.05
    Cma: float = -1.2
    Cmde: float = -1.3
    zT: float = 1.5             # m, thrust line below the centre of gravity (positive: nose-up thrust moment)

    @property
    def W(self) -> float:
        return self.mass * G0


def isa_density(h: float) -> float:
    """ISA troposphere density [kg/m^3] at geopotential altitude h [m] (h <= 11 km)."""
    T = 288.15 - 0.0065 * h
    p = 101325.0 * (T / 288.15) ** (G0 / (0.0065 * 287.05287))
    return p / (287.05287 * T)


def dynamic_pressure(rho: float, V: float) -> float:
    """q = rho V^2 / 2 [Pa]."""
    return 0.5 * rho * V * V


def trim_residual(x: np.ndarray, ac: Aircraft, V: float, rho: float) -> np.ndarray:
    """Steady level flight (flight-path angle zero, so the angle of attack equals the pitch attitude):

        F1 = (q S CL + T sin(alpha) - W) / W            vertical force balance
        F2 = (T cos(alpha) - q S CD) / W                 horizontal force balance
        F3 = (q S c Cm + T zT) / (W c)                   pitching-moment balance about the CG

    with x = (alpha, delta_e, T/W), CL = CL0 + CLa alpha + CLde delta_e, CD = CD0 + k CL^2,
    Cm = Cm0 + Cma alpha + Cmde delta_e. Works with complex x (for the complex-step Jacobian).
    """
    a, de, t = x
    W = ac.W
    T = t * W
    qS = dynamic_pressure(rho, V) * ac.S
    CL = ac.CL0 + ac.CLa * a + ac.CLde * de
    CD = ac.CD0 + ac.k * CL ** 2
    Cm = ac.Cm0 + ac.Cma * a + ac.Cmde * de
    return np.array([(qS * CL + T * np.sin(a) - W) / W,
                     (T * np.cos(a) - qS * CD) / W,
                     (qS * ac.c * Cm + T * ac.zT) / (W * ac.c)])


def trim_jacobian(x: np.ndarray, ac: Aircraft, V: float, rho: float) -> np.ndarray:
    """Analytic Jacobian dF/dx of ``trim_residual`` (3 x 3)."""
    a, de, t = x
    W = ac.W
    T = t * W
    qS = dynamic_pressure(rho, V) * ac.S
    CL = ac.CL0 + ac.CLa * a + ac.CLde * de
    dCD_dCL = 2.0 * ac.k * CL
    J = np.empty((3, 3))
    J[0] = [(qS * ac.CLa + T * np.cos(a)) / W, qS * ac.CLde / W, np.sin(a)]
    J[1] = [(-T * np.sin(a) - qS * dCD_dCL * ac.CLa) / W, -qS * dCD_dCL * ac.CLde / W, np.cos(a)]
    J[2] = [qS * ac.Cma / W, qS * ac.Cmde / W, ac.zT / ac.c]
    return J


def trim_initial_guess(ac: Aircraft, V: float, rho: float, crude: bool = False) -> np.ndarray:
    """A physically motivated start: the linear lift equation with delta_e = 0 for alpha, T/W = CD/CL.
    With crude=True, the uninformed start (alpha, delta_e, T/W) = (0, 0, 0) that every sweep point of the
    cold-start study and the convergence-order study use."""
    if crude:
        return np.zeros(3)
    qS = dynamic_pressure(rho, V) * ac.S
    CL = ac.W / qS
    a = (CL - ac.CL0) / ac.CLa
    CD = ac.CD0 + ac.k * CL ** 2
    return np.array([a, 0.0, CD / CL])


# ---------------------------------------------------------------- Jacobians by differencing

def fd_jacobian(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray, h: float | None = None) -> np.ndarray:
    """Forward-difference Jacobian, column j = (f(x + h_j e_j) - f(x)) / h_j with h_j = sqrt(eps) max(1, |x_j|)
    unless h is given. Costs n + 1 residual evaluations; accuracy O(h) ~ 1e-8 relative at best."""
    x = np.asarray(x, dtype=float)
    f0 = np.asarray(f(x), dtype=float)
    J = np.empty((f0.size, x.size))
    for j in range(x.size):
        hj = np.sqrt(np.finfo(float).eps) * max(1.0, abs(x[j])) if h is None else h
        xp = x.copy()
        xp[j] += hj
        J[:, j] = (np.asarray(f(xp), dtype=float) - f0) / hj
    return J


def complex_step_jacobian(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray, h: float = 1e-30) -> np.ndarray:
    """Complex-step Jacobian: column j = Im f(x + i h e_j) / h. No subtraction, so no cancellation:
    exact to round-off for any analytic f written with NumPy functions."""
    x = np.asarray(x, dtype=float)
    n = x.size
    J = None
    for j in range(n):
        xc = x.astype(complex)
        xc[j] += 1j * h
        col = np.imag(np.asarray(f(xc))) / h
        if J is None:
            J = np.empty((col.size, n))
        J[:, j] = col
    return J


# ---------------------------------------------------------------- Newton-type solvers

@dataclass
class SolveResult:
    x: np.ndarray
    converged: bool
    iterations: int
    residual_norms: np.ndarray            # ||F(x_k)||, k = 0 .. iterations
    xs: list = field(default_factory=list, repr=False)   # the iterates x_0 .. x_k
    nfev: int = 0
    step_lengths: list = field(default_factory=list, repr=False)   # damping factors used (damped Newton)


def newton(f: Callable, x0: np.ndarray, jac: Callable | None = None, tol: float = 1e-10,
           max_iter: int = 50) -> SolveResult:
    """Undamped Newton-Raphson: solve J(x_k) dx = -F(x_k), x_{k+1} = x_k + dx.

    ``jac(x)`` returns the Jacobian; None means a forward-difference Jacobian (n extra evaluations).
    Stops when ||F|| <= tol or when ||dx|| <= tol * max(1, ||x||); raises nothing on divergence but
    reports converged=False (and stops early if the iterate becomes non-finite or huge).
    """
    x = np.array(x0, dtype=float)
    F = np.asarray(f(x), dtype=float)
    nfev = 1
    norms = [np.linalg.norm(F)]
    xs = [x.copy()]
    for k in range(1, max_iter + 1):
        if norms[-1] <= tol:
            return SolveResult(x, True, k - 1, np.array(norms), xs, nfev)
        if jac is None:
            J = fd_jacobian(f, x)
            nfev += x.size
        else:
            J = np.asarray(jac(x), dtype=float)
        dx = np.linalg.solve(J, -F)
        x = x + dx
        F = np.asarray(f(x), dtype=float)
        nfev += 1
        norms.append(np.linalg.norm(F))
        xs.append(x.copy())
        if not np.all(np.isfinite(x)) or np.linalg.norm(x) > 1e12:
            return SolveResult(x, False, k, np.array(norms), xs, nfev)
        if np.linalg.norm(dx) <= tol * max(1.0, np.linalg.norm(x)) and norms[-1] <= np.sqrt(tol):
            return SolveResult(x, True, k, np.array(norms), xs, nfev)
    return SolveResult(x, bool(norms[-1] <= tol), max_iter, np.array(norms), xs, nfev)


def damped_newton(f: Callable, x0: np.ndarray, jac: Callable | None = None, tol: float = 1e-10,
                  max_iter: int = 100, c1: float = 1e-4, shrink: float = 0.5, min_step: float = 1e-6) -> SolveResult:
    """Newton with an Armijo backtracking line search on the merit function m(x) = ||F(x)||^2 / 2.

    The Newton direction is a descent direction for m (grad m = J^T F, and dx = -J^{-1} F gives
    grad m . dx = -||F||^2 < 0), so some step lambda in (0, 1] satisfies
        m(x + lambda dx) <= (1 - 2 c1 lambda) m(x).
    Start at lambda = 1 and halve until it holds. Far from the root the step is shortened; near the
    root lambda = 1 is accepted and the quadratic convergence of Newton is recovered.
    """
    x = np.array(x0, dtype=float)
    F = np.asarray(f(x), dtype=float)
    nfev = 1
    norms = [np.linalg.norm(F)]
    xs = [x.copy()]
    lambdas = []
    for k in range(1, max_iter + 1):
        if norms[-1] <= tol:
            return SolveResult(x, True, k - 1, np.array(norms), xs, nfev, lambdas)
        if jac is None:
            J = fd_jacobian(f, x)
            nfev += x.size
        else:
            J = np.asarray(jac(x), dtype=float)
        dx = np.linalg.solve(J, -F)
        m0 = 0.5 * norms[-1] ** 2
        lam = 1.0
        while True:
            x_try = x + lam * dx
            F_try = np.asarray(f(x_try), dtype=float)
            nfev += 1
            if 0.5 * np.dot(F_try, F_try) <= (1.0 - 2.0 * c1 * lam) * m0 or lam < min_step:
                break
            lam *= shrink
        x, F = x_try, F_try
        lambdas.append(lam)
        norms.append(np.linalg.norm(F))
        xs.append(x.copy())
        if np.linalg.norm(lam * dx) <= tol * max(1.0, np.linalg.norm(x)) and norms[-1] <= np.sqrt(tol):
            return SolveResult(x, True, k, np.array(norms), xs, nfev, lambdas)
    return SolveResult(x, bool(norms[-1] <= tol), max_iter, np.array(norms), xs, nfev, lambdas)


def broyden(f: Callable, x0: np.ndarray, J0: np.ndarray | None = None, tol: float = 1e-10,
            max_iter: int = 200) -> SolveResult:
    """Broyden's method ("good" update): keep an approximate Jacobian B and after each step
    s = x_new - x, y = F_new - F, update B += (y - B s) s^T / (s^T s), the least change that satisfies
    the secant condition B s = y. One residual evaluation per iteration, superlinear convergence.
    J0 defaults to a finite-difference Jacobian at x0."""
    x = np.array(x0, dtype=float)
    F = np.asarray(f(x), dtype=float)
    nfev = 1
    if J0 is None:
        B = fd_jacobian(f, x)
        nfev += x.size
    else:
        B = np.array(J0, dtype=float)
    norms = [np.linalg.norm(F)]
    xs = [x.copy()]
    for k in range(1, max_iter + 1):
        if norms[-1] <= tol:
            return SolveResult(x, True, k - 1, np.array(norms), xs, nfev)
        s = np.linalg.solve(B, -F)
        x_new = x + s
        F_new = np.asarray(f(x_new), dtype=float)
        nfev += 1
        y = F_new - F
        B = B + np.outer(y - B @ s, s) / np.dot(s, s)
        x, F = x_new, F_new
        norms.append(np.linalg.norm(F))
        xs.append(x.copy())
        if not np.all(np.isfinite(x)):
            return SolveResult(x, False, k, np.array(norms), xs, nfev)
        if np.linalg.norm(s) <= tol * max(1.0, np.linalg.norm(x)) and norms[-1] <= np.sqrt(tol):
            return SolveResult(x, True, k, np.array(norms), xs, nfev)
    return SolveResult(x, bool(norms[-1] <= tol), max_iter, np.array(norms), xs, nfev)


def iterate_errors(xs: list, x_star: np.ndarray) -> np.ndarray:
    """||x_k - x*|| for each iterate (x* a converged solution), used for the order estimates."""
    return np.array([np.linalg.norm(np.asarray(x) - x_star) for x in xs])


def convergence_order(errors: np.ndarray, floor: float = 1e-13) -> float:
    """Observed order p from e_{k+1} ~ C e_k^p: least-squares slope of log e_{k+1} against log e_k over
    the consecutive pairs whose errors are both above ``floor`` (round-off) and below 1."""
    e = np.asarray(errors, dtype=float)
    pairs = [(e[k], e[k + 1]) for k in range(len(e) - 1) if floor < e[k + 1] and e[k] < 1.0 and e[k] > floor]
    if len(pairs) < 2:
        raise ValueError("need at least two usable error pairs to estimate an order")
    a = np.log([p[0] for p in pairs])
    b = np.log([p[1] for p in pairs])
    return float(np.polyfit(a, b, 1)[0])


# ---------------------------------------------------------------- trim: solving and sweeping

def solve_trim(ac: Aircraft, V: float, h: float, method: str = "newton", x0: np.ndarray | None = None,
               tol: float = 1e-12, crude: bool = False) -> SolveResult:
    """Trim at speed V [m/s] and altitude h [m]. method: 'newton' (analytic J), 'newton-fd', 'broyden', 'damped'.
    x0 defaults to ``trim_initial_guess`` (the crude zero start if crude=True)."""
    rho = isa_density(h)
    f = lambda x: trim_residual(x, ac, V, rho)
    jac = lambda x: trim_jacobian(x, ac, V, rho)
    if x0 is None:
        x0 = trim_initial_guess(ac, V, rho, crude)
    if method == "newton":
        return newton(f, x0, jac, tol=tol)
    if method == "newton-fd":
        return newton(f, x0, None, tol=tol)
    if method == "broyden":
        return broyden(f, x0, tol=tol)
    if method == "damped":
        return damped_newton(f, x0, jac, tol=tol)
    raise ValueError(f"unknown method {method!r}")


def trim_sweep(ac: Aircraft, Vs: np.ndarray, h: float, warm_start: bool = True, tol: float = 1e-12,
               crude: bool = True) -> tuple[np.ndarray, np.ndarray]:
    """Trim over a range of speeds. With warm_start the previous trim is the next initial guess
    (continuation); otherwise each point starts from ``trim_initial_guess`` (crude zero start by default).
    Returns (X (len(Vs), 3), iterations)."""
    X, its = [], []
    x_prev = None
    for V in Vs:
        res = solve_trim(ac, V, h, "newton", x0=x_prev if warm_start else None, tol=tol, crude=crude)
        if not res.converged:
            raise RuntimeError(f"trim did not converge at V = {V}")
        X.append(res.x)
        its.append(res.iterations)
        x_prev = res.x
    return np.array(X), np.array(its)


def stall_speed(ac: Aircraft, h: float, CL_max: float = 1.6) -> float:
    """V at which trimmed flight would need CL_max with delta_e = 0 and T sin(alpha) neglected: sqrt(2W/(rho S CL_max))."""
    return float(np.sqrt(2.0 * ac.W / (isa_density(h) * ac.S * CL_max)))


# ---------------------------------------------------------------- a case where undamped Newton fails

def arctan_residual(x: np.ndarray) -> np.ndarray:
    """F(x) = arctan(x): a single root at 0; Newton from |x0| > 1.3917 diverges (the tangent overshoots)."""
    return np.arctan(np.asarray(x, dtype=float))


def arctan_jacobian(x: np.ndarray) -> np.ndarray:
    return np.array([[1.0 / (1.0 + float(x[0]) ** 2)]])


# ---------------------------------------------------------------- nonlinear heat conduction (256 unknowns)

def face_conductances(kfield: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Harmonic-mean conductances on the vertical faces kx[J, I] (between (J, I) and (J, I+1)) and the
    horizontal faces ky[J, I] (between (J, I) and (J+1, I)) of the (n+2) x (n+2) grid."""
    kx = 2.0 * kfield[:, :-1] * kfield[:, 1:] / (kfield[:, :-1] + kfield[:, 1:])
    ky = 2.0 * kfield[:-1, :] * kfield[1:, :] / (kfield[:-1, :] + kfield[1:, :])
    return kx, ky


def conduction_matrix(kfield: np.ndarray, n: int, h: float) -> np.ndarray:
    """Dense five-point matrix of -div(k grad T) on the n x n interior grid for a conductivity field on
    the (n+2) x (n+2) grid (harmonic means at faces); a small-size twin of P03.09's sparse assembly."""
    kx, ky = face_conductances(kfield)
    west, east = kx[1:-1, :-1], kx[1:-1, 1:]          # (n, n): faces left and right of each interior point
    south, north = ky[:-1, 1:-1], ky[1:, 1:-1]
    N = n * n
    A = np.zeros((N, N))
    idx = np.arange(N).reshape(n, n)                   # idx[j, i] = i + n j
    A[idx, idx] = (west + east + south + north) / h ** 2
    A[idx[:, :-1], idx[:, 1:]] = -east[:, :-1] / h ** 2
    A[idx[:, 1:], idx[:, :-1]] = -west[:, 1:] / h ** 2
    A[idx[:-1, :], idx[1:, :]] = -north[:-1, :] / h ** 2
    A[idx[1:, :], idx[:-1, :]] = -south[1:, :] / h ** 2
    return A


def conduction_boundary_rhs(kfield: np.ndarray, n: int, h: float, T_wall: float) -> np.ndarray:
    """Wall-temperature contribution k_face T_wall / h^2 for the boundary neighbours."""
    kx, ky = face_conductances(kfield)
    b = np.zeros((n, n))
    b[:, 0] += kx[1:-1, 0]
    b[:, -1] += kx[1:-1, -1]
    b[0, :] += ky[0, 1:-1]
    b[-1, :] += ky[-1, 1:-1]
    return b.reshape(n * n) * T_wall / h ** 2


@dataclass
class HotPlate:
    """Square plate whose conductance depends on temperature: k(T) = k0 (1 + beta (T - T_wall))."""
    n: int = 16
    L: float = 0.2
    k0: float = 0.34          # W/K (conductivity x thickness)
    beta: float = -0.004      # 1/K (conductance falls with temperature: a stiffer nonlinearity than aluminium's)
    q: float = 5000.0         # W/m^2 absorbed
    T_wall: float = 300.0

    @property
    def h(self) -> float:
        return self.L / (self.n + 1)


def plate_kfield(T: np.ndarray, pl: HotPlate) -> np.ndarray:
    """Conductance on the (n+2)^2 grid from the interior temperatures T (n^2,), wall ring at T_wall."""
    n = pl.n
    Tfull = np.full((n + 2, n + 2), pl.T_wall)
    Tfull[1:-1, 1:-1] = np.asarray(T).reshape(n, n)
    return pl.k0 * (1.0 + pl.beta * (Tfull - pl.T_wall))


def plate_residual(T: np.ndarray, pl: HotPlate) -> np.ndarray:
    """F(T) = A(T) T - b(T): the discrete -div(k(T) grad T) - q with Dirichlet walls."""
    kf = plate_kfield(T, pl)
    A = conduction_matrix(kf, pl.n, pl.h)
    b = conduction_boundary_rhs(kf, pl.n, pl.h, pl.T_wall) + pl.q
    return A @ np.asarray(T, dtype=float) - b


def picard_plate(pl: HotPlate, tol: float = 1e-8, max_iter: int = 200) -> SolveResult:
    """Successive substitution: freeze k at the current T, solve the linear problem A(T_k) T_{k+1} = b(T_k).
    Converges linearly (order 1) when the nonlinearity is mild; fails when it is strong."""
    n = pl.n
    T = np.full(n * n, pl.T_wall)
    norms = [np.linalg.norm(plate_residual(T, pl))]
    xs = [T.copy()]
    nfev = 1
    for k in range(1, max_iter + 1):
        if norms[-1] <= tol:
            return SolveResult(T, True, k - 1, np.array(norms), xs, nfev)
        kf = plate_kfield(T, pl)
        A = conduction_matrix(kf, n, pl.h)
        b = conduction_boundary_rhs(kf, n, pl.h, pl.T_wall) + pl.q
        T = np.linalg.solve(A, b)
        norms.append(np.linalg.norm(plate_residual(T, pl)))
        nfev += 1
        xs.append(T.copy())
    return SolveResult(T, bool(norms[-1] <= tol), max_iter, np.array(norms), xs, nfev)


def newton_plate(pl: HotPlate, tol: float = 1e-8, damped: bool = True, max_iter: int = 50) -> SolveResult:
    """Newton on the plate with a forward-difference Jacobian (n^2 + 1 residual evaluations per
    iteration, each an O(n^2) assembly) and a dense solve; damped by default."""
    f = lambda T: plate_residual(T, pl)
    T0 = np.full(pl.n * pl.n, pl.T_wall)
    return damped_newton(f, T0, None, tol=tol, max_iter=max_iter) if damped else newton(f, T0, None, tol=tol, max_iter=max_iter)


def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's experiments, print the tables and save ``nonlinear_systems.png``."""
    import time
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import root, newton_krylov

    out_dir = Path(out_dir)
    ac = Aircraft()
    h_alt, V = 5000.0, 200.0
    rho = isa_density(h_alt)
    print(f"regional jet, W = {ac.W / 1e3:.1f} kN, h = {h_alt / 1e3:.0f} km (rho = {rho:.4f} kg/m^3), V = {V} m/s, q S = {dynamic_pressure(rho, V) * ac.S / 1e3:.1f} kN")

    # ---- Jacobians
    x0 = trim_initial_guess(ac, V, rho)
    f = lambda x: trim_residual(x, ac, V, rho)
    Ja, Jfd, Jcs = trim_jacobian(x0, ac, V, rho), fd_jacobian(f, x0), complex_step_jacobian(f, x0)
    print(f"initial guess alpha = {np.degrees(x0[0]):.3f} deg, delta_e = 0, T/W = {x0[2]:.4f}; ||F(x0)|| = {np.linalg.norm(f(x0)):.3e}")
    print(f"Jacobian check: |J_fd - J| max = {np.max(np.abs(Jfd - Ja)):.1e}, |J_cs - J| max = {np.max(np.abs(Jcs - Ja)):.1e}; kappa(J) = {np.linalg.cond(Ja):.1f}")
    for hh in (1e-2, 1e-4, 1e-6, 1e-8, 1e-10, 1e-12):
        print(f"  forward difference h = {hh:.0e}: max error {np.max(np.abs(fd_jacobian(f, x0, hh) - Ja)):.2e}")

    # ---- Newton, Broyden, scipy
    res_smart = solve_trim(ac, V, h_alt, "newton")
    res = solve_trim(ac, V, h_alt, "newton", crude=True)
    x_star = res.x
    print(f"\nNewton from the physical guess: {res_smart.iterations} iterations; from (0, 0, 0):")
    print(f"Newton (analytic J): {res.iterations} iterations, alpha = {np.degrees(x_star[0]):.4f} deg, delta_e = {np.degrees(x_star[1]):.4f} deg, "
          f"T = {x_star[2] * ac.W / 1e3:.3f} kN (T/W = {x_star[2]:.5f}); residual norms {np.array2string(res.residual_norms, precision=2)}")
    errs = iterate_errors(res.xs, x_star)
    print(f"  errors ||x_k - x*||: {np.array2string(errs, precision=2)}; observed order {convergence_order(errs):.2f}")
    for method in ("newton-fd", "broyden", "damped"):
        r = solve_trim(ac, V, h_alt, method, crude=True)
        e = iterate_errors(r.xs, x_star)
        try:
            order = f"{convergence_order(e):.2f}"
        except ValueError:
            order = "n/a"
        print(f"{method:10s}: {r.iterations} iterations, {r.nfev} residual evaluations, |x - x*| = {np.linalg.norm(r.x - x_star):.1e}, "
              f"errors {np.array2string(e, precision=1)}, order {order}")
    sol = root(f, np.zeros(3), jac=lambda x: trim_jacobian(x, ac, V, rho), method="hybr", tol=1e-12)
    print(f"scipy.optimize.root (hybr): {sol.nfev} evaluations, |x - x*| = {np.linalg.norm(sol.x - x_star):.1e}, success {sol.success}")
    # a hand check: thrust = drag at the trimmed CL
    qS = dynamic_pressure(rho, V) * ac.S
    CL = ac.W / qS
    CD = ac.CD0 + ac.k * CL ** 2
    print(f"hand check: CL = W/(qS) = {CL:.4f} (code CL = {ac.CL0 + ac.CLa * x_star[0] + ac.CLde * x_star[1]:.4f}), "
          f"T ~ qS CD = {qS * CD / 1e3:.3f} kN, L/D = {CL / CD:.2f}")

    # ---- sweep over speed: continuation
    Vs = np.arange(110.0, 261.0, 5.0)
    t0 = time.perf_counter(); Xw, itw = trim_sweep(ac, Vs, h_alt, warm_start=True); tw = time.perf_counter() - t0
    t0 = time.perf_counter(); Xc, itc = trim_sweep(ac, Vs, h_alt, warm_start=False); tc = time.perf_counter() - t0
    print(f"\nsweep V = {Vs[0]:.0f}..{Vs[-1]:.0f} m/s ({Vs.size} points): warm start {itw.sum()} Newton iterations ({tw * 1e3:.1f} ms), "
          f"cold start {itc.sum()} ({tc * 1e3:.1f} ms); max |x_warm - x_cold| = {np.max(np.abs(Xw - Xc)):.1e}")
    i_min = np.argmin(Xw[:, 2])
    print(f"minimum thrust {Xw[i_min, 2] * ac.W / 1e3:.2f} kN at V = {Vs[i_min]:.0f} m/s (max L/D speed ~ {np.sqrt(2 * ac.W / (rho * ac.S)) * (ac.k / ac.CD0) ** 0.25:.0f} m/s "
          f"from CD0 = k CL^2); alpha from {np.degrees(Xw[0, 0]):.2f} to {np.degrees(Xw[-1, 0]):.2f} deg, delta_e from {np.degrees(Xw[0, 1]):.2f} to {np.degrees(Xw[-1, 1]):.2f} deg; "
          f"stall speed (CL_max 1.6) {stall_speed(ac, h_alt):.1f} m/s")

    # ---- Newton failing and the fix
    for x0s in (1.0, 1.5, 3.0):
        ru = newton(arctan_residual, np.array([x0s]), arctan_jacobian, tol=1e-12, max_iter=20)
        rd = damped_newton(arctan_residual, np.array([x0s]), arctan_jacobian, tol=1e-12)
        print(f"arctan from x0 = {x0s}: undamped {'converged' if ru.converged else 'diverged'} after {ru.iterations} (|x| = {abs(ru.x[0]):.2e}); "
              f"damped converged in {rd.iterations} with steps {np.array2string(np.array(rd.step_lengths), precision=3)}")

    # ---- the hot plate: Picard vs Newton vs scipy newton_krylov
    pl = HotPlate()
    t0 = time.perf_counter(); rp = picard_plate(pl); tp = time.perf_counter() - t0
    t0 = time.perf_counter(); rn = newton_plate(pl); tn = time.perf_counter() - t0
    t0 = time.perf_counter(); Tk = newton_krylov(lambda T: plate_residual(T, pl), np.full(pl.n ** 2, pl.T_wall), f_tol=1e-8); tk = time.perf_counter() - t0
    T_star = rn.x
    ep, en = iterate_errors(rp.xs, T_star), iterate_errors(rn.xs, T_star)
    print(f"\nhot plate n = {pl.n} ({pl.n ** 2} unknowns), k(T) = k0(1 + beta (T - T_wall)), beta = {pl.beta} 1/K, q = {pl.q} W/m^2: "
          f"max rise {T_star.max() - pl.T_wall:.2f} K (linear k0 would give {np.max(np.linalg.solve(conduction_matrix(plate_kfield(np.full(pl.n ** 2, pl.T_wall), pl), pl.n, pl.h), conduction_boundary_rhs(plate_kfield(np.full(pl.n ** 2, pl.T_wall), pl), pl.n, pl.h, pl.T_wall) + pl.q)) - pl.T_wall:.2f} K); "
          f"k at the hottest point {pl.k0 * (1 + pl.beta * (T_star.max() - pl.T_wall)):.3f} W/K")
    print(f"Picard: {rp.iterations} iterations ({tp * 1e3:.0f} ms), residuals {np.array2string(rp.residual_norms, precision=1)}, order {convergence_order(ep / ep[0]):.2f}; "
          f"ratio of successive errors {np.array2string(ep[2:6] / ep[1:5], precision=3)}")
    print(f"Newton (FD Jacobian, damped): {rn.iterations} iterations, {rn.nfev} residual evaluations ({tn * 1e3:.0f} ms), residuals {np.array2string(rn.residual_norms, precision=1)}, "
          f"steps {rn.step_lengths}, order {convergence_order(en / en[0]):.2f}")
    print(f"scipy newton_krylov: {tk * 1e3:.0f} ms, max |T - T_newton| = {np.max(np.abs(Tk - T_star)):.1e} K; Picard vs Newton {np.max(np.abs(rp.x - T_star)):.1e} K")
    for beta in (-0.002, -0.004, -0.006, -0.0075):
        plb = HotPlate(beta=beta)
        rpb = picard_plate(plb, max_iter=400)
        rnb = newton_plate(plb)
        print(f"  beta = {beta}: Picard {'converged' if rpb.converged else 'FAILED'} in {rpb.iterations}, Newton {rnb.iterations} iterations (steps {rnb.step_lengths}), max rise {rnb.x.max() - 300:.1f} K, "
              f"min k/k0 {1 + beta * (rnb.x.max() - 300):.3f}")
    plx = HotPlate(q=10000.0, beta=-0.0075)
    rpx = picard_plate(plx, max_iter=15)
    rnx = newton_plate(plx, max_iter=15)
    print(f"  q = {plx.q:.0f}, beta = {plx.beta}: Picard {'converged' if rpx.converged else 'FAILED'} (residual {rpx.residual_norms[-1]:.1e} after {rpx.iterations}), "
          f"damped Newton {'converged' if rnx.converged else 'FAILED'} (residual {rnx.residual_norms[-1]:.1e}, last step {rnx.step_lengths[-1]:.1e}); "
          f"the linear k(T) would have to reach {1 + plx.beta * (rnx.x.max() - 300):.2f} k0: no physical solution exists")

    # ---- figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    ax.plot(Vs, np.degrees(Xw[:, 0]), label="alpha [deg]")
    ax.plot(Vs, np.degrees(Xw[:, 1]), label="delta_e [deg]")
    ax.plot(Vs, Xw[:, 2] * ac.W / 1e3 / 10, label="thrust [kN] / 10")
    ax.set_xlabel("V [m/s]"); ax.set_title("trim sweep at 5 km (warm-started Newton)"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[0, 1]
    for method, mk in (("newton", "o-"), ("newton-fd", "s--"), ("broyden", "^-"), ("damped", "x:")):
        r = solve_trim(ac, V, h_alt, method, crude=True)
        e = iterate_errors(r.xs, x_star)
        ax.semilogy(np.arange(e.size), np.maximum(e, 1e-17), mk, label=method)
    ax.set_xlabel("iteration"); ax.set_ylabel("||x_k - x*||"); ax.set_title("trim at 200 m/s: convergence"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1, 0]
    ax.semilogy(np.arange(ep.size), np.maximum(ep, 1e-13), "o-", label=f"Picard (order {convergence_order(ep / ep[0]):.2f})")
    ax.semilogy(np.arange(en.size), np.maximum(en, 1e-13), "s-", label=f"Newton (order {convergence_order(en / en[0]):.2f})")
    ax.set_xlabel("iteration"); ax.set_ylabel("||T_k - T*|| [K]"); ax.set_title("hot plate, 256 unknowns"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    ax = axes[1, 1]
    xx = np.linspace(-4, 4, 400)
    ax.plot(xx, np.arctan(xx), "k")
    ru = newton(arctan_residual, np.array([1.5]), arctan_jacobian, max_iter=3)
    for a_, b_ in zip(ru.xs[:-1], ru.xs[1:]):
        ax.plot([a_[0], b_[0]], [np.arctan(a_[0]), 0.0], "r-", lw=0.8)
    ax.plot([x[0] for x in ru.xs], np.arctan([x[0] for x in ru.xs]), "ro", ms=4, label="undamped from 1.5")
    rd = damped_newton(arctan_residual, np.array([1.5]), arctan_jacobian)
    ax.plot([x[0] for x in rd.xs], np.arctan([x[0] for x in rd.xs]), "b.-", ms=6, label="damped from 1.5")
    ax.set_ylim(-1.7, 1.7); ax.set_title("F(x) = arctan x: why Newton needs damping"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
    fig.suptitle("P03.10  Newton-Raphson for trim and for nonlinear conduction")
    fig.tight_layout()
    fig.savefig(out_dir / "nonlinear_systems.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'nonlinear_systems.png'}")


if __name__ == "__main__":
    main()
