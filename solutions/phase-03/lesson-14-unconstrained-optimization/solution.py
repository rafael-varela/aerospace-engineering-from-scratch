"""P03.14 Unconstrained Optimization: reference solution.

Minimisers written from scratch: steepest descent, Newton's method and the BFGS quasi-Newton method, each
with a line search (Armijo backtracking, or a strong-Wolfe bracketing/zoom search), plus a central-difference
gradient. They are exercised on the Rosenbrock function and a quadratic (where every property has a known
answer) and then on an aerospace design problem: the cruise Mach number and altitude that minimise the fuel
burned per kilometre by a jet transport, with a compressible drag polar and an ISA atmosphere.
``scipy.optimize.minimize`` is the oracle.

Conventions: x is a 1-D NumPy array; f(x) returns a float; grad(x) returns the gradient array; results are
dictionaries with keys 'x', 'f', 'grad_norm', 'iterations', 'path' (array of iterates, one per row),
'converged' and, for BFGS, 'H' (the final inverse-Hessian approximation). The design variables of the cruise
problem are x = (M, h / 10 km) so that both are of order one.

Run:  python solution.py          (prints the verification tables, writes bfgs_cruise.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-14-unconstrained-optimization
"""
from __future__ import annotations

from pathlib import Path
from typing import Callable

import numpy as np

Func = Callable[[np.ndarray], float]
Grad = Callable[[np.ndarray], np.ndarray]


# ---------------------------------------------------------------- test functions

def rosenbrock(x: np.ndarray) -> float:
    """Return f = 100 (x2 - x1^2)^2 + (1 - x1)^2: the banana valley, minimum 0 at (1, 1)."""
    x = np.asarray(x, dtype=float)
    return float(100.0 * (x[1] - x[0] ** 2) ** 2 + (1.0 - x[0]) ** 2)


def rosenbrock_grad(x: np.ndarray) -> np.ndarray:
    """Return the gradient of ``rosenbrock``."""
    x = np.asarray(x, dtype=float)
    return np.array([-400.0 * x[0] * (x[1] - x[0] ** 2) - 2.0 * (1.0 - x[0]), 200.0 * (x[1] - x[0] ** 2)])


def rosenbrock_hess(x: np.ndarray) -> np.ndarray:
    """Return the Hessian of ``rosenbrock``."""
    x = np.asarray(x, dtype=float)
    return np.array([[1200.0 * x[0] ** 2 - 400.0 * x[1] + 2.0, -400.0 * x[0]], [-400.0 * x[0], 200.0]])


def quadratic(Q: np.ndarray, b: np.ndarray) -> tuple[Func, Grad]:
    """Return (f, grad) for f(x) = 1/2 x^T Q x - b^T x, whose minimiser is Q^{-1} b (Q symmetric positive definite)."""
    Q = np.asarray(Q, dtype=float)
    b = np.asarray(b, dtype=float)

    def f(x: np.ndarray) -> float:
        x = np.asarray(x, dtype=float)
        return float(0.5 * x @ Q @ x - b @ x)

    def g(x: np.ndarray) -> np.ndarray:
        return Q @ np.asarray(x, dtype=float) - b

    return f, g


# ---------------------------------------------------------------- derivatives

def gradient_fd(f: Func, x: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """Return the central-difference gradient (f(x + h e_i) - f(x - h e_i)) / (2h): error O(h^2) + eps/h.

    2n function evaluations. The truncation error is h^2 f'''/6 and the rounding error about eps |f| / h,
    so the best h is near (eps)^(1/3) ~ 1e-5 to 1e-6 for well-scaled x; the lesson measures both.
    """
    x = np.asarray(x, dtype=float)
    g = np.empty(x.size)
    for i in range(x.size):
        e = np.zeros(x.size)
        e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2.0 * h)
    return g


def hessian_fd(f: Func, x: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Return a symmetric central-difference Hessian, O(h^2), from 1 + 2n + 2n(n-1) function evaluations."""
    x = np.asarray(x, dtype=float)
    n = x.size
    H = np.empty((n, n))
    f0 = f(x)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h
        H[i, i] = (f(x + ei) - 2.0 * f0 + f(x - ei)) / h ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = h
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4.0 * h ** 2)
            H[j, i] = H[i, j]
    return H


# ---------------------------------------------------------------- line searches

def line_search_backtracking(f: Func, x: np.ndarray, p: np.ndarray, g: np.ndarray, alpha0: float = 1.0,
                             rho: float = 0.5, c1: float = 1e-4, max_iter: int = 60) -> float:
    """Return a step alpha satisfying the Armijo (sufficient decrease) condition f(x + a p) <= f(x) + c1 a g.p.

    Start from alpha0 and shrink by rho until the condition holds. Requires a descent direction (g.p < 0);
    raise ``ValueError`` otherwise.
    """
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=float)
    slope = float(g @ p)
    if slope >= 0.0:
        raise ValueError("not a descent direction")
    f0 = f(x)
    alpha = alpha0
    for _ in range(max_iter):
        if f(x + alpha * p) <= f0 + c1 * alpha * slope:
            return alpha
        alpha *= rho
    return alpha


def line_search_wolfe(f: Func, grad: Grad, x: np.ndarray, p: np.ndarray, alpha0: float = 1.0, c1: float = 1e-4,
                      c2: float = 0.9, alpha_max: float = 50.0, max_iter: int = 40) -> float:
    """Return a step satisfying the strong Wolfe conditions by bracketing and zooming (bisection/interpolation).

    Sufficient decrease: phi(a) <= phi(0) + c1 a phi'(0); curvature: |phi'(a)| <= c2 |phi'(0)|, with
    phi(a) = f(x + a p). The curvature condition guarantees s.y > 0 for the BFGS update, which keeps the
    inverse-Hessian approximation positive definite. Raise ``ValueError`` for a non-descent direction.
    """
    x = np.asarray(x, dtype=float)
    p = np.asarray(p, dtype=float)
    phi0 = f(x)
    dphi0 = float(grad(x) @ p)
    if dphi0 >= 0.0:
        raise ValueError("not a descent direction")

    def phi(a: float) -> float:
        return f(x + a * p)

    def dphi(a: float) -> float:
        return float(grad(x + a * p) @ p)

    def zoom(lo: float, hi: float, phi_lo: float) -> float:
        for _ in range(max_iter):
            a = 0.5 * (lo + hi)                      # bisection inside the bracket
            pa = phi(a)
            if pa > phi0 + c1 * a * dphi0 or pa >= phi_lo:
                hi = a
            else:
                da = dphi(a)
                if abs(da) <= -c2 * dphi0:
                    return a
                if da * (hi - lo) >= 0.0:
                    hi = lo
                lo, phi_lo = a, pa
        return lo

    a_prev, phi_prev = 0.0, phi0
    a = alpha0
    for i in range(max_iter):
        pa = phi(a)
        if pa > phi0 + c1 * a * dphi0 or (i > 0 and pa >= phi_prev):
            return zoom(a_prev, a, phi_prev)
        da = dphi(a)
        if abs(da) <= -c2 * dphi0:
            return a
        if da >= 0.0:
            return zoom(a, a_prev, pa)
        a_prev, phi_prev = a, pa
        a = min(2.0 * a, alpha_max)
    return a


# ---------------------------------------------------------------- minimisers

def _finish(f: Func, path: list[np.ndarray], g: np.ndarray, converged: bool, extra: dict | None = None) -> dict:
    """Pack the result dictionary (x, f, grad_norm, iterations, path, converged, plus ``extra``) from the iterate list."""
    out = {"x": path[-1], "f": f(path[-1]), "grad_norm": float(np.linalg.norm(g)), "iterations": len(path) - 1,
           "path": np.array(path), "converged": converged}
    if extra:
        out.update(extra)
    return out


def steepest_descent(f: Func, grad: Grad, x0: np.ndarray, tol: float = 1e-6, max_iter: int = 5000,
                     exact_quadratic: np.ndarray | None = None) -> dict:
    """Minimise f from x0 along -grad with a Wolfe line search; stop when ||grad|| <= tol.

    With ``exact_quadratic=Q`` the exact step for f = 1/2 x^T Q x - b^T x, alpha = g.g / (g^T Q g), is used
    instead, which reproduces the textbook linear convergence rate ((kappa - 1)/(kappa + 1))^2 in f.
    """
    x = np.asarray(x0, dtype=float).copy()
    path = [x.copy()]
    g = grad(x)
    for _ in range(max_iter):
        if np.linalg.norm(g) <= tol:
            return _finish(f, path, g, True)
        p = -g
        if exact_quadratic is not None:
            alpha = float(g @ g / (g @ exact_quadratic @ g))
        else:
            alpha = line_search_wolfe(f, grad, x, p, alpha0=min(1.0, 1.0 / max(np.linalg.norm(g), 1e-12)))
        x = x + alpha * p
        g = grad(x)
        path.append(x.copy())
    return _finish(f, path, g, bool(np.linalg.norm(g) <= tol))


def newton_method(f: Func, grad: Grad, hess: Callable[[np.ndarray], np.ndarray], x0: np.ndarray, tol: float = 1e-8,
                  max_iter: int = 100) -> dict:
    """Minimise f by damped Newton: p = -H^{-1} g (falling back to -g if H is not positive definite), Armijo step.

    Quadratic convergence near a minimiser with a positive-definite Hessian; one linear solve per iteration.
    """
    x = np.asarray(x0, dtype=float).copy()
    path = [x.copy()]
    g = grad(x)
    for _ in range(max_iter):
        if np.linalg.norm(g) <= tol:
            return _finish(f, path, g, True)
        H = hess(x)
        try:
            np.linalg.cholesky(H)
            p = -np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            p = -g
        alpha = line_search_backtracking(f, x, p, g)
        x = x + alpha * p
        g = grad(x)
        path.append(x.copy())
    return _finish(f, path, g, bool(np.linalg.norm(g) <= tol))


def bfgs_update(H: np.ndarray, s: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the BFGS update of the inverse-Hessian approximation H from the step s and gradient change y.

    H+ = (I - rho s y^T) H (I - rho y s^T) + rho s s^T, rho = 1 / (y.s). It is the closest matrix to H (in a
    weighted norm) that satisfies the secant equation H+ y = s, and it stays symmetric positive definite
    whenever y.s > 0. Raise ``ValueError`` if y.s <= 0.
    """
    s = np.asarray(s, dtype=float)
    y = np.asarray(y, dtype=float)
    ys = float(y @ s)
    if ys <= 0.0:
        raise ValueError("curvature condition y.s > 0 violated")
    rho = 1.0 / ys
    I = np.eye(s.size)
    V = I - rho * np.outer(s, y)
    return V @ np.asarray(H, dtype=float) @ V.T + rho * np.outer(s, s)


def bfgs(f: Func, grad: Grad, x0: np.ndarray, tol: float = 1e-8, max_iter: int = 500, H0: np.ndarray | None = None,
         exact_quadratic: np.ndarray | None = None) -> dict:
    """Minimise f by BFGS: p = -H g, strong-Wolfe line search, secant update of H; stop when ||grad|| <= tol.

    H0 defaults to the identity, rescaled after the first step by (y.s)/(y.y) (Nocedal-Wright's heuristic).
    With ``exact_quadratic=Q`` the exact step alpha = -g.p / (p^T Q p) is used, in which case BFGS on an
    n-dimensional quadratic terminates in at most n iterations with H = Q^{-1}. Updates with y.s <= 0
    (possible only when the line search fails) are skipped.
    """
    x = np.asarray(x0, dtype=float).copy()
    n = x.size
    H = np.eye(n) if H0 is None else np.asarray(H0, dtype=float).copy()
    path = [x.copy()]
    g = grad(x)
    for k in range(max_iter):
        if np.linalg.norm(g) <= tol:
            return _finish(f, path, g, True, {"H": H})
        p = -H @ g
        if g @ p >= 0.0:                        # lost positive definiteness through rounding: reset
            H = np.eye(n)
            p = -g
        if exact_quadratic is not None:
            alpha = float(-(g @ p) / (p @ exact_quadratic @ p))
        else:
            alpha = line_search_wolfe(f, grad, x, p)
        s = alpha * p
        x_new = x + s
        g_new = grad(x_new)
        y = g_new - g
        ys = float(y @ s)
        if ys > 1e-12 * np.linalg.norm(y) * np.linalg.norm(s):
            if k == 0 and H0 is None:
                H = (ys / (y @ y)) * np.eye(n)
            H = bfgs_update(H, s, y)
        x, g = x_new, g_new
        path.append(x.copy())
    return _finish(f, path, g, bool(np.linalg.norm(g) <= tol), {"H": H})


# ---------------------------------------------------------------- convergence diagnostics

def error_sequence(path: np.ndarray, x_star: np.ndarray) -> np.ndarray:
    """Return ||x_k - x*|| for every iterate in ``path``."""
    return np.linalg.norm(np.asarray(path, dtype=float) - np.asarray(x_star, dtype=float), axis=1)


def convergence_order(errors: np.ndarray, floor: float = 1e-13) -> np.ndarray:
    """Return the estimated order p_k = log(e_{k+1}/e_k) / log(e_k/e_{k-1}) from consecutive errors above ``floor``.

    1 means linear, 2 quadratic; BFGS gives values that drift upward (superlinear: the ratio e_{k+1}/e_k
    itself tends to zero). Errors at or below ``floor`` (rounding) are dropped.
    """
    e = np.asarray(errors, dtype=float)
    e = e[e > floor]
    if e.size < 3:
        return np.array([])
    r = np.log(e[1:] / e[:-1])
    return r[1:] / r[:-1]


def linear_rate_bound(Q: np.ndarray) -> float:
    """Return ((kappa - 1)/(kappa + 1))^2 with kappa = cond(Q): the steepest-descent bound on f_{k+1}/f_k (exact line search)."""
    lam = np.linalg.eigvalsh(np.asarray(Q, dtype=float))
    kappa = lam[-1] / lam[0]
    return float(((kappa - 1.0) / (kappa + 1.0)) ** 2)


# ---------------------------------------------------------------- the cruise problem

# GIVEN
def isa(h: float) -> tuple[float, float, float, float]:
    """Return (T [K], p [Pa], rho [kg/m^3], a [m/s]) of the International Standard Atmosphere at height h [m] (P03.03).

    Linear lapse of 6.5 K/km to 11 km, isothermal 216.65 K above; gamma = 1.4, R = 287.05 J/(kg K).
    """
    R, g0 = 287.05, 9.80665
    if h < 11000.0:
        T = 288.15 - 0.0065 * h
        p = 101325.0 * (T / 288.15) ** (g0 / (R * 0.0065))
    else:
        T = 216.65
        p = 22632.06 * np.exp(-g0 * (h - 11000.0) / (R * T))
    rho = p / (R * T)
    return float(T), float(p), float(rho), float(np.sqrt(1.4 * R * T))


# GIVEN
def cruise_aircraft() -> dict[str, float]:
    """Return a narrow-body jet transport in mid-cruise, in SI units.

    W [N] weight, S [m^2] wing area, CD0 zero-lift drag, K induced-drag factor (CD = CD0 + K CL^2), M_cr0 and
    dM_dCL the drag-divergence Mach number model M_dd = M_cr0 - dM_dCL * CL, k_wave the wave-drag coefficient
    (CD_wave = k_wave (M - M_dd)^4 for M > M_dd), cT0 [kg/(N s)] the thrust-specific fuel consumption at sea
    level, scaled by sqrt(T/288.15) with altitude. Representative numbers, not a specific type.
    """
    return {"W": 600.0e3, "S": 125.0, "CD0": 0.020, "K": 0.045, "M_cr0": 0.82, "dM_dCL": 0.20, "k_wave": 30.0,
            "cT0": 1.6e-5}


def fuel_per_km(x: np.ndarray, ac: dict[str, float]) -> float:
    """Return the fuel burned per kilometre [kg/km] in steady level cruise at x = (Mach, altitude / 10 km).

    Lift = weight fixes CL = W / (q S) with q = rho V^2 / 2 and V = M a; drag D = q S (CD0 + K CL^2 + CD_wave);
    thrust = drag; fuel flow = cT D; fuel per distance = cT D / V. Raise ``ValueError`` for M <= 0 or h < 0.
    """
    x = np.asarray(x, dtype=float)
    M, h = float(x[0]), float(x[1]) * 1.0e4
    if M <= 0.0 or h < 0.0:
        raise ValueError("Mach must be positive and altitude non-negative")
    T, _, rho, a = isa(h)
    V = M * a
    q = 0.5 * rho * V ** 2
    CL = ac["W"] / (q * ac["S"])
    M_dd = ac["M_cr0"] - ac["dM_dCL"] * CL
    CD = ac["CD0"] + ac["K"] * CL ** 2 + ac["k_wave"] * max(M - M_dd, 0.0) ** 4
    D = q * ac["S"] * CD
    cT = ac["cT0"] * np.sqrt(T / 288.15)
    return float(cT * D / V * 1000.0)


def cruise_breakdown(x: np.ndarray, ac: dict[str, float]) -> dict[str, float]:
    """Return the physical quantities at a cruise point: V, CL, L/D, drag, drag-divergence Mach, fuel flow, altitude."""
    x = np.asarray(x, dtype=float)
    M, h = float(x[0]), float(x[1]) * 1.0e4
    T, _, rho, a = isa(h)
    V = M * a
    q = 0.5 * rho * V ** 2
    CL = ac["W"] / (q * ac["S"])
    M_dd = ac["M_cr0"] - ac["dM_dCL"] * CL
    CD = ac["CD0"] + ac["K"] * CL ** 2 + ac["k_wave"] * max(M - M_dd, 0.0) ** 4
    D = q * ac["S"] * CD
    return {"M": M, "h_m": h, "V_mps": V, "CL": CL, "L_over_D": CL / CD, "drag_N": D, "M_dd": M_dd,
            "fuel_flow_kgps": ac["cT0"] * np.sqrt(T / 288.15) * D, "fuel_per_km": fuel_per_km(x, ac)}


def optimal_cruise(ac: dict[str, float], x0: np.ndarray | None = None, tol: float = 1e-7) -> dict:
    """Return the BFGS result for the minimum-fuel cruise point, using the central-difference gradient (h = 1e-6).

    Points outside the model's domain (M <= 0 or h < 0), which a long trial step can reach, are given f = inf
    so that the line search rejects them instead of raising.
    """
    x0 = np.array([0.70, 0.80]) if x0 is None else np.asarray(x0, dtype=float)

    def f(x: np.ndarray) -> float:
        return fuel_per_km(x, ac) if (x[0] > 0.0 and x[1] >= 0.0) else np.inf

    g = lambda x: gradient_fd(f, x, 1e-6)                             # noqa: E731
    return bfgs(f, g, x0, tol=tol)


# ---------------------------------------------------------------- main

# GIVEN
def main(out_dir: str | Path = ".") -> None:
    """Compare the three minimisers on Rosenbrock, check BFGS on a quadratic, optimise the cruise point, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize

    out_dir = Path(out_dir)
    x0 = np.array([-1.2, 1.0])
    x_star = np.array([1.0, 1.0])
    print("Rosenbrock from (-1.2, 1): iterations, f, ||g||, final convergence-order estimates")
    results = {}
    for name, res in (("steepest descent", steepest_descent(rosenbrock, rosenbrock_grad, x0, tol=1e-6, max_iter=20000)),
                      ("Newton", newton_method(rosenbrock, rosenbrock_grad, rosenbrock_hess, x0, tol=1e-10)),
                      ("BFGS", bfgs(rosenbrock, rosenbrock_grad, x0, tol=1e-10))):
        results[name] = res
        e = error_sequence(res["path"], x_star)
        orders = convergence_order(e)
        tail = ", ".join(f"{o:.2f}" for o in orders[-3:]) if orders.size else "n/a"
        print(f"  {name:17s} {res['iterations']:6d} iterations, f = {res['f']:.2e}, ||g|| = {res['grad_norm']:.1e}, orders {tail}")
    e_b = error_sequence(results["BFGS"]["path"], x_star)
    print("  BFGS error ratios e_{k+1}/e_k over the last 6 steps:", np.round((e_b[1:] / e_b[:-1])[-6:], 3))
    sp = minimize(rosenbrock, x0, jac=rosenbrock_grad, method="BFGS", options={"gtol": 1e-10})
    print(f"  scipy BFGS: {sp.nit} iterations, f = {sp.fun:.2e}; our final inverse Hessian vs exact H^-1 at x*:"
          f" max |diff| = {np.max(np.abs(results['BFGS']['H'] - np.linalg.inv(rosenbrock_hess(x_star)))):.2e}")

    Q = np.array([[10.0, 2.0, 0.0], [2.0, 4.0, 1.0], [0.0, 1.0, 1.0]])
    b = np.array([1.0, -2.0, 0.5])
    fq, gq = quadratic(Q, b)
    res = bfgs(fq, gq, np.zeros(3), tol=1e-12, exact_quadratic=Q)
    print(f"\nquadratic, 3 variables, exact line search: BFGS finished in {res['iterations']} iterations,"
          f" |H - Q^-1| = {np.max(np.abs(res['H'] - np.linalg.inv(Q))):.1e}, |x - Q^-1 b| = {np.linalg.norm(res['x'] - np.linalg.solve(Q, b)):.1e}")
    sd = steepest_descent(fq, gq, np.zeros(3), tol=1e-10, exact_quadratic=Q)
    fs = np.array([fq(xk) for xk in sd["path"]]) - fq(np.linalg.solve(Q, b))
    fs = fs[fs > 1e-12]                                      # drop the rounding floor
    print(f"steepest descent on the same quadratic: {sd['iterations']} iterations; worst observed f-ratio {np.max(fs[1:] / fs[:-1]):.4f}"
          f" against the bound ((k-1)/(k+1))^2 = {linear_rate_bound(Q):.4f} (kappa = {np.linalg.cond(Q):.1f})")

    f = lambda x: rosenbrock(x)                                        # noqa: E731
    xg = np.array([0.3, -0.7])
    print("\ncentral-difference gradient at (0.3, -0.7): |error| against the analytic gradient")
    for h in (1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-8, 1e-10):
        print(f"  h = {h:.0e}: {np.max(np.abs(gradient_fd(f, xg, h) - rosenbrock_grad(xg))):.2e}")

    ac = cruise_aircraft()
    fc = lambda x: fuel_per_km(x, ac)                                  # noqa: E731
    res = optimal_cruise(ac)
    sp = minimize(fc, np.array([0.70, 0.80]), method="BFGS", options={"gtol": 1e-9})
    bd = cruise_breakdown(res["x"], ac)
    print(f"\nminimum-fuel cruise (BFGS with FD gradient, {res['iterations']} iterations, ||g|| = {res['grad_norm']:.1e}):")
    print(f"  M = {bd['M']:.4f}, h = {bd['h_m']:.0f} m, V = {bd['V_mps']:.1f} m/s, CL = {bd['CL']:.3f}, L/D = {bd['L_over_D']:.2f},"
          f" M_dd = {bd['M_dd']:.3f}, drag = {bd['drag_N']/1e3:.1f} kN, fuel = {bd['fuel_per_km']:.4f} kg/km")
    print(f"  scipy BFGS: M = {sp.x[0]:.4f}, h = {sp.x[1]*1e4:.0f} m, fuel = {sp.fun:.4f} kg/km, {sp.nit} iterations")
    Hc = hessian_fd(fc, res["x"], 1e-3)
    print(f"  FD Hessian at the optimum: eigenvalues {np.round(np.linalg.eigvalsh(Hc), 3)} (both positive: a minimum),"
          f" condition number {np.linalg.cond(Hc):.1f}")
    for lab, x in (("M = 0.78, h = 10 km", [0.78, 1.0]), ("M = 0.80, h = 11 km", [0.80, 1.1]), ("optimum", res["x"])):
        print(f"  {lab:22s}: {fuel_per_km(np.asarray(x), ac):.4f} kg/km")

    # Figure.
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    ax = axes[0]
    X, Y = np.meshgrid(np.linspace(-1.5, 1.5, 300), np.linspace(-0.5, 1.5, 300))
    Z = 100 * (Y - X ** 2) ** 2 + (1 - X) ** 2
    ax.contour(X, Y, np.log10(Z + 1e-3), levels=20, cmap="Greys")
    for name, style in (("steepest descent", "r-"), ("Newton", "g.-"), ("BFGS", "b.-")):
        p = results[name]["path"]
        ax.plot(p[:, 0], p[:, 1], style, lw=0.8, ms=3, label=f"{name} ({results[name]['iterations']} it.)")
    ax.plot(1, 1, "k*", ms=10)
    ax.set_xlabel("x1"); ax.set_ylabel("x2"); ax.set_title("Rosenbrock: iterates"); ax.legend(fontsize=8)
    ax = axes[1]
    Ms = np.linspace(0.65, 0.86, 120)
    hs = np.linspace(7000, 13500, 120)
    F = np.array([[fuel_per_km(np.array([m, h / 1e4]), ac) for m in Ms] for h in hs])
    cs = ax.contour(Ms, hs / 1e3, F, levels=np.arange(2.2, 3.6, 0.1))
    ax.clabel(cs, fmt="%.1f", fontsize=7)
    p = res["path"]
    ax.plot(p[:, 0], p[:, 1] * 10, "b.-", lw=1, label="BFGS path")
    ax.plot(res["x"][0], res["x"][1] * 10, "r*", ms=12, label="optimum")
    ax.set_xlabel("Mach"); ax.set_ylabel("altitude [km]"); ax.set_title("fuel per km [kg/km]"); ax.legend(fontsize=8)
    fig.suptitle("P03.14  Unconstrained optimisation: BFGS from scratch")
    fig.tight_layout()
    fig.savefig(out_dir / "bfgs_cruise.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'bfgs_cruise.png'}")


if __name__ == "__main__":
    main()
