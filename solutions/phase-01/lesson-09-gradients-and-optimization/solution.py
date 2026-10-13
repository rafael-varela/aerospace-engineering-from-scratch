"""P01.09 Gradients, Directional Derivatives & Optimization: reference solution.

Drag of an aircraft in level flight as a function of two design variables, cruise speed V and
wing aspect ratio AR:

    D(V, AR) = 1/2 rho V**2 S C_D0 + 2 K(AR) W(AR)**2 / (rho V**2 S),
    K(AR) = 1 / (pi e AR),            W(AR) = W0 + c_w AR**2   (wing-weight surrogate).

The gradient (analytic and by central differences), directional derivatives, the Hessian and
the second-derivative test, steepest-descent and Newton searches (scaled variables), and the
Lagrange-multiplier solution of "minimum drag subject to lift = weight" in the (V, C_L) plane,
which reproduces the classical minimum-drag speed V_md = sqrt(2W/(rho S)) (K/C_D0)**(1/4).

Run:  python solution.py          (prints the optimisation table, writes drag_landscape.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-09-gradients-and-optimization
"""
from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import sympy as sp


@dataclass(frozen=True)
class Aircraft:
    """Aircraft and air parameters in SI units.

    rho : air density, kg/m**3
    S : wing area, m**2
    cd0 : zero-lift drag coefficient
    e : Oswald span-efficiency factor
    W0 : weight of everything except the aspect-ratio penalty, N
    c_w : wing-weight surrogate coefficient, N per unit AR**2 (W = W0 + c_w AR**2)
    """

    rho: float = 1.225
    S: float = 16.0
    cd0: float = 0.025
    e: float = 0.8
    W0: float = 13000.0
    c_w: float = 40.0


DEFAULT_AIRCRAFT = Aircraft()


def induced_factor(AR: np.ndarray | float, ac: Aircraft = DEFAULT_AIRCRAFT) -> np.ndarray | float:
    """Return K = 1 / (pi e AR). Raise ``ValueError`` for AR <= 0."""
    AR = np.asarray(AR, dtype=float)
    if np.any(AR <= 0.0):
        raise ValueError("aspect ratio must be positive")
    out = 1.0 / (np.pi * ac.e * AR)
    return float(out) if out.ndim == 0 else out


def weight(AR: np.ndarray | float, ac: Aircraft = DEFAULT_AIRCRAFT) -> np.ndarray | float:
    """Return the weight surrogate W = W0 + c_w AR**2 in newtons."""
    AR = np.asarray(AR, dtype=float)
    out = ac.W0 + ac.c_w * AR**2
    return float(out) if out.ndim == 0 else out


def drag(V: np.ndarray | float, AR: np.ndarray | float, ac: Aircraft = DEFAULT_AIRCRAFT) -> np.ndarray | float:
    """Return the level-flight drag D(V, AR) in newtons; V and AR broadcast. Raise ``ValueError`` for V <= 0."""
    V = np.asarray(V, dtype=float)
    if np.any(V <= 0.0):
        raise ValueError("speed must be positive")
    q = 0.5 * ac.rho * V**2
    out = q * ac.S * ac.cd0 + induced_factor(AR, ac) * weight(AR, ac) ** 2 / (q * ac.S)
    return float(out) if np.ndim(out) == 0 else out


def drag_gradient(V: float, AR: float, ac: Aircraft = DEFAULT_AIRCRAFT) -> np.ndarray:
    """Return grad D = (dD/dV, dD/dAR) analytically, units N/(m/s) and N per unit AR.

    dD/dV = rho V S C_D0 - 4 K W**2 / (rho V**3 S);
    dD/dAR = (2/(rho V**2 S)) [K' W**2 + 2 K W W'] with K' = -K/AR, W' = 2 c_w AR.
    """
    K, W = induced_factor(AR, ac), weight(AR, ac)
    dK, dW = -K / AR, 2.0 * ac.c_w * AR
    q = 0.5 * ac.rho * V**2
    dD_dV = ac.rho * V * ac.S * ac.cd0 - 4.0 * K * W**2 / (ac.rho * V**3 * ac.S)
    dD_dAR = (dK * W**2 + 2.0 * K * W * dW) / (q * ac.S)
    return np.array([dD_dV, dD_dAR])


def gradient_fd(f: Callable[[np.ndarray], float], x: np.ndarray, h: np.ndarray | float) -> np.ndarray:
    """Return the central-difference gradient of f at x; ``h`` is a scalar or one step per coordinate.

    Raise ``ValueError`` if any step is not positive.
    """
    x = np.asarray(x, dtype=float)
    h = np.broadcast_to(np.asarray(h, dtype=float), x.shape)
    if np.any(h <= 0.0):
        raise ValueError("steps must be positive")
    g = np.empty_like(x)
    for i in range(x.size):
        e = np.zeros_like(x)
        e[i] = h[i]
        g[i] = (f(x + e) - f(x - e)) / (2.0 * h[i])
    return g


def hessian_fd(f: Callable[[np.ndarray], float], x: np.ndarray, h: np.ndarray | float) -> np.ndarray:
    """Return the central-difference Hessian (matrix of second partials) of f at x, symmetrised."""
    x = np.asarray(x, dtype=float)
    h = np.broadcast_to(np.asarray(h, dtype=float), x.shape)
    if np.any(h <= 0.0):
        raise ValueError("steps must be positive")
    n = x.size
    H = np.empty((n, n))
    f0 = f(x)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = h[i]
        H[i, i] = (f(x + ei) - 2.0 * f0 + f(x - ei)) / h[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = h[j]
            H[i, j] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4.0 * h[i] * h[j])
            H[j, i] = H[i, j]
    return H


def directional_derivative(grad: np.ndarray, direction: np.ndarray) -> float:
    """Return D_u f = grad . u with u = direction / |direction|. Raise ``ValueError`` for a zero direction."""
    d = np.asarray(direction, dtype=float)
    n = np.linalg.norm(d)
    if n == 0.0:
        raise ValueError("direction must be non-zero")
    return float(np.dot(np.asarray(grad, dtype=float), d / n))


def classify_critical_point(H: np.ndarray) -> str:
    """Return 'minimum', 'maximum', 'saddle' or 'degenerate' from the Hessian's eigenvalues (second-derivative test)."""
    H = np.asarray(H, dtype=float)
    w = np.linalg.eigvalsh(0.5 * (H + H.T))
    tol = 1e-12 * max(1.0, np.max(np.abs(w)))
    if np.all(w > tol):
        return "minimum"
    if np.all(w < -tol):
        return "maximum"
    if np.any(np.abs(w) <= tol):
        return "degenerate"
    return "saddle"


def gradient_descent(f: Callable[[np.ndarray], float], grad: Callable[[np.ndarray], np.ndarray], x0: np.ndarray,
                     step: float = 0.1, tol: float = 1e-8, max_iter: int = 10000) -> tuple[np.ndarray, np.ndarray]:
    """Steepest descent with a fixed step and backtracking: x <- x - step * grad, halving step while f does not decrease.

    Return (x_opt, path) with ``path`` the (n_iter + 1, dim) array of iterates. Stop when |grad| < tol
    or after ``max_iter`` steps. Raise ``ValueError`` if step <= 0.
    """
    if step <= 0.0:
        raise ValueError("step must be positive")
    x = np.asarray(x0, dtype=float).copy()
    path = [x.copy()]
    fx = f(x)
    for _ in range(max_iter):
        g = grad(x)
        if np.linalg.norm(g) < tol:
            break
        s = step
        while True:
            x_new = x - s * g
            try:
                f_new = f(x_new)
            except ValueError:          # stepped outside the function's domain (e.g. V <= 0): shorten
                f_new = np.inf
            if np.isfinite(f_new) and f_new < fx:
                break
            s *= 0.5
            if s < 1e-16:
                return x, np.array(path)
        x, fx = x_new, f_new
        path.append(x.copy())
    return x, np.array(path)


def newton_minimize(f: Callable[[np.ndarray], float], grad: Callable[[np.ndarray], np.ndarray],
                    hess: Callable[[np.ndarray], np.ndarray], x0: np.ndarray, tol: float = 1e-10,
                    max_iter: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Newton's method for a minimum: x <- x - H^{-1} grad. Return (x_opt, path). Stop when |grad| < tol."""
    x = np.asarray(x0, dtype=float).copy()
    path = [x.copy()]
    for _ in range(max_iter):
        g = grad(x)
        if np.linalg.norm(g) < tol:
            break
        x = x - np.linalg.solve(hess(x), g)
        path.append(x.copy())
    return x, np.array(path)


def scaled_drag(ac: Aircraft = DEFAULT_AIRCRAFT, V_ref: float = 50.0, AR_ref: float = 10.0, D_ref: float = 1000.0) -> tuple[Callable[[np.ndarray], float], Callable[[np.ndarray], np.ndarray]]:
    """Return (f, grad) of the dimensionless drag f(z) = D(V_ref z0, AR_ref z1) / D_ref and its gradient in z."""
    def f(z: np.ndarray) -> float:
        return float(drag(V_ref * z[0], AR_ref * z[1], ac) / D_ref)

    def g(z: np.ndarray) -> np.ndarray:
        gr = drag_gradient(V_ref * z[0], AR_ref * z[1], ac)
        return np.array([gr[0] * V_ref, gr[1] * AR_ref]) / D_ref

    return f, g


def minimum_drag_speed(W: float, AR: float, ac: Aircraft = DEFAULT_AIRCRAFT) -> float:
    """Return V_md = sqrt(2 W / (rho S)) (K / C_D0)**(1/4): the speed of minimum drag at fixed weight and AR."""
    if W <= 0.0:
        raise ValueError("weight must be positive")
    K = induced_factor(AR, ac)
    return float(np.sqrt(2.0 * W / (ac.rho * ac.S)) * (K / ac.cd0) ** 0.25)


def optimum_analytic(ac: Aircraft = DEFAULT_AIRCRAFT) -> tuple[float, float, float]:
    """Return (V*, AR*, D*) of the unconstrained minimum of D(V, AR).

    Minimising over V first gives D_min(AR) = 2 W(AR) sqrt(K(AR) C_D0); then dD_min/dAR = 0 gives
    AR* = sqrt(W0 / (3 c_w)) for W = W0 + c_w AR**2 and K = 1/(pi e AR).
    """
    AR = np.sqrt(ac.W0 / (3.0 * ac.c_w))
    V = minimum_drag_speed(weight(AR, ac), AR, ac)
    return V, float(AR), float(drag(V, AR, ac))


def lagrange_min_drag(W: float, ac: Aircraft = DEFAULT_AIRCRAFT) -> tuple[float, float, float]:
    """Solve min D(V, C_L) subject to L(V, C_L) = W with Lagrange multipliers, in the (V, C_L) plane.

    D = q S (C_D0 + K C_L**2), L = q S C_L, q = 1/2 rho V**2, with K at the aircraft's optimum AR.
    grad D = lambda grad L gives C_L* = sqrt(C_D0 / K); the constraint then gives V*.
    Return (V*, C_L*, lambda*). lambda* = dD/dW is the marginal drag cost of weight.
    """
    _, AR, _ = optimum_analytic(ac)
    K = induced_factor(AR, ac)
    cl = np.sqrt(ac.cd0 / K)
    V = np.sqrt(2.0 * W / (ac.rho * ac.S * cl))
    lam = 2.0 * K * cl
    return float(V), float(cl), float(lam)


def lagrange_sympy() -> dict[str, sp.Expr]:
    """Return the symbolic Lagrange solution {'C_L': ..., 'V': ..., 'lambda': ...} of min D s.t. L = W.

    Symbols: V, C_L, lambda, rho, S, K, C_D0, W (all positive).
    """
    V, CL, lam, rho, S, K, cd0, W = sp.symbols("V C_L lambda rho S K C_D0 W", positive=True)
    q = sp.Rational(1, 2) * rho * V**2
    D = q * S * (cd0 + K * CL**2)
    L = q * S * CL
    eqs = [sp.diff(D, V) - lam * sp.diff(L, V), sp.diff(D, CL) - lam * sp.diff(L, CL), L - W]
    sols = sp.solve(eqs, [V, CL, lam], dict=True)
    (sol,) = [s for s in sols if all(v.is_positive for v in s.values())]
    return {"C_L": sp.simplify(sol[CL]), "V": sp.simplify(sol[V]), "lambda": sp.simplify(sol[lam])}


def main(out_dir: str | Path = ".") -> None:
    """Print the optimisation table and save drag_landscape.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize

    out_dir = Path(out_dir)
    ac = DEFAULT_AIRCRAFT
    V_s, AR_s, D_s = optimum_analytic(ac)
    print(f"Aircraft: S = {ac.S} m^2, C_D0 = {ac.cd0}, e = {ac.e}, W = {ac.W0} + {ac.c_w} AR^2 N")
    print(f"Analytic optimum: V* = {V_s:.3f} m/s, AR* = {AR_s:.4f}, D* = {D_s:.3f} N (W = {weight(AR_s, ac):.0f} N, L/D = {weight(AR_s, ac)/D_s:.2f})")
    x0 = np.array([35.0, 6.0])
    g_an = drag_gradient(*x0, ac)
    g_fd = gradient_fd(lambda x: drag(x[0], x[1], ac), x0, np.array([1e-3, 1e-4]))
    print(f"At V = {x0[0]} m/s, AR = {x0[1]}: D = {drag(*x0, ac):.2f} N, grad = {g_an} (FD {g_fd})")
    u = -g_an / np.linalg.norm(g_an)
    print(f"  steepest descent direction {np.round(u, 4)}, slope {directional_derivative(g_an, u):.3f}; along +V only {directional_derivative(g_an, [1, 0]):.3f}; along +AR only {directional_derivative(g_an, [0, 1]):.3f}")

    f_raw = lambda x: drag(x[0], x[1], ac)
    g_raw = lambda x: drag_gradient(x[0], x[1], ac)
    for step in (1e-3, 1e-1):
        x_raw, p_raw = gradient_descent(f_raw, g_raw, x0, step=step, tol=1e-6, max_iter=20000)
        print(f"  raw-variable descent, step {step:g}: {len(p_raw)-1} iterations -> V = {x_raw[0]:.3f}, AR = {x_raw[1]:.4f}, |grad| = {np.linalg.norm(g_raw(x_raw)):.1e}")
    f_s, g_s = scaled_drag(ac)
    z0 = np.array([x0[0] / 50.0, x0[1] / 10.0])
    z_opt, path = gradient_descent(f_s, g_s, z0, step=0.2, tol=1e-9, max_iter=20000)
    print(f"  scaled descent: {len(path)-1} iterations -> V = {50*z_opt[0]:.3f}, AR = {10*z_opt[1]:.4f}")
    hz = lambda z: hessian_fd(f_s, z, 1e-4)
    z_n, path_n = newton_minimize(f_s, g_s, hz, z0)
    H = hessian_fd(f_raw, np.array([V_s, AR_s]), np.array([1e-2, 1e-3]))
    print(f"  Newton: {len(path_n)-1} iterations -> V = {50*z_n[0]:.3f}, AR = {10*z_n[1]:.4f}")
    print(f"  Hessian at optimum (N/(m/s)^2, N/(m/s), N): {np.round(H, 4).tolist()}, eigenvalues {np.round(np.linalg.eigvalsh(H), 5)} -> {classify_critical_point(H)}")
    res = minimize(f_raw, x0, method="BFGS", jac=g_raw, options={"gtol": 1e-10})
    print(f"  scipy.optimize.minimize (BFGS): V = {res.x[0]:.3f}, AR = {res.x[1]:.4f}, {res.nit} iterations")
    W = weight(AR_s, ac)
    Vl, cl, lam = lagrange_min_drag(W, ac)
    print(f"Lagrange (min D s.t. L = W = {W:.0f} N at AR* = {AR_s:.3f}): C_L* = {cl:.4f}, V* = {Vl:.3f} m/s, lambda = dD/dW = {lam:.4f}")
    print(f"  V_md formula: {minimum_drag_speed(W, AR_s, ac):.3f} m/s; sympy: {lagrange_sympy()}")

    Vg = np.linspace(25.0, 75.0, 201)
    ARg = np.linspace(3.0, 20.0, 171)
    VV, AA = np.meshgrid(Vg, ARg, indexing="ij")
    fig, ax = plt.subplots(figsize=(7, 5))
    cs = ax.contour(VV, AA, drag(VV, AA, ac), levels=[700, 750, 800, 900, 1000, 1200, 1500, 2000, 3000], cmap="viridis")
    ax.clabel(cs, fmt="%.0f N", fontsize=7)
    ax.plot(50 * path[:, 0], 10 * path[:, 1], "r.-", ms=3, lw=0.7, label=f"scaled steepest descent ({len(path)-1} steps)")
    ax.plot(50 * path_n[:, 0], 10 * path_n[:, 1], "bs-", ms=4, lw=0.8, label=f"Newton ({len(path_n)-1} steps)")
    ax.plot(V_s, AR_s, "k*", ms=12, label="analytic optimum")
    ARline = np.linspace(3.0, 20.0, 50)
    ax.plot([minimum_drag_speed(weight(a, ac), a, ac) for a in ARline], ARline, "k--", lw=0.8, label="V_md(AR): valley floor")
    ax.set_xlabel("cruise speed V [m/s]")
    ax.set_ylabel("aspect ratio AR")
    ax.set_title("Drag landscape D(V, AR) with descent paths")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(out_dir / "drag_landscape.png", dpi=150)
    print(f"saved {out_dir / 'drag_landscape.png'}")


if __name__ == "__main__":
    main()
