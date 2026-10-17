"""P03.15 Constrained Optimization: starter code.

Sequential quadratic programming (SQP) written from scratch, applied to the sizing of a thin-walled
aluminium strut (minimum mass under stress, Euler buckling, local shell buckling and gauge limits),
with ``scipy.optimize.minimize`` (SLSQP, trust-constr) as the oracle.

Conventions:
    * A problem is  minimise f(x)  subject to  g_i(x) >= 0, i = 1..m.  The Lagrangian is
      L(x, lam) = f(x) - lam . g(x)  with lam >= 0, so the KKT conditions read
      grad f = J^T lam,  g >= 0,  lam >= 0,  lam_i g_i = 0.
    * Constraints are *normalised*: each g_i is written as 1 - demand/capacity, so every g_i is O(1)
      and a violation of 0.01 means "1 % over the limit" whatever the physical units.
    * Design variables are scaled: x = (r / r_ref, t / t_ref) with r_ref = 50 mm, t_ref = 1 mm, so both
      components are O(1). ``unscale`` returns (r, t) in metres.
    * The QP subproblem  min 1/2 p^T B p + grad^T p  s.t.  g + J p >= 0  is solved by enumerating active
      sets (exact for small m, exponential in m: the lesson explains why libraries do not do this).

Given: the dataclasses, ``unscale``, ``strut_problem``, the two textbook problems and ``main``;
everything else is yours.

Run:  python main.py          (prints the sizing table, writes strut_sqp.png)
Test: pytest curriculum/phase-03/lesson-15-constrained-optimization
"""
from __future__ import annotations

import itertools
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

# ---------------------------------------------------------------- the strut sizing problem (GIVEN)


@dataclass(frozen=True)
class StrutParams:
    """Parameters of the tubular strut. SI units.

    A straight aluminium tube of length L carries a compressive axial load P with pinned ends. The
    local shell-buckling stress uses the classical cylinder result E t / (r sqrt(3 (1 - nu^2)))
    multiplied by a knock-down factor ``gamma`` for imperfections (a common design practice: thin
    shells buckle well below the classical load).
    """

    P: float = 40.0e3            # axial compressive load [N]
    L: float = 2.5               # length between pinned ends [m]
    E: float = 68.9e9            # Young's modulus, 6061-T6 aluminium [Pa]
    sigma_allow: float = 240.0e6  # allowable compressive stress [Pa]
    rho: float = 2700.0          # density [kg/m^3]
    nu: float = 0.33             # Poisson's ratio
    gamma: float = 0.5           # knock-down factor on classical shell buckling
    t_min: float = 0.5e-3        # minimum manufacturable wall thickness [m]
    r_min: float = 10.0e-3       # minimum mean radius (fittings) [m]
    r_ref: float = 50.0e-3       # scaling of the radius design variable [m]
    t_ref: float = 1.0e-3        # scaling of the thickness design variable [m]

    @property
    def k_local(self) -> float:
        """Coefficient k in the local buckling stress sigma_cr = k E t / r."""
        return self.gamma / np.sqrt(3.0 * (1.0 - self.nu ** 2))


@dataclass
class Problem:
    """A smooth inequality-constrained problem: minimise f(x) s.t. g(x) >= 0.

    ``hess_lagr(x, lam)`` returns the Hessian of the Lagrangian f - lam . g; it may be None (then SQP
    must use a quasi-Newton approximation).
    """

    f: Callable[[np.ndarray], float]
    grad: Callable[[np.ndarray], np.ndarray]
    cons: Callable[[np.ndarray], np.ndarray]
    jac: Callable[[np.ndarray], np.ndarray]
    n: int
    m: int
    hess_lagr: Callable[[np.ndarray, np.ndarray], np.ndarray] | None = None
    names: tuple[str, ...] = field(default_factory=tuple)


@dataclass
class SQPResult:
    """What ``sqp`` returns: the point, multipliers, objective, history of iterates and KKT residuals."""

    x: np.ndarray
    lam: np.ndarray
    fun: float
    n_iter: int
    converged: bool
    history: list[np.ndarray]
    kkt_history: list[float]


def unscale(x: np.ndarray, params: StrutParams) -> tuple[float, float]:
    """Return (r, t) in metres from the scaled design vector x = (r / r_ref, t / t_ref).

    Both are floored at 1 micrometre so that a line search which overshoots to r <= 0 or t <= 0 sees a
    huge (finite) constraint violation instead of a division by zero or a negative mass.
    """
    return max(float(x[0] * params.r_ref), 1e-6), max(float(x[1] * params.t_ref), 1e-6)


# ---------------------------------------------------------------- the strut model


def strut_mass(x: np.ndarray, params: StrutParams) -> float:
    """Return the strut mass [kg], rho * L * 2 pi r t, for the scaled design vector x."""
    raise NotImplementedError


def strut_mass_grad(x: np.ndarray, params: StrutParams) -> np.ndarray:
    """Return d(mass)/dx for the scaled variables (chain rule through r_ref, t_ref)."""
    raise NotImplementedError


def strut_constraints(x: np.ndarray, params: StrutParams) -> np.ndarray:
    """Return the five normalised constraints g(x) >= 0 of the strut, in this order:

    0 stress:        1 - P / (2 pi r t sigma_allow)
    1 Euler:         1 - P L^2 / (pi^3 E r^3 t)          (P <= pi^2 E I / L^2 with I = pi r^3 t)
    2 local shell:   1 - P / (2 pi k E t^2)              (P / (2 pi r t) <= k E t / r)
    3 gauge:         t / t_min - 1
    4 radius:        r / r_min - 1
    """
    raise NotImplementedError


def strut_jacobian(x: np.ndarray, params: StrutParams) -> np.ndarray:
    """Return the 5 x 2 Jacobian dg/dx of ``strut_constraints`` for the scaled variables (analytic)."""
    raise NotImplementedError


def strut_hessian_lagrangian(x: np.ndarray, lam: np.ndarray, params: StrutParams) -> np.ndarray:
    """Return the Hessian of L = mass - lam . g in the scaled variables (analytic)."""
    raise NotImplementedError


def strut_problem(params: StrutParams | None = None) -> Problem:
    """Return the strut sizing problem as a ``Problem`` in the scaled variables."""
    p = params or StrutParams()
    return Problem(
        f=lambda x: strut_mass(x, p),
        grad=lambda x: strut_mass_grad(x, p),
        cons=lambda x: strut_constraints(x, p),
        jac=lambda x: strut_jacobian(x, p),
        n=2, m=5,
        hess_lagr=lambda x, lam: strut_hessian_lagrangian(x, lam, p),
        names=("stress", "Euler", "local", "gauge", "radius"),
    )


def analytic_optimum(params: StrutParams | None = None) -> tuple[np.ndarray, tuple[int, int], float]:
    """Return (x*, active pair, mass*) by enumerating vertices of the problem in log variables.

    Every constraint of the strut is a monomial inequality r^a t^b >= C, which in u = ln r, v = ln t is
    linear, and the objective ln(r t) = u + v is linear too: the problem is a linear program in (u, v)
    and its minimum sits at a vertex where two constraints are active. Solve each 2 x 2 pair in
    closed form, discard infeasible vertices, keep the lightest. Raises if no vertex is feasible.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- KKT machinery


def kkt_residual(problem: Problem, x: np.ndarray, lam: np.ndarray) -> float:
    """Return the infinity-norm KKT residual: max of stationarity |grad f - J^T lam|,
    primal infeasibility |min(g, 0)|, dual infeasibility |min(lam, 0)| and complementarity |lam g|."""
    raise NotImplementedError


def multipliers_least_squares(problem: Problem, x: np.ndarray, active: list[int] | np.ndarray) -> np.ndarray:
    """Return the multiplier vector (length m) that best satisfies grad f = J_A^T lam_A in least squares,
    with zeros on the inactive constraints. At a KKT point with LICQ the fit is exact."""
    raise NotImplementedError


def solve_qp(H: np.ndarray, c: np.ndarray, A: np.ndarray, b: np.ndarray, tol: float = 1e-10) -> tuple[np.ndarray, np.ndarray]:
    """Solve  min 1/2 p^T H p + c^T p  s.t.  A p >= b  by enumerating active sets.

    For each candidate active set W (at most n constraints, so that A_W can have full row rank) solve
    the KKT system [[H, -A_W^T], [A_W, 0]] [p; lam_W] = [-c; b_W], and accept the first p with
    lam_W >= 0 and A p >= b (all others get lam = 0). For H positive definite that KKT point is the
    unique global minimiser. Cost is exponential in m: fine for m <= 12, which is all this lesson needs;
    production codes use active-set or interior-point iterations instead. Raises ``ValueError`` if H is
    not positive definite and ``RuntimeError`` if no KKT point exists (infeasible linearisation).
    """
    raise NotImplementedError


def damped_bfgs_update(B: np.ndarray, s: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return the BFGS update of B with Powell's damping, which keeps B positive definite even when
    s^T y <= 0 (common in constrained problems, where the Lagrangian need not be convex)."""
    raise NotImplementedError


def merit_l1(problem: Problem, x: np.ndarray, mu: float) -> float:
    """Return the l1 merit function f(x) + mu * sum of constraint violations max(0, -g_i)."""
    raise NotImplementedError


def _make_positive_definite(H: np.ndarray, floor: float = 1e-6) -> np.ndarray:
    """Return H with eigenvalues below ``floor`` raised to ``floor`` (modified Newton)."""
    w, V = np.linalg.eigh(0.5 * (H + H.T))
    return (V * np.maximum(w, floor)) @ V.T


def sqp(problem: Problem, x0: np.ndarray, hessian: str = "bfgs", tol: float = 1e-8,
        max_iter: int = 100) -> SQPResult:
    """Minimise ``problem`` by sequential quadratic programming from x0.

    Each iteration linearises the constraints, builds a quadratic model of the Lagrangian with B
    (``hessian='bfgs'``: damped BFGS from the identity; ``'exact'``: problem.hess_lagr made positive
    definite), solves the QP for the step p and the new multipliers, and takes x + alpha p with a
    backtracking (Armijo) line search on the l1 merit function whose penalty mu exceeds the largest
    multiplier. Stops when the KKT residual is below ``tol``.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- penalty and augmented Lagrangian


def penalty_minimize(problem: Problem, mu: float, x0: np.ndarray, gtol: float = 1e-10) -> np.ndarray:
    """Return the minimiser of the quadratic penalty function f + (mu/2) sum min(g_i, 0)^2.

    The inner unconstrained minimisation uses BFGS (``scipy.optimize.minimize``), the method of
    P03.14. The penalty minimiser approaches the constrained optimum as mu -> infinity with error
    O(1/mu), while the Hessian's condition number grows like mu: that trade-off is the lesson.
    """
    raise NotImplementedError


def augmented_lagrangian(problem: Problem, x0: np.ndarray, mu: float = 10.0, n_outer: int = 30,
                         tol: float = 1e-8) -> tuple[np.ndarray, np.ndarray, int]:
    """Minimise ``problem`` by the augmented Lagrangian (method of multipliers) for inequalities.

    Outer iteration k: minimise  L_A(x) = f(x) + (mu/2) sum_i [max(0, lam_i/mu - g_i(x))^2 - (lam_i/mu)^2]
    (BFGS inner solve), then update lam_i <- max(0, lam_i - mu g_i(x)). With a *fixed* mu the multipliers
    converge linearly and the Hessian stays well conditioned, which is why this method, not the plain
    penalty, is used in practice. Returns (x, lam, outer iterations).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- analysis helpers


def observed_order(errors: np.ndarray) -> float:
    """Return the convergence order q estimated from three successive errors e1 > e2 > e3 of an
    iteration: q = ln(e3 / e2) / ln(e2 / e1). q = 1 is linear, 2 is quadratic (Newton)."""
    raise NotImplementedError


def observed_order_from_resolutions(errors: np.ndarray, h: np.ndarray) -> float:
    """Return the slope of log(error) against log(h) by least squares: the observed order in h."""
    raise NotImplementedError


def shadow_price_check(params: StrutParams | None = None, delta: float = 1e-3) -> tuple[float, float]:
    """Return (predicted, finite-difference) values of d(mass*)/dP in kg/N.

    Prediction: at a KKT point the sensitivity of the optimum to a parameter is that of the Lagrangian,
    d f*/dP = dL/dP = -sum_i lam_i dg_i/dP; with the normalised constraints dg_i/dP = -(1 - g_i)/P, which
    is -1/P on an active constraint, so d(mass*)/dP = (sum of the active multipliers) / P.
    The finite difference recomputes the analytic optimum at P (1 +- delta).
    """
    raise NotImplementedError


def design_chart(loads: np.ndarray, params: StrutParams | None = None) -> list[tuple[float, float, tuple[int, int]]]:
    """Return [(P, mass*, active pair)] for each load in ``loads``: the sizing chart of part 4."""
    raise NotImplementedError


def textbook_problem() -> Problem:
    """Return  min x^2 + y^2  s.t.  x + y >= 1,  x >= -1.  Solution (1/2, 1/2), lam = (1, 0)."""
    return Problem(
        f=lambda x: float(x[0] ** 2 + x[1] ** 2),
        grad=lambda x: 2.0 * np.asarray(x, dtype=float),
        cons=lambda x: np.array([x[0] + x[1] - 1.0, x[0] + 1.0]),
        jac=lambda x: np.array([[1.0, 1.0], [1.0, 0.0]]),
        n=2, m=2,
        hess_lagr=lambda x, lam: 2.0 * np.eye(2),
        names=("sum", "box"),
    )


def rosenbrock_disk_problem(radius2: float = 1.5) -> Problem:
    """Return the Rosenbrock function minimised inside the disk x^2 + y^2 <= radius2 (one constraint)."""
    def f(x: np.ndarray) -> float:
        return float((1.0 - x[0]) ** 2 + 100.0 * (x[1] - x[0] ** 2) ** 2)

    def grad(x: np.ndarray) -> np.ndarray:
        return np.array([-2.0 * (1.0 - x[0]) - 400.0 * x[0] * (x[1] - x[0] ** 2), 200.0 * (x[1] - x[0] ** 2)])

    def hess_lagr(x: np.ndarray, lam: np.ndarray) -> np.ndarray:
        Hf = np.array([[2.0 - 400.0 * (x[1] - x[0] ** 2) + 800.0 * x[0] ** 2, -400.0 * x[0]], [-400.0 * x[0], 200.0]])
        return Hf + lam[0] * 2.0 * np.eye(2)             # g = R - x^2 - y^2, so -lam * hess(g) = +2 lam I

    return Problem(f=f, grad=grad, cons=lambda x: np.array([radius2 - x[0] ** 2 - x[1] ** 2]),
                   jac=lambda x: np.array([[-2.0 * x[0], -2.0 * x[1]]]), n=2, m=1, hess_lagr=hess_lagr,
                   names=("disk",))


# ---------------------------------------------------------------- main


def main(out_dir: str | Path = ".") -> None:
    """Size the strut by SQP, check against the analytic vertex and SciPy, and plot the design space."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.optimize import minimize

    out_dir = Path(out_dir)
    params = StrutParams()
    prob = strut_problem(params)
    x0 = np.array([1.0, 2.0])                              # r = 50 mm, t = 2 mm: feasible, heavy
    print(f"strut: P = {params.P / 1e3:.0f} kN, L = {params.L} m, 6061-T6, k_local = {params.k_local:.4f}")
    print(f"start: r = {x0[0] * params.r_ref * 1e3:.1f} mm, t = {x0[1] * params.t_ref * 1e3:.2f} mm,"
          f" mass {strut_mass(x0, params):.3f} kg, g = {np.round(strut_constraints(x0, params), 3)}")

    x_star, active, m_star = analytic_optimum(params)
    r_s, t_s = unscale(x_star, params)
    lam_star = multipliers_least_squares(prob, x_star, list(active))
    print(f"\nanalytic vertex: active = {[prob.names[i] for i in active]}, r = {r_s * 1e3:.3f} mm,"
          f" t = {t_s * 1e3:.4f} mm, mass = {m_star:.5f} kg")
    print(f"  stress at optimum {params.P / (2 * np.pi * r_s * t_s) / 1e6:.1f} MPa (allowable {params.sigma_allow / 1e6:.0f}),"
          f" r/t = {r_s / t_s:.1f}, multipliers = {np.round(lam_star, 5)} kg per unit g")

    for mode in ("bfgs", "exact"):
        res = sqp(prob, x0, hessian=mode, tol=1e-10)
        err = [np.linalg.norm(h - x_star) for h in res.history]
        print(f"\nSQP ({mode} Hessian): {res.n_iter} iterations, converged = {res.converged}, mass = {res.fun:.6f} kg,"
              f" |x - x*| = {err[-1]:.2e}, KKT = {res.kkt_history[-1]:.2e}")
        print("  |x_k - x*| per iteration:", " ".join(f"{e:.1e}" for e in err))
        print("  multipliers:", np.round(res.lam, 5))

    cons = [{"type": "ineq", "fun": prob.cons, "jac": prob.jac}]
    r1 = minimize(prob.f, x0, jac=prob.grad, constraints=cons, method="SLSQP", options={"ftol": 1e-14, "maxiter": 200})
    print(f"\nscipy SLSQP: {r1.nit} iterations, mass = {r1.fun:.6f} kg, |x - x*| = {np.linalg.norm(r1.x - x_star):.2e}")
    from scipy.optimize import NonlinearConstraint
    r2 = minimize(prob.f, x0, jac=prob.grad, method="trust-constr",
                  constraints=[NonlinearConstraint(prob.cons, 0.0, np.inf, jac=prob.jac)],
                  options={"gtol": 1e-10, "xtol": 1e-12, "maxiter": 1000})
    print(f"scipy trust-constr: {r2.nit} iterations, mass = {r2.fun:.6f} kg, |x - x*| = {np.linalg.norm(r2.x - x_star):.2e}")

    pred, fd = shadow_price_check(params)
    print(f"\nshadow price d(mass*)/dP: predicted from multipliers {pred * 1e3:.5f} g/N, finite difference {fd * 1e3:.5f} g/N")

    print("\nquadratic penalty: error |x_mu - x*| against mu")
    mus = np.array([1e1, 1e2, 1e3, 1e4])
    errs = []
    for mu in mus:
        xm = penalty_minimize(prob, mu, x0)
        errs.append(np.linalg.norm(xm - x_star))
        print(f"  mu = {mu:8.0f}: |x_mu - x*| = {errs[-1]:.2e}, max violation = {-np.min(np.minimum(prob.cons(xm), 0)):.2e}")
    print(f"  observed order in 1/mu: {observed_order_from_resolutions(errs, 1.0 / mus):.3f} (theory 1)")
    x_al, lam_al, k_al = augmented_lagrangian(prob, x0, mu=10.0)
    print(f"augmented Lagrangian (mu = 10): {k_al} outer iterations, |x - x*| = {np.linalg.norm(x_al - x_star):.2e},"
          f" multipliers {np.round(lam_al, 5)}")

    print("\nsizing chart: load -> optimal mass and active constraints")
    for P, mass, act in design_chart(np.array([5e3, 10e3, 20e3, 40e3, 80e3, 160e3]), params):
        print(f"  P = {P / 1e3:6.1f} kN: mass {mass:.4f} kg, active {[prob.names[i] for i in act]}")

    # Figure: design space with constraint boundaries and the SQP path.
    res = sqp(prob, x0, hessian="bfgs", tol=1e-10)
    r_mm = np.linspace(15.0, 100.0, 300)
    t_mm = np.linspace(0.3, 2.5, 300)
    R, T = np.meshgrid(r_mm, t_mm)
    X = np.stack([R / (params.r_ref * 1e3), T / (params.t_ref * 1e3)], axis=-1)
    mass = params.rho * params.L * 2 * np.pi * (R * 1e-3) * (T * 1e-3)
    G = np.stack([strut_constraints(X[i, j], params) for i in range(R.shape[0]) for j in range(R.shape[1])])
    G = G.reshape(R.shape + (5,))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6))
    ax = axes[0]
    cs = ax.contour(R, T, mass, levels=12, colors="0.6", linewidths=0.7)
    ax.clabel(cs, fmt="%.1f kg", fontsize=7)
    feasible = np.all(G >= 0, axis=-1)
    ax.contourf(R, T, feasible, levels=[0.5, 1.5], colors=["#cfe8cf"], alpha=0.6)
    for i, (name, col) in enumerate(zip(prob.names, ["C0", "C1", "C2", "C3", "C4"])):
        ax.contour(R, T, G[..., i], levels=[0.0], colors=[col], linewidths=1.5)
        ax.plot([], [], color=col, label=f"{name} = 0")
    path = np.array(res.history)
    ax.plot(path[:, 0] * params.r_ref * 1e3, path[:, 1] * params.t_ref * 1e3, "ko-", ms=4, lw=1, label="SQP iterates")
    ax.plot(r_s * 1e3, t_s * 1e3, "r*", ms=12, label="analytic optimum")
    ax.set_xlabel("mean radius r [mm]"); ax.set_ylabel("wall thickness t [mm]")
    ax.set_title("strut design space (green = feasible)"); ax.legend(fontsize=7, loc="upper right")
    ax = axes[1]
    for mode, mk in (("bfgs", "o-"), ("exact", "s-")):
        rr = sqp(prob, x0, hessian=mode, tol=1e-12)
        e = [np.linalg.norm(h - x_star) + 1e-17 for h in rr.history]
        ax.semilogy(e, mk, ms=4, label=f"SQP, {mode} Hessian")
    ax.set_xlabel("iteration"); ax.set_ylabel("|x_k - x*| (scaled variables)")
    ax.set_title("convergence of SQP"); ax.grid(True, which="both", lw=0.3); ax.legend()
    fig.suptitle("P03.15  SQP sizing of a tubular strut")
    fig.tight_layout()
    fig.savefig(out_dir / "strut_sqp.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'strut_sqp.png'}")


if __name__ == "__main__":
    main()
