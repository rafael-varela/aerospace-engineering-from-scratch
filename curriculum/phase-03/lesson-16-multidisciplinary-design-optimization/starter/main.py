"""P03.16 Multidisciplinary Design Optimization: starter code.

A two-discipline wing model (aerodynamics + mission, structures) coupled through the root bending
moment and the wing weight, solved as a fixed point (Gauss-Seidel and Newton multidisciplinary
analysis), differentiated by finite differences, by the direct (forward) coupled method and by the
coupled adjoint, and optimised for minimum block fuel with ``scipy.optimize.minimize`` (SLSQP) as the
driver in two architectures: MDF (the analysis converged inside every function evaluation) and IDF
(the coupling variables become design variables with consistency constraints).

Everything is NumPy/SciPy. OpenMDAO, which automates exactly this bookkeeping, is an optional
extension of the lesson and is never imported here.

Conventions:
    * Design variables x = (b, S): span [m] and reference area [m^2].
    * Coupling variables y = (y1, y2) = (M_root [N m], W_wing [N]).  Discipline 1 (aero + mission)
      computes y1 from (x, y2); discipline 2 (structures) computes y2 from (x, y1).
    * Coupled residual R(x, y) = y - G(x, y) with G = (G1(x, y2), G2(x, y1)); the MDA solves R = 0.
    * Functions of interest: block fuel W_fuel [N] (objective) and the stall constraint g_stall;
      the deflection and gate constraints depend on x only.
    * Weights are forces in newtons (W = m g).

Given: the dataclasses, ``design_constraints_x_only``, ``mdf_optimize``, ``idf_optimize`` and ``main``;
everything else is yours.

Run:  python main.py          (prints the MDA, derivative and optimisation tables, writes wing_mdo.png)
Test: pytest curriculum/phase-03/lesson-16-multidisciplinary-design-optimization
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

G0 = 9.80665  # standard gravity [m/s^2]


# ---------------------------------------------------------------- parameters (GIVEN)


@dataclass(frozen=True)
class WingParams:
    """Fixed parameters of the wing aero-structural model. SI units; weights in newtons.

    The aircraft is a light twin-engine utility aeroplane of the author's invention: realistic in
    magnitude, not a specific type. Nothing here was measured.
    """

    W0: float = 25.0e3          # weight of everything except wing structure and fuel [N]
    rho: float = 0.909          # cruise air density, ISA 3000 m [kg/m^3]
    V: float = 80.0             # cruise true airspeed [m/s]
    f_rest: float = 0.25        # equivalent flat-plate drag area of fuselage, tail, nacelles [m^2]
    c_dw: float = 0.008         # wing profile drag coefficient (per unit S)
    e: float = 0.85             # Oswald span-efficiency factor
    R: float = 1.5e6            # design range [m]
    c_T: float = 1.4e-5         # thrust-specific fuel consumption [kg/(N s)]
    tau: float = 0.12           # spar depth / chord
    sigma_a: float = 120.0e6    # allowable (fatigue) stress in the spar caps [Pa]
    rho_m: float = 2700.0       # cap material density [kg/m^3]
    E: float = 70.0e9           # cap material Young's modulus [Pa]
    n_ult: float = 3.75         # ultimate load factor (2.5 g limit x 1.5 safety)
    w_skin: float = 60.0        # areal weight of skins, ribs, controls, systems [N/m^2]
    q_stall: float = 0.5 * 1.225 * 45.0 ** 2   # dynamic pressure at the 45 m/s sea-level stall speed [Pa]
    CL_max: float = 1.8         # maximum lift coefficient with flaps
    defl_max: float = 0.06      # allowed tip deflection / semi-span at limit load
    b_max: float = 20.0         # hangar/gate limit on span [m]

    @property
    def q(self) -> float:
        """Cruise dynamic pressure [Pa]."""
        return 0.5 * self.rho * self.V ** 2

    @property
    def kappa(self) -> float:
        """Breguet exponent factor R g c_T / V (dimensionless): W_fuel = W_zf (exp(kappa D/L) - 1)."""
        return self.R * G0 * self.c_T / self.V


@dataclass
class MDAResult:
    """What an MDA solver returns: the converged coupling vector, the residual history and the count."""

    y: np.ndarray
    residuals: list[float]
    n_iter: int
    converged: bool


# ---------------------------------------------------------------- the two disciplines


def aero_mission(x: np.ndarray, W_wing: float, params: WingParams) -> dict[str, float]:
    """Discipline 1: aerodynamics and mission. Return a dict with keys
    'L' (cruise lift = zero-fuel weight) [N], 'D' (cruise drag) [N], 'W_fuel' [N], 'W_TO' [N],
    'M_root' (ultimate root bending moment) [N m], and 'LD' (lift-to-drag ratio).

    Lift equals the zero-fuel weight W_zf = W0 + W_wing (the lift-to-drag ratio is evaluated at the
    end-of-cruise weight, a stated simplification). Drag is the parabolic polar in dimensional form,
    D = q (f_rest + c_dw S) + L^2 / (q pi e b^2). Block fuel follows Breguet,
    W_fuel = W_zf (exp(kappa D / L) - 1), take-off weight is W_zf + W_fuel, and the ultimate root
    bending moment of an elliptically loaded wing is M_root = n_ult W_TO b / (3 pi).
    """
    raise NotImplementedError


def aero_mission_partials(x: np.ndarray, W_wing: float, params: WingParams) -> dict[str, np.ndarray]:
    """Return the analytic partial derivatives of discipline 1 with respect to (b, S, W_wing), as a dict
    mapping 'M_root', 'W_fuel', 'W_TO' to length-3 arrays [d/db, d/dS, d/dW_wing].

    Derived by hand from the formulas in ``aero_mission`` with phi = D / L:
    dphi/db = -2 L / (q pi e b^3),  dphi/dS = q c_dw / L,  dphi/dL = -q (f_rest + c_dw S) / L^2 + 1 / (q pi e b^2),
    W_TO = L exp(kappa phi),  W_fuel = W_TO - L,  M_root = n b W_TO / (3 pi).
    """
    raise NotImplementedError


def structures(x: np.ndarray, M_root: float, params: WingParams) -> dict[str, float]:
    """Discipline 2: structures. Return a dict with 'A_cap' (spar cap area at the root) [m^2],
    'W_caps' [N], 'W_wing' [N] and 'defl_ratio' (tip deflection / semi-span at limit load).

    The spar box has depth h = tau c with c = S / b. Two caps of area A carry the bending moment at
    the allowable stress: A = M_root / (h sigma_a). Sizing every section to its own moment, with
    M(y) = M_root (1 - 2|y|/b)^2, the cap weight is rho_m g (2/3) A b. Skins, ribs and systems add
    w_skin S. A fully stressed beam has constant curvature 2 sigma_a / (E h) at ultimate load, so the
    limit-load tip deflection ratio is sigma_a b / (3 E h) = sigma_a b^2 / (3 E tau S): it depends on the
    geometry only, not on the load.
    """
    raise NotImplementedError


def structures_partials(x: np.ndarray, M_root: float, params: WingParams) -> np.ndarray:
    """Return the analytic partials of W_wing with respect to (b, S, M_root) as a length-3 array.

    W_wing = (2/3) rho_m g b^2 M_root / (tau S sigma_a) + w_skin S.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the coupled system


def coupled_residual(x: np.ndarray, y: np.ndarray, params: WingParams) -> np.ndarray:
    """Return R(x, y) = y - G(x, y) with G1 = M_root(x, y2) from aero and G2 = W_wing(x, y1) from structures."""
    raise NotImplementedError


def coupled_jacobians(x: np.ndarray, y: np.ndarray, params: WingParams) -> tuple[np.ndarray, np.ndarray]:
    """Return (dR/dy, dR/dx): the 2 x 2 and 2 x 2 Jacobians of the coupled residual at (x, y).

    dR/dy = I - dG/dy, where dG/dy has the off-diagonal coupling terms dM_root/dW_wing and
    dW_wing/dM_root; dR/dx = -dG/dx.
    """
    raise NotImplementedError


def mda_gauss_seidel(x: np.ndarray, params: WingParams, y0: np.ndarray | None = None, tol: float = 1e-10,
                     max_iter: int = 200) -> MDAResult:
    """Solve R(x, y) = 0 by block Gauss-Seidel: y1 <- G1(x, y2), then y2 <- G2(x, y1), repeated.

    The residual norm |R| is recorded after every sweep. Converges linearly with rate equal to the
    coupling strength k = (dM_root/dW_wing)(dW_wing/dM_root) when |k| < 1 (the fixed point is a
    contraction), and diverges otherwise.
    """
    raise NotImplementedError


def mda_newton(x: np.ndarray, params: WingParams, y0: np.ndarray | None = None, tol: float = 1e-12,
               max_iter: int = 50) -> MDAResult:
    """Solve R(x, y) = 0 by Newton's method with the analytic dR/dy: quadratic convergence."""
    raise NotImplementedError


def coupling_strength(x: np.ndarray, y: np.ndarray, params: WingParams) -> float:
    """Return k = (dM_root/dW_wing)(dW_wing/dM_root) at (x, y): the Gauss-Seidel contraction factor."""
    raise NotImplementedError


# ---------------------------------------------------------------- functions of interest and their derivatives


def functions_of_interest(x: np.ndarray, y: np.ndarray, params: WingParams) -> np.ndarray:
    """Return F(x, y) = [W_fuel, g_stall] at a (not necessarily converged) coupling state y.

    g_stall = 1 - W_TO / (q_stall CL_max S) >= 0 says the wing is large enough to fly at the stall
    speed at take-off weight.
    """
    raise NotImplementedError


def functions_partials(x: np.ndarray, y: np.ndarray, params: WingParams) -> tuple[np.ndarray, np.ndarray]:
    """Return (dF/dx, dF/dy), each 2 x 2, the *partial* derivatives of F holding the other arguments fixed."""
    raise NotImplementedError


def total_derivatives_direct(x: np.ndarray, params: WingParams) -> tuple[np.ndarray, np.ndarray]:
    """Return (F, dF/dx) at the converged MDA by the direct (forward) coupled method.

    Solve (dR/dy) (dy/dx) = -(dR/dx) once per design variable (here: one 2 x 2 solve with two
    right-hand sides), then dF/dx = dF/dx|_partial + (dF/dy)(dy/dx).
    """
    raise NotImplementedError


def total_derivatives_adjoint(x: np.ndarray, params: WingParams) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (F, dF/dx, psi) at the converged MDA by the coupled adjoint method.

    For each function F_i solve (dR/dy)^T psi_i = (dF_i/dy)^T, then dF_i/dx = dF_i/dx|_partial - psi_i^T dR/dx.
    One linear solve per function instead of one per design variable: the method of choice when there
    are many design variables and few functions. psi has one row per function.
    """
    raise NotImplementedError


def total_derivatives_fd(x: np.ndarray, params: WingParams, h_rel: float = 1e-6, mda_tol: float = 1e-12,
                         mda_solver: str = "newton") -> np.ndarray:
    """Return dF/dx by central differences of the converged MDA, with relative step h_rel in each
    design variable and the MDA converged to ``mda_tol`` by ``mda_solver`` ('newton' or 'gs').
    Costs two MDA solves per design variable and carries an O(h^2) truncation error plus an error of
    order (MDA error in y) x |dF/dy| / h: a loosely converged loop poisons small steps."""
    raise NotImplementedError


def fd_convergence_study(x: np.ndarray, params: WingParams, h_rels: np.ndarray, mda_tol: float = 1e-12,
                         mda_solver: str = "newton") -> np.ndarray:
    """Return the error |dW_fuel/db (FD) - dW_fuel/db (adjoint)| for each relative step in ``h_rels``,
    with the MDA converged to ``mda_tol`` by ``mda_solver`` inside the differences: the observed order
    should be 2 until the floor set by the MDA error (or rounding) takes over."""
    raise NotImplementedError


# ---------------------------------------------------------------- optimisation


def design_constraints_x_only(x: np.ndarray, params: WingParams) -> np.ndarray:
    """Return the constraints that depend on x alone: [g_defl, g_gate] = [1 - defl_ratio/defl_max, 1 - b/b_max]."""
    b, S = float(x[0]), float(x[1])
    return np.array([1.0 - params.sigma_a * b ** 2 / (3.0 * params.E * params.tau * S) / params.defl_max,
                     1.0 - b / params.b_max])


def design_constraints_x_only_jac(x: np.ndarray, params: WingParams) -> np.ndarray:
    """Return the 2 x 2 Jacobian of ``design_constraints_x_only``."""
    raise NotImplementedError


def mdf_optimize(params: WingParams, x0: np.ndarray, use_adjoint: bool = True) -> tuple[np.ndarray, float, dict]:
    """Minimise block fuel over (b, S) in the MDF architecture: every function evaluation converges the
    MDA (Newton) first. Gradients come from the coupled adjoint when ``use_adjoint`` else SLSQP's own
    finite differences. Returns (x*, W_fuel*, info) with info holding the SciPy result and counts."""
    from scipy.optimize import minimize

    counts = {"mda": 0}
    scale = 1.0e-3                                                   # objective in kN for conditioning

    def analysis(x):
        counts["mda"] += 1
        if use_adjoint:
            F, dF, _ = total_derivatives_adjoint(x, params)
        else:
            y = mda_newton(x, params).y
            F, dF = functions_of_interest(x, y, params), None
        return F, dF

    def obj(x):
        return analysis(x)[0][0] * scale

    def obj_jac(x):
        return analysis(x)[1][0] * scale

    def cons(x):
        return np.concatenate([[analysis(x)[0][1]], design_constraints_x_only(x, params)])

    def cons_jac(x):
        return np.vstack([analysis(x)[1][1], design_constraints_x_only_jac(x, params)])

    constraints = {"type": "ineq", "fun": cons}
    if use_adjoint:
        constraints["jac"] = cons_jac
    res = minimize(obj, np.asarray(x0, dtype=float), jac=obj_jac if use_adjoint else None, method="SLSQP",
                   constraints=[constraints], bounds=[(5.0, 30.0), (5.0, 60.0)],
                   options={"ftol": 1e-12, "maxiter": 200})
    return res.x, float(res.fun / scale), {"result": res, "n_mda": counts["mda"], "n_iter": res.nit}


def idf_optimize(params: WingParams, x0: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, dict]:
    """Minimise block fuel in the IDF architecture: z = (b, S, M_root/1e5, W_wing/1e3) are all design
    variables, the disciplines are evaluated once each per function call (no MDA), and two equality
    constraints R(x, y) = 0 enforce consistency at the optimum. Returns (x*, y*, W_fuel*, info)."""
    from scipy.optimize import minimize

    ys = np.array([1.0e5, 1.0e3])                                    # scaling of the coupling variables
    scale = 1.0e-3

    def split(z):
        return z[:2], z[2:] * ys

    def obj(z):
        x, y = split(z)
        return functions_of_interest(x, y, params)[0] * scale

    def obj_jac(z):
        x, y = split(z)
        dFdx, dFdy = functions_partials(x, y, params)
        return np.concatenate([dFdx[0], dFdy[0] * ys]) * scale

    def ineq(z):
        x, y = split(z)
        return np.concatenate([[functions_of_interest(x, y, params)[1]], design_constraints_x_only(x, params)])

    def ineq_jac(z):
        x, y = split(z)
        dFdx, dFdy = functions_partials(x, y, params)
        row0 = np.concatenate([dFdx[1], dFdy[1] * ys])
        rows = np.hstack([design_constraints_x_only_jac(x, params), np.zeros((2, 2))])
        return np.vstack([row0, rows])

    def eq(z):
        x, y = split(z)
        return coupled_residual(x, y, params) / ys

    def eq_jac(z):
        x, y = split(z)
        dRdy, dRdx = coupled_jacobians(x, y, params)
        return np.hstack([dRdx, dRdy * ys[None, :]]) / ys[:, None]

    y0 = mda_newton(np.asarray(x0, dtype=float), params).y
    z0 = np.concatenate([x0, y0 / ys])
    res = minimize(obj, z0, jac=obj_jac, method="SLSQP",
                   constraints=[{"type": "ineq", "fun": ineq, "jac": ineq_jac}, {"type": "eq", "fun": eq, "jac": eq_jac}],
                   bounds=[(5.0, 30.0), (5.0, 60.0), (0.1, 100.0), (0.1, 100.0)],
                   options={"ftol": 1e-12, "maxiter": 300})
    x, y = split(res.x)
    return x, y, float(res.fun / scale), {"result": res, "n_iter": res.nit}


# ---------------------------------------------------------------- main


def main(out_dir: str | Path = ".") -> None:
    """Run the MDA, compare the three derivative routes, optimise in MDF and IDF, and plot the trade."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = WingParams()
    x0 = np.array([12.0, 16.8])
    print(f"baseline: b = {x0[0]} m, S = {x0[1]} m^2, AR = {x0[0] ** 2 / x0[1]:.2f}, q = {p.q:.0f} Pa, kappa = {p.kappa:.4f}")

    gs = mda_gauss_seidel(x0, p)
    nw = mda_newton(x0, p)
    y = nw.y
    a = aero_mission(x0, y[1], p)
    s = structures(x0, y[0], p)
    print(f"\nMDA: Gauss-Seidel {gs.n_iter} sweeps, Newton {nw.n_iter} iterations, |y_GS - y_N| = {np.linalg.norm(gs.y - y):.2e}")
    print("  GS residuals:", " ".join(f"{r:.2e}" for r in gs.residuals))
    print("  Newton residuals:", " ".join(f"{r:.2e}" for r in nw.residuals))
    print(f"  coupling strength k = {coupling_strength(x0, y, p):.4f}; observed GS ratio {gs.residuals[-2] / gs.residuals[-3]:.4f}")
    print(f"  converged: M_root = {y[0] / 1e3:.2f} kN m, W_wing = {y[1]:.1f} N ({y[1] / G0:.1f} kg; caps {s['W_caps']:.1f} N,"
          f" cap area {s['A_cap'] * 1e4:.1f} cm^2)")
    print(f"  W_zf = {a['L']:.0f} N, L/D = {a['LD']:.2f}, D = {a['D']:.0f} N, W_fuel = {a['W_fuel']:.0f} N"
          f" ({a['W_fuel'] / G0:.0f} kg), W_TO = {a['W_TO']:.0f} N, defl ratio {s['defl_ratio']:.4f}")

    F, dF_dir = total_derivatives_direct(x0, p)
    _, dF_adj, psi = total_derivatives_adjoint(x0, p)
    dF_fd = total_derivatives_fd(x0, p)
    print("\ntotal derivatives dW_fuel/d(b, S) [N/m, N/m^2] and dg_stall/d(b, S):")
    print("  direct :", np.round(dF_dir, 6).tolist())
    print("  adjoint:", np.round(dF_adj, 6).tolist())
    print("  FD     :", np.round(dF_fd, 6).tolist())
    print(f"  adjoint variables psi (fuel row): {psi[0]}")
    dFdx_partial, _ = functions_partials(x0, y, p)
    print(f"  partial dW_fuel/db (coupling frozen) = {dFdx_partial[0, 0]:.3f} N/m against total {dF_adj[0, 0]:.3f} N/m")
    hs = np.array([1e-1, 1e-2, 1e-3, 1e-4])
    errs = fd_convergence_study(x0, p, hs)
    order = np.polyfit(np.log(hs[:3]), np.log(errs[:3]), 1)[0]
    print("  FD error vs relative step:", " ".join(f"{h:.0e}:{e:.2e}" for h, e in zip(hs, errs)), f" order {order:.2f}")

    x_mdf, f_mdf, info_mdf = mdf_optimize(p, x0)
    x_mdf_fd, f_mdf_fd, info_fd = mdf_optimize(p, x0, use_adjoint=False)
    x_idf, y_idf, f_idf, info_idf = idf_optimize(p, x0)
    y_mdf = mda_newton(x_mdf, p).y
    print(f"\nMDF (adjoint gradients): b = {x_mdf[0]:.3f} m, S = {x_mdf[1]:.3f} m^2, AR = {x_mdf[0] ** 2 / x_mdf[1]:.2f},"
          f" W_fuel = {f_mdf:.1f} N ({f_mdf / G0:.1f} kg), {info_mdf['n_iter']} SLSQP iterations, {info_mdf['n_mda']} MDA solves")
    print(f"MDF (FD gradients)     : b = {x_mdf_fd[0]:.3f} m, S = {x_mdf_fd[1]:.3f} m^2, W_fuel = {f_mdf_fd:.1f} N,"
          f" {info_fd['n_iter']} iterations, {info_fd['n_mda']} MDA solves")
    print(f"IDF                    : b = {x_idf[0]:.3f} m, S = {x_idf[1]:.3f} m^2, W_fuel = {f_idf:.1f} N, {info_idf['n_iter']} iterations,"
          f" consistency residual {np.linalg.norm(coupled_residual(x_idf, y_idf, p)):.2e}")
    g = np.concatenate([[functions_of_interest(x_mdf, y_mdf, p)[1]], design_constraints_x_only(x_mdf, p)])
    print(f"  constraints at the MDF optimum (stall, deflection, gate): {np.round(g, 5)}")
    a_opt = aero_mission(x_mdf, y_mdf[1], p)
    print(f"  optimum: W_wing = {y_mdf[1]:.0f} N, L/D = {a_opt['LD']:.2f}, fuel saving vs baseline {(a['W_fuel'] - f_mdf) / G0:.1f} kg"
          f" ({100 * (a['W_fuel'] - f_mdf) / a['W_fuel']:.1f} %)")

    # Figure: fuel against span at the stall-limited area, and the MDA convergence.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    bs = np.linspace(8.0, 20.0, 49)
    for S_fixed, ls in ((14.0, "--"), (x_mdf[1], "-"), (22.0, ":")):
        fuel = [aero_mission(np.array([b, S_fixed]), mda_newton(np.array([b, S_fixed]), p).y[1], p)["W_fuel"] / G0 for b in bs]
        ax.plot(bs, fuel, ls, label=f"S = {S_fixed:.1f} m$^2$")
    ax.plot(x_mdf[0], f_mdf / G0, "r*", ms=12, label="MDF optimum")
    ax.plot(x0[0], a["W_fuel"] / G0, "ko", label="baseline")
    ax.set_xlabel("span b [m]"); ax.set_ylabel("block fuel [kg]"); ax.set_title("span against weight: the coupled trade")
    ax.legend(fontsize=8); ax.grid(lw=0.3)
    ax = axes[1]
    ax.semilogy(range(1, len(gs.residuals) + 1), gs.residuals, "o-", label="Gauss-Seidel")
    ax.semilogy(range(len(nw.residuals)), np.array(nw.residuals) + 1e-16, "s-", label="Newton")
    ax.set_xlabel("iteration"); ax.set_ylabel("|R(x, y)|"); ax.set_title("MDA convergence at the baseline")
    ax.legend(); ax.grid(lw=0.3, which="both")
    fig.suptitle("P03.16  Wing aero-structural MDO (NumPy/SciPy)")
    fig.tight_layout()
    fig.savefig(out_dir / "wing_mdo.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'wing_mdo.png'}")


if __name__ == "__main__":
    main()
