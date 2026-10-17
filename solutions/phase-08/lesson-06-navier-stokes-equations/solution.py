"""P08.06 Navier-Stokes Equations: reference solution.

How much hydraulic fluid leaks past a flight-control actuator's piston, how much viscous friction
does the piston feel, and how fast does the gap flow settle? The answer comes from the exact
solutions of the incompressible Navier-Stokes equations for a Newtonian fluid, each checked twice:

* closed forms: Couette, plane Poiseuille and Couette-Poiseuille flow between parallel plates,
  pipe (Hagen-) Poiseuille flow, flow in a concentric annulus, Stokes's first problem (a plate
  started impulsively), and the start-up of plane Poiseuille flow (a Fourier series);
* substitution: SymPy puts each closed form back into the full Navier-Stokes equations (all the
  terms, including the nonlinear ones) and the residual must simplify to zero;
* finite differences: a steady two-point boundary-value solver for the channel and the pipe, and a
  Crank-Nicolson solver for the unsteady problems, whose observed order of convergence is measured.

Coordinates: x along the flow, y across a gap of height h (wall y = 0 at rest, wall y = h moving
at speed U in +x). G = -dp/dx is the driving pressure gradient (Pa/m, positive when pressure falls
in +x). In the pipe, r is the radius and the axial pressure gradient is again G = -dp/dz.

Units: SI (m, s, Pa, Pa s, m^2/s).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes navier_stokes.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-06-navier-stokes-equations
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np
import sympy as sp
from scipy.linalg import solve_banded
from scipy.special import erfc

# Typical properties of a MIL-PRF-83282-type synthetic hydrocarbon hydraulic fluid near 40 degC
# (kinematic viscosity about 14 mm^2/s, density about 845 kg/m^3); typical values, not a data sheet.
RHO_HYD = 845.0                 # kg/m^3
NU_HYD = 14.0e-6                # m^2/s
MU_HYD = RHO_HYD * NU_HYD       # Pa s (about 0.0118)

# SymPy symbols used by the substitution checks
X, Y, T = sp.symbols("x y t", real=True)
R_SYM, Z_SYM = sp.symbols("r z", positive=True)


# ---------------------------------------------------------------- closed forms

def couette_poiseuille(y: np.ndarray, h: float, U: float, G: float, mu: float) -> np.ndarray:
    """Steady fully developed flow between plates y = 0 (at rest) and y = h (moving at U), driven by
    the pressure gradient G = -dp/dx:  u(y) = U y/h + G y (h - y) / (2 mu).
    G = 0 is Couette flow; U = 0 is plane Poiseuille flow."""
    y = np.asarray(y, dtype=float)
    return U * y / h + G * y * (h - y) / (2.0 * mu)


def plane_flow_rate(h: float, U: float, G: float, mu: float) -> float:
    """Volume flow per unit width (m^2/s) of Couette-Poiseuille flow: q = U h / 2 + G h^3 / (12 mu)."""
    return U * h / 2.0 + G * h ** 3 / (12.0 * mu)


def wall_shear_stress(h: float, U: float, G: float, mu: float) -> tuple[float, float]:
    """The shear stress tau_xy = mu du/dy at y = 0 and at y = h for Couette-Poiseuille flow:
    mu U/h + G h/2 and mu U/h - G h/2. The force per unit area the fluid exerts on the lower wall is
    +tau_xy(0) in x, on the upper wall -tau_xy(h) (the fluid lies on the other side)."""
    return mu * U / h + G * h / 2.0, mu * U / h - G * h / 2.0


def pipe_poiseuille(r: np.ndarray, R: float, G: float, mu: float) -> np.ndarray:
    """Hagen-Poiseuille flow in a round pipe of radius R: w(r) = G (R^2 - r^2) / (4 mu)."""
    r = np.asarray(r, dtype=float)
    return G * (R ** 2 - r ** 2) / (4.0 * mu)


def pipe_flow_rate(R: float, G: float, mu: float) -> float:
    """Volume flow of Hagen-Poiseuille flow: Q = pi R^4 G / (8 mu)."""
    return np.pi * R ** 4 * G / (8.0 * mu)


def annulus_flow_rate(Ri: float, Ro: float, G: float, mu: float) -> float:
    """Pressure-driven volume flow in a concentric annulus Ri < r < Ro (both walls at rest):
    Q = pi G / (8 mu) [Ro^4 - Ri^4 - (Ro^2 - Ri^2)^2 / ln(Ro/Ri)].
    For a thin gap the bracket is a small difference of O(1) terms, so it is factored with
    delta = (Ro - Ri)/Ro and k = 1 - delta as delta (1 + k) [1 + k^2 - delta (1 + k) / ln(1/k)], with
    ln(1/k) = -log1p(-delta); the result tends to the plane slot pi (Ri + Ro) (Ro - Ri)^3 G / (12 mu)."""
    delta = (Ro - Ri) / Ro
    k = Ri / Ro
    bracket = delta * (1.0 + k) * ((1.0 + k ** 2) - delta * (1.0 + k) / (-np.log1p(-delta)))
    return np.pi * G * Ro ** 4 * bracket / (8.0 * mu)


def stokes_first_problem(y: np.ndarray, t: float, U: float, nu: float) -> np.ndarray:
    """A plate at y = 0 started impulsively to speed U at t = 0 under fluid at rest (y > 0):
    u = U erfc(y / (2 sqrt(nu t))). The 99 % thickness is 3.64 sqrt(nu t)."""
    y = np.asarray(y, dtype=float)
    return U * erfc(y / (2.0 * np.sqrt(nu * t)))


def startup_channel_series(y: np.ndarray, t: float, h: float, G: float, rho: float, mu: float,
                           n_terms: int = 200) -> np.ndarray:
    """Start-up of plane Poiseuille flow: fluid at rest between fixed plates, a pressure gradient G
    switched on at t = 0. u = u_s(y) - (G/(2 mu)) sum over odd n of 8 h^2/(n pi)^3 sin(n pi y/h)
    exp(-n^2 pi^2 nu t / h^2), with u_s the steady parabola (the Fourier sine series of u_s)."""
    y = np.asarray(y, dtype=float)
    nu = mu / rho
    u = G * y * (h - y) / (2.0 * mu)
    for n in range(1, 2 * n_terms, 2):
        k = n * np.pi / h
        u = u - G / (2.0 * mu) * 8.0 * h ** 2 / (n * np.pi) ** 3 * np.sin(k * y) * np.exp(-nu * k ** 2 * t)
    return u


def startup_time(h: float, nu: float, fraction: float = 0.99) -> float:
    """Time for the centreline speed of start-up plane Poiseuille flow to reach ``fraction`` of its
    steady value, from the slowest (n = 1) term: t = h^2 / (pi^2 nu) ln(32 / (pi^3 (1 - fraction)))."""
    return h ** 2 / (np.pi ** 2 * nu) * np.log(32.0 / (np.pi ** 3 * (1.0 - fraction)))


# ---------------------------------------------------------------- substitution into the equations

def ns_residuals(u: sp.Expr, v: sp.Expr, p: sp.Expr, rho: sp.Expr, mu: sp.Expr,
                 gx: sp.Expr = 0, gy: sp.Expr = 0) -> tuple[sp.Expr, sp.Expr, sp.Expr]:
    """Residuals of the 2-D incompressible Navier-Stokes equations for a Newtonian fluid with constant
    rho and mu, for velocity (u, v) and pressure p given as SymPy expressions in X, Y, T:

        continuity:  u_x + v_y
        x-momentum:  rho (u_t + u u_x + v u_y) + p_x - mu (u_xx + u_yy) - rho gx
        y-momentum:  rho (v_t + u v_x + v v_y) + p_y - mu (v_xx + v_yy) - rho gy

    Each is returned simplified; an exact solution gives (0, 0, 0)."""
    d = sp.diff
    cont = d(u, X) + d(v, Y)
    xm = rho * (d(u, T) + u * d(u, X) + v * d(u, Y)) + d(p, X) - mu * (d(u, X, 2) + d(u, Y, 2)) - rho * gx
    ym = rho * (d(v, T) + u * d(v, X) + v * d(v, Y)) + d(p, Y) - mu * (d(v, X, 2) + d(v, Y, 2)) - rho * gy
    return sp.simplify(cont), sp.simplify(xm), sp.simplify(ym)


def pipe_axial_residual(w: sp.Expr, p: sp.Expr, mu: sp.Expr) -> sp.Expr:
    """Residual of the axial momentum equation for steady, fully developed, axisymmetric pipe flow
    w(r) (continuity and the convective terms vanish identically for w = w(r), v_r = v_theta = 0):
    mu (1/r) d/dr (r dw/dr) - dp/dz, in the symbols R_SYM and Z_SYM."""
    return sp.simplify(mu / R_SYM * sp.diff(R_SYM * sp.diff(w, R_SYM), R_SYM) - sp.diff(p, Z_SYM))


# ---------------------------------------------------------------- finite differences

def _tridiag_solve(lower: np.ndarray, diag: np.ndarray, upper: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """Solve a tridiagonal system (lower[i] multiplies x[i-1], upper[i] multiplies x[i+1])."""
    n = diag.size
    ab = np.zeros((3, n))
    ab[0, 1:] = upper[:-1]
    ab[1] = diag
    ab[2, :-1] = lower[1:]
    return solve_banded((1, 1), ab, rhs)


def fd_steady_channel(n: int, h: float, U: float, f: float | Callable[[np.ndarray], np.ndarray],
                      mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Finite-difference solution of mu u'' = -f(y) on 0 < y < h with u(0) = 0, u(h) = U, on n equal
    intervals, by the central difference (u_{i-1} - 2 u_i + u_{i+1}) / dy^2 (truncation error
    dy^2 u''''/12). ``f`` is the pressure gradient G (a constant) or any body-force density f(y).
    Returns (y, u) at the n + 1 nodes."""
    y = np.linspace(0.0, h, n + 1)
    dy = h / n
    fy = f(y) if callable(f) else np.full(n + 1, float(f))
    m = n - 1
    rhs = -fy[1:-1] * dy ** 2 / mu
    rhs[-1] -= U
    u_int = _tridiag_solve(np.ones(m), -2.0 * np.ones(m), np.ones(m), rhs)
    return y, np.concatenate([[0.0], u_int, [U]])


def fd_steady_pipe(n: int, R: float, G: float, mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Finite-volume-style solution of mu (1/r)(r w')' = -G on 0 <= r < R with w(R) = 0 and symmetry
    on the axis, on n equal intervals. Interior rows use face radii r_{i +- 1/2}; the axis row uses
    (1/r)(r w')' -> 2 w'' with the mirror node w_{-1} = w_1. Returns (r, w) at the n + 1 nodes."""
    r = np.linspace(0.0, R, n + 1)
    dr = R / n
    N = n                                   # unknowns w_0 .. w_{n-1}; w_n = 0
    lower, diag, upper = np.zeros(N), np.zeros(N), np.zeros(N)
    rhs = np.full(N, -G / mu * dr ** 2)
    diag[0], upper[0] = -4.0, 4.0
    for i in range(1, N):
        rp, rm = r[i] + dr / 2, r[i] - dr / 2
        lower[i], diag[i], upper[i] = rm / r[i], -(rp + rm) / r[i], rp / r[i]
    w = _tridiag_solve(lower, diag, upper, rhs)
    return r, np.concatenate([w, [0.0]])


def crank_nicolson(y: np.ndarray, u0: np.ndarray, t0: float, t1: float, nt: int, nu: float,
                   source: float = 0.0, u_lo: Callable[[float], float] | float = 0.0,
                   u_hi: Callable[[float], float] | float = 0.0) -> np.ndarray:
    """March u_t = nu u_yy + source from t0 to t1 in nt equal steps by Crank-Nicolson (second order
    in time and, with the central difference, in space), with Dirichlet values u_lo(t), u_hi(t) at the
    ends of the uniform grid y. ``source`` is a constant acceleration (G/rho for a pressure gradient).
    Returns u at t1."""
    y = np.asarray(y, dtype=float)
    u = np.array(u0, dtype=float)
    dy = y[1] - y[0]
    dt = (t1 - t0) / nt
    lam = nu * dt / dy ** 2
    m = y.size - 2
    lo = u_lo if callable(u_lo) else (lambda _t, c=float(u_lo): c)
    hi = u_hi if callable(u_hi) else (lambda _t, c=float(u_hi): c)
    a = -0.5 * lam * np.ones(m)
    d = (1.0 + lam) * np.ones(m)
    t = t0
    for _ in range(nt):
        tn = t + dt
        u[0], u[-1] = lo(t), hi(t)
        rhs = u[1:-1] + 0.5 * lam * (u[:-2] - 2 * u[1:-1] + u[2:]) + dt * source
        b0, b1 = lo(tn), hi(tn)
        rhs[0] += 0.5 * lam * b0
        rhs[-1] += 0.5 * lam * b1
        u[1:-1] = _tridiag_solve(a, d, a, rhs)
        u[0], u[-1] = b0, b1
        t = tn
    return u


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) of errors on grids refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- the actuator

def actuator_leakage(D: float, h: float, L_land: float, dp: float, U: float, mu: float,
                     rho: float) -> dict[str, float]:
    """Leakage and friction of a piston of diameter D sliding at U (positive in the direction of the
    leak, from high to low pressure) in a bore with radial clearance h, across a land of length
    L_land holding a pressure difference dp. The gap is treated as a plane slot of width pi D (valid
    for h << D). Returns flow rate Q (m^3/s), mean gap speed, gap Reynolds number rho V h/mu, the
    shear force on the piston (N, in x), its pressure-driven and wall-driven parts, and the start-up
    time to 99 %."""
    G = dp / L_land
    width = np.pi * D
    q = plane_flow_rate(h, U, G, mu)
    _, tau_h = wall_shear_stress(h, U, G, mu)
    force = -tau_h * width * L_land            # fluid on the moving (upper) wall
    return {
        "G": G,
        "Q": q * width,
        "V_mean": q / h,
        "Re_gap": rho * (q / h) * h / mu,
        "force": force,
        "force_pressure": G * h / 2.0 * width * L_land,
        "force_couette": -mu * U / h * width * L_land,
        "t99": startup_time(h, mu / rho),
    }


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save navier_stokes.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    D, h, L_land, dp, U = 25e-3, 10e-6, 12e-3, 20.7e6, 0.1
    mu, rho = MU_HYD, RHO_HYD
    a = actuator_leakage(D, h, L_land, dp, U, mu, rho)
    a0 = actuator_leakage(D, h, L_land, dp, 0.0, mu, rho)
    print(f"hydraulic fluid (typical): rho = {rho:.0f} kg/m^3, mu = {mu:.4g} Pa s, nu = {NU_HYD * 1e6:.1f} mm^2/s")
    print(f"actuator: D = {D * 1e3:.0f} mm, clearance h = {h * 1e6:.0f} um, land {L_land * 1e3:.0f} mm, "
          f"dp = {dp / 1e6:.1f} MPa (G = {a['G']:.4g} Pa/m), piston speed {U} m/s")
    print(f"  leakage, piston at rest:   Q = {a0['Q']:.4e} m^3/s = {a0['Q'] * 6e4:.4f} L/min")
    print(f"  leakage, piston moving:    Q = {a['Q']:.4e} m^3/s = {a['Q'] * 6e4:.4f} L/min "
          f"(Couette part {U * h / 2 * np.pi * D:.3e} m^3/s)")
    print(f"  mean gap speed {a['V_mean']:.3f} m/s, gap Reynolds number {a['Re_gap']:.3f}")
    print(f"  shear force on piston {a['force']:.3f} N = pressure-driven {a['force_pressure']:.3f} N "
          f"+ wall-driven {a['force_couette']:.4f} N")
    print(f"  pressure force on the land annulus dp * pi D h = {dp * np.pi * D * h:.3f} N")
    print(f"  start-up to 99 %: {a['t99'] * 1e6:.3f} us (h^2/nu = {h ** 2 / NU_HYD * 1e6:.3f} us)")
    Ri, Ro = D / 2, D / 2 + h
    Qa = annulus_flow_rate(Ri, Ro, a['G'], mu)
    Qs = np.pi * (Ri + Ro) * h ** 3 * a['G'] / (12 * mu)
    print(f"  concentric annulus Q = {Qa:.5e} m^3/s, plane slot pi (Ri+Ro) h^3 G/(12 mu) = {Qs:.5e}: "
          f"|ratio - 1| < 1e-7 is {abs(Qa / Qs - 1) < 1e-7}")
    for hh in (0.1e-3, 1e-3):
        r_ = annulus_flow_rate(Ri, Ri + hh, 1.0, 1.0) / (np.pi * (2 * Ri + hh) * hh ** 3 / 12)
        print(f"  annulus / slot for a {hh * 1e3:.1f} mm gap on a {D * 1e3:.0f} mm piston: {r_:.6f}")
    for hh in (5e-6, 10e-6, 20e-6):
        print(f"  clearance {hh * 1e6:4.0f} um: leakage {actuator_leakage(D, hh, L_land, dp, 0.0, mu, rho)['Q'] * 6e4:.4f} L/min")

    # verification: finite differences
    print("\nfinite differences, steady Couette-Poiseuille (max |u_FD - u_exact| / max u):")
    for n in (4, 8, 16):
        yy, uu = fd_steady_channel(n, h, U, a['G'], mu)
        ex = couette_poiseuille(yy, h, U, a['G'], mu)
        print(f"  n = {n:3d}: {np.max(np.abs(uu - ex)) / np.max(ex):.2e}  (exact: u'''' = 0)")
    hm = 1.0
    f_m = lambda yy: np.pi ** 2 * np.sin(np.pi * yy / hm)          # manufactured: u = sin(pi y)
    errs = []
    for n in (10, 20, 40, 80):
        yy, uu = fd_steady_channel(n, hm, 0.0, f_m, 1.0)
        errs.append(np.max(np.abs(uu - np.sin(np.pi * yy))))
    print(f"  manufactured u = sin(pi y): errors {np.array2string(np.array(errs), precision=3, floatmode='fixed', formatter={'float': lambda v: f'{v:.3e}'})}, "
          f"orders {np.array2string(observed_order(errs), precision=3)}")
    errs_s = []
    nu1, t0, t1, Yd = 1.0, 0.01, 0.05, 2.0
    for n in (20, 40, 80, 160):
        yy = np.linspace(0, Yd, n + 1)
        u = crank_nicolson(yy, stokes_first_problem(yy, t0, 1.0, nu1), t0, t1, n // 2, nu1,
                           u_lo=1.0, u_hi=lambda tt: float(stokes_first_problem(Yd, tt, 1.0, nu1)))
        errs_s.append(np.max(np.abs(u - stokes_first_problem(yy, t1, 1.0, nu1))))
    print(f"  Stokes's first problem, Crank-Nicolson (dt ~ dy): errors {np.array2string(np.array(errs_s), precision=3)}, "
          f"orders {np.array2string(observed_order(errs_s), precision=3)}")
    errs_u = []
    for n in (20, 40, 80, 160):
        yy = np.linspace(0, 1.0, n + 1)
        t_end = 0.05
        u = crank_nicolson(yy, np.zeros_like(yy), 0.0, t_end, n, 1.0, source=1.0)
        errs_u.append(np.max(np.abs(u - startup_channel_series(yy, t_end, 1.0, 1.0, 1.0, 1.0))))
    print(f"  start-up Poiseuille from rest: errors {np.array2string(np.array(errs_u), precision=3)}, "
          f"orders {np.array2string(observed_order(errs_u), precision=3)}")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    yy = np.linspace(0, h, 101)
    for UU, lab in ((0.0, "piston at rest (Poiseuille)"), (U, f"piston at {U} m/s (Couette-Poiseuille)"),
                    (-U, f"piston at {-U} m/s")):
        ax[0].plot(couette_poiseuille(yy, h, UU, a['G'], mu), yy * 1e6, label=lab)
    yf, uf = fd_steady_channel(8, h, U, a['G'], mu)
    ax[0].plot(uf, yf * 1e6, "ko", ms=4, label="finite differences, n = 8")
    ax[0].set_xlabel("u (m/s)")
    ax[0].set_ylabel("y (um)")
    ax[0].set_title("actuator clearance flow")
    ax[0].legend(fontsize=8)
    eta = np.linspace(0, 4, 200)
    for tt in (0.25, 1.0, 4.0):
        ax[1].plot(stokes_first_problem(eta, tt, 1.0, 1.0), eta, label=f"t nu / l^2 = {tt}")
    ax[1].set_xlabel("u / U")
    ax[1].set_ylabel("y / l")
    ax[1].set_title("Stokes's first problem")
    ax[1].legend(fontsize=8)
    ns_ = np.array([20, 40, 80, 160])
    ax[2].loglog(1 / ns_, errs_s, "o-", label="Stokes's first problem (CN)")
    ax[2].loglog(1 / ns_, errs_u, "s-", label="start-up Poiseuille (CN)")
    ax[2].loglog(1 / np.array([10, 20, 40, 80]), errs, "^-", label="manufactured steady")
    ax[2].loglog(1 / ns_, errs_s[0] * (ns_[0] / ns_) ** 2, "k--", lw=0.8, label="slope 2")
    ax[2].set_xlabel("dy / l")
    ax[2].set_ylabel("max error")
    ax[2].set_title("finite-difference convergence")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "navier_stokes.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'navier_stokes.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
