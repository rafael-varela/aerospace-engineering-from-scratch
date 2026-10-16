"""P03.06 ODE Integrators II: Adaptive & Stiff: starter code.

Two integrators written from scratch and one question each:

* an embedded Dormand-Prince 5(4) Runge-Kutta pair with automatic step-size control
  (``rk45_step``, ``solve_rk45``), applied to a ballistic re-entry whose deceleration peak
  lasts a few seconds inside a flight of several minutes: where should the steps go?
* implicit (backward) Euler with a Newton solve at every step (``implicit_euler_step``,
  ``solve_implicit_euler``), applied to a two-node skin/structure thermal model whose time
  constants differ by five orders of magnitude: why does every explicit method crawl, and how
  does an implicit one step over the fast mode?

Conventions:
    * a right-hand side is ``f(t, y) -> dy/dt`` with ``y`` a 1-D array (the ``solve_ivp`` convention);
    * a Jacobian is ``jac(t, y) -> d f_i / d y_j`` (n x n);
    * ``rtol`` and ``atol`` enter through the scale ``atol + rtol * max(|y_n|, |y_{n+1}|)``, as in SciPy;
    * the atmosphere is the exponential model rho = rho0 exp(-h/H); it is a verification model,
      not a prediction of any real entry.

Run:  python main.py          (prints the verification tables, writes adaptive_stiff.png)
Test: pytest curriculum/phase-03/lesson-06-ode-integrators-adaptive-stiff
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np

RHS = Callable[[float, np.ndarray], np.ndarray]
JAC = Callable[[float, np.ndarray], np.ndarray]

# ---------------------------------------------------------------- Dormand-Prince 5(4) coefficients
# c: stage times; A: stage weights (lower triangular); B5: fifth-order weights; B4: fourth-order
# weights. E = B5 - B4 gives the embedded error estimate h * sum(E_i k_i) at no extra cost.
DP_C = np.array([0.0, 1 / 5, 3 / 10, 4 / 5, 8 / 9, 1.0, 1.0])
DP_A = np.array([
    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [1 / 5, 0.0, 0.0, 0.0, 0.0, 0.0],
    [3 / 40, 9 / 40, 0.0, 0.0, 0.0, 0.0],
    [44 / 45, -56 / 15, 32 / 9, 0.0, 0.0, 0.0],
    [19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729, 0.0, 0.0],
    [9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656, 0.0],
    [35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84],
])
DP_B5 = np.array([35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0])
DP_B4 = np.array([5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40])
DP_E = DP_B5 - DP_B4


@dataclass
class IVPResult:
    """What an integrator returns: the accepted times and states plus the work it took."""
    t: np.ndarray
    y: np.ndarray                      # shape (n_states, n_times)
    n_accepted: int = 0
    n_rejected: int = 0
    nfev: int = 0
    n_newton: int = 0
    h: np.ndarray = field(default_factory=lambda: np.zeros(0))   # accepted step sizes


# ---------------------------------------------------------------- one-step methods

def explicit_euler_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """One forward-Euler step y + h f(t, y): first order, stability interval (-2, 0) on the real axis."""
    raise NotImplementedError


def rk4_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """One classical fourth-order Runge-Kutta step (the P03.05 method), for comparison."""
    raise NotImplementedError


def rk45_step(f: RHS, t: float, y: np.ndarray, h: float) -> tuple[np.ndarray, np.ndarray]:
    """One Dormand-Prince 5(4) step.

    Returns ``(y_new, err)``: the fifth-order solution at t + h and the embedded estimate of the
    local error of the *fourth*-order solution, err = y5 - y4 = h * sum_i (b5_i - b4_i) k_i.
    Seven stages, of which the seventh equals the first stage of the next step (FSAL): a careful
    implementation reuses it; this one recomputes it for clarity and counts seven evaluations.
    """
    raise NotImplementedError


def fd_jacobian(f: RHS, t: float, y: np.ndarray, eps: float = 1e-7) -> np.ndarray:
    """Forward-difference Jacobian d f_i / d y_j with step eps * max(1, |y_j|) in each column."""
    raise NotImplementedError


def implicit_euler_step(f: RHS, t: float, y: np.ndarray, h: float, jac: JAC | None = None,
                        tol: float = 1e-10, max_iter: int = 20) -> tuple[np.ndarray, int]:
    """One backward-Euler step: solve z = y + h f(t + h, z) for z by Newton's method.

    G(z) = z - y - h f(t + h, z), G'(z) = I - h J(t + h, z). Newton starts from z = y (the
    explicit-Euler predictor is a bad guess when h is large compared with the fast time
    constant, which is exactly when an implicit method is used). Stop when the Newton update is below ``tol`` (absolute + relative to |z|).
    Returns ``(z, n_newton_iterations)``; raises ``RuntimeError`` if Newton does not converge.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- fixed-step drivers

def solve_fixed(stepper: Callable, f: RHS, t_span: tuple[float, float], y0: np.ndarray, n_steps: int,
                **kwargs) -> IVPResult:
    """Integrate with ``n_steps`` equal steps of a one-step method ``stepper(f, t, y, h, **kwargs)``.

    Accepts ``rk45_step`` (whose first return value is used) and ``implicit_euler_step`` (whose
    second return value is the Newton count) as well as plain steppers returning ``y_new``.
    """
    raise NotImplementedError


def solve_implicit_euler(f: RHS, t_span: tuple[float, float], y0: np.ndarray, n_steps: int,
                         jac: JAC | None = None) -> IVPResult:
    """Backward Euler with ``n_steps`` equal steps; ``nfev`` counts right-hand-side calls inside Newton."""
    raise NotImplementedError


# ---------------------------------------------------------------- adaptive driver

def error_norm(err: np.ndarray, y_old: np.ndarray, y_new: np.ndarray, rtol: float, atol: float) -> float:
    """RMS of the error estimate scaled by atol + rtol * max(|y_old|, |y_new|): accept the step when <= 1."""
    raise NotImplementedError


def initial_step(f: RHS, t0: float, y0: np.ndarray, rtol: float, atol: float, order: int = 5) -> float:
    """A safe first step, chosen from the size of f and of its first difference (the standard recipe).

    With the tolerance scale s = atol + rtol |y0|: d0 = ||y0 / s||, d1 = ||f0 / s||, h0 = 0.01 d0 / d1
    (one percent of the solution per step); then one explicit Euler step to estimate the second
    derivative d2 = ||(f(t0 + h0, y0 + h0 f0) - f0) / s|| / h0 and h1 = (0.01 / max(d1, d2))^(1/order).
    Returns min(100 h0, h1). The RMS norm is used throughout, as in ``error_norm``.
    """
    raise NotImplementedError


def solve_rk45(f: RHS, t_span: tuple[float, float], y0: np.ndarray, rtol: float = 1e-6, atol: float = 1e-9,
               h0: float | None = None, h_max: float = np.inf, max_steps: int = 100_000,
               stop: Callable[[float, np.ndarray], bool] | None = None) -> IVPResult:
    """Dormand-Prince 5(4) with step-size control.

    Each attempted step returns y5 and the error estimate err. With
    E = error_norm(err, y, y5, rtol, atol):
        * E <= 1: accept, then h <- h * min(5, max(0.2, 0.9 E^(-1/5)));
        * E >  1: reject, h <- h * max(0.2, 0.9 E^(-1/5)), and try again from the same point.
    The exponent -1/5 is because the estimated local error is O(h^5). The last step is clipped
    to land on t1 exactly. If ``stop(t, y)`` returns True after an accepted step the integration
    ends there (no root-finding on the event). Raises ``RuntimeError`` after ``max_steps`` attempts.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- linear stability

def stability_function(method: str, z: complex) -> complex:
    """R(z) for the test equation y' = lambda y with z = h lambda: one step from y = 1 with h = 1.

    Computed by *running* the method on the complex scalar equation, not by typing in a formula,
    so it is the stability function of the code as written. ``method`` is one of
    'explicit_euler', 'rk4', 'rk45', 'implicit_euler'.
    """
    raise NotImplementedError


def stability_interval(method: str, x_max: float = 50.0) -> float:
    """Length L of the real stability interval (-L, 0): the largest L with |R(-x)| <= 1 for 0 < x < L.

    Found by scanning in steps of 0.01 and bisecting the first crossing of |R| = 1. Returns
    ``inf`` if |R(-x)| <= 1 all the way to -x_max (the method is A-stable on the real axis).
    """
    raise NotImplementedError


def stiffness_ratio(J: np.ndarray) -> float:
    """max |Re lambda| / min |Re lambda| over the eigenvalues of J (the Jacobian); inf if a real part is zero."""
    raise NotImplementedError


# ---------------------------------------------------------------- problem 1: ballistic re-entry

@dataclass(frozen=True)
class Capsule:
    """A blunt ballistic capsule in an exponential atmosphere (SI)."""
    mass: float = 1200.0          # kg
    cd: float = 1.2               # drag coefficient
    area: float = 4.5             # m^2 reference area
    rho0: float = 1.225           # kg/m^3 sea-level density
    H: float = 7200.0             # m scale height
    R: float = 6.371e6            # m Earth radius
    g0: float = 9.80665           # m/s^2

    @property
    def beta(self) -> float:
        """Ballistic coefficient m / (C_D A) in kg/m^2."""
        return self.mass / (self.cd * self.area)


def density(h: np.ndarray | float, cap: Capsule = Capsule()) -> np.ndarray | float:
    """Exponential atmosphere rho0 exp(-h/H)."""
    return cap.rho0 * np.exp(-np.asarray(h) / cap.H)


def reentry_rhs(t: float, y: np.ndarray, cap: Capsule = Capsule(), gravity: bool = True) -> np.ndarray:
    """Planar point-mass re-entry, state y = [h, v, gamma, s] (altitude m, speed m/s, flight-path
    angle rad, positive upward, ground range m).

        dh/dt = v sin(gamma)
        dv/dt = -rho v^2 / (2 beta) - g sin(gamma)
        dgamma/dt = (v / (R + h) - g / v) cos(gamma)
        ds/dt = v cos(gamma) R / (R + h)

    with g = g0 (R / (R + h))^2. With ``gravity=False`` the g terms *and* the curvature term
    are dropped, which is the straight-line, drag-only model whose solution is the Allen-Eggers
    formula (``allen_eggers_velocity``); use it to verify the integrator, not to fly.
    """
    h, v, gamma, _ = y
    rho = cap.rho0 * np.exp(-h / cap.H)
    drag = rho * v * v / (2.0 * cap.beta)
    if gravity:
        g = cap.g0 * (cap.R / (cap.R + h)) ** 2
        dv = -drag - g * np.sin(gamma)
        dgamma = (v / (cap.R + h) - g / v) * np.cos(gamma)
    else:
        dv = -drag
        dgamma = 0.0
    return np.array([v * np.sin(gamma), dv, dgamma, v * np.cos(gamma) * cap.R / (cap.R + h)])


def allen_eggers_velocity(h: np.ndarray | float, h_e: float, v_e: float, gamma_e: float,
                          cap: Capsule = Capsule()) -> np.ndarray | float:
    """Speed against altitude for straight-line, drag-only entry (Allen-Eggers):

        v(h) = v_e exp( -(rho0 H / (2 beta |sin gamma_e|)) (e^{-h/H} - e^{-h_e/H}) ).

    Derivation: dv/dh = (dv/dt)/(dh/dt) = -rho v^2/(2 beta) / (v sin gamma) is separable.
    """
    raise NotImplementedError


def allen_eggers_peak(v_e: float, gamma_e: float, cap: Capsule = Capsule()) -> tuple[float, float]:
    """Peak deceleration and the altitude at which it occurs, for entry from very high altitude:

        a_max = v_e^2 |sin gamma_e| / (2 e H),   h* = H ln( rho0 H / (beta |sin gamma_e|) ).

    a_max does not depend on the ballistic coefficient; h* does. Both follow from maximising
    a(h) = rho v^2 / (2 beta) with v from ``allen_eggers_velocity`` (set u = rho0 H e^{-h/H} /
    (beta |sin gamma_e|); then a = (|sin gamma_e| v_e^2 / 2H) u e^{-u}, maximal at u = 1).
    """
    raise NotImplementedError


def deceleration(res: IVPResult, cap: Capsule = Capsule()) -> np.ndarray:
    """Aerodynamic deceleration rho v^2 / (2 beta) [m/s^2] along a re-entry result."""
    h, v = res.y[0], res.y[1]
    return density(h, cap) * v ** 2 / (2.0 * cap.beta)


def entry_state(h_e: float = 120e3, v_e: float = 7800.0, gamma_e_deg: float = -6.0) -> np.ndarray:
    """Initial state [h, v, gamma, s] for the lesson's entry case."""
    return np.array([h_e, v_e, np.radians(gamma_e_deg), 0.0])


# ---------------------------------------------------------------- problem 2: stiff two-node thermal model

@dataclass(frozen=True)
class SkinStructure:
    """A thin skin node coupled to a thick structure node (SI).

    dT_s/dt = [ q_in - k_c (T_s - T_b) - eps sigma (T_s^4 - T_inf^4) ] / C_s
    dT_b/dt = [ k_c (T_s - T_b) - k_b (T_b - T_inf) ] / C_b
    """
    C_s: float = 500.0            # J/K heat capacity of the skin node
    C_b: float = 2.0e5            # J/K heat capacity of the structure node
    k_c: float = 2500.0           # W/K skin-to-structure conductance
    k_b: float = 20.0             # W/K structure-to-sink conductance
    q_in: float = 2000.0          # W absorbed by the skin
    eps: float = 0.8              # skin emissivity (0 gives the linear model)
    T_inf: float = 300.0          # K sink temperature
    sigma: float = 5.670374419e-8


def thermal_rhs(t: float, y: np.ndarray, p: SkinStructure = SkinStructure()) -> np.ndarray:
    """Right-hand side of the two-node model, y = [T_s, T_b] in kelvin."""
    Ts, Tb = y
    rad = p.eps * p.sigma * (Ts ** 4 - p.T_inf ** 4)
    return np.array([(p.q_in - p.k_c * (Ts - Tb) - rad) / p.C_s,
                     (p.k_c * (Ts - Tb) - p.k_b * (Tb - p.T_inf)) / p.C_b])


def thermal_jac(t: float, y: np.ndarray, p: SkinStructure = SkinStructure()) -> np.ndarray:
    """Analytic Jacobian of ``thermal_rhs``."""
    Ts = y[0]
    return np.array([[(-p.k_c - 4.0 * p.eps * p.sigma * Ts ** 3) / p.C_s, p.k_c / p.C_s],
                     [p.k_c / p.C_b, (-p.k_c - p.k_b) / p.C_b]])


def thermal_linear_system(p: SkinStructure = SkinStructure()) -> tuple[np.ndarray, np.ndarray]:
    """For eps = 0 the model is y' = A y + c; return (A, c)."""
    A = np.array([[-p.k_c / p.C_s, p.k_c / p.C_s],
                  [p.k_c / p.C_b, -(p.k_c + p.k_b) / p.C_b]])
    c = np.array([p.q_in / p.C_s, p.k_b * p.T_inf / p.C_b])
    return A, c


def thermal_exact_linear(t: np.ndarray | float, y0: np.ndarray, p: SkinStructure = SkinStructure()) -> np.ndarray:
    """Exact solution of the linear (eps = 0) model via the matrix exponential (P02.05):
    y(t) = y_eq + e^{At} (y0 - y_eq), y_eq = -A^{-1} c. Returns shape (2, len(t))."""
    from scipy.linalg import expm
    A, c = thermal_linear_system(p)
    y_eq = -np.linalg.solve(A, c)
    t = np.atleast_1d(np.asarray(t, dtype=float))
    out = np.zeros((2, t.size))
    for i, ti in enumerate(t):
        out[:, i] = y_eq + expm(A * ti) @ (np.asarray(y0) - y_eq)
    return out


# ---------------------------------------------------------------- studies

def convergence_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Least-squares slope of log(err) against log(h): the observed order of accuracy."""
    hs, errs = np.asarray(hs, dtype=float), np.asarray(errs, dtype=float)
    if np.any(errs <= 0):
        raise ValueError("errors must be positive to take logarithms")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


def fixed_step_errors(stepper: Callable, f: RHS, exact: Callable[[float], np.ndarray],
                      t_span: tuple[float, float], y0: np.ndarray, n_steps_list: list[int],
                      component: int | None = None, **kwargs) -> np.ndarray:
    """Global error at t1 for each step count, for a fixed-step convergence study.

    The error is the max norm over all components, or |error| of one ``component`` (use this when
    the components have different units, as in the re-entry state).
    """
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's experiments, print the tables and save ``adaptive_stiff.png``."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    cap = Capsule()
    print(f"capsule: beta = {cap.beta:.1f} kg/m^2")

    # ---- stability intervals of the four methods
    print("\nreal stability interval (-L, 0):")
    for m in ("explicit_euler", "rk4", "rk45", "implicit_euler"):
        print(f"  {m:15s} L = {stability_interval(m):.4f}")

    # ---- re-entry with step control
    y0 = entry_state()
    a_max, h_star = allen_eggers_peak(y0[1], y0[2], cap)
    print(f"\nAllen-Eggers (no gravity): a_max = {a_max:.1f} m/s^2 = {a_max / cap.g0:.2f} g at h* = {h_star / 1e3:.2f} km")
    stop = lambda t, y: y[0] < 20e3
    for rtol in (1e-4, 1e-6, 1e-8):
        res = solve_rk45(lambda t, y: reentry_rhs(t, y, cap), (0.0, 1500.0), y0, rtol=rtol, atol=1e-3, stop=stop)
        a = deceleration(res, cap)
        i = int(np.argmax(a))
        ref = solve_ivp(lambda t, y: reentry_rhs(t, y, cap), (0.0, res.t[-1]), y0, method="DOP853",
                        rtol=1e-12, atol=1e-9, t_eval=res.t)
        err = np.max(np.abs(res.y[1] - ref.y[1]))
        print(f"rtol = {rtol:.0e}: {res.n_accepted} accepted, {res.n_rejected} rejected, nfev = {res.nfev}, "
              f"h in [{res.h.min():.3f}, {res.h.max():.1f}] s, peak {a[i] / cap.g0:.2f} g at t = {res.t[i]:.1f} s, "
              f"h = {res.y[0, i] / 1e3:.2f} km, v = {res.y[1, i]:.0f} m/s; max |v - DOP853| = {err:.2e} m/s")
    res_g = solve_rk45(lambda t, y: reentry_rhs(t, y, cap, gravity=False), (0.0, 1500.0), y0, rtol=1e-8, atol=1e-3, stop=stop)
    v_ae = allen_eggers_velocity(res_g.y[0], y0[0], y0[1], y0[2], cap)
    a_g = deceleration(res_g, cap)
    print(f"gravity off: max |v - Allen-Eggers| = {np.max(np.abs(res_g.y[1] - v_ae)):.2e} m/s, "
          f"peak {np.max(a_g) / cap.g0:.3f} g vs formula {a_max / cap.g0:.3f} g at h = {res_g.y[0, np.argmax(a_g)] / 1e3:.2f} km")

    # ---- fixed-step convergence on the entry problem (short arc through the peak)
    f_entry = lambda t, y: reentry_rhs(t, y, cap)
    t_arc = (0.0, 300.0)
    ref = solve_ivp(f_entry, t_arc, y0, method="DOP853", rtol=1e-13, atol=1e-10).y[:, -1]
    ns = [100, 200, 400, 800]
    for name, st in (("rk4", rk4_step), ("rk45 (fixed h)", rk45_step), ("implicit Euler", implicit_euler_step)):
        errs = fixed_step_errors(st, f_entry, lambda t: ref, t_arc, y0, ns, component=1)
        print(f"{name:15s} error in v at t = {t_arc[1]:.0f} s [m/s] for n = {ns}: {np.array2string(errs, precision=2)}"
              f"  observed order {convergence_order(300.0 / np.array(ns), errs):.2f}")

    # ---- the stiff thermal problem
    p = SkinStructure()
    T0 = np.array([300.0, 300.0])
    J = thermal_jac(0.0, T0, p)
    lam = np.linalg.eigvals(J).real
    print(f"\nthermal Jacobian eigenvalues {lam[0]:.4e}, {lam[1]:.4e} 1/s; stiffness ratio {stiffness_ratio(J):.3e}; "
          f"explicit Euler needs h < {2 / np.max(np.abs(lam)):.3f} s")
    f_th = lambda t, y: thermal_rhs(t, y, p)
    j_th = lambda t, y: thermal_jac(t, y, p)
    T_end = 3600.0
    ref = solve_ivp(f_th, (0.0, T_end), T0, method="Radau", jac=j_th, rtol=1e-12, atol=1e-10)
    Tb_ref = ref.y[1, -1]
    print(f"Radau reference at t = {T_end:.0f} s: T_s = {ref.y[0, -1]:.3f} K, T_b = {Tb_ref:.3f} K ({ref.nfev} evaluations)")
    for h in (0.3, 0.5):
        with np.errstate(over="ignore", invalid="ignore"):
            res = solve_fixed(explicit_euler_step, f_th, (0.0, 60.0), T0, int(60.0 / h))
        print(f"explicit Euler h = {h} s over 60 s: |T_s| at end = {np.abs(res.y[0, -1]):.3e} K  "
              f"({int(T_end / h)} steps per hour would be needed)")
    res45 = solve_rk45(f_th, (0.0, T_end), T0, rtol=1e-6, atol=1e-6)
    print(f"adaptive RK45 rtol 1e-6: {res45.n_accepted} steps, nfev = {res45.nfev}, max h = {res45.h.max():.3f} s, "
          f"error in T_b = {abs(res45.y[1, -1] - Tb_ref):.2e} K")
    for n in (60, 600, 6000):
        res = solve_implicit_euler(f_th, (0.0, T_end), T0, n, jac=j_th)
        print(f"implicit Euler h = {T_end / n:5.1f} s: {n} steps, {res.n_newton} Newton iterations, nfev = {res.nfev}, "
              f"error in T_b = {abs(res.y[1, -1] - Tb_ref):.2e} K, T_s = {res.y[0, -1]:.3f} K")
    errs_ie = []
    for n in (60, 120, 240):
        res = solve_implicit_euler(f_th, (0.0, T_end), T0, n, jac=j_th)
        errs_ie.append(np.max(np.abs(res.y[:, -1] - ref.y[:, -1])))
    print(f"implicit Euler observed order on the thermal problem: {convergence_order(T_end / np.array([60, 120, 240]), errs_ie):.2f}")
    # linear case against the matrix exponential
    p_lin = SkinStructure(eps=0.0)
    ex = thermal_exact_linear(T_end, T0, p_lin)[:, 0]
    res_lin = solve_implicit_euler(lambda t, y: thermal_rhs(t, y, p_lin), (0.0, T_end), T0, 600,
                                   jac=lambda t, y: thermal_jac(t, y, p_lin))
    print(f"linear model (eps = 0) at {T_end:.0f} s: exact T_b = {ex[1]:.3f} K, implicit Euler h = 6 s T_b = {res_lin.y[1, -1]:.3f} K")

    # ---- figure
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    ax = axes[0, 0]
    res = solve_rk45(f_entry, (0.0, 1500.0), y0, rtol=1e-6, atol=1e-3, stop=stop)
    ax.plot(res.t, deceleration(res, cap) / cap.g0, lw=1.2)
    ax.set_xlabel("time [s]"); ax.set_ylabel("deceleration [g]")
    ax2 = ax.twinx()
    ax2.semilogy(res.t[1:], res.h, ".", ms=3, color="tab:red")
    ax2.set_ylabel("accepted step [s]", color="tab:red")
    ax.set_title(f"re-entry, RK45 rtol 1e-6: {res.n_accepted} steps")
    ax = axes[0, 1]
    xs = np.linspace(-5.0, 1.0, 121); ys_ = np.linspace(-4.0, 4.0, 121)
    X, Y = np.meshgrid(xs, ys_)
    for m, col in (("explicit_euler", "tab:blue"), ("rk4", "tab:green"), ("rk45", "tab:orange")):
        Rm = np.array([[abs(stability_function(m, complex(x, yy))) for x in xs] for yy in ys_])
        ax.contour(X, Y, Rm, levels=[1.0], colors=[col])
        ax.plot([], [], color=col, label=m)
    ax.axhline(0, color="k", lw=0.5); ax.axvline(0, color="k", lw=0.5)
    ax.set_xlabel("Re(h lambda)"); ax.set_ylabel("Im(h lambda)"); ax.set_title("|R(z)| = 1 (stable inside); implicit Euler: outside |z - 1| = 1")
    ax.legend(loc="upper left", fontsize=8); ax.set_aspect("equal")
    ax = axes[1, 0]
    res_ie = solve_implicit_euler(f_th, (0.0, T_end), T0, 60, jac=j_th)
    ax.plot(ref.t, ref.y[0], "k-", lw=1, label="Radau (reference) T_s")
    ax.plot(ref.t, ref.y[1], "k--", lw=1, label="Radau T_b")
    ax.plot(res_ie.t, res_ie.y[0], "o", ms=3, label="implicit Euler h = 60 s, T_s")
    ax.plot(res_ie.t, res_ie.y[1], "s", ms=3, label="implicit Euler h = 60 s, T_b")
    ax.set_xlabel("time [s]"); ax.set_ylabel("temperature [K]"); ax.set_title("stiff two-node thermal model"); ax.legend(fontsize=8)
    ax = axes[1, 1]
    rtols = np.logspace(-3, -10, 8)
    nsteps, errs = [], []
    ref_e = solve_ivp(f_entry, (0.0, 300.0), y0, method="DOP853", rtol=1e-13, atol=1e-10).y[1, -1]
    for rt in rtols:
        r = solve_rk45(f_entry, (0.0, 300.0), y0, rtol=rt, atol=rt * 1e-3)
        nsteps.append(r.n_accepted); errs.append(abs(r.y[1, -1] - ref_e))
    ax.loglog(rtols, errs, "o-", label="|error in v(300 s)| [m/s]")
    ax.loglog(rtols, nsteps, "s--", label="accepted steps")
    ax.set_xlabel("rtol"); ax.set_title("re-entry: error and work against tolerance"); ax.legend(fontsize=8); ax.invert_xaxis()
    fig.suptitle("P03.06  Adaptive and stiff integration")
    fig.tight_layout()
    fig.savefig(out_dir / "adaptive_stiff.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'adaptive_stiff.png'}")


if __name__ == "__main__":
    main()
