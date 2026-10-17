"""P03.05 ODE Integrators I: Euler to RK4: reference solution.

Explicit one-step integrators written from scratch: forward Euler, the midpoint and Heun (RK2) methods,
Kutta's third-order rule and the classical fourth-order Runge-Kutta method, both as hand-written steps and
as a generic explicit Runge-Kutta step driven by a Butcher tableau. They are applied to a planar two-body
orbit whose exact solution comes from Kepler's equation (P03.02), so that the global error, the energy and
angular-momentum drift, the observed convergence order and the cost of reaching a position tolerance can all
be measured against the truth. ``scipy.integrate.solve_ivp`` is the oracle for problems without a closed form.

Conventions:
    * An ODE is y' = f(t, y) with y a 1-D NumPy array; f(t, y) returns an array of the same shape.
    * A method of order p has local (one-step) error O(h^(p+1)) and global error O(h^p) on a fixed interval.
    * ``integrate`` takes a fixed nominal step h and adjusts it slightly so that it lands exactly on t1.
    * The stability function of an s-stage explicit RK method of order p = s <= 4 is the truncated
      exponential series R(z) = sum_{k<=p} z^k / k!, z = h * lambda; the method is stable where |R(z)| <= 1.
    * The planar two-body state is y = [x, y, vx, vy] in m and m/s, with mu = G M of the central body.
    * Orbital elements: a [m], e, and the time of perigee passage t_p = 0, perigee on the +x axis.
    * Convergence order is the least-squares slope of log(error) against log(h).

Run:  python solution.py          (prints the order, stability, orbit and work-precision tables, writes rk4_orbit.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-05-ode-integrators-euler-rk4
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np

MU_EARTH: float = 3.986004418e14   # m^3 / s^2
R_EARTH: float = 6.371e6           # m

RHS = Callable[[float, np.ndarray], np.ndarray]


# ---------------------------------------------------------------- Butcher tableaus

class ButcherTableau(NamedTuple):
    """An explicit Runge-Kutta method: stage matrix ``a`` (strictly lower triangular), weights ``b``,
    nodes ``c`` and the method's classical order ``order``."""

    name: str
    a: np.ndarray
    b: np.ndarray
    c: np.ndarray
    order: int


TABLEAUS: dict[str, ButcherTableau] = {
    "euler": ButcherTableau("euler", np.zeros((1, 1)), np.array([1.0]), np.array([0.0]), 1),
    "midpoint": ButcherTableau("midpoint", np.array([[0.0, 0.0], [0.5, 0.0]]), np.array([0.0, 1.0]), np.array([0.0, 0.5]), 2),
    "heun": ButcherTableau("heun", np.array([[0.0, 0.0], [1.0, 0.0]]), np.array([0.5, 0.5]), np.array([0.0, 1.0]), 2),
    "kutta3": ButcherTableau("kutta3", np.array([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [-1.0, 2.0, 0.0]]),
                             np.array([1.0, 4.0, 1.0]) / 6.0, np.array([0.0, 0.5, 1.0]), 3),
    "rk4": ButcherTableau("rk4", np.array([[0.0, 0.0, 0.0, 0.0], [0.5, 0.0, 0.0, 0.0], [0.0, 0.5, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]),
                          np.array([1.0, 2.0, 2.0, 1.0]) / 6.0, np.array([0.0, 0.5, 0.5, 1.0]), 4),
}


def explicit_rk_step(f: RHS, t: float, y: np.ndarray, h: float, tableau: ButcherTableau) -> np.ndarray:
    """One step of the explicit Runge-Kutta method given by ``tableau``:
    k_i = f(t + c_i h, y + h sum_j a_ij k_j), y_next = y + h sum_i b_i k_i."""
    y = np.asarray(y, dtype=float)
    s = tableau.b.size
    k = np.zeros((s, y.size))
    for i in range(s):
        yi = y + h * (tableau.a[i, :i] @ k[:i]) if i > 0 else y
        k[i] = f(t + tableau.c[i] * h, yi)
    return y + h * (tableau.b @ k)


def euler_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Forward Euler: y + h f(t, y)."""
    return np.asarray(y, dtype=float) + h * f(t, y)


def midpoint_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Explicit midpoint (RK2): a half Euler step to the midpoint, then a full step with the midpoint slope."""
    y = np.asarray(y, dtype=float)
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    return y + h * k2


def heun_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Heun (improved Euler, RK2): average of the slopes at the start and at the Euler predictor."""
    y = np.asarray(y, dtype=float)
    k1 = f(t, y)
    k2 = f(t + h, y + h * k1)
    return y + 0.5 * h * (k1 + k2)


def rk4_step(f: RHS, t: float, y: np.ndarray, h: float) -> np.ndarray:
    """The classical fourth-order Runge-Kutta step, Simpson's rule on the slopes."""
    y = np.asarray(y, dtype=float)
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


STEPPERS: dict[str, Callable[[RHS, float, np.ndarray, float], np.ndarray]] = {
    "euler": euler_step,
    "midpoint": midpoint_step,
    "heun": heun_step,
    "kutta3": lambda f, t, y, h: explicit_rk_step(f, t, y, h, TABLEAUS["kutta3"]),
    "rk4": rk4_step,
}


def integrate(f: RHS, t_span: tuple[float, float], y0: np.ndarray, h: float, method: str = "rk4") -> tuple[np.ndarray, np.ndarray]:
    """Integrate y' = f(t, y) from t0 to t1 with n = ceil((t1 - t0) / h) equal steps of the named method.

    Returns ``(t, Y)`` with ``t`` of shape (n + 1,) ending exactly at t1 and ``Y`` of shape (n + 1, len(y0)).
    Raise ``ValueError`` for an unknown method, h <= 0 or t1 <= t0.
    """
    if method not in STEPPERS:
        raise ValueError(f"unknown method {method!r}; choose from {sorted(STEPPERS)}")
    t0, t1 = float(t_span[0]), float(t_span[1])
    if h <= 0.0 or t1 <= t0:
        raise ValueError("need h > 0 and t1 > t0")
    n = max(1, int(math.ceil((t1 - t0) / h - 1e-12)))
    hh = (t1 - t0) / n
    step = STEPPERS[method]
    y = np.array(y0, dtype=float)
    t = t0 + hh * np.arange(n + 1)
    Y = np.empty((n + 1, y.size))
    Y[0] = y
    for i in range(n):
        y = step(f, t[i], y, hh)
        Y[i + 1] = y
    t[-1] = t1
    return t, Y


def stages(method: str) -> int:
    """Number of f evaluations per step of the named method."""
    return TABLEAUS[method].b.size


# ---------------------------------------------------------------- error measurement

def observed_order(hs: np.ndarray, errs: np.ndarray) -> float:
    """Return the least-squares slope of log(err) against log(h). Raise ``ValueError`` for < 2 points or non-positive entries."""
    hs = np.asarray(hs, dtype=float)
    errs = np.asarray(errs, dtype=float)
    if hs.size < 2 or hs.size != errs.size or np.any(hs <= 0.0) or np.any(errs <= 0.0):
        raise ValueError("need >= 2 positive (h, err) pairs")
    return float(np.polyfit(np.log(hs), np.log(errs), 1)[0])


def local_error(f: RHS, exact: Callable[[float], np.ndarray], t0: float, h: float, method: str) -> float:
    """Error of one step of ``method`` started from the exact solution at t0, measured against exact(t0 + h)
    (2-norm). For a method of order p this is O(h^(p+1))."""
    y1 = STEPPERS[method](f, t0, exact(t0), h)
    return float(np.linalg.norm(y1 - exact(t0 + h)))


def global_error(f: RHS, exact: Callable[[float], np.ndarray], t_span: tuple[float, float], h: float, method: str) -> float:
    """2-norm error at t1 of ``integrate`` with nominal step h, against ``exact(t1)``."""
    _, Y = integrate(f, t_span, exact(t_span[0]), h, method)
    return float(np.linalg.norm(Y[-1] - exact(t_span[1])))


def step_for_tolerance(h_ref: float, err_ref: float, order: int, tol: float) -> float:
    """Step predicted to bring the global error from ``err_ref`` (observed at ``h_ref``) down to ``tol`` for a
    method of the given order: h = h_ref (tol / err_ref)^(1 / order)."""
    if h_ref <= 0.0 or err_ref <= 0.0 or tol <= 0.0 or order < 1:
        raise ValueError("need positive h_ref, err_ref, tol and order >= 1")
    return h_ref * (tol / err_ref) ** (1.0 / order)


# ---------------------------------------------------------------- linear stability

def stability_function(z: complex | np.ndarray, order: int) -> complex | np.ndarray:
    """R(z) = sum_{k=0}^{p} z^k / k!: the amplification factor per step of an explicit RK method with s = p
    stages and order p (Euler p = 1, midpoint/Heun p = 2, Kutta p = 3, RK4 p = 4) applied to y' = lambda y
    with z = h lambda. Raise ``ValueError`` unless 1 <= order <= 4."""
    if not 1 <= order <= 4:
        raise ValueError("order must be 1..4 (s = p explicit RK)")
    z = np.asarray(z)
    out = np.ones_like(z, dtype=complex)
    term = np.ones_like(z, dtype=complex)
    for k in range(1, order + 1):
        term = term * z / k
        out = out + term
    return out if out.ndim else complex(out)


def stability_limit_real(order: int, tol: float = 1e-10) -> float:
    """Largest x > 0 such that |R(-x')| <= 1 for all 0 <= x' <= x: the real-axis stability limit
    (2 for Euler and the RK2 methods, 2.5127 for Kutta3, 2.7853 for RK4), by bisection on a grid scan."""
    xs = np.linspace(0.0, 4.0, 4001)
    unstable = np.abs(stability_function(-xs, order)) > 1.0 + 1e-13
    first = int(np.argmax(unstable))
    lo, hi = xs[first - 1], xs[first]
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if abs(stability_function(-mid, order)) <= 1.0:
            lo = mid
        else:
            hi = mid
    return float(0.5 * (lo + hi))


def max_stable_step(order: int, lam: float) -> float:
    """Largest step for which the method of the given order is stable on y' = lam y with lam < 0:
    stability_limit_real(order) / |lam|."""
    if lam >= 0.0:
        raise ValueError("need lam < 0")
    return stability_limit_real(order) / abs(lam)


# ---------------------------------------------------------------- the two-body problem

def two_body_rhs(t: float, y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Planar two-body right-hand side for y = [x, y, vx, vy]: r'' = -mu r / |r|^3."""
    r3 = (y[0] * y[0] + y[1] * y[1]) ** 1.5
    return np.array([y[2], y[3], -mu * y[0] / r3, -mu * y[1] / r3])


def solve_kepler(M: float, e: float, tol: float = 1e-14) -> float:
    """Eccentric anomaly E from Kepler's equation M = E - e sin E by Newton's method (P03.02), 0 <= e < 1."""
    if not 0.0 <= e < 1.0:
        raise ValueError("elliptic orbits only")
    E = M + e * math.sin(M)
    for _ in range(50):
        d = (E - e * math.sin(E) - M) / (1.0 - e * math.cos(E))
        E -= d
        if abs(d) < tol:
            break
    return E


def orbit_period(a: float, mu: float = MU_EARTH) -> float:
    """T = 2 pi sqrt(a^3 / mu)."""
    return 2.0 * math.pi * math.sqrt(a ** 3 / mu)


def kepler_state(t: float, a: float, e: float, mu: float = MU_EARTH) -> np.ndarray:
    """Exact planar state [x, y, vx, vy] at time t of an orbit with semi-major axis a and eccentricity e,
    perigee on +x at t = 0, counter-clockwise motion."""
    n = math.sqrt(mu / a ** 3)
    E = solve_kepler(math.remainder(n * t, 2.0 * math.pi), e)
    b = a * math.sqrt(1.0 - e * e)
    x = a * (math.cos(E) - e)
    yy = b * math.sin(E)
    r = a * (1.0 - e * math.cos(E))
    Edot = n * a / r
    return np.array([x, yy, -a * math.sin(E) * Edot, b * math.cos(E) * Edot])


def specific_energy(y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray | float:
    """v^2 / 2 - mu / r for one state (shape (4,)) or a stack of states (shape (n, 4))."""
    y = np.asarray(y, dtype=float)
    r = np.hypot(y[..., 0], y[..., 1])
    v2 = y[..., 2] ** 2 + y[..., 3] ** 2
    return 0.5 * v2 - mu / r


def angular_momentum(y: np.ndarray) -> np.ndarray | float:
    """x vy - y vx for one state or a stack of states."""
    y = np.asarray(y, dtype=float)
    return y[..., 0] * y[..., 3] - y[..., 1] * y[..., 2]


def propagate_orbit(a: float, e: float, n_orbits: float, steps_per_orbit: int, method: str = "rk4",
                    mu: float = MU_EARTH) -> tuple[np.ndarray, np.ndarray]:
    """Integrate the two-body problem from perigee for ``n_orbits`` periods with ``steps_per_orbit`` steps
    per period; returns (t, Y)."""
    T = orbit_period(a, mu)
    f = lambda t, y: two_body_rhs(t, y, mu)
    return integrate(f, (0.0, n_orbits * T), kepler_state(0.0, a, e, mu), T / steps_per_orbit, method)


def orbit_position_error(a: float, e: float, n_orbits: float, steps_per_orbit: int, method: str = "rk4",
                         mu: float = MU_EARTH) -> float:
    """Final position error [m] of ``propagate_orbit`` against ``kepler_state``."""
    t, Y = propagate_orbit(a, e, n_orbits, steps_per_orbit, method, mu)
    ref = kepler_state(t[-1], a, e, mu)
    return float(np.hypot(Y[-1, 0] - ref[0], Y[-1, 1] - ref[1]))


def energy_drift(Y: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Relative energy error (eps(t) - eps(0)) / |eps(0)| along a trajectory."""
    eps = specific_energy(Y, mu)
    return (eps - eps[0]) / abs(eps[0])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Measure local and global orders, stability limits, orbit errors and cost, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    methods = ["euler", "midpoint", "heun", "kutta3", "rk4"]

    # --- orders on the harmonic oscillator y'' = -y, y = [cos t, -sin t]
    f_h = lambda t, y: np.array([y[1], -y[0]])
    ex_h = lambda t: np.array([math.cos(t), -math.sin(t)])
    hs_loc = np.array([0.2, 0.1, 0.05])
    hs_glob = 2.0 * math.pi / np.array([50, 100, 200, 400])
    print("harmonic oscillator over one period: local order (one step) and global order (error at t = 2 pi)")
    print(f"{'method':>9}{'stages':>7}{'local':>8}{'global':>8}{'err h=2pi/400':>16}")
    for m in methods:
        el = [local_error(f_h, ex_h, 0.3, h, m) for h in hs_loc]
        eg = [global_error(f_h, ex_h, (0.0, 2.0 * math.pi), h, m) for h in hs_glob]
        print(f"{m:>9}{stages(m):>7}{observed_order(hs_loc, el):>8.2f}{observed_order(hs_glob, eg):>8.2f}{eg[-1]:>16.2e}")

    # --- the RK4 step on y' = y is the truncated exponential
    f_e = lambda t, y: y
    y1 = rk4_step(f_e, 0.0, np.array([1.0]), 0.5)[0]
    print(f"RK4 on y' = y, h = 0.5: {y1:.12f}; 1 + h + h^2/2 + h^3/6 + h^4/24 = {sum(0.5 ** k / math.factorial(k) for k in range(5)):.12f}; "
          f"exp(0.5) = {math.exp(0.5):.12f}; local error {y1 - math.exp(0.5):.2e} vs h^5/120 = {0.5 ** 5 / 120:.2e}")

    # --- stability
    print("real-axis stability limits |h lambda|:", ", ".join(f"order {p}: {stability_limit_real(p):.4f}" for p in range(1, 5)))
    lam = -1000.0
    for p, m in [(1, "euler"), (4, "rk4")]:
        hmax = max_stable_step(p, lam)
        for h in (0.95 * hmax, 1.05 * hmax):
            _, Y = integrate(lambda t, y: lam * y, (0.0, 0.1), np.array([1.0]), h, m)
            print(f"  y' = -1000 y, {m:>5}, h = {h:.5f} ({'0.95' if h < hmax else '1.05'} h_max): |y(0.1)| = {abs(Y[-1, 0]):.2e}")

    # --- two-body orbit
    a, e = 7.5e6, 0.1
    T = orbit_period(a)
    print(f"orbit a = {a / 1e3:.0f} km, e = {e}: perigee altitude {(a * (1 - e) - R_EARTH) / 1e3:.0f} km, period {T:.1f} s")
    y0 = kepler_state(0.0, a, e)
    sol = solve_ivp(two_body_rhs, (0.0, T), y0, method="DOP853", rtol=1e-13, atol=1e-6)
    ref = kepler_state(T, a, e)
    print(f"Kepler reference vs solve_ivp(DOP853, rtol 1e-13) after one orbit: {np.hypot(*(sol.y[:2, -1] - ref[:2])):.1e} m")
    spo = [100, 200, 400, 800]
    print(f"{'steps/orbit':>12}" + "".join(f"{m:>14}" for m in methods) + "   (position error after 1 orbit, m)")
    errs = {m: [] for m in methods}
    for n in spo:
        row = f"{n:>12}"
        for m in methods:
            errs[m].append(orbit_position_error(a, e, 1.0, n, m))
            row += f"{errs[m][-1]:>14.3e}"
        print(row)
    hs = T / np.array(spo, dtype=float)
    print("observed orders:", ", ".join(f"{m} {observed_order(hs, errs[m]):.2f}" for m in methods))
    fine = [3200, 6400, 12800]
    e_fine = [orbit_position_error(a, e, 1.0, n, "euler") for n in fine]
    print("Euler is not yet asymptotic above: at", fine, "steps/orbit its errors are",
          ", ".join(f"{x:.3e}" for x in e_fine), f"m, order {observed_order(T / np.array(fine, dtype=float), e_fine):.2f}")
    # drift over ten orbits
    print("ten orbits at 400 steps/orbit: final position error and relative energy / angular-momentum drift")
    drift = {}
    for m in ["euler", "heun", "rk4"]:
        t, Y = propagate_orbit(a, e, 10.0, 400, m)
        de = energy_drift(Y)
        dh = (angular_momentum(Y) - angular_momentum(Y[0])) / abs(angular_momentum(Y[0]))
        pe = np.hypot(Y[-1, 0] - kepler_state(t[-1], a, e)[0], Y[-1, 1] - kepler_state(t[-1], a, e)[1])
        drift[m] = (t, Y, de)
        print(f"  {m:>5}: position {pe:.3e} m, energy {de[-1]:+.3e}, angular momentum {dh[-1]:+.3e}")
    # cost for a tolerance
    print("steps per orbit and f evaluations predicted for a 1 m position error after one orbit:")
    for m in methods:
        p = TABLEAUS[m].order
        h_req = step_for_tolerance(hs[-1], errs[m][-1], p, 1.0)
        n_req = T / h_req
        print(f"  {m:>9}: {n_req:12.0f} steps, {n_req * stages(m):14.0f} evaluations")
    n_rk4 = int(math.ceil(T / step_for_tolerance(hs[-1], errs["rk4"][-1], 4, 1.0)))
    print(f"check: RK4 with {n_rk4} steps/orbit gives {orbit_position_error(a, e, 1.0, n_rk4):.3f} m")

    # --- figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ax = axes[0]
    for m in methods:
        ax.loglog(spo, errs[m], "o-", label=f"{m} ({observed_order(hs, errs[m]):.2f})")
    ax.set_xlabel("steps per orbit")
    ax.set_ylabel("position error after one orbit [m]")
    ax.set_title("Global error against step")
    ax.legend(fontsize=8)
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    for m, (t, Y, de) in drift.items():
        ax.plot(t / T, de, label=m)
    ax.set_yscale("symlog", linthresh=1e-10)
    ax.set_xlabel("orbits")
    ax.set_ylabel("relative energy error")
    ax.set_title("Energy drift, 400 steps/orbit")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax = axes[2]
    xg, yg = np.meshgrid(np.linspace(-3.2, 0.8, 321), np.linspace(-3.2, 3.2, 321))
    zg = xg + 1j * yg
    for p, col in zip(range(1, 5), ["C0", "C1", "C2", "C3"]):
        ax.contour(xg, yg, np.abs(stability_function(zg, p)), levels=[1.0], colors=col)
        ax.plot([], [], color=col, label=f"order {p}: real limit {stability_limit_real(p):.3f}")
    ax.set_aspect("equal")
    ax.set_xlabel("Re(h lambda)")
    ax.set_ylabel("Im(h lambda)")
    ax.set_title("Stability regions |R(z)| = 1")
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "rk4_orbit.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
