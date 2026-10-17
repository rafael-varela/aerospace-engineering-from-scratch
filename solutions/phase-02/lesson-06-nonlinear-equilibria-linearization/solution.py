"""P02.06 Nonlinear ODEs, Equilibria & Linearization: reference solution.

The compound pendulum of an aircraft inertia swing test, treated as a nonlinear autonomous
system: equilibria, the Jacobian, classification by the linear theory, an RK4 integrator
written from scratch, the exact large-amplitude period, and the error the small-angle
(linearised) model makes when it is used to infer a moment of inertia.

Conventions:
    * State y = (theta, theta_dot) [rad, rad/s]. The pendulum equation is written as
          theta'' + 2 zeta omega0 theta' + omega0**2 sin(theta) = 0,
      omega0 = sqrt(g / L_eff) the small-angle natural frequency [rad/s], zeta the damping ratio.
    * For a compound pendulum (a rigid body swung about a pivot a distance d from its centre of
      mass) L_eff = I_pivot / (m d), with I_pivot = I_cg + m d**2.
    * All right-hand sides have the SciPy signature f(t, y, ...) so solve_ivp can be the oracle.

Run:  python solution.py          (prints the swing-test tables, writes pendulum_linearization.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-06-nonlinear-equilibria-linearization
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.special import ellipk

G0 = 9.80665  # standard gravity [m/s^2]


def uav_swing_test() -> dict[str, float]:
    """Return the nominal swing-test case: a 4.2 kg fixed-wing UAV hung from a knife edge.

    m [kg], d [m] from the pivot to the centre of mass, I_cg [kg m^2] about the pitch axis through
    the centre of mass (the value the test is meant to recover), zeta the damping ratio of the swing.
    """
    return {"m": 4.2, "d": 0.35, "I_cg": 0.21, "zeta": 0.01}


# ---------------------------------------------------------------- the compound pendulum

def effective_length(I_cg: float, m: float, d: float) -> float:
    """Return L_eff = (I_cg + m d^2) / (m d), the simple-pendulum length with the same small-angle period.

    Raise ``ValueError`` unless m, d and I_cg are positive.
    """
    if m <= 0.0 or d <= 0.0 or I_cg <= 0.0:
        raise ValueError("mass, pivot distance and inertia must be positive")
    return (I_cg + m * d * d) / (m * d)


def natural_frequency(L_eff: float, g: float = G0) -> float:
    """Return omega0 = sqrt(g / L_eff) [rad/s], the small-angle natural frequency."""
    if L_eff <= 0.0 or g <= 0.0:
        raise ValueError("length and gravity must be positive")
    return float(np.sqrt(g / L_eff))


def inertia_from_period(T: float, m: float, d: float, g: float = G0) -> float:
    """Return I_cg inferred from a measured period with the small-angle formula.

    I_pivot = m g d T^2 / (4 pi^2); I_cg = I_pivot - m d^2. This is what a swing test reports; it is
    exact only in the limit of vanishing amplitude.
    """
    if T <= 0.0:
        raise ValueError("period must be positive")
    return m * g * d * T * T / (4.0 * np.pi ** 2) - m * d * d


# ---------------------------------------------------------------- the nonlinear system

def pendulum_rhs(t: float, y: np.ndarray, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return dy/dt for y = (theta, theta_dot): theta'' = -2 zeta omega0 theta' - omega0^2 sin(theta)."""
    theta, theta_dot = y
    return np.array([theta_dot, -2.0 * zeta * omega0 * theta_dot - omega0 ** 2 * np.sin(theta)])


def linear_pendulum_rhs(t: float, y: np.ndarray, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return dy/dt for the linearised pendulum, sin(theta) replaced by theta."""
    theta, theta_dot = y
    return np.array([theta_dot, -2.0 * zeta * omega0 * theta_dot - omega0 ** 2 * theta])


def pendulum_equilibria() -> list[np.ndarray]:
    """Return the two distinct equilibria of the pendulum in [-pi, pi): hanging (0, 0) and inverted (pi, 0)."""
    return [np.array([0.0, 0.0]), np.array([np.pi, 0.0])]


def jacobian_numeric(f: Callable[[np.ndarray], np.ndarray], x: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Return the Jacobian matrix of a vector field f at x by central differences, J[i, j] = df_i/dx_j.

    ``f`` takes the state alone; wrap a (t, y, ...) right-hand side with a lambda. The step ``eps`` is
    absolute; the error is O(eps^2) plus rounding of order 1e-16 / eps.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    f0 = np.asarray(f(x), dtype=float)
    J = np.zeros((f0.size, n))
    for j in range(n):
        dx = np.zeros(n)
        dx[j] = eps
        J[:, j] = (np.asarray(f(x + dx), dtype=float) - np.asarray(f(x - dx), dtype=float)) / (2.0 * eps)
    return J


def pendulum_jacobian(theta: float, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return the analytic Jacobian [[0, 1], [-omega0^2 cos(theta), -2 zeta omega0]] of pendulum_rhs at (theta, any theta_dot)."""
    return np.array([[0.0, 1.0], [-omega0 ** 2 * np.cos(theta), -2.0 * zeta * omega0]])


def pendulum_jacobian_sympy(theta: float, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return the same Jacobian computed symbolically with SymPy (``Matrix.jacobian``), then evaluated.

    An independent check that the hand derivative is right; the lesson uses it once and then keeps
    the hand version because it is a thousand times faster to evaluate.
    """
    import sympy as sp

    th, thd, w, z = sp.symbols("theta theta_dot omega0 zeta", real=True)
    F = sp.Matrix([thd, -2 * z * w * thd - w ** 2 * sp.sin(th)])
    J = F.jacobian(sp.Matrix([th, thd]))
    return np.array(J.subs({th: theta, w: omega0, z: zeta}).evalf(), dtype=float)


def classify_equilibrium(J: np.ndarray, tol: float = 1e-12) -> str:
    """Classify the linear system y' = J y for a 2x2 Jacobian by its trace and determinant.

    Returns one of 'saddle', 'center', 'stable spiral', 'unstable spiral', 'stable node',
    'unstable node', 'degenerate' (zero determinant or trace^2 = 4 det within ``tol``). A center is
    reported only when the trace is zero within ``tol``; the lesson explains why that verdict is the
    one the linearisation cannot be trusted to deliver for the nonlinear system.
    """
    J = np.asarray(J, dtype=float)
    if J.shape != (2, 2):
        raise ValueError("classify_equilibrium is written for 2x2 Jacobians")
    tr = J[0, 0] + J[1, 1]
    det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    scale = max(1.0, np.abs(J).max() ** 2)
    if abs(det) <= tol * scale:
        return "degenerate"
    if det < 0.0:
        return "saddle"
    disc = tr * tr - 4.0 * det
    if abs(tr) <= tol * np.sqrt(scale):
        return "center"
    if abs(disc) <= tol * scale:
        return "degenerate"
    if disc < 0.0:
        return "stable spiral" if tr < 0.0 else "unstable spiral"
    return "stable node" if tr < 0.0 else "unstable node"


def pendulum_eigenvalues(theta: float, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return the two eigenvalues of the pendulum Jacobian at ``theta`` (numpy.linalg.eigvals, complex dtype)."""
    return np.linalg.eigvals(pendulum_jacobian(theta, omega0, zeta))


# ---------------------------------------------------------------- integration from scratch

def rk4_step(f: Callable[[float, np.ndarray], np.ndarray], t: float, y: np.ndarray, h: float) -> np.ndarray:
    """Advance y' = f(t, y) by one classical fourth-order Runge-Kutta step of size h."""
    k1 = np.asarray(f(t, y), dtype=float)
    k2 = np.asarray(f(t + 0.5 * h, y + 0.5 * h * k1), dtype=float)
    k3 = np.asarray(f(t + 0.5 * h, y + 0.5 * h * k2), dtype=float)
    k4 = np.asarray(f(t + h, y + h * k3), dtype=float)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_rk4(f: Callable[[float, np.ndarray], np.ndarray], y0: np.ndarray, t_end: float,
                  h: float, t0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Integrate y' = f(t, y) from t0 to t_end with fixed-step RK4.

    Returns (t, Y) with t of length N + 1 and Y of shape (N + 1, len(y0)); the step count N is
    ceil((t_end - t0) / h) and the last step is shortened to land exactly on t_end.
    Raise ``ValueError`` for a non-positive step or t_end <= t0.
    """
    if h <= 0.0 or t_end <= t0:
        raise ValueError("need h > 0 and t_end > t0")
    n = int(np.ceil((t_end - t0) / h - 1e-12))
    t = np.empty(n + 1)
    Y = np.empty((n + 1, np.size(y0)))
    t[0] = t0
    Y[0] = np.asarray(y0, dtype=float)
    for k in range(n):
        step = min(h, t_end - t[k])
        Y[k + 1] = rk4_step(f, t[k], Y[k], step)
        t[k + 1] = t[k] + step
    t[-1] = t_end
    return t, Y


def linear_response(t: np.ndarray, theta0: float, omega0: float, zeta: float = 0.0) -> np.ndarray:
    """Return theta(t) of the linearised pendulum released from rest at theta0 (closed form, 0 <= zeta < 1).

    theta = theta0 e^{-zeta omega0 t} [cos(omega_d t) + (zeta / sqrt(1 - zeta^2)) sin(omega_d t)],
    omega_d = omega0 sqrt(1 - zeta^2). Raise ``ValueError`` for zeta outside [0, 1).
    """
    if not 0.0 <= zeta < 1.0:
        raise ValueError("closed form written for the underdamped case 0 <= zeta < 1")
    t = np.asarray(t, dtype=float)
    wd = omega0 * np.sqrt(1.0 - zeta ** 2)
    return theta0 * np.exp(-zeta * omega0 * t) * (np.cos(wd * t) + zeta / np.sqrt(1.0 - zeta ** 2) * np.sin(wd * t))


# ---------------------------------------------------------------- energy, separatrix, period

def energy(theta: np.ndarray, theta_dot: np.ndarray, omega0: float) -> np.ndarray:
    """Return the specific energy E = theta_dot^2 / 2 + omega0^2 (1 - cos theta) [rad^2/s^2]; conserved when zeta = 0."""
    return 0.5 * np.asarray(theta_dot, dtype=float) ** 2 + omega0 ** 2 * (1.0 - np.cos(np.asarray(theta, dtype=float)))


def separatrix_energy(omega0: float) -> float:
    """Return 2 omega0^2, the energy of the orbit through the inverted equilibrium that separates swinging from rotating."""
    return 2.0 * omega0 ** 2


def is_rotating(theta: float, theta_dot: float, omega0: float) -> bool:
    """Return True if the undamped pendulum with this state goes over the top instead of swinging back."""
    return bool(energy(theta, theta_dot, omega0) > separatrix_energy(omega0))


def period_small_angle(omega0: float) -> float:
    """Return T0 = 2 pi / omega0, the period of the linearised pendulum."""
    return 2.0 * np.pi / omega0


def period_exact(theta0: float, omega0: float) -> float:
    """Return the exact period of the undamped pendulum released from rest at amplitude theta0 (0 <= theta0 < pi).

    T = (4 / omega0) K(k), k = sin(theta0 / 2), with K the complete elliptic integral of the first kind
    (scipy.special.ellipk takes the parameter m = k^2). Raise ``ValueError`` outside [0, pi).
    """
    if not 0.0 <= theta0 < np.pi:
        raise ValueError("amplitude must lie in [0, pi)")
    k = np.sin(0.5 * theta0)
    return float(4.0 / omega0 * ellipk(k * k))


def period_series(theta0: float, omega0: float) -> float:
    """Return the three-term series T0 (1 + theta0^2/16 + 11 theta0^4/3072) for the period at amplitude theta0."""
    return period_small_angle(omega0) * (1.0 + theta0 ** 2 / 16.0 + 11.0 * theta0 ** 4 / 3072.0)


def period_from_trajectory(t: np.ndarray, theta: np.ndarray) -> float:
    """Return the mean period from the upward zero crossings of theta(t), each located by linear interpolation.

    Needs at least two upward crossings; raise ``ValueError`` otherwise.
    """
    t = np.asarray(t, dtype=float)
    th = np.asarray(theta, dtype=float)
    idx = np.where((th[:-1] < 0.0) & (th[1:] >= 0.0))[0]
    if idx.size < 2:
        raise ValueError("fewer than two upward zero crossings in the trajectory")
    tc = t[idx] - th[idx] * (t[idx + 1] - t[idx]) / (th[idx + 1] - th[idx])
    return float((tc[-1] - tc[0]) / (tc.size - 1))


def max_amplitude_for_period_error(tol: float, lo: float = 0.0, hi: float = 3.0) -> float:
    """Return the largest release amplitude [rad] at which the exact period exceeds T0 by at most ``tol`` (relative).

    Found by bisection on period_exact / T0 - 1, which is monotone in the amplitude. Independent of omega0.
    """
    if tol <= 0.0:
        raise ValueError("tolerance must be positive")
    f = lambda a: period_exact(a, 1.0) / period_small_angle(1.0) - 1.0 - tol
    if f(hi) < 0.0:
        return hi
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- a limit cycle, for contrast

def van_der_pol_rhs(t: float, y: np.ndarray, mu: float) -> np.ndarray:
    """Return dy/dt for the Van der Pol oscillator x'' - mu (1 - x^2) x' + x = 0, y = (x, x')."""
    x, v = y
    return np.array([v, mu * (1.0 - x * x) * v - x])


def limit_cycle_amplitude(mu: float, y0: np.ndarray, t_end: float = 60.0, h: float = 0.01) -> float:
    """Return the amplitude max|x| over the last 20 % of an RK4 run of the Van der Pol oscillator from y0.

    For every non-zero start the trajectory settles on the same closed orbit (amplitude 2.0 as mu -> 0);
    the origin is an unstable spiral for 0 < mu < 2 and no linear analysis predicts the closed orbit.
    """
    t, Y = integrate_rk4(lambda tt, yy: van_der_pol_rhs(tt, yy, mu), np.asarray(y0, dtype=float), t_end, h)
    tail = Y[int(0.8 * Y.shape[0]):, 0]
    return float(np.max(np.abs(tail)))


# ---------------------------------------------------------------- the swing-test study

def swing_test_error(theta0: float, case: dict[str, float], g: float = G0) -> dict[str, float]:
    """Return the inertia error a small-angle swing test makes when released at amplitude theta0.

    Keys: 'T0' (small-angle period), 'T_exact', 'I_cg_inferred' from T_exact with the small-angle formula,
    'rel_error' = I_cg_inferred / I_cg - 1. Undamped; light damping changes the period only at second
    order in zeta.
    """
    L = effective_length(case["I_cg"], case["m"], case["d"])
    w0 = natural_frequency(L, g)
    T0 = period_small_angle(w0)
    T = period_exact(theta0, w0)
    I_inf = inertia_from_period(T, case["m"], case["d"], g)
    return {"T0": T0, "T_exact": T, "I_cg_inferred": I_inf, "rel_error": I_inf / case["I_cg"] - 1.0}


def max_amplitude_for_inertia_error(tol: float, case: dict[str, float], g: float = G0) -> float:
    """Return the largest release amplitude [rad] at which the swing test's inertia error stays below ``tol`` (relative).

    Bisection on swing_test_error(...)['rel_error'] - tol. The error in I_cg is larger than the error in
    the period by a factor of about 2 I_pivot / I_cg, because the test subtracts m d^2 from I_pivot.
    """
    if tol <= 0.0:
        raise ValueError("tolerance must be positive")
    lo, hi = 0.0, 3.0
    f = lambda a: swing_test_error(a, case, g)["rel_error"] - tol
    if f(hi) < 0.0:
        return hi
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


def main(out_dir: str | Path = ".") -> None:
    """Print the swing-test analysis and save the phase portrait and amplitude-error figure."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.integrate import solve_ivp

    out_dir = Path(out_dir)
    np.set_printoptions(precision=5, suppress=True)
    case = uav_swing_test()
    L = effective_length(case["I_cg"], case["m"], case["d"])
    w0 = natural_frequency(L)
    zeta = case["zeta"]
    print(f"UAV swing test: m = {case['m']} kg, d = {case['d']} m, I_cg = {case['I_cg']} kg m^2")
    print(f"I_pivot = {case['I_cg'] + case['m'] * case['d'] ** 2:.4f} kg m^2, L_eff = {L:.4f} m,"
          f" omega0 = {w0:.4f} rad/s, T0 = {period_small_angle(w0):.4f} s")

    print("\nequilibria and their linearisations (zeta = {:.2f}):".format(zeta))
    for eq in pendulum_equilibria():
        J = pendulum_jacobian(eq[0], w0, zeta)
        Jn = jacobian_numeric(lambda x: pendulum_rhs(0.0, x, w0, zeta), eq)
        print(f"  theta = {eq[0]:.4f}: J = {J.tolist()}, |J - J_numeric| = {np.abs(J - Jn).max():.1e},"
              f" eigenvalues {pendulum_eigenvalues(eq[0], w0, zeta)}, class: {classify_equilibrium(J)}")
    print("  undamped hanging equilibrium:", classify_equilibrium(pendulum_jacobian(0.0, w0, 0.0)))

    print("\nsmall-angle swing test: inertia error against release amplitude")
    print("  amp [deg]   T_exact/T0     I_cg inferred [kg m^2]   error")
    for deg in (5.0, 10.0, 20.0, 30.0, 45.0, 60.0, 90.0):
        r = swing_test_error(np.radians(deg), case)
        print(f"  {deg:6.1f}      {r['T_exact'] / r['T0']:.5f}        {r['I_cg_inferred']:.5f}            {100 * r['rel_error']:+.2f} %")
    lever = (case["I_cg"] + case["m"] * case["d"] ** 2) / case["I_cg"]
    print(f"  amplification I_pivot / I_cg = {lever:.3f}: the inertia error is about {2 * lever:.2f} x the period error")
    for tol in (0.005, 0.01, 0.02):
        a_T = max_amplitude_for_period_error(tol)
        a_I = max_amplitude_for_inertia_error(tol, case)
        print(f"  {100 * tol:.1f} % error: period allows {np.degrees(a_T):.1f} deg, inferred I_cg allows {np.degrees(a_I):.1f} deg")

    # RK4 against the closed form (linear) and against solve_ivp (nonlinear), with convergence.
    theta0 = np.radians(20.0)
    t_end = 3.0 * period_small_angle(w0)
    print("\nRK4 convergence on the linearised pendulum (max |theta - closed form| over three periods):")
    for h in (0.02, 0.01, 0.005, 0.0025):
        t, Y = integrate_rk4(lambda tt, yy: linear_pendulum_rhs(tt, yy, w0, zeta), [theta0, 0.0], t_end, h)
        err = np.max(np.abs(Y[:, 0] - linear_response(t, theta0, w0, zeta)))
        print(f"  h = {h:.4f} s: error {err:.3e}")
    f_nl = lambda tt, yy: pendulum_rhs(tt, yy, w0, 0.0)
    t, Y = integrate_rk4(f_nl, [theta0, 0.0], 5.0 * period_small_angle(w0), 0.002)
    sol = solve_ivp(f_nl, (0.0, t[-1]), [theta0, 0.0], t_eval=t, rtol=1e-10, atol=1e-12)
    E = energy(Y[:, 0], Y[:, 1], w0)
    print(f"RK4 vs solve_ivp (undamped, 20 deg, h = 2 ms): max |dtheta| = {np.max(np.abs(Y[:, 0] - sol.y[0])):.2e} rad;"
          f" energy drift {np.max(np.abs(E - E[0])) / E[0]:.2e};"
          f" period from zero crossings {period_from_trajectory(t, Y[:, 0]):.5f} s, exact {period_exact(theta0, w0):.5f} s")
    print(f"Van der Pol (mu = 1) limit-cycle amplitude from (0.1, 0): {limit_cycle_amplitude(1.0, [0.1, 0.0]):.4f},"
          f" from (4, 0): {limit_cycle_amplitude(1.0, [4.0, 0.0]):.4f}")

    # Figure: phase portrait with separatrix, and the inertia error against amplitude.
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    th = np.linspace(-2.0 * np.pi, 2.0 * np.pi, 41)
    om = np.linspace(-2.5 * w0, 2.5 * w0, 25)
    TH, OM = np.meshgrid(th, om)
    U = OM
    V = -w0 ** 2 * np.sin(TH)
    ax1.streamplot(TH, OM, U, V, density=1.0, linewidth=0.6, color="0.6")
    for E0 in (0.25, 1.0, 2.0, 3.0):
        Es = E0 * w0 ** 2
        ths = np.linspace(-2.0 * np.pi, 2.0 * np.pi, 1001)
        arg = 2.0 * (Es - w0 ** 2 * (1.0 - np.cos(ths)))
        om_c = np.sqrt(np.clip(arg, 0.0, None))
        om_c[arg < 0.0] = np.nan
        style = "r-" if E0 == 2.0 else "b-"
        ax1.plot(ths, om_c, style, lw=1.2)
        ax1.plot(ths, -om_c, style, lw=1.2)
    ax1.plot([0.0, 2 * np.pi, -2 * np.pi], [0.0, 0.0, 0.0], "ko", ms=5)
    ax1.plot([np.pi, -np.pi], [0.0, 0.0], "kx", ms=7, mew=2)
    ax1.set_xlabel("theta [rad]")
    ax1.set_ylabel("theta_dot [rad/s]")
    ax1.set_title("undamped pendulum: centers (o), saddles (x), separatrix (red)")
    degs = np.linspace(1.0, 60.0, 120)
    errs = [100.0 * swing_test_error(np.radians(dd), case)["rel_error"] for dd in degs]
    ax2.plot(degs, errs, "b-", label="exact period in the small-angle formula")
    lever = (case["I_cg"] + case["m"] * case["d"] ** 2) / case["I_cg"]
    ax2.plot(degs, 100.0 * lever * ((1.0 + np.radians(degs) ** 2 / 16.0) ** 2 - 1.0), "k--",
             label=r"$(I_p/I_{cg})\,[(1 + \theta_0^2/16)^2 - 1]$")
    ax2.axhline(1.0, color="0.5", lw=0.8)
    ax2.set_xlabel("release amplitude [deg]")
    ax2.set_ylabel("error in inferred I_cg [%]")
    ax2.set_ylim(0, 55)
    ax2.grid(True, lw=0.4)
    ax2.legend()
    ax2.set_title("swing-test error from the linearisation")
    fig.tight_layout()
    fig.savefig(out_dir / "pendulum_linearization.png", dpi=150)
    print(f"saved {out_dir / 'pendulum_linearization.png'}")


if __name__ == "__main__":
    main()
