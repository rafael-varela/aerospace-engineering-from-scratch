"""P03.02 Root Finding: reference solution.

Bisection, Newton, secant and fixed-point iteration written from scratch, each returning its whole
iterate history so the convergence order can be measured; applied to Kepler's equation
M = E - e sin E (the position of a satellite on an elliptical orbit at a given time) and to the
Rayleigh pitot formula (Mach number from a supersonic pitot pressure ratio). ``scipy.optimize`` is
the oracle.

Conventions:
    * Every solver returns a ``RootResult(root, iterates, n_evals)``: ``iterates`` is the array of
      successive approximations x_0, x_1, ..., x_k (the last entry is ``root``); ``n_evals`` counts
      calls of f (and of f' for Newton).
    * Stopping: |x_{k+1} - x_k| <= tol (Newton, secant, fixed point); half-interval <= tol (bisection).
      ``max_iter`` exceeded raises ``RuntimeError``; a zero derivative raises ``ZeroDivisionError``;
      a bracket without a sign change raises ``ValueError``.
    * Angles in radians. Kepler: M mean anomaly, E eccentric anomaly, nu true anomaly, 0 <= e < 1.
      The solver reduces M to [0, 2 pi) and adds the revolutions back to E.
    * Orbit state is given in the perifocal frame (x toward periapsis, y along the velocity at
      periapsis), SI units: metres, seconds, mu in m^3/s^2.

Run:  python solution.py          (prints the convergence tables, writes root_finding.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-03/lesson-02-root-finding
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np

MU_EARTH: float = 3.986004418e14        # m^3/s^2
R_EARTH: float = 6.378137e6             # m


class RootResult(NamedTuple):
    """Result of a scalar root finder: the root, every iterate on the way, and the function-call count."""

    root: float
    iterates: np.ndarray
    n_evals: int


# ---------------------------------------------------------------- the four iterations

def bisection(f: Callable[[float], float], a: float, b: float, tol: float = 1e-12, max_iter: int = 200) -> RootResult:
    """Bracketing root finder: halve [a, b] keeping the sign change, until the half-width is <= tol.

    ``iterates`` holds the successive midpoints. The error after k steps is at most (b - a) / 2^(k+1):
    linear convergence with rate exactly 1/2, one function evaluation per step, guaranteed for any
    continuous f with f(a) f(b) < 0. Raise ``ValueError`` if there is no sign change, ``RuntimeError``
    if ``max_iter`` steps do not reach ``tol``.
    """
    fa, fb = f(a), f(b)
    n_evals = 2
    if fa == 0.0:
        return RootResult(float(a), np.array([float(a)]), n_evals)
    if fb == 0.0:
        return RootResult(float(b), np.array([float(b)]), n_evals)
    if fa * fb > 0.0:
        raise ValueError("bisection needs f(a) and f(b) of opposite sign")
    mids: list[float] = []
    for _ in range(max_iter):
        m = 0.5 * (a + b)
        mids.append(m)
        fm = f(m)
        n_evals += 1
        if fm == 0.0 or 0.5 * (b - a) <= tol:
            return RootResult(m, np.array(mids), n_evals)
        if fa * fm < 0.0:
            b, fb = m, fm
        else:
            a, fa = m, fm
    raise RuntimeError(f"bisection did not reach tol={tol} in {max_iter} steps")


def bisection_steps_needed(a: float, b: float, tol: float) -> int:
    """Return the number of midpoints k after which the bracket half-width (b - a) / 2^k is <= tol: ceil(log2((b - a) / tol)), at least 1.

    The k-th midpoint is within tol of the root, so this is the number of f evaluations (beyond the
    two endpoint checks) that ``bisection`` makes. Raise ``ValueError`` if b <= a or tol <= 0.
    """
    if tol <= 0.0 or b <= a:
        raise ValueError("need b > a and tol > 0")
    return max(1, math.ceil(math.log2((b - a) / tol)))


def newton(f: Callable[[float], float], dfdx: Callable[[float], float], x0: float, tol: float = 1e-12,
           max_iter: int = 50) -> RootResult:
    """Newton's method x_{k+1} = x_k - f(x_k) / f'(x_k), stopped when |x_{k+1} - x_k| <= tol.

    Quadratic convergence near a simple root (error squares each step); linear with rate 1/2 at a
    double root; may diverge or cycle from a poor start, which surfaces as ``RuntimeError`` after
    ``max_iter`` steps. Raise ``ZeroDivisionError`` if f'(x_k) = 0. Counts f and f' calls in n_evals.
    """
    x = float(x0)
    xs = [x]
    n_evals = 0
    for _ in range(max_iter):
        fx, dfx = f(x), dfdx(x)
        n_evals += 2
        if dfx == 0.0:
            raise ZeroDivisionError("Newton hit a zero derivative")
        x_new = x - fx / dfx
        if not math.isfinite(x_new) or abs(x_new) > 1e150:
            raise RuntimeError("Newton diverged (iterate not finite or beyond 1e150)")
        xs.append(x_new)
        if abs(x_new - x) <= tol:
            return RootResult(x_new, np.array(xs), n_evals)
        x = x_new
    raise RuntimeError(f"Newton did not converge in {max_iter} steps")


def secant(f: Callable[[float], float], x0: float, x1: float, tol: float = 1e-12, max_iter: int = 100) -> RootResult:
    """Secant method: Newton with f' replaced by the slope through the last two iterates.

    Order (1 + sqrt 5) / 2 = 1.618 at a simple root, one new f evaluation per step. Raise
    ``ZeroDivisionError`` if f(x_k) = f(x_{k-1}) before convergence, ``RuntimeError`` on ``max_iter``.
    """
    x_prev, x = float(x0), float(x1)
    f_prev, fx = f(x_prev), f(x)
    n_evals = 2
    xs = [x_prev, x]
    for _ in range(max_iter):
        if fx == f_prev:
            if fx == 0.0:
                return RootResult(x, np.array(xs), n_evals)
            raise ZeroDivisionError("secant slope is zero")
        x_new = x - fx * (x - x_prev) / (fx - f_prev)
        if not math.isfinite(x_new):
            raise RuntimeError("secant iterate is not finite")
        xs.append(x_new)
        if abs(x_new - x) <= tol:
            return RootResult(x_new, np.array(xs), n_evals)
        x_prev, f_prev = x, fx
        x, fx = x_new, f(x_new)
        n_evals += 1
    raise RuntimeError(f"secant did not converge in {max_iter} steps")


def fixed_point(g: Callable[[float], float], x0: float, tol: float = 1e-12, max_iter: int = 10_000) -> RootResult:
    """Iterate x_{k+1} = g(x_k) until |x_{k+1} - x_k| <= tol.

    Converges to a fixed point x* = g(x*) when |g'(x*)| < 1, linearly with rate |g'(x*)|; diverges
    when |g'(x*)| > 1 (``RuntimeError`` after ``max_iter`` or when an iterate is not finite).
    """
    x = float(x0)
    xs = [x]
    for k in range(max_iter):
        x_new = g(x)
        if not math.isfinite(x_new):
            raise RuntimeError("fixed-point iterate is not finite")
        xs.append(x_new)
        if abs(x_new - x) <= tol:
            return RootResult(x_new, np.array(xs), k + 1)
        x = x_new
    raise RuntimeError(f"fixed-point iteration did not converge in {max_iter} steps")


# ---------------------------------------------------------------- measuring convergence

def errors_from_iterates(iterates: np.ndarray, root: float) -> np.ndarray:
    """Return |x_k - root| for every iterate, dropping entries that are exactly zero (rounding floor)."""
    e = np.abs(np.asarray(iterates, dtype=float) - root)
    return e[e > 0.0]


def observed_convergence_order(errors: np.ndarray, floor: float = 1e-13) -> np.ndarray:
    """Return the per-step estimates q_k = log(e_{k+1}/e_k) / log(e_k/e_{k-1}) of the order in e_{k+1} ~ C e_k^q.

    Only triples of consecutive errors all above ``floor`` and strictly decreasing are used (below the
    floor rounding, not the method, sets the error). Returns an empty array if no triple qualifies.
    """
    e = np.asarray(errors, dtype=float)
    out = []
    for k in range(1, e.size - 1):
        if e[k - 1] > floor and e[k] > floor and e[k + 1] > floor and e[k + 1] < e[k] < e[k - 1]:
            out.append(math.log(e[k + 1] / e[k]) / math.log(e[k] / e[k - 1]))
    return np.array(out)


def fitted_convergence_order(errors: np.ndarray, floor: float = 1e-13) -> float:
    """Return the least-squares slope of log e_{k+1} against log e_k over all consecutive pairs above ``floor``.

    A steadier estimate of q than the per-step ratios of ``observed_convergence_order``: for the
    secant method the per-step estimates oscillate around 1.618 while this fit lands close to it.
    Raise ``ValueError`` with fewer than three usable pairs.
    """
    e = np.asarray(errors, dtype=float)
    e = e[e > floor]
    if e.size < 4:
        raise ValueError("need at least four errors above the floor")
    slope, _ = np.polyfit(np.log(e[:-1]), np.log(e[1:]), 1)
    return float(slope)


def observed_linear_rate(errors: np.ndarray, floor: float = 1e-13) -> float:
    """Return the ratio e_{k+1}/e_k averaged over the last few steps above ``floor`` (the linear rate)."""
    e = np.asarray(errors, dtype=float)
    e = e[e > floor]
    if e.size < 2:
        raise ValueError("need at least two errors above the floor")
    ratios = e[1:] / e[:-1]
    return float(np.mean(ratios[-5:]))


# ---------------------------------------------------------------- Kepler's equation

def kepler_residual(E: float, M: float, e: float) -> float:
    """Return f(E) = E - e sin E - M, whose root is the eccentric anomaly."""
    return E - e * math.sin(E) - M


def kepler_derivative(E: float, e: float) -> float:
    """Return f'(E) = 1 - e cos E, which lies in [1 - e, 1 + e] and never vanishes for e < 1."""
    return 1.0 - e * math.cos(E)


def solve_kepler(M: float, e: float, method: str = "newton", tol: float = 1e-14) -> RootResult:
    """Solve Kepler's equation for E given M (rad) and e in [0, 1) by 'newton', 'secant', 'bisection' or 'fixed_point'.

    M is reduced to [0, 2 pi) and the whole revolutions are added back to E. Starting values:
    Newton and fixed point from E0 = M + e sin M (the first-order fixed-point step), secant from
    (M, M + e), bisection on [M - e, M + e] (f is monotone, so this bracket always works).
    Raise ``ValueError`` for e outside [0, 1) or an unknown method.
    """
    if not 0.0 <= e < 1.0:
        raise ValueError("eccentricity must satisfy 0 <= e < 1 for an ellipse")
    revs = math.floor(M / (2.0 * math.pi))
    m = M - 2.0 * math.pi * revs
    f = lambda E: kepler_residual(E, m, e)                       # noqa: E731
    df = lambda E: kepler_derivative(E, e)                       # noqa: E731
    E0 = m + e * math.sin(m)
    if method == "newton":
        res = newton(f, df, E0, tol=tol)
    elif method == "secant":
        res = secant(f, m, m + max(e, 1e-3), tol=tol)
    elif method == "bisection":
        res = bisection(f, m - e, m + e, tol=tol)
    elif method == "fixed_point":
        res = fixed_point(lambda E: m + e * math.sin(E), E0, tol=tol)
    else:
        raise ValueError(f"unknown method {method!r}")
    shift = 2.0 * math.pi * revs
    return RootResult(res.root + shift, res.iterates + shift, res.n_evals)


def true_anomaly(E: float, e: float) -> float:
    """Return the true anomaly nu from E: tan(nu/2) = sqrt((1 + e)/(1 - e)) tan(E/2), via atan2 (no quadrant trouble)."""
    return 2.0 * math.atan2(math.sqrt(1.0 + e) * math.sin(E / 2.0), math.sqrt(1.0 - e) * math.cos(E / 2.0))


def mean_motion(a: float, mu: float = MU_EARTH) -> float:
    """Return n = sqrt(mu / a^3) [rad/s]."""
    if a <= 0.0:
        raise ValueError("semi-major axis must be positive")
    return math.sqrt(mu / a ** 3)


def orbit_state(a: float, e: float, t: float, mu: float = MU_EARTH, M0: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Return position [m] and velocity [m/s] in the perifocal frame at time t after M = M0.

    M = M0 + n t; E from ``solve_kepler``; r = a (cos E - e, sqrt(1 - e^2) sin E);
    v = (sqrt(mu a) / |r|) (-sin E, sqrt(1 - e^2) cos E). Satisfies vis-viva v^2 = mu (2/|r| - 1/a)
    and |r x v| = sqrt(mu a (1 - e^2)) exactly (to rounding).
    """
    n = mean_motion(a, mu)
    E = solve_kepler(M0 + n * t, e).root
    b_over_a = math.sqrt(1.0 - e * e)
    r = a * np.array([math.cos(E) - e, b_over_a * math.sin(E)])
    rnorm = a * (1.0 - e * math.cos(E))
    v = (math.sqrt(mu * a) / rnorm) * np.array([-math.sin(E), b_over_a * math.cos(E)])
    return r, v


def specific_energy(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> float:
    """Return v^2/2 - mu/|r| [J/kg]; equals -mu / (2a) on every point of the ellipse."""
    return 0.5 * float(np.dot(v, v)) - mu / float(np.linalg.norm(r))


# ---------------------------------------------------------------- a second application: supersonic pitot

def rayleigh_pitot_ratio(M: float, gamma: float = 1.4) -> float:
    """Return p02 / p1 for a pitot tube behind a normal shock at Mach M > 1 (Rayleigh pitot formula).

    p02/p1 = [((gamma+1)^2 M^2 / (4 gamma M^2 - 2 (gamma-1)))]^(gamma/(gamma-1)) * (1 - gamma + 2 gamma M^2)/(gamma+1).
    Monotone increasing in M; equals 1.8929 at M = 1 for gamma = 1.4. Raise ``ValueError`` for M < 1.
    """
    if M < 1.0:
        raise ValueError("the Rayleigh pitot formula needs M >= 1 (a shock in front of the probe)")
    g = gamma
    term1 = ((g + 1.0) ** 2 * M * M / (4.0 * g * M * M - 2.0 * (g - 1.0))) ** (g / (g - 1.0))
    term2 = (1.0 - g + 2.0 * g * M * M) / (g + 1.0)
    return term1 * term2


def mach_from_pitot_ratio(ratio: float, gamma: float = 1.4, tol: float = 1e-12) -> RootResult:
    """Invert the Rayleigh pitot formula for M by bisection on [1, 20] followed by Newton polish.

    There is no closed form, so this is a root-finding problem in flight-test data reduction.
    Raise ``ValueError`` if ``ratio`` is below the M = 1 value (the flow is subsonic: use the isentropic formula).
    """
    r1 = rayleigh_pitot_ratio(1.0, gamma)
    if ratio < r1:
        raise ValueError("pressure ratio below the M = 1 value: subsonic, no shock")
    f = lambda M: rayleigh_pitot_ratio(M, gamma) - ratio      # noqa: E731
    coarse = bisection(f, 1.0, 20.0, tol=1e-3)
    df = lambda M: (f(M + 1e-6) - f(M - 1e-6)) / 2e-6          # noqa: E731
    fine = newton(f, df, coarse.root, tol=tol)
    return RootResult(fine.root, np.concatenate([coarse.iterates, fine.iterates[1:]]), coarse.n_evals + fine.n_evals)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Solve Kepler's equation four ways, measure each method's order, place a GTO satellite, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import optimize

    out_dir = Path(out_dir)
    e, M = 0.73, math.radians(60.0)                           # a geostationary transfer orbit, one hour-ish past perigee
    ref = optimize.brentq(lambda E: kepler_residual(E, M, e), M - e, M + e, xtol=1e-15)
    print(f"Kepler's equation, e = {e}, M = {math.degrees(M):.0f} deg: scipy brentq E = {ref:.15f} rad")
    print(f"{'method':<12}{'E (rad)':>20}{'|f(E)|':>12}{'steps':>7}{'f evals':>9}{'order estimates':>40}")
    for method in ("bisection", "fixed_point", "secant", "newton"):
        res = solve_kepler(M, e, method)
        errs = errors_from_iterates(res.iterates, ref)
        q = observed_convergence_order(errs)
        qs = ", ".join(f"{v:.2f}" for v in q[-3:]) if q.size else "-"
        if method == "bisection":
            qs = "n/a: bound halves each step (rate 0.5)"
        print(f"{method:<12}{res.root:>20.15f}{abs(kepler_residual(res.root, M, e)):>12.1e}"
              f"{res.iterates.size - 1:>7}{res.n_evals:>9}{qs:>40}")
    fp = solve_kepler(M, e, "fixed_point")
    rate = observed_linear_rate(errors_from_iterates(fp.iterates, ref))
    print(f"fixed-point linear rate observed {rate:.4f}, predicted |g'(E*)| = e |cos E*| = {e * abs(math.cos(ref)):.4f}")
    bs = solve_kepler(M, e, "bisection")
    print(f"bisection: {bs.iterates.size} midpoints for tol 1e-14 on a bracket of width {2 * e:.2f}; "
          f"predicted {bisection_steps_needed(M - e, M + e, 1e-14)}")

    # --- place the satellite
    a = (R_EARTH + 250e3 + R_EARTH + 35_786e3) / 2.0             # GTO: 250 km x 35 786 km
    e_gto = 1.0 - (R_EARTH + 250e3) / a
    n = mean_motion(a)
    T = 2.0 * math.pi / n
    print(f"\nGTO: a = {a / 1e3:.1f} km, e = {e_gto:.4f}, period {T / 3600:.3f} h")
    for t in (0.0, 1800.0, 3600.0, T / 2, 3 * T / 4):
        r, v = orbit_state(a, e_gto, t)
        E = solve_kepler(n * t, e_gto).root
        nu = true_anomaly(E, e_gto)
        eps = specific_energy(r, v)
        print(f"  t = {t:8.0f} s: M = {math.degrees(n * t):7.2f} deg, E = {math.degrees(E):7.2f} deg, nu = {math.degrees(nu) % 360:7.2f} deg, "
              f"|r| = {np.linalg.norm(r) / 1e3:9.1f} km, |v| = {np.linalg.norm(v):8.1f} m/s, energy {eps / 1e6:.4f} MJ/kg "
              f"(-mu/2a = {-MU_EARTH / (2 * a) / 1e6:.4f})")

    # --- Newton versus eccentricity near e -> 1
    print("\nNewton steps from E0 = M + e sin M to 1e-14, M = 5 deg:")
    for ecc in (0.1, 0.5, 0.9, 0.99, 0.999):
        res = solve_kepler(math.radians(5.0), ecc)
        print(f"  e = {ecc:<6} steps = {res.iterates.size - 1}")

    # --- pitot
    ratio = rayleigh_pitot_ratio(2.0)
    inv = mach_from_pitot_ratio(ratio)
    print(f"\nRayleigh pitot: p02/p1 at M = 2.0 is {ratio:.5f}; inverted M = {inv.root:.12f} in {inv.iterates.size - 1} steps")

    # --- Newton failure
    try:
        newton(math.atan, lambda x: 1.0 / (1.0 + x * x), 1.5)
    except (RuntimeError, ZeroDivisionError) as exc:
        print(f"Newton on atan from x0 = 1.5: {exc}")
    dbl = newton(lambda x: (x - 1.0) ** 2, lambda x: 2.0 * (x - 1.0), 2.0, tol=1e-10)
    print(f"Newton on a double root (x-1)^2 from 2: {dbl.iterates.size - 1} steps, rate "
          f"{observed_linear_rate(errors_from_iterates(dbl.iterates, 1.0), floor=1e-8):.3f} (expected 0.5)")

    # --- figure
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    for method, marker in (("bisection", "."), ("fixed_point", "x"), ("secant", "s"), ("newton", "o")):
        errs = np.abs(solve_kepler(M, e, method).iterates - ref)
        ax.semilogy(np.arange(errs.size), np.maximum(errs, 1e-17), marker=marker, label=method)
    ax.set_xlabel("iteration k")
    ax.set_ylabel("|E_k - E*| [rad]")
    ax.set_title(f"Kepler's equation, e = {e}, M = {math.degrees(M):.0f}°")
    ax.set_xlim(0, 40)
    ax.legend()
    ax.grid(True, which="both", alpha=0.3)
    ax = axes[1]
    Ms = np.linspace(0.0, 2.0 * math.pi, 361)
    for ecc in (0.0, 0.3, 0.73, 0.95):
        Es = np.array([solve_kepler(m, ecc).root for m in Ms])
        ax.plot(np.degrees(Ms), np.degrees(Es), label=f"e = {ecc}")
    ax.set_xlabel("mean anomaly M [deg]")
    ax.set_ylabel("eccentric anomaly E [deg]")
    ax.set_title("E(M) from Newton")
    ax.legend()
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "root_finding.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main()
