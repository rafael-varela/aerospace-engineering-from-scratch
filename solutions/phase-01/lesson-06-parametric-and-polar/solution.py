"""P01.06 Parametric Curves & Polar Coordinates: reference solution.

A conic orbit in polar form about the attracting focus,

    r(theta) = p / (1 + e cos theta),        p = semi-latus rectum [m], e = eccentricity [-],

is traced, differentiated and integrated: arc length and swept area from the polar formulas,
time of flight from Kepler's second law (dA/dt = h/2, h = sqrt(mu p)), and the speed along the
curve from the parametric derivatives (dx/dt, dy/dt), checked against the vis-viva equation.

Run:  python solution.py          (prints the orbit table, writes conic_orbits.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-06-parametric-and-polar
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

MU_EARTH = 3.986004418e14  # m**3/s**2, Earth's gravitational parameter G*M


def _check_conic(p: float, e: float) -> None:
    if not np.isfinite(p) or p <= 0.0:
        raise ValueError(f"semi-latus rectum p must be positive and finite, got {p!r}")
    if not np.isfinite(e) or e < 0.0:
        raise ValueError(f"eccentricity e must be non-negative and finite, got {e!r}")


def true_anomaly_limit(e: float) -> float:
    """Return the largest true anomaly (rad) at which r(theta) is finite and positive.

    For e < 1 every angle is allowed and the limit is pi (a closed curve). For e >= 1 the
    denominator 1 + e cos(theta) vanishes at theta_inf = arccos(-1/e) and the branch is open.
    """
    _check_conic(1.0, e)
    if e < 1.0:
        return float(np.pi)
    return float(np.arccos(-1.0 / e))


def conic_radius(theta: np.ndarray | float, p: float, e: float) -> np.ndarray | float:
    """Return r = p / (1 + e cos theta) in metres, vectorized over ``theta`` (rad).

    Raise ``ValueError`` if p <= 0, e < 0, or any theta lies on or beyond the asymptote of an
    open conic (1 + e cos theta <= 0).
    """
    _check_conic(p, e)
    th = np.asarray(theta, dtype=float)
    denom = 1.0 + e * np.cos(th)
    if np.any(denom <= 0.0):
        raise ValueError("theta lies on or beyond the asymptote of an open conic")
    r = p / denom
    return float(r) if r.ndim == 0 else r


def conic_radius_derivative(theta: np.ndarray | float, p: float, e: float) -> np.ndarray | float:
    """Return dr/dtheta = p e sin(theta) / (1 + e cos theta)**2 in metres per radian."""
    _check_conic(p, e)
    th = np.asarray(theta, dtype=float)
    denom = 1.0 + e * np.cos(th)
    if np.any(denom <= 0.0):
        raise ValueError("theta lies on or beyond the asymptote of an open conic")
    d = p * e * np.sin(th) / denom**2
    return float(d) if d.ndim == 0 else d


def polar_to_cartesian(r: np.ndarray | float, theta: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Return (x, y) = (r cos theta, r sin theta) with the focus at the origin."""
    r = np.asarray(r, dtype=float)
    th = np.asarray(theta, dtype=float)
    return r * np.cos(th), r * np.sin(th)


def periapsis_apoapsis(p: float, e: float) -> tuple[float, float]:
    """Return (r_p, r_a) in metres: p/(1+e) and p/(1-e); r_a is ``inf`` for e >= 1."""
    _check_conic(p, e)
    rp = p / (1.0 + e)
    ra = p / (1.0 - e) if e < 1.0 else float("inf")
    return float(rp), float(ra)


def semi_major_axis(p: float, e: float) -> float:
    """Return a = p / (1 - e**2): positive for an ellipse, ``inf`` for a parabola, negative for a hyperbola."""
    _check_conic(p, e)
    if e == 1.0:
        return float("inf")
    return float(p / (1.0 - e**2))


def conic_from_apsides(r_p: float, r_a: float) -> tuple[float, float]:
    """Return (p, e) of the ellipse with periapsis ``r_p`` and apoapsis ``r_a`` (both in metres).

    e = (r_a - r_p)/(r_a + r_p) and p = r_p (1 + e). Raise ``ValueError`` unless 0 < r_p <= r_a.
    """
    if not (0.0 < r_p <= r_a):
        raise ValueError("need 0 < r_p <= r_a")
    e = (r_a - r_p) / (r_a + r_p)
    return float(r_p * (1.0 + e)), float(e)


def orbit_points(p: float, e: float, n: int = 361, theta_max: float | None = None) -> tuple[np.ndarray, np.ndarray]:
    """Return arrays (x, y) of ``n`` points on the conic for theta in [-theta_max, theta_max].

    ``theta_max`` defaults to the full closed curve (pi) for e < 1 and to 95 % of the asymptote
    angle for e >= 1 so the open branch stays finite.
    """
    _check_conic(p, e)
    if n < 2:
        raise ValueError("n must be at least 2")
    if theta_max is None:
        theta_max = np.pi if e < 1.0 else 0.95 * true_anomaly_limit(e)
    th = np.linspace(-theta_max, theta_max, n)
    return polar_to_cartesian(conic_radius(th, p, e), th)


def simpson(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int) -> float:
    """Composite Simpson's rule with ``n`` (even, >= 2) panels on [a, b]. Raise ``ValueError`` otherwise."""
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer >= 2")
    x = np.linspace(a, b, n + 1)
    y = np.asarray(f(x), dtype=float)
    h = (b - a) / n
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))


def polar_arc_length(p: float, e: float, theta0: float, theta1: float, n: int = 400) -> float:
    """Return the arc length of the conic between ``theta0`` and ``theta1`` (rad), metres.

    L = integral of sqrt(r**2 + (dr/dtheta)**2) dtheta, evaluated with composite Simpson.
    """
    def integrand(th: np.ndarray) -> np.ndarray:
        r = conic_radius(th, p, e)
        return np.sqrt(r**2 + conic_radius_derivative(th, p, e) ** 2)

    return simpson(integrand, theta0, theta1, n)


def sector_area(p: float, e: float, theta0: float, theta1: float, n: int = 400) -> float:
    """Return the area swept by the radius vector between ``theta0`` and ``theta1``: 1/2 integral r**2 dtheta (m**2)."""
    return simpson(lambda th: 0.5 * conic_radius(th, p, e) ** 2, theta0, theta1, n)


def specific_angular_momentum(p: float, mu: float = MU_EARTH) -> float:
    """Return h = sqrt(mu p) in m**2/s, the constant r**2 dtheta/dt of the orbit."""
    _check_conic(p, 0.0)
    if mu <= 0.0:
        raise ValueError("mu must be positive")
    return float(np.sqrt(mu * p))


def time_of_flight(p: float, e: float, theta0: float, theta1: float, mu: float = MU_EARTH, n: int = 400) -> float:
    """Return the time (s) to move from ``theta0`` to ``theta1`` along the conic.

    Kepler's second law: dA/dt = h/2, so dt = (r**2 / h) dtheta and t = (2/h) * swept area.
    """
    return 2.0 * sector_area(p, e, theta0, theta1, n) / specific_angular_momentum(p, mu)


def orbital_period(p: float, e: float, mu: float = MU_EARTH) -> float:
    """Return the period T = 2 pi sqrt(a**3 / mu) of a closed orbit (e < 1); raise ``ValueError`` otherwise."""
    _check_conic(p, e)
    if e >= 1.0:
        raise ValueError("only a closed orbit (e < 1) has a period")
    a = semi_major_axis(p, e)
    return float(2.0 * np.pi * np.sqrt(a**3 / mu))


def orbital_speed(theta: np.ndarray | float, p: float, e: float, mu: float = MU_EARTH) -> np.ndarray | float:
    """Return the speed |(dx/dt, dy/dt)| along the conic at true anomaly ``theta`` (m/s).

    From the polar velocity components r' = (h/p) e sin(theta) and r theta' = (h/p)(1 + e cos theta):
    v = (h/p) sqrt(1 + 2 e cos theta + e**2).
    """
    _check_conic(p, e)
    th = np.asarray(theta, dtype=float)
    v = specific_angular_momentum(p, mu) / p * np.sqrt(1.0 + 2.0 * e * np.cos(th) + e**2)
    return float(v) if v.ndim == 0 else v


def vis_viva_speed(r: np.ndarray | float, a: float, mu: float = MU_EARTH) -> np.ndarray | float:
    """Return the speed from energy conservation, v = sqrt(mu (2/r - 1/a)); use a = inf for a parabola."""
    r = np.asarray(r, dtype=float)
    inv_a = 0.0 if np.isinf(a) else 1.0 / a
    v = np.sqrt(mu * (2.0 / r - inv_a))
    return float(v) if v.ndim == 0 else v


def flight_path_angle(theta: np.ndarray | float, e: float) -> np.ndarray | float:
    """Return the flight-path angle gamma (rad) between the velocity and the local horizontal.

    tan(gamma) = r' / (r theta') = e sin(theta) / (1 + e cos theta); gamma = 0 at the apsides.
    """
    _check_conic(1.0, e)
    th = np.asarray(theta, dtype=float)
    g = np.arctan2(e * np.sin(th), 1.0 + e * np.cos(th))
    return float(g) if g.ndim == 0 else g


def ellipse_parametric(a: float, e: float, E: np.ndarray | float) -> tuple[np.ndarray, np.ndarray]:
    """Return focus-centred (x, y) of the ellipse from the eccentric anomaly ``E`` (rad).

    x = a (cos E - e), y = a sqrt(1 - e**2) sin E. Raise ``ValueError`` unless a > 0 and 0 <= e < 1.
    """
    if not (a > 0.0) or not (0.0 <= e < 1.0):
        raise ValueError("ellipse needs a > 0 and 0 <= e < 1")
    E = np.asarray(E, dtype=float)
    return a * (np.cos(E) - e), a * np.sqrt(1.0 - e**2) * np.sin(E)


def arc_length_convergence(p: float, e: float, ns: list[int], theta1: float = np.pi / 2) -> np.ndarray:
    """Return |Simpson arc length(n) - reference| over [0, theta1] for each ``n`` in ``ns``.

    The reference is ``scipy.integrate.quad`` at tight tolerance; Simpson should converge at fourth order.
    """
    from scipy.integrate import quad

    def integrand(th: float) -> float:
        r = conic_radius(th, p, e)
        return float(np.sqrt(r**2 + conic_radius_derivative(th, p, e) ** 2))

    ref, _ = quad(integrand, 0.0, theta1, epsabs=0.0, epsrel=1e-13, limit=200)
    return np.array([abs(polar_arc_length(p, e, 0.0, theta1, n) - ref) for n in ns])


def main(out_dir: str | Path = ".") -> None:
    """Print the transfer-orbit table and save conic_orbits.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import ellipe

    out_dir = Path(out_dir)
    r_earth = 6378.137e3
    r_p, r_a = r_earth + 300e3, 42164e3               # LEO parking orbit to geostationary radius
    p, e = conic_from_apsides(r_p, r_a)
    a = semi_major_axis(p, e)
    b = a * np.sqrt(1.0 - e**2)
    h = specific_angular_momentum(p)
    print("Geostationary transfer ellipse (focus at Earth's centre)")
    print(f"  r_p = {r_p/1e3:.1f} km, r_a = {r_a/1e3:.1f} km  ->  e = {e:.5f}, p = {p/1e3:.1f} km, a = {a/1e3:.1f} km")
    print(f"  h = sqrt(mu p) = {h:.4e} m^2/s")
    for name, th in (("periapsis", 0.0), ("theta = 90 deg", np.pi / 2), ("apoapsis", np.pi)):
        r = conic_radius(th, p, e)
        print(f"  {name:14s} r = {r/1e3:9.1f} km, v(parametric) = {orbital_speed(th, p, e):8.1f} m/s,"
              f" v(vis-viva) = {vis_viva_speed(r, a):8.1f} m/s, gamma = {np.degrees(flight_path_angle(th, e)):6.2f} deg")
    t_half = time_of_flight(p, e, 0.0, np.pi)
    T = orbital_period(p, e)
    print(f"  time periapsis -> apoapsis (Simpson, n=400) = {t_half:.2f} s; Kepler III period / 2 = {T/2:.2f} s")
    area_full = sector_area(p, e, -np.pi, np.pi)
    print(f"  full swept area = {area_full:.6e} m^2; pi a b = {np.pi*a*b:.6e} m^2")
    per = polar_arc_length(p, e, -np.pi, np.pi, 400)
    print(f"  perimeter (Simpson) = {per/1e3:.3f} km; 4 a E(e^2) (scipy.special.ellipe) = {4*a*ellipe(e**2)/1e3:.3f} km")
    ns = [8, 16, 32, 64, 128]
    errs = arc_length_convergence(p, e, ns)
    print("  Simpson arc-length error over [0, pi/2]:", " ".join(f"n={n}: {er:.2e}" for n, er in zip(ns, errs)))
    print("  observed orders:", np.round(np.log2(errs[:-1] / errs[1:]), 2))

    fig, ax = plt.subplots(figsize=(6, 6))
    for ecc, label in ((0.0, "circle e = 0"), (e, f"ellipse e = {e:.3f}"), (1.0, "parabola e = 1"), (1.5, "hyperbola e = 1.5")):
        x, y = orbit_points(p, ecc, 721)
        ax.plot(x / 1e3, y / 1e3, label=label)
    ax.plot(0.0, 0.0, "ko", label="focus (Earth's centre)")
    ax.add_patch(plt.Circle((0.0, 0.0), r_earth / 1e3, color="0.8"))
    ax.set_aspect("equal")
    ax.set_xlim(-60e3, 25e3)
    ax.set_ylim(-40e3, 40e3)
    ax.set_xlabel("x [km]")
    ax.set_ylabel("y [km]")
    ax.set_title(f"Conics with the same semi-latus rectum p = {p/1e3:.0f} km")
    ax.grid(True)
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "conic_orbits.png", dpi=150)
    print(f"saved {out_dir / 'conic_orbits.png'}")


if __name__ == "__main__":
    main()
