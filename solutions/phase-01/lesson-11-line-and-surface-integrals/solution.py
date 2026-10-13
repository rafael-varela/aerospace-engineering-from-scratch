"""P01.11 Vector Fields, Line & Surface Integrals: reference solution.

Circulation around a vortex, work in a gravitational field and mass flow through a surface,
all computed as line and surface integrals of vector fields and checked against closed forms.

Conventions (hold them in your head, they are the source of every sign error):
    * 2D fields are callables ``field(x, y) -> (u, v)`` that accept NumPy arrays.
    * 3D fields are callables ``field(x, y, z) -> (u, v, w)``.
    * Closed paths are (n, 2) arrays of vertices listed counter-clockwise; the last vertex
      joins back to the first. Circulation is positive counter-clockwise.
    * Parametric surfaces are callables ``surface(s, t) -> (x, y, z)`` with the normal
      taken as r_s x r_t; the helpers here orient that normal outward (sphere) or +z (disc).

Run:  python solution.py          (prints the verification table, writes circulation_vortex.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-11-line-and-surface-integrals
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np

Field2D = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]
Field3D = Callable[[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]
Surface = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]

RHO_SL = 1.225          # kg/m**3, ISA sea-level density
RHO_10KM = 0.4135       # kg/m**3, ISA density at 10 km geopotential


# ----------------------------------------------------------------------------- 2D fields

def vortex_velocity(x: np.ndarray, y: np.ndarray, gamma: float,
                    center: tuple[float, float] = (0.0, 0.0)) -> tuple[np.ndarray, np.ndarray]:
    """Velocity (u, v) of a 2D point vortex of circulation ``gamma`` (m**2/s, positive counter-clockwise).

    v_theta = gamma / (2 pi r); the field is irrotational everywhere except at the centre.
    """
    dx = np.asarray(x, dtype=float) - center[0]
    dy = np.asarray(y, dtype=float) - center[1]
    r2 = dx**2 + dy**2
    return -gamma * dy / (2.0 * np.pi * r2), gamma * dx / (2.0 * np.pi * r2)


def source_velocity(x: np.ndarray, y: np.ndarray, q: float,
                    center: tuple[float, float] = (0.0, 0.0)) -> tuple[np.ndarray, np.ndarray]:
    """Velocity of a 2D point source of volume flow rate ``q`` per unit depth (m**2/s): v_r = q / (2 pi r)."""
    dx = np.asarray(x, dtype=float) - center[0]
    dy = np.asarray(y, dtype=float) - center[1]
    r2 = dx**2 + dy**2
    return q * dx / (2.0 * np.pi * r2), q * dy / (2.0 * np.pi * r2)


def uniform_velocity(x: np.ndarray, y: np.ndarray, v_inf: float, alpha: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
    """Uniform stream of speed ``v_inf`` at angle ``alpha`` (rad) above the x axis, broadcast to the shape of x."""
    x = np.asarray(x, dtype=float)
    shape = np.broadcast(x, np.asarray(y, dtype=float)).shape
    return np.full(shape, v_inf * np.cos(alpha)), np.full(shape, v_inf * np.sin(alpha))


def cylinder_flow_velocity(x: np.ndarray, y: np.ndarray, v_inf: float, radius: float,
                           gamma: float) -> tuple[np.ndarray, np.ndarray]:
    """Potential flow past a circular cylinder of radius ``radius`` with circulation ``gamma``.

    Uniform stream along +x plus a doublet of strength 2 pi v_inf radius**2 plus a point vortex at the
    origin. Valid for r >= radius only; the field inside the cylinder is fictitious.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r2 = x**2 + y**2
    r4 = r2**2
    u_d = -v_inf * radius**2 * (x**2 - y**2) / r4
    v_d = -2.0 * v_inf * radius**2 * x * y / r4
    u_v, v_v = vortex_velocity(x, y, gamma)
    return v_inf + u_d + u_v, v_d + v_v


def rankine_vortex_velocity(x: np.ndarray, y: np.ndarray, gamma: float, core_radius: float) -> tuple[np.ndarray, np.ndarray]:
    """Rankine vortex: solid-body rotation inside ``core_radius``, potential vortex outside.

    v_theta = gamma r / (2 pi a**2) for r < a and gamma / (2 pi r) for r >= a; continuous at r = a.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r = np.hypot(x, y)
    r_safe = np.where(r == 0.0, 1.0, r)
    v_theta = np.where(r < core_radius, gamma * r / (2.0 * np.pi * core_radius**2), gamma / (2.0 * np.pi * r_safe))
    u = np.where(r == 0.0, 0.0, -v_theta * y / r_safe)
    v = np.where(r == 0.0, 0.0, v_theta * x / r_safe)
    return u, v


def gravity_field(x: np.ndarray, y: np.ndarray, mu: float) -> tuple[np.ndarray, np.ndarray]:
    """Gravitational acceleration -mu r / |r|**3 (m/s**2) of a point mass with mu = G M at the origin."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r3 = (x**2 + y**2) ** 1.5
    return -mu * x / r3, -mu * y / r3


def gravity_potential(x: np.ndarray, y: np.ndarray, mu: float) -> np.ndarray:
    """Potential phi = -mu / r (J/kg) with F = -grad(phi); zero at infinity."""
    return -mu / np.hypot(np.asarray(x, dtype=float), np.asarray(y, dtype=float))


# ----------------------------------------------------------------------------- paths and line integrals

def circle_path(radius: float, n: int, center: tuple[float, float] = (0.0, 0.0)) -> np.ndarray:
    """``n`` vertices of a counter-clockwise circle, shape (n, 2), first vertex at angle 0, not repeated."""
    if n < 3 or radius <= 0.0:
        raise ValueError("need n >= 3 vertices and radius > 0")
    th = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([center[0] + radius * np.cos(th), center[1] + radius * np.sin(th)])


def rectangle_path(xmin: float, xmax: float, ymin: float, ymax: float, n_per_side: int) -> np.ndarray:
    """4 * n_per_side vertices of a counter-clockwise rectangle, shape (4n, 2), corners included once."""
    if n_per_side < 1 or xmax <= xmin or ymax <= ymin:
        raise ValueError("need n_per_side >= 1 and xmax > xmin, ymax > ymin")
    s = np.arange(n_per_side) / n_per_side
    bottom = np.column_stack([xmin + (xmax - xmin) * s, np.full(n_per_side, ymin)])
    right = np.column_stack([np.full(n_per_side, xmax), ymin + (ymax - ymin) * s])
    top = np.column_stack([xmax - (xmax - xmin) * s, np.full(n_per_side, ymax)])
    left = np.column_stack([np.full(n_per_side, xmin), ymax - (ymax - ymin) * s])
    return np.vstack([bottom, right, top, left])


def _segments(points: np.ndarray, closed: bool) -> tuple[np.ndarray, np.ndarray]:
    """Return (midpoints, displacement vectors) of the polygon's segments, each of shape (n_segments, 2).

    A closed polygon has as many segments as vertices (the last joins back to the first);
    an open one has one fewer. Raise ``ValueError`` unless ``points`` is (n, 2) with n >= 2.
    """
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2 or pts.shape[0] < 2:
        raise ValueError("points must be an (n, 2) array with n >= 2")
    nxt = np.roll(pts, -1, axis=0) if closed else pts[1:]
    cur = pts if closed else pts[:-1]
    return 0.5 * (cur + nxt), nxt - cur


def line_integral(field: Field2D, points: np.ndarray, closed: bool = True) -> float:
    """Midpoint-rule estimate of the line integral of F . dr along the polygon through ``points``.

    Each segment contributes F(midpoint) . (p_next - p). Second-order accurate for a smooth field.
    """
    mid, d = _segments(points, closed)
    u, v = field(mid[:, 0], mid[:, 1])
    return float(np.sum(u * d[:, 0] + v * d[:, 1]))


def circulation(field: Field2D, points: np.ndarray) -> float:
    """Circulation: the line integral of velocity around the closed path, counter-clockwise positive (m**2/s)."""
    return line_integral(field, points, closed=True)


def flux_2d(field: Field2D, points: np.ndarray) -> float:
    """Outward flux of a 2D field through the closed counter-clockwise path: the integral of F . n ds.

    With dr = (dx, dy) along the path, the outward normal element is n ds = (dy, -dx).
    """
    mid, d = _segments(points, closed=True)
    u, v = field(mid[:, 0], mid[:, 1])
    return float(np.sum(u * d[:, 1] - v * d[:, 0]))


def work_along_path(field: Field2D, points: np.ndarray) -> float:
    """Work per unit mass done by a force field along the open polygon through ``points`` (J/kg)."""
    return line_integral(field, points, closed=False)


def kutta_joukowski_lift(rho: float, v_inf: float, gamma: float) -> float:
    """Lift per unit span L' (N/m) in the +y direction for a stream along +x with circulation ``gamma``.

    With circulation positive counter-clockwise, L'_y = -rho v_inf gamma: a clockwise vortex lifts upward.
    """
    return -rho * v_inf * gamma


# ----------------------------------------------------------------------------- 3D fields and surfaces

def point_source_3d(x: np.ndarray, y: np.ndarray, z: np.ndarray, q: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Velocity of a 3D point source of volume flow rate ``q`` (m**3/s): v = q r / (4 pi |r|**3)."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    r3 = (x**2 + y**2 + z**2) ** 1.5
    c = q / (4.0 * np.pi)
    return c * x / r3, c * y / r3, c * z / r3


def sphere_surface(radius: float, center: tuple[float, float, float] = (0.0, 0.0, 0.0)) -> Surface:
    """Parametric sphere r(s, t) with s = polar angle in [0, pi], t = azimuth in [0, 2 pi]; r_s x r_t points outward."""
    def r(s: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return (center[0] + radius * np.sin(s) * np.cos(t),
                center[1] + radius * np.sin(s) * np.sin(t),
                center[2] + radius * np.cos(s))
    return r


def disc_surface(radius: float, z0: float = 0.0) -> Surface:
    """Parametric disc r(s, t) = (s cos t, s sin t, z0), s in [0, radius], t in [0, 2 pi]; r_s x r_t points along +z."""
    def r(s: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return s * np.cos(t), s * np.sin(t), np.full(np.broadcast(s, t).shape, z0)
    return r


def surface_normal(surface: Surface, s: np.ndarray, t: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """The (unnormalised) normal r_s x r_t at (s, t), with the tangents by central differences, shape (..., 3)."""
    xs1, ys1, zs1 = surface(s + h, t)
    xs0, ys0, zs0 = surface(s - h, t)
    xt1, yt1, zt1 = surface(s, t + h)
    xt0, yt0, zt0 = surface(s, t - h)
    r_s = np.stack([xs1 - xs0, ys1 - ys0, zs1 - zs0], axis=-1) / (2.0 * h)
    r_t = np.stack([xt1 - xt0, yt1 - yt0, zt1 - zt0], axis=-1) / (2.0 * h)
    return np.cross(r_s, r_t)


def surface_flux(field: Field3D, surface: Surface, s_range: tuple[float, float], t_range: tuple[float, float],
                 ns: int, nt: int) -> float:
    """Midpoint-rule estimate of the flux integral of F . (r_s x r_t) ds dt over the parameter rectangle."""
    if ns < 1 or nt < 1:
        raise ValueError("need at least one cell in each parameter direction")
    ds = (s_range[1] - s_range[0]) / ns
    dt = (t_range[1] - t_range[0]) / nt
    s = s_range[0] + ds * (np.arange(ns) + 0.5)
    t = t_range[0] + dt * (np.arange(nt) + 0.5)
    S, T = np.meshgrid(s, t, indexing="ij")
    x, y, z = surface(S, T)
    n = surface_normal(surface, S, T)
    u, v, w = field(x, y, z)
    return float(np.sum(u * n[..., 0] + v * n[..., 1] + w * n[..., 2]) * ds * dt)


def mass_flow_rate(rho: float, v: float, radius: float) -> float:
    """Mass flow rate rho v A (kg/s) through a circular inlet of radius ``radius`` with uniform normal speed v."""
    if rho <= 0.0 or radius <= 0.0:
        raise ValueError("density and radius must be positive")
    return rho * v * np.pi * radius**2


# ----------------------------------------------------------------------------- studies

def circulation_convergence(field: Field2D, gamma_exact: float, ns: list[int], radius: float = 1.0) -> np.ndarray:
    """|circulation - gamma_exact| around a circle and around a square of half-width ``radius`` for each n."""
    out = np.empty((len(ns), 2))
    for i, n in enumerate(ns):
        out[i, 0] = abs(circulation(field, circle_path(radius, n)) - gamma_exact)
        out[i, 1] = abs(circulation(field, rectangle_path(-radius, radius, -radius, radius, n // 4)) - gamma_exact)
    return out


def main(out_dir: str | Path = ".") -> None:
    """Print the verification table and save the vortex and convergence figures."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # Lifting cylinder: v_inf = 50 m/s, radius 0.5 m, clockwise circulation of 20 m**2/s.
    v_inf, a, gamma = 50.0, 0.5, -20.0
    cyl = lambda x, y: cylinder_flow_velocity(x, y, v_inf, a, gamma)  # noqa: E731
    print("Circulation around the lifting cylinder (exact: %.1f m^2/s)" % gamma)
    for name, path in [("circle r = 1 m", circle_path(1.0, 400)), ("circle r = 3 m", circle_path(3.0, 400)),
                       ("square 4 m x 4 m", rectangle_path(-2, 2, -2, 2, 100)),
                       ("circle r = 0.3 m at (2, 2), not enclosing", circle_path(0.3, 400, (2.0, 2.0)))]:
        print(f"  {name:42s} {circulation(cyl, path):+10.4f}")
    lift = kutta_joukowski_lift(RHO_SL, v_inf, gamma)
    print(f"Kutta-Joukowski lift per span  L' = -rho V Gamma = {lift:.1f} N/m")

    # Wake vortex of a transport aircraft on approach, modelled as a Rankine vortex.
    weight, span, v_app, core = 686700.0, 34.0, 70.0, 1.5
    gamma_w = 4.0 * weight / (np.pi * RHO_SL * v_app * span)
    rank = lambda x, y: rankine_vortex_velocity(x, y, gamma_w, core)  # noqa: E731
    print(f"\nWake vortex: Gamma = {gamma_w:.1f} m^2/s, peak swirl {gamma_w / (2 * np.pi * core):.1f} m/s at r = {core} m")
    radii = np.array([0.5, 1.0, 1.5, 3.0, 10.0])
    for r in radii:
        print(f"  circulation around r = {r:5.1f} m: {circulation(rank, circle_path(r, 2000)):9.2f}"
              f"   (solid-body prediction {gamma_w * min(1.0, (r / core) ** 2):9.2f})")

    # Work in a gravitational field: path independence.
    mu_earth = 3.986004418e14
    p1 = np.array([6.771e6, 0.0])
    p2 = np.array([0.0, 7.171e6])
    straight = np.linspace(p1, p2, 2001)
    th = np.linspace(0.0, np.pi / 2, 2001)
    arc = np.column_stack([np.interp(th, [0, np.pi / 2], [p1[0], p2[1]]) * np.cos(th),
                           np.interp(th, [0, np.pi / 2], [p1[0], p2[1]]) * np.sin(th)])
    grav = lambda x, y: gravity_field(x, y, mu_earth)  # noqa: E731
    w1, w2 = work_along_path(grav, straight), work_along_path(grav, arc)
    dphi = gravity_potential(*p1, mu_earth) - gravity_potential(*p2, mu_earth)
    print(f"\nGravity work 400 km -> 800 km altitude: straight {w1:.6e}, arc {w2:.6e}, -delta(phi) {dphi:.6e} J/kg")

    # Mass flow through a turbofan inlet at 10 km, 250 m/s, radius 1 m.
    rho, v, r_in = RHO_10KM, 250.0, 1.0
    mflux = lambda x, y, z: (np.zeros_like(x), np.zeros_like(x), np.full_like(x, rho * v))  # noqa: E731
    mdot = surface_flux(mflux, disc_surface(r_in), (0.0, r_in), (0.0, 2 * np.pi), 50, 50)
    print(f"Inlet mass flow: surface integral {mdot:.3f} kg/s, rho V A = {mass_flow_rate(rho, v, r_in):.3f} kg/s")
    src = lambda x, y, z: point_source_3d(x, y, z, 2.0)  # noqa: E731
    for R in (0.5, 2.0):
        print(f"Source q = 2 m^3/s through sphere R = {R}: {surface_flux(src, sphere_surface(R), (0, np.pi), (0, 2 * np.pi), 60, 60):.6f}")

    # Figures.
    ns = [8, 16, 32, 64, 128, 256, 512]
    errs = circulation_convergence(lambda x, y: vortex_velocity(x, y, 1.0), 1.0, ns)
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    X, Y = np.meshgrid(np.linspace(-2, 2, 24), np.linspace(-2, 2, 24))
    U, V = cyl(X, Y)
    mask = np.hypot(X, Y) < a
    U[mask], V[mask] = np.nan, np.nan
    ax[0].quiver(X, Y, U, V, scale=1500, width=0.003)
    ax[0].add_patch(plt.Circle((0, 0), a, color="0.6"))
    ax[0].set_aspect("equal")
    ax[0].set_title("Lifting cylinder: V = 50 m/s, Gamma = -20 m^2/s")
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("y [m]")
    ax[1].loglog(ns, errs[:, 0], "o-", label="circle")
    ax[1].loglog(ns, errs[:, 1], "s-", label="square")
    ax[1].loglog(ns, errs[0, 0] * (ns[0] / np.array(ns, dtype=float)) ** 2, "k:", label="n^-2")
    ax[1].set_xlabel("segments n")
    ax[1].set_ylabel("|circulation error| / Gamma")
    ax[1].legend()
    ax[1].grid(True, which="both")
    fig.tight_layout()
    fig.savefig(out_dir / "circulation_vortex.png", dpi=150)
    print(f"saved {out_dir / 'circulation_vortex.png'}")


if __name__ == "__main__":
    main()
