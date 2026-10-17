"""P01.12 Divergence, Stokes & Green's Theorems: reference solution.

Verify the divergence theorem numerically on a control volume, Green's theorem on a rectangle
and Stokes' theorem on a hemisphere, with divergence and curl computed both symbolically
(SymPy) and by central differences.

Conventions:
    * 3D fields are callables ``field(x, y, z) -> (u, v, w)`` accepting NumPy arrays.
    * 2D fields are callables ``field(x, y) -> (u, v)``.
    * A box is ``bounds = ((x0, x1), (y0, y1), (z0, z1))``; its outward normal is used on all six faces.
    * Closed 2D paths are (n, 2) arrays of vertices, counter-clockwise, last vertex joining the first.
    * Symbolic fields are lists of three SymPy expressions in the symbols returned by ``coords()``.

Run:  python solution.py          (prints the verification table, writes divergence_theorem.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-12-divergence-stokes-green
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np
import sympy as sp

Field3D = Callable[[np.ndarray, np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]
Field2D = Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]
Scalar3D = Callable[[np.ndarray, np.ndarray, np.ndarray], np.ndarray]
Bounds = tuple[tuple[float, float], tuple[float, float], tuple[float, float]]

X, Y, Z = sp.symbols("x y z", real=True)
G_NEWTON = 6.674e-11   # m**3 kg**-1 s**-2


def coords() -> tuple[sp.Symbol, sp.Symbol, sp.Symbol]:
    """Return the shared Cartesian symbols (x, y, z) used by every symbolic field in this module."""
    return X, Y, Z


# ----------------------------------------------------------------------------- symbolic operators

def gradient_sym(phi: sp.Expr) -> list[sp.Expr]:
    """Return grad(phi) = [d phi/dx, d phi/dy, d phi/dz] as SymPy expressions."""
    return [sp.diff(phi, X), sp.diff(phi, Y), sp.diff(phi, Z)]


def divergence_sym(field: list[sp.Expr]) -> sp.Expr:
    """Return div(F) = dF_x/dx + dF_y/dy + dF_z/dz for a three-component symbolic field."""
    if len(field) != 3:
        raise ValueError("a 3D field has exactly three components")
    return sp.diff(field[0], X) + sp.diff(field[1], Y) + sp.diff(field[2], Z)


def curl_sym(field: list[sp.Expr]) -> list[sp.Expr]:
    """Return curl(F) = (dF_z/dy - dF_y/dz, dF_x/dz - dF_z/dx, dF_y/dx - dF_x/dy)."""
    if len(field) != 3:
        raise ValueError("a 3D field has exactly three components")
    fx, fy, fz = field
    return [sp.diff(fz, Y) - sp.diff(fy, Z), sp.diff(fx, Z) - sp.diff(fz, X), sp.diff(fy, X) - sp.diff(fx, Y)]


def lambdify_field(field: list[sp.Expr]) -> Field3D:
    """Turn a symbolic 3-component field into a NumPy callable ``f(x, y, z) -> (u, v, w)`` that broadcasts."""
    fns = [sp.lambdify((X, Y, Z), comp, modules="numpy") for comp in field]

    def f(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
        return tuple(np.broadcast_to(np.asarray(g(x, y, z), dtype=float), x.shape) for g in fns)  # type: ignore[return-value]

    return f


def lambdify_scalar(expr: sp.Expr) -> Scalar3D:
    """Turn a symbolic scalar into a NumPy callable ``f(x, y, z)`` that broadcasts (constants included)."""
    g = sp.lambdify((X, Y, Z), expr, modules="numpy")

    def f(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> np.ndarray:
        x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
        return np.broadcast_to(np.asarray(g(x, y, z), dtype=float), x.shape)

    return f


# ----------------------------------------------------------------------------- finite-difference operators

def divergence_fd(field: Field3D, x: np.ndarray, y: np.ndarray, z: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Central-difference divergence of a 3D field at the points (x, y, z); error O(h**2). Raise if h <= 0."""
    if h <= 0.0:
        raise ValueError("step h must be positive")
    x, y, z = (np.asarray(a, dtype=float) for a in (x, y, z))
    dudx = (field(x + h, y, z)[0] - field(x - h, y, z)[0]) / (2.0 * h)
    dvdy = (field(x, y + h, z)[1] - field(x, y - h, z)[1]) / (2.0 * h)
    dwdz = (field(x, y, z + h)[2] - field(x, y, z - h)[2]) / (2.0 * h)
    return dudx + dvdy + dwdz


def curl_fd(field: Field3D, x: np.ndarray, y: np.ndarray, z: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Central-difference curl of a 3D field at the points (x, y, z); returns shape (..., 3). Raise if h <= 0."""
    if h <= 0.0:
        raise ValueError("step h must be positive")
    x, y, z = (np.asarray(a, dtype=float) for a in (x, y, z))
    fxp, fxm = field(x + h, y, z), field(x - h, y, z)
    fyp, fym = field(x, y + h, z), field(x, y - h, z)
    fzp, fzm = field(x, y, z + h), field(x, y, z - h)
    d_dx = [(fxp[i] - fxm[i]) / (2.0 * h) for i in range(3)]
    d_dy = [(fyp[i] - fym[i]) / (2.0 * h) for i in range(3)]
    d_dz = [(fzp[i] - fzm[i]) / (2.0 * h) for i in range(3)]
    return np.stack([d_dy[2] - d_dz[1], d_dz[0] - d_dx[2], d_dx[1] - d_dy[0]], axis=-1)


# ----------------------------------------------------------------------------- example fields

def contraction_flow(x: np.ndarray, y: np.ndarray, z: np.ndarray, k: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Axisymmetric converging flow v = (-k x / 2, -k y / 2, k z): the core of a wind-tunnel contraction.

    Incompressible (divergence zero) and irrotational (curl zero); k has units 1/s.
    """
    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
    return -0.5 * k * x, -0.5 * k * y, k * z


def rigid_rotation(x: np.ndarray, y: np.ndarray, z: np.ndarray, omega: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Velocity omega x r of a rigid rotation about the z axis: (-omega y, omega x, 0); curl = (0, 0, 2 omega)."""
    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
    return -omega * y, omega * x, np.zeros_like(z)


def uniform_expansion(x: np.ndarray, y: np.ndarray, z: np.ndarray, k: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Velocity k r of a uniform expansion: divergence 3k everywhere (volume grows as exp(3 k t))."""
    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
    return k * x, k * y, k * z


def rankine_vortex_2d(x: np.ndarray, y: np.ndarray, gamma: float, core_radius: float) -> tuple[np.ndarray, np.ndarray]:
    """Rankine vortex velocity (u, v): solid-body core of radius ``core_radius``, potential vortex outside."""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    r = np.hypot(x, y)
    r_safe = np.where(r == 0.0, 1.0, r)
    v_theta = np.where(r < core_radius, gamma * r / (2.0 * np.pi * core_radius**2), gamma / (2.0 * np.pi * r_safe))
    return np.where(r == 0.0, 0.0, -v_theta * y / r_safe), np.where(r == 0.0, 0.0, v_theta * x / r_safe)


def rankine_vorticity(x: np.ndarray, y: np.ndarray, gamma: float, core_radius: float) -> np.ndarray:
    """Vorticity (z component of curl) of the Rankine vortex: gamma / (pi a**2) inside the core, 0 outside."""
    x, y = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float))
    return np.where(np.hypot(x, y) < core_radius, gamma / (np.pi * core_radius**2), 0.0)


def uniform_sphere_gravity(x: np.ndarray, y: np.ndarray, z: np.ndarray, rho_body: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Gravitational acceleration inside a uniform sphere of density ``rho_body``: g = -(4/3) pi G rho r.

    Valid for |r| below the body's radius; div g = -4 pi G rho (Gauss's law for gravity).
    """
    c = -4.0 * np.pi * G_NEWTON * rho_body / 3.0
    x, y, z = np.broadcast_arrays(np.asarray(x, dtype=float), np.asarray(y, dtype=float), np.asarray(z, dtype=float))
    return c * x, c * y, c * z


# ----------------------------------------------------------------------------- box integrals (divergence theorem)

def _check_bounds(bounds: Bounds, n: int) -> None:
    if n < 1:
        raise ValueError("need at least one cell per edge")
    for lo, hi in bounds:
        if hi <= lo:
            raise ValueError("each bound must satisfy hi > lo")


def box_face_fluxes(field: Field3D, bounds: Bounds, n: int) -> dict[str, float]:
    """Outward flux through each of the six faces of the box, midpoint rule with n x n cells per face.

    Keys: 'x0', 'x1', 'y0', 'y1', 'z0', 'z1' (the face where that coordinate is at its lower/upper bound).
    """
    _check_bounds(bounds, n)
    (x0, x1), (y0, y1), (z0, z1) = bounds
    mids = [lo + (hi - lo) * (np.arange(n) + 0.5) / n for lo, hi in bounds]
    cell = [(hi - lo) / n for lo, hi in bounds]
    out: dict[str, float] = {}
    # x faces: integrate u over (y, z)
    Ym, Zm = np.meshgrid(mids[1], mids[2], indexing="ij")
    da = cell[1] * cell[2]
    out["x0"] = float(-np.sum(field(np.full_like(Ym, x0), Ym, Zm)[0]) * da)
    out["x1"] = float(np.sum(field(np.full_like(Ym, x1), Ym, Zm)[0]) * da)
    Xm, Zm = np.meshgrid(mids[0], mids[2], indexing="ij")
    da = cell[0] * cell[2]
    out["y0"] = float(-np.sum(field(Xm, np.full_like(Xm, y0), Zm)[1]) * da)
    out["y1"] = float(np.sum(field(Xm, np.full_like(Xm, y1), Zm)[1]) * da)
    Xm, Ym = np.meshgrid(mids[0], mids[1], indexing="ij")
    da = cell[0] * cell[1]
    out["z0"] = float(-np.sum(field(Xm, Ym, np.full_like(Xm, z0))[2]) * da)
    out["z1"] = float(np.sum(field(Xm, Ym, np.full_like(Xm, z1))[2]) * da)
    return out


def box_flux(field: Field3D, bounds: Bounds, n: int) -> float:
    """Net outward flux of ``field`` through the closed surface of the box (sum of the six face fluxes)."""
    return float(sum(box_face_fluxes(field, bounds, n).values()))


def box_volume_integral(f: Scalar3D, bounds: Bounds, n: int) -> float:
    """Midpoint-rule triple integral of the scalar ``f(x, y, z)`` over the box with n**3 cells."""
    _check_bounds(bounds, n)
    mids = [lo + (hi - lo) * (np.arange(n) + 0.5) / n for lo, hi in bounds]
    dv = np.prod([(hi - lo) / n for lo, hi in bounds])
    Xm, Ym, Zm = np.meshgrid(*mids, indexing="ij")
    return float(np.sum(f(Xm, Ym, Zm)) * dv)


def divergence_theorem_check(field: Field3D, div: Scalar3D, bounds: Bounds, n: int) -> tuple[float, float, float]:
    """Return (surface flux, volume integral of div, |difference|) for the box: the two sides of Gauss's theorem."""
    flux = box_flux(field, bounds, n)
    vol = box_volume_integral(div, bounds, n)
    return flux, vol, abs(flux - vol)


# ----------------------------------------------------------------------------- Green's theorem

def polygon_circulation(field: Field2D, points: np.ndarray) -> float:
    """Midpoint-rule circulation of a 2D field around the closed polygon ``points`` (counter-clockwise)."""
    pts = np.asarray(points, dtype=float)
    if pts.ndim != 2 or pts.shape[1] != 2 or pts.shape[0] < 3:
        raise ValueError("a closed polygon needs an (n, 2) array with n >= 3")
    nxt = np.roll(pts, -1, axis=0)
    mid, d = 0.5 * (pts + nxt), nxt - pts
    u, v = field(mid[:, 0], mid[:, 1])
    return float(np.sum(u * d[:, 0] + v * d[:, 1]))


def rectangle_boundary(xmin: float, xmax: float, ymin: float, ymax: float, n_per_side: int) -> np.ndarray:
    """4 * n_per_side vertices of the counter-clockwise rectangle boundary, corners included once."""
    if n_per_side < 1 or xmax <= xmin or ymax <= ymin:
        raise ValueError("need n_per_side >= 1 and xmax > xmin, ymax > ymin")
    s = np.arange(n_per_side) / n_per_side
    bottom = np.column_stack([xmin + (xmax - xmin) * s, np.full(n_per_side, ymin)])
    right = np.column_stack([np.full(n_per_side, xmax), ymin + (ymax - ymin) * s])
    top = np.column_stack([xmax - (xmax - xmin) * s, np.full(n_per_side, ymax)])
    left = np.column_stack([np.full(n_per_side, xmin), ymax - (ymax - ymin) * s])
    return np.vstack([bottom, right, top, left])


def rectangle_area_integral(f: Callable[[np.ndarray, np.ndarray], np.ndarray], xmin: float, xmax: float,
                            ymin: float, ymax: float, n: int) -> float:
    """Midpoint-rule double integral of ``f(x, y)`` over the rectangle with n x n cells."""
    if n < 1 or xmax <= xmin or ymax <= ymin:
        raise ValueError("need n >= 1 and xmax > xmin, ymax > ymin")
    xm = xmin + (xmax - xmin) * (np.arange(n) + 0.5) / n
    ym = ymin + (ymax - ymin) * (np.arange(n) + 0.5) / n
    Xm, Ym = np.meshgrid(xm, ym, indexing="ij")
    return float(np.sum(f(Xm, Ym)) * (xmax - xmin) * (ymax - ymin) / n**2)


def greens_theorem_check(field: Field2D, curl_z: Callable[[np.ndarray, np.ndarray], np.ndarray],
                         xmin: float, xmax: float, ymin: float, ymax: float, n: int) -> tuple[float, float, float]:
    """Return (circulation around the rectangle, area integral of curl_z, |difference|): Green's theorem."""
    circ = polygon_circulation(field, rectangle_boundary(xmin, xmax, ymin, ymax, n))
    area = rectangle_area_integral(curl_z, xmin, xmax, ymin, ymax, n)
    return circ, area, abs(circ - area)


# ----------------------------------------------------------------------------- Stokes' theorem on a hemisphere

def hemisphere_surface(radius: float) -> Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Upper hemisphere r(s, t) = R (sin s cos t, sin s sin t, cos s), s in [0, pi/2], t in [0, 2 pi]; normal outward (+z at the pole)."""
    def r(s: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return radius * np.sin(s) * np.cos(t), radius * np.sin(s) * np.sin(t), radius * np.cos(s)
    return r


def disc_surface(radius: float, z0: float = 0.0) -> Callable[[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    """Flat disc r(s, t) = (s cos t, s sin t, z0), s in [0, R], t in [0, 2 pi]; normal along +z."""
    def r(s: np.ndarray, t: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        s, t = np.broadcast_arrays(s, t)
        return s * np.cos(t), s * np.sin(t), np.full(s.shape, z0)
    return r


def surface_flux(field: Field3D, surface: Callable, s_range: tuple[float, float], t_range: tuple[float, float],
                 ns: int, nt: int, h: float = 1e-6) -> float:
    """Midpoint-rule flux of ``field`` through the parametric surface, normal r_s x r_t by central differences."""
    if ns < 1 or nt < 1:
        raise ValueError("need at least one cell in each parameter direction")
    ds = (s_range[1] - s_range[0]) / ns
    dt = (t_range[1] - t_range[0]) / nt
    s = s_range[0] + ds * (np.arange(ns) + 0.5)
    t = t_range[0] + dt * (np.arange(nt) + 0.5)
    S, T = np.meshgrid(s, t, indexing="ij")
    p = [np.stack(surface(S + a, T + b), axis=-1) for a, b in ((h, 0), (-h, 0), (0, h), (0, -h))]
    normal = np.cross((p[0] - p[1]) / (2 * h), (p[2] - p[3]) / (2 * h))
    x, y, z = surface(S, T)
    u, v, w = field(x, y, z)
    return float(np.sum(u * normal[..., 0] + v * normal[..., 1] + w * normal[..., 2]) * ds * dt)


def rim_circulation(field: Field3D, radius: float, n: int, z0: float = 0.0) -> float:
    """Circulation of a 3D field around the counter-clockwise circle of radius R in the plane z = z0 (midpoint rule)."""
    if n < 3:
        raise ValueError("need at least three segments")
    th = 2.0 * np.pi * np.arange(n) / n
    pts = np.column_stack([radius * np.cos(th), radius * np.sin(th), np.full(n, z0)])
    nxt = np.roll(pts, -1, axis=0)
    mid, d = 0.5 * (pts + nxt), nxt - pts
    u, v, w = field(mid[:, 0], mid[:, 1], mid[:, 2])
    return float(np.sum(u * d[:, 0] + v * d[:, 1] + w * d[:, 2]))


def stokes_theorem_check(field: Field3D, curl: Field3D, radius: float, n: int) -> tuple[float, float, float]:
    """Return (rim circulation, flux of curl through the hemisphere, flux of curl through the flat disc) for radius R.

    Stokes' theorem says all three agree: the rim is the boundary of both surfaces.
    """
    circ = rim_circulation(field, radius, 4 * n)
    cap = surface_flux(curl, hemisphere_surface(radius), (0.0, np.pi / 2), (0.0, 2 * np.pi), n, n)
    disc = surface_flux(curl, disc_surface(radius), (0.0, radius), (0.0, 2 * np.pi), n, n)
    return circ, cap, disc


# ----------------------------------------------------------------------------- studies

def divergence_theorem_convergence(field: Field3D, div: Scalar3D, bounds: Bounds, ns: list[int]) -> np.ndarray:
    """|flux - volume integral| of the divergence theorem for each resolution in ``ns``."""
    return np.array([divergence_theorem_check(field, div, bounds, n)[2] for n in ns])


def main(out_dir: str | Path = ".") -> None:
    """Print the verification table and save divergence_theorem.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    x, y, z = coords()

    # 1. Wind-tunnel contraction core: incompressible, so what enters the control volume leaves it.
    k = 10.0  # 1/s
    box: Bounds = ((-0.2, 0.2), (-0.2, 0.2), (1.0, 2.0))
    flow = lambda X_, Y_, Z_: contraction_flow(X_, Y_, Z_, k)  # noqa: E731
    faces = box_face_fluxes(flow, box, 40)
    print("Contraction flow v = (-k x/2, -k y/2, k z), k = 10 1/s, box 0.4 x 0.4 x 1 m:")
    for name, val in faces.items():
        print(f"  face {name}: {val:+10.5f} m^3/s")
    print(f"  net outward flux: {sum(faces.values()):+.2e} m^3/s   (divergence theorem with div v = 0)")
    sym = [-k * x / 2, -k * y / 2, k * z]
    print(f"  symbolic div = {divergence_sym(sym)}, curl = {curl_sym(sym)}")

    # 2. Uniform expansion: div = 3k, flux = 3k * volume.
    ke = 2.0
    expand = lambda X_, Y_, Z_: uniform_expansion(X_, Y_, Z_, ke)  # noqa: E731
    flux, vol, diff = divergence_theorem_check(expand, lambda X_, Y_, Z_: np.full(np.shape(X_), 3 * ke), box, 10)
    print(f"\nUniform expansion k = {ke} 1/s: flux {flux:.6f}, volume integral {vol:.6f}, |diff| {diff:.1e} (exact 3 k V = {3 * ke * 0.16:.6f})")

    # 3. A field whose divergence is not constant: F = (y e^x, z sin y, x cos z); convergence of the residual.
    F = [y * sp.exp(x), z * sp.sin(y), x * sp.cos(z)]
    f_num, div_num = lambdify_field(F), lambdify_scalar(divergence_sym(F))
    ns = [2, 4, 8, 16, 32, 64]
    res = divergence_theorem_convergence(f_num, div_num, box, ns)
    print(f"\nF = {F}, div F = {divergence_sym(F)}")
    print("  n      |flux - volume|")
    for n, r in zip(ns, res):
        print(f"  {n:3d}    {r:.3e}")
    print("  observed orders:", np.round(np.log2(res[:-1] / res[1:]), 2))
    Q = [x**2 * y, y * z**2, sp.sin(x) * z]
    print(f"  and for F = {Q}, quadratic in each component's own coordinate, the residual at n = 4 is "
          f"{divergence_theorem_check(lambdify_field(Q), lambdify_scalar(divergence_sym(Q)), box, 4)[2]:.1e}: the discrete theorem is exact")

    # 4. Green's theorem: a Rankine wake vortex inside and outside its core.
    gamma, a = 300.0, 1.5
    rv = lambda X_, Y_: rankine_vortex_2d(X_, Y_, gamma, a)  # noqa: E731
    wz = lambda X_, Y_: rankine_vorticity(X_, Y_, gamma, a)  # noqa: E731
    for name, rect in (("1 m square inside the core", (-0.5, 0.5, -0.5, 0.5)), ("6 m square enclosing the core", (-3, 3, -3, 3)),
                       ("1 m square outside the core", (2.0, 3.0, 2.0, 3.0))):
        circ, area, d = greens_theorem_check(rv, wz, *rect, 200)
        print(f"\nGreen, {name}: circulation {circ:9.4f}, integral of vorticity {area:9.4f}, |diff| {d:.2e} m^2/s")

    # 5. Stokes' theorem for rigid rotation and for a non-trivial field.
    omega, R = 0.5, 2.0
    rot = lambda X_, Y_, Z_: rigid_rotation(X_, Y_, Z_, omega)  # noqa: E731
    curl_rot = lambda X_, Y_, Z_: (np.zeros_like(X_), np.zeros_like(X_), np.full_like(X_, 2 * omega))  # noqa: E731
    circ, cap, disc = stokes_theorem_check(rot, curl_rot, R, 40)
    print(f"\nStokes, rigid rotation omega = {omega} rad/s, R = {R} m: rim {circ:.6f}, hemisphere {cap:.6f}, disc {disc:.6f}, exact 2 pi omega R^2 = {2 * np.pi * omega * R**2:.6f}")
    G = [z, x, y]
    circ, cap, disc = stokes_theorem_check(lambdify_field(G), lambdify_field(curl_sym(G)), R, 40)
    print(f"Stokes, F = {G}, curl = {curl_sym(G)}: rim {circ:.6f}, hemisphere {cap:.6f}, disc {disc:.6f}, exact pi R^2 = {np.pi * R**2:.6f}")

    # Figure: contraction streamlines-ish quiver and the convergence plot.
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    Xg, Zg = np.meshgrid(np.linspace(-0.3, 0.3, 13), np.linspace(0.8, 2.2, 15))
    U, _, W = contraction_flow(Xg, 0.0, Zg, k)
    ax[0].quiver(Xg, Zg, U, W, scale=300, width=0.003)
    ax[0].add_patch(plt.Rectangle((box[0][0], box[2][0]), 0.4, 1.0, fill=False, color="C3", lw=1.5))
    ax[0].set_xlabel("x [m]")
    ax[0].set_ylabel("z [m]")
    ax[0].set_title("Contraction core flow and control volume")
    ax[1].loglog(ns, res, "o-", label="|flux - volume integral|")
    ax[1].loglog(ns, res[0] * (ns[0] / np.array(ns, dtype=float)) ** 2, "k:", label="n^-2")
    ax[1].set_xlabel("cells per edge n")
    ax[1].set_ylabel("residual")
    ax[1].set_title("Divergence theorem residual")
    ax[1].grid(True, which="both")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "divergence_theorem.png", dpi=150)
    print(f"\nsaved {out_dir / 'divergence_theorem.png'}")


if __name__ == "__main__":
    main()
