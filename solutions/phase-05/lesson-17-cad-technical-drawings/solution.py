"""P05.17 CAD & Technical Drawings I: reference solution.

The CAD work of this lesson is done in FreeCAD (or any free CAD tool). This module is the independent,
tested check on what the CAD model and its drawing claim. It does not read the CAD file; it works from
the bracket's parameters and the dimension scheme on the drawing:

* mass properties of an extruded plate with holes: area, centroid and second moments of the outline by
  Green's theorem (the shoelace formulas), exact circle formulas for the holes, then the prism's volume,
  mass, centre of mass and inertia tensor, checked against a raster integration and against a polygon
  approximation of each hole;
* tolerance stack-up: worst case, root-sum-square, and a seeded Monte Carlo that confirms the RSS
  statistics; baseline against chain dimensioning; the floating-fastener rule for hole positions;
* the dimension scheme as a geometric constraint system: rank of its Jacobian gives the remaining
  degrees of freedom and the redundant dimensions (what the FreeCAD sketcher reports as "fully
  constrained" or "redundant constraints"), and a least-squares solve rebuilds the geometry from the
  dimensions alone.

Units: millimetres and kilograms for the part (the drawing's units), SI elsewhere.

Run:  python main.py [out_dir]     (prints the checks, writes bracket_check.png)
Test: pytest curriculum/phase-05/lesson-17-cad-technical-drawings
"""
from __future__ import annotations

import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares


# ---------------------------------------------------------------- the bracket

@dataclass
class Bracket:
    """A flat mounting bracket: an outline polygon (counter-clockwise vertices, mm) with circular holes,
    extruded to thickness t (mm), of density rho (kg/m^3).

    The default is the lesson's actuator-mounting bracket in 2024-T3 plate (typical density 2780 kg/m^3):
    a 60 x 40 mm plate with the upper-right corner cut by a sloping edge from (60, 15) to (20, 40), two
    6.6 mm clearance holes for M6 bolts and one 12 mm lightening hole, 4 mm thick."""
    outline: np.ndarray = field(default_factory=lambda: np.array(
        [[0.0, 0.0], [60.0, 0.0], [60.0, 15.0], [20.0, 40.0], [0.0, 40.0]]))
    holes: list[tuple[float, float, float]] = field(default_factory=lambda: [
        (10.0, 10.0, 3.3), (10.0, 30.0, 3.3), (35.0, 15.0, 6.0)])     # (x, y, radius) in mm
    t: float = 4.0
    rho: float = 2780.0


def polygon_properties(P: np.ndarray) -> dict[str, float]:
    """Area, first moments and second moments of a simple polygon by Green's theorem.

    With cross_i = x_i y_{i+1} - x_{i+1} y_i:  A = sum cross / 2;  Q_y = int x dA = sum (x_i + x_{i+1}) cross / 6;
    Q_x = int y dA = sum (y_i + y_{i+1}) cross / 6;  I_xx = int y^2 dA = sum (y_i^2 + y_i y_{i+1} + y_{i+1}^2) cross / 12;
    I_yy = int x^2 dA likewise in x;  I_xy = int x y dA = sum (x_i y_{i+1} + 2 x_i y_i + 2 x_{i+1} y_{i+1} + x_{i+1} y_i) cross / 24.
    Moments are about the origin. A clockwise polygon is reversed first, so the area is positive."""
    P = np.asarray(P, dtype=float)
    x, y = P[:, 0], P[:, 1]
    xn, yn = np.roll(x, -1), np.roll(y, -1)
    cr = x * yn - xn * y
    A = cr.sum() / 2.0
    sgn = 1.0 if A >= 0 else -1.0
    cr = cr * sgn
    return {
        "A": abs(A),
        "Qy": float(np.sum((x + xn) * cr) / 6.0),
        "Qx": float(np.sum((y + yn) * cr) / 6.0),
        "Ixx": float(np.sum((y ** 2 + y * yn + yn ** 2) * cr) / 12.0),
        "Iyy": float(np.sum((x ** 2 + x * xn + xn ** 2) * cr) / 12.0),
        "Ixy": float(np.sum((x * yn + 2 * x * y + 2 * xn * yn + xn * y) * cr) / 24.0),
    }


def circle_properties(xc: float, yc: float, r: float) -> dict[str, float]:
    """Same quantities as polygon_properties for a full circle, exactly: A = pi r^2, first moments A x_c,
    A y_c, second moments pi r^4/4 + A y_c^2 (parallel axis), product A x_c y_c."""
    A = np.pi * r ** 2
    return {"A": A, "Qy": A * xc, "Qx": A * yc,
            "Ixx": np.pi * r ** 4 / 4 + A * yc ** 2, "Iyy": np.pi * r ** 4 / 4 + A * xc ** 2, "Ixy": A * xc * yc}


def regular_polygon(xc: float, yc: float, r: float, n: int) -> np.ndarray:
    """n-gon inscribed in a circle, counter-clockwise: what a mesh or an STL export does to a hole."""
    a = 2.0 * np.pi * np.arange(n) / n
    return np.column_stack([xc + r * np.cos(a), yc + r * np.sin(a)])


def section_properties(b: Bracket, hole_polygon_sides: int | None = None) -> dict[str, float]:
    """Area properties of the plate (outline minus holes), about the origin and about the centroid.
    Holes are exact circles, or n-gons when hole_polygon_sides is given (to check convergence)."""
    tot = polygon_properties(b.outline)
    for (xc, yc, r) in b.holes:
        h = circle_properties(xc, yc, r) if hole_polygon_sides is None else polygon_properties(
            regular_polygon(xc, yc, r, hole_polygon_sides))
        for k in tot:
            tot[k] -= h[k]
    A = tot["A"]
    xc, yc = tot["Qy"] / A, tot["Qx"] / A
    tot.update({"xc": xc, "yc": yc,
                "Ixx_c": tot["Ixx"] - A * yc ** 2, "Iyy_c": tot["Iyy"] - A * xc ** 2, "Ixy_c": tot["Ixy"] - A * xc * yc})
    return tot


def mass_properties(b: Bracket) -> dict:
    """Volume (mm^3), mass (kg), centre of mass (mm, z = t/2) and the inertia tensor about the centre of
    mass (kg mm^2) of the extruded plate. For a prism of area properties I_xx_c, I_yy_c, I_xy_c and
    thickness t, with density rho' = rho per mm^3:
        J_xx = rho' t I_xx_c + m t^2/12,  J_yy = rho' t I_yy_c + m t^2/12,  J_zz = rho' t (I_xx_c + I_yy_c),
        J_xy = -rho' t I_xy_c,  J_xz = J_yz = 0 (mid-plane symmetry)."""
    s = section_properties(b)
    rho_mm = b.rho * 1e-9                                   # kg per mm^3
    V = s["A"] * b.t
    m = rho_mm * V
    Jxx = rho_mm * b.t * s["Ixx_c"] + m * b.t ** 2 / 12.0
    Jyy = rho_mm * b.t * s["Iyy_c"] + m * b.t ** 2 / 12.0
    Jzz = rho_mm * b.t * (s["Ixx_c"] + s["Iyy_c"])
    Jxy = -rho_mm * b.t * s["Ixy_c"]
    J = np.array([[Jxx, Jxy, 0.0], [Jxy, Jyy, 0.0], [0.0, 0.0, Jzz]])
    return {"V": V, "m": m, "com": np.array([s["xc"], s["yc"], b.t / 2.0]), "J": J, "area": s["A"]}


def point_in_polygon(px: np.ndarray, py: np.ndarray, P: np.ndarray) -> np.ndarray:
    """Even-odd ray-casting test, vectorised over points."""
    inside = np.zeros(np.shape(px), dtype=bool)
    n = len(P)
    for i in range(n):
        x1, y1 = P[i]
        x2, y2 = P[(i + 1) % n]
        cond = (y1 > py) != (y2 > py)
        with np.errstate(divide="ignore", invalid="ignore"):
            xint = x1 + (py - y1) * (x2 - x1) / (y2 - y1)
        inside ^= cond & (px < xint)
    return inside


def raster_properties(b: Bracket, h: float = 0.1) -> dict[str, float]:
    """Independent check: area, centroid and centroidal second moments by midpoint sampling of the part
    on a square grid of spacing h (mm). Error O(h) from the boundary cells, about 1e-3 relative at
    h = 0.1 mm for this part."""
    xmin, ymin = b.outline.min(axis=0)
    xmax, ymax = b.outline.max(axis=0)
    xs = np.arange(xmin + h / 2, xmax, h)
    ys = np.arange(ymin + h / 2, ymax, h)
    X, Y = np.meshgrid(xs, ys)
    inside = point_in_polygon(X, Y, b.outline)
    for (xc, yc, r) in b.holes:
        inside &= (X - xc) ** 2 + (Y - yc) ** 2 > r ** 2
    dA = h * h
    A = inside.sum() * dA
    xc, yc = X[inside].sum() * dA / A, Y[inside].sum() * dA / A
    return {"A": A, "xc": xc, "yc": yc,
            "Ixx_c": ((Y[inside] - yc) ** 2).sum() * dA + inside.sum() * h ** 4 / 12,
            "Iyy_c": ((X[inside] - xc) ** 2).sum() * dA + inside.sum() * h ** 4 / 12,
            "Ixy_c": ((X[inside] - xc) * (Y[inside] - yc)).sum() * dA}


def principal_axes(Ixx: float, Iyy: float, Ixy: float) -> tuple[float, float, float]:
    """Principal second moments (I1 >= I2) and the angle (rad) of the I1 axis from x, by eigh."""
    w, v = np.linalg.eigh(np.array([[Ixx, -Ixy], [-Ixy, Iyy]]))
    ang = float(np.arctan2(v[1, 1], v[0, 1]))
    return float(w[1]), float(w[0]), ang


# ---------------------------------------------------------------- tolerance stack-up

def worst_case_stack(sens: np.ndarray, tol: np.ndarray) -> float:
    """Worst-case half-range of a linear stack y = sum a_i x_i: sum |a_i| t_i."""
    return float(np.sum(np.abs(sens) * np.abs(tol)))


def rss_stack(sens: np.ndarray, tol: np.ndarray) -> float:
    """Root-sum-square half-range: sqrt(sum (a_i t_i)^2). If each tolerance is a +/-3 sigma band of an
    independent normal variable, this is the +/-3 sigma band of the result."""
    return float(np.sqrt(np.sum((np.asarray(sens) * np.asarray(tol)) ** 2)))


def monte_carlo_stack(nominal: np.ndarray, sens: np.ndarray, tol: np.ndarray, n: int = 200_000,
                      seed: int = 17, dist: str = "normal") -> np.ndarray:
    """Samples of y = sum a_i x_i with each x_i drawn about its nominal: 'normal' with sigma = t_i/3, or
    'uniform' on [-t_i, t_i]. Seeded with numpy.random.default_rng(seed)."""
    rng = np.random.default_rng(seed)
    nominal = np.asarray(nominal, dtype=float)
    tol = np.asarray(tol, dtype=float)
    if dist == "normal":
        dev = rng.normal(0.0, tol / 3.0, size=(n, len(tol)))
    elif dist == "uniform":
        dev = rng.uniform(-tol, tol, size=(n, len(tol)))
    else:
        raise ValueError(dist)
    return (nominal + dev) @ np.asarray(sens, dtype=float)


def floating_fastener_position_tolerance(hole_dia: float, fastener_dia: float) -> float:
    """Floating-fastener rule (two parts, clearance holes in both, a bolt and nut): each part's hole may
    be off its true position by a diametral tolerance T = H - F with H the smallest hole and F the largest
    bolt. Returned as the radial allowance (H - F)/2 per part."""
    return 0.5 * (hole_dia - fastener_dia)


# ---------------------------------------------------------------- the dimension scheme as constraints

# A sketch: points (x, y) in a flat vector, radii appended at the end. Each constraint is a tuple:
#   ("fixx", i, value)  x_i = value          ("fixy", i, value)   y_i = value
#   ("dx", i, j, d)     x_j - x_i = d        ("dy", i, j, d)      y_j - y_i = d
#   ("dist", i, j, d)   |p_j - p_i| = d      ("radius", k, r)     radius k = r
#   ("horizontal", i, j)  y_j = y_i          ("vertical", i, j)   x_j = x_i
#   ("angle", i, j, k, a)  angle at j between the lines j->i and j->k equals a (rad)
#   ("on_segment", c, i, j)  point c is the point (i + j)/2  -- not used by default

def constraint_residuals(z: np.ndarray, cons: list[tuple], n_pts: int) -> np.ndarray:
    """Residual of every constraint at the geometry z = (x0, y0, x1, y1, ..., r0, r1, ...)."""
    P = z[:2 * n_pts].reshape(n_pts, 2)
    R = z[2 * n_pts:]
    out = []
    for c in cons:
        k = c[0]
        if k == "fixx":
            out.append(P[c[1], 0] - c[2])
        elif k == "fixy":
            out.append(P[c[1], 1] - c[2])
        elif k == "dx":
            out.append(P[c[2], 0] - P[c[1], 0] - c[3])
        elif k == "dy":
            out.append(P[c[2], 1] - P[c[1], 1] - c[3])
        elif k == "dist":
            out.append(np.linalg.norm(P[c[2]] - P[c[1]]) - c[3])
        elif k == "radius":
            out.append(R[c[1]] - c[2])
        elif k == "horizontal":
            out.append(P[c[2], 1] - P[c[1], 1])
        elif k == "vertical":
            out.append(P[c[2], 0] - P[c[1], 0])
        elif k == "angle":
            u = P[c[1]] - P[c[2]]
            v = P[c[3]] - P[c[2]]
            out.append(np.arctan2(u[0] * v[1] - u[1] * v[0], u @ v) - c[4])
        else:
            raise ValueError(k)
    return np.array(out, dtype=float)


def constraint_jacobian(z: np.ndarray, cons: list[tuple], n_pts: int, h: float = 1e-6) -> np.ndarray:
    """Jacobian of the residuals by central differences (exact for the linear constraints)."""
    J = np.empty((len(cons), len(z)))
    for k in range(len(z)):
        e = np.zeros_like(z)
        e[k] = h
        J[:, k] = (constraint_residuals(z + e, cons, n_pts) - constraint_residuals(z - e, cons, n_pts)) / (2 * h)
    return J


def constraint_status(z: np.ndarray, cons: list[tuple], n_pts: int) -> dict[str, int]:
    """Degrees of freedom left and redundant constraints, from the rank of the Jacobian at z:
    dof = n_unknowns - rank, redundant = n_constraints - rank. Fully constrained without redundancy
    means dof = 0 and redundant = 0: then the dimensions fix the geometry (locally) uniquely."""
    J = constraint_jacobian(z, cons, n_pts)
    sv = np.linalg.svd(J, compute_uv=False)
    rank = int(np.sum(sv > 1e-8 * max(sv[0], 1.0)))
    return {"unknowns": len(z), "constraints": len(cons), "rank": rank,
            "dof": len(z) - rank, "redundant": len(cons) - rank}


def bracket_sketch(b: Bracket) -> tuple[np.ndarray, list[tuple], int]:
    """The bracket's geometry as a constraint system with the lesson's dimension scheme: datum at the
    lower-left corner (point 0 fixed, like a sketch attached to the origin), every other feature located
    by baseline dimensions from it. Points 0-4 are the outline corners, 5-7 the hole centres.
    Returns (z_nominal, constraints, n_points)."""
    P = np.vstack([b.outline, np.array([[h[0], h[1]] for h in b.holes])])
    R = np.array([h[2] for h in b.holes])
    z = np.concatenate([P.ravel(), R])
    o = np.asarray(b.outline, dtype=float).tolist()          # plain floats: readable constraint tuples
    cons = [("fixx", 0, o[0][0]), ("fixy", 0, o[0][1]),
            ("horizontal", 0, 1), ("dx", 0, 1, o[1][0] - o[0][0]),          # bottom edge, overall width
            ("vertical", 1, 2), ("dy", 1, 2, o[2][1] - o[1][1]),            # right edge height
            ("dx", 0, 3, o[3][0] - o[0][0]), ("dy", 0, 3, o[3][1] - o[0][1]),  # top corner of the sloping edge
            ("vertical", 0, 4), ("horizontal", 4, 3)]                       # left edge, top edge
    for k, (xc, yc, r) in enumerate(b.holes):
        cons += [("dx", 0, 5 + k, float(xc - o[0][0])), ("dy", 0, 5 + k, float(yc - o[0][1])), ("radius", k, float(r))]
    return z, cons, len(P)


def solve_sketch(z0: np.ndarray, cons: list[tuple], n_pts: int) -> np.ndarray:
    """Rebuild the geometry from the dimensions alone, as a CAD sketcher does: least-squares solve of
    the constraint residuals from a starting guess z0."""
    sol = least_squares(constraint_residuals, z0, args=(cons, n_pts), xtol=1e-14, ftol=1e-14, gtol=1e-14)
    return sol.x


def edge_distance_ratio(b: Bracket) -> list[float]:
    """Edge distance over diameter, e/D, for each hole: the distance from the hole's centre to the nearest
    outline edge divided by the hole's diameter. Bolted aluminium plates commonly ask for e/D >= 1.5 so
    that a loaded fastener cannot tear out through the edge; a drawing check reports it."""
    out = []
    P = b.outline
    for (xc, yc, r) in b.holes:
        best = np.inf
        for i in range(len(P)):
            a, c = P[i], P[(i + 1) % len(P)]
            d = c - a
            s = np.clip(((xc - a[0]) * d[0] + (yc - a[1]) * d[1]) / (d @ d), 0.0, 1.0)
            best = min(best, float(np.hypot(xc - a[0] - s * d[0], yc - a[1] - s * d[1])))
        out.append(best / (2 * r))
    return out


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Print the mass-property, tolerance and constraint checks; save bracket_check.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    b = Bracket()

    # ---- mass properties
    s = section_properties(b)
    mp = mass_properties(b)
    ras = raster_properties(b, 0.1)
    print(f"bracket outline area {polygon_properties(b.outline)['A']:.2f} mm^2, holes {polygon_properties(b.outline)['A'] - s['A']:.2f} mm^2, "
          f"net {s['A']:.2f} mm^2")
    print(f"centroid ({s['xc']:.3f}, {s['yc']:.3f}) mm; raster h = 0.1 mm: ({ras['xc']:.3f}, {ras['yc']:.3f}) mm, area {ras['A']:.2f} mm^2")
    print(f"centroidal I_xx {s['Ixx_c']:.0f}, I_yy {s['Iyy_c']:.0f}, I_xy {s['Ixy_c']:.0f} mm^4; raster "
          f"{ras['Ixx_c']:.0f}, {ras['Iyy_c']:.0f}, {ras['Ixy_c']:.0f} mm^4")
    I1, I2, ang = principal_axes(s["Ixx_c"], s["Iyy_c"], s["Ixy_c"])
    print(f"principal I1 {I1:.0f}, I2 {I2:.0f} mm^4, I1 axis at {np.degrees(ang):.1f} deg from x")
    print(f"volume {mp['V']:.1f} mm^3, mass {mp['m']*1e3:.2f} g (2024-T3, {b.rho:.0f} kg/m^3), centre of mass "
          f"({mp['com'][0]:.3f}, {mp['com'][1]:.3f}, {mp['com'][2]:.3f}) mm")
    J = mp["J"]
    print(f"inertia about the centre of mass [kg mm^2]: Jxx {J[0,0]:.4f}, Jyy {J[1,1]:.4f}, Jzz {J[2,2]:.4f}, Jxy {J[0,1]:.4f}")
    print("holes as n-gons (what a mesh export does): n, relative area error")
    for n in (8, 16, 32, 64):
        sn = section_properties(b, n)
        print(f"  {n:3d}  {(sn['A'] - s['A']) / s['A']:+.2e}")
    print("edge distance / hole diameter: " + ", ".join(f"{v:.2f}" for v in edge_distance_ratio(b)))

    # ---- tolerance stack: hole-to-hole spacing, baseline vs chain
    print("\nhole spacing x (hole 1 to lightening hole) from two baseline dimensions, each +/-0.1 mm:")
    sens, tol = np.array([-1.0, 1.0]), np.array([0.1, 0.1])
    print(f"  worst case +/-{worst_case_stack(sens, tol):.3f} mm, RSS +/-{rss_stack(sens, tol):.3f} mm; "
          f"a direct (chain) dimension would give +/-0.100 mm")
    # an assembly gap: housing slot 62.0 +/-0.10 minus bracket width 60.0 +/-0.05 minus two shims 0.50 +/-0.02
    nom = np.array([62.0, 60.0, 0.5, 0.5])
    sens = np.array([1.0, -1.0, -1.0, -1.0])
    tol = np.array([0.10, 0.05, 0.02, 0.02])
    y = monte_carlo_stack(nom, sens, tol, seed=17)
    wc, rs = worst_case_stack(sens, tol), rss_stack(sens, tol)
    print(f"assembly gap nominal {nom @ sens:.3f} mm: worst case +/-{wc:.3f}, RSS +/-{rs:.3f}; Monte Carlo (normal, "
          f"sigma = t/3, 200000 samples) 3 sigma = {3*np.std(y):.4f}, fraction inside RSS band {np.mean(np.abs(y - nom @ sens) <= rs):.4f}")
    yu = monte_carlo_stack(nom, sens, tol, seed=17, dist="uniform")
    print(f"  uniform tolerances instead: 3 sigma = {3*np.std(yu):.4f} mm (theory sqrt(3) x RSS = {np.sqrt(3)*rs:.4f}); "
          f"fraction inside the RSS band {np.mean(np.abs(yu - nom @ sens) <= rs):.4f}")
    print(f"floating fastener: 6.6 mm holes, 6.0 mm bolt -> each hole within {floating_fastener_position_tolerance(6.6, 6.0):.2f} mm "
          f"of true position (diameter {2*floating_fastener_position_tolerance(6.6, 6.0):.2f} mm)")

    # ---- constraint check of the dimension scheme
    z, cons, npt_ = bracket_sketch(b)
    st = constraint_status(z, cons, npt_)
    print(f"\ndimension scheme: {st['unknowns']} unknowns, {st['constraints']} dimensions/constraints, rank {st['rank']}: "
          f"{st['dof']} DOF left, {st['redundant']} redundant")
    st2 = constraint_status(z, cons + [("dy", 2, 3, 25.0)], npt_)
    print(f"  plus the height of the sloping edge (25 mm): {st2['dof']} DOF, {st2['redundant']} redundant")
    st3 = constraint_status(z, [c for c in cons if c != ("dx", 0, 3, 20.0)], npt_)
    print(f"  without the 20 mm location of the slope's top corner: {st3['dof']} DOF, {st3['redundant']} redundant")
    rng = np.random.default_rng(1717)
    zs = solve_sketch(z + rng.normal(0.0, 2.0, z.size), cons, npt_)
    print(f"  solved from a guess perturbed by 2 mm: max |z - nominal| = {np.max(np.abs(zs - z)):.1e} mm")

    # ---- figure: a check drawing
    fig, ax = plt.subplots(figsize=(7, 5))
    P = np.vstack([b.outline, b.outline[:1]])
    ax.plot(P[:, 0], P[:, 1], "k-", lw=1.5)
    for (xc, yc, r) in b.holes:
        ax.add_patch(plt.Circle((xc, yc), r, fill=False, color="k"))
        ax.plot([xc - r - 2, xc + r + 2], [yc, yc], "k-.", lw=0.5)
        ax.plot([xc, xc], [yc - r - 2, yc + r + 2], "k-.", lw=0.5)
    ax.plot(s["xc"], s["yc"], "r+", ms=12, label="centroid")
    for (lab, xa, ya, xb, yb, off) in (("60", 0, 0, 60, 0, -6), ("40", 0, 0, 0, 40, -6), ("20", 0, 40, 20, 40, 5)):
        if ya == yb:
            ax.annotate("", (xa, ya + off), (xb, yb + off), arrowprops=dict(arrowstyle="<->", lw=0.7))
            ax.text((xa + xb) / 2, ya + off + 0.8, lab, ha="center", fontsize=8)
        else:
            ax.annotate("", (xa + off, ya), (xb + off, yb), arrowprops=dict(arrowstyle="<->", lw=0.7))
            ax.text(xa + off - 0.8, (ya + yb) / 2, lab, va="center", ha="right", fontsize=8)
    ax.set_aspect("equal"); ax.set_xlim(-12, 68); ax.set_ylim(-10, 46)
    ax.set_title("bracket: check drawing from the parameters (mm)"); ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "bracket_check.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'bracket_check.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
