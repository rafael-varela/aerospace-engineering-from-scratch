"""P05.05 Centroids & Moments of Inertia: reference solution.

How stiff a beam is in bending, and how much it weighs, depend on its cross-section through a handful of
area integrals: the area A, the centroid (where the neutral axis of a beam in pure bending passes), the
second moments Ix, Iy and the product Ixy. This module computes them two independent ways and checks
one against the other:

* composite sections: rectangles, circles and point areas (stringer "booms") added or subtracted, each
  with its own centroidal properties moved to the section centroid by the parallel-axis theorem
  I = I_c + A d^2 (and Ixy = Ixy_c + A dx dy);
* polygon sections: exact integrals over any simple polygon by Green's theorem, one sum over the edges,
  with holes as polygons of negative weight;
* rotation of axes, principal second moments and the principal angle (Mohr's circle for areas), checked
  against numpy.linalg.eigh of the 2 x 2 second-moment matrix;
* thin-walled approximations (terms in t^3 dropped) and their error, and the boom idealisation of a
  stiffened wing box;
* the centroid of a set of masses: an aircraft's centre of gravity in % of the mean aerodynamic chord.

Axes: x horizontal, y vertical, in the plane of the section. Ix = integral y^2 dA (bending about x),
Iy = integral x^2 dA, Ixy = integral x y dA. Units: SI (m, m^2, m^4); main() prints mm and mm^4.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes section_properties.png)
Test: pytest curriculum/phase-05/lesson-05-centroids-moments-of-inertia
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Typical density of 7075-T6 and 2024-T3 aluminium (see the lesson for the source and caveats).
RHO_7075 = 2810.0                 # kg/m^3
RHO_2024 = 2780.0                 # kg/m^3


# ---------------------------------------------------------------- polygon integrals (Green's theorem)

def polygon_integrals(xy: np.ndarray) -> np.ndarray:
    """Area integrals of a simple polygon about the ORIGIN: (A, Qx, Qy, Ixx, Iyy, Ixy) with
    Qx = int y dA, Qy = int x dA, Ixx = int y^2 dA, Iyy = int x^2 dA, Ixy = int x y dA.

    Green's theorem turns each area integral into a sum over the edges (x_i, y_i) -> (x_{i+1}, y_{i+1})
    with c_i = x_i y_{i+1} - x_{i+1} y_i:
        A   = 1/2  sum c_i
        Qy  = 1/6  sum (x_i + x_{i+1}) c_i            Qx  = 1/6  sum (y_i + y_{i+1}) c_i
        Iyy = 1/12 sum (x_i^2 + x_i x_{i+1} + x_{i+1}^2) c_i
        Ixx = 1/12 sum (y_i^2 + y_i y_{i+1} + y_{i+1}^2) c_i
        Ixy = 1/24 sum (x_i y_{i+1} + 2 x_i y_i + 2 x_{i+1} y_{i+1} + x_{i+1} y_i) c_i
    The sums are signed (negative for a clockwise polygon); the result is normalised to positive area
    so the vertex order does not matter. Exact for straight edges."""
    p = np.asarray(xy, dtype=float)
    x, y = p[:, 0], p[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    c = x * y1 - x1 * y
    out = np.array([c.sum() / 2,
                    ((y + y1) * c).sum() / 6,
                    ((x + x1) * c).sum() / 6,
                    ((y * y + y * y1 + y1 * y1) * c).sum() / 12,
                    ((x * x + x * x1 + x1 * x1) * c).sum() / 12,
                    ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * c).sum() / 24])
    return out if out[0] >= 0 else -out


@dataclass
class Section:
    """Section properties: area, centroid and the centroidal second moments (m, m^2, m^4)."""
    A: float
    cx: float
    cy: float
    Ix: float
    Iy: float
    Ixy: float

    def as_tuple(self) -> tuple:
        return (self.A, self.cx, self.cy, self.Ix, self.Iy, self.Ixy)


def section_from_integrals(I: np.ndarray) -> Section:
    """Centroidal properties from integrals about the origin, by the parallel-axis theorem used backwards:
    Ix = Ixx - A cy^2, Iy = Iyy - A cx^2, Ixy = Ixy_O - A cx cy."""
    A, Qx, Qy, Ixx, Iyy, Ixy = I
    if A <= 0:
        raise ValueError("section has no area")
    cx, cy = Qy / A, Qx / A
    return Section(A, cx, cy, Ixx - A * cy ** 2, Iyy - A * cx ** 2, Ixy - A * cx * cy)


def polygon_section(outer: np.ndarray, holes: list[np.ndarray] | None = None) -> Section:
    """Properties of a polygon with optional polygonal holes (holes must lie inside ``outer``)."""
    I = polygon_integrals(outer)
    for h in holes or []:
        I = I - polygon_integrals(h)
    return section_from_integrals(I)


# ---------------------------------------------------------------- composite sections (parallel axis)

@dataclass
class Part:
    """A piece of a composite section with its own centroidal properties. sign = -1 for a hole."""
    A: float
    cx: float
    cy: float
    Ix: float = 0.0
    Iy: float = 0.0
    Ixy: float = 0.0
    sign: float = 1.0


def rect(b: float, h: float, cx: float, cy: float, sign: float = 1.0) -> Part:
    """Rectangle of width b (along x) and height h centred at (cx, cy): Ix = b h^3 / 12, Iy = h b^3 / 12."""
    return Part(b * h, cx, cy, b * h ** 3 / 12.0, h * b ** 3 / 12.0, 0.0, sign)


def circle(r: float, cx: float, cy: float, sign: float = 1.0) -> Part:
    """Solid circle of radius r: Ix = Iy = pi r^4 / 4."""
    return Part(np.pi * r ** 2, cx, cy, np.pi * r ** 4 / 4.0, np.pi * r ** 4 / 4.0, 0.0, sign)


def boom(area: float, cx: float, cy: float) -> Part:
    """Point area (a stringer or spar cap idealised as a boom): its own second moments are neglected."""
    return Part(area, cx, cy)


def composite_section(parts: list[Part]) -> Section:
    """Section properties of a composite by the parallel-axis theorem.

    Centroid: c = sum s_i A_i c_i / sum s_i A_i. Second moments about it:
    Ix = sum s_i (Ix_i + A_i (cy_i - cy)^2), Iy likewise, Ixy = sum s_i (Ixy_i + A_i dx_i dy_i)."""
    s = np.array([p.sign for p in parts])
    A = np.array([p.A for p in parts]) * s
    At = A.sum()
    if At <= 0:
        raise ValueError("composite has no net area")
    cx = float(np.sum(A * [p.cx for p in parts]) / At)
    cy = float(np.sum(A * [p.cy for p in parts]) / At)
    Ix = sum(p.sign * (p.Ix + p.A * (p.cy - cy) ** 2) for p in parts)
    Iy = sum(p.sign * (p.Iy + p.A * (p.cx - cx) ** 2) for p in parts)
    Ixy = sum(p.sign * (p.Ixy + p.A * (p.cx - cx) * (p.cy - cy)) for p in parts)
    return Section(float(At), cx, cy, float(Ix), float(Iy), float(Ixy))


# ---------------------------------------------------------------- rotation and principal axes

def rotate_axes(Ix: float, Iy: float, Ixy: float, theta: float) -> tuple[float, float, float]:
    """Second moments about axes u, v turned by theta (anticlockwise) from x, y through the same point:
        Iu  = Ix cos^2 - 2 Ixy sin cos + Iy sin^2
        Iv  = Ix sin^2 + 2 Ixy sin cos + Iy cos^2
        Iuv = (Ix - Iy) sin cos + Ixy (cos^2 - sin^2)
    Iu + Iv = Ix + Iy and Iu Iv - Iuv^2 = Ix Iy - Ixy^2 for every theta (the invariants)."""
    c, s = np.cos(theta), np.sin(theta)
    Iu = Ix * c * c - 2 * Ixy * s * c + Iy * s * s
    Iv = Ix * s * s + 2 * Ixy * s * c + Iy * c * c
    Iuv = (Ix - Iy) * s * c + Ixy * (c * c - s * s)
    return float(Iu), float(Iv), float(Iuv)


def principal_axes(Ix: float, Iy: float, Ixy: float) -> tuple[float, float, float]:
    """(I_max, I_min, theta_p): principal second moments from Mohr's circle, centre (Ix + Iy)/2 and radius
    sqrt(((Ix - Iy)/2)^2 + Ixy^2), and the angle theta_p (rad, in (-pi/2, pi/2]) of the axis about which
    I = I_max, from tan 2 theta_p = -2 Ixy / (Ix - Iy). Iuv = 0 on the principal axes."""
    c = 0.5 * (Ix + Iy)
    R = float(np.hypot(0.5 * (Ix - Iy), Ixy))
    th = 0.5 * np.arctan2(-2.0 * Ixy, Ix - Iy)
    return c + R, c - R, float(th)


def principal_axes_eig(Ix: float, Iy: float, Ixy: float) -> tuple[np.ndarray, np.ndarray]:
    """Independent route: eigenvalues (ascending) and eigenvectors of C = [[Iy, Ixy], [Ixy, Ix]].
    For an axis along the unit vector e, I_axis = n . C n with n = (-e_y, e_x) the normal to the axis, so
    each eigenvector of C is the NORMAL of a principal axis."""
    return np.linalg.eigh(np.array([[Iy, Ixy], [Ixy, Ix]]))


def radius_of_gyration(I: float, A: float) -> float:
    """r = sqrt(I / A): the distance at which the whole area would give the same second moment."""
    return float(np.sqrt(I / A))


def section_moduli(sec: Section, outline: np.ndarray) -> tuple[float, float]:
    """Elastic section moduli about the centroidal x axis to the top and bottom fibres: Ix / c, with c the
    distance from the centroid to the extreme y of ``outline``. Used in P05.10 for sigma = M / S."""
    y = np.asarray(outline, float)[:, 1]
    return sec.Ix / (y.max() - sec.cy), sec.Ix / (sec.cy - y.min())


# ---------------------------------------------------------------- shapes

def rect_polygon(b: float, h: float, cx: float = 0.0, cy: float = 0.0) -> np.ndarray:
    return np.array([[cx - b / 2, cy - h / 2], [cx + b / 2, cy - h / 2], [cx + b / 2, cy + h / 2], [cx - b / 2, cy + h / 2]])


def circle_polygon(r: float, n: int, cx: float = 0.0, cy: float = 0.0) -> np.ndarray:
    """Regular n-gon inscribed in a circle of radius r."""
    t = 2 * np.pi * np.arange(n) / n
    return np.column_stack([cx + r * np.cos(t), cy + r * np.sin(t)])


def i_section_parts(h: float, b: float, tf: float, tw: float) -> list[Part]:
    """Symmetric I-section (spar): overall height h, flange width b and thickness tf, web thickness tw,
    centred at the origin. Two flange rectangles and the web between them."""
    return [rect(b, tf, 0.0, (h - tf) / 2), rect(b, tf, 0.0, -(h - tf) / 2), rect(tw, h - 2 * tf, 0.0, 0.0)]


def i_section_polygon(h: float, b: float, tf: float, tw: float) -> np.ndarray:
    """The same I-section as one 12-vertex polygon (anticlockwise)."""
    hb, hw, hh = b / 2, tw / 2, h / 2
    return np.array([[-hb, -hh], [hb, -hh], [hb, -hh + tf], [hw, -hh + tf], [hw, hh - tf], [hb, hh - tf],
                     [hb, hh], [-hb, hh], [-hb, hh - tf], [-hw, hh - tf], [-hw, -hh + tf], [-hb, -hh + tf]])


def z_section_polygon(h: float, b: float, t: float) -> np.ndarray:
    """Z-section stringer: web of height h (outside to outside) and thickness t along y, top flange of
    width b (web included) running to +x, bottom flange running to -x, flange thickness t. The web's
    centre line is x = 0. Eight vertices, anticlockwise."""
    ht, hh = t / 2, h / 2
    return np.array([[-(b - ht), -hh], [ht, -hh], [ht, hh - t], [b - ht, hh - t], [b - ht, hh], [-ht, hh],
                     [-ht, -hh + t], [-(b - ht), -hh + t]])


def z_section_parts(h: float, b: float, t: float) -> list[Part]:
    """The Z-section as three non-overlapping rectangles: web (full height) and two flanges (b - t wide)."""
    ht = t / 2
    fw = b - t
    return [rect(t, h, 0.0, 0.0), rect(fw, t, ht + fw / 2, h / 2 - ht), rect(fw, t, -(ht + fw / 2), -(h / 2 - ht))]


def z_section_thin(h: float, b: float, t: float) -> Section:
    """Thin-walled approximation of the Z: the wall is a line along its mid-line (web of length h - t, flanges
    of length b - t/2 measured from the web's centre line), and terms in t^3 are dropped:
    Ix = t d^3 / 12 + 2 f t (d/2)^2, Iy = 2 t f^3 / 3, Ixy = 2 f t (f/2)(d/2) with d = h - t, f = b - t/2
    (both flanges lie in the quadrants where x y > 0, so Ixy > 0)."""
    d, f = h - t, b - t / 2
    A = (d + 2 * f) * t
    return Section(A, 0.0, 0.0, t * d ** 3 / 12 + 2 * f * t * (d / 2) ** 2, 2 * t * f ** 3 / 3, 2 * f * t * (f / 2) * (d / 2))


def angle_polygon(a: float, b: float, t: float) -> np.ndarray:
    """Unequal angle: leg of length a along y and leg of length b along x, thickness t, heel at the origin."""
    return np.array([[0.0, 0.0], [b, 0.0], [b, t], [t, t], [t, a], [0.0, a]])


def wing_box_parts(width: float, height: float, t_skin: float, t_web: float, cap_area: float,
                   n_stringers: int, stringer_area: float, idealisation: str = "full") -> list[Part]:
    """Single-cell wing box (symmetric about x and y): two spar webs of thickness t_web, upper and lower skins
    of thickness t_skin between the webs, four spar caps (booms of cap_area at the corners) and
    n_stringers stringers (booms) evenly spaced along each skin.

    idealisation: 'full' (skins and webs as thin rectangles plus the booms), 'booms' (caps and stringers
    only: skins and webs dropped), or 'lumped' (booms only, but each skin panel of width w between two
    booms adds t w / 6 (2 + s2/s1) to each boom, with s2/s1 = 1 for two booms at the same height and -1
    for the two ends of a spar web under bending about x: the classic boom idealisation)."""
    xs = np.linspace(-width / 2, width / 2, n_stringers + 2)
    areas = np.full(len(xs), stringer_area)
    areas[0] = areas[-1] = cap_area
    if idealisation == "lumped":
        w = xs[1] - xs[0]
        add = np.full(len(xs), 2 * t_skin * w / 6 * (2 + 1))      # two panels meet at an interior boom
        add[0] = add[-1] = t_skin * w / 6 * (2 + 1) + t_web * height / 6 * (2 - 1)
        areas = areas + add
    elif idealisation not in ("full", "booms"):
        raise ValueError(f"unknown idealisation {idealisation!r}")
    parts = [boom(a, x, sy * height / 2) for x, a in zip(xs, areas) for sy in (-1, 1)]
    if idealisation == "full":
        parts += [rect(width - t_web, t_skin, 0.0, sy * (height - t_skin) / 2) for sy in (-1, 1)]
        parts += [rect(t_web, height, sx * (width - t_web) / 2, 0.0) for sx in (-1, 1)]
    return parts


# ---------------------------------------------------------------- centroid of masses: weight and balance

def composite_centroid(masses: np.ndarray, positions: np.ndarray) -> np.ndarray:
    """Centre of mass sum m_i r_i / sum m_i of point masses (positions of shape (n,) or (n, d))."""
    m = np.asarray(masses, dtype=float)
    r = np.asarray(positions, dtype=float)
    return np.tensordot(m, r, axes=(0, 0)) / m.sum()


def cg_percent_mac(x_cg: float, x_lemac: float, mac: float) -> float:
    """CG position in percent of the mean aerodynamic chord, measured aft of its leading edge."""
    return 100.0 * (x_cg - x_lemac) / mac


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save section_properties.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    mm, mm2, mm4 = 1e3, 1e6, 1e12

    def show(label: str, s: Section) -> None:
        print(f"  {label:26s} A = {s.A * mm2:9.2f} mm^2, c = ({s.cx * mm:8.4f}, {s.cy * mm:8.4f}) mm, "
              f"Ix = {s.Ix * mm4:12.1f}, Iy = {s.Iy * mm4:11.1f}, Ixy = {s.Ixy * mm4:10.1f} mm^4")

    # ---- I spar
    h, b, tf, tw = 0.150, 0.050, 0.006, 0.003
    sc = composite_section(i_section_parts(h, b, tf, tw))
    sp = polygon_section(i_section_polygon(h, b, tf, tw))
    print(f"I spar {h * mm:.0f} x {b * mm:.0f} mm, flanges {tf * mm:.0f} mm, web {tw * mm:.0f} mm:")
    show("composite (parallel axis)", sc)
    show("polygon (Green)", sp)
    flange_share = 2 * (b * tf * ((h - tf) / 2) ** 2 + b * tf ** 3 / 12) / sc.Ix
    St, Sb = section_moduli(sc, i_section_polygon(h, b, tf, tw))
    print(f"  flanges carry {100 * flange_share:.1f} % of Ix with {100 * 2 * b * tf / sc.A:.1f} % of the area; "
          f"S = {St * 1e9:.0f} mm^3, r_x = {radius_of_gyration(sc.Ix, sc.A) * mm:.2f} mm, "
          f"r_y = {radius_of_gyration(sc.Iy, sc.A) * mm:.2f} mm; mass {RHO_7075 * sc.A:.3f} kg/m in 7075-T6")
    solid = rect(b, h, 0.0, 0.0)
    print(f"  a solid {b * mm:.0f} x {h * mm:.0f} mm bar of the same Ix would need height "
          f"{(12 * sc.Ix / b) ** (1 / 3) * mm:.1f} mm and area {b * (12 * sc.Ix / b) ** (1 / 3) * mm2:.0f} mm^2 "
          f"(the full 150 mm bar has Ix = {solid.Ix * mm4:.0f} mm^4)")

    # ---- Z stringer
    hz, bz, tz = 0.025, 0.015, 0.0016
    zp = polygon_section(z_section_polygon(hz, bz, tz))
    zc = composite_section(z_section_parts(hz, bz, tz))
    zt = z_section_thin(hz, bz, tz)
    print(f"\nZ stringer h = {hz * mm:.0f} mm, b = {bz * mm:.0f} mm, t = {tz * mm:.1f} mm:")
    show("polygon (Green)", zp)
    show("composite (parallel axis)", zc)
    show("thin-walled", zt)
    I1, I2, th = principal_axes(zp.Ix, zp.Iy, zp.Ixy)
    w, V = principal_axes_eig(zp.Ix, zp.Iy, zp.Ixy)
    print(f"  principal: I_max = {I1 * mm4:.1f}, I_min = {I2 * mm4:.1f} mm^4 at theta_p = {np.degrees(th):.3f} deg; "
          f"eigh: {w[::-1] * mm4} mm^4")
    print(f"  rotated by theta_p: Iuv = {rotate_axes(zp.Ix, zp.Iy, zp.Ixy, th)[2] * mm4:.2e} mm^4")
    print(f"  thin-walled errors: Ix {100 * (zt.Ix / zp.Ix - 1):.2f} %, Iy {100 * (zt.Iy / zp.Iy - 1):.2f} %, "
          f"Ixy {100 * (zt.Ixy / zp.Ixy - 1):.2f} %")

    # ---- angle
    La = angle_polygon(0.030, 0.020, 0.002)
    sa = polygon_section(La)
    print("\nangle 30 x 20 x 2 mm:")
    show("polygon (Green)", sa)
    I1a, I2a, tha = principal_axes(sa.Ix, sa.Iy, sa.Ixy)
    print(f"  principal I = {I1a * mm4:.1f}, {I2a * mm4:.1f} mm^4 at {np.degrees(tha):.3f} deg")

    # ---- tube by polygon convergence
    ro, ri = 0.010, 0.009
    exact = np.pi / 4 * (ro ** 4 - ri ** 4)
    print(f"\ntube {2 * ro * mm:.0f} x {(ro - ri) * mm:.0f} mm: exact Ix = {exact * mm4:.3f} mm^4")
    errs = []
    ns = [8, 16, 32, 64, 128, 256]
    for n in ns:
        s = polygon_section(circle_polygon(ro, n), [circle_polygon(ri, n)])
        errs.append(abs(s.Ix - exact) / exact)
        print(f"  {n:4d}-gon: Ix = {s.Ix * mm4:.4f} mm^4, relative error {errs[-1]:.3e}")
    orders = np.diff(np.log(errs)) / np.log(0.5)
    print(f"  observed orders {np.array2string(orders, precision=3)}")

    # ---- wing box: full composite against booms only
    wb = dict(width=0.40, height=0.12, t_skin=0.0016, t_web=0.002, cap_area=300e-6, n_stringers=4, stringer_area=60e-6)
    full = composite_section(wing_box_parts(**wb))
    booms = composite_section(wing_box_parts(**wb, idealisation="booms"))
    lump = composite_section(wing_box_parts(**wb, idealisation="lumped"))
    print(f"\nwing box 400 x 120 mm, skins 1.6 mm, webs 2 mm, caps 300 mm^2, 4 stringers of 60 mm^2 per skin:")
    print(f"  full Ix = {full.Ix * mm4:.0f} mm^4; booms only {booms.Ix * mm4:.0f} mm^4 "
          f"({100 * (booms.Ix / full.Ix - 1):.1f} %); skin lumped into booms {lump.Ix * mm4:.0f} mm^4 "
          f"({100 * (lump.Ix / full.Ix - 1):.2f} %); mass {RHO_2024 * full.A:.2f} kg/m in 2024-T3")

    # ---- the experiment: thin-walled error vs t/b, principal angle vs flange ratio
    print("\nZ section h = 25 mm, b = 15 mm: thin-walled error vs t")
    for t_ in (0.0005, 0.001, 0.0016, 0.0025, 0.004):
        p_ = polygon_section(z_section_polygon(hz, bz, t_))
        q_ = z_section_thin(hz, bz, t_)
        print(f"  t = {t_ * mm:.1f} mm (t/b = {t_ / bz:.3f}): Ix {100 * (q_.Ix / p_.Ix - 1):+.3f} %, "
              f"Iy {100 * (q_.Iy / p_.Iy - 1):+.3f} %, Ixy {100 * (q_.Ixy / p_.Ixy - 1):+.3f} %")
    print("Z section h = 25 mm, t = 1.6 mm: principal angle vs flange width b")
    bs = np.array([0.005, 0.010, 0.015, 0.020, 0.025, 0.030])
    angs = []
    for b_ in bs:
        p_ = polygon_section(z_section_polygon(hz, b_, tz))
        angs.append(np.degrees(principal_axes(p_.Ix, p_.Iy, p_.Ixy)[2]))
        print(f"  b = {b_ * mm:4.0f} mm: theta_p = {angs[-1]:7.3f} deg, I_min/I_max = "
              f"{principal_axes(p_.Ix, p_.Iy, p_.Ixy)[1] / principal_axes(p_.Ix, p_.Iy, p_.Ixy)[0]:.4f}")

    # ---- weight and balance
    items = {"empty aircraft": (700.0, 2.10), "pilot": (80.0, 2.00), "passenger": (75.0, 2.00),
             "fuel 100 L": (72.0, 2.25), "baggage": (20.0, 3.00)}
    m = np.array([v[0] for v in items.values()])
    x = np.array([v[1] for v in items.values()])
    xcg = float(composite_centroid(m, x))
    print(f"\nweight and balance (illustrative): total {m.sum():.0f} kg, x_cg = {xcg:.4f} m aft of datum, "
          f"{cg_percent_mac(xcg, 1.80, 1.50):.2f} % MAC (LEMAC 1.80 m, MAC 1.50 m)")

    # ---- figure
    fig, axs = plt.subplots(1, 3, figsize=(14, 4.4))
    ax = axs[0]
    for poly, lab in ((i_section_polygon(h, b, tf, tw), "I spar"),):
        P = np.vstack([poly, poly[:1]])
        ax.plot(P[:, 0] * mm, P[:, 1] * mm, "k-")
    ax.plot(sc.cx * mm, sc.cy * mm, "r+", ms=12)
    zpoly = z_section_polygon(hz, bz, tz) + np.array([0.06, 0.0])
    P = np.vstack([zpoly, zpoly[:1]])
    ax.plot(P[:, 0] * mm, P[:, 1] * mm, "b-")
    for ang, ls in ((th, "r--"), (th + np.pi / 2, "g--")):
        e = np.array([np.cos(ang), np.sin(ang)]) * 0.02
        ax.plot((0.06 + np.array([-e[0], e[0]])) * mm, np.array([-e[1], e[1]]) * mm, ls)
    ax.set_aspect("equal"); ax.set_title("I spar and Z stringer (principal axes)"); ax.set_xlabel("mm")
    ax = axs[1]
    ax.loglog(ns, errs, "o-", label="tube n-gon")
    ax.loglog(ns, errs[0] * (np.array(ns) / ns[0]) ** -2.0, "k:", label="slope -2")
    ax.set_xlabel("vertices n"); ax.set_ylabel("relative error in Ix"); ax.legend(); ax.set_title("polygon convergence")
    ax = axs[2]
    ax.plot(bs * mm, angs, "o-")
    ax.set_xlabel("Z flange width b (mm)"); ax.set_ylabel("principal angle (deg)"); ax.set_title("Z stringer")
    fig.tight_layout()
    fig.savefig(out_dir / "section_properties.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
