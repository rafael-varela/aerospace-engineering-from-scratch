"""P01.07 Vectors, Dot & Cross Products: reference solution.

A small vector library for three-dimensional mechanics, written from components so that every
identity can be checked against NumPy's own ``dot``, ``cross`` and ``linalg.norm``. Every
function accepts a single vector of shape (3,) or a stack of vectors of shape (N, 3) and works
row by row.

Applications: the orbital plane, inclination and specific angular momentum h = r x v of a
satellite from one position/velocity pair; the moment (torque) of distributed wing loads about
the root; tangential and normal components of acceleration along a flight path.

Run:  python solution.py          (prints the satellite and wing tables, writes vector_mechanics.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-01/lesson-07-vectors-dot-cross
"""
from __future__ import annotations

from pathlib import Path

import numpy as np

MU_EARTH = 3.986004418e14  # m**3/s**2


def _as_vectors(a: np.ndarray) -> np.ndarray:
    """Return ``a`` as a float array of shape (N, 3); raise ``ValueError`` if the last axis is not 3."""
    arr = np.asarray(a, dtype=float)
    if arr.ndim == 0 or arr.shape[-1] != 3:
        raise ValueError(f"expected vectors with 3 components, got shape {arr.shape}")
    return arr.reshape(-1, 3)


def _pair(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool]:
    """Return both inputs as broadcast (N, 3) arrays plus a flag: were both single vectors?"""
    A, B = np.broadcast_arrays(_as_vectors(a), _as_vectors(b))
    single = np.asarray(a).ndim == 1 and np.asarray(b).ndim == 1
    return A, B, single


def _out(out: np.ndarray, single: bool) -> np.ndarray | float:
    """Drop the leading axis for single-vector input; return a plain float for a single scalar."""
    if single:
        out = out[0]
        return float(out) if np.ndim(out) == 0 else out
    return out


def dot(a: np.ndarray, b: np.ndarray) -> np.ndarray | float:
    """Return the dot product a . b = a1 b1 + a2 b2 + a3 b3 (row by row for (N, 3) inputs)."""
    A, B, single = _pair(a, b)
    return _out(A[:, 0] * B[:, 0] + A[:, 1] * B[:, 1] + A[:, 2] * B[:, 2], single)


def cross(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Return the cross product a x b from the component formula (right-handed axes)."""
    A, B, single = _pair(a, b)
    out = np.column_stack([
        A[:, 1] * B[:, 2] - A[:, 2] * B[:, 1],
        A[:, 2] * B[:, 0] - A[:, 0] * B[:, 2],
        A[:, 0] * B[:, 1] - A[:, 1] * B[:, 0],
    ])
    return _out(out, single)


def norm(a: np.ndarray) -> np.ndarray | float:
    """Return the Euclidean length |a| = sqrt(a . a)."""
    A = _as_vectors(a)
    return _out(np.sqrt(np.sum(A * A, axis=1)), np.asarray(a).ndim == 1)


def unit(a: np.ndarray) -> np.ndarray:
    """Return a / |a|. Raise ``ValueError`` if any vector has zero length."""
    A = _as_vectors(a)
    n = np.sqrt(np.sum(A * A, axis=1))
    if np.any(n == 0.0):
        raise ValueError("cannot normalise a zero vector")
    return _out(A / n[:, None], np.asarray(a).ndim == 1)


def angle_between(a: np.ndarray, b: np.ndarray) -> np.ndarray | float:
    """Return the angle (rad, in [0, pi]) between ``a`` and ``b`` using atan2(|a x b|, a . b).

    The atan2 form is accurate for nearly parallel vectors, where arccos(a . b / (|a||b|))
    loses digits. Raise ``ValueError`` if either vector is zero.
    """
    A, B, single = _pair(a, b)
    if np.any(np.sum(A * A, axis=1) == 0.0) or np.any(np.sum(B * B, axis=1) == 0.0):
        raise ValueError("angle with a zero vector is undefined")
    c = _as_vectors(cross(A, B))
    return _out(np.arctan2(np.sqrt(np.sum(c * c, axis=1)), np.sum(A * B, axis=1)), single)


def projection(a: np.ndarray, onto: np.ndarray) -> np.ndarray:
    """Return the vector projection of ``a`` onto ``onto``: (a . u) u with u = onto/|onto|."""
    A, U, single = _pair(a, unit(onto))
    return _out(np.sum(A * U, axis=1)[:, None] * U, single)


def rejection(a: np.ndarray, from_vec: np.ndarray) -> np.ndarray:
    """Return the component of ``a`` perpendicular to ``from_vec``: a - projection(a, from_vec)."""
    A, P, single = _pair(a, projection(a, from_vec))
    return _out(A - P, single)


def scalar_triple(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray | float:
    """Return a . (b x c): the signed volume of the parallelepiped spanned by a, b, c."""
    return dot(a, cross(b, c))


def moment(r: np.ndarray, force: np.ndarray) -> np.ndarray:
    """Return the moment M = r x F (N m) of force ``force`` applied at position ``r`` about the origin."""
    return cross(r, force)


def resultant(points: np.ndarray, forces: np.ndarray, about: np.ndarray = (0.0, 0.0, 0.0)) -> tuple[np.ndarray, np.ndarray]:
    """Return (F_total, M_total) of forces ``forces`` (N, 3) applied at ``points`` (N, 3), moments about ``about``.

    M_total = sum_i (r_i - about) x F_i. Raise ``ValueError`` if the two arrays differ in length.
    """
    P, F = _as_vectors(points), _as_vectors(forces)
    if P.shape[0] != F.shape[0]:
        raise ValueError("points and forces must have the same number of rows")
    arm = P - np.asarray(about, dtype=float)
    return F.sum(axis=0), _as_vectors(cross(arm, F)).sum(axis=0)


def specific_angular_momentum(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Return h = r x v (m**2/s), the specific angular momentum of a point mass at ``r`` moving with ``v``."""
    return cross(r, v)


def orbital_plane(r: np.ndarray, v: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (unit normal to the orbital plane, inclination in rad) from one state vector.

    The normal is h/|h| with h = r x v; the inclination is the angle between the normal and the
    +z axis (Earth's spin axis), arccos(h_z/|h|). Raise ``ValueError`` if r and v are parallel.
    """
    h = specific_angular_momentum(r, v)
    hn = norm(h)
    if hn == 0.0:
        raise ValueError("r and v are parallel: the orbital plane is undefined")
    n = h / hn
    return n, float(np.arctan2(np.hypot(n[0], n[1]), n[2]))


def semi_latus_rectum(r: np.ndarray, v: np.ndarray, mu: float = MU_EARTH) -> float:
    """Return p = |h|**2 / mu (m), linking this lesson to the conic r = p/(1 + e cos theta) of P01.06."""
    return float(norm(specific_angular_momentum(r, v)) ** 2 / mu)


def flight_path_angle(r: np.ndarray, v: np.ndarray) -> float:
    """Return the angle (rad) between ``v`` and the local horizontal (the plane perpendicular to ``r``).

    gamma = atan2(v . r_hat, |v - (v . r_hat) r_hat|); positive when climbing.
    """
    rhat = unit(r)
    vr = dot(v, rhat)
    vt = norm(rejection(v, r))
    return float(np.arctan2(vr, vt))


def tangential_normal_components(v: np.ndarray, a: np.ndarray) -> tuple[float, float]:
    """Return (a_t, a_n): acceleration components along and perpendicular to the velocity.

    a_t = a . v_hat (rate of change of speed), a_n = |a x v_hat| (centripetal, v**2/R). Raise
    ``ValueError`` for zero velocity.
    """
    vhat = unit(v)
    return float(dot(a, vhat)), float(norm(cross(a, vhat)))


def point_line_distance(point: np.ndarray, line_point: np.ndarray, direction: np.ndarray) -> float:
    """Return the shortest distance from ``point`` to the line through ``line_point`` with ``direction``."""
    d = unit(direction)
    return float(norm(cross(np.asarray(point, float) - np.asarray(line_point, float), d)))


def plane_from_points(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (unit normal n, offset d) of the plane n . x = d through three points. Raise ``ValueError`` if collinear."""
    p1, p2, p3 = (np.asarray(p, dtype=float) for p in (p1, p2, p3))
    n = cross(p2 - p1, p3 - p1)
    if norm(n) == 0.0:
        raise ValueError("the three points are collinear")
    n = unit(n)
    return n, float(dot(n, p1))


def point_plane_distance(point: np.ndarray, normal: np.ndarray, d: float) -> float:
    """Return the signed distance from ``point`` to the plane n . x = d (positive on the side ``normal`` points to)."""
    n = unit(normal)
    return float(dot(point, n) - d / norm(normal))


def elliptic_lift_distribution(y: np.ndarray, lift_total: float, semi_span: float) -> np.ndarray:
    """Return the lift per unit span (N/m) of an elliptic distribution at stations ``y`` on a half-wing.

    l(y) = (4 L / (pi b)) sqrt(1 - (2y/b)**2) with b = 2 * semi_span, so that the integral over the
    whole span [-b/2, b/2] equals ``lift_total``. Raise ``ValueError`` if any |y| > semi_span.
    """
    y = np.asarray(y, dtype=float)
    if np.any(np.abs(y) > semi_span) or semi_span <= 0.0:
        raise ValueError("stations must lie within the semi-span")
    b = 2.0 * semi_span
    return 4.0 * lift_total / (np.pi * b) * np.sqrt(1.0 - (2.0 * y / b) ** 2)


def wing_root_loads(lift_total: float, semi_span: float, sweep_rad: float, n_strips: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Return (F, M) at the root of a half-wing carrying half of ``lift_total`` with an elliptic distribution.

    Axes: x aft, y along the right wing, z up. Strip loads act at the quarter-chord line
    x = y tan(sweep). M_x is the root bending moment, M_y the root torsion. The half-wing is split
    into ``n_strips`` equal strips with the load applied at each strip's mid-span station.
    """
    if n_strips < 1:
        raise ValueError("need at least one strip")
    edges = np.linspace(0.0, semi_span, n_strips + 1)
    y = 0.5 * (edges[:-1] + edges[1:])
    dy = np.diff(edges)
    lift = elliptic_lift_distribution(y, lift_total, semi_span) * dy
    points = np.column_stack([y * np.tan(sweep_rad), y, np.zeros_like(y)])
    forces = np.column_stack([np.zeros_like(y), np.zeros_like(y), lift])
    return resultant(points, forces)


def main(out_dir: str | Path = ".") -> None:
    """Print the satellite and wing results and save vector_mechanics.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # One tracked state of a low Earth orbit satellite (Earth-centred axes, z = spin axis), SI units.
    r = np.array([3500.0e3, -5000.0e3, 3000.0e3])
    v = np.array([6.020e3, 4.670e3, 0.730e3])
    h = specific_angular_momentum(r, v)
    n, inc = orbital_plane(r, v)
    print("Satellite state vector")
    print(f"  |r| = {norm(r)/1e3:.1f} km (altitude {(norm(r) - 6378.137e3)/1e3:.1f} km), |v| = {norm(v):.1f} m/s")
    print(f"  h = r x v = {h/1e9} x 1e9 m^2/s, |h| = {norm(h):.4e} m^2/s")
    print(f"  orbital plane normal = {np.round(n, 5)}, inclination = {np.degrees(inc):.3f} deg")
    print(f"  r . h = {dot(r, h):.3e}, v . h = {dot(v, h):.3e}  (both zero: r and v lie in the plane)")
    print(f"  flight-path angle = {np.degrees(flight_path_angle(r, v)):.3f} deg, p = |h|^2/mu = {semi_latus_rectum(r, v)/1e3:.1f} km")
    a_grav = -MU_EARTH / norm(r) ** 3 * r
    at, an = tangential_normal_components(v, a_grav)
    print(f"  gravity {norm(a_grav):.4f} m/s^2: tangential {at:+.4f}, normal {an:.4f} m/s^2  -> turning radius v^2/a_n = {norm(v)**2/an/1e3:.1f} km")
    print(f"  numpy check: max |cross - np.cross| = {np.max(np.abs(h - np.cross(r, v))):.1e}, |dot - np.dot| = {abs(dot(r, v) - np.dot(r, v)):.1e}")

    # A 16 m span wing lifting 15.5 kN with 25 deg of quarter-chord sweep: root loads of the right half.
    L_total, semi_span, sweep = 15.5e3, 8.0, np.radians(25.0)
    F, M = wing_root_loads(L_total, semi_span, sweep, 50)
    print("\nRight half-wing root loads (elliptic lift, 50 strips)")
    print(f"  F = {np.round(F, 1)} N  (should be 0, 0, {L_total/2:.0f})")
    print(f"  M = {np.round(M, 1)} N m: bending M_x = {M[0]:.1f}, torsion M_y = {M[1]:.1f}")
    print(f"  analytic bending (L/2)(4 b/2 / (3 pi)) = {L_total/2 * 4*semi_span/(3*np.pi):.1f} N m; "
          f"torsion = -tan(sweep) * bending = {-np.tan(sweep) * L_total/2 * 4*semi_span/(3*np.pi):.1f} N m")

    fig = plt.figure(figsize=(10, 4.5))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    s = 1e-6
    for vec, col, lab in ((r * s, "k", "r"), (v * 1e-3, "b", "v (km/s)"), (n * 5, "r", "h / |h| x 5")):
        ax.quiver(0, 0, 0, *vec, color=col, label=lab)
    th = np.linspace(0, 2 * np.pi, 200)
    e1 = unit(r)
    e2 = unit(rejection(v, r))
    circ = norm(r) * s * (np.outer(np.cos(th), e1) + np.outer(np.sin(th), e2))
    ax.plot(circ[:, 0], circ[:, 1], circ[:, 2], "g--", lw=0.8, label="orbital plane (circle |r|)")
    ax.set_xlabel("x [1000 km]")
    ax.set_ylabel("y [1000 km]")
    ax.set_zlabel("z [1000 km]")
    ax.set_title(f"Orbital plane, i = {np.degrees(inc):.1f} deg")
    ax.legend(fontsize=7)
    ax2 = fig.add_subplot(1, 2, 2)
    y = np.linspace(0, semi_span, 200)
    ax2.plot(y, elliptic_lift_distribution(y, L_total, semi_span), label="l(y) [N/m]")
    ax2.set_xlabel("spanwise station y [m]")
    ax2.set_ylabel("lift per unit span [N/m]")
    ax2.set_title(f"Elliptic lift; root bending {M[0]/1e3:.1f} kN m, torsion {M[1]/1e3:.1f} kN m")
    ax2.grid(True)
    ax2.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "vector_mechanics.png", dpi=150)
    print(f"saved {out_dir / 'vector_mechanics.png'}")


if __name__ == "__main__":
    main()
