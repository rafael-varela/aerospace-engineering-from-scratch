"""P04.09 Magnetism & Induction: reference solution.

A magnetorquer is a current loop (or a coil on a ferromagnetic rod) that a small satellite uses to
turn itself by pushing against the Earth's magnetic field. This module builds the physics of that
device from the Lorentz force up:

* the force on a moving charge and on a current element (Lorentz force, F = q v x B);
* the field of a current: Biot-Savart summed over a polyline of straight segments, the closed form on
  the axis of a circular loop (the limit every discretisation must reproduce), the on-axis field of a
  finite solenoid and its infinite-solenoid limit mu0 n I;
* the Earth's field as a centred, untilted magnetic dipole, and that field along a circular orbit;
* the torque on a current loop, tau = m x B, with m = N I A, and why one coil cannot torque about
  the field direction;
* Faraday's law and Lenz's law: the emf of a rotating coil, the back-emf of a motor, and the RL time
  constant that limits how fast a torquer coil can switch;
* a one-axis "B-dot" detumble simulation whose unsaturated phase has the closed form
  omega(t) = omega0 exp(-t / T_d), T_d = I / (k B^2).

Units: SI throughout. Field in tesla (1 T = 10^4 gauss); the Earth's surface field is a few tens of
microtesla. Coordinates: right-handed, dipole axis along z.

Run:  python solution.py [out_dir]     (prints the verification tables, writes magnetorquer.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-04/lesson-09-magnetism-induction
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import quad, solve_ivp

MU0 = 4.0e-7 * np.pi                 # T m / A, vacuum permeability (exact to 1e-10 relative)
R_EARTH = 6.371e6                    # m, mean Earth radius
B0_EQUATOR = 3.1e-5                  # T, surface equatorial field of the untilted centred-dipole model
EARTH_DIPOLE_MOMENT = B0_EQUATOR * R_EARTH ** 3 / (MU0 / (4.0 * np.pi))   # A m^2, about 8.0e22
RHO_COPPER = 1.68e-8                 # ohm m, copper resistivity at 20 C
Q_PROTON = 1.602176634e-19           # C
M_PROTON = 1.67262192e-27            # kg


# ---------------------------------------------------------------- the Lorentz force

def lorentz_force(q: float, v: np.ndarray, E: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Force on a charge q moving at v in fields E and B: F = q (E + v x B). All vectors shape (3,)."""
    v, E, B = (np.asarray(a, dtype=float) for a in (v, E, B))
    return q * (E + np.cross(v, B))


def gyroradius(m: float, v_perp: float, q: float, B: float) -> float:
    """Radius of the circle a charge of mass m and charge q describes with speed v_perp across a
    uniform field B: r = m v_perp / (|q| B), from q v B = m v^2 / r."""
    return m * v_perp / (abs(q) * B)


def cyclotron_frequency(q: float, B: float, m: float) -> float:
    """Angular frequency of that circular motion, omega_c = |q| B / m [rad/s]: independent of speed."""
    return abs(q) * B / m


# ---------------------------------------------------------------- Biot-Savart

def wire_field(points: np.ndarray, I: float, r: np.ndarray, closed: bool = False) -> np.ndarray:
    """Magnetic field at r of a current I flowing along the polyline through ``points`` (shape (n, 3)),
    by the Biot-Savart law summed over straight segments with the midpoint rule:

        B(r) = (mu0 I / 4 pi) sum_k  dl_k x (r - r_k) / |r - r_k|^3,

    dl_k the segment vector and r_k its midpoint. ``closed=True`` adds the segment from the last point
    back to the first. ``r`` has shape (3,) or (m, 3); the result has the same shape. The error of the
    midpoint rule is O(h^2) in the segment length h away from the wire; it blows up on the wire itself,
    where the true field of a thin wire is infinite anyway.
    """
    P = np.asarray(points, dtype=float)
    if closed:
        P = np.vstack([P, P[:1]])
    dl = P[1:] - P[:-1]                                    # (n, 3)
    mid = 0.5 * (P[1:] + P[:-1])                           # (n, 3)
    r = np.asarray(r, dtype=float)
    single = r.ndim == 1
    R = np.atleast_2d(r)[:, None, :] - mid[None, :, :]     # (m, n, 3)
    d3 = np.linalg.norm(R, axis=-1) ** 3                   # (m, n)
    B = (MU0 * I / (4.0 * np.pi)) * np.sum(np.cross(dl[None, :, :], R) / d3[..., None], axis=1)
    return B[0] if single else B


def loop_points(R: float, n: int, center: np.ndarray | None = None, normal: np.ndarray | None = None) -> np.ndarray:
    """Vertices of a regular n-gon inscribed in a circle of radius R, centred at ``center`` (default the
    origin) with its axis along ``normal`` (default +z); current flows anticlockwise seen from +normal,
    so the loop's magnetic moment is along +normal."""
    phi = 2.0 * np.pi * np.arange(n) / n
    pts = np.column_stack([R * np.cos(phi), R * np.sin(phi), np.zeros(n)])
    if normal is not None:
        nrm = np.asarray(normal, dtype=float)
        nrm = nrm / np.linalg.norm(nrm)
        z = np.array([0.0, 0.0, 1.0])
        v = np.cross(z, nrm)
        s, c = np.linalg.norm(v), float(np.dot(z, nrm))
        if s > 1e-12:
            vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
            Rm = np.eye(3) + vx + vx @ vx * ((1 - c) / s ** 2)
            pts = pts @ Rm.T
        elif c < 0:
            pts = pts * np.array([1.0, -1.0, -1.0])
    if center is not None:
        pts = pts + np.asarray(center, dtype=float)
    return pts


def loop_field(r: np.ndarray, R: float, I: float, n_seg: int = 360) -> np.ndarray:
    """Field at r (shape (3,) or (m, 3)) of a circular loop of radius R in the z = 0 plane centred on the
    origin carrying I anticlockwise (moment +z), by ``wire_field`` on an n_seg-gon."""
    return wire_field(loop_points(R, n_seg), I, r, closed=True)


def loop_field_axis(z: np.ndarray | float, R: float, I: float) -> np.ndarray | float:
    """Closed form on the axis of the loop: B_z(z) = mu0 I R^2 / (2 (R^2 + z^2)^(3/2)).

    At z = 0 this is mu0 I / (2 R); for |z| >> R it is mu0 m / (2 pi |z|^3) with m = I pi R^2, the
    dipole far field on the axis."""
    z = np.asarray(z, dtype=float)
    return MU0 * I * R ** 2 / (2.0 * (R ** 2 + z ** 2) ** 1.5)


def loop_field_quad(rho: float, z: float, R: float, I: float) -> tuple[float, float]:
    """(B_rho, B_z) of the circular loop at cylindrical (rho, z), by ``scipy.integrate.quad`` on the exact
    Biot-Savart integrand over the loop angle phi':

        B_z   = (mu0 I / 4 pi) int R (R - rho cos phi') / d^3 dphi',
        B_rho = (mu0 I / 4 pi) int R z cos phi' / d^3 dphi',    d^2 = rho^2 + R^2 - 2 rho R cos phi' + z^2.

    The field point is taken at (rho, 0, z); by symmetry B_phi = 0."""
    def d3(p):
        return (rho ** 2 + R ** 2 - 2.0 * rho * R * np.cos(p) + z ** 2) ** 1.5
    bz = quad(lambda p: R * (R - rho * np.cos(p)) / d3(p), 0.0, 2.0 * np.pi, limit=200)[0]
    br = quad(lambda p: R * z * np.cos(p) / d3(p), 0.0, 2.0 * np.pi, limit=200)[0]
    k = MU0 * I / (4.0 * np.pi)
    return k * br, k * bz


def straight_wire_field(d: float, I: float) -> float:
    """Magnitude of the field at distance d from an infinite straight wire: mu0 I / (2 pi d), from
    Ampere's law (or from Biot-Savart integrated along the wire)."""
    return MU0 * I / (2.0 * np.pi * d)


def solenoid_field_axis(z: np.ndarray | float, R: float, I: float, n_per_m: float, L: float) -> np.ndarray | float:
    """On-axis field of a finite solenoid of radius R, length L (centred at z = 0), n_per_m turns per
    metre, obtained by integrating ``loop_field_axis`` over the windings:

        B_z(z) = (mu0 n I / 2) [ (z + L/2) / sqrt((z + L/2)^2 + R^2) - (z - L/2) / sqrt((z - L/2)^2 + R^2) ].

    At the centre of a long solenoid (L >> R) this tends to mu0 n I; at either end to half of it."""
    z = np.asarray(z, dtype=float)
    a, b = z + 0.5 * L, z - 0.5 * L
    return 0.5 * MU0 * n_per_m * I * (a / np.sqrt(a ** 2 + R ** 2) - b / np.sqrt(b ** 2 + R ** 2))


def solenoid_inductance(N: int, A: float, L: float) -> float:
    """Self-inductance of a long air-core solenoid, L_ind = mu0 N^2 A / L [H] (flux linkage N Phi per
    ampere with B = mu0 N I / L inside). An overestimate for a short coil, since the end field is weaker."""
    return MU0 * N ** 2 * A / L


# ---------------------------------------------------------------- the Earth's dipole field

def dipole_field(r: np.ndarray, m: np.ndarray) -> np.ndarray:
    """Field of a point magnetic dipole m at position r (both shape (3,), or r of shape (n, 3)):

        B = (mu0 / 4 pi) [ 3 (m . r_hat) r_hat - m ] / |r|^3.
    """
    r = np.asarray(r, dtype=float)
    m = np.asarray(m, dtype=float)
    single = r.ndim == 1
    rr = np.atleast_2d(r)
    dist = np.linalg.norm(rr, axis=1, keepdims=True)
    rhat = rr / dist
    B = (MU0 / (4.0 * np.pi)) * (3.0 * (rhat @ m)[:, None] * rhat - m) / dist ** 3
    return B[0] if single else B


def earth_dipole_field(r: np.ndarray) -> np.ndarray:
    """The Earth's field in the untilted centred-dipole model: ``dipole_field`` with
    m = -EARTH_DIPOLE_MOMENT z_hat (the moment points south, so the field at the north pole points
    down into the ground and the field at the equator points north, +z). Surface equatorial magnitude
    B0_EQUATOR, polar magnitude 2 B0_EQUATOR, falling as 1/r^3.

    Model limits: the real dipole is tilted about 11 degrees from the spin axis and offset from the
    centre, and the non-dipole part of the field is of order 10 % at the surface; use this for
    magnitudes and torque budgets, not for navigation."""
    return dipole_field(r, np.array([0.0, 0.0, -EARTH_DIPOLE_MOMENT]))


def orbit_positions(altitude: float, inclination_deg: float, n: int = 361) -> np.ndarray:
    """Positions (n, 3) along one circular orbit of the given altitude and inclination in the dipole
    frame, parametrised by the argument of latitude u in [0, 2 pi]: r = a [cos u, cos i sin u, sin i sin u]
    with a = R_EARTH + altitude. The ascending node is on the x axis."""
    a = R_EARTH + altitude
    i = np.radians(inclination_deg)
    u = np.linspace(0.0, 2.0 * np.pi, n)
    return a * np.column_stack([np.cos(u), np.cos(i) * np.sin(u), np.sin(i) * np.sin(u)])


def field_along_orbit(altitude: float, inclination_deg: float, n: int = 361) -> tuple[np.ndarray, np.ndarray]:
    """(u, B) with u the argument of latitude (n,) and B the dipole field (n, 3) along the orbit of
    ``orbit_positions``."""
    pos = orbit_positions(altitude, inclination_deg, n)
    return np.linspace(0.0, 2.0 * np.pi, n), earth_dipole_field(pos)


# ---------------------------------------------------------------- torque on a coil

def magnetic_moment(N: int, I: float, A: float) -> float:
    """Magnetic moment of a flat coil of N turns, current I, area A: m = N I A [A m^2]."""
    return N * I * A


def torque_on_dipole(m: np.ndarray, B: np.ndarray) -> np.ndarray:
    """Torque on a magnetic moment m in a field B: tau = m x B [N m]. Its component along B is zero,
    so a magnetorquer can never produce torque about the local field direction."""
    return np.cross(np.asarray(m, dtype=float), np.asarray(B, dtype=float))


def coil_resistance(N: int, side: float, wire_diameter: float, rho: float = RHO_COPPER) -> float:
    """Resistance of N square turns of side ``side`` wound from wire of the given diameter:
    R = rho (4 N side) / (pi d^2 / 4)."""
    length = 4.0 * N * side
    area = np.pi * wire_diameter ** 2 / 4.0
    return rho * length / area


@dataclass
class Torquer:
    """An air-core magnetorquer coil on one face of a CubeSat."""
    N: int = 200                   # turns
    side: float = 0.08             # m, square coil
    wire_diameter: float = 0.25e-3 # m
    I: float = 0.12                # A, nominal current

    @property
    def area(self) -> float:
        return self.side ** 2

    @property
    def moment(self) -> float:
        return magnetic_moment(self.N, self.I, self.area)

    @property
    def resistance(self) -> float:
        return coil_resistance(self.N, self.side, self.wire_diameter)

    @property
    def power(self) -> float:
        return self.I ** 2 * self.resistance

    @property
    def inductance(self) -> float:
        # a flat coil is not a long solenoid; use the solenoid formula with an effective length equal
        # to the side as an order-of-magnitude estimate (the lesson says so)
        return solenoid_inductance(self.N, self.area, self.side)

    @property
    def time_constant(self) -> float:
        return self.inductance / self.resistance


# ---------------------------------------------------------------- Faraday, Lenz, RL

def faraday_emf(t: np.ndarray, flux: np.ndarray) -> np.ndarray:
    """emf = -dPhi/dt by central differences on sampled flux linkage (Faraday's law with Lenz's sign:
    the emf drives a current whose own field opposes the change in flux)."""
    return -np.gradient(np.asarray(flux, dtype=float), np.asarray(t, dtype=float))


def generator_emf(t: np.ndarray | float, N: int, A: float, B: float, omega: float) -> np.ndarray | float:
    """emf of a coil of N turns and area A rotating at omega in a uniform field B, from
    Phi = N A B cos(omega t): emf = N A B omega sin(omega t). Peak value N A B omega."""
    return N * A * B * omega * np.sin(omega * np.asarray(t, dtype=float))


def motor_back_emf(k_e: float, omega: float) -> float:
    """Back-emf of a motor with voltage constant k_e [V s/rad] turning at omega: e = k_e omega. The
    same constant gives the torque per ampere (k_t = k_e in SI), so the no-load speed is V / k_e."""
    return k_e * omega


def rl_current(t: np.ndarray | float, V: float, R: float, L: float, i0: float = 0.0) -> np.ndarray | float:
    """Current in a series RL circuit after a step V is applied at t = 0 with initial current i0:
    i(t) = V/R + (i0 - V/R) exp(-t / tau), tau = L / R. Initial slope (V - i0 R) / L."""
    tau = L / R
    return V / R + (i0 - V / R) * np.exp(-np.asarray(t, dtype=float) / tau)


def rl_simulate(t_end: float, V_fn: Callable[[float], float], R: float, L: float, i0: float = 0.0,
                n: int = 401) -> tuple[np.ndarray, np.ndarray]:
    """Integrate L di/dt = V(t) - R i with ``scipy.integrate.solve_ivp`` (RK45, tight tolerances) for an
    arbitrary source V(t). Returns (t, i) on n points."""
    t = np.linspace(0.0, t_end, n)
    sol = solve_ivp(lambda tt, y: [(V_fn(tt) - R * y[0]) / L], (0.0, t_end), [i0], t_eval=t,
                    rtol=1e-10, atol=1e-13)
    return t, sol.y[0]


# ---------------------------------------------------------------- one-axis B-dot detumble

def detumble_time_constant(I_sc: float, k: float, B: float) -> float:
    """Time constant of the unsaturated B-dot law about an axis perpendicular to B:
    with m = k omega B the torque is -k B^2 omega, so I omega' = -k B^2 omega and T_d = I / (k B^2)."""
    return I_sc / (k * B ** 2)


def saturation_rate(m_max: float, k: float, B: float) -> float:
    """Spin rate below which the commanded moment k omega B stays under the coil limit m_max."""
    return m_max / (k * B)


def bdot_detumble(I_sc: float, omega0: float, B: float, k: float, m_max: float, t_end: float,
                  n: int = 2001) -> tuple[np.ndarray, np.ndarray]:
    """One-axis detumble: a body with moment of inertia I_sc spins at omega about an axis
    perpendicular to a fixed field of magnitude B; the coil moment is m = min(k omega B, m_max) in the
    direction that makes tau = m x B oppose the spin. Then

        I domega/dt = - min(k omega B, m_max) B.

    Saturated phase: linear decay at m_max B / I. Unsaturated: omega = omega_s exp(-(t - t_s) / T_d).
    Returns (t, omega) on n points from solve_ivp."""
    def rhs(tt, y):
        w = y[0]
        m = min(k * abs(w) * B, m_max)
        return [-np.sign(w) * m * B / I_sc]
    t = np.linspace(0.0, t_end, n)
    sol = solve_ivp(rhs, (0.0, t_end), [omega0], t_eval=t, rtol=1e-9, atol=1e-12, max_step=t_end / 400)
    return t, sol.y[0]


def momentum_dump_time(I_sc: float, omega0: float, tau_max: float) -> float:
    """Time to remove the angular momentum I omega0 at a constant torque tau_max: I omega0 / tau_max.
    The lower bound for any detumble that cannot exceed tau_max."""
    return I_sc * omega0 / tau_max


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's numbers, print them and save magnetorquer.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- Lorentz force: a radiation-belt proton
    B_leo = 3.0e-5
    v = 1.0e6
    print(f"proton at {v:.0e} m/s across {B_leo * 1e6:.0f} uT: force {Q_PROTON * v * B_leo:.2e} N, "
          f"gyroradius {gyroradius(M_PROTON, v, Q_PROTON, B_leo):.0f} m, "
          f"cyclotron frequency {cyclotron_frequency(Q_PROTON, B_leo, M_PROTON) / (2 * np.pi):.0f} Hz")

    # ---- Biot-Savart: loop on axis, convergence with segments
    R, I = 0.04, 0.12
    print(f"\nloop R = {R * 100:.0f} cm, I = {I} A: centre field {loop_field_axis(0.0, R, I) * 1e6:.3f} uT, "
          f"far field at z = 10 R: {loop_field_axis(10 * R, R, I) * 1e9:.3f} nT (dipole: "
          f"{MU0 * I * np.pi * R ** 2 / (2 * np.pi * (10 * R) ** 3) * 1e9:.3f} nT)")
    print("segments   |B_z(0,0,R/2) - closed form| / closed form")
    exact = loop_field_axis(0.5 * R, R, I)
    errs = []
    for n in (8, 16, 32, 64, 128):
        b = loop_field(np.array([0.0, 0.0, 0.5 * R]), R, I, n)[2]
        errs.append(abs(b - exact) / exact)
        print(f"{n:8d}   {errs[-1]:.2e}")
    orders = np.diff(np.log(errs)) / np.log(0.5)
    print(f"observed orders: {np.array2string(orders, precision=2)}")
    br, bz = loop_field_quad(0.5 * R, 0.3 * R, R, I)
    bseg = loop_field(np.array([0.5 * R, 0.0, 0.3 * R]), R, I, 360)
    print(f"off axis (rho, z) = (R/2, 0.3 R): quad (B_rho, B_z) = ({br * 1e6:.4f}, {bz * 1e6:.4f}) uT, "
          f"360-gon ({bseg[0] * 1e6:.4f}, {bseg[2] * 1e6:.4f}) uT, B_phi = {bseg[1]:.1e} T")

    # straight wire against Ampere
    pts = np.column_stack([np.zeros(8001), np.zeros(8001), np.linspace(-20.0, 20.0, 8001)])
    bw = wire_field(pts, 1.0, np.array([0.05, 0.0, 0.0]))
    print(f"wire 40 m long in 5 mm segments, 1 A, at 5 cm: {np.linalg.norm(bw) * 1e6:.4f} uT against mu0 I / 2 pi d = "
          f"{straight_wire_field(0.05, 1.0) * 1e6:.4f} uT")

    # solenoid
    n_pm, Ls, Rs = 2000.0, 0.10, 0.005
    print(f"solenoid {n_pm:.0f} turns/m, L = {Ls * 100:.0f} cm, R = {Rs * 1000:.0f} mm, 1 A: centre "
          f"{solenoid_field_axis(0.0, Rs, 1.0, n_pm, Ls) * 1e3:.4f} mT, end {solenoid_field_axis(0.5 * Ls, Rs, 1.0, n_pm, Ls) * 1e3:.4f} mT, "
          f"mu0 n I = {MU0 * n_pm * 1e3:.4f} mT")

    # ---- the Earth's dipole and the orbit
    print(f"\nEarth dipole: m = {EARTH_DIPOLE_MOMENT:.3e} A m^2; surface equator {B0_EQUATOR * 1e6:.1f} uT, pole "
          f"{np.linalg.norm(earth_dipole_field(np.array([0, 0, R_EARTH]))) * 1e6:.1f} uT")
    alt = 500e3
    for inc in (0.0, 51.6, 98.0):
        u, Bo = field_along_orbit(alt, inc)
        mag = np.linalg.norm(Bo, axis=1)
        print(f"  {alt / 1e3:.0f} km, i = {inc:5.1f} deg: |B| min {mag.min() * 1e6:.2f}, max {mag.max() * 1e6:.2f}, "
              f"mean {mag.mean() * 1e6:.2f} uT")

    # ---- the torquer
    tq = Torquer()
    print(f"\ntorquer: N = {tq.N}, side {tq.side * 100:.0f} cm, wire {tq.wire_diameter * 1e3:.2f} mm, I = {tq.I} A: "
          f"m = {tq.moment:.3f} A m^2, R = {tq.resistance:.1f} ohm, P = {tq.power:.2f} W, "
          f"L ~ {tq.inductance * 1e3:.1f} mH, tau = L/R ~ {tq.time_constant * 1e3:.2f} ms")
    u, Bo = field_along_orbit(alt, 51.6)
    m_body = np.array([0.0, 0.0, tq.moment])          # coil axis along the dipole axis z, for the demonstration
    tau = np.cross(m_body, Bo)
    print(f"torque with m along z, 500 km, i = 51.6: max {np.linalg.norm(tau, axis=1).max() * 1e6:.2f} uN m, "
          f"min {np.linalg.norm(tau, axis=1).min() * 1e6:.2f} uN m; component along B: "
          f"{np.max(np.abs(np.sum(tau * Bo, axis=1) / np.linalg.norm(Bo, axis=1))):.1e} N m")

    # ---- Faraday and RL
    omega = 2 * np.pi * 50.0
    tt = np.linspace(0.0, 0.04, 4001)
    flux = tq.N * tq.area * 0.5 * np.cos(omega * tt)
    emf = faraday_emf(tt, flux)
    print(f"\nrotating the coil at 50 Hz in 0.5 T: peak emf {emf.max():.3f} V (N A B omega = {tq.N * tq.area * 0.5 * omega:.3f} V)")
    tau_rl = tq.time_constant
    t_rl, i_rl = rl_simulate(5 * tau_rl, lambda t: tq.I * tq.resistance, tq.resistance, tq.inductance)
    i_cf = rl_current(t_rl, tq.I * tq.resistance, tq.resistance, tq.inductance)
    k63 = np.argmin(np.abs(t_rl - tau_rl))
    print(f"RL step to {tq.I} A: at t = tau the current is {i_rl[k63] / tq.I * 100:.1f} % of final; solve_ivp vs closed form "
          f"max |diff| {np.max(np.abs(i_rl - i_cf)):.1e} A")

    # ---- detumble
    I_sc, w0, Bf, k, m_max = 0.04, np.radians(10.0), 3.0e-5, 5.0e4, tq.moment
    T_d = detumble_time_constant(I_sc, k, Bf)
    w_s = saturation_rate(m_max, k, Bf)
    t_dt, w_dt = bdot_detumble(I_sc, w0, Bf, k, m_max, 6000.0)
    t_sat = (w0 - w_s) * I_sc / (m_max * Bf)
    print(f"\ndetumble: I = {I_sc} kg m^2, omega0 = 10 deg/s, B = {Bf * 1e6:.0f} uT, k = {k:.0e}, m_max = {m_max:.3f} A m^2")
    print(f"  saturated until omega_s = {np.degrees(w_s):.2f} deg/s at t = {t_sat:.0f} s (slope {np.degrees(m_max * Bf / I_sc) * 60:.3f} deg/s per min); "
          f"then T_d = {T_d:.0f} s; lower bound I omega0 / tau_max = {momentum_dump_time(I_sc, w0, m_max * Bf):.0f} s")
    for tk in (1000, 2000, 4000, 6000):
        print(f"  t = {tk:4d} s: omega = {np.degrees(w_dt[np.argmin(np.abs(t_dt - tk))]):.3f} deg/s")

    # ---- figure
    fig, axs = plt.subplots(2, 2, figsize=(11, 8))
    rho = np.linspace(-3 * R, 3 * R, 121)
    zz = np.linspace(-3 * R, 3 * R, 121)
    PR, PZ = np.meshgrid(rho, zz)
    pts3 = np.column_stack([PR.ravel(), np.zeros(PR.size), PZ.ravel()])
    Bm = loop_field(pts3, R, I, 180)
    Bx, Bz = Bm[:, 0].reshape(PR.shape), Bm[:, 2].reshape(PR.shape)
    mag = np.hypot(Bx, Bz)
    ax = axs[0, 0]
    ax.streamplot(rho / R, zz / R, Bx, Bz, color=np.log10(mag * 1e6), cmap="viridis", density=1.2)
    ax.plot([-1, 1], [0, 0], "ro")
    ax.set_xlabel("x / R"); ax.set_ylabel("z / R"); ax.set_title("field of a current loop (Biot-Savart, 180 segments)")
    ax = axs[0, 1]
    for inc in (0.0, 51.6, 98.0):
        u, Bo = field_along_orbit(alt, inc)
        ax.plot(np.degrees(u), np.linalg.norm(Bo, axis=1) * 1e6, label=f"i = {inc} deg")
    ax.set_xlabel("argument of latitude (deg)"); ax.set_ylabel("|B| (uT)"); ax.set_title("dipole field along a 500 km orbit"); ax.legend()
    ax = axs[1, 0]
    ax.plot(t_rl * 1e3, i_rl * 1e3, label="solve_ivp"); ax.plot(t_rl * 1e3, i_cf * 1e3, "--", label="closed form")
    ax.axvline(tau_rl * 1e3, color="gray", lw=0.8); ax.set_xlabel("t (ms)"); ax.set_ylabel("i (mA)"); ax.set_title("RL rise of the coil current"); ax.legend()
    ax = axs[1, 1]
    ax.plot(t_dt, np.degrees(w_dt), label="B-dot, saturated coil")
    mask = t_dt >= t_sat
    ax.plot(t_dt[mask], np.degrees(w_s * np.exp(-(t_dt[mask] - t_sat) / T_d)), "--", label="exp(-(t - t_s)/T_d)")
    ax.set_xlabel("t (s)"); ax.set_ylabel("spin rate (deg/s)"); ax.set_title("one-axis detumble"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "magnetorquer.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
