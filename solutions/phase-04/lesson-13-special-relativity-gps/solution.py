"""P04.13 Special Relativity for GPS: reference solution.

Special relativity (Lorentz factor, time dilation, length contraction, the Lorentz transformation,
velocity addition, the relativistic Doppler effect, relativistic energy and momentum), the weak-field
gravitational time dilation, and the clock-correction budget of a GNSS satellite:

* the velocity term  -v^2 / (2 c^2)                  (about -7.2 us/day for a GPS orbit),
* the gravitational term  (GM/c^2)(1/R_E - 1/a)      (about +45.7 us/day),
* the rotation and J2 terms of the ground clock      (about +0.1 us/day),
* the net factory offset  Delta f / f = -4.465e-10    (about +38.6 us/day),
* the eccentricity correction  -2 sqrt(GM a) e sin E / c^2  (46 ns amplitude at e = 0.02),
* the Sagnac term  omega_E (x_s y_r - y_s x_r) / c^2       (up to 133 ns),
* the pseudorange error that accumulates without the correction (11.6 km per day).

Oracle: the net frequency offset of a circular GPS orbit is -3GM/(2 a c^2) - Phi_0/c^2 with Phi_0 the
geoid potential; both sides are computed independently and compared in the tests.

Run:  python solution.py          (prints the clock budget tables, writes gps_relativity.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-04/lesson-13-special-relativity-gps
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

# ---------------------------------------------------------------- constants (SI)

C = 299_792_458.0                 # m/s, speed of light (exact)
GM_EARTH = 3.986004418e14         # m^3/s^2, geocentric gravitational constant
R_EARTH = 6.378137e6              # m, equatorial radius (WGS 84)
OMEGA_EARTH = 7.292115e-5         # rad/s, Earth rotation rate
J2_EARTH = 1.0826e-3              # Earth's oblateness coefficient
GPS_SEMI_MAJOR_AXIS = 26_561.75e3 # m, nominal GPS orbit (period half a sidereal day)
GPS_CLOCK_HZ = 10.23e6            # Hz, nominal satellite clock frequency
GPS_L1_HZ = 1575.42e6             # Hz, L1 carrier
SECONDS_PER_DAY = 86_400.0


# ---------------------------------------------------------------- kinematics of special relativity

def _check_speed(v: np.ndarray | float) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    if np.any(np.abs(v) >= C):
        raise ValueError("speed must be below c")
    return v


def gamma(v: np.ndarray | float) -> np.ndarray | float:
    """Lorentz factor 1 / sqrt(1 - v^2/c^2). Raises ValueError for |v| >= c. Accepts arrays."""
    v = _check_speed(v)
    return 1.0 / np.sqrt(1.0 - (v / C) ** 2)


def gamma_minus_one(v: np.ndarray | float) -> np.ndarray | float:
    """gamma - 1 computed without cancellation: beta^2 / (sqrt(1 - beta^2) (1 + sqrt(1 - beta^2))).

    At v = 4 km/s, gamma - 1 = 9e-11 and the naive gamma(v) - 1 keeps only six digits; this form keeps
    all sixteen. It is the quantity every clock-rate calculation needs.
    """
    v = _check_speed(v)
    b2 = (v / C) ** 2
    s = np.sqrt(1.0 - b2)
    return b2 / (s * (1.0 + s))


def time_dilation_rate(v: np.ndarray | float) -> np.ndarray | float:
    """Fractional rate of a clock moving at speed v relative to one at rest: 1/gamma - 1 (negative).

    To second order in v/c this is -v^2 / (2 c^2): a moving clock loses v^2/(2c^2) seconds per second.
    """
    g1 = gamma_minus_one(v)
    return -g1 / (1.0 + g1)


def length_contraction(L0: float, v: float) -> float:
    """Length L0 / gamma of an object of rest length L0 moving at speed v along its length."""
    return L0 / gamma(v)


def lorentz_transform(t: np.ndarray | float, x: np.ndarray | float, v: float) -> tuple[np.ndarray, np.ndarray]:
    """Event (t, x) in frame S expressed in frame S' moving at velocity v along +x: (t', x').

    t' = gamma (t - v x / c^2),  x' = gamma (x - v t).
    """
    g = gamma(v)
    t = np.asarray(t, dtype=float)
    x = np.asarray(x, dtype=float)
    return g * (t - v * x / C ** 2), g * (x - v * t)


def lorentz_inverse(tp: np.ndarray | float, xp: np.ndarray | float, v: float) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of ``lorentz_transform``: the same boost with -v."""
    return lorentz_transform(tp, xp, -v)


def spacetime_interval(dt: np.ndarray | float, dx: np.ndarray | float) -> np.ndarray | float:
    """The invariant interval s^2 = c^2 dt^2 - dx^2 (same value in every inertial frame)."""
    return (C * np.asarray(dt, dtype=float)) ** 2 - np.asarray(dx, dtype=float) ** 2


def velocity_addition(u: float, v: float) -> float:
    """Velocity of an object moving at u in S', where S' moves at v relative to S: (u + v) / (1 + u v / c^2)."""
    _check_speed(v)
    if abs(u) > C:
        raise ValueError("speed must not exceed c")
    return (u + v) / (1.0 + u * v / C ** 2)


def relativistic_doppler(f0: float, v: float, theta: float = 0.0) -> float:
    """Frequency received from a source of rest frequency f0 moving at speed v.

    ``theta`` is the angle, in the receiver's frame, between the source velocity and the line from the
    source to the receiver: 0 approaching, pi receding, pi/2 transverse.
    f = f0 / (gamma (1 - (v/c) cos theta)).
    """
    return f0 / (gamma(v) * (1.0 - (v / C) * np.cos(theta)))


# ---------------------------------------------------------------- relativistic dynamics

def relativistic_energy(m: float, v: float) -> float:
    """Total energy gamma m c^2 [J]."""
    return gamma(v) * m * C ** 2


def kinetic_energy(m: float, v: float) -> float:
    """Kinetic energy (gamma - 1) m c^2 [J], computed without cancellation at low speed."""
    return gamma_minus_one(v) * m * C ** 2


def momentum(m: float, v: float) -> float:
    """Relativistic momentum gamma m v [kg m/s]."""
    return gamma(v) * m * v


# ---------------------------------------------------------------- gravity and the GPS budget

def circular_orbit_speed(a: float, mu: float = GM_EARTH) -> float:
    """Speed sqrt(mu / a) of a circular orbit of radius a [m]."""
    return float(np.sqrt(mu / a))


def gravitational_rate(r: np.ndarray | float, mu: float = GM_EARTH) -> np.ndarray | float:
    """Weak-field fractional rate of a clock at rest at radius r relative to one far away: -GM / (r c^2).

    Valid when GM/(r c^2) << 1 (here ~7e-10). A clock deeper in the potential runs slower.
    """
    return -mu / (np.asarray(r, dtype=float) * C ** 2)


def geoid_potential_ratio(mu: float = GM_EARTH, R: float = R_EARTH, omega: float = OMEGA_EARTH,
                          J2: float = J2_EARTH) -> float:
    """Phi_0 / c^2 for a clock on the rotating geoid at the equator: -(GM/R)(1 + J2/2) - omega^2 R^2 / 2, over c^2.

    The same value applies anywhere on the geoid (that is what defines it), which is why every ground
    clock at sea level ticks at the same rate.
    """
    phi = -(mu / R) * (1.0 + J2 / 2.0) - 0.5 * omega ** 2 * R ** 2
    return phi / C ** 2


@dataclass
class ClockBudget:
    """Fractional clock-rate terms of a satellite clock relative to a clock on the geoid (dimensionless),
    and the same terms in microseconds per day."""
    velocity: float            # -v^2 / (2 c^2)
    gravity: float             # (GM/c^2)(1/R - 1/a)
    ground_rotation: float     # + omega^2 R^2 / (2 c^2)  (the ground clock also moves)
    ground_j2: float           # + (GM/(R c^2)) J2 / 2     (the geoid is not a sphere)

    @property
    def net(self) -> float:
        return self.velocity + self.gravity + self.ground_rotation + self.ground_j2

    def us_per_day(self) -> dict[str, float]:
        f = 1e6 * SECONDS_PER_DAY
        return {"velocity": self.velocity * f, "gravity": self.gravity * f,
                "ground_rotation": self.ground_rotation * f, "ground_j2": self.ground_j2 * f,
                "net": self.net * f}


def gps_clock_budget(a: float = GPS_SEMI_MAJOR_AXIS, mu: float = GM_EARTH, R: float = R_EARTH,
                     omega: float = OMEGA_EARTH, J2: float = J2_EARTH) -> ClockBudget:
    """Clock-rate budget of a satellite on a circular orbit of radius a relative to a clock on the geoid.

    velocity:  the satellite clock runs slow by v^2/(2c^2) with v = sqrt(mu/a);
    gravity:   it runs fast because it is higher in the potential, by (mu/c^2)(1/R - 1/a);
    ground:    the ground clock itself is moving at omega R and sits on an oblate geoid; both make the
               ground clock slightly slower, so both add to the satellite's apparent gain.
    """
    v = circular_orbit_speed(a, mu)
    return ClockBudget(velocity=-v ** 2 / (2.0 * C ** 2),
                       gravity=(mu / C ** 2) * (1.0 / R - 1.0 / a),
                       ground_rotation=0.5 * (omega * R) ** 2 / C ** 2,
                       ground_j2=(mu / (R * C ** 2)) * J2 / 2.0)


def gps_frequency_offset_oracle(a: float = GPS_SEMI_MAJOR_AXIS) -> float:
    """Independent form of the net offset: -3GM/(2 a c^2) - Phi_0/c^2 (the orbit's total energy per unit
    mass over c^2, minus the geoid potential). Equals ClockBudget.net to round-off."""
    return -1.5 * GM_EARTH / (a * C ** 2) - geoid_potential_ratio()


def factory_clock_frequency(nominal_hz: float = GPS_CLOCK_HZ, offset: float | None = None) -> float:
    """Frequency the satellite clock is set to before launch so that it is seen at ``nominal_hz`` from the ground:
    nominal (1 - offset). With the GPS budget this is 10 229 999.995 43 Hz."""
    if offset is None:
        offset = gps_clock_budget().net
    return nominal_hz * (1.0 - offset)


def crossover_radius(mu: float = GM_EARTH, R: float = R_EARTH) -> float:
    """Orbit radius where the velocity and gravitational terms cancel in the simple (spherical, non-rotating)
    model: GM/(2a) = GM (1/R - 1/a)  =>  a = 3R/2. Found numerically by bisection as a check."""
    f = lambda a: -mu / (2.0 * a) + mu * (1.0 / R - 1.0 / a)
    return brentq(f, 1.01 * R, 10.0 * R, xtol=1e-6)


def eccentricity_correction(e: float, E: np.ndarray | float, a: float = GPS_SEMI_MAJOR_AXIS,
                            mu: float = GM_EARTH) -> np.ndarray | float:
    """Periodic clock correction -2 sqrt(mu a) e sin E / c^2 [s] at eccentric anomaly E [rad].

    The orbit-average rate is already in the factory offset; this is the part that varies around the
    orbit because speed and radius vary together. Amplitude 2 sqrt(mu a) e / c^2: 45.8 ns at e = 0.02.
    """
    return -2.0 * np.sqrt(mu * a) * e * np.sin(np.asarray(E, dtype=float)) / C ** 2


def sagnac_correction(r_sat: np.ndarray, r_rec: np.ndarray, omega: float = OMEGA_EARTH) -> float:
    """Sagnac term omega (x_s y_r - y_s x_r) / c^2 [s] for satellite and receiver positions in an
    Earth-centred frame with z along the rotation axis.

    The receiver moves while the signal is in flight; the extra path is omega * (z-component of r_s x r_r) / c.
    Positive when the receiver moves away from the satellite (satellite to the west).
    """
    rs = np.asarray(r_sat, dtype=float)
    rr = np.asarray(r_rec, dtype=float)
    return float(omega * (rs[0] * rr[1] - rs[1] * rr[0]) / C ** 2)


def max_sagnac(a: float = GPS_SEMI_MAJOR_AXIS, R: float = R_EARTH, omega: float = OMEGA_EARTH) -> float:
    """Largest Sagnac term: satellite on the equatorial horizon of an equatorial receiver, omega a R sin(theta)/c^2
    with cos(theta) = R/a."""
    theta = np.arccos(R / a)
    return float(omega * a * R * np.sin(theta) / C ** 2)


def position_error_without_correction(t: np.ndarray | float, offset: float | None = None) -> np.ndarray | float:
    """Pseudorange error c * offset * t [m] that accumulates if the clock rate is not corrected
    (11.6 km after one day for the GPS budget)."""
    if offset is None:
        offset = gps_clock_budget().net
    return C * offset * np.asarray(t, dtype=float)


def pass_doppler(f0: float, v_radial: np.ndarray, v_total: float) -> np.ndarray:
    """Received frequency during a pass: first-order Doppler from the radial velocity (positive = approaching)
    plus the second-order (time-dilation) shift -f0 v_total^2/(2c^2), which is the same at every point of the pass."""
    vr = np.asarray(v_radial, dtype=float)
    return f0 * (1.0 + vr / C) * (1.0 - gamma_minus_one(v_total) / gamma(v_total))


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the GPS clock budget and write gps_relativity.png (offset against altitude; accumulated error;
    eccentricity and Sagnac terms)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    b = gps_clock_budget()
    v = circular_orbit_speed(GPS_SEMI_MAJOR_AXIS)
    print(f"GPS orbit: a = {GPS_SEMI_MAJOR_AXIS/1e3:.2f} km, v = {v:.1f} m/s, gamma - 1 = {gamma_minus_one(v):.3e}")
    u = b.us_per_day()
    for k in ("velocity", "gravity", "ground_rotation", "ground_j2", "net"):
        print(f"  {k:16s} {getattr(b, k) if k != 'net' else b.net:+.4e}   {u[k]:+8.3f} us/day")
    print(f"oracle -3GM/(2ac^2) - Phi0/c^2 = {gps_frequency_offset_oracle():+.4e}")
    print(f"factory clock frequency: {factory_clock_frequency():.5f} Hz  (nominal {GPS_CLOCK_HZ:.0f} Hz)")
    print(f"crossover radius (simple model): {crossover_radius()/1e3:.1f} km = {crossover_radius()/R_EARTH:.4f} R_E, "
          f"altitude {(crossover_radius()-R_EARTH)/1e3:.0f} km")
    for name, a in (("ISS (h = 400 km)", R_EARTH + 400e3), ("GPS", GPS_SEMI_MAJOR_AXIS), ("GEO", 42_164e3)):
        bb = gps_clock_budget(a).us_per_day()
        print(f"  {name:18s} velocity {bb['velocity']:+7.2f}  gravity {bb['gravity']:+7.2f}  net {bb['net']:+7.2f} us/day")
    print(f"eccentricity correction amplitude at e = 0.02: {abs(eccentricity_correction(0.02, np.pi/2))*1e9:.1f} ns")
    print(f"maximum Sagnac term: {max_sagnac()*1e9:.1f} ns")
    print(f"pseudorange error after one day without correction: {position_error_without_correction(SECONDS_PER_DAY)/1e3:.2f} km")
    df1 = GPS_L1_HZ * 800.0 / C
    df2 = GPS_L1_HZ * gamma_minus_one(v)
    print(f"L1 Doppler at 800 m/s radial: {df1:+.0f} Hz; second-order (time-dilation) shift: {-df2:+.3f} Hz")

    alts = np.linspace(200e3, 40_000e3, 400)
    terms = np.array([[*gps_clock_budget(R_EARTH + h).us_per_day().values()] for h in alts])
    t_days = np.linspace(0.0, 1.0, 200)
    E = np.linspace(0.0, 2 * np.pi, 361)
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    ax[0, 0].plot(alts / 1e3, terms[:, 0], label="velocity")
    ax[0, 0].plot(alts / 1e3, terms[:, 1], label="gravity")
    ax[0, 0].plot(alts / 1e3, terms[:, 4], "k", label="net")
    ax[0, 0].axhline(0, color="gray", lw=0.5)
    ax[0, 0].axvline((GPS_SEMI_MAJOR_AXIS - R_EARTH) / 1e3, color="gray", ls="--", lw=0.8)
    ax[0, 0].set_xlabel("altitude (km)"); ax[0, 0].set_ylabel("clock offset (us/day)"); ax[0, 0].legend()
    ax[0, 0].set_title("satellite clock relative to the geoid")
    ax[0, 1].plot(t_days * 24, position_error_without_correction(t_days * SECONDS_PER_DAY) / 1e3)
    ax[0, 1].set_xlabel("time since last correction (h)"); ax[0, 1].set_ylabel("pseudorange error (km)")
    ax[0, 1].set_title("error without the relativistic correction")
    ax[1, 0].plot(np.degrees(E), eccentricity_correction(0.02, E) * 1e9)
    ax[1, 0].set_xlabel("eccentric anomaly (deg)"); ax[1, 0].set_ylabel("eccentricity correction (ns)")
    ax[1, 0].set_title("periodic term, e = 0.02")
    az = np.linspace(-np.pi, np.pi, 361)
    theta = np.arccos(R_EARTH / GPS_SEMI_MAJOR_AXIS)
    sag = [sagnac_correction(GPS_SEMI_MAJOR_AXIS * np.array([np.cos(theta) * 1.0, np.sin(theta) * np.sin(p), 0.0]),
                             np.array([R_EARTH, 0.0, 0.0])) for p in az]
    ax[1, 1].plot(np.degrees(az), np.array(sag) * 1e9)
    ax[1, 1].set_xlabel("satellite east-west position angle (deg)"); ax[1, 1].set_ylabel("Sagnac term (ns)")
    ax[1, 1].set_title("Sagnac term for an equatorial receiver")
    fig.tight_layout()
    fig.savefig(out_dir / "gps_relativity.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
