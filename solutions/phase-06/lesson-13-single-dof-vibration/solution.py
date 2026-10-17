"""P06.13 Single-DOF Vibration: reference solution.

A camera on a cantilever bracket, as a single-degree-of-freedom vibrating system:

* the model: stiffness of a cantilever bracket 3EI/L^3, an effective mass that includes 33/140 of the
  bracket's own mass (Rayleigh), natural frequency and damping ratio;
* free response in all three damping regimes in closed form; damping identified from a ring-down
  ("tap test") record by the logarithmic decrement over n cycles;
* frequency response functions: receptance, mobility and accelerance, and damping identified from the
  half-power bandwidth of a sampled FRF;
* base excitation, as on a shaker table or a vibrating spacecraft panel: absolute transmissibility
  (acceleration and force), relative-motion transmissibility (bracket deflection and stress), and the
  full transient response to a base sine that starts at t = 0, in closed form;
* Coulomb (dry-friction) damping: the closed-form free decay with its linear envelope, and the
  equivalent viscous damping ratio;
* every closed form checked against `scipy.integrate.solve_ivp`.

Numbers are illustrative (a typical 6061-T6 bracket and a 0.6 kg camera); the tap-test record is
SYNTHETIC, generated from a seed, and never a measurement.

Conventions: x absolute displacement of the mass, y base displacement, z = x - y relative displacement;
SI units; omega in rad/s; zeta dimensionless.

Run:  python solution.py [out_dir]     (prints the tables, writes sdof_vibration.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-06/lesson-13-single-dof-vibration
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

G0 = 9.80665
E_6061 = 68.9e9          # Pa, typical Young's modulus of 6061-T6 aluminium
RHO_6061 = 2700.0        # kg/m^3, typical density
SY_6061 = 276e6          # Pa, typical tensile yield strength of 6061-T6 (not a design allowable)


def camera_bracket() -> dict[str, float]:
    """Illustrative case: a 0.6 kg camera at the tip of a 6061-T6 cantilever bracket 100 mm long,
    40 mm wide and 6 mm thick (bending about the thin direction), and a 2 g (0-peak) base sine."""
    return dict(m_tip=0.6, L=0.100, b=0.040, t=0.006, E=E_6061, rho=RHO_6061, zeta=0.02, base_g=2.0)


# ---------------------------------------------------------------- the model

# LEARNER
def cantilever_stiffness(E: float, I: float, L: float) -> float:
    """Tip stiffness of a uniform cantilever, k = 3 E I / L^3 (N/m), from w_tip = P L^3 / (3 E I) (P05.12)."""
    return 3.0 * E * I / L ** 3


# LEARNER
def effective_mass(m_tip: float, m_beam: float) -> float:
    """Tip mass plus the Rayleigh share of a uniform cantilever's own mass, m_tip + (33/140) m_beam (kg),
    from the kinetic energy of the static deflection shape moving at the tip velocity."""
    return m_tip + 33.0 / 140.0 * m_beam


# LEARNER
def natural_frequency(k: float, m: float) -> float:
    """Undamped natural frequency w_n = sqrt(k/m) (rad/s)."""
    return float(np.sqrt(k / m))


# LEARNER
def damping_ratio(c: float, m: float, k: float) -> float:
    """zeta = c / (2 sqrt(k m)): the damping over the critical damping."""
    return c / (2.0 * np.sqrt(k * m))


# LEARNER
def free_response(t: np.ndarray, x0: float, v0: float, wn: float, zeta: float) -> np.ndarray:
    """Closed-form free response of x'' + 2 zeta wn x' + wn^2 x = 0 with x(0) = x0, x'(0) = v0, in all
    three regimes: underdamped (zeta < 1), critically damped (zeta = 1, to within 1e-9) and overdamped."""
    t = np.asarray(t, dtype=float)
    if abs(zeta - 1.0) < 1e-9:
        return (x0 + (v0 + wn * x0) * t) * np.exp(-wn * t)
    if zeta < 1.0:
        wd = wn * np.sqrt(1.0 - zeta ** 2)
        s = zeta * wn
        return np.exp(-s * t) * (x0 * np.cos(wd * t) + (v0 + s * x0) / wd * np.sin(wd * t))
    r = wn * np.sqrt(zeta ** 2 - 1.0)
    s1, s2 = -zeta * wn + r, -zeta * wn - r
    c1 = (v0 - s2 * x0) / (s1 - s2)
    c2 = x0 - c1
    return c1 * np.exp(s1 * t) + c2 * np.exp(s2 * t)


def find_peaks_refined(t: np.ndarray, x: np.ndarray, f_guess: float) -> tuple[np.ndarray, np.ndarray]:
    """Positive peaks of a sampled oscillation, at most one per 0.7 / f_guess seconds (so that noise
    ripples are not counted as cycles), each refined by a parabola through three samples (given)."""
    from scipy.signal import find_peaks
    dt = t[1] - t[0]
    i, _ = find_peaks(x, height=0.0, distance=max(1, int(0.7 / (f_guess * dt))))
    i = i[(i > 0) & (i < len(x) - 1)]
    a, b, c = x[i - 1], x[i], x[i + 1]
    den = a - 2 * b + c
    d = np.where(den != 0, 0.5 * (a - c) / np.where(den != 0, den, 1.0), 0.0)
    return t[i] + d * dt, b - 0.25 * (a - c) * d


# LEARNER
def log_decrement(peaks: np.ndarray, n: int) -> float:
    """Logarithmic decrement from successive positive peaks, averaged over n cycles:
    delta = ln(x_0 / x_n) / n."""
    peaks = np.asarray(peaks, dtype=float)
    return float(np.log(peaks[0] / peaks[n]) / n)


# LEARNER
def zeta_from_log_decrement(delta: float) -> float:
    """Exact inversion of delta = 2 pi zeta / sqrt(1 - zeta^2): zeta = delta / sqrt(4 pi^2 + delta^2)."""
    return delta / np.sqrt(4.0 * np.pi ** 2 + delta ** 2)


def synthetic_tap_test(seed: int, fn: float = 78.0, zeta: float = 0.02, x0: float = 1e-4,
                       duration: float = 0.5, fs: float = 5000.0, noise: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """SYNTHETIC ring-down record of a tapped bracket (not a measurement): a free decay with the given
    natural frequency (Hz) and damping, sampled at fs (Hz), plus seeded Gaussian noise of standard
    deviation ``noise`` (m). Returns (t, x) (given)."""
    rng = np.random.default_rng(seed)
    t = np.arange(0.0, duration, 1.0 / fs)
    x = free_response(t, x0, 0.0, 2 * np.pi * fn, zeta)
    return t, x + rng.normal(0.0, noise, t.size)


# ---------------------------------------------------------------- frequency response

# LEARNER
def frf(omega: np.ndarray, m: float, c: float, k: float, kind: str = "receptance") -> np.ndarray:
    """Complex frequency response of m x'' + c x' + k x = F e^{i w t}: receptance X/F = 1/(k - m w^2 + i c w)
    (m/N), mobility V/F = i w X/F (m/s/N), accelerance A/F = -w^2 X/F (m/s^2/N)."""
    w = np.asarray(omega, dtype=float)
    H = 1.0 / (k - m * w ** 2 + 1j * c * w)
    if kind == "receptance":
        return H
    if kind == "mobility":
        return 1j * w * H
    if kind == "accelerance":
        return -(w ** 2) * H
    raise ValueError(f"unknown FRF kind {kind!r}")


# LEARNER
def half_power_zeta(omega: np.ndarray, H: np.ndarray) -> tuple[float, float]:
    """Natural frequency and damping ratio from a sampled receptance magnitude by the half-power method:
    the peak |H|max, the two frequencies w1 < w2 where |H| = |H|max / sqrt(2) (linear interpolation),
    w_n ~ (w1 + w2)/2 and zeta ~ (w2 - w1) / (2 w_n). Accurate to O(zeta^2)."""
    w = np.asarray(omega, dtype=float)
    a = np.abs(H)
    k = int(np.argmax(a))
    level = a[k] / np.sqrt(2.0)
    i = k
    while a[i] > level:
        i -= 1
    w1 = w[i] + (level - a[i]) * (w[i + 1] - w[i]) / (a[i + 1] - a[i])
    j = k
    while a[j] > level:
        j += 1
    w2 = w[j - 1] + (level - a[j - 1]) * (w[j] - w[j - 1]) / (a[j] - a[j - 1])
    wn = 0.5 * (w1 + w2)
    return wn, (w2 - w1) / (2.0 * wn)


# ---------------------------------------------------------------- base excitation

# LEARNER
def transmissibility(r: np.ndarray, zeta: float) -> np.ndarray:
    """Absolute transmissibility |X/Y| = |A_mass/A_base| = |F_base/(k Y)| for base excitation:
    sqrt((1 + (2 zeta r)^2) / ((1 - r^2)^2 + (2 zeta r)^2)), r = w / w_n."""
    r = np.asarray(r, dtype=float)
    return np.sqrt((1 + (2 * zeta * r) ** 2) / ((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2))


# LEARNER
def relative_transmissibility(r: np.ndarray, zeta: float) -> np.ndarray:
    """Relative-motion transmissibility |Z/Y| = r^2 / sqrt((1 - r^2)^2 + (2 zeta r)^2): the bracket's
    deflection per unit base displacement. Per unit base acceleration, |Z| = |A_base| / w_n^2 times
    1/sqrt((1 - r^2)^2 + (2 zeta r)^2)."""
    r = np.asarray(r, dtype=float)
    return r ** 2 / np.sqrt((1 - r ** 2) ** 2 + (2 * zeta * r) ** 2)


# LEARNER
def base_sine_response(t: np.ndarray, A: float, omega: float, wn: float, zeta: float) -> tuple[np.ndarray, np.ndarray]:
    """Closed-form response to a base acceleration y'' = A sin(w t) switched on at t = 0, with the mass
    initially moving with the base (z(0) = z'(0) = 0, zeta < 1). The relative motion obeys
    z'' + 2 zeta wn z' + wn^2 z = -A sin(w t): the steady part is Im(Z e^{i w t}) with
    Z = -A / (wn^2 - w^2 + 2 i zeta wn w), plus the free response that cancels its initial values.
    Returns (z, x_ddot) with the absolute acceleration x'' = -(2 zeta wn z' + wn^2 z)."""
    t = np.asarray(t, dtype=float)
    Z = -A / (wn ** 2 - omega ** 2 + 2j * zeta * wn * omega)
    zs = np.imag(Z * np.exp(1j * omega * t))
    zsd = np.imag(1j * omega * Z * np.exp(1j * omega * t))
    z0s, v0s = np.imag(Z), np.real(omega * Z)               # steady part's values at t = 0
    zh = free_response(t, -z0s, -v0s, wn, zeta)
    # derivative of the homogeneous part, by the same closed form (underdamped)
    wd = wn * np.sqrt(1 - zeta ** 2)
    s = zeta * wn
    x0, v0 = -z0s, -v0s
    zhd = np.exp(-s * t) * (v0 * np.cos(wd * t) - (s * v0 + wn ** 2 * x0) / wd * np.sin(wd * t))
    z = zs + zh
    zd = zsd + zhd
    return z, -(2 * zeta * wn * zd + wn ** 2 * z)


# LEARNER
def simulate_base_sine(t: np.ndarray, A: float, omega: float, wn: float, zeta: float,
                       rtol: float = 1e-10, atol: float = 1e-14) -> tuple[np.ndarray, np.ndarray]:
    """The same problem integrated with solve_ivp (DOP853, max_step a tenth of the shorter period).
    Returns (z, x_ddot) at t."""
    f = lambda tt, u: [u[1], -A * np.sin(omega * tt) - 2 * zeta * wn * u[1] - wn ** 2 * u[0]]
    hmax = 0.1 * 2 * np.pi / max(wn, omega)
    sol = solve_ivp(f, (t[0], t[-1]), [0.0, 0.0], method="DOP853", t_eval=t, rtol=rtol, atol=atol, max_step=hmax)
    z, zd = sol.y
    return z, -(2 * zeta * wn * zd + wn ** 2 * z)


# LEARNER
def bracket_root_stress(z: np.ndarray, E: float, t: float, L: float) -> np.ndarray:
    """Root bending stress of a cantilever of thickness t whose tip deflects z under a tip load:
    M_root = (3 E I / L^3) z L, sigma = M (t/2) / I = 3 E t z / (2 L^2)."""
    return 3.0 * E * t * np.asarray(z) / (2.0 * L ** 2)


# ---------------------------------------------------------------- Coulomb damping

# LEARNER
def coulomb_free_response(t: np.ndarray, x0: float, wn: float, d: float) -> np.ndarray:
    """Closed-form free response of x'' + wn^2 x = -wn^2 d sign(x'), released from rest at x0 > 0, with
    d = F/k the friction force over the stiffness. Each half cycle is a cosine of the undamped frequency
    about the offset +d or -d; the turning-point amplitude falls by 2d per half cycle (4d per cycle),
    linearly, and the mass sticks at the first turning point with |x| <= d."""
    t = np.asarray(t, dtype=float)
    out = np.empty_like(t)
    half = np.pi / wn
    xs, t_start = float(x0), 0.0               # position at the turning point that starts the half cycle
    while True:
        if abs(xs) <= d:                       # the spring cannot overcome static friction: it sticks
            out[t >= t_start] = xs
            break
        centre = np.sign(xs) * d               # moving towards -sign(xs): friction pushes along +sign(xs)
        sel = (t >= t_start) & (t < t_start + half)
        out[sel] = centre + (xs - centre) * np.cos(wn * (t[sel] - t_start))
        xs = 2.0 * centre - xs                 # the next turning point: |x| falls by 2d
        t_start += half
        if t_start > t[-1]:
            break
    return out


# LEARNER
def simulate_coulomb(t: np.ndarray, x0: float, wn: float, d: float) -> np.ndarray:
    """Coulomb free decay integrated with solve_ivp, one half cycle at a time: integrate with the friction
    sign fixed until the velocity crosses zero (an event), then flip the sign; stop when the spring force
    at a turning point no longer exceeds the friction (|x| <= d)."""
    t = np.asarray(t, dtype=float)
    out = np.full_like(t, np.nan)
    y, t0, s = np.array([x0, 0.0]), t[0], -1.0         # s = sign of the velocity in this half cycle
    while t0 < t[-1]:
        if abs(y[0]) <= d:
            out[t >= t0] = y[0]
            break
        f = lambda tt, u, s=s: [u[1], -wn ** 2 * u[0] - wn ** 2 * d * s]
        ev = lambda tt, u, s=s: u[1]
        ev.terminal, ev.direction = True, -s
        sel = t >= t0
        sol = solve_ivp(f, (t0, t[-1] + 1.0), y, method="DOP853", rtol=1e-11, atol=1e-14, events=ev, dense_output=True)
        t1 = sol.t_events[0][0] if sol.t_events[0].size else t[-1] + 1.0
        seg = sel & (t < t1)
        out[seg] = sol.sol(t[seg])[0]
        y = sol.sol(t1)
        y[1] = 0.0
        t0, s = t1, -s
    return out


# LEARNER
def coulomb_equivalent_zeta(d: float, X: float) -> float:
    """Equivalent viscous damping ratio of Coulomb friction at amplitude X near resonance, equating the
    energy lost per cycle 4 F X to that of a viscous damper pi c w X^2 at w = w_n:
    c_eq = 4 F / (pi w_n X), zeta_eq = 2 d / (pi X) with d = F/k."""
    return 2.0 * d / (np.pi * X)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save sdof_vibration.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    c = camera_bracket()
    I = c["b"] * c["t"] ** 3 / 12
    k = cantilever_stiffness(c["E"], I, c["L"])
    m_beam = c["rho"] * c["b"] * c["t"] * c["L"]
    m = effective_mass(c["m_tip"], m_beam)
    wn = natural_frequency(k, m)
    print(f"bracket: I = {I:.4e} m^4, k = {k:.4e} N/m, beam mass {m_beam * 1e3:.1f} g, effective mass {m:.4f} kg")
    print(f"natural frequency {wn:.2f} rad/s = {wn / 2 / np.pi:.2f} Hz (tip mass only: {np.sqrt(k / c['m_tip']) / 2 / np.pi:.2f} Hz)")

    # tap test (SYNTHETIC)
    tt, xx = synthetic_tap_test(seed=13, fn=wn / 2 / np.pi, zeta=c["zeta"])
    tp, xp = find_peaks_refined(tt, xx, wn / 2 / np.pi)
    print("\nSYNTHETIC tap test (seed 13, zeta 0.02 used to generate it):")
    for n in (1, 5, 10, 20):
        z_est = zeta_from_log_decrement(log_decrement(xp, n))
        print(f"   log decrement over {n:2d} cycles: zeta = {z_est:.5f}")
    fd = 1.0 / np.mean(np.diff(tp[:21]))
    print(f"   damped frequency from 20 peak intervals: {fd:.3f} Hz")

    # FRF and half power
    zeta = c["zeta"]
    cc = 2 * zeta * np.sqrt(k * m)
    w = np.linspace(0.5 * wn, 1.5 * wn, 4001)
    wn_hp, z_hp = half_power_zeta(w, frf(w, m, cc, k))
    print(f"\nhalf-power method on the receptance (4001 points, 0.5-1.5 w_n): w_n {wn_hp:.3f} rad/s, zeta {z_hp:.5f}")

    # sine test
    A = c["base_g"] * G0
    Tpk = float(np.max(transmissibility(np.linspace(0.5, 1.5, 200001), zeta)))
    zpk = A / wn ** 2 / (2 * zeta * np.sqrt(1 - zeta ** 2))
    sig = float(bracket_root_stress(zpk, c["E"], c["t"], c["L"]))
    print(f"\nsine test {c['base_g']:.1f} g (0-peak) base: peak transmissibility {Tpk:.3f} -> {Tpk * c['base_g']:.1f} g at the camera; "
          f"peak relative deflection {zpk * 1e3:.3f} mm; root stress {sig / 1e6:.1f} MPa "
          f"(typical 6061-T6 yield {SY_6061 / 1e6:.0f} MPa)")
    t = np.linspace(0.0, 0.5, 20001)
    zc, ac = base_sine_response(t, A, wn, wn, zeta)
    zn, an = simulate_base_sine(t, A, wn, wn, zeta)
    print(f"   resonant dwell from rest, 0.5 s: closed form vs solve_ivp max |dz| {np.max(np.abs(zc - zn)):.1e} m "
          f"(peak {np.max(np.abs(zc)) * 1e3:.3f} mm); max |d accel| {np.max(np.abs(ac - an)):.1e} m/s^2")
    t1 = np.linspace(0.0, 0.2, 8001)
    z1, a1 = base_sine_response(t1, A, 2 * np.pi * 20.0, wn, zeta)
    z2, a2 = simulate_base_sine(t1, A, 2 * np.pi * 20.0, wn, zeta)
    print(f"   20 Hz from rest: max |dz| {np.max(np.abs(z1 - z2)):.1e} m")

    # Coulomb
    d = 2e-6
    tc = np.linspace(0, 0.2, 4001)
    xc = coulomb_free_response(tc, 1e-4, wn, d)
    xs = simulate_coulomb(tc, 1e-4, wn, d)
    print(f"\nCoulomb decay (d = F/k = {d * 1e6:.1f} um, x0 = 100 um): closed form vs solve_ivp max diff {np.nanmax(np.abs(xc - xs)):.1e} m; "
          f"equivalent zeta at 100 um {coulomb_equivalent_zeta(d, 1e-4):.4f}, at 20 um {coulomb_equivalent_zeta(d, 2e-5):.4f}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].plot(tt * 1e3, xx * 1e6, lw=0.6)
    ax[0].plot(tp * 1e3, xp * 1e6, "r.", ms=3)
    ax[0].set_xlabel("t (ms)")
    ax[0].set_ylabel("tip deflection (um)")
    ax[0].set_title("SYNTHETIC tap test (seed 13)")
    r = np.logspace(-1, 1, 600)
    for zz in (0.01, 0.02, 0.05, 0.2):
        ax[1].loglog(r, transmissibility(r, zz), label=f"absolute, zeta = {zz}")
    ax[1].loglog(r, relative_transmissibility(r, 0.02), "k--", label="relative, zeta = 0.02")
    ax[1].axvline(np.sqrt(2), color="0.6", lw=0.8)
    ax[1].set_xlabel("r = w / w_n")
    ax[1].set_ylabel("transmissibility")
    ax[1].set_title("base excitation")
    ax[1].legend(fontsize=7)
    ax[2].plot(tc * 1e3, xc * 1e6, label="Coulomb, closed form")
    ax[2].plot(tc[::40] * 1e3, xs[::40] * 1e6, "o", ms=3, label="Coulomb, solve_ivp")
    ax[2].plot(tc * 1e3, free_response(tc, 1e-4, 0, wn, coulomb_equivalent_zeta(d, 1e-4)) * 1e6, lw=0.7,
               label="viscous, same initial loss")
    ax[2].set_xlabel("t (ms)")
    ax[2].set_ylabel("x (um)")
    ax[2].set_title("linear vs exponential decay")
    ax[2].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "sdof_vibration.png", dpi=110)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'sdof_vibration.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
