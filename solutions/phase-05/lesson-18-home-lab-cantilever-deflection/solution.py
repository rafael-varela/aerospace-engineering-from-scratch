"""P05.18 Home Lab: Cantilever Beam Deflection: reference solution.

A strip clamped to a table edge is a cantilever. Hang a known mass at distance a from the clamp and
Euler-Bernoulli beam theory predicts the deflection everywhere along it:

    w(x) = P x^2 (3a - x) / (6 E I)       for 0 <= x <= a
    w(x) = P a^2 (3x - a) / (6 E I)       for a <= x <= L     (straight beyond the load)

This module turns a measurement of that into evidence:

* the prediction with its uncertainty, from E (typical handbook value with a stated spread) and from
  the measured width, thickness and load position, by first-order propagation and checked by a seeded
  Monte Carlo;
* the safe load from the bending stress at the clamp, so the strip stays elastic;
* a weighted straight-line fit of deflection against load: the slope gives E with its uncertainty, the
  intercept shows a zero offset, chi^2 shows whether the scatter matches the reading uncertainty;
* a fit of the deflected shape along the beam to the cubic plus a rigid rotation at the clamp, which
  detects a compliant clamp;
* two independent checks of the theory itself: a finite-difference solution of E I w'' = M, and the
  large-deflection elastica (solve_bvp), which says how far the linear theory can be trusted;
* a SYNTHETIC dataset generator for learners without a kit. Synthetic data are not measurements.

Units: SI in the functions (N, m, Pa); the driver prints millimetres and grams.

Run:  python main.py [out_dir]     (prints prediction, fit and checks, writes cantilever_lab.png)
Test: pytest curriculum/phase-05/lesson-18-home-lab-cantilever-deflection
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_bvp

G = 9.81                     # m/s^2; use your local g if you know it (P00.14)
E_AL_TYPICAL = 69.0e9        # Pa, typical for 6061/6063 aluminium; +/- 2 GPa covers alloy and handbook spread
U_E_AL = 2.0e9
FTY_6063_T5 = 145e6          # Pa, typical yield of the weakest common hardware-store alloy (6063-T5)
RHO_AL = 2700.0              # kg/m^3


# ---------------------------------------------------------------- section and theory

def rect_second_moment(b: float, h: float, u_b: float = 0.0, u_h: float = 0.0) -> tuple[float, float]:
    """I = b h^3 / 12 of a rectangular strip bent about its weak axis, and its standard uncertainty
    u_I = I sqrt((u_b/b)^2 + (3 u_h/h)^2): the thickness enters cubed."""
    I = b * h ** 3 / 12.0
    return I, I * np.hypot(u_b / b, 3.0 * u_h / h)


def deflection_point_load(x: np.ndarray, P: float, a: float, EI: float) -> np.ndarray:
    """Downward-positive deflection of a cantilever (clamp at x = 0) under a downward point load P at x = a: P x^2 (3a - x)/(6EI) for x <= a, P a^2 (3x - a)/(6EI) beyond."""
    x = np.asarray(x, dtype=float)
    return np.where(x <= a, P * x ** 2 * (3 * a - x), P * a ** 2 * (3 * x - a)) / (6.0 * EI)


def deflection_self_weight(x: np.ndarray, q: float, L: float, EI: float) -> np.ndarray:
    """Downward deflection under the strip's own weight q per unit length: q x^2 (6L^2 - 4Lx + x^2)/(24EI).
    It is present before any load is hung, so a measurement of the CHANGE in reading removes it."""
    x = np.asarray(x, dtype=float)
    return q * x ** 2 * (6 * L ** 2 - 4 * L * x + x ** 2) / (24.0 * EI)


def compliance_shape(x: float, a: float) -> float:
    """g(x, a) such that w(x) = P g / (E I): x^2(3a - x)/6 for x <= a, a^2(3x - a)/6 beyond."""
    return float(deflection_point_load(x, 1.0, a, 1.0))


def clamp_stress(P: float, a: float, b: float, h: float) -> float:
    """Largest bending stress, at the clamp's edge: sigma = M c / I = P a (h/2) / (b h^3/12) = 6 P a / (b h^2)."""
    return 6.0 * P * a / (b * h ** 2)


def max_safe_load(F_allow: float, a: float, b: float, h: float) -> float:
    """Load P [N] at which the clamp stress reaches F_allow. The lesson uses F_allow = F_ty / 2."""
    return F_allow * b * h ** 2 / (6.0 * a)


# ---------------------------------------------------------------- prediction with uncertainty

def predict_deflection(P: float, a: float, E: float, u_E: float, b: float, u_b: float, h: float, u_h: float,
                       u_a: float, x: float | None = None, u_P: float = 0.0) -> tuple[float, float]:
    """Predicted deflection and its standard uncertainty by first-order (GUM) propagation, with numerical
    partial derivatives (central differences, 1e-6 relative) and independent inputs.

    x is the reading position; None means the reading is taken at the load point, so x moves with a
    (w = P a^3 / 3EI, and a enters cubed). A fixed x (a mark on the strip) gives w = P g(x, a) / EI."""
    def f(v):
        P_, a_, E_, b_, h_ = v
        xx = a_ if x is None else x
        return float(deflection_point_load(xx, P_, a_, E_ * b_ * h_ ** 3 / 12.0))

    v = np.array([P, a, E, b, h], dtype=float)
    u = np.array([u_P, u_a, u_E, u_b, u_h], dtype=float)
    d0 = f(v)
    var = 0.0
    for i in range(len(v)):
        dv = np.zeros_like(v)
        dv[i] = 1e-6 * abs(v[i])
        c = (f(v + dv) - f(v - dv)) / (2 * dv[i])
        var += (c * u[i]) ** 2
    return d0, float(np.sqrt(var))


def predict_deflection_mc(P: float, a: float, E: float, u_E: float, b: float, u_b: float, h: float, u_h: float,
                          u_a: float, x: float | None = None, n: int = 200_000, seed: int = 18) -> tuple[float, float]:
    """Independent check of predict_deflection: Monte Carlo with normal inputs (seeded with
    numpy.random.default_rng(seed)). Returns the sample mean and standard deviation of the deflection."""
    rng = np.random.default_rng(seed)
    Es = rng.normal(E, u_E, n)
    bs = rng.normal(b, u_b, n)
    hs = rng.normal(h, u_h, n)
    as_ = rng.normal(a, u_a, n)
    xs = as_ if x is None else np.full(n, float(x))
    EI = Es * bs * hs ** 3 / 12.0
    d = np.where(xs <= as_, P * xs ** 2 * (3 * as_ - xs), P * as_ ** 2 * (3 * xs - as_)) / (6.0 * EI)
    return float(np.mean(d)), float(np.std(d))


# ---------------------------------------------------------------- fitting the measurement

def fit_line(P: np.ndarray, d: np.ndarray, u_d: float | np.ndarray) -> dict[str, float]:
    """Weighted least-squares line d = s P + c with standard uncertainties of s and c from the
    covariance (u_d taken as known), chi^2 and the degrees of freedom."""
    P = np.asarray(P, dtype=float)
    d = np.asarray(d, dtype=float)
    w = 1.0 / np.broadcast_to(np.asarray(u_d, dtype=float), P.shape) ** 2
    A = np.column_stack([P, np.ones_like(P)])
    cov = np.linalg.inv(A.T @ (w[:, None] * A))
    s, c = cov @ (A.T @ (w * d))
    r = d - (s * P + c)
    return {"slope": float(s), "u_slope": float(np.sqrt(cov[0, 0])), "intercept": float(c),
            "u_intercept": float(np.sqrt(cov[1, 1])), "chi2": float(np.sum(w * r ** 2)), "dof": len(P) - 2,
            "residuals": r}


def modulus_from_slope(slope: float, u_slope: float, a: float, u_a: float, b: float, u_b: float,
                       h: float, u_h: float, x: float | None = None) -> tuple[float, float]:
    """E from the measured compliance slope = d(deflection)/dP: E = g / (slope I), with g = a^3/3 when the
    reading is taken at the load point (x None) or g(x, a) at a fixed mark x. The relative uncertainty
    combines slope, a, b and h in quadrature (first order); with x None, a contributes 3 u_a / a."""
    I = b * h ** 3 / 12.0
    g = (lambda aa: aa ** 3 / 3.0) if x is None else (lambda aa: compliance_shape(x, aa))
    E = g(a) / (slope * I)
    da = 1e-6 * a
    dlng_da = (np.log(g(a + da)) - np.log(g(a - da))) / (2 * da)
    rel = np.sqrt((u_slope / slope) ** 2 + (dlng_da * u_a) ** 2 + (u_b / b) ** 2 + (3 * u_h / h) ** 2)
    return float(E), float(E * rel)


def fit_profile(x: np.ndarray, w: np.ndarray, P: float, a: float, u_w: float) -> dict[str, float]:
    """Fit the measured shape w(x) under load P at a to  w = P g(x, a) / EI + theta0 x : the cubic of beam
    theory plus a rigid rotation theta0 of a compliant clamp. Returns EI, theta0, their uncertainties
    and chi^2. A clamp that is truly rigid gives theta0 consistent with zero."""
    x = np.asarray(x, dtype=float)
    g = np.array([compliance_shape(xi, a) for xi in x])
    A = np.column_stack([P * g, x])
    cov = np.linalg.inv(A.T @ A) * u_w ** 2
    (inv_EI, th), *_ = np.linalg.lstsq(A, np.asarray(w, dtype=float), rcond=None)
    r = w - A @ np.array([inv_EI, th])
    EI = 1.0 / inv_EI
    return {"EI": float(EI), "u_EI": float(EI ** 2 * np.sqrt(cov[0, 0])), "theta0": float(th),
            "u_theta0": float(np.sqrt(cov[1, 1])), "chi2": float(np.sum(r ** 2) / u_w ** 2), "dof": len(x) - 2}


def en_number(x1: float, u1: float, x2: float, u2: float, k: float = 2.0) -> float:
    """E_n = |x1 - x2| / sqrt(U1^2 + U2^2) with expanded uncertainties U = k u. E_n < 1: consistent."""
    return abs(x1 - x2) / (k * np.hypot(u1, u2))


# ---------------------------------------------------------------- independent checks of the theory

def fd_cantilever(P: float, a: float, L: float, EI: float, n: int = 200) -> tuple[np.ndarray, np.ndarray]:
    """Finite differences on E I w'' = M(x), M(x) = P (a - x) for x < a and 0 beyond (downward w positive),
    with w(0) = 0 and w'(0) = 0 imposed through a ghost node (w_{-1} = w_1). Second order in h = L/n.
    Returns (x, w)."""
    x = np.linspace(0.0, L, n + 1)
    h = L / n
    rhs = np.where(x < a, P * (a - x), 0.0) / EI * h ** 2
    A = np.zeros((n + 1, n + 1))
    A[0, 0] = 1.0                                   # w_0 = 0
    b = rhs.copy()
    b[0] = 0.0
    # at node 0 the curvature equation with the ghost w_{-1} = w_1 (zero slope) reads 2 w_1 - 2 w_0 = h^2 M_0/EI;
    # use it as the row for node 1's unknown, and the ordinary three-point rows for nodes 1..n-1
    A[1, 0], A[1, 1] = -2.0, 2.0
    b[1] = rhs[0]
    for i in range(1, n):
        A[i + 1, i - 1], A[i + 1, i], A[i + 1, i + 1] = 1.0, -2.0, 1.0
        b[i + 1] = rhs[i]
    return x, np.linalg.solve(A, b)


def elastica_deflection(P: float, a: float, EI: float, n: int = 400) -> float:
    """Large-deflection (elastica) vertical deflection at the load point of a cantilever, by solve_bvp on
    the arc length s in [0, a]: E I theta'' = -P cos(theta), theta(0) = 0, theta'(a) = 0, then
    w(a) = integral sin(theta) ds. The load stays vertical; the strip does not stretch. For small loads
    it tends to P a^3 / (3 E I); for large ones it is smaller, because the load's lever arm shrinks."""
    s = np.linspace(0.0, a, n)
    lin = P / EI * (a * s - s ** 2 / 2)                     # small-deflection slope as the first guess
    y0 = np.vstack([lin, P / EI * (a - s)])

    def rhs(_, y):
        return np.vstack([y[1], -P / EI * np.cos(y[0])])

    def bc(ya, yb):
        return np.array([ya[0], yb[1]])

    sol = solve_bvp(rhs, bc, s, y0, tol=1e-10, max_nodes=100_000)
    if not sol.success:
        raise RuntimeError(sol.message)
    ss = np.linspace(0.0, a, 4001)
    th = sol.sol(ss)[0]
    return float(cumulative_trapezoid(np.sin(th), ss, initial=0.0)[-1])


# ---------------------------------------------------------------- SYNTHETIC data

def make_synthetic_cantilever_data(seed: int = 518) -> dict:
    """SYNTHETIC measurements of the lesson's strip (25 x 3 mm aluminium, 300 mm free length, load at
    290 mm, deflection read at the load point with a ruler photographed by a phone). Not a measurement.

    The generator's "truth" differs from the nominal inputs on purpose: E = 67.5 GPa, thickness 3.01 mm,
    width 25.02 mm, a clamp that rotates by M / k_theta with k_theta = 2.0e3 N m/rad, and a reading noise
    of 0.30 mm (standard deviation) on each ruler reading. Seeded: numpy.random.default_rng(seed).

    Returns a dict with 'loads' (mass_g, reading_unloaded_mm, reading_loaded_mm) rows and 'profile'
    (x_mm, reading_unloaded_mm, reading_loaded_mm) rows at 700 g, plus the true parameters."""
    rng = np.random.default_rng(seed)
    E, b, h, k_th = 67.5e9, 0.02502, 0.00301, 2.0e3
    L, a = 0.300, 0.290
    EI = E * b * h ** 3 / 12.0
    noise = 0.30e-3

    def true_w(x, m):
        P = m * G
        return deflection_point_load(x, P, a, EI) + P * a / k_th * x

    masses = np.array([0, 100, 200, 300, 400, 500, 600, 700, 800], dtype=float) / 1e3
    loads = []
    base = 0.1523                                           # ruler reading of the unloaded strip [m]
    for m in masses:
        r0 = base + rng.normal(0, noise)
        r1 = base + true_w(a, m) + rng.normal(0, noise)
        loads.append((m * 1e3, round(r0 * 1e3, 1), round(r1 * 1e3, 1)))
    xs = np.array([50, 100, 150, 200, 250, 290, 300], dtype=float) / 1e3
    prof = []
    for x in xs:
        r0 = base + rng.normal(0, noise)
        r1 = base + true_w(x, 0.700) + rng.normal(0, noise)
        prof.append((x * 1e3, round(r0 * 1e3, 1), round(r1 * 1e3, 1)))
    return {"loads": loads, "profile": prof,
            "truth": {"E": E, "b": b, "h": h, "k_theta": k_th, "L": L, "a": a, "noise": noise}}


def write_synthetic_csv(out_dir: str | Path, seed: int = 518) -> list[Path]:
    """Write the SYNTHETIC dataset as two CSV files whose names and first lines say SYNTHETIC."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = make_synthetic_cantilever_data(seed)
    paths = []
    for key, cols in (("loads", "mass_g,reading_unloaded_mm,reading_loaded_mm"),
                      ("profile", "x_mm,reading_unloaded_mm,reading_loaded_mm")):
        p = out_dir / f"SYNTHETIC_seed{seed}_{key}.csv"
        with open(p, "w", newline="") as fh:
            fh.write(f"# SYNTHETIC data from make_synthetic_cantilever_data(seed={seed}) in P05.18 solution.py: "
                     "NOT a measurement\n")
            fh.write("# strip 25 x 3 mm aluminium (nominal), free length 300 mm, load at a = 290 mm; "
                     + ("deflection read at the load point" if key == "loads" else "load 700 g, readings along the strip") + "\n")
            fh.write(cols + "\n")
            wr = csv.writer(fh)
            for row in data[key]:
                wr.writerow([f"{row[0]:g}", f"{row[1]:.1f}", f"{row[2]:.1f}"])
        paths.append(p)
    return paths


def read_csv(path: str | Path) -> np.ndarray:
    """Read one of the data CSV files: lines starting with # are comments, then one header line, then
    numbers. Returns a 2-D array of the numeric rows."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#") and ln.strip()]
    return np.loadtxt(lines[1:], delimiter=",", ndmin=2)


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Prediction, SYNTHETIC reduction and theory checks; writes cantilever_lab.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    b, u_b = 0.0250, 0.05e-3            # caliper, but the strip's width varies along its length
    h, u_h = 0.00300, 0.02e-3           # caliper
    a, u_a = 0.290, 1.0e-3
    L = 0.300
    I, u_I = rect_second_moment(b, h, u_b, u_h)
    EI = E_AL_TYPICAL * I
    print(f"strip {b*1e3:.1f} x {h*1e3:.2f} mm: I = {I*1e12:.2f} +/- {u_I*1e12:.2f} mm^4, E = {E_AL_TYPICAL/1e9:.0f} +/- "
          f"{U_E_AL/1e9:.0f} GPa (typical), EI = {EI:.4f} N m^2")
    Pmax = max_safe_load(FTY_6063_T5 / 2, a, b, h)
    print(f"safe load (clamp stress <= F_ty/2 = {FTY_6063_T5/2e6:.1f} MPa): {Pmax:.2f} N = {Pmax/G*1e3:.0f} g; "
          f"stress at 800 g: {clamp_stress(0.8*G, a, b, h)/1e6:.1f} MPa")
    q = RHO_AL * G * b * h
    print(f"self weight {q:.2f} N/m: tip sag {deflection_self_weight(L, q, L, EI)*1e3:.2f} mm before loading (removed by differencing)")

    print("\nPREDICTION (write it down and date it before loading): deflection at the load point x = a")
    print("mass [g]   delta [mm]   u (GUM) [mm]   Monte Carlo mean, sd [mm]")
    for m in (100, 400, 800):
        d, u = predict_deflection(m / 1e3 * G, a, E_AL_TYPICAL, U_E_AL, b, u_b, h, u_h, u_a)
        dm, um = predict_deflection_mc(m / 1e3 * G, a, E_AL_TYPICAL, U_E_AL, b, u_b, h, u_h, u_a)
        print(f"{m:8d}   {d*1e3:10.2f}   {u*1e3:12.2f}   {dm*1e3:8.2f}, {um*1e3:.2f}")
    s_pred = compliance_shape(a, a) / EI
    _, u_rel = predict_deflection(1.0, a, E_AL_TYPICAL, U_E_AL, b, u_b, h, u_h, u_a)
    print(f"predicted compliance at x = a: {s_pred*1e3:.3f} mm/N +/- {u_rel/(s_pred)*100:.1f} %")

    # ---- the SYNTHETIC dataset
    data = make_synthetic_cantilever_data(518)
    Ld = np.array(data["loads"])
    P = Ld[:, 0] / 1e3 * G
    d = (Ld[:, 2] - Ld[:, 1]) / 1e3
    u_d = np.sqrt(2) * 0.30e-3
    fit = fit_line(P, d, u_d)
    E_fit, u_Efit = modulus_from_slope(fit["slope"], fit["u_slope"], a, u_a, b, u_b, h, u_h)
    print("\nSYNTHETIC dataset seed 518 (not a measurement): load-deflection at x = a")
    for row, dd in zip(data["loads"], d):
        print(f"  {row[0]:5.0f} g  readings {row[1]:6.1f} -> {row[2]:6.1f} mm   delta {dd*1e3:5.1f} mm")
    print(f"slope {fit['slope']*1e3:.3f} +/- {fit['u_slope']*1e3:.3f} mm/N, intercept {fit['intercept']*1e3:+.2f} +/- "
          f"{fit['u_intercept']*1e3:.2f} mm, chi^2 = {fit['chi2']:.1f} for {fit['dof']} dof")
    print(f"E from the slope: {E_fit/1e9:.1f} +/- {u_Efit/1e9:.1f} GPa (k = 1); E_n against the typical "
          f"{E_AL_TYPICAL/1e9:.0f} +/- {U_E_AL/1e9:.0f} GPa: {en_number(E_fit, u_Efit, E_AL_TYPICAL, U_E_AL):.2f}")
    print(f"measured compliance {fit['slope']*1e3:.3f} mm/N against predicted {s_pred*1e3:.3f} mm/N: "
          f"{100*(fit['slope']/s_pred-1):+.1f} %, E_n = {en_number(fit['slope'], fit['u_slope'], s_pred, u_rel):.2f}")
    Pr = np.array(data["profile"])
    xp = Pr[:, 0] / 1e3
    wp = (Pr[:, 2] - Pr[:, 1]) / 1e3
    pf = fit_profile(xp, wp, 0.7 * G, a, u_d)
    print(f"profile at 700 g: EI = {pf['EI']:.4f} +/- {pf['u_EI']:.4f} N m^2, clamp rotation theta0 = "
          f"{pf['theta0']*1e3:.2f} +/- {pf['u_theta0']*1e3:.2f} mrad, chi^2 = {pf['chi2']:.1f} for {pf['dof']} dof")
    E_prof = pf["EI"] / I
    print(f"  E from the profile, clamp rotation removed: {E_prof/1e9:.1f} +/- {E_prof*np.hypot(pf['u_EI']/pf['EI'], u_I/I)/1e9:.1f} GPa")

    # ---- theory checks
    xf, wf = fd_cantilever(0.7 * G, a, L, EI, 200)
    print(f"\nfinite differences (n = 200) against the closed form at 700 g: max |diff| = "
          f"{np.max(np.abs(wf - deflection_point_load(xf, 0.7*G, a, EI)))*1e6:.3f} um")
    print("elastica against linear theory at the load point:")
    for m in (0.2, 0.5, 0.8, 2.0):                     # 2 kg is computed only: it would overstress this strip
        lin = float(deflection_point_load(a, m * G, a, EI))
        ela = elastica_deflection(m * G, a, EI)
        print(f"  {m*1e3:5.0f} g: linear {lin*1e3:6.2f} mm, elastica {ela*1e3:6.2f} mm, linear high by {100*(lin/ela-1):.2f} %"
              f"  (delta/a = {ela/a:.3f})")

    # ---- figure
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
    ax[0].errorbar(P, d * 1e3, yerr=2 * u_d * 1e3, fmt="o", label="SYNTHETIC data (k = 2 bars)")
    Pl = np.linspace(0, P.max(), 50)
    ax[0].plot(Pl, (fit["slope"] * Pl + fit["intercept"]) * 1e3, label="weighted fit")
    ax[0].fill_between(Pl, (s_pred - 2 * u_rel) * Pl * 1e3, (s_pred + 2 * u_rel) * Pl * 1e3, alpha=0.3,
                       label="prediction, k = 2 band")
    ax[0].set_xlabel("load P [N]"); ax[0].set_ylabel("deflection at the load [mm]")
    ax[0].set_title("SYNTHETIC: load-deflection"); ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
    xx = np.linspace(0, L, 200)
    ax[1].errorbar(xp * 1e3, wp * 1e3, yerr=2 * u_d * 1e3, fmt="o", label="SYNTHETIC profile, 700 g")
    ax[1].plot(xx * 1e3, deflection_point_load(xx, 0.7 * G, a, EI) * 1e3, label="prediction (typical E)")
    ax[1].plot(xx * 1e3, (deflection_point_load(xx, 0.7 * G, a, pf["EI"]) + pf["theta0"] * xx) * 1e3, "--",
               label="fit: cubic + clamp rotation")
    ax[1].set_xlabel("x from the clamp [mm]"); ax[1].set_ylabel("deflection [mm]")
    ax[1].set_title("SYNTHETIC: deflected shape"); ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_dir / "cantilever_lab.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'cantilever_lab.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
