"""P07.15 Home Lab: Cooling Curve: reference solution.

A small aluminium bar, warmed in tap-hot water (never above 60 C) and hung in still room air, cools
by natural convection and radiation. If the bar is small and conducts well (Biot number Bi = h Lc / k
much less than 0.1) its temperature stays uniform and follows the lumped-capacitance law

    rho c V dT/dt = -h A (T - T_inf)   =>   theta = T - T_inf = theta_0 exp(-t / tau),  tau = rho c V / (h A).

This module turns that into a prediction and a test of it:

* the prediction of tau with an uncertainty, from h = h_conv + h_rad: Churchill-Chu natural
  convection on a horizontal cylinder and its two vertical end faces (lesson P07.13, air properties
  from Lienhard Table A.6) plus the linearised radiation coefficient (lesson P07.14), by first-order
  propagation and a seeded Monte Carlo;
* two fits of a measured cooling curve: a weighted straight line through ln(theta) (log-linear) and a
  nonlinear least-squares exponential with T_inf free; a third fit includes the thermometer's own
  first-order lag; fits over successive time windows show whether h is really constant;
* the Biot number and a one-dimensional distributed model (the series solutions for an infinite
  cylinder and a slab, combined by the product rule into the short cylinder) against the lumped model;
* a nonlinear lumped model whose h follows the correlations as the bar cools (solve_ivp), with its
  energy balance closing to the solver tolerance;
* a SYNTHETIC dataset generator for learners without a kit. Synthetic data are not measurements.

Properties: 6061-T6 aluminium rho = 2700 kg/m^3, c = 896 J/(kg K), k = 167 W/(m K) (Lienhard Table
A.1, at 20 C); emittance of aluminium from Lienhard Table 10.1 (commercial sheet 0.09; polished
0.04-0.06; heavily oxidised 0.20-0.33). Units: SI in the functions (K, m, s, W); the driver prints
degrees Celsius and minutes.

Run:  python solution.py [out_dir]     (prints prediction, fits and checks, writes cooling_lab.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-07/lesson-15-home-lab-cooling-curve
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad, solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.special import j0, j1

G0 = 9.80665                  # m/s^2
SIGMA = 5.670374419e-8        # W/(m^2 K^4), CODATA 2018
# 6061-T6 aluminium at 20 C, Lienhard & Lienhard, A Heat Transfer Textbook 6th ed., Table A.1
RHO_AL, C_AL, K_AL = 2700.0, 896.0, 167.0
# AISI 1010 carbon steel and AISI 304 stainless steel at 20 C, Table A.1 (for comparisons)
RHO_1010, C_1010, K_1010 = 7830.0, 434.0, 64.0
RHO_304, C_304, K_304 = 8000.0, 400.0, 13.8
EPS_AL_NOMINAL = 0.09         # Lienhard Table 10.1, aluminium commercial sheet (90 C)
U_EPS_AL = 0.04               # our judgment for an as-bought bar: spans polished (0.05) to lightly oxidised (0.13)
U_HCONV_REL = 0.10            # assumed relative standard uncertainty of a natural-convection correlation
# Air at 1 atm, Lienhard Table A.6 (subset): T [K], nu [1e-5 m^2/s], k [W/(m K)], alpha [1e-5 m^2/s], Pr
AIR_A6 = np.array([
    [270.0, 1.304, 0.0241, 1.83, 0.711],
    [280.0, 1.392, 0.0249, 1.96, 0.710],
    [290.0, 1.482, 0.0256, 2.09, 0.708],
    [300.0, 1.575, 0.0264, 2.23, 0.707],
    [310.0, 1.670, 0.0271, 2.37, 0.706],
    [320.0, 1.766, 0.0279, 2.51, 0.705],
    [330.0, 1.865, 0.0286, 2.65, 0.704],
    [340.0, 1.966, 0.0293, 2.80, 0.703],
    [350.0, 2.069, 0.0300, 2.95, 0.702],
])
# The lesson's bar: 1 inch (25.4 mm) round 6061 bar cut 50 mm long, hung horizontally in a string sling.
D_BAR, L_BAR = 0.0254, 0.0500
T_WATER_MAX_C = 60.0          # hard limit from projects/p07-cooling-curve/safety.md


# ---------------------------------------------------------------- the object

def bar_geometry(D: float, L: float) -> dict[str, float]:
    """Volume, total area (side plus both ends), side and end areas, and the conduction length
    Lc = V/A of a solid cylinder of diameter D and length L."""
    V = np.pi * D ** 2 / 4 * L
    A_side, A_end = np.pi * D * L, 2 * np.pi * D ** 2 / 4
    A = A_side + A_end
    return {"V": V, "A": A, "A_side": A_side, "A_end": A_end, "Lc": V / A}


def air_props(T: float) -> dict[str, float]:
    """nu, k, alpha, Pr of air at 1 atm from the Table A.6 subset (270-350 K, linear interpolation),
    and beta = 1/T."""
    if not (AIR_A6[0, 0] <= T <= AIR_A6[-1, 0]):
        raise ValueError("air_props: T outside 270-350 K")
    c = lambda j: float(np.interp(T, AIR_A6[:, 0], AIR_A6[:, j]))
    return {"nu": c(1) * 1e-5, "k": c(2), "alpha": c(3) * 1e-5, "Pr": c(4), "beta": 1.0 / T}


def h_conv_bar(Ts: float, Tinf: float, D: float, L: float) -> dict[str, float]:
    """Natural convection from a horizontal bar at Ts in still air at Tinf, properties at the film
    temperature: the curved side by Churchill-Chu for a horizontal cylinder (Lienhard 8.28, length D),
    each end as a vertical plate of height D (8.13a). Returns h_side, h_end and the area-weighted h."""
    if Ts == Tinf:
        return {"h_side": 0.0, "h_end": 0.0, "h": 0.0}
    p = air_props(0.5 * (Ts + Tinf))
    Ra = G0 * p["beta"] * abs(Ts - Tinf) * D ** 3 / (p["nu"] * p["alpha"])
    Pr = p["Pr"]
    Nu_side = 0.36 + 0.518 * Ra ** 0.25 / (1 + (0.559 / Pr) ** (9 / 16)) ** (4 / 9)
    Nu_end = 0.68 + 0.670 * Ra ** 0.25 / (1 + (0.492 / Pr) ** (9 / 16)) ** (4 / 9)
    h_side, h_end = Nu_side * p["k"] / D, Nu_end * p["k"] / D
    g = bar_geometry(D, L)
    return {"h_side": h_side, "h_end": h_end, "h": (h_side * g["A_side"] + h_end * g["A_end"]) / g["A"], "Ra": Ra}


def h_rad_linear(eps: float, Ts: float, Tsur: float) -> float:
    """eps sigma (Ts^2 + Tsur^2)(Ts + Tsur): exact linearisation of eps sigma (Ts^4 - Tsur^4) (P07.14)."""
    return eps * SIGMA * (Ts ** 2 + Tsur ** 2) * (Ts + Tsur)


def biot_number(h: float, k: float, Lc: float) -> float:
    """Bi = h Lc / k (Lc = V/A for the lumped criterion Bi < 0.1)."""
    return h * Lc / k


def time_constant(rho: float, c: float, V: float, h: float, A: float) -> float:
    """tau = rho c V / (h A)."""
    return rho * c * V / (h * A)


def predict_tau(theta_rep: float, Tinf: float, D: float = D_BAR, L: float = L_BAR, eps: float = EPS_AL_NOMINAL,
                rho: float = RHO_AL, c: float = C_AL, f_conv: float = 1.0, Tsur: float | None = None) -> dict[str, float]:
    """Predicted lumped time constant with h = f_conv h_conv + h_rad evaluated at a representative excess
    temperature theta_rep (K above Tinf). Returns tau, h_conv, h_rad, h, the Biot number and rho c V."""
    Tsur = Tinf if Tsur is None else Tsur
    g = bar_geometry(D, L)
    Ts = Tinf + theta_rep
    hc = f_conv * h_conv_bar(Ts, Tinf, D, L)["h"]
    hr = h_rad_linear(eps, Ts, Tsur)
    h = hc + hr
    return {"tau": time_constant(rho, c, g["V"], h, g["A"]), "h_conv": hc, "h_rad": hr, "h": h,
            "Bi": biot_number(h, K_AL, g["Lc"]), "C": rho * c * g["V"]}


def predict_tau_uncertainty(theta_rep: float, Tinf: float, u: dict[str, float], n: int = 20000,
                            seed: int = 2026) -> dict[str, float]:
    """Uncertainty of the predicted tau, two ways. ``u`` holds standard uncertainties of D, L, eps, rho, c,
    Tinf (absolute) and f_conv (relative scale on h_conv, nominal 1). First-order (GUM): central-difference
    sensitivities, terms added in quadrature, returned per input. Monte Carlo: n normal draws
    (eps clipped to [0.01, 1]) with numpy.random.default_rng(seed); returns mean, sd."""
    nominal = {"D": D_BAR, "L": L_BAR, "eps": EPS_AL_NOMINAL, "rho": RHO_AL, "c": C_AL, "Tinf": Tinf, "f_conv": 1.0}

    def tau_of(v):
        return predict_tau(theta_rep, v["Tinf"], v["D"], v["L"], v["eps"], v["rho"], v["c"], v["f_conv"])["tau"]

    tau0 = tau_of(nominal)
    terms = {}
    for key, uk in u.items():
        if uk == 0:
            terms[key] = 0.0
            continue
        hi, lo = dict(nominal), dict(nominal)
        step = 1e-4 * (abs(nominal[key]) if nominal[key] != 0 else 1.0)
        hi[key] += step
        lo[key] -= step
        terms[key] = (tau_of(hi) - tau_of(lo)) / (2 * step) * uk
    u_gum = float(np.sqrt(sum(t ** 2 for t in terms.values())))
    rng = np.random.default_rng(seed)
    draws = {k: nominal[k] + rng.normal(0.0, u.get(k, 0.0), n) for k in nominal}
    draws["eps"] = np.clip(draws["eps"], 0.01, 1.0)
    taus = np.array([tau_of({k: draws[k][i] for k in nominal}) for i in range(n)])
    return {"tau": tau0, "u_gum": u_gum, "terms": terms, "mc_mean": float(taus.mean()), "mc_sd": float(taus.std(ddof=1))}


# ---------------------------------------------------------------- fits

def fit_loglinear(t: np.ndarray, T: np.ndarray, Tinf: float, u_T: float, t_min: float = 0.0,
                  theta_min: float = 0.5) -> dict[str, float]:
    """Weighted straight-line fit of ln(theta) = ln(theta_0) - t/tau, theta = T - Tinf, over the points with
    t >= t_min and theta >= theta_min. Each ln(theta) has standard uncertainty u_T/theta, so late, small-theta
    points get less weight. Returns tau, u_tau, theta0, u_theta0, chi2, dof, n."""
    t, T = np.asarray(t, float), np.asarray(T, float)
    th = T - Tinf
    m = (t >= t_min) & (th >= theta_min)
    tt, y, s = t[m], np.log(th[m]), u_T / th[m]
    A = np.column_stack([np.ones_like(tt), tt])
    W = 1 / s ** 2
    cov = np.linalg.inv(A.T @ (A * W[:, None]))
    beta = cov @ (A.T @ (W * y))
    r = y - A @ beta
    slope, u_slope = beta[1], np.sqrt(cov[1, 1])
    tau = -1 / slope
    return {"tau": tau, "u_tau": u_slope / slope ** 2, "theta0": float(np.exp(beta[0])),
            "u_theta0": float(np.exp(beta[0]) * np.sqrt(cov[0, 0])), "chi2": float(np.sum(W * r ** 2)),
            "dof": int(m.sum() - 2), "n": int(m.sum())}


def exp_model(t, theta0, tau, Tinf):
    """T(t) = Tinf + theta0 exp(-t / tau)."""
    return Tinf + theta0 * np.exp(-np.asarray(t, float) / tau)


def fit_exponential(t: np.ndarray, T: np.ndarray, u_T: float, t_min: float = 0.0,
                    Tinf: float | None = None) -> dict[str, float]:
    """Nonlinear least squares (scipy.optimize.curve_fit, absolute_sigma) of T = Tinf + theta0 exp(-t/tau)
    to the points with t >= t_min, with Tinf a third free parameter (Tinf None) or fixed. Returns the
    parameters, their standard uncertainties, the tau-Tinf correlation and chi2/dof."""
    t, T = np.asarray(t, float), np.asarray(T, float)
    m = t >= t_min
    tt, TT = t[m], T[m]
    T_end = TT[-1]
    p0_tau = max((tt[-1] - tt[0]) / 2, 1.0)
    if Tinf is None:
        f = exp_model
        p0 = [TT[0] - T_end + 1.0, p0_tau, T_end - 1.0]
    else:
        f = lambda x, a, b: exp_model(x, a, b, Tinf)
        p0 = [TT[0] - Tinf, p0_tau]
    popt, pcov = curve_fit(f, tt, TT, p0=p0, sigma=np.full_like(TT, u_T), absolute_sigma=True, maxfev=20000)
    r = TT - f(tt, *popt)
    perr = np.sqrt(np.diag(pcov))
    out = {"theta0": popt[0], "tau": popt[1], "u_theta0": perr[0], "u_tau": perr[1],
           "chi2": float(np.sum((r / u_T) ** 2)), "dof": int(m.sum() - len(popt))}
    if Tinf is None:
        out.update(Tinf=popt[2], u_Tinf=perr[2], corr_tau_Tinf=float(pcov[1, 2] / (perr[1] * perr[2])))
    else:
        out["Tinf"] = Tinf
    return out


def lagged_model(t, theta0, tau, tau_p, Tinf):
    """Reading of a first-order probe (time constant tau_p) attached at t = 0 at Tinf to a body whose excess
    temperature is theta0 exp(-t/tau): tau_p dTm/dt = T - Tm gives
    Tm = Tinf + theta0 tau / (tau - tau_p) (exp(-t/tau) - exp(-t/tau_p))."""
    t = np.asarray(t, float)
    return Tinf + theta0 * tau / (tau - tau_p) * (np.exp(-t / tau) - np.exp(-t / tau_p))


def fit_lagged(t: np.ndarray, T: np.ndarray, u_T: float, Tinf: float) -> dict[str, float]:
    """Fit lagged_model (theta0, tau, tau_p free; Tinf fixed) to the whole record including the first
    minutes, where the probe is still catching up. Returns parameters, uncertainties and chi2/dof."""
    t, T = np.asarray(t, float), np.asarray(T, float)
    i_pk = int(np.argmax(T))
    p0 = [T[i_pk] - Tinf, max(t[-1] / 2, 10.0), max(t[i_pk] / 3, 1.0)]
    f = lambda x, a, b, c: lagged_model(x, a, b, c, Tinf)
    popt, pcov = curve_fit(f, t, T, p0=p0, sigma=np.full_like(T, u_T), absolute_sigma=True, maxfev=20000,
                           bounds=([0, 1, 0.01], [200, 1e6, 600]))
    perr = np.sqrt(np.diag(pcov))
    r = T - f(t, *popt)
    return {"theta0": popt[0], "tau": popt[1], "tau_p": popt[2], "u_theta0": perr[0], "u_tau": perr[1],
            "u_tau_p": perr[2], "chi2": float(np.sum((r / u_T) ** 2)), "dof": int(len(t) - 3)}


def h_from_tau(tau: float, u_tau: float, rho: float, c: float, V: float, A: float,
               u_rel_rhocV_over_A: float = 0.0) -> tuple[float, float]:
    """Measured h = rho c V / (tau A) and its standard uncertainty from tau and from rho c V / A."""
    h = rho * c * V / (tau * A)
    return h, h * np.hypot(u_tau / tau, u_rel_rhocV_over_A)


def en_number(x1: float, u1: float, x2: float, u2: float, k: float = 2.0) -> float:
    """Normalised error |x1 - x2| / sqrt(U1^2 + U2^2) with expanded uncertainties U = k u (P00.15)."""
    return abs(x1 - x2) / np.hypot(k * u1, k * u2)


# ---------------------------------------------------------------- distributed model (series solutions)

def cylinder_eigenvalues(Bi: float, n: int) -> np.ndarray:
    """First n roots of lambda J1(lambda) = Bi J0(lambda) (infinite cylinder, Bi = h R / k): one root in each
    interval between successive zeros of J1 and J0, found with brentq."""
    from scipy.special import jn_zeros
    z0, z1 = jn_zeros(0, n), np.concatenate([[0.0], jn_zeros(1, n)])
    f = lambda lam: lam * j1(lam) - Bi * j0(lam)
    return np.array([brentq(f, max(z1[i], 1e-12), z0[i] - 1e-12) for i in range(n)])


def slab_eigenvalues(Bi: float, n: int) -> np.ndarray:
    """First n roots of lambda tan(lambda) = Bi (slab of half-thickness l, Bi = h l / k): one root in each
    interval (k pi, k pi + pi/2), found with brentq on lambda sin - Bi cos."""
    f = lambda lam: lam * np.sin(lam) - Bi * np.cos(lam)
    return np.array([brentq(f, k * np.pi + 1e-12, k * np.pi + np.pi / 2 - 1e-12) for k in range(n)])


def cylinder_series(rho_r: np.ndarray, Bi: float, Fo: np.ndarray, n_terms: int = 40) -> np.ndarray:
    """theta/theta_i in an infinite cylinder, uniform initial temperature, convection h at r = R:
    sum_n C_n exp(-lambda_n^2 Fo) J0(lambda_n r/R), C_n = 2 J1(l) / (l (J0(l)^2 + J1(l)^2)), Fo = alpha t / R^2.
    Returns an array of shape broadcast(Fo, rho_r)."""
    lam = cylinder_eigenvalues(Bi, n_terms)
    C = 2 * j1(lam) / (lam * (j0(lam) ** 2 + j1(lam) ** 2))
    Fo = np.asarray(Fo, float)[..., None]
    r = np.asarray(rho_r, float)[..., None]
    return np.sum(C * np.exp(-lam ** 2 * Fo) * j0(lam * r), axis=-1)


def cylinder_series_mean(Bi: float, Fo: np.ndarray, n_terms: int = 40) -> np.ndarray:
    """Area-mean theta/theta_i of the infinite cylinder: sum_n C_n exp(-lambda_n^2 Fo) 2 J1(l)/l (energy content)."""
    lam = cylinder_eigenvalues(Bi, n_terms)
    C = 2 * j1(lam) / (lam * (j0(lam) ** 2 + j1(lam) ** 2))
    Fo = np.asarray(Fo, float)[..., None]
    return np.sum(C * np.exp(-lam ** 2 * Fo) * 2 * j1(lam) / lam, axis=-1)


def slab_series(xi: np.ndarray, Bi: float, Fo: np.ndarray, n_terms: int = 40) -> np.ndarray:
    """theta/theta_i in a slab of half-thickness l (both faces convecting, xi = x/l from the mid-plane):
    sum_n C_n exp(-lambda_n^2 Fo) cos(lambda_n xi), C_n = 4 sin(l) / (2 l + sin(2 l)), Fo = alpha t / l^2."""
    lam = slab_eigenvalues(Bi, n_terms)
    C = 4 * np.sin(lam) / (2 * lam + np.sin(2 * lam))
    Fo = np.asarray(Fo, float)[..., None]
    x = np.asarray(xi, float)[..., None]
    return np.sum(C * np.exp(-lam ** 2 * Fo) * np.cos(lam * x), axis=-1)


def short_cylinder(r_over_R: float, z_over_l: float, t: np.ndarray, h: float, D: float, L: float,
                   k: float, alpha: float, n_terms: int = 40) -> np.ndarray:
    """theta/theta_i in a finite cylinder (diameter D, length L = 2 l) with the same h on every face, by the
    product rule: (infinite cylinder of radius R) x (slab of half-thickness l)."""
    R, l = D / 2, L / 2
    t = np.asarray(t, float)
    return cylinder_series(r_over_R, h * R / k, alpha * t / R ** 2, n_terms) * slab_series(z_over_l, h * l / k, alpha * t / l ** 2, n_terms)


def lumped_theta(t: np.ndarray, h: float, rho: float, c: float, V: float, A: float) -> np.ndarray:
    """theta/theta_i = exp(-t / tau) of the lumped model."""
    return np.exp(-np.asarray(t, float) / time_constant(rho, c, V, h, A))


# ---------------------------------------------------------------- nonlinear lumped model (h varies as the bar cools)

def simulate_cooling(t_eval: np.ndarray, T0: float, Tinf: float, D: float = D_BAR, L: float = L_BAR,
                     eps: float = EPS_AL_NOMINAL, rho: float = RHO_AL, c: float = C_AL, f_conv: float = 1.0,
                     Tsur: float | None = None) -> dict[str, np.ndarray]:
    """Lumped bar with h_conv(T) from the correlations (scaled by f_conv) and exact gray radiation:
    rho c V dT/dt = -f_conv h_conv(T) A (T - Tinf) - eps sigma A (T^4 - Tsur^4). solve_ivp (LSODA,
    rtol 1e-10) also integrates the heat lost by each mode, so the energy balance
    rho c V (T0 - T(t)) = Q_conv(t) + Q_rad(t) can be checked. Returns t, T, Q_conv, Q_rad."""
    Tsur = Tinf if Tsur is None else Tsur
    g = bar_geometry(D, L)
    C = rho * c * g["V"]

    def rhs(_, y):
        T = y[0]
        qc = f_conv * h_conv_bar(T, Tinf, D, L)["h"] * g["A"] * (T - Tinf)
        qr = eps * SIGMA * g["A"] * (T ** 4 - Tsur ** 4)
        return [-(qc + qr) / C, qc, qr]

    sol = solve_ivp(rhs, (0.0, float(np.max(t_eval))), [T0, 0.0, 0.0], t_eval=t_eval, method="LSODA",
                    rtol=1e-10, atol=1e-10)
    return {"t": sol.t, "T": sol.y[0], "Q_conv": sol.y[1], "Q_rad": sol.y[2], "C": C}


def fit_variable_h(t: np.ndarray, T: np.ndarray, u_T: float, Tinf: float, eps: float = EPS_AL_NOMINAL,
                   t_min: float = 120.0, D: float = D_BAR, L: float = L_BAR) -> dict[str, float]:
    """Physics-based nonlinear fit: simulate_cooling started at t_min with two free parameters, the
    temperature T_start at t_min and the scale factor f_conv on the correlation's h_conv (eps, Tinf, the
    geometry and the material fixed), by scipy.optimize.least_squares. f_conv = 1 means the correlation
    was right; its uncertainty comes from the Jacobian, scaled by sqrt(chi2/dof) if that exceeds one.
    Returns T_start, f_conv, u_f_conv, chi2, dof."""
    from scipy.optimize import least_squares
    t, T = np.asarray(t, float), np.asarray(T, float)
    m = t >= t_min
    tt, TT = t[m], T[m]

    def resid(p):
        sim = simulate_cooling(tt - tt[0], p[0], Tinf, D, L, eps, f_conv=p[1])
        return (sim["T"] - TT) / u_T

    sol = least_squares(resid, [TT[0], 1.0], x_scale=[1.0, 0.1], method="lm")
    chi2 = float(np.sum(sol.fun ** 2))
    dof = int(m.sum() - 2)
    cov = np.linalg.inv(sol.jac.T @ sol.jac) * max(1.0, chi2 / dof)
    return {"T_start": sol.x[0], "f_conv": sol.x[1], "u_f_conv": float(np.sqrt(cov[1, 1])), "chi2": chi2, "dof": dof}


def fitted_model_energy_balance(theta0: float, tau: float, t_end: float, C: float) -> tuple[float, float]:
    """Energy balance of a fitted exponential: the change in stored energy C (theta(0) - theta(t_end)) and the
    heat lost, int_0^t_end h A theta dt with h A = C / tau, by adaptive quadrature. They must agree."""
    stored = C * theta0 * (1 - np.exp(-t_end / tau))
    lost = quad(lambda s: (C / tau) * theta0 * np.exp(-s / tau), 0.0, t_end, epsabs=0.0, epsrel=1e-13)[0]
    return stored, lost


# ---------------------------------------------------------------- SYNTHETIC data

def make_synthetic_cooling_data(seed: int = 715) -> dict:
    """SYNTHETIC cooling curve of the lesson's bar (25.4 mm x 50 mm 6061 aluminium, hung horizontally in
    still air), read every 15 s for 60 min by a probe thermometer of 0.1 C resolution taped to its surface.
    Not a measurement.

    The generator's hidden "truth" differs from the nominal inputs on purpose: convection 12 % stronger
    than the correlation (f_conv = 1.12, a room with gentle draughts), emittance 0.12 (nominal 0.09), a bar
    0.3 mm longer (50.3 mm), room air at 21.3 C, the bar at 49.6 C when the probe is attached at t = 0, and
    a probe that lags with a 20 s time constant starting from room temperature. Each reading adds normal
    noise (sd 0.04 C) and is rounded to 0.1 C (the thermometer's quantization). The room thermometer reads
    the air at 0, 30 and 60 min with the same noise and rounding. Seeded: numpy.random.default_rng(seed).

    Returns {'cooling': rows (t_s, T_C), 'ambient': rows (t_s, T_air_C), 'truth': dict}."""
    rng = np.random.default_rng(seed)
    truth = {"f_conv": 1.12, "eps": 0.12, "D": D_BAR, "L": 0.0503, "Tinf_C": 21.3, "T0_C": 49.6,
             "tau_p": 20.0, "noise_C": 0.04, "resolution_C": 0.1}
    Tinf = truth["Tinf_C"] + 273.15
    t_fine = np.linspace(0.0, 3600.0, 3601)
    sim = simulate_cooling(t_fine, truth["T0_C"] + 273.15, Tinf, truth["D"], truth["L"], truth["eps"],
                           f_conv=truth["f_conv"])
    # probe lag: tau_p dTm/dt = T - Tm, Tm(0) = Tinf, integrated exactly on the fine grid (T piecewise linear)
    T = sim["T"]
    Tm = np.empty_like(T)
    Tm[0] = Tinf
    dt = t_fine[1] - t_fine[0]
    a = np.exp(-dt / truth["tau_p"])
    for i in range(1, len(T)):
        # exact update for linear T over the step
        slope = (T[i] - T[i - 1]) / dt
        Tm[i] = T[i] - slope * truth["tau_p"] + (Tm[i - 1] - T[i - 1] + slope * truth["tau_p"]) * a
    t_read = np.arange(0.0, 3600.0 + 1e-9, 15.0)
    Tm_read = np.interp(t_read, t_fine, Tm) - 273.15
    q = truth["resolution_C"]
    rows = [(float(tt), round(float(np.round((x + rng.normal(0, truth["noise_C"])) / q) * q), 1))
            for tt, x in zip(t_read, Tm_read)]
    amb = [(float(tt), round(float(np.round((truth["Tinf_C"] + rng.normal(0, truth["noise_C"])) / q) * q), 1))
           for tt in (0.0, 1800.0, 3600.0)]
    return {"cooling": rows, "ambient": amb, "truth": truth}


def write_synthetic_csv(out_dir: str | Path, seed: int = 715) -> list[Path]:
    """Write the SYNTHETIC dataset as two CSV files whose names and first lines say SYNTHETIC."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = make_synthetic_cooling_data(seed)
    paths = []
    for key, cols, what in (("cooling", "t_s,T_C", "bar surface temperature, probe taped on, read every 15 s"),
                            ("ambient", "t_s,T_air_C", "room air temperature, second thermometer, 1 m from the bar")):
        p = out_dir / f"SYNTHETIC_seed{seed}_{key}.csv"
        with open(p, "w", newline="") as fh:
            fh.write(f"# SYNTHETIC data from make_synthetic_cooling_data(seed={seed}) in P07.15 solution.py: "
                     "NOT a measurement\n")
            fh.write(f"# 25.4 mm x 50 mm 6061 aluminium bar (nominal), hung horizontally in still air; {what}; "
                     "0.1 C resolution\n")
            fh.write(cols + "\n")
            wr = csv.writer(fh)
            for row in data[key]:
                wr.writerow([f"{row[0]:g}", f"{row[1]:.1f}"])
        paths.append(p)
    return paths


def read_csv(path: str | Path) -> np.ndarray:
    """Read a data CSV: lines starting with # are comments, then one header line, then numbers."""
    with open(path, encoding="utf-8") as fh:
        lines = [ln for ln in fh if not ln.startswith("#") and ln.strip()]
    return np.loadtxt(lines[1:], delimiter=",", ndmin=2)


# ---------------------------------------------------------------- driver

def main(out_dir: str | Path = ".") -> None:
    """Prediction, SYNTHETIC reduction, distributed-model check; writes cooling_lab.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    g = bar_geometry(D_BAR, L_BAR)
    C = RHO_AL * C_AL * g["V"]
    print(f"bar {D_BAR * 1e3:.1f} mm x {L_BAR * 1e3:.1f} mm 6061-T6: V = {g['V'] * 1e6:.3f} cm^3, mass "
          f"{RHO_AL * g['V'] * 1e3:.1f} g, A = {g['A'] * 1e4:.2f} cm^2, Lc = V/A = {g['Lc'] * 1e3:.3f} mm, "
          f"rho c V = {C:.2f} J/K")

    # ---- prediction (before measuring)
    Tinf_nom = 294.15                          # 21.0 C, the room as read before the run
    theta_rep = 14.0                           # K: mean excess over a run from about 50 C down to 25 C
    pr = predict_tau(theta_rep, Tinf_nom)
    print(f"\nPREDICTION (write it down and date it before measuring), at theta = {theta_rep:.0f} K, air 21.0 C:")
    hc = h_conv_bar(Tinf_nom + theta_rep, Tinf_nom, D_BAR, L_BAR)
    print(f"  Ra_D = {hc['Ra']:.3e}; h_side {hc['h_side']:.2f}, h_end {hc['h_end']:.2f}, h_conv {pr['h_conv']:.2f} W/m^2 K; "
          f"h_rad (eps {EPS_AL_NOMINAL}) {pr['h_rad']:.3f}; h {pr['h']:.2f} W/m^2 K")
    print(f"  Bi = {pr['Bi']:.2e} (lumped if < 0.1); tau = {pr['tau']:.0f} s = {pr['tau'] / 60:.1f} min")
    for th in (5.0, 14.0, 28.0):
        p_ = predict_tau(th, Tinf_nom)
        print(f"  at theta = {th:4.0f} K: h_conv {p_['h_conv']:.2f}, h_rad {p_['h_rad']:.3f}, tau {p_['tau'] / 60:.1f} min")
    u_in = {"D": 0.05e-3, "L": 0.1e-3, "eps": U_EPS_AL, "rho": 0.005 * RHO_AL, "c": 0.01 * C_AL, "Tinf": 0.5,
            "f_conv": U_HCONV_REL}
    pu = predict_tau_uncertainty(theta_rep, Tinf_nom, u_in, n=4000)
    print(f"  tau = {pu['tau']:.0f} +/- {pu['u_gum']:.0f} s (GUM, k = 1); Monte Carlo {pu['mc_mean']:.0f}, sd {pu['mc_sd']:.0f} s")
    print("  budget [s]: " + ", ".join(f"{k} {abs(v):.1f}" for k, v in pu["terms"].items()))
    eps_hi = predict_tau(theta_rep, Tinf_nom, eps=0.9)
    print(f"  painted or taped black (eps 0.9, illustrative): h_rad {eps_hi['h_rad']:.2f}, tau {eps_hi['tau'] / 60:.1f} min")

    # ---- SYNTHETIC reduction
    data = make_synthetic_cooling_data(715)
    cool = np.array(data["cooling"])
    amb = np.array(data["ambient"])
    t, T = cool[:, 0], cool[:, 1] + 273.15
    Tinf_meas = float(np.mean(amb[:, 1])) + 273.15
    u_T = float(np.hypot(0.1 / np.sqrt(12), 0.04))   # quantization (resolution / sqrt 12) and reading noise
    print(f"\nSYNTHETIC dataset seed 715 (not a measurement): {len(t)} readings every 15 s; first readings "
          f"{', '.join(f'{x:.1f}' for x in cool[:6, 1])} C; room {', '.join(f'{x:.1f}' for x in amb[:, 1])} C "
          f"(mean {Tinf_meas - 273.15:.2f} C); reading u = {u_T:.3f} C")
    ll = fit_loglinear(t, T, Tinf_meas, u_T, t_min=120.0, theta_min=1.0)
    print(f"log-linear fit (t >= 2 min, theta >= 1 K, Tinf = room mean): tau = {ll['tau']:.1f} +/- {ll['u_tau']:.1f} s, "
          f"theta0 = {ll['theta0']:.2f} K, chi^2 = {ll['chi2']:.0f} for {ll['dof']} dof")
    ne = fit_exponential(t, T, u_T, t_min=120.0)
    print(f"nonlinear fit, Tinf free: tau = {ne['tau']:.1f} +/- {ne['u_tau']:.1f} s, Tinf = {ne['Tinf'] - 273.15:.2f} +/- "
          f"{ne['u_Tinf']:.2f} C, corr(tau, Tinf) = {ne['corr_tau_Tinf']:.2f}, chi^2 = {ne['chi2']:.0f} for {ne['dof']} dof")
    nf = fit_exponential(t, T, u_T, t_min=120.0, Tinf=Tinf_meas)
    print(f"nonlinear fit, Tinf fixed at room mean: tau = {nf['tau']:.1f} +/- {nf['u_tau']:.1f} s, chi^2 = {nf['chi2']:.0f} for {nf['dof']} dof")
    lg = fit_lagged(t, T, u_T, Tinf_meas)
    print(f"lagged fit, whole record: tau = {lg['tau']:.1f} +/- {lg['u_tau']:.1f} s, probe tau_p = {lg['tau_p']:.1f} +/- "
          f"{lg['u_tau_p']:.1f} s, theta0 = {lg['theta0']:.2f} K, chi^2 = {lg['chi2']:.0f} for {lg['dof']} dof")
    print("windowed log-linear fits (is h constant?):")
    for t0, t1 in ((120, 1200), (1200, 2400), (2400, 3600)):
        m = (t >= t0) & (t <= t1)
        w = fit_loglinear(t[m], T[m], Tinf_meas, u_T, theta_min=1.0)
        th_mid = float(np.interp(0.5 * (t0 + t1), t, T)) - Tinf_meas
        print(f"  {t0 / 60:4.0f}-{t1 / 60:4.0f} min (theta ~ {th_mid:4.1f} K): tau = {w['tau']:.0f} +/- {w['u_tau']:.0f} s, "
              f"predicted {predict_tau(th_mid, Tinf_meas)['tau']:.0f} s")
    birge = np.sqrt(ll["chi2"] / ll["dof"])
    print(f"chi^2/dof = {ll['chi2'] / ll['dof']:.1f}: one exponential does not fit to the reading noise; "
          f"Birge-scaled u(tau) = {ll['u_tau'] * birge:.1f} s")
    for eps_f in (EPS_AL_NOMINAL, 0.05, 0.13):
        vh = fit_variable_h(t, T, u_T, Tinf_meas, eps=eps_f)
        print(f"variable-h fit (correlation x f_conv, eps fixed at {eps_f:.2f}): f_conv = {vh['f_conv']:.3f} +/- "
              f"{vh['u_f_conv']:.3f}, chi^2 = {vh['chi2']:.0f} for {vh['dof']} dof")
    u_rel_geom = np.hypot(np.hypot(0.005, 0.01), 0.05e-3 / D_BAR)
    h_m, u_hm = h_from_tau(ll["tau"], ll["u_tau"], RHO_AL, C_AL, g["V"], g["A"], u_rel_geom)
    th_bar = 14.0
    pm = predict_tau(th_bar, Tinf_meas)
    u_hp = np.hypot(U_HCONV_REL * pm["h_conv"], pm["h_rad"] * U_EPS_AL / EPS_AL_NOMINAL)
    print(f"measured h = rho c V/(tau A) = {h_m:.2f} +/- {u_hm:.2f} W/m^2 K; predicted h_conv + h_rad = {pm['h']:.2f} +/- "
          f"{u_hp:.2f}; E_n = {en_number(h_m, u_hm, pm['h'], u_hp):.2f}")
    print(f"tau measured {ll['tau']:.0f} +/- {ll['u_tau']:.0f} s against predicted {pu['tau']:.0f} +/- {pu['u_gum']:.0f} s: "
          f"E_n = {en_number(ll['tau'], ll['u_tau'], pu['tau'], pu['u_gum']):.2f}")
    print(f"Biot number with the measured h: {biot_number(h_m, K_AL, g['Lc']):.2e}")
    st, lo = fitted_model_energy_balance(ll["theta0"], ll["tau"], 3600.0, C)
    print(f"energy balance of the fitted model over 60 min: stored change {st:.4f} J, heat lost {lo:.4f} J, "
          f"difference {st - lo:.1e} J")
    tr = data["truth"]
    sim = simulate_cooling(np.linspace(0, 3600, 361), tr["T0_C"] + 273.15, tr["Tinf_C"] + 273.15, tr["D"], tr["L"],
                           tr["eps"], f_conv=tr["f_conv"])
    print(f"generator truth (hidden in a real run): f_conv {tr['f_conv']}, eps {tr['eps']}, L {tr['L'] * 1e3:.1f} mm, "
          f"room {tr['Tinf_C']} C, probe tau_p {tr['tau_p']:.0f} s; at theta = 14 K its tau is "
          f"{predict_tau(14.0, tr['Tinf_C'] + 273.15, tr['D'], tr['L'], tr['eps'], f_conv=tr['f_conv'])['tau']:.0f} s")
    print(f"variable-h model energy balance: rho c V dT {sim['C'] * (sim['T'][0] - sim['T'][-1]):.4f} J, "
          f"Q_conv + Q_rad {sim['Q_conv'][-1] + sim['Q_rad'][-1]:.4f} J (radiation {100 * sim['Q_rad'][-1] / (sim['Q_conv'][-1] + sim['Q_rad'][-1]):.1f} %)")

    # ---- distributed model against lumped
    alpha = K_AL / (RHO_AL * C_AL)
    tt = np.array([60.0, 600.0, 1800.0, 3600.0])
    lump = lumped_theta(tt, pr["h"], RHO_AL, C_AL, g["V"], g["A"])
    cen = short_cylinder(0.0, 0.0, tt, pr["h"], D_BAR, L_BAR, K_AL, alpha)
    cor = short_cylinder(1.0, 1.0, tt, pr["h"], D_BAR, L_BAR, K_AL, alpha)
    print(f"\ndistributed (short-cylinder product series) against lumped, h = {pr['h']:.2f} W/m^2 K:")
    for i, ti in enumerate(tt):
        print(f"  t = {ti / 60:4.0f} min: lumped {lump[i]:.6f}, centre {cen[i]:.6f}, edge/corner {cor[i]:.6f}, "
              f"centre - corner {cen[i] - cor[i]:.1e}")
    print("lumped error at the surface (infinite cylinder, theta = 0.1 time):")
    for Bi in (0.001, 0.01, 0.1, 1.0):
        Fo = brentq(lambda f: float(cylinder_series_mean(Bi, f)) - 0.1, 1e-4, 1e5)
        surf = float(cylinder_series(1.0, Bi, Fo))
        lumpv = np.exp(-2 * Bi * Fo)
        print(f"  Bi = {Bi:5.3f}: mean 0.1 at Fo = {Fo:.4g}; surface {surf:.4f}, centre {float(cylinder_series(0.0, Bi, Fo)):.4f}, "
              f"lumped {lumpv:.4f}")

    # ---- figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    ax[0].plot(t / 60, T - 273.15, ".", ms=3, label="SYNTHETIC readings")
    ax[0].plot(t / 60, exp_model(t, ne["theta0"], ne["tau"], ne["Tinf"]) - 273.15, label="exponential fit, Tinf free")
    ax[0].plot(t / 60, lagged_model(t, lg["theta0"], lg["tau"], lg["tau_p"], Tinf_meas) - 273.15, "--", label="fit with probe lag")
    ax[0].set_xlabel("time (min)"); ax[0].set_ylabel("temperature (C)"); ax[0].set_title("SYNTHETIC: cooling curve")
    ax[0].legend(fontsize=8)
    th = T - Tinf_meas
    m = th > 0.05
    ax[1].semilogy(t[m] / 60, th[m], ".", ms=3, label="SYNTHETIC theta")
    ax[1].semilogy(t / 60, ll["theta0"] * np.exp(-t / ll["tau"]), label=f"log-linear fit, tau {ll['tau']:.0f} s")
    ax[1].semilogy(t / 60, ll["theta0"] * np.exp(-t / pu["tau"]), ":", label=f"prediction, tau {pu['tau']:.0f} s")
    ax[1].set_xlabel("time (min)"); ax[1].set_ylabel("T - Tinf (K)"); ax[1].set_title("SYNTHETIC: semi-log")
    ax[1].legend(fontsize=8)
    Fo = np.linspace(0, 3, 300)
    for Bi in (0.01, 0.1, 1.0):
        ax[2].plot(Fo, cylinder_series(1.0, Bi, Fo), label=f"surface, Bi {Bi}")
        ax[2].plot(Fo, np.exp(-2 * Bi * Fo), "k:", lw=0.8)
    ax[2].set_xlabel("Fo = alpha t / R^2"); ax[2].set_ylabel("theta / theta_i")
    ax[2].set_title("infinite cylinder: series (solid) vs lumped (dotted)"); ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "cooling_lab.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'cooling_lab.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
