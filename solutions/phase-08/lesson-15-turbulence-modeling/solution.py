"""P08.15 Turbulence Modeling: reference solution.

What does a turbulence model close, and how much can two constants predict? Averaging the
Navier-Stokes equations (RANS) leaves one unknown, the Reynolds shear stress -rho <u'v'>. This
module closes it two ways for fully developed flow between parallel plates (a plane channel, the
model of a flat cabin-air duct):

* Prandtl's mixing length with van Driest's near-wall damping, l+ = kappa y+ (1 - exp(-y+/A+)),
  optionally with Nikuradse's outer-region formula. The mean momentum balance gives the total shear
  stress exactly, tau+ = 1 - y+/Re_tau, so the velocity gradient is the positive root of a quadratic
  at every point and the profile follows by quadrature (P03.04). The model is then compared with the
  log law (by a least-squares fit of kappa and B) and with friction laws: the log-law friction law and
  Dean's correlation c_f = 0.073 Re_b^(-1/4) (a correlation, not a measurement).
* The standard k-epsilon model (Launder & Spalding constants C_mu = 0.09, C_eps1 = 1.44,
  C_eps2 = 1.92, sigma_k = 1.0, sigma_eps = 1.3): in the log layer, where production equals
  dissipation, k = u_tau^2/sqrt(C_mu), eps = u_tau^3/(kappa y) and nu_t = kappa u_tau y. The eps
  equation is satisfied there only if kappa^2 = sigma_eps sqrt(C_mu) (C_eps2 - C_eps1); the module
  checks that relation by substituting the profiles into the transport equations with finite
  differences.

Wall units: u_tau = sqrt(tau_w/rho), y+ = y u_tau/nu, u+ = U/u_tau, Re_tau = h u_tau/nu with h the
half-height. Bulk Reynolds number Re_b = U_b (2h)/nu (full height), c_f = tau_w/(rho U_b^2/2).

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes turbulence_models.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-15-turbulence-modeling
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

KAPPA = 0.41          # von Karman constant (typical value; published fits range about 0.38-0.42)
B_LOG = 5.0           # log-law additive constant (typical; 5.0-5.5 in the literature)
A_PLUS = 26.0         # van Driest damping constant (his 1956 value)
C_MU, C_EPS1, C_EPS2, SIGMA_K, SIGMA_EPS = 0.09, 1.44, 1.92, 1.0, 1.3   # standard k-epsilon constants
RHO_AIR, MU_AIR = 1.204, 1.81e-5          # air at 20 C, 101 325 Pa (kg/m^3, Pa s)


# ---------------------------------------------------------------- given helpers

def channel_grid(re_tau: float, n: int = 2001, beta: float = 2.5) -> np.ndarray:
    """y+ from the wall (0) to the centreline (Re_tau) on n points clustered towards the wall with a
    tanh stretching of strength beta (beta -> 0 gives a uniform grid)."""
    s = np.linspace(0.0, 1.0, n)
    if beta == 0:
        return re_tau * s
    return re_tau * (1.0 + np.tanh(beta * (s - 1.0)) / np.tanh(beta))


def cumtrapz0(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoid integral of f over x, starting at 0."""
    out = np.zeros_like(f, dtype=float)
    out[1:] = np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))
    return out


def log_law(y_plus: np.ndarray, kappa: float = KAPPA, B: float = B_LOG) -> np.ndarray:
    """The log law u+ = ln(y+)/kappa + B (valid roughly for 30 < y+ < 0.2 Re_tau)."""
    return np.log(y_plus) / kappa + B


def dean_cf(re_b: float | np.ndarray) -> float | np.ndarray:
    """Dean's (1978) correlation for fully developed turbulent channel flow, c_f = 0.073 Re_b^(-1/4),
    Re_b = U_b (2h)/nu, fitted to experiments for roughly 6e3 < Re_b < 6e5. A correlation."""
    return 0.073 * np.asarray(re_b, dtype=float) ** -0.25


# ---------------------------------------------------------------- mixing-length model

def total_stress(y_plus: np.ndarray, re_tau: float) -> np.ndarray:
    """Total (viscous + turbulent) shear stress in wall units from the mean x-momentum balance of a
    fully developed channel: tau+ = 1 - y+/Re_tau (linear, 1 at the wall, 0 on the centreline)."""
    return 1.0 - np.asarray(y_plus, dtype=float) / re_tau


def mixing_length(y_plus: np.ndarray, re_tau: float, kappa: float = KAPPA, A_plus: float = A_PLUS,
                  outer: str = "prandtl") -> np.ndarray:
    """Mixing length in wall units with van Driest damping D = 1 - exp(-y+/A+).
    outer = "prandtl": l+ = kappa y+ D.  outer = "nikuradse": l+ = Re_tau [0.14 - 0.08 (1 - eta)^2
    - 0.06 (1 - eta)^4] D with eta = y+/Re_tau (Nikuradse's pipe fit, which tends to 0.4 y+ near the
    wall and caps l at 0.14 h on the centreline). A_plus = 0 switches the damping off."""
    yp = np.asarray(y_plus, dtype=float)
    D = 1.0 - np.exp(-yp / A_plus) if A_plus > 0 else np.ones_like(yp)
    if outer == "prandtl":
        return kappa * yp * D
    if outer == "nikuradse":
        e = 1.0 - yp / re_tau
        return re_tau * (0.14 - 0.08 * e ** 2 - 0.06 * e ** 4) * D
    raise ValueError("outer must be 'prandtl' or 'nikuradse'")


def velocity_gradient(tau_plus: np.ndarray, l_plus: np.ndarray) -> np.ndarray:
    """du+/dy+ from tau+ = (1 + l+^2 |du+/dy+|) du+/dy+, the positive root of the quadratic
    l+^2 g^2 + g - tau+ = 0, written in the cancellation-free form g = 2 tau+/(1 + sqrt(1 + 4 l+^2 tau+)).
    Reduces to g = tau+ (laminar) when l+ = 0."""
    tau = np.asarray(tau_plus, dtype=float)
    l2 = np.asarray(l_plus, dtype=float) ** 2
    return 2.0 * tau / (1.0 + np.sqrt(1.0 + 4.0 * l2 * np.clip(tau, 0.0, None)))


def mixing_length_profile(re_tau: float, n: int = 2001, kappa: float = KAPPA, A_plus: float = A_PLUS,
                          outer: str = "prandtl", beta: float = 2.5) -> dict:
    """Integrate the mixing-length model across the half-channel. Returns a dict with y_plus, u_plus,
    dudy (du+/dy+), nu_t (nu_t/nu = l+^2 du+/dy+), reynolds_stress (-<u'v'>+ = nu_t du+/dy+),
    tau (total stress) and l_plus."""
    yp = channel_grid(re_tau, n, beta)
    tau = total_stress(yp, re_tau)
    lp = mixing_length(yp, re_tau, kappa, A_plus, outer)
    g = velocity_gradient(tau, lp)
    nu_t = lp ** 2 * g
    return {"y_plus": yp, "u_plus": cumtrapz0(g, yp), "dudy": g, "nu_t": nu_t,
            "reynolds_stress": nu_t * g, "tau": tau, "l_plus": lp, "re_tau": re_tau}


def laminar_profile(y_plus: np.ndarray, re_tau: float) -> np.ndarray:
    """Plane Poiseuille flow in wall units: u+ = y+ - y+^2/(2 Re_tau) (the kappa -> 0 limit)."""
    yp = np.asarray(y_plus, dtype=float)
    return yp - yp ** 2 / (2.0 * re_tau)


def bulk_velocity(prof: dict) -> float:
    """Bulk (mean) velocity in wall units, U_b+ = (1/Re_tau) integral of u+ dy+ over the half-channel."""
    return float(np.trapezoid(prof["u_plus"], prof["y_plus"]) / prof["re_tau"])


def friction(prof: dict) -> tuple[float, float]:
    """Bulk Reynolds number Re_b = 2 U_b+ Re_tau and skin-friction coefficient c_f = 2/U_b+^2 of a
    profile (both follow from the definitions of the wall units)."""
    ub = bulk_velocity(prof)
    return 2.0 * ub * prof["re_tau"], 2.0 / ub ** 2


def fit_log_law(y_plus: np.ndarray, u_plus: np.ndarray, lo: float = 30.0,
                hi: float | None = None) -> tuple[float, float]:
    """Least-squares fit of u+ = ln(y+)/kappa + B over lo <= y+ <= hi (default hi = 0.2 y+_max, the
    usual outer limit of the log layer). Returns (kappa, B)."""
    yp, up = np.asarray(y_plus), np.asarray(u_plus)
    hi = 0.2 * yp.max() if hi is None else hi
    m = (yp >= lo) & (yp <= hi)
    A = np.column_stack([np.log(yp[m]), np.ones(m.sum())])
    slope, B = np.linalg.lstsq(A, up[m], rcond=None)[0]
    return float(1.0 / slope), float(B)


def log_law_friction(re_tau: float, kappa: float = KAPPA, B: float = B_LOG) -> tuple[float, float]:
    """Friction law from integrating the log law across the whole half-channel (viscous layer and wake
    neglected): U_b+ = ln(Re_tau)/kappa + B - 1/kappa. Returns (Re_b, c_f)."""
    ub = np.log(re_tau) / kappa + B - 1.0 / kappa
    return 2.0 * ub * re_tau, 2.0 / ub ** 2


def re_tau_for_re_b(re_b: float, **model) -> float:
    """Invert the mixing-length model: the Re_tau whose profile has bulk Reynolds number re_b
    (brentq on Re_b(Re_tau) - re_b; Re_b is monotonic in Re_tau). Keyword arguments go to
    ``mixing_length_profile``."""
    f = lambda rt: friction(mixing_length_profile(rt, **model))[0] - re_b
    return float(brentq(f, 20.0, 1e5, xtol=1e-10, rtol=1e-12))


def duct_flow(V_bulk: float, gap: float, rho: float = RHO_AIR, mu: float = MU_AIR, **model) -> dict:
    """A wide flat duct of height `gap` (= 2h) with bulk velocity V_bulk, by the mixing-length model.
    Returns Re_b, Re_tau, u_tau (m/s), tau_w (Pa), dp/dx (Pa/m, = -tau_w/h), c_f, the Dean c_f, and the
    first-cell heights (m) for y+ = 1 and y+ = 30 (the two classic CFD near-wall targets)."""
    nu = mu / rho
    re_b = V_bulk * gap / nu
    rt = re_tau_for_re_b(re_b, **model)
    u_tau = rt * nu / (gap / 2.0)
    tau_w = rho * u_tau ** 2
    return {"re_b": re_b, "re_tau": rt, "u_tau": u_tau, "tau_w": tau_w, "dpdx": -tau_w / (gap / 2.0),
            "cf": tau_w / (0.5 * rho * V_bulk ** 2), "cf_dean": float(dean_cf(re_b)),
            "y1": nu / u_tau, "y30": 30.0 * nu / u_tau}


# ---------------------------------------------------------------- k-epsilon in the log layer

def keps_log_layer(y: np.ndarray, u_tau: float, kappa: float = KAPPA, c_mu: float = C_MU) -> dict:
    """The equilibrium log-layer solution of the k-epsilon model: constant stress u_tau^2, production
    = dissipation. k = u_tau^2/sqrt(C_mu), eps = u_tau^3/(kappa y), dU/dy = u_tau/(kappa y),
    nu_t = C_mu k^2/eps = kappa u_tau y. Returns a dict of arrays k, eps, nu_t, dudy."""
    y = np.asarray(y, dtype=float)
    k = np.full_like(y, u_tau ** 2 / np.sqrt(c_mu))
    eps = u_tau ** 3 / (kappa * y)
    return {"k": k, "eps": eps, "nu_t": c_mu * k ** 2 / eps, "dudy": u_tau / (kappa * y)}


def keps_residuals(y: np.ndarray, u_tau: float, kappa: float = KAPPA, c_mu: float = C_MU,
                   c1: float = C_EPS1, c2: float = C_EPS2, sigma_k: float = SIGMA_K,
                   sigma_eps: float = SIGMA_EPS) -> tuple[np.ndarray, np.ndarray]:
    """Substitute the log-layer profiles into the steady, fully developed k and eps equations,
        0 = d/dy(nu_t/sigma_k dk/dy) + P - eps,
        0 = d/dy(nu_t/sigma_eps deps/dy) + (C1 P - C2 eps) eps/k,   P = nu_t (dU/dy)^2,
    with the diffusion terms by central differences on the given (uniform or not) y grid, and return
    the two residuals at the interior points, each divided by its largest source term (eps and
    C2 eps^2/k). Both vanish only if kappa^2 = sigma_eps sqrt(C_mu) (C2 - C1)."""
    y = np.asarray(y, dtype=float)
    p = keps_log_layer(y, u_tau, kappa, c_mu)
    k, eps, nu_t = p["k"], p["eps"], p["nu_t"]
    P = nu_t * p["dudy"] ** 2

    def diffusion(phi, sigma):
        ym = 0.5 * (y[1:] + y[:-1])
        num = 0.5 * (nu_t[1:] + nu_t[:-1])
        flux = num / sigma * np.diff(phi) / np.diff(y)
        return np.diff(flux) / np.diff(ym)

    rk = (diffusion(k, sigma_k) + (P - eps)[1:-1]) / eps[1:-1]
    re = (diffusion(eps, sigma_eps) + ((c1 * P - c2 * eps) * eps / k)[1:-1]) / (c2 * eps ** 2 / k)[1:-1]
    return rk, re


def kappa_from_keps(c_mu: float = C_MU, c1: float = C_EPS1, c2: float = C_EPS2,
                    sigma_eps: float = SIGMA_EPS) -> float:
    """The von Karman constant implied by the k-epsilon constants: kappa = sqrt(sigma_eps sqrt(C_mu) (C2 - C1))."""
    return float(np.sqrt(sigma_eps * np.sqrt(c_mu) * (c2 - c1)))


def c_mu_from_stress_ratio(ratio: float) -> float:
    """C_mu from the measured log-layer ratio -<u'v'>/k: since -<u'v'> = u_tau^2 and k = u_tau^2/sqrt(C_mu),
    the ratio is sqrt(C_mu), so C_mu = ratio^2 (a ratio of about 0.3 gives the standard 0.09)."""
    return float(ratio ** 2)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save turbulence_models.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # the engineering question: a flat cabin-air duct
    d = duct_flow(5.0, 0.05)
    print("flat cabin-air duct, gap 50 mm, bulk velocity 5 m/s, air at 20 C (mixing length, van Driest A+ = 26):")
    print(f"  Re_b = {d['re_b']:.0f}, Re_tau = {d['re_tau']:.1f}, u_tau = {d['u_tau']:.4f} m/s, "
          f"tau_w = {d['tau_w']:.4f} Pa, dp/dx = {d['dpdx']:.3f} Pa/m")
    print(f"  c_f = {d['cf']:.5f} (model) vs {d['cf_dean']:.5f} (Dean correlation): "
          f"{(d['cf'] / d['cf_dean'] - 1) * 100:+.1f} %")
    dn = duct_flow(5.0, 0.05, outer="nikuradse")
    print(f"  with Nikuradse's outer mixing length: Re_tau = {dn['re_tau']:.1f}, tau_w = {dn['tau_w']:.4f} Pa, "
          f"dp/dx = {dn['dpdx']:.3f} Pa/m, c_f = {dn['cf']:.5f} ({(dn['cf'] / dn['cf_dean'] - 1) * 100:+.1f} %)")
    tw_dean = d["cf_dean"] * 0.5 * RHO_AIR * 5.0 ** 2
    print(f"  Dean correlation: tau_w = {tw_dean:.4f} Pa, dp/dx = {-tw_dean / 0.025:.3f} Pa/m")
    rb, cf_log = log_law_friction(d["re_tau"])
    print(f"  log-law friction law at the same Re_tau: Re_b = {rb:.0f}, c_f = {cf_log:.5f}")
    print(f"  first-cell height for y+ = 1: {d['y1'] * 1e6:.1f} um; for y+ = 30 (wall functions): {d['y30'] * 1e3:.2f} mm")
    lam = 12.0 / d["re_b"]
    print(f"  (laminar flow at this Re_b would have c_f = 12/Re_b = {lam:.5f}, {d['cf'] / lam:.1f} times less)")

    print("\nRe_tau sweep (mixing length, prandtl outer):")
    print("  Re_tau     Re_b      c_f model   c_f Dean   diff %   kappa_fit  B_fit   U_b+   u+_cl")
    rows = []
    for rt in (180, 550, 1000, 2000, 5200):
        prof = mixing_length_profile(rt)
        reb, cf = friction(prof)
        kf, Bf = fit_log_law(prof["y_plus"], prof["u_plus"])
        rows.append((rt, reb, cf))
        print(f"  {rt:6d}  {reb:9.0f}   {cf:.5f}    {dean_cf(reb):.5f}   {(cf / dean_cf(reb) - 1) * 100:+5.1f}"
              f"     {kf:.4f}   {Bf:.3f}  {bulk_velocity(prof):.3f}  {prof['u_plus'][-1]:.3f}")
    prof = mixing_length_profile(1000, outer="nikuradse")
    reb, cf = friction(prof)
    print(f"  Nikuradse outer length at Re_tau = 1000: Re_b = {reb:.0f}, c_f = {cf:.5f} "
          f"({(cf / dean_cf(reb) - 1) * 100:+.1f} % against Dean), u+_cl = {prof['u_plus'][-1]:.3f}")

    print("\nk-epsilon log layer:")
    kk = kappa_from_keps()
    print(f"  kappa implied by the standard constants: sqrt(1.3 sqrt(0.09) (1.92 - 1.44)) = {kk:.4f}")
    y = np.linspace(0.01, 0.1, 201)
    rk, re = keps_residuals(y, 1.0, kappa=kk)
    print(f"  residuals with kappa = {kk:.4f}: k eqn {np.max(np.abs(rk)):.1e}, eps eqn {np.max(np.abs(re)):.1e}")
    rk, re = keps_residuals(y, 1.0, kappa=0.41)
    print(f"  residuals with kappa = 0.41:   k eqn {np.max(np.abs(rk)):.1e}, eps eqn {np.max(np.abs(re)):.3f}")
    print(f"  C_eps2 needed for kappa = 0.41 with the other constants: "
          f"{C_EPS1 + 0.41 ** 2 / (SIGMA_EPS * np.sqrt(C_MU)):.4f}; sigma_eps needed: "
          f"{0.41 ** 2 / (np.sqrt(C_MU) * (C_EPS2 - C_EPS1)):.4f}")
    print(f"  C_mu from a measured -<u'v'>/k = 0.30: {c_mu_from_stress_ratio(0.30):.3f}")

    prof = mixing_length_profile(1000)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    yp, up = prof["y_plus"], prof["u_plus"]
    ax[0].semilogx(yp[1:], up[1:], label="mixing length + van Driest (model)")
    yl = np.logspace(0, np.log10(1000), 100)
    ax[0].semilogx(yl, log_law(yl), "k--", lw=0.8, label="log law, kappa 0.41, B 5.0")
    ax[0].semilogx(yl[yl < 15], yl[yl < 15], "k:", lw=0.8, label="u+ = y+")
    ax[0].set_xlabel("y+")
    ax[0].set_ylabel("u+")
    ax[0].set_title("channel, Re_tau = 1000 (model, not data)")
    ax[0].legend(fontsize=8)
    rb = np.logspace(3.6, 5.6, 50)
    ax[1].loglog(rb, dean_cf(rb), "k--", label="Dean 0.073 Re_b^-1/4 (correlation)")
    r = np.array(rows)
    ax[1].loglog(r[:, 1], r[:, 2], "o", label="mixing-length model")
    rts = np.logspace(2, 3.8, 30)
    lf = np.array([log_law_friction(t) for t in rts])
    ax[1].loglog(lf[:, 0], lf[:, 1], ":", label="log-law friction law")
    ax[1].set_xlabel("Re_b = U_b 2h / nu")
    ax[1].set_ylabel("c_f")
    ax[1].set_title("channel skin friction")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "turbulence_models.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'turbulence_models.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
