"""P08.14 Turbulent Boundary Layers: starter code.

How much skin-friction drag does the light aeroplane's wing carry when its boundary layer is
turbulent, compared with the all-laminar layer of P08.10 and with transition at the e^N point of
P08.12? This module builds the pieces:

* the law of the wall: the viscous sublayer u+ = y+, the log law u+ = ln(y+)/kappa + B, and
  Spalding's single composite formula y+(u+) that joins them through the buffer layer, inverted for
  u+(y+) by Newton's method; Coles's wake function for the outer part;
* a Clauser fit, which recovers the friction velocity u_tau from a measured (here SYNTHETIC) mean
  profile by fitting the log law, and therefore the wall shear without measuring it;
* the 1/7th power law and the momentum integral (with Blasius's pipe-flow wall-shear law), giving
  delta/x = 0.37 Re_x^(-1/5) and c_f = 0.0577 Re_x^(-1/5);
* flat-plate skin-friction correlations (labelled correlations, not measurements), local and
  integrated, and the mixed laminar-turbulent plate two ways (Prandtl's subtraction and a
  momentum-thickness match at transition);
* an independent friction law built from the log-wake profile and the momentum integral alone, with
  no empirical skin-friction formula, to check the correlations against.

Constants: kappa = 0.41, B = 5.0 (classical values; other fits such as kappa = 0.384, B = 4.17 are
discussed in the lesson), Coles wake parameter Pi = 0.55 for a flat plate at high Reynolds number.

Units: SI (m, s, kg, N).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes turbulent_bl.png)
Test: pytest curriculum/phase-08/lesson-14-turbulent-boundary-layers
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

KAPPA, B_LOG, PI_COLES = 0.41, 5.0, 0.55
BLASIUS_FPP0 = 0.332057          # published; P08.10

# ISA constants (shared Phase 08 numbers)
T0, P0, R_AIR, G0, LAPSE = 288.15, 101325.0, 287.05287, 9.80665, -0.0065


# ---------------------------------------------------------------- given helpers

def isa_troposphere(h: float) -> tuple[float, float, float, float]:
    """ISA temperature (K), pressure (Pa), density (kg/m^3) and Sutherland viscosity (Pa s) below
    11 km (P08.01)."""
    T = T0 + LAPSE * h
    p = P0 * (T / T0) ** (-G0 / (R_AIR * LAPSE))
    return T, p, p / (R_AIR * T), 1.458e-6 * T ** 1.5 / (T + 110.4)


def coles_wake(eta: np.ndarray, Pi: float = PI_COLES, kappa: float = KAPPA) -> np.ndarray:
    """Coles's wake function added to the law of the wall in the outer layer:
    (2 Pi / kappa) sin^2(pi eta / 2), eta = y/delta (clipped to [0, 1])."""
    e = np.clip(np.asarray(eta, dtype=float), 0.0, 1.0)
    return 2 * Pi / kappa * np.sin(np.pi * e / 2) ** 2


def laminar_CF(Re: np.ndarray) -> np.ndarray:
    """Blasius flat-plate friction coefficient, both laminar results of P08.10: integrated
    C_F = 4 f''(0) / sqrt(Re_L) = 1.328 / sqrt(Re_L)."""
    return 4 * BLASIUS_FPP0 / np.sqrt(np.asarray(Re, dtype=float))


# ---------------------------------------------------------------- the law of the wall

def log_law(yplus: np.ndarray, kappa: float = KAPPA, B: float = B_LOG) -> np.ndarray:
    """Log law of the wall: u+ = ln(y+)/kappa + B."""
    raise NotImplementedError


def spalding_yplus(uplus: np.ndarray, kappa: float = KAPPA, B: float = B_LOG) -> np.ndarray:
    """Spalding's (1961) composite law of the wall, explicit in u+:
    y+ = u+ + exp(-kappa B) [exp(kappa u+) - 1 - kappa u+ - (kappa u+)^2/2 - (kappa u+)^3/6].
    It tends to y+ = u+ at the wall and to the log law far from it."""
    raise NotImplementedError


def spalding_uplus(yplus: np.ndarray, kappa: float = KAPPA, B: float = B_LOG, tol: float = 1e-12) -> np.ndarray:
    """Invert Spalding's law for u+(y+) by Newton's method (vectorised), starting from the smaller of
    y+ and the log law. y+(u+) is increasing and convex, so after the first step every iterate is at
    or above the root and Newton converges monotonically from above (the guess is already above the
    root in the sublayer and buffer layer). Stops when every |step| < tol (u+ units). y+ = 0 gives u+ = 0."""
    raise NotImplementedError


def synthetic_profile(y: np.ndarray, u_tau: float, nu: float, delta: float, Pi: float = PI_COLES,
                      noise: float = 0.01, seed: int = 0) -> np.ndarray:
    """SYNTHETIC mean-velocity profile (not a measurement): Spalding's law plus Coles's wake,
    U = u_tau [u+(y u_tau/nu) + w(y/delta)], with multiplicative Gaussian noise of relative standard
    deviation ``noise`` (numpy.random.default_rng(seed)) standing in for probe scatter. Points above
    delta get the edge value."""
    yp = np.asarray(y, dtype=float) * u_tau / nu
    U = u_tau * (spalding_uplus(np.minimum(yp, delta * u_tau / nu)) + coles_wake(np.asarray(y) / delta, Pi))
    rng = np.random.default_rng(seed)
    return U * (1 + noise * rng.standard_normal(U.shape))


def clauser_fit(y: np.ndarray, U: np.ndarray, nu: float, delta: float, kappa: float = KAPPA, B: float = B_LOG,
                yplus_min: float = 30.0, eta_max: float = 0.2) -> tuple[float, int]:
    """Clauser fit: the friction velocity u_tau that makes the measured points of the log region
    (y+ >= yplus_min and y/delta <= eta_max) fit U/u_tau = ln(y u_tau/nu)/kappa + B in least squares.
    Because the log region depends on u_tau, the point selection and a bounded scalar minimisation
    of sum [U - u_tau (ln(y u_tau/nu)/kappa + B)]^2 over u_tau in [1e-4, 0.2] max(U) are repeated
    (starting from u_tau = 0.04 max(U)) until the selection stops changing.
    Returns (u_tau, number of points used). Raises ValueError if fewer than 3 points qualify."""
    raise NotImplementedError


# ---------------------------------------------------------------- power law and correlations

def power_law_bl(x: np.ndarray, U: float, nu: float) -> dict[str, np.ndarray]:
    """Turbulent flat-plate layer from the 1/7th power law u/U = (y/delta)^(1/7), turbulent from x = 0,
    and Blasius's pipe-flow wall shear tau_w = 0.0225 rho U^2 (nu/(U delta))^(1/4) in the momentum
    integral d(theta)/dx = c_f/2 with theta = 7 delta / 72: delta/x = 0.3707 Re_x^(-1/5),
    c_f = 0.0577 Re_x^(-1/5) (exact coefficients computed, not rounded), delta* = delta/8,
    theta = 7 delta/72, H = 9/7. Returns a dict with keys Re_x, delta, delta_star, theta, H, cf."""
    raise NotImplementedError


def cf_local(Re_x: np.ndarray, kind: str = "white") -> np.ndarray:
    """Local turbulent skin-friction CORRELATIONS for a smooth flat plate:
    'power'           0.0592 Re_x^(-1/5)                 (Prandtl's 1/7th-law fit; Re_x < 1e7)
    'schultz-grunow'  0.370 (log10 Re_x)^(-2.584)
    'white'           0.455 / ln^2(0.06 Re_x)."""
    raise NotImplementedError


def CF_plate(Re_L: np.ndarray, kind: str = "prandtl-schlichting") -> np.ndarray:
    """Integrated (average) turbulent skin-friction CORRELATIONS, turbulent from the leading edge:
    'power'                0.074 Re_L^(-1/5)
    'prandtl-schlichting'  0.455 (log10 Re_L)^(-2.58)
    'karman-schoenherr'    0.242/sqrt(C_F) = log10(Re_L C_F), solved by Newton's method
    'ittc'                 0.075 / (log10 Re_L - 2)^2   (the ITTC-1957 ship model-correlation line)."""
    raise NotImplementedError


def integrate_cf(Re_L: float, kind: str = "white", Re_start: float = 1e3, n: int = 4001) -> float:
    """Average of a local correlation over a plate turbulent from (nearly) the leading edge,
    C_F = (1/Re_L) int c_f dRe_x, by the trapezoid rule on n log-spaced points from Re_start to Re_L.
    (The piece below Re_start is dropped: the correlations are meaningless there and it is ~0.1 %.)"""
    raise NotImplementedError


def mixed_plate_prandtl(Re_L: float, Re_tr: float, kind: str = "prandtl-schlichting") -> float:
    """Mixed laminar-turbulent plate by Prandtl's subtraction: the turbulent plate of length L, minus
    the turbulent drag it would have had up to x_tr, plus the laminar drag up to x_tr:
    C_F = C_F,turb(Re_L) - (Re_tr/Re_L) [C_F,turb(Re_tr) - C_F,lam(Re_tr)]."""
    raise NotImplementedError


def mixed_plate_momentum(Re_L: float, Re_tr: float) -> tuple[float, float]:
    """Mixed plate by matching momentum thickness at transition: laminar (Blasius) theta up to x_tr,
    theta_tr = 2 f''(0) x_tr / sqrt(Re_tr); beyond it a 1/7th-law turbulent layer grown from a virtual
    origin x0 chosen so that its theta equals theta_tr at x_tr (theta = 0.0360 (x - x0) Re_{x-x0}^(-1/5)
    from ``power_law_bl``). The drag of one side up to L is rho U^2 theta(L) (momentum integral), so
    C_F = 2 theta(L)/L. Works in Reynolds numbers (lengths in units of nu/U).
    Returns (C_F, Re_x0) with Re_x0 the virtual origin's Reynolds number (may be negative)."""
    raise NotImplementedError


# ---------------------------------------------------------------- an independent friction law

def log_wake_state(delta_plus: float, Pi: float = PI_COLES, kappa: float = KAPPA, B: float = B_LOG,
                   n: int = 4000) -> dict[str, float]:
    """Integral properties of the composite profile u+ = Spalding(y+) + Coles wake(y/delta) for a
    layer of thickness delta+ = delta u_tau / nu: the edge velocity Ue+ = u+(delta+), the friction
    coefficient c_f = 2 / Ue+^2, Re_theta = Ue+ theta+ and Re_delta* = Ue+ delta*+ (theta+ and
    delta*+ by the trapezoid rule on n points log-spaced in y+ from 1e-3 to delta+, plus the
    sublayer piece below), and H. Returns a dict with keys Ue_plus, cf, Re_theta, Re_dstar, H."""
    raise NotImplementedError


def log_wake_plate(Re_x: np.ndarray, Pi: float = PI_COLES, kappa: float = KAPPA, B: float = B_LOG) -> np.ndarray:
    """Local c_f(Re_x) of a flat plate from the log-wake profile and the momentum integral alone:
    along a sweep of delta+ (log-spaced 200 to 2e5), dRe_x/dRe_theta = 2/c_f, so
    Re_x = int (2/c_f) dRe_theta (trapezoid rule, starting from Re_x = Re_theta / 0.0360^... chosen as
    the 1/7th-law value at the first point), then c_f is interpolated in log Re_x. No skin-friction
    correlation is used."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save turbulent_bl.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    h, V, chord, semi = 3000.0, 60.0, 1.5, 5.0
    T, p, rho, mu = isa_troposphere(h)
    nu = mu / rho
    q = 0.5 * rho * V ** 2
    S_wet = 2 * 2 * semi * chord
    ReL = V * chord / nu
    print(f"wing at {h:.0f} m, {V} m/s: rho {rho:.4f} kg/m^3, nu {nu:.4e} m^2/s, q {q:.1f} Pa, Re_c {ReL:.4g}, "
          f"wetted area {S_wet:.0f} m^2")

    print("\nintegrated skin friction C_F at Re_c (correlations, turbulent from the leading edge):")
    for kind in ("power", "prandtl-schlichting", "karman-schoenherr", "ittc"):
        print(f"  {kind:20s} C_F = {CF_plate(ReL, kind):.5f}, drag {q * CF_plate(ReL, kind) * S_wet:.1f} N")
    for kind in ("power", "schultz-grunow", "white"):
        print(f"  integral of local '{kind}': C_F = {integrate_cf(ReL, kind):.5f}; local c_f at the trailing edge "
              f"{float(cf_local(ReL, kind)):.5f}")
    CFl = float(laminar_CF(ReL))
    print(f"  all-laminar (P08.10): C_F = {CFl:.5f}, drag {q * CFl * S_wet:.1f} N")
    for Rtr, label in ((5e5, "Re_tr = 5e5 (a common textbook default)"), (2.785e6, "e^9 point of P08.12 (Re_tr = 2.785e6)")):
        cfp = mixed_plate_prandtl(ReL, Rtr)
        cfm, rx0 = mixed_plate_momentum(ReL, Rtr)
        print(f"  mixed plate, {label}: x_tr = {Rtr * nu / V:.3f} m; C_F Prandtl {cfp:.5f} ({q * cfp * S_wet:.1f} N), "
              f"momentum match {cfm:.5f} ({q * cfm * S_wet:.1f} N), virtual origin {rx0 * nu / V * 1e3:.0f} mm")

    pl = power_law_bl(chord, V, nu)
    cfw = float(cf_local(ReL, "white"))
    ut = V * np.sqrt(cfw / 2)
    print(f"\nturbulent trailing edge (from the leading edge): delta {float(pl['delta']) * 1e3:.1f} mm (1/7th law), "
          f"theta {float(pl['theta']) * 1e3:.2f} mm, c_f {cfw:.5f} (White), u_tau {ut:.3f} m/s, "
          f"tau_w {rho * ut ** 2:.2f} Pa, viscous length nu/u_tau {nu / ut * 1e6:.2f} um, delta+ {float(pl['delta']) * ut / nu:.0f}")
    print(f"  first CFD cell for y+ = 1: {nu / ut * 1e6:.1f} um; log layer starts (y+ = 30) at {30 * nu / ut * 1e3:.3f} mm")

    # Clauser fit on a SYNTHETIC profile
    delta = float(pl["delta"])
    y = np.geomspace(2 * nu / ut, delta, 40)
    Uy = synthetic_profile(y, ut, nu, delta, noise=0.01, seed=14)
    ut_fit, npts = clauser_fit(y, Uy, nu, delta)
    print(f"\nClauser fit to a SYNTHETIC trailing-edge profile (Spalding + Coles wake, Pi = {PI_COLES}, 1 % noise, seed 14):")
    print(f"  u_tau used {ut:.4f} m/s, recovered {ut_fit:.4f} m/s ({(ut_fit / ut - 1) * 100:+.2f} %) from {npts} points")
    ut_alt, _ = clauser_fit(y, Uy, nu, delta, kappa=0.384, B=4.17)
    print(f"  with kappa = 0.384, B = 4.17 instead: {ut_alt:.4f} m/s ({(ut_alt / ut - 1) * 100:+.2f} %)")

    # independent friction law
    Rx = np.array([1e6, ReL, 1e7, 1e8])
    cflw = log_wake_plate(Rx)
    print("\nlog-wake profile + momentum integral (no correlation) against correlations, local c_f:")
    for R, c in zip(Rx, cflw):
        print(f"  Re_x {R:.3g}: log-wake {c:.5f}, White {float(cf_local(R, 'white')):.5f}, "
              f"Schultz-Grunow {float(cf_local(R, 'schultz-grunow')):.5f}, power {float(cf_local(R, 'power')):.5f}")
    st = log_wake_state(1e4)
    print(f"  at delta+ = 1e4: Re_theta {st['Re_theta']:.0f}, H {st['H']:.3f}, c_f {st['cf']:.5f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    yp = np.geomspace(0.1, 1e4, 400)
    ax[0].semilogx(yp, spalding_uplus(yp), label="Spalding composite")
    ax[0].semilogx(yp[yp < 15], yp[yp < 15], "--", label="u⁺ = y⁺")
    ax[0].semilogx(yp[yp > 5], log_law(yp[yp > 5]), "--", label="log law κ = 0.41, B = 5.0")
    ax[0].semilogx(y * ut / nu, Uy / ut, "k.", ms=4, label="SYNTHETIC profile (u_τ known)")
    ax[0].set_xlabel("y⁺")
    ax[0].set_ylabel("u⁺")
    ax[0].set_title("law of the wall (synthetic data)")
    ax[0].legend(fontsize=8)
    R = np.geomspace(1e5, 1e8, 200)
    for kind in ("power", "schultz-grunow", "white"):
        ax[1].loglog(R, cf_local(R, kind), label=f"local c_f, {kind} (correlation)")
    ax[1].loglog(R[R > 3e5], log_wake_plate(R[R > 3e5]), "k--", label="log-wake + momentum integral")
    ax[1].loglog(R, 2 * BLASIUS_FPP0 / np.sqrt(R), ":", label="laminar (Blasius)")
    ax[1].set_xlabel("Re_x")
    ax[1].set_ylabel("c_f")
    ax[1].set_title("flat-plate skin friction")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out_dir / "turbulent_bl.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'turbulent_bl.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
