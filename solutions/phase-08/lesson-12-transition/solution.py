"""P08.12 Transition: reference solution.

How far back can the boundary layer on the light aeroplane's wing stay laminar? This module builds
the two halves of the e^N method:

* linear stability theory: the Orr-Sommerfeld equation for small wavy disturbances
  v(y) exp(i alpha (x - c t)) on a parallel shear flow U(y),

      (U - c)(D^2 - alpha^2) v - U'' v = (D^2 - alpha^2)^2 v / (i alpha Re),

  solved for plane Poiseuille flow U = 1 - y^2 (y in [-1, 1], v = Dv = 0 at the walls) as a
  generalised matrix eigenvalue problem A v = c B v by Chebyshev collocation. A mode grows when
  Im(c) > 0. Its checks are published (labelled published): the critical Reynolds number 5772.22 at
  alpha = 1.02056 (Orszag 1971) and Orszag's eigenvalue c = 0.23752649 + 0.00373967 i at alpha = 1,
  Re = 10 000, with spectral (faster than any power of 1/N) convergence;

* the e^N envelope for a boundary layer: the amplification exponent N of the most amplified
  disturbance grows almost linearly with the momentum-thickness Reynolds number Re_theta beyond a
  critical value, and transition is predicted where N reaches N_crit (about 9 in quiet air). The
  growth rate dN/dRe_theta and the critical Re_theta are taken from Drela and Giles's (1987)
  *approximate fit* to Orr-Sommerfeld results for Falkner-Skan profiles, as functions of the shape
  factor H. The fit is supplied, not derived; it is labelled an approximate fit everywhere.

Units: SI (m, s) where dimensional; the stability problem is dimensionless.

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes transition.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-08/lesson-12-transition
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eig
from scipy.optimize import brentq, minimize_scalar

RE_CRIT_PUBLISHED, ALPHA_CRIT_PUBLISHED = 5772.22, 1.02056          # Orszag (1971), published
C_ORSZAG = 0.23752649 + 0.00373967j                                  # alpha = 1, Re = 10 000, published
BLASIUS = {"theta": 0.664115, "H": 2.591}                            # Blasius constants (published; P08.10)

# ISA constants (shared Phase 08 numbers)
T0, P0, R_AIR, G0, LAPSE = 288.15, 101325.0, 287.05287, 9.80665, -0.0065


# ---------------------------------------------------------------- given helpers

def isa_nu(h: float) -> float:
    """Kinematic viscosity (m^2/s) of the ISA troposphere at geopotential altitude h (m) (P08.01)."""
    T = T0 + LAPSE * h
    p = P0 * (T / T0) ** (-G0 / (R_AIR * LAPSE))
    return 1.458e-6 * T ** 1.5 / (T + 110.4) / (p / (R_AIR * T))


def drela_giles_dn_dretheta(H: np.ndarray) -> np.ndarray:
    """APPROXIMATE FIT (Drela and Giles 1987) to the slope of the e^N envelope of Falkner-Skan
    profiles, dN/dRe_theta = 0.01 sqrt([2.4 H - 3.7 + 2.5 tanh(1.5 H - 4.65)]^2 + 0.25).
    Supplied as a stability chart in formula form; not derived in this lesson."""
    H = np.asarray(H, dtype=float)
    return 0.01 * np.sqrt((2.4 * H - 3.7 + 2.5 * np.tanh(1.5 * H - 4.65)) ** 2 + 0.25)


def drela_giles_retheta0(H: np.ndarray) -> np.ndarray:
    """APPROXIMATE FIT (Drela and Giles 1987) to the critical momentum-thickness Reynolds number at
    which the envelope starts to grow: log10 Re_theta0 = (1.415/(H-1) - 0.489) tanh(20/(H-1) - 12.9)
    + 3.295/(H-1) + 0.44."""
    Hm = np.asarray(H, dtype=float) - 1.0
    return 10.0 ** ((1.415 / Hm - 0.489) * np.tanh(20.0 / Hm - 12.9) + 3.295 / Hm + 0.44)


def thwaites_wedge(m: float) -> tuple[float, float]:
    """Thwaites's method (P08.11) for a wedge flow Ue = C x^m: Re_theta^2 = 0.45 Re_x / (5m + 1)
    (Re_x = Ue x / nu) and H from Cebeci and Bradshaw's fit at lambda = 0.45 m / (5m + 1).
    Returns (Re_theta / sqrt(Re_x), H)."""
    lam = 0.45 * m / (5 * m + 1)
    H = 2.61 - 3.75 * lam + 5.24 * lam ** 2 if lam >= 0 else 2.088 + 0.0731 / (lam + 0.14)
    return float(np.sqrt(0.45 / (5 * m + 1))), float(H)


def mack_ncrit(Tu: float) -> float:
    """Mack's (1977) correlation between free-stream turbulence intensity Tu (a fraction, e.g. 0.001
    for 0.1 %) and the critical amplification exponent: N_crit = -8.43 - 2.4 ln(Tu). A correlation,
    valid roughly for 0.0007 < Tu < 0.03."""
    return -8.43 - 2.4 * np.log(Tu)


def observed_decay(errs: np.ndarray) -> np.ndarray:
    """Ratio of successive errors e_k / e_{k+1} (constant or growing for spectral convergence)."""
    e = np.asarray(errs, dtype=float)
    return e[:-1] / e[1:]


# ---------------------------------------------------------------- Chebyshev collocation

def cheb(N: int) -> tuple[np.ndarray, np.ndarray]:
    """Chebyshev differentiation matrix on the N + 1 Gauss-Lobatto points x_j = cos(pi j / N),
    j = 0..N (so x runs from 1 down to -1), after Trefethen, Spectral Methods in MATLAB:
    D_ij = (c_i / c_j) (-1)^(i+j) / (x_i - x_j) off the diagonal, c_0 = c_N = 2, c_j = 1 otherwise,
    and the diagonal from the negative row sums (so D @ ones = 0 exactly). Returns (D, x)."""
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1)
    c[0] = c[-1] = 2.0
    c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    dX = X - X.T
    D = np.outer(c, 1.0 / c) / (dX + np.eye(N + 1))
    D -= np.diag(D.sum(axis=1))
    return D, x


def os_matrices(N: int, alpha: float, Re: float) -> tuple[np.ndarray, np.ndarray]:
    """Orr-Sommerfeld matrices for plane Poiseuille flow U = 1 - y^2 on the N - 1 interior
    Chebyshev points, A v = c B v with

        A = U (D2 - alpha^2) - U'' - (D4 - 2 alpha^2 D2 + alpha^4) / (i alpha Re),   B = D2 - alpha^2,

    where v = 0 at the walls is imposed by deleting the boundary rows and columns, and Dv = 0 by
    building the fourth derivative on v = (1 - y^2) g (Trefethen's trick):
    D4 = (diag(1 - y^2) D^4 - 8 diag(y) D^3 - 12 D^2) diag(0, 1/(1 - y^2), 0), interior block.
    Returns (A, B), complex arrays of shape (N - 1, N - 1)."""
    D, y = cheb(N)
    D2full = D @ D
    S = np.zeros(N + 1)
    S[1:-1] = 1.0 / (1.0 - y[1:-1] ** 2)
    D3full = D2full @ D
    D4full = (np.diag(1 - y ** 2) @ D3full @ D - 8 * np.diag(y) @ D3full - 12 * D2full) @ np.diag(S)
    D2 = D2full[1:-1, 1:-1]
    D4 = D4full[1:-1, 1:-1]
    yi = y[1:-1]
    I = np.eye(N - 1)
    U, Upp = 1.0 - yi ** 2, -2.0
    L2 = D2 - alpha ** 2 * I
    A = np.diag(U) @ L2 - Upp * I - (D4 - 2 * alpha ** 2 * D2 + alpha ** 4 * I) / (1j * alpha * Re)
    return A.astype(complex), L2.astype(complex)


def os_eigenvalues(N: int, alpha: float, Re: float) -> np.ndarray:
    """All finite eigenvalues c of the Orr-Sommerfeld problem (``os_matrices``), sorted by decreasing
    imaginary part (the least stable first)."""
    A, B = os_matrices(N, alpha, Re)
    c = eig(A, B, right=False)
    c = c[np.isfinite(c)]
    return c[np.argsort(-c.imag)]


def least_stable(alpha: float, Re: float, N: int = 60) -> complex:
    """The eigenvalue c with the largest imaginary part: the mode that decays slowest or grows."""
    return complex(os_eigenvalues(N, alpha, Re)[0])


def neutral_reynolds(alpha: float, N: int = 60, Re_lo: float = 4000.0, Re_hi: float = 10000.0) -> float:
    """Reynolds number on the lower branch of the neutral curve at wavenumber alpha: the root of
    Im c(alpha, Re) = 0 for the least stable mode, by Brent's method on [Re_lo, Re_hi] (relative
    tolerance 1e-12). Raises ValueError if the bracket does not contain a sign change."""
    f = lambda R: least_stable(alpha, R, N).imag
    return float(brentq(f, Re_lo, Re_hi, rtol=1e-12))


def critical_point(N: int = 60, bracket: tuple[float, float] = (0.95, 1.08)) -> tuple[float, float]:
    """The critical point (Re_c, alpha_c): the minimum over alpha of the neutral Reynolds number, by a
    bounded scalar minimisation (xatol 1e-6) on ``bracket``. Returns (Re_c, alpha_c)."""
    res = minimize_scalar(lambda a: neutral_reynolds(a, N), bounds=bracket, method="bounded",
                          options={"xatol": 1e-6})
    return float(res.fun), float(res.x)


# ---------------------------------------------------------------- the e^N envelope

def en_envelope(Re_theta: np.ndarray, H: float) -> np.ndarray:
    """Amplification exponent N of the most amplified disturbance at momentum-thickness Reynolds
    numbers Re_theta for a boundary layer of constant shape factor H, from the approximate fit:
    N = dN/dRe_theta (Re_theta - Re_theta0) beyond Re_theta0, and 0 before it."""
    Rt = np.asarray(Re_theta, dtype=float)
    return np.where(Rt > drela_giles_retheta0(H), drela_giles_dn_dretheta(H) * (Rt - drela_giles_retheta0(H)), 0.0)


def transition_reynolds(H: float, retheta_per_rootrex: float, N_crit: float = 9.0) -> dict[str, float]:
    """Instability and transition points of a similar boundary layer (constant H, Re_theta =
    k sqrt(Re_x) with k = ``retheta_per_rootrex``; k = 0.664 for Blasius): the critical
    Re_theta0 and Re_x0 = (Re_theta0 / k)^2 where the envelope starts, and the transition
    Re_theta_tr = Re_theta0 + N_crit / (dN/dRe_theta) and Re_x_tr = (Re_theta_tr / k)^2.
    Returns a dict with keys Re_theta0, Re_x0, Re_theta_tr, Re_x_tr."""
    R0 = float(drela_giles_retheta0(H))
    Rtr = R0 + N_crit / float(drela_giles_dn_dretheta(H))
    k = retheta_per_rootrex
    return {"Re_theta0": R0, "Re_x0": (R0 / k) ** 2, "Re_theta_tr": Rtr, "Re_x_tr": (Rtr / k) ** 2}


def flat_plate_transition(U: float, nu: float, N_crit: float = 9.0) -> dict[str, float]:
    """Instability point x0 and e^N transition point x_tr (m) on a flat plate at zero pressure
    gradient in a stream U (m/s) of kinematic viscosity nu (m^2/s), with the Blasius layer
    (H = 2.591, Re_theta = 0.664 sqrt(Re_x)). Returns the dict of ``transition_reynolds`` plus
    x0 and x_tr."""
    r = transition_reynolds(BLASIUS["H"], BLASIUS["theta"], N_crit)
    r["x0"] = r["Re_x0"] * nu / U
    r["x_tr"] = r["Re_x_tr"] * nu / U
    return r


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save transition.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("Orr-Sommerfeld, plane Poiseuille flow, alpha = 1, Re = 10 000 (published c = 0.23752649 + 0.00373967i):")
    Ns = [20, 30, 40, 50, 60, 80]
    errs = []
    for N in Ns:
        c = least_stable(1.0, 1e4, N)
        errs.append(abs(c - C_ORSZAG))
        print(f"  N = {N:3d}: c = {c.real:.8f} + {c.imag:.8f}i, |c - published| = {errs[-1]:.2e}")

    Re_c, a_c = critical_point(60)
    print(f"critical point (N = 60): Re_c = {Re_c:.2f}, alpha_c = {a_c:.5f}; published {RE_CRIT_PUBLISHED}, {ALPHA_CRIT_PUBLISHED}")
    print(f"  relative error in Re_c {(Re_c / RE_CRIT_PUBLISHED - 1):.1e}")
    for N in (30, 40, 50):
        print(f"  N = {N}: neutral Re at alpha = 1.02056: {neutral_reynolds(1.02056, N):.3f}")

    # the e^N envelope for the light aeroplane's wing (flat plate, zero pressure gradient)
    nu, U, chord = isa_nu(3000.0), 60.0, 1.5
    fp = flat_plate_transition(U, nu, 9.0)
    print(f"\nflat plate, Blasius H = {BLASIUS['H']} (approximate-fit stability chart): dN/dRe_theta = "
          f"{float(drela_giles_dn_dretheta(BLASIUS['H'])):.5f}, Re_theta0 = {fp['Re_theta0']:.1f}")
    print(f"  instability: Re_x0 = {fp['Re_x0']:.4g}; transition (N = 9): Re_theta = {fp['Re_theta_tr']:.0f}, "
          f"Re_x = {fp['Re_x_tr']:.4g}")
    print(f"  wing at 3000 m, {U} m/s (nu = {nu:.4e}): instability at x = {fp['x0'] * 1e3:.1f} mm, transition at "
          f"x = {fp['x_tr']:.3f} m = {fp['x_tr'] / chord * 100:.0f} % of the {chord} m chord")
    print("  sensitivity to N_crit (Mack: Tu -> N_crit):")
    for Tu in (0.0007, 0.001, 0.003, 0.01, 0.02):
        Nc = mack_ncrit(Tu)
        r = flat_plate_transition(U, nu, Nc)
        print(f"    Tu = {Tu * 100:4.2f} %: N_crit = {Nc:5.2f}, Re_x_tr = {r['Re_x_tr']:.3g}, x_tr = {r['x_tr']:.3f} m "
              f"({min(r['x_tr'] / chord, 9.99) * 100:.0f} % chord)")
    print("  pressure gradient (wedge flows Ue ~ x^m, Thwaites H; N_crit = 9):")
    for m in (-0.05, -0.02, 0.0, 0.05, 0.1, 0.2):
        k, H = thwaites_wedge(m)
        r = transition_reynolds(H, k, 9.0)
        print(f"    m = {m:+.2f}: H = {H:.3f}, Re_theta0 = {r['Re_theta0']:.0f}, Re_x_tr = {r['Re_x_tr']:.3g}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].semilogy(Ns, errs, "o-")
    ax[0].set_xlabel("Chebyshev points N")
    ax[0].set_ylabel("|c − published|")
    ax[0].set_title("Orr–Sommerfeld: spectral convergence (α = 1, Re = 10⁴)")
    x = np.linspace(0.0, chord, 400)
    Rt = BLASIUS["theta"] * np.sqrt(U * x / nu)
    ax[1].plot(x / chord, en_envelope(Rt, BLASIUS["H"]), label="flat plate (approximate-fit chart)")
    for m in (-0.02, 0.05):
        k, H = thwaites_wedge(m)
        ax[1].plot(x / chord, en_envelope(k * np.sqrt(U * x / nu), H), label=f"wedge flow m = {m:+.2f}")
    ax[1].axhline(9.0, color="r", ls="--", lw=0.8, label="N = 9")
    ax[1].set_xlabel("x / chord (U = 60 m/s, 3000 m)")
    ax[1].set_ylabel("N = ln(A/A₀)")
    ax[1].set_title("e^N envelope on the wing")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "transition.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'transition.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
