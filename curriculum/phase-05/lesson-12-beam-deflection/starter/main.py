"""P05.12 Beam Deflection: starter code.

How far does a wing tip rise at limit load, and does the deflection stay acceptable? This module
answers with three independent routes to the deflection of an Euler-Bernoulli beam:

* closed forms from integrating EI w'' = M by hand for the standard cases (cantilever and simply
  supported beams under point loads, uniform loads and end moments), combined by superposition,
  which also solves a statically indeterminate propped cantilever;
* numerical double integration of the curvature M/EI with the boundary conditions, for any load and
  any spanwise stiffness EI(x) (a tapered spar);
* a finite-difference solution of (EI w'')'' = q written as two second-order equations,
  m = EI w'' and m'' = q, with second-order boundary rows: it needs no bending-moment diagram, so it
  checks the other two and also handles indeterminate supports.

Sign convention: x along the beam from the left end (the root of a cantilever), w is deflection
(positive up), q is distributed load (positive up, N/m), point loads P positive up, the bending
moment M is positive when it compresses the top fibres (sagging). Then EI w'' = M and M'' = q.
Lift on a wing (cantilever) gives M > 0 and an upward tip deflection.

Units: SI (N, m, Pa, N m^2).

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes wing_deflection.png)
Test: pytest curriculum/phase-05/lesson-12-beam-deflection
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665
E_7075 = 71.7e9                       # Pa, typical Young's modulus of 7075-T6 (ASM Aerospace Specification Metals data sheet)
I_SPAR_ROOT = 9.704e-6                # m^4, the P05.10 root spar section (180 mm I-section)


# ---------------------------------------------------------------- closed forms

def cantilever_point(x: np.ndarray, L: float, a: float, P: float, EI: float) -> np.ndarray:
    """Cantilever clamped at x = 0, point load P (up) at x = a:
    w = P x^2 (3a - x) / (6 EI) for x <= a, and P a^2 (3x - a) / (6 EI) beyond (straight line)."""
    raise NotImplementedError


def cantilever_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Cantilever clamped at x = 0 under uniform load q (up): w = q x^2 (6L^2 - 4Lx + x^2) / (24 EI);
    tip deflection q L^4 / (8 EI)."""
    raise NotImplementedError


def cantilever_end_moment(x: np.ndarray, L: float, M0: float, EI: float) -> np.ndarray:
    """Cantilever with a sagging-positive moment M0 applied at the free end: constant curvature,
    w = M0 x^2 / (2 EI)."""
    raise NotImplementedError


def simply_supported_point(x: np.ndarray, L: float, a: float, P: float, EI: float) -> np.ndarray:
    """Simply supported beam (pins at 0 and L), point load P (up) at x = a, b = L - a:
    w = P b x (L^2 - b^2 - x^2) / (6 L EI) for x <= a, and the mirror image beyond."""
    raise NotImplementedError


def simply_supported_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Simply supported beam under uniform load q (up): w = q x (L^3 - 2 L x^2 + x^3) / (24 EI);
    midspan deflection 5 q L^4 / (384 EI)."""
    raise NotImplementedError


def propped_cantilever_reaction(L: float, q: float) -> float:
    """Prop force at x = L of a cantilever (clamped at 0) under uniform q, by superposition: the prop
    must cancel the free cantilever's tip deflection, R L^3 / (3 EI) = q L^4 / (8 EI), so R = 3 q L / 8
    (acting against the load; returned as a positive magnitude)."""
    raise NotImplementedError


def propped_cantilever_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Deflection of the propped cantilever (clamped at 0, pinned at L) under uniform q, by
    superposition of ``cantilever_udl`` and a tip point load -R: w = q x^2 (3L^2 - 5Lx + 2x^2) / (48 EI)."""
    raise NotImplementedError


# ---------------------------------------------------------------- bending moment and integration

def _cumtrapz(f: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(x))])


def bending_moment(x: np.ndarray, q: np.ndarray, point_loads: list[tuple[float, float]] = (),
                   support: str = "cantilever") -> np.ndarray:
    """Bending moment of a statically determinate beam from M'' = q (plus point loads, which add a
    kink P (x - a)_+), with the two moment boundary conditions of the support:
    'cantilever' (clamped at x[0], free at x[-1]: M = M' = 0 at the tip) or 'simply_supported'
    (M = 0 at both ends). ``q`` holds the distributed load at the points x; integration is by the
    trapezoid rule."""
    raise NotImplementedError


def deflection_by_integration(x: np.ndarray, M: np.ndarray, EI: np.ndarray | float,
                              support: str = "cantilever") -> tuple[np.ndarray, np.ndarray]:
    """Slope and deflection from the curvature kappa = M / EI by two cumulative trapezoid integrals,
    with the support's kinematic conditions: 'cantilever' (w = w' = 0 at x[0]) or 'simply_supported'
    (w = 0 at both ends, the slope constant fixed by the second condition). EI may vary along x.
    Returns (theta, w)."""
    raise NotImplementedError


# ---------------------------------------------------------------- independent check: finite differences

def fd_beam(L: float, EI: float | Callable[[np.ndarray], np.ndarray], q: Callable[[np.ndarray], np.ndarray],
            n: int, bc: str = "clamped-free", point_loads: list[tuple[float, float]] = ()) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Finite-difference solution of (EI w'')'' = q on n equal intervals, as the coupled pair
    m'' = q, w'' = m / EI (m is the bending moment), with central differences at interior nodes.

    Boundary rows: w = 0 and m = 0 are imposed directly. A slope condition uses the last interval
    plus the curvature there, w'(0) ~ (w1 - w0)/h - (h/2) m0/EI0 (and its mirror at x = L); a shear
    condition uses m'(0) ~ (m1 - m0)/h - (h/2) q0 (and its mirror). Both are second-order accurate
    and need no ghost nodes or wide stencils, so a point load near an end does not spoil them.

    ``bc`` is one of 'clamped-free' (cantilever), 'pinned-pinned', 'clamped-clamped' or
    'clamped-pinned' (propped cantilever; indeterminate). ``point_loads`` (a, P) are snapped to the
    nearest node: at an interior node they enter as P/h in m'' = q, which makes m exact at the
    nodes; at a free end they enter the shear condition (m' = -P at x = L, +P at x = 0); at a
    supported end they go straight into the reaction and do not deflect the beam.
    Returns (x, w, m) with x the n + 1 nodes."""
    raise NotImplementedError


# ---------------------------------------------------------------- the wing

def elliptic_lift(y: np.ndarray, L_semi: float, s: float) -> np.ndarray:
    """Elliptic lift per unit span (N/m) on a semi-span s carrying L_semi: 4 L_semi/(pi s) sqrt(1 - (y/s)^2)."""
    u = np.clip(np.asarray(y, dtype=float) / s, 0.0, 1.0)
    return 4.0 * L_semi / (np.pi * s) * np.sqrt(1.0 - u ** 2)


def elliptic_tip(L_semi: float, s: float, EI: float) -> tuple[float, float]:
    """Tip slope and tip deflection of a uniform cantilever under elliptic lift, by superposing the
    tip response of a point load at eta, eta^2/(2EI) and eta^2 (3s - eta)/(6EI), over the load
    (int t^2 sqrt(1-t^2) = pi/16, int t^3 sqrt(1-t^2) = 2/15 on [0, 1]):

        theta_tip = L s^2 / (8 EI),   w_tip = L s^3 / EI * (1/8 - 4/(45 pi)).

    Compare the uniform load: L s^2 / (6 EI) and L s^3 / (8 EI)."""
    raise NotImplementedError


def tapered_EI(EI_root: float, taper: float, s: float) -> Callable[[np.ndarray], np.ndarray]:
    """Spanwise stiffness falling linearly from EI_root at the root to (1 - taper) EI_root at the tip."""
    return lambda y: EI_root * (1.0 - taper * np.asarray(y, dtype=float) / s)


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed convergence orders log(e_k / e_{k+1}) / log(ratio) of errors on grids refined by ``ratio``."""
    e = np.asarray(errors, dtype=float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save wing_deflection.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    s, L = 5.0, 3.8 * 750.0 * G0 / 2.0
    EI = E_7075 * I_SPAR_ROOT
    print(f"spar EI = {EI:.4g} N m^2 (7075-T6, I = {I_SPAR_ROOT * 1e12:.4g} mm^4), semi-span {s} m, lift {L:.1f} N")

    # closed forms for the bounding loads
    w_udl = float(cantilever_udl(s, s, L / s, EI))
    w_pt = float(cantilever_point(s, s, 4 * s / (3 * np.pi), L, EI))
    print(f"uniform lift: tip {w_udl:.4f} m; all lift at the elliptic centre of lift 4s/(3 pi): tip {w_pt:.4f} m")

    # elliptic lift: integration and FD
    y = np.linspace(0.0, s, 4001)
    M = bending_moment(y, elliptic_lift(y, L, s))
    th, w = deflection_by_integration(y, M, EI)
    print(f"elliptic lift: root moment {M[0]:.1f} N m, tip slope {np.degrees(th[-1]):.3f} deg, "
          f"tip deflection {w[-1]:.4f} m by integration ({w[-1] / s * 100:.2f} % of semi-span)")
    th_ex, w_ex = elliptic_tip(L, s, EI)
    print(f"closed form (superposed point loads): tip slope {np.degrees(th_ex):.3f} deg, tip {w_ex:.5f} m; "
          f"integration error {w[-1] - w_ex:.1e} m")
    print("finite differences, uniform EI: |FD tip - closed form| (m)")
    print("     n    uniform lift    elliptic lift")
    ns = np.array([25, 50, 100, 200, 400])
    errs_u, errs = [], []
    for n in ns:
        _, wu, _ = fd_beam(s, EI, lambda yy: np.full_like(yy, L / s), int(n))
        _, wf, mf = fd_beam(s, EI, lambda yy: elliptic_lift(yy, L, s), int(n))
        errs_u.append(abs(wu[-1] - w_udl))
        errs.append(abs(wf[-1] - w_ex))
        print(f"  {n:4d}    {errs_u[-1]:.3e}       {errs[-1]:.3e}")
    print(f"  observed orders, uniform:  {np.array2string(observed_order(errs_u), precision=2)}")
    print(f"  observed orders, elliptic: {np.array2string(observed_order(errs), precision=2)} "
          f"(load ~ sqrt(s - y) at the tip: order tends to 1.5)")

    # superposition: lift, wing weight, a fuel tank
    m_wing, m_tank, y_tank, n_lim = 45.0, 40.0, 2.0, 3.8
    q_wing = -n_lim * m_wing * G0 / s
    P_tank = -n_lim * m_tank * G0
    w_lift = w[-1]
    w_wt = float(cantilever_udl(s, s, q_wing, EI))
    w_tk = float(cantilever_point(s, s, y_tank, P_tank, EI))
    x, wf, mf = fd_beam(s, EI, lambda yy: elliptic_lift(yy, L, s) + q_wing, 200, point_loads=[(y_tank, P_tank)])
    print(f"\nsuperposition at n = {n_lim}: lift {w_lift:.4f} + wing weight {w_wt:.4f} ({m_wing:.0f} kg uniform) + "
          f"tank {w_tk:.4f} ({m_tank:.0f} kg at {y_tank} m) = {w_lift + w_wt + w_tk:.4f} m; FD with all loads {wf[-1]:.4f} m")
    print(f"root moment relieved from {M[0]:.0f} to {M[0] + q_wing * s ** 2 / 2 + P_tank * y_tank:.0f} N m "
          f"by superposition (FD {mf[0]:.0f} N m)")

    # tapered spar
    EIt = tapered_EI(EI, 0.7, s)
    th_t, w_t = deflection_by_integration(y, M, EIt(y))
    x, wf_t, _ = fd_beam(s, EIt, lambda yy: elliptic_lift(yy, L, s), 400)
    print(f"\ntapered spar (EI tip = 0.3 EI root): tip {w_t[-1]:.4f} m by integration, {wf_t[-1]:.4f} m by FD")

    # indeterminate: propped cantilever
    Lp, qp, EIp = 1.2, 2000.0, 2.0e4
    xp, wp, mp = fd_beam(Lp, EIp, lambda xx: np.full_like(xx, qp), 400, bc="clamped-pinned")
    print(f"\npropped cantilever L = {Lp} m, q = {qp:.0f} N/m: prop force 3qL/8 = {propped_cantilever_reaction(Lp, qp):.1f} N, "
          f"clamp moment qL^2/8 = {qp * Lp ** 2 / 8:.1f} N m (FD {mp[0]:.1f}), "
          f"max deflection {np.max(np.abs(propped_cantilever_udl(xp, Lp, qp, EIp))) * 1e3:.3f} mm "
          f"(FD {np.max(np.abs(wp)) * 1e3:.3f} mm)")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].plot(y, w * 1e3, label="elliptic lift, integration")
    x50, w50, _ = fd_beam(s, EI, lambda yy: elliptic_lift(yy, L, s), 25)
    ax[0].plot(x50, w50 * 1e3, "o", ms=4, label="elliptic lift, FD n = 25")
    ax[0].plot(y, cantilever_udl(y, s, L / s, EI) * 1e3, "--", label="uniform lift, closed form")
    ax[0].plot(y, w_t * 1e3, label="elliptic lift, tapered EI")
    ax[0].set_xlabel("span station y (m)")
    ax[0].set_ylabel("deflection w (mm)")
    ax[0].set_title("wing spar at limit load")
    ax[0].legend()
    ax[1].loglog(s / ns, errs_u, "o-", label="uniform lift")
    ax[1].loglog(s / ns, errs, "s-", label="elliptic lift")
    ax[1].loglog(s / ns, errs_u[0] * (ns[0] / ns) ** 2, "k--", lw=0.8, label="slope 2")
    ax[1].loglog(s / ns, errs[0] * (ns[0] / ns) ** 1.5, "k:", lw=0.8, label="slope 1.5")
    ax[1].set_xlabel("grid spacing h (m)")
    ax[1].set_ylabel("|FD tip - closed form| (m)")
    ax[1].set_title("finite-difference convergence")
    ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "wing_deflection.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'wing_deflection.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
