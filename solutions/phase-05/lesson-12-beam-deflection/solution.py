"""P05.12 Beam Deflection: reference solution.

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

Run:  python solution.py [out_dir]     (prints the lesson's numbers, writes wing_deflection.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-05/lesson-12-beam-deflection
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
    x = np.asarray(x, dtype=float)
    return np.where(x <= a, P * x ** 2 * (3 * a - x), P * a ** 2 * (3 * x - a)) / (6 * EI)


def cantilever_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Cantilever clamped at x = 0 under uniform load q (up): w = q x^2 (6L^2 - 4Lx + x^2) / (24 EI);
    tip deflection q L^4 / (8 EI)."""
    x = np.asarray(x, dtype=float)
    return q * x ** 2 * (6 * L ** 2 - 4 * L * x + x ** 2) / (24 * EI)


def cantilever_end_moment(x: np.ndarray, L: float, M0: float, EI: float) -> np.ndarray:
    """Cantilever with a sagging-positive moment M0 applied at the free end: constant curvature,
    w = M0 x^2 / (2 EI)."""
    x = np.asarray(x, dtype=float)
    return M0 * x ** 2 / (2 * EI)


def simply_supported_point(x: np.ndarray, L: float, a: float, P: float, EI: float) -> np.ndarray:
    """Simply supported beam (pins at 0 and L), point load P (up) at x = a, b = L - a:
    w = P b x (L^2 - b^2 - x^2) / (6 L EI) for x <= a, and the mirror image beyond."""
    x = np.asarray(x, dtype=float)
    b = L - a
    left = P * b * x * (L ** 2 - b ** 2 - x ** 2) / (6 * L * EI)
    xr = L - x
    right = P * a * xr * (L ** 2 - a ** 2 - xr ** 2) / (6 * L * EI)
    return np.where(x <= a, left, right)


def simply_supported_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Simply supported beam under uniform load q (up): w = q x (L^3 - 2 L x^2 + x^3) / (24 EI);
    midspan deflection 5 q L^4 / (384 EI)."""
    x = np.asarray(x, dtype=float)
    return q * x * (L ** 3 - 2 * L * x ** 2 + x ** 3) / (24 * EI)


def propped_cantilever_reaction(L: float, q: float) -> float:
    """Prop force at x = L of a cantilever (clamped at 0) under uniform q, by superposition: the prop
    must cancel the free cantilever's tip deflection, R L^3 / (3 EI) = q L^4 / (8 EI), so R = 3 q L / 8
    (acting against the load; returned as a positive magnitude)."""
    return 3.0 * q * L / 8.0


def propped_cantilever_udl(x: np.ndarray, L: float, q: float, EI: float) -> np.ndarray:
    """Deflection of the propped cantilever (clamped at 0, pinned at L) under uniform q, by
    superposition of ``cantilever_udl`` and a tip point load -R: w = q x^2 (3L^2 - 5Lx + 2x^2) / (48 EI)."""
    R = propped_cantilever_reaction(L, q)
    return cantilever_udl(x, L, q, EI) + cantilever_point(x, L, L, -R, EI)


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
    x = np.asarray(x, dtype=float)
    q = np.broadcast_to(np.asarray(q, dtype=float), x.shape)
    I1 = _cumtrapz(q, x)
    M = _cumtrapz(I1, x)
    slope_end = I1[-1]
    for a, P in point_loads:
        M = M + P * np.clip(x - a, 0.0, None)
        slope_end += P
    L = x[-1] - x[0]
    if support == "cantilever":
        C1 = -slope_end
        C0 = -(M[-1] + C1 * L)
    elif support == "simply_supported":
        C0 = 0.0
        C1 = -M[-1] / L
    else:
        raise ValueError("support must be 'cantilever' or 'simply_supported'")
    return M + C1 * (x - x[0]) + C0


def deflection_by_integration(x: np.ndarray, M: np.ndarray, EI: np.ndarray | float,
                              support: str = "cantilever") -> tuple[np.ndarray, np.ndarray]:
    """Slope and deflection from the curvature kappa = M / EI by two cumulative trapezoid integrals,
    with the support's kinematic conditions: 'cantilever' (w = w' = 0 at x[0]) or 'simply_supported'
    (w = 0 at both ends, the slope constant fixed by the second condition). EI may vary along x.
    Returns (theta, w)."""
    x = np.asarray(x, dtype=float)
    kappa = np.asarray(M, dtype=float) / np.broadcast_to(np.asarray(EI, dtype=float), x.shape)
    theta = _cumtrapz(kappa, x)
    w = _cumtrapz(theta, x)
    if support == "simply_supported":
        c = -w[-1] / (x[-1] - x[0])
        theta = theta + c
        w = w + c * (x - x[0])
    elif support != "cantilever":
        raise ValueError("support must be 'cantilever' or 'simply_supported'")
    return theta, w


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
    x = np.linspace(0.0, L, n + 1)
    h = L / n
    N = n + 1
    EIx = EI(x) if callable(EI) else np.full(N, float(EI))
    qx = np.asarray(q(x), dtype=float) * np.ones(N)
    q_dist = qx.copy()                # distributed load only, for the end shear rows
    P_end = {"L": 0.0, "R": 0.0}
    for a, P in point_loads:
        j = int(round(a / h))
        if j == 0:
            P_end["L"] += P
        elif j == N - 1:
            P_end["R"] += P
        else:
            qx[j] += P / h
    A = np.zeros((2 * N, 2 * N))
    b = np.zeros(2 * N)
    W = lambda i: i                   # column of w_i
    Mc = lambda i: N + i              # column of m_i
    row = 0
    for i in range(1, N - 1):         # interior: m'' = q and w'' = m / EI
        A[row, Mc(i - 1)], A[row, Mc(i)], A[row, Mc(i + 1)] = 1 / h ** 2, -2 / h ** 2, 1 / h ** 2
        b[row] = qx[i]
        row += 1
        A[row, W(i - 1)], A[row, W(i)], A[row, W(i + 1)] = 1 / h ** 2, -2 / h ** 2, 1 / h ** 2
        A[row, Mc(i)] = -1.0 / EIx[i]
        row += 1

    def set_row(cols, vals, rhs=0.0):
        nonlocal row
        for c, v in zip(cols, vals):
            A[row, c] += v
        b[row] = rhs
        row += 1

    left, right = bc.split("-")
    for end, kind in (("L", left), ("R", right)):
        i0, i1 = (0, 1) if end == "L" else (N - 1, N - 2)
        sgn = 1.0 if end == "L" else -1.0
        # Taylor: f'(end) = (f_L1 - f_L0)/h - (h/2) f''(0) at the left end, and
        #         f'(end) = (f_R0 - f_R1)/h + (h/2) f''(L) at the right end, both to O(h^2).
        if kind == "clamped":                  # w = 0, w' = 0 (w'' = m / EI)
            set_row([W(i0)], [1.0])
            set_row([W(i1), W(i0), Mc(i0)], [sgn / h, -sgn / h, -sgn * 0.5 * h / EIx[i0]])
        elif kind == "pinned":                 # w = 0, m = 0
            set_row([W(i0)], [1.0])
            set_row([Mc(i0)], [1.0])
        elif kind == "free":                   # m = 0, m' = V_end
            set_row([Mc(i0)], [1.0])
            V_end = P_end["L"] if end == "L" else -P_end["R"]
            set_row([Mc(i1), Mc(i0)], [sgn / h, -sgn / h], V_end + sgn * 0.5 * h * q_dist[i0])
        else:
            raise ValueError(f"unknown end condition {kind!r}")
    if row != 2 * N:
        raise ValueError("boundary conditions do not close the system")
    sol = np.linalg.solve(A, b)
    return x, sol[:N], sol[N:]


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
    return L_semi * s ** 2 / (8 * EI), L_semi * s ** 3 / EI * (1.0 / 8.0 - 4.0 / (45.0 * np.pi))


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
