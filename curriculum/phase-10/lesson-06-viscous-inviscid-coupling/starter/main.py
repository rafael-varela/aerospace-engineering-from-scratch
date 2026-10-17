"""P10.06 Viscous-Inviscid Coupling: starter code.

A panel method (`P10.05`) gives the inviscid surface velocity u_e(s). A boundary layer grows under
that velocity, and its momentum deficit at the trailing edge is the section's profile drag. This
module couples them one way (inviscid -> viscous) and optionally feeds the displacement thickness
back once:

* the panel solution is split at the stagnation point into upper and lower surfaces, with arc
  length s measured from the stagnation point;
* laminar boundary layer by Thwaites' method (`P08.11`): theta^2 u_e^6 = 0.45 nu int u_e^5 ds,
  lambda = theta^2 u_e'/nu, with H(lambda) and l(lambda) from the Cebeci-Bradshaw fits;
  laminar separation at lambda < -0.09;
* transition by Michel's criterion (`P08.12`):
  Re_theta >= 1.174 (1 + 22400/Re_x) Re_x^0.46, or at laminar separation if that comes first;
* turbulent boundary layer by Head's entrainment method (`P08.14`), two ODEs for theta and the
  entrainment shape factor H1, closed by the Ludwieg-Tillmann skin friction, integrated with
  scipy.integrate.solve_ivp (`P03.05`);
* profile drag by the Squire-Young formula cd = 2 theta_TE/c (u_e,TE/V)^((H_TE + 5)/2), summed
  over both surfaces;
* checks: Falkner-Skan similarity solutions (Blasius, Hiemenz) by solve_bvp against Thwaites;
  Head against the 1/7-power-law flat plate; the NACA 0012 prediction against the authentic
  Ladson (NASA TM-4074) free-transition drag data in curriculum/phase-10/data/.

Units: chord c = 1 and V = 1 inside the boundary-layer code, so nu = 1/Re; SI elsewhere.

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes viscous_drag.png)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import integrate, interpolate

NU_ISA = 1.461e-5
TRAINER = dict(mass=750.0, S=12.0, c=1.2, V=55.0, rho=1.225, g=9.81)
LADSON_FILE = "ladson_naca0012_M015.csv"


# ---------------------------------------------------------------- given: geometry and panel method (P10.02, P10.05)

def naca4(code: str, n: int = 100) -> tuple[np.ndarray, np.ndarray]:
    """Closed-trailing-edge NACA 4-digit coordinates, Selig order, n panels per surface (P10.02, given)."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    if m > 0:
        yc = np.where(x <= p, m / p ** 2 * (2 * p * x - x ** 2), m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x ** 2))
        dyc = np.where(x <= p, 2 * m / p ** 2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))
    else:
        yc, dyc = np.zeros_like(x), np.zeros_like(x)
    th = np.arctan(dyc)
    xu, yu, xl, yl = x - yt * np.sin(th), yc + yt * np.cos(th), x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xu[::-1], xl[1:]]), np.concatenate([yu[::-1], yl[1:]])


def panel_solve(x: np.ndarray, y: np.ndarray, alpha: float) -> dict:
    """Compact linear-strength vortex panel method of P10.05 (given). Returns node gamma (clockwise
    positive), the signed surface velocity at the nodes along the node order, Vt = -gamma, cl from
    the circulation, and the panel normals. V = 1, chord ~ 1."""
    dx, dy = np.diff(x), np.diff(y)
    L = np.hypot(dx, dy)
    tx, ty = dx / L, dy / L
    nx, ny = ty, -tx
    xc, yc = 0.5 * (x[:-1] + x[1:]), 0.5 * (y[:-1] + y[1:])
    N = L.size
    rx, ry = xc[:, None] - x[None, :-1], yc[:, None] - y[None, :-1]
    xi, eta = rx * tx + ry * ty, -rx * ty + ry * tx
    k = np.arange(N)
    eta[k, k] = -0.0
    beta = np.arctan2(eta, xi - L) - np.arctan2(eta, xi)
    beta[k, k] = -np.pi
    lnr = 0.5 * np.log((xi ** 2 + eta ** 2) / ((xi - L) ** 2 + eta ** 2))
    I3, I4 = -eta * lnr + xi * beta, -L + eta * beta + xi * lnr
    c = 1 / (2 * np.pi)
    ua, ub, va, vb = c * (beta - I3 / L), c * I3 / L, -c * (lnr - I4 / L), -c * I4 / L
    rot = lambda ul, vl: (ul * tx - vl * ty, ul * ty + vl * tx)
    (uA, vA), (uB, vB) = rot(ua, va), rot(ub, vb)
    U, V = np.zeros((N, N + 1)), np.zeros((N, N + 1))
    U[:, :N] += uA; U[:, 1:] += uB; V[:, :N] += vA; V[:, 1:] += vB
    A = np.zeros((N + 1, N + 1))
    A[:N] = U * nx[:, None] + V * ny[:, None]
    A[N, 0] = A[N, N] = 1.0
    rhs = np.concatenate([-(np.cos(alpha) * nx + np.sin(alpha) * ny), [0.0]])
    g = np.linalg.solve(A, rhs)
    chord = float(np.max(np.hypot(x - 0.5 * (x[0] + x[-1]), y - 0.5 * (y[0] + y[-1]))))
    return dict(gamma=g, Vt=-g, cl=float(2 * np.sum(L * 0.5 * (g[:-1] + g[1:])) / chord), nx=nx, ny=ny, L=L,
                chord=chord)


# ---------------------------------------------------------------- given: data access

def find_data_file(name: str = LADSON_FILE, start: str | Path | None = None) -> Path | None:
    """Search upward from this file (or start) for curriculum/phase-10/data/<name>; works from both
    the curriculum/ and the solutions/ copies of the lesson. Returns None if not found (given)."""
    here = Path(start or __file__).resolve()
    for d in [here] + list(here.parents):
        cand = d / "curriculum" / "phase-10" / "data" / name
        if cand.exists():
            return cand
    return None


def read_ladson(path: str | Path) -> dict:
    """Read the Ladson NACA 0012 CSV (columns table, transition, mach, reynolds, alpha_deg, cd, cl,
    cm, l_over_d; '#' lines are provenance; blank cd = post-stall, returned as NaN) (given)."""
    lines = [ln for ln in Path(path).read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    cols = lines[0].split(",")
    rows = [ln.split(",") for ln in lines[1:]]
    out = {}
    for i, cname in enumerate(cols):
        vals = [r[i] for r in rows]
        if cname in ("table", "transition"):
            out[cname] = np.array(vals)
        else:
            out[cname] = np.array([float(v) if v != "" else np.nan for v in vals])
    return out


def ladson_zero_lift_drag(data: dict, table: str = "I", max_abs_alpha: float = 0.15) -> tuple[np.ndarray, np.ndarray]:
    """Reynolds numbers and cd near alpha = 0 (|alpha| <= max_abs_alpha deg) from one table, averaging
    repeated points at the same Re (given; values are never altered)."""
    sel = (data["table"] == table) & (np.abs(data["alpha_deg"]) <= max_abs_alpha) & np.isfinite(data["cd"])
    Re = np.unique(data["reynolds"][sel])
    cd = np.array([data["cd"][sel & (data["reynolds"] == r)].mean() for r in Re])
    return Re, cd


# ---------------------------------------------------------------- similarity solutions (independent check)

def falkner_skan(beta: float, eta_max: float = 10.0) -> dict:
    """Falkner-Skan similarity solution f''' + f f'' + beta (1 - f'^2) = 0, f(0) = f'(0) = 0,
    f'(inf) = 1, by scipy.integrate.solve_bvp. With eta = y sqrt((m+1) u_e / (2 nu x)) and
    beta = 2m/(m+1): beta = 0 is Blasius, beta = 1 the Hiemenz stagnation flow.
    Returns f''(0) and the integrals I_theta = int f'(1 - f') d eta, I_delta = int (1 - f') d eta,
    H = I_delta / I_theta, so theta = I_theta sqrt(2 nu x / ((m+1) u_e))."""
    raise NotImplementedError


# ---------------------------------------------------------------- laminar: Thwaites

def thwaites_H_l(lam: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Shape factor H and shear function l = (theta/u_e) du/dy|_wall as functions of Thwaites'
    lambda (Cebeci-Bradshaw fits of Thwaites' tables):
        0 <= lambda <= 0.1:  l = 0.22 + 1.57 lambda - 1.8 lambda^2,  H = 2.61 - 3.75 lambda + 5.24 lambda^2;
        -0.1 <= lambda < 0:  l = 0.22 + 1.402 lambda + 0.018 lambda/(lambda + 0.107),  H = 2.088 + 0.0731/(lambda + 0.14).
    lambda is clipped to [-0.1, 0.1]."""
    raise NotImplementedError


def thwaites(s: np.ndarray, ue: np.ndarray, nu: float) -> dict:
    """Thwaites' laminar boundary layer along arc length s (s[0] = 0) with edge velocity ue (>= 0).

    theta^2(s) = 0.45 nu / ue^6 int_0^s ue^5 ds, the integral taken exactly for ue linear between
    stations, int ue^5 ds = ds (u0^5 + u0^4 u1 + ... + u1^5)/6 (the trapezoid rule is three times
    too large in the first interval after a stagnation point, where ue ~ s). If ue[0] = 0 (a stagnation
    point) the limit theta^2 = 0.075 nu / ue'(0) is used there. due/ds by np.gradient (second order on
    a non-uniform grid). Returns theta, lam, H, cf = 2 nu l /(ue theta), Re_theta = ue theta/nu,
    due_ds and i_sep: the first index with lambda < -0.09 (laminar separation), or None."""
    raise NotImplementedError


def michel_transition(Re_x: np.ndarray, Re_theta: np.ndarray) -> np.ndarray:
    """Michel's criterion: True where Re_theta >= 1.174 (1 + 22400/Re_x) Re_x^0.46."""
    raise NotImplementedError


# ---------------------------------------------------------------- turbulent: Head

def head_H1(H: np.ndarray) -> np.ndarray:
    """Head's entrainment shape factor H1 = (delta - delta*)/theta as a function of H (Cebeci-Bradshaw fit):
    H1 = 3.3 + 0.8234 (H - 1.1)^-1.287 for H <= 1.6, 3.3 + 1.5501 (H - 0.6778)^-3.064 for H > 1.6."""
    raise NotImplementedError


def head_H_of_H1(H1: np.ndarray) -> np.ndarray:
    """Inverse of ``head_H1``: H = 1.1 + ((H1 - 3.3)/0.8234)^(-1/1.287) for H1 >= 5.3,
    0.6778 + ((H1 - 3.3)/1.5501)^(-1/3.064) for H1 < 5.3 (H1 is floored just above 3.3)."""
    raise NotImplementedError


def ludwieg_tillmann_cf(H: np.ndarray, Re_theta: np.ndarray) -> np.ndarray:
    """Turbulent skin friction cf = 0.246 10^(-0.678 H) Re_theta^-0.268."""
    raise NotImplementedError


def head_turbulent(s: np.ndarray, ue: np.ndarray, nu: float, s0: float, theta0: float, H0: float = 1.4,
                   H_sep: float = 2.4) -> dict:
    """Head's method from s0 (transition) to the end of s, with ue interpolated by a cubic spline:

        d theta/ds = cf/2 - (H + 2) (theta/ue) due/ds,
        d(ue theta H1)/ds = ue * 0.0306 (H1 - 3)^-0.6169,

    integrated for y = [theta, ue theta H1] with solve_ivp (RK45, rtol 1e-8). The event H >= H_sep
    stops the integration (turbulent separation). Returns s, theta, H, cf at the output stations
    (the grid points beyond s0) and s_sep (None if attached)."""
    raise NotImplementedError


# ---------------------------------------------------------------- coupling

def split_at_stagnation(x: np.ndarray, y: np.ndarray, Vt: np.ndarray) -> tuple[dict, dict]:
    """Split a Selig-ordered panel solution at the stagnation point. Vt is the signed surface
    velocity at the nodes along the node order: negative on the upper surface (flow runs against
    the order, towards the trailing edge at node 0), positive on the lower. The stagnation point
    lies where Vt changes sign from - to +, linearly interpolated in arc length. Returns
    (upper, lower), each with s (from the stagnation point, starting with 0), ue = |Vt| (starting
    with 0), x (chordwise position) and the node indices used."""
    raise NotImplementedError


def surface_boundary_layer(s: np.ndarray, ue: np.ndarray, nu: float, x: np.ndarray | None = None,
                           force_transition_s: float | None = None) -> dict:
    """Boundary layer along one surface: Thwaites from the stagnation point; transition at the first
    station where Michel's criterion holds, or at laminar separation if that comes first (a short
    bubble assumed to reattach turbulent), or at force_transition_s (a trip); then Head with theta
    continuous and H = 1.4. Returns the trailing-edge theta, H, ue, the transition arc length and
    chordwise position, the mechanism ('michel', 'laminar separation', 'forced', 'none'), and
    s_sep (turbulent separation, None if attached)."""
    raise NotImplementedError


def squire_young(theta_te: float, H_te: float, ue_te: float, V: float = 1.0, chord: float = 1.0) -> float:
    """Squire-Young drag of one surface: cd = 2 (theta_TE/c) (ue_TE/V)^((H_TE + 5)/2)."""
    raise NotImplementedError


def trim_branch(br: dict, x_end: float) -> dict:
    """Keep the stations of a surface branch up to the last one with x <= x_end (given). The inviscid
    flow decelerates sharply into the closed trailing-edge corner over the last ~1 % of chord, an
    artefact the real (viscous) flow does not see; the boundary layer is stopped just ahead of it."""
    k = int(np.where(br["x"] <= x_end)[0][-1]) + 1
    return dict(s=br["s"][:k], ue=br["ue"][:k], x=br["x"][:k], idx=br["idx"][:k - 1])


def section_drag(code: str, alpha: float, Re: float, n: int = 100, trip_x: float | None = None,
                 x_end: float = 0.99) -> dict:
    """Profile drag of a NACA 4-digit section: panel solution (n panels per surface), split at the
    stagnation point, each branch trimmed at x_end, boundary layer on each surface, Squire-Young on
    each (at x_end), summed. trip_x forces transition at that chordwise position on both surfaces
    (fixed transition). Returns cd, the upper and lower results, cl (inviscid), the panel solution
    and the branches."""
    raise NotImplementedError


def displacement_iteration(code: str, alpha: float, Re: float, n: int = 100) -> dict:
    """One viscous-inviscid iteration: add the displacement thickness delta* = H theta (laminar and
    turbulent parts) to the surface along the outward normal at each node, re-run the panel method
    on the thickened body, and recompute the drag. Returns cl and cd before and after."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed orders log(e_k/e_{k+1})/log(ratio) (given)."""
    e = np.asarray(errors, float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".", n: int = 100) -> None:
    """Print the lesson's numbers and save viscous_drag.png (n panels per surface; the tests pass a
    smaller n to run quickly)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    # similarity checks
    bl = falkner_skan(0.0)
    hz = falkner_skan(1.0)
    print(f"Blasius (Falkner-Skan beta = 0): f''(0) = {bl['fpp0']:.5f}, theta sqrt(Re_x)/x = {bl['I_theta'] * np.sqrt(2):.5f}, "
          f"H = {bl['H']:.4f}; cf sqrt(Re_x) = {2 * bl['fpp0'] / np.sqrt(2):.5f}")
    s = np.linspace(0, 1, 2001)
    nu = 1e-6
    fp = thwaites(s, np.ones_like(s), nu)
    print(f"Thwaites flat plate: theta sqrt(Re_x)/x = {fp['theta'][-1] * np.sqrt(1 / nu) / 1:.5f} (sqrt 0.45 = {np.sqrt(0.45):.5f}), "
          f"H = {fp['H'][-1]:.3f}, cf sqrt(Re_x) = {fp['cf'][-1] * np.sqrt(1 / nu):.4f}")
    print(f"Hiemenz (beta = 1): theta sqrt(C/nu) = {hz['I_theta']:.5f}, H = {hz['H']:.4f}, f''(0) = {hz['fpp0']:.5f}")
    st = thwaites(s, 3.0 * s, nu)
    print(f"Thwaites stagnation flow: theta sqrt(C/nu) = {st['theta'][1000] * np.sqrt(3.0 / nu):.5f} (sqrt 0.075 = "
          f"{np.sqrt(0.075):.5f}), lambda = {st['lam'][1000]:.4f}, H = {st['H'][1000]:.3f}")
    # Head on a flat plate from a laminar start, against the 1/7-power law
    Re_L = 1e7
    xs = np.linspace(0, 1, 4001)
    lamfp = thwaites(xs, np.ones_like(xs), 1 / Re_L)
    i0 = 40
    hd = head_turbulent(xs, np.ones_like(xs), 1 / Re_L, xs[i0], lamfp["theta"][i0])
    th_pl = 0.036 * hd["s"] * (hd["s"] * Re_L) ** -0.2
    print(f"Head flat plate (Re_L = 1e7, start x = {xs[i0]:.3f}): theta(L) = {hd['theta'][-1]:.6f}, 1/7-power law "
          f"0.036 L Re_L^-0.2 = {th_pl[-1]:.6f}; H(L) = {hd['H'][-1]:.3f}; cf(L) = {hd['cf'][-1]:.5f} vs "
          f"0.0592 Re_x^-0.2 = {0.0592 * Re_L ** -0.2:.5f}")

    # NACA 0012 drag vs Re, against Ladson
    print(f"\nNACA 0012, alpha = 0, free transition (Michel), N = {2 * n} panels:")
    Res = [2e6, 3.94e6, 5.97e6, 8.86e6, 11.9e6, 18.9e6]
    cds = []
    for Re in Res:
        r = section_drag("0012", 0.0, Re, n)
        cds.append(r["cd"])
        print(f"   Re = {Re:.3g}: cd = {r['cd']:.5f}; transition x/c = {r['upper']['x_tr']:.3f} ({r['upper']['mechanism']})")
    trip = [section_drag("0012", 0.0, Re, n, trip_x=0.05)["cd"] for Re in Res]
    print("   tripped at 5 % chord: " + ", ".join(f"{c:.5f}" for c in trip))
    path = find_data_file()
    if path is not None:
        data = read_ladson(path)
        ReL, cdL = ladson_zero_lift_drag(data, "I")
        ReF, cdF = ladson_zero_lift_drag(data, "VII")
        print(f"Ladson (NASA TM-4074) Table I, free transition, |alpha| <= 0.15 deg ({path.name}):")
        for R, c in zip(ReL, cdL):
            pred = section_drag("0012", 0.0, R, n)["cd"]
            print(f"   Re = {R:.3g}: measured cd = {c:.5f}; predicted {pred:.5f} ({(pred / c - 1) * 100:+.1f} %)")
        print("Table VII (wraparound grit at 5 % chord): " + ", ".join(f"Re {R:.3g}: {c:.5f}" for R, c in zip(ReF, cdF)))
    else:
        ReL = cdL = ReF = cdF = None
        print("Ladson data file not found; comparison skipped")

    # convergence with panel count
    print("\npanel-count convergence, NACA 0012, Re = 5.97e6:")
    cdn = []
    for nn in (n // 2, n, 2 * n):
        cdn.append(section_drag("0012", 0.0, 5.97e6, nn)["cd"])
        print(f"   {2 * nn} panels: cd = {cdn[-1]:.6f}")

    # the trainer
    t = TRAINER
    Re_tr = t["V"] * t["c"] / NU_ISA
    cl_need = t["mass"] * t["g"] / (0.5 * t["rho"] * t["V"] ** 2 * t["S"])
    a_lo, a_hi = 0.0, np.radians(4.0)
    c_lo, c_hi = section_drag("2412", a_lo, Re_tr, n)["cl"], section_drag("2412", a_hi, Re_tr, n)["cl"]
    a_cr = a_lo + (cl_need - c_lo) / (c_hi - c_lo) * (a_hi - a_lo)
    r = section_drag("2412", a_cr, Re_tr, n)
    print(f"\ntrainer NACA 2412 at cruise: Re = {Re_tr:.3g}, inviscid cl = {r['cl']:.4f} at alpha = {np.degrees(a_cr):.3f} deg: "
          f"cd = {r['cd']:.5f} (upper {r['upper']['cd']:.5f}, transition x/c {r['upper']['x_tr']:.3f} {r['upper']['mechanism']}; "
          f"lower {r['lower']['cd']:.5f}, x/c {r['lower']['x_tr']:.3f} {r['lower']['mechanism']})")
    rt = section_drag("2412", a_cr, Re_tr, n, trip_x=0.05)
    print(f"   tripped at 5 % chord (insects, rain, rivets): cd = {rt['cd']:.5f}")
    D = rt["cd"] * 0.5 * t["rho"] * t["V"] ** 2 * t["S"]
    D0 = r["cd"] * 0.5 * t["rho"] * t["V"] ** 2 * t["S"]
    print(f"   wing profile drag at cruise: {D0:.0f} N clean, {D:.0f} N tripped (section cd x q x S)")
    di = displacement_iteration("2412", a_cr, Re_tr, n)
    print(f"   one displacement-thickness iteration: cl {di['cl0']:.4f} -> {di['cl1']:.4f} "
          f"({(di['cl1'] / di['cl0'] - 1) * 100:+.1f} %)")
    print("alpha sweep, NACA 2412, Re = trainer cruise:")
    for a_deg in (0.0, 2.0, 4.0, 6.0, 8.0):
        ra = section_drag("2412", np.radians(a_deg), Re_tr, n)
        print(f"   alpha = {a_deg:3.0f} deg: cl = {ra['cl']:.3f}, cd = {ra['cd']:.5f}, x_tr upper {ra['upper']['x_tr']:.3f} "
              f"({ra['upper']['mechanism']}), lower {ra['lower']['x_tr']:.3f}; turbulent separation upper: "
              f"{'none' if ra['upper']['s_sep'] is None else 'yes'}")

    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    ax[0].semilogx(Res, cds, "o-", label="this model, free transition (Michel)")
    ax[0].semilogx(Res, trip, "s--", label="this model, tripped at 0.05c")
    if ReL is not None:
        ax[0].semilogx(ReL, cdL, "k^", label="Ladson TM-4074, free transition (measured)")
        ax[0].semilogx(ReF, cdF, "kv", mfc="none", label="Ladson, wraparound grit (measured)")
    ax[0].set_xlabel("Re"); ax[0].set_ylabel("cd"); ax[0].legend(fontsize=8); ax[0].set_title("NACA 0012, alpha = 0")
    r = section_drag("0012", 0.0, 5.97e6, n)
    up = r["branches"][0]
    lamu = r["upper"]["laminar"]
    tu = r["upper"]["turbulent"]
    itr = r["upper"]["i_tr"]
    ax[1].plot(up["x"][:itr + 1], lamu["theta"][:itr + 1] * 1e3, label="theta laminar (Thwaites)")
    xt = np.interp(tu["s"], up["s"], up["x"])
    ax[1].plot(xt, tu["theta"] * 1e3, label="theta turbulent (Head)")
    ax[1].plot(xt, tu["H"], label="H turbulent")
    ax[1].plot(up["x"][:itr + 1], lamu["H"][:itr + 1], label="H laminar")
    ax[1].set_xlabel("x/c"); ax[1].set_ylabel("theta/c x 1000, H"); ax[1].legend(fontsize=8)
    ax[1].set_title("NACA 0012 upper surface, Re = 5.97e6")
    fig.tight_layout()
    fig.savefig(out_dir / "viscous_drag.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'viscous_drag.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
