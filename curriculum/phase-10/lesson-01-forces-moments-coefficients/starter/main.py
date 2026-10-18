"""P10.01 Forces, Moments & Coefficients: starter code.

A wind-tunnel balance does not measure lift, drag and pitching moment. It measures forces along its
own axes and a moment about its own centre. This module turns balance readings into the
coefficients an aerodynamicist uses, and says how uncertain they are:

* dynamic pressure, Reynolds number and Mach number (the similarity parameters of `P08.05`);
* body (balance) axes -> wind axes: normal and axial force -> lift and drag;
* moment transfer from the balance moment centre to any reference point (the quarter chord);
* cl, cd, cm from forces, the dynamic pressure and the reference area and chord;
* the aerodynamic centre from the slope of cm against the normal-force coefficient;
* uncertainty propagation (`P00.14`): first-order (Jacobian) and Monte Carlo.

The dataset is SYNTHETIC: `synthetic_balance_data(seed)` invents a plausible NACA 2412 section
model in a 2-D tunnel, adds noise, and records what an imaginary balance would read. It is not a
measurement; every table and plot that uses it says so.

Axes and signs (Phase 10 conventions): x along the chord from the leading edge (aft positive),
z up; alpha positive nose-up; N (normal force) along +z, A (axial force) along +x (aft); lift
perpendicular and drag parallel to the free stream; pitching moment positive nose-up.

Units: SI (N, m, Pa, rad inside the functions; degrees only where the name says _deg).

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes balance_coefficients.png
                                      and SYNTHETIC_seed042_balance.csv into out_dir)
"""
from __future__ import annotations

import sys
from pathlib import Path
from typing import Callable

import numpy as np

# ISA sea level (Phase 10 shared conventions)
RHO_ISA = 1.225          # kg/m^3
MU_ISA = 1.789e-5        # Pa s
NU_ISA = 1.461e-5        # m^2/s
A_ISA = 340.3            # m/s
G0 = 9.81                # m/s^2
# room air, 20 C, 101.3 kPa (used for the synthetic tunnel)
RHO_ROOM = 1.204
MU_ROOM = 1.81e-5
A_ROOM = 343.2           # m/s, speed of sound at 20 C

# the reference trainer
TRAINER = dict(mass=750.0, S=12.0, b=10.0, c=1.2, V_cruise=55.0, V_approach=30.0)

# the SYNTHETIC tunnel model: a 2-D NACA 2412 section spanning the test section
MODEL = dict(c=0.30, span=0.60, V=35.0, x_bal=0.120, z_bal=-0.020)  # balance centre 0.40c aft, 20 mm below chord
# invented "true" section behaviour used to generate the synthetic readings (plausible, not measured)
TRUE_SECTION = dict(a0=6.0, alpha_L0_deg=-2.1, cd0=0.0062, k_cd=0.0060, cl_cd_min=0.25, x_ac=0.248, cm_ac=-0.050)
# invented standard uncertainties of the imaginary instruments (1 sigma)
NOISE = dict(N=0.20, A=0.05, M=0.010, q_rel=0.005, alpha_deg=0.05)


# ---------------------------------------------------------------- similarity parameters

def dynamic_pressure(rho: float | np.ndarray, V: float | np.ndarray) -> float | np.ndarray:
    """Dynamic pressure q = rho V^2 / 2 (Pa)."""
    raise NotImplementedError


def reynolds_number(V: float | np.ndarray, length: float, nu: float = NU_ISA) -> float | np.ndarray:
    """Reynolds number Re = V l / nu based on the reference length l (the chord for a section)."""
    raise NotImplementedError


def mach_number(V: float | np.ndarray, a: float = A_ISA) -> float | np.ndarray:
    """Mach number M = V / a."""
    raise NotImplementedError


# ---------------------------------------------------------------- axes and moments

def body_to_wind(N: np.ndarray, A: np.ndarray, alpha: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rotate body-axis normal force N (up, perpendicular to the chord) and axial force A (aft, along
    the chord) into lift L (perpendicular to the stream) and drag D (along the stream):

        L = N cos(alpha) - A sin(alpha),   D = N sin(alpha) + A cos(alpha),   alpha in radians."""
    raise NotImplementedError


def wind_to_body(L: np.ndarray, D: np.ndarray, alpha: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Inverse of ``body_to_wind``: N = L cos(alpha) + D sin(alpha), A = -L sin(alpha) + D cos(alpha)."""
    raise NotImplementedError


def transfer_moment(M: np.ndarray, N: np.ndarray, A: np.ndarray, x_from: float, z_from: float,
                    x_to: float, z_to: float) -> np.ndarray:
    """Pitching moment (nose-up positive) about (x_to, z_to) of a force system given as N, A acting
    through (x_from, z_from) plus a couple M. With x aft and z up, a force (F_x, F_z) at r has the
    nose-up moment (r_z - z_P) F_x - (r_x - x_P) F_z about P, so

        M_to = M + (x_to - x_from) N - (z_to - z_from) A."""
    raise NotImplementedError


def coefficients(L: np.ndarray, D: np.ndarray, M: np.ndarray, q: np.ndarray, S: float,
                 c: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Non-dimensional coefficients: C_L = L/(q S), C_D = D/(q S), C_m = M/(q S c)."""
    raise NotImplementedError


def reduce_balance(N: np.ndarray, A: np.ndarray, M_bal: np.ndarray, alpha_deg: np.ndarray, q: np.ndarray,
                   S: float, c: float, x_bal: float, z_bal: float, x_ref: float = None,
                   z_ref: float = 0.0) -> dict:
    """Balance readings -> coefficients. Moments are transferred from the balance centre
    (x_bal, z_bal) to the reference point (x_ref, z_ref) (default the quarter chord on the chord
    line), then N, A are rotated to L, D and everything is divided by q S (and c).

    Returns a dict with cl, cd, cm (about the reference point), cn, ca (body-axis coefficients),
    L, D, M_ref (N, N, N m)."""
    raise NotImplementedError


def fit_aerodynamic_centre(cn: np.ndarray, cm: np.ndarray, x_ref: float, c: float) -> tuple[float, float, float]:
    """Aerodynamic centre from a straight-line fit cm_ref = cm_ac + k cn. Transferring the moment from
    the a.c. to x_ref gives cm_ref = cm_ac + (x_ref - x_ac)/c * cn, so x_ac = x_ref - k c.

    Returns (x_ac, cm_ac, k). Use cn (exact); cl in place of cn is the usual small-angle shortcut."""
    raise NotImplementedError


# ---------------------------------------------------------------- uncertainty (P00.14)

def propagate_uncertainty(func: Callable[[np.ndarray], np.ndarray], x: np.ndarray, u: np.ndarray,
                          rel_step: float = 1e-6) -> tuple[np.ndarray, np.ndarray]:
    """First-order (GUM) propagation of independent standard uncertainties u through y = func(x).

    The sensitivity coefficients J_ij = d y_i / d x_j are central differences with step
    h_j = rel_step * max(|x_j|, u_j). Returns (u_y, contributions) where contributions[i, j] =
    (J_ij u_j)^2 is input j's share of the variance of output i, and u_y = sqrt(sum_j contributions)."""
    raise NotImplementedError


def monte_carlo_uncertainty(func: Callable[[np.ndarray], np.ndarray], x: np.ndarray, u: np.ndarray,
                            n: int = 20000, seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Monte Carlo propagation: draw n samples of x ~ Normal(x, u) (independent), evaluate func on
    each, and return (mean, standard deviation) of the outputs."""
    raise NotImplementedError


def point_coefficients(inputs: np.ndarray, S: float, c: float, x_bal: float, z_bal: float,
                       x_ref: float) -> np.ndarray:
    """(cl, cd, cm) of one balance point from inputs = [N, A, M_bal, alpha_deg, q]; the function the
    uncertainty routines differentiate."""
    N, A, M, a_deg, q = inputs
    r = reduce_balance(N, A, M, a_deg, q, S, c, x_bal, z_bal, x_ref)
    return np.array([float(r["cl"]), float(r["cd"]), float(r["cm"])])


# ---------------------------------------------------------------- the SYNTHETIC dataset (given)

def true_section(alpha_deg: np.ndarray, x_ref_over_c: float = 0.25) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """The invented 'true' section used to make the synthetic data (given; not a measurement):
    cl = a0 (alpha - alpha_L0), cd = cd0 + k (cl - cl*)^2, and a pitching moment that is constant,
    cm_ac, about x_ac; returned cm is about x_ref_over_c on the chord line (transferred with cn)."""
    p = TRUE_SECTION
    a = np.radians(np.asarray(alpha_deg, float))
    cl = p["a0"] * (a - np.radians(p["alpha_L0_deg"]))
    cd = p["cd0"] + p["k_cd"] * (cl - p["cl_cd_min"]) ** 2
    cn = cl * np.cos(a) + cd * np.sin(a)
    cm = p["cm_ac"] + (x_ref_over_c - p["x_ac"]) * cn
    return cl, cd, cm


def synthetic_balance_data(seed: int = 42, alpha_deg: np.ndarray | None = None, repeats: int = 2) -> dict:
    """SYNTHETIC balance readings for the section model in MODEL (given).

    Builds the true loads from ``true_section`` at the nominal tunnel q, converts them to what the
    balance would read (N, A in chord axes; moment about the balance centre), then adds independent
    Gaussian noise with the NOISE standard deviations to N, A, M, to the recorded q (relative) and to
    the recorded alpha. The recorded alpha and q are what the learner sees; the true ones are
    returned too, for checking. Seeded with numpy.random.default_rng(seed). Not a measurement."""
    if alpha_deg is None:
        alpha_deg = np.arange(-4.0, 10.01, 1.0)
    rng = np.random.default_rng(seed)
    m = MODEL
    S, c = m["c"] * m["span"], m["c"]
    q_true = float(dynamic_pressure(RHO_ROOM, m["V"]))
    a_true = np.repeat(np.asarray(alpha_deg, float), repeats)
    run = np.tile(np.arange(1, repeats + 1), len(alpha_deg))
    cl, cd, cm_qc = true_section(a_true, 0.25)
    L, D = cl * q_true * S, cd * q_true * S
    N, A = wind_to_body(L, D, np.radians(a_true))
    M_qc = cm_qc * q_true * S * c
    M_bal = transfer_moment(M_qc, N, A, 0.25 * c, 0.0, m["x_bal"], m["z_bal"])
    n = a_true.size
    return dict(
        run=run,
        alpha_deg=a_true + NOISE["alpha_deg"] * rng.standard_normal(n),
        q=q_true * (1.0 + NOISE["q_rel"] * rng.standard_normal(n)),
        N=N + NOISE["N"] * rng.standard_normal(n),
        A=A + NOISE["A"] * rng.standard_normal(n),
        M_bal=M_bal + NOISE["M"] * rng.standard_normal(n),
        alpha_true_deg=a_true, q_true=np.full(n, q_true), N_true=N, A_true=A, M_bal_true=M_bal,
        seed=seed,
    )


def write_balance_csv(data: dict, path: str | Path) -> None:
    """Write synthetic readings to CSV with the SYNTHETIC provenance header (given)."""
    m = MODEL
    header = [
        f"# SYNTHETIC: generated by synthetic_balance_data(seed={data['seed']}) in solution.py; not a measurement",
        "# Invented 2-D NACA 2412 section model in an imaginary tunnel: chord 0.30 m, span 0.60 m, V = 35 m/s, room air",
        f"# Balance moment centre at x = {m['x_bal']} m aft of the leading edge, z = {m['z_bal']} m (below the chord line)",
        "# Columns: run, alpha_deg (recorded), q_Pa (recorded), N_N (normal force, chord axes, up), "
        "A_N (axial force, aft), M_bal_Nm (nose-up, about the balance centre)",
        f"# Noise (1 sigma): N {NOISE['N']} N, A {NOISE['A']} N, M {NOISE['M']} N m, q {NOISE['q_rel'] * 100} %, "
        f"alpha {NOISE['alpha_deg']} deg",
        "run,alpha_deg,q_Pa,N_N,A_N,M_bal_Nm",
    ]
    rows = [f"{r:d},{a:.6f},{q:.6f},{N:.6f},{A:.6f},{M:.6f}"
            for r, a, q, N, A, M in zip(data["run"], data["alpha_deg"], data["q"], data["N"], data["A"], data["M_bal"])]
    Path(path).write_text("\n".join(header + rows) + "\n")


def read_balance_csv(path: str | Path) -> dict:
    """Read a balance CSV written by ``write_balance_csv`` (given). '#' lines are provenance."""
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    cols = lines[0].split(",")
    arr = np.array([[float(v) for v in ln.split(",")] for ln in lines[1:]])
    d = {name: arr[:, i] for i, name in enumerate(cols)}
    return dict(run=d["run"].astype(int), alpha_deg=d["alpha_deg"], q=d["q_Pa"], N=d["N_N"], A=d["A_N"],
                M_bal=d["M_bal_Nm"])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save balance_coefficients.png and the synthetic CSV."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    t = TRAINER
    W = t["mass"] * G0
    q_cr = float(dynamic_pressure(RHO_ISA, t["V_cruise"]))
    q_ap = float(dynamic_pressure(RHO_ISA, t["V_approach"]))
    print(f"trainer: W = {W:.1f} N, cruise q = {q_cr:.1f} Pa, CL_cruise = W/(qS) = {W / (q_cr * t['S']):.4f}; "
          f"approach q = {q_ap:.1f} Pa, CL_approach = {W / (q_ap * t['S']):.4f}")
    Re_tr = float(reynolds_number(t["V_cruise"], t["c"]))
    print(f"trainer cruise: Re_c = {Re_tr:.4g}, M = {float(mach_number(t['V_cruise'])):.4f}")

    m = MODEL
    S, c = m["c"] * m["span"], m["c"]
    nu_room = MU_ROOM / RHO_ROOM
    q_t = float(dynamic_pressure(RHO_ROOM, m["V"]))
    Re_m = float(reynolds_number(m["V"], c, nu_room))
    print(f"SYNTHETIC tunnel model: c = {c} m, span {m['span']} m, S = {S:.3f} m^2, V = {m['V']} m/s, q = {q_t:.1f} Pa, "
          f"Re_c = {Re_m:.4g}, M = {m['V'] / A_ROOM:.4f}")
    V_match = Re_tr * nu_room / c
    print(f"  to match the trainer's Re with this model: V = {V_match:.1f} m/s, M = {V_match / A_ROOM:.3f}")

    data = synthetic_balance_data(42)
    write_balance_csv(data, out_dir / "SYNTHETIC_seed042_balance.csv")
    r = reduce_balance(data["N"], data["A"], data["M_bal"], data["alpha_deg"], data["q"], S, c, m["x_bal"], m["z_bal"])
    x_ac, cm_ac, k = fit_aerodynamic_centre(r["cn"], r["cm"], 0.25 * c, c)
    a = data["alpha_deg"]
    slope, icpt = np.polyfit(np.radians(a), r["cl"], 1)
    print(f"\nSYNTHETIC seed 42 ({a.size} points, alpha -4..10 deg): lift slope {slope:.3f} /rad, "
          f"alpha_L0 = {np.degrees(-icpt / slope):.3f} deg (true {TRUE_SECTION['a0']}, {TRUE_SECTION['alpha_L0_deg']})")
    cl_need = W / (q_cr * t["S"])
    print(f"section alpha for the trainer's cruise cl {cl_need:.4f} from the fit: "
          f"{np.degrees((cl_need - icpt) / slope):.2f} deg (2-D; finite-wing effects come in P10.09-P10.10)")
    (p, cov) = np.polyfit(r["cn"], r["cm"], 1, cov=True)
    u_k = np.sqrt(cov[0, 0])
    print(f"aerodynamic centre: x_ac/c = {x_ac / c:.4f} +/- {u_k:.4f} (1 sigma, fit), cm_ac = {cm_ac:.4f} "
          f"(true {TRUE_SECTION['x_ac']}, {TRUE_SECTION['cm_ac']})")
    # wrong reductions
    r_noz = reduce_balance(data["N"], data["A"], data["M_bal"], data["alpha_deg"], data["q"], S, c, m["x_bal"], 0.0)
    x_noz, _, _ = fit_aerodynamic_centre(r_noz["cn"], r_noz["cm"], 0.25 * c, c)
    cl_noaxis = data["N"] / (data["q"] * S)
    cd_noaxis = data["A"] / (data["q"] * S)
    i4 = np.where(np.isclose(data["alpha_true_deg"], 4.0))[0][0]
    print(f"ignoring the balance z offset: x_ac/c = {x_noz / c:.4f}; at alpha = 4 deg cm error "
          f"{r_noz['cm'][i4] - r['cm'][i4]:+.5f}")
    print(f"skipping the axis rotation at alpha = 4 deg: cl {cl_noaxis[i4]:.4f} vs {r['cl'][i4]:.4f}; "
          f"cd {cd_noaxis[i4]:.5f} vs {r['cd'][i4]:.5f}")
    # moment reference shift
    print("cm-vs-cn slope for different reference points (expected (x_ref - x_ac)/c):")
    for xr in (0.0, 0.25, 0.40, 0.50):
        rr = reduce_balance(data["N"], data["A"], data["M_bal"], data["alpha_deg"], data["q"], S, c, m["x_bal"],
                            m["z_bal"], xr * c)
        _, _, kk = fit_aerodynamic_centre(rr["cn"], rr["cm"], xr * c, c)
        print(f"  x_ref/c = {xr:.2f}: slope {kk:+.4f} (expected {xr - x_ac / c:+.4f})")

    # uncertainty at alpha = 4 deg, run 1
    x0 = np.array([data["N"][i4], data["A"][i4], data["M_bal"][i4], data["alpha_deg"][i4], data["q"][i4]])
    u0 = np.array([NOISE["N"], NOISE["A"], NOISE["M"], NOISE["alpha_deg"], NOISE["q_rel"] * data["q"][i4]])
    f = lambda xx: point_coefficients(xx, S, c, m["x_bal"], m["z_bal"], 0.25 * c)
    y0 = f(x0)
    uy, contrib = propagate_uncertainty(f, x0, u0)
    mc_mean, mc_std = monte_carlo_uncertainty(f, x0, u0, n=20000, seed=1)
    names = ["N", "A", "M", "alpha", "q"]
    print(f"\nalpha = 4 deg (run 1): cl = {y0[0]:.4f} +/- {2 * uy[0]:.4f}, cd = {y0[1]:.5f} +/- {2 * uy[1]:.5f}, "
          f"cm_c/4 = {y0[2]:.4f} +/- {2 * uy[2]:.4f}  (k = 2)")
    print(f"  Monte Carlo (20000, seed 1) 1-sigma: cl {mc_std[0]:.5f}, cd {mc_std[1]:.6f}, cm {mc_std[2]:.6f}; "
          f"linear: {uy[0]:.5f}, {uy[1]:.6f}, {uy[2]:.6f}")
    for i, out in enumerate(["cl", "cd", "cm"]):
        share = contrib[i] / contrib[i].sum() * 100
        print(f"  variance budget {out}: " + ", ".join(f"{n} {s:.1f} %" for n, s in zip(names, share)))
    print(f"  alpha uncertainty that equalises the alpha and A shares of u(cd) at 4 deg: "
          f"{NOISE['alpha_deg'] * np.sqrt(contrib[1][1] / contrib[1][3]):.3f} deg")
    cl_true, cd_true, cm_true = true_section(np.array([4.0]))
    print(f"  true values (synthetic generator): cl {cl_true[0]:.4f}, cd {cd_true[0]:.5f}, cm {cm_true[0]:.4f}")
    dcd = f(x0 + np.array([0, 0, 0, 0.1, 0]))[1] - y0[1]
    print(f"  an alpha error of +0.1 deg changes cd by {dcd:+.5f} ({dcd / y0[1] * 100:+.1f} %); cn * 0.1 deg = "
          f"{r['cn'][i4] * np.radians(0.1):.5f}")

    # alpha sweep of the uncertainty of cd
    u_cd = []
    for i in range(a.size):
        xi = np.array([data["N"][i], data["A"][i], data["M_bal"][i], data["alpha_deg"][i], data["q"][i]])
        ui = np.array([NOISE["N"], NOISE["A"], NOISE["M"], NOISE["alpha_deg"], NOISE["q_rel"] * data["q"][i]])
        u_cd.append(propagate_uncertainty(f, xi, ui)[0])
    u_cd = np.array(u_cd)

    fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
    ax[0].errorbar(a, r["cl"], yerr=2 * u_cd[:, 0], fmt="o", ms=3, label="SYNTHETIC balance, k = 2")
    aa = np.linspace(-4, 10, 50)
    ax[0].plot(aa, true_section(aa)[0], "k--", lw=0.8, label="synthetic generator (truth)")
    ax[0].set_xlabel("alpha (deg)"); ax[0].set_ylabel("cl"); ax[0].set_title("SYNTHETIC data: lift"); ax[0].legend()
    ax[1].errorbar(r["cl"], r["cd"], yerr=2 * u_cd[:, 1], fmt="o", ms=3, label="SYNTHETIC balance, k = 2")
    ax[1].plot(true_section(aa)[0], true_section(aa)[1], "k--", lw=0.8, label="truth")
    ax[1].set_xlabel("cl"); ax[1].set_ylabel("cd"); ax[1].set_title("SYNTHETIC data: drag polar"); ax[1].legend()
    for xr in (0.0, 0.25, 0.5):
        rr = reduce_balance(data["N"], data["A"], data["M_bal"], data["alpha_deg"], data["q"], S, c, m["x_bal"],
                            m["z_bal"], xr * c)
        ax[2].plot(rr["cn"], rr["cm"], "o", ms=3, label=f"about x/c = {xr}")
    ax[2].set_xlabel("cn"); ax[2].set_ylabel("cm"); ax[2].set_title("SYNTHETIC data: moment reference shift")
    ax[2].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "balance_coefficients.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'balance_coefficients.png'} and {out_dir / 'SYNTHETIC_seed042_balance.csv'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
