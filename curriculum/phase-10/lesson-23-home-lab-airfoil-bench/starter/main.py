"""P10.23 Home Lab: Airfoil Lift & Drag on a Bench Rig: starter code.

A guarded desk fan with a drinking-straw flow straightener blows on a 100 mm chord foam NACA 0012
section between end plates. A kitchen scale under the stand reads the lift (the reading falls when
the wing lifts); a pivot arm (bell-crank) turns the drag into a push on a second, finer scale. An
inclined water manometer on a pitot tube gives the dynamic pressure. This module turns those
readings into lift and drag coefficients with uncertainties and compares them with theory,
computation and published data:

* room air density and viscosity from temperature, pressure and humidity;
* dynamic pressure from the inclined-tube manometer (q = (rho_w - rho) g dl sin(theta)), or from a
  vane anemometer, and the airspeed;
* balance reduction: tares (fan on with the sting alone; the model at zero lift), the pivot-arm
  ratio from a calibration with known horizontal loads, a spring balance as an alternative;
* CL and CD with first-order uncertainty propagation (P00.14), the Reynolds number, and a blockage
  estimate carried as an uncertainty on q;
* predictions: thin-airfoil and compact panel cl, a finite-span correction for the end-plated model,
  and a low-Reynolds-number drag estimate; and the comparison with Ladson's NACA 0012 data (NASA
  TM-4074) at a Reynolds number about 70 times higher;
* a seeded SYNTHETIC dataset generator, make_synthetic_bench_data(seed=...), for learners without a
  kit. Synthetic data are not measurements; every file, table and plot made from them says SYNTHETIC.

Units: SI inside the functions; scale readings in grams and manometer readings in millimetres as
they are read from the instruments.

Run:  python main.py [out_dir]     (prediction, SYNTHETIC reduction, comparisons; writes airfoil_bench.png)
Test: pytest curriculum/phase-10/lesson-23-home-lab-airfoil-bench
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

G = 9.81                    # m/s^2
RHO_W = 998.0               # kg/m^3, water at about 20 C
R_DRY = 287.05              # J/(kg K)
R_VAP = 461.5               # J/(kg K)
RHO_ROOM = 1.204            # kg/m^3, the brief's room air (20 C, 101.3 kPa)
MU_ROOM = 1.81e-5           # Pa s
# The lesson's model and rig (nominal; measure your own)
CHORD = 0.100               # m
SPAN = 0.200                # m between the end plates
THICK = 0.012               # m (NACA 0012 at 100 mm chord)
JET_AREA = 0.30 * 0.30      # m^2, straightener exit
ENDPLATE_H = 0.12           # m, end-plate height (chordwise length 0.15 m)
ARM_NOMINAL = 2.0           # l_D / l_s of the bell-crank (120 mm / 60 mm)
SYNTH_SEED = 230


# ---------------------------------------------------------------- given: data location

def find_repo_root(start: str | Path | None = None) -> Path:
    """Given: the repository root (the folder that contains curriculum/phase-10/data), searched upward
    from `start` (default this file), so the code runs from curriculum/... and from solutions/...."""
    p = Path(start if start is not None else __file__).resolve()
    for parent in [p] + list(p.parents):
        if (parent / "curriculum" / "phase-10" / "data" / "ladson_naca0012_M015.csv").exists():
            return parent
    raise FileNotFoundError("repository root not found above " + str(p))


def load_ladson(table: str = "I", reynolds: float = 2.0e6) -> dict:
    """Given: one Reynolds-number block of Ladson (NASA TM-4074), as published (alpha_deg, cl, cd, cm)."""
    path = find_repo_root() / "curriculum" / "phase-10" / "data" / "ladson_naca0012_M015.csv"
    rows = list(csv.DictReader(l for l in path.open() if not l.startswith("#")))
    sel = [r for r in rows if r["table"] == table and abs(float(r["reynolds"]) - reynolds) < 0.01 * reynolds]
    f = lambda k: np.array([float(r[k]) if r[k] != "" else np.nan for r in sel])
    return {"alpha_deg": f("alpha_deg"), "cl": f("cl"), "cd": f("cd"), "reynolds": float(sel[0]["reynolds"])}


# ---------------------------------------------------------------- given: compact panel method and boundary layer

def naca4(code: str, n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    """Given (from P10.02): symmetric or cambered NACA 4-digit section, closed trailing edge, cosine
    spacing, clockwise from the trailing edge (lower surface first), chord 1."""
    m, p, t = int(code[0]) / 100, int(code[1]) / 10, int(code[2:]) / 100
    x = 0.5 * (1 - np.cos(np.linspace(0, np.pi, n + 1)))
    yt = 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2 + 0.2843 * x ** 3 - 0.1036 * x ** 4)
    yc, dyc = np.zeros_like(x), np.zeros_like(x)
    if m > 0:
        f = x < p
        yc[f] = m / p ** 2 * (2 * p * x[f] - x[f] ** 2); dyc[f] = 2 * m / p ** 2 * (p - x[f])
        yc[~f] = m / (1 - p) ** 2 * (1 - 2 * p + 2 * p * x[~f] - x[~f] ** 2); dyc[~f] = 2 * m / (1 - p) ** 2 * (p - x[~f])
    th = np.arctan(dyc)
    xu, yu, xl, yl = x - yt * np.sin(th), yc + yt * np.cos(th), x + yt * np.sin(th), yc - yt * np.cos(th)
    return np.concatenate([xl[::-1], xu[1:]]), np.concatenate([yl[::-1], yu[1:]])


def _panel(code: str, alpha_deg, n: int = 80) -> dict:
    """Given (compact version of P10.05): linear-strength vortex panel method (Kuethe-Chow form, Kutta
    condition gamma_1 + gamma_N+1 = 0) for a NACA 4-digit section at one or more angles [deg].
    Returns cl per angle and, per angle, the tangential surface velocity Vt/V at the panel midpoints
    (positive along the clockwise panel direction), with the midpoints xc and panel lengths S."""
    X, Y = naca4(code, n)
    N = len(X) - 1
    xc, yc = 0.5 * (X[:-1] + X[1:]), 0.5 * (Y[:-1] + Y[1:])
    S = np.hypot(np.diff(X), np.diff(Y))
    th = np.arctan2(np.diff(Y), np.diff(X))
    xi, eta = xc[:, None] - X[None, :-1], yc[:, None] - Y[None, :-1]
    thi, thj, Sj = th[:, None], th[None, :], S[None, :]
    A = -xi * np.cos(thj) - eta * np.sin(thj)
    Bq = xi ** 2 + eta ** 2
    Cc, D = np.sin(thi - thj), np.cos(thi - thj)
    E = xi * np.sin(thj) - eta * np.cos(thj)
    F = np.log(1 + Sj * (Sj + 2 * A) / Bq)
    Gm = np.arctan2(E * Sj, Bq + A * Sj)
    P = xi * np.sin(thi - 2 * thj) + eta * np.cos(thi - 2 * thj)
    Q = xi * np.cos(thi - 2 * thj) - eta * np.sin(thi - 2 * thj)
    Cn2 = D + 0.5 * Q * F / Sj - (A * Cc + D * E) * Gm / Sj
    Cn1 = 0.5 * D * F + Cc * Gm - Cn2
    Ct2 = Cc + 0.5 * P * F / Sj + (A * D - Cc * E) * Gm / Sj
    Ct1 = 0.5 * Cc * F - D * Gm - Ct2
    i = np.arange(N)
    Cn1[i, i], Cn2[i, i], Ct1[i, i], Ct2[i, i] = -1.0, 1.0, 0.5 * np.pi, 0.5 * np.pi
    An = np.zeros((N + 1, N + 1)); At = np.zeros((N, N + 1))
    An[:N, :N] += Cn1; An[:N, 1:] += Cn2; At[:, :N] += Ct1; At[:, 1:] += Ct2
    An[N, 0] = An[N, N] = 1.0
    al = np.radians(np.atleast_1d(np.asarray(alpha_deg, dtype=float)))
    rhs = np.vstack([np.sin(th[:, None] - al[None, :]), np.zeros((1, len(al)))])
    g = np.linalg.solve(An, rhs)
    Vt = np.cos(th[:, None] - al[None, :]) + At @ g
    Gamma = 2 * np.pi * np.sum(0.5 * (g[:-1] + g[1:]) * S[:, None], axis=0)
    return {"cl": 2 * Gamma / (X.max() - X.min()), "Vt": Vt, "xc": xc, "S": S}


def panel_cl(code: str, alpha_deg, n: int = 80) -> np.ndarray:
    """Given: inviscid panel-method cl of a NACA 4-digit section at each angle [deg]."""
    return _panel(code, alpha_deg, n)["cl"]


def upper_surface_velocity(code: str, alpha_deg: float, n: int = 80) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Given: the panel method's upper-surface edge velocity Ue/V against arc length s/c from the
    stagnation point (s = 0, Ue = 0 prepended) and the chordwise position x/c of each station."""
    r = _panel(code, alpha_deg, n)
    Vt, S, xc = r["Vt"][:, 0], r["S"], r["xc"]
    i0 = int(np.where((Vt[:-1] < 0) & (Vt[1:] >= 0))[0][0])
    sc = np.concatenate([[0.0], np.cumsum(0.5 * (S[:-1] + S[1:]))])
    s_stag = sc[i0] + (-Vt[i0] / (Vt[i0 + 1] - Vt[i0])) * (sc[i0 + 1] - sc[i0])
    idx = np.arange(i0 + 1, len(Vt))
    return (np.concatenate([[0.0], sc[idx] - s_stag]), np.concatenate([[0.0], Vt[idx]]),
            np.concatenate([[xc[idx[0]]], xc[idx]]))


# ---------------------------------------------------------------- the learner's functions: air and airspeed

def air_density(T_C: float, p_Pa: float, RH: float = 0.0) -> float:
    """Moist-air density: dry-air and water-vapour partial densities, with the saturation pressure
    from the Magnus-Tetens formula p_sat = 610.94 exp(17.625 T / (T + 243.04)) Pa (T in C)."""
    raise NotImplementedError


def air_viscosity(T_C: float) -> float:
    """Sutherland's law for air: mu = 1.716e-5 (T/273.15)^1.5 (273.15 + 110.4)/(T + 110.4) Pa s."""
    raise NotImplementedError


def manometer_q(dl_mm: float, incline_deg: float, rho_air: float, u_dl_mm: float = 0.0,
                u_incline_deg: float = 0.0, rho_liquid: float = RHO_W) -> tuple[float, float]:
    """Dynamic pressure from an inclined-tube manometer on a pitot-static probe (P08.19):
    q = (rho_liquid - rho_air) g dl sin(theta), with dl the reading change along the tube.
    Returns q [Pa] and its standard uncertainty from dl and theta."""
    raise NotImplementedError


def airspeed(q: float, rho: float, u_q: float = 0.0, u_rho: float = 0.0) -> tuple[float, float]:
    """V = sqrt(2 q / rho) and its standard uncertainty: u_V/V = (1/2) sqrt((u_q/q)^2 + (u_rho/rho)^2)."""
    raise NotImplementedError


def anemometer_q(V: float, rho: float, rel_spec: float = 0.03, abs_spec: float = 0.1) -> tuple[float, float]:
    """q = rho V^2 / 2 from a vane anemometer whose accuracy is +/-(rel_spec V + abs_spec); the
    specification is read as the half-width of a uniform band, so u_V = (rel_spec V + abs_spec)/sqrt 3."""
    raise NotImplementedError


def reynolds_number(V: float, c: float, rho: float = RHO_ROOM, mu: float = MU_ROOM) -> float:
    """Re = rho V c / mu."""
    raise NotImplementedError


# ---------------------------------------------------------------- the learner's functions: balance

def fit_line(x: np.ndarray, y: np.ndarray) -> dict:
    """Ordinary least squares y = m x + b with standard uncertainties from the residual scatter."""
    raise NotImplementedError


def pivot_calibration(applied_g: np.ndarray, reading_g: np.ndarray) -> tuple[float, float]:
    """Arm ratio R = l_D / l_s of the drag bell-crank from known horizontal loads hung over a pulley at
    the model's drag point: the slope of scale reading against applied load. Returns R and u_R."""
    raise NotImplementedError


def spring_calibration(mass_g: np.ndarray, extension_mm: np.ndarray) -> tuple[float, float]:
    """Spring rate k [N/m] and u_k from hung masses and measured extensions (slope of force on extension)."""
    raise NotImplementedError


def lift_from_scale(off_g, on_g, tare_g: float = 0.0) -> np.ndarray:
    """Lift [N] from the lift scale: the reading FALLS by the lift. tare_g is the change (on - off) of
    the tare run (fan on, sting without model), which is removed: L = -((on - off) - tare) g."""
    raise NotImplementedError


def drag_from_pivot(off_g, on_g, ratio: float, tare_g: float = 0.0) -> np.ndarray:
    """Drag [N] from the bell-crank scale: the reading RISES by D R (R = l_D / l_s); the tare change
    (sting alone, fan on) is removed before dividing by R."""
    raise NotImplementedError


def drag_from_spring(x_off_mm, x_on_mm, k: float) -> np.ndarray:
    """Drag [N] from a calibrated spring: D = k (x_on - x_off)."""
    raise NotImplementedError


def coefficient(F: float, u_F: float, q: float, u_q: float, S: float, u_S: float = 0.0) -> tuple[float, float]:
    """C = F / (q S) and its standard uncertainty by first-order propagation:
    u_C^2 = (u_F/(q S))^2 + (C u_q/q)^2 + (C u_S/S)^2 (valid also when F is near zero)."""
    raise NotImplementedError


def blockage_fraction(alpha_deg: float, c: float = CHORD, s: float = SPAN, t: float = THICK,
                      endplate_t: float = 0.003, endplate_h: float = ENDPLATE_H, jet_area: float = JET_AREA) -> float:
    """Rough blockage estimate epsilon = A_frontal / (4 A_jet) (the simple solid-blockage rule for a
    closed section), with the frontal area of the model at alpha and of its end plates edge-on. In an
    open jet the true correction is smaller and can change sign, so the lesson carries epsilon as an
    uncertainty on the velocity (2 epsilon on q), not as a correction."""
    raise NotImplementedError


def en_number(x1: float, u1: float, x2: float, u2: float, k: float = 2.0) -> float:
    """Normalised error |x1 - x2| / sqrt(U1^2 + U2^2) with expanded uncertainties U = k u (P00.15)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the learner's functions: predictions

def thin_airfoil_cl(alpha_deg) -> np.ndarray:
    """Thin-airfoil theory for a symmetric section: cl = 2 pi alpha."""
    raise NotImplementedError


def finite_span_slope(a0: float, AR_eff: float, e: float = 1.0) -> float:
    """Lift-curve slope of a finite wing from the section slope a0 (lifting-line form,
    a = a0 / (1 + a0 / (pi e AR_eff))). An end-plated model has an effective aspect ratio above its
    geometric one; how much above is uncertain, so the lesson brackets it."""
    raise NotImplementedError


def low_re_cd(Re: float, cl, alpha_deg, t_over_c: float = 0.12, suction_loss: float = 0.5,
              AR_eff: float | None = None, e: float = 0.9) -> np.ndarray:
    """A low-Reynolds-number drag estimate: laminar (Blasius) skin friction on both surfaces times a
    thickness form factor, 2 (1.328/sqrt(Re)) (1 + 2 t/c + 60 (t/c)^4); plus the pressure drag of a lost
    fraction of leading-edge suction, suction_loss * cl tan(alpha) (0 = full suction, as in attached
    potential flow; 1 = none, as on a separated flat plate); plus induced drag cl^2/(pi e AR_eff) if
    AR_eff is given. It ignores laminar separation bubbles, so it is a lower-bound style estimate."""
    raise NotImplementedError


def laminar_upper_surface(code: str, alpha_deg: float, Re: float, n: int = 80) -> dict:
    """With a boundary-layer estimate (P08.11, P08.12 recaps in the lesson): Thwaites' method on the
    panel method's upper surface, theta^2 = 0.45 (nu/V c) Ue^-6 int Ue^5 ds, lambda = theta^2 Re dUe/ds.
    Returns the laminar-separation station x_sep/c (first lambda < -0.09), the station where Michel's
    criterion Re_theta > 1.174 (1 + 22400/Re_x) Re_x^0.46 is first met (x_tr/c, or None), Re_theta at
    separation, and whether transition comes first. lambda, hence x_sep, does not depend on Re."""
    raise NotImplementedError


def predict_bench(alpha_deg, V: float, rho: float = RHO_ROOM, mu: float = MU_ROOM, c: float = CHORD,
                  s: float = SPAN, AR_eff: tuple = (2 * SPAN / CHORD, 3 * SPAN / CHORD),
                  ratio: float = ARM_NOMINAL) -> dict:
    """The prediction to write down before measuring: Re; 2-D thin-airfoil and panel cl; the finite-span
    CL band for AR_eff between AR_eff[0] and AR_eff[1] (the end plates' effect is uncertain); the
    CD band from low_re_cd with suction_loss 0 to 1; and the scale changes these imply [g]."""
    raise NotImplementedError


# ---------------------------------------------------------------- SYNTHETIC data (given)

def _synthetic_truth_CL(alpha_true_deg):
    """Generator assumption (not a measured fact): a low-Re, end-plated model with slope 0.060 /deg,
    a soft stall to CL_max ~ 0.6 near 9-10 deg and a gentle fall beyond 11 deg."""
    a = np.asarray(alpha_true_deg, float)
    lin = 0.060 * a
    CLmax = 0.62
    cl = lin / (1 + (np.abs(lin) / CLmax) ** 4) ** 0.25
    return cl * (1 - 0.03 * np.clip(np.abs(a) - 11, 0, None))


def _synthetic_truth_CD(alpha_true_deg):
    """Generator assumption (not a measured fact): CD = 0.035 + 0.10 CL^2 + 6e-4 max(|alpha| - 7, 0)^2."""
    a = np.asarray(alpha_true_deg, float)
    cl = _synthetic_truth_CL(a)
    return 0.035 + 0.10 * cl ** 2 + 6e-4 * np.clip(np.abs(a) - 7, 0, None) ** 2


def make_synthetic_bench_data(seed: int = SYNTH_SEED) -> dict:
    """SYNTHETIC readings of the lesson's rig (foam NACA 0012, c = 100 mm, s = 200 mm, end plates;
    fan about 4 m/s; lift on a 0.1 g scale; drag on a 0.01 g scale via a bell-crank). NOT a measurement.

    The generator's assumptions (stated, not measured): CL and CD from _synthetic_truth_CL/_CD
    (reduced lift slope, CD 0.035-0.1, early soft stall), a flow angularity of +0.6 deg in the jet,
    a fan speed that drifts from 4.05 to 3.97 m/s with 1.5 % run-to-run unsteadiness, a sting with
    0.005 N of drag and no lift, a true arm ratio of 2.04 against the nominal 2.0, a true manometer
    incline of 10.15 deg recorded as 10.0, reading noise of 0.3 mm on the manometer, 0.12 g on the
    lift scale and 0.04 g on the drag scale, and room air at 21.0 C, 100.8 kPa, 45 % RH.
    Seeded with numpy.random.default_rng(seed). Returns a dict of tables (lists of row tuples)."""
    rng = np.random.default_rng(seed)
    T_C, p, RH = 21.0, 100.8e3, 0.45
    p_v = RH * 610.94 * np.exp(17.625 * T_C / (T_C + 243.04))       # same physics as air_density, kept
    rho = (p - p_v) / (R_DRY * (T_C + 273.15)) + p_v / (R_VAP * (T_C + 273.15))   # self-contained here
    flow_angle, R_true, incl_true, D_sting = 0.6, 2.04, 10.15, 0.005
    S = CHORD * SPAN
    base_lift, base_drag = 352.6, 18.40

    def rnd(x, r):
        return np.round(np.asarray(x) / r) * r

    conditions = [("temperature", 21.2, "C", 0.5), ("pressure", 100.7e3, "Pa", 300.0),
                  ("relative_humidity", 0.45, "1", 0.05), ("chord", 0.1002, "m", 0.0005),
                  ("span", 0.1995, "m", 0.0005), ("incline", 10.0, "deg", 0.2),
                  ("lift_scale_resolution", 0.1, "g", 0.0), ("drag_scale_resolution", 0.01, "g", 0.0)]
    # manometer: 3 readings before and 3 after the sweep
    airspeed_rows = []
    for run, V in enumerate([4.05, 4.04, 4.06, 3.98, 3.97, 3.98], start=1):
        q = 0.5 * rho * V ** 2
        dl = q / ((RHO_W - rho) * G * np.sin(np.radians(incl_true))) * 1000
        z = 50.0 + rng.normal(0, 0.3)
        airspeed_rows.append((run, 10.0, float(rnd(z, 0.5)), float(rnd(z + dl + rng.normal(0, 0.3), 0.5))))
    # pivot calibration (horizontal loads over a pulley) and a lift-scale check
    cal_rows = []
    for m in (0.0, 1.0, 2.0, 5.0, 10.0, 5.0, 2.0, 0.0):
        cal_rows.append(("drag_pivot", m, float(rnd(base_drag + R_true * m + rng.normal(0, 0.03), 0.01))))
    for m in (0.0, 20.0, 50.0):
        cal_rows.append(("lift_scale_check", m, float(rnd(base_lift + m + rng.normal(0, 0.05), 0.1))))

    def reading_pair(alpha_set, model=True, V=4.0):
        qq = 0.5 * rho * (V * (1 + rng.normal(0, 0.015))) ** 2
        at = alpha_set + flow_angle
        L = _synthetic_truth_CL(at) * qq * S if model else 0.0
        D = (_synthetic_truth_CD(at) * qq * S if model else 0.0) + D_sting * qq / (0.5 * rho * 4.0 ** 2)
        lo = base_lift + rng.normal(0, 0.05)
        do = base_drag + rng.normal(0, 0.02)
        lon = lo - L / G * 1000 + rng.normal(0, 0.12)
        don = do + R_true * D / G * 1000 + rng.normal(0, 0.04)
        return float(rnd(lo, 0.1)), float(rnd(lon, 0.1)), float(rnd(do, 0.01)), float(rnd(don, 0.01))

    alphas = [-4, -2, 0, 2, 4, 6, 8, 10, 12, 14]
    sweep_rows = []
    for run, seq in ((1, alphas), (2, alphas[::-1])):
        for k, a in enumerate(seq):
            V = 4.02 - 0.04 * (k / (len(seq) - 1) + (run - 1))
            sweep_rows.append((run, a) + reading_pair(a, True, V))
    tare_rows = [("sting_only",) + reading_pair(0, False) for _ in range(3)]
    tare_rows += [("model_alpha0",) + reading_pair(0, True) for _ in range(3)]
    return {"conditions": conditions, "airspeed": airspeed_rows, "calibration": cal_rows, "sweep": sweep_rows,
            "tare": tare_rows,
            "truth": {"rho": rho, "flow_angle_deg": flow_angle, "ratio": R_true, "incline_deg": incl_true,
                      "sting_drag_N_at_4ms": D_sting, "CL": _synthetic_truth_CL, "CD": _synthetic_truth_CD}}


SYNTH_COLUMNS = {
    "conditions": "quantity,value,unit,standard_uncertainty",
    "airspeed": "run,incline_deg,reading_zero_mm,reading_fan_mm",
    "calibration": "kind,applied_g,reading_g",
    "sweep": "run,alpha_set_deg,lift_off_g,lift_on_g,drag_off_g,drag_on_g",
    "tare": "config,lift_off_g,lift_on_g,drag_off_g,drag_on_g",
}


def write_synthetic_csv(out_dir: str | Path, seed: int = SYNTH_SEED) -> list[Path]:
    """Write the SYNTHETIC dataset as five CSV files named SYNTHETIC_seed<NNN>_<table>.csv whose first
    line says SYNTHETIC."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = make_synthetic_bench_data(seed)
    paths = []
    for key, cols in SYNTH_COLUMNS.items():
        p = out_dir / f"SYNTHETIC_seed{seed:03d}_{key}.csv"
        with open(p, "w", newline="") as fh:
            fh.write(f"# SYNTHETIC: generated by make_synthetic_bench_data(seed={seed}) in solution.py; not a measurement\n")
            fh.write("# P10.23 bench rig: foam NACA 0012, chord 100 mm, span 200 mm between end plates, fan ~4 m/s; "
                     "see projects/p10-airfoil-bench/data/README.md\n")
            fh.write(cols + "\n")
            wr = csv.writer(fh)
            for row in data[key]:
                wr.writerow([f"{v:g}" if isinstance(v, float) else v for v in row])
        paths.append(p)
    return paths


def read_bench_csv(path: str | Path) -> dict:
    """Given: read a bench CSV (lines starting with # are comments) into a dict of columns; numeric
    columns become float arrays, text columns stay lists."""
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(l for l in fh if not l.startswith("#") and l.strip()))
    out = {}
    for k in rows[0]:
        vals = [r[k] for r in rows]
        try:
            out[k] = np.array([float(v) for v in vals])
        except ValueError:
            out[k] = vals
    return out


def load_bench_dataset(folder: str | Path, seed: int = SYNTH_SEED, prefix: str | None = None) -> dict:
    """Given: read the five tables of one dataset from a folder (by default the SYNTHETIC set
    SYNTHETIC_seed<NNN>_*.csv; pass prefix for your own files, e.g. '2026-11-02_recommended')."""
    folder = Path(folder)
    pre = prefix if prefix is not None else f"SYNTHETIC_seed{seed:03d}"
    return {k: read_bench_csv(folder / f"{pre}_{k}.csv") for k in SYNTH_COLUMNS}


# ---------------------------------------------------------------- reduction (the learner's)

def reduce_bench(ds: dict) -> dict:
    """Reduce one dataset (as returned by load_bench_dataset) to CL and CD against alpha with
    uncertainties. Steps: room air; q and V from the manometer (mean of all runs, u from the scatter of
    the mean, the reading resolution and the incline); the arm ratio from the calibration; tares from
    the sting-only runs; lift and drag per reading pair; mean and standard uncertainty at each alpha
    over the repeated runs; CL, CD with u including the blockage band on q; the flow angle from the
    zero-lift intercept of the linear part, and the slope there. The fit's u_slope reflects only the
    scatter about the line; u_slope_total adds the manometer's uncertainty of q, which scales every CL
    together and so cannot show up in the scatter."""
    raise NotImplementedError


def compare_with_ladson(alpha_deg, CL, CD, reynolds: float = 2.0e6) -> dict:
    """Interpolate Ladson's published NACA 0012 cl and cd (Table I, free transition) at the bench
    angles and return them with the ratios bench / published (CL over cl, CD over cd)."""
    raise NotImplementedError


# ---------------------------------------------------------------- main

def synthetic_folder() -> Path:
    """Given: projects/p10-airfoil-bench/data/synthetic in the repository."""
    return find_repo_root() / "projects" / "p10-airfoil-bench" / "data" / "synthetic"


def main(out_dir: str | Path = ".") -> None:
    """Prediction, reduction of the SYNTHETIC dataset, comparisons; writes airfoil_bench.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    V0 = 4.0
    q0 = 0.5 * RHO_ROOM * V0 ** 2
    print(f"room air {RHO_ROOM} kg/m^3, {MU_ROOM} Pa s; at {V0} m/s: q = {q0:.2f} Pa = {q0/(RHO_W*G)*1000:.3f} mm of water"
          f" ({q0/(RHO_W*G)*1000/np.sin(np.radians(10)):.2f} mm along a tube inclined at 10 deg)")
    print(f"Re = {reynolds_number(V0, CHORD):.3g} on the bench vs 2.0e6-18.9e6 in Ladson's tunnel ({2e6/reynolds_number(V0, CHORD):.0f} x)")
    qm, uqm = manometer_q(5.6, 10.0, RHO_ROOM, 0.5 / np.sqrt(3) * np.sqrt(2), 0.2)
    print(f"one manometer reading: dl = 5.6 mm at 10.0 +/- 0.2 deg, reading +/- 0.5 mm (rectangular, two readings):"
          f" q = {qm:.2f} +/- {uqm:.2f} Pa ({100*uqm/qm:.1f} %), V = {airspeed(qm, RHO_ROOM, uqm)[0]:.2f} +/- {airspeed(qm, RHO_ROOM, uqm)[1]:.2f} m/s")
    qa, uqa = anemometer_q(V0, RHO_ROOM)
    print(f"vane anemometer +/-(3 % + 0.1 m/s): q = {qa:.2f} +/- {uqa:.2f} Pa ({100*uqa/qa:.1f} %)")
    print(f"blockage estimate A_f/(4 A_jet): {blockage_fraction(0):.4f} at 0 deg, {blockage_fraction(10):.4f} at 10 deg")

    print("\nPREDICTION (write it down and date it before measuring), V = 4.0 m/s, AR_eff 4-6 (end plates), suction loss 0-1")
    al = np.array([-4, -2, 0, 2, 4, 6, 8, 10, 12, 14], float)
    pr = predict_bench(al, V0)
    print(" alpha  cl_thin  cl_panel   CL band        CD band          lift change [g]   drag-scale change [g]")
    for i, a in enumerate(al):
        print(f" {a:5.1f}  {pr['cl_thin'][i]:7.3f}  {pr['cl_panel'][i]:7.3f}  {pr['CL_lo'][i]:6.3f} to {pr['CL_hi'][i]:6.3f}"
              f"  {pr['CD_lo'][i]:.4f} to {pr['CD_hi'][i]:.4f}   {pr['lift_g_lo'][i]:6.2f} to {pr['lift_g_hi'][i]:6.2f}"
              f"      {pr['drag_scale_g_lo'][i]:5.2f} to {pr['drag_scale_g_hi'][i]:5.2f}")
    print(f" finite-span slope: {np.radians(finite_span_slope(2*np.pi, 4.0)):.4f} to {np.radians(finite_span_slope(2*np.pi, 6.0)):.4f} /deg "
          f"(2-D thin airfoil {np.radians(2*np.pi):.4f}, panel {pr['cl_panel'][4]/4:.4f})")

    print("\nwith a boundary-layer estimate: Thwaites on the panel upper surface, laminar separation vs Michel transition")
    print(" alpha    x_sep/c    Re 2.6e4: x_tr/c, Re_theta at sep     Re 2.0e6: x_tr/c")
    for a in (0.0, 2.0, 4.0, 6.0, 8.0):
        b1 = laminar_upper_surface("0012", a, reynolds_number(V0, CHORD))
        b2 = laminar_upper_surface("0012", a, 2.0e6)
        f1 = "none" if b1["x_tr"] is None else f"{b1['x_tr']:.3f}"
        print(f" {a:5.1f}    {b1['x_sep']:.3f}      {f1:>6s}, {b1['Re_theta_sep']:6.0f} ({'transition first' if b1['transition_first'] else 'separates laminar'})"
              f"     {b2['x_tr']:.3f} ({'transition first' if b2['transition_first'] else 'separates laminar'})")

    # ---- the SYNTHETIC dataset
    folder = synthetic_folder()
    if (folder / f"SYNTHETIC_seed{SYNTH_SEED:03d}_sweep.csv").exists():
        ds = load_bench_dataset(folder)
        src = str(folder.relative_to(find_repo_root()))
    else:
        import tempfile
        tmp = Path(tempfile.mkdtemp())
        write_synthetic_csv(tmp)
        ds = load_bench_dataset(tmp)
        src = "regenerated in a temporary folder"
    r = reduce_bench(ds)
    print(f"\nSYNTHETIC dataset seed {SYNTH_SEED} ({src}); NOT a measurement")
    print(f"  air {r['rho']:.4f} +/- {r['u_rho']:.4f} kg/m^3; q = {r['q']:.2f} +/- {r['u_q_manometer']:.2f} Pa (manometer);"
          f" V = {r['V']:.2f} +/- {r['u_V']:.2f} m/s; drift during the sweep {100*r['drift']:+.1f} % in q; Re = {r['Re']:.3g}")
    print(f"  arm ratio R = {r['ratio']:.3f} +/- {r['u_ratio']:.3f} (nominal {ARM_NOMINAL}); tares: sting lift {r['tare_lift_g']:+.2f} g,"
          f" sting drag {r['tare_drag_g']:.2f} +/- {r['u_tare_drag_g']:.2f} g on the drag scale")
    print("  alpha_set   CL +/- u        CD +/- u        eps")
    for row in r["rows"]:
        print(f"   {row['alpha_set']:5.1f}   {row['CL']:6.3f} +/- {row['u_CL']:.3f}   {row['CD']:.4f} +/- {row['u_CD']:.4f}   {row['eps']:.4f}")
    print(f"  slope {r['slope_per_deg']:.4f} +/- {r['u_slope']:.4f} (fit) +/- {r['u_slope_total']:.4f} (with q) /deg (-4..4 deg); zero-lift at alpha_set = {r['flow_angle_deg']:+.2f} +/- {r['u_flow_angle_deg']:.2f} deg"
          f" (flow angularity); CL_max = {r['CL_max']:.3f} +/- {r['u_CL_max']:.3f} at {r['alpha_CL_max']:.0f} deg")
    lo, hi = np.radians(finite_span_slope(2 * np.pi, 4.0)), np.radians(finite_span_slope(2 * np.pi, 6.0))
    mid, u_mid = 0.5 * (lo + hi), (hi - lo) / np.sqrt(12)
    print(f"  slope vs predicted finite-span band {lo:.4f}-{hi:.4f}: ratio to mid {r['slope_per_deg']/mid:.3f},"
          f" E_n = {en_number(r['slope_per_deg'], r['u_slope_total'], mid, u_mid):.2f}")
    a_corr = np.array([row["alpha_set"] for row in r["rows"]]) - r["flow_angle_deg"]
    CLs = np.array([row["CL"] for row in r["rows"]]); CDs = np.array([row["CD"] for row in r["rows"]])
    cmp_ = compare_with_ladson(a_corr, CLs, CDs)
    d2 = load_ladson("I", 2.0e6)
    sl_l = np.polyfit(d2["alpha_deg"][:7], d2["cl"][:7], 1)[0]
    print(f"  Ladson Re 2e6 slope {sl_l:.4f} /deg; SYNTHETIC bench slope / Ladson = {r['slope_per_deg']/sl_l:.3f}")
    print("  alpha (flow-angle corrected), CL bench/Ladson, CD bench/Ladson, low-Re cd estimate (suction loss 0.5, AR_eff 5):")
    for a, row, rc, rd in zip(a_corr, r["rows"], cmp_["CL_ratio"], cmp_["CD_ratio"]):
        est = float(low_re_cd(r["Re"], row["CL"], a, AR_eff=5.0))
        print(f"   {a:6.2f}   {rc:6.3f}   {rd:6.2f}   {est:.4f} (SYNTHETIC CD {row['CD']:.4f})")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.6))
    ax[0].errorbar(a_corr, CLs, yerr=2 * np.array([x["u_CL"] for x in r["rows"]]), fmt="o", color="C3",
                   label="SYNTHETIC bench data (k = 2), Re 2.6e4")
    aa = np.linspace(-5, 15, 81)
    ax[0].plot(aa, thin_airfoil_cl(aa), "k--", lw=1, label="thin airfoil 2 pi alpha (2-D)")
    ax[0].plot(aa, panel_cl("0012", aa, 60), "k:", lw=1, label="panel NACA 0012 (2-D, inviscid)")
    ax[0].fill_between(aa, np.radians(aa) * finite_span_slope(2 * np.pi, 4.0), np.radians(aa) * finite_span_slope(2 * np.pi, 6.0),
                       color="0.8", label="finite span, AR_eff 4-6")
    ax[0].plot(d2["alpha_deg"], d2["cl"], "s", color="C0", ms=4, label="Ladson TM-4074, Re 2e6")
    ax[0].set_xlim(-5, 16); ax[0].set_ylim(-0.7, 1.6)
    ax[0].set_xlabel("alpha (flow-angle corrected) [deg]"); ax[0].set_ylabel("CL, cl")
    ax[0].set_title("SYNTHETIC bench lift vs theory, computation and NACA data", fontsize=9)
    ax[0].legend(fontsize=7); ax[0].grid(alpha=0.3)
    ax[1].errorbar(a_corr, CDs, yerr=2 * np.array([x["u_CD"] for x in r["rows"]]), fmt="o", color="C3", label="SYNTHETIC bench data (k = 2)")
    ok = np.isfinite(d2["cd"])
    ax[1].plot(d2["alpha_deg"][ok], d2["cd"][ok], "s", color="C0", ms=4, label="Ladson, Re 2e6")
    ax[1].plot(aa, low_re_cd(r["Re"], np.radians(aa) * finite_span_slope(2 * np.pi, 5.0), aa, AR_eff=5.0), "k--", lw=1,
               label="low-Re estimate (suction loss 0.5, AR_eff 5)")
    ax[1].set_xlabel("alpha [deg]"); ax[1].set_ylabel("CD, cd"); ax[1].set_ylim(0, 0.16)
    ax[1].set_title("SYNTHETIC bench drag", fontsize=9); ax[1].legend(fontsize=7); ax[1].grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "airfoil_bench.png", dpi=120)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
