"""P08.19 Home Lab: Flow & Pressure Measurement: starter code.

Two low-energy experiments with water at room temperature (never more than 2 L) and a small guarded fan
or a hair dryer on its cool setting:

1. Draining bottle (Torricelli). Water leaves a small round hole near the base of a bottle. Bernoulli
   from the free surface to the vena contracta, with a discharge coefficient C_d for the contraction and
   the small losses, and mass conservation for the bottle, give the quasi-steady model
       A dh/dt = -C_d a sqrt(2 g h)   =>   sqrt(h) = sqrt(h0) - k t,   k = C_d (a/A) sqrt(g/2),
   so sqrt(h) is a straight line in time and its slope gives C_d. The prediction of the draining time,
   with C_d from the typical sharp-edged-orifice range, carries an uncertainty budget.
2. U-tube manometer and airspeed. One leg of a water manometer is connected to a small tube facing the
   stream (a pitot probe); the other is open to the still room air, which in a free jet is at the jet's
   static pressure. Then dp = (rho_w - rho_air) g dh and V = sqrt(2 dp/rho_air). Tilting the whole
   U-tube on a board at an angle theta multiplies the reading along the tube by 1/sin(theta).

Also: weighted straight-line fits, E_n, the resolution limit of the manometer, Reynolds numbers of the
jet, the probe and the fan, and a SYNTHETIC dataset generator for learners without a kit.
SYNTHETIC DATA ARE NOT MEASUREMENTS.

Units: SI in the functions (m, s, Pa, kg/m^3); the driver prints mm, s and m/s.

Run:  python main.py [out_dir]     (prediction, SYNTHETIC reduction, checks; writes flow_lab.png)
Test: pytest curriculum/phase-08/lesson-19-home-lab-flow-pressure-measurement
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

G = 9.81                       # m/s^2; use your local g if you know it (P00.14)
RHO_WATER, MU_WATER = 998.2, 1.002e-3      # water at 20 C
MU_AIR = 1.81e-5               # air at 20 C (Pa s); depends on T only weakly
R_AIR = 287.05287              # J/(kg K)
CD_ORIFICE_TYPICAL = 0.61      # typical discharge coefficient of a sharp-edged round orifice (textbook range 0.60-0.62)
U_CD_TYPICAL = 0.03            # standard uncertainty we allow for a hole in a thin plastic wall (not truly sharp-edged)
SYNTHETIC_SEED = 819


# ---------------------------------------------------------------- experiment 1: the draining bottle

def orifice_area(d: float) -> float:
    """Area of a round hole of diameter d: pi d^2/4."""
    raise NotImplementedError


def drain_rate_constant(cd: float, d_hole: float, D_bottle: float, g: float = G) -> float:
    """k = C_d (a/A) sqrt(g/2) in sqrt(h) = sqrt(h0) - k t (units m^0.5/s); a hole of diameter d_hole in
    a straight-walled bottle of inner diameter D_bottle."""
    raise NotImplementedError


def torricelli_level(t: np.ndarray, h0: float, cd: float, d_hole: float, D_bottle: float,
                     g: float = G) -> np.ndarray:
    """Water level above the hole centre, h(t) = (sqrt(h0) - k t)^2, and zero once the bottle has drained
    to the hole (t >= sqrt(h0)/k)."""
    raise NotImplementedError


def drain_time(h_start: float, h_end: float, cd: float | np.ndarray, d_hole: float | np.ndarray,
               D_bottle: float | np.ndarray, g: float = G) -> float | np.ndarray:
    """Time for the level to fall from h_start to h_end (both above the hole centre):
    T = (sqrt(h_start) - sqrt(h_end))/k. Vectorised over its arguments (for a Monte Carlo); a float
    when all arguments are scalars."""
    raise NotImplementedError


def predict_drain_time(h_start: float, h_end: float, cd: float, u_cd: float, d_hole: float, u_d: float,
                       D_bottle: float, u_D: float, u_h: float = 0.0) -> dict:
    """The prediction written down before measuring: T with its standard uncertainty by first-order
    propagation (central-difference sensitivities, relative step 1e-6) and the budget term by term
    (keys 'cd', 'd', 'D', 'h_start', 'h_end', each a standard uncertainty in seconds)."""
    raise NotImplementedError


def fit_line(x: np.ndarray, y: np.ndarray, u_y: float | np.ndarray) -> dict:
    """Weighted least squares for y = slope x + intercept with standard uncertainties u_y. Returns slope,
    intercept, their standard uncertainties (from (A^T W A)^-1), their covariance, chi2 and dof."""
    raise NotImplementedError


def fit_discharge_coefficient(t: np.ndarray, h: np.ndarray, u_h: float, u_t: float, d_hole: float, u_d: float,
                              D_bottle: float, u_D: float, g: float = G) -> dict:
    """Fit sqrt(h) = c - k t to readings (t_i, h_i) and return C_d = k (A/a) sqrt(2/g) with its uncertainty.
    The uncertainty of each sqrt(h_i) combines the level reading, u_h/(2 sqrt(h)), and the timing,
    k u_t (effective variance; k from an unweighted first pass). u(C_d) adds the slope's uncertainty
    and the geometry in quadrature: (u_Cd/Cd)^2 = (u_k/k)^2 + (2 u_d/d)^2 + (2 u_D/D)^2.
    Returns cd, u_cd, its budget (relative terms), the line fit, and the implied initial level."""
    raise NotImplementedError


def jet_reynolds(h: float | np.ndarray, d_hole: float, rho: float = RHO_WATER, mu: float = MU_WATER,
                 g: float = G) -> np.ndarray:
    """Reynolds number of the jet, rho V d/mu with the ideal Torricelli speed V = sqrt(2 g h)."""
    raise NotImplementedError


# ---------------------------------------------------------------- experiment 2: manometer and airspeed

def air_density(p: float, T: float) -> float:
    """Ideal-gas density of dry air, p/(R T) (P07.01): p in Pa, T in K."""
    raise NotImplementedError


def manometer_dp(dh: float | np.ndarray, rho_air: float, rho_liquid: float = RHO_WATER,
                 g: float = G) -> np.ndarray:
    """Pressure difference across a U-tube from the vertical height difference dh of the two menisci:
    dp = (rho_liquid - rho_air) g dh (the air columns above the menisci, which also differ by dh, are
    the small correction)."""
    raise NotImplementedError


def vertical_head(reading_along_tube: float | np.ndarray, angle_deg: float) -> np.ndarray:
    """Vertical height difference from the difference of the two meniscus positions measured ALONG a
    U-tube lying on a board tilted at angle_deg to the horizontal (90 = upright): dh = L sin(theta)."""
    raise NotImplementedError


def airspeed(dp: float | np.ndarray, rho_air: float) -> np.ndarray:
    """Incompressible pitot airspeed V = sqrt(2 dp/rho_air) (lesson P08.04); NaN for dp < 0."""
    raise NotImplementedError


def airspeed_with_uncertainty(L_reading: float, u_L: float, angle_deg: float, u_angle_deg: float,
                              rho_air: float, u_rho_air: float) -> tuple[float, float]:
    """V and u_V from an along-tube reading difference L (m) with standard uncertainty u_L, the tilt and
    its uncertainty, and the air density. First order: (u_V/V)^2 = (1/4)[(u_L/L)^2 + (u_theta/tan theta)^2
    + (u_rho/rho)^2] (u_theta in radians; the water density and g are known far better)."""
    raise NotImplementedError


def min_speed_for_uncertainty(u_dh: float, rel_target: float, rho_air: float,
                              rho_liquid: float = RHO_WATER, g: float = G) -> float:
    """The lowest airspeed at which a head uncertainty u_dh (vertical, m) still gives a relative airspeed
    uncertainty of rel_target: u_V/V = u_dh/(2 dh), so dh_min = u_dh/(2 rel_target) and
    V_min = sqrt(2 (rho_liquid - rho_air) g dh_min/rho_air)."""
    raise NotImplementedError


def reynolds(rho: float, V: float | np.ndarray, d: float, mu: float) -> np.ndarray:
    """Reynolds number rho V d/mu."""
    raise NotImplementedError


def fan_mean_velocity(Q_m3_per_h: float, blade_d: float, hub_d: float = 0.0) -> float:
    """Mean axial velocity through a fan's swept annulus from a rated volume flow (a MANUFACTURER'S
    FIGURE, usually measured in a test rig, not at your table): V = Q/(pi (D^2 - d_hub^2)/4)."""
    raise NotImplementedError


def en_number(x1: float, u1: float, x2: float, u2: float, k: float = 2.0) -> float:
    """Normalised error E_n = |x1 - x2| / sqrt(U1^2 + U2^2) with expanded uncertainties U = k u (P00.15)."""
    raise NotImplementedError


# ---------------------------------------------------------------- SYNTHETIC data (given)

BOTTLE_NOMINAL = {"d_hole": 4.0e-3, "u_d": 0.1e-3, "D_bottle": 90.0e-3, "u_D": 0.5e-3}
ROOM_NOMINAL = {"p": 101_000.0, "T": 294.15}


def make_synthetic_flow_data(seed: int = SYNTHETIC_SEED) -> dict:
    """SYNTHETIC measurements for both experiments (NOT a measurement; for learners without a kit).

    Draining bottle: a straight-walled bottle, nominal inner diameter 90.0 mm and hole 4.0 mm; the hidden
    truth is D = 90.4 mm, d = 4.05 mm, C_d = 0.62, level 205 mm above the hole centre when the tape comes
    off at t = 0. Three runs; for each, the time at which the level crosses the marks 200, 190, ..., 30 mm,
    read from a phone video with normal noise of 0.08 s on the time and 0.5 mm on where the level really
    was when it 'crossed' (meniscus and parallax).

    Manometer: room at 101.0 kPa and 21 C. True speeds 12.0 and 8.5 m/s (hair dryer, cool, high and low,
    50 mm from the nozzle) and 3.6 and 2.8 m/s (guarded fan, high and low, 300 mm in front). Each is read
    on the upright U-tube (90 deg) and tilted at 10 deg, twice. The along-tube position of each meniscus is
    read to a normal noise of 0.3 mm; the stream's unsteadiness adds a normal 2 % to each dp. An
    'anemometer_mps' column imitates an enhanced-tier vane anemometer (3 % + 0.1 m/s noise).

    Returns {'drain': array (run, mark_mm, t_s), 'manometer': list of row dicts}."""
    rng = np.random.default_rng(seed)
    D_t, d_t, cd_t, h0_t = 90.4e-3, 4.05e-3, 0.62, 0.205
    k = drain_rate_constant(cd_t, d_t, D_t)
    marks = np.arange(200, 29, -10) * 1e-3
    drain = []
    for run in (1, 2, 3):
        h_true = marks + rng.normal(0.0, 0.5e-3, marks.size)
        t = (np.sqrt(h0_t) - np.sqrt(h_true)) / k + rng.normal(0.0, 0.08, marks.size)
        drain += [(run, round(m * 1e3), round(float(tt), 2)) for m, tt in zip(marks, t)]
    rho_a = air_density(ROOM_NOMINAL["p"], ROOM_NOMINAL["T"])
    rows = []
    cases = [("hair dryer cool", "high", 50, 12.0), ("hair dryer cool", "low", 50, 8.5),
             ("fan", "high", 300, 3.6), ("fan", "low", 300, 2.8)]
    for source, setting, dist, V in cases:
        for angle in (90.0, 10.0):
            for rep in (1, 2):
                dp = 0.5 * rho_a * V ** 2 * (1.0 + rng.normal(0.0, 0.02))
                L = float(dp / ((RHO_WATER - rho_a) * G) / np.sin(np.radians(angle)))
                zl, zr = 100.0 + rng.normal(0, 0.3), 100.0 + rng.normal(0, 0.3)
                left = 100.0 - 0.5 * L * 1e3 + rng.normal(0, 0.3)
                right = 100.0 + 0.5 * L * 1e3 + rng.normal(0, 0.3)
                anemo = V * (1 + rng.normal(0, 0.03)) + rng.normal(0, 0.1)
                rows.append({"source": source, "setting": setting, "distance_mm": dist, "angle_deg": angle,
                             "repeat": rep, "zero_left_mm": round(zl, 1), "zero_right_mm": round(zr, 1),
                             "left_mm": round(left, 1), "right_mm": round(right, 1),
                             "anemometer_mps": round(anemo, 2)})
    return {"drain": np.array(drain, dtype=float), "manometer": rows}


def write_synthetic_csv(out_dir: str | Path, seed: int = SYNTHETIC_SEED) -> list[Path]:
    """Write the SYNTHETIC dataset as two CSV files whose names and first lines say SYNTHETIC."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    data = make_synthetic_flow_data(seed)
    head = f"# SYNTHETIC data from make_synthetic_flow_data(seed={seed}) in P08.19 solution.py: NOT a measurement\n"
    p1 = out_dir / f"SYNTHETIC_seed{seed}_drain.csv"
    with open(p1, "w", newline="\n") as fh:
        fh.write(head)
        fh.write("# bottle inner diameter 90.0 +/- 0.5 mm (nominal), hole 4.0 +/- 0.1 mm (nominal); "
                 "level read above the hole centre\n")
        fh.write("run,mark_mm,t_s\n")
        for run, mark, t in data["drain"]:
            fh.write(f"{int(run)},{int(mark)},{t:.2f}\n")
    p2 = out_dir / f"SYNTHETIC_seed{seed}_manometer.csv"
    keys = ["source", "setting", "distance_mm", "angle_deg", "repeat", "zero_left_mm", "zero_right_mm",
            "left_mm", "right_mm", "anemometer_mps"]
    with open(p2, "w", newline="\n") as fh:
        fh.write(head)
        fh.write("# room 101.0 kPa, 21.0 C; readings are meniscus positions ALONG the tube (mm); "
                 "pitot on the left leg; angle 90 = upright\n")
        fh.write(",".join(keys) + "\n")
        for r in data["manometer"]:
            fh.write(",".join(f"{r[k]:g}" if isinstance(r[k], float) else str(r[k]) for k in keys) + "\n")
    return [p1, p2]


def read_csv(path: str | Path) -> dict[str, np.ndarray]:
    """Read a CSV with '#' comment lines into {column: array}; numeric columns become floats."""
    lines = [ln for ln in Path(path).read_text().splitlines() if ln and not ln.startswith("#")]
    keys = lines[0].split(",")
    cols = list(zip(*[ln.split(",") for ln in lines[1:]]))
    out = {}
    for k, c in zip(keys, cols):
        try:
            out[k] = np.array(c, dtype=float)
        except ValueError:
            out[k] = np.array(c)
    return out


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Prediction, SYNTHETIC reduction and checks; writes flow_lab.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    nb = BOTTLE_NOMINAL
    print("EXPERIMENT 1: draining bottle (prediction, before any measurement)")
    pr = predict_drain_time(0.200, 0.030, CD_ORIFICE_TYPICAL, U_CD_TYPICAL, nb["d_hole"], nb["u_d"],
                            nb["D_bottle"], nb["u_D"], u_h=0.5e-3)
    print(f"  bottle 90.0 +/- 0.5 mm, hole 4.0 +/- 0.1 mm, C_d = {CD_ORIFICE_TYPICAL} +/- {U_CD_TYPICAL} (typical)")
    print(f"  predicted time from 200 mm to 30 mm: {pr['T']:.1f} +/- {pr['u_T']:.1f} s "
          f"({pr['u_T'] / pr['T'] * 100:.1f} %); budget (s): "
          + ", ".join(f"{k} {v:.2f}" for k, v in pr["budget"].items()))
    print(f"  to empty from 200 mm: {drain_time(0.2, 0.0, 0.61, nb['d_hole'], nb['D_bottle']):.1f} s; "
          f"water used up to 200 mm above a hole 30 mm from the base: "
          f"{np.pi * nb['D_bottle'] ** 2 / 4 * 0.230 * 1e3:.2f} L")
    re = jet_reynolds(np.array([0.2, 0.03]), nb["d_hole"])
    print(f"  jet Reynolds number: {re[0]:.0f} at 200 mm, {re[1]:.0f} at 30 mm; area ratio (a/A)^2 = "
          f"{(nb['d_hole'] / nb['D_bottle']) ** 4:.1e}")

    data = make_synthetic_flow_data()
    dr = data["drain"]
    print(f"\nSYNTHETIC dataset seed {SYNTHETIC_SEED} (not a measurement): C_d by run")
    res = []
    for run in (1, 2, 3):
        m = dr[:, 0] == run
        r = fit_discharge_coefficient(dr[m, 2], dr[m, 1] * 1e-3, 0.5e-3, 0.08, nb["d_hole"], nb["u_d"],
                                      nb["D_bottle"], nb["u_D"])
        res.append(r)
        print(f"  run {run}: C_d = {r['cd']:.4f} +/- {r['u_cd']:.4f}  (slope term {r['rel_budget']['slope'] * 100:.2f} %, "
              f"hole {r['rel_budget']['d'] * 100:.1f} %, bottle {r['rel_budget']['D'] * 100:.1f} %); "
              f"chi2/dof = {r['fit']['chi2']:.1f}/{r['fit']['dof']}; h0 = {r['h0'] * 1e3:.1f} mm")
    cds = np.array([r["cd"] for r in res])
    sd = cds.std(ddof=1)
    print(f"  mean of 3 runs {cds.mean():.4f}, run-to-run sd {sd:.4f}; E_n against the typical "
          f"{CD_ORIFICE_TYPICAL} +/- {U_CD_TYPICAL}: {en_number(cds.mean(), res[0]['u_cd'], CD_ORIFICE_TYPICAL, U_CD_TYPICAL):.2f}")
    m1 = dr[:, 0] == 1
    T_meas = dr[m1, 2][-1] - dr[m1, 2][0]
    print(f"  run 1 measured 200 -> 30 mm: {T_meas:.1f} s against the prediction {pr['T']:.1f} +/- {pr['u_T']:.1f} s: "
          f"E_n = {en_number(T_meas, 0.08 * np.sqrt(2), pr['T'], pr['u_T']):.2f}")

    print("\nEXPERIMENT 2: U-tube manometer (prediction)")
    rho_a = air_density(ROOM_NOMINAL["p"], ROOM_NOMINAL["T"])
    print(f"  room 101.0 kPa, 21 C: rho_air = {rho_a:.4f} kg/m^3")
    for V in (3.0, 5.0, 10.0, 15.0):
        dh = 0.5 * rho_a * V ** 2 / ((RHO_WATER - rho_a) * G)
        print(f"  V = {V:4.1f} m/s: dp = {0.5 * rho_a * V ** 2:6.1f} Pa, dh = {dh * 1e3:5.2f} mm upright, "
              f"{dh / np.sin(np.radians(10)) * 1e3:5.1f} mm along a tube tilted at 10 deg")
    u_L = 0.3e-3 * 2.0                  # four readings (two menisci, flow on and zero), each to 0.3 mm
    for ang in (90.0, 30.0, 10.0):
        v10 = min_speed_for_uncertainty(u_L * np.sin(np.radians(ang)), 0.10, rho_a)
        print(f"  lowest speed read to 10 % (four readings of 0.3 mm), tilt {ang:.0f} deg: {v10:.2f} m/s")
    for name, V, d in (("probe (6 mm straw) at 3 m/s", 3.0, 0.006), ("probe at 12 m/s", 12.0, 0.006),
                       ("fan, 280 mm blades, at 3.6 m/s", 3.6, 0.28), ("dryer nozzle 40 mm at 12 m/s", 12.0, 0.04)):
        print(f"  Reynolds number, {name}: {float(reynolds(rho_a, V, d, MU_AIR)):.0f}")
    print(f"  fan rated 1100 m^3/h (manufacturer's figure, SYNTHETIC example), 280 mm blades, 80 mm hub: "
          f"mean {fan_mean_velocity(1100, 0.28, 0.08):.2f} m/s through the annulus")

    print(f"\nSYNTHETIC manometer readings, seed {SYNTHETIC_SEED} (not a measurement):")
    print("  source            setting angle  L_along (mm)  dh (mm)   V (m/s)    anemometer (m/s)")
    rows = data["manometer"]
    table = []
    for r in rows:
        L = ((r["right_mm"] - r["left_mm"]) - (r["zero_right_mm"] - r["zero_left_mm"])) * 1e-3
        V, uV = airspeed_with_uncertainty(L, 0.3e-3 * 2, r["angle_deg"], 0.5, rho_a, 0.01 * rho_a)
        table.append((r["source"], r["setting"], r["angle_deg"], L, V, uV, r["anemometer_mps"]))
        print(f"  {r['source']:17s} {r['setting']:5s} {r['angle_deg']:5.0f}   {L * 1e3:7.1f}     "
              f"{float(vertical_head(L, r['angle_deg'])) * 1e3:6.2f}   {V:5.2f} +/- {uV:4.2f}   {r['anemometer_mps']:.2f}")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for run, mk in zip((1, 2, 3), ("o", "s", "^")):
        m = dr[:, 0] == run
        ax[0].plot(dr[m, 2], np.sqrt(dr[m, 1] * 1e-3), mk, mfc="none", label=f"SYNTHETIC run {run}")
    tt = np.linspace(0, dr[:, 2].max(), 50)
    ax[0].plot(tt, np.sqrt(0.2) + 0 * tt - drain_rate_constant(0.61, nb["d_hole"], nb["D_bottle"]) * (tt - dr[m1, 2][0]),
               "k--", lw=0.8, label="prediction, C_d = 0.61 (typical)")
    ax[0].set_xlabel("t (s)")
    ax[0].set_ylabel("sqrt(h) (m^0.5)")
    ax[0].set_title("SYNTHETIC: draining bottle")
    ax[0].legend(fontsize=8)
    for ang, mk in ((90.0, "o"), (10.0, "s")):
        sel = [x for x in table if x[2] == ang]
        ax[1].errorbar([x[6] for x in sel], [x[4] for x in sel], yerr=[2 * x[5] for x in sel], fmt=mk,
                       label=f"SYNTHETIC, tilt {ang:.0f} deg (k = 2 bars)")
    ax[1].plot([0, 14], [0, 14], "k:", lw=0.8)
    ax[1].set_xlabel("anemometer (m/s)")
    ax[1].set_ylabel("manometer airspeed (m/s)")
    ax[1].set_title("SYNTHETIC: pitot + U-tube")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "flow_lab.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'flow_lab.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
