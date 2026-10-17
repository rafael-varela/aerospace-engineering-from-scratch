"""P04.15 Electronics Bench Lab: starter code.

Complete every function marked NotImplementedError. When finished, running this file produces the lesson's output and figure.

Analysis code for the bench lab in projects/p04-electronics-bench/:

* design: voltage divider, LED series resistor, nearest standard (E12) values, resistor power;
* instruments: the loading of a circuit by a voltmeter's input resistance and by an ammeter's shunt;
* the ADC: counts to volts (n-bit and MicroPython's 16-bit ``read_u16`` scale), the LSB, the
  quantization noise LSB/sqrt(12), an ideal quantizer, and a two-parameter calibration against a multimeter;
* the thermistor: beta model and Steinhart-Hart, both fitted by linear least squares, the ratiometric
  divider read by the ADC (the reference voltage cancels), and temperature from counts;
* uncertainty: a multimeter's accuracy specification as a Type B term and first-order propagation with a
  finite-difference Jacobian (as in P00.14).

A SYNTHETIC dataset (generated here from a known thermistor, labelled synthetic in every file and column
note) lets a learner without the kit run the whole analysis. It is a simulation, not a measurement.

Run:  python main.py          (prints the design and analysis tables, writes bench_lab.png)
Test: pytest curriculum/phase-04/lesson-15-electronics-bench-lab
"""
from __future__ import annotations

import csv
from pathlib import Path
from typing import Callable

import numpy as np

KELVIN = 273.15
E12 = np.array([1.0, 1.2, 1.5, 1.8, 2.2, 2.7, 3.3, 3.9, 4.7, 5.6, 6.8, 8.2])
U16_FULL_SCALE = 65535          # MicroPython ADC.read_u16() returns 0..65535

# the "true" thermistor used to generate the SYNTHETIC dataset: Steinhart-Hart coefficients of a
# generic 10 kOhm NTC (about 10.0 kOhm at 25 C, beta about 3950 K)
SYNTH_SH = (1.125e-3, 2.347e-4, 8.6e-8)
SYNTH_R_FIXED = 10_000.0        # ohm, the bottom resistor of the synthetic divider
SYNTH_VREF = 3.3                # V


# ---------------------------------------------------------------- design

def nearest_e_series(value: float, series: np.ndarray = E12, direction: str = "nearest") -> float:
    """Standard resistor value from a decade series closest to ``value`` ('nearest'), or the smallest
    one not below it ('up'), or the largest not above it ('down')."""
    raise NotImplementedError


def divider_output(vin: float, r_top: float, r_bottom: float, r_load: float = np.inf) -> float:
    """Output of a two-resistor divider, optionally loaded by r_load in parallel with r_bottom:
    vin * Rb' / (r_top + Rb'),  Rb' = r_bottom || r_load."""
    raise NotImplementedError


def divider_design(vin: float, vout: float, r_total: float) -> tuple[float, float]:
    """Exact (r_top, r_bottom) giving vout from vin with r_top + r_bottom = r_total."""
    raise NotImplementedError


def led_resistor(v_supply: float, v_forward: float, i_target: float) -> float:
    """Series resistor (exact, before rounding) for an LED: (V_supply - V_f) / I."""
    raise NotImplementedError


def led_current(v_supply: float, v_forward: float, r: float, r_extra: float = 0.0) -> float:
    """LED current with a fixed-forward-voltage LED model: (V_supply - V_f) / (r + r_extra)."""
    raise NotImplementedError


def resistor_power(current: float, r: float) -> float:
    """Power I^2 R dissipated in a resistor [W]."""
    raise NotImplementedError


# ---------------------------------------------------------------- instruments

def loading_error(r_top: float, r_bottom: float, r_meter: float) -> float:
    """Fractional error (reading - true)/true of a voltmeter of input resistance r_meter across r_bottom
    of a divider. Equals -R_th / (R_th + r_meter) with R_th = r_top || r_bottom (always negative)."""
    raise NotImplementedError


def thevenin_resistance(r_top: float, r_bottom: float) -> float:
    """Thevenin (source) resistance seen at the divider output: r_top || r_bottom."""
    raise NotImplementedError


# ---------------------------------------------------------------- the ADC

def lsb(v_ref: float, n_bits: int) -> float:
    """Voltage of one least-significant bit: v_ref / 2^n."""
    raise NotImplementedError


def quantization_noise_rms(v_ref: float, n_bits: int) -> float:
    """RMS quantization error of an ideal rounding quantizer: LSB / sqrt(12)."""
    raise NotImplementedError


def volts_to_counts(v: np.ndarray | float, v_ref: float, n_bits: int) -> np.ndarray:
    """Ideal n-bit quantizer: round(v / LSB), clipped to 0 .. 2^n - 1."""
    raise NotImplementedError


def counts_to_volts(counts: np.ndarray | float, v_ref: float, n_bits: int) -> np.ndarray | float:
    """Inverse of the ideal quantizer: counts x LSB."""
    raise NotImplementedError


def u16_to_volts(raw: np.ndarray | float, v_ref: float = 3.3) -> np.ndarray | float:
    """MicroPython ``ADC.read_u16()`` value (0..65535) to volts: raw / 65535 x v_ref.

    Most ports scale a 10- or 12-bit conversion up to 16 bits, so the true resolution is that of the
    hardware, not 16 bits; the low bits are either zero or repeat the high bits."""
    raise NotImplementedError


def adc_calibration_fit(raw: np.ndarray, volts: np.ndarray) -> tuple[float, float, float]:
    """Straight-line calibration volts = gain x raw + offset against a multimeter, by least squares.
    Returns (gain [V/count], offset [V], residual standard deviation [V])."""
    raise NotImplementedError


# ---------------------------------------------------------------- the thermistor

def thermistor_beta_r(T_K: np.ndarray | float, r0: float, beta: float, T0: float = 298.15) -> np.ndarray | float:
    """Beta model: R(T) = R0 exp(beta (1/T - 1/T0)), temperatures in kelvin."""
    raise NotImplementedError


def thermistor_beta_t(R: np.ndarray | float, r0: float, beta: float, T0: float = 298.15) -> np.ndarray | float:
    """Inverse beta model: T = 1 / (1/T0 + ln(R/R0)/beta) [K]."""
    raise NotImplementedError


def fit_beta(R: np.ndarray, T_K: np.ndarray, T0: float = 298.15) -> tuple[float, float]:
    """Least-squares fit of ln R = ln R0 + beta (1/T - 1/T0). Returns (R0, beta)."""
    raise NotImplementedError


def steinhart_hart_t(R: np.ndarray | float, a: float, b: float, c: float) -> np.ndarray | float:
    """Steinhart-Hart: 1/T = a + b ln R + c (ln R)^3 [T in K, R in ohm]."""
    raise NotImplementedError


def steinhart_hart_r(T_K: np.ndarray | float, a: float, b: float, c: float) -> np.ndarray | float:
    """Inverse Steinhart-Hart by Cardano's formula for the depressed cubic c L^3 + b L + (a - 1/T) = 0."""
    raise NotImplementedError


def fit_steinhart_hart(R: np.ndarray, T_K: np.ndarray) -> tuple[float, float, float]:
    """Least-squares fit of 1/T = a + b ln R + c (ln R)^3 (linear in a, b, c). Needs >= 3 points."""
    raise NotImplementedError


def thermistor_r_from_ratio(ratio: np.ndarray | float, r_fixed: float) -> np.ndarray | float:
    """Thermistor on top, fixed resistor to ground, ADC at the junction: ratio = V_out/V_ref
    = R_fixed/(R_th + R_fixed), so R_th = R_fixed (1/ratio - 1). V_ref cancels when the divider is
    powered from the ADC's own reference (ratiometric measurement)."""
    raise NotImplementedError


def temperature_from_u16(raw: np.ndarray | float, r_fixed: float, sh: tuple[float, float, float]) -> np.ndarray | float:
    """Temperature [deg C] from a read_u16 value of the ratiometric thermistor divider."""
    raise NotImplementedError


def best_fixed_resistor(T_lo_K: float, T_hi_K: float, r0: float, beta: float, T0: float = 298.15) -> float:
    """Fixed resistor that maximises the voltage swing of the divider between two temperatures:
    sqrt(R(T_lo) R(T_hi)) (the geometric mean)."""
    raise NotImplementedError


# ---------------------------------------------------------------- uncertainty

def dmm_standard_uncertainty(reading: float, pct: float, digits: int, resolution: float) -> float:
    """Standard uncertainty of a multimeter reading from an accuracy spec +-(pct % of reading + digits x
    resolution), treated as the half-width of a rectangular distribution (Type B, divide by sqrt 3)."""
    raise NotImplementedError


def propagate(f: Callable[[np.ndarray], float], x: np.ndarray, u: np.ndarray,
              rel_step: float = 1e-6) -> tuple[float, float]:
    """First-order propagation of independent standard uncertainties u through f at x, with a central
    finite-difference gradient (P00.14). Returns (f(x), u_f)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the SYNTHETIC dataset

def make_synthetic_bench_data(seed: int = 2026) -> dict[str, np.ndarray]:
    """SYNTHETIC (simulated) bench data with a known truth, for learners without the kit and for tests.

    Thermistor sweep: reference temperatures 5..45 C read with a +-0.5 C thermometer (Gaussian sigma
    0.25 C), a 12-bit ADC scaled to 16 bits (as many MicroPython ports do), 16 samples averaged, an ADC
    noise of 1.5 LSB(12) rms, and a 10 kOhm fixed resistor that is really 9.95 kOhm (0.5 % low).
    ADC calibration: five potentiometer settings read by the ADC (with 0.8 % gain error and 6 mV offset)
    and by a multimeter (+-(0.5 % + 2 digits) at 1 mV resolution).
    Nothing here is a measurement.
    """
    rng = np.random.default_rng(seed)
    T_ref = np.array([5.0, 15.0, 25.0, 35.0, 45.0])
    T_true = T_ref + rng.normal(0.0, 0.25, T_ref.size)
    R_true = steinhart_hart_r(T_true + KELVIN, *SYNTH_SH)
    r_fixed_true = 9_950.0
    ratio = r_fixed_true / (R_true + r_fixed_true)
    n_samp = 16
    samples = ratio[:, None] * 4095 + rng.normal(0.0, 1.5, (T_ref.size, n_samp))
    c12 = np.clip(np.round(samples), 0, 4095)
    u16 = c12 * 16 + c12 // 256                       # 12-bit to 16-bit scaling with bit replication
    pot_true = np.array([0.30, 0.95, 1.65, 2.40, 3.05])
    adc_v = pot_true * 1.008 + 0.006
    adc_raw = np.round(adc_v / SYNTH_VREF * 4095) * 16
    dmm = np.round(pot_true + rng.normal(0.0, 0.002, pot_true.size), 3)
    return {"ref_temp_C": T_ref, "u16_mean": u16.mean(axis=1), "u16_std": u16.std(axis=1, ddof=1),
            "n_samples": np.full(T_ref.size, n_samp), "pot_dmm_V": dmm, "pot_u16": adc_raw}


def write_synthetic_csv(out_dir: str | Path, seed: int = 2026) -> tuple[Path, Path]:
    """Write the SYNTHETIC dataset as two CSV files whose first line says it is synthetic.
    Returns the two paths."""
    out_dir = Path(out_dir)
    d = make_synthetic_bench_data(seed)
    p1 = out_dir / "synthetic_thermistor_sweep.csv"
    p2 = out_dir / "synthetic_adc_calibration.csv"
    head = f"# SYNTHETIC DATA generated by make_synthetic_bench_data(seed={seed}); a simulation, not a measurement\n"
    with p1.open("w", newline="") as fh:
        fh.write(head)
        w = csv.writer(fh)
        w.writerow(["ref_temp_C", "u16_mean", "u16_std", "n_samples", "r_fixed_nominal_ohm", "note"])
        for i in range(d["ref_temp_C"].size):
            w.writerow([f"{d['ref_temp_C'][i]:.1f}", f"{d['u16_mean'][i]:.2f}", f"{d['u16_std'][i]:.2f}",
                        int(d["n_samples"][i]), f"{SYNTH_R_FIXED:.0f}", "synthetic"])
    with p2.open("w", newline="") as fh:
        fh.write(head)
        w = csv.writer(fh)
        w.writerow(["setting", "dmm_V", "u16_mean", "note"])
        for i in range(d["pot_dmm_V"].size):
            w.writerow([i + 1, f"{d['pot_dmm_V'][i]:.3f}", f"{d['pot_u16'][i]:.0f}", "synthetic"])
    return p1, p2


def read_bench_csv(path: str | Path) -> dict[str, np.ndarray]:
    """Read a bench CSV (comment lines starting with '#' skipped) into numeric columns; non-numeric
    columns are returned as string arrays."""
    with Path(path).open() as fh:
        rows = [r for r in csv.reader(line for line in fh if not line.startswith("#"))]
    header, body = rows[0], rows[1:]
    out: dict[str, np.ndarray] = {}
    for j, name in enumerate(header):
        col = [r[j] for r in body]
        try:
            out[name] = np.array([float(v) for v in col])
        except ValueError:
            out[name] = np.array(col)
    return out


# ---------------------------------------------------------------- main

def analyse_thermistor(data: dict[str, np.ndarray], r_fixed: float = SYNTH_R_FIXED) -> dict[str, object]:
    """Reduce a thermistor sweep: resistances, beta and Steinhart-Hart fits, residuals in deg C."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the design and analysis tables for the SYNTHETIC dataset and write bench_lab.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print("DESIGN")
    rt, rb = divider_design(5.0, 3.0, 20e3)
    print(f"  divider 5 V -> 3.0 V, 20 kOhm total: exact {rt:.0f} / {rb:.0f} ohm; E12 "
          f"{nearest_e_series(rt):.0f} / {nearest_e_series(rb):.0f} -> {divider_output(5.0, nearest_e_series(rt), nearest_e_series(rb)):.3f} V")
    for vs in (5.0, 3.3):
        r = led_resistor(vs, 2.0, 10e-3)
        ru = nearest_e_series(r, direction="up")
        i = led_current(vs, 2.0, ru)
        print(f"  red LED (Vf 2.0 V, 10 mA) from {vs} V: exact {r:.0f} ohm, E12 up {ru:.0f} ohm -> "
              f"{i * 1e3:.2f} mA, {resistor_power(i, ru) * 1e3:.1f} mW in the resistor")
    print("METER LOADING (10 MOhm voltmeter across the bottom resistor of an equal divider)")
    for r in (1e3, 10e3, 100e3, 1e6, 10e6):
        print(f"  {r:9.0f} ohm + {r:9.0f} ohm: R_th {thevenin_resistance(r, r):9.0f} ohm, error {loading_error(r, r, 10e6) * 100:+.4f} %")
    i0 = led_current(5.0, 2.0, 330.0)
    i1 = led_current(5.0, 2.0, 330.0, r_extra=10.0)
    print(f"AMMETER BURDEN: 330 ohm LED circuit at 5 V: {i0 * 1e3:.3f} mA without, {i1 * 1e3:.3f} mA with a 10 ohm shunt "
          f"({(i1 / i0 - 1) * 100:+.2f} %)")
    print("ADC")
    for n in (10, 12, 16):
        print(f"  {n:2d} bit at 3.3 V: LSB {lsb(3.3, n) * 1e3:7.4f} mV, quantization noise {quantization_noise_rms(3.3, n) * 1e3:7.4f} mV rms")
    print("SYNTHETIC DATASET (simulation, not a measurement)")
    d = make_synthetic_bench_data()
    g, o, s = adc_calibration_fit(d["pot_u16"], d["pot_dmm_V"])
    print(f"  ADC calibration vs multimeter: gain {g * 1e6:.4f} uV/count (ideal {3.3 / 65535 * 1e6:.4f}), "
          f"offset {o * 1e3:+.2f} mV, residual sd {s * 1e3:.2f} mV")
    a = analyse_thermistor(d)
    r0, beta = a["beta"]
    print(f"  beta fit: R25 = {r0:.0f} ohm, beta = {beta:.0f} K; residuals {np.round(a['res_beta_C'], 3)} C")
    print(f"  Steinhart-Hart fit: a = {a['sh'][0]:.4e}, b = {a['sh'][1]:.4e}, c = {a['sh'][2]:.3e}; "
          f"residuals {np.round(a['res_sh_C'], 3)} C")
    raw25 = d["u16_mean"][2]

    def t_of(p: np.ndarray) -> float:
        return float(temperature_from_u16(p[0], p[1], a["sh"]))
    u_raw = np.hypot(d["u16_std"][2] / np.sqrt(16), 16 / np.sqrt(12))
    T25, uT = propagate(t_of, np.array([raw25, 10_000.0]), np.array([u_raw, 10_000.0 * 0.01 / np.sqrt(3)]))
    print(f"  25 C point: raw {raw25:.1f} +- {u_raw:.1f} counts, R_fixed 10 kOhm +-1 % (rectangular): "
          f"T = {T25:.2f} C, u(T) = {uT:.3f} C")
    print(f"  best fixed resistor for 0-50 C (beta model): {best_fixed_resistor(273.15, 323.15, r0, beta):.0f} ohm")

    fig, ax = plt.subplots(1, 3, figsize=(14, 4))
    Tg = np.linspace(0, 50, 200) + KELVIN
    ax[0].semilogy(d["ref_temp_C"], a["R"], "o", label="SYNTHETIC data")
    ax[0].semilogy(Tg - KELVIN, steinhart_hart_r(Tg, *a["sh"]), label="Steinhart-Hart fit")
    ax[0].semilogy(Tg - KELVIN, thermistor_beta_r(Tg, r0, beta), "--", label="beta fit")
    ax[0].set_xlabel("temperature (C)"); ax[0].set_ylabel("resistance (ohm)"); ax[0].legend()
    ax[0].set_title("thermistor (synthetic data)")
    ax[1].plot(d["ref_temp_C"], a["res_beta_C"], "s-", label="beta")
    ax[1].plot(d["ref_temp_C"], a["res_sh_C"], "o-", label="Steinhart-Hart")
    ax[1].axhline(0, color="gray", lw=0.5)
    ax[1].set_xlabel("reference temperature (C)"); ax[1].set_ylabel("fit residual (C)"); ax[1].legend()
    ax[1].set_title("residuals (synthetic data)")
    rr = np.logspace(2, 7.3, 100)
    ax[2].semilogx(rr, -100 * np.array([loading_error(r, r, 10e6) for r in rr]))
    ax[2].set_xlabel("each divider resistor (ohm)"); ax[2].set_ylabel("voltmeter reads low by (%)")
    ax[2].set_title("loading by a 10 MOhm voltmeter")
    fig.tight_layout()
    fig.savefig(out_dir / "bench_lab.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
