"""P04.10 AC Circuits & Electronics Basics: starter code.

The signal chain of every digital sensor on an aircraft or a spacecraft: a sensor voltage, an analog
filter, a buffer or an amplifier, an analog-to-digital converter (ADC), and a microcontroller that
reads integers. This module builds the analog half with phasors and the digital half with sampling
and quantisation:

* impedances of R, L and C as complex numbers; series and parallel combination;
* the series RLC circuit: current phasor, resonance, quality factor, bandwidth;
* RMS values and AC power on a 115 V, 400 Hz aircraft bus (power factor);
* the RC low-pass filter: transfer function, corner frequency, phase, group delay, step response;
  cascaded buffered stages and the loaded (unbuffered) two-stage ladder;
* the ideal op-amp: buffer, non-inverting and inverting amplifiers, and the closed-loop gain with a
  finite open-loop gain A (the limit A -> infinity recovers the ideal result);
* the diode equation and a transistor used as a switch, at the conceptual level;
* logic gates, a half adder, and what an ADC returns: aliasing, the LSB, quantisation;
* the design of an anti-aliasing filter for a sensor sampled at 100 Hz by a 12-bit ADC, checked in
  the time domain with ``scipy.signal.lsim`` on SYNTHETIC test signals.

Units: SI throughout (ohm, farad, henry, hertz, volt, ampere). Phasors use the e^{j omega t}
convention, so Z_L = j omega L and Z_C = 1 / (j omega C).

Run:  python main.py [out_dir]     (prints the verification tables, writes rc_filter.png)
Test: pytest curriculum/phase-04/lesson-10-ac-circuits-electronics-basics
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import signal

K_BOLTZMANN = 1.380649e-23     # J/K
Q_E = 1.602176634e-19          # C


# ---------------------------------------------------------------- impedance and phasors

def impedance_R(R: float, f: np.ndarray | float) -> np.ndarray | complex:
    """Impedance of a resistor: R, at every frequency (returned as complex, broadcast to f)."""
    raise NotImplementedError


def impedance_L(L: float, f: np.ndarray | float) -> np.ndarray | complex:
    """Impedance of an inductor: j omega L, omega = 2 pi f."""
    raise NotImplementedError


def impedance_C(C: float, f: np.ndarray | float) -> np.ndarray | complex:
    """Impedance of a capacitor: 1 / (j omega C). Infinite at f = 0 (an open circuit at DC)."""
    raise NotImplementedError


def series(*Z) -> np.ndarray | complex:
    """Impedances in series add."""
    raise NotImplementedError


def parallel(*Z) -> np.ndarray | complex:
    """Impedances in parallel: admittances add, 1/Z = sum 1/Z_k."""
    raise NotImplementedError


def rlc_series_current(V: complex, R: float, L: float, C: float, f: np.ndarray | float) -> np.ndarray | complex:
    """Current phasor I = V / (R + j omega L + 1/(j omega C)) of a series RLC driven by the phasor V."""
    raise NotImplementedError


def resonant_frequency(L: float, C: float) -> float:
    """f0 = 1 / (2 pi sqrt(L C)): where X_L = X_C and the series impedance is purely R."""
    raise NotImplementedError


def quality_factor(R: float, L: float, C: float) -> float:
    """Series-RLC quality factor Q = (1/R) sqrt(L/C) = omega0 L / R; the -3 dB bandwidth is f0 / Q."""
    raise NotImplementedError


# ---------------------------------------------------------------- RMS and AC power

def rms(x: np.ndarray) -> float:
    """Root-mean-square of a sampled waveform (over a whole number of periods for a periodic signal)."""
    raise NotImplementedError


def ac_power(V_rms: float, Z: complex) -> tuple[float, float, float]:
    """(P, Q, pf) for a load Z on a sinusoidal source of RMS voltage V_rms:
    S = V I* = |V|^2 / Z*, P = Re S (watts), Q = Im S (vars), power factor pf = cos(arg Z) = P / |S|."""
    raise NotImplementedError


# ---------------------------------------------------------------- the RC low-pass filter

def rc_corner_frequency(R: float, C: float) -> float:
    """Corner (-3 dB) frequency of a first-order RC filter: f_c = 1 / (2 pi R C)."""
    raise NotImplementedError


def rc_lowpass(f: np.ndarray | float, R: float, C: float, stages: int = 1) -> np.ndarray | complex:
    """Transfer function V_out / V_in of ``stages`` identical RC low-pass sections separated by ideal
    buffers (so no stage loads the previous one): H = (1 / (1 + j f / f_c))^stages."""
    raise NotImplementedError


def rc_highpass(f: np.ndarray | float, R: float, C: float) -> np.ndarray | complex:
    """First-order RC high-pass (series C, shunt R): H = (j f / f_c) / (1 + j f / f_c)."""
    raise NotImplementedError


def rc2_unbuffered(f: np.ndarray | float, R: float, C: float) -> np.ndarray | complex:
    """Two identical RC sections connected directly (the second loads the first). Nodal analysis with
    s = j omega and x = s R C gives H = 1 / (1 + 3 x + x^2): not the buffered 1 / (1 + x)^2 = 1/(1 + 2x + x^2)."""
    raise NotImplementedError


def rc_loaded(f: np.ndarray | float, R: float, C: float, R_load: float) -> np.ndarray | complex:
    """One RC low-pass driving a resistive load R_load (for example an amplifier or an ADC input):
    H = Z_p / (R + Z_p) with Z_p = C || R_load. DC gain R_load / (R + R_load) and a higher corner."""
    raise NotImplementedError


def gain_db(H: np.ndarray | complex) -> np.ndarray | float:
    """20 log10 |H|."""
    raise NotImplementedError


def group_delay_dc(R: float, C: float, stages: int = 1) -> float:
    """Low-frequency group delay -d(phase)/d(omega) of ``stages`` buffered RC sections: stages * R C."""
    raise NotImplementedError


def rc_step_response(t: np.ndarray, R: float, C: float, V: float = 1.0) -> np.ndarray:
    """Capacitor voltage after a step V at t = 0 from rest: V (1 - exp(-t / RC)) (the RC time constant of P04.08)."""
    raise NotImplementedError


def rc_lti(R: float, C: float, stages: int = 1) -> signal.lti:
    """The buffered RC cascade as a ``scipy.signal.lti`` object: H(s) = 1 / (RC s + 1)^stages."""
    raise NotImplementedError


def required_stages_and_corner(f_signal: float, f_reject: float, attenuation: float,
                               max_droop_db: float, max_stages: int = 4) -> tuple[int, float]:
    """Smallest number n of identical buffered RC stages, and the highest corner f_c, such that
    |H(f_reject)| <= 1 / attenuation and the passband droop at f_signal is at most max_droop_db.

    For n stages |H(f)| = (1 + (f / f_c)^2)^(-n/2); the attenuation requirement gives the largest
    allowed corner f_c = f_reject / sqrt(attenuation^(2/n) - 1), and the droop is evaluated there.
    Raises ValueError if no n <= max_stages works."""
    raise NotImplementedError


# ---------------------------------------------------------------- the ideal op-amp

def noninverting_gain(Rf: float, Rg: float) -> float:
    """Ideal non-inverting amplifier: G = 1 + Rf / Rg (Rf from output to the inverting input, Rg from
    there to ground). Rf = 0 or Rg = infinity is the unity-gain buffer."""
    raise NotImplementedError


def inverting_gain(Rf: float, Rin: float) -> float:
    """Ideal inverting amplifier: G = -Rf / Rin (virtual ground at the inverting input)."""
    raise NotImplementedError


def closed_loop_gain(A: float, beta: float) -> float:
    """Negative-feedback amplifier with open-loop gain A and feedback fraction beta: G = A / (1 + A beta).
    For the non-inverting amplifier beta = Rg / (Rf + Rg); as A -> infinity, G -> 1 / beta."""
    raise NotImplementedError


# ---------------------------------------------------------------- diode and transistor (conceptual)

def thermal_voltage(T: float = 300.0) -> float:
    """V_T = k T / q, about 25.9 mV at 300 K."""
    raise NotImplementedError


def diode_current(V: np.ndarray | float, I_s: float = 1e-12, n: float = 1.0, T: float = 300.0) -> np.ndarray | float:
    """Shockley diode equation I = I_s (exp(V / (n V_T)) - 1): exponential in forward bias, -I_s in reverse."""
    raise NotImplementedError


def bjt_base_resistor(V_drive: float, I_c: float, beta: float, V_be: float = 0.7, overdrive: float = 5.0) -> float:
    """Base resistor for an NPN transistor used as a saturated switch: the base current is overdrive
    times the minimum I_c / beta, R_b = (V_drive - V_be) / (overdrive I_c / beta)."""
    raise NotImplementedError


# ---------------------------------------------------------------- logic and the ADC

def gate(name: str, a: int | np.ndarray, b: int | np.ndarray | None = None) -> int | np.ndarray:
    """The basic gates on 0/1 integers (or arrays): NOT (one input), AND, OR, NAND, NOR, XOR."""
    raise NotImplementedError


def half_adder(a: int, b: int) -> tuple[int, int]:
    """(sum, carry) = (a XOR b, a AND b): the one-bit adder built from two gates."""
    raise NotImplementedError


def aliased_frequency(f: np.ndarray | float, fs: float) -> np.ndarray | float:
    """Apparent frequency, in [0, fs/2], of a sinusoid of frequency f sampled at fs: |f - fs round(f / fs)|."""
    raise NotImplementedError


def adc_lsb(v_ref: float, bits: int) -> float:
    """Voltage of one least-significant bit of an ADC with full scale v_ref: v_ref / 2^bits."""
    raise NotImplementedError


def adc_convert(v: np.ndarray | float, v_ref: float = 3.3, bits: int = 12) -> np.ndarray:
    """What a microcontroller ADC returns: integer counts floor(v / LSB), clipped to [0, 2^bits - 1]
    (a mid-tread rounding converter differs by half an LSB; datasheets say which)."""
    raise NotImplementedError


def counts_to_volts(counts: np.ndarray | int, v_ref: float = 3.3, bits: int = 12) -> np.ndarray | float:
    """Back to volts at the centre of each code: (counts + 1/2) LSB."""
    raise NotImplementedError


def quantisation_noise_rms(v_ref: float, bits: int) -> float:
    """RMS error of an ideal quantiser for a busy signal: LSB / sqrt(12)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the anti-aliasing experiment (SYNTHETIC signals)

def synthetic_sensor_signal(t: np.ndarray, f_sig: float = 2.0, a_sig: float = 0.5, offset: float = 1.65,
                            f_vib: float = 430.0, a_vib: float = 0.05) -> np.ndarray:
    """SYNTHETIC sensor voltage for testing the chain, not a measurement: a slow signal of interest
    (offset + a_sig sin 2 pi f_sig t) plus a vibration pickup a_vib sin 2 pi f_vib t."""
    t = np.asarray(t, dtype=float)
    return offset + a_sig * np.sin(2 * np.pi * f_sig * t) + a_vib * np.sin(2 * np.pi * f_vib * t)


def filter_and_sample(t: np.ndarray, v: np.ndarray, R: float, C: float, stages: int, fs: float,
                      v_ref: float = 3.3, bits: int = 12) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Pass the finely sampled analog waveform v(t) through the buffered RC cascade (``scipy.signal.lsim``,
    starting from the steady DC state of v[0]; stages = 0 means no filter), then sample at fs and
    quantise. Returns (t_samples, v_filtered_at_samples, counts)."""
    raise NotImplementedError


def tone_amplitude(x: np.ndarray, fs: float, f: float) -> float:
    """Amplitude of the sinusoidal component of x at frequency f (least-squares fit of a sin and cos
    at f plus a constant): the measured alias tone in the samples."""
    n = np.arange(x.size)
    M = np.column_stack([np.ones(x.size), np.cos(2 * np.pi * f * n / fs), np.sin(2 * np.pi * f * n / fs)])
    c, *_ = np.linalg.lstsq(M, np.asarray(x, dtype=float), rcond=None)
    return float(np.hypot(c[1], c[2]))


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Run the lesson's numbers, print them and save rc_filter.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- the 400 Hz aircraft bus
    V, f_bus = 115.0, 400.0
    Rl, Ll = 10.0, 2e-3
    Z = series(impedance_R(Rl, f_bus), impedance_L(Ll, f_bus))
    P, Qv, pf = ac_power(V, Z)
    print(f"400 Hz bus, 115 V rms, load R = {Rl} ohm + L = {Ll * 1e3} mH: X_L = {Z.imag:.3f} ohm, |Z| = {abs(Z):.3f} ohm, "
          f"phase {np.degrees(np.angle(Z)):.1f} deg, I = {V / abs(Z):.3f} A rms, P = {P:.0f} W, Q = {Qv:.0f} var, pf = {pf:.3f}")
    tt = np.linspace(0, 1 / f_bus, 10001)[:-1]
    print(f"rms of a 115 V rms sine sampled over one period: {rms(np.sqrt(2) * V * np.sin(2 * np.pi * f_bus * tt)):.4f} V")

    # ---- series RLC resonance
    R, L, C = 10.0, 10e-3, 1e-6
    f0, Qf = resonant_frequency(L, C), quality_factor(R, L, C)
    ff = np.linspace(100.0, 10000.0, 990001)                      # 10 mHz grid
    Ir = np.abs(rlc_series_current(1.0, R, L, C, ff))
    above = ff[Ir >= Ir.max() / np.sqrt(2)]
    print(f"\nseries RLC R = {R} ohm, L = {L * 1e3} mH, C = {C * 1e6} uF: f0 = {f0:.1f} Hz, Q = {Qf:.2f}, f0/Q = {f0 / Qf:.1f} Hz; "
          f"swept peak at {ff[np.argmax(Ir)]:.1f} Hz, -3 dB width {above[-1] - above[0]:.1f} Hz")

    # ---- the anti-aliasing design
    fs, bits, vref = 100.0, 12, 3.3
    f_sig, f_vib, a_vib = 5.0, 430.0, 0.05
    lsb = adc_lsb(vref, bits)
    att = a_vib / lsb
    print(f"\nADC: {bits} bits, {vref} V: LSB = {lsb * 1e3:.4f} mV, quantisation noise {quantisation_noise_rms(vref, bits) * 1e3:.4f} mV rms")
    print(f"vibration {f_vib:.0f} Hz sampled at {fs:.0f} Hz appears at {aliased_frequency(f_vib, fs):.0f} Hz; "
          f"{a_vib * 1e3:.0f} mV = {att:.1f} LSB, so the filter must attenuate by {att:.1f} ({20 * np.log10(att):.1f} dB)")
    for n in (1, 2, 3):
        fc = f_vib / np.sqrt(att ** (2.0 / n) - 1.0)
        droop = -gain_db(rc_lowpass(f_sig, 1.0, 1.0 / (2 * np.pi * fc), n))
        print(f"  {n} stage(s): largest corner {fc:6.2f} Hz, droop at {f_sig:.0f} Hz {droop:.3f} dB")
    n_req, fc_req = required_stages_and_corner(f_sig, f_vib, att, 0.5)
    Rr, Cr = 3.3e3, 1e-6
    fc = rc_corner_frequency(Rr, Cr)
    H430 = abs(rc_lowpass(f_vib, Rr, Cr, 2))
    print(f"design: {n_req} buffered stages, corner <= {fc_req:.2f} Hz; chosen R = 3.3 kohm, C = 1 uF: f_c = {fc:.2f} Hz, "
          f"|H(430)| = {H430:.5f} ({gain_db(H430):.1f} dB), residual {a_vib * H430 * 1e3:.3f} mV = {a_vib * H430 / lsb:.2f} LSB, "
          f"droop at 5 Hz {-gain_db(rc_lowpass(f_sig, Rr, Cr, 2)):.3f} dB, delay {group_delay_dc(Rr, Cr, 2) * 1e3:.1f} ms")
    w = 2 * np.pi * np.array([f_sig, f_vib])
    _, Hs = signal.freqs(*[rc_lti(Rr, Cr, 2).num, rc_lti(Rr, Cr, 2).den], worN=w)
    print(f"scipy.signal.freqs check: |H| = {np.abs(Hs)[0]:.6f}, {np.abs(Hs)[1]:.6f}; ours {abs(rc_lowpass(f_sig, Rr, Cr, 2)):.6f}, {H430:.6f}")
    Hu = abs(rc2_unbuffered(f_vib, Rr, Cr))
    print(f"the same two sections without a buffer: |H(430)| = {Hu:.5f}, droop at 5 Hz {-gain_db(rc2_unbuffered(f_sig, Rr, Cr)):.3f} dB, "
          f"-3 dB at {np.interp(-3.0, -gain_db(rc2_unbuffered(np.linspace(1, 100, 10000), Rr, Cr))[::-1] * -1, np.linspace(1, 100, 10000)[::-1]):.2f} Hz")

    # ---- op-amp
    print(f"\nnon-inverting amplifier Rf = 27 k, Rg = 10 k: ideal G = {noninverting_gain(27e3, 10e3):.2f}; "
          f"with A = 1e5: {closed_loop_gain(1e5, 10e3 / 37e3):.5f}; A = 100: {closed_loop_gain(100.0, 10e3 / 37e3):.4f}")
    print(f"diode at 0.6 V, I_s = 1e-12 A: {diode_current(0.6) * 1e3:.3f} mA; at 0.66 V: {diode_current(0.66) * 1e3:.3f} mA "
          f"(V_T = {thermal_voltage() * 1e3:.2f} mV)")
    print(f"NPN switch for the 0.12 A torquer coil of P04.09 from a 3.3 V pin, beta = 100: R_b = {bjt_base_resistor(3.3, 0.12, 100):.0f} ohm")

    # ---- SYNTHETIC time-domain check
    t = np.arange(0.0, 4.0, 1.0 / 20000.0)
    v = synthetic_sensor_signal(t, f_sig=2.0, a_sig=0.5, f_vib=f_vib, a_vib=a_vib)
    print("\nSYNTHETIC test signal (not a measurement): 1.65 V + 0.5 V at 2 Hz + 50 mV at 430 Hz, sampled at 100 Hz, 12 bits")
    results = {}
    for n in (0, 1, 2):
        ts, vs, counts = filter_and_sample(t, v, Rr, Cr, n, fs)
        keep = ts >= 1.0                                            # skip the filter start-up
        amp = tone_amplitude(counts[keep].astype(float), fs, aliased_frequency(f_vib, fs))
        results[n] = (ts, vs, counts)
        pred = a_vib * (abs(rc_lowpass(f_vib, Rr, Cr, n)) if n else 1.0) / lsb
        print(f"  {n} stage(s): 30 Hz alias in the ADC counts {amp:7.2f} LSB (predicted {pred:6.2f} LSB)")

    # ---- figure
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.5))
    fr = np.logspace(0, 3.5, 400)
    ax = axs[0]
    for n in (1, 2, 3):
        ax.semilogx(fr, gain_db(rc_lowpass(fr, Rr, Cr, n)), label=f"{n} buffered stage(s)")
    ax.semilogx(fr, gain_db(rc2_unbuffered(fr, Rr, Cr)), "k--", label="2 stages, unbuffered")
    ax.axhline(-20 * np.log10(att), color="r", lw=0.8); ax.axvline(f_vib, color="r", lw=0.8)
    ax.set_xlabel("f (Hz)"); ax.set_ylabel("gain (dB)"); ax.set_title("RC low-pass, R = 3.3 kohm, C = 1 uF"); ax.legend()
    ax = axs[1]
    ax.semilogx(ff[::100], Ir[::100], label="|I| for 1 V")
    ax.axvline(f0, color="gray", lw=0.8)
    ax.set_xlabel("f (Hz)"); ax.set_ylabel("|I| (A)"); ax.set_title(f"series RLC, Q = {Qf:.0f}")
    ax = axs[2]
    for n in (0, 2):
        ts, vs, counts = results[n]
        ax.plot(ts, counts, ".-", ms=3, label=f"{n} stages")
    ax.set_xlim(1.0, 2.0); ax.set_xlabel("t (s)"); ax.set_ylabel("ADC counts")
    ax.set_title("SYNTHETIC signal sampled at 100 Hz"); ax.legend()
    fig.tight_layout()
    fig.savefig(out_dir / "rc_filter.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
