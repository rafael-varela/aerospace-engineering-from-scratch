# P09.05  Frequency Response & Bode Plots

> Drive a stable linear system with a sine and, once the transient has died, a sine comes out at the same frequency, scaled by |G(jω)| and shifted by ∠G(jω). A Bode plot is those two numbers drawn for every frequency, and a sine test is how you check it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the steady-state sinusoidal response of a stable LTI system from convolution, and explain why gains multiply (dB add) and phases add in a cascade. 2. Implement and test a Bode plotter from scratch: exact magnitude and phase from G(jω), the phase unwrapped and, independently, built from the angles of the pole and zero factors; checked against `scipy.signal.bode`. 3. Build the asymptotic (straight-line) Bode diagram from the Bode form, and quantify its error at real and complex corners (−3.01 dB, −20 log₁₀ 2ζ, 5.7° of phase). 4. Verify a frequency response by a simulated sine test (time-domain simulation, least-squares fit of the steady-state amplitude and phase), and read resonance, bandwidth, phase lag and gain crossover off the plot for an elevator servo, a gyro anti-aliasing filter and a pitch-attitude loop. |
| **Prerequisites** | `P09.04` (poles and time response, the servo design), `P09.03` (transfer functions, the pitch loop), `P09.02` (convolution, the RC front end), `P02.10` (Fourier transform; the frequency response as G on the imaginary axis), `P02.07` (phasors, complex gain), `P04.10` (RC corner frequency, dB). Skill check: you can compute the magnitude and angle of (1 + jω/a)⁻¹ and convert a ratio to dB. |
| **Engineering question** | The elevator servo of `P09.04` (designed for 10 % overshoot) must follow attitude commands from a pilot or an autopilot. Up to what frequency does it follow them, with what gain error and phase lag at 1 Hz and 2 Hz? How much does the gyro's RC front end attenuate a 1200 Hz vibration line, and what does it cost in phase at the pitch loop's frequencies? Where does the pitch-attitude loop of `P09.03` cross 0 dB, and with what phase? And does a simulated sine test reproduce the computed Bode plot? |
| **Assumptions and limits** | Linear, time-invariant, stable systems (the steady-state sine response exists only if every pole is in the open left half plane; for a marginally stable open loop such as a type-1 servo, Bode is read from G(jω) but cannot be measured open-loop). The sine test is a SIMULATION of a frequency-response test on the model, not a measurement; a real test sees noise, nonlinearity (backlash, saturation, rate limits) and amplitude dependence. The straight-line diagram is an approximation with known errors at the corners. |
| **Analytical work** | The steady-state response A|G(jω)| cos(ωt + φ + ∠G) from the convolution integral; the Bode form K_B s^{−type} Π(1 ± s/ω_b)^{±1} Π(1 + 2ζs/ω_n + s²/ω_n²)^{±1}; asymptotes and slopes (±20, ±40 dB/decade); corner errors −3.01 dB and −20 log₁₀(2ζ); the decade rule for phase and its 5.7° maximum error; resonance ω_r = ω_n√(1 − 2ζ²), M_r = 1/(2ζ√(1 − ζ²)); the −3 dB bandwidth of the second order; the gain crossover of K/(s(τs + 1)); minimum and non-minimum phase (an all-pass factor); phase as a time delay. |
| **Computational deliverable** | `starter/main.py`: `freq_response`, `mag_db`, `unwrap_deg`, `phase_from_factors`, `bode`, `bode_form`, `asymptotic_mag_db`, `asymptotic_phase_deg`, `bandwidth`, `gain_crossover`, `resonance`, `fit_sinusoid`, `sine_test` (`tf`, `tf_to_ss`, `simulate_tf`, `servo_closed_loop`, `rc_filter`, `pitch_open_loop` and `main` are given). Tests prove: magnitude and phase match `scipy.signal.bode` to 10⁻⁸ dB and 10⁻⁷ deg and G(jω) matches `scipy.signal.freqs` for five systems (servo, two RC stages, the type-1 pitch loop, a lightly damped mode, a non-minimum-phase system); the factor-sum phase equals the unwrapped phase; unwrapping removes only 360° jumps; the simulated sine test reproduces |G| to 10⁻⁵ and ∠G to 10⁻³ deg at three servo frequencies, at the resonance of a ζ = 0.05 mode and on the RC filter at 1200 Hz; the least-squares fit is exact without noise and within statistical bounds with seeded noise; low- and high-frequency limits and slopes; the corner corrections −20 log₁₀ 2ζ, −3.0103 dB, −0.0432 dB a decade away and 5.71° maximum phase error; straight-line slopes from the Bode form; an all-pass has 0 dB and −2 arctan ω; resonance and bandwidth formulas against numerical peaks; the servo's gain crossover against its closed form. |
| **Practical activity** | **computational experiment**: a SIMULATED sine-sweep test of the elevator servo compared with the computed Bode plot; the gyro RC filter at the pitch-loop frequencies and at the vibration line; straight-line against exact Bode for the pitch loop; the transient-length study of the sine test. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `bode.png`, the derivation of the steady-state response and of the corner corrections, the table of computed against SIMULATED sine-test gain and phase (labelled), the straight-line error figures, and a short note answering the engineering question with numbers. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the steady-state response and the corner corrections are derived; all public tests pass; the computed and simulated responses agree and the residual is explained; the SIMULATED test is labelled as a simulation everywhere. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 9 "Frequency Response" (Eigenfunctions; Complex Exponentials; Rational System Functions; Vector Diagrams; Example: Mass, Spring, and Dashpot) and Lecture 11 "Continuous-Time (CT) Frequency Response and Bode Plot" (Frequency Response: H(s)\|s←jω; Asymptotic Behavior: Isolated Zero; Asymptotic Behavior: Isolated Pole; Bode Plot: Adding Instead of Multiplying; Bode Plot: dB; Bode Plot: Accuracy; Frequency Response of a High-Q System) (required-free): G(jω) as the steady-state complex gain on a sinusoid, why dB and log frequency make factors add, break-frequency asymptotes and their worst-case error, and the resonant peak of a lightly damped pair. `@mit-ocw-16-06-principles-of-automatic-control` Lecture 16 "Frequency response design" (The Frequency Response; Bode Diagram), Lecture 17 "Bode plot problems" (PDF "Frequency Response Design": Bode Plot Construction) and Lecture 18 "Bode plots (cont.)" (required-free): building a Bode diagram by hand, factor by factor. `@mit-ocw-16-30-feedback-control-systems` Topic 3 "Frequency response methods": Frequency Response Function; Bode Example (required-free): one more worked example. `@scipy-userguide` API reference "Signal processing (scipy.signal)": "Continuous-time linear systems" (`bode`, `freqs`) (required-free, oracle only). |
| **Extensions** | ★ Add a pure time delay e^{−sT} to the plotter (magnitude 1, phase −ωT) and compare with its first-order Padé approximation (Problem 6). ★ Replace the single-frequency sine test by a SYNTHETIC chirp and estimate the frequency response from the ratio of the FFTs of output and input (`P02.10`); compare its accuracy and test time with the sine sweep. ★ Implement the Nichols chart (magnitude against phase) and plot the pitch loop on it. |
| **Time** | unmeasured |

## 1. The engineering problem

Every actuator and sensor in a flight-control loop is specified in the frequency domain. An elevator servo is said to have a "bandwidth of 3 Hz"; a gyro front end "attenuates 1200 Hz by 22 dB"; a pitch loop "crosses over at 5.7 rad/s with 50° of phase margin". The numbers come from the **frequency response**: how a stable linear system scales and delays a sinusoid at each frequency. On the ground the frequency response is *measured* by a sine test: drive the actuator with a sinusoid at one frequency, wait for the transient to die, record the amplitude and phase of the output, step to the next frequency. In analysis it is *computed* from the transfer function, as G(jω).

This lesson builds the computation from scratch and checks it against the test, on the model. The examples continue the earlier lessons: the elevator servo K/(τs² + s + K) with K = 14.31 s⁻¹ and τ = 50 ms (the 10 % overshoot design of `P09.04`); the gyro's RC front end with τ = 1.6 ms (`P09.02`); and the return ratio of the pitch-attitude loop (`P09.03`), whose 0 dB crossing sets up the stability margins of `P09.10`.

Success criteria, written first: the computed magnitude and phase match `scipy.signal.bode` to 10⁻⁸ dB and 10⁻⁷ deg; the phase is continuous and starts at its physical low-frequency value (−90° for a type-1 loop); a simulated sine test reproduces the computed gain and phase to better than 10⁻³ dB and 10⁻² deg at every test frequency; and the straight-line diagram's errors at the corners are the ones derived by hand.

## 2. Understand and derive

### 2.1 The steady-state response to a sinusoid

Let a stable LTI system with impulse response h be driven by x(t) = Re{A e^{jωt}} from t = 0. By convolution (`P09.02`),

y(t) = Re{A ∫₀ᵗ h(τ) e^{jω(t − τ)} dτ} = Re{A e^{jωt} [G(jω) − ∫_t^∞ h(τ) e^{−jωτ} dτ]},

with **G(jω) = ∫₀^∞ h(τ) e^{−jωτ} dτ**, the Laplace transform on the imaginary axis (`P02.10`; 6.003 Lecture 11, "Frequency Response: H(s)|s←jω"). The last integral is the transient: for a stable system h decays exponentially, so it decays too, at the rate of the slowest pole. What remains is

**y_ss(t) = A |G(jω)| cos(ωt + ∠G(jω))**.   (1)

The sinusoid keeps its frequency; its amplitude is multiplied by |G(jω)| and its phase is shifted by ∠G(jω). For an unstable system the transient grows and (1) is never reached, which is why a frequency response is measured on a stable system or inside a stabilising loop.

In a cascade G = G₁G₂, |G| = |G₁||G₂| and ∠G = ∠G₁ + ∠G₂. In **decibels**, 20 log₁₀|G| = 20 log₁₀|G₁| + 20 log₁₀|G₂|: magnitudes add, and so do phases. That is why a Bode plot (magnitude in dB and phase in degrees, both against log ω) can be drawn factor by factor.

### 2.2 Phase: unwrapping and the factor sum

The phase of a complex number from `np.angle` lies in (−180°, 180°], so a phase that really falls through −180° (any system with three or more poles) jumps by 360° on a sweep. Two cures, implemented independently and tested against each other:

- **Unwrapping**: walk along the sweep and remove every jump larger than 180° between neighbours by adding a multiple of 360°. This needs a sweep dense enough that the true phase never changes by 180° between two points, and it leaves a 360° ambiguity at the start.
- **Factor sum**: write G(s) = C s^{n₀} Π(1 − s/z_i)/Π(1 − s/p_i) over the non-zero zeros and poles, with n₀ the net number of zeros at the origin and C real. Then **∠G = (180° if C < 0) + 90° n₀ + Σ∠(1 − jω/z_i) − Σ∠(1 − jω/p_i)**. Each factor starts at 0° at ω = 0 and moves continuously, so the sum is the physical phase with no ambiguity. `bode` unwraps the angle of G(jω) and then shifts it by the multiple of 360° that makes its first value agree with the factor sum; this reproduces `scipy.signal.bode` exactly for all five test systems.

### 2.3 The Bode form and the straight-line diagram

Collect the factors in **Bode form** (6.003 Lecture 11; 16.06 Lecture 17, "Bode Plot Construction"):

G(s) = K_B s^{−type} Π(1 + s/ω_b)^{±1} Π(1 + 2ζs/ω_n + s²/ω_n²)^{±1},   (2)

where K_B is the gain of the low-frequency asymptote (the DC gain for type 0; lim s^{type}G(s) otherwise). For ω much below a corner a factor is ≈ 1 (0 dB, 0°); much above, (1 + jω/ω_b) ≈ jω/ω_b, a line rising at **+20 dB/decade** with +90°, and a quadratic factor ≈ −ω²/ω_n², ±40 dB/decade with ±180°. The **straight-line magnitude** is the sum of these asymptotes, starting at 20 log₁₀|K_B| − 20·type·log₁₀ω. Its error is largest at the corners:

- a real pole at its corner: |1/(1 + j)| = 1/√2, **−3.01 dB**; one decade away, −0.043 dB;
- a complex pair at ω_n: |1/(2jζ)| = 1/(2ζ), an error of **−20 log₁₀(2ζ)** dB: +20 dB at ζ = 0.05, 0 dB at ζ = 0.5.

The **straight-line phase** of a real factor is 0° below ω_b/10, ±90° above 10ω_b and a straight line in log ω between (45° per decade); its largest error is **5.71°**, a decade either side of the corner (arctan 0.1). For a complex pair the phase changes by 180° over a band whose width shrinks with ζ; the code uses ω_n 10^{±ζ}, a step as ζ → 0. A right-half-plane zero has the magnitude of its mirror image but the phase of a pole: magnitude alone cannot tell them apart.

### 2.4 Resonance and bandwidth of the second order

For ω_n²/(s² + 2ζω_n s + ω_n²), |G|² = 1/((1 − u²)² + 4ζ²u²) with u = ω/ω_n. Setting its derivative to zero gives the **resonant peak**

**ω_r = ω_n√(1 − 2ζ²),  M_r = 1/(2ζ√(1 − ζ²))**   (ζ < 1/√2; no peak above),   (3)

and setting |G|² = ½ gives the **−3 dB bandwidth** ω_B = ω_n[1 − 2ζ² + √(4ζ⁴ − 4ζ² + 2)]^{1/2}. For the servo (ζ = 0.5911, ω_n = 16.92 rad/s): ω_r = 9.284 rad/s, M_r = 1.0487 (**+0.41 dB**), bandwidth **19.62 rad/s = 3.12 Hz**. The resonant peak is the frequency-domain face of the 10 % overshoot: small damping means a tall peak and a ringing step.

### 2.5 Phase as delay, and the servo's lag

A phase lag φ at frequency ω is a time delay φ/ω. For the servo at 1 Hz, ∠G = −27.0°, a delay of 75 ms. At low frequency, −dφ/dω is the group delay of `P04.10`, which for a first-order lag is τ and for the servo is 1/K = 70 ms (the error of a type-1 loop following a ramp, `P09.08`). A command that arrives 75 ms late inside an outer loop costs that loop phase margin (`P09.10`).

### 2.6 The gain crossover of a loop

For a loop transfer function L(s), the **gain crossover** ω_c is where |L(jω)| = 1 (0 dB). For the servo's open loop K/(s(τs + 1)), |L|² = 1 gives ω_c² (1 + ω_c²τ²) = K², hence ω_c = √[(−1 + √(1 + 4K²τ²))/(2τ²)] = 12.21 rad/s, with phase −90° − arctan(ω_cτ) = −121.4°. The distance of that phase from −180° is the **phase margin** (58.6° here), the subject of `P09.10`; this lesson only computes the crossover.

### 2.7 Where the model stops being valid

The frequency response is a property of a linear model. A real servo tested with a small sine sees backlash and friction (the output may not move at all below some amplitude); with a large sine it hits its rate limit, and the "gain" falls with amplitude. A real test therefore states its amplitude, and a careful one repeats at two amplitudes. Above the frequencies the model was built for, unmodelled dynamics take over: the motor's electrical time constant, structural modes of the surface and its linkage. The straight-line diagram is a sketching tool; for numbers, compute G(jω).

## 3. Build the computational model

Open `starter/main.py`. `tf`, `tf_to_ss`, `simulate_tf` (from `P09.03`), the three example systems and `main` are given.

### Step 1: the exact response

```python
def freq_response(G, w):
    s = 1j * np.asarray(w, dtype=float)
    return np.polyval(G[0], s) / np.polyval(G[1], s)

def unwrap_deg(p):
    d = np.diff(p)
    p = p.copy()
    p[1:] += np.cumsum(-360.0 * np.round(d / 360.0))
    return p
```

`mag_db` is 20 log₁₀|H|. `phase_from_factors` builds the factor sum of 2.2 from `np.roots` of the numerator and denominator, counting roots at the origin separately and accumulating the sign of C. `bode` combines them as described.

### Step 2: the Bode form and the straight lines

`bode_form` sorts the roots into origin roots (the type), real corners (with a flag for the half plane) and conjugate pairs (ω_n, ζ), and accumulates K_B: a real root r contributes (−r)^{±1}, a pair |r|^{±2}. `asymptotic_mag_db` and `asymptotic_phase_deg` sum the asymptotes of 2.3.

### Step 3: features

`resonance` is (3); `bandwidth` scans a log grid for the first frequency 3.01 dB below the DC gain and bisects; `gain_crossover` does the same for 0 dB and returns the factor-sum phase there.

### Step 4: the simulated sine test

```python
def fit_sinusoid(t, y, w):
    M = np.column_stack([np.cos(w * t), np.sin(w * t), np.ones_like(t)])
    (a, b, _), *_ = np.linalg.lstsq(M, y, rcond=None)
    return np.hypot(a, b), np.degrees(np.arctan2(-b, a))
```

`sine_test(G, w)` simulates G from rest with input cos ωt for a settling time of 10/|Re p_slow| (the transient decays by e^{−10}) plus five periods, with at least 200 points per period and a step no larger than 0.5/|p_fast| for RK4 accuracy, and fits the last five periods. The constant column of the fit absorbs any residual offset.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-05-frequency-response-bode`. The tests:

1. **Against SciPy**: `bode` matches `scipy.signal.bode` (magnitude 10⁻⁸ dB, phase 10⁻⁷ deg) and `freq_response` matches `scipy.signal.freqs` on 3001 frequencies from 10⁻² to 10⁴ rad/s, for the servo, two RC stages, the type-1 pitch loop, a ζ = 0.05 mode and a non-minimum-phase system.
2. **Factor-sum phase** equals the unwrapped phase on the same five systems.
3. **Unwrapping** recovers a wrapped ramp to −700° and leaves a smooth phase alone.
4. **Simulated sine test**: the servo at 0.5, 2 and 6 Hz (gain 10⁻⁵ relative, phase 10⁻³ deg); the ζ = 0.05 mode at its resonance (gain M_r = 10.01 to 10⁻⁴); the RC filter at 1200 Hz.
5. **Fit**: exact for a noise-free sine with offset; within five standard deviations with seeded noise.
6. **Limits**: 0 dB and 0° at low frequency for the servo, −40 dB/decade and −180° at high frequency; the pitch loop's −20 dB/decade and −90° at low frequency and its K_B = 108 000/42 800.
7. **Corners**: −20 log₁₀(2ζ) at ω_n for ζ = 0.05, 0.2, 0.5; −0.0432, −3.0103, −0.0432 dB around a real corner; 5.71° maximum phase error.
8. **Slopes** of the straight lines for the pitch loop: −20 and −60 dB/decade, phases −90° and −270°.
9. **All-pass**: (1 − s)/(1 + s) is 0 dB with phase −2 arctan ω; an RHP zero has the magnitude of its mirror and more lag.
10. **Resonance and bandwidth** against a fine numerical search for ζ = 0.1, 0.3, 0.6; the first-order bandwidth 1/τ; the servo's crossover against its closed form.
11. **main** writes the figure.

### Compare

`scipy.signal.bode` evaluates G(jω) with `freqs` and unwraps the angle with `numpy.unwrap`; for every test system our phase and SciPy's agree to 10⁻⁷ deg, because both start from the low-frequency value in (−180°, 180°]. The factor sum is the more robust of the two methods: it needs no dense grid and has no 360° ambiguity, which matters for lightly damped systems where the phase drops by 180° within a few percent of ω_n. The simulated sine test agrees with both to 10⁻⁴ dB or better; its residual is the RK4 error and the e^{−10} transient that has not quite died, and it can be made smaller by waiting longer (part 4). A laboratory frequency-response analyser does the same fit on a measured signal, usually by correlating with sin and cos over whole periods, which is the least-squares fit of Step 4 in disguise.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Every sine test here is a SIMULATION on the model, not a measurement. Write each prediction first.

1. **Servo sine sweep (SIMULATED).** Predict |G| and ∠G of the servo at 0.5, 1, 2, 4 and 8 Hz from (3) and the Bode plot; then run `sine_test` at each and tabulate. Predict the frequency of the resonant peak, the −3 dB bandwidth and the phase lag at 1 Hz as a time delay.
2. **The gyro RC filter.** For one and two buffered stages with τ = 1.6 ms, predict the magnitude and phase at 10 Hz (inside the pitch loop's band), at the 99.5 Hz corner and at the 1200 Hz vibration line, exactly and from the straight lines. Check one stage at 1200 Hz with a SIMULATED sine test.
3. **Pitch loop.** From the Bode form of L(s) of `P09.03` list the corners and K_B; predict where the straight-line magnitude crosses 0 dB; then compute the exact crossover and its phase, and the largest straight-line errors in magnitude and phase from 0.1 to 1000 rad/s.
4. **How long to wait.** For the closed pitch loop (slowest pole −1.068 s⁻¹) at ω = 1 + s/10 rad/s, run `sine_test` with settle = 0, 1, 3 and 10 s and three fitted periods; record the gain and phase errors. Predict how the error falls with the settling time.

Raw data format: CSV tables `f_Hz, mag_exact_dB, phase_exact_deg, mag_sim_dB, phase_sim_deg` with SIMULATED in the header (1); `stages, f_Hz, mag_dB, phase_deg, mag_straight_dB` (2); `corner_rad_s, kind` and `w_c, phase_at_w_c, max_mag_err_dB, max_phase_err_deg` (3); `settle_s, gain_err_dB, phase_err_deg` (4).

What you should find. (1)

| f (Hz) | \|G\| exact (dB) | ∠G exact (deg) | SIMULATED test (dB, deg) |
|---|---|---|---|
| 0.5 | +0.0859 | −12.810 | +0.0859, −12.810 |
| 1 | +0.2875 | −26.991 | +0.2875, −26.991 |
| 2 | +0.1230 | −62.959 | +0.1230, −62.959 |
| 4 | −6.5720 | −124.499 | −6.5721, −124.499 |
| 8 | −18.6699 | −155.834 | −18.6699, −155.834 |

The servo *amplifies* slightly up to its resonance at 9.28 rad/s (1.48 Hz, +0.41 dB), has a −3 dB bandwidth of 3.12 Hz, and lags by 27.0° at 1 Hz, a 75 ms delay. It follows 1 Hz commands within 0.3 dB but not in phase. (2) One stage: −0.044 dB and −5.7° at 10 Hz; −3.010 dB and −45.0° at the corner; −21.66 dB and −85.3° at 1200 Hz (the straight line says −21.63 dB). Two stages: −0.087 dB and −11.5° at 10 Hz, −6.02 dB at the corner, −43.3 dB at 1200 Hz. The SIMULATED 1200 Hz test gives −21.659 dB, −85.26°. The 1200 Hz line falls by a factor of 12 per stage; the price at 10 Hz is 5.7° of phase per stage, which the pitch loop pays. (3) K_B = 2.523, type 1; corners: zeros at 1.5 and 100 rad/s, poles at 13.59 and 100.9 rad/s, and a pole pair at ω_n = 5.586 rad/s with ζ = 0.780. The exact crossover is 5.736 rad/s with phase −129.45°, so 180° + phase = 50.6°. The straight lines are off by up to 4.22 dB and 22.6°, the large phase error coming from the heavily damped pair whose real phase changes over a much wider band than any straight-line rule assumes. (4) At ω = 2.5 rad/s the exact response is −1.603 dB, −26.96°; fitting three periods with no settling time gives errors of −0.32 dB and −0.84°; settling 1 s: +0.006 dB, +0.048°; 3 s: +0.0013 dB, +0.0099°; 10 s: below 10⁻⁴. The error falls like e^{−1.068·t_settle}, the slowest pole's decay.

## 5. Validate and improve

| quantity | theory (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| servo \|G\|, ∠G at 1 Hz | +0.2875 dB, −26.991° | same | SIMULATED sine test: same to 10⁻⁴ | < 10⁻⁴ | transient e^{−10}, RK4 |
| servo resonance | (3): 9.284 rad/s, +0.413 dB | numerical peak | test at ζ = 0.05 in the suite | < 10⁻³ rad/s | grid |
| servo bandwidth | closed form 19.62 rad/s | `bandwidth` 19.624 rad/s | test against the ζ formula | < 10⁻⁷ relative | — |
| RC at its corner | −3.0103 dB, −45° | −3.010 dB, −45.00° | `scipy.signal.bode` | 10⁻⁸ | — |
| RC at 1200 Hz | −21.659 dB | −21.659 dB | SIMULATED test −21.659 dB | < 10⁻³ | — |
| pitch-loop crossover | ω_c of \|L\| = 1 | 5.736 rad/s, −129.45° | `scipy.signal.bode` crossing | < 10⁻⁶ | bisection |
| straight lines, pitch loop | — | max errors 4.22 dB, 22.6° | exact Bode | large near the damped pair | approximation |

Computation, oracle and simulated test agree to the precision each claims. The disagreement worth reporting is the straight-line diagram's for the pitch loop: 22.6° of phase is the difference between a comfortable and a marginal loop, so the asymptotes are for sketching and checking the shape, never for margins. The improvement for a real test campaign is the one 2.7 names: run the sine test at two amplitudes and look for amplitude dependence before trusting any linear model; and choose the settling time from the slowest pole, which part 4 measured.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1), of the corner errors and of (3); `bode.png`; the four experiment tables (the sine-test ones labelled SIMULATED); and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the servo's bandwidth and its gain and phase at 1 and 2 Hz, the gyro filter's attenuation at 1200 Hz and its phase cost at 10 Hz, and the pitch loop's crossover.

### Adversarial self-checks

- **What happens as ω → 0 and ω → ∞?** Type 0: |G| → |G(0)|, phase → 0° (or 180° for negative gain). Type 1: |G| grows as 1/ω and the phase is −90°. At high frequency the slope is −20 dB/decade per unit of relative degree and, for a minimum-phase system, the phase tends to −90° per unit.
- **What if the servo gain doubles?** ζ falls to 0.418, the resonant peak rises to M_r = 1/(2ζ√(1 − ζ²)) = 1.32 (+2.4 dB) and the bandwidth grows; more overshoot in time is more peak in frequency.
- **Which assumption fails first in a real sine test?** Amplitude: backlash dominates small amplitudes, and the rate limit large amplitudes at high frequency (the surface rate is amplitude × ω × |G|). Then noise, which the least-squares fit averages down only as 1/√N.
- **What measurement could falsify the model?** A sine test on the real servo whose phase at 4 Hz is far below −124.5° would reveal unmodelled lag (the motor's electrical time constant, a filter); a gain that changes with amplitude reveals a nonlinearity.
- **Does the result survive a change of method?** G(jω), `scipy.signal.bode`, the factor sum and the SIMULATED sine test agree to 10⁻⁴ dB or better. The straight lines do not, and are not meant to.
- **Which decision changes if the gyro filter's capacitor is 20 % high?** The 1200 Hz attenuation improves by about 1.6 dB, but the phase lag at 10 Hz grows from 5.7° to 6.9° per stage, phase taken from the pitch loop's margin.

## Check yourself

1. A stable system has G(j2) = 0.5 e^{−j60°}. What is its steady-state output for the input 3 cos(2t)?
2. Why do magnitudes in dB and phases add in a cascade?
3. What are the magnitude and phase of a single real pole at its corner frequency, and how far is the straight line off there?
4. A lightly damped mode has ζ = 0.05. How far above the straight line is the true magnitude at ω_n?
5. What is the low-frequency slope and phase of a type-1 loop, and what is its Bode gain K_B?
6. Why can a frequency response not be measured open-loop on the servo K/(s(τs + 1)), and how do engineers measure it anyway?
7. The servo lags by 27° at 1 Hz. What time delay is that?
8. A sine test on a system with slowest pole −1 s⁻¹ fits the output after only 1 s of settling. Roughly what fraction of the transient remains?
