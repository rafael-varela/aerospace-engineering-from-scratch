# P06.17  Vibration Isolation

> An isolator only starts to work above √2 times its own natural frequency, it makes things worse below that, and the softer it is the further the box swings. Isolator design is choosing where to put that resonance and how much to damp it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the base-excitation transmissibility of a spring–damper isolator, its phase, the √2 crossover, the exact resonant peak and its high-frequency roll-off, and compare viscous with hysteretic (elastomer-like) damping. 2. Size an isolator: the stiffness for a target natural frequency, the damping for an allowed resonant peak, and the corner frequency needed for a required attenuation at a given frequency. 3. Implement and test the transmissibility against a `solve_ivp` steady state, and the random response of the isolated box (transmitted Grms and sway) by the integral of `P06.16`, checked against Miles' equation and a seeded SYNTHETIC time-domain simulation. 4. Run a design search over natural frequency and damping against four competing requirements, and decide which requirement sets the design and how sensitive the choice is to the mass and the elastomer's stiffness. |
| **Prerequisites** | `P06.16` (PSD, |H|²G, Miles' equation, the time-domain check), `P06.13` (SDOF base excitation, damping ratio), `P02.07` (phasors and complex gain), `P03.02` (root finding with `brentq`), `P02.19` (the variance of an RMS estimate). Skill check: you can write the steady-state amplitude of a damped oscillator driven through its base. |
| **Engineering question** | A 0.35 kg MEMS inertial measurement unit (IMU) sits on four elastomer mounts on a launch vehicle's avionics shelf, which vibrates at the NASA-STD-7001C minimum-workmanship level (6.78 Grms). The mount mode must stay at or above 60 Hz, the static plus 3σ sway must stay within 0.8 mm, the resonant peak must not exceed 4 and the transmissibility at 500 Hz must not exceed 0.06 (illustrative requirements). Which stiffness and damping minimise the vibration the IMU feels, and which requirement decides? The answer: fₙ = 65 Hz, ζ = 0.20 (58.4 kN/m and 57.2 N·s/m in total), 2.88 Grms transmitted (42 % of the input), 0.757 mm sway, a peak of 2.73 and T(500 Hz) = 0.056; the sway and the 500 Hz requirement bind. |
| **Assumptions and limits** | The IMU is a rigid mass on a linear isolator moving along one axis (real mounts couple translation and rotation, giving six modes); the shelf is rigid and its motion is the specified PSD; viscous damping, or hysteretic damping with a constant loss factor; stiffness independent of frequency, amplitude and temperature (real elastomers stiffen with frequency and with cold, by tens of per cent to a factor of two or more); Gaussian stationary input. The IMU mass and the four requirements are illustrative, chosen to make the trade visible. Random records are SYNTHETIC. |
| **Analytical work** | Transmissibility and phase from the equation of motion; the √2 crossover for every damping; the peak r_p² = (√(1 + 8ζ²) − 1)/(4ζ²); asymptotes 2ζ/r (viscous) and 1/r² (hysteretic); stiffness from frequency and from static sag, fₙ = (1/2π)√(g/δ); the corner frequency for a target attenuation; the random response and sway by the `P06.16` integrals; the sway budget (static plus 3σ). |
| **Computational deliverable** | `starter/main.py`: `natural_frequency`, `stiffness_for_frequency`, `damping_coefficient`, `static_deflection`, `transmissibility`, `transmissibility_phase`, `transmissibility_hysteretic`, `peak_transmissibility`, `zeta_for_peak`, `fn_for_attenuation`, `steady_state_ivp`, `random_response`, `design_sweep`, `choose_isolator` (the `P06.16` helpers `psd_loglog`, `psd_grms`, `freq_grid`, `synthesize_random`, `simulate_sdof_base`, `miles_grms`, and `main`, are given). Tests (17) prove: stiffness, frequency, damping and static sag are consistent; the transmissibility is 1 at r → 0 and at r = √2 for every ζ, below 1 beyond, 1/|1 − r²| undamped; the phase limits; −20 dB/decade viscous and −40 dB/decade hysteretic roll-off; damping lowers the peak and raises high-frequency transmission; the exact peak against a brute-force search and the 1/(2ζ) limit; the inverse design functions; the closed form against `solve_ivp` at five (r, ζ) pairs to 10⁻⁶; the random response's rigid and soft limits, Miles' limit, the white-noise sway formula and a seeded time-domain simulation; the design sweep's trends; and a feasible optimum that is the best feasible point on the grid. |
| **Practical activity** | **computational experiment**: the transmissibility against `solve_ivp`; the damping trade at resonance and at r = 10, viscous against hysteretic; a sweep of natural frequency and damping against the four requirements; the sensitivity of the chosen design to mass, dynamic stiffening and temperature. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vibration_isolation.png`, the derivation of the transmissibility and its peak, the design table, the sensitivity table, and a half-page design note stating the chosen k and c, which requirement binds, and what you would test first. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: transmissibility and its peak derived; all public tests pass; the closed form is checked against `solve_ivp` and the random response against the time domain with a stated tolerance; the note names the binding requirements and the sensitivity to elastomer stiffness. |
| **References** | `@mit-ocw-2-003sc-engineering-dynamics` Lecture 21: Vibration Isolation (the transfer function with complex numbers, isolation system design, natural frequency from static deflection √(g/δ), isolating the source); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response (required-free). `@nasa-std-7001c-payload-vibroacoustic-test-criteria` 4.3.3 Workmanship (Table 2, the input spectrum for an IMU-sized component); 4.4.6 Test Tailoring Methods (required-free). `@mit-ocw-18-03-differential-equations` notes 14 "Complex gain, dashpot drive" and 16 "Frequency response" (required-free): the complex-gain method used in 2.2. `@mit-ocw-6-011-signals-systems-and-inference` Lecture 17: LTI Filtering of WSS Processes (supplementary-free): the random response. Sway space and IMU isolator sizing are original. Optional textbooks for deeper study (paid, never required; try a library): Inman, *Engineering Vibration*, the chapter on vibration isolation. |
| **Extensions** | ★ Add the rotational modes: put the four mounts at the corners of the IMU's base with the centre of mass above them, and find the coupled rocking and translation frequencies (`P06.14`); explain why isolator layouts aim to put the elastic centre at the centre of mass. ★ Model the elastomer's dynamic stiffening as k(f) = k₀(1 + 0.3 log₁₀(f/10 Hz)) and redo the design. ★ Design a two-stage isolator (an intermediate mass between two mount sets) and show its −80 dB/decade roll-off and its extra resonance. |
| **Time** | unmeasured |

## 1. The engineering problem

An inertial measurement unit measures the vehicle's acceleration and rotation rate. Its gyroscopes and accelerometers are small vibrating or suspended structures with their own resonances, and the electronics that read them assume the input changes slowly. Shake an IMU hard at high frequency and two things go wrong: its sensors *rectify* the vibration into an apparent steady bias (a navigation error that grows with time), and in the worst case a sensor's own resonance is excited and the reading saturates. So IMUs on launch vehicles, missiles and helicopters often sit on **isolators**: small elastomer or wire-rope mounts that put a soft spring between the vehicle and the box.

An isolator is the base-excited oscillator of `P06.13` and `P06.16` used on purpose. Above √2 times its natural frequency it transmits less than it receives, and the higher the frequency the better it isolates. But it creates a resonance, which amplifies whatever arrives near fₙ, and it lets the box move relative to the shelf: cables and connectors must allow that sway, and the box must not hit its neighbours. It also changes what the IMU measures: below fₙ the box follows the vehicle, above it the box does not, so fₙ must sit well above the frequencies the guidance system needs to see. Every choice of stiffness and damping is a trade among these.

This lesson builds the tools for that trade and applies them to a 0.35 kg IMU on the NASA-STD-7001C workmanship input. Success criteria, written first: the closed-form transmissibility agrees with a `solve_ivp` steady state to 10⁻⁶; the random response agrees with a seeded time-domain simulation within four standard deviations of its statistical scatter; the design search returns the best design that meets every requirement, and the note says which requirements bind.

## 2. Understand and derive

### 2.1 The isolator and its two numbers

Let the shelf move with displacement y(t) and the IMU, of mass m, with x(t). Between them are the mounts, with total stiffness k and viscous damping c (four mounts in parallel add their stiffnesses and dampings). Newton's second law for the IMU in an inertial frame:

**m ẍ = −k(x − y) − c(ẋ − ẏ).**

Two numbers describe the mount: the natural frequency fₙ = (1/2π)√(k/m) and the damping ratio ζ = c/(2√(km)). Inverting, **k = m(2πfₙ)²** and **c = 2ζ√(km)**. A 0.35 kg IMU at fₙ = 65 Hz needs k = 0.35 × (2π × 65)² = 58.4 kN/m in total, 14.6 kN/m per mount.

A second route to fₙ is the static sag (2.003SC Lecture 21). Under 1 g the mounts deflect δ = mg/k, so ωₙ² = k/m = g/δ and **fₙ = (1/2π)√(g/δ)**: a 65 Hz mount sags 0.059 mm under the IMU's weight. Under a steady acceleration of n g (the launch vehicle's thrust, 5 g here) the sag is n g/ωₙ² = 0.294 mm. That is part of the sway budget.

### 2.2 Transmissibility

Put y = Y e^{iωt} and x = X e^{iωt} (the complex-gain method, 18.03 notes 14 and 16). Then (−mω² + icω + k)X = (icω + k)Y, and dividing by k with r = ω/ωₙ:

X/Y = (1 + 2iζr)/(1 − r² + 2iζr).

Its magnitude is the **transmissibility**

**T(r) = √[ (1 + (2ζr)²) / ((1 − r²)² + (2ζr)²) ],**

and its phase lag is φ = −arg(X/Y). Because ẍ/ÿ = X/Y as well, the same T relates accelerations, and the same T relates the force transmitted to a rigid base to an applied force (isolation of the source, 2.003SC Lecture 21). T² is the |H|² of `P06.16`.

Read it:

- **r → 0:** T → 1, φ → 0. The box follows the shelf.
- **r = 1:** T = √(1 + 4ζ²)/(2ζ) ≈ 1/(2ζ) = Q, φ = arctan(1/(2ζ)). With ζ = 0.1, T = 5.10 and the box lags 78.7°.
- **r = √2:** put r² = 2: the numerator and denominator are both 1 + 8ζ², so **T = 1 for every damping**. Below this frequency the isolator amplifies; above it, it isolates.
- **r ≫ 1:** T ≈ √(1 + 4ζ²r²)/r² → 2ζ/r. The viscous isolator rolls off at only −20 dB per decade, because the damper connects the box directly to the shelf and transmits a force cω Y that grows with frequency.

The undamped isolator, T = 1/|1 − r²|, rolls off at −40 dB/decade but is infinite at resonance. That is the **damping trade**: damping limits the peak and spoils the high-frequency isolation.

### 2.3 The exact peak and the damping for an allowed peak

Setting dT²/dr² = 0 (write u = r²) gives a quadratic, 2ζ²u² + u − 1 = 0, whose positive root is

**r_p² = (√(1 + 8ζ²) − 1)/(4ζ²).**

The peak sits slightly below r = 1 (r_p = 0.9903 at ζ = 0.1, 0.9647 at ζ = 0.2) and tends to 1 as ζ → 0, where T_peak → 1/(2ζ). Inverting T_peak(ζ) numerically with `brentq`, a peak of 4 needs **ζ = 0.130**.

| ζ | 0.02 | 0.05 | 0.10 | 0.20 | 0.40 |
|---|---|---|---|---|---|
| T_peak | 25.02 | 10.06 | 5.12 | 2.73 | 1.66 |
| T(r = 10), viscous | 0.0109 | 0.0143 | 0.0226 | 0.0416 | 0.0812 |
| T(r = 10), hysteretic η = 2ζ | 0.0101 | 0.0102 | 0.0103 | 0.0109 | 0.0129 |

Going from ζ = 0.02 to 0.4 cuts the peak by a factor of 15 and multiplies the transmission at r = 10 by 7.5.

### 2.4 Hysteretic damping: why elastomers are used

Rubber does not behave like a dashpot. Over a wide band its energy loss per cycle is roughly independent of frequency, which is modelled by a complex stiffness k(1 + iη) with a constant **loss factor** η. Then (−mω² + k(1 + iη))X = k(1 + iη)Y and

**T_h(r) = √[ (1 + η²) / ((1 − r²)² + η²) ].**

At resonance T_h = √(1 + η²)/η (2.69 for η = 0.4, close to the viscous 2.73 at ζ = 0.2, since η ≈ 2ζ at resonance). At high frequency T_h ≈ √(1 + η²)/r², a −40 dB/decade roll-off. The table shows the consequence: with the same resonant peak, an elastomer isolates almost as well as an undamped spring at r = 10. Real elastomers are not this simple (their stiffness rises with frequency and falls with temperature, and η varies), but the lesson stands: *where* the damping force acts matters as much as *how much* there is.

### 2.5 The corner frequency for a required attenuation

Often the requirement is a transmissibility T* at some frequency f*: "no more than 0.06 at 500 Hz and above". Solve T(r, ζ) = T* for r > √2 (T falls monotonically there), and fₙ = f*/r is the highest natural frequency that meets it. Undamped, T = 1/(r² − 1) gives r = √(1 + 1/T*): for T* = 0.05 that is r = 4.58 and fₙ ≤ 109 Hz at 500 Hz. With ζ = 0.1 the code gives fₙ ≤ 90.0 Hz; with ζ = 0.2, fₙ ≤ 59.2 Hz. Damping costs corner frequency.

### 2.6 The random response, from P06.16

On a random input with one-sided PSD G(f), the transmitted acceleration and the sway are the integrals of `P06.16`:

**Grms_out² = ∫ T(f/fₙ)² G(f) df,  z_rms² = g₀² ∫ G(f) / [(ωₙ² − ω²)² + (2ζωₙω)²] df.**

The sway budget is the static sag plus three standard deviations of z: δ_static + 3z_rms. For an IMU, fₙ must also stay above the guidance and control bandwidth with margin: 60 Hz here, an illustrative value.

The design question is now a small optimisation: over fₙ and ζ, minimise Grms_out subject to fₙ ≥ 60 Hz, δ_static + 3z_rms ≤ 0.8 mm, T_peak ≤ 4 and T(500 Hz) ≤ 0.06.

### 2.7 Where the model fails

- **Six degrees of freedom.** Four mounts under a box give three translations and three rotations; if the mounts' elastic centre is not at the box's centre of mass, translation and rocking couple and there are two or more resonances in the band, and a translation input produces a rotation the gyros see.
- **Elastomer stiffness is not constant.** Dynamic stiffness above the static value by 20–50 % is common, and cold can double it; both raise fₙ and spoil the high-frequency isolation (part 4, experiment 4).
- **The shelf is not rigid,** and the IMU's resonance modifies the shelf's motion near fₙ (`P06.16`, 2.9).
- **Nonlinearity.** Mounts stiffen at large deflection and snubbers limit travel; a bottoming mount transmits shocks.
- **Cables.** The harness is a parallel spring and damper that the isolator analysis often forgets; a stiff cable bundle can short-circuit an isolator.

## 3. Build the computational model

Open `starter/main.py`. The `P06.16` helpers and `main` are given.

### Step 1: the isolator in closed form

`natural_frequency`, `stiffness_for_frequency`, `damping_coefficient` and `static_deflection` are 2.1. `transmissibility`, `transmissibility_phase` (use `np.angle` of the complex ratio, negated) and `transmissibility_hysteretic` are 2.2 and 2.4. `peak_transmissibility` is 2.3; `zeta_for_peak` and `fn_for_attenuation` invert T with `scipy.optimize.brentq` (bracket r between √2 and a large number for the second).

### Step 2: the independent check

`steady_state_ivp(r, zeta)` integrates the equation of 2.1 from rest with y = sin ωt using `solve_ivp` (DOP853, rtol 10⁻¹⁰), long enough for the transient to fall below 10⁻⁸ (time ln(10⁸)/(ζωₙ), at least 30 cycles), then fits x = a sin ωt + b cos ωt to the last five cycles by least squares. The amplitude is √(a² + b²) and the lag is atan2(−b, a). Fitting both quadratures is far more accurate than picking peaks.

### Step 3: the random response and the design search

`random_response` evaluates both integrals of 2.6 on `freq_grid`. `design_sweep` loops over a grid of fₙ and ζ and returns four arrays: transmitted Grms, static + 3σ sway, peak T, T at f_hf. `choose_isolator` masks the infeasible points and returns the feasible point with the lowest Grms, with its k and c; it raises `ValueError` when nothing is feasible.

Running `main()` prints

```
base input: NASA-STD-7001C Table 2 workmanship level, 6.784 Grms; IMU mass 0.35 kg
rigidly mounted IMU: 6.784 Grms (the input itself); on a stiff bracket with its mode at 3000 Hz, zeta = 0.05: 7.722 Grms

transmissibility: closed form vs solve_ivp steady state (zeta = 0.1)
     r     closed      ivp        |diff|     phase closed  phase ivp [deg]
  0.500   1.328229   1.328229   2.0e-12       1.884        1.884
  1.000   5.099020   5.099020   5.0e-09      78.690       78.690
  1.414   1.000000   1.000000   6.8e-10     148.414      148.414
  3.000   0.145366   0.145366   2.8e-10     144.747      144.747
  10.000   0.022582   0.022582   4.9e-11     115.408      115.408
...
requirements: fn >= 60.0 Hz, static (5.0 g) + 3 sigma sway <= 0.8 mm, peak T <= 4.0, T(500 Hz) <= 0.06
chosen isolator: fn = 65 Hz, zeta = 0.2, k = 58.4 kN/m total (14.59 kN/m per mount of 4), c = 57.18 N s/m
  transmitted 2.876 Grms (42 % of the input), sway 0.757 mm, peak T 2.73, T(500 Hz) 0.0555
  an elastomer with nearly the same resonant peak (eta = 2 zeta = 0.4) has T(500 Hz) = 0.0185
  integral 2.876 Grms vs Miles 2.880 Grms; sway: static 0.294 mm + 3 x 0.1544 mm

SYNTHETIC time-domain check (seed 17, 32 s): transmitted 2.881 Grms (+0.2 %), sway 0.1549 mm rms (+0.3 %); expected scatter about +/- 1.0 % (1 sigma)
```

and writes `vibration_isolation.png`: the transmissibility curves with the `solve_ivp` points, the sway–Grms trade for six dampings with the chosen design, and the input and isolated PSDs.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-17-vibration-isolation`. The 17 tests check, in order: the consistency of k, fₙ, c and static sag; the limits and the √2 crossover for ζ from 0 to 1; the phase at r → 0, r = 1 (arctan(1/(2ζ))) and r → ∞ (90° for viscous damping); the −1 and −2 log-log slopes of the viscous and hysteretic models between r = 100 and 1000; that more damping lowers the peak and raises T at r = 5 and 10; the exact peak against a 400 001-point search for four dampings, and T_peak → 1/(2ζ); `zeta_for_peak` and `fn_for_attenuation` inverting their functions, including the undamped closed form; the closed form against `solve_ivp` at (r, ζ) = (0.5, 0.1), (1, 0.1), (√2, 0.3), (3, 0.05), (8, 0.2); the random response tending to the input Grms for a rigid mount and to almost nothing for a 2 Hz mount, rising in Grms and falling in sway as fₙ rises; Miles within 1 % for light damping on the plateau and the white-noise sway formula of `P06.16`; the seeded time-domain check; the sweep's four trends; and `choose_isolator` against a brute-force search of the same grid, with an infeasible specification raising `ValueError`.

### Compare

The first block of output is the comparison with an established tool: `solve_ivp`'s DOP853 integrator, which knows nothing of the closed form, agrees with T to 5 × 10⁻⁹ or better and with the phase to the digits shown, from r = 0.5 to r = 10. The largest difference is at resonance, where the transient decays slowest. Miles' equation (2.880 Grms) agrees with the integral (2.876 Grms) to 0.15 % here: fₙ = 65 Hz sits on the +3 dB/oct slope, whose error (Miles high) is offset by the heavier damping's wider resonance.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Closed form against `solve_ivp`.** The table of part 3, extended to the test cases. Predict where the agreement is worst. *Finding:* worst at r = 1 (5 × 10⁻⁹), where the transient decays slowest relative to the period and the integration runs longest; 7 × 10⁻¹⁰ or less at the other four ratios.

2. **The damping trade.** The table of 2.3. Predict the factor by which ζ = 0.4 raises T(r = 10) over ζ = 0.02. *Finding:* 7.5 for viscous damping, 1.28 for hysteretic damping with the same peaks. The viscous damper is the problem, not damping as such.

3. **Natural frequency at ζ = 0.2.** Sweep fₙ from 40 to 200 Hz. Predict which requirements bind at the low and high ends. *Finding:*

   | fₙ [Hz] | 40 | 60 | 80 | 100 | 120 | 160 | 200 |
   |---|---|---|---|---|---|---|---|
   | transmitted Grms | 1.80 | 2.67 | 3.43 | 3.99 | 4.44 | 5.21 | 5.86 |
   | static + 3σ sway [mm] | 1.53 | 0.85 | 0.56 | 0.40 | 0.30 | 0.19 | 0.13 |
   | T(500 Hz) | 0.033 | 0.051 | 0.071 | 0.093 | 0.118 | 0.181 | 0.265 |

   Softer is quieter, but below 62.6 Hz the sway exceeds 0.8 mm, and above 69.6 Hz T(500 Hz) exceeds 0.06, *at this damping*: a feasible band only 7 Hz wide. Every isolator is caught between sway (too soft) and transmission (too stiff).

4. **The design search and its sensitivity.** Grid fₙ = 40–200 Hz in 5 Hz steps, ζ = 0.05–0.30 in steps of 0.025. Predict the binding requirements. *Finding:* fₙ = 65 Hz, ζ = 0.20 (2.88 Grms, 0.757 mm, peak 2.73, T(500 Hz) 0.0555). The sway and the 500 Hz requirement bind: at 60 Hz the sway is 0.849 mm, at ζ = 0.225 T(500 Hz) would exceed 0.06; the peak requirement (ζ ≥ 0.130) and fₙ ≥ 60 Hz are slack. A finer grid (0.5 Hz, ζ in steps of 0.005) moves the optimum to fₙ = 60.5 Hz, ζ = 0.235, with both binding constraints met almost exactly (0.794 mm, T(500 Hz) = 0.0595) and 2.53 Grms: the 5 Hz grid's answer is 14 % worse. Then perturb the built isolator (ζ = 0.2 fixed):

   | case | fₙ [Hz] | Grms | sway [mm] | T(500 Hz) |
   |---|---|---|---|---|
   | as designed | 65.0 | 2.88 | 0.757 | 0.056 |
   | IMU + harness 20 % heavier | 59.3 | 2.64 | 0.863 | 0.050 |
   | elastomer 30 % stiffer dynamically | 74.1 | 3.22 | 0.627 | 0.065 |
   | elastomer twice as stiff (cold) | 91.9 | 3.78 | 0.454 | 0.084 |
   | elastomer 30 % softer (hot) | 54.4 | 2.43 | 0.978 | 0.046 |

   Each plausible perturbation breaks one requirement: heavier or hotter breaks sway; stiffer or colder breaks the 500 Hz requirement. A design on the edge of two constraints has no margin for an elastomer whose stiffness is uncertain by ±30 %. Changing damping alone does not rescue it either: ζ = 0.15 gives 0.841 mm, ζ = 0.25 gives T(500 Hz) = 0.068.

5. **Hysteretic instead of viscous.** The same 65 Hz mount with η = 0.4 instead of ζ = 0.2. *Finding:* transmitted 2.81 Grms (against 2.88), peak 2.69, T(500 Hz) 0.0185 (against 0.0555). Only 1.8 % of the viscous design's transmitted mean square lies above 200 Hz, so the better roll-off barely changes the Grms; it changes the 500 Hz requirement by a factor of three, which is what frees the design.

Raw data format: CSV tables `r, T_closed, T_ivp, phase_closed_deg, phase_ivp_deg` (1); `zeta, T_peak, T10_viscous, T10_hysteretic` (2); `fn_Hz, grms, sway3_mm, T500` (3, 4); `case, fn_Hz, grms, sway3_mm, T500` (4); plus `vibration_isolation.png`.

## 5. Validate and improve

Two independent checks support the design.

| quantity | closed form / integral | independent route | difference | expected |
|---|---|---|---|---|
| T at r = 0.5, 1, √2, 3, 10 (ζ = 0.1) | 1.328229, 5.099020, 1, 0.145366, 0.022582 | `solve_ivp` steady state | ≤ 5 × 10⁻⁹ | integrator tolerance |
| phase at the same r | 1.884°, 78.690°, 148.414°, 144.747°, 115.408° | `solve_ivp` fit | < 10⁻⁶ rad (test) | — |
| transmitted Grms, chosen design | 2.876 | SYNTHETIC simulation, seed 17, 32 s: 2.881 | +0.2 % | ± 1.0 % (1σ) |
| sway RMS | 0.1544 mm | SYNTHETIC simulation: 0.1549 mm | +0.3 % | ± 1.0 % |
| transmitted Grms, 20 seeds × 8 s | 2.876 | mean ratio 0.9990, sd 2.28 % | −0.1 % | sd predicted 1.96 % |
| transmitted Grms | 2.876 | Miles 2.880 | +0.15 % | approximation |

The closed form and the integrator agree to their tolerances; the integral and the seeded simulation agree well inside one standard deviation; the observed seed-to-seed scatter (2.28 %) matches the 1/(2√(ζωₙT)) law (1.96 %) within the scatter of a 20-sample standard deviation.

The engineering answer: **mounts of 58.4 kN/m and 57.2 N·s/m in total (fₙ = 65 Hz, ζ = 0.20) cut the IMU's vibration from 6.78 to 2.88 Grms while keeping the sway to 0.76 mm, the resonant peak to 2.7 and the 500 Hz transmissibility to 0.056.** The design is bound by the sway and the high-frequency requirement together, and experiment 4 shows it has no margin against a ±30 % elastomer stiffness. Two improvements follow from the numbers. First, prefer a mount whose damping is hysteretic (an elastomer chosen for its loss factor rather than a viscous damper in parallel): it meets the 500 Hz requirement with a factor of three to spare, which leaves room to stiffen the mount and recover sway margin. Second, measure the mount's dynamic stiffness at temperature before freezing the design; the stiffness uncertainty, not the analysis, decides whether the design works.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of T, the √2 crossover and r_p; `vibration_isolation.png`; the experiment tables; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) with the chosen k and c, the binding requirements, the sensitivity to elastomer stiffness and the first test you would run.

### Adversarial self-checks

- **What if the input level doubles (+3 dB)?** Transmitted Grms and dynamic sway rise by √2 (2.88 → 4.07 Grms, 3z_rms 0.46 → 0.65 mm), and the sway budget is broken (0.29 + 0.65 = 0.95 mm). Linear isolators do not care about level; their sway budgets do.
- **What if the damping goes to zero?** The peak becomes infinite and T(500 Hz) improves; on a broadband input the transmitted Grms and the sway grow without bound as 1/√ζ, because the resonance sits on an energetic part of the spectrum.
- **Which assumption fails first?** Constant elastomer stiffness (experiment 4), then single-axis motion (rocking modes), then the forgotten cable harness.
- **What measurement could falsify the model?** A sine sweep of the isolated IMU on a shaker at low level: the measured T(f) should cross 1 at √2 fₙ, peak at about 2.7 and fall at −20 dB/decade if the damping is viscous or −40 if it is hysteretic. A second peak means a rocking mode; a peak that moves with level means a nonlinear mount.
- **Does the result survive a change of method?** The closed form, `solve_ivp`, the random-response integral, Miles and the time-domain simulation agree to their stated errors. The grid does not survive refinement: the 5 Hz, 0.025 grid gives 65 Hz and ζ = 0.20; a 0.5 Hz, 0.005 grid gives 60.5 Hz and ζ = 0.235, 12 % lower in Grms. Report the optimum with the grid it was found on, or refine until it stops moving.
- **Which decision changes if the uncertainty doubles?** With the elastomer stiffness known to ±15 %, viscous designs exist that meet all four requirements at both ends of the band (for example fₙ = 70 Hz, ζ = 0.18). At ±30 % none does, on a grid of fₙ = 40–120 Hz and ζ = 0.05–0.40. Doubling the uncertainty changes the decision from "choose k and c" to "change the damping mechanism or relax a requirement". Hysteretic mounts (with the sway estimated by ζ = η/2) do have designs robust to ±30 %, for example 72 Hz with η = 0.5.

## Check yourself

1. Above what frequency ratio does any linear spring–damper isolator start to isolate, and why does that not depend on damping?
2. What is the transmissibility at resonance for ζ = 0.1, approximately and exactly?
3. Why does a viscously damped isolator roll off at only −20 dB/decade, while a hysteretic one rolls off at −40?
4. An IMU mount sags 0.1 mm under the IMU's weight. What is its natural frequency?
5. You must halve the transmissibility at 500 Hz without changing the damping ratio. Which way do you move fₙ, and what happens to the sway?
6. What two contributions make up the sway budget in this lesson, and which dominates for the chosen design?
7. Why does more damping *reduce* the transmitted Grms on the workmanship spectrum here, even though it worsens isolation at high frequency?
8. Name two reasons the as-built isolator's natural frequency may differ from the design value.
