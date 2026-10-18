# P06.19  Home Lab: Pendulum System Identification 🔬

> Pluck a structure, record the decay, and fit a model: that is how every aircraft's modes are measured. Do it on a string pendulum, and find out how much your error bars can be trusted before you trust the numbers inside them.

| | |
|---|---|
| **Learning outcomes** | 1. Predict, before the swing and with uncertainties, a pendulum's small-angle frequency from L and g, its large-amplitude correction from the exact law, and an order-of-magnitude damping ratio from air drag on the bob. 2. Log a free decay with a phone video or a gyroscope, and extract the damped frequency (zero crossings and FFT peak) and the damping ratio (logarithmic decrement over n cycles and a fit of ln(peak) against cycle number). 3. Fit θ(t) = A e^{−ζωₙt} cos(ω_d t + φ) + c by nonlinear least squares with `curve_fit`, report parameter uncertainties from the covariance, and check them against a seeded bootstrap and the scatter of repeat runs; recognise from the residuals when the covariance understates the uncertainty. 4. Tell viscous, quadratic (air-drag) and Coulomb (pivot-friction) damping apart by the shape of the decay envelope, and measure the amplitude-dependent frequency against the exact law. 5. Show, on SYNTHETIC data with a known truth, that the analysis recovers the generator's parameters within its stated uncertainty. |
| **Prerequisites** | `P04.06` (the pendulum, the exact period, the `p04-pendulum` kit and its frame-timing uncertainty), `P06.13` (damping ratio, logarithmic decrement), `P02.18` (least squares, nonlinear fits, parameter covariance), `P02.17` (the bootstrap), `P02.10` (FFT, windowing, zero padding), `P00.13` (safety), `P00.14` (measurement and uncertainty), `P00.15` (the first experiment, normalised error). Skill check: you can write the logarithmic decrement of x(t) = A e^{−ζωₙt} cos ω_d t over n cycles, and you know what the diagonal of a covariance matrix means. |
| **Engineering question** | A string pendulum with L = 0.800 ± 0.003 m is released from 20° and logged for 120 s. What are its small-angle natural frequency and its damping ratio, with uncertainties you can defend; does the frequency change with amplitude as predicted; and is the damping viscous, quadratic or Coulomb? On the SYNTHETIC record of this lesson (seed 619) the answers are: f₀ = 0.55732 Hz against a predicted 0.5573 ± 0.0011 Hz; ζ = 9.21 × 10⁻⁴ ± 0.48 × 10⁻⁴ (block bootstrap; the covariance claims ± 0.10 × 10⁻⁴); the frequency rises as the swing decays, at 0.86 of the exact small-amplitude rate; and the envelope exponent p = 1.38 ± 0.29 lies between viscous and quadratic, which is what the generator put in. |
| **Assumptions and limits** | A simple pendulum (point bob, massless inextensible string, rigid pivot) moving in one plane; damping light (ζ ≪ 1) and made of viscous, quadratic and Coulomb parts; the damped-cosine fit assumes a constant frequency and damping ratio, which is only approximately true at 20° (the frequency changes by 0.5 % and ζ_eq by 20 % over the record); angle readings with independent noise. Fails first: a support that rocks, the bob's pouch or a logging cable adding a second mode or extra damping, a swing that is not planar (an ellipse), and amplitude errors in the frequency correction. All data shown in this lesson are SYNTHETIC; your own record replaces them. |
| **Analytical work** | The small-angle frequency and its uncertainty; the exact amplitude correction ω(A)/ω₀ = π/(2K(sin²(A/2))) and its series; the cycle-averaged energy balance for each damping law, giving dA/dt = −(b₁/2)A − (4b₂ω/3π)A² − 2b₀/(πω) and the equivalent damping ratio; the closed-form envelopes for each law; the logarithmic decrement and its uncertainty; the drag coefficient b₂ = ρC_dSL/(2m); the linear-interpolation zero crossing; the parabolic peak refinement; the covariance of a least-squares fit and what makes it wrong. |
| **Computational deliverable** | `starter/main.py`: `amplitude_frequency_ratio`, `predict_frequency`, `quadratic_drag_coefficient`, `envelope_rate`, `equivalent_zeta`, `pendulum_rhs`, `zero_crossings`, `period_from_crossings`, `fft_peak_frequency`, `find_peaks_interp`, `log_decrement`, `log_decrement_fit`, `damped_cosine`, `fit_damped_cosine`, `bootstrap_fit`, `envelope_power_law`, `fit_envelope`, `envelope_viscous_quadratic`, `fit_envelope_vq`, `frequency_vs_amplitude` (the SYNTHETIC tools `simulate_swing`, `make_synthetic_swing`, `write_synthetic_csv`, `read_swing_csv`, and `initial_guess` and `main` are given). Tests (20) prove: the exact correction against its series and a `solve_ivp` period; the prediction's propagation; the three damping signatures; the cycle-averaged decay against simulation for each law; energy conservation and decay; the generator's reproducibility, noise level and labelled CSV round trip; crossings, FFT and peaks on known signals; the decrement on an exact exponential; recovery of a known truth by NLS and log decrement within 3σ; the covariance against a seeded bootstrap and against the scatter of 12 repeat runs; the envelope exponent separating Coulomb, viscous and quadratic damping; recovery of b₁ and b₂ together; the frequency–amplitude law; and the committed SYNTHETIC CSV. |
| **Practical activity** | **physical experiment**: log a swing of the `p04-pendulum` string pendulum (string ≤ 2 m, bob ≤ 0.2 kg, release ≤ 45°) with a phone video, the phone's gyroscope or an IMU, after writing a dated prediction; fit frequency and damping with uncertainties; diagnose the damping. See `projects/p06-pendulum-sysid/`. A learner without the kit completes the lesson with the SYNTHETIC record `projects/p06-pendulum-sysid/data/SYNTHETIC_seed619_swing.csv`, labelled as synthetic. |
| **Equipment** | Minimum (Level 1): the `p04-pendulum` kit, a phone camera on a stand, a printed protractor or grid, a contrast marker. Recommended: the phone's gyroscope logged by a free sensor-logging app, with the phone as the bob in a padded pouch (total ≤ 0.2 kg), or slow-motion video; a kitchen scale (`p00-calibration`). Enhanced: an IMU breakout on the `p04-electronics-bench` microcontroller. See `projects/p06-pendulum-sysid/bom.yaml` (reuses `p04-pendulum`, `p00-calibration` and `p04-electronics-bench`). |
| **Safety** | See `projects/p06-pendulum-sysid/safety.md`. Summary box at the start of part 4. |
| **Evidence** | `data/prediction.md` (dated before the first swing), `data/calibration.csv`, `data/raw/*.csv`, the three plots (data with fit, decay envelope, frequency against amplitude), the results table with every uncertainty estimate, the damping diagnosis, `pendulum_sysid.png`, the report. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; all four competencies. Threshold: a prediction with uncertainty recorded before measurement; all public tests pass; f₀ and ζ reported with uncertainties checked by at least two methods and E_n for f₀ against the prediction; the damping mechanism argued from the envelope with its uncertainty; the safety box completed. |
| **References** | `@openstax-university-physics-vol1` 15.4 Pendulums (The Simple Pendulum); 15.5 Damped Oscillations (required-free). `@mit-ocw-16-07-dynamics` L24 "Pendulums" (Simple Pendulum; Small amplitude approximation) (required-free). `@mit-ocw-2-003sc-engineering-dynamics` Lecture 20: Linear System Modeling a Single Degree of Freedom Oscillator (the experimental determination of the damping ratio) (required-free). `@jcgm-100-2008` GUM: Type A and Type B evaluation, combined standard uncertainty, coverage factor (required-free). `@openstax-introductory-statistics-2e` 12.3 The Regression Equation (residuals) (required-free). `@scipy-userguide` "Optimization (scipy.optimize)" (required-free): `curve_fit`. Logging tools, one option each, never required (supplementary-free): `@phyphox-app` phyphox (the phone's gyroscope, CSV export, FAQ "Export formats"); `@tracker-video-analysis` Tracker video analysis and modeling tool (frame-by-frame tracking of a marker). Telling damping laws apart by the envelope, and bootstrap uncertainty for fit parameters, are original instruction. |
| **Extensions** | ★ Fit the full nonlinear ODE (`simulate_swing` with b₁, b₂, b₀ and L as parameters) directly to θ(t) with `scipy.optimize.least_squares`, which removes the constant-frequency assumption; compare the residual RMS with the damped-cosine fit's. ★ Change the bob's shape (a flat card edge-on, then face-on) at the same mass and measure b₂ against the drag prediction. ★ Hang the phone in its pouch from a single point and find the second (pouch-swinging) mode in the gyro spectrum: a double pendulum (`P06.11`). |
| **Time** | unmeasured |

## 1. The engineering problem

Before a new aircraft flies, it is shaken on the ground. A ground vibration test excites the structure (with shakers, or by a sharp tap) and records how it responds. From the records, engineers identify each mode's natural frequency and damping ratio, and those numbers go straight into the flutter analysis, where a damping ratio of 0.02 instead of 0.01 can move the predicted flutter speed. A spacecraft is shaken the same way before launch (`P06.16`). The identification step, fitting a model to a measured decay and saying how well its parameters are known, is the subject of this lesson.

A string pendulum is the cheapest structure with a mode you can test. It has one dominant mode whose frequency you can predict to a fraction of a per cent from its length, a damping ratio around 10⁻³ (as light as a deployed solar array's), and three damping mechanisms you can tell apart: air drag on the bob (quadratic in velocity), drag on the string and losses in the clip (roughly viscous), and friction at a rubbing pivot (Coulomb). Its frequency also depends slightly on amplitude, the nonlinearity of `P04.06`, which a fit that assumes a constant frequency will absorb into its errors. That is the realistic difficulty: the model you fit is never quite the system you measured.

The success criteria, written first: a dated prediction with uncertainties exists before the swing; f₀ agrees with it (E_n < 2) or a cause is found and tested; ζ falls within the predicted order of magnitude; every uncertainty is checked by a second method; on SYNTHETIC data the analysis recovers the known truth within its stated uncertainty; and the damping mechanism is argued from the envelope, with an honest statement of how decisive the evidence is.

## 2. Understand and derive

### 2.1 The model

A bob on a string of length L (to the bob's centre), at angle θ from the vertical, with three damping torques per unit inertia:

**θ̈ + (g/L) sin θ + b₁θ̇ + b₂θ̇|θ̇| + b₀ sgn θ̇ = 0.**

- **Viscous**, b₁ [s⁻¹]: a torque proportional to angular velocity, from drag on the thin string and losses in the clip at low Reynolds number. On its own, with sin θ ≈ θ, it is the damped oscillator of `P06.13` with ζ = b₁/(2ω₀).
- **Quadratic**, b₂ [rad⁻¹]: air drag on the bob, ½ρC_dS v² with S the bob's frontal area and v = Lθ̇. Its torque about the pivot, divided by the inertia mL², gives **b₂ = ρC_dSL/(2m)**. For a 50 mm sphere (C_d ≈ 0.47) of 150 g on 0.8 m in air (ρ = 1.2 kg/m³), b₂ = 0.0030 rad⁻¹.
- **Coulomb**, b₀ [rad/s²]: a constant friction torque at a pivot that rubs (a string over a rod, a loose hinge), always opposing the motion. The clip of the `p04-pendulum` kit makes it small.

The code integrates this equation with `solve_ivp` (the sign function smoothed as tanh(θ̇/ε) with ε = 10⁻³ rad/s so the integrator does not stall), and the SYNTHETIC generator samples it like a phone video.

### 2.2 Frequency: the prediction and the amplitude effect

The small-angle frequency is **f₀ = (1/2π)√(g/L)**. Its uncertainty follows from the power rule (`P00.14`): u_f/f = ½√((u_g/g)² + (u_L/L)²). With L = 0.800 ± 0.003 m and g = 9.81 ± 0.01 m/s², **f₀ = 0.5573 ± 0.0011 Hz**. The length dominates: 3 mm on 0.8 m is 0.19 % in f, the 0.01 in g only 0.05 %.

At amplitude A the exact period of the undamped pendulum is T = 4√(L/g)K(m) with m = sin²(A/2) (`P04.06`), so

**ω(A)/ω₀ = π/(2K(sin²(A/2))) ≈ 1 − A²/16.**

At 20° the frequency is 0.76 % low (0.5531 Hz). As a damped swing decays, its frequency therefore *rises*. Over the SYNTHETIC record the amplitude falls from 20.0° to 13.5°, so the instantaneous frequency rises by 0.4 %. A fit with one constant frequency returns an amplitude-averaged value; to compare with f₀ it must be divided by the mean of ω(A)/ω₀ over the record.

The damping also lowers the frequency, ω_d = ωₙ√(1 − ζ²), but at ζ = 10⁻³ by only 5 × 10⁻⁷: invisible beside the amplitude effect.

### 2.3 Damping: the energy balance over one cycle

When the damping is light, the swing over one cycle is very nearly θ = A cos ωt with A changing slowly. The energy per unit inertia is E = ½ω²A², so dE/dt = ω²A dA/dt, and over one cycle the damping removes:

- viscous: ∮ b₁θ̇² dt = b₁A²ω²(T/2) = πb₁A²ω;
- quadratic: ∮ b₂|θ̇|³ dt = b₂A³ω³ × (8/3)/ω = (8/3)b₂A³ω²;
- Coulomb: ∮ b₀|θ̇| dt = 4b₀A.

Dividing by the period 2π/ω and by ω²A gives the cycle-averaged envelope law:

**dA/dt = −(b₁/2)A − (4b₂ω/3π)A² − 2b₀/(πω),**

and, writing dA/dt = −ζ_eq ωA, the **equivalent damping ratio**

**ζ_eq(A) = b₁/(2ω) + (4b₂/3π)A + 2b₀/(πω²A).**

The three mechanisms have three signatures. Viscous damping gives a ζ independent of amplitude and an exponential decay. Quadratic drag gives ζ proportional to A: large swings die fast, small ones linger, and 1/A grows linearly in time. Coulomb friction gives ζ proportional to 1/A: the amplitude falls by the same 4b₀/ω² every cycle, linearly in time, and stops in finite time. For the predicted b₂ = 0.0030 rad⁻¹, ζ_eq = 4.4 × 10⁻⁴ at 20° and 1.1 × 10⁻⁴ at 5°.

A single law dA/dt = −κA^p covers the three pure cases with p = 0, 1, 2. Its solution is

A(t) = A₀[1 − (1 − p)κA₀^{p−1}t]^{1/(1−p)}  (and A₀e^{−κt} at p = 1),

and fitting it to the peaks gives p with an uncertainty: a **test that distinguishes the mechanisms**. For a mixture of viscous and quadratic damping, dA/dt = −αA − βA² has the closed form

A(t) = αA₀e^{−αt} / (α + βA₀(1 − e^{−αt})),  with α = b₁/2 and β = 4b₂ω/(3π),

which recovers b₁ and b₂ separately, if the amplitude range is wide enough to tell a straight line in ζ_eq(A) from a constant.

### 2.4 The logarithmic decrement and its uncertainty

For the damped cosine, successive peaks fall by e^{−ζωₙT_d} per cycle, so

**δ = (1/n) ln(A₀/Aₙ),  ζ = δ/√(4π² + δ²) ≈ δ/(2π).**

Over n cycles the decrement uses only two peaks. If each peak height carries a relative uncertainty s, the decrement's is √2 s/n, and ζ's is √2 s/(2πn). Fitting a straight line to ln A_k against k uses every peak, and its slope's standard error (from the residual scatter) is the better estimate. Neither knows about the amplitude dependence of ζ_eq: they report an average over the cycles used.

### 2.5 Nonlinear least squares and what its covariance means

Fit θ(t) = A e^{−ζωₙt} cos(ω_d t + φ) + c by minimising the sum of squared residuals over five parameters with `scipy.optimize.curve_fit` (`P02.18`). The offset c absorbs a tilted protractor or a gyro bias. Near the minimum, the parameter covariance is s²(JᵀJ)⁻¹, with J the Jacobian of the model and s² the residual variance; the square roots of its diagonal are the standard errors. That estimate is right when the residuals are independent, of constant variance, and *the model is correct*. Two checks:

- **Bootstrap** (`P02.17`): add resampled residuals to the fitted curve, refit, repeat (100 times), and take the standard deviation of each parameter. Resampling single residuals assumes they are independent; resampling **blocks** (2 s here, about one swing) keeps their correlation over a swing.
- **Repeat runs**: release the pendulum again, fit again, and compare the scatter of the estimates with the covariance.

If the residual RMS is much larger than the known measurement noise, the residuals are not noise: they are model error, smooth in time, and the covariance (and the single-residual bootstrap) will understate the uncertainty of anything the model error biases.

### 2.6 Extracting the frequency without a fit

- **Zero crossings**: interpolate linearly between the samples on either side of each upward crossing of the mean, then fit crossing time against crossing number; the slope is T_d and its standard error follows from the residuals. Sixty crossings timed to a few milliseconds give T_d to about 10⁻⁴ relative.
- **FFT peak** (`P02.10`): resample the irregular timestamps onto a uniform grid, apply a Hann window, zero-pad by 16 and refine the peak bin with a parabola through the log magnitudes of three bins.
- **Peaks**: for each positive lobe between crossings, fit a parabola through the largest sample and two neighbours on each side. This finds the peak to a fraction of a frame. Peak heights are measured from the record's mean, which is not exactly the equilibrium for a decaying record; a constant offset cancels in the decrement only to first order.

### 2.7 Where the model fails

In order of how likely they are at home: the pivot moves (a support that flexes adds a second mode and drains energy); the swing is not planar (released with a sideways push, the bob traces an ellipse and the video sees only one component); a phone pouch swings relative to the string (a double pendulum); a logging cable adds stiffness and damping; the angle read from a video is distorted by perspective if the camera is not square to the swing; air currents from heating or a fan. The damped-cosine model itself fails in a known way: constant frequency and constant ζ are approximations whose errors are visible in the residuals.

## 3. Build the computational model

Open `starter/main.py`. The SYNTHETIC tools, `initial_guess` and `main` are given; the generator runs once `pendulum_rhs` is written.

### Step 1: prediction

`amplitude_frequency_ratio` (with `scipy.special.ellipk`, whose argument is the parameter m = k²), `predict_frequency`, `quadratic_drag_coefficient`, `envelope_rate` and `equivalent_zeta` are 2.2 and 2.3.

### Step 2: the model

`pendulum_rhs` is the equation of 2.1 with the smoothed sign. `simulate_swing` (given) integrates it with DOP853 at rtol 10⁻¹⁰ and dense output.

### Step 3: frequency and peaks

`zero_crossings`, `period_from_crossings`, `fft_peak_frequency` and `find_peaks_interp` are 2.6. `log_decrement` and `log_decrement_fit` are 2.4.

### Step 4: the fit and its uncertainty

`damped_cosine` is the model of 2.5; `fit_damped_cosine` calls `curve_fit` from `initial_guess` (FFT frequency, decrement damping, first-peak amplitude and phase); `bootstrap_fit` resamples residuals singly or in blocks with `numpy.random.default_rng(seed)`.

### Step 5: diagnosis

`envelope_power_law`, `fit_envelope`, `envelope_viscous_quadratic`, `fit_envelope_vq` and `frequency_vs_amplitude` (zero-crossing periods in windows of six cycles, with the mean peak in each) are 2.3 and 2.2.

Running `main()` on the SYNTHETIC record prints

```
prediction (L = 0.800 +/- 0.003 m, g = 9.81 +/- 0.01 m/s^2): f0 = 0.5573 +/- 0.0011 Hz; at 20 deg 0.5531 Hz (-0.76 %)
  air drag on a 50 mm, 150 g bob (Cd 0.47): b2 = 0.0030 /rad, zeta_eq at 20 deg = 4.37e-04, at 5 deg = 1.09e-04

data: SYNTHETIC CSV SYNTHETIC_seed619_swing.csv: 3600 samples over 120.0 s, mean interval 33.33 ms (sd 2.87 ms)
damped frequency: zero crossings 0.55443 +/- 0.00004 Hz (66 crossings); FFT peak 0.55445 Hz
peaks: 66, first 20.00 deg, last 13.47 deg
log decrement over 50 cycles: delta = 0.00575, zeta = 9.145e-04; ln(peak) line over all peaks: zeta = 9.113e-04 +/- 1.3e-05
NLS fit: A = 19.820 +/- 0.041 deg, zeta = 9.2078e-04 +/- 9.6e-06, omega_n = 3.48333 +/- 0.00003 rad/s (f_n = 0.55439 Hz), c = 0.485 +/- 0.013 deg; residual RMS 0.800 deg
  bootstrap sd (100 resamples): zeta 9.0e-06 (blocks of 2 s: 4.8e-05), omega_n 3.4e-05 (blocks: 1.5e-04) rad/s; covariance: 9.6e-06, 3.3e-05
envelope dA/dt = -kappa A^p: p = 1.383 +/- 0.294 (0 Coulomb, 1 viscous, 2 quadratic)
frequency against amplitude (10 windows of 6 cycles): f = 0.55694 Hz x (1 -0.0540 A^2); exact small-amplitude law: 1 - 0.0625 A^2
envelope dA/dt = -alpha A - beta A^2: b1 = 3.90e-03 +/- 1.9e-03 1/s, b2 = 2.88e-03 +/- 2.2e-03 1/rad (correlation -1.00)
small-angle frequency: NLS f_n / mean amplitude ratio = 0.55732 Hz; windowed intercept 0.55694 +/- 0.00037 Hz

SYNTHETIC truth (generator): omega_0 = 3.50179 rad/s, f0 = 0.55733 Hz; zeta_eq from 1.021e-03 at the first peak to 8.275e-04 at the last (mean 9.178e-04)
```

(`python main.py . <path to SYNTHETIC_seed619_swing.csv>` reads the committed file; with no path, `main` regenerates the same record from the seed) and writes `pendulum_sysid.png`, every panel titled SYNTHETIC: the first 20 s with the fit, the decay envelope with the power-law and exponential fits, and the windowed frequency against amplitude with the exact law.

### Verify

`pytest curriculum/phase-06/lesson-19-home-lab-pendulum-system-identification` runs the 20 tests listed in the metadata. Three are statistical, and each tolerance is derived in a comment: the bootstrap against the covariance (60 resamples, ratio 0.67–1.5), the repeat-run scatter (12 runs, the 99 % interval of a χ-distributed ratio with 11 degrees of freedom), and recovery within 3σ.

### Compare

The exact correction agrees with the period of a `solve_ivp` integration at 30° to 10⁻⁹. On a known damped cosine with 2 ms jitter, the zero crossings give the frequency to 2 × 10⁻⁵, the FFT peak to 2 × 10⁻⁴, and the peaks match the envelope to 2 × 10⁻⁴. `curve_fit` is the established tool here; the comparison that matters is with the truth, in part 5.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the data formats are in `projects/p06-pendulum-sysid/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). The apparatus is the `p04-pendulum` kit. Read the safety box first, every session. If you have no kit, run the whole analysis on `projects/p06-pendulum-sysid/data/SYNTHETIC_seed619_swing.csv` (or on `make_synthetic_swing` with your own seed) and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p06-pendulum-sysid/safety.md`).** Hazards: a swinging mass (bob or phone) striking a person, a pet, a screen or a window; the phone or bob falling from a failed pouch or clip; the support falling; tripping over the string or a logging cable. Energy limits: **string length ≤ 2.0 m, bob mass ≤ 0.2 kg** (phone, pouch and clip included when the phone is the bob), **release amplitude ≤ 45°** (at the limits: 3.4 m/s and 1.15 J at the bottom, 3.1 N peak string tension). Controls: a soft bob, or the phone zipped into a padded pouch hung from two points; **the swing arc plus 0.5 m kept clear of people, pets, screens and glass**, its ends marked on the floor, a blanket below when the phone swings; the support tested with twice the bob's weight; the string pinched in a clip; the camera on a stand outside the arc; any cable light and slack; release from beside the swing plane, stop by catching the string from the side. Supervision: none for an adult; children only with an adult handling release and catch. Stop criteria: anyone enters the arc; the support moves or creaks; the string frays or slips; the bob or phone loosens; a cable tugs; the phone gets warm. Not permitted: heavier or longer pendulums, larger amplitudes or loops, swinging over people, beds or screens, batteries on the bob, and **any powered drive of the pendulum**.

The setup, seen from the camera:

![The pendulum rig as seen from the camera: a string pinched in a clip at the pivot on a rigid support hangs a bob with a contrast marker; L is the length from the pivot to the centre of the bob; a protractor or grid is taped behind the pivot in the swing plane; the swing arc plus 0.5 m is kept clear.](figures/pendulum-rig.svg)
*Figure 1. The pendulum as the camera sees it. L is measured to the centre of the bob, and the swing plane is square to the camera so that the angle reads true against the protractor.*

**Prediction first (in `data/prediction.md`, dated).** Measure L three times. Compute f₀ ± u with `predict_frequency`, the frequency at your release amplitude and at half of it, and whether the difference is detectable from about 60 crossings (between 20° and 10° it is 0.57 %, about 80 times the 7 × 10⁻⁵ relative uncertainty of the crossing frequency). Weigh the bob, measure its diameter, and compute b₂ and ζ_eq at the release amplitude; write ζ as a range of a factor of three. Predict the amplitude after 120 s. Freeze the file.

**Procedure in outline** (full version in `procedure.md`). Calibrate the time base against a stopwatch, L to the bob's centre, the angle-reading scatter on a still frame (video) or the bias at rest (gyro). Release from rest at 20° and record at least 120 s untouched; three such runs; then one run at 40° and one at 10°, which separate the damping mechanisms and show the amplitude effect. Save every record raw.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_run<k>.csv`, one file per run, never edited). The two example rows below are SYNTHETIC: they are the first rows of `SYNTHETIC_seed619_swing.csv`, not a measurement:

```
t_s,theta_deg
0.000,20.362
0.030,20.145
```

or `t_s,omega_deg_s` for a gyroscope. The angular rate of a damped cosine is itself a damped cosine with the same ζ and ωₙ, so the same fit applies, with c the gyro bias.

**What the logging tier buys** (SYNTHETIC records of 120 s at 10°, NLS covariance uncertainties): video at 30 frames/s and 0.3° gives u_ζ = 7.0 × 10⁻⁶ and u_ωₙ = 2.4 × 10⁻⁵ rad/s; 60 frames/s gives 4.9 × 10⁻⁶ and 1.7 × 10⁻⁵ (√2 better, more samples); a gyroscope-like record at 100 Hz with 0.05° equivalent noise gives 1.2 × 10⁻⁶ and 4.1 × 10⁻⁶; a careless 1.0° reading at 30 frames/s gives 2.3 × 10⁻⁵ and 7.9 × 10⁻⁵. These are noise-only figures: part 5 shows that model error, not noise, sets the real uncertainty at 20°.

## 5. Validate and improve

**The independent check: recover a known truth.** The SYNTHETIC record was generated with L = 0.80 m, g = 9.81 m/s², b₁ = 0.003 s⁻¹ and b₂ = 0.004 rad⁻¹, so every quantity the analysis estimates has a known true value. (The analysis was not told them.)

| quantity (SYNTHETIC, seed 619) | estimate | its uncertainty | truth | difference |
|---|---|---|---|---|
| damped frequency, zero crossings | 0.55443 Hz | ± 0.00004 | 0.55440 Hz (exact law averaged over the peaks) | 0.75σ |
| damped frequency, FFT peak | 0.55445 Hz | — | 0.55440 Hz | 9 × 10⁻⁵ relative |
| ζ, log decrement over 50 cycles | 9.145 × 10⁻⁴ | ± 0.58 × 10⁻⁴ (two-peak propagation) | 9.38 × 10⁻⁴ (mean ζ_eq over those cycles) | 0.4σ |
| ζ, ln(peak) line over 66 peaks | 9.113 × 10⁻⁴ | ± 0.13 × 10⁻⁴ | 9.178 × 10⁻⁴ (mean ζ_eq) | 0.5σ |
| ζ, NLS | 9.208 × 10⁻⁴ | ± 0.096 × 10⁻⁴ (covariance); ± 0.48 × 10⁻⁴ (block bootstrap) | 9.178 × 10⁻⁴ | 0.3σ (cov.) |
| ωₙ, NLS (amplitude-averaged) | 3.48333 rad/s | ± 0.00003 (cov.); ± 0.00015 (block) | 3.48338 rad/s | 1.6σ (cov.), 0.4σ (block) |
| f₀, NLS ÷ mean ω(A)/ω₀ | 0.55732 Hz | ± 0.00002 (block) | 0.55733 Hz | 0.4σ |
| f₀, windowed intercept against A² | 0.55694 Hz | ± 0.00037 | 0.55733 Hz | 1.05σ |
| offset c | 0.485° | ± 0.013° | 0.5° | 1.2σ |
| envelope exponent p | 1.38 | ± 0.29 | between 1 and 2 (a mixture) | consistent |
| b₁, b₂ (viscous + quadratic envelope) | 3.9 × 10⁻³ s⁻¹, 2.9 × 10⁻³ rad⁻¹ | ± 1.9 × 10⁻³, ± 2.2 × 10⁻³ (correlation −1.00) | 3.0 × 10⁻³, 4.0 × 10⁻³ | 0.5σ, 0.5σ |

Every estimate recovers its truth within its stated uncertainty, but the table also shows which uncertainty to state. The fit's residual RMS is **0.80°**, against a known noise of 0.30°: the residuals are model error, because the damped cosine has one frequency and one ζ while the true swing's frequency drifts by 0.4 % and its ζ_eq by 20 % over the record. The covariance and the single-residual bootstrap (9.0 × 10⁻⁶ for ζ) treat that model error as white noise and understate the uncertainty; the block bootstrap, which keeps the error's correlation over a swing, gives **ζ ± 4.8 × 10⁻⁵ and ωₙ ± 1.5 × 10⁻⁴ rad/s**, five times larger. Repeat runs do not reveal the problem: six SYNTHETIC repeats with fresh noise scatter by only 5.0 × 10⁻⁶ in ζ and 1.5 × 10⁻⁵ in ωₙ, because every repeat has the same model error. Real repeats, released by hand at slightly different amplitudes, would scatter more; but a bias common to every run is invisible to repeat-run scatter. On a small-amplitude record where the damped cosine *is* the right model (the tests: 4°, purely viscous), covariance, bootstrap and repeat runs agree, and the residual RMS equals the noise.

**Prediction against measurement.** The SYNTHETIC f₀ = 0.55732 Hz against the prediction 0.5573 ± 0.0011 Hz: E_n ≈ 0. The measurement is fifty times more precise than the prediction, whose uncertainty is all in L, so the honest conclusion runs the other way: the swing measures the effective length, L = g/(2πf₀)² = 0.8000 m. The SYNTHETIC ζ ≈ 9.2 × 10⁻⁴ is twice the air-drag prediction for a 150 g, 50 mm bob (4.4 × 10⁻⁴): the generator's bob drags more (b₂ = 0.004 against 0.003 rad⁻¹) and it adds viscous loss. With your own pendulum, a measured ζ far above the drag prediction points at the pivot or the support.

**The damping mechanism.** The exponent p = 1.38 ± 0.29 rules out Coulomb friction (p = 0, at 4.7σ) and sits between viscous and quadratic, but a 20° → 13.5° decay cannot say more: the viscous-plus-quadratic fit recovers b₁ and b₂ only with a −1.00 correlation, because over a 1.5:1 amplitude range ζ_eq(A) = b₁/(2ω) + (4b₂/3π)A cannot tell a slope from an offset. The improvement is a wider amplitude range. On SYNTHETIC records released at 30° with lighter noise, the same test returns p = 0.06 ± 0.05 for pure Coulomb, 1.02 ± 0.05 for pure viscous and 1.97 ± 0.02 for pure quadratic damping; from 40° decaying to 6° over 150 s, the viscous-plus-quadratic fit recovers b₁ = 0.00989 ± 0.00020 s⁻¹ (truth 0.01) and b₂ = 0.02003 ± 0.00022 rad⁻¹ (truth 0.02). That is why the procedure asks for a 40° run and a long record.

**The amplitude effect.** The windowed frequency rises as the swing decays, with a fitted coefficient −0.054 on A² against the exact −0.0625: the effect is detected, at 0.86 of the predicted rate, within the scatter of ten windows (the intercept's uncertainty is 0.07 %).

| pattern in your data | likely cause | test |
|---|---|---|
| residual RMS ≫ reading noise | constant-frequency model at large amplitude | refit a 5° run; fit in windows; block bootstrap |
| f₀ low, E_n > 2, ζ high | support flexing or rocking | clamp the support stiffer; repeat |
| f₀ off, ζ normal | L measured to the wrong point | re-measure to the bob's centre |
| amplitude falls linearly, p ≈ 0 | string rubbing on a rod, loose clip | pinch the string in a clip |
| beating in the envelope, second spectral peak | elliptical swing, or a pouch swinging on its string | release more carefully; hang the pouch from two points |
| ζ different with and without the logging cable | the cable's stiffness and damping | run on video without the cable |

## 6. Defend your work

Submit: `data/prediction.md` (dated before the first swing), `data/calibration.csv`, `data/raw/*.csv`, your analysis script, the three plots, the results table with every uncertainty estimate and the one you report, the damping diagnosis, `pendulum_sysid.png`, and a report following `projects/p06-pendulum-sysid/report-template.md`. If you used the SYNTHETIC dataset, every table and figure says so in its caption.

### Adversarial self-checks

- **What if the damping doubles?** The decay over 120 s doubles in log terms, the amplitude range widens, and the envelope exponent becomes easier to determine; ζ's relative uncertainty from noise halves. The frequency is unaffected (ω_d/ωₙ changes by 10⁻⁶).
- **What if the amplitude goes to zero?** The amplitude effect and the quadratic drag vanish (ζ_eq → b₁/(2ω) plus any Coulomb part, which grows as 1/A and stops the swing); the damped cosine becomes the right model and the covariance becomes trustworthy, but the signal-to-noise ratio falls. A 5° run is the cleanest check of f₀.
- **Which assumption fails first?** The constant-frequency, constant-ζ fit at 20°: the residual RMS (0.80° against 0.30° of noise) says so directly. Then the rigid support, then a planar swing.
- **What measurement could falsify the model?** A swing whose amplitude falls linearly (p ≈ 0) falsifies viscous and drag damping and points at pivot friction; a frequency that falls as the amplitude decays (the opposite of the exact law) would mean the support is softening or the string stretching; a second spectral peak falsifies the single-mode model.
- **Does the result survive a change of method?** Zero crossings, the FFT peak and NLS agree on the frequency to 10⁻⁴; decrement, ln-peak line and NLS agree on ζ within about 1 %. The *uncertainty* does not survive a change of method: covariance and block bootstrap differ by a factor of five, and the report must say which it uses and why.
- **Which decision changes if the uncertainty doubles?** None for f₀, whose uncertainty is a hundredth of the prediction's. For the damping mechanism, a doubled u_p (0.59) would exclude Coulomb friction only at 2.3σ and could not separate viscous from quadratic damping at all; a pivot or bob redesign could not be justified from this run, and the 40° run becomes essential.

## Check yourself

1. What is the logarithmic decrement over n cycles, and how is it converted to a damping ratio?
2. Why does a free swing's frequency rise as it decays? By how much is a 20° swing slower than a small one?
3. How does the equivalent damping ratio depend on amplitude for viscous, quadratic and Coulomb damping?
4. A decay envelope is a straight line on linear axes. Which damping mechanism is that, and where on the apparatus would you look for it?
5. The covariance of a fit gives ζ ± 1 × 10⁻⁵ and the residual RMS is 2.7 times the known noise. What do you conclude, and what do you do?
6. Why can repeat runs with fresh noise fail to reveal a bias common to all of them?
7. Your f₀ measurement is fifty times more precise than your prediction. What does the comparison then measure?
8. Why must the prediction be written and dated before the first swing?
