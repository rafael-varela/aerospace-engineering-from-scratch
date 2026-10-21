# P04.06  Oscillations & Waves 🔬

> A restoring force proportional to displacement gives a frequency that does not depend on amplitude, and that one fact explains clocks, resonance, and why a solar array can shake itself apart. Every real oscillator departs from it; the skill is knowing by how much, and measuring it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive simple harmonic motion from a linear restoring force, its energy, and the free response of a damped oscillator in all three regimes; extract the damping ratio from a logarithmic decrement. 2. Derive the forced steady-state amplitude and phase, the quality factor and the half-power bandwidth, and use them to decide whether a structural mode will be excited by a disturbance frequency. 3. Derive the pendulum's small-angle period, state the amplitude correction (exact elliptic-integral result and its series), and predict, before measuring, whether your timing can detect it. 4. Derive the wave equation for a string, its speed √(T/μ), d'Alembert's solution and the standing-wave modes; implement and test a finite-difference wave solver with its stability (CFL) limit and energy conservation. 5. Measure the period of a string pendulum against length and amplitude, reduce the timings to g with an uncertainty, and compare with the prediction written beforehand. |
| **Prerequisites** | `P02.02` (second-order linear ODEs: the three damping regimes), `P02.03` (forced response and resonance), `P02.11` (the wave PDE and separation of variables), `P03.11` (finite differences for PDEs), `P04.02` (energy), `P00.13` (safety), `P00.14` (measurement and uncertainty), `P00.15` (frame timing, normalized error). Skill check: you can solve x″ + 2ζω₀x′ + ω₀²x = 0 for ζ < 1 and propagate an uncertainty through a power law. |
| **Engineering question** | Does a string pendulum's period follow 2π√(L/g) within the uncertainty of a phone-video measurement, at what amplitude does the departure become visible, and, as the same physics at a larger scale, which disturbance frequencies must a lightly damped spacecraft appendage be kept away from? |
| **Assumptions and limits** | Linear restoring force (small amplitude); point mass on a massless, inextensible string; viscous damping proportional to velocity; ideal string with uniform tension and small slopes. The first to fail on the bench is small amplitude (1.7 % at 30°); then the length reference (L must run to the bob's centre); the point-mass assumption costs only (1/5)(r/L)² of the period for a bob of radius r, 0.003 % for a 2 cm bob on 0.8 m; the viscous damping law is wrong in form (air drag on a bob is closer to quadratic) but barely touches the period. For waves: the finite-difference scheme is exact only at Courant number 1 and disperses below it. |
| **Analytical work** | SHM and its energy; the damped free response and log decrement; the forced steady state, Q and bandwidth; the pendulum equation, small-angle period and the amplitude correction; the string wave equation, wave speed, d'Alembert's solution with fixed ends, standing waves; the uncertainty of g from L and T. |
| **Computational deliverable** | `starter/main.py`: `shm`, `damped_response`, `log_decrement_to_zeta`, the pendulum periods (small-angle, exact, series, simulated with `solve_ivp`), the resonance curve, Q, bandwidth and peak frequency, `wave_equation_fd`, `dalembert`, `standing_wave`, `wave_energy`, and the data reduction (`period_from_timings`, `g_from_pendulum`, `fit_power_law`), plus `make_synthetic_pendulum_data` (SYNTHETIC, for tests only). Tests (19) check the equation of motion and energy of SHM, the damped closed form against `solve_ivp` in all regimes, the zero-damping limit, the small-angle limit and √L scaling, exact = series = simulated periods, the resonance curve's limits, peak and bandwidth, the phase, the wave speed and harmonics, the standing wave against the PDE, exactness and energy conservation of the FD solver at Courant 1, its second-order convergence below 1 and blow-up above, the reduction formulas, and recovery of g from synthetic data. |
| **Practical activity** | **physical experiment**: time a string pendulum (string ≤ 2 m, bob ≤ 0.2 kg) at five or more lengths and at several amplitudes with phone video, after writing the prediction down; reduce to g and to the amplitude dependence. See `projects/p04-pendulum/`. A learner without the kit completes the lesson with the SYNTHETIC dataset from `make_synthetic_pendulum_data`, labelled as synthetic. |
| **Equipment** | Level 1: string or thread, a compact bob of 50–200 g, a tape measure, a phone camera, a rigid support (door frame clamp or a stick weighed down on a table), a protractor printout. Optional: a phone accelerometer app, a clamp stand. See `projects/p04-pendulum/bom.yaml` (reuses `p00-first-experiment` items). |
| **Safety** | See `projects/p04-pendulum/safety.md`. Summary box at the start of part 4. |
| **Evidence** | `data/prediction.md` (dated before any timing), `data/raw/*.csv`, the T-against-L plot with k = 2 error bars and the 2π√(L/g) curve, the T/T₀-against-amplitude plot with the exact curve, g with its uncertainty, the E_n table, `oscillations_waves.png`, the report. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; all four competencies. Threshold: a prediction with uncertainty recorded before measurement; all public tests pass; g reported with an uncertainty budget and E_n against the local value; the amplitude effect either detected and compared with the exact curve or shown to be below the measurement's resolution, with the numbers; the safety box completed. |
| **References** | `@openstax-university-physics-vol1` Ch. 15 "Oscillations": 15.1 "Simple Harmonic Motion", 15.2 "Energy in Simple Harmonic Motion", 15.4 "Pendulums", 15.5 "Damped Oscillations", 15.6 "Forced Oscillations" (required-free; https://openstax.org/books/university-physics-volume-1/pages/15-6-forced-oscillations). Ch. 16 "Waves": 16.1 "Traveling Waves", 16.2 "Mathematics of Waves", 16.4 "Energy and Power of a Wave", 16.5 "Interference of Waves", 16.6 "Standing Waves and Resonance" (required-free). |
| **Extensions** | ★ Measure the damping: film 100 swings and fit the decay of the amplitude; decide whether it is exponential (viscous) or linear (dry friction at the pivot) or something in between (quadratic air drag). ★ Add a driven pendulum (move the support by hand at a steady rhythm) and find the resonance by ear. ★ Pluck a stretched string of known tension and mass per length, record the sound with a phone, and find the harmonics with an FFT (`P02.10`). |
| **Time** | unmeasured |

## 1. The engineering problem

Every structure that carries load vibrates. A wing bends and twists at a few hertz; a satellite's solar array, a long, light panel on a hinge, can have a first bending mode below 1 Hz; a launch vehicle's fuel slosh, an antenna boom, a control cable, all have natural frequencies. Each is a lightly damped oscillator, and a lightly damped oscillator driven near its natural frequency responds tens of times more strongly than it would to the same force applied slowly. That is resonance. A reaction wheel spinning at 48 revolutions per minute with a few grams of imbalance shakes the spacecraft at 0.8 Hz, and if the solar array's first mode is there, the pointing jitter can ruin every image the spacecraft takes. The engineering response is to know the natural frequency and the damping, and keep disturbances away from the resonance band. Flutter, where the airflow itself supplies negative damping, is the same mathematics with the sign of ζ reversed, and is the subject of later phases.

You will learn the physics on the simplest oscillator you can build at home: a weight on a string. The prediction is the oldest in the subject, T = 2π√(L/g), and it is exact only in a limit (small swings, point mass, rigid support, no air). This lesson derives it, derives how it fails as the swing grows, and asks you to predict, before timing anything, whether your phone camera can see that failure. Then it moves from one oscillator to a continuum of them, a stretched string, where the same restoring force produces travelling waves, reflections and standing-wave modes, and builds a finite-difference solver that propagates a pulse exactly.

Success criteria. Before timing: a dated table of predicted periods with uncertainties at each length and amplitude. After: g from each length with its uncertainty, E_n < 2 against your local g (or a named and tested cause where not), and a statement, with numbers, of whether the amplitude effect at 20–30° was detected. On the computer: the FD solver reproduces d'Alembert's solution to 10⁻¹² at Courant number 1 and conserves its discrete energy to 10⁻¹².

## 2. Understand and derive

### 2.1 Simple harmonic motion

A mass m on a spring of stiffness k, displaced by x, feels F = −kx (Hooke's law: a linear restoring force). Newton's second law gives mẍ = −kx, or

**ẍ + ω₀² x = 0**,  ω₀ = √(k/m).

From `P02.02`, the solution is **x(t) = x₀ cos ω₀t + (v₀/ω₀) sin ω₀t**, with period T = 2π/ω₀, *independent of the amplitude*. That independence is the defining property of SHM and the reason pendulum clocks work. The energy E = ½mv² + ½kx² is constant: differentiate it and substitute the equation of motion, dE/dt = v(mẍ + kx) = 0. The energy sloshes between kinetic and potential twice per period, and its total is ½kA² for amplitude A.

### 2.2 Damping

Add a viscous damper, a force −cẋ. Writing c/m = 2ζω₀ defines the **damping ratio** ζ:

**ẍ + 2ζω₀ẋ + ω₀²x = 0**.

The characteristic roots are −ζω₀ ± ω₀√(ζ² − 1). For ζ < 1 (underdamped) the motion oscillates at the damped frequency ω_d = ω₀√(1 − ζ²) inside an envelope e^{−ζω₀t}; for ζ = 1 (critical) it returns fastest without overshoot; for ζ > 1 (overdamped) it creeps back. The closed form for all three, with initial displacement x₀ and velocity v₀, is in `damped_response`, and a test checks each regime against `solve_ivp` to 10⁻⁸.

Two successive peaks of an underdamped response are one damped period apart, so their ratio is e^{ζω₀T_d}. Its logarithm, the **logarithmic decrement** δ = ln(x_n/x_{n+1}) = 2πζ/√(1 − ζ²), inverts to

**ζ = δ / √(4π² + δ²)**.

This is how damping is measured in practice: pluck the structure and read the decay. For a 1.5 Hz oscillator with ζ = 0.05, the code's closed form agrees with DOP853 to 2.1 × 10⁻¹³ m, and the decrement of its simulated peaks, 0.3146, returns ζ = 0.0500.

Energy is no longer conserved: dE/dt = −cẋ² ≤ 0. Light structures are lightly damped. A bolted aluminium structure has ζ of about 0.01–0.05; a deployed solar array in vacuum can be below 0.01, because there is no air to help.

### 2.3 Forced oscillation and resonance

Drive the oscillator with a force F₀cos ωt (per unit mass f₀ = F₀/m). After the transient decays, the response is x = X cos(ωt − φ) with

**X = f₀ / √((ω₀² − ω²)² + (2ζω₀ω)²)**,  **φ = atan2(2ζω₀ω, ω₀² − ω²)**.

Read the limits. At ω → 0, X = f₀/ω₀², the static deflection: a slow force just stretches the spring. At ω = ω₀, X = f₀/(2ζω₀²): the static deflection multiplied by the **quality factor Q = 1/(2ζ)**. At ω ≫ ω₀, X ≈ f₀/ω²: the mass is too sluggish to follow, and the response falls with the square of frequency. The phase lag runs from 0 through exactly 90° at ω₀ to 180°. The peak of X sits at ω₀√(1 − 2ζ²), slightly below ω₀, and vanishes for ζ ≥ 1/√2. The **half-power bandwidth**, the width of the band where X exceeds its peak divided by √2, is 2ζω₀ = ω₀/Q for small ζ. For the 1.5 Hz, ζ = 0.05 oscillator: Q = 10.0, amplitude at resonance 10.0 times static, bandwidth 0.150 Hz.

The engineering reading: with Q = 50 (ζ = 0.01), a disturbance anywhere inside a band of about ±5 % around the natural frequency is amplified by more than ten. Problem 3 does this for a solar-array boom and a reaction wheel.

### 2.4 The pendulum

A bob of mass m on a string of length L, at angle θ from the vertical. The tension acts along the string and does no torque about the pivot; gravity does, −mgL sin θ. With moment of inertia mL² about the pivot:

**θ̈ + (g/L) sin θ = 0**.

This is *not* SHM: the restoring torque is proportional to sin θ, not θ. For small angles sin θ ≈ θ, and the equation becomes SHM with ω₀ = √(g/L):

**T₀ = 2π √(L/g)**,

independent of mass and (in this limit) amplitude. For L = 1.000 m and g = 9.80665 m/s², T₀ = 2.0064 s.

At finite amplitude the period lengthens, because sin θ < θ and the restoring torque is weaker than the linear model assumes. Energy conservation gives the exact period as an integral (Problem 2), which is a complete elliptic integral of the first kind:

**T = 4√(L/g) K(k²)**,  k = sin(θ₀/2),

and its expansion for small amplitude is

**T ≈ T₀ (1 + θ₀²/16 + 11θ₀⁴/3072)**.

The relative correction at a few amplitudes, from the code's exact, series and `solve_ivp` periods (the three agree to the digits shown up to 30°):

| amplitude | 5° | 10° | 20° | 30° | 45° | 60° | 90° |
|---|---|---|---|---|---|---|---|
| T/T₀ − 1 | 0.048 % | 0.19 % | 0.77 % | 1.74 % | 4.00 % | 7.32 % | 18.0 % |

So "small angle" means: below about 10° the correction is under 0.2 %; at 30° it is 1.7 %. The series to θ₀⁴ is good to 10⁻⁵ at 30° but 0.4 % low at 90°.

### 2.5 Can your camera see the amplitude effect?

From `P00.15`, timing one event pair from video at frame rate f carries a standard uncertainty u = √2 × 0.5/(√3 f) = 0.41/f: 13.6 ms at 30 fps. Time **N swings** at once (from one passage through the lowest point in a given direction to the N-th later one) and divide: the frame term per period falls to 13.6 ms/N. For L = 1 m and N = 10, that is 1.36 ms on 2.006 s, 0.07 %. The 30° correction is 1.74 %, about 25 times larger: easily visible. The 10° correction, 0.19 %, is about 3 times the frame term: visible with repeats. The 5° correction is not. **Write this down before you measure.** It is the same question as `P00.15` asked about drag: not "is it negligible?" but "can my instrument see it?"

Counting many swings has a cost: the amplitude decays during them, so the period you measure is an average over a slightly decreasing amplitude. Record the amplitude at the start and the end of each run and use the mean in the correction.

### 2.6 From timings to g, with uncertainty

With T measured and L measured from the pivot to the **centre of the bob**,

**g = 4π² L / T²**,  **u_g/g = √((u_L/L)² + (2u_T/T)²)**

(the power rule from `P00.14`; the period's uncertainty counts twice because T is squared). For L = 1.000 ± 0.002 m and T = 2.006 ± 0.002 s, the length term is 0.2 % and the period term 0.2 %, so g is good to about 0.28 %, ±0.03 m/s². Two systematic corrections must be applied before comparing with a reference: the amplitude correction (divide T by T/T₀ from 2.4, so g rises by (T/T₀)²; at 0.15 rad this is a factor 1.0028) and the length reference (the centre of the bob, not its top: a 2 cm bob measured to the knot gives L about 1 cm short, a 1 % error in g).

A second, independent reduction uses all lengths at once: fit T = A L^p in log–log coordinates (`fit_power_law`). The model predicts p = 1/2 and A = 2π/√g. A fitted p far from 0.5 is evidence against the model (or against the length reference: a constant offset in L bends the log–log line).

### 2.7 Waves on a string

Take a string under tension T_s with mass per unit length μ, displaced by a small amount u(x, t). A short element dx has its two ends pulled by the tension along the local slopes; the net vertical force is T_s(u_x(x + dx) − u_x(x)) ≈ T_s u_xx dx (small slopes: assumption). Its mass is μ dx. Newton's second law:

**u_tt = c² u_xx**,  **c = √(T_s/μ)**,

the **wave equation** (`P02.11`). Any function of x − ct (a shape moving right at speed c) or of x + ct (moving left) satisfies it. For a string released from rest with shape f(x), **d'Alembert's solution** is u = ½[f(x − ct) + f(x + ct)]: the initial shape splits into two halves of half the height, moving apart at c. With fixed ends at x = 0 and L, extend f as an odd function about each end (period 2L), and the same formula holds; the odd extension is what makes a pulse come back *inverted* after reflecting from a fixed end.

Separation of variables (`P02.11`) gives the **standing waves**: u_n = sin(nπx/L) cos(nπct/L), with frequencies

**f_n = n c / (2L)**,  n = 1, 2, 3, …

the harmonic series. Each standing wave is two equal travelling waves going opposite ways (OpenStax 16.6). A guitar-like string with T_s = 50 N, μ = 2.0 g/m and L = 0.65 m has c = 158.1 m/s and modes 121.6, 243.3, 364.9, 486.5 Hz. A wave carries energy: per unit length, ½μu_t² kinetic plus ½T_s u_x² potential (OpenStax 16.4); for a string with fixed ends the total is conserved.

### 2.8 The finite-difference wave solver

Discretise x with spacing Δx and t with step Δt, and replace both second derivatives by central differences (`P03.11`):

**u_i^{n+1} = 2u_i^n − u_i^{n−1} + C² (u_{i+1}^n − 2u_i^n + u_{i−1}^n)**,  **C = cΔt/Δx**,

the **Courant number**. The first step, from rest, uses u_t = 0 to write u^{−1} = u^1, giving u^1 = u^0 + (C²/2)(Δ²u^0). The scheme is second order in both Δx and Δt and has three properties the tests check.

- **Stability (CFL condition): C ≤ 1.** The numerical domain of dependence must contain the physical one: information travels one cell per step in the scheme and c Δt in reality. For C > 1 the highest-frequency grid mode grows each step: on the 200-cell string of 2.7, a unit pulse at C = 1.02 has grown to 3.5 × 10⁴⁷ after one round trip.
- **Exactness at C = 1.** Substitute C = 1: u_i^{n+1} = u_{i+1}^n + u_{i−1}^n − u_i^{n−1}, which is exactly satisfied by any f(x − ct) + g(x + ct) sampled on the grid. The scheme then reproduces d'Alembert to rounding: 9.8 × 10⁻¹⁵ after a full round trip of a Gaussian pulse on 200 cells.
- **Energy.** A discrete energy (kinetic from the difference between levels, potential from the product of slopes at two levels) is conserved exactly by the scheme at every C ≤ 1: spread 9.5 × 10⁻¹⁵ at C = 1 and 3.7 × 10⁻¹⁵ at C = 0.8.

Below C = 1 the scheme is still stable and energy-conserving but **dispersive**: short wavelengths travel slower than c, and a sharp pulse sheds a train of ripples behind it. At C = 0.8 the maximum error after one round trip is 2.8 × 10⁻² of the pulse height. On a smooth standing wave (mode 3, C = 0.9, three periods) the error is 3.1 × 10⁻⁴: a phase error that converges at second order as the grid is refined.

### 2.9 Where the models fail

- **Pendulum:** finite amplitude (2.4); a bob of finite size (a physical pendulum with I = mL² + (2/5)mr²: T rises by (1/5)(r/L)² relative, negligible for a compact bob); a stretchy string (L grows under load and at the bottom of the swing); a support that moves (a door frame flexes; a stand rocks); air drag, which damps the amplitude (shifting the mean amplitude during a run) and very slightly lowers the frequency.
- **Damping model:** viscous damping is a convenient fiction. Air drag on a bob is closer to quadratic in speed (amplitude decays linearly in the number of swings at large amplitude); friction at a pivot is dry (constant decrement per swing). The ★ extension distinguishes them.
- **String:** small slopes, uniform tension, no stiffness (a real string has bending stiffness that raises the high harmonics: inharmonicity), no damping.

## 3. Build the computational model

Open `starter/main.py`. `G_STANDARD` is given.

### Step 1: SHM, damping and the log decrement

`shm` is one line. `damped_response` has three branches on ζ (below, equal to, above 1); use `np.cosh`/`np.sinh` in the overdamped branch. `log_decrement_to_zeta` is the inversion of 2.2. `shm_energy` returns ½mv² + ½kx² for arrays.

### Step 2: the pendulum

```python
def pendulum_period_exact(L, theta0, g=G_STANDARD):
    k2 = np.sin(0.5 * theta0) ** 2
    return float(4.0 * np.sqrt(L / g) * ellipk(k2))     # scipy.special.ellipk takes m = k^2
```

`scipy.special.ellipk` takes the *parameter* m = k², not the modulus k: passing sin(θ₀/2) instead of its square is the classic bug, and the small-angle test catches it only weakly (both give π/2 at θ₀ = 0), so the 30° and 90° tests are there to catch it. `pendulum_period_simulated` integrates `pendulum_rhs` with `solve_ivp` and a zero-crossing event. Released from rest, the bob crosses θ = 0 at T/4, 3T/4, 5T/4, …, so the first and third crossings are exactly one period apart. `amplitude_correction` is T/T₀.

### Step 3: the resonance curve

`steady_state_amplitude`, `steady_state_phase` (use `np.arctan2` so the phase passes smoothly through 90°), `quality_factor`, `half_power_bandwidth`, `resonance_peak_frequency` (zero when ζ ≥ 1/√2).

### Step 4: the wave solver

```python
C2 = (c * dt / dx) ** 2
U[1, 1:-1] = U[0, 1:-1] + 0.5 * C2 * (U[0, 2:] - 2 * U[0, 1:-1] + U[0, :-2])
for n in range(1, n_t):
    U[n + 1, 1:-1] = 2 * U[n, 1:-1] - U[n - 1, 1:-1] + C2 * (U[n, 2:] - 2 * U[n, 1:-1] + U[n, :-2])
```

The fixed ends are the zero first and last columns, never updated. `wave_equation_fd` adjusts Δt so that the run ends exactly at `t_end`. `odd_periodic_extension` and `dalembert` are the analytic oracle; `standing_wave` is the mode; `wave_energy` is the discrete energy of 2.8.

### Step 5: reducing pendulum data

`period_from_timings(times, n_swings, u_single)` divides each timing by the number of swings, takes the mean and the Type A term s/√n, and adds the Type B frame term u_single/n_swings in quadrature *after* averaging (it is the same for every repeat, so it does not average down; a test checks this with fifty identical timings). `g_from_pendulum` and `fit_power_law` are 2.6.

`make_synthetic_pendulum_data(seed)` is **SYNTHETIC**: it computes the exact period at θ₀ = 0.15 rad and g = 9.81 m/s² for six lengths, adds a 0.05 % release scatter, and quantises the ten-swing timings to whole frames at 30 fps with a random phase. It exists so that the reduction can be tested on data whose truth is known. It is not a measurement, and nothing from it belongs in your report as data.

Running `main()` prints the damped-oscillator check, the pendulum table of 2.4, the frame-term comparison of 2.5, the reduction of a SYNTHETIC dataset (seed 4061), the string modes and the four FD checks, and writes `oscillations_waves.png` (resonance curves, T/T₀ against amplitude, a pulse reflecting on a fixed–fixed string, and the synthetic T against L).

### Verify

`python -m pytest curriculum/phase-04/lesson-06-oscillations-waves` runs 19 tests: SHM satisfies its ODE and conserves energy; the damped closed form matches `solve_ivp` for ζ = 0, 0.1, 1, 2.5; ζ = 0 reproduces SHM and the log decrement returns ζ; the small-angle limit and T ∝ √L; exact, series and simulated periods agree (and 30° gives 1.74 %, 90° gives 1.180); T grows monotonically with amplitude; the pendulum right-hand side; the resonance curve's static limit, Q, peak position, bandwidth, roll-off and the undamped singularity; the phase; wave speed and harmonics; the standing wave satisfies the PDE and the boundary conditions; the FD solver at C = 1 matches d'Alembert to 10⁻¹², returns the pulse after a round trip and inverts it after half, and conserves energy to 10⁻¹²; d'Alembert splits the pulse into two halves at ±ct; second-order convergence at C = 0.5 and blow-up at C = 1.1; the reduction formulas; g recovered from four synthetic datasets within 3u_g after the amplitude correction.

### Compare

Two independent oracles. For the pendulum, `solve_ivp` with DOP853 at rtol 10⁻¹⁰ and a zero-crossing event returns the period to 10⁻⁷ of the elliptic-integral value at every amplitude tested, so the closed form and the simulation confirm each other. For the string, d'Alembert's formula is exact and needs no grid; the FD solver at C = 1 reproduces it to rounding, which is the strongest possible check that the update, the first step and the boundaries are all right. Below C = 1 the comparison measures the scheme's dispersion, a property of the method rather than a bug.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the raw-data format are in `projects/p04-pendulum/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). Read the safety box first, every session. If you have no kit, run the whole reduction on `make_synthetic_pendulum_data` with your own seed and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p04-pendulum/safety.md`).** Hazards: a swinging mass striking a person, a pet, a screen or a window; the support or its clamp falling; a snapping string releasing the bob; tripping over the string or the tape. Energy limits: bob mass ≤ 0.2 kg, string length ≤ 2.0 m, amplitude ≤ 45° (bob speed at the bottom ≤ 3.4 m/s, kinetic energy ≤ 1.2 J). Controls: a soft, compact bob (a bag of rice, a rubber ball, a padded nut); the full swing arc plus 0.5 m kept clear of people, screens and glass; the support tested with twice the bob's weight before use; the phone on a stand outside the arc; release by hand from beside the arc, never by standing in it. Supervision: none needed for an adult; children only with an adult handling release. Location: indoors, a doorway or a table edge with the swing plane clear. Stop criteria: anyone or anything enters the arc; the support moves or creaks; the string frays; the bob loosens. Not permitted: heavier bobs, glass or metal bobs without padding, swinging over people, beds, desks with screens, or stairs; standing on chairs to rig the support.

The setup, seen from the camera:

![The pendulum setup seen from the camera. A string is clamped at a pivot on a rigid support, with a printed protractor behind the pivot. The length L runs from the pivot to the centre of the bob, displaced by θ from the vertical. The swing arc, plus 0.5 m, is kept clear, and the phone stands 1.5–2 m in front of the swing plane with its lens level with the bob.](figures/pendulum-setup.svg)
*Figure 1. Measure L to the centre of the bob, not to its top, and read the starting and ending amplitude θ on the protractor behind the pivot.*

**Prediction first (in `data/prediction.md`, dated).** For L = 0.40, 0.60, 0.80, 1.00, 1.20, 1.50 m (or whatever your support allows, up to 2 m), compute T₀ = 2π√(L/g) with your local g (or 9.81 ± 0.03 m/s²), its uncertainty from u_L (typically 2 mm, mostly the bob-centre reference), and the frame term per period for N = 10 swings at your frame rate. At L = 1.00 m add the predicted periods at 10°, 20° and 30° from the exact formula, and state, using 2.5, whether each should be detectable. Then freeze the file.

**Procedure in outline** (full version in `procedure.md`): calibrate the frame rate by filming a stopwatch; measure L to the centre of the bob three times; release from rest at about 10° (small-angle series) and time ten full swings from the video, five repeats per length; then at L = 1.00 m repeat at 10°, 20° and 30°, reading the starting and ending amplitude from the protractor; record everything raw. Optional: run a phone accelerometer app with the phone *as* the bob only if the phone is inside a padded pouch and the total mass is within the 0.2 kg limit; otherwise do not.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_timings.csv`, one row per run, never edited):

```
length_id,L_reading_m,run,amplitude_start_deg,amplitude_end_deg,n_swings,start_frame,end_frame,fps_nominal,note
L100,1.003,1,10,8,10,145,749,30,
L100,1.003,2,10,8,10,1021,1624,30,
```

## 5. Validate and improve

Run the reduction on your data (or the SYNTHETIC set) and put three comparisons side by side.

**Prediction against measurement, length by length.** For each L: T ± u_T, g ± u_g (after the amplitude correction), and E_n = |g − g_ref|/√(u_g² + u_ref²). With honest timings at 30 fps and ten swings, u_T is about 1.5 ms and u_g about 0.02–0.05 m/s², larger at short lengths. The SYNTHETIC dataset (seed 4061) shows what consistent data look like: T from 1.2707 s at 0.4 m to 2.4607 s at 1.5 m, u_T 1.5–1.6 ms, and g between 9.776 and 9.791 m/s² *before* the amplitude correction. Those g values are all low by about 0.3 %, consistently, because the synthetic timings were made at θ₀ = 0.15 rad; multiplying by (T/T₀)² = 1.0028 brings the power-law estimate from 9.781 to 9.809 m/s², against the 9.81 used to generate it. That is the pattern to recognise in real data: a *consistent* offset in g across lengths points at a systematic effect that does not depend on L (amplitude, a timing bias) or that scales with L (a stretchy string).

| Pattern in g across lengths | Likely cause | Test |
|---|---|---|
| all low by ≈ 0.3–1 %, no trend | amplitude correction not applied, or larger amplitude than recorded | recompute with the measured mean amplitude; repeat one length at 5° |
| low at short L, fine at long L | L measured to the knot or the top of the bob (constant offset) | re-measure to the bob's centre; fit T² = (4π²/g)(L + ΔL) for ΔL |
| scatter far above the Type A estimate | inconsistent frame rule for "passes the lowest point" | re-read frames with one fixed rule; use the passage through the vertical, not the turning point |
| g drifts within a session | support creeping, string stretching | re-measure L at the end of the session |

**Amplitude against theory.** At L = 1 m, plot T/T₀ against the mean amplitude with the exact curve of 2.4. At 30°, the predicted 1.74 % (35 ms on 2.006 s) is about 25 times the frame term per period with ten swings, so it must show; if it does not, your amplitude reading or swing count is wrong. At 10° the predicted 0.19 % is about three times the frame term: expect it to be marginal on a single run and visible on the mean of five.

**Code against closed form.** The FD wave solver and the pendulum ODE have been checked against exact solutions (part 3). Hand-check one length: T, u_T, g and u_g with a calculator from your raw frames.

**Improve and re-run.** The dominant terms at 30 fps are the frame term (cut it by counting more swings or using 120 or 240 fps) and the length reference (cut it by measuring to the bob's centre with a ruler across the bob). Re-run one length after the improvement and put both results in the report.

## 6. Defend your work

Submit: `data/prediction.md` (dated before any timing), `data/calibration.csv`, `data/raw/*.csv`, your analysis script, the two plots with k = 2 error bars, the g table with E_n, the amplitude comparison, `oscillations_waves.png`, and a report following `projects/p04-pendulum/report-template.md`. If you used the SYNTHETIC dataset, every table and figure says so in its caption.

### Adversarial self-checks

- **What if the bob's mass doubles?** The period does not change (mass cancels in θ̈ = −(g/L)sin θ), but a heavier bob stretches the string more and loads the support more: check that L was measured *with* the bob hanging.
- **What if the amplitude goes to zero?** T → T₀ exactly, and the measurement gets harder (the bob moves slowly at the turning points, but its passage through the vertical stays sharp). What happens to damping-induced amplitude drift?
- **Which assumption fails first?** On the bench: small amplitude above about 10°, then the length reference. On the spacecraft: viscous damping (structural damping is not proportional to velocity), and linearity at large deflection.
- **What measurement could falsify the model?** A fitted exponent p far from 0.5 across lengths (with L measured correctly), or a period that does not lengthen with amplitude at 30°.
- **Does the result survive a change in method?** Compare video timing with a phone accelerometer or with a stopwatch over 20 swings: the stopwatch's reaction bias is a constant per run, so its effect per period shrinks as 1/N.
- **Which decision changes if the uncertainty doubles?** At twice the frame term, the 10° amplitude effect is no longer detectable, and a claim to have seen it must be withdrawn; the 30° effect still is.

## Check yourself

1. Write the small-angle period of a pendulum. Which two quantities does it not depend on?
2. Define the damping ratio ζ and the quality factor Q, and give the amplitude at resonance relative to the static deflection.
3. A lightly damped oscillator's successive peaks have ratio 1.37. What is ζ?
4. By what percentage does a pendulum's period exceed T₀ at 30° amplitude? At 5°?
5. Give the wave speed on a string and the frequencies of its standing waves with both ends fixed.
6. What is the CFL condition for the explicit wave scheme, and what happens at Courant number exactly 1?
7. Why does a pulse come back inverted after reflecting from a fixed end?
8. You time ten swings at 30 fps. What is the frame-quantisation uncertainty of the period, and does it shrink if you repeat the run five times?
