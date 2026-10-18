# P06.13  Single-DOF Vibration

> One mass, one spring, one damper: tap it to learn its frequency and damping, then shake its base and those two numbers tell you how hard it will ring and how much its bracket will bend.

| | |
|---|---|
| **Learning outcomes** | 1. Reduce a component on a flexible bracket to a single-degree-of-freedom model: stiffness from beam theory, an effective mass that includes the bracket's own mass, natural frequency and damping ratio. 2. Identify the damping ratio from a ring-down record by the logarithmic decrement over n cycles and from a frequency response function by the half-power bandwidth, and say how noise and frequency resolution bias each. 3. Derive the absolute and relative transmissibilities of a base-excited system and the full transient response to a base sine switched on at t = 0, and check the closed form against `solve_ivp`. 4. Predict the response acceleration, the bracket deflection and the root stress of a component in a sine-vibration test, and decide whether stiffening or damping is the better fix. 5. Distinguish viscous from Coulomb (dry-friction) damping by the shape of the decay, and compute the amplitude-dependent equivalent viscous damping of friction. |
| **Prerequisites** | `P02.02` (free response, damping regimes, logarithmic decrement), `P02.03` (forced response, complex gain, magnification, transmissibility), `P02.07` (phasors), `P04.06` (damped and driven oscillators), `P05.12` (cantilever deflection PL³/3EI), `P05.10` (bending stress Mc/I), `P06.11` (Lagrange; the Rayleigh energy argument), `P03.06` (`solve_ivp`, events). Skill check: you can write the steady response of m ẍ + c ẋ + k x = F₀ cos ωt as a complex amplitude, and compute the tip stiffness of a cantilever. |
| **Engineering question** | A 0.6 kg camera sits at the tip of a 6061-T6 aluminium cantilever bracket, 100 × 40 × 6 mm (illustrative). Before a 2 g (0-peak) sine-vibration test of the spacecraft panel it is mounted on, what is the bracket's natural frequency, what damping does a tap test give, how many g does the camera see at resonance, how much does the bracket bend and how close does its root stress come to yield, and is it better to stiffen the bracket or to add damping? The answers: 78.3 Hz; ζ ≈ 0.020 from a (SYNTHETIC) ring-down; 50 g; 2.03 mm and 126 MPa, a factor of 2.2 on a typical 276 MPa yield; at ζ = 0.01 the factor falls to 1.1, and an 8 mm bracket (120 Hz) cuts the stress to 71 MPa. |
| **Assumptions and limits** | One mode (the bracket's first bending mode; the next is several times higher); a linear spring and linear viscous damping with constant coefficients (Coulomb friction treated separately); small deflections; the base motion is prescribed and is not affected by the component (a stiff, heavy shaker); steady sinusoidal input (a slow sweep); the root stress from the tip-load bending moment, ignoring stress concentrations at the bolts. The camera is a rigid point mass. Damping in real bolted joints depends on amplitude, and the value from a gentle tap may differ from the value at test level. |
| **Analytical work** | The cantilever's stiffness 3EI/L³ and the Rayleigh effective mass m_tip + (33/140)m_beam; the free response in all regimes; the logarithmic decrement and its exact inversion; receptance, mobility and accelerance and their limits; the half-power estimate of ζ; base excitation in relative coordinates; absolute and relative transmissibility; the transient response to a base sine switched on at t = 0 and the build-up time; deflection to root stress; Coulomb decay (4F/k per cycle) and its equivalent viscous damping 2(F/k)/(πX). |
| **Computational deliverable** | `starter/main.py`: `cantilever_stiffness`, `effective_mass`, `natural_frequency`, `damping_ratio`, `free_response`, `log_decrement`, `zeta_from_log_decrement`, `frf`, `half_power_zeta`, `transmissibility`, `relative_transmissibility`, `base_sine_response`, `simulate_base_sine`, `bracket_root_stress`, `coulomb_free_response`, `simulate_coulomb`, `coulomb_equivalent_zeta` (`camera_bracket`, `find_peaks_refined`, `synthetic_tap_test` and `main` are given). Tests prove: the stiffness scales as L⁻³ and the effective mass as 33/140; the free response in five regimes meets its initial conditions and matches `solve_ivp` to 10⁻¹⁰; the log decrement recovers ζ to 10⁻⁴ from clean records and within 0.0015 from a noisy SYNTHETIC one; the FRFs have their static, resonant (1/(2ζk), −90°) and mass-line limits and the mobility peaks at 1/c; the half-power method recovers ω_n and ζ to O(ζ²); transmissibility is 1 at r = √2, peaks at ≈ 1/(2ζ) and reverses its damping trend above √2; the closed-form base-sine response from rest matches `solve_ivp` to 10⁻⁸ and tends to the steady amplitudes; the root stress is Mc/I; the Coulomb decay falls 2F/k per half period at the undamped frequency, sticks, and matches an event-driven `solve_ivp` to 10⁻¹¹ m; the equivalent viscous ζ of friction is inversely proportional to amplitude. |
| **Practical activity** | **computational experiment**: a seeded study of the logarithmic decrement's bias and scatter against the number of cycles used (on SYNTHETIC records); the half-power method against frequency resolution; a damping sweep and a stiffness sweep of the sine-test response; and a test that tells Coulomb from viscous decay. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Lesson `P06.19` measures a real pendulum's damping with the same tools.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `sdof_vibration.png`, the derivations of the effective mass, the transmissibilities and the base-sine transient, the experiment tables with your predictions, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the base-excitation equations are derived (not quoted) and the closed form is shown to agree with `solve_ivp`; all public tests pass; ζ is identified two ways with a stated uncertainty; the note gives the resonant g level, the deflection and the stress margin, and recommends stiffening or damping with numbers. |
| **References** | `@mit-ocw-2-003sc-engineering-dynamics` Lecture 19: Introduction to Mechanical Vibration; Lecture 20: Linear System Modeling a Single Degree of Freedom Oscillator (experimental determination of damping ratio; steady-state response to harmonic input); Recitation 10 Notes: Natural Frequency From Deflection & Frequency Response (required-free). `@mit-ocw-16-07-dynamics` L19 "Vibration, normal modes, natural frequencies, instability" (supplementary-free). `@libretexts-cline-variational-principles` 3.5 Linearly-damped Free Linear Oscillator; 3.6 Sinusoidally-driven, linearly-damped, linear oscillator (required-free). `@openstax-university-physics-vol1` 15.5 Damped Oscillations; 15.6 Forced Oscillations (supplementary-free). `@mit-ocw-18-03-differential-equations` notes 14 "Complex gain, dashpot drive"; 16 "Frequency response" (supplementary-free). `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (required-free). |
| **Extensions** | ★ Add hysteretic (structural) damping, k(1 + iη), to `frf` and show that the half-power method then returns ζ ≈ η/2, and that the free decay of a hysteretic model is not defined in the time domain. ★ Fit ω_n, ζ and an amplitude to the whole SYNTHETIC ring-down with `scipy.optimize.curve_fit` and compare the scatter over 20 seeds with the log decrement's (preview of `P06.19`). ★ Model a sine sweep at 2 octaves per minute through the resonance and find how much of the steady peak the camera actually sees. |
| **Time** | unmeasured |

## 1. The engineering problem

Every component that flies on a launch vehicle or a spacecraft is shaken on the ground first. In a **sine-vibration test** the component, or the panel it is bolted to, sits on a shaker table that moves sinusoidally, sweeping slowly through a frequency band at a specified base acceleration. The test is meant to show that the hardware survives the low-frequency vibration of launch. It can also break hardware that would have flown happily, if a resonance of the test article sits in the band and is barely damped: the response at resonance can be 20 to 100 times the input. Test engineers therefore predict each component's response before the test, decide where the input must be reduced ("notched") to protect it, and size brackets so that resonances are either above the band or tolerable within it.

The prediction starts with the simplest model in vibration: one mass on one spring with one damper, driven through its base. `P02.03` derived the forced response of that system and its transmissibility; this lesson uses them as an engineer does. It builds the model from the hardware (beam theory for the spring, an energy argument for the mass), identifies the damping from data (a tap test and a frequency response), predicts the test response in acceleration, deflection and stress, and checks the closed-form transient against `solve_ivp`. It ends with the damping mechanism that most often breaks the linear model, dry friction in joints.

The hardware (illustrative): a 0.6 kg camera at the tip of a 6061-T6 cantilever bracket, 100 mm long, 40 mm wide and 6 mm thick, bending about its thin direction, with a 2 g (0-peak) sine input at the bracket root. Success criteria, written first: the closed-form free and forced responses agree with `solve_ivp` to 10⁻⁸ of the response scale; the log decrement returns the damping that generated a clean record to 10⁻⁴; the half-power method returns it to O(ζ²); and the Coulomb closed form agrees with an event-driven integration to 10⁻¹¹ m.

## 2. Understand and derive

### 2.1 The model from the hardware

**Stiffness.** A tip load P deflects a cantilever by PL³/(3EI) (`P05.12`), so the bracket is a spring of

**k = 3EI/L³.**   (1)

With I = bt³/12 = 40 × 6³/12 = 720 mm⁴ = 7.2 × 10⁻¹⁰ m⁴ and E = 68.9 GPa (typical), k = 1.4882 × 10⁵ N/m. Because k ∝ t³/L³, a 33 % thicker bracket is 2.4 times stiffer, and a bracket twice as long is eight times softer.

**Effective mass.** The bracket's own mass (64.8 g) moves too, but not all of it at the tip velocity. Rayleigh's idea: assume the vibrating shape is the static deflection shape under a tip load, w(x) = w_tip (3x²L − x³)/(2L³), compute the kinetic energy of the beam moving in that shape, ½∫ρA ẇ² dx, and write it as ½ m_eq ẇ_tip². The integral gives m_eq = (33/140) m_beam, so

**m = m_tip + (33/140) m_beam.**   (2)

Here m = 0.6 + 0.2357 × 0.0648 = 0.6153 kg. Then

**ω_n = √(k/m) = 491.82 rad/s, f_n = 78.27 Hz** (79.26 Hz if the bracket's mass is ignored; 1.3 % high).

Lesson `P06.15` computes the same frequency from the continuous beam and shows how good (2) is.

**Damping.** Write the damper as c = 2ζ√(km) with the **damping ratio** ζ. A bolted aluminium bracket typically has ζ between 0.005 and 0.05; it is not predictable from drawings and must be measured.

### 2.2 Identifying ζ from a ring-down

Tap the camera and let it ring. For ζ < 1 the free response is x = e^{−ζω_n t}(x₀ cos ω_d t + ((v₀ + ζω_n x₀)/ω_d) sin ω_d t) with ω_d = ω_n√(1 − ζ²) (`P02.02`; Cline 3.5; the critically damped and overdamped forms are in `free_response` too). Successive positive peaks are one damped period T_d = 2π/ω_d apart and their ratio is e^{ζω_n T_d}. The **logarithmic decrement** over n cycles is

**δ = (1/n) ln(x₀/x_n) = 2πζ/√(1 − ζ²),  so  ζ = δ/√(4π² + δ²)**   (3)

exactly (2.003SC Lecture 20 treats the experimental determination of the damping ratio). Using n cycles rather than one averages out noise on individual peaks, but if n is so large that the last peak is buried in noise, the measured peak is biased upwards by the noise and ζ comes out low. Part 4 measures both effects on SYNTHETIC records.

### 2.3 Identifying ζ from a frequency response

A shaker or an instrumented hammer measures a **frequency response function** (FRF): the response per unit force at each frequency. For m ẍ + cẋ + kx = F e^{iωt} the three common forms are

**receptance X/F = 1/(k − mω² + icω), mobility V/F = iω X/F, accelerance A/F = −ω² X/F.**   (4)

Their limits are worth knowing by heart: the receptance is 1/k at low frequency (static), 1/(2ζk) at ω_n with a −90° phase, and falls as 1/(mω²); the accelerance tends to 1/m at high frequency (the "mass line"); the mobility peaks exactly at ω_n with the value 1/c. The **half-power method**: find the two frequencies ω₁ < ω₂ at which |H| = |H|_max/√2; for light damping

**ζ ≈ (ω₂ − ω₁)/(2ω_n),**   (5)

with an error of order ζ². The half-power bandwidth of the camera bracket at ζ = 0.02 is 2ζf_n = 3.13 Hz, so an FRF must be sampled much more finely than that to resolve it (part 4).

### 2.4 Base excitation

On the shaker the base moves, y(t), and the camera moves x(t). The spring and damper act on the **relative** motion z = x − y:

m ẍ = −c(ẋ − ẏ) − k(x − y)  ⇒  **m z̈ + c ż + k z = −m ÿ.**   (6)

The base acceleration acts as an inertial force on the mass. For a harmonic base, ÿ = A e^{iωt}, the steady relative motion is z = Z e^{iωt} with Z = −A/(ω_n² − ω² + 2iζω_nω). Two ratios matter.

**The relative transmissibility** gives the bracket's deflection, which sets its stress. Per unit base displacement (Y = −A/ω²),

**|Z/Y| = r²/√((1 − r²)² + (2ζr)²),  r = ω/ω_n,**   (7)

which is 0 at low frequency (the bracket carries the camera rigidly), 1/(2ζ) at resonance, and 1 at high frequency (the camera stands still while the base moves under it). Per unit base acceleration, |Z| = (A/ω_n²)/√((1 − r²)² + (2ζr)²): the static deflection under the base acceleration, magnified.

**The absolute transmissibility** gives the camera's acceleration (and the force the bracket puts into the base). The mass's absolute acceleration is ẍ = −(cż + kz)/m, so

**|X/Y| = |ẍ/ÿ| = √((1 + (2ζr)²)/((1 − r²)² + (2ζr)²)),**   (8)

the result of `P02.03` §2.7: 1 at low frequency, about 1/(2ζ) at resonance, exactly 1 at r = √2 for every ζ, and falling as 2ζ/r above it, where *more* damping transmits *more*.

**The camera bracket in the test.** At resonance, with A = 2g = 19.61 m/s² and ζ = 0.02: the peak absolute transmissibility is 25.02, so the camera sees **50.0 g**; the peak relative deflection is A/(ω_n² × 2ζ√(1 − ζ²)) = **2.028 mm** (the static deflection under 1 g is only 0.041 mm).

### 2.5 From deflection to stress

A tip deflection z of the cantilever needs a tip force kz and puts a root moment kzL into the bracket. With M = kzL = 3EIz/L² and σ = M(t/2)/I (`P05.10`):

**σ_root = 3Etz/(2L²).**   (9)

For z = 2.028 mm: **σ = 125.7 MPa**, against a typical 6061-T6 yield of 276 MPa (a typical value, not a design allowable): a factor of 2.2. Equation (9) says what to change: at a fixed base acceleration, z ∝ 1/(ζω_n²) and σ ∝ tz, so stiffening helps through ω_n² ∝ t³ and damping through 1/ζ.

### 2.6 The transient: switching the sine on

A real test does not start in steady state. With the base sine switched on at t = 0 and the camera moving with the base (z(0) = ż(0) = 0), the solution of (6) is the steady part plus the free response that cancels its initial values:

**z(t) = Im(Z e^{iωt}) + z_h(t),  z_h = free response with z_h(0) = −Im Z, ż_h(0) = −Re(ωZ).**   (10)

At resonance the amplitude builds up as (1 − e^{−ζω_n t}), reaching 99 % after ln(100)/(ζω_n) = 0.468 s, or 37 cycles. A sweep that passes through the half-power band (3.13 Hz wide) faster than that never reaches the steady peak; that is why sine-sweep rates are specified, and extension 3 quantifies it. Equation (10) is the closed form that part 5 checks against `solve_ivp`.

### 2.7 Coulomb damping

Bolted joints, cable bundles and pivots dissipate energy by dry friction: a force of constant size F opposing the velocity. The free motion m ẍ + kx = −F sign(ẋ) is linear within each half cycle: an undamped oscillation about the offset +F/k or −F/k. Released from rest at x₀, the mass reaches −(x₀ − 2F/k) after half a period, then x₀ − 4F/k, and so on:

**the amplitude falls by 4F/k per cycle, linearly, at the undamped frequency ω_n,**   (11)

and the motion stops at the first turning point inside |x| ≤ F/k, where the spring can no longer overcome friction. Contrast viscous damping, whose amplitude falls by a constant *ratio* per cycle. Plot the peaks: a straight line in linear axes is friction, a straight line in log axes is viscous. Equating the energy lost per cycle, 4FX for friction and πcω X² for a viscous damper, gives the **equivalent viscous damping** at amplitude X near resonance,

**ζ_eq = 2(F/k)/(πX),**   (12)

inversely proportional to amplitude: friction damps small vibrations strongly and large ones weakly, the opposite of what you want in a test. A damping ratio measured with a gentle tap can overestimate the damping at test level.

### 2.8 Where the model stops being valid

- **One mode.** The bracket's second bending mode and the camera's own modes are higher, but a base input that reaches them excites them; `P06.14` handles several.
- **Linear damping.** Real joints mix viscous-like material damping with friction; ζ varies with amplitude. Identify ζ at the level of the test, or bracket it.
- **A stiff, heavy base.** The model prescribes ÿ. A light panel is itself a spring–mass, and the camera's reaction changes the panel's motion; the two must then be modelled together (`P06.14`).
- **Linear elasticity.** At 126 MPa the bracket is elastic; at ζ = 0.005 the predicted 503 MPa is beyond yield, where the linear model has already failed (the bracket would yield or the joint slip first).
- **Stress concentrations.** Bolt holes and fillets at the root raise the local stress by a factor of 2–3; (9) is the nominal stress.

## 3. Build the computational model

Open `starter/main.py`. `camera_bracket`, `find_peaks_refined` (one peak per 0.7 period, parabola-refined), `synthetic_tap_test` (SYNTHETIC, seeded) and `main` are given.

### Step 1: the model and the free response

`cantilever_stiffness` is (1), `effective_mass` is (2), `natural_frequency` and `damping_ratio` are one line each. `free_response` covers the three regimes, with |ζ − 1| < 10⁻⁹ treated as critical.

### Step 2: identification

`log_decrement(peaks, n)` and `zeta_from_log_decrement` are (3). `frf` is (4) with a `kind` switch. `half_power_zeta` walks down from the peak to the first samples below |H|_max/√2 on each side, interpolates linearly, and applies (5).

### Step 3: base excitation

`transmissibility` (8) and `relative_transmissibility` (7). `base_sine_response(t, A, ω, ω_n, ζ)` is (10): the steady complex amplitude, then `free_response` for the homogeneous part, and its derivative in closed form; it returns z and the absolute acceleration −(2ζω_n ż + ω_n² z). `simulate_base_sine` integrates (6) with DOP853 (rtol 10⁻¹⁰, a maximum step of a tenth of the shorter period). `bracket_root_stress` is (9).

```python
Z = -A / (wn**2 - omega**2 + 2j * zeta * wn * omega)
zs = np.imag(Z * np.exp(1j * omega * t))                 # steady part
zh = free_response(t, -np.imag(Z), -np.real(omega * Z), wn, zeta)   # cancels z(0), z'(0)
```

### Step 4: Coulomb damping

`coulomb_free_response` steps from turning point to turning point: each half cycle is a cosine about ±F/k, the next turning point is 2(±F/k) − x, and the mass sticks when |x| ≤ F/k. `simulate_coulomb` integrates with the friction sign fixed until an event finds ẋ = 0, flips the sign and restarts. `coulomb_equivalent_zeta` is (12).

Running `main()` prints

```
bracket: I = 7.2000e-10 m^4, k = 1.4882e+05 N/m, beam mass 64.8 g, effective mass 0.6153 kg
natural frequency 491.82 rad/s = 78.27 Hz (tip mass only: 79.26 Hz)

SYNTHETIC tap test (seed 13, zeta 0.02 used to generate it):
   log decrement over  1 cycles: zeta = 0.01685
   log decrement over  5 cycles: zeta = 0.01938
   log decrement over 10 cycles: zeta = 0.01957
   log decrement over 20 cycles: zeta = 0.01862
   damped frequency from 20 peak intervals: 78.387 Hz

half-power method on the receptance (4001 points, 0.5-1.5 w_n): w_n 491.520 rad/s, zeta 0.02002

sine test 2.0 g (0-peak) base: peak transmissibility 25.025 -> 50.0 g at the camera; peak relative deflection 2.028 mm; root stress 125.7 MPa (typical 6061-T6 yield 276 MPa)
   resonant dwell from rest, 0.5 s: closed form vs solve_ivp max |dz| 7.9e-13 m (peak 2.012 mm); max |d accel| 1.9e-07 m/s^2
   20 Hz from rest: max |dz| 5.3e-14 m

Coulomb decay (d = F/k = 2.0 um, x0 = 100 um): closed form vs solve_ivp max diff 9.2e-15 m; equivalent zeta at 100 um 0.0127, at 20 um 0.0637
```

and writes `sdof_vibration.png`: the SYNTHETIC tap-test record with its detected peaks; the absolute transmissibility for four damping ratios with the relative transmissibility and the r = √2 line; and a Coulomb decay (closed form and `solve_ivp`) against a viscous decay with the same initial energy loss per cycle.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-13-single-dof-vibration`. The tests:

1. **Model**: (1) and its L⁻³ scaling; the 33/140 share; ω_n and ζ for simple numbers.
2. **Free response** at ζ = 0, 0.05, 0.7, 1 and 2.5: initial position and velocity, and agreement with `solve_ivp` to 10⁻¹⁰ m.
3. **Log decrement**: ζ = 0.005, 0.02, 0.1 recovered to 10⁻⁴ from clean records; the exact inversion of (3); ζ = 0.02 within ±0.0015 from a SYNTHETIC noisy record (seed 4).
4. **FRF**: static 1/k; 1/(2ζk) at −90° at resonance; the mass line 1/m; mobility = iω × receptance, peaking at 1/c; an unknown kind is refused. The half-power method recovers ω_n and ζ to O(ζ²) for ζ = 0.01, 0.02, 0.05.
5. **Transmissibility**: 1 at r → 0 and at r = √2 for three ζ; peak ≈ 1/(2ζ); the damping trend reverses across √2; the relative transmissibility's three limits.
6. **Base sine**: the closed form (10) against `solve_ivp` at r = 0.3, 1 and 2.2; after 10–12 time constants the amplitudes of z and ẍ equal the steady values from (7) and (8) to 10⁻³.
7. **Stress**: (9) equals kzL(t/2)/I.
8. **Coulomb**: turning points at multiples of half the undamped period, falling 2F/k each; agreement with the event-driven integration; sticking; ζ_eq ∝ 1/X and the energy balance behind (12).
9. `main` writes the figure.

### Compare

The closed form (10) and `solve_ivp` agree to 7.9 × 10⁻¹³ m in a 0.5 s resonant dwell (a response of 2 mm: relative agreement 4 × 10⁻¹⁰) and to 5.3 × 10⁻¹⁴ m at 20 Hz; the Coulomb closed form and the event-driven integration agree to 9.2 × 10⁻¹⁵ m. The half-power method applied to the exact receptance returns ζ = 0.02002 and ω_n 0.06 % low, the O(ζ²) error that (5) predicts. On the identification side the established tool is `scipy.signal.find_peaks`, which `find_peaks_refined` wraps: without its `distance` argument the first version of this lesson's peak finder counted noise ripples as cycles and returned a "damped frequency" of 197 Hz from a 78 Hz record. A peak finder must know roughly how often a peak can occur.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All records here are SYNTHETIC (generated by `synthetic_tap_test` from a seed, with ζ = 0.02 and 1 μm of Gaussian noise on a 100 μm release), and every response is a simulation, not a measurement. Write each prediction before running it.

1. **The log decrement against n.** For 20 seeds (0–19), estimate ζ from n = 1, 3, 5, 10, 20, 30 cycles. Predict which n is best. *Finding:* mean (standard deviation) 0.01881 (0.00285), 0.01928 (0.00093), 0.01934 (0.00071), 0.01960 (0.00035), 0.01875 (0.00037), 0.01701 (0.00051). The scatter falls as n grows, but beyond 10 cycles the bias grows: the amplitude after 20 cycles is 8.1 μm and after 30 only 2.3 μm, comparable with the 1 μm noise, and a noisy peak is more often too high than too low (the peak finder picks the largest sample). Use the cycles in which the signal is at least ten times the noise. The noiseless record gives 0.020000 at any n, and a damped frequency of 78.2592 Hz against the exact 78.2592 Hz.
2. **The half-power method against resolution.** Sample the receptance from 60 to 100 Hz in steps of 0.05, 0.25, 0.5, 1 and 2 Hz (as a stepped-sine test would) and apply (5). *Finding:* ζ = 0.02002, 0.02005, 0.02059, 0.02085, 0.02035, and f_n = 78.228, 78.227, 78.228, 78.212, 78.180 Hz. With a 3.13 Hz half-power bandwidth, steps of 0.25 Hz or less (about 12 per bandwidth) give ζ to 0.3 %; at 1–2 Hz the estimate is off by several per cent and unpredictable in sign, because the true peak falls between samples.
3. **Damping sweep.** The 2 g resonant response for ζ = 0.005, 0.01, 0.02, 0.05, 0.1. Predict the stress at ζ = 0.01. *Finding:*

   | ζ | camera (g) | deflection (mm) | root stress (MPa) | factor on 276 MPa | 99 % build-up |
   |---|---|---|---|---|---|
   | 0.005 | 200.0 | 8.109 | 502.8 | 0.55 | 1.87 s (147 cycles) |
   | 0.01 | 100.0 | 4.055 | 251.4 | 1.10 | 0.94 s (73 cycles) |
   | 0.02 | 50.0 | 2.028 | 125.7 | 2.20 | 0.47 s (37 cycles) |
   | 0.05 | 20.1 | 0.812 | 50.3 | 5.48 | 0.19 s (15 cycles) |
   | 0.1 | 10.2 | 0.407 | 25.3 | 10.9 | 0.09 s (7 cycles) |

   The answer depends on a number nobody can predict from the drawing: if the joint gives 1 % instead of 2 %, the bracket is at yield.
4. **Stiffness sweep.** Thickness 5, 6, 8 and 10 mm at ζ = 0.02. *Finding:* f_n 59.67, 78.27, 120.02, 167.04 Hz; deflection 3.488, 2.027, 0.862, 0.445 mm; stress 180.3, 125.7, 71.3, 46.0 MPa. Thickening to 8 mm (+33 % bracket mass, 21.6 g) halves the stress and moves the resonance to 120 Hz; if the test band ends at 100 Hz, it leaves the band altogether.
5. **Friction or viscous?** A Coulomb decay with F/k = 2 μm from 100 μm, and a viscous decay with the same first-cycle loss (ζ = ζ_eq = 0.0127). Tabulate the peaks, their differences and their ratios. *Finding:* Coulomb peaks 100, 92, 84, …, 12 μm: constant differences of 8.000 μm (4F/k), ratios falling from 0.920 to 0.600; the motion stops after 12.25 cycles. Viscous peaks 100, 92.31, 85.21, …, 41.48 μm: constant ratios of 0.9231, differences falling from 7.69 to 3.46 μm. ζ_eq of the friction is 0.0127, 0.0255, 0.0637 and 0.127 at 100, 50, 20 and 10 μm.

Raw data format: CSV tables `n_cycles, zeta_mean, zeta_std` (1); `step_Hz, fn_Hz, zeta` (2); `zeta, camera_g, z_mm, sigma_MPa, factor, buildup_s` (3); `t_mm, fn_Hz, z_mm, sigma_MPa` (4); `cycle, coulomb_um, viscous_um` (5); plus `sdof_vibration.png`. Label every table derived from `synthetic_tap_test` SYNTHETIC.

## 5. Validate and improve

| quantity | closed form (§2) | numerical | difference | cause |
|---|---|---|---|---|
| free response, five regimes | `free_response` | `solve_ivp` (tests) | < 10⁻¹⁰ m | solver tolerance |
| resonant dwell from rest, 0.5 s | (10), peak 2.012 mm | `solve_ivp` | 7.9 × 10⁻¹³ m | solver tolerance |
| steady deflection at resonance | A/(2ζ√(1 − ζ²)ω_n²) = 2.028 mm | (10) after 12 time constants (tests) | < 10⁻³ | residual transient |
| peak transmissibility | ≈ 1/(2ζ) = 25 | 25.025 (sampled maximum of (8)) | 0.1 % | the (2ζr)² term in the numerator |
| ζ from the receptance | 0.02 | half-power 0.02002 | 0.1 % | O(ζ²) of (5) |
| ζ from a SYNTHETIC ring-down (10 cycles, 20 seeds) | 0.02 | 0.01960 ± 0.00035 | −2 % | noise bias on late peaks |
| Coulomb decay | (11) | event-driven `solve_ivp` | 9.2 × 10⁻¹⁵ m | — |

The required independent check of this lesson, closed form against `solve_ivp`, holds to the solver tolerance for the free response, the base-sine transient and the Coulomb decay. The engineering answer: **the camera bracket resonates at 78.3 Hz; with ζ = 0.02 a 2 g sine input gives 50 g at the camera, a 2.03 mm bracket deflection and a 126 MPa root stress, a factor of 2.2 on yield.** The weak link is ζ, which the model cannot supply: at ζ = 0.01 the factor is 1.1. Stiffening (8 mm: 71 MPa, 120 Hz) is the more robust fix, because it does not rely on a damping value that may change with amplitude, temperature and bolt torque; adding damping (a viscoelastic washer) is attractive only if its ζ is measured at test level.

What to improve first: measure ζ at the test amplitude, not with a tap (2.7 shows friction damping falls with amplitude), and include the stress concentration at the root. What not to improve: the closed form, which is exact for the model.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (6)–(10) and (12); `sdof_vibration.png`; the five experiment tables with your predictions (SYNTHETIC data labelled); and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving f_n, ζ with its uncertainty, the resonant g level, deflection and stress, and your recommendation.

### Adversarial self-checks

- **What if the input level doubles?** In the linear model everything doubles: 100 g, 4.06 mm, 251 MPa, a factor of 1.1. In reality a bolted joint behaves differently at higher level: the equivalent damping of pure friction falls with amplitude (12), while microslip over a larger contact area adds dissipation, so ζ can move either way. The linear prediction is a guide, not a guarantee.
- **What if ζ goes to zero?** The steady response at resonance is infinite and the build-up time ln(100)/(ζω_n) infinite: the response grows linearly in time (`P02.03`'s exceptional case), and the test duration, not the damping, limits it.
- **Which assumption fails first in the real world?** Linear, amplitude-independent damping. Then the rigid, heavy base: a light panel and the camera bracket form a two-degree-of-freedom system (`P06.14`).
- **What measurement could falsify the model?** A low-level sine sweep with an accelerometer on the camera: the measured f_n should be within a few per cent of 78 Hz (the clamp's compliance at the root lowers it) and the peak transmissibility should equal 1/(2ζ) with the ζ from a tap test at the same level. A peak that rises less than proportionally with input level means amplitude-dependent damping.
- **Does the result survive a change of method?** The closed form and `solve_ivp` agree to 10⁻¹⁰ relative; ζ from the ring-down and from the FRF agree within the ring-down's bias (experiment 1).
- **Which decision changes if the uncertainty doubles?** With ζ known only as 0.02 ± 0.01 (not ± 0.005), the stress range spans 84–251 MPa: the 6 mm bracket can no longer be cleared for test without a notch, while the 8 mm bracket (48–143 MPa) still can.

## Check yourself

1. Why does the effective mass of a cantilever with a tip mass include only 33/140 of the bracket's own mass?
2. Write the logarithmic decrement over n cycles and its exact relation to ζ. Why not always use as many cycles as possible?
3. What are the low-frequency, resonant and high-frequency limits of the receptance and of the accelerance?
4. Write the equation of motion of a base-excited mass in relative coordinates. What plays the role of the force?
5. At what frequency ratio is the absolute transmissibility 1 for every ζ, and what does more damping do above it?
6. Which transmissibility sets the bracket's stress, and which sets the camera's acceleration?
7. How long does a resonant response take to build up to 99 % of its steady value, and why does that matter for a sine sweep?
8. How can you tell Coulomb damping from viscous damping in a decay record, and how does the equivalent ζ of friction depend on amplitude?
