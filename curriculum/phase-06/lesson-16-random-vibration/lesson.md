# P06.16  Random Vibration

> A random test has no amplitude, only a spectrum. Multiply the spectrum by |H|², integrate, take the square root, and you have the one number the box must survive, and then remember that "3σ" is not the largest peak it will see.

| | |
|---|---|
| **Learning outcomes** | 1. Read a random-vibration test specification given as PSD breakpoints and dB/octave slopes, and integrate it exactly to its overall Grms. 2. Derive the response PSD of a base-excited single-degree-of-freedom component, |H(f)|²G(f), the white-noise integral that gives Miles' equation, and the relative-displacement RMS; say exactly what Miles' equation leaves out. 3. Implement and test a seeded synthesizer of Gaussian random acceleration with a target PSD and an exact (first-order-hold) discretization of the SDOF model, and use them to check the analytic Grms by an independent time-domain simulation with a derived statistical tolerance. 4. Compute the 3σ load, the Gaussian exceedance fraction and Davenport's expected largest peak, and decide what load a mount must be designed for in a test of given duration. |
| **Prerequisites** | `P02.19` (stationary random processes, one-sided PSD, Welch estimates, the variance of a variance estimate), `P02.10` (FFT and sampling), `P02.16` (Gaussian random variables), `P02.08` (transfer functions, |H(iω)|), `P02.05` (matrix exponential and discretization), `P06.13` (SDOF base excitation and transmissibility), `P06.14` (modes; a box is many SDOF systems). Skill check: you can say why the variance of the output of a linear filter driven by a stationary input is ∫|H|²S df. |
| **Engineering question** | A 2 kg avionics box whose first mode is 150 Hz with Q = 10 is bolted to a spacecraft panel and must pass the NASA-STD-7001C minimum-workmanship random test (Table 2, overall 6.78 Grms) for 60 s. What is its RMS response, what 3σ load must its mounts take, how far does it move on its mounts, and is Miles' equation good enough? The answers: 9.66 Grms by the full integral (Miles 9.71), 29.0 g and 568 N at 3σ, 0.106 mm RMS relative displacement, and an expected largest peak of 4.40σ (42.5 g) in 60 s. |
| **Assumptions and limits** | The base acceleration is a stationary, zero-mean Gaussian process with a known one-sided PSD; the component is a linear SDOF system with viscous damping (Q constant); the base is rigid and its motion does not depend on the component (no mechanical-impedance coupling, which is what force-limited testing corrects); one mode dominates. Miles' equation adds two more: the input PSD is flat around fn, and the resonance dominates the integral. Real boxes have several modes, amplitude-dependent damping and non-Gaussian peaks from nonlinear joints. The input spectrum is the NASA-STD-7001C Table 2 workmanship level; the box (2 kg, 150 Hz, Q = 10) is illustrative. All time histories here are SYNTHETIC. |
| **Analytical work** | dB/octave slopes and the exact area of a log-log segment; the base-excited SDOF and its absolute-acceleration transmissibility; the input–output PSD relation; the white-noise integral ∫|H|²df = (π/2) fn Q (1 + 4ζ²) and Miles' equation; the relative-displacement RMS; Gaussian exceedance, Rice's upcrossing rate from spectral moments and Davenport's peak factor; the variance of a time-domain RMS estimate, 1/√(ζωₙT). |
| **Computational deliverable** | `starter/main.py`: `segment_slopes_db_per_octave`, `psd_loglog`, `psd_area`, `grms`, `transmissibility_sq`, `miles_grms`, `response_grms`, `relative_displacement_rms`, `spectral_moments`, `upcrossing_rate`, `expected_peak_factor`, `synthesize_random`, `sdof_discrete`, `simulate_sdof_base` (`_log_grid`, `gaussian_exceedance`, `variance_estimate_relative_sd` and `main` are given). Tests (16) prove: exact segment areas against hand closed forms and adaptive quadrature; dB/octave slopes; √-scaling with level; transmissibility limits (1 at f = 0 and at r = √2 for every ζ, 1 + 1/(4ζ²) at resonance); the full integral reproduces the exact white-noise result and tends to Miles as ζ → 0; grid convergence; the relative displacement against its closed form and its fn^(−3/2) scaling; the upcrossing rate near fn; the synthesized record's variance, mean, Gaussian tails, reproducibility and band-by-band spectrum; the discretization against sinusoidal steady state, a static input and `scipy.signal.lsim`; and the analytic Grms against a seeded time-domain simulation within four derived standard deviations. |
| **Practical activity** | **computational experiment**: Miles against the full integral over natural frequency and damping; the scatter of a time-domain Grms against record length over 20 seeds, compared with the derived law; peaks per second against Davenport. All random data in this lesson are SYNTHETIC. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `random_vibration.png`, the hand derivation of the white-noise integral and of the relative-displacement RMS, the experiment tables, and a half-page note stating the design load for the box's mounts with the assumption you trust least. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the exact segment area and the white-noise integral are derived; all public tests pass; the analytic and time-domain Grms are compared with a stated statistical tolerance; the note gives a design load that accounts for test duration, not only "3σ". |
| **References** | `@nasa-std-7001c-payload-vibroacoustic-test-criteria` NASA-STD-7001C *Payload Vibroacoustic Test Criteria*: 4.1 Summary of Verification Test Requirements (Table 1); 4.3 Test Levels, 4.3.3 Workmanship (Table 2—Component Minimum Workmanship Random Vibration Test Levels); 4.4.2 Random Vibration Tests; 4.4.3 Test Duration; 4.4.4 Test Control Tolerances (required-free): the input spectrum, the +3 dB qualification margin and the test durations. `@mit-ocw-6-011-signals-systems-and-inference` Lecture 16: Wide-sense Stationary Processes, LTI Filtering of WSS Processes; Lecture 17: LTI Filtering of WSS Processes; Lecture 18: Power Spectral Density; Lecture 19: Einstein-Wiener-Khinchin Theorem, PSD Applications, Modeling Filters (required-free): S_yy = |H|²S_xx and the PSD. `@mit-ocw-18-03-differential-equations` notes 16 "Frequency response" (required-free): the gain of a second-order system. `@mit-ocw-6-003-signals-and-systems` Lecture 21: Sampling (supplementary-free refresher). Miles' equation, the 3σ convention and Davenport's peak factor are not in these sources; part 2 derives them. Optional textbooks for deeper study (paid, never required; try a library): Wirsching, Paez and Ortiz, *Random Vibrations: Theory and Practice*, for SDOF and MDOF random response and Miles' equation. |
| **Extensions** | ★ Give the box a second mode (a 2-DOF model, `P06.14`) and compare the response Grms with the root-sum-square of two Miles estimates. ★ Replace viscous damping by a constant Q that varies with frequency (Q = √fn is a common rule of thumb) and repeat the fn sweep. ★ Count stress cycles with the Gaussian bands (68.3 % of the time within 1σ, 27.1 % between 1σ and 2σ, 4.3 % between 2σ and 3σ) and estimate fatigue damage with a Miner's sum (`P05.14`). |
| **Time** | unmeasured |

## 1. The engineering problem

A launch vehicle shakes everything it carries. Rocket-engine noise, turbulent flow over the fairing and combustion roughness reach the payload through its structure as broadband vibration, from about 20 Hz to 2000 Hz, with no repeating waveform. It is not a sine wave with an amplitude you can read off; it is a random process with a spectrum. Every electronics box, star tracker and propellant valve on a spacecraft is therefore tested on a shaker driven with random noise shaped to a specified **acceleration power spectral density** (PSD): a curve in g²/Hz, given as a few breakpoints joined by straight lines on log-log axes. NASA-STD-7001C sets a *minimum workmanship* level (its 4.3.3 Workmanship and Table 2), a spectrum every flight component must pass to show that it was built properly (no loose screws, cracked solder joints or rubbing cables). Its Table 1 puts qualification levels 3 dB above flight acceptance, with minimum durations of 1 min for acceptance and protoflight tests and 2 min for qualification.

The engineer's job before the test is to predict what the box will feel. A box mounted on its brackets is a mass on a spring: at its natural frequency it amplifies the input by roughly Q, the quality factor, so a component that sees 6.8 Grms at its feet can see 10 Grms or more at its centre of mass. The mounts, the circuit boards and the connectors must carry that load, and the box must not hit its neighbours. Get the prediction wrong and either the box breaks on the shaker (a programme delay) or it is over-designed (mass that could have been payload).

This lesson builds that prediction three ways: Miles' equation, the one-line formula every structures engineer uses; the full integral of |H|²·G, which Miles approximates; and a time-domain simulation driven by a synthesized random signal, which checks both by an independent route. Success criteria, written first: the segment integrals are exact; the full integral reproduces the exact white-noise result; the simulated RMS agrees with the full integral within four standard deviations of its derived statistical scatter; and the design load is stated with the test duration in it.

## 2. Understand and derive

### 2.1 A spectrum in g²/Hz and its overall level

Let a(t) be the base acceleration in g. Its **one-sided PSD** G(f) (`P02.19`) is defined for f ≥ 0 so that the mean square is the area under it:

**a²_rms = ∫₀^∞ G(f) df**, and the overall level is **Grms = √(∫G df)**.

Units: (g²/Hz) × Hz = g². A test specification gives G at breakpoints and joins them by straight lines on log-log axes, so on each segment G is a power law, G = G₁(f/f₁)ⁿ. The slope is quoted in **dB per octave**: one octave is a doubling of frequency, and a change of 10 log₁₀(G₂/G₁) dB over log₂(f₂/f₁) octaves gives

m = 10 log₁₀(G₂/G₁) / log₂(f₂/f₁)  dB/oct,  n = m / (10 log₁₀ 2) = m/3.0103.

So +3 dB/oct is n = 1 (G proportional to f), and −6 dB/oct is n = −2. The area of one segment is exact:

**∫_{f₁}^{f₂} G df = G₁f₁[(f₂/f₁)^{n+1} − 1]/(n + 1),** and **G₁f₁ ln(f₂/f₁)** when n = −1.

The second case is not a curiosity: −3 dB/oct is a common roll-off, and the general formula divides by zero there. Integrating a log-log spectrum with the trapezoid rule on its breakpoints is a classic error (it joins the points by straight lines on *linear* axes and overestimates every sloped segment).

**Worked example: the workmanship spectrum.** The NASA-STD-7001C Table 2 spectrum (components under 50 kg) rises at +3 dB/oct from 0.01 g²/Hz at 20 Hz to 0.04 g²/Hz at 80 Hz, stays flat to 500 Hz, and falls at −3 dB/oct to 0.01 g²/Hz at 2000 Hz. The three areas are 0.01/20 × (80² − 20²)/2 = 1.500 g², 0.04 × 420 = 16.800 g², and 0.04 × 500 × ln 4 = 27.726 g². The total is 46.026 g², so **Grms = 6.784**, the 6.8 grms the standard states for it: a built-in check of the integration. 60 % of the energy sits above 500 Hz, where the spectrum is falling: the slope is not a detail.

### 2.2 The component as a base-excited SDOF system

Model the box as a mass m on a spring k and a viscous damper c whose other ends move with the base displacement y(t) (`P06.13`). Let x be the absolute displacement of the mass and z = x − y the relative displacement, the deflection of the mounts. Newton's second law for the mass, in an inertial frame:

m ẍ = −k z − c ż.

Write ẍ = z̈ + ÿ, divide by m, and use ωₙ² = k/m, 2ζωₙ = c/m:

**z̈ + 2ζωₙ ż + ωₙ² z = −ÿ,**  and the absolute acceleration is **ẍ = −(2ζωₙ ż + ωₙ² z).**

For a harmonic base acceleration ÿ = A e^{iωt}, z = Z e^{iωt} with Z = −A/(ωₙ² − ω² + 2iζωₙω), and ẍ = (ωₙ² + 2iζωₙω)Z·(−1). With r = ω/ωₙ, the ratio of absolute accelerations is the **transmissibility**

**|H(f)|² = [1 + (2ζr)²] / [(1 − r²)² + (2ζr)²].**

Read its limits. At r → 0 it is 1: a slow shake moves the box with the base. At r = 1 it is 1 + 1/(4ζ²) ≈ Q², with Q = 1/(2ζ). At r = √2 it is exactly 1 for every damping (put r² = 2: numerator and denominator are both 1 + 8ζ²); above that the box is isolated, which is lesson `P06.17`. At r ≫ 1 it falls as (2ζ/r)².

A linear, time-invariant system driven by a stationary random input has a stationary output whose PSD is the input PSD multiplied by the squared gain at each frequency (`P02.19`):

**S_out(f) = |H(f)|² G(f),  Grms_out² = ∫₀^∞ |H(f)|² G(f) df.**   (1)

That is the *full integral*. It needs only the transfer function and the spectrum; no time history.

### 2.3 Miles' equation, and what it leaves out

If G were constant (white) over the whole axis, (1) could be done by hand. Substitute f = r fₙ:

∫₀^∞ |H|² df = fₙ ∫₀^∞ (1 + 4ζ²r²) dr / [(1 − r²)² + 4ζ²r²].

Call I₀ = ∫₀^∞ dr/D(r) and I₂ = ∫₀^∞ r² dr/D(r), with D = (1 − r²)² + 4ζ²r². The substitution r → 1/r maps D(r) to D(1/r) = D(r)/r⁴ and dr to −dr/r², so I₂ = I₀: the two integrals are equal. A standard contour integral (or the table result for ∫dω/|ωₙ² − ω² + 2iζωₙω|², `P02.19`) gives I₀ = π/(4ζ). Hence

**∫₀^∞ |H|² df = fₙ π(1 + 4ζ²)/(4ζ) = (π/2) fₙ Q (1 + 4ζ²).**

For a white input of level G this is exact. Dropping the 4ζ² (it is 1 % at ζ = 0.05) and reading G at fₙ gives **Miles' equation**:

**Grms_out ≈ √( (π/2) fₙ Q G(fₙ) ).**   (2)

Miles therefore makes three approximations: the input is flat *across the whole band that matters*, the factor (1 + 4ζ²) is 1, and nothing outside the band is missing or extra (a real specification stops at 20 Hz and 2000 Hz, while the white-noise integral runs from 0 to ∞). Below resonance |H|² ≈ 1, so the input itself passes through; Miles ignores that *quasi-static* contribution, which matters when fₙ is near the top of the spectrum. For the box, fₙ = 150 Hz sits on the flat plateau: G(fₙ) = 0.04 g²/Hz, and (2) gives √(π/2 × 150 × 10 × 0.04) = **9.708 Grms**. The full integral (1) gives **9.658 Grms**. Miles is 0.5 % high here; part 4 maps where it is not.

### 2.4 How far the box moves

The mounts deflect by z. From 2.2, |Z/Ÿ|² = 1/[(ωₙ² − ω²)² + (2ζωₙω)²], so with G in g²/Hz the relative-displacement RMS is

z_rms² = g₀² ∫ G(f) / [(ωₙ² − ω²)² + (2ζωₙω)²] df,

with g₀ = 9.80665 m/s² per g. For white G the same integral I₀ gives ∫df/[…] = 1/(64π³ζfₙ³), so

**z_rms ≈ g₀ √( G(fₙ)/(64π³ζfₙ³) ) = g₀ Grms_Miles / ωₙ².**

For the box, 0.1061 mm RMS by the full integral (0.1072 mm by the Miles-type estimate), 0.318 mm at 3σ. A stiffer mount moves less: for a flat input, z_rms falls as fₙ^(−3/2) while the acceleration rises as fₙ^(1/2). That trade is the subject of `P06.17`.

### 2.5 Peaks: what "3σ" does and does not mean

A linear system driven by a Gaussian input has a Gaussian output (a weighted sum of Gaussians). Its standard deviation σ is the Grms of (1), and the fraction of time it spends beyond ±kσ is erfc(k/√2): 31.7 % beyond 1σ, 4.55 % beyond 2σ, **0.270 % beyond 3σ**. That is why loads are often quoted as "3σ": 29.0 g for the box, and 3σ × m g₀ = **568 N** through its mounts if the whole 2 kg moves as one mass.

But 0.27 % of a long test is a lot of time. The response of a lightly damped system is narrow-band: it oscillates near fₙ with a slowly varying envelope, and it crosses zero upward at Rice's rate

**ν₀⁺ = √(m₂/m₀),**  mₖ = ∫ fᵏ S_out(f) df,

which is 151.0 per second for the box, close to fₙ as it should be. In a test of duration T there are N = ν₀⁺T cycles, and Davenport's result for the expected largest peak of a stationary Gaussian narrow-band process is

**peak/σ ≈ √(2 ln N) + 0.5772/√(2 ln N).**

For 60 s, N = 9061 and the expected largest peak is **4.40σ = 42.5 g**, half again the 3σ value. A part that must not yield anywhere in the test should be checked at the expected largest peak (or a higher percentile), not at 3σ. A part limited by fatigue is the other way round: what matters is how many cycles it sees at each level, and most are below 2σ.

### 2.6 Making a random signal with a given PSD

To check (1) independently we need a time history whose PSD is G. Put N samples at rate f_s over a duration T = N/f_s; the FFT frequencies are f_k = k Δf with Δf = 1/T. Build

a(t_m) = Σ_k √(G(f_k)Δf) [α_k cos(2πf_k t_m) − β_k sin(2πf_k t_m)],  α_k, β_k independent standard normal.

Each term has variance G(f_k)Δf, the terms are independent, so the variance is Σ G(f_k)Δf, the Riemann sum of ∫G df, and the PSD is G at every bin. Because the coefficients are Gaussian, so is a(t). One inverse real FFT computes the whole sum (`np.fft.irfft` returns (1/N)[X₀ + 2ΣRe(X_k e^{2πikm/N})], so the result is scaled by N/2). The record is periodic with period T, which is harmless here. Using random *phases* with fixed amplitudes instead gives a record whose variance is exactly the target but whose values are not quite Gaussian; the version here is the honest model of a random process. Everything it produces is SYNTHETIC: a model of a shaker input, not a measurement.

### 2.7 Simulating the box exactly

Write the SDOF equation of 2.2 in state form, ẋ = A x + B u with x = [z, ż], u = ÿ, A = [[0, 1], [−ωₙ², −2ζωₙ]], B = [0, −1]ᵀ, output y = ẍ = C x with C = [−ωₙ², −2ζωₙ]. If the input varies linearly between samples (a **first-order hold**), the state update over one step is exact (`P02.05`):

x_{k+1} = Φ x_k + (Γ₁ − Γ₂) u_k + Γ₂ u_{k+1},

where Φ, Γ₁ and Γ₂ are blocks of the exponential of the augmented matrix [[AΔt, BΔt, 0], [0, 0, 1], [0, 0, 0]]. Eliminating the state gives a second-order recursive filter, which `scipy.signal.lfilter` runs in compiled code. With f_s = 8192 Hz the discretization is exact for piecewise-linear input and the remaining error is the hold itself, which matters only near the Nyquist frequency where |H|² is tiny.

### 2.8 How long a simulation must be

The simulated RMS is itself a random estimate. For the response of 2.2, the autocorrelation is close to σ² e^{−ζωₙ|τ|} cos ωₙτ, and the variance of a variance estimate over a record of length T is (2/T)∫R²dτ (`P02.19`). The integral is σ⁴/(2ζωₙ), so

**sd(σ̂²)/σ² ≈ 1/√(ζωₙT),**  and the RMS scatters by half that.

For the box (ζωₙ = 47.1 s⁻¹) a 2 s record scatters by 5.2 %, a 32 s record by 1.3 %. Lightly damped systems need long records, because their response "remembers" for about 1/(ζωₙ) seconds and each memory length is one independent sample.

### 2.9 Where the model fails

- **The base is not rigid.** On a shaker the fixture is very stiff and the input spectrum is imposed regardless of the box; in flight the panel is flexible and the box's own resonance pulls the panel's motion down near fₙ (a *notch* in the interface acceleration). Shaker tests therefore over-test at resonance; force-limited vibration testing corrects this by limiting the interface force.
- **One mode.** A box has many modes; each contributes its own Miles-type term, and the board-level modes (often 200–600 Hz) matter more for solder joints than the box's mounting mode.
- **Linear, viscous damping.** Q depends on amplitude in bolted joints, so a response measured at a low level overpredicts Q at full level; a rule of thumb Q ≈ √fₙ is common for electronics. Nonlinear joints also make the peaks non-Gaussian.
- **Stationarity.** Real launch vibration varies through flight; the test spectrum is an envelope.

## 3. Build the computational model

Open `starter/main.py`. The constants, `_log_grid` (an integration grid packed around fₙ), `gaussian_exceedance`, `variance_estimate_relative_sd` and `main` are given.

### Step 1: the input spectrum

`segment_slopes_db_per_octave` is the formula of 2.1. `psd_loglog(f, f_bp, G_bp)` interpolates log G against log f with `np.interp` and returns zero outside the band. `psd_area` loops over segments with the exact area, handling n = −1 separately. `grms` is its square root.

### Step 2: the SDOF response

`transmissibility_sq` and `miles_grms` are 2.2 and (2). `response_grms` evaluates |H|²G on `_log_grid(f_bp, fn, zeta, n_per_decade)` and integrates with `np.trapezoid`; the grid matters because at ζ = 0.01 the resonance is only 3 Hz wide. `relative_displacement_rms` is 2.4 on the same grid. `spectral_moments` returns m₀ and m₂, `upcrossing_rate` is √(m₂/m₀), and `expected_peak_factor` is Davenport's formula.

### Step 3: the time domain

`synthesize_random(f_bp, G_bp, fs, duration, seed)` builds the complex coefficients √(GΔf)(α + iβ) on the rfft grid with `np.random.default_rng(seed)`, zeroes the DC and Nyquist bins, and calls `np.fft.irfft(X, n) * n / 2`. `sdof_discrete(fn, zeta, dt)` builds the augmented matrix, takes `scipy.linalg.expm`, and converts the state update to filter coefficients (b, a) through the characteristic polynomial z² − tr(Φ) z + det Φ and the adjugate of zI − Φ. `simulate_sdof_base` calls `lfilter(b, a, a_base)`.

Running `main()` prints

```
input spectrum: NASA-STD-7001C Table 2 component minimum workmanship level:
  breakpoints [20.0, 80.0, 500.0, 2000.0] Hz; [0.01, 0.04, 0.04, 0.01] g^2/Hz; slopes [3.01, 0.0, -3.01] dB/oct
  overall 6.784 Grms (area 46.026 g^2)

avionics box: fn = 150.0 Hz, zeta = 0.05 (Q = 10), mass 2.0 kg
  Miles: 9.708 Grms; full integral: 9.658 Grms (Miles/full = 1.0052)
  3-sigma acceleration 28.97 g; 3-sigma mount force 568 N
  relative displacement 0.1061 mm rms, 3-sigma 0.318 mm (Miles-type estimate g Grms/wn^2 = 0.1072 mm)
  upcrossing rate 151.0 /s; expected largest peak in a 60 s test: 4.40 sigma = 42.5 g; fraction of time beyond 3 sigma 0.270 %

SYNTHETIC time-domain check (seed 16, 32 s at 8192 Hz):
  input 6.758 Grms (target 6.784); response 9.492 Grms (integral 9.658; expected scatter +/- 1.3 % at 1 sigma)
  fraction beyond 3 sigma 0.275 % (Gaussian 0.270 %); largest |peak| 4.25 sigma (Davenport for 31.5 s: 4.26)
```

followed by the Miles table of part 4, and writes `random_vibration.png`: the input and response PSDs with their Welch estimates from the SYNTHETIC record, a quarter-second of the SYNTHETIC time history with the analytic ±3σ lines, and Miles' accuracy against damping.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-16-random-vibration`. The 16 tests:

1. **Spectrum**: the workmanship area by hand (including the logarithmic −3 dB/oct case) to 10⁻¹²; an arbitrary five-breakpoint spectrum against adaptive quadrature to 10⁻⁹; interpolation through the breakpoints; ±3.01 dB/oct slopes; a +3 dB/oct segment doubling G per octave; Grms doubling when G quadruples.
2. **Transmissibility**: 1 at f = 0; exactly 1 at r = √2 for three dampings; 1 + 1/(4ζ²) at resonance; monotone fall above resonance, lower for less damping.
3. **Integrals**: the full integral against the exact white-noise result (π/2)fₙQG(1 + 4ζ²) minus the computed band tails, to 2 × 10⁻⁴; Miles within 1 % on the plateau and closer as ζ falls; grid convergence; the relative displacement against 1/(64π³ζfₙ³) and its fₙ^(−3/2) scaling; ν₀⁺ within 10 % of fₙ.
4. **Statistics**: exceedance fractions; Davenport's factor between 3 and 5 and growing with T.
5. **Synthesis**: variance within 4 % (four standard deviations of the ≈ 1 % bin scatter), zero mean, the Gaussian 2σ fraction, reproducibility by seed, and Welch band averages within 5 % of G.
6. **Discretization**: sinusoidal steady state against |H| at 40, 150 and 400 Hz within the sample-peak-picking error; a constant base acceleration passed unchanged; zero in, zero out; agreement with `scipy.signal.lsim`.
7. **The independent check** (part 5): the analytic Grms against the seeded simulation, and the simulated 3σ exceedance against 0.27 %.
8. **main** writes the figure and says SYNTHETIC.

### Compare

`scipy.signal.lsim` integrates the same state-space model with its own first-order-hold discretization. After the start-up transient (the two codes assume different input before t = 0, and the difference decays as e^{−ζωₙt}), they agree to **4.3 × 10⁻¹¹ g** on a response of 10 g RMS: the same discretization, implemented twice. `scipy.signal.welch` (`P02.19`) estimates the PSD of the SYNTHETIC response; averaged over 140–160 Hz it gives 2.860 g²/Hz against 2.885 g²/Hz from |H|²G, a 0.9 % difference consistent with the Welch estimate's scatter over 19 bins. The library tools agree with the code; the independent check of the *physics* is part 5.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All random records are SYNTHETIC. Write each prediction before running it.

1. **Where is Miles good enough?** Miles / full integral for fₙ from 30 to 1800 Hz and ζ = 0.01, 0.05, 0.10. Predict where Miles over- and under-estimates. *Finding:*

   | fₙ [Hz] | ζ = 0.01 | ζ = 0.05 | ζ = 0.10 |
   |---|---|---|---|
   | 30 | 1.0057 | 1.0243 | 1.0349 |
   | 60 | 1.0043 | 1.0173 | 1.0227 |
   | 150 | 1.0018 | 1.0052 | 1.0006 |
   | 400 | 1.0011 | 1.0014 | 0.9938 |
   | 1000 | 0.9936 | 0.9672 | 0.9330 |
   | 1800 | 0.9975 | 0.9827 | 0.9514 |

   On the rising +3 dB/oct segment Miles is high: G(fₙ) overstates the spectrum below fₙ, where half the resonance's area sits. On the falling segment Miles is low, by up to 6.7 % at ζ = 0.1 and 1000 Hz, because the quasi-static part (|H|² ≈ 1 below resonance, over the whole energetic plateau) is not in Miles at all. Heavier damping widens the resonance and makes both effects bigger. On the plateau Miles is within 0.6 %.

2. **Where should the box's frequency be?** Response Grms and z_rms against fₙ at ζ = 0.05. Predict the shape. *Finding:* 2.60, 4.36, 6.85, 9.66, 13.72, 17.47, 18.18, 18.42, 18.04 Grms at 30, 50, 80, 150, 300, 500, 800, 1200, 1800 Hz; z_rms 0.712, 0.431, 0.265, 0.106, 0.038, 0.017, 0.007, 0.003, 0.001 mm. Acceleration rises as √(fₙG(fₙ)) and flattens on the −3 dB/oct segment (where fₙG is constant); displacement falls steeply everywhere. A stiff box is not a quiet box: it trades displacement for acceleration.

3. **Damping.** At fₙ = 150 Hz, ζ = 0.01–0.2. *Finding:* 21.67, 15.30, 9.66, 6.86, 5.03 Grms at ζ = 0.01, 0.02, 0.05, 0.10, 0.20 (Q = 50 to 2.5), following Q^(1/2); Miles within 0.5 % up to ζ = 0.1 and 3.5 % low at ζ = 0.2. Damping pays as its square root: halving the response needs four times the damping.

4. **How long must a simulation run?** 20 seeds at each of T = 2, 8, 32 s; RMS / full integral. Predict the scatter from 2.8. *Finding:* mean 1.0009, 0.9942, 0.9953; standard deviation 4.57 %, 2.38 %, 1.64 % against the predicted 5.15 %, 2.58 %, 1.29 %. The law holds to within the scatter of 20-sample standard deviations (which is itself about 16 %). The extremes at T = 2 s were 0.929 and 1.139: a single short simulation can be 14 % off and look like a disagreement with theory when it is only noise.

5. **Peaks.** Split the 32 s SYNTHETIC response (seed 16) into 1 s windows; largest |a|/σ in each. *Finding:* mean 3.38 (standard deviation 0.37) against Davenport's 3.35 for 1 s; over the whole 32 s, 4.25σ against Davenport's 4.26. Fractions beyond 1σ, 2σ, 3σ: 31.9 %, 4.47 %, 0.278 % against 31.7 %, 4.55 %, 0.270 % for a Gaussian. For a 60 s test the expected largest peak is 4.40σ; for 180 s, 4.65σ.

Raw data format: CSV tables `fn_Hz, ratio_z001, ratio_z005, ratio_z010` (1); `fn_Hz, G_fn, grms_miles, grms_full, z_rms_mm` (2); `zeta, Q, grms_miles, grms_full` (3); `T_s, seed, rms_ratio` (4); `window, peak_over_sigma` (5); all labelled SYNTHETIC where a time history was used.

## 5. Validate and improve

The independent check: the analytic Grms of (1) against a time-domain simulation of a box driven by a synthesized random signal. They share only the PSD and the SDOF parameters.

| quantity | analytic | SYNTHETIC simulation (seed 16, 32 s) | difference | expected scatter (1σ) |
|---|---|---|---|---|
| input Grms | 6.784 (exact area) | 6.758 | −0.4 % | 0.22 % (half the 0.44 % bin scatter of the variance) |
| response Grms | 9.658 (full integral) | 9.492 | −1.7 % | 1.3 % (2.8) |
| response Grms, mean of 20 seeds × 32 s | 9.658 | ratio 0.9953 | −0.47 % | 0.37 % (1.64 %/√20) |
| fraction beyond 3σ | 0.270 % | 0.275 % | +0.005 % | — |
| largest peak in 31.5 s | 4.26σ (Davenport) | 4.25σ | −0.01σ | ≈ 0.31σ (Gumbel spread π/√(12 ln N)) |
| PSD at 140–160 Hz | 2.885 g²/Hz | 2.860 g²/Hz (Welch) | −0.9 % | ≈ 1 % |
| Miles' equation | 9.708 | — | +0.5 % vs full | approximation |

Every difference is inside two standard deviations of its predicted scatter (the input is 1.7, the response 1.3). The single 32 s run is 1.7 % low, which would be alarming without the scatter estimate and is unremarkable with it. The mean of 20 seeds is closer, as it must be. Miles' equation, an approximation, is the largest systematic difference in the table, and still only 0.5 %.

The engineering answer: **the box responds at 9.66 Grms; its mounts see 29.0 g (568 N) at 3σ and an expected largest peak of 42.5 g (4.40σ) in a 60 s test; it moves 0.11 mm RMS on its mounts.** For a strength check use the expected largest peak (or a 99th-percentile peak, which is higher); for clearance use 3σ to 5σ of z; for fatigue count cycles by level. Miles is good to 0.6 % for this box because 150 Hz sits on the flat plateau; a box with fₙ near 1000 Hz on this spectrum needs the full integral (Miles is 3.3 % low there at ζ = 0.05 and 6.7 % low at ζ = 0.1).

What to improve: the model, not the arithmetic. The largest uncertainty is Q (an estimate of 10 could be 5 or 20 in a bolted box, changing the response by a factor √2 either way), then the extra modes, then the rigid-base assumption. A low-level sine sweep on the shaker before the random test measures fₙ and Q; re-running (1) with the measured values is the improvement that matters.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of the segment area, of I₂ = I₀ and of the white-noise integral, and of z_rms; `random_vibration.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the mount design load with its duration and the assumption you trust least.

### Adversarial self-checks

- **What if the input level doubles (+3 dB, the qualification margin of NASA-STD-7001C Table 1)?** Every RMS rises by √2: 6.78 → 9.59 Grms input, 9.66 → 13.66 Grms response. Linear theory says nothing new; a real box with amplitude-dependent damping responds a little less than √2 more.
- **What if the damping goes to zero?** Miles and the full integral both grow as 1/√ζ without bound, and the simulation needs ever longer records (2.8): T must exceed many times 1/(ζωₙ). An undamped resonance on a broadband spectrum has infinite RMS response, which is why no real test meets one.
- **Which assumption fails first?** The rigid base: in flight, the panel's motion dips at the box's resonance and the true response is lower than the shaker test's. Then Q, then single-mode behaviour.
- **What measurement could falsify the model?** A low-level random run on a shaker with an accelerometer on the box: the measured response PSD divided by the input PSD should be |H|² with the fitted fₙ and ζ. A response that grows less than linearly with input level falsifies linearity; a kurtosis well above 3 falsifies the Gaussian peak statistics.
- **Does the result survive a change of method?** Exact segment area, quadrature, the full integral, Miles and the seeded simulation agree to within their stated errors; quadrupling the integration grid changes the full integral by about 10⁻⁶.
- **Which decision changes if the uncertainty doubles?** If Q is 10 within a factor of 2 instead of √2, the 3σ load spans 21–41 g rather than 24–34 g, and a mount sized for 3σ at Q = 10 has no margin left at Q = 20. A low-level sine survey to measure Q becomes mandatory.

## Check yourself

1. What are the units of an acceleration PSD, and what is the area under it?
2. A spectrum falls at −6 dB/oct. What power of f is that, and by what factor does G change over one octave?
3. Why does integrating a log-log spectrum with the trapezoid rule on its breakpoints overestimate Grms?
4. State Miles' equation and its three assumptions.
5. Where on a +3 dB/oct / flat / −3 dB/oct spectrum is Miles' equation high, and where is it low? Why?
6. What fraction of the time does a Gaussian response spend beyond 3σ, and why is 3σ not the largest peak in a 60 s test?
7. A simulated RMS differs from the analytic Grms by 2 %. What else must you know before calling it a disagreement?
8. Doubling a box's natural frequency on a flat spectrum changes its response acceleration and its relative displacement by what factors?
