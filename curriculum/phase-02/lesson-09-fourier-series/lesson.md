# P02.09  Fourier Series

> Any periodic signal is a sum of sinusoids, and a linear structure answers each sinusoid separately. Find the harmonic that lands on a resonance and you have found the vibration.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Fourier coefficient formulas from orthogonality, the coefficients of the square and triangle waves, and Parseval's identity, and verify each by quadrature and with SymPy. 2. Implement and test, from scratch, a trapezoid-rule coefficient calculator, partial sums, the complex form, the Gibbs overshoot and the rms truncation error, and show that the convergence rate is set by the smoothness of the waveform. 3. Compute the steady-state response of a lightly damped structural mode to a periodic thruster pulse train harmonic by harmonic, and check it against `solve_ivp`. 4. Decide which pulse periods a reaction-control schedule must avoid so that no harmonic of the pulse train lands on a solar-array mode. |
| **Prerequisites** | `P02.03` (forced response, superposition), `P02.07` (phasors and the complex gain), `P01.04` (integration of products of trigonometric functions), `P01.05` (series and convergence), `P00.09` (SymPy). Skill check: you can integrate sin(mx) sin(nx) over a period, sum a geometric series, and state what superposition means for a linear ODE. |
| **Engineering question** | A spacecraft's reaction-control thrusters fire in a 1 s on-off pattern (a square-wave torque); the solar-array boom has a bending mode at 5.0 Hz with damping ratio 0.01: how much does the pulse train excite the boom compared with a smooth sinusoidal torque of the same amplitude, which harmonic is responsible, and which pulse periods must the schedule avoid? |
| **Assumptions and limits** | The torque is exactly periodic with a sharp square waveform; the boom mode is linear with unit modal mass and viscous damping; the response is the steady state after the transient e^{−ζω_n t} has died (about 150 s here); every harmonic is treated independently (superposition), which needs linearity. The Fourier series of a signal with a jump converges in the mean square but not uniformly: the Gibbs overshoot is a property of the partial sum, not of the thruster. The method says nothing about non-periodic signals (P02.10) or about a mode that shifts frequency with amplitude. |
| **Analytical work** | Orthogonality of sines and cosines on a period; the coefficient formulas; the coefficients of the square wave (4/(πn), odd n) and the triangle wave (8/(π²n²), odd n); even and odd symmetry; the complex form c_n; Parseval's identity and the mean square; the rms truncation error from the tail of Parseval's sum and the N^{−1/2} and N^{−3/2} convergence rates; the Gibbs constant (2/π) Si(π); the steady-state response to a periodic force as the harmonic sum with the complex gain. |
| **Computational deliverable** | `starter/main.py`: `square_wave`, `triangle_wave`, `square_wave_coefficients`, `triangle_wave_coefficients`, `fourier_coefficients`, `partial_sum`, `complex_coefficients`, `from_complex_coefficients`, `gibbs_overshoot`, `parseval_power`, `mean_square`, `truncation_error_rms`, `convergence_rate`, `complex_gain`, `response_coefficients`, `harmonic_amplitudes`, `dominant_harmonic`, `forced_mode_rhs`, `sympy_square_wave_coefficient`. Tests prove: the waveforms have the stated symmetry and period; the closed-form coefficients match `quad` to 10⁻¹⁰ and SymPy exactly; the trapezoid coefficients converge to them and are exact to rounding for a band-limited signal; the partial sum reproduces a finite series and its rms error falls monotonically and matches Parseval's tail; the convergence exponents are 0.5 (jump) and 1.5 (kink); the Gibbs overshoot tends to 1.17898 from above and never shrinks; Parseval holds for every waveform; the complex form round-trips with conjugate symmetry; the response coefficients have the static and resonant limits and the 90° lag at resonance; the harmonic sum matches `solve_ivp` to 0.1 % in steady state. |
| **Practical activity** | **computational experiment**: a sweep of the boom response against pulse period (which periods resonate), of the Gibbs overshoot against the number of terms, and of the truncation error against smoothness. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `fourier_series.png`, the hand derivation of the square-wave coefficients and of Parseval's identity, the pulse-period sweep table, and a short note stating the pulse periods to avoid and the margin. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the coefficient formulas are derived from orthogonality; all public tests pass; the hand value of the fifth-harmonic response agrees with the code to three figures; the note names the forbidden pulse periods with the harmonic responsible for each. |
| **References** | `@mit-ocw-18-03-differential-equations` L15 "Introduction to Fourier Series", L16 "Continuation: More General Periods", L17 "Finding Particular Solutions via Fourier Series", notes 20 "Fourier series", 21 "Operations on Fourier series", 22 "Periodic solutions; resonance" (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 24 "Markov matrices; fourier series" (supplementary-free, the second half: Fourier series as projection onto an orthogonal basis). |
| **Extensions** | ★ Replace the square wave by a trapezoidal pulse with 50 ms rise and fall times and show the coefficients fall as 1/n² beyond n ≈ T/(π × rise time), which is why real thruster valves excite high modes less than the model says. ★ Compute the response of a mode whose frequency is *between* harmonics as the pulse duty cycle varies from 10 % to 90 %. ★ Implement Fejér (Cesàro) averaging of the partial sums and show it removes the Gibbs overshoot at the cost of a slower approach to the jump. |
| **Time** | unmeasured |

## 1. The engineering problem

A spacecraft holds its attitude with small thrusters that are either on or off. A typical limit-cycle schedule fires them in a regular pattern: for the case here, a torque of +τ₀ for half a second, −τ₀ for the next half second, repeated, a square wave with period T = 1 s. The solar arrays are long, light and lightly damped; the first bending mode of a boom is at f = 5.0 Hz with damping ratio ζ = 0.01, so its resonant magnification is 1/(2ζ) = 50. A 1 Hz forcing seems far from a 5 Hz mode. It is not: a square wave is not one sinusoid but a sum of them at 1, 3, 5, 7, … Hz with amplitudes 4/(πn), and the fifth harmonic lands exactly on the boom's mode.

The design question has three parts. How much does the pulse train excite the boom compared with a pure 1 Hz sinusoid of the same amplitude? Which harmonic is responsible, and what fraction of the response does it carry? And, since the pulse period is a free parameter of the control software, which periods must be avoided? Each is answered by the same two steps: decompose the periodic torque into its harmonics (the **Fourier series**, `@mit-ocw-18-03-differential-equations` notes 20) and, because the boom is linear, add up the response to each harmonic computed with the complex gain of P02.07 (notes 22). The second step is P02.03's superposition principle used fifty times at once.

Success criteria, written down first: our trapezoid-rule coefficients agree with the closed forms to 10⁻⁵ at 2048 samples and with `scipy.integrate.quad` and SymPy where those apply; the Gibbs overshoot of the partial sum is reproduced to 3 × 10⁻⁴ of its limit 1.17898; the rms truncation error falls as N^{−0.5 ± 0.06} for the square wave and N^{−1.5 ± 0.12} for the triangle wave; and the harmonic-sum response of the boom agrees with a direct integration of the modal equation to 0.1 % of its peak.

## 2. Understand and derive

### 2.1 The series and the coefficient formulas

A function f with period T can be written, with ω₀ = 2π/T,

**f(t) = a₀/2 + Σ_{n≥1} [a_n cos(nω₀t) + b_n sin(nω₀t)]**.   (1)

The basis functions are **orthogonal** over a period: ∫₀^T cos(mω₀t) cos(nω₀t) dt = (T/2) δ_mn for m, n ≥ 1, the same for sines, and every sine–cosine product integrates to zero (write the products as sums of cosines with the P01.04 identities and integrate whole cycles). Multiply (1) by cos(nω₀t), integrate over a period, and every term but one vanishes:

**a_n = (2/T) ∫₀^T f(t) cos(nω₀t) dt,   b_n = (2/T) ∫₀^T f(t) sin(nω₀t) dt,   a₀ = (2/T) ∫₀^T f dt**,   (2)

so a₀/2 is the mean of f. This is exactly the projection of a vector onto an orthogonal basis from P01.14 and `@mit-ocw-18-06-linear-algebra` Lecture 24, with ∫ f g dt as the inner product and the period as the space. Units: a_n and b_n have the units of f. Two symmetries save half the work: an **odd** function (f(−t) = −f(t)) has only sines, an **even** one only cosines.

### 2.2 The square and triangle waves

The odd square wave of amplitude A, +A on (0, T/2) and −A on (T/2, T), has a_n = 0 and

b_n = (2A/T)[∫₀^{T/2} sin(nω₀t) dt − ∫_{T/2}^T sin(nω₀t) dt] = (2A/(nπ))(1 − cos nπ) = **4A/(πn) for odd n, 0 for even n**.   (3)

SymPy's `integrate` returns 4/π, 0, 4/(3π), … for n = 1, 2, 3, … and the test checks it. The even triangle wave (+A at t = 0, −A at t = T/2, straight lines between) has b_n = 0 and, by two integrations by parts,

**a_n = 8A/(π²n²) for odd n, 0 for even n**.   (4)

The pattern to remember: a **jump** gives coefficients that fall as 1/n; a **kink** (a jump in the derivative) gives 1/n²; each further derivative that is continuous buys another power of n. The decay rate of the coefficients is the smoothness of the signal, and it will decide the convergence rate in 2.4.

```
   square (jump)                 triangle (kink)
    +A ┌────┐    ┌────          +A  /\        /\
       │    │    │                 /  \      /  \
    ───┼────┼────┼────  t        ─/────\────/────\─  t
       │    │    │                       \  /
    -A ┘    └────┘               -A       \/
   b_n = 4A/(pi n), odd n        a_n = 8A/(pi^2 n^2), odd n
```

### 2.3 Parseval's identity: power lives in the coefficients

Square (1) and integrate over a period; orthogonality kills every cross term:

**(1/T) ∫₀^T f² dt = a₀²/4 + ½ Σ_{n≥1} (a_n² + b_n²)**.   (5)

The mean square of the signal is the sum of the mean squares of its harmonics (each sinusoid of amplitude c has mean square c²/2). For the unit square wave the left side is 1, so (5) says (8/π²)(1 + 1/9 + 1/25 + …) = 1, i.e. Σ_{odd} 1/n² = π²/8, a series P01.05 could not sum. The code's `parseval_power` and `mean_square` check (5) for every waveform, and the tail of (5) gives the truncation error for free.

### 2.4 Convergence: how many terms, and the Gibbs overshoot

Let S_N be the sum of the first N harmonics. By (5) applied to f − S_N, the rms truncation error is

**rms(f − S_N) = √(½ Σ_{n>N} (a_n² + b_n²))**.   (6)

For the square wave the tail is (8/π²) Σ_{odd n>N} 1/n² ≈ 4/(π²N), so the rms error falls as **2/(π√N)**: N^{−1/2}, slow. The code measures exponents 0.48 for the square wave and 1.42 for the triangle wave, whose tail (32/π⁴) Σ 1/n⁴ ≈ 32/(3π⁴N³) gives N^{−3/2}. To halve the error of the square-wave series you need four times the terms; for the triangle wave, 1.6 times.

The mean-square convergence hides something. Near the jump the partial sum of the square wave overshoots by a fixed amount, **8.95 % of the jump height**, no matter how many terms are kept: the peak of S_N on (0, T/2) tends to (2/π) Si(π) = 1.17898 for a unit square wave (jump 2), approaching from above (1.2004 at N = 3, 1.1823 at N = 9, 1.1794 at N = 29, 1.1790 at N = 99) and moving ever closer to the jump, at t = T/(2(N + 1)). This is the **Gibbs phenomenon** (problem 6 derives the constant). The overshoot's *width* shrinks as 1/N so its energy vanishes and (6) is unaffected, but its *height* does not. Engineering consequence: a Fourier synthesis of a square wave always shows ears near the edges; they are a property of truncating the series, not of the thruster, and a plot that shows them is not showing a measurement.

### 2.5 The complex form

With Euler's formula (P02.07), cos(nω₀t) = (e^{inω₀t} + e^{−inω₀t})/2 and sin likewise, (1) becomes

**f(t) = Σ_{n=−∞}^{∞} c_n e^{inω₀t}**,   c₀ = a₀/2,   c_n = (a_n − i b_n)/2,   c_{−n} = c̄_n.   (7)

One formula replaces two, c_n = (1/T)∫₀^T f e^{−inω₀t} dt, and the harmonic n carries the phasor 2c_n = a_n − i b_n in the sense of P02.07 (Re(2c_n e^{inω₀t}) = a_n cos + b_n sin). The conjugate symmetry is what keeps f real. This form is the one that becomes the Fourier *transform* in P02.10 and the DFT in the same lesson.

### 2.6 Response of a mode to a periodic force

The boom mode obeys q̈ + 2ζω_n q̇ + ω_n² q = f(t), with q the modal coordinate per unit modal mass; for a unit constant force its static deflection is 1/ω_n². Each harmonic of f is a sinusoid at ω = nω₀, so by P02.07 its steady-state response is the input phasor times the complex gain,

**Q_n = H(r_n) (a_n − i b_n)/ω_n²,   H(r) = 1/(1 − r² + 2iζr),   r_n = nω₀/ω_n**,   (8)

and by superposition (the equation is linear, `@mit-ocw-18-03-differential-equations` notes 22) the steady-state response is the series with coefficients A_n = Re Q_n, B_n = −Im Q_n and A₀ = a₀/ω_n². For the thruster case ω₀ = 2π rad/s and ω_n = 10π rad/s, so r_n = n/5 and the fifth harmonic sits at r = 1 with |H| = 50:

| harmonic n | forcing amplitude 4/(πn) | r_n | |H| | response amplitude (static deflections) |
|---|---|---|---|---|
| 1 | 1.2732 | 0.2 | 1.042 | 1.326 |
| 3 | 0.4244 | 0.6 | 1.562 | 0.663 |
| 5 | 0.2546 | 1.0 | 50 | **12.732** |
| 7 | 0.1819 | 1.4 | 1.041 | 0.189 |
| 9 | 0.1415 | 1.8 | 0.446 | 0.063 |

The fifth harmonic, only a fifth as strong as the fundamental in the forcing, produces 9.6 times the response. The peak of the summed response is 14.33 static deflections, against 1.04 for a unit sinusoid at 1 Hz: **the pulse train shakes the boom 14 times harder than a smooth torque of the same amplitude would**, and 89 % of that is one harmonic. At exact resonance the response of that harmonic lags the forcing by 90°, so the sin(5ω₀t) forcing produces a −cos(5ω₀t) response; the test checks the sign.

This is the general rule for periodic forcing: **resonance occurs whenever any harmonic of the forcing, nω₀, coincides with ω_n**, not only when ω₀ does. For a square wave the odd harmonics are the dangerous ones; for a pulse train with another duty cycle (problem 1), every harmonic is present. Part 4 sweeps the pulse period and finds the forbidden values T = n/f_mode = 0.2, 0.6, 1.0, 1.4 s (odd n) for this boom.

### 2.7 Where the model stops being valid

- **Periodicity.** A thruster schedule that changes every few seconds is not periodic; its spectrum is continuous and P02.10's transform is the tool. The series still applies piecewise if the pattern persists for many periods compared with 1/(ζω_n) = 3.2 s here.
- **Sharp edges.** A real valve opens in 20–50 ms; the coefficients then fall as 1/n² beyond n ≈ T/(π t_rise), and the model's 1/n overestimates the high harmonics (extension).
- **Linearity.** Superposition needs it. A boom with a loose joint or large-amplitude stiffening answers the harmonics *jointly*, and the fifth-harmonic resonance can shift off 5 Hz as the amplitude grows.
- **Steady state.** With ζ = 0.01 the transient lasts about 5/(ζω_n) = 160 s; a pattern held for ten seconds never reaches the amplitudes of the table, which is itself a design option (change the pattern often).
- **Damping known.** The resonant response is 1/(2ζ); if ζ is 0.005 rather than 0.01 the fifth-harmonic response doubles. Nothing in the series is as uncertain as ζ.

## 3. Build the computational model

Open `starter/main.py`. `thruster_boom_case()` holds T, f_mode and ζ. Coefficient arrays are indexed from harmonic 1; a₀ is returned separately.

### Step 1: coefficients by the trapezoid rule

```python
def fourier_coefficients(f, T, n_max, n_samples=4096):
    if n_samples <= 2 * n_max: raise ValueError("aliasing")
    t = (np.arange(n_samples) + 0.5) * T / n_samples        # midpoints: never sample a jump exactly
    w0 = 2 * np.pi / T
    ft = f(t)
    a0 = 2 * np.mean(ft)
    n = np.arange(1, n_max + 1)
    an = 2 * (np.cos(np.outer(n, w0 * t)) @ ft) / n_samples
    bn = 2 * (np.sin(np.outer(n, w0 * t)) @ ft) / n_samples
    return a0, an, bn
```

For a periodic integrand the trapezoid (here midpoint) rule on a whole period is remarkable: for a band-limited f it is *exact* once n_samples exceeds twice the highest harmonic present, which the test checks to 10⁻¹³ with 64 samples, and for a signal with jumps its error falls as 1/n_samples² with midpoint sampling. The guard against n_samples ≤ 2 n_max is the sampling theorem in its simplest form: with fewer than two samples per cycle, harmonic n_max is indistinguishable from a lower one (**aliasing**, P02.10).

### Step 2: partial sums, the complex form, Parseval

`partial_sum` evaluates (1) up to N terms with one matrix product; `complex_coefficients` builds (7) as an array c_{−N..N} and `from_complex_coefficients` inverts it; `parseval_power` is the right side of (5) and `mean_square` the left, by the midpoint rule; `truncation_error_rms` is the rms of f − S_N on a fine grid and `convergence_rate` fits log(error) against log(N).

### Step 3: Gibbs

`gibbs_overshoot(N)` evaluates the N-term square-wave sum on a grid of 200 points per half-cycle of the highest harmonic over (0, T/4) and returns its maximum. The grid density matters: the overshoot peak is only T/(2(N + 1)) from the jump and a coarse grid misses it.

### Step 4: the response

```python
def response_coefficients(a0, an, bn, T, wn, zeta):
    n = np.arange(1, an.size + 1)
    r = n * (2 * np.pi / T) / wn
    out = complex_gain(r, zeta) * (an - 1j * bn) / wn**2             # (8)
    return a0 / wn**2, out.real, -out.imag
```

`harmonic_amplitudes` and `dominant_harmonic` read the table of 2.6; `forced_mode_rhs` is the modal equation for `solve_ivp`.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-09-fourier-series`. The tests:

1. **Analytic cases.** Square and triangle waves have the stated period, parity and values; (3) and (4) match `quad` to 10⁻¹⁰ and SymPy exactly; the trapezoid coefficients converge to (3) at better than second order and are exact to 10⁻¹³ for a three-term signal; the aliasing guard raises.
2. **Series algebra.** A finite series is reproduced exactly; the rms error falls monotonically with N and matches the Parseval tail to 2 %; the complex form round-trips, has conjugate symmetry and sums to the real series.
3. **Convergence and Gibbs.** Exponents 0.5 ± 0.06 and 1.5 ± 0.12; the overshoot decreases toward 1.17898 (3 × 10⁻⁴ at N = 199) and never below 1.17; one term gives 4/π.
4. **Conservation.** Parseval for the square wave (mean square 1), the triangle wave (1/3) and a trapezoid-computed signal (to 10⁻¹²).
5. **Response.** Amplitudes equal forcing × |H| to 10⁻¹²; harmonic 5 dominates with gain 1/(2ζ) and a 90° lag; a stiff mode responds statically; the harmonic sum matches `solve_ivp` after 40 s to 0.1 % of the peak; the square wave excites the mode more than 13 times a unit sinusoid.

### Compare

`scipy.integrate.quad` with a breakpoint at the jump is the oracle for the coefficients and agrees with (3) to 10⁻¹⁰; our midpoint rule gets 10⁻⁵ at 2048 samples and 10⁻⁷ at 8192, and for a smooth signal it beats `quad` outright (exact to rounding). `numpy.fft.rfft` would compute the same sums in O(N log N) and is what you use for a sampled signal (P02.10); we wrote the matrix version to see that a Fourier coefficient is an inner product and nothing more. `solve_ivp` on the modal equation reproduces the harmonic-sum response to 2.4 × 10⁻³ static deflections out of 14.3 after 35 s; the residual is the solver locating the jumps to within its 2 ms step and the series truncated at n = 51. The direct integration costs 20 000 steps and a wait for the transient; the harmonic sum is fifty complex divisions, and it also tells you *which* harmonic to blame.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three sweeps with `main.py`; all records are synthetic.

1. **Pulse-period sweep.** For T from 0.5 to 1.5 s in 0.01 s steps, the peak steady-state response of the 5 Hz boom (static deflections) to the unit square wave. **Prediction to write down first:** sharp peaks at T = 0.6, 1.0 and 1.4 s (harmonics 3, 5, 7 on the mode), with heights falling as 1/n, and a floor near 2 between them.
2. **Gibbs sweep.** The overshoot for N = 1, 3, 5, …, 199 and the position of the peak; predict the limit and the 1/(N + 1) law for the position.
3. **Smoothness sweep.** rms truncation error against N for the square wave, the triangle wave, and the square wave passed through a 50 ms rise-time ramp (make the waveform yourself); predict exponents 0.5, 1.5 and a transition from 0.5 to 1.5 near N ≈ T/(π × 0.05) = 6.

Raw data format: three CSV tables, `T, peak_response, dominant_harmonic`; `N, overshoot, t_peak`; `N, err_square, err_triangle, err_ramped`, each with a header naming the fixed parameters.

What you should find: resonance peaks of 14.3, 22.6 and 10.8 static deflections at T = 1.0, 0.6 and 1.4 s (the third harmonic is strongest because its forcing coefficient is larger), a floor of about 2.0–2.5 between them; the overshoot 1.2732, 1.2004, … → 1.1790 with t_peak = T/(2(N + 1)) to four figures; the ramped square wave's error following the square wave's for N below 6 and the triangle wave's slope beyond.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| b₁, b₃, b₅ of the square wave | 4/π, 4/(3π), 4/(5π) = 1.27324, 0.42441, 0.25465 | trapezoid 4096: 1.273240, 0.424414, 0.254649 | 10⁻⁶ | midpoint rule at a jump, O(1/n_samples²) |
| Parseval, 99 harmonics | 1 − 4/(π²·99)... tail ≈ 0.0041 → 0.9959 | 0.99595 | 10⁻⁵ | tail estimate is an integral approximation |
| rms error, N = 81 | 2/(π√81) = 0.0707 | 0.0703 | 0.6 % | same |
| Gibbs limit | (2/π) Si(π) = 1.17898 | 1.17901 at N = 99 | 3 × 10⁻⁵ | finite N: approaches as 1/N² |
| 5th-harmonic response | 0.25465 × 50 = 12.732 | 12.7324 | 0 | — |
| peak response | ≈ 12.73 + 1.33 + … ≈ 14.3 | 14.332 | — | phases of the other harmonics |
| harmonic sum vs `solve_ivp` | — | 2.4 × 10⁻³ of 14.3 | 0.02 % | solver step at the jumps; n ≤ 51 |

Every number is reproduced and every difference has a named, small cause. What the agreement does not show: that the boom's response is really 14 static deflections. The resonant term is 1/(2ζ), and a measured ζ of a deployed array has an uncertainty of a factor two either way; the peak response could be 7 or 28. The series has done its job by pointing at the fifth harmonic; the number it attaches to it is only as good as ζ, which is the quantity a ground vibration test or an in-orbit ring-down measures.

Improvements worth making: a trapezoidal pulse with the valve's real rise time (extension), which lowers the high harmonics and makes the 1/n² decay honest; sweeping duty cycle as well as period; reporting response with a ζ band instead of a single value. Improvement worth not making: adding harmonics beyond 51 to "fix" the 0.02 % disagreement with `solve_ivp`; the next harmonics change the peak by 10⁻⁴ and the uncertainty in ζ changes it by a factor two.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (2) from orthogonality and of (3) and (5); `fourier_series.png`; the three CSV tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the forbidden pulse periods for the 5 Hz boom, the harmonic responsible for each, the response at the chosen period, and the assumption (ζ) that limits the number's trustworthiness.

### Adversarial self-checks

- **What if the torque amplitude doubles?** Every coefficient and every response doubles; the dominant harmonic and the forbidden periods do not move. What if ζ halves? The fifth-harmonic response doubles to 25.5; everything else barely changes.
- **What if T → 0 or T → ∞?** Short periods push every harmonic above the mode (the boom sees the mean torque, a₀/2 = 0 here, plus small ripples); long periods put many harmonics below the mode and the response tends to the quasi-static square wave, amplitude 1 static deflection, plus a resonant burst at each edge.
- **Which assumption fails first in the real world?** The sharp edge: a 30 ms valve already halves the 15th harmonic. Second: the value of ζ.
- **What measurement could falsify the model?** A strain gauge on the boom during a 1 s pulse pattern: the model predicts a response at 5 Hz 9.6 times the 1 Hz component. A measured ratio of 3 means ζ is three times larger than assumed; a response at 5 Hz with a 1.2 s pattern (where no harmonic lands on the mode) means the boom is nonlinear or the mode is not at 5 Hz.
- **Does the result survive a change in the number of terms, samples or solver?** The coefficients change by 10⁻⁶ between 2048 and 4096 samples; the response peak by 10⁻⁴ between 51 and 101 harmonics; `solve_ivp` at rtol 10⁻⁹ and 10⁻⁶ by 10⁻⁴. The Gibbs overshoot does *not* go away with more terms, and that is the point.
- **Which decision changes if the uncertainty doubles?** If the mode were at 5.0 ± 0.5 Hz rather than ± 0.1 Hz, the forbidden periods would become bands (0.9–1.1 s instead of 0.98–1.02 s) and the schedule's free choice of T would shrink to the gaps; at ± 1 Hz the 1 s pattern would have to be abandoned in favour of patterns changed faster than the 160 s transient.

## Check yourself

1. Write the coefficient formulas (2) and explain in one sentence why the factor is 2/T and not 1/T.
2. What are b₁ and b₃ of the unit odd square wave, and why are the even b_n zero?
3. A signal has a kink but no jump. How fast do its coefficients fall, and how does the rms error of its N-term partial sum scale with N?
4. State Parseval's identity and use it with the square wave to evaluate 1 + 1/9 + 1/25 + ….
5. What is the Gibbs overshoot as a fraction of the jump, and what happens to it as N → ∞?
6. Write the complex coefficient c_n in terms of a_n and b_n, and explain why c_{−n} = c̄_n.
7. A mode at 7 Hz with ζ = 0.02 is driven by a 1 Hz square wave. Which harmonic resonates, what is its response amplitude in static deflections, and what is the phase of the response relative to that harmonic's forcing?
8. Give two pulse periods between 0.5 and 1.5 s that a 5 Hz boom's schedule must avoid, and one it may use.
