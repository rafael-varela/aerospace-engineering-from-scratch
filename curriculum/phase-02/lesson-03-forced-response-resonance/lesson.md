# P02.03  Forced Response & Resonance

> Drive a linear system with a sinusoid and it answers with the same sinusoid, scaled and delayed. The scale and the delay are one complex number, and near resonance that number is enormous.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the steady-state response of m ẍ + c ẋ + k x = F₀ cos ωt by undetermined coefficients and by the exponential response formula, and show the two agree through the complex gain H(iω). 2. Define magnification, phase lag, the resonance peak, half-power bandwidth and Q, and derive the exceptional resonant case (t sin ω_n t) and beats. 3. Implement and test a resonance-explorer module: complex gain, steady and full responses, peak and bandwidth, beats, step response, transmissibility, and an RK4 integrator checked against `scipy.integrate.solve_ivp`. 4. Decide, for an avionics bracket driven near its natural frequency, whether to stiffen it or to damp it, with numbers, and state where the linear forced model fails. |
| **Prerequisites** | `P02.02` (free response, ζ, ω_n, ω_d, the homogeneous solution), `P02.01` (RK4), `P01.05` (Euler's formula, complex exponentials), `P01.16` (the damped oscillator as a 2×2 system), `P00.03`, `P00.07`. Skill check: you can divide complex numbers, convert a + ib to R e^{iφ}, and differentiate t sin ωt. |
| **Engineering question** | A 2.5 kg avionics box sits on a bracket with a 55 Hz natural frequency and 2 % damping; the engine shakes it with a 20 N force whose frequency sweeps 50–60 Hz. How large is the steady vibration across the band, how long does it take to build up, and is stiffening the bracket to 80 Hz or adding damping to ζ = 0.1 the better fix? |
| **Assumptions and limits** | Linear spring and viscous damper with constant coefficients; a single degree of freedom; a sinusoidal force of fixed amplitude (an unbalance force actually grows as ω²: part 6); steady state reached (the transient takes 1/(ζω_n) × ln 100 ≈ 0.67 s here, so a sweep faster than a few Hz/s never reaches it: part 4). Real brackets have several modes, amplitude-dependent damping, and joints that stiffen or slip at large amplitude; the model is right for small motions and one mode at a time. |
| **Analytical work** | Undetermined coefficients for cos forcing; the exponential response formula and complex gain; magnification M(r, ζ) and phase lag; peak at r = √(1 − 2ζ²) by differentiation; half-power points and Q ≈ 1/(2ζ); the exceptional case by the limit of the beats formula and by direct substitution; the step response and its overshoot e^{−πζ/√(1 − ζ²)}; transmissibility and the universal crossing at r = √2. |
| **Computational deliverable** | `starter/main.py`: `natural_frequency`, `damping_ratio`, `complex_gain`, `magnification`, `phase_lag`, `undetermined_coefficients`, `steady_state_response`, `peak_frequency`, `peak_magnification`, `half_power_bandwidth`, `quality_factor`, `free_response`, `full_response`, `undamped_resonance`, `beats`, `step_response`, `step_overshoot`, `transient_decay_time`, `transmissibility`, `rk4_forced`. Tests prove: every closed form satisfies its differential equation by second differences; undetermined coefficients equal the real and imaginary parts of F₀H; the gain has the static, resonant (M = 1/(2ζ), lag 90°) and mass-controlled (1/r²) limits; the analytic peak beats 20 000 sampled frequencies and the half-power width is 2ζω_n to 1 %; the full response meets its initial conditions, decays to the steady state and agrees with RK4 (order 4) and `solve_ivp` to 10⁻¹⁰ m; the resonant solution grows linearly and the beats stay inside their envelope and tend to it; the step overshoot matches e^{−πζ/√(1 − ζ²)}; transmissibility is 1 at r = √2 for every ζ. |
| **Practical activity** | **computational experiment**: a frequency sweep of magnification and phase for five damping ratios; the build-up of the resonant amplitude against time; a chirp through resonance at four sweep rates to measure how much of the steady peak a fast sweep sees. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `resonance_explorer.png`, the hand derivation of H(iω) by both methods and of the exceptional case, the sweep-rate table with your predictions, and a short note recommending stiffening or damping for the bracket with the numbers behind it. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: both derivations of the steady state are complete and shown to agree; all public tests pass; hand values of M and φ at 50, 55 and 60 Hz agree with the code to three figures; the note gives a recommendation with the amplitude it achieves and the assumption it trusts least. |
| **References** | `@openstax-calculus-volume-3` §7.2 "Nonhomogeneous Linear Equations" (required-free). `@mit-ocw-18-03-differential-equations` Lecture 12 "Continuation: General Theory for Inhomogeneous ODEs", Lecture 13 "Finding Particular Solutions to Inhomogeneous ODEs", Lecture 14 "Interpretation of the Exceptional Case: Resonance" (required-free); notes 13 "Exponential response formula, spring drive", 14 "Complex gain, dashpot drive", 15 "Operators, undetermined coefficients, resonance", 16 "Frequency response" (supplementary-free). |
| **Extensions** | ★ Replace F₀ by the unbalance force m_e e ω² and show the high-frequency amplitude tends to m_e e/m. ★ Add a tuned absorber (a second mass-spring on the box) and show the original resonance splits into two. ★ Drive the system with a square wave built from its first five Fourier terms (preview of `P02.09`) and find which harmonic resonates. |
| **Time** | unmeasured |

## 1. The engineering problem

Every piece of equipment on an aircraft or launch vehicle sits on a structure that vibrates, and every structure has natural frequencies. When the frequency of a persistent disturbance, an engine's rotation, a propeller blade passage, a pump, lands on a natural frequency, the response can be tens of times larger than the static deflection the same force would cause. That is **resonance**, and it breaks brackets, fatigues solder joints, blurs sensors and, in the extreme case of aeroelastic flutter, destroys wings. The engineering job is to know where the resonances are, how sharp they are, and what the steady vibration will be across the frequency range the vehicle actually sees.

The model is the forced spring-mass-damper, **m ẍ + c ẋ + k x = F(t)**, the equation of `P02.02` with a right-hand side. `P02.02` solved the homogeneous equation and found that its solutions die away; this lesson finds the **particular solution** that persists, shows that for a sinusoidal force it is a sinusoid of the same frequency with a frequency-dependent amplitude and phase, and extracts the whole of vibration engineering's vocabulary from that one function: magnification, phase lag, bandwidth, Q, transmissibility.

The concrete problem: a 2.5 kg avionics box on a bracket with f_n = 55 Hz (k = 298.6 kN/m) and ζ = 0.02 is excited by a 20 N engine-induced force between 50 and 60 Hz. The static deflection F₀/k is 67 μm. Find the steady amplitude across the band, how fast it builds up, and compare two fixes: stiffening the bracket to 80 Hz, or adding damping to ζ = 0.1. Success criteria, written first: the steady-state closed form satisfies the equation to 10⁻⁴ relative; undetermined coefficients and the complex gain agree to 10⁻¹⁵; the peak magnification 1/(2ζ√(1 − ζ²)) exceeds every sampled frequency's; the full transient agrees with RK4 and `solve_ivp` to 10⁻⁹ m; and the undamped resonant solution grows linearly inside the predicted envelope.

## 2. Understand and derive

### 2.1 Structure of the solution: homogeneous plus particular

For a linear equation L[x] = F, if x_p is any one solution (a **particular solution**) and x_h is the general solution of L[x] = 0, then x = x_h + x_p is the general solution of L[x] = F (`@openstax-calculus-volume-3` §7.2; 18.03 Lecture 12): L[x_h + x_p] = 0 + F, and any two solutions differ by a solution of the homogeneous equation. For c > 0, `P02.02` showed every x_h decays as e^{−ζω_n t}; so after a time of order 1/(ζω_n) only x_p is left. x_h is the **transient**, x_p the **steady state**, and the initial conditions affect only the transient. For the bracket, ζω_n = 6.91 s⁻¹: the transient falls to 1 % in 0.666 s, 37 cycles. Everything in 2.2–2.5 is about x_p; 2.6 puts the transient back.

### 2.2 Undetermined coefficients

For F = F₀ cos ωt, guess x_p = A cos ωt + B sin ωt (Lecture 13). Substitute: ẋ_p = −Aω sin ωt + Bω cos ωt, ẍ_p = −ω²x_p, so

(k − mω²)(A cos ωt + B sin ωt) + cω(B cos ωt − A sin ωt) = F₀ cos ωt.

Match the cos and sin coefficients:

**(k − mω²) A + cω B = F₀,  −cω A + (k − mω²) B = 0.**   (1)

A 2×2 linear system (`P01.13`); the solution is A = F₀(k − mω²)/D, B = F₀cω/D with D = (k − mω²)² + (cω)². Then x_p = X cos(ωt − φ) with **X = F₀/√D** and **tan φ = cω/(k − mω²)**. The method works because the derivative of a sinusoid at ω is a sinusoid at ω: the guess is closed under the operator. It fails in exactly one case, when the guess is itself a homogeneous solution (c = 0 and ω = ω_n), and 2.5 is about that.

### 2.3 The exponential response formula and the complex gain

The same answer with less algebra (notes 13, 14): since cos ωt = Re e^{iωt}, solve the complex equation m z̈ + c ż + k z = F₀ e^{iωt} and take the real part at the end (the equation has real coefficients, so real and imaginary parts separate). Guess z = Z e^{iωt}: (−mω² + icω + k) Z = F₀, so

**z_p = F₀ H(ω) e^{iωt},   H(ω) = 1/(k − mω² + icω)**.   (2)

H is the **complex gain** (units m/N): one complex number per frequency that holds both the amplitude ratio |H| = 1/√D and the phase lag φ = −arg H. Then x_p = Re[F₀ H e^{iωt}] = F₀|H| cos(ωt − φ), which is 2.2 again: A = F₀ Re H, B = −F₀ Im H (the test checks these identities to 10⁻¹⁵). The exponential response formula is the method you will use for the rest of your career, because it turns a differential equation into division: every linear time-invariant system has an H(ω), `P02.08` will call it a transfer function evaluated at s = iω, and `P02.10` will measure it from data.

### 2.4 Magnification and phase: the frequency response

Write r = ω/ω_n (the **frequency ratio**) and divide numerator and denominator of (2) by k. The **magnification** M = X/(F₀/k) = k|H| is the steady amplitude in units of the static deflection:

**M(r, ζ) = 1/√((1 − r²)² + (2ζr)²),   tan φ = 2ζr/(1 − r²)**.   (3)

Three regimes, visible in `resonance_explorer.png`:

- **Stiffness-controlled, r ≪ 1:** M → 1, φ → 0. The mass follows the force quasi-statically; x ≈ F/k.
- **Resonance, r ≈ 1:** at r = 1 exactly, M = 1/(2ζ) and φ = 90° whatever ζ. The spring and inertia forces cancel and only the damper resists: X = F₀/(cω_n). For the bracket M = 25.0: 20 N produces 1.675 mm, 25 times the static 67 μm.
- **Mass-controlled, r ≫ 1:** M → 1/r², φ → 180°. The force reverses too fast for the spring to matter; x ≈ −F/(mω²), opposite in sign to the force.

The peak of M is not exactly at r = 1. Differentiate the denominator (1 − r²)² + 4ζ²r² with respect to r² and set it to zero: −2(1 − r²) + 4ζ² = 0, so

**r_peak = √(1 − 2ζ²),   M_peak = 1/(2ζ√(1 − ζ²))**,   (4)

and there is no peak at all for ζ ≥ 1/√2 (M then decreases monotonically from 1). For the bracket r_peak = 0.9996, 54.98 Hz, M_peak = 25.01: indistinguishable from the r = 1 values at this damping, which is why the light-damping shorthand M_peak ≈ 1/(2ζ) is universal.

**Half-power bandwidth and Q.** The two frequencies at which M² has fallen to half its peak (−3 dB) bracket the resonance. For light damping, expand (3) near r = 1 with r = 1 + ε: (1 − r²)² ≈ 4ε², so M² ≈ 1/(4ε² + 4ζ²), which is half its maximum at ε = ±ζ. The band is **Δω ≈ 2ζω_n**, and the **quality factor** Q = ω_peak/Δω ≈ 1/(2ζ). For the bracket: 53.87–56.07 Hz, Δf = 2.20 Hz, Q = 25.0. The code finds the half-power points with `brentq` on (3) rather than the approximation, and the test confirms the approximation to 1 %. Q is what a vibration test measures directly: sweep, read the peak and its width, and you have ζ without ever seeing a decay record.

The bracket across its drive band (`main()` output):

```
 f [Hz]       M   X [um]  phase lag [deg]  T (base)
   40.0    2.12    141.9              3.5      2.12
   50.0    5.64    377.8             11.8      5.64
   53.0   12.32    825.6             28.4     12.33
   55.0   25.00   1674.7             90.0     25.02
   57.0   11.78    789.4            150.8     11.79
   60.0    5.13    343.5            167.1      5.13
   80.0    0.90     60.0            177.0      0.90
  110.0    0.33     22.3            178.5      0.33
```

The phase column is as diagnostic as the amplitude column: a response that lags the force by 90° *is* at resonance, and a phase that swings from near 0 to near 180° across a few hertz is the signature of light damping. The two fixes: stiffening to 80 Hz puts the 55 Hz drive at r = 0.69 and M = 1.90 (127 μm, a 13-fold reduction); raising ζ to 0.10 leaves the resonance in the band and caps M at 5.0 (335 μm, a 5-fold reduction). Stiffening wins if the drive band is fixed; damping wins if the drive frequency can be anywhere, because a stiffened bracket still has a resonance somewhere and Q = 25 at it.

### 2.5 The exceptional case: resonance without damping, and beats

With c = 0 and ω = ω_n, the guess A cos ω_n t + B sin ω_n t is a homogeneous solution, so it cannot produce the right-hand side: D = 0 in (1) and H in (2) is infinite. The fix (Lecture 14, notes 15) is the usual one for a repeated root: multiply the guess by t. Try x_p = t(A cos ω_n t + B sin ω_n t); substituting into m ẍ + k x = F₀ cos ω_n t, the t-terms cancel because of ω_n² = k/m and what remains is 2mω_n(B cos ω_n t − A sin ω_n t) = F₀ cos ω_n t, so A = 0, B = F₀/(2mω_n):

**x_p = (F₀/(2mω_n)) t sin ω_n t**.   (5)

The amplitude grows linearly for ever: 11.6 mm per second for the bracket if it had no damping. Real structures have damping, so the growth saturates at F₀/(cω_n) after a time of order 1/(ζω_n); but a lightly damped structure driven at resonance spends that time at ever-increasing amplitude, and for ζ = 0.02 the first 0.1 s already reaches 0.81 mm, half the final value.

**Beats** are the undamped response just *off* resonance, from rest: with x_p = F₀ cos ωt/(m(ω_n² − ω²)) and the homogeneous part chosen to make x(0) = ẋ(0) = 0,

**x = F₀/(m(ω_n² − ω²)) (cos ωt − cos ω_n t) = (2F₀/(m(ω_n² − ω²))) sin((ω_n − ω)t/2) sin((ω_n + ω)t/2)**.   (6)

A fast oscillation at the mean frequency inside a slow envelope at |ω_n − ω|/2: the energy flows into the oscillator and back out with period 2π/|ω_n − ω|. For a 52 Hz drive on the 55 Hz bracket the envelope repeats every 1/3 s. Let ω → ω_n in (6): sin((ω_n − ω)t/2) ≈ (ω_n − ω)t/2 and ω_n² − ω² ≈ 2ω_n(ω_n − ω), so x → (F₀/(2mω_n)) t sin ω_n t, which is (5). Resonance is the limit of beats whose period has become infinite; the test checks that limit numerically at ω = 1.0001 ω_n.

### 2.6 The full response and the step response

With damping, the complete solution from x(0) = x₀, ẋ(0) = v₀ is x = x_h + x_p where x_h is the `P02.02` free response started from (x₀ − x_p(0), v₀ − ẋ_p(0)), so that the *sum* meets the initial conditions. At resonance from rest the amplitude of the bracket grows to 0.49, 0.81, 1.24, 1.62 mm after 0.05, 0.1, 0.2 and 0.5 s, approaching 1.675 mm as 1 − e^{−ζω_n t}: the build-up time constant is the same 1/(ζω_n) = 0.145 s that governs the decay. This matters for a sweep: if the drive passes through the resonance faster than the amplitude can build, the structure never sees the full peak (part 4).

A constant force F₀ switched on at t = 0 (a **step**) has the particular solution x_p = F₀/k and, from rest,

**x = (F₀/k)[1 − e^{−ζω_n t}(cos ω_d t + (ζω_n/ω_d) sin ω_d t)]**,   (7)

which overshoots its final value by the fraction e^{−πζ/√(1 − ζ²)} at the first peak (t = π/ω_d): 93.9 % for the bracket, 4.6 % for ζ = 0.7, 0 for ζ ≥ 1. The step response is the bridge to control theory (`P02.08`), where overshoot and settling time are the specifications.

### 2.7 Base excitation and transmissibility

Often the force is not applied to the mass; the *base* moves, y = Y cos ωt, and the question is how much motion reaches the equipment. With relative coordinate z = x − y the equation is m z̈ + c ż + k z = −m ÿ, and the absolute motion of the mass is x = y + z. Working it through in complex form, x/y = (k + icω)H(ω), so the **transmissibility**

**T = |x|/|y| = √((1 + (2ζr)²)/((1 − r²)² + (2ζr)²))**.   (8)

Two facts to carry away. Every T curve, whatever ζ, passes through T = 1 at r = √2: below √2 the mount amplifies, above it isolates. And above r = √2, *more* damping transmits *more* (the damper is a second path for force at high frequency). A vibration isolator is therefore designed soft, so that the lowest disturbance frequency is well above √2 ω_n, and only as damped as needed to survive passing through resonance at start-up. For the bracket, T and M are indistinguishable below 60 Hz (the (2ζr)² term is 0.0016) and separate only at high r.

### 2.8 Where the model stops being valid

- **Constant force amplitude.** A rotating unbalance produces F₀ = m_e e ω², so the high-frequency response tends to the constant m_e e/m rather than falling as 1/r²; extension 1.
- **One mode.** A real bracket has a second bending mode and a torsion mode within a factor of three of the first; each has its own H, and they add. `P02.10` separates them in the frequency domain.
- **Linear damping.** Joint friction and material damping depend on amplitude; Q measured at small amplitude is usually higher than at large.
- **Steady state.** The sweep-rate experiment in part 4 shows a 20 Hz/s sweep sees only 78 % of the steady peak.
- **Small motion.** A 1.7 mm steady amplitude at 55 Hz is a peak acceleration of 200 m/s² (20 g) on the box; bolted joints slip and the spring is no longer linear long before that.

## 3. Build the computational model

Open `starter/main.py`. `avionics_bracket()` is given; everything else is yours.

### Step 1: the complex gain and what it unpacks into

```python
def complex_gain(omega, m, c, k):
    omega = np.asarray(omega, dtype=float)
    denom = k - m * omega ** 2 + 1j * c * omega
    if np.any(denom == 0.0):
        raise ValueError("undamped system driven exactly at resonance")
    return 1.0 / denom                                                 # (2)

def magnification(omega, m, c, k):  return k * np.abs(complex_gain(omega, m, c, k))
def phase_lag(omega, m, c, k):      return np.arctan2(c * omega, k - m * omega ** 2)

def steady_state_response(t, omega, F0, m, c, k):
    return np.real(F0 * complex_gain(omega, m, c, k) * np.exp(1j * omega * t))
```

Use `arctan2` for the phase so that it runs 0 → π through 90° rather than wrapping. `undetermined_coefficients` solves the 2×2 system (1) with `np.linalg.solve`; the test compares it with Re and −Im of F₀H to 10⁻¹⁵, which is the whole point of writing both.

### Step 2: the peak

`peak_frequency` and `peak_magnification` are (4) with the ζ ≥ 1/√2 branch returning (0, 1). `half_power_bandwidth` finds the two roots of M² − ½M_peak² with `brentq`, one on each side of the peak; `quality_factor` divides.

### Step 3: transients

`free_response` is the underdamped formula from `P02.02`; `full_response` adds the steady state with the initial conditions corrected by x_p(0) = Re(F₀H) and ẋ_p(0) = Re(iωF₀H). `undamped_resonance` is (5), `beats` is (6), `step_response` is (7), `step_overshoot` is e^{−πζ/√(1 − ζ²)}, `transient_decay_time` is −ln(fraction)/(ζω_n).

### Step 4: transmissibility and the RK4 oracle

`transmissibility` is |(k + icω)H|. `rk4_forced(m, c, k, force, x0, v0, t_end, n)` integrates d/dt(x, v) = (v, (force(t) − cv − kx)/m) with the stepper of `P02.01`; it takes any callable `force(t)`, which is what makes the chirp experiment in part 4 possible.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-03-forced-response-resonance`. The tests:

1. **Closed forms satisfy their equations.** Steady state at 20, 54, 55 and 70 Hz; the full response; the resonant solution; the beats; the step response, each by second differences to 10⁻³–10⁻⁴ of F₀.
2. **Two methods agree.** A = F₀ Re H and B = −F₀ Im H to 10⁻¹⁵; A cos + B sin equals the steady-state function.
3. **Limits.** M(0) = 1, φ(0) = 0; M → 1/r² and φ → π at 50 ω_n; M(ω_n) = 1/(2ζ) and φ(ω_n) = π/2 exactly; the undamped gain at ω_n raises; T(√2 ω_n) = 1 for four damping ratios; T(0) = 1.
4. **The peak.** Analytic peak ≥ the maximum over 20 000 samples and equal to it to 10⁻⁶; half-power points have M² = ½M_peak² to 10⁻⁸; Δω = 2ζω_n and Q = 1/(2ζ) to 1 %; ζ = 0.8 has no peak and `half_power_bandwidth` raises.
5. **Transients.** Initial conditions met; steady state reached after 1.5 decay times; RK4 and `solve_ivp` within 10⁻⁹ and 10⁻¹⁰ m; RK4 order 4.0 ± 0.4; two drives superpose; the resonant amplitude grows linearly (peak ratios 0.25/0.15, 0.35/0.25, 0.45/0.35 within 10 %); beats vanish at multiples of 2π/|ω_n − ω| and tend to the resonant solution at ω = 1.0001 ω_n; step overshoot equals the formula to 10⁻³.

### Compare

`solve_ivp` (RK45, rtol 10⁻¹⁰) agrees with the closed-form full response at resonance to 1.6 × 10⁻¹² m over one second; our RK4 with 20 000 steps agrees to 6 × 10⁻¹¹ m. Both are irrelevant as *solvers* for this problem, since the closed form is exact and instant, and both are essential as *checks*: the closed form in 2.6 has four places to drop a sign (the two initial-condition corrections and the two Im parts), and the integrator catches every one. In practice the frequency response of a linear system is computed exactly as `complex_gain` does, by evaluating H at each frequency; `scipy.signal.freqs` does the same for any transfer function, and `P02.08` connects the two. The integrator is reserved for what H cannot do: a nonlinear spring, amplitude-dependent damping, or a force that is not a sum of sinusoids (though `P02.09` will show that almost everything is).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Frequency sweep.** M and φ from 20 to 120 Hz for ζ = 0.02, 0.05, 0.1, 0.2, 0.5 (`resonance_explorer.png`, top row). Predict: the peak heights (1/(2ζ)), the frequency at which all phase curves cross (ω_n, 90°), and whether the curves cross anywhere in amplitude (they do not; more damping is lower everywhere below r = √2 for M; check T separately).
2. **Build-up at resonance.** From rest at 55 Hz, record the envelope every 0.05 s and fit 1 − e^{−t/τ}. Predict τ = 1/(ζω_n) = 0.145 s.
3. **Chirp through resonance.** Drive the bracket with F₀ cos(2π(45t + 10t²/T)) (a linear sweep from 45 to 65 Hz over T seconds) using `rk4_forced` with 20 000 steps per second, for T = 1, 3, 15, 60 s. Record the peak |x| and the instantaneous frequency at which it occurs. Predict first: the fraction of the steady peak (1.675 mm) seen at each rate, and whether the apparent resonance shifts up or down in frequency.
4. **Amplitude limit band.** If the box tolerates 0.5 mm, find the frequency band in which the steady response exceeds it for ζ = 0.02 and ζ = 0.05, and compare each band with the half-power bandwidth.

Raw data format: CSV tables (`f_hz, M, phase_deg` per ζ; `t, envelope_m`; `T_sweep, rate_hz_per_s, peak_m, f_at_peak_hz`; `zeta, f_lo, f_hi`), each with a header naming the fixed parameters.

What you should find: sweeps of 20, 6.7, 1.3 and 0.33 Hz/s reach 1.306, 1.523, 1.654 and 1.673 mm (78 %, 91 %, 99 %, 99.9 % of the steady peak), with the apparent peak at 57.6, 56.2, 55.3 and 55.1 Hz: a fast upward sweep *under-reads* the peak and *over-reads* its frequency, because the amplitude is still building when the drive has already passed ω_n. The rule of thumb that falls out is that the sweep must spend several build-up times 1/(ζω_n) inside the half-power band 2ζf_n, i.e. a rate well below 2ζ²f_n²/… in practice below about ζ²f_n² = 1.2 Hz/s here. The 0.5 mm band is 51.3–58.4 Hz (7.1 Hz) at ζ = 0.02, three times the half-power width, and 52.3–57.3 Hz (4.9 Hz) at ζ = 0.05, where the peak itself is only 0.67 mm.

## 5. Validate and improve

| quantity (bracket) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| M at 55 Hz | 1/(2ζ) = 25.0 | 25.00 | 0 | exact |
| X at 55 Hz (mm) | 67.0 μm × 25 = 1.675 | 1.6747 | 0 | exact |
| M at 50 Hz | 1/√((1 − 0.826)² + (0.0364)²) = 5.64 | 5.64 | 0 | exact |
| φ at 55 Hz | 90° | 90.0° | 0 | exact |
| peak frequency (Hz) | 55√(1 − 0.0008) = 54.98 | 54.978 | 0 | exact |
| half-power width (Hz) | 2ζf_n = 2.20 | 2.20 (brentq) | 0.04 % | light-damping approximation |
| Q | 25 | 25.0 | 0 | same |
| 1 % decay time (s) | ln 100/(ζω_n) = 0.666 | 0.666 | 0 | exact |
| step overshoot | e^{−πζ/√(1 − ζ²)} = 93.9 % | 93.9 % (max of (7)) | < 0.1 % | grid |
| RK4 vs closed form (m) | — | 6 × 10⁻¹¹ | — | h = 50 μs, order 4 |
| `solve_ivp` vs closed form (m) | — | 1.6 × 10⁻¹² | — | rtol 10⁻¹⁰ |

The numbers are self-consistent to the last digit, and that is the least important thing about them. What they rest on is three inputs nobody knows to better than 10–20 %: the damping ratio (bolted aluminium brackets range from 0.5 % to 5 %), the natural frequency (a 5 % error in k is 2.5 % in f_n, and at Q = 25 a 2.5 % frequency error changes M at a fixed drive frequency by a factor 2), and the force amplitude (engine-induced loads are measured, not calculated). The sharper the resonance, the less a prediction at a single frequency means; the honest deliverable is the whole curve M(f) with a band of uncertainty, and the design rule that follows from it: keep every drive frequency out of the band f_n(1 ± 3ζ) or add damping until the band is wide and low.

What would validate the model is a shaker test: mount the box on an electrodynamic shaker, sweep slowly (part 4 says how slowly: below 1 Hz/s for this bracket), measure the transfer function with an accelerometer, and read f_n, Q and the amplitude at the drive frequencies directly. The model's two predictable failures are a second peak within a factor of three of the first (a mode the single-degree-of-freedom model does not have) and a Q that falls as the drive level rises (amplitude-dependent damping). Both show up in the shaker data, and neither can be fixed by tuning m, c and k.

Improvements worth making: replace F₀ by m_e e ω² and redo the sweep (extension 1); add the tuned absorber (extension 2) and see the single resonance become two, with a notch between them at the original f_n, which is how helicopter cabins are quietened. Improvement worth *not* making: a finer time step. The integrator is already seven orders of magnitude below the uncertainty in ζ.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of H(ω) by undetermined coefficients and by the exponential response formula, with the identities A = F₀ Re H, B = −F₀ Im H; the derivation of (5) by substitution; `resonance_explorer.png`; the sweep-rate table with predictions; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending stiffening or damping for the bracket, giving the amplitude each achieves and the assumption you trust least.

### Adversarial self-checks

- **What if the force doubles?** Every amplitude doubles and nothing else changes: the response is linear in F₀. What if the *mass* doubles? f_n falls to 38.9 Hz and ζ falls to 0.014 (same c); the 50–60 Hz drive band is now at r = 1.29–1.54, M = 1.53–0.72. Check with the code; a heavier box on the same bracket is *better* here, which is the mass-controlled regime at work.
- **What if ω → 0? ω → ∞?** M → 1 (static) and M → 1/r² → 0 (the mass cannot follow); φ → 0 and φ → 180°. The 1/r² limit is why isolation works.
- **Which assumption fails first in the real world?** Single degree of freedom: the bracket's second mode. Second: amplitude-independent damping.
- **What measurement could falsify the model?** A shaker sweep at two force levels: if Q differs between them, the damping is not viscous; if there is a second peak below 165 Hz, the model is missing a mode.
- **Does the result survive a change of solver or sweep rate?** Closed form, RK4 and `solve_ivp` agree to 10⁻¹⁰ m. The *measured* peak does not survive a change of sweep rate: 78 % at 20 Hz/s against 99.9 % at 0.33 Hz/s. Which decision changes if a test engineer swept at 10 Hz/s and reported Q = 20?
- **Which decision changes if the uncertainty doubles?** If ζ is known only as 0.01–0.04 rather than 0.015–0.025, the peak amplitude is anywhere between 0.84 and 3.35 mm; stiffening to 80 Hz gives 127 μm *almost independently of ζ* (M = 1.90 at r = 0.69 is stiffness-controlled), while adding damping to a nominal ζ = 0.1 could be adding to 0.05 or 0.15. Under doubled uncertainty, stiffening is the robust choice.

## Check yourself

1. State why the general solution of m ẍ + c ẋ + k x = F(t) is x_h + x_p, and which part the initial conditions determine.
2. Write H(ω) for the forced oscillator and give |H| and the phase lag at ω = ω_n in terms of c and k.
3. What is the magnification M(r, ζ)? Evaluate it at r = 0.5, 1 and 2 for ζ = 0.02.
4. Where exactly is the peak of M, and for which ζ is there no peak?
5. The half-power bandwidth of a resonance is 2.2 Hz at 55 Hz. What are Q and ζ?
6. Why does the guess A cos ω_n t + B sin ω_n t fail for the undamped system driven at ω_n, and what replaces it?
7. A 52 Hz drive on an undamped 55 Hz oscillator from rest gives beats. What is the beat period and what is the envelope frequency?
8. At what frequency ratio does every transmissibility curve pass through 1, and what does more damping do above it?
