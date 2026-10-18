# P09.17  Digital Control

> A flight computer does not see the aircraft; it sees samples, thinks for a while, and holds its answer until the next sample. Each of those three acts costs phase, and you can say in advance exactly how much.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the zero-order-hold (ZOH) equivalent of a state-space plant from the augmented matrix exponential, extend it to a computational delay 0 ≤ τ ≤ T, and show it is exact at the sample instants (step invariance). 2. Derive the Tustin (bilinear) and matched pole-zero discretisations, the frequency-warping law ω_d = (2/T)arctan(ω_aT/2) and prewarping, and implement all three from scratch, checked against `scipy.signal.cont2discrete` and `scipy.linalg.expm`. 3. Use the z-transform tools (difference equations, inverse z-transform by long division, frequency response on the unit circle, the z-plane stability test |z| < 1) on a digital pitch-attitude loop. 4. Predict the phase margin a sample rate and a computational delay cost, ΔPM ≈ ω_c(T/2 + τ), verify it against the margin of the exact sampled loop, and choose a sample rate. |
| **Prerequisites** | `P09.13` (state-space models), `P09.10` (Nyquist, phase margin, delay margin), `P09.07` (digital filters, the difference equation), `P09.06` (sampling and aliasing), `P09.11` (pitch-attitude loop), `P02.05` (matrix exponential, the augmented-matrix identity), `P02.08` (Laplace transform). Skill check: you can compute e^{AT} for a 2 × 2 matrix and evaluate a polynomial at a complex number. |
| **Engineering question** | A small UAV's pitch-attitude loop, designed in continuous time with a 66° phase margin at a 7.1 rad/s crossover, is to run on a flight computer that samples the attitude, computes for up to one sample period and then holds the elevator command. At 100, 50, 20, 10 and 7 Hz, with no delay, half a sample and a full sample of delay, how much margin and damping survive, which rates are unsafe, and how well does a one-line rule predict the result? |
| **Assumptions and limits** | Linear time-invariant plant (a representative short-period pitch model of the same form as `P09.11`, with a stiffer M_α, and a first-order elevator servo); ideal, synchronous sampling at a fixed period T; the controller output is held constant between updates (ZOH) and appears after a fixed delay τ ≤ T; no quantisation, no sensor noise, no anti-aliasing filter (its phase lag is added in Problem 5); a delay longer than one sample needs more augmented states. The equivalent-delay rule is a first-order approximation, accurate when ω_cT ≪ 1. |
| **Analytical work** | The ZOH equivalent and the augmented matrix; step invariance; the delay formulas G₀, G₁ and the augmented state (x_k, u_{k−1}); the z-transform of the difference equation and the z-plane stability test; Tustin's map, its warping law and prewarping; matched pole-zero with zeros at infinity mapped to z = −1 and DC-gain matching; the ZOH frequency response T e^{−jωT/2} sinc(ωT/2); the phase-margin rule PM_d ≈ PM_c − ω_c(T/2 + τ). |
| **Computational deliverable** | `starter/main.py`: `c2d_zoh`, `c2d_tustin`, `c2d_zoh_delay`, `tustin_tf`, `matched_tf`, `dlsim`, `impulse_by_long_division`, `freqresp_z`, `spectral_radius`, `z_to_s`, `discrete_pitch_loop` (`expm`, `pitch_model`, `pitch_gains`, `freqresp_s`, `step_metrics`, `hybrid_simulation`, `discrete_loop_margin`, `sample_rate_study` and `main` are given). Tests prove: ZOH matches `cont2discrete(method="zoh")` and `expm` to 10⁻¹⁰ on random systems and the closed form A⁻¹(e^{AT} − I)B; the double-integrator ZOH is exact; the discrete step response equals `scipy.signal.step` at the samples to 10⁻⁹ while Tustin's does not; Tustin in state space and as a transfer function matches `cont2discrete(method="bilinear")`; the warping law holds and prewarping is exact at its frequency; matched pole-zero maps poles to e^{pT}, keeps the DC gain, and puts zeros at infinity at z = −1; long division equals `scipy.signal.dimpulse`; `dlsim` is linear and time invariant; ZOH and Tustin keep stable systems stable; the delay discretisation reduces to ZOH at τ = 0 and to a pure one-sample delay at τ = T, has G₀ + G₁ = B_d, and matches a simulation resolved between samples to 10⁻¹² for four delays; fast sampling recovers the continuous poles; slower sampling and longer delay cost damping and margin as the rule predicts; 7 Hz with a full-sample delay is unstable. |
| **Practical activity** | **computational experiment**: a sample-rate and delay sweep of the pitch loop (spectral radius, equivalent damping, measured phase margin against the rule); a comparison of Tustin, prewarped Tustin and matched pole-zero for a lead compensator and (Problem 3) a bending-mode notch filter; step-invariance and between-sample checks. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `digital_pitch.png`, the derivations of the augmented-matrix ZOH, the delay formulas and the Tustin warping law, the sample-rate table with the rule column, and a half-page note recommending a sample rate and a maximum computational delay for the pitch loop. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the ZOH and delay formulas are derived; all public tests pass; the margin table is reproduced and the rule's error is explained; the note states a rate and a delay budget with the margin they leave. |
| **References** | `@mit-ocw-16-06-principles-of-automatic-control` Lecture 29 "Digital control, the z-transform" (https://ocw.mit.edu/courses/16-06-principles-of-automatic-control-fall-2012/resources/mit16_06f12_lecture_29/): Discrete Time Control Time line; Z-Transform (required-free): the sample-compute-hold cycle and z-plane stability. Lecture 30 "The z‐transform, design by emulation, the Tustin transform" (PDF heading "Z-transform Inversion"): Design by Emulation; Tustin's approximation; matched pole-zero (MPZ) (required-free): the controller discretisations of §2.5. Lecture 31 "Compensator design examples, time delay of ZOH"; Lecture 32 "Discrete design"; Lecture 34 "Design examples, pre‐warping, direct design" (required-free): the ZOH's half-sample delay and its phase lag at crossover, and prewarping. `@mit-ocw-16-30-feedback-control-systems` Topic 20 "Digital control basics" (https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec20/): Digital Control Mechanics; Sample and Hold Analysis ("Effective Delay"); "So why not just make Ts small?"; Control Law; Computer Code Layout (required-free): the effective delay of sampling and computation as lost phase margin. `@mit-ocw-6-003-signals-and-systems` Lecture 7 "Discrete Approximation of Continuous-Time Systems": Forward Euler, Backward Euler and Trapezoidal Rule: Mapping CT poles to DT poles; Lecture 5 "Z Transform" (supplementary-free): where each approximation sends the s-plane. `@nasa-tn-d-7843-f8-digital-fly-by-wire` (supplementary-free; US Government work, historical): Deets, paper 2, FLIGHT CONTROL SOFTWARE → Control Law Modes and DIGITAL FLY-BY-WIRE DESIGN → Design Synthesis and Analysis (https://ntrs.nasa.gov/citations/19750010175); Szalai, paper 6, CONTROL SYSTEM PERFORMANCE → Direct Mode and Stability Augmentation System Mode (https://ntrs.nasa.gov/citations/19750010179): how a real flight-control computer chose its sample rate and what quantisation and sampling did in flight. `@scipy-userguide` API reference "Signal processing (scipy.signal)": "Discrete-time linear systems" and `cont2discrete`; "Linear algebra (scipy.linalg)": `expm` (supplementary-free): the test oracles only. Partial GAP (original instruction): the augmented-matrix (Van Loan) identity (2) for the ZOH equivalent is not derived in the loaded notes; it is derived here, building on `P02.05`. |
| **Extensions** | ★ Add a delay of 1.5 samples by augmenting two past inputs and check against `hybrid_simulation` extended to it. ★ Design the pitch controller directly in the z-domain by placing the closed-loop poles at e^{s_iT} with Ackermann's formula on (Φ, Γ) (`P09.15`) and compare its margin at 10 Hz with the emulated design. ★ Model timing jitter by drawing τ uniformly in [0, T] at each sample (seeded) and estimate the spread of the overshoot. |
| **Time** | unmeasured |

## 1. The engineering problem

Every modern autopilot is a computer program. Once per sample period T it reads the sensors, runs the control law, and writes a new command to the actuator, which then holds that command until the next update. The continuous-time designs of `P09.11`–`P09.16` assume the controller acts instantly and continuously. A flight computer does neither: it acts every T seconds, its output is a staircase, and the computation itself takes time, so the command for the state measured at t_k reaches the servo at t_k + τ.

The questions for the engineer are concrete. What is the slowest sample rate at which the continuous design still works? How much does a computational delay of half a sample or a whole sample cost? And how should a continuous filter, such as a lead compensator or a notch for a structural mode, be turned into a difference equation the computer can run? Get these wrong and an autopilot that was stable in every analysis oscillates in flight; the mechanism is lost phase margin, which is exactly computable.

The example is the pitch-attitude loop of a small UAV: a short-period model of the same form as `P09.11` (α̇ = −L_αα + q, q̇ = M_αα + M_qq + M_δδ_e, θ̇ = q, with L_α = 1.2 s⁻¹, M_q = −2.5 s⁻¹ and M_δ = 12 s⁻² as there but a stiffer M_α = −9 s⁻² instead of −6 s⁻²) with a 0.05 s elevator servo, and the attitude law δ_c = K_p(θ_ref − θ) − K_dq with K_p = 2.0 and K_d = 0.5 s. Success criteria, written first: our ZOH matches `scipy.signal.cont2discrete` and `scipy.linalg.expm` to 10⁻¹⁰; the discrete step response equals the continuous one at the samples to 10⁻⁹; the delay model equals a simulation resolved between samples to 10⁻¹²; and the rule PM_d ≈ PM_c − ω_c(T/2 + τ) predicts the exact sampled margin within 2° at 20 Hz and above.

## 2. Understand and derive

### 2.1 Sampling, holding and the ZOH equivalent

Let the plant be ẋ = Ax + Bu. The computer applies u(t) = u_k for kT ≤ t < (k + 1)T. Over one period the variation-of-parameters formula (`P02.05`) gives exactly

**x_{k+1} = A_d x_k + B_d u_k,  A_d = e^{AT},  B_d = ∫₀ᵀ e^{As} ds B.**   (1)

Both come from one exponential of the **augmented matrix**:

exp([[A, B], [0, 0]] T) = [[A_d, B_d], [0, I]],   (2)

because the series of the left side has upper-right block Σ_{k≥1} A^{k−1}BT^k/k! = (∫₀ᵀ e^{As}ds)B (Problem 2a). No quadrature, no approximation: (1) is exact for a held input. (16.30 Topic 20 and 16.06 Lecture 29 set up the same sample-and-hold picture; the augmented-matrix route to it is this lesson's.) That is **step invariance**: a step is a held input, so the discrete step response equals the continuous step response at every sample instant. The test checks this against `scipy.signal.step` to 10⁻⁹ (the measured difference is 5 × 10⁻¹⁵).

For small T, A_d ≈ I + AT and B_d ≈ BT, which is forward Euler; (1) is its exact replacement. The eigenvalues map as **z = e^{sT}**: each continuous pole s becomes a discrete pole e^{sT}. For the pitch plant at 50 Hz the open-loop poles 0, −1.85 ± 2.929j and −20 s⁻¹ become 1, 0.962 ± 0.0564j and 0.6703, equal to e^{λT} to 1.1 × 10⁻¹⁶.

### 2.2 The z-transform and stability

Take the z-transform of x_{k+1} = A_dx_k + B_du_k, y_k = C_dx_k + D_du_k with zero initial state: zX(z) = A_dX(z) + B_dU(z), so G_d(z) = C_d(zI − A_d)⁻¹B_d + D_d, a ratio of polynomials in z. Three consequences used all lesson:

- **Stability.** The solution of x_{k+1} = Mx_k is Σ c_i λ_iᵏ v_i, so it decays iff every eigenvalue has |λ_i| < 1: the discrete system is stable iff its **spectral radius** max|λ_i| is below 1. The map z = e^{sT} sends the left half plane into the unit disk (|e^{sT}| = e^{Re(s)T}), so ZOH preserves stability.
- **Frequency response.** A sinusoid sampled at frequency ω is z^k with z = e^{jωT}, so the frequency response of G_d is G_d(e^{jωT}), evaluated on the unit circle, periodic in ω with period 2π/T and meaningful up to the Nyquist frequency π/T.
- **Inverse transform by long division.** Dividing numerator by denominator in powers of z⁻¹ gives h₀ + h₁z⁻¹ + h₂z⁻² + …, and the coefficients are the impulse response. For 1/(z − 0.6) this gives 0, 1, 0.6, 0.36, 0.216, … (`impulse_by_long_division`, checked against `scipy.signal.dimpulse`).

### 2.3 The computational delay

Let the controller output computed from x_k reach the servo at t_k + τ, with 0 ≤ τ ≤ T. During [t_k, t_k + τ) the plant still sees u_{k−1}; during [t_k + τ, t_{k+1}) it sees u_k. Integrating each piece exactly:

**x_{k+1} = e^{AT}x_k + G₁u_{k−1} + G₀u_k,  G₀ = ∫₀^{T−τ} e^{As}ds B,  G₁ = e^{A(T−τ)}∫₀^τ e^{As}ds B.**   (3)

The old input u_{k−1} is now part of the state. With z_k = (x_k, u_{k−1}),

z_{k+1} = Φ z_k + Γ u_k,  Φ = [[e^{AT}, G₁], [0, 0]],  Γ = [[G₀], [I]].   (4)

Limits: τ = 0 gives G₁ = 0, G₀ = B_d, plain ZOH; τ = T gives G₀ = 0, G₁ = B_d, a pure one-sample delay, which adds a pole at z = 0. For any τ, G₀ + G₁ = B_d (a constant input does not care when it switched). The tests check all three, and check (4) against `hybrid_simulation`, which advances the plant on a grid of T/20 and switches the input at t_k + τ without ever forming G₀ or G₁: the two agree to 10⁻¹² at every sample for τ = 0, T/4, T/2 and T.

### 2.4 Where the phase goes: the rule

Why does sampling cost margin? A ZOH turns the impulse train of samples into a staircase. Its frequency response, relative to an ideal sampler, is

H_ZOH(jω) = (1 − e^{−jωT})/(jω) = T e^{−jωT/2} sinc(ωT/2),  sinc x = sin x / x,   (5)

(Problem 2d): unit gain at low frequency and a **phase lag of ωT/2**, the same as a delay of half a sample (the "Effective Delay" of 16.30 Topic 20; 16.06 Lecture 31). The staircase is, on average, half a sample late. A computational delay τ adds ωτ. At the crossover ω_c, the phase margin therefore falls by about

**ΔPM ≈ ω_c (T/2 + τ)** [rad].   (6)

For the pitch loop, the continuous loop broken at the servo input has ω_c = 7.12 rad/s and PM = 65.7°. At 20 Hz (T = 0.05 s) the hold costs 7.12 × 0.025 rad = 10.2° and a full-sample delay another 20.4°, predicting 35.1°. At 10 Hz the prediction is 4.5°: a loop on the edge of instability, from a design that had 66°.

The rule is checked against the exact margin of the sampled loop, L(e^{jωT}) = [K, 0](e^{jωT}I − Φ)⁻¹Γ, computed by `discrete_loop_margin` (part 4). The rule is a first-order expansion in ω_cT. The data of part 4 show it is always slightly pessimistic, and more so at low rates where ω_cT is no longer small; it also ignores the gain of the sinc factor and the small shift of the crossover frequency.

### 2.5 Discretising a controller: Tustin and matched pole-zero

A continuous controller C(s) must become a difference equation. Three common routes:

- **ZOH equivalent** (2.1): exact for a held input, natural for the *plant*, less so for a controller whose input (the sampled error) is not held.
- **Tustin (bilinear)** (16.06 Lecture 30; 6.003 Lecture 7 calls it the trapezoidal rule): substitute **s = (2/T)(z − 1)/(z + 1)**. It comes from the trapezoidal rule for integration, 1/s ≈ (T/2)(z + 1)/(z − 1). In state space, with W = I − AT/2: A_d = W⁻¹(I + AT/2), B_d = W⁻¹BT, C_d = CW⁻¹, D_d = D + CW⁻¹BT/2 (the realisation `cont2discrete(method="bilinear")` returns). The map sends the imaginary axis onto the unit circle exactly (Problem 2c), so a stable C(s) gives a stable C_d(z) at any T, and the gain at every frequency is reproduced, but at the wrong frequency: the continuous response at ω_a appears at **ω_d = (2/T)arctan(ω_aT/2)**, compressed toward the Nyquist frequency. **Prewarping** replaces 2/T by ω₀/tan(ω₀T/2), which makes the discrete response exact at one chosen frequency ω₀ (the crossover, or a notch).
- **Matched pole-zero**: map every pole and finite zero by z = e^{sT}, put each zero at infinity at z = −1 (the highest frequency the discrete system can represent), and scale so the DC gains agree. It preserves the dynamics' time constants exactly but not the frequency response.

For the lead compensator C(s) = (0.5s + 2)/(s/20 + 1) at 50 Hz, Tustin gives C_d(z) = (8.667z − 8.000)/(z − 0.6667) and matched pole-zero gives (8.576z − 7.917)/(z − 0.6703). At the crossover both are close to the continuous response: Tustin's gain is high by 0.11 % with a 0.01° phase error, matched pole-zero's phase is off by 0.22°; at 20 Hz the errors grow to 0.69 % and 1.34°. Prewarped Tustin is exact at ω_c at both rates. For a lead at a crossover this far below the Nyquist frequency the choice barely matters; for a notch near a structural mode it decides whether the filter works (Problem 3).

### 2.6 Where the model stops being valid

- **Fixed, synchronous timing.** Real schedulers jitter; τ varies from sample to sample. Design for the worst case τ = T (Problem 5).
- **No anti-aliasing filter.** The gyro and attitude signals pass through analog and digital filters (`P09.06`, `P09.07`), each adding phase lag at ω_c. A second-order Butterworth at 40 Hz costs 2.3° at 7.12 rad/s.
- **No quantisation.** Fixed-point arithmetic and sensor resolution add noise and, in loops with integrators, limit cycles.
- **Linear plant.** The elevator has rate and travel limits (`P09.11`, `P09.20`).
- **The rule (6)** assumes ω_cT ≪ 1. At 7 Hz with a full-sample delay, ω_cT = 1.02, the rule is 3.5° pessimistic and the loop is unstable anyway.

## 3. Build the computational model

Open `starter/main.py`. A matrix exponential, the pitch model and gains, the continuous frequency response, the hybrid simulation, the discrete margin and the sweep are given.

### Step 1: ZOH and delay

```python
def c2d_zoh(A, B, T):
    n, m = B.shape
    M = np.zeros((n + m, n + m))
    M[:n, :n], M[:n, n:] = A, B
    E = expm(M * T)
    return E[:n, :n], E[:n, n:]
```

`c2d_zoh_delay(A, B, T, tau)` calls `c2d_zoh` three times (over T, over τ and over T − τ, with the identity and zero for an empty interval), forms G₀ and G₁ from (3), and assembles (4). It raises `ValueError` for τ outside [0, T].

### Step 2: Tustin and matched pole-zero

`c2d_tustin` is four lines with `np.linalg.solve(W, ·)`; use `np.linalg.solve(W.T, C.T).T` for C W⁻¹. `tustin_tf(num, den, T, w_prewarp=None)` substitutes s = c(z − 1)/(z + 1) term by term: the coefficient a of s^k becomes a c^k (z − 1)^k (z + 1)^{n−k}, built with `np.poly` (which returns the monic polynomial with given roots) and `np.polymul`; normalise so the leading denominator coefficient is 1. `matched_tf` uses `np.roots`, maps, appends −1 for each zero at infinity, rebuilds with `np.poly` and scales by G(0)/G_d(1); it raises for a pole at s = 0.

### Step 3: z-domain tools and the loop

`dlsim` loops the difference equation and returns (y, x); accept a 1-D input for one channel. `impulse_by_long_division` pads the numerator to the denominator's length and runs h_k = b_k − Σ a_ih_{k−i}. `freqresp_z` is `np.polyval` at e^{jωT}; `spectral_radius` is `max(abs(eigvals))`; `z_to_s` is `np.log(z)/T`. `discrete_pitch_loop(T, tau, K)` closes u_k = −Kx_k + K_pθ_ref around (4): A_cl = Φ − Γ[K, 0], B_cl = ΓK_p.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-17-digital-control`. The tests:

1. **ZOH, the independent check**: against `cont2discrete(method="zoh")` and `scipy.linalg.expm` on random (2,1), (3,2), (4,1), (5,2) systems at four periods to 10⁻¹⁰, and against A⁻¹(e^{AT} − I)B; the double integrator gives [[1, T], [0, 1]] and (T²/2, T) exactly; the small-T limit is Euler.
2. **Step invariance**: the closed-loop pitch step by ZOH equals `scipy.signal.step` at the samples to 10⁻⁹; Tustin's differs by more than 10⁻⁴ (here 0.048 rad per rad of reference).
3. **Tustin**: state space against `cont2discrete(method="bilinear")` at three periods; transfer function against the same; G_d(e^{jω_dT}) = G(jω_a) with the warping law; prewarping exact at ω₀ and plain Tustin not.
4. **Matched pole-zero**: poles e^{pT}, zeros e^{qT}, DC gain, zeros at infinity at −1, refusal for an integrator.
5. **z tools**: long division equals `dimpulse` and the geometric series; `dlsim` superposition, zero input, time invariance; stable maps to |z| < 1 for ZOH and Tustin at every T; ln(z)/T recovers the continuous eigenvalues.
6. **Delay**: the τ = 0 and τ = T limits, G₀ + G₁ = B_d, refusal of τ > T; (4) equals `hybrid_simulation` to 10⁻¹² for τ/T = 0, 0.25, 0.5, 1.
7. **Pitch loop**: the discrete closed-loop poles approach the continuous ones as T → 0; damping and margin fall with T and τ; the rule within 2° at 20 Hz and 50 Hz; 7 Hz with a full-sample delay unstable, 50 Hz stable.

### Compare

`scipy.signal.cont2discrete` implements ZOH with `scipy.linalg.expm` of the same augmented matrix, and agrees with ours to 1.1 × 10⁻¹⁶ (A_d) and 7 × 10⁻¹⁸ (B_d) on the pitch plant at 50 Hz: they are the same formula. Its "bilinear" method is the generalised bilinear transform with α = ½ in the realisation of 2.5. It has no matched pole-zero method and no computational delay: the delay model (4) is something a flight-control engineer writes; The z-plane checks in the tests use `dlti` and `dimpulse` only as oracles.

Running `main()` prints (abridged):

```
continuous loop at the servo input: crossover 7.12 rad/s, phase margin 65.7 deg
  fs =   20 Hz: ZOH phase lag wc T/2 =  10.2 deg, one-sample delay adds  20.4 deg -> predicted margin  35.1 deg
  fs =   10 Hz: ZOH phase lag wc T/2 =  20.4 deg, one-sample delay adds  40.8 deg -> predicted margin   4.5 deg
ZOH at 50 Hz: eig(Ad) = [1.    +0.j     0.962 +0.0564j 0.962 -0.0564j 0.6703+0.j    ] = e^(lambda T): max error 1.1e-16
fs = 1 kHz, no delay: ln(z)/T = [-12.33+0.j   -5.25-5.5j  -5.25+5.5j  -0.81+0.j ] vs continuous [-12.41+0.j    -5.24-5.47j  -5.24+5.47j  -0.81+0.j  ]
20 Hz, tau = T/2: hybrid simulation vs difference equation at the samples, max |dx| = 2.5e-14
```

and the sample-rate table of part 4, and writes `digital_pitch.png` (θ step responses resolved between samples at 50, 20 and 10 Hz with a full-sample delay against the continuous one; the closed-loop poles in the z-plane at four rates).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Sample rate and delay sweep.** For f_s = 100, 50, 20, 10 and 7 Hz and τ/T = 0, 0.5, 1, record the spectral radius, the equivalent damping ratio of the least-damped complex pole (from ln z/T), the θ step overshoot and the phase margin of the exact sampled loop, and put the rule (6) beside it. Predict from (6) before running: the rates at which the margin falls below 45° and below 0.
2. **Continuous limit.** At f_s = 1 kHz with no delay, compare ln(z)/T of the closed-loop poles with the continuous poles. Predict the size of the difference from the half-sample delay of the hold.
3. **Controller discretisation.** For the lead compensator, compare Tustin, prewarped Tustin and matched pole-zero at 50 Hz and 20 Hz by the gain ratio and phase error at ω_c. Predict which errs most at 20 Hz.
4. **Step invariance and between samples.** Compare the ZOH and Tustin discretisations of the closed loop against `scipy.signal.step` at the samples; then run `hybrid_simulation` at 20 Hz with τ = T/2 and compare its samples with the difference equation.

Raw data format: CSV with header `fs_Hz, tau_over_T, spectral_radius, zeta_eq, overshoot, pm_deg, pm_rule_deg`; `method, fs_Hz, gain_ratio, phase_error_deg`; `check, max_abs_difference`. Each header line names the plant and gains.

What you should find:

| f_s (Hz) | τ/T | spectral radius | ζ_eq | overshoot | PM measured (°) | PM rule (°) |
|---|---|---|---|---|---|---|
| 100 | 0 / 0.5 / 1 | 0.9919 | 0.679 / 0.663 / 0.642 | 0 | 63.8 / 61.7 / 59.7 | 63.7 / 61.7 / 59.6 |
| 50 | 0 / 0.5 / 1 | 0.9839 | 0.662 / 0.617 / 0.559 | 0 | 61.7 / 57.6 / 53.6 | 61.7 / 57.6 / 53.5 |
| 20 | 0 / 0.5 / 1 | 0.960 | 0.585 / 0.432 / 0.293 | 0 | 55.8 / 45.6 / 35.5 | 55.5 / 45.3 / 35.1 |
| 10 | 0 / 0.5 / 1 | 0.921 / 0.920 / 0.976 | 0.427 / 0.184 / 0.033 | 0 / 0 / 0.211 | 46.2 / 25.9 / 5.9 | 45.3 / 24.9 / 4.5 |
| 7 | 0 / 0.5 / 1 | 0.889 / 0.945 / 1.076 | 0.319 / 0.054 / −0.082 | 0 / 0.132 / unstable | 38.7 / 9.6 / −18.2 | 36.6 / 7.4 / −21.7 |

Three readings. First, the **spectral radius is a poor guide**: it is set by the slow real pole near −0.81 s⁻¹ (e^{−0.81T}) until the loop is nearly unstable, while the damping and the margin fall steadily from the first row. Second, the **rule (6) is excellent**: within 0.4° at 20 Hz and above, and within 1.4° even at 10 Hz with a full-sample delay; it is always slightly pessimistic. Third, **delay is worse than slow sampling**: 20 Hz with a full-sample delay (35.5°) is worse than 10 Hz with none (46.2°). The overshoot column shows nothing until 10 Hz because the slow pole (the "θ zero" at −L_α) dominates the step; the damping column is the honest early warning. At 1 kHz the discrete poles are within 0.08 s⁻¹ of the continuous ones (−12.33 against −12.41 for the servo pole, −5.25 ± 5.50j against −5.24 ± 5.47j). Tustin's step response is off by up to 0.048 at the samples while ZOH is exact to 5 × 10⁻¹⁵, and the between-sample simulation agrees with the difference equation to 2.5 × 10⁻¹⁴. For a flown example of the same trade, read how the NASA F-8 digital fly-by-wire team chose its sample rates (Deets, TN D-7843 paper 2, "Design Synthesis and Analysis" and "Control Law Modes").

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| A_d at 50 Hz | e^{AT} | `c2d_zoh` | `cont2discrete`, `expm` | 1.1 × 10⁻¹⁶ | rounding |
| eig(A_d) | e^{λT}: 1, 0.962 ± 0.0564j, 0.6703 | same | — | 1.1 × 10⁻¹⁶ | rounding |
| step at the samples | step invariance: 0 | ZOH | `scipy.signal.step` | 5 × 10⁻¹⁵ | rounding |
| delay model, τ = T/2, 20 Hz | (3)–(4) | `discrete_pitch_loop` | `hybrid_simulation` | 2.5 × 10⁻¹⁴ | rounding |
| PM, 20 Hz, τ = T | 65.7 − 30.6 = 35.1° | 35.5° | — | 0.4° | first-order rule |
| PM, 10 Hz, τ = T | 4.5° | 5.9° | — | 1.4° | ω_cT = 0.71, not ≪ 1 |
| 7 Hz, τ = T | PM < 0 | spectral radius 1.076 | — | — | unstable |
| lead, Tustin at 20 Hz | warping | gain +0.69 %, phase +0.07° | prewarped: exact | — | (2/T)arctan law |

The discretisations are exact where the theory says they are and wrong by the predicted amount where it says they are not. The engineering conclusion: for this loop, **50 Hz with at most half a sample of computational delay keeps 57.6° of the 65.7°**; 20 Hz is acceptable only if the computation finishes early in the period; 10 Hz is not acceptable with any realistic delay. A common rule of thumb asks for a sample rate 20 to 40 times the closed-loop bandwidth; here ω_s/ω_c = 2π·50/7.12 = 44.

What to improve: add the anti-aliasing filter's phase (Problem 5) and the elevator rate limit; check jitter (Extension 3). What not to improve: the discretisation of the plant. ZOH is exact, so a finer integrator only adds cost.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3)–(4), the Tustin warping law and (5); `digital_pitch.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a sample rate and a computational-delay budget for the pitch loop, with the margin that remains and the effects (anti-aliasing filter, jitter) that will eat into it.

### Adversarial self-checks

- **What if the flight computer is twice as slow?** Halving f_s from 50 to 25 Hz with a full-sample delay takes the margin from 53.6° to 41.5° (Problem 4) and the damping ratio from 0.56 to 0.37: still stable, no longer comfortable.
- **What if T → 0?** The discrete poles approach e^{s_iT} → 1 and ln(z)/T approaches the continuous poles; the tests show the error shrinking more than three times from T = 4 ms to 1 ms.
- **Which assumption fails first?** Fixed timing. If the computation sometimes overruns into the next period, τ exceeds T and the augmented model of (4) no longer describes the loop; the flight software must then either skip an update (a jitter of a full sample) or delay all of them.
- **What measurement could falsify the model?** Log the elevator command and the attitude at a high rate in a ground test with the actual flight computer and inject a small sinusoid at ω_c: the measured phase of the loop at ω_c should match L(e^{jω_cT}) within a degree or two; a larger lag reveals an unmodelled filter or delay.
- **Does the result survive a change of method?** ZOH by our `expm`, by SciPy, and by a between-sample simulation agree to 10⁻¹⁴; the margin by the exact sampled loop and by the rule agree to 0.4° at 20 Hz.
- **Which decision changes if the delay were twice as uncertain?** Designing for τ = T instead of τ = T/2 at 50 Hz costs 4° (57.6° to 53.6°); at 20 Hz it costs 10° (45.6° to 35.5°). The faster computer makes the delay budget almost irrelevant, which is the strongest argument for it.

## Check yourself

1. Write the ZOH equivalent (A_d, B_d) of ẋ = Ax + Bu and the augmented matrix that produces both.
2. What does step invariance mean, and which discretisation has it?
3. State the stability condition for x_{k+1} = Mx_k and the map that takes continuous poles to discrete ones under ZOH.
4. Why does a zero-order hold behave like a delay of half a sample, and what phase does it cost at frequency ω?
5. A loop has crossover 10 rad/s and runs at 50 Hz with a full-sample computational delay. Estimate the phase-margin loss.
6. Write Tustin's substitution and the warping law. What does prewarping fix, and what does it not fix?
7. In matched pole-zero discretisation, where do zeros at infinity go, and why?
8. Why is the spectral radius of the pitch loop a poor indicator of its robustness?

## References

The sources, with the exact lectures, topics and sections to read, are in the **References** row at the top of this lesson: 16.06 Lectures 29–34 and 16.30 Topic 20 for sampling, emulation and the effective delay; 6.003 Lectures 5 and 7 for the z-transform and the pole mappings; the NASA F-8 digital fly-by-wire report for flight context; SciPy for the test oracles. The augmented-matrix derivation of the ZOH equivalent is original instruction (Partial GAP in `references/phase-09.md`).
