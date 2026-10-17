# P09.19  Nonlinear Control Basics

> Linear analysis tells you what happens near the operating point. A Lyapunov function tells you how far "near" reaches, a describing function tells you how big the oscillation of a saturated loop will be, and feedback linearisation tells you what it costs to make the nonlinearity disappear.

| | |
|---|---|
| **Learning outcomes** | 1. Solve the Lyapunov equation AᵀP + PA + Q = 0 from scratch, verify it against `scipy.linalg.solve_continuous_lyapunov` and its residual, and use V = xᵀPx to estimate the region of attraction of a nonlinear, torque-limited loop, checked by simulation. 2. Derive the describing functions of a saturation and of a relay with hysteresis, predict the amplitude and frequency of a limit cycle by harmonic balance G(jω)N(A) = −1, and compare the prediction with a simulation of the nonlinear loop. 3. Derive and implement a feedback-linearising control law for an inverted pendulum, show that the closed loop is exactly linear, and quantify what it costs in torque and how it degrades with model error. 4. Decide which of these tools to trust for a given loop and explain where each one fails. |
| **Prerequisites** | `P09.16` (Lyapunov equations in LQR), `P09.15` (state feedback), `P09.10` (Nyquist criterion), `P09.14` (Gramians by Lyapunov solves), `P02.06` (linearisation of nonlinear ODEs), `P02.09` (Fourier series), `P03.05` (RK4 and `solve_ivp`). Skill check: you can linearise ẍ = sin x about x = 0 and compute the first Fourier sine coefficient of a square wave. |
| **Engineering question** | An inverted-pendulum test stand (a vehicle balanced on its pivot) has a 6 N·m motor, too weak to hold it beyond 37.7°: from how large an initial tilt can a linear controller be *proven* to recover it, and how large in simulation? A fin servo under relay feedback, or with an over-driven saturating amplifier, oscillates: can the amplitude and frequency be predicted before the test? And can a controller that cancels gravity make the pendulum behave exactly linearly at 60°, and at what price? |
| **Assumptions and limits** | Pendulum: a rigid body on a frictionless-plus-viscous pivot, m l²θ̈ = mgl sin θ − bθ̇ + u, with representative m = 2 kg, l = 0.5 m, b = 0.1 N·m·s and a hard torque limit of 6 N·m; state measured exactly. Servo: the linear plant 30/(s(0.04s + 1)(0.01s + 1)) (the `P09.18` servo with a 10 ms amplifier lag) with one static nonlinearity in the loop. The describing function keeps only the first harmonic: it assumes the plant filters the higher harmonics (low-pass), and is an approximation, not a theorem. Lyapunov estimates are sufficient conditions and conservative. Feedback linearisation assumes the model is exact and the torque is unlimited. |
| **Analytical work** | The Lyapunov equation and its integral solution; V̇ = −xᵀQx along linear trajectories; Lyapunov's indirect method and why it is local; the sublevel-set estimate of the region of attraction; the first-harmonic describing functions of saturation and relay with hysteresis; harmonic balance and the stability of the predicted limit cycle; the feedback-linearising law u = ml²v − mgl sin θ + bθ̇ and its torque requirement. |
| **Computational deliverable** | `starter/main.py`: `lyapunov_solve`, `jacobian_fd`, `lyapunov_derivative`, `roa_estimate`, `df_saturation`, `df_relay`, `phase_crossover`, `predict_relay_limit_cycle`, `predict_relay_hysteresis_limit_cycle`, `predict_saturation_limit_cycle`, `oscillation_amplitude_frequency`, `feedback_linearising_torque` (the models, `tf_to_ss`, the two loop simulations, `first_harmonic_gain`, `converges`, `largest_recoverable_angle`, `simulate_pendulum` and `main` are given). Tests prove: the Lyapunov solution has residual below 10⁻¹⁰ and matches `solve_continuous_lyapunov` for n = 2, 3, 6, is positive definite for Hurwitz A and indefinite for an unstable A; dV/dt = −xᵀQx along linear trajectories; the finite-difference Jacobian matches the analytic linearisation at three angles; the upright pendulum is a saddle and the controller places its poles; dV/dt < 0 throughout the estimated region, 25 random starts inside it converge, and the estimate lies inside the static limit; more torque gives a larger region; both describing functions match a numerical Fourier coefficient; the relay, relay-with-hysteresis and saturation limit-cycle predictions match simulation within 6 % in amplitude and 4 % in frequency; below the critical gain no cycle is predicted and none appears; feedback linearisation cancels the nonlinearity exactly and the 60° response equals the linear closed form to 10⁻⁸. |
| **Practical activity** | **computational experiment**: region-of-attraction estimates for several Lyapunov weights against the simulated recovery limit; limit-cycle prediction against simulation for relay, relay with hysteresis and saturation at several gains; feedback linearisation against linear control from large angles, with model error and with the torque limit. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `nonlinear.png`, the derivations of the Lyapunov integral solution, of both describing functions and of the feedback-linearising law, the three experiment tables, and a half-page note stating which tool you would use to certify the test stand and the servo, with the conservatism or error of each. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the Lyapunov and describing-function results are derived; all public tests pass; the gap between the estimated and simulated regions of attraction and the describing-function error are quantified and explained; the note separates what is proven from what is predicted. |
| **References** | `@mit-ocw-16-30-feedback-control-systems` Topic 22 "Analysis of nonlinear systems" (https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec22/): Lyapunov Stability Analysis; Lyapunov Stability Theorem; Example 1: Pendulum; Example 2: Linear System; Example 4: Local Region; Example 5: Saturation; Invariance Principle (required-free): V > 0 with V̇ < 0, the Lyapunov equation AᵀP + PA = −Q, and a sublevel set of V as a region-of-attraction estimate (§2.1–2.2). Topic 21 "Systems with nonlinear functions" (.../mit16_30f10_lec21/): Describing Function Analysis; Saturation Nonlinearity; Odd Nonlinearities with Memory; Limit Cycle Analysis; DF Analysis Example; Limit Cycle Stability (required-free): N(A) of saturation and of a relay with hysteresis, and limit cycles where G(jω) = −1/N(A) (§2.3). Topic 5 "Introduction to state-space models": Linearization; Stability of LTI Systems (supplementary-free review): Lyapunov's indirect method. `@mit-ocw-6-243j-dynamics-of-nonlinear-systems` Lectures 5 "Lyapunov Functions and Storage Functions", 7 "Finding Lyapunov Functions", 8 "Local Behavior at Equilibria", 13 "Feedback Linearization" (https://ocw.mit.edu/courses/6-243j-dynamics-of-nonlinear-systems-fall-2003/resources/lec13_6243_2003/) (supplementary-free; graduate and proof-oriented, beyond this lesson's depth). `@scipy-userguide` API reference "Linear algebra (scipy.linalg)": Matrix Equation Solvers (`solve_continuous_lyapunov`) (supplementary-free): the test oracle; `scipy.integrate.solve_ivp` is the simulator. Partial GAP (original instruction): the feedback linearisation of the inverted pendulum at this phase's level (§2.4) is written here; the 6.243J lecture on it is graduate level and its PDF was not checked. |
| **Extensions** | ★ Replace the quadratic Lyapunov function by xᵀPx + κ(1 − cos θ) and search κ for a larger proven region. ★ Add the third harmonic to the relay describing function analysis and estimate how much of the 4 % amplitude error it explains. ★ Add a dead zone (gear backlash) to the servo and predict its limit cycle with the dead-zone describing function. |
| **Time** | unmeasured |

## 1. The engineering problem

Every real control loop is nonlinear somewhere. Actuators saturate: a fin servo's amplifier cannot exceed its supply voltage, a reaction wheel its torque. Some actuators are nonlinear by design: an on-off valve on a cold-gas thruster is a relay. And some plants are nonlinear in their physics: gravity acts on a pendulum through sin θ, not θ. Linear tools, used inside their range, are excellent; the danger is using them outside it without knowing where the range ends.

This lesson takes three loops and asks three questions. First, an inverted-pendulum test stand, a 2 kg body on a pivot 0.5 m below its centre of mass (a stand-in for a vehicle balanced on its thrust vector), held by a 6 N·m motor. A linear state-feedback controller stabilises it near upright, but gravity's torque mgl sin θ reaches 6 N·m at 37.7°, beyond which no controller can hold it. *From which initial states can we prove that the controller recovers?* Second, the fin servo of `P09.18`: under relay (bang-bang) feedback, or with a saturating amplifier and too much gain, it settles into a sustained oscillation, a **limit cycle**. *What amplitude and frequency, before the test?* Third: *can the controller simply cancel gravity*, and what does that cost?

Success criteria, written first: the Lyapunov solver agrees with SciPy to 10⁻⁹ with residual below 10⁻¹⁰; every initial state inside the proven region converges in simulation; the predicted limit cycles agree with simulation within 6 % in amplitude and 4 % in frequency; the feedback-linearised response matches the linear closed form to 10⁻⁸ rad at 60°.

## 2. Understand and derive

### 2.1 Lyapunov functions and the Lyapunov equation

A **Lyapunov function** for an equilibrium x = 0 of ẋ = f(x) is a function V(x) that is positive except at 0 and decreases along every trajectory, V̇ = ∇V·f < 0 for x ≠ 0. If one exists on a region, trajectories there can only go "downhill" in V and must end at 0. No solution of the ODE is needed: that is the point.

For a linear system ẋ = Ax, try V = xᵀPx with P symmetric. Then V̇ = xᵀ(AᵀP + PA)x, and asking for V̇ = −xᵀQx with a chosen Q > 0 gives the **Lyapunov equation** (16.30 Topic 22, "Example 2: Linear System")

**AᵀP + PA + Q = 0.**   (1)

If A is Hurwitz, (1) has the unique solution P = ∫₀^∞ e^{Aᵀt}Qe^{At}dt, which is positive definite (Problem 2a); conversely, P > 0 solving (1) with Q > 0 proves A Hurwitz. `lyapunov_solve` solves (1) as the n² × n² linear system (I ⊗ Aᵀ + Aᵀ ⊗ I) vec P = −vec Q, as in `P09.16`.

For a nonlinear system with Jacobian A = ∂f/∂x at 0 (Lyapunov's indirect method, `P02.06`), the same V proves *local* stability: near 0 the nonlinear terms are small compared with −xᵀQx. How near is the question of 2.2.

**The pendulum.** With θ from upright, m l²θ̈ = mgl sin θ − bθ̇ + u. The Jacobian at upright is [[0, 1], [g/l, −b/(ml²)]] = [[0, 1], [19.62, −0.2]], eigenvalues +4.331 and −4.531 s⁻¹: a saddle, unstable. The state feedback u = −k₁θ − k₂θ̇ with k₁ = ml²(ω_n² + g/l) = 17.810 N·m/rad and k₂ = 2ζω_nml² − b = 2.700 N·m·s/rad places the linearised poles at ω_n = 4 rad/s, ζ = 0.7: −2.8 ± 2.857j s⁻¹. With Q = diag(1, 0.1), (1) gives P = [[0.4071, 0.0312], [0.0312, 0.0145]] with residual 2.8 × 10⁻¹⁷.

### 2.2 The region of attraction

The real loop is u = sat(−Kx, 6 N·m) and sin θ: nonlinear twice. Along its trajectories V̇(x) = 2xᵀPf(x) is negative near 0 and becomes positive somewhere further out. Let c* be the largest level such that V̇ < 0 at every x ≠ 0 with V(x) < c*. Then the ellipse {V < c*} is **positively invariant** (a trajectory cannot cross a level set of V outward while V decreases) and every trajectory starting in it converges to 0: it is a proven inner estimate of the **region of attraction** (16.30 Topic 22, "Example 4: Local Region" and "Example 5: Saturation", make the same argument).

`roa_estimate` finds c* by scanning: in coordinates where V = r² (x = rP^{−1/2}(cos φ, sin φ)), march outward along 360 directions and record the first radius at which V̇ ≥ 0; c* is the square of the smallest. For the stand, **c* = 0.0689**: the ellipse reaches ±25.8° in θ and ±23.6° when starting at rest. A hundred random starts inside it (seed 19) all converge. Simulation from rest recovers up to **37.71°**, which is the physical limit arcsin(6/9.81) = 37.7°: beyond it gravity's torque exceeds the motor's at rest. The guarantee is valid and conservative: it proves 23.6° of the 37.7° that are actually achievable. A quadratic V is an ellipse, and the true region is not an ellipse.

### 2.3 Describing functions

Put a static nonlinearity y = φ(e) in a loop with a linear plant G, e = −Gy. Suppose the loop oscillates with e ≈ A sin ωt. Then y is periodic with the same period; keep only its fundamental, b₁ sin ωt + a₁ cos ωt, and define the **describing function**

N(A) = (b₁ + ja₁)/A,   (2)

the "gain" of the nonlinearity for a sinusoid of amplitude A (16.30 Topic 21, "Describing Function Analysis"). If the plant filters the higher harmonics (|G(3jω)| ≪ |G(jω)|), the loop closes on the fundamental alone, and a self-sustained oscillation requires **G(jω)N(A) = −1** (harmonic balance): G(jω) = −1/N(A).

- **Saturation** of slope k and linear range δ (output k·clip(e, −δ, δ)): N = k for A ≤ δ, and for A > δ (Problem 2b)

  N(A) = (2k/π)[arcsin(δ/A) + (δ/A)√(1 − (δ/A)²)],   (3)

  real (no phase), falling from k to 4kδ/(πA) as A grows: the bigger the swing, the smaller the effective gain.
- **Relay** ±M: N = 4M/(πA) (the first harmonic of a square wave). With **hysteresis** h (switch up when e rises through +h, down when it falls through −h), N(A) = (4M/(πA))[√(1 − (h/A)²) − jh/A], and −1/N(A) = −(π/(4M))(√(A² − h²) + jh) is a horizontal line at Im = −πh/(4M) (Problem 2c).

Both describing functions match a direct numerical Fourier coefficient in the tests (`first_harmonic_gain`).

**The servo.** G(s) = 30/(s(0.04s + 1)(0.01s + 1)) has its phase crossover at ω₁₈₀ = 1/√(0.04 × 0.01) = 50 rad/s, where |G| = 30 × 0.0004/0.05 = 0.24 rad/V.

- Relay ±1 V: −1/N is the negative real axis, so ω = 50 rad/s and A = 4M|G(jω₁₈₀)|/π = **0.3056 rad**. Simulation (`simulate_relay_loop`, switching instants located exactly by event detection): **0.3194 rad at 48.55 rad/s**. The prediction is 4.3 % low in amplitude and 3.0 % high in frequency, the price of neglecting the square wave's harmonics.
- Relay with h = 0.02 rad: the intersection with the horizontal line gives **0.3517 rad at 46.57 rad/s**; simulation **0.3655 rad at 45.36 rad/s**. Hysteresis makes the cycle bigger and slower, as predicted.
- Saturating amplifier, limit 1 V, gain twice the critical gain 1/0.24 = 4.167 V/rad: the linear loop is unstable, but N(A) falls as A grows until N(A)|G(jω₁₈₀)| = 1. Prediction **0.2970 rad at 50.00 rad/s**; simulation **0.3046 rad at 49.18 rad/s**. Below the critical gain no solution exists and the simulated loop decays.

A predicted limit cycle is **stable** (it attracts nearby oscillations) if a slightly larger amplitude makes the loop "more stable": for these monotonically falling N(A), a larger A means a smaller gain, the oscillation shrinks back. All three cycles here are stable, which is why simulation finds them.

### 2.4 Feedback linearisation

If the model is known, the controller can cancel the nonlinearity. For the pendulum choose

**u = ml²v − mgl sin θ + bθ̇,**   (4)

then the plant becomes θ̈ = v exactly, for every θ. Choosing v = −ω_n²θ − 2ζω_nθ̇ makes the closed loop *exactly* the linear system with poles −2.8 ± 2.857j, at 60° as at 1°. The simulation from 60° matches the linear closed form to 1.6 × 10⁻¹¹ rad. Three costs:

- **Torque.** (4) must supply the full gravity torque: 16.87 N·m peak from 60° (the plain linear controller with the same poles asks 18.65 N·m and overshoots 3.30° against 2.76°). With a 6 N·m motor, (4) is clipped and the exact linearity is gone (Problem 6).
- **Model knowledge.** The cancellation is exact only if m, l and b are. A 10 % mass error leaves 10 % of gravity uncancelled (Problem 6: a 4.5° deviation from the linear response from 60°).
- **The physical limit stays.** No controller holds the stand beyond 37.7° at rest with 6 N·m.

### 2.5 Where the models stop being valid

- **Describing functions** assume one dominant frequency and a low-pass plant. A plant with a resonance near 3ω, a nonlinearity with strong even harmonics (an asymmetric one), or two interacting nonlinearities break it; so does a plant whose phase never reaches the −1/N locus (Problem 3a: an ideal relay on a plant with phase above −180° chatters instead).
- **Quadratic Lyapunov estimates** are always ellipses; the true region of attraction of the stand is a strip bounded by the stable manifolds of the saddles near ±37.7°.
- **The pendulum model** has no motor dynamics, no sensor delay and no pivot friction beyond viscous.
- **The ideal relay** switches instantly; a real valve has a minimum on-time, which changes the cycle at small amplitudes.

## 3. Build the computational model

Open `starter/main.py`. The models, the two loop simulations (event-driven for the relay, `solve_ivp` for saturation), the numerical Fourier coefficient and the pendulum simulator are given.

### Step 1: Lyapunov and the region of attraction

`lyapunov_solve` is the Kronecker solve of `P09.16`. `jacobian_fd` uses central differences with a step of 10⁻⁶ max(1, |x_i|). `lyapunov_derivative` is 2xᵀPf(x). `roa_estimate`:

```python
lam, U = np.linalg.eigh(P)
Pmh = U @ np.diag(1 / np.sqrt(lam)) @ U.T            # P^(-1/2): V(r Pmh d) = r^2
rs = np.linspace(r_max / n_r, r_max, n_r)
r_first = r_max
for phi in np.linspace(0, 2 * np.pi, n_dir, endpoint=False):
    d = Pmh @ np.array([np.cos(phi), np.sin(phi)])
    for j, r in enumerate(rs):
        if r >= r_first:
            break                                       # cannot lower the minimum further
        if lyapunov_derivative(f, P, r * d) >= 0:
            r_first = rs[j - 1] if j else 0.0           # the last radius that was still decreasing
            break
return r_first ** 2
```

### Step 2: describing functions and harmonic balance

`df_saturation` and `df_relay` are (3) and the relay formula, vectorised. `phase_crossover` brackets the −180° crossing of the unwrapped phase on a log grid and bisects. `predict_relay_limit_cycle` returns 4M|G(jω₁₈₀)|/π and ω₁₈₀. `predict_relay_hysteresis_limit_cycle` bisects Im G(jω) = −πh/(4M) on the bracket with Re G < 0 and returns A = √(h² + (4M Re G/π)²). `predict_saturation_limit_cycle` raises if gain·|G(jω₁₈₀)| ≤ 1 and otherwise bisects N(A)|G| = 1. `oscillation_amplitude_frequency` takes half the peak-to-peak after a settling time and 2π over the mean spacing of upward crossings of the mid-level.

### Step 3: feedback linearisation

`feedback_linearising_torque(x, v, p)` is (4), one line.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-19-nonlinear-control-basics`. The tests:

1. **Lyapunov, the independent check**: residual below 10⁻¹⁰‖Q‖ and agreement with `solve_continuous_lyapunov(A.T, -Q)` for n = 2, 3, 6; P > 0 for Hurwitz A; a singular case raises; an unstable A gives an indefinite P; dV/dt by central differences along e^{At}x₀ equals −xᵀQx to 10⁻⁶.
2. **Linearisation**: the Jacobian against the analytic one at θ = 0, 0.7 and π; the upright saddle; the closed-loop poles at the roots of s² + 5.6s + 16.
3. **Region of attraction**: V̇ < 0 at 216 points inside the estimate; 25 random starts inside converge; the at-rest extent is below the static limit; c* grows with the torque limit (3, 6, 12 N·m).
4. **Describing functions**: saturation at five amplitudes and relay with h = 0, 0.1, 0.3 against the numerical Fourier coefficient; monotonicity and the large-A limit 4kδ/(πA).
5. **Limit cycles, the independent check**: the relay, relay-with-hysteresis and saturation predictions against event-driven or adaptive simulation (amplitude within 6 %, frequency within 4 %); harmonic balance holds exactly at the predicted point; hysteresis gives a larger, slower cycle; at 0.8 of the critical gain no cycle is predicted, and at 0.5 the simulated loop decays.
6. **Feedback linearisation**: θ̈ = v for random states and inputs; the 60° response equals e^{A_lin t}x₀ to 10⁻⁸; the linear controller's response does not.

### Compare

`scipy.linalg.solve_continuous_lyapunov` uses the Bartels–Stewart algorithm (Schur forms, O(n³)); on these 2 × 2 to 6 × 6 problems it agrees with the Kronecker solve to 10⁻⁹. There is no library describing-function tool in SciPy; the oracle is the simulation, run with `solve_ivp` (RK45, relative tolerance 10⁻⁹–10⁻¹⁰) and, for the relay, with exact event location at each switch so that the discontinuity is never smoothed. The 3–4 % disagreement between prediction and simulation is far larger than the integration error at these tolerances, so it is the method's own error: the neglected harmonics.

`main()` prints:

```
pendulum: m = 2.0 kg, l = 0.5 m, b = 0.1 N m s; m g l = 9.81 N m, torque limit 6.0 N m -> largest holdable angle 37.7 deg
region-of-attraction estimate: V < 0.0689, theta extent +-25.8 deg; theta at rest inside it: +-23.6 deg
   100 random starts inside the estimate (seed 19): 100 converge
   from rest, simulation recovers up to 37.71 deg (bisection): the estimate is conservative
relay +-1 V on G = 30/(s(0.04s+1)(0.01s+1)): predicted A = 0.3056 rad at 50.00 rad/s; simulated A = 0.3194 rad at 48.55 rad/s (errors -4.3 %, +3.0 %)
saturating amplifier, gain 8.333 V/rad (= 2 x critical 4.167), limit 1.0 V: predicted A = 0.2970 rad at 50.00 rad/s; simulated A = 0.3046 rad at 49.18 rad/s
feedback linearisation from 60 deg (no torque limit): max |theta - linear closed form| = 1.6e-11 rad; peak torque 16.87 N m
```

and writes `nonlinear.png` (the phase portrait of the torque-limited stand with the proven ellipse; the simulated relay and saturation cycles against the predicted amplitude; feedback-linearised, linear-closed-form and linear-controller responses from 60°).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Proven versus simulated region.** For Q = diag(1, 0.1), diag(1, 0.15) and diag(1, 1), compute c* and the angle at rest inside the ellipse; seed 100 random starts inside each (seed 19) and count convergences; find the simulated recovery limit from rest by bisection. Predict: every start inside converges; the proven angle is well below 37.7°; changing Q changes the ellipse's shape more than its reach.
2. **Limit cycles.** Relay with h = 0 and 0.02 rad; saturating amplifier at 1.5, 2 and 2.25 times the critical gain. Record predicted and simulated amplitude and frequency. Predict the sign of the error (the harmonics are neglected).
3. **Feedback linearisation.** From 60° at rest, with the exact model, with 10 % and 20 % more mass than modelled, and with the 6 N·m limit; record the largest deviation from the linear response and the largest recoverable angle.

Raw data format: CSV with headers `Q11, Q22, c_star, theta_rest_deg, inside_converged_of_100, sim_limit_deg`; `case, A_pred, w_pred, A_sim, w_sim, err_A_pct, err_w_pct`; `case, max_dev_deg, peak_torque_Nm, recover_limit_deg`.

What you should find. The three Lyapunov weights give c* = 0.0689, 0.0812 and 0.2525 with at-rest reach 23.6°, 23.6° and 22.1°: the level c* changes a lot but the reach hardly at all, because the ellipse is limited by the saturated corner of the state space whatever its shape. All 300 starts inside converge; the simulated limit is 37.7° in every case. The limit cycles are all predicted within 5 %: relay 0.3056 against 0.3194 rad (−4.3 %), with hysteresis 0.3517 against 0.3655 rad, saturation at twice critical 0.2970 against 0.3046 rad, and at 2.25 times critical 0.2990 against 0.3075 rad (−2.8 %, frequency +1.9 %). The prediction is always low in amplitude and high in frequency. Feedback linearisation is exact with the exact model; a 10 % mass error causes a 4.5° deviation and 20 % causes 9.3°, but the loop still converges (the feedback v keeps it stable). With the 6 N·m limit both the clipped feedback-linearising controller and the clipped linear controller recover from 37.70°: the physical limit, not the control law, decides.

## 5. Validate and improve

| quantity | theory (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| Lyapunov P, residual | 0 | 2.8 × 10⁻¹⁷ | `solve_continuous_lyapunov`: 10⁻⁹ rel | — | rounding |
| region at rest | ≤ 37.7° (physics) | proven 23.6° | simulated 37.71° | 14° | ellipse vs true region |
| relay amplitude | DF 0.3056 rad | — | simulated 0.3194 rad | −4.3 % | harmonics |
| relay frequency | DF 50.00 rad/s | — | 48.55 rad/s | +3.0 % | harmonics |
| saturation, 2 × critical | DF 0.2970 rad, 50.00 rad/s | — | 0.3046 rad, 49.18 rad/s | −2.5 %, +1.7 % | harmonics |
| FL response from 60° | e^{A_lin t}x₀ | `solve_ivp` | — | 1.6 × 10⁻¹¹ rad | integrator |
| static limit | arcsin(6/9.81) = 37.71° | bisection 37.71° | — | < 0.02° | bisection |

Each tool does what its theory claims. The Lyapunov estimate is a proof and is 37 % short of the truth; the describing function is a prediction and is within 5 %; feedback linearisation is exact and fragile. For certifying the stand, the honest statement is "proven to recover from 23.6° at rest, demonstrated in simulation to 37.7°". For the servo, the describing function is good enough to size a relay test or to set an amplifier limit, and the simulation is the check.

What to improve: a non-quadratic Lyapunov function (Extension 1) or sum-of-squares methods enlarge the guaranteed region; the third harmonic (Extension 2) shrinks the describing-function error. What not to improve: the region-of-attraction scan resolution. The estimate is limited by the shape of V, not by the grid.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the Lyapunov integral solution, (3), the hysteresis relay describing function and (4); `nonlinear.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating what you can certify about the stand and predict about the servo, with the conservatism and error of each tool.

### Adversarial self-checks

- **What if the motor is twice as strong?** The physical limit becomes arcsin(12/9.81): no limit at all for a static hold (12 > 9.81 N·m); the proven region grows (test: c* increases from 6 to 12 N·m).
- **What if the relay has no hysteresis on a second-order plant?** The phase of K/(s(τs + 1)) never reaches −180°, so no limit cycle is predicted; the simulation chatters instead (more than 2000 switches in 1.04 s, Problem 3). The describing function is silent about chattering.
- **Which assumption fails first?** For the servo: the ideal relay (valves and drivers have minimum on-times and delays, which add phase and enlarge the cycle). For the stand: the exactly known model in feedback linearisation.
- **What measurement could falsify the prediction?** A relay test on the real servo: an amplitude more than about 10 % from 0.31 rad, or a frequency far from 49 rad/s, means the linear model (or the amplifier lag) is wrong; this is exactly how relay autotuning identifies ω₁₈₀.
- **Does the result survive a change of method?** The Kronecker and Bartels–Stewart Lyapunov solutions agree to 10⁻⁹; the describing-function prediction and the simulation agree to 5 % on every cycle tried, and the direction of the error (amplitude low, frequency high) is the same each time, which is what a systematic, harmonic-neglect error should look like.
- **Which decision changes if the mass is uncertain by ±20 %?** Feedback linearisation then deviates by up to 9.3° from its design response; a robust linear design (`P09.18`) with the same poles may be preferable, since it does not depend on cancelling gravity exactly.

## Check yourself

1. Write the Lyapunov equation and say what its solution proves when Q > 0 and P > 0.
2. Why does a Lyapunov function for the linearisation prove only local stability of the nonlinear system?
3. What makes the set {V < c*} a valid estimate of the region of attraction?
4. Define the describing function and state the harmonic-balance condition.
5. Give the describing function of an ideal relay ±M. For the servo with |G(jω₁₈₀)| = 0.24 rad/V and M = 2 V, what amplitude is predicted?
6. Why does a saturating amplifier turn an unstable linear loop into a limit cycle rather than divergence?
7. Write the feedback-linearising torque for the pendulum and name two things it requires.
8. Name one situation in which the describing-function prediction should not be trusted.

## References

The sources, with the exact topics and sections to read, are in the **References** row at the top of this lesson: 16.30 Topic 22 for Lyapunov functions, the Lyapunov equation and regions of attraction; 16.30 Topic 21 for describing functions and limit cycles; 6.243J as graduate-level further reading; SciPy for the oracle and the simulator. Feedback linearisation at this level is original instruction (Partial GAP in `references/phase-09.md`).
