# P09.11  PID Design & Tuning

> A PID controller is three numbers and a clip function. The three numbers come from a model; the clip function decides whether they still mean anything when the elevator hits its stop.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the transfer function of a PID with a filtered derivative, Kp + Ki/s + Kd s/(Tf s + 1), and of the two-degree-of-freedom closed loop T = P C_r/(1 + P C) with setpoint weight b and derivative on the measurement; show that b and the derivative input move only the zeros. 2. Implement and test a discrete PID (backward-Euler derivative filter, forward-Euler integrator, output limit) and two anti-windup schemes, clamping and back-calculation, on a pitch-attitude loop with an elevator servo. 3. Tune it two ways, by the Ziegler–Nichols ultimate-gain rules (Ku and Pu from the j-omega crossing) and by assigning three closed-loop poles exactly, and verify the sampled linear loop against the transfer-function prediction. 4. Decide which design to fly from overshoot, elevator demand, phase margin and behaviour under saturation, and say where the linear model stops being trustworthy. |
| **Prerequisites** | `P09.03` (transfer functions, block-diagram reduction), `P09.04` (step-response metrics, effect of a zero), `P09.08` (system type, steady-state error, disturbance rejection), `P09.09` (root locus, j-omega crossing by Routh), `P09.10` (gain and phase margins), `P02.05` (matrix exponential, zero-order-hold discretisation), `P02.08` (Laplace transform). Skill check: you can multiply polynomials with `np.polymul` and find the closed-loop poles of 1 + K G(s) = 0. |
| **Engineering question** | A small fixed-wing aircraft needs a pitch-attitude hold: elevator command from pitch error, through a servo with a 0.05 s time constant and a ±15° limit. Which PID gains give a well-damped 10° pitch step without demanding more elevator than exists, and what happens to the integrator during a 20° pull-up when the elevator saturates? Answer for the representative model of this lesson: the Ziegler–Nichols gains overshoot 64.7 % and ask for 32° of elevator; gains that place closed-loop poles at −4 ± 4j and −1.5 with setpoint weight b = 0 overshoot 10.6 % with 8.5° peak elevator; under saturation the unprotected integrator turns a 41 % overshoot into 68 %, and clamping anti-windup brings it to 13 %. |
| **Assumptions and limits** | Linear short-period dynamics about one trimmed flight condition (constant airspeed, small angles), Z_δ = 0, a first-order elevator servo with a pure magnitude limit (no rate limit), a perfect pitch measurement sampled at 1 kHz, no turbulence. The numbers are a representative small-aircraft model, not any specific aircraft's data. The model is wrong at large angle of attack (stall), at large pitch excursions where airspeed changes (the phugoid is ignored), and whenever the servo rate limit, not its position limit, is the binding constraint. |
| **Analytical work** | θ/δ_cmd from the short-period equations; the PID polynomial pair; T = P C_r/(1 + P C) and its invariant denominator; unit DC gain from integral action; the ultimate-gain condition Im[den(jω) num(−jω)] = 0; the Ziegler–Nichols rules; the affine dependence of the characteristic polynomial on (Kp, Ki, Kd) and the 3 × 3 placement solve; the PD steady error −d/Kp under a trim moment; the discrete PID update and the two anti-windup laws. |
| **Computational deliverable** | `starter/main.py`: `pitch_plant_tf`, `pid_tf`, `closed_loop_tf`, `ultimate_gain`, `ziegler_nichols`, `pid_place_three_poles`, `pid_step`, `simulate_pitch_loop` (`pitch_plant`, `pitch_state_space`, `expm`, `discretize_zoh`, `tf_step`, `step_metrics`, `main` are given). Tests prove: the plant TF equals C(sI − A)⁻¹B pointwise and `scipy.signal`'s realisation; the PID and closed-loop TFs equal their defining expressions for four (b, derivative-input) choices with unit DC gain; **the sampled linear loop without saturation matches `scipy.signal.step` of the closed-loop TF** to 0.6 % of the step, with error halving as dt halves; superposition; Ku puts roots on the j-omega axis (and Ku = 6 for 1/(s(s+1)(s+2))); the ZN rules; three-pole placement makes the chosen poles exact roots; the discrete PID's integration, derivative kick and clamping; anti-windup changes nothing without saturation and **reduces overshoot under saturation** by more than 40 %; integral action removes a trim disturbance. |
| **Practical activity** | **computational experiment**: tuning-method comparison (Ziegler–Nichols against pole placement, b = 1 against b = 0), sample-time convergence of the discrete loop, and an anti-windup study with the elevator limit binding. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pid_pitch.png`, the hand derivation of T(s) and of the placement equations, the tuning and anti-windup tables, and a half-page recommendation of a gain set and anti-windup scheme with its overshoot, peak elevator and phase margin. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: T(s) and the three-pole equations are derived; all public tests pass; the note explains why b = 0 cuts the overshoot without moving a single closed-loop pole, and why the unprotected integrator overshoots under saturation, with numbers. |
| **References** | `@mit-ocw-6-003-signals-and-systems` course lecture notes (required-free, course level: feedback and the Laplace-domain view of LTI systems). `@scipy-userguide` "Linear Algebra (scipy.linalg)" (supplementary-free: the oracle routines used in the tests). `@mit-ocw-18-06-linear-algebra` Lecture 23 "Differential equations and exp(At)" (supplementary-free: the exponential behind the zero-order-hold simulation). The Phase 09 reference map (`references/phase-09.md`) did not exist when this lesson was written, so no section of a control text is mapped yet: the PID, compensation, state-space and observer material here is original instruction (GAP), to be linked to exact sections when the map is published. |
| **Extensions** | ★ Add a servo rate limit (for example 60°/s) and show that it, not the position limit, binds for fast steps; add rate-limit-aware anti-windup. ★ Replace the backward-Euler derivative filter by the Tustin (bilinear) form and compare the sample-time convergence. ★ Close an altitude-hold outer loop around the tuned pitch loop and find how slow the outer loop must be for the inner loop to look like unity. |
| **Time** | unmeasured |

## 1. The engineering problem

Most flight-control loops flying today are PID loops or close relatives. An autopilot's pitch-attitude hold takes the difference between commanded and measured pitch, multiplies it by a gain (proportional), adds its running integral (so that a constant trim moment is eventually cancelled) and its rate of change (to add damping), and sends the sum to the elevator servo. Three gains. The difficulty is not the formula but the numbers in it, and what happens at the edges of the linear world: the elevator can only deflect so far, and once it is against its stop the integrator keeps integrating an error the actuator cannot remove. When the aircraft finally reaches the commanded attitude, the integrator holds a large stored command that drives the aircraft well past it. This is **integrator windup**, and it has caused real overshoots in real aircraft and spacecraft.

This lesson writes a pitch-attitude PID as flight software would: a discrete update every millisecond, with a filtered derivative, an output limit and an anti-windup rule. The plant is a representative short-period model of a light aircraft or small UAV, an integrator from pitch rate to pitch angle, and a first-order elevator servo. The gains are chosen two ways: by the Ziegler–Nichols rules, a 1942 recipe from process control that needs only one experiment (or one root-locus calculation), and by **placing three closed-loop poles exactly**, which needs the model. Both designs are then checked: first in the linear regime against the transfer-function prediction (the code and the algebra must agree), then with the elevator limit binding (where only simulation can say what happens).

Success criteria, written first: the sampled linear loop agrees with the closed-loop transfer function to better than 1 % of the step at a 1 ms sample time, and the discrepancy halves when the sample time halves; the placed gains put −4 ± 4j and −1.5 among the closed-loop roots to 10⁻⁶; for a 10° step the recommended design overshoots less than 15 % with less than 15° of elevator; under a 20° step with a 15° elevator limit, anti-windup cuts the overshoot of the unprotected loop by more than 40 %.

## 2. Understand and derive

### 2.1 The plant

Write the short-period equations (`P02.06` linearisation, `P09.03` block diagrams) for angle of attack α, pitch rate q, pitch angle θ and elevator δ, with δ > 0 trailing edge up (nose-up moment):

α̇ = Z_α α + q,  q̇ = M_α α + M_q q + M_δ δ,  θ̇ = q,  δ̇ = ω_s(δ_cmd − δ).   (1)

Z_α (1/s) is the lift-curve slope term L_α/V with its sign; M_α (1/s²) is the static stability (negative is stable); M_q (1/s) is pitch damping; M_δ (1/s²) is elevator power; ω_s is the servo bandwidth. The representative numbers are Z_α = −1.2 s⁻¹, M_α = −6 s⁻², M_q = −2.5 s⁻¹, M_δ = 12 s⁻², ω_s = 20 rad/s (a 0.05 s servo), elevator limit ±15°. Eliminating α: from the first equation q = (s − Z_α)α, so the second gives α[s² − (Z_α + M_q)s + Z_αM_q − M_α] = M_δ δ. Call the bracket Δ(s). Then q/δ = M_δ(s − Z_α)/Δ, θ = q/s, and with the servo

**θ/δ_cmd = M_δ ω_s (s − Z_α) / [s Δ(s) (s + ω_s)] = 240(s + 1.2) / [s (s² + 3.7s + 9)(s + 20)].**   (2)

Δ has ω_n = 3.000 rad/s and ζ = 0.617 (poles −1.85 ± 2.362j): the short-period mode. The integrator at s = 0 is the attitude kinematics. The zero at s = Z_α = −1.2 is physical: a nose-up elevator first raises the nose (q), and only as α builds does lift curve the flight path. Every pitch-attitude loop carries this slow zero, and §2.5 shows it costs overshoot.

### 2.2 The PID and the two-degree-of-freedom loop

The ideal derivative Kd s has unbounded gain at high frequency and amplifies sensor noise without limit, so it is always filtered: D(s) = Kd s/(Tf s + 1). The controller is

**C(s) = Kp + Ki/s + Kd s/(Tf s + 1) = [(Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki] / [s(Tf s + 1)].**   (3)

Its gain tends to Kp + Kd/Tf at high frequency, instead of infinity. Real implementations also feed the reference and the measurement differently. With a **setpoint weight** b on the proportional term and the derivative acting on the measurement only,

u = Kp(b r − y) + (Ki/s)(r − y) − D(s) y.

Collect: u = C_r(s) r − C(s) y with C_r = b Kp + Ki/s (plus D(s) if the derivative acts on the error). With y = P u, y(1 + PC) = P C_r r, so

**T(s) = y/r = P C_r / (1 + P C).**   (4)

The denominator 1 + PC, the **characteristic polynomial** den_P · den_C + num_P · num_C, contains C only. Therefore **b and the choice of derivative input do not move a single closed-loop pole**; they change only the zeros of T, i.e. how hard the reference "kicks" the loop. And T(0) = 1 whenever Ki > 0 and P has no zero at the origin, because C(s) → ∞ as s → 0 and C_r/C → 1: integral action gives zero steady-state error to a step, whatever b is (`P09.08`).

### 2.3 Ziegler–Nichols: the ultimate gain

With proportional control only, the closed-loop poles are the roots of den(s) + K num(s) = 0 (`P09.09`). As K grows, a pair crosses the imaginary axis at K = Ku, the **ultimate gain**, where the loop oscillates steadily at ω_u with period Pu = 2π/ω_u. On s = jω, K = −den(jω)/num(jω) must be real and positive, which needs Im[den(jω) num(−jω)] = 0: a real polynomial in ω whose positive roots are the candidate crossings (`ultimate_gain`). For (2): **Ku = 5.004, ω_u = 7.633 rad/s, Pu = 0.823 s**. The classic Ziegler–Nichols rules then set Kp = 0.6Ku, Ti = Pu/2, Td = Pu/8 (so Ki = Kp/Ti, Kd = Kp Td): **Kp = 3.002, Ki = 7.295 s⁻¹, Kd = 0.309 s**. The rules were tuned empirically for a quarter-amplitude decay, which is about 50 % overshoot: aggressive by flight-control standards. Note how they are obtained in a test: by raising the gain on the real system until it oscillates. On an aircraft that is a flight-test experiment nobody should run at low altitude (Problem 6).

### 2.4 Placing three poles exactly

The closed-loop polynomial of (2) with (3) is of degree 6 (four plant poles, the integrator and the filter pole of C), and the PID has three gains. It cannot place all six roots. It can place three: from (3), num_C = (Kp Tf + Kd)s² + (Kp + Ki Tf)s + Ki is **linear in (Kp, Ki, Kd)** for a fixed Tf, so the characteristic polynomial

χ(s) = den_P(s)(Tf s² + s) + num_P(s) num_C(s)

is an affine function of the gains at every s. Requiring χ(λ) = 0 for a real λ is one linear equation; for a complex pair λ, λ* the real and imaginary parts of χ(λ) = 0 are two. A complex pair plus one real pole give three equations in three unknowns (`pid_place_three_poles`). The remaining three roots fall where they fall, and **must be checked**.

Choosing a dominant pair at −4 ± 4j (ζ = 0.707, ω_n = 5.66 rad/s) and a real pole at −1.5 near the plant zero −1.2 (to keep the slow mode that the zero creates short-lived) gives **Kp = 2.417, Ki = 5.848 s⁻¹, Kd = 0.463 s** with Tf = 0.02 s; the other roots are −5.45 ± 1.79j and −53.3, all faster than the placed ones. For comparison the ZN gains give closed-loop poles −1.33 ± 6.13j (ζ = 0.21), −1.47, −2.33, −14.9 and −52.3. If the plant were the simple design model b/(s(s + a)) and Tf → 0, the polynomial would be cubic, s³ + (a + bKd)s² + bKp s + bKi, and matching s³ + c₁s² + c₂s + c₃ gives Kd = (c₁ − a)/b, Kp = c₂/b, Ki = c₃/b: the textbook result, which the test checks as a limiting case of the general solve.

### 2.5 Why the zeros matter: setpoint weighting

With b = 1 the placed design overshoots **41.1 %** on a 10° step although its dominant poles have ζ = 0.707 (4.3 % for a pure second-order system, `P09.04`). The extra overshoot comes from the zeros of T: the plant zero at −1.2 and the PID zeros near −Ki/Kp ≈ −2.4, both slower than the poles. With **b = 0** the proportional term acts on −y only, the reference enters through the integrator alone, the PID zeros leave T, and the overshoot falls to **10.6 %** with the same poles, the same disturbance rejection and the same stability margins; the rise time grows from 0.242 to 0.493 s. The remaining overshoot is the plant's own zero at −1.2, which no feedback can remove (only a reference prefilter can). Setpoint weighting is a free improvement whenever the reference is a step.

### 2.6 The discrete controller and anti-windup

Flight software runs the controller at a sample time dt (1 ms here, that is 1 kHz) and holds the command between samples. One sample of `pid_step`:

e = r − y;  x_d = −y (derivative on measurement);  D_k = (Tf D_{k−1} + Kd(x_d,k − x_d,k−1))/(Tf + dt);
v = Kp(b r − y) + I + D;  u = clip(v, ±δ_max);  I ← I + Ki dt e.

The derivative line is the backward-Euler discretisation of Tf Ḋ + D = Kd ẋ_d; the integrator is forward Euler. When u ≠ v the actuator is saturated and the integrator is integrating an error the actuator cannot act on: **windup**. Two standard fixes:

- **Clamping (conditional integration)**: skip the integrator update when the output is saturated *and* the error pushes it further into saturation, e(v − u) > 0.
- **Back-calculation**: I ← I + Ki dt e + (dt/Tt)(u − v). While saturated, the extra term drains the integrator toward the value that would make v equal to the limit, with time constant Tt; a common choice is Tt = √(Ti Td) (0.281 s here). Without saturation u = v and the term vanishes.

Both are exactly inactive when the limit is never reached: the test checks bit-for-bit equality.

### 2.7 Where the model stops being valid

- **Sampling.** The zero-order hold adds an effective delay of about dt/2; at 1 ms against a 6.4 rad/s crossover that costs 0.2° of phase. It is negligible here, not in general (`P09.17`).
- **Servo rate limit.** Real elevator servos are limited in rate (typically 40–100°/s) as well as position. A 20° step demanding 51° of elevator in the linear loop would hit the rate limit first; this model does not have one.
- **Airspeed and the phugoid.** Equation (1) holds airspeed constant. A sustained 20° pitch-up bleeds airspeed over tens of seconds and changes every coefficient.
- **Angle of attack.** Linear lift ends near stall. The 20° step with clamping anti-windup drives α to 14.0° in this model, close to the stall angle of many light aircraft: the saturated runs are at the edge of where (1) can be believed.

## 3. Build the computational model

Open `starter/main.py`. `pitch_plant`, `pitch_state_space`, `expm`, `discretize_zoh`, `tf_step` (exact step response of a TF through its companion realisation, `P09.13` explains why that works) and `step_metrics` are given.

### Step 1: transfer functions

`pitch_plant_tf` is (2): build Δ(s) and multiply polynomials with `np.polymul`. `pid_tf` is (3). `closed_loop_tf(num_p, den_p, gains, b, d_on_measurement)` returns (num_P · num_Cr, den_P · den_C + num_P · num_C), with num_Cr = (b Kp s + Ki)(Tf s + 1) when the derivative is on the measurement (sharing the denominator s(Tf s + 1) of C) and (b Kp Tf + Kd)s² + (b Kp + Ki Tf)s + Ki when it is on the error.

### Step 2: tuning

```python
def ultimate_gain(num, den):
    num_neg = num * (-1.0) ** np.arange(len(num) - 1, -1, -1)   # num(-s)
    q = np.polymul(den, num_neg)                                # Im q(jw) = 0 at a crossing
    ...  # keep odd powers with sign (-1)^((k-1)/2), np.roots, keep real w > 0 with real K > 0, take the smallest K
```

`ziegler_nichols` is three lines. `pid_place_three_poles` builds χ for the gains (0, 0, 0) and for each unit gain vector, evaluates them at the chosen poles, stacks real and imaginary parts and calls `np.linalg.solve`.

### Step 3: the discrete loop

`pid_step` is §2.6 line by line, returning (u, v, new_state) with state (I, D, x_prev); zero state means all signals were zero before t = 0, which is the transfer function's convention too. `simulate_pitch_loop` discretises the plant once with `discretize_zoh` (with a second input column for a constant disturbance moment), then loops: measure θ = Cx, call `pid_step`, apply u for one dt.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-11-pid-design-tuning`. The tests:

1. **Plant**: (2) equals C(sI − A)⁻¹B of the state-space model at four complex points and SciPy's realisation of it on a frequency grid; its poles are eig(A).
2. **PID and loop TFs** equal Kp + Ki/s + Kd s/(Tf s + 1) and P C_r/(1 + PC) pointwise for b ∈ {1, 0.5, 0} and both derivative inputs; DC gain 1.
3. **Independent check**: the sampled loop without a limit against `scipy.signal.step` of `closed_loop_tf`, for three configurations, within 0.6 % of the step; the error ratio between dt = 4, 2, 1 ms is 2.0 ± 15 % (first order); superposition; zero in, zero out.
4. **Tuning**: Ku puts a root at jω_u, |Ku P(jω_u)| = 1, the loop is stable at 0.98Ku and unstable at 1.02Ku; Ku = 6 and ω_u = √2 for 1/(s(s + 1)(s + 2)); the ZN ratios; three pole sets placed to 10⁻⁶; the b/(s(s + a)) closed form.
5. **Controller**: first-sample output Kp + Kd/(Tf + dt) with derivative on the error and Kp without; the integrator step Ki dt; the clamp freezes the integrator; the output never exceeds the limit.
6. **Independent check**: both anti-windup schemes are inactive without saturation and cut the 20° overshoot below 0.6 of the unprotected value while still converging; integral action removes a 2° trim disturbance that leaves a PD loop with a steady error.

### Compare

The comparison tool is `scipy.signal`: `step(lti(num, den))` integrates the closed-loop transfer function directly, with nothing in common with our sample-by-sample loop except the polynomials. Running `main()` prints, for the 10° step and dt = 1 ms, the largest difference between the two:

```
  ZN PID  b = 1: overshoot  64.7 %, rise 0.194 s, settle 2.82 s, peak elevator   32.3 deg; max |sim - TF| = 0.0336 deg
  ZN PID  b = 0: overshoot   8.6 %, rise 0.388 s, settle 2.10 s, peak elevator   11.3 deg; max |sim - TF| = 0.0176 deg
  placed  b = 1: overshoot  41.1 %, rise 0.242 s, settle 1.54 s, peak elevator   25.7 deg; max |sim - TF| = 0.0228 deg
  placed  b = 0: overshoot  10.6 %, rise 0.493 s, settle 2.29 s, peak elevator    8.5 deg; max |sim - TF| = 0.0097 deg
  convergence: dt = 4 ms, max |sim - TF| = 0.0921 deg
  convergence: dt = 2 ms, max |sim - TF| = 0.0457 deg
  convergence: dt = 1 ms, max |sim - TF| = 0.0228 deg
```

The residual is the sampling itself: it halves with dt, as a first-order method should (forward-Euler integrator, backward-Euler filter, half-sample hold delay). It is largest for the ZN design with b = 1, the one with the most violent initial transient.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Tuning-method comparison.** For the ZN and placed gains, with b = 1 and b = 0, record overshoot, rise time, 2 % settling time and peak elevator for a 10° step, and the phase margin of L = PC. Predict: b changes overshoot and peak elevator but not the poles or the margins; ZN is less damped.
2. **Sample-time convergence.** Repeat the TF comparison at dt = 4, 2, 1 ms and plot the error against dt on log axes. Predict a slope of 1.
3. **Anti-windup.** Command a 20° step with the elevator limited to ±15°, with no anti-windup, clamping, and back-calculation (Tt = 0.1, 0.3 s and √(Ti Td)). Record overshoot, settling time, the time spent saturated and the peak integrator value. Predict which scheme cuts the overshoot most.
4. **Disturbance.** Hold 5° pitch and apply a −2° equivalent-elevator trim moment at t = 3 s (a CG shift or flap deployment), with the full PID and with Ki = 0. Predict the PD steady error from §2.2 before running.

Raw data format: CSV tables `design, b, overshoot_pct, rise_s, settle_s, peak_elevator_deg, pm_deg`; `dt_s, max_err_deg`; `scheme, Tt_s, overshoot_pct, settle_s, saturated_s, peak_integrator_deg`; `controller, steady_error_deg`, each with a header stating the plant parameters and that every number is from simulation of a representative model.

What you should find: ZN with b = 1 overshoots 64.7 % (peak elevator 32.3°, more than twice the limit), placed with b = 1 41.1 % (25.7°), placed with b = 0 10.6 % (8.5°), ZN with b = 0 8.6 % (11.3°); the phase margins of L are 22.1° for ZN and 43.8° for the placed gains, at crossovers of 6.33 and 6.42 rad/s, unchanged by b. The gain margins are nearly equal (19.2 and 19.5 dB): ZN's weakness is phase, not gain. Under the limit: without anti-windup the overshoot is **67.9 %** (0.86 s saturated, integrator peaking at the equivalent of 51.1° of elevator), with clamping **13.2 %** (0.38 s saturated, integrator peak 13.2°), with back-calculation at Tt = √(Ti Td) = 0.281 s **23.4 %** (Problem 5 shows that a smaller Tt approaches clamping). The linear loop with no limit would overshoot 41.1 % and command 51.3° of elevator. The PD loop settles with a 0.827° error under the trim disturbance, exactly 2°/Kp; the PID removes it.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | cause of difference |
|---|---|---|---|---|
| short-period ω_n, ζ | √9 = 3.000, 3.7/6 = 0.617 | 3.000, 0.617 | eig(A) | — |
| Ku, Pu | Im[den num(−jω)] = 0 | 5.004, 0.823 s | root at jω_u to 10⁻⁶ | — |
| placed poles | −4 ± 4j, −1.5 | roots of χ | to 10⁻⁶ (test) | — |
| 10° step, placed, b = 1 | T(s) | sampled loop | `scipy.signal.step`: 0.0228° max diff | sample-and-hold, O(dt) |
| error ratio dt 4 → 2 → 1 ms | 2, 2 (first order) | 2.02, 2.00 | — | — |
| PD steady error under 2° trim | 2°/Kp = 0.827° | 0.827° | — | — |
| overshoot under saturation | (no linear prediction) | 67.9 / 13.2 / 23.4 % | linear loop 41.1 % | windup |

The linear results agree three ways (hand algebra, sampled loop, SciPy) to the accuracy the sample time allows, which validates the code. The saturation results have no linear prediction: the integrator dynamics change discontinuously at the limit. What validates them is the sanity of their mechanism: with clamping the integrator peaks at 13.2° instead of 51.1°, and the overshoot falls in proportion. They are still simulations of a model, not flight data; a flight test would see turbulence, a rate-limited servo, airspeed changes and sensor noise.

Improvements worth making: a servo rate limit and rate-limit-aware anti-windup (extension 1); a reference prefilter to cancel the effect of the plant zero; gain scheduling with airspeed, since M_δ scales with dynamic pressure. Improvement worth not making: driving Tt to zero in the back-calculation law. It approaches clamping and loses the smooth recovery that back-calculation offers with sensor noise.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (2), (4), the ultimate-gain condition and the three-pole equations; `pid_pitch.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending gains, b, Tf and an anti-windup scheme for the pitch hold, with the overshoot, peak elevator, phase margin and saturated overshoot of your choice, and one sentence on the first thing a flight test would have to check.

### Adversarial self-checks

- **What if the servo is twice as slow (ω_s = 10 rad/s)?** The servo lag at 6.4 rad/s grows by 14.9°. With the gains unchanged, the phase margin of the placed design falls from 43.8° to 30.0° and that of ZN from 22.1° to 8.0°. Re-solving the three-pole placement still puts roots at −4 ± 4j and −1.5 (Kp = 1.421, Ki = 4.327, Kd = 0.359), but an unplaced pair moves to −1.56 ± 3.20j and becomes the dominant, less damped mode. Always look at all six roots.
- **What if the aircraft is at higher speed?** M_δ and M_α grow roughly with dynamic pressure. A fixed gain set then has a higher crossover and a smaller phase margin. This is why autopilots schedule gains with airspeed.
- **Which assumption fails first in the real world?** The servo rate limit, for any step larger than a few degrees; then airspeed decay in a sustained pull-up; then turbulence and sensor noise through the derivative.
- **What measurement could falsify the model?** A small doublet in flight with the pitch rate recorded: if the measured short-period frequency or damping differs from 3.0 rad/s and 0.62, the plant (1) is wrong and so are the gains.
- **Does the result survive a change of sample time?** Yes in the linear regime (error ∝ dt, 0.02° at 1 ms); test it under saturation too, because a slower controller saturates for slightly longer.
- **Which decision changes if the uncertainty doubles?** An error in elevator power alone matters little: with M_δ 30 % high the crossovers rise to 7.84 and 7.35 rad/s but the phase margins barely move (42.1° placed, 22.0° ZN). An error in the servo bandwidth matters a lot (first check above). If the servo bandwidth is known only to a factor of two, ZN's margin can fall to 8°. Only the placed design keeps a usable margin, and only if it is re-checked on the slow-servo model.

## Check yourself

1. Write the transfer function of a PID with a filtered derivative and give its high-frequency gain.
2. Why does setpoint weighting (b = 0) reduce overshoot without changing any closed-loop pole?
3. What is the ultimate gain, and what condition on den(jω) and num(jω) locates it?
4. A PID has three gains. Why can it place only three of the six closed-loop poles of the pitch loop, and what must you do about the other three?
5. Explain integrator windup in two sentences, and say what clamping does differently from back-calculation.
6. A PD pitch loop with Kp = 2.5 holds attitude against a trim moment equivalent to 1° of elevator. What is the steady pitch error, and what removes it?
7. The sampled loop differs from the transfer-function prediction by 0.023° at dt = 1 ms. What do you expect at 0.5 ms, and why?
8. Why is the Ziegler–Nichols experiment (raise the gain until the loop oscillates) a poor idea on an aircraft in flight?
