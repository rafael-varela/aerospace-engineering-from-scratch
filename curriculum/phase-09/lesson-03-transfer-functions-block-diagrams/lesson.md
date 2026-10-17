# P09.03  Transfer Functions & Block Diagrams

> A block diagram is a set of linear equations drawn as a picture. Reduce it with polynomial algebra, then check the reduction the hard way: simulate every block of the unreduced diagram and see if the answers agree.

| | |
|---|---|
| **Learning outcomes** | 1. Write transfer functions for the blocks of a feedback loop (an elevator servo, short-period pitch dynamics, an integrator, a rate-gyro filter) and reduce series, parallel and feedback connections by hand, including the general G/(1 ± GH). 2. Implement and test a polynomial transfer-function algebra (series, parallel, feedback, poles, zeros, DC gain, pole-zero cancellation) and a controllable-canonical state-space realisation. 3. Check a reduction by an independent route: simulate the unreduced block diagram in the time domain, block by block, and compare with the step response of the reduced transfer function and with `scipy.signal`. 4. Read stability, steady state and closed-loop zeros off the reduced transfer function of a pitch-attitude autopilot, and judge when a reduced transfer function hides an unstable internal mode. |
| **Prerequisites** | `P09.01` (LTI systems), `P09.02` (impulse response, convolution), `P02.08` (Laplace transform, transfer functions, poles), `P02.05` (state-space form ẋ = Ax + Bu), `P03.05` (RK4), `P01.13` (linear systems). Skill check: you can multiply two polynomials, find the roots of a quadratic, and write a second-order ODE as two first-order ones. |
| **Engineering question** | A pitch-attitude autopilot for a light aircraft or small UAV closes an inner pitch-rate loop (a rate gyro with a 10 ms filter, gain K_q = 0.3 s) and an outer attitude loop (K_θ = 3) around a 50 ms elevator servo and a representative short-period model. What is the closed-loop transfer function θ/θ_c, is the loop stable, what are its steady-state gain and closed-loop zeros, how much elevator does a 5 deg step command need, and does a simulation of the unreduced diagram confirm the reduction? |
| **Assumptions and limits** | Every block is linear and time-invariant: small perturbations about one trimmed flight condition (fixed speed and altitude), a short-period approximation that ignores the slow phugoid mode, an elevator servo without rate or deflection limits, an ideal rate gyro behind a first-order filter. The numbers are a representative model, not any specific aircraft's data. Polynomial algebra never cancels a pole against a zero by itself, and numerical roots of high-degree expanded polynomials lose accuracy. The diagram simulator requires strictly proper blocks (no algebraic loops). |
| **Analytical work** | The transfer functions of the servo, the short period (from its two ODEs), the integrator and the gyro filter; series, parallel and feedback rules from the signal equations; the closed-loop characteristic equation 1 + GH = 0; limits K → 0 and K → ∞ (closed loop → 1/H); poles of H becoming zeros of the closed loop; the controllable canonical realisation and C(sI − A)⁻¹B + D = G; the two-loop pitch reduction by hand; pole-zero cancellation and internal stability. |
| **Computational deliverable** | `starter/main.py`: `tf`, `evaluate`, `series`, `parallel`, `feedback`, `poles`, `zeros`, `dc_gain`, `minreal`, `tf_to_ss`, `simulate_tf`, `reduce_pitch_loop` (`trim`, `poly_add`, `gain`, `simulate_diagram`, `pitch_model`, `pitch_blocks`, `pitch_connect` and `main` are given). Tests prove: `tf` normalises and rejects improper or empty transfer functions; series, parallel and both signs of feedback equal the block formulas evaluated at random complex s; K → 0 gives the open loop and K → ∞ gives 1/H; closed-loop poles satisfy 1 + GH = 0; DC gains including integrators; `minreal` cancels only coincident pairs; the realisation reproduces G pointwise, its eigenvalues are the poles and `scipy.signal.lti` recovers G; the RK4 step response matches `scipy.signal.step`; the pitch reduction equals the hand reduction, is fifth order with unit DC gain and has zeros at −1/T_θ2 and −1/τ_g; the unreduced diagram simulation matches the reduced TF and `scipy.signal.step` to 10⁻⁷; the loop is stable with the damper and unstable without it at K_θ = 4; a servo position loop reduces to s² + 20s + 400. |
| **Practical activity** | **computational experiment**: reduce the pitch loop by hand and by code, simulate the unreduced diagram, sweep the damper gain K_q and watch the closed-loop poles move, and find the attitude gain at which the loop without a damper goes unstable. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `block_diagram.png`, the hand reduction of the pitch loop, the table of closed-loop poles against K_q, the comparison of the reduced and unreduced simulations, and a short note on what the reduction hides (internal signals, cancellations, limits). |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the feedback formula and the two-loop reduction are derived by hand; all public tests pass; the reduced and unreduced responses agree and the agreement is explained; the note identifies at least one thing the reduced transfer function cannot show. |
| **References** | `@mit-ocw-16-06-principles-of-automatic-control` Lecture 2 "Block diagrams, the effect of feedback" (plant, actuator, sensor and sensor noise; open-loop vs. closed-loop control) and Lecture 4 "Block diagram manipulations, Mason' rule" [sic] (PDF "Block Diagram Manipulations": the single-loop gain rule; Mason's Rule) (required-free): the blocks of a control loop and the series, parallel and feedback rules. `@mit-ocw-16-30-feedback-control-systems` Topic 2 "Basic root locus: analysis and examples": Aircraft Longitudinal Control; First Short Period Autopilot (required-free): the pitch-attitude loop with an elevator servo that this lesson reduces. `@mit-ocw-6-003-signals-and-systems` Lecture 6 "Laplace Transform" (Solving Differential Equations with Laplace Transforms; Properties of Laplace Transforms) and Lecture 10 "Feedback and Control" (Structure of a Control Problem; Analysis of wallFinder System: Block Diagram; Analysis of wallFinder System: System Function) (required-free): G(s) as the transform of the impulse response, and L/(1 + L) for a whole loop. `@scipy-userguide` API reference "Signal processing (scipy.signal)": "Continuous-time linear systems" (`lti`, `lsim`, `step`) (required-free, oracle only). The polynomial-algebra reducer and its check by simulating the unreduced diagram are original instruction in this lesson; the pitch-loop numbers are a representative model written for it. |
| **Extensions** | ★ Implement Mason's gain formula for a signal-flow graph and check it on the pitch loop. ★ Add a feed-forward path from θ_c to the servo command and show that it adds a zero but leaves the closed-loop poles unchanged. ★ Extend `simulate_diagram` to blocks with a direct term (D ≠ 0) by ordering the evaluation, and detect algebraic loops. |
| **Time** | unmeasured |

## 1. The engineering problem

An autopilot is drawn as a block diagram: boxes for the airframe, the actuator, the sensors and the control law, joined by arrows (signals) that meet at summing junctions and split at pick-off points. Every box is a transfer function (`P02.08`), and the diagram as a whole is a set of linear equations. Reducing it to a single closed-loop transfer function answers the questions that matter first: is the loop stable, what is its steady-state gain, how fast and how well damped is it?

The vehicle for this lesson is the pitch-attitude loop of a light aircraft or small UAV. The elevator servo is a first-order lag, 1/(0.05s + 1). The short-period pitch dynamics take elevator deflection δ_e to pitch rate q with a representative model, q/δ_e = 12(s + 1.5)/(s² + 3.2s + 16), a natural frequency of 4 rad/s, a damping ratio of 0.4 and a zero at −1/T_θ2 = −1.5 s⁻¹. Pitch rate integrates to pitch angle, θ = q/s. A rate gyro measures q behind a 10 ms filter, and the autopilot commands δ_e = K_θ(θ_c − θ) − K_q q_measured: an inner "pitch damper" loop and an outer attitude loop. With K_q = 0.3 s and K_θ = 3, what does the aircraft do when the pilot or an outer navigation loop asks for a 5 deg pitch change?

Reductions are easy to get wrong. A sign at a summing junction, a block on the wrong side of a pick-off point, a cancelled factor that should not have been cancelled: each gives a transfer function that looks plausible. So the reduction is checked by a route that shares none of its algebra: every block of the *unreduced* diagram is realised in state-space form and the whole interconnection is integrated in time. Success criteria, written first: the coded reduction equals the hand reduction at random complex points to 10⁻⁹; the unreduced simulation matches the reduced transfer function's step response and `scipy.signal.step` to 10⁻⁷; and the steady-state gain is exactly 1.

## 2. Understand and derive

### 2.1 The blocks

**Servo.** A position-controlled actuator that follows its command with a lag: τ_s δ̇_e + δ_e = u, so δ_e/u = 1/(τ_s s + 1) with τ_s = 0.05 s.

**Short period.** For small perturbations from trim at constant speed, the angle of attack α and pitch rate q obey two linear ODEs (Phase 06 derives them): α̇ = Z_α α + q + Z_δ δ_e and q̇ = M_α α + M_q q + M_δ δ_e. Eliminating α and neglecting Z_δ gives a second-order transfer function with one zero:

q/δ_e = M_δ (s + 1/T_θ2)/(s² + 2ζ_sp ω_sp s + ω_sp²).   (1)

The lesson takes ω_sp = 4 rad/s, ζ_sp = 0.4, 1/T_θ2 = 1.5 s⁻¹ and M_δ = 12 s⁻², with positive δ_e defined as nose-up. These are representative of a light aircraft or small UAV at moderate speed, not data for a specific airframe. The loop they form is the one 16.30 Topic 2 sketches as its "First Short Period Autopilot".

**Integrator.** θ̇ = q, so θ/q = 1/s: one pole at the origin. **Gyro filter.** q_m/q = 1/(τ_g s + 1), τ_g = 0.01 s.

### 2.2 Series, parallel, feedback

Write each signal as the output of its block and solve (16.06 Lecture 4). **Series**: Y = G₂(G₁U) gives G₁G₂; with G_i = n_i/d_i, the numerator is n₁n₂ and the denominator d₁d₂. **Parallel**: Y = G₁U + G₂U gives (n₁d₂ + n₂d₁)/(d₁d₂). **Feedback**: with forward block G, feedback block H and error E = R ∓ HY (minus for negative feedback), Y = GE gives Y(1 ± GH) = GR:

**T = G/(1 + GH) = n_G d_H/(d_G d_H + n_G n_H)**   (negative feedback).   (2)

Three consequences follow from (2) alone. The closed-loop poles are the roots of **d_G d_H + n_G n_H**, equivalently of the **characteristic equation 1 + G(s)H(s) = 0**; they are not the poles of G or of H. The closed-loop zeros are the zeros of G **and the poles of H**: a sensor filter in the feedback path appears as a zero of the closed loop. And the limits: as the loop gain K → 0, T → KG (the loop is open); as K → ∞, T → **1/H**: a high-gain loop makes the output follow the inverse of the sensor, which is why feedback is only as good as its sensor.

The **DC gain** is G(0) = n(0)/d(0). A pole at s = 0 (an integrator) makes the open-loop DC gain infinite, and then unity feedback gives T(0) = G(0)/(1 + G(0)) = 1 exactly: an integrator in the forward path removes steady-state error to a step (`P09.08` will make this precise).

### 2.3 The pitch loop, reduced by hand

*Inner loop.* Forward path G_f = servo × short period = 12(s + 1.5)/[(0.05s + 1)(s² + 3.2s + 16)], normalised to 240(s + 1.5)/(s³ + 23.2s² + 80s + 320), poles −20 and −1.6 ± 3.666j, DC gain 1.125. Feedback H_q = K_q/(0.01s + 1). By (2):

q/u = G_f/(1 + K_q G_f/(0.01s + 1)) = 240(s + 1.5)(s + 100)/(s⁴ + 123.2s³ + 2400s² + 15 520s + 42 800).

The damper has moved the short-period poles from −1.6 ± 3.666j to −4.356 ± 3.498j (damping ratio 0.78), and the gyro filter pole −100 has become a zero.

*Outer loop.* θ/u = (q/u)/s, the open-loop transfer function from servo command to attitude, and with K_θ in the forward path and unity feedback:

**θ/θ_c = 720(s + 1.5)(s + 100)/(s⁵ + 123.2s⁴ + 2400s³ + 16 240s² + 115 880s + 108 000).**   (3)

Fifth order (servo, two short-period, gyro, integrator). DC gain 108 000/108 000 = 1. Zeros at −1.5 (the airframe's 1/T_θ2) and −100 (the gyro-filter pole). Poles −100.897, −17.241, −1.997 ± 7.358j and −1.068, all in the left half plane: **stable**. The slow real pole at −1.068 lies close to the zero at −1.5, a near-cancellation whose consequences `P09.04` explores.

### 2.4 From transfer function to simulation

A proper transfer function G = n/d with monic d = sⁿ + a_{n−1}sⁿ⁻¹ + ⋯ + a₀ has the **controllable canonical realisation**: the state x = (z, ż, …, z⁽ⁿ⁻¹⁾) of the auxiliary equation d(D)z = u, with

A = [[0, 1, 0, …], …, [−a₀, −a₁, …, −a_{n−1}]], B = (0, …, 0, 1)ᵀ, C = (c₀, …, c_{n−1}), D = b_n,

where b_n is the sⁿ coefficient of n and c_k are the coefficients of the strictly proper remainder n − b_n d. One checks C(sI − A)⁻¹B + D = G (Problem 2) and that the eigenvalues of A are the roots of d. `simulate_tf` integrates ẋ = Ax + Bu with fixed-step RK4 (`P03.05`).

`simulate_diagram` does the same for a whole diagram: each block gets its own realisation, the joint state is the stack of block states, and at every RK4 stage a wiring function computes each block's input from the current block outputs and the reference. Because every block is strictly proper, its output depends only on its state, so the wiring is explicit (no algebraic loop). The gains and summing junctions live in the wiring function; for the pitch loop it is the single line u = K_θ(θ_c − θ) − K_q q_m.

### 2.5 What a reduced transfer function hides

Polynomial algebra does not cancel. series((s + 2)/(s + 1), 3/(s + 2)) is 3(s + 2)/((s + 1)(s + 2)), and only `minreal` (which removes pole-zero pairs that coincide to a tolerance) gives 3/(s + 1). When the cancelled pole is stable this is bookkeeping. When it is **unstable**, it is a hidden failure: a controller zero placed on an unstable plant pole gives a reduced transfer function that looks stable, but the internal signal of the cancelled mode grows from any disturbance (Problem 5). A transfer function describes the response from one input to one output; internal stability needs every closed-loop pole, which the unreduced diagram keeps.

Two numerical traps complete the picture. Roots of an expanded polynomial are sensitive to its coefficients: np.roots of (s + 1)¹⁰ expanded is wrong by 0.05 and of (s + 1)¹⁶ by 0.23, from rounding alone. Keep loops in factored form when you can and reduce numerically only when you must. And a high-order reduced transfer function gives no physical insight on its own; the block diagram does.

### 2.6 Where the model stops being valid

The short-period model holds for small perturbations about one flight condition; its coefficients change with speed, altitude and centre of gravity, and the phugoid (a slow speed-altitude exchange) is absent. The servo saturates: the elevator needed for a 5 deg attitude step is 11.5 deg at its peak (2.295 rad per rad of command), and an elevator limited to about 25 deg allows linear steps up to roughly 11 deg only. The servo also has a rate limit. The gyro has noise and bias (`P09.01`), and the flight computer samples everything (`P09.17`).

## 3. Build the computational model

Open `starter/main.py`. `trim`, `poly_add`, `gain`, `simulate_diagram`, the pitch-model functions and `main` are given.

### Step 1: the algebra

```python
def tf(num, den):
    num, den = trim(num), trim(den)          # drop leading zeros
    ...                                      # raise on a zero den or deg num > deg den
    return num / den[0], den / den[0]        # monic denominator

def feedback(G, H=None, sign=-1):
    nG, dG = G
    nH, dH = (np.array([1.0]), np.array([1.0])) if H is None else H
    return tf(np.polymul(nG, dH), poly_add(np.polymul(dG, dH), -sign * np.polymul(nG, nH)))
```

`series` and `parallel` follow 2.2 with `np.polymul` and `poly_add`. `poles` and `zeros` are `np.roots` sorted by real part; `dc_gain` returns `inf` for a pole at the origin; `minreal` pairs zeros with poles within a tolerance, rebuilds both polynomials with `np.poly`, and fixes the gain by matching G at a test point away from every root.

### Step 2: realisation and simulation

`tf_to_ss` builds (A, B, C, D) of 2.4. `simulate_tf(G, u, t)` integrates ẋ = Ax + Bu(t) with RK4, `substeps` steps per grid interval, and returns y = Cx + Du.

### Step 3: the pitch loop

`reduce_pitch_loop(p)` returns the stages of 2.3: `forward = series(servo, short_period)`, `inner = feedback(forward, series(gain(K_q), gyro))`, `theta_u = series(inner, integrator)`, `closed = feedback(series(gain(K_θ), theta_u))`.

Running `main()` prints each stage (numerator, denominator, poles, zeros, DC gain), the comparison of the two simulations, and the pole table of part 4, and writes `block_diagram.png`.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-03-transfer-functions-block-diagrams`. The tests:

1. **Normalisation**: monic denominator, leading zeros removed, improper and zero-denominator inputs rejected.
2. **Block formulas** at five random complex points for three pairs of random transfer functions: series = G₁G₂, parallel = G₁ + G₂, feedback = G₁/(1 ± G₁G₂) and G₁/(1 + G₁).
3. **Limits**: K = 10⁻⁹ gives KG, K = 10⁹ gives 1/H.
4. **Characteristic equation**: |1 + GH| < 10⁻⁷ at every closed-loop pole; series poles are the union.
5. **DC gain**: plain, integrator (∞), parallel (sum), and unit DC gain of a unity-feedback type-1 loop for three gains.
6. **Cancellation**: `minreal` removes (s + 2) from (s + 2)/((s + 1)(s + 2)), keeps the values, and leaves (s + 2.1)/((s + 1)(s + 2)) alone.
7. **Realisation**: C(sI − A)⁻¹B + D = G at five points for four transfer functions, eig A = poles, `scipy.signal.lti(A, B, C, D)` recovers G.
8. **Simulation**: the RK4 step response matches `scipy.signal.step` to 10⁻⁸ for three systems, one with D ≠ 0.
9. **Pitch loop**: reduction = hand formula (built independently in the test), order 5, DC gain 1, zeros −100 and −1.5; the unreduced diagram matches `scipy.signal.step` of the reduced TF to 10⁻⁷ and `simulate_tf` to 10⁻⁹; the measured rate decays and θ → θ_c.
10. **Stability**: all poles in the left half plane with the damper; K_q = 0, K_θ = 4 has a right-half-plane pair and its simulated error grows.
11. **Servo loop**: K/(s(τs + 1)) with unity feedback reduces to s² + 20s + 400, and the diagram simulation agrees; a block with D ≠ 0 is refused.
12. **main** writes the figure.

### Compare

`scipy.signal` has the same objects (`lti`, `TransferFunction`, `StateSpace`) but no block-diagram algebra; its `tf2ss` returns a realisation in a different (but similar) companion form. Our realisation and SciPy's give the same transfer function to 10⁻⁸ (test 7), and the RK4 simulation of the unreduced pitch diagram matches `scipy.signal.step` (which uses the exact matrix exponential) to 1.1 × 10⁻¹⁰ on a 5 ms grid with four sub-steps. The unreduced and reduced RK4 simulations agree with each other to 1.7 × 10⁻¹⁵: not a coincidence but a theorem, since RK4 applied to two similar realisations of the same system gives identical results in exact arithmetic. So the comparison between the two RK4 runs checks the algebra, and the comparison with `scipy.signal.step` checks the integration.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Reduce and simulate.** Run `main()`. Before looking, predict from (3): the DC gain, the order, the closed-loop zeros, whether the unit-step response overshoots (the dominant complex pair has ζ = 0.26, but there is a zero at −1.5 and a pole at −1.07 near it). Record the peak, the peak time, the 2 % settling time and the peak elevator per radian of command.
2. **Damper sweep.** For K_θ = 3 and K_q = 0, 0.1, 0.3, 0.5, tabulate the closed-loop poles and the largest real part. Predict how the complex pair moves as K_q rises.
3. **No damper.** Set K_q = 0 and increase K_θ to 4. Predict whether the loop is stable, and confirm with both the poles and the unreduced simulation over 60 s.
4. **Cancellation.** Reduce series((s + 2)/(s + 1), 3/(s + 2)) and (s + 2)/((s + 1)(s + 2)) with and without `minreal`, and compute the roots of the expanded (s + 1)ⁿ for n = 4, 10, 16.

Raw data format: CSV tables `stage, num, den, poles, zeros, dc_gain` (1); `K_q, pole_1 … pole_5, max_real` (2); `K_theta, K_q, max_real, sim_growth` (3); `n, max_root_error` (4).

What you should find. (1) The unreduced and reduced responses agree to 1.7 × 10⁻¹⁵ rad. The response peaks at 1.078 at t = 0.51 s (7.8 % overshoot), settles to within 2 % in 2.75 s and reaches 0.9995 at 6 s; the elevator peaks at 2.295 rad per rad of command, 11.5 deg for a 5 deg step. The complex pair alone would predict 43 % overshoot; the zero and the slow pole cut it to 7.8 % (`P09.04`). (2) With K_q = 0 the largest real part is −0.239 (the pair at −0.24 ± 6.72j, barely damped); K_q = 0.1 moves it to −0.763, K_q = 0.3 to −2.00 ± 7.36j, after which the slow real pole at −1.07 is the slowest and more damping no longer helps (−1.042 at K_q = 0.5). (3) K_q = 0, K_θ = 4 puts a pair at +0.027 ± 7.42j: unstable; the simulated attitude error over the last 10 s is more than twice that between 20 and 30 s. The pitch damper is what makes a stiff attitude loop possible. (4) Without `minreal` the denominators keep the cancelled factor; the roots of (s + 1)ⁿ are wrong by 2.2 × 10⁻⁴, 0.051 and 0.23.

## 5. Validate and improve

| quantity | hand (§2.3) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| inner-loop denominator | s⁴ + 123.2s³ + 2400s² + 15 520s + 42 800 | same | — | 0 | exact |
| closed-loop θ/θ_c | (3) | same at 5 complex points | hand formula in the test | < 10⁻⁹ relative | rounding |
| DC gain | 1 | 1 | θ(6 s) = 0.9995 in simulation | slow pole still decaying | e^{−1.07·6} |
| closed-loop zeros | −1.5, −100 | −1.5, −100 | poles of H, zeros of G | 0 | — |
| step response, reduced vs unreduced (RK4) | equal | 1.7 × 10⁻¹⁵ rad | — | — | similar realisations |
| step response vs `scipy.signal.step` | equal | 1.1 × 10⁻¹⁰ rad | exact matrix exponential | — | RK4 error at h = 1.25 ms |
| stability, K_q = 0.3 | all Re < 0 | max Re = −1.068 | simulation decays | — | — |
| stability, K_q = 0, K_θ = 4 | unstable | max Re = +0.027 | simulation grows | — | — |

The reduction is right, by two routes that share no algebra. What the comparison cannot validate is the model: the short-period coefficients are representative, not measured; the servo has limits; and the attitude loop's apparent stiffness depends on a gyro that has noise and bias. The improvement most worth making is to keep the diagram, not just the reduced transfer function, in every later analysis: the internal signals (elevator, measured rate) are what saturate and what a flight test records.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand reduction of 2.3 and the derivation of (2); `block_diagram.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the closed-loop transfer function, its poles and zeros, the step-response figures, the peak elevator for a 5 deg command, and two things the reduced transfer function cannot tell you.

### Adversarial self-checks

- **What if the loop gain goes to zero or to infinity?** K_θ → 0 leaves θ/θ_c → 0 (no command authority) and the closed-loop poles tend to the open-loop ones; raising K_θ with K_q = 0.3 pushes the complex pair toward the imaginary axis, and at K_θ = 9.62 it crosses at ±11.03j: the loop has a gain margin of about 3.2 over the design value (`P09.09` draws the locus, `P09.10` names the margin).
- **Which assumption fails first?** Elevator saturation: with about 25 deg of travel, steps larger than about 11 deg of attitude command exceed it; then the servo rate limit. The linear reduction says nothing about either.
- **What if the gyro filter is ten times slower (τ_g = 0.1 s)?** Its pole at −10 is no longer far from the loop's poles: the closed-loop zero moves to −10, the complex pair moves to −0.94 ± 8.02j (damping ratio 0.12 instead of 0.26) and the overshoot doubles from 7.8 % to 15.7 %. Re-run the sweep before believing the K_q = 0.3 design.
- **What measurement could falsify the model?** A small attitude step in flight (or on a simulator with a validated model): the overshoot and the period of the residual oscillation should match the pair −2.0 ± 7.36j (about 0.85 s period); a slow tail with the −1.07 time constant should be visible.
- **Does the result survive a change of tool?** Hand algebra, the polynomial reducer, the unreduced simulation and `scipy.signal` agree to 10⁻⁷ or better. They would all agree on a wrong diagram, which is why the wiring function is written from the drawing, not from the reduced formula.
- **Which decision changes if the short-period damping is half as large?** ζ_sp = 0.2 moves the open-loop pair toward the axis; the inner damper must do more, and the K_q sweep should be repeated to keep the outer loop's pair damped.

## Check yourself

1. Write the closed-loop transfer function of a negative-feedback loop with forward block G and feedback block H, and the equation whose roots are the closed-loop poles.
2. A feedback sensor has transfer function 1/(0.01s + 1). Where does its pole appear in the closed-loop transfer function?
3. A loop has open-loop transfer function 5/(s(s + 2)) with unity feedback. What are the closed-loop poles and DC gain?
4. What does a negative-feedback loop tend to as its gain becomes very large, and what does this say about the sensor?
5. Why does polynomial algebra give (s + 2)/((s + 1)(s + 2)) instead of 1/(s + 1), and when is cancelling the factor dangerous?
6. What are A, B, C, D of the controllable canonical realisation of (2s + 1)/(s² + 3s + 5)?
7. Two RK4 simulations of the same loop, one reduced and one unreduced, agree to 10⁻¹⁵. What does that check, and what does it not check?
