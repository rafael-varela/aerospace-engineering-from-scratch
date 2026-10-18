# P09.13  State-Space Models

> A transfer function says what comes out. A state-space model also says what is going on inside. The same input-output behaviour has infinitely many inside stories, and they must all tell the same transfer function.

| | |
|---|---|
| **Learning outcomes** | 1. Write a state-space model ẋ = Ax + Bu, y = Cx + Du directly from physical ODEs (the short-period pitch dynamics extended with pitch angle and altitude) and from a transfer function (controllable and observable canonical forms). 2. Derive and implement the transfer function of any (A, B, C, D) by the Faddeev–LeVerrier algorithm, and show it is unchanged by any similarity transform x = Tz, including the real modal form. 3. Implement and test simulation by the zero-order-hold matrix exponential, and verify a model three independent ways: C(sI − A)⁻¹B + D evaluated pointwise against a hand-derived polynomial TF, simulation against `scipy.signal.lti`, and different realisations against each other. 4. Judge which coordinates to use: physical for meaning and limits, modal for insight and numerical robustness, companion only for small orders. |
| **Prerequisites** | `P09.03` (transfer functions), `P09.04` (poles and time response), `P02.04` (systems of ODEs), `P02.05` (the matrix exponential and zero-order-hold discretisation), `P01.16` (eigenvalues, diagonalisation), `P01.13` (matrix inverse and products). Skill check: you can turn ÿ + 3ẏ + 2y = u into two first-order equations. |
| **Engineering question** | A small UAV's longitudinal model must feed an altitude-hold design (`P09.15`, `P09.20`). Written from the physics, the states are angle of attack, pitch rate, pitch angle and altitude. Is that 4-state model the same system as the transfer function h/δ = −V M_δ Z_α/(s²Δ(s)) derived by hand, and does it simulate the same as SciPy? What does a 1° elevator step do in 6 s? Answer: C(sI − A)⁻¹B + D equals 432/(s²(s² + 3.7s + 9)) to 6 × 10⁻¹⁶ relative; the physical, companion and observable realisations agree with each other and with `scipy.signal.step` to 10⁻¹³ m; after 6 s the aircraft has climbed 13.06 m with pitch 10.28°, angle of attack 1.33° and a pitch rate that, in this constant-airspeed model, never stops. |
| **Assumptions and limits** | Linear, time-invariant, single input and single output (the methods extend to MIMO); for the aircraft, the short period with constant airspeed (no phugoid), small angles, Z_δ = 0, still air, altitude rate V(θ − α). The representative numbers are those of `P09.11`, V = 30 m/s, not a specific aircraft. The model predicts a steady pitch rate under a held elevator, which a real aircraft cannot sustain because it loses speed. Companion forms are numerically fragile beyond order 6–8; the modal form needs a diagonalisable A. |
| **Analytical work** | The pitch-altitude ODEs as (A, B, C, D); h/δ by elimination; the controllable canonical form from the ODE with an internal variable; the observable form as its dual; TF invariance under x = Tz; the Faddeev–LeVerrier recursion from (sI − A) adj(sI − A) = det(sI − A) I; the real 2 × 2 modal block of a complex pair; the DC gain −CA⁻¹B + D; the steady short-period state under a held elevator. |
| **Computational deliverable** | `starter/main.py`: `tf_eval`, `ss_eval`, `tf_to_controllable`, `tf_to_observable`, `faddeev_leverrier`, `ss_to_tf`, `similarity`, `modal_form`, `discretize_zoh`, `simulate`, `dc_gain`, `pitch_altitude_ss`, `pitch_altitude_tf_hand` (`pitch_altitude_params`, `expm`, `round_sig`, `main` are given). Tests prove: both canonical forms reproduce three TFs pointwise (including a biproper one) and the companion form is exactly `scipy.signal.tf2ss`'s; Faddeev–LeVerrier equals `np.poly`, satisfies Cayley–Hamilton and gives adj(sI − A); `ss_to_tf` matches `scipy.signal`'s and round-trips; similarity and modal transforms preserve eigenvalues and the TF, the modal form is real block diagonal; **the physical pitch-altitude model equals the hand TF pointwise** (the lesson's independent check) with two poles at the origin; **its step response equals `scipy.signal.step`** and a held random input equals `scipy.signal.lsim(interp=False)` to 10⁻⁹; ZOH matrices equal `scipy.linalg.expm` and a quadrature; superposition, time invariance, the free response and the DC gain. |
| **Practical activity** | **computational experiment**: compare realisations of the same system (physical, companion, observable, modal) in simulation, and a coefficient-sensitivity study of companion versus modal forms as the order grows. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `state_space.png`, the hand derivations of h/δ and of the canonical forms, the realisation-comparison and sensitivity tables, and a short note on which coordinates you would put in flight software and why. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: h/δ and the companion form are derived by hand; all public tests pass; the note explains the coefficient-sensitivity result with numbers and states the constant-airspeed limit of the aircraft model. |
| **References** | `@mit-ocw-16-30-feedback-control-systems` Topic 5 "Introduction to state-space models" (https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec05/): Basic Definitions; Creating State-Space Models; Equilibrium Points; Linearization; Stability of LTI Systems (required-free): states from the physics. Topic 6 "Developing state-space models based on transfer functions": TF's to State-Space Models; Modal Form; SS to TF; State-Space Transformations (required-free): companion and modal forms, G(s) = C(sI − A)⁻¹B + D and its invariance under x = Tz. Topic 7 "State-space models: basic properties": Time Response; SS: Forced Solution; Eigenvalues and Eigenvectors; Jordan Form; Diagonalization with Complex Roots (required-free): the forced solution and the real 2 × 2 block of a complex pair. Topic 8 "System zeros and transfer function matrices" (supplementary-free). `@mit-ocw-6-011-signals-systems-and-inference` Lecture 4 "State-Space Models"; Lecture 5 "State-Space Models, Equilibrium, Linearization"; Lecture 6 "Modal Solution of Undriven CT LTI State-Space Models"; Lecture 8 "Matrix Exponential, ZIR+ZSR, Transfer Function, Hidden Modes, Reaching Target States" (https://ocw.mit.edu/courses/6-011-signals-systems-and-inference-spring-2018/resources/mit6_011s18lec8/) (supplementary-free; sparse slides): a second view of modal solutions and the matrix exponential. `@mit-ocw-18-06-linear-algebra` Lecture 23 "Differential equations and exp(At)", Lecture 28 "Similar matrices and Jordan form" (required-free): e^{At}, similarity and the defective case; Lecture 21 "Eigenvalues and eigenvectors", Lecture 22 "Diagonalization and powers of A" (supplementary-free review). `@scipy-userguide` API reference "Signal processing (scipy.signal)": "LTI representations" (`tf2ss`, `ss2tf`), "Continuous-time linear systems" (`lti`, `lsim`, `step`); "Linear algebra (scipy.linalg)": Matrix Functions (`expm`) (supplementary-free): the test oracles only. |
| **Extensions** | ★ Add the phugoid (airspeed as a state, the X_u and Z_u derivatives) and show that the steady pitch rate disappears. ★ Extend `ss_to_tf` to MIMO (a matrix of numerators) and get the transfer matrix from elevator and throttle to pitch and altitude. ★ Implement a balanced realisation (`P09.14` Gramians) and compare its coefficient sensitivity with the modal form. |
| **Time** | unmeasured |

## 1. The engineering problem

Up to `P09.12` a system was a transfer function: a ratio of polynomials that maps input to output. That is enough to design a single loop, but it hides three things a flight-control engineer needs. The first is the inside of the system. A transfer function from elevator to altitude says nothing about the angle of attack on the way, and the angle of attack is what stalls the wing. The second is the initial conditions: a transfer function assumes the system starts at rest, but an aircraft in a gust does not. The third is many inputs and many outputs: an aircraft has elevator, throttle and flaps, and pitch, speed and altitude, and a pile of scalar transfer functions is a clumsy way to hold that.

The **state-space model** ẋ = Ax + Bu, y = Cx + Du keeps all three. Its state x is a list of numbers that, together with the future input, determines the future completely. Here that is angle of attack, pitch rate, pitch angle and altitude. The price is that the description is not unique: any invertible change of coordinates x = Tz gives another valid (A, B, C, D) with the same input-output behaviour. This lesson builds state-space models both ways, from the physics and from a transfer function, converts back by an algorithm that needs only matrix products and traces, and checks that every route tells the same story.

Success criteria, written first: the physical 4-state pitch-altitude model and the hand-derived TF agree pointwise in the complex plane to 10⁻¹² relative; the physical, companion and observable realisations simulate a 1° elevator step identically to 10⁻¹² m and agree with `scipy.signal.step`; similarity and modal transforms leave the TF unchanged to 10⁻⁹.

## 2. Understand and derive

### 2.1 State-space from the physics

The short-period equations of `P09.11`, without the servo, for angle of attack α, pitch rate q and pitch angle θ, plus the altitude h: climb rate is airspeed times the flight-path angle γ = θ − α (small angles, still air), so ḣ = V(θ − α). With the elevator δ as the input and h as the output:

```
d/dt [α]   [Z_α   1   0   0] [α]   [ 0 ]
     [q] = [M_α  M_q  0   0] [q] + [M_δ] δ,     y = h = [0 0 0 1] x,  D = 0.   (1)
     [θ]   [ 0    1   0   0] [θ]   [ 0 ]
     [h]   [−V    0   V   0] [h]   [ 0 ]
```

Every entry has a meaning and a unit. Rows 1 and 2 are the short period, row 3 is kinematics and row 4 is geometry. The eigenvalues of A are the short-period pair, −1.85 ± 2.362j for the representative numbers (Z_α = −1.2 s⁻¹, M_α = −6 s⁻², M_q = −2.5 s⁻¹, M_δ = 12 s⁻², V = 30 m/s), and two zeros: θ and h are integrals that nothing pulls back.

**The TF by hand.** In Laplace form, from row 1, q = (s − Z_α)α. Row 2 then gives α Δ(s) = M_δ δ with Δ = s² − (Z_α + M_q)s + Z_αM_q − M_α, as in `P09.11`. So α/δ = M_δ/Δ and θ/δ = q/(sδ) = M_δ(s − Z_α)/(sΔ). The flight path is γ/δ = (θ − α)/δ = M_δ[(s − Z_α) − s]/(sΔ) = −M_δZ_α/(sΔ), and h = Vγ/s:

**h/δ = −V M_δ Z_α / (s² Δ(s)) = 432 / (s² (s² + 3.7s + 9)).**   (2)

The numerator is positive (Z_α < 0) and has no zero: altitude responds to elevator only through the lift change that α creates.

**The steady state under a held elevator.** Set α̇ = q̇ = 0 in rows 1 and 2: q = −Z_αα and (M_α − M_qZ_α)α = −M_δδ, so α_ss = M_δδ/(Z_αM_q − M_α) = 12/9 δ = **1.333° per degree** and q_ss = −Z_αα_ss = **1.600°/s per degree**. In this model the aircraft keeps pitching up at 1.6°/s for as long as the elevator is held: it loops. That is not a bug in the algebra. It is the constant-airspeed assumption, and it tells you exactly where the model stops (§2.7).

### 2.2 State-space from a transfer function: the canonical forms

For G(s) = (b₀sⁿ + b₁sⁿ⁻¹ + ⋯ + b_n)/(sⁿ + a₁sⁿ⁻¹ + ⋯ + a_n), introduce an internal variable w with w⁽ⁿ⁾ + a₁w⁽ⁿ⁻¹⁾ + ⋯ + a_nw = u. Its derivatives (w⁽ⁿ⁻¹⁾, …, w) are a state. The first state's derivative is w⁽ⁿ⁾ = −a₁w⁽ⁿ⁻¹⁾ − ⋯ − a_nw + u, and each other state is the derivative of the next. The output y = b₀w⁽ⁿ⁾ + ⋯ + b_nw, after substituting w⁽ⁿ⁾, is

**A = [[−a₁, −a₂, …, −a_n], [1, 0, …, 0], …, [0, …, 1, 0]],  B = e₁,  C = [b₁ − b₀a₁, …, b_n − b₀a_n],  D = b₀**:   (3)

the **controllable canonical (companion) form** (`tf_to_controllable`; it is also what `scipy.signal.tf2ss` returns; 16.30 Topic 6, "TF's to State-Space Models", builds it the same way). Because G is a scalar, G = Gᵀ = Bᵀ(sI − Aᵀ)⁻¹Cᵀ + D, so (Aᵀ, Cᵀ, Bᵀ, D) is another realisation: the **observable canonical form**. The names will make sense in `P09.14`.

### 2.3 Back from state space: Faddeev–LeVerrier

G(s) = C(sI − A)⁻¹B + D = C adj(sI − A) B/det(sI − A) + D. Write det(sI − A) = sⁿ + c₁sⁿ⁻¹ + ⋯ + c_n and adj(sI − A) = N₁sⁿ⁻¹ + N₂sⁿ⁻² + ⋯ + N_n, matrices of numbers. Expand (sI − A) adj(sI − A) = det(sI − A) I and compare powers of s: N₁ = I, N_k = AN_{k−1} + c_{k−1}I, and the constant term gives AN_n + c_nI = 0 (Cayley–Hamilton). Taking the trace of each step, using tr(adj) = d(det)/ds, gives c_k = −tr(AN_k)/k. So

**N₁ = I, c₁ = −tr A;  N_k = AN_{k−1} + c_{k−1}I, c_k = −tr(AN_k)/k;  num = Σ C N_k B s^{n−k} + D · den.**   (4)

The algorithm uses only products and traces, and Cayley–Hamilton is a built-in residual check. Applied to (1) it returns den = s⁴ + 3.7s³ + 9s² + (4 × 10⁻¹⁶)s − 8 × 10⁻¹⁶ and num = 432: equation (2), with rounding where the two zeros at the origin should be. Its weakness is accuracy at high order (Problem 6).

### 2.4 Changes of coordinates

Put x = Tz with T invertible. Then ż = T⁻¹ATz + T⁻¹Bu and y = CTz + Du, and

C T (sI − T⁻¹AT)⁻¹ T⁻¹ B + D = C T [T⁻¹(sI − A)T]⁻¹ T⁻¹B + D = C(sI − A)⁻¹B + D.   (5)

The TF and the eigenvalues are invariant (16.30 Topic 6, "State-Space Transformations"). Everything else (the matrix entries, the meaning of the states, the numerical conditioning) changes. The **modal form** picks T from the eigenvectors. A real eigenvalue contributes its eigenvector. A complex pair σ ± jω with eigenvector v = v_r + jv_i contributes the two real columns (v_r, v_i). Since Av = (σ + jω)v gives Av_r = σv_r − ωv_i and Av_i = ωv_r + σv_i, the block of T⁻¹AT is **[[σ, ω], [−ω, σ]]**. For the short period, `modal_form` returns [[−1.85, 2.3617], [−2.3617, −1.85]]: the damping and the damped frequency, read straight off the matrix. The full model (1) cannot be put in modal form: its double eigenvalue at 0 has only one eigenvector (θ and h form a chain of integrators, a defective matrix as in `P02.05`).

### 2.5 Simulation by the matrix exponential

With u held constant over each step of length dt (zero-order hold), the exact update is x_{k+1} = A_d x_k + B_d u_k with A_d = e^{A dt} and B_d = ∫₀^dt e^{As}ds B, both from one exponential of the augmented matrix [[A, B], [0, 0]] (`P02.05`). `simulate` loops it and forms y_k = Cx_k + Du_k. For a step input it is exact at the sample instants, whatever dt is.

### 2.6 What the 1° step does

A 1° trailing-edge-up elevator step, simulated for 6 s with dt = 10 ms: **h(6 s) = 13.062 m, θ(6 s) = 10.276°, α(6 s) = 1.3333°**. The short period has settled (α at its steady value 4/3°), the pitch rate is 1.6°/s, and the climb is accelerating as γ = θ − α grows.

### 2.7 Where the model stops being valid

- **Constant airspeed.** The model holds V constant. A real aircraft climbing at a growing flight-path angle slows down within seconds to tens of seconds (the phugoid), and lift, M_δ and the pitch rate all fall. The predicted 13 m climb in 6 s is plausible. The predicted loop is not.
- **Small angles.** ḣ = V(θ − α) assumes sin γ ≈ γ. At γ = 9° the error is 0.4 %; at 30° it is 4.5 %.
- **Still air.** Wind enters as an extra input on α and h; the state-space form takes it as another column of B.
- **Numerical form.** The companion form of a high-order system is so sensitive to its coefficients that rounding them to five figures can move poles by more than one unit (§4).

## 3. Build the computational model

Open `starter/main.py`. `pitch_altitude_params`, `expm` and `round_sig` are given.

### Step 1: evaluate and realise

`tf_eval` is `np.polyval(num, s)/np.polyval(den, s)`. `ss_eval` solves (sI − A)X = B for each s instead of inverting. `tf_to_controllable` normalises den to monic, pads num to n + 1 coefficients, and builds (3). `tf_to_observable` transposes it.

### Step 2: Faddeev–LeVerrier

```python
c = np.zeros(n + 1); c[0] = 1.0
N = np.eye(n); Ns = [N]; c[1] = -np.trace(A @ N)
for k in range(2, n + 1):
    N = A @ N + c[k - 1] * np.eye(n); Ns.append(N)
    c[k] = -np.trace(A @ N) / k
```

`ss_to_tf` then forms num[k] = C N_k B and adds D times den.

### Step 3: coordinates and simulation

`similarity` is (5). `modal_form` sorts the eigenvalues and pairs each complex one with its conjugate before building T. `discretize_zoh` builds the augmented matrix, and `simulate` loops the update. `dc_gain` is D − CA⁻¹B with a solve. `pitch_altitude_ss` is (1) and `pitch_altitude_tf_hand` is (2).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-13-state-space-models`. The tests:

1. **Canonical forms** reproduce three TFs (second order, a non-monic third order, a biproper one with D = 3) at five complex points to 10⁻¹²; the companion form equals `scipy.signal.lti(num, den).to_ss()` entry by entry; an improper TF raises.
2. **Faddeev–LeVerrier** equals `np.poly(A)` for random 3 × 3 and 5 × 5 matrices, satisfies AN_n + c_nI = 0, and its adjugate equals det · inverse at a test point; `ss_to_tf` equals SciPy's conversion and round-trips.
3. **Coordinates**: a random similarity keeps the eigenvalues and the TF; the modal form of a random 5-state system is real, block diagonal, with σ on the diagonal and ±ω off it, and the same TF; the short-period block is [[σ, ω], [−ω, σ]] with σ = (Z_α + M_q)/2.
4. **Independent check, part 1**: the physical model (1) equals the hand TF (2) at five points to 10⁻¹², both directly and through Faddeev–LeVerrier; A has exactly two eigenvalues at the origin; the TF doubles when V doubles.
5. **Independent check, part 2**: the step response of (1) equals `scipy.signal.step(lti(A, B, C, D))` to 10⁻⁹ relative, and so does the companion realisation of (2); a held random input with a random initial state equals `scipy.signal.lsim(..., interp=False)` in output and state.
6. **Simulation laws**: A_d = `scipy.linalg.expm`, B_d equals a 4001-point quadrature; superposition; a 10-sample input delay delays the output by 10 samples; the free response is e^{At}x₀; the step's final value equals the DC gain.

### Compare

`scipy.signal.lti` is the oracle: its `step` integrates the state-space model with its own discretisation, and `to_ss`/`to_tf` convert with its own algorithms. `main()` prints:

```
max relative |C(sI-A)^-1 B + D - hand TF| at 4 test points = 6.23e-16
1 deg elevator step: h(6 s) = 13.062 m, theta(6 s) = 10.276 deg, alpha(6 s) = 1.3333 deg
max |physical - companion| = 5.51e-14 m, |physical - observable| = 5.86e-14 m, |physical - scipy.signal.step| = 8.70e-14 m
```

Three realisations with different A matrices (the physical one has entries like −V = −30, the companion one has the polynomial coefficients −3.7 and −9 in its first row and 432 in C) produce the same altitude history to 10⁻¹³ m. SciPy agrees to the same level. What the hand derivation adds is the meaning of every coefficient. What the state space adds is α, θ and h along the way, and the starting point for `P09.14`–`P09.16`.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Realisations side by side.** For the 1° elevator step, simulate the physical, companion and observable realisations and SciPy's step for 6 s at dt = 10 ms. Record the largest pairwise differences in h and the final α, θ, h. Predict: identical to rounding, because (5) is exact.
2. **Coefficient sensitivity.** For G = 1/((s + 1)(s + 2)⋯(s + n)) with n = 4, 6, 8, 10, build the companion matrix from the coefficients rounded to 5 significant figures (as if typed into flight software from a printout), and the modal (here diagonal) matrix with its entries rounded the same way. Record the largest pole error and the largest spurious imaginary part. Predict how the error grows with n.
3. **Steady short period.** From the simulation, read α and q after the short period has settled and compare with α_ss and q_ss from §2.1.

Raw data format: CSV tables `realisation, h_6s_m, max_diff_vs_physical_m`; `n, largest_coefficient, pole_error_companion, max_imag_companion, pole_error_modal`; `quantity, simulated, hand`, each with a header naming the model and stating that the values are from simulation.

What you should find: the realisations agree to 5.5 × 10⁻¹⁴ m (companion), 5.9 × 10⁻¹⁴ m (observable) and 8.7 × 10⁻¹⁴ m (SciPy). In the sensitivity study the companion form's pole error is 3.9 × 10⁻¹⁴ at n = 4 and 5.5 × 10⁻¹³ at n = 6, where rounding to five figures happens to change nothing because the coefficients are small integers. At n = 8 (largest coefficient 1.18 × 10⁵) it is **0.64**, and two real poles have become a complex pair with imaginary part 0.30. At n = 10 (largest coefficient 1.28 × 10⁷) it is **2.75**, with imaginary parts up to 2.42. The modal form's error is exactly zero, because rounding −1, −2, … to five figures changes nothing, and in general a rounding error of ε moves a modal pole by about ε. The simulated α(6 s) = 1.3333° equals α_ss = 4/3°.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference |
|---|---|---|---|---|
| h/δ at 4 points | 432/(s²(s² + 3.7s + 9)) | C(sI − A)⁻¹B | relative 6.2 × 10⁻¹⁶ | rounding |
| den from Faddeev–LeVerrier | s⁴ + 3.7s³ + 9s² | [1, 3.7, 9, 4 × 10⁻¹⁶, −8 × 10⁻¹⁶] | `np.poly` (test) | rounding at the double root |
| short-period block | σ = −1.85, ω = 2.3617 | [[−1.85, 2.3617], [−2.3617, −1.85]] | eig(A) | — |
| h(6 s), 1° step | — | 13.062 m | `scipy.signal.step`: 8.7 × 10⁻¹⁴ m | rounding |
| α at 6 s | 4/3° = 1.3333° | 1.3333° | — | short period settled |
| q steady | 1.6°/s | 1.6000°/s (from the ODE steady state) | — | — |

Every route agrees to rounding, so the representation is not in doubt. The model is: the steady pitch rate under a held elevator is the clearest prediction the model gets wrong, and it is wrong for a reason the model states (constant V). Extension 1 adds airspeed and removes it. The second lesson of the validation is numerical. A correct companion form can still be a bad thing to put in software at order 8 or more, because its coefficients span many decades and its poles depend on them violently.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(5); `state_space.png`; the three experiment tables; and a short note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending the coordinates in which an altitude-hold design and its flight software should hold this model, with one number from the sensitivity study and one limit of the physics.

### Adversarial self-checks

- **What if V doubles?** Only row 4 of A changes; h/δ doubles (the test checks it) because the same flight-path angle gives twice the climb rate. In reality M_δ, M_α and Z_α change with dynamic pressure too, so the short period moves. Which entries of A would you rescale?
- **What if the elevator input were replaced by a vertical gust?** A gust changes α directly: a new column of B with Z_α in row 1 (and M_α in row 2). The same A, a different B, and a different TF with the same poles.
- **Which assumption fails first in the real world?** Constant airspeed, within seconds of a held elevator; then small angles; then the linear lift curve.
- **What measurement could falsify the model?** A logged elevator doublet with α, q, θ and h: the short-period frequency and damping in q, and the ratio of climb rate to γ, test rows 1–2 and row 4 separately.
- **Does the result survive a change of realisation or solver?** Yes to 10⁻¹³ m, which is the point of §2.4; it does not survive rounding the coefficients of a high-order companion form.
- **Which decision changes if the parameter uncertainty doubles?** If M_δ is uncertain by ±30 %, h/δ scales by the same factor and the predicted 13 m climb ranges from about 9 to 17 m. An altitude-hold loop must be designed with that range in mind (`P09.18`).

## Check yourself

1. Write ẋ = Ax + Bu, y = Cx + Du for ÿ + 3ẏ + 2y = u, y the output, in controllable canonical form.
2. Why are there infinitely many state-space models with the same transfer function? Give the transformation.
3. Show in one line that the transfer function does not change under x = Tz.
4. What two quantities can you read directly off a 2 × 2 real modal block [[σ, ω], [−ω, σ]]?
5. In the pitch-altitude model, which eigenvalues are at the origin and what physical integrations do they represent?
6. Why does a held 1° elevator produce a pitch rate that never stops in this model, and what is missing?
7. What does Cayley–Hamilton say, and how is it used as a check in Faddeev–LeVerrier?
8. Why is the companion form a poor choice for a 10th-order model in flight software?
