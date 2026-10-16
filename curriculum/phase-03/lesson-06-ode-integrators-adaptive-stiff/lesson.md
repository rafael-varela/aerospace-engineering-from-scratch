# P03.06  ODE Integrators II: Adaptive & Stiff

> An integrator that cannot measure its own error cannot choose its own step; and no explicit method, however accurate, can take a step longer than the fastest time constant in the problem allows.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the embedded error estimate of a Runge–Kutta pair and the step-size rule h_new = h·(tol/err)^{1/5} from the order of the local error, and explain the roles of `rtol`, `atol`, the safety factor and step rejection. 2. Derive the stability function R(z) of a one-step method, locate the real stability interval of explicit Euler, RK4 and the Dormand–Prince pair, and show that implicit Euler is stable for every step on a decaying mode. 3. Implement and test, from scratch, a Dormand–Prince 5(4) integrator with step control and a backward-Euler integrator with a Newton solve, against `scipy.integrate.solve_ivp` (DOP853, Radau) and `scipy.linalg.expm` as oracles, with an observed-order study for each. 4. Diagnose stiffness from a Jacobian's eigenvalues and from the behaviour of an adaptive explicit solver, and decide which solver a given problem needs. |
| **Prerequisites** | `P03.05` (Euler, RK4, local and global error, order of accuracy), `P03.02` (Newton's method), `P02.05` (e^{At}, stability iff Re λ < 0), `P02.04` (eigenvalues of a linear system), `P01.16`. Skill check: you can write the RK4 Butcher table from memory and state why its global error is O(h⁴). |
| **Engineering question** | A 1200 kg ballistic capsule enters at 120 km, 7.8 km/s, −6°: what peak deceleration does it see, at what altitude, and how many right-hand-side evaluations does a tolerance of 10⁻⁶ cost if the integrator chooses its own steps? Then: a sensor skin on a vehicle responds to heating in 0.2 s while the structure under it warms over hours; what is the structure temperature after one hour, and why does the same adaptive integrator need 5 000 steps for that when 60 would do? |
| **Assumptions and limits** | Point-mass planar entry in an exponential atmosphere (ρ = ρ₀e^{−h/H}, H = 7.2 km) with constant drag coefficient, no lift, inverse-square gravity and a round non-rotating Earth; the drag-only straight-line version is the Allen–Eggers model, used as an analytic oracle and not as a prediction. The thermal model is two lumped nodes with constant properties, a constant absorbed heat flow and grey-body radiation from the skin; it says nothing about where the heat comes from. Step control is "error per step" with a scaled RMS norm; the stopping condition is detected after a step, not located by root-finding. |
| **Analytical work** | The embedded error estimate and the step-size law; the stability function of each method by hand; the stability interval of explicit Euler (2) and RK4 (2.785); A-stability of implicit Euler; the Allen–Eggers velocity law and its peak-deceleration formula; the eigenvalues and stiffness ratio of the two-node thermal model; the Newton iteration for the backward-Euler step. |
| **Computational deliverable** | `starter/main.py`: `rk45_step`, `solve_rk45` (with `error_norm`, `initial_step`), `explicit_euler_step`, `rk4_step`, `implicit_euler_step`, `solve_fixed`, `solve_implicit_euler`, `fd_jacobian`, `stability_function`, `stability_interval`, `stiffness_ratio`, `allen_eggers_velocity`, `allen_eggers_peak`, `fixed_step_errors` (the physical models and `convergence_order` are given). Tests prove: the DP5 step is fifth-order and its embedded estimate scales as h⁵; fixed-step DP5, RK4 and backward Euler show observed orders 5, 4 and 1; the stability functions match the hand formulas and the real stability intervals are 2.000, 2.785, 3.307 and ∞; the adaptive driver lands on t₁ exactly, rejects steps at the deceleration peak, concentrates its steps there, and matches DOP853 to 10⁻² m/s on the entry; the drag-only entry reproduces the Allen–Eggers velocity to 10⁻² m/s and the peak formula to 0.5 %; backward Euler matches the matrix exponential on the linear thermal model and the Radau oracle on the nonlinear one, with order 1.00; on the stiff problem the adaptive RK45 is capped at the stability step and costs more than 50 times the work of backward Euler. |
| **Practical activity** | **computational experiment**: a tolerance sweep on the entry (error and work against `rtol` from 10⁻³ to 10⁻¹⁰), a fixed-step convergence study of three methods through the deceleration peak, and a cost comparison of explicit, adaptive-explicit and implicit integration on the stiff thermal model against the five `solve_ivp` methods. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `adaptive_stiff.png`, the hand derivations of the step-size law and of the four stability functions, the tolerance-sweep table, the stiff-cost table, and a half-page note recommending an integrator for a coupled trajectory-plus-thermal simulation with a reason and a cost estimate. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the step-size law and the stability intervals are derived; all public tests pass; the peak-deceleration numbers (Allen–Eggers and with gravity) and the stiffness ratio are reproduced by hand; the note names a solver and gives a step-count argument for it. |
| **References** | `@ocw-16-90` Unit 1 "1.6 Systems of ODE's and Eigenvalue Stability", "1.7 Stiffness and Implicit Methods" (1.7.1 Stiffness, 1.7.3 Implicit Methods for Linear Systems of ODEs) (required-free); `@langtangen-fd-computing-exponential-decay` Chapter 3 "Generalizations": "3.2.10 Example: Runge–Kutta Methods", "3.2.11 Example: Adaptive Runge–Kutta Methods" (required-free: the 4th-minus-5th-order difference as the error indicator and a tolerance sweep); `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" and "Solving a system with a banded Jacobian matrix" (required-free: RK45/DOP853 against Radau/BDF/LSODA, `rtol`/`atol`); `@ocw-16-90` Unit 1 "1.8 Multi-Step Methods" (supplementary-free: Adams and BDF families). |
| **Extensions** | ★ Add the FSAL (first-same-as-last) trick to `rk45_step` and verify the evaluation count drops from 7 to 6 per step. ★ Implement BDF2 with a variable step and compare its cost with backward Euler at equal error on the thermal model. ★ Add dense output (the DP5 quartic interpolant) and locate the 20 km crossing by root-finding instead of stopping after the step. |
| **Time** | unmeasured |

## 1. The engineering problem

`P03.05` built RK4 and showed that its global error falls as h⁴. It left two questions open. First, *which* h? A re-entry flight is quiet for a minute, violent for twenty seconds around peak deceleration, and quiet again; a step small enough for the peak wastes most of its evaluations on the quiet parts, and a step chosen for the quiet parts misses the peak. An integrator must therefore estimate its own error at every step and change the step accordingly. Second, there are problems on which *no* explicit method, adaptive or not, can be made cheap: those with widely separated time scales, called **stiff**. A sensor skin that follows its heating in 0.2 s is bolted to a structure that warms over hours. The skin's transient is over almost at once and the solution is smooth afterwards, yet every explicit method is forced to keep its step below about 0.4 s for the whole hour, not for accuracy but for *stability*. That is 9 000 steps where 60 would do.

This lesson builds both remedies from scratch. The first is an **embedded Runge–Kutta pair** (Dormand–Prince 5(4), the method behind `solve_ivp`'s default `RK45`): two solutions of different order from the same stage evaluations, whose difference estimates the error and drives a step-size controller. The second is **implicit (backward) Euler**: a method whose stability region contains the whole left half-plane, bought at the price of a nonlinear solve at every step. Each is applied to the problem that motivates it: the capsule entry, where the step must follow the physics, and the skin-structure model, where the step must ignore it.

Success criteria, written first: the adaptive integrator reproduces `solve_ivp`'s DOP853 to 10⁻² m/s on the entry at `rtol` = 10⁻⁸ with fewer than 60 steps; the drag-only entry matches the closed-form Allen–Eggers velocity to 10⁻² m/s and the predicted peak deceleration to 0.5 %; on the thermal model backward Euler with 60 s steps is within 0.1 K of a Radau reference after an hour, and the adaptive explicit method is shown to be stability-limited (step capped near 3.3/|λ_fast| regardless of tolerance) and at least 50 times more expensive.

## 2. Understand and derive

### 2.1 Local error, global error, and what a controller can measure

For a one-step method of order p the **local error** (the error of one step started from the exact solution) is C h^{p+1} and the **global error** after T/h steps is O(h^p) (`P03.05`, `@ocw-16-90` 1.3–1.4). A controller cannot measure the global error, which depends on everything that happened before; it can only estimate the local error of the step it has just taken. The strategy that works is to keep each local error below a tolerance and let the global error take care of itself, which it does for a stable method: the global error is a sum of propagated local errors.

The tolerance is applied component-wise through a **scale** s_i = atol + rtol·max(|y_i|, |y_i,new|): `rtol` is the relative accuracy you want, and `atol` is the absolute error below which you stop caring, which matters for any component that passes through or starts at zero. The entry state starts with ground range s = 0; with atol = 10⁻⁹ m the controller tries to resolve the first metre of range to a nanometre and takes a first step of 0.2 μs. Setting atol = 10⁻³ m is the fix, and the lesson uses it throughout. The step is accepted when

**E = sqrt( (1/n) Σ_i (err_i / s_i)² ) ≤ 1**,   (1)

an RMS over the components, which is what `solve_ivp` does.

### 2.2 The embedded estimate and the step-size law

An s-stage explicit Runge–Kutta method computes k_i = f(t + c_i h, y + h Σ_{j<i} a_{ij} k_j) and y_new = y + h Σ b_i k_i. The Dormand–Prince pair has seven stages and **two sets of weights**: b (fifth order) and b̂ (fourth order). The same seven k_i give two answers, y₅ and y₄, and

**err = y₅ − y₄ = h Σ_i (b_i − b̂_i) k_i**   (2)

costs nothing beyond a dot product. Since y₅ is accurate to O(h⁶) and y₄ to O(h⁵), the difference *is* the local error of y₄ up to higher-order terms. The method advances with y₅ (local extrapolation) and uses err, which overestimates y₅'s own error, as a conservative control signal; the tests check that |err| exceeds the true fifth-order error on y' = y and scales as h⁵. The seventh stage is evaluated at t + h with the fifth-order weights, so it equals the first stage of the next step and can be reused (FSAL, extension 1); our implementation recomputes it for clarity and counts seven evaluations per attempt against `solve_ivp`'s six.

Now the controller. If the local error estimate at step h is E (in the scaled norm of (1)) and the estimate is O(h⁵), then at a step h_new it would be E·(h_new/h)⁵. Demanding that it equal 1 gives h_new = h·E^{−1/5}. Three engineering corrections make this usable:

**h_new = h · min(5, max(0.2, 0.9 · E^{−1/5})).**   (3)

The **safety factor** 0.9 aims below the tolerance, so that the next step is accepted more often than not; the **clamps** stop the step growing more than 5× (the estimate E was taken at h and says nothing reliable about 5h) or shrinking more than 5× in one go. If E > 1 the step is **rejected**: nothing is stored, the same formula shrinks h, and the step is retried from the same point. Rejections are not failures; they are the controller discovering that the solution has become harder, which on the entry happens where the drag rises steeply (8 rejections among 28 accepted steps at rtol = 10⁻⁶). The exponent −1/5 is set by the order of the *estimate*: a 2(3) pair would use −1/3. This is the whole of step control; everything else in a production code (PI control, step-size smoothing, dense output) refines it.

Two consequences follow directly from (3). The error of the result is roughly proportional to the tolerance (the sweep in part 4 shows error ∝ rtol^{1.1} over seven decades). And the *work* scales as rtol^{−1/5}: tightening the tolerance by 100 costs only 100^{1/5} = 2.5 times the steps, which is why high-order adaptive methods are cheap at high accuracy and why `rtol` = 10⁻⁸ is not an extravagance.

### 2.3 Stability functions

Apply a one-step method to the scalar **test equation** y' = λy with λ complex. Every explicit RK method gives y_{n+1} = R(z) y_n with z = hλ and R a polynomial; the method is stable (the numerical solution does not grow when the exact one decays) only if **|R(z)| ≤ 1**. The set of z where that holds is the **stability region** (`@ocw-16-90` 1.6). By hand:

- explicit Euler: y + h λ y, **R(z) = 1 + z**; on the negative real axis |R| ≤ 1 for −2 ≤ z ≤ 0: the stability interval has length L = 2;
- RK4: R(z) = 1 + z + z²/2 + z³/6 + z⁴/24 (the Taylor polynomial of e^z, which is what fourth-order accuracy means on this equation); |R(−x)| = 1 again at x = 2.7853, so L = 2.785;
- Dormand–Prince 5(4), fifth-order weights: a sixth-degree polynomial equal to e^z through z⁵ (seven stages give one more coefficient, 1/600 instead of 1/720, for z⁶); L = 3.307;
- backward (implicit) Euler: y_{n+1} = y_n + hλ y_{n+1}, so **R(z) = 1/(1 − z)**, and |R| ≤ 1 for every z outside the disc |z − 1| < 1, which contains the entire left half-plane. A method stable for every Re z < 0 is **A-stable**; implicit Euler also has R → 0 as z → −∞ (**L-stable**), so a fast decaying mode is damped out in one step rather than merely kept bounded.

`stability_function` does not type these in: it *runs* each stepper on y' = zy with h = 1, so the number tested is the stability function of the code as written; `stability_interval` scans the negative axis and bisects the crossing of |R| = 1. The figure in `main()` draws the |R| = 1 contours. For a linear system y' = Ay the same analysis applies mode by mode: with A = VΛV⁻¹ the method acts on each eigen-component with z = hλ_k, so **the step must keep every hλ_k inside the stability region**. That one sentence is the whole theory of stiffness.

### 2.4 Stiffness

A system is **stiff** when the step needed for stability is far smaller than the step needed for accuracy on the part of the solution you care about (`@ocw-16-90` 1.7.1). The usual measure is the **stiffness ratio** max|Re λ|/min|Re λ| of the Jacobian's eigenvalues, though the ratio alone misleads: a problem is stiff only on the interval where the fast modes have already decayed and you still have to integrate.

The thermal model is two lumped nodes: a skin of heat capacity C_s = 500 J/K absorbing q = 2 kW and radiating with emissivity 0.8, bolted through a conductance k_c = 2500 W/K to a structure of C_b = 2 × 10⁵ J/K that leaks to a 300 K sink through k_b = 20 W/K:

**C_s Ṫ_s = q − k_c(T_s − T_b) − εσ(T_s⁴ − T_∞⁴),  C_b Ṫ_b = k_c(T_s − T_b) − k_b(T_b − T_∞).**   (4)

With radiation off it is linear, y' = Ay + c, and `P02.05` gives the exact solution y(t) = y_eq + e^{At}(y₀ − y_eq) with y_eq = −A⁻¹c: the structure settles at T_∞ + q/k_b = 400 K and the skin q/k_c = 0.8 K above it. The Jacobian has eigenvalues −5.022 s⁻¹ and −1.241 × 10⁻⁴ s⁻¹ (close to −k_c/C_s and −k_b/C_b, the time constants of 0.20 s and 2.2 h): stiffness ratio **4.05 × 10⁴**. Explicit Euler is stable only for h < 2/5.022 = 0.398 s, the DP5 pair for h < 3.307/5.022 = 0.658 s. After the first second the fast mode has decayed by e⁻⁵ and the solution is as smooth as any curve can be, but the eigenvalue is still there in the Jacobian and the stability limit still applies: an hour costs at least 9 000 explicit Euler steps or 5 500 RK45 steps. The run confirms it: the adaptive RK45 at rtol = 10⁻⁶ takes **5 472 steps** with a largest step of 0.766 s (it may stray slightly past the stability boundary for a step or two before the error estimate drags it back), 39 307 evaluations, and delivers T_b to 10⁻⁷ K: a thousand times more accuracy than anyone asked for, because the step was never limited by accuracy at all.

### 2.5 Backward Euler with Newton

Backward Euler advances by **y_{n+1} = y_n + h f(t_{n+1}, y_{n+1})**: the slope is taken at the *end* of the step, so y_{n+1} appears on both sides and each step is a system of nonlinear equations G(z) = z − y_n − h f(t_{n+1}, z) = 0. Newton's method (`P03.02`, generalised to systems; `P03.10` does it properly) iterates

**(I − h J(t_{n+1}, z)) Δz = −G(z),  z ← z + Δz,**   (5)

with J = ∂f/∂y, from the start z = y_n, until ‖Δz‖ is below a tolerance. On a linear problem one iteration is exact; on (4) it takes two or three. The Jacobian is supplied analytically (`thermal_jac`) or, failing that, by forward differences (`fd_jacobian`, a column at a time). The cost per step is one Jacobian, one linear solve and a few evaluations, against seven evaluations for DP5, but the step is free of any stability limit: at h = 60 s, z_fast = −301 and R = 1/302, so the skin transient is crushed in one step while the slow mode, with z_slow = −0.0074, is integrated with a relative error of z²/2 per step, i.e. first-order accuracy. The price is accuracy: backward Euler is only first order, and the error in T_b after an hour is 0.087 K at h = 60 s, 0.0087 K at 6 s, 0.00087 K at 0.6 s (observed order 1.00). The remedies are higher-order implicit methods: the trapezoidal rule (second order, A-stable but not L-stable; Problem 5 shows what that costs), the BDF family (`@ocw-16-90` 1.8; `solve_ivp`'s `BDF`) and implicit Runge–Kutta (`Radau`).

### 2.6 The entry model and its analytic oracle

With state (h, v, γ, s) (altitude, speed, flight-path angle positive upward, ground range), ballistic coefficient β = m/(C_D A) and g = g₀(R/(R + h))²:

**ḣ = v sin γ,  v̇ = −ρv²/(2β) − g sin γ,  γ̇ = (v/(R + h) − g/v) cos γ,  ṡ = v cos γ · R/(R + h).**   (6)

Drop gravity and curvature and γ is constant; then dv/dh = v̇/ḣ = −ρ₀e^{−h/H} v/(2β sin γ) separates into the **Allen–Eggers** law v(h) = v_E exp[−(ρ₀H/(2β|sin γ_E|))(e^{−h/H} − e^{−h_E/H})]. The deceleration a = ρv²/(2β), written in u = ρ₀He^{−h/H}/(β|sin γ_E|), is (|sin γ_E| v_E²/2H)·u e^{−u} (for entry from high altitude), maximal at u = 1:

**a_max = v_E² |sin γ_E| / (2eH),  h* = H ln[ρ₀H/(β|sin γ_E|)].**   (7)

a_max does not depend on the ballistic coefficient; h* does. For the lesson's capsule (m = 1200 kg, C_D = 1.2, A = 4.5 m², so β = 222.2 kg/m²) entering at 120 km, 7 800 m/s, −6°: **a_max = 162.5 m/s² = 16.57 g at h* = 42.76 km**. The integrator with `gravity=False` returns 16.557 g at the nearest stored point, 43.0 km, and matches v(h) to 1.8 × 10⁻⁴ m/s. With gravity and curvature on, the peak is **18.12 g at 41.7 km, t = 99.4 s, v = 4 609 m/s** (finely resolved reference): gravity adds to the deceleration along the path and steepens γ, which brings the peak lower and higher. The capsule reaches 20 km after 172 s, 821 km down-range, at 281 m/s. Nothing here is a prediction for a real vehicle: the exponential atmosphere is 20 % off the ISA at 40 km, the capsule has no lift, and the Allen–Eggers assumptions (constant γ, no gravity) are precisely the ones that hold only for steep entries. It is an oracle for the *integrator*, which is what this lesson needs.

### 2.7 Where the methods stop being valid

- **Error estimate.** (2) is asymptotic: it is trustworthy only when h is small enough that the next Taylor term is smaller than the one estimated. A controller started with a huge step produces an estimate that means nothing, which is why `initial_step` begins from the size of f and its first difference, and why the clamps in (3) exist.
- **Tolerance and output.** Step control bounds the local error, not the global one, and only at the step points. The rtol = 10⁻⁴ entry run reports a "peak" of 16.9 g at t = 104.6 s, 1.2 g low and 5 s late, not because the integration is wrong (its velocity error is 5.7 m/s) but because with 29 s steps no stored point lands on the peak. Dense output (extension 3) or a finer `t_eval` fixes that; a tighter tolerance fixes it by accident.
- **Stopping.** The 20 km stop is detected after the step that crosses it; the stored end point is up to one step past the altitude. A root-finder on the interpolant is the proper tool.
- **Stiffness detection.** `stiffness_ratio` needs the Jacobian, which a black-box right-hand side does not expose; the practical symptom is an adaptive explicit solver whose step size stops responding to `rtol`.
- **Newton.** (5) converges from z = y_n when h‖J‖ is not absurdly large relative to the nonlinearity; the radiation term T⁴ can break it for a bad start (Problem 5 shows the explicit predictor failing at h = 60 s from 400 K). Production codes reduce h when Newton fails; ours raises.

## 3. Build the computational model

Open `starter/main.py`. The physical models (`Capsule`, `reentry_rhs`, `SkinStructure`, `thermal_rhs`, `thermal_jac`, `thermal_linear_system`, `thermal_exact_linear`), `entry_state`, `deceleration`, `convergence_order` and `main` are given.

### Step 1: one Dormand–Prince step

The tableau is in the module as `DP_C`, `DP_A`, `DP_B5`, `DP_B4`, with `DP_E = DP_B5 − DP_B4`.

```python
def rk45_step(f, t, y, h):
    k = np.zeros((7, y.size))
    k[0] = f(t, y)
    for i in range(1, 7):
        k[i] = f(t + DP_C[i] * h, y + h * (DP_A[i, :i] @ k[:i]))
    return y + h * (DP_B5 @ k), h * (DP_E @ k)        # y5, and y5 - y4
```

`rk4_step` and `explicit_euler_step` are the `P03.05` steps, three and one lines. Keep the arrays' dtype: `stability_function` feeds complex values through these same functions.

### Step 2: the controller

`error_norm` is (1). `initial_step` follows the standard recipe in its docstring (one percent of the solution per step, then corrected by an estimate of the second derivative). `solve_rk45` is the loop of 2.2: clip h to land on t₁, attempt a step, compute E, accept or reject, update h by (3), store accepted steps in `h`, count `n_accepted`, `n_rejected`, `nfev`, and stop after an accepted step if `stop(t, y)` is true. Return an `IVPResult`.

### Step 3: backward Euler

`implicit_euler_step` is (5) with `np.linalg.solve`, starting from z = y, and returns the Newton count; `fd_jacobian` perturbs one component at a time by 10⁻⁷ max(1, |y_j|). `solve_fixed` drives any stepper with equal steps and understands the two return conventions; `solve_implicit_euler` wraps it and counts evaluations.

### Step 4: stability and stiffness

`stability_function` runs the named stepper once on f(t, y) = zy with h = 1 from y = 1; for implicit Euler return 1/(1 − z) directly (Newton on a linear problem is exact in one step, so this is what the code computes). `stability_interval` scans x = 0.01, 0.02, … and bisects the first x with |R(−x)| > 1; `inf` if there is none up to 50. `stiffness_ratio` uses `np.linalg.eigvals`.

### Step 5: the oracle and the studies

`allen_eggers_velocity` and `allen_eggers_peak` are (7) and the law above it. `fixed_step_errors` runs `solve_fixed` at each step count and returns the error at t₁ against an `exact` callable, for all components or one.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-06-ode-integrators-adaptive-stiff`. The tests:

1. **One step.** On y' = y the DP5 step is within 2h⁶ of e^h and its estimate within 2h⁵, larger than the true error; on the oscillator the estimate falls by 32 ± 15 % per halving of h.
2. **Backward Euler step** equals y/(1 − hλ) on a linear problem in at most two Newton iterations, with the analytic or the finite-difference Jacobian; it raises when Newton is denied enough iterations; it stays *above* a decaying exponential (1/(1 + h) > e^{−h}).
3. **Observed orders** on y'' = −y over one period: DP5 fixed-step 4.95, RK4 3.98; backward Euler on y' = −y: 1.0 within ±0.1.
4. **Stability.** R(z) of explicit Euler, RK4 and implicit Euler equal the hand formulas at z = −0.7 + 0.4i, DP5 agrees with e^z through z⁵; intervals 2.000, 2.7853 ± 0.002, (3.2, 3.4), ∞. Explicit Euler on the thermal model survives 200 steps at 0.8 × 0.398 s and blows up at 1.2 × 0.398 s.
5. **Driver.** Lands on t₁ to 10⁻¹²; `h` equals `diff(t)`; the error on y' = −y falls and the step count rises monotonically with rtol, with steps ∝ rtol^{−0.12…−0.30}; the entry forces rejections and stops just below 20 km; more than `max_steps` attempts raises.
6. **Entry.** Our RK45 at rtol = 10⁻⁸ agrees with DOP853 (rtol 10⁻¹²) to 10⁻² m/s and 1 m at every stored time; the peak is 17.5–18.5 g between 40 and 44 km; after the start-up steps the smallest step lies within 40 s of the peak and the largest exceeds three times the smallest. Drag-only: velocity to 10⁻² m/s of Allen–Eggers, peak within 0.5 % of (7) and 1.5 km of h*, a_max unchanged and h* lowered when the mass doubles.
7. **Thermal.** `thermal_linear_system` agrees with `thermal_rhs` and `thermal_jac` at ε = 0; backward Euler with 6 s steps is within 0.02 K of e^{At} after an hour; the equilibrium is T_∞ + q/k_b and q/k_c above it. Nonlinear: within 0.2 K of Radau (rtol 10⁻¹²) with 60 s steps, order 1.0 ± 0.1 across 60, 120, 240 steps. Adaptive RK45 over 1 800 s never exceeds 1.2 × 3.307/|λ_fast| per step and takes more than 50 times the steps and evaluations of backward Euler with 30 steps, for the same T_b to 0.1 K.
8. `main` writes the figure.

### Compare

`solve_ivp` with `RK45` is the same Dormand–Prince pair with the same error norm; the differences are FSAL (six evaluations per step, not seven), a slightly different initial step and clamps, and a quartic interpolant for dense output. On the 300 s entry arc at rtol = 10⁻⁶ it takes 63 steps and 422 evaluations against our 66 steps and 527 evaluations; at 10⁻⁸, 138 against our 142. `DOP853` is an eighth-order pair that costs twelve stages but wins at tight tolerances (2 832 steps against RK45's 5 469 on the stiff problem, both stability-bound). For the stiff thermal model the implicit methods are in a different league at rtol = 10⁻⁶: `Radau` 17 steps and 127 evaluations, `BDF` 50 and 102, `LSODA` 50 and 108, against 38 000 evaluations for `RK45`. Our backward Euler with 60 steps and 180 evaluations sits between them in cost and far below them in accuracy (0.087 K against 10⁻⁴ K), because it is first order; what the library adds is order (Radau is fifth order, BDF up to fifth), variable step and order selection, and Newton with a reused, factored Jacobian. We still call the library in practice for exactly those reasons, and because its failure modes are documented. We write our own once so that "the solver took 5 000 steps on a smooth curve" reads as a diagnosis and not a mystery.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Tolerance sweep.** Integrate the entry from 0 to 300 s with rtol = 10⁻³ … 10⁻¹⁰ (atol = 10⁻³ rtol) and record accepted steps, rejections and the error in v(300 s) against a DOP853 reference at rtol 10⁻¹³. Predict from 2.2: error roughly ∝ rtol, steps ∝ rtol^{−1/5} (a factor 2.5 per two decades).
2. **Fixed-step convergence through the peak.** With `fixed_step_errors` and the same reference, run RK4, fixed-step DP5 and backward Euler with 100, 200, 400 and 800 steps over 0–300 s. Predict orders 4, 5 and 1; predict which will be furthest from its asymptotic order at these step sizes (the peak is 20 s wide and h = 3 s at n = 100).
3. **Stiff cost table.** On the thermal model over one hour: explicit Euler at h = 0.3 s and 0.5 s (predict: stable and unstable); adaptive RK45 at rtol = 10⁻⁶ (predict the step count from 3 600 × 5.022/3.307 ≈ 5 500); backward Euler with h = 60, 6, 0.6 s (predict error ∝ h); the five `solve_ivp` methods at rtol = 10⁻⁶. Record steps, evaluations and the error in T_b against Radau at 10⁻¹².
4. **Where does the step go?** Plot the accepted step size against time for the entry at rtol = 10⁻⁶ (the top-left panel of `adaptive_stiff.png`). Predict: growth by 5× per step from the start-up value, a minimum near t ≈ 100 s, then slow recovery as the deceleration decays.

Raw data format: CSV tables (`rtol, n_accepted, n_rejected, nfev, err_v`; `method, n_steps, err_v`; `method, h_or_rtol, n_steps, nfev, err_Tb`), each with a header naming the fixed parameters.

What you should find: the sweep gives 29, 36, 48, 66, 95, 142, 214, 331 steps and errors 0.28, 1.6 × 10⁻², 8.5 × 10⁻⁴, 8.8 × 10⁻⁵, 6.7 × 10⁻⁶, 5.8 × 10⁻⁷, 4.7 × 10⁻⁸, 3.9 × 10⁻⁹ m/s for rtol = 10⁻³ … 10⁻¹⁰, so steps grow by 1.2–1.6 per decade (the −1/5 law predicts 1.58) and the error falls by about a decade per decade of tolerance; the fixed-step study gives observed orders 4.53 (RK4, errors 7.9 × 10⁻⁵ to 6.4 × 10⁻⁹ m/s), 5.52 (DP5, 1.1 × 10⁻⁵ to 1.2 × 10⁻¹⁰) and 1.17 (backward Euler, 0.13 to 0.01 m/s), all above the asymptotic value because the error constant is dominated by the peak and the coarsest runs are not yet asymptotic (the same methods on the oscillator in the tests give 3.98, 4.95 and 1.0); explicit Euler at 0.3 s holds 301.4 K after 60 s and at 0.5 s overflows to NaN; adaptive RK45 takes 5 472 steps and 39 307 evaluations; backward Euler errors are 8.66 × 10⁻², 8.69 × 10⁻³, 8.69 × 10⁻⁴ K (order 1.00) at 180, 1 202 and 12 002 evaluations.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| stability interval, explicit Euler / RK4 / DP5 | 2 / 2.7853 / — | 2.0000 / 2.7853 / 3.3066 | < 10⁻⁴ | bisection |
| thermal eigenvalues (s⁻¹) | ≈ −k_c/C_s = −5.0, −k_b/C_b = −1.0 × 10⁻⁴ | −5.0223, −1.2414 × 10⁻⁴ | 0.4 %, 24 % | the linearised radiation 4εσT³ = 4.9 W/K adds to k_c for the fast mode and adds a quarter to the 20 W/K leak for the slow one |
| explicit Euler limit (s) | 2/5.022 = 0.398 | stable at 0.3, NaN at 0.5 | — | consistent |
| RK45 step on the stiff model (s) | ≤ 3.307/5.022 = 0.658 | max 0.766 | one step 16 % over | error control, not stability, rejects the next step |
| a_max, drag only (g) | 16.567 | 16.557 | 0.06 % | stored points miss the exact peak |
| h* (km) | 42.76 | 43.02 | 0.26 km | same |
| v(h), drag only (m/s) | Allen–Eggers | max dev. 1.8 × 10⁻⁴ | — | rtol 10⁻⁸ |
| peak with gravity (g, km) | — | 18.12, 41.7 | — | gravity adds g sin γ and steepens the path |
| v vs DOP853, rtol 10⁻⁸ (m/s) | — | 3.4 × 10⁻⁴ | — | tolerance |
| T_b after 1 h, ε = 0 (K) | e^{At}: 330.168 | backward Euler h = 6 s: 330.161 | 0.007 | first-order lag |
| T_b after 1 h, ε = 0.8 (K) | Radau: 328.765 | backward Euler h = 60 s: 328.678 | 0.087 | first-order lag ∝ h |
| backward Euler order | 1 | 1.00 | 0 | asymptotic |
| RK45 steps on the stiff model | ≈ 3600/0.658 = 5 470 | 5 472 | 0.04 % | stability-bound |

Two validations are worth separating. The *integrators* are validated against closed forms (Allen–Eggers, e^{At}) and against SciPy's independent implementations, and the agreement is at the tolerance level in every row: the code does what the theory says. The *models* are not validated at all by this, and the table says so: the entry peak of 18 g belongs to an exponential atmosphere and a liftless capsule, and the thermal model's 0.087 K error is small against the uncertainty in q, k_c or ε, which a real design would carry at 10–30 %. That comparison, numerical error against model uncertainty, is the one that decides how tight a tolerance to pay for: for the thermal problem h = 60 s and first order are already better than the inputs deserve; for the entry, a guidance study that differentiates the trajectory with respect to γ_E needs rtol = 10⁻⁸ or better so that the differences are not noise.

Improvements worth making: FSAL (extension 1), which is free; a second-order L-stable implicit method (BDF2, extension 2) for the thermal model, which buys two orders of magnitude in accuracy at the same step; dense output and event location for the entry (extension 3). Improvement worth *not* making: pushing rtol below 10⁻¹⁰ on the entry, where the 331-step run is already within 4 × 10⁻⁹ m/s of the reference and rounding in the 10⁶ m range coordinate starts to show.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (2), (3), the four stability functions and (7); `adaptive_stiff.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending an integrator for a simulation that couples the entry trajectory to the skin-structure thermal model, with a step-count estimate for each candidate.

### Adversarial self-checks

- **What if the entry angle doubles to −12°?** a_max doubles (it is ∝ |sin γ_E|) to about 33 g and h* falls by H ln 2 ≈ 5 km; the peak gets narrower in time, so the controller takes more rejections and smaller steps there. Check with `entry_state(gamma_e_deg=-12)`.
- **What if the skin capacity goes to zero?** λ_fast → −∞ and every explicit method's stable step → 0, while backward Euler is untouched (R → 0). In the limit the skin equation becomes algebraic, q = k_c(T_s − T_b) + radiation, and the system is a differential-algebraic equation: the stiff limit of an ODE *is* a DAE, and implicit methods are the ones that survive the limit.
- **Which assumption fails first in the real world?** For the entry, the exponential atmosphere (the real density at 40 km differs by 20 %), then the absence of lift (any real capsule has some); for the thermal model, the constant q and the lumping of the structure into one node (a thick structure has its own internal time constants, which add more stiff modes, not fewer).
- **What measurement could falsify the model?** An onboard accelerometer trace: a peak 20 % away from 18 g or 5 s late at the predicted altitude indicates the atmosphere; a skin thermocouple reaching equilibrium with the structure in seconds rather than the predicted 0.2 s indicates k_c.
- **Does the result survive a change of solver?** Our RK45, DOP853, RK4 at 800 steps and fixed-step DP5 agree on v(300 s) to 10⁻⁸ m/s; Radau, BDF, LSODA and backward Euler agree on T_b to 0.1 K. Explicit Euler at h = 0.5 s does not survive at all, and that is the point of the lesson.
- **Which decision changes if the tolerance doubles?** None: rtol from 10⁻⁶ to 2 × 10⁻⁶ changes the entry velocity by 10⁻⁴ m/s and the step count by 15 %. If the *model* uncertainty doubles (say ε from 0.8 ± 0.1 to ± 0.2) the thermal answer moves by kelvins, so the decision to use a first-order method with 60 s steps, already far below the model error, is unchanged, and the decision that *does* change is how much to spend on measuring ε.

## Check yourself

1. Why does the Dormand–Prince pair advance with the fifth-order solution if the error estimate belongs to the fourth-order one? What is this choice called?
2. A controller uses h_new = h·E^{−1/5}. Why −1/5 and not −1/4 or −1/6, and what changes for a 2(3) pair?
3. Write R(z) for explicit Euler, RK4 and implicit Euler, and give the real stability interval of each.
4. What is the stiffness ratio of the thermal model, and what step would explicit Euler need over an hour? How many steps is that?
5. State the two Allen–Eggers results (7) and say which of them depends on the ballistic coefficient.
6. Give one symptom, visible from an adaptive explicit solver's output alone, that a problem is stiff.
7. Why does the backward-Euler Newton iteration start from z = y_n rather than from the explicit Euler predictor?
8. Why is `atol` needed at all, and what went wrong in this lesson when it was set to 10⁻⁹ m?
