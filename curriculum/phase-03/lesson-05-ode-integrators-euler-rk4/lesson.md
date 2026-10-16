# P03.05  ODE Integrators I: Euler to RK4

> Every trajectory you will ever compute is a sequence of small steps, each one a guess at where the slope points. The art is in how many slopes you sample per step, and in knowing what the guess costs.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the forward Euler, midpoint, Heun and classical Runge–Kutta methods from Taylor expansions, state their local and global orders, and show why a one-step method's global order is one less than its local order. 2. Derive the stability function R(z) = Σz^k/k! of the explicit RK methods, locate their real-axis limits (2, 2, 2.51, 2.79) and the imaginary-axis limit 2√2 of RK4, and predict from them when a fixed-step integration blows up. 3. Implement the four methods from scratch, both as hand-written steps and through a Butcher tableau, verify the order conditions, and measure local and global orders against exact solutions and `scipy.integrate.solve_ivp`. 4. Propagate a two-body orbit with each method, measure position error, energy and angular-momentum drift against the Kepler solution of `P03.02`, and size the step for a stated position tolerance from one measured error and the method's order. |
| **Prerequisites** | `P02.04` (first-order systems, the state vector), `P02.01` (what an initial-value problem is; exact solutions used as references), `P03.01` (truncation against rounding, observed order), `P03.02` (Newton's method for Kepler's equation, reused here as the exact orbit), `P03.04` (Simpson's rule: RK4 reduces to it; Richardson-style error estimation), `P01.05` (Taylor series in two variables). Skill check: you can write a second-order ODE as a first-order system and expand f(t + h, y + k) to second order. |
| **Engineering question** | A satellite in a 379 × 1879 km orbit (a = 7500 km, e = 0.1, period 6464 s) must be propagated numerically for one orbit with a position error below 1 m, and for ten orbits without an energy drift that turns the orbit into a spiral. Which fixed-step method, and how many steps per orbit, does that take? Forward Euler, with one slope per step, needs about 3 × 10⁹ steps per orbit; RK4, with four, needs 210. The question is how to know that before running, from a single coarse run and the method's order, and how to recognise, from the energy history, which way the numerical orbit is drifting. |
| **Assumptions and limits** | Explicit one-step methods with a fixed step on a smooth right-hand side (f Lipschitz with bounded derivatives up to the method's order): the error laws hold only then, and fail at discontinuities (thrust on/off) or where the solution is stiff (part 2.6). Planar two-body dynamics with a point-mass Earth and no perturbations: the exact Kepler orbit is the reference, so every error measured here is numerical. Fixed steps are the wrong tool for high-eccentricity orbits and stiff problems, which is the subject of `P03.06`; conservation of energy over long times is the subject of `P03.07`. |
| **Analytical work** | Taylor derivation of Euler's local error h²y″/2 and of the order-2 conditions for two-stage methods; the stability function and its limits; RK4 on y′ = g(t) as Simpson's rule; the local-to-global error argument (n = T/h steps each with error Ch^{p+1}, amplified by at most e^{LT}); the exact two-body state from Kepler's equation; the step h = h₀(tol/e₀)^{1/p} for a tolerance. |
| **Computational deliverable** | `starter/main.py`: `TABLEAUS` (given) and `explicit_rk_step`, `euler_step`, `midpoint_step`, `heun_step`, `rk4_step`, `integrate`, `stages`, `observed_order`, `local_error`, `global_error`, `step_for_tolerance`, `stability_function`, `stability_limit_real`, `max_stable_step`, `two_body_rhs`, `specific_energy`, `angular_momentum`, `propagate_orbit`, `orbit_position_error`, `energy_drift`; `solve_kepler`, `kepler_state`, `orbit_period` are given. Tests prove: every tableau satisfies the Butcher conditions through its order; the hand-written steps equal the tableau step to 10⁻¹⁴; one step on y′ = y is the truncated exponential and RK4's local error is h⁵/120 within 15 %; RK4 is exact for quartic y(t) and not quintic; `integrate` lands on t1 with the right shapes and rejects bad input; RK4 at h = 10⁻³ matches `solve_ivp(DOP853, rtol=10⁻¹²)` on the Van der Pol oscillator to 10⁻⁸; local orders p + 1 within 0.1 and global orders p within 0.12 for all five methods on the harmonic oscillator; `step_for_tolerance` hits a 10⁻⁸ target within 30 %; the stability limits 2, 2, 2.5127, 2.7853 within 10⁻³ and the imaginary-axis limit 2√2; a decay with λ = −1000 decays at 0.95h_max and blows up at 1.05h_max for Euler, Heun and RK4; the Kepler state conserves energy and angular momentum to 10⁻¹² and matches DOP853 to 1 mm after half an orbit; RK4's orbit error converges at order 4 ± 0.3 and is below 1 cm at 800 steps per orbit while Euler's is above 1000 km; RK4's energy error is negative and grows linearly, Euler's is positive and growing; the predicted step for 1 m gives an error between 0.5 and 2 m and beats Heun at equal cost by 10³. |
| **Practical activity** | **computational experiment**: a convergence and work–precision study of five explicit methods on the two-body problem, a stability experiment on a stiff decay and on the undamped oscillator, and a ten-orbit drift study of energy and angular momentum. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `rk4_orbit.png`, the order and work–precision tables, the derivations of the Euler and RK2 error terms and of the stability function, and a half-page note recommending a method and step for a one-day orbit propagation with 10 m accuracy, with the evidence. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the order-2 conditions and the stability function are derived; all public tests pass; the convergence study reports observed orders with the pre-asymptotic Euler case explained; the note states a method, a step and the measured basis for both. |
| **References** | `@ocw-16-90` Unit 1 "Numerical Integration of ODEs": "1.2 Discretizing ODEs" (1.2.4 The Forward Euler Method, 1.2.5 The Midpoint Method), "1.3 Order of Accuracy", "1.4 Convergence", "1.9 Runge-Kutta Methods" (required-free): local against global error, consistency, the RK4 Butcher table; "1.5 Zero Stability and the Dahlquist Equivalence Theorem" (supplementary-free): why consistency plus zero-stability gives convergence. `@ocw-18-330` "Lecture notes 11–15: Methods for Ordinary Differential Equations" (required-free): derivations of Euler, RK2 and RK4 and the stability analysis. `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (required-free): the library interface and its methods. |
| **Extensions** | ★ Implement Ralston's RK2 (c₂ = 2/3) and show it has the smallest local error constant among two-stage methods on y′ = −2ty. ★ Plot the full stability regions |R(z)| = 1 for orders 1–8 (the series grows to s = p only up to 4: look up the 6-stage order-5 tableau). ★ Propagate the orbit with a J₂ perturbation and compare RK4's secular drift of the argument of perigee with the analytic rate. |
| **Time** | unmeasured |

## 1. The engineering problem

Every dynamics problem from here on, orbits, re-entries, attitude motion, flexible structures, chemical kinetics in a combustor, ends in an initial-value problem y′ = f(t, y), y(t₀) = y₀, with y a state vector of anything from 2 to 10⁶ components. `P02.01`–`P02.05` solved the linear ones exactly. The real ones have no closed form, and a computer advances them in steps: from (tₙ, yₙ) to (tₙ₊₁, yₙ₊₁) using evaluations of f. Forward Euler uses one slope per step; the classical fourth-order Runge–Kutta method uses four, combined so that the step matches the Taylor series of the true solution through h⁴. This lesson builds both, and the two-stage methods between them, and asks the only question that matters in practice: for a required accuracy, how many evaluations of f does each method cost, and how do you know your answer is right when there is no exact solution to compare with?

The test problem is a planar two-body orbit, chosen because `P03.02` gives its exact solution at any time through Kepler's equation. A 7500 km semi-major axis and eccentricity 0.1 put the perigee at 379 km altitude with a period of 6464 s. Propagating it one orbit with 400 RK4 steps (1600 evaluations of the gravity law) gives a position error of 7.9 cm; the same 1600 evaluations spent on Euler give 2300 km, which is to say no orbit at all. The ratio is not a matter of constants but of exponents: Euler's error falls like h, RK4's like h⁴, so each halving of the step buys Euler one bit and RK4 four. The lesson's job is to make that sentence quantitative and to give you the tools (local and global order, the stability function, the drift of conserved quantities, the step-for-tolerance formula) to decide a method and a step before the expensive run.

Success criteria, written first: each method shows its theoretical local order p + 1 and global order p within 0.1 on a problem with an exact solution; the stability limits are reproduced to three figures and predict blow-up to within 5 % in step; the orbit error converges at order 4 for RK4 and the predicted step for a 1 m error delivers between 0.5 and 2 m; the note at the end sizes a one-day propagation.

## 2. Understand and derive

### 2.1 Forward Euler and its two errors

The Taylor expansion y(t + h) = y + hy′ + h²y″/2 + … with y′ = f(t, y) gives, after dropping everything beyond the first order, **yₙ₊₁ = yₙ + hf(tₙ, yₙ)**. Started from the exact solution, one step is off by the dropped terms: the **local error** is h²y″(ξ)/2 = O(h²). After n = T/h steps the errors accumulate; each is carried forward and amplified by at most (1 + hL) per step, L the Lipschitz constant of f, so the **global error** at T is bounded by (e^{LT} − 1)/L × max|y″|h/2 = O(h): one order lower than the local error, because the number of steps is proportional to 1/h. This is the general rule for one-step methods, *local order p + 1, global order p* (16.90 §1.3–1.4), and both are measured in `main()`: on the harmonic oscillator y″ = −y, one Euler step has observed order 2.00 and the error at t = 2π has order 1.08 (5.1 × 10⁻² with 400 steps).

The bound e^{LT} is pessimistic for most problems but exact in spirit: errors made early are propagated by the dynamics, and an unstable system (L large and positive) amplifies them exponentially while a dissipative one damps them. For an orbit, L is of order the mean motion and the errors grow roughly linearly in time (part 4, item 4).

### 2.2 Two slopes: midpoint and Heun

A general two-stage explicit method is k₁ = f(t, y), k₂ = f(t + c₂h, y + c₂hk₁), yₙ₊₁ = y + h(b₁k₁ + b₂k₂). Expanding k₂ in two variables, k₂ = f + c₂h(f_t + f_y f) + O(h²), and comparing with y(t + h) = y + hf + h²(f_t + f_y f)/2 + O(h³), the step matches through h² exactly when

**b₁ + b₂ = 1, b₂c₂ = ½**   (1)

(Problem 7 does this in full). Every solution is second order: c₂ = ½, b = (0, 1) is the **midpoint method** (16.90 §1.2.5); c₂ = 1, b = (½, ½) is **Heun's method** (the trapezoid rule on the slopes); c₂ = ⅔, b = (¼, ¾) is Ralston's. No choice reaches third order: the h³ terms contain four elementary differentials and two free parameters cannot match them all (Problem 7c). Observed: local order 3.00, global 2.00 for both midpoint and Heun, with identical errors 2.6 × 10⁻⁴ at 400 steps on the oscillator (their error constants differ in general, not on this problem).

### 2.3 Four slopes: the classical Runge–Kutta method

k₁ = f(t, y), k₂ = f(t + h/2, y + hk₁/2), k₃ = f(t + h/2, y + hk₂/2), k₄ = f(t + h, y + hk₃), **yₙ₊₁ = y + (h/6)(k₁ + 2k₂ + 2k₃ + k₄)**.

The weights are Simpson's: on y′ = g(t) the stages are g at t, t + h/2 (twice) and t + h and the step is exactly Simpson's rule, exact for cubic g, i.e. quartic y (the tests check this). On y′ = λy the step is y(1 + z + z²/2 + z³/6 + z⁴/24), z = hλ, the exponential series through fourth order, so the local error is h⁵λ⁵y/120 and the method is fourth order (local 5.00, global 4.00 observed; 3.2 × 10⁻⁹ at 400 steps). In general a method's coefficients are collected in a **Butcher tableau** (c | A ; b) and the order conditions are polynomial identities in them: Σbᵢ = 1 (order 1), Σbᵢcᵢ = ½ (order 2), Σbᵢcᵢ² = ⅓ and Σbᵢaᵢⱼcⱼ = ⅙ (order 3), and four more for order 4 (16.90 §1.9); `explicit_rk_step` runs any tableau and the tests verify all eight conditions for RK4 and the appropriate subsets for the others. Four stages for order four is the last free lunch: order 5 needs six stages, order 8 eleven (the Butcher barriers), which is one reason RK4 is the workhorse and why `P03.06`'s RK45 pairs a fourth- and a fifth-order method sharing six stages.

### 2.4 From one measured error to the step you need

If the global error is e₀ = Ch₀ᵖ at step h₀, then the step for a tolerance tol is

**h = h₀(tol/e₀)^{1/p}**,   (2)

and the cost is s(T/h) evaluations with s stages. For the orbit at h₀ = T/800: RK4 e₀ = 4.7 mm, so 1 m needs h = h₀(1/0.0047)^{1/4} = 3.8h₀: 209 steps per orbit (836 evaluations); a check run at 210 steps gives 1.15 m. Heun at the same h₀ has e₀ = 4.8 km, so 1 m needs h₀/69: 55 000 steps, 110 000 evaluations. Euler: 3 × 10⁹ steps. Formula (2) is only as good as the assumption that the error is already in its hᵖ regime, which the Euler row of the orbit table violates (its apparent order over 100–800 steps is 0.77, because 7 × 10⁶ m errors mean the integration never found the orbit; at 3200–12 800 steps the order is 0.99). Always check the order before trusting the extrapolation.

### 2.5 Linear stability

Apply a method to y′ = λy, Re λ < 0, whose solution decays. Each step multiplies y by the **stability function** R(z), z = hλ; the numerical solution decays if and only if |R(z)| ≤ 1. For the s = p explicit methods, R(z) = Σ_{k≤p}z^k/k!, the truncated exponential, so stability is the set where the truncated series stays inside the unit disc. On the negative real axis the limits are |z| ≤ 2 for Euler (1 + z between −1 and 1) and for both RK2 methods (1 + z + z²/2 touches 1 at z = −2), 2.5127 for Kutta's third-order rule and **2.7853 for RK4**; on the imaginary axis Euler is unstable everywhere (|1 + iθ| > 1: an undamped oscillator integrated with Euler spirals outward, by a factor √(1 + θ²) per step), RK2 likewise, and RK4 is stable up to |z| = 2√2 (16.90 §1.5 treats the related zero-stability of multistep methods). `main()` confirms the limits by integrating y′ = −1000y: Euler with h = 0.0019 gives |y(0.1)| = 1.7 × 10⁻³ and with h = 0.0021 gives 47; RK4 with 0.00265 gives 1.4 × 10⁻⁴ and with 0.00292 gives 44. The transition is sharp and it is set by the eigenvalue of largest magnitude, not by the accuracy you want.

### 2.6 Stiffness, in one paragraph

A problem is **stiff** when the stability limit, not the accuracy requirement, dictates the step: y′ = −1000(y − cos t) has a solution that, after a 10⁻³ s transient, is cos t to 10⁻⁶ and would be integrated to that accuracy by RK4 with h = 0.1, but stability forces h < 0.0028 and 3600 steps per unit time, a factor 36 wasted. Chemical kinetics, structural dynamics with stiff modes, boundary-layer equations and any problem with widely separated time scales are stiff; the cure is an *implicit* method whose stability region contains the whole left half-plane, which is `P03.06`. The diagnosis tool is already here: estimate the largest |λ| of the Jacobian and compare h_max = 2.79/|λ| with the step accuracy would allow.

### 2.7 The two-body test problem

With state y = [x, y, vx, vy], the right-hand side is ẍ = −μx/r³, ÿ = −μy/r³. The exact solution from perigee at t = 0 follows from the eccentric anomaly E, solved from M = nt = E − e sin E by Newton's method (`P03.02`): x = a(cos E − e), y = b sin E with b = a√(1 − e²), and the velocities from Ė = na/r. Energy ε = v²/2 − μ/r = −μ/2a and angular momentum ℓ = xvy − yvx = √(μa(1 − e²)) are conserved exactly; a numerical method conserves neither, and the sign and growth of their drift are a fingerprint of the method. RK4, applied to a Hamiltonian system, is slightly **dissipative**: its energy error is negative and grows linearly (−2.3 × 10⁻⁹ after ten orbits at 400 steps per orbit; the first-orbit value times ten). Euler is anti-dissipative: the energy grows (+49 % after ten orbits, the satellite spiralling out), which is the stability picture of 2.5 in disguise (the orbit's eigenvalues are imaginary and Euler amplifies every imaginary mode). Heun sits between, +1.4 × 10⁻⁴. None of the three is **symplectic**: `P03.07` builds the integrators whose energy error stays bounded for ever.

### 2.8 Where the error laws fail

- A discontinuous f (a thrust switch, a stage separation, a contact) breaks every order to 1 unless the step lands on the event; the fix is event location (`P03.06`).
- Stiff eigenvalues make the step a stability matter (2.6).
- Rounding: at very small h the global error of RK4 reaches a floor near ε(T/h) and starts to grow; for a one-orbit integration with 10⁶ steps the floor is 10⁻¹⁰ relative, well below anything here, but a 10-year propagation at 1 s steps meets it.
- The step-for-tolerance formula assumes the asymptotic regime (2.4).
- Energy drift is a *separate* failure from position error: a method can have a tiny position error after one orbit and still spiral after a thousand.

## 3. Build the computational model

Open `starter/main.py`. `TABLEAUS`, `ButcherTableau`, `solve_kepler`, `kepler_state` and `orbit_period` are given; the integrators are yours.

### Step 1: the generic step and the four hand-written ones

`explicit_rk_step(f, t, y, h, tableau)` loops over the stages: kᵢ = f(t + cᵢh, y + hΣⱼ<ᵢ aᵢⱼkⱼ), then y + hΣbᵢkᵢ. Then write `euler_step`, `midpoint_step`, `heun_step`, `rk4_step` by hand with their k's; the tests compare them with the tableau version to 10⁻¹⁴. `STEPPERS` maps the names to the functions (Kutta's third-order rule goes through the tableau).

### Step 2: integrate

`integrate(f, t_span, y0, h, method)` computes n = ⌈(t₁ − t₀)/h⌉ and the adjusted step (t₁ − t₀)/n, allocates `Y` of shape (n + 1, dim), steps, and sets `t[-1] = t1` exactly. Raise `ValueError` for h ≤ 0, t₁ ≤ t₀ or an unknown method.

### Step 3: errors and orders

`local_error` steps once from `exact(t0)` and compares with `exact(t0 + h)`; `global_error` integrates and compares at t₁; `observed_order` is the log-log slope; `step_for_tolerance` is (2).

### Step 4: stability

`stability_function(z, order)` sums the series (works on arrays). `stability_limit_real(order)` scans x ∈ [0, 4] for the first |R(−x)| > 1 and bisects; `max_stable_step(order, lam)` divides by |λ|.

### Step 5: the orbit

`two_body_rhs`, `specific_energy`, `angular_momentum` (both must accept a stack of states, shape (n, 4)), `propagate_orbit(a, e, n_orbits, steps_per_orbit, method)`, `orbit_position_error` against `kepler_state(t[-1], a, e)`, `energy_drift`.

Running `main()` prints

```
harmonic oscillator over one period: local order (one step) and global order (error at t = 2 pi)
   method stages   local  global   err h=2pi/400
    euler      1    2.00    1.08        5.06e-02
 midpoint      2    3.00    2.00        2.58e-04
     heun      2    3.00    2.00        2.58e-04
   kutta3      3    4.00    3.00        1.01e-06
      rk4      4    5.00    4.00        3.19e-09
RK4 on y' = y, h = 0.5: 1.648437500000; 1 + h + h^2/2 + h^3/6 + h^4/24 = 1.648437500000; exp(0.5) = 1.648721270700; local error -2.84e-04 vs h^5/120 = 2.60e-04
real-axis stability limits |h lambda|: order 1: 2.0000, order 2: 2.0000, order 3: 2.5127, order 4: 2.7853
orbit a = 7500 km, e = 0.1: perigee altitude 379 km, period 6464.0 s
Kepler reference vs solve_ivp(DOP853, rtol 1e-13) after one orbit: 1.5e-05 m
 steps/orbit         euler      midpoint          heun        kutta3           rk4   (position error after 1 orbit, m)
         100     1.855e+07     6.460e+04     3.348e+05     1.192e+04     2.718e+01
         200     1.200e+07     1.512e+04     7.992e+04     1.481e+03     1.409e+00
         400     6.980e+06     3.640e+03     1.946e+04     1.848e+02     7.913e-02
         800     3.763e+06     8.920e+02     4.799e+03     2.309e+01     4.667e-03
observed orders: euler 0.77, midpoint 2.06, heun 2.04, kutta3 3.00, rk4 4.17
ten orbits at 400 steps/orbit: final position error and relative energy / angular-momentum drift
  euler: position 2.303e+07 m, energy +4.876e-01, angular momentum +3.826e-01
   heun: position 2.441e+05 m, energy +1.407e-04, angular momentum +6.547e-05
    rk4: position 1.591e+00 m, energy -2.275e-09, angular momentum -1.027e-09
```

followed by the stiff-decay lines, the step-for-tolerance table and `check: RK4 with 210 steps/orbit gives 1.148 m`, and writes `rk4_orbit.png` (orbit error against steps per orbit for the five methods; energy drift over ten orbits; the stability regions of orders 1–4).

### Verify

Run `python -m pytest curriculum/phase-03/lesson-05-ode-integrators-euler-rk4`. The tests:

1. **Tableaus**: explicit, c = row sums of A, the Butcher conditions through each method's order; stage counts 1, 2, 2, 3, 4.
2. **Steps**: hand-written steps equal the tableau step to 10⁻¹⁴ on a nonlinear 4-state system; one step on y′ = y is the truncated exponential of the right degree for every method; RK4's local error on it is h⁵/120 within 15 %; RK4 is exact for quartic y(t) (random cubic g) and not quintic; Euler exact for linear y.
3. **integrate**: shapes, exact landing on t₁, h adjustment, the three `ValueError`s; exp(2) to 10⁻⁸ with RK4 at h = 0.01; Van der Pol (μ = 2) to 10⁻⁸ against DOP853.
4. **Orders**: local p + 1 within 0.1 on the oscillator (and within 0.35 on y′ = −2ty for the four lower methods: RK4's leading term changes sign there, a reminder that order estimates need a non-vanishing constant); global p within 0.12; the ranking Euler > Heun > Kutta3 > RK4 at 400 steps; `step_for_tolerance` reproduces 10⁻⁸ within 30 %.
5. **Stability**: R(z) for scalars and arrays; |R(−2)| = 1 for orders 1 and 2; the four limits to 10⁻³, with |R| > 1 just outside and < 1 just inside; RK4's imaginary limit between 2.8 and 2.85; order 5 raises; λ = −1000 decays at 0.95h_max and blows up at 1.05h_max for Euler, Heun and RK4; λ > 0 raises.
6. **Two-body**: the circular-orbit acceleration; perigee state and speed √(μ(1 + e)/(a(1 − e))); return to the start after one period; energy −μ/2a and angular momentum √(μa(1 − e²)) to 10⁻¹² at 20 random times; DOP853 agrees with `kepler_state` to 1 mm after half an orbit; Kepler's equation solved to 10⁻¹³.
7. **Orbit convergence**: orders 4, 3, 2 ± 0.3 for RK4, Kutta3, Heun over 200–800 steps; RK4 below 1 cm at 800 steps; Euler above 1000 km at 400.
8. **Drift**: RK4 energy and angular momentum within 10⁻⁸ over three orbits, the energy error negative and three times the one-orbit value within 25 %; Euler's energy rising monotonically.
9. **Cost**: the step predicted for 1 m from the 400-step error lies in 150–300 steps per orbit and delivers 0.5–2 m; Heun with twice the steps is worse by 10³.
10. **main** writes the figure.

### Compare

`scipy.integrate.solve_ivp` offers RK23, RK45 (Dormand–Prince), DOP853, and the implicit Radau, BDF and LSODA, all with adaptive steps and error control. Called with `method="DOP853", rtol=1e-13` it reproduces the Kepler orbit to 15 µm after one revolution and is the oracle for the Van der Pol test, which has no closed form. Nothing in this lesson is faster: `solve_ivp` picks its own step, lengthens it at apogee and shortens it at perigee, and reports an error estimate. What the hand-written methods bought: the knowledge that a fixed-step RK4 at 210 steps per orbit meets 1 m and at 400 meets 8 cm, and *why* (fourth order, and the measured constant); the stability limit 2.79/|λ| that explains a blow-up before you see it; and the energy-drift fingerprint, which the library's error estimate does not report and which is the first thing to plot after any long propagation. `P03.06` adds the step control; `P03.07` the conservation.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Convergence.** Tabulate the one-orbit position error of the five methods at 100, 200, …, 6400 steps per orbit and fit the orders. Predict 1, 2, 2, 3, 4; predict that Euler needs more than 3200 steps before its order appears (the orbit must first be resolved) and that RK4 reaches 10⁻⁴ m near 2000 steps.
2. **Work–precision.** Plot error against evaluations (stages × steps). Predict that for any tolerance below 1 km RK4 is cheapest, and that Heun beats Euler everywhere.
3. **Stability.** For y′ = −1000y, scan h from 0.9 to 1.1 h_max in 2 % steps for each method and mark the first blow-up. Predict the limits of 2.5 within one 2 % step.
4. **Error growth in time.** With RK4 at 800 steps per orbit, record the position error after 1, 2, 4, 8 orbits. Predict linear growth (the error is dominated by a secular phase error: the numerical orbit has a slightly wrong period) and compare with the energy drift, which should also be linear and negative.
5. **Eccentricity.** Repeat the convergence study for e = 0.3 and e = 0.6 at the same a. Predict that the RK4 error at fixed steps per orbit grows steeply with e (the perigee speed grows and the fixed step wastes work at apogee: roughly as the fourth power of the perigee-to-apogee speed ratio), motivating `P03.06`.

Raw data format: CSV tables `method, steps_per_orbit, n_evals, pos_err_m` (1–2), `method, h_over_hmax, final_abs_y` (3), `orbits, pos_err_m, energy_drift` (4), `e, steps_per_orbit, pos_err_m` (5); plus `rk4_orbit.png`.

What you should find: orders 0.99 (Euler, 3200–12 800), 2.0, 2.0, 3.0, 4.1; RK4 at 4.7 mm with 800 steps and 1.15 × 10⁻⁴ m at 2000; 1.6 m after ten orbits at 400 steps (twenty times the one-orbit 8 cm: the growth is linear once the error vector has settled into the along-track direction, and the first orbit is not yet there); RK4 cheapest for every tolerance below 10 km; blow-ups between 1.00 and 1.02 h_max; position error 4.7, 9.9, 22.0, 52.9 mm after 1, 2, 4, 8 orbits at 800 steps, i.e. linear growth of about 6 mm per orbit once settled; at e = 0.6 the RK4 error at 800 steps is 3.9 m instead of 4.7 mm (e = 0.3: 3.1 cm; e = 0: 4.3 mm).

## 5. Validate and improve

| quantity | theory (§2) | observed (code) | difference | cause |
|---|---|---|---|---|
| local orders | 2, 3, 3, 4, 5 | 2.00, 3.00, 3.00, 4.00, 5.00 | 0 | oscillator: constants do not vanish |
| global orders (oscillator) | 1, 2, 2, 3, 4 | 1.08, 2.00, 2.00, 3.00, 4.00 | 0.08 for Euler | 5 % error at 400 steps: still a visible h² term |
| RK4 local error on y′ = y, h = 0.5 | h⁵/120 = 2.60 × 10⁻⁴ | 2.84 × 10⁻⁴ | 9 % | h⁶/720 term |
| real-axis limits | 2, 2, 2.5127, 2.7853 | same to 10⁻⁴ | 0 | — |
| blow-up of y′ = −1000y | h_max = 0.002, 0.00279 | decays at 0.95h_max, grows at 1.05h_max | within 5 % | — |
| orbit orders (100–800) | 1, 2, 2, 3, 4 | 0.77, 2.06, 2.04, 3.00, 4.17 | Euler 0.23 | Euler not asymptotic: 10⁷ m errors |
| Euler order (3200–12 800) | 1 | 0.99 | 0.01 | asymptotic at last |
| RK4 step for 1 m from (2) | 209 steps | 210 steps give 1.15 m | 15 % | order 4.17 ≠ 4 over the range |
| Kepler reference vs DOP853 | 0 | 1.5 × 10⁻⁵ m | — | DOP853 at rtol 10⁻¹³ |
| RK4 energy drift, 10 orbits | linear, negative | −2.3 × 10⁻⁹, 10 × the one-orbit value | — | dissipative RK4 |
| Euler energy drift, 10 orbits | positive, growing | +49 % | — | Euler amplifies imaginary modes |

Everything with a theoretical value is met to the stated figure. The two rows that depart, Euler's orbit order and the 15 % miss of the step prediction, are the same lesson: an order is a property of the asymptotic regime, and before you are in it neither the order nor any extrapolation based on it means anything. Resolve the solution first (a few dozen steps per period of the fastest motion), measure the order, *then* extrapolate.

What to improve: never use Euler for anything but a first sanity check; prefer RK4 for smooth problems with mild eccentricity at fixed step; before any fixed-step run, estimate |λ|_max and check h against 2.79/|λ|. What not to improve here: the position error of a 400-step RK4 orbit (8 cm) against the model error of a point-mass Earth (J₂ moves the perigee by kilometres per day). The next errors to attack are adaptive step control for eccentric orbits (`P03.06`) and long-term energy behaviour (`P03.07`), not more digits from RK4.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of Euler's local and global error, of conditions (1), of the stability function and of RK4-as-Simpson; `rk4_orbit.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a method and a step for propagating a satellite in the lesson's orbit for one day with a 10 m position error, with the evaluation count and the evidence (measured error at one step, observed order, the growth law in time).

### Adversarial self-checks

- **What if the step doubles?** RK4's one-orbit error rises 16-fold (8 cm → 1.4 m at 200 steps), Euler's about doubles, and if the doubled step crosses 2.79/|λ| for any mode the solution blows up regardless of accuracy; the energy drift of RK4 rises 16-fold too.
- **What if the orbit is nearly circular (e = 0)?** Every method does better at fixed steps per orbit (the speed is constant and the fixed step is efficient); the order is unchanged. At e = 0.6 RK4's error at 800 steps is 3.9 m, a thousand times the e = 0.1 value: the fixed step is the wrong tool and adaptivity (`P03.06`) is the fix.
- **Which assumption fails first in the real world?** The point-mass Earth: J₂ alone perturbs the orbit by kilometres per day, 10⁴ times the RK4 error at 400 steps. Numerically, the first failure is a thrust event inside a step, which collapses the order to 1.
- **What measurement could falsify the integrator?** The energy history: a conservative system whose numerical energy drifts secularly is being integrated wrongly, whatever the position error looks like; the conserved quantity is the cheapest and most sensitive diagnostic available.
- **Does the result survive a change in method or step?** RK4 at 400 and 800 steps and DOP853 agree on the one-orbit state to 8 cm, and Kutta3 at 800 to 23 m; the spread between methods is the error estimate when there is no Kepler to compare with.
- **Which decision changes if the requirement tightens from 1 m to 1 mm per orbit?** Formula (2) gives RK4 1180 steps per orbit (4700 evaluations; a run at 1200 gives 0.90 mm), still cheap; Heun would need 5.5 × 10⁶; and at 1 mm the point-mass model error dominates by 10⁷, so the real decision is to add J₂ and the Moon, not to refine the step.

## Check yourself

1. State the local and global error orders of Euler, Heun and RK4, and explain in one sentence why the global order is one less than the local.
2. Write the two conditions a two-stage explicit RK method must satisfy to be second order, and give two methods that satisfy them.
3. What is the stability function of RK4, and what is its real-axis stability limit? What does the limit mean for y′ = −1000y?
4. Why is forward Euler unstable on an undamped oscillator at every step size, and what does that do to a numerical orbit?
5. You measure a 4.7 mm one-orbit error with RK4 at 800 steps per orbit. How many steps per orbit give 1 m, and what assumption does the answer rest on?
6. On y′ = g(t), which quadrature rule is RK4, and for which polynomial degree of g is it exact?
7. RK4's energy error on an orbit is negative and grows linearly; Euler's is positive. Explain both from the stability function.
8. A problem is stiff. In one sentence, what decides the step, and what is the remedy?
