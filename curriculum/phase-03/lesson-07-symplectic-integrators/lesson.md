# P03.07  Symplectic Integrators

> A fourth-order method gets each step more nearly right; a symplectic method gets the *kind* of motion right, so its errors stop adding up. For a million orbits the second matters more.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Störmer–Verlet (leapfrog, velocity Verlet) step from the centred difference of q'' = a(q), show it is second order, symmetric in time, area-preserving on the harmonic oscillator, and that it conserves a *modified* energy exactly. 2. State what a symplectic map is, test the Jacobian of a one-step map for symplecticity, and explain through the modified Hamiltonian why the energy error of a symplectic method stays bounded while that of Runge–Kutta drifts. 3. Implement and test, from scratch, explicit Euler, symplectic Euler, Störmer–Verlet, Yoshida's fourth-order composition and RK4 on the two-body problem, against the closed-form Kepler solution and `solve_ivp`, with observed orders 1, 2, 4, 4. 4. Decide, for a given propagation length, eccentricity and accuracy, whether a symplectic fixed-step method or an adaptive Runge–Kutta method is the right tool. |
| **Prerequisites** | `P03.05` (RK4, order of accuracy), `P03.06` (stability functions, convergence studies), `P03.02` (Newton's method for Kepler's equation), `P02.04` (phase-plane view of a second-order system), `P01.16` (determinants and Jacobians). Skill check: you can write the energy of a two-body orbit and compute its period from the semi-major axis. |
| **Engineering question** | A satellite in a 7 000 km, e = 0.1 orbit is to be propagated for hundreds of orbits with a fixed-step integrator at 200 steps per orbit: which methods keep the orbital energy and angular momentum, which let the orbit spiral, by how much per orbit, and at what cost in force evaluations? |
| **Assumptions and limits** | Two-body point-mass gravity, planar motion, no perturbations; fixed time step; the oscillator and Kepler problems as the test cases. Symplectic methods are stated for separable Hamiltonians H = p²/2 + V(q) with unit mass; drag, thrust and time-dependent forces are not Hamiltonian and the conservation arguments do not apply to them. Fixed steps are a poor fit for high eccentricity; the lesson quantifies why and leaves the remedy (time regularisation) to the extensions. |
| **Analytical work** | Verlet from the centred difference; its 2×2 map on the oscillator, its determinant, its numerical frequency and the exactly conserved modified energy; the composition rule and Yoshida's coefficients; conservation of angular momentum by kick-drift-kick maps; symmetry and reversibility; symplectic Euler as a conjugate of Verlet; the closed-form two-body solution via Kepler's equation. |
| **Computational deliverable** | `starter/main.py`: `euler_step`, `symplectic_euler_step`, `stormer_verlet_step`, `yoshida4_step`, `rk4_step`, `propagate`, `solve_kepler`, `kepler_exact`, `step_jacobian`, `symplectic_defect`, `reversibility_defect`, `modified_energy_oscillator`, `energy_drift`, `drift_slope`, `position_errors` (the force, invariants, `perigee_state` and `convergence_order` are given). Tests prove: the Kepler closed form matches DOP853 to 10⁻³ m and returns to perigee after one period; each step matches its hand formula on the oscillator; observed orders over 1.5 orbits are 1, 2, 4, 4 (and symplectic Euler looks second order over exactly one period); the symplectic maps have unit Jacobian determinant and zero symplectic defect to 10⁻⁷ while RK4 shrinks phase area by 2.4 % at hω = 1; Verlet and Yoshida are reversible to 10⁻¹²; Verlet conserves the modified energy to 10⁻¹² while the true energy swings by (hω)²/4; the kick-drift methods conserve angular momentum to 10⁻¹² and RK4 does not; over 40 orbits Verlet's energy error has the same amplitude at the end as at the start and no trend, RK4's grows linearly, and explicit Euler's grows monotonically. |
| **Practical activity** | **computational experiment**: a 60-orbit energy and angular-momentum study of six integrator/step combinations at matched cost, a convergence study over 1.5 orbits, and an eccentricity sweep that shows where fixed steps fail. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `symplectic_orbit.png`, the hand derivation of the Verlet map's determinant and modified energy, the 60-orbit table, and a half-page note recommending an integrator for a 10 000-orbit debris propagation with a cost estimate. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the modified energy of Verlet on the oscillator is derived; all public tests pass; the drift rates of RK4 and the bounded swing of Verlet are reproduced and explained; the note states a method, a step and a cost. |
| **References** | `@langtangen-linge-fd-computing-pdes` Chapter 1 "Vibration ODEs": "1.1.2 A Centered Finite Difference Scheme", "1.4 Analysis of the Numerical Scheme" (1.4.2 The Error in the Numerical Frequency, 1.4.7 Stability), "1.6 Energy Considerations" (1.6.1 Derivation of the Energy Expression, 1.6.2 An Error Measure Based on Energy), "1.7 The Euler-Cromer Method" (1.7.1 Forward-Backward Discretization, 1.7.2 Equivalence with the Scheme for the Second-Order ODE, 1.7.4 The Störmer-Verlet Algorithm) (required-free); `@langtangen-fd-computing-exponential-decay` Chapter 3 "Generalizations": "3.2.4 Leapfrog Schemes" (supplementary-free); `@ocw-18-085` "Lecture 10: Finite Differences in Time" (supplementary-free: leapfrog and its energy behaviour); `@langtangen-linge-fd-computing-pdes` "1.5 Alternative Schemes Based on 1st-Order Equations" (1.5.1 The Forward Euler Scheme, 1.5.5 Runge-Kutta Methods) and "1.9 Exercises and Problems" Problem 1.6 "Derive and investigate the velocity Verlet method" (supplementary-free). |
| **Extensions** | ★ Implement the Sundman time transformation dt = r dτ (steps in fictitious time τ, so that the physical step shrinks at perigee) and repeat the eccentricity sweep. ★ Add the J₂ perturbation (a conservative force, so still Hamiltonian) and show the energy stays bounded while the orbit precesses. ★ Add a drag term, show that no method conserves energy any more, and explain why a symplectic method is then no better than RK4 of the same order. |
| **Time** | unmeasured |

## 1. The engineering problem

Propagating an orbit is integrating q'' = −μ q/|q|³. `P03.05` and `P03.06` would reach for RK4 or an adaptive Dormand–Prince method, and for a few orbits they are excellent. But orbital mechanics often asks for *long* propagations: a debris cloud over decades, a constellation's drift over years, a planetary system over 10⁸ years, a Monte Carlo of 10⁴ trajectories each for a thousand orbits. Over such spans the question changes from "how accurate is each step?" to "what does the accumulated error *do*?" An RK4 integration of a two-body orbit does not merely acquire position error; its energy drifts in one direction, so the orbit slowly spirals in or out and its period changes. Position error after N orbits then grows as N², not N, and no quantity of the simulation can be trusted qualitatively. A **symplectic** integrator makes an error of a different kind: it solves *exactly* (to round-off) a slightly different Hamiltonian problem. Its energy error oscillates within a bound set by h² and never drifts; the orbit is a slightly wrong ellipse rather than a slowly decaying spiral.

This lesson builds the family: symplectic Euler, Störmer–Verlet (the leapfrog of `@ocw-18-085` Lecture 10 and `@langtangen-linge-fd-computing-pdes` 1.1.2, the velocity Verlet of molecular dynamics, the Euler–Cromer of 1.7, all the same scheme), and Yoshida's fourth-order composition, with explicit Euler and RK4 as the non-symplectic foils. It tests them against the closed-form Kepler solution and quantifies what "bounded" and "secular" mean in numbers. The test orbit is a = 7 000 km, e = 0.1 (perigee 6 300 km, period 97.1 min).

Success criteria, written first: the from-scratch Kepler solution agrees with DOP853 at rtol 10⁻¹³ to 10⁻³ m; the observed orders are 1, 2, 4 and 4; Verlet and Yoshida conserve angular momentum to round-off and show no energy trend over 60 orbits (slope below 10⁻⁹ per orbit) while RK4's energy error grows linearly at about 10⁻⁸ per orbit; and the one-step Jacobian of every symplectic method satisfies MᵀJM = J to 10⁻⁷.

## 2. Understand and derive

### 2.1 Störmer–Verlet from the centred difference

For q'' = a(q), the centred second difference (`@langtangen-linge-fd-computing-pdes` 1.1.2) gives

**q_{n+1} − 2q_n + q_{n−1} = h² a(q_n).**   (1)

It is second order (the odd Taylor terms cancel by symmetry) and symmetric: swapping n + 1 and n − 1 leaves it unchanged, so running it backwards retraces the same points. Introduce the velocity p_{n+½} = (q_{n+1} − q_n)/h and (1) becomes the **leapfrog**: p_{n+½} = p_{n−½} + h a(q_n), q_{n+1} = q_n + h p_{n+½}. Shifting the velocity to the whole steps by p_n = (p_{n−½} + p_{n+½})/2 gives the **velocity Verlet** form used in the code, **kick–drift–kick**:

**p_{n+½} = p_n + (h/2) a(q_n),  q_{n+1} = q_n + h p_{n+½},  p_{n+1} = p_{n+½} + (h/2) a(q_{n+1}).**   (2)

One force evaluation per step (the last kick reuses the next step's first), second order in both q and p, symmetric. The **symplectic Euler** (Euler–Cromer, 1.7.1) method is the half-way house: p_{n+1} = p_n + h a(q_n), q_{n+1} = q_n + h p_{n+1}; the velocity update is explicit and the position update uses the *new* velocity. It is first order, but 1.7.2 shows it is (1) in disguise: two symplectic Euler steps compose into a Verlet step shifted by half a step. **Explicit Euler** updates both with old values, q_{n+1} = q_n + h p_n, p_{n+1} = p_n + h a(q_n), and that one difference, old versus new velocity, decides everything that follows.

### 2.2 The oscillator: a 2×2 map and what it preserves

On q'' = −ω²q with unit mass the Verlet step (2) is linear, (q, p)_{n+1} = M (q, p)_n with

**M = [[1 − h²ω²/2, h], [−hω²(1 − h²ω²/4), 1 − h²ω²/2]],  det M = 1.**   (3)

(Check: (1 − x/2)² + x(1 − x/4) = 1 − x + x²/4 + x − x²/4 = 1 with x = h²ω².) det M = 1 means the map preserves area in the (q, p) plane; for explicit Euler det = 1 + h²ω² > 1 (spirals out, energy grows, `@langtangen-linge-fd-computing-pdes` 1.6 "Growth of energy in the Forward Euler scheme") and for RK4 det = |R(iωh)|² < 1 (spirals in; 0.9758 at hω = 1). The eigenvalues of M lie on the unit circle, e^{±iω_h h} with cos(ω_h h) = 1 − h²ω²/2, for hω < 2: the **numerical frequency** is ω_h = (2/h) asin(hω/2) = ω(1 + h²ω²/24 + …) (1.4.2), the method is stable for hω < 2 (1.4.7), that is for more than π steps per period, and above it the eigenvalues become real and the solution grows exponentially (the test checks hω = 2.1).

A 2×2 matrix with det M = 1 preserves a quadratic form: with M = [[α, β], [γ, α]], the form S = diag(−γ, β) satisfies MᵀSM = S exactly when α² − βγ = 1. Dividing by β = h, Verlet conserves **exactly**

**H̃(q, p) = p²/2 + (ω²/2)(1 − h²ω²/4) q²,**   (4)

a *modified* energy that differs from the true H = p²/2 + ω²q²/2 by (h²ω⁴/8) q². The true energy therefore oscillates as q² does, with relative amplitude h²ω²/4 (the test: 6.25 % at hω = 0.5, while H̃ is constant to 5 × 10⁻¹⁵), and never drifts, because H̃ is a constant of the numerical motion and the two differ by a bounded amount. This is the whole phenomenon in one formula, and it has nothing to do with accuracy: at hω = 0.5 the energy error is 6 %, but it is the same 6 % after a million steps.

### 2.3 Symplectic maps and the modified Hamiltonian

The exact flow of any Hamiltonian system (q, p) → (q(t), p(t)) is a **symplectic map**: its Jacobian M satisfies MᵀJM = J with J = [[0, I], [−I, 0]]. In one degree of freedom this is det M = 1; in more it preserves the sum of the projected areas and, in particular, phase-space volume. A numerical one-step method is called symplectic when its map has the same property for every h. Symplectic Euler and Verlet are symplectic because each of their sub-steps is: a kick p ← p + h a(q) has Jacobian [[I, 0], [h ∂a/∂q, I]], a drift q ← q + h p has [[I, hI], [0, I]], and both satisfy MᵀJM = J whenever ∂a/∂q is symmetric, which it is for any force derived from a potential (`P01.08`: mixed partials are equal). Compositions of symplectic maps are symplectic, so **Yoshida's** fourth-order method, three Verlet steps of sizes w₁h, w₀h, w₁h with w₁ = 1/(2 − 2^{1/3}) = 1.3512 and w₀ = 1 − 2w₁ = −1.7024 (a backward step in the middle, which cancels the h³ error term: the conditions are w₀ + 2w₁ = 1 and w₀³ + 2w₁³ = 0), is symplectic and fourth order at three force evaluations per step. RK4 is not: its map is a polynomial in the Jacobian, and no explicit Runge–Kutta method satisfies MᵀJM = J for general h.

Why does symplecticity bound the energy error? The result of (4) generalises: for a symplectic method of order p applied to H there is a **modified Hamiltonian** H̃ = H + h^p H_p + h^{p+1} H_{p+1} + … whose *exact* flow the numerical map follows (to all orders in h; the series is asymptotic, and the statement holds up to an exponentially small remainder for times of order e^{c/h}). The numerical trajectory therefore conserves H̃ and, since H − H̃ = O(h^p) is bounded along a bounded trajectory, the error in H oscillates within O(h^p) for exponentially long times. A non-symplectic method has no such H̃: its energy error accumulates like any other error, at a rate O(h^p) per unit time, and the drift is **secular**, proportional to t. Both behaviours are measured in part 4: Verlet's energy error at 200 steps per orbit swings by 1.06 × 10⁻⁴ every orbit and has a fitted slope of 10⁻⁹ per orbit (noise); RK4's at the same step drifts at −7.3 × 10⁻⁹ per orbit, 60 times smaller than Verlet's swing after 60 orbits but growing without limit.

Two more properties of the kick-drift-kick form matter for orbits. **Angular momentum**: a kick changes p by a vector parallel to q (a central force), which leaves q × p unchanged; a drift changes q by a vector parallel to p, which also leaves q × p unchanged; so symplectic Euler, Verlet and Yoshida conserve L = q × p *exactly*, to round-off, at any step size, and RK4 does not (2 × 10⁻⁷ relative after 60 orbits). **Reversibility**: Verlet is symmetric, so a step of +h followed by a step of −h returns to the starting point to round-off (3 × 10⁻¹⁸ in the test); Yoshida inherits this; symplectic Euler, explicit Euler and RK4 do not (defects of 10⁻³, 10⁻³ and 4 × 10⁻¹⁰ in the nondimensional Kepler test). Symmetric methods have only even powers of h in their error expansion, which is why Verlet is second order although each half-kick is first order.

### 2.4 A curiosity that is also a diagnostic

Symplectic Euler is Verlet conjugated by a half-step drift (1.7.2). Write a kick as K_h and a drift as D_h, so that SE_h = D_h ∘ K_h and the drift–kick–drift Verlet step is Ṽ_h = D_{h/2} ∘ K_h ∘ D_{h/2}. Then, exactly, **SE_h^N = D_{h/2} ∘ Ṽ_h^N ∘ D_{−h/2}** (insert D_{h/2}D_{−h/2} = I between every pair of steps): the first-order error of symplectic Euler is a *fixed* change of variables of size h/2 × p, not an accumulating one. Measure its position error after exactly one period, when V_h^N has returned near the starting point, and the O(h) parts cancel: the observed order is 2.00 (errors 46 484, 11 617, 2 904, 726 m at 200–1 600 steps per orbit). Measure it after 1.5 periods and the order is 0.94. The test checks both, and the moral is general: an order measured at one special time can lie, so measure at a generic one.

### 2.5 The two-body oracle

The closed-form solution (`P03.02` solved Kepler's equation; here it is used in anger) takes the state at t = 0, forms a = −μ/(2E), the angular momentum h_z, the eccentricity vector e = (v × h)/μ − q/|q| and the eccentric anomaly E₀ from |q| = a(1 − e cos E₀), then for each t solves M = E − e sin E for E with M = E₀ − e sin E₀ + nt and n = √(μ/a³), and writes the position in the perifocal frame, a(cos E − e, √(1 − e²) sin E), rotated by the argument of perigee. `kepler_exact` agrees with DOP853 at rtol 10⁻¹³ to 1.4 × 10⁻⁵ m in position and 1.5 × 10⁻⁸ m/s in velocity over 2.5 orbits, and returns to the initial point after one period to 9 × 10⁻¹⁰ m: an oracle good to a nanometre on a 7 000 km orbit, so that even the 1 mm errors of RK4 at 2 400 steps per orbit are measured, not guessed.

### 2.6 Where the methods stop being valid

- **Fixed step and eccentricity.** The perigee angular rate ω_p = v_p/r_p exceeds the apogee rate by ((1 + e)/(1 − e))²: 1.5 at e = 0.1 but 45 at e = 0.74 (a Molniya orbit). A step chosen for perigee wastes most of the orbit; one chosen for apogee misses perigee and the energy swing grows as (hω_p)². Problem 3 quantifies it; the Sundman transformation of extension 1 is the symplectic-friendly cure (adaptive step control, which changes h from the error estimate, breaks symplecticity).
- **Non-Hamiltonian forces.** Drag, thrust and any velocity-dependent or time-dependent force have no H̃ to conserve; the energy *should* change, and a symplectic method offers nothing beyond its order (extension 3).
- **Round-off.** The 10⁻¹⁴ angular-momentum conservation is round-off; over 10⁹ steps even round-off accumulates as a random walk, and compensated summation or extended precision becomes the next topic.
- **Accuracy is not conservation.** Verlet at 200 steps per orbit conserves energy to 10⁻⁴ and angular momentum to 10⁻¹⁴ and is 1 140 km off in position after 60 orbits (phase error: the numerical period is wrong by h²ω²/24 per 2.2, and phase errors grow linearly even for symplectic methods). The right ellipse, traversed at slightly the wrong speed.

## 3. Build the computational model

Open `starter/main.py`. `kepler_acceleration`, `energy`, `angular_momentum`, `orbital_period`, `perigee_state`, `convergence_order` and `main` are given.

### Step 1: five one-step maps

Each takes `(accel, q, p, h)` and returns `(q_new, p_new)`:

```python
def stormer_verlet_step(accel, q, p, h):
    p_half = p + 0.5 * h * accel(q)
    q_new = q + h * p_half
    return q_new, p_half + 0.5 * h * accel(q_new)

def yoshida4_step(accel, q, p, h):
    for w in (YOSHIDA_W1, YOSHIDA_W0, YOSHIDA_W1):
        q, p = stormer_verlet_step(accel, q, p, w * h)
    return q, p
```

`euler_step` evaluates the force once at the old q and updates both with old values; `symplectic_euler_step` kicks first and drifts with the new p; `rk4_step` concatenates (q, p) and applies the `P03.05` stages to y' = (p, a(q)). `propagate` loops and stores.

### Step 2: the oracle

`solve_kepler` is Newton on E − e sin E − M, started at M (or π for e ≥ 0.8), to 10⁻¹⁴. `kepler_exact` follows 2.5; it raises for E ≥ 0.

### Step 3: symplecticity and reversibility

`step_jacobian` perturbs each of the four phase-space components by ±ε max(1, |y|) and takes central differences of the one-step map; on the oscillator the map is linear and the result is exact. `symplectic_defect` returns ‖MᵀJM − J‖_F. `reversibility_defect` steps +h then −h. `modified_energy_oscillator` is (4).

### Step 4: studies

`energy_drift` propagates for a number of orbits at a given number of steps per orbit and returns (t/P, (E − E₀)/|E₀|); `drift_slope` fits a line; `position_errors` compares q(t_end) with `kepler_exact` at each step count.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-07-symplectic-integrators`. The tests:

1. **Bookkeeping.** The perigee state has E = −μ/2a, L = √(μa(1 − e²)) and period 2π√(a³/μ); the force points inward with magnitude μ/r²; unbound and e ≥ 1 raise.
2. **Oracle.** Kepler's equation solved to 10⁻¹² for e up to 0.9; `kepler_exact` against DOP853 to 10⁻³ m and 10⁻⁶ m/s at five times, back at perigee after P to 10⁻⁶ m, at apogee a(1 + e) after P/2, a circular orbit keeps its radius.
3. **One step** of each method on the oscillator matches its hand formula; RK4 reproduces the degree-4 Taylor polynomial of the rotation; Yoshida's weights satisfy both order conditions.
4. **Orders** over 1.5 orbits: 1 ± 0.2, 2 ± 0.15, 4 ± 0.3, 4 ± 0.3 at 300, 600, 1 200 steps; symplectic Euler over exactly one orbit: 2 ± 0.15.
5. **Symplecticity** on the oscillator at hω = 1: det M = 1 to 10⁻⁸ and defect < 10⁻⁷ for the three symplectic maps; RK4's det < 0.99, Euler's > 1.01; on the nondimensional Kepler problem at h = 0.2 the defects are < 10⁻⁷ against > 10⁻⁷ for RK4.
6. **Reversibility**: Verlet and Yoshida < 10⁻¹², the others > 10⁻¹¹.
7. **Modified energy** constant to 10⁻¹² while the true energy swings by (hω)²/4 at hω = 0.5; instability beyond hω = 2.
8. **Angular momentum** to 10⁻¹² for the kick-drift methods over 10 orbits; RK4 above 10⁻⁹.
9. **Drift**: over 40 orbits Verlet's swing is below 2 × 10⁻⁴, its last-orbit amplitude equals its first to 5 % and its slope times 40 orbits is below 2 % of the swing; RK4's final error equals slope × 40 to 10 % and exceeds ten times its first-orbit swing; Euler's error exceeds 0.1 in three orbits and never decreases.
10. `main` writes the figure.

### Compare

`solve_ivp` offers no symplectic method; `DOP853` at rtol 10⁻¹³ is the oracle for the oracle. The industry tools for long propagations do: planetary ephemerides and asteroid-family studies use symplectic mixed-variable methods (Wisdom–Holman), molecular dynamics uses velocity Verlet almost exclusively, and particle accelerator codes use symplectic maps by construction, all because a bounded energy error over 10⁹ steps is worth more than a smaller error over 10³. For a few orbits, or for any problem with drag or thrust, an adaptive Dormand–Prince method remains the right default, and the 60-orbit table of part 4 is honest about it: at matched cost (800 force evaluations per orbit) RK4's energy error after 60 orbits, 4.4 × 10⁻⁷, is 15 times smaller than Verlet's swing of 6.6 × 10⁻⁶; the crossover, where RK4's linear drift at 7.3 × 10⁻⁹ per orbit overtakes Verlet's bound, is at about 900 orbits, and Yoshida at 600 evaluations per orbit beats both with 2.3 × 10⁻⁷ bounded. The decision rule, in numbers, is in Problem 5.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Sixty orbits.** For the e = 0.1 orbit, run explicit Euler (5 orbits only), symplectic Euler, Verlet at 200 and at 800 steps per orbit, RK4 and Yoshida at 200. Record the maximum |ΔE/E|, the final ΔE/E, the fitted slope per orbit, the maximum |ΔL/L|, the final position error against `kepler_exact`, and the force evaluations per orbit. Predict: Verlet's swing falls by 16 when the step is quartered; RK4's slope is negative (R(iωh) has modulus below 1); symplectic Euler's swing is an order of magnitude or more larger than Verlet's (first order against second: the ratio is of order 1/(hω_p) ≈ 26).
2. **Orders.** Position error after 1.5 orbits at 300, 600, 1 200, 2 400 steps; then symplectic Euler after exactly 1.0 orbits. Predict 1, 2, 4, 4 and the one-period surprise.
3. **Eccentricity sweep.** Keep the perigee at 6 300 km and 200 steps per orbit; vary e from 0 to 0.5 and record Verlet's energy swing. Predict growth with (hω_p)² = (2π v_p/(200 r_p n))² and, since the perigee passage also sharpens, somewhat faster.
4. **Matched cost.** Plot (E − E₀)/|E₀| for Verlet at 800 and RK4 at 200 steps per orbit on the same axes and extrapolate RK4's line to the orbit count where it crosses Verlet's bound.

Raw data format: CSV tables (`method, steps_per_orbit, max_dE, final_dE, slope_per_orbit, max_dL, pos_err_m, evals_per_orbit`; `method, n_steps, err_m`; `e, h_omega_p, max_dE`), each with a header naming the fixed parameters.

What you should find (60 orbits, e = 0.1): explicit Euler gains 48 % in energy in five orbits and spirals out (position error 19 900 km); symplectic Euler swings by 4.1 × 10⁻³, slope −1.8 × 10⁻⁶ per orbit (its swing is centred away from zero, so a short fit shows a spurious slope: look at the amplitude, which is constant), L to 4 × 10⁻¹⁵; Verlet at 200 swings by 1.06 × 10⁻⁴, final +3 × 10⁻⁷, slope +1.2 × 10⁻⁹ per orbit, L to 10⁻¹⁴, position error 1 140 km; Verlet at 800 swings by 6.59 × 10⁻⁶ (ratio 16.1), position error 71 km; RK4 at 200: maximum 4.4 × 10⁻⁷ equal to its final value, slope −7.28 × 10⁻⁹ per orbit, L to 2 × 10⁻⁷, position error 1.0 km; Yoshida at 200: swing 2.28 × 10⁻⁷, final 7 × 10⁻¹⁵, slope 10⁻¹⁴ per orbit, L to 10⁻¹⁴, position error 2.5 km. Orders 0.94, 2.00, 4.07, 4.00 and 2.00 for the one-period symplectic Euler. In the sweep, Verlet's swing at e = 0.25 is 3.7 × 10⁻⁴, 3.5 times the e = 0.1 value against 2.0 predicted from (hω_p)² alone.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| det M, Verlet, hω = 1 | 1 | 1.000000000292 | 3 × 10⁻¹⁰ | finite-difference Jacobian |
| det M, RK4, hω = 1 | \|R(i)\|² = 0.9758 | 0.975842 | 0 | exact |
| Verlet energy swing, hω = 0.5 | (hω)²/4 = 0.0625 | 0.06250 | < 10⁻⁵ | sampling of the extremes |
| modified energy H̃ | constant | varies by 5 × 10⁻¹⁵ | — | round-off |
| orders (1.5 orbits) | 1, 2, 4, 4 | 0.94, 2.00, 4.07, 4.00 | ≤ 0.07 | pre-asymptotic at 300 steps |
| Kepler closed form vs DOP853 (m) | 0 | 1.4 × 10⁻⁵ | — | rtol 10⁻¹³ on a 7 × 10⁶ m orbit |
| Verlet swing, 200 → 800 steps | ×1/16 | ×1/16.1 | 0.6 % | h² law |
| angular momentum, kick-drift methods | exact | 10⁻¹⁴ | — | round-off |
| RK4 energy after 60 orbits | slope × 60 = −4.37 × 10⁻⁷ | −4.37 × 10⁻⁷ | 0 | linear drift |
| Verlet trend over 60 orbits | 0 | 1.2 × 10⁻⁹ per orbit | — | fit noise on a 10⁻⁴ oscillation |

The integrators are validated against each other, against a closed form good to 10⁻⁵ m and against the hand formulas of 2.2 to round-off, so the numerical claims are settled. The engineering claim, "use a symplectic method for long propagations", is validated *in the regime tested*: a conservative two-body force, moderate eccentricity, 60 to a few hundred orbits, where the table shows exactly where RK4 is still better and where it stops being. Beyond that regime the lesson has not measured anything, and the honest extension is to measure: extension 2 (J₂) checks that a perturbed but conservative problem behaves the same way; extension 3 (drag) checks that a non-conservative one does not, and that the "advantage" vanishes.

Improvements worth making: force reuse in Verlet (the second kick's force is the next step's first; halves the cost as written); the Sundman transformation for eccentric orbits (extension 1); Yoshida's sixth- and eighth-order compositions for ephemeris-grade work. Improvement worth *not* making: adaptive step control on a symplectic method; it destroys the modified-Hamiltonian structure and the energy starts to drift like RK4's.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (3) and (4), of angular-momentum conservation by kick-drift maps, and of Yoshida's two order conditions; `symplectic_orbit.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending an integrator and a step for a 10 000-orbit propagation of a debris object in the e = 0.1 orbit with a 1 km position requirement, with a cost in force evaluations and a statement of what the requirement cannot be met by.

### Adversarial self-checks

- **What if the step doubles?** Verlet's energy swing quadruples (h²), Yoshida's and RK4's error grow 16-fold, RK4's drift rate 16-fold too; symplectic Euler's swing doubles. Beyond hω_p = 2 at perigee Verlet is unstable; for the e = 0.1 orbit that is about 4 steps per orbit, so 200 is far from the cliff.
- **What if the eccentricity goes to 0.74?** ω_p rises 45-fold relative to ω_a; at 200 steps per orbit hω_p = 0.31 and Verlet's swing is of order 5 % (it is 4.7 × 10⁻⁴ at 2 000 steps per orbit), while RK4 (fourth order in hω_p) is still fine: fixed-step symplectic integration *loses* to RK4 on a Molniya orbit unless the step is regularised.
- **Which assumption fails first in the real world?** Conservative forces: atmospheric drag below 800 km makes the energy change physically, so energy conservation is no longer the test of the integrator, and the first thing to check becomes the drag model, not the integrator.
- **What measurement could falsify the model?** Tracking data: a real orbit's semi-major axis decays by drag at a rate the two-body propagation cannot show; an integrator-induced decay (RK4's) looks the same in the output, which is exactly why the integrator must be shown drift-free *before* the drag model is tuned to the data.
- **Does the result survive a change of solver?** Verlet, Yoshida and DOP853 agree on the orbit to the accuracy of each; the *energy behaviour* does not survive a change of solver, by design: that is the result.
- **Which decision changes if the requirement doubles?** With a 2 km instead of 1 km position tolerance after 10 000 orbits, phase error dominates, and it grows linearly with time for the symplectic methods. Extrapolating the 60-orbit table (19 km per orbit for Verlet at 200 steps, 42 m per orbit for Yoshida, both ∝ h^p): Verlet would need about 60 000 steps per orbit, Yoshida about 800; for 1 km, 90 000 and 950. RK4 at 200 steps has 17 m per orbit *today*, but its energy drift of 7.3 × 10⁻⁹ per orbit changes the period progressively, so its phase error grows as N², and after 10 000 orbits the fractional energy error of 7 × 10⁻⁵ alone shifts the semi-major axis by 500 m. The symplectic *property* buys the qualitative trustworthiness of the 10 000-orbit answer; the *order* buys its accuracy, and the requirement is met by order.

## Check yourself

1. Write the velocity Verlet step and say which quantity is evaluated at the half step.
2. What is det M for the Verlet map on the oscillator, and what does it say about areas in the (q, p) plane?
3. State the modified energy that Verlet conserves on the oscillator and the size of the true energy's oscillation.
4. Why does the energy error of RK4 drift while that of Verlet stays bounded? Name the object that makes the difference.
5. Why do kick-drift-kick methods conserve angular momentum exactly for a central force?
6. What are Yoshida's three sub-step weights, and why must one of them be negative?
7. Symplectic Euler is first order. Why does its position error after exactly one period scale as h²?
8. Give one orbital situation where a symplectic method is the wrong choice and say why.
