# P03.10  Nonlinear Systems

> Trim is three equations in three unknowns that are not linear: the drag depends on the square of the lift, the thrust tilts with the angle of attack. Newton's method linearises, solves, and repeats; close to the answer the number of correct digits doubles every step. Far from it, it can fly off to infinity, and the lesson is as much about that as about the doubling.

| | |
|---|---|
| **Learning outcomes** | 1. Write a nonlinear engineering problem as F(x) = 0 with scaled unknowns and residuals, derive the Newton–Raphson iteration from a Taylor expansion, and prove its quadratic convergence and its conditions. 2. Build the Jacobian three ways (analytic, forward differences, complex step) and state the error of each; derive the Armijo line search and Broyden's secant update. 3. Implement and test, from scratch, Newton, damped Newton, Broyden, continuation over a parameter, and Picard iteration, against `scipy.optimize.root` and `newton_krylov` as oracles, measuring the observed orders 2, 1 and the superlinear case. 4. Diagnose a failing Newton iteration (bad start, singular Jacobian, no solution) and choose among Newton, quasi-Newton, Picard and Newton–Krylov for a given problem size and cost of F. |
| **Prerequisites** | `P03.02` (Newton in one variable, convergence orders), `P03.08` (LU, condition number), `P03.09` (sparse SPD systems, CG), `P01.14` (Jacobian matrices, Taylor expansion in several variables), `P02.11` (lift, drag polar, pitching moment). Skill check: you can write Newton's iteration for a scalar f(x) = 0 and say why it converges quadratically. |
| **Engineering question** | A 40 t regional jet (S = 80 m², c = 3 m) cruises at 200 m/s at 5 km: what angle of attack, elevator deflection and thrust hold it in steady level flight, how do they change from 110 to 260 m/s, and how many iterations and residual evaluations does each answer cost with Newton, with a differenced Jacobian, with Broyden, and with the previous speed's answer as the starting point? |
| **Assumptions and limits** | Steady level flight (flight-path angle zero, so the angle of attack is the pitch attitude), linear lift and moment coefficients in α and δe, parabolic drag polar, thrust along the body axis offset z_T below the centre of gravity, ISA troposphere, no compressibility and no stall model (the linear C_L is used up to 9°, which is at the edge of validity). The conduction plate has a linear k(T), which goes negative at high temperature and then has no solution; that failure is part of the lesson. Dense Jacobians throughout (n ≤ 256). |
| **Analytical work** | Taylor expansion and the Newton step; the quadratic-convergence theorem with its constant ½‖J⁻¹‖‖F''‖; the forward-difference error h F''/2 + 2ε/h and the optimal h ≈ √ε; the complex-step formula; the Newton direction as a descent direction for ‖F‖²/2 and the Armijo condition; Broyden's least-change secant update; the divergence threshold of Newton on arctan x; hand trim of the jet from C_L = W/(qS). |
| **Computational deliverable** | `starter/main.py`: `trim_residual`, `trim_jacobian`, `fd_jacobian`, `complex_step_jacobian`, `newton`, `damped_newton`, `broyden`, `solve_trim`, `trim_sweep`, `conduction_matrix`, `conduction_boundary_rhs`, `plate_residual`, `picard_plate`, `newton_plate` (the aircraft data, the atmosphere, the initial guesses, `iterate_errors`, `convergence_order`, `stall_speed`, the arctan test function, `face_conductances`, `HotPlate`, `plate_kfield` and `main` are given). Tests prove: the residual vanishes on a hand-built state; the analytic Jacobian matches forward differences to 10⁻⁶ and the complex step to 10⁻¹³; the forward-difference error is first order in h (slope 1.0 ± 0.1); Newton on x² − 2 is the Babylonian iteration with observed order 2 ± 0.2; the differenced Newton counts n + 1 evaluations per step; undamped Newton diverges on arctan x from 1.5 and 3 while the damped version converges with one shortened step then full steps; the order estimator returns 2 and 1 on synthetic sequences; Newton trim from (0, 0, 0) converges in ≤ 6 steps to the `scipy.optimize.root` solution with order 1.6–2.7 and C_L = W/(qS) to 0.5 %; differenced Newton, Broyden and damped Newton reach the same trim, Broyden with fewer evaluations and more steps, with improving contraction; the physical initial guess converges no slower than the crude one; the warm-started sweep equals the cold-started one to 10⁻⁸ with no more iterations, α falls and δe rises with speed, and the minimum thrust lies within 10 m/s of the minimum-drag speed; the conduction matrix equals the Kronecker Laplacian for constant k, is SPD, and gives the wall temperature without a source; the linear plate converges in one Picard and one Newton step; on the nonlinear plate Picard has order 1 ± 0.15 and Newton above 1.6 in ≤ 6 steps, both agree to 10⁻⁹ and with `newton_krylov` to 10⁻⁸, and the plate runs hotter than its linear-k model. |
| **Practical activity** | **computational experiment**: a Jacobian-accuracy study against the step h, a convergence race (Newton analytic, Newton differenced, Broyden, damped Newton, `scipy.optimize.root`) from a crude start, a trim sweep with and without continuation, the arctan divergence demonstration, and a Picard-against-Newton study on a 256-unknown plate with increasing nonlinearity up to the point where no solution exists. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `nonlinear_systems.png`, the hand derivations of the Newton step, the quadratic-convergence constant, the optimal difference step and the Armijo condition, the four experiment tables, and a half-page note stating the trim at 200 m/s with the evaluation cost of each method and the start you would use in a sweep. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the Newton step and its convergence are derived; all public tests pass; the observed orders are reported and explained; the note names the cause of each failure mode seen and the remedy. |
| **References** | `@ocw-16-90` Unit 1 "1.7.4 Newton-Raphson Implement Implicit Methods on Nonlinear Problems", "1.7.5 Apply Newton-Rhapson" (required-free: Newton for systems arising from implicit time stepping, titles as published); `@ocw-18-335j` Week 11 "Lecture 30: Quasi-Newton Methods and the BFGS Algorithm" (required-free: secant updates, the least-change idea behind Broyden); `@scipy-userguide` "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations" (required-free: `root`, the `hybr` and `lm` methods), "Root finding for large problems", "Still too slow? Preconditioning." (supplementary-free: `newton_krylov` and preconditioned Newton–Krylov for the plate-sized problems of `P03.09` and `P03.11`). |
| **Extensions** | ★ Write a Jacobian-free Newton–Krylov solver: Newton with the inner linear system solved by your `P03.09` CG or a GMRES you write, with J·v approximated by (F(x + εv) − F(x))/ε; compare with `newton_krylov` on the plate at n = 64. ★ Add a stall model C_L = C_Lmax sin(2α)/... or a tabulated lift curve to the trim equations and show where Newton needs damping. ★ Implement pseudo-arclength continuation and trace the trim curve through the minimum-thrust point as a function of thrust rather than speed, where there are two solutions. |
| **Time** | unmeasured |

## 1. The engineering problem

Every lesson so far ended in a linear system. Trim does not: the lift balance qS C_L(α, δe) + T sin α = W, the drag balance T cos α = qS [C_D0 + k C_L²] and the moment balance qS c C_m(α, δe) + T z_T = 0 are three equations in three unknowns with a square, a sine and a cosine in them. The implicit integrators of `P03.06` solved exactly such a system at every step (`@ocw-16-90`, 1.7.4), and the finite-element and optimisation lessons ahead will too. The method is always Newton–Raphson: linearise about the current guess, solve the linear system, correct, repeat. This lesson writes it for systems, measures the quadratic convergence that makes it the default, and then spends as much effort on the three ways it fails: a bad start, a Jacobian that is wrong or singular, and a problem with no solution.

The vehicle is a regional jet, W = 392.3 kN, cruising at 200 m/s at 5 km (ρ = 0.7361 kg/m³, qS = 1 177.8 kN). The hand answer is immediate for the lift: C_L = W/(qS) = 0.3331, α ≈ 0.86° with the elevator neutral; the thrust is drag, qS(C_D0 + k C_L²) = 35.32 kN; the moment equation then wants an elevator deflection that changes C_L, and so on round. Newton does all three at once: **α = 0.700°, δe = 2.218°, T = 35.31 kN**, in 2 iterations from the hand guess and 3 from (0, 0, 0), with residual norms 3.0 × 10⁻¹, 1.4 × 10⁻³, 7.5 × 10⁻⁸, 1.8 × 10⁻¹⁶: four, then eight, then sixteen digits. The second vehicle is the plate of `P03.09` with a conductance that falls with temperature, 256 unknowns, where Newton meets Picard iteration and, when the heating is raised, meets a problem that has no answer.

Success criteria, written first: the Jacobian from code agrees with the analytic one to 10⁻⁶ (differences) and 10⁻¹³ (complex step); Newton from a crude start reaches the `scipy.optimize.root` solution to 10⁻¹⁰ in at most 6 iterations with observed order near 2; Broyden reaches it with fewer residual evaluations than differenced Newton; the trim sweep with warm starts needs no more iterations than with cold starts and gives the same answers to 10⁻⁸; undamped Newton on arctan x from 1.5 diverges and damped Newton converges; Picard on the plate shows order 1.0 ± 0.15 and Newton above 1.6, agreeing to 10⁻⁹.

## 2. Understand and derive

### 2.1 The Newton step and its convergence

Expand F about the current iterate: F(x + δ) = F(x) + J(x) δ + O(‖δ‖²), where J_ij = ∂F_i/∂x_j is the Jacobian. Setting the linear part to zero gives the **Newton step**

**J(x_k) δ_k = −F(x_k),  x_{k+1} = x_k + δ_k.**   (1)

Each iteration costs one Jacobian, one residual and one linear solve (`P03.08`: 2n³/3 flops, trivial here at n = 3, the dominant cost at n = 256). Subtract the root x* and expand F(x*) = 0 about x_k to second order: 0 = F_k + J_k (x* − x_k) + ½ F''_k [e_k, e_k] + …, with e_k = x_k − x*. Combine with (1):

**e_{k+1} = ½ J_k⁻¹ F''_k [e_k, e_k] + O(‖e_k‖³),  ‖e_{k+1}‖ ≤ C ‖e_k‖²,  C = ½ ‖J⁻¹‖ max ‖F''‖.**   (2)

Three conditions carry this result: J(x*) nonsingular, F twice differentiable near x*, and ‖e₀‖ small enough that C‖e₀‖ < 1. Then the error squares every step: 10⁻², 10⁻⁴, 10⁻⁸, 10⁻¹⁶. The code measures this with `convergence_order`, the slope of log e_{k+1} against log e_k: 2.30 for the trim from (0, 0, 0) (errors 9.9 × 10⁻², 9.3 × 10⁻⁴, 2.0 × 10⁻⁸; only two usable pairs, so the estimate is rough), 1.93 on the plate, 2.00 on x² − 2 (the test). Near a root where J is singular the convergence drops to linear; when F has a double root, Newton halves the error per step like bisection.

For the trim, the unknowns are scaled, x = (α [rad], δe [rad], T/W), and the residuals divided by W (and by c for the moment), so every entry of x and F is O(0.1–1) and the Jacobian is well conditioned, κ(J) = 17.3. Without the scaling, T in newtons and α in radians give κ(J) ≈ 10⁶ and the tolerance tests become meaningless; **scale first** is the rule whenever unknowns differ by orders of magnitude.

### 2.2 Three Jacobians

**Analytic.** Differentiate by hand: for the trim, ∂F₁/∂α = (qS C_Lα + T cos α)/W, ∂F₂/∂α = (−T sin α − 2 qS k C_L C_Lα)/W, ∂F₃/∂(T/W) = z_T/c, and so on; nine entries. It is exact, it is the fastest, and it is where the bugs are: the code checks it against the two numerical Jacobians at every development step, and so should you.

**Forward differences.** Column j is [F(x + h e_j) − F(x)]/h: n extra residual evaluations, with the error h‖F''‖/2 from truncation plus 2ε‖F‖/h from cancellation. The sum is minimised at **h ≈ √ε max(1, |x_j|) ≈ 10⁻⁸**, with error about √ε ≈ 10⁻⁸. The run shows it: errors 4.1 × 10⁻², 4.1 × 10⁻⁴, 4.1 × 10⁻⁶, 4.0 × 10⁻⁸ at h = 10⁻², 10⁻⁴, 10⁻⁶, 10⁻⁸ (first order, slope 1.0, as the test confirms), then back up to 1.9 × 10⁻⁶ at 10⁻¹⁰ and 1.2 × 10⁻⁴ at 10⁻¹²: the V-shape of `P03.04`. A Jacobian wrong by 10⁻⁸ barely slows Newton (the convergence becomes linear with factor 10⁻⁸, indistinguishable from quadratic until the error reaches 10⁻⁸), so differenced Jacobians are the practical default when F is cheap; the cost is n + 1 evaluations per step, 13 for the trim, 1 029 for the plate.

**Complex step.** If F is analytic and written with functions that accept complex arguments, F(x + ih e_j) = F(x) + ih J e_j − h² F''/2 − …, so **J e_j = Im F(x + ih e_j)/h** with no subtraction and no cancellation: h = 10⁻³⁰ gives the derivative to round-off, 8.9 × 10⁻¹⁶ in the run. The price is a complex evaluation (about four times a real one) and code that does not branch on comparisons or use `abs`. The tests use it as the referee between the analytic and differenced Jacobians.

### 2.3 When Newton fails, and damping

Condition three of (2) is the one that bites. For F(x) = arctan x, with the single root 0, the Newton step is x − (1 + x²) arctan x. From x₀ = 1 the iterates go 1, −0.570, 0.117, −0.001, 10⁻⁹: converged. From 1.5 they go 1.5, −1.69, 2.32, −5.11, 32.3, …, 2 × 10¹³: the tangent at a point where the function has flattened crosses zero far beyond the root, and each step overshoots more than the last. The threshold, where the step lands at −x₀ exactly, solves 2x = (1 + x²) arctan x: **x₀ = 1.3917** (Problem 2). Any Newton code used on a problem it has not seen before will meet this.

The remedy is to treat the Newton direction as a direction, not a step. The merit function m(x) = ½‖F(x)‖² has gradient ∇m = JᵀF, and along δ = −J⁻¹F the directional derivative is ∇m·δ = −FᵀJJ⁻¹F = −‖F‖² < 0: **the Newton direction always descends m**, whatever the start. So take x + λδ and shorten λ from 1 until the **Armijo condition**

**m(x + λδ) ≤ (1 − 2c₁λ) m(x),  c₁ = 10⁻⁴,**   (3)

holds (halving each time): a step that reduces ‖F‖ by at least a fraction of what the linear model promised. Near the root, λ = 1 satisfies (3) and the quadratic convergence of (2) is untouched; far away, the step is shortened and m decreases monotonically. In the run, damped Newton from 1.5 takes λ = 0.5 once, then 1, 1, 1 and converges in 4 iterations; from 3.0 it takes λ = 0.25 then full steps. On the trim from (0, 0, 0) every λ is 1: the crude start is already inside Newton's basin, so damping cost one extra residual evaluation (4 against 3) and changed nothing else. Damping cannot create a solution: for the overheated plate of part 4, where k(T) would have to go negative, the line search shrinks λ to 10⁻⁶ and the residual stalls at 3.7 × 10⁴, which is the correct report that no solution exists.

### 2.4 Quasi-Newton: Broyden's method

When the Jacobian is expensive and there is no analytic form, replace it by an approximation B that is corrected from the residual changes already computed. After a step s = x_{k+1} − x_k with residual change y = F_{k+1} − F_k, the true Jacobian satisfies approximately J s ≈ y (the **secant condition**). Broyden's "good" update chooses the B_{k+1} that satisfies B_{k+1} s = y exactly and is otherwise as close as possible to B_k in the Frobenius norm (`@ocw-18-335j` Lecture 30; Problem 2 proves the least-change property):

**B_{k+1} = B_k + (y − B_k s) sᵀ/(sᵀs).**   (4)

A rank-one change, n² flops, and one residual evaluation per iteration instead of n + 1. The price is **superlinear rather than quadratic** convergence: the trim from (0, 0, 0) took 5 Broyden steps and 9 evaluations (errors 9.9 × 10⁻², 9.3 × 10⁻⁴, 1.7 × 10⁻⁵, 1.0 × 10⁻⁷, 3.3 × 10⁻¹⁰, 3.7 × 10⁻¹⁴) against 3 Newton steps and 13 evaluations with a differenced Jacobian. The error ratios fall from 10⁻² to 10⁻⁴, the signature of superlinear convergence, but the least-squares order over so few points reads 1.09; the test checks the ratio trend, not a number. BFGS, the symmetric version for optimisation, is the subject of `P03.15`. `scipy.optimize.root` with `method="hybr"` (Powell's hybrid, the MINPACK code) solved the trim in 11 evaluations: it is a trust-region Newton with a Broyden-updated Jacobian, both ideas of this lesson in one routine (`@scipy-userguide`, "Sets of equations").

### 2.5 Continuation: the previous answer as the next start

A sweep over a parameter, speed here, is a sequence of related problems; the trim at 115 m/s is a few per cent away from the trim at 110 m/s. **Continuation** feeds each solution to the next problem as its initial guess. On the 31-point sweep from 110 to 260 m/s the warm start took 94 Newton iterations against 101 from the crude (0, 0, 0) start, the same answers to 8.6 × 10⁻¹³. The saving is modest here because the problem is nearly linear and the crude start is already inside the basin; it becomes decisive when the start matters (stall models, transonic drag rise, structural buckling) and when the curve folds back on itself (two trims for one thrust, Extension 3), where a cold start lands on the wrong branch or nowhere. The sweep itself is aircraft performance in miniature: α from 9.06° at 110 m/s (C_L ≈ 1.1) to −0.80° at 260 m/s, δe from −4.41° (nose-up) to +3.54°, and a thrust minimum of **26.27 kN at 135 m/s**, against the hand value of the minimum-drag speed, where C_D0 = kC_L², of 134 m/s. The stall speed at C_Lmax = 1.6 is 91.2 m/s, so the sweep stays within the linear lift range, barely.

### 2.6 Picard against Newton on the plate

The plate of `P03.09` with k(T) = k₀(1 + β(T − T_wall)) gives F(T) = A(T) T − b(T) = 0 with A the five-point matrix assembled from the current temperatures, 256 unknowns at n = 16. **Picard** (successive substitution) freezes the coefficient: solve A(T_k) T_{k+1} = b(T_k), a linear SPD system per iteration, no Jacobian. Its error contracts by a factor proportional to |β| ΔT, so it converges **linearly**, order 1: 10 iterations for 10⁻⁸ at β = −0.004 K⁻¹, with error ratios 0.062, 0.050, 0.042, 0.037 and measured order 1.03. **Newton** on the same F, with a differenced Jacobian (257 evaluations per step), takes 4 iterations with residuals 8.0 × 10⁴, 7.1 × 10³, 91, 1.5 × 10⁻², 9 × 10⁻⁹, order 1.93; the first Newton step *is* the first Picard step (same linear system) and the two then part ways. Both give a maximum rise of 47.48 K where the constant-k₀ model gives 42.97 K: the plate is 10 % hotter because its conductance has fallen to 0.275 W/K at the hot spot. As |β| grows the Picard count grows (9, 10, 12, 13 at β = −0.002 … −0.0075) while Newton stays at 4–5; at q = 10 kW/m² and β = −0.0075 K⁻¹ the linear k(T) would have to fall to 0.03 k₀, and both methods fail because there is nothing to find, which is the model's fault, not the solver's. `scipy.optimize.newton_krylov` agrees with our Newton to 7 × 10⁻¹³ K in 30 ms, against 220 ms for our dense differenced Jacobian: it never forms J, approximating J·v by a difference of residuals and solving the Newton system by GMRES (`@scipy-userguide`, "Root finding for large problems"); at n = 1 000 that is the only way, and the preconditioner of `P03.09` is what makes it fast ("Still too slow? Preconditioning.").

### 2.7 Where the methods stop being valid

Newton needs J(x*) nonsingular: at a fold of the solution curve (the minimum-thrust point, if thrust rather than speed were the parameter) J is singular and the convergence degrades to linear or the iteration jumps branches. It needs a start in the basin: the arctan example shows that "close" can mean |x₀| < 1.39 for a function as tame as arctan. It needs a solution to exist: the overheated plate has none. Differenced Jacobians need F smooth at the scale of h: a tabulated drag polar with linear interpolation has a discontinuous derivative and a differenced J that jumps between table intervals, which stalls Newton at the table's resolution. Picard needs a contraction: its factor is of order |β|ΔT and it fails before Newton does. And the aircraft model is linear in α to 9°: beyond about 12° the real C_L bends over, two trims can exist at the same speed (the "back side" of the thrust curve), and the model, not the method, is what must change.

## 3. Build the computational model

Open `starter/main.py`. `Aircraft`, `isa_density`, `dynamic_pressure`, `trim_initial_guess`, `SolveResult`, `iterate_errors`, `convergence_order`, `stall_speed`, `arctan_residual`, `arctan_jacobian`, `face_conductances`, `HotPlate`, `plate_kfield` and `main` are given.

### Step 1: the trim equations

`trim_residual(x, ac, V, rho)` returns the three scaled residuals of 2.1 for x = (α, δe, T/W); write it with `np.sin`/`np.cos` so that complex arguments pass through. `trim_jacobian` is the 3 × 3 matrix of partial derivatives; check it with the two numerical Jacobians before using it.

### Step 2: Jacobians by differencing

```python
def fd_jacobian(f, x, h=None):
    f0 = f(x); J = np.empty((f0.size, x.size))
    for j in range(x.size):
        hj = np.sqrt(np.finfo(float).eps) * max(1.0, abs(x[j])) if h is None else h
        xp = x.copy(); xp[j] += hj
        J[:, j] = (f(xp) - f0) / hj
    return J
```

`complex_step_jacobian` perturbs `x.astype(complex)` by `1j * h` and takes `np.imag(f(xc)) / h`.

### Step 3: the solvers

`newton` is (1) in a loop with `np.linalg.solve`, recording ‖F‖ and the iterate at every step; stop at ‖F‖ ≤ tol or when the step is below tol relative to ‖x‖ (and ‖F‖ is already small); report `converged=False` if the iterate becomes non-finite or huge. `damped_newton` wraps the step in the backtracking loop of (3), recording λ. `broyden` keeps `B`, starts it from `fd_jacobian` unless `J0` is given, and applies (4). All three return a `SolveResult` with `nfev` counted honestly (the differenced Jacobian costs n evaluations).

### Step 4: trim

`solve_trim` builds the closures for the chosen method; `trim_sweep` loops over speeds, passing the previous solution as `x0` when `warm_start` is true and raising if any point fails.

### Step 5: the plate

`conduction_matrix` assembles the dense five-point matrix from the face conductances (`face_conductances` gives the harmonic means on vertical and horizontal faces; index the interior points with `np.arange(n*n).reshape(n, n)` and fill the five diagonals by fancy indexing); `conduction_boundary_rhs` collects k_face T_wall/h² on the four edges. `plate_residual` is A(T)T − b(T) − q. `picard_plate` solves the frozen linear system repeatedly; `newton_plate` calls `damped_newton` (or `newton`) with a differenced Jacobian.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-10-nonlinear-systems`. The tests:

1. **Atmosphere**: ρ(0) = 1.225, ρ(5 km) = 0.7364; q = ρV²/2.
2. **Residual** on a state with T = qS C_D/cos α: F₂ = 0 exactly, F₁ and F₃ as defined.
3. **Jacobians**: analytic = differenced to 10⁻⁶, = complex step to 10⁻¹³; on a 2 × 2 test map the differenced error has slope 1.0 ± 0.1 in h, the default step is accurate to 10⁻⁶, the complex step to 10⁻¹⁴.
4. **Newton on x² − 2**: converged to √2, Babylonian first iterate, order 2 ± 0.2; the differenced version counts 1 + 2k evaluations.
5. **arctan**: undamped diverges from 1.5 and 3.0 (|x| > 10, `converged=False`), converges from 1.0; damped converges from 1.5 and 3.0 with a shortened first step then unit steps.
6. **Order estimator**: 2.0 and 1.0 on synthetic sequences; raises with too few points.
7. **Trim**: Newton from (0, 0, 0) converged in ≤ 6, equal to `root(method="hybr")` to 10⁻¹⁰, residual < 10⁻¹², order 1.6–2.7, C_L = W/(qS) to 0.5 %, T = qS C_D to 1 %, 0 < δe < 5°; differenced Newton, Broyden and damped reach the same trim to 10⁻⁹; `nfev` = 1 + 4k for the differenced version; Broyden fewer evaluations, more steps, improving error ratios; damped takes unit steps; an unknown method raises.
8. **Initial guesses**: physical guess has δe = 0 and C_L = W/(qS); crude is zero; physical converges no slower.
9. **Sweep** 120–250 m/s: warm = cold to 10⁻⁸, no more iterations, α decreasing, δe increasing, minimum thrust within 10 m/s of the minimum-drag speed; `stall_speed` formula.
10. **Plate**: constant k gives the Kronecker Laplacian; random k gives a symmetric SPD matrix; no source gives T_wall; β = 0 converges in one Picard and one Newton step; at n = 8 Picard order 1 ± 0.15, Newton order > 1.6 in ≤ 6 steps, agreement 10⁻⁹, hotter than the linear model, rise 47 ± 3 K; equal to `newton_krylov` to 10⁻⁸.
11. `main` writes the figure.

### Compare

`scipy.optimize.root` (`hybr`) is MINPACK's `hybrd`/`hybrj`: Newton's direction inside a trust region (the dogleg step) with Broyden's update of the Jacobian between re-evaluations; it solved the trim in 11 evaluations against our 4 (Newton with the analytic Jacobian), 13 (differenced) and 9 (Broyden), with the same answer to 7 × 10⁻¹⁷. `method="lm"` is Levenberg–Marquardt, the right call when there are more equations than unknowns. `newton_krylov` is for large sparse systems: Jacobian-free, GMRES inside, optionally with a preconditioner (`@scipy-userguide`, "Root finding for large problems"), and it beat our dense differenced Newton on the plate by a factor 7 at n = 16 and would beat it by n² at larger n. We use the library for the trust region and for Krylov; we write Newton ourselves because every implicit integrator, finite-element nonlinearity and optimiser in the phases ahead contains one, and when it stalls, the diagnosis (start, Jacobian, existence) is the one this lesson practised.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Jacobian accuracy.** Max error of `fd_jacobian` against `trim_jacobian` at h = 10⁻², 10⁻⁴, …, 10⁻¹²; and of `complex_step_jacobian`. Predict: error ∝ h down to about 10⁻⁸, then rising as ε/h; complex step at round-off.
2. **The race.** From (0, 0, 0) at 200 m/s, 5 km: iterations, residual evaluations and final error for Newton (analytic), Newton (differenced), Broyden, damped Newton, `scipy.optimize.root`. Predict: 3–4 Newton steps with quadratic residual decay; Broyden more steps, fewer evaluations than differenced Newton.
3. **Continuation.** Sweep 110–260 m/s in 5 m/s steps, warm and cold; total iterations, agreement, the thrust minimum and its speed. Predict the minimum-drag speed from C_D0 = kC_L².
4. **Picard against Newton.** The plate at β = −0.002, −0.004, −0.006, −0.0075 K⁻¹, q = 5 kW/m²: iterations of each, maximum rise, k/k₀ at the hot spot; then q = 10 kW/m², β = −0.0075. Predict: Picard count grows with |β|, Newton's does not; the last case fails.

Raw data format: CSV tables (`h, fd_error, cs_error`; `method, iterations, nfev, final_error`; `V, alpha_deg, de_deg, T_kN, it_warm, it_cold`; `beta, q, it_picard, it_newton, max_rise_K, k_ratio`), each with a header naming the fixed parameters.

What you should find: errors 4.1 × 10⁻², 4.1 × 10⁻⁴, 4.1 × 10⁻⁶, 4.0 × 10⁻⁸, 1.9 × 10⁻⁶, 1.2 × 10⁻⁴ and 8.9 × 10⁻¹⁶; Newton 3 steps/4 evaluations, differenced 3/13, Broyden 5/9, damped 3/4, `root` 11 evaluations, all within 4 × 10⁻¹⁴; 94 against 101 iterations, agreement 8.6 × 10⁻¹³, minimum thrust 26.27 kN at 135 m/s against 134 m/s predicted; Picard 9, 10, 12, 13 against Newton 4, 4, 5, 5, rises 45.0, 47.5, 50.7, 53.9 K, k/k₀ 0.91, 0.81, 0.70, 0.60; the last case: Picard residual 2.3 × 10⁵ and damped Newton stalled at λ ≈ 10⁻⁶ with residual 3.7 × 10⁴.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| C_L at trim | W/(qS) = 0.3331 | 0.3327 | 0.1 % | T sin α carries 0.1 % of the weight |
| thrust at trim (kN) | qS C_D = 35.32 | 35.31 | 0.03 % | cos α and the C_L difference |
| α with δe = 0 (hand guess) | 0.865° | 0.700° at trim | 0.17° | the 2.2° elevator adds 0.015 to C_L |
| J (analytic vs complex step) | equal | 8.9 × 10⁻¹⁶ | — | round-off |
| optimal difference step | √ε ≈ 1.5 × 10⁻⁸ | best at 10⁻⁸ (4.0 × 10⁻⁸ error) | — | V-shape of truncation + cancellation |
| Newton order (trim) | 2 | 2.30 | 0.3 | two usable pairs only |
| Newton order (plate) | 2 | 1.93 | 0.07 | differenced Jacobian, finite tolerance |
| Picard order | 1 | 1.03 | 0.03 | asymptotic |
| arctan threshold | 1.3917 | diverges from 1.5, converges from 1.0 | — | consistent |
| Armijo steps from 1.5 / 3.0 | < 1 then 1 | 0.5, 1, 1, 1 / 0.25, 1, 1, 1 | — | as derived |
| minimum-drag speed (m/s) | 134 | 135 (5 m/s grid) | 1 | grid resolution; T sin α |
| warm vs cold sweep | equal | 8.6 × 10⁻¹³ | — | same root |
| plate rise, β = −0.004 (K) | > 42.97 (linear) | 47.48 | +10 % | k falls to 0.81 k₀ |
| Newton vs `newton_krylov` (K) | equal | 7 × 10⁻¹³ | — | same root, different linear solver |

Every row is a check against an exact statement or a limiting model, except the two orders measured from three or four points, which are the weakest numbers in the lesson and are reported with that caveat. The engineering answer is that the jet trims at 200 m/s with **α = 0.70°, δe = 2.22° (trailing edge down), T = 35.3 kN**, L/D = 11.1; the solver's contribution to its uncertainty is 10⁻¹⁶ and the model's (linear aerodynamics, ±10 % on C_D0) is 10 % on the thrust. The decision the numerics *do* inform is the cost of a sweep: 3 analytic-Jacobian Newton steps per point against 13 evaluations with differences, which matters when F is a CFD run.

Improvements worth making: an analytic Jacobian for the plate (the five-point pattern gives it in O(n²) instead of 257 residual evaluations), which turns 220 ms into a few; a Jacobian-free Newton–Krylov with the CG of `P03.09` (Extension 1); a stall model in C_L, which is where damping earns its place. Improvement worth *not* making: a tighter tolerance than 10⁻¹² on the trim; the aerodynamic coefficients are known to two digits.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (1)–(4) and of the optimal difference step; `nonlinear_systems.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the trim at 200 m/s, the evaluation cost of each method, the start you would use in a sweep, and the three failure modes you observed with their remedies.

### Adversarial self-checks

- **What if the weight doubles?** C_L doubles to 0.67, α to 4.3°, induced drag quadruples (thrust 53 kN by the hand estimate, 52.9 kN from the code); Newton still converges in 4 steps from (0, 0, 0) because the problem is still nearly linear, and the physical guess is still good. What changes is the model's validity margin: 4.6° is well within the linear range, but at 110 m/s the required C_L of 2.2 exceeds C_Lmax and the trim equations have a solution the aircraft does not.
- **What if the Jacobian is wrong by 10 %?** Newton becomes a linearly convergent iteration with factor 0.091 (measured with J scaled by 1.1): 12 iterations for 10⁻¹² instead of 3, still converging. Wrong by a sign in one entry, it diverges. The complex-step check catches both before the first solve.
- **Which assumption fails first in the real world?** Linear C_L(α) at the low-speed end of the sweep (α = 9°), then the constant C_D0 at the high-speed end (compressibility above about 230 m/s at 5 km, M ≈ 0.72); the thrust model (constant with speed) is the next.
- **What measurement could falsify the model?** Flight-test trim curves: δe against speed at fixed CG. A measured slope of δe/V that differs from the model's (7.9° over 150 m/s) by more than the ±10 % of the coefficients indicates a wrong C_mα or z_T; a measured thrust at 200 m/s of 40 kN rather than 35 points at C_D0.
- **Does the result survive a change of solver?** Newton (analytic), Newton (differenced), Broyden, damped Newton and `root` agree to 4 × 10⁻¹⁴ in the scaled unknowns; Picard, Newton and `newton_krylov` on the plate to 10⁻¹² K. Disagreement at that level would mean a wrong Jacobian, not a wrong method.
- **Which decision changes if the uncertainty doubles?** None of the numerics: with ±20 % on C_D0 the thrust is known to ±15 % and solving it to 10⁻¹² is still free. What changes is how the sweep should be run in an expensive model: at 13 evaluations per point for a differenced Jacobian against 1 per point for Broyden after the first, the quasi-Newton choice is what a CFD-based trim loop makes, and continuation is what keeps it inside the basin.

## Check yourself

1. Derive the Newton step for a system from the Taylor expansion and state the three conditions for quadratic convergence.
2. Why is the forward-difference step h ≈ √ε optimal, and what accuracy does it give? Why does the complex step do better?
3. Show that the Newton direction is a descent direction for ½‖F‖² and write the Armijo condition.
4. For F(x) = arctan x, why does Newton diverge from x₀ = 1.5 and converge from 1.0?
5. State Broyden's update and the two properties that define it. What convergence order does it give?
6. What is continuation, and when does it matter more than it did in the trim sweep?
7. Compare Picard and Newton on A(T)T = b(T): cost per iteration, order, and which fails first as the nonlinearity grows.
8. A Newton iteration stalls with the line search taking λ = 10⁻⁶. Name the three possible causes and how you would tell them apart.
