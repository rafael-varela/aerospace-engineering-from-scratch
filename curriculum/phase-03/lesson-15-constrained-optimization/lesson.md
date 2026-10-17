# P03.15  Constrained Optimization

> A constraint that is active at the optimum is not an obstacle but a price: its multiplier says exactly how much the objective would improve if you could relax it by one unit. SQP finds the point where every price is right.

| | |
|---|---|
| **Learning outcomes** | 1. Write the Lagrangian of an inequality-constrained problem, derive the Karush–Kuhn–Tucker (KKT) conditions, and interpret the multipliers as shadow prices (sensitivities of the optimum to the constraint limits). 2. Derive the sequential quadratic programming (SQP) step as Newton's method on the KKT system, state why the quadratic subproblem needs a positive definite Hessian and why a merit function is needed for global convergence. 3. Implement and test, from scratch, a small active-set QP solver, SQP with damped BFGS and with the exact Lagrangian Hessian, the quadratic penalty method and the augmented Lagrangian, with `scipy.optimize.minimize` (SLSQP, trust-constr) as the oracle. 4. Size a thin-walled compression strut for minimum mass under stress, Euler-buckling, shell-buckling and gauge limits, confirm the optimum by hand from the two active constraints, and decide from the multipliers which constraint is worth relaxing. |
| **Prerequisites** | `P03.14` (BFGS, line searches, convergence orders of Newton and quasi-Newton), `P03.10` (Newton for nonlinear systems), `P03.08` (solving the small saddle-point system), `P01.09` (gradients, Lagrange multipliers for equality constraints), `P01.13` (block matrices). Skill check: you can write the Lagrange condition ∇f = λ∇g for one equality constraint and solve a 2×2 linear system by hand. |
| **Engineering question** | A straight 6061-T6 aluminium tube 2.5 m long must carry 40 kN in compression between pinned ends. What mean radius and wall thickness give the lightest tube that neither yields, nor bows as an Euler column, nor wrinkles by local shell buckling, nor falls below a 0.5 mm manufacturing gauge; which of those limits actually set the design; and how many grams does each extra newton of load cost? |
| **Assumptions and limits** | Thin-walled circular tube (t ≪ r), so the area is 2πrt and the second moment πr³t; pinned ends with no eccentricity; a linear-elastic material with one allowable stress; Euler buckling of the whole column and the classical shell-buckling stress E t/(r√(3(1 − ν²))) multiplied by a knock-down factor of 0.5 for imperfections; no interaction between the two buckling modes and no fatigue, corrosion or joint design. The optimiser's assumptions: smooth functions, a regular active set (constraint gradients independent) and second-order sufficiency at the solution. |
| **Analytical work** | The Lagrangian and KKT conditions; the shadow-price interpretation d f*/d b_i = λ_i; the SQP step as Newton on the KKT system and its equivalence to a QP; the l₁ merit function and its directional derivative; the quadratic penalty's O(1/μ) bias and O(μ) conditioning; the augmented-Lagrangian multiplier update; the strut problem rewritten in log variables as a linear program, its vertex solution t² = P/(2πkE), r³t = PL²/(π³E), and the multipliers λ_Euler = λ_local = m*/3 by hand. |
| **Computational deliverable** | `starter/main.py`: `strut_mass`, `strut_mass_grad`, `strut_constraints`, `strut_jacobian`, `strut_hessian_lagrangian`, `analytic_optimum`, `kkt_residual`, `multipliers_least_squares`, `solve_qp`, `damped_bfgs_update`, `merit_l1`, `sqp`, `penalty_minimize`, `augmented_lagrangian`, `observed_order`, `observed_order_from_resolutions`, `shadow_price_check`, `design_chart` (the dataclasses, `unscale`, `strut_problem`, two textbook problems and `main` are given). Tests prove: gradients, Jacobian and Lagrangian Hessian match finite differences; the analytic vertex is feasible, has exactly two active constraints, satisfies the hand formulas and beats 200 random feasible designs; the QP solver reproduces hand solutions and `trust-constr` on random convex QPs and raises on indefinite or infeasible input; damped BFGS stays positive definite and satisfies the secant equation when undamped; SQP with BFGS and with the exact Hessian reaches the vertex from feasible and infeasible starts with the right multipliers, agrees with SLSQP to 10⁻⁷ and with trust-constr to 10⁻⁴, and shows an observed convergence order above 1.5 on its last three steps; the multipliers predict d m*/dP to 10⁻⁴ relative; the penalty error is first order in 1/μ (observed order 1.0 ± 0.1); the augmented Lagrangian with fixed μ = 10 reaches the same point and multipliers; the active set changes from Euler + gauge to Euler + local to stress + Euler as the load grows. |
| **Practical activity** | **computational experiment**: a sizing chart of optimal mass and active set against load from 5 to 160 kN; the penalty-parameter sweep and its observed order; SQP iteration counts from four starting points and in three variable scalings; the shadow price against a finite-difference re-optimisation. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `strut_sqp.png`, the hand derivation of the KKT conditions and of the vertex solution with its multipliers, the sizing-chart table, the penalty table with its fitted order, and a half-page sizing note stating the design, the governing constraints and the price of load. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: KKT conditions derived and the multipliers computed by hand to three figures; all public tests pass; the penalty order and the SQP order are measured and explained; the note names which constraint to relax and what it would buy. |
| **References** | `@ocw-16-90` Unit 3 "3.6 Introduction to Design Optimization" (required-free: problem statement, constraints, gradient-based search). `@ocw-18-335j` Week 11 "Lecture 29: Lagrange Dual Problems" (required-free: Lagrangian, KKT and duality); Week 9 "Lecture 26: Adjoint Methods" (supplementary-free: the multiplier as an adjoint variable, used again in `P03.16`). `@scipy-userguide` "Optimization (scipy.optimize)" > "Constrained minimization": "Sequential Least SQuares Programming (SLSQP) Algorithm", "Trust-Region Constrained Algorithm (method='trust-constr')" (required-free: the two oracles and their options). |
| **Extensions** | ★ Replace the enumerating QP solver by a primal active-set method that adds and drops one constraint per iteration and compare operation counts as m grows. ★ Add the interaction between Euler and local buckling through an effective modulus and see how the active set and the multipliers move. ★ Solve the strut as a linear program in log variables with `scipy.optimize.linprog` and show the dual variables equal the multipliers divided by m*. |
| **Time** | unmeasured |

## 1. The engineering problem

A compression strut (a wing-brace, a landing-gear drag link, a truss member in a launch-tower) is the simplest structural element that is governed by *several* failure modes at once. A thin tube is light, but a tube that is too thin for its radius wrinkles (local shell buckling), one that is too slender bows as a whole (Euler buckling), and one that is too small in cross-section yields. The designer's job is to pick the radius r and wall thickness t that minimise mass while every one of these modes stays inside its limit, and to know which limits are *governing*: those are the ones a better material or a shorter span would relieve, while the others would not change the design at all.

The lesson's strut is 2.5 m of 6061-T6 aluminium (E = 68.9 GPa, allowable 240 MPa, density 2700 kg/m³) carrying P = 40 kN between pinned ends, with a 0.5 mm manufacturing gauge and a 10 mm minimum radius for the end fittings. Five inequality constraints and a bilinear objective m = ρL·2πrt: two variables, so you can draw the design space, and simple enough that the optimum can be found by hand from the two constraints that turn out to be active. That hand solution is the oracle for the algorithm this lesson builds, **sequential quadratic programming**, which is Newton's method applied to the optimality conditions of a constrained problem and is what `scipy.optimize.minimize(method="SLSQP")` and most aircraft-sizing codes run inside.

Success criteria, written first: SQP from the heavy starting design (r = 50 mm, t = 2 mm, 4.24 kg) reaches the hand optimum to 10⁻⁸ in the scaled variables in fewer than 30 iterations from feasible and infeasible starts; its multipliers equal the hand multipliers m*/3; the multipliers predict the mass change per newton of load to 10⁻⁴ relative against a re-optimisation; the quadratic penalty method shows the first-order bias the theory predicts; and SLSQP and trust-constr land on the same point.

## 2. Understand and derive

### 2.1 The problem and its Lagrangian

Write every design problem in one standard form:

**minimise f(x) subject to g_i(x) ≥ 0, i = 1 … m.**   (1)

An equality h(x) = 0 is two inequalities, and a "≤" constraint is a "≥" after a sign change, so (1) loses nothing. Each g_i in this lesson is **normalised**, g_i = 1 − demand/capacity, so that g_i = −0.01 means "1 % over the limit" whether the limit is a stress in pascals or a thickness in metres; this matters because every algorithm below compares constraint values with each other and with the objective. The strut's five constraints are, with k = γ/√(3(1 − ν²)) = 0.3058 the knocked-down shell-buckling coefficient,

| i | mode | g_i(r, t) | form r^a t^b ≥ C |
|---|---|---|---|
| 0 | yield | 1 − P/(2πrt σ_a) | rt ≥ a_s = P/(2πσ_a) = 2.653 × 10⁻⁵ m² |
| 1 | Euler | 1 − PL²/(π³E r³t) | r³t ≥ a_e = PL²/(π³E) = 1.170 × 10⁻⁷ m⁴ |
| 2 | local shell | 1 − P/(2πkE t²) | t² ≥ a_l = P/(2πkE) = 3.021 × 10⁻⁷ m² |
| 3 | gauge | t/t_min − 1 | t ≥ 0.5 mm |
| 4 | fittings | r/r_min − 1 | r ≥ 10 mm |

(Euler: P ≤ π²EI/L² with I = πr³t. Local: the stress P/(2πrt) must not exceed kE t/r, and r cancels.) The objective is m = ρL 2πrt.

The **Lagrangian** attaches a price λ_i ≥ 0 to each constraint:

**L(x, λ) = f(x) − Σ_i λ_i g_i(x).**   (2)

The sign convention (minus, with λ ≥ 0 for g ≥ 0) makes the multipliers positive when a constraint pushes against the objective. At a local minimum x* of (1) at which the active constraint gradients are linearly independent, there exist multipliers λ* such that the **Karush–Kuhn–Tucker conditions** hold (`@ocw-18-335j` Lecture 29):

- stationarity: ∇f(x*) = Σ λ_i* ∇g_i(x*) = J(x*)ᵀλ*, where J is the m × n Jacobian of g;
- primal feasibility: g(x*) ≥ 0;
- dual feasibility: λ* ≥ 0;
- complementarity: λ_i* g_i(x*) = 0 for every i.

Complementarity is the heart of it: a constraint is either **active** (g_i = 0, it holds the design back, λ_i may be positive) or **inactive** (g_i > 0, it plays no role, λ_i = 0). The unconstrained condition ∇f = 0 of `P03.14` is the special case with no active constraints. Stationarity says that at the optimum the objective's gradient is a non-negative combination of the active constraint gradients: you cannot move to lower f without crossing a constraint. `kkt_residual` returns the largest violation of the four conditions in the infinity norm; it is zero at a solution and is the stopping criterion of every method in this lesson.

### 2.2 Multipliers are prices

Suppose the limit of constraint i is loosened: g_i(x) ≥ −ε instead of ≥ 0. Solve the perturbed problem and call the optimum f*(ε). Differentiating the Lagrangian along the path of optima, and using stationarity and complementarity to kill every term but one,

**d f*/dε = −λ_i*.**   (3)

So λ_i* is the **shadow price** of constraint i: how much objective one unit of constraint is worth. For a parameter p that appears inside the constraints (the load P, say) the same argument gives d f*/dp = ∂L/∂p = −Σ λ_i ∂g_i/∂p evaluated at the optimum, and for our normalised constraints ∂g_i/∂P = −(1 − g_i)/P, which is −1/P on an active constraint, so

**d m*/dP = (Σ_active λ_i)/P.**   (4)

`shadow_price_check` compares (4) with a central difference of two re-optimisations at P(1 ± 10⁻³): 0.023200 g/N both ways. This is the single most useful output of a constrained optimiser: the design is one number per variable, the multipliers are one number per limit, and the second set tells you where to spend engineering effort.

### 2.3 SQP: Newton's method on the KKT system

Pretend for a moment that the active set A is known and that all other constraints can be ignored. Then x* and λ_A* solve the square nonlinear system

∇f(x) − J_A(x)ᵀλ_A = 0, g_A(x) = 0,

and Newton's method (`P03.10`) on it gives, with W = ∇²_xx L = ∇²f − Σλ_i ∇²g_i the **Hessian of the Lagrangian**,

[[W, −J_Aᵀ], [J_A, 0]] [[p], [λ_A⁺]] = [[−∇f], [−g_A]].   (5)

Look at (5) differently: it is exactly the KKT system of the quadratic program

**minimise ½ pᵀWp + ∇fᵀp subject to g + Jp ≥ 0,**   (6)

the objective being a quadratic model of the Lagrangian and the constraints their linearisations at x. Solving (6) *with all m constraints as inequalities* lets the QP choose the active set for you, which is what makes SQP a method rather than a guess: at each iteration linearise, solve the QP for the step p and the new multipliers λ⁺, update x ← x + αp. Near a regular solution with second-order sufficiency the QP's active set settles to the true one, (6) reduces to (5), and the iteration converges **quadratically** when W is exact and **superlinearly** when W is a BFGS approximation, exactly as for Newton and quasi-Newton in `P03.14`. The run confirms it: the errors ‖x_k − x*‖ of the last three BFGS steps are 6.7 × 10⁻³, 1.2 × 10⁻⁴, 3.8 × 10⁻⁸, an observed order ln(3.8 × 10⁻⁸/1.2 × 10⁻⁴)/ln(1.2 × 10⁻⁴/6.7 × 10⁻³) = 2.0; with the exact Hessian 6.2 × 10⁻³, 1.0 × 10⁻⁴, 2.8 × 10⁻⁸, order 1.98. The multipliers must enter the Hessian: W is the Hessian of the Lagrangian, not of f, because the curvature of the constraints is what keeps the step on the feasible boundary. For the strut, ∇²m is the indefinite [[0, c], [c, 0]] (mass is bilinear), yet W at the solution has eigenvalues 2.13 and 14.06 in the scaled variables: the constraint curvature, weighted by λ, makes the Lagrangian convex where it matters.

Two things can go wrong, and the implementation handles both:

**W must be positive definite for (6) to have a unique bounded solution.** With `hessian="exact"` the Lagrangian Hessian is computed analytically and any eigenvalue below 10⁻⁶ is raised to it (a modified Newton step, inactive at the solution here). With `hessian="bfgs"` the matrix B starts as the identity and is updated with the Lagrangian gradient difference y = ∇L(x⁺, λ⁺) − ∇L(x, λ⁺); because the Lagrangian need not be convex, sᵀy can be negative, and **Powell's damping** replaces y by θy + (1 − θ)Bs whenever sᵀy < 0.2 sᵀBs so that the update keeps B positive definite (`damped_bfgs_update`; the tests check positive definiteness after twenty random updates and the secant equation Bs = y whenever damping is not triggered).

**The full step can be bad far from the solution.** A step that reduces f may increase constraint violation, or vice versa, so progress is measured by a **merit function** that combines both, the l₁ merit

**φ_μ(x) = f(x) + μ Σ_i max(0, −g_i(x)),**   (7)

with the penalty μ larger than the largest multiplier (we use μ ≥ 2 max λ⁺ + 10⁻³). The QP step is a descent direction for φ_μ, with directional derivative D = ∇fᵀp − μ Σ max(0, −g_i) ≤ 0, and a backtracking Armijo search accepts the first α ∈ {1, ½, ¼, …} with φ_μ(x + αp) ≤ φ_μ(x) + 10⁻⁴ α D. In the run from (50 mm, 2 mm) the merit function halves the first four steps (the errors go 1.5, 0.78, 0.39, 0.18, 0.068: the full QP step would overshoot the curved Euler boundary), and then full steps are taken and the superlinear phase begins. From the infeasible start (25 mm, 0.4 mm), which violates Euler, local buckling and the gauge at once, the same code converges in 8 iterations.

### 2.4 The QP subproblem by active-set enumeration

The QP (6) is itself a constrained problem, but a special one: convex (W positive definite) with linear constraints, so its KKT conditions are necessary *and* sufficient, and there are only finitely many possible active sets. For an active set 𝒲 the KKT system is linear, [[W, −J_𝒲ᵀ], [J_𝒲, 0]][p; λ_𝒲] = [−∇f; −g_𝒲], and a candidate is the answer if λ_𝒲 ≥ 0 and the inactive constraints are satisfied, g + Jp ≥ 0. `solve_qp` simply tries every subset of at most n constraints (a vertex in n dimensions cannot have more than n independent active constraints) and returns the first that passes. This is exact and transparent, and it is exponential in m: 2⁵ = 32 subsets for the strut, 2³⁰ for a thirty-constraint wing box. Production QP solvers do what `P03.09` did for linear systems, iterate: a primal active-set method adds or drops one constraint per step, an interior-point method follows a central path with all constraints slightly inactive. The lesson's version exists so that you can see that the QP's answer is a KKT point of a linear system, nothing more; extension 1 is the iterative version.

### 2.5 Two alternatives and why they lose

**Quadratic penalty.** Replace (1) by the unconstrained problem minimise φ(x) = f + (μ/2) Σ min(g_i, 0)² and solve it with BFGS (`penalty_minimize` calls `scipy.optimize.minimize(method="BFGS")`, the method of `P03.14`). The minimiser x_μ is always slightly infeasible, because the penalty only bites once a constraint is violated: setting the gradient to zero, ∇f = μ Σ(−g_i)∇g_i, and comparing with stationarity shows μ(−g_i) → λ_i, so the violation of an active constraint is −g_i ≈ λ_i/μ and the error in x is O(1/μ). The run gives ‖x_μ − x*‖ = 1.42 × 10⁻², 1.56 × 10⁻³, 1.57 × 10⁻⁴, 1.57 × 10⁻⁵ for μ = 10, 10², 10³, 10⁴, a fitted order 0.987 in 1/μ, and violations 4.6 × 10⁻³ (0.46 % over the buckling load) at μ = 100. To get the constraint to 10⁻⁸ you need μ = 10⁸, and the Hessian of φ then has eigenvalues of order μ along the constraint normals and of order 1 along them: a condition number of 10⁸ that makes the inner BFGS crawl. Bias and conditioning pull in opposite directions; that is the method's flaw.

**Augmented Lagrangian** (method of multipliers). Add the multiplier estimate to the penalty: for inequalities, L_A(x; λ, μ) = f + (μ/2) Σ [max(0, λ_i/μ − g_i)² − (λ_i/μ)²], minimise in x, then update λ_i ← max(0, λ_i − μ g_i(x)). The update is the identity μ(−g_i) → λ_i of the penalty method turned into an iteration, and it converges *linearly* at a rate that shrinks as μ grows, with μ held **fixed** so the conditioning never degrades. With μ = 10 the strut takes 7 outer iterations to a KKT residual of 10⁻⁸ and recovers the multipliers (0, 0.46401, 0.46401, 0, 0) to six figures. It is the method of choice when the QP of SQP is too large to solve, and it is what `P03.16`'s coupled problems fall back on when SQP's linearisations are unreliable.

### 2.6 The strut by hand: a linear program in disguise

Every row of the table in 2.1 is a **monomial** inequality r^a t^b ≥ C, and the objective is a monomial too. Take logarithms, u = ln r, v = ln t: the constraints become *linear*, au + bv ≥ ln C, and the objective ln(m/(2πρL)) = u + v is linear. The strut is a linear program in (u, v) (a **geometric program**, the class behind many sizing problems), and a linear program in two variables attains its minimum at a **vertex**, where two constraints are active. `analytic_optimum` solves each of the ten 2×2 pairs, discards the infeasible vertices and keeps the lightest.

For P = 40 kN the winning pair is Euler + local shell buckling:

t* = √a_l = √(3.021 × 10⁻⁷) = **0.5497 mm**, r* = (a_e/t*)^{1/3} = (1.170 × 10⁻⁷/5.497 × 10⁻⁴)^{1/3} = **59.71 mm**, m* = 2700 × 2.5 × 2π × 0.05971 × 5.497 × 10⁻⁴ = **1.3920 kg**.

Check the inactive constraints: the stress is 40 000/(2π × 0.05971 × 5.497 × 10⁻⁴) = 194.0 MPa, 19 % below the allowable (g₀ = 0.192); t is 9.9 % above the gauge (g₃ = 0.099); r/t = 108.6, thin-walled as assumed. The vertex's multipliers follow from stationarity in log variables: ∇(u + v) = (1, 1) = λ'_E (3, 1) + λ'_l (0, 2) gives λ'_E = λ'_l = 1/3, and converting from the log objective to the mass (dm = m d ln m, and on an active constraint dg_i = d ln(r^a t^b)) multiplies both by m*: **λ_Euler = λ_local = m*/3 = 0.46401 kg**. The code's least-squares multipliers agree to six figures, and (4) gives d m*/dP = (2/3)m*/P = 2.3200 × 10⁻⁵ kg/N, consistent with m* ∝ P^{2/3} on this vertex (t ∝ P^{1/2}, r ∝ P^{1/6}).

The active set changes with the load. From the vertex formulas, t* = √(P/(2πkE)) falls to the 0.5 mm gauge at P = 2πkE t_min² = 33.1 kN, below which Euler + gauge govern; and the stress P/(2πr*t*) ∝ P^{1/3} reaches 240 MPa at P = 40 × (240/194.0)³ = 75.8 kN, above which stress + Euler govern and m* = ρL a_s ∝ P exactly (2.25 kg at 80 kN, 4.50 kg at 160 kN). Three regimes, three different "what would help" answers: a thicker gauge limit is pure cost at 5 kN, a stiffer material (E) helps everywhere Euler is active, and a stronger alloy (σ_a) helps only above 76 kN.

### 2.7 Scaling

SQP with BFGS starts from B = I, so it assumes that a unit step in each variable is equally reasonable. In metres the strut variables are 0.06 and 0.0005: the identity is a terrible Hessian guess and the first steps are wasted learning the scale. The lesson therefore works in **scaled variables** x = (r/50 mm, t/1 mm), both O(1). The experiment in part 4 measures the effect here (small for two variables, decisive for twenty), and the same reasoning gave the normalised constraints: the merit function (7) adds a mass in kilograms to constraint violations, which only makes sense if the violations are dimensionless fractions.

### 2.8 Where the model stops being valid

- **Interaction of buckling modes.** A real tube near the Euler + local vertex fails below both loads: the modes couple through imperfections. Design codes handle this with interaction curves or lower knock-downs (γ = 0.25 moves the optimum to 0.777 mm and 1.754 kg; γ = 1 to the gauge at 1.307 kg).
- **Thin-wall formulas.** At r/t = 109 the area 2πrt is within 1 % of the exact annulus; below r/t ≈ 10 the formulas and the local-buckling model both fail.
- **Pinned ends, pure compression.** Any eccentricity or lateral load adds bending stress and makes the yield constraint active earlier.
- **Regularity.** If two active constraint gradients were parallel (a degenerate vertex) the multipliers would not be unique and SQP's QP would be singular. At our vertex the gradients (2.51, 1.82) and (0, 3.64) are far from parallel.
- **Local optimality.** KKT conditions are necessary, not sufficient, for a nonconvex problem. Here the log-variable form proves the solution is global; in general SQP finds *a* KKT point and you check others by restarting.

## 3. Build the computational model

Open `starter/main.py`. The dataclasses, `unscale`, `strut_problem`, `textbook_problem`, `rosenbrock_disk_problem` and `main` are given; everything else is yours.

### Step 1: the strut model

```python
def strut_constraints(x, params):
    r, t = unscale(x, params)
    p = params
    return np.array([
        1.0 - p.P / (2.0 * np.pi * r * t * p.sigma_allow),              # yield
        1.0 - p.P * p.L ** 2 / (np.pi ** 3 * p.E * r ** 3 * t),         # Euler
        1.0 - p.P / (2.0 * np.pi * p.k_local * p.E * t ** 2),           # local shell buckling
        t / p.t_min - 1.0,                                               # gauge
        r / p.r_min - 1.0,                                               # fittings
    ])
```

`strut_jacobian` differentiates each row with respect to r and t and multiplies the columns by (r_ref, t_ref) for the chain rule into the scaled variables; `strut_hessian_lagrangian` assembles ∇²m − Σλ_i∇²g_i the same way (the linear rows contribute nothing). Check both against finite differences before going on: every optimiser failure in this lesson that is not a bug in `sqp` is a wrong derivative.

### Step 2: KKT residual, multipliers, QP

```python
def solve_qp(H, c, A, b, tol=1e-10):
    n, m = H.shape[0], A.shape[0]
    for k in range(0, min(m, n) + 1):
        for W in itertools.combinations(range(m), k):
            W = list(W)
            K = np.block([[H, -A[W].T], [A[W], np.zeros((k, k))]])
            sol = np.linalg.solve(K, np.concatenate([-c, b[W]]))        # skip singular K
            p, lam_w = sol[:n], sol[n:]
            if np.all(lam_w >= -tol) and np.all(A @ p - b >= -tol * (1 + np.abs(b))):
                lam = np.zeros(m); lam[W] = lam_w
                return p, lam
    raise RuntimeError("infeasible linearisation")
```

Raise `ValueError` first if H is not positive definite (check the smallest eigenvalue). `kkt_residual` is the maximum of the four infinity norms in 2.1. `multipliers_least_squares(problem, x, active)` solves J_Aᵀλ_A ≈ ∇f with `np.linalg.lstsq` and puts zeros elsewhere.

### Step 3: SQP

```python
p, lam_new = solve_qp(B, g, J, -c)                       # c + J p >= 0
mu = max(mu, 2.0 * np.max(lam_new) + 1e-3)
phi0 = merit_l1(problem, x, mu)
d_phi = g @ p - mu * np.sum(np.maximum(0.0, -c))
alpha = 1.0
while merit_l1(problem, x + alpha * p, mu) > phi0 + 1e-4 * alpha * d_phi:
    alpha *= 0.5
x_new = x + alpha * p
s = x_new - x
y = (problem.grad(x_new) - problem.jac(x_new).T @ lam_new) - (g - J.T @ lam_new)
B = damped_bfgs_update(B, s, y)
```

Record `x` in `history` and the KKT residual in `kkt_history` every iteration, stop when the residual is below `tol`, and return an `SQPResult`. With `hessian="exact"` replace B by the problem's Lagrangian Hessian with its eigenvalues floored at 10⁻⁶.

### Step 4: penalty, augmented Lagrangian, analysis

`penalty_minimize` minimises f + (μ/2)Σ min(g, 0)² with `scipy.optimize.minimize(method="BFGS")` and an analytic gradient ∇f + μ Jᵀ min(g, 0). `augmented_lagrangian` loops the inner minimisation and the update λ ← max(0, λ − μg). `analytic_optimum` enumerates the ten constraint pairs in log variables. `observed_order` returns ln(e₃/e₂)/ln(e₂/e₁) from three successive iteration errors; `observed_order_from_resolutions` fits the slope of log error against log h. `shadow_price_check` and `design_chart` wrap `analytic_optimum` for part 4.

Running `main()` prints

```
analytic vertex: active = ['Euler', 'local'], r = 59.711 mm, t = 0.5497 mm, mass = 1.39202 kg
  stress at optimum 194.0 MPa (allowable 240), r/t = 108.6, multipliers = [0. 0.46401 0.46401 0. 0.] kg per unit g
SQP (bfgs Hessian): 9 iterations, converged = True, mass = 1.392020 kg, |x - x*| = 4.58e-16, KKT = 6.66e-16
  |x_k - x*| per iteration: 1.5e+00 7.8e-01 3.9e-01 1.8e-01 6.8e-02 6.7e-03 1.2e-04 3.8e-08 3.6e-15 4.6e-16
SQP (exact Hessian): 8 iterations, ...
scipy SLSQP: 11 iterations, mass = 1.392020 kg, |x - x*| = 6.75e-16
scipy trust-constr: 20 iterations, mass = 1.392034 kg, |x - x*| = 5.44e-06
shadow price d(mass*)/dP: predicted from multipliers 0.02320 g/N, finite difference 0.02320 g/N
```

followed by the penalty table, the augmented-Lagrangian line, the sizing chart and the figure.

### Verify

- **Analytic cases:** the textbook problem min x² + y² s.t. x + y ≥ 1 has the solution (½, ½) with λ = 1 from three starting points; the QP projection of the origin onto x + y ≥ 1 gives (½, ½) with λ = ½; the strut vertex satisfies t² = a_l and r³t = a_e to 10⁻¹⁰.
- **Limiting and failure cases:** an inactive QP constraint gets λ = 0 and the unconstrained minimiser; an indefinite H raises; contradictory linear constraints raise; `observed_order` of (10⁻¹, 10⁻², 10⁻⁴) is exactly 2.
- **Conservation and structure:** complementarity λ_i g_i = 0 and stationarity to 10⁻⁸ on random convex QPs; damped BFGS stays positive definite; the active set is Euler + gauge at 5 kN, Euler + local at 40 kN, stress + Euler at 160 kN, and mass doubles exactly between 160 and 320 kN.
- **Convergence:** the last three SQP errors give an observed order above 1.5 and a final ratio below 0.05 (superlinear) for both Hessians; the penalty error has fitted order 1.0 ± 0.1 in 1/μ over μ = 10², 10³, 10⁴.
- **Agreement with the oracle:** SLSQP to 10⁻⁷, trust-constr to 10⁻⁴ (its own tolerance) on the strut and on Rosenbrock inside a disk, where the constraint is active with a positive multiplier.

### Compare

`scipy.optimize.minimize(method="SLSQP")` is an SQP: a least-squares formulation of the QP, a BFGS Lagrangian Hessian, an l₁ merit function; it takes 11 iterations where ours takes 9, and reaches the same point to 10⁻¹⁵. `trust-constr` is an interior-point method with a trust region: it approaches the vertex from inside the feasible set, takes 20 iterations, and stops 5 × 10⁻⁶ short because its barrier parameter stops at its own tolerance; this is normal for interior-point methods and why their answers sit a hair inside the constraints. We call the library in practice because it handles equality constraints, bounds, degenerate active sets and infeasible QP subproblems (SLSQP's least-squares formulation never fails to return a step; ours raises), and because its QP solver is iterative. What writing our own bought is the knowledge that an "SLSQP converged" message means a KKT residual below a tolerance in *scaled* variables with *normalised* constraints: feed it radii in metres and stresses in pascals and the same 11-iteration problem can take hundreds or stop early at a false optimum.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Sizing chart.** Run `design_chart` for P = 5, 10, 20, 40, 80, 160 kN. Predict the two load thresholds 2πkE t_min² and 40 (σ_a/194.0)³ kN at which the active set changes, and the local exponents of m*(P): 1/3 on Euler + gauge (t fixed, r ∝ P^{1/3}), 2/3 on Euler + local, 1 on stress + Euler. Record `P_kN, mass_kg, active_pair`.
2. **Penalty order.** Run `penalty_minimize` for μ = 10, 10², 10³, 10⁴ and tabulate ‖x_μ − x*‖ and the largest violation. Predict the slope 1 on log-log axes and the violation λ/μ ≈ 0.464/μ in normalised units (4.6 × 10⁻³ at μ = 100).
3. **Starting points and scaling.** Count SQP iterations to KKT < 10⁻¹⁰ from (50 mm, 2 mm), the infeasible (25 mm, 0.4 mm), (100 mm, 1 mm) and (75 mm, 0.6 mm), for both Hessians; then rebuild the problem with `r_ref = t_ref = 1` (metres) and with `r_ref = t_ref = 10⁻³` (millimetres) and repeat from the first point with BFGS. Predict that the exact Hessian never needs more iterations than BFGS and that the unscaled runs need more.
4. **Shadow price.** Compare `shadow_price_check` with a finite difference of `analytic_optimum` at δ = 10⁻², 10⁻³, 10⁻⁴; predict agreement improving as δ² until rounding.

Raw data format: CSV tables (`P_kN, mass_kg, active`; `mu, error, violation`; `start, hessian, scaling, iterations`), each with a header naming the fixed parameters; `strut_sqp.png`.

What you should find: thresholds at 33.1 and 75.8 kN; masses 0.653, 0.823, 1.037, 1.392, 2.250, 4.500 kg (the 5 → 10 → 20 kN ratios are 2^{1/3} = 1.26 as predicted, 80 → 160 kN exactly 2); the penalty errors 1.42 × 10⁻², 1.56 × 10⁻³, 1.57 × 10⁻⁴, 1.57 × 10⁻⁵ (slope 0.987) with violations 4.3 × 10⁻², 4.6 × 10⁻³, 4.6 × 10⁻⁴, 4.6 × 10⁻⁵; iteration counts 9/8, 8/8, 7/7, 6/5 (BFGS/exact) for the four starts, and 10 (metres) and 13 (millimetres) against 9 (scaled) for BFGS: the effect of scaling is real but modest with two variables; the shadow price agrees to 7 × 10⁻⁶, 7 × 10⁻⁸ and 7 × 10⁻¹⁰ relative at δ = 10⁻², 10⁻³, 10⁻⁴ (the δ² law of a central difference).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| t*, r* | 0.5497 mm, 59.71 mm | 0.5497 mm, 59.711 mm | < 10⁻¹⁰ relative | same formulas |
| m* | 1.3920 kg | 1.392020 kg | 0 | — |
| λ_Euler, λ_local | m*/3 = 0.46401 kg | 0.46401, 0.46401 | < 10⁻⁷ | least squares at the vertex |
| d m*/dP | (2/3) m*/P = 2.3200 × 10⁻⁵ kg/N | 2.32003 × 10⁻⁵ (FD) | 10⁻⁷ relative | central difference O(δ²) |
| SQP vs vertex | 0 | 4.6 × 10⁻¹⁶ | rounding | — |
| SQP vs SLSQP | 0 | < 10⁻¹⁵ | rounding | same algorithm family |
| SQP vs trust-constr | 0 | 5.4 × 10⁻⁶ | barrier parameter | interior-point tolerance |
| SQP order (last three steps) | 2 (exact), superlinear (BFGS) | 1.98, 2.0 | — | asymptotic regime reached |
| penalty order in 1/μ | 1 | 0.987 | 1 % | finite range of μ |
| stress at optimum | 194 MPa < 240 | 194.0 MPa | 0 | constraint inactive, as the vertex test requires |

The numerical side of this lesson validates completely: five independent routes (hand vertex, SQP with two Hessians, augmented Lagrangian, SLSQP, trust-constr) agree to the tolerances each one claims. The engineering validation is a different matter, and the multipliers tell you where to look. Both active limits are *buckling* limits, whose formulas carry the largest modelling uncertainty in the whole problem: the knock-down γ = 0.5 is an engineering judgement, and the interaction between the two modes is ignored. The honest statement of the result is therefore "1.39 kg with r = 60 mm, t = 0.55 mm, governed by stability, with the mass sensitive to the shell-buckling knock-down roughly as m* ∝ γ^{−1/3}" (γ = 0.25 gives 1.754 kg), and the next step is not a better optimiser but a buckling test or a validated interaction curve. Note also what the optimum is *not* sensitive to: the allowable stress could drop to 200 MPa without changing the design, and a tube made of a 20 % stronger alloy would weigh exactly the same.

Improvements worth making: an iterative active-set QP (extension 1) so that m can grow; a mode-interaction constraint (extension 2); the LP view of the problem (extension 3), which also shows that sizing problems of this monomial type are globally solvable. Improvement worth *not* making: pushing the penalty parameter to 10⁸ to make the penalty method "accurate". The augmented Lagrangian gets the same accuracy at μ = 10.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of the KKT conditions, of (3)–(4), of the QP equivalence (5)–(6), and of the vertex solution with its multipliers; `strut_sqp.png`; the four tables of part 4; and a half-page sizing note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the design, the two governing constraints, the price of load in g/N, the change that would *not* help (a stronger alloy) and the assumption you trust least.

### Adversarial self-checks

- **What if the load doubles?** From 40 to 80 kN the active set changes: stress joins Euler, the thickness jumps from 0.55 to 0.80 mm and the mass from 1.39 to 2.25 kg (×1.62, not ×2 and not ×2^{2/3} = 1.59, because the vertex moved). Multipliers computed at 40 kN predict the first few per cent of the change and nothing beyond the regime boundary at 75.8 kN: a shadow price is a derivative, not a forecast.
- **What if the length doubles?** a_e grows four-fold, r* by 4^{1/3} = 1.59 to 94.8 mm, the mass to 4.42 kg; t* does not move at all, because the local-buckling constraint does not contain L. Which of the two active constraints would you attack first for a long strut?
- **Which assumption fails first in the real world?** The independence of the two buckling modes, then the pinned-end idealisation (a 1 % eccentricity at r = 60 mm adds a bending stress of order P e/(πr²t), tens of megapascals here), then the thin-wall formulas; never the optimiser.
- **What measurement could falsify the model?** A compression test of a tube built to the design: failure below 40 kN, or failure by a combined wrinkle-and-bow mode, falsifies the no-interaction assumption; the lesson's numbers are predictions, not test results.
- **Does the result survive a change of solver or scaling?** SQP (two Hessians), augmented Lagrangian, SLSQP and trust-constr agree; the penalty method lands 10⁻³ away at μ = 100 and *says so* through its visible constraint violation; scaling in metres or millimetres changes the iteration count (10, 13 against 9), not the answer.
- **Which decision changes if the uncertainty doubles?** If the shell-buckling knock-down is uncertain by a factor 2 (γ between 0.25 and 0.5), the mass is uncertain between 1.39 and 1.75 kg and the thickness between 0.55 and 0.78 mm. A designer who must freeze the drawing today picks the thicker wall; one who can afford a buckling test on three tubes buys back up to 0.36 kg per strut. The multipliers say that a test of the *Euler* load would buy nothing: its formula is exact for an ideal pinned column.

## Check yourself

1. Write the Lagrangian and the four KKT conditions for minimise f(x) subject to g(x) ≥ 0, and say what complementarity means for an inactive constraint.
2. A constraint at the optimum has multiplier λ = 0.46 kg. By how much does the minimum mass change if the constraint limit is loosened by 1 % of its normalised value, and in which direction?
3. Why is the Hessian in the SQP subproblem that of the Lagrangian and not that of the objective? Give the strut as an example.
4. Why does the quadratic penalty method always return a slightly infeasible point, and how does the violation scale with μ?
5. What does Powell's damping protect, and when is it triggered?
6. The strut's optimum sits where Euler and local buckling are active. What happens to the design if the allowable stress is raised by 20 %? If the gauge limit is raised to 0.6 mm?
7. State the observed convergence order of SQP with an exact Hessian and with BFGS, and how you measure it from three successive errors.
8. Why does rewriting the strut in logarithmic variables prove that the vertex found is the global optimum?
