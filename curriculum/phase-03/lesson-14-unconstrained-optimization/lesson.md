# P03.14  Unconstrained Optimization

> An optimiser is only as good as its model of the function's curvature. Steepest descent assumes none, Newton computes all of it, and BFGS learns it from the gradients it has already paid for. That is the whole lesson; the rest is making each claim measurable.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the first- and second-order optimality conditions, the exact-line-search steepest-descent rate ((κ − 1)/(κ + 1))² on a quadratic, the Wolfe conditions, the secant equation and the BFGS inverse-Hessian update, and show that the update preserves symmetry and positive definiteness when yᵀs > 0. 2. Explain why steepest descent converges linearly, Newton quadratically and BFGS superlinearly, and why BFGS on an n-dimensional quadratic with exact line search terminates in n steps with the exact inverse Hessian. 3. Implement and test, from scratch, a central-difference gradient, Armijo and strong-Wolfe line searches, steepest descent, damped Newton and BFGS, against closed-form answers and `scipy.optimize.minimize`, and measure the observed orders (2 for the finite-difference gradient, linear, quadratic and superlinear for the three methods). 4. Find the minimum-fuel cruise Mach number and altitude of a jet transport, explain the optimum from the physics of the drag polar and the isothermal stratosphere, and state what an unconstrained optimum does and does not tell an operator. |
| **Prerequisites** | `P01.09` (gradients and the idea of a gradient search), `P03.02` (Newton's method for root finding, convergence order and its measurement), `P03.10` (Newton's method for systems), `P03.08` (symmetric positive-definite matrices, Cholesky), `P03.04` (finite-difference derivatives and the truncation/rounding trade-off), `P00.03` (the ISA atmosphere), `P01.02` (the drag polar and L/D), `P01.16` and `P01.18` (eigenvalues, quadratic forms, condition number). Skill check: you can write the Taylor expansion of a function of several variables to second order and you know what a positive-definite matrix is. |
| **Engineering question** | A 60-tonne narrow-body airliner (W = 600 kN, S = 125 m², C_D = 0.020 + 0.045 C_L² plus wave drag above a drag-divergence Mach number that falls with lift coefficient, TSFC 1.6 × 10⁻⁵ kg/(N s) at sea level scaling with √T) is in mid-cruise. At what Mach number and altitude does it burn the least fuel per kilometre, how much does it save against "M 0.78 at 10 km", how flat is the optimum, and how does the answer move as fuel burns off? Alongside: which minimiser finds that point in the fewest function evaluations, and how do you know it has found it? |
| **Assumptions and limits** | Steady level cruise with lift = weight and thrust = drag; a parabolic polar with a smooth quartic wave-drag rise (a textbook model, not a measured polar); TSFC ∝ √T and independent of Mach and throttle; ISA with no wind; the optimum sought without constraints (no ceiling, no buffet margin, no M_MO, no airway or ATC level). The optimisers assume a twice-differentiable objective with a positive-definite Hessian at the minimum; they find a local minimum, which is the global one only when the function has a single basin (true here, not in general). |
| **Analytical work** | Optimality conditions; the quadratic model and the linear rate of steepest descent; Newton's method and quadratic convergence; the Wolfe conditions; the secant equation, the BFGS update and its properties; finite termination on quadratics; the truncation–rounding analysis of the central-difference gradient; the cruise objective reduced to (M, C_L) in the isothermal stratosphere and the resulting cruise-climb law. |
| **Computational deliverable** | `starter/main.py`: `rosenbrock`, `rosenbrock_grad`, `rosenbrock_hess`, `quadratic`, `gradient_fd`, `hessian_fd`, `line_search_backtracking`, `line_search_wolfe`, `steepest_descent`, `newton_method`, `bfgs_update`, `bfgs`, `error_sequence`, `convergence_order`, `linear_rate_bound`, `fuel_per_km`, `cruise_breakdown`, `optimal_cruise` (`isa`, `cruise_aircraft` and `main` are given). Tests prove: the central-difference gradient has observed order 2.00 ± 0.05 and its error turns up again below h ≈ 10⁻⁸; the finite-difference Hessian matches the analytic one; both line searches satisfy their conditions and reject ascent directions; the BFGS update satisfies the secant equation to 10⁻¹² and keeps H symmetric positive definite; BFGS with exact line search finishes a 3-variable quadratic in 3 steps with H = Q⁻¹ to 10⁻¹⁰; steepest descent's f-ratios never exceed ((κ − 1)/(κ + 1))²; Newton's error is squared at each late step on Rosenbrock; BFGS's error ratios fall below 0.2 and its answer matches `scipy.optimize.minimize(method='BFGS')` to 10⁻⁶; the cruise optimum matches SciPy's to 10⁻⁴ in Mach and 10 m in altitude and has a positive-definite Hessian. |
| **Practical activity** | **computational experiment**: a step-size study for the finite-difference gradient; a method comparison on Rosenbrock counting function and gradient evaluations; a scaling study on the cruise problem (altitude in metres against tens of kilometres); a cruise-climb study with weight. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `bfgs_cruise.png`, the hand derivations of the steepest-descent rate, the BFGS update and its positive-definiteness, the four tables of §4, and a half-page note stating the optimum cruise point with its sensitivity and the constraints that would move it. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the BFGS update is derived and its properties proved; all public tests pass; the three convergence rates are measured and explained; the note gives the cruise optimum with the 1 % window and names the constraint most likely to bind. |
| **References** | `@ocw-18-335j` Week 9 "Lecture 25: Overview of Optimization Algorithms" (required-free: the landscape of methods, local against global, why gradients matter); Week 11 "Lecture 30: Quasi-Newton Methods and the BFGS Algorithm", "Lecture 31: Derivation of the BFGS Update" (required-free: the secant condition and the update used here); Week 10 "Lecture 28: Trust-Regions Methods and the CCSA Algorithm" and Week 12 "Lecture 32: Derivative-Free Optimization by Linear and Quadratic Approximations" (supplementary-free: what to use when line searches or gradients are unavailable). `@scipy-userguide` "Optimization (scipy.optimize)" > "Unconstrained minimization": "Broyden-Fletcher-Goldfarb-Shanno algorithm" (required-free: the oracle's options and stopping rules), "Nelder-Mead Simplex algorithm", "Newton-Conjugate-Gradient algorithm" (supplementary-free: the alternatives, for the comparison in §3). |
| **Extensions** | ★ Add the L-BFGS two-loop recursion (`@ocw-18-335j` Lecture 30) and minimise a 1000-variable discretised energy without ever forming H. ★ Add a bounded cruise problem (ceiling 12.5 km, M ≤ 0.80) through a penalty or by `scipy.optimize.minimize(method='L-BFGS-B')` and show which constraint binds for the heavy and the light aircraft. ★ Replace the finite-difference gradient by complex-step differentiation and show the step-size trade-off disappears. |
| **Time** | unmeasured |

## 1. The engineering problem

Every cruise of every jet transport is an optimisation problem solved several times an hour: for the current weight, which Mach number and flight level burn the least fuel per kilometre? The answer sits on a ridge between two drag mechanisms. Fly too slow or too low and the induced drag K C_L² grows; fly too fast and the wave drag rises as the fourth power of the Mach excess over drag divergence, and drag divergence itself retreats as C_L rises. Fuel per kilometre is c_T D/V: the specific fuel consumption times drag over speed, with c_T falling as √T up to the tropopause, V = Ma falling with the speed of sound, and D following the polar. Nothing in that function is quadratic, nothing is separable, and its minimum has no closed form. It is the ordinary situation of design: a smooth objective of a few variables, cheap to evaluate, whose gradient is not available analytically.

The **BFGS method** (`@ocw-18-335j` Lecture 30) is the workhorse for exactly this situation, and `scipy.optimize.minimize` uses it by default when no Hessian is given. This lesson builds it from the ground up, together with the two methods it sits between: **steepest descent**, which uses the gradient alone and crawls down ill-conditioned valleys, and **Newton's method**, which uses the Hessian and converges quadratically but needs second derivatives and a safeguard far from the minimum. Each is tested where the answer is known (quadratics, the Rosenbrock function) and then applied to the cruise problem, where the oracle is SciPy and the physics is the check.

Success criteria, written down first: on a quadratic with condition number κ, steepest descent's energy ratio per step is bounded by ((κ − 1)/(κ + 1))² (a theorem; the measured 0.7801 against 0.7805 for κ = 16.2 must not exceed it); BFGS with exact line search finishes a 3-variable quadratic in 3 steps with the exact inverse Hessian; on Rosenbrock from (−1.2, 1) Newton takes about 20 iterations, BFGS about 35 and steepest descent more than 10 000; the cruise optimum agrees with SciPy to 10⁻⁴ in Mach and 10 m in altitude and makes physical sense (C_L a little below the maximum-L/D value, M a little above drag divergence).

## 2. Understand and derive

### 2.1 What a minimum looks like

For smooth f: ℝⁿ → ℝ, the Taylor expansion about x is f(x + p) = f(x) + gᵀp + ½pᵀHp + O(|p|³) with gradient g = ∇f and Hessian H = ∇²f. At a local minimiser x* the **first-order condition** g(x*) = 0 is necessary (otherwise p = −εg lowers f); the **second-order condition** H(x*) positive semi-definite is necessary and, if positive definite, sufficient. A direction p with gᵀp < 0 is a **descent direction**: f falls along it for small enough steps. All three methods of this lesson have the form x_{k+1} = x_k + α_k p_k with p_k a descent direction and α_k > 0 chosen by a **line search**; they differ only in p_k.

### 2.2 Steepest descent and the price of ignoring curvature

p = −g is the direction of fastest local decrease. On the quadratic f = ½xᵀQx − bᵀx (Q symmetric positive definite, minimiser x* = Q⁻¹b) the exact step along −g is α = gᵀg/(gᵀQg), and Kantorovich's inequality gives (`@ocw-18-335j` Lecture 25)

**(f_{k+1} − f*)/(f_k − f*) ≤ ((κ − 1)/(κ + 1))²,  κ = λ_max(Q)/λ_min(Q).**   (1)

The rate is **linear** (a fixed factor per step) and the factor approaches 1 as κ grows: for κ = 100 it is 0.961, so 10⁻¹⁰ accuracy costs about 580 iterations, and the iterates zig-zag across the valley because each step is orthogonal to the previous one (g_{k+1}ᵀg_k = 0 at an exact line search). `main()` measures 0.7801 for a 3 × 3 matrix with κ = 16.2 (bound 0.7805), 189 iterations to 10⁻¹⁰; on Rosenbrock, whose Hessian at the minimum has κ ≈ 2500, steepest descent needs 14 786 iterations from (−1.2, 1). Scaling the variables changes Q, hence κ, hence everything: that is why the cruise variables are (M, h/10 km) and not (M, h in metres), for which the Hessian's condition number is 1.4 × 10⁹.

### 2.3 Newton's method and quadratic convergence

Minimising the quadratic model gives the **Newton step** p = −H⁻¹g, which lands on the minimiser of a quadratic in one step and, near a minimiser with H(x*) positive definite, converges **quadratically**: |e_{k+1}| ≤ C|e_k|² (`P03.02` for the one-variable case; here e_k = x_k − x*). Two things go wrong far away. H may not be positive definite, so the Newton step may point uphill; the solver checks with a Cholesky attempt and falls back to −g when it fails. And the unit step may overshoot; an **Armijo backtracking** search, halving α until f(x + αp) ≤ f(x) + c₁αgᵀp with c₁ = 10⁻⁴, protects against that. On Rosenbrock from (−1.2, 1) the damped Newton method takes 22 iterations, of which the last four reduce the error from 9 × 10⁻³ to 1 × 10⁻³, 3 × 10⁻⁶ and 1 × 10⁻¹⁰: each error roughly the square of the previous once the full step is accepted. Pure Newton without the line search takes 6 iterations from (−1.3, 1) but passes through f = 2197 on the way, forty times its starting value: it converges here by luck, and in general it does not.

### 2.4 The Wolfe conditions

A line search need not be exact; it needs to guarantee progress and to supply usable curvature information. The **strong Wolfe conditions** (`@ocw-18-335j` Lecture 30) on φ(α) = f(x + αp) are

sufficient decrease φ(α) ≤ φ(0) + c₁αφ'(0),  curvature |φ'(α)| ≤ c₂|φ'(0)|,  0 < c₁ < c₂ < 1,   (2)

with c₁ = 10⁻⁴ and c₂ = 0.9 for quasi-Newton methods. The first rules out steps too long to help; the second rules out steps too short to matter and, because it forces φ'(α) ≥ c₂φ'(0) > φ'(0), it implies (g_{k+1} − g_k)ᵀp > 0, which is the **curvature condition yᵀs > 0** that BFGS needs. The search here brackets (expanding α until the conditions fail or the slope turns positive) and then zooms by bisection; it is the structure of every production line search, with polynomial interpolation left as an improvement.

### 2.5 BFGS: learning the curvature from the gradients

A quasi-Newton method keeps an approximation H_k to the inverse Hessian, takes p = −H_k g_k, and updates H after each step from the pair s = x_{k+1} − x_k, y = g_{k+1} − g_k. For a quadratic, y = Qs exactly, so the natural requirement is the **secant equation** H_{k+1}y = s. In n > 1 dimensions it under-determines H_{k+1}; the BFGS choice is the symmetric matrix closest to H_k in a weighted Frobenius norm that satisfies it (`@ocw-18-335j` Lecture 31), and it has the closed form

**H_{k+1} = (I − ρsyᵀ) H_k (I − ρysᵀ) + ρssᵀ,  ρ = 1/(yᵀs).**   (3)

Three properties, each tested. It satisfies the secant equation: multiply (3) by y, the first term's (I − ρysᵀ)y = y − y = 0 and the second gives ρs(yᵀs) = s. It preserves symmetry, by inspection. And it preserves positive definiteness when yᵀs > 0: for any z ≠ 0, zᵀH_{k+1}z = wᵀH_kw + ρ(sᵀz)² with w = (I − ρysᵀ)z; the first term is non-negative, and if it is zero then w = 0, so z = ρ(sᵀz)y, which forces sᵀz ≠ 0 (otherwise z = 0) and makes the second term positive. The Wolfe line search guarantees yᵀs > 0, so H stays a valid metric and p = −Hg stays a descent direction. The first H₀ is the identity, rescaled after the first step by yᵀs/yᵀy so that its size matches the true inverse curvature.

Two consequences are measurable. On an n-dimensional quadratic with exact line searches, BFGS generates conjugate directions and terminates in at most n steps with H_n = Q⁻¹ exactly (`main()`: 3 steps, |H − Q⁻¹| = 6 × 10⁻¹⁷). On a general function it converges **superlinearly**: |e_{k+1}|/|e_k| → 0, without ever reaching a fixed order 2. On Rosenbrock from (−1.2, 1) it takes 36 iterations (SciPy: 34) and its last error ratios are 0.60, 0.07, 0.023, 0.011, 0.037: tending to zero but irregularly, not squaring. Its final H differs from the exact inverse Hessian [[0.5, 1], [1, 2.005]] by 3 × 10⁻⁵, a reminder that BFGS learns the curvature along the directions it has travelled and nowhere else.

### 2.6 Gradients by finite differences

When no analytic gradient exists the central difference (f(x + he_i) − f(x − he_i))/(2h) has truncation error h²f'''/6 and rounding error about ε|f|/h, so the error is smallest near h ≈ (3ε|f|/|f'''|)^{1/3} ~ 10⁻⁵ to 10⁻⁶ for a well-scaled function (`P03.04`). `main()` measures it on Rosenbrock at (0.3, −0.7): 1.2 × 10⁻², 1.2 × 10⁻⁴, 1.2 × 10⁻⁶ for h = 10⁻², 10⁻³, 10⁻⁴ (order 2 exactly), 10⁻⁸ at h = 10⁻⁵ and 10⁻⁶, then 6 × 10⁻⁷ at h = 10⁻⁸ and 5 × 10⁻⁵ at h = 10⁻¹⁰ as rounding takes over. Each gradient costs 2n function evaluations; the cruise optimisation needs 9 BFGS iterations and 135 evaluations of the model, which is why the function evaluation count, not the iteration count, is the cost that matters.

### 2.7 The cruise problem reduced by hand

In the isothermal stratosphere (above 11 km, T = 216.65 K) a and c_T are constants, so fuel per km = c_T D/V = c_T W/(a M (L/D)), and L/D = C_L/C_D depends only on (M, C_L). The minimum over (M, h) is therefore a minimum over (M, C_L), **independent of weight**, and the altitude follows from lift = weight: p = 2W/(γM²SC_L). Hence as fuel burns off the optimum Mach and C_L stay fixed and the optimal altitude rises along a **cruise-climb**, dh*/dW = −RT/(g₀W) = −10.6 m per kN at 600 kN. Without wave drag the optimum C_L would be √(C_D0/K) = 0.667 (maximum L/D = 16.67) and the fuel per km would fall without limit as M rose. The wave-drag term k(M − M_dd(C_L))⁴ with M_dd = 0.82 − 0.20 C_L stops the Mach rise and pulls C_L down: at the optimum C_L = 0.585, L/D = 16.14, M = 0.776 with M − M_dd = 0.073 and C_D,wave = 8.5 × 10⁻⁴ (2.4 % of the total C_D of 0.0362), which is the balance point where a further Mach increase costs more in wave drag than it saves in time. The altitude from the pressure condition is p = 19 469 Pa, h = 11 955 m; BFGS finds 11 955 m. Within 1 % of the minimum fuel the window is ±0.023 in Mach and ±715 m in altitude (from the Hessian: Δf/f ≈ ½H₁₁ΔM²/f).

### 2.8 Where the model stops being valid

The polar is a textbook model: a real aircraft's drag rise is measured, not quartic, and its TSFC depends on Mach and throttle. The unconstrained optimum ignores the ceiling (thrust and buffet limits), M_MO, the flight levels ATC offers (a step-climb, not a cruise-climb), wind (which shifts the ground-speed optimum) and the cost of time, which airlines weigh through a cost index that moves the optimum Mach up. The methods have their own limits: all three find the nearest local minimum, none guarantees a global one, and a finite-difference gradient at a kink (the wave-drag term is C³ at M = M_dd, smooth enough here) or at a noisy objective can send BFGS astray.

## 3. Build the computational model

Open `starter/main.py`. `isa`, `cruise_aircraft` and `main` are given.

### Step 1: derivatives and line searches

```python
def gradient_fd(f, x, h=1e-6):
    g = np.empty(x.size)
    for i in range(x.size):
        e = np.zeros(x.size); e[i] = h
        g[i] = (f(x + e) - f(x - e)) / (2 * h)
    return g

def line_search_backtracking(f, x, p, g, alpha0=1.0, rho=0.5, c1=1e-4, max_iter=60):
    slope = g @ p                              # must be negative: raise ValueError otherwise
    f0, alpha = f(x), alpha0
    for _ in range(max_iter):
        if f(x + alpha * p) <= f0 + c1 * alpha * slope:
            return alpha
        alpha *= rho
    return alpha
```

`line_search_wolfe` implements the bracketing–zoom structure of §2.4: expand α (doubling, up to `alpha_max`) until sufficient decrease fails, φ stops falling or φ' turns non-negative, then zoom by bisection on the bracket, keeping the end with the lower φ. Return the step as soon as both conditions in (2) hold.

### Step 2: the three minimisers

Each takes `f`, `grad`, `x0`, `tol`, `max_iter`, stops when ‖g‖ ≤ tol, and returns the result dictionary of the module docstring. `steepest_descent` uses p = −g with the Wolfe search (or the exact quadratic step when `exact_quadratic=Q` is given, to test (1)). `newton_method` solves Hp = −g after a Cholesky check and backtracks. `bfgs` keeps H, computes p = −Hg, takes a Wolfe step, forms s and y, rescales H after the first step and applies `bfgs_update` whenever yᵀs > 0:

```python
def bfgs_update(H, s, y):
    rho = 1.0 / (y @ s)                        # raise ValueError if y.s <= 0
    V = np.eye(s.size) - rho * np.outer(s, y)
    return V @ H @ V.T + rho * np.outer(s, s)  # (3)
```

### Step 3: diagnostics and the cruise model

`error_sequence` and `convergence_order` turn a path into |e_k| and the order estimates log(e_{k+1}/e_k)/log(e_k/e_{k−1}) (`P03.02`), dropping errors at the rounding floor. `fuel_per_km(x, ac)` evaluates §2.7's chain (ISA, V, q, C_L, M_dd, C_D, D, c_T) for x = (M, h/10 km); `cruise_breakdown` returns the intermediate quantities; `optimal_cruise` runs BFGS with the central-difference gradient from (0.70, 0.80), with f = ∞ outside the model's domain so that a long trial step is rejected rather than crashing.

Running `main()` prints

```
Rosenbrock from (-1.2, 1): iterations, f, ||g||, final convergence-order estimates
  steepest descent   14786 iterations, f = 6.32e-13, ||g|| = 9.9e-07, orders 1.00, 1.00, 1.00
  Newton                22 iterations, f = 0.00e+00, ||g|| = 0.0e+00, orders 1.00, 2.95, 1.63
  BFGS                  36 iterations, f = 4.93e-32, ||g|| = 4.4e-16, orders 1.42, 1.19, 0.74
  BFGS error ratios e_{k+1}/e_k over the last 6 steps: [0.601 0.07  0.023 0.011 0.037 0.   ]
  scipy BFGS: 34 iterations, f = 4.35e-25; our final inverse Hessian vs exact H^-1 at x*: max |diff| = 3.42e-05

quadratic, 3 variables, exact line search: BFGS finished in 3 iterations, |H - Q^-1| = 5.6e-17, |x - Q^-1 b| = 3.2e-16
steepest descent on the same quadratic: 189 iterations; worst observed f-ratio 0.7801 against the bound ((k-1)/(k+1))^2 = 0.7805 (kappa = 16.2)
...
minimum-fuel cruise (BFGS with FD gradient, 9 iterations, ||g|| = 5.0e-10):
  M = 0.7760, h = 11955 m, V = 229.0 m/s, CL = 0.585, L/D = 16.14, M_dd = 0.703, drag = 37.2 kN, fuel = 2.2529 kg/km
  scipy BFGS: M = 0.7760, h = 11955 m, fuel = 2.2529 kg/km, 9 iterations
  FD Hessian at the optimum: eigenvalues [ 8.698 86.002] (both positive: a minimum), condition number 9.9
```

followed by the comparison points (2.4012 kg/km at M 0.78, 10 km; 2.3105 at M 0.80, 11 km) and the figure. The order estimates for Newton and BFGS are erratic at the end because the last errors are at the rounding floor; the error ratios are the better diagnostic, and §4 asks for them.

### Verify

- **Analytic cases:** Rosenbrock's gradient vanishes at (1, 1) and its Hessian is positive definite there; the central-difference gradient has observed order 2.00 ± 0.05; the FD Hessian matches the analytic one to 10⁻⁴; a 1-D update gives H = s/y.
- **Limiting cases:** both line searches raise on an ascent direction; the Wolfe search accepts the exact quadratic step unchanged; a step with yᵀs ≤ 0 is refused by `bfgs_update`; below drag divergence the fuel per km does not depend on the wave-drag coefficient at all.
- **Structure:** after each update H satisfies Hy = s to 10⁻¹², is symmetric and positive definite; BFGS with exact line search finishes a 3-variable quadratic in ≤ 3 steps with H = Q⁻¹ to 10⁻¹⁰.
- **Convergence (observed):** steepest-descent ratios ≤ (1) and > 0.3 (linear, not better); Newton's late errors satisfy e_{k+1} < 20 e_k²; BFGS's last ratios average below 0.2 while its first five average above 0.5.
- **Agreement with the oracle:** Rosenbrock's minimiser to 10⁻⁶ and the cruise optimum to 10⁻⁴ in M and 10 m in h against `scipy.optimize.minimize(method='BFGS')`; the cruise optimum beats a 3 × 3 grid of neighbours and has a positive-definite FD Hessian.

### Compare

`scipy.optimize.minimize(method='BFGS')` (`@scipy-userguide` "Broyden-Fletcher-Goldfarb-Shanno algorithm") is the same algorithm with a polynomial-interpolating Wolfe search and a forward-difference gradient by default (`jac=None`), which is why its iteration counts differ slightly (34 against 36 on Rosenbrock, 9 against 9 on the cruise problem). Use it as the production tool; use this lesson's version to see the matrix H, the step lengths and the yᵀs values it hides. `method='Newton-CG'` and the trust-region methods want a Hessian or Hessian-vector product; `method='Nelder-Mead'` wants nothing but f and is the right choice for a noisy or non-smooth objective, at the price of a linear rate and many evaluations.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before the run.

1. **Step size.** For the cruise objective at (0.715, 1.0), tabulate the central- and forward-difference gradient errors for h = 10⁻¹ … 10⁻¹⁰ against a Richardson-extrapolated reference. Predict: central error ∝ h² down to about 10⁻⁵, forward ∝ h down to about 10⁻⁸, then both rise as ε/h.
2. **Method comparison.** On Rosenbrock from (−1.3, 1), count iterations, function evaluations and gradient evaluations for the three methods, and tabulate the last six error ratios of each. Predict a constant ratio (steepest descent), ratios collapsing like e_k itself (Newton) and ratios falling toward zero irregularly (BFGS).
3. **Scaling.** Run BFGS on the cruise objective with x = (M, h in metres) and with x = (M, h/10 km), from the same physical start, with the same tolerance. Predict which converges and estimate the Hessian condition numbers of the two parameterisations.
4. **Cruise-climb.** Optimise for W = 600, 570, 540, 500 kN. Predict from §2.7 that M* and C_L* do not move and h* rises by RT/(g₀) × ln(W₀/W) (325 m from 600 to 570 kN).

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `bfgs_cruise.png`.

What you should find: central-difference errors 2.4 × 10⁻³, 2.5 × 10⁻⁵, 2.5 × 10⁻⁷, 2.5 × 10⁻⁹ in ∂f/∂M for h = 10⁻², 10⁻³, 10⁻⁴, 10⁻⁵, a minimum of 3 × 10⁻¹⁰ at h = 10⁻⁶, then 2 × 10⁻⁷ at 10⁻⁹; on Rosenbrock from (−1.3, 1), 14 546 / 23 / 25 iterations and about 159 000 / 54 / 63 function evaluations for steepest descent, Newton and BFGS; with altitude in metres the Hessian has condition number 1.4 × 10⁹ and BFGS does not meet a gradient tolerance of 10⁻⁷ in 500 iterations although it reaches the same point in f, against 9 iterations with the scaled variables (condition number 10); optima at h* = 11 955, 12 280, 12 623, 13 111 m with M* = 0.7760 and C_L = 0.585 throughout, 2.2529, 2.1402, 2.0276, 1.8774 kg/km.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| steepest-descent f-ratio, κ = 16.2 | ≤ 0.7805 | 0.7801 (worst) | within bound | (1) is sharp for the worst starting direction |
| BFGS steps on a 3-variable quadratic | 3, H = Q⁻¹ | 3, 6 × 10⁻¹⁷ | rounding | finite termination with exact line search |
| inverse Hessian of Rosenbrock at (1, 1) | [[0.5, 1], [1, 2.005]] | BFGS H to 3 × 10⁻⁵ | 3 × 10⁻⁵ | H is learned only along the directions travelled |
| FD gradient error at h = 10⁻³, 10⁻⁴ | ratio 100 | 1.2 × 10⁻⁴ → 1.2 × 10⁻⁶ | 0 | order 2 |
| cruise C_L without wave drag | 0.667, L/D 16.67 | 0.585, L/D 16.14 | −12 % | wave drag pulls C_L down (M_dd falls with C_L) |
| cruise altitude from p = 2W/(γM²SC_L) | 11 955 m | 11 955 m | 0 | consistent ISA |
| cruise-climb slope | −10.6 m/kN | −10.8 m/kN (600 → 570 kN) | logarithm vs. tangent | exact law h ∝ ln(1/W) |
| 1 % window | ±0.023 M, ±715 m | grid in test | — | quadratic model of the optimum |
| iterations, cruise (BFGS / SciPy) | — | 9 / 9 | — | same algorithm, different line search |

Every discrepancy is explained by the derivation, and two deserve a sentence in a report. The optimum C_L is 12 % below the maximum-L/D value and the fuel saving against "M 0.78 at 10 km" is 6.2 %, almost all of it from altitude (at 10 km the aircraft is at C_L = 0.43 and L/D 15.1, and the air is warmer so c_T is higher); the saving against "M 0.80 at 11 km" is 2.5 %, from wave drag. And the optimum is flat: a ±0.02 change in Mach costs under 1 %, which is why airlines fly a cost-index Mach rather than the fuel-optimal one with no loss worth measuring, while an altitude 2 km below the optimum costs several per cent.

What agreement with SciPy does *not* validate: the model. Both codes minimise the same function; a wrong polar gives the same wrong optimum twice. Validation of the cruise point would need the manufacturer's measured drag rise and engine deck, and the note in §6 must say so. Improvements worth making: a cubic-interpolating zoom in the Wolfe search (fewer evaluations per step), L-BFGS for large n (extension 1), and constraints (extension 2), which is where the real operating problem lives. Improvement worth *not* making: tightening `tol` below about 10⁻⁷ with a finite-difference gradient, whose own error is 10⁻⁹ to 10⁻¹⁰ here; the iterations spent beyond that chase noise.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (1), the Newton step, the implication "curvature condition ⇒ yᵀs > 0", the secant and positive-definiteness properties of (3) and the cruise-climb law; `bfgs_cruise.png`; the four tables of §4; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the optimum cruise point, its 1 % window, the saving against two reference points, and the constraint most likely to bind on the real aircraft.

### Adversarial self-checks

- **What if the weight rises 10 %?** M* and C_L* are unchanged (§2.7); h* falls by 6341 × ln 1.1 = 604 m to 11 350 m and fuel per km rises exactly 10 % to 2.478 kg/km, because f ∝ W at fixed (M, C_L). A 10 % rise in c_T0 gives the same fuel and the same point; a 10 % rise in C_D0 moves the optimum to C_L = 0.609, M = 0.773, 2.375 kg/km.
- **What if the drag-divergence model is wrong?** Doubling k_wave moves the optimum only to M = 0.762 and 11 682 m (+1.5 % fuel); raising M_cr0 to 0.84 moves it to M = 0.795, 12 281 m (−2.5 %). The optimum Mach tracks the drag-divergence Mach almost one for one, so it is only as good as that one measured number.
- **Which assumption fails first in the real world?** The absence of constraints. At 600 kN an optimum of 11 955 m is plausible; at 500 kN it is 13 111 m, which is above most narrow-bodies' ceilings. The real answer is the ceiling, and the method of this lesson cannot say so: a constrained optimiser (extension 2) can.
- **What measurement could falsify the model?** Fuel flow at the optimum: the model predicts 0.516 kg/s (1857 kg/h) at 229 m/s. A flight-data value 10 % off with the same weight and conditions means the polar or the engine model is wrong by that much, and the optimum has to be recomputed, not trusted.
- **Does the result survive a change of starting point or tolerance?** From (0.85, 1.3), (0.6, 0.5) and (0.9, 0.9) BFGS reaches the same point in 13, 9 and 12 iterations; with tol = 10⁻³ it stops 3 iterations earlier at a point 0.25 m away in altitude and 7 × 10⁻⁹ kg/km higher. The function has one basin.
- **Which decision changes if the uncertainty doubles?** If the drag-divergence Mach were known to ±0.02 instead of ±0.01, the optimum Mach would be uncertain by about the same amount, which is still inside the 1 % window: the flight-planning decision (M ≈ 0.78, as high as the ceiling allows) would not change. If the wing area or weight were uncertain by 10 %, the altitude would move by 600 m, which does change the flight level requested.

## Check yourself

1. State the first- and second-order conditions for a local minimum and explain why a saddle point satisfies one but not the other.
2. Derive the exact line-search step for steepest descent on a quadratic and state the bound (1). Why does scaling the variables change the rate?
3. Why does Newton's method need a line search and a positive-definiteness check far from the minimum, and why does it need neither near it?
4. Write the strong Wolfe conditions and show that the curvature condition implies yᵀs > 0.
5. Write the BFGS update (3) and prove that it satisfies the secant equation.
6. Why does BFGS on an n-variable quadratic with exact line search terminate in n steps, and what does H equal then?
7. Steepest descent converges linearly, Newton quadratically, BFGS superlinearly. Define each in terms of the error ratios and say how you would measure it from a run.
8. The optimum cruise Mach and C_L of §2.7 do not depend on weight, but the altitude does. Explain both from fuel per km = c_T W/(a M (L/D)), and say what breaks the argument below the tropopause.
