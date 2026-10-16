# Phase 03 Exam: Numerical Methods & Scientific Computing

Integrated assessment of the whole phase. Suggested self-assessment mode: closed reference, a time limit chosen and recorded in advance (the exam has not been piloted, so no duration is suggested), no code except in Part C, no computer algebra in Part B; a scientific calculator (or Python used only as a calculator) is allowed for the arithmetic in Part B, since several questions need logarithms, square roots and a few Newton iterations by hand. There is no enforced lock; the point is to find out what you can do unaided. Solutions and pass thresholds are in `solutions/phase-03/exam/`.

Report your result per competency: (1) mathematical and physical reasoning, (2) analytical problem solving, (3) computational competence, (4) engineering interpretation. Questions are tagged with the competencies they assess. Every question is written over an aerospace situation; the numerical analysis is what is being tested, the framing is there so that you also practise deciding which method applies and what its error model says. The organising idea of the phase is *error*: in every question say which kind (truncation, rounding, stability, convergence order) is in play. Angles are in radians unless a degree sign is written. SI units throughout. Every data set and every convergence table in this exam is invented for the exam and is labelled synthetic; nothing here is a measurement, and no table was produced by a solver you have seen.

## Part A: Concepts

**1.** (1, 4) An orbit-geometry routine evaluates the versine 1 − cos θ for the small angle θ = 10⁻⁴ rad between two position vectors. (a) Machine epsilon in IEEE double precision is 2.2 × 10⁻¹⁶; estimate the *relative* error of computing 1 − cos θ directly, and say how many significant figures survive. (b) Say what the same computation returns in single precision (ε = 1.2 × 10⁻⁷) and why. (c) Name the phenomenon, rewrite the expression in a form that loses nothing, and say why the rewritten form is accurate. (d) Compute the condition number κ = |θ f′(θ)/f(θ)| of f(θ) = 1 − cos θ at small θ and use it to state, in two sentences, the difference between an ill-conditioned *problem* and an unstable *algorithm*. Which is this?

**2.** (1) A re-entry thermal model is integrated with an explicit Runge–Kutta method. (a) Define local truncation error and global error for a one-step method and state how their orders are related; name the two properties whose combination gives convergence (the Dahlquist equivalence theorem). (b) Define the stability function R(z) of a Runge–Kutta method applied to y′ = λy with z = hλ, and say what "the step is inside the stability region" means and what happens to the computed solution when it is not. (c) Say what makes a system *stiff* (in terms of its eigenvalues and the time span of interest) and why an adaptive explicit method such as RK45 does not fail on a stiff problem but becomes uselessly expensive: which limit sets its step, accuracy or stability? (d) An embedded pair gives an error estimate err at step h for a tolerance tol. Write the standard step-size update for a method of order p and say why the exponent is 1/(p + 1), not 1/p.

**3.** (1, 4) A satellite orbit is propagated for a year, first with RK4 at a fixed step, then with the leapfrog (Störmer–Verlet) method at the same step. The RK4 energy falls steadily; the leapfrog energy oscillates with a bounded amplitude and no drift. (a) Explain both behaviours: what property does leapfrog have that RK4 lacks, and what geometric quantity in phase space does it preserve? (b) Explain why the leapfrog energy error is *bounded* but not zero (what quantity *is* conserved to round-off?). (c) RK4 is fourth order and leapfrog second: say why RK4 is nevertheless the wrong tool for a long propagation at fixed step, and when it is the right one. (d) Leapfrog is time-reversible. Say what that means operationally and name one way to test it in code.

**4.** (1) A potential-flow solver discretises Laplace's equation on an N × N interior grid with the five-point stencil and solves the resulting system with the conjugate gradient (CG) method. (a) State the condition CG needs on the matrix and why the five-point Laplacian with Dirichlet boundaries satisfies it. (b) State the two orthogonality properties that define the CG iterates (of the residuals and of the search directions) and the finite-termination property that follows in exact arithmetic; say why nobody relies on finite termination in practice. (c) State the error bound in the A-norm in terms of the condition number κ, and, given that κ ≈ 4(N + 1)²/π² for this matrix, say how the iteration count to a fixed tolerance scales with N and how the cost per iteration scales with N. Evaluate the bound's iteration count for N = 100 and a relative tolerance of 10⁻⁶. (d) Say why Jacobi (diagonal) preconditioning does nothing for this matrix and name a preconditioner that does help.

## Part B: Analysis

Work every derivation by hand and show the steps. Numbers to three significant figures unless more are needed to show a check or a convergence ratio.

**5.** (2) Root finding: Kepler's equation. The eccentric anomaly E of an orbit with eccentricity e = 0.2 at mean anomaly M = 0.5 rad satisfies f(E) = E − e sin E − M = 0.
(a) Write Newton's iteration for this f and perform two iterations from E₀ = M, carrying six decimals. Give f(E₀), f′(E₀), E₁, f(E₁), f′(E₁) and E₂.
(b) The root is E* = 0.615468. Give the errors |E₀ − E*|, |E₁ − E*| and |E₂ − E*|, and show that the ratio |E_{k+1} − E*|/|E_k − E*|² is approximately constant; derive what constant Newton's method predicts, |f″(E*)/(2f′(E*))|, and evaluate it.
(c) How many bisection steps are needed to reduce the bracket [0, 2π] to 10⁻¹⁰? Compare the number of function evaluations with Newton's and state the convergence orders of bisection, secant and Newton.
(d) Say why E₀ = M is a safe start for small e and what can go wrong with Newton from that start when e = 0.99 and M = 0.05 (think about f′ near E = 0). Name a robust alternative start or method.

**6.** (2) Interpolation: an atmospheric table. An exponential atmosphere p(h) = p₀ e^{−h/H} with scale height H = 8500 m is tabulated every Δh = 500 m and read back by interpolation.
(a) Write the error formula for linear interpolation between two nodes and show that its maximum relative error on this table is Δh²/(8H²); evaluate it.
(b) Check (a) directly: evaluate the relative error of linear interpolation at h = 250 m from the nodes h = 0 and 500 m.
(c) What table spacing gives a relative error below 10⁻⁶ with linear interpolation? Write the error bound for a cubic spline, (5/384) h⁴ max |f⁗|, and show that the 500 m table already meets 10⁻⁶ with a spline.
(d) Say why a single degree-40 polynomial through a 41-point table is a bad idea (name the phenomenon and the point distribution that cures it), and say in one sentence what the atmosphere table has that makes a spline the right choice.

**7.** (2) Quadrature: the rocket-equation integral. The velocity gained by a rocket burn is ∫₀^{t_b} ṁv_e/(m₀ − ṁt) dt with m₀ = 250 kg, ṁ = 7.5 kg/s, v_e = 2200 m/s and t_b = 20 s; the exact value is v_e ln 2.5 = 2015.84 m/s.
(a) Show that the two-point Gauss–Legendre rule on [−1, 1] with nodes ±1/√3 and weights 1 integrates every cubic exactly, by requiring exactness on 1, x, x² and x³ (use symmetry). State the degree of exactness of an n-point Gauss rule and of Simpson's rule.
(b) Apply the two-point Gauss rule (mapped to [0, 20]) and Simpson's rule with two intervals to the integral. Give both results and both errors.
(c) The two rules use two and three evaluations. Explain, using the fourth derivative of the integrand f(t) = v_e/(m₀/ṁ − t), why both errors are so large here and where in the interval the error comes from. What would you do about it with the same number of evaluations?
(d) Numerical differentiation. For the central difference (f(x + h) − f(x − h))/(2h), write the two competing error terms, derive the optimal h for a function with f ≈ f‴ ≈ 1, and give the smallest attainable error in double precision. Say why h = 10⁻⁸ (the "small is good" instinct) is worse than h = 10⁻⁵.

**8.** (2) Runge–Kutta order and stability: a pitch-damping mode. The short-period pitch rate obeys q̇ = λq with λ = −8 s⁻¹.
(a) Apply the explicit midpoint method (k₁ = f(t, y), k₂ = f(t + h/2, y + hk₁/2), y_{n+1} = y_n + hk₂) to y′ = λy and show that y_{n+1} = R(z) y_n with R(z) = 1 + z + z²/2, z = hλ. Compare with e^z and deduce the local error order and the global order. Show that the real-axis stability interval is −2 ≤ z ≤ 0.
(b) Write R(z) for classical RK4 and, given that its real-axis stability limit is z = −2.785, give the largest stable step for this mode with forward Euler, with the midpoint method and with RK4.
(c) RK4 on this problem (synthetic run, q(0) = 1) gives errors in q(1 s) of 1.84 × 10⁻⁵ at h = 0.1 s, 8.01 × 10⁻⁷ at h = 0.05 s and 4.23 × 10⁻⁸ at h = 0.025 s. Compute the observed order from each pair, compare with the theoretical order, and explain the direction of the discrepancy. Evaluate R(−0.8) and e^{−0.8} and show the per-step error is of the size −z⁵/120 predicts.
(d) The mode is stable, yet forward Euler with h = 0.3 s produces a growing oscillation. Give the amplification factor and say in one sentence why this is a stability error, not a truncation error, and how the two are distinguished in a convergence study.

**9.** (2) Stiffness and implicit Euler: a two-layer heat-shield model. The temperature excesses of a thin skin and the structure beneath obey ẏ = Ay with

```
A = [ −500.5   499.5 ]  (1/s)
    [  499.5  −500.5 ]
```

(a) Show that the eigenvalues are −1 and −1000 s⁻¹ with eigenvectors (1, 1) and (1, −1), and give the stiffness ratio. Say which mode you care about after the first few milliseconds.
(b) Give the largest stable step for forward Euler and the number of steps to reach t = 10 s. Then write the implicit Euler step (I − hA) y_{n+1} = y_n, show that its amplification factor on each mode is 1/(1 − hλ), and show that |1/(1 − hλ)| < 1 for every h > 0 when λ < 0 (A-stability).
(c) With implicit Euler at h = 0.1 s, give the amplification factors of both modes per step, the slow-mode value after 100 steps and the exact e^{−10}; give the relative error and its order. Repeat the slow-mode comparison at h = 0.01 s and confirm the error falls by the expected factor.
(d) Say what each implicit step costs that an explicit step does not (for a linear system; for a nonlinear one), and why that cost is still worth it here. An adaptive RK45 run on this system to t = 10 s (synthetic run) took 3045 steps of mean length 0.0033 s: say what set that step length.

**10.** (2) Symplectic integration: a circular orbit. In normalised units a circular orbit obeys ẍ = −ω²x, ÿ = −ω²y with ω = 1; consider one coordinate. The Störmer–Verlet (leapfrog) step is x_{n+1} = x_n + hv_n + (h²/2)a_n, v_{n+1} = v_n + (h/2)(a_n + a_{n+1}) with a = −ω²x.
(a) Write one step as (x_{n+1}, v_{n+1})ᵀ = M (x_n, v_n)ᵀ and give the 2 × 2 matrix M.
(b) Show that det M = 1 and give tr M. Say what det M = 1 means for areas in the (x, v) plane and why that is the discrete version of what the exact flow does.
(c) Show that the method is stable (|eigenvalues| = 1) exactly when hω < 2, and that the eigenvalues are e^{±iω̃h} with cos(ω̃h) = 1 − h²ω²/2. For h = 0.1, evaluate ω̃/ω − 1 and the phase error per orbit; show the leading term is h²ω²/24.
(d) The energy E = (v² + ω²x²)/2 is not exactly conserved, but Ẽ = (v² + ω²x²(1 − h²ω²/4))/2 is. Verify that Ẽ is unchanged by one step (use M). Deduce the maximum relative energy error over any number of orbits for h = 0.1 and say why it does not grow.

**11.** (2) Direct solvers: a three-node structural model. The stiffness equations of a three-degree-of-freedom model are Kx = b with

```
K = [ 4  2   2 ]        b = [ 10 ]
    [ 2  5   3 ]            [ 13 ]
    [ 2  3  11 ]            [ 27 ]
```

(a) Compute the Cholesky factor L (K = LLᵀ) by hand, writing the formula for each entry as you go. Why is K guaranteed to have one (what property of a stiffness matrix with the supports applied)?
(b) Solve Ly = b by forward substitution and Lᵀx = y by back substitution; check Kx = b.
(c) Give the leading-order flop counts of LU with partial pivoting and of Cholesky for an n × n matrix, and of one triangular solve. For n = 10⁴ evaluate them and say what the ratio is. The model is then solved for 500 load cases: say how the work divides between factorisation and solves, and why you never form K⁻¹.
(d) Say why partial pivoting is unnecessary for an SPD matrix and what the Cholesky algorithm does when the matrix is *not* positive definite (what goes wrong, and what you learn from it about the structure).

**12.** (2) Nonlinear systems: trimming an aircraft. Longitudinal trim of a light aircraft is the pair of equations F(α, V) = 0,

F₁ = ½ρV²S(C_{L0} + C_{Lα}α) − W,  F₂ = C_{m0} + C_{mα}α + C_{mδ}δ_e,

with ρ = 1.225 kg/m³, S = 16 m², W = 12 000 N, C_{L0} = 0.3, C_{Lα} = 5.0 /rad, C_{m0} = 0.05, C_{mα} = −0.8 /rad, C_{mδ} = −1.2 /rad and the elevator fixed at δ_e = 0.
(a) Write the Jacobian J(α, V) analytically.
(b) From (α₀, V₀) = (0.05 rad, 40 m/s), evaluate F and J, solve the 2 × 2 Newton system J Δz = −F by hand, and give (α₁, V₁).
(c) Evaluate F at (α₁, V₁). One component is now exactly zero and stays so: say why. Perform the second Newton step and give V₂; the exact trim speed is 44.712 m/s. Show the V-errors of the three iterates behave as quadratic convergence predicts.
(d) If the Jacobian were computed by forward finite differences, say how to choose the step (in terms of machine epsilon) and why, and state the cost of each Newton iteration in function evaluations for n unknowns. Say what safeguard you add when the starting guess is far from the trim point.

**13.** (2) Finite differences: FTCS for the re-entry tile. The slab u_t = αu_xx with α = 5 × 10⁻⁷ m²/s, L = 0.05 m and u = 0 at both faces is discretised with N = 50 cells (Δx = 1 mm) and the forward-time centred-space scheme u_j^{n+1} = u_j^n + r(u_{j+1}^n − 2u_j^n + u_{j−1}^n), r = αΔt/Δx².
(a) Substitute the von Neumann mode u_j^n = Gⁿe^{ikjΔx} and derive G(k) = 1 − 4r sin²(kΔx/2). Deduce the stability condition r ≤ ½ and say which wavelength goes unstable first.
(b) Give the largest stable Δt and the number of steps to reach t = 1167 s (the time at which the fundamental mode has decayed to 10% of its initial value, from the Phase 02 exam). If Δx is halved, by what factor do the step count and the total work change?
(c) State the truncation error of FTCS in Δt and Δx. At r = ½ the fundamental mode (k = π/L) is damped by G per step: evaluate G and Gⁿ after the 1167 steps and compare with e^{−t/τ₁}, τ₁ = L²/(απ²) = 506.6 s. Give the relative error.
(d) Write the BTCS and Crank–Nicolson amplification factors, say which is unconditionally stable and which is second order in time, and say what each costs per step compared with FTCS. Which would you choose to run this slab to steady state, and why?

**14.** (2) Finite volumes: an acoustic pulse in a duct. The linear advection equation u_t + au_x = 0 with a = 340 m/s is solved on cells of width Δx = 0.01 m with ν = aΔt/Δx = 0.8.
(a) Write the cell-average update in conservation form u_j^{n+1} = u_j^n − (Δt/Δx)(F_{j+½} − F_{j−½}) and give the upwind flux. Show in two lines that the scheme conserves Σ_j u_j exactly on a periodic domain.
(b) Derive the upwind amplification factor G = 1 − ν + νe^{−iθ} (θ = kΔx) and show |G|² = 1 − 2ν(1 − ν)(1 − cos θ), hence stability for 0 ≤ ν ≤ 1 (the CFL condition). Give Δt here.
(c) Derive the amplification factor of the central (FTCS) scheme and show it is unstable for every ν > 0. Evaluate |G| for both schemes at θ = π/2 (the four-cell wavelength) with ν = 0.8.
(d) By Taylor-expanding the upwind scheme, show it solves u_t + au_x = Du_xx + ... with D = aΔx(1 − ν)/2, and evaluate D. Over the 1/340 s the pulse takes to travel 1 m, how far does numerical diffusion spread a sharp front (take √(2Dt))? Say why ν = 1 would be perfect and why it cannot be used in practice.

**15.** (2) Finite elements: a strut under launch load. An aluminium strut (E = 70 GPa, A = 2 × 10⁻⁴ m², L = 2 m) is fixed at x = 0 and carries an axial tip load P = 20 kN. Model it with two equal linear elements.
(a) Starting from the weak form ∫ EA u′ w′ dx = P w(L), derive the element stiffness matrix k_e = (EA/L_e)[[1, −1], [−1, 1]] using the linear shape functions.
(b) Assemble the 3 × 3 global stiffness matrix. Show that it is singular, give the null vector and say what it means physically. Apply the boundary condition and give the reduced system.
(c) Solve for the nodal displacements and compare with the exact u(x) = Px/(EA). Explain why the nodal values are exact here and what the error in the stress is.
(d) Now add the inertial body load of a launch acceleration a = 50 m/s² (ρ = 2700 kg/m³, so q = ρaA = 27 N/m, synthetic case; P removed). Give the consistent nodal loads, solve, and compare with the exact u(x) = (q/(EA))(Lx − x²/2) at the nodes. State the convergence rates in the energy norm and in the displacement for linear elements as h → 0.

**16.** (2) Unconstrained optimisation: a two-variable drag model. A scaled wing-design drag model is f(x) = ½xᵀHx with H = diag(2, 10) (synthetic; minimum at the origin), starting at x₀ = (5, 1).
(a) Give the gradient g₀, the exact line-search step α₀ = g₀ᵀg₀/(g₀ᵀHg₀) for steepest descent, x₁ and g₁. Give f(x₀) and f(x₁). Show that the Armijo condition f(x₁) ≤ f(x₀) + c₁α₀g₀ᵀp₀ holds for c₁ = 10⁻⁴.
(b) With s₀ = x₁ − x₀ and y₀ = g₁ − g₀, apply the BFGS update to the inverse Hessian approximation H₀ = I and give H₁. Verify the secant condition H₁y₀ = s₀. Why must s₀ᵀy₀ > 0, and what guarantees it here?
(c) Compute the BFGS direction p₁ = −H₁g₁ and the exact line-search step along it. Where does x₂ land, and why is this not a coincidence (which theorem)?
(d) Steepest descent with exact line search on this f reduces f by the factor ((κ − 1)/(κ + 1))² per iteration with κ = 5: evaluate it and the number of iterations for a 10⁶ reduction. Say what this teaches about scaling the design variables before calling `scipy.optimize.minimize`.

**17.** (2) Constrained optimisation: a propellant tank in a fairing. A cylindrical tank of radius r and height h must hold volume V = 2 m³ with the least wall area (mass) f = 2πr² + 2πrh.
(a) Write the problem as min f subject to g(r, h) = V − πr²h ≤ 0, form the Lagrangian with multiplier λ, and write the KKT conditions (stationarity, primal feasibility, dual feasibility, complementary slackness).
(b) Solve them: show h = 2r, give r, h, f* and λ, and check complementary slackness.
(c) Show that λ = df*/dV (the shadow price), by differentiating f*(V) = 6π(V/2π)^{2/3}. Give its units and its engineering meaning.
(d) The fairing limits r ≤ 0.6 m. Add the constraint c = r − 0.6 ≤ 0 with multiplier μ, re-solve the KKT system, and give r, h, f*, λ and μ. Say what the sign of μ tells you, what the fairing costs in area, and in one sentence what the SQP method does with these conditions (what it solves at each iteration).

**18.** (2) Multidisciplinary analysis: coupled derivatives. A synthetic two-discipline wing model couples a lift variable y₁ and a tip-deflection variable y₂ through a span parameter x: aerodynamics y₁ = 2x + 0.3y₂, structure y₂ = 0.5y₁ + x, and the objective is f = y₁ − 0.2y₂².
(a) Solve the coupled system at x = 1 exactly. Then run Gauss–Seidel from (0, 0) (aero first, then structure) for three iterations and give the errors in y₁; show that they fall by the factor 0.15 per iteration and say where 0.15 comes from. How many iterations to 10⁻⁸?
(b) Write the residual form R(x, y) = 0 and derive the direct method for the total derivative: (∂R/∂y) dy/dx = −∂R/∂x. Solve the 2 × 2 system for dy₁/dx and dy₂/dx at x = 1.
(c) Assemble df/dx = ∂f/∂y · dy/dx and evaluate it. Then derive the adjoint method: solve (∂R/∂y)ᵀψ = (∂f/∂y)ᵀ and form df/dx = −ψᵀ∂R/∂x; evaluate ψ and show the same number results.
(d) Say when direct is cheaper and when adjoint is (in terms of the numbers of design variables and objectives), and why a finite-difference total derivative is both expensive and inaccurate here (what must be re-converged inside each perturbation, and how its tolerance interacts with the difference step).

**19.** (2) Uncertainty quantification: Δv of a stage. The burnout velocity gain Δv = v_e ln(m₀/m_f) has m₀ = 250 kg fixed, v_e ~ N(2200, 30²) m/s and m_f ~ N(100, 2²) kg, independent.
(a) Linearise and give the two contributions to σ_Δv and their root-sum-square. Give the nominal Δv.
(b) A Monte Carlo estimate of the mean has standard error σ/√N: give N for a standard error of 1 m/s. The relative standard error of the *standard deviation* estimate is about 1/√(2N): give N for 1%. Say how the error of the 1st-percentile estimate compares and why.
(c) Hold v_e at 2200 and expand only in ξ = (m_f − 100)/2 ~ N(0, 1). Taylor-expand ln(100 + 2ξ) to second order and rewrite the result in probabilists' Hermite polynomials He₀ = 1, He₁ = ξ, He₂ = ξ² − 1 to get the polynomial-chaos coefficients c₀, c₁, c₂ of Δv. Use the orthogonality ⟨He_m He_n⟩ = n! δ_mn to give the mean and the standard deviation from the coefficients. Compare with (a) and say where the mean shift comes from.
(d) State how the polynomial-chaos error decreases with the degree for a smooth response and how the Monte Carlo error decreases with N, give the number of terms of a total-degree-3 expansion in 2 and in 10 random variables, and name the two situations in which Monte Carlo is the method to use.

**20.** (2) Grid convergence and Richardson extrapolation. A skin-friction coefficient from a boundary-layer solver (second-order scheme, synthetic run) on three nested grids with spacing h, h/2 and h/4 is C_f = 0.030928, 0.030216 and 0.030052.
(a) Derive the observed order of accuracy p from the three values assuming C_f(h) = C* + Ch^p, and evaluate it.
(b) Derive the Richardson extrapolate C* from the two finest grids and evaluate it with the observed p. Repeat with the theoretical p = 2 and compare the two.
(c) The exact value (known to the examiner) is C* = 0.030000. Give the error of each raw grid value and of each extrapolate. Say why the observed order is above 2 and what that tells you about the asymptotic range.
(d) Define the grid convergence index GCI = 1.25|ε|/(r^p − 1) with ε the relative change between the two finest grids, evaluate it, and say what it is an estimate of and what it is *not* (what Richardson extrapolation can never detect). What is the difference between code verification and solution verification, and which one is this?

## Part C: Computation

**21.** (3) Write a Python function `conjugate_gradient(A, b, tol, maxiter)` that solves Ax = b for a symmetric positive definite A from x₀ = 0 and returns the solution and the history of residual norms, and a helper `poisson_1d(n)` that returns the n × n matrix of −u″ with u = 0 at both ends (scaled by h²). Write three pytest tests: one that solves −u″ = 1 on (0, 1) and compares with the closed form u = x(1 − x)/2, which the second-difference scheme reproduces exactly at the nodes; one that checks finite termination (no more than n iterations to 10⁻¹² on a random SPD matrix) and agreement with `scipy.sparse.linalg.cg` or `numpy.linalg.solve`; and one that checks that the error decreases monotonically in the A-norm ‖e‖_A = √(eᵀAe) from one iteration to the next. Then explain in two sentences why the third test uses the A-norm and not the residual 2-norm, and in one sentence why `scipy.sparse.linalg.cg` with a preconditioner is still what you would call in production.

## Part D: Interpretation

**22.** (4) A learner's lesson 05 RK4 integrator passes its tests, including an observed-order test on y′ = −8y that reports 3.98. Applied to the forced mode y′ = −y + sin t, y(0) = 0, on [0, 2] and compared with the exact solution, it gives the following convergence table (synthetic):

| h | error in y(2) |
|---|---|
| 0.2 | 1.45 × 10⁻² |
| 0.1 | 6.59 × 10⁻³ |
| 0.05 | 3.14 × 10⁻³ |
| 0.025 | 1.53 × 10⁻³ |
| 0.0125 | 7.55 × 10⁻⁴ |

(a) Compute the observed order from each successive pair. Is the method converging? To what order? Does it converge to the right answer?
(b) List at least five distinct hypotheses that produce a lower observed order than the theory promises (think of: the implementation, the problem, the test, the arithmetic), and for each give the one-line diagnostic that would confirm or rule it out.
(c) Here the autonomous test problem reports order 4 and the non-autonomous one reports order 1. Say which hypothesis this pattern points to and what single change to the test suite would have caught it.
(d) At the other end: the correct RK4 on this problem has error 1.4 × 10⁻¹⁵ at h = 2 × 10⁻⁴ but 2.1 × 10⁻¹² at h = 2 × 10⁻⁶. Explain the reversal and give the error model that predicts where the floor is.
(e) Can Richardson extrapolation be applied to the table above, and with which p? Write the two-sentence statement about this solver that belongs in a verification report, and say what *validation* would require that no amount of this work provides.
