# P03.16  Multidisciplinary Design Optimization

> A wing is not an aerodynamic object plus a structural object: its weight depends on its loads and its loads depend on its weight. Converge that loop, differentiate through it, and only then hand the design to an optimiser.

| | |
|---|---|
| **Learning outcomes** | 1. Write a two-discipline wing model as a coupled residual R(x, y) = 0 between design variables x and coupling variables y, derive the Gauss–Seidel fixed-point iteration and its linear convergence rate from the coupling strength, and the Newton multidisciplinary analysis with its quadratic rate. 2. Derive the total derivative of a function of interest through the coupling by the direct (forward) method and by the coupled adjoint, and state when each is cheaper. 3. Implement and test, from scratch in NumPy, the disciplines with analytic partials, both MDA solvers, the three derivative routes (finite difference, direct, adjoint) with an observed second-order finite-difference convergence, and the MDF and IDF optimisation architectures with `scipy.optimize.minimize` (SLSQP) as the driver and `scipy.optimize.root` as the MDA oracle. 4. Minimise block fuel over span and area, read which constraint sets the span, and judge where the model stops being believable. |
| **Prerequisites** | `P03.15` (KKT, SLSQP, multipliers as prices), `P03.10` (Newton for nonlinear systems, fixed-point contraction), `P03.14` (gradient-based search, why gradient accuracy matters), `P03.08` (small linear solves and their transposes), `P01.09` (chain rule in several variables, implicit differentiation), `P01.13` (block matrices). Skill check: you can apply the implicit function theorem to R(x, y) = 0 to get dy/dx = −(∂R/∂y)⁻¹ ∂R/∂x, and you can write the Breguet range equation. |
| **Engineering question** | A light twin of 25 kN empty weight (less wing) must fly 1500 km at 80 m/s at 3000 m. Its wing's span and area are free; its spar caps are sized to the ultimate root bending moment, which depends on the take-off weight, which includes the fuel, which depends on the drag, which depends on the span and on the weight of the spar caps. What span and area burn the least fuel, which limit stops the span from growing, and how much would a stiffer spar material be worth in kilograms of fuel? |
| **Assumptions and limits** | Aerodynamics: a parabolic drag polar with fixed parasite area and Oswald factor, elliptic lift distribution, lift-to-drag ratio evaluated at the zero-fuel weight, Breguet cruise with constant specific fuel consumption. Structures: a two-cap spar of depth 12 % chord sized to a fatigue allowable at every section for a moment distribution M_root(1 − 2|y|/b)², skins and systems at a fixed areal weight, a fully stressed beam (so the tip deflection is load-independent), no flutter, no torsion, no fuel in the wing. All weights are the author's invented but realistic numbers, not a specific aircraft. The optimiser assumes smooth functions and a regular active set. |
| **Analytical work** | The two disciplines and their analytic partials; the coupled residual and its Jacobian; the Gauss–Seidel contraction factor k = (∂M_root/∂W_wing)(∂W_wing/∂M_root) and the Newton step; the implicit-function derivative dy/dx; the direct and adjoint total-derivative formulas and their costs; the fully stressed beam's constant curvature and the resulting aspect-ratio limit AR ≤ 3Eτδ_max/σ_a = 12.6; the first Gauss–Seidel sweep by hand. |
| **Computational deliverable** | `starter/main.py`: `aero_mission`, `aero_mission_partials`, `structures`, `structures_partials`, `coupled_residual`, `coupled_jacobians`, `mda_gauss_seidel`, `mda_newton`, `coupling_strength`, `functions_of_interest`, `functions_partials`, `total_derivatives_direct`, `total_derivatives_adjoint`, `total_derivatives_fd`, `fd_convergence_study`, `design_constraints_x_only_jac` (the dataclasses, `design_constraints_x_only`, `mdf_optimize`, `idf_optimize` and `main` are given). Tests prove: each discipline reproduces its formulas and limits (zero range burns no fuel, doubling the span removes three quarters of the induced drag, a zero moment leaves the skin weight, the deflection ratio is load-free); every analytic partial matches central differences to 10⁻⁶; Gauss–Seidel converges with an observed ratio equal to k to 3 %; Newton converges with observed order above 1.7, agrees with Gauss–Seidel and with `scipy.optimize.root` to 10⁻⁹; direct and adjoint totals agree to 10⁻¹⁰ and with finite differences to 10⁻⁵; the adjoint variables satisfy their linear system; the finite-difference error falls as h² (slope 2.0 ± 0.1); the MDF optimum is feasible, 10 % better than the baseline, sits exactly on the aspect-ratio limit, is reached from two starts, equals the IDF optimum to 10⁻⁵, and the FD-gradient run costs more MDA solves; a stiffer material buys span and fuel, a longer range costs fuel. |
| **Practical activity** | **computational experiment**: the fuel-against-span curve at fixed area with and without the deflection limit; the coupling strength and Gauss–Seidel behaviour as the span grows until the fixed point stops being a contraction; the finite-difference step study; MDA-solve counts for adjoint against finite-difference gradients; the optimum against range and spar stiffness. |
| **Equipment** | None. A computer with the environment from `P00.01`. OpenMDAO (part 6, extension) is optional and separately installed; nothing here requires it. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `wing_mdo.png`, the hand derivation of the partials, of k and of the adjoint equations, the first Gauss–Seidel sweep by hand, the five experiment tables, and a one-page design note stating the optimum, the active constraint, the value of the deflection allowance in kilograms of fuel, and the span beyond which you would not trust the model. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the adjoint equations are derived and the direct/adjoint cost comparison stated; all public tests pass; the coupling strength and the FD order are measured and explained; the design note names the active constraint and prices it. |
| **References** | `@ocw-16-90` Unit 3 "3.6 Introduction to Design Optimization" (required-free: design variables, objective, constraints, gradient-based optimisation). `@ocw-18-335j` Week 9 "Lecture 26: Adjoint Methods" (required-free: the adjoint trick, one solve per output instead of one per input). `@scipy-userguide` "Optimization (scipy.optimize)" > "Constrained minimization": "Sequential Least SQuares Programming (SLSQP) Algorithm" (required-free: the driver); "Root finding" > "Sets of equations" (supplementary-free: `root` as the MDA oracle). `@openmdao-docs` "Getting Started"; "Basic User Guide": "Single-Discipline Optimization", "Multidisciplinary Optimization", "Using the OpenMDAO Command Line Utility", "Recording and Reading Data" (supplementary-free, optional: the Sellar problem, Groups and Components, the N2 diagram, total and partial derivatives; read it after part 6 if you install the package). |
| **Extensions** | ★ Rebuild the model as OpenMDAO components with `declare_partials`, let its Newton solver and coupled adjoint reproduce this lesson's numbers, and draw the N2 diagram (optional; requires `pip install openmdao`, not part of the tested code). ★ Add wing fuel (fuel in the wing relieves the root moment, a third coupling) and see the coupling strength change sign. ★ Replace the fatigue allowable by a flutter-speed constraint that scales with torsional stiffness and recover a different active constraint. |
| **Time** | unmeasured |

## 1. The engineering problem

Aircraft design is coupled. The aerodynamicist wants span, because induced drag falls as 1/b²; the structural engineer pays for it, because the root bending moment grows with b and the spar caps that carry it grow with b² and more; and the mission analyst closes the loop, because every newton of structure must be lifted, which costs drag, which costs fuel, which must be lifted too. No single discipline can answer "what span?", and the conventional way of designing (freeze the weight, size the wing, resize the structure, repeat by hand) is a fixed-point iteration carried out by meetings. **Multidisciplinary design optimization** (MDO) makes the loop explicit, converges it numerically, differentiates through it, and lets an optimiser move the design.

This lesson builds the smallest version with all of the structure of the real thing: two disciplines (aerodynamics-plus-mission and structures), two design variables (span b and area S), two coupling variables (the ultimate root bending moment M_root the aerodynamics sends to the structure, and the wing weight W_wing the structure sends back), one objective (block fuel for a 1500 km mission), and three constraints (stall speed, tip deflection, hangar span). The aircraft is an invented light twin: 25 kN of airframe, engines, payload and systems other than the wing and the fuel, cruising at 80 m/s at 3000 m. Everything is NumPy; `scipy.optimize.minimize` drives the design and `scipy.optimize.root` checks the coupled analysis. Frameworks such as OpenMDAO exist to automate the bookkeeping for dozens of disciplines (part 6); you should write it once by hand to know what they do.

Success criteria, written first: the coupled analysis converges to a residual of 10⁻¹⁰ relative; its Gauss–Seidel rate equals the hand-derived coupling strength to 3 %; the total derivative of fuel with respect to span computed three independent ways (finite differences through the converged loop, the direct method, the coupled adjoint) agrees to 10⁻⁵, and the finite-difference error falls as the square of the step; the MDF and IDF architectures find the same optimum to 10⁻⁵; and the optimum is explained by one active constraint with a stated price in kilograms of fuel.

## 2. Understand and derive

### 2.1 Discipline 1: aerodynamics and mission

Inputs: b, S and the wing weight W_wing. In cruise, lift equals weight; we evaluate the lift-to-drag ratio at the **zero-fuel weight** W_zf = W₀ + W_wing (a stated simplification: the end-of-cruise condition), so L = W_zf. The drag in dimensional form, with q = ½ρV² = 2909 Pa,

**D = q (f_rest + c_dw S) + L²/(q π e b²),**   (1)

has a parasite part (an equivalent flat-plate area f_rest = 0.25 m² for fuselage, nacelles and tail, plus the wing's profile drag c_dw S with c_dw = 0.008) and the induced drag of an elliptically loaded wing with Oswald factor e = 0.85. At the baseline (b = 12 m, S = 16.8 m²) the converged loop gives D = 1118 + 683 = 1802 N and L/D = 15.35 at C_L = 0.566. **Breguet** with thrust-specific fuel consumption c_T = 1.4 × 10⁻⁵ kg/(N s) over the range R gives the block fuel

**W_fuel = W_zf (e^{κ D/L} − 1), κ = R g c_T/V = 2.574,**   (2)

the take-off weight W_TO = W_zf + W_fuel, and the **ultimate root bending moment** of one wing under elliptic loading at the ultimate load factor n = 2.5 × 1.5 = 3.75 (the lift centroid of a half-ellipse is at 2b/(3π) from the root),

**M_root = n W_TO b/(3π).**   (3)

Outputs: W_fuel (the objective), M_root (to the structure), and W_TO (for the stall constraint, S ≥ W_TO/(q_stall C_Lmax) with a 45 m/s sea-level stall speed and C_Lmax = 1.8: 14.65 m² at the baseline).

### 2.2 Discipline 2: structures

Inputs: b, S and M_root. The spar box has depth h = τc = τS/b with τ = 0.12 (0.168 m at the baseline). Two caps of area A carry M at the allowable stress σ_a = 120 MPa (a fatigue allowable for aluminium, not its yield strength): σ = M/(Ah), so A = M_root/(hσ_a), 77.4 cm² at the baseline. Size every section to its own moment, M(y) = M_root(1 − 2|y|/b)² (the moment of a triangular-ish loading, a stated approximation), integrate 2A(y) over both halves, and the cap weight is ρ_m g (2/3) A_root b. Skins, ribs, controls and systems add w_skin S with w_skin = 60 N/m²:

**W_wing = (2/3) ρ_m g b² M_root/(τ S σ_a) + w_skin S.**   (4)

The b²/S makes span expensive: at fixed area a longer wing has a shallower spar *and* a longer lever arm. A fully stressed beam has the same curvature everywhere, 2σ_a/(Eh) at ultimate load, so the tip deflection at limit load (one 1.5th of ultimate) is δ = σ_a (b/2)²/(1.5 E h) and the **deflection ratio** δ/(b/2) = σ_a b²/(3Eτ S) does not depend on the load at all: it is a function of the aspect ratio b²/S alone, 0.0408 at the baseline. The deflection limit δ_max = 0.06 therefore becomes an aspect-ratio limit, **AR ≤ 3Eτδ_max/σ_a = 12.6**, and it is this constraint, not stress, that will stop the optimiser.

### 2.3 The coupled system

Collect the coupling variables y = (y₁, y₂) = (M_root, W_wing). Discipline 1 is a map G₁(x, y₂) → y₁; discipline 2 is G₂(x, y₁) → y₂. A consistent design satisfies the **coupled residual**

**R(x, y) = y − G(x, y) = 0.**   (5)

Solving (5) for y at fixed x is the **multidisciplinary analysis** (MDA). Two ways:

**Block Gauss–Seidel.** Start with W_wing = 0, run discipline 1 (M_root = 141.86 kN m), run discipline 2 with that moment (W_wing = 2498.6 N), and repeat. Near the solution the error e in y₂ obeys e⁺ = k e with

**k = (∂G₁/∂y₂)(∂G₂/∂y₁) = (∂M_root/∂W_wing)(∂W_wing/∂M_root),**   (6)

the **coupling strength**: the product of the two off-diagonal entries of ∂G/∂y, which is also the spectral radius of the Gauss–Seidel iteration matrix for a 2 × 2 system. The iteration converges linearly at rate |k| when |k| < 1 and diverges otherwise. At the baseline ∂M_root/∂W_wing = 5.418 N m/N and ∂W_wing/∂M_root = 0.01051 N/(N m), so k = 0.0569, and the residuals of successive sweeps are 1.34 × 10⁴, 765, 43.6, 2.48, 0.141, 8.0 × 10⁻³, … : each is 0.0569 times the last, as the test demands to 3 %. Nine sweeps reach 10⁻¹⁰ relative. The loop converges to M_root = 156.12 kN m and W_wing = 2648.4 N (270 kg, 8 % of the take-off weight), W_fuel = 5049 N (515 kg), W_TO = 32 698 N.

**Newton.** Linearise (5): (∂R/∂y) Δy = −R with ∂R/∂y = I − ∂G/∂y = [[1, −5.418], [−0.01051, 1]], solve the 2 × 2 system, update. Two iterations from a rough guess take the residual from 5.3 × 10⁴ to 5.5 to 5.0 × 10⁻⁸: order ln(5.0 × 10⁻⁸/5.5)/ln(5.5/5.3 × 10⁴) = 2.0, the quadratic convergence of `P03.10`. Newton needs the coupling derivatives; Gauss–Seidel does not, which is why loosely coupled, derivative-free disciplines (a commercial CFD code and a commercial FE code) are usually run Gauss–Seidel with relaxation, and tightly coupled ones with Newton. With k = 0.06 the choice hardly matters here; part 4 pushes the span until k passes 1 and Gauss–Seidel fails while Newton does not.

### 2.4 Derivatives through the coupling

The optimiser needs dF/dx for the functions of interest F = (W_fuel, g_stall). F depends on x directly and through y(x), and y(x) is defined implicitly by (5). Differentiate R(x, y(x)) = 0 (the implicit function theorem, `P01.09`):

**(∂R/∂y)(dy/dx) = −∂R/∂x, dF/dx = ∂F/∂x + (∂F/∂y)(dy/dx).**   (7)

Here ∂F/∂x and ∂F/∂y are **partial** derivatives with the other argument frozen, cheap and local to one discipline, while dy/dx is the **total** response of the converged loop. Two ways to organise (7):

- **Direct (forward) method:** solve (∂R/∂y) dy/dx = −∂R/∂x, one linear solve per design variable (n_x right-hand sides), then multiply. At the baseline dy/dx = [[13 610, 142.4], [416.4, −36.1]]: a metre of span adds 13.6 kN m of root moment and 416 N of wing weight *after* the loop has re-converged.
- **Adjoint method** (`@ocw-18-335j` Lecture 26): for each function F_i solve the transposed system (∂R/∂y)ᵀ ψ_i = (∂F_i/∂y)ᵀ, then

**dF_i/dx = ∂F_i/∂x − ψ_iᵀ ∂R/∂x.**   (8)

One solve per *function* (n_F of them) instead of one per *variable*. For a wing with 200 thickness and twist variables and 3 functions, the adjoint is 70 times cheaper; for this lesson's 2 × 2 the two are the same work and, as the test checks, give identical numbers to 10⁻¹⁰ because they are the same linear algebra transposed. The fuel adjoint at the baseline is ψ = (0.00150, 0.1429): the second entry says a newton of residual in the wing-weight equation is worth 0.143 N of fuel, which is also dW_fuel/dW_wing through the loop.

The result: **dW_fuel/db = −290.65 N/m** and dW_fuel/dS = +65.97 N/m² at the baseline. Compare the partial with the coupling frozen, ∂W_fuel/∂b = −346.77 N/m: ignoring the weight feedback overstates the benefit of span by 19 %. A designer who optimises the wing with a frozen structural weight is pulled toward too much span; this gap is the whole reason MDO exists.

**Finite differences** are the third route: converge the MDA at x ± h and difference. They cost two MDAs per variable and carry a truncation error O(h²) plus the MDA tolerance divided by h. With the MDA converged to 10⁻¹³ the errors against the adjoint are 6.78, 6.66 × 10⁻², 6.66 × 10⁻⁴, 6.66 × 10⁻⁶ N/m for relative steps 10⁻¹, 10⁻², 10⁻³, 10⁻⁴: slope 2.00, down to a rounding floor near 3 × 10⁻⁸ N/m. Replace the Newton MDA by Gauss–Seidel stopped at a relative residual of 10⁻⁶ and the error stops falling at 2.4 × 10⁻⁴ N/m (10⁻³ gives 0.04 N/m): a loosely converged loop puts a floor on every finite-difference gradient at about the loop's tolerance times the derivative, and hides the h² law below h = 10⁻³. That floor is the standard way an MDO run goes wrong quietly.

### 2.5 Two architectures

**MDF (multidisciplinary feasible).** The optimiser sees only x; every function evaluation runs a full MDA, and the gradient comes from (8). Every iterate is a consistent design, so the run can be stopped at any time with a usable answer. Cost: an MDA per evaluation, which is why the MDA must be fast and the derivatives must not be finite differences.

**IDF (individual discipline feasible).** Promote y to design variables, z = (x, y), evaluate each discipline once per function call with *no* inner loop, and add the consistency equations R(x, y) = 0 as equality constraints. The disciplines are consistent only at the optimum; in between, the optimiser is free to explore inconsistent designs, which is what lets it run the disciplines in parallel. Both architectures must find the same optimum, and the test checks that they do to 10⁻⁵; the `mdf_optimize` and `idf_optimize` wrappers are given, scaled (objective in kilonewtons, coupling variables by 10⁵ N m and 10³ N) exactly as `P03.15` taught.

### 2.6 The optimum

From the baseline, SLSQP with adjoint gradients takes 12 iterations and 49 MDA solves to **b = 14.869 m, S = 17.546 m², AR = 12.60, W_fuel = 4511.9 N = 460.1 kg**, a saving of 54.8 kg (10.6 %) over the baseline's 515 kg. The wing weight grows from 270 to 421 kg and the cruise L/D from 15.35 to 17.87. The deflection constraint is active to 10⁻¹² (AR = 12.60 exactly), the stall constraint has 14 % margin and the gate limit 26 %. With SLSQP's own finite-difference gradients the same point costs 85 MDA solves; IDF takes 13 iterations with a final consistency residual of 10⁻¹². Note what the optimiser did with the area: it *grew* S from 16.8 to 17.5 m², not to reduce wing loading (stall is slack) but because at the aspect-ratio limit more area allows more span, b = √(12.6 S), and the induced-drag gain of span outweighs the profile drag and skin weight of area. Only a coupled, constrained view produces that answer.

The multiplier of the deflection constraint, by least squares on the KKT condition of `P03.15`, is 1024 N of fuel per unit of the normalised constraint: allowing 10 % more deflection would save **10.4 kg of fuel** per flight, and a spar material 50 % stiffer moves the optimum to b = 17.4 m, AR = 18.9 and 427 kg of fuel (−33 kg). That is the number a programme manager needs to weigh a carbon spar against an aluminium one.

### 2.7 Where the model stops being valid

- **Fully stressed, fatigue-limited caps with no stiffness design.** Real spars are often stiffness- or flutter-designed above AR ≈ 12; the deflection limit is a crude proxy for that, and the multiplier above is only as good as the 6 % number.
- **No fuel in the wing, no inertia relief.** Wing fuel and engine mass relieve the root moment; adding them (extension 2) lowers W_wing and weakens k.
- **Elliptic loading at every span.** High aspect ratios with no twist do not load elliptically, and e = 0.85 is not a constant.
- **The unconstrained limit is nonsense.** Remove the deflection limit and the optimiser goes to b = 19.6 m, AR = 22.4, a wing of 946 kg (37 % of W₀) and 423 kg of fuel: the model's fuel benefit of span never runs out because nothing in it penalises span faster than b³. A result that drives a sub-model to a regime its assumptions exclude is a result about the model.
- **Contraction.** Gauss–Seidel needs |k| < 1: k reaches 0.46 at b = 25 m and 1.02 at b = 35 m with S = 17.5 m², where the fixed point is no longer attracting and the physical statement is that the wing cannot carry its own weight growth.

## 3. Build the computational model

Open `starter/main.py`. `WingParams`, `MDAResult`, `design_constraints_x_only`, `mdf_optimize`, `idf_optimize` and `main` are given.

### Step 1: the disciplines and their partials

```python
def aero_mission(x, W_wing, params):
    b, S = float(x[0]), float(x[1])
    p = params
    L = p.W0 + W_wing
    D = p.q * (p.f_rest + p.c_dw * S) + L ** 2 / (p.q * np.pi * p.e * b ** 2)      # (1)
    W_fuel = L * (np.exp(p.kappa * D / L) - 1.0)                                  # (2)
    W_TO = L + W_fuel
    M_root = p.n_ult * W_TO * b / (3.0 * np.pi)                                   # (3)
    return {"L": L, "D": D, "W_fuel": W_fuel, "W_TO": W_TO, "M_root": M_root, "LD": L / D}
```

`aero_mission_partials` returns d(M_root, W_fuel, W_TO)/d(b, S, W_wing) from the hand derivation: with φ = D/L, ∂φ/∂b = −2L/(qπeb³), ∂φ/∂S = qc_dw/L, ∂φ/∂L = −q(f_rest + c_dw S)/L² + 1/(qπeb²), W_TO = L e^{κφ}, so ∂W_TO = L e^{κφ} κ ∂φ plus e^{κφ} on the W_wing component. `structures` implements (4) and the deflection ratio; `structures_partials` differentiates (4). Check every partial against central differences before writing anything else: a wrong partial makes Newton converge linearly and the adjoint wrong by a per cent, both of which look like success.

### Step 2: residual, Jacobians, two MDAs

```python
def coupled_jacobians(x, y, params):
    da = aero_mission_partials(x, y[1], params)["M_root"]      # [d/db, d/dS, d/dW_wing]
    ds = structures_partials(x, y[0], params)                   # [d/db, d/dS, d/dM_root]
    dGdy = np.array([[0.0, da[2]], [ds[2], 0.0]])
    dGdx = np.array([[da[0], da[1]], [ds[0], ds[1]]])
    return np.eye(2) - dGdy, -dGdx
```

`mda_gauss_seidel` alternates the two disciplines and records ‖R‖ after each sweep, stopping at `tol (1 + ‖y‖)`; `mda_newton` records ‖R‖ before each step and solves with the analytic ∂R/∂y. `coupling_strength` returns the product of the off-diagonal entries.

### Step 3: functions of interest and three derivative routes

`functions_of_interest` returns (W_fuel, g_stall) with g_stall = 1 − W_TO/(q_stall C_Lmax S); `functions_partials` returns ∂F/∂x and ∂F/∂y (the first column of ∂F/∂y is zero: neither function sees M_root directly). Then

```python
def total_derivatives_adjoint(x, params):
    y = mda_newton(x, params).y
    dRdy, dRdx = coupled_jacobians(x, y, params)
    dFdx, dFdy = functions_partials(x, y, params)
    psi = np.linalg.solve(dRdy.T, dFdy.T).T          # row i: adjoint of F_i
    return functions_of_interest(x, y, params), dFdx - psi @ dRdx, psi
```

`total_derivatives_direct` solves `dydx = solve(dRdy, -dRdx)` and returns `dFdx + dFdy @ dydx`; `total_derivatives_fd` differences the converged loop; `fd_convergence_study` tabulates the FD error against the adjoint for a list of steps.

Running `main()` prints

```
MDA: Gauss-Seidel 9 sweeps, Newton 2 iterations, |y_GS - y_N| = 1.52e-06
  GS residuals: 1.34e+04 7.65e+02 4.36e+01 2.48e+00 1.41e-01 8.04e-03 4.57e-04 2.60e-05 1.48e-06
  coupling strength k = 0.0569; observed GS ratio 0.0569
  converged: M_root = 156.12 kN m, W_wing = 2648.4 N (270.1 kg; caps 1640.4 N, cap area 77.4 cm^2)
  W_zf = 27648 N, L/D = 15.35, D = 1802 N, W_fuel = 5049 N (515 kg), W_TO = 32698 N, defl ratio 0.0408
total derivatives dW_fuel/d(b, S) [N/m, N/m^2] and dg_stall/d(b, S):
  direct : [[-290.65465, 65.972861], [-0.003353, 0.051096]]
  adjoint: [[-290.65465, 65.972861], [-0.003353, 0.051096]]
  FD     : [[-290.65465, 65.972861], [-0.003353, 0.051096]]
MDF (adjoint gradients): b = 14.869 m, S = 17.546 m^2, AR = 12.60, W_fuel = 4511.9 N (460.1 kg), 12 SLSQP iterations, 49 MDA solves
```

and the IDF line, the constraint values and the figure.

### Verify

- **Analytic cases:** each discipline against its formula; zero range gives zero fuel; doubling b removes ¾ of the induced drag; zero moment leaves w_skin S; the deflection ratio is unchanged when the moment doubles; the diagonal of ∂R/∂y is exactly 1.
- **Limiting and failure cases:** Newton converges from y = 0 and from a bad guess; the coupling strength at b = 18 m is more than 1.5 times that at 12 m.
- **Structure:** direct = adjoint to 10⁻¹⁰; (∂R/∂y)ᵀψ = (∂F/∂y)ᵀ; the frozen partial differs from the total by more than 5 %; signs (span saves fuel, area costs fuel and helps stall).
- **Convergence:** the Gauss–Seidel ratio equals k to 3 % over four consecutive sweeps; Newton's observed order exceeds 1.7; the FD error has slope 2.0 ± 0.1 and falls by 100 per decade of step to 5 %.
- **Agreement with the oracle:** `scipy.optimize.root` to 10⁻⁹ for the MDA; MDF with adjoint gradients, MDF with SLSQP's finite differences and IDF agree to 10⁻⁵ in x and 10⁻⁷ in fuel, and the two MDF starts agree.

### Compare

`scipy.optimize.root` (hybrid Powell) solves the MDA to the same y as our Newton without being told any derivative, at the cost of extra residual evaluations; it is the right tool for a loop whose Jacobian you cannot write, and ours is the right one when you can. SLSQP with its own finite-difference gradients reaches the same optimum here, with 85 MDA solves instead of 49; with twenty variables it would be 10 times worse and with a loosely converged MDA it would stall. We call the library for the optimisation and the oracle root solve, and write the MDA and the derivatives ourselves, because those are where the physics and the cost live. OpenMDAO (part 6) packages exactly this division of labour: you write the disciplines and their partials, it assembles R, runs Newton or Gauss–Seidel, and forms the direct or adjoint totals for any number of disciplines, and it is what you would reach for at three disciplines or more.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Predict before running.

1. **The span trade.** At S = 17.55 m², converge the MDA and compute the fuel for b = 10, 12, 14, 14.87, 16, 18, 20, 22 m. Predict a curve falling steeply at first (induced drag ∝ 1/b²) and flattening as the b³ cap weight bites, with its unconstrained minimum beyond 18 m, and mark where the AR = 12.6 limit cuts it.
2. **Coupling strength.** Record k at each span of experiment 1 and at b = 25, 30, 35 m, and whether `mda_gauss_seidel` converges in 50 sweeps. Predict k ∝ b³/S roughly (from (3), (4) and the b-dependence of W_TO) and failure near k = 1.
3. **Finite-difference step.** Run `fd_convergence_study` with h = 10⁻¹ … 10⁻⁷ with the Newton MDA at 10⁻¹³, then with `mda_solver="gs"` at tolerances 10⁻¹², 10⁻⁶ and 10⁻³. Predict a slope of 2 until a floor; predict that the floor of the loose runs scales with the tolerance (the error in y is the same on both sides of the difference, so it does not grow as 1/h but does not vanish either).
4. **Cost of gradients.** Count MDA solves for `mdf_optimize` with `use_adjoint=True` and `False`. Predict the ratio from (2 n_x + 1)/1 evaluations per iteration against 1.
5. **Parameters.** Re-optimise with R = 1000, 1500, 2500 km and with E × 1.5. Predict that range moves the fuel almost proportionally and the design hardly at all (the active constraint is geometric), while stiffness moves the design.

Raw data format: CSV tables (`b_m, fuel_kg, W_wing_N, k, gs_converged`; `h_rel, fd_error_tight, fd_error_loose`; `gradients, iterations, mda_solves`; `case, b, S, fuel_kg`) with a header naming the fixed parameters; `wing_mdo.png`.

What you should find: fuel 601, 520, 474, 460, 446, 430, 423, 424 kg for the eight spans, with the unconstrained minimum near 20–21 m where the wing weighs over 900 kg; k = 0.033, 0.054, 0.084, 0.100, 0.123, 0.173, 0.235, 0.312 and then 0.46, 0.83, 1.02 at 25, 30, 35 m, Gauss–Seidel failing at 30 m within 50 sweeps and diverging at 35 m while Newton still converges; FD errors 6.78, 6.66 × 10⁻², 6.66 × 10⁻⁴, 6.66 × 10⁻⁶ (slope 2.00) for the tight MDA with floors of 3 × 10⁻⁸ (Newton or Gauss–Seidel at 10⁻¹²), 2.4 × 10⁻⁴ (Gauss–Seidel at 10⁻⁶) and 0.04 N/m (10⁻³); 49 against 85 MDA solves; fuel 298, 460, 811 kg at the three ranges with b moving only from 14.82 to 14.98 m, and b = 17.44 m, 427 kg at 1.5 E.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| first Gauss–Seidel sweep, M_root / W_wing | 141.86 kN m / 2498.6 N | same | 0 | one evaluation each |
| coupling strength k | 5.418 × 0.01051 = 0.0569 | 0.0569 (ratio of residuals) | < 3 % | asymptotic regime |
| Newton order | 2 | 2.0 | — | exact Jacobian |
| dW_fuel/db: direct / adjoint / FD | — | −290.65465 / −290.65465 / −290.65465 N/m | 10⁻¹⁰, 10⁻⁵ | transposed algebra; O(h²) |
| frozen partial ∂W_fuel/∂b | — | −346.77 N/m | 19 % from the total | weight feedback |
| FD order | 2 | 2.00 | — | tight MDA |
| MDF vs IDF optimum | same | b 14.869 / 14.869 m, fuel 4511.9 / 4511.9 N | < 10⁻⁷ | both satisfy the same KKT system |
| active constraint | AR = 12.6 | AR = 12.60, g_defl = 0 to 10⁻¹² | 0 | — |
| MDA solves, adjoint vs FD gradients | fewer | 49 vs 85 | — | 2n_x + 1 evaluations per FD gradient |

The numerics validate completely and the physics is self-consistent; the thing to doubt is the sub-models. The model says the span should be as large as stiffness allows, which is true of real aircraft, and it says the fuel benefit of span never saturates on its own, which is not, because real wing weight grows faster than b³ once the spar becomes stiffness-designed and because flutter, gust loads, taxiing clearance and gate boxes intervene. A validated statement from this lesson is therefore: "the optimum is on the stiffness limit, each per cent of deflection allowance is worth about 1 kg of fuel per flight, and the model is not to be used beyond AR ≈ 15 without a stiffness and flutter model". The variance between architectures and derivative routes is zero to the tolerances claimed; the variance between modelling choices (fuel in the wing, elliptic loading, fatigue allowable) is tens of kilograms, and that is where the next hour of work belongs.

Improvements worth making: wing fuel with inertia relief (extension 2), which also makes k smaller and the design heavier-wing-friendly; a flutter or torsional-stiffness constraint (extension 3); the OpenMDAO rebuild (extension 1) if you will go on to more disciplines. Improvement worth *not* making: tightening the SLSQP tolerance below 10⁻¹² or the MDA below 10⁻¹³. The optimum is known to six figures and the inputs to two.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of the partials of (1)–(4), of k, of (7) and (8), and the first Gauss–Seidel sweep; `wing_mdo.png`; the five experiment tables; and a one-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the optimum, the active constraint and its price, the span beyond which you would not trust the model, and what one additional sub-model would change the answer most.

**On OpenMDAO (optional).** The framework from NASA Glenn does for N disciplines what this lesson did for two: you write each discipline as a component with `compute` and `compute_partials`, connect them in a group, choose a nonlinear solver (NonlinearBlockGS is our Gauss–Seidel, Newton is our Newton), and ask for total derivatives in forward or reverse (adjoint) mode; its N2 diagram is the picture of ∂R/∂y. Its "Multidisciplinary Optimization" guide solves the Sellar problem, a two-discipline algebraic model with the same structure as this wing. Install it separately if you want to rebuild this lesson in it; nothing in the tests depends on it, and the lesson's numbers are the ones to reproduce.

### Adversarial self-checks

- **What if the empty weight doubles?** W₀ = 50 kN: W_TO, M_root and the cap weight roughly double (W_wing 4260 N at the baseline geometry), but k rises only from 0.057 to 0.064, because both coupling derivatives are set by the geometry b²/S and not by the load. The optimum runs to the 20 m gate limit with S = 31.7 m² (AR still 12.6: two active constraints) and the stall margin shrinks to 5 %. The aspect-ratio limit does not move because it is geometric; the gate limit is now the price-setting constraint.
- **What if the range goes to zero?** No fuel, W_TO = W_zf, and the objective is identically zero: the optimiser has nothing to minimise and any feasible design is optimal. The lesson's objective is only meaningful because a mission is specified.
- **Which assumption fails first in the real world?** The fully stressed, fatigue-allowable spar with no stiffness design at AR above about 12; then the elliptic loading and constant e; then the frozen L/D at the zero-fuel weight (a 15 % error in the fuel fraction's exponent, which moves the fuel by a few per cent but hardly moves the design).
- **What measurement could falsify the model?** A static test of the spar to limit load: a tip deflection far from σ_a b²/(3Eτ S) × (b/2) falsifies the fully-stressed assumption; a flight-test L/D at the design C_L far from 17.9 falsifies the polar. Both are routine; neither is a simulation.
- **Does the result survive a change of solver, architecture or derivative route?** Gauss–Seidel and Newton agree to 10⁻⁹; direct, adjoint and FD gradients agree to 10⁻⁵; MDF and IDF agree to 10⁻⁵; SLSQP with FD gradients agrees at higher cost. A disagreement among these would be a bug, never physics.
- **Which decision changes if the uncertainty doubles?** If the allowable stress is uncertain by ±20 % instead of ±10 %, the aspect-ratio limit 3Eτδ_max/σ_a is uncertain by the same fraction, the optimum span by ±10 %, and the fuel by about ±10 kg: still a 10 % better design than the baseline either way, so the decision "go to the stiffness limit" survives while the decision "carbon or aluminium spar" (worth 33 kg per flight against its cost) does not without a better stiffness model.

## Check yourself

1. Write the coupled residual R(x, y) = y − G(x, y) for the wing, naming x, y and the two maps G₁ and G₂ with their inputs.
2. What is the coupling strength k, how do you compute it from the discipline partials, and what does k > 1 mean for Gauss–Seidel and for the physics?
3. State the implicit-function formula for dy/dx and explain why dW_fuel/db computed with W_wing frozen is wrong, and in which direction.
4. Write the adjoint equation and the total-derivative formula (8). For 200 design variables and 3 functions, how many linear solves does each method need?
5. Why does the finite-difference derivative of a converged loop have an error floor, and what sets it?
6. What is the difference between the MDF and IDF architectures, and why must they agree at the optimum?
7. The tip deflection ratio of a fully stressed beam is σ_a b²/(3Eτ S). Why does the load not appear, and what kind of constraint does a deflection limit therefore become?
8. The optimiser increased the wing area even though the stall constraint had 14 % margin. Explain why in one sentence.
