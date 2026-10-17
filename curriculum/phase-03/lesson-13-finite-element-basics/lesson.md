# P03.13  Finite Element Basics

> A finite element method does not approximate the differential equation; it approximates the solution, by the best piecewise-polynomial the mesh can hold, in the energy the problem itself defines. Everything else, stiffness matrices, assembly, convergence rates, follows from that one sentence.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the weak form of the axially loaded bar from the strong form by integration by parts, show that the Galerkin method with hat functions produces the element stiffness matrix (EA/h)[[1, −1], [−1, 1]] and the consistent load vector, and prove that the global matrix is symmetric, positive semi-definite with one rigid-body mode, and tridiagonal. 2. Explain why the finite-element solution minimises the potential energy over the trial space and why, for a uniform bar, it is exact at the nodes while the stress converges one order more slowly. 3. Implement and test, from scratch, element matrices, Gauss-quadrature loads, assembly, constraints by partitioning, stresses, reactions, energies and the axial eigenproblem, against two closed-form solutions and `scipy.integrate.quad`, `scipy.integrate.solve_bvp`, `scipy.linalg.solve_banded` and `scipy.linalg.eigh`, and measure the observed orders (2 in L2, 1 in energy, 2 for the first frequency). 4. Decide how many elements a rotor-blade spar and a tapered strut need for a stated stress accuracy, and say where the one-dimensional bar model of a rotor blade stops being physics. |
| **Prerequisites** | `P03.08` (symmetric positive-definite systems, Cholesky), `P03.04` (Gauss–Legendre quadrature), `P03.09` (banded and sparse matrices; why K is tridiagonal), `P02.13` (calculus of variations: a minimum of a functional satisfies a differential equation), `P02.11` (the fixed–free bar's axial modes by separation of variables), `P01.16` (the generalised eigenproblem), `P01.13` (partitioned linear systems). Skill check: you can integrate by parts and you know that a symmetric matrix with non-negative eigenvalues is positive semi-definite. |
| **Engineering question** | The aluminium spar of a light helicopter's rotor blade is 5 m long with a 600 mm² section and carries a 2 kg tip mass at 400 rpm. What are the root stress, the tip extension and the hub reaction under centrifugal load, how many elements does the finite-element model need to give the root stress to 1 %, and what are the spar's first three axial natural frequencies? Alongside: how does the same code handle a tapered strut, where the exact stress varies as 1/A(x)? |
| **Assumptions and limits** | Linear elasticity (stresses far below yield, 88 MPa against about 270 MPa for 2024 alloy); a one-dimensional bar (uniaxial stress, no bending, no centrifugal stiffening of bending modes, no Coriolis coupling); the spar modelled from the hub centre with a fixed root (a real blade attaches through a grip and pins at a cut-out radius, where the stress concentration lives); the tip mass as a point load; the section constant along the span; the natural frequencies computed for the non-spinning bar. Linear elements carry a constant stress per element, so a varying stress is always under- or over-estimated inside an element; the consistent-mass frequencies are upper bounds that converge from above. |
| **Analytical work** | Strong form, weak form and the principle of virtual work; Galerkin with hat functions; the element stiffness and consistent mass matrices and the consistent load vector; assembly; the properties of K; the minimum-energy characterisation and Clapeyron's theorem; nodal exactness for a uniform bar; the error orders in L2 and energy; the closed-form solutions of the rotating bar and the tapered bar; the exact axial frequencies of a fixed–free bar. |
| **Computational deliverable** | `starter/main.py`: `exact_rotating_bar`, `exact_tapered_bar`, `shape_functions`, `element_stiffness`, `element_mass`, `element_load`, `assemble`, `apply_point_load`, `solve_bar`, `element_stress`, `interpolate`, `strain_energy`, `potential_energy`, `l2_error`, `energy_error`, `observed_order`, `solve_rotor_blade`, `convergence_study`, `axial_frequencies`, `exact_axial_frequencies` (`rotor_blade`, `tapered_bar_bvp_oracle` and `main` are given). Tests prove: both closed forms satisfy their equations and boundary conditions; the shape functions sum to one; element matrices have the textbook form and the right rigid-body and total-mass properties; the element load equals `quad` to 10⁻¹² and preserves the total load; K is symmetric, tridiagonal, rank-deficient by one before constraint and SPD after; the solution equals `solve_banded` to 10⁻¹⁰ and the reaction balances the load to 10⁻¹⁰; a linear field passes the patch test on a random mesh to 10⁻¹²; the uniform bar is exact at the nodes and its element stresses are the exact element averages; the solution minimises the potential energy and satisfies Clapeyron's theorem; the tapered bar matches `solve_bvp` to 10⁻⁷ and converges at order 2.0 ± 0.1 at the nodes; observed orders are 2.00 ± 0.05 in L2 and 1.00 ± 0.05 in energy; the frequencies equal `eigh`'s to 10⁻⁹, bound the exact ones from above, fall monotonically with refinement and converge at order 2. |
| **Practical activity** | **computational experiment**: a mesh-refinement study for displacement, stress and frequency; a comparison of consistent and lumped mass; a stress-recovery study (element centre versus nodal extrapolation) on the tapered strut; a non-uniform-mesh study grading elements toward the root. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `bar_fem.png`, the hand derivation of the weak form, the element matrices and the nodal-exactness argument, the convergence and frequency tables, and a half-page note stating the spar's root stress, tip extension and first frequency with the mesh that justifies them and the model assumption you trust least. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the weak form and element stiffness are derived; all public tests pass; the observed orders in L2, energy and frequency are explained; the note gives each number with its mesh and an error estimate. |
| **References** | `@ocw-16-90` Unit 2 "2.8 Method of Weighted Residuals", "2.9 Introduction to Finite Elements", "2.10 More on Finite Element Methods" (required-free: weighted residuals, the weak form, linear elements and assembly); "2.11 The Finite Element Method for Two-Dimensional Diffusion" (supplementary-free: the same ideas on triangles). `@ocw-18-085` "Lecture 17: Finite Elements in 1D (part 1)", "Lecture 18: Finite Elements in 1D (part 2)" (required-free: hat functions, the stiffness matrix, nodal exactness); "Lecture 19: Quadratic/Cubic Elements", "Lecture 20: Element Matrices; 4th Order Bending Equations" (supplementary-free: higher elements and the beam, for extension 1). `@ocw-2-092` SES 4 "The principle of virtual work", SES 5 "The finite element formulation", SES 6–7 "Finite element solution process" (supplementary-free: the structural-mechanics view of the same formulation). |
| **Extensions** | ★ Add three-node quadratic elements and show the uniform spar's stress becomes exact element by element and the tapered strut's error falls at order 3 in L2. ★ Grade the mesh (elements shrinking toward the tip of the tapered strut, where the stress gradient is largest) and show the same accuracy with half the elements. ★ Add the Euler–Bernoulli beam element (`@ocw-18-085` Lecture 20) and compute the blade's first flapping frequency, then include centrifugal stiffening through a geometric stiffness matrix. |
| **Time** | unmeasured |

## 1. The engineering problem

The spar of a rotor blade is the member that keeps the blade attached to the aircraft. At 400 rpm every kilogram of the 5 m blade pulls outward with ρAΩ²x per metre, and the 2 kg tip cap and balance weight add a point load of m_tipΩ²L = 17.55 kN at the end. The questions are the ordinary ones of stress analysis (the root stress, the tip extension, the hub reaction) and the extraordinary one that governs rotorcraft design: the spar's natural frequencies, which must stay clear of the rotor's harmonics. For a uniform bar every one of these has a closed form, which is exactly why it is the right first finite-element problem: the method can be checked at every step.

The **finite element method** (`@ocw-16-90` 2.9, `@ocw-18-085` Lecture 17) is the general way of turning a boundary-value problem into a linear system that honours the physics (equilibrium of each piece, compatibility at the joints) on any geometry. The bar is its simplest instance, and the one-dimensional solver built here has every part of a structural code: element matrices, quadrature loads, assembly, constraints, stress recovery, an eigenproblem. The second case, a tapered strut with a tip load, has a stress that varies as 1/A(x) and a logarithmic displacement, so that the solver's approximation of a varying section can be tested against a closed form and against `scipy.integrate.solve_bvp`.

Success criteria, written down first: the uniform spar's nodal displacements agree with the closed form to rounding (a theorem for this problem, so anything less is a bug); the element stresses equal the exact element-average stress; the L2 displacement error falls at order 2 and the energy-norm stress error at order 1 on refinement; the tapered strut's tip displacement is within 0.1 % of the closed form with 16 elements; the first axial frequency is within 0.1 % with 16 elements, and the finite-element frequencies are upper bounds that fall monotonically as the mesh is refined.

## 2. Understand and derive

### 2.1 Strong form and weak form

A bar of section A(x) and modulus E, loaded by a body force b(x) per unit length and a tip force P, obeys equilibrium of an element dx: d(Aσ)/dx + b = 0 with σ = Eu'. The **strong form** is

−(EAu')' = b on 0 < x < L,  u(0) = 0,  EAu'(L) = P.   (1)

Multiply (1) by any **test function** v(x) with v(0) = 0 (a virtual displacement that respects the support), integrate over the bar, and integrate the left side by parts:

∫₀ᴸ EAu'v' dx − [EAu'v]₀ᴸ = ∫₀ᴸ bv dx  ⇒  **∫₀ᴸ EAu'v' dx = ∫₀ᴸ bv dx + Pv(L)**  for all admissible v.   (2)

This is the **weak form**: the internal virtual work equals the external virtual work for every virtual displacement (the principle of virtual work of `@ocw-2-092` SES 4). Two things happened. The second derivative became a first derivative on each side, so u need only be once differentiable, piecewise: a kinked function is admissible. And the tip force P entered through the boundary term, by itself: a force boundary condition is **natural**, it appears in the equation; the displacement condition u(0) = 0 is **essential**, it must be built into the space of u and v. If (2) holds for all v and u is smooth, integrating by parts backwards recovers (1) and EAu'(L) = P, so the two forms are equivalent; the weak form merely asks less of u.

### 2.2 Galerkin with hat functions

Divide the bar into n elements by nodes x₀ = 0 < x₁ < … < x_n = L and let N_j(x) be the **hat function** that is 1 at node j, 0 at every other node, and linear in between. Approximate u by u_h = Σ_j u_j N_j, so that the coefficients u_j are the nodal displacements (N_j(x_i) = δ_ij), and take the test functions from the same space, v = N_i (**Galerkin**). The weak form becomes n + 1 equations:

Σ_j [∫ EAN_i'N_j' dx] u_j = ∫ bN_i dx + Pδ_{in}, i.e. **Ku = f**,   (3)

with **K_ij = ∫EAN_i'N_j' dx** (the stiffness matrix) and **f_i = ∫bN_i dx + Pδ_{in}** (the consistent load vector). Because N_i and N_j overlap only when the nodes are the same or adjacent, K is tridiagonal, and because each integrand is a product of two piecewise-constant derivatives it can be evaluated element by element.

On element e, spanning [x_e, x_{e+1}] of length h, the two hat functions restricted to it are the **shape functions** N₁ = (x_{e+1} − x)/h, N₂ = (x − x_e)/h, or in the parent coordinate ξ ∈ [−1, 1], N₁ = (1 − ξ)/2, N₂ = (1 + ξ)/2, with derivatives ∓1/h. The element's contribution to K is

**k_e = (EA/h) [[1, −1], [−1, 1]],**   (4)

exact for constant EA, and exact with A replaced by its element mean for linearly varying A (the integrand is then linear and the mean integrates it exactly). Its contribution to f is **f_e = ∫N^T b dx**, computed by Gauss–Legendre quadrature (`P03.04`) in ξ: for uniform b the load splits half and half, for b ∝ x it splits 1/3 and 2/3 toward the loaded end (the tests check both). **Assembly** adds k_e into rows and columns e, e + 1 of K and f_e into entries e, e + 1 of f: the scatter-add that every finite-element code performs.

For vibration, d'Alembert's inertia force −ρAü replaces b, giving Mü + Ku = 0 with the **consistent mass matrix** m_e = ∫ρAN^TN dx = (ρAh/6)[[2, 1], [1, 2]], and the harmonic ansatz u = φ cos ωt gives the generalised eigenproblem **Kφ = ω²Mφ** (`P01.16`).

### 2.3 What K and the solution satisfy

- **Symmetry**: K_ij = K_ji from the definition. **Positive semi-definiteness**: vᵀKv = ∫EA(v_h')² dx ≥ 0, with equality only for a constant v_h: the vector of ones is the single null vector (rigid translation). Fixing one node removes it and K becomes positive definite (`P03.08`), so the partitioned solve K_ff u_f = f_f − K_fc u_c has a unique answer and the Cholesky factor exists. The reactions are r = Ku − f, non-zero only at constrained nodes, and Σr = −Σf: equilibrium of the whole bar, which the test checks to 10⁻¹⁰.
- **Minimum potential energy**: (3) is the stationarity condition of Π(u) = ½uᵀKu − fᵀu over the trial space, and since K is positive definite on it the stationary point is the minimum (`P02.13` in finite dimensions). Any other admissible nodal vector has a larger Π; the test perturbs the solution randomly and checks. At the minimum, **Clapeyron's theorem** ½uᵀKu = ½fᵀu holds (strain energy equals half the work of the loads applied gradually), and Π = −U. For the spar: U = 107.88 J on 64 elements.
- **Galerkin orthogonality and the best approximation**: subtracting (3) from (2) with v in the trial space gives ∫EA(u − u_h)'v_h' dx = 0: the error is orthogonal, in the energy inner product, to every trial function. Hence u_h is the **best approximation** of u in the energy norm ‖w‖_E = (∫EAw'² dx)^{1/2} among all piecewise-linear functions on the mesh. The best piecewise-linear approximation of a smooth function has slope error O(h) and value error O(h²), which gives the orders **1 in energy (stress)** and **2 in L2 (displacement)** that the convergence study measures (0.99–1.00 and 1.99–2.00).
- **Nodal exactness for a uniform bar** (`@ocw-18-085` Lecture 17): the Green's function of −EAu'' with these boundary conditions is piecewise linear with its kink at the source point, so for a source at a node it *is* a hat function and lies in the trial space; Galerkin orthogonality then makes u_h(x_j) = u(x_j) exactly, for any load b. The spar's nodal displacements agree with the closed form to 10⁻¹⁸ m on 4 elements. The element stress E(u_{e+1} − u_e)/h is then the exact stress **averaged over the element**, which differs from the exact stress at the element centre by h²σ''/24 = −ρΩ²h²/24 (0.31 MPa at n = 4, 0.02 MPa at n = 16; tested) and from the root value by ρΩ²h²/6 (1.23 MPa at n = 4), because the exact stress has zero slope at the root. For a varying section the nodes are no longer exact (the Green's function is no longer piecewise linear), and the tapered strut converges at order 2 at the nodes instead.

### 2.4 The two closed forms

**Rotating uniform bar.** With b = ρAΩ²x and tip load P, (1) integrates to EAu' = ρAΩ²(L² − x²)/2 + P, so

u = (ρΩ²/2E)(L²x − x³/3) + Px/(EA),  σ = ρΩ²(L² − x²)/2 + P/A.   (5)

The root stress is the whole centrifugal load of the spar (ρAΩ²L²/2 = 35.53 kN) plus P, over A: 88.46 MPa, a third of the alloy's yield; the tip extension is 4.909 mm; the hub reaction is 53.08 kN, the weight of 5.4 tonnes hanging on a 600 mm² rod. These are what the spar's root fitting must carry, forever, 400 times a minute.

**Tapered bar.** With A = A₀(1 − cx/L), b = 0 and tip load P, EAu' = P is constant, so σ = P/A(x) and u = (PL/(EA₀c)) ln(1/(1 − cx/L)): 2.0453 mm at the tip for A₀ = 800 mm², c = 0.6, L = 1.5 m, P = 50 kN, with the tip stress 156.2 MPa against 62.5 MPa at the root. The finite-element solution with the mid-element area gives 2.0346, 2.0446 and 2.0452 mm on 4, 16 and 64 elements (errors 10.7, 0.69 and 0.043 µm: order 2), and `solve_bvp` gives 2.0453 mm. The last element's stress, P/A(mid-element), is 131.6, 149.3 and 154.4 MPa: always below the true tip value, because the element reports the stress half an element inboard, and the difference is 1/A' × h/2 relative: the first-order stress error of §2.3 made visible.

**Axial frequencies.** Separation of variables for ρü = Eu'' with u(0) = 0, u'(L) = 0 (`P02.11`) gives ω_k = (2k − 1)(π/2L)√(E/ρ): 254.59, 763.76, 1272.94 Hz for the spar (√(E/ρ) = 5092 m/s is the bar wave speed). The consistent-mass model gives 256.23, 808.30, 1468.31 Hz on 4 elements and 254.69, 766.53, 1285.75 Hz on 16: **upper bounds** (the trial space constrains the structure, which stiffens it, a Rayleigh-quotient argument), converging from above at order 2 for the first mode, and the k-th mode needs about 10k elements for 1 %.

### 2.5 Where the model stops being valid

- **A bar is not a blade.** The spar bends (flapping, lead-lag) and twists; its bending frequencies, not the axial ones computed here, are the ones that must avoid the rotor harmonics, and they rise with rotor speed through centrifugal stiffening (extension 3). The axial model is right for the root tension and the tip extension and for nothing dynamic.
- **The root is not a clamp.** The blade attaches through a grip at a cut-out radius of perhaps 0.3 m with pin holes; the stress concentration at the holes, not the 88 MPa nominal, governs fatigue life. Nothing in a bar model sees it.
- **The section is not constant.** Real spars taper and change shape; the tapered-strut machinery handles A(x), and the error analysis says the stress near a rapid change needs small elements there (extension 2).
- **Linear elasticity.** Fine at 88 MPa; the bar formulation itself would survive plasticity with a tangent modulus and a Newton loop (`P03.10`), the closed forms would not.

## 3. Build the computational model

Open `starter/main.py`. `rotor_blade`, `tapered_bar_bvp_oracle` and `main` are given.

### Step 1: element level

```python
def shape_functions(xi):
    return np.column_stack([0.5 * (1 - xi), 0.5 * (1 + xi)])         # N1, N2 on [-1, 1]

def element_stiffness(E, A, h):
    return E * A / h * np.array([[1.0, -1.0], [-1.0, 1.0]])          # (4)

def element_load(b, x1, x2, n_gauss=3):
    xi, w = np.polynomial.legendre.leggauss(n_gauss)
    x = 0.5 * (x1 + x2) + 0.5 * (x2 - x1) * xi                        # map parent to physical
    return 0.5 * (x2 - x1) * (shape_functions(xi).T @ (w * b(x)))     # int N^T b dx
```

`element_mass` is (ρAh/6)[[2, 1], [1, 2]]. Raise `ValueError` for h ≤ 0.

### Step 2: assembly and solution

`assemble(nodes, E, A, b, rho)` loops over elements, evaluates A at the element mid-point if it is callable, and scatter-adds `k_e`, `f_e` and (if ρ is given) `m_e` into dense global arrays with `K[e:e+2, e:e+2] += k_e`. Dense is fine here (a few hundred nodes); `P03.09` showed how the same tridiagonal matrix is stored sparsely. `apply_point_load` adds P to one entry of f. `solve_bar(K, f, fixed)` partitions the dofs into constrained and free, solves K_ff u_f = f_f − K_fc u_c with `np.linalg.solve`, and returns u and the reactions r = Ku − f. `element_stress` is E Δu/Δx per element; `interpolate` is `np.interp`; `strain_energy` and `potential_energy` are the quadratic forms of §2.3.

### Step 3: errors, convergence and the eigenproblem

`l2_error` integrates (u_h − u)² element by element with Gauss points (u_h from the shape functions); `energy_error` integrates (σ_h − σ)²/E. `solve_rotor_blade(case, n)` builds the spar on n elements with b = ρAΩ²x and the tip load. `convergence_study` returns the two error sequences. `axial_frequencies(K, M, fixed)` removes the fixed dofs, reduces Kφ = λMφ to a standard symmetric problem with the Cholesky factor of M (L⁻¹KL⁻ᵀ, `P03.08`) and returns √λ ascending; `exact_axial_frequencies` is the formula of §2.4.

Running `main()` prints

```
rotor-blade spar: L = 5.0 m, A = 600 mm^2, omega = 41.89 rad/s (400 rpm), tip load P = 17.55 kN, total centrifugal load of the spar 35.53 kN
  n =   4: tip u = 4.9087 mm (exact 4.9087), max |nodal error| = 2.2e-19 m, root stress = 87.23 MPa (exact at root 88.46, at element centre 87.54), reaction = -53.077 kN, strain energy = 107.36 J
  n =  16: tip u = 4.9087 mm (exact 4.9087), max |nodal error| = 8.7e-19 m, root stress = 88.38 MPa (exact at root 88.46, at element centre 88.40), reaction = -53.077 kN, strain energy = 107.85 J
  n =  64: tip u = 4.9087 mm (exact 4.9087), max |nodal error| = 3.1e-17 m, root stress = 88.46 MPa (exact at root 88.46, at element centre 88.46), reaction = -53.077 kN, strain energy = 107.88 J

convergence (n = [4, 8, 16, 32, 64]):
  L2 displacement error  [6.18e-05 1.55e-05 3.89e-06 9.74e-07 2.43e-07] orders [1.99 2.   2.   2.  ]
  energy-norm stress err [4.14e+01 2.08e+01 1.04e+01 5.21e+00 2.61e+00] orders [0.99 1.   1.   1.  ]
```

followed by the tapered strut, the frequency table and the figure.

### Verify

- **Analytic cases:** (5) satisfies σ = Eu' and (Aσ)' = −b to 10⁻⁴ and its boundary conditions exactly; the tapered closed form has constant axial force to 10⁻¹²; the uniform bar is exact at the nodes to 10⁻¹⁸ m for n = 3, 10, 25 and its element stresses are the exact element averages.
- **Limiting cases:** a rigid translation produces no element force; a uniform load splits evenly and a linear load 1 : 2; a prescribed end displacement with no load gives uniform strain and the reaction EAu/L; the patch test on a random 12-node mesh reproduces a linear field to 10⁻¹².
- **Structure:** K symmetric, tridiagonal, rank n, SPD after one constraint; M symmetric positive definite with total mass ρAL; reactions only at constrained dofs and summing to minus the load.
- **Energy:** the solution minimises Π among five random perturbations; Clapeyron's theorem and Π = −U hold to 10⁻¹².
- **Convergence (observed order):** 2.00 ± 0.05 in L2, 1.00 ± 0.05 in energy, 2.0 ± 0.1 at the nodes for the tapered bar, 2.0 ± 0.1 for the first frequency.
- **Agreement with the oracles:** `quad` for the element load to 10⁻¹²; `solve_banded` for the solution to 10⁻¹⁰; `solve_bvp` for the tapered bar to 10⁻⁷; `eigh` for the frequencies to 10⁻⁹.

### Compare

SciPy has no finite-element package; what it offers is the numerical linear algebra and quadrature a finite-element code is built from, and the lesson uses those as oracles: `solve_banded` (the tridiagonal solve, as in `P03.11`), `eigh` (the generalised symmetric eigenproblem, which does internally the Cholesky reduction `axial_frequencies` does by hand), `quad` (the element integrals) and `solve_bvp` (a collocation solver for the same boundary-value problem, which is the *other* way to solve (1) numerically and knows nothing of elements). In practice a structural problem goes to a finite-element code (there are open-source ones, and the bar element here is the first entry in every one of their element libraries) and the engineer's job is the one this lesson trains: choosing the element, the mesh and the constraints, and reading the stress with its known first-order error. Call `eigh` for the eigenproblem; write the assembly yourself at least once, because the assembly is where the physics enters and where the mistakes (a missing boundary term, a lumped load that should be consistent, a node fixed that should be free) live.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before the run.

1. **Refinement.** For the spar on n = 4, 8, 16, 32, 64, record the root element stress, its error against the exact root stress, the L2 and energy errors, and the first three frequencies. Predict: the root element stress error falls as 1/n² (it is the difference between the root value and the element average of a stress with zero slope at the root, ρΩ²h²/6 from (5)), the L2 error as 1/n², the energy error as 1/n, and the k-th frequency error as (k/n)².
2. **Consistent against lumped mass.** Replace M by the diagonal matrix of its row sums (each node gets half of each adjacent element's mass) and repeat the frequency study. Predict that the lumped frequencies are *below* the exact ones and the consistent ones above, so that the two bracket the answer.
3. **Stress recovery on the tapered strut.** On n = 4 … 128, record the last element's stress against the exact tip stress and against the exact stress at the last element's centre. Predict: the second agrees to rounding (the element stress is P/A at the element mid-point), the first is low by about (c/(1 − c))·(h/2L) relative. Then extrapolate the last two element stresses linearly to the tip and record the improvement.
4. **A graded mesh.** For the tapered strut, put the nodes at x_j = L(1 − (1 − j/n)^p) for p = 1, 1.5, 2 (elements shrinking toward the tip) and record the tip-stress error at n = 16. Predict which p wins and why, from where σ' is largest.

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `bar_fem.png`.

What you should find: root element stresses 87.23, 88.15, 88.38, 88.44, 88.46 MPa for n = 4 … 64 against 88.46 (error quartering each time, 1/n²); L2 orders 1.99, 2.00, 2.00, 2.00 and energy orders 0.99, 1.00, 1.00, 1.00; first-mode frequency ratios 1.0064, 1.0016, 1.0004, 1.0001 for n = 4, 8, 16, 32 (consistent, from above) and 0.9936, 0.9984, 0.9996, 0.9999 (lumped, from below); the tapered strut's last-element stress equal to the exact centre value to 10⁻³ % on every mesh while the tip error falls from −9.3 % (n = 4, c = 0.45) to −0.3 % (n = 128).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| tip load P = m_tipΩ²L | 2 × 41.89² × 5 = 17.55 kN | 17.55 kN | 0 | — |
| spar centrifugal load ρAΩ²L²/2 | 35.53 kN | 35.53 kN | 0 | — |
| root stress (35.53 + 17.55)/600 mm² | 88.46 MPa | 88.46 (n = 64 element), 87.23 (n = 4) | 0 / 1.4 % | element average: ρΩ²h²/6 below the root value |
| tip extension (5) | 4.9087 mm | 4.9087 mm on 4 elements | 2 × 10⁻¹⁹ m | nodal exactness |
| hub reaction | −53.08 kN | −53.077 kN | 10⁻¹⁰ rel. | equilibrium |
| strain energy | ½fᵀu | 107.36 (n = 4) → 107.88 J (n = 64) | 0.5 % | U_h ≤ U: the trial space is stiffer |
| L2 / energy order | 2 / 1 | 2.00 / 1.00 | — | best approximation in energy |
| tapered tip u (16 elements) | 2.0453 mm | 2.0446 mm | 0.03 % | mid-element area, O(h²) |
| tapered tip stress (64 elements) | 156.2 MPa | 154.4 MPa | 1.2 % | element stress is P/A(mid), half an element inboard |
| first frequency (16 elements) | 254.59 Hz | 254.69 Hz | 0.04 %, above | consistent mass is an upper bound |
| third frequency (16 elements) | 1272.94 Hz | 1285.75 Hz | 1.0 %, above | 5 elements per half-wave |

Every discrepancy is a derived property of linear elements, not noise, and two deserve a sentence in any report. The strain energy of the finite-element solution is *below* the true one (107.36 against 107.88 J): the trial space cannot bend as the true solution does, so it stores less energy under the same loads, which is the same fact as "frequencies are upper bounds" seen from the static side. And the element stress is a first-order quantity that always sits half an element away from where you want it: for the strut's tip stress that is a 1.2 % under-read even on 64 elements, and the fix is not more elements but a better recovery (extrapolate from the element centres, or use the natural boundary condition σ(L) = P/A(L) directly), which is what structural codes do under the name of stress recovery.

What agreement with the closed forms does *not* validate: the blade. The root stress of 88 MPa is the nominal tension at the hub centre of an idealised bar; the real spar's life is set by the stress at the pin holes of its grip, by bending under lift, and by the fatigue cycles of each revolution. The finite-element machinery scales to all of that; the one-dimensional bar does not, and the note in §6 must say which of its numbers survive (the tip extension and the hub tension, to the accuracy of A and ρ) and which do not (any statement about frequencies that matter).

Improvements worth making: quadratic elements (extension 1), which make the uniform spar's stress exact in every element and the strut converge at order 3; a graded mesh toward the strut's tip (extension 2); a stress-recovery step. Improvement worth *not* making: lumping the loads (half the element's body force to each node) to save the quadrature. It costs an order of accuracy in the displacement for a linear load and loses the nodal exactness.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (2), (4), the consistent mass matrix, the properties of K and the nodal-exactness argument; `bar_fem.png`; the four tables of §4; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the spar's root stress, tip extension and first axial frequency, each with the mesh used and the error that mesh leaves, and the assumption of the bar model you trust least.

### Adversarial self-checks

- **What if the rotor speed rises 10 %?** Every centrifugal load scales with Ω², so stresses, extension and reaction rise 21 %: root stress 107 MPa. The axial frequencies of the non-spinning bar do not change at all, which is the clearest sign that the axial model is not the one that governs the blade's dynamics.
- **What if the tip mass goes to zero?** The root stress falls to 59.2 MPa (the spar's own load) and the tip extension to 2.82 mm; nothing about the convergence or exactness changes. What if the section doubles? Stress halves, extension halves, reaction unchanged.
- **Which assumption fails first in the real world?** The fixed root at x = 0. A real blade's root fitting starts at a cut-out radius with pinned lugs; the nominal stress there is close to the bar's value but the peak at the pin holes is two to three times higher.
- **What measurement could falsify the model?** A strain gauge at mid-span on a whirl-tower run: the model predicts ε = σ/E = (ρΩ²(L² − x²)/2 + P/A)/E = 1.05 × 10⁻³ at x = 2.5 m. A reading 20 % off with the same Ω points to a wrong A, ρ or tip mass; a reading that varies with azimuth points to bending the bar model cannot see.
- **Does the result survive a change of mesh?** The nodal displacements and the reaction do not change at all with n (nodal exactness); the element stresses change by 1/n and the frequencies by 1/n². Report the first set as exact for the model and the second with its mesh.
- **Which decision changes if the uncertainty doubles?** If the alloy's allowable were known to ±20 % instead of ±10 %, the 1.4 % discretisation error of a 4-element root stress would still be irrelevant and the decision (the spar is at a third of yield) unchanged; a 1 % error on a frequency near a rotor harmonic would not be.

## Check yourself

1. Write the strong form and the weak form of the axially loaded bar and say which boundary condition is essential and which is natural.
2. Derive the element stiffness matrix of a two-node bar element from K_ij = ∫EAN_i'N_j' dx.
3. Why is the global stiffness matrix singular before any displacement is fixed, and what is its null vector?
4. State the consistent load vector for a uniform body force b₀ on an element of length h, and for b = b₁x/h measured from the element's first node.
5. Why is the finite-element solution of a uniform bar exact at the nodes, and why does that not hold for a tapered bar?
6. The L2 displacement error falls at order 2 and the energy-norm stress error at order 1. Explain both from the best-approximation property.
7. The consistent-mass frequencies converge to the exact ones from above. Why, and from which side do lumped-mass frequencies converge?
8. On 64 elements the tapered strut's last-element stress is 1.2 % below the exact tip stress. Is that a convergence failure? What would you do about it?
