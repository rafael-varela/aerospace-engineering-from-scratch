# Phase 01 Exam: Mathematics I: Calculus & Linear Algebra

Integrated assessment of the whole phase. Suggested self-assessment mode: closed reference, a time limit chosen and recorded in advance (the exam has not been piloted, so no duration is suggested), no code except in Part C, no computer algebra in Part B. There is no enforced lock; the point is to find out what you can do unaided. Solutions and pass thresholds are in `solutions/phase-01/exam/`.

Report your result per competency: (1) mathematical and physical reasoning, (2) analytical problem solving, (3) computational competence, (4) engineering interpretation. Questions are tagged with the competencies they assess. Every question is written over an aerospace situation; the mathematics is what is being tested, the framing is there so that you also practise deciding which mathematics applies. Angles are in radians unless a degree sign is written. SI units throughout.

## Part A: Concepts

**1.** (1, 4) A colleague reports that their forward-difference derivative "matches the SymPy derivative exactly, so I am going to shrink the step to h = 10⁻¹² to make it even more exact." Explain, in six sentences or fewer, why "exactly" cannot be literally true and why h = 10⁻¹² will make the result *worse*. Name the two error sources, say how each one scales with h, and give the order of magnitude of the step that balances them in double precision for a forward difference and for a central difference.

**2.** (1) An orbit in polar form is r(θ) = p / (1 + e cos θ). State what shape the curve has for e = 0, 0 < e < 1, e = 1 and e > 1, and at which θ the radius is smallest. Then show in two lines why the vector h = r × v is constant in time when the acceleration is directed along −r (a central force), and say what that constancy implies for the rate at which the radius vector sweeps area.

**3.** (4) A planar pin-jointed truss has j = 4 joints, m = 4 members and r = 4 reaction components, so the counting rule m + r = 2j is satisfied and a textbook would call it statically determinate. You assemble the 8 × 8 equilibrium matrix A (unknowns: 4 member forces and 4 reactions) and find rank A = 7. State what the nullspace of A contains and what it means physically, state what the left nullspace (A^T y = 0) tells you about which loads the structure can carry, and say what the designer should look for in the drawing.

**4.** (1) A direction cosine matrix C that rotates vector components from an inertial frame to a body frame satisfies two algebraic properties. State them. Say what a computed C with det C = −1 would represent and why it can never be a physical attitude. Finally, explain in two sentences why a C obtained by multiplying thousands of small incremental rotations in floating point stops satisfying the properties, and name one method from lessons 17 or 18 that restores them.

## Part B: Analysis

Work every derivation by hand and show the steps. Numbers to three significant figures unless more are needed to show a check.

**5.** (1, 2) Limits and the derivative.
(a) Using only the definition of the derivative and the angle-sum identity, show that d(sin α)/dα = cos α. You may use lim_{h→0} (sin h)/h = 1 and lim_{h→0} (1 − cos h)/h = 0.
(b) Hence state the lift-curve slope dC_L/dα at α = 0 for the thin-aerofoil model C_L = 2π sin α.
(c) Evaluate lim_{α→0} (1 − cos α)/α² without L'Hôpital's rule, and check your answer by evaluating the quotient at α = 0.01.

**6.** (2) A glider has the drag polar C_D = C_D0 + K C_L² with C_D0 = 0.020 and K = 0.050. Derive, by differentiating L/D = C_L / C_D with respect to C_L, the lift coefficient that maximises L/D, show that it is a maximum, and evaluate (L/D)_max. State the glide range from 1000 m height in still air. Show also that at this condition the drag is split equally between the C_D0 term and the K C_L² term. (A glider of mass 300 kg and wing area 12 m² flying this C_L at sea level: what airspeed? Use ρ = 1.225 kg/m³, g = 9.81 m/s².)

**7.** (2) Integrals.
(a) A wing of span b carries an elliptic lift distribution L'(y) = L₀ √(1 − (2y/b)²) per unit span, for −b/2 ≤ y ≤ b/2. Derive the total lift ∫L' dy (a trigonometric substitution does it) and the bending moment at the root from one half-wing, M = ∫₀^{b/2} y L'(y) dy (a u-substitution does it). Evaluate both for b = 10 m and L₀ = 400 N/m, and say what mass of aircraft this lift supports in level flight.
(b) A thin triangular fin of uniform areal density σ = 4 kg/m² has vertices at (0, 0), (c, 0) and (c, s), with c = 1.2 m and s = 0.6 m. Set up the double integrals for the area and the first moments and derive the centroid (x̄, ȳ) in terms of c and s; then give the numbers and the fin mass.

**8.** (1, 2) Taylor remainder. Write the Maclaurin series of sin α to the α³ term. Using the Lagrange form of the remainder, bound the absolute error of the small-angle approximation sin α ≈ α for |α| ≤ 10°, and compare the bound with the actual error at 10°. Give the relative error at 10°. Then find the angle at which the bound first allows a 1% relative error. Do the same bound for cos α ≈ 1 − α²/2 at 10° in one line.

**9.** (2, 4) Gradient and sensitivity. The lift of a wing in the linear range is L(α, V) = ½ ρ V² S C_Lα α. Take ρ = 1.225 kg/m³, S = 16 m², C_Lα = 5.7 per radian, α = 4°, V = 50 m/s.
(a) Derive the gradient ∇L = (∂L/∂α, ∂L/∂V) and evaluate it, giving ∂L/∂α both per radian and per degree.
(b) Show that the relative differential is dL/L = dα/α + 2 dV/V. A 2% increase in α and a 2% increase in V: which changes lift more, and by what ratio?
(c) Both increase by 2% together. Give the linearised relative change in L and the exact relative change, and comment on the difference.

**10.** (1, 2) Circulation. A two-dimensional point vortex at the origin has velocity field u = −Γy / (2π(x² + y²)), v = Γx / (2π(x² + y²)), with Γ = 12 m²/s.
(a) Show that the speed on a circle of radius R is Γ/(2πR) and give its value at R = 0.5 m.
(b) Compute the line integral ∮ (u dx + v dy) around the circle of radius R by parametrising it, and show the result does not depend on R.
(c) Show that ∂v/∂x − ∂u/∂y = 0 everywhere except the origin. Reconcile this with Green's (or Stokes') theorem and your answer to (b), and state without integrating what the circulation is around the square 1 ≤ x ≤ 2, 1 ≤ y ≤ 2.

**11.** (2) Elimination. Three unknown forces F₁, F₂, F₃ (kN) in a strut assembly satisfy the two force balances and one moment balance

```
 2F₁ +  F₂ −  F₃ =   8
−3F₁ −  F₂ + 2F₃ = −11
−2F₁ +  F₂ + 2F₃ =  −3
```

Solve by Gaussian elimination without row exchanges, showing the multipliers and the three pivots, then back-substitute. State the determinant from the pivots. Verify by substituting into the original equations.

**12.** (1, 2) Eigenvalues. A flat satellite panel has, in body axes, the planar inertia tensor

```
I = [ 5  −2 ]  kg·m²
    [−2   2 ]
```

Find the principal moments of inertia and the principal axes (unit eigenvectors). Check with the trace and the determinant. Show that the two eigenvectors are orthogonal and state why that was guaranteed. Give the angle of the axis of maximum inertia measured from the body x-axis.

**13.** (2, 3) Least squares. Four wind-tunnel points for a wing: C_L = 0.2, 0.4, 0.6, 0.8 with C_D = 0.023, 0.027, 0.039, 0.051. Fit C_D = C_D0 + K C_L² by least squares. Write the design matrix with columns (1, C_L²), form the normal equations A^T A θ = A^T b by hand (show the four sums), solve the 2 × 2 system, and report C_D0 and K. Compute the residual vector and show that it is orthogonal to both columns. Give (L/D)_max implied by the fit.

## Part C: Computation

**14.** (3) Write a Python function `power_iteration(A, n_iter=200)` that returns the dominant eigenvalue and a unit eigenvector of a real symmetric matrix `A` by repeated multiplication and normalisation, with the eigenvalue taken as the Rayleigh quotient at the end. Write three pytest tests for it: one that compares against `numpy.linalg.eigh` on a 3 × 3 inertia tensor, one that asserts the eigen-residual ‖A v − λ v‖ is small, and one limiting case (the identity matrix). Then explain in two sentences why the tests assert a residual and an oracle comparison rather than a hard-coded eigenvalue, and in one sentence what property of the eigenvalues sets the convergence rate.

## Part D: Interpretation

**15.** (4) A 2-D frame model gives a stiffness matrix K with condition number κ₂(K) = 5 × 10⁴. The load vector comes from strain-gauge-calibrated actuators with a relative uncertainty of 10⁻⁴. A colleague reports the solved displacements to six significant figures and says "the solver is double precision so the numbers are good to 15 digits."
(a) Using the bound ‖δx‖/‖x‖ ≤ κ ‖δb‖/‖b‖, state the worst-case relative error in the displacements due to the load uncertainty, and the worst-case contribution from double-precision rounding (ε ≈ 2.2 × 10⁻¹⁶). Which dominates, and is six significant figures justified?
(b) The bound is a worst case. Say what direction of load vector attains it and what that direction means physically for a structure (what the smallest singular value of K is telling you).
(c) Name two common causes of a large κ in a stiffness matrix that are *modelling* choices rather than physics, and say what you would check first.
