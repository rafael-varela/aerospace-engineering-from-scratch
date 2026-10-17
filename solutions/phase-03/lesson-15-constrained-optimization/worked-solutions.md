# P03.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. KKT by hand ($s = 15$: right-hand side 3.5)

**Derivation.** g(x, y) = 3.5 − x − y ≥ 0, L = (x − 3)² + (y − 2)² − λ(3.5 − x − y). KKT: 2(x − 3) = −λ, 2(y − 2) = −λ (stationarity, ∇f = λ∇g with ∇g = (−1, −1)); 3.5 − x − y ≥ 0; λ ≥ 0; λ(3.5 − x − y) = 0.

(b) The unconstrained minimum (3, 2) has x + y = 5 > 3.5, so the constraint is active. From stationarity x − 3 = y − 2 = −λ/2; substituting into x + y = 3.5: 5 − λ = 3.5, **λ* = 1.5**, **(x*, y*) = (2.25, 1.25)**, f* = 2 × 0.75² = **1.125**. Shadow price: loosening the right-hand side by 0.1 changes f* by −λ × 0.1 = **−0.15** (exactly −0.1425 by recomputation: the derivative is exact only to first order).

(c) Inactive when 5 ≤ 2 + s/10, i.e. **s ≥ 30**: then (x*, y*) = (3, 2), λ* = 0, f* = 0.

**Check.** (2.25, 1.25) is the foot of the perpendicular from (3, 2) to the line x + y = 3.5: the displacement (−0.75, −0.75) is parallel to the normal (1, 1) ✓. Dimensions: λ has the units of f per unit of g.

**Common mistakes.** Writing the constraint as g = x + y − 3.5 ≥ 0 (the wrong side) and getting λ < 0; forgetting to test whether the unconstrained minimum is feasible before assuming the constraint is active.

## 2. Prices, Newton, penalty

(a) On the path of optima each active constraint satisfies g_i(x*(ε)) + ε = 0 and each inactive one has λ_i = 0, so f*(ε) = f(x*) = f(x*) − Σλ_i(g_i + ε) = L(x*, λ*; ε). Then df*/dε = ∇_xL · dx*/dε − Σ (g_i + ε) dλ_i/dε + ∂L/∂ε. The first term vanishes by **stationarity** (∇_xL = 0), the second by **complementarity** (g_i + ε = 0 wherever dλ_i/dε can be nonzero), and ∂L/∂ε = −λ_i for the perturbed constraint. Hence df*/dε = −λ_i*.

(b) The QP's Lagrangian is ½pᵀWp + ∇fᵀp − λᵀ(g_𝒲 + J_𝒲p). Its stationarity Wp + ∇f − J_𝒲ᵀλ = 0 together with the active linearised constraints g_𝒲 + J_𝒲p = 0 is [[W, −J_𝒲ᵀ], [J_𝒲, 0]][p; λ] = [−∇f; −g_𝒲]. Newton on F(x, λ) = (∇f − J_𝒲ᵀλ, g_𝒲) has Jacobian [[∇²f − Σλ_i∇²g_i, −J_𝒲ᵀ], [J_𝒲, 0]] = the same matrix with W the Lagrangian Hessian, and right-hand side −F = (−∇f + J_𝒲ᵀλ_old, −g_𝒲); writing λ_new = λ_old + Δλ turns it into the QP system. The block matrix is nonsingular iff J_𝒲 has full row rank (independent active gradients, LICQ) and W is positive definite on the null space of J_𝒲.

(c) ∇φ = ∇f + μ Σ_i min(g_i, 0)∇g_i = 0. On a violated constraint min(g_i, 0) = g_i < 0, so ∇f = Σ μ(−g_i)∇g_i; matching with ∇f = Σλ_i∇g_i gives **μ(−g_i) → λ_i**, a violation of λ_i/μ. Prediction for the lesson: 0.464/100 = **4.6 × 10⁻³**; the run gives 4.60 × 10⁻³ ✓.

(d) With u = ln r, v = ln t: minimise u + v subject to u + v ≥ ln a_s, 3u + v ≥ ln a_e, 2v ≥ ln a_l, v ≥ ln t_min, u ≥ ln r_min: a linear program. At the Euler + local vertex (1, 1) = λ'_E(3, 1) + λ'_l(0, 2) gives 3λ'_E = 1 and λ'_E + 2λ'_l = 1, so **λ'_E = λ'_l = 1/3**. Near activity g_E = 1 − a_e/(r³t) ≈ ln(r³t/a_e), so dg_E = 3du + dv (the same gradient as the log constraint), and dm = m(du + dv); the multiplier that prices dm in kilograms is m* times the one that prices d ln m: **λ_E = λ_l = m*/3 = 1.3920/3 = 0.46401 kg**, matching the code's least-squares multipliers.

**Common mistakes.** In (a), differentiating f directly and forgetting that x* moves; in (c), using the sign convention g ≤ 0; in (d), forgetting the factor m* when converting.

## 3. A propellant tank ($s = 15$: V = 3.5 m³)

**Assumptions.** Thin wall, uniform thickness, hemispherical ends, hoop stress pR/t in the cylinder governs (the hemisphere's pR/(2t) is half of it), no end fittings or welds.

(a) L = (V − (4/3)πR³)/(πR²), so m = ρt(2πRL + 4πR²) = **ρt(2V/R + (4/3)πR²)**. Constraints: g₁ = 1 − pR/(σ_a t), g₂ = t/t_min − 1, g₃ = 1 − R/R_max, g₄ = 1 − (4/3)πR³/V, all ≥ 0.

(b) At fixed t, dm/dR = ρt(−2V/R² + (8/3)πR) = 0 gives R³ = 3V/(4π): the volume of a sphere of radius R equals V, so L = 0. A **sphere** has the least surface for its volume. For V = 3.5 m³, R_sphere = **0.942 m** > 0.6 m: it does not fit.

(c) The radius limit is active, R = **0.6 m**, so L = (3.5 − 0.9048)/1.1310 = **2.295 m**. Hoop stress at the gauge: pR/t_min = 0.3 × 10⁶ × 0.6/1.2 × 10⁻³ = 150 MPa < 180 MPa, so the stress constraint is inactive and the **gauge is active**, t = **1.2 mm**. Mass m = 2700 × 1.2 × 10⁻³ × (2 × 3.5/0.6 + (4/3)π × 0.36) = 3.24 × (11.667 + 1.508) = **42.69 kg**. Multiplier of the radius limit: λ_R = −dm/dR = −3.24 × (−19.44 + 5.03) = **46.7 kg/m**: a fairing 10 mm wider would save about **0.47 kg**. The gauge multiplier is ∂m/∂t × t_min = 35.57 kg per mm of gauge (`sqp` reports 42.69 per unit of the normalised constraint t/t_min − 1).

(d) With x = (R/0.6, t/1.2 mm) and the four normalised constraints, `sqp` from (0.3 m, 3.6 mm) converges in **2 iterations** (the problem is nearly linear in the scaled variables) to R = 0.6000 m, t = 1.2000 mm, m = 42.6858 kg, multipliers (0, 42.69, 28.03, 0) in normalised units; 28.03/0.6 = 46.7 kg/m ✓. SLSQP: the same point to 10⁻¹⁰ in 10 iterations.

**Check.** The sphere limit: with R_max = 1 m the optimum would be the sphere, m = ρ t_min 4πR² = 3.24 × 4π × 0.887 = 36.1 kg, lighter than 42.7 kg as it must be, since a relaxed constraint cannot raise the optimum.

**Common mistakes.** Keeping L as a variable with the volume as an equality constraint in an inequality-only solver; forgetting that the hemispherical ends also use material; reporting the multiplier in normalised units and calling it kg/m.

## 4. Your own strut ($s = 15$: P = 45 kN, L = 2.75 m)

(a) Gauge threshold 2πkE t_min² = 33.1 kN < 45 kN, so the gauge is inactive. t* = √(45 000/(2π × 0.3058 × 68.9 × 10⁹)) = **0.5830 mm**; a_e = PL²/(π³E) = 45 000 × 7.5625/(31.006 × 68.9 × 10⁹) = 1.593 × 10⁻⁷ m⁴; r* = (a_e/t*)^{1/3} = **64.89 mm**; stress 45 000/(2π × 0.06489 × 5.830 × 10⁻⁴) = **189.3 MPa** < 240 MPa, so the active pair is **Euler + local** and m* = 2π × 2700 × 2.75 × 0.06489 × 5.830 × 10⁻⁴ = **1.7650 kg**.

(b) and (d):

| quantity | hand | SQP-BFGS | SQP-exact | SLSQP |
|---|---|---|---|---|
| iterations | — | 10 | 8 | 8 |
| final KKT residual | 0 | < 10⁻¹² | < 10⁻¹² | — |
| ‖x − x*‖ (scaled) | 0 | 1.6 × 10⁻¹⁰ (then rounding) | 3.3 × 10⁻⁹ (then rounding) | 1.8 × 10⁻¹⁵ |
| λ_Euler, λ_local [kg] | m*/3 = 0.58832 | 0.58832, 0.58832 | 0.58832, 0.58832 | — |
| observed order, last three errors | 2 | 2.01 | 2.01 | — |

The BFGS errors end 3.6 × 10⁻², 1.8 × 10⁻³, 8.0 × 10⁻⁶, 1.6 × 10⁻¹⁰; the exact-Hessian ones 3.7 × 10⁻³, 3.7 × 10⁻⁵, 3.3 × 10⁻⁹.

(c) d m*/dP = (λ_E + λ_l)/P = 1.17664/45 000 = **2.6148 × 10⁻⁵ kg/N** = 26.1 g/kN; `shadow_price_check` gives 2.61476 × 10⁻⁵ (prediction) and 2.61476 × 10⁻⁵ (finite difference). Every kilonewton of load removed saves 26 g, to first order.

**Check.** (2/3)m*/P = (2/3) × 1.7650/45 000 = 2.6148 × 10⁻⁵ ✓ (m* ∝ P^{2/3} on this vertex).

**Common mistakes.** Using the unscaled history to compute the order (the numbers change but the order does not; mixing the two does); keeping errors below 10⁻¹² in the order estimate, where rounding gives nonsense.

## 5. The price of a penalty (default strut, start (50 mm, 2 mm))

(a)

| μ | ‖x_μ − x*‖ | max violation | objective evaluations |
|---|---|---|---|
| 10 | 1.42 × 10⁻² | 4.32 × 10⁻² | 64 |
| 10² | 1.56 × 10⁻³ | 4.60 × 10⁻³ | 47 |
| 10³ | 1.57 × 10⁻⁴ | 4.64 × 10⁻⁴ | 66 |
| 10⁴ | 1.57 × 10⁻⁵ | 4.64 × 10⁻⁵ | 89 |
| 10⁵ | 1.57 × 10⁻⁶ | 4.64 × 10⁻⁶ | 185 |
| 10⁶ | 1.57 × 10⁻⁷ | 4.64 × 10⁻⁷ | 243 |

Fitted order in 1/μ over 10²–10⁴: **0.987**; over the whole range essentially 1. The violation is λ/μ = 0.464/μ to three figures from μ = 10³ on. The evaluation count grows because the penalty Hessian W + μ Σ ∇g_i∇g_iᵀ has two eigenvalues of order μ (two active constraints in two variables here, so the growth is mild) and BFGS's line search needs more trials per step as the valley narrows; with more variables than active constraints the condition number grows as μ and the count grows faster.

(b) Augmented Lagrangian: μ = 1: **27** outer iterations, error 1.1 × 10⁻⁹, 399 evaluations; μ = 10: **7**, 1.9 × 10⁻⁹, 156; μ = 100: **4**, 7.1 × 10⁻¹⁰, 77. Each outer iteration shrinks the multiplier error by a factor of order 1/μ (the linear rate of the method of multipliers), so μ = 100 needs fewer outer passes; at μ = 1 the first subproblem is barely penalised and the inner BFGS wanders toward thin walls before the multipliers build up.

(c) SQP: **23** objective evaluations in 9 iterations. For a violation of 10⁻⁸ the penalty needs 0.464/μ ≤ 10⁻⁸, μ ≥ 4.6 × 10⁷; extrapolating the roughly ×1.3–2 growth per decade of (a) gives of order 500–1000 evaluations, twenty to forty times SQP's, for a point that is still infeasible.

**Common mistakes.** Fitting the order including μ = 10, where the asymptotic law has not set in; comparing evaluation counts without the same stopping tolerance.

## 6. Reading an optimiser's output

(a) **Two scaling problems.** Variables in metres are O(10⁻³) while BFGS starts from the identity, so its first steps are of order 1 m: the line search must shrink them by a factor of 10³ before anything sensible happens, and the Hessian approximation spends its first dozen updates learning the scale. Constraints in pascals are O(10⁸) while the objective is O(10²) kg: the l₁ merit function f + μΣ max(0, −g) is dominated by any constraint violation unless μ is tiny, and the QP's multipliers are O(10⁻⁶) kg/Pa, below the solver's tolerances. The "positive directional derivative" message is the typical symptom: the QP step is not a descent direction for the merit function at the precision available.

(b) 3 × 10⁴ Pa on a 300 MPa limit is a violation of **10⁻⁴** of the allowable: within any realistic knowledge of the allowable itself, and far below the 10 % imperfection knock-down on a buckling limit. With g = 1 − σ/300 MPa the output would have read g = −10⁻⁴ and the question would not arise.

(c) A multiplier of 2 × 10⁻⁹ kg/m on a 400 kg wing box is a relative sensitivity of 5 × 10⁻¹² per metre: the bound is active only in the sense that the optimiser landed on it, and relaxing it by a metre would save 2 ng. Do not relax it; check instead that the constraint is not degenerate (two constraints coinciding). Compare with the lesson's strut, where each active multiplier is one third of the mass.

(d) With both active constraints being buckling limits with ±30 % knock-down uncertainty, the lesson's law m* ∝ γ^{−1/3} gives a mass uncertainty of about (1.3^{1/3} − 1) ≈ **+9 %**/−11 % (γ = 0.7 raises m* by 12.6 %; γ = 1.3 lowers it by 8.4 %): 412.3 kg is 410 ± 40 kg, and only the first two figures are defensible. The quoted precision is that of the arithmetic, not of the physics.

## Check-yourself answers

1. L = f − Σλ_i g_i; stationarity ∇f = Σλ_i∇g_i, primal feasibility g ≥ 0, dual feasibility λ ≥ 0, complementarity λ_i g_i = 0. For an inactive constraint (g_i > 0) complementarity forces λ_i = 0: it does not influence the optimum at all.
2. Δm* ≈ −λ Δg = −0.46 × 0.01 = −4.6 g: the mass *decreases* when a binding limit is loosened.
3. Because the constraints are curved: moving along a curved boundary changes f at second order through the constraint's curvature, and only ∇²f − Σλ_i∇²g_i captures it. The strut's mass Hessian [[0, c], [c, 0]] is indefinite, yet the Lagrangian Hessian at the optimum has eigenvalues 2.13 and 14.06, positive, because the buckling constraints are convex and their multipliers are 0.464 kg.
4. The penalty only acts once a constraint is violated, so the minimiser balances the objective's pull against a finite restoring force; the violation of an active constraint is ≈ λ_i/μ, first order in 1/μ.
5. Positive definiteness of the BFGS matrix B, which the QP needs; it is triggered when sᵀy < 0.2 sᵀBs, that is when the Lagrangian shows little or negative curvature along the step.
6. Raising the allowable stress by 20 % changes nothing: the stress constraint is inactive (λ = 0) and stays inactive (194 MPa < 288 MPa). Raising the gauge to 0.6 mm makes it active (t* = 0.55 mm < 0.6 mm): the design moves to Euler + gauge, r = 58.0 mm, t = 0.60 mm, m = 1.476 kg, 6 % heavier.
7. Quadratic (order 2) with the exact Hessian, superlinear (order between 1 and 2, in practice measured near 2 on the last steps) with BFGS; q = ln(e₃/e₂)/ln(e₂/e₁) from three successive errors above the rounding floor.
8. Because in u = ln r, v = ln t the objective and all constraints are linear, so the feasible set is a convex polygon and the objective has no local minima other than the global one at a vertex; the vertex enumeration finds it exactly.
