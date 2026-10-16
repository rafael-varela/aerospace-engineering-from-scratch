# P03.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. One Newton step by hand ($s = 15$: c = 2.5)

(a) J = [[2x, 2y], [y, x]]. At (1.5, 0.5): F = (2.25 + 0.25 − 2.5, 0.75 − 0.3) = **(0, 0.45)**, ‖F‖ = 0.45; **J = [[3, 1], [0.5, 1.5]]**, det J = 4. Solve Jδ = −F: 3δ₁ + δ₂ = 0, 0.5δ₁ + 1.5δ₂ = −0.45 → δ₂ = −3δ₁, 0.5δ₁ − 4.5δ₁ = −0.45, **δ = (0.1125, −0.3375)**, **x₁ = (1.6125, 0.1625)**, F(x₁) = (0.12656, −0.03797), **‖F₁‖ = 0.1321**.

(b) `newton` from (1.5, 0.5): residual norms 0.45, 0.132, 2.80 × 10⁻³, 1.75 × 10⁻⁶, 7.2 × 10⁻¹³, 4 × 10⁻¹⁶; iterates converge to **(1.569543, 0.191138)** in 5 steps; errors 0.317, 0.0516, 1.25 × 10⁻³, 8.0 × 10⁻⁷, 3.3 × 10⁻¹³; **observed order 2.00**.

(c) (x + y)² = x² + y² + 2xy = c + 0.6 = 3.1 and (x − y)² = c − 0.6 = 1.9, so x + y = ±1.76068, x − y = ±1.37840: **x = 1.569543, y = 0.191138** (and the other three sign combinations: (0.191, 1.570), (−1.570, −0.191), (−0.191, −1.570)). Newton found the root in whose basin it started. At (1, 1), det J = 2x² − 2y² = 0: the Jacobian is singular (the two level curves x² + y² = c and xy = 0.3 have parallel normals there) and the step is undefined; `np.linalg.solve` raises.

**Check.** F(1.569543, 0.191138) = (2.4636 + 0.0365 − 2.5, 0.3000 − 0.3) = (0, 0) to five digits.

**Common mistakes.** Forgetting the factor 2 in ∂(x²)/∂x; expecting Newton to find "the" root of a system with four; reading a `LinAlgError` as a bug rather than as det J = 0.

## 2. Convergence, damping, Broyden, and the arctan threshold

(a) Taylor with Lagrange remainder about x_k, evaluated at x*: 0 = f(x*) = f(x_k) + f'(x_k)(x* − x_k) + ½f''(ξ)(x* − x_k)². Divide by f'(x_k): 0 = f/f' − e_k + f''(ξ)e_k²/(2f'), and since x_{k+1} = x_k − f/f', e_{k+1} = e_k − f/f' = **f''(ξ)e_k²/(2f'(x_k))**, so |e_{k+1}| ≤ C e_k² with C = max|f''|/(2 min|f'|) in a neighbourhood where f' does not vanish. At a double root f'(x_k) ≈ f''(x*)e_k, so e_{k+1} ≈ e_k/2: Newton becomes linear with factor ½, the Babylonian iteration on x² has it.

(b) m = ½Σ_i F_i², so ∂m/∂x_j = Σ_i F_i ∂F_i/∂x_j = (JᵀF)_j: **∇m = JᵀF**. With δ = −J⁻¹F: ∇m·δ = −FᵀJJ⁻¹F = **−‖F‖² < 0** whenever F ≠ 0: a descent direction regardless of where x is. Expanding, F(x + λδ) = F + λJδ + O(λ²) = (1 − λ)F + O(λ²), so m(x + λδ) = ½(1 − λ)²‖F‖² + O(λ²) = **(1 − 2λ)m(x) + O(λ²)**, which lies below (1 − 2c₁λ)m(x) for all sufficiently small λ when c₁ < 1. Backtracking from λ = 1 by halving therefore always terminates.

(c) Write B = B_k + E. The secant condition Bs = y becomes Es = r := y − B_k s. For any matrix, ‖Es‖₂ ≤ ‖E‖_F‖s‖₂ (apply Cauchy–Schwarz to each row and sum squares), so every admissible E has ‖E‖_F ≥ ‖r‖₂/‖s‖₂. The rank-one E* = rsᵀ/(sᵀs) satisfies E*s = r and ‖E*‖_F = ‖r‖₂‖s‖₂/(sᵀs) = ‖r‖₂/‖s‖₂: it attains the lower bound, so **B_{k+1} = B_k + (y − B_k s)sᵀ/(sᵀs)** is the least change in the Frobenius norm. (It is also unique, since the bound is attained only when every row of E is parallel to sᵀ.)

(d) The Newton iterate is x₁ = x₀ − arctan(x₀)/(1/(1 + x₀²)) = x₀ − (1 + x₀²) arctan x₀. Setting x₁ = −x₀: **2x₀ = (1 + x₀²) arctan x₀**, which `newton` on g(x) = 2x − (1 + x²) arctan x from 1.5 solves as **x₀ = 1.39175**. The map x ↦ |x₁| is increasing in |x₀| (its derivative, 2|x| arctan|x| for |x| > 0 after simplification of the sign, is positive), so above the threshold |x₁| > |x₀| and each step is larger than the last; below it |x₁| < |x₀| and the iterates shrink toward the root, where (a) takes over with f'' = −2x/(1 + x²)², f' = 1/(1 + x²): C = 2x₀ at most, and quadratic convergence.

**Check.** The code's run: from 1.0 the iterates are 1, −0.571, 0.117, −1.06 × 10⁻³, 8 × 10⁻¹⁰ (quadratic from the second step); from 1.5 they are −1.69, 2.32, −5.11, 32.3, −1 575, 3.9 × 10⁶.

**Common mistakes.** In (a), dropping the remainder's ξ and writing f''(x*) (true only asymptotically); in (b), believing a descent direction guarantees convergence (it guarantees monotone decrease of m, which can stall at a local minimum of m where F ≠ 0); in (c), proving E* works without proving nothing smaller does.

## 3. A two-node radiator ($s = 15$: Q₁ = 35 W)

(a) **F₁ = Q₁ − G(T₁ − T₂) − εσA₁T₁⁴ = 0, F₂ = G(T₁ − T₂) − εσA₂T₂⁴ = 0**, with εσA₁ = 2.410 × 10⁻⁹ W/K⁴ and εσA₂ = 9.640 × 10⁻⁹ W/K⁴. **J = [[−G − 4εσA₁T₁³, G], [G, −G − 4εσA₂T₂³]]**; at (300, 300) K: [[−0.7602, 0.5], [0.5, −1.5411]] W/K, and `fd_jacobian` agrees to 2 × 10⁻⁸.

(b) Newton from (300, 300) K: **T₁ = 266.40 K, T₂ = 220.68 K**, 5 iterations, residual norms 79.6, 14.4, 1.39, 2.2 × 10⁻², 6.4 × 10⁻⁶, 5.6 × 10⁻¹³ W; errors 86.1, 28.9, 3.66, 6.3 × 10⁻², 1.9 × 10⁻⁵ K, **observed order 1.98**. Energy: radiated 12.14 + 22.86 = **35.00 W = Q₁**; conducted 22.86 W = the radiator's emission, as F₂ demands.

(c) At (1, 1) K the T⁴ terms are 10⁻⁹ W and their derivatives 10⁻⁸ W/K, so J ≈ 0.5[[−1, 1], [1, −1]] with det J = 2.4 × 10⁻⁸: nearly singular, because conduction alone cannot fix the temperature level (adding a constant to both T's changes nothing). The first step is **(7.3 × 10⁸, 7.3 × 10⁸) K**. From there the problem is essentially εσ(A₁ + A₂)T⁴ = Q₁, a quartic, for which Newton's step is T/4 − Q₁/(…): the iterate shrinks by 3/4 per step, and the undamped iteration needs **58 steps** (log(7.3 × 10⁸/266)/log(4/3) ≈ 52, plus the quadratic finish). `damped_newton` accepts only λ = 10⁻⁶ of that first step (landing at 693 K), then takes full steps: **9 iterations**, 30 residual evaluations.

(d) Broyden from (300, 300): 9 iterations, **12 evaluations** against Newton's 6 evaluations plus 5 Jacobians (or 16 evaluations with `fd_jacobian`). The isothermal node of area 0.25 m² sits at (35/(0.85 × 5.67 × 10⁻⁸ × 0.25))^{1/4} = **232.2 K**. The finite conductance needs a 45.7 K drop to carry 22.9 W to the radiator, so the box runs hotter (266.4 K) and the radiator, which receives only 22.9 of the 35 W, cooler (220.7 K) than the isothermal case. The weighted mean A₁T₁⁴ + A₂T₂⁴ is fixed by Q₁, so pulling T₂ down must push T₁ up.

**Check.** Dimensional: εσAT⁴ at 266 K and 0.05 m² = 12.1 W, consistent with a few hundred W/m² for a warm radiating surface.

**Common mistakes.** Using Celsius; dropping the factor 4 in the Jacobian (Newton then converges linearly); starting from (0, 0), where J is exactly singular.

## 4. A more stable aircraft and its ceiling ($s = 15$: C_mα = −1.5 rad⁻¹, T₀ = 85 kN)

(a)

| h (km) | α (°) | δe (°) | T (kN) | iterations |
|---|---|---|---|---|
| 0 | −0.785 | 3.700 | 52.55 | 2 |
| 5 | 0.713 | 2.042 | 35.31 | 2 |
| 10 | 3.637 | −1.093 | 26.96 | 2 |

At 5 km, δe = 2.042° against 2.218° for C_mα = −1.2: the more stable aircraft at the same positive α has a more nose-down aerodynamic moment, which is balanced by the thrust offset and C_m0, so the elevator has to add *less* nose-down (trailing-edge-down) deflection. α itself barely changes (0.713° against 0.700°) because the lift equation is nearly independent of C_mα. Thrust is unchanged to four digits.

(b)

| h (km) | T_required (kN) | T_avail = T₀ρ/ρ₀ (kN) |
|---|---|---|
| 0 | 52.55 | 85.00 |
| 5 | 35.31 | 51.08 |
| 8 | 29.22 | 36.44 |
| 10 | 26.96 | 28.64 |
| 11 | 26.40 | 25.25 |

They cross between 10 and 11 km.

(c) Appending F₄ = (T − T₀ρ(h)/ρ₀)/W and using `newton` with a differenced 4 × 4 Jacobian from (α, δe, T/W, h) = (5 km trim, 5): **h = 10.60 km, α = 4.149°, δe = −1.628°, T = 26.58 kN**, 5 iterations. ρ(h) is a power law of (1 − 0.0065h/288.15); its derivative is a few lines of algebra that would be written once, used once, and be the most likely place for a bug; the differenced Jacobian costs 5 evaluations of a trivially cheap function and is accurate to 10⁻⁸, which Newton does not notice.

(d) The code's `isa_density` is the troposphere formula; above 11 km the real atmosphere is isothermal and ρ falls exponentially, faster than the formula says, so a ceiling computed above 11 km would be too high. Check `h <= 11 km` on the converged value (10.60 km: valid here); for a heavier aircraft the check would fail, and the fix is the stratosphere branch in `isa_density`, not in the solver.

**Common mistakes.** Mixing h in metres into a vector whose other entries are O(0.1) (κ(J) jumps by 10⁴; use km); taking the ceiling as the altitude where T_avail = T_required at *minimum-drag* speed (the absolute ceiling) rather than at the specified 200 m/s; forgetting that T/W, not T, is the unknown.

## 5. The chord method ($s = 15$: β = −0.006 K⁻¹, n = 8)

(a)

| m | iterations | residual evaluations | observed order |
|---|---|---|---|
| 1 (Newton) | 6 | 391 | 1.94 |
| 2 | 6 | 199 | 1.94 |
| 3 | 8 | 201 | 1.24 |
| 5 | 9 | 138 | 1.53 |
| ∞ (chord) | 27 | 92 | 1.00 |

(b) **m = ∞** minimises the evaluations here: 92 against 391 for Newton. The Jacobian costs 64 evaluations and each iteration one; recomputing it buys quadratic convergence, but on a problem with a mild nonlinearity (error ratio 0.28 per chord step) the 21 extra cheap iterations cost less than one Jacobian. Had the ratio been 0.9, the chord method would have needed 200 iterations and m = 3–5 would have won; the optimum depends on how nonlinear the problem is, which the first few error ratios reveal.

(c) The chord iteration is x_{k+1} = g(x_k) with g(x) = x − J₀⁻¹F(x), a fixed-point iteration; at the root, g'(x*) = I − J₀⁻¹J(x*), and e_{k+1} ≈ g'(x*)e_k, so the asymptotic error ratio is **ρ(I − J₀⁻¹J(x*))**, computed as **0.298**. The measured ratios over the first five chord steps are 0.118, 0.242, 0.265, 0.275, 0.281, approaching the spectral radius from below as the error aligns with the dominant eigenvector. For m = 5 the ratios are the same until the Jacobian is refreshed, then the error drops by orders of magnitude.

(d) Picard to 10⁻⁹: **25 iterations**, each an assembly of A(T) and one 64 × 64 solve, 26 residual-sized evaluations in all; about the same iteration count as the chord method (27) with comparable cost, and no Jacobian code at all. At n = 64 (4 096 unknowns) a differenced Jacobian costs 4 097 assemblies per Newton step, which is prohibitive; Picard, or Newton with the analytic Jacobian of A(T)T (a five-point stencil plus the dk/dT terms), or Newton–Krylov with a differenced Jacobian–vector product, are the choices, and Picard is the one you can write in ten lines.

**Common mistakes.** Measuring the order from the m = 2 run and calling it quadratic (the refresh steps are quadratic, the others linear; the least-squares slope mixes them); counting the Jacobian as one evaluation.

## 6. Choosing the nonlinear solver

(a) Newton with the analytic Jacobian: nine derivatives by hand, 3 iterations, done.
(b) Newton from the previous frame's trim (continuation): one step per frame is usually enough because the state moved only slightly; the Jacobian can even be reused for several frames (chord).
(c) Broyden: one 10-minute evaluation per iteration; a differenced Jacobian would cost n + 1 runs and quadratic convergence is not worth an hour per step. Start B from a cheap model's Jacobian.
(d) Newton–Krylov with a preconditioner (incomplete factorisation or multigrid from `P03.09`), Jacobian-free; a dense or differenced Jacobian is out of the question at 10⁶.
(e) Pseudo-arclength continuation: at the fold the Jacobian with respect to the load is singular and plain Newton in displacement jumps or fails; parameterise by arc length along the curve.
(f) Damped Newton (the line search refuses steps that increase ‖F‖), or continuation from a solvable parameter value, or a change of variables (solve for ln T so negative temperatures cannot be represented).
(g) Picard: A(u) is cheap to assemble, the Jacobian of A(u)u is not available, and a mild nonlinearity gives a small contraction factor; switch to Newton only if the Picard ratios approach 1.
