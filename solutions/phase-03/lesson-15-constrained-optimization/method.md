# P03.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. g(x, y) = 2 + s/10 − x − y ≥ 0; L = f − λg; KKT: ∇f = λ∇g, g ≥ 0, λ ≥ 0, λg = 0. 2. Test (3, 2): if 5 ≤ 2 + s/10 the constraint is inactive. Otherwise set g = 0 and ∇f = λ(−1, −1): 2(x − 3) = 2(y − 2) = −λ, so x − 3 = y − 2 = −λ/2; substitute into x + y = 2 + s/10 to get λ, then x*, y*, f*. 3. Shadow price: Δf* ≈ −λ × Δ(right-hand side). For (c) solve 2 + s/10 ≥ 5.

## 2.
(a) f*(ε) = f(x*(ε)) = L(x*, λ*; ε) because λ_i(g_i + ε) = 0 on the path. Differentiate: dL/dε = ∇_x L · dx*/dε + ∂L/∂λ · dλ*/dε + ∂L/∂ε; the first term is zero by stationarity, the second by complementarity (either λ_i = 0 or g_i + ε = 0 stays zero), leaving ∂L/∂ε = −λ_i. (b) QP Lagrangian ½pᵀWp + ∇fᵀp − λᵀ(g + Jp); ∇_p = Wp + ∇f − Jᵀλ = 0 and g_𝒲 + J_𝒲 p = 0; stack as a block system and compare with Newton on F(x, λ) = (∇f − J_𝒲ᵀλ, g_𝒲). Nonsingular iff W is positive definite on the null space of J_𝒲 and J_𝒲 has full row rank (LICQ). (c) ∇φ = ∇f − Σ μ(−g_i)∇g_i for violated constraints; compare with ∇f = Σλ_i∇g_i: μ(−g_i) ≈ λ_i. (d) Constraints a_i u + b_i v ≥ ln C_i, objective u + v: an LP. At the vertex, (1, 1) = λ'_E(3, 1) + λ'_l(0, 2). Convert: the physical constraint near activity is g_i ≈ ln(r^a t^b/C_i), so dg_i = a du + b dv, and dm = m(du + dv); hence λ_i = m* λ'_i.

## 3.
(a) L = (V − (4/3)πR³)/(πR²) ≥ 0 gives (4/3)πR³ ≤ V. Substitute into m: m = ρt(2V/R + (4/3)πR²). Constraints: 1 − pR/(σ_a t) ≥ 0, t/t_min − 1 ≥ 0, 1 − R/R_max ≥ 0, 1 − (4/3)πR³/V ≥ 0. (b) With t = t_min, dm/dR = ρt(−2V/R² + (8/3)πR) = 0 gives R³ = 3V/(4π), i.e. L = 0; compare with R_max. (c) If the sphere does not fit, R = R_max; check whether pR/σ_a exceeds t_min to decide between gauge and hoop stress; compute L and m. Multiplier of the radius limit: λ_R = −dm/dR × R_max (normalised units) or −dm/dR in kg/m. (d) Scale x = (R/R_max, t/t_min); code f, grad, cons, jac; call `sqp` and SLSQP from a feasible start such as (0.3 m, 3.6 mm).

## 4.
(a) Gauge threshold 2πkE t_min² is independent of L; compute t* = √(P/(2πkE)), r* = (PL²/(π³E t*))^{1/3}, stress P/(2πr*t*) against 240 MPa. m* = 2πρL r*t*. (b) `res = sqp(prob, x0, hessian=...)`; errors from `res.history` against `analytic_optimum(params)[0]`; `observed_order` on the last three above 10⁻¹². (c) d m*/dP = (λ_E + λ_l)/P = (2/3)m*/P; `shadow_price_check(params)`; multiply by 1000 for kg/kN. (d) `minimize(..., method="SLSQP", options={"ftol": 1e-14})`.

## 5.
(a) For each μ: `penalty_minimize(counted_problem, mu, x0)`, error against x*, violation `-min(min(cons, 0))`, count. Fit `observed_order_from_resolutions(errors, 1/mus)`. Conditioning: penalty Hessian ≈ W + μ Σ ∇g_i∇g_iᵀ over violated constraints. (b) `augmented_lagrangian(counted_problem, x0, mu=...)` returns the outer count; relate the per-iteration contraction of the multiplier error to 1/μ. (c) Count SQP evaluations; the violation is ≈ λ/μ, so solve 0.464/μ ≤ 10⁻⁸; extrapolate the evaluation count from the trend of (a).

## 6.
(a) Variable scale (metres, 10⁻³) and constraint scale (pascals, 10⁸): identity initial Hessian wrong by many decades; merit function adds kilograms to pascals. (b) Divide the violation by the limit. (c) Compare the multiplier with the objective's magnitude and with a finite-difference re-optimisation; a multiplier at rounding level means the constraint is active only by coincidence or weakly. (d) Propagate ±30 % in the knock-down through m* ∝ γ^{−1/3}.
