# P01.05  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. θ = (5 + s/2)° × π/180. 2. Errors θ²/6, θ²/2, θ⁴/24, θ²/3. 3. Compare each with 0.01, or compare θ with the 1 % angles 13.99°, 8.07°, 37.9°, 9.91°. 4. Terms θ^(2k+1)/(2k+1)! until one is below 10⁻⁸; the count of terms *kept* is the index of that first small one. 5. `small_angle_error(theta, a)` for each `a`; `terms_needed(theta, "sin", 1e-8)`.

## 2.
(a) Sₙ − rSₙ = 1 − rⁿ⁺¹. Tail after n + 1 terms: rⁿ⁺¹ + rⁿ⁺² + … = rⁿ⁺¹/(1 − r). For Σfᵏ/k: tail = Σ_{k>n} fᵏ/k ≤ (1/(n + 1)) Σ_{k>n} fᵏ = fⁿ⁺¹/((n + 1)(1 − f)). (b) f = sin, x₀ = 0, polynomial θ − θ³/6 = T₃ = T₄; use n = 4: |R₄| ≤ max|sin⁽⁵⁾| θ⁵/5! = θ⁵/120. Evaluate at θ = π/6 and compare with |θ − θ³/6 − ½|. (c) C(−½, k)(−M²)ᵏ for k = 0, 1, 2: 1, ½M², ⅜M⁴. Ratio of consecutive terms (2k + 1)M²/(2k + 2) → M², so convergence needs M² < 1. (d) e^(iθ) = Σ(iθ)ᵏ/k!; k even gives (−1)^(k/2)θᵏ/k! (the cos series), k odd gives i(−1)^((k−1)/2)θᵏ/k! (i times the sin series). Multiply e^(iθ) by e^(−iθ) = cos θ − i sin θ.

## 3.
(a) 1/(1 − t) = Σtᵏ for |t| < 1; integrate 0 to f: Σfᵏ⁺¹/(k + 1) = Σ_{k≥1} fᵏ/k, R = 1. (b) `v_e * sum(f**k / k for k in range(1, n + 1))` against `-v_e * np.log(1 - f)`; count terms until within 1 %. (c) Thrust × burn time / initial mass; the mass falls during the burn, so the true acceleration is higher. (d) Bound fⁿ⁺¹/((n + 1)(1 − f)) at f = 0.5, n = 5; actual tail = closed form − partial sum.

## 4.
(a) Loop over θ; `series_terms(theta, "sin", 120)` and `np.abs(...).max()`. Achievable absolute accuracy ≈ 2.2 × 10⁻¹⁶ × (largest term). (b) `maclaurin_partial_sum(40.0, "sin", n)` for n = 60, 80, 120; all three differ from `np.sin(40.0)` by a number of order eps × 1.5 × 10¹⁶. (c) θ′ = 40 − 12 × 2π; `maclaurin_partial_sum(theta_r, "sin", 12)`; `terms_needed(theta_r, "sin", 1e-15)`. (d) `maclaurin_partial_sum(-40.0, "exp", 120)` against `1 / maclaurin_partial_sum(40.0, "exp", 120)` and `np.exp(-40.0)`; for the positive argument all terms are positive and there is no cancellation.

## 5.
(a) u = 0.2M²: (1 + u)^3.5 − 1 = 3.5u + 4.375u² + 2.1875u³ + 0.2734u⁴ + …; divide by 0.7M² = 3.5u: 1 + 1.25u + 0.625u² + 0.078125u³ + …; substitute u. (b) `binomial_coefficients(3.5, 4)` → [1, 3.5, 4.375, 2.1875, 0.2734375]. (c) Exact ratio ((1 + 0.2M²)^3.5 − 1)/(0.7M²); pressure error 1 − q/q_c; airspeed error 1 − √(q/q_c). (d) Binomial series in u converges for |u| < 1, i.e. M < √5 ≈ 2.24; the isentropic pitot formula itself is wrong above M = 1 because a shock forms ahead of the probe (Rayleigh pitot formula, Phase 05).

## 6.
(a) Integral test: Hₙ ≈ ln n + 0.577, unbounded; the per-term change 1/n is the derivative of ln n. (b) Three terms at M = 0.9: compare with the exact factor 2.294; the series converges but needs 16 terms for 1 %; use the closed form. (c) φ = 60° = 1.047 rad is outside the 1 % range of cos ≈ 1 − φ²/2 (37.9°); use cos φ exactly: n = 2. (d) The alternating bound needs decreasing terms; at θ = 10 the terms grow until k = 4 (θ⁹/9!); the bound applies only after the maximum; use `terms_needed` or reduce the argument.

## 7.
(a) (1 − k² sin²φ)^(−1/2) = Σ C(−½, n)(−1)ⁿ k²ⁿ sin²ⁿφ = Σ [(2n)!/(4ⁿ n!²)] k²ⁿ sin²ⁿφ; integrate with the given formula: K = (π/2) Σ [(2n)!/(4ⁿ n!²)]² k²ⁿ; T/T₀ = (2/π)K. Re-expand 1 + k²/4 with k² = sin²(θ₀/2) ≈ θ₀²/4: 1 + θ₀²/16. (b) `pendulum_period_ratio_exact(np.radians(deg))`; inertia error 1 − 1/ratio²; bisection or `scipy.optimize.brentq` on ratio² − 1.01 and − 1.02. (c) θ₀ = 10 + s/2 degrees; error and correction (1/ratio²); `pendulum_period_ratio(theta0, n)` for n = 1, 2, 3 against exact. (d) Timing error in I is 2 × 0.2 % = 0.4 %; compare with the amplitude error; correct if the amplitude error exceeds the timing error and is known to better than that.
