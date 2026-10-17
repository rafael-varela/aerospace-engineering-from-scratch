# P01.05  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Convert θ_s to radians first; every error law in lesson §2.4 is in radians. The alternating bound is the first *omitted* term: list the sin terms θ, θ³/6, θ⁵/120, θ⁷/5040 at your angle and find the first below 10⁻⁸.
2. For (a) multiply the partial sum by r and subtract; for the logarithm tail, factor out fⁿ⁺¹/(n + 1) and bound what remains by a geometric series. For (b) the θ⁴ coefficient of sin is zero, so take n = 4 and the fifth derivative. For (c) C(−½, k) = (−½)(−3/2)…/k!; for the ratio test compute |cₖ₊₁ M²ᵏ⁺²/(cₖ M²ᵏ)|. For (d) i² = −1, i³ = −i, i⁴ = 1: sort the powers.
3. ln(1/(1 − f)) = ∫₀ᶠ dt/(1 − t); integrate Σtᵏ term by term. Impulse divided by the *initial* mass is v_e f.
4. The largest term of the sin series at θ is near k ≈ θ/2 and is about e^θ/√(2πθ) in size; machine epsilon is 2.2 × 10⁻¹⁶. Argument reduction uses the period 2π; e⁻ˣ = 1/eˣ exactly.
5. (1 + u)^(3.5) = 1 + 3.5u + 4.375u² + 2.1875u³ + … with u = 0.2M²; subtract 1 and divide by 0.7M². A binomial series in u converges for |u| < 1, and u = 0.2M² < 1 means M < √5.
6. Each statement misuses one of: the integral test, the radius of convergence (and the rate of convergence near it), the range of validity of cos θ ≈ 1 − θ²/2 (φ = 60° is 1.05 rad), the hypotheses of the alternating-series bound (decreasing terms).
7. (1 − k² sin²φ)^(−1/2) = Σ C(−½, n)(−k² sin²φ)ⁿ; the sign (−1)ⁿ cancels the one in C(−½, n). For the re-expansion, sin²(θ₀/2) = θ₀²/4 − θ₀⁴/48 + …. The inertia goes as T², so its relative error is about twice the period's.
