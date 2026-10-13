# P01.04  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything is in (1) and the two formulas beside it: mass σ c_r s(1 + τ), centroid s(1 + 2τ)/(3(1 + τ)). A rectangular wing of the same mass has chord c̄ = (mass)/(2σ s); its inertia is 2σ c̄ s³/3.
2. For (a), always differentiate the power of y and integrate the exponential; two rounds remove y². For (b), with 0 ≤ θ ≤ π/2 the cosine is non-negative, which is why the square root is cos θ and not |cos θ|. For (c), expand ln((1 + x)/(1 − x)) = 2(x + x³/3 + …) with x = v/v_t. For (d), the midpoint error is +(h²/24)Δf′ and the trapezoid's is −(h²/12)Δf′; the weights 2 and 1 cancel them.
3. The transform of 1 is the simplest improper integral there is; the others follow from parts, and sin ωt from two rounds of parts that bring back the original integral (solve for it). Convergence needs the integrand to decay: compare the growth rate of f with s.
4. For (b), the θ-form integrand is (1 − cos 4θ)/8, and [0, π/2] is exactly one period of cos 4θ; an equal-weight rule on a full period integrates the k-th harmonic exactly unless k is a multiple of n, and here k = 1. Simpson's unequal weights break that. For (d), set the bound equal to 10⁻⁴ × I_xx and solve for n.
5. With y = s sin θ, y L′ dy becomes L₀ s² sin θ cos²θ dθ, and d(cos³θ)/dθ = −3 cos²θ sin θ. The root moment is W/2 times the centroid of the load distribution.
6. Each item attacks one of: constant σ (a point mass adds m y²), the planform near the root (a polynomial you can add), tabulated data (the tip entry hides a square root), the take-off model (the denominator T − k v² changes, and so does the partial-fraction split).
7. For the linear taper expand [1 − (1 − τ)u]² before integrating; for the ellipse c² = c_r²(1 − y²/s²) and the integrand is y² − y⁴/s². Simpson's error carries f⁽⁴⁾, which for a quartic is a non-zero constant: the error is exactly (b − a)⁵ f⁽⁴⁾/(180 n⁴).
