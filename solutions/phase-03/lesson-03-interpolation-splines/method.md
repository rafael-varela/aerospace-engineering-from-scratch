# P03.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. ρ = ρ₀uᵏ, u = 1 − h/H, so ρ′ = −ρ₀k uᵏ⁻¹/H and ρ″ = ρ₀k(k − 1)uᵏ⁻²/H². 2. Bound: (2000)² |ρ″|/8. 3. Pointwise: ½|ρ″|(h − hᵢ)(hᵢ₊₁ − h) with hᵢ = 2 km ⌊h/2 km⌋. 4. `interpolate_isa(isa_table(2000.0), h, "rho", method)` for the three methods, minus `isa_atmosphere(h)[2]`.

## 2.
(a) Let p, q both interpolate; r = p − q has degree ≤ n and r(xᵢ) = 0 for n + 1 distinct xᵢ; a nonzero polynomial of degree ≤ n has at most n roots. (b) g(t) has n + 2 distinct zeros (the nodes and x); Rolle gives g′ n + 1 zeros, …, g⁽ⁿ⁺¹⁾ one zero ξ; compute g⁽ⁿ⁺¹⁾ using p⁽ⁿ⁺¹⁾ = 0 and (Π(t − xᵢ))⁽ⁿ⁺¹⁾ = (n + 1)!; solve for f(x) − p(x). (c) With n = 1, |f − p| = ½|f″(ξ)||(x − xᵢ)(x − xᵢ₊₁)|; the product is maximal at the midpoint with value h²/4. (d) Differentiate (4): S′ = (−Mᵢa² + Mᵢ₊₁b²)/(2hᵢ) + dᵢ − (Mᵢ₊₁ − Mᵢ)hᵢ/6; evaluate at x = xᵢ from the right (a = hᵢ, b = 0) and at x = xᵢ from the left on the previous interval (a = 0, b = hᵢ₋₁); set equal; multiply by 6.

## 3.
(a) α = 0, 2, …, 20; C_L by the two-branch formula; evaluate each interpolant at the three angles; subtract the formula. (b) Fine grid, `max` of each interpolant, location of the maximum; compare with 0.1α_s. (c) Choose on the basis of overshoot near the peak and cost; the table change is about where the samples are taken.

## 4.
(a) For each method and spacing: `max_interpolation_error(f_rho, 20000, 32000, int(12000/step), method)`; `observed_order(steps, errs)`. (b) Same with (9000, 13000) and the finer spacings. (c) The spline bound 5h⁴f⁗/384 assumes f ∈ C⁴; a jump in f′ gives the spline an O(h) error near the knot (its C² pieces must bend through the corner). The linear and log-linear interpolants are local to one interval and exact for a kink that sits on a node.

## 5.
(a) `lagrange_interpolate` on `numpy.linspace(-1, 1, n)` and on `chebyshev_nodes(n, -1, 1)`; maximum absolute error on a 2001-point grid (Chebyshev: restrict to [x₀, xₙ₋₁]). (b) Loop n upward (Chebyshev) until error < 10⁻³; loop n upward (equispaced) until error > 1. (c) `spline_evaluate(cubic_spline(xs, runge(xs), "not-a-knot"), xf)` for 9, 17, 33 nodes; `observed_order` with h = 2/(n − 1).

## 6.
(a) Spacing from the lesson's convergence tables; boundaries as rows; rows = 86/spacing + 1 + 7 (minus coincidences). (b) Density: log-linear or per-layer spline, with the lesson's measured maxima; pressure: log-linear (exact in isothermal layers); temperature: linear (exact). (c) The global spline across kinks: order 1 there, worst global maximum at every spacing. (d) Use the part 5 numbers for 0.01 %: which method reaches it at which spacing, and the row counts.

## 7.
(a) Algebra: (g″)² = (S″)² + (g″ − S″)² + 2S″(g″ − S″). (b) On [xᵢ, xᵢ₊₁]: ∫S″(g − S)″ = [S″(g − S)′]ᵢⁱ⁺¹ − ∫S‴(g − S)′ = [S″(g − S)′]ᵢⁱ⁺¹ − S‴ᵢ[(g − S)]ᵢⁱ⁺¹; the last bracket is zero (interpolation); sum over intervals; the boundary terms telescope (S″ and (g − S)′ are continuous) to S″(xₙ)(g − S)′(xₙ) − S″(x₀)(g − S)′(x₀) = 0 (natural ends). (c) Drop the non-negative middle term; equality iff g″ = S″ everywhere, which with equal values at two or more points forces g = S. The natural spline is the flattest-at-the-ends interpolant, which is wrong when the function is not flat there.
