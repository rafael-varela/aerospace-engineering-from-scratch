# P03.03  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Write ρ = ρ₀ uᵏ with u = 1 − h/H and k = 4.2559; then ρ″ = ρ₀ k(k − 1) uᵏ⁻²/H². The bound uses the interval width (2000 m); the pointwise formula uses the two distances from h to the interval ends.
2. (a) Two interpolants differ by a polynomial of degree ≤ n with n + 1 roots. (b) g has zeros at the n + 1 nodes and at x; each application of Rolle removes one zero and one degree; after n + 1 applications g⁽ⁿ⁺¹⁾(ξ) = 0 and p⁽ⁿ⁺¹⁾ = 0. (c) Maximise |(x − xᵢ)(x − xᵢ₊₁)| on the interval. (d) Differentiate (4) once; at x = xᵢ the factor a or b vanishes on each side.
3. Build C_L with `numpy.where`. The spline and the polynomial are smooth curves forced through a corner; the linear interpolant cannot exceed its data. A table that samples at the stall angle itself changes what every interpolant sees.
4. `max_interpolation_error(lambda h: isa_atmosphere(h)[2], a, b, n_intervals, method)` with n_intervals = (b − a)/step; `observed_order(steps, errs)`. The spline bound needs f ∈ C⁴; across 11 km, ρ′ jumps.
5. The lesson's `main()` already prints the equispaced and Chebyshev errors for n = 5, 9, 13, 17, 21; loop over n for the thresholds. For the spline use `cubic_spline(xs, runge(xs), "not-a-knot")` on equispaced nodes; the function is smooth, so the spline does not care about Runge's phenomenon.
6. The numbers are in the lesson's part 4 and part 5 tables: the plain 2 km table fails every method; adding the boundaries rescues log-linear and the per-layer spline but not the global spline; linear needs 500 m. Temperature is piecewise linear: pick accordingly.
7. Expand (g″)² = ((g″ − S″) + S″)². For the cross term, ∫S″(g − S)″ = [S″(g − S)′] − ∫S‴(g − S)′ on each interval, and S‴ is a constant there, so the second integral is S‴[(g − S)] between knots, which is zero. The boundary terms telescope and vanish at the ends because S″ = 0 there.
