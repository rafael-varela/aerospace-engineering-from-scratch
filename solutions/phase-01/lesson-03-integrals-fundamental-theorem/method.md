# P01.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. t = 0.08 + s/400, c = 1 + s/20. 2. A = 0.685083 t c². 3. x̄ = 0.42044 c. 4. `section_area(t, c, n=1000)`; subtract A. 5. Predicted |error| = 0.022 (t/0.12) c² × 1000^(−1.5). 6. `section_area_substituted(t, c, n=1000)`; subtract A.

## 2.
(a) A = 2∫₀ᶜ y_t dx; put ξ = x/c, dx = c dξ; the integral of each term c · 5tc · aₖ ∫₀¹ ξ^(eₖ) dξ by the power rule and FTC part 2 (G(1) − G(0) with G = ξ^(e+1)/(e+1)); sum. For x̄ the extra factor x = cξ raises each exponent by one. (b) Linear: the midpoint value times h equals the trapezoid area of a straight segment. Error: integrate the Taylor expansion over one panel; the (x − x_m)² term gives f″(x_m)h³/24; sum over panels; replace Σ f″(x_m)h by ∫f″ = f′(b) − f′(a) (needs f″ continuous, f′ finite at the ends). (c) Substitute, expand, multiply by 2u; g′ is a polynomial, evaluate at 0 and 1; error = (h²/24)[g′(1) − g′(0)] with h = 1/1280.

## 3.
(a) Trapezoid areas: ½(F_i + F_{i+1})(t_{i+1} − t_i) over the four segments; average = I/1.5. (b) `Fi = lambda x: np.interp(x, tp, Fp)`; `riemann_sum(Fi, 0.0, 1.5, n, rule)`. With n = 15 the edges are multiples of 0.1 s, so every breakpoint is an edge; midpoint is exact per linear piece; left and right sums equal the trapezoid sum because they differ only by h[F(1.5) − F(0)] = 0. (c) Δv₁ = I/m − g·1.5. Δv₂ = `riemann_sum(lambda x: Fi(x)/(0.25 - 0.02*x), 0, 1.5, 3000)` − g·1.5. Δv₃ = I/m_avg − g·1.5. (d) Drag at these speeds, and the thrust curve's own tolerance.

## 4.
(a) `convergence_study` for each rule on both windows; `observed_order`. (b) Compute `riemann_sum(...) - A` (signed) for left and right; tabulate; average the two sums and fit the order. (c) `cumulative_integral(f, 0, [x+h, x-h], n=2000)`, central difference, subtract `f(x)`; the leading-edge point has the largest f‴ and the sum's own error enters through the h = 10⁻⁴ subtraction. (d) `ax.loglog` for all series.

## 5.
(a) `section_centroid(0.12, 1.0, m, 0.4, n=4000)[1]` for the three m; divide by m. (b) Both numerator and denominator of ȳ are integrals; only the numerator is proportional to m. A weighted mean of y_c with positive weights lies between min and max of y_c. (c) θ_max = arctan(max dy_c/dx); for the four-digit camber line the slope is largest at the leading edge, 2m/p. Area change from the horizontal shift: upper surface moves back by y_t sin θ, lower forward by the same, so the net area change is second order in θ. (d) Build the two surface polylines, close the loop, apply ½|Σ(x_i y_{i+1} − x_{i+1} y_i)|.

## 6.
For each: name the property that fails (bounded, finite integral, continuity, resolution, noise); give the expected order (0.5 for x^(−1/2); none for 1/x; 1 for a step; 2 once h ≪ period, but garbage before; limited by noise); the observable (error ratio per doubling; sum growing like log n; error jumping with n; error oscillating with n; error floor independent of n).

## 7.
(a) `yt = 5*t*c*sum(...)`; `sp.integrate(sp.expand(yt**3), (x, 0, c))` times 2/3; read off κ. (b) `convergence_study(lambda x: (2/3)*naca4_thickness(x,0.12)**3, I, 0, 1, ns)`; `observed_order`. Compare (h²/24)[f′(c) − f′(0)] (f′(0) = 0 because y_t³ ∝ x^(3/2); f′(c) = 2y_t(c)² y_t′(c), tiny) with the x^(3/2) singular term ∝ h^(5/2). (c) √(I/A); rectangle I = c(tc)³/12; the section's material sits near the chord line where it contributes little to I. (d) I ∝ t³ c⁴: (15/12)³ and (0.30/0.25)⁴; mass ∝ A ∝ t c².
