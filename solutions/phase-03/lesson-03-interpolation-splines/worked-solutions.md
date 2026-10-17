# P03.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. The h²/8 rule on the troposphere (s = 15: h = 4.5 km)

**Derivation.** ρ = ρ₀uᵏ, u = 1 − h/H, k = 4.2559, H = 44 331 m. ρ″ = ρ₀k(k − 1)uᵏ⁻²/H². At h = 4500 m, u = 0.89849, uᵏ⁻² = u^2.2559 = 0.7846: ρ″ = 1.225 × 4.2559 × 3.2559 × 0.7846/(44 331)² = **6.78 × 10⁻⁹ kg/m³/m²**.

(a) as above. (b) Bound on [4000, 6000] m: h²|ρ″|/8 = 4 × 10⁶ × 6.78 × 10⁻⁹/8 = **3.39 × 10⁻³ kg/m³** (0.44 % of ρ = 0.7768). Pointwise: ½ × 6.78 × 10⁻⁹ × 500 × 1500 = **2.54 × 10⁻³ kg/m³**. (c) From the 2 km table: linear 0.779271, error **+2.50 × 10⁻³** (the chord lies above a convex curve; the pointwise prediction is right to 2 %); log-linear 0.775980, error **−7.9 × 10⁻⁴** (ln ρ is concave, so the chord lies below; (ln ρ)″ is a third of ρ″/ρ, and 2.5/3 ≈ 0.8 ✓); natural spline 0.776910, error **+1.4 × 10⁻⁴** (the spline's fourth-order term is small, but the 10–12 km kink interval two intervals away still leaks in).

**Check.** Relative linear error 3.2 × 10⁻³ against the lesson's 7.5 km value 3.9 × 10⁻³: larger at 7.5 km because ρ″/ρ grows with altitude (u shrinks) ✓.

**Common mistakes.** Using the interval width for the pointwise error (gives the bound, 3.4 × 10⁻³, not the 2.5 × 10⁻³ at that point); forgetting the 1/H² in ρ″ (an error of 2 × 10⁹).

## 2. Error formula, h²/8 bound and the spline equations

(a) If p and q both have degree ≤ n and interpolate, r = p − q has degree ≤ n and vanishes at n + 1 distinct points. A nonzero polynomial of degree ≤ n has at most n roots, so r ≡ 0. ∎

(b) Fix x not a node and set w(t) = Π(t − xᵢ), g(t) = f(t) − p(t) − [f(x) − p(x)] w(t)/w(x). Then g vanishes at the n + 1 nodes and at t = x: n + 2 distinct zeros. Rolle: g′ has n + 1 zeros, g″ has n, …, g⁽ⁿ⁺¹⁾ has one, call it ξ. But g⁽ⁿ⁺¹⁾(t) = f⁽ⁿ⁺¹⁾(t) − 0 − [f(x) − p(x)](n + 1)!/w(x), since p has degree ≤ n and w is monic of degree n + 1. Setting t = ξ: **f(x) − p(x) = f⁽ⁿ⁺¹⁾(ξ) w(x)/(n + 1)!**. ∎

(c) n = 1 on [xᵢ, xᵢ₊₁]: |f − p| = ½|f″(ξ)||(x − xᵢ)(x − xᵢ₊₁)|. The product of distances is a downward parabola with maximum (h/2)² at the midpoint, so |f − p| ≤ ½ max|f″| h²/4 = **h² max|f″|/8**. ∎

(d) From (4), S′(x) = (−Mᵢa² + Mᵢ₊₁b²)/(2hᵢ) + dᵢ − (Mᵢ₊₁ − Mᵢ)hᵢ/6. At the left end of [xᵢ, xᵢ₊₁] (a = hᵢ, b = 0): S′(xᵢ⁺) = −Mᵢhᵢ/2 + dᵢ − (Mᵢ₊₁ − Mᵢ)hᵢ/6 = dᵢ − hᵢ(2Mᵢ + Mᵢ₊₁)/6. At the right end of [xᵢ₋₁, xᵢ] (a = 0, b = hᵢ₋₁): S′(xᵢ⁻) = Mᵢhᵢ₋₁/2 + dᵢ₋₁ − (Mᵢ − Mᵢ₋₁)hᵢ₋₁/6 = dᵢ₋₁ + hᵢ₋₁(Mᵢ₋₁ + 2Mᵢ)/6. Equating and multiplying by 6: **hᵢ₋₁Mᵢ₋₁ + 2(hᵢ₋₁ + hᵢ)Mᵢ + hᵢMᵢ₊₁ = 6(dᵢ − dᵢ₋₁)**, which is (5). ∎

**Check.** For data on a straight line, dᵢ = dᵢ₋₁ and M = 0 solves (5): the spline is the line ✓.

**Common mistakes.** In (b), forgetting that x itself is a zero of g (then only n + 1 zeros, and the argument fails); in (d), evaluating S′ on the wrong piece at the knot.

## 3. A stall in a lift table (s = 15: α_s = 15°; SYNTHETIC data)

Table (α = 0, 2, …, 20°): C_L = 0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.45, 1.05, 0.25. The stall at 15° falls *between* the rows at 14° and 16°.

| α | true | linear | natural spline | not-a-knot | Lagrange (degree 10) |
|---|---|---|---|---|---|
| 14.5° | 1.4500 | 1.4125 (−0.038) | 1.4387 (−0.011) | 1.4392 (−0.011) | 1.4377 (−0.012) |
| 15.5° | 1.4875 | 1.4375 (−0.050) | 1.4705 (−0.017) | 1.4713 (−0.016) | 1.4684 (−0.019) |
| 17.5° | 1.1875 | 1.1500 (−0.038) | 1.1964 (+0.009) | 1.1934 (+0.006) | 1.2021 (+0.015) |

(b) True C_Lmax = 1.500 at 15°. Maxima: linear **1.45 at 16°** (it can never exceed a data point: undershoots by 0.05); natural spline **1.4712 at 15.38°** (undershoots by 0.029 and misplaces the peak by 0.4°); not-a-knot similar; Lagrange **1.4690**. None overshoots *here* because the peak is a corner between two samples that both lie below it; the smooth interpolants overshoot *beyond* the stall instead (17.5°: +0.006 to +0.015, the ringing of a smooth curve through a corner), and the maximum errors over the whole range are 0.075 (linear), 0.036 (spline), 0.038 (Lagrange).

(c) Near stall at 400 Hz: **piecewise linear** (or a monotone cubic), because it never produces a lift larger than the table and never rings; its 3 % undershoot at the peak is a known, one-sided bias that a control law can tolerate, while a spline's ringing changes sign with angle and feeds a derivative. The one table change: **sample at the stall angle itself** (add the row α = 15°, C_L = 1.5). Then every interpolant passes through the peak, the linear error at 14.5° and 15.5° drops to zero on the linear branch and to 0.0125 on the quadratic one, and the kink is on a node, where the lesson showed linear interpolation is exact and the spline's ringing is smallest.

**Check.** On the linear branch (α ≤ 14°) all four interpolants reproduce 0.1α exactly for the linear one and to 10⁻³ for the others ✓.

**Common mistakes.** Calling the Lagrange polynomial "most accurate because highest degree" (it is comparable to the spline here and would diverge with more rows); reading C_Lmax from the table (1.45) as the wing's maximum.

## 4. Convergence inside a layer and across a kink

(a) 20–32 km (one lapse rate, +1 K/km), max |error| in kg/m³:

| step m | linear | log-linear | natural | not-a-knot | clamped |
|---|---|---|---|---|---|
| 2000 | 1.01 × 10⁻³ | 2.80 × 10⁻⁵ | 4.61 × 10⁻⁴ | 2.25 × 10⁻⁵ | 2.79 × 10⁻⁶ |
| 1000 | 2.74 × 10⁻⁴ | 7.58 × 10⁻⁶ | 1.17 × 10⁻⁴ | 1.69 × 10⁻⁶ | 1.82 × 10⁻⁷ |
| 500 | 7.14 × 10⁻⁵ | 1.98 × 10⁻⁶ | 2.92 × 10⁻⁵ | 1.16 × 10⁻⁷ | 1.15 × 10⁻⁸ |
| 250 | 1.82 × 10⁻⁵ | 5.04 × 10⁻⁷ | 7.32 × 10⁻⁶ | 7.59 × 10⁻⁹ | 7.27 × 10⁻¹⁰ |
| order | **1.93** | **1.93** | **1.99** | **3.85** | **3.97** |

Log-linear is 36 times better than linear at every spacing (the ratio ρ″/(ρ(ln ρ)″) in this layer); the natural spline is second order (f″ ≠ 0 at 20 and 32 km); not-a-knot and clamped are fourth order, the clamped one eight times better because it is given the true end slopes.

(b) 9–13 km (kink at 11 km, always a node):

| step m | linear | log-linear | natural | not-a-knot |
|---|---|---|---|---|
| 1000 | 1.05 × 10⁻³ | 1.93 × 10⁻⁴ | 1.00 × 10⁻³ | 1.19 × 10⁻³ |
| 500 | 2.72 × 10⁻⁴ | 4.89 × 10⁻⁵ | 4.91 × 10⁻⁴ | 4.92 × 10⁻⁴ |
| 250 | 6.93 × 10⁻⁵ | 1.23 × 10⁻⁵ | 2.39 × 10⁻⁴ | 2.39 × 10⁻⁴ |
| 125 | 1.75 × 10⁻⁵ | 3.09 × 10⁻⁶ | 1.18 × 10⁻⁴ | 1.18 × 10⁻⁴ |
| order | **1.97** | **1.99** | **1.03** | **1.11** |

(c) The spline bound 5h⁴ max|f⁗|/384 requires f ∈ C⁴ across every interval, including the knot at 11 km. There ρ′ jumps (the lapse rate changes from −6.5 to 0 K/km), so f⁗ does not exist and the bound says nothing. What actually happens: the C² condition forces the two cubics meeting at 11 km to share a slope, while the data demand a corner; the spline splits the slope jump between the two sides, producing an O(h × Δslope) error that decays over a few intervals on each side: first order, and 14 times worse than linear interpolation at 125 m. The linear and log-linear interpolants use only the two bracketing nodes, never cross the kink when it is a node, and keep their second order from each layer's own curvature.

**Check.** At 250 m the two splines give the same error to three figures: the end condition no longer matters when the kink dominates ✓.

**Common mistakes.** Fitting the order over all four points when the first lies in a pre-asymptotic regime (fine here, but check the ratios step by step); concluding that "splines are bad" rather than "splines across kinks are bad".

## 5. Runge and its cure (s = 15: n = 13 nodes)

(a) Equispaced, 13 nodes: max error **3.66**; Chebyshev, 13 nodes: **0.069** (inside the outermost nodes ±0.9927).

(b) Chebyshev: the error first drops below 10⁻³ at **n = 35** nodes (9.5 × 10⁻⁴). Equispaced: the error first exceeds 1 at **n = 9** (1.045; at n = 7 it is 0.62, at n = 8, 0.25: even counts put a node at the centre and do better for a while, but n = 11 gives 1.9 and n = 21 gives 60).

(c) Not-a-knot spline on 9, 17, 33 equispaced nodes: **0.056, 3.7 × 10⁻³, 6.6 × 10⁻⁴**, observed order 3.2 over this range (4.05 from 33 to 129 nodes, where the asymptotic regime is reached). For a table of this function: **the spline on equispaced nodes**. It converges at fourth order with no node placement freedom needed (a table usually comes on a fixed grid), while the Chebyshev polynomial requires choosing where to sample and the equispaced polynomial is unusable beyond 8 nodes.

**Check.** The spline at 33 nodes (6.6 × 10⁻⁴) is comparable to Chebyshev at 35 nodes (9.5 × 10⁻⁴); the spline wins by a wide margin at 129 nodes (2.4 × 10⁻⁶) where the polynomial would need ~70 Chebyshev nodes.

**Common mistakes.** Measuring the Chebyshev error on the full [−1, 1] (that includes extrapolation beyond the outer nodes, which the code refuses); taking the equispaced error at n = 8 (0.25) as evidence that the phenomenon is mild.

## 6. Specifying the simulator's atmosphere table

Specification (numbers from the lesson's parts 4 and 5):

1. Rows every 1000 m from 0 to 86 km **plus the seven layer boundaries** (11, 20, 32, 47, 51, 71, 84.852 km; six coincide with the grid): **88 rows**, 4 columns.
2. Density: **log-linear in ρ**; measured maximum relative error 4.7 × 10⁻⁴ from 0 to 86 km (2 km rows with boundaries: 1.7 × 10⁻³, fails).
3. Pressure: **log-linear in p**; exact in the two isothermal layers, 2.1 × 10⁻³ at 2 km and below 10⁻³ at 1 km (the lesson's factor 7.5 over linear).
4. Temperature: **piecewise linear**; exact (10⁻¹⁴ K) because T is piecewise linear with every kink on a row.
5. Evaluation: one `searchsorted` and four flops per call; 10⁵ calls/s is 10⁻³ of a core.
6. Alternative for density: **one not-a-knot spline per layer** on the 2 km-plus-boundaries table (49 rows): 3.7 × 10⁻⁴, fewer rows, more code (eight spline objects and a layer lookup).
7. Rejected: **a single cubic spline across all rows**, although it is the only fourth-order method here, because across the kinks it is first order: 3.2 × 10⁻³ at 2 km with boundaries, 3.0 × 10⁻³ at 1 km, 6.7 × 10⁻⁴ at 250 m, 3.3 × 10⁻⁴ at 125 m: it never meets 0.1 % with fewer than 345 rows and never meets 0.01 %.
8. For 0.01 % on density: log-linear needs **250 m rows (345 rows, 3.0 × 10⁻⁵)**; the per-layer spline needs 500 m (174 rows, 1.7 × 10⁻⁶); linear would need 125 m (689 rows). The spacing changes, the method does not; the cost is 4× the rows for log-linear or 2× for the per-layer spline.

**Check.** Each quoted error is from the convergence runs of `main()` and part 4; the orders (2 for log-linear) predict that 1 km → 250 m cuts the error by 16: 4.7 × 10⁻⁴/16 = 2.9 × 10⁻⁵ against the measured 3.0 × 10⁻⁵ ✓.

**Common mistakes.** Choosing "the spline" for its order without checking the kinks; forgetting 84.852 km, which no round-number grid contains; quoting an absolute error in kg/m³ where the requirement is relative (the 86 km density is 7 × 10⁻⁶ kg/m³).

## 7. The natural spline has minimum curvature

(a) (g″)² = ((g″ − S″) + S″)² = (g″ − S″)² + 2S″(g″ − S″) + (S″)². Integrate over [x₀, xₙ]. ∎

(b) Let e = g − S, which is C², vanishes at every knot, and has e″ = g″ − S″. On [xᵢ, xᵢ₊₁], integrate by parts: ∫S″e″ = [S″e′]ᵢⁱ⁺¹ − ∫S‴e′. On this interval S‴ is a constant cᵢ, so ∫S‴e′ = cᵢ[e]ᵢⁱ⁺¹ = cᵢ(e(xᵢ₊₁) − e(xᵢ)) = 0. Summing the boundary terms over all intervals, the interior contributions cancel in pairs because S″ and e′ are continuous at every knot (S is C², g is C²), leaving S″(xₙ)e′(xₙ) − S″(x₀)e′(x₀), which is **zero because S″ = 0 at both ends** (natural condition). Hence ∫S″(g″ − S″) = 0. ∎

(c) With the cross term gone, ∫(g″)² = ∫(S″)² + ∫(g″ − S″)² ≥ ∫(S″)², with equality iff g″ = S″ everywhere; then g − S is linear, and since it vanishes at two or more knots it is zero: g = S. ∎ The natural spline is the interpolant with the least total curvature, which it achieves by going *straight* into the ends; when the true function curves there (f″(x₀) ≠ 0), the spline is forced flat where it should bend and the end intervals carry an O(h²) error, which is why its global order is 2 unless f″ = 0 at the ends (then 4, as measured on sin over [0, π]).

**Check.** The lesson's test: natural order 4.0 on sin (f″ = 0 at 0 and π), 2.0 on e^{−x/2}cos 2x ✓.

**Common mistakes.** Dropping the boundary term without invoking S″ = 0 at the ends (a clamped spline does *not* minimise curvature among all C² interpolants, only among those with the same end slopes); forgetting that S‴ is constant per piece, not globally.

## Check-yourself answers

1. Exactly one: the difference of two would be a polynomial of degree ≤ n with n + 1 roots.
2. f(x) − p(x) = f⁽ⁿ⁺¹⁾(ξ)/(n + 1)! × Π(x − xᵢ): the derivative factor (the function's smoothness) and the node polynomial (the node placement).
3. The derivatives of 1/(1 + 25x²) grow like 5ⁿ⁺¹(n + 1)!, cancelling the factorial, and the equispaced node polynomial is large near the ends; Chebyshev nodes xₖ = cos((2k + 1)π/2n) minimise the node polynomial's maximum and give 0.015.
4. Error ≤ h² max|f″|/8 = (2000)² × 6.8 × 10⁻⁹/8 = 3.4 × 10⁻³ kg/m³ (0.44 %).
5. S, S′ and S″ are continuous at every interior knot. Interpolating both ends of each of the n pieces (2n conditions) already makes S continuous; continuity of S′ and S″ at the n − 1 interior knots adds 2(n − 1), for 4n − 2 conditions on 4n coefficients, leaving two for the ends.
6. Natural (M₀ = Mₙ = 0), clamped (end slopes given), not-a-knot (S‴ continuous at x₁ and xₙ₋₁). Natural is second order in general and fourth order only when f″ = 0 at both ends.
7. Temperature is piecewise linear with its kink at 11 km on a table row, so the linear interpolant is exact; the C² spline cannot represent a corner and bends through it, ringing by 0.5 K over the neighbouring intervals.
8. Because p falls nearly exponentially, ln p is nearly linear and the linear interpolant of ln p has a far smaller second-derivative factor; it is exact in the isothermal layers 11–20 km and 47–51 km.
