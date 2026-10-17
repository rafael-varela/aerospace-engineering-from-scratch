# P03.03  Interpolation & Splines

> A table is a finite set of facts. Everything you read between its rows is a model, and the model's error is as much yours to know as the measurement's.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the uniqueness of the interpolating polynomial, its error formula f(x) − p(x) = f⁽ⁿ⁺¹⁾(ξ) Π(x − xᵢ)/(n + 1)!, the h²/8 bound for piecewise-linear interpolation and the second-derivative equations of the cubic spline. 2. Implement and test from scratch the Lagrange and Newton forms, Chebyshev nodes, piecewise-linear, log-linear and C² cubic-spline interpolation with natural, clamped and not-a-knot ends, against `scipy.interpolate` as oracle. 3. Measure the observed convergence order of each interpolant on a smooth function (2, 2, 4) and on a function with a kink (1 for the spline), and explain both. 4. Decide, for an atmosphere table in a flight simulator, the spacing, the node placement, the variable to interpolate and the method that meet a 0.1 % density requirement. |
| **Prerequisites** | `P00.03` (the ISA atmosphere, which this lesson treats as the truth), `P03.01` (observed order, rounding floor), `P01.05` (Taylor remainder), `P01.13` (solving a small linear system; the spline needs one). Skill check: you can write the straight line through two points and expand f(x + h) to second order. |
| **Engineering question** | A trajectory simulator carries the standard atmosphere as a table every 2 km and must deliver density to 0.1 % at any altitude below 86 km, thousands of times per second. Which interpolant (linear, log-linear, cubic spline), which table spacing and which node placement meet the requirement, how large is the error at 7.5 km today, and where in the atmosphere does each method fail worst and why? |
| **Assumptions and limits** | The ISA of `P00.03` is the truth: seven lapse-rate layers to 86 km geopotential, hydrostatic, ideal gas, dry air; the tabulated values are exact (no measurement noise); interpolants are evaluated only inside the table range (extrapolation raises). Polynomial error bounds need f ∈ Cⁿ⁺¹ on the interval; the spline bounds need C⁴ and, for the natural spline's fourth-order rate, f″ = 0 at the ends. The ISA's density is only C⁰ in its first derivative at layer boundaries (temperature has kinks), which is exactly where the smooth-function bounds stop applying. |
| **Analytical work** | Uniqueness of the interpolant through n + 1 nodes; the Lagrange basis and the Newton divided-difference form; the error formula and the node polynomial; why Chebyshev nodes minimise its maximum; the h²/8 linear bound; the spline's piecewise form in terms of knot second derivatives Mᵢ and the tridiagonal continuity equations; the three end conditions; the 5h⁴/384 clamped bound; log-linear exactness for exponentials; the condition number view of why equispaced high-degree interpolation diverges (Runge). |
| **Computational deliverable** | `starter/main.py`: `isa_table`, `lagrange_interpolate`, `newton_divided_differences`, `newton_evaluate`, `chebyshev_nodes`, `runge`, `piecewise_linear`, `log_linear_interpolate`, `cubic_spline`, `spline_evaluate`, `linear_error_bound`, `spline_error_bound`, `observed_order`, `max_interpolation_error`, `interpolate_isa` (`isa_atmosphere` is given). Tests prove: the ISA hits its defining values; Lagrange reproduces its nodes, equals `BarycentricInterpolator` to 10⁻¹², is exact for a quartic, and refuses extrapolation and repeated nodes; the Newton form equals the Lagrange form and its coefficients are incremental; Runge's error grows monotonically to 60 on 21 equispaced nodes and falls monotonically to 0.015 on Chebyshev nodes; piecewise linear equals `numpy.interp`; log-linear is exact for the isothermal-layer pressure to 10⁻¹²; all three spline end conditions equal `CubicSpline` to 10⁻¹⁰ in value, slope and curvature, are C² at every knot, and reproduce a cubic exactly (clamped, not-a-knot) or not (natural); the linear bound holds and is tight within a factor 2, the clamped bound holds; observed orders are 2.0 ± 0.15 (linear), 2.0 ± 0.3 (natural with f″ ≠ 0 at the ends), 4.0 ± 0.4 (not-a-knot), 4.0 ± 0.3 (clamped) and 4.0 ± 0.3 (natural on sin, where f″ = 0 at the ends); the ISA table errors fall by more than 2.5× when the spacing halves, log-linear beats linear by 4× for pressure, and the spline rings by 0.2–2 K around the tropopause kink while linear interpolation of the piecewise-linear temperature is exact. |
| **Practical activity** | **computational experiment**: interpolate density from the 2 km ISA table four ways and map the relative error against altitude; measure the convergence order of each method against table spacing inside a smooth layer and across the tropopause; find the smallest spacing that meets 0.1 % with and without nodes at the layer boundaries; reproduce the Runge phenomenon and its Chebyshev cure. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `atmosphere_interpolation.png`, the error-against-altitude table, the convergence table with observed orders, the Runge table, the hand derivations of the error formula and the spline equations, and a half-page note specifying the simulator's atmosphere table (spacing, nodes, variable, method) with the error it guarantees. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the error formula and the spline interior equation are derived; all public tests pass; the convergence orders are reported with the altitude range they were measured on and the kink explained; the note states a spacing and a method with a number. |
| **References** | `@ocw-18-330` "Lecture notes 5–8: Interpolation" (required-free): Lagrange and Newton forms, the error formula, Runge, splines. `@ocw-18-330` "Lecture notes 21–25: Spectral Interpolation, Differentiation, Quadrature" (supplementary-free): Chebyshev nodes and why they work. `@ocw-18-335j` Week 13 "Lecture 35: Chebyshev Approximation" (supplementary-free). `@scipy-userguide` "Interpolation (scipy.interpolate)" (required-free): `interp1d`, `CubicSpline`, `BarycentricInterpolator`, `PchipInterpolator`. |
| **Extensions** | ★ Implement the barycentric form of Lagrange interpolation (O(n) per point after O(n²) weights) and show it is stable where the direct form loses digits at n = 40. ★ Implement a monotone cubic (PCHIP) interpolant and show it removes the spline's overshoot at the tropopause and at an airfoil's stall. ★ Fit one not-a-knot spline per ISA layer on ln ρ from a table that includes the layer boundaries and find the coarsest spacing meeting 10⁻⁴ (the lesson's per-layer spline on ρ gives 3.7 × 10⁻⁴ at 2 km and 1.7 × 10⁻⁶ at 500 m). |
| **Time** | unmeasured |

## 1. The engineering problem

Every trajectory code, flight simulator and performance tool needs the atmosphere: density for drag and lift, pressure for altimetry and engine thrust, temperature for the speed of sound. The International Standard Atmosphere (`P00.03`) defines all three analytically, layer by layer, but the analytic form costs a logarithm, a power and a layer search per call, and a real atmosphere (a sounding, a forecast grid, a planetary entry profile) has no formula at all. So codes carry a **table**, sampled at some altitude spacing, and interpolate. The table every 2 km from 0 to 86 km has 44 rows; a six-degree-of-freedom simulation asks it for density perhaps a hundred thousand times a second. The question is what happens between the rows.

Interpolation is the model you use when you have values at points and need values in between. The choice of model decides the error. Draw a straight line between neighbouring rows and the error is set by the curvature of the true profile and the square of the spacing. Fit one polynomial through all 44 rows and the error is enormous, not small: high-degree polynomials through equispaced data oscillate wildly near the ends, the **Runge phenomenon**. Fit a **cubic spline** (a chain of cubics with continuous slope and curvature) and the error falls as the fourth power of the spacing, but the spline *rings* near any kink in the data, and the ISA temperature has a kink at every layer boundary. Interpolate ln p instead of p and a straight line becomes exact where the atmosphere is isothermal. None of these facts is visible from the table; all of them follow from the error formulas of part 2.

Success criteria, written first: every interpolant written here equals its `scipy.interpolate` counterpart to 10⁻¹⁰; its observed convergence order on a smooth stretch of atmosphere matches the theory (2 for linear and log-linear, 4 for a clamped or not-a-knot spline); the density at 7.5 km from the 2 km table is known with an error estimate before it is computed; and the lesson ends with a table specification (spacing, node placement, variable, method) that meets 0.1 % with a stated margin.

## 2. Understand and derive

### 2.1 One polynomial through n + 1 points

Given nodes x₀ < x₁ < … < xₙ with values yᵢ, there is exactly one polynomial of degree ≤ n with p(xᵢ) = yᵢ. *Uniqueness*: the difference of two such polynomials has degree ≤ n and n + 1 roots, so it is zero. *Existence*: the **Lagrange form** writes it down,

**p(x) = Σᵢ yᵢ ℓᵢ(x), ℓᵢ(x) = Π_{j≠i} (x − xⱼ)/(xᵢ − xⱼ)**,   (1)

where ℓᵢ is 1 at xᵢ and 0 at every other node. It costs O(n²) per evaluation point and must be rebuilt if a node is added. The **Newton form** p(x) = c₀ + c₁(x − x₀) + c₂(x − x₀)(x − x₁) + … uses the **divided differences** cₖ = f[x₀, …, xₖ], built recursively (f[xᵢ, …, xⱼ] = (f[xᵢ₊₁, …, xⱼ] − f[xᵢ, …, xⱼ₋₁])/(xⱼ − xᵢ)) in O(n²) once and then evaluated by nested multiplication in O(n); adding a node adds one coefficient and leaves the others unchanged. The two forms are the same polynomial (the tests check this to 10⁻¹²), and the Newton form is also the one that reveals the error.

### 2.2 The error formula and what it says about nodes

If f ∈ Cⁿ⁺¹ on an interval containing the nodes and x, then for some ξ in that interval

**f(x) − p(x) = f⁽ⁿ⁺¹⁾(ξ)/(n + 1)! · Π_{i=0}^{n} (x − xᵢ)**.   (2)

(Problem 2b derives it by applying Rolle's theorem n + 1 times to an auxiliary function.) Two factors, two lessons. The derivative factor f⁽ⁿ⁺¹⁾/(n + 1)! is the function's business: for sin or exp it falls like 1/(n + 1)! and the interpolant converges fast for any node choice; for f(x) = 1/(1 + 25x²) the derivatives grow like 5ⁿ⁺¹ (n + 1)!, which cancels the factorial, and convergence depends entirely on the second factor. The **node polynomial** Π(x − xᵢ) is the node choice's business: on equispaced nodes it is tiny in the middle and large near the ends (the product of distances to the far nodes), so the error concentrates at the ends and, for Runge's function, grows without bound: `main()` measures 0.44, 1.05, 3.7, 14, 60 for n = 5, 9, 13, 17, 21 nodes. **Chebyshev nodes** xₖ = cos((2k + 1)π/2n) (mapped to the interval) cluster near the ends and make max|Π(x − xᵢ)| as small as any choice of n nodes can, 2¹⁻ⁿ on [−1, 1]; the same function on Chebyshev nodes gives 0.40, 0.17, 0.069, 0.033, 0.015, falling steadily. The rule for engineering tables: **never put one high-degree polynomial through many equispaced rows**. Either choose the nodes (Chebyshev, when you control where the data are taken) or go piecewise.

### 2.3 Piecewise linear: the h²/8 rule

Between neighbouring nodes use the straight line. On [xᵢ, xᵢ₊₁] of width h, (2) with n = 1 gives |f(x) − p(x)| = |f″(ξ)|/2 · |(x − xᵢ)(x − xᵢ₊₁)| ≤ |f″| h²/8, the product being largest at the midpoint:

**max error ≤ (h²/8) max|f″|**.   (3)

Second order: halve the spacing, quarter the error. For ISA density in the troposphere, ρ = ρ₀(T/T₀)^4.2559 has ρ″ ≈ 6.8 × 10⁻⁹ kg/m³/m² at 4.5 km, so a 2 km table predicts an error up to 3.4 × 10⁻³ kg/m³ (0.44 %), and the code measures 2.5 × 10⁻³ at 4.5 km (Problem 1): the bound is tight within a factor 1.4. The method never overshoots the data, is exact for piecewise-linear functions with knots on the nodes (the ISA temperature from a 1 km table: error 3 × 10⁻¹⁴ K), and is what `numpy.interp` does.

**Log-linear interpolation** applies (3) to ln f: interpolate ln y linearly and exponentiate. It is *exact* when y = A e^{Bx}, which is the ISA pressure in an isothermal layer (11–20 km, 47–51 km), and nearly so everywhere else in the atmosphere, where p and ρ fall roughly exponentially: the relative error becomes (h²/8)|(ln f)″|, and (ln ρ)″ = −4.26/(H − h)² in the troposphere is a third of ρ″/ρ. From the 2 km table at 7.5 km the density error falls from 0.39 % (linear) to 0.12 % (log-linear). For pressure over the whole table the maximum relative error falls from 1.6 % to 0.21 %, a factor 7.5. Physics tells you what to interpolate.

### 2.4 Cubic splines: smooth and fourth order

A **cubic spline** S on knots x₀ < … < xₙ is a cubic on each [xᵢ, xᵢ₊₁] with S, S′ and S″ continuous at every interior knot (C²). Write Mᵢ = S″(xᵢ). Because S″ is linear on each piece, S″(x) = (Mᵢ(xᵢ₊₁ − x) + Mᵢ₊₁(x − xᵢ))/hᵢ, and integrating twice with S(xᵢ) = yᵢ, S(xᵢ₊₁) = yᵢ₊₁ gives, with a = xᵢ₊₁ − x, b = x − xᵢ,

**S(x) = Mᵢ a³/(6hᵢ) + Mᵢ₊₁ b³/(6hᵢ) + (yᵢ − Mᵢhᵢ²/6) a/hᵢ + (yᵢ₊₁ − Mᵢ₊₁hᵢ²/6) b/hᵢ**.   (4)

S and S″ are continuous by construction. Demanding that S′ match from the left and the right at each interior knot gives (Problem 2d) one equation per interior knot:

**hᵢ₋₁ Mᵢ₋₁ + 2(hᵢ₋₁ + hᵢ) Mᵢ + hᵢ Mᵢ₊₁ = 6(dᵢ − dᵢ₋₁)**, dᵢ = (yᵢ₊₁ − yᵢ)/hᵢ,   (5)

n − 1 equations for n + 1 unknowns. Two **end conditions** close the system:

- **natural**: M₀ = Mₙ = 0 (zero curvature at the ends; minimises ∫S″² over all C² interpolants, Problem 7);
- **clamped**: prescribed end slopes S′(x₀) = s₀, S′(xₙ) = sₙ, giving 2h₀M₀ + h₀M₁ = 6(d₀ − s₀) and its mirror;
- **not-a-knot**: S‴ continuous at x₁ and xₙ₋₁, so the first two and last two pieces are each one cubic: (M₁ − M₀)/h₀ = (M₂ − M₁)/h₁ and its mirror. This is `scipy.interpolate.CubicSpline`'s default and needs no information beyond the data.

The system is tridiagonal (plus two rows), diagonally dominant, and solved here with `numpy.linalg.solve`; `P03.08` builds the O(n) tridiagonal solver. The error: for f ∈ C⁴, the clamped spline satisfies **max|f − S| ≤ (5/384) h⁴ max|f⁗|**, fourth order, and not-a-knot matches it in the interior. The natural spline forces S″ = 0 at the ends where f″ generally is not, which costs an O(h²) error in the end intervals; its order is 4 only when f″ vanishes at the ends (sin on [0, π]: measured 4.0) and 2 otherwise (measured 2.0 on e^{−x/2} cos 2x). On the smooth troposphere density from 2 km to 250 m spacing `main()` measures orders 1.91 (linear), 1.91 (log-linear), 1.92 (natural), 3.84 (not-a-knot), 3.85 (clamped), the last two from errors of 3.9 × 10⁻⁶ and 3.7 × 10⁻⁷ kg/m³ at 2 km spacing: *three orders of magnitude* better than linear at the same spacing. The not-a-knot order in the max norm comes out at 3.7 on some test functions (3.68 on e^{−x/2} cos 2x from 20 to 160 intervals) because its end pieces converge a little slower than the interior, which is at 4.01.

### 2.5 Kinks: where the smooth theory stops

Every bound above assumes derivatives that exist. The ISA temperature is piecewise linear with a kink at 11 km (lapse −6.5 K/km to 0), so its density has a jump in ρ′ there. Two consequences, both measured. First, the **spline rings**: from a 1 km table (which has a node at 11 km) the natural spline's temperature error peaks at 0.55 K near the tropopause and decays away from it, while linear interpolation is exact to 10⁻¹⁴. The C² requirement forces the cubic to bend smoothly through a corner it cannot represent, and the bend spreads over several intervals on each side. Second, the **order collapses**: on 9–13 km the spline's density error falls only as h¹ (observed 1.03 for natural, 1.11 for not-a-knot), linear interpolation keeps its h² (1.97) and log-linear too (1.99). Third, if the kink is *not* a node (a 2 km table has rows at 10 and 12 km, none at 11), every method's error is set by the kink interval and the 0–86 km maxima become 1.4 % (linear), 1.6 % (log-linear), 1.0 % (natural), 1.0 % (not-a-knot): all failing a 0.1 % requirement. The first fix costs nothing: **put the seven layer boundaries in the table** (49 rows instead of 44). With them the log-linear error drops to 0.17 %, the global spline only to 0.32 % (it still rings across the kinks, now from nodes) and linear to 1.4 % (its worst interval is now a smooth but strongly curved one near 33 km). The second fix is to **stop the spline at every kink**: fit one not-a-knot spline per layer, and the 49-row table delivers 3.7 × 10⁻⁴ everywhere, 500 m rows 1.7 × 10⁻⁶. This is the central design lesson of the phase for tables: node placement and respecting the data's smoothness matter more than the order of the interpolant, and the smoothest method is the worst one across a corner.

### 2.6 Conditioning and cost

Interpolation is a linear map from data to values, and its condition number, the **Lebesgue constant** Λ = max Σ|ℓᵢ(x)|, bounds how much a data error is amplified: 1 for piecewise linear (the interpolant stays between its neighbours), about 2 for a cubic spline, 2ⁿ/(e n ln n) for equispaced polynomial interpolation (10⁴ at n = 20: a 0.01 % error in one table entry becomes 100 % somewhere) and (2/π) ln n for Chebyshev. Runge's phenomenon is this conditioning made visible with exact data. Cost per evaluation: a binary search plus a few flops for the piecewise methods, O(n²) for the direct Lagrange form; a 44-row spline evaluated 10⁵ times a second is cheap, a 44-node polynomial is both expensive and wrong.

## 3. Build the computational model

Open `starter/main.py`. `isa_atmosphere` (and its layer-base helper) is given; everything else is yours.

### Step 1: polynomials

`lagrange_interpolate(xs, ys, x)` validates the nodes (1-D, strictly increasing, ≥ 2), refuses points outside [x₀, xₙ], and sums (1) with two nested loops over nodes, vectorised over x. `newton_divided_differences` builds the table in place:

```python
table = ys.copy(); coeffs[0] = table[0]
for k in range(1, n):
    table[k:] = (table[k:] - table[k - 1:-1]) / (xs[k:] - xs[:-k]); coeffs[k] = table[k]
```

`newton_evaluate` runs Horner from the last coefficient: `out = out * (x - xs[k]) + coeffs[k]`. `chebyshev_nodes(n, a, b)` is the mapped cosine formula, sorted. `runge` is 1/(1 + 25x²).

### Step 2: piecewise methods

`piecewise_linear` locates the interval with `np.searchsorted(xs, x, side="right") - 1`, clipped to [0, n − 1], and blends. `log_linear_interpolate` rejects non-positive values, interpolates ln y and exponentiates. `cubic_spline(xs, ys, bc, end_slopes)` assembles the (n + 1) × (n + 1) system of (5) with the chosen boundary rows and returns `Spline(xs, ys, M)`; `spline_evaluate(spline, x, deriv)` applies (4) or its first or second derivative:

```python
S  = (M[i]*a**3 + M[i+1]*b**3)/(6h) + (y[i] - M[i]*h*h/6)*a/h + (y[i+1] - M[i+1]*h*h/6)*b/h
S' = (-M[i]*a**2 + M[i+1]*b**2)/(2h) + (y[i+1] - y[i])/h - (M[i+1] - M[i])*h/6
S''= (M[i]*a + M[i+1]*b)/h
```

### Step 3: error models and studies

`linear_error_bound(h, f2)` = h²f₂/8; `spline_error_bound(h, f4)` = 5h⁴f₄/384; `observed_order` is the log-log slope. `max_interpolation_error(f, a, b, n_intervals, method)` samples f on equispaced nodes, interpolates onto a 2001-point grid and returns the maximum absolute error (for 'clamped' it estimates the end slopes with one-sided second-order differences). `interpolate_isa(table, h, column, method)` dispatches to the right interpolant for a column of `isa_table`.

Running `main()` prints

```
ISA check: 11 km T = 216.65 K (216.65), p = 22632 Pa (22632), rho = 0.3639 kg/m^3 (0.3639)

2 km table (44 rows): density at 7.5 km, truth 0.556623 kg/m^3
method         rho(7.5 km)  rel. error    max rel. error 0-86 km    max rel. error 0-11 km
linear            0.558800    3.91e-03                  1.43e-02                  5.91e-03
loglinear         0.555980    1.16e-03                  1.58e-02                  1.58e-02
natural           0.556431    3.45e-04                  1.00e-02                  1.00e-02
not-a-knot        0.556447    3.16e-04                  9.99e-03                  9.99e-03
natural spline against scipy.interpolate.CubicSpline: max |diff| = 2.22e-16
```

then the convergence table of 2.4, the tropopause ringing of 2.5, the Runge table of 2.2, and `atmosphere_interpolation.png` (error against altitude with the layer boundaries marked; convergence against spacing; Runge's function with 13 equispaced and 13 Chebyshev nodes).

### Verify

Run `python -m pytest curriculum/phase-03/lesson-03-interpolation-splines`. The tests:

1. **ISA**: sea level, 11, 20 and 47 km match the defining values; the 2 km table has 44 rows; negative altitude raises.
2. **Lagrange**: reproduces nodes to 10⁻¹³, equals `BarycentricInterpolator` to 10⁻¹², is exact for a quartic through five nodes, raises outside the range and for a repeated node.
3. **Newton**: equals Lagrange; the first five coefficients are unchanged when a sixth node is added; c₀ = y₀, c₁ = first divided difference.
4. **Runge**: equispaced errors increase monotonically past 10 at 21 nodes; Chebyshev errors decrease monotonically below 0.02; the node formula, count, ordering and symmetry are right.
5. **Piecewise linear and log-linear**: equal to `numpy.interp`; exact for the isothermal-layer pressure (where linear is off by more than 10 Pa); raise on negative values and outside the range.
6. **Splines** (each end condition): value, slope and curvature equal `CubicSpline` to 10⁻¹⁰/10⁻⁹/10⁻⁸ on random knots; C² at every interior knot; natural has M₀ = Mₙ = 0; clamped reproduces its end slopes; not-a-knot has continuous S‴ at x₁; clamped and not-a-knot reproduce a cubic exactly while natural does not; unknown end conditions and extrapolation raise.
7. **Bounds**: on sin over [0, π] the linear error is within (½, 1] of h²/8 and the clamped error under 5h⁴/384 for three spacings.
8. **Convergence orders**: linear 2.0 ± 0.15; natural 2.0 ± 0.3 on a function with f″ ≠ 0 at the ends and 4.0 ± 0.3 on sin; not-a-knot 4.0 ± 0.4; clamped 4.0 ± 0.3; `observed_order` exact on a synthetic law.
9. **The table**: halving the spacing cuts every method's pressure and density error by more than 2.5×; log-linear beats linear by 4× for pressure over the whole range and by 2× for density in the troposphere interior; log-linear is exact for pressure in 11–20 km; linear is exact for temperature from a 1 km table; the spline overshoots the tropopause by 0.2–2 K with the error at 5.5 km below 5 % of that at 10.5 km.

### Compare

`scipy.interpolate.CubicSpline` builds the same spline (our natural spline differs from it by 2 × 10⁻¹⁶) with a banded solver and offers `periodic` and derivative-specified ends; `interp1d(kind="linear")` and `numpy.interp` are the linear method; `BarycentricInterpolator` is the Lagrange polynomial in the numerically stable barycentric form (extension 1); `PchipInterpolator` is the monotone cubic that does not overshoot at kinks (extension 2); `Akima1DInterpolator` is another. Use them: they handle the degenerate cases (two knots, repeated x, huge n) and are compiled. What building ours bought is the ability to read an error: to look at a 0.5 K overshoot near the tropopause and know it is the C² constraint meeting a kink, not a data error; to predict from (3) that a 2 km table cannot give 0.1 % by linear interpolation before writing a line; and to know that `CubicSpline`'s default is not-a-knot and `natural` is one keyword away, with the order-2 end penalty that comes with it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Error against altitude.** From the 2 km table, map the relative density error of linear, log-linear and natural-spline interpolation on a 10 m grid from 0 to 86 km (the left panel of the figure). Predict from (3) the troposphere error of linear interpolation at the midpoint of the 6–8 km interval, and predict that every method peaks in the intervals that straddle 11, 20, 32, 47, 51 and 71 km.
2. **Convergence inside a layer.** On 20–32 km (one lapse rate, no kink), tabulate the maximum density error against spacing 2000, 1000, 500, 250 m for linear, log-linear, natural, not-a-knot and clamped. Predict orders 2, 2, 2, 4, 4 and that log-linear beats linear by the ratio |(ln ρ)″| ρ/|ρ″|.
3. **Convergence across a kink.** Repeat on 9–13 km with spacings 1000, 500, 250, 125 m (the kink at 11 km is always a node). Predict that the splines drop to order 1 and linear stays at 2. Then shift the range to 9.5–13.5 km so the kink falls between nodes at the 1000 m spacing and record what happens to the first point of each curve.
4. **Node placement.** Build a 2 km table with the seven layer boundaries inserted as extra rows (49 rows) and repeat study 1 for linear, log-linear and the global not-a-knot spline; then fit one not-a-knot spline per layer (log-linear in the two-row top segment) on the same table. Predict that log-linear improves most (its only large errors were in kink intervals), that the global spline improves least (it still rings from the kink nodes, 2.5), and that the per-layer spline reaches the smooth-layer level of a few 10⁻⁴ at 2 km.

Raw data format: CSV tables `h_m, rel_err_linear, rel_err_loglinear, rel_err_natural` (1 and 4), `step_m, method, max_abs_err` (2 and 3), each with a header naming the table spacing and the altitude range; plus `atmosphere_interpolation.png`.

What you should find: in study 1 the troposphere error of linear interpolation at 7 km is 0.50 % against the (3) bound of 0.58 % and the global maxima sit exactly in the kink intervals, 1.4 %, 1.6 % and 1.0 %; in study 2 the orders come out 1.93, 1.93, 1.99, 3.85 and 3.97 and the log-linear error is 36 times below the linear one at every spacing; in study 3 the splines show 1.03 and 1.11 against 1.97 for linear, and in the shifted range the coarsest point is 2 to 15 times worse than its neighbours' trend for every method; in study 4 the maxima become 1.4 % (linear, now set by the curved 32–47 km layer), 0.17 % (log-linear), 0.32 % (global spline) and 0.037 % (per-layer spline), and halving the spacing to 1 km with boundaries (88 rows) gives 0.35 %, 0.047 %, 0.30 % and 0.023 %.

## 5. Validate and improve

| quantity | theory (§2) | observed (code) | ratio | cause |
|---|---|---|---|---|
| ρ at 7.5 km, linear, 2 km table | ≤ h²ρ″/8 ≈ 0.40 % | 0.39 % | 0.98 | (3), tight at a midpoint |
| ρ at 7.5 km, log-linear | ≤ h²|(ln ρ)″|/8 ≈ 0.12 % | 0.12 % | 1.0 | exact for exponentials; ρ nearly is |
| ρ at 7.5 km, not-a-knot | O(h⁴): ~10⁻⁵ in-layer | 0.03 % | — | the 10–12 km kink interval contaminates its neighbours |
| ρ, 0–86 km, 49-row table: log-linear / global spline / per-layer spline | — | 0.17 % / 0.32 % / 0.037 % | — | kinks at the nodes: order 1 for a spline that crosses them, order 4 for one that stops |
| our natural spline vs `CubicSpline` | 0 | 2.2 × 10⁻¹⁶ | — | rounding |
| orders, 0–11 km | 2, 2, 2, 4, 4 | 1.91, 1.91, 1.92, 3.84, 3.85 | — | pre-asymptotic at 2 km; errors at the coarsest step |
| natural spline order, f″ = 0 at ends / not | 4 / 2 | 4.0 / 2.0 | — | end-condition penalty |
| not-a-knot order, max norm | 4 | 3.68–3.92 (function-dependent) | — | end pieces; interior measures 4.01 |
| tropopause ringing, 1 km table | O(h × slope jump) ≈ 0.5 K | 0.55 K | 1.1 | C² across a kink |
| linear T from 1 km table | exact | 3 × 10⁻¹⁴ K | — | piecewise-linear data on its knots |
| Runge, 21 equispaced / Chebyshev | grows / shrinks | 60 / 0.015 | — | node polynomial; Lebesgue constant |
| pressure, 2 km table, linear → log-linear | factor ~ ρ″/(ρ (ln ρ)″) ≈ 3–8 | 1.6 % → 0.21 %, factor 7.5 | — | p is closer to exponential than ρ |

The bounds hold and are tight where the theory applies (smooth layers) and silent where it does not (kinks), and the measured orders say which regime you are in without knowing the function: a spline reporting order 1 is telling you there is a kink in the data. The 0.1 % density requirement of the engineering question is met, from the numbers above, by a 1 km table with the layer boundaries (88 rows) and log-linear interpolation (0.047 %), or by a 2 km table with the boundaries (49 rows) and one not-a-knot spline per layer (0.037 %); it is not met by any method on the plain 2 km table, nor by a single spline across the kinks at any spacing down to 125 m (0.033 %, still above the per-layer spline at 2 km). That is the answer a specification should carry: not "use a spline" but "rows at 0, 1, …, 86 km plus the seven layer boundaries; log-linear in ρ; relative error below 5 × 10⁻⁴ by the convergence study of part 4; or the 49-row table with a per-layer spline for 4 × 10⁻⁴".

Improvements worth making: the barycentric Lagrange form (stable to n = 100s), a monotone interpolant for kinked data, and the per-layer spline on ln ρ (exact for the isothermal layers, order 4 elsewhere: 3.4 × 10⁻⁸ at 500 m rows). Improvement worth not making: a finer table with the kinks still between rows, or a finer table under a single global spline. Halving the plain table to 1 km moves six of the seven kinks *onto* nodes (which is why 1 km looks so much better than 2 km) but leaves 84.852 km between rows; the gain is accidental, not designed, and the global spline's 0.30 % at 1 km shows that a kink on a node is still a kink.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of the error formula (2), the h²/8 bound (3) and the spline equations (4)–(5); `atmosphere_interpolation.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) specifying the simulator's atmosphere table: spacing, node placement, variable interpolated, method, guaranteed relative error, and the evidence from part 4 for each.

### Adversarial self-checks

- **What if the spacing doubles to 4 km?** Linear error ×4, spline error ×16 in smooth layers; the kink intervals widen and every method's global maximum roughly doubles. Which requirement would you drop first, 0.1 % or the 4 km table?
- **What if the table came from a radiosonde instead of the ISA?** Then the data carry noise (0.5 K, 1 hPa), the spline interpolates the noise with Lebesgue constant ~2 and ringing between every pair of points, and the right tool is a smoothing spline or a fit (`P02.18`), not an interpolant. Interpolation assumes the data are exact.
- **Which assumption fails first in the real world?** Smoothness: real soundings have inversions and the tropopause is a sharp fold, so every spline rings somewhere; a monotone interpolant or linear-in-ln p is the robust default for atmospheres.
- **What measurement could falsify the model?** None directly, since the ISA is a definition; but comparing the interpolated table against the analytic ISA on a fine grid *is* the test, and a spline reporting a 0.5 K error where the analytic form is piecewise linear falsifies the assumption that the data are C².
- **Does the result survive a change of interpolant?** Linear, log-linear and spline agree to 0.03 % at 7.5 km, which is the honest uncertainty of the 2 km table there; where they disagree by 1 % (the 10–12 km interval) the table, not the method, is at fault.
- **Which decision changes if the requirement tightens to 0.01 %?** Log-linear needs 250 m rows (345 rows, 3.0 × 10⁻⁵); linear needs 125 m (689 rows, 5.6 × 10⁻⁵); a single global spline never gets there (3.3 × 10⁻⁴ at 125 m); the per-layer spline passes at 500 m (174 rows, 1.7 × 10⁻⁶). The decision to respect the layer boundaries does not change; the decision between log-linear and per-layer spline now turns on memory against code complexity.

## Check yourself

1. How many polynomials of degree ≤ n pass through n + 1 points with distinct abscissae, and why?
2. Write the interpolation error formula and name the two factors that control it.
3. Why does a degree-20 polynomial through 21 equispaced samples of 1/(1 + 25x²) have an error of 60, and what node choice brings it to 0.015?
4. State the maximum error of piecewise-linear interpolation with spacing h and evaluate it for ISA density at 4.5 km from a 2 km table (ρ″ ≈ 6.8 × 10⁻⁹ kg/m³/m²).
5. What does "C²" mean for a cubic spline, how many conditions does it impose, and how many are left for the ends?
6. Name the three end conditions and say which one has only second-order accuracy in general and when it recovers fourth order.
7. A spline through ISA temperatures from a 1 km table shows a 0.5 K error near 11 km while linear interpolation is exact there. Explain both facts.
8. Why interpolate ln p rather than p, and in which ISA layers is the log-linear interpolant exact?
