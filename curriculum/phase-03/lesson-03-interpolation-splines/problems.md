# P03.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are bounds and which are observed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the h²/8 rule on the troposphere)

The ISA density in the troposphere is ρ(h) = ρ₀ (1 − h/H)^4.2559 with ρ₀ = 1.225 kg/m³ and H = T₀/L = 44 331 m. A table lists ρ every 2 km.

(a) Differentiate twice and evaluate ρ″ at h = (3 + 0.1 s) km.
(b) Predict the maximum error of linear interpolation on the 2 km interval containing h from the h²/8 bound, and the error at h itself from the exact formula ½|ρ″| (h − hᵢ)(hᵢ₊₁ − h).
(c) Compute the linear, log-linear and natural-spline values from `isa_table(2000.0)` with `interpolate_isa` and compare their errors against `isa_atmosphere`.

Answer form: ρ″ with units, two predicted errors, three observed errors with a one-line comment each.

## 2. (Derivation: the error formula, the h²/8 bound and the spline equations)

(a) Prove that the polynomial of degree ≤ n through n + 1 points with distinct abscissae is unique.
(b) Derive the error formula f(x) − p(x) = f⁽ⁿ⁺¹⁾(ξ) Π(x − xᵢ)/(n + 1)!. Hint: for fixed x not a node, consider g(t) = f(t) − p(t) − [f(x) − p(x)] Π(t − xᵢ)/Π(x − xᵢ), which has n + 2 zeros, and apply Rolle's theorem repeatedly.
(c) From (b) with n = 1 derive the bound max|f − p| ≤ h² max|f″|/8 on an interval of width h.
(d) Starting from the piecewise form (4) of the lesson, compute S′ at the right end of [xᵢ₋₁, xᵢ] and at the left end of [xᵢ, xᵢ₊₁], equate them, and obtain the interior equation (5).

Answer form: four derivations.

## 3. (Unfamiliar scenario: a stall in a lift table)

A wind-tunnel-style table (SYNTHETIC: generated from a formula, not measured) gives the lift coefficient of a wing every 2° of angle of attack from 0° to 20°: C_L = 0.1 α for α ≤ α_s and C_L = 0.1 α_s − 0.05 (α − α_s)² beyond, with the stall angle α_s = 12 + (s mod 4) degrees (so 12°, 13°, 14° or 15°).

(a) Tabulate the eleven values. Interpolate C_L at α_s − 0.5°, α_s + 0.5° and α_s + 2.5° with piecewise-linear, natural-spline, not-a-knot-spline and the degree-10 Lagrange polynomial, and compare with the formula.
(b) Find the maximum of each interpolant on a fine grid and compare with the true C_Lmax = 0.1 α_s at α_s. Which interpolants overshoot or undershoot the peak, and why?
(c) A flight-control law reads C_L from this table at 400 Hz near stall. Which interpolant would you specify, and what one change to the *table* would improve every interpolant at once?

Answer form: a table of 3 × 4 values with errors, four maxima, a two-sentence recommendation.

## 4. (Numerical, uses your code: convergence inside a layer and across a kink)

(a) With `max_interpolation_error` on the ISA density over 20–32 km (one lapse rate, no kink) and spacings 2000, 1000, 500, 250 m, tabulate the maximum error for linear, log-linear, natural, not-a-knot and clamped interpolation and fit the observed orders.
(b) Repeat over 9–13 km with spacings 1000, 500, 250, 125 m (the 11 km kink is always a node). Report the orders.
(c) Explain the difference between (a) and (b) for the splines using the smoothness the spline error bound requires, and say why linear and log-linear are unaffected.

Answer form: two 5 × 4 tables with orders, one paragraph.

## 5. (Numerical, uses your code: Runge and its cure)

Let n = 8 + (s mod 10) nodes.

(a) Interpolate 1/(1 + 25x²) on [−1, 1] through n equispaced nodes and through n Chebyshev nodes; report the maximum error of each on a fine grid (for Chebyshev, inside the outermost nodes).
(b) Find the smallest Chebyshev node count that brings the error below 10⁻³, and the equispaced count at which the error first exceeds 1.
(c) Interpolate the same function with a not-a-knot spline on 9, 17 and 33 equispaced nodes and report the errors and the observed order. Which of the three approaches would you use for a table of this function, and why?

Answer form: two errors, two node counts, three spline errors with an order, one sentence.

## 6. (Judgment: specifying the simulator's atmosphere table)

A six-degree-of-freedom simulator needs ISA density to 0.1 % everywhere from 0 to 86 km, pressure to 0.1 %, and temperature to 0.1 K, from one table with one interpolation routine per variable, evaluated 10⁵ times per second.

(a) Specify the table: spacing, whether the seven layer boundaries are rows, total rows.
(b) Specify the interpolant for each variable (linear, log-linear, global spline, per-layer spline) with the error you expect from the lesson's measurements.
(c) Say which combination you would reject even though it has the highest order, and why.
(d) The requirement tightens to 0.01 % on density. Which part of your specification changes, and what does it cost in rows?

Answer form: a specification of about ten lines with a number on every line.

## 7. (Proof: the natural spline has minimum curvature)

Let S be the natural cubic spline through (xᵢ, yᵢ), i = 0..n, and let g be any C² function on [x₀, xₙ] with g(xᵢ) = yᵢ.

(a) Show that ∫(g″)² dx = ∫(S″)² dx + ∫(g″ − S″)² dx + 2∫S″(g″ − S″) dx.
(b) Integrate the cross term by parts twice on each interval, use S‴ = constant on each piece, g − S = 0 at every knot and S″ = 0 at the ends, and show the cross term vanishes.
(c) Conclude that ∫(S″)² ≤ ∫(g″)² with equality only if g = S, and explain in one sentence what this says about the natural spline's end behaviour and why it costs accuracy when f″ ≠ 0 at the ends.

Answer form: a proof in three steps and one sentence.
