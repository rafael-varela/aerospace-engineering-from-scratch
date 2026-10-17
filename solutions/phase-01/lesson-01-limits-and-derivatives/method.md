# P01.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. C_Lmax = 1.2 + s/50; α = (s/2)° converted to radians. 2. u = a₀α/C_Lmax with a₀ = 2π. 3. C_L = C_Lmax tanh u. 4. dC_L/dα = a₀(1 − tanh²u) in rad⁻¹; multiply by π/180 for deg⁻¹. 5. Fraction = 1 − tanh²u. 6. `central_difference(lambda x: lift_coefficient(x, cl_max=clmax), alpha, 1e-4)`; the difference should be of order (h²/6)|C_L‴| ≈ 10⁻⁸.

## 2.
(a) Forward: subtract f(x) from the + expansion, divide by h; the leftover starts at (h/2)f″. Central: subtract the − expansion from the + expansion; even-order terms cancel, divide by 2h; leftover starts at (h²/6)f‴. (b) D(2h) = f′ + (4h²/6)f‴ + (16h⁴/120)f⁽⁵⁾ + …; form [4D(h) − D(2h)]/3 and collect the function values: weights −1, 8, −8, 1 over 12h; the leftover is −(h⁴/30)f⁽⁵⁾. (c) Replace h by h/2 in (b)'s combination and compare with the five-point stencil at step h/2; the same four points with the same weights. (d) A formula is exact when every leftover term contains a derivative that vanishes: degree ≤ 1 for forward, ≤ 2 for central, ≤ 4 for five-point.

## 3.
(a) Chain rule through T: dM/dt = (dM/dT)(dT/dh)ḣ with dM/dT = −V/(2a T) = −M/(2T). (b) T from the troposphere formula at each altitude, a = √(γRT), M = V/a, dT/dh = −0.0065 K/m; substitute. The sign follows from dT/dh < 0: a falls as the air cools, so M rises at constant V. (c) Define `M = lambda h: V / np.sqrt(gamma * R * isa_temperature_two_layer(h))`, then `central_difference(M, h, 1.0) * hdot`. (d) `one_sided_derivatives(M, 11000.0, 10.0)`: the left value is dM/dh in the troposphere, the right is zero because T is constant above the tropopause.

## 4.
(a) `error_vs_step` for both functions; `observed_order` on h = 10⁻² … 10⁻⁴ for forward and central and on 10⁻² … 10⁻³ for five-point (it reaches rounding sooner). (b) `sp.diff(expr, x, 3)` at the point; `optimal_step_central(|f|, |f‴|)`; compare with `hs[np.argmin(err)]`. (c) Write h_opt ∝ (|f|/|f‴|)^(1/3) for both and take the ratio. (d) `np.imag(lift_coefficient(x + 1j * h)) / h`; the only operations are a complex tanh and a division, no cancellation.

## 5.
(a) Build the table with the seed; for spacing k degrees use (cl_tab[8 + k] − cl_tab[8 − k])/(2h) with h = k·π/180. (b) Truncation (h²/6)·12.85; noise 0.005/h; compare. (c) Set (h²/6)·12.85 = 0.005/h and solve for h; convert to degrees. (d) Same stencil on `lift_coefficient(np.radians(deg))` without noise; now only truncation is left and the smallest spacing wins.

## 6.
For each case name the violated assumption (smooth across the stencil; enough digits; the derivative exists at all), then one observable: the value depends on h (kink), the error scales as 1/h (digits), the slope changes sign within one step (stall), the estimate grows without bound as h shrinks (jump), the error is larger than at h = 10⁻⁶ (rounding).

## 7.
(a) dα/dC_L = 1/(dC_L/dα) = cosh²(u)/a₀; at C_L = 1.0, tanh u = 1/C_Lmax so sech²u = 1 − 1/C_Lmax². (b) Exact: α = (C_Lmax/a₀) artanh(C_L/C_Lmax). Bisection: C_L(8°) < 1 < C_L(12°), halve until the bracket is below 10⁻⁶ rad: n = ⌈log₂(bracket/10⁻⁶)⌉. (c) Δα ≈ (dα/dC_L)·0.1; exact Δα = α(1.1) − α(1.0); the difference is the second-order term, which grows as the curve bends towards the plateau. (d) The measured table has noise; use Problem 5's balance to choose the spacing.
