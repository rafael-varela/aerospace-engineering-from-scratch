# P01.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; symbolic results are exact, numerical ones are given to four significant figures and come from `solution.py`.

## 1. Slope of your own lift curve

**Assumptions.** The lesson's tanh model; a₀ = 2π rad⁻¹; angles in radians inside the model.

**Derivation** (s = 15: C_Lmax = 1.5, α = 7.5° = 0.130900 rad). (a) u = 6.28319 × 0.130900/1.5 = 0.5483. (b) tanh u = 0.49925, C_L = 1.5 × 0.49925 = 0.7489. (c) dC_L/dα = a₀(1 − tanh²u) = 6.28319 × (1 − 0.24926) = 4.717 rad⁻¹ = 4.717 × π/180 = 0.08233 deg⁻¹. (d) Fraction remaining = 1 − tanh²u = 0.7507.

**Result.** u = 0.5483, C_L = 0.7489, dC_L/dα = 4.717 rad⁻¹ (0.08233 deg⁻¹), 75.1 % of the thin-airfoil slope.

**Check.** `central_difference(lambda x: lift_coefficient(x, cl_max=1.5), np.radians(7.5), 1e-4)` returns 4.717076, differing from the exact 4.717076 by 7 × 10⁻⁸, which is the (h²/6)|C_L‴| order expected. Limiting case: with α → 0 the fraction → 1, and for s = 31 (C_Lmax = 1.82, α = 15.5°) the fraction is 0.46: twice the angle costs far more than twice the slope, even with a higher C_Lmax.

**Common mistakes.** Feeding degrees into tanh; multiplying by 180/π instead of π/180 for the per-degree slope; reporting the fraction as C_L/(a₀α) (that is the secant ratio, 0.91 here, not the tangent ratio).

## 2. The error of each formula

**Derivation.** (a) Forward: f(x+h) − f(x) = h f′ + (h²/2)f″ + (h³/6)f‴ + …; divide by h: f′ + (h/2)f″ + O(h²). First order, leading error (h/2)f″. Central: f(x+h) − f(x−h) = 2h f′ + (h³/3)f‴ + (h⁵/60)f⁽⁵⁾ + … (all even-order terms cancel); divide by 2h: f′ + (h²/6)f‴ + (h⁴/120)f⁽⁵⁾ + …. Second order, leading error (h²/6)f‴.

(b) D(2h) = f′ + (4h²/6)f‴ + (16h⁴/120)f⁽⁵⁾ + …. Then [4D(h) − D(2h)]/3 = f′ + [(4 − 16)/120](h⁴/3)f⁽⁵⁾ + … = f′ − (h⁴/30)f⁽⁵⁾ + …. Collecting function values: 4D(h) = [2f(x+h) − 2f(x−h)]/h and D(2h) = [f(x+2h) − f(x−2h)]/(4h); so [4D(h) − D(2h)]/3 = [−f(x+2h) + 8f(x+h) − 8f(x−h) + f(x−2h)]/(12h). The five-point stencil, leading error −(h⁴/30)f⁽⁵⁾, fourth order.

(c) Replace h by h/2 in (b): [4D(h/2) − D(h)]/3 = [−f(x+h) + 8f(x+h/2) − 8f(x−h/2) + f(x−h)]/(6h), which is the five-point stencil with step h/2 (its denominator is 12 × h/2 = 6h). Identical, not merely equal to O(h⁴).

(d) Forward is exact when f″ = f‴ = … = 0: polynomials of degree ≤ 1. Central is exact when f‴ and higher odd derivatives vanish: degree ≤ 2. Five-point is exact when f⁽⁵⁾ and higher vanish: degree ≤ 4 (the f‴ term was cancelled by construction).

**Check.** The tests `test_central_difference_is_exact_for_quadratics_and_five_point_for_quartics` and `test_richardson_equals_five_point_with_half_the_step` are (d) and (c) in code. Dimensional check: every error term has the units of f′ because hⁿ f⁽ⁿ⁺¹⁾ has the units of f/x.

**Common mistakes.** Forgetting that the even terms cancel only when the two points are *symmetric* about x; dividing D(2h) by 2h instead of 4h; claiming the central difference is exact for cubics (it is not: f‴ ≠ 0).

## 3. Mach number in a climb

**Assumptions.** ISA troposphere, γ = 1.4, R = 287.05287 J/(kg K); constant V = 150 m/s and ḣ = 10 m/s; geopotential altitude.

**Derivation.** (a) M = V/√(γRT), so dM/dT = −V/(2T√(γRT)) = −M/(2T). Chain rule: dM/dt = (dM/dT)(dT/dh)(dh/dt) = −(M/(2T))(dT/dh)ḣ.

(b) At 5000 m: T = 255.65 K, a = 320.53 m/s, M = 0.4680; dM/dt = −(0.4680/(2 × 255.65))(−0.0065)(10) = 5.949 × 10⁻⁵ s⁻¹. At 10 000 m: T = 223.15 K, a = 299.46 m/s, M = 0.5009; dM/dt = 7.295 × 10⁻⁵ s⁻¹. M rises (about 0.0036 per minute at 5 km) because the air cools with altitude, the speed of sound falls, and the same true airspeed is a larger fraction of it.

(c) `central_difference(M, 5000.0, 1.0) * 10.0` = 5.949 × 10⁻⁵; at 10 000 m, 7.295 × 10⁻⁵. Agreement to four figures: M(h) is nearly linear over 1 m, so the h² term is negligible.

(d) `one_sided_derivatives(M, 11000.0, 10.0)` returns (7.626 × 10⁻⁶ m⁻¹, 0). Times ḣ: dM/dt jumps from 7.63 × 10⁻⁵ s⁻¹ to 0 at the tropopause, because T, hence a, is constant above 11 km and a constant-V climb no longer changes M (M = 0.5084 from there up).

**Check.** Units: (1/K)(K/m)(m/s) = 1/s. ✓ Limiting case: ḣ = 0 gives dM/dt = 0; a level flight at constant V has constant M. ✓

**Common mistakes.** Differentiating a = √(γRT) as if a ∝ T (factor 2 error); using the lapse rate in K/km with ḣ in m/s; expecting the central difference *at* 11 km to return either one-sided value (it returns the average, 3.8 × 10⁻⁶).

## 4. The error curve for two functions

**Result** (author's run; rounding-dominated entries vary by a factor of 2–3 between machines).

sin x at x = 1 (f = 0.8415, f‴ = −0.5403):

| h | forward | central | five-point |
|---|---|---|---|
| 10⁻¹ | 4.3 × 10⁻² | 9.0 × 10⁻⁴ | 1.8 × 10⁻⁶ |
| 10⁻² | 4.2 × 10⁻³ | 9.0 × 10⁻⁶ | 1.8 × 10⁻¹⁰ |
| 10⁻³ | 4.2 × 10⁻⁴ | 9.0 × 10⁻⁸ | 7.6 × 10⁻¹⁴ |
| 10⁻⁴ | 4.2 × 10⁻⁵ | 9.0 × 10⁻¹⁰ | 4.2 × 10⁻¹³ |
| 10⁻⁵ | 4.2 × 10⁻⁶ | 1.1 × 10⁻¹¹ | 3.9 × 10⁻¹⁴ |
| 10⁻⁸ | 3.0 × 10⁻⁹ | 2.6 × 10⁻⁹ | 2.6 × 10⁻⁹ |
| 10⁻¹² | 4.3 × 10⁻⁵ | 1.2 × 10⁻⁵ | 1.5 × 10⁻⁵ |

Lift model at α = 16° (f = 1.1889, f‴ = 82.12): central errors 1.3 × 10⁻¹, 1.4 × 10⁻³, 1.4 × 10⁻⁵, 1.4 × 10⁻⁷, 1.4 × 10⁻⁹ for h = 10⁻¹ … 10⁻⁵, minimum 6.6 × 10⁻¹¹ at h = 10⁻⁶, then 2.9 × 10⁻¹⁰, 9.2 × 10⁻⁹, …, 6.4 × 10⁻⁵ at 10⁻¹².

(a) Fitted orders for sin on h = 10⁻² … 10⁻⁴: forward 1.00, central 2.00; five-point on 10⁻¹ … 10⁻² (it is at 10⁻¹³ by 10⁻³): 4.00. For the lift model on the same windows: 1.00, 2.00, 4.0.

(b) sin: h_opt = (3 × 2.2 × 10⁻¹⁶ × 0.8415/0.5403)^(1/3) = 1.0 × 10⁻⁵; measured minimum at 10⁻⁵ (1.1 × 10⁻¹¹). Lift at 16°: h_opt = (3ε × 1.1889/82.12)^(1/3) = 2.1 × 10⁻⁶; measured minimum at 10⁻⁶. Both within the decade.

(c) h_opt ∝ (f/f‴)^(1/3). For sin, f/f‴ = 1.56; for the lift model at 16°, 0.0145: a ratio of 108, whose cube root is 4.8, which is the ratio of the two h_opt values (1.0 × 10⁻⁵/2.1 × 10⁻⁶ = 4.8). The lift curve bends hard near the plateau (large f‴), so truncation stays significant down to a smaller h.

(d) Complex step at 16°: error 1.4 × 10⁻⁷ at h = 10⁻⁴ (the same h² truncation as the central difference, since Im f(x+ih)/h = f′ − (h²/6)f‴ + …), then 1.6 × 10⁻¹⁵ at 10⁻⁸ and 2 × 10⁻¹⁶ at 10⁻¹² and at 10⁻²⁰. There is no subtraction of nearly equal numbers: the imaginary part is built from terms of size h f′ directly, so making h tiny loses no digits; the only limit is the precision of f′ itself.

**Check.** The central error at h = 10⁻² for sin is 9.0 × 10⁻⁶ and (h²/6)|f‴| = 10⁻⁴ × 0.5403/6 = 9.0 × 10⁻⁶. ✓

**Common mistakes.** Fitting the five-point order on a window that includes its rounding floor (the order comes out near 3.4 on 10⁻² … 10⁻⁴); claiming a complex-step error of exactly zero (it is a few ε).

## 5. Derivative from a table with noise

**Assumptions.** Noise uniform in ±0.005, independent between points; exact slope at 8° is 4.3427 rad⁻¹; C_L‴ = −12.85 rad⁻³.

**Result** (s = 15):

| spacing | h (rad) | estimate | error | truncation (h²/6)·12.85 | noise bound 0.005/h |
|---|---|---|---|---|---|
| 1° | 0.01745 | 4.523 | +0.180 | 6.5 × 10⁻⁴ | 0.286 |
| 2° | 0.03491 | 4.349 | +0.006 | 2.6 × 10⁻³ | 0.143 |
| 4° | 0.06981 | 4.328 | −0.015 | 1.0 × 10⁻² | 0.072 |

(b) Noise dominates at every spacing in the table; at 1° it is 440 times the truncation error. The 2° estimate happens to be good for this seed; another seed can put it anywhere within ±0.14.

(c) Set (h²/6)·12.85 = 0.005/h: h³ = 0.03/12.85, h = 0.133 rad = 7.6°. That is 10⁵ times larger than the smooth-function optimum of about 10⁻⁶ rad. With real data the "small h" instinct is exactly wrong; the noise, not machine epsilon, sets the floor, and it sets it 10 orders of magnitude higher.

(d) Without noise: errors −6.7 × 10⁻⁴ (1°), −2.8 × 10⁻³ (2°), −1.4 × 10⁻² (4°): now the smallest spacing wins and the errors follow h² (factors 4.2 and 4.8 per doubling, close to 4; the excess is the h⁴ term growing with h).

**Check.** The noise-free errors agree with (h²/6)·12.85 to 3 % at 1°, confirming C_L‴. The sign is negative because C_L‴ < 0 at 8°.

**Common mistakes.** Using 0.005/(2h) as the noise bound (two independent errors can add); treating the lucky 2° result as evidence that 2° is the right spacing; differentiating the table at 1° and then smoothing the *derivative* instead of the data.

## 6. When a numerical derivative lies

(a) **Kink.** The function is not differentiable, so the stencil returns the average of the two one-sided slopes: −0.00325 K/m at 11 km, independent of h. Sign: the answer does not change as h shrinks, yet `one_sided_derivatives` gives two different values. Do: report both one-sided derivatives, or move the evaluation point off the kink.

(b) **Three significant figures.** The table has a rounding error of about 5 × 10⁻⁴ per entry; the derivative error is about 5 × 10⁻⁴/h, so with h = 0.01 the slope is uncertain by ±0.05 rad⁻¹ (0.8 % of a₀), and with h = 10⁻³ by ±0.5. Sign: the error scales as 1/h as you refine. Do: use a spacing balanced as in Problem 5 (about 4° for this precision), or fit a smooth curve first.

(c) **Just below a sharp stall.** The stencil straddles a slope discontinuity, so it mixes a positive pre-stall slope with a negative post-stall one and can return almost any number, including zero. Sign: the value swings with h, and the one-sided derivatives differ in sign. Do: use the backward (one-sided) difference from the pre-stall side and say so.

(d) **Shock.** The function jumps by Δp; the central difference returns Δp/(2h), which grows without bound as h → 0. This is not a derivative problem: the derivative does not exist and the physical quantity is the jump itself (Rankine–Hugoniot, Phase 11). Sign: the estimate doubles every time h halves. Do: measure Δp and the shock location, not a slope.

(e) **h = 10⁻¹² rad.** Rounding: ε·f/h ≈ 2.2 × 10⁻¹⁶ × 0.78/10⁻¹² = 1.7 × 10⁻⁴, and the measured error at 8° is 7.3 × 10⁻⁵, far above the 6 × 10⁻¹¹ reachable at 10⁻⁶. Sign: the error at 10⁻¹² is larger than at 10⁻⁶, and all four formulas return the same wrong value. Do: use h near (3ε f/f‴)^(1/3), or the complex step.

**Common mistakes.** Treating a stable-looking answer at a kink as correct; "fixing" a shock by smoothing it; choosing h by how small it looks rather than by an error estimate.

## 7. The inverse problem

**Assumptions.** Lesson model with C_Lmax = 1.4, a₀ = 2π; the inverse α(C_L) exists because C_L(α) is strictly increasing.

**Derivation.** (a) dα/dC_L = 1/(dC_L/dα) = cosh²(u)/a₀. At C_L = 1.0, tanh u = 1/1.4 = 0.7143, sech²u = 1 − 0.5102 = 0.4898, dC_L/dα = 6.2832 × 0.4898 = 3.077 rad⁻¹, so dα/dC_L = 0.3249 rad per unit C_L = 18.62° per unit C_L.

(b) Exact: α = (C_Lmax/a₀) artanh(C_L/C_Lmax) = (1.4/6.2832) × artanh(0.7143) = 0.22282 × 0.89588 = 0.19962 rad = 11.44°. Bisection: C_L(8°) = 0.7780 < 1 < C_L(12°) = 1.0293, so by the intermediate value theorem (C_L continuous) a root lies in (8°, 12°), a bracket of 4° = 0.0698 rad. Halving to 10⁻⁶ rad needs ⌈log₂(0.0698/10⁻⁶)⌉ = 17 steps (18 if you start the count from the 0.2 rad bracket (0.1, 0.3)); the code converges to 0.199617 rad.

(c) Linear: Δα ≈ 0.3249 × 0.1 = 0.03249 rad = 1.862°. Exact: α(1.1) − α(1.0) = (0.22282 × artanh(0.78571)) − 0.19962 = 0.23622 − 0.19962 = 0.03660 rad = 2.097°. The linear estimate is 11 % low. Near C_L = 1.0 the slope dC_L/dα is already falling fast (d²C_L/dα² is large and negative), so the inverse slope dα/dC_L rises steeply over the interval; at C_L = 0.5 the curve is nearly straight (slope 5.48 rad⁻¹ against 6.28 at zero) and the same ΔC_L = 0.1 is linearised to within 3 %.

(d) The certifying engineer wants dα/dC_L to convert a required lift margin into an angle-of-attack margin before stall; on a measured curve with ±0.005 noise they should use the spacing from Problem 5(c), about 7°, or fit a smooth curve through the data and differentiate the fit.

**Check.** Inverse-function consistency: 3.077 × 0.3249 = 1.000. ✓ Limiting case: at C_L → 0, dα/dC_L → 1/a₀ = 0.159 rad = 9.1° per unit C_L, which is the thin-airfoil "9° per unit C_L" rule of thumb. ✓

**Common mistakes.** Inverting the slope at the wrong point (using dC_L/dα at α = 0 instead of at C_L = 1.0, which gives 9.1° instead of 18.6°); forgetting that bisection needs the IVT, i.e. continuity, to guarantee a root; reporting Δα for ΔC_L = +0.1 as the same size as for ΔC_L = −0.1 (the curve is not symmetric about C_L = 1.0).

## Check-yourself answers

1. f′(x) = lim_{h→0} [f(x+h) − f(x)]/h. It never evaluates anything *at* h = 0; the value of the quotient at h = 0 is undefined (0/0), and the limit is what the quotient approaches.
2. No. The one-sided derivatives are −0.0065 K/m from below and 0 K/m from above; they differ, so the derivative does not exist, although the function is continuous (216.65 K from both sides).
3. dC_L/dα = a₀ sech²(a₀α/C_Lmax) = a₀(1 − tanh²(a₀α/C_Lmax)). At α = 0 the slope is a₀ = 2π = 6.283 rad⁻¹ = 0.1097 deg⁻¹.
4. In the central difference the two points are symmetric about x, so f(x+h) − f(x−h) cancels every even-order term of the Taylor expansion, including (h²/2)f″; the forward difference uses f(x) itself and nothing cancels the f″ term.
5. Order = log₁₀(2.2 × 10⁻⁴/2.1 × 10⁻⁶) = 2.0, so the formula is second order and the error at h = 10⁻⁴ should be another factor of 100 smaller: about 2 × 10⁻⁸ (the author's run gives 2.14 × 10⁻⁸).
6. Each function value is rounded to about 16 digits; subtracting two values that agree to 12 digits leaves 4, and dividing by 10⁻¹² amplifies the rounding to ε f/h ≈ 10⁻⁴. At h = 10⁻⁶ the rounding contribution is 10⁻¹⁰ and the truncation 10⁻¹²; the boundary is set by machine epsilon ε = 2⁻⁵² ≈ 2.2 × 10⁻¹⁶ through h_opt = (3εf/f‴)^(1/3).
7. The five-point (fourth-order) stencil with step h/2: same four function values, same weights, identically.
8. The linear estimate is the tangent line; its error is about ½(d²C_L/dα²)Δα². At α = 0 the second derivative is zero (the curve is odd), so the error is third order and small (7 %); at 8° the second derivative is −21.7 rad⁻² and the gust carries the wing to 13.7°, where the curve is bending hard towards the plateau, so the tangent overshoots by 32 %.
