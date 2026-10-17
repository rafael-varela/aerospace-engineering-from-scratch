# P03.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. Two Newton steps by hand (s = 15: M = 150° = 2.617994 rad, e = 0.74)

**Derivation.** E₀ = M + e sin M = 2.617994 + 0.74 × 0.5 = **2.987994 rad**.
Step 1: f(E₀) = E₀ − e sin E₀ − M = 2.987994 − 0.74 × 0.15306 − 2.617994 = 0.256783; f′(E₀) = 1 − 0.74 cos E₀ = 1 − 0.74 × (−0.98822) = 1.731288; step = −0.148319; **E₁ = 2.839675 rad**.
Step 2: f(E₁) = 1.640 × 10⁻³; f′(E₁) = 1.706528; step = −9.612 × 10⁻⁴; **E₂ = 2.8387135 rad**.

(b) Residual f(E₂) = **1.02 × 10⁻⁷**. Quadratic rule: e₂ ≈ (f″/2f′) e₁² with f″/2f′ = e sin E/(2f′) = 0.0647 and e₁ ≈ |step 2| = 9.6 × 10⁻⁴: e₂ ≈ 0.0647 × 9.2 × 10⁻⁷ = **6.0 × 10⁻⁸ rad**. The code's converged value gives the true error of E₂ as 5.96 × 10⁻⁸: the estimate is right to two figures.

(c) `solve_kepler(2.617994, 0.74)`: E = 2.838713404161357 rad = **162.646°** in **4 steps**; ν = **173.248°**.

(d) `bisection_steps_needed(M − e, M + e, 1e-14)` = ⌈log₂(1.48/10⁻¹⁴)⌉ = **48** midpoints.

**Check.** f(E)/f′(E) ≈ error: 1.02 × 10⁻⁷/1.7 = 6.0 × 10⁻⁸ ✓, the same number from the residual route.

**Common mistakes.** Working in degrees (sin of a degree value is wrong, and E − e sin E mixes units); using E₀ = M (two extra steps); reporting ν = E.

## 2. The four error recurrences

(a) Expand f about x_k: 0 = f(x*) = f(x_k) + (x* − x_k)f′(x_k) + ½(x* − x_k)²f″(ξ). Divide by f′(x_k): 0 = f(x_k)/f′(x_k) + (x* − x_k) + ½(x* − x_k)²f″(ξ)/f′(x_k). Since x_{k+1} = x_k − f(x_k)/f′(x_k), the first two terms are x* − x_{k+1}, so x* − x_{k+1} = −½(x* − x_k)²f″(ξ)/f′(x_k) and **e_{k+1} = |f″(ξ)/(2f′(x_k))| e_k²**. Assumed: f ∈ C² near x*, f′(x*) ≠ 0 so that f′(x_k) is bounded away from zero for x_k near x*, and x_k close enough that the constant stays bounded.

(b) With f = (x − x*)²h: f′ = 2(x − x*)h + (x − x*)²h′, so f/f′ = (x − x*)h/(2h + (x − x*)h′) = (x − x*)/2 × 1/(1 + (x − x*)h′/2h) = (x − x*)/2 + O((x − x*)²). Then x_{k+1} − x* = (x_k − x*) − f/f′ = (x_k − x*)/2 + O(e_k²): **e_{k+1} = e_k/2 + O(e_k²)**, linear with rate ½. For h ≡ 1, f/f′ = (x − x*)²/(2(x − x*)) = (x − x*)/2 exactly, so e_{k+1} = e_k/2 with no remainder; `main()` measures the rate 0.500.

(c) x_{k+1} − x* = g(x_k) − g(x*) = g′(ξ)(x_k − x*) by the mean value theorem, ξ between x_k and x*. If |g′| ≤ ρ < 1 on an interval around x* that contains x₀, every iterate stays in it and e_k ≤ ρᵏe₀ → 0: **linear convergence at rate |g′(x*)|**. If |g′(x*)| > 1 the fixed point repels.

(d) The secant line p through (x_{k−1}, f_{k−1}) and (x_k, f_k) is the linear interpolant; the interpolation error is f(x) − p(x) = ½f″(ξ)(x − x_{k−1})(x − x_k). At x = x*, f = 0, so p(x*) = −½f″(ξ)(x* − x_{k−1})(x* − x_k). The secant iterate is where p = 0, and p is linear with slope f[x_{k−1}, x_k] ≈ f′, so p(x*) = f′ (x* − x_{k+1}). Hence x* − x_{k+1} = −(f″/2f′)(x* − x_{k−1})(x* − x_k), i.e. **e_{k+1} ≈ |f″/2f′| e_k e_{k−1}**. Ansatz e_{k+1} = C e_k^q, so e_k = C e_{k−1}^q and e_{k−1} = (e_k/C)^{1/q}. Substituting: C e_k^q = K e_k (e_k/C)^{1/q}, so q = 1 + 1/q, **q² − q − 1 = 0, q = (1 + √5)/2 = 1.6180**.

**Check.** The fitted orders in the lesson (2.00, 0.500 rate, 1.00 at rate e|cos E*|, 1.53–1.71) match all four.

**Common mistakes.** In (a) dropping the second-order remainder *and* then claiming quadratic convergence (the remainder is the convergence); in (d) assuming e_{k+1} ∝ e_k² because "two points are used".

## 3. A hyperbolic flyby (s = 15: e = 1.35, M_h = 2.0)

(a) f(H) = e sinh H − H − M_h has f′(H) = e cosh H − 1 ≥ e − 1 = 0.35 > 0, so f is strictly increasing on ℝ and has at most one root. f(0) = −M_h < 0; for H ≥ 0, sinh H ≥ H gives f(H) ≥ (e − 1)H − M_h > 0 at H = M_h/(e − 1) = 5.71. Bracket **[0, 5.71]** (any upper end beyond the root works). Exactly one root by the intermediate value theorem.

(b) From H₀ = asinh(M_h/e) = asinh(1.4815) = 1.1844: **H = 1.744296296284** in **6 steps**. From H₀ = M_h = 2.0: 5 steps (the start is closer: f(2) = 2.9 against f(1.18) = −1.2). `brentq` on [0, 10] gives 1.7442962962835 ✓.

(c) ν = 2 atan(√(2.35/0.35) tanh(0.8721)) = 2 atan(2.5912 × 0.7026) = **122.43°**. Asymptote ν_∞ = arccos(−1/1.35) = **137.79°**: the spacecraft at M_h = 2 is 15° short of its outbound asymptote.

**Check.** e sinh H − H = 1.35 × 2.7465 − 1.7443 = 2.0000 ✓. ν < ν_∞ as required for a finite time ✓.

**Common mistakes.** Using sin and cos (the elliptic formulas) with e > 1 and getting a complex √(1 − e²); forgetting the sign change M_h = e sinh H − H (hyperbolic) against M = E − e sin E (elliptic).

## 4. Measured orders (s = 15: e = 0.65, M = 20°; E* = 0.8277197 rad from `brentq`)

| method | steps | evaluations | fitted order (observed) | theory | digits per evaluation |
|---|---|---|---|---|---|
| Newton | 5 | 10 | 2.00 | 2 | 1.51 |
| secant | 8 | 8 | 1.62 | 1.618 | 1.48 |
| fixed point | 39 | 39 | 1.00; rate 0.4398 | 1; e|cos E*| = 0.4398 | 0.35 |
| bisection | 46 | 49 | 0.98 (bound halves: rate ½) | bound rate ½ | 0.28 |

(a), (b) as tabulated. The fixed-point rate matches e|cos E*| to four figures; bisection's 46 midpoints are within one of the bound's 47 (the stopping test is on the half-width, which happened to fall below 10⁻¹⁴ a step early).

(c) **Digits per evaluation** = log₁₀(e₀/e_final)/evaluations: Newton 1.51 and secant 1.48 are tied; fixed point 0.35 and bisection 0.28 are five times worse. If f′ cost ten times f, a Newton step would cost 11 units against the secant's 1: Newton's digits per unit cost would fall to 0.27 and the **secant method would win** outright. For Kepler's equation f′ is one cosine, so Newton and secant are equivalent and Newton's simpler start decides it.

**Check.** 10 evaluations × 1.51 digits = 15.1 digits = log₁₀(0.256/1.1 × 10⁻¹⁶) ✓.

**Common mistakes.** Comparing steps instead of evaluations (Newton "wins" by steps but costs two per step); fitting the order over errors at the rounding floor (the slope collapses toward 0).

## 5. GTO radius crossing (s = 15: R = 22 500 km; a = 24 396.1 km, e = 0.7283, T = 37 922 s)

(a) g(t) = ‖r(t)‖ − R. Bisection on [0, T/2] to 10⁻⁶ s: **t = 4464.594507 s** (74.41 min) in 35 midpoints = 37 Kepler solves. Secant from (1000, 2000): **t = 4464.594507 s** in 7 steps = 8 solves. Difference 2 × 10⁻⁸ s, below the tolerance. The secant method used a fifth of the solves.

(b) cos E = (1 − R/a)/e = (1 − 0.92228)/0.7283 = 0.10672, E = 83.874°; M = E − e sin E = 1.46384 − 0.72414 = 0.73970 rad; t = M/n = 0.73970/1.6569 × 10⁻⁴ = **4464.5945 s**. The root finders agree with the closed form to 2 × 10⁻⁸ s (bisection) and 10⁻⁹ s (secant).

(c) The second crossing is at T − t = 37 922.1 − 4464.6 = **33 457.5 s** (557.6 min) by the symmetry of the orbit about the apse line. On [0, T] both endpoints are perigee, where g = a(1 − e) − R < 0 at both ends: no sign change, bisection raises `ValueError`, and a bracket that did contain both roots (an even number) hides them from any sign-change test. Bracket each monotone half separately.

**Check.** r at t = 4464.6 s from `orbit_state` is 22 500.000 km ✓; ν there is 132.4° against E = 83.9°.

**Common mistakes.** Bracketing over a full period; using geometric altitude instead of radius for R.

## 6. Specifying the flight-computer solver

(a) **Newton from E₀ = M + e sin M**, cap 10. Per step it costs one sin and one cos (2 evaluations); at e = 0.95 the worst case in the lesson's sweep is 8 steps = 16 evaluations, over budget, so pair it with the start E₀ = π for e > 0.9 (8 steps at M = 5°, 4 to 7 for M from 150° down to 30°, but 10 at M = 1°) or accept a 12-evaluation cap with the fallback below. Rejected: bisection (48 evaluations, 4× budget); fixed point (rate e|cos E*| → 0.95 as M → 0: about 600 steps); secant (7 to 11 evaluations, within budget, but needs two starts and can step outside the bracket on the steep part of f near perigee at high e; a safeguarded secant is Brent's method, which is the second choice).

(b) Stop when |E_{k+1} − E_k| ≤ 4 ulp(E) ≈ 10⁻¹⁵ rad (for Newton the error is then ≈ 10⁻³⁰, the floor); cap at 10 steps; if the cap is hit, run bisection on [M − e, M + e] (always valid, 48 steps) and raise a flag so the event is logged. A solver that "must never return a wrong E" needs the bracketing fallback, not a bigger cap.

(c) At e = 0.999 the standard start overshoots to 316° and Newton needs 13 steps (26 evaluations); from E₀ = π it needs 8 (16 evaluations). Cost: the budget rises to 16 evaluations per solve, or the rate drops to 7 Hz, or a better initial guess (a cubic in M near perigee) brings it back to 10. Also f′ = 1 − e cos E ≥ 10⁻³ near perigee, so the step f/f′ amplifies the residual's rounding error by 10³: the stopping test must tolerate step sizes of 10⁻¹³ rad there.

(d) |f(E)| = |f′(ξ)| × |E − E*| with f′ ∈ [1 − e, 1 + e] = [0.05, 1.95] at e = 0.95. A residual of 10⁻¹² therefore guarantees an error of at most 2 × 10⁻¹¹ rad (near perigee) and at least 5 × 10⁻¹³ rad (near apogee): **worse** as a guarantee than the step test by a factor 20 at perigee, and not uniform around the orbit. Equivalent only at e = 0. Keep the step test; report the residual as a diagnostic.

**Common mistakes.** Choosing the fastest method without a fallback; a cap so high that a stalled solve eats the 100 ms frame.

## 7. The Kepler map is a contraction

(a) g(E₁) − g(E₂) = e(sin E₁ − sin E₂) = e cos ξ (E₁ − E₂) for some ξ between E₁ and E₂ (mean value theorem), so |g(E₁) − g(E₂)| ≤ e |E₁ − E₂| since |cos ξ| ≤ 1. ∎

(b) *Uniqueness*: if E₁ and E₂ are both fixed points, |E₁ − E₂| = |g(E₁) − g(E₂)| ≤ e|E₁ − E₂|, which for e < 1 forces E₁ = E₂. *Existence and convergence*: |E_{k+1} − E_k| = |g(E_k) − g(E_{k−1})| ≤ e|E_k − E_{k−1}| ≤ eᵏ|E₁ − E₀|, so Σ|E_{k+1} − E_k| converges (geometric series) and (E_k) is Cauchy; its limit E* satisfies E* = g(E*) by continuity. Error bound by induction: |E_{k+1} − E*| = |g(E_k) − g(E*)| ≤ e|E_k − E*|, hence **|E_k − E*| ≤ eᵏ|E₀ − E*|**, from any E₀. ∎

(c) For the orbit of Problem 1 (e = 0.74, M = 150°): |E₀ − E*| = 0.1493 rad, so k ≥ log(10⁻¹²/0.1493)/log 0.74 = **85.4 → 86 steps** guaranteed. The code (`solve_kepler(M, 0.74, "fixed_point", tol=1e-12)`) takes **77 steps**. The bound uses the worst rate e at every step; the actual rate is e|cos E*| = 0.74 × |cos 162.6°| = 0.706, which predicts log(10⁻¹²/0.1493)/log 0.706 = 79 steps, within two of the observation (the stopping test on the step rather than the error accounts for the rest).

**Check.** The lesson's e = 0.73, M = 60° case has rate 0.1399 and took 17 steps; the e-rate bound would have promised 108. Both bounds hold; the e|cos E*| estimate is the useful one.

**Common mistakes.** Claiming convergence "because e < 1" without the mean value step; forgetting that the bound needs the error at E₀, not the first step size (they differ by the factor 1/(1 − e)).

## Check-yourself answers

1. e_k ≤ (b − a)/2ᵏ after k midpoints; k = ⌈log₂(w/tol)⌉ midpoints guarantee e_k ≤ tol.
2. g(E) = M + e sin E has |g′| = e|cos E| ≤ e < 1 everywhere, so it is a contraction: convergence from any start at rate e|cos E*| (asymptotically), bounded by e.
3. The tangent at x_k is y = f(x_k) + f′(x_k)(x − x_k); its zero is x_{k+1} = x_k − f/f′. Order 2 at a simple root (e_{k+1} ≈ |f″/2f′| e_k²); linear with rate ½ at a double root.
4. e_{k+1} ≈ C e_k e_{k−1}; with e_{k+1} = C e_k^q this requires q = 1 + 1/q, so q = (1 + √5)/2. Per evaluation the secant method gains 1.618 orders per step against Newton's 2 per two evaluations (√2 = 1.414 per evaluation), so it is cheaper whenever f′ costs about as much as f or more.
5. Zero derivative (x² − 1 from 0); divergence (atan x from 1.5); a cycle (x³ − 2x + 2 from 0); and, as a fourth, degraded linear convergence at a multiple root ((x − 1)²).
6. [M − e, M + e], because |E − M| = e|sin E| ≤ e; f′ = 1 − e cos E ≥ 1 − e > 0.
7. tan(ν/2) = √(1.5/0.5) tan 45° = √3, so ν = 120°; r = a(1 − e cos 90°) = a.
8. Vis-viva v² = μ(2/r − 1/a); |r × v| = √(μa(1 − e²)); ε = v²/2 − μ/r = −μ/2a. All three should hold to about 10⁻¹² relative (a few ulps through the square roots and divisions); a larger violation means a formula error, not rounding.
