# P01.05  Worked solutions

Seeded problems are worked for **s = 15**. Every number was produced by `solution.py`; if yours differs, run the code before trusting either.

## 1. Your own angle (s = 15)

θ_s = 5 + 7.5 = 12.5° = 0.21817 rad; θ² = 0.047596.

(a) Leading-term errors: sin θ ≈ θ: θ²/6 = **0.79 %**; cos θ ≈ 1: θ²/2 = **2.38 %**; cos θ ≈ 1 − θ²/2: θ⁴/24 = **9.4 × 10⁻⁵**; tan θ ≈ θ: θ²/3 = **1.59 %**. The code's exact errors: 0.798 %, 2.43 %, 9.65 × 10⁻⁵, 1.59 %; the leading term is within 2 % of each.
(b) 12.5° is inside the 1 % range for sin (13.99°) and cos₂ (37.9°), outside for cos (8.07°) and tan (9.91°): **yes, no, yes, no**.
(c) Terms: θ = 0.21817, θ³/6 = 1.73 × 10⁻³, θ⁵/120 = 4.12 × 10⁻⁶, θ⁷/5040 = 4.67 × 10⁻⁹ < 10⁻⁸. The first omitted term must be the fourth, so **3 terms** (θ − θ³/6 + θ⁵/120) give an error below 4.7 × 10⁻⁹. `terms_needed(0.21817, "sin", 1e-8)` returns 3.

## 2. The tools

(a) Sₙ = 1 + r + … + rⁿ; rSₙ = r + … + rⁿ⁺¹; subtract: (1 − r)Sₙ = 1 − rⁿ⁺¹. For |r| < 1, Sₙ → 1/(1 − r), and the tail S − Sₙ = rⁿ⁺¹/(1 − r) exactly (it is itself a geometric series starting at rⁿ⁺¹). For the logarithm: Σ_{k≥n+1} fᵏ/k ≤ (1/(n + 1)) Σ_{k≥n+1} fᵏ = fⁿ⁺¹/((n + 1)(1 − f)), because 1/k ≤ 1/(n + 1) for every k in the tail.

(b) Taylor's theorem: f(x) = Σₖ₌₀ⁿ f⁽ᵏ⁾(x₀)(x − x₀)ᵏ/k! + f⁽ⁿ⁺¹⁾(ξ)(x − x₀)ⁿ⁺¹/(n + 1)!, ξ between x₀ and x. For sin about 0, θ − θ³/6 is T₃ and also T₄ (the θ⁴ coefficient is sin⁽⁴⁾(0)/24 = 0), so the theorem with n = 4 gives |R₄| ≤ max|sin⁽⁵⁾| θ⁵/120 = θ⁵/120 (every derivative of sin is bounded by 1). At θ = π/6 = 0.52360: bound **3.28 × 10⁻⁴**; actual |0.52360 − 0.02392 − 0.5| = **3.26 × 10⁻⁴**. With n = 3 the bound θ⁴/24 = 3.1 × 10⁻³ is true but ten times looser.

(c) C(−½, 0) = 1; C(−½, 1) = −½; C(−½, 2) = (−½)(−3/2)/2 = 3/8; C(−½, 3) = (−½)(−3/2)(−5/2)/6 = −5/16. With x = −M²: (1 − M²)^(−1/2) = 1 + M²/2 + 3M⁴/8 + 5M⁶/16 + …. Ratio test: the k-th term is cₖM²ᵏ with cₖ = (2k)!/(4ᵏk!²), and cₖ₊₁/cₖ = (2k + 1)/(2k + 2); the term ratio is (2k + 1)M²/(2k + 2) → M², below 1 if and only if M < 1: **R = 1**. (The code's estimate from 200 coefficients is 1.0025, approaching 1 from above as (2k + 2)/(2k + 1).)

(d) e^(iθ) = Σ(iθ)ᵏ/k!. Since i⁰, i¹, i², i³ = 1, i, −1, −i and repeat: the even terms are (−1)^(k/2)θᵏ/k! = 1 − θ²/2 + θ⁴/24 − … = cos θ; the odd terms are i(−1)^((k−1)/2)θᵏ/k! = i(θ − θ³/6 + …) = i sin θ. So e^(iθ) = cos θ + i sin θ. Replacing θ by −θ: e^(−iθ) = cos θ − i sin θ. Their product is e⁰ = 1 and also (cos θ + i sin θ)(cos θ − i sin θ) = cos²θ + sin²θ. `euler_formula_residual(1.0, 10)` returns 0.0: with matched truncations the identity holds term by term, not only in the limit.

## 3. The rocket equation as a series

(a) 1/(1 − t) = Σ_{k≥0} tᵏ for |t| < 1; integrating from 0 to f term by term (allowed inside the radius): ln(1/(1 − f)) = Σ_{k≥1} fᵏ/k = f + f²/2 + f³/3 + …, radius **R = 1**. Physically f = 1 means the whole rocket is propellant; the velocity gain diverges logarithmically as f → 1 (and the series' terms stop shrinking), and f > 1 is meaningless.

(b) v_e = 3000 m/s, Δv in m/s:

| f | 1 term | 2 | 3 | 5 | 10 | exact | terms for 1 % |
|---|---|---|---|---|---|---|---|
| 0.1 | 300.0 | 315.0 | 316.0 | 316.1 | 316.1 | 316.1 | 2 |
| 0.5 | 1500 | 1875 | 2000 | 2065.6 | 2079.2 | 2079.4 | 5 |
| 0.9 | 2700 | 3915 | 4644 | 5490 | 6356 | 6907.8 | 24 |

(c) v_e f = (ṁ v_e)(t_burn)/m₀ = (thrust × burn time)/(initial mass): the velocity a body of *constant* mass m₀ would gain from the same impulse. The real rocket gets lighter as it burns, so the same thrust accelerates it more and more; every higher term in the series is positive, so the linear estimate is always low (by 5 % at f = 0.1, 28 % at 0.5, 61 % at 0.9).

(d) Tail bound at f = 0.5 after n = 5 terms: 0.5⁶/(6 × 0.5) = **5.21 × 10⁻³** (times v_e: 15.6 m/s); actual tail 0.69315 − 0.68854 = **4.61 × 10⁻³** (13.8 m/s). The bound is honest and 13 % loose because 1/k keeps falling inside the tail.

## 4. Convergence against cancellation

(a)

| θ | terms for 10⁻¹⁰ | largest term |
|---|---|---|
| 0.25 | 4 | 0.25 |
| 0.5 | 5 | 0.5 |
| 1 | 6 | 1 |
| 2 | 9 | 2 |
| 4 | 13 | 10.7 (k = 1: θ³/6) |
| 8 | 19 | 4.2 × 10² |
| 16 | 31 | 4.1 × 10⁵ |
| 40 | 63 | 1.5 × 10¹⁶ |

The largest term exceeds 10¹⁰ somewhere between θ = 16 and 40 (at θ ≈ 27). Float64 carries 16 digits, so the sum cannot be better than about 2.2 × 10⁻¹⁶ × (largest term) in absolute terms: 10⁻¹⁰ at θ = 16 is just reachable, and at θ = 40 the attainable accuracy is of order 3, i.e. none.

(b) `maclaurin_partial_sum(40.0, "sin", n)` for n = 60, 80, 120 all return **2.0004** (the series has converged in the sense that further terms are below 10⁻¹⁶ of the largest), against `np.sin(40.0)` = 0.74511: an error of **1.26**. The largest term is 40¹⁹/19! = 1.48 × 10¹⁶, and 2.2 × 10⁻¹⁶ × 1.48 × 10¹⁶ = **3.3**: the error is the rounding error of the largest terms, and no number of further terms repairs it because the damage was done when the big terms were rounded.

(c) round(40/2π) = round(6.37) = 6, so θ′ = 40 − 6 × 2π = 40 − 37.699 = **2.3009 rad**. Twelve terms give 0.7451131604793478 against `np.sin(2.3009)` = 0.7451131604793478 and `np.sin(40.0)` = 0.7451131604793488: equal to 10⁻¹⁵ (the 10⁻¹⁵ difference is the rounding of 40 − 37.699 itself). `terms_needed(2.3009, "sin", 1e-15)` guarantees **12** terms.

(d) The exp series at x = −40 with 120 terms returns **−8.0** (negative!) while e⁻⁴⁰ = 4.25 × 10⁻¹⁸; the largest term, 40³⁹/39! = 1.5 × 10¹⁶, has a rounding error 10³³ times the answer. The series for e⁺⁴⁰ has all terms positive, no cancellation, and a relative error of 10⁻¹⁶; its reciprocal is 4.248354255291589 × 10⁻¹⁸, equal to `np.exp(-40.0)` to every digit. Rule: never sum an alternating series whose terms are much larger than its sum; find an identity that gives the same number from a series of positive terms, or reduce the argument.

## 5. Compressible dynamic pressure

(a) With u = 0.2M², (1 + u)^(3.5) = 1 + 3.5u + (3.5 × 2.5/2)u² + (3.5 × 2.5 × 1.5/6)u³ + (3.5 × 2.5 × 1.5 × 0.5/24)u⁴ + … = 1 + 3.5u + 4.375u² + 2.1875u³ + 0.27344u⁴ + …. Subtract 1 and divide by 0.7M² = 3.5u: q_c/q = 1 + 1.25u + 0.625u² + 0.078125u³ + …; with u = 0.2M²: **1 + M²/4 + M⁴/40 + M⁶/1600 + …** (the next term is 0.078125 × 0.0016 M⁸ = M⁸/8000).

(b) `binomial_coefficients(3.5, 4)` = [1, 3.5, 4.375, 2.1875, 0.2734375]. Multiplying by 0.2ᵏ/0.7 for k = 1..4: 1, 0.25, 0.025, 0.000625 = 1, 1/4, 1/40, 1/1600. Confirmed.

(c)

| M | exact q_c/q | 2 terms | 3 terms | 4 terms | pressure error 1 − q/q_c | airspeed error 1 − √(q/q_c) |
|---|---|---|---|---|---|---|
| 0.2 | 1.010040 | 1.0100 | 1.010040 | 1.010040 | 0.99 % | 0.50 % |
| 0.3 | 1.022703 | 1.0225 | 1.022703 | 1.022703 | 2.22 % | 1.12 % |
| 0.5 | 1.064072 | 1.0625 | 1.064063 | 1.064072 | 6.02 % | 3.06 % |
| 0.7 | 1.128575 | 1.1225 | 1.128503 | 1.128576 | 11.4 % | 5.87 % |
| 0.9 | 1.219229 | 1.2025 | 1.218903 | 1.219235 | 18.0 % | 9.44 % |

Two terms are enough to 0.5 % up to M = 0.5; three to 10⁻⁴ up to M = 0.7. An indicator that ignores compressibility over-reads the airspeed by 1 % at M = 0.3 and 3 % at M = 0.5 (the "compressibility correction" of the airspeed tables).

(d) The binomial series in u converges for |u| < 1, i.e. 0.2M² < 1, **M < √5 = 2.24**; in M² the radius is 5. The series is not what stops the formula: at M = 1 a normal shock forms ahead of the pitot probe, the flow into it is no longer isentropic, and q_c is given by the Rayleigh pitot formula instead (Phase 05). A convergent series can describe a physically wrong formula perfectly well.

## 6. Four series that mislead

(a) **Integral test ignored.** Hₙ ≈ ln n + 0.5772 has no limit: H₁₀₆ = 14.39, H₁₀⁷ = 16.70, H₁₀⁸ = 19.00, each decade adding 2.30. A per-term change of 10⁻⁶ at n = 10⁶ is just 1/n, the derivative of ln n, and says nothing about the sum. Correction: compare with ∫dx/x; the series diverges.
(b) **Inside the radius but far too few terms.** Three terms give 1.65; the exact factor at M = 0.9 is 2.294, so the slope is 28 % low; the series needs 16 terms for 1 %. Worse, Prandtl–Glauert itself (a linearised small-perturbation theory) is unreliable above the critical Mach number, 0.7–0.8 for most sections. Correction: use the closed form, and do not trust it at M = 0.9.
(c) **Outside the 1 % range of the approximation.** φ = 60° = 1.047 rad; cos φ ≈ 1 − φ²/2 = 0.452 against the true 0.5, so the "load factor" 2.21 is 11 % high. The series is fine (it converges), the truncation is not: the 1 % range of cos ≈ 1 − φ²/2 ends at 37.9°. Correction: n = 1/cos 60° = 2.00 exactly; use the function.
(d) **Alternating bound without decreasing terms.** At θ = 10 the terms θ, θ³/6, θ⁵/120, θ⁷/5040, θ⁹/9! are 10, 167, 833, 1984, 2756, still growing until k = 4; the bound |S − Sₖ| ≤ |aₖ₊₁| holds only once the terms decrease. The actual error after three terms is |10 − 166.7 + 833.3 − (−0.544)| = 677, not below 1984 by luck of the sign but not guaranteed by anything. Correction: reduce the argument (sin 10 = sin(10 − 4π) = sin(−2.566)) or wait for the terms to turn, as `terms_needed` does (22 terms for 10⁻¹⁰).

## 7. The inertia swing test

(a) (1 − k² sin²φ)^(−1/2) = Σₙ C(−½, n)(−k² sin²φ)ⁿ, and C(−½, n)(−1)ⁿ = (2n)!/(4ⁿ n!²) (the same coefficients as Prandtl–Glauert), so the integrand is Σ [(2n)!/(4ⁿ n!²)] k²ⁿ sin²ⁿφ. Term-by-term integration is legitimate inside |k| < 1: K(k) = Σ [(2n)!/(4ⁿ n!²)] k²ⁿ × (π/2)(2n)!/(4ⁿ n!²) = (π/2) Σ [(2n)!/(4ⁿ n!²)]² k²ⁿ, and T/T₀ = (2/π)K = 1 + k²/4 + 9k⁴/64 + 25k⁶/256 + …. With k² = sin²(θ₀/2) = θ₀²/4 − θ₀⁴/48 + …: 1 + θ₀²/16 + (9/64)(θ₀⁴/16) − θ₀⁴/192 + … = 1 + θ₀²/16 + 11θ₀⁴/3072 + … (SymPy confirms, and the next term is 173θ₀⁶/737280).

(b)

| θ₀ | T/T₀ exact | 2 terms | 3 terms | inertia error 1 − (T₀/T)² |
|---|---|---|---|---|
| 10° | 1.001907 | 1.001899 | 1.001907 | 0.38 % |
| 15° | 1.004301 | 1.004259 | 1.004300 | 0.85 % |
| 23° | 1.010165 | 1.009937 | 1.010159 | 2.00 % |
| 30° | 1.017409 | 1.016747 | 1.017378 | 3.39 % |
| 45° | 1.039973 | 1.036612 | 1.039628 | 7.54 % |

The inertia error reaches **1 % at θ₀ = 16.2°** and **2 % at 23.0°** (the period error reaches 1 % at 22.8°; the inertia error is 2 × (T/T₀ − 1) to first order, which is the "twice" of the lesson).

(c) s = 15: θ₀ = 17.5° = 0.30543 rad. T/T₀ = 1.005858 (exact, by the elliptic integral); assuming T₀ gives I too high by **1.16 %** (1 − 1/1.005858² = 0.01162). Correction factor (T₀/T)² = **0.98838**. Series: 2 terms give 1.005787 (error 7 × 10⁻⁵), 3 terms 1.005857 (1.4 × 10⁻⁶ error, within 10⁻⁵): **3 terms**.

(d) Timing uncertainty 0.2 % in T gives 0.4 % in I. The amplitude bias at 17.5° is 1.16 %, three times larger and *systematic* (always the same sign), so it must be corrected: apply 0.9884 and the residual error of the correction (from an amplitude known to ±1°, dT/dθ₀ ≈ θ₀/8 per radian → ±0.07 % in T, ±0.13 % in I) is well below the timing noise. Reducing the amplitude to below 8° would bring the bias under the timing error but would shorten the measurable swing and worsen the timing; ignoring it is wrong by 1.2 %, which for a flight simulator's roll response is noticeable. What would change the decision: an amplitude the rig cannot measure to ±2° (then reduce the swing instead), or a timing uncertainty above 1 % (then the bias is lost in the noise and correcting it is pointless until the timing improves).
