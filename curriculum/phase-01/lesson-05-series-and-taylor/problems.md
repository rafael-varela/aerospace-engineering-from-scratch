# P01.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units where there are any, and for every series give the number of terms used and the bound that justifies it.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your own angle)

Let θ_s = 5 + s/2 degrees. By hand, using the leading-term error laws of lesson §2.4: (a) the relative error of sin θ ≈ θ, cos θ ≈ 1, cos θ ≈ 1 − θ²/2 and tan θ ≈ θ at θ_s; (b) whether θ_s is inside the 1 % range of each; (c) the number of non-zero terms of the sin series needed for an absolute error below 10⁻⁸ at θ_s, from the alternating-series bound. Check (a) with `small_angle_error` and (c) with `terms_needed`. Answer form: four errors, four yes/no, one term count, and the four exact errors from the code beside the hand ones.

## 2. (Derivation: the tools)

(a) Derive the sum of the geometric series and the exact tail after n + 1 terms. Then show that the tail of ln(1/(1 − f)) = Σ fᵏ/k after n terms is below fⁿ⁺¹/((n + 1)(1 − f)).
(b) Write Taylor's theorem with the Lagrange remainder and use it to bound the error of sin θ ≈ θ − θ³/6 at θ = 30°, choosing the degree n so that the bound is as tight as the theorem allows. Compare with the actual error.
(c) Derive the first three terms of the binomial series of (1 − M²)^(−1/2) from the definition of C(α, k), and show by the ratio test that its radius of convergence in M is 1.
(d) From the Maclaurin series of eˣ with x = iθ, separate real and imaginary parts and obtain Euler's formula; then show that e^(iθ)e^(−iθ) = 1 gives cos²θ + sin²θ = 1.

Answer form: four derivations; (b) with two numbers.

## 3. (Unfamiliar scenario: the rocket equation as a series)

The ideal rocket equation gives the velocity gain Δv = v_e ln(m₀/m₁) = v_e ln(1/(1 − f)), with v_e the exhaust velocity and f = 1 − m₁/m₀ the propellant mass fraction.

(a) Expand ln(1/(1 − f)) as a power series in f (integrate the geometric series term by term) and state its radius of convergence. What does the radius mean for a rocket?
(b) With v_e = 3000 m/s, compute Δv by the series with 1, 2, 3, 5 and 10 terms for f = 0.1, 0.5 and 0.9, and compare with the closed form. How many terms does each need for 1 %?
(c) The one-term result Δv ≈ v_e f is the "linear" rocket: what is it physically (think impulse = thrust × time, with constant mass), and why is it an *under*estimate?
(d) Use the tail bound from Problem 2(a) to check the 5-term result at f = 0.5 and compare with the actual tail.

Answer form: a series with its radius and a sentence; a 3 × 5 table plus three counts; two sentences; two numbers.

## 4. (Numerical, uses your code: convergence against cancellation)

(a) Tabulate `terms_needed(θ, "sin", 1e-10)` and the size of the largest term for θ = 0.25, 0.5, 1, 2, 4, 8, 16, 40. At what θ does the largest term exceed 10¹⁰, and what does that imply for the achievable accuracy in float64?
(b) Evaluate the sin series at θ = 40 with 60, 80 and 120 terms and compare with `np.sin(40.0)`. Explain the error with one number (the largest term times machine epsilon).
(c) Reduce the argument: θ′ = 40 − 2π·round(40/2π). Evaluate the series at θ′ with 12 terms and compare. How many terms does `terms_needed` guarantee for 10⁻¹⁵ at θ′?
(d) Repeat for e⁻⁴⁰ by the exp series (120 terms) and by 1/(series for e⁴⁰). State which is right and why the other is not.

Answer form: a table and a sentence; three numbers and an explanation; two numbers and a count; two numbers and a sentence.

## 5. (Compressible dynamic pressure: a binomial series you will meet in Phase 05)

An airspeed indicator measures the impact pressure q_c = p[(1 + 0.2M²)^(3.5) − 1] (isentropic flow, γ = 1.4), while the incompressible dynamic pressure is q = ½ρV² = 0.7 p M².

(a) Expand (1 + 0.2M²)^(3.5) by the binomial series to the M⁶ term and show that q_c/q = 1 + M²/4 + M⁴/40 + M⁶/1600 + ….
(b) Use `binomial_coefficients(3.5, 4)` to confirm the four coefficients.
(c) At M = 0.2, 0.3, 0.5, 0.7 and 0.9 compute q_c/q exactly and by two, three and four terms, and the error an instrument makes if it assumes q_c = q: in pressure and in the indicated airspeed (which goes as √q).
(d) Why does this series converge for every M (what is its radius in M²), unlike Prandtl–Glauert's? What stops the formula itself at M = 1?

Answer form: one derivation, one check, a table with two error columns, two sentences.

## 6. (Judgment: four series that mislead)

For each case, say what is wrong with the reasoning, give one number, and say what to do instead:

(a) "The partial sums of Σ1/k have reached 14.39 after a million terms and change by only 10⁻⁶ per term, so the sum is about 14.4."
(b) "The lift-curve slope at M = 0.9 is the incompressible slope times 1 + M²/2 + 3M⁴/8 = 1.65."
(c) "A 60° banked turn has load factor 1/cos φ ≈ 1/(1 − φ²/2) = 2.21, so the wing must carry 2.21 W."
(d) "The sin series alternates, so at θ = 10 the error after three terms is below the fourth term, 10⁷/5040 ≈ 1984."

Answer form: four short paragraphs, each with a number and a correction.

## 7. (Hard, uses your code: the inertia swing test and the pendulum series)

The swing test measures a body's moment of inertia from the period of a compound pendulum, I = T²·(m g d)/(4π²), where d is the distance from the pivot to the centre of mass; the formula assumes T = T₀, the small-amplitude period.

(a) Expand the integrand of K(k) = ∫₀^(π/2) dφ/√(1 − k² sin²φ) by the binomial series and integrate term by term using ∫₀^(π/2) sin²ⁿφ dφ = (π/2)(2n)!/(4ⁿ n!²); show that T/T₀ = (2/π)K(k) = Σ [(2n)!/(4ⁿ n!²)]² k²ⁿ with k = sin(θ₀/2), and re-expand the first two terms to get 1 + θ₀²/16 + ….
(b) With `pendulum_period_ratio` and `pendulum_period_ratio_exact`, tabulate T/T₀ and the resulting *inertia* error 1 − (T₀/T)² for amplitudes 10°, 15°, 23°, 30°, 45°. At what amplitude does the inertia error reach 1 %, and at what amplitude 2 %?
(c) Your test rig can hold the amplitude to θ₀ = 10 + s/2 degrees. Compute the inertia error if T₀ is assumed, and the amplitude correction factor (T₀/T)² you would apply instead. How many series terms does the correction need for 10⁻⁵?
(d) The timing uncertainty of the rig is 0.2 % of the period. Compare it with the amplitude correction; decide whether to correct, to reduce the amplitude, or to ignore, and say what would change your decision.

Answer form: one derivation; a table and two amplitudes; two numbers and a count; a decision with its reasoning.
