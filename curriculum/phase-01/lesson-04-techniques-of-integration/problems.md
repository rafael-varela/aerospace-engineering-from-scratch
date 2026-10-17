# P01.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, name the integration technique used for every closed form, and for every numerical integral state the rule, the panel count and an error estimate or bound.

Several problems use your personal seed **s**: the day of the month you were born (1–31). "The lesson's wing" means b = 11 m, c_r = 1.6 m, τ = 0.6, σ = 12 kg/m².

## 1. (Warm-up: your own wing)

Your wing has span b = 9 + s/5 m, root chord c_r = 1.2 + s/50 m, taper ratio 0.5 and areal density σ = 10 + s/10 kg/m², with mass per unit span λ = σ c(y). By hand: (a) the mass of the wing pair; (b) the spanwise centroid of one half-wing; (c) the roll moment of inertia I_xx from (1). Check (c) with `roll_inertia_numeric(..., n=2, rule="simpson")` and say why two panels suffice. (d) Compute I_xx for a rectangular wing of the same mass and span and the radius of gyration √(I/m) of your wing. Answer form: five numbers with units, four significant figures.

## 2. (Derivation: the four techniques)

(a) Derive (3) by integrating ∫₀ˢ y² e^(−y/ℓ) dy by parts twice, showing the boundary terms at each step.
(b) Derive (2) by the substitution y = s sin θ, including the double-angle reduction of sin²θ cos²θ, and state why the substitution is valid on [0, s] (which branch of the square root, which range of θ).
(c) Derive (4) by partial fractions, showing the decomposition of 1/(T − k v²) and the two logarithms; then show that for k v² ≪ T the result reduces to t ≈ m v/T + (m k/(3T²)) v³, and interpret the correction.
(d) Show that the Simpson error law contains f⁽⁴⁾ by writing Simpson as (2M_n + T_n)/3 and using the midpoint and trapezoid error laws; conclude that Simpson is exact for cubics.

Answer form: four derivations with each technique named where used.

## 3. (Unfamiliar scenario: the Laplace transform as an improper integral)

Phase 02 defines the Laplace transform F(s) = ∫₀^∞ e^(−st) f(t) dt, which needs this lesson's improper integrals and integration by parts.

(a) Compute F(s) by hand for f = 1, f = t (parts once), f = e^(−at), and f = sin ωt (parts twice, or Euler's formula), stating for each the condition on s that makes the improper integral converge.
(b) With s = 2, a = 1, ω = 3, use `improper_integral` with cutoffs 2, 5, 10, 20 to build the four cutoff sequences and compare the last term with your closed forms.
(c) For f = e^(+t) and s = 2, does the integral converge? For s = 0.5? Show the cutoff sequence in both cases and explain the pattern.
(d) In one sentence, why does the convergence condition on s matter for an engineer using the transform on a system that grows exponentially?

Answer form: four closed forms with conditions, four comparisons, two sequences with an explanation, one sentence.

## 4. (Numerical, uses your code: rules, bounds and the tip)

(a) Reproduce the elliptic-planform convergence table of lesson part 4 and fit the order of each raw rule with a log-log fit on n ≥ 16. Explain why Simpson's extra accuracy does not survive.
(b) Replace the elliptic integrand by the θ-form and run the *midpoint* and *Simpson* rules on it (not only the trapezoid). Which of the three is exact at n = 4, and why does the trapezoid rule do better than Simpson here?
(c) For the exponential wing, tabulate error and bound for the trapezoid rule at n = 8, 16, 32, 64, 128 and for Simpson at n = 8, 16, 32. Compute the ratio bound/error for each and explain why it is nearly constant.
(d) Find, from the bounds alone, the smallest n for which the trapezoid rule and Simpson's rule each guarantee four significant figures in I_xx of the exponential wing. Then verify by running them.

Answer form: a table and three orders; one sentence; a table of ratios with an explanation; two panel counts with verification.

## 5. (Trigonometric substitution in the load: root bending moment)

A wing of semispan s = 5.5 m carries an elliptic lift distribution L′(y) = L₀ √(1 − (y/s)²) N/m, and each half-wing lifts half the aircraft's weight, W = 10 kN.

(a) Find L₀ from ∫₀ˢ L′ dy = W/2 (trigonometric substitution).
(b) Derive the root bending moment M₀ = ∫₀ˢ y L′(y) dy in closed form (the substitution gives ∫ sin θ cos²θ dθ) and show it equals 2Ws/(3π). Evaluate it.
(c) Compare with the root bending moment for a lift distribution proportional to the lesson's tapered chord (M₀ = (W/2) ȳ with ȳ from P01.04 §2.2) and for a uniform distribution. Which is kindest to the spar, and why?
(d) Compute the shear force and bending moment at y = 2 m for the elliptic case numerically (Simpson), stating n and an error estimate; the integrals are ∫₂ˢ L′ dy and ∫₂ˢ (y − 2) L′ dy.

Answer form: L₀ with units, one derivation and a number, three moments and a sentence, two numbers with error estimates.

## 6. (Judgment: where the beam model and the rules stop)

For each statement, name the assumption that fails, estimate the size of the effect with one number, and give a measurable sign:

(a) The wing carries 60 kg of fuel in tanks centred at y = 3.5 m on each side.
(b) The chord is constant over the inner 1.5 m (the fuselage), and the fuselage itself contributes to I_xx.
(c) The tabulated planform has 12 chord values at 0.5 m spacing and the tip chord is listed as 0.
(d) The take-off is on grass, with rolling friction coefficient 0.08, and the propeller's thrust falls from 2500 N static to 1700 N at 30 m/s.

Answer form: four short paragraphs, each with a number and an observable.

## 7. (Hard, uses your code: a heavier-at-the-root mass law)

A better mass model for a wing of constant thickness-to-chord ratio is λ = σ′ c(y)², because the spar depth and the skin area per unit span both scale with chord. Take σ′ = σ/c_r so the root mass per unit span is unchanged (19.2 kg/m for the lesson's wing).

(a) Derive I_xx = 2σ′ ∫₀ˢ y² c² dy in closed form for the linear taper (the integrand is a quartic) and for the elliptic planform (c² is a polynomial: no substitution needed).
(b) Evaluate both for the lesson's root chord and span, and the wing masses and radii of gyration. Compare with the λ ∝ c results (1490.7 kg m², 169.0 kg; 1254.4 kg m², 165.9 kg).
(c) Simpson's rule with n = 2 is no longer exact for the linear taper. By how much is it off, and why? What n makes it exact to 10⁻⁶?
(d) Judgment: which mass law would you use for (i) a foam-core model aircraft wing, (ii) a metal light-aircraft wing with a single spar, (iii) a sailplane wing with water ballast in the leading edge? Give a one-line reason each.

Answer form: two closed forms, six numbers with units, one error with explanation and one n, three choices with reasons.
