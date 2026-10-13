# P01.03  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which results are exact and which are numerical, and for every numerical integral state the rule, the number of panels and an error estimate.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your own section)

A symmetric NACA four-digit section has thickness ratio t = 0.08 + s/400 and chord c = 1 + s/20 m (open trailing edge). By hand, with the closed forms of lesson §2.4 and §2.5: (a) the area A in m² and as a multiple of t c²; (b) the chordwise centroid x̄ in m and as a fraction of c. Then with the code: (c) `section_area(t, c, n=1000)` and its error against (a); (d) predict that error from the law |error| ≈ 0.022 (t/0.12) c² n^(−1.5) and compare; (e) `section_area_substituted(t, c, n=1000)` and its error. Answer form: five numbers with units, four significant figures, and one sentence on (d).

## 2. (Derivation: the closed forms and the error law)

(a) Starting from (1), derive A = 10 t c² Σ aₖ/(eₖ + 1) and x̄ = c [Σ aₖ/(eₖ + 2)]/[Σ aₖ/(eₖ + 1)], naming the Fundamental Theorem (part 2) and the power rule for integrals where you use them, and showing the substitution ξ = x/c explicitly.
(b) Prove that the midpoint rule is exact for every linear function f(x) = αx + β, and derive the error law ∫ₐᵇ f − M_n ≈ (h²/24)[f′(b) − f′(a)] from a Taylor expansion about each panel midpoint, stating what you assume about f.
(c) Show that the substitution x = c u² turns 2∫₀ᶜ y_t dx into ∫₀¹ g(u) du with g a polynomial, and write out g(u) for t = 0.12, c = 1. Evaluate g′(1) − g′(0) and predict the midpoint error at n = 1280.

Answer form: derivation with each rule named; one polynomial; one predicted number.

## 3. (Unfamiliar scenario: total impulse of a model-rocket motor)

A small solid motor's thrust curve is piecewise linear through the points (t, F) = (0, 0), (0.1 s, F_peak), (0.3 s, 18 N), (1.2 s, 12 N), (1.5 s, 0), with F_peak = 20 + s/3 N.

(a) The **total impulse** is I = ∫₀^1.5 F dt (the net change of the rocket's momentum from thrust alone). Compute it exactly by adding trapezoids, and the average thrust I/1.5.
(b) Compute I with `riemann_sum` using `np.interp` for F, with n = 15, 150 and 1500 and all three rules. Explain why all nine numbers agree, with reference to where the panel edges fall and to the exactness of the midpoint rule.
(c) A rocket of mass 0.25 kg leaves the pad; neglecting drag and mass change, its speed at burnout is Δv = I/m − g × 1.5 s. Compute it. Then let the mass fall linearly from 0.25 to 0.22 kg during the burn, compute Δv = ∫(F/m(t)) dt − g × 1.5 numerically, and compare with the constant-mass value and with using the average mass.
(d) Which of the three Δv estimates would you quote, and what would make all three wrong by more than 10 %?

Answer form: two numbers with units, one explanation, three speeds with units, two sentences.

## 4. (Numerical, uses your code: convergence and the Fundamental Theorem)

(a) Reproduce the convergence table of lesson part 4 for n = 10 … 1280 and fit the orders on n ≥ 40. Then fit them on n ≥ 320 only: which rules' orders change, and in which direction?
(b) Use the signed errors (not absolute values) of the left and right sums to show that the first-order terms have opposite signs and the n^(−1.5) terms the same sign, and compute (L_n + R_n)/2 − A for each n. What order does this average converge with, and why?
(c) Check the Fundamental Theorem part 1 at x = 0.05 c, 0.5 c and 0.95 c with `cumulative_integral` (n = 2000) and a central difference with h = 10⁻³, 10⁻⁴. Where is the error largest, and why there?
(d) Plot everything on one figure and save it as `riemann_convergence_full.png`.

Answer form: a table, eight fitted orders, one signed-error table with a conclusion, three FTC residuals, one figure.

## 5. (Centroid of a cambered section, and the vertical-thickness assumption)

(a) With `section_centroid`, find ȳ for NACA 2412, 4412 and 6412 (m = 0.02, 0.04, 0.06; p = 0.4; t = 0.12). Show ȳ is linear in m and give ȳ/m.
(b) Derive why ȳ/m is independent of m under the vertical-thickness assumption, and why it is less than 1 (the maximum camber) even though the camber line's own mean height is used.
(c) Estimate, by an order-of-magnitude argument, the error in A caused by adding the thickness vertically instead of perpendicular to the camber line: the surface point moves horizontally by y_t sin θ with θ = arctan(dy_c/dx). Bound the maximum of θ for NACA 2412 and state whether the area error is first or second order in θ.
(d) ★ Implement the perpendicular construction (upper surface x_u = x − y_t sin θ, y_u = y_c + y_t cos θ; lower with opposite signs) and compute the enclosed area with the shoelace formula on 2000 points. Report the difference from the vertical-thickness area for NACA 2412 and 6412.

Answer form: three numbers, one ratio, a two-sentence argument, a bound in degrees, an order; ★ two percentages.

## 6. (Judgment: when a Riemann sum misleads)

For each integrand on [0, 1], say what goes wrong with a midpoint sum, what order you expect, what sign in the convergence plot reveals the problem, and what you would do instead:

(a) f(x) = x^(−1/2) (integrable, value 2).
(b) f(x) = 1/x (not integrable).
(c) f(x) = a step function: 0 for x < 0.37, 1 for x ≥ 0.37 (a deflected flap's extra drag switching on).
(d) f(x) = sin(200x) (an oscillatory pressure signal).
(e) f given only at 11 measured points with ±1 % noise.

Answer form: five short paragraphs, each with an order (or "none") and an observable.

## 7. (Hard: second moment of area of the section)

For bending about the chord line, a solid symmetric section has second moment of area I = ∫₀ᶜ (2/3) y_t³ dx (each vertical strip of height 2y_t contributes (2y_t)³/12).

(a) Expand y_t³ symbolically (SymPy) and integrate exactly to show I = κ t³ c⁴ and find κ to five figures.
(b) Compute I for the NACA 0012 (c = 1) with midpoint sums for n = 40 … 1280 and fit the order. It is not 1.5 and not 2; explain the value you get by considering the behaviour of y_t³ near the leading edge *and* the size of f′(c) − f′(0) in the error law.
(c) Compute the radius of gyration √(I/A) as a fraction of c, and compare I with that of a rectangle c × tc (I = c (tc)³/12). Why is the section's I less than half the rectangle's when its area is 68 % of the rectangle's?
(d) A foam wing of chord 0.25 m is loaded in bending. By what factor does I change if the section is changed from NACA 0012 to NACA 0015, and from chord 0.25 m to 0.30 m? Which change is cheaper in mass?

Answer form: one constant, one order with a two-part explanation, two numbers and one sentence, two factors and a judgment.
