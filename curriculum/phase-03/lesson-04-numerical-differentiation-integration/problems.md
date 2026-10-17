# P03.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are bounds and which are observed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: deceleration from a difference of the velocity profile)

Use the lesson's capsule (`EntryParams()` defaults: V_e = 7800 m/s, γ = 3°, β = 300 kg/m², H = 7.2 km) at altitude h = (40 + s/2) km.

(a) Compute the deceleration a = ρV²/(2β) from the closed form and from `deceleration_from_velocity` with dh = 100 m, 50 m and 25 m. Report the three errors and the observed order.
(b) The central-difference error is dh²V‴/6 times V sin γ. Estimate V‴ at h by a fourth-order stencil of your own, or by differentiating (4) symbolically, and predict the dh = 100 m error. Compare.
(c) Predict the optimal step h* = (3εF₀/F₃)^{1/3} with F₀ = V(h), F₃ = |V‴|, and sweep dh = 10, 1, 0.1, 0.01 m to find where the error stops falling.

Answer form: a table of dh against error, the order, the predicted and observed dh = 100 m errors, h* and the observed turning point.

## 2. (Derivation: the error laws and the Allen–Eggers closed forms)

(a) From Taylor expansions derive the leading error terms of the forward, central and second-difference stencils.
(b) Derive the one-panel trapezoid error −h³f″(ξ)/12 by integrating the linear interpolation error from `P03.03`, and the composite law (2).
(c) Show that Simpson's rule is exact for every cubic on a symmetric pair of panels, and deduce from the composite error of a quartic that the constant in (3) is 1/180.
(d) Starting from m dV/dt = −½ρV²C_DA and dh/dt = −V sin γ, derive the velocity profile (4), the deceleration a = ρV²/2β, the altitude h_a = H ln 2B of its peak and the value a_max = V_e² sin γ/(2eH), and the peak-heating altitude h_q = H ln 6B for q̇ ∝ √ρ V³.
(e) Derive the heat-load closed form (5) by the substitution u = e^{−h/H}.

Answer form: five derivations; (d) and (e) should end with the lesson's 45.6 km, 81.3 m/s², 53.5 km and 119.87 MJ/m² for the default parameters.

## 3. (Unfamiliar scenario: total impulse of a thrust curve)

A solid motor's thrust (SYNTHETIC: a formula standing in for a test-stand record) rises linearly from 0 to 2000 N over the first 0.2 s and then decays as F(t) = 2000 exp(−(t − 0.2)/τ) N with τ = 1 + s/15 seconds. The record runs from t = 0 to 6 s. The total impulse is the integral of F.

(a) Compute the exact total impulse from the formula.
(b) Sample the record at Δt = 0.1 s (60 panels, so the kink at 0.2 s falls on a sample) and integrate the samples with `trapezoid` and `simpson`. Report both errors. Then sample at Δt = 0.15 s (40 panels, kink between samples) and repeat. Explain the two Simpson results.
(c) Halve the step from 0.1 s to 0.05 s. By what factor should each rule's error fall if the kink is on a node and the integrand is piecewise smooth? What do you observe? Predict the trapezoid error at Δt = 0.1 s quantitatively from the Euler–Maclaurin end terms of the two smooth pieces.

Answer form: the exact impulse in N·s, four errors for (b), two for (c), and a one-paragraph explanation of the Simpson pair.

## 4. (Uses the code: a heavier capsule)

Repeat the lesson's heat-load study for β = (200 + 10 s) kg/m², other parameters default.

(a) Compute B, the peak-heating altitude, the peak deceleration and its altitude, and the analytic heat load from 120 km. Check the scaling Q ∝ 1/√B and h_q = H ln 6B against the lesson's β = 300 values.
(b) Tabulate the relative error of `heat_load_numeric` for trapezoid, Simpson and Gauss at n = 16, 32, 64, 128, 256, and the observed orders of trapezoid and Simpson over 32–128.
(c) Find the smallest even n at which each rule meets 10⁻⁶ relative error. State which rule you would use when the integrand is a function and which when it is a set of samples.

Answer form: five numbers for (a) with the scaling check, the 3 × 5 table and two orders for (b), three counts and a two-line recommendation for (c).

## 5. (Uses the code: when the laws fail)

(a) The lesson's Simpson column on ∫₀² e^{−x}cos 3x dx reaches 1.2 × 10⁻⁷ at n = 64. Find the smallest n at which Simpson reaches 10⁻¹⁰, and the smallest Gauss point count that does. Compare the cost in evaluations.
(b) Print the Romberg diagonal for the same integral to level 6 and explain why the level-3 entry (three extrapolations over 1, 2, 4 and 8 panels) is *worse* than the level-2 entry (Boole's rule on 4 panels).
(c) Integrate x^a on [0, 1] with a = 0.5 + s/100 by trapezoid, Simpson and Gauss for n = 8, 16, …, 128 and report the observed orders. Explain each in terms of the smoothness of the integrand at x = 0, and say what a general-purpose integrator does about it (hint: `scipy.integrate.quad`).

Answer form: two counts for (a), the six-entry diagonal with an explanation for (b), three orders with explanations for (c).

## 6. (Judgment: a quadrature rule for a trajectory tool)

A 6-DOF simulation delivers the entry trajectory as samples of time, altitude, velocity and density every 0.1 s (about 40 m in altitude near peak heating, 400 m near entry interface). The thermal group asks for the stagnation heat load with a relative error below 10⁻⁶ and an error estimate attached to every number.

(a) Which rule from this lesson can run directly on the samples, and which cannot? Why?
(b) Using the lesson's resolution study (or your Problem 4 table), what sample spacing in altitude meets 10⁻⁶ for Simpson, and does the delivered spacing meet it?
(c) Propose the error estimate the tool should report, built from the samples alone (hint: Richardson with every other sample), and state what it would read for the lesson's default capsule at 1 km spacing.
(d) The analyst notices the model's own error is about 10 % (part 2.8). Is the 10⁻⁶ requirement sensible? Write three sentences for the thermal group.

Answer form: four short paragraphs; (c) must contain a number.

## 7. (Proof: Gauss–Legendre is exact to degree 2n − 1 and has positive weights)

Let x₁, …, xₙ be the roots of the Legendre polynomial Pₙ on [−1, 1] and let wᵢ = ∫₋₁¹ ℓᵢ(x) dx where ℓᵢ is the Lagrange basis polynomial for the nodes.

(a) Show the rule Σwᵢf(xᵢ) is exact for every polynomial of degree ≤ n − 1.
(b) For a polynomial p of degree ≤ 2n − 1 write p = qPₙ + r with deg q, deg r ≤ n − 1 and show the rule integrates p exactly. Which property of Pₙ is used, and where?
(c) Show that no n-point rule can be exact for degree 2n (hint: consider Π(x − xᵢ)²).
(d) Show that every weight is positive (hint: apply the rule to ℓᵢ²).
(e) Check (a)–(d) numerically with `gauss_legendre_nodes(5)`: integrate x⁹ and x¹⁰ on [−1, 1], and show the weights sum to 2.

Answer form: four short proofs and the two numerical integrals with their exact values.
