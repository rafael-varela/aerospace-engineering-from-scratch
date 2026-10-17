# P03.04  Numerical Differentiation & Integration

> A derivative is a difference that has forgotten its step; an integral is a sum that has forgotten its panels. Numerical methods put the step and the panels back and tell you exactly what remembering costs.

| | |
|---|---|
| **Learning outcomes** | 1. Derive from Taylor expansions the forward, central and second-difference stencils with their leading error terms, the Richardson extrapolation that cancels them, and the optimal step where truncation meets rounding. 2. Derive the trapezoid and Simpson error laws from the interpolation error, show why Simpson is exact for cubics and why n-point Gauss–Legendre is exact to degree 2n − 1, and explain when the trapezoid rule is spectrally accurate. 3. Implement and test, from scratch, the stencils, Richardson, trapezoid, Simpson, Romberg, Gauss–Legendre nodes and weights by Newton on the Legendre recurrence, and adaptive Simpson, against `scipy.integrate` and `numpy.polynomial.legendre.leggauss`. 4. Compute the deceleration and total heat load of a ballistic re-entry from the Allen–Eggers profile, check both against closed forms, measure the observed orders, and decide which rule and how many panels a trajectory tool needs for a 10⁻⁶ heat-load requirement. |
| **Prerequisites** | `P03.01` (truncation against rounding, observed order), `P03.03` (interpolation error formula: every quadrature rule here integrates an interpolant), `P01.05` (Taylor series), `P02.01` (the Allen–Eggers velocity profile is a first-order ODE solved by separation; here it is taken as given). Skill check: you can expand f(x ± h) to fourth order and integrate a parabola through three points. |
| **Engineering question** | A ballistic capsule (ballistic coefficient 300 kg/m², nose radius 1 m) enters at 7.8 km/s and 3° below the horizontal. What peak deceleration and peak heating rate does it meet, at what altitudes, and what total heat load per square metre must its heat shield absorb from 120 km to the ground? The deceleration is a derivative of the velocity profile and the heat load an integral of the heating rate; both have closed forms in the Allen–Eggers model, so the question becomes: which numerical rule, at what resolution, reproduces them to 10⁻⁶, and how do you know the resolution is enough when there is no closed form? |
| **Assumptions and limits** | Allen–Eggers ballistic entry: no gravity, no lift, constant flight-path angle γ and ballistic coefficient β, exponential atmosphere ρ = ρ₀e^{−h/H} with H = 7.2 km; the model is good from entry to peak heating and poor below about 30 km, where it brings the vehicle to rest instead of to a terminal velocity (gravity is missing). Convective stagnation heating in Sutton–Graves form q̇ = k√(ρ/Rₙ)V³ with k = 1.74 × 10⁻⁴ kg^½/m, which omits radiative heating (small below 9 km/s). Numerically: all integrands and functions are smooth enough for the stated error laws (C² for trapezoid, C⁴ for Simpson, analytic for the spectral claims); a kink or an endpoint singularity reduces every order, as Problems 3 and 5 show. |
| **Analytical work** | Stencil derivations with leading constants hf″/2, h²f‴/6, h²f⁗/12; the Richardson formula and its order gain; the optimal step h* = (3εf/f‴)^{1/3} for the central difference; the trapezoid error −(b − a)h²f″/12 and the Euler–Maclaurin view (the error is −h²[f′(b) − f′(a)]/12 + …, which vanishes for periodic integrands); Simpson's exactness for cubics by symmetry and its −(b − a)h⁴f⁗/180 law; Gauss–Legendre exactness to degree 2n − 1 by division by Pₙ; the Allen–Eggers velocity V = V_e exp(−B e^{−h/H}), the deceleration ρV²/2β, its maximum at h = H ln 2B with value V_e² sin γ/(2eH), the peak heating altitude H ln 6B, and the heat load in closed form through erf. |
| **Computational deliverable** | `starter/main.py`: `forward_diff`, `central_diff`, `second_diff`, `richardson`, `optimal_step_central`, `trapezoid`, `simpson`, `romberg`, `gauss_legendre_nodes`, `gauss_legendre`, `adaptive_simpson`, `observed_order`, `entry_B`, `density_exp`, `entry_velocity`, `deceleration_analytic`, `deceleration_from_velocity`, `peak_deceleration`, `heating_rate`, `peak_heating_altitude`, `heat_load_integrand`, `heat_load_numeric`, `heat_load_analytic`. Tests prove: the stencils have orders 1, 2, 2 within 0.02 and their leading constants within 1 %; forward is exact for lines, central and second for quadratics; Richardson lifts central to order 4 and forward to 2 and is exact on a cubic; the central optimum lies within a decade of h* with error below 10⁻¹⁰; trapezoid is exact for lines and Simpson for cubics but not quartics; orders 2.0 ± 0.1 and 4.0 ± 0.2 on e^{−x}cos 3x with the Euler–Maclaurin constant within 2 %; both equal `scipy.integrate.trapezoid`/`simpson` on the same grid; the trapezoid rule reaches 10⁻¹⁴ with 16 panels on a periodic integrand; the Romberg table's columns are the trapezoid and Simpson rules and its diagonal converges monotonically to 10⁻¹⁰; Gauss nodes and weights equal `leggauss` to 10⁻¹⁴ for n = 1 to 16, are symmetric, positive, sum to 2, integrate every monomial to degree 2n − 1 exactly and fail at degree 2n; Gauss gains more than five digits from 4 to 8 points; adaptive Simpson meets 10⁻⁶ and 10⁻¹⁰, costs more at the tighter tolerance, handles a near-singular peak to 10⁻⁶ in under 5000 evaluations, and raises at `max_depth`; the entry profile has the right limits, its analytic peaks are the maxima of the sampled curves, the numerical deceleration converges at order 2 to the analytic one, and the numerical heat load matches the erf closed form and `quad` with trapezoid order 2.0 ± 0.15, Simpson at least 3.5, Gauss to 10⁻¹⁰ at 64 points. |
| **Practical activity** | **computational experiment**: a resolution study of the stencils and of trapezoid, Simpson, Romberg and Gauss on the heat-load integrand, locating the pre-asymptotic regime set by the 7.2 km width of the heating pulse, the asymptotic orders, and the rounding floors; a cost study (function evaluations against accuracy) across the five rules. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `reentry_heating.png`, the stencil and quadrature tables with observed orders, the work–precision table, the hand derivations of the stencil constants, the trapezoid and Simpson error laws and the Allen–Eggers results, and a half-page note recommending a quadrature rule and panel count for the heat load of a sampled trajectory with an error estimate. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the central-difference constant, the trapezoid law and Gauss exactness are derived; all public tests pass; the resolution study identifies the pre-asymptotic regime and explains it by the pulse width; the note states a rule, a count and the evidence. |
| **References** | `@ocw-18-330` "Lecture notes 3–4: Integrals as Sums and Derivatives as Differences" (required-free): stencils, the truncation–roundoff trade-off, Newton–Cotes rules. `@ocw-18-335j` Week 12 "Lecture 33: Numerical Integration and the Convergence of the Trapezoidal Rule" (required-free): why the trapezoid rule is spectral for periodic and analytic integrands; Week 13 "Lecture 36: Integration with Weight Functions, and Gaussian Quadrature" (required-free): orthogonal polynomials and Gauss rules; Week 12 "Lecture 34: Clenshaw-Curtis Quadrature" and Week 13 "Lecture 37: Adaptive and Multidimensional Quadrature" (supplementary-free). `@scipy-userguide` "Integration (scipy.integrate)" > "General integration (quad)", "Gaussian quadrature", "Integrating using Samples" (required-free). |
| **Extensions** | ★ Implement Clenshaw–Curtis quadrature (Chebyshev nodes, weights by FFT) and compare its convergence with Gauss on the heat-load integrand. ★ Implement the complex-step derivative Im f(x + ih)/h and show it has no rounding floor. ★ Replace the exponential atmosphere by the ISA interpolant of `P03.03` and integrate the heat load numerically where no closed form exists, with a Richardson error estimate. |
| **Time** | unmeasured |

## 1. The engineering problem

A capsule returning from orbit enters the atmosphere at V_e = 7.8 km/s on a shallow ballistic path, γ = 3° below the horizontal. Three numbers size its design: the **peak deceleration**, which sets the structure and what the crew can take; the **peak heating rate**, which sets the heat-shield material; and the **total heat load**, the time integral of the heating rate, which sets the heat shield's thickness. In the classical Allen–Eggers model all three follow from one velocity profile V(h), and the deceleration is its derivative while the heat load is an integral along the path. The model is simple enough to have closed forms (part 2), which makes it the perfect test bed: every numerical derivative and integral in this lesson is checked against an exact answer, and the methods that pass are the ones you then trust on the real trajectory, which comes as samples from a six-degree-of-freedom simulation and has no closed form at all.

Numerical differentiation and integration are the two oldest tasks of numerical analysis and they are opposites in character. **Differentiation amplifies error**: a difference quotient divides a small difference of rounded numbers by a small step, so there is an optimal step and a floor of accuracy that no stencil can pass (`P03.01` met this with the forward difference; here the central and second differences, and Richardson extrapolation, move the floor). **Integration smooths error**: a sum of many rounded values has relative error of order ε, and the only question is how fast the truncation error falls with the number of panels, which ranges from h² for the trapezoid rule through h⁴ for Simpson to faster than any power for Gauss quadrature on smooth integrands. Every rule here is the exact integral of an interpolant from `P03.03`, and its error law follows from the interpolation error formula.

Success criteria, written first: each stencil and rule shows its theoretical order to within 0.1 on a smooth test function and its leading error constant to within a few per cent; the Gauss nodes and weights match `leggauss` to 10⁻¹⁴; the peak deceleration and the heat load match their closed forms to 10⁻⁶ with a stated number of evaluations; and the lesson ends with a resolution rule for a trajectory tool.

## 2. Understand and derive

### 2.1 Stencils from Taylor expansions

Write f(x ± h) = f ± hf′ + h²f″/2 ± h³f‴/6 + h⁴f⁗/24 ± …. Then

- **forward**: (f(x + h) − f(x))/h = f′ + hf″/2 + O(h²): first order;
- **central**: (f(x + h) − f(x − h))/2h = f′ + h²f‴/6 + O(h⁴): second order, the odd terms cancel;
- **second difference**: (f(x + h) − 2f(x) + f(x − h))/h² = f″ + h²f⁗/12 + O(h⁴): second order.

`main()` measures them on sin at x = 1: at h = 10⁻³ the errors are 4.21 × 10⁻⁴, 9.01 × 10⁻⁸ and 7.01 × 10⁻⁸, and the predictions hf″/2 = 4.21 × 10⁻⁴, h²f‴/6 = 9.00 × 10⁻⁸, h²f⁗/12 = 7.01 × 10⁻⁸ agree to three figures; the observed orders over h = 10⁻¹ to 10⁻³ are 1.00, 2.00, 2.00. The central difference is exact for quadratics and the second difference too, which the tests check with exact arithmetic.

**Rounding** adds ε|f|/h to a first-derivative stencil and about 4ε|f|/h² to the second difference. For the central difference the total h²F₃/6 + εF₀/h is minimised at **h* = (3εF₀/F₃)^{1/3}** ≈ 8.7 × 10⁻⁶ for F₀ = F₃ = 1, with best error ≈ 10⁻¹¹: two thirds of the digits, against half for the forward difference (`P03.01`). Observed: the best central error is 1.1 × 10⁻¹¹ at h = 10⁻⁵, and the error at h = 10⁻⁸ has climbed back to 2.6 × 10⁻⁹. The second difference is worse off: its optimum near h ≈ ε^{1/4} ≈ 10⁻⁴ gives only 3 × 10⁻⁹ and at h = 10⁻⁸ it returns 0.27 for a quantity of size 0.84: garbage. A second derivative by differences is a last resort.

### 2.2 Richardson extrapolation

If an estimate satisfies D(h) = exact + Chᵖ + O(hᵖ⁺¹), then D(h/2) = exact + Chᵖ/2ᵖ + …, and the combination

**(2ᵖD(h/2) − D(h))/(2ᵖ − 1) = exact + O(hᵖ⁺¹)**   (1)

cancels the leading term using nothing but two evaluations of the same formula. When the expansion contains only even powers (central differences, trapezoid), the gain is two orders. Applied to the central difference with p = 2: errors 1.13 × 10⁻⁷ at h = 0.1 and 1.13 × 10⁻¹¹ at h = 0.01, order 4.00; by h = 10⁻³ it has hit the rounding floor at 3.5 × 10⁻¹⁴. Applied to the trapezoid rule it gives Simpson; applied again, Boole's rule; applied repeatedly to a sequence of halvings it is **Romberg integration** (2.4). Richardson is also the standard **error estimator**: D(h/2) − D(h) is about (1 − 2⁻ᵖ) times the error of D(h), which is how adaptive quadrature and adaptive ODE solvers (`P03.06`) decide when to refine.

### 2.3 Trapezoid and Simpson: integrating an interpolant

Replace f on [xᵢ, xᵢ₊₁] by the chord and integrate: the **trapezoid rule** h(fᵢ + fᵢ₊₁)/2. Its error on one panel is the integral of the interpolation error ½f″(ξ)(x − xᵢ)(x − xᵢ₊₁), which is −h³f″/12; over n panels of total length b − a,

**E_trap = −(b − a)h²f″(ξ)/12**.   (2)

Second order. The Euler–Maclaurin formula sharpens this: E_trap = −h²[f′(b) − f′(a)]/12 + h⁴[f‴(b) − f‴(a)]/720 − …, every term a difference of odd derivatives at the ends. Two consequences. First, the tests check the constant: for ∫₀² e^{−x}cos 3x with 128 panels the error is 2.00 × 10⁻⁵ against h²[f′(2) − f′(0)]/12 = 2.00 × 10⁻⁵. Second, if all the end differences vanish, because the integrand is **periodic over the interval** or decays to zero with all its derivatives at both ends, the trapezoid rule has no power-law error at all and converges exponentially (18.335J Lecture 33): ∫₀²π e^{cos x}dx is right to 10⁻⁵ with 8 panels and 10⁻¹⁴ with 16, a factor 10⁹ for one doubling instead of 4. This is why the FFT-based methods of `P02.10` work, and it is why, on the heat-load integrand below (a pulse that nearly vanishes at both ends), the trapezoid rule beats Simpson at coarse resolution.

**Simpson's rule** integrates the parabola through three points: (h/3)(f₀ + 4f₁ + f₂) per pair of panels. Its error should be O(h⁴) from the quadratic interpolant, but it is better: by symmetry about the midpoint the cubic term integrates to zero, so Simpson is **exact for cubics** and its error is

**E_Simp = −(b − a)h⁴f⁗(ξ)/180**,   (3)

fourth order. On e^{−x}cos 3x: 1.15 × 10⁻², 5.5 × 10⁻⁴, 3.2 × 10⁻⁵, 2.0 × 10⁻⁶, 1.2 × 10⁻⁷ for 4 to 64 panels, order 4.11 (trapezoid 2.03 on the same runs). Both rules use n + 1 equally spaced samples, so they are the rules for **data** (a sampled trajectory, a thrust curve), and `scipy.integrate.trapezoid` and `simpson` take exactly those samples.

### 2.4 Romberg: Richardson on the trapezoid rule

Compute the trapezoid rule with 1, 2, 4, …, 2ⁱ panels (each halving reuses the previous samples) and extrapolate across the table: R[i, j] = (4ʲR[i, j − 1] − R[i − 1, j − 1])/(4ʲ − 1). Column 0 is the trapezoid rule, column 1 is Simpson (the tests confirm this to 10⁻¹³), column 2 is Boole, and the diagonal converges as fast as the integrand's smoothness allows: for ∫₀¹ dx/(1 + x²) = π/4 the diagonal reaches 10⁻¹¹ at level 6 (33 evaluations) while the trapezoid column is still at 10⁻⁵. Romberg is Gauss's accuracy with Simpson's samples, and the difference of the last two diagonal entries is its error estimate.

### 2.5 Gauss–Legendre: choosing the nodes

An interpolatory rule on n fixed nodes is exact for polynomials of degree n − 1. If the nodes are free, n nodes and n weights are 2n parameters and one can hope for degree 2n − 1, and the **Gauss–Legendre** rule achieves it: the nodes are the roots of the Legendre polynomial Pₙ and the weights wᵢ = 2/((1 − xᵢ²)Pₙ′(xᵢ)²). The proof (Problem 7): divide a polynomial of degree ≤ 2n − 1 as p = qPₙ + r with deg q, deg r ≤ n − 1; ∫qPₙ = 0 by orthogonality; the rule integrates r exactly and gives zero on qPₙ because Pₙ vanishes at every node. The weights are positive (wᵢ = ∫ℓᵢ² > 0) so there is no cancellation, and for an integrand analytic in a neighbourhood of the interval the error falls like ρ⁻²ⁿ for some ρ > 1: **spectral** convergence. On e^{−x}cos 3x: 5.7 × 10⁻⁴ with 4 points, 1.5 × 10⁻¹¹ with 8, rounding floor at 16. Five digits per four extra points. The cost: the nodes are irrational, so Gauss needs the integrand as a *function*, not as samples. `gauss_legendre_nodes` finds the roots by Newton's method (`P03.02`) on the three-term recurrence (k + 1)Pₖ₊₁ = (2k + 1)xPₖ − kPₖ₋₁ from the asymptotic guess cos(π(i + ¾)/(n + ½)); the result equals `numpy.polynomial.legendre.leggauss` to 10⁻¹⁴ for every n tested.

### 2.6 Adaptive quadrature

When the integrand has a sharp feature somewhere (a peak, a near-singularity), a uniform grid wastes panels where f is flat. **Adaptive Simpson** compares Simpson on an interval with Simpson on its two halves; by (3) the difference is 15 times the error of the finer estimate, so the interval is accepted when |S₂ − S₁| ≤ 15 tol and split otherwise, with the tolerance shared between halves. The accepted value includes the Richardson correction (S₂ − S₁)/15. For the smooth e^{−x}cos 3x it is wasteful (645 evaluations to reach 2.8 × 10⁻¹⁵ against Romberg's 33 for 10⁻⁹: the tolerance splitting is conservative); for 1/(10⁻³ + (x − 0.3)²), a peak of width 0.03 on [0, 1], it matches `quad` to 10⁻⁶ relative with fewer than 5000 evaluations by concentrating them near 0.3. `scipy.integrate.quad` is the production version of this idea (QUADPACK: adaptive Gauss–Kronrod with extrapolation), and it is the oracle for every heat load in the tests.

### 2.7 The ballistic entry

Along a straight path at angle γ below the horizontal, with drag as the only force, m dV/dt = −½ρV²C_DA and dh/dt = −V sin γ. Dividing, dV/dh = ρV/(2β sin γ) with β = m/(C_DA), and with ρ = ρ₀e^{−h/H} this separates:

**V(h) = V_e exp(−B e^{−h/H}), B = ρ₀H/(2β sin γ)**.   (4)

For the example B = 281: the exponent is negligible above 80 km, V has fallen to e⁻¹V_e where e^{−h/H} = 1/B (h = 40.6 km), and the model drives V to zero at sea level because gravity, which would provide a terminal velocity, is absent. Below about 30 km (V = 100 m/s) the model is wrong and the heat load there is negligible anyway.

**Deceleration** a = −dV/dt = (dV/dh)V sin γ = **ρV²/(2β)**: drag per unit mass, as it must be. It peaks where d(ρV²)/dh = 0, i.e. e^{−h/H} = 1/(2B): **h_a = H ln 2B = 45.6 km**, with **a_max = V_e² sin γ/(2eH) = 81.3 m/s² = 8.3 g**, independent of β (a heavier capsule decelerates lower, not harder). The test computes a from a central difference of (4) and gets 81.3 to 10⁻⁶ relative with dh = 10 m, with observed order 2.0 as dh is halved.

**Heating.** q̇ = k√(ρ/Rₙ)V³ ∝ e^{−h/2H}exp(−3Be^{−h/H}), maximal where e^{−h/H} = 1/(6B): **h_q = H ln 6B = 53.5 km**, above the deceleration peak (V is still 6.6 km/s there), with q̇_max = 135 W/cm² for Rₙ = 1 m. The **heat load** is Q = ∫q̇ dt = ∫q̇ dh/(V sin γ) = (k√(ρ₀/Rₙ)/sin γ)∫V² e^{−h/2H}dh. Substituting u = e^{−h/H}, dh = −H du/u, V² = V_e²e^{−2Bu}:

**Q = C√(π/2B)[erf(√(2B)) − erf(√(2Bu_top))], C = k√(ρ₀/Rₙ)V_e²H/sin γ**,   (5)

since ∫u^{−1/2}e^{−2Bu}du = √(π/2B) erf(√(2Bu)). From 120 km: **Q = 119.87 MJ/m² = 12.0 kJ/cm²**; 90 % of it accrues below 80 km and 61 % below 60 km. The integrand q̇/(V sin γ) in h is a pulse of width about H centred near 50 km, nearly zero at both ends of [0, 120 km]: smooth, but narrow compared with the range, which is what makes the resolution study of part 4 interesting.

### 2.8 Where the model stops being valid

- No gravity: fine above ~35 km (drag ≫ g), wrong below 30 km; for shallow entries (γ < 2°) the omission matters even at peak heating.
- Constant γ: real ballistic entries steepen slightly; lifting entries (Apollo, Soyuz nominal) violate it entirely and halve the g-load.
- Sutton–Graves convective heating: ±15 % against detailed CFD; radiative heating adds significantly above 9 km/s (lunar return), not at 7.8.
- Exponential atmosphere: H = 7.2 km fits 20–80 km to about 20 %; `P03.03`'s ISA interpolant replaces it in extension 3 at the cost of the closed form.
- Numerically, the error laws (2), (3) and the spectral claims need the stated smoothness; a thrust curve with an ignition kink (Problem 3) or an integrand like √x (Problem 5) drops Simpson to second or 1.5th order, and no rule is spectral across a kink.

## 3. Build the computational model

Open `starter/main.py`. `EntryParams` is given; every function is yours.

### Step 1: stencils and Richardson

The three stencils are one line each. `richardson(estimate, h, order, ratio)` evaluates `estimate(h)` and `estimate(h / ratio)` and applies (1) with 2ᵖ replaced by ratioᵖ; validate order ≥ 1 and ratio > 1. `optimal_step_central` is (3εF₀/F₃)^{1/3}.

### Step 2: Newton–Cotes and Romberg

`trapezoid(f, a, b, n)` and `simpson(f, a, b, n)` build `np.linspace(a, b, n + 1)`, evaluate f once on the array, and apply the weights (Simpson: `y[1:-1:2]` get 4, `y[2:-1:2]` get 2; raise for odd n). `romberg(f, a, b, levels)` fills the lower triangle:

```python
for i in range(levels):
    R[i, 0] = trapezoid(f, a, b, 2 ** i)
    for j in range(1, i + 1):
        R[i, j] = (4.0 ** j * R[i, j - 1] - R[i - 1, j - 1]) / (4.0 ** j - 1.0)
```

### Step 3: Gauss–Legendre

`gauss_legendre_nodes(n)` starts from x = cos(π(i + ¾)/(n + ½)), and in each Newton sweep runs the recurrence to get Pₙ(x) and Pₙ₋₁(x), forms Pₙ′ = n(xPₙ − Pₙ₋₁)/(x² − 1), and updates x ← x − Pₙ/Pₙ′ until the step is below 10⁻¹⁵; then w = 2/((1 − x²)Pₙ′²); sort by x. `gauss_legendre(f, a, b, n)` maps the nodes to [a, b] and scales the weights by (b − a)/2.

### Step 4: adaptive Simpson

Recursive: given f at a, m, b and the whole-interval Simpson S₁, evaluate f at the quarter points, form the two half Simpsons, accept when |S₂ − S₁| ≤ 15 tol returning S₂ + (S₂ − S₁)/15, otherwise recurse on each half with tol/2. Count evaluations through a closure; raise `RuntimeError` beyond `max_depth`.

### Step 5: the entry model

`entry_B`, `density_exp`, `entry_velocity`, `deceleration_analytic`, `heating_rate`, `heat_load_integrand`, `peak_deceleration`, `peak_heating_altitude` and `heat_load_analytic` type in 2.7; `deceleration_from_velocity(h, p, dh)` is `central_diff` of `entry_velocity` times V sin γ; `heat_load_numeric(h_top, p, method, n)` dispatches to the three rules on `heat_load_integrand`.

Running `main()` prints the stencil table of 2.1, the quadrature table of 2.3–2.6, and

```
ballistic entry: V_e = 7800 m/s, gamma = 3 deg, beta = 300 kg/m^2, R_n = 1.0 m, B = 280.9
peak heating at h = 53.5 km: qdot = 135.1 W/cm^2, V = 6603 m/s
peak deceleration at h = 45.6 km: 81.3 m/s^2 = 8.3 g; from a central difference of V: 81.3 m/s^2
heat load from 120 km to sea level: analytic 119.8708 MJ/m^2 (erf closed form)
    n     trapezoid       Simpson         Gauss   (relative errors)
    8      2.42e-02      7.38e-02      9.00e-02
   16      2.28e-04      8.37e-03      3.38e-03
   32      3.65e-05      2.73e-05      2.05e-06
   64      9.14e-06      8.78e-09      3.93e-13
  128      2.28e-06      6.45e-10      1.24e-16
scipy.integrate.simpson on 1001 samples: relative error 1.7e-13; quad: 0.0e+00
peak heating altitude by bisection on the central difference: 53.494 km (analytic 53.494 km)
```

and writes `reentry_heating.png` (stencil errors against h; quadrature errors against n; heating rate and deceleration against altitude with both peaks marked).

### Verify

Run `python -m pytest curriculum/phase-03/lesson-04-numerical-differentiation-integration`. The tests:

1. **Stencils**: orders 1, 2, 2 within 0.02 over h = 10⁻¹ to 10⁻³ and leading constants within 1 % at 10⁻³; exactness on lines and quadratics.
2. **Richardson**: central → order 4.0 ± 0.1, forward → 2.0 ± 0.1, exact on a cubic and on a pure h² error; order 0 raises.
3. **Optimal step**: h* formula; the best observed central step within a decade of it with error below 10⁻¹⁰; the rounding regime at 10⁻¹² is 100× worse.
4. **Trapezoid and Simpson**: exact classes; orders 2.0 ± 0.1 and 4.0 ± 0.2; Euler–Maclaurin constant within 2 %; equality with SciPy's sample rules to 10⁻¹³; odd n and n = 0 raise; spectral accuracy on e^{cos x} (10⁻¹⁴ with 16 panels).
5. **Romberg**: column 0 is the trapezoid rule, column 1 is Simpson to 10⁻¹³, the diagonal is monotone and reaches 10⁻¹⁰ while the trapezoid column is at 10⁻⁵.
6. **Gauss–Legendre** for n = 1, 2, 3, 5, 8, 16: nodes and weights equal `leggauss` to 10⁻¹⁴, symmetric, positive, summing to 2; exact for every monomial to degree 2n − 1 and not at 2n; more than five digits gained from 4 to 8 points; a mapped interval integrates x² exactly.
7. **Adaptive Simpson**: meets 10⁻⁶ and 10⁻¹⁰, costs more at the tighter tolerance, handles the 0.03-wide peak to 10⁻⁶ in under 5000 evaluations, raises at `max_depth` on |x − ½|^{−0.9}.
8. **Entry**: V → V_e at 300 km and V_e e^{−B} at sea level; h_a and a_max formulas; the analytic peaks are the maxima of the sampled profiles to 20 m; 5 g < a_max < 12 g; the numerical deceleration converges at order 2.0 ± 0.1 and matches to 10⁻⁶; the heat load matches `quad` to 10⁻⁹ and the erf form, with trapezoid order 2.0 ± 0.15 over 32–256 panels, Simpson at least 3.5 over 32–128, Gauss 10⁻¹⁰ at 64 points; 90 % of Q lies below 80 km.

### Compare

`scipy.integrate.quad` (QUADPACK's adaptive 21-point Gauss–Kronrod with Wynn extrapolation) returns the heat load to the last digit with an error estimate, and it is what you call when the integrand is a function. `scipy.integrate.trapezoid` and `simpson` are the rules for samples and agree with ours to 10⁻¹³ on the same grid; `scipy.integrate.fixed_quad` is Gauss–Legendre; `numpy.polynomial.legendre.leggauss` gives the nodes by the Golub–Welsch eigenvalue method rather than Newton, to the same 10⁻¹⁴. Nothing here is faster or more accurate than the library. What writing the rules bought: knowing that Simpson on 16 panels of a 7 km pulse spread over 120 km is *worse* than the trapezoid rule (7.4 × 10⁻² against 2.3 × 10⁻⁴), because its fourth-order law has not started while the trapezoid rule is already enjoying its Euler–Maclaurin cancellation at the nearly-flat ends; knowing from the Romberg table that two refinements that agree to 10⁻⁸ mean the answer is good to about 10⁻⁸; and knowing that a derivative by differences cannot give more than two thirds of the machine's digits, so a reported deceleration to ten figures from a sampled velocity is a reporting error.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Resolution study.** Tabulate the relative error of the heat load from 120 km for trapezoid, Simpson and Gauss at n = 8, 16, …, 512. Predict where the asymptotic orders 2 and 4 begin: when the panel width h = 120 km/n is well below the pulse width H = 7.2 km, i.e. n ≳ 64; predict that Gauss reaches the floor by n = 64; and predict the crossover where Simpson overtakes the trapezoid rule.
2. **Work–precision.** For each rule and for Romberg and adaptive Simpson, plot error against number of function evaluations. Predict the ordering at 10⁻⁶: Gauss < Romberg < Simpson < adaptive Simpson < trapezoid in evaluations.
3. **Range and peaks.** Vary the upper limit from 80 to 200 km and the parameters β ∈ {200, 300, 400} kg/m², γ ∈ {2°, 3°, 5°}. Predict from (5) that Q grows like 1/√B ∝ √(β sin γ) and that h_q = H ln 6B falls as β rises; confirm the peak altitudes by bisection on the central difference of q̇, as `main()` does.
4. **The derivative floor.** For the deceleration at 47.5 km from a central difference of V, sweep dh from 10⁴ m down to 10⁻³ m and find the step where the error turns up; predict it from h* = (3εF₀/F₃)^{1/3} with F₀ = V ≈ 5300 m/s and F₃ = |V‴| ≈ 1.5 × 10⁻¹¹ m/s per m³ (Problem 1 computes it): h* ≈ 0.6 m.

Raw data format: CSV tables `n, rel_err_trap, rel_err_simpson, rel_err_gauss` (1), `method, n_evals, rel_err` (2), `beta, gamma_deg, h_top_km, Q_MJ_m2, h_q_km, h_a_km` (3), `dh_m, abs_err` (4); plus `reentry_heating.png`.

What you should find: the trapezoid rule is already at 2.3 × 10⁻⁴ with 16 panels and settles to order 2 beyond 32 while Simpson is pre-asymptotic until 64 (apparent order near 8 over 32–128 as it "catches up") and overtakes the trapezoid rule between 32 and 64 panels; Gauss is at 4 × 10⁻¹³ with 64 points; the 10⁻⁶ threshold costs about 30 Gauss points, 44 Simpson panels, 188 trapezoid panels; Q scales as √(β sin γ) within 1 %, and the peak altitudes move by H ln(ratio) exactly; the central-difference error at 47.5 km falls as dh² down to about 1 m (8 × 10⁻¹⁰ m/s², 10⁻¹¹ relative) and rises again below 0.1 m.

## 5. Validate and improve

| quantity | theory (§2) | observed (code) | ratio or difference | cause |
|---|---|---|---|---|
| forward / central / second constants at h = 10⁻³ | 4.21 × 10⁻⁴ / 9.00 × 10⁻⁸ / 7.01 × 10⁻⁸ | 4.21 × 10⁻⁴ / 9.01 × 10⁻⁸ / 7.01 × 10⁻⁸ | < 1 % | next Taylor term |
| central optimum | h* = 8.7 × 10⁻⁶, error ~10⁻¹¹ | 1.1 × 10⁻¹¹ at h = 10⁻⁵ | — | model is a bound |
| Richardson of central | order 4 | 4.00 (h = 0.1 → 0.01) | 0 | even expansion |
| trapezoid, 128 panels | h²Δf′/12 = 2.00 × 10⁻⁵ | 2.00 × 10⁻⁵ | < 1 % | Euler–Maclaurin |
| Simpson order | 4 | 4.11 | 0.11 | h⁶ term still visible at n = 4 |
| trapezoid on e^{cos x}, 16 panels | exponential | 10⁻¹⁴ | — | periodic: all end terms vanish |
| Gauss 5-point nodes | ±0.9061798459, ±0.5384693101, 0 | same to 10⁻¹⁴ | — | Newton on the recurrence |
| a_max, h_a | 81.3 m/s², 45.6 km | 81.3 (central difference), 45.6 (argmax) | 10⁻⁶ relative, 20 m | dh = 10 m; grid 10 m |
| h_q | 53.494 km | 53.494 km (bisection on central difference) | < 1 m | — |
| Q from 120 km | 119.8708 MJ/m² (erf) | `quad` 0; Gauss(64) 4 × 10⁻¹³; Simpson(128) 6 × 10⁻¹⁰; trapezoid(128) 2.3 × 10⁻⁶ | — | spectral / h⁴ / h² |
| trapezoid vs Simpson, 16 panels | h ≫ H: pre-asymptotic | 2.3 × 10⁻⁴ vs 8.4 × 10⁻³ | Simpson 37× worse | Simpson's correction assumes an h⁴ law not yet valid; the trapezoid rule's end terms nearly vanish |

Every row where the theory makes a quantitative claim is met to the stated figure; the row that *looks* like a failure (Simpson losing to the trapezoid rule by a factor 37) is the one to understand. An error law such as (3) holds asymptotically, when the panel width is small compared with the scale of the integrand; before that the "fourth-order" correction (4S_{2h} − S_h)/3 that Simpson adds to the trapezoid rule is a correction of the wrong size. A resolution study is the only way to know which regime you are in, which is why every integration in this course is run at two resolutions and the difference reported. For a trajectory sampled by a simulation (part 6), the tool is Simpson on the samples plus a Richardson estimate from every other sample, and the requirement 10⁻⁶ on Q is met once the sample spacing is below about 2 km (the 64-panel line above), which any 6-DOF run at 0.1 s already provides.

What to improve: use Gauss or Romberg whenever the integrand is a function; use the trapezoid rule, not Simpson, on integrands that decay to zero at both ends; never estimate a second derivative from noisy samples by (2.1) without smoothing. What not to improve: the Allen–Eggers constants. A better atmosphere or a gravity term changes Q by tens of per cent (extension 3); a better quadrature rule changes it by 10⁻¹⁰. Spend the effort where the error is.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the three stencil constants, of (2) and (3), of Gauss exactness and of (4)–(5); `reentry_heating.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending, for a trajectory tool that receives V(t), h(t) as samples, the quadrature rule, the sample spacing and the error estimate it should report for Q, with the evidence.

### Adversarial self-checks

- **What if γ doubles to 6°?** B halves (141), h_q and h_a drop by H ln 2 = 5.0 km, a_max doubles to 16.6 g (it is proportional to sin γ), and Q falls by 1/√2 (shorter exposure wins over the higher rate). The pulse narrows in time but not in altitude, so the quadrature resolution requirement in h is unchanged.
- **What if β → ∞ (a cannonball)?** B → 0, V stays at V_e all the way down, a_max is unchanged (independent of β) but h_a → −∞, meaning the peak is below the ground: the vehicle hits at full speed. The model's silence about the ground is the first thing a real trajectory code adds.
- **Which assumption fails first in the real world?** Gravity and the constant γ, below about 35 km; then the exponential atmosphere (20 % density error at 60–80 km); the heating correlation last (15 %). None of them is a numerical error.
- **What measurement could falsify the model?** Accelerometer data from the capsule: a peak at the predicted altitude but 30 % lower would indicate lift or a shallower path; a heat-shield recession larger than Q/(heat of ablation) would indicate radiative heating or a wrong k.
- **Does the result survive a change in rule or step?** Q from trapezoid, Simpson, Gauss, Romberg, adaptive Simpson and `quad` agree to 10⁻⁹ at adequate resolution; the deceleration from differences agrees with the closed form to 10⁻⁶; change dh by a factor 2 and the central-difference result moves by 10⁻⁷ relative (its order-2 error), which is the uncertainty to report.
- **Which decision changes if the requirement on Q tightens to 10⁻⁹?** Trapezoid on samples cannot deliver it below 10⁴ panels; Simpson needs about 128; Gauss 64; the decision is to run Simpson at 1 km sample spacing with a Richardson check, and no change to the model, whose own error is 10⁻¹.

## Check yourself

1. Write the leading error term of the forward, central and second-difference stencils, and the order of each.
2. Why does the central difference have an optimal step near ε^{1/3} and what accuracy does that give? What is the corresponding figure for the second difference?
3. State the Richardson extrapolation formula for an estimate with error Chᵖ and say what it yields when applied to the trapezoid rule.
4. Give the trapezoid and Simpson error laws and explain why Simpson is exact for cubics although it integrates a parabola.
5. Under what conditions on the integrand does the trapezoid rule converge faster than any power of h? Give an example.
6. Why is n-point Gauss–Legendre exact for degree 2n − 1, and what prevents its use on sampled data?
7. A capsule with B = 281 and H = 7.2 km: at what altitudes do deceleration and heating peak, and which comes first?
8. The heat-load integrand is a pulse of width ~7 km on a 120 km range. Why does Simpson with 16 panels do worse than the trapezoid rule, and what single test tells you whether a rule is in its asymptotic regime?
