# P01.03  Integrals & the Fundamental Theorem

> An integral is a sum you can refine until it stops changing; the Fundamental Theorem tells you where it is heading.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the area and chordwise centroid of a NACA four-digit section in closed form by the Fundamental Theorem, and state where the "thickness applied vertically" simplification enters. 2. Implement and test left, right and midpoint Riemann sums, an accumulation function, and a substitution that removes an integrable singularity; measure the order of convergence of each. 3. Verify both parts of the Fundamental Theorem numerically: the derivative of the accumulated area is the local thickness, and the antiderivative at the ends reproduces the sum. 4. Decide, from an error-versus-panel-count plot, whether an integrand's singularity or its smoothness limits the accuracy, and what to do about it. |
| **Prerequisites** | `P01.01` (derivatives, convergence order, log-log fits), `P01.02` (optimisation vocabulary), `P00.03` (NumPy), `P00.09` (SymPy `integrate`). Skill check: you can integrate xⁿ for any real n ≠ −1, sketch √x near zero, and evaluate a polynomial at a point. |
| **Engineering question** | What are the cross-sectional area and the centroid of a NACA 0012 wing section (for its mass per unit span and the location of its bending neutral axis), and how many panels does a numerical integration need before the answer is good to five figures? |
| **Assumptions and limits** | The NACA four-digit thickness polynomial as given (an open trailing edge of 0.21 % c unless the closed-edge coefficient is used); for cambered sections the thickness is added vertically to the camber line, not perpendicular to it as in the NACA construction, which changes the area and centroid at second order in the camber slope; the section is solid (a real wing is a skin, spars and ribs, so the solid centroid is a reference, not the structural neutral axis). Riemann sums assume a bounded integrand on a finite interval; the leading-edge √x slope is infinite but integrable, and it slows convergence. |
| **Analytical work** | Riemann sums as limits; the definite integral and its properties; both parts of the Fundamental Theorem; antiderivatives of powers including x^(1/2); the closed-form area and centroid; the error law of the midpoint rule from Taylor's expansion; the substitution x = c u² and why it restores second-order convergence; the net change theorem on a vertical-speed record. |
| **Computational deliverable** | `starter/main.py`: `naca4_thickness`, `naca4_camber`, `riemann_sum`, `cumulative_integral`, `section_area_exact`, `section_centroid_x_exact`, `symbolic_section_area`, `definite_from_antiderivative`, `section_area`, `section_area_substituted`, `section_centroid`, `net_change`, `total_variation`, `convergence_study`, `observed_order`. Tests prove: left and right sums bracket an increasing integrand, midpoint is exact for linear functions, orders 1 and 2 on a smooth integrand, order 1.5 on the leading-edge singularity and 2 after substitution, the accumulation function's derivative equals the integrand, the antiderivative route and the sum agree, the symbolic area equals the closed form, area and centroid scale correctly with t and c, and the camber lifts the centroid linearly. |
| **Practical activity** | **computational experiment**: a convergence study of three Riemann rules on the section area over seven doublings of n, with and without the leading-edge substitution, and a numerical check of the Fundamental Theorem part 1. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `airfoil_section.png`, `riemann_convergence.png` with its table and fitted orders, the hand derivation of A and x̄, and a short note explaining the 1.5 and the 2. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the closed-form derivation is complete with the x^(1/2) antiderivative shown; all public tests pass; the hand value of A agrees with the code to five figures; the convergence note names the cause of each observed order with a number. |
| **References** | `@openstax-calculus-volume-1` §4.10 "Antiderivatives", §5.1 "Approximating Areas", §5.2 "The Definite Integral", §5.3 "The Fundamental Theorem of Calculus", §5.4 "Integration Formulas and the Net Change Theorem", §5.5 "Substitution" (all required-free). The same material appears as `@openstax-calculus-volume-2` Chapter 1 "Integration", §1.1–§1.5 (supplementary-free; either volume may be read). |
| **Extensions** | ★ Implement the true NACA construction (thickness perpendicular to the camber line) as a parametric curve and measure the change in area and centroid for NACA 2412 and 6412. ★ Compute the second moment of area of the section about its centroid and the radius of gyration (Problem 7 starts this). |
| **Time** | unmeasured |

## 1. The engineering problem

A wing is a beam. Before any structural calculation you need the geometry of its cross-section: the **area** (for mass per unit span, and for a first estimate of the stiffness of a solid or foam-filled section), and the **centroid** (the point about which bending is computed; for a symmetric section it lies on the chord line, but *where* along the chord is not obvious because the nose is fat and the tail is thin). The NACA four-digit sections, published in the 1930s and still used for tails, struts and model aircraft, are defined by a formula: the half-thickness is a polynomial in x/c with one awkward term, 0.2969 √(x/c), which gives the rounded nose. The area is the integral of that formula; the centroid is the integral of x times it, divided by the area.

So the engineering question is: what are A and x̄ for a NACA 0012 of chord c, as formulas and as numbers, and if the formula were *not* available (a measured section, a CFD mesh, a table from a report), how many sample points would a numerical integral need to get the same answer to five figures? The second question has a surprise in it: the square root that draws the nose also slows every simple numerical integrator from second order to order 1.5, and knowing why, and how to fix it with a change of variable, is worth more than the airfoil itself.

Success criteria, written first: the closed-form area agrees with the published value for the NACA 0012 (about 0.082 c²) and with SymPy exactly; a midpoint sum with 1280 panels agrees with the closed form to 5 × 10⁻⁷ c²; the observed orders are 1.5 for the midpoint rule on the raw integrand and 2.0 after the substitution x = c u²; and the derivative of the accumulated area matches the local thickness to 10⁻⁶.

## 2. Understand and derive

### 2.1 The section

The NACA four-digit half-thickness distribution, for 0 ≤ x ≤ c, is

y_t(x) = 5 t c [0.2969 √(x/c) − 0.1260 (x/c) − 0.3516 (x/c)² + 0.2843 (x/c)³ − 0.1015 (x/c)⁴],   (1)

where t is the maximum thickness as a fraction of chord (0.12 for a NACA 0012). The factor 5t makes the maximum of 2y_t equal to t c, at x = 0.30 c. The coefficients sum to 0.0021, so the trailing edge is not quite closed: 2y_t(c) = 0.0021 × 5 t c = 0.00252 c for t = 0.12. Replacing −0.1015 by −0.1036 closes it; the code offers both.

```
  y/c
 0.06 ┤      _______
      │    /         \______
      │   /                 \_______
 0    ┼──┼────────┼────────┼────────┼──▶ x/c       NACA 0012: thickness t = 0.12 c at x = 0.30 c
      │   \                 /‾‾‾‾‾‾‾
      │    \_________/‾‾‾‾‾‾
-0.06 ┤
      0        0.3      0.5      0.75     1
```

**Assumption 1 (solid section).** Area and centroid are computed for the full region between the upper and lower surfaces. For a cambered section (camber line y_c(x), maximum camber m at position p), **Assumption 2 (vertical thickness):** upper surface y_c + y_t, lower y_c − y_t. The NACA definition adds the thickness perpendicular to the camber line; the difference in the area is of order (y_t · dy_c/dx)² and is below 0.1 % for 2 % camber, but it is an assumption and it is stated.

### 2.2 Riemann sums: the area as a limit of rectangles

Split [a, b] into n panels of width h = (b − a)/n. On each panel approximate f by a constant and add the rectangles (`@openstax-calculus-volume-1` §5.1):

- **left sum** L_n = h Σ f(x_i), i = 0…n−1;
- **right sum** R_n = h Σ f(x_i), i = 1…n;
- **midpoint sum** M_n = h Σ f(x_i + h/2), i = 0…n−1.

The **definite integral** is the common limit, ∫ₐᵇ f dx = lim_{n→∞} of any of them (§5.2), and it exists for every continuous f and for many discontinuous ones. Three properties you will use: linearity (∫(αf + βg) = α∫f + β∫g), additivity over intervals (∫ₐᶜ = ∫ₐᵇ + ∫ᵇᶜ), and the sign flip ∫ᵇₐ = −∫ₐᵇ. For an increasing f, L_n ≤ ∫ ≤ R_n: the two sums *bracket* the answer, which is a test.

Worked by hand, n = 4 midpoint sum of x² on [0, 1]: h = 0.25, midpoints 0.125, 0.375, 0.625, 0.875, squares 0.015625, 0.140625, 0.390625, 0.765625, sum 1.3125, times h: **0.328125**. The exact value is 1/3 = 0.333333; the error is 0.005208.

### 2.3 The Fundamental Theorem, both parts

Define the **accumulation function** F(x) = ∫ₐˣ f(s) ds: the area from a up to a moving right end. **Part 1** (§5.3): if f is continuous, F is differentiable and F′(x) = f(x). The picture: moving the right end by dx adds a sliver of height f(x) and width dx. **Part 2**: if G is *any* antiderivative of f (G′ = f), then ∫ₐᵇ f dx = G(b) − G(a). Part 1 says integrals undo derivatives; part 2 says derivatives undo integrals, and it turns an infinite sum into two evaluations.

Antiderivatives you need here (§4.10): ∫xⁿ dx = xⁿ⁺¹/(n + 1) for every real n ≠ −1, including n = ½, so ∫√x dx = (2/3) x^(3/2). The constant of integration cancels in G(b) − G(a) and is irrelevant for a definite integral.

### 2.4 The area, in closed form

A = 2 ∫₀ᶜ y_t dx (both surfaces). Substitute ξ = x/c, dx = c dξ, and integrate (1) term by term with the power rule:

∫₀¹ ξ^(1/2) dξ = 2/3,  ∫₀¹ ξ dξ = 1/2,  ∫₀¹ ξ² dξ = 1/3,  ∫₀¹ ξ³ dξ = 1/4,  ∫₀¹ ξ⁴ dξ = 1/5.

So A = 10 t c² [0.2969 × 2/3 − 0.1260/2 − 0.3516/3 + 0.2843/4 − 0.1015/5] = 10 t c² [0.197933 − 0.063000 − 0.117200 + 0.071075 − 0.020300] = 10 t c² × 0.0685083 = **0.685083 t c²**.   (2)

For the NACA 0012 with c = 1 m: A = 0.082210 m². SymPy gives the same as a rational, 8221/12000 t c² (exact because the coefficients are decimals with four places), and the test asks for agreement to 10⁻¹². Units: y_t is a length, dx a length, A an area; t is dimensionless. Scaling: A ∝ t and A ∝ c², both tested. With the closed trailing edge A = 0.081706 c², 0.61 % less.

### 2.5 The centroid

By symmetry ȳ = 0 for a symmetric section. The chordwise centroid is the first moment divided by the area:

x̄ = (1/A) ∫ x · 2y_t dx = c [Σ aₖ/(eₖ + 2)] / [Σ aₖ/(eₖ + 1)],

where aₖ are the coefficients and eₖ the exponents in (1), because ∫₀¹ ξ · ξ^e dξ = 1/(e + 2). The numerator sum is 0.2969 × 2/5 − 0.1260/3 − 0.3516/4 + 0.2843/5 − 0.1015/6 = 0.118760 − 0.042000 − 0.087900 + 0.056860 − 0.016917 = 0.0288033, so

**x̄ = 0.0288033/0.0685083 c = 0.42044 c**.   (3)

The centroid sits ahead of mid-chord, at 42 % c, because the front of the section is thicker. For a cambered section with vertical thickness, x̄ is unchanged (the thickness distribution is the same) and ȳ = (1/A) ∫ y_c · 2y_t dx, which is linear in the camber m: for a NACA 2412, ȳ = 0.01543 c, 77 % of the maximum camber.

### 2.6 How wrong is a midpoint sum? The error law

Expand f about the panel midpoint x_m: f(x) = f(x_m) + f′(x_m)(x − x_m) + ½ f″(x_m)(x − x_m)² + …. Integrate over the panel: the linear term integrates to zero by symmetry (which is why the midpoint rule is exact for straight lines), and the quadratic term gives f″(x_m) h³/24. Summing n panels and writing Σ f″(x_m) h ≈ ∫ f″ = f′(b) − f′(a):

**∫ₐᵇ f dx − M_n ≈ (h²/24) [f′(b) − f′(a)]**.   (4)

Second order: doubling n divides the error by four. Check on the hand example: f = x², f′(1) − f′(0) = 2, h = 0.25: (0.0625/24) × 2 = 0.005208, the error found in §2.2 to all the digits shown. For left and right sums the same expansion gives an error of ∓(h/2)[f(b) − f(a)] + O(h²): first order, and opposite in sign, so (L_n + R_n)/2 is second order (that is the trapezoid rule, P01.04).

### 2.7 The leading-edge singularity, and the substitution that removes it

Law (4) needs f′ to be finite at the ends. At the leading edge y_t ∝ √x, so f′(0) = ∞ and (4) breaks down. Look at the first panel alone: ∫₀ʰ √x dx = (2/3)h^(3/2), while the midpoint rule gives h √(h/2) = 0.7071 h^(3/2); the error of that one panel is 0.040 h^(3/2), and the later panels add smaller contributions of the same type. The total error goes as h^(3/2): **order 1.5**, not 2. The code measures |A − M_n| ≈ 0.022 h^(3/2) c² for the NACA 0012, falling from 7.4 × 10⁻⁴ at n = 10 to 4.8 × 10⁻⁷ at n = 1280, a factor of 2.83 = 2^(1.5) per doubling.

The cure is **substitution** (§5.5), the integral's version of the chain rule: if x = g(u), then ∫ f(x) dx = ∫ f(g(u)) g′(u) du. Take x = c u², dx = 2cu du. Then √(x/c) = u and every term of (1) becomes a power of u times 2cu:

A = 2 ∫₀¹ y_t(cu²) · 2cu du = 10 t c² ∫₀¹ [0.2969 u − 0.1260 u² − 0.3516 u⁴ + 0.2843 u⁶ − 0.1015 u⁸] 2u du,

a polynomial in u with nothing singular anywhere. For t = 0.12, c = 1 the integrand is g(u) = 0.71256 u² − 0.3024 u³ − 0.84384 u⁵ + 0.68232 u⁷ − 0.2436 u⁹, and (4) applies with g′(1) − g′(0) = −1.1174: the predicted error at n = 1280 is 1.1174/(24 × 1280²) = 2.84 × 10⁻⁸, and the code measures 2.84 × 10⁻⁸. Same area, same n, seventeen times smaller error, and the order is back to 2.00.

### 2.8 Net change: what the integral of a rate means

If a quantity Q has rate dQ/dt = r(t), then Q(b) − Q(a) = ∫ₐᵇ r dt (§5.4, the **net change theorem**; it is part 2 of the Fundamental Theorem read as physics). Two quantities with the same rate differ by a constant, so the integral gives the *change*, not the value; the initial value must be supplied. Applied to a glider whose vertical speed in a weak thermal is v(t) = 2 sin(πt/60) − 0.8 m/s over 120 s: the net height change is ∫v dt = −96.0 m (the sine averages to zero over its 120 s period; only the −0.8 m/s remains), while the total vertical distance ∫|v| dt = 165.2 m. A variometer integrates v; an altimeter reads the net change; the pilot's workload is closer to the total variation.

### 2.9 Where this stops being valid

- **The polynomial is the section.** Real wings are built to tolerance; a measured section differs from (1) by 0.1–0.5 mm on a 1 m chord, which changes A by about 0.5 %.
- **Solid section.** The structural neutral axis of a built-up wing is set by spar caps and skin, not by the solid area. The solid centroid is the right reference for a foam-core model wing and a first guess for anything else.
- **Vertical thickness.** For camber above 4 % the perpendicular construction shifts the nose forward and the area by more than 0.3 %; use the extension.
- **Riemann sums on singular integrands.** A √x endpoint costs half an order; a 1/√x endpoint (still integrable) costs more, and a 1/x endpoint is not integrable at all. The substitution in §2.7 handles the first two; the third needs the limit process of P01.04.

## 3. Build the computational model

Open `starter/main.py`. The NACA coefficients are at the top as arrays.

### Step 1: the geometry

```python
import numpy as np

NACA_COEFFS = np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1015])
NACA_EXPONENTS = np.array([0.5, 1.0, 2.0, 3.0, 4.0])

def naca4_thickness(x, t, c=1.0, closed_te=False):
    x = np.asarray(x, dtype=float)
    if t <= 0.0 or c <= 0.0:
        raise ValueError("thickness ratio and chord must be positive")
    if np.any(x < 0.0) or np.any(x > c):
        raise ValueError("x must lie in [0, c]")
    xi = x / c
    coeffs = NACA_COEFFS_CLOSED if closed_te else NACA_COEFFS
    out = 5.0 * t * c * sum(a * xi**e for a, e in zip(coeffs, NACA_EXPONENTS))      # (1)
    return float(out) if np.ndim(out) == 0 else out
```

Keeping the coefficients and exponents as parallel arrays means the closed-form area and centroid are one line each: `10 t c² Σ a/(e+1)` and `c Σ a/(e+2) / Σ a/(e+1)`. `naca4_camber` is the two-piece parabola with `np.where`, continuous at x = pc (a test checks it).

### Step 2: Riemann sums

```python
def riemann_sum(f, a, b, n, rule="midpoint"):
    if n < 1:
        raise ValueError("n must be at least 1")
    if rule not in ("left", "right", "midpoint"):
        raise ValueError("rule must be left, right or midpoint")
    h = (b - a) / n
    edges = a + h * np.arange(n + 1)
    pts = {"left": edges[:-1], "right": edges[1:], "midpoint": 0.5 * (edges[:-1] + edges[1:])}[rule]
    return float(h * np.sum(f(pts)))
```

One vectorised call to `f`, no Python loop over panels. Note `b < a` is allowed and gives the negative, as the integral does.

### Step 3: the accumulation function and the antiderivative route

```python
def cumulative_integral(f, a, xs, n=400):
    xs = np.atleast_1d(np.asarray(xs, dtype=float))
    return np.array([riemann_sum(f, a, x, n, "midpoint") if x != a else 0.0 for x in xs])

def definite_from_antiderivative(expr, x, a, b):
    big_f = sp.integrate(expr, x)
    return float(big_f.subs(x, b) - big_f.subs(x, a))
```

The first is part 1 of the Fundamental Theorem as a function (the test differentiates it and compares with f); the second is part 2 (the test compares it with a Riemann sum). `symbolic_section_area` integrates (1) with SymPy using `sp.Rational("0.2969")` so the result is exact.

### Step 4: area and centroid, numerically

```python
def section_centroid(t, c=1.0, m=0.0, p=0.4, n=1000, closed_te=False):
    yt = lambda x: naca4_thickness(x, t, c, closed_te)
    area = 2.0 * riemann_sum(yt, 0.0, c, n)
    sx = 2.0 * riemann_sum(lambda x: x * yt(x), 0.0, c, n)
    sy = 2.0 * riemann_sum(lambda x: naca4_camber(x, m, p, c) * yt(x), 0.0, c, n)
    return sx / area, sy / area
```

`section_area_substituted` integrates `naca4_thickness(c*u**2, ...) * 2*c*u` over u ∈ [0, 1]: the substitution of §2.7 in one lambda.

### Step 5: convergence tools

`convergence_study(f, exact, a, b, ns, rule)` returns |M_n − exact| for a list of n, and `observed_order(ns, errs)` fits log(err) against log(n) and returns the negative slope, as in P01.01. `net_change` and `total_variation` are midpoint sums of r and |r|.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-03-integrals-fundamental-theorem`. The tests:

1. **Geometry.** y_t(0) = 0; maximum 2y_t = t c to 0.2 % at x = 0.30 c; the closed-edge coefficient gives y_t(c) = 0; y_t scales with c and with t; camber is continuous at p, peaks at m c, vanishes at both ends.
2. **Riemann sums.** Left < exact < right for eˣ; midpoint exact for 3x + 1; reversed limits flip the sign; n = 0 and an unknown rule raise.
3. **Orders.** Left 1.0 and midpoint 2.0 on cos x; midpoint 1.5 on the raw thickness and 2.0 after substitution, with every substituted error smaller than the raw one.
4. **Fundamental Theorem.** The central difference of `cumulative_integral` matches 2y_t to 0.2 % at two step sizes; F increases; the antiderivative of x³ − 2x + 1 gives 2 on [0, 2] and matches a 4000-panel sum to 10⁻⁶; the SymPy area is a number times t c² and equals (2) to 10⁻¹² for both trailing edges.
5. **Area and centroid.** Numeric area matches closed form to 10⁻⁶ (raw, n = 20 000) and 10⁻⁷ (substituted, n = 4000); A ∝ c² and ∝ t; closing the trailing edge reduces A; ȳ = 0 for the symmetric section and x̄ ∈ (0.40, 0.45) c; camber leaves x̄ unchanged and lifts ȳ linearly.
6. **Net change.** ∫v = −96 m exactly to 10⁻⁹ (the sine integrates to zero); total variation exceeds the net; ∫cos over a period is zero.

### Compare

Three routes to A for the NACA 0012, c = 1 m: the closed form (2) gives 0.0822100; SymPy's exact integration gives 8221/12000 × 0.12 = 0.0822100; the midpoint sum with n = 1280 gives 0.0822095 (raw) and 0.08221003 (substituted). `scipy.integrate.quad`, the library's adaptive integrator, returns 0.0822100 with an error estimate of 10⁻¹⁶ after 189 function evaluations; it handles the √x endpoint by subdividing the interval where the error estimate is largest, which is the automated version of what the substitution does by hand, and it needs nine times as many evaluations as the substituted midpoint rule at n = 20 to reach the same accuracy. When a closed form exists it is the reference; when it does not (a measured section), `quad` on an interpolant or the midpoint rule with the substitution are the tools, and the convergence plot is the evidence that they have converged.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** A convergence study with a prediction before each run.

1. For n = 10, 20, …, 1280 compute |A − S_n| for the left, right and midpoint rules on f(x) = 2y_t(x) (NACA 0012, c = 1), and for the midpoint rule on the substituted integrand g(u).
2. **Before running, write down** the expected order for each: 1 for left and right (first order), 1.5 for midpoint on f (the singularity), 2 for midpoint on g. Then run and fit `observed_order` on n ≥ 40.
3. Plot the four error curves on log-log axes with reference lines n^(−1.5) and n^(−2); save as `riemann_convergence.png`.
4. Check the Fundamental Theorem part 1: differentiate `cumulative_integral` by a central difference with h = 10⁻³ at nine points in [0.1, 0.9] and record the largest |dF/dx − 2y_t|.

Raw data format: CSV with columns `n, err_left, err_right, err_midpoint, err_substituted`, and one line for the FTC check.

The author's run:

| n | left | right | midpoint | midpoint, u-subst |
|---|---|---|---|---|
| 10 | 2.6 × 10⁻³ | 2.3 × 10⁻³ | 7.4 × 10⁻⁴ | 4.6 × 10⁻⁴ |
| 40 | 3.3 × 10⁻⁴ | 2.7 × 10⁻⁴ | 8.9 × 10⁻⁵ | 2.9 × 10⁻⁵ |
| 160 | 4.5 × 10⁻⁵ | 2.9 × 10⁻⁵ | 1.1 × 10⁻⁵ | 1.8 × 10⁻⁶ |
| 640 | 6.6 × 10⁻⁶ | 2.6 × 10⁻⁶ | 1.4 × 10⁻⁶ | 1.1 × 10⁻⁷ |
| 1280 | 2.6 × 10⁻⁶ | 6.4 × 10⁻⁷ | 4.8 × 10⁻⁷ | 2.8 × 10⁻⁸ |

Fitted orders on n ≥ 40: left 1.40, right 1.73, midpoint 1.51, substituted 2.00. The FTC check gives 4.6 × 10⁻⁷. Two of these numbers will not match your prediction, and part 5 explains them.

## 5. Validate and improve

| quantity | predicted | measured | cause of the difference |
|---|---|---|---|
| A, NACA 0012 (c = 1 m) | 0.0822100 (closed form, SymPy) | 0.0822095 (midpoint, n = 1280) | the 0.022 h^(3/2) singular error, 4.8 × 10⁻⁷ |
| midpoint order, raw | 1.5 | 1.51 | the h² smooth term still contributes a little |
| midpoint error, substituted, n = 1280 | 2.84 × 10⁻⁸ from (4) | 2.84 × 10⁻⁸ | none: (4) is exact to the digits shown |
| left/right order | 1.0 | 1.40 / 1.73 | see below |
| x̄ | 0.42044 c | 0.42044 c (n = 4000) | — |
| dF/dx against 2y_t | 0 | 4.6 × 10⁻⁷ | central-difference truncation, h²/6 f‴, plus the sum's own error |

The left and right orders are the discrepancy worth understanding. Law (4)'s cousin for one-sided sums says the first-order error is ∓(h/2)[f(b) − f(a)]. Here f(c) − f(0) = 2y_t(c) = 0.00252 c, tiny because the trailing edge is almost closed, so the first-order term is only 0.0013 h and the singular 0.04 h^(3/2) term competes with it at every n in the table: the left sum (where both terms have the same sign) shows a blend near 1.4, and the right sum (where they partly cancel) shows 1.73 and is *more* accurate than the midpoint rule at n = 1280. A fit that mixes two power laws gives a slope that is neither. The lesson for practice: an observed order that is not an integer or a half-integer usually means two error terms of comparable size, and the fix is to look at the signs, not to refine further.

What the agreement does not show: that the NACA 0012 of a real model wing has this area. Published values for the section (about 0.0822 c²) agree with (2) because they come from the same formula; only a measurement of a built wing validates the geometry, with the tolerances named in §2.9.

Improvements worth making: use the substituted integrand for every airfoil integral in the rest of the course (P01.04 adds Simpson's rule, which on g(u) reaches 10⁻⁸ with 64 panels and 6 × 10⁻¹⁰ with 128); compute the centroid of the *skin* (a line integral along the surface, P01.11) for a built-up wing. Improvement worth *not* making: pushing n to 10⁶ on the raw integrand to recover the digits the singularity costs; the substitution gets them with 1000.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (2) and (3) with every power-rule step; `airfoil_section.png` and `riemann_convergence.png` with the CSV; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) explaining the three observed orders and stating which integrand you would hand to a numerical integrator and why.

### Adversarial self-checks

- **What if the thickness doubles?** A doubles, x̄ is unchanged (the shape of the distribution is the same), ȳ of a cambered section is unchanged too (it is a weighted mean of y_c with weights ∝ y_t). What if the chord doubles? A quadruples, x̄ doubles. Check all four with the closed forms.
- **What if n goes to 1? To 10⁷?** With n = 1 the midpoint sum is c · 2y_t(c/2) = 0.1059 c², 29 % high (the mid-chord thickness is far above the average); with 10⁷ panels the raw midpoint error is 7 × 10⁻¹³ c², below the rounding error of summing 10⁷ terms (about 10⁻¹⁵ × 10⁷ × 0.08 = 10⁻⁹): a limit of the arithmetic, not the method.
- **Which assumption fails first in the real world?** Solid section: a built wing's mass and bending centroid are set by its structure, not its outline. The outline's centroid is still needed, as the reference for aerodynamic centre and structural axis positions.
- **What measurement could falsify the model?** Weigh a foam-core wing of known span and density: the mass per span should be ρ_foam × 0.0822 c². A 2 % disagreement points at the foam density; a 10 % one at the section.
- **Does the result survive a change in rule or in n?** Midpoint, left, right and substituted all converge to the same 0.08221 c²; they differ in how fast. Simpson (P01.04) and `quad` agree to 10⁻¹⁰.
- **Which decision would change if the uncertainty doubled?** If a built section's profile were uncertain by 1 mm instead of 0.5 mm on a 1 m chord, the area would be uncertain by 1 % instead of 0.5 %, and the mass budget of a 20-wing production run would need twice the margin. The integral's own error (10⁻⁷) is irrelevant next to either.

## Check yourself

1. State both parts of the Fundamental Theorem in one sentence each.
2. Why is the midpoint rule exact for a straight line but the left rule not?
3. The midpoint sum of x² on [0, 1] with four panels is 0.328125. Use the error law (h²/24)[f′(b) − f′(a)] to predict the error, and compare with 1/3 − 0.328125.
4. Why does the midpoint rule converge as n^(−1.5) on the NACA thickness instead of n^(−2), and which single term of (1) is responsible?
5. What substitution turns √(x/c) into a polynomial, and what happens to dx?
6. The area of a NACA 0012 is 0.0822 c². Without recomputing, give the area of a NACA 0018 and of a NACA 0012 with twice the chord.
7. The chordwise centroid of the symmetric section is at 0.42 c. Why is it ahead of mid-chord, and does adding 2 % camber (vertical thickness) move it?
8. A glider's vertical speed integrates to −96 m over two minutes while ∫|v| dt = 165 m. What does each number mean physically, and which instrument reads which?
