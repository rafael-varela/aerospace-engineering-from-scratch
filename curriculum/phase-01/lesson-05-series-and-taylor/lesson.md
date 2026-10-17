# P01.05  Sequences, Series & Taylor Expansions

> Every "≈" in an engineering formula is a truncated series; the remainder term is the engineer's promise about how wrong it is.

| | |
|---|---|
| **Learning outcomes** | 1. Decide whether a series converges (geometric tail, ratio test, alternating series, comparison with the harmonic series) and bound how far a partial sum is from the limit. 2. Build Taylor polynomials by hand and with SymPy, state Taylor's theorem with the Lagrange remainder, and use it to bound the error of sin θ ≈ θ, cos θ ≈ 1 − θ²/2, tan θ ≈ θ and 1/√(1 − M²) ≈ 1 + M²/2. 3. Implement partial sums, remainder bounds and a small-angle error explorer, and measure how the error grows with angle and how many terms a series needs at a given argument. 4. Recognise when a convergent series is useless in floating point (cancellation at large arguments) and what to do instead (argument reduction, reciprocal), and when a divergent series pretends to converge (the harmonic series). |
| **Prerequisites** | `P01.01` (derivatives of all orders, orders of convergence), `P01.03` (the integral test needs ∫dx/x and ∫dx/x²), `P00.09` (SymPy `diff`, `series`). Skill check: you can differentiate sin, cos, eˣ, (1 + x)^α repeatedly; you know n! and can evaluate 7!; you can sum 1 + ½ + ¼ + … in your head. |
| **Engineering question** | The linearised equations of flight dynamics replace sin α by α and cos α by 1; the Prandtl–Glauert rule replaces 1/√(1 − M²) by a few terms; the inertia swing test replaces a pendulum's period by its small-amplitude value. At what angle, Mach number or amplitude does each replacement cost 1 %, and how is that number computed *before* the exact answer is known? |
| **Assumptions and limits** | The functions expanded are analytic where they are used (sin, cos, exp everywhere; 1/√(1 − M²) only for M < 1; ln(1 − f) only for f < 1); the remainder bounds assume a bound on the next derivative, which is easy for sin and cos (≤ 1) and must be computed for the others. Floating-point partial sums are exact in the mathematics and wrong in the arithmetic once the terms are far larger than the sum. The engineering models the approximations serve (linearised dynamics, Prandtl–Glauert, the pendulum) carry their own assumptions, which are Phase 05 and Phase 15 material. |
| **Analytical work** | Geometric series and its tail; the harmonic series and the integral test; ratio test and the radius of convergence; the alternating-series bound; Taylor's theorem with the Lagrange remainder; Maclaurin series of sin, cos, exp, 1/(1 − x), (1 + x)^α; the small-angle error laws θ²/6, θ²/2, θ⁴/24, θ²/3; the binomial series of the Prandtl–Glauert factor; the series for the pendulum period; Euler's formula from the exp series. |
| **Computational deliverable** | `starter/main.py`: `partial_sums`, `geometric_terms`, `geometric_sum`, `harmonic_partial_sums`, `ratio_test_estimate`, `radius_of_convergence`, `alternating_series_bound`, `taylor_polynomial`, `taylor_coefficients`, `series_terms`, `maclaurin_partial_sum`, `lagrange_remainder_bound`, `terms_needed`, `small_angle_error`, `small_angle_error_leading`, `max_angle_for_error`, `binomial_coefficients`, `prandtl_glauert`, `prandtl_glauert_series`, `pendulum_period_ratio`, `pendulum_period_ratio_exact`, `euler_formula_residual`. Tests prove: geometric partial sums obey the tail bound and diverge outside |r| < 1; the harmonic sums grow like ln n; ratio and radius estimates reproduce known radii (½, 1, ∞, 1); the alternating bound holds for sin and cos; Taylor polynomials match SymPy and have the right value, slope and curvature at the centre; Maclaurin sums converge to NumPy for moderate arguments and fail by cancellation at 40; the Lagrange bound holds and `terms_needed` is minimal; the small-angle errors follow their leading Taylor term with orders 2 and 4; the inverse angles are ordered cos < tan < sin < cos₂; the binomial series reproduces Prandtl–Glauert monotonically inside M < 1 and slowly near it; the pendulum series converges to the elliptic integral; Euler's formula holds term by term. |
| **Practical activity** | **computational experiment**: the small-angle error explorer (error against angle for four approximations, with the leading Taylor term overlaid), a terms-needed study across arguments, and a cancellation experiment at large arguments. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `small_angle_error.png` with its table of 1 % and 0.1 % angles, the hand derivations of the four small-angle error laws from Taylor's theorem, and a short note on when a convergent series should not be summed. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: Taylor's theorem is stated with the remainder and applied to at least two approximations with a numerical bound; all public tests pass; the 1 % angles agree with the code to three figures; the note names the cancellation mechanism with a number. |
| **References** | `@openstax-calculus-volume-2` §5.1 "Sequences", §5.2 "Infinite Series", §5.3 "The Divergence and Integral Tests", §5.4 "Comparison Tests", §5.5 "Alternating Series", §5.6 "Ratio and Root Tests", §6.1 "Power Series and Functions", §6.2 "Properties of Power Series", §6.3 "Taylor and Maclaurin Series", §6.4 "Working with Taylor Series" (all required-free). The must-knows are the ratio test, the alternating-series bound and Taylor's theorem with remainder. |
| **Extensions** | ★ Implement argument reduction for sin (subtract the nearest multiple of 2π, then use symmetry to land in [0, π/4]) and show the series then needs at most 8 terms for 10⁻¹⁶ at any argument. ★ Derive the series for the pendulum period from the binomial expansion of the elliptic integrand (Problem 7) and extend it to the amplitude at which it stops converging. |
| **Time** | unmeasured |

## 1. The engineering problem

Open any flight-dynamics text at the linearised equations of motion and you will find sin α replaced by α, cos α by 1, and tan γ by γ. Open an aerodynamics text at compressibility and you will find the lift-curve slope multiplied by 1/√(1 − M²), or by its first two terms 1 + M²/2. Open a test report on the measurement of an aircraft's moment of inertia (the swing test of P01.04) and you will find the period of a pendulum taken as 2π√(L/g), which is only true for a vanishing amplitude. All three are **truncated Taylor series**, and all three come with an unstated promise: the terms dropped are small. This lesson makes the promise explicit. Taylor's theorem says exactly what is dropped, and the Lagrange remainder bounds it with the next derivative; for an alternating series the bound is simply the next term. The practical output is a number: sin θ ≈ θ is within 1 % up to 14.0°, cos θ ≈ 1 only up to 8.1°, cos θ ≈ 1 − θ²/2 up to 37.9°, tan θ ≈ θ up to 9.9°; 1 + M²/2 is within 1 % of the Prandtl–Glauert factor up to M = 0.42 and useless at M = 0.9; the pendulum period is within 1 % of 2π√(L/g) up to an amplitude of 22.8°, and because the inertia goes as the period squared, that is 2 % in I_xx.

Behind the approximations is the question of whether an infinite sum means anything at all. 1 + ½ + ¼ + … means 2; 1 + ½ + ⅓ + … means nothing, although after a million terms it has only reached 14.4 and looks as settled as anything; and the series for 1/√(1 − M²) means 1/√(1 − M²) for M < 1 and nothing at all for M > 1, which is the mathematics' way of saying that subsonic theory stops at the speed of sound. Deciding convergence, and estimating the radius inside which a series is trustworthy, is the first half of the lesson; the error of a truncation is the second; and what floating point does to a convergent series at a large argument is the warning at the end.

Success criteria, written first: the partial-sum tail of the geometric series is bounded by r^n/(1 − r) and the harmonic sums track ln n + 0.577; the ratio test reproduces the radii ½, 1 and ∞ for 1/(1 − 2x), 1/(1 + x²) and eˣ; every Taylor polynomial agrees with SymPy; the four small-angle error curves lie on their leading-term lines θ²/6, θ²/2, θ⁴/24 and θ²/3 to 2 % below 4°; the Lagrange and alternating bounds are never violated; and the 60-term series for sin 40 is wrong by more than 1 for a reason you can name.

## 2. Understand and derive

### 2.1 Sequences, series and the two series that teach everything

A **sequence** a₀, a₁, … converges to L when its terms eventually stay within any ε of L. A **series** Σaₖ converges when its sequence of partial sums Sₙ = a₀ + … + aₙ converges; the sum *is* that limit, and nothing else.

**Geometric series.** For Sₙ = 1 + r + … + rⁿ, multiply by r and subtract: Sₙ(1 − r) = 1 − rⁿ⁺¹, so Sₙ = (1 − rⁿ⁺¹)/(1 − r). If |r| < 1, rⁿ⁺¹ → 0 and the sum is **1/(1 − r)**; the tail after n + 1 terms is exactly rⁿ⁺¹/(1 − r), which is the first remainder bound of the lesson and the model for all the others. If |r| ≥ 1 the terms do not go to zero and the series cannot converge (the **divergence test**: terms → 0 is necessary, not sufficient).

**Harmonic series.** 1 + ½ + ⅓ + … has terms that go to zero and still diverges. The **integral test** (`@openstax-calculus-volume-2` §5.3) compares the sum with ∫₁ⁿ dx/x = ln n, which P01.04 showed grows without bound: Hₙ ≈ ln n + 0.5772. H₁₀ = 2.93, H₁₀₀₀ = 7.49, H₁₀₀₀₀₀ = 12.09, H₁₀₀₀₀₀₀ = 14.39: each factor of ten adds 2.30, forever. The lesson's test asks for exactly that. By the same comparison Σ1/k² converges (∫dx/x² is finite); the boundary is p = 1 in Σ1/kᵖ, and the **comparison tests** (§5.4) let you decide most series by holding them against a geometric or a p-series.

### 2.2 The ratio test and the radius of convergence

For a series with terms aₖ, let ρ = lim |aₖ₊₁/aₖ|. If ρ < 1 the series converges absolutely (its tail is eventually dominated by a geometric series with ratio just above ρ); if ρ > 1 it diverges; if ρ = 1 the test says nothing (the harmonic series and Σ1/k² both have ρ = 1). For a **power series** Σcₖxᵏ the ratio is |cₖ₊₁/cₖ| |x|, so convergence holds for |x| < R = lim |cₖ/cₖ₊₁|, the **radius of convergence** (§6.1). Examples the code checks: 1/(1 − 2x) = Σ2ᵏxᵏ has R = ½; 1/(1 + x²) = Σ(−1)ᵏx²ᵏ has R = 1, although the function is perfectly smooth on the whole real line (the trouble is at x = ±i, which a real plot never shows); eˣ = Σxᵏ/k! has cₖ/cₖ₊₁ = k + 1 → ∞, so R = ∞. The Prandtl–Glauert series (§2.6) has coefficient ratio (2k + 1)/(2k + 2) → 1, so R = 1 in M: the series knows about the speed of sound.

Inside the radius a power series can be differentiated and integrated term by term (§6.2), which is how the series for arctan x = x − x³/3 + x⁵/5 − … comes from integrating 1/(1 + x²), and how the pendulum series comes from integrating a binomial series (Problem 7).

### 2.3 Taylor's theorem and the Lagrange remainder

The degree-n **Taylor polynomial** of f about x₀ is Tₙ(x) = Σₖ₌₀ⁿ f⁽ᵏ⁾(x₀)(x − x₀)ᵏ/k!: the unique polynomial of degree n with the same value and first n derivatives as f at x₀. **Taylor's theorem with the Lagrange remainder** (§6.3): if f⁽ⁿ⁺¹⁾ exists between x₀ and x,

f(x) = Tₙ(x) + Rₙ(x),  Rₙ(x) = f⁽ⁿ⁺¹⁾(ξ)(x − x₀)ⁿ⁺¹/(n + 1)!  for some ξ between x₀ and x,

so that **|Rₙ(x)| ≤ M |x − x₀|ⁿ⁺¹/(n + 1)!** whenever |f⁽ⁿ⁺¹⁾| ≤ M on that interval. For sin and cos every derivative is bounded by 1, so M = 1 always; for eˣ on [0, x], M = eˣ; for other functions you must bound the derivative yourself, which is the work. When the series alternates with terms decreasing in size, the **alternating-series bound** (§5.5) is simpler and usually tighter: |S − Sₖ| ≤ |aₖ₊₁|, the first omitted term.

The Maclaurin series (x₀ = 0) you must know (§6.4): sin x = x − x³/3! + x⁵/5! − …, cos x = 1 − x²/2! + x⁴/4! − …, eˣ = Σxᵏ/k!, 1/(1 − x) = Σxᵏ (|x| < 1), ln(1 + x) = x − x²/2 + x³/3 − … (−1 < x ≤ 1), and the **binomial series** (1 + x)^α = Σ C(α, k) xᵏ with C(α, k) = α(α − 1)…(α − k + 1)/k! (|x| < 1; terminates when α is a non-negative integer). Euler's formula e^(iθ) = cos θ + i sin θ is the exp series with x = iθ: the even powers of i are ±1 and give the cos series, the odd ones ±i and give the sin series. The code checks it term by term.

### 2.4 The small-angle approximations, with their error laws

Each approximation is a Taylor polynomial, and its relative error at small θ is the first omitted term divided by the function:

| approximation | polynomial | first omitted term | relative error ≈ | 1 % at | 0.1 % at |
|---|---|---|---|---|---|
| sin θ ≈ θ | T₁ = T₂ | −θ³/6 | θ²/6 | 14.0° | 4.44° |
| cos θ ≈ 1 | T₀ = T₁ | −θ²/2 | θ²/2 | 8.07° | 2.56° |
| cos θ ≈ 1 − θ²/2 | T₂ = T₃ | +θ⁴/24 | θ⁴/24 | 37.9° | 22.2° |
| tan θ ≈ θ | T₁ = T₂ | +θ³/3 | θ²/3 | 9.91° | 3.14° |

The 1 % angles come from inverting the exact error numerically (`max_angle_for_error`); the leading-term estimate √(0.06) = 14.03° for sin is within 0.3 % of the exact 13.99°. Two lessons sit in the table. First, the next *omitted* term is what matters, not the degree: sin θ ≈ θ is a degree-1 polynomial but its error is cubic because the θ² coefficient happens to be zero, which is why the Lagrange bound for sin should be written with n = 2 (bound θ³/6) and not n = 1 (bound θ²/2, which is true but four times too loose at 10°). Second, cos θ ≈ 1 is by far the worst of the four, and it is the one hidden in "the lift balances the weight" (W cos γ ≈ W) and in the load factor n = 1/cos φ in a turn; at 30° of bank the error of cos φ ≈ 1 is 15 %, and at 60° the "small-angle" load factor 1/(1 − φ²/2) = 2.21 instead of the true 2.00.

For the remainder bound at a definite angle, say sin 30°: θ = 0.5236, the degree-4 Lagrange bound (next derivative bounded by 1) is θ⁵/120 = 3.28 × 10⁻⁴, and the actual error of θ − θ³/6 is 3.26 × 10⁻⁴: the bound is honest to 1 %, because sin⁽⁵⁾ = cos ξ ≈ 1 near zero.

### 2.5 How many terms, and when the series must not be used

The sin series alternates with decreasing terms once θ² < (2k + 2)(2k + 3), so the first omitted term bounds the error. At θ = 0.5 the series needs 4 non-zero terms for 10⁻⁶ and 5 for 10⁻¹⁰; at θ = 2, 7 and 9; at θ = 10, 22 for 10⁻¹⁰ (the terms first grow to 10⁴/9! ≈ 2.8 × 10⁴ before they shrink). The series converges for every θ, so in exact arithmetic any argument is fine. In float64 it is not: at θ = 40 the largest term is 40¹⁹/19! = 1.5 × 10¹⁶ while the sum is 0.745, and each term carries a rounding error of about 2 × 10⁻¹⁶ of itself, so the sum carries an error of order 3: the code returns 2.0004 for sin 40. Likewise the exp series at x = −40 returns −8.0 instead of 4.2 × 10⁻¹⁸. The cure is not more terms (the tests check that 60 and 120 do not help) but **argument reduction**: sin 40 = sin(40 − 12π) = sin 2.3009, which 12 terms give to 10⁻¹⁵; and e⁻⁴⁰ = 1/e⁴⁰, where the series for e⁴⁰ has no cancellation and its reciprocal is correct to 16 digits. A series that converges is a mathematical object; evaluating it is a numerical decision.

### 2.6 Prandtl–Glauert as a binomial series

The subsonic compressibility correction multiplies incompressible pressure coefficients and lift-curve slopes by β⁻¹ = 1/√(1 − M²) = (1 − M²)^(−1/2). The binomial series with α = −½ and x = −M²:

1/√(1 − M²) = 1 + M²/2 + 3M⁴/8 + 5M⁶/16 + 35M⁸/128 + …,  coefficients C(−½, k)(−1)ᵏ = (2k)!/(4ᵏ k!²).

The coefficient ratio (2k + 1)/(2k + 2) tends to 1, so R = 1: the series converges for M < 1 and diverges for M ≥ 1, where the factor itself is imaginary or infinite. Inside the disc the convergence slows as M → 1 because the terms fall like M²ᵏ/√k: two terms are 0.3 % low at M = 0.3 and 5.6 % low at M = 0.6; ten terms are 0.3 % low at M = 0.8 and 4.3 % low at M = 0.9, where 16 terms are needed for 1 %. The "1 + M²/2" rule of thumb is a 1 % rule only up to M ≈ 0.42; the full factor costs nothing to compute and should be used, but the series is what tells you how the correction *scales* (as M², which is why compressibility is invisible below M = 0.3 and dominant above 0.7).

### 2.7 The pendulum swing test

The period of a pendulum released from amplitude θ₀ is T = T₀ (2/π) K(k) with k = sin(θ₀/2) and K the complete elliptic integral; expanding its integrand (1 − k² sin²φ)^(−1/2) by the binomial series and integrating term by term gives

T/T₀ = 1 + k²/4 + 9k⁴/64 + 25k⁶/256 + …,  coefficients [(2n)!/(4ⁿ n!²)]²,

or, re-expanded in θ₀, T/T₀ = 1 + θ₀²/16 + 11θ₀⁴/3072 + …. The inertia swing test of P01.04 measures I from T² and uses T₀; the period error at 10° of amplitude is 0.19 %, at 15° 0.43 %, at 22.8° 1.0 % and at 45° 4.0 %, and the error in I is twice that. The series with two terms is within 10⁻⁴ of the exact ratio up to 20° and within 0.4 % at 45°; at 90° (where k² = ½) it needs five terms for 0.3 % and the exact ratio is 1.1803. A swing test should be done below 10° of amplitude, or the measured period corrected by the series; either way the series gives the number.

### 2.8 Where this stops being valid

- **Outside the radius.** The Prandtl–Glauert series at M ≥ 1, ln(1 − f) at f ≥ 1 (a rocket cannot burn more than its mass), 1/(1 + x²) at |x| > 1: the terms grow and the partial sums are noise. The ratio test tells you the boundary; no number of terms crosses it.
- **At ρ = 1.** The ratio test is silent; the harmonic series looks convergent to any finite computation. Use the integral or comparison test.
- **Non-alternating or non-decreasing terms.** The alternating bound needs both; for sin at θ = 10 the first terms grow and the bound applies only from the ninth term on. `terms_needed` checks this.
- **Floating point at large arguments.** Cancellation (§2.5). Reduce the argument first.
- **The physics behind the approximation.** sin α ≈ α is the least of the linearisation's sins; the lift curve itself is only linear to 10–12° (P01.01), and that limit arrives first.

## 3. Build the computational model

Open `starter/main.py`.

### Step 1: series and their partial sums

```python
import numpy as np

def geometric_terms(r, n):
    return float(r) ** np.arange(n, dtype=float)

def partial_sums(terms):
    return np.cumsum(np.asarray(terms, dtype=float))

def series_terms(x, kind, n_terms):
    x = np.asarray(x, dtype=float)
    out = np.empty((n_terms,) + x.shape)
    if kind == "sin":
        term = x.copy()
        for k in range(n_terms):
            out[k] = term
            term = -term * x * x / ((2 * k + 2) * (2 * k + 3))   # ratio of consecutive terms
    ...
```

Each term is the previous one times a ratio; this never forms 40¹⁹ or 19! separately (both would overflow long before the ratio does) and it makes the ratio test visible in the code: the series converges because that ratio goes to zero. `maclaurin_partial_sum` sums the terms along axis 0.

### Step 2: Taylor polynomials with SymPy

```python
def taylor_polynomial(expr, x0, n):
    c = _centre(x0)                      # recognises pi/3, pi/4 ... in a float
    poly, d = sp.Integer(0), expr
    for k in range(n + 1):
        poly += sp.simplify(d.subs(x_s, c)) * (x_s - c) ** k / sp.factorial(k)
        d = sp.diff(d, x_s)
    return poly
```

The loop *is* Taylor's formula; `sp.series` would give the same result at x₀ = 0 and the test compares them. `taylor_coefficients` returns the numbers f⁽ᵏ⁾(x₀)/k!, which `radius_of_convergence` turns into R by the ratio of the last two non-zero coefficients (with the exponent 1/(k − j) when the non-zero coefficients are spaced, as in a series in x²).

### Step 3: remainders and terms needed

```python
def lagrange_remainder_bound(m_next, x, x0, n):
    return m_next * abs(x - x0) ** (n + 1) / factorial(n + 1)

def terms_needed(x, kind, tol, max_terms=200):
    t = series_terms(x, kind, max_terms + 1)
    for k in range(1, max_terms + 1):
        nxt = abs(t[k])
        if kind == "exp":
            nxt = np.exp(abs(x)) * abs(x) ** k / factorial(k)     # Lagrange, M = e^|x|
        elif abs(t[k - 1]) < nxt:
            continue                                               # terms still growing: no alternating bound yet
        if nxt < tol:
            return k
```

The `continue` is the part learners forget: the alternating bound is a theorem with hypotheses, and until the terms decrease it says nothing.

### Step 4: the small-angle explorer

```python
def small_angle_error(theta, approx):
    th = np.asarray(theta, dtype=float)
    exact, appr = {"sin": (np.sin(th), th), "cos": (np.cos(th), np.ones_like(th)),
                   "cos2": (np.cos(th), 1 - th * th / 2), "tan": (np.tan(th), th)}[approx]
    return np.abs(appr - exact) / np.abs(exact)

def max_angle_for_error(approx, tol):
    lo, hi = 0.0, 0.5 * np.pi - 1e-9
    for _ in range(100):                     # bisection on a monotone error
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if small_angle_error(mid, approx) < tol else (lo, mid)
    return 0.5 * (lo + hi)
```

`small_angle_error_leading` returns θ²/6, θ²/2, θ⁴/24, θ²/3 so that the plot can overlay the Taylor prediction on the measured curve. `prandtl_glauert_series` and `pendulum_period_ratio` are both built on `binomial_coefficients(-0.5, n)`, the first with x = −M² and the second with the coefficients squared.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-05-series-and-taylor`. Thirteen tests, grouped as in the metadata table. Two that fail for a wrong reason are worth knowing: `test_maclaurin_partial_sums…` *requires* the 60-term sin series at 40 to be wrong by more than 10⁻³ (an implementation that secretly reduces the argument is doing something the lesson has not asked for), and `test_lagrange_remainder_bound…` checks that `terms_needed` returns the smallest count the bound guarantees, so returning "a safe large number" fails.

### Compare

`sp.series(sp.sin(x), x, 0, 8).removeO()` is the Taylor polynomial of Step 2 in one call, and `sp.series` also handles centres and compositions (`sp.series(sp.sin(sp.exp(x)), x, 0, 5)`) that would be tedious by the loop. `np.sin(40.0)` is right to the last bit because the library reduces the argument modulo π/4 with an extended-precision constant before using a polynomial of about ten terms, which is Step 1 plus the extension. `scipy.special.ellipk(k**2)` gives K(k) directly, and `(2/np.pi) * ellipk(np.sin(theta0/2)**2)` reproduces `pendulum_period_ratio_exact` to 10⁻¹⁵; the lesson's midpoint rule on the elliptic integrand is there so that you see the series and the integral agree, not because SciPy lacks the function.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Prediction first, then the run.

1. **The error explorer.** Predict, from the table in §2.4, the ordering of the four curves and their slopes on a log-log plot (2, 2, 4, 2). Run `main()` and read off the 1 % and 0.1 % angles; compare with the leading-term estimates √(6 × tol), √(2 × tol), (24 × tol)^(1/4), √(3 × tol) in degrees. Record both columns; the leading-term estimate is within 1 % of the truth for sin and tan at 1 %, and 5 % off for cos₂ (37.9° against 40.1°) because θ⁴/24 is not small enough there to be the whole story.

| approximation | tol | leading-term angle | exact angle (code) |
|---|---|---|---|
| sin θ ≈ θ | 10⁻² | 14.03° | 13.99° |
| sin θ ≈ θ | 10⁻³ | 4.44° | 4.44° |
| cos θ ≈ 1 | 10⁻² | 8.10° | 8.07° |
| cos θ ≈ 1 − θ²/2 | 10⁻² | 40.1° | 37.9° |
| cos θ ≈ 1 − θ²/2 | 10⁻³ | 22.5° | 22.2° |
| tan θ ≈ θ | 10⁻² | 9.92° | 9.91° |

2. **Terms needed.** For θ = 0.25, 0.5, 1, 2, 4, 8 tabulate `terms_needed(θ, "sin", 1e-10)` and the size of the largest term. Predict: the count grows roughly linearly in θ once θ > 2 (the terms peak near k ≈ θ/2 and then need about e·θ/2 more to fall to the tolerance).

3. **Cancellation.** For θ = 5, 10, 20, 40 compute the sin series with 80 terms and the error against `np.sin`; alongside, the largest term times 2.2 × 10⁻¹⁶. The two columns track each other: the error of the sum is the rounding error of its largest term. Then reduce the argument by the nearest multiple of 2π and repeat.

Raw data format: three CSV tables (`approx, tol, angle_deg_exact, angle_deg_leading`; `theta, terms, max_term`; `theta, error, eps_times_max_term, error_reduced`).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| 1 % angle, sin θ ≈ θ | √0.06 = 14.03° | 13.99° | 0.3 % | next term θ⁴/120 in the ratio |
| 1 % angle, cos θ ≈ 1 − θ²/2 | (0.24)^(1/4) = 40.1° | 37.9° | 5 % | θ⁶/720 and the shrinking cos in the denominator |
| sin 30° by θ − θ³/6 | error ≤ 3.28 × 10⁻⁴ (Lagrange) | 3.26 × 10⁻⁴ | 1 % | sin⁽⁵⁾ = cos ξ < 1 |
| sin 40, 60 terms | ≈ 2 × 10⁻¹⁶ × 1.5 × 10¹⁶ ≈ 3 (rounding) | off by 1.26 | — | same order; the sign and size of the rounding errors are not predictable |
| 1/√(1 − 0.6²) by 2 terms | 1.18 against 1.25: 5.6 % low | 5.6 % | 0 | — |
| T/T₀ at 45° | 1 + θ²/16 + 11θ⁴/3072 = 1.03992 | 1.03997 (exact) | 5 × 10⁻⁵ | θ⁶ term |

What the agreement does not show: that sin α ≈ α is a good model of an aircraft at α. The linearisation of the equations of motion is a Taylor expansion of the whole dynamics about a trim point, and its validity is set by the first nonlinearity to bite, which is the lift curve near stall (P01.01), not the sine. The series here tell you the *mathematical* error of each replacement; the physical error is another lesson. What the agreement does show is that every "≈" has a number attached, and that the number was available before the exact answer was.

Improvements worth making: argument reduction for sin and cos; a `terms_needed` for the Prandtl–Glauert series with the geometric tail bound (terms fall like M²ᵏ, so the tail after k terms is below the next term divided by 1 − M²). Not worth making: summing the raw sin series at θ = 40 in higher precision; the right fix is cheaper and exact.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of the four small-angle error laws from Taylor's theorem, with the Lagrange bound for one of them at a stated angle; `small_angle_error.png` with its table; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5 and 7) on the cancellation experiment: the mechanism, the number, the fix.

### Adversarial self-checks

- **What if the angle doubles?** The sin and tan errors quadruple, the cos₂ error grows sixteen-fold; the 1 % angles in §2.4 scale as √tol or tol^(1/4). What if the tolerance is 10⁻⁴? sin: 1.40°, cos: 0.81°, cos₂: 12.6°, tan: 0.99°.
- **What if θ → 0 or θ → π/2?** All errors → 0 as their leading power; at π/2 the relative error of sin θ ≈ θ is 57 %, of cos θ ≈ 1 it is infinite (cos → 0), which is why `max_angle_for_error("cos2", 0.5)` stops at 76° and not at 90°.
- **Which assumption fails first in the real world?** For the linearised equations, the lift curve's linearity (10–12°), before the sine's 14°. For Prandtl–Glauert, the theory itself (it is a linearised small-perturbation theory) fails above the critical Mach number, typically 0.7–0.8, before its series does at 1.
- **What measurement could falsify the model?** The period of a real pendulum at 45° of amplitude against the series: 1.040 T₀. A real pendulum also has air drag and a non-rigid suspension; the series predicts the amplitude dependence, the measurement adds the rest.
- **Does the result survive a change in method?** The 1 % angles from bisection, from the leading-term formula and from reading the plot agree to 0.3 % (sin, tan, cos) and 5 % (cos₂); `pendulum_period_ratio_exact` by the midpoint rule, the series with 40 terms and SciPy's `ellipk` agree to 10⁻⁹.
- **Which decision would change if the uncertainty doubled?** If the acceptable error of the swing test were 2 % in I instead of 1 %, the amplitude limit would rise from 16° to 23°; a test that cannot keep the amplitude below 16° should correct its period by the series rather than reduce the swing further and lose timing accuracy.

## Check yourself

1. State the geometric series sum and its tail bound after n + 1 terms. Why is the tail bound exact for this series and only an inequality for others?
2. The terms of the harmonic series go to zero. Why is that not enough for convergence, and which test shows it diverges?
3. State the ratio test. What does it say about Σxᵏ/k!, and what does that mean for the exp series at x = −40 in exact arithmetic?
4. Write Taylor's theorem with the Lagrange remainder. What is M for sin, and why is the bound for sin θ ≈ θ best written with n = 2 rather than n = 1?
5. Derive the leading relative error of tan θ ≈ θ from the series tan θ = θ + θ³/3 + ….
6. Why does cos θ ≈ 1 fail at a smaller angle than sin θ ≈ θ? Name one flight-mechanics formula in which each is hidden.
7. The binomial series for 1/√(1 − M²) has radius 1. What happens to the partial sums at M = 1.2, and what does the radius mean physically?
8. The 60-term sin series at θ = 40 returns 2.0004. The series converges. Explain the number, and the fix.
