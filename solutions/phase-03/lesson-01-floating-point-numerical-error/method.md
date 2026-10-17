# P03.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r in metres; find e = ⌊log₂r⌋ (it is 22 for every seed). 2. ulp₆₄ = 2^(e − 52), ulp₃₂ = 2^(e − 23). 3. digits = −log₁₀(ulp/r) for each. 4. Compare each ulp with 0.01 m: the subtraction error is the ulp of the stored absolute positions. 5. The fix stores positions relative to a nearby reference so that the ulp at the *stored* magnitude is below the requirement.

## 2.
(a) fl(x̂ − ŷ) = (x(1 + δ₁) − y(1 + δ₂))(1 + δ₃); subtract x − y, divide by |x − y|, bound each δ by u. (b) Put x = 1 (exact), y = cos θ with absolute error ≤ u, x − y = θ²/2 + O(θ⁴); for the stable form, sin(θ/2) has relative error ≤ u, squaring doubles it, the factor 2 is exact. (c) E[x²] = m² + σ², each of E[x²] and m² is computed with relative error ~ε, so absolute error ~ε m²; divide by σ². The two-pass form subtracts m from each x (exact to u relative on numbers of size σ after the shift), then squares and sums positives. (d) Taylor: f(x + h) = f + h f′ + h² f″(ξ)/2. Rounding: f(x), f(x + h) each ≤ u F₀ absolute; x + h stored with error ≤ u|x|, which moves f by ≤ u|x f′|. Collect, divide by h, bound, differentiate E(h), solve, substitute back.

## 3.
1. Binade of 604 800 (between 2¹⁹ and 2²⁰) and of 3600 (between 2¹¹ and 2¹²). 2. ulp₃₂ = 2^(e − 23) in each; multiply by c. 3. ulp₆₄ = 2^(e − 52) at the end of the week; multiply by c. 4. Compare each with 10 ns × c = 3 m. 5. Represent time as an integer count of milliseconds (or whole seconds) plus a single-precision fraction, so the float only ever holds a number below 1.

## 4.
(a) `x = np.full(n, 0.1, np.float32)`; exact reference n × float(np.float32(0.1)); run the three sums with dtype float32; take absolute differences. (b) `observed_order(ns, errs_naive)`; explain through the grid spacing at the running sum S: spacing ∝ S ∝ n, bias per step ∝ spacing, n steps. (c) n × 0.1 − n × float32(0.1) in double; Kahan compensates the additions, not the stored value.

## 5.
(a) Build the data in double; call the four variance routines (two formulas × two dtypes); `correct_digits` against `numpy.var`. (b) Compare `numpy.spacing(numpy.float32(4e5))` with σ; the readings quantise to a grid coarser than their scatter. (c) Subtract 400 000 before casting; repeat (a).

## 6.
For each item: magnitude → ulp in binary32 → compare with the requirement; then classify the fix (none needed / formula / representation). The rule concerns what a float32 variable may hold (a quantity whose required absolute precision exceeds ~10⁻⁷ of its magnitude) and what must be stored as integer-plus-fraction or relative to a reference.

## 7.
(a) Assume x ≥ y > 0. Show 0 ≤ x − y ≤ y. Write x = X·ulp(y), y = Y·ulp(y) with integers X, Y (x's ulp is a multiple of y's because x ≥ y); then x − y = (X − Y)·ulp(y) with 0 ≤ X − Y ≤ Y < 2⁵³, hence representable. (b) With |y| ≤ |s|, t = fl(s + y) lies between s/2 and 2s (check both signs of y), so t − s is exact by (a); then c = fl((t − s) − y) is one rounding of the exact rounding error. (c) The compensated error per step is bounded by u times the step's rounding error, i.e. u² scale, plus the final rounding: summing n of these gives O(n u²) + O(u) Σ|xᵢ|, and n u² ≪ u for any n you can afford.
