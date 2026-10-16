# P03.01  Floating Point & Numerical Error

> A double carries sixteen digits. Every formula you write decides how many of them survive, and the ones that die are lost in subtraction, not in multiplication.

| | |
|---|---|
| **Learning outcomes** | 1. Decode an IEEE 754 double into sign, exponent and fraction, state machine epsilon and the unit roundoff, and compute the spacing (ulp) of doubles and singles at any magnitude. 2. Derive the rounding model fl(x op y) = (x op y)(1 + δ), |δ| ≤ u, verify it with exact rational arithmetic, and use it to predict the relative error of a cancelling subtraction and the n-dependence of a long sum. 3. Implement and test, from scratch, a cancelling and a stable form of 1 − cos θ, the quadratic formula, e^x by series, summation (naive, pairwise, Kahan) and the variance, and measure the digits each loses. 4. Derive the truncation-plus-rounding error model of a forward difference, predict the optimal step h* = 2√(ε f/f″), confirm the observed orders +1 and −1 on either side of it, and decide from the relative condition number whether a loss of digits is the problem's fault or the algorithm's. |
| **Prerequisites** | `P00.03` (NumPy arrays, dtypes), `P00.07` (Taylor series), `P01.05` (series, small-angle approximations and their remainders), `P02.05` (the Clohessy–Wiltshire coefficient 2(1 − cos nt)/n used as the running example). Skill check: you can write 0.1 in binary to six places and expand cos θ to the θ⁴ term. |
| **Engineering question** | A relative-navigation filter recovers the 1 km separation between a chaser and a station by subtracting two Earth-centred position vectors of magnitude 6.778 × 10⁶ m; a timing loop accumulates 0.1 s ticks for 100 hours in a 24-bit fixed-point register; and a rendezvous propagator evaluates 2(1 − cos nt)/n for a 1 s step at n = 1.13 × 10⁻³ rad/s. In single and double precision, how many correct digits does each computation deliver, and which formula should the propagator use? |
| **Assumptions and limits** | IEEE 754 binary64 and binary32 with round-to-nearest-even; no extended-precision intermediates (NumPy scalars and arrays do not use x87 80-bit registers); `math.sin`, `math.cos`, `math.exp` correctly rounded to within one ulp (true of every libm this course runs on, not guaranteed by the standard); the error models are worst-case bounds, and observed errors are typically 3 to 10 times smaller; nothing here covers decimal floating point, fixed point beyond one demonstration, or interval arithmetic. |
| **Analytical work** | The binary64 layout and the value (1 + f/2⁵²)·2^(e − 1023); ε = 2⁻⁵², u = ε/2, ulp(x) = 2^(⌊log₂x⌋ − 52); the rounding model and the propagation of relative error through +, −, ×, ÷; the relative error ε/(θ²/2) of 1 − cos θ computed naively; the Vieta rearrangement of the quadratic formula; the size of the largest term of the alternating series for e^(−x); the n u Σ|xᵢ| bound for naive summation and the 2u Σ|xᵢ| bound for compensated summation; the condition number ε m²/σ² of the one-pass variance; the forward-difference error h f″/2 + 2ε f/h and its minimiser. |
| **Computational deliverable** | `starter/main.py`: `float_components`, `machine_epsilon`, `ulp`, `relative_error`, `correct_digits`, `exact_rounding_error`, `fixed_point_truncation_error`, `one_minus_cos_naive`/`_stable`, `cw_position_coefficient`, `quadratic_roots_naive`/`_stable`, `exp_series`, `exp_stable`, `subtract_positions`, `naive_sum`, `kahan_sum`, `pairwise_sum`, `variance_one_pass`/`_two_pass`, `forward_difference`, `fd_total_error_bound`, `fd_optimal_step`, `relative_condition_number`, `observed_order`. Tests prove: ε and ulp agree with `numpy.finfo` and `numpy.spacing`; the decoded fields rebuild the double; every one of 800 random +, −, ×, ÷ results has |δ| ≤ u against exact `Fraction` arithmetic; the stable 1 − cos θ holds 4ε while the naive one follows ε/(θ²/2) and is exactly 0 at θ = 10⁻⁸; the stable quadratic roots satisfy Vieta's relations to 10⁻¹⁵ while the textbook small root is 25 % off; the series e^(−20) has no correct digits and e^(−30) is wrong by 10⁷; subtracting absolute positions costs the ulp of 6.778 × 10⁶ m, 0.5 m in single precision; the three sums agree with `math.fsum` to 10⁻¹² in double and in single Kahan ≤ pairwise ≪ naive with the naive error growing with n; the two-pass variance matches `numpy.var` while the one-pass form loses every digit for a 400 km mean; the forward difference has observed order 1.00 ± 0.05 in the truncation regime and −1.0 ± 0.4 in the rounding regime, stays under the error model everywhere, and reaches its minimum within a decade of h*. |
| **Practical activity** | **computational experiment**: measure digits lost against θ for 1 − cos θ, against n for three summation algorithms in single precision, against the mean-to-sigma ratio for the one-pass variance, and against h for the forward difference, each with the prediction written first. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `floating_point_error.png`, the four experiment tables, the hand derivations of the rounding model consequences (2.3), the quadratic rearrangement (2.4), the forward-difference optimum (2.6), and a half-page note stating for each of the three engineering-question computations the digits delivered in single and double precision and the formula or representation you recommend. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the ε/(θ²/2) law, the Vieta rearrangement and h* are derived; all public tests pass; each experiment table has a prediction column and an observed column with the discrepancy explained; the note names the representation change for the single-precision navigation case. |
| **References** | `@ocw-18-335j` Week 1 "Lecture 2: Floating-Point Arithmetic" (handout "Floating-Point Arithmetic"), Week 2 "Lecture 3: Floating-Point Summation and Backwards Stability", "Lecture 5: Condition Numbers" (required-free): formats, ε, backward error, conditioning. `@goldberg-1991` "Rounding Error" > "Relative Error and Ulps", "Cancellation"; "The IEEE Standard" > "Special Quantities" (required-free): the rounding model, benign versus catastrophic cancellation, inf/nan/subnormals. `@ocw-18-330` "Lecture notes 1–2: Series and Sequences" 1.2 The big-O notation (supplementary-free): the O(h) notation used in 2.6. |
| **Extensions** | ★ Implement Welford's single-pass variance update and show it matches the two-pass result to ε even for the 400 km data. ★ Implement `fma`-free double-double addition (two doubles carrying a 106-bit sum) and recompute the 10⁶-term sum. ★ Replace the forward difference by the complex-step derivative Im f(x + ih)/h and show its error has no rounding regime at all. |
| **Time** | unmeasured |

## 1. The engineering problem

Three computations from earlier lessons are about to run on a flight computer, and each hides a subtraction of nearly equal numbers.

The first is relative navigation. A chaser and a station each know their Earth-centred position to a few centimetres. The station is at radius 6.778 × 10⁶ m; the chaser is 50 m below and 1 km behind it (`P02.05`). The guidance law needs their *separation*, and the obvious code stores both absolute positions and subtracts. In double precision the spacing of representable numbers at 6.778 × 10⁶ m is 9.3 × 10⁻¹⁰ m, so the subtraction costs a nanometre. In single precision, which many flight processors and most GPUs prefer, the spacing at that magnitude is **0.5 m**: the centimetre requirement is dead before the first filter update, and no amount of clever filtering afterwards brings the digits back.

The second is a clock. A real-time loop counts ticks of 0.1 s and multiplies the count by 0.1 to get elapsed time. 0.1 has no finite binary expansion, so a 24-bit fixed-point register stores it short by 3.6 × 10⁻⁸ s. Over 100 h that is 3.6 × 10⁶ ticks and 0.13 s of drift; at 1676 m/s, the closing speed of a ballistic target, 0.13 s is 216 m. A widely documented 1991 air-defence failure traced to exactly this mechanism, and the lesson of that incident is not "use more bits" but "know what the representation can and cannot hold".

The third is the Clohessy–Wiltshire coefficient 2(1 − cos nt)/n from `P02.05`. At n = 1.13 × 10⁻³ rad/s and a 1 s step the argument is 10⁻³ rad, cos nt = 0.99999949…, and 1 − cos nt subtracts two numbers that agree in their first six digits. Six of the sixteen digits of the result are gone, and a propagator that takes 10⁵ such steps per orbit accumulates whatever bias they leave. The identity 1 − cos θ = 2 sin²(θ/2) costs one more function call and loses nothing.

Success criteria, written first: every number in this lesson is predicted by an error model before it is computed, and the model must bound the observed error without being more than a factor 10 pessimistic. A learner who finishes can look at a formula, say which subtraction will cancel, estimate the digits it loses from the sizes of the operands alone, and rearrange it or change the representation so that the loss does not happen.

## 2. Understand and derive

### 2.1 What a double is

An IEEE 754 binary64 number is 64 bits: 1 sign bit s, an 11-bit exponent field e (0 to 2047) and a 52-bit fraction field f. For 0 < e < 2047 the value is

**x = (−1)^s · (1 + f/2⁵²) · 2^(e − 1023)**.   (1)

The "1 +" is the implicit leading bit: every normal number has a 53-bit significand with the first bit always 1, so it is not stored. e = 0 encodes zero and the **subnormal** numbers (no implicit 1, exponent −1022), which fill the gap between 0 and 2⁻¹⁰²² evenly; e = 2047 encodes ±inf (f = 0) and nan (f ≠ 0). binary32 has 8 exponent bits and 23 fraction bits, so a 24-bit significand.

`float_components(0.1)` returns sign 0, exponent 1019 (so 2⁻⁴) and fraction 0x999999999999a: 0.1 is stored as 3602879701896397/36028797018963968 = 0.1000000000000000055511151231257827…, the nearest double, 5.55 × 10⁻¹⁸ too large. Every decimal that is not a dyadic rational is stored with such an error, which is why 0.1 + 0.2 ≠ 0.3 and why a tick of 0.1 s drifts.

Two numbers describe the grid. **Machine epsilon** ε = 2⁻⁵² = 2.22 × 10⁻¹⁶ is the spacing of doubles between 1 and 2: the smallest power of two that changes 1.0 when added to it. The spacing at any x, called an **ulp** (unit in the last place), scales with the binade: **ulp(x) = 2^(⌊log₂|x|⌋ − 52)**, so ulp(6.778 × 10⁶) = 2^(22 − 52) = 9.3 × 10⁻¹⁰ and the relative spacing ulp(x)/x is always between ε/2 and ε. In binary32, ε = 2⁻²³ = 1.19 × 10⁻⁷ and ulp(6.778 × 10⁶) = 2^(22 − 23) = 0.5. Sixteen digits in double, seven in single, *relative* to the number held: that is the whole story of the navigation example.

### 2.2 The rounding model

Round-to-nearest places the result of any single operation within half an ulp of the exact result, so with the **unit roundoff** u = ε/2 = 1.11 × 10⁻¹⁶,

**fl(x op y) = (x op y)(1 + δ), |δ| ≤ u**, op ∈ {+, −, ×, ÷}.   (2)

This is the only fact about arithmetic the rest of the lesson uses. `exact_rounding_error` checks it directly: compute x op y exactly with `fractions.Fraction`, compare with the double result, and read δ. The tests do this for 800 random operations and never find |δ| > u; for 0.1 + 0.2 the model is also tight, δ = 9.25 × 10⁻¹⁷.

Three consequences. **Multiplication and division are safe**: a relative error of δ₁ in x and δ₂ in y gives xy(1 + δ₁)(1 + δ₂)(1 + δ₃) ≈ xy(1 + δ₁ + δ₂ + δ₃), so relative errors *add*, and a product of twenty rounded numbers is still good to 20u ≈ 2 × 10⁻¹⁵. **Addition of like signs is safe** for the same reason. **Subtraction of nearly equal numbers is not**: if x and y each carry relative error u and x − y is much smaller than x, the *absolute* errors u|x| and u|y| are unchanged while the result shrank, so the relative error of the difference is about u(|x| + |y|)/|x − y|. The subtraction itself is exact (Sterbenz: if y/2 ≤ x ≤ 2y, x − y is a double); what it does is *expose* the rounding errors already present in x and y. This is **catastrophic cancellation**, and it is where digits die.

### 2.3 Four cancellations, predicted

**1 − cos θ.** For small θ, cos θ = 1 − θ²/2 + θ⁴/24 − …, and the double nearest cos θ carries absolute error up to u. Subtracting from 1 leaves that absolute error on a result of size θ²/2, so

relative error of 1 − cos θ (naive) ≈ u/(θ²/2) to ε/(θ²/2).   (3)

At θ = 10⁻³: ε/(5 × 10⁻⁷) = 4.4 × 10⁻¹⁰, about six digits lost. At θ = 10⁻⁶: 4.4 × 10⁻⁴, twelve digits lost. At θ = 10⁻⁸, θ²/2 = 5 × 10⁻¹⁷ is below u, cos θ rounds to exactly 1.0, and the naive result is exactly 0: every digit lost. The identity 1 − cos θ = **2 sin²(θ/2)** evaluates sin to relative error u and squares and doubles it (safe operations), so its relative error is about 3u for every θ. The problem is well conditioned (2.5); only the formula is bad. For the Clohessy–Wiltshire coefficient at n t = 1.13 × 10⁻³ the two forms give 1.131399879407 × 10⁻³ and 1.131399879311 × 10⁻³ s: they disagree at the tenth digit, and the stable one agrees with the series n t²(1 − (nt)²/12 + …) to 10⁻¹⁴.

**The quadratic formula.** For a x² + b x + c = 0 with b² ≫ |4ac|, √(b² − 4ac) ≈ |b| and the root −b + √D (for b > 0) cancels. Rearrange: let q = −(b + sign(b)√D)/2, which adds like signs. Then x₁ = q/a is the large root and, by Vieta's product x₁x₂ = c/a, **x₂ = c/q** is the small one, obtained by a division instead of a subtraction. For x² + 10⁸x + 1 = 0 the small root is −1.0000000000000001 × 10⁻⁸; the textbook formula returns −7.45 × 10⁻⁹ (25 % error: D = 10¹⁶ − 4 is representable but √D = 10⁸ − 2 × 10⁻⁸ rounds to the nearest double at spacing 1.5 × 10⁻⁸) and the Vieta form returns −1.0 × 10⁻⁸ to 10⁻¹⁵. The aerospace version: a chaser closing at v = 1 m/s on a target 10 m away under a thrust acceleration a = 10⁻⁶ m/s² reaches it when ½at² + vt − 10 = 0; the textbook form gives 9.999950000505 s, the stable form 9.999950000500 s, and only the second has zero residual.

**e^x by series.** Σ xᵏ/k! converges for every x, but for x = −20 the terms alternate in sign and the largest, 20²⁰/20! = 4.3 × 10⁷, is 2 × 10¹⁶ times the answer e⁻²⁰ = 2.06 × 10⁻⁹. Each term carries rounding error of order u × (its size), so the sum carries absolute error ≈ ε × 4.3 × 10⁷ ≈ 10⁻⁸: larger than the answer. `exp_series(−20.0)` returns 6.15 × 10⁻⁹ (no correct digits) and `exp_series(−30.0)` returns 6.1 × 10⁻⁶ against 9.4 × 10⁻¹⁴. The fix is again algebraic: e⁻ˣ = 1/eˣ, and the series for e⁺ˣ has only positive terms, so its relative error is a few u.

**The one-pass variance.** σ² = E[x²] − (E[x])² is one subtraction of two numbers of size m² (m the mean) to get a result of size σ². The relative error is therefore about ε m²/σ²: the **condition number of the formula**. For altimeter readings with m = 400 km and σ = 1 cm, m²/σ² = 1.6 × 10¹⁵ and ε m²/σ² = 0.36: the one-pass variance has no correct digits in double (`main()` prints −2.44 × 10⁻⁴ for a true 1.01 × 10⁻⁴; a negative variance). In single precision it returns 2.7 × 10⁷. The two-pass formula mean((x − m)²) subtracts first, while the numbers are still exact to u, and then squares and adds (safe), giving relative error of a few ε. Note what single precision does to the *data*: ulp(4 × 10⁵) in binary32 is 0.03 m, more than σ, so the readings cannot even be stored, and no formula helps. That is a representation problem, not an algorithm problem, and the fix is to store x − 400 000 m.

### 2.4 Summation: how error grows with n

Adding n numbers left to right applies (2) n − 1 times. Each rounding error is at most u times the running partial sum, so the worst-case error is about **n u Σ|xᵢ|** (18.335J Lecture 3). Two things improve on it. **Pairwise summation** (add halves recursively, which is what `numpy.sum` does in blocks) applies only log₂n roundings to any one term, giving ≈ u log₂n Σ|xᵢ|. **Kahan compensated summation** computes the rounding error of each addition exactly (c = (t − s) − y, which is exact by Sterbenz) and feeds it into the next addition, giving ≈ 2u Σ|xᵢ| independent of n. In double these are invisible on 10⁶ numbers of order 1 (all three agree with `math.fsum`, the exactly rounded sum, to 10⁻¹²), so `main()` runs them in single: for 10⁶ uniform numbers the naive sum is off by 8.4 against a total of 499 897, pairwise by 0.015 and Kahan by 0.015, both at the half-ulp floor of the binary32 result itself (ulp(5 × 10⁵) = 0.03). For equal summands the naive error is systematic rather than random and grows like n² (every addition rounds the same way once the sum exceeds the summand by enough): Problem 4 measures that order.

### 2.5 Conditioning: the problem's fault or the algorithm's?

The **relative condition number** of a function at x is the factor by which a relative error in the input is amplified in the output:

**κ(x) = |x f′(x)/f(x)|**.   (4)

κ(√x) = ½, κ(eˣ) = |x|, κ(1 − cos x) → 2 as x → 0, κ(ln x) = 1/|ln x| → ∞ as x → 1. A problem with κ ≫ 1 loses log₁₀κ digits *whatever algorithm is used*, because the input itself is only known to u; a problem with κ ~ 1 that loses digits is being solved by an **unstable algorithm**, and a stable one exists. 1 − cos θ is the clean case: κ = 2, so the six lost digits in the naive formula are entirely the formula's fault, and 2 sin²(θ/2) recovers them. ln(1 + 10⁻⁶) is the other case: κ = 10⁶, six digits are lost by anyone who only knows 1 + 10⁻⁶ to u; the only cure is to keep the small quantity 10⁻⁶ itself (which is what `numpy.log1p` asks you to pass). An algorithm is **backward stable** if its output is the exact answer to a nearby question (relative perturbation of the input by a few u); backward stability plus a condition number gives the forward error bound κ · few u, and that factorisation is how every result in Phase 03 will be judged.

### 2.6 Truncation against rounding: the forward difference

Numerical methods have a second error source that has nothing to do with floating point: **truncation**, the error of replacing a limit by a finite step. The forward difference D(h) = (f(x + h) − f(x))/h estimates f′(x) with truncation error h f″(ξ)/2 (Taylor with remainder, `P01.05`). But f(x + h) and f(x) each carry rounding error up to u|f|, their difference carries up to ε|f|, and x + h itself is rounded with error up to u|x| that enters as u|x f′|. Divided by h:

**|D(h) − f′(x)| ≤ (h/2)·F₂ + 2ε F₀/h**, F₂ ≥ |f″|, F₀ ≥ |f|, |x f′|.   (5)

The first term falls with h, the second rises, and the sum is minimised where they are equal: **h* = 2√(ε F₀/F₂)**, with minimum error 2√(ε F₀ F₂). For f = sin at x = 1, h* = 2.98 × 10⁻⁸ and the best achievable error is 3 × 10⁻⁸: *half* the digits of the machine, which no choice of h can beat. `main()` confirms it: the observed error is 4.2 × 10⁻³ at h = 10⁻², falls by a factor 10 per decade (observed order 1.000) to 2.97 × 10⁻⁹ at h = 10⁻⁸, then rises again by about a factor 10 per decade (observed order −1.1 over h = 10⁻¹⁰ to 10⁻¹⁵) to 0.54 at h = 10⁻¹⁶, where x + h = x and D = 0. The bound (5) sits above every observed point and is within a factor 10 of it. Lesson `P03.04` builds better stencils; the shape of this curve, a V on log-log axes with its vertex at √ε, is what every difference quotient looks like and what Richardson extrapolation and the complex-step method (extension) are designed to beat.

**Observed order** is the tool for reading such a curve: fit log(error) against log(h) by least squares and read the slope p of the law error ≈ C hᵖ. It is +1 for a first-order method in its truncation regime, −1 where rounding dominates a quotient, and will be 2, 4 and higher for the methods of later lessons. Every lesson in this phase ends with a measurement of it.

### 2.7 Where the model stops being valid

- (2) assumes no overflow (|x| > 1.8 × 10³⁰⁸ gives inf) and no underflow into the subnormal range, where the relative spacing grows. An orbit propagator in SI units never approaches either, but a probability of 10⁻³²⁰ does.
- Library functions are not guaranteed correctly rounded by the standard; modern libms are within 1 ulp, and the models here allow a factor 2 for that.
- Fixed point is a different arithmetic: the error of a tick is absolute, not relative, and accumulates linearly with the count (2.1). The lesson treats it only through `fixed_point_truncation_error`.
- Worst-case bounds (n u for sums) are rarely attained with random data, where errors partially cancel and grow like √n; they *are* attained with systematic data (Problem 4). Design to the bound, expect the typical.

## 3. Build the computational model

Open `starter/main.py`. Every function body is yours; `main()` is given.

### Step 1: the representation

`float_components` packs the double big-endian with `struct.pack(">d", x)`, unpacks it as an unsigned 64-bit integer and slices the fields:

```python
(bits,) = struct.unpack(">Q", struct.pack(">d", float(x)))
sign, exponent, fraction = bits >> 63, (bits >> 52) & 0x7FF, bits & ((1 << 52) - 1)
```

`machine_epsilon(dtype)` starts at eps = 1 and halves while `one + eps/2 != one`, doing every operation in the requested dtype (cast with `dtype(...)` at each step, or NumPy promotes to double). `ulp(x)` is 2^(⌊log₂|x|⌋ − 52), clamped at exponent −1022 for subnormals, 2⁻¹⁰⁷⁴ at zero. `exact_rounding_error` converts both operands to `Fraction`, performs the exact operation and the double operation, and returns (fl − exact)/exact as a float. `fixed_point_truncation_error(value, bits)` is value − ⌊value·2^bits⌋/2^bits in `Fraction` arithmetic (use `limit_denominator` so that 0.1 means one tenth, not the double nearest it).

### Step 2: the cancelling pairs

Type in the formulas of 2.3 exactly as written: `1.0 - np.cos(theta)` against `2.0 * np.sin(theta / 2) ** 2`; the textbook (−b ± √D)/2a against q = −(b + copysign(√D, b))/2, roots q/a and c/q; the series with the recurrence term *= x/k stopped when |term| ≤ ε|sum| against 1/series(−x) for negative x. Both members of each pair return the same ordering and raise the same `ValueError` for a = 0 or D < 0, so the tests can compare them blind. `cw_position_coefficient(n, t, stable)` is 2 · (one of the two forms)(n t)/n. `subtract_positions` casts the station position and the station-plus-separation to the requested dtype and subtracts; the result's error is the ulp of the absolute position.

### Step 3: sums and variances in a chosen dtype

Loop explicitly (a NumPy reduction would silently do pairwise summation and hide the effect) and cast every intermediate to the dtype. Kahan's four lines are

```python
y = dtype(v - comp); t = dtype(total + y); comp = dtype(dtype(t - total) - y); total = t
```

`pairwise_sum` recurses on halves and adds the two results in the dtype. The two variance functions follow 2.3; both loop in the dtype.

### Step 4: the forward difference and its model

`forward_difference(f, x, h)` is one line. `fd_total_error_bound` is (5) and must accept an array of h; `fd_optimal_step` is h*. `relative_condition_number` is (4) and raises where f(x) = 0. `observed_order(hs, errs)` is `np.polyfit(log hs, log errs, 1)[0]` after validating that there are at least two positive pairs.

Running `main()` prints

```
machine epsilon (double) 2.220e-16 = 2^-52; (single) 1.192e-07 = 2^-23
0.1 stored as sign=0 exponent=1019 (2^-4) fraction=0x999999999999a; exactly 3602879701896397/36028797018963968
ulp(1.0) = 2.220e-16, ulp(6.778e6 m) = 9.313e-10 m, float32 ulp at 6.778e6 m = 0.500 m
separation recovered in float64: error 4.053e-10 m
separation recovered in float32: error 1.230e-02 m
```

followed by the 1 − cos θ table (naive relative error 1.1 × 10⁻¹⁴ at θ = 0.1, 1.6 × 10⁻¹¹ at 10⁻³, 8.9 × 10⁻⁵ at 10⁻⁶, exactly 1 at 10⁻⁸; stable 0 at every θ against the exact-rational reference), the quadratic, series, summation and variance lines of 2.3–2.4, the forward-difference summary (best h = 10⁻⁸, error 2.97 × 10⁻⁹, model h* = 2.98 × 10⁻⁸; orders 1.000 and −1.107), the clock drift (3.576 × 10⁻⁸ s per tick, 0.129 s in 100 h, 216 m at 1676 m/s), and the figure.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-01-floating-point-numerical-error`. The tests:

1. **Representation.** ε equals `numpy.finfo(...).eps` for both dtypes; the fields of 1.0, −2.0, 0.1 and 0.0 are right and (1) rebuilds 0.1 exactly; `ulp` equals `numpy.spacing` at six magnitudes from 3 × 10⁻³¹⁰ (subnormal) to 10³⁰⁰, and raises for inf.
2. **The rounding model.** 800 random operations in all four ops have |δ| ≤ u by exact rational comparison; 0.1 + 0.2 has the predicted nonzero δ; the fixed-point error of 0.1 is below one LSB at 8, 16, 24 and 32 bits and zero for 0.5.
3. **Cancellation pairs.** The stable 1 − cos θ holds 4ε at seven magnitudes while the naive one follows (3) within a factor 10, exceeds 10⁻⁴ at θ = 10⁻⁷ and is exactly 0 at 10⁻⁸; the Clohessy–Wiltshire coefficient matches its series only in stable form; Vieta's product and sum hold to 10⁻¹⁵ for the stable roots and the textbook small root of x² + 10⁸x + 1 is more than 10 % off; both root finders agree with `numpy.roots` for twenty random well-conditioned quadratics and raise on a = 0 and D < 0; the series e^x matches `math.exp` to 10⁻¹⁴ for x ≥ 0 and fails by more than 10 % at −20 and 10³ at −30.
4. **Representation of positions.** Subtraction error ≤ ulp(6.778 × 10⁶) in double and between 1 mm and 0.5 m in single.
5. **Sums and variances.** All three sums match `math.fsum` to 10⁻¹² in double; in single, Kahan ≤ 1.5 × pairwise < naive/5, naive > 10⁻³, and Kahan satisfies its 2u Σ|x| bound; the naive error grows more than tenfold from n = 2000 to 200 000 while Kahan's relative error stays below 10⁻⁶; the two-pass variance matches `numpy.var` to 10⁻¹⁰, the one-pass loses more than 10 % in double and a factor 10³ in single, and both agree on centred data.
6. **Convergence order.** The forward difference of sin at 1 has observed order 1.00 ± 0.05 over h = 10⁻², 10⁻³, 10⁻⁴ and −1.0 ± 0.4 over 10⁻¹⁰ to 10⁻¹⁵; the bound (5) exceeds every observed error; the best h lies within a decade of h*; `observed_order` recovers slopes 2 and −1 exactly from synthetic power laws.
7. **Conditioning.** κ = ½ for √x at 4, 3 for eˣ at 3, 2 for 1 − cos x near 0, > 10⁵ for ln x at 1 + 10⁻⁶, and a raise at f(x) = 0.

### Compare

The oracles here are not SciPy but the exact tools Python already ships: `fractions.Fraction` for exact rational arithmetic, `math.fsum` for the correctly rounded sum of a list, `numpy.spacing` and `numpy.finfo` for the grid. In practice you call `numpy.sum` (pairwise, in C, and fast), `numpy.var` (two-pass), `numpy.log1p`, `numpy.expm1` and `numpy.hypot` (each a stable form of a cancelling expression), and `numpy.roots` (an eigenvalue method that is stable for all coefficient sizes). What writing the naive versions bought is the ability to recognise, in code you did not write, the shape `big - big`, `1 - cos`, `exp(-large)` by series, `E[x²] - E[x]²`, `(f(x+h) - f(x))/h` with h = 10⁻¹², and to estimate the damage from the magnitudes alone before running anything.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Digits against θ.** For θ = 10⁻¹ … 10⁻⁸, tabulate the relative error of `one_minus_cos_naive` against the stable form. Predict the digits lost from (3): log₁₀(2/θ²) = 2 log₁₀(1/θ) + 0.3, i.e. 2.3 at 0.1, 6.3 at 10⁻³, 12.3 at 10⁻⁶. Record `theta, rel_err_naive, rel_err_stable, digits_lost_predicted, digits_lost_observed`.
2. **Summation against n in single precision.** For n = 10³, 10⁴, 10⁵, 10⁶ uniform random numbers (`default_rng(301)`), tabulate the absolute error of the three sums against `math.fsum`. Predict: naive error ∝ between √n and n, pairwise ∝ log n, Kahan flat at the half-ulp of the result. Then repeat with every summand equal to 0.1 and fit the observed order of the naive error against n (Problem 4 asks for the number and the reason).
3. **Variance against m/σ.** Fix σ = 1 and sweep the mean m = 1, 10², 10⁴, 10⁶, 10⁸ (10⁴ samples, seeded). Predict the relative error of the one-pass form as ε m²/σ² (a scale, not a bound: the sum of 10⁴ squares adds its own rounding) and record where it first exceeds 10⁻⁶, 10⁻², 1.
4. **The V-curve.** For f = exp at x = 2 (F₀ = F₂ = e²), tabulate the forward-difference error for h = 10⁻¹ … 10⁻¹⁶, fit the observed order on each side, and compare the best h with h* = 2√ε = 3 × 10⁻⁸. Then repeat for f = sin at x = 100 (F₀ must now cover |x f′| = 100) and predict how h* moves.

Raw data format: four CSV tables with a one-line header naming the dtype, seed and fixed parameters; `floating_point_error.png` with the two panels of `main()`.

What you should find: the digits-lost column sits one to two digits below the prediction of (3) (1.7, 4.8 and 11.6 observed), and the observed error is 2 to 30 times *below* the bound because actual rounding errors are smaller than u; in single precision the naive uniform sum error grows from about 10⁻² at n = 10⁴ to 8.4 at 10⁶ while Kahan and pairwise sit at 0.015; the one-pass variance error runs about 10 times above ε m²/σ² (the 10⁴ squared terms accumulate their own rounding, 2.4) and passes relative error 1 between m/σ = 10⁶ and 10⁸ in double; the V-curve has slopes +1.00 and about −1 with its vertex at h = 10⁻⁸ for exp at 2, but for sin at x = 100 the vertex moves to h = 10⁻⁷ (model h* = 3 × 10⁻⁷) because the argument's own rounding error u·100 now dominates F₀.

## 5. Validate and improve

| quantity | model (§2) | observed (code) | ratio | cause of the gap |
|---|---|---|---|---|
| ε, ulp(6.778 × 10⁶) double / single | 2⁻⁵², 9.31 × 10⁻¹⁰ m / 0.5 m | same | 1 | exact |
| max |δ| over 800 operations | ≤ 1.11 × 10⁻¹⁶ | ≤ 1.11 × 10⁻¹⁶ | ≤ 1 | bound attained within 10 % |
| 1 − cos θ naive, θ = 10⁻³ | ≤ 4.4 × 10⁻¹⁰ | 1.57 × 10⁻¹¹ | 0.04 | actual rounding of cos θ was far below u |
| 1 − cos θ naive, θ = 10⁻⁶ | ≤ 4.4 × 10⁻⁴ | 8.9 × 10⁻⁵ | 0.2 | same |
| CW coefficient, two forms | differ at 10⁻¹⁰ | differ by 8.5 × 10⁻¹¹ relative | — | six digits lost, as predicted |
| small root of x² + 10⁸x + 1, textbook | O(ulp(10⁸)/10⁻⁸) = O(1) | 25 % | — | √D rounded to a grid of 1.5 × 10⁻⁸ |
| e⁻²⁰ by series, absolute error | ≈ 10⁻⁸ | 4.1 × 10⁻⁹ | 0.4 | partial cancellation of term errors |
| naive float32 sum, n = 10⁶ | ≤ n u Σ|x| = 3 × 10⁴ (worst), ~√n u Σ|x| = 30 (typical) | 8.4 | 0.3 of typical | random rounding, partial cancellation |
| one-pass variance, 400 km mean | relative error ~0.36 | −3.4 (negative result) | — | the "0.36" is the scale of the noise; sign is luck |
| forward difference, best error | 3 × 10⁻⁸ at h* = 3 × 10⁻⁸ | 3 × 10⁻⁹ at h = 10⁻⁸ | 0.1 | (5) is a worst-case bound; both terms were below their bounds |
| observed orders | +1, −1 | +1.000, −1.11 | — | rounding regime is noisy; ±0.4 is the honest tolerance |

Every row has the model above the observation, and never by more than a factor 30: the error model is doing its job, which is to tell you *before* you run whether a computation can meet its requirement. The discrepancies all have the same cause, that (2) is a bound and most roundings are smaller, which is why engineering practice designs to the bound and reports the observed. The one row where the model is not a bound but a scale, the one-pass variance, is also the one where the sign of the answer came out wrong: when the predicted relative error is of order 1, the result is noise and should be reported as such, not as "−2.4 × 10⁻⁴".

What to change: nothing in the stable forms. In the navigation example, change the *representation*: store the chaser state relative to the station (or relative to a reference point updated each orbit), so that single precision holds centimetres on a 1 km quantity (ulp(1000) in binary32 is 6 × 10⁻⁵ m). In the clock, count ticks in an integer and convert once, or use a tick that is a power of two (1/16 s). In the propagator, use 2 sin²(θ/2), or the small-angle series when nt < 10⁻⁴, and test both against each other.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), the Vieta rearrangement, the ε m²/σ² condition number and h*; `floating_point_error.png`; the four experiment tables with prediction and observation columns; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that, for each of the three engineering-question computations, states the digits delivered in double and single precision and recommends a formula or a representation.

### Adversarial self-checks

- **What if the separation is 10 m instead of 1 km?** The absolute positions are unchanged, so the subtraction error is unchanged (0.5 m in single, 10⁻⁹ m in double) while the quantity shrank by 100: the relative error of the single-precision separation is now 5 %. Which representation would you choose if the mission's final approach ends at 10 m?
- **What if the step shrinks to zero?** The naive 1 − cos θ returns exactly 0 for θ < 1.5 × 10⁻⁸, the forward difference returns exactly 0 for h < ulp(x), and the one-pass variance may return a negative number. A result of exactly 0 from a formula that should be small is the signature of cancellation; add a test for it.
- **Which assumption fails first in the real world?** That the libm functions are within an ulp: on some embedded targets `cosf` is a polynomial good to 10⁻⁶, and the error model of 2.3 must use that error instead of u. Measure the library before trusting the model.
- **What measurement could falsify the model?** A single operation with |δ| > u found by `exact_rounding_error` would mean the arithmetic is not round-to-nearest (a flush-to-zero or directed-rounding mode has been set, which some compilers do for speed); a forward-difference error above (5) would mean f is not evaluated to within the assumed ulp.
- **Does the result survive a change of dtype, order or compiler?** The stable forms give the same digits in single and double up to each dtype's ε; the naive ones change their *values* when the summation order changes (`naive_sum` versus `pairwise_sum` versus `numpy.sum` differ in the last digits), which is why a reproducibility requirement in flight software fixes the order of every reduction.
- **Which decision changes if the requirement tightens by a factor 2?** A 5 mm separation requirement is still met in double by subtraction (10⁻⁹ m) and still missed in single by any formula; a 5 × 10⁻¹¹ s clock requirement rules out any fixed-point representation of 0.1 s with fewer than 35 bits. The decisions are representation decisions and do not change with the formula.

## Check yourself

1. Write the value of a binary64 number in terms of its sign, exponent field and fraction field. What do e = 0 and e = 2047 encode?
2. What are machine epsilon, the unit roundoff and the ulp of 6.778 × 10⁶ in double and in single precision?
3. State the rounding model for one operation and say which of +, −, ×, ÷ can lose many digits in one step, and under what condition.
4. Predict the relative error of 1 − cos θ computed naively at θ = 10⁻⁴ and name the stable form.
5. Why does the textbook quadratic formula fail for x² + 10⁸x + 1 = 0 and what replaces the failing subtraction?
6. The one-pass variance of readings with mean 400 km and standard deviation 1 cm: estimate its relative error in double and say whether a two-pass formula in single precision would fix it.
7. Write the forward-difference error model, give h* for f ≈ f″ ≈ 1 and the best error achievable, and say what the observed order is on each side of h*.
8. What is the relative condition number, what is it for 1 − cos x near 0, and what does that tell you about who is to blame for the lost digits?
