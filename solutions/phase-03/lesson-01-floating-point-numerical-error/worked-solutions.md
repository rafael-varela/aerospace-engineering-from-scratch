# P03.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. The grid at orbital altitude (s = 15: r = 6650 km = 6.65 × 10⁶ m)

**Derivation.** 2²² = 4 194 304 ≤ 6 650 000 < 2²³, so ⌊log₂r⌋ = 22. (a) ulp₆₄ = 2^(22 − 52) = 2⁻³⁰ = **9.31 × 10⁻¹⁰ m**; ulp₃₂ = 2^(22 − 23) = **0.5 m**. `ulp(6.65e6)` and `numpy.spacing(numpy.float32(6.65e6))` return the same. (b) digits = −log₁₀(ulp/r): **15.9** in double, **7.1** in single. (c) The subtraction of two stored positions is exact (Sterbenz), so its error is the rounding already in the operands: up to one ulp of 6.65 × 10⁶ m. Double: 10⁻⁹ m ≪ 1 cm, **met**. Single: up to 0.5 m, **not met** by a factor 50. Fix: store positions *relative to the station* (or to a reference point refreshed every orbit); at 1 km magnitude the binary32 ulp is 6 × 10⁻⁵ m and the centimetre requirement is met with 150 to spare.

**Check.** 9.31 × 10⁻¹⁰ × 2²⁹ = 0.5 ✓ (the two ulps differ by exactly 2²⁹ because the significands differ by 29 bits).

**Common mistakes.** Using the ulp of the separation (1 km) instead of the absolute position; quoting 16 and 7 digits without the −log₁₀ of the *relative* spacing (the number of digits depends on where in the binade r sits: 15.9 here, 16.2 just above a power of two).

## 2. The error models

(a) fl(x̂ − ŷ) = [x(1 + δ₁) − y(1 + δ₂)](1 + δ₃) = (x − y)(1 + δ₃) + (xδ₁ − yδ₂)(1 + δ₃). Subtracting x − y and dividing by |x − y|: relative error ≤ |δ₃| + (|x||δ₁| + |y||δ₂|)(1 + u)/|x − y| ≤ u + (|x| + |y|)u/|x − y| + O(u²). The second term is large exactly when |x − y| ≪ |x| + |y|: cancellation. The subtraction itself (δ₃) contributes only u; it *reveals* δ₁ and δ₂.

(b) x = 1 is exact (δ₁ = 0); y = fl(cos θ) has |δ₂| ≤ u and |y| ≤ 1, so absolute error ≤ u; x − y = 1 − cos θ = θ²/2 − θ⁴/24 + … The bound of (a) gives relative error ≤ u + u/(θ²/2) ≈ **ε/θ²** = ε/(θ²/2)·½; the lesson's (3) quotes the range u/(θ²/2) to ε/(θ²/2) to allow a libm error of up to one ulp. For 2 sin²(θ/2): θ/2 is exact (a power-of-two scaling), fl(sin) = sin(1 + δ₁), squaring gives (1 + δ₁)²(1 + δ₂) ≈ 1 + 2δ₁ + δ₂, the factor 2 is exact: relative error ≤ **3u** independent of θ.

(c) E[x²] = m² + σ². Computed as a mean of n squares, it carries relative error of order ε (times a modest factor for the summation), so absolute error ~ε(m² + σ²) ≈ ε m². (E[x])² likewise has absolute error ~ε m². The difference σ² therefore has absolute error ~ε m² and **relative error ~ε m²/σ²**. The two-pass form computes d = x − m first; m is known to relative u, so d has absolute error ~u m, which is relative u m/σ on d, then squaring and summing positives adds O(ε): relative error ~ε m/σ, smaller by the factor m/σ than the one-pass form (for m/σ = 4 × 10⁷ that is 10⁻⁸ against 10 for the one pass: still worth caring about, which is why `numpy.var` uses it and Welford's update does even better).

(d) D(h) = (f(x + h) − f(x))/h with exact arithmetic equals f′(x) + h f″(ξ)/2, |truncation| ≤ h F₂/2. Each function value is rounded: absolute error ≤ u F₀ each, so ≤ 2u F₀ = ε F₀ on the difference. x + h is rounded with absolute error ≤ u|x|, which shifts f(x + h) by ≤ u|x f′| ≤ u F₀. Total rounding ≤ (ε + u) F₀ ≤ 2ε F₀ (allowing a further ulp for the libm), divided by h. So E(h) = h F₂/2 + 2ε F₀/h, (5). dE/dh = F₂/2 − 2ε F₀/h² = 0 at **h* = 2√(ε F₀/F₂)**, and E(h*) = √(ε F₀ F₂) + √(ε F₀ F₂) = **2√(ε F₀ F₂)**. For F₀ = F₂ = 1: h* = 2√(2.22 × 10⁻¹⁶) = **2.98 × 10⁻⁸**, E(h*) = **2.98 × 10⁻⁸**.

**Check.** `fd_optimal_step(1, 1)` = 2.98 × 10⁻⁸; `main()` finds the observed minimum 2.97 × 10⁻⁹ at h = 10⁻⁸, a factor 10 under the bound and within a decade of h*.

**Common mistakes.** In (a) forgetting δ₃ (it is harmless but belongs in the bound); in (b) claiming the stable form has error u/θ² "because sin also rounds" (the rounding is relative, and the result is not a difference); in (d) omitting the x + h rounding (it matters when |x| ≫ 1: for sin at x = 100 it moves h* to 3 × 10⁻⁷).

## 3. GPS time of week in single precision

**Derivation.** 2¹⁹ = 524 288 ≤ 604 800 < 2²⁰, so e = 19. (a) ulp₃₂ = 2^(19 − 23) = 2⁻⁴ = **0.0625 s**; × c = **1.87 × 10⁷ m** (18 700 km, more than the Earth's diameter). At 1 h = 3600 s, 2¹¹ ≤ 3600 < 2¹², ulp₃₂ = 2^(11 − 23) = **2.44 × 10⁻⁴ s**, × c = **73 km**. (b) ulp₆₄ at 604 800 = 2^(19 − 52) = **1.16 × 10⁻¹⁰ s**, × c = **3.5 cm**. (c) A 10 ns clock error is 3 m of range; the representation error of single precision is 73 km to 18 700 km, four to seven orders worse, and it is not noise but a systematic quantisation that changes as the week advances (the ulp doubles at every power of two). Representation: an integer count of whole seconds (or milliseconds) plus a single-precision fraction in [0, 1); the fraction has ulp 6 × 10⁻⁸ s at most, i.e. 18 m, still worse than the clock, so a double for the fraction (ulp 10⁻¹⁶ s) is the honest choice, and GPS receivers do keep time in an integer-plus-fraction form.

**Check.** The 1 h result scales as expected: 3600/604 800 ≈ 2⁻⁷·(1.1), and 0.0625/2.44 × 10⁻⁴ = 256 = 2⁸ (the extra factor 2 is where in the binade each sits) ✓.

**Common mistakes.** Using the relative precision 1.2 × 10⁻⁷ × 604 800 = 0.072 s (close, but the grid is a power of two, 0.0625 s); forgetting that the fraction alone in single precision is still 18 m.

## 4. The order of summation error (single precision, summands float32(0.1) = 0.100000001490116…)

| n | exact n × float32(0.1) | naive error | pairwise error | Kahan error |
|---|---|---|---|---|
| 10³ | 100.0000015 | 9.55 × 10⁻⁴ | 1.49 × 10⁻⁶ | 1.49 × 10⁻⁶ |
| 10⁴ | 1 000.0000149 | 9.71 × 10⁻² | 1.49 × 10⁻⁵ | 1.49 × 10⁻⁵ |
| 10⁵ | 10 000.000149 | 1.44 | 1.49 × 10⁻⁴ | 1.49 × 10⁻⁴ |
| 10⁶ | 100 000.00149 | 958 | 9.30 × 10⁻³ | 1.49 × 10⁻³ |

(a) As tabulated (`naive_sum`, `pairwise_sum`, `kahan_sum` with `dtype=np.float32`). The Kahan and pairwise errors are the half-ulp of the float32 *result* (ulp(100) = 7.6 × 10⁻⁶, ulp(10⁵) = 7.8 × 10⁻³): the last rounding, which no algorithm can avoid; pairwise at 10⁶ shows one extra ulp.

(b) `observed_order([1e3, 1e4, 1e5, 1e6], naive errors)` = **1.92**. Once the running sum S is much larger than 0.1, adding 0.1 lands between grid points spaced ulp(S) ≈ 1.2 × 10⁻⁷ S apart and rounds to the nearest; because every summand is identical, the rounding goes the *same way* at every step within a binade instead of averaging out, so the error per step is a fixed fraction of ulp(S) ∝ S ∝ n, and the sum over n steps is ∝ n². The naive worst-case bound n u Σ|x| = n² u × 0.1 is actually attained here (at n = 10⁶: 6 × 10⁻⁸ × 10¹² × 0.1 = 6000, against 958 observed: the bound within a factor 6), which is the danger of equal summands such as fixed time steps.

(c) n × float32(0.1) − n × 0.1 = 10⁶ × 1.49 × 10⁻⁹ = **1.49 × 10⁻³**: the stored constant is 1.49 × 10⁻⁸ relative too large and that bias is simply multiplied by n. Kahan summation removes the *accumulation* error (958 → 0.0015, which is the final rounding) and cannot touch the *representation* error of the summand: a loop that adds float32(0.1) a million times with perfect summation still reaches 100 000.0015, not 100 000.

**Check.** The naive error ratio between successive decades is 102, 15, 665: not uniform (the n² law holds within binades and jumps at binade boundaries) but averaging to 10^1.92 per decade ✓.

**Common mistakes.** Comparing against n × 0.1 instead of n × float32(0.1) and so blaming the algorithm for the representation; expecting √n growth (that is for *random* rounding errors).

## 5. An altimeter record single precision cannot hold (s = 15: σ = 0.025 m, 10⁴ samples, `default_rng(15)`; SYNTHETIC)

(a) `numpy.var` = 6.245386 × 10⁻⁴ m². One-pass double: **0.0** (zero correct digits: E[x²] and m² agree to every stored bit); two-pass double: 6.245386 × 10⁻⁴ (**14.4 digits**); one-pass single: 8.36 × 10⁵ (**no digits**, off by 10⁹); two-pass single: 185.6 (**no digits**, off by 3 × 10⁵).

(b) ulp(4 × 10⁵) in binary32 = 2^(18 − 23) = **0.03125 m** > σ = 0.025 m. Cast to single, every reading is one of a handful of grid values (400 000.00, 400 000.03, 399 999.97, …); the first sample's true deviation of −0.0358 m is stored as −0.03125. The two-pass formula then computes the variance of the *quantised* data correctly, and since the quantisation error of ±0.016 m is comparable to σ itself the answer is wrong in the first digit, and the mean in single precision carries the sum of 10⁴ terms of 4 × 10⁵ (grid spacing at 4 × 10⁹ is 256) so it is off by 13.6 m, which the squared deviations then turn into the 185 m² reported. Single precision would work only if ulp(mean) ≪ σ, i.e. the mean below about 10⁻³ × σ/ε₃₂ ≈ 2 × 10² m for σ = 2.5 cm.

(c) Shifted data d = h − 400 000 m (magnitude ~0.025, binary32 ulp 1.9 × 10⁻⁹ m): one-pass double 6.245386 × 10⁻⁴ (14.4 digits); two-pass double unchanged; one-pass single 6.245386 × 10⁻⁴ (**6.9 digits**); two-pass single 6.245374 × 10⁻⁴ (**5.7 digits**). Every result is now within the dtype's precision; single precision holds the shifted data to 7 digits as it should.

**Check.** The variance is shift-invariant mathematically, so (c) must equal (a)'s correct value; it does to 14 digits in double ✓.

**Common mistakes.** Blaming the algorithm in (b) (the two-pass form is stable; the *data* were destroyed by the cast); forgetting that the shift must be applied before the cast to single, not after.

## 6. A proposal to go single precision

| item | binary32 ulp at magnitude | requirement | verdict | fix |
|---|---|---|---|---|
| ECI position 10⁷ m | 1.0 m | 1 cm | fails by 100 | representation: store relative to a reference point or as integer metres + float fraction |
| velocity 10⁴ m/s | 10⁻³ m/s | 1 mm/s | marginal (one ulp = the requirement) | representation: store relative to the reference orbit's velocity, or double |
| quaternion, components ≤ 1 | 1.2 × 10⁻⁷ (at 1), finer below | 10⁻⁶ | meets, with 8× margin | none; renormalise periodically |
| CW propagation, 1 km state, 1 s steps | 6 × 10⁻⁵ m at 1 km; but 1 − cos(10⁻³) in single = 4.77 × 10⁻⁷ against the true 5.0 × 10⁻⁷ (5 % error) | 1 cm after 10⁴ steps | fails via the *formula* | formula: 2 sin²(θ/2) or the series; the state itself fits |
| mean of 1000 samples of ~10 m/s² | sum ~10⁴, ulp 10⁻³ m/s², mean error ~10⁻⁶ m/s² | typical 10⁻⁵ m/s² | meets, barely; grows with sample count | pairwise or Kahan summation if the count rises |
| mission time 10⁶ s | 0.0625 s | 1 ms | fails by 60 | representation: integer seconds + float fraction, or integer milliseconds |

Rule for the coding standard: *a single-precision variable may hold a quantity only if 10⁻⁷ × its largest magnitude is at least ten times below its required absolute precision; otherwise store it as an integer part plus a fractional part, or relative to a reference, and never subtract two large stored values to get a small one.* The team's processor-load saving is real for the quaternion and sensor averaging and illusory for the state vector, which is where a navigation computer spends little of its time anyway.

**Common mistakes.** Judging by relative precision ("single has 7 digits, that is plenty") without asking what magnitude the variable *holds*; treating the CW item as a representation problem when the 1 km state fits fine and the formula is at fault.

## 7. Why Kahan's correction is exact

(a) **Sterbenz.** Let x ≥ y > 0 be doubles with y ≤ x ≤ 2y (the other half of the hypothesis is the symmetric case). Then 0 ≤ x − y ≤ y. Write y = Y·2^(e_y − 52) with integer Y < 2⁵³ and e_y the exponent of y. Since x ≥ y, e_x ≥ e_y, so ulp(x) is an integer multiple of ulp(y) and x = X·2^(e_y − 52) for an integer X. Then x − y = (X − Y)·2^(e_y − 52) with 0 ≤ X − Y ≤ Y < 2⁵³: an integer of at most 53 bits times a power of two, which is exactly a double (or zero). No rounding occurs. ∎

(b) In the Kahan step with |y| ≤ |s|: t = fl(s + y) satisfies |t − s| ≤ |y|(1 + u) ≤ |s|(1 + u), so t lies within a factor 2 of s (for s + y and s of the same sign, which |y| ≤ |s| guarantees), and by (a) t − s is computed **exactly**: it equals the part of y that was actually absorbed into the sum. Then (t − s) − y is the exact rounding error of the addition, up to the single rounding of that subtraction, which is of relative size u on a quantity already of size u|s|: an error of order u²|s|. So c captures the lost low-order bits of y to within u² relative, and feeding −c into the next addition restores them.

(c) In the naive sum each of the n additions contributes an uncompensated error of up to u|S_k|, and these add to ≤ n u Σ|xᵢ| in the worst case. In the compensated sum each addition's error is recovered to within u² |S_k| and re-injected; what remains per step is O(u²)|S_k|, summing to O(n u²) Σ|xᵢ|, plus the final result's own rounding O(u) Σ|xᵢ|. Hence the bound ≈ (2u + O(n u²)) Σ|xᵢ|, and n u² < u for every n < 1/u = 9 × 10¹⁵: for any sum you will ever compute, the bound is independent of n. That is the sense in which Kahan summation "beats" the arithmetic: it does not have more bits, it carries the lost ones in a second variable.

**Common mistakes.** Claiming the compensation is exact with no rounding at all (the final c has one rounding; what is exact is t − s); applying Sterbenz without checking the signs (if s and y have opposite signs and |y| > |s|, the lemma can fail and a variant, Neumaier's, swaps the roles of s and y).

## Check-yourself answers

1. x = (−1)^s (1 + f/2⁵²) 2^(e − 1023) for 0 < e < 2047; e = 0 encodes ±0 and subnormals (no implicit 1, exponent −1022); e = 2047 encodes ±inf (f = 0) and nan (f ≠ 0).
2. ε = 2⁻⁵² = 2.22 × 10⁻¹⁶ (double), 2⁻²³ = 1.19 × 10⁻⁷ (single); u = ε/2; ulp(6.778 × 10⁶) = 2⁻³⁰ = 9.3 × 10⁻¹⁰ in double, 2⁻¹ = 0.5 in single.
3. fl(x op y) = (x op y)(1 + δ), |δ| ≤ u. Only subtraction of nearly equal numbers (or addition of nearly opposite ones) loses many digits in one step, when |x − y| ≪ |x| + |y|; the loss is the exposure of the operands' prior rounding errors.
4. ε/(θ²/2) = 2.2 × 10⁻¹⁶/(5 × 10⁻⁹) ≈ 4 × 10⁻⁸ (about eight digits lost; observed 5 × 10⁻⁹). Stable form: 2 sin²(θ/2).
5. √(b² − 4ac) ≈ b = 10⁸ to within 2 × 10⁻⁸, which is below the double grid spacing 1.5 × 10⁻⁸ at 10⁸, so −b + √D is a difference of two numbers equal to the grid: 25 % error. Replace it by Vieta: x₂ = c/q with q = −(b + sign(b)√D)/2.
6. ε m²/σ² = 2.2 × 10⁻¹⁶ × (4 × 10⁵)²/(10⁻²)² = 0.36: no correct digits. A two-pass formula in single precision would not help, because binary32 cannot store 400 km to better than 3 cm; the data must be shifted (stored as deviations from 400 km) before any formula is applied.
7. E(h) = h F₂/2 + 2ε F₀/h; h* = 2√ε ≈ 3 × 10⁻⁸ for F₀ = F₂ = 1, with best error 2√ε ≈ 3 × 10⁻⁸ (half the digits); observed order +1 for h > h* (truncation), −1 for h < h* (rounding).
8. κ = |x f′(x)/f(x)|, the amplification of relative input error. For 1 − cos x near 0, κ → 2: the problem is well conditioned, so the digits lost by the naive formula are the algorithm's fault, and a stable algorithm (2 sin²(x/2)) recovers them.
