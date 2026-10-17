# P03.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which are bounds, and which are observed.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the grid at orbital altitude)

A satellite's Earth-centred position has magnitude r = (6500 + 10 s) km.

(a) Give ulp(r) in double and in single precision, by hand from 2^(⌊log₂r⌋ − 52) and 2^(⌊log₂r⌋ − 23), and check with `ulp` and `numpy.spacing(numpy.float32(r))`.
(b) How many decimal digits of r does each format hold?
(c) A relative-navigation requirement is 1 cm on the separation of two vehicles whose absolute positions are both stored at this magnitude and subtracted. Is it met in double? In single? What representation change meets it in single?

Answer form: four numbers with units, two yes/no answers with reasons, one sentence.

## 2. (Derivation: the error models of the lesson)

(a) From the rounding model fl(x op y) = (x op y)(1 + δ), |δ| ≤ u, show that if x̂ = x(1 + δ₁) and ŷ = y(1 + δ₂) then the relative error of fl(x̂ − ŷ) is bounded by u + (|x| + |y|)u/|x − y|, and say when that bound is large.
(b) Hence derive equation (3) of the lesson, the relative error ε/(θ²/2) of 1 − cos θ computed naively, and show that 2 sin²(θ/2) has relative error at most about 3u for every θ.
(c) Show that the one-pass variance E[x²] − (E[x])² has relative error of order ε m²/σ² when the data have mean m and standard deviation σ, and that the two-pass formula does not.
(d) Starting from the Taylor remainder and the rounding of f(x), f(x + h) and x + h, derive the forward-difference bound (5), its minimiser h* = 2√(ε F₀/F₂) and the minimum error 2√(ε F₀ F₂). Evaluate both for F₀ = F₂ = 1.

Answer form: four short derivations; two numbers in (d).

## 3. (Unfamiliar scenario: GPS time of week in single precision)

A receiver stores the GPS time of week (seconds since Sunday 00:00, 0 to 604 800 s) as a single-precision float and computes pseudoranges by multiplying time differences by c = 299 792 458 m/s.

(a) What is the spacing of representable times at the end of the week, and at 1 h into the week? Convert both to a range error at the speed of light.
(b) Repeat for double precision at the end of the week.
(c) A colleague says "single is fine, our clock is only good to 10 ns anyway". Explain in two sentences why the representation error is worse than the clock error, and propose a representation that keeps single precision for the fractional part.

Answer form: four numbers with units, two sentences, one representation.

## 4. (Numerical, uses your code: the order of summation error)

Sum n copies of the single-precision number 0.1 for n = 10³, 10⁴, 10⁵, 10⁶ with `naive_sum`, `pairwise_sum` and `kahan_sum` (all in `numpy.float32`). The exact sum of the stored values is n × float32(0.1), which you can compute in double.

(a) Tabulate the absolute error of each method against n.
(b) Fit the observed order of the naive error against n with `observed_order`. It is not 1. Explain, from the rounding model, why equal summands make the naive error grow faster than linearly in n.
(c) Separately, how far is n × float32(0.1) from n × 0.1 at n = 10⁶? Which of the two errors (representation or accumulation) would a Kahan sum remove, and which would it not?

Answer form: a 4 × 3 table, one observed order with a two-sentence explanation, two numbers.

## 5. (Numerical, uses your code: an altimeter record that single precision cannot hold)

Generate 10⁴ altimeter readings h = 400 000 m + σ·z with z standard normal from `numpy.random.default_rng(s)` and σ = 0.01 + 0.001 s metres (synthetic data).

(a) Compute the variance with `variance_one_pass` and `variance_two_pass` in double and in single, and against `numpy.var`. Report correct digits for each of the four.
(b) The two-pass formula fails in single precision although it is the stable algorithm. Explain why from ulp(4 × 10⁵) in binary32 and state what the data would have to look like for single precision to work.
(c) Shift the data by 400 000 m before converting to single and repeat. Which results change?

Answer form: four digit counts, one explanation, one short table.

## 6. (Judgment: a proposal to go single precision)

A flight-software team proposes computing the whole navigation and guidance chain in single precision to halve the processor load. The chain contains: Earth-centred positions (10⁷ m, 1 cm requirement), velocities (10⁴ m/s, 1 mm/s requirement), quaternion attitude (components of order 1, 10⁻⁶ requirement), the Clohessy–Wiltshire propagation of a 1 km relative state over 1 s steps, the mean of 1000 accelerometer samples per second, and a mission elapsed time counter (10⁶ s, 1 ms requirement).

For each item say whether single precision meets the requirement as stated, and if not, whether a change of *formula* or a change of *representation* rescues it. One sentence each. Then state the single rule you would write into the coding standard.

Answer form: six verdicts with reasons, one rule.

## 7. (Proof: why Kahan's correction is exact)

(a) Prove Sterbenz's lemma: if x and y are doubles with y/2 ≤ x ≤ 2y, then x − y is exactly representable (no rounding occurs). Hint: compare the exponents of x, y and x − y, and count significand bits.
(b) In the Kahan loop, t = fl(s + y) and c = fl(fl(t − s) − y). Use (a) to argue that t − s is computed exactly when |y| ≤ |s|, and hence that c is the rounding error of the addition up to one further rounding.
(c) Explain why this makes the Kahan error bound independent of n while the naive bound grows as n.

Answer form: a proof, an argument, one paragraph.
