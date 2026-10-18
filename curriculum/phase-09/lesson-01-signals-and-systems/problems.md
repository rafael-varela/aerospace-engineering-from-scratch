# P09.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: energy, power, even and odd)

(a) Compute by hand the energy of x(t) = (s/10) e^{−2t} u(t). (b) Compute the average power of y(t) = 3 cos(2π s t) + 1, and say why the cross term 6 cos(2π s t) contributes nothing. (c) Write the even and odd parts of z(t) = e^{−t}u(t) and give the energy of each. (d) Check all three with `energy`, `average_power`, `even_part` and `odd_part` on a fine grid, and say how the Riemann-sum error depends on dt.

Answer form: three numbers, two functions with their energies, one sentence on the error.

## 2. (Derivation: four short proofs)

(a) Prove the BIBO theorem of 2.4 in both directions, writing out the worst-case input for the necessity part. (b) Show that y[n] = x[n] − b (b ≠ 0) is not linear, but that the difference of the outputs for two inputs depends linearly on the difference of the inputs ("incrementally linear"). (c) Show by hand, with a two-line counterexample each, that the decimator y[n] = x[2n] is neither time-invariant nor causal. (d) Show that x(at − b) is obtained by shifting by b then scaling by a, or by scaling by a then shifting by b/a, and that scaling then shifting by b gives x(at − ab).

Answer form: four derivations.

## 3. (Unfamiliar scenario: a star tracker's output filter)

A star tracker reports a spacecraft attitude angle at 10 Hz. Its software applies, in turn, (i) a three-sample **median filter** y[n] = median(x[n−2], x[n−1], x[n]) to reject single-sample outliers from stray light, and (ii) a **deadband** y[n] = 0 if |x[n]| < d, else x[n], with d = 0.5 arcsec, so the attitude controller ignores errors smaller than the sensor noise.

(a) Classify both blocks by hand for the five properties, with one line of reasoning per entry. (b) Write both as Python callables in the zoo format and run `classify`; compare. (c) Which of the two could you replace by a transfer function in a frequency-domain analysis, and over what input range?

Answer form: a 2 × 5 table with reasons, the code's table, one paragraph.

## 4. (Numerical, uses your code: bias drift and the record-length trap)

(a) A gyro has an uncorrected bias b = 0.002 s deg/s. Using `accumulator` at 100 Hz, compute the attitude error after 10 minutes and after 1 hour, and check by hand. (b) The low-pass y[n] = a y[n−1] + (1 − a)x[n] with a = 1 − 1/(10 s) has a time constant of 10 s samples. Run `is_bibo_stable` with n = 10, 20, 40, …, 1280 (doubling) and report the smallest n that gives the correct verdict. Express it in time constants. (c) Explain in two sentences why the test needs a record several time constants long.

Answer form: two angles in deg, a table of eight verdicts, a ratio, two sentences.

## 5. (Judgment: can a property test be trusted?)

For each of the five testers, describe one system for which the tester (with the lesson's default settings) returns the wrong verdict, and say what change to the test would expose it. Then answer: a reviewer proposes certifying a flight-software filter as "linear and stable" because `classify` says so. Write the paragraph you would send back.

Answer form: five short counterexamples with fixes, one paragraph.

## 6. (Hard: a marginally stable resonator)

The second-order recursion y[n] = 2r cos θ · y[n−1] − r² y[n−2] + x[n] models a lightly damped structural mode sampled in a flight computer.

(a) Show that its impulse response is h[n] = rⁿ sin((n + 1)θ)/sin θ (verify the recursion and the first two values). (b) Show that Σ|h| is finite if and only if r < 1, and compute Σ|h| over 5000 samples for r = 0.95, θ = 2π/50. (c) For r = 1, show that the bounded input x[n] = cos(θn) produces an output that grows like n, while a unit step produces a bounded output. (d) Run `is_bibo_stable` for (r, θ) = (0.95, 2π/50), (1, 2π/50) and (1, 2π/37). Which of the four bounded test inputs exposes the instability in each unstable case, and what does that tell you about choosing test inputs?

Answer form: a verification, a proof and a number, a short derivation, a 3-row table with an explanation.
