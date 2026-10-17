# P09.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Sketch every locus by hand before computing it.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your altitude loop)

Take the lesson's altitude plant with the flight-path lag T = 0.8 + s/50 s (V = 50 m/s, τ = 0.5 s).

(a) Give the open-loop poles, the real-axis segments, the asymptote angles and centroid.
(b) Find the breakaway point and the gain there.
(c) Find the imaginary-axis crossing (gain and frequency) by hand with Routh, and confirm with `jw_crossings`.
(d) Find the gain for ζ = 0.5 with `gain_for_damping`, the third pole, and compare the predicted and simulated overshoot and settling time.

Answer form: three poles, two segments, angles and a centroid; a point and a gain; a gain and a frequency; a gain, a pole and two pairs of metrics.

## 2. (Derivation: the rules)

(a) Show that a point s is on the locus for some K > 0 if and only if ∠G(s) = 180° (mod 360°), and that the gain there is 1/\|G(s)\|.
(b) Derive the real-axis rule from the angle condition.
(c) Show that for large \|s\| the locus approaches straight lines at angles (2k + 1)180°/(n − m) through σ_a = (Σp − Σz)/(n − m). (Hint: expand D/N for large s to two terms.)
(d) Show that breakaway points satisfy N D′ − N′D = 0, and that when n − m ≥ 2 the sum of the closed-loop poles does not depend on K.

Answer form: four derivations.

## 3. (Unfamiliar scenario: an unstable airframe)

A launch vehicle under thrust-vector control, or a statically unstable airframe, has a pitch plant with a real pole in the right half plane. Take L = K(s + 1)/((s² − 4)(0.1s + 1)): aerodynamic poles at ±2 s⁻¹, an actuator pole at −10 s⁻¹ and a controller zero at −1 s⁻¹.

(a) Sketch the locus: segments, asymptotes, breakaway.
(b) Use Routh to find the range of K for which the closed loop is stable. Where does the locus cross the imaginary axis?
(c) Find the gain for ζ = 0.5 and the remaining pole. Is the ζ = 0.5 pair dominant?
(d) What does it mean physically that the loop is unstable for small K?

Answer form: a sketch with numbers; a range and a crossing point; a gain, two poles and a verdict; one sentence.

## 4. (Numerical, uses your code: the climb-rate feedback)

(a) Run `pd_design_table` for T_d = 0.5 + s/10 s, 1 s and 3 s, and tabulate K, the pair, the third pole, the predicted and simulated overshoot and settling.
(b) Explain each mismatch between prediction and simulation from the positions of the third pole and the zero.
(c) Which T_d gives the fastest settling? Is it the one with the largest ω_n?

Answer form: a table, three explanations, an answer with a reason.

## 5. (Routh's special cases)

For each polynomial build the Routh array by hand, count the right-half-plane roots, and check with `routh_rhp_count` and `np.roots`:

(a) s⁴ + 2s³ + 3s² + 4s + 5
(b) s³ + 3s + 2
(c) s⁵ + s⁴ + 4s³ + 4s² + 3s + 3
(d) 2s⁴ + s³ + 2s² + s + 1
(e) s³ + 3s² + 3s + 1

Answer form: five arrays (first columns are enough), five counts, and a note on which special case each one triggers.

## 6. (Judgment: is ζ = 0.5 enough?)

(a) With the lesson's design gain K = 0.01037 rad/m fixed, compute ζ and the critical gain for τ = 0.4, 0.5, 0.7 s and for V = 40, 50, 60 m/s. Which uncertainty matters more?
(b) The T_d = 1 s design uses a measured climb rate. What does noise on that signal do, and why is the closed form "exact" only in the model?
(c) A pilot asks for the altitude to settle "within 5 s, with no more than 10 % overshoot" with P control alone. Is it possible on this plant? Use the locus to argue.

Answer form: a table and a sentence; two sentences; a yes or no with a locus-based argument.
