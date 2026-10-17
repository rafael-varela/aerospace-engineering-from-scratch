# P02.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: your own chute)

A chute runs from (0, 0) to x_f = 3 + s/10 m along and y_f = 1.5 m down.

(a) Solve (θ_f − sin θ_f)/(1 − cos θ_f) = x_f/y_f by hand to three figures (bisection with a calculator is fine), then with `cycloid_parameter`. Give k.
(b) Give the cycloid time from (7) and the ramp time, and the percentage saving.
(c) Does the track dip below the destination? If so, how deep and where?
(d) Confirm (b) with `travel_time` on the 401-point cycloid and state the discretisation error.

Answer form: θ_f, k, two times, one percentage, a yes/no with two numbers.

## 2. (Derivation: from the functional to the cycloid, and Snell)

(a) Derive the Euler–Lagrange equation from δJ = 0, stating where the end conditions on η are used.
(b) Prove Beltrami's identity (3) and explain why it fails when F depends explicitly on x.
(c) Apply (3) to the brachistochrone integrand and obtain y(1 + y′²) = 2k. Then substitute y = k(1 − cos θ) and integrate to get x = k(θ − sin θ).
(d) Show that along the cycloid sin θ_path/v = 1/√(4gk), where θ_path is the angle of the track from the vertical, and identify this with Snell's law.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a sound ray in a temperature inversion)

On a still night the air near a runway is colder at the ground than above it, so the speed of sound increases with height: c(y) = c₀ + a·y with c₀ = 340 m/s, a = s/1000 s⁻¹, y the height. Fermat's principle says sound takes the path of least time, and for a medium stratified in y the extremals are circular arcs.

(a) Show from (8) that a ray that is horizontal at its highest point h leaves the ground at an angle θ₀ from the vertical with sin θ₀ = c₀/(c₀ + ah), and that the ray is a circle of radius (c₀ + ah)/a centred on the line y = −c₀/a.
(b) A source and a receiver are both on the ground 2000 m apart. Give the radius of the arc joining them, its peak height, and the launch elevation above the horizontal.
(c) Build the arc as a 2001-point path and compute its travel time with `path_time` (speed field c₀ + ay) and compare with the straight ground path. Which is faster, and by how many microseconds? Why is the difference so small although the arc is longer?
(d) Check with `snell_invariant` that sin θ/c is constant on the arc to 10⁻¹² and explain why the ground path is *not* an extremal even though it is straight.

Answer form: a derivation, three numbers, two times with a difference, one check.

## 4. (Numerical, uses your code: Ritz on your chute)

For the chute of problem 1:

(a) Run `ritz_minimum_time` with 1, 2, 4 and 8 sine modes. Tabulate T and the excess over the cycloid time.
(b) Is the sequence monotone? Must it be? What is the upper-bound property that guarantees every entry is above the cycloid time?
(c) Why does the excess fall so slowly? Look at the slope of the Ritz path at x = 0 and compare with the cycloid's.
(d) Design a two-parameter family that includes the right behaviour at the start (hint: x^{2/3}), minimise over it with `scipy.optimize.minimize` and `travel_time`, and report the excess.

Answer form: a four-row table, three short answers, one number.

## 5. (Numerical, uses your code: stationary or not?)

(a) For the chute of problem 1, compute `perturbed_time_increase` for ε = 0.1, 0.05, 0.025 m on the cycloid and on the circular arc that starts vertically (radius (x_f² + y_f²)/(2x_f)). Fit the exponent for each.
(b) Explain, from the expansion J[y + εη] = J + εδJ + ½ε²δ²J + …, what the two exponents mean.
(c) For y_f = 2 m and x_f = 2 + s/10 m, decide without the code whether the fastest track dips below its end, then confirm with `cycloid_parameter` and give the depth of the lowest point.
(d) For the dipping case, the part climbs at the end. Explain in two sentences why climbing can still be time-optimal.

Answer form: six numbers and two exponents, two short paragraphs, one decision with two numbers.

## 6. (Judgment: building the chute)

(a) Rank the assumptions (frictionless, point mass, released from rest, planar track) by how much each would change the predicted 25 % saving of the lesson's chute, with a direction and a rough number for the first.
(b) The cycloid starts vertically. Propose a practical modification for a part with friction coefficient 0.2 and say what it costs in time, in words.
(c) Describe a bench test with two tracks and a steel ball that would validate the model, and say what measured ratio of times would make you distrust it.
(d) A colleague says "the Euler–Lagrange equation proves the cycloid is fastest". Reply in three sentences, naming what the equation proves and what else the lesson needed.

Answer form: four short paragraphs, each with at least one number.
