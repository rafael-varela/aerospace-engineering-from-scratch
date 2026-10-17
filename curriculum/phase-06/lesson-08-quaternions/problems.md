# P06.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: two turns in two orders)

An aircraft yaws 90° (about its axis 3) and then pitches 90° (about its new axis 2). A second aircraft pitches 90° first and then yaws 90°.

(a) By hand with (2) and (6), the quaternion of each final attitude. Check with `quat_from_axis_angle` and `quat_multiply`.
(b) The single axis and angle of each. What do they have in common, and how do they differ?
(c) Where does each nose point (body axis 1 in north–east–down components)? Use a row of C(q).
(d) The body components of the north direction for each.

Answer form: two quaternions; two axes and two angles with a sentence; two vectors; two vectors.

## 2. (Derivation: the quaternion as a rotation)

(a) Expand the product of two quaternions using Hamilton's rules (1) and obtain the scalar–vector form (2).
(b) Prove (p ⊗ q)* = q* ⊗ p* and then |p ⊗ q| = |p||q|.
(c) Expand (4) with (2) and derive the DCM (5). Show that it equals the passive Rodrigues matrix for q = (cos(a/2), sin(a/2) n).
(d) Derive the kinematic equation (8) from the composition rule (6), and show that q · q̇ = 0.
(e) Show that SLERP (10) has unit norm for all u, and that the attitude turns from q₀ to slerp(q₀, q₁, u) by the angle 2uΩ, u times the whole turn 2Ω from q₀ to q₁.

Answer form: five derivations.

## 3. (Unfamiliar scenario: averaging star-tracker samples)

A star tracker reports ten samples of a fixed attitude, each with a random error of about 10⁻³ rad, and with an arbitrary sign: some samples have q₀ < 0. Generate them with `rng = np.random.default_rng(s)`: the true attitude q_t = `quat_from_euler_321` of yaw (20 + s)°, pitch 5°, roll −10°; for each sample draw e = `rng.normal(size=3) * 1e-3`, set q = q_t ⊗ `quat_from_axis_angle(e, |e|)`, and flip its sign if `rng.random() < 0.5`.

(a) The naive average (sum of the ten arrays divided by ten): its norm, and the error angle of its normalised value from q_t.
(b) A fix: flip each sample so that it has a positive dot product with the first, then average and normalise. Its error, and the mean error of a single sample.
(c) The eigenvector method: the eigenvector of Σ qqᵀ with the largest eigenvalue. Its error. Why is it immune to the signs?
(d) Why does averaging improve the error, and by roughly what factor should it for ten independent samples?

Answer form: two numbers; two numbers; a number and a sentence; a factor and a sentence.

## 4. (Numerical, uses your code: quaternion against DCM)

A spacecraft cones with ψ̇ = 0.2 + s/100 rad/s, φ̇ = 1.5 + s/10 rad/s and β = (10 + s/2)° (the case of `P06.07` Problem 4).

(a) Integrate 60 s by `integrate_quat` and by `integrate_dcm`. The largest misalignment between the two routes, and the largest | |q| − 1 |.
(b) RK4 on the quaternion at Δt = 0.04, 0.02 and 0.01 s: the final error angle and the norm drift, with their observed orders.
(c) Compare (b) with your `P06.07` Problem 4 RK4 errors for the DCM at the same steps. Explain the ratio.
(d) A designer proposes renormalising only once a second at Δt = 0.02 s. What largest |q| − 1 and DCM orthogonality error should she expect between renormalisations?

Answer form: two numbers; a table with two orders; three ratios and a sentence; two numbers.

## 5. (Numerical, uses your code: interpolating between star-tracker reports)

Use the lesson's coning spacecraft (ψ̇ = 0.4, φ̇ = 2.0 rad/s, β = 20°). Sample its attitude every h = (0.1 + s/100) s and every h/2.

(a) The largest SLERP and NLERP interpolation errors over 30 s, evaluated every 0.01 s.
(b) A colleague instead interpolates the 3-2-1 angles linearly (after unwrapping them). Its largest error at the same two intervals.
(c) The observed order of each method.
(d) Which method would you put in flight software, and why?

Answer form: four numbers; two numbers; three orders; two sentences.

## 6. (Judgment: quaternions in practice)

For each, give a short argued answer (one to three sentences).

(a) A ground tool reads a star tracker's quaternion as (x, y, z, w) while the attitude software writes (w, x, y, z). What attitude error results, and which test in this lesson would have caught it?
(b) Why is it safe to renormalise a quaternion by dividing by its norm, while the DCM needed a singular-value decomposition?
(c) A telemetry plot of q₀ jumps from +0.7 to −0.7 between two samples while the spacecraft turns smoothly. Is something wrong?
(d) Henderson's NASA report gives quaternion formulas for all twelve Euler sequences. What two conventions must you check before using them, and how will you check them numerically?
(e) When is it acceptable to use NLERP instead of SLERP?
(f) A colleague says quaternions have "no singularities, so they are always better". Name one thing Euler angles do better and one cost of the quaternion.

Answer form: six labelled items.
