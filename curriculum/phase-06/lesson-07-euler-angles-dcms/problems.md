# P06.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an orbit as a 3-1-3 rotation)

The perifocal frame of an orbit (axis 1 towards periapsis, axis 3 along the orbit normal) is reached from the inertial frame by the 3-1-3 sequence (Ω, i, ω). Take Ω = (10 + s)°, i = 51.6° and ω = 30°.

(a) The DCM C_p/i with `dcm_from_euler`, its orthogonality error and its determinant.
(b) The inertial components of the orbit normal, from the matrix and from the formula (sin i sin Ω, −sin i cos Ω, cos i). Which row or column of C is it, and why?
(c) The 3-2-1 angles of the same attitude, and the 3-1-3 angles extracted back.

Answer form: a matrix and two numbers; a vector and a sentence; two angle triples.

## 2. (Derivation: rates, determinants and the kinematic equation)

(a) From (2), derive the 3-2-1 body rates (3), then invert them to (4).
(b) Show that |det B| = |cos a₂| for an asymmetric sequence and |sin a₂| for a symmetric sequence, using the geometry of the three columns.
(c) Derive Ċ = −[ω×]C from CCᵀ = E, and show that the exact solution keeps CᵀC = E.
(d) Derive the coning body rates (6) from (2) for the 3-1-3 sequence with constant ψ̇, β and φ̇.
(e) For a constant ω, show that one RK4 step of (5) multiplies C by the first five terms of the series for exp(−[ω×]Δt), and estimate the orthogonality error that this leaves after one step.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a rocket on the pad)

A launch vehicle stands vertically: its body x axis points up (−down), its body z axis points north.

(a) Write its DCM relative to north–east–down by thinking about rows (row i is body axis i in NED components). Confirm with `dcm_from_euler([0, 90°, 0], "321")`.
(b) For each of the twelve sequences, extract the angles and compute |det B|. Which sequences are singular at this attitude?
(c) During the ascent the vehicle pitches over from vertical to horizontal in the north–down plane. Which of the non-singular sequences of (b) stays non-singular through the whole pitch-over? Check your answer at pitch 90°, 60°, 30° and 0°.
(d) Why do many launch-vehicle guidance systems carry a quaternion instead of any Euler sequence?

Answer form: a matrix; a list of four sequences; a sequence with four determinants; one sentence.

## 4. (Numerical, uses your code: your own coning spacecraft)

A spacecraft cones with ψ̇ = 0.2 + s/100 rad/s, φ̇ = 1.5 + s/10 rad/s and β = (10 + s/2)°.

(a) Its body rate magnitude. Integrate 60 s with `integrate_dcm` and report the largest misalignment from `coning_dcm` and the largest orthogonality error.
(b) RK4 with Δt = 0.04, 0.02, 0.01 and 0.005 s: the final misalignment and orthogonality error, and the observed order of each.
(c) The same with `renormalize_every=1`. What changed and what did not?
(d) Your flight computer runs at 50 Hz. Does RK4 meet a 10⁻⁴ rad attitude budget over 60 s?

Answer form: a rate and two numbers; a table and two orders; a sentence; yes or no with the number.

## 5. (Numerical, uses your code: how fast must the computer run?)

Fly the hammerhead with q = r = 0.5 rad/s and δ = (0.05 + s/100)°.

(a) The largest ψ̇ on the path, estimated from (4) as r/sin δ.
(b) The final misalignment of Euler-angle RK4 and of DCM RK4 (`integrate_segments_rk4`) at 50, 100, 200 and 400 Hz.
(c) The observed order of the Euler-angle error. Is it the RK4 order? Why might it not be at the coarsest rates?
(d) The lowest of these rates that keeps the Euler-angle error below 10⁻⁴ rad. What rate would the DCM need?

Answer form: a number; a table; three orders and a sentence; two rates.

## 6. (Judgment: attitude representations in practice)

For each, give a short argued answer (one to three sentences).

(a) A pilot's attitude display shows yaw, pitch and roll. Should the flight computer integrate those three angles? What should it integrate, and when should it compute the angles?
(b) A colleague's code passes every test for 3-2-1 but extracts wrong angles for 1-2-3. What is the likely bug?
(c) A colleague re-orthonormalises the DCM after every step and concludes that the attitude is now exact. What is wrong with the conclusion?
(d) Henderson's 1977 NASA report tabulates all twelve sequences. Before comparing your matrices with its tables, what must you check, and what do you do about it?
(e) Is gimbal lock a physical phenomenon? Explain the difference between a mechanical gimbal platform and a strapdown computer here.
(f) A spacecraft's mission profile never passes within 30° of the singular attitude of 3-2-1. Is it safe to integrate 3-2-1 angles? What could still make you choose the DCM?

Answer form: six labelled items.
