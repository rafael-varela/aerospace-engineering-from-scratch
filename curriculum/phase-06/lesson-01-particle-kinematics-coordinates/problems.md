# P06.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the unit vectors before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one radar plot, three descriptions)

A surveillance radar at the origin reports an aircraft at horizontal range r = 10 + s/5 km and azimuth θ = 30° (from east toward north), with ṙ = −150 m/s, θ̇ = 0.010 rad/s, r̈ = 2.0 m/s² and θ̈ = −1.0 × 10⁻⁴ rad/s². The aircraft is level (z constant).

(a) Compute v_r, v_θ and the speed.
(b) Compute a_r and a_θ, and the size of each of the four terms of (3). Which term dominates?
(c) Convert v and a to east–north components with `cylindrical_to_cartesian`.
(d) With `path_kinematics`, find a_t, a_n, the radius of curvature and the load factor contribution a_n/g. Is the aircraft speeding up or slowing down, and is it turning left or right?

Answer form: three numbers; two numbers and four terms; two vectors; four numbers and two sentences.

## 2. (Derivation: turning unit vectors, straight lines and helices)

(a) Starting from e_r = (cos θ, sin θ) and e_θ = (−sin θ, cos θ), derive (1), then (2) and (3). Say where the product rule and where the chain rule are used.
(b) Prove that r²θ̇ is constant for a particle moving in a straight line at constant velocity. For an aircraft at speed V passing the radar at miss distance d, show that θ̇ peaks at V/d and that |θ̈| peaks at (3√3/8)(V/d)². Where on the track does |θ̈| peak?
(c) From (6), derive the curvature of the helix x = R cos ωt, y = R sin ωt, z = wt and show it equals R/(R² + b²) with b = w/ω. What are its limits as w → 0 and as R → 0?

Answer form: three derivations, two maxima with their locations.

## 3. (Unfamiliar scenario: a tracking antenna for a sounding rocket)

A tracking antenna stands on flat ground D = 2.0 km from a launch pad. A sounding rocket rises vertically from rest with constant acceleration a = 30 m/s² (ignore its burnout). The antenna describes the rocket in polar coordinates in the vertical plane: slant range R and elevation angle e.

(a) Write R(t) and e(t), and show that the elevation rate is ė = D ż/(D² + z²).
(b) Find the largest elevation rate the antenna drive must provide, the time at which it occurs, the rocket's height then and the elevation angle. Give the exact height in terms of D.
(c) At that instant, compute the polar acceleration components a_R = R̈ − Rė² and a_e = Rë + 2Ṙė from your R(t), e(t). Show that they equal a sin e and a cos e, and explain why.
(d) What is the largest elevation acceleration the drive must give, and when?

Answer form: two formulas and a derivation; four numbers; two numbers and a sentence; a number and a time.

## 4. (Numerical, uses your code: your own noisy radar)

Use the lesson's climbing turn sampled at 10 Hz for 120 s.

(a) With `accel_error_vs_stride`, seed s, σ_r = 5 m, σ_θ = 1 mrad and σ_z = 5 m, tabulate the RMS acceleration error for strides 10–140 (h = 1–14 s). Find the best step and its error and compare with `optimal_step` using the per-axis noise of the lesson.
(b) Repeat with σ_θ = 2 mrad. By what factor did the best step and the best error change? Predict both factors from (9) first, using the new per-axis noise.
(c) Sample at 20 Hz instead (same σ, seed s). Does the best error improve? Explain with the structure of the three-point stencil.

Answer form: a table, a step, an error and a comparison; two factors predicted and measured; a number and an explanation.

## 5. (Judgment: is that 6 g real?)

A range-safety display differences the same radar's samples with a 0.5 s step and shows the aircraft of the lesson pulling 6 g for a moment.

(a) Using the lesson's noise study, how large an acceleration error should you expect at h = 0.5 s? Is a 6 g reading evidence of a 6 g manoeuvre?
(b) What step would you use instead, and what does it cost you in the timing of a manoeuvre's onset?
(c) Name two things that could produce a false transverse acceleration even with a long step (think about which term of (3) a slowly drifting azimuth bias enters, and about an overhead pass).
(d) What independent measurement would settle it?

Answer form: a number with units and a sentence; a step and a trade-off; two mechanisms; one measurement.
