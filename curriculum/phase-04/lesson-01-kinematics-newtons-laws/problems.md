# P04.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Where a problem uses the lesson's code, say which function and which tolerance.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a toss by hand and by code)

A ping-pong ball (k = 0.1396 m⁻¹) is tossed at v₀ = (5 + s/10) m/s at 40° above the horizontal from floor level.

(a) By hand, without drag: time of flight, range and apex height. Check with `vacuum_state` at t = T.
(b) With drag, using `simulate_projectile`: range, time of flight, apex height, impact speed and impact angle. Give the range as a fraction of the vacuum range.
(c) Which of the five quantities in (b) changed most, in relative terms, and why is the time of flight *shorter* with drag although the ball is slowed?

Answer form: three numbers with units; five numbers with units and a fraction; two sentences.

## 2. (Derivation: the closed forms and their expansions)

(a) From m dv/dt = mg − ½ρC_D A v², derive v(t) = v_t tanh(gt/v_t) and s(t) = (v_t²/g) ln cosh(gt/v_t) for a fall from rest, stating where "from rest" enters.
(b) Using ln cosh u = u²/2 − u⁴/12 + O(u⁶), show that the time to fall h is t ≈ t₀(1 + gh/(6v_t²)) with t₀ = √(2h/g), to first order in gh/v_t². Evaluate the correction for a ping-pong ball at h = 2 m and compare with the exact value from `vertical_drop_time`.
(c) Derive x(t) = ln(1 + kv₀t)/k for horizontal motion with quadratic drag and no gravity, and show that the distance over which the speed halves is ln 2/k independent of v₀.
(d) For the coupled two-dimensional case, prove from (1) that ẋ(t) is strictly decreasing while ẋ > 0, and deduce that the impact angle below the horizontal is steeper than the launch angle above it whenever the launch and landing heights are equal. (Hint: compare |ż| at launch and at impact using the energy argument of §2.6, or compare the horizontal speeds alone.)

Answer form: four derivations; one number compared with the code.

## 3. (Unfamiliar scenario: a supply drop under a parachute)

A drone hovering at 30 m releases a (20 + s/10) kg package under a round parachute of 1.2 m diameter with C_D = 1.5 (based on the canopy's projected area). Air density 1.20 kg/m³.

(a) Terminal speed of the package, time to reach the ground, impact speed, and the fraction of v_t reached at impact. Compare the time and impact speed with a vacuum drop.
(b) The parachute fails to open and the package falls as a 0.4 m lump with C_D = 1.0. Repeat (a). Why is the parachute's effect on the *time* smaller than its effect on the *impact speed*?
(c) The drone is moving horizontally at 10 m/s when it releases the package (parachute open). Using `simulate_projectile` with a horizontal launch from z₀ = 30 m, find where the package lands relative to the release point and its impact velocity components. Compare the horizontal distance with the vacuum value 10 × √(2h/g) and with the pure-horizontal-decay distance ln(1 + kv₀t)/k over the same flight time, and explain the ordering of the three.
(d) Now there is a 4 m/s wind from the north and the drone flies east. Without any new simulation, say where the package lands and why (relative motion).

Answer form: four numbers; four numbers and a sentence; three distances with a two-sentence explanation; one sentence.

## 4. (Numerical, uses your code: designing the drop experiment)

The experiment of part 4 drops the ping-pong ball from h = 0.6, 0.9, 1.2, 1.5, 1.8 m and fits v_t. Your camera runs at f = 30, 60 or 120 fps (choose the two you have access to, or 30 and 60).

(a) For each frame rate, tabulate the drag excess t − t₀ and the normalized error against the vacuum model at each height, taking the measurement uncertainty as the frame term 0.41/f alone (the most optimistic case). At which heights is the vacuum model rejected (E_n > 2)?
(b) Using `make_synthetic_drop_data` (SYNTHETIC: generated from the drag model, for design only) with seeds 0–19, n = 5 and n = 10 drops per height, at each frame rate: the mean and the standard deviation over seeds of the fitted v_t, and the mean reported uncertainty u_v. Does the reported u_v agree with the scatter over seeds? Convert to a relative uncertainty on C_D.
(c) Choose a frame rate and a number of drops that give C_D to ±20 % and justify the choice in three sentences, including why doubling the number of drops helps less than doubling the frame rate.
(d) Suppose your release adds a constant bias of +(s) ms to every timing (`release_bias_s` in the generator). Rerun (b) at 60 fps, n = 10. What happens to the fitted v_t, and how would the bean-bag control reveal the bias?

Answer form: a table per frame rate; a table; three sentences; two numbers and two sentences.

## 5. (Relative motion: a UAV survey square)

A survey UAV flies at a true airspeed of 25 m/s around a square of 2 km sides with tracks 000°, 090°, 180°, 270°, in a wind of W = (5 + s/5) m/s from 040°.

(a) For each leg: the heading to fly, the ground speed and the leg time (`heading_for_track`). Total time around the square against 320 s in still air.
(b) Prove that for a straight out-and-back flight along the wind, the round trip takes longer than in still air by the factor 1/(1 − W²/V²), and that the same holds for an out-and-back across the wind with factor 1/√(1 − W²/V²). Which is worse?
(c) At what wind speed does the UAV become unable to hold the 270° leg, and what does `heading_for_track` do there?
(d) The UAV's battery gives 20 minutes of flight. How many squares can it fly in this wind, and in still air?

Answer form: a four-row table and a total; two derivations and a comparison; a number and a sentence; two integers.

## 6. (Judgment: where the model stops)

For each situation, say which assumption of the lesson's model fails first, what the sign of the error in the predicted range or fall time would be, and what single measurement would show it, in one or two sentences each.

(a) A ping-pong ball tossed with heavy backspin.
(b) A golf ball driven at 70 m/s.
(c) A feather dropped from 2 m.
(d) A ping-pong ball dropped in a stairwell with the front door open.
(e) The same experiment repeated in a mountain town at 2500 m altitude with the sea-level ρ in the code.
(f) A steel ball bearing (8 mm) dropped from 2 m: can your camera distinguish its fall from free fall, and should you bother modelling the drag?

Answer form: six labelled items.
