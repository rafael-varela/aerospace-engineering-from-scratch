# P00.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer. Where a problem says "your data", use the drops you recorded; otherwise use the seeded synthetic generator and label the result as synthetic.

## 1. (Warm-up: one prediction)

h = 1.500 ± 0.010 m, g = 9.81 ± 0.03 m/s² (standard uncertainties).

Find t_pred and u_pred, report them correctly, and say which input dominates the budget. Then: by how much would you have to improve u_h for the two inputs to contribute equally?

## 2. (Derivation: a hand that pushes)

The ball leaves the hand with a small downward speed v₀ instead of from rest.

(a) Starting from h = v₀t + ½gt², derive t(v₀, h).
(b) Show that dt/dv₀ at v₀ = 0 equals −1/g, independent of h, and compute the time shift for v₀ = 0.1 m/s.
(c) Drag lengthens the fall by Δt_drag ≈ t₀ gh/(6v_t²). Show how Δt_drag and the push shift Δt_push scale with h, and explain how five heights let you tell them apart even when both are below the frame resolution of a single drop.

## 3. (Unfamiliar scenario: the same experiment on Mars)

A rover team repeats your drop experiment on Mars: g = 3.71 m/s², air density 0.020 kg/m³, same tennis ball (m = 0.058 kg, d = 67 mm, C_D = 0.5), same 30 fps camera, same heights.

(a) Compute t_pred at 0.6 m and 1.8 m.
(b) Compute the terminal speed and the first-order drag correction at 1.8 m. Compare with the Earth value from the lesson.
(c) The frame term u_frame is the same in seconds. Is the Mars comparison sharper or blunter in *relative* terms than on Earth, and why?
(d) An engineer proposes a 1.00 s fall on Earth to make the timing easier. What height does that need, and is it permitted under the `P00.13` envelope?

## 4. (Numerical: uses your code)

Use your `data/raw/drops.csv`, or `make_synthetic_timings(HEIGHTS, n_rep=6, fps=30, seed=SEED)` with your own seed.

(a) Run `compare` with u_h = 0.005 m and the frame term from `frame_timing(1, fps)`. Report the table (h, t_pred ± u_pred, t_meas ± u_meas, E_n) and the plot.
(b) Repeat at 120 fps (re-film, or regenerate with `fps=120`). How do u_meas and E_n change? Which part of u_meas changed and which did not?
(c) Simulate a stopwatch: subtract 0.15 s from every timing (a late start) and recompute. Describe the E_n pattern and explain why more drops would not help.
(d) Run `fit_g` on your t_meas and report g ± u_g. Is it consistent with 9.81 m/s²?

## 5. (Drag by design: the ping-pong ball)

A ping-pong ball has m = 2.7 g and d = 40 mm; take C_D = 0.5 and ρ = 1.2 kg/m³.

(a) Compute its terminal speed.
(b) Compute the drag correction to the fall time at 0.6, 1.2 and 1.8 m, first with the first-order formula and then with `fall_time_with_drag`. How good is the first-order formula at 1.8 m?
(c) With u_meas ≈ 15 ms at 30 fps and ≈ 4.5 ms at 120 fps, compute E_n of the *drag-free* prediction against a measurement that follows the drag model exactly, at each height and frame rate. At which frame rate and heights is drag detectable (E_n > 2)?
(d) Design a side-by-side drop of the tennis ball and the ping-pong ball from 1.8 m that makes the drag effect visible in a single video. Predict the arrival gap in milliseconds and in frames at 120 fps.

## 6. (Judgment: what E_n does and does not tell you)

(a) A learner reports E_n < 2 at all five heights, with u_meas = 50 ms from a stopwatch. What has been learned about the model, and what has not? What single change would make the test more informative?
(b) Another learner reports E_n = 0.3, 0.4, 0.5, 1.1, 2.6 at 0.6 to 1.8 m, with t_meas longer than t_pred at every height. Rank three candidate causes, say what pattern each would produce, and name the one follow-up measurement you would make first.
(c) A third reports E_n = 3.1, 1.9, 1.2, 0.8, 0.6, with t_meas *shorter* than t_pred and the discrepancy largest at 0.6 m. What is the most likely cause, and what is the fix?

## 7. (Where is h = 0? The ball's radius)

A learner measures each height to the *centre* of the tennis ball (radius r = 33.5 mm) but the model needs the distance the ball's *bottom* falls.

(a) Derive the ratio t_pred(wrong)/t_pred(right) as a function of h and r.
(b) Evaluate the error in milliseconds at 0.6 m and 1.8 m.
(c) Compute E_n at both heights for u_meas = 15 ms (30 fps) and 4.5 ms (120 fps). At which frame rate would this mistake be caught, and what is its signature across heights?
