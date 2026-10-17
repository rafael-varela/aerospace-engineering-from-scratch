# Procedure: Drop and toss light objects, time them on video, and measure a drag coefficient

**Engineering question.** Does the quadratic-drag model F_D = ½ρC_D A v² with C_D = 0.5 predict the fall time of a ping-pong ball from 0.6–2.0 m and the range of a gentle toss to within the uncertainty of a phone-video measurement, and what value of C_D, with what uncertainty, do the drop times give?

This project repeats the method of `p00-first-experiment` with one new ingredient: a second object light enough for drag to be visible, and a control object dense enough for drag to be invisible. The difference between them, released the same way and timed by the same rule, is the drag. The analysis code is `curriculum/phase-04/lesson-01-kinematics-newtons-laws`.

## 0. Predict first

In `data/prediction.md`, before any drop or toss, dated:

1. Weigh the ball (0.1 g scale; or state the nominal 2.7 g ± 0.1 g) and measure its diameter (callipers, or a ruler between two books; 0.5 mm). Weigh and measure the control object. Compute today's air density from pressure and temperature, ρ = p/(287.05 T) with p in Pa and T in K (sea level, 20 °C: 1.20 kg/m³); if you cannot, use 1.20 ± 0.05.
2. Compute k = ρC_D A/(2m) and v_t = √(g/k) for both objects with C_D = 0.5 ± 0.1 for the ball and 1.0 ± 0.5 for the bag (a lumpy soft object; its C_D is a guess, and it does not matter, since the bag is a control).
3. Choose five drop heights between 0.6 and 2.0 m (the lesson uses 0.6, 0.9, 1.2, 1.5, 1.8 m). For each, tabulate the vacuum time √(2h/g) and the drag time from `vertical_drop_time` for both objects, with the uncertainty propagated (finite differences, as in P00.14) from u_h = 5 mm, u_g = 0.03 m/s², u_CD = 0.1, u_ρ = 0.05 kg/m³, u_m = 0.1 g, u_d = 0.5 mm.
4. Write your frame rate and the Type B frame term 0.41/f, and the normalized error you expect at each height between the drag model and the vacuum model if the frame term were the only uncertainty. Decide from this table, and from Problem 4 of the lesson, your frame rate and number of drops; write the decision and the reason.
5. For the toss: predict the range band of `simulate_projectile` for launch speeds 5–7 m/s and angles 40–50° with your k, and state the expected effect of spin (longer with backspin, shorter with topspin). This is a prediction to compare with, not a measurement of C_D.

## 1. Set up

Fix the phone on its stand, side-on, level, 1.5–2 m from the drop line with the whole fall in frame (drops) and later 1.5 m to the side of the toss lane with the whole lane in frame (tosses). Put the reference marks on the wall (every 0.3 m, measured from the floor) and on the floor along the lane (every 0.5 m from the release mark). Set the camera to its highest frame rate and record the setting. Lay the mat. Clear the circle and the lane. Close windows and doors. Read the safety box once more.

## 2. Calibrate

Film a running stopwatch display for 10 s and count frames: record the measured frame rate in `data/calibration.csv`. Check the wall marks against the tape (three readings each). Decide and write your frame rule: release = first frame in which a gap is visible between fingers and object; impact = first frame in which the object touches the mat. Apply the same rule to both objects, every time.

## 3. Drop

For each height, alternate the two objects: ball, bag, ball, bag, … until each has at least 8 drops (10 if time allows). Release from rest by opening two fingers from below, or by opening the clip, at the marked height measured to the *bottom* of the object. Record `height_m, object, drop, video_file, release_frame, impact_frame, fps_nominal, note` in `data/raw/YYYY-MM-DD_<tier>_drops.csv` as you go; never edit it. A drop you think you fumbled is recorded with a note, not deleted.

## 4. Toss

From the release mark, with the ball resting on the open fingers (least spin), toss gently underarm at about 45° and about 6 m/s so that the ball lands on the mat 2–3 m away. Ten tosses. From the video record, for each toss: the release frame; the frame and the position (against the floor marks) of the apex; the landing frame and position; and the position two frames after release (for the launch speed and angle). Record in `data/raw/YYYY-MM-DD_<tier>_tosses.csv` with the columns of `data/README.md`. If you mark the ball with a pen cross, note in the `spin` column whether you can see it rotating and roughly how fast (turns per 10 frames).

## 5. Analyze

With the lesson's code:

- Per height and object: `reduce_timings` on the frame-count times with the frame term as `u_single`; E_n against the vacuum time and against the drag time (with the predicted uncertainty from step 0).
- `fit_terminal_speed(heights, means, u_means)` on the ball; convert to C_D with `drag_coefficient_from_terminal_speed` and the relative uncertainty u_CD/C_D = 2u_v/v_t. Report the fit residuals at each height.
- The same fit on the bag, as a control: a sensible result is a large v_t with a large uncertainty (undetermined); a small v_t for the bag means a bias in the release or the frame rule that affects both objects.
- Per toss: launch speed and angle from the two-frame displacement (uncertainty from the 1 cm reading over 2/f seconds); `simulate_projectile` with your fitted k; E_n of the measured range against the predicted one with the launch-speed uncertainty propagated through the simulation (run it at v₀ ± u).
- The plot: fall time against height for both objects with k = 2 error bars, the vacuum curve and the two drag curves, as in `projectile_drag.png` with your data in place of the SYNTHETIC points.

## 6. Compare and explain

Where E_n > 2, rank the causes and test one. The patterns of P00.15 apply (a constant offset at every height: release; an error shrinking with height: the height reference) with two new ones: the control object above the vacuum curve by more than a frame (a release bias, or a bag with more drag than assumed: tell them apart by the height dependence, which is flat for a bias and growing for drag); tosses consistently long or short at the measured launch state (spin: confirm from the pen mark). If the fitted C_D differs from 0.5 by more than 2u, say whether the sphere-drag literature's range (0.4–0.5 at these Reynolds numbers) covers it before suspecting the apparatus.

If E_n is comfortably below 1 everywhere, say which term of your uncertainty budget was too pessimistic.

## 7. Report

Fill `report-template.md`. Submit: the prediction sheet, the calibration, the raw frame tables for drops and tosses, the reduced tables, the fit with residuals and C_D ± u, the plot, the ranked causes with the follow-up result, limitations, and what you would redesign.
