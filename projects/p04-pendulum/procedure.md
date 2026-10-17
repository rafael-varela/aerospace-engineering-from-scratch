# Procedure: Period of a string pendulum against length and amplitude

**Engineering question.** Does a string pendulum's period follow T₀ = 2π√(L/g) within the uncertainty of a phone-video measurement at lengths from about 0.4 m to 1.5 m, and at what amplitude does the lengthening predicted by the exact (elliptic-integral) period become visible?

This is the bench version of a question every structural engineer asks of a lightly damped mode: what is its natural frequency, how does it depend on the geometry, and where does the linear model stop being right. Read `safety.md` first. If you have no kit, run every analysis step on `make_synthetic_pendulum_data(seed=<your seed>)` from the lesson code and label every result **SYNTHETIC**.

## 0. Predict first

In `data/prediction.md`, before any timing, dated:

1. Choose five or six lengths between 0.4 m and the longest your support allows (never above 2.0 m). Write the targets.
2. For each, compute T₀ = 2π√(L/g) with your local g if you know it (write the source) or 9.81 ± 0.03 m/s², and propagate u_L (typically 2 mm, mostly from locating the centre of the bob) and u_g into u_T₀.
3. Write your timing uncertainty before measuring: for video at frame rate f, one timing carries 0.41/f (P00.15); timing N swings at once divides it by N. Choose N (10 is a good start) and write the frame term per period.
4. At L ≈ 1.0 m, compute the exact period at 10°, 20° and 30° with `pendulum_period_exact`, and state for each whether the excess over T₀ is larger than twice your frame term per period. That is your prediction of what you will and will not detect.
5. Freeze the file. The comparison in step 5 is against these numbers.

## 1. Build

Minimum-cost tier: clamp or weight a stick on a sturdy table so 0.3 m overhangs, pinch the string in a binder clip at the end of the stick (the clip's jaws are the pivot), and tie the bob on. Tape the printed protractor behind the pivot so its 0° line hangs vertical (use the resting string as the plumb line). Mark the ends of the swing arc plus 0.5 m on the floor. Put the phone on its stand 1.5–2 m away, outside the arc, lens level with the bob, the whole arc and the protractor in frame.

Recommended tier: add the ruler or caliper for the bob's diameter, the scale for its mass, and the slow-motion mode. Enhanced tier: mount the break-beam at the bottom of the arc and log each interruption time over USB serial; the procedure is otherwise the same.

## 2. Calibrate

1. **Frame rate.** Film a running stopwatch for 10 s at each frame rate you will use; count frames; record the measured rate in `data/calibration.csv`.
2. **Length.** Measure from the clip's jaws to the top of the bob and add half the bob's height (or measure directly to its centre). Three readings, with the bob hanging. Record all three.
3. **Frame rule.** Decide which event you time: the frame in which the string passes the protractor's 0° line moving left to right. Write the rule down and use it every time; the passage through the vertical is sharp, the turning point is not.

## 3. Measure

Envelope: L ≤ 2.0 m, bob ≤ 0.2 kg, amplitude ≤ 45°. Check the stop criteria before each run.

1. **Length series.** At each length, release from rest at about 10° (read the protractor), wait one swing, then film at least N + 2 swings. Five runs per length. Shuffle the order of lengths if you can.
2. **Amplitude series.** At L ≈ 1.0 m, five runs each at 10°, 20° and 30°. For each run, read the amplitude at the start and at the end of the timed swings from the video.
3. **Raw data.** Record every run in `data/raw/YYYY-MM-DD_<tier>_timings.csv` (format in `data/README.md`): start and end frames, not just their difference; note any doubtful run in the `note` column instead of deleting it.

## 4. Analyze

Reduce the data with the lesson's code (`curriculum/phase-04/lesson-06-oscillations-waves`):

- For each length: T ± u_T with `period_from_timings(times, N, u_single=0.41/f)`.
- The amplitude correction from the mean of the start and end amplitudes; corrected T₀ = T / `amplitude_correction(θ_mean)`.
- g ± u_g at each length with `g_from_pendulum`, and their weighted mean.
- E_n against your reference g: |g − g_ref|/√(u_g² + u_ref²).
- `fit_power_law` on T against L: the exponent p (model: 0.5) and g from the coefficient.
- Plots: T against L with k = 2 error bars and the 2π√(L/g) curve; T/T₀ against mean amplitude at L ≈ 1 m with the exact curve.

## 5. Compare and explain

Where E_n exceeds 2, or p differs from 0.5 by more than its uncertainty, list candidate causes and rank them: amplitude correction missing or wrong amplitude; length measured to the wrong point; a frame rule applied inconsistently; a flexing or rocking support; a stretchy string; the true frame rate different from the setting. Use the pattern table in the lesson's part 5. Test at least one cause with a follow-up run (re-measure L to the centre; repeat one length at 5°; repeat with a stiffer support) and say what it showed.

For the amplitude series: did you detect the effects you predicted in step 0, and not the ones you predicted you would not? If a predicted detection failed, find out why before claiming anything.

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, calibration file, raw timing files, analysis script, both plots, the g table with E_n, the amplitude comparison, ranked causes with the follow-up result, the completed safety box, limitations, and what you would redesign.
