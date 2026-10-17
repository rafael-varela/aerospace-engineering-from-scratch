# Procedure: Frequency and damping of a string pendulum, identified from a logged swing

**Engineering question.** From a two-minute logged swing of the `p04-pendulum` string pendulum, what are its small-angle natural frequency and its damping ratio, each with an uncertainty you can defend; does the frequency fall with amplitude as the exact pendulum law predicts; and is the damping viscous, quadratic (air drag) or Coulomb (pivot friction)?

This is the bench version of a ground vibration test: excite a structure, record its free decay, and identify its modal frequency and damping. Read `safety.md` first. If you have no kit, run every analysis step on `data/SYNTHETIC_seed619_swing.csv` (or on `make_synthetic_swing(seed=<your seed>)` from the lesson code), label every result **SYNTHETIC**, and say so in the report.

## 0. Predict first

In `data/prediction.md`, before the first swing, **dated**:

1. Measure L from the pivot (the clip's jaws) to the bob's centre, three readings, and write L ± u_L. Write the g you will use and its source (your local value, or 9.81 ± 0.01 m/s²).
2. Predict the small-angle frequency f₀ = (1/2π)√(g/L) ± u with `predict_frequency`, and the frequency at your release amplitude (20° is the default) and at half of it, from the exact law (`amplitude_frequency_ratio`). Is the difference larger than the uncertainty you expect from 60 crossings (about 10⁻⁴ relative)? Write detectable / not detectable.
3. Predict an order of magnitude for the damping: weigh the bob, measure its diameter, and compute the air-drag coefficient b₂ with `quadratic_drag_coefficient` (C_d ≈ 0.47 for a sphere; say what you assumed for another shape) and the equivalent damping ratio at your release amplitude with `equivalent_zeta`. Write it as a range (a factor of 3 either way is honest; the string and pivot add to it).
4. Predict how far the amplitude will fall in 120 s, from the predicted ζ: A(120 s) ≈ A₀ e^{−ζω₀·120}.
5. Freeze the file. The comparison in step 5 is against these numbers.

## 1. Build

Rig the `p04-pendulum` support, string and bob (see `projects/p04-pendulum/procedure.md`, step 1), with L between 0.5 m and 1.5 m (never more than 2.0 m). Mark the arc's ends plus 0.5 m on the floor.

- **Minimum-cost (video).** Stick a high-contrast marker on the bob. Tape the printed protractor or grid flat behind the pivot, with its 0° line on the resting string. Put the phone on its stand 1.5–2 m away, square to the swing plane, lens level with the bob's lowest point, the whole arc and the reference in frame.
- **Recommended (phone gyroscope).** Weigh the phone in its padded pouch with the clip: **total ≤ 0.2 kg, or do not use the phone as the bob.** Hang the pouch from two points so the phone's screen faces the swing direction and it cannot spin; L runs to the phone's centre. Start the sensor-logging app on the gyroscope at its highest rate. Or: film a slow-motion video of the marker.
- **Enhanced (IMU).** Tape the IMU breakout flat to the bob, with one gyro axis along the pivot axis; run four light wires along the string to the `p04-electronics-bench` board; log angular rate over USB serial at 100 Hz or more. Weigh bob + IMU + tape.

## 2. Calibrate

1. **Time base.** Video: film a running stopwatch for 10 s and count frames (the p04-pendulum check); use the timestamps your tracking tool exports, not the nominal rate. Gyro: log for 60 s against a stopwatch and compare elapsed times.
2. **Length.** Three readings of L to the bob's (or phone's) centre with the bob hanging.
3. **Angle reference (video).** Photograph the protractor or grid with a known angle (a ruler held at a measured slope); the angle-reading uncertainty of your method is the scatter of 10 repeated readings of one still frame.
4. **Gyro bias (rate tiers).** Log 30 s with the pendulum at rest; the mean is the bias, the standard deviation the noise.
5. Record all of it in `data/calibration.csv`.

## 3. Measure

Envelope: L ≤ 2.0 m, bob ≤ 0.2 kg, amplitude ≤ 45°, released from rest by hand, no powered drive. Check the stop criteria before each run.

1. Start the recording, then release the bob from rest at about 20° (read the protractor), standing beside the plane of the swing. Do not touch anything for at least 120 s (longer is better: 300 s separates the damping mechanisms far better). Stop the swing by catching the string from the side.
2. **Repeat runs:** at least three runs at 20°, so the scatter of repeat runs can be compared with the fit's covariance.
3. **Amplitude runs:** one run released at 40° (still ≤ 45°) and one at 10°. A wide amplitude range is what tells viscous from quadratic damping, and it shows the frequency change.
4. Save each record untouched as `data/raw/YYYY-MM-DD_<tier>_run<k>.csv` (format in `data/README.md`), with the release amplitude and anything unusual in `data/raw/runs.csv`. Keep the videos and app exports; they are raw data too.

## 4. Analyze

Reduce every run with the lesson's code (`curriculum/phase-06/lesson-19-home-lab-pendulum-system-identification`):

- Damped frequency from zero crossings (`zero_crossings`, `period_from_crossings`) and from the FFT peak (`fft_peak_frequency`); they should agree.
- Peaks (`find_peaks_interp`); damping by the logarithmic decrement over n cycles (`log_decrement`) and by the ln(peak) line (`log_decrement_fit`).
- Nonlinear least squares of the damped cosine (`fit_damped_cosine`): A, ζ, ωₙ, φ, c with covariance uncertainties; the residual RMS compared with your angle-reading (or gyro) noise.
- Uncertainty checks: `bootstrap_fit` with block = 1 and with blocks of about 2 s; the scatter of ζ and ωₙ across your repeat runs.
- Diagnosis: `fit_envelope` (the exponent p) and `fit_envelope_vq` (b₁ and b₂, with their correlation) on each run, and on the 40° run in particular; `frequency_vs_amplitude` against the exact law.
- Small-angle frequency: the NLS ωₙ divided by the mean of `amplitude_frequency_ratio` over the peaks, and the intercept of the windowed frequency against A².

Rate data (gyroscope): the angular rate of a damped cosine is itself a damped cosine with a shifted phase, so the same fit applies directly; c is then the gyro bias.

## 5. Compare and explain

Compare the measured f₀ with the prediction by the normalised error E_n = |f_meas − f_pred|/√(U_meas² + U_pred²) with k = 2. Compare ζ with your predicted range. Where they disagree, rank the candidate causes and test at least one: L measured to the wrong point (E_n large in frequency, see `P04.06`); a support that flexes or rocks (the frequency is low and the damping high); the string's own stretch; a logging cable adding stiffness and damping (repeat without it, on video); the phone's pouch swinging relative to the string (a double-pendulum mode in the gyro spectrum); air currents. For the damping mechanism: does the exponent p sit near 0, 1 or 2, or between? Does the 40° run give a different ζ from the 10° run, as quadratic drag predicts? If the residual RMS is much larger than your noise, say which model error causes it (amplitude-dependent frequency and damping) and how you allowed for it in the uncertainty. "Experimental error" is not an explanation.

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, calibration file, raw records, analysis script, the plots (data with fit, decay envelope with the power-law and exponential fits, frequency against amplitude with the exact law), the results table with all uncertainty estimates, the damping diagnosis, the completed safety box, limitations, and what you would redesign. If any part used SYNTHETIC data, say which and give the seed in every caption.
