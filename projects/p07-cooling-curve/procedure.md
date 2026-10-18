# Procedure: Cooling curve of a small aluminium bar

**Engineering question.** Does a small aluminium bar, warmed in tap-hot water no hotter than 60 °C and hung in still room air, cool with the time constant τ = ρcV/(hA) predicted from natural convection plus radiation, and is a lumped (uniform-temperature) model valid for it?

This is the bench version of a thermal-balance test: predict a unit's thermal response, measure it, and correct the model. Read `safety.md` first. If you have no kit, run every analysis step on the SYNTHETIC files in `data/synthetic/` (or on `make_synthetic_cooling_data(seed=<your seed>)` from the lesson code) and label every result **SYNTHETIC**.

## 0. Predict first

In `data/prediction.md`, before the bar goes into the water, dated:

1. Measure the bar: D at three places and L at three places (caliper if you have one; record the instrument and its resolution), and weigh it. Compare the mass with ρπD²L/4 (2700 kg/m³ for 6061): a mass more than 3 % off means a different alloy or a hollow bar.
2. Record the room temperature where the bar will hang.
3. Compute V, A, L_c = V/A, ρcV (Lienhard Table A.1 values from the lesson), h_conv from `h_conv_bar` and h_rad from `h_rad_linear` at θ = 14 K, the Biot number, and τ.
4. Compute u(τ) with `predict_tau_uncertainty`, with the budget term by term; say which term dominates.
5. Predict τ in each third of the run (2–20, 20–40, 40–60 min, at the mean excess in each) from `predict_tau`.
6. Freeze the file. The comparison in step 5 is against these numbers.

## 1. Build

Minimum-cost tier. Deburr the bar's ends. Build the stand: two equal stacks of books about 20 cm high, 15 cm apart, with a pencil or ruler across them, at the dry end of the table, on a towel. Tie two loops of thin cotton string from the pencil, 3 cm apart, so that the bar will hang horizontally about 10 cm below the pencil and well clear of the books and the table. Cut two squares of foil tape (or masking tape) about 15 × 15 mm. Put the mug on a tray at the opposite end of the table. Set the room thermometer 1 m from the stand at the same height. Switch on the probe thermometer and check it reads the room within 0.5 °C of the room thermometer.

Recommended tier. Build the thermistor divider and logging script of `p04-electronics-bench` (USB power from a laptop or power bank), with the board and laptop at the dry end; route the thermistor leads to the stand. Log one reading per second.

Enhanced tier. Wrap the second bar in one layer of black tape and hang it in a second sling 15 cm from the first, so both cool together in the same air. Have the cardboard draught shield ready.

## 2. Calibrate

1. **Thermometers against each other.** Put the probe (or thermistor in a small zip bag) and the room thermometer side by side for 10 minutes and record both every minute: the mean difference is an offset to correct, the scatter is your Type A reading uncertainty (if below the resolution divided by √12, use that).
2. **Probe response time (kitchen probe only).** Move the kitchen probe from room air into a mug of lukewarm tap water (about 35 °C) and read every 2 s. The thermistor never goes in the water: it and its leads stay dry, and its lag is left to the lagged fit. For the kitchen probe, the time to reach 63 % of the step is a rough τ_p in water; taped to the bar it will be different, and the lagged fit measures it.
3. Record everything in `data/calibration.csv`.

## 3. Measure

Envelope: water ≤ 60 °C, measured; bar handled with tongs; battery or USB power only; stop at the stop criteria.

1. Fill the mug with tap-hot water at the sink, carry it on the tray, and **measure it**. Above 60 °C, add cold water until it reads 60 °C or less. Record the reading.
2. With the tongs, lower the bar into the water. Leave it 5 minutes. Measure the water again (it will have cooled a little; record it).
3. Start reading the room thermometer (t = 0, 30, 60 min).
4. With the tongs, lift the bar out, dry it on the towel (tongs plus a folded cloth), and place it in the sling. Tape the probe tip (or thermistor) to the middle of its underside. **Start the clock when the probe touches the bar.**
5. Read every 15 s for 60 minutes (or let the logger run). Do not move, breathe on, or stand within 1 m of the bar; note anything that happens (a door, a person, the heating switching on) in the `note` column with its time.
6. At 60 min, record the room again. Leave the bar to cool fully before handling it.
7. Repeat the run on another day or at another time of day. Enhanced tier: one run with the draught shield over the stand, one run with both bars.

Record every reading in `data/raw/YYYY-MM-DD_<tier>_cooling.csv` and the room readings in `..._ambient.csv` (format in `data/README.md`). Never edit a raw file.

## 4. Analyze

Reduce the data with the lesson's code (`curriculum/phase-07/lesson-15-home-lab-cooling-curve`):

- Convert to kelvin; T∞ = the mean room reading (corrected for any offset found in calibration); u_T from calibration.
- `fit_loglinear(t, T, Tinf, u_T, t_min=120, theta_min=1.0)`: τ, u(τ), χ²/dof.
- `fit_exponential` with T∞ fixed and with T∞ free; `fit_lagged` on the whole record (τ_p).
- Windowed log-linear fits over 2–20, 20–40 and 40–60 min.
- `fit_variable_h` with ε = 0.09 (and 0.05, 0.13): f_conv.
- h = ρcV/(τA) with `h_from_tau`; E_n of τ and of h against the prediction (`en_number`).
- The Biot number with the measured h, and `short_cylinder` against `lumped_theta` at t = τ.
- Plots: T against t with the fits; ln θ against t with the log-linear fit and the predicted slope; the windowed τ against the mean excess with the predicted curve.

## 5. Compare and explain

Where E_n exceeds 2, or χ²/dof is far from 1, list candidate causes and rank them: draughts; T∞ wrong or drifting; the emittance; the probe's mass and stem; the sling touching too much of the bar; sunlight or a lamp; the alloy (last). Use the pattern table in the lesson's part 5. Test at least one cause with a follow-up (the draught shield; a second room thermometer; the black-tape bar; the thermistor instead of the kitchen probe) and say what it showed.

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, calibration file, raw data files, analysis script, the three plots, the τ, h and E_n table, the windowed τ, the Biot number and distributed comparison, ranked causes with the follow-up result, the completed safety box, limitations, and what you would redesign.
