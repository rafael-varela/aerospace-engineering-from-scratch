# Procedure: Discharge coefficient of a draining bottle, and airspeed from a water U-tube manometer

**Engineering question.** Does a gravity-driven jet from a small hole follow Torricelli's law with a discharge coefficient in the typical sharp-edged range (0.60–0.62), and how well can a water U-tube read the dynamic pressure of a hair dryer on cool and of a small guarded fan, upright and tilted? This is the bench version of two things every flow facility does: calibrate a flow element (an orifice) and measure airspeed from a pitot pressure.

Read `safety.md` first. If you have no kit, run every analysis step on the SYNTHETIC files in `data/synthetic/` (or on `make_synthetic_flow_data(seed=<your seed>)` from the lesson code) and label every result **SYNTHETIC**. Synthetic data are not measurements.

## 0. Predict first

In `data/prediction.md`, before any water flows or any fan runs, dated:

1. **Bottle.** Measure the inner diameter D of the straight part (outer diameter with a caliper minus twice the wall, or the circumference with a paper strip divided by π, minus twice the wall) and the hole diameter d (caliper inside jaws or the largest drill shank that enters). Write each with its standard uncertainty and the instrument. Mark the hole centre as the zero of the ruler taped to the bottle.
2. **Drain-time prediction.** With C_d = 0.61 ± 0.03 (typical for a sharp-edged round orifice; a hole in a thin plastic wall is close to that, but not exactly), predict the time for the level to fall from 200 mm to 30 mm above the hole centre with `predict_drain_time`, and the budget term by term. Say which term dominates. For the lesson's bottle (D = 90.0 ± 0.5 mm, d = 4.0 ± 0.1 mm) the prediction is **102.7 ± 7.3 s**, with C_d and d contributing about 5 s each.
3. **Jet Reynolds number** at 200 mm and at 30 mm (`jet_reynolds`), and one sentence on whether viscosity should matter.
4. **Air density** from the room's pressure and temperature (`air_density`), with its uncertainty.
5. **Manometer resolution.** With your reading uncertainty per meniscus (0.3 mm from a zoomed photo of a 1 mm scale is a fair start; four readings make each difference), compute with `min_speed_for_uncertainty` the lowest airspeed you can read to 10 % upright and at your planned tilt. For 0.3 mm per reading: 7.0 m/s upright, 4.9 m/s at 30°, 2.9 m/s at 10°.
6. **Airspeed predictions.** For the fan: its rated volume flow divided by its swept annulus (`fan_mean_velocity`; a manufacturer's figure, measured in a test rig, and an upper bound for the speed a metre away). For the hair dryer: write a range (a household dryer's cool jet near the nozzle is typically 8–15 m/s) and the Δh that range would give upright and tilted. Then decide, from step 5, which source needs the tilted tube.
7. Freeze the file. The comparisons in step 5 are against these numbers.

## 1. Build

**Bottle.** Make the 4 mm hole 20–40 mm above the base before filling (safety box: cutting board, cut away from the holding hand). Deburr the edge with fine sandpaper so it is clean and square; do not flare it. Tape the ruler vertically on the far side of the bottle with its zero level with the hole centre. Seal the hole with a strip of tape folded back on itself so it pulls off in one movement. Stand the bottle in the sink or basin (on the floor if a basin), on a level surface, with the hole pointing into the basin and the jet's landing point inside it. Set the phone on its stand level with the 100 mm mark, about 300 mm away, with the ruler and the water surface in frame.

**Manometer.** Tape the clear tube into a U on the board, legs parallel and about 30 mm apart, with a ruler between them. Fill it with about 20 mL of room-temperature water with a syringe (or by dipping one end into a cup and lowering the U), so the menisci sit near the middle of the ruler; tap out any bubbles. Push the straw 10–20 mm into the end of the left leg and tape it airtight; cut its open end square. The right leg stays open to the room. Stand the board upright on a towel on the manometer table. For the tilted reading, lay the board on the ramp at 10–30° and measure the angle from rise over run (to about 0.5°).

**Air source.** Place the fan or the hair dryer on a separate dry table or hold it at the board's edge; never above the basin. Hair dryer: heat switch on cool before plugging in. Fan: grille intact, rating plate ≤ 30 W.

## 2. Calibrate

1. **Hole and bottle.** Measure d three times (rotate the gauge) and D at three heights; record means and scatter.
2. **Timing.** Film a phone stopwatch for 60 s and check the video's frame times against it; the frame interval (1/30 s) and your judgement of the crossing set u_t (0.05–0.1 s).
3. **Level reading.** With the hole sealed, read the level from ten photos, moving the phone away and back; the scatter is u_h (0.5 mm is typical with a meniscus and parallax).
4. **Manometer zero and reading.** With no flow, photograph the menisci ten times; the scatter is the reading uncertainty. The two zero readings differ by the meniscus and the tube's bore; they are subtracted, so they need not agree.
5. **Tilt.** Measure the board's angle twice and record the mean and half the difference as its uncertainty.
6. **"Calibrated pressure".** The U-tube is a primary instrument: Δp = (ρ_w − ρ_a) g Δh needs only the water density at room temperature, g, and a ruler. The tilted reading is calibrated against it: read the hair dryer's high setting both upright and tilted, at the same position, and check that L sin θ equals the upright Δh within its uncertainty.

## 3. Measure

Envelope: ≤ 2 L of water at room temperature in the session; hair dryer on cool only; fans guarded, ≤ 30 W; water and mains on separate tables; mop any spill before switching anything on. Stop at the stop criteria.

1. **Draining bottle.** Fill to about 205 mm above the hole with the measuring jug, keeping count of the water used. Start the video, pull the tape off in one movement, and let the bottle drain to below 30 mm. From the video, read the time at which the level crosses 200, 190, …, 30 mm. Repeat three times (refill from the basin if it stays clean, still within 2 L in total).
2. **Manometer, hair dryer.** With the dryer on cool, hold the straw's open end 50 mm in front of the nozzle on its axis, facing into the jet. Wait 5 s for the menisci to settle; photograph; take the straw out of the jet and photograph the zero. Do it twice for each setting, upright and tilted.
3. **Manometer, fan.** The same with the straw 300 mm in front of the fan's hub, on its axis, at each speed setting, upright and tilted.
4. **Enhanced tier.** Hold the vane anemometer at the same point as the straw and record its reading after the straw's.
5. **Raw data.** Record every reading in `data/raw/YYYY-MM-DD_<tier>_drain.csv` and `..._manometer.csv` (format in `data/README.md`) with photo or video file names; note any doubtful reading instead of deleting it.

## 4. Analyze

Reduce the data with the lesson's code (`curriculum/phase-08/lesson-19-home-lab-flow-pressure-measurement`):

- `fit_discharge_coefficient(t, h, u_h, u_t, d, u_d, D, u_D)` for each run: C_d ± u, the budget (slope, hole, bottle), χ²/dof, and the implied starting level. The mean and scatter of the three runs.
- E_n of C_d against 0.61 ± 0.03, and of the measured 200 → 30 mm time against your frozen prediction.
- For each manometer reading: L = (right − left) − (zero right − zero left) along the tube, Δh = L sin θ (`vertical_head`), Δp (`manometer_dp`), V ± u (`airspeed_with_uncertainty`).
- Reynolds numbers of the probe and the source at the measured speeds.
- Plots: √h against t with the prediction line; manometer airspeed against the anemometer (enhanced) or against the setting, with k = 2 error bars, upright and tilted.

## 5. Compare and explain

Where E_n exceeds 2, or χ² is far from the dof, list candidate causes and rank them: the hole diameter (it enters squared and is the hardest dimension to measure); a burr or a flared edge (changes C_d itself); the bottle not cylindrical over the range used; the ruler's zero not at the hole centre; the start of the jet disturbed by the tape; the straw not aligned with the jet (a pitot is insensitive to about ±15°, but not more); the straw outside the jet core; the jet's unsteadiness; the tilt angle's uncertainty at shallow angles. Test at least one cause with a follow-up (re-gauge the hole with drill shanks; redo one run with the ruler zero re-checked; move the straw 10 mm off axis and see the change) and say what it showed.

For the fan: did the upright tube resolve it, as predicted in step 0? If not, report "below resolution" with the upper bound, not a number.

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, the calibration file, the raw data files, the analysis script, both plots, the C_d table with E_n, the airspeed table with uncertainties and Reynolds numbers, ranked causes with the follow-up result, the completed safety box, limitations, and what you would redesign.
