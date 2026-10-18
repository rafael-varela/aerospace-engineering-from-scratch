# Procedure: Lift and drag of a foam airfoil section on a guarded-fan bench rig

**Engineering question.** At the bench rig's Reynolds number (about 3 × 10⁴), how do the lift and drag coefficients of a NACA 0012 section compare with thin-airfoil theory, your own panel code with and without a boundary-layer estimate, and the published NACA data at Reynolds numbers of millions (Ladson, NASA TM-4074), and can you explain the gap?

This is the bench version of a wind-tunnel test: calibrate the instruments, write the prediction down, measure with tares and repeats, reduce with uncertainties, and compare. Read `safety.md` first, every session. If you have no kit, run every analysis step on the SYNTHETIC files in `data/synthetic/` (or on `make_synthetic_bench_data(seed=<your seed>)` from the lesson code) and label every result **SYNTHETIC**.

## 0. Predict first

In `data/prediction.md`, before the fan is ever switched on with the model in place, dated:

1. Measure the model: chord c and span s between the end plates (calipers or ruler, three readings each), thickness at the 30 % chord point, end-plate size. Record the instrument and its resolution.
2. Measure the jet: the straightener exit size (jet area), and decide the model's position (centred, about 50–100 mm from the straightener exit). Compute the blockage estimate A_frontal/(4A_jet) at 0° and at your largest angle.
3. Record the room temperature, pressure and humidity; compute ρ and μ (`air_density`, `air_viscosity`).
4. Predict the airspeed: from the fan's datasheet if it gives one, or from one anemometer reading taken with **no model** in place (the fan unplugged while you set the anemometer). Write V with its uncertainty, q, and the manometer reading you expect at your incline (Δl = q/((ρ_w − ρ)g sin θ)).
5. Compute Re = ρVc/μ and write it next to Ladson's lowest, 2.0 × 10⁶.
6. For each planned α (−4° to 14° in 2° steps): thin-airfoil cl, panel cl, the finite-span CL band for an effective aspect ratio you choose and justify (the lesson uses 4–6 for its end-plated model), the CD band from the low-Re estimate, and the scale changes they imply in grams on each scale (`predict_bench`). Mark which predicted changes are smaller than five times your scale's resolution: those points will be dominated by resolution.
7. Write one sentence on what you expect the bench data to show relative to the predictions, and why (Reynolds number). Freeze the file. The comparison in step 5 is against these numbers.

## 1. Build (fan unplugged throughout)

1. Shape the model: two NACA 0012 card templates at 100 mm chord (from the `P10.02` generator) glued to the ends of a foam blank; sand down to the templates. Mark the chord line and the quarter chord on the end plates. Craft knife and sanding only with the fan unplugged; hot-wire cutting only with an adult and ventilation.
2. Glue the end plates square to the span. Drill a small hole at the quarter chord through one end plate for the sting pin.
3. Build the straightener: straws cut to equal length, packed tightly in the square duct at the fan end and taped in; leave a settling length of at least one straw length after them.
4. Lift balance: the sting stands on a heavy base on the lift scale; the model sits on the sting at the quarter chord with a pin, so α can be set and locked. The protractor template is fixed to the sting with its zero along the chord line at α = 0 (set with a ruler square to the jet axis).
5. Drag balance: the sting is mounted on the vertical arm of the bell-crank, which pivots on a nail; the horizontal arm rests on the second scale (or the same scale in a second session). Measure l_D (pivot to the model's quarter chord) and l_s (pivot to the scale pad); the nominal ratio is l_D/l_s.
6. Airspeed: the pitot-static tube in the jet beside the model position (or in its place for a no-model run), the inclined manometer below and away from the fan.
7. Check every safety control in `safety.md`, then plug in.

## 2. Calibrate

1. **Lift scale.** Place known masses (0, 20, 50 g) on the model base; record readings. The slope should be 1.00; record any offset.
2. **Pivot arm.** With the fan unplugged, tie a thread to the model at the quarter chord, run it horizontally downstream over the pencil pulley, and hang known masses (0, 1, 2, 5, 10 g, then down again). Record the drag-scale reading each time. The slope is the arm ratio R (`pivot_calibration`); the up and down readings show friction.
3. **Manometer.** Measure the incline as rise over run along the board (±0.5 mm on each); compute θ and its uncertainty. Read the zero with the fan off, at eye level square to the scale.
4. **Airspeed check.** Fan on, no model: read the manometer (or anemometer) at least three times, with a minute between readings.

## 3. Measure

Envelope: the fan at its chosen speed setting only; nothing changed while it runs. **Unplug before every change of angle.**

1. **Tare runs.** (a) Sting and balance without the model: fan off reading, fan on reading, for both scales, three times (the sting's own drag and lift). (b) Model at α = 0: three fan-off/fan-on pairs (the zero-lift check; a non-zero lift here is flow angularity or misalignment).
2. **The α sweep.** For each α from −4° to 14° in 2° steps: unplug, set α, lock it, photograph the protractor; plug in, fan off: read both scales; fan on: wait 20 s, read both scales (photograph the displays; if they fluctuate, record the middle and the range); fan off. Then repeat the whole sweep in the opposite order (from 14° down to −4°).
3. **Airspeed during the sweep.** Read the manometer at the start, middle and end of each sweep.
4. **Raw data.** Record every reading in `data/raw/YYYY-MM-DD_<tier>_*.csv` (formats in `data/README.md`), with photo file names; note any doubtful reading instead of deleting it.

## 4. Analyze

Reduce with the lesson's code (`curriculum/phase-10/lesson-23-home-lab-airfoil-bench`), `load_bench_dataset(folder, prefix=...)` then `reduce_bench`:

- ρ, μ; q and V from the manometer (mean, with the uncertainty from the scatter, the reading resolution and the incline); the drift between the start and the end.
- R from the calibration; the tares.
- Lift and drag per reading pair; mean and uncertainty at each α (two runs); CL and CD with their uncertainty, including the blockage band on q.
- The slope dCL/dα over −4° to 4°, the zero-lift angle (flow angularity), CL_max and the angle where it occurs.
- Plots: CL and CD against α (corrected for the flow angle) with k = 2 bars, the prediction bands, the panel curve and Ladson's data at 2.0 × 10⁶.

## 5. Compare and explain

1. E_n of the measured lift slope against your predicted band (`en_number`). Expect disagreement with inviscid theory; explain it with the Reynolds-number argument (the Thwaites/Michel comparison of the lesson: at Re ≈ 3 × 10⁴ the laminar boundary layer separates before it can become turbulent).
2. CL and CD against Ladson's at the same α: ratios, and the gap in CL_max and the stall angle.
3. Rank the alternative explanations for any remaining disagreement: blockage and jet non-uniformity, flow angularity, the end plates' leakage (finite span), balance friction, fan unsteadiness, model accuracy (thickness, leading-edge shape), sting interference. Test at least one: repeat a point with the model moved 30 mm off-centre (non-uniformity), or reverse the model to measure the flow angle (a symmetric section at +α and −α).

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, calibration records, raw data files, analysis script, plots, the CL/CD table with uncertainties, E_n values, the ranked explanations with the follow-up result, the completed safety box, limitations, and what you would redesign.
