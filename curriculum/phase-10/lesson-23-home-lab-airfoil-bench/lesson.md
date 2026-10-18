# P10.23  Home Lab: Airfoil Lift & Drag on a Bench Rig 🔬

> On a kitchen table a 100 mm wing at 4 m/s makes about ten grams of lift and one gram of drag at a Reynolds number seventy-five times below any NACA tunnel. Predict it, measure it with its uncertainty, and explain why it is not the textbook polar.

| | |
|---|---|
| **Learning outcomes** | 1. Build and calibrate a safe guarded-fan bench rig with a straw-honeycomb flow straightener, a kitchen-scale lift balance and a bell-crank drag balance, and measure the airspeed with an inclined water manometer on a pitot-static tube (or a vane anemometer). 2. Write, before measuring, a dated prediction of CL and CD at each angle, with an uncertainty band, from thin-airfoil theory, your panel code, a finite-span correction and a low-Reynolds-number drag estimate, including the scale readings it implies. 3. Reduce raw scale and manometer readings to CL and CD with tares, the calibrated arm ratio and first-order uncertainty propagation, and state the Reynolds number. 4. Compare measurement, computation and published data (Ladson's NACA 0012 at 2 × 10⁶) and explain the gap with the Reynolds-number argument (laminar separation before transition), separating it from rig effects (blockage, flow angularity, the fan's unsteadiness). |
| **Prerequisites** | `P10.01` (balance data to coefficients), `P10.02` (NACA generator), `P10.04` (thin-airfoil theory), `P10.05` (panel method), `P10.06` (Thwaites and Michel), `P10.09` (finite-wing lift slope), `P10.14` (low-Reynolds-number aerodynamics), `P08.19` (the flow and pressure home lab: pitot tube and manometer; recapped in §2.2), `P00.13` (safety), `P00.14` (measurement and uncertainty), `P00.15` (first experiment, normalised error). Skill check: you can propagate the uncertainty of C = F/(qS) and convert a ± half-width to a standard uncertainty. |
| **Engineering question** | A foam NACA 0012 section, chord 100 mm and span 200 mm between end plates, sits in the 4 m/s jet of a guarded desk fan (Re ≈ 2.6 × 10⁴). What lift and drag coefficients does it give against angle of attack, with what uncertainty, and how and why do they differ from thin-airfoil theory, your panel code and the NACA data at 2 × 10⁶? The prediction, before measuring: q = 9.6 Pa (0.98 mm of water), CL slope 0.073–0.082 per degree for an effective aspect ratio of 4–6 (inviscid), CD at least 0.020 (laminar friction). The lesson's SYNTHETIC dataset (not a measurement) gives a slope of 0.057 ± 0.002 per degree (55 % of Ladson's 0.104 at 2 × 10⁶), CL_max 0.53 ± 0.04 near 12°, CD 0.034 at zero lift (5.4 times Ladson's), and a flow angularity of 0.65°: the pattern low-Reynolds-number physics predicts, because the laminar boundary layer separates before it can become turbulent. |
| **Assumptions and limits** | Steady mean forces (the fan's flow is unsteady; readings are time-averaged by eye or by repeats); a model at the jet centre in a reasonably uniform core; the end plates give partial two-dimensionality (an effective aspect ratio between the geometric 2 and infinity, bracketed as 4–6); incompressible room air; small angles in the balance geometry; a blockage estimate from closed-section practice, carried as an uncertainty. Thin-airfoil and panel predictions are inviscid; the drag estimate ignores laminar separation bubbles. The rig cannot reach the Reynolds numbers of the NACA data, and that is the point of the comparison. |
| **Analytical work** | The manometer relation and its uncertainty (the cot θ term); the pitot-static principle; the balance equations with tares and the bell-crank ratio; C = F/(qS) and its absolute uncertainty; Reynolds number; the finite-span slope a₀/(1 + a₀/(πAR_eff)); the low-Re drag estimate; Thwaites' separation point (independent of Re) against Michel's transition (dependent on Re); the blockage estimate. |
| **Computational deliverable** | `starter/main.py`: `air_density`, `air_viscosity`, `manometer_q`, `airspeed`, `anemometer_q`, `reynolds_number`, `fit_line`, `pivot_calibration`, `spring_calibration`, `lift_from_scale`, `drag_from_pivot`, `drag_from_spring`, `coefficient`, `blockage_fraction`, `en_number`, `thin_airfoil_cl`, `finite_span_slope`, `low_re_cd`, `laminar_upper_surface`, `predict_bench`, `reduce_bench`, `compare_with_ladson` (the repository and Ladson loaders, the NACA generator, the compact panel method, the SYNTHETIC generator `make_synthetic_bench_data`, `write_synthetic_csv`, `read_bench_csv`, `load_bench_dataset` and `main` are given). Tests (16) prove: room-air density and viscosity reproduce the standard values and respond correctly to humidity and temperature; the manometer relation, its incline amplification and its uncertainty terms; airspeed and anemometer propagation; line fits and calibrations recover exact data and report honest slope uncertainties on seeded noise; the balance reduction inverts its forward model with tares; coefficient uncertainties stay finite at zero load; blockage grows with α; thin-airfoil, panel and finite-span predictions have the right limits; the drag estimate scales as Re^−½; Thwaites' separation point is independent of Re while transition is not; the prediction table is consistent; the Ladson comparison returns unity on Ladson's own data; the SYNTHETIC files are labelled, seeded and identical to the shipped ones, and reduce back to the generator's stated truth. |
| **Practical activity** | **physical experiment**: build the guarded-fan rig, write a dated prediction, calibrate the balances and the manometer, sweep α up and down with tares, and reduce to CL and CD with uncertainties; compare with theory, your panel code and the NACA data. See `projects/p10-airfoil-bench/`. A learner without a kit completes the lesson with the SYNTHETIC dataset in `projects/p10-airfoil-bench/data/synthetic/`, labelled as synthetic everywhere it appears. |
| **Equipment** | Minimum: a desk fan of at most 50 W with intact front and back grilles, drinking straws and a cardboard duct, XPS or EPS foam board, card templates and end plates, a skewer or dowel sting, a 1 g kitchen scale, a bell-crank arm, coins, a protractor, a ruler, a thermometer and a pressure source, safety glasses; a vane anemometer if you have one. Recommended: a 0.1 g pocket scale for lift, a 0.01 g jewellery scale for drag, a pitot-static tube and an inclined water manometer. Enhanced: a load cell with a microcontroller logger and a digital differential-pressure sensor. See `projects/p10-airfoil-bench/bom.yaml` (reuses items from `p00-calibration`, `p00-first-experiment`, `p04-drag-drop`, `p04-electronics-bench` and `p05-cantilever`). |
| **Safety** | See `projects/p10-airfoil-bench/safety.md`. Summary box at the start of part 4. |
| **Evidence** | `data/prediction.md` (dated before measuring), the calibration records, `data/raw/*.csv`, photographs of readings, the CL and CD table with uncertainties, the plots with the prediction bands, the panel curve and Ladson's data, E_n values, the ranked explanations of the gap with a follow-up test, `airfoil_bench.png`, and the report following `projects/p10-airfoil-bench/report-template.md`. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; all four competencies. Threshold: a dated prediction with uncertainty; the safety box completed; all public tests pass; CL and CD reported with uncertainty budgets and the Reynolds number; the comparison with theory, computation and the NACA data explains the gap with the Reynolds-number argument and tests at least one rig effect; anything synthetic is labelled SYNTHETIC. |
| **References** | `@nasa-glenn-beginners-guide-aeronautics` NASA Glenn "Build Your Own Wind Tunnel" (https://www.grc.nasa.gov/www/k-12/WindTunnel/build.html); "Wind Tunnel Parts"; "Open Return Wind Tunnel"; "Wind Tunnel Testing"; "Reynolds Number" (required-free): classroom tunnels and their parts, a model for this rig, not its procedure. `@leishman-intro-aerospace-flight-vehicles` "Wind Tunnel Testing" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/wind-tunnel-testing/): Components of a Low-Speed Wind Tunnel; Flow Quality; Freestream Speed Measurement; Wall Interference; Instrumentation & Measurement Techniques; Data Reduction (required-free). `@mit-ocw-16-100-aerodynamics` Lecture 35 "Drag Tare Due to Mount"; Lecture 34 "Three-Dimensional Wall Effects" (required-free): tares and wall interference. `@ntrs-nasa-tm-4074-ladson-naca0012` Ladson, NASA TM-4074 (required-free): the comparison data at 2–18.9 × 10⁶. `@ntrs-nasa-tm-4062-eppler-387-low-re` McGhee, Walker & Millard, NASA TM-4062: Introduction; Presentation of Results; Discussion of Results (required-free): low-Re airfoil data at 60 000–460 000 that bridge toward the bench. `@ntrs-naca-tr-824-summary-of-airfoil-data` (required-free): the high-Re section data. `@uiuc-selig-summary-low-speed-airfoil-data` Selig et al., Vol. 1, §2 "Wind-Tunnel Facility and Measurement Techniques" (supplementary-free, link only, all rights reserved): how a small university tunnel measures and corrects lift and drag. `@jcgm-100-2008` (required-free): uncertainty propagation as in `P00.14`. |
| **Extensions** | ★ Repeat the sweep at three fan speeds (Re ≈ 2.0, 2.7 and 3.3 × 10⁴) and plot the slope, CL_max and CD at zero lift against Re. ★ Add a trip strip (a thin thread or a strip of tape at 10 % chord on both surfaces) and see whether the lift slope rises and the drag falls, as tripping a laminar separation bubble can do at low Re. ★ Build the enhanced-tier load cell and pressure logger and measure the spectrum of the lift fluctuation: how much of your scatter is the fan? |
| **Time** | unmeasured |

## 1. The engineering problem

Every number in Phase 10's earlier lessons came from a model or a published tunnel. This lesson makes one. A wind-tunnel test has the same parts at any scale: a flow of known speed and quality, a model of known geometry, a balance that separates lift from drag, tares that remove what is not the model, calibrations that turn readings into forces, and a prediction written down beforehand so the result is a test and not a description. The bench rig is a guarded desk fan, a pack of drinking straws, a foam wing and two kitchen scales. It can be built in an afternoon and run safely, and it produces lift and drag coefficients with honest uncertainties.

It also produces numbers that disagree with the textbook, for a reason worth understanding. The rig's Reynolds number is about 2.6 × 10⁴; Ladson's NACA 0012 data (`@ntrs-nasa-tm-4074-ladson-naca0012`) start at 2 × 10⁶. At the bench's Reynolds number the boundary layer stays laminar and separates before it can become turbulent; the section effectively loses camber and suction, the lift slope falls, stall comes early and softly, and the drag is several times the high-Re value. That is the regime of small drones, model aircraft and insects (`P10.14`), and it is why a small UAV's designer does not scale NACA data down. The experiment measures the size of that effect on your table and asks you to show it is the Reynolds number, not the rig.

Success criteria, written first. Before measuring: a dated prediction table with uncertainties, including the predicted scale readings, and a statement of which points will be resolution-limited. After: CL and CD at every angle with uncertainty budgets; the lift slope, CL_max and the flow angle with uncertainties; E_n against the prediction; a comparison with the panel code and Ladson's data with each difference explained; and one follow-up test of a rig effect. In code: every reduction step inverts its forward model and the SYNTHETIC dataset reduces back to its stated truth.

## 2. Understand and derive

### 2.1 How small the forces are

At V = 4 m/s in room air (ρ = 1.204 kg/m³), the dynamic pressure is q = ½ρV² = **9.63 Pa**: the weight of a 0.98 mm column of water. On the lesson's model (S = cs = 0.100 × 0.200 = 0.020 m²) a lift coefficient of 0.5 is a force of 0.096 N, **9.8 g** on a kitchen scale, and a drag coefficient of 0.05 is **0.98 g**. A 1 g kitchen scale resolves the lift to about 10 % and cannot resolve the drag; that is why the recommended tier uses a 0.1 g pocket scale for lift and a 0.01 g jewellery scale behind a bell-crank of ratio 2 for drag. The Reynolds number is Re = ρVc/μ = **2.66 × 10⁴**, 75 times below Ladson's lowest.

### 2.2 Airspeed: the pitot-static tube and the inclined manometer (recap of `P08.19`)

A pitot tube facing the flow brings the air to rest at its mouth, where the pressure is the stagnation pressure p₀ + q (Bernoulli, incompressible); side holes on a static tube read the free-stream pressure p₀. Connect the two to the legs of a U-tube of water. The water stands higher on the static side by a vertical height h with (ρ_w − ρ)g h = q. A height of 1 mm is hard to read, so one leg is laid along a board inclined at a small angle θ: the meniscus then moves Δl = h/sin θ along the tube, and

**q = (ρ_w − ρ) g Δl sin θ.**   (1)

At 10° a 0.98 mm column becomes **5.67 mm** along the tube. The incline's own uncertainty enters as u_q/q = (cot θ) u_θ: at 10° ± 0.2°, 2.0 %. A single reading pair read to ±0.5 mm (a uniform band) gives u_q/q = 7.6 % and V = 3.97 ± 0.15 m/s; averaging six pairs brings the manometer term to 3.7 %. A vane anemometer specified to ±(3 % + 0.1 m/s) gives 6.4 % on q. Then V = √(2q/ρ), with half q's relative uncertainty.

### 2.3 The balances

**Lift.** The model sits on a vertical sting standing on the lift scale. With the fan on, the air pushes the model up by L, and the scale reading **falls** by L/g. Read the scale fan off and fan on and subtract; the sting alone, without the model, gives a tare (its own lift, usually tiny) that is subtracted too:

**L = −g[(m_on − m_off) − (m_on − m_off)_tare].**   (2)

**Drag.** The sting is mounted on the vertical arm of a bell-crank pivoting on a nail; drag D at height l_D above the pivot pushes the horizontal arm, of length l_s, down on the drag scale with a force R D, R = l_D/l_s (moments about the pivot). The sting's own drag is not small (the lesson's synthetic sting has about 0.005 N, half the model's drag at zero lift), so the sting-only tare matters:

**D = g[(m_on − m_off) − (m_on − m_off)_tare]/R.**   (3)

R is **calibrated**, not measured with a ruler: a thread tied to the model at the quarter chord runs horizontally over a pencil pulley to known hanging masses, and the slope of the drag-scale reading against the applied load is R (with friction showing as the difference between loading and unloading). A calibrated spring read against a ruler is the alternative, D = kΔx.

### 2.4 Coefficients, uncertainty and blockage

C = F/(qS). Because F can be near zero (lift at α = 0), propagate absolute uncertainties (`P00.14`):

**u_C² = (u_F/(qS))² + (C u_q/q)² + (C u_S/S)².**   (4)

u_F combines the scale resolution (a reading difference of two readings, each uniform over the resolution: √2 r/√12), the repeatability of the two sweeps, the tare's uncertainty and, for drag, R's. u_q combines the manometer, the fan's run-to-run unsteadiness and **blockage**: the model obstructs part of the jet. A simple closed-section estimate of the effect on the velocity is ε = A_frontal/(4A_jet) (for the lesson's model, 0.0087 at 0° and 0.018 at 10° in a 300 × 300 mm jet). In an open jet the correction is smaller and can change sign, so this lesson does not apply it; it carries ±2ε on q as an uncertainty.

### 2.5 The predictions

**Thin airfoil** (`P10.04`): cl = 2πα for a symmetric section, 0.110 per degree. **Panel** (`P10.05`): the 12 % thickness raises it to 0.121 per degree, inviscid. **Finite span** (`P10.09`): the end plates do not make the flow two-dimensional; with an effective aspect ratio AR_eff, a = a₀/(1 + a₀/(πAR_eff)). The geometric aspect ratio is 2 and the plates raise the effective value by an uncertain amount; the lesson brackets it as 4–6, which gives **0.073–0.082 per degree**. **Drag**: laminar (Blasius) friction on both surfaces, 2 × 1.328/√Re, times a thickness form factor 1 + 2t/c + 60(t/c)⁴, plus a lost fraction (0 to 1) of leading-edge suction, |cl| tan α, plus induced drag cl²/(πeAR_eff): at Re = 2.66 × 10⁴, **CD ≥ 0.020** at zero lift.

**With a boundary-layer estimate** (`P10.06`, recapping `P08.11` and `P08.12`). Thwaites' method gives θ² = (0.45ν/Ue⁶)∫Ue⁵ds and λ = (θ²/ν)dUe/ds; laminar separation is at λ ≈ −0.09. Because θ² ∝ ν, λ and hence the **separation point do not depend on Re**. Michel's criterion for transition, Re_θ > 1.174(1 + 22 400/Re_x)Re_x^0.46, **does**. On the panel method's upper surface of the NACA 0012:

| α | laminar separation x/c | Re = 2.6 × 10⁴: Michel x_tr/c, Re_θ at separation | Re = 2.0 × 10⁶: Michel x_tr/c |
|---|---|---|---|
| 0° | 0.626 | 1.000 (none), 105: separates laminar | 0.393: transition first |
| 2° | 0.451 | 0.997, 96: separates laminar | 0.230: transition first |
| 4° | 0.282 | 0.984, 85: separates laminar | 0.126: transition first |
| 6° | 0.059 | 0.941, 47: separates laminar | 0.069: separates laminar (near the leading edge) |
| 8° | 0.022 | 0.817, 34: separates laminar | 0.050: separates laminar |

At 2 × 10⁶ the boundary layer transitions before it can separate up to about 4°, and beyond that a short leading-edge bubble re-attaches turbulent. At the bench's Reynolds number it separates laminar at every angle, with a momentum-thickness Reynolds number of order 100: the separated layer may not re-attach at all. That is the prediction of the gap, made before measuring.

## 3. Build the computational model

Open `starter/main.py`. The data loaders, the NACA generator, the compact panel method (from `P10.05`), the SYNTHETIC generator and its writer and reader, and `main` are given.

### Step 1: air and airspeed

`air_density(T, p, RH)` is the moist-air ideal gas (dry-air plus vapour partial pressures, Magnus–Tetens saturation pressure); `air_viscosity` is Sutherland's law; `manometer_q` is (1) with both uncertainty terms; `airspeed` and `anemometer_q` propagate to V and q; `reynolds_number` is ρVc/μ.

### Step 2: balances and coefficients

`fit_line` is least squares with uncertainties from the residual scatter; `pivot_calibration` and `spring_calibration` use it; `lift_from_scale` is (2), `drag_from_pivot` (3), `drag_from_spring` D = kΔx; `coefficient` is (4); `blockage_fraction` is §2.4; `en_number` is |x₁ − x₂|/(2√(u₁² + u₂²)).

### Step 3: predictions and the reduction

`thin_airfoil_cl`, `finite_span_slope`, `low_re_cd` and `laminar_upper_surface` are §2.5; `predict_bench` assembles the prediction table with the scale readings it implies. `reduce_bench(ds)` takes the five tables of a dataset (conditions, airspeed, calibration, sweep, tare) and returns CL and CD with uncertainties at every angle, the slope over −4° to 4° (with the fit's uncertainty and with q's added, since an error in q scales every CL together and is invisible in the fit's scatter), the zero-lift angle (the jet's flow angularity) and CL_max. `compare_with_ladson` interpolates Ladson's Table I at the bench angles.

### Verify

`pytest curriculum/phase-10/lesson-23-home-lab-airfoil-bench` runs 16 property tests. Then `python main.py` prints the prediction:

| α | cl thin | cl panel | CL band (AR_eff 4–6) | CD band | lift change [g] | drag-scale change [g] |
|---|---|---|---|---|---|---|
| 0° | 0 | 0 | 0 | 0.0204 | 0 | 0.80 |
| 2° | 0.219 | 0.241 | 0.146–0.164 | 0.0217–0.0285 | 2.87–3.23 | 0.85–1.12 |
| 4° | 0.439 | 0.483 | 0.292–0.329 | 0.0254–0.0530 | 5.74–6.46 | 1.00–2.08 |
| 6° | 0.658 | 0.723 | 0.439–0.493 | 0.0317–0.0938 | 8.61–9.69 | 1.25–3.68 |
| 8° | 0.877 | 0.963 | 0.585–0.658 | 0.0406–0.1511 | 11.49–12.92 | 1.59–5.94 |
| 10° | 1.097 | 1.201 | 0.731–0.822 | 0.0519–0.2252 | 14.36–16.15 | 2.04–8.85 |

(The drag band widens with α because the leading-edge-suction loss is unknown between 0 and 1.) On a 1 g lift scale, the 2° points are resolution-limited; on a 1 g drag scale, all of them are.

### Compare

The panel code converges at second order (`P10.05`); its cl at 4° agrees with the thin-airfoil value times 1 + 0.77t/c within 2 %, a test. The SYNTHETIC generator writes exactly the shipped files (a test compares them byte for byte), and reducing them recovers the generator's arm ratio, flow angle and lift slope within the computed uncertainties (a test).

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the data formats are in `projects/p10-airfoil-bench/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). Read the safety box first, every session. If you have no kit, run the whole reduction on the SYNTHETIC files in `projects/p10-airfoil-bench/data/synthetic/` (or on `make_synthetic_bench_data` with your own seed) and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p10-airfoil-bench/safety.md`).** Hazards: fingers, hair, scarves, string, paper or straws drawn into a fan; a fan that tips or walks; parts blown out of the jet; a cable tripped over or pulled into the fan; an overheated motor; knife cuts and hot-wire fumes; water near a mains appliance. **Hard limits: the fan is guarded on both sides (front and back grilles intact, finger-proof), rated at most 50 W, and stands on a stable surface where it cannot tip or walk. No propellers, drones or unguarded blades, ever. No cutting tools near the running fan. Hot-wire foam cutting only with adult supervision and ventilation; otherwise a craft knife, with the fan unplugged. Unplug the fan before adjusting the model** and before reaching anywhere near the fan. Energy limits: 50 W, a 3–5 m/s jet, mains only inside the unmodified fan. Controls: the cable routed so it cannot be tripped over or pulled in; hair tied back, loose paper, scarves and string away from the intake; the straws taped into their frame and the model and balance secured so nothing can be blown into the fan; safety glasses while the fan runs; the manometer below and away from the fan on a tray. Supervision: an adult may work alone; a child only with an adult present whenever the fan is plugged in. Location: a sturdy table, space around the fan, no pets or small children loose. **Stop criteria: the fan rattles, wobbles or walks; a grille is damaged or anything touches one; the motor smells hot.** Not permitted: any propeller, drone or unguarded blade, fans above 50 W or with damaged grilles, opening or modifying the fan, reaching near a plugged-in fan.

The setup, seen from the side:

```
   guarded fan      straw honeycomb   settling    model (foam NACA 0012, end plates)
   (≤ 50 W,        ┌───────────────┬─────────┐       ╭───╮  ← pitot-static tube beside it
   grilles both     │|||||||||||||||│         │  →   ╰─┬─╯     to the inclined manometer
   sides)           │|||||||||||||||│         │        │ sting (pin at c/4, protractor)
     ║█║  →→→→→→    │|||||||||||||||│         │        │
     ║█║            └───────────────┴─────────┘   pivot ●───┐  bell-crank: l_D above,
   ══╩═╩══ table ═══════════════════════════════════   ▲    └── l_s → [drag scale 0.01 g]
                                                  [lift scale 0.1 g under the base]
```

**Prediction first (in `data/prediction.md`, dated).** Measure the model (chord, span, thickness) and the jet; record the room air; predict V (from one no-model reading or the fan's datasheet) and q, the manometer reading, Re; and for each planned angle the CL and CD bands and the scale changes (`predict_bench`), marking the resolution-limited points. Write one sentence on what you expect relative to inviscid theory and Ladson's data, and why. Date the file and do not edit it afterwards.

**Procedure in outline** (full version in `procedure.md`). Calibrate the lift scale with known masses, the bell-crank with horizontal loads over the pulley (loading and unloading), and the manometer incline by rise over run. Tare: sting only, fan off and on, three times; model at 0°, three times. Sweep α from −4° to 14° in 2° steps, **unplugging the fan before every change of angle**, reading both scales fan off and fan on (20 s after switching on), photographing each display; then sweep back down. Read the manometer at the start, middle and end of each sweep.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_sweep.csv`, one row per angle per sweep, never edited; the other four tables are described in `data/README.md`):

```
run,alpha_set_deg,lift_off_g,lift_on_g,drag_off_g,drag_on_g
1,-4,352.7,356.3,18.41,20.9
1,0,352.6,351.8,18.4,20.81
```

## 5. Validate and improve

Run the reduction on your data (or the SYNTHETIC files) and put measurement, computation and published data side by side. The SYNTHETIC dataset (seed 230, `projects/p10-airfoil-bench/data/synthetic/`) was generated with stated assumptions (its README lists them): a low-Re lift slope of 0.060 per degree with a soft stall, CD from 0.035 up to about 0.1, a flow angularity of +0.6°, a drifting fan, a sting with drag, an arm ratio of 2.04 against the nominal 2.0, and a manometer incline recorded as 10.0° though the generator used 10.15°.

**SYNTHETIC results (not a measurement).** Air 1.1868 ± 0.0041 kg/m³; q = 9.48 ± 0.35 Pa from six manometer pairs; V = 4.00 ± 0.07 m/s; q drifted −2.9 % between the readings before and after the sweep; Re = 2.61 × 10⁴. Arm ratio R = 2.038 ± 0.004 (nominal 2.0); sting tares −0.03 g of lift and 1.05 ± 0.03 g on the drag scale.

| α_set | CL ± u | CD ± u | ε (blockage) |
|---|---|---|---|
| −4° | −0.185 ± 0.010 | 0.0359 ± 0.0019 | 0.0125 |
| −2° | −0.087 ± 0.006 | 0.0356 ± 0.0031 | 0.0106 |
| 0° | 0.037 ± 0.005 | 0.0335 ± 0.0019 | 0.0087 |
| 2° | 0.154 ± 0.009 | 0.0367 ± 0.0018 | 0.0106 |
| 4° | 0.267 ± 0.014 | 0.0424 ± 0.0023 | 0.0125 |
| 6° | 0.386 ± 0.022 | 0.0508 ± 0.0039 | 0.0145 |
| 8° | 0.459 ± 0.024 | 0.0564 ± 0.0031 | 0.0164 |
| 10° | 0.529 ± 0.031 | 0.0708 ± 0.0041 | 0.0183 |
| 12° | 0.531 ± 0.037 | 0.0823 ± 0.0080 | 0.0202 |
| 14° | 0.508 ± 0.033 | 0.0951 ± 0.0059 | 0.0221 |

Slope over −4° to 4°: **0.0573 ± 0.0011 per degree from the fit's scatter, ± 0.0024 with q's uncertainty**; zero lift at α_set = **−0.65 ± 0.05°** (the jet's flow angularity: the generator's +0.6°, found); CL_max **0.531 ± 0.037 at 12°**.

**Three comparisons.**

| quantity | inviscid prediction | SYNTHETIC bench | Ladson, 2 × 10⁶ | verdict |
|---|---|---|---|---|
| lift slope [/deg] | 0.073–0.082 (AR_eff 4–6) | 0.057 ± 0.002 | 0.104 (2-D) | E_n = 2.86 against the band's middle: disagrees; bench/Ladson = 0.55 |
| CL at 4.65° (flow-angle corrected) | 0.34–0.38 | 0.267 | 0.484 (interpolated) | ratio to Ladson 0.55 |
| CL_max | — | 0.53 near 12° | 1.41 at 15.3° | early, soft stall |
| CD at zero lift | ≥ 0.020 (laminar friction) | 0.034 | 0.0062 | 5.4 × Ladson |

Each difference has a named cause. The slope is below the inviscid finite-span band and far below Ladson's: the laminar separation predicted in §2.5 takes away camber and suction, exactly the low-Re signature. Part of the gap to Ladson is also three-dimensional (the end-plated model is not a 2-D section, which the band allows for). The drag is above the laminar-friction estimate by 0.013 at zero lift: the estimate ignores the pressure drag of laminar separation, so it is a lower bound by construction. The generator's planted errors also show: the slope recovered from this seed is 4.5 % below its stated truth of 0.060, a gap of 1.1u of the total uncertainty (2.5u of the fit's alone). The planted incline error does not explain it: recording the incline as 10.0° when it was 10.15° biases q low by 1.5 %, which raises every CL and so the slope; correcting it moves the slope to 0.0565 per degree, further from the truth. It is nearly cancelled by a second planted effect, the fan's drift: the sweep ran at a mean speed about 0.8 % below that of the manometer readings, so its true q was about 1.6 % lower and every CL comes out low by the same amount. What remains is the soft stall, which bends the truth to 0.0597 per degree over the fitted angles, and this seed's scatter. A fit uncertainty that ignores q would have called a 2.5u gap a disagreement; with q's uncertainty it is consistent at k = 2.

| pattern in your data | likely cause | test |
|---|---|---|
| lift at α_set = 0 not zero | flow angularity or a misaligned protractor | turn the model over and repeat ±α; the curves shift oppositely |
| slope rises with fan speed | Reynolds-number effect | three fan speeds, same rig |
| slope the same at all speeds but low | end-plate leakage, jet non-uniformity, sting compliance | move the model 30 mm off-centre; check α under load in a photo |
| drag reading drifts while the fan runs | the fan warming up; the pivot sticking | read later; tap the bench gently and re-read; repeat calibration |
| large scatter between sweeps | fan unsteadiness, readings taken too early | wait 20 s; take the middle of the fluctuation; more repeats |
| CD negative at small α | tare larger than the model's drag change; R wrong | repeat the sting tare; recalibrate R |

**Improve and re-run.** For CL at 4° the largest term is q (4.6 %: the manometer 3.7 %, blockage 2.5 %); a digital pressure sensor or more manometer readings with a better-measured incline helps most. For CD at zero lift, the force (4.0 %) and q (4.2 %) terms are similar. Re-run one sweep after the improvement and report both.

## 6. Defend your work

Submit: `data/prediction.md` (dated), the calibration records, `data/raw/*.csv`, photographs, your analysis script, the CL and CD plots with k = 2 bars, the prediction bands, the panel curve and Ladson's data, the table of three comparisons with E_n, the ranked explanations with your follow-up test, `airfoil_bench.png`, and a report following `projects/p10-airfoil-bench/report-template.md`, with the completed safety box. If you used the SYNTHETIC dataset, every table and figure says so.

### Adversarial self-checks

- **What if the airspeed doubles?** q quadruples (forces 4 ×, resolution terms 4 × smaller relative), Re doubles to about 5 × 10⁴: the slope should rise toward the high-Re value if the gap is a Reynolds-number effect. The fan's 50 W limit is fixed, so this is a thought experiment unless your fan has a higher speed setting within its rating.
- **What if the angle goes to zero?** A symmetric section must give zero lift; any lift is flow angularity or misalignment (the SYNTHETIC set shows 0.65°), which is why the slope, not the lift at one angle, is the robust result.
- **Which assumption fails first?** A uniform, steady jet: a desk fan's flow has a hub wake, swirl and unsteadiness that the straws reduce but do not remove.
- **What measurement could falsify the Reynolds-number explanation?** A slope that does not change across fan speeds, or a trip strip that does not change the slope or the drag; either would point to a rig effect.
- **Does the result survive a change of method?** Lift from the scale and lift inferred from the manometer and a second balance (or a load cell) should agree; drag from the bell-crank and from a calibrated spring should agree within their uncertainties.
- **Which decision changes if the uncertainty doubles?** At ± 0.005 per degree the slope still differs from Ladson's by an order of magnitude more than its uncertainty, so the qualitative conclusion stands; a claim that the bench slope agrees with the inviscid finite-span band (E_n would fall to about 1.9) would become untestable.

## Check yourself

1. What is the dynamic pressure at 4 m/s in room air, and how high a water column is it? Why is the manometer inclined?
2. What is the bench rig's Reynolds number, and how many times lower is it than Ladson's lowest?
3. Which way does the lift scale move when the model lifts, and what do the two tare runs remove?
4. Why does a NACA 0012 at Re ≈ 3 × 10⁴ have a lower lift slope, an earlier and softer stall and a much higher drag than at 2 × 10⁶?
5. Why must the prediction be written and dated before measuring?
6. State the hard safety limits for the fan and for cutting the foam.
7. Why does the slope's uncertainty from the line fit understate the real uncertainty?
