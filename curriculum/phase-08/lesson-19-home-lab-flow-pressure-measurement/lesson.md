# P08.19  Home Lab: Flow & Pressure Measurement 🔬

> Bernoulli says a hole 4 mm across drains a 90 mm bottle from 200 mm to 30 mm in 103 seconds, if its discharge coefficient is the textbook 0.61. Write that down with its uncertainty, pull the tape, and time it. Then let a column of water tell you how fast your hair dryer blows.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the draining law of a bottle from Bernoulli's equation and mass conservation, with a discharge coefficient, show that √h falls linearly in time, and predict a draining time with an uncertainty budget before measuring. 2. Measure the level against time from a phone video, fit C_d with a weighted straight line, and compare it with the prediction and with the typical sharp-edged value by a normalised error E_n, identifying the dominant uncertainty (the hole diameter, which enters squared). 3. Build a water U-tube manometer with a drinking-straw pitot probe, derive Δp = (ρ_w − ρ_a) g Δh and V = √(2Δp/ρ_a), and calculate, before measuring, the lowest airspeed it can read to 10 % upright and tilted. 4. Measure the airspeed of a hair dryer on cool and of a small guarded fan, calibrate the tilted manometer against the upright one, and judge each result against its resolution. 5. Use Reynolds numbers to decide whether viscosity matters in each experiment. |
| **Prerequisites** | `P08.04` (Bernoulli's equation, the pitot tube), `P08.03` (mass conservation for a control volume), `P08.01` (hydrostatics and the manometer, air density), `P08.05` (Reynolds number), `P00.13` (home-experiment safety), `P00.14` (measurement, calibration and uncertainty), `P00.15` (the first experiment, normalised error), `P01.17` (least squares), `P07.01` Systems, Properties & State (the ideal-gas law p = ρRT for air: density from the room's pressure and temperature, with R = 287.05 J/(kg·K); restated here because Phase 07 may not be published yet). Skill check: you can propagate the uncertainty of y = c x²/z by the relative-uncertainty rule. |
| **Engineering question** | Can a kitchen-table rig measure a discharge coefficient to about 5 % and an airspeed of a few metres per second to 10 %? The predictions, before anything flows: a 90.0 ± 0.5 mm bottle with a 4.0 ± 0.1 mm hole drains from 200 to 30 mm above the hole in 102.7 ± 7.3 s (7.1 %, half from the hole diameter, half from C_d = 0.61 ± 0.03); the jet's Reynolds number is 3000–8000, so viscosity barely matters. A water U-tube read to 0.3 mm per meniscus resolves airspeeds to 10 % only above 7.0 m/s upright (a hair dryer on cool), 4.9 m/s tilted at 30° and 2.9 m/s tilted at 10° (a small fan). |
| **Assumptions and limits** | Quasi-steady, inviscid flow from the free surface to the vena contracta, with all real-fluid effects lumped into C_d; a straight-walled bottle; level read above the hole centre; (a/A)² ≪ 1 (here 4 × 10⁻⁶). Below a few millimetres of head the jet is held back by surface tension and dribbles (2σ/d to 4σ/d is 4–7.5 mm of water), so the fit stops at 30 mm. Pitot: the probe faces the stream within about ±15°, its Reynolds number is above about 100, the room air is at the jet's static pressure (true in a free jet), incompressible flow (M < 0.05). The manometer: the same bore in both legs, water at room temperature, no bubbles. |
| **Analytical work** | Bernoulli plus continuity for the draining bottle; the √h law and the draining time; C_d = C_c C_v; the exact correction for (a/A)²; the uncertainty budget of the time and of C_d; the manometer relation with the air-column correction; the tilted tube and its gain 1/sin θ; the uncertainty of V and the resolution limit; Reynolds numbers of the jet, the probe and the source; the capillary head that stops the jet. |
| **Computational deliverable** | `starter/main.py`: `orifice_area`, `drain_rate_constant`, `torricelli_level`, `drain_time`, `predict_drain_time`, `fit_line`, `fit_discharge_coefficient`, `jet_reynolds`, `air_density`, `manometer_dp`, `vertical_head`, `airspeed`, `airspeed_with_uncertainty`, `min_speed_for_uncertainty`, `reynolds`, `fan_mean_velocity`, `en_number` (the constants, the SYNTHETIC data generator and writer, `read_csv` and `main` are given). Tests (14, under 2 s) prove: the Torricelli level satisfies the draining ODE and empties at the predicted time; the draining time scales as 1/C_d, (D/d)² and 1/√g; first-order propagation agrees with the power-law rule and a seeded Monte Carlo; the weighted fit recovers an exact line and has mean χ² = dof on seeded noise; C_d is recovered from noise-free levels and from the SYNTHETIC data within its uncertainty; the air density reproduces ISA sea level; the manometer, tilt and airspeed formulas invert each other; the resolution limit gives exactly the target uncertainty; the shipped SYNTHETIC files match the seeded generator byte for byte. |
| **Practical activity** | **physical experiment**: drain a bottle through a small hole and fit its discharge coefficient; read a water U-tube with a straw pitot probe in the stream of a hair dryer on cool and of a small guarded fan, upright and tilted, with a dated prediction before each. See `projects/p08-flow-pressure/`. A learner without a kit completes the lesson with the SYNTHETIC dataset in `projects/p08-flow-pressure/data/synthetic/`, labelled as synthetic. |
| **Equipment** | Minimum: a 1–1.5 L PET bottle with a 4 mm hole, tape, two rulers, a phone and stand, a basin, a measuring jug, 1.5 m of 4–6 mm clear tubing on a board, a straw, a hair dryer used on cool or a guarded fan of at most 30 W, a thermometer and the phone's barometer. Recommended: a caliper and drill bits as hole gauges, a ramp to tilt the board. Enhanced: a vane anemometer. See `projects/p08-flow-pressure/bom.yaml` (reuses `p00-calibration`, `p00-first-experiment`, `p04-pendulum` and `p04-drag-drop` items). |
| **Safety** | See `projects/p08-flow-pressure/safety.md`. Summary box at the start of part 4: water at room temperature, no more than 2 L drawn (repeat runs reuse it); fans guarded and at most 30 W; hair dryer on cool only and away from water; no compressed gas, no pressurised containers, no mains appliance near spilled water. |
| **Evidence** | `data/prediction.md` (dated before anything flows), `data/calibration.csv`, `data/raw/*.csv` with the video and photo names, the √h–t plot with the prediction, the C_d table with its budget and E_n, the airspeed table with uncertainties and Reynolds numbers, `flow_lab.png`, the report. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; all four competencies. Threshold: a dated prediction with uncertainty before each experiment; all public tests pass; C_d reported with its budget and E_n; airspeeds reported with uncertainties, and any reading below the predicted resolution reported as a bound, not a number; the safety box completed. |
| **References** | `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual` 1.6 "Experiment #6: Orifice and Free Jet Flow": 7. Theory (7.1 Determination of the Coefficient of Velocity; 7.2 Determination of the Coefficient of Discharge) (required-free: the laboratory version of experiment 1; the home procedure here is different). `@openstax-university-physics-vol1` 14.2 Measuring Pressure (Manometers; Example "Fluid Heights in an Open U-Tube"); 14.6 Bernoulli's Equation; 14.7 Viscosity and Turbulence (Measuring Turbulence) (required-free). `@leishman-intro-aerospace-flight-vehicles` chapter "Fluid Statics & the Hydrostatic Equation": U-Tube Manometer; chapter "Applications of the Conservation Laws": Pitot Tube & Pitot-Static Tube (required-free; https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/). `@nasa-glenn-beginners-guide-aeronautics` "Pitot – Static Tube – Speedometer"; "Reynolds Number" (supplementary-free). `@mit-ocw-1-060-engineering-mechanics-ii` L4 "Manometry, Stability of Dams and Gates" (supplementary-free). |
| **Extensions** | ★ Drill a second bottle with a 2 mm hole and measure C_d at jet Reynolds numbers of 1000–4000: does it change? ★ Measure the hair dryer's velocity profile across its jet at 50, 100 and 200 mm from the nozzle and integrate it to a volume flow; compare with the dryer's rating plate. ★ Flare one hole with a countersink (a rounded entry) and compare its C_d with the sharp hole's: the vena contracta disappears and C_d rises towards 0.9. |
| **Time** | unmeasured |

## 1. The engineering problem

Two measurements sit under almost every flow test in aerospace. The first is the **discharge coefficient** of a flow element: the injector orifices of a rocket engine, the restrictor in a fuel line, the metering hole of a carburettor, the nozzle of a flow bench. Ideal theory gives a flow rate; the real element passes less, and the ratio, C_d, is measured, because nothing else gives it to a few per cent. The second is **airspeed from a pressure difference**: the pitot tube on every aeroplane (`P08.04`) and the pitot-static rakes of every wind tunnel turn a pressure into a speed through Bernoulli's equation, and the pressure is read by a transducer that was calibrated, ultimately, against a column of liquid.

This lab does both with household things. A plastic bottle with a small hole near its base is an orifice under a falling head: Torricelli's law, corrected by C_d, predicts how its level falls with time, and a phone video measures it. A clear tube bent into a U and half-filled with water is a manometer, the primary standard of low pressure; a drinking straw facing the stream of a hair dryer on its cool setting, or of a small guarded fan, is a pitot probe. The physics is a single equation each time. The engineering is in the prediction you write before you measure, in knowing which uncertainty dominates, and in recognising when your instrument cannot see what you are trying to measure.

Success criteria, written first. Before measuring: a dated prediction with uncertainties for the draining time and for the smallest airspeed each manometer arrangement can read. After: C_d with its uncertainty budget and E_n against the typical sharp-edged value and against the prediction; χ² of the fit consistent with the reading uncertainties; airspeeds with uncertainties, a tilted reading calibrated against the upright one, and an honest "below resolution" where the prediction says so.

## 2. Understand and derive

### 2.1 The draining bottle

Let the bottle have a straight section of area A = πD²/4 and a round hole of area a = πd²/4, with the water level a height h above the hole centre. Take a streamline from the free surface (pressure atmospheric, speed V_s) to the narrowest section of the jet just outside the hole, the **vena contracta**, where the streamlines are parallel and the pressure is again atmospheric. For steady (quasi-steady: h changes slowly), inviscid, incompressible flow, Bernoulli's equation (`P08.04`; `@openstax-university-physics-vol1` 14.6) gives

V_j² − V_s² = 2gh.

Mass conservation (`P08.03`): A V_s = a_c V_j, where a_c is the contracted jet area. If a_c ≪ A, V_s² is negligible and V_j = √(2gh), Torricelli's result. Two corrections are lumped into the **discharge coefficient**: the jet contracts to a_c = C_c a (C_c ≈ 0.61–0.64 for a sharp-edged hole, because the flow arriving from the side cannot turn the corner sharply), and friction slows it slightly, V_j = C_v √(2gh) (C_v ≈ 0.97–0.99). The volume flow is

Q = C_d a √(2gh),  C_d = C_c C_v.

The textbook range for a sharp-edged round orifice is about 0.60–0.62 (typical values; `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual` 1.6, 7.2 shows the laboratory method). A hole drilled in a thin plastic wall is close to sharp-edged but rarely perfect; the lesson takes C_d = 0.61 ± 0.03.

The level falls as the bottle loses volume: A dh/dt = −Q = −C_d a √(2gh). Divide by √h:

d(√h)/dt = −C_d (a/A) √(g/2) ≡ −k,  so  **√h = √h₀ − k t,  k = C_d (d/D)² √(g/2).**   (1)

**√h is a straight line in time**, and its slope gives C_d. The time to fall from h₁ to h₂ is

**T = (√h₁ − √h₂)/k.**   (2)

The neglected surface velocity: keeping V_s gives V_j² = 2gh/(1 − (C_d a/A)²), a relative change of (C_d a/A)²/2 = 7 × 10⁻⁷ for the lesson's bottle. The quasi-steady assumption is as good: the level changes by less than a millimetre in the time a fluid particle takes to cross the bottle.

**The lesson's bottle.** D = 90.0 ± 0.5 mm, d = 4.0 ± 0.1 mm, hole centre 30 mm above the base, filled to 200 mm above the hole: 1.46 L of water. With C_d = 0.61, (2) gives 102.7 s from 200 to 30 mm and 167.6 s to empty.

### 2.2 The uncertainty of the prediction

T ∝ C_d⁻¹ d⁻² D², so by the relative-uncertainty rule (`P00.14`)

(u_T/T)² = (u_Cd/C_d)² + (2u_d/d)² + (2u_D/D)² + (small level terms).

For the lesson's bottle: C_d contributes 4.9 % (5.05 s), the hole 5.0 % (5.13 s), the bottle 1.1 % (1.14 s), the level readings 0.5 mm each under 0.6 s. Total **102.7 ± 7.3 s (7.1 %)**. `predict_drain_time` does this with numerical sensitivities, and the tests check it against a seeded Monte Carlo. Read the budget before measuring: the hole diameter is as important as the unknown C_d. Gauging the hole with drill shanks to ± 0.05 mm halves its term.

### 2.3 Fitting C_d

From the video, record the time t_i at which the level crosses each ruler mark h_i (200, 190, …, 30 mm). Fit √h_i = c − k t_i by weighted least squares (`P01.17`). Each √h_i is uncertain from the level reading, u_h/(2√h), and from the timing, k u_t; adding them in quadrature is the effective-variance method. Then

**C_d = k (D/d)² √(2/g),  (u_Cd/C_d)² = (u_k/k)² + (2u_d/d)² + (2u_D/D)².**   (3)

The fit's slope is typically known to a fraction of a per cent; the geometry dominates again. **Repeated runs agree with each other far better than this**, because they share the same hole: their scatter measures the timing and the reading, not the hole. The χ² of the fit should be close to its degrees of freedom (marks − 2) if the reading uncertainties are honest.

### 2.4 Does viscosity matter? The jet's Reynolds number

Re_jet = ρ_w √(2gh) d/μ_w is **7900 at 200 mm and 3100 at 30 mm** for the 4 mm hole (`@openstax-university-physics-vol1` 14.7). Viscosity acts only in a thin boundary layer at the hole's edge; C_d of a sharp orifice changes little over this range. With a 1 mm hole Re falls to 760–2000 and C_d begins to depend on it (Extension 1). Surface tension matters only at the end: the pressure inside a jet of diameter d exceeds the outside by about 2σ/d, and forming a drop at the hole needs 4σ/d; with σ = 0.073 N/m that is 3.7–7.5 mm of water, so below a few millimetres of head the jet dribbles and Torricelli no longer applies. That is why the marks stop at 30 mm.

### 2.5 The manometer and the pitot probe

A straw whose open end faces the stream brings the flow at its tip to rest: the pressure inside it is the **total pressure** p₀ = p + ½ρ_aV² (`P08.04`; `@leishman-intro-aerospace-flight-vehicles` "Pitot Tube & Pitot-Static Tube"). The other leg of the manometer is open to the still room. In a **free jet** (the stream of a dryer or a fan, open to the room) the static pressure equals the surrounding pressure, because a jet cannot sustain a pressure difference across its boundary; so the open leg reads the jet's static pressure, and the manometer reads p₀ − p = ½ρ_aV².

The U-tube (`P08.01`; `@openstax-university-physics-vol1` 14.2, "Fluid Heights in an Open U-Tube"): the pressure difference pushes the water down one leg and up the other until the hydrostatic heads balance. Above the lower meniscus the air columns differ in height by Δh too, so

**Δp = (ρ_w − ρ_a) g Δh,  V = √(2Δp/ρ_a).**   (4)

The air correction is 0.12 %. ρ_a = p/(RT) from the room's pressure and temperature (`P07.01`): 1.196 kg/m³ at 101.0 kPa and 21 °C.

**The resolution problem.** At 3 m/s, Δp = ½ × 1.196 × 9 = 5.4 Pa and Δh = **0.55 mm**. Each reading of a meniscus from a zoomed photo of a 1 mm scale is uncertain by about 0.3 mm, and Δh is the difference of four readings (two menisci, with and without flow), so u_Δh ≈ 0.6 mm. Since V ∝ √Δh, u_V/V = ½ u_Δh/Δh: a 10 % airspeed needs Δh ≥ 3 mm, which is **7.0 m/s**. A hair dryer on cool, at 8–15 m/s near the nozzle, is readable; a desk fan at 3 m/s is not.

**The tilted tube.** Tape the U-tube to a board tilted at θ to the horizontal, with both legs along the board. A vertical head Δh moves the menisci apart by L = Δh/sin θ along the tube: the reading is magnified by 1/sin θ (2 at 30°, 5.8 at 10°) while the reading uncertainty stays 0.6 mm along the tube. The 10 % threshold drops to **4.9 m/s at 30° and 2.9 m/s at 10°**. The price is the angle: (u_V/V)² gains a term ¼(u_θ/tan θ)², 2.5 % at 10° for u_θ = 0.5°. Below about 5° the angle, not the reading, dominates.

**Calibrated pressure.** The upright U-tube is a primary instrument: (4) needs only ρ_w, g and a ruler. The tilted tube's gain depends on an angle you measured; calibrate it by reading the same source (the dryer's high setting, at the same spot) both ways, and check L sin θ against the upright Δh.

### 2.6 Reynolds numbers in experiment 2

The 6 mm straw at 3 m/s has Re = 1190, at 12 m/s 4760. A pitot probe reads ½ρV² to within a per cent as long as its Reynolds number (on its diameter) is above about a hundred; below that, viscous effects make it read high. The straw is safe; a 1 mm hypodermic tube at 3 m/s (Re = 200) is near the edge. The fan (280 mm blades at 3.6 m/s, Re ≈ 67 000) and the dryer's 40 mm nozzle at 12 m/s (Re ≈ 32 000) are turbulent jets: their speed fluctuates, and the manometer, whose water column responds in about a second, averages it.

### 2.7 Where the models fail

In order of likelihood at home: the hole diameter is not what you think (burr, out-of-round, measured on the wrong side); the hole edge is not sharp (a flared or punched edge raises C_d); the bottle tapers; the ruler's zero is not at the hole centre; the jet's start is disturbed by the tape. For the manometer: the straw is not aligned with the stream or is outside the jet's core; the dryer's jet decays with distance from the nozzle, so "the airspeed" depends on where you put the straw; bubbles in the tube; the board's angle. Bernoulli itself comes last.

## 3. Build the computational model

Open `starter/main.py`. The constants, `BOTTLE_NOMINAL` and `ROOM_NOMINAL`, the SYNTHETIC generator `make_synthetic_flow_data`, `write_synthetic_csv`, `read_csv` and `main` are given. The generator calls `drain_rate_constant` and `air_density`, so it runs once you have written them.

### Step 1: the draining law

`orifice_area`, `drain_rate_constant` (k of (1)), `torricelli_level` (the level, clipped at zero once empty), `drain_time` (2), and `jet_reynolds`.

### Step 2: prediction and fitting

`predict_drain_time` takes central-difference sensitivities with a relative step of 10⁻⁶ and returns T, u_T and the budget term by term. `fit_line` is weighted least squares with covariance (AᵀWA)⁻¹, χ² and the degrees of freedom. `fit_discharge_coefficient` is (3), with a first unweighted pass to estimate k for the timing term:

```python
k0 = -np.polyfit(t, np.sqrt(h), 1)[0]
u_s = np.sqrt((u_h / (2 * np.sqrt(h)))**2 + (k0 * u_t)**2)
fit = fit_line(t, np.sqrt(h), u_s)
cd = -fit["slope"] * (D / d)**2 * np.sqrt(2 / g)
```

### Step 3: the manometer

`air_density`, `manometer_dp` (4), `vertical_head` (L sin θ), `airspeed`, `airspeed_with_uncertainty`, and `min_speed_for_uncertainty`, which inverts u_V/V = u_Δh/(2Δh) for the lowest readable speed. `fan_mean_velocity` turns a rated volume flow (a manufacturer's figure) into a mean speed through the fan's swept annulus; `en_number` is the normalised error of `P00.15`.

### Verify

`pytest curriculum/phase-08/lesson-19-home-lab-flow-pressure-measurement` runs 14 property tests in under 2 s, listed in the metadata. One of them regenerates the SYNTHETIC files in a temporary folder and checks them byte for byte against `projects/p08-flow-pressure/data/synthetic/`. Then run `python main.py`; the output starts:

```
EXPERIMENT 1: draining bottle (prediction, before any measurement)
  bottle 90.0 +/- 0.5 mm, hole 4.0 +/- 0.1 mm, C_d = 0.61 +/- 0.03 (typical)
  predicted time from 200 mm to 30 mm: 102.7 +/- 7.3 s (7.1 %); budget (s): cd 5.05, d 5.13, D 1.14, h_start 0.21, h_end 0.54
  to empty from 200 mm: 167.6 s; water used up to 200 mm above a hole 30 mm from the base: 1.46 L
  jet Reynolds number: 7894 at 200 mm, 3057 at 30 mm; area ratio (a/A)^2 = 3.9e-06
```

### Compare

The fit is checked by the tests against noise-free levels from (1) (C_d recovered to 10⁻⁹) and the draining ODE against a finite difference of the level (10⁻⁴). The manometer has no computational check worth the name: it is one formula. Its check is the experiment's own: the upright and tilted readings of the same jet must agree.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the data formats are in `projects/p08-flow-pressure/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). Read the safety box first, every session. If you have no kit, run the whole reduction on the SYNTHETIC files in `projects/p08-flow-pressure/data/synthetic/` (or on `make_synthetic_flow_data` with your own seed) and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p08-flow-pressure/safety.md`).** Hazards: water and mains electricity together (a dryer or fan in a basin, a plug in a puddle); slipping on a wet floor; fingers in a fan; hot air from a dryer on a heated setting; cuts while making the hole; a squeezed or pressurised bottle squirting. **Hard limits: water at room temperature, no more than 2 L drawn in the whole session (the repeat runs reuse the same water, caught in a clean basin); fans guarded (finger-proof grille) and drawing at most 30 W; a hair dryer on its cool setting only, kept away from water; no compressed gas and no pressurised containers (never squeeze, cap, pump or blow into the bottle or the manometer); no mains appliance near spilled water.** Controls: **the water and the mains appliance are never on the same table**; the bottle drains into a clean basin in the sink or on the floor; **mop any spill before switching anything on**, with dry hands; cables on the dry side, a socket above floor level, an RCD/GFCI socket where you have one; the hole made before filling, on a cutting board, with a drill bit turned by hand (never a heated nail); the manometer filled with a syringe or by dipping, never by mouth. Supervision: none for an adult; a child only with an adult handling the appliance and the hole. Stop criteria: **water on or near a cable, plug, socket or appliance**; any spill on the floor; a smell of hot plastic or a dryer blowing warm on cool; a damaged fan grille; the bottle leaking elsewhere or bulging; water reaching the end of the manometer tube. On a stop: switch off at the wall and unplug with dry hands, move the appliance away, then mop.

The two set-ups:

```
  EXPERIMENT 1 (in the sink or basin)          EXPERIMENT 2 (separate dry table)
                                                                       board, upright or tilted
     ruler, zero at hole centre                      straw pitot            ┌──────────────┐
     │ 200 ─ ─ ~~~~~ water surface                  ══════╗  ← jet ← ◄)) │ ║    ║ ruler │
     │ 190    │       │                                   ║   dryer on    │ ║    ║       │
     │  …     │ 90 mm │                                   ╚═══════════════╪═╝    ║       │
     │  30    │       │                                                   │  ╚══╝ water  │
     │   0 ───┤ o ────┼──►  jet (4 mm hole, tape pulled at t = 0)          │  open leg ↑  │
     │        └───────┘     into the basin                                └──────────────┘
                                                                   towel under the board
```

**Prediction first (in `data/prediction.md`, dated).** Experiment 1: measure D (inside) and d with their uncertainties and instruments; compute k, the 200 → 30 mm draining time with `predict_drain_time` and its budget term by term, and the jet's Reynolds numbers. Experiment 2: the room's pressure and temperature and ρ_a; your reading uncertainty; `min_speed_for_uncertainty` upright and at your tilt; the fan's mean speed from its rating plate or box (a manufacturer's figure); the Δh you expect for each source; and a decision, written down, on which source needs the tilted tube. For the lesson's kit the prediction is the table in the metadata. Date the file and do not edit it afterwards.

**Procedure in outline** (full version in `procedure.md`). Experiment 1: seal the hole, fill to about 205 mm with the measuring jug (counting the water), start the video, pull the tape, and let it drain below 30 mm; read the crossing times of the 10 mm marks from the video; three runs. Experiment 2: fill the U-tube with about 20 mL of water, tape the straw to one leg, stand the board upright on a towel on a separate dry table; dryer on cool (switched to cool before plugging in), straw 50 mm from the nozzle on its axis; wait 5 s, photograph, remove the straw from the jet, photograph the zero; twice per setting; then the fan, straw 300 mm in front of the hub; then repeat tilted. Enhanced tier: an anemometer at the same spot.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_drain.csv` and `..._manometer.csv`, never edited):

```
run,mark_mm,t_s,video,note
1,200,2.00,VID_0201,
source,setting,distance_mm,angle_deg,repeat,zero_left_mm,zero_right_mm,left_mm,right_mm,anemometer_mps,photo_zero,photo_flow,note
hair dryer cool,high,50,90,1,99.8,99.9,95.6,104.9,,IMG_0311,IMG_0312,
```

## 5. Validate and improve

Run the reduction on your data (or the SYNTHETIC files) and put the comparisons side by side. The SYNTHETIC dataset (seed 819) was generated with a hidden "truth" that differs from the nominal inputs on purpose: C_d = 0.62, a hole of 4.05 mm, a bottle of 90.4 mm, 0.08 s timing noise and 0.5 mm level noise; on the manometer, true speeds of 12.0 and 8.5 m/s (dryer) and 3.6 and 2.8 m/s (fan), 0.3 mm reading noise and 2 % unsteadiness.

**Experiment 1, SYNTHETIC.**

| quantity | prediction (before measuring) | SYNTHETIC result (seed 819) | E_n (k = 2) |
|---|---|---|---|
| time 200 → 30 mm, run 1 | 102.7 ± 7.3 s | 98.7 ± 0.1 s | 0.27 |
| C_d, runs 1, 2, 3 | 0.61 ± 0.03 (typical) | 0.6313, 0.6324, 0.6312, each ± 0.0324 | 0.25 (mean) |
| run-to-run scatter of C_d | — | sd 0.0007 | — |
| budget of u(C_d) | — | slope 0.26 %, hole 5.0 %, bottle 1.1 % | — |
| χ²/dof | about 1 | 14.9/16, 14.4/16, 15.0/16 | consistent |

The fit is consistent with the prediction, and it would still be wrong to say "C_d = 0.632". The generator's true C_d is 0.62; read through the nominal hole and bottle it looks like 0.62 × (4.05/4.0)² × (90.0/90.4)² = 0.630, and that is what all three runs found. The three runs agree to 0.0007 because they share the same wrong hole diameter: **repeatability is not accuracy**. Only a better measurement of d (drill-shank gauges, a caliper) can reduce the 5 % hole term, and it is the term that the experiment exists to beat.

**Experiment 2, SYNTHETIC** (mean of two repeats, with the uncertainty of one reading):

| source | upright: Δh, V | tilted 10°: Δh, V | anemometer (SYNTHETIC, enhanced) |
|---|---|---|---|
| dryer, high | 9.10 mm, 12.2 ± 0.4 m/s | 8.78 mm, 12.0 ± 0.3 m/s | 11.6 m/s |
| dryer, low | 4.65 mm, 8.7 ± 0.6 m/s | 4.51 mm, 8.6 ± 0.2 m/s | 8.6 m/s |
| fan, high | 0.25 mm: below resolution (u_Δh = 0.6 mm) | 0.73 mm, 3.45 ± 0.27 m/s | 3.7 m/s |
| fan, low | 0.35 mm: below resolution | 0.55 mm, 2.99 ± 0.29 m/s | 2.8 m/s |

The calibration check (dryer, high): 9.10 ± 0.42 mm upright against 8.78 ± 0.44 mm tilted, E_n = 0.26; the upright tube confirms the tilted tube's gain to about 5 %. The fan behaves as predicted: unreadable upright, readable to about 8–10 % at 10°. The fan's rated flow (1100 m³/h in the SYNTHETIC example) gives a mean of 5.4 m/s through its annulus, 1.5 times what the manometer reads 300 mm downstream on the axis: the rating is measured in a test rig at the fan face, and the jet has spread and slowed by the time it reaches the straw. It is a weak independent check; the anemometer is a good one.

| pattern in your data | likely cause | test |
|---|---|---|
| C_d above 0.65, χ² fine | hole larger than measured, or a flared / punched edge | gauge with drill shanks; inspect the edge with a magnifier |
| C_d below 0.58 | hole smaller (burr), level zero above the hole centre | deburr and re-gauge; re-check the ruler zero |
| √h–t line curves at the end | surface tension below a few mm; a tapering bottle base | drop the last marks; measure D near the bottom |
| χ² ≫ dof | timing or level uncertainty underestimated | re-read a run frame by frame; film level with the surface |
| manometer reading drifts | water creeping, a leak at the straw joint, temperature | reseal the joint; read the zero before and after |
| tilted and upright disagree (E_n > 2) | angle wrong, board flexing | re-measure the angle by rise over run |

**Improve and re-run.** For experiment 1, gauge the hole to ± 0.05 mm (drill shanks in 0.1 mm steps): the hole term falls to 2.5 % and u(C_d) to about 3 %. For experiment 2, the tilt and a steadier straw position bring the fan within reach; below 5° the angle's uncertainty dominates, so measure it better before tilting further (Problem 5).

## 6. Defend your work

Submit: `data/prediction.md` (dated before anything flowed), `data/calibration.csv`, `data/raw/*.csv`, your analysis script, the √h–t plot with the prediction line and the airspeed plot with k = 2 error bars, the C_d table with budget and E_n, the airspeed table with uncertainties and Reynolds numbers, `flow_lab.png`, and a report following `projects/p08-flow-pressure/report-template.md` with the completed safety box. If you used the SYNTHETIC dataset, every table and figure says so in its caption.

### Adversarial self-checks

- **What if the hole is 3.9 mm, not 4.0?** The draining time rises by (4.0/3.9)² − 1 = 5.2 %, and a C_d fitted with the nominal 4.0 mm comes out 5 % low. The hole is measured, at least two ways, before the prediction.
- **What if the head goes to zero?** √(2gh) → 0, the jet's Reynolds number falls, and below a few millimetres surface tension holds the water back: the model has no business there, and the fit stops at 30 mm.
- **Which assumption fails first?** In experiment 1, a sharp edge (C_d depends on it more than on anything else in the physics). In experiment 2, "the airspeed": the jet of a dryer decays and spreads, so the straw's position matters as much as the reading.
- **What measurement could falsify the model?** A √h–t plot that is clearly curved between 200 and 30 mm (a tapering bottle, or a C_d that changes with Re); upright and tilted manometer readings of the same jet that disagree by more than their uncertainties.
- **Does the result survive a change in method?** C_d from the slope of all marks and from the single 200 → 30 mm time agree; the anemometer and the manometer use different physics (a rotating vane against a pressure).
- **Which decision changes if the uncertainty doubles?** At twice the reading uncertainty the upright tube needs 9.9 m/s for 10 %, so even the dryer's low setting is at the edge, and the fan needs a tilt of about 5°, where the angle's own uncertainty dominates: you would buy the anemometer (enhanced tier) rather than tilt further.

## Check yourself

1. Write the draining law of a bottle with a discharge coefficient and say why √h, not h, is linear in time.
2. What two effects does C_d combine, and what is its typical value for a sharp-edged round hole?
3. Which measured dimension dominates the uncertainty of C_d in this lab, and why?
4. A pitot probe in a free jet is connected to one leg of a U-tube; the other is open to the room. Why does the open leg read the jet's static pressure?
5. At 5 m/s in air of 1.2 kg/m³, what is the dynamic pressure, and what Δh of water does it give?
6. Why does tilting the U-tube improve the resolution, and what limits how far you can tilt it?
7. Three runs give C_d = 0.6313, 0.6324 and 0.6312, each ± 0.032. Why is the run-to-run scatter so much smaller than the stated uncertainty, and which one do you report?
8. Why must the prediction be written and dated before the tape comes off the hole?
