# P07.15  Home Lab: Cooling Curve 🔬

> A warm aluminium bar in still air is a spacecraft thermal model you can hold (with tongs). Predict its time constant from the correlations, with an uncertainty, write it down and date it, and then let the bar tell you whether the correlations were right.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the lumped-capacitance model ρcV dT/dt = −hA(T − T∞), its exponential solution and time constant τ = ρcV/(hA), and the Biot-number condition Bi = hL_c/k < 0.1 under which it holds. 2. Predict τ for a small metal object before measuring, from h = h_conv + h_rad (natural convection from `P07.13`, linearised radiation from `P07.14`) with an uncertainty budget by first-order propagation and a seeded Monte Carlo, and name the dominant term. 3. Measure a cooling curve safely at home (water ≤ 60 °C, tongs, battery or USB power only) and fit it three ways: a weighted line through ln θ, a nonlinear exponential, and a model that includes the thermometer's own lag; read χ² to decide whether one exponential is enough. 4. Compare the lumped model with a one-dimensional distributed model (series solutions for the cylinder and the slab, combined by the product rule) and show where the lumped model fails. 5. Extract h from the fit, compare it with the correlations by a normalised error, and identify which of drafts, emittance, ambient temperature and probe effects explains a discrepancy. |
| **Prerequisites** | `P07.13` (natural convection on a horizontal cylinder and a vertical plate; film temperature), `P07.14` (radiation and its linearisation), `P07.12` (conduction, the heat equation in a solid), `P02.11` (the heat equation and separation of variables), `P02.01` (first-order ODEs), `P02.18` (regression and curve fitting), `P03.05` (ODE integrators, `solve_ivp`), `P04.15` (the 10 kΩ NTC thermistor on a microcontroller), `P00.13` (safety), `P00.14` (measurement, calibration and uncertainty), `P00.15` (the first experiment, normalised error). Skill check: you can fit a weighted straight line with its covariance and run `scipy.optimize.curve_fit` with `absolute_sigma=True`. |
| **Engineering question** | A 25.4 mm diameter, 50 mm long 6061 aluminium bar is warmed in tap-hot water (never above 60 °C), lifted out with tongs, dried and hung in a string sling in still 21 °C room air. Predict before measuring: its time constant from natural convection plus radiation, and whether a lumped model is valid. Then measure and compare. The prediction: τ = 1968 ± 197 s (32.8 min) at a 14 K excess, h = 6.22 W/(m² K), Bi = 1.9 × 10⁻⁴, lumped valid. The SYNTHETIC test (seed 715): τ = 1772 s from the log-linear fit, h = 6.91 ± 0.08 W/(m² K), E_n = 0.55 against the prediction; the distributed model agrees with the lumped one to 0.1 %; but χ²/dof = 34 shows h falls as the bar cools (τ rises from 1616 to 2077 s across the run), and a fit with the correlation's own temperature dependence says the room's convection is 1.15 times the correlation. |
| **Assumptions and limits** | The bar is isothermal (Bi ≪ 0.1, checked) with uniform h on all faces; still air (no drafts, no person nearby, no fan, no sun on the bar); surroundings at the air temperature; ambient constant during the run; the sling and the probe take negligible heat (checked as an extension); properties of 6061-T6 from Lienhard Table A.1 at 20 °C (the alloy of a hardware-store bar is usually 6061 or 6063 and is not always marked); the emittance of bare aluminium is uncertain (Lienhard Table 10.1 spans 0.04 polished to 0.33 heavily oxidised). The correlations reproduce their data to roughly ±10 %, which dominates the prediction. The probe is a first-order lag. |
| **Analytical work** | The lumped energy balance and its solution; the Biot number from comparing internal conduction with surface convection; the linearised radiation coefficient; the time constant's uncertainty budget; the probe-lag solution and why it does not bias τ once its transient has died; the slab and cylinder series solutions with their eigenvalue equations, the product rule for a short cylinder, and the small-Bi limit λ₁² → 2Bi; the energy balance of the fitted model. |
| **Computational deliverable** | `starter/main.py`: `bar_geometry`, `h_conv_bar`, `h_rad_linear`, `biot_number`, `time_constant`, `predict_tau`, `predict_tau_uncertainty`, `fit_loglinear`, `fit_exponential`, `lagged_model`, `fit_lagged`, `h_from_tau`, `en_number`, `cylinder_eigenvalues`, `slab_eigenvalues`, `cylinder_series`, `cylinder_series_mean`, `slab_series`, `short_cylinder`, `lumped_theta`, `simulate_cooling`, `fit_variable_h`, `fitted_model_energy_balance` (`air_props`, `exp_model`, the SYNTHETIC generator `make_synthetic_cooling_data`, `write_synthetic_csv`, `read_csv` and `main` are given; the generator calls your `simulate_cooling`). Tests (18) prove: geometry and τ definitions; the air table nodes; the correlation's zero, monotone growth and ΔT^(1/4)-like scaling; Bi ≪ 0.1 and τ responding to emittance and material; GUM against Monte Carlo; exact recovery of planted exponentials by all three fits; honest uncertainties (pull sd ≈ 1, χ²/dof ≈ 1 on seeded noise); the lag ODE, its fit and the unbiased τ after the transient; a wrong T∞ biasing τ; the variable-h fit recovering a planted factor; eigenvalues solving their equations with the right limits; the series conserving energy and reducing to the lumped model; the bar being lumped to 0.1 %; the nonlinear model's energy balance; the fitted model's energy balance; the SYNTHETIC files labelled, seeded, quantized and read back. |
| **Practical activity** | **physical experiment**: warm a small aluminium bar in tap-hot water no hotter than 60 °C, lift it out with tongs, dry it, hang it in still air, and log its temperature every 15 s for 60 min with a kitchen probe thermometer or the `P04.15` thermistor on a USB- or battery-powered microcontroller, after writing a dated prediction. See `projects/p07-cooling-curve/`. A learner without a kit completes the lesson with the SYNTHETIC dataset in `projects/p07-cooling-curve/data/synthetic/`, labelled as synthetic. |
| **Equipment** | Minimum: an aluminium bar (25 mm round, 50 mm long, or similar), kitchen tongs, a heat-safe plastic or metal (not glass) mug or jug, a kitchen probe thermometer with 0.1 °C resolution, string, tape, a towel, a stand made of books and a pencil or ruler, a phone timer. Recommended: the `p04-electronics-bench` microcontroller and 10 kΩ NTC thermistor for logging, aluminium foil tape, a second thermometer for the room air, a caliper (from `p00-calibration`). Enhanced: a second, black-taped bar for a radiation comparison and a 0.1 °C-accuracy reference thermometer. See `projects/p07-cooling-curve/bom.yaml` (reuses `p00-calibration`, `p00-first-experiment` and `p04-electronics-bench` items). |
| **Safety** | See `projects/p07-cooling-curve/safety.md`. Summary box at the start of part 4. |
| **Evidence** | `data/prediction.md` (dated before measuring), `data/calibration.csv`, `data/raw/*.csv`, the semi-log plot with both fits and the prediction band, the windowed-τ table, h with its uncertainty and E_n, the Biot number, the distributed-model comparison, `cooling_lab.png`, the report. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; all four competencies. Threshold: a dated prediction with uncertainty recorded before the bar went into the water; all public tests pass; τ and h reported with uncertainties and E_n against the prediction; χ² interpreted (one exponential or not); the Biot number and the distributed comparison reported; the safety box completed. |
| **References** | `@lienhard-heat-transfer-textbook` (required-free; https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf): §1.3 Modes of heat transfer, subsection "Heat Convection" (the lumped-capacity solution and its Biot-number condition); 5.2 Lumped-capacity solutions; 5.3 Transient conduction in a one-dimensional slab (Table 5.1 Terms of series solutions for slabs, cylinders, and spheres); 5.5 One-term solutions (Table 5.2); 8.3 Laminar natural convection on a vertical isothermal surface; 8.4 Natural convection in other situations (Natural convection from horizontal isothermal cylinders); §2.3, subsection "Resistance for thermal radiation"; Table 10.1 Total emittances for a variety of surfaces; Table A.1 Properties of metallic solids; Table A.6 Thermophysical properties of gases at atmospheric pressure (101325 Pa). `@mit-16-unified-thermodynamics-propulsion-notes` 18.3 Transient Heat Transfer (Convective Cooling or Heating); 17.3 Dimensionless Numbers and Analysis of Results (required-free; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node129.html; link only). `@mit-ocw-16-050-thermal-energy` Part 3: 6.0 Transient Heat Transfer (Convective Cooling or Heating) (supplementary-free). `@openstax-university-physics-vol2` 1.4 Heat Transfer, Specific Heat, and Calorimetry; 1.6 Mechanisms of Heat Transfer (supplementary-free; https://openstax.org/books/university-physics-volume-2/pages/1-4-heat-transfer-specific-heat-and-calorimetry). |
| **Extensions** | ★ Wrap a second identical bar in matte black electrical tape (emittance about 0.9, typical, illustrative and untested here) and cool both together: the difference in τ measures radiation directly (predicted 18.2 against 32.8 min). ★ Hang the bar vertically and repeat: the side becomes a vertical plate; predict the change first. ★ Estimate the heat drawn by the probe's stem as a fin (`P07.12`) and the probe's own heat capacity, and correct τ for both. |
| **Time** | unmeasured |

## 1. The engineering problem

Every spacecraft and aircraft equipment item has a thermal model, and every thermal model is tested. Before flight, units are heated and cooled in chambers while thermocouples record how fast they respond; the analysts' predicted curves are laid over the measured ones, and the model's conductances and coefficients are corrected until they match. A thermal model that has never been compared with a measurement is a guess with many decimal places.

This lab is the smallest honest version of that test. A short aluminium bar is warmed in tap-hot water, hung in still room air, and its temperature recorded as it cools. If heat flows through the metal much more easily than it leaves the surface, the bar stays nearly isothermal, and a single number, its time constant τ, describes the whole cooling curve. The number can be predicted from what you already have: the bar's dimensions and material, the natural-convection correlation for a horizontal cylinder (`P07.13`), and the radiation coefficient for a gray surface in a room (`P07.14`). The correlations come with uncertainties of order ten per cent, and the emittance of bare aluminium is uncertain by a factor of several, so the prediction carries an honest error bar too.

Then the bar answers back, in three ways. The value of τ says whether the correlations were right for your room. The *shape* of the curve says whether h is really constant (it is not: natural convection weakens as the bar approaches room temperature). And the first minute of the record says how fast your thermometer is. Success criteria, written first: a dated prediction with uncertainty; a τ and an h with uncertainties and E_n against the prediction; a χ² that you can explain; a Biot number below 0.01 and a distributed model that agrees with the lumped one; a completed safety box.

## 2. Understand and derive

### 2.1 The lumped model

Take the bar as the control volume. No work is done; heat leaves through its surface by convection and radiation. The first law with a uniform temperature T:

**ρ c V dT/dt = −h A (T − T∞)**,

with ρ the density, c the specific heat, V the volume, A the surface area and h = h_conv + h_rad the combined coefficient. With θ = T − T∞ and h constant,

**θ(t) = θ₀ exp(−t/τ),  τ = ρ c V/(h A)**.

The **time constant** τ is the heat capacity ρcV divided by the surface conductance hA (`P07.12`'s thermal resistance 1/(hA)): an RC circuit, with the bar as the capacitor (Lienhard 5.2). After one τ the excess has fallen to 36.8 %, after three to 5.0 %. A straight line through ln θ against t has slope −1/τ.

**When is the bar uniform?** Heat must cross the metal, a distance of order L_c = V/A, through a resistance L_c/(kA), before it can cross the surface resistance 1/(hA). Their ratio is the **Biot number**

**Bi = h L_c / k**.

If Bi ≪ 1 the internal temperature differences are a small fraction Bi of the surface-to-air difference, and the lumped model holds; the usual rule is Bi < 0.1 (Lienhard §1.3, "Heat Convection"). For the lesson's bar, L_c = 5.064 mm and k = 167 W/(m K) (6061-T6, Lienhard Table A.1): with h ≈ 6 W/(m² K), Bi = 1.9 × 10⁻⁴. The bar is uniform to two parts in ten thousand.

### 2.2 The bar and its prediction

A 25.4 mm (1 inch) round bar, 50.0 mm long: V = 25.335 cm³, mass 68.4 g, A = 50.03 cm² (the side 39.90 cm² and both ends), ρcV = 61.29 J/K with ρ = 2700 kg/m³ and c = 896 J/(kg K) (Table A.1).

**Convection.** Hung horizontally, the curved side is a horizontal cylinder (Churchill–Chu, Lienhard 8.4, eq. 8.28) and each end is a small vertical plate of height D (8.13a), both with air properties at the film temperature from Table A.6. At a 14 K excess in 21 °C air, Ra_D = 2.097 × 10⁴ (laminar), h_side = 5.29 and h_end = 7.16 W/(m² K), area-weighted **h_conv = 5.67 W/(m² K)**.

**Radiation.** With the exact linearisation of `P07.14`, h_rad = εσ(T² + T∞²)(T + T∞). Bare aluminium is the uncertain input: Lienhard's Table 10.1 gives 0.09 for commercial sheet, 0.04–0.06 polished, 0.20–0.33 heavily oxidised. Take **ε = 0.09 ± 0.04** for an as-bought, unpolished bar (our judgment, which the experiment can test): **h_rad = 0.558 W/(m² K)**, 9 % of the total. A black-taped bar (extension) would have h_rad = 5.58 and τ = 18.2 min.

**The prediction.** h = 6.22 W/(m² K) and **τ = 1968 s = 32.8 min**. Because h_conv depends on ΔT^(1/4), τ depends on where in the run you look: 28.3 min at θ = 28 K, 32.8 min at 14 K, 40.6 min at 5 K. The prediction is made at 14 K, about the mean excess over a run from 50 °C down to 25 °C, and the lesson tests the drift too.

### 2.3 The uncertainty budget

τ is a product and quotient of inputs, but h_conv and h_rad enter as a sum, so `predict_tau_uncertainty` takes central-difference sensitivities rather than exponent rules. Standard uncertainties: D ± 0.05 mm, L ± 0.1 mm (caliper), ε ± 0.04, ρ ± 0.5 %, c ± 1 %, T∞ ± 0.5 K, and **± 10 % on the convection correlation** (an assumed value, typical of the scatter of natural-convection data about Churchill–Chu's curves). The budget, in seconds:

| input | contribution to u(τ) |
|---|---|
| convection correlation (10 %) | 179.2 |
| emittance (± 0.04) | 78.4 |
| specific heat (1 %) | 19.7 |
| density (0.5 %) | 9.8 |
| diameter | 4.0 |
| length | 1.0 |
| ambient temperature | 0.0 |
| **total** | **197 s (10.0 %)** |

A seeded Monte Carlo (4000 draws) gives a mean of 1989 s (1 % above the nominal) and a standard deviation of 207 s (5 % above the first-order value), because τ ∝ 1/h is convex. **Prediction: τ = 1968 ± 197 s (k = 1).** The correlation dominates; the dimensions do not matter; the emittance matters only because it is so poorly known. (The ambient temperature contributes nothing at first order to τ at a given θ; it matters through the fit, 2.6.)

### 2.4 The thermometer's lag

A probe taped to the bar is itself a small body with its own time constant τ_p (seconds to tens of seconds, depending on contact). It starts at room temperature when you attach it, and obeys τ_p dT_m/dt = T − T_m. With T = T∞ + θ₀e^(−t/τ) the solution is

**T_m = T∞ + θ₀ τ/(τ − τ_p) [e^(−t/τ) − e^(−t/τ_p)]**.

The reading rises for the first minute or two, peaks, and then follows the bar with the *same* τ, scaled in amplitude by τ/(τ − τ_p) (here about 1.01). So: **discard the first 5τ_p (about 2 min), and the lag does not bias τ**; or fit the whole record with the lagged model, which also measures τ_p. The tests check both.

### 2.5 Is one exponential enough?

The lumped model with constant h is a straight line on a semi-log plot. With h_conv ∝ θ^(1/4), the line curves: the slope −1/τ shrinks as the bar cools. Between 28 K and 5 K of excess the predicted local τ changes by 43 % (28.3 to 40.6 min). With 0.1 °C readings every 15 s that curvature is far larger than the noise, and χ² will say so. Two remedies: fit windows of the record separately and report τ(θ); or fit a **variable-h model**, ρcV dT/dt = −f h_conv(T)A(T − T∞) − εσA(T⁴ − T∞⁴), solved with `solve_ivp`, whose one free factor f says how much stronger or weaker the room's convection is than the correlation.

### 2.6 What a free ambient temperature does

The nonlinear fit T = T∞ + θ₀e^(−t/τ) can take T∞ as a free parameter. That is useful when you did not measure the room, but T∞ and τ are strongly correlated (a higher T∞ makes the late part of the curve look flatter), and if the model is wrong T∞ absorbs the error. Measure the room air with a second thermometer, a metre from the bar at the bar's height, at the start, middle and end. Setting T∞ 0.5 K too high in the fit shortens τ (the tests check the sign).

### 2.7 The distributed model

To check the lumped model rather than assume it, solve the heat equation in the bar (`P02.11`, `P07.12`). For an infinite cylinder of radius R with convection h at its surface and a uniform initial temperature, separation of variables gives (Lienhard 5.3, Table 5.1)

θ/θ_i = Σ_n C_n exp(−λ_n² Fo) J₀(λ_n r/R),  λ_n J₁(λ_n) = Bi J₀(λ_n),  C_n = 2J₁(λ_n)/[λ_n (J₀² + J₁²)],

with Fo = αt/R² and Bi = hR/k; for a slab of half-thickness l, θ/θ_i = Σ C_n exp(−λ_n² Fo) cos(λ_n x/l) with λ_n tan λ_n = Bi and C_n = 4 sin λ_n/(2λ_n + sin 2λ_n). A short cylinder of length 2l is the intersection of the two, and its solution is their **product** (Lienhard 5.8, "Transient multidimensional heat conduction—the tactic of superposition"). As Bi → 0, λ₁² → 2Bi for the cylinder and Bi for the slab, every higher term dies quickly, and the product becomes exp(−h(2/R + 1/l)t/(ρc)) = exp(−t/τ): the lumped model. The series also conserves energy: d(mean θ)/dFo = −2Bi θ_surface for the cylinder, which the tests check by differentiating the series.

For the bar at h = 6.22 W/(m² K), centre and corner differ by 6.8 × 10⁻⁴ of θ_i after 1 min and 1.1 × 10⁻⁴ after 60 min; the centre and the lumped solution differ by less than 3 × 10⁻⁴. That is 0.01 K: below the thermometer's resolution. The distributed model would matter at Bi = 0.1, where the surface is 2.5 % below the mean and the lumped curve is 5.6 % off by the time the bar has cooled to a mean of θ/θ_i = 0.1 (Problem 5 and `main()`).

### 2.8 Where the model fails

In order of likelihood at home: drafts (a door, a heating vent, a person walking past) raise h_conv by tens of per cent and make it vary; the room's temperature drifts; the emittance is not 0.09; the probe's stem conducts heat out of the bar like a fin, and its mass adds to ρcV (a 5 g steel probe tip adds about 4 %); the sling touches the bar over a significant area; sunlight or a warm lamp falls on the bar. The lumped model and the correlations' functional forms come last.

## 3. Build the computational model

Open `starter/main.py`. The constants (Table A.1 properties, the Table A.6 subset, the nominal emittance and its uncertainty, the 60 °C water limit), `air_props`, `exp_model`, the SYNTHETIC generator and writer, `read_csv` and `main` are given.

### Step 1: the prediction

`bar_geometry`, `h_conv_bar` (both correlations of 2.2 at the film temperature), `h_rad_linear`, `biot_number`, `time_constant`, `predict_tau` (h = f_conv h_conv + h_rad at a representative excess θ) and `predict_tau_uncertainty` (central differences and a seeded Monte Carlo).

### Step 2: the fits

`fit_loglinear(t, T, Tinf, u_T, t_min, theta_min)`: weighted least squares on ln θ, with u(ln θ) = u_T/θ, so the late, small-θ points (where 0.1 °C is a large fraction of θ) get little weight; return τ, its uncertainty from the covariance, χ² and the degrees of freedom. `fit_exponential` uses `curve_fit` with `absolute_sigma=True`, T∞ free or fixed. `lagged_model` is 2.4; `fit_lagged` fits it to the whole record with bounds that keep τ_p positive. `fit_variable_h` wraps `simulate_cooling` in `scipy.optimize.least_squares`. `h_from_tau` and `en_number` are one line each.

### Step 3: the distributed model and the energy balances

`cylinder_eigenvalues` and `slab_eigenvalues` find one root per interval with `brentq` (between successive zeros of J₁ and J₀; between kπ and kπ + π/2). `cylinder_series`, `cylinder_series_mean`, `slab_series`, `short_cylinder` (the product) and `lumped_theta` follow 2.7. `simulate_cooling` integrates the variable-h model with `solve_ivp` (LSODA, rtol 10⁻¹⁰) and two extra states, the heat lost by convection and by radiation, so that ρcV(T₀ − T) = Q_conv + Q_rad can be checked. `fitted_model_energy_balance` compares the stored-energy change of the fitted exponential with the integral of hAθ.

### Verify

`pytest curriculum/phase-07/lesson-15-home-lab-cooling-curve` runs the 18 tests listed in the metadata. Then `python main.py`; it starts

```
bar 25.4 mm x 50.0 mm 6061-T6: V = 25.335 cm^3, mass 68.4 g, A = 50.03 cm^2, Lc = V/A = 5.064 mm, rho c V = 61.29 J/K
PREDICTION (write it down and date it before measuring), at theta = 14 K, air 21.0 C:
  Ra_D = 2.097e+04; h_side 5.29, h_end 7.16, h_conv 5.67 W/m^2 K; h_rad (eps 0.09) 0.558; h 6.22 W/m^2 K
  Bi = 1.89e-04 (lumped if < 0.1); tau = 1968 s = 32.8 min
  tau = 1968 +/- 197 s (GUM, k = 1); Monte Carlo 1989, sd 207 s
```

### Compare

Two independent checks of the theory, before any data. The short-cylinder series agrees with the lumped exponential to 3 × 10⁻⁴ (part 2.7), so the lumped model is not an assumption here but a verified approximation. And the variable-h model closes its energy balance: over 60 min the SYNTHETIC truth's bar gives up 1497.9111 J, and the integrated convective and radiative losses add to 1497.9111 J, radiation carrying 10.6 % of it. The fitted exponential's stored-energy change and integrated loss over 60 min agree to 4.5 × 10⁻¹³ J out of 1464.16 J.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure, the report template and the data formats are in `projects/p07-cooling-curve/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). Read the safety box first, every session. If you have no kit, run the whole reduction on the SYNTHETIC files in `projects/p07-cooling-curve/data/synthetic/` (or on `make_synthetic_cooling_data` with your own seed) and label every table and plot **SYNTHETIC**; you still complete the lesson, but you do not get to claim a measurement.

> **Safety (summary; the authoritative box is `projects/p07-cooling-curve/safety.md`).** Hazards: scalds from hot water and from the hot bar; a dropped bar; water reaching electronics or a socket; a broken container. Hard limits: **tap-hot water no hotter than 60 °C, checked with the thermometer before the bar goes in** (add cold water if it is hotter); **no boiling water, no kettle at the boil, no stove, oven, flame, heat gun, hot plate, microwave heating of metal, or soldering iron**; **the bar is handled with tongs** while it is warm; **battery or USB power only** for any electronics (no mains-powered heaters or circuits); **no glass containers**; **keep water away from the microcontroller and the computer** (water on a tray at the far end of the table). Supervision: none for an adult; children only with an adult handling the water and the bar. Stop criteria: **the bar is too hot to touch through a folded cloth; the water reads above 60 °C; any electronics get wet**; the sling or stand slips; anyone is scalded (cool the skin under running cold water for 20 minutes).

The setup, seen from the side:

```
   pencil or ruler across two stacks of books (stand)
   ══════════════════════════════════════════════
            │  │   cotton string sling (two loops)
            │  │
          ┌─┴──┴──────┐  aluminium bar, 25.4 mm x 50 mm, horizontal
          └─────■─────┘  ■ probe tip or thermistor under foil tape
                │
                └── lead to the thermometer / microcontroller (dry end of the table)
   room thermometer 1 m away at the same height, out of any draught
```

**Prediction first (in `data/prediction.md`, dated).** Measure D and L with a caliper (or a ruler, with its larger uncertainty), weigh the bar (a check on ρ: 68.4 g expected), record the room temperature, and compute h_conv, h_rad, τ ± u(τ) with the budget term by term, the Biot number, and the τ you expect in each third of the run (part 2.5). Say whether the lumped model will hold. Date the file and do not edit it afterwards.

**Procedure in outline** (full version in `procedure.md`). Fill a plastic or metal mug with tap-hot water and measure it: if it is above 60 °C, add cold water until it is not. Put the bar in for five minutes. Hang the sling from the stand and start the room thermometer. Lift the bar out **with tongs**, dry it with a towel (still with tongs and a folded cloth), place it in the sling, tape the probe to its underside with a small square of foil tape, start the clock at the moment the probe touches, and read every 15 s for 60 min without moving, breathing on or standing near the bar. Read the room at 0, 30 and 60 min. Repeat the run once more on another day or at another time of day.

**Raw data format** (`data/raw/YYYY-MM-DD_<tier>_cooling.csv`, one row per reading, never edited):

```
run,t_s,T_C,note
1,0,21.4,probe attached
1,15,36.0,
```

and `data/raw/YYYY-MM-DD_<tier>_ambient.csv` with `run,t_s,T_air_C,note`.

## 5. Validate and improve

Run the reduction on your data (or the SYNTHETIC files) and put the comparisons side by side. The SYNTHETIC dataset (seed 715) was generated with a hidden "truth" that differs from the nominal inputs on purpose: convection 12 % stronger than the correlation (a room with gentle draughts), ε = 0.12, a bar 0.3 mm longer, room air at 21.3 °C, the bar at 49.6 °C when the probe is attached, a probe lag of 20 s, reading noise of 0.04 °C and 0.1 °C quantization. Its first readings, 21.3, 36.1, 43.0, 46.1, 47.4, 47.9 °C, show the lag; the room reads 21.3 °C three times.

| quantity | prediction (before measuring) | SYNTHETIC result (seed 715) | E_n (k = 2) |
|---|---|---|---|
| τ, log-linear, t ≥ 2 min, θ ≥ 1 K [s] | 1968 ± 197 | 1771.7 ± 1.0 (statistical; Birge-scaled ± 6.1) | 0.50 |
| h = ρcV/(τA) [W/(m² K)] | 6.22 ± 0.62 | 6.91 ± 0.08 | 0.55 |
| χ²/dof of one exponential | about 1 if h constant | 7922/231 = 34.3 | model too simple |
| τ in windows 2–20, 20–40, 40–60 min [s] | 1851, 2121, 2399 | 1616 ± 3, 1853 ± 5, 2077 ± 13 | — |
| τ, nonlinear, T∞ fixed at the room [s] | — | 1768.2 ± 1.1 | — |
| τ and T∞, nonlinear, T∞ free | — | 1490.0 ± 2.9 s, 23.09 ± 0.02 °C (correlation −0.94) | T∞ wrong by 1.8 K |
| probe lag τ_p, lagged fit [s] | 5–30 (guess) | 18.2 ± 0.1 (truth 20) | — |
| f_conv, variable-h fit, ε = 0.05 / 0.09 / 0.13 | 1 | 1.196 / 1.151 / 1.106 (χ² 246, 229, 224 for 231 dof) | — |
| Biot number | 1.9 × 10⁻⁴ | 2.1 × 10⁻⁴ | lumped valid |

How to read it. **The prediction passes** (E_n 0.50): the measured τ lies 10 % below the prediction, inside the 10 % uncertainty the correlation carries. **One exponential fails** (χ²/dof = 34): the windows show τ rising by 29 % as the bar cools, in the ratio the correlation predicts (1851 : 2399 predicted, 1616 : 2077 measured, both a rise of about 29 %); the room is consistently more effective than the correlation by a factor of about 0.87 in τ. **A free T∞ is a trap**: it reaches χ² = 556 by moving T∞ 1.8 K above the measured room and shortening τ by 16 %. **The physics-based fit works**: with the correlation's own θ^(1/4) dependence, one factor f fits the record to the noise (χ² ≈ dof); it cannot separate f from ε (they trade along f = 1.196 to 1.106 for ε = 0.05 to 0.13, with χ² slightly preferring the higher ε). The generator's truth, f = 1.12 and ε = 0.12, sits on that line. To separate them, change ε and keep the room (the black-taped bar, Extension 1).

| pattern in your data | likely cause | test |
|---|---|---|
| τ much shorter than predicted, scatter large | drafts, someone nearby | repeat in a closed room; shield the bar with a cardboard box open at the bottom, not touching it |
| τ shorter, curve clean | emittance higher than 0.09 (oxidised or anodised bar), or convection stronger | the black-tape comparison |
| semi-log line bends down at the end | T∞ in the fit too high | re-read the room thermometer; fit with T∞ fixed |
| semi-log line bends up at the end | T∞ too low, or the room warming | ambient readings at 0, 30, 60 min |
| reading rises for a minute at the start | probe lag | discard 5τ_p, or fit the lagged model |
| τ longer than predicted | more heat capacity than assumed (a thicker or longer bar, a steel bar), or air stiller than the correlation's data; probe and sling losses would shorten τ, not lengthen it | weigh the bar; check it with a magnet |

**Improve and re-run.** The biggest term in the prediction budget was the correlation (10 %); the measurement has now replaced it for your room. The next is the emittance: the black-taped bar measures the radiation share directly. Re-run with both bars together and report τ for each.

## 6. Defend your work

Submit: `data/prediction.md` (dated before measuring), `data/calibration.csv`, `data/raw/*.csv`, your analysis script, the semi-log plot with both fits and the prediction, the windowed-τ table, h with its uncertainty and E_n, the Biot number and the distributed comparison, `cooling_lab.png`, and a report following `projects/p07-cooling-curve/report-template.md`. If you used the SYNTHETIC dataset, every table and figure says so in its caption.

### Adversarial self-checks

- **What if the water is only 40 °C?** θ₀ halves; τ rises (h_conv falls as ΔT^(1/4)), and the late part of the record falls into the region where 0.1 °C is a large fraction of θ. The fit becomes noisier, not biased; the water limit costs you signal, not correctness.
- **What if the bar were stainless steel?** k = 13.8 W/(m K) (Table A.1): Bi = 2.5 × 10⁻³ in air, still lumped; ρc is 32 % larger, so τ is longer. In the warm-up water (h of hundreds) Bi approaches 0.2 and the bar is *not* uniform when you lift it out: wait 5 min in the water.
- **Which assumption fails first?** Still air. A door opening during the run is visible as a kink in ln θ; mark it in the `note` column.
- **What measurement could falsify the model?** A semi-log plot whose curvature has the wrong sign (τ falling as the bar cools), or windowed τ that do not follow the correlation's θ^(1/4) trend even after f is fitted, would mean the heat-loss model is wrong (for example, the bar touching a conducting support).
- **Does the result survive a change in method?** τ from the log-linear fit, the nonlinear fit with T∞ fixed and the lagged fit agree within 11 s (1761–1772 s); the free-T∞ fit disagrees by 280 s, which is a statement about that method, not about the bar.
- **Which decision changes if the uncertainty doubles?** With 20 % on the prediction, E_n halves and any room would pass; the experiment would no longer be able to tell a draughty room from a still one. The windowed τ ratio, which does not depend on the correlation's absolute level, would still test the θ^(1/4) law.

## Check yourself

1. Write the lumped energy balance for a cooling body and its solution. What is the time constant in terms of ρ, c, V, h and A?
2. Define the Biot number. What is it for the lesson's bar, and what does that number mean physically?
3. Why is the prediction's uncertainty dominated by the convection correlation, not by the bar's dimensions?
4. Why must the water be checked with the thermometer before the bar goes in, and what is the limit?
5. A probe with τ_p = 20 s is attached to the bar at room temperature. Which readings do you discard, and why does the lag not bias τ after that?
6. Your log-linear fit has χ²/dof = 34. What does that tell you, and what do you do next?
7. Why is leaving T∞ free in the nonlinear fit risky?
8. Why must the prediction be written and dated before the bar goes into the water?
