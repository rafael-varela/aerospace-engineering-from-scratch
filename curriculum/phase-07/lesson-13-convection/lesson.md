# P07.13  Heat Transfer: Convection

> A correlation is a measurement in disguise: it tells you h to ten or twenty per cent, inside the range where it was measured, and nothing outside it. Know its range, check it two ways, and design with margin.

| | |
|---|---|
| **Learning outcomes** | 1. Define the Reynolds, Prandtl, Nusselt, Stanton, Grashof and Rayleigh numbers, say what each compares, and evaluate them with air and water properties interpolated from Lienhard's Tables A.6 and A.3 at the film temperature. 2. Use the flat-plate (laminar, turbulent, mixed), cylinder cross-flow (Churchill–Bernstein), pipe-flow (fully developed laminar, Gnielinski, Dittus–Boelter) and natural-convection (Churchill–Chu, Raithby–Hollands, Fujii–Imura) correlations inside their stated ranges. 3. Derive the laminar flat-plate results average Nu = 2 × local Nu at the trailing edge and St·Pr^(2/3) = c_f/2, and show by computation that the turbulent power laws obey the analogy only to the factor Pr^(0.6 − 1/3). 4. Implement and test the correlations, check the mixed laminar–turbulent average against a numerical average of the local coefficient, and find limits and discontinuities. 5. Decide whether an avionics box needs a fan and how much coolant a cold plate needs, with the uncertainty that the correlations carry. |
| **Prerequisites** | `P07.12` (Fourier's law, thermal resistance, the convection boundary condition), `P07.01` (properties and the ideal gas), `P03.02` (root finding), `P03.03` (interpolation), `P03.04` (numerical integration), `P01.05` (integration), `P00.14` (uncertainty). Skill check: you can solve f(T) = 0 with `scipy.optimize.brentq` and interpolate a table with `numpy.interp`. |
| **Engineering question** | A 60 W avionics box, 0.30 × 0.12 × 0.20 m, painted black (ε = 0.90) and standing on a rack tray in still 25 °C air, has a 50 °C surface limit: does it need a fan, and how fast must the air move? A 400 W cold plate has a 4 mm, 1 m water channel with a 20 °C inlet: what flow keeps its hottest wall below 40 °C? The answers: 51.8 °C without a fan (radiation carries 56 % of the heat; unpainted, 77.5 °C), so it needs 1.06 m/s of air; the cold plate needs at least 8.8 g/s (Re 3171, in the transitional range), and a design at 15 g/s (Re 5134, wall 31.6 °C) is the one to defend. |
| **Assumptions and limits** | Steady state; isothermal box faces; constant properties evaluated at the film temperature (Tw + T∞)/2; ideal gas (β = 1/T); the correlations' stated ranges (Pr ≥ 0.6 for the laminar plate; 2300 ≤ Re ≤ 5 × 10⁶ for Gnielinski; Ra ≤ 10⁹ for the laminar Churchill–Chu forms; 1 < Ra* < 10⁷ for the Raithby–Hollands upper-face form); a sharp transition at Re = 5 × 10⁵ on a plate (real transition depends on free-stream turbulence and is not sharp); fully developed pipe flow (entry effects raise h, so ignoring them is conservative); no mixed convection (valid when Gr/Re² ≪ 1); radiation from a small gray body to large surroundings. Correlations reproduce their data to roughly ±10–20 %, and up to a factor of two in transition. |
| **Analytical work** | The dimensionless groups from the energy balance on a surface; the laminar flat plate: local and average Nu and c_f, the factor of two, the Reynolds–Colburn analogy; the mixed average by integrating the local law in two pieces; natural convection scalings Ra^(1/4) (laminar) and Ra^(1/3) (turbulent, h independent of size); the bulk-temperature energy balance and the wall temperature of a uniformly heated channel; the linearised radiation coefficient. |
| **Computational deliverable** | `starter/main.py`: `air_properties`, `water_properties`, `reynolds`, `rayleigh`, `h_from_nu`, `stanton`, the plate correlations (`plate_laminar_local_nu`, `plate_laminar_average_nu`, `plate_laminar_cf_local`, `plate_turbulent_local_nu`, `plate_turbulent_cf_local`, `plate_local_nu`, `plate_average_nu`, `plate_average_nu_numerical`), `cylinder_crossflow_nu`, `filonenko_f`, `gnielinski_nu`, `dittus_boelter_nu`, `pipe_nu`, `vertical_plate_nu`, `vertical_plate_nu_allra`, `horizontal_cylinder_nu`, `horizontal_plate_up_nu`, `horizontal_plate_down_nu`, `h_natural`, `h_radiation_linear`, `box_heat_loss_natural`, `box_surface_temperature`, `box_heat_loss_forced`, `box_required_speed`, `cold_plate` (the property tables, `box_faces` and `main` are given). Tests (18) prove: table nodes reproduced and columns consistent; the laminar Reynolds–Colburn analogy exact for three Prandtl numbers; the turbulent analogy off by exactly Pr^(0.6 − 1/3); average = 2 × local; the mixed average equal to numerical integration to 10⁻⁷ and continuous at transition; √Re scaling; Churchill–Bernstein's limit 0.3 and its low-Re form; pipe limits, range checks and the Dittus–Boelter/Gnielinski agreement for gases; natural-convection limits and scalings; hot/cold symmetry; the radiation linearisation; and the energy balances of the box, the fan and the cold plate. |
| **Practical activity** | **computational experiment**: a power sweep and a fan-speed requirement for the avionics box; the sensitivity of the answer to the property temperature, the correlation's ±20 % and the paint; a flow sweep of the cold plate across the laminar–turbulent jump; a correlation comparison (Gnielinski against Dittus–Boelter for air and water, Churchill–Chu laminar against all-Ra); and a water-cooled thrust-chamber channel at 2 MW/m². |
| **Equipment** | None. A computer with the environment from `P00.01`. (Lesson `P07.15` measures natural convection plus radiation on a small bar at home.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `convection.png`, the derivations of the factor of two and the laminar analogy, the experiment tables, and a half-page note answering the engineering question with the margin you would carry and why. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the dimensionless groups defined and the laminar results derived; all public tests pass; the box temperature, fan speed and cold-plate flow reproduced; the note names the correlation range used for each answer and the uncertainty it carries. |
| **References** | `@lienhard-heat-transfer-textbook` (required-free; https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf): Ch. 6 "Laminar and turbulent boundary layers": 6.5 Heat transfer coefficient for laminar, incompressible flow over a flat surface; 6.6 The Reynolds-Colburn analogy; 6.8 Heat transfer in turbulent boundary layers (Power-law correlations for local turbulent heat transfer); 6.9 Turbulent transition and overall heat transfer (Average Nusselt number for the entire plate); Ch. 7 "Forced convection in a variety of configurations": 7.2 Heat transfer to or from laminar flows in pipes; 7.3 Turbulent pipe flow; 7.6 Heat transfer during cross flow over cylinders; 7.7 Finding and assessing correlations for other configurations; Ch. 8 "Natural convection in single-phase fluids and during film condensation": 8.3 Laminar natural convection on a vertical isothermal surface; 8.4 Natural convection in other situations (Natural convection from horizontal isothermal cylinders; Laminar heat transfer from inclined and horizontal plates); Table A.6 Thermophysical properties of gases at atmospheric pressure (101325 Pa); Table A.3 Thermophysical properties of saturated liquids. `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 17 "Convective Heat Transfer": 17.1 The Reynolds Analogy; 17.3 Dimensionless Numbers and Analysis of Results (required-free; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node122.html; link only). `@mit-ocw-16-050-thermal-energy` Part 3: 3.0 Convective Heat Transfer; 3.1 The Reynolds Analogy (supplementary-free). `@openstax-university-physics-vol2` 1.6 Mechanisms of Heat Transfer (Convection) (supplementary-free; the physical picture). |
| **Extensions** | ★ Replace the sharp plate transition with Lienhard's blended local formula (6.117) and integrate it numerically with `plate_average_nu_numerical`; compare with the sharp model at Re_L = 10⁶. ★ Add mixed convection to the box: combine forced and natural Nu as (Nu_F³ + Nu_N³)^(1/3) and find the speed below which the fan is pointless. ★ Give the cold plate an entry-length correction and a wall-conduction resistance (`P07.12`) between the heat source and the channel. |
| **Time** | unmeasured |

## 1. The engineering problem

Every watt an avionics box dissipates has to leave through its skin, and the last step is almost always convection: heat crossing from a solid surface into a moving fluid. In a pressurised equipment bay with no fan, the fluid moves only because warm air is lighter than cool air; that is natural convection, and it is weak. With a fan, or in a cold plate where water flows through channels, it is forced convection, and it can be a thousand times stronger. Rocket engines push this to its limit: the propellant flowing through the channels of a regeneratively cooled chamber carries away megawatts per square metre, and the wall survives only because the coolant's heat transfer coefficient is enormous.

The quantity that matters is the heat transfer coefficient h in q = h (T_w − T∞), in W/(m² K). It is not a property of the fluid. It depends on the flow, the shape, the size and the temperature difference, and for the flows engineers meet it cannot be computed from first principles in closed form. It comes from **correlations**: formulas fitted to careful experiments, written in dimensionless groups so that one experiment on a small plate in air can predict a large plate in water. Using them well means three things: computing the groups with the right properties, staying inside the range where each correlation was fitted, and knowing how far to trust the answer.

This lesson builds a library of the standard correlations, tests them against each other and against two exact results for the laminar flat plate, and uses them on three aerospace problems: an avionics box cooled by still air (and by radiation, which turns out to matter as much), the same box with a fan, and a cold plate with a water channel. Success criteria, written first: every node of the property tables is reproduced exactly; the laminar Reynolds–Colburn analogy holds to round-off; the mixed-flow plate average equals a numerical integral of the local coefficient to 10⁻⁷; every application closes its energy balance; and each design answer comes with the correlation range it used.

## 2. Understand and derive

### 2.1 What h depends on, in dimensionless form

Next to a cooled wall the fluid is at rest (it sticks to the surface), so the heat leaves the wall by conduction into the fluid: q = −k_f ∂T/∂y at y = 0. Writing that equal to h(T_w − T∞) and making y dimensionless with a length L gives the **Nusselt number**

**Nu = h L / k_f**,

the wall temperature gradient in units of (T_w − T∞)/L. Nu = 1 means pure conduction through a fluid layer of thickness L; large Nu means the flow has thinned the layer. Nu depends on two other groups. For forced flow:

- **Reynolds number** Re = U L / ν, inertia against viscosity (ν = μ/ρ, the kinematic viscosity): it sets whether the flow is laminar or turbulent and how thin the boundary layer is.
- **Prandtl number** Pr = ν/α = μ c_p / k, momentum diffusivity against thermal diffusivity (α = k/(ρ c_p)): a fluid property, about 0.71 for air and 5.85 for water at 300 K. When Pr > 1 the thermal boundary layer is thinner than the velocity boundary layer.

For natural convection the velocity is not given; buoyancy makes it. A fluid of thermal expansion coefficient β (= 1/T for an ideal gas) heated by ΔT is lighter by ρβΔT, and the ratio of buoyancy to viscous forces is the **Grashof number** Gr = g β ΔT L³ / ν². Its product with Pr is the **Rayleigh number**

**Ra = g β ΔT L³ / (ν α)**,

the group natural-convection correlations use. Finally the **Stanton number** St = Nu/(Re Pr) = h/(ρ c_p U) compares the heat actually transferred with the enthalpy the stream carries past.

Where are the properties evaluated? They change across the boundary layer. The usual compromise is the **film temperature** T_f = (T_w + T∞)/2. The code interpolates Lienhard's Table A.6 (air at 1 atm) and Table A.3 (saturated liquid water) linearly in T; between 290 and 340 K the air's properties change by a few per cent per 10 K, so the interpolation error is far below the correlations' own scatter. At 315 K, halfway between nodes, it gives ν = 1.7180 × 10⁻⁵ m²/s, k = 0.02750 W/(m K) and Pr = 0.7055.

### 2.2 The laminar flat plate: two exact results

For a flat plate with a uniform wall temperature in laminar flow, boundary-layer theory (Lienhard 6.5) gives the **local** Nusselt number at distance x from the leading edge,

**Nu_x = h_x x / k = 0.332 Re_x^(1/2) Pr^(1/3)**,  Pr ≥ 0.6,

and the local skin-friction coefficient c_f = τ_w/(½ρU²) = 0.664 Re_x^(−1/2) (the Blasius solution). Two results follow by algebra.

**The factor of two.** h_x = 0.332 k Pr^(1/3) (U/ν)^(1/2) x^(−1/2) falls as x^(−1/2): the boundary layer thickens downstream. The average over a plate of length L is

h̄ = (1/L) ∫₀ᴸ h_x dx = (1/L) · 0.332 k Pr^(1/3) (U/ν)^(1/2) · 2L^(1/2),

so **Nu_L = h̄ L/k = 0.664 Re_L^(1/2) Pr^(1/3) = 2 Nu_x(x = L)**. The average is twice the trailing-edge value because x^(−1/2) integrates to 2x^(1/2). At Re = 10⁵ in air the code gives Nu_x = 93.53 and Nu_L = 187.06.

**The Reynolds–Colburn analogy** (Lienhard 6.6). Divide Nu_x by Re_x Pr^(1/3):

St_x Pr^(2/3) = Nu_x / (Re_x Pr^(1/3)) = 0.332 Re_x^(−1/2) = c_f/2.

**St Pr^(2/3) = c_f/2** links heat transfer to wall friction: measure one and you know the other. For the laminar plate it is exact, because both results come from the same similarity solution: at Re = 10⁵ both sides are 1.049876 × 10⁻³. It fails for a uniform-heat-flux wall (Lienhard's constant there is 0.4587 instead of 0.332), with pressure gradients, and for liquid metals.

### 2.3 Turbulent and mixed flow on a plate

Beyond a transition Reynolds number the boundary layer becomes turbulent and h jumps. For gases Lienhard (6.8) gives the power law **Nu_x = 0.0296 Re_x^0.8 Pr^0.6**, and the matching one-fifth-power friction law is c_f = 0.0592 Re_x^(−0.2). The analogy then gives St Pr^(2/3)/(c_f/2) = Pr^(0.6 − 1/3) exactly: 1 at Pr = 1, 0.9117 for air. A turbulent boundary layer carries heat and momentum by the same eddies, so the analogy holds well, but the Prandtl exponents fitted to data are not the laminar 1/3.

Take transition as sharp at Re_c = 5 × 10⁵ (a convention; a plate in a turbulent free stream can trip much earlier). At Re_x = Re_c the local Nu jumps from 209.1 (laminar) to 871.2 (turbulent): a factor 4.2. The average over a plate longer than the laminar run is the laminar integral up to x_c plus the turbulent one beyond it:

**Nu_L = 0.037 Pr^0.6 (Re_L^0.8 − Re_c^0.8) + 0.664 Re_c^(1/2) Pr^(1/3)**,  Re_L > Re_c,

(0.037 = 0.0296/0.8 from integrating x^(−0.2)). This is Lienhard's (6.120) with the transition zone shrunk to a point. It is continuous at Re_L = Re_c even though the local value jumps: the average is an integral, and integrals of step functions are continuous. The code checks it against `plate_average_nu_numerical`, which integrates the local law with adaptive quadrature: 264.54, 418.27 and 2630.51 at Re_L = 2 × 10⁵, 5 × 10⁵ and 2 × 10⁶, identical to the digits printed.

### 2.4 Cylinders and pipes

**A cylinder in cross-flow** (a probe, a wire, a tube in a heat exchanger) has a separated wake, and no simple theory. Churchill and Bernstein's single formula covers the whole range Re Pr > 0.2 (Lienhard 7.6, eq. 7.65):

Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3) / [1 + (0.4/Pr)^(2/3)]^(1/4) · [1 + (Re/282 000)^(5/8)]^(4/5).

As Re → 0 it tends to 0.3; below Re ≈ 4000 the last bracket is close to 1 (5.5 % at 4000) and the formula is the simpler (7.66). For a 6 mm tube in 300 K air it gives h = 43.5, 140.1 and 338.7 W/(m² K) at 1, 10 and 50 m/s.

**Pipes.** Far from the inlet, laminar flow has a fixed velocity and temperature profile and Nu_D = hD/k is a constant: **4.364 for a uniform wall heat flux, 3.657 for a uniform wall temperature** (Lienhard 7.2). Turbulent flow (Lienhard 7.3) uses **Gnielinski's correlation**

Nu_D = (f/8)(Re − 1000) Pr / [1 + 12.7 (f/8)^(1/2) (Pr^(2/3) − 1)],  f = (1.82 log₁₀ Re − 1.64)^(−2),

valid for 2300 ≤ Re ≤ 5 × 10⁶ and 0.6 ≤ Pr ≤ 10⁵ (f is Filonenko's smooth-pipe friction factor), or the older **Dittus–Boelter** form Nu = 0.023 Re^0.8 Pr^0.4. The code refuses Gnielinski outside its range. Note what happens at Re = 2300: laminar gives 4.364, Gnielinski 7.22 for air and 14.57 for water. Between them is the transition range, where the flow can be either and h is uncertain by a factor of two; designers keep out of it.

### 2.5 Natural convection

On a **vertical plate** of height L, Churchill and Chu's laminar correlation (Lienhard 8.3, eq. 8.13a) is

**Nu_L = 0.68 + 0.670 Ra_L^(1/4) / [1 + (0.492/Pr)^(9/16)]^(4/9)**,  Ra_L ≤ 10⁹,

and their all-Ra form (8.13b), Nu = {0.825 + 0.387 Ra^(1/6)/[1 + (0.492/Pr)^(9/16)]^(8/27)}², is used above. The Ra^(1/4) is the laminar boundary-layer scaling: Nu − 0.68 doubles when Ra grows sixteen-fold, so h ∝ ΔT^(1/4) L^(−1/4). The two forms do not meet at 10⁹: in air they give 92.1 and 122.9, a 33 % jump. That is not a bug; it is the transition to turbulence, where the data scatter and the laminar form falls low. A **horizontal cylinder** of diameter D uses (8.28) Nu_D = 0.36 + 0.518 Ra_D^(1/4)/[1 + (0.559/Pr)^(9/16)]^(4/9).

**Horizontal plates** behave differently on their two faces (Lienhard 8.4, "Laminar heat transfer from inclined and horizontal plates"). Above a hot plate the warm air must rise from the middle as a plume; the flow is unstable and vigorous. With the length scale L* = A/P (area over perimeter), Raithby and Hollands give Nu = 0.560 Ra*^(1/4)/[1 + (0.492/Pr)^(9/16)]^(4/9) for 1 < Ra* < 10⁷, corrected to 1.4/ln(1 + 1.4/Nu) when Nu < 10 (thick boundary layers), and Nu = 0.14 Ra^(1/3)(1 + 0.0107 Pr)/(1 + 0.01 Pr) above 10⁷. With Nu ∝ Ra^(1/3), h ∝ (L³)^(1/3)/L does not depend on the plate's size: turbulent natural convection is local. Under a hot plate the warm air is trapped and must creep to the edges: Fujii and Imura's Nu = 0.58 Ra^(1/5) (for 10⁶ ≤ Ra ≤ 10¹¹) is much smaller. A cold surface swaps the two faces.

### 2.6 Radiation alongside

A box in still air also radiates. For a small gray body of emittance ε in large surroundings at T_sur (lesson `P07.14` derives it),

q_rad = ε σ (T_w⁴ − T_sur⁴) = h_rad (T_w − T_sur),  h_rad = ε σ (T_w² + T_sur²)(T_w + T_sur).

For black paint (ε = 0.90, Lienhard Table 10.1, "black gloss") near 40 °C, h_rad is 6.18 W/(m² K): larger than the natural-convection coefficient. Natural convection and radiation are the same size near room temperature, which is why spacecraft engineers and avionics engineers both care about paint.

### 2.7 Worked problem: the avionics box

The box is 0.30 m deep, 0.12 m wide and 0.20 m high, its bottom on a rack tray treated as adiabatic (heat lost into the tray is left out, which is conservative). Its four vertical faces (2 × 0.20 × 0.30 m² and 2 × 0.20 × 0.12 m², height 0.20 m) use the vertical-plate correlation; the 0.036 m² top uses the upper-face correlation with L* = 0.036/0.84 = 0.0429 m. Total area 0.2040 m². For a surface temperature T_w, the heat it loses is

Q(T_w) = Σ_faces h_i(T_w) A_i (T_w − T∞) + h_rad(T_w) A (T_w − T∞),

and the operating temperature is the root of Q(T_w) = 60 W (`brentq`). With black paint the root is **T_w = 51.8 °C**: the sides see Ra = 1.673 × 10⁷ (laminar) and h = 4.57 W/(m² K), the top Ra* = 1.646 × 10⁵ and h = 5.93 W/(m² K), the area-weighted h_conv = 4.81 and h_rad = 6.18 W/(m² K), and radiation carries 56 % of the heat. Unpainted (ε = 0), the same box would run at 77.5 °C. Check by hand: (4.81 + 6.18) × 0.204 × (51.8 − 25) = 60 W.

The box misses its 50 °C limit by 1.8 K. A fan blowing along its 0.30 m depth turns the sides and top (0.132 m²) into flat plates in forced flow (the ends and the tray are left out, again conservatively). The required speed is the root of Q_forced(U) = 60 W at T_w = 50 °C: **U = 1.06 m/s**, Re_L = 18 965 (laminar), h = 7.37 W/(m² K). Is natural convection still important at that speed? Gr_H/Re_L² = 0.062: buoyancy is a small correction, so treating the flow as forced is consistent.

### 2.8 Worked problem: the cold plate

400 W enter water flowing at ṁ through a 4 mm channel 1 m long, uniformly along its length, with a 20 °C inlet. The **first law on the coolant** fixes the outlet: T_out = T_in + Q/(ṁ c_p). The wall is hottest at the outlet, where the water is warmest: T_w,out = T_out + q″/h with q″ = Q/(πDL) = 31 831 W/m². The code iterates c_p on the mean bulk temperature and evaluates h from `pipe_nu`. At 5 g/s the flow is laminar (Re 1973, Nu 4.36, h 671 W/(m² K)) and the wall reaches 86.6 °C; at 10 g/s it is turbulent (Re 3549, h 4011 W/(m² K)) and the wall is at 37.5 °C. The minimum flow for a 40 °C wall is **8.80 g/s** (Re 3171, outlet 30.87 °C), squarely in the range where Gnielinski is least reliable and where a slight change could drop the flow back to laminar. A defensible design runs at 15 g/s (Re 5134, wall 31.6 °C).

### 2.9 Where the model stops

- **Correlation scatter.** Gnielinski reproduces 90 % of liquid data within ±20 %; the natural-convection correlations a few per cent for clean laminar cases, much worse near transition. Experiment 3 shows ±20 % on h_conv moves the box by about ±2 K.
- **Transition.** Plate transition is not at a fixed Re; pipe flow between 2300 and about 10⁴ is unpredictable; natural convection near Ra = 10⁹ has a 33 % jump between forms.
- **Mixed convection.** At low fan speeds forced and natural convection add (Extension 2).
- **Real boxes.** Fins, vents, internal circulation, conduction into the tray, neighbouring boxes that block radiation (`P07.14`).
- **Variable properties.** Large temperature differences (a rocket coolant channel) need property-ratio corrections the film-temperature rule does not supply.

## 3. Build the computational model

Open `starter/main.py`. The property tables (copied from Lienhard Tables A.6 and A.3 with their source in a comment), `box_faces` and `main` are given.

### Step 1: properties and groups

`air_properties(T)` and `water_properties(T)` interpolate each column with `np.interp` and raise `ValueError` outside the table: extrapolating a property table silently is how a correlation gets used where nobody measured. Return β = 1/T for air (ideal gas) and the table's β for water. `reynolds`, `rayleigh`, `h_from_nu` and `stanton` are one line each.

### Step 2: the plate

The laminar and turbulent local and average laws of 2.2 and 2.3 are one line each; write them with `np.sqrt`, `np.cbrt` and `**` so they accept arrays. `plate_local_nu` uses `np.where` for the sharp transition; `plate_average_nu` is the mixed formula. `plate_average_nu_numerical` is the independent check: in s = x/L the average is ∫₀¹ Nu_x(Re_L s)/s ds; give `quad` the transition point as a breakpoint, or it will step over the jump.

```python
f = lambda s: float(plate_local_nu(Re_L * s, Pr, Re_crit)) / s
s_tr = Re_crit / Re_L
return quad(f, 0, s_tr)[0] + quad(f, s_tr, 1)[0]      # when s_tr < 1
```

### Step 3: cylinders, pipes, natural convection

Each correlation of 2.4 and 2.5 is a formula; `gnielinski_nu` raises outside 2300–5 × 10⁶. `h_natural(geometry, Tw, Tinf, L)` evaluates properties at the film temperature, computes Ra with β = 1/T_f, swaps 'up' and 'down' for a cold surface, and returns h.

### Step 4: the applications

`box_heat_loss_natural` sums `h_natural × area × ΔT` over `box_faces` and adds radiation with `h_radiation_linear`; `box_surface_temperature` is `brentq` on it. `box_heat_loss_forced` uses `plate_average_nu` on the sides and top; `box_required_speed` is `brentq` in U, returning 0 if natural convection already suffices. `cold_plate` is 2.8.

Running `main()` prints (abridged)

```
flat plate, Re_x = 1e+05, Pr = 0.707: laminar Nu_x 93.53, average Nu_L 187.06 (= 2 x local)
Reynolds-Colburn, laminar: St Pr^(2/3) = 1.049876e-03, c_f/2 = 1.049876e-03
Reynolds-Colburn, turbulent Re_x = 1e+06: St Pr^(2/3) / (c_f/2) = 0.9117 (= Pr^(0.6 - 1/3) = 0.9117)
plate average Re_L = 2e+06: closed form 2630.51, numerical average of local 2630.51
  area 0.2040 m^2; convection only: Tw = 77.5 C (h_conv 5.61 W/m^2 K)
  black paint eps = 0.9: Tw = 51.8 C, h_conv 4.81, h_rad 6.18 W/m^2 K, radiation carries 56 %
  to hold Tw = 50 C: air speed 1.060 m/s (Re_L 18965, h 7.37 W/m^2 K)
```

and writes `convection.png`: the plate laws, the box temperature against power with and without radiation, and the two pipe correlations.

### Verify

`python -m pytest curriculum/phase-07/lesson-13-convection` runs 18 tests: property nodes reproduced exactly and the columns consistent (ν = μ/ρ to 0.3 %, Pr = μc_p/k to 1 %); the laminar analogy to 10⁻¹² at Pr = 0.7, 1 and 5.85; the turbulent analogy off by exactly Pr^(0.6 − 1/3); average = 2 × local; the mixed average equal to numerical integration to 10⁻⁷ at five Re_L, continuous at Re_c while the local value jumps; √Re scaling of h; Churchill–Bernstein → 0.3 and its low-Re form; the pipe constants, the range check, Filonenko against Blasius' 0.316 Re^(−1/4) within 6 %; Dittus–Boelter within 20 % of Gnielinski for air; Churchill–Chu → 0.68, the cylinder → 0.36, the Ra^(1/4) and Ra^(1/3) scalings; the hot/cold face swap; the exact radiation linearisation; and the box, fan and cold-plate energy balances.

### Compare

Two independent routes agree on every check the lesson relies on. The mixed-flow average comes from a closed-form integral and from adaptive quadrature of the local law; they agree to 10⁻⁷, which validates the algebra of (6.120) and the code that implements it. The laminar Reynolds–Colburn identity checks the heat-transfer law against the friction law. For pipes, two independent correlations fitted to different data sets bracket the answer: Dittus–Boelter is 6–12 % above Gnielinski for air but 16–28 % below it for water at Re ≥ 10⁵ (Experiment 4). That spread is the real uncertainty of turbulent pipe-flow h, and no amount of arithmetic removes it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **The box against its power.** Q = 20–100 W, black paint. *Predict* whether doubling the power doubles the temperature rise. *Finding:* 35.3, 43.9, 51.8, 59.1, 66.0 °C at 20, 40, 60, 80, 100 W. The rise (10.3 K at 20 W, 41.0 K at 100 W) grows less than linearly, because both h_conv (∝ ΔT^(1/4)) and h_rad (∝ T³) grow as the box warms.
2. **How much fan?** The speed that holds the black box at 40, 45, 50 and 60 °C at 60 W. *Predict* the shape of the curve. *Finding:* 6.31, 2.54 and 1.06 m/s, and none at 60 °C (natural convection and radiation suffice). Each 5 K of margin costs a factor of about 2.5 in speed near the limit, because forced convection must take over what radiation and natural convection lose as the box cools.
3. **What the answer is sensitive to.** At 60 W, evaluate the properties at the ambient and the wall temperature instead of the film; scale h_conv by 0.8 and 1.2; change ε. *Finding:* 51.63, 51.76 and 51.86 °C for ambient, film and wall properties (convection only: 77.04, 77.46, 78.51 °C): negligible. ±20 % on h_conv gives 53.93 and 49.92 °C: about ±2 K. The paint: 77.46, 72.64, 59.03, 51.76 and 50.80 °C for ε = 0, 0.09 (bare aluminium sheet, Table 10.1), 0.5, 0.90 and 0.97. The emittance is the biggest lever on this box, and a bare aluminium box is 21 K hotter than a black one.
4. **Correlation against correlation.** Dittus–Boelter over Gnielinski at Re = 10⁴, 3 × 10⁴, 10⁵, 10⁶: 1.061, 1.083, 1.115, 1.109 for air (Pr 0.707); 0.996, 0.906, 0.843, 0.724 for water (Pr 5.85). Churchill–Chu at Ra = 10⁹: 92.1 laminar, 122.9 all-Ra. *Predict* first which correlation you would trust for a water-cooled chamber.
5. **The cold plate across transition.** ṁ = 5–30 g/s. *Finding:* hottest wall 86.60, 37.50, 31.61, 28.75, 25.92 °C at 5, 10, 15, 20, 30 g/s (Re 1973 → 9903); minimum flow for 40 °C, 8.80 g/s at Re 3171. The 49 K drop between 5 and 10 g/s is the laminar–turbulent jump in h (671 to 4011 W/(m² K)); the real curve between them is uncertain.
6. **A water-cooled thrust-chamber channel.** Test stands often cool a chamber with water rather than propellant (regenerative-style channels, without the propellant). A 3 mm channel, water at 300 K, a heat flux of 2 MW/m² (typical, illustrative for a small chamber's throat region). *Predict* the velocity that keeps the wall-to-water difference below 50 K. *Finding:* film ΔT = q″/h of 185.4, 80.1, 43.8, 30.9 and 24.1 K at 2, 5, 10, 15 and 20 m/s; 50 K needs 8.58 m/s (Re 30 052, Gnielinski), where Dittus–Boelter would predict 55.2 K. A real chamber channel must also stay below the coolant's local boiling point at its pressure, which this model does not check.

Raw data format: CSV tables `Q_W, Tw_C` (1); `T_limit_C, U_mps` (2); `case, Tw_C` (3); `Re, Pr, Nu_Gn, Nu_DB` (4); `mdot_gps, Re, Nu, h, T_out_C, T_wall_C` (5); `U_mps, Re, h, dT_film_K` (6); plus `convection.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| laminar St Pr^(2/3) vs c_f/2, Re = 10⁵ | equal | 1.049876 × 10⁻³ both | < 10⁻¹² | identity |
| turbulent ratio, air | Pr^(0.6 − 1/3) = 0.9117 | 0.9117 | 0 | exponent mismatch, by construction |
| average/local, laminar | 2 | 2 | < 10⁻¹² | — |
| mixed average, Re_L = 2 × 10⁶ | closed form 2630.51 | quadrature 2630.51 | < 10⁻⁷ relative | — |
| box energy balance at T_w | 60 W | 60 W | < 10⁻⁹ | root-finding tolerance |
| cold plate first law | ṁ c_p ΔT = 400 W | 400 W | < 10⁻⁹ | — |
| pipe Nu, water, Re = 10⁵ | Gnielinski | Dittus–Boelter 16 % lower | 16 % | two correlations, two data sets |

The arithmetic is checked to round-off; the physics is checked only to the scatter of the correlations, and that scatter decides the design margin. The engineering answers: **the black box runs at 51.8 °C in still air, 1.8 K over its limit; a 1.06 m/s air stream holds it at 50 °C, but given ±2 K from the correlation alone, a designer would ask for 2.5 m/s (45 °C) or a lower-power design. The cold plate needs 8.8 g/s in principle and 15 g/s in practice, to stay clear of transition.** The unpainted box (77.5 °C) shows that the cheapest cooling improvement is the paint.

What to improve: the box model leaves out the tray, internal circulation and neighbours; the next step is a thermal network with a conduction path to the tray (`P07.12`) and a view-factor model of the bay (`P07.14`). What not to improve: the property interpolation, which moves the answer by 0.2 K.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of 2.2 (the factor of two and the analogy) and of the mixed average; `convection.png`; the six experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question, naming the correlation and its range for each number, and stating the margin you would carry.

### Adversarial self-checks

- **What if the box power doubles to 120 W?** It runs at 72.7 °C in still air; the vertical faces stay laminar (Ra grows only as ΔT), so the same correlations apply. Holding 50 °C then takes 10.1 m/s of air, and a cold plate becomes the better answer.
- **What if the velocity goes to zero?** Forced-convection h → 0 as U^(1/2), which is wrong: natural convection takes over. At 0.5 m/s the forced model already gives less convection (19.7 W at 50 °C) than still air does (24.1 W), so at low speeds neither correlation alone is right (mixed convection, Extension 2). `box_required_speed` returns 0 when natural convection suffices.
- **Which assumption fails first in the real world?** For the box, the still-air environment and the adiabatic tray: real bays have stray airflow and conduction paths. For the cold plate, the transition range.
- **What measurement could falsify the model?** Thermocouples on the box faces at 60 W in still air: a mean within ±2 K of 51.8 °C, and a top face slightly cooler than the sides if the case conducts poorly (its h is 5.93 against 4.57 W/(m² K)). A box much hotter than predicted would point to blocked radiation (a neighbour at the same temperature).
- **Does the result survive a change of method or input?** Property temperature: ±0.2 K. Correlation choice: Dittus–Boelter against Gnielinski changes the chamber channel's film ΔT from 50 to 55 K. The design numbers survive only with the stated margins.
- **Which decision changes if the uncertainty doubles?** At ±40 % on h_conv the box spans 48.3–56.5 °C; the decision "fan needed" would not change, but the fan sizing would have to be done by test.

## Check yourself

1. Define Nu, Re, Pr, Ra and St. Which one is a property of the fluid alone?
2. Why is the average Nusselt number of a laminar flat plate twice the local value at the trailing edge?
3. State the Reynolds–Colburn analogy. For which plate boundary condition is it exact, and why does it fail for uniform heat flux?
4. At what film temperature should properties be evaluated for a 50 °C wall in 20 °C air, and why does it barely matter for the box?
5. Why does h on the upper face of a large hot plate in turbulent natural convection not depend on the plate's size?
6. A designer uses Gnielinski at Re = 1500. What is wrong, and what Nu would you use instead for a uniformly heated tube?
7. Why does a black box in still air run more than 20 K cooler than a bare aluminium one?
8. The cold plate needs 8.8 g/s for a 40 °C wall. Why would you not design to 8.8 g/s?
