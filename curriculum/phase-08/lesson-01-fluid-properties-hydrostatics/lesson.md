# P08.01  Fluid Properties & Hydrostatics

> A fluid at rest carries no shear, so its pressure must hold up the weight above it: dp/dh = −ρg. Add the gas law and a temperature profile, and the whole standard atmosphere follows.

| | |
|---|---|
| **Learning outcomes** | 1. Define a fluid, the continuum, density, pressure and dynamic and kinematic viscosity, and compute the viscosity of air from Sutherland's law. 2. Derive the hydrostatic equation dp/dh = −ρg from a force balance and integrate it by hand, with the ideal-gas law, for an isothermal and a constant-lapse-rate layer. 3. Implement and test the International Standard Atmosphere to 47 km by two independent routes, the layer-by-layer closed forms and a fourth-order Runge–Kutta integration of the hydrostatic equation, and check both against the published anchor values. 4. Apply hydrostatics to manometers, an accelerating propellant tank and buoyancy, and verify Archimedes' principle by integrating surface pressure over a body. 5. Decide when a constant-density, isothermal or geometric-altitude simplification is good enough, and quantify density altitude on a non-standard day. |
| **Prerequisites** | `P07.01` (Systems, Properties & State: the ideal-gas law p = ρRT with R the specific gas constant, state variables), `P00.03` (the vectorised ISA table you built there, now derived), `P00.05` (units), `P03.05` (RK4 and observed order), `P03.04` (quadrature), `P04.01` (Newton's laws, free-body diagrams), `P04.07` (kinetic theory, for why gas viscosity rises with temperature). Skill check: you can integrate dp/p = −k dh and solve a one-variable root-finding problem with `scipy.optimize.brentq`. |
| **Engineering question** | What are the temperature, pressure, density and viscosity of the air the light aeroplane cruises in at 3000 m, and of the air a sounding rocket climbs through up to 47 km, computed from the hydrostatic equation and checked against the published ISA? The answer at 3000 m: 268.65 K, 70 108.5 Pa, 0.9091 kg/m³ (74.2 % of sea-level density), μ = 1.694 × 10⁻⁵ Pa·s and ν = 1.863 × 10⁻⁵ m²/s; the closed forms and the RK4 integration agree to 3 × 10⁻⁹ and both match the four published pressure anchors to within 4 × 10⁻⁶. |
| **Assumptions and limits** | Continuum (the molecular mean free path, tens of nanometres at sea level and growing as the density falls, ≪ every length of interest); fluid at rest (no shear stress); dry air as an ideal gas with constant R = 287.05287 J/(kg·K); the ISA's prescribed piecewise-linear temperature in geopotential altitude (a model day, not a real one); g₀ constant in geopotential altitude; Sutherland's law for air at ordinary temperatures. Fails above 86 km (dissociation, varying molar mass), in a real atmosphere with humidity and weather, and in any flowing fluid where shear matters. |
| **Analytical work** | The slab force balance and the hydrostatic equation; the closed-form pressure in gradient and isothermal layers and the a → 0 limit; the scale height; the constant-density approximation as the first Taylor term; geopotential altitude; manometer and accelerating-tank pressures; Archimedes' principle from the divergence theorem; density altitude. |
| **Computational deliverable** | `starter/main.py`: `layer_base_pressures`, `isa_temperature`, `isa_pressure`, `isa_density`, `integrate_hydrostatic`, `incompressible_pressure`, `scale_height`, `sutherland_viscosity`, `kinematic_viscosity`, `manometer_delta_p`, `tank_bottom_pressure`, `pressure_force_on_sphere`, `displaced_weight` (constants, `_layer_index`, `_base_temperatures`, `geopotential_altitude`, `observed_order` and `main` are given). Tests prove: the temperature profile and lapse rates; the published pressure anchors at 11, 20, 32 and 47 km to 10⁻⁵ and ρ(0) = 1.225, ρ(11 km) = 0.36392 kg/m³; continuity at the layer boundaries; the closed forms satisfy dp/dh = −ρg₀ and p = ρRT; the isothermal layer is exponential with the scale height; RK4 matches the closed forms to 10⁻⁸ and converges at order 4.0 ± 0.3; the constant-density error grows as h²; Sutherland's sea-level value and its √T limit; manometer, tank and buoyancy limits; Archimedes' principle by surface integration in a liquid and in the compressible ISA, converging at order 2. |
| **Practical activity** | **computational experiment**: a step-size study of the RK4 atmosphere with layer boundaries on and off the grid; the altitude at which the constant-density and isothermal approximations fail; viscosity and kinematic viscosity with altitude; the geopotential–geometric correction; hot and cold days and density altitude. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Lesson `P08.19` reads a water manometer at home.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `isa_atmosphere.png`, the derivations of the two layer formulas and of the scale height, the five experiment tables, and a half-page note giving the 3000 m properties with the two routes' agreement and the simplification that fails first. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: dp/dh = −ρg is derived from a force balance and integrated for both layer types; all public tests pass; the observed RK4 order is reported with and without aligned layer boundaries and the difference explained; the note states the 3000 m properties and when the constant-density model fails. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapters "Fundamental Properties of Fluids" (What is a Fluid?; Concept of a Continuum; Fluid Pressure; Equation of State; Viscosity of Fluids), "Fluid Statics & the Hydrostatic Equation" (Aerostatic/Hydrostatic Pressure Field; Solutions to the Aero-Hydrostatic Equation; Archimedes' Principle) and "Atmospheric Properties" (Geopotential Altitude vs Geometric Altitude; Equations of the ISA in the Troposphere; Pressure Altitude and Density Altitude), https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/fluid-statics/ and https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/international-standard-atmosphere-isa/ (required-free; link only, CC BY-NC-ND). `@ntrs-us-standard-atmosphere-1976` Part 1, 1.2.1 Adopted constants, 1.2.3 Gravity and geopotential altitude, 1.3.1 Pressure, 1.3.11 Dynamic viscosity, 1.4.1 Sea-level values, and Part 4 Table I (required-free: the defining constants, the anchor values the tests use and Sutherland's constants). `@libretexts-bar-meir-fluid-mechanics` Ch. 4 "Fluids Statics": 4.2 The Hydrostatic Equation; 4.3 Pressure and Density in a Gravitational Field; 4.6 Buoyancy and Stability (supplementary-free). `@openstax-university-physics-vol1` Ch. 14: 14.2 Measuring Pressure; 14.4 Archimedes' Principle and Buoyancy (supplementary-free: manometers and buoyancy at first-year level). `@libretexts-marchman-aerodynamics-aircraft-performance` Ch. 1: 1.4 The Standard Atmosphere; 1.7 Altitude Measurement (supplementary-free). `@ntrs-naca-tr-1235-standard-atmosphere` (supplementary-free: the historical table format). |
| **Extensions** | ★ Extend the model to 86 km with the two mesosphere layers of USSA-1976 1.2.5 and check the 71 km anchor. ★ Add humidity: replace R by the virtual-temperature form and compute the density altitude of a humid tropical day. ★ Integrate the hydrostatic equation in geometric altitude with g(z) = g₀(r₀/(r₀ + z))² and show it reproduces the geopotential result to round-off. |
| **Time** | unmeasured |

## 1. The engineering problem

Every aerodynamic number the rest of this course computes needs the air's state first. Lift and drag scale with density; the Reynolds number needs density and viscosity; the speed of sound needs temperature; an altimeter converts pressure to height; an engine's power falls with density. The light aeroplane of the earlier phases cruises at 3000 m at 60 m/s; a sounding rocket in `P08.03` and `P08.16` climbs through the stratosphere. Before either can be analysed we need p, T, ρ and μ as functions of altitude, and we need to trust them.

The International Standard Atmosphere (ISA) is that table. It is not a measurement of any real day; it is an agreed model, defined by a sea-level state, a temperature profile and one law of physics: a fluid at rest supports its own weight by a pressure gradient. Built from that, the ISA is the reference for altimeter calibration, aircraft performance charts and the airspeed definitions of `P08.04`. In `P00.03` you tabulated it from given formulas; here you derive those formulas, integrate the governing equation numerically as an independent check, and find where the simplifications people use (constant density, an isothermal atmosphere, geometric altitude) stop being good enough.

Success criteria, written first: the closed forms reproduce the published ISA anchor pressures (22 632.06 Pa at 11 km, 5474.889 Pa at 20 km, 868.019 Pa at 32 km, 110.906 Pa at 47 km) to 10⁻⁵ and ρ(11 km) = 0.36392 kg/m³; an RK4 integration that never sees the closed forms agrees with them to 10⁻⁸ at a 100 m step and converges at fourth order; Sutherland's law gives μ₀ = 1.7894 × 10⁻⁵ Pa·s; and surface integration of pressure over a balloon reproduces the weight of the air it displaces.

## 2. Understand and derive

### 2.1 What a fluid is, and the properties we need

A **fluid** is a substance that deforms continuously under a shear stress, however small; a solid resists it with a finite deformation. A fluid can therefore be at rest only if it carries no shear stress at all. We treat it as a **continuum**: properties such as density are defined at a point as averages over a volume small compared with the flow but large compared with the molecular spacing. That holds when the mean free path (tens of nanometres in sea-level air) is far below every length in the problem; it fails in the upper atmosphere and in micro-devices (Leishman, "Fundamental Properties of Fluids", Concept of a Continuum).

- **Density** ρ (kg/m³): mass per unit volume.
- **Pressure** p (Pa = N/m²): the normal force per unit area that the fluid exerts on any surface, the same in every direction at a point of a fluid at rest (Pascal's law).
- **Equation of state**: for air, the ideal-gas law **p = ρRT** with R = 287.05287 J/(kg·K), the universal gas constant divided by the molar mass of dry air (`P07.01`).
- **Dynamic viscosity** μ (Pa·s): the ratio of shear stress to shear rate in a Newtonian fluid, τ = μ du/dy (`P08.06` builds on it). **Kinematic viscosity** ν = μ/ρ (m²/s) is the diffusivity of momentum.

In a gas, viscosity comes from molecules carrying momentum between layers moving at different speeds (`P04.07`); faster molecules carry it faster, so μ rises with temperature. **Sutherland's law** fits air well at ordinary temperatures (USSA-1976 1.3.11):

**μ = C₁ T^{3/2}/(T + S),  C₁ = 1.458 × 10⁻⁶ kg/(m·s·K^{1/2}),  S = 110.4 K.**   (1)

At 288.15 K it gives 1.7894 × 10⁻⁵ Pa·s. For T ≫ S it tends to C₁√T, the kinetic-theory result for hard spheres. A liquid behaves the other way: its viscosity comes from intermolecular attraction and falls with temperature.

### 2.2 The hydrostatic equation

Take a thin horizontal slab of fluid at rest, of area A and thickness dh, with h measured upward.

![A thin horizontal slab of fluid at rest, of area A, between heights h and h + dh: pressure p + dp pushes down on its top face, pressure p pushes up on its bottom face, and its weight ρ g A dh acts downward. Pressure on the side faces cancels horizontally.](figures/hydrostatic-slab.svg)
*Figure 1. Free-body diagram of a fluid slab at rest. Only the pressures on the top and bottom faces and the weight act vertically, which gives dp/dh = −ρg.*

Vertical equilibrium: pA − (p + dp)A − ρgA dh = 0, so

**dp/dh = −ρg.**   (2)

The no-shear property of a fluid at rest is what makes the side faces irrelevant (Leishman, "Fluid Statics & the Hydrostatic Equation", Aerostatic/Hydrostatic Pressure Field; Bar-Meir 4.2). Equation (2) is the momentum equation with zero velocity and the hydrostatic stress model: the stress is −p times the identity.

**A liquid.** With ρ constant, p = p_s + ρg·depth. That is all a manometer, a fuel-tank gauge or a dam needs. In a U-tube with liquid of density ρ_m and legs filled above it by a fluid of density ρ_f, following the pressure down one leg and up the other gives **Δp = (ρ_m − ρ_f) g Δh**. Inclining the tube at θ to the horizontal spreads the same vertical Δh over Δh/sin θ of tube, improving the reading resolution by 1/sin θ (5.76 at 10°) without changing the pressure. A 25 mm water column under air is 244.43 Pa.

**An accelerating tank.** In the frame of a propellant tank accelerating upward at a, the liquid feels an effective gravity g + a, and the outlet pressure is p_ullage + ρ(g + a)h. For liquid oxygen (1141 kg/m³, a typical value) in a 2 m column under 3.0 bar, the outlet sees 3.224 bar on the pad and 3.671 bar at a total of 3g. In a coast (a = −g) the head vanishes and the liquid is no longer held over the outlet, which is why upper stages fire small ullage motors before relighting.

### 2.3 The atmosphere: two layer types

For a gas, ρ = p/(RT), and (2) becomes

**dp/p = −g₀ dh/(R T(h)).**   (3)

The ISA prescribes T(h) as straight segments in **geopotential altitude** h (section 2.4): −6.5 K/km from 288.15 K to 11 km, isothermal at 216.65 K to 20 km, +1.0 K/km to 32 km and +2.8 K/km to 47 km (USSA-1976 1.2.5).

**Gradient layer** (T = T_b + a(h − h_b), a ≠ 0). With dh = dT/a, (3) integrates to

**p = p_b (T/T_b)^{−g₀/(aR)}.**   (4)

In the troposphere −g₀/(aR) = 5.2559, so p ∝ T^{5.256} and ρ ∝ T^{4.256}.

**Isothermal layer** (a = 0). (3) integrates directly to

**p = p_b exp(−(h − h_b)/H),  H = RT_b/g₀,**   (5)

with the **scale height** H: the height over which pressure falls by e. It is 8435 m at 288.15 K and 6342 m at 216.65 K. Equation (5) is the a → 0 limit of (4) (Problem 2). Stepping (4) and (5) up from p₀ = 101 325 Pa gives each layer's base pressure; that is `layer_base_pressures`, and `isa_pressure` then evaluates the right formula in each layer (Leishman, "Fluid Statics & the Hydrostatic Equation", Solutions to the Aero-Hydrostatic Equation).

**The constant-density shortcut.** Expanding (4) about sea level, p = p₀ − ρ₀g₀h + O(h²): the liquid formula is the first-order term. Its error grows as h²; it is 1 % low at 1249 m, 5 % low at 2609 m, and it reaches zero pressure at p₀/(ρ₀g₀) = RT₀/g₀ = 8434.5 m, exactly one scale height.

### 2.4 Geopotential altitude

Gravity weakens with geometric height z as g₀(r₀/(r₀ + z))². The ISA absorbs that by defining **geopotential altitude** h, the height that has the same potential energy per unit mass at constant g₀: g₀ dh = g dz, which integrates to **h = r₀z/(r₀ + z)** with r₀ = 6 356 766 m (USSA-1976 1.2.3; Leishman, "Atmospheric Properties", Geopotential Altitude vs Geometric Altitude). In h, equation (2) holds with g₀ exactly. The difference is 1.4 m at 3 km, 19 m at 11 km and 345 m at 47 km.

### 2.5 Buoyancy

Integrate the pressure over the closed surface S of a body: F = −∮p n dA = −∫∇p dV (divergence theorem) = ∫ρg dV upward by (2). The buoyant force is the weight of the displaced fluid, **Archimedes' principle**, and it holds in a compressible atmosphere too, with ρ the local density averaged over the body (Leishman, "Fluid Statics & the Hydrostatic Equation", Archimedes' Principle). For a small body ρ at the centre times gV is accurate to second order in its size. The code computes the left-hand side by integrating surface pressure over a sphere and the right-hand side by integrating density over its volume: two routes that share only the atmosphere.

**A worked example: a balloon at the cruise altitude.** A sphere of radius 2 m centred at 3000 m: the surface pressure integral (400 bands) gives 298.762 N, the volume integral of ρg gives 298.759 N, and ρgV at the centre 298.759 N. Filled with helium (R = 2077.1 J/(kg·K)) at ambient pressure and temperature, its net lift is 26.25 kg; at sea level the same volume would lift 35.4 kg.

### 2.6 Where the model stops

- **A real day.** Temperatures depart from the ISA by ±20 K or more, and humid air is slightly less dense than dry air at the same p and T (water vapour is lighter than air). Part 4 shows what a 15 K warm day does.
- **High altitude.** Above 86 km the gas dissociates and diffuses; the molar mass is no longer constant and the continuum assumption itself weakens.
- **Motion.** The moment the air moves with shear, (2) is replaced by the momentum equation with viscous stresses (`P08.06`); the hydrostatic pressure distribution survives only across thin layers such as a boundary layer (`P08.10`).

## 3. Build the computational model

Open `starter/main.py`. The constants, `_layer_index`, `_base_temperatures`, `geopotential_altitude`, `observed_order` and `main` are given.

### Step 1: the temperature profile and the base pressures

`isa_temperature(h)` finds each altitude's layer with `_layer_index` and returns T_b + a(h − h_b). `layer_base_pressures()` steps (4) or (5) up the four layers from P0; keep it separate so you can print and check it (22 632.04, 5474.877, 868.016, 110.906 Pa).

### Step 2: pressure, density, viscosity

`isa_pressure(h)` evaluates (4) or (5) in the right layer, vectorised: compute both forms with `np.where`, guarding the division by a = 0. `isa_density` is p/(RT). `sutherland_viscosity` is (1), refusing T ≤ 0; `kinematic_viscosity` divides by the ISA density.

### Step 3: the independent route

`integrate_hydrostatic(h_max, dh, p_start, T_of_h)` integrates (3) with classical RK4 (`P03.05`) on a fixed step, using only T(h):

```python
def f(hh, pp):
    return -G0 * pp / (R_AIR * float(T_of_h(hh)))
for i in range(n):
    k1 = f(h[i], p[i]); k2 = f(h[i] + dh/2, p[i] + dh/2*k1)
    k3 = f(h[i] + dh/2, p[i] + dh/2*k2); k4 = f(h[i] + dh, p[i] + dh*k3)
    p[i+1] = p[i] + dh/6*(k1 + 2*k2 + 2*k3 + k4)
```

The `T_of_h` argument lets the same integrator compute a non-standard day. T(h) has kinks at 11, 20 and 32 km; RK4's error expansion needs a smooth right-hand side within each step, so choose dh to put the kinks on step boundaries.

### Step 4: manometers, tanks and buoyancy

`manometer_delta_p`, `tank_bottom_pressure` and `incompressible_pressure` are one line each from 2.2–2.3. `pressure_force_on_sphere(p_of_z, z_c, R, n)` sums −p cos θ · 2πR² sin θ dθ over n polar bands at their midpoints; `displaced_weight` sums ρ(z)g over n horizontal discs.

Running `main()` prints

```
ISA anchors (published data) against the closed forms and RK4 (dh = 100 m):
   h (m)    T (K)    published p (Pa)  closed form (Pa)   RK4 (Pa)       rel. diff closed/published
   11000   216.65       22632.060       22632.040       22632.040   -8.8e-07
   20000   216.65        5474.889        5474.877        5474.877   -2.1e-06
   32000   228.65         868.019         868.016         868.016   -3.7e-06
   47000   270.65         110.906         110.906         110.906   -2.0e-06
rho(0) = 1.22500 kg/m^3, rho(11 km) = 0.36392 kg/m^3 (published 1.225, 0.36392)
mu(T0) = 1.78938e-05 Pa s (published 1.7894e-5)

light aeroplane at 3000 m: T = 268.65 K, p = 70108.5 Pa, rho = 0.9091 kg/m^3, sigma = rho/rho0 = 0.7421, mu = 1.6937e-05 Pa s, nu = 1.8630e-05 m^2/s
  incompressible estimate p0 - rho0 g h = 65285.6 Pa (-6.88%)

RK4 step study, error at 47 km against the closed form:
  dh = 1000.0 m: relative error 3.018e-05
  dh =  500.0 m: relative error 1.773e-06
  dh =  250.0 m: relative error 1.074e-07
  dh =  125.0 m: relative error 6.611e-09
  observed orders [4.09 4.04 4.02]

scale height at T0: 8435 m; at 216.65 K: 6342 m

U-tube, 25.0 mm of water under air: 244.43 Pa; inclined at 10 deg, the same pressure gives 144.0 mm along the tube
LOX tank (1141 kg/m^3 typical), 2 m column, 3.0 bar ullage, 3 g total: outlet 3.671 bar (on the pad 3.224 bar)
balloon r = 2.0 m at 3000 m: surface-pressure force 298.762 N, displaced weight 298.759 N, rho g V = 298.759 N
```

and writes `isa_atmosphere.png`: temperature, pressure by both routes with the published anchors, and μ and ν against altitude.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-01-fluid-properties-hydrostatics`. The 20 tests:

1. **Profile**: T at 0, 11, 15, 20, 32, 47 km and the four lapse rates by finite difference.
2. **Published anchors**: p at 11, 20, 32, 47 km to 10⁻⁵; ρ(0) and ρ(11 km) to 10⁻⁴.
3. **Consistency**: p and T continuous at each layer boundary; the closed form satisfies dp/dh = −ρg₀ (central differences, 10⁻⁷) at eight heights in all four layers; p = ρRT; the isothermal layer equals p₁₁ exp(−Δh/H).
4. **Interface**: vectorised, monotonic, and a `ValueError` above 47 km.
5. **RK4**: matches the closed form to 10⁻⁸ at dh = 100 m and the anchors to 10⁻⁵; observed order 3.7–4.3 for dh = 1000, 500, 250 m; with a constant temperature it is the exponential.
6. **Limits**: the constant-density error quadruples when h doubles; geopotential altitude is below geometric and equal near the ground.
7. **Viscosity**: Sutherland's sea-level value to 10⁻⁴, monotonic in T, the √T limit, rejection of T ≤ 0; ν rises with altitude.
8. **Statics**: no manometer reading between identical fluids; the inclined-tube identity; the tank's free-fall and linear-head limits.
9. **Archimedes**: surface integration equals ρgV in a liquid at two depths, a uniform pressure gives no force, and in the ISA the surface and volume integrals agree to 10⁻⁴ with midpoint order above 1.8.
10. `main` writes the figure.

### Compare

The closed forms and the RK4 integration are independent routes to the same pressure: one is the integral of (3) done by hand, the other a numerical integration that knows only T(h). They agree to 2.7 × 10⁻⁹ everywhere at dh = 100 m, and the step study shows the error falling as dh⁴ (orders 4.09, 4.04, 4.02). Both differ from the published anchors by 0.9–3.7 × 10⁻⁶; that difference is the same for both routes, so it is not a numerical error but a constant: the published values were computed with slightly different rounding of the defining constants and base pressures. 

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Step size and kinks.** Run `integrate_hydrostatic` to 47 km with dh = 1000, 500, 250, 125 m (layer boundaries on the grid) and with dh = 47 000/60, /120, /240, /480 m (boundaries between grid points). Predict both orders. *Finding:* aligned: errors 3.0 × 10⁻⁵, 1.8 × 10⁻⁶, 1.07 × 10⁻⁷, 6.6 × 10⁻⁹, orders 4.09, 4.04, 4.02. Unaligned: errors 2.5 × 10⁻⁵, 1.0 × 10⁻⁵, 2.3 × 10⁻⁶, 1.7 × 10⁻⁷, orders 1.28, 2.17, 3.75: erratic, because a step that straddles a kink in T(h) loses RK4's error cancellation and how much it loses depends on where the kink falls in the step. Put the kinks on the grid.
2. **When does constant density fail?** Compare p₀ − ρ₀g₀h and the isothermal formula p₀e^{−h/H₀} with the ISA at 1, 3, 5, 8, 11 km. Predict the signs. *Finding:* constant density −0.63, −6.88, −23.6, −85.3, −236 %; isothermal (T₀) +0.14, +1.27, +3.7, +10.2, +21.5 %. Constant density is good to 1 % only below 1249 m (5 % below 2609 m): fine for a cabin or a fuel tank, useless for flight. The isothermal model errs the other way because real air is colder and denser aloft.
3. **Viscosity with altitude.** Tabulate μ and ν at 0, 3, 11, 20, 32, 47 km. Predict which changes more. *Finding:* μ = 1.789, 1.694, 1.422, 1.422, 1.487, 1.704 × 10⁻⁵ Pa·s (it follows T only); ν = 1.461 × 10⁻⁵ m²/s at sea level, 1.275, 2.674, 11.06, 77.0 and 817 times that. A power law μ ∝ T^n fitted over 200–300 K gives n = 0.81. Because ν grows so fast, a rocket at 47 km at the same speed and size as at sea level has a Reynolds number 800 times smaller (`P08.05`).
4. **Geometric versus geopotential.** At geometric 3, 11 and 47 km compute h and the pressure error of using z in place of h. *Finding:* 1.4 m and 0.018 %; 19.0 m and 0.30 %; 345 m and 4.27 %. Negligible for the aeroplane, not for the rocket's pressure altimeter.
5. **A non-standard day.** At a pressure altitude of 3000 m (p = 70 108.5 Pa), with the whole atmosphere offset by ΔT = −15, +15, +25 K, compute the density altitude and, with `integrate_hydrostatic(..., T_of_h=lambda h: isa_temperature(h) + ΔT)`, the true height. *Finding:* density altitude 2438, 3524, 3855 m; true height 2844, 3156, 3260 m. On a 15 K hot day the aeroplane performs as if 524 m higher and is in fact 156 m higher than its altimeter reads. Cold days are the dangerous ones for terrain clearance: the altimeter reads high.

Raw data format: CSV tables `dh_m, aligned, rel_err_47km` (1); `h_m, p_isa_Pa, err_const_rho, err_isothermal` (2); `h_m, T_K, mu_Pa_s, nu_m2_s` (3); `z_m, h_m, p_err` (4); `dT_K, density_alt_m, true_height_m` (5); plus `isa_atmosphere.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| p at 11, 20, 32, 47 km | published anchors | closed form and RK4 | 0.9–3.7 × 10⁻⁶ | rounding of the published constants |
| closed form vs RK4 | identical in exact arithmetic | 2.7 × 10⁻⁹ at dh = 100 m | 2.7 × 10⁻⁹ | RK4 truncation, order 4 |
| ρ(0), ρ(11 km) | 1.225, 0.36392 (published) | 1.22500, 0.36392 | < 10⁻⁵ | — |
| μ(288.15 K) | 1.7894 × 10⁻⁵ (published) | 1.78938 × 10⁻⁵ | 10⁻⁵ | rounding |
| buoyancy at 3000 m, r = 2 m | ρgV = 298.759 N | surface integral 298.762 N | 1 × 10⁻⁵ | midpoint rule, order 2 |
| scale height | RT₀/g₀ | 8434.5 m = p₀/(ρ₀g₀) | 10⁻⁸ | identity |

Every result has been reproduced by an independent route. The engineering answer: **at 3000 m the ISA air is at 268.65 K and 70 108.5 Pa, with density 0.9091 kg/m³ (σ = 0.742), viscosity 1.694 × 10⁻⁵ Pa·s and kinematic viscosity 1.863 × 10⁻⁵ m²/s**, and the numerical and closed-form routes agree to nine digits. The uncertainty that matters is not numerical: a real day differs from the ISA by tens of kelvin, which moves the density by several per cent (Experiment 5) and the rocket's pressure altitude by a kilometre at 30 km (Problem 6).

What to improve: when the actual temperature profile is known (a radiosonde, or the aircraft's own sensor), integrate (3) with it through `T_of_h` rather than using the ISA; that is a one-line change because the integrator takes T(h) as an argument. What not to improve: the RK4 step below 100 m.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (4), (5) and the scale height; `isa_atmosphere.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the 3000 m properties, the two routes' agreement, and the simplification that fails first for the aeroplane and for the rocket.

### Adversarial self-checks

- **What if the temperature goes to zero lapse everywhere?** The atmosphere becomes isothermal and (4) must reduce to (5): the test with constant temperature checks exactly that limit in the integrator.
- **What if the altitude is 50 km?** The model refuses: the code raises rather than extrapolating a layer whose temperature profile it does not know.
- **Which assumption fails first in the real world?** The ISA temperature profile: a 15 K warm day lowers the density at a pressure altitude of 3000 m by 5.3 % (ratio 268.65/283.65), more than any numerical or geopotential effect at that height.
- **What measurement could falsify the model?** A radiosonde profile of p and T: integrating (3) with its T(h) must reproduce its p(h) to the sensor's accuracy; a systematic drift would point to humidity (R changes) or to sensor lag, not to hydrostatics.
- **Does the result survive a change of method or step?** Closed form and RK4 agree to 3 × 10⁻⁹; halving the step reduces the RK4 error sixteenfold when the layer boundaries are on the grid.
- **Which decision changes if the uncertainty doubles?** A take-off-distance calculation at a hot, high airfield: take-off distance grows faster than 1/σ, so doubling the density uncertainty more than doubles the uncertainty of the distance, enough to move a marginal take-off from allowed to forbidden.

## Check yourself

1. What is a fluid, and why can a fluid at rest carry no shear stress?
2. Write the hydrostatic equation and say what the sign means.
3. Why does pressure fall exponentially in an isothermal layer, and what is the scale height at 216.65 K?
4. Why does the viscosity of air rise with temperature while that of water falls?
5. What is the difference between dynamic and kinematic viscosity, and which grows faster with altitude?
6. Why does an inclined manometer improve resolution, and does it change the pressure it measures?
7. Is Archimedes' principle still "ρgV" in a compressible atmosphere?
8. What is geopotential altitude, and why is the ISA defined in it?
