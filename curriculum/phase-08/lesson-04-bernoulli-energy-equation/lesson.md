# P08.04  Bernoulli & Energy Equation

> Along a streamline in steady, frictionless, incompressible flow, pressure and speed trade one for the other: p + ½ρV² + ρgz stays fixed. A pitot tube turns that trade into an airspeed, and every condition of the trade is a way for the reading to go wrong.

| | |
|---|---|
| **Learning outcomes** | 1. Derive Bernoulli's equation from the momentum balance along a streamline, naming each of its four conditions where it is used, and derive the pitot-static relation and the venturi formula. 2. Define TAS, EAS, CAS and IAS, and compute each from the others for the light aeroplane, including the effect of a static-port pressure error on airspeed and altitude. 3. Implement and test the airspeed chain with every conversion inverted to round-off, the compressible calibration law tending to Bernoulli as M → 0, and an independent numerical integration of the Euler equation through a wind-tunnel contraction converging to Bernoulli at second order. 4. Apply the steady-flow energy equation with head losses and shaft work to a pump and a wind tunnel. 5. Decide when Bernoulli's equation may be used and what replaces it when it may not. |
| **Prerequisites** | `P08.01` (ISA, hydrostatics, manometers), `P08.03` (control-volume momentum), `P08.02` (streamlines, the material derivative), `P04.02` (work and energy), `P07.01` (Systems, Properties & State: the first law for an open system as energy in = energy out plus storage), `P03.02` (root finding), `P03.04` (quadrature). Skill check: you can solve ½ρV² = Δp for V and invert a monotonic function with `brentq`. |
| **Engineering question** | The light aeroplane cruises at 3000 m (ISA) at a true airspeed of 60 m/s. What do its pitot and static ports sense, what does the airspeed indicator read, how big an error does a static port in slightly fast air cause, and how wrong is plain Bernoulli? The answers: static 70 108.5 Pa, pitot 71 758.6 Pa (impact pressure 1650.1 Pa); CAS 51.75 m/s (100.6 kt), EAS 51.69 m/s (100.5 kt), TAS 116.6 kt; a port at c_p = −0.02 makes the dial read 0.99 kt high and the altimeter 3.7 m high; Bernoulli with the true density overestimates TAS by 0.42 %. |
| **Assumptions and limits** | Steady, inviscid, incompressible flow along a streamline (Bernoulli); subsonic isentropic deceleration for the standard calibration law (used as given here, derived in `P08.17`); ISA day at the pressure altitude; the static port senses p + c_p q with c_p a constant property of its location; no instrument lag. The energy equation assumes steady flow between two sections with uniform or corrected profiles. Bernoulli fails across fans and propellers (work), in boundary layers and wakes (viscosity), above M ≈ 0.3 (compressibility), in rapidly changing flows (unsteadiness) and between streamlines of a rotational flow. |
| **Analytical work** | Euler's equation along a streamline and its integration; stagnation, static and dynamic pressure; the pitot-static and venturi relations; the low-Mach expansion q_c = q(1 + M²/4 + …); EAS, CAS and IAS; the static-port error in airspeed and altitude; the steady-flow energy equation with losses and shaft work, and Bernoulli as its lossless limit; the energy budget of an open-circuit wind tunnel. |
| **Computational deliverable** | `starter/main.py`: `dynamic_pressure`, `bernoulli_airspeed`, `cas_from_impact_pressure`, `impact_pressure_from_cas`, `airspeeds_from_tas`, `tas_from_cas`, `static_error_effects`, `venturi_flow_rate`, `euler_pressure_along_duct`, `bernoulli_pressure_along_duct`, `shaft_head`, `shaft_power`, `open_tunnel_fan_power` (`isa_troposphere`, `impact_pressure`, `contraction_area`, `observed_order` and `main` are given). Tests prove: the troposphere matches the published ISA; Bernoulli inverts the dynamic pressure; the compressible law tends to Bernoulli with the M²/4 term; at sea level CAS = EAS = TAS; CAS ↔ q_c and TAS ↔ CAS round trips close to 10⁻¹² at four altitudes; EAS keeps the dynamic pressure and CAS exceeds it by less than 0.3 %; the static-port error has the right signs and its small-error linearisation; the altimeter error agrees with the ISA's hydrostatic gradient; the Euler integration converges to Bernoulli at order 2.0 ± 0.2; total pressure is constant along the duct; the venturi formula inverts; the energy equation reduces to Bernoulli, adds losses linearly and gives negative head for a turbine; tunnel power scales with V³ and falls with diffuser area. |
| **Practical activity** | **computational experiment**: the airspeed chain across altitude, a static-port error sweep, the Mach number at which Bernoulli fails, a grid study of the Euler integration through a contraction, and the energy budget of an open-circuit tunnel against diffuser size and losses. |
| **Equipment** | None. A computer with the environment from `P00.01`. (Lesson `P08.19` measures an airspeed with a home-made pitot probe and a water manometer.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pitot_static.png`, the derivations of Bernoulli, the pitot relation and the M²/4 term, the five experiment tables, and a half-page note on the aeroplane's airspeeds and how a static-port error propagates. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: Bernoulli is derived with its four conditions named; all public tests pass; the IAS-to-TAS conversion for your own numbers is correct and inverted; the note states where Bernoulli stops being good enough and what replaces it. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapters "Energy Equation & Bernoulli's Equation" (Bernoulli's Equation from Euler's Equation; Other Forms of the Bernoulli Equation), "Applications of the Conservation Laws" (Flow in a Venturi; Total, Static, & Dynamic Pressures; Pitot Tube & Pitot-Static Tube; Flow Speed in a Wind Tunnel (Take 2)) and "Airspeed Definitions & Measurement" (Definitions of Airspeeds; Pitot-Static System), https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/determination-of-airspeed/ (required-free; link only). `@nasa-glenn-beginners-guide-aeronautics` "Bernoulli's Equation", "Pitot – Static Tube – Speedometer", "Dynamic Pressure" (required-free). `@libretexts-marchman-aerodynamics-aircraft-performance` Ch. 1: 1.8 Bernoulli's Equation; 1.9 Airspeed Measurement; 1.10 Bernoulli's Equation for Compressible Flow (supplementary-free: the pilot's-eye view of the airspeeds). `@mit-ocw-2-25-advanced-fluid-mechanics` "Inviscid Flow and Bernoulli" class notes: Equation of Motion in Streamline Coordinates; Unsteady Bernoulli Equation (supplementary-free). `@libretexts-bar-meir-fluid-mechanics` Ch. 7 "Energy Conservation": 7.3 Approximation of Energy Equation; 7.4.5 Energy Losses in Incompressible Flow (supplementary-free). `@openstax-university-physics-vol1` 14.6 Bernoulli's Equation (supplementary-free). |
| **Extensions** | ★ Add instrument lag: model the pitot line as a first-order system with a time constant and find the IAS error during a 2 kt/s acceleration. ★ Fly a "tower fly-by" in code: generate pitot and static pressures with a position error that varies with angle of attack, and recover c_p(α) from passes at known TAS. ★ Use the unsteady Bernoulli equation (with ∂φ/∂t) for a pitot tube in an oscillating stream and find the frequency at which the steady formula errs by 1 %. |
| **Time** | unmeasured |

## 1. The engineering problem

An aeroplane's airspeed indicator is a pressure gauge. A forward-facing **pitot tube** brings the oncoming air to rest and senses the resulting total pressure; a **static port**, a flush hole on the side of the fuselage, senses the pressure of the passing air. The instrument converts the difference into a speed through a calibration law and shows it on a dial. Everything a pilot does with airspeed (rotating for take-off, avoiding the stall, staying below the structural limit speed, navigating) rests on that one pressure difference and on what it means.

It means less than it seems. The difference depends on the air density, which falls with height, so the dial does not show the true airspeed; it is built to show the speed that would produce the same pressure difference at sea level. The static port senses the true static pressure only if it sits where the local flow is at free-stream pressure; a port in slightly fast air reads low, the airspeed reads high and the altimeter, fed by the same port, reads high too. And the conversion assumes Bernoulli's equation, whose four conditions are each violated somewhere on an aeroplane.

This lesson works the chain for the light aeroplane at 3000 m and 60 m/s true airspeed, checks every conversion by inverting it, verifies Bernoulli itself by integrating the Euler equation through a wind-tunnel contraction, and then extends it to the energy equation with losses and shaft work that sizes fans and pumps. Success criteria, written first: every conversion inverts to 10⁻¹²; the calibration law tends to Bernoulli with a relative difference of M²/4; the Euler integration converges to Bernoulli at second order; a static-port error moves the altimeter by exactly the height over which the ISA pressure changes by that error.

## 2. Understand and derive

### 2.1 Bernoulli's equation from the momentum balance

Take a small fluid element of cross-section δA and length δs along a streamline (coordinate s, unit tangent along the velocity V; Leishman, "Energy Equation & Bernoulli's Equation", Bernoulli's Equation from Euler's Equation; 2.25 "Equation of Motion in Streamline Coordinates"). The forces along s are the pressure difference across its ends and the component of its weight; there are **no shear forces (inviscid)**. Its acceleration along s is ∂V/∂t + V ∂V/∂s, and in **steady** flow the first term is zero. Dividing by the volume,

**∂p/∂s + ρV ∂V/∂s + ρg ∂z/∂s = 0.**   (1)

This is the Euler equation along a streamline. For **constant density** it is the derivative along s of one quantity, so

**p + ½ρV² + ρgz = constant along a streamline,**   (2)

provided nothing between the two points adds or removes energy (**no shaft work**). These are the four conditions: steady, inviscid, incompressible, no work, and the two points on the same streamline. In an irrotational flow the constant is the same on every streamline (`P08.08`).

**Total, static and dynamic pressure.** p is the **static pressure**, the pressure the moving fluid exerts on a surface moving with it; q = ½ρV² is the **dynamic pressure**; their sum p₀ = p + q is the **total (stagnation) pressure**, reached when the flow is brought to rest without loss.

### 2.2 The pitot-static tube

Follow the streamline that ends on the nose of a forward-facing tube. Far upstream: p, V. At the nose (stagnation point): p₀, 0. Equation (2) gives

**p₀ − p = ½ρV²,  V = √(2(p₀ − p)/ρ).**   (3)

The static holes, on the side of the tube or the fuselage where the flow has recovered to free-stream speed, sense p. In the lesson's tunnel a water manometer reading 45 mm gives 27.03 m/s (Problem 1).

**Compressibility.** In a gas the deceleration to rest compresses the air, and the pitot senses more than ½ρV². For subsonic isentropic deceleration (derived in `P08.17` and used here as given), the **impact pressure** is

**q_c = p₀ − p = p[(1 + ((γ − 1)/2)M²)^{γ/(γ−1)} − 1],  M = V/√(γRT).**   (4)

Expanding with x = 0.2M², (1 + x)^{3.5} − 1 = 3.5x + 4.375x² + …, and since 0.7pM² = ½ρV²,

**q_c = ½ρV²(1 + M²/4 + …).**   (5)

At the cruise point M = 0.1826, so q_c/q − 1 = 0.00836 (M²/4 = 0.00834). Applying (3) to q_c with the true density gives a speed 0.42 % too high. The error depends on M alone and reaches 1 % at M = 0.282 and 5 % at M = 0.628 (Problem 5): "Bernoulli is fine below M ≈ 0.3" is a 1 % statement.

### 2.3 The airspeeds

The instrument cannot know the density at altitude. It is built to a standard law: the **calibrated airspeed** is the speed that would produce the measured q_c at sea-level standard conditions (p₀ = 101 325 Pa, a₀ = 340.29 m/s), the inverse of (4) at sea level:

**CAS = a₀ √(5[(q_c/p₀ + 1)^{2/7} − 1]).**   (6)

The other speeds (Leishman, "Airspeed Definitions & Measurement", Definitions of Airspeeds; Marchman 1.9):

- **TAS**, the true airspeed relative to the air, what navigation needs.
- **EAS = TAS √(ρ/ρ₀)**: the sea-level speed with the same dynamic pressure, ½ρ₀EAS² = ½ρTAS². Lift, drag and the stall depend on q, so the stall happens at a fixed EAS at every height. At low Mach number CAS ≈ EAS; at altitude CAS exceeds EAS slightly because q_c/q grows with M.
- **IAS**: what the dial shows, CAS plus the **position error** of the ports and the instrument error.

**At the cruise point** (3000 m ISA, TAS 60 m/s = 116.6 kt): p = 70 108.5 Pa, ρ = 0.9091 kg/m³, q = 1636.42 Pa, q_c = 1650.11 Pa, so the pitot senses 71 758.6 Pa; EAS = 51.689 m/s (100.5 kt), CAS = 51.754 m/s (100.6 kt). The dial reads 16 kt less than the aeroplane is flying through the air. The inverse, TAS from CAS at a pressure altitude (`tas_from_cas`), finds q_c from (6), the Mach number from q_c/p by inverting (4), and TAS = M·a; it returns 60.000000000000 m/s.

### 2.4 The static-port error

A static port sits on the fuselage where the local pressure coefficient c_p = (p_local − p)/q is close to zero, but not exactly. If it senses p + c_p q:

- the airspeed indicator sees q_c − c_p q: with c_p < 0 (the port in slightly fast air) IAS reads high;
- the altimeter, which converts pressure to height with the ISA (dp/dh = −ρg, `P08.01`), reads high by **Δh = −c_p q/(ρg)**.

At the cruise point, c_p = −0.02 gives Δp = −32.73 Pa, IAS 0.99 kt high and the altimeter 3.67 m high. For small errors, ΔCAS ≈ −c_p q/(ρ₀ CAS): the speed error is about c_p/2 relative, 1 % for c_p = −0.02 whatever the speed (Experiment 2). Position errors are calibrated in flight test for every aircraft type, and the stall speeds in the flight manual are published as IAS for that reason.

### 2.5 A venturi and a contraction

In a duct, continuity gives V₁A₁ = V₂A₂ (incompressible, `P08.03`), and (2) between the two sections gives the **venturi** relation

**Q = A₂ √(2Δp/(ρ(1 − (A₂/A₁)²))),**   (7)

with a discharge coefficient C_d ≤ 1 for real losses (Leishman, "Applications of the Conservation Laws", Flow in a Venturi; Flow Speed in a Wind Tunnel (Take 2)). A wind tunnel sets its speed exactly this way, from the pressure drop across the contraction measured by wall rings: Δp = q_test(1 − 1/CR²) with CR the contraction ratio. The independent check of Bernoulli is to integrate (1) itself numerically along the contraction's centreline, dp/dx = −ρV dV/dx with V = Q/A(x), and compare.

### 2.6 The energy equation with losses and work

Bernoulli is a mechanical energy balance with nothing lost. For steady flow between two sections of a pipe or duct, the first law for an open system (`P07.01`) per unit weight of fluid, with no heat transfer, is

**p₁/(ρg) + α₁V₁²/(2g) + z₁ + h_pump = p₂/(ρg) + α₂V₂²/(2g) + z₂ + h_turbine + h_L**   (8)

(Bar-Meir 7.3, 7.4.5). Each term is a **head** in metres; α is the kinetic-energy correction factor (1 for a uniform profile, 2 for laminar pipe flow); h_L, the **head loss**, is mechanical energy turned into heat by friction, the drop in total head when no machine acts. The shaft power is ρgQh/η. With no machine and no loss, (8) is (2).

**A worked example: an open-circuit wind tunnel.** Air at 1.204 kg/m³ is drawn from still room air, accelerated to 30 m/s in a 0.5 m² test section and discharged through a diffuser exit of 1.0 m². The fan must supply the kinetic energy thrown away at the exit, V_exit²/(2g) with V_exit = 15 m/s, plus losses written as 0.6 test-section dynamic heads. With a 70 % fan, (8) gives a shaft power of **9869 W**. With no losses and no diffuser the fan would supply exactly the kinetic-energy flux ½ρQV², 8127 W at 100 % efficiency.

### 2.7 Where Bernoulli stops

- **Work**: across a propeller or fan the total pressure jumps; use (2) on each side and a momentum or energy balance across it (`P08.03`, Problem 3).
- **Viscosity**: inside a boundary layer total pressure falls towards the wall (`P08.10`); in a separated wake it is lost to turbulence (`P08.16`).
- **Compressibility**: above M ≈ 0.3 the 1 % error of (5) grows quickly; use (4) (`P08.17`).
- **Unsteadiness**: in a gust with reduced frequency ωc/U of order 1 the ∂V/∂t term of the momentum balance is not negligible (the unsteady Bernoulli equation, 2.25 "Unsteady Bernoulli Equation").
- **Different streamlines**: in a rotational flow the constant in (2) changes from streamline to streamline.

## 3. Build the computational model

Open `starter/main.py`. `isa_troposphere`, `impact_pressure` (equation (4)), `contraction_area` and `observed_order` are given, with `main`.

### Step 1: Bernoulli and the calibration law

`dynamic_pressure` and `bernoulli_airspeed` are (3), refusing a negative difference. `cas_from_impact_pressure` is (6); `impact_pressure_from_cas` is (4) evaluated at sea-level standard conditions.

### Step 2: the airspeed chain

`airspeeds_from_tas(tas, h, cp_static, instrument_error)` gets T, p, ρ from `isa_troposphere`, computes q and q_c, applies the port error to both the sensed static pressure and the sensed impact pressure (q_c − c_p q), and returns TAS, EAS, CAS, IAS and the pressures. `tas_from_cas` inverts: q_c from CAS, M from q_c/p,

```python
M = np.sqrt(2/(GAMMA - 1) * ((qc/p + 1)**((GAMMA - 1)/GAMMA) - 1))
```

then TAS = M a. `static_error_effects` returns Δp = c_p q and Δh = −Δp/(ρg).

### Step 3: ducts and the energy equation

`venturi_flow_rate` is (7). `bernoulli_pressure_along_duct` is (2) with continuity; `euler_pressure_along_duct` integrates (1) by the trapezoid rule, with dV/dx from `np.gradient`. `shaft_head`, `shaft_power` and `open_tunnel_fan_power` are (8) and its use in 2.6.

Running `main()` prints

```
light aeroplane, 3000 m ISA, TAS 60.0 m/s (116.6 kt), M = 0.1826
  static p = 70108.5 Pa, rho = 0.9091 kg/m^3; q = 1636.42 Pa; impact pressure q_c = 1650.11 Pa (q_c/q - 1 = 0.00836, M^2/4 = 0.00834)
  EAS 51.689 m/s (100.5 kt)   CAS 51.754 m/s (100.6 kt)   IAS (perfect ports) 51.754 m/s
  inversion: tas_from_cas(CAS) = 60.000000000000 m/s
  Bernoulli with the true density on q_c: 60.250 m/s (+0.417% vs TAS)

static-port error (port in air of pressure coefficient cp):
    cp      dp (Pa)   altimeter error (m)   IAS (m/s)   IAS - CAS (kt)
   -0.04     -65.46         +7.34          52.765      +1.96
   -0.02     -32.73         +3.67          52.262      +0.99
   +0.00      +0.00         -0.00          51.754      +0.00
   +0.02     +32.73         -3.67          51.242      -1.00

Bernoulli vs the compressible law: (Bernoulli V - TAS)/TAS from q_c at 3000 m
   TAS    30 m/s  M = 0.091  error +0.1042%
   TAS    60 m/s  M = 0.183  error +0.4173%
   TAS   100 m/s  M = 0.304  error +1.1618%
   TAS   150 m/s  M = 0.457  error +2.6252%
   TAS   200 m/s  M = 0.609  error +4.6942%

contraction 4.0 -> 0.5 m^2 over 2.0 m, 30 m/s in the test section, air 1.204 kg/m^3:
   n =   25: max |Euler - Bernoulli| = 8.899e+00 Pa (pressure drop 533.33 Pa)
   n =   50: max |Euler - Bernoulli| = 2.292e+00 Pa (pressure drop 533.33 Pa)
   n =  100: max |Euler - Bernoulli| = 5.790e-01 Pa (pressure drop 533.33 Pa)
   n =  200: max |Euler - Bernoulli| = 1.454e-01 Pa (pressure drop 533.33 Pa)
   observed orders [1.96 1.98 1.99]
   venturi formula inverted from that drop: Q = 15.000000 m^3/s (true 15.0)

open-circuit tunnel, 30 m/s in 0.5 m^2, diffuser exit 1.0 m^2, losses 0.6 q_test, fan 70 %: shaft power 9869 W
```

and writes `pitot_static.png`: EAS, CAS and TAS against altitude at a constant 60 m/s TAS, and the contraction's pressure by Bernoulli and by the Euler integration.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-04-bernoulli-energy-equation`. The 19 tests:

1. **Atmosphere**: ρ(0) = 1.225, p(11 km) = 22 632.06 Pa (published), T(11 km) = 216.65 K; refusal above 11 km.
2. **Bernoulli**: (3) inverts q; four times the pressure gives twice the speed; negative differences refused.
3. **Compressibility**: q_c/q − 1 = M²/4 within 2 % at M = 0.01, 0.05, 0.1; zero at rest; refusal at M > 1.
4. **The chain**: CAS = EAS = TAS at sea level; CAS ↔ q_c to 10⁻¹³; TAS → CAS → TAS to 10⁻¹² at 0, 1500, 3000, 9000 m; EAS keeps q; EAS < CAS < TAS above sea level with CAS − EAS below 0.3 % and growing with height.
5. **Static port**: c_p < 0 raises IAS and lowers the sensed pressure; the altimeter reads high; ΔIAS = 0.02q/(ρ₀CAS) within 2 %; the instrument error adds; and the ISA pressure at the indicated height differs from the true static pressure by c_pq to 0.2 %.
6. **Ducts**: Euler converges to Bernoulli at order 1.8–2.2; total pressure constant to 10⁻¹⁴; a constant-area duct keeps its pressure; the contraction's ends and midpoint; the venturi inversion to 10⁻¹⁴ and C_d scaling.
7. **Energy**: (8) is zero for a Bernoulli pair, equals h_L with a loss, is symmetric, and negative for a turbine; α₂ = 2 doubles the outlet kinetic head; power with efficiency; tunnel power ∝ V³, lower with a larger exit, and equal to the kinetic-energy flux with no losses.
8. `main` writes the figure.

### Compare

Two routes give the contraction's pressure: Bernoulli's algebra and the numerical integration of the differential equation (1) it comes from. They agree to 0.145 Pa on a 533.33 Pa drop at n = 200, and the difference falls as h² (orders 1.96, 1.98, 1.99): the trapezoid rule and the central-difference dV/dx are both second order. The airspeed conversions are checked the other way round: every function has an inverse written independently ((6) against (4) at sea level; `tas_from_cas` against `airspeeds_from_tas`; (7) against Bernoulli), and every round trip closes to 10⁻¹². The static-port altimeter error is checked against `isa_troposphere`'s pressure profile, a third route.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **The chain across altitude.** At a constant 60 m/s TAS, EAS and CAS at 0, 1500, 3000, 6000, 9000 m. Predict whether CAS − EAS grows. *Finding:* EAS = 60.000, 55.762, 51.689, 44.031, 37.020 m/s; CAS = 60.000, 55.799, 51.754, 44.136, 37.146 m/s (116.6, 108.5, 100.6, 85.8, 72.2 kt); CAS − EAS = 0, 0.037, 0.066, 0.105, 0.125 m/s. The dial falls by 44 kt between sea level and 9 km at the same true speed; the compressibility correction is a fraction of a knot for this aeroplane and matters only for fast ones.
2. **Static-port errors.** c_p = −0.04, −0.02, +0.02 at 30 m/s and 60 m/s TAS at sea level. Predict whether the relative error depends on speed. *Finding:* IAS − CAS = +1.15, +0.58, −0.58 kt at 30 m/s and +2.27, +1.14, −1.15 kt at 60 m/s: relative errors +1.97, +0.99, −1.00 % and +1.95, +0.98, −0.99 %, essentially −c_p/2 at both speeds. A port error is a percentage of airspeed; it matters most near the stall, where the margin is a few knots.
3. **Where Bernoulli fails.** The table in part 3 at 3000 m, and root finds for 1 % and 5 % (Problem 5). *Finding:* 0.10, 0.42, 1.16, 2.63, 4.69 % at 30–200 m/s; 1 % at M = 0.2824 and 5 % at M = 0.6279 at any altitude.
4. **Grid study.** The Euler integration at n = 25–200 (the table in part 3). Predict the order. *Finding:* 1.96, 1.98, 1.99, and the venturi inversion of the computed drop returns 15.000000 m³/s.
5. **Tunnel energy budget.** The open-circuit tunnel at 30 m/s with a 70 % fan: diffuser exits of 0.5, 0.75, 1.0, 1.5, 2.0 m² (losses 0.6 q), then losses K = 0, 0.3, 0.6, 1.0 (exit 1.0 m²). Predict which matters more. *Finding:* 18 576, 12 126, 9869, 8256, 7692 W against exit area; 2903, 6386, 9869, 14 513 W against K. Without a diffuser the fan throws away the whole test-section kinetic energy; doubling the exit area quarters that part. Already at an exit of 1.0 m² the losses (6966 W) exceed the discarded kinetic energy (2903 W), so a still larger diffuser buys little and the losses become the target; closed-circuit tunnels go further and recover the exit energy instead of discarding it.

Raw data format: CSV tables `h_m, TAS, EAS, CAS` (1); `TAS, cp, IAS_minus_CAS_kt, rel` (2); `TAS, M, bernoulli_err` (3); `n, max_err_Pa` (4); `A_exit_m2, K, P_W` (5); plus `pitot_static.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| q_c/q − 1 at M = 0.1826 | M²/4 = 0.00834 | 0.00836 | 2 × 10⁻⁵ | O(M⁴) terms |
| TAS from CAS | 60 m/s | 60.000000000000 | < 10⁻¹² | exact inversion |
| contraction drop | q(1 − 1/CR²) = 533.33 Pa | Euler, n = 200: 533.33 Pa | max 0.145 Pa along the duct | trapezoid, order 2 |
| venturi inversion | Q = 15 m³/s | 15.000000 | < 10⁻¹² | exact |
| port error c_p = −0.02 | ΔCAS ≈ 0.02q/(ρ₀CAS) | 0.99 kt | within 2 % of the linear estimate | O(c_p²) |
| altimeter error | −c_pq/(ρg) = 3.67 m | ISA pressure at that height | 0.2 % | ρ varies over 3.67 m |
| tunnel power, no losses, no diffuser | ½ρQV² = 8127 W | 8127 W | 0 | — |

The engineering answer: **at 3000 m and 60 m/s TAS the static port senses 70 108.5 Pa and the pitot 71 758.6 Pa; the indicator reads a CAS of 51.75 m/s (100.6 kt), 16 kt below the true airspeed, and EAS is 51.69 m/s; a static port in air with c_p = −0.02 adds 0.99 kt to the reading and 3.7 m to the altimeter; Bernoulli applied to the impact pressure is 0.42 % fast.** For this aeroplane the position error of the ports is the largest error in the chain, about twice Bernoulli's compressibility error, and it is the one flight test must measure.

What to improve: a position-error calibration c_p(α, flap) from flight test (Extension 2), then lag. What not to improve: the compressibility correction for this aeroplane; it is a tenth of a knot in CAS − EAS.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), (5) and (7); `pitot_static.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the aeroplane's airspeeds at cruise, the effect of a static-port error, and the condition of Bernoulli that fails first for a pitot tube and for the whole aeroplane.

### Adversarial self-checks

- **What if the speed goes to zero?** q and q_c → 0 and every airspeed → 0; `bernoulli_airspeed` returns 0, and a negative difference (a blocked or leaking line) is refused rather than turned into an imaginary speed.
- **What if the pitot tube ices over and its drain hole stays open?** The pitot pressure bleeds down to static: the indicator reads zero. If the drain also blocks, the trapped pressure stays constant and the indicator behaves like an altimeter, reading higher as the aeroplane climbs: a known accident chain.
- **Which assumption fails first in the real world?** The static pressure: the port's c_p depends on angle of attack and flap setting, so the position error changes in the approach. Then compressibility for faster aircraft.
- **What measurement could falsify the analysis?** A tower fly-by at known height and ground speed in still air: the altimeter error gives c_p directly (Δh = −c_pq/(ρg)), and the airspeed error must then match −c_p/2 of CAS; a mismatch points to a pitot error instead.
- **Does the result survive a change of route?** Bernoulli's algebra and the Euler integration agree to 0.03 % at n = 200; every conversion inverts to round-off.
- **Which decision changes if the uncertainty doubles?** If c_p is known to ±0.04 instead of ±0.02, the stall-speed uncertainty doubles to about ±2 % of IAS, roughly ±1 kt at the light aeroplane's approach speed: enough to change the approach speed a careful pilot flies.

## Check yourself

1. State the conditions under which Bernoulli's equation holds.
2. Define static, dynamic and total pressure.
3. Define TAS, EAS, CAS and IAS.
4. Why does a wing stall at the same EAS at every altitude, and why does that make IAS useful to pilots?
5. What does a static port sense if it sits where c_p = −0.02, and what do the airspeed indicator and the altimeter then show?
6. At what Mach number does incompressible Bernoulli make about a 1 % error in airspeed from a pitot tube?
7. What is head loss, and how do a pump and a turbine enter the energy equation?
8. Where does the power of an open-circuit wind tunnel's fan go, and why does a longer diffuser lower it?
