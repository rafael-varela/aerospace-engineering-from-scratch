# P08.17  Intro to Compressible Flow

> Air is incompressible only as long as it moves slowly compared with the speed of sound. Measure that ratio, the Mach number, and the total-to-static relations tell you how much Bernoulli is off. Above Mach 1 a shock stands in front of every pitot probe.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the speed of sound a = √(γRT) and the Mach number through a standard atmosphere, and derive the total-to-static relations T₀/T, p₀/p and ρ₀/ρ from the energy equation and the isentropic law. 2. Show that the stagnation pressure rise tends to Bernoulli's ½ρV² as M → 0, with the correction 1 + M²/4, and find the Mach number at which the incompressible airspeed of `P08.04` is wrong by a given amount. 3. Implement and test the isentropic relations, their inversion, the subsonic compressible pitot formula, the normal-shock relations and Rayleigh's supersonic pitot formula, each checked by an independent route (conservation laws, numerical integration of Euler's equation, inversion). 4. Decide when a flow must be treated as compressible for a given purpose, and interpret a pitot reading on a supersonic vehicle. |
| **Prerequisites** | `P08.01` (the ISA), `P08.04` (Bernoulli, pitot-static airspeed, IAS/CAS/EAS/TAS), `P08.03` (control-volume balances), `P07.05` (Isentropic Processes & Efficiencies; restated in 2.2), `P07.01` (Systems, Properties & State: the ideal gas), `P03.02` (root finding), `P03.06` (ODE integration with `solve_ivp`). Skill check: you can write p = ρRT and the steady-flow energy equation h + V²/2 = constant. |
| **Engineering question** | The light aeroplane cruises at 60 m/s at 3000 m: what is its Mach number, how wrong is the incompressible Bernoulli airspeed there, at what Mach number does that error reach 1 %, and what does a pitot probe read on a sounding rocket at Mach 2 at 11 km? The answers: M = 0.183 and an error of 0.42 % (0.25 m/s); 1 % at M = 0.282 (93 m/s at 3000 m); on the rocket a normal shock stands ahead of the probe, the pitot reads 5.640 times the static pressure (127.7 kPa), the Rayleigh formula recovers M = 2 exactly, and the subsonic formula would say M = 1.79 (−10.6 %). |
| **Assumptions and limits** | Air as a calorically perfect gas (γ = 1.4, R = 287.05287 J/(kg·K)); steady flow; adiabatic and reversible (isentropic) deceleration into a subsonic probe; a thin, normal, adiabatic shock ahead of a supersonic probe, followed by isentropic compression; the ISA to 47 km. Fails at high temperatures (γ varies above roughly 600 K, and air dissociates far above that), for oblique or detached shocks on probes at incidence, for probe Reynolds numbers so low that viscosity matters, and in unsteady flow. |
| **Analytical work** | The speed of sound; Mach number; the energy equation and T₀/T; the isentropic relations and p₀/p, ρ₀/ρ; the low-Mach series of p₀ − p and the M²/8 airspeed error; the subsonic pitot inversion; the normal-shock relations from mass, momentum and energy; Rayleigh's pitot formula; the total-pressure loss across a shock. |
| **Computational deliverable** | `starter/main.py`: `speed_of_sound`, `T0_over_T`, `p0_over_p`, `rho0_over_rho`, `mach_from_p0_over_p`, `pressure_coefficient_ratio`, `bernoulli_airspeed`, `subsonic_pitot_airspeed`, `bernoulli_error`, `mach_for_bernoulli_error`, `stagnation_temperature`, `normal_shock`, `rayleigh_pitot`, `mach_from_pitot` (the small ISA `isa` and `main` are given). Tests prove: the ISA reproduces its published anchor pressures, the hydrostatic equation and continuity at the layer boundaries; a² = γp/ρ and a ∝ √T; the three isentropic ratios are mutually consistent and invert exactly; numerical integration of Euler's equation along an isentropic streamline gives p₀/p to 10⁻⁸; the low-Mach limit is Bernoulli with the M²/4 and M⁴ corrections; the 1 % Mach number; T₀ from the energy equation; the compressible pitot returns the true airspeed while Bernoulli reads high; the normal shock conserves mass, momentum and energy to 10⁻¹², is the identity at M = 1, has the strong-shock limits and loses total pressure; Rayleigh's formula equals shock-then-isentropic compression and joins the subsonic formula at M = 1; the pitot inversion works on both branches; the results do not depend on the pressure unit. |
| **Practical activity** | **computational experiment**: the Bernoulli error as a function of Mach number against its M²/8 estimate; a check of the isentropic relation by integrating Euler's equation; the speed of sound through the ISA; a Mach sweep of the supersonic pitot reading, the shock loss and the error of the subsonic formula. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `compressible.png`, the derivations of T₀/T, p₀/p, the M → 0 limit and the normal-shock pressure ratio, the four experiment tables, and a half-page note on when the pitot-static formula must be compressible and how to read a supersonic pitot. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the total-to-static relations are derived with the assumption each needs; the M → 0 limit is shown; all public tests pass; the note gives the light aeroplane's Bernoulli error, the 1 % Mach number and the rocket's pitot interpretation with the shock loss. |
| **References** | `@nasa-glenn-beginners-guide-aeronautics` "Speed of Sound" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/speed-of-sound-interactive/), "Mach Number" (https://www.grc.nasa.gov/www/k-12/airplane/mach.html), "Isentropic Flow Equations" (https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html) (required-free: short statements of the relations). `@mit-ocw-16-01-unified-engineering-fluids` 16.03–16.04 notes F12 "Energy Conservation", F13 "Stagnation Quantities, Introduction to Shock Waves", F14 "Normal Shock Waves, Speed of Sound", F15 "Mach Number Relations, Normal Shock Properties", F16 "Shock Losses, Compressible Flow Pitot-Tube" (required-free: the derivations and the Rayleigh pitot formula). `@leishman-intro-aerospace-flight-vehicles` chapters "Mach Number & Reynolds Number": Definition of Mach Number; Significance of the Speed of Sound; and "Airspeed Definitions & Measurement": Airspeed Measurement at Higher Mach Numbers (required-free; link only). `@mit-ocw-16-100-aerodynamics` Lecture 39 "Derivation of Sound Wave Properties" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre39_cj/) (supplementary-free). `@libretexts-bar-meir-fluid-mechanics` Ch. 11 "Compressible Flow One Dimensional": 11.3.2 Speed of Sound in Ideal and Perfect Gases; 11.4.1 Stagnation State for Ideal Gas Model; 11.4.3.1 The pressure Mach number relationship (supplementary-free). `@ntrs-naca-tr-1135-compressible-flow` (https://ntrs.nasa.gov/citations/19930091059): the tables of isentropic and normal-shock ratios for γ = 7/5 (supplementary-free: a check table for your code). `@ntrs-us-standard-atmosphere-1976` 1.3.10 Speed of sound and Table III (supplementary-free). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics* (paid; library), the compressible-flow chapters. |
| **Extensions** | ★ Add the area–Mach relation A/A* and design the converging–diverging nozzle that runs a supersonic wind tunnel at M = 2. ★ Model a temperature probe with a recovery factor r (T_probe = T + rV²/(2c_p)) and find the static temperature from a measured probe temperature on the airliner. ★ Replace the perfect gas with γ(T) from a polynomial c_p(T) and find the Mach number at which the T₀/T error of the perfect-gas formula reaches 1 %. |
| **Time** | unmeasured |

## 1. The engineering problem

`P08.04` measured the light aeroplane's airspeed with a pitot-static tube and Bernoulli's equation, V = √(2(p₀ − p)/ρ). That equation assumes the density stays constant as the air is brought to rest in the probe. Air is compressible: slowing it down raises its pressure, and with the pressure its density and temperature. At 60 m/s the effect is small. At an airliner's cruise speed it is not, and above the speed of sound the probe no longer sees the free stream at all: a shock wave stands in front of it and destroys part of the total pressure the probe is meant to measure.

This lesson asks how large the effect is, when it can be ignored, and what replaces Bernoulli when it cannot. The vehicles are the light aeroplane at 60 m/s and 3000 m, an airliner at Mach 0.83, and a sounding rocket at Mach 2 at 11 km whose nose carries a pitot probe and whose skin carries a static port. Success criteria, written first: the speed of sound and the ISA reproduce their published anchor values; the isentropic relations are checked against an independent numerical integration of Euler's equation and recover Bernoulli as M → 0; the normal-shock relations conserve mass, momentum and energy to round-off; every inversion (pitot to Mach number, subsonic and supersonic) returns its input.

## 2. Understand and derive

### 2.1 The speed of sound and the Mach number

A weak pressure disturbance travels through a gas at the **speed of sound** a. A control volume moving with a weak wave (mass and momentum balances of `P08.03`) gives a² = (∂p/∂ρ) at constant entropy, since a sound wave is too fast and too weak to exchange heat or dissipate (MIT 16.100 Lecture 39; NASA Glenn "Speed of Sound"). For a perfect gas with p ∝ ρ^γ along an isentrope and p = ρRT,

**a = √(γp/ρ) = √(γRT).**   (1)

It depends only on temperature: **340.294 m/s** at sea level (288.15 K), 328.58 m/s at 3000 m, 295.07 m/s in the isothermal layer from 11 to 20 km, 303.13 m/s at 32 km and 329.80 m/s at 47 km, where the ISA warms again (US Standard Atmosphere 1976, 1.3.10 Speed of sound). The **Mach number** M = V/a compares the flow speed with the speed at which the flow can "warn" the air ahead of a change (NASA Glenn "Mach Number"; Leishman, "Mach Number & Reynolds Number"). The light aeroplane flies at M = 60/328.58 = **0.1826**.

### 2.2 Total conditions

Restated from `P07.05` (Isentropic Processes & Efficiencies): a process is **isentropic** when it is adiabatic and reversible; for a perfect gas along an isentrope, p/ρ^γ is constant, so p ∝ T^{γ/(γ−1)} and ρ ∝ T^{1/(γ−1)}. Restated from `P07.01`: the enthalpy of a calorically perfect gas is h = c_pT with c_p = γR/(γ − 1) = 1004.7 J/(kg·K) for air.

The steady-flow energy equation for an adiabatic flow with no work (`P08.04`) is h + V²/2 = h₀, constant along a streamline. The **total** (stagnation) state is the state the gas would reach if brought to rest. Dividing by c_pT and using (1):

**T₀/T = 1 + ((γ − 1)/2) M².**   (2)

This needs only an adiabatic flow; it holds even across a shock. If the deceleration is also reversible, the isentropic power laws give

**p₀/p = (1 + ((γ − 1)/2) M²)^{γ/(γ−1)},   ρ₀/ρ = (1 + ((γ − 1)/2) M²)^{1/(γ−1)}.**   (3)

(NASA Glenn "Isentropic Flow Equations"; Unified Engineering F13 "Stagnation Quantities"; Bar-Meir 11.4.1.) At M = 1, p₀/p = 1.8929 and T₀/T = 1.2. The total temperature rise is T₀ − T = V²/(2c_p): **1.79 K** at 60 m/s, **173 K** at Mach 2 at 11 km, where a probe tip reaches 390 K. Relation (3) inverts in closed form, M = √((2/(γ − 1))((p₀/p)^{(γ−1)/γ} − 1)), and that is the **compressible (subsonic) pitot formula**: from the impact pressure q_c = p₀ − p, the static pressure p and the static temperature T, M follows from (3) and V = Ma.

### 2.3 Where Bernoulli fails

Expanding (3) for small M with the binomial series and using ½ρV² = ½γpM²:

**p₀ − p = ½ρV² (1 + M²/4 + (2 − γ)M⁴/24 + …).**   (4)

As M → 0 the bracket tends to 1: **Bernoulli is the low-Mach limit of the isentropic relation**. The first correction is M²/4 in pressure. The incompressible airspeed, √(2q_c/ρ) with the true q_c and the true static density, is therefore too high by about **M²/8**. The exact error, `bernoulli_error`, is 0.125 %, 1.129 %, 3.154 %, 8.19 % and 12.94 % at M = 0.1, 0.3, 0.5, 0.8 and 1.0 (against M²/8 = 0.125, 1.125, 3.125, 8.0 and 12.5 %). For the light aeroplane: q_c = **1650.11 Pa** against ½ρV² = 1636.42 Pa, and the Bernoulli airspeed is **60.250 m/s, 0.417 % high**; the compressible formula returns 60.000000 m/s. The error reaches **0.1 % at M = 0.0894, 1 % at M = 0.2824 and 5 % at M = 0.6279**, which at 3000 m means 29, 93 and 206 m/s. That is the origin of the "M < 0.3 is incompressible" rule. The density change it implies is larger, ρ₀/ρ − 1 = 4.6 % at M = 0.3 (Problem 6): "incompressible" is a statement about the accuracy wanted, not about the air.

### 2.4 The normal shock

Above M = 1 the air ahead of a probe cannot be warned that the probe is there, and the flow adjusts through a **shock**: a thin layer (a few mean free paths) across which the pressure, density and temperature jump. A control volume enclosing a normal shock, with no area change, gives three balances (Unified Engineering F14, F15):

ρ₁u₁ = ρ₂u₂,   p₁ + ρ₁u₁² = p₂ + ρ₂u₂²,   c_pT₁ + u₁²/2 = c_pT₂ + u₂²/2.   (5)

Solving with ρu² = γpM²:

**M₂² = (1 + ((γ − 1)/2)M₁²)/(γM₁² − (γ − 1)/2),  p₂/p₁ = 1 + 2γ(M₁² − 1)/(γ + 1),  ρ₂/ρ₁ = (γ + 1)M₁²/((γ − 1)M₁² + 2).**   (6)

The flow behind a normal shock is subsonic. T₀ is unchanged (the shock is adiabatic) but the shock is irreversible: entropy rises and the total pressure falls, p₀₂/p₀₁ = 0.9928, 0.9298, 0.7209, 0.3283 and 0.0617 at M₁ = 1.2, 1.5, 2, 3 and 5. As M₁ → ∞ the density ratio tends to (γ + 1)/(γ − 1) = 6 and M₂ to √((γ − 1)/(2γ)) = 0.378.

### 2.5 The supersonic pitot

A pitot probe in supersonic flow sits behind its own (detached, but near the opening nearly normal) bow shock. The air crosses the shock, slows to M₂, and is then brought to rest isentropically in the probe. So the pitot reads p₀₂ = (p₂/p₁)(p₀₂/p₂)p₁, which with (6) and (3) is **Rayleigh's pitot formula** (Unified Engineering F16 "Shock Losses, Compressible Flow Pitot-Tube"):

**p₀₂/p₁ = [(γ + 1)²M₁²/(4γM₁² − 2(γ − 1))]^{γ/(γ−1)} (1 − γ + 2γM₁²)/(γ + 1).**   (7)

At M₁ = 1 it equals the isentropic 1.8929, and the two branches join. Above that it rises monotonically, so it can be inverted for M₁ by a bracketing root-finder. For the rocket at M = 2 at 11 km (T = 216.65 K, p = 22 632.0 Pa, a = 295.07 m/s, V = 590.14 m/s): **p₀₂/p₁ = 5.6404 and the pitot reads 127.65 kPa**, while isentropic compression without a shock would have given 7.8245. Read with the subsonic formula, that pitot ratio gives **M = 1.788 (−10.6 %)**; read with Bernoulli, 759.7 m/s (+28.7 %). Read with (7), it gives M = 2.0000000000.

### 2.6 Where the model stops

- **The perfect gas.** γ = 1.4 holds while the vibrational modes of N₂ and O₂ are frozen. Above roughly 600 K c_p rises and γ falls; at M = 3 the probe tip reaches 625 K at 10 km (Problem 4), where the perfect-gas formulas begin to drift (Extension 3).
- **The probe.** A probe at incidence sees an oblique, not a normal, shock; a probe with a tiny opening at low density sees viscous effects; a probe in a boundary layer or wake measures the wrong stream. And the static port has a position error that `P08.04` showed must be calibrated.
- **Unsteady flow.** The relations are steady. A gust or a manoeuvre faster than the pneumatic lag of the probe and its tubing makes the reading lag the flow.
- **The ISA.** The real atmosphere differs from the standard by tens of kelvin; a in (1) must use the actual static temperature, which is measured with a total-temperature probe and corrected with (2).

## 3. Build the computational model

Open `starter/main.py`. `isa` (the layers of `P08.01`, reimplemented here so the module stands alone) and `main` are given.

### Step 1: sound speed and the isentropic relations

`speed_of_sound` is (1). `T0_over_T`, `p0_over_p` and `rho0_over_rho` are (2)–(3), vectorised; `mach_from_p0_over_p` is the closed-form inverse. `pressure_coefficient_ratio(M)` is (p₀ − p)/(½ρV²) = 2(p₀/p − 1)/(γM²). Write p₀/p − 1 as `expm1(γ/(γ − 1) · log1p((γ − 1)M²/2))`. The obvious `p0_over_p(M) - 1` loses half its digits at M = 10⁻⁴, which the low-Mach test found in the first version.

### Step 2: airspeeds

`bernoulli_airspeed` is √(2q_c/ρ); `subsonic_pitot_airspeed(qc, p, T)` inverts (3) with 1 + q_c/p and multiplies by a(T); `bernoulli_error` is √(pressure_coefficient_ratio) − 1; `mach_for_bernoulli_error` finds the M at which the error equals a given value with `brentq` (`P03.02`); `stagnation_temperature` is T + V²/(2c_p).

### Step 3: shock and supersonic pitot

`normal_shock(M1)` returns the dictionary of (6) plus p₀₂/p₀₁ = (p₂/p₁)(p₀/p)(M₂)/(p₀/p)(M₁), refusing M₁ < 1. `rayleigh_pitot` is (7). `mach_from_pitot(r)` uses the isentropic inverse for r ≤ 1.8929 and otherwise solves (7) = r on [1, 50] by `brentq`:

```python
if r <= float(p0_over_p(1.0, gamma)):
    return float(mach_from_p0_over_p(r, gamma))
return brentq(lambda M: float(rayleigh_pitot(M, gamma)) - r, 1.0, 50.0, xtol=1e-14)
```

Running `main()` prints

```
ISA 3000 m: T = 268.650 K, p = 70108.5 Pa, rho = 0.90912 kg/m^3, a = 328.578 m/s
light aeroplane at 60.0 m/s: M = 0.1826; impact pressure qc = 1650.11 Pa, rho V^2/2 = 1636.42 Pa
  Bernoulli airspeed 60.2504 m/s (error 0.4173 %, M^2/8 = 0.4168 %); compressible formula 60.000000 m/s
  total temperature rise V^2/(2 cp) = 1.792 K
  Bernoulli error reaches 0.1 % at M = 0.0894 (29.4 m/s at 3000 m, 30.4 m/s at sea level)
  Bernoulli error reaches 1 % at M = 0.2824 (92.8 m/s at 3000 m, 96.1 m/s at sea level)
  Bernoulli error reaches 5 % at M = 0.6279 (206.3 m/s at 3000 m, 213.7 m/s at sea level)
  sea level: a = 340.294 m/s; p0/p at M = 1: 1.89293

sounding rocket at M = 2.0 at 11 km: T = 216.65 K, p = 22632.0 Pa, a = 295.07 m/s, V = 590.14 m/s
  normal shock: M2 = 0.57735, p2/p1 = 4.5000, T2/T1 = 1.68750, p02/p01 = 0.72087
  Rayleigh pitot p02/p1 = 5.64044 -> pitot reads 127654.7 Pa; isentropic (no shock) p01/p1 would be 7.82445
  inverted: M = 2.0000000000
  subsonic formula on that reading: M = 1.7879 (-10.61 %); Bernoulli: V = 759.7 m/s (28.74 %)
  total temperature 389.97 K
```

and writes `compressible.png`: the Bernoulli error against M with M²/8 and the light aeroplane's cruise marked; the pitot-to-static ratio for the subsonic, the shockless and the Rayleigh cases; the ISA speed of sound to 47 km.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-17-intro-compressible-flow`. The tests:

1. **ISA**: sea level; the published anchor pressures at 11, 20, 32 and 47 km to 5 × 10⁻⁶ and ρ(11 km) to 10⁻⁴ (published data, US Standard Atmosphere 1976); dp/dh = −ρg₀ by differences in four layers; continuity at the layer boundaries; refusal above 47 km.
2. **Sound**: 340.294 m/s at sea level (the published ISA value); a² = γp/ρ; a ∝ √T; another gas (helium).
3. **Isentropic**: p₀/p = (ρ₀/ρ)(T₀/T) and the power law; the inverse returns M from 0 to 3.5; integrating dp/dV = −ρ(p)V with p/ρ^γ constant (`solve_ivp`, 10⁻¹¹) from the free stream to rest gives p₀/p to 10⁻⁸ at M = 0.3, 0.8 and 1.5.
4. **Low Mach**: the ratio equals 1 + M²/4 + (2 − γ)M⁴/24 to 10⁻⁸ for M ≤ 0.05 and 1 at M = 10⁻⁴; the error is monotone and starts at M²/8.
5. **Airspeed**: the 1 % Mach number is between 0.25 and 0.32; T₀ from the energy equation; the compressible pitot returns the true speed to 10⁻¹² at three altitudes while Bernoulli reads high.
6. **Shock**: mass, momentum, energy and the equation of state across the shock to 10⁻¹² at M = 1.2, 2, 4; subsonic downstream, total pressure lost; identity at M = 1; strong-shock limits; total-pressure loss growing with M.
7. **Rayleigh**: equals (p₂/p₁)·(p₀/p)(M₂) at four Mach numbers; equals the isentropic value at M = 1; reads below isentropic compression at M = 2; `mach_from_pitot` returns M on both branches (including 0.999, 1 and 1.001); the reading is monotone; a reading below static is refused.
8. **Units**: kPa and Pa give the same airspeed. **main** writes the figure.

### Compare

The established check is a numerical integration that knows nothing of the closed form: Euler's equation along a streamline, dp = −ρV dV, with ρ(p) from the isentrope, integrated by `solve_ivp` from the free stream to rest. It reproduces (3) to 6 × 10⁻¹³, 2 × 10⁻¹² and 2 × 10⁻¹² at M = 0.3, 0.8 and 1.5 (the last is the shockless isentropic value, a valid integral even though a probe would not see it). For the shock, the check is the conservation laws (5) evaluated with the dimensional states the function returns. For a table check, the NACA TR 1135 tables of isentropic and normal-shock ratios for γ = 7/5 list exactly these functions of M; comparing a few rows with your code is a good exercise. The ISA anchor values are published data.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Where Bernoulli fails.** The Bernoulli airspeed error from M = 0.01 to 1, against M²/8. Predict the Mach number for 1 %. *Finding:* 0.125 %, 1.129 %, 3.154 %, 8.19 % and 12.94 % at M = 0.1, 0.3, 0.5, 0.8, 1.0, so M²/8 is good to 1 % of itself up to M = 0.5; 1 % at M = 0.2824. The light aeroplane, at M = 0.183, sits at 0.42 %.
2. **The isentrope by integration.** Euler's equation integrated along an isentropic streamline at M = 0.3, 0.8, 1.5 (the Compare section). *Finding:* p₀/p = 1.064430, 1.524340, 3.671031, equal to (3) to 2 × 10⁻¹². The relation is not an extra assumption beyond Euler's equation and the isentrope; it is their exact integral.
3. **Speed of sound through the ISA.** a at 0, 11, 20, 32, 47 km. Predict whether a falls monotonically. *Finding:* 340.294, 295.069, 295.069, 303.131, 329.799 m/s: it follows the temperature, constant in the isothermal layer and rising above 20 km. The same TAS is a higher Mach number at 11 km than at sea level by 15 %.
4. **The supersonic pitot.** M₁ = 1.0–5.0: the shock's M₂ and p₀₂/p₀₁, the Rayleigh ratio, the shockless isentropic ratio and the Mach number the subsonic formula would report. *Finding:*

| M₁ | M₂ | p₂/p₁ | p₀₂/p₀₁ | p₀₂/p₁ (Rayleigh) | p₀/p (no shock) | subsonic formula reads |
|---|---|---|---|---|---|---|
| 1.0 | 1.0 | 1.0 | 1.0 | 1.8929 | 1.8929 | 1.0 |
| 1.2 | 0.84217 | 1.5133 | 0.99280 | 2.4075 | 2.4250 | 1.1945 |
| 1.5 | 0.70109 | 2.4583 | 0.92979 | 3.4133 | 3.6710 | 1.4494 |
| 2.0 | 0.57735 | 4.5 | 0.72087 | 5.6404 | 7.8244 | 1.7879 |
| 3.0 | 0.47519 | 10.333 | 0.32834 | 12.061 | 36.733 | 2.2769 |
| 5.0 | 0.41523 | 29.0 | 0.06172 | 32.654 | 529.09 | 2.9218 |

Just above M = 1 the shock is weak and the subsonic formula is nearly right (−0.5 % at 1.2); at M = 5 it under-reads by 42 %, because the probe recovers only 6 % of the free-stream total pressure.

Raw data format: CSV tables `M, err_bernoulli_pct, M2_over_8_pct` (1); `M, p0_p_formula, p0_p_integrated` (2); `h_m, T_K, a_m_s` (3); `M1, M2, p2_p1, p02_p01, p02_p1, p0_p_isentropic, M_subsonic_formula` (4); plus `compressible.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| ISA p at 11, 20, 32, 47 km | 22 632.06, 5 474.889, 868.019, 110.906 Pa (published) | 22 632.04, 5 474.877, 868.016, 110.906 | ≤ 4 × 10⁻⁶ | R and g₀ rounding |
| a at sea level | 340.294 m/s (published) | 340.294 | < 10⁻³ | — |
| p₀/p | Euler integrated along the isentrope | formula (3) | ≤ 2 × 10⁻¹² | integrator tolerance |
| low-Mach limit | 1 + M²/4 + (2 − γ)M⁴/24 | to 10⁻⁸ (M ≤ 0.05) | — | truncated series |
| compressible pitot | V = 60 m/s | 60.000000 m/s | < 10⁻⁶ | — |
| normal shock at M = 2 | M₂ = 0.57735, p₂/p₁ = 4.5 (exact) | 0.57735, 4.5000; balances to 10⁻¹² | — | — |
| Rayleigh inversion | M = 2 | 2.0000000000 | < 10⁻¹⁰ | brentq tolerance |

The relations are exact for their assumptions and the code reproduces them by independent routes. The engineering answer: **the light aeroplane flies at M = 0.18, where Bernoulli over-reads by 0.42 % (0.25 m/s); the error reaches 1 % at M = 0.28, about 93 m/s at 3000 m; and the sounding rocket's pitot at Mach 2 reads 5.64 times the static pressure behind a bow shock, which the Rayleigh formula turns back into M = 2 while the subsonic formula would report 1.79.** The uncertainty of a real air-data system lies in the static pressure (the static-port position error of `P08.04`) and in the static temperature, not in these relations.

What to improve: measure the static temperature (with a total-temperature probe and (2)) instead of taking it from the ISA; then account for γ(T) at high Mach number. What not to improve: the formulas for M < 0.3 when 1 % is good enough.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), (4) and (6); `compressible.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the light aeroplane's Mach number and Bernoulli error, the Mach number at which compressibility must be included for 1 % accuracy, and how to read the rocket's pitot with the shock loss.

### Adversarial self-checks

- **What if the speed doubles?** The light aeroplane at 120 m/s is at M = 0.37: the Bernoulli error grows fourfold (it goes as M²), to about 1.7 %, past the 1 % threshold.
- **What if the velocity goes to zero?** T₀ = T, p₀ = p, and (4) is exactly Bernoulli; the tests check the limit to 10⁻⁸.
- **Which assumption fails first in the real world?** The static pressure and temperature the formulas need: the static port's position error and the difference between the real and the standard atmosphere. On the rocket, the probe's angle of attack (an oblique shock) and the tip heating.
- **What measurement could falsify the model?** A calibrated wind-tunnel run at known M with a pitot and a static tap: the ratio should follow (3) below M = 1 and (7) above it. A ratio below (7) at supersonic speed means the probe is at incidence or its opening sits in a separated or viscous region.
- **Does the result survive a change of method?** The closed form, the integrated Euler equation, the conservation laws across the shock and both inversions agree to 10⁻¹⁰ or better.
- **Which decision changes if the uncertainty doubles?** If the static pressure were uncertain by 2 % instead of 1 %, the rocket's Mach number from the Rayleigh formula (whose ratio grows roughly as M² at M = 2) would be uncertain by about 1 % instead of 0.5 %, still adequate for trajectory reconstruction. The decision that does change is the subsonic one: for the light aeroplane a 2 % error in q_c is a 1 % airspeed error, as large as the compressibility correction itself.

## Check yourself

1. What is the speed of sound at 288.15 K, at 3000 m and at 216.65 K, and why does it depend only on temperature?
2. Write the relations for T₀/T, p₀/p and ρ₀/ρ and say which need the flow to be isentropic.
3. At what Mach number is the incompressible pitot airspeed 1 % too high, and how does the error scale with M?
4. Why does the compressible pitot formula reduce to Bernoulli as M → 0?
5. What does a normal shock in front of a supersonic pitot probe do, and what does Rayleigh's formula account for?
6. What is the pitot-to-static ratio at M = 1 from each formula, and why do they agree there?
7. Which quantities are conserved across a normal shock and which are not?
8. Which flight quantity does the steady energy equation predict without any isentropic assumption, and what is its value for the light aeroplane and for the rocket?
