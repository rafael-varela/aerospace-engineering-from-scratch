# P04.07  Kinetic Theory of Gases

> Pressure is molecules hitting a wall; temperature is how fast they move; the speed of sound is a fixed fraction of that speed. The air an aircraft flies through is 10²⁵ molecules per cubic metre doing nothing more complicated than bouncing.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the pressure of an ideal gas from the momentum molecules transfer to a wall, p = nm⟨v²⟩/3, and combine it with p = nkT to identify temperature with mean molecular kinetic energy. 2. Use the Maxwell–Boltzmann speed distribution: compute its most probable, mean and rms speeds, the fraction of molecules above a given speed, and decide which gases a planet can hold. 3. Derive the molar heat capacities and γ = (f + 2)/f from equipartition, and the speed of sound a = √(γRT/M) with its ratio √(γ/3) to the rms speed. 4. Compute the mean free path and Knudsen number and decide where the continuum model of the air stops. 5. Implement the International Standard Atmosphere to 86 km from the hydrostatic equation and the ideal-gas law, test it against standard values, and run two seeded simulations (free molecules in a box; random elastic collisions) that reproduce p = nkT and relax any speed distribution to Maxwell–Boltzmann. |
| **Prerequisites** | `P04.03` (momentum and impulse), `P04.02` (kinetic energy), `P04.05` (escape speed), `P02.16` (probability densities, moments, sampling), `P02.15` (probability fundamentals), `P01.07` (vectors). Skill check: you can compute the mean and variance of a density by integration and draw normal samples with `numpy.random.default_rng`. |
| **Engineering question** | Starting from molecules, what are the pressure, density and speed of sound of the air at an airliner's cruise altitude, how does the speed of sound follow from the molecular speeds, and above what altitude does a 1 m vehicle stop seeing the air as a continuum? |
| **Assumptions and limits** | Ideal gas: point molecules (their volume is negligible), no forces between them except brief elastic collisions, classical equipartition. Air is treated as a single species of molar mass 28.96 g/mol with γ = 1.4 (rotation excited, vibration frozen). The ISA is a standard, not a forecast: real days differ by tens of kelvin. The ideal-gas law fails at high pressure and low temperature (near condensation); γ = 1.4 drifts as vibration wakes up (1.383 at 600 K, 1.343 at 1000 K for N₂ with a 3390 K vibrational temperature) and fails completely when molecules dissociate (re-entry); the continuum model fails when the mean free path approaches the body size (Kn above about 0.01). |
| **Analytical work** | Momentum transfer to a wall and p = nm⟨v²⟩/3; ⟨KE⟩ = (3/2)kT; the Maxwell–Boltzmann density and its three speeds; the escape fraction; equipartition, C_v, C_p and γ; a = √(γRT/M); the mean free path; the hydrostatic equation and the ISA layers. |
| **Computational deliverable** | `starter/main.py`: pressure and number density, the three speeds, `maxwell_boltzmann_pdf`, `fraction_above`, `sample_speeds`, `gamma_from_dof`, `molar_heat_capacities`, `speed_of_sound`, `mean_free_path`, `knudsen_number`, the seven-layer ISA (`isa_temperature`, `isa_pressure`, `isa_density`, `isa_speed_of_sound`, `hydrostatic_residual`), and two SIMULATIONS (`simulate_box_pressure`, `relax_to_maxwellian`). Tests (13) check momentum-transfer pressure = nkT, the order and scaling of the speeds, normalisation and moments of the density, the escape fractions, sampled moments, γ and the heat capacities, the speed of sound, the mean free path scaling, ISA standard values and continuity, the hydrostatic residual, the box simulation against nkT, and the collision simulation's exact conservation and its relaxation to Maxwell–Boltzmann. |
| **Practical activity** | **computational experiment**: two seeded simulations (pressure from free molecules in a box; relaxation of a single-speed gas by elastic collisions) compared with the analytic results, and a sweep of the ISA from 0 to 86 km computing the speed of sound, mean free path and Knudsen number. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `kinetic_theory.png`; the simulation tables (pressure ratio against N with its scatter; kurtosis against collisions per molecule); the ISA table to 86 km; a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: p = nm⟨v²⟩/3 and γ = (f + 2)/f derived; all public tests pass; sea-level and 11 km ISA values and the speeds of sound reproduced by hand to three figures; the note gives the cruise-altitude state and the continuum limit with the Knudsen number that sets it, and says what the simulations do and do not show. |
| **References** | `@openstax-university-physics-vol2` Ch. 2 "The Kinetic Theory of Gases": 2.1 "Molecular Model of an Ideal Gas", 2.2 "Pressure, Temperature, and RMS Speed", 2.3 "Heat Capacity and Equipartition of Energy", 2.4 "Distribution of Molecular Speeds" (required-free; https://openstax.org/books/university-physics-volume-2/pages/2-2-pressure-temperature-and-rms-speed). `@leishman-intro-aerospace-flight-vehicles` chapter "Atmospheric Properties & the ISA", https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/international-standard-atmosphere-isa/ (required-free; the ISA as aeronautical engineers use it). |
| **Extensions** | ★ Add vibrational excitation: give N₂ a vibrational temperature of 3390 K and compute C_v(T) and γ(T) from 300 to 3000 K; find where γ = 1.4 is 5 % wrong. ★ Replace the random-pair collisions by a hard-sphere molecular-dynamics simulation in 2D and measure the mean free path directly. ★ Build a "non-standard day" ISA (ISA + 15 K) and compute how much it changes an aircraft's true airspeed at Mach 0.85 and 11 km. |
| **Time** | unmeasured |

## 1. The engineering problem

An airliner cruises at 11 km at Mach 0.85. Its pitot tube reads a pressure, its static port another, and its air-data computer turns them into airspeed, altitude and Mach number using a model of the atmosphere: the International Standard Atmosphere (ISA). Mach number is speed divided by the speed of sound, and the speed of sound at 11 km is not the 340 m/s of sea level; it is 295 m/s, because the air there is 71 K colder. Why should a colder gas carry sound more slowly? Why does pressure fall by a factor of 4.5 between the ground and 11 km while temperature falls by only 25 %? And why does a re-entering capsule at 80 km, in air whose pressure is one ten-thousandth of sea level, still need a heat shield, while above about 100 km aerodynamics as an aircraft designer knows it stops working?

All of these answers come from treating the air as what it is: an enormous number of molecules (2.5 × 10²⁵ per cubic metre at sea level) flying in straight lines between collisions. This lesson derives pressure from the momentum they carry to a wall, temperature from their kinetic energy, the distribution of their speeds, the heat capacities and γ from the number of ways a molecule can hold energy, and the speed of sound from all of it. It ends at the ISA, built from the hydrostatic equation and the ideal-gas law, and at the Knudsen number, which says where the continuum picture fails.

Success criteria, written first: momentum-transfer pressure equals nkT to 10⁻¹⁰; the speed of sound at 288.15 K is 340.3 m/s and equals √(1.4/3) v_rms; the ISA matches its standard values at every layer boundary to 10⁻³ and satisfies the hydrostatic equation to 10⁻⁵; a simulation of free molecules reproduces nkT within its statistical scatter, and a collision simulation conserves energy and momentum to rounding while relaxing to Maxwell–Boltzmann.

## 2. Understand and derive

### 2.1 Pressure from momentum transfer

Put N molecules of mass m in a cubic box of side L (volume V = L³, number density n = N/V). **Assumptions:** the molecules are points; they collide elastically with the walls; between wall collisions they move in straight lines (collisions with each other only shuffle velocities among molecules, which does not change the averages used below).

A molecule with x-velocity v_x hits the wall at x = L, reverses v_x, and transfers momentum 2m|v_x| to it. It comes back to the same wall after travelling 2L in x, a time 2L/|v_x|. The average force it exerts on that wall is momentum per hit divided by time between hits: 2m|v_x| × |v_x|/(2L) = mv_x²/L. Summing over all molecules and dividing by the wall's area L²:

p = Σ m v_x² / L³ = n m ⟨v_x²⟩.

There is nothing special about x: the gas is **isotropic**, so ⟨v_x²⟩ = ⟨v_y²⟩ = ⟨v_z²⟩ = ⟨v²⟩/3. Hence

**p = (1/3) n m ⟨v²⟩ = (1/3) n m v_rms²**,

with v_rms = √⟨v²⟩. Check units: m⁻³ × kg × m²/s² = kg m⁻¹ s⁻² = Pa.

### 2.2 Temperature is kinetic energy

The ideal-gas law, measured long before molecules were accepted, says **p = nkT**, with Boltzmann's constant k = 1.380649 × 10⁻²³ J/K. Equate it to 2.1: nkT = (2/3) n ⟨½mv²⟩, so

**⟨½ m v²⟩ = (3/2) k T**.

The absolute temperature is a measure of the mean translational kinetic energy of the molecules, and nothing else. Per mole (multiply by Avogadro's number N_A, and use R = kN_A = 8.314 J/(mol K), M = mN_A):

**v_rms = √(3RT/M)**.

For air (M = 28.96 g/mol) at 288.15 K: v_rms = 498.1 m/s. At sea level the number density is n = p/kT = 2.547 × 10²⁵ m⁻³, each molecule weighs 4.810 × 10⁻²⁶ kg, and (1/3)nm v_rms² = 101 325.0 Pa, identical to nkT (the code checks this). The rate at which molecules strike each square metre of wall, n⟨v⟩/4 (Problem 2), is 2.9 × 10²⁷ per second.

### 2.3 The Maxwell–Boltzmann distribution

Not every molecule moves at v_rms. In equilibrium each velocity component is normally distributed with mean zero and variance kT/m (the derivation, from the requirement that the distribution be isotropic and that the three components be independent, is in OpenStax 2.4). The speed is the length of a three-dimensional normal vector, and its density is

**f(v) = 4π (m/2πkT)^{3/2} v² exp(−mv²/2kT)**.

The v² factor counts the volume of a thin spherical shell in velocity space (there are more ways to have speed v as v grows); the exponential is the Boltzmann factor, which makes high energies rare. Three characteristic speeds follow by calculus (Problem 2):

- **most probable** (the peak): v_p = √(2RT/M) = 406.7 m/s for air at 288.15 K;
- **mean**: ⟨v⟩ = √(8RT/πM) = 458.9 m/s;
- **rms**: v_rms = √(3RT/M) = 498.1 m/s.

They stand in the fixed ratios 1 : 1.128 : 1.225, at every temperature and for every gas. All three scale as √(T/M): four times the temperature doubles them; a molecule sixteen times lighter (H₂ against O₂) is four times faster.

The density integrates to 1 (the code does it with `scipy.integrate.quad`: 1.000000), and its tail is what matters for many applications. At sea level 0.71 % of air molecules move faster than 1000 m/s.

### 2.4 Why the Earth kept its nitrogen and lost its hydrogen

At the exobase, about 500 km up, the thermosphere's temperature is around 1000 K, molecules rarely collide, and one moving upward faster than the escape speed (`P04.05`: √(2μ/r) = 10.77 km/s at 500 km altitude; take 10.8 km/s) leaves the planet. The fraction of molecules above 10.8 km/s at 1000 K:

| gas | v_rms at 1000 K | fraction above 10.8 km/s |
|---|---|---|
| H₂ | 3518 m/s | 3.2 × 10⁻⁶ |
| He | 2496 m/s | 4.4 × 10⁻¹² |
| N₂ | 944 m/s | 7.4 × 10⁻⁸⁵ |
| O₂ | 883 m/s | 5.7 × 10⁻⁹⁷ |

Three parts per million sounds small, but the escaping molecules are continuously replaced from the Maxwell–Boltzmann tail, collisions happen constantly below the exobase, and the Earth has 4.5 × 10⁹ years: hydrogen and helium leak away; nitrogen and oxygen are retained forever. The exponential tail turns a factor of 3.7 in rms speed into 79 orders of magnitude in escape fraction. Problem 3 asks the same question of Mars.

### 2.5 Equipartition, heat capacity and γ

The result ⟨½mv_x²⟩ = ½kT generalises: in thermal equilibrium every independent quadratic term in a molecule's energy (a "degree of freedom") holds on average **½kT**. This is the **equipartition theorem** (OpenStax 2.3). A monatomic gas (He, Ar) has three translational degrees of freedom: U = (3/2)RT per mole. A diatomic molecule (N₂, O₂) also rotates about the two axes perpendicular to its bond (rotation about the bond itself stores negligible energy), adding two: f = 5, U = (5/2)RT. Its vibration would add two more (kinetic and potential), but at room temperature the vibrational quantum is so large (a characteristic temperature of about 3390 K for N₂, used as an input in the extension) that vibration is "frozen out" (it adds 0.0016R to C_v at 300 K in the Einstein model): a quantum effect classical equipartition cannot see.

The molar heat capacity at constant volume is **C_v = dU/dT = fR/2**. At constant pressure the gas also does work pdV = RdT per mole as it expands, so **C_p = C_v + R**, and the ratio

**γ = C_p/C_v = (f + 2)/f**:

5/3 = 1.6667 for a monatomic gas, 7/5 = 1.400 for air (C_v = 20.79 and C_p = 29.10 J/(mol K)), and 9/7 = 1.29 if vibration were fully excited.

### 2.6 The speed of sound

Sound is a small pressure wave. It travels at **a = √(γp/ρ) = √(γRT/M)**, where the γ (rather than 1) appears because the compressions are too fast for heat to flow out of them: they are adiabatic. That formula comes from fluid mechanics; kinetic theory adds the insight. Divide by v_rms = √(3RT/M):

**a / v_rms = √(γ/3)**.

For air, √(1.4/3) = 0.6831: sound travels at 68 % of the rms molecular speed, and for exactly the reason one would hope, because the molecules themselves carry the disturbance. At 288.15 K, a = 340.3 m/s; at 216.65 K (11 km), 295.1 m/s. The speed of sound depends on temperature and composition only, not on pressure or density. A gas cannot carry a pressure signal faster than its molecules move, which is why a = 0.68 v_rms < ⟨v⟩.

### 2.7 Mean free path and the Knudsen number

Model molecules as hard spheres of diameter d. A molecule moving through the gas sweeps a cylinder of cross-section πd²; it collides once per volume 1/n swept, and accounting for the motion of the others adds a factor √2. The **mean free path** is

**λ = 1/(√2 π d² n) = kT/(√2 π d² p)**.

With d = 3.7 × 10⁻¹⁰ m for N₂, λ = 64.6 nm at sea level, and each molecule collides 7.1 × 10⁹ times per second. The **Knudsen number** Kn = λ/ℓ compares it with a body size ℓ. Below Kn ≈ 0.01 the gas is a continuum: pressure, density and velocity are smooth fields, and the Navier–Stokes equations of later phases apply. Above Kn ≈ 10 the molecules hit the body without meeting each other first (free-molecular flow: satellite drag in low orbit is computed this way). In between lies the transition regime that re-entry vehicles cross.

### 2.8 The International Standard Atmosphere

The ISA (Leishman's "Atmospheric Properties & the ISA" chapter) is a defined model: a temperature profile made of straight-line segments in geopotential altitude, sea-level values T₀ = 288.15 K and p₀ = 101 325 Pa, air as an ideal gas with R_air = R/M = 287.05 J/(kg K), and g₀ = 9.80665 m/s². Pressure follows from two equations. **Hydrostatic balance** of a thin slab of air (its weight is supported by the pressure difference across it): dp/dh = −ρg₀. **The ideal-gas law:** ρ = p/(R_air T). Together, dp/p = −g₀ dh/(R_air T). In a layer where T = T_b + L(h − h_b) (lapse rate L), this integrates to

**p = p_b (T_b/T)^{g₀/(L R_air)}**,

and in an isothermal layer (L = 0) to

**p = p_b exp(−g₀(h − h_b)/(R_air T_b))**:

exponential decay with the **scale height** H = R_air T/g₀ (8.4 km at 288.15 K). The ISA layers to 86 km:

| from (km) | to (km) | lapse rate (K/km) | name |
|---|---|---|---|
| 0 | 11 | −6.5 | troposphere |
| 11 | 20 | 0 | tropopause |
| 20 | 32 | +1.0 | stratosphere |
| 32 | 47 | +2.8 | stratosphere |
| 47 | 51 | 0 | stratopause |
| 51 | 71 | −2.8 | mesosphere |
| 71 | 86 | −2.0 | mesosphere |

Each layer starts from the temperature and pressure at the top of the one below, so the model is continuous. The code's table, with kinetic-theory quantities alongside:

```
 h (km)    T (K)     p (Pa)  rho (kg/m^3)  a (m/s)  lambda (m)   Kn (1 m)
      0   288.15  101325.00    1.2250e+00    340.3    6.46e-08   6.46e-08
      5   255.65   54020.14    7.3611e-01    320.5    1.07e-07   1.07e-07
     11   216.65   22632.29    3.6392e-01    295.1    2.17e-07   2.17e-07
     20   216.65    5475.00    8.8036e-02    295.1    8.98e-07   8.98e-07
     32   228.65     868.05    1.3225e-02    303.1    5.98e-06   5.98e-06
     47   270.65     110.91    1.4276e-03    329.8    5.54e-05   5.54e-05
     71   214.65       3.96    6.4215e-05    293.7    1.23e-03   1.23e-03
     86   184.65       0.30    5.7043e-06    272.4    1.39e-02   1.39e-02
```

Two things to notice. Pressure falls by 4.5 between 0 and 11 km while T falls by 25 %: pressure is exponential in altitude with a scale height of a few kilometres, temperature is linear. And for a 1 m body the Knudsen number reaches 0.01 at 84.1 km, at the top of the ISA's range: above it, the continuum model of this table is itself the wrong model, which is one reason the ISA stops at 86 km.

### 2.9 Where the model fails

- **Dense or cold gas.** Molecules occupy volume and attract each other, so the ideal-gas law drifts at the pressures of high-pressure gas tanks and fails near condensation. Not an issue in flight; an issue in pressurant and propellant tanks.
- **Hot gas.** Vibration excites gradually: with an Einstein model and a 3390 K vibrational temperature, γ of N₂ is 1.383 at 600 K, 1.343 at 1000 K and 1.294 at 3000 K, and it is 5 % below 1.4 at 1191 K; behind a re-entry shock O₂ and N₂ then dissociate. The extension computes this.
- **Moist air.** Water vapour (M = 18 g/mol) makes humid air *less* dense than dry air at the same p and T.
- **Real atmosphere.** The ISA is a standard day; real temperatures at 11 km vary by ±20 K or more with latitude, season and weather, and aircraft performance is corrected to ISA conditions for that reason.
- **Rarefied gas.** Above Kn ≈ 0.01 the continuum equations need slip corrections; above 10, free-molecular theory.

## 3. Build the computational model

Open `starter/main.py`. The constants (k, N_A, R, the molar masses, the ISA layer table) are given.

### Step 1: pressure and the three speeds

`pressure_from_momentum_transfer(n, m, v_rms)`, `ideal_gas_pressure`, `number_density`, `rms_speed`, `mean_speed`, `most_probable_speed`: one line each.

### Step 2: the distribution

```python
def maxwell_boltzmann_pdf(v, T, M):
    v = np.asarray(v, dtype=float)
    a = (M / N_A) / (2.0 * K_B * T)
    return 4.0 * np.pi * (a / np.pi) ** 1.5 * v ** 2 * np.exp(-a * v ** 2)
```

`fraction_above(v_min, T, M)` integrates the density from v_min to ∞ with `scipy.integrate.quad`; with `limit=200` it resolves tails down to 10⁻¹⁰⁰. `sample_speeds(n, T, M, rng)` draws three normal components of variance kT/m and takes the norm: that *is* the Maxwell–Boltzmann distribution, by construction.

### Step 3: equipartition and sound

`gamma_from_dof(f)`, `molar_heat_capacities(f)`, `speed_of_sound(T, gamma, M)`, `sound_to_rms_ratio(gamma)`, `mean_free_path(T, p, d)`, `collision_frequency`, `knudsen_number`.

### Step 4: the ISA

Loop over the seven layers carrying the base temperature and pressure upward; inside each, apply the lapse-rate or isothermal formula of 2.8 to the altitudes in that layer. Return NaN outside 0–86 km rather than extrapolating. `hydrostatic_residual(h)` differentiates `isa_pressure` numerically and returns (dp/dh + ρg₀)/(ρg₀), which should be zero to the finite-difference error.

### Step 5: two simulations

Both are **SIMULATIONS**: computer experiments on model molecules, labelled as such wherever their output appears. Neither is a measurement of a real gas.

`simulate_box_pressure(N, T, M, L, t_total, rng)`: draw |v_x| for N molecules from the normal distribution; each hits the wall at x = L every 2L/|v_x|; count the hits in t_total from a random starting phase; sum 2m|v_x| per hit; divide by t_total and L². With L = 1 μm and t_total = 0.1 μs each molecule hits about 24 times. The result is the *time-averaged* momentum flux, and it should equal nkT within a statistical scatter that falls as 1/√N.

`relax_to_maxwellian(velocities, n_collisions, rng)`: pick a random pair (i, j); keep their centre-of-mass velocity and the magnitude g of their relative velocity; give the relative velocity a uniformly random new direction. That is an elastic collision between equal masses (momentum and kinetic energy conserved exactly, the scattering angle random). Start every molecule at the same speed in random directions, the most non-Maxwellian gas with a given energy, and watch the distribution of speeds become Maxwell–Boltzmann at the temperature (3/2)kT = ⟨½mv²⟩.

Running `main()` prints the sea-level kinetic quantities, the escape table, the ISA table, the two simulations, and writes `kinetic_theory.png` (Maxwell–Boltzmann densities at three temperatures and for H₂; the relaxed simulation's histogram against the density; T, a and p against altitude).

### Verify

`python -m pytest curriculum/phase-04/lesson-07-kinetic-theory-of-gases` runs 13 tests: momentum-transfer pressure equals nkT and quadruples with v_rms doubling; v_p < ⟨v⟩ < v_rms in the exact ratios, √T and 1/√M scaling, air near 500 m/s; the density's normalisation, first and second moments, peak and ⟨KE⟩ = (3/2)kT; the escape fractions (H₂ between 10⁻⁷ and 10⁻⁵, N₂ below 10⁻⁶⁰); 2 × 10⁵ samples reproduce ⟨v⟩, v_rms and a tail fraction (tolerances justified in the test as 5σ); γ and C_p − C_v = R; a = 340.3 and 295.1 m/s, a/v_rms = √(γ/3), helium faster; λ ∝ T/(pd²); ISA values at 0, 11, 20, 32, 47, 71, 86 km, monotone p and ρ, continuity at the boundaries, NaN above 86 km; the hydrostatic residual below 10⁻⁵ and the isothermal layer's exponential; the box simulation within 2 % of nkT at 10⁵ molecules, scattering by several percent at 100; the collision simulation conserving energy to 10⁻¹² and momentum to rounding, keeping T, and reaching the Maxwell–Boltzmann ratios ⟨v⟩/v_rms = √(8/3π) and ⟨v⁴⟩/⟨v²⟩² = 5/3 within 2 %, with a histogram check.

### Compare

Two independent routes to the same numbers. The ISA's pressure at the layer boundaries matches the published standard values (22 632 Pa at 11 km, 5475 Pa at 20 km, 868 Pa at 32 km) to the four or five figures tested, and its hydrostatic residual is below 10⁻⁵ everywhere: the closed-form layer integrals and the differential equation agree. The simulations reproduce the analytic results by an entirely different path (counting impacts, colliding pairs), within statistical tolerances that are themselves predicted. When the two routes agree, both are probably right.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it. Both simulations use `numpy.random.default_rng(seed)`; record the seed.

1. **Pressure from impacts.** Run `simulate_box_pressure` for N₂ at 300 K with N = 10³, 10⁴, 10⁵, forty seeds each, and tabulate the mean and the standard deviation of p_sim/(nkT). *Predict:* mean 1, scatter ∝ 1/√N.
2. **Relaxation.** Start 2 × 10⁴ molecules at one speed and apply 0.25, 0.5, 1, 2, 4 and 8 pair collisions per molecule. Track ⟨v⁴⟩/⟨v²⟩² (1 for a single speed, 5/3 for Maxwell–Boltzmann) and ⟨v⟩/v_rms (1, then √(8/3π) = 0.921). *Predict* how many collisions per molecule it takes to come within 1 % of 5/3, then compare with the sampling scatter of the same statistic for an exact Maxwell–Boltzmann sample of 2 × 10⁴ molecules.
3. **Escape.** Tabulate the escape fraction of H, H₂, He, N₂ at 1000 K (Earth) and 200 K (Mars, with its own escape speed from `P04.05`). *Predict* which gases each planet retains.
4. **The ISA sweep.** From 0 to 86 km in 1 km steps compute T, p, ρ, a, λ and Kn for ℓ = 1 m and 4 m. *Predict* the altitude at which Kn = 0.01 for each, and where the speed of sound is lowest.

Raw data format: CSV files `N, seed, p_ratio` (1), `collisions_per_molecule, kurtosis, mean_over_rms, T_kin` (2), `gas, T_K, v_esc_mps, fraction` (3), `h_m, T_K, p_Pa, rho, a_mps, lambda_m, Kn_1m, Kn_4m` (4); plus `kinetic_theory.png`. Every file from items 1–2 has `SIMULATION` in its name.

What you should find for items 1 and 2 (seeded, from the reference code): the pressure ratio's mean over forty seeds is 0.994, 0.998 and 1.000 with standard deviations 4.6 %, 1.3 % and 0.47 % at N = 10³, 10⁴, 10⁵ (the single 10⁴ run printed by `main()` happens to be 2.6 % high, a 2σ draw); the kurtosis reaches 1.6689 after ten pair collisions per molecule (`main()`), against the sampling scatter of an exact Maxwell–Boltzmann sample of the same size, which is about 0.009 (Problem 5).

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| p from momentum transfer at sea level | nkT = 101 325 Pa | 101 325.0 Pa | < 10⁻¹⁰ | identity |
| a at 288.15 K | √(γRT/M) | 340.3 m/s | — | — |
| a / v_rms | √(1.4/3) = 0.6831 | 0.6831 | 0 | — |
| γ (diatomic, f = 5) | 7/5 | 1.400 | 0 | — |
| ISA p at 11 km | 22 632 Pa (standard) | 22 632.29 Pa | 10⁻⁵ | rounding of the standard |
| ISA hydrostatic residual | 0 | ≤ 7.5 × 10⁻⁶ | finite difference across a layer kink at 11 km | numerical derivative |
| box simulation, N = 10⁵ | 1 | 1.005 (one seed) | 0.5 % | statistical, σ = 0.47 % |
| collision simulation energy | conserved | ratio 1.000000000000 | < 10⁻¹² | exact by construction |
| collision simulation momentum | conserved | 4.75 × 10⁻¹⁰ m/s change in the sum | rounding | — |
| relaxed ⟨v⁴⟩/⟨v²⟩² | 5/3 = 1.6667 | 1.6689 | 0.13 % | sampling (σ ≈ 0.5 %) |
| relaxed most probable speed | v_p = 422.0 m/s at 300 K | histogram peak in the 400 m/s bin | one bin width | histogram resolution |

Every analytic quantity is reproduced to the precision of the arithmetic; every simulated one to within its statistical scatter, which is predicted in advance and is the right yardstick. A simulation that agreed with theory *better* than its scatter allows would be suspicious, not reassuring.

What to improve. The air model's weakest point for an aircraft is not kinetic theory but the standard day: on a hot day at 11 km the air can be 15 K warmer than ISA, raising a by 3.4 % and lowering ρ at the same pressure by 6.5 %. For a re-entry vehicle, it is γ: above roughly 600 K vibration is excited and γ falls, and above a few thousand kelvin the gas dissociates. The next model to build is one with a temperature-dependent γ (★ extension; 5 % low already at 1191 K), not a better ISA.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of p = nm⟨v²⟩/3, ⟨KE⟩ = (3/2)kT, γ = (f + 2)/f and a/v_rms = √(γ/3); `kinetic_theory.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) answering the engineering question for 11 km: T, p, ρ, a, the true airspeed at Mach 0.85, the mean free path, and the altitude above which a 1 m vehicle leaves the continuum regime, with the evidence.

### Adversarial self-checks

- **What if the temperature doubles?** All molecular speeds and the speed of sound rise by √2; at fixed n the pressure doubles; at fixed p the density halves and the mean free path doubles.
- **What if the molecular mass doubles at the same temperature?** Speeds fall by √2, the speed of sound too (same γ); the pressure at fixed n is unchanged (⟨½mv²⟩ = (3/2)kT does not depend on m).
- **Which assumption fails first in the real world?** For aircraft, the standard day; for re-entry, constant γ, then the continuum; for a satellite at 400 km, everything continuum (Kn ≫ 10).
- **What measurement could falsify the model?** A speed of sound that depends on pressure at fixed temperature (it does not, to a fraction of a percent, in air near ambient conditions), or a C_v of air near 7R/2 at room temperature (vibration frozen out says 5R/2).
- **Does the result survive a change in method?** Pressure from impacts and from p = nkT; speeds from the integral, from samples, and from a collision simulation that never assumed the distribution: all agree within their stated tolerances.
- **Which decision changes if T at 11 km is 20 K higher than ISA?** The Mach number for a given true airspeed falls by 4.3 %; the engine's thrust and the wing's lift at a given indicated airspeed change through ρ; the ISA-based performance charts must be corrected.

## Check yourself

1. Write the pressure of an ideal gas in terms of n, m and v_rms. Where does the factor 1/3 come from?
2. What is the mean translational kinetic energy of a molecule at temperature T?
3. Put v_p, ⟨v⟩ and v_rms in order and give their ratios.
4. Why does the Earth retain nitrogen but not hydrogen?
5. Give γ for a monatomic gas and for air near room temperature, and explain the difference in terms of degrees of freedom.
6. What is the speed of sound in air at 11 km in the ISA, and why is it lower than at sea level?
7. Define the Knudsen number. Roughly what value separates continuum from transitional flow?
8. In the collision simulation, which quantities are conserved exactly by each collision, and which quantity is not but relaxes to a fixed value?
