# P08.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Propellant and Mars values are typical, not specifications.

## 1. A manometer and a propellant tank (s = 15: Δh = 35 mm)

**Assumptions.** Fluids at rest, uniform g, incompressible liquids, the manometer's legs filled with air of density 1.204 kg/m³.

(a) Δp = (998.2 − 1.204) × 9.80665 × 0.035 = **342.20 Pa**. Ignoring the air gives 342.61 Pa, **0.12 % high** (the ratio ρ_air/ρ_water). For an air-over-water manometer the correction is below anything you can read; for a mercury-under-water manometer it would not be.

(b) Along a tube at 10°, the reading is 35/sin 10° = **201.6 mm**. A 1 mm reading error is **9.78 Pa** on the vertical tube and **1.70 Pa** on the inclined one: the same pressure is spread over 5.76 times more scale. The inclination changes the reading, not the pressure.

(c) Pad: 2.5 × 10⁵ + 810 × 9.80665 × 1.5 = **2.619 bar**. Burn at 4g total: 2.5 × 10⁵ + 810 × 4 × 9.80665 × 1.5 = **2.977 bar**. Coast: the effective gravity in the tank's frame is zero, so the outlet sees only the ullage, **2.500 bar**, and nothing holds the liquid over the outlet: it can float away from it and the pump would ingest gas. Ullage motors give a small acceleration that settles the propellant over the outlet before the main engine relights.

**Check.** Dimensions: kg/m³ × m/s² × m = N/m² ✓. Free-fall limit: the hydrostatic head vanishes, as it must in an orbiting tank.

**Common mistakes.** Using the reading along an inclined tube as the vertical height; adding the thrust acceleration without gravity on the pad, or forgetting it is gravity *plus* thrust acceleration during the burn.

## 2. The standard atmosphere by hand

(a) A slab of area A between h and h + dh: the pressure below pushes up with p A, the pressure above pushes down with (p + dp) A, gravity pulls down with ρ g A dh. At rest the side faces carry only pressure, which cancels horizontally, because a fluid cannot sustain a shear stress without deforming: that is the definition of a fluid. Equilibrium: **dp/dh = −ρg**.

(b) With ρ = p/(RT): dp/p = −g₀ dh/(RT). In a gradient layer T = T_b + a(h − h_b), dh = dT/a, so ln(p/p_b) = −(g₀/(aR)) ln(T/T_b), giving **p = p_b (T/T_b)^(−g₀/(aR))**. In an isothermal layer T = T_b and ln(p/p_b) = −g₀(h − h_b)/(RT_b), giving **p = p_b exp(−g₀(h − h_b)/(RT_b))**. Limit: (T/T_b)^(−g₀/(aR)) = exp[−(g₀/(aR)) ln(1 + aΔh/T_b)] = exp[−(g₀/(aR))(aΔh/T_b − a²Δh²/(2T_b²) + …)] → exp(−g₀Δh/(RT_b)) as a → 0 ✓.

(c) p′(0) = −ρ₀g₀, so p ≈ p₀ − ρ₀g₀h + O(h²): the constant-density formula is the linear part. It reaches zero at h = p₀/(ρ₀g₀) = RT₀/g₀ = **8434.5 m**, which is exactly the scale height at T₀ (since p₀/ρ₀ = RT₀). The constant-density atmosphere is an atmosphere one scale height thick.

(d) ρ = p/(RT) ∝ (T/T₀)^(−g₀/(aR)) / T ∝ **(T/T₀)^(−g₀/(aR) − 1)**; with a = −0.0065 K/m the exponent is 5.2559 − 1 = **4.2559**.

**Check.** The code's `isa_pressure` satisfies dp/dh = −ρg₀ to 10⁻⁷ by central differences at eight heights (test `test_closed_form_satisfies_the_hydrostatic_equation`), and RK4 on the ODE reproduces it to 2.7 × 10⁻⁹ at dh = 100 m.

**Common mistakes.** Dropping the minus sign of the lapse rate (a = −6.5 × 10⁻³ K/m in the troposphere); integrating in geometric altitude with constant g and calling it the ISA.

## 3. A helicopter on Mars

(a) **R = 188.92 J/(kg·K)**; ρ = 610/(188.92 × 210) = **0.01538 kg/m³** (1.26 % of Earth's sea-level 1.225); **H = RT/g = 10 694 m**, against 8434.5 m on Earth at 288.15 K: lower gravity more than offsets the heavier gas and the colder air.

(b) At 10 km: 610 e^(−10000/10694) = **239 Pa**; at 21.9 km: **79 Pa**. An isothermal model is crude over 22 km (the real temperature falls with height), so treat these as estimates.

(c) `isa_density(h) = 0.01538` at **h = 31.0 km**: flying at the Martian surface is like flying at 31 km on Earth, far above any helicopter's ceiling.

(d) Re = ρVL/μ: at a similar tip speed (similar Mach number; CO₂'s speed of sound at 210 K is about 230 m/s) and a viscosity of the same order as Earth air's, a density 80 times lower gives a Reynolds number roughly 80 times lower, so the blade sections work in a low-Reynolds-number regime where laminar separation dominates (`P08.10`–`P08.12`).

**Common mistakes.** Using Earth's R or g; quoting the pressure without saying it rests on an isothermal assumption.

## 4. A hot day and density altitude (s = 15: ΔT = 17.5 K)

**Assumptions.** ISA + 17.5 K at every height (a uniform offset; the lapse rate is unchanged), dry air, hydrostatic.

(a) p = 70 108.5 Pa (the ISA value at 3000 m, fixed by the altimeter setting); T = 268.65 + 17.5 = **286.15 K**; ρ = p/(RT) = **0.8535 kg/m³**, **6.1 % below** the ISA value 0.9091 kg/m³ (σ = 0.6968).

(b) Density altitude: **3608 m**. The aeroplane takes off, climbs and produces engine power as if it were 600 m higher than its altimeter says.

(c) RK4 (dh = 10 m) from 101 325 Pa with T = T_ISA + 17.5 K reaches 70 108.5 Pa at **3182.2 m**; the closed form h = ((T₀ + ΔT)/a)[(p/p₀)^(−aR/g₀) − 1] gives **3182.2 m** (the two differ by under 1 mm, the interpolation error). The warm column is less dense, so the same pressure drop needs more height: the aeroplane is 182 m higher than it reads.

(d) Terrain clearance: the true height (here the altimeter under-reads, which is the safe direction on a hot day and the dangerous one on a cold day). Take-off and climb performance: density altitude. Air-traffic separation: pressure altitude, because every aeroplane's altimeter makes the same error at the same place.

**Check.** The true height exceeds 3000 m by roughly ΔT/T̄ × 3000 m = 17.5/278 × 3000 ≈ 189 m: the thickness of a layer is proportional to its mean temperature (the hypsometric equation, from integrating dp/p = −g dh/(RT)).

**Common mistakes.** Recomputing the pressure at 3000 m with the hot temperature (the altimeter fixes the pressure, not the height); confusing density altitude with true altitude.

## 5. A stratospheric balloon (s = 15: h = 21 500 m)

**Assumptions.** Gas at ambient pressure and temperature, envelope spherical, ISA.

(a) ρ_air = **0.06907 kg/m³**, ρ_He = p/(R_He T) = **0.009545 kg/m³**; V = 50/(0.06907 − 0.009545) = **840.0 m³**; r = **5.853 m**.

(b) Surface integration (n = 400): **568.967 N**; displaced weight by volume integration: **568.963 N**; ρgV at the centre: **568.962 N**. Errors of the surface integral against the volume integral for n = 25, 50, 100, 200: 1.13, 0.281, 0.0701, 0.0175 N, orders **2.00, 2.00, 2.00** (midpoint rule).

(c) The helium mass is ρ_He V = **8.02 kg**. At sea level it occupies **47.36 m³** (r = 2.24 m) and lifts (1.225 − ρ_He,0) × 47.36 = **50.0 kg**: the same. Both densities scale as p/T, so (ρ_air − ρ_He)V = m_He(ρ_air/ρ_He − 1) = m_He(R_He/R_air − 1), a constant. A zero-pressure balloon's lift does not change as it rises; it stops rising when the envelope is full and the gas starts spilling out of the bottom.

(d) ρ(h + r)/ρ(h − r) = **0.9981**: 0.19 % over the balloon's height. Because the density varies almost linearly across it, ρgV at the centre is accurate to 10⁻⁵ here (Archimedes' principle holds exactly as the weight of displaced fluid; using the centre density is the approximation).

**Common mistakes.** Treating the helium at sea-level density at float altitude; forgetting that the payload mass includes the envelope.

## 6. Which simplification fails first?

(a) At z = 47 km geometric, h = r₀z/(r₀ + z) = **46 655 m**, **345 m** less. Reading the ISA at 47 000 m instead of 46 655 m gives a pressure 4.27 % too low (110.906 Pa against 115.851 Pa).

(b) The constant-density pressure is 1 % low at **1249 m** and 5 % low at **2609 m**; it reaches zero at 8434.5 m.

(c) The error grows with height: a warmer layer has a larger scale height, so the pressure at a given height is higher than the ISA's, and the altimeter, which converts pressure to height with the ISA, reads low by an amount that accumulates layer by layer. For ISA + 10 K the pressure-altitude error is −101 m at 3 km, −369 m at 11 km, −764 m at 20 km and **−1238 m at 30 km**.

(d) At 30 km: the constant-density model is useless (it predicts negative pressure above 8.4 km); a 10 K warm stratosphere is a **1.2 km** error; the geopotential-geometric confusion is **141 m** (2.1 % in pressure). Correct the temperature first (from a radiosonde or the rocket's own temperature sensor), then convert geopotential to geometric height.

**Common mistakes.** Assuming the geopotential correction is negligible at all heights (it is 0.05 % at 3 km but 0.7 % at 47 km); treating the ISA as the actual atmosphere.

## Check-yourself answers

1. A fluid is a substance that deforms continuously under any shear stress, however small. So a fluid at rest can carry no shear stress; only pressure, which acts normal to every surface and is the same in all directions at a point.
2. dp/dh = −ρg, with h measured upward; pressure falls as you go up because each layer carries the weight of the fluid above it.
3. With T constant, dp/p = −g dh/(RT): the fractional pressure drop per metre is constant, which is the definition of exponential decay. H = RT/g = 6342 m at 216.65 K.
4. In a gas, viscosity comes from molecules carrying momentum between layers; they move faster when hotter, so μ rises with T. In a liquid, viscosity comes from intermolecular attraction, which thermal motion overcomes more easily when hotter, so μ falls.
5. Dynamic viscosity μ (Pa·s) relates shear stress to shear rate; kinematic viscosity ν = μ/ρ (m²/s) is momentum diffusivity. ν grows much faster with altitude (2.7 times at 11 km, 11 at 20 km) because ρ falls while μ changes little.
6. The pressure difference sets the vertical height of the column; a tube inclined at θ spreads that height over 1/sin θ more length, so a reading error of 1 mm along the tube is a smaller pressure error. The pressure itself is unchanged.
7. Yes, when ρ is the density of the displaced fluid averaged over the body's volume: the buoyancy is the weight of the displaced fluid, ∫ρg dV, which follows from integrating −p n over the surface and the divergence theorem with ∇p = ρg. Over a small body ρgV at the centre is accurate to second order.
8. Geopotential altitude is the height in a uniform gravity g₀ with the same potential energy per unit mass as the real height in the inverse-square field: h = r₀z/(r₀ + z). The ISA is defined in it so that dp/dh = −ρg₀ holds exactly.
