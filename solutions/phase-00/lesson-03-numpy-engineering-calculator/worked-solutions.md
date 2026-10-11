# P00.03  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Numbers below are for the example parameters stated; yours will differ if your inputs differ.

## 1. The air where you live (example: z_home = 500 m)

**Assumptions.** ISA troposphere; dry air; geopotential altitude from a spherical Earth, R_E = 6 356 766 m.
**Units.** z, h in m; T in K; p in Pa; ρ in kg/m³; a in m/s.
**Derivation.** h = 6 356 766 × 500/(6 357 266) = 499.96 m (a 4 cm correction, negligible here but not at 11 km). T = 288.15 − 0.0065 × 499.96 = 284.90 K. p = 101 325 × (284.90/288.15)^5.2559 = 101 325 × 0.94213 = 95 461 Pa. ρ = 95 461/(287.05287 × 284.90) = 1.1673 kg/m³. a = √(1.4 × 287.05287 × 284.90) = 338.37 m/s.
**Result.** Code: T = 284.900 K, p = 95 460.8 Pa, ρ = 1.16727 kg/m³, a = 338.369 m/s. Relative differences all below 1 × 10⁻⁵, i.e. agreement to five figures; the hand value is limited by rounding the exponent and the base, not by physics.
**Check.** Rule of thumb: pressure falls about 1.2 % per 100 m near sea level; 500 m gives −5.8 %, and 95 461/101 325 = 0.942 is −5.8 %. Dimensional check on a: √(J/(kg K) × K) = √(m²/s²) = m/s.
**Barometer (optional).** The two leading causes of a ±3 % disagreement are weather (sea-level pressure varies by roughly 950–1050 hPa, i.e. ±5 %) and the phone sensor's own offset (typically ±1 hPa, 0.1 %). Weather dominates.
**Common mistakes.** Using z instead of h and then wondering why the "exact" answer is 0.01 % off at high altitude; forgetting that L is negative (temperature must *fall*); using R = 8.314 (the universal constant, per mole) with mass density.

## 2. Isothermal layer and scale height

**Assumptions.** Hydrostatic equilibrium (used to write dp = −ρ g₀ dh); ideal gas (used to replace ρ by p/(R T)); constant T = T_b in the layer (used to pull 1/T out of the integral); g = g₀ constant, which is the geopotential convention.
**Derivation.** dp/p = −(g₀/(R T_b)) dh. Integrate from h_b to h: ∫_{p_b}^{p} dp′/p′ = −(g₀/(R T_b)) ∫_{h_b}^{h} dh′, so ln(p/p_b) = −(g₀/(R T_b))(h − h_b). Exponentiate: p = p_b exp(−(h − h_b)/H) with **H = R T_b/g₀**.
**Units.** H: (J/(kg K))(K)/(m/s²) = (m²/s²)/(m/s²) = m. A length, as required for the exponent to be dimensionless.
**Result.** H = 287.05287 × 216.65/9.80665 = **6 341.6 m**. p(20 km) = 22 632.0 × exp(−9000/6341.6) = 22 632.0 × 0.24191 = **5 474.9 Pa**. Code: `layer_bases()[1][2]` = 5 474.877 Pa; relative difference 5 × 10⁻⁶.
**Check.** Limiting case: as h → h_b, p → p_b. Each scale height divides pressure by e ≈ 2.718; 9000 m is 1.42 scale heights, and e^−1.42 ≈ 0.242. Alternative: the gradient formula (6) with L → 0 has the limit (1 + L Δh/T_b)^(−g₀/(R L)) → exp(−g₀ Δh/(R T_b)), the same result, which is why the two branches join continuously.
**Common mistakes.** Writing H = R T/g with the universal R (gives H in metres per mole, which is nonsense); using the sea-level T₀ for a layer that is at 216.65 K; forgetting the minus sign and getting pressure *rising*.

## 3. Mars

**Assumptions.** The idealized single gradient layer stated in the problem; CO₂ as an ideal gas with the given R and γ; hydrostatic equilibrium. Geopotential conversion skipped (declared), since the exercise is about the structure of the code; at 10 km on Mars (R_M = 3 389 500 m) the correction is 0.3 %, so this choice is visible in the fourth figure.
**Derivation.** Exponent −g/(R L) = 3.71/(188.92 × 0.0025) = 7.855. T(10 km) = 210 − 25 = 185.0 K. p = 610 × (185/210)^7.855 = 225.4 Pa. ρ = 225.4/(188.92 × 185.0) = 6.449 × 10⁻³ kg/m³. a = √(1.3 × 188.92 × 185.0) = 213.2 m/s.
**Result.**

| at 10 km | T (K) | p (Pa) | ρ (kg/m³) | a (m/s) |
|---|---|---|---|---|
| Earth (ISA) | 223.15 | 26 436 | 0.4127 | 299.5 |
| Mars (idealized) | 185.0 | 225.4 | 0.00645 | 213.2 |

Scale heights at the surface: Earth 287.05 × 288.15/9.80665 = **8 435 m**; Mars 188.92 × 210/3.71 = **10 694 m**. Mars's atmosphere is thinner in *amount* (surface density 0.0154 kg/m³, 1.3 % of Earth's) but *thicker* in extent, because low gravity lets it spread upward.
**Argument.** A parachute decelerates by drag, D = ½ ρ V² C_D S, and must roughly balance weight m g at terminal speed. On Mars g is 38 % of Earth's but ρ at the altitudes where parachutes open is about 1 % of Earth's, so for the same C_D S the terminal speed is √(0.38/0.01) ≈ 6 times higher. To reach the same terminal speed the area S must be about 40 times larger, which is why Mars landers use parachutes tens of metres across and still need retro-rockets.
**Check.** Dimensional: ρ = p/(R T) in kg/m³ on both planets. Limiting case: at h = 0 the code must return exactly 610 Pa and 210 K. The untouched functions were `_prepare`, `temperature`, `pressure`, `density`, `speed_of_sound`, `isa`, `isa_table`: every planet-specific number lived in a named constant, which is what Step 1 of the lesson was for. If you had to edit a formula, a constant was hiding in it.
**Common mistakes.** Keeping Earth's R_AIR (gives ρ 1.5 times too small); keeping H_TOP at 84 852 m so the range check lets through altitudes the single-layer model was never defined for; keeping γ = 1.4 for CO₂.

## 4. Density altitude (example: ΔT = 20 K)

**Assumptions.** ISA pressure at the airport elevation (the "ISA + ΔT" convention); dry air.
**Derivation.** At h = 1 655 m: T_ISA = 288.15 − 10.76 = 277.39 K; p_ISA = 82 960 Pa. Hot-day T = 297.39 K; ρ = 82 960/(287.05287 × 297.39) = **0.9718 kg/m³**. Invert the ISA density on a grid of 2 × 10⁶ points from 0 to 20 km (1 cm spacing): `np.interp(0.9718, rho_grid[::-1], grid[::-1])`.
**Result.** Density altitude = **2 347 m ± 0.01 m** (grid resolution; the physical uncertainty is set by how well you know ΔT, and 1 K in ΔT moves the result by about 33 m, which is the number that matters). Confirmation: `density(2347.43)` = 0.971800 kg/m³, equal to ρ to seven figures. For other seeds: ΔT = 11 K gives 2 043 m; 25 K gives 2 512 m; 30 K gives 2 672 m.
**Check.** Sanity: a 20 K excess at about 288 K is a 7 % density drop; near 1.6 km ISA density falls about 1 % per 100 m, so about 700 m of extra altitude. 1 655 + 700 ≈ 2 350 m. Consistent.
**Interpretation.** The aircraft takes off as if from a 2 350 m runway elevation: less lift at a given speed (so a higher ground speed is needed), less thrust from a normally aspirated engine, and typically a 15–20 % longer takeoff roll than at ISA. Pilots compute density altitude for exactly this reason.
**Common mistakes.** Passing a decreasing x array to `np.interp` (it returns garbage without complaint); forgetting to add ΔT to the ISA temperature *at the airport altitude* rather than at sea level; reporting the grid resolution as if it were the physical uncertainty.

## 5. Dynamic pressure and Mach tables

**Assumptions.** ISA; incompressible dynamic pressure ½ ρ V² (fine as a loading measure at any Mach; the *compressible* impact pressure differs above M ≈ 0.3, but q itself is defined this way).
**Derivation.** `q = 0.5 * density(hs)[:, None] * V[None, :]**2`; `M = V[None, :] / speed_of_sound(temperature(hs))[:, None]`.
**Result.** (a) and (b): both shape (21, 6). (c) q > 20 kPa in 20 cells: rows 0–2 (0–2 km) at columns 3–5 (200–300 m/s), rows 3–6 at columns 4–5, rows 7–9 at column 5 only. The region is the bottom-left corner of the table (low altitude, high speed) with a boundary sloping so that higher altitude needs higher speed: q is fixed when ρ V² is fixed. M > 0.8 in 35 cells: column 5 (300 m/s) at every altitude from 0 km (M = 0.88), column 4 (250 m/s) from 7 km upward (M = 0.847 at 11 km), and column 3 (200 m/s) nowhere below 20 km. The Mach region is the right-hand side and widens with altitude, because a falls from 340 to 295 m/s as temperature falls to the tropopause, then stays constant.
**Hand check.** Cell (h = 0, V = 300): q = 0.5 × 1.225 × 300² = 55 125 Pa; table gives 55 125.0. M = 300/340.29 = 0.8816; table 0.8816. Cell (20 km, 300 m/s): ρ = 0.08803, q = 3 961 Pa; table 3 961.6.
**Check.** Dimensional: (kg/m³)(m²/s²) = kg/(m s²) = Pa. Limiting: the q table at V = 0 would be all zeros; M is dimensionless and increases monotonically along each row.
**Common mistakes.** Using `density(hs) * V**2` without the `None` indices (shape error, or worse, silent broadcasting if the lengths happen to match); computing M with density instead of temperature; reading `argwhere` pairs as (speed, altitude).

## 6. When the table lies

**(a)** Measurable sign: the reported surface temperature differs from the ISA value at that elevation by more than about 10 K, or a radiosonde shows the tropopause more than 2 km from 11 km. Decision that changes: takeoff. A +25 K day at a 1 600 m airport is a density altitude near 2 500 m, and for a light piston aircraft the takeoff roll grows by a quarter or more; the go/no-go on a short runway flips. For a rocket, a cold day raises density and hence maximum dynamic pressure at the same trajectory, which can force a throttle-down to protect the structure.
**(b)** Convert: h = 6 356 766 × 100 000/6 456 766 = 98 451 m. Extrapolate layer 7: T = 214.65 − 0.002 × (98 451 − 71 000) = **159.7 K**; exponent −g₀/(R L) = 9.80665/(287.05287 × 0.002) = 17.08; p = 3.9564 × (159.7/214.65)^17.08 = **0.025 Pa**. The upper-atmosphere models give roughly 0.03 Pa and 190–200 K at 100 km. The pressure is within a factor of 1.3 of reality, which is *why this extrapolation is dangerous*: it produces a plausible number from a wrong temperature (35 K low, and falling where the real atmosphere starts warming at the mesopause near 86–90 km). A user who saw 10⁻⁹ Pa would know something was broken; a user who sees 0.025 Pa will believe it, and then compute a density and a speed of sound that are 20 % wrong, in a regime where the ideal-gas and fixed-composition assumptions are themselves failing. Models should refuse outside their domain, as ours does, rather than produce believable fiction.
**Common mistakes.** Forgetting the geopotential conversion (1.5 km error at 100 km); arguing from "it is not in the standard" rather than from the physics (temperature trend reverses; composition changes).

## Check-yourself answers

1. K/m in the code (−6.5 × 10⁻³); −6.5 K/km in the troposphere.
2. About 22 632/e ≈ 8 300 Pa (code: 8 326 Pa).
3. Because a = √(γ R T) depends only on temperature, and temperature rises through the stratosphere (ozone absorbs ultraviolet light). Density has nothing to do with it.
4. Everywhere. Both branches are computed for every element; `np.where` only chooses afterwards. Hence the `lapse_safe` trick.
5. (21, 6): shape (21, 1) broadcast against (1, 6).
6. Always smaller for z > 0, since R_E/(R_E + z) < 1; equal only at z = 0.
7. The test compares an exact derivative with a finite-difference estimate whose truncation error is of order dh² × p‴/p, about 10⁻⁸ relative here, plus rounding error from subtracting two nearly equal pressures. 10⁻⁵ is loose enough to pass a correct implementation on any machine and tight enough that any formula error, which produces residuals of order 1, fails.
8. Derived. The ISA defines T₀ = 288.15 K and p₀ = 101 325 Pa (and the gas constant through the composition); ρ₀ = p₀/(R T₀) = 1.225 kg/m³ follows from the ideal gas law.
