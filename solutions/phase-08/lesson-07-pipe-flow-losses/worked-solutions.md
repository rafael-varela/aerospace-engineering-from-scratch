# P08.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Properties, roughnesses and loss coefficients are typical values; Colebrook, Haaland and Blasius are correlations of measured data.

## 1. A cabin air duct (s = 15: V = 8 m/s)

**Assumptions.** Incompressible (the drop is a few hundred pascals against 101 kPa), fully developed, isothermal air at 20 °C, ε = 0.15 mm.

(a) Re = 1.204 × 8 × 0.10/1.81 × 10⁻⁵ = **53 215**: turbulent (well above 4000).

(b) ε/D = 1.5 × 10⁻³. Colebrook by Newton: **f = 0.02515**, in **3** iterations from the Haaland start; Haaland: **f = 0.02490**, **−0.99 %**.

(c) Friction drop f(L/D)ρV²/2 = 0.02515 × 80 × 38.53 = **77.50 Pa**; Q = 8 × π × 0.1²/4 = **0.0628 m³/s**; fan power for the friction alone Q dp = **4.87 W**. With K = 4 × 0.3 + 0.5 + 1.0 = 2.7 the fittings add 104.03 Pa: total **181.5 Pa** (fan power 11.4 W before the fan's own efficiency).

(d) Smooth pipe at the same Re: f = 0.02060; fully rough limit: 0.02173. The duct's f is above the fully rough value at this ε/D. This is not a contradiction: Colebrook adds the smooth and rough terms inside one logarithm, so in the transitionally rough range f exceeds both limits. The duct is in that range, and roughness raises f by 22 % over a smooth duct.

**Independent check.** Haaland, an explicit formula fitted independently, agrees to 1 %.

**Common mistakes.** Using ε in millimetres with D in metres; adding minor losses with the wrong velocity.

## 2. Derivations

(a) With V = 4Q/(πD²) and f = 64μ/(ρVD): dp = (64μ/(ρVD))(L/D)(ρV²/2) = 32μLV/D² = 32μL(4Q/(πD²))/D² = **128μLQ/(πD⁴)** ✓, the Hagen–Poiseuille law. The density cancels: laminar pipe flow does not depend on inertia at all, which is why roughness does not matter either.

(b) g(x) = x + 2 log₁₀(a + bx), a = ε/(3.7D), b = 2.51/Re > 0. **g′(x) = 1 + 2b/((a + bx) ln 10) > 0** and **g″(x) = −2b²/((a + bx)² ln 10) < 0**: increasing and concave. Newton: **x_{k+1} = x_k − g(x_k)/g′(x_k)**. For a concave function the tangent lies above the curve, so the tangent's zero lies at or to the right of the root: after one step every iterate is at or above x*, and from there the iterates decrease monotonically to x* (each tangent zero is between x* and the previous point). Hence convergence from any x₀ > 0. The fixed point x_{k+1} = φ(x_k) = −2 log₁₀(a + bx_k) converges linearly with rate **|φ′(x*)| = 2b/((a + bx*) ln 10)**, which is small (about 0.1) because b is small.

(c) V = 4ṁ/(ρπD²) at fixed ṁ, so ρV²/2 ∝ D⁻⁴ and dp ∝ f D⁻⁵. Fully rough: f constant (at fixed ε/D), **dp ∝ D⁻⁵**. Blasius: Re ∝ VD ∝ D⁻¹, f ∝ Re^{−1/4} ∝ D^{1/4}, **dp ∝ D^{−4.75}**. With `segment_loss` at 10, 20, 40 mm: smooth slopes −4.78 and −4.77 (Blasius's −4.75, with the smooth Colebrook curve a little steeper than its fit); at fixed ε/D = 10⁻³, −4.89 and −4.84, approaching −5 as the flow becomes fully rough.

**Common mistakes.** Forgetting that a fixed ε with a changing D changes ε/D; expecting the exponent to be exactly −5 at moderate Re.

## 3. The oxidiser side

**Assumptions.** As in the lesson; LOX properties typical near 90 K; no cavitation (the tank is at 2.5 MPa, far above the vapour pressure of about one atmosphere near the normal boiling point).

(a) V = **6.570 m/s**, Re = **961 133**, f = **0.01288** (smooth would be 0.01172). Re is nineteen times the fuel line's because μ_LOX is ten times smaller than μ_RP-1 while ρV D is larger; at that Re the flow is nearly fully rough even for drawn tubing, and f depends mostly on ε/D.

(b) Injector area 174.0 mm²; friction 19.03 kPa; minor (K = 5.20) 128.07 kPa; injector 400.0 kPa; head 20.14 kPa: **tank 2527.0 kPa**. The fittings dominate the line loss by a factor of nearly seven: the dynamic pressure ρV²/2 = 24.6 kPa is multiplied by K = 5.2.

(c) Fuel tank (lesson) 2513.6 kPa. With both tank pressures held and p_c raised to 2.1 MPa: ṁ_O = 3.3261 kg/s, ṁ_F = 1.4422 kg/s, **O/F = 2.306** (from 2.300). Both flows fall by about 10 %, but almost in proportion, because both feed systems are dominated by quadratic losses (orifice and fittings): the ratio moves by 0.3 %. The mixture ratio is robust; the thrust is not.

**Common mistakes.** Treating the LOX line as laminar because the flow "looks slow"; forgetting to keep the injector area fixed when the flow changes.

## 4. Sizing your own feed line (s = 15: ṁ = 1.95 kg/s)

(a) **D_min = 20.04 mm** for 10 % (200 kPa); **23.61 mm** for 5 % (100 kPa).

(b)

| tube (OD × wall, mm) | ID (mm) | V (m/s) | Re | line loss (kPa) | % of p_c |
|---|---|---|---|---|---|
| 15.875 × 0.889 | 14.097 | 15.42 | 88 062 | 897.9 | 44.9 |
| 19.05 × 0.889 | 17.272 | 10.28 | 71 874 | 376.2 | 18.8 |
| 22.225 × 1.245 | 19.735 | 7.87 | 62 904 | 213.4 | 10.7 |
| 25.4 × 1.245 | 22.910 | 5.84 | 54 186 | 113.5 | 5.7 |

The 22.225 mm tube misses the 10 % limit (10.7 %); choose **25.4 × 1.245 mm**, which gives 5.7 %, a margin of 86.5 kPa against the 200 kPa allowance.

(c) A larger line weighs more, traps more propellant (unusable residual mass and a longer chill-down for cryogens) and responds more slowly at start-up and shut-down; its lower velocity also changes the line's dynamic coupling with the combustion chamber. The aim is the smallest line that meets the loss limit with margin.

**Common mistakes.** Using the outside diameter; forgetting that the minor losses scale with D⁻⁴ too.

## 5. Newton against fixed point, and a laminar line

(a) a = 10⁻⁴/3.7, b = 2.51 × 10⁻⁵. Fixed point: **10 iterations** from Haaland, **11** from x = 7, with rate |φ′(x*)| = 0.103 per iteration. Newton: **4** from either start. Newton's errors from x = 7: 0.3494, 7.133 × 10⁻⁴, 2.82 × 10⁻⁹, then zero to round-off; e_{k+1}/e_k² = 0.0058 and 0.0055: constant, i.e. quadratic convergence.

(b) ṁ = 0.005 kg/s: V = 0.1761 m/s, **Re = 1179** (laminar), f = 0.05429, **dp = 211.36 Pa**; Hagen–Poiseuille 128μLQ/(πD⁴) = **211.36 Pa** ✓ (identical, as 2(a) proves). At 0.01 kg/s, Re = 2358 (f = 0.02825, dp = 439.9 Pa), and at 0.02 kg/s, Re = 4716 (f = 0.03830, dp = 2386 Pa). The 0.01 kg/s result is the one not to trust: Re = 2358 is in the transition range, where the code interpolates between 64/Re and Colebrook. The real flow may be laminar or intermittently turbulent depending on disturbances, and f can be anywhere between the laminar 0.0271 and the turbulent (Colebrook) 0.0471 at that Re. A line designed to run there has an unpredictable drop.

**Common mistakes.** Counting the initial guess as an iteration; trusting the third significant figure of a transitional friction factor.

## 6. What controls the tank pressure?

Baseline tank pressure 2513.6 kPa (friction 42.18 kPa, minor 95.27 kPa, injector 400 kPa).

| uncertainty | change in p_tank (kPa) |
|---|---|
| μ × 0.7 / × 1.3 | −3.0 / +2.5 |
| ε × 10 / × 30 | +4.4 / +11.9 |
| every K × 0.5 / × 1.5 | −47.6 / +47.6 |
| Cd 0.665 / 0.735 (area fixed) | +43.2 / −37.2 |

Ranking: **fitting loss coefficients ≈ injector discharge coefficient ≫ roughness > viscosity**. Viscosity barely matters because f varies slowly with Re in turbulent flow (as Re^{−0.2} or less), and roughness only matters once it is large. Test time should go into a water-flow test of the assembled line and injector, which measures the K's and Cd directly. Refining the friction factor is wasted effort when the fittings carry more than twice the friction loss.

**Common mistakes.** Ranking by the size of the input uncertainty (±30 % on μ sounds large) rather than by its effect on the output.

## Check-yourself answers

1. f = 8τ_w/(ρV²), the wall shear stress made dimensionless with the dynamic pressure; dp = f(L/D)ρV²/2.
2. 64/Re follows from the exact Hagen–Poiseuille solution of the Navier–Stokes equations. Turbulent flow has no exact solution: Colebrook is a fit to measurements (with Prandtl's and von Kármán's log-law forms in its two limits).
3. ε/D is the typical height of the wall roughness relative to the diameter. When the roughness elements stick out through the viscous sublayer (large Re), the drag is form drag on the roughness and f depends on ε/D alone: the fully rough regime.
4. Haaland is within about 1.5 % of Colebrook, so the start is close, and Newton then doubles the number of correct digits each step; the function is increasing and concave, so the iteration cannot overshoot into trouble.
5. dp ∝ D^{−4.75} to D⁻⁵: a 20 % smaller line has 0.8^{−4.75} to 0.8⁻⁵ ≈ 2.9 to 3.1 times the loss.
6. p_tank = p_chamber + Δp_injector + Δp_friction + Δp_minor − ρ a H: the hydrostatic term is subtracted because the propellant column above the injector, under the vehicle's acceleration a, pushes the propellant along.
7. Between about 2300 and 4000 the flow is transitional: turbulent puffs and slugs come and go depending on inlet disturbances, so f is not a single function of Re.
8. The line is short (L/D = 125) and has many fittings (K = 5.95), so K ρV²/2 exceeds f(L/D)ρV²/2 = 2.6 ρV²/2. Fewer fittings, smoother bends or a rounded tank outlet save more than a smoother tube would.
