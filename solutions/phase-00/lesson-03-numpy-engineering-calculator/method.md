# P00.03  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. h = R_E z/(R_E + z). 2. T = T₀ + L h with L = −0.0065 K/m. 3. p = p₀ (T/T₀)^5.2559. 4. ρ = p/(R T). 5. a = √(γ R T). 6. Call `isa(h)`; compute (hand − code)/code for each. 7. Expect agreement at the fourth figure, limited by the rounded exponent and by your calculator's rounding, not by physics.

## 2.
1. With T = T_b constant, dp/p = −(g₀/(R T_b)) dh has a constant coefficient. 2. Integrate both sides from h_b to h: ln(p/p_b) = −(g₀/(R T_b))(h − h_b). 3. Exponentiate; define H = R T_b/g₀. 4. Dimensional check on H. 5. Evaluate H at 216.65 K. 6. p(20 km) = p(11 km) exp(−9000/H). 7. Compare with the code's third layer base.

## 3.
1. Copy the module; set T0 = 210, P0 = 610, G0 = 3.71, R_AIR = 188.92, GAMMA = 1.3, H_BASE = [0], LAPSE = [−2.5e-3]; set H_TOP to 40 000 and R_EARTH to Mars's 3 389 500 m (or skip the geopotential step for this exercise and say so). 2. Run `isa(10000.0)`. 3. Scale heights: R T₀/g for each planet. 4. For the parachute: drag force is ½ ρ V² C_D S; set it against weight m g on each planet and compare the ρ that is available. 5. Note which functions were untouched.

## 4.
1. ΔT from your seed. 2. ρ_actual = pressure(1655)/(R (temperature(1655) + ΔT)). 3. grid = np.linspace(0, 20000, N) with N large (10⁶ gives 2 cm resolution); rho_grid = density(grid). 4. h_d = np.interp(ρ_actual, rho_grid[::-1], grid[::-1]). 5. Check density(h_d) ≈ ρ_actual. 6. Interpret: the aircraft "sees" altitude h_d; takeoff distance grows roughly with 1/ρ² for a fixed-pitch piston aircraft, somewhat less for a jet.

## 5.
1. `rho = density(hs)`; `q = 0.5 * rho[:, None] * V[None, :]**2`. 2. `a = speed_of_sound(temperature(hs))`; `M = V[None, :] / a[:, None]`. 3. `np.argwhere(q > 20000)`, `np.argwhere(M > 0.8)`; read the index pairs as (row = altitude, column = speed). 4. The q region is low altitude and high speed (dense air, fast); the M region is high speed at any altitude, growing at higher altitude because a falls with T until 11 km. 5. Hand check: pick a cell, compute ρ and a from the lesson's equations.

## 6.
(a) Pick a measured quantity (surface temperature, radiosonde temperature at altitude, measured sea-level pressure) and a threshold; tie it to a performance decision. (b) Convert 100 km geometric to geopotential; apply equation (6) from the 71 km base with L = −0.002; compare with the stated real range; argue about plausibility and the behaviour of the temperature trend above the mesopause.
