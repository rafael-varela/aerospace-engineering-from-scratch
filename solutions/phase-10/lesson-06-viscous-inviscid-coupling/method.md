# P10.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Re_x = Vx/ν. 2. `bl = falkner_skan(0)`: θ = I_θ √2 x/√Re_x, δ* = Hθ, c_f = 2f″(0)/√(2Re_x). 3. Thwaites λ = 0: θ = √(0.45νx/V), H = 2.61, c_f = 2νℓ(0)/(Vθ) with ℓ(0) = 0.22 (or call `thwaites` with a constant u_e). 4. Michel on the Blasius layer: scan Re_x with `michel_transition(Re_x, 0.664√Re_x)`. 5. Favourable gradients delay transition; adverse ones advance it.

## 2.
(a) Multiply the momentum integral by 2u_eθ/ν; use λ and c_f = 2νℓ/(u_eθ); write F(λ) = 2[ℓ − (2 + H)λ]; replace F by 0.45 − 6λ; the equation becomes d(θ²u_e⁶)/ds = 0.45νu_e⁵. (b) Integrate with u_e = Cs; λ = θ²C/ν. Compare with `falkner_skan(1)["I_theta"]`. (c) Set the ratio to 1; argue that the wake's momentum deficit falls as the wake velocity recovers.

## 3.
1. Re = 12 × 0.20/ν. 2. `section_drag("2412", radians(α), Re)` for α = 0, 4; read `x_tr`, `mechanism`, `cd`. 3. Compare with the physics of a laminar separation bubble; decide direction of error.

## 4.
1. `d = read_ladson(find_data_file())`; `Re, cd = ladson_zero_lift_drag(d, "I")`. 2. Loop: `section_drag("0012", 0, R, 100)["cd"]`; percentage differences. 3. Repeat for "VII" with `trip_x = 0.05`. 4. Plot both on a log Re axis, measured points with distinct markers and labels.

## 5.
1. Loop α; `section_drag("2412", radians(α), 4.52e6)`. 2. `displacement_iteration` for four α; percentage change of c_l. 3. D = c_d q S; P = DV for clean and tripped; difference.

## 6.
1. `section_drag("0012", 0, 1.19e7, 100, trip_x=x)` for x = 0.05, 0.1, 0.2; interpolate to the measured c_d. 2. Reasons: free-stream turbulence and acoustic disturbances, surface roughness, the empirical Michel correlation's spread; low-biased turbulent closure or Squire–Young at high Re. 3. Choose an amplification (e^N) transition method or a two-equation turbulent closure, and say what data (transition location, wake profiles) would validate it.
