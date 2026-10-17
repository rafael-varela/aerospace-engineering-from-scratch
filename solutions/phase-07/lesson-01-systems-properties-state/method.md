# P07.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T = 300 + 20s. 2. `cp_molar("N2", T)`; divide by `molar_mass("N2")` for per kilogram and by R_U for cp/R; γ = cp/(cp − R_U). 3. `h_molar`, `s0_molar`; s = s° − R_U ln 2; divide by M. 4. `specific_volume(2e5, T, M)`, ρ = 1/v. 5. Central difference of `h_molar` over ±10 K; compare with (a).

## 2.
(a) ∫ aₖ Tᵏ dT for k = −2, −1, 0, 1, 2, 3, 4 (the k = −1 term gives a logarithm); add b₁R; divide by RT. ∫ aₖ Tᵏ⁻¹ dT similarly for s°/R (k = 0 gives ln T), add b₂. (b) h/(RT) = 2.5 + b₁/T = 0 at T = 298.15. (c) h = u + RT; dh/dT = du/dT + R; both are total derivatives for an ideal gas; R > 0 so cp > cv. (d) Gibbs relation per species at partial pressure p_i; sum with mole-fraction weights; separate Σ x_i ln x_i; evaluate with the four mole fractions; divide by M_air.

## 3.
(a) `cp_molar("CO2", 210)` and at 300; γ = cp/(cp − R_U). (b) f = 2/(γ − 1); subtract 5. (c) Compare with 1.29; ratio of √γ for the speed of sound. (d) `cp_molar("H2", 300)`; compare the H₂ rotational and vibrational temperature scales qualitatively.

## 4.
(a) Two calls of `exit_temperature`. (b) q/(T₂ − T₁). (c) `air_state(p=2e6, h=float(air_h(800) + q))`. (d) Compare both exit temperatures with 1650 K.

## 5.
(a) s₁ = `air_s(288.15, 101325)`; `air_state(p=pi*101325, s=float(s1))`. (b) 288.15 × π^(0.4/1.4). (c) `air_h(T2) − air_h(T1)` and 1005 × (T₂,const − T₁). (d) State postulate: p and s are independent; T and h are not.

## 6.
For each: estimate how close to condensation (Z), whether vibration or dissociation is active, and whether the composition is air; name the replacement model (lesson 03 tables, lesson 10 products, lesson 11 equilibrium, a real-gas equation of state).
