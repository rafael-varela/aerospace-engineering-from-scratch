# P08.17  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. `isa(10668.0)`; `speed_of_sound(T)`; M = V/a. 2. T₀ = T × `T0_over_T(M)`; p₀ = p × `p0_over_p(M)`; q_c = p₀ − p; ½ρV². 3. `bernoulli_airspeed(qc, rho)`; error V_B/V − 1; `rho0_over_rho(M)`.

## 2.
(a) Energy per unit mass c_pT + V²/2 = c_pT₀; divide by c_pT; use a² = γRT. Then the isentropic power laws. T₀/T needs only adiabatic flow; p₀/p and ρ₀/ρ need isentropic. (b) Series in x = (γ − 1)M²/2 up to x³; subtract 1; multiply by p and use p = ρa²/γ; divide by ½ρV² = ½γpM²; airspeed ∝ √(q_c) gives half the relative change. (c) ρ₁u₁ = ρ₂u₂; p₁ + ρ₁u₁² = p₂ + ρ₂u₂²; h₁ + u₁²/2 = h₂ + u₂²/2. Solve for p₂/p₁ and M₂; Rayleigh: p₀₂/p₁ = (p₂/p₁)(p₀₂/p₂) with p₀₂/p₂ = `p0_over_p(M2)`; simplify.

## 3.
(a) T = T₀/`T0_over_T(0.8)`, p = p₀/`p0_over_p(0.8)`, ρ = p/(RT), V = M√(γRT), ρV. (b) √(2(p₀ − p)/ρ) and √(2(p₀ − p)/ρ₀) with ρ₀ = p₀/(RT₀); compare with V. (c) Solve for the density that makes Bernoulli exact; compare with ρ and ρ₀.

## 4.
1. `isa(h)`; for each M, `rayleigh_pitot(M)`, reading = ratio × p, `mach_from_pitot(ratio)`. 2. `mach_from_p0_over_p(ratio)`; `bernoulli_airspeed(reading - p, rho)`; errors against M and M a. 3. `normal_shock(M)["p02_p01"]`. 4. `stagnation_temperature(T, 3 a)`.

## 5.
cas(qc) = a_SL√(5((q_c/p_SL + 1)^{2/7} − 1)) with a_SL = `speed_of_sound(T_SL)`; q_c from the isentropic relation at the flight condition; EAS = V√(ρ/ρ_SL); convert with 0.514444 m/s per knot; then √(2q_c/ρ_SL).

## 6.
(a) `rho0_over_rho(M) - 1` and `bernoulli_error(M)` at the three Mach numbers. (b) Compute the Mach number of each case and compare its error with the accuracy needed. (c) Think about the probe's angle to the flow, its Reynolds number, and the time the flow has to come to rest in the probe.
