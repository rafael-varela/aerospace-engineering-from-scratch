# P04.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `wire_resistance(2.0, 0.644e-3)`; per metre /2; drop IR; power I²R. (b) R₂ = (10 + s) kΩ; `voltage_divider`; I = 12/(R₁ + R₂); P = 12 I; R_th = `parallel([R1, R2])`.

## 2.
(a) Rows: KCL at nodes 1, 2, 3 and the source constraint V₁ = V_s; g's of R₁…R_g; the source current I_s in rows 1 and 4. (b) KCL at 2 and 3 with I_g = (V₂ − V₃)/R_g; I_g = 0 ⇒ V₂ = V₃ = V_s R₂/(R₁ + R₂) = V_s R_x/(R₃ + R_x) ⇒ R₁/R₂ = R₃/R_x; conversely, if the ratios are equal, V₂ = V₃ is a solution with I_g = 0, and the solution is unique. (c) Loops: (V_s, R₁, R₂), (R₁, R_g, R₃), (R₂, R_x, R_g); write the 3 × 3 symmetric system; I_g is the difference of the two loop currents sharing R_g. (d) ∫₀^∞ (V_s/R)² e^{−2t/RC} R dt; or: charge Q = CV_s passes through the source, which does work QV_s = CV_s²; the capacitor stores ½CV_s²; the rest is dissipated.

## 3.
(a) C = 4πε₀R; Q = CV; N = |Q|/e. (b) |E| = k|Q|/r² at r = R and r = 2 m; V(2 m) = kQ/2; `electric_field([Q], [[0,0,0]], [2,0,0])` and `electric_potential`. (c) ½QV; compare the energy with a microcircuit's tolerance and the time scale with its response.

## 4.
(a) P = (100 + s) + 100; I = P/28; `wire_diameter_for_drop(I, 8, 28, 0.03)` with each resistivity; mass = density × 2L × πd²/4. (b) `bus_circuit(28, wire_resistance(8, d), loads)`; `solve_dc`; V_load = V₂ − V₃; current = `source_currents[0]`; loss 2I²R_w. (c) Fixed-point iteration V ← 28 − 2R_w P/V from V = 28. (d) Mass, diameter, and what drives the choice: connectors, bend radius, creep and galvanic corrosion with aluminium, temperature.

## 5.
(a) Build the divider; `voltmeter_reading(div, 2, 0, R_m)` for 10 MΩ, 1 MΩ and infinity. (b) V_read = V_true R_m/(R_m + R_th); solve for V_true; R_th = 0.5 MΩ. (c) `ammeter_reading(c, 0, R_a)`. (d) The error is R_th/R_m to first order: state the ratio needed for a given accuracy.

## 6.
(a) τ = 10 × 1000 × 10⁻⁶; I₀ = 28/10; t₉₉ = τ ln 100. (b) ½CV² each; `rc_energy_accounting(28, 10, 1e-3, 25 τ)`. (c) `solve_ivp(rc_rhs, …)` against `rc_capacitor_voltage`. (d) 0.1 Ω: I₀ = 280 A, τ = 0.1 ms; the energy lost is unchanged; the peak current is what trips protection and stresses switches and connectors.

## 7.
For each, name the broken assumption and its observable: temperature coefficient of copper; meter resistance comparable to the circuit's; unmodelled series resistance; negative incremental resistance of constant-power loads; a singular matrix with no reference node.
