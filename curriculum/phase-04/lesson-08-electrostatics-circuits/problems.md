# P04.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which numerical, and (if you did the optional bench work) which measured.

Several problems use your personal seed **s**: the day of the month you were born (1–31). ρ(copper) = 1.68 × 10⁻⁸ Ω m, ρ(aluminium) = 2.65 × 10⁻⁸ Ω m, densities 8960 and 2700 kg/m³.

## 1. (Warm-up: a sensor lead and a divider)

(a) A sensor lead is 2.0 m of 0.644 mm diameter copper (AWG 22) and carries 0.5 A. Compute its resistance, its resistance per metre, the voltage drop and the power dissipated.
(b) A divider of R₁ = 10 kΩ and R₂ = (10 + s) kΩ sits across 12 V. Compute the output voltage, the current, the total power and the Thévenin resistance seen at the output.

Answer form: eight numbers with units.

## 2. (Derivation: nodal analysis, bridge balance, the half-energy rule)

(a) For the Wheatstone bridge of the lesson (V_s from node 1 to ground; R₁: 1–2; R₂: 2–0; R₃: 1–3; R_x: 3–0; R_g: 2–3), write the MNA system with unknowns (V₁, V₂, V₃, I_s): the 4 × 4 matrix and the right-hand side. Show the matrix is symmetric.
(b) Show from the node equations that the galvanometer current is zero for every R_g if and only if R₁/R₂ = R₃/R_x.
(c) Write the same circuit with three loop currents and KVL, and show that the galvanometer current is the difference of two loop currents.
(d) Prove that charging a capacitor C from 0 to V_s through any resistance R dissipates ½CV_s², by integrating I²R. Then prove it without any integral, from the charge that passes through the source.

Answer form: a matrix and a vector, two proofs, three equations, two proofs.

## 3. (Unfamiliar scenario: a charged spacecraft)

In geostationary orbit during eclipse, a spacecraft can charge negatively relative to the surrounding plasma. Model it as an isolated conducting sphere of radius 1.0 m at a potential V = −(1 + s/10) kV.

(a) Compute its capacitance 4πε₀R, the charge, and the number of excess electrons.
(b) Compute the field at its surface and at 2 m from its centre, and the potential at 2 m. Check the field against `electric_field` and the potential against `electric_potential` for a point charge at the centre.
(c) Compute the stored electrostatic energy ½QV. A discharge (an arc) releases a fraction of it in nanoseconds. Why is that dangerous for electronics even though the energy is small?

Answer form: three numbers, three numbers with two code checks, one number and two sentences.

## 4. (Numerical, uses your code: sizing the bus harness)

The bus is 28 V; the loads are (100 + s) W, 60 W and 40 W at 28 V; the harness is 8 m each way.

(a) With `wire_diameter_for_drop`, find the copper and aluminium diameters for a 3 % drop at the nominal current, and the harness mass of each (supply plus return).
(b) Build the circuit with `bus_circuit` at the copper diameter and solve it. Report the load voltage, drop, current, load power and harness loss. Why is the drop below 3 %?
(c) Real loads behind DC–DC converters draw constant power. Iterate V_load = 28 − 2R_w P/V_load to convergence and report the drop and current. Which model is the conservative one for sizing?
(d) Which metal would you choose, and what else would you check before committing?

Answer form: four numbers, five numbers and a sentence, two numbers and a sentence, two sentences.

## 5. (Numerical, uses your code: trust your multimeter)

(a) A 1 MΩ/1 MΩ divider sits on 3.3 V. Compute what voltmeters of 10 MΩ and 1 MΩ read at the midpoint with `voltmeter_reading`, and the true value.
(b) Derive the correction V_true = V_read(1 + R_th/R_m) and apply it to the 10 MΩ reading.
(c) 5 V drives 100 Ω. Compute the current read by ammeters of 0, 1 and 10 Ω.
(d) Write a rule of thumb, with a number, for when a reading can be used without correction.

Answer form: three numbers, a derivation and one number, three numbers, one sentence.

## 6. (Numerical: inrush into a capacitor bank)

A 1000 μF capacitor bank is switched onto a 28 V bus through a 10 Ω inrush-limiting resistor.

(a) Compute τ, the peak current, and the time to reach 99 % of 28 V.
(b) Compute the energy stored and the energy dissipated in the resistor; check with `rc_energy_accounting`.
(c) Verify the closed form against `solve_ivp` with `rc_rhs` over 5τ.
(d) Without the limiter, the only resistance is 0.1 Ω of wiring. What are the peak current and τ, and what changes in the energy lost? Why do designers fit inrush limiters anyway?

Answer form: three numbers, two numbers with a check, one error, two numbers and two sentences.

## 7. (Judgment: limits of the DC model)

For each, say which assumption fails, what you would observe, and what to do. Two sentences each.

(a) The bus harness in direct sunlight on the outside of the spacecraft.
(b) Measuring a 10 MΩ sensor resistance network with a 10 MΩ multimeter.
(c) A crimp connector that has corroded.
(d) Ten constant-power loads switched on at once on a long, thin harness.
(e) A solver that returns node voltages for a circuit you forgot to ground.

Answer form: five short paragraphs.
