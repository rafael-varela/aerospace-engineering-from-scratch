# P10.23  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, give uncertainties as standard uncertainties unless a problem asks for k = 2, and label every number that comes from the synthetic dataset **SYNTHETIC**. Room air: 20 °C, 101.3 kPa, ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s unless you measured your own.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: how small are the numbers?)

The lesson's model has chord 100 mm and span 200 mm. For jet speeds of 3, 4 and 5 m/s compute:

(a) the dynamic pressure q, the water-column height it supports, and the reading along a manometer tube inclined at 10°;
(b) the Reynolds number, and its ratio to Ladson's lowest, 2.0 × 10⁶;
(c) the lift in grams at CL = 0.5 and the drag in grams at CD = 0.05;
(d) the number of resolution steps those forces span on a 1 g kitchen scale and a 0.01 g jewellery scale (the latter through a bell-crank of ratio 2). Which quantity cannot be measured at the minimum-cost tier?

Answer form: a table of nine numbers; three numbers; six numbers; a short table and a sentence.

## 2. (Derivation: the instruments)

(a) Derive the inclined-manometer relation q = (ρ_w − ρ)g Δl sin θ from the hydrostatics of a U-tube with one leg inclined, assuming the reservoir level does not move. Show that the relative uncertainty of q from the incline alone is (cot θ) u_θ, and evaluate it for θ = 10° measured to ±0.2°. Why is a shallow incline both helpful and harmful?
(b) Draw the bell-crank: drag D at height l_D above the pivot, the scale under the horizontal arm at distance l_s. Derive the scale's reading change ΔF = D l_D/l_s by moments, and say why calibrating with known horizontal loads is better than measuring l_D and l_s.
(c) Derive the uncertainty of C = F/(qS) when F can be near zero, and explain why the relative form (u_C/C)² = Σ(u_i/x_i)² fails there.

Answer form: three derivations with a number and two sentences.

## 3. (Unfamiliar scenario: where does a small drone wing sit?)

The reference small UAV (Phase 10 conventions) has a 0.20 m chord and cruises at 12 m/s at sea level.

(a) Compute its Reynolds number and place it between the bench rig (about 2.6 × 10⁴) and Ladson's tunnel.
(b) With `laminar_upper_surface`, find the laminar separation point and the Michel transition point at α = 0°, 2° and 4° for the UAV. Does the boundary layer transition before it separates?
(c) Which published dataset in the lesson's references bridges the gap between the bench and Ladson, and in what Reynolds-number range?
(d) What would you expect of the UAV's section cd and lift slope compared with the bench and with Ladson, and why does it matter for choosing an airfoil (`P10.14`)?

Answer form: a number and a sentence; a table and a yes/no; a citation and a range; a paragraph.

## 4. (Numerical, uses your code: a synthetic dataset of your own)

Generate the SYNTHETIC dataset with seed 200 + s (`write_synthetic_csv(folder, seed=200 + s)`), reduce it with `reduce_bench`, and label everything SYNTHETIC.

(a) Report q, V, Re, the arm ratio R, the sting tares, the lift slope over −4° to 4°, the zero-lift angle and CL_max with their uncertainties.
(b) Compare each with the generator's stated truth (`make_synthetic_bench_data(...)["truth"]` and its docstring). Which recovered values are within 2u of the truth, and which are not? Explain any that are not.
(c) Compute E_n of your lift slope against the inviscid finite-span prediction band (AR_eff 4–6) and against Ladson's slope at 2.0 × 10⁶. What does each say?

Answer form: nine numbers with uncertainties; a table of verdicts with reasons; two E_n values and two sentences.

## 5. (Numerical, uses your code: the uncertainty budget)

For the lesson's SYNTHETIC dataset (seed 230):

(a) At α = 0°, 4° and 10°, split the relative uncertainty of CL into the force term and the q term, and split the q term into the manometer, the blockage band and the fan's unsteadiness.
(b) Which single instrument upgrade would most reduce u(CL) at 4°? And u(CD) at 0°?
(c) Estimate u(CL) and u(CD) at 4° if the lift and drag had been read on a 1 g kitchen scale (minimum-cost tier, bell-crank ratio 2, two runs). Is CD still measurable?

Answer form: a table; two upgrades with reasons; two uncertainties and a sentence.

## 6. (Judgment: is it the Reynolds number or the rig?)

The bench lift slope is about half of Ladson's and the drag about five times larger.

(a) List four rig effects (not Reynolds-number effects) that could lower the measured lift slope, and the sign of each.
(b) Design a test that separates the Reynolds-number effect from the rig effects using only the equipment of the recommended tier (hint: what can you change without changing the rig?). What result would point to the Reynolds number?
(c) The blockage estimate A_frontal/(4A_jet) comes from closed test sections. Why is the sign of the correction uncertain in an open jet, and why does the lesson carry it as an uncertainty rather than apply it?
(d) Write the two sentences of your report's conclusion that a reviewer could not fault, if your measured slope were 0.057 ± 0.001 per degree.

Answer form: four effects with signs; a test design and a criterion; a short argument; two sentences.
