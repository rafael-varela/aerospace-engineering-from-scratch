# P10.23  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. q = ½ρV², h = q/((ρ_w − ρ)g), Δl = h/sin 10°. 2. `reynolds_number(V, 0.1)`. 3. L = CL q S, D = CD q S, in grams. 4. Divide by the resolution (times R = 2 for the drag scale).

## 2.
(a) Pressure balance across the liquid: Δp = (ρ_w − ρ)g × vertical height = (ρ_w − ρ)g Δl sin θ; differentiate with respect to θ. (b) Moments about the pivot with the scale force; calibration includes friction, compliance and the true lever arms. (c) First-order propagation with absolute terms.

## 3.
1. Re for c = 0.20 m, V = 12 m/s, ISA sea level. 2. `laminar_upper_surface("0012", α, Re)` for three α. 3. The Eppler 387 report in the references. 4. Argue from laminar separation bubbles.

## 4.
1. `write_synthetic_csv(tmp, seed=200 + s)`; `load_bench_dataset(tmp, seed=200 + s)`; `reduce_bench`. 2. Compare with the truth dict and docstring (the flow angle appears as −truth). 3. `en_number` with the predicted band's midpoint and its half-width/√3; Ladson slope from a fit to Table I at 2.0 × 10⁶.

## 5.
1. Rebuild u_q and u_L from the reduction's pieces at each α. 2. Compare the largest terms. 3. Resolution term √2 r/√12 per reading pair, divided by √n and by R for drag.

## 6.
(a) Blockage, jet non-uniformity and swirl, end-plate leakage, flow angularity, balance friction, sting interference. (b) Two or three fan speeds; plot slope against Re; reverse the model. (c) Open-jet boundary condition; sign. (d) Numbers, scope, and the word SYNTHETIC if applicable.
