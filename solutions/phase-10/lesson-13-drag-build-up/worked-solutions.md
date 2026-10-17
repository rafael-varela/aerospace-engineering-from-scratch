# P10.13  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results come from `solution.py` and a probe script that calls it, given to the figures shown. Values from the Ladson file are authentic measurements (NASA TM-4074); geometry, interference factors and add-ons are stated assumptions.

## 1. Friction and transition (s = 15: V = 55 m/s)

(a) Re = 55 × 1.2/1.461 × 10⁻⁵ = **4.52 × 10⁶**. C_f = **0.003422** (x_tr = 0), **0.003099** (0.1), **0.002496** (0.3), **0.001936** (0.5), **0.000625** (1.0, all laminar).

(b) Laminar to mid-chord saves **43 %** of the fully turbulent friction. A sailplane's drag at its best glide is about half friction, so a long laminar run is worth a large fraction of its performance; that is why sailplane wings are built to tight waviness tolerances and cleaned of insects before competition.

**Independent check.** The saving grows linearly with x_tr only roughly: from 0 to 0.1 it is 9.4 %, from 0.1 to 0.3 a further 17.6 %; the laminar friction itself rises as √x_tr.

**Common mistakes.** Using Re based on x_tr in the turbulent term; forgetting that C_f is per side.

## 2. Derivations

(a) Per unit width the drag of one side is qC_fL. A fully turbulent plate gives qC_f,t(Re)L. Over the laminar length x_tr, that turbulent model would have contributed qC_f,t(Re_tr)x_tr; the laminar layer actually contributes qC_f,l(Re_tr)x_tr. Replacing one by the other and dividing by qL: **C_f = C_f,t(Re) − (x_tr/L)[C_f,t(Re_tr) − C_f,l(Re_tr)]**. The approximation: the turbulent layer behind the transition point is treated as if it had grown from the leading edge (its "virtual origin" is at x = 0), whereas in reality it starts with the laminar layer's momentum thickness at x_tr, i.e. from a virtual origin somewhat behind the leading edge.

(b) L/D = C_L/(C_D0 + kC_L²). Setting d/dC_L = 0: (C_D0 + kC_L²) − 2kC_L² = 0 ⇒ **kC_L² = C_D0**: the induced drag equals the zero-lift drag. Then **C_L* = √(C_D0/k)** and **(L/D)max = C_L*/(2C_D0) = 1/(2√(kC_D0))**.

(c) **V* = √(2(W/S)/(ρC_L*))** with C_L* ∝ C_D0^{1/2}: V* ∝ (W/S)^{1/2}C_D0^{−1/4}. Retracting the gear lowers C_D0 (to 0.0187 in Problem 4) and so lowers C_L* and raises V* (by about (0.0318/0.0187)^{1/4} ≈ 1.14, partly offset by the change of e).

**Common mistakes.** Maximising C_L/C_D with respect to V instead of C_L (same answer, more algebra); forgetting that e itself depends on C_D0 in (4).

## 3. A sailplane

**Assumptions.** As stated; c = 11/15 = 0.733 m, AR = 20.45; exposed wing area 11 − 0.6 × 0.733 m²; the wing's wetted area uses the NACA 2412 perimeter ratio as an approximation for a 14 % section.

(a) At 30 m/s:

| component | Re | C_f | FF | Q | S_wet (m²) | D/q (m²) |
|---|---|---|---|---|---|---|
| wing (laminar to 0.5c) | 1.51 × 10⁶ | 0.00255 | 1.303 | 1.00 | 21.56 | 0.0717 |
| fuselage | 1.33 × 10⁷ | 0.00287 | 1.074 | 1.00 | 7.92 | 0.0244 |
| tails | 1.23 × 10⁶ | 0.00430 | 1.206 | 1.05 | 4.46 | 0.0243 |

Sum of D/q = 0.1205 m²; on 11 m² with 5 % miscellaneous, **C_D0 = 0.0105**.

(b) e_LL(AR 20.45, taper 0.4) times the fuselage factor (1 − 2(0.6/15)²), then (4): **e = 0.770**. **(L/D)max = 34.3** at **C_L* = 0.722**, **V* = 30.1 m/s**, minimum-glide sink rate **V*/(L/D) = 0.88 m/s**.

(c) The landing gear: the trainer's fixed gear is 41 % of its C_D0; the sailplane has none in flight. Then the long laminar run on a clean wing, and the high aspect ratio for the induced drag.

**Independent check.** A best glide ratio in the mid-thirties at about 30 m/s is the right order for a club sailplane of this span; the build-up is internally consistent (at C_L*, C_Di = C_D0 = 0.0105).

**Common mistakes.** Leaving the 12 m² of the trainer in the denominator; using the trainer's AR in e.

## 4. The price of fixed gear

(a)

| variant | C_D0 | e | (L/D)max |
|---|---|---|---|
| retractable (D/q = 0) | 0.0187 | 0.779 | 16.5 |
| baseline fairings (0.15 m²) | 0.0318 | 0.707 | 12.1 |
| no fairings (0.25 m²) | 0.0406 | 0.666 | 10.4 |

(b) At the baseline's 46.13 kW of drag power: retractable **66.7 m/s**, no fairings **50.2 m/s** (against 55 m/s).

(c) Retracting the gear would add 21 % to the cruise speed for the same power, but a trainer flies short legs, does many landings and is flown by students who may forget to lower the gear; 20 kg and the maintenance of a retraction system cost more than the speed is worth. Wheel fairings, on the other hand, are cheap and are worth nearly 5 m/s (or about 20 % less fuel per hour at the same speed). The numbers rest on assumed gear drag areas, which are the least certain inputs.

**Common mistakes.** Keeping e fixed while C_D0 changes.

## 5. Drag with lift (authentic data)

(a) and (b), Ladson Table I:

| Re | α (°) | cl | cd | cd at α ≈ 0 | x_tr at 4° | x_tr at 0° |
|---|---|---|---|---|---|---|
| 2.00 × 10⁶ | 4.25 | 0.445 | 0.00670 | 0.00620 | 0.398 | 0.464 |
| 3.94 × 10⁶ | 3.96 | 0.420 | 0.00631 | 0.00659 | 0.317 | 0.279 |
| 5.97 × 10⁶ | 4.18 | 0.452 | 0.00760 | 0.00650 | 0.072 | 0.217 |
| 8.86 × 10⁶ | 4.13 | 0.450 | 0.00700 | 0.00630 | 0.084 | 0.179 |
| 1.19 × 10⁷ | 4.17 | 0.458 | 0.00687 | 0.00634 | 0.056 | 0.128 |
| 1.89 × 10⁷ | 4.20 | 0.464 | 0.00716 | 0.00655 | NaN | 0.032 |

With lift the effective (average of both surfaces) transition moves forward at most Reynolds numbers: the upper surface's suction peak is followed by a stronger adverse gradient, which destabilises the laminar layer. The 3.94 × 10⁶ row goes the other way, and the README flags free-transition drag points at this Re (cd at 10° below cd at 8°) as laminar-bucket behaviour: the lower surface's laminar run lengthens with α while the upper's shortens, and the balance varies.

(c) At 1.89 × 10⁷ and 4.2° the measured cd (0.00716) exceeds even the fully turbulent zero-lift estimate (0.00681): no transition location can reproduce it, because the drag increase with lift (thicker, pressure-gradient-loaded boundary layer, pressure drag) is not in a zero-lift form factor.

(d) The zero-lift estimate is a good C_D0 ingredient only if the transition location is right; at cruise lift coefficients the section drag at α ≈ 4° differs from its zero-lift value by −4 % to +17 % in these data (+8 % to +17 % in five of the six blocks), which is exactly the lift-dependent profile drag represented in the lesson by the K term of (4).

**Common mistakes.** Treating the effective transition as a measured transition point: it is the location that makes a flat-plate model match, averaged over both surfaces.

## 6. How good is C_D0 = 0.032?

(a) Ranked by effect on C_D0 over a plausible range:
1. gear drag area ±0.05 m²: **±0.0044**;
2. cooling drag area ±0.02 m²: **±0.0018**;
3. interference factors ±0.1 on wing and fuselage: **±0.0013**;
4. wing transition 0–0.3c: **0.0304–0.0326** (about ±0.0011);
5. miscellaneous allowance ±3 % of the base: **±0.0009**.

(b) Treating them as independent and adding in quadrature: **C_D0 = 0.032 ± 0.005** (about ±16 %).

(c) The section estimate's error is a friction-model error on the wing (and, less, the tails and fuselage), which make up 48 % of C_D0; the build-up already chose a transition location (0.1c) instead of fully turbulent flow, and the free-transition data at the trainer's Re suggest 0.2–0.3c, a change of 0.0014 in C_D0. The 59 % applies to a fully turbulent model of a smooth 2-D section at Re = 2 × 10⁶, not to the trainer's estimate.

(d) Power-off glides (engine at idle, propeller drag measured or feathered) at a range of steady airspeeds at known weight, recording sink rate: C_D = (W/(qS))(V_sink/V). Plot C_D against C_L²: the intercept is C_D0 and the slope is 1/(πeAR).

**Common mistakes.** Adding uncertainties linearly and then calling the result "conservative" without saying so.

## Check-yourself answers

1. Laminar C_f = 1.328/√Re, turbulent 0.455/(log₁₀Re)^{2.58}; at Re = 5 × 10⁵: 0.00188 and 0.00511, a factor of 2.7.
2. The extra drag of a thick body over a flat plate of the same wetted area: the air speeds up over the thickness (more friction) and the pressure recovery adds pressure drag; both grow with thickness.
3. C_D0 = Σ C_f,i FF_i Q_i S_wet,i/S + (add-on drag areas)/S, times (1 + miscellaneous): friction coefficient, form factor, interference factor, wetted area, reference area.
4. They are equal.
5. As Re rises the friction coefficient falls, but transition moves forward (from 46 % to 3 % of the chord in the inversion), adding turbulent friction; the two effects almost cancel.
6. The grit adds its own pressure drag and produces a thicker, rougher turbulent boundary layer than a smooth plate's; a smooth-plate C_f contains neither.
7. The landing gear's drag area (41 % of C_D0 here): it is a bluff, partly separated body whose drag is pressure drag, not skin friction, so it must come from measurements of similar installations.
