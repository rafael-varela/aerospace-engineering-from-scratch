# P07.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, say which correlation you used and check that you are inside its stated range. Evaluate properties at the film temperature unless told otherwise.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Properties: air from Lienhard Table A.6 and water from Table A.3 (your `air_properties` and `water_properties`).

## 1. (Warm-up: a heated probe)

An anti-iced air-data probe is a cylinder of diameter D = 6 + s/10 mm held across a stream of air at T∞ = 260 K and U = 20 + s m/s; its heater keeps the surface at 330 K. Treat the air as at 1 atm.

(a) Compute the film temperature, the properties, Re_D and Pr.
(b) Compute Nu_D by Churchill–Bernstein and by the simplified low-Re form (7.66). Which applies, and why do they differ?
(c) Compute h and the heater power per metre of probe.
(d) By what factor does the heater power change if the airspeed doubles? Compare with the factor √2 a pure Re^(1/2) law would give, and explain the difference.

Answer form: four numbers; two Nusselt numbers and a sentence; h and a power per metre; a factor and a sentence.

## 2. (Derivation: the flat plate)

(a) From the local laminar law Nu_x = 0.332 Re_x^(1/2) Pr^(1/3), derive the average Nu_L and show it is twice the local value at x = L.
(b) Show that St_x Pr^(2/3) = c_f/2 for the laminar plate with c_f = 0.664 Re_x^(−1/2). Then show that for Lienhard's turbulent power laws (Nu_x = 0.0296 Re_x^0.8 Pr^0.6, c_f = 0.0592 Re_x^(−0.2)) the ratio St Pr^(2/3)/(c_f/2) is Pr^(0.6 − 1/3), and evaluate it for air and for water at 300 K.
(c) Derive the mixed-flow average Nu_L = 0.037 Pr^0.6 (Re_L^0.8 − Re_c^0.8) + 0.664 Re_c^(1/2) Pr^(1/3) for a sharp transition at Re_c, and prove that it is continuous at Re_L = Re_c.
(d) For a plate tripped to turbulence at its leading edge (Re_c → 0), what is the ratio of the average to the trailing-edge local Nu? Why is it smaller than 2?

Answer form: four derivations and two numbers.

## 3. (Unfamiliar scenario: an unpressurised equipment bay)

The lesson's 60 W black box (0.30 × 0.12 × 0.20 m, ε = 0.90, bottom on an adiabatic tray) is moved to an unpressurised bay of an aircraft cruising at 12 km, where the ISA pressure is 19 330 Pa (`P04.07`). For an ideal gas μ, k and c_p do not depend on pressure, but ρ ∝ p.

(a) Show that ν and α scale as 1/p, that Ra scales as p², and that a laminar natural-convection h scales roughly as p^(1/2). Evaluate p^(1/2) for this bay.
(b) Write a wrapper around `air_properties` that scales ρ, ν and α by the pressure ratio, and find the box's surface temperature with bay air and walls at 25 °C, at 1 atm and at 19 330 Pa.
(c) Repeat with bay air and walls at 250 K. What share of the heat does radiation carry in each of the four cases?
(d) Which assumption of the radiation model is most doubtful in a real bay, and in which direction would it move the answer?

Answer form: a derivation and a number; two temperatures; two temperatures and four percentages; a sentence.

## 4. (Numerical, uses your code: designing the cold plate)

A cold plate takes Q = 300 + 10s W into water in a single channel of length 1.2 m, inlet 20 °C. The hottest wall must stay below 40 °C.

(a) For D = 3, 4 and 5 mm, find the minimum flow that meets the wall limit (search upward from laminar flow; the laminar branch never meets it) and its Re.
(b) For each D, find the flow that puts Re at 5000, and the hottest wall there.
(c) Compute the pressure drop Δp = f (L/D) ρU²/2 (Filonenko's f) and the pumping power ṁΔp/ρ at the flows of (a) and (b).
(d) Choose D and ṁ and defend the choice in three sentences, using Re, wall margin and pumping power.

Answer form: a table of six flows with Re, wall temperatures, pressure drops and powers; a choice with reasons.

## 5. (The Reynolds–Colburn analogy in a wind tunnel)

A smooth flat plate 0.40 m long (in the flow direction) and 0.30 m wide sits edge-on in a wind tunnel at 10 m/s, air at 300 K. A drag balance measures the skin-friction drag on both faces as 0.0383 N.

(a) Compute Re_L and say whether the boundary layer should be laminar.
(b) Compute the measured average C_F and compare it with the laminar value 1.328 Re_L^(−1/2).
(c) Use the analogy (in its average form, St̄ Pr^(2/3) = C_F/2) to predict the average h, and compare with `plate_average_nu`.
(d) If the plate were tripped to turbulence at its leading edge, what drag and what h would you expect (C_F = 0.074 Re_L^(−0.2), Nu_L = 0.037 Re_L^0.8 Pr^0.6)? Could the drag balance tell you whether the plate was tripped?

Answer form: a number and a sentence; two coefficients and a ratio; two values of h; two numbers and a sentence.

## 6. (Judgment: review a thermal analysis)

A colleague's analysis of a water-cooled box contains these steps. For each, say what is wrong, estimate its effect with your code, and say what you would do instead.

(a) "Water at Re = 1800 in a 4 mm tube, uniform heat flux: Gnielinski gives Nu = 9.5, so h = 1420 W/(m² K)."
(b) "The 1.0 m tall vertical side of a hot cabinet in air has Ra = 3 × 10⁹; with (8.13a), Nu = 121."
(c) "The exposed bottom face of the lesson's box at 50 °C (25 °C air) has the same h as the top face, 5.8 W/(m² K)."
(d) "Radiation is small at these temperatures and is neglected." (for the lesson's black box at 60 W)
(e) "The correlation is accurate to 1 %, so the predicted 51.8 °C is good to ±0.1 K."

Answer form: five short paragraphs, each with a number.
