# P07.14  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (closed forms, rules) and which are numerical. Draw every enclosure and label its surfaces before writing an equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Surface data: Lienhard Tables 10.1 and 10.5 (in the lesson's constants). σ = 5.670374419 × 10⁻⁸ W/(m² K⁴).

## 1. (Warm-up: a sunshade tube)

An instrument looks out through a cylindrical sunshade of radius r = 0.05 m and depth h = 0.1 + s/100 m. Its three surfaces are the bottom disk (1), the side wall (2) and the open mouth (3), treated as a disk.

(a) Use `vf_coaxial_disks` to find F₁₃.
(b) Use the summation rule and reciprocity, and nothing else, to find F₁₂, F₂₁, F₂₃ and F₂₂.
(c) Check that every row of the 3 × 3 matrix sums to one and that A_i F_ij = A_j F_ji.
(d) Why is F₂₂ not zero, while F₁₁ is?

Answer form: one view factor from a formula, four from the rules, a check, a sentence.

## 2. (Derivation: the radiosity method and its limits)

(a) Starting from the view-factor integral, prove reciprocity in one line, and prove that Σ_j Q_j = 0 for any enclosure solved by the radiosity method.
(b) Derive the two radiosity equations (known T, known Q) from J = εE_b + (1 − ε)G and the definition of Q.
(c) For two surfaces, eliminate J₁ and J₂ to obtain the three-resistance formula. Reduce it to (i) large parallel plates and (ii) a small convex body in a large enclosure (A₁/A₂ → 0). What does each limit say about the role of the enclosure's emittance?
(d) n thin shields of emittance ε_s on both faces are placed between two large plates of emittances ε₁ and ε₂. Show that the flux is σ(T₁⁴ − T₂⁴)/[(1/ε₁ + 1/ε₂ − 1) + n(2/ε_s − 1)].

Answer form: four derivations.

## 3. (Unfamiliar scenario: a thermocouple that reads low)

An engine test cell measures the temperature of the exhaust air in a duct with a bare thermocouple bead. The gas is at T_g = 500 + 5s K, the duct walls are 150 K cooler, the bead has ε = 0.8, and the convection coefficient from the gas to the bead is 100 W/(m² K) (from a cylinder or sphere correlation of `P07.13`). The bead is small compared with the duct.

(a) Write the bead's energy balance and find its temperature with `brentq`. What is the measurement error?
(b) Repeat with (i) a bead coated to ε = 0.2 and (ii) twice the convection coefficient (a faster gas or a smaller bead).
(c) A thin cylindrical radiation shield of polished metal (ε = 0.1) surrounds the bead; gas flows over both its faces with h = 80 W/(m² K), and it sees the wall on its outside. Neglecting the bead's small effect on the shield, find the shield's and then the bead's temperature and the new error.
(d) Which of the three fixes would you choose, and why might each one fail in practice?

Answer form: a temperature and an error; two errors; two temperatures and an error; a short argument.

## 4. (Numerical, uses your code: a tilted solar array)

The lesson's 1.0 × 0.6 m radiator (its 1.0 m edge along the hinge) has a 1.0 m wide array hinged on that edge, but the array is opened to an angle φ = 100 + 2s degrees from the radiator's plane (φ = 90° is the lesson's perpendicular case; larger φ opens the pair up).

(a) Set up the two parallelograms for `vf_numerical` and compute F (radiator → array) for n = 10, 20, 40. Check reciprocity at each n.
(b) No closed form is given for this angle. Use the order you measured in the lesson for a shared edge to extrapolate (Richardson) to n → ∞, and check the same extrapolation on the 90° case against `vf_perpendicular_rectangles`.
(c) Build the three-surface enclosure (radiator Q = 200 W, ε = 0.93; array face at 330 K, ε = 0.83; space black at 3 K) and find the radiator temperature. Compare with the 90° case.

Answer form: a convergence table; an extrapolated value and a check; two temperatures.

## 5. (Numerical, uses your code: why a hole looks black)

A cubic cavity of side 1 m has five walls at 400 K with emittance ε_w = 0.1 + s/50 and one open face (a black surface at 0 K, i.e. the outside).

(a) Use `box_view_factors` and `solve_enclosure` to find the heat leaving through the opening, and the cavity's effective emittance ε_eff = Q_out/(A_open σT⁴).
(b) Repeat for ε_w = 0.1 and 0.9, and for a deeper 1 × 1 × 3 m cavity open at one end.
(c) Explain the trend, and say how this effect is used in a black-body calibration source and how it affects a spacecraft's deep, open-ended bays.

Answer form: two numbers; a table; a paragraph.

## 6. (Judgment: review a radiator design)

A colleague's radiator analysis contains these steps. For each, say what is wrong, compute its effect with your code, and say what you would do instead.

(a) "We linearised radiation about the sink: q = 4εσT_sink³(T − T_sink)."
(b) "We used α = ε = 0.93 for sunlight on the white paint, because Kirchhoff's law says α = ε."
(c) "The radiator is sized at 1.079 m² for normal sun with the paint's beginning-of-life α_s = 0.15." (What happens at end of life if α_s has risen to 0.30?)
(d) "The deployed array is neglected; it is only 0.3 of the view."
(e) "Twenty layers of aluminised film reduce the loss to 0.54 W/m², so the bay needs no heater."

Answer form: five short paragraphs, each with a number.
