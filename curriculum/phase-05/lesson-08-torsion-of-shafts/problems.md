# P05.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Material values are the typical values of `MATERIALS` in `starter/main.py`, and shear strengths come from `shear_allowables`.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a tail-rotor drive tube)

A helicopter tail-rotor drive shaft is a 2024-T3 tube of outside diameter 50 mm and wall 1.5 mm, 2 m between couplings, transmitting 100 + 5s kW at 6000 rpm.

(a) Find the torque, the shear stress at the outer and inner surfaces, and the twist over 2 m.
(b) Compare the outer stress with the shear yield estimate Fty/√3. What is the ratio?
(c) Treat the tube as thin and compute the shear stress from Bredt–Batho, τ = T/(2πr_m²t). How far is it from (a), and which of the two numbers would you report?

Answer form: four numbers with units, a ratio, a number and a sentence.

## 2. (Derivation: circular shafts, thin closed sections and open sections)

(a) Derive τ = Tr/J and φ = TL/GJ, naming where compatibility, Hooke's law and equilibrium each enter.
(b) Derive q = T/(2A_m) for a single-cell thin-walled section from equilibrium of a wall element and of the whole section. Then derive the twist rate (4) by equating the work of the torque to the strain energy in the walls.
(c) For a thin circular tube, show that J_Bredt/J_exact = 1/(1 + (t/2r_m)²) exactly.
(d) Show that a circular tube slit along its length is softer in torsion than the closed tube by the factor 3(r_m/t)².

Answer form: four derivations.

## 3. (Unfamiliar scenario: the D-spar of a rotor blade)

The leading edge of a helicopter rotor blade is often a closed D-shaped spar: a semicircular nose of mid-line radius 40 mm (wall 1.2 mm) closed by a flat vertical web (wall 2.0 mm). Take it as 7075-T6 and apply a torque of 400 + 10s N·m.

(a) With `bredt_batho`, using a polygon of at least 400 segments for the nose, find A_m, the shear flow, the shear stress in the nose and in the web, J and the twist rate in degrees per metre.
(b) Check A_m and J against the exact semicircle (A_m = πR²/2, ∮ds/t = πR/1.2 + 2R/2.0).
(c) If the web were omitted (an open semicircular section plus a separate flat web, connected only at one edge), what would J be? Give the ratio of closed to open stiffness.

Answer form: six numbers with units, two checks, one number and a ratio.

## 4. (Numerical, uses your code: your own turboshaft)

Redo the lesson's sizing for an engine of 1000 + 20s kW at 20 900 rpm, with the same factors, length and twist requirement, in Ti-6Al-4V.

(a) For d_i/d_o = 0 and 0.8, tabulate the diameters required by each criterion, the governing one, the mass per metre and the limit shear stress.
(b) Find the twist requirement (in degrees at limit torque over 600 mm) at which the ultimate-strength criterion would take over from the twist criterion, for both ratios.
(c) What does (b) tell you about how the shaft should be specified?

Answer form: a 2-row table, two angles, two sentences.

## 5. (Numerical, uses your code: Prandtl FD as the check)

(a) Compute J of a solid 10 × 1 strip with `prandtl_fd` at h = 0.1, 0.05 and 0.025, and compare with the series (6) and with bt³/3. Report the observed order.
(b) Compute J of a hollow 2 × 1 rectangle with a wall of 0.1 on all sides at h = t/2, t/4, t/8 and t/16. Report the observed order and explain why it is not 2.
(c) Compare the finest result of (b) with Bredt–Batho on the mid-line rectangle. Is the difference consistent with the lesson's square-tube results?

Answer form: three numbers, two references and an order; four numbers, an order and a sentence; a ratio and a sentence.

## 6. (Judgment: torsion in practice)

Answer each in two or three sentences, with a number where one is asked for.

(a) Why does the lesson's design case put a factor of 1.25 on rated torque before applying the ultimate factor of 1.5, and what would justify a different number for a real engine?
(b) A designer proposes a 0.9 diameter-ratio titanium shaft to save 31 % of the mass. Name two checks that the lesson did not make and that this shaft must pass.
(c) An access hole must be cut in the upper skin of the wing box between two ribs. Using 2.5, explain what happens to the torsional stiffness of that bay and what structural measure restores it.
(d) Why are splined and keyed shaft ends the usual sites of fatigue cracks, even when the shaft body is lightly stressed?

Answer form: four short paragraphs.
