# P05.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical. Draw the free-body diagram before writing any equation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Material data (typical room-temperature values, not design allowables): 2024-T3 E = 73.1 GPa, F_ty = 345 MPa, F_tu = 483 MPa; 7075-T6 E = 71.7 GPa, F_ty = 503 MPa, F_tu = 572 MPa, ρ = 2810 kg/m³.

## 1. (Warm-up: a bracket arm on edge and on its side)

A control-system bracket arm is a solid 2024-T3 bar 20 mm wide and h = (30 + s) mm deep, built in at one end, with a downward 800 N limit load at the free end 0.25 m away.

(a) Draw the free-body diagram of the arm and find the bending moment at the support.
(b) Compute I, the section modulus and the largest bending stress at the support. Which edge is in tension?
(c) Compute the margins of safety on yield (factor 1.0) and on ultimate (factor 1.5).
(d) The arm is mounted on its side by mistake (20 mm deep, h wide). By what factor does the stress change, and what is the new yield margin?

Answer form: a moment; three numbers and an edge; two margins; a factor and a margin.

## 2. (Derivation: the flexure formula, the general formula and the two-flange limit)

(a) From plane sections, ε = −κy, and Hooke's law, show that zero net axial force puts the neutral axis through the centroid and that the moment resultant gives κ = M/(EI). State which step uses equilibrium, which compatibility and which the constitutive law.
(b) Let ε = a + bx + cy with x, y from the centroid. Impose N = 0, M = −∫σy dA and M_h = −∫σx dA and derive equation (6) of the lesson. Show that it reduces to σ = −My/I when I_xy = 0 and M_h = 0, and find the neutral-axis angle.
(c) An idealised spar has two caps of area A, thickness t and width b = A/t whose centroids are h_c apart, and no web, so the overall depth is h = h_c + t. Show that the ratio of the two-flange stress M/(A h_c) to the flexure stress at the outer fibre is (h_c² + t²/3)/(h_c h). Evaluate it for h_c = 170 mm and t = 10, 5 and 2 mm, and find its limit as t → 0.

Answer form: three derivations, a ratio with three values and a limit.

## 3. (Unfamiliar scenario: a launch-vehicle tank at maximum dynamic pressure)

A launch vehicle's propellant tank is a thin aluminium cylinder of radius R = 1.8 m and wall thickness t = (4 + s/10) mm. At maximum dynamic pressure a crosswind gives a bending moment of 2.0 MN·m at a station where the axial force is 1.5 MN in compression (thrust pushing against the inertia of the stages above). Ignore internal pressure.

(a) Show that a thin-walled tube has A ≈ 2πRt and I ≈ πR³t, and find the relative error of that I against the exact π/4((R + t/2)⁴ − (R − t/2)⁴).
(b) Find the axial stress from the axial force alone, the largest bending stress, and the largest compressive and tensile total stresses. Where around the circumference do they act?
(c) Below what bending moment is the whole wall in compression? Show that the answer is |N|R/2, independent of t.
(d) The compressive side is the critical one even though the tensile allowable and compressive yield strength of the alloy are similar. Why? (One sentence; the subject of `P05.15`.)

Answer form: two formulas and an error; four stresses and two locations; a moment and a proof; one sentence.

## 4. (Numerical, uses your code: size the spar of your own aeroplane)

Your aeroplane has a maximum mass of m = 550 + 12s kg, a limit load factor of 3.8 and a semi-span of 4.5 + s/20 m, with elliptic lift. The root spar is a 7075-T6 I-section of depth h = 0.14 + s/500 m with 10 mm thick caps and a 3 mm web.

(a) Compute the semi-span lift and the root bending moment.
(b) With `section_properties`, `i_section` and `scipy.optimize.brentq`, find the cap width b_f that makes MS_ultimate = 0 using the flexure formula at the outer fibre. Compare the cap area with the two-flange estimate from `cap_area_required`.
(c) What is the limit stress at that size, and the yield margin? Explain why the limit stress does not depend on your seed.
(d) Check your spar with `cantilever_fem` using 20 elements: compare E κ c at the first element midpoint with the flexure formula at the same station.

Answer form: two numbers; a width, two areas; a stress, a margin and a sentence; two stresses and their relative difference.

## 5. (Numerical, uses your code: a Z-stringer, free and riveted)

A wing upper-skin stringer is a 2024-T3 Z-section 30 mm deep with 15 mm flanges, 1.6 mm thick (`z_section(0.030, 0.015, 0.0016)`), carrying a bending moment M = (40 + s) N·m about the horizontal axis (positive: compression at the top).

(a) Compute I_xx, I_yy and I_xy, and the neutral-axis angle when the stringer is free to bend in any direction.
(b) Compute the stress at the two flange tips and the two web–flange corners, at mid-thickness: (±14.2, ±14.2) mm and (∓0.8, ±14.2) mm, by the general formula and by the symmetric formula.
(c) When the stringer is riveted to the skin, the skin prevents it from bending sideways, so its curvature is purely vertical (b = 0 in ε = a + bx + cy). Show that the stress is then exactly −My/I_xx and that the skin must supply a lateral moment M_h = M I_xy/I_xx. Evaluate M_h and confirm with `bending_stress_general` that this M_h gives the symmetric stress.
(d) Which analysis would you use for a stringer between two widely spaced ribs where it is riveted to a thin, flexible skin? Say what you would check.

Answer form: three properties and an angle; a table of eight stresses; a proof, a moment and a check; a short argument.

## 6. (Judgment: is the spar margin real?)

A colleague reports MS_ultimate = +0.39 for the lesson's spar root and wants to remove material from the caps. For each assumption below, say whether removing it raises or lowers the root cap stress, estimate the size of the effect with your code where you can, and rank them by how much they matter.

(a) The lift is elliptic. (Repeat with the same total lift distributed uniformly.)
(b) The wing's own weight is ignored. (Take 45 kg per wing, spread uniformly, at the same load factor.)
(c) The tail load is ignored. (Take a tail download equal to 10 % of the aircraft weight, carried by extra wing lift.)
(d) The spar caps are solid. (Take two 6 mm fastener holes across the 60 mm bottom cap; compute the net-section properties and stress.)
(e) The skins carry no bending.
(f) The compression cap cannot buckle.

Then answer: should material be removed? What test or analysis would you ask for first?

Answer form: six signed effects with numbers where computable, a ranking, a recommendation.
