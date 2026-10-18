# P10.02  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Coordinates are in chords unless a chord is given. Say which results are exact (closed form) and which are numerical, and with how many points.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one station by hand)

Your section is the NACA 44TT with TT = 10 + (s mod 9) (for s = 15, NACA 4416).

(a) At x = 0.25, compute by hand the half-thickness y_t, the camber y_c, the camber slope dy_c/dx and its angle θ, and the upper and lower surface points (x_u, y_u) and (x_l, y_l). Note that x_u ≠ x_l.
(b) Compute the leading-edge radius and the trailing-edge thickness of the open (standard) form.
(c) Check (a) against `naca_coordinates` by evaluating the formulas in the code at the same station (the generated points are not at x = 0.25, so do not interpolate; call the building blocks).

Answer form: eight numbers; two numbers; a check to 10⁻¹².

## 2. (Derivation: the constants behind the families)

(a) Integrate 2y_t from 0 to 1 and show that the area of the 4-digit thickness form is 0.68508 t c² (open) and 0.68088 t c² (closed). Explain why the 0.685 often quoted is a property of the thickness form only, and estimate how much camber changes the area of a 2412.
(b) Show that near the leading edge y_t² ≈ 2rx with r = (5 × 0.2969 t)²/2 = 1.1019 t².
(c) Show that the closed-trailing-edge coefficient −0.1036 makes the five thickness coefficients sum to zero, and compute the open trailing-edge thickness.
(d) For the 5-digit camber line, show that the camber is greatest at x = p exactly when r³ − 3r² + 6pr − 3p² = 0, and find r for p = 0.15. Then explain why the ideal lift coefficient of the line is proportional to k₁, and use thin-airfoil theory (c_l,ideal = πA₁, A₁ = (2/π)∫₀^π dy_c/dx cos θ dθ, quoted here, derived in `P10.04`) to set up the equation that fixes k₁.

Answer form: four derivations; r to four figures.

## 3. (Unfamiliar scenario: a battery pocket in a foam UAV wing)

The reference small UAV has a rectangular wing of chord 0.20 m and span 1.5 m, cut from extruded polystyrene foam (take 30 kg/m³).

(a) For a NACA 2412 and a NACA 4412 section, compute the cross-section area, the volume of foam in the wing and its mass.
(b) A 20 mm thick flat battery pack must sit inside the wing with at least its full thickness available at every chordwise station it covers. Over what chordwise range is the local vertical thickness at least 20 mm? Does the camber change the answer?
(c) Compute the volume between 20 % and 60 % chord over the full span: how much foam does a pocket of that extent remove, as a fraction of the total?

Answer form: three numbers per section; a range; two numbers and a sentence.

## 4. (Numerical, uses your code: does spacing matter?)

(a) For the NACA 2412, compute the polygon area with n = 10, 20, …, 320 panels per surface on cosine spacing and on uniform spacing (build the uniform version from the same formulas), against a reference with 20 000 panels per surface. Report the observed orders.
(b) Explain the order 1.5 on uniform spacing from the behaviour of y_t near x = 0. Show by a substitution why cosine spacing restores order 2.
(c) A panel method (`P10.05`) on this section needs the leading edge resolved. With uniform spacing and n = 100, what is the x-coordinate of the first point after the nose, compared with the leading-edge radius? With cosine spacing?

Answer form: a table and two orders; a short argument; two comparisons.

## 5. (Numerical, uses your code: the 5-digit constants)

(a) Use `naca5_constants` to derive (r, k₁) for the 210, 220, 230, 240 and 250 camber lines (design c_l 0.3). Compare with the commonly quoted table values (0.0580, 361.4), (0.1260, 51.64), (0.2025, 15.957), (0.2900, 6.643), (0.3910, 3.230) and compute the ideal c_l that each quoted pair actually gives.
(b) Which quoted pair is furthest from its design condition? By how much, and is the discrepancy in r or in k₁?
(c) What are r and k₁ for a "430" line (design c_l 0.6, p = 0.15)? Explain why r does not change.

Answer form: a table; a sentence with a number; two numbers and a sentence.

## 6. (Judgment: what does "12 % thick" mean?)

(a) Measure the maximum vertical thickness (upper minus lower at the same x) and the maximum of the mean line (y_u + y_l)/2 of the generated NACA 0012, 2412, 6412 and 9415 sections, and their chordwise positions. Compare with the designation.
(b) The 23012's camber line has its maximum at x = 0.15 by construction. Where does the mean line of the generated surface peak, and why are they different?
(c) A colleague builds a wind-tunnel model from coordinates generated with thickness added *vertically* to the camber line, not perpendicular to it. For which sections, and where on the section, does that matter? Would a panel method notice?

Answer form: a table; a number and a sentence; a short argument.
