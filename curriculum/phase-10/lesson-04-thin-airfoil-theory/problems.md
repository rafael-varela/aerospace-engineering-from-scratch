# P10.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Angles in radians inside formulas; quote degrees in answers. x in chords from the leading edge; cm about the quarter chord unless stated; flap deflection δ positive trailing-edge down.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a parabolic camber line by hand)

A thin airfoil has the camber line z = 4h x(1 − x) with h = 0.01 + s/1000.

(a) Show that dz/dx = 4h cos θ under the Glauert substitution, and write down A₀, A₁, A₂ at α = 3° without integrating anything else.
(b) Find α_L0, c_l at 3°, c_m,c/4, the ideal lift coefficient (the c_l at which A₀ = 0) and the centre of pressure at 3°.
(c) Check every number with `thin_airfoil`.

Answer form: three coefficients; five numbers; a check.

## 2. (Derivation: from the vortex sheet to the Fourier coefficients)

(a) Starting from the flow-tangency condition for a vortex sheet γ(x) on the chord, (1/2π)∫₀^c γ(ξ) dξ/(x − ξ) = V(α − dz/dx), substitute x = (c/2)(1 − cos θ) and show that γ(θ) = 2V[A₀(1 + cos θ)/sin θ + Σ Aₙ sin nθ] satisfies it, using the Glauert integral ∫₀^π cos nφ dφ/(cos φ − cos θ) = π sin nθ/sin θ. Derive the formulas for A₀ and Aₙ.
(b) Integrate γ for the lift and the moment and derive c_l = π(2A₀ + A₁) and c_m,c/4 = (π/4)(A₂ − A₁). Show that c_m,c/4 does not depend on α.
(c) Derive α_L0 = −(1/π)∫₀^π (dz/dx)(cos θ − 1) dθ.
(d) For a plain flap of chord ratio E on a flat plate, derive τ = 1 − (θ_h − sin θ_h)/π and dc_m,c/4/dδ = −½ sin θ_h(1 − cos θ_h), with cos θ_h = 1 − 2(1 − E). Show that τ → (4/π)√E for a small tab.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a tailless glider needs a reflexed section)

A flying wing has no tail to balance a nose-down pitching moment, so its section needs c_m,c/4 ≈ 0. Start from the NACA 2412 camber line and add a reflexed trailing-edge tab of chord ratio E = 0.15 (a plain flap deflected trailing-edge *up*).

(a) What tab deflection makes c_m,c/4 = 0?
(b) What do the tab and the 2412 camber give for α_L0 together? At what α does the reflexed section reach c_l = 0.5, compared with the plain 2412?
(c) What does the designer pay for a zero pitching moment, and why do purpose-designed reflexed sections usually put the reflex in the camber line rather than a tab?

Answer form: an angle; two angles and a comparison; two sentences.

## 4. (Numerical, uses your code: camber position)

(a) Tabulate α_L0 and c_m,c/4 of the NACA 4-digit camber lines with m = 0.02, 0.04, 0.06 and p = 0.2, 0.3, 0.4, 0.5, 0.6.
(b) Show from the table that both are proportional to m, and describe how they depend on p.
(c) Compare the 230 camber line (r = 0.2025, k₁ = 15.957) with the 2412's. Why did NACA move the camber forward in the 5-digit series?

Answer form: a table; two sentences; two numbers per section and a sentence.

## 5. (Numerical, uses your code: the discrete-vortex check)

(a) For the NACA 2412 camber line at α = 4°, compute the error of the lumped-vortex c_l against `thin_airfoil` for n = 10, 20, …, 320 on uniform and on cosine panel spacing, and the observed orders.
(b) Explain why uniform spacing gives a clean second order while cosine spacing gives erratic errors, though both converge. (Look at where the panel edges fall relative to the kink at x = 0.4.)
(c) For the parabolic camber line of Problem 1 the lumped-vortex c_l is exact for every n and every spacing, but c_m,c/4 is not. Report the c_m error for n = 10, 20, 40, 80 (uniform) and its order.

Answer form: a table and orders; a short argument; four numbers and an order.

## 6. (Judgment: how far can you trust a thin-airfoil flap?)

Thin-airfoil theory gives the 25 % plain flap on the 2412 a lift increment Δc_l = 2πτδ, so 0.334 at δ = 5°, 1.336 at 20° and 2.004 at 30°.

(a) List the three assumptions of the theory that a 30° plain flap violates, and say which one fails first.
(b) In `P10.05` the panel method gives the 12 % thick 2412 a lift 11 % above thin-airfoil theory at 4°. Does that mean a thick section's flap is also more effective than τ predicts? What else changes for a thick section with a flap?
(c) For the trainer's approach (c_l = 1.112) the lesson finds δ = 6.7° at α = 4°. Would you design the real flap for that deflection? What would you do instead?

Answer form: three assumptions and a sentence; a short argument; a decision with a reason.
