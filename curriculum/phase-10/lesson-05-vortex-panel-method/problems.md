# P10.05  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. γ and Γ are positive clockwise; the surface is listed in Selig order (counter-clockwise); V = 1 and chord-scaled geometry unless stated. Report every convergence result with its observed order.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one panel by hand)

A single panel runs from (0, 0) to (1, 0) and carries a uniform clockwise vortex sheet γ = 1 (that is, g_a = g_b = 1 in the linear form).

(a) Using I₁ = β (the angle the panel subtends) and I₂ = ln(r₁/r₂), compute by hand the velocity (u, v) it induces at (0.5, 0.5), (0.5, −0.5) and (2, 0).
(b) Compare each with a point vortex of strength Γ = 1 at the panel midpoint. Where is the difference largest, and why?
(c) Check (a) with `linear_vortex_influence` (add the a and b parts). At (0.5, 0.5), split the result into the parts from g_a and g_b and explain why their u parts are equal but their v parts are opposite.

Answer form: three velocity pairs; three comparisons and a sentence; a check and a sentence.

## 2. (Derivation: the influence integrals and the system)

(a) In panel coordinates (ξ along the panel, η to its left), derive I₃ = ∫₀^L s η/D ds = η ln(r₂/r₁) + ξβ and I₄ = ∫₀^L s(ξ − s)/D ds = −L + ηβ + ξ ln(r₁/r₂), with D = (ξ − s)² + η².
(b) Show how the velocity of a linear sheet γ(s) = g_a(1 − s/L) + g_b s/L follows from I₁–I₄, and why a point exactly on the panel needs the side to be specified (β → ±π).
(c) Count unknowns and equations for N panels and show why the Kutta condition is needed. Explain why γ₀ + γ_N = 0 expresses "equal speeds leaving the trailing edge" in this sign convention.
(d) Show that with the fluid inside the body at rest, the tangential velocity just outside the sheet equals the sheet strength (up to sign), and hence that c_l = 2Γ/(Vc) with Γ = Σ L_j(γ_j + γ_{j+1})/2.

Answer form: four derivations.

## 3. (Unfamiliar scenario: a "Kutta condition" on a circle)

A circle has no sharp edge, so potential flow does not fix its circulation. Make a polygon of N panels on the unit circle starting and ending at θ = 0 (the rearmost point), and run `vortex_panel` at α = 0° and 5° for N = 40, 80 and 160.

(a) What circulation does the method produce? Compare with 4πVa sin α.
(b) Explain the result with `P10.03`: what did the condition γ₀ + γ_N = 0 impose on the circle flow?
(c) What does this tell you about the Kutta condition's role on a real airfoil, and about a panel method's lift on a body with a rounded trailing edge?

Answer form: a table; a sentence; two sentences.

## 4. (Numerical, uses your code: your own exact benchmark)

Take the Karman–Trefftz section with b = 1, μ = (−0.05 − s/1000) + (0.03 + s/1000)i and τ = 5° + s/2 degrees, at α = 3°.

(a) Run `kt_convergence` for N = 20, 40, 80, 160, 320. Tabulate the errors in c_l (from Γ and from Cp), in c_m,c/4, and the largest Cp error on 5–95 % of the chord, with observed orders.
(b) How many panels does c_l need for 0.1 % accuracy? For 0.01 %?
(c) Why are the orders irregular at small N and closer to 2 at large N?

Answer form: a table with orders; two panel counts; a short argument.

## 5. (Numerical, uses your code: three NACA sections)

For the NACA 0012, 2412 and 4412 with 160 panels:

(a) Find the lift slope (as a multiple of 2π), α_L0, c_m,c/4 at α = 4°, and the minimum Cp at 4° with its chordwise position.
(b) Compare α_L0 and c_m,c/4 with thin-airfoil theory (`thin_airfoil_naca4`). Which differs more, and why?
(c) The 0012 has c_m,c/4 ≠ 0 at 4°. What does that say about the position of its aerodynamic centre in inviscid flow?

Answer form: a table; two comparisons and a sentence; a sentence.

## 6. (Judgment: where not to trust the panel Cp)

In the lesson's convergence study the largest Cp error over the whole surface stays near 0.5 as N grows, while the error on 5–95 % of the chord falls at second order, and the c_m error stops falling at about 1.5 × 10⁻⁵.

(a) Explain the stalled whole-surface error from the trailing-edge corner flow of a finite-angle edge (`P10.03`, §2.4).
(b) A boundary-layer code (`P10.06`) needs u_e(s) right up to the trailing edge. What should it do about the last few panels, and does the choice matter for the drag?
(c) Rank the following as sources of error in a predicted c_l for the trainer's 2412 at cruise: panel count (N = 160), the corner at the trailing edge, the boundary layer, the finite wing. Justify each rank with a number from `P10.03`–`P10.06`.

Answer form: a short argument; a decision and a sentence; a ranked list.
