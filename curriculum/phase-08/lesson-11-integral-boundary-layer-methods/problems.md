# P08.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results come from Thwaites's method (an integral closure), which from an exact similarity solution, and which from the marching PDE solver.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Air: at 3000 m ISA, ν = 1.8630 × 10⁻⁵ m²/s and ρ = 0.9091 kg/m³; at sea level, ν = 1.4607 × 10⁻⁵ m²/s.

## 1. (Warm-up: the wing's leading edge)

Near the stagnation point of a rounded leading edge of radius r, the potential-flow edge velocity is close to that of a cylinder, Ue ≈ 2U sin(x/r) ≈ (2U/r)x = a x, with x the arc length from the stagnation point. Take the light aeroplane at 3000 m and U = 60 m/s with r = (15 + s/2) mm.

(a) Compute a, and Thwaites's θ, λ, H, δ* and l at the stagnation point.
(b) Compare θ and δ* with the exact Hiemenz values from `falkner_skan(1.0)` (θ√(a/ν) and δ*√(a/ν)). Which does Thwaites get better, and why?
(c) Compute the wall shear stress at x = 2 mm by Thwaites and from the exact f″(0). How good is the linear approximation of sin(x/r) there?

Answer form: six numbers with units; two comparisons with percentages and a sentence; two stresses and a percentage.

## 2. (Derivation: from the momentum integral to a quadrature)

(a) Integrate the boundary-layer momentum equation u u_x + v u_y = Ue Ue′ + ν u_yy across the layer and derive von Kármán's momentum integral equation dθ/dx + (2 + H)(θ/Ue)Ue′ = c_f/2. State where you use continuity and the wall conditions.
(b) Multiply by Ue θ/ν and show that Ue d(θ²)/dx = 2[l − (2 + H)λ] with λ = θ²Ue′/ν and l = τ_w θ/(μUe). Thwaites observed that the right-hand side of exact solutions is close to 0.45 − 6λ. Use that line to derive θ²Ue⁶ = θ₀²Ue₀⁶ + 0.45ν∫Ue⁵dx.
(c) For a wedge flow Ue = Cx^m show that λ = 0.45m/(5m + 1) everywhere, and for Howarth's flow Ue = U₀(1 − x/L) show that λ = −0.075[(1 − x/L)^(−6) − 1]. Find Thwaites's separation point x_s/L in closed form.
(d) Show that the separation angle Thwaites predicts on a circular cylinder does not depend on U, R or ν.

Answer form: four derivations; one number for (c).

## 3. (Unfamiliar scenario: a cabin air-duct diffuser)

A cabin air duct at sea level carries air at U₀ = 10 m/s along a straight section of length l₀ = 0.2 + s/100 m and then into a gentle diffuser in which the core velocity falls as Ue = U₀/(1 + ξ/L_d) with L_d = 3 m (ξ is the distance from the diffuser entrance). Treat the wall boundary layer as two-dimensional and laminar.

(a) Ignore the straight section (the layer starts at the diffuser entrance). Show that Thwaites gives λ = −0.1125[(1 + ξ/L_d)⁴ − 1], and find where the layer separates and by how much Ue has fallen there.
(b) Now include the straight section: start the diffuser with Thwaites's flat-plate θ₀ at l₀. Find λ at the entrance and the separation point with your code, and the fall in Ue. Explain the difference from (a) in one sentence.
(c) Check (b) with `march_boundary_layer`: build Ue(x) for the straight section followed by the diffuser and march from the start of the straight section. Compare the separation points.
(d) What would you change in the duct to avoid laminar separation, and what does a separated layer do to the diffuser's pressure recovery?

Answer form: an expression, a length and a percentage; a λ, a length, a percentage and a sentence; two lengths and a percentage; two sentences.

## 4. (Numerical, uses your code: how good is the closure?)

Let m₁ = s/100 − 0.08, and take also m = −0.05 and m = 0.5.

(a) For each m, compute the exact Falkner–Skan solution with `falkner_skan(2m/(m + 1))` and Thwaites's closed form with `wedge_thwaites(m)`. Tabulate θ√Re_x/x, c_f√Re_x and H from both (Re_x = Ue x/ν). Hint: in the Falkner–Skan scaling, θ√Re_x/x = θ_η√(2/(m + 1)) and c_f√Re_x = 2f″(0)√((m + 1)/2).
(b) Plot the percentage error of Thwaites in θ and c_f against m. Where is the method best, and where is it worst?
(c) The exact solutions lose their wall shear at m = −0.0904 (β = −0.1988), where λ_exact = −0.068. Thwaites's criterion is λ = −0.09. Explain why a method whose separation λ is "wrong" by 30 % still predicts the separation *point* of Howarth's flow within 3 %.

Answer form: a table of 18 numbers; a plot and two sentences; a paragraph.

## 5. (Numerical, uses your code: two routes to separation)

(a) Run `march_boundary_layer` on Howarth's flow (U = 10 m/s, sea-level ν, L = 1 m, y_max = 20√(0.12νL/U), ny = 300) on linear x-grids from 10⁻⁴ m to 0.2 m with nx = 1000, 2000, 4000, 8000 steps. Report x_s/L from `march_separation`, and estimate the converged value. Compare with the published 0.1199 and with Thwaites.
(b) Repeat for the circular cylinder of radius 0.05 m in a 20 m/s stream with ν = 1.5 × 10⁻⁵ m²/s, marching from 10⁻³R to 0.7πR, with (nx, ny) = (2000, 150), (4000, 300), (8000, 300), (16000, 600). Report the separation angle and compare with Thwaites.
(c) Why does `march_boundary_layer` stop at separation instead of continuing into the reversed flow? What would a solver need to go further?

Answer form: four numbers and an estimate with two comparisons; four angles and a comparison; two sentences.

## 6. (Judgment: a laminar-flow wing)

The lesson's experiment 5 puts the suction peak of a model wing section at 0.1c, 0.3c, 0.5c and 0.7c with the same peak velocity (1.25U) and trailing-edge velocity (0.95U), and finds laminar separation at x/c = 0.506, 0.508, 0.540 and 0.700.

(a) Explain why moving the peak back from 0.5c to 0.7c makes the layer separate *at* the peak instead of 60 mm behind it.
(b) Real laminar-flow aerofoils do not shed their boundary layer at x/c ≈ 0.5. What usually happens after laminar separation at flight Reynolds numbers, and what does it cost?
(c) Name the assumption of this lesson that fails first near separation, and the measurement on a wind-tunnel model that would show it.

Answer form: three short paragraphs.
