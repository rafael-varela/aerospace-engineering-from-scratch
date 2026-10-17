# P08.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact (similarity solution), which are numerical, and which are assumptions about the flow (laminar, zero pressure gradient).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Air at sea level (ISA): ρ = 1.225 kg/m³, μ = 1.7894 × 10⁻⁵ Pa·s, ν = 1.4607 × 10⁻⁵ m²/s. Use the Blasius constants your code computes (or the published f″(0) = 0.332057 if you have not finished it).

## 1. (Warm-up: a wind-tunnel fin)

A thin flat fin of chord 0.20 m and span 0.6 m sits at zero incidence in a wind tunnel running at U = 20 + s/2 m/s with sea-level air. Treat each side as a flat plate with a laminar boundary layer from the leading edge.

(a) Compute the Reynolds number at the trailing edge, and δ99, δ* and θ there.
(b) Compute c_f and the wall shear stress at mid-chord and at the trailing edge. By what factor do they differ, and why that factor?
(c) Compute C_F and the friction drag of the fin (both sides).
(d) Is the laminar assumption plausible at this Reynolds number? Say what lesson `P08.12` will add to this judgment.

Answer form: four numbers with units; four numbers and a ratio with a reason; a coefficient and a force; two sentences.

## 2. (Derivation: similarity and a trick that removes the shooting)

(a) Starting from Prandtl's boundary-layer equations for a flat plate at zero pressure gradient, u u_x + v u_y = ν u_yy and u_x + v_y = 0, introduce the stream function ψ = √(νUx) f(η) with η = y√(U/(νx)). Show that u = U f′(η), that v = ½√(νU/x)(η f′ − f), and that the momentum equation becomes f‴ + ½ f f″ = 0 with f(0) = f′(0) = 0 and f′(∞) = 1. State where each assumption enters.
(b) Show that if F(ζ) solves F‴ + ½ F F″ = 0 with F(0) = F′(0) = 0 and F″(0) = 1, then f(η) = a F(aη) also solves the Blasius equation, with f″(0) = a³ and f′(∞) = a² F′(∞). Deduce f″(0) = [F′(∞)]^(−3/2): **one** forward integration and no root finding. (This is Töpfer's transformation.)
(c) Integrate the momentum equation across the layer to obtain the von Kármán momentum integral for a flat plate, dθ/dx = c_f/2. Then show from the similarity form that θ = 2 f″(0) x/√Re_x, so θ√Re_x/x = 2 f″(0) exactly.

Answer form: three derivations, each final result boxed.

## 3. (Unfamiliar scenario: why wind-tunnel test sections diverge)

A small low-speed tunnel has a square test section 0.5 m × 0.5 m and 2 m long, entered at U = 15 + s/5 m/s with sea-level air. Assume a laminar boundary layer starts on all four walls at the test-section entrance.

(a) Find δ99 and δ* on each wall at the end of the test section.
(b) The displacement thickness acts like a thickening of the walls. By mass conservation, by what fraction does the core velocity rise between the entrance and the exit if the walls are parallel? Compare with the approximation 4δ*/W.
(c) By what angle must each wall diverge at the exit so that the core velocity stays constant? Compare it with the outward velocity v at the edge of the layer, v/U = (η f′ − f)/(2√Re_x).
(d) A real tunnel's walls are turbulent (lesson `P08.14`). Would the divergence need to be larger or smaller? Answer without a calculation.

Answer form: two lengths; a fraction and a comparison; an angle and a velocity ratio; one sentence.

## 4. (Numerical, uses your code: shooting accuracy)

Let η_max = 6 + s/5.

(a) With your `shoot_blasius`, compute f″(0) with n = 10 η_max, 20 η_max, 40 η_max and 80 η_max RK4 steps. As the reference, use `blasius_bvp(η_max, tol=1e-10)` at the same η_max. Report the errors and the observed orders. Why must the reference use the same η_max?
(b) With 100 steps per unit η, compute f″(0) for η_max = 4, 5, …, 10, 12 and the difference from 0.332057336. How large must η_max be for six correct digits? Why does too small an η_max make f″(0) too large?
(c) Use Problem 2(b): integrate once with F″(0) = 1 to ζ = 10 with 2000 steps and compute f″(0) = [F′(10)]^(−3/2). How many digits does it get right, and what did it cost compared with your secant shooting?

Answer form: an error table, orders and a sentence; a table, a number and a sentence; a number, a digit count and a cost comparison.

## 5. (Numerical, uses your code: the PDE forgets its starting profile)

Run `march_boundary_layer` with U = 5 + s/3 m/s, sea-level ν, L = 1 m, y_max = 3 × 5L/√Re_L, ny = 300, and nx = 200, 400, 800, 1600.

(a) Report c_f√Re_x at x = L for each nx (use `wall_cf_from_profile`), its error against 2f″(0), and the observed order. Which discretisation error dominates, and how do you know?
(b) Richardson-extrapolate the two finest values with the order you found. How close is the result to the Blasius value?
(c) At nx = 1600, plot u/U against η at x = 0.25 m and x = 1 m with the Blasius f′(η). Report the largest difference. The march started from a uniform profile, not from Blasius: what does the collapse tell you about the boundary layer's memory?

Answer form: a table and an order with a reason; a number; a plot, a number and a sentence.

## 6. (Judgment: where Blasius stops describing the wing)

The lesson's engineering answer treats the light aeroplane's wing (3000 m, 60 m/s, chord 1.5 m) as a flat plate with a laminar boundary layer from the leading edge to the trailing edge.

(a) Near the leading edge the boundary-layer approximation itself fails because δ is not small compared with x. Find the Reynolds number and the distance from the leading edge at which δ99/x = 0.1.
(b) List three further ways the real wing differs from the Blasius plate. For each, say whether it raises or lowers the friction drag compared with 29.7 N, and which later lesson deals with it.
(c) Give one measurement on a real wing (in flight or in a tunnel) that would show the laminar model has failed, and what reading you would expect if it had.

Answer form: a Reynolds number and a length; three items with a direction and a lesson each; one measurement and an expected reading.
