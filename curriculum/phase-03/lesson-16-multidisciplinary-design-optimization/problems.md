# P03.16  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical (and to what tolerance).

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: one Gauss–Seidel sweep by hand)

Take the lesson's aircraft with S = 16.8 m² and span b = 10 + s/10 m.

(a) Starting from W_wing = 0, compute by hand (calculator, the formulas of 2.1–2.2) the first discipline-1 outputs L, D, W_fuel, W_TO and M_root, then the discipline-2 output W_wing. Check against `aero_mission` and `structures`.
(b) Estimate the coupling strength k from the two partials (6) at your point and predict how many further sweeps bring the relative error in W_wing below 10⁻⁶. Check with `mda_gauss_seidel`.
(c) Would you bother with Newton for this k? Say why in one sentence.

Answer form: six numbers with units, k, a sweep count, one sentence.

## 2. (Derivation: coupling, adjoint, deflection)

(a) From (1)–(3) derive ∂M_root/∂W_wing at fixed (b, S), including the Breguet factor; and from (4) derive ∂W_wing/∂M_root. Write k and show it does not depend on W₀ except through the exponential e^{κφ}.
(b) Show that the Gauss–Seidel error in W_wing obeys e⁺ = k e to first order and that the iteration is a contraction iff |k| < 1.
(c) Starting from the direct formula (7), derive the adjoint formula (8). Count linear solves for n_x design variables and n_F functions for each method, and state the break-even.
(d) A beam with the bending stress equal to σ_a at every section has curvature 2σ_a/(Eh). Integrate twice over the semi-span to get the tip deflection at ultimate load, divide by 1.5 for limit load, and show that δ/(b/2) = σ_a b²/(3Eτ S). Evaluate the aspect-ratio limit for δ_max = 0.06 and for a carbon spar with E = 140 GPa, σ_a = 400 MPa.

Answer form: two derivations ending in formulas, a solve count, two aspect-ratio numbers.

## 3. (Unfamiliar scenario: an electric aircraft that must carry its battery)

A small electric aircraft has W₀ = 8 kN without battery, cruises at L/D = 15 with a propulsive-chain efficiency η = 0.75, and carries a battery of specific energy e = 250 Wh/kg. The cruise energy for range R is E = W R/(η L/D) and the battery weight is W_b = g E/e. The design must fly R = 200 + 10s km.

(a) Write the coupled fixed point W = W₀ + W_b(W) and solve it exactly. What is the coupling strength k, and what does k = 1 mean physically?
(b) Starting from W = W₀, how many Gauss–Seidel sweeps bring the relative error below 10⁻⁶? How many Newton steps?
(c) Compute dW/dR two ways: the total derivative through the loop and the partial with W frozen at its converged value. By what factor do they differ, and what is that factor in terms of k?
(d) At what range does no design close? What would the frozen-weight analysis say at that range?

Answer form: W and W_b in newtons, k, two iteration counts, two derivatives in N/km, one range.

## 4. (Numerical, uses your code: your own aircraft)

Set W₀ = (20 + s/2) kN and R = (1000 + 50s) km, everything else as in the lesson.

(a) Converge the MDA at the baseline (12 m, 16.8 m²) by Gauss–Seidel and Newton; report y, k and the observed Gauss–Seidel ratio; compare the Newton residual history with the quadratic law.
(b) Compute dW_fuel/d(b, S) by the adjoint and by finite differences (h = 10⁻⁶ relative) and report the relative difference; report also the frozen partial ∂W_fuel/∂b and the percentage by which it misjudges the span benefit.
(c) Optimise in MDF (adjoint) and IDF; report b, S, AR, fuel, the active constraints and the saving over the baseline.
(d) Estimate the multiplier of the deflection constraint by least squares on the KKT condition and state the fuel saved per 10 % of extra deflection allowance.

Answer form: a table (quantity, Gauss–Seidel, Newton), a derivative table (adjoint, FD, partial), a design table (MDF, IDF), one multiplier with its interpretation.

## 5. (Numerical, uses your code: the finite-difference floor)

For the lesson's baseline, run `fd_convergence_study` with h = 10⁻¹, 10⁻², …, 10⁻⁷ in four settings: Newton at 10⁻¹³; Gauss–Seidel at 10⁻¹²; Gauss–Seidel at 10⁻⁶; Gauss–Seidel at 10⁻³.

(a) Tabulate the errors and fit the slope over the steps where the error still falls.
(b) Explain why the floor of the loose Gauss–Seidel runs is roughly constant in h rather than growing like 1/h (think about which errors are common to both sides of the central difference).
(c) An SLSQP run with finite-difference gradients uses h ≈ 1.5 × 10⁻⁸ relative by default. From your table, what gradient accuracy would it get with each MDA setting, and at which setting would you expect it to stall or wander?

Answer form: a 7 × 4 table, a fitted slope, two short paragraphs.

## 6. (Judgment: scaling up)

A team proposes an MDO of a regional aircraft with four disciplines (aerodynamics by a panel code, structures by a finite-element model, propulsion by an engine deck, mission), 150 design variables, 12 constraints and one objective. The aerodynamics and structures codes are commercial and give no derivatives; the engine deck is a table.

(a) Which MDA solver do you recommend and what would you measure first to justify it? What would you do if k turned out to be 0.9?
(b) Compare the cost of finite-difference, direct and adjoint total derivatives per optimiser iteration for this problem. Which can the team actually use, given the codes they have, and what would it take to use the cheapest?
(c) MDF or IDF? Give one argument for each and the one that decides it here.
(d) The first result shows the span at the gate limit and a 7 % fuel saving. List, in order, the three model assumptions you would check before believing the 7 %, and say what evidence each check needs.

Answer form: four short paragraphs with at least one number each.
