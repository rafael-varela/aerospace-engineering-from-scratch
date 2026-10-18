# P08.07  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results come from exact solutions (laminar flow), which from correlations (Colebrook, Haaland, Blasius, loss coefficients K) and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Properties (typical values): air at 20 °C, 101 325 Pa: ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s; RP-1 ρ = 810 kg/m³, μ = 2.0 × 10⁻³ Pa·s; liquid oxygen near 90 K ρ = 1141 kg/m³, μ = 1.95 × 10⁻⁴ Pa·s; hydrazine near 20 °C ρ = 1004 kg/m³, μ = 0.9 × 10⁻³ Pa·s. Roughness (typical): drawn tubing 1.5 µm, galvanised sheet-metal duct 0.15 mm.

## 1. (Warm-up: a cabin air duct)

A cabin air-distribution duct of diameter 0.10 m and length 8 m (galvanised, ε = 0.15 mm) carries air at 20 °C at V = 5 + s/5 m/s.

(a) Compute Re and say which regime the flow is in.
(b) Find the Darcy friction factor from Colebrook (by your Newton solver) and from Haaland, and their difference in per cent. How many Newton iterations did Colebrook need from the Haaland start?
(c) Compute the friction pressure drop, the volume flow and the fan power needed to overcome it. Add four 90° bends (K = 0.3 each), an inlet (K = 0.5) and an outlet (K = 1.0): what is the total drop?
(d) Compare f with the smooth-pipe value at the same Re and with the fully rough limit. Is the duct closer to smooth or to fully rough?

Answer form: a number and a regime; two friction factors, a percentage and a count; three numbers and a fourth; two friction factors and a sentence.

## 2. (Derivation: the laminar friction factor, Newton for Colebrook, and the D⁻⁵ law)

(a) From the Hagen–Poiseuille profile of `P08.06`, show that the Darcy–Weisbach form dp = f(L/D)ρV²/2 with f = 64/Re is exactly the Hagen–Poiseuille law dp = 128μLQ/(πD⁴).
(b) Write the Colebrook equation as g(x) = x + 2 log₁₀(a + bx) = 0 with x = 1/√f. Derive the Newton iteration, and show that g is increasing and concave for x > 0. Explain why, after the first step, every Newton iterate lies at or to the left of the root and the iterates increase monotonically to it, so that Newton converges from any start between 0 and the root. What can go wrong from a start far to the right of the root? Derive the fixed-point iteration x_{k+1} = −2 log₁₀(a + bx_k) and its linear convergence rate.
(c) At fixed mass flow, show that the friction drop scales as D⁻⁵ when f is constant (fully rough) and as D^{−4.75} when f follows Blasius's 0.3164 Re^{−1/4}. Check both exponents with your `segment_loss`.

Answer form: three derivations; two exponents from your code.

## 3. (Unfamiliar scenario: the oxidiser side)

The same engine as the lesson burns at a mixture ratio O/F = 2.3, so its liquid-oxygen line carries 3.68 kg/s. The line is 1.5 m of 25 mm drawn tubing with a rounded tank outlet (K = 0.05), two bends (0.3 each), a ball valve (0.05), a filter (1.5), a check valve (2.0) and the exit into the manifold (1.0). The injector's oxidiser orifices are sized for a 0.4 MPa drop at Cd = 0.7; the LOX column above the injector is 0.6 m at 3 g; the chamber is at 2.0 MPa.

(a) Find V, Re and f. Why is Re almost twenty times the fuel line's, and what does it do to f?
(b) Build the budget with `feed_line_budget` and give the tank pressure. Which term dominates the line loss?
(c) Compare with the fuel tank pressure of the lesson. If both tanks are held at their pressures and the chamber pressure rises by 5 %, both flows fall. Find the new flows and the new O/F (use a root finder on the budget). Is the engine's mixture ratio robust to a chamber-pressure excursion?

Answer form: three numbers and two sentences; a pressure and a term; two flows, a ratio and a sentence.

## 4. (Numerical, uses your code: sizing your own feed line)

Your engine needs ṁ = 1.2 + s/20 kg/s of RP-1 through the lesson's line (2.5 m, the same fittings, chamber 2.0 MPa).

(a) Find the smallest inside diameter that keeps the line loss below 10 % of the chamber pressure, and below 5 %.
(b) Stainless tubing comes in standard outside diameters; consider 15.875 × 0.889, 19.05 × 0.889, 22.225 × 1.245 and 25.4 × 1.245 mm (OD × wall). For each, give the inside diameter, the velocity, Re and the line loss as a percentage of chamber pressure. Which tube do you choose for the 10 % limit, and what margin does it leave?
(c) Name one reason not to choose a much larger tube than needed.

Answer form: two diameters; a four-row table and a choice with its margin; a sentence.

## 5. (Numerical, extends your code: Newton against fixed point, and a laminar line)

(a) At Re = 10⁵, ε/D = 10⁻⁴, count the iterations to a relative error of 10⁻¹² for the fixed-point iteration of 2(b) and for Newton, both from the Haaland start and from x = 7. Tabulate Newton's errors and show that e_{k+1}/e_k² is roughly constant.
(b) A reaction-control thruster's hydrazine line is 1.5 m of 6 mm bore. For ṁ = 0.005 kg/s compute Re, f and the drop, and check the drop against Hagen–Poiseuille. Repeat at 0.01 and 0.02 kg/s and say which result you would not trust, and why.

Answer form: four counts, a table of errors and a ratio; three drops with their regimes and a sentence.

## 6. (Judgment: what controls the tank pressure?)

For the lesson's fuel line, rank the following uncertainties by their effect on the tank pressure needed, and say which you would reduce first: (i) RP-1 viscosity ±30 % (a temperature range); (ii) wall roughness 10 or 30 times the drawn-tubing value (a poorly finished or corroded line); (iii) the fitting loss coefficients K ±50 %; (iv) the injector discharge coefficient ±5 %. Give each effect in kPa. What does the ranking say about where a feed-system designer should spend test time?

Answer form: four changes in kPa, a ranking and two sentences.
