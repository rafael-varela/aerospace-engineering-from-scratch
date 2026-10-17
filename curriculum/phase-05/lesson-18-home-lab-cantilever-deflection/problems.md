# P05.18  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Material values are typical handbook values, as in the lesson: aluminium E = 69 ± 2 GPa, F_ty = 145 MPa (6063-T5), ρ = 2700 kg/m³; g = 9.81 m/s².

## 1. (Warm-up: predict before you load)

An aluminium strip 20.0 mm wide and 2.00 mm thick is clamped with L = 300 mm free. A bag of mass m = (50 + 5s) g hangs at a = 250 mm.

(a) Compute I and EI.
(b) Compute the deflection at the load point, at x = 150 mm and at the free end x = 300 mm.
(c) Compute the bending stress at the clamp and the safe load (in grams) for an allowable of F_ty/2. Is your mass safe?
(d) With u_b = 0.05 mm, u_h = 0.02 mm, u_a = 1 mm and u_E = 2 GPa, give the predicted deflection at the load point with its standard uncertainty, by the relative-uncertainty rule and with `predict_deflection`.
(e) How much does the strip sag at x = a under its own weight before you hang anything, and why does it not matter?

Answer form: I, EI, three deflections, a stress, a safe mass, δ ± u, a sag.

## 2. (Derivations)

(a) Integrate E I w″ = P(a − x) with w(0) = w′(0) = 0 to obtain w(x) for x ≤ a. Then show that w is a straight line for x ≥ a, find its slope, and write w(x) there.
(b) Show that for a reading at the load point the relative uncertainty of δ = 4Pa³/(Ebh³) is √((u_E/E)² + (u_b/b)² + (3u_h/h)² + (3u_a/a)²). Now suppose the reading is taken at a fixed mark x on the strip and only the load's position a is uncertain. Show that the sensitivity d ln w/d ln a is 3a/(3a − x) for x ≤ a, and evaluate it at x = a. Why is it smaller than 3?
(c) A clamp with rotational stiffness k_θ rotates by θ₀ = Pa/k_θ. Show that the measured compliance at the load point is s = a³/(3EI)·(1 + 3EI/(k_θ a)), so that E computed from s by the rigid-clamp formula is E/(1 + 3EI/(k_θ a)). Evaluate the factor for the lesson's strip (EI = 3.881 N·m², a = 290 mm) with k_θ = 2000 N·m/rad.
(d) Prove that taking the difference of a loaded and an unloaded reading removes the self-weight deflection exactly, whatever its shape, and say which assumption of the beam model this relies on.

## 3. (Unfamiliar scenario: a thin steel rule)

A stainless-steel rule, 25 mm wide and 0.80 mm thick, is clamped so that the load hangs at a = 240 mm. Typical values for annealed 304 stainless: E = 193 GPa, F_ty = 215 MPa (a hardened rule is stronger; use the annealed value to be safe). A bag of mass m = (20 + 2s) g hangs at a.

(a) Compute the linear deflection, the elastica deflection with `elastica_deflection`, the percentage by which linear theory over-predicts, δ/a, and the clamp stress.
(b) Compute the safe load for F_ty/2, and repeat (a) at that load. At what mass does linear theory over-predict by 1 %? Compare δ/a there with the lesson's rule.
(c) With u_h = 0.01 mm, u_E = 5 % (the alloy of a rule is rarely marked) and the other uncertainties as in Problem 1, what is the relative uncertainty of the prediction? Could the 1 % nonlinearity be detected by comparing with the prediction? Could it be detected from the *shape* of the load–deflection curve, and how big is the departure from the straight line at the safe load, in millimetres?
(d) Which extra safety control does this specimen need?

Answer form: two deflections, a percentage, δ/a and a stress at your mass; the same at the safe load; a mass; an uncertainty and two yes/no answers with reasons; a control.

## 4. (Code: reduce a SYNTHETIC dataset)

Generate `make_synthetic_cantilever_data(seed=s)`. Every result in this problem is SYNTHETIC and must be labelled so.

(a) Reduce the load series as the lesson's `main` does (difference the readings, u_d = √2 × 0.30 mm): fit the line and report the slope, intercept, χ² and dof; compute E with `modulus_from_slope` and E_n against 69 ± 2 GPa.
(b) Fit the profile at 700 g with `fit_profile`. Report EI, θ₀ and their uncertainties and χ².
(c) The generator's truth is in `data["truth"]`. Compute the true compliance at x = a including the clamp, a³/(3EI) + a²/k_θ, and the apparent E it implies through the nominal b = 25.0 mm and h = 3.00 mm. Is your E from (a) consistent with it?
(d) Repeat the profile fit on noise-free predicted data with u_w = √2 × 0.01 mm (a dial gauge). What is u(θ₀) now, and would the generator's clamp rotation at 700 g be detected?

## 5. (Design: a 1 % measurement of E)

You want E of an aluminium strip to ± 1 % (k = 1), from the load–deflection slope, using the lesson's strip and reading at the load point. Start from the lesson's SYNTHETIC budget: slope 2.67 %, thickness 3u_h/h with u_h = 0.02 mm, load position 3u_a/a with u_a = 1 mm, width u_b/b with u_b = 0.05 mm.

(a) Tabulate the four terms and the total.
(b) A dial gauge cuts the slope term to 0.09 %. Which terms must you improve next, and to what? Give a combination of u_h and u_a that reaches 1 %, and the resulting total.
(c) Even with that budget, a C-clamp with k_θ = 2000 N·m/rad biases E by about 2 % (Problem 2c). Is that inside or outside your 1 % target? Describe how to measure and remove it, and what it does to your uncertainty.
(d) Would reading at a fixed mark instead of at the load point help with the load-position term? Use Problem 2(b).

## 6. (Judgment: a surprising modulus)

A classmate clamped a strip sold as "aluminium flat bar 25 × 3 mm", 310 mm long, with about 20 mm under the clamp, and hung loads at "a = 290 mm, measured from the end of the strip under the clamp". They read the deflection with a ruler, fitted the slope, and report E = 56 ± 2 GPa (k = 1), with E_n = 2.3 against 69 ± 2 GPa. Their conclusion: "the strip is not aluminium, or the handbook value is wrong for this alloy."

Write a review of at most one page. Rank at least three candidate causes. For each, compute how large it would have to be to explain the result (an error in a, an error in h, a clamp stiffness k_θ), and say how plausible that is. Propose the cheapest test that discriminates between them, and say what result of each test would confirm each cause. Say whether the conclusion as written is acceptable.
