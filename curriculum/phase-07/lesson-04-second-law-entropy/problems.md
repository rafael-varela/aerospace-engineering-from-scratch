# P07.04  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which inputs are typical, illustrative values. A negative entropy generation is never an answer: if you get one, find the mistake.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Dry air: R = 287.045 J/(kg·K), NASA variable cp (the lesson's helpers). Dead state T₀ = 288.15 K.

## 1. (Warm-up: three generators by hand)

(a) An avionics box at 340 K rejects Q̇ = 100 + 10s W to cabin air at 290 K. Compute Ṡ_gen.
(b) Two 1 kg aluminium blocks (c = 900 J/(kg·K)), one at 300 + 2s K and one at 300 K, touch inside an insulated box. Find T_f and S_gen, and compare with the small-ΔT formula C ΔT²/(4T_m²) with T_m the mean temperature.
(c) Air is throttled from (1 + 0.1s) MPa to 0.1 MPa. Find s_gen and the lost work T₀ s_gen per kilogram.

Answer form: one number; three numbers; two numbers.

## 2. (Derivation: the second law as a balance)

(a) From the Clausius inequality, show that S₂ − S₁ = ∫δQ/T + S_gen with S_gen ≥ 0, and write the steady-flow entropy balance for a control volume.
(b) Derive T ds = du + p dv and T ds = dh − v dp, and from them ds of an ideal gas with variable cp.
(c) Prove that two bodies equilibrating always give S_gen ≥ 0. (Hint: ln y ≤ y − 1, or the concavity of ln.)
(d) Show that for a two-reservoir cycle S_gen ≥ 0 is equivalent to η ≤ 1 − T_C/T_H.

Answer form: four derivations.

## 3. (Unfamiliar scenario: where does a satellite's entropy go?)

A satellite's electronics turn 200 W of electrical work into heat at 350 K; the heat is conducted to a radiator at 290 K and radiated to deep space, taken as a reservoir at 3 K (typical, illustrative).

(a) Compute the entropy generated in each of the three steps (dissipation of work at 350 K, conduction 350 → 290 K, radiation 290 → 3 K) and the total.
(b) Show that the total is W/T_space, independent of the intermediate temperatures. Why?
(c) Which step dominates? Does that mean the radiator is "badly designed"? What does the choice of dead-state temperature have to do with your answer?

Answer form: four numbers, a one-line proof, a paragraph.

## 4. (Numerical, uses your code: entropy is a property)

Take air from (300 K, 100 kPa) to T₂ = 500 + 10s K, p₂ = (2 + s/5) bar.

(a) Compute Δs with `ds_ideal_gas`.
(b) Integrate along an isobaric-then-isothermal path, an isothermal-then-isobaric path, and the path T = T₁ + (T₂ − T₁)λ + 100 sin(2πλ), p = p₁(p₂/p₁)^λ (1 + 0.3 sin(πλ)), with `ds_path_integral`.
(c) For the third path, compute the error at n = 101, 201, 401, 801 and the observed order.
(d) Compare Δs with two constant-cp estimates: cp = 1005 J/(kg·K) and the mean cp = Δh/ΔT. Which is better, and why?

Answer form: one number, three numbers, four errors and an order, two numbers and a sentence.

## 5. (Numerical, uses your code: redesigning the bleed chain)

Copy `BLEED` and change the valve outlet to p_prv = (0.2 + 0.01s) MPa and the precooler outlet to 450 + s K, keeping the duct outlet at 0.28 MPa.

(a) Build the ledger with `bleed_air_ledger` and print it ranked, with shares and lost work.
(b) Compare the total with the lesson's 535.8 W/K. Explain why moving the valve's outlet pressure barely changes the total.
(c) With your precooler outlet, compute the precooler's S_gen for fan-air flows of 0.8, 1.2 and 2.0 kg/s. Which would you choose, and what limits how far you can go?
(d) Propose one change that reduces the total substantially, and estimate the reduction with your code.

Answer form: a table, two numbers and a sentence, three numbers and a sentence, a proposal with a number.

## 6. (Judgment: claims, choices and limits)

(a) A supplier claims a heat engine that takes 300 J at 900 K, rejects 90 J at 300 K and delivers 210 J per cycle. Check it with `check_cycle` and explain the verdict in one sentence.
(b) Another engine takes 500 J at 1200 K and 200 J at 800 K and rejects 450 J at 350 K. Is it possible? What is its efficiency against the Carnot efficiency between 1200 K and 350 K?
(c) A colleague says, "the mixer's S_gen is tiny, so mixing streams is always harmless." When is that wrong?
(d) Your ledger uses T₀ = 288.15 K. At cruise the outside air is at about 220 K. Which numbers change and which conclusion does not?
(e) Name one process in this phase where the ideal-gas throttling formula s_gen = R ln(p₁/p₂) would be wrong, and why.

Answer form: five short paragraphs.
