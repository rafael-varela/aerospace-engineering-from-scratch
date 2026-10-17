# P07.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are predictions, which are fits to data, and whether the data are measured or **SYNTHETIC**.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Properties: 6061-T6 and AISI 304 from Lienhard Table A.1; air from Table A.6; σ = 5.670374419 × 10⁻⁸ W/(m² K⁴). No problem here asks for water hotter than 60 °C or for any heater other than tap-hot water.

## 1. (Warm-up: predict your own bar)

Your bar is the lesson's 25.4 mm 6061 round bar cut to L = 40 + s mm, hung horizontally in 21.0 °C air, ε = 0.09.

(a) Compute V, A, L_c and ρcV.
(b) With `predict_tau` at θ = 14 K, find h_conv, h_rad, τ and Bi.
(c) Will τ be longer or shorter than the lesson's 50 mm bar? Explain from V/A, not from the code.

Answer form: four numbers; four numbers; a sentence.

## 2. (Derivation: the lumped model, the probe and the small-Bi limit)

(a) Derive θ = θ₀e^(−t/τ) from the first law for an isothermal body. Then, by comparing the conduction resistance L_c/(kA) with the surface resistance 1/(hA), show why the internal temperature differences are a fraction of order Bi of θ.
(b) Solve τ_p dT_m/dt = T − T_m with T = T∞ + θ₀e^(−t/τ) and T_m(0) = T∞. Show that for t ≫ τ_p the reading is an exponential with the same τ and amplitude θ₀τ/(τ − τ_p).
(c) Show that 4εσT_m³ with T_m = (T + T∞)/2 underestimates the exact h_rad by a relative amount of order ((T − T∞)/(T + T∞))², and evaluate the exact and approximate h_rad for ε = 0.09 at θ = 5, 14 and 28 K above 21 °C.
(d) From λJ₁(λ) = Bi J₀(λ), show that λ₁² → 2Bi as Bi → 0, and hence that the cylinder's series solution reduces to exp(−2Bi Fo), the lumped model for an infinite cylinder.

Answer form: four derivations and a short table.

## 3. (Unfamiliar scenario: the bar in a thermal-vacuum chamber)

Spacecraft units are tested in a vacuum chamber, where there is no convection. Put the lesson's bar (ε = 0.09) in such a chamber with walls at 21 °C.

(a) Compute h_rad at θ = 14 K and the linearised time constant. How does it compare with the 32.8 min in room air?
(b) Find the exact time to cool from 50 °C to 30 °C by integrating ρcV dT/(εσA(T⁴ − T_w⁴)), and compare with the linearised estimate τ ln(29/9).
(c) Repeat (a) and (b) for the bar wrapped in black tape (ε = 0.9, typical, illustrative).
(d) If the chamber walls were an ideal cold shroud at 0 K, derive t = ρcV/(3εσA)(1/T₂³ − 1/T₁³) for cooling from T₁ to T₂ and evaluate it for both emittances.

Answer form: two numbers and a ratio; two times; four numbers; a formula and two times.

## 4. (Numerical, uses your code: your own SYNTHETIC run)

Generate a SYNTHETIC dataset with `make_synthetic_cooling_data(seed = 715 + s)` (label everything SYNTHETIC).

(a) Fit it with `fit_loglinear` (t ≥ 2 min, θ ≥ 1 K, T∞ = the mean room reading), `fit_exponential` with T∞ free, and `fit_lagged`. Tabulate τ, its uncertainty and χ²/dof for each.
(b) Compute E_n of the log-linear τ against the prediction (re-run `predict_tau_uncertainty` at your measured room temperature).
(c) Fit the three windows 2–20, 20–40 and 40–60 min. Is h constant?
(d) Fit `fit_variable_h` with ε = 0.09 and report f_conv. The generator's truth is in the returned dict: how close did you get, and why can you not recover f_conv and ε separately?

Answer form: a table; a number; three τ and a sentence; a number and two sentences.

## 5. (Numerical, uses your code: when does the lumped model fail?)

Take the lesson's bar shape (25.4 × 50 mm) and compare four cases: aluminium in air (h = 6.9 W/(m² K)), aluminium in a stirred water bath (h = 500 W/(m² K), typical, illustrative), AISI 304 stainless steel in air, and 304 in the water bath.

(a) For each, compute Bi and τ.
(b) Use `short_cylinder` at t = τ to compute the centre and corner temperatures, and compare with the lumped value e⁻¹ = 0.3679.
(c) In the experiment the bar is warmed in water and then cooled in air. Which of the four cases describe the warm-up, and what does your answer imply for how long the bar must stay in the water?

Answer form: a table; a sentence on each case; a recommendation with a time.

## 6. (Judgment: why is the measured τ 20 % short?)

A learner's measured τ is 20 % shorter than the prediction, with χ²/dof = 60 on a single exponential.

(a) Estimate the τ the prediction gives if convection is 30 % stronger than the correlation (a draughty room), and if the bar's emittance is 0.3 rather than 0.09. Which is more plausible, and how would you tell them apart?
(b) Re-fit the SYNTHETIC seed-730 data with T∞ set 1 K too high. What τ results, and in which direction?
(c) A kitchen probe has a steel tip of a few grams. Estimate the fraction it adds to ρcV for a 2 g and a 5 g tip (take c = 500 J/(kg K), typical, illustrative). In which direction does the probe's stem, acting as a fin, move τ?
(d) Rank the four explanations and propose one test for each.

Answer form: two times and a sentence; a time and a direction; two percentages and a direction; a ranked list.
