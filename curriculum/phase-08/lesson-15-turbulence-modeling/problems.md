# P08.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer, and say which results are exact, which come from the model, and which from a correlation.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Air at 20 °C: ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s.

## 1. (Warm-up: wall units of a duct)

A flat cabin-air duct has a gap of 40 + s mm (treat it as a channel of height 2h) and carries air at a bulk speed of 4 m/s.

(a) Compute Re_b. With `duct_flow`, find Re_τ, u_τ, τ_w and the pressure drop over 2 m of straight duct, for both outer mixing lengths.
(b) Compare τ_w with Dean's correlation.
(c) Give the first-cell height for y⁺ = 1 and for y⁺ = 30, in micrometres and millimetres.

Answer form: six numbers per closure; two percentages; two lengths.

## 2. (Derivation: the channel model and the k–ε constraint)

(a) Starting from the Reynolds-averaged x-momentum equation for fully developed channel flow, derive τ⁺ = 1 − y⁺/Re_τ, and show the same result from a momentum balance on a control volume one metre long spanning the channel.
(b) Derive the velocity-gradient formula g = 2τ⁺/(1 + √(1 + 4l⁺²τ⁺)) from the mixing-length closure. Show its two limits and say why this form is preferred to the usual root formula when l⁺ is small.
(c) For the standard k–ε model in the log layer, derive k, ε and ν_t, and then the constraint κ² = σ_ε √C_μ (C_ε2 − C_ε1).

Answer form: three derivations.

## 3. (Derivation and numbers: the log-law friction law)

(a) Integrate the log law across the half-channel to obtain U_b⁺ = (1/κ) ln Re_τ + B − 1/κ. Which two parts of the real profile does this ignore, and in which direction does each push U_b⁺?
(b) Solve the resulting implicit relation for c_f at Re_b = 2 × 10⁴ and 10⁵ (κ = 0.41, B = 5.0), and compare with Dean's correlation.
(c) Repeat (b) with κ = 0.4327 (the k–ε value) and the same B. How much does the friction change? What does that say about a k–ε code with wall functions that assume κ = 0.41?

Answer form: a derivation and two sentences; two c_f values and two percentages; two c_f values and a sentence.

## 4. (Numerical, uses your code: which constant matters?)

At Re_τ = 1000:

(a) Tabulate c_f, its difference from Dean's correlation, and the fitted B for κ = 0.38, 0.41, 0.42; A⁺ = 22, 26, 30; and both outer mixing lengths.
(b) Rank the three choices by their effect on c_f. Which part of the profile does each one change?
(c) Is the model's uncertainty in c_f larger or smaller than the scatter of Dean's correlation (a few per cent)? What does that imply for a duct design?

Answer form: a table; a ranking with reasons; a comparison and a sentence.

## 5. (Unfamiliar scenario: a propellant feed line is a pipe)

The feed line of a small liquid rocket is a round pipe, not a channel. The mixing-length model carries over with two changes: the total stress is still linear, τ⁺ = 1 − y⁺/R⁺ (R⁺ = R u_τ/ν is the pipe's friction Reynolds number), but the bulk velocity weights each annulus by its area, U_b⁺ = (2/R⁺²)∫₀^{R⁺} u⁺ (R⁺ − y⁺) dy⁺, and Re_D = 2U_b⁺R⁺.

(a) Derive the area weighting and the linear stress in a pipe from a force balance on a cylinder of radius r.
(b) Compute the Fanning friction factor c_f = 2/U_b⁺² for R⁺ = 500 and 2000 with both outer mixing lengths, and compare with the Blasius correlation c_f = 0.0791 Re_D^(−1/4).
(c) Which outer closure works better in the pipe, and why might that be expected?
(d) Give the ratio of centreline to bulk velocity. Why does it matter for a flowmeter placed on the axis?

Answer form: a derivation; a table of four c_f values against the correlation; a sentence; two ratios and a sentence.

## 6. (Judgment: meshing the duct for CFD)

You must set up a RANS simulation of the duct of Problem 1 with a bend.

(a) For a low-Reynolds-number k–ε or k–ω model, how thin must the first cell be, and how many cells would you put between y⁺ = 1 and the centreline with a growth ratio of 1.2?
(b) For wall functions, where must the first cell centre lie, and why is that hard to meet at Re_τ ≈ 180?
(c) Which flow features in the bend will an eddy-viscosity model get wrong, whichever wall treatment you choose?
(d) Your manager asks for "the pressure drop to 2 %". Using the sensitivity of Problem 4, answer in two sentences.

Answer form: a length and a count; a range and a reason; a list; two sentences.
