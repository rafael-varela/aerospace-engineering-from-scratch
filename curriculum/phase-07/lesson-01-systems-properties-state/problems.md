# P07.01  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact, which come from the NASA fits, and which are estimates.

Several problems use your personal seed **s**: the day of the month you were born (1–31). Constants: R_U = 8.314 462 618 J/(mol·K), p° = 1 bar = 10⁵ Pa. Use your `main.py` functions where a problem says so; otherwise work by hand with a calculator and check with the code.

## 1. (Warm-up: nitrogen at your temperature)

Take N₂ at T = 300 + 20s K and p = 2 bar.

(a) Compute cp (per mole and per kilogram), cp/R and γ from the NASA fit.
(b) Compute h and s°, then s at 2 bar, per mole and per kilogram.
(c) Compute v and ρ from the ideal-gas law.
(d) Estimate cp by hand from [h(T + 10) − h(T − 10)]/20 using `h_molar`, and give its relative difference from (a).

Answer form: eight numbers with units and one relative difference.

## 2. (Derivation: where the NASA forms and the ideal-gas relations come from)

(a) Starting from the cp/R polynomial, derive the h/(RT) and s°/R forms of the lesson by integration, identifying b₁ and b₂ as integration constants.
(b) For argon in the 200–1000 K interval only a₃ = 2.5 is non-zero. Show that requiring h(298.15 K) = 0 gives b₁ = −745.375, the value in the database.
(c) Show that cp − cv = R and γ > 1 for any ideal gas, starting from h = u + p v.
(d) Derive the entropy of an ideal-gas mixture, s = Σ x_i [s°_i − R_U ln(x_i p/p°)] per mole, and compute the mixing part −R_U Σ x_i ln x_i for the Phase 07 dry air per mole and per kilogram.

Answer form: four derivations, one check and two numbers.

## 3. (Unfamiliar scenario: is γ = 1.29 right for Martian air?)

The Martian atmosphere is mostly CO₂. Problem 3 of `P04.07` took γ ≈ 1.29 near the surface, where T ≈ 210 K.

(a) From the NASA CO₂ fit, compute cp/R and γ at 210 K and at 300 K.
(b) How many active quadratic degrees of freedom f = 2/(γ − 1) does each value imply? A linear triatomic molecule has 3 translational and 2 rotational ones; what does the excess mean physically?
(c) At which of the two temperatures is γ ≈ 1.29? What is the consequence for a Martian speed-of-sound estimate at 210 K (a ∝ √γ)?
(d) For comparison compute cp/R and γ of H₂ at 300 K and say why H₂ is close to 7/2 R even though it is the lightest molecule.

Answer form: four numbers, two counts and a sentence, a judgment with a percentage, two numbers and a sentence.

## 4. (Numerical, uses your code: a combustor without constant cp)

Air enters a combustor at 800 K and 2 MPa and receives q = (600 + 20s) kJ/kg at constant pressure.

(a) Find the exit temperature with `exit_temperature` (variable cp) and with constant cp = 1005 J/(kg·K). Give the difference in kelvin.
(b) Find the mean cp over the process.
(c) Use `air_state(p=2e6, h=...)` to complete the exit state: v, u and s.
(d) The engine's turbine is rated for 1650 K. With your q, does the constant-cp designer think it is safe? Is it?

Answer form: three temperatures, one cp, three properties, two sentences.

## 5. (State postulate in practice: an isentropic compressor)

Air at 288.15 K and 101 325 Pa is compressed isentropically (s₂ = s₁) to a pressure ratio π = 10 + s/3.

(a) Compute s₁ with `air_s`, then find T₂ with `air_state(p=π p₁, s=s₁)`.
(b) Compare with the constant-γ formula T₂ = T₁ π^((γ−1)/γ) with γ = 1.4.
(c) Compute the work per kilogram h₂ − h₁ both ways (the second with cp = 1005 J/(kg·K)).
(d) Which two properties fixed state 2? Could you have used (T₂, h₂) instead? Explain.

Answer form: three numbers, one comparison, two numbers, a short argument. (Lesson 05 develops this calculation into a compressor-efficiency calculator.)

## 6. (Judgment: which model, and why?)

For each case say whether this lesson's model (ideal gas, NASA cp(T), fixed dry-air composition) is adequate, which assumption fails first if not, and what you would use instead. Two sentences each.

(a) Cabin air at 0.8 bar and 22 °C in an airliner's air-conditioning duct.
(b) Gas leaving a gas-turbine combustor at 1900 K.
(c) Nitrogen in a 30 MPa pressurant tank at 300 K.
(d) Air behind the bow shock of a capsule re-entering at 7 km/s, where the equilibrium temperature is several thousand kelvin.
(e) Air in a wind tunnel's settling chamber at 10 bar and 300 K, used to compute the test-section speed of sound.

Answer form: five short paragraphs.
