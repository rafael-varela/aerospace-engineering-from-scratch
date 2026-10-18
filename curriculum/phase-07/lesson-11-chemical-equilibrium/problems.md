# P07.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which numbers are equilibrium values and which complete-combustion values. **Everything here is computed: there is no combustion, flame, burning, fuel handling or fuel of any kind.**

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Data: the module's NASA polynomials and species list. Chamber pressures, mixture ratios and combustor conditions are typical, illustrative values.

## 1. (Warm-up: water falling apart, by hand)

Pure water vapour at T = 2500 + 20 s K and 1 bar.

(a) Compute ΔG° and K_p for H₂O ⇌ H₂ + ½O₂ with `g_species`.
(b) Assuming this is the only reaction, solve by hand (or with `brentq`) for the fraction α dissociated and the mole fractions.
(c) Run `equilibrium_tp({"H": 2, "O": 1}, T, 1e5)` and compare. Which species did (b) leave out, and how much do they matter?

Answer form: two numbers; a number and three mole fractions; a table and a sentence.

## 2. (Derivation)

(a) Form the Lagrangian of min G/RT subject to A n = b and derive μ_j/RT = Σ_i a_ij π_i. Show that ∂(Σ n_k ln(n_k/n))/∂n_j = ln(n_j/n).
(b) From (a), derive the law of mass action for any atom-balanced reaction, and K_p = exp(−ΔG°/RT).
(c) Starting from the three linearised equations of RP-1311 (lesson §2.3), eliminate Δln n_j and derive the reduced system (6).
(d) Explain why G/RT is convex in n, and why that makes the answer independent of the starting guess.

Answer form: four derivations.

## 3. (Unfamiliar scenario: nitric oxide in a combustor)

Jet-A/air with 800 K air at 30 bar (typical, illustrative).

(a) Find the equilibrium flame temperature and the equilibrium NO mole fraction at the turbine-inlet mixture φ = 0.396 and at φ = 1.
(b) At φ = 0.396, find the equilibrium x_NO at 1500, 1700 and 1900 K. Why does it rise so steeply?
(c) K_p for ½N₂ + ½O₂ ⇌ NO at 1700 K, and a hand estimate of x_NO from x_N₂ and x_O₂.
(d) Real combustors make far less NO than these equilibrium values at the exit. Why, and is equilibrium an upper or a lower bound?

Answer form: two temperatures and two fractions; three fractions and a sentence; a constant and an estimate; a paragraph.

## 4. (Numerical, uses your code: your own rocket chamber)

H₂/O₂ at O/F = 4 + s/10 (mass), gaseous propellants at 298.15 K, chamber pressure p_c = (5 + s/3) MPa (illustrative).

(a) Find φ, the equilibrium chamber temperature, the mean molar mass and the composition at p_c, at p_c/10 and at 2 p_c, and the complete-combustion value.
(b) Sweep O/F from 3 to 8 at p_c and tabulate T_c, M̄ and √(T_c/M̄). The ideal exhaust velocity of a rocket grows roughly like √(T_c/M̄): which O/F maximises it, and why is it richer than stoichiometric (7.94)?

Answer form: a table for three pressures; a table and two sentences.

## 5. (Numerical, uses your code: testing the solver)

(a) For stoichiometric CH₄/air at 2200 K and 1 atm compare every mole fraction from `equilibrium_tp` and `equilibrium_scipy`, and their G/RT. Which species differ most in relative terms, and why does that not matter?
(b) Count Newton iterations for H₂ + ½O₂ at 3000 K and at 1000 K. Explain the difference using the conditioning of the element matrix in (6).
(c) Restart `equilibrium_tp` from three random log-amount vectors and show the answer does not change.

Answer form: a table and a sentence; two counts and an explanation; a check.

## 6. (Judgment: the species list)

(a) Recompute the stoichiometric H₂/O₂ flame at 1 bar with OH removed from the species list, and with only H₂O, H₂ and O₂. How much does each change T_ad?
(b) Which species missing from this module (for example HO₂, H₂O₂, NO₂, solid carbon) could matter in a combustor or a fuel-rich preburner, and in which direction would including them move the answer?
(c) A colleague wants to use this solver for the composition at the exit of a rocket nozzle. What would you tell them?

Answer form: two temperatures; a short list with directions; a paragraph.
