# P07.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers are from `solution.py` and the probe scripts that call it, given to the figures shown. Chamber and combustor conditions are typical, illustrative values; propellants are gases at 298.15 K. Everything is computed: there is no combustion, flame or fuel of any kind.

## 1. Water falling apart (s = 15: T = 2800 K)

**Assumptions.** Ideal gases, 1 bar, equilibrium.

(a) ΔG° = g_H₂ + ½ g_O₂ − g_H₂O = **89.12 kJ/mol** at 2800 K; K_p = exp(−89 122/(8.31446 × 2800)) = **0.02175**.

(b) With (1 − α) H₂O, α H₂, α/2 O₂ (total 1 + α/2) at p = p°: α(α/2)^½/[(1 − α)(1 + α/2)^½] = 0.02175 ⇒ **α = 0.0934**; x_H₂O = **0.866**, x_H₂ = **0.0892**, x_O₂ = **0.0446**.

(c) Full solver: H₂O **0.7755**, H₂ **0.0937**, O₂ **0.0324**, OH **0.0636**, H **0.0251**, O **0.0098**. The one-reaction model left out OH, H and O, which together hold 9.8 % of the moles; OH is as abundant as H₂. Water at 2800 K is 22 % dissociated (by mole fraction), not 13 %.

**Independent check.** The full solution satisfies H₂O ⇌ H₂ + ½O₂ too: 0.0937 × √0.0324/0.7755 = 0.0217 = K_p. Both satisfy the same constant; they differ because the full solution also uses oxygen and hydrogen atoms in OH, H and O.

**Common mistakes.** Forgetting the (1 + α/2) total in the mole fractions; using ΔG in J and R in kJ.

## 2. Derivations

(a) L = G/RT − Σ_i π_i(Σ_j a_ij n_j − b_i). With G/RT = Σ_k n_k (g°_k/RT + ln p/p°) + Σ_k n_k ln n_k − n ln n: ∂(n_j ln n_j)/∂n_j = ln n_j + 1 and ∂(n ln n)/∂n_j = ln n + 1, so ∂(Σ n_k ln(n_k/n))/∂n_j = **ln(n_j/n)**. Hence ∂L/∂n_j = μ_j/RT − Σ_i a_ij π_i = 0: **μ_j/RT = Σ_i a_ij π_i**.

(b) For a reaction with Σ_j ν_j a_ij = 0 for every element: Σ_j ν_j μ_j/RT = Σ_i π_i (Σ_j ν_j a_ij) = 0. With μ_j = g°_j + RT ln(x_j p/p°): Σ ν_j g°_j/RT + Σ ν_j ln(x_j p/p°) = 0, so **Π (x_j p/p°)^ν_j = exp(−ΔG°/RT) = K_p**.

(c) From the first equation Δln n_j = −μ_j/RT + Σ_i a_ij π_i + Δln n. Multiply by a_kj n_j and sum over j: Σ_j a_kj n_j Δln n_j = −Σ_j a_kj n_j μ_j/RT + Σ_i (Σ_j a_kj a_ij n_j) π_i + b_k Δln n, and set it equal to b⁰_k − b_k: the first row of (6). Multiply instead by n_j and sum: Σ_j n_j Δln n_j = −Σ n_j μ_j/RT + Σ_i b_i π_i + (Σ n_j) Δln n, and set it equal to n − Σ n_j + n Δln n: the second row.

(d) The terms linear in n are convex trivially; Σ n_k ln(n_k/n) is n times the negative ideal entropy of mixing per mole, a convex (positively homogeneous) function whose Hessian is positive semidefinite. With linear constraints the feasible set is convex, so any stationary point is the global minimum and in the interior (n_j > 0) it is unique: every start converges to the same answer.

## 3. Nitric oxide in a combustor

(a) φ = 0.396: T_ad = **1694.3 K**, x_NO = **0.00212**. φ = 1: T_ad = **2600.3 K**, x_NO = **0.00447**.

(b) φ = 0.396 at fixed T: x_NO = **0.00092, 0.00217, 0.00427** at 1500, 1700, 1900 K. NO formation is strongly endothermic (about 90 kJ/mol), so K_p rises steeply with T: 0.00301, 0.00713, 0.01410 at the same temperatures.

(c) K_p(1700 K) = **0.00713**. With x_N₂ = 0.7593 and x_O₂ = 0.1221 from the solve, x_NO ≈ K_p √(x_N₂ x_O₂) = 0.00713 × √0.0927 = **0.00217**, matching the solver to every printed digit (Δν = 0, so pressure drops out).

(d) NO forms by slow chain reactions (the Zeldovich mechanism, initiated by O atoms attacking N₂) with time constants of milliseconds to seconds at these temperatures, comparable to or longer than the combustor residence time. So the exit holds less NO than equilibrium: equilibrium is an **upper bound** at the exit (and, at the hot primary zone, the measure of what slow chemistry is heading towards). It also explains why designers fight peak temperature and residence time to cut NOx.

## 4. Your own rocket chamber (s = 15: O/F = 5.5, p_c = 10 MPa)

(a) φ = 7.937/5.5 = **1.4430**.

| p (MPa) | T_c (K) | M̄ (g/mol) | H₂O | H₂ | OH | H | O | O₂ |
|---|---|---|---|---|---|---|---|---|
| 1 | 3290.8 | 12.218 | 0.5786 | 0.2968 | 0.0525 | 0.0616 | 0.0060 | 0.0045 |
| 10 | **3561.1** | **12.545** | 0.6180 | 0.3014 | 0.0388 | 0.0372 | 0.0026 | 0.0020 |
| 20 | 3634.2 | 12.636 | 0.6296 | 0.3024 | 0.0339 | 0.0308 | 0.0019 | 0.0015 |

Complete combustion: **3997.1 K**.

(b) At 10 MPa:

| O/F | T_c (K) | M̄ (g/mol) | √(T_c/M̄) |
|---|---|---|---|
| 3 | 2699.0 | 8.042 | 579.3 |
| 4 | 3155.6 | 9.954 | 563.0 |
| 5 | 3458.2 | 11.724 | 543.1 |
| 6 | 3636.1 | 13.318 | 522.5 |
| 7 | 3716.5 | 14.718 | 502.5 |
| 8 | 3730.1 | 15.928 | 483.9 |

√(T_c/M̄) keeps rising as the mixture gets richer over this range: excess H₂ lowers M̄ faster than it lowers T_c. (Real engines run near O/F 5–6 rather than 3 because exhaust velocity is not the only criterion: liquid hydrogen is bulky, so tank mass and volume push toward more oxygen; this simple index also ignores γ.)

**Common mistakes.** Using the stoichiometric O/F of 8 instead of 7.937; forgetting that p_c changes T_c only through dissociation.

## 5. Testing the solver

(a) CH₄/air, 2200 K, 1 atm: G/RT = −350.1196695612 by both; difference 1.1 × 10⁻¹³. Major species agree to about 10⁻⁸ (CO₂ 0.0864838 against 0.0864838, H₂O 0.1835604, N₂ 0.7012120 …). The largest relative differences are in the trace species: N (1.0357 × 10⁻⁸ against 1.0361 × 10⁻⁸, 0.03 %) and CH₄ (2.37 × 10⁻¹⁷ against 3.11 × 10⁻¹⁷). SLSQP's tolerance is on G, and a trace species barely changes G, so it is determined loosely; the element-potential solution, which satisfies (4) for every species, is the one to trust for trace species.

(b) **9** iterations at 3000 K, **27** at 1000 K. At 1000 K nearly all H and O sit in H₂O, so the element block Σ_j a_kj a_ij n_j is almost the rank-one matrix n_H₂O (2,1)(2,1)ᵀ: nearly singular. Newton steps for the trace species are then poorly determined and the damping keeps them small, so more iterations are needed.

(c) Three random starts drawn from ln n_j ∈ [−12, 1] (seeded `default_rng(1)`) give mole fractions within **1.6 × 10⁻¹⁵** of the default start.

## 6. The species list

(a) Stoichiometric H₂/O₂, 1 bar: full list **3072.8 K**; without OH **3145.5 K** (+73 K); with only H₂O, H₂, O₂ **3507.4 K** (+435 K). Each species left out is a way the mixture could have absorbed energy; leaving it out makes the flame hotter.

(b) HO₂ and H₂O₂: matter at moderate temperature and high pressure (fuel-rich preburners, ignition), little at flame peaks; including them would lower T slightly. NO₂: small at flame temperatures. Solid carbon: decisive in fuel-rich hydrocarbon mixtures (it forms when O is scarce), and changes the gas composition completely; including it moves the rich-side answer substantially. Ions: only at much higher temperatures.

(c) The chamber composition is a fair equilibrium answer. In the nozzle the gas expands in milliseconds: recombination (H + OH → H₂O and so on) cannot keep up, especially downstream of the throat. The truth lies between "frozen" (composition fixed at the chamber value) and "shifting equilibrium" (re-equilibrated at every pressure). Run both (Extension 2) and quote the range; for performance, frozen flow gives the lower bound.

## Check-yourself answers

1. At fixed T and p the second law for a closed system becomes dG ≤ 0 (the surroundings' entropy change is −dH/T, so the total entropy rises exactly when G falls); H and S separately are not the right potentials when heat is exchanged at fixed T.
2. π_O is the chemical potential of oxygen atoms (divided by RT): the change in G/RT per extra mole of O atoms made available, whatever molecule they end up in. Because at equilibrium moving atoms between molecules cannot lower G, each molecule's μ must be the sum of its atoms'.
3. 13 − 5 = **8** independent reactions.
4. K_p = x_H₂ x_O₂^½/x_H₂O · (p/p°)^½. Raising p makes the factor (p/p°)^½ larger, so the mole-fraction ratio must fall: less dissociation.
5. To keep n_j positive automatically and to handle amounts spanning 30 orders of magnitude with steps of the same scale.
6. Dissociation absorbs part of the released energy, so less goes into temperature. At φ = 0.4 the gas reaches only 1700 K, where K for dissociation is tiny (5.7 K difference).
7. For rich mixtures `P07.10` splits the products by a fixed rule; equilibrium chooses a lower-G split of CO, CO₂, H₂ and H₂O that may also release a little more energy.
8. In a nozzle (frozen chemistry: equilibrium overestimates recombination heat and so exhaust velocity) or NO in a combustor (equilibrium overestimates exit NO).
