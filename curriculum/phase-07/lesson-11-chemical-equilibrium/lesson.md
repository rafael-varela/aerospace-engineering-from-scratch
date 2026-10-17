# P07.11  Chemical Equilibrium

> At flame temperatures the products are not finished molecules but a negotiation: every species takes the share that makes the Gibbs energy smallest, and the atoms are the only thing nobody can trade.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the equilibrium condition of a reacting ideal-gas mixture at fixed T and p as a constrained minimum of the Gibbs energy, identify the Lagrange multipliers as element potentials, and show that they imply the law of mass action with K_p = exp(−ΔG°/RT) for every reaction. 2. Implement the element-potential (reduced Gibbs iteration) method of NASA RP-1311 for H/O and C/H/O/N/Ar systems and test it against a general-purpose constrained minimiser and the law of mass action. 3. Compute equilibrium adiabatic flame temperatures and compare them with the complete-combustion values of `P07.10`. 4. Explain with numbers how temperature and pressure shift dissociation (Le Chatelier), and compute an illustrative rocket-chamber composition. 5. Decide when equilibrium, rather than complete combustion, must be used. |
| **Prerequisites** | `P03.15` (constrained optimization, Lagrange multipliers, `scipy.optimize.minimize`), `P07.10` (stoichiometry, enthalpies from the elements, complete-combustion T_ad), `P07.01` (NASA h and s°), `P07.04` (entropy and the second law), `P03.10` (Newton's method for nonlinear systems), `P03.02` (`brentq`), Bridge B2 or equivalent general chemistry (equilibrium constants). Skill check: you can write the Lagrangian of a function minimised under linear equality constraints and state its stationarity conditions. |
| **Engineering question** | How much hotter does complete combustion say a flame is than it really can be? For Jet-A/air with 800 K air at 30 bar (typical, illustrative combustor), and for hydrogen–oxygen, what are the equilibrium flame temperature and composition, and what does an illustrative rocket chamber (H₂/O₂ at O/F = 6, 20 MPa) contain? The answers: Jet-A/air peaks at 2612.9 K near φ = 1.05 (stoichiometric 2600.3 K, 154 K below complete combustion), but at the lean turbine-inlet mixture equilibrium is only 5.7 K cooler; H₂/O₂ at 1 bar reaches 3072.8 K, not 4930.6 K; the chamber reaches 3721.4 K with 66.0 % H₂O, 25.0 % H₂, 4.9 % OH and 3.4 % H by mole. |
| **Assumptions and limits** | Ideal-gas mixture; gas-phase species only (CO₂, H₂O, CO, H₂, O₂, N₂, Ar, OH, H, O, NO, N, CH₄; no soot, no liquid, no ions, no HO₂, H₂O₂, NO₂); chemical equilibrium reached (infinitely fast chemistry); adiabatic, steady, no work; NASA polynomials to 6000 K; fuels as gases at 298.15 K (the rocket's propellants would really enter as cryogenic liquids, whose enthalpies are not in this course's data, so its chamber temperature is illustrative and high). Fails: where chemistry is slower than the flow (NO formation, CO burnout, nozzle expansion: frozen or kinetic models), with condensed phases (soot, alumina), and at very low pressure or very high temperature (ionisation). **Computational only: there is no at-home combustion, flame, burning, fuel handling or fuel of any kind in this lesson.** |
| **Analytical work** | G = Σ n_j μ_j and its minimisation under A n = b; the Lagrange conditions μ_j/RT = Σ_i a_ij π_i; from them the law of mass action and K_p = exp(−ΔG°/RT); the reduced Newton system of RP-1311; the one-reaction dissociation estimate for water; Le Chatelier in T and p from ΔH° and Δν; the equilibrium flame temperature. |
| **Computational deliverable** | `starter/main.py`: `gibbs_rt`, `equilibrium_tp` (element-potential method), `equilibrium_scipy` (SLSQP check), `equilibrium_constant`, `formation_reactions`, `mass_action_residuals`, `equilibrium_flame_temperature`, `oxidizer_fuel_to_phi`; the NASA functions, `element_matrix`, `species_for`, the complete-combustion helpers from `P07.10` and `main` are given. Tests (24) prove Gibbs-energy identities, K composing like reactions, balanced formation reactions, element conservation to 10⁻¹² and the Lagrange conditions at four states, the same G from SLSQP at three states, no feasible perturbation lowering G, the law of mass action for every independent reaction at four states, independence from the starting guess, Le Chatelier in T and p, the low-temperature complete-combustion limit, the energy balance at the equilibrium flame temperature and its position below complete combustion, the pressure and φ trends, and the rocket mixture ratio. |
| **Practical activity** | **computational experiment** (computational only: no combustion, flame, burning, fuel handling or fuel of any kind): a method comparison (element potentials against SLSQP against the law of mass action), temperature and pressure sweeps of H₂/O₂ dissociation, equilibrium against complete-combustion flame temperatures for five flames and a φ sweep, and an illustrative rocket chamber at two pressures. |
| **Equipment** | None (a computer with the course's Python environment). There is no fuel of any kind in this lesson; nothing is burned. |
| **Safety** | None beyond normal computer use. This lesson is computational only: no at-home combustion, flame, burning, fuel handling or fuel of any kind. |
| **Evidence** | `equilibrium.png`; the method-comparison and mass-action tables; the T–p table of x_H₂O; the flame-temperature comparison; the chamber compositions; a half-page note answering the engineering question and stating when complete combustion is good enough. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the Lagrange conditions and the law of mass action derived; all public tests pass; the one-reaction water dissociation at 2800 K (or your seed) solved by hand; the note gives the equilibrium–complete gap for a lean and a stoichiometric flame with the reason they differ. |
| **References** | `@nasa-rp-1311-cea-analysis` Gordon & McBride, NASA RP-1311 Part I, Ch. 2 "Equations Describing Chemical Equilibrium": 2.2 Equation of State; 2.3 Minimization of Gibbs Energy (2.3.1 Gibbs Iteration Equations; 2.3.2 Reduced Gibbs Iteration Equations); Ch. 3 "Procedure for Obtaining Equilibrium Compositions": 3.1 Initial Estimates; 3.2 Magnitude of Species Used During Iteration; 3.3 Convergence; Ch. 6: 6.3.1 Combustion Conditions (required-free; https://ntrs.nasa.gov/citations/19950013764). `@libretexts-ellgen-thermodynamics-chemical-equilibrium` Ch. 13 "Equilibria in Reactions of Ideal Gases": 13.1, 13.2, 13.3, 13.6 The Standard Gibbs Free Energy Change and Equilibrium in Ideal Gas Reactions, 13.7 (supplementary-free; https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Thermodynamics_and_Chemical_Equilibrium_(Ellgen)/13%3A_Equilibria_in_Reactions_of_Ideal_Gases). `@openstax-chemistry-2e` 13.2 Equilibrium Constants; 13.3 Shifting Equilibria: Le Châtelier's Principle; 16.4 Free Energy (Free Energy and Equilibrium) (required-free; https://openstax.org/books/chemistry-2e/pages/13-2-equilibrium-constants). `@nasa-tp-2002-211556-glenn-coefficients` "Heats of Formation and Equilibrium Constants" (required-free; https://ntrs.nasa.gov/citations/20020085330). `@nasa-cea-github` the reference implementation, for comparison only (supplementary-free; https://github.com/nasa/cea). |
| **Extensions** | ★ Add the energy equation to the Newton system (RP-1311's hp problem, with Δln T as an unknown) and solve the flame temperature in one iteration loop instead of `brentq` around a T–p solver. ★ Compute frozen and equilibrium isentropic expansion of the chamber gas through a nozzle (fixed composition against re-equilibrating at each pressure) and compare the exit temperatures. ★ Add solid carbon as a condensed species (its chemical potential has no mixing term) and find the φ at which soot first appears in CH₄/air at 1500 K. |
| **Time** | unmeasured |

## 1. The engineering problem

Lesson `P07.10` gave a stoichiometric hydrogen–oxygen flame 4931 K and methane–oxygen 5166 K. No rocket chamber runs that hot, and the reason is not heat loss. At those temperatures water and carbon dioxide do not stay together: some H₂O breaks into H₂ and O₂, into OH and H, into H and O; some CO₂ into CO and O₂. Breaking bonds absorbs energy, so the gas settles at a lower temperature, with a mixture of fragments whose proportions depend on temperature and pressure. The same thing, less dramatically, caps the hottest zone of a jet-engine combustor and sets how much nitric oxide forms there.

Which mixture does the gas settle at? The second law answers: at fixed temperature and pressure a system reaches equilibrium where its Gibbs energy G = H − TS is smallest. So equilibrium is an optimization problem: minimise G over the amounts of every species, subject to conservation of every element. This lesson solves it two ways: with the element-potential method that NASA's Chemical Equilibrium with Applications (CEA) program uses (Gordon & McBride, RP-1311), and with a general-purpose constrained minimiser from `P03.15`. Then it checks both against the law of mass action, and wraps the solver in an energy balance to get the real (equilibrium) flame temperature.

**This lesson is computational only. There is no at-home combustion, flame, burning, fuel handling or fuel of any kind; all results are computed.**

Success criteria, written first: atoms conserved to 10⁻¹²; the second minimiser reaches the same G to 10⁻⁹ relative; ln Q = ln K to 10⁻⁶ for every independent reaction; at 1000 K the composition is complete combustion; the equilibrium flame temperature closes its energy balance and lies below the complete-combustion value for lean and stoichiometric mixtures.

## 2. Understand and derive

### 2.1 Equilibrium is a minimum of G

For a closed system at fixed T and p, the second law (`P07.04`) with dG = −S dT + V dp + Σ μ_j dn_j says any spontaneous change has dG ≤ 0: G falls until it cannot fall further. For an ideal-gas mixture of amounts n_j (total n), the chemical potential of species j is

μ_j = g°_j(T) + R T ln(n_j/n) + R T ln(p/p°),   (1)

with g°_j = h_j − T s°_j from the NASA polynomials at p° = 1 bar. Then G = Σ n_j μ_j, or in dimensionless form

**G/RT = Σ_j n_j [g°_j/RT + ln(n_j/n) + ln(p/p°)]**.   (2)

The amounts are not free: element i (C, H, O, N, Ar) has a fixed total b_i, and a_ij atoms of it sit in each molecule of j:

**Σ_j a_ij n_j = b_i**  (A n = b).   (3)

Equilibrium is the minimum of (2) subject to (3) and n_j ≥ 0. G/RT is convex in n (the ideal mixing terms make it so), so the minimum is unique and any method that finds a stationary point has found it.

### 2.2 Lagrange multipliers are element potentials

Form the Lagrangian L = G/RT − Σ_i π_i (Σ_j a_ij n_j − b_i). Since ∂(G/RT)/∂n_j = μ_j/RT (the ln n terms differentiate to zero in total), stationarity gives

**μ_j/RT = Σ_i a_ij π_i**  for every species j.   (4)

Read it as: each element has a chemical potential π_i (in units of RT), and a molecule's chemical potential is the sum of its atoms'. Water's is 2π_H + π_O; OH's is π_H + π_O. The code finds π_H = −11.4246 and π_O = −16.6928 for H₂ + ½O₂ at 3000 K and 1 bar.

**The law of mass action follows.** Take any reaction Σ ν_j X_j = 0 that balances atoms, Σ_j ν_j a_ij = 0. Multiply (4) by ν_j and sum: Σ ν_j μ_j = Σ_i π_i Σ_j ν_j a_ij = 0. Substituting (1),

Σ ν_j g°_j + RT Σ ν_j ln(x_j p/p°) = 0  ⇒  **Π (x_j p/p°)^ν_j = K_p = exp(−ΔG°/RT)**,   (5)

with ΔG° = Σ ν_j g°_j. Every reaction you can write among the species is in equilibrium at once; only (number of species − number of elements) of them are independent. Ellgen §13.6 derives (5) one reaction at a time; the element potentials give all of them at once.

**A one-reaction estimate by hand.** Water alone at 2800 K and 1 bar: H₂O ⇌ H₂ + ½O₂ with ΔG° = 89.12 kJ/mol, K_p = 0.02175. If a fraction α dissociates, the mole fractions are (1−α), α, α/2 over (1 + α/2), and (5) is α(α/2)^½/(1−α) · (p/p°/(1+α/2))^½ = K_p, solved by α = 0.0934: x_H₂O = 0.866. The full solver, which also allows OH, H and O, gives x_H₂O = 0.775: at 2800 K OH (6.4 %) matters as much as H₂, so the one-reaction picture is a first estimate only (Problem 1).

### 2.3 The element-potential method (RP-1311)

Newton's method on (3) and (4) in the variables ln n_j is efficient and keeps every n_j positive. RP-1311 §2.3 carries the total amount n as a separate unknown and writes the corrections Δln n_j, Δln n and the multipliers π_i (they appear linearly). Linearising:

Δln n_j − Σ_i a_ij π_i − Δln n = −μ_j/RT   (one per species)
Σ_j a_kj n_j Δln n_j = b_k − Σ_j a_kj n_j   (one per element)
Σ_j n_j Δln n_j − n Δln n = n − Σ_j n_j   (total amount)

The first line gives Δln n_j explicitly; substituting it into the other two leaves a small square system, one row per element plus one:

Σ_i [Σ_j a_kj a_ij n_j] π_i + b_k Δln n = b⁰_k − b_k + Σ_j a_kj n_j μ_j/RT,
Σ_i b_i π_i + (Σ_j n_j − n) Δln n = n − Σ_j n_j + Σ_j n_j μ_j/RT.   (6)

For C/H/O/N/Ar that is 6 × 6, whatever the number of species. Solve, recover Δln n_j, and update with a damping factor λ: λ = min(1, 2/max(5|Δln n|, |Δln n_j| of major species)), plus a limit that keeps a trace species (n_j/n < 10⁻⁸) from jumping above 10⁻⁴ in one step (RP-1311 §3.2–3.3). Start from n_j = 0.1 b_total/N_species each (RP-1311 §3.1 uses 0.1/N). Converge when every (n_j/n)|Δln n_j| and |Δln n| is below 10⁻¹² and every |Δln n_j| below 10⁻⁶.

One numerical wrinkle, found by the first version of the code: at 1000 K almost all hydrogen and oxygen sit in H₂O, the 2 × 2 element block of (6) becomes nearly singular (two elements, one molecule), and Δln n_j of the trace species stalls near 10⁻⁸ instead of reaching round-off. A criterion on |Δln n_j| alone never triggered; the mole-fraction-weighted criterion of RP-1311 does, and the law of mass action still holds to 10⁻¹⁴ there.

### 2.4 Le Chatelier, with numbers

Dissociation is endothermic (ΔH° > 0), so ln K grows with T (van 't Hoff, d ln K/dT = ΔH°/RT²): hotter means more fragments. It also makes more molecules (Δν > 0), so by (5) a higher pressure p/p° pushes it back. The code's x_H₂O for stoichiometric H₂/O₂:

| T (K) | 0.1 bar | 1 bar | 10 bar | 100 bar |
|---|---|---|---|---|
| 2000 | 0.9766 | 0.9892 | 0.9950 | 0.9977 |
| 2500 | 0.8022 | 0.9090 | 0.9581 | 0.9806 |
| 3000 | 0.3194 | 0.6391 | 0.8270 | 0.9192 |
| 3500 | 0.0227 | 0.2312 | 0.5597 | 0.7812 |

That is why rocket engines like high chamber pressure: less dissociation, more of the chemical energy turned into temperature.

### 2.5 The equilibrium flame temperature

The adiabatic, constant-pressure flame satisfies the same energy balance as in `P07.10`, but the products are now the equilibrium mixture at the unknown temperature:

**Σ_j n_j^eq(T, p) h_j(T) = H_reactants**.   (7)

Each trial T needs an equilibrium solve; `brentq` on T, warm-starting each solve from the last, converges in a few milliseconds. Because dissociation absorbs energy, the equilibrium T_ad is below the complete-combustion value for lean and stoichiometric mixtures. For rich mixtures `P07.10`'s fixed order of CO, H₂O and CO₂ is not the equilibrium split, and equilibrium can come out slightly hotter (Jet-A/air at 30 bar: 2465.9 K against 2462.8 K at φ = 1.3).

### 2.6 What generates entropy, and where the model stops

At equilibrium nothing generates entropy: G is at its minimum. The entropy of a combustor is generated in getting there, by the reaction itself (`P07.10` §2.4). The model stops where chemistry is too slow for the flow: NO forms over milliseconds (its equilibrium value is an upper bound in a combustor), CO burns out slowly at the combustor exit, and in a rocket nozzle the expanding gas cannot keep re-equilibrating (frozen flow, Extension 2). It also stops where condensed phases appear (soot) and where species it omits matter (HO₂ and H₂O₂ at moderate temperature and high pressure).

## 3. Build the computational model

Open `starter/main.py`. The NASA coefficients for CO₂, H₂O, CO, H₂, O₂, N₂, Ar, OH, H, O, NO, N, CH₄ and Jet-A(g), `h_species`, `s_species`, `g_species`, `molar_mass`, `element_totals`, `element_matrix`, `species_for`, the complete-combustion helpers of `P07.10` (`o2_stoich`, `reactants`, `reactant_enthalpy`, `complete_products`, `complete_flame_temperature`) and `main` are given.

### Step 1: Gibbs energy and K

`gibbs_rt(n, species, T, p)` is (2), skipping species with n_j = 0. `equilibrium_constant(reaction, T)` is exp(−Σ ν g°/RT).

### Step 2: the element-potential solver

`equilibrium_tp(b0, T, p, species=None, ln_n0=None)` follows 2.3. The heart of each iteration:

```python
mu = g0 + ln_n - ln_N + np.log(p / P_REF)          # mu_j / RT with the carried total N
M[:ne, :ne] = (A * nj) @ A.T;  M[:ne, ne] = bk;  M[ne, :ne] = bk;  M[ne, ne] = nj.sum() - N
rhs[:ne] = b - bk + A @ (nj * mu);  rhs[ne] = N - nj.sum() + np.sum(nj * mu)
pi, dlnN = np.split(np.linalg.solve(M, rhs), [ne]);  dln = -mu + A.T @ pi + dlnN
```

followed by the damping and the convergence test. Return the species, amounts, logarithms, mole fractions, π, the iteration count, G/RT and the element residual.

### Step 3: the independent checks

`equilibrium_scipy` minimises (2) with `scipy.optimize.minimize(method="SLSQP")`, the equality constraints A n = b, bounds n_j ≥ 10⁻³⁰ and the analytic gradient μ_j/RT; variables scaled by the total atoms. `formation_reactions(species, base)` writes one reaction per non-base species, solving A_base ν = a_j; `mass_action_residuals` returns (reaction, ln K, ln Q) for each, from the solver's ln n_j (so trace species at 10⁻¹⁴ are checked too).

### Step 4: the flame

`equilibrium_flame_temperature(fuel, phi, oxidizer, p, T_fuel, T_ox)` is (7) with `brentq` between 1000 K and the complete-combustion value plus 300 K, warm-started. `oxidizer_fuel_to_phi(fuel, of)` converts a rocket's mass mixture ratio: φ = (O/F)_st/(O/F), with (O/F)_st = 7.94 for H₂/O₂.

Running `main()` prints

```
H2 + O2/2 at 3000 K, 1 bar: 9 Newton iterations, element residual 1.4e-15
  mole fractions: H2O 0.6391, H2 0.1347, O2 0.0451, OH 0.0991, H 0.0580, O 0.0240
  element potentials pi_H = -11.4246, pi_O = -16.6928
  G/RT: element potentials -39.541969632401, SLSQP -39.541969632401 (difference -5.0e-14)
  law of mass action (base species H2, O2):
    1 H2 + 0.5 O2 -> 1 H2O       K_p = 2.2348e+01, ln K = +3.106738, ln Q = +3.106738
    0.5 H2 + 0.5 O2 -> 1 OH      K_p = 1.2715e+00, ln K = +0.240233, ln Q = +0.240233
    0.5 H2 -> 1 H                K_p = 1.5815e-01, ln K = -1.844199, ln Q = -1.844199
    0.5 O2 -> 1 O                K_p = 1.1315e-01, ln K = -2.179007, ln Q = -2.179007
```

then the tables of parts 2.4 and 4, and writes `equilibrium.png`: the equilibrium mole fractions of stoichiometric H₂/O₂ at 1 bar from 1500 to 4500 K; and T_ad against φ for Jet-A/air with 800 K air, complete combustion against equilibrium at 1 and 30 bar.

### Verify

Run `python -m pytest curriculum/phase-07/lesson-11-chemical-equilibrium`. The tests:

1. **Gibbs energy and K**: G/RT of a pure gas; mixing lowers G by n R ln 2 for equal parts; reversing a reaction inverts K; K of a sum of reactions is the product; dissociation's K grows with T; formation reactions balance atoms and number N_species − N_elements.
2. **Solver**: element residual below 10⁻¹² and the Lagrange conditions (4) to 10⁻⁶ for H/O at 3000 K, 1 bar and 1500 K, 200 bar, and C/H/O/N/Ar at 2200 K and 2500 K; SLSQP reaches the same G to 10⁻⁹ relative, never lower, with the same major-species mole fractions to 10⁻⁴; random feasible perturbations (along the null space of A) all raise G; mass action to 10⁻⁶ for every reaction with four different base sets; a random starting guess gives the same answer.
3. **Physics**: x_H₂O rises with p and falls with T; at 1000 K the composition is complete combustion (to 10⁻⁶ for H/O, 10⁻⁴ for CH₄/air).
4. **Flame**: the energy balance (7) closes; T_ad is below complete combustion for CH₄/air, H₂/O₂ and lean Jet-A/air; higher p gives a hotter H₂/O₂ flame by more than 100 K; the gap grows from φ = 0.4 to 0.7 to 1.0; the rocket's φ follows from O/F, its T_c lies between 3000 K and the complete value, and it holds more than 20 % H₂.
5. **main** writes the figure.

### Compare

Three independent routes reach the same answer: the element-potential solution, SLSQP (difference in G/RT 5 × 10⁻¹⁴ at 3000 K, 1.1 × 10⁻¹³ for CH₄/air at 2200 K) and the law of mass action (ln Q − ln K below 10⁻¹³ at every printed state). NASA CEA (`@nasa-cea-github`) implements the same method with a larger species list and condensed phases; this lesson's solver must not call it in tested code. If you run CEA yourself for CH₄/air at 1 atm and H₂/O₂ at 1 bar, compare with this module's 2224.3 K and 3072.8 K; small differences come from CEA's longer species list.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment. There is no at-home combustion, flame, burning, fuel handling or fuel of any kind in this lesson; every result below is computed.** Write each prediction before running it.

1. **Method comparison.** H₂ + ½O₂ at 3000 K, 1 bar: element potentials (9 Newton iterations, element residual 1.4 × 10⁻¹⁵) against SLSQP and the four mass-action reactions. *Finding:* the table above; identical to 13 digits.
2. **Dissociation against T and p.** The table of 2.4, plus the low-temperature limit: at 1000 K and 1 bar 1 − x_H₂O = 3.85 × 10⁻⁷, and CH₄/air at φ = 0.8 and 1000 K differs from complete combustion by at most 6.45 × 10⁻⁶ in mole fraction (NO 1.28 × 10⁻⁵). Predict the pressure at which 3000 K water is 90 % intact. *Finding:* between 10 and 100 bar.
3. **Flame temperatures, complete against equilibrium.** Predict which gap is largest. *Finding:*

| flame | complete (K) | equilibrium (K) | gap (K) |
|---|---|---|---|
| CH₄/air, φ = 1, 1 atm, 298 K | 2326.4 | 2224.3 | 102.1 |
| Jet-A/air, φ = 1, 1 atm, 298 K | 2418.1 | 2278.7 | 139.4 |
| Jet-A/air, φ = 1, 30 bar, 800 K air | 2754.7 | 2600.3 | 154.4 |
| Jet-A/air, φ = 0.396, 30 bar, 800 K air | 1700.0 | 1694.3 | 5.7 |
| H₂/O₂, φ = 1, 1 bar | 4930.6 | 3072.8 | 1857.7 |
| CH₄/O₂, φ = 1, 1 bar | 5166.5 | 3048.5 | 2118.0 |

4. **The φ sweep**, Jet-A/air with 800 K air at 30 bar: 2389.7, 2457.5, 2517.6, 2566.8, 2600.3, **2612.9**, 2603.5, 2578.0, 2543.8, 2505.7, 2465.9, 2425.5, 2385.2 K at φ = 0.80 to 1.40 in steps of 0.05. Predict where the peak moves. *Finding:* near φ = 1.05, slightly rich (at stoichiometric the products still hold spare O₂ and CO from dissociation; a little extra fuel uses the oxygen without much more dissociation), and the corner of `P07.10` is rounded off. At 1 bar the whole curve is lower (the figure).
5. **Illustrative rocket chamber**, H₂/O₂ at O/F = 6 (φ = 1.3228), gaseous propellants at 298.15 K. Predict the effect of going from 0.1 to 20 MPa. *Finding:* 3040.1 K and M̄ = 12.567 g/mol at 0.1 MPa (H₂O 0.5586, H₂ 0.2504, OH 0.0739, H 0.0894, O 0.0148, O₂ 0.0130); **3721.4 K and 13.432 g/mol at 20 MPa** (H₂O 0.6604, H₂ 0.2503, OH 0.0487, H 0.0336, O 0.0035, O₂ 0.0035), against 4208.2 K for complete combustion. SLSQP agrees on G/RT to 10 digits at the chamber state. Real engines inject liquid hydrogen and oxygen, whose lower enthalpy would lower T_c; this is a gas-propellant illustration.

Raw data format: CSV tables `species, x_ep, x_slsqp` and `reaction, lnK, lnQ` (1); `T_K, p_bar, x_H2O` (2); `flame, T_complete_K, T_eq_K` (3); `phi, T_eq_K, T_complete_K` (4); `p_c_MPa, T_c_K, M_gmol, x_...` (5); plus `equilibrium.png`.

## 5. Validate and improve

| check | route 1 | route 2 | difference |
|---|---|---|---|
| G/RT, H/O, 3000 K | element potentials −39.541969632401 | SLSQP −39.541969632401 | 5 × 10⁻¹⁴ |
| G/RT, CH₄/air, 2200 K | −350.1196695612 | −350.1196695612 | 1 × 10⁻¹³ |
| ln K vs ln Q, 4 reactions | ΔG° | composition | < 10⁻¹³ |
| element balance | b | A n | 1.4 × 10⁻¹⁵ relative |
| low-T limit | complete combustion | equilibrium at 1000 K | 4 × 10⁻⁷ |

The engineering answer: **complete combustion overstates a stoichiometric air flame by 100–155 K (Jet-A/air at 30 bar: 2754.7 → 2600.3 K, peaking at 2612.9 K near φ = 1.05) and an oxygen flame by about 1860–2120 K (H₂/O₂ at 1 bar: 3072.8 K), but at a lean turbine-inlet mixture the difference is 5.7 K.** Use equilibrium whenever the products are above about 2000 K: in primary zones, afterburners and every rocket chamber. The illustrative H₂/O₂ chamber at 20 MPa reaches 3721 K with a quarter of its moles still free H₂.

What to improve: liquid-propellant enthalpies for the rocket, and more species (HO₂, H₂O₂, condensed carbon). What not to improve: the convergence tolerance; three methods already agree to 10⁻¹³.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (4) and (5); the hand solution of the one-reaction water dissociation; `equilibrium.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the pressure is ten times higher?** Dissociation falls (Le Chatelier): stoichiometric H₂/O₂ at 3000 K goes from 63.9 % to 82.7 % H₂O, and the flame gets hotter (test: more than 100 K from 1 to 100 bar).
- **What if T falls to 1000 K?** Complete combustion is recovered to 4 × 10⁻⁷; the solver's matrix becomes ill-conditioned, which is why the convergence test is weighted.
- **Which assumption fails first?** Equilibrium itself, where chemistry is slower than the flow: NO in a combustor, and the expanding gas in a nozzle. Then the species list.
- **What measurement could falsify the model?** Spectroscopic temperature of a premixed flat flame at known φ and p: equilibrium should match it within heat losses; a flame measurably above the equilibrium T_ad would mean the enthalpy data or the species list are wrong.
- **Does the result survive a change of method?** Three independent routes agree to 10⁻¹³; dropping OH from the species list raises the H₂/O₂ flame from 3072.8 K to 3145.5 K, and keeping only H₂O, H₂ and O₂ raises it to 3507.4 K (Problem 6): the species list matters more than the solver.
- **Which decision changes if the uncertainty doubles?** A liner designed for 2755 K from complete combustion would be over-designed by 150 K; one designed for 2600 K with no margin would fail at φ = 1.05 (2613 K).

## Check yourself

1. Why is equilibrium at fixed T and p a minimum of G and not of H or of −S?
2. What does the element potential π_O mean physically, and why does a molecule's μ/RT equal the sum of its atoms' π?
3. How many independent reactions are there among the 13 species of a C/H/O/N/Ar mixture?
4. Write K_p for H₂O ⇌ H₂ + ½O₂ in terms of mole fractions and p. Does raising p increase or decrease dissociation?
5. Why does RP-1311 iterate on ln n_j rather than n_j?
6. Why does the equilibrium flame temperature fall below the complete-combustion one, and why is the gap tiny at φ = 0.4?
7. Why might equilibrium come out *hotter* than `P07.10`'s value for a rich mixture?
8. Name a situation in an engine where the equilibrium assumption fails, and say which way the error goes.
