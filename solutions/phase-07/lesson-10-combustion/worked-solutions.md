# P07.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers are from `solution.py` and the probe scripts that call it, given to the figures shown. All flame temperatures are complete-combustion values. Everything is computed: there is no combustion, flame or fuel of any kind.

## 1. Jet-A stoichiometry and the exhaust's dew point (s = 15: φ = 0.6)

**Assumptions.** Complete combustion; the phase's dry air; ideal gases.

(a) a = 12 + 23/4 = **17.75**. With air: C₁₂H₂₃ + 17.75 O₂ + 66.154 N₂ + 0.788 Ar + 0.034 CO₂ → 12.034 CO₂ + 11.5 H₂O + 66.154 N₂ + 0.788 Ar. Air mass per mole of fuel: 17.75/0.2095 × 28.9657 = 2454.1 g; **f_st = 167.311/2454.1 = 0.06818, AFR 14.67**.

(b) f = 0.6 × 0.06818 = **0.04091** (AFR 24.45). Per mole of fuel, reactants: Jet-A 1, O₂ 29.583, N₂ 110.256, Ar 1.313, CO₂ 0.056. Products: CO₂ 12.056, H₂O 11.5, O₂ 11.833, N₂ 110.256, Ar 1.313 (146.959 mol). Mole fractions: **CO₂ 0.0820, H₂O 0.0783, O₂ 0.0805, N₂ 0.7503, Ar 0.0089**.

(c) p_v = 0.07825 × 101 325 = **7929 Pa**; dew point **41.3 °C** (from the table). Exhaust carries a lot of water: when it mixes with air at −50 °C the mixture passes through saturation (over ice) on its way to ambient, and the water freezes into a contrail. The mixing line of `P07.09` is the tool for that question.

**Independent check.** Atom count: O on the right 2 × 12.034 + 11.5 = 35.568 = 2 × 17.75 + 2 × 0.034. ✓

**Common mistakes.** Forgetting the 0.0004 CO₂ in air (tiny, but the atom balance then fails at 10⁻³); using the mass fraction of water instead of the mole fraction for p_v.

## 2. Derivations

(a) C: x CO₂. H: y/2 H₂O. O: z + 2a = 2x + y/2 ⇒ **a = x + y/4 − z/2**.

(b) LHV and HHV differ only by the water term: (HHV − LHV) M_fuel = (y/2)[Hf(H₂O, g) − Hf(H₂O, l)] = (y/2) M_w h_fg(298.15 K). For methane: 2 × 44 004 J/mol/16.04246 g/mol = **5.486 MJ/kg**; with the steam table's 2441.68 kJ/kg, 2 × 18.01528/16.04246 × 2441.68 = 5.484 MJ/kg.

(c) Insert the products at T₀: H_p(T_ad) − H_p(T₀) = H_r(T₀) − H_p(T₀). The right side is the heat a reaction at T₀ would release, m_fuel × LHV (the spare O₂ and inerts appear on both sides at T₀ and cancel). The left side is m_p c̄_p (T_ad − T₀) by definition of c̄_p. Hence **T_ad − T₀ = m_fuel LHV/(m_p c̄_p)**. (For reactants not at T₀, add their sensible enthalpy to the numerator.)

(d) For ideal gases h = h(T) only, and the complete-combustion products do not depend on p, so neither side of H_p(T_ad) = H_r depends on p. The assumption is "complete combustion": with equilibrium products the composition depends on p (`P07.11`), and so does T_ad.

**Common mistakes.** Using the HHV in (c) with water vapour in the products.

## 3. A hydrogen-fuelled turbofan

(a) Hydrogen: **φ = 0.3518, f = 0.01026, 1.026 kg per 100 kg of air**. Jet-A: **φ = 0.3960, f = 0.02700, 2.700 kg per 100 kg of air**.

(b) f × LHV = **1230.6 kJ/kg air** (H₂) and **1170.4 kJ/kg air** (Jet-A). Both must heat the same air from 800 to 1700 K, so they are close. They differ because the products differ: hydrogen's products contain more water, whose cp is higher than CO₂-and-air at these temperatures, so it takes more heat to reach 1700 K; and the fuel enters at 298.15 K and must be heated too.

(c) Water **0.0917 kg/kg air** (H₂) against **0.0334 kg/kg air** (Jet-A): 2.7 times more. Hydrogen engines make more water per unit of thrust energy, so contrails form over a wider range of ambient conditions (the ice crystals may be fewer and larger, since there is no soot to nucleate on).

(d) Stoichiometric from an 800 K inlet: **2811.1 K** (H₂) and **2754.7 K** (Jet-A). Hydrogen's primary zone is hotter (and hydrogen burns faster), which raises the NOx and liner-cooling challenge.

**Common mistakes.** Comparing fuel mass flows without noting the factor of 2.8 in LHV.

## 4. Your own combustor (s = 15: T3 = 850 K, T4 = 1750 K)

(a) **φ = 0.4003, f = 0.02729**, fuel flow **1.364 kg/s** for 50 kg/s of core air.

(b) Stoichiometric T_ad = **2789.9 K**.

(c) dT4/dφ ≈ **1992 K** per unit φ (2.92 × 10⁴ K per unit f); 1 % more fuel (Δφ = 0.0040) adds **8.0 K**.

(d) Residual H_p(T4) − H_r = **−1.5 × 10⁻⁷ J** per mole of fuel, round-off.

**Independent check.** The lesson's case (T3 = 800, T4 = 1700) gave φ = 0.3960; adding 50 K to both moved φ up by 0.004, consistent with the products' cp exceeding the inlet air's.

## 5. The heat of reaction is not constant

(a) LHV/HHV (MJ/kg): CH₄ 50.027/55.513, H₂ 119.960/141.788, Jet-A 43.353/46.377; HHV − LHV = 5.4859, 21.8286, 3.0246 against 5.4839, 21.8205, 3.0235 from the table (+0.037 % each).

(b) Heat released with reactants and products at 1000 K: CH₄ **49.946 MJ/kg** (−0.08 MJ/kg, −0.16 %), H₂ **122.951 MJ/kg** (+2.99 MJ/kg, +2.5 %), Jet-A **43.655 MJ/kg** (+0.30 MJ/kg, +0.7 %).

(c) Kirchhoff: d(heat released)/dT = Σcp(reactants) − Σcp(products). For H₂ + ½O₂ the reactants' cp (44.8 J/K at 500 K) exceeds one H₂O's (35.2 J/K): fewer product molecules, less heat capacity, so the release grows with T. For CH₄ + 2O₂ → CO₂ + 2H₂O the two sides are close (108.8 vs 115.1 J/K at 500 K, 131.5 vs 128.9 J/K at 800 K), so the change is tiny and of mixed sign. Jet-A, with 17.75 O₂ against 23.5 product molecules plus a large fuel molecule, sits between. The heating value is a 298.15 K number; a first-law balance with h measured from the elements never needs it at any other temperature, which is why the module never uses it to compute T_ad.

## 6. Where complete combustion stops

(a) Turbine-inlet φ (1750 K): **use it**: below the dissociation range; `P07.11` shows the equilibrium correction at the lesson's turbine inlet is about 6 K. Stoichiometric Jet-A at 800 K: **only as an upper bound**, equilibrium is about 154 K lower at 30 bar. CH₄/O₂ at 5166 K: **do not use**: equilibrium at 1 bar is 3049 K; the products at 5166 K would be mostly atoms.

(b) `complete_products("CH4", 3, "O2")`: **CO 1.000, H₂ 1.667, H₂O 0.333** per mole of CH₄ at **1511 K**. Wrong: (1) at 1500 K and fuel-rich, equilibrium puts the oxygen differently (the water-gas balance, and methane itself survives at low temperature), so the split between CO, H₂ and H₂O is arbitrary here; (2) real fuel-rich methane combustion forms soot and leaves unburned CH₄ governed by kinetics, which no equilibrium or fixed-order model captures.

(c) A liquid fuel enters with less enthalpy (its heat of vaporization), so the same φ gives a slightly lower T4: the turbine-inlet φ must rise a little.

## Check-yourself answers

1. 17.75 mol of O₂; AFR 14.67 in the phase's dry air.
2. φ = f/f_st (actual over stoichiometric fuel–air ratio). Overall lean (φ ≈ 0.4) because the turbine cannot accept stoichiometric temperatures; the extra air cools.
3. Because every species' enthalpy is then measured from the same reference (the elements at 298.15 K), so the chemical energy is already in the numbers and a reaction's heat is just a difference of sums.
4. In the HHV the product water is liquid (its latent heat is recovered); in the LHV it is vapour. Engines exhaust vapour, so they are rated on the LHV.
5. Ideal-gas enthalpies depend on T only and complete products do not depend on p.
6. The products' cp rises strongly with T (vibrational modes of CO₂ and H₂O); a room-temperature cp is too small, so the temperature rise comes out too large.
7. Lean, the spare air dilutes; rich, there is not enough oxygen to release all the fuel's energy and the extra fuel dilutes. Both sides lower T_ad.
8. Most: the oxygen flames (CH₄/O₂, H₂/O₂, by about 2000 K) and stoichiometric air flames (100–150 K). Hardly at all: the lean turbine-inlet solve and the heating values.
