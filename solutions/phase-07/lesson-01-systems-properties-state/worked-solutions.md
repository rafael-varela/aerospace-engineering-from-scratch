# P07.01  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numbers come from `solution.py` and the probe scripts that call it. R_U = 8.314 462 618 J/(mol·K); p° = 1 bar.

## 1. Nitrogen at 600 K and 2 bar (s = 15)

**Assumptions.** Ideal gas; NASA N₂ fit, lower interval (600 K < 1000 K); M = 28.0134 g/mol.

(a) cp = **30.109 J/(mol·K)** = **1074.81 J/(kg·K)**; cp/R = **3.6213**; γ = cp/(cp − R_U) = **1.3815**.

(b) h = **8893.96 J/mol** = **317.49 kJ/kg** (above the 298.15 K datum where h = 0 for an element); s° = **212.176 J/(mol·K)**; s(2 bar) = s° − R_U ln 2 = **206.413 J/(mol·K)** = **7368.36 J/(kg·K)**.

(c) v = R T/p = (8.314 4626/0.028 0134) × 600/2 × 10⁵ = **0.89041 m³/kg**; ρ = **1.1231 kg/m³**.

(d) [h(610) − h(590)]/20 = 30.10925 J/(mol·K), **6.5 × 10⁻⁶** above cp: the central difference's truncation error, (ΔT²/6) times the second derivative of cp.

**Check.** cp/R = 3.62 lies between 7/2 (vibration frozen) and 9/2 (fully excited): at 600 K the N₂ vibration (characteristic temperature 3390 K, `P04.07`) is partly active. γ = 1.38 matches `P04.07`'s Einstein-model estimate of 1.383 at 600 K.

**Common mistakes.** Mixing per-mole and per-kilogram values; using ln(p) instead of ln(p/p°); applying the pressure correction to h (an ideal gas's h does not depend on p).

## 2. Derivations

(a) h = ∫ cp dT. With cp/R = Σ aₖ Tᵏ over k = −2…4: ∫a₋₂T⁻² dT = −a₁T⁻¹; ∫a₋₁T⁻¹ dT = a₂ ln T; ∫a₀ dT = a₃T; ∫aₖTᵏ dT = aₖTᵏ⁺¹/(k + 1) for k = 1…4. Adding the constant b₁R and dividing by RT gives **h/(RT) = −a₁T⁻² + a₂ ln(T)/T + a₃ + a₄T/2 + a₅T²/3 + a₆T³/4 + a₇T⁴/5 + b₁/T**. For s° = ∫ cp/T dT: ∫a₁T⁻³ dT = −a₁T⁻²/2; ∫a₂T⁻² dT = −a₂T⁻¹; ∫a₃T⁻¹ dT = a₃ ln T; ∫aₖTᵏ⁻¹ dT = aₖTᵏ/k for k = 1…4; plus b₂: **s°/R = −a₁T⁻²/2 − a₂T⁻¹ + a₃ ln T + a₄T + a₅T²/2 + a₆T³/3 + a₇T⁴/4 + b₂**. b₁ sets the enthalpy datum, b₂ the entropy datum (third law).

(b) For Ar, h/(RT) = 2.5 + b₁/T. h(298.15) = 0 requires b₁ = −2.5 × 298.15 = **−745.375**, exactly the tabulated b₁.

(c) h = u + p v = u + R T for an ideal gas, and u = u(T) only. Differentiating, dh/dT = du/dT + R, i.e. **cp = cv + R**. Since R > 0 and cv > 0, **γ = cp/cv = 1 + R/cv > 1**.

(d) In Dalton's model species i occupies the whole volume at partial pressure p_i = x_i p, so its molar entropy is s°_i(T) − R_U ln(p_i/p°). Entropy is extensive, so the mixture's molar entropy is the mole-fraction average: **s = Σ x_i [s°_i − R_U ln(x_i p/p°)] = Σ x_i s°_i − R_U ln(p/p°) − R_U Σ x_i ln x_i**. For dry air, −R_U Σ x_i ln x_i = **4.7167 J/(mol·K)** = **162.84 J/(kg·K)**, the value the lesson quotes.

**Common mistakes.** Forgetting that the T⁻¹ term of cp integrates to a logarithm (the a₂ term); dividing the s° integrand by T twice; omitting the −R_U Σ x_i ln x_i term and then wondering why air's s is not the average of its species'.

## 3. Martian air

**Assumptions.** Ideal-gas CO₂, NASA fit (lower interval, valid from 200 K).

(a) At 210 K: cp/R = **3.9497**, γ = **1.3390**. At 300 K: cp/R = **4.4765**, γ = **1.2876**.

(b) f = 2/(γ − 1) = **5.90** at 210 K and **6.95** at 300 K. Translation and rotation give 5; the excess (0.9 at 210 K, about 2 at 300 K) is CO₂'s doubly degenerate bending vibration, soft enough to be partly excited even at Martian temperatures.

(c) γ ≈ 1.29 belongs to **300 K**, not 210 K. At 210 K the ideal-gas value is 1.339, so a speed-of-sound estimate with 1.29 is low by √(1.29/1.339) − 1 = **−1.9 %** (the speed of sound itself is about 1.9 % higher than that estimate). A small error, but a systematic one, and exactly the kind of number that should come from a property fit, not from memory.

(d) H₂ at 300 K: cp/R = **3.4697**, γ = **1.4049**. Light as it is, H₂ is a diatomic molecule with rotation active (5 quadratic terms) and a very stiff bond, whose vibration is frozen; cp/R is a little below 7/2 because hydrogen's rotational levels are so widely spaced that rotation is not quite fully classical at 300 K. Molecular mass sets the speeds of `P04.07`, not the number of active degrees of freedom.

**Check.** At 300 K the NASA γ of CO₂ (1.288) is the textbook room-temperature value; the 210 K value is higher, as it must be, because fewer vibrational states are populated.

**Common mistakes.** Treating γ as a property of the molecule rather than of the molecule at a temperature; using the speed-of-sound ratio γ instead of √γ.

## 4. A combustor (s = 15: q = 900 kJ/kg)

**Assumptions.** Steady flow, constant pressure, q = h₂ − h₁ (lesson 02), air composition unchanged (in reality the gas becomes combustion products).

(a) Variable cp: **T₂ = 1570.77 K**. Constant cp: 800 + 900 000/1005 = **1695.52 K**. Difference **124.75 K**.

(b) Mean cp = 900 000/(1570.77 − 800) = **1167.67 J/(kg·K)**.

(c) h₂ = h(800) + q = 518.36 + 900.00 = 1418.36 kJ/kg (NASA datum). `air_state(p=2e6, h=1.41836e6)` gives T = 1570.77 K, **v = 0.22544 m³/kg**, **u = 967.48 kJ/kg**, **s = 7810.18 J/(kg·K)**.

(d) The constant-cp designer predicts 1695.5 K, above the 1650 K rating, and would throttle the engine back or add cooling air. The variable-cp answer, 1570.8 K, is 79 K inside the rating. Here constant cp errs on the conservative side, but at the price of performance; in Problem 5 and in the turbine (lesson 06) the same error can be unconservative.

**Check.** The mean cp (1167.7) lies between cp(800 K) = 1099 and cp(1570 K) ≈ 1218, as it must.

**Common mistakes.** Using constant cp to "check" the variable-cp answer; forgetting the NASA datum when computing h₂ (any datum works, as long as it is the same one for h₁ and h₂).

## 5. An isentropic compressor (s = 15: π = 15)

**Assumptions.** Reversible adiabatic compression (s₂ = s₁), ideal gas, variable cp.

(a) s₁ = air_s(288.15, 101 325) = 6826.23 J/(kg·K). With p₂ = 15 × 101 325 Pa and s₂ = s₁: **T₂ = 617.03 K**.

(b) Constant γ = 1.4: 288.15 × 15^0.2857 = **624.66 K**, **7.63 K** higher, because the true γ falls from 1.400 at 288 K to 1.374 at 617 K as the air heats.

(c) Variable cp: h₂ − h₁ = **336.78 kJ/kg**. Constant cp and γ: 1005 × (624.66 − 288.15) = **338.19 kJ/kg**, 0.4 % high. The two errors (too high a T₂, too low a cp) partly cancel.

(d) **p and s**: independent, so the state postulate says they fix the state. (T₂, h₂) could not have: both are functions of T alone for an ideal gas, so they fix T but not p; and you did not know T₂ in advance anyway.

**Check.** At constant p, s increases with T, and at constant T, s decreases with p; compression at constant s therefore must raise T. Also the code's (p, s) state reproduces T₁ when p = p₁.

**Common mistakes.** Using the constant-γ formula and calling it "isentropic" without saying it assumes constant cp; passing a NumPy array where a float is expected.

## 6. Judgment

(a) **Adequate.** Cabin air at 0.8 bar and 295 K is far from condensation, so the ideal-gas law is excellent, and it is in the range where cp is flat (cp(295 K) = 1004.6 J/(kg·K), so 1005 is right to 0.1 %); only its humidity, which changes M and cp by under 1 % (Experiment 3), might need lesson 09's moist-air model.

(b) **Composition fails first.** At 1900 K the gas is combustion products (N₂, H₂O, CO₂, excess O₂), with a larger cp than air, and dissociation is beginning; use the species fits in a product mixture (lesson 10) and equilibrium composition if accuracy matters (lesson 11).

(c) **The ideal gas fails first.** At 30 MPa nitrogen's molecules are close enough that Z departs noticeably from 1, which changes the mass the tank holds for a given p and T; use a real-gas equation of state or tables (lesson 03's method) for tank sizing, and keep the ideal gas only for trend studies (lesson 02).

(d) **Fixed composition fails.** Several thousand kelvin dissociates O₂ and then N₂; the gas is a reacting mixture whose composition depends on p and T, so you need chemical equilibrium (lesson 11) or non-equilibrium chemistry; the species fits remain valid, the mixture does not.

(e) **Adequate, but check Z once.** At 10 bar and 300 K air is far above its condensation temperatures, so the ideal gas is a good model, though at this pressure it is worth checking Z against a real-air table once rather than assuming it; the test-section speed of sound depends on the static temperature after expansion, where the gas is even more ideal; if the expansion is cold enough for oxygen to approach condensation, that changes, which is why hypersonic tunnels heat their air.

## Check-yourself answers

1. A **control volume**: propellants flow in and exhaust flows out. A closed system would follow one parcel of gas through it, which is harder.
2. Extensive: V, U, m. Intensive: v, T, h, p.
3. Two independent intensive properties fix the state of a simple compressible substance. For an ideal gas h depends on T only, so (T, h) is one property counted twice.
4. h = u + RT; dh/dT = du/dT + R; cp = cv + R.
5. Vibration of the N–N bond is progressively excited (`P04.07`); cp/R rises from 7/2 towards 9/2 (4.45 at 3000 K).
6. b₁ is the integration constant chosen so that h(298.15 K) equals the enthalpy of formation. Air's h is negative at 288 K because the datum is 298.15 K (N₂, O₂, Ar have h = 0 there, so they are slightly negative below it) and its 0.04 % of CO₂ carries CO₂'s large negative enthalpy of formation.
7. s = s°(T) − R ln(p/p°); doubling p lowers s by R ln 2 (199.0 J/(kg·K) for air).
8. About 620 K (618 K for the Phase 07 dry air).
