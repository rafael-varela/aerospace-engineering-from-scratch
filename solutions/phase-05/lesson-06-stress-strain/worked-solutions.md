# P05.06  Worked solutions

Each solution states assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Every tensile record here is **synthetic**.

## 1. A control-rod end in tension (s = 15: d = 7.5 mm)

**Assumptions.** Uniform axial stress away from the ends; typical 7075-T6 values (E = 71.7 GPa, ν = 0.33, Fty = 503 MPa); elastic.

(a) A = π(7.5)²/4 = 44.179 mm². σ = 10 000/44.179 = **226.4 MPa**. ε = 226.4/71 700 = **3.157 × 10⁻³**. ΔL = 3.157 × 10⁻³ × 300 = **0.947 mm**. Δd = −0.33 × 3.157 × 10⁻³ × 7.5 mm = **−7.81 μm**.

(b) σ_t = σ(1 + ε) = **227.07 MPa**, ε_t = ln(1 + ε) = **3.152 × 10⁻³**; ratio **1.0032**. A 0.3 % difference is smaller than the uncertainty of a typical modulus (about 1–2 %); in the elastic range the distinction is not worth keeping.

(c) F_y = 503 × 44.179 = **22.22 kN**.

**Check.** Units: N/mm² = MPa. The yield load over the applied load, 2.22, equals Fty/σ = 503/226.4.

**Common mistakes.** Using the diameter instead of the radius in πr²; reporting the diameter change as −νε without multiplying by d; quoting true stress as if it mattered at 0.3 % strain.

## 2. True measures, Considère and the elastic constants

(a) Constant volume: A L = A₀ L₀, so A = A₀ L₀/L = A₀/(1 + ε). Then **σ_t = F/A = (F/A₀)(1 + ε) = σ(1 + ε)**. The true strain sums the increments relative to the current length: **ε_t = ∫_{L₀}^{L} dL/L = ln(L/L₀) = ln(1 + ε)**. Both need the deformation to be uniform along the gauge length; after the maximum load a neck forms and the area there is not A₀/(1 + ε).

(b) F = σ_t A. dF = A dσ_t + σ_t dA, and dA/A = −dL/L = −dε_t at constant volume, so dF = A(dσ_t − σ_t dε_t), which is zero when **dσ_t/dε_t = σ_t**. For σ_t = Kε_tᵐ: mKε_t^(m−1) = Kε_tᵐ ⇒ **ε_t = m**. For the lesson's Ramberg–Osgood material (n = 9), `considere_point` gives ε_u = 0.1249, true **ln(1.1249) = 0.1177**, against **1/n = 0.1111**. The plastic part alone is a power law σ_t = σ_y(ε_p/0.002)^(1/n), for which Considère gives ε_p = 1/n; the total true strain also contains the elastic strain σ_t/E (about 0.007 at the 539 MPa true stress of the maximum), and the criterion is applied to the total, which moves the point.

(c) In (6) the plastic strain is 0.002(σ/σ_y)ⁿ, equal to 0.002 at σ = σ_y. The 0.2 % offset line is σ = E(ε − 0.002), which meets (6) where σ/E + 0.002(σ/σ_y)ⁿ − 0.002 = σ/E, i.e. at σ = σ_y. For offset o: 0.002(σ/σ_y)ⁿ = o ⇒ **σ = σ_y (o/0.002)^(1/n)**. For the lesson's material the 0.05 % proof stress is 345 × 0.25^(1/9) = **295.7 MPa**; the test `test_vanishing_offset_and_monotone_order_of_points` checks this formula.

(d) Principal stresses +τ along one diagonal and −τ along the other. Strain along the +τ diagonal: ε_d = τ/E − ν(−τ)/E = (1 + ν)τ/E; along the other, −(1 + ν)τ/E. A square whose diagonals stretch by +ε_d and −ε_d becomes a rhombus whose corner angle changes by γ = 2ε_d (for small strains: the diagonal at 45° stretches by ε_d = γ/2). So γ = 2(1 + ν)τ/E and **τ/γ = G = E/2(1 + ν)**.

**Check.** 2024-T3: 73.1/(2 × 1.33) = 27.5 GPa, the value in the lesson's table. `test_pure_shear_and_hydrostatic_pressure_give_G_and_K` checks the same identity through (5).

**Common mistakes.** Writing ε_t = ln(ε); applying the true-stress formulas after the maximum; confusing the tensor shear strain ε₁₂ with the engineering shear strain γ = 2ε₁₂ (that factor of 2 is exactly the "2" in G = E/2(1 + ν)).

## 3. A titanium pressurant sphere (s = 15: p = 2.0 MPa)

**Assumptions.** Thin wall (t/R = 0.004); membrane stresses; radial stress (between −p and 0) neglected against 250 MPa; typical Ti-6Al-4V (E = 113.8 GPa, ν = 0.342, Fty = 880 MPa).

(a) Cut through the centre. The pressure acts on the projected area πR² and is balanced by the wall stress on the ring 2πRt: p πR² = σ 2πRt ⇒ **σ = pR/2t = 2.0 × 300/2.4 = 250.0 MPa**, the same in every tangential direction by symmetry.

(b) `hooke_strain(np.diag([250, 250, 0]), 113.8e3, 0.342)`: tangential strain **1.4455 × 10⁻³** (= σ(1 − ν)/E); thickness strain **−1.5026 × 10⁻³** (= −2νσ/E). ΔR = εR = **0.434 mm**; Δt = **−1.803 μm**; ΔV/V = 3ε = **4.337 × 10⁻³**; with V = 113.10 L, **ΔV = 490.5 cm³**.

(c) Ignoring the Poisson effect: ΔR = σR/E = **0.659 mm**, an overestimate by 1/(1 − ν) − 1 = **52 %**.

(d) σ/Fty = **0.284**. A real tank is sized by a proof factor and a burst factor on the maximum expected operating pressure, by fatigue from pressure cycles, by the knocked-down strength of welds, and by temperature; the membrane stress alone is the start of that calculation.

**Check.** The tangential strain by hand: σ(1 − ν)/E = 250 × 0.658/113 800 = 1.4455 × 10⁻³, as `hooke_strain` gives; and the thickness strain is −2ν/(1 − ν) = −1.04 times it, as the two printed strains show. (ΔV/V = 3ε is the small-strain linearisation; the exact (1 + ε)³ − 1 is 0.14 % larger.)

**Common mistakes.** Using the cylinder's hoop formula pR/t (twice too large) for a sphere; forgetting that both tangential stresses reduce each other's strain through ν.

## 4. Reduce a titanium coupon (synthetic, s = 15)

(a) `analyze_tensile(..., window=21)` on seed 15:

| property | value |
|---|---|
| E | 112.37 GPa |
| Fty (0.2 % offset) | 874.1 MPa |
| Ftu | 974.6 MPa |
| ε_u | 0.0572 |
| elongation after fracture | 0.1322 |
| toughness | 128.4 MJ/m³ |
| resilience | 3.400 MJ/m³ |

(b) Inputs and independent values: E = 113.8 GPa; σ_y = 880 MPa (true-stress law); Considère **ε_u = 0.0604, σ_u = 974.4 MPa**. The clean record (no noise) reduces to E = 113.06 GPa, Fty = 869.8 MPa, Ftu = 974.4 MPa at 0.0604, elongation 0.1323, toughness 128.44 MJ/m³. Differences: Fty clean is 1.2 % below 880 MPa because the engineering stress is σ_t/(1 + ε) with ε ≈ 0.0097 at yield; the noisy-record Fty (874.1) is above the clean one by 4.3 MPa, about two standard deviations of the noise study (d); E is 0.6 % low (noise in strain, see (d)); Ftu agrees within the noise; **ε_u moves from 0.0604 to 0.0572**: the top of the curve is so flat that the *location* of the maximum is poorly determined by noisy data, even when its *value* is well determined. Compared with the typical Ftu of 950 MPa, the synthetic material is 2.6 % strong: the generator's n = 20 was not calibrated to it.

(c) Windows 5–20 %, 10–40 %, 20–60 %: E = 113.37, 112.37, 112.63 GPa; Fty = 869.9, 874.1, 873.7 MPa. The spread (0.9 %) is noise, not curvature: with n = 20, the plastic strain at 60 % of the maximum is 0.002 × (0.664)²⁰ ≈ 6 × 10⁻⁷, against 2024-T3 with n = 9 where it is large enough to pull the 20–60 % modulus down by 10 %. A sharp knee makes the elastic range clean.

(d) 100 seeds: **E = 113.02 ± 0.64 GPa, Fty = 870.1 ± 2.08 MPa**. Seed 15 is −1.0σ in E and +1.9σ in Fty: **yes**, within 2σ, though only just in Fty.

**Common mistakes.** Comparing Fty of the engineering record with the true-curve σ_y and calling the difference an analyzer error; reporting ε_u to four figures from a noisy record.

## 5. Diagnose a crosshead-strain record (synthetic, s = 15: c = 6 × 10⁻⁶ per MPa)

(a) `analyze_tensile` gives **E_app = 50.52 GPa**. No aluminium alloy has a modulus near 50 GPa (they lie around 69–73 GPa); a value this low says that the strain was not measured on the gauge length.

(b) c_est = 1/50.52 GPa − 1/72.49 GPa = **6.000 × 10⁻⁶ per MPa**, the value used. (Against the handbook 73.1 GPa it would be 6.115 × 10⁻⁶: the 0.8 % difference between the clean-record and handbook moduli is a 2 % error in c.)

(c) Recorded strain ε_rec = σ/E + cσ + ε_p = σ/E_app + ε_p. The offset line with the fitted modulus is σ = E_app(ε_rec − ε₀ − 0.002). Substituting, σ = σ + E_app(ε_p − ε₀ − 0.002), so the crossing is at ε_p = 0.002 + ε₀ regardless of c. Numerically **Fty = 342.3 MPa**, the same as the clean record and as c = 0, 2, 5, 10 × 10⁻⁶ in the lesson.

(d) With E = 73.1 GPa: **Fty = 278.2 MPa**, **64.1 MPa (18.7 %) too low**. The handbook line is steeper than the record's elastic slope, so it meets the curve too early.

**Check.** The lesson's compliance table brackets this result: 294.0 MPa at c = 5 × 10⁻⁶ and 196.8 MPa at 10 × 10⁻⁶.

**Common mistakes.** "Correcting" the modulus by replacing it with the handbook value without correcting the strain; concluding that the material is weak.

## 6. What the numbers can and cannot tell you

(a) Ask (1) where the strain was measured (extensometer on the gauge length, or crosshead travel) and (2) which modulus was used to draw the offset line. Distrust **E = 52 GPa first**: it is physically impossible for an aluminium alloy and is the signature of compliance. If the offset line then used a handbook modulus with crosshead strain, Fty comes out low (the lesson's table gives 294 MPa at a compliance that produces 53 GPa), so the 300 MPa may be a reduction error, not a weak lot. Ask for the raw record and re-reduce it with the fitted modulus before rejecting anything.

(b) KE = ½ × 90 × 4² = 720 J. Estimated toughness (Fty + Ftu)/2 × elongation and the strap mass m = KE ρ/U:

| alloy | U (MJ/m³) | U/ρ (kJ/kg) | minimum mass |
|---|---|---|---|
| 2024-T3 | 74.5 | 26.81 | **26.9 g** |
| 7075-T6 | 59.1 | 21.04 | **34.2 g** |
| Ti-6Al-4V | 128.1 | 28.92 | **24.9 g** |
| 4130 | 140.9 | 17.95 | **40.1 g** |

Ti-6Al-4V wins narrowly per kilogram, with 2024-T3 close behind; 7075-T6 loses because it is less ductile. These are lower bounds: a real absorber does not strain all its material uniformly to fracture. The modulus of resilience σ_y²/2E (293 J/kg for 2024-T3) is the energy an elastic spring stores and *returns*; a crash absorber must dissipate the energy, not give it back to the occupant, so toughness is the relevant property.

(c) A typical value is an average: about half of all material is weaker. An allowable is a statistically based lower bound (a value that a stated high fraction of production is expected to exceed with stated confidence), so that a structure sized to it is safe with the material actually delivered. More tests make the typical value more precise; they do not make it a lower bound.

(d) Stretch-forming a fuselage skin or a leading edge (or deep drawing, or crash simulation): strains of tens of per cent, where the material's response is the true stress against true strain, and the engineering curve, referred to the original area and shaped by the neck of one coupon, describes the coupon, not the material.

## Check-yourself answers

1. A₀ = 37.5 mm²; σ = 15 000/37.5 = 400 MPa; true stress 400 × 1.1 = 440 MPa.
2. True strain adds increments relative to the current length, so it is additive over successive stretches; the difference from ε is about ε²/2, under 1 % of ε below 2 % strain and significant in forming.
3. dσ_t/dε_t = σ_t. Beyond it the hardening no longer makes up for the shrinking area, so the force (and the engineering stress, force over the original area) falls while the neck thins.
4. G = E/2(1 + ν), K = E/3(1 − 2ν). K → ∞: the material cannot change its volume (incompressible), as plastic flow of metals nearly is.
5. 500 MPa: the plastic strain of (6) at σ = σ_y is exactly 0.002.
6. Not E. Fty with the fitted modulus and the elongation are both immune to a linear compliance.
7. MJ/m³ (= MPa × strain); σ_y²/2E = 345²/(2 × 73 100) = 0.814 MJ/m³.
8. Low: strain is the x variable of the fit, and noise in x flattens a least-squares slope by var(x)/(var(x) + σ²); noise in y only scatters the slope, it does not bias it.
