# P05.06  Stress & Strain

> Force over area, stretch over length, and one line between them. Everything else in a tensile test is deciding where that line stops being true, and reading the record honestly.

| | |
|---|---|
| **Learning outcomes** | 1. Convert a force–elongation record into engineering and true stress and strain, derive the constant-volume link between them, and say where it fails (necking). 2. Apply the isotropic Hooke's law in three dimensions (E, ν, G, K) to uniaxial, biaxial and hydrostatic stress states, including a pressurised fuselage skin. 3. Implement and test a stress–strain curve analyzer (modulus by least squares, 0.2 % offset yield, ultimate strength, elongation after fracture, toughness, resilience) and verify it against two independent results: the Ramberg–Osgood curve, whose offset yield is known by construction, and Considère's criterion, which says where the engineering curve peaks. 4. Quantify how the modulus window, measurement noise and machine compliance move the reported properties, and decide whether a reported property can be trusted. |
| **Prerequisites** | `P05.01` and `P05.02` (free-body diagrams and equilibrium: a cut through a coupon), `P02.14` (the stress tensor and index notation), `P01.17` and `P02.18` (least squares and regression), `P03.02` (root finding), `P03.04` (the trapezoid rule), `P02.16` (random variables), `P00.14` (reporting a number with its uncertainty). Skill check: you can fit a straight line by least squares and say what its residuals mean. |
| **Engineering question** | A 2024-T3 sheet coupon has been pulled to fracture (a **synthetic** record: 2000 points, load-cell noise 1 MPa, extensometer noise 2 × 10⁻⁵, seed 2024). What are its Young's modulus, 0.2 % offset yield, ultimate strength, elongation and toughness; how far can each number be trusted; and what goes wrong if strain is taken from crosshead travel instead of an extensometer? The answers: E = 71.4 GPa (72.5 ± 0.6 GPa over 200 seeds), Fty = 344 MPa (342.3 ± 1.2), Ftu = 479 MPa, elongation 17.4 %, toughness 81 MJ/m³; a machine compliance of 5 × 10⁻⁶ per MPa drops the apparent modulus to 53 GPa and, if the handbook modulus is then used for the offset line, the reported yield to 294 MPa. |
| **Assumptions and limits** | Uniaxial stress in the gauge section, uniform up to the maximum load; small strains for Hooke's law; isotropic, homogeneous material (rolled sheet is mildly anisotropic, which a single coupon cannot show); quasi-static loading at room temperature (no strain-rate or temperature effects); the synthetic material follows a Ramberg–Osgood true stress–strain law up to the maximum and a smooth imitation of necking after it, which is a model of a test, not a test. Material values are typical, not design allowables. |
| **Analytical work** | Engineering stress and strain from a free-body diagram and a gauge length; true stress and strain at constant volume; Hooke's law for a general stress tensor and its inverse with the Lamé constants; G = E/2(1 + ν) and K = E/3(1 − 2ν); strains in a pressurised cylinder; the Ramberg–Osgood curve and why its 0.2 % offset yield is σ_y; Considère's criterion dσ_t/dε_t = σ_t; toughness and resilience as areas; springs in series for machine compliance and why the offset construction survives it. |
| **Computational deliverable** | `starter/main.py`: `engineering_stress`, `engineering_strain`, `true_stress_strain`, `shear_modulus`, `bulk_modulus`, `hooke_strain`, `hooke_stress`, `ramberg_osgood_strain`, `considere_point`, `fit_modulus`, `offset_yield`, `ultimate_strength`, `elongation_after_fracture`, `toughness`, `resilience`, `remove_compliance`, `analyze_tensile` (`MATERIALS`, `moving_average`, `synthetic_tensile_test`, `TensileResult` and `main` are given). Tests prove: constant volume links the true and engineering measures; Hooke's law inverts, gives −ν contraction, G in pure shear and K under pressure, and no volume change at ν = 0.5; the analyzer returns σ_y of a Ramberg–Osgood curve to 10⁻⁶ and the modulus and toe of a line exactly; the maximum sits at Considère's point found independently by root finding; the toughness of an elastic–perfectly-plastic curve is exact and the trapezoid converges at second order; results scale with stress; a linear compliance changes E by springs in series and leaves the offset yield unchanged; seeded noise moves the results by statistically bounded amounts. |
| **Practical activity** | **computational experiment** on a synthetic tensile dataset (labelled synthetic wherever it appears): reduce one record, then study the modulus window, noise (200 seeds) and machine compliance, and compare the analyzer with Ramberg–Osgood and Considère. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `stress_strain.png`, the derivations of the true-stress relations, Considère's criterion and G = E/2(1 + ν), the three experiment tables, and a half-page data-reduction note for the coupon stating each property with its uncertainty and the strain source. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: the constant-volume relations and Considère's criterion are derived; all public tests pass; every reported property carries an uncertainty and a statement of how it was obtained (window, smoothing, strain source); the compliance experiment is explained, not just run. |
| **References** | `@openstax-university-physics-vol1` Ch. 12 "Static Equilibrium and Elasticity": 12.3 Stress, Strain, and Elastic Modulus (tensile, bulk and shear stress, strain and modulus); 12.4 Elasticity and Plasticity (required-free): the definitions and the shape of a stress–strain curve. `@mit-ocw-3-11-mechanics-of-materials` Roylance modules "Introduction to Elastic Response" (Tensile strength and tensile stress; Stiffness), "Stress-Strain Curves", and "Constitutive Relations" (Isotropic elastic materials; Hydrostatic and distortional components) (required-free): the engineering reading of the curve and Hooke's law in 3-D. `@mit-ocw-16-001-unified-materials-structures` lecture notes "Concept and Uses of Stress at a Point", "Strain; Transformation of Strain Components", "Constitutive Material Response", and Lab 01 §2 "Material Tensile Testing Lab" with Appendix B "Analyzing Material Tensile Strength Testing Data" (required-free): how a real test record is reduced. `@tudelft-alderliesten-aerospace-structures-materials` Ch. 1 "Material physics & properties": 1.2 Stress-strain, 1.8 Stiffness and apparent stiffness, 1.10 Toughness (required-free, link only: https://books.open.tudelft.nl/home/catalog/book/197): the aerospace framing. |
| **Extensions** | ★ Fit the Ramberg–Osgood exponent n to the plastic part of the record by least squares on log(ε_p) against log σ and report it with a confidence interval. ★ Add a strain-rate step to the synthetic generator (a sudden 2 % rise in flow stress) and teach the analyzer to flag it. ★ Find a public tensile dataset with a stated licence, reduce it with your analyzer, and compare with the published properties. |
| **Time** | unmeasured |

## 1. The engineering problem

Every structural number in this phase starts with a coupon. A strip of the sheet that will become a fuselage skin is cut to a dog-bone shape, clamped in a testing machine and pulled at a slow constant rate while a load cell records the force and an extensometer, two knife edges clipped 50 mm apart on the gauge section, records the stretch. The record is a few thousand pairs of numbers. From it a test laboratory reports five or six properties, and a stress engineer uses them for decades: the modulus E sets every deflection in `P05.12` and every buckling load in `P05.15`; the 0.2 % offset yield strength Fty is the limit-load allowable of `P05.14`; the ultimate strength Ftu is the ultimate-load allowable; the elongation says whether the material will warn before it breaks.

The reduction looks trivial, and it is where good engineers get quietly wrong numbers. Which points are "the straight part"? A modulus fitted over the wrong window can be 10 % low. Where was strain measured? A modulus from crosshead travel includes the stretch of the grips and the machine frame and can be 30 % low; if that low number is then mixed with a handbook modulus, the yield comes out wrong. How does noise move the maximum? The largest of 2000 noisy readings is biased high. A certificate that reports Fty = 294 MPa for 2024-T3 is not evidence of a bad batch until you know how it was reduced.

This lesson builds the analyzer and the physics behind it. Part 2 defines stress and strain from a free-body diagram, derives the true measures, writes Hooke's law in three dimensions and derives the two results the analyzer is checked against. Part 3 implements it. Part 4 runs it on a **synthetic** 2024-T3 record, labelled as such everywhere, because a synthetic record has known answers: the analyzer must give back what the generator was told. Success criteria, written first: the analyzer returns the σ_y of a Ramberg–Osgood curve to 10⁻⁶; the ultimate point agrees with Considère's criterion; the noise study shows how much each property moves; the compliance study shows which property survives a bad strain measurement.

The organising idea of the phase is *equilibrium, compatibility and constitutive law*. This lesson is about the third: the constitutive law, the material's own relation between stress and strain. Equilibrium enters once (the force on a cut through the coupon) and compatibility once (strain is the change in length over the length).

## 2. Understand and derive

### 2.1 Stress and strain from a free-body diagram

Cut the coupon across its gauge section and draw the free-body diagram of the lower half. The grip pulls down with force F; equilibrium needs an equal upward force on the cut face, spread over the original area A₀. The **engineering (nominal) stress** is

**σ = F / A₀,**   (1)

in N/mm² = MPa. If the gauge length L₀ grows by ΔL, the **engineering strain** is

**ε = ΔL / L₀**   (2)

(dimensionless). A 2024-T3 strap 25 mm × 3 mm carrying 20 kN has σ = 266.7 MPa; with the typical E = 73.1 GPa its strain is 3.648 × 10⁻³, so a 200 mm length stretches 0.730 mm, and the 25 mm width shrinks by ν ε b = 30.1 μm. The stress is 0.773 of the typical yield strength: a working stress, not a failure.

```
  stress
   ^             ultimate (Ftu) o
   |                   ___..--''   ''--.._        necking
   |           ..--''                     ''-.
   |  yield .-'                                x  fracture
   |      .'| 0.2 % offset line
   |     /  |/
   |    /   /
   |   /   /    slope E on both lines
   |  /   /
   | /   /
   |/___/__________________________________________> strain
   0  0.002                     eps_u           eps_f
```

### 2.2 True stress and strain

The engineering measures refer to the original geometry. As the coupon stretches its cross-section shrinks, so the stress on the material is F/A, not F/A₀. Metals deform plastically at nearly constant volume, A L = A₀ L₀, so A = A₀ L₀/L = A₀/(1 + ε) and the **true stress** is

**σ_t = F / A = σ (1 + ε).**   (3)

The **true (logarithmic) strain** adds up the increments dL/L: ε_t = ∫ dL/L = ln(L/L₀) = **ln(1 + ε)**. For small strains the two measures agree (at ε = 10⁻⁵ to one part in 10⁵); at 10 % strain the true stress is 10 % higher. Equations (3) hold only while the deformation is uniform: once a neck forms, the stress in the neck is not F over the average area. That is where the next result comes from.

### 2.3 Considère's criterion: where the curve peaks

The load F = σ_t A. As the coupon stretches, σ_t rises (the metal work-hardens) and A falls. The load peaks when the two balance, dF = 0: A dσ_t + σ_t dA = 0. At constant volume dA/A = −dL/L = −dε_t, so

**dσ_t / dε_t = σ_t   at the maximum load.**   (4)

This is **Considère's criterion**: necking starts when the hardening rate falls to the current true stress. Beyond it, any section that happens to be a little thinner carries a little more stress, hardens too slowly to compensate, and thins faster: the instability that ends the uniform part of the test. The engineering stress at the maximum is the **ultimate tensile strength** Ftu, and the engineering strain there, ε_u, the **uniform elongation**. Equation (4) gives the analyzer an independent check: if the material's true curve is known, the maximum of the record must be where (4) holds.

### 2.4 Hooke's law in three dimensions

In the elastic range a uniaxial stress σ produces the axial strain σ/E and a lateral strain −νσ/E in each perpendicular direction; ν is **Poisson's ratio**, 0.29–0.34 for the metals here. For a general stress tensor σ (`P02.14`), superposing the three normal stresses and treating shear separately gives the **isotropic Hooke's law**

**ε = [(1 + ν) σ − ν tr(σ) I] / E,**   (5)

whose inverse is σ = λ tr(ε) I + 2μ ε with the Lamé constants λ = Eν/((1 + ν)(1 − 2ν)) and μ = G. Two special states give the other two moduli.

*Pure shear.* With only σ₁₂ = τ, (5) gives ε₁₂ = (1 + ν)τ/E. The engineering shear strain is γ = 2ε₁₂, so τ = Gγ with **G = E / 2(1 + ν)**: 27.5 GPa for 2024-T3 and 42.4 GPa for Ti-6Al-4V. An isotropic material has only two independent elastic constants.

*Hydrostatic pressure.* With σ = −pI, tr(ε) = −3p(1 − 2ν)/E is the volume strain, and the **bulk modulus** is **K = E / 3(1 − 2ν)**. At ν = 0.5 the material is incompressible: (5) gives tr(ε) = 0 for any stress, which is the constant-volume assumption of 2.2 applied to elastic strain.

**A fuselage skin.** A cylinder of radius R = 1.975 m (a narrow-body fuselage) with skin t = 1.6 mm, pressurised to Δp = 55.2 kPa, carries a hoop stress pR/t = 68.1 MPa and an axial stress pR/2t = 34.1 MPa. Both are equilibrium results: cut the cylinder in half along its length and the pressure force 2pRL on the half is balanced by the wall force 2σtL; cut it across and pπR² is balanced by σ 2πRt. Equation (5) with typical 2024-T3 values gives a hoop strain 7.783 × 10⁻⁴, an axial strain 1.585 × 10⁻⁴ and a thickness strain −4.614 × 10⁻⁴. The radius grows ε_hoop R = 1.54 mm; ignoring the Poisson effect of the axial stress would give 1.84 mm, 20 % too much. Biaxial stress is the normal state of an aircraft skin, and a uniaxial modulus alone does not describe it.

### 2.5 The Ramberg–Osgood curve

Aluminium and titanium alloys have no sharp yield point: the curve bends over gradually. A two-parameter description that fits them well is the **Ramberg–Osgood** curve

**ε = σ/E + 0.002 (σ/σ_y)ⁿ,**   (6)

an elastic part plus a plastic part that grows as a power n of the stress. At σ = σ_y the plastic strain is exactly 0.002, so **the 0.2 % offset yield of (6) is σ_y by construction**. That is the analyzer's first independent check: hand it a sampled curve (6) and it must return σ_y. A larger n means a sharper knee; n of about 10 to 30 covers the alloys of this phase.

The lesson's synthetic 2024-T3 material uses (6) as its *true* stress–strain law with E = 73.1 GPa, σ_y = 345 MPa and n = 9, converted to engineering values by (3). Solving (4) for it (with (6), dε_t/dσ_t = 1/E + 0.002 n σ_t^(n−1)/σ_yⁿ, so (4) becomes σ_t/E + 0.002 n (σ_t/σ_y)ⁿ = 1, one equation for σ_t solved by bisection) gives **ε_u = 0.1249 and Ftu = 478.9 MPa**, close to the typical 483 MPa for this alloy. After the maximum the generator imitates necking with a smooth 10 % fall of engineering stress to fracture at ε_f = 0.18.

### 2.6 Reading the curve

The five reduction steps, each a decision:

1. **Modulus.** Fit a straight line σ = E(ε − ε₀) by least squares to the points between 10 % and 40 % of the maximum stress, before the maximum. The intercept ε₀ is the *toe*: slack in the grips, seating of the extensometer, a zero offset. Every later strain is measured from it. (A Ramberg–Osgood curve has no exactly straight part, so the fitted E depends slightly on the window; part 4 measures how much.)
2. **Offset yield.** Draw a line of slope E through (ε₀ + 0.002, 0) and find where the curve first crosses it; interpolate linearly between the two straddling points.
3. **Ultimate.** The largest engineering stress, after a short centred moving average so that load-cell noise does not bias the maximum upward (on the noisy coupon: 481.4 MPa unsmoothed, 479.3 MPa after 21 points).
4. **Elongation after fracture.** When the coupon breaks the elastic strain springs back; what is left, and what a technician measures by fitting the halves together, is ε_f − σ_f/E − ε₀.
5. **Toughness and resilience.** The area under the whole curve is the energy absorbed per unit volume up to fracture, MPa × strain = MJ/m³; the trapezoid rule is exact for the piecewise-linear record and second-order accurate for the underlying curve. The **modulus of resilience** σ_y²/2E is the elastic energy stored at yield, the energy a part can absorb and still spring back.

### 2.7 Machine compliance

If strain is computed from crosshead travel, the measured stretch includes the grips, the load train and the frame, which act as a spring in series with the coupon. If that extra stretch is proportional to the force, it adds c σ to the strain (c in 1/MPa) and the apparent modulus is

**1/E_app = 1/E + c,**   (7)

springs in series. Now look at the offset construction. The recorded strain is ε_rec = σ/E + cσ + ε_p = σ/E_app + ε_p, with ε_p the plastic strain. The offset line with slope E_app through (0.002, 0) is σ = E_app(ε_rec − 0.002); substituting ε_rec, it meets the record where E_app(ε_p − 0.002) = 0, that is where ε_p = 0.002, whatever c is. **The 0.2 % offset yield is immune to a linear compliance, provided the offset line uses the modulus fitted to the same record.** Mixing a handbook E with crosshead strain breaks this: the line is too steep, crosses too early, and reports a yield that is too low. The elongation after fracture is immune for the same reason.

### 2.8 Where the model stops

- Equations (3) and (4) assume uniform deformation at constant volume; after the maximum, the true stress in the neck needs a correction for the triaxial stress there, which this lesson does not make.
- Hooke's law (5) is for small strains of an isotropic material; rolled sheet has slightly different properties along and across the rolling direction, and composites are strongly anisotropic.
- The values in `MATERIALS` (2024-T3: E = 73.1 GPa, ν = 0.33, Fty = 345 MPa, Ftu = 483 MPa, elongation 18 %; 7075-T6, Ti-6Al-4V and 4130 in the printed table) are typical room-temperature values, rounded from alloy producers' published data sheets; you can cross-check the moduli against the "Material Properties" appendix of `@mit-ocw-3-11-mechanics-of-materials`. They are not design allowables: those are statistically based minima from a materials-properties handbook, specific to product form, thickness and grain direction. A coupon above the typical value says nothing about the next batch.
- The synthetic record has no strain-rate effects, no serrations, no extensometer slip and no grip failure; a real record can have all four.

## 3. Build the computational model

Open `starter/main.py`. `MATERIALS` (typical values: see 2.8), `moving_average`, the generator `synthetic_tensile_test`, the `TensileResult` dataclass and `main` are given.

### Step 1: the measures

`engineering_stress` and `engineering_strain` are (1) and (2), vectorised; `true_stress_strain` is (3) with `np.log1p`.

### Step 2: Hooke's law

`shear_modulus` and `bulk_modulus` are the formulas of 2.4 (return `np.inf` for K at ν = 0.5). `hooke_strain` is (5) on a 3 × 3 array; `hooke_stress` is its inverse with the Lamé constants. The tests check that one undoes the other on random symmetric tensors.

### Step 3: the material models

`ramberg_osgood_strain` is (6). `considere_point` solves σ_t/E + 0.002 n (σ_t/σ_y)ⁿ = 1 with `scipy.optimize.brentq` (bracket it between a tiny stress and 50 σ_y), then returns the engineering strain e^(ε_t) − 1 and stress σ_t e^(−ε_t) at that point.

### Step 4: the analyzer

```python
def fit_modulus(strain, stress, lo=0.10, hi=0.40):
    i_max = np.argmax(stress); smax = stress[i_max]
    sel = (np.arange(stress.size) <= i_max) & (stress >= lo * smax) & (stress <= hi * smax)
    slope, intercept = np.polyfit(strain[sel], stress[sel], 1)
    return slope, -intercept / slope          # E and the toe eps0
```

`offset_yield` forms f = σ − E(ε − ε₀ − 0.002) before the maximum, finds the first index where f goes from positive to non-positive, and interpolates. `ultimate_strength` takes the maximum of `moving_average(stress, window)`. `elongation_after_fracture`, `toughness` (`np.trapezoid`), `resilience` and `remove_compliance` (subtract cσ from the strain) are one line each. `analyze_tensile` chains them into a `TensileResult`.

Running `main()` prints

```
SYNTHETIC 2024-T3 coupon     E (GPa)  Fty (MPa)  Ftu (MPa)  eps_u   elong   toughness  resilience
  clean                        72.49      342.3      478.9  0.1248  0.1741    81.0     0.808
  noisy, seed 2024             71.43      344.4      479.3  0.1255  0.1740    81.0     0.830
```

after the material table, the strap and the fuselage-skin lines of 2.1 and 2.4, and the Ramberg–Osgood check (`offset yield of the exact curve 345.0000 MPa`), followed by the three experiment tables of part 4, and writes `stress_strain.png` (the synthetic record with its true-stress curve, offset line, yield and ultimate points; a zoom on the elastic region).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-06-stress-strain`. The tests:

1. **Measures**: (1) and (2) are linear; true and engineering measures satisfy σ_t = σ e^(ε_t) and agree as ε → 0.
2. **Hooke's law**: −ν contraction in uniaxial tension, no shear coupling; `hooke_stress` inverts `hooke_strain` to 10⁻⁹ MPa for three materials; τ/γ = G in pure shear, −p/tr(ε) = K under pressure, tr(ε) = 0 at ν = 0.5.
3. **Ramberg–Osgood**: offset yield of three sampled curves (n = 6, 9, 25) equal to σ_y to 10⁻⁶; a smaller offset gives σ_y(o/0.002)^(1/n).
4. **Modulus**: slope and toe of a line with a toe recovered to 10⁻⁹.
5. **Considère**: `considere_point` against an independent root of dσ_t/dε_t = σ_t to 10⁻⁸, and the maximum of the generated record at that point.
6. **Reduction of the clean coupon**: E within 1 %, Fty within 1 % (and below σ_y, because engineering stress is σ_t/(1 + ε)), elongation = ε_f − σ_f/E − ε₀.
7. **Toughness**: σ_y(ε_f − ε_y/2) for an elastic–perfectly-plastic curve, exactly; second-order convergence of the trapezoid on a Ramberg–Osgood curve.
8. **Scaling**: multiplying stress by 2.5 multiplies E, Fty, Ftu and toughness by 2.5 and leaves strains unchanged.
9. **Compliance**: apparent E from (7); offset yield and elongation unchanged; corrected strain restores E; the handbook modulus with crosshead strain under-reads the yield by more than 10 %.
10. **Noise**: 40 seeds keep the mean E within 0.4 GPa and the mean Fty within 1 MPa of the clean values; the same seed gives the same record.
11. **main** writes the figure.

### Compare

There is no library tensile analyzer in the scientific stack to compare with, which is the point of building one; the comparisons are with mathematics. `scipy.optimize.brentq` solves Considère's equation independently in the tests, and the analyzer's maximum lands on it to within one strain step (9 × 10⁻⁵). `np.polyfit` is the least-squares fit of `P01.17`; its slope is the same as the normal-equations solution you wrote there. For a real record, the standard test methods for metals specify the modulus window and the strain measurement; the Lab 01 appendix of `@mit-ocw-16-001-unified-materials-structures` shows a reduction of real data that you can repeat with this code.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment on a synthetic dataset.** Every record here is generated by `synthetic_tensile_test` and is labelled synthetic. Write each prediction before running.

1. **Modulus window.** Reduce the noisy coupon (seed 2024) with windows 5–20 %, 10–40 %, 20–60 % and 30–80 % of the maximum. Predict which window gives the highest E and which one breaks.
2. **Noise.** For three noise levels (0.5 MPa and 10⁻⁵; 1 MPa and 2 × 10⁻⁵; 3 MPa and 5 × 10⁻⁵), reduce 200 seeds and record the mean and standard deviation of E and Fty. Predict how the spread scales with the noise, and whether the mean moves.
3. **Compliance.** Add c = 0, 2, 5 and 10 × 10⁻⁶ per MPa. Record the apparent E, Fty with the fitted E, Fty with the handbook E = 73.1 GPa, the elongation, and E after `remove_compliance`. Predict which numbers move.
4. **Analyzer against the two closed forms.** Run `offset_yield` on a sampled Ramberg–Osgood curve and `ultimate_strength` on the clean record; compare with σ_y and with `considere_point`.

Raw data format: CSV tables `lo, hi, E_GPa, Fty_MPa` (1); `noise_MPa, noise_strain, mean_E, sd_E, mean_Fty, sd_Fty` (2); `c_per_MPa, E_app_GPa, Fty_fitted, Fty_handbook, elongation, E_corrected` (3); plus `stress_strain.png`, each file headed "SYNTHETIC".

What you should find. (1) E = 72.89, 71.43, 64.92 and 21.34 GPa, Fty = 342.2, 344.4, 347.1 and 387.6 MPa: the low windows agree within 2 %; 20–60 % already reaches into the knee (the Ramberg–Osgood plastic strain at 60 % of the maximum is no longer negligible) and 30–80 % fits a line through plastic data and is meaningless. (2) Spreads of 0.31, 0.61 and 1.60 GPa in E and 0.60, 1.17 and 3.25 MPa in Fty, roughly proportional to the noise; the mean Fty stays at 342.3 MPa, but the mean E falls to 71.67 GPa at the highest noise. That bias is *regression dilution*: noise in the strain (the x variable of the fit) flattens a least-squares slope by the factor var(x)/(var(x) + σ²_noise), here about 1 %. (3) Apparent E = 72.49, 63.31, 53.21 and 42.03 GPa, exactly (7); Fty with the fitted E stays at 342.3 MPa and the elongation at 0.1741 for every c, as 2.7 predicts; Fty with the handbook E falls to 341.5, 326.8, 294.0 and 196.8 MPa; `remove_compliance` restores 72.49 GPa. (4) 345.0000 MPa and 478.9 MPa at ε_u = 0.1249.

## 5. Validate and improve

| quantity | independent value | analyzer | difference | cause |
|---|---|---|---|---|
| offset yield of exact Ramberg–Osgood curve | σ_y = 345 MPa | 345.0000 MPa | < 10⁻⁶ | — |
| ultimate of the clean record | Considère: 478.9 MPa at ε_u = 0.1249 | 478.9 MPa at 0.1248 | one strain step | sampling every 9 × 10⁻⁵ |
| modulus, clean record | generator E = 73.1 GPa | 72.49 GPa | −0.8 % | Ramberg–Osgood curvature in the 10–40 % window, plus the engineering conversion |
| offset yield, clean record | true-curve σ_y = 345 MPa | 342.3 MPa | −0.8 % | engineering stress is σ_t/(1 + ε), 0.67 % lower at ε = 0.0067, and the fitted E |
| elongation | ε_f − σ_f/E − ε₀ | 0.1741 | 0 | — |
| apparent E with c = 5 × 10⁻⁶ /MPa | (1/72.49 GPa + c)⁻¹ = 53.2 GPa | 53.21 GPa | < 10⁻⁵ | — |
| Fty with crosshead strain, fitted E | unchanged | 342.3 MPa | 0 | the offset line moves with the data |

Every analyzer result agrees with its independent value or departs by an amount that 2.5 and 2.7 predict. The two departures in the clean record are not errors of the code: they are the difference between a material's true-stress law and the engineering-stress record a laboratory reports, and the fact that "the straight part" of a rounded curve depends on where you look. A report that gives E without the window, or Fty without the strain source, is incomplete.

**The engineering answer** for the synthetic coupon: E = 71.4 GPa from this record, with a single-test standard deviation of 0.6 GPa from the noise study, so 71.4 ± 1.2 GPa (2σ), lower than the 72.5 GPa clean value by 1.7σ, an unlucky but unremarkable draw; Fty = 344.4 ± 2.3 MPa (2σ); Ftu = 479 MPa; elongation 17.4 %; toughness 81 MJ/m³; resilience 0.83 MJ/m³. If the laboratory took strain from the crosshead, E is meaningless but Fty is still trustworthy as long as it was reduced with the fitted modulus.

What to improve: report E from several windows and quote the spread; correct for compliance by measuring the machine's own stiffness with a stiff dummy specimen; reduce strain noise rather than load noise, because strain noise biases E. What not to improve: more digits of Fty than the batch-to-batch scatter (several per cent for a real alloy) can support.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; derivations of (3), (4), G = E/2(1 + ν) and (7) with the offset-immunity argument; `stress_strain.png`; the four experiment tables, each marked synthetic; and a half-page data-reduction note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving each property with its uncertainty, the window, the smoothing and the strain source.

### Adversarial self-checks

- **What if the load doubles?** In the elastic range everything is linear: stress and strain double, E does not change. In the plastic range a doubled stress is not possible for this material (it fractures first); the test 2.5 × stress in the tests is a *different material* of the same shape, and every stress property scales while every strain property stays.
- **What if the noise goes to zero?** The analyzer must return the clean values (72.49 GPa, 342.3 MPa, 478.9 MPa) exactly; it does, because the window and smoothing are deterministic. With very high strain noise, E is biased low, not just scattered.
- **Which assumption fails first in the real world?** The uniform deformation after the maximum: the engineering curve beyond Ftu depends on the gauge length and the coupon geometry, and the true-stress formulas (3) are wrong there. Before that, the assumption that the modulus window is elastic fails for materials with a rounded knee.
- **What measurement could falsify the model?** A strain gauge or extensometer reading at a known load in the elastic range: if strain/stress differs from the fitted 1/E by more than the noise, the fit window or the strain source is wrong. A record whose maximum does not satisfy Considère's criterion on its true curve means the material is not deforming uniformly (a flaw, a grip failure, or a material that work-softens).
- **Does the result survive a change of method?** Changing the window from 5–20 % to 10–40 % moves E by 2 % and Fty by 0.6 %; that is the method uncertainty to quote. The unsmoothed maximum of the noisy record is 481.4 MPa, 2.5 MPa above the clean 478.9 MPa (the largest of many noisy readings is biased high); a 21-point average gives 479.3 MPa and a 51-point one 479.0 MPa.
- **Which decision changes if the uncertainty doubles?** Accepting a batch against a specified minimum Fty a little below the typical value: with ±2.3 MPa the coupon clearly passes; with an unknown strain source and a handbook modulus, the same record could be reported anywhere from 197 to 342 MPa (the compliance table), and a good batch could be wrongly rejected.

## Check yourself

1. A coupon of 12.5 mm × 3 mm carries 15 kN. What is the engineering stress, and what is the true stress if the strain is 10 %?
2. Why is the true strain ln(1 + ε) rather than ε, and when does the difference matter?
3. State Considère's criterion and explain in one sentence why the engineering stress falls after it is met.
4. Write G and K in terms of E and ν. What happens to K as ν → 0.5, and what does that mean physically?
5. What is the 0.2 % offset yield of a Ramberg–Osgood curve with σ_y = 500 MPa, and why?
6. Strain is taken from crosshead travel and the fitted modulus is 53 GPa for an aluminium alloy. Which of E, Fty (offset with the fitted E) and elongation can you still trust?
7. What are the units of toughness, and what is the modulus of resilience of a material with Fty = 345 MPa and E = 73.1 GPa?
8. Noise in the strain channel biases the fitted modulus. Which way, and why does noise in the load channel not do the same?
