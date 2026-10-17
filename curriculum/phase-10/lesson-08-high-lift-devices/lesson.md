# P10.08  High-Lift Devices

> A flap adds camber; a slat shields a leading edge; a slot lets one element help another. Each buys lift for landing, and each is limited by where the flow on some element first gives up.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the thin-airfoil plain-flap results: the hinge angle, the flap effectiveness τ = 1 − (θ_h − sin θ_h)/π, the increments Δcl = 2πτδ, Δα_L0 = −τδ and Δcm_c/4 = −½ sin θ_h (1 − cos θ_h) δ, and the hinge moment by integrating the load over the flap. 2. Implement and test them, and check them by an independent route: the panel method with a geometrically deflected flap as the thickness goes to zero. 3. Build a two- and three-element panel model and show A. M. O. Smith's slat effect (a slat lowers the main element's leading-edge suction peak) and circulation effect (a slotted flap raises the main element's lift). 4. Estimate Δcl_max with a calibrated suction-peak criterion, state where it fails, and turn CL_max into stall and approach speeds. 5. Decide which device the reference trainer needs for its approach. |
| **Prerequisites** | `P10.04` (thin-airfoil theory, Glauert coefficients), `P10.05` (panel method, Kutta condition), `P10.07` (camber lines and the panel check; the multi-element panel solver), `P10.02` (NACA sections), `P01.13` (linear systems), `P02.09` (Fourier series). Skill check: you can evaluate ∫_{θ_h}^{π} cos nθ dθ by hand. |
| **Engineering question** | Can the reference trainer (750 kg, S = 12 m²) fly its 30 m/s approach at the customary margin of at least 1.3 times the stall speed, and which device makes that possible? Answer: clean, CL_max ≈ 1.44 gives V_s = 26.4 m/s and a margin of only 1.14; the 1.3 margin needs CL_max = 1.88. A 25 % plain flap at 20° gives V_s = 23.6 m/s (1.27, just short); a 30 % slotted flap at 30° gives V_s = 21.7 m/s (1.38), and stays above 1.3 for any assumed clean section cl_max from 1.4 to 1.8. |
| **Assumptions and limits** | Thin-airfoil flap theory: small deflections, a sharp hinge, inviscid flow, no gap at the hinge. Panel model: inviscid, incompressible, steady, one Kutta condition per element, no wakes, no boundary layers. Δcl_max criterion: the section stalls when the highest leading-edge suction peak of a forward element reaches the clean section's peak at an *assumed* clean cl_max of 1.6; it knows nothing of flap separation, confluent boundary layers or slat stall at the slat's own Reynolds number. 3-D: CL_max = 0.9 × (cl_max + 0.6 Δcl_max) (both factors are stated assumptions). The 1.3 V_s margin is customary practice for approach speeds, used here as the success criterion. |
| **Analytical work** | The kinked camber line's Glauert coefficients A0, A1, A2; τ, Δα_L0, Δcm and their limits (τ → (4/π)√E for small flaps, the largest |Δcm| at E = 0.25); the hinge moment integral; the stall speed V_s = √(2W/(ρ S CL_max)) and the CL_max needed for a given approach margin. |
| **Computational deliverable** | `starter/main.py`: `hinge_angle`, `flap_effectiveness`, `flap_increments`, `flap_slope`, `hinge_moment`, `flapped_section`, `panel_flap_derivative`, `element_peaks`, `critical_peak`, `clmax_by_peak`, `stall_speed`, `wing_clmax` (the NACA generator, `thin_airfoil`, the multi-element panel solver `panel_basis`/`panel_at`/`panel_solve`, `place_element`, `slat_configuration`, `slotted_flap_configuration` and `main` are given). Tests prove: the hinge sits where cos θ_h = 2E − 1; τ runs from 0 to 1 with the small-flap asymptote; the closed forms equal quadrature of the kinked camber line; flap increments are linear and superpose on camber; the largest nose-down Δcm per radian is 3√3/8 at E = 0.25; the hinge moment is restoring, linear, and has the whole-plate limit; the deflected-flap panel model approaches 2πτ as t → 0; far-apart elements do not interact; the slat and circulation effects; the cl_max criterion is self-consistent; stall speed scales correctly. |
| **Practical activity** | **computational experiment**: thin-airfoil flap theory against the panel method as the thickness goes to zero; a slat and a slotted flap in a multi-element panel model; a cl_max estimate for six configurations with a sensitivity study; the landing-speed payoff for the trainer. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `high_lift.png`, the derivation of τ and Δcm, the thickness-convergence table, the slat-effect table, the cl_max table with its assumptions, and a half-page note on the trainer's approach speed. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: τ and Δcm are derived; all public tests pass; the panel-vs-thin-airfoil comparison is reported with the trend in thickness; the note states the approach margin with and without flaps and the assumption that matters most. |
| **References** | `@ntrs-naca-tr-938-trailing-edge-high-lift` (Cahill, *Summary of Section Data on Trailing-Edge High-Lift Devices*, https://ntrs.nasa.gov/citations/19930092003): Flap Theory; Discussion of Experimental Results: Plain Flaps (Flap Loads and Moments; Summary of Plain-Flap Data); Single Slotted Flaps (required-free): measured flap effectiveness and Δcl_max against which this lesson's theory must be judged. `@leishman-intro-aerospace-flight-vehicles` "Airfoil Shapes" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-geometries/): Airfoils with Flaps; and "Stalling & Spinning" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/maximum-lift-stalling-spinning/): Importance of High-Lift Devices; Stall Control Devices (required-free). `@ntrs-naca-tr-824-summary-of-airfoil-data` "Experimental Characteristics": High-Lift Devices (https://ntrs.nasa.gov/citations/19930090976) (required-free). `@nasa-grc-beginners-guide-aeronautics` "Flaps and Slats Interactive" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/flaps-and-slats/) (required-free, qualitative). `@ntrs-nasa-cr-4746-high-lift-systems` (Rudolph, https://ntrs.nasa.gov/citations/19960052267) Chapter 1 "Types of High-Lift Systems: Their Geometry, Functions, and Design Criteria": 1.1 Types of High-Lift Systems; 1.3 Geometric Parameters of High-Lift Devices; 1.5 Characteristics and Constraints of High-Lift Devices; 1.6 Interactions between Leading- and Trailing-Edge Devices (supplementary-free): airliner high-lift systems. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Abbott & von Doenhoff, *Theory of Wing Sections*, for flap and slat section data; Anderson, *Fundamentals of Aerodynamics*, for thin-airfoil flap theory. |
| **Extensions** | ★ Add a gap and a rounded flap nose to the plain-flap geometry and see how the panel result changes. ★ Optimise the slat's position and deflection to maximise the criterion's cl_max, with the slat's own peak as a constraint. ★ Compare the thin-airfoil τ with the measured plain-flap effectiveness in TR-938's summary figures, and fit an empirical efficiency factor against deflection. |
| **Time** | unmeasured |

## 1. The engineering problem

A wing sized for cruise is too small for landing. The trainer cruises at 55 m/s at CL = 0.33; at its 30 m/s approach it needs CL = 1.11, and it must fly that approach with a margin above the stall, customarily at least 1.3 V_s, so that a gust or a mishandled flare does not stall it near the ground. Making the wing big enough for landing would make it draggy in cruise; instead the wing changes shape. A **trailing-edge flap** adds camber (more lift at a given angle). A **leading-edge slat** is a small airfoil ahead of the main element that protects the main element's nose (more angle before it stalls). A **slot** between elements lets each element help the others. Airliners carry all three, in several elements each; light aircraft usually carry only a plain or slotted flap.

This lesson quantifies those effects with three models: thin-airfoil flap theory (a flap is a kinked camber line), a multi-element panel method (the given Hess–Smith solver of `P10.07`, standing in for `P10.05`), and a calibrated criterion that turns inviscid pressures into a cl_max estimate. Each model is checked against the next one: the thin-airfoil flap against the panel flap as the thickness goes to zero; the panel model's two-element interactions against the explanations A. M. O. Smith gave for why slots work.

Success criteria, written first: the closed forms equal numerical quadrature of the kinked camber line to 10⁻⁴; the panel flap derivative at 2 % thickness is within 1 % of 2πτ; the slat lowers the main element's suction peak at fixed α; the cl_max criterion returns exactly the assumed clean cl_max on the clean section; the approach-margin answer is reported with its sensitivity to the assumed clean cl_max.

## 2. Understand and derive

### 2.1 A plain flap is a camber line with a kink

Hinge the aft fraction E = c_f/c of a flat plate at x_h = 1 − E and deflect it by δ (trailing edge down positive). The camber slope is zero ahead of the hinge and −δ behind it (small angles). In Glauert's variable x = (1 − cos θ)/2 the hinge sits at

**cos θ_h = 1 − 2x_h = 2E − 1.**

From `P10.04`, A0 = α − (1/π)∫₀^π (dz/dx) dθ and Aₙ = (2/π)∫₀^π (dz/dx) cos nθ dθ:

A0 = α + δ(π − θ_h)/π,  Aₙ = (2δ/π) sin(nθ_h)/n,  so A1 = 2δ sin θ_h/π, A2 = δ sin 2θ_h/π.

Then cl = 2π(A0 + A1/2) = 2πα + 2δ(π − θ_h + sin θ_h), and

**Δcl = 2πτδ,  τ = 1 − (θ_h − sin θ_h)/π,  Δα_L0 = −τδ.**   (1)

τ is the **flap effectiveness**: one degree of flap is worth τ degrees of angle of attack. Limits: E → 1 (the whole plate turns) gives θ_h → 0 and τ → 1; E → 0 gives θ_h → π − 2√E and τ → (4/π)√E. Small flaps are worth far more than their chord: a 10 % flap has τ = 0.396, a 25 % flap 0.609. That square root is why ailerons and elevators can be narrow.

The moment: cm_c/4 = (π/4)(A2 − A1) gives

**Δcm_c/4 = −½ sin θ_h (1 − cos θ_h) δ.**   (2)

A trailing-edge-down flap always pitches the section nose down. Differentiating, the largest |Δcm| per radian is at cos θ_h = −½, θ_h = 120°, which is **E = 0.25 exactly**, where it is 3√3/8 = 0.6495. For E = 0.25 and δ = 20°: Δcl = 1.336, Δα_L0 = −12.18°, Δcm_c/4 = −0.227.

### 2.2 The hinge moment

The pilot (or an actuator) holds the flap against the moment of its own load about the hinge. With H positive when it would push the trailing edge down, and the load Δcp = 4[A0 cot(θ/2) + Σ Aₙ sin nθ] of `P10.07`,

**c_h = H/(q c_f²) = −(1/E²) ∫_{x_h}^{1} Δcp (x − x_h) dx.**   (3)

Both ∂c_h/∂α and ∂c_h/∂δ are negative: an upward load behind the hinge pushes the trailing edge up. For E = 0.25 the code gives c_hα = −0.565 and c_hδ = −0.943 per radian. A check: as E → 1 the "flap" is the whole plate turning about its leading edge, and c_h must become the leading-edge moment of a flat plate, −cl/4 = −(π/2)(α + δ).

### 2.3 Why slots work: Smith's effects

A slat or slotted flap is a separate airfoil. In potential flow every element sees the velocity field of all the others, and A. M. O. Smith's classic explanation of multi-element lift names the consequences: (i) the **slat effect**: the slat's circulation induces a velocity on the main element's nose that opposes the local flow, lowering its leading-edge suction peak, so the main element can reach a higher angle before its boundary layer separates; (ii) the **circulation effect**: an element downstream (a flap) puts the trailing edge of the element ahead of it in a region of high velocity, which raises that element's circulation through its Kutta condition; (iii) the **dumping effect**: the upstream element's boundary layer leaves at that higher velocity, so it has less pressure to recover; (iv) off-surface pressure recovery and (v) a fresh boundary layer on each element, which are viscous effects outside this lesson's inviscid model. The multi-element panel method shows (i)–(iii) directly. Leishman's "Stalling & Spinning", Importance of High-Lift Devices, describes the devices (`@leishman-intro-aerospace-flight-vehicles`); Rudolph's CR-4746, 1.6 Interactions between Leading- and Trailing-Edge Devices, discusses them on airliners (`@ntrs-nasa-cr-4746-high-lift-systems`).

### 2.4 Estimating cl_max, honestly

No inviscid model predicts stall. What it can do is compare *pressure distributions*: if a section stalls when its leading-edge suction peak reaches some value (a common behaviour of thin and moderately thick sections whose stall starts at the leading edge), then a device that delays that peak delays the stall. The **suction-peak criterion** used here:

1. Assume the clean section's cl_max (here 1.6 for the NACA 2412: an assumption, not a measurement; replacing it with TR-824's measured value is the first improvement, and `P10.15` fits stall to data).
2. Find the angle at which the clean panel model reaches that cl (11.24°) and record its leading-edge peak: **cp_crit = −7.234**.
3. For each configuration, increase α until the leading-edge peak of a forward element (main or slat) reaches cp_crit; the total cl there is the estimate.

The criterion is self-consistent on the clean section by construction (a test). Its limits are as important as its result: it ignores separation on the flap itself (a plain flap at 20–30° stalls at its hinge, which the criterion does not see); it treats a slat's small nose at its lower Reynolds number like the main element's; it excludes the slotted flap, whose fresh boundary layer tolerates much higher peaks; and the 1.6 is an assumption, tested in part 4.

### 2.5 From cl_max to stall speed

At the stall in level flight, W = ½ρV_s² S CL_max, so

**V_s = √(2W/(ρ S CL_max)).**   (4)

For the trainer, W = 750 × 9.81 = 7357.5 N. Two stated assumptions turn section values into a wing CL_max: CL_max = 0.9 cl_max for an unswept untwisted wing (the 3-D wing reaches its peak local cl at a lower CL; `P10.12`), and flaps over 60 % of the wing area add 0.6 Δcl_max. The approach margin V_app/V_s ≥ 1.3 at V_app = 30 m/s needs CL_max ≥ 2W/(ρS(30/1.3)²) = **1.880**.

## 3. Build the computational model

Open `starter/main.py`. Given: the NACA generator, `thin_airfoil`, the multi-element panel solver (`panel_basis` solves the system once for a stream along x and once along y; `panel_at` superposes them at any α), `place_element`, `slat_configuration` (a NACA 6418-shaped slat of chord 0.15 deflected 30° nose-down, trailing edge at (0.03, 0.045) above the main element's nose), `slotted_flap_configuration` (a NACA 4412-shaped flap of chord 0.30 with its leading edge at (0.95, −0.03)), and `main`.

### Step 1: flap theory

`hinge_angle`, `flap_effectiveness` and `flap_increments` are (1) and (2). `flap_slope(E, δ, base)` returns the kinked slope, optionally added to a base camber line, so that the given `thin_airfoil` can check the closed forms by quadrature. `hinge_moment` evaluates (3) in θ from the Glauert coefficients (300 terms, 8000 points; the load is log-singular at the hinge, and the θ-midpoint rule handles it).

### Step 2: a flap with thickness

`flapped_section(code, E, δ)` turns every node aft of the hinge about the hinge point on the mean line. On the lower surface the rotation can push nodes forward of their neighbours and fold the contour; drop those nodes. `panel_flap_derivative(t, E)` differences cl at δ = ±2° for a NACA 00tt section.

### Step 3: elements and the criterion

`element_peaks(basis, α)` returns each element's cl (by pressure integration over that element) and its upper-surface leading-edge peak. `critical_peak` calibrates cp_crit with `brentq`; `clmax_by_peak` scans α in 0.5° steps and refines the first crossing with `brentq`, checking only the elements listed. `stall_speed` and `wing_clmax` are (4) and the 3-D rule.

Running `main()` prints (abridged):

```
thin-airfoil plain flap: E = cf/c, theta_h, tau, dcl/ddelta (per rad), dcm/ddelta, ch_alpha, ch_delta (per rad)
  0.10  143.13 deg  0.3958  2.4870  -0.5400  -0.3449  -0.8840
  0.25  120.00 deg  0.6090  3.8264  -0.6495  -0.5653  -0.9434
  0.50   90.00 deg  0.8183  5.1416  -0.5000  -0.8584  -1.0658
E = 0.25, delta = 20 deg: dcl = 1.3357 (quadrature 1.3357), dalpha_L0 = -12.180 deg, dcm_c4 = -0.2267 (quadrature -0.2267)
panel dcl/ddelta (alpha = 0, +-2 deg) against thin theory 2 pi tau = 3.8264:
  t/c = 0.02: 3.8312  (ratio 1.0013)
  t/c = 0.04: 3.9286  (ratio 1.0267)
  t/c = 0.08: 4.0858  (ratio 1.0678)
  t/c = 0.12: 4.2339  (ratio 1.1065)
clean NACA 2412: assumed cl_max = 1.6, reached (inviscid) at alpha = 11.24 deg; critical LE peak cp = -7.234
```

and writes `high_lift.png`: τ against flap chord with the panel point; the clean and slatted pressure distributions at 10°; and the three-element geometry.

### Verify

`python -m pytest curriculum/phase-10/lesson-08-high-lift-devices` runs 16 tests (listed in the metadata). Two deserve comment. The far-apart test places a second NACA 0012 a thousand chords away: each element's cl then matches the isolated value to 0.2 %. At fifty chords it did not: the other element's circulation induces an upwash Γ/(2πd) that changes cl by about 0.8 %, which is physics, not error, and a useful reminder of how far circulation reaches. And the small-flap limit τ(10⁻⁶) = (4/π) × 10⁻³ was found by a test written first with the wrong expectation (τ → 0 "fast"); the square root is the result.

### Compare

The independent route: the panel method knows nothing of Glauert's series. Its flap derivative at 2 % thickness is 3.8312 per radian against 2πτ = 3.8264 (+0.13 %); it rises with thickness to +10.7 % at 12 %, as thickness raises the lift slope (`P10.07`). Grid uncertainty matters here: between 80 and 160 nodes per surface the 2 % value moves from 3.809 to 3.846 (1 %), the 12 % value from 4.224 to 4.239 (0.4 %); the thin section is the harder case, because the kink is sharp compared with the thickness. The two models agree in the limit where both are valid, within that 1 %.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Flap chord.** τ, Δcl/δ, Δcm/δ, c_hα, c_hδ for E = 0.1–0.5 (table above, also 0.2: 0.5498, 3.4546, −0.6400, −0.4994, −0.9231; 0.3: 0.6607, 4.1516, −0.6416, −0.6274, −0.9655; 0.4: 0.7478, 4.6985, −0.5879, −0.7448, −1.0128). Predict where |Δcm| peaks. *Finding:* at E = 0.25. Doubling a 25 % flap to 50 % adds only 34 % to τ, but c_hδ grows by 13 % and H scales with c_f², so the hinge moment for the same q and δ grows 4.5 times: actuator size grows much faster than effectiveness.
2. **Thickness.** The panel flap derivative for t = 2, 4, 8, 12 %. *Finding:* ratios to 2πτ of 1.0013, 1.0267, 1.0678, 1.1065, almost linear in t.
3. **The slat effect.** Main element alone and with the slat, same α:

| α | clean cl | clean peak cp | slat + main cl | main cl | main peak cp | slat peak cp |
|---|---|---|---|---|---|---|
| 0° | 0.256 | −0.572 | 0.240 | 0.284 | −0.572 | +0.038 |
| 5° | 0.858 | −1.916 | 0.883 | 0.726 | −1.079 | −1.026 |
| 10° | 1.453 | −5.879 | 1.521 | 1.145 | −1.711 | −3.022 |
| 15° | 2.038 | −12.384 | 2.147 | 1.523 | −2.413 | −10.591 |

*Finding:* at 10° the slat cuts the main element's peak from −5.88 to −1.71 and unloads the main element (cl 1.145 instead of 1.453), while the pair lifts slightly more (1.521). The slat does not add much lift at a given angle; it moves the stall to a higher angle, and the slat itself becomes the element that runs out of margin (−10.6 at 15°).

4. **The circulation effect.** Slotted flap at 30°. *Finding:* at α = 0 the main element's own cl rises from 0.256 to 1.864 (the flap adds 0.793 more), and the pressure at the main trailing edge falls from +0.510 to −0.888: the main element now sheds its flow into a fast stream (dumping). The price is a main-element nose peak of −5.40 already at α = 0.
5. **cl_max and Δcl_max** by the criterion (cp_crit = −7.234):

| configuration | α_max | cl_max | Δcl_max | limiting element | thin-airfoil lift shift |
|---|---|---|---|---|---|
| clean | 11.24° | 1.600 | — | main | — |
| plain flap 25 %, 10° | 7.99° | 1.932 | +0.332 | main | 0.668 |
| plain flap 25 %, 20° | 4.80° | 2.258 | +0.658 | main | 1.336 |
| slat 15 %, 30° | 13.15° | 1.916 | +0.316 | slat | — |
| slotted flap 30 %, 20° | 3.78° | 2.567 | +0.967 | main | — |
| slotted flap 30 %, 30° | 1.58° | 2.872 | +1.272 | main | — |
| slat + slotted flap 30° | 3.86° | 3.264 | +1.664 | slat | — |

Predict the plain-flap Δcl_max as a fraction of its lift shift. *Finding:* about one half: the flap shifts the lift curve up but also brings the critical peak to a lower angle. The slat's gain comes the other way, through a higher stall angle. With the slotted flap at 30°, the flap's own nose peak (−16.5) would exceed the criterion many times over if it applied; it does not, because the flap's boundary layer starts fresh in the slot flow. That exclusion is a judgment the criterion cannot make for you.

6. **The trainer** (W = 7357.5 N, S = 12 m², needed CL_max 1.880):

| configuration | CL_max | V_s (m/s) | V_app/V_s |
|---|---|---|---|
| clean | 1.440 | 26.37 | 1.138 |
| plain flap 20° | 1.795 | 23.61 | 1.270 |
| slotted flap 30° | 2.127 | 21.69 | 1.383 |
| slat + slotted flap | 2.339 | 20.69 | 1.450 |

Sensitivity to the assumed clean cl_max (the criterion recalibrated each time): 1.4 gives cp_crit = −5.415, clean V_app/V_s = 1.064, slotted 1.321; 1.8 gives −9.420, 1.207 and 1.443.

Raw data format: CSV tables `E, tau, dcl_ddelta, dcm_ddelta, ch_alpha, ch_delta` (1); `t_c, panel_dcl_ddelta, ratio` (2); `alpha, clean_cl, clean_peak, total_cl, main_cl, main_peak, slat_peak` (3); `config, alpha_max, cl_max, dcl_max, limiting_element` (5); `config, CL_max, V_s, margin, assumed_clmax` (6); plus `high_lift.png`.

## 5. Validate and improve

| quantity | theory | code | difference | cause |
|---|---|---|---|---|
| Δcl, E = 0.25, δ = 20° | 2πτδ = 1.3357 | quadrature 1.3357 | < 10⁻⁴ | — |
| Δcm, E = 0.25, δ = 20° | −0.2267 | quadrature −0.2267 | < 10⁻⁴ | — |
| dcl/dδ, E = 0.25 | 3.8264 | panel, t = 2 %: 3.8312 | +0.13 % | thickness, discretisation |
| dcl/dδ, E = 0.25 | 3.8264 | panel, t = 12 %: 4.2339 | +10.7 % | thickness (not in thin theory) |
| clean cl_max | assumed 1.600 | criterion 1.600 | 0 | by construction |

The flap theory is verified twice; its *validity* is limited by viscosity, which thickens and eventually separates the boundary layer behind the hinge as the deflection grows. TR-938's plain-flap summaries (`@ntrs-naca-tr-938-trailing-edge-high-lift`) are the measured data to hold it against (Extension 3). The cl_max criterion is a model with a calibration constant and no viscous physics; treat its Δcl_max as an estimate whose error is not known from this lesson alone.

The engineering answer: **clean, the trainer's stall speed of 26.4 m/s leaves a 30 m/s approach only 14 % above the stall; it needs CL_max 1.88. A 30 % slotted flap at 30° brings V_s to 21.7 m/s (margin 1.38), and the conclusion survives clean cl_max assumptions from 1.4 to 1.8 (margin 1.32–1.44). A 20° plain flap gets close (1.27) but not there, and its estimate is the least trustworthy, because plain-flap stall starts at the hinge.** To improve: replace the assumed 1.6 by the TR-824 data for the 2412, calibrate the criterion against TR-938's measured slotted-flap Δcl_max, and add a boundary layer (`P10.06`) on the flap.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(3) and of the E = 0.25 extremum; `high_lift.png`; the six experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the trainer's approach margin with each device and the assumptions that set it.

### Adversarial self-checks

- **What if the flap deflection doubles?** Thin theory doubles Δcl; the criterion gives Δcl_max of 0.975 at 30° and 1.279 at 40°, still rising. Real plain flaps lose effectiveness long before that because the flow separates at the hinge; the inviscid model cannot see it, and its large-deflection numbers should not be used.
- **What if the clean cl_max assumption is 20 % lower?** The calibration constant changes with it; the slotted-flap margin falls from 1.38 to 1.32 at cl_max 1.4, still above 1.3. The *clean* margin falls to 1.06.
- **Which assumption fails first?** The 3-D factors 0.9 and 0.6: the margin scales with their square root, and both are guesses until `P10.12`'s lifting line (and a flapped-span calculation) replaces them.
- **What measurement could falsify the model?** A pressure-tapped two-element section showing the main element's peak *not* falling with the slat at fixed α would contradict the potential-flow mechanism; a measured Δcl_max far below the criterion's would show it is viscosity-limited.
- **Does the result survive a change of method or grid?** The flap derivative changes by 1 % at 2 % thickness and 0.4 % at 12 % between 80 and 160 nodes per surface; thin theory and panel agree as t → 0 within that.
- **Which decision changes if the uncertainty doubles?** If Δcl_max were uncertain by ±0.4, the 20° plain flap's margin would span 1.21 to 1.34, either side of 1.3, while the slotted flap's would stay at or above 1.31; the choice of the slotted flap survives, the plain flap would need a flight test of the stall speed.

## Check yourself

1. Where is the hinge of a 30 % flap in Glauert's variable θ?
2. What is the flap effectiveness τ, and what are its limits for a very small and a full-chord flap?
3. Why does a trailing-edge-down flap always pitch the section nose down, and for which flap chord is the effect per degree largest?
4. What sign do c_hα and c_hδ have, and why does that make a flap "want" to float trailing edge up?
5. State Smith's slat effect and circulation effect in one sentence each.
6. Why does a slat add little lift at a given angle of attack but raise cl_max?
7. What does the stall speed become if CL_max rises by 44 %?
8. Name two things the suction-peak cl_max criterion cannot account for.
