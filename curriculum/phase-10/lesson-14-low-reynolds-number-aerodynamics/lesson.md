# P10.14  Low Reynolds Number Aerodynamics

> At a small drone's Reynolds number the laminar boundary layer separates before it can become turbulent; what decides the section's performance is whether the separated layer comes back, and how much pressure it still has to recover if it does not.

| | |
|---|---|
| **Learning outcomes** | 1. Compute and compare the chord Reynolds numbers of insects, birds, small UAVs, sailplanes and light aircraft, and explain which flow features change across that range. 2. Implement Thwaites' method on an arbitrary edge-velocity distribution, verify it against its Falkner–Skan closed forms, Blasius, the stagnation point and Howarth's retarded flow, and use λ = −0.09 to locate laminar separation; show that the separation point does not depend on Re. 3. Combine Thwaites with Michel's transition criterion and an empirical bubble-length correlation to decide whether a laminar separation bubble forms, how long it is, and below which Reynolds number it cannot reattach on the airfoil. 4. Explain, with computed pressure-recovery numbers, why a thin cambered plate copes better than a 12 % thick section near Re = 5 × 10⁴, and where it does not. 5. Read the Reynolds-number trends of cd and cl_max from authentic NACA 0012 data at the high-Re end and say what may and may not be extrapolated down to a drone. |
| **Prerequisites** | `P10.05` (panel method: inviscid surface speed), `P10.06` (Thwaites, Michel, viscous drag), `P10.13` (flat-plate friction laws), `P03.04` (numerical integration), `P00.14`. From Phase 08 (recapped in §2.2): `P08.11` (Thwaites' integral method), `P08.12` (transition). Skill check: you can explain what an adverse pressure gradient does to a boundary layer. |
| **Engineering question** | At the reference UAV's cruise (V = 12 m/s, chord 0.20 m, Re = 1.64 × 10⁵, section lift coefficient 0.445), where does the laminar boundary layer on the upper surface of a 12 % thick NACA 0012 separate, does natural transition get there first, how long is the separation bubble, and below which Reynolds number would the bubble no longer close? The answers: laminar separation at x/c = 0.29, before natural transition (which would need Re > 4.2 × 10⁵); a bubble about 0.24c long that closes; and closure fails below Re ≈ 5.5 × 10⁴ (3.3–7.7 × 10⁴ over the correlation's uncertainty). A thin cambered plate (NACA 4402) separates later (x/c = 0.76) with only a third of the pressure recovery left (ΔCp = 0.36 against 1.16), which is why it degrades more gently near Re = 5 × 10⁴. |
| **Assumptions and limits** | Inviscid edge velocity from a panel method (no displacement-thickness feedback); Thwaites' one-parameter laminar method (good to a few per cent in attached flow, unreliable close to separation); Michel's empirical transition criterion (calibrated on airfoils in low-turbulence flow; it knows nothing of free-stream turbulence, roughness or noise); a bubble length from one empirical constant, Ue_s ℓ/ν ≈ 5 × 10⁴, uncertain by a factor of about 1.5; no model of the bubble's own pressure distribution, of bursting, or of hysteresis. The friction-drag curve is a flat-plate model with a thickness factor; it ignores bubbles entirely. |
| **Analytical work** | Reynolds-number scaling; Thwaites' quadrature and its closed forms for Ue = s^m; the Re-independence of the laminar separation point; Michel's criterion; the bubble length and the minimum Reynolds number for closure Re_min = 5 × 10⁴ /((Ue_s/V)(s_TE − s_s)/c); the pressure recovery left after separation; laminar, mixed and turbulent friction drag against Re. |
| **Computational deliverable** | `starter/main.py`: `reynolds_number`, `thwaites`, `laminar_separation`, `michel_transition`, `bubble_length`, `bubble_analysis`, `natural_transition_re`, `cf_laminar`, `cf_turbulent`, `cf_mixed`, `section_cd_friction`, `clmax_vs_re`, `cd_at_zero_lift` (the data readers, `naca4`, `hess_smith`, `alpha_for_cl`, `upper_surface_speed`, `surface_for_cl` and `main` are given). Tests prove: Thwaites reproduces θ² = 0.45νs^{1−m}/(5m + 1) and λ = 0.45m/(5m + 1) for four power laws, Blasius within 1.2 %, the stagnation value λ = 0.075 and Howarth's separation within 3 %; separation is Re-independent and θ ∝ √ν; Michel's flat-plate transition Re_x is speed-independent; the bubble length scales as 1/Re and closure switches at Re_min; the thin cambered plate separates later with less than half the recovery left; the natural-transition Re separates the two regimes; the friction laws join continuously; the authentic data load completely and show the grit strip's penalty. |
| **Practical activity** | **computational experiment**: Thwaites against closed forms; laminar separation and bubbles on three sections across Re = 2 × 10⁴–10⁶; sensitivity to the bubble correlation; the friction-drag model against the authentic NACA 0012 data (Ladson, NASA TM-4074) at Re = 2 × 10⁶–1.9 × 10⁷; cl_max against Re from the same data. |
| **Equipment** | None. A computer with the environment from `P00.01`. (The bench rig of `P10.23` measures a section at Re ≈ 3 × 10⁴, where this lesson's physics dominates.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `low_re.png`, the derivations of §2, the experiment tables of §4 and a half-page note recommending a section type for the UAV's wing (and for a smaller, slower one), with the assumption the recommendation rests on. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: Thwaites is verified against closed forms; all public tests pass; the bubble analysis is reported with its correlation uncertainty; the note explains the role of pressure recovery after separation and states what the high-Re data can and cannot say about a drone. |
| **References** | `@mit-ocw-16-100-aerodynamics` Lecture 23 "Correlation Methods for Integral Boundary Layers" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre23/) (required-free: Thwaites' method). `@leishman-intro-aerospace-flight-vehicles` "Boundary Layer Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/): Laminar Separation Bubbles; "Gliders & Sailplanes" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/gliders-sailplanes/): Sailplane Aerodynamics; Airfoils for Sailplanes (required-free, link only). `@ntrs-nasa-tm-4062-eppler-387-low-re` (McGhee, Walker & Millard, NASA TM-4062, https://ntrs.nasa.gov/citations/19890001471): Introduction; Presentation of Results; Discussion of Results (required-free: measured low-Re section data with oil-flow separation and reattachment lines, Re 60 000–460 000). `@ntrs-nasa-tm-4074-ladson-naca0012` (Ladson, NASA TM-4074, https://ntrs.nasa.gov/citations/19880019495; transcribed in `curriculum/phase-10/data/`): Tests and Procedures; Presentation of Results (required-free: the data used in §4). `@uiuc-selig-summary-low-speed-airfoil-data` *Summary of Low-Speed Airfoil Data* Vol. 1 (https://m-selig.ae.illinois.edu/uiuc_lsat.html): 3 "Summary of Airfoil Data" (supplementary-free, link only: measured drone and sailplane sections). `@drela-xfoil` XFOIL User Primer: General Description; "Re, Mach dependence" (supplementary-free: a code that models transitional bubbles, for Extension 2). `@mit-ocw-16-13-aerodynamics-viscous-fluids` Lecture 23 (transition prediction) (supplementary-free). Optional textbooks for deeper study: M. Drela, *Flight Vehicle Aerodynamics* (transition and bubbles inside a coupled integral boundary-layer method). |
| **Extensions** | ★ Replace the bubble-length constant by a model that turns the laminar shear layer turbulent at the bubble's transition point and closes it with Head's method (`P10.06`), and compute the drag increment of the bubble. ★ Run XFOIL (optional, outside the tested code) on the NACA 0012 and a thin cambered section at Re = 5 × 10⁴–2 × 10⁵ and compare its predicted bubble positions with the Thwaites separation points here. ★ Compare the separation and reattachment lines of the Eppler 387 oil-flow results in NASA TM-4062 with a Thwaites analysis of that section (its coordinates are tabulated in the report). |
| **Time** | unmeasured |

## 1. The engineering problem

Scale matters in aerodynamics through the **Reynolds number** Re = Vc/ν, the ratio of inertial to viscous forces at the scale of the chord. The reference trainer cruises at Re = 4.5 × 10⁶; the reference UAV (chord 0.20 m, 12 m/s) at 1.64 × 10⁵, nearly thirty times lower; a bumblebee near 800. The published airfoil data most engineers use (`P10.06`, `P10.13`) were measured at millions. Below about 5 × 10⁵ a section does not simply get a little more drag: its boundary layer stays laminar long enough to meet the adverse pressure gradient of the aft surface while still laminar, and laminar layers separate easily. The separated layer may become turbulent in the air and reattach, enclosing a **laminar separation bubble**, or it may never reattach. The section's drag can double, its maximum lift fall, and its lift curve develop kinks and hysteresis. Small-drone, model-aircraft and sailplane designers choose airfoils with this in mind; a choice made from high-Re data can fail badly.

This lesson asks a concrete question about the reference UAV's wing at cruise (section cl = 0.445): where does the laminar layer separate, does a bubble form, does it close, and at which Reynolds number would it stop closing? It answers with a chain of simple models (panel method → Thwaites → Michel → a bubble-length correlation), compares a 12 % thick section with a thin cambered one, and reads the high-Re end of the trends from authentic NACA 0012 data.

Success criteria, written first: Thwaites reproduces its Falkner–Skan closed forms to 0.2 % and Howarth's separation point to 3 %; the separation point is the same at every Re to round-off; the bubble length scales exactly as 1/Re; the results are reported with the bubble correlation's uncertainty, and nothing below the data's Re range is presented as data.

## 2. Understand and derive

### 2.1 Reynolds numbers across flight

| flyer (illustrative round values) | V (m/s) | c (m) | Re = Vc/ν |
|---|---|---|---|
| fruit fly | 2.0 | 0.0008 | 110 |
| bumblebee | 4.0 | 0.003 | 821 |
| pigeon | 15 | 0.10 | 1.03 × 10⁵ |
| reference UAV | 12 | 0.20 | 1.64 × 10⁵ |
| model sailplane | 10 | 0.25 | 1.71 × 10⁵ |
| sailplane | 30 | 0.80 | 1.64 × 10⁶ |
| reference trainer, approach | 30 | 1.2 | 2.46 × 10⁶ |
| reference trainer, cruise | 55 | 1.2 | 4.52 × 10⁶ |

(Sea-level ν = 1.461 × 10⁻⁵ m²/s; the speeds and chords are round illustrative numbers, not measurements.) Insects live where viscosity dominates and flat, unsteady flapping wings work best. Birds, drones and model aircraft sit in the awkward band 10⁴–5 × 10⁵ where bubbles dominate. Sailplanes and light aircraft are above it, where the boundary layer transitions naturally before it can separate laminar.

### 2.2 Recap: Thwaites' method (`P08.11`) and transition (`P08.12`)

For a laminar boundary layer with edge velocity Ue(s) along the surface arc length s, the momentum-integral equation with Thwaites' empirical closure gives the momentum thickness θ in closed form,

**θ²(s) = (0.45ν/Ue⁶) ∫₀ˢ Ue⁵ ds,  λ = (θ²/ν) dUe/ds,**   (1)

and laminar separation where **λ = −0.09**. At a stagnation point (Ue = as) the limit is λ = 0.075 and θ is constant. For a power law Ue = s^m, (1) gives exactly θ² = 0.45νs^{1−m}/(5m + 1) and λ = 0.45m/(5m + 1); the flat plate (m = 0) gives θ = 0.671√(νs/U), within 1 % of the exact Blasius value 0.664. For Howarth's linearly retarded flow Ue = U(1 − s/L), (1) predicts separation at s/L = 0.123 against the similarity-solution value 0.120.

**Michel's criterion** (empirical, from airfoil transition data in low-turbulence flow): transition occurs where

Re_θ ≥ 1.174 (1 + 22 400/Re_s) Re_s^{0.46},  Re_θ = Ueθ/ν, Re_s = Ue s/ν.   (2)

### 2.3 The laminar separation point does not depend on Re

In (1), θ² is proportional to ν: halve ν (double Re) and θ falls by √2 everywhere. But λ = θ²Ue′/ν is independent of ν. So **the laminar separation point depends only on the shape of Ue(s), that is, on the section and α**. What Re decides is whether the layer becomes turbulent first: Re_θ grows as √Re at a fixed station and Re_s as Re, so at high enough Re Michel's criterion is met before λ reaches −0.09 and the turbulent layer, which tolerates adverse gradients far better, carries the flow to the trailing edge. At low Re the laminar layer separates first.

### 2.4 The separation bubble

After laminar separation the free shear layer is very unstable; it becomes turbulent after a distance ℓ, spreads, and may reattach. Bubble experiments show that ℓ, measured in units of ν/Ue_s, is roughly constant. This lesson uses

**ℓ = N_ℓ ν/Ue_s,  N_ℓ ≈ 5 × 10⁴**   (3)

as a stated empirical correlation (reported values span roughly 3–7 × 10⁴; treat it as uncertain by a factor of 1.5). In chord units ℓ/c = N_ℓ/(Re·Ue_s/V): the bubble shortens as 1/Re. If separation at arc length s_s plus ℓ lies beyond the trailing edge, the bubble cannot close on the airfoil. The **minimum Reynolds number for closure** is therefore

**Re_min = N_ℓ / ((Ue_s/V)(s_TE − s_s)/c).**   (4)

**Pressure recovery left after separation.** If the layer does not reattach, the pressure it would have recovered between separation and the trailing edge is lost, which is the source of the large drag and lift loss of an open separation. Its size is

**ΔCp_rec = (Ue_s² − Ue_TE²)/V².**   (5)

A section that separates late, with little recovery left, loses little; one that separates early, with a large recovery still to make, loses a lot.

### 2.5 Why thin cambered plates do well at Re ≈ 5 × 10⁴

At its design angle a thin cambered plate has no leading-edge suction peak: the flow arrives smoothly and the upper-surface speed rises gently to mid-chord and falls towards the trailing edge (Figure, left panel). A thick symmetric section at the same cl has a suction peak near the front and must recover more pressure. At cl = 0.445 the panel method gives peak speeds of 1.55V (NACA 0012) and 1.19V (NACA 4402, 4 % camber, 2 % thickness, the closest NACA 4-digit shape to a cambered plate). Thwaites then places separation at x/c = 0.29 on the 0012 with ΔCp_rec = 1.16 left, and at x/c = 0.76 on the 4402 with ΔCp_rec = 0.36 left. When the Reynolds number is too low for either bubble to close, the cambered plate has far less to lose. The flip side (Problem 4): at an angle above its design point the thin plate's sharp leading edge produces a large suction spike and separates almost at once.

### 2.6 Friction drag against Re

Ignoring bubbles, a section's zero-lift drag is roughly twice a flat plate's skin friction times a thickness form factor (`P10.13`):

cd ≈ 2C_f(Re)(1 + 2t/c + 60(t/c)⁴),   (6)

with C_f = 1.328/√Re laminar (Blasius), 0.455/(log₁₀Re)^{2.58} turbulent (Prandtl–Schlichting), or mixed with transition at Re_x = Re_tr. Laminar friction falls as Re^{−1/2}: at Re = 5 × 10⁴ the laminar estimate for a 12 % section is 0.0149, at 1.64 × 10⁵ 0.0082. A bubble adds drag on top of the laminar friction, so below 5 × 10⁵ (6) with a laminar C_f is a floor, not a prediction; measured low-Re data such as the Eppler 387 results of TM-4062 are the test.

### 2.7 Where the models stop

- Thwaites is least accurate just before separation, exactly where it is used; the separation point is good to a few per cent of chord.
- The inviscid Ue ignores the bubble's own displacement effect, which flattens the pressure plateau and moves separation.
- Michel's criterion assumes a quiet free stream; turbulence, noise, vibration and roughness all cause earlier transition (which is why turbulators work).
- Bubble bursting, the abrupt change from a short to a long bubble at higher α, and hysteresis are not modelled.

## 3. Build the computational model

Open `starter/main.py`. Given: `find_data_file`, `read_ladson`, `ladson_blocks` (the authentic data), `naca4` (from `P10.02`), `hess_smith` (a compact panel method; `P10.05` builds one from scratch), `alpha_for_cl`, `upper_surface_speed` (finds the stagnation point where the panel speed changes sign and returns s, x and Ue from there to the trailing edge), `surface_for_cl` (a cached version) and `main`.

### Step 1: Thwaites

`thwaites(s, ue, nu)` evaluates (1). The integral ∫Ue⁵ds must be accurate near the stagnation point, where Ue⁵ ~ s⁵ and the trapezoid rule is badly wrong on the first intervals (it was the first bug found by the stagnation-point test). Take Ue linear on each interval and integrate exactly: ∫Ue⁵ds = h(u₁⁶ − u₀⁶)/(6(u₁ − u₀)). Where Ue = 0 use the limit θ² = 0.075ν/Ue′. Then `laminar_separation(s, lam, start)` interpolates the first downward crossing of −0.09 after the speed peak.

### Step 2: transition and the bubble

`michel_transition(s, ue, theta, nu)` interpolates the first station that satisfies (2). `bubble_length` is (3). `bubble_analysis(code, cl, Re)` works in chord units (c = 1, V = 1, ν = 1/Re): it finds α for the cl, runs Thwaites and Michel, and returns the separation point, whether a bubble forms (separation before transition), its length, whether it closes, ΔCp_rec of (5) and Re_min of (4). `natural_transition_re` finds by Brent's method the Re at which Michel's transition reaches the separation point.

### Step 3: friction drag and the data

`cf_laminar`, `cf_turbulent`, `cf_mixed` (transition at Re_x = Re_tr), `section_cd_friction` (6), `clmax_vs_re` (largest measured cl per block) and `cd_at_zero_lift` (the row nearest α = 0 per block).

Running `main()` prints (abridged):

```
flat plate: theta sqrt(U/(nu x)) at x = 1: 0.6708 (Blasius 0.664)
Howarth linearly retarded flow Ue = 1 - x: Thwaites separation at x = 0.1231 (exact similarity-solution value about 0.120)
upper-surface laminar analysis at cl = 0.445, Re = 1.643e+05 (bubble Re_l = 5e+04):
  NACA 0012: alpha = 3.68 deg, peak Ue = 1.546 V, laminar separation x/c = 0.2875 (Ue_s = 1.268 V, theta_s = 1.013e-3 c), Michel transition s = 0.521, bubble yes, length 0.240 c, closes True, recovery left dCp = 1.162, Re_min = 5.51e+04, no-bubble above Re = 4.18e+05
  NACA 2412: alpha = 1.53 deg, peak Ue = 1.326 V, laminar separation x/c = 0.4040 (Ue_s = 1.249 V, theta_s = 1.047e-3 c), Michel transition s = 0.711, bubble yes, length 0.244 c, closes True, recovery left dCp = 1.114, Re_min = 6.65e+04, no-bubble above Re = 9.59e+05
  NACA 4402: alpha = -0.18 deg, peak Ue = 1.188 V, laminar separation x/c = 0.7564 (Ue_s = 1.111 V, theta_s = 1.477e-3 c), Michel transition s = 0.944, bubble yes, length 0.274 c, closes False, recovery left dCp = 0.363, Re_min = 1.83e+05, no-bubble above Re = 6.19e+05
```

and writes `low_re.png`: the upper-surface speeds of the 0012 and 4402 at cl = 0.445 with their separation points; the friction-drag curves with Ladson's data; and Ladson's cl_max against Re.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-14-low-reynolds-number-aerodynamics`. The tests: θ and λ against the power-law closed forms for m = 0, 0.1, 0.3, 1; θ ∝ √ν and Blasius within 1.2 %; λ = 0.075 and constant θ at a stagnation point; Howarth's separation within 3 % and identical for ν = 10⁻⁷–10⁻³; no separation in a favourable gradient; Michel's flat-plate transition at the same Re_x for three speeds and between 5 × 10⁵ and 5 × 10⁶; ℓ ∝ ν; x_s independent of Re, ℓ ∝ 1/Re, θ_s ∝ Re^{−1/2}; closure failing at 0.95 Re_min and succeeding at 1.05 Re_min; the cambered plate's later separation, gentler peak and less than half the recovery; bubbles below and none above the natural-transition Re; the friction laws, their continuity and the constant laminar-run correction; the data's 161 rows and 6 + 5 blocks, and the grit's higher cd and lower cl_max.

### Compare

Thwaites on the flat plate gives 0.6708 against Blasius' 0.664 (1.0 % high, the known offset of its closure); on Howarth's flow 0.1231 against 0.120 (2.6 % late). Those are the accuracies to carry into the airfoil results: separation points good to a few per cent of chord.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment**, with a comparison against an authentic public dataset (Ladson, NASA TM-4074; US Government work, public domain; NACA 0012 in the Langley Low-Turbulence Pressure Tunnel at M = 0.15, cl and cm from surface pressures, cd from a wake rake; transcribed with its provenance in `curriculum/phase-10/data/`). Write each prediction before running it.

1. **Three sections at the UAV's cruise.** cl = 0.445, Re = 1.64 × 10⁵. Predict which separates first. *Finding:* the table in Step 3's output. Laminar separation at x/c = 0.29 (0012), 0.40 (2412), 0.76 (4402); all three separate before natural transition; the bubbles are 0.24–0.27c long; the 0012's and 2412's close, the 4402's (starting so far aft) does not, but it leaves only ΔCp = 0.36 unrecovered against 1.16 and 1.11.

2. **Reynolds-number sweep.** 0012 against 4402 at cl = 0.445, Re = 2 × 10⁴–10⁶. Predict where each bubble opens. *Finding:* bubble lengths 1.97/2.25c at 2 × 10⁴ (both open), 0.79/0.90c at 5 × 10⁴ (both open), 0.39/0.45c at 10⁵ (0012 closes, 4402 open), 0.13/0.15c at 3 × 10⁵ (both close), and at 10⁶ natural transition comes first on both. Below 5.5 × 10⁴ the 0012 has an open separation from x/c = 0.29 with ΔCp = 1.16 to recover; the 4402 has an open separation from 0.76 with 0.36 to recover. That is the computed version of "thin cambered plates beat thick sections at Re ≈ 5 × 10⁴".

3. **Sensitivity to the bubble correlation.** N_ℓ = 3, 5, 7 × 10⁴. *Finding:* Re_min = 3.31, 5.51, 7.71 × 10⁴ (0012) and 1.10, 1.83, 2.57 × 10⁵ (4402): proportional to N_ℓ. The separation points do not depend on N_ℓ, so the qualitative conclusion of Experiment 2 survives; the Reynolds number at which it bites is uncertain by ±40 %.

4. **Friction drag against the authentic data.** NACA 0012 near α = 0 at Re = 2 × 10⁶–1.9 × 10⁷. Predict where the measurements sit between the laminar and turbulent curves. *Finding:* free transition 0.00620, 0.00659, 0.00650, 0.00630, 0.00634, 0.00655 at Re = 2.0, 3.94, 5.97, 8.86, 11.9, 18.9 × 10⁶, almost independent of Re; fully turbulent estimate 0.00987 → 0.00681; with transition at Re_x = 5 × 10⁵, 0.00785 → 0.00659. The measured free-transition drag is below both: in that low-turbulence tunnel the laminar run is longer than Re_x = 5 × 10⁵ (`P10.13` infers it). With the wraparound grit strip the measurements (0.01100 → 0.00842) lie 9–15 % above the fully turbulent estimate. At the UAV's Re the same model gives 0.0082 laminar and 0.0161 turbulent; no data in this file reach that Re.

5. **cl_max against Re from the data.** *Finding:* free transition 1.408, 1.484, 1.660, 1.693, 1.617, 1.659 at Re = 2.0, 3.94, 5.97, 8.86, 11.9, 18.9 × 10⁶; with the grit 1.165, 1.207, 1.217, 1.229, 1.121 at 2.0–12.1 × 10⁶. A straight-line fit in log Re to the free-transition values rises 0.273 per decade; extrapolated to 1.64 × 10⁵ it would give 1.15. That extrapolation is not valid: below about 5 × 10⁵ the mechanism changes (bubbles, bursting) and the data stop at 2 × 10⁶. The lesson's models do not predict cl_max at low Re; measured low-Re data such as TM-4062 or the UIUC volumes are needed.

Raw data format: CSV tables `section, alpha_deg, x_sep, ue_sep, s_tr, bubble, length_c, closes, dCp_rec, Re_min, Re_natural` (1); `Re, section, length_c, closes, bubble` (2); `N_l, Re_min_0012, Re_min_4402` (3); `table, Re, cd_measured, cd_lam, cd_tr5e5, cd_turb` (4); `table, Re, cl_max` (5).

## 5. Validate and improve

| quantity | reference | this lesson | difference | cause |
|---|---|---|---|---|
| flat-plate θ√(U/νx) | Blasius 0.664 | Thwaites 0.6708 | +1.0 % | Thwaites' closure |
| Howarth separation s/L | similarity solution ≈ 0.120 | 0.1231 | +2.6 % | Thwaites near separation |
| separation point vs Re | Re-independent (§2.3) | identical to round-off | 0 | — |
| 0012 cd, free transition, Re 2 × 10⁶–1.9 × 10⁷ | Ladson 0.0062–0.0066 | turbulent model 0.0099–0.0068 | +59 % to +4 % | laminar runs in the measurement |
| 0012 cd, wraparound grit | Ladson 0.0110–0.0084 | turbulent model 0.0099–0.0073 | −8 % to −13 % | grit roughness drag |

The engineering answer: **at the UAV's cruise the 0012's upper-surface laminar layer separates at x/c ≈ 0.29, natural transition would need Re > 4.2 × 10⁵, the bubble is about 0.24c long and closes, and closure is lost below Re ≈ 5.5 × 10⁴ (3.3–7.7 × 10⁴).** The UAV at 12 m/s has margin; the same wing at half the chord, or a smaller drone at 6 m/s, would not. The recommendation follows from §2.5: for a small, slow drone, a thin cambered section operated near its design cl; for the reference UAV, a moderately thick cambered section (or a turbulator ahead of the separation point) is acceptable. What would change it: measured low-Re data for the actual section (TM-4062 is the free model of such a measurement), and a bubble model that predicts its drag (Extension 1).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of §2.2–2.4; `low_re.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending section types for the reference UAV and for a smaller drone, with the assumption each rests on.

### Adversarial self-checks

- **What if the speed halves?** Re halves to 8.2 × 10⁴: the separation point does not move, the bubble doubles to about 0.48c and still closes on the 0012 (Re_min = 5.5 × 10⁴), but the margin is gone.
- **What if Re → ∞?** Transition precedes separation; bubbles vanish (Experiment 2 at 10⁶) and the turbulent friction law takes over.
- **Which assumption fails first in the real world?** Michel's quiet-flow transition: free-stream turbulence, propeller wash and surface roughness trigger transition earlier, usually helping at low Re.
- **What measurement could falsify the model?** Oil-flow or surface-pressure measurements of the separation and reattachment lines across Re, as in TM-4062: separation should not move with Re, reattachment should move forward as 1/Re.
- **Does the result survive a change of method or grid?** The panel resolution and the Thwaites quadrature change x_s in the fourth figure; the bubble correlation changes Re_min by ±40 %.
- **Which decision changes if the uncertainty doubles?** If N_ℓ were uncertain by a factor of 2, Re_min for the 0012 could reach 1.1 × 10⁵, close to the UAV's slow-flight Re: the recommendation would need a turbulator or a test.

## Check yourself

1. What is the chord Reynolds number of the reference UAV at cruise, and why is it in an awkward range?
2. State Thwaites' formula for θ² and the separation criterion.
3. Why does the laminar separation point not depend on the Reynolds number in Thwaites' method?
4. What is a laminar separation bubble, and what decides whether it closes on the airfoil?
5. How does the bubble length scale with Re, and what is Re_min?
6. Why does a thin cambered plate near its design angle lose less when its separation does not reattach?
7. Why is cd = 2C_f·FF with a laminar C_f a floor rather than a prediction below Re ≈ 5 × 10⁵?
8. Why may cl_max not be extrapolated from Ladson's 2 × 10⁶–1.9 × 10⁷ data down to the UAV's Re?
