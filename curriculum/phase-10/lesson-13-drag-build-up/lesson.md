# P10.13  Drag Build-Up

> An aeroplane's zero-lift drag is the sum of its parts' skin friction, scaled by how much faster the air must run over their thickness, plus the bits that are not streamlined at all; the hard part is knowing which correlation to trust.

| | |
|---|---|
| **Learning outcomes** | 1. Compute flat-plate skin friction for laminar, turbulent and transitional boundary layers and explain how the transition location changes it. 2. Build up an aeroplane's zero-lift drag coefficient from wetted areas, skin friction, form factors, interference factors and stated add-ons, and say which component dominates. 3. Combine C_D0 with a lifting-line span efficiency and an Oswald estimate into a drag polar, and find the best lift-to-drag ratio, its speed and the cruise drag power. 4. Test the section estimate cd = 2C_f(1 + 2t/c + 60(t/c)⁴) against authentic NACA 0012 drag data with free and fixed transition, quantify the discrepancy and explain it physically. 5. Decide which inputs of a drag estimate deserve the most care. |
| **Prerequisites** | `P10.09` (induced drag, the polar, best L/D), `P10.10` (lifting line, span efficiency), `P10.06` (boundary layers on airfoils, transition), `P10.16` (fuselage geometry; may be taken in either order), `P03.02` (root finding), `P00.14`. From Phase 08 (recapped in §2.1): `P08.14` (turbulent boundary layers and skin friction). Skill check: you can compute the surface area of a solid of revolution. |
| **Engineering question** | What is the reference trainer's zero-lift drag coefficient, its best lift-to-drag ratio and the speed for it, and how much power does its drag absorb in the 55 m/s cruise? And how well does the skin-friction-times-form-factor estimate it rests on reproduce measured section drag? The answers: C_D0 = 0.0318 (fixed landing gear 41 % of it), (L/D)max = 12.1 at C_L = 0.77 and 36.1 m/s, and 46.1 kW of drag power at 55 m/s; the fully turbulent section estimate is 4–59 % above the NACA 0012's free-transition drag (the measurement has long laminar runs) and 8–13 % below its wraparound-grit drag (the grit adds roughness drag). |
| **Assumptions and limits** | Incompressible, attached flow; each component's friction from a flat plate at its own Reynolds number; empirical form factors and interference factors (stated, not derived); wetted areas from simplified geometry; landing gear, engine cooling and miscellaneous drag as stated drag areas; a parabolic polar with an Oswald factor that includes an empirical lift-dependent profile-drag term; no trim drag, no propeller slipstream, sea level. The estimate is good to perhaps ±15 % on the airframe and worse on the add-ons. |
| **Analytical work** | Blasius and Prandtl–Schlichting friction laws; the mixed laminar–turbulent plate; form factors as velocity-increase corrections; wetted areas of a section and a body of revolution; C_D0 = ΣC_f FF Q S_wet/S + add-ons; the parabolic polar and its optimum C_L* = √(πeAR C_D0), (L/D)max = ½√(πeAR/C_D0); inverting the section estimate for an effective transition location. |
| **Computational deliverable** | `starter/main.py`: `cf_laminar`, `cf_turbulent`, `cf_mixed`, `form_factor_section`, `form_factor_body`, `section_cd`, `perimeter_ratio`, `body_wetted_area`, `component`, `trainer_components`, `cd0_build_up`, `oswald_e`, `polar`, `best_ld`, `cruise_drag`, `ladson_cd_near_zero_lift`, `effective_transition` (`find_data_file`, `read_ladson`, `naca4`, `trainer_fuselage_radius`, `lifting_line_e` and `main` are given). Tests prove: the friction laws and their mixed form reduce to the pure limits; later transition lowers friction; the laminar-run correction depends only on the transition Reynolds number; form factors tend to 1; perimeters and wetted areas match closed forms; drag scales with wetted area; C_D0 adds up; lifting line gives e < 1 for a rectangle and > 0.98 at taper 0.4; the polar's best L/D matches a brute-force search and has C_Di = C_D0 there; lift equals weight; the fully turbulent estimate lies above all free-transition data and below all grit data, the gap closing with Re; the effective transition moves forward with Re. |
| **Practical activity** | **data-based investigation**: the 2-D estimate cd = 2C_f FF against the authentic NACA 0012 section drag of Ladson (NASA TM-4074; free transition at Re = 2.0–18.9 × 10⁶ and a wraparound grit strip at 2.0–12.1 × 10⁶), with the discrepancy quantified and explained; plus the trainer's build-up and its sensitivity to the transition location and the gear drag. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `drag_build_up.png`, the component table, the polar numbers, the data comparison with its two discrepancy tables, and a half-page note stating C_D0 with an uncertainty and the input that dominates it. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the friction laws and the build-up are derived and stated with their assumptions; all public tests pass; the comparison with Ladson's data reports the discrepancy at every Re and explains its sign; the note ranks the inputs. |
| **References** | `@ntrs-nasa-tm-4074-ladson-naca0012` (Ladson, NASA TM-4074, https://ntrs.nasa.gov/citations/19880019495; transcribed in `curriculum/phase-10/data/`): Introduction; Wake Survey Rake; Tests and Procedures; Presentation of Results; Conclusions (required-free: the data of §4 and how they were measured). `@mit-ocw-16-unified-fluids` 16.01 Fluids (spring) Lecture 10 "Aircraft Performance Analysis" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f10_sp/): Parasite Drag Estimation (required-free). `@openvsp-docs` "Parasite Drag Tool" (https://openvsp.org/wiki/doku.php?id=parasitedrag): Coefficient of Friction Equations (Laminar; Turbulent); Form Factor Equations (supplementary-free: a catalogue of the correlations a design tool offers). `@ntrs-naca-tr-824-summary-of-airfoil-data` "Drag Characteristics of Smooth Airfoils"; "Effect of Surface Irregularities on Drag" (Drag with fixed transition) (supplementary-free). `@mit-ocw-16-100-aerodynamics` study material "CD NACA" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/8cdnaca/) (supplementary-free: C_D against Re for NACA 0004 and 0012). `@leishman-intro-aerospace-flight-vehicles` "Aerodynamics of Finite Wings": Drag Polar for an Airplane; Generalization of the Airplane Drag Polar (supplementary-free, link only). `@nasa-grc-beginners-guide-aeronautics` "Shape Effects on Drag" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/shape-effects-on-drag/) (supplementary-free). Optional textbooks for deeper study: S. F. Hoerner, *Fluid-Dynamic Drag* (the component drag data behind the section and body form factors used here); D. P. Raymer, *Aircraft Design: A Conceptual Approach* (the component build-up method with interference factors and add-ons as used in conceptual design). |
| **Extensions** | ★ Replace the fixed wing transition location by Michel's criterion on a panel-method pressure distribution (`P10.06`) at cruise C_L and see how C_D0 changes. ★ Rebuild the estimate with the OpenVSP Parasite Drag Tool's alternative friction and form-factor equations (free; optional, outside the tested code) and compare the spread with the uncertainty stated here. ★ Add a retractable-gear variant and a cowl-flap (cooling) variant and compute the change in cruise speed for the same power. |
| **Time** | unmeasured |

## 1. The engineering problem

Every performance number of an aeroplane (cruise speed for a given power, range, climb rate, glide angle) depends on its drag polar, and the polar's first number is the **zero-lift drag coefficient** C_D0. Long before anything is built, C_D0 must be estimated, and the standard way is a **component build-up**: treat each part (wing, fuselage, tails) as a streamlined body whose drag is mostly skin friction, correct that friction for the part's thickness with a **form factor**, add an **interference** penalty where parts meet, and add the bluff bits (landing gear, cooling inlets, antennas) as drag areas from experience. The method is quick and transparent, but most of its numbers are empirical correlations whose authors fitted them to data from other aircraft. A designer has to know which numbers matter and how far to trust them.

This lesson builds C_D0 for the reference trainer (750 kg, wing 12.0 m², span 10.0 m, rectangular NACA 2412 wing of chord 1.2 m, fuselage 7.0 m × 1.1 m, horizontal tail 2.0 m², cruise 55 m/s at sea level), turns it into a polar and a best L/D, and then tests its key ingredient, the 2-D section estimate cd = 2C_f(1 + 2t/c + 60(t/c)⁴), against authentic wind-tunnel data for the NACA 0012.

Success criteria, written first: every component's numbers are traceable to a stated formula or assumption; the polar's best L/D matches a brute-force search to 10⁻⁸; the comparison with Ladson's data reports the discrepancy at every Reynolds number, with a physical explanation of its sign and trend rather than "experimental error".

## 2. Understand and derive

### 2.1 Recap: skin friction on a flat plate (`P08.14`)

A smooth flat plate of length L in a stream V has a mean skin-friction coefficient C_f = (friction force)/(q × wetted area) that depends on Re_L = VL/ν and the state of the boundary layer:

- laminar (Blasius' exact similarity solution): **C_f = 1.328/√Re**;
- turbulent (Prandtl–Schlichting's fit to measurements): **C_f = 0.455/(log₁₀Re)^{2.58}**.

At Re = 4.52 × 10⁶ (the trainer's wing) these are 0.00062 and 0.00342: a turbulent layer has five times the friction. A real surface is laminar from the leading edge to a transition point x_tr and turbulent after it. Counting the turbulent layer as if it had started at the leading edge, removing its share over the laminar part and putting the laminar friction back,

**C_f = C_f,t(Re) − (x_tr/L)[C_f,t(Re_tr) − C_f,l(Re_tr)],  Re_tr = (x_tr/L) Re.**   (1)

The correction term Re_tr[C_f,t(Re_tr) − C_f,l(Re_tr)] depends only on Re_tr (for Re_tr = 5 × 10⁵ it is 1614 with these two laws; textbooks quote values near 1700 for slightly different turbulent laws). Equation (1) is itself approximate: the real turbulent layer starts with the laminar layer's momentum thickness, not from zero.

### 2.2 Form factors, interference and wetted area

A thick body speeds the air up over its surface and adds pressure drag, so its drag exceeds that of a plate with the same wetted area by a **form factor** FF. The correlations used here are empirical fits to drag data (their origins are in the optional textbooks of the References row; the OpenVSP Parasite Drag Tool lists several alternatives):

- sections (wings, tails): **FF = 1 + 2(t/c) + 60(t/c)⁴**, 1.252 for t/c = 0.12. The 2(t/c) term is roughly twice the mean velocity increase over a symmetric section (friction scales as velocity squared); the (t/c)⁴ term is pressure drag that grows quickly for thick sections;
- bodies: **FF = 1 + 60/f³ + f/400** with f = L/D the fineness ratio (`P10.16`), 1.249 for the trainer's f = 6.36.

The **interference factor** Q accounts for the extra drag where components meet (junction flows, mutual velocity increases). Values are judgement: here Q = 1.0 for the wing (a high or well-faired mid wing) and the fuselage, 1.05 for the tails.

The **wetted area** is computed from the geometry: for a wing, twice the exposed planform area times the section's perimeter divided by twice its chord (1.0208 for NACA 2412, 1.0146 for NACA 0010); for the fuselage, the surface of revolution ∫2πr√(1 + r′²)dx of an equivalent round body (the profile of `P10.16`).

### 2.3 The build-up

For each component with reference length ℓ (chord or body length),

**ΔC_D = C_f(Vℓ/ν, x_tr) · FF · Q · S_wet/S,**   (2)

and

**C_D0 = [ΣΔC_D + (D/q)_gear/S + (D/q)_cooling/S](1 + misc).**   (3)

The add-ons are stated assumptions for a light trainer: a fixed tricycle landing gear with wheel fairings and struts, (D/q)_gear = 0.15 m²; engine cooling and induction, 0.04 m²; and a miscellaneous 5 % for leakage, antennas, steps and gaps. Neither the gear nor the cooling drag is predicted by skin friction; they are bluff-body drag areas, and they are the least certain inputs.

### 2.4 Induced drag and the polar

Lifting-line theory (`P10.10`) gives the span efficiency of the untwisted rectangular wing: e_LL = 0.934 at AR = 8.33 (δ = 0.071). The fuselage interrupts the span loading; an empirical factor (1 − 2(d/b)²) = 0.976 for d/b = 0.11 reduces it to e_inv = 0.911. Profile drag also grows with C_L (the boundary layer thickens as the pressure gradients steepen); an empirical way to include it in a parabolic polar is

**1/e = 1/e_inv + K π AR C_D0,  K ≈ 0.38,**   (4)

which gives the aeroplane's **Oswald efficiency factor e = 0.707**. The polar is C_D = C_D0 + C_L²/(πeAR), and its best lift-to-drag ratio is where induced and zero-lift drag are equal:

**C_L* = √(πeAR C_D0),  (L/D)max = ½√(πeAR/C_D0),  V* = √(2W/(ρSC_L*)).**   (5)

### 2.5 The section estimate and what it leaves out

For a 2-D section at zero lift, (2) becomes **cd = 2C_f(Re, x_tr)(1 + 2t/c + 60(t/c)⁴)**, the factor 2 for the two surfaces. It ignores: the real velocity distribution (the friction is not that of a plate even with the form factor); the effect of pressure gradient on transition; surface roughness and, for a grit strip, the drag of the grit itself; wall corrections and measurement uncertainty in the experiment. §4 tests it.

## 3. Build the computational model

Open `starter/main.py`. Given: `find_data_file` and `read_ladson` (the authentic data), `naca4` (from `P10.02`), `trainer_fuselage_radius` (from `P10.16`), `lifting_line_e` (a compact Fourier lifting line from `P10.10`) and `main`.

### Step 1: friction

`cf_laminar`, `cf_turbulent` and `cf_mixed(Re, x_tr)`, (1), with the pure limits at x_tr = 0 and 1.

### Step 2: geometry and components

`form_factor_section`, `form_factor_body`, `section_cd`; `perimeter_ratio(code)` sums the panel lengths of the NACA loop; `body_wetted_area(x, r)` sums frustums. `component(name, L_ref, S_wet, FF, Q, V, x_tr)` is (2). `trainer_components(V, x_tr_wing)` assembles the wing (exposed area 12.0 − 1.1 × 1.2 = 10.68 m², laminar to 10 % chord by assumption), the fuselage and the two tails (assumed vertical tail 1.2 m², tail chord 0.8 m, NACA 0010, turbulent).

### Step 3: the aeroplane

`cd0_build_up` is (3); `oswald_e(CD0)` is (4) with the lifting-line value; `polar`, `best_ld` (5) and `cruise_drag(V)` (C_L from L = W, drag and drag power).

### Step 4: the data

`ladson_cd_near_zero_lift(rows, table)` picks, in each Reynolds-number block, the row with α closest to zero (all within 0.1°). `effective_transition(Re, cd)` solves section_cd(Re, 0.12, x_tr) = cd for x_tr with `brentq` (NaN if no x_tr in [0, 1) fits).

Running `main()` prints (abridged):

```
component build-up at 55 m/s:
  component         Re        Cf       FF      Q     S_wet m^2   D/q m^2   Delta C_D
  wing              4.52e+06  0.00310  1.2524  1.00    21.805    0.0846  0.00705
  fuselage          2.64e+07  0.00258  1.2487  1.00    18.214    0.0588  0.00490
  horizontal tail   3.01e+06  0.00367  1.2060  1.05     4.058    0.0188  0.00157
  vertical tail     3.01e+06  0.00367  1.2060  1.05     2.435    0.0113  0.00094
  components 0.01446 + gear 0.01250 (D/q = 0.15 m^2) + cooling 0.00333 (D/q = 0.04 m^2) + misc 0.00151 (5%) = C_D0 0.03181; equivalent flat-plate area 0.382 m^2
span efficiency: lifting line (rectangular, AR 8.333) 0.9339; with fuselage factor 0.9113; Oswald e with lift-dependent profile drag 0.7073
best L/D = 12.06 at C_L = 0.768, V = 36.11 m/s
cruise 55 m/s: C_L = 0.3309, C_D = 0.03772 (induced 0.00591), L/D = 8.77, drag 838.7 N, drag power 46.13 kW
```

and writes `drag_build_up.png`: the polar with its tangent and cruise point, and the section estimate against Ladson's data.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-13-drag-build-up`. The tests check the friction formulas, their ordering and monotonicity; (1)'s limits, its fall with later transition and the Re_tr-only correction; the form factors; the perimeter of a 1 % thick section (≈ 1) and its growth with thickness; a cylinder's and a sphere's surface; drag proportional to S_wet·Q; the C_D0 sum and the effect of the gear term; e of a rectangle (0.9–0.99) and a taper-0.4 wing (> 0.98), and its fall with AR for a rectangle; the Oswald limits; (5) against a brute-force search with C_D = 2C_D0 at the optimum; L = W in cruise; the data blocks; the sign of the discrepancies; the inversion for x_tr.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** The data: Ladson, *Effects of Independent Variation of Mach and Reynolds Numbers on the Low-Speed Aerodynamic Characteristics of the NACA 0012 Airfoil Section*, NASA TM-4074 (1988), NTRS 19880019495; US Government work, public domain. A 60.1 cm chord NACA 0012 model spanning the Langley Low-Turbulence Pressure Tunnel at M = 0.15; cd from a wake-rake survey about one chord behind the trailing edge, corrected by the report for wall effects (about 2 %); repeatability of cd about ±0.0002. Table I is free transition; Table VII has a NACA-type wraparound strip of No. 60 grit from 0.05c on the upper surface round the leading edge to 0.05c on the lower surface. The transcription and its checks are in `curriculum/phase-10/data/README.md`. The values are used exactly as printed.

**Prediction, written before looking at the data:** with free transition in a low-turbulence tunnel the measured drag should lie *below* the fully turbulent estimate (laminar runs), and the gap should shrink as Re rises (transition moves forward); with the grit strip the drag should lie *above* even a turbulent estimate (roughness and grit drag).

1. **Free transition.** cd near α = 0 against the fully turbulent estimate (FF = 1.2524):

| Re | measured | turbulent estimate | difference | x_tr that matches |
|---|---|---|---|---|
| 2.00 × 10⁶ | 0.00620 | 0.00987 | +59.2 % | 0.464 |
| 3.94 × 10⁶ | 0.00659 | 0.00877 | +33.1 % | 0.279 |
| 5.97 × 10⁶ | 0.00650 | 0.00818 | +25.9 % | 0.217 |
| 8.86 × 10⁶ | 0.00630 | 0.00767 | +21.8 % | 0.179 |
| 1.19 × 10⁷ | 0.00634 | 0.00732 | +15.4 % | 0.128 |
| 1.89 × 10⁷ | 0.00655 | 0.00681 | +3.9 % | 0.032 |

The measured drag barely changes with Re (0.0062–0.0066) while the turbulent estimate falls by 31 %. The explanation is transition: inverting the estimate, the laminar run that matches the data shrinks from 46 % of the chord at 2 × 10⁶ to 3 % at 1.89 × 10⁷, i.e. transition at Re_x ≈ 0.6–1.6 × 10⁶. Lower friction from a shorter laminar run would lower the drag; the forward movement of transition raises it; the two nearly cancel. The rms discrepancy of the fully turbulent estimate is 31.6 %.

2. **Wraparound grit strip.** cd near α = 0:

| Re | measured | turbulent from LE | turbulent from 0.05c | grit increment |
|---|---|---|---|---|
| 2.00 × 10⁶ | 0.01100 | 0.00987 (−10.3 %) | 0.00950 (−13.6 %) | 0.00150 |
| 4.00 × 10⁶ | 0.00954 | 0.00875 (−8.3 %) | 0.00835 (−12.5 %) | 0.00119 |
| 6.00 × 10⁶ | 0.00895 | 0.00818 (−8.6 %) | 0.00777 (−13.2 %) | 0.00118 |
| 8.95 × 10⁶ | 0.00862 | 0.00766 (−11.1 %) | 0.00725 (−15.8 %) | 0.00137 |
| 1.21 × 10⁷ | 0.00842 | 0.00730 (−13.3 %) | 0.00690 (−18.1 %) | 0.00152 |

Every measurement is above the estimate, by 8–13 % for turbulence from the leading edge and 13–18 % for turbulence from the strip. The difference (0.0012–0.0015) is the drag of the roughness: the grit's own pressure drag and the thicker, rougher turbulent layer it produces, which a smooth-plate C_f does not contain. The rms discrepancy is 10.5 %.

3. **What it means for the build-up.** The trainer's wing (Re = 4.5 × 10⁶) used turbulent flow from 10 % chord. The data say that a smooth, clean wing could keep a laminar run of 20–30 % of the chord at that Re, and a rough or contaminated one (insects, rain, rivets) behaves more like the grit. Moving the wing's transition from 0.10 to 0.30 c changes C_D0 from 0.0318 to 0.0304 and (L/D)max from 12.1 to 12.4; fully turbulent gives 0.0326 and 11.9.

4. **The dominant input.** The gear's drag area: (D/q)_gear = 0, 0.15, 0.30 m² gives C_D0 = 0.0187, 0.0318, 0.0449 and (L/D)max = 16.5, 12.1, 9.7. A plausible ±0.05 m² on the gear moves C_D0 by about ±0.0044, about four times the spread of the whole transition question (0.0304–0.0326).

Raw data format: CSV tables `Re, alpha_deg, cd_measured, cd_turbulent, diff_pct, x_tr_match` (1); `Re, alpha_deg, cd_measured, cd_turb_LE, cd_turb_005, grit_increment` (2); `x_tr_wing, CD0, LD_max` (3); `gear_DQ, CD0, LD_max` (4).

## 5. Validate and improve

| quantity | estimate | data or independent route | difference | cause |
|---|---|---|---|---|
| 0012 cd, free transition, Re 2–19 × 10⁶ | turbulent 2C_fFF | Ladson Table I | +4 % to +59 % (rms 32 %) | laminar runs of 3–46 % c in the measurement |
| 0012 cd, wraparound grit, Re 2–12 × 10⁶ | turbulent 2C_fFF | Ladson Table VII | −8 % to −13 % (rms 10.5 %) | roughness and grit drag |
| best L/D | (5): 12.06 | brute-force search | < 10⁻⁸ | — |
| e of the wing alone | lifting line 0.934 | elliptic limit 1 | δ = 0.071 | rectangular planform |

The correlation cd = 2C_f(1 + 2t/c + 60(t/c)⁴) is consistent with the data *once the boundary-layer state is right*: with the matching transition locations it reproduces the free-transition drag exactly by construction, and the matching locations are physically sensible (decreasing with Re, Re_x,tr ≈ 10⁶). Its form factor is not what limits it; the transition location and the surface condition are. The engineering answer: **C_D0 ≈ 0.032 (±0.005, dominated by the gear and cooling add-ons), (L/D)max ≈ 12 at about 36 m/s, and about 46 kW of drag power at 55 m/s**; with a propeller efficiency of 0.8 that is about 58 kW of shaft power. To improve the estimate: measure or look up the gear's drag area for the actual installation, then predict the wing's transition from its pressure distribution (Extension 1).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1) and (5); the component table; `drag_build_up.png`; the two data tables of §4 with their explanation; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating C_D0 with its uncertainty and the input that dominates it.

### Adversarial self-checks

- **What if the speed doubles?** Every Re doubles and the turbulent C_f falls by about 10 %: the airframe part of C_D0 falls from 0.0145 to 0.0129 and C_D0 from 0.0318 to 0.0302; the add-ons, bluff bodies, do not change in coefficient. The drag rises roughly as V².
- **What if the wing were laminar to 30 %?** C_D0 falls from 0.0318 to 0.0304 (Experiment 3); the data show this needs a clean, smooth wing.
- **Which assumption fails first in the real world?** The gear and cooling drag areas; then the wing's transition location (insects, rain, paint).
- **What measurement could falsify the model?** A flight-test polar from glides or power-off descents (C_D0 and e from the intercept and slope of C_D against C_L²); the build-up should fall within its ±15 % if the add-ons are right.
- **Does the result survive a change of correlation?** The OpenVSP tool lists several friction and form-factor alternatives; they differ by a few per cent on the airframe part, small compared with the gear question (Extension 2).
- **Which decision changes if the uncertainty doubles?** Engine sizing for a 55 m/s cruise: ±10 kW of drag power instead of ±5 kW would push the choice to the next engine size.

## Check yourself

1. State the laminar and turbulent flat-plate friction laws and compare them at Re = 5 × 10⁵.
2. What does a form factor represent, and why does it grow with thickness?
3. Write the component build-up formula for C_D0 and name each term.
4. At the best L/D of a parabolic polar, how do the induced and zero-lift drags compare?
5. Why does the NACA 0012's free-transition drag stay nearly constant from Re = 2 × 10⁶ to 1.9 × 10⁷ while the turbulent estimate falls?
6. Why is the grit-strip drag higher than a turbulent smooth-plate estimate?
7. Which single input of the trainer's build-up has the largest effect on C_D0, and why can skin-friction theory not predict it?
