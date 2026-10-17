# P10.06  Viscous-Inviscid Coupling

> The panel method gives the velocity at the edge of the boundary layer; the boundary layer grows under it; its momentum deficit at the trailing edge is the drag. Couple the two and test the result against a wind tunnel.

| | |
|---|---|
| **Learning outcomes** | 1. Split a panel solution at the stagnation point and march a laminar boundary layer along each surface by Thwaites' method, detect laminar separation, and predict transition with Michel's criterion. 2. Integrate Head's turbulent entrainment method with an ODE solver and obtain the section's profile drag from the trailing-edge state by the Squire–Young formula. 3. Verify each piece against an independent solution (Falkner–Skan similarity solutions by a boundary-value solver; the 1/7-power-law flat plate) and validate the whole against Ladson's authentic NACA 0012 drag data. 4. Estimate the trainer's section drag at cruise, clean and tripped, and judge which part of the method limits its accuracy. |
| **Prerequisites** | `P10.05` (the vortex panel method, copied here as a given helper), `P08.11` (integral boundary layers and Thwaites' method; recapped in §2.2), `P08.12` (transition; recapped in §2.3), `P08.14` (turbulent boundary layers; recapped in §2.4), `P03.05` and `P03.06` (ODE integrators, `solve_ivp`), `P03.03` (cubic splines). Skill check: you can state the momentum integral equation and define θ, δ* and H. |
| **Engineering question** | What is the profile drag coefficient of the trainer's NACA 2412 section at cruise (Re = 4.52 × 10⁶, c_l = 0.331), clean and with transition fixed near the leading edge by insects or rain, and how well does the same method predict the measured drag of a NACA 0012? Answers: c_d = 0.0064 clean (transition at 36 % chord on the upper surface, 26 % on the lower), 0.0087 tripped at 5 % chord: 142 N and 193 N of wing profile drag, 2.8 kW of extra power. For the NACA 0012 at zero lift the prediction is within −6 % to +11 % of Ladson's free-transition data for Re = 2–9 × 10⁶, and 19 % low at 1.9 × 10⁷. |
| **Assumptions and limits** | One-way coupling (inviscid u_e drives the boundary layer; an optional single displacement-thickness iteration feeds it back); attached flow, thin boundary layers; Thwaites' one-parameter laminar closure; Michel's empirical transition correlation, which knows nothing of free-stream turbulence or roughness; instant transition at laminar separation (no bubble model); Head's two-equation turbulent closure with Ludwieg–Tillmann friction; Squire–Young for the wake. Fails at low Re (bubbles, `P10.14`), near stall (separated flow, `P10.15`), and gives only profile drag (no induced drag, `P10.09`). |
| **Analytical work** | The momentum integral equation; Thwaites' quadrature and its stagnation-point limit θ² = 0.075ν/u_e′; λ, H(λ), ℓ(λ); Michel's criterion on a flat plate; Head's entrainment equations; Squire–Young and its limit 2θ/c; the flat-plate values θ√Re_x/x = 0.664 (Blasius) and √0.45 = 0.671 (Thwaites). |
| **Computational deliverable** | `starter/main.py`: `falkner_skan`, `thwaites_H_l`, `thwaites`, `michel_transition`, `head_H1`, `head_H_of_H1`, `ludwieg_tillmann_cf`, `head_turbulent`, `split_at_stagnation`, `surface_boundary_layer`, `squire_young`, `section_drag`, `displacement_iteration` (`naca4`, `panel_solve` from `P10.05`, `find_data_file`, `read_ladson`, `ladson_zero_lift_drag`, `trim_branch`, `observed_order` and `main` are given). Tests prove: Falkner–Skan reproduces Blasius (f″(0) = 0.4696, θ√Re_x/x = 0.664, H = 2.59) and Hiemenz; Thwaites gives √0.45 on a flat plate, within 1.5 % of Blasius, and a constant θ within 7 % of Hiemenz in stagnation flow; the correlations are continuous and flag separation; Michel puts flat-plate transition between Re_x = 2 × 10⁶ and 10⁷; Head matches the 1/7-power law to 5 % and separates in a strong adverse gradient; Squire–Young limits; a symmetric section splits symmetrically and its drag falls with Re; a trip raises drag; the displacement thickness reduces lift by 0–6 %; the NACA 0012 prediction is within 25 % of Ladson's free-transition data for Re ≤ 1.25 × 10⁷. |
| **Practical activity** | **computational experiment, compared with the authentic NACA 0012 data** of Ladson (NASA TM-4074) in `curriculum/phase-10/data/ladson_naca0012_M015.csv`: zero-lift drag against Re, free and fixed transition; panel-count sensitivity; the trainer's polar; the displacement-thickness iteration. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `viscous_drag.png` with measured points labelled as measured, the derivations of §2.2–2.5, the experiment tables, and a one-page note comparing prediction and measurement with an explanation of the discrepancy. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: Thwaites' quadrature and the stagnation limit derived; all public tests pass; each closure verified against its independent check; the Ladson comparison tabulated with percentage differences and the high-Re discrepancy explained by transition. |
| **References** | `@mit-ocw-16-100-aerodynamics` Lecture 22 "Integral Boundary Layer Equations", Lecture 23 "Correlation Methods for Integral Boundary Layers", Lecture 25 "Falkner-Skan Flows" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre22_cg/, …/16100lectre23/, …/16100lectre25/) (required-free). `@mit-ocw-16-13-aerodynamics-viscous-fluids` handouts "Thwaites' Method", "Integral Thicknesses", "Displacement Effects of Boundary Layer on Potential Flow"; Lecture 17 "2D Interaction Models: Displacement Body, Transpiration. Form Drag, Stall Mechanisms."; Lecture 23 "Transition Mechanisms. Transition Prediction: Local Correlations, Amplification Methods." (https://ocw.mit.edu/courses/16-13-aerodynamics-of-viscous-fluids-fall-2003/resources/thwaites_method/, …/lecture_17/, …/lecture_23/) (supplementary-free). `@leishman-intro-aerospace-flight-vehicles` "Boundary Layer Flows": Pressure Gradients & Flow Separation; von Kármán Momentum Integral Equation; Skin Friction Drag on an Airfoil (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/) (required-free, link only). `@drela-xfoil` User Primer "Viscous Formulation" (https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) (supplementary-free: a strongly coupled production method). Data: `@ntrs-nasa-tm-4074-ladson-naca0012` (https://ntrs.nasa.gov/citations/19880019495), Tables I and VII, transcribed in `curriculum/phase-10/data/` (required-free; read that folder's README). `@scipy-userguide` "Integration (scipy.integrate)" (supplementary-free: `solve_ivp`, `solve_bvp`). Optional textbooks for deeper study: M. Drela, *Flight Vehicle Aerodynamics* (MIT Press; paid, library), integral boundary layers and viscous–inviscid coupling at graduate depth. |
| **Extensions** | ★ Replace Michel's criterion by a simple e^N envelope method and fit N to Ladson's free-transition drag. ★ Iterate the displacement-thickness coupling to convergence with under-relaxation and watch it fail near stall. ★ Add a wake to the panel model and compute the drag from the far-wake momentum deficit instead of Squire–Young. |
| **Time** | unmeasured |

## 1. The engineering problem

A panel method (`P10.05`) predicts lift well and drag not at all: d'Alembert's paradox makes the inviscid drag zero. Yet section drag decides the cruise power of an aeroplane, the glide ratio of a sailplane and the efficiency of a propeller. It comes from a layer of air a few millimetres thick on the surface, the **boundary layer**, where viscosity slows the flow. Prandtl's insight was that for a streamlined body at high Reynolds number this layer is thin, so the outer flow is nearly the inviscid one, and the boundary layer can be computed afterwards, driven by the inviscid edge velocity u_e(s). That is **viscous-inviscid coupling**, and its simplest form, one-way, is what this lesson builds.

The chain is: panel solution → split at the stagnation point → laminar boundary layer (Thwaites) → transition (Michel) → turbulent boundary layer (Head) → drag from the trailing-edge state (Squire–Young). Each link is a correlation or a model, and each is checked against an independent solution before the chain is trusted. Then the chain is tested against a measurement: the NACA 0012 drag measured by Ladson in the NASA Langley Low-Turbulence Pressure Tunnel (NASA TM-4074, 1988; authentic, US Government work, public domain; transcribed in `curriculum/phase-10/data/` with its provenance).

The trainer's question is the drag of its NACA 2412 section at cruise (Re = 4.52 × 10⁶): clean, and with the leading edge contaminated by insects or rain so that the boundary layer goes turbulent near the nose.

Success criteria, written first: Thwaites reproduces √0.45 on a flat plate and agrees with Blasius to 1.5 %; Head reproduces the turbulent flat plate to 5 %; a symmetric section at zero lift has equal upper and lower drag; and the NACA 0012 zero-lift drag is predicted within 25 % of the measured values for Re = 2–12 × 10⁶, with the trend explained.

## 2. Understand and derive

### 2.1 The quantities

At a station s along the surface, the velocity rises from zero at the wall to u_e at the edge. Two integral thicknesses summarise the profile: the **displacement thickness** δ* = ∫(1 − u/u_e) dy, how far the wall is effectively pushed out, and the **momentum thickness** θ = ∫(u/u_e)(1 − u/u_e) dy, the momentum deficit. Their ratio is the **shape factor** H = δ*/θ: about 2.6 for a laminar flat plate, 1.3–1.4 for a turbulent one, rising towards separation. Integrating the boundary-layer equations across the layer gives the **momentum integral equation** (`P08.11`):

**dθ/ds + (2 + H)(θ/u_e) du_e/ds = c_f/2.**   (1)

It has three unknowns (θ, H, c_f): every integral method is a choice of closure.

### 2.2 Laminar: Thwaites' method (recap of `P08.11`)

Multiplying (1) by 2u_eθ/ν and defining **λ = θ²u_e′/ν** gives (u_e/ν) d(θ²)/ds = F(λ). Thwaites observed that for many exact laminar solutions F ≈ 0.45 − 6λ. Then (Problem 2) the equation integrates in closed form:

**θ² u_e⁶ = 0.45 ν ∫₀^s u_e⁵ ds,**   (2)

with H(λ) and ℓ(λ) = c_f u_eθ/(2ν) from fits to Thwaites' tables, and **laminar separation at λ ≈ −0.09**. At a stagnation point (u_e = Cs) the integral gives a constant **θ² = 0.075ν/C**; the code uses that limit at s = 0. On a flat plate (λ = 0), θ√Re_x/x = √0.45 = **0.6708**, H = 2.61. The independent check is the Falkner–Skan similarity equation f‴ + ff″ + β(1 − f′²) = 0, solved here by `scipy.integrate.solve_bvp`: Blasius (β = 0) gives f″(0) = 0.46960 and θ√Re_x/x = **0.66411**, H = 2.5911; Hiemenz (β = 1) gives θ = 0.29234√(ν/C), H = 2.2162. Thwaites is 1.0 % high on the plate and 6.3 % low at a stagnation point.

The integral in (2) is taken exactly for u_e linear between stations, ∫u_e⁵ ds = Δs(u₀⁵ + u₀⁴u₁ + … + u₁⁵)/6; the trapezoid rule would be three times too large in the first interval after the stagnation point, where u_e ∝ s. (The first version of the code used the trapezoid rule, and the stagnation-flow test found it.)

### 2.3 Transition: Michel's criterion (recap of `P08.12`)

A laminar layer becomes turbulent when small disturbances grow; the local correlation of Michel says this happens when

**Re_θ ≥ 1.174 (1 + 22 400/Re_x) Re_x^{0.46},**   (3)

with Re_x = u_e s/ν measured from the stagnation point. On a flat plate (Blasius) it is met at Re_x = 2.02 × 10⁶: 0.54 m from the leading edge at 55 m/s. A favourable pressure gradient (accelerating flow near the nose) keeps Re_θ low and delays transition; an adverse one advances it. If λ reaches −0.09 first, the code puts transition at the laminar separation point (it assumes a short bubble that reattaches turbulent).

### 2.4 Turbulent: Head's method (recap of `P08.14`)

Head's method adds an equation for the rate at which the turbulent layer entrains outer fluid. With the entrainment shape factor H₁ = (δ − δ*)/θ,

**d(u_e θ H₁)/ds = u_e · 0.0306 (H₁ − 3)^{−0.6169},**   (4)

closed with the fits H₁(H) and the **Ludwieg–Tillmann** friction law c_f = 0.246 × 10^{−0.678H} Re_θ^{−0.268}. Equations (1) and (4) are integrated for y = [θ, u_eθH₁] with `solve_ivp` (RK45, `P03.05`) from the transition point, with θ continuous and H = 1.4, and u_e and u_e′ from a cubic spline. Turbulent separation is taken at H = 2.4 (an event stops the integration). The independent check is the 1/7-power-law flat plate, θ = 0.036 x Re_x^{−0.2} and c_f = 0.0592 Re_x^{−0.2}: at Re_L = 10⁷ Head gives θ = 0.001409 L against 0.001433 L and c_f = 0.00234 against 0.00236, with H = 1.343.

### 2.5 Drag: the Squire–Young formula

The drag is the momentum deficit far downstream, θ_∞: D′ = ρV²θ_∞. The wake between the trailing edge and "far downstream" still accelerates from u_e,TE to V, and its θ falls as it does. Squire and Young integrated that with an assumed H–u_e relation:

**c_d = 2 (θ_TE/c) (u_e,TE/V)^{(H_TE + 5)/2},**   (5)

per surface, summed. With u_e,TE = V it is just 2θ/c.

**A closed-trailing-edge detail.** The inviscid flow decelerates sharply into a finite-angle trailing edge over the last 1 % of chord (`P10.05`, §2.6), an artefact the real, viscous flow does not see. The code stops the boundary layer at x/c = 0.99 and applies (5) there (the given `trim_branch`). Without that, H blew up in the last two panels and the method reported turbulent separation at α = 0°.

### 2.6 One step back: displacement thickness

The boundary layer displaces the outer flow by δ*. Adding δ* to the surface along its normal and re-running the panel method once (the "displacement body") shows how much the viscous layer decambers the section: it takes 2.1 % of the trainer's cruise lift.

### 2.7 Where the model stops

- **Transition** is the weakest link: Michel's criterion is a correlation, blind to free-stream turbulence, noise and roughness (§4 shows what it costs).
- **Low Re**: laminar separation bubbles are not modelled (`P10.14`).
- **High α**: turbulent separation near the trailing edge, then stall: one-way coupling cannot predict the lift loss (`P10.15`).
- **Profile drag only**: no induced drag (`P10.09`), no compressibility (Phase 11).

## 3. Build the computational model

Open `starter/main.py`. The NACA generator, the compact panel method from `P10.05` (it returns the signed surface velocity at the nodes, V_t = −γ), the Ladson reader and `trim_branch` are given.

### Step 1: the independent checks first

`falkner_skan(beta)` solves the similarity equation with `solve_bvp` on η ∈ [0, 10] and returns f″(0), I_θ, I_δ and H. Do this before Thwaites: it is your reference.

### Step 2: laminar and transition

`thwaites_H_l(lam)` is the pair of fits; `thwaites(s, ue, nu)` is (2) with the exact-for-linear integral and the stagnation limit, returning θ, λ, H, c_f, Re_θ and the separation index. `michel_transition(Re_x, Re_theta)` is (3).

### Step 3: turbulent and drag

`head_H1`, `head_H_of_H1` and `ludwieg_tillmann_cf` are the closures; `head_turbulent(s, ue, nu, s0, theta0)` integrates (1) and (4):

```python
def rhs(x, y):
    th, Q = y
    u = spl(x); H1 = Q / (u * th); H = head_H_of_H1(H1)
    cf = ludwieg_tillmann_cf(H, max(u * th / nu, 10.0))
    return [cf / 2 - (H + 2) * th / u * dspl(x), u * 0.0306 * max(H1 - 3.0, 1e-6) ** -0.6169]
```

`squire_young` is (5).

### Step 4: coupling

`split_at_stagnation(x, y, Vt)` finds the node interval where V_t changes sign from − to +, interpolates the stagnation point, and returns two branches with s from it. `surface_boundary_layer` chains Thwaites, Michel (or laminar separation, or a forced trip) and Head. `section_drag(code, alpha, Re, n, trip_x)` runs everything for both surfaces; `displacement_iteration` thickens the body once.

Running `main()` prints (abridged)

```
Blasius (Falkner-Skan beta = 0): f''(0) = 0.46960, theta sqrt(Re_x)/x = 0.66411, H = 2.5911; cf sqrt(Re_x) = 0.66411
Thwaites flat plate: theta sqrt(Re_x)/x = 0.67082 (sqrt 0.45 = 0.67082), H = 2.610, cf sqrt(Re_x) = 0.6559
Hiemenz (beta = 1): theta sqrt(C/nu) = 0.29234, H = 2.2162, f''(0) = 1.23259
Thwaites stagnation flow: theta sqrt(C/nu) = 0.27386 (sqrt 0.075 = 0.27386), lambda = 0.0750, H = 2.358
Head flat plate (Re_L = 1e7, start x = 0.010): theta(L) = 0.001409, 1/7-power law 0.036 L Re_L^-0.2 = 0.001433; H(L) = 1.343; cf(L) = 0.00234 vs 0.0592 Re_x^-0.2 = 0.00236
Ladson (NASA TM-4074) Table I, free transition, |alpha| <= 0.15 deg (ladson_naca0012_M015.csv):
   Re = 5.97e+06: measured cd = 0.00650; predicted 0.00617 (-5.1 %)
trainer NACA 2412 at cruise: Re = 4.52e+06, inviscid cl = 0.3311 at alpha = 0.592 deg: cd = 0.00637 (upper 0.00367, transition x/c 0.360 michel; lower 0.00270, x/c 0.261 michel)
   tripped at 5 % chord (insects, rain, rivets): cd = 0.00866
   one displacement-thickness iteration: cl 0.3311 -> 0.3242 (-2.1 %)
```

and writes `viscous_drag.png`: c_d against Re for the model (free and tripped) and Ladson's measured points (marked "measured"), and θ and H along the upper surface of the NACA 0012.

### Verify

`python -m pytest curriculum/phase-10/lesson-06-viscous-inviscid-coupling` runs 15 tests (listed in the metadata). The Ladson test reads the data file only if it is present (it is, read-only) and is skipped otherwise.

### Compare

Each closure has its own independent check (Falkner–Skan, the power-law plate); the coupled method is compared with measurement. XFOIL's viscous mode (optional, `P10.21`) solves the same problem with strong coupling and an e^N transition model.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment compared with authentic data.** The measured values are Ladson's (NASA TM-4074, Langley LTPT, 2-D model of 0.601 m chord, M = 0.15, c_d from a wake rake, data wall-corrected by the report); read `curriculum/phase-10/data/README.md` before using them. Every other number here is computed. Write your predictions first.

1. **NACA 0012 at zero lift, free transition.** Predict whether the model's c_d will rise or fall with Re. *Finding:*

| Re | measured c_d | predicted c_d | difference | predicted transition x/c |
|---|---|---|---|---|
| 2.00 × 10⁶ | 0.00620 | 0.00690 | +11.2 % | 0.391 |
| 3.94 × 10⁶ | 0.00659 | 0.00649 | −1.5 % | 0.316 |
| 5.97 × 10⁶ | 0.00650 | 0.00617 | −5.1 % | 0.287 |
| 8.86 × 10⁶ | 0.00630 | 0.00590 | −6.3 % | 0.259 |
| 1.19 × 10⁷ | 0.00635 | 0.00567 | −10.8 % | 0.245 |
| 1.89 × 10⁷ | 0.00660 | 0.00536 | −18.7 % | 0.219 |

The model falls steadily with Re (thinner turbulent layer, slowly forward-moving transition); the measurement stays at 0.0062–0.0066.
2. **Fixed transition.** Ladson's Table VII used a NACA wraparound grit strip from 5 % upper round the nose to 5 % lower. *Finding:* measured 0.01096, 0.00953, 0.00895, 0.00855, 0.00835 for Re = 2.0–12.1 × 10⁶; the model tripped at 5 % chord gives 0.00974, 0.00860, 0.00799, 0.00743, 0.00704, 10–16 % low: the grit makes the layer turbulent from the nose and adds roughness drag the model does not have (Problem 4).
3. **Where must transition be?** At Re = 1.19 × 10⁷, trips at 5, 10 and 20 % give 0.00706, 0.00670, 0.00596 against the measured 0.00635: transition near 14 % chord, ahead of Michel's 24.5 % (Problem 6).
4. **Panel count.** NACA 0012, Re = 5.97 × 10⁶: c_d = 0.006146, 0.006171, 0.006182 with 100, 200, 400 panels; the numerical error (0.6 %) is far below the modelling error.
5. **The trainer's polar** at Re = 4.52 × 10⁶: c_d = 0.00633, 0.00664, 0.00731, 0.00828, 0.01029 at α = 0°, 2°, 4°, 6°, 8° (c_l = 0.260, 0.501, 0.742, 0.981, 1.220 inviscid), upper transition moving from 0.391c to 0.022c (by laminar separation at 8°, where Head also predicts turbulent separation before the trailing edge). One displacement iteration takes 1.7 %, 3.1 %, 4.8 %, 7.9 % off c_l at 0°–6°.

Raw data format: CSV `Re, cd_measured, cd_predicted, diff_percent, x_tr` (1); `Re, cd_measured_grit, cd_predicted_trip05` (2); `trip_x, cd` (3); `alpha_deg, cl, cd, x_tr_upper, x_tr_lower, separation` (5); every measured column headed with its source.

## 5. Validate and improve

| check | reference | method | difference |
|---|---|---|---|
| flat plate θ√Re_x/x | Blasius (solve_bvp): 0.66411 | Thwaites: 0.67082 | +1.0 % |
| stagnation θ√(C/ν) | Hiemenz: 0.29234 | Thwaites: 0.27386 | −6.3 % |
| turbulent plate θ(L), Re_L = 10⁷ | 1/7 power law: 0.001433 | Head: 0.001409 | −1.7 % |
| NACA 0012 c_d, Re = 2–8.9 × 10⁶ | Ladson, measured | predicted | −6.3 % to +11.2 % |
| NACA 0012 c_d, Re = 1.89 × 10⁷ | Ladson: 0.00660 | 0.00536 | −18.7 % |
| tripped c_d, Re = 2–12 × 10⁶ | Ladson grit: 0.0110–0.0083 | 0.0097–0.0070 | −10 to −16 % |

The pieces are verified to a few per cent; the chain is validated to about 10 % up to Re ≈ 10⁷, with a discrepancy that grows with Re and is explained by transition: matching the data at 1.19 × 10⁷ needs transition at 14 % chord, not 24.5 %. The engineering answer: **the trainer's NACA 2412 has c_d ≈ 0.0064 at cruise with natural transition (36 % and 26 % chord), 0.0087 with the boundary layer tripped at 5 %; over 12 m² at 55 m/s that is 142 N or 193 N of profile drag, and the 51 N difference costs 2.8 kW.** Given the comparison with Ladson, quote the clean value as 0.0064 ± 10 %.

What to improve: transition (an e^N method calibrated to the tunnel, Extension 1); then a strongly coupled iteration for high α. What not to improve: the panel count.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), the stagnation limit and (5); `viscous_drag.png` with measured points labelled as measured and their source; the five experiment tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1–7) on the trainer's section drag with its uncertainty, and the comparison with Ladson's data.

### Adversarial self-checks

- **What if Re doubles?** Laminar θ falls as Re^{−1/2}, turbulent as about Re^{−0.2}, and transition moves forward: the model's c_d falls by 8 % from 5.97 to 11.9 × 10⁶; the measured one does not, which is the model's weakness, not a numerical error.
- **What if V → 0 (Re → small)?** Laminar separation comes first everywhere (Problem 3: the UAV at 1.6 × 10⁵) and the instant-transition assumption is no longer defensible.
- **Which assumption fails first?** Michel's criterion in a real environment; then attached flow near the trailing edge at high α.
- **What measurement could falsify the model?** Measured transition positions (hot films, oil flow) at the same Re: if the measured transition at 1.19 × 10⁷ is near 25 % and the drag is still 0.00635, the turbulent closure, not the transition, is wrong.
- **Does the result survive a change of grid or solver?** Panel count changes c_d by 0.6 % from 100 to 400 panels; `solve_ivp`'s tolerance (10⁻⁸) by far less.
- **Which decision changes if the drag uncertainty doubles to ±20 %?** The cruise-power estimate (±1.5 kW on the wing's profile drag alone) would no longer separate a clean laminar wing from a tripped one; a flight test or a tunnel test of the real section would be needed before claiming the laminar benefit.

## Check yourself

1. Define the momentum thickness and the shape factor.
2. State Thwaites' integral and the criterion for laminar separation.
3. State Michel's transition criterion. At what Re_x does it put transition on a flat plate?
4. What does Head's method add to the momentum integral, and what closes it?
5. Write the Squire–Young formula.
6. Compare θ√Re_x/x on a flat plate from Thwaites and from Blasius.
7. How well does the method predict Ladson's NACA 0012 drag, and why does the agreement worsen at high Re?
