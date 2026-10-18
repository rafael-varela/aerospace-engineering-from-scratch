# P10.07  Airfoil Design & Selection

> Thin-airfoil theory runs both ways: give it a camber line and it returns the load; give it the load you want and it returns the camber line. Choosing a section is a different problem, a weighed judgment, and it deserves a sensitivity study as much as any calculation.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the thin-airfoil inverse problem: from a target chordwise load Δcp(x) to the Glauert sine coefficients, the camber slope as the matching cosine series, and the ideal angle that closes the camber line. 2. Show by hand and in code that a uniform load gives the NACA a = 1.0 mean line, and that a target (cl, cm_c/4) gives a cubic camber line. 3. Implement and test the inversion, verify every recovered camber line by a forward thin-airfoil analysis, wrap a thickness distribution around it and check it with a panel method. 4. Build a weighted decision matrix over candidate sections from criteria your own models compute, and test the choice with a Monte Carlo over the weights and the exact weights at which the winner changes. 5. State what thin-airfoil inverse design cannot do. |
| **Prerequisites** | `P10.04` (Glauert's substitution, A0, A1, A2, cl and cm from the coefficients), `P10.02` (NACA mean lines and thickness, cosine spacing), `P10.05` (panel method), `P02.09` (Fourier series and orthogonality), `P01.09` (optimisation and trade-offs), `P03.14` (sensitivity of an optimum). Skill check: you can compute (2/π)∫₀^π f(θ) sin nθ dθ for a simple f by hand. |
| **Engineering question** | What camber line carries a uniform load Δcp = 0.4 along the whole chord, how far does adding 12 % thickness move it from its design, and which of seven candidate sections best suits the reference trainer (and, for contrast, a sailplane), and how robust is that choice? Answer: the inversion recovers the NACA a = 1.0 mean line to 1.2 × 10⁻⁶ c (maximum camber 2.206 % at mid-chord, ideal angle 0°); with 12 % thickness the panel method gives cl = 0.433 at the design angle, 8 % above the design, and a load that differs from uniform by 0.053 on average. On the lesson's criteria the NACA 4415 scores 0.818 against 0.623 for the trainer's NACA 2412, and wins 99.4 % of random weightings; the 2412 wins only if the pitching-moment weight exceeds 0.54. |
| **Assumptions and limits** | Thin-airfoil theory: small camber and angle, zero thickness, inviscid, incompressible, Kutta condition at the trailing edge; the target load must be finite at the leading edge (the design point is then the ideal angle). Thickness enters only through the panel check. No boundary layer anywhere: the selection criteria are inviscid proxies (no cd, no cl_max, no transition); the decision-matrix weights are the author's judgment, stated and tested, not data. |
| **Analytical work** | The inverse thin-airfoil relation dz/dx = α_i + Σ A_n cos nθ; the closure integral I_n = (1 + cos nπ)/(1 − n²); the camber line term by term; the uniform-load series and its closed form z = −(cl_i/4π)[x ln x + (1 − x) ln(1 − x)]; the two-term design A1 = cl/π, A2 = A1 + 4cm/π, α_i = A2/3; the min–max normalisation and the linear dependence of the totals on one weight. |
| **Computational deliverable** | `starter/main.py`: `glauert_load`, `load_coefficients`, `ideal_angle`, `camber_from_coefficients`, `naca_a1_camber`, `camber_from_cl_cm`, `load_cl_cm`, `inverse_design`, `thicken`, `panel_load`, `section_metrics`, `decision_matrix`, `weight_monte_carlo`, `weight_flip_points` (the NACA generator, `thin_airfoil`, the panel solver, `alpha_for_cl`, `mission_cls`, `build_candidates`, `select_sections` and `main` are given). Tests prove: the load series round-trips; the ideal angle closes the camber line; a uniform load recovers the a = 1.0 line at observed order 2; the forward analysis of each recovered line reproduces its load, cl and cm; the inversion is linear; `thicken` reproduces the NACA 2412; the panel cl tends to the thin-airfoil value as t → 0; the decision matrix is normalised, scale-invariant and never picks a dominated section; flip points match brute force. |
| **Practical activity** | **computational experiment**: inverse design of three target loads with forward verification, a series-truncation convergence study, a thickness sweep with the panel method, and a selection study for two missions with a weight sensitivity analysis. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `airfoil_design.png`, the derivations of the inverse relation, the closure integral and the a = 1.0 line, the convergence and thickness tables, the two decision matrices with their sensitivity results, and a half-page selection note for the trainer. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the inverse relation is derived; all public tests pass; the a = 1.0 recovery is verified forward and its observed order reported; the selection note gives the winner, the weight at which it changes, and two things the criteria leave out. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` "Classic Airfoil Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/classic-airfoil-theory/): Thin Airfoil Theory; Inverse Airfoil Design; and "Aerodynamics of Airfoil Sections" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airfoil-characteristics/): Inverse Airfoil Design and Shape Optimization; Drag Polars (required-free). `@mit-ocw-16-01-unified-engineering-fluids` Fluids (fall) Lecture 21 "Airfoil Polar Relations" and Lecture 20 "Airfoils: Detailed Look" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f21_fall/) (required-free): how a mission's cl range is read against a polar. `@ntrs-naca-tr-824-summary-of-airfoil-data` "Application to Wing Design": Application of Section Data; Selection of Root Section; Selection of Tip Section (https://ntrs.nasa.gov/citations/19930090976) (required-free): the NACA engineers' own selection logic. `@drela-xfoil` User Primer "Inverse Formulation", "Full-Inverse Design Routine (MDES)", "Mixed-Inverse Design Routine (QDES)" (https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) (supplementary-free): inverse design with thickness and viscosity, the step beyond this lesson. Optional textbooks for deeper study (optional-paid, not required; library or used copy): Abbott & von Doenhoff, *Theory of Wing Sections*, for the NACA a-series mean lines in tabulated form; Anderson, *Fundamentals of Aerodynamics*, for thin-airfoil theory worked at length. |
| **Extensions** | ★ Add the NACA 6-series thickness and test your a = 0.6 roof-top camber in the panel method at several angles: where does the roof-top survive? ★ Replace the sine series by a least-squares fit to a target Δcp given only at a few stations, with a smoothness penalty, and see how the camber line responds to noise. ★ Add a viscous criterion from `P10.06` (cd at the design cl) to the decision matrix and redo the sensitivity study. |
| **Time** | unmeasured |

## 1. The engineering problem

An airfoil designer usually knows what the section must *do* before knowing what it looks like: carry a given lift at cruise with the flow entering the leading edge smoothly, keep the pitching moment small so the tail does not have to fight it, and spread the load so that the pressure recovery near the trailing edge is gentle. Turning such a wish list into a shape is **inverse design**. The analysis problem of `P10.04` (camber line → load) has a clean inverse in thin-airfoil theory, and working it out shows what a camber line *is*: the shape that makes a vortex sheet with a chosen strength satisfy flow tangency.

The other half of the job is **selection**: picking, from sections that already exist, the one that best serves a mission. Real selection weighs many criteria that do not reduce to one number: thickness (spar depth and fuel volume), the cl at which the section is happiest (its ideal cl), the pitching moment (trim drag and wing torsion), and stall behaviour. Engineers often formalise this as a **weighted decision matrix**. The danger is that the weights carry the decision while looking objective, so the matrix is only as good as its sensitivity study.

This lesson does both for the reference trainer (750 kg, S = 12 m², cruise 55 m/s, approach 30 m/s, so cl ≈ 0.331 in cruise and 1.112 on approach) and, for contrast, a sailplane whose design cl is 0.70 and thermalling cl 1.20 (invented but typical values). Success criteria, written first: the inversion of a uniform load matches the NACA a = 1.0 closed form to better than 2 × 10⁻⁶ c with 200 terms and converges at a measured order; every recovered camber line, analysed forward, reproduces its target load, cl and cm; the panel method agrees with thin-airfoil cl to 1 % at 2 % thickness; the decision matrix reports a winner *and* the weight at which it would change.

## 2. Understand and derive

### 2.1 Thin-airfoil theory in both directions

Recall `P10.04`. Replace the section by a vortex sheet γ(x) on the chord (0 ≤ x ≤ c, take c = 1). Flow tangency on the camber line z(x) requires the sheet's induced downwash to cancel the free stream's normal component:

(1/2π) ∫₀¹ γ(ξ) dξ/(x − ξ) = V∞ (α − dz/dx).

With x = (1 − cos θ)/2, Glauert's solution is

γ(θ)/(2V∞) = A0 (1 + cos θ)/sin θ + Σₙ Aₙ sin nθ,   dz/dx = (α − A0) + Σₙ Aₙ cos nθ.   (1)

The pressure difference across the sheet is the load: **Δcp = cp,lower − cp,upper = 2γ/V∞ = 4[A0 cot(θ/2) + Σ Aₙ sin nθ]** (positive when the section lifts). Integrating, cl = 2π(A0 + A1/2) and cm_c/4 = (π/4)(A2 − A1).

Equation (1) is the whole inverse method. The load is a **sine** series in θ; the camber slope is the **cosine** series with the same coefficients. In the analysis direction we compute Aₙ from dz/dx by cosine projection; in the design direction we compute Aₙ from Δcp by sine projection:

**Aₙ = (2/π) ∫₀^π (Δcp/4) sin nθ dθ.**   (2)

The A0 term is the leading-edge singularity (γ ~ 1/√x). A designer wants the flow to enter the nose smoothly, so we design for loads that are finite at the leading edge: A0 = 0. That fixes the angle at which the design load occurs, the **ideal angle** α_i, and (1) becomes

**dz/dx = α_i + Σ Aₙ cos nθ.**   (3)

### 2.2 Closing the camber line

The camber line must start and end on the chord: z(0) = z(1) = 0, so ∫₀¹ dz/dx dx = 0. With dx = ½ sin θ dθ,

0 = ∫₀^π (α_i + Σ Aₙ cos nθ) sin θ dθ = 2α_i + Σ Aₙ Iₙ,  Iₙ = ∫₀^π cos nθ sin θ dθ = (1 + cos nπ)/(1 − n²) (I₁ = 0).

So **α_i = −½ Σ Aₙ Iₙ**. Only even n contribute (Iₙ = −2/(n² − 1) for even n, 0 for odd n): a load symmetric about mid-chord (odd n only) is carried at α_i = 0. The camber itself follows term by term from sin θ cos nθ = ½[sin(n + 1)θ − sin(n − 1)θ]:

z(θ) = ½[α_i(1 − cos θ) + Σ Aₙ Jₙ(θ)],  Jₙ = ½[(1 − cos(n + 1)θ)/(n + 1) − (1 − cos(n − 1)θ)/(n − 1)],  J₁ = ½ sin²θ.   (4)

Units: Δcp, Aₙ, α_i and z/c are all dimensionless (angles in radians).

### 2.3 Worked example: the uniform load gives the a = 1.0 mean line

Let Δcp = cl_i everywhere. Then (2) gives Aₙ = (2/π)(cl_i/4)(1 − cos nπ)/n = **cl_i/(πn) for odd n**, zero for even n, so α_i = 0 and

dz/dx = (cl_i/π) Σ_{n odd} cos nθ / n = (cl_i/π)(−½ ln tan(θ/2)) = (cl_i/4π) ln((1 − x)/x),

using the classical series Σ_{n odd} cos nθ/n = −½ ln tan(θ/2) and tan²(θ/2) = x/(1 − x). Integrating with z(0) = 0:

**z = −(cl_i/4π)[x ln x + (1 − x) ln(1 − x)],**   (5)

the NACA a = 1.0 mean line. Its maximum, at x = ½, is (cl_i/4π) ln 2: **0.02206 c for cl_i = 0.4**. Its slope is infinite at both ends (logarithmically), the price of a load that stays finite right up to the leading and trailing edges. The centre of pressure is at mid-chord, so cm_c/4 = −cl_i/4 = −0.1, a large nose-down moment: uniform loading is great for keeping the suction peak low and bad for trim.

The truncation error of the series is easy to predict: the tail Σ_{n > N} Aₙ Jₙ has terms of size 1/n × 1/n, so the error in z falls as **1/N²**. The code measures it.

### 2.4 Designing for cl and cm

The simplest load with a chosen cl and cm_c/4 at the ideal angle has two terms: A0 = 0, **A1 = cl/π, A2 = A1 + 4cm_c/4/π**, and α_i = −½A2 I₂ = **A2/3**. Equation (4) then gives a cubic camber line. For cl = 0.30 and cm = −0.053 (the NACA 2412's thin-airfoil moment): A1 = 0.0955, A2 = 0.0280, α_i = 0.535°, maximum camber 2.47 % at 41 % chord, close to the 2412's 2 % at 40 %. The 2412's own ideal cl is 0.256: the cubic line is a smooth relative of the 2412 mean line, not a copy, because the 2412's line has a kink in its curvature at 40 % chord that two terms cannot represent.

### 2.5 What the thin-airfoil inverse cannot do

- **Thickness.** The theory has none. Thickness adds lift (about 8 % at t/c = 0.12 here), and it adds a symmetric velocity increment that changes the surface Cp everywhere, so the design *load* survives better than the design *pressures*. A real roof-top pressure distribution (for laminar flow) is a statement about the upper-surface speed, which thin-airfoil theory does not see. XFOIL's full-inverse and mixed-inverse routines design the surface speed of a thick section (`@drela-xfoil`, "Inverse Formulation").
- **Viscosity.** No boundary layer, so no drag, no transition, no separation, no cl_max. A load that keeps the suction peak low is a good proxy for a gentle pressure recovery, but only a boundary-layer calculation (`P10.06`) tells you whether the flow stays attached.
- **Large camber and angle.** The linearisation fails as camber exceeds a few per cent or near the logarithmic ends of the a = 1.0 line.

### 2.6 The decision matrix

Let r_aj be the raw value of criterion j for candidate a. **Min–max normalisation** maps each criterion to [0, 1] over the candidate set, with 1 the best: s_aj = (r_aj − min_a r)/(max − min), or one minus that if smaller is better. With weights w_j summing to 1, the total is T_a = Σ_j w_j s_aj, and the highest total wins. Two consequences follow from linearity. If the weight of criterion k varies from 0 to 1 with the others rescaled in proportion, every total is a straight line, T_a(w) = w s_ak + (1 − w)R_a, so the winner can change only where two lines cross: the **flip points** are exact, not searched for. And the normalisation depends on the candidate set: adding or removing a weak candidate changes the scale of a criterion and can reorder the others. Neither fact is a flaw to hide; both belong in the report.

The lesson's four criteria are computed by its own models: thickness t/c (larger is better, structure); |cl_design − cl_ideal| (smaller is better: near its ideal cl a section has no leading-edge suction spike, a proxy for low drag); |cm_c/4| from thin-airfoil theory (smaller is better, trim and torsion); and the panel method's suction peak −cp,min at the mission's high-lift cl (smaller is better, a proxy for a gentle adverse gradient and a benign stall). The weights are judgments: trainer [0.30, 0.15, 0.20, 0.35] (stall manners and structure first), sailplane [0.15, 0.40, 0.20, 0.25] (operating near the ideal cl first). NACA TR-824's "Selection of Root Section" discusses the same trade-offs with measured data (`@ntrs-naca-tr-824-summary-of-airfoil-data`).

## 3. Build the computational model

Open `starter/main.py`. Given: `naca4`, `naca_thickness`, `naca4_camber`, `naca4_slope` (from `P10.02`), `thin_airfoil` (from `P10.04`, midpoint rule in θ), `panel_solve` (a Hess–Smith source–vortex panel method standing in for your `P10.05` solver; it handles several elements, which `P10.08` needs), `alpha_for_cl`, `mission_cls`, `build_candidates`, `select_sections`, `main`.

### Step 1: the load and its series

`glauert_load(A0, An, x)` evaluates Δcp from (1). `load_coefficients(dcp, n_terms, m)` is (2) by the midpoint rule on m points in θ; the midpoint rule is spectrally accurate for smooth periodic integrands and never evaluates the ends, where a target may be awkward. `load_cl_cm(dcp)` integrates cl and cm_c/4 of a load directly, in θ.

### Step 2: the camber line

`ideal_angle(An)` is the closure sum; `camber_from_coefficients(An, x)` returns z, dz/dx and α_i from (3) and (4), vectorised over x and n. `inverse_design(dcp)` chains them and returns a callable slope for the forward check. `naca_a1_camber(x, cl_i)` is (5); `camber_from_cl_cm(cl, cm)` is section 2.4.

```python
n = np.arange(1, len(An) + 1)
dzdx = ai + np.cos(np.multiply.outer(th, n)) @ An
```

### Step 3: thickness and the panel check

`thicken(camber, slope, t, n)` lays the NACA 00tt half-thickness perpendicular to the camber line (as the NACA sections do), with cosine spacing and nodes clockwise from the trailing edge; with the 2412 mean line it must reproduce `naca4("2412")` exactly (a test). The a = 1.0 slope is infinite at the ends, so evaluate it at clipped stations. `panel_load(X, Y, alpha, x)` interpolates cp on the two surfaces and returns cp,lower − cp,upper.

### Step 4: selection

`section_metrics` computes the four criteria; `decision_matrix` normalises and totals; `weight_monte_carlo` draws weight vectors from a Dirichlet distribution centred on the nominal weights (concentration 20, seed 7) and counts wins; `weight_flip_points` finds the exact crossings of section 2.6.

Running `main()` prints (abridged):

```
uniform load Delta_cp = 0.4: alpha_i = -1.07e-16 deg, max camber 0.02206 c at x = 0.500 c (closed form 0.02206 c at 0.5 c)
  N terms   max |z - z_a1| / c
     25     7.048e-05
     50     1.909e-05
    100     4.773e-06
    200     1.193e-06
  observed orders [1.8846 1.9995 2.0006]
forward thin-airfoil check of the recovered camber line at alpha_i: cl = 0.4000 (target 0.4000), cm_c/4 = -0.1000 (target -0.1000), A0 = 2.7e-17
  forward load at x = 0.05..0.95: min 0.3951, max 0.4050 (target 0.4)
roof-top load (a = 0.6, cl = 0.500): alpha_i = 1.292 deg, cm_c/4 = -0.0792, max camber 0.0369 c at x = 0.476
  forward check: cl = 0.5000, cm_c/4 = -0.0792, max |load error| on 0.02-0.98 = 0.0066
a = 1.0, cl_i = 0.4 with NACA 00tt thickness, panel method at alpha_i (thin theory cl = 0.4):
   t/c     cl_panel   mean |Delta_cp - 0.4| on 0.1-0.9
   0.02    0.3979     0.0069
   0.06    0.4119     0.0256
   0.12    0.4327     0.0530
   0.18    0.4537     0.0785
```

and writes `airfoil_design.png`: the recovered camber lines with the a = 1.0 closed form and the 2412 mean line; the panel load of the thickened design at 2 % and 12 % thickness against the design load; and the decision-matrix totals for both missions.

### Verify

`python -m pytest curriculum/phase-10/lesson-07-airfoil-design-selection` runs 18 tests: the load of random coefficients integrates to the thin-airfoil cl and cm; `load_coefficients` recovers the coefficients of a load built from them to 10⁻⁸; a symmetric load needs no ideal angle and the camber line closes; the slope is the derivative of the camber; the a = 1.0 closed form is symmetric, has the right maximum and slope; the uniform-load recovery converges at order 1.8–2.2 with error below 2 × 10⁻⁶ at N = 200; the forward analysis reproduces the uniform and roof-top loads; doubling the load doubles the camber and α_i; the (cl, cm) design hits both targets; `thicken` reproduces the 2412; the panel cl approaches 0.4 as t → 0 and the load departs monotonically as t grows; the panel load integrates to the panel cl; NACA 0012 has α_L0 = cm0 = cl_ideal = 0 and α_L0 and cm0 are linear in camber; the decision-matrix properties; the Monte Carlo on a symmetric problem; the flip points against a brute-force scan.

### Compare

The independent route is the **forward** analysis: the given `thin_airfoil` computes Aₙ from dz/dx by cosine projection, a different integral from the sine projection that produced them. For the uniform load it returns cl = 0.4000, cm_c/4 = −0.1000 and A0 = 2.7 × 10⁻¹⁷; the reconstructed load lies between 0.3951 and 0.4050 on 5–95 % chord (the residual is the Gibbs-like ripple of a 60-term series of a load that jumps to zero off the ends). For the roof-top it reproduces the load to 0.0066. The second independent route is the panel method, which knows nothing of Fourier series: at 2 % thickness its cl is 0.3979, 0.5 % from the design.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Truncation convergence.** Recover the a = 1.0 line with N = 25–200 terms. Predict the order from section 2.3. *Finding:* errors 7.05 × 10⁻⁵, 1.91 × 10⁻⁵, 4.77 × 10⁻⁶, 1.19 × 10⁻⁶ c; observed orders 1.88, 2.00, 2.00.
2. **Three design loads.** (a) Uniform cl_i = 0.4: α_i = 0, maximum camber 2.206 % at 50 %. (b) Roof-top: uniform to 60 % chord, then falling linearly to zero, cl = 0.5: α_i = 1.292° (unloading the rear needs a positive angle), cm_c/4 = −0.0792, maximum camber 3.69 % at 47.6 %. (c) cl = 0.30, cm_c/4 = −0.053: the cubic line of section 2.4. Predict before running which load needs the most camber per unit cl and why. *Finding:* the roof-top (7.4 % camber per unit cl against 5.5 % for the uniform load): moving load forward needs more curvature ahead of mid-chord and a positive ideal angle.
3. **Thickness.** Wrap NACA 00tt thickness around the cl_i = 0.4 line, t = 2–18 %, and run the panel method at α_i. Predict the panel cl at 12 %. *Finding:* cl = 0.3979, 0.4119, 0.4327, 0.4537; the mean load error grows almost linearly with thickness (0.0069, 0.0256, 0.0530, 0.0785). Thickness raises cl by about 0.7 × t/c relative, and at 12 % the load is no longer uniform: it rises to 0.46–0.49 over 10–30 % chord and falls to 0.37 at 90 %.
4. **Selection for two missions.** Compute the criteria for NACA 0012, 2412, 4412, 4415, 6409 and the a = 1.0 lines with cl_i = 0.4 and 0.8 (each thickened to 12 %).

| section | t/c | α_L0 (°) | cm_c/4 | cl_ideal | −cp,min at cl = 1.112 | trainer total | sailplane total |
|---|---|---|---|---|---|---|---|
| NACA 0012 | 0.120 | 0.000 | 0.0000 | 0.000 | 5.441 | 0.401 | 0.274 |
| NACA 2412 | 0.120 | −2.077 | −0.0531 | 0.256 | 3.285 | 0.623 | 0.504 |
| NACA 4412 | 0.120 | −4.154 | −0.1062 | 0.512 | 1.791 | 0.654 | 0.701 |
| NACA 4415 | 0.150 | −4.154 | −0.1062 | 0.512 | 1.631 | **0.818** | **0.792** |
| NACA 6409 | 0.090 | −6.232 | −0.1594 | 0.768 | 1.226 | 0.403 | 0.691 |
| a = 1.0, cl_i = 0.4 | 0.120 | −3.647 | −0.1000 | 0.400 | 2.673 | 0.629 | 0.584 |
| a = 1.0, cl_i = 0.8 | 0.120 | −7.294 | −0.2000 | 0.800 | 1.299 | 0.494 | 0.698 |

(The sailplane's suction peaks are at cl = 1.20: 6.224, 3.874, 2.171, 1.862, 1.360, 3.184, 1.496.) Predict the winners first. *Finding:* the NACA 4415 wins both.

5. **Weight sensitivity.** Monte Carlo over the weights and exact flip points. *Finding:* trainer: the 4415 wins 99.4 % of draws (the 2412 the other 0.6 %); the winner changes only if the cl-matching weight exceeds 0.493 (to the cl_i = 0.4 design), the moment weight exceeds 0.539 (to the 2412; above 0.564 the symmetric 0012), or the suction-peak weight exceeds 0.873; no thickness weight changes it. Sailplane: the 4415 wins 87.3 %, the 6409 11.8 %, the cl_i = 0.8 design 0.9 %; the winner changes at a cl-matching weight of 0.609 or a suction-peak weight of 0.622 (both to the 6409). The sailplane's choice is the fragile one.

Raw data format: CSV tables `N, err_z` (1); `load, alpha_i_deg, cm_c4, zmax, x_zmax` (2); `t_c, cl_panel, mean_load_error` (3); `section, criterion..., total, mission` (4); `criterion, flip_weight, new_winner, mission` (5); plus `airfoil_design.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| a = 1.0 camber, max | (cl_i/4π) ln 2 = 0.02206 c | inverse, N = 200: 0.02206 c | 1.2 × 10⁻⁶ c | series truncation, order 2 |
| uniform load forward | cl 0.4, cm −0.1, A0 = 0 | 0.4000, −0.1000, 2.7 × 10⁻¹⁷ | < 10⁻⁴ | quadrature of a log-singular slope |
| two-term design | cl 0.30, cm −0.053 | forward 0.3000, −0.0530 | < 10⁻⁶ | — |
| panel at t → 0 | cl 0.4 | 0.3979 at t = 0.02 | −0.5 % | thickness, panel discretisation |
| panel at t = 0.12 | cl 0.4 | 0.4327 | +8.2 % | thickness (not in the theory) |

The inverse method is verified to the accuracy its series predicts. Its *validity* is a different matter: section 5's last row is the model's error, not the code's, and it grows with thickness. A designer would therefore iterate: invert, thicken, analyse with the panel method (and `P10.06`), adjust the target, repeat. That loop is what XFOIL's mixed-inverse mode automates.

The selection answer: **for the trainer, on these inviscid criteria, the NACA 4415 (0.818) beats the 2412 (0.623) by a wide and robust margin; the 2412 would win only with a pitching-moment weight above 0.54.** That the trainer's actual 2412 loses is information, not an embarrassment: the matrix has no drag, no cl_max and no measured stall behaviour, and the 2412's real advantages (low cm, low drag, docile stall in the NACA data) live in exactly those missing criteria. What to improve: add cd at the cruise cl and a measured cl_max (`P10.06`, `P10.15`, TR-824's data) as criteria, and re-run the flip-point study.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(5) and of α_i = A2/3; `airfoil_design.png`; the five experiment tables; and a half-page selection note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) naming the winner, its margin, the weight at which it changes, and the criteria the matrix leaves out.

### Adversarial self-checks

- **What if the design cl doubles?** Every Aₙ, the camber and α_i double (the inversion is linear; a test checks it). Thin-airfoil theory does not care, but a section with twice the camber is further outside its small-camber assumption, and the panel check will show a larger departure.
- **What if the target load is infinite at the leading edge?** Then A0 ≠ 0 and the design is not at the ideal angle; the sine projection (2) silently absorbs part of the singularity into a slowly converging series. Check Δcp sin θ → 0 at θ → 0 before trusting the result.
- **Which assumption fails first?** Zero thickness: an 8 % lift error at 12 % thickness, far larger than any numerical error in the lesson.
- **What could falsify the selection?** A wind-tunnel polar showing the 4415 with a markedly higher cd at cl 0.33, or an abrupt stall: both are outside the criteria. Measured data beat proxies.
- **Does the result survive a change of method?** The inverse result survives a change of N (order 2) and a forward check by a different projection. The selection result survives 99.4 % of weight draws for the trainer but only 87.3 % for the sailplane.
- **Which decision changes if the uncertainty doubles?** With a Dirichlet concentration of 10 instead of 20 (roughly twice the variance in the weights), the trainer's 4415 still wins 95.6 % of draws, but the sailplane's falls to 76.8 % against 20.1 % for the 6409: the sailplane decision needs better-founded weights or a viscous criterion before it is made.

## Check yourself

1. In thin-airfoil theory, which series represents the load and which the camber slope, and what links them?
2. Why do we design for A0 = 0, and what is the angle called at which the design load occurs?
3. A target load is symmetric about mid-chord. What is its ideal angle, and why?
4. What camber line carries a uniform load, and what is its maximum camber for cl_i = 0.6?
5. Why does the error of the truncated a = 1.0 series fall as 1/N²?
6. Name two things the thin-airfoil inverse method cannot predict and the tool that can.
7. In a min–max-normalised decision matrix, why can removing the worst candidate change the ranking of the others?
8. Why are the flip points of a weighted-sum matrix exact, and what do they tell a reviewer that the totals do not?
