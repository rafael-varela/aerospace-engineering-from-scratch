# P10.15  Stall & Post-Stall

> A wing stalls where its sections first run out of lift, and the planform, not the section, decides where that is.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Kirchhoff/Helmholtz separation-point lift model cl = a₀α((1 + √f)/2)² and its limits, and state the flat-plate post-stall limit used for large angles. 2. Fit the model to authentic NACA 0012 lift curves by nonlinear least squares, report parameter uncertainties and residuals, and explain where the model's physics does not match the data. 3. Implement a nonlinear lifting line with a full-range section curve, solved by a relaxed iteration, and prove that it reduces to the linear lifting line below stall and converges in the number of Fourier terms. 4. Predict where along the span the reference trainer's wing stalls first and its maximum lift coefficient for rectangular, tapered and twisted planforms. 5. Explain why lifting line past stall has multiple solutions, what that has to do with hysteresis, and why those post-stall solutions must not be trusted quantitatively. |
| **Prerequisites** | `P10.04` (lift-curve slope 2π, thin-airfoil limit), `P10.10` (lifting line by Fourier collocation), `P10.12` (stall-onset location from local cl), `P02.18` (regression and least squares), `P03.10` (nonlinear systems and fixed-point iteration), `P00.14`. Skill check: you can explain why an iteration x ← g(x) converges when |g′| < 1. |
| **Engineering question** | Where along the span does the reference trainer's wing (span 10 m, area 12 m²) stall first, and at what wing lift coefficient, if it is rectangular, tapered to 0.45, or rectangular with 4° of washout, all with the NACA 0012 section fitted to Ladson's Re = 2 × 10⁶ data? The answers: rectangular, at the root, at α = 16.7° and C_L = 1.31 (0.90 of the section's fitted cl_max of 1.456); tapered, at 56 % of the semi-span, at 17.1° and C_L = 1.38; with washout, at the root, at 17.2° and C_L = 1.22. |
| **Assumptions and limits** | Section data are 2-D, steady and taken at one Reynolds number (2 × 10⁶; the trainer's approach is 2.46 × 10⁶ and its tips are not at a different Re in this model); the Kirchhoff model describes trailing-edge separation, while the NACA 0012 at this Re stalls abruptly from the leading edge, so the fitted sigmoid is nearly a step; the flat-plate limit c_n = 2 sin α is an empirical approximation; lifting line assumes a straight, unswept, high-aspect-ratio wing and a planar wake; beyond the first local stall the nonlinear lifting line is ill-posed (its solutions are not unique and depend on the discretisation). |
| **Analytical work** | The Kirchhoff lift and its limits f = 1 and f = 0; the flat-plate normal-force model and its cl and cd; the Gauss–Newton (least-squares) fit and its parameter covariance; the nonlinear lifting-line equation; the "decambering" iteration and why it is stable; the stability limit of the direct Fourier fixed point. |
| **Computational deliverable** | `starter/main.py`: `kirchhoff_cl`, `separation_point`, `fit_kirchhoff`, `flat_plate`, `full_range_cl`, `lifting_line_linear`, `lifting_line_nonlinear`, `fixed_point_relaxation_limit`, `lifting_line_fixed_point`, `stall_onset`, `section_stall_angle`, `stall_onset_alpha`, `wing_cl_max` (`find_data_file`, `read_ladson_block`, `ladson_reynolds_numbers`, the three planforms, `_ll_grid`, `_ll_out` and `main` are given). Tests prove: Kirchhoff reduces to 2πα at f = 1 and to a quarter at f = 0; the sigmoid's properties; the fit recovers known parameters from model data and fits the authentic data to rms < 0.05; the flat plate's maximum at 45° and cd = 2 at 90°; the blend; the linear lifting line equals the elliptic closed form; the nonlinear lifting line equals the linear one for a straight section and below stall; the iteration converges and agrees with the textbook fixed point, which fails beyond its computed stability limit; the Fourier series converges; rectangular wings stall at the root and tapered ones outboard. |
| **Practical activity** | **computational experiment** fitted to authentic data: the Kirchhoff model fitted to the six free-transition NACA 0012 lift curves of Ladson (NASA TM-4074); the stall pattern of three planforms with the nonlinear lifting line; iteration and Fourier convergence studies; an up-and-down α sweep that exposes the model's non-unique post-stall solutions. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `stall.png`, the fit table with uncertainties and residuals, the stall-pattern table for the three planforms, the convergence studies, and a half-page note recommending a planform change for benign stall, with the model's limits stated. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the Kirchhoff limits and the iteration are derived; all public tests pass; the fit is reported with residuals and a physical comment on their pattern; the stall locations are explained; the post-stall non-uniqueness is reported as a limitation, not a result. |
| **References** | `@ntrs-naca-tn-2502-airfoil-stall-types` (McCullough & Gault, NACA TN-2502, https://ntrs.nasa.gov/citations/19930083422): Summary; Introduction; Concluding remarks (required-free: trailing-edge, leading-edge and thin-airfoil stall). `@leishman-intro-aerospace-flight-vehicles` "Stalling & Spinning" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/maximum-lift-stalling-spinning/): General Flow Physics of Stall; Two-Dimensional Stall Characteristics; Three-Dimensional Stall; Stall Control Devices; and "Lifting Line & Finite Wing Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lifting-line-theory/): Numerical Solution to the Monoplane Wing Equation (required-free, link only). `@ntrs-nasa-tm-4074-ladson-naca0012` (https://ntrs.nasa.gov/citations/19880019495; transcribed in `curriculum/phase-10/data/`): Tests and Procedures; Presentation of Results (required-free: the data fitted here). `@ntrs-naca-tr-824-summary-of-airfoil-data` Table II "Maximum Lift and Stalling Characteristics of Models Tested in the NACA 19-Foot Pressure Tunnel"; "Lift Characteristics of Smooth Airfoils" (supplementary-free). `@mit-ocw-16-01-unified-engineering-fluids` 16.01 Fluids (spring) Lecture 9 "General Wings" (supplementary-free: lifting line for general planforms). `@ntrs-naca-tn-1422-sivells-three-wings` (supplementary-free: measured and calculated characteristics of tapered wings with and without washout, for Extension 2). `@mit-ocw-16-13-aerodynamics-viscous-fluids` Lecture 17 "... Form Drag, Stall Mechanisms." (supplementary-free). `@scipy-userguide` "Optimization (scipy.optimize)" > "Least-squares minimization (least_squares)" (supplementary-free). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics* (the nonlinear lifting-line iteration as a worked numerical method). |
| **Extensions** | ★ Refit with the Re = 8.86 × 10⁶ block, whose stall is gentler, and see how the stall pattern and CL_max change. ★ Compare the lifting-line CL_max of a taper-0.4 wing with and without 2° washout against the trends reported in NACA TN-1422 (different sections, so compare trends, not values). ★ Add a spanwise Reynolds-number variation (tip chord smaller) by interpolating between Ladson blocks, and see whether the tapered wing's stall moves outboard. |
| **Time** | unmeasured |

## 1. The engineering problem

How an aeroplane stalls matters as much as when. A wing whose root stalls first keeps its ailerons in attached flow and gives the pilot buffet from the separated wake on the tail as a warning; a wing whose tips stall first can drop a wing abruptly and enter a spin. Trainers are designed for the first behaviour. The designer controls it through the planform (taper puts relatively more lift outboard), twist (washout lowers the tip's angle of attack), and the section. To predict it you need two things: a section lift curve that goes through and beyond its maximum, and a wing model that turns it into a spanwise lift distribution.

This lesson builds both. The section curve is the Kirchhoff/Helmholtz separation-point model, a physical formula with three parameters, fitted to the authentic NACA 0012 data of Ladson at Re = 2 × 10⁶ and extended to 90° with a flat-plate limit. The wing model is Prandtl's lifting line made nonlinear: the straight-line section of `P10.10` is replaced by the fitted curve, and the equations are solved iteratively. Applied to the reference trainer's wing (span 10.0 m, area 12.0 m², the section varied to the NACA 0012 for this study, as the shared conventions allow), it predicts where each of three planforms stalls first.

Success criteria, written first: the fit reports its parameters with standard errors and an rms residual below 0.05 in cl; the nonlinear lifting line reproduces the linear solution to 10⁻⁹ for a straight section and to 10⁻⁵ below stall; its iteration converges to 10⁻¹⁰ and agrees with the textbook fixed point to 10⁻⁸ where both converge; the Fourier series converges (changes falling by more than 4 per doubling); every post-stall claim is accompanied by its dependence on the discretisation.

## 2. Understand and derive

### 2.1 The Kirchhoff/Helmholtz separation-point model

Kirchhoff and Helmholtz modelled a flat plate whose flow separates at a point along its chord, with a free streamline behind it. Let f be the fraction of the chord over which the flow is attached (f = 1: attached everywhere; f = 0: separated from the leading edge). The theory gives the normal force, and for small angles the lift,

**cl = a₀ α ((1 + √f)/2)².**   (1)

Limits: f = 1 gives **cl = a₀α**, thin-airfoil theory's 2πα when a₀ = 2π (`P10.04`); f = 0 gives **a₀α/4**, the fully separated plate. Between them, the separation point moving forward lowers the lift progressively: that is trailing-edge stall. The model needs f(α). Here it is a sigmoid,

**f(α) = 1/(1 + exp((|α| − α₁)/S)),**   (2)

equal to ½ at α₁, with S setting how abruptly the separation moves forward. Three parameters, (a₀, α₁, S), are fitted.

### 2.2 The flat-plate post-stall limit

Far beyond stall any section behaves like a flat plate whose force is nearly perpendicular to it. A stated empirical model: the normal-force coefficient c_n = 2 sin α (close to the c_n ≈ 2 of a long flat plate broadside to the flow), with no chordwise force, so

**cl = c_n cos α = sin 2α,  cd = c_n sin α = 2 sin²α.**   (3)

cl peaks at 1 at 45°, and cd reaches 2 at 90°. (Equation (1) is only valid for moderate α: with f → 0 it keeps growing linearly, 2.26 at 90°.) The full-range curve blends (1) into (3) with a smoothstep weight between 20° and 30°: exactly the fitted Kirchhoff curve below 20°, exactly the plate above 30°.

### 2.3 Fitting by least squares (`P02.18`)

Given the measured pairs (α_k, cl_k), choose p = (a₀, α₁, S) to minimise Σ r_k², r_k = cl_model(α_k; p) − cl_k. The model is nonlinear in α₁ and S, so the minimum is found iteratively (`scipy.optimize.least_squares`, a trust-region Gauss–Newton method) from a start at a₀ = 2π and α₁ just above the angle of the largest measured cl. At the solution, with J the Jacobian of r and σ² = Σr²/(n − 3), the parameter covariance is σ²(JᵀJ)⁻¹: the standard errors in the fit table come from it. They describe the fit's statistical precision, not whether the model is right; the residual pattern says that.

### 2.4 The nonlinear lifting line

Lifting line (`P10.10`) writes the circulation as Γ(θ) = 2bV Σ A_n sin nθ with y = (b/2) cos θ, which gives the induced angle α_i(θ) = Σ nA_n sin nθ/sin θ, and requires at every station that the section's lift match the circulation: ½ c cl = Γ/V. With a straight-line section this is the linear monoplane equation. With a curved section,

**2b Σ A_n sin nθ = ½ c(θ) · cl_sec(α + τ(θ) − α_i(θ)),**   (4)

τ the local twist, is nonlinear in A.

**The textbook iteration** relaxes the coefficients directly: compute α_i from A, cl_sec from the effective angle, new coefficients from (4), and A ← A + ω(A_new − A). It is explicit in α_i, and the high harmonics make it stiff: for a straight-line section of slope a₀ its iteration matrix J has real negative eigenvalues, and it converges only if **ω < 2/(1 − λ_min(J))**, which the code computes as **0.062** for the trainer's rectangular wing with 16 terms (a₀ = 5.74).

**The decambering iteration** used here keeps the linear lifting line as the solver and iterates only on a spanwise zero-lift-angle shift δ(y): solve the linear system with slope a₀ and zero-lift angle δ(y); read off α_eff; the section says cl_sec(α_eff); the shift that would put the straight line through that point is δ_new = α_eff − cl_sec(α_eff)/a₀; relax δ ← δ + ω(δ_new − δ). At convergence a₀(α_eff − δ) = cl_sec(α_eff) at every station, so (4) holds with the true section curve. The stiff part (the induced angle) is inside the linear solve, so ω close to 1 works below stall.

**Stall onset.** The critical-section estimate of the wing's maximum lift is the C_L at which the first station reaches the section's stall angle (the angle of its maximum cl); `stall_onset_alpha` finds that α by bisection.

### 2.5 Why the model breaks past stall

Beyond its maximum the section curve has a negative slope. Then (4) can have several solutions at the same α: one with the root stalled and the tips attached, one with more of the span stalled, and so on. Which one an iteration finds depends on where it starts. A real wing behaves similarly (stall cells, hysteresis between increasing and decreasing α), but lifting line has no physics to choose among them: no spanwise flow, no unsteady wake, a 2-D section curve that itself comes from steady 2-D measurements. Results past the first local stall are therefore qualitative at best.

### 2.6 Where the model stops

- The NACA 0012 at Re = 2 × 10⁶ stalls from the leading edge (TN-2502's leading-edge stall type): the lift collapses within a degree. Kirchhoff's model describes trailing-edge separation; fitting it to an abrupt stall gives S ≈ 0.28°, a step, and the fitted curve overshoots the last pre-stall points.
- 2-D data on a 3-D wing: near the root (fuselage junction) and the tips the flow is three-dimensional.
- Unsteadiness: dynamic stall (`P10.19`) changes cl_max when α changes quickly.

## 3. Build the computational model

Open `starter/main.py`. Given: the data readers, the three planforms (`rectangular_wing`, `tapered_wing(0.45)` with the same area and span, `twisted_wing(4.0)` with linear washout), the lifting-line grid and output helpers `_ll_grid` and `_ll_out`, and `main`.

### Step 1: the section

`kirchhoff_cl` (1), `separation_point` (2), `flat_plate` (3), `full_range_cl` (the blend), `section_stall_angle` (dense search for the maximum of the fitted curve). `fit_kirchhoff(alpha_deg, cl)` wraps `least_squares` with bounds a₀ ∈ [3, 8], α₁ ∈ [5°, 30°], S ∈ [0.05°, 10°] and returns the parameters, the residuals and the standard errors.

### Step 2: the lifting lines

`lifting_line_linear` is `P10.10`'s collocation with odd terms on the half span. `lifting_line_nonlinear` is the decambering iteration:

```python
for it in range(1, max_iter + 1):
    A = np.linalg.solve(M, mu * (alpha + tw - delta) * np.sin(th))
    alpha_i = np.sin(np.outer(th, n)) @ (n * A) / np.sin(th)
    a_eff = alpha + tw - alpha_i
    cl_sec = section_cl(a_eff)
    r = np.max(np.abs(a0 * (a_eff - delta) - cl_sec))
    if r < tol: break
    delta = delta + relax * ((a_eff - cl_sec / a0) - delta)
```

`lifting_line_fixed_point` is the textbook iteration, kept for comparison; `fixed_point_relaxation_limit` computes the eigenvalue bound of §2.4.

### Step 3: stall

`stall_onset` (the station with the largest cl/cl_max), `stall_onset_alpha` (bisection on α for the first station reaching the section stall angle) and `wing_cl_max` (an α sweep with continuation).

Running `main()` prints (abridged):

```
Kirchhoff fits to Ladson NACA 0012, free transition (NASA TM-4074):
  Re = 2e+06: a0 = 5.741 +- 0.051 /rad, alpha1 = 15.72 +- 0.05 deg, S = 0.281 deg, rms residual 0.0235, max |residual| 0.0620; model cl_max 1.456, measured max cl 1.408
  Re = 8.86e+06: a0 = 6.284 +- 0.018 /rad, alpha1 = 20.89 +- 0.09 deg, S = 3.056 deg, rms residual 0.0034, max |residual| 0.0078; model cl_max 1.700, measured max cl 1.693
linear section, alpha = 5 deg: linear LL CL = 0.395788, nonlinear iteration CL = 0.395788 (1 iterations, converged True), e = 0.9277
elliptic wing, 2 pi section: CL = 0.442187, closed form 2 pi alpha AR/(AR+2) = 0.442187, e = 1.000000
stall pattern (section: cl_max = 1.456 at 14.81 deg):
  rectangular                 : e = 0.9277, CL_alpha = 4.535/rad; first local stall at alpha = 16.67 deg, 2|y|/b = 0.000, CL = 1.313 (CL/cl_max = 0.902)
  tapered 0.45                : e = 0.9830, CL_alpha = 4.661/rad; first local stall at alpha = 17.10 deg, 2|y|/b = 0.556, CL = 1.384 (CL/cl_max = 0.951)
  rectangular, 4 deg washout  : e = 0.9440, CL_alpha = 4.535/rad; first local stall at alpha = 17.24 deg, 2|y|/b = 0.000, CL = 1.216 (CL/cl_max = 0.836)
```

and writes `stall.png`: the six data sets with their fits; the full-range section curve with the flat plate and the data; and local cl/cl_max along the span at the first local stall for the three wings.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-15-stall-post-stall`. The tests: (1) at f = 1 and f = 0 and its monotonicity in f; (2) at α₁, its symmetry and limits; the fit recovering (5.9, 16°, 1.2°) from noiseless model data with rms < 10⁻⁸; the fit to Re = 2 × 10⁶ (14 points, rms < 0.05, a₀ within 0.85–1.02 × 2π, cl_max within 0.1 of the data); (3) at 45° and 90° and its resultant 2 sin α; the blend's two ends; the elliptic closed form CL = 2πα AR/(AR + 2), e = 1 and C_Di = C_L²/(πAR); nonlinear = linear for a straight section on a twisted wing and below stall, with the lifting-line equation satisfied to 10⁻⁸; convergence near stall from a deliberately wrong slope, agreement with the fixed point to 10⁻⁸, failure of the fixed point at 2.5 times its limit; Fourier convergence; root stall for the rectangle, outboard for the taper, root with washout, and C_L below the section cl_max.

### Compare

Linear and nonlinear lifting line agree to 10⁻⁹ for a straight section (one iteration), and the linear solver reproduces the elliptic closed form CL = 0.442187 at 5° exactly. Below stall, at 12°, the decambering iteration and the textbook fixed point give the same C_L = 0.949890; the first needs 2–13 iterations (ω = 1.0–0.3), the second 376 at ω = 0.05 and fails at ω = 0.1. The CL change per doubling of Fourier terms at 12° is 9.1 × 10⁻⁴, 5.2 × 10⁻⁵, 3.1 × 10⁻⁶, 1.9 × 10⁻⁷ (a factor of about 17 per doubling).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment fitted to authentic data** (Ladson, NASA TM-4074, NTRS 19880019495, US Government work, public domain; NACA 0012 in the Langley Low-Turbulence Pressure Tunnel at M = 0.15, cl from surface-pressure integration, wall-corrected by the report; free transition, Table I; see `curriculum/phase-10/data/README.md`). Write each prediction before running it.

1. **The fits.** All six free-transition blocks. Predict how S changes with Re. *Finding:*

| Re | a₀ (/rad) | α₁ (°) | S (°) | rms | max |r| | model cl_max | measured max cl |
|---|---|---|---|---|---|---|---|
| 2.00 × 10⁶ | 5.741 ± 0.051 | 15.72 ± 0.05 | 0.281 | 0.0235 | 0.062 | 1.456 | 1.408 |
| 3.94 × 10⁶ | 5.987 ± 0.035 | 17.07 ± 0.08 | 1.044 | 0.0122 | 0.022 | 1.472 | 1.484 |
| 5.97 × 10⁶ | 5.867 ± 0.057 | 18.14 ± 0.06 | 0.386 | 0.0352 | 0.079 | 1.697 | 1.660 |
| 8.86 × 10⁶ | 6.284 ± 0.018 | 20.89 ± 0.09 | 3.056 | 0.0034 | 0.008 | 1.700 | 1.693 |
| 1.19 × 10⁷ | 6.144 ± 0.146 | 18.14 ± 0.20 | 0.953 | 0.0635 | 0.177 | 1.636 | 1.617 |
| 1.89 × 10⁷ | 6.214 ± 0.165 | 18.11 ± 0.30 | 1.146 | 0.0765 | 0.231 | 1.614 | 1.659 |

S scatters rather than trending: where the data contain one post-stall point after an abrupt drop (2.0, 5.97 × 10⁶), the fit is a near-step (S < 0.4°) and overshoots the largest measured cl by 0.037–0.048; at 8.86 × 10⁶ the data end at 18.4° with the lift curve still bending over gently, and the fit is good to 0.008. At the two highest Re the lift falls by 0.40 and 0.57 between two angles 0.2° and 0.85° apart, and the largest residuals (0.18 and 0.23) are at that pair of points: the drop is too sudden for even a narrow sigmoid to place. The fitted a₀ rises from 5.74 to 6.2–6.3 /rad with Re (91–100 % of 2π), the expected thinning of the boundary layer.

2. **The flat-plate limit.** cl at 20, 30, 45, 60 and 90°. *Finding:* Kirchhoff (f → 0) 0.502, 0.752, 1.127, 1.503, 2.255; plate 0.643, 0.866, 1.000, 0.866, 0.000 (cd 0.234, 0.500, 1.000, 1.500, 2.000); the blended curve 0.502, 0.866, 1.000, 0.866, 0.000. Equation (1) alone would be absurd beyond 30°.

3. **Iteration study.** At 12°, decambering with ω = 0.3–1.0 and the fixed point with ω = 0.02–0.2. *Finding:* above; the fixed point's computed limit is 0.062, it converges at 0.05 and fails at 0.1.

4. **Stall pattern.** Three planforms, first local stall. Predict the location for each. *Finding:* rectangular: root (2|y|/b = 0.000), 16.67°, C_L = 1.313; tapered 0.45: 2|y|/b = 0.556, 17.10°, C_L = 1.384; 4° washout: root, 17.24°, C_L = 1.216. The rectangle's local cl is largest at the root because its outer sections are unloaded by the tip vortices; taper moves load outboard and the critical station with it (and raises e to 0.983); washout reduces the tips' angles further, keeping stall at the root at the price of 7 % lower C_L at stall onset. The largest C_L found in the continuation sweep (1.302, 1.378, 1.233) is close to the critical-section value: with an abrupt section stall the wing loses lift as soon as one station stalls.

5. **Up-and-down sweep.** α from 12° to 22° and back, each solve starting from the previous one. *Finding:* going up, C_L = 0.950, 1.108, 1.266 at 12, 14, 16°; coming down, 0.645, 0.696, 0.810: the iteration stays on a branch with the root stalled down to 12°, 4.6 tenths of C_L lower. A second check shows this is a property of the model's ill-posedness, not of the wing: with 8 Fourier terms the down-sweep reattaches by 12° (C_L 0.950), with 32 it does not (0.661), and the number of stalled stations depends on the resolution. Starting from zero at 16°, the two iterations converge to different solutions (C_L = 1.266 and 0.720). Real wings do show stall hysteresis, but this model cannot predict its size.

Raw data format: CSV tables `Re, a0, a0_se, alpha1_deg, alpha1_se_deg, S_deg, rms, max_resid, clmax_model, clmax_data` (1); `alpha_deg, cl_kirchhoff, cl_plate, cd_plate, cl_blend` (2); `method, relax, converged, iterations, CL` (3); `planform, e, CL_alpha, alpha_onset_deg, eta_onset, CL_onset, CL_max_sweep` (4); `alpha_deg, CL_up, CL_down, n_terms` (5).

## 5. Validate and improve

| quantity | reference | model | difference | cause |
|---|---|---|---|---|
| cl, Re = 2 × 10⁶ | Ladson data (14 points) | Kirchhoff fit | rms 0.024, max 0.062 | abrupt leading-edge stall |
| cl_max, Re = 2 × 10⁶ | measured 1.408 | fit 1.456 | +3.4 % | step-like sigmoid overshoot |
| linear LL, elliptic wing | 2παAR/(AR + 2) = 0.442187 | 0.442187 | < 10⁻⁹ | — |
| nonlinear vs linear, straight section | 0.395788 | 0.395788 | < 10⁻⁹ | — |
| decambering vs fixed point, 12° | 0.949890 | 0.949890 | < 10⁻⁸ | — |

The engineering answer: **the rectangular trainer wing stalls first at the root, at α ≈ 16.7° and C_L ≈ 1.31; a 0.45 taper raises C_L at stall onset to 1.38 but moves the first stall to 56 % of the semi-span, towards the ailerons; 4° of washout keeps the rectangle's root stall with a lower C_L of 1.22.** For a trainer the rectangular wing's root-first stall is the desirable pattern; if a taper is wanted for its higher e, it needs washout or a root-stall strip to restore root-first stall (Problem 4 tests how much). The numbers rest on a 2-D section at one Re and on lifting line's critical-section reasoning: the validation route is a wind-tunnel test of the wing with tufts, the kind of finite-wing data NACA TN-1422 reports for tapered wings with and without washout (Extension 2), and the comparisons of `P10.22`.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(4) and the two iterations; `stall.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on the trainer's stall pattern and a planform recommendation, with the model's limits.

### Adversarial self-checks

- **What if the section stalled gently?** With the Re = 8.86 × 10⁶ fit (S = 3.1°) the wing's lift would round over instead of dropping, and the continuation sweep would carry C_L beyond the critical-section value (Extension 1).
- **What if AR → ∞?** The induced angle vanishes; every station sees α + τ, the rectangle stalls along the whole span at once, and C_L,max → cl_max.
- **Which assumption fails first in the real world?** The 2-D section data at the root (fuselage junction) and the abrupt-stall section near its maximum, where small differences in Re or surface finish move cl_max by 0.1.
- **What measurement could falsify the model?** Tufts or a surface-flow pattern on a wind-tunnel wing at increasing α: the first tufts to reverse should be at the root for the rectangle and near mid-semi-span for the taper.
- **Does the result survive a change of method or grid?** Up to the first local stall, yes: the two iterations agree to 10⁻⁸ and the Fourier truncation error is 5 × 10⁻⁵ at 16 terms. Past it, no (Experiment 5).
- **Which decision changes if the uncertainty doubles?** The taper decision: if the predicted onset station of the tapered wing were uncertain by ±0.2 of the semi-span instead of ±0.1, it might be over the ailerons, and the design would need washout regardless.

## Check yourself

1. Write the Kirchhoff lift formula and its values at f = 1 and f = 0.
2. What does the parameter S in the sigmoid f(α) control, and why did the fit at Re = 2 × 10⁶ give a very small S?
3. What are cl and cd of the flat-plate model at 45° and 90°?
4. Why does the textbook fixed-point iteration of the nonlinear lifting line need a small relaxation factor, and how does the decambering iteration avoid it?
5. Why does a rectangular wing stall first at the root?
6. What does washout do to the stall pattern and to C_L at stall onset?
7. Why are nonlinear lifting-line solutions past the first local stall not unique, and what does that mean for predicting hysteresis?
