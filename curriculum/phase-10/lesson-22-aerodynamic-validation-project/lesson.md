# P10.22  Aerodynamic Validation Project

> Agreement with a wind tunnel is not a yes or a no. It is an error, an uncertainty on each side, and a statement of the conditions under which the model can be trusted.

| | |
|---|---|
| **Learning outcomes** | 1. Organise a verification-and-validation (V&V) study of an aerodynamic code into code verification, solution verification and validation, and say what evidence belongs to each. 2. Measure the observed order of a vortex lattice's convergence, extrapolate by Richardson, and estimate the lattice-size uncertainty with Roache's grid convergence index (GCI). 3. Derive and apply a section-lift-slope correction to a flat-plate lattice from lifting-line theory, with the section slope taken from published wind-tunnel data. 4. Compare predictions with twelve measured planforms, combining reading uncertainty, run-to-run scatter, numerical and input uncertainty into a validation uncertainty, and report a validation metric. 5. Identify where the linear model fails (low aspect ratio, taper) and write a validation conclusion with its domain of validity. |
| **Prerequisites** | `P10.11` (the vortex lattice method), `P10.10` (lifting line), `P10.12` (planform effects), `P10.13` and `P10.21` (Ladson's NACA 0012 data), `P03.18` (verification and validation), `P03.11` (observed order), `P00.14` (uncertainty, uniform distributions), `P02.18` (least squares with covariance). Skill check: you can combine independent standard uncertainties in quadrature and convert a ± half-width of a uniform band to a standard uncertainty. |
| **Engineering question** | How well does your vortex lattice, with a section-lift-slope correction from Ladson's NACA 0012 data, predict the measured lift-curve slopes of twelve low-speed wing planforms (NACA TN-1669, TN-1671, TN-1422), and where does it fail? The answers: mean absolute error 4.3 %, mean signed error +0.2 %, largest error +16.2 % (the AR 1.34 unswept wing, whose measured lift is nonlinear); 12 of 15 comparisons (80 %) agree within the k = 2 combined uncertainty; the other failures are the two tapered 45° wings of TN-1671 (+5.8 % and +7.8 %); without the section correction the mean signed error is +3.7 %. |
| **Assumptions and limits** | Flat vortex lattice (inviscid, incompressible, attached flow, flat wake), uniform spacing, half-wing with a symmetry image; the section correction is exact for lifting line and approximate for the lattice; NACA 0012 section slope from Ladson at Re = 2.0 × 10⁶ used for all wings (the tests ran at about 0.8–1.6 × 10⁶, and TN-1422's sections are 64-210 and 65-210); dihedral and washout ignored (neither changes a linear lift slope much); compressibility at M = 0.13–0.17 ignored (about 1 %). Data: TN-1669 and TN-1671 values digitized from figures; the CL range of the TN-1669 slopes is not stated; jet-boundary corrections applied by the reports, no blockage corrections. |
| **Analytical work** | The three parts of V&V; observed order, Richardson extrapolation and the GCI from three grids; the lifting-line chord-scaling identity CL_α(AR, a₀) = κ CL_α,2π(AR/κ); the validation comparison error E = S − D and the combined uncertainty u_val = √(u_num² + u_input² + u_D²); uniform-band standard uncertainties. |
| **Computational deliverable** | `starter/main.py`: `vlm_cla`, `observed_order_3`, `richardson`, `gci_fine`, `grid_study`, `section_slope_from_ladson`, `load_finite_wing_data`, `validation_uncertainty`, `validation_table`, `validation_metrics`, `symmetry_residual` (the data finder, the Biot–Savart and horseshoe helpers, `Lattice`, `trapezoidal_lattice`, `solve_vlm`, `lift_coefficient` and `main` are given, copied in compact form from `P10.11`). Tests (14) prove: a symmetric wing's span load is symmetric to 10⁻¹²; the half-wing with its image equals the full wing; lift is linear in sin α; a very long unswept wing tends to 2π and to κ2π with the correction; the lattice converges monotonically with an observed order between 0.5 and 1.5; the three-grid formulas recover a manufactured solution exactly and the GCI covers its error; the grid study gives a small positive GCI; Ladson's section slope is physically sensible and rises with Re; the data file is read with its stated uncertainties; the correction lowers the slope less than proportionally; the uncertainty combination and metrics are right; and against the twelve planforms the mean error is below 8 %, at least 60 % agree, and the worst case is the AR 1.34 unswept wing. |
| **Practical activity** | **data-based investigation**: a V&V report comparing the learner's lattice with the authentic measured lift-curve slopes in `curriculum/phase-10/data/finite_wing_cla.csv` (NACA TN-1669, TN-1671, TN-1422; US Government works, public domain), with the section slope from Ladson's NACA 0012 data. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vlm_validation.png`, the V&V report (structure in part 6) with the code-verification checks, the grid study with GCI, the validation table for all 15 rows with uncertainties, the validation metric, and the conclusion with its domain of validity. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: each piece of evidence is assigned to the right part of V&V; all public tests pass; the GCI and the observed order are reported; the validation table carries every uncertainty term and the metric is stated with the assumptions it depends on; the failures are explained physically, not as "experimental error". |
| **References** | `@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep` NACA TN-1669: Introduction; Apparatus and Tests; Corrections; Results and Discussion; Conclusions (required-free). `@ntrs-naca-tn-1671-swept-wings-taper-ratio` NACA TN-1671 (required-free). `@ntrs-naca-tn-1422-sivells-three-wings` NACA TN-1422 (required-free): the measured slopes, transcribed in the shared data file with its README. `@ntrs-nasa-tm-4074-ladson-naca0012` Ladson, NASA TM-4074 (required-free): the section slope. `@ntrs-nasa-tn-d-6142-vortex-lattice-program` Margason & Lamar, NASA TN D-6142: "Effect of Vortex-Lattice Arrangement on Computed Aerodynamic Characteristics"; "Sample Cases" (supplementary-free): lattice arrangement and convergence. `@leishman-intro-aerospace-flight-vehicles` "Wind Tunnel Testing" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/wind-tunnel-testing/): Model Scaling & Similarity Parameters; Wall Interference; Data Reduction (required-free). `@mit-ocw-16-100-aerodynamics` Lecture 31 "Similarity in Wind Tunnel Testing" (required-free). `@nasa-nparc-cfd-vv-tutorial` (required-free): verification, validation and the GCI as used in `P03.18`. `@ntrs-naca-tr-921-span-loading` (supplementary-free): loading characteristics for comparison. Optional textbooks for deeper study (paid, never required): Katz & Plotkin, *Low-Speed Aerodynamics*, for lifting-surface theory and lattice convergence; library or used copies are fine. |
| **Extensions** | ★ Replace uniform spanwise spacing by cosine spacing (with the control points placed accordingly) and show that the observed order and the GCI improve. ★ Add a Polhamus-type leading- and side-edge vortex-lift term and see whether the AR 1.34 wing comes into agreement. ★ Digitize one more planform from one of the reports that the data README names as leads (TN-1862 lists 29 swept planforms and their sources) and add it to your own copy of the validation, with its reading uncertainty. |
| **Time** | unmeasured |

## 1. The engineering problem

A vortex lattice code (`P10.11`) gives lift-curve slopes, span loads, induced drag and stability derivatives in seconds, and designers use those numbers to size tails, set control-surface throws and predict handling. Before anyone builds on them, someone has to ask the question this lesson answers: how wrong are they, for which wings, and how do you know? "The code agrees with the wind tunnel" is not an answer. The code's own numerical error may be larger than the difference; the measurement has its own uncertainty, including the figure-reading error of digitized data and the scatter between two reports on the same wing; and the code needs an input, the section lift slope, that is itself uncertain.

The discipline that sorts this out is verification and validation (`P03.18`). **Code verification** shows the code solves its equations correctly. **Solution verification** estimates the numerical error of each answer. **Validation** compares the answers with measurements, with all the uncertainties on the table, and states where the model can be trusted. This lesson does all three for the learner's lattice against authentic data: twelve planforms (fifteen rows) measured at low speed in the Langley stability tunnel and the Langley 19-foot pressure tunnel, transcribed with their provenance in `curriculum/phase-10/data/finite_wing_cla.csv` (read its README first).

Success criteria, written first: the code-verification checks pass to round-off; the lattice's observed order and GCI are reported for every planform; each comparison carries u_num, u_input and u_D; the validation metric is stated with the assumptions it depends on; and each failure has a physical explanation.

## 2. Understand and derive

### 2.1 What each part of V&V can show

| part | question | evidence here | it cannot show |
|---|---|---|---|
| code verification | are the equations solved correctly? | symmetry, linearity, the 2-D limit, half-wing = full wing, observed order | that the equations describe the flow |
| solution verification | how large is this answer's numerical error? | three-lattice grid study, Richardson value, GCI | anything about modelling error |
| validation | does the model represent the flow for this purpose? | prediction against measurement with u_val | anything outside the tested conditions |

### 2.2 Observed order, Richardson and the GCI

Solve on three lattices refined by r = 2 in both directions, f₃ (coarse), f₂, f₁ (fine). If the error behaves as Chᵖ,

**p = ln((f₃ − f₂)/(f₂ − f₁)) / ln r,  f* ≈ f₁ + (f₁ − f₂)/(rᵖ − 1),**   (1)

provided the three converge monotonically. Roache's grid convergence index on the fine lattice,

**GCI = Fs |(f₂ − f₁)/f₁| / (rᵖ − 1),  Fs = 1.25,**   (2)

is meant as a band of roughly 95 % on the fine-lattice error; the lesson takes **u_num = GCI·|f₁|/2** as a standard uncertainty. When p is outside 0.5–2 it is clipped for the estimate, as common GCI practice recommends.

### 2.3 The section-slope correction

A flat lattice panel has the thin-airfoil slope 2π. A real section has a₀ = κ2π with κ < 1, because the boundary layer reduces the circulation (`P10.21` measured the gap: 12.7 % at 6 × 10⁶). In lifting-line theory (`P10.10`) the circulation obeys Γ(y) = ½Vc(y)a₀[α − α_i(y)], in which the chord appears only in the product c a₀. So a wing of chord c with slope κ2π has exactly the circulation of a wing of chord κc with slope 2π, at the same span. That wing has area κS and aspect ratio AR/κ, so

**CL_α(AR, a₀) = κ · CL_α,2π(AR/κ).**   (3)

The lesson applies (3) to the lattice: scale every chord by κ about the quarter-chord line (span and quarter-chord sweep unchanged), run the 2π lattice, multiply by κ. It is exact for lifting line and an approximation for the lattice, in which the chord also sets the chordwise vorticity distribution. (The DATCOM-type slope formula has the same scaling built in.)

κ comes from Ladson's NACA 0012 data (`@ntrs-nasa-tm-4074-ladson-naca0012`, Table I, free transition): a line through cl(α) for −4.5° ≤ α ≤ 8.5° gives a₀ = 5.9329 ± 0.0276 per radian at Re = 2.0 × 10⁶ and 6.1236 ± 0.0191 at 5.97 × 10⁶. The wings were tested at 0.8–1.6 × 10⁶, below Ladson's lowest Reynolds number; the lesson uses the 2.0 × 10⁶ value, **κ = 0.9443**, and takes the change from 2.0 to 5.97 × 10⁶ as the input uncertainty for extrapolating downward, **u_κ = 0.0303**. Its effect on each prediction is u_input = |∂CL_α/∂κ| u_κ, with the derivative by central difference.

### 2.4 The validation comparison

For each row, the comparison error is **E = S − D** (simulation minus data). As in the comparison-error approach of `P03.18`, the validation uncertainty combines the independent standard uncertainties of the simulation's numerics, its inputs and the data:

**u_val² = u_num² + u_input² + u_D²,  u_D² = (u_read/√3)² + (0.03 D/√3)²,**   (4)

where u_read is the stated figure-reading uncertainty (±0.0005 per degree for TN-1669, ±0.0003 for TN-1671, and ±0.0005, half the last printed digit, for TN-1422's tabulated values) and 3 % is the documented scatter between TN-1669 and TN-1671 on the same wing (0.0463 against 0.0449 per degree). Both are stated as half-widths, so each is divided by √3 (a uniform band, `P00.14`). A row is "within" if |E| ≤ 2u_val. The **validation metric** is the mean and largest |E/D|, the mean signed E/D (bias) and the fraction within. That fraction depends on the uncertainties assumed (Problem 4), so it is never reported alone.

### 2.5 Where the linear model should fail

The data README warns of two things. For the low-aspect-ratio wings (AR 1.34, and 2.61 unswept) TN-1669 reports that the lift-curve slope *increases* with CL: vorticity rolls up along the side edges into vortices above the wing (vortex lift), which a flat-wake linear model cannot produce, and which makes a measured "slope" depend on the CL range used (not stated in the report). And the two reports differ by 3 % on the same wing. High sweep (60°) is a third place to look: thick swept sections behave differently from the flat-panel model, and the jet-boundary corrections were based on unswept-wing theory.

## 3. Build the computational model

Open `starter/main.py`. The data finder (it walks up to the repository root, so the code runs from `curriculum/` and from `solutions/`), the Biot–Savart helpers, the trapezoidal planform builder (aspect ratio, taper, quarter-chord sweep, streamwise tips, optional chord scale), the half-wing solver with its symmetry image and the lift coefficient are given.

### Step 1: code verification

`symmetry_residual` solves the full wing without the image and compares Γ(y) with Γ(−y). `vlm_cla(AR, taper, sweep, nc, ns, kappa, half)` returns κ times the 2π slope of the chord-scaled planform (3), from a solution at 1°.

### Step 2: solution verification

`observed_order_3`, `richardson` and `gci_fine` are (1) and (2). `grid_study` runs the three lattices 2 × 6, 4 × 12 and 8 × 24 (chordwise × spanwise panels per half; 8 × 48 for the full span) and returns f, p, the Richardson value, the GCI and u_num.

### Step 3: inputs, data, validation

`section_slope_from_ladson` fits a₀ with `np.polyfit(..., cov=True)`; `load_finite_wing_data` reads the 15 rows and assigns the reading uncertainty from the source note; `validation_uncertainty` is (4); `validation_table` runs the grid study once per distinct planform, computes the sensitivity to κ, and fills E, E/D, u_num, u_input, u_val and the verdict; `validation_metrics` summarises.

### Verify

```
CODE VERIFICATION
  symmetry residual, AR 2.61 taper 0.25 sweep 45 (full wing): 7.89e-16
  half-wing with image vs full wing, AR 5.16 sweep 45: 3.3437765039 vs 3.3437765039
  2-D limit: AR 1000 unswept, CL_alpha = 6.2633 /rad (2 pi = 6.2832)
  linearity: CL(6 deg)/CL(2 deg) = 2.9951281005, sin ratio 2.9951281005
```

and the convergence of the AR 5.16 unswept wing with κ = 1:

| lattice (full span) | 1 × 6 | 2 × 12 | 4 × 24 | 8 × 48 | 16 × 96 |
|---|---|---|---|---|---|
| CL_α [/rad] | 4.364588 | 4.203807 | 4.107908 | 4.055608 | 4.028321 |

Observed orders 0.75, 0.88, 0.94: approaching first order from below. The span load falls to zero at the tips like a square root, which a uniform lattice resolves only to first order (cosine spacing would do better; Extension 1). The 14 tests check everything listed in the metadata, including that the three-grid formulas recover a manufactured solution f* + Chᵖ exactly.

### Compare: is the GCI honest?

For three planforms, the 8 × 48 value against a fourth, 16 × 96 lattice: actual differences 0.88 % (AR 1.34), 0.43 % (AR 9) and 0.68 % (AR 2.61, 60°), against GCIs of 2.5 %, 1.5 % and 1.6 %. The GCI covers the observed error with a margin, as intended.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** The data: fifteen measured lift-curve slopes of twelve planforms, from NACA TN-1669 (`@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep`; nine untapered wings, AR 1.34, 2.61, 5.16 at 0°, 45° and 60° sweep, NACA 0012 normal to the leading edge, Langley stability tunnel, about 1 × 10⁶, M = 0.13, digitized from its Figure 5 with ±0.0005 per degree), TN-1671 (`@ntrs-naca-tn-1671-swept-wings-taper-ratio`; 45° sweep, AR 2.61, taper 1.0, 0.5, 0.25, digitized from Figure 6 for CL < 0.6 with ±0.0003) and TN-1422 (`@ntrs-naca-tn-1422-sivells-three-wings`; AR 9, taper 0.4, three wings, printed in its Table I). All are US Government works in the public domain; the reports applied jet-boundary corrections and no blockage corrections. These values are the reports' measurements, not this curriculum's.

**Prediction first.** Before running the comparison, write down which rows you expect to fail and why (§2.5).

**The validation table** (fine lattice; CL_α per degree; uncertainties standard, per degree):

| report | AR | taper | Λ_c/4 | measured | predicted | E/D | p | GCI | u_num | u_input | u_val | within k = 2 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| TN-1669 | 1.34 | 1.00 | 0° | 0.0280 | 0.0325 | +16.2 % | 0.88 | 2.49 % | 0.00041 | 0.00025 | 0.00074 | no |
| TN-1669 | 2.61 | 1.00 | 0° | 0.0486 | 0.0504 | +3.7 % | 0.88 | 2.19 % | 0.00055 | 0.00073 | 0.00128 | yes |
| TN-1669 | 5.16 | 1.00 | 0° | 0.0689 | 0.0682 | −1.0 % | 0.87 | 1.92 % | 0.00065 | 0.00143 | 0.00199 | yes |
| TN-1669 | 1.34 | 1.00 | 45° | 0.0315 | 0.0310 | −1.6 % | 1.01 | 1.71 % | 0.00026 | 0.00034 | 0.00075 | yes |
| TN-1669 | 2.61 | 1.00 | 45° | 0.0463 | 0.0440 | −5.0 % | 1.00 | 1.71 % | 0.00038 | 0.00085 | 0.00126 | yes |
| TN-1669 | 5.16 | 1.00 | 45° | 0.0560 | 0.0544 | −2.9 % | 0.96 | 1.82 % | 0.00049 | 0.00133 | 0.00174 | yes |
| TN-1669 | 1.34 | 1.00 | 60° | 0.0286 | 0.0279 | −2.4 % | 1.13 | 1.24 % | 0.00017 | 0.00046 | 0.00075 | yes |
| TN-1669 | 2.61 | 1.00 | 60° | 0.0360 | 0.0357 | −0.8 % | 1.06 | 1.64 % | 0.00029 | 0.00083 | 0.00112 | yes |
| TN-1669 | 5.16 | 1.00 | 60° | 0.0434 | 0.0415 | −4.3 % | 0.95 | 1.85 % | 0.00038 | 0.00110 | 0.00142 | yes |
| TN-1671 | 2.61 | 0.25 | 45° | 0.0440 | 0.0474 | +7.8 % | 0.90 | 1.42 % | 0.00034 | 0.00087 | 0.00122 | no |
| TN-1671 | 2.61 | 0.50 | 45° | 0.0439 | 0.0465 | +5.8 % | 0.97 | 1.54 % | 0.00036 | 0.00087 | 0.00122 | no |
| TN-1671 | 2.61 | 1.00 | 45° | 0.0449 | 0.0440 | −2.1 % | 1.00 | 1.71 % | 0.00038 | 0.00085 | 0.00122 | yes |
| TN-1422 (64-210, 2° washout) | 9 | 0.40 | 0° | 0.086 | 0.0823 | −4.3 % | 0.71 | 1.51 % | 0.00062 | 0.00206 | 0.00263 | yes |
| TN-1422 (65-210, 2° washout) | 9 | 0.40 | 0° | 0.085 | 0.0823 | −3.1 % | 0.71 | 1.51 % | 0.00062 | 0.00206 | 0.00262 | yes |
| TN-1422 (65-210, no washout) | 9 | 0.40 | 0° | 0.085 | 0.0823 | −3.1 % | 0.71 | 1.51 % | 0.00062 | 0.00206 | 0.00262 | yes |

**The validation metric.** All 15 rows: mean |E/D| = **4.3 %**, largest **16.2 %** (wing 1, AR 1.34 unswept), mean signed **+0.2 %**, **80 %** within k = 2. Excluding AR 1.34 and 60° sweep (10 rows): mean 3.9 %, largest 7.8 %, 80 % within. Without the section correction (κ = 1): mean |E/D| 3.9 %, mean signed **+3.7 %**, 73 % within: the correction removes the bias without reducing the scatter, which is set by the failures.

**TN-1422.** The lattice gives 0.0823 per degree against the measured 0.085–0.086 and the report's own calculated 0.083–0.085. Using NACA 0012's a₀ for 6-series sections is the main suspect (u_input is the largest term in this row); the lattice's own error is in the same direction (the Richardson value, 0.0813, is lower still).

Raw data format: the table above as CSV with columns `report,planform,AR,taper,sweep_deg,measured_per_deg,u_read_per_deg,predicted_per_deg,E_over_D,p_obs,gci,u_num,u_input,u_val,within`; the grid-study table `planform,nc,ns_half,cla_per_rad`; `vlm_validation.png`.

## 5. Validate and improve

Every row is now explained or flagged:

- **AR 1.34 unswept (+16 %)**: the measured lift is nonlinear (vortex lift along the side edges grows faster than α), the slope's CL range is unstated, and the same effect makes the 45° wing of the same aspect ratio read higher than the unswept one in the data. The linear lattice is outside its domain here, as §2.5 predicted.
- **TN-1671 taper 0.25 and 0.5 (+8 % and +6 %)**: the lattice says taper raises the slope of these 45° wings; the report found little effect below CL = 0.6. Candidate causes: the approximate jet-boundary correction (based on unswept-wing theory), the narrow, highly loaded tips of the tapered planforms, and the 3 % scatter between the two reports. The data cannot separate them.
- **The rest** (11 rows): within 5 %, all within u_val, with a slight negative bias for the larger aspect ratios that is consistent with the lattice still converging downward and with the section-slope proxy.

The engineering answer: **the vortex lattice with a section-slope correction from Ladson's data predicts the lift-curve slope of these low-speed wings with a mean absolute error of 4.3 % and no significant bias (+0.2 %); 80 % of the comparisons agree within their k = 2 combined uncertainty. It is validated, at this level, for attached-flow lift slopes of wings of aspect ratio above about 2.5 with sweep up to 60°; it fails for the AR 1.34 wing (+16 %), where the measured lift is nonlinear, and overpredicts the effect of taper on 45° wings by 6–8 %.**

What to improve first: the section slope (the largest uncertainty term for most rows) by using sections' own data at the test Reynolds numbers; then cosine spacing to cut u_num. What would not help: finer uniform lattices for the AR 1.34 wing, whose error is physical.

## 6. Defend your work

Submit a V&V report (template in `projects/_template/report-template.md`) with these sections: (1) purpose and the quantity validated (CL_α, attached flow); (2) the model and its assumptions; (3) code verification (the four checks with their results); (4) solution verification (the grid study, observed orders, GCI for every planform); (5) the data and their uncertainty (sources, reading uncertainty, scatter, corrections applied by the reports); (6) the input uncertainty (a₀ and κ); (7) the validation table and metric; (8) the failures and their explanations; (9) the domain of validity and the conclusion. Include `vlm_validation.png`, the completed `main.py` and the passing tests.

### Adversarial self-checks

- **What if the aspect ratio goes to infinity?** CL_α → κ2π = 5.93 per radian, Ladson's section slope; the test checks the 2-D limit with and without κ.
- **What if the uncertainty I assumed for the data were halved?** The fraction within falls from 80 % to 73 % (1.5 % scatter); the metric's mean and largest errors do not change. That is why both are reported.
- **Which assumption fails first?** Linearity, for the low-aspect-ratio wings. Then the section-slope proxy (Re and section family).
- **What measurement could falsify the conclusion?** A new planform inside the claimed domain (say AR 4, 30° sweep, taper 0.5) whose measured slope differs from the prediction by more than 2u_val.
- **Does the result survive a change of lattice?** The 8 × 48 values move by 0.4–0.9 % on refinement to 16 × 96, inside their GCIs; no verdict changes.
- **Which decision changes if the uncertainty doubles?** If u_κ doubled, every row would pass and the comparison would show nothing: a validation is only as informative as its uncertainty is small. The tapered-wing disagreement would hide inside the band.

## Check yourself

1. What are code verification, solution verification and validation, in one sentence each?
2. How do you get the observed order, a Richardson value and a GCI from three lattices? What order does the uniform lattice show here?
3. How does the lesson correct a 2π lattice for a real section's lift slope?
4. Which terms make up the validation uncertainty in this lesson?
5. State the validation metric for the twelve planforms.
6. Which planform fails worst, and why can a linear lattice not fix it?
7. Why is "the fraction within uncertainty" never reported alone?
