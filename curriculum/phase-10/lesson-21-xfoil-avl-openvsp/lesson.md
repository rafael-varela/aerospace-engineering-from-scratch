# P10.21  Using XFOIL, AVL & OpenVSP

> A tool's polar without its settings is a rumour. Record the tool, the version and every setting; compare with numbers that mean something; and remember that the wind tunnel, not the tool, is the reference.

| | |
|---|---|
| **Learning outcomes** | 1. Write the exact XFOIL command sequence and AVL input file for a stated case (panels, Ncrit, Re, Mach, iteration limit; lattice density, reference quantities), and say what each setting changes. 2. Implement and test a polar reader that refuses any file whose header does not record tool, version and settings, and a parser that extracts the settings XFOIL writes into its own polar file. 3. Compute and interpret comparison metrics between any two polars: lift-curve slope and zero-lift angle from a fit over a stated linear range, the RMS cl difference at common angles, and the minimum drag. 4. Validate your own panel and boundary-layer codes against published wind-tunnel data (Ladson, NASA TM-4074) and attribute each difference to a missing physical effect or a setting. 5. Decide, when a tool, your code and the tunnel disagree, which to trust and what to write. |
| **Prerequisites** | `P10.02` (NACA generator), `P10.05` (vortex panel method), `P10.06` (Thwaites, Michel, turbulent momentum integral, Squire–Young), `P10.11` (vortex lattice method), `P10.13` (the Ladson data), `P03.18` (verification and validation), `P02.18` (least-squares fitting and its covariance). Skill check: you can fit a straight line and state the standard error of its slope. |
| **Engineering question** | How well do your own codes reproduce the published NACA 0012 data (Ladson, NASA TM-4074, M = 0.15), and is a polar handed to you by a tool trustworthy enough to use? The answers: inviscid panel lift is 12.7 % too steep at Re = 5.97 × 10⁶ (0.1205 against 0.1069 per degree; 16.4 % at 2 × 10⁶), the zero-lift angle agrees within 0.05°, and the one-way drag estimate is within 3.5 % of the measured minimum drag at 5.97 × 10⁶ (18.5 % high at 2 × 10⁶); a tool polar is usable only when its header records the tool, version and settings, and the reader written here refuses it otherwise. |
| **Assumptions and limits** | The own codes are two-dimensional, incompressible, attached-flow models: a linear-strength vortex panel method and a one-way (uncoupled) boundary layer with Thwaites' laminar method, Michel's transition criterion, a fixed-shape-factor turbulent momentum integral and the Squire–Young formula applied at x/c = 0.99. They do not feed the boundary layer back into the inviscid flow, so they overpredict the lift slope and say nothing about cl_max. The published data include the report's two-dimensional wall corrections. The XFOIL file layout assumed by the parser is the one in XFOIL's documentation; check it against your own file. |
| **Analytical work** | The meaning of each tool setting; the line fit for cl_α and α_L0 with the uncertainty of α_L0 from the fit covariance; why RMS differences must be computed at common angles and quoted with their range; the expected size of the viscous lift-slope reduction. |
| **Computational deliverable** | `starter/main.py`: `parse_header`, `read_tool_polar`, `parse_xfoil_polar`, `xfoil_polar_to_labelled_csv`, `polar_metrics`, `compare_polars`, `own_polar`, `validate_against_ladson` (`find_phase10_data`, `load_ladson`, `naca4`, `panel_system`, `vortex_panel`, the boundary-layer helpers, `bl_drag`, `write_tool_polar` and `main` are given). Tests (14) prove: the reader refuses a file with no header, with no version or settings, and an XFOIL or AVL polar missing a required setting, and round-trips a labelled one; the XFOIL-layout parser extracts version, Re, Mach, Ncrit and xtr from a synthetic fixture and refuses a file that is not one; converted files carry the hand-added settings; the metrics recover a known line, α_L0 and cd_min exactly; a polar compared with itself gives zero and with a shifted copy gives the shift; the Ladson file is found and read as published; the own panel code gives cl = 0 for a symmetric section at α = 0, tends to 2πα as thickness → 0 and converges at order 2; the drag estimate falls with Re and rises when tripped; the own codes sit where the physics says against Ladson; an own-code polar round-trips with its own label. Nothing in the tests imports or runs XFOIL, AVL or OpenVSP; the polar fixtures are SYNTHETIC files written by the tests. |
| **Practical activity** | **computational experiment**: the own codes against Ladson's free- and fixed-transition data at 2.0 and 5.97 × 10⁶ (and the learner's own Reynolds number), with a sensitivity study of the settings that matter; optional practice for learners who install the free tools, storing their output as labelled CSVs and running the same comparison. |
| **Equipment** | None required. Optional: XFOIL 6.99 (or 6.996 from source) and AVL 3.52 from Drela's pages, OpenVSP 3.55 with VSPAERO; all free. |
| **Safety** | None beyond normal computer use. Download the tools only from their official pages. |
| **Evidence** | Completed `main.py`, passing tests, `tool_validation.png`, the validation tables at two Reynolds numbers, the settings-sensitivity table, and a half-page validation statement in the form of Problem 6(d) (or, if you ran the tools, your labelled CSVs and their comparison). |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the reader refuses unlabelled files; all public tests pass; the own-code validation is reported with cl_α, α_L0, RMS cl over a stated range and cd_min, each difference attributed; tool output, if any, is labelled and never confused with computation or measurement. |
| **References** | `@drela-xfoil` XFOIL page (https://web.mit.edu/drela/Public/web/xfoil/): Release Conditions; Software; Documents; User Primer `xfoil_doc.txt`: Program Execution; Airfoil file formats; Analysis; Numerical accuracy; Polar calculations and plotting; Re, Mach dependence; Caveats; "Viscous Formulation" (supplementary-free; GNU GPL). `@drela-avl` AVL page (https://web.mit.edu/drela/Public/web/avl/) and User Primer `avl_doc.txt`: Input Files; Geometry Input File -- xxx.avl (Coordinate system; File format); Program Execution; OPER Routine -- Flow Analysis (supplementary-free; GNU GPL). `@openvsp-docs` OpenVSP (https://openvsp.org/): Download; License; wiki "Analysis Guides" (VSPAERO VLM & Panel Code; Definitions of Aerodynamic Coefficients); "vspaerotutorial" (supplementary-free; NASA Open Source Agreement 1.3, wiki CC BY-SA 4.0). `@ntrs-nasa-tm-4074-ladson-naca0012` Ladson, NASA TM-4074: Wake Survey Rake; Tests and Procedures; Presentation of Results (required-free): the reference data. `@nasa-nparc-cfd-vv-tutorial` (required-free): the vocabulary of validation used in `P03.18`. Optional textbooks for deeper study (paid, never required): Drela, *Flight Vehicle Aerodynamics*, for the strongly coupled viscous–inviscid formulation XFOIL uses; Katz & Plotkin, *Low-Speed Aerodynamics*, for panel and lattice methods; library or used copies are fine. |
| **Extensions** | ★ Add displacement-thickness feedback (transpiration velocities v = d(Ue δ*)/ds on the panels) to the own code and iterate; show that the lift-slope error against Ladson falls. ★ If you have AVL, run the `P10.11` wing and compare CL_α, CDi and e with your lattice in the same labelled-CSV format, with the AVL settings in the header. ★ Write a reader for VSPAERO's `.polar` output that produces the same labelled CSV. |
| **Time** | unmeasured |

## 1. The engineering problem

XFOIL, AVL and OpenVSP are free, fast and good, and engineers use them every day: XFOIL for airfoil polars with transition and drag, AVL for vortex-lattice loads and stability derivatives, OpenVSP for building geometry and running its VSPAERO lattice and panel solvers. They are excellent checks on your own codes, and your own codes are excellent checks on them. They also fail in quiet ways. The same airfoil gives a different polar at a different transition parameter, a different panel count or an unconverged point, and a CSV file of numbers says none of that. A drag coefficient lifted from an unlabelled polar into an aircraft sizing is a common source of optimistic designs.

This lesson does three things. It gives the exact commands and settings to run the tools, for learners who choose to install them. It builds the bookkeeping (a reader that refuses any polar whose header does not say which tool, which version and which settings produced it) and comparison metrics that mean something. And it validates the learner's own codes from `P10.05` and `P10.06` against the published wind-tunnel data of Ladson (`@ntrs-nasa-tm-4074-ladson-naca0012`), which is the reference both your code and the tools answer to. **The tools were not available in the environment in which this lesson was written, and no tool output ships with it.** Nothing in the tests imports or runs them.

Success criteria, written first: the reader refuses unlabelled files and accepts labelled ones; the metrics recover known lines exactly; the own codes reproduce the symmetric-section and thin-airfoil limits and converge; against Ladson, every difference in cl_α, α_L0 and cd is quantified and given a cause.

## 2. Understand and derive

### 2.1 What the tools compute, and what their settings do

**XFOIL** (`@drela-xfoil`) is an inviscid linear-vorticity panel method strongly coupled to a two-equation integral boundary layer and wake, with transition predicted by an e^N envelope method. The settings that change a polar:

| setting | XFOIL command | what it changes |
|---|---|---|
| Reynolds number | `VISC <Re>` in OPER | boundary-layer thickness, transition, cd, cl_max |
| Mach number | `MACH <M>` | a compressibility correction to the inviscid solution (about +1 % lift slope at 0.15) |
| Ncrit | `VPAR`, then `N <Ncrit>` | the amplification exponent at which transition is predicted; 9 is the usual default for a quiet tunnel |
| forced transition | `VPAR`, then `XTR <top> <bot>` | trips the boundary layer at those x/c |
| panel count | `PPAR`, then `N <n>` | resolution of the leading-edge peak and the boundary-layer march |
| iteration limit | `ITER <n>` in OPER | how many Newton iterations before a point is declared unconverged |

**AVL** (`@drela-avl`) is a vortex-lattice code (like `P10.11`) with slender-body fuselages and a stability-derivative output; its settings are the lattice (chordwise and spanwise counts and spacing), the reference area, chord and span, the moment reference point, and the Mach number (a Prandtl–Glauert correction). **OpenVSP** (`@openvsp-docs`) builds parametric geometry and runs VSPAERO, a lattice or panel solver; its settings include the method, the tessellation, the Mach and Reynolds numbers and the number of wake iterations.

### 2.2 A polar is data only with its header

A labelled polar CSV, as this lesson defines it, starts with `#` lines:

```
# tool: XFOIL
# version: 6.99
# settings: re=6e6; mach=0.15; ncrit=9; panels=160; iter=100; xtr_top=1.0; xtr_bot=1.0
# airfoil: NACA 0012
# produced_by: <your name>, <date>, <machine and OS>
alpha_deg,cl,cd,cm
```

`read_tool_polar` refuses a file without `tool`, `version` and `settings`, and refuses an XFOIL polar whose settings lack any of re, mach, ncrit, panels, iter (AVL: mach, nchord, nspan, sref, cref, bref; OpenVSP: mach, re, method, wake_iter, tessellation). XFOIL's own polar file records the version, the airfoil name, xtrf, Mach, Re and Ncrit; it does not record the panel count or the iteration limit, so the converter asks for them.

### 2.3 Metrics

Over a stated linear range, fit cl = m α + c by least squares (`P02.18`): **cl_α = m**, **α_L0 = −c/m**. To first order u²(α_L0) = [u²(c) + (c/m)²u²(m) + 2(c/m)Cov(m, c)]/m², so α_L0 is well determined only when the range straddles it. Compare two polars at common angles: interpolate the second onto the first's angles inside the linear range and form **RMS Δcl**. A constant slope ratio 1 + ε gives RMS Δcl = εm√⟨α²⟩, which grows with the range: always quote the range. Add the relative differences in cl_α and in **cd_min**.

### 2.4 What to expect from the own codes

An inviscid panel method gives a NACA 0012 lift slope about 2π(1 + 0.77 t/c) ≈ 6.87 per radian, 0.120 per degree. A real section has less, because the boundary layer's displacement thickness, thicker on the suction side, acts like a loss of camber and circulation. One-way coupling computes the boundary layer from the inviscid velocity but never feeds it back, so it cannot capture that loss: expect the own code's cl_α to be 10–15 % high. The drag estimate depends on where transition happens; recap of the two correlations it uses (`P08.11`, `P08.12`, written here because Phase 08 is not yet available):

- **Thwaites** (laminar): θ² = (0.45ν/Ue⁶)∫₀ˢ Ue⁵ ds, λ = (θ²/ν) dUe/ds, laminar separation at λ ≈ −0.09.
- **Michel** (transition): Re_θ > 1.174(1 + 22 400/Re_x) Re_x^0.46.
- Turbulent: dθ/ds = cf/2 − (H + 2)(θ/Ue) dUe/ds with H = 1.4 and cf = 0.246 × 10^(−0.678H) Re_θ^(−0.268); drag by **Squire–Young**, cd = Σ 2θ (Ue/V)^((H+5)/2), evaluated here at x/c = 0.99 because the inviscid velocity on the last panel plunges toward a trailing-edge stagnation point that the real, viscous flow does not have (moving the station between 0.98 and 0.995 changes cd by about ±1.4 %).

## 3. Build the computational model

Open `starter/main.py`. The data loader, the NACA generator, the panel method (factorised once per airfoil and reused for every angle), the boundary-layer march, `bl_drag`, `write_tool_polar` and `main` are given; their docstrings say which lesson each comes from.

### Step 1: the reader and the parser

`parse_header` turns `# key: value` lines into a dict and splits `settings` into key = value pairs. `read_tool_polar` applies the refusals of §2.2, then reads the columns with `csv.DictReader` (blank cells become NaN). `parse_xfoil_polar` uses regular expressions on XFOIL's header lines (`Version`, `Calculated polar for:`, `xtrf = … (top) … (bottom)`, `Mach = … Re = … e … Ncrit = …`) and reads the seven-column table after the dashed line. `xfoil_polar_to_labelled_csv` adds the panel count, the iteration limit and who ran it.

### Step 2: metrics and comparison

`polar_metrics` fits the line with `np.polyfit` over the linear range and finds cd_min; `compare_polars` interpolates and differences as in §2.3.

### Step 3: the own codes as a labelled polar

`own_polar` runs the panel method and `bl_drag` at each angle and returns a polar dict *with its own header*: tool "own code (P10.05 panel + P10.06 boundary layer)", version, and settings. `validate_against_ladson` runs it at the published angles of one Ladson block and compares.

### Verify

The tests (14) check everything listed in the metadata; the polar fixtures they use are SYNTHETIC files they write into `tmp_path`, labelled as such in their headers. The own panel code's convergence at α = 4°:

| panels | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| cl | 0.481970 | 0.482466 | 0.482591 | 0.482622 |

Successive differences 4.96 × 10⁻⁴, 1.25 × 10⁻⁴, 3.10 × 10⁻⁵: **observed order 1.99, 2.01**.

### Compare: the own codes against the published data

Ladson's free-transition data first, fitted over −4.5° to 8.5°: cl_α = 0.10355 per degree (5.933 per radian), α_L0 = −0.041°, cd_min = 0.00620 at Re = 2.0 × 10⁶; cl_α = 0.10688 per degree (6.124 per radian), α_L0 = −0.008°, cd_min = 0.00650 at 5.97 × 10⁶. The own codes (200 panels) against them:

| Re | cl_α own / published [/deg] | Δcl_α | Δα_L0 | RMS Δcl (max) | cd_min own / published | Δcd_min |
|---|---|---|---|---|---|---|
| 2.0 × 10⁶ | 0.12049 / 0.10355 | +16.4 % | +0.042° | 0.0760 (0.126) | 0.00734 / 0.00620 | +18.5 % |
| 5.97 × 10⁶ | 0.12049 / 0.10688 | +12.7 % | +0.009° | 0.0622 (0.109) | 0.00627 / 0.00650 | −3.5 % |

Drag at 5.97 × 10⁶, angle by angle (own code with Michel transition):

| α [deg] | −4.05 | −2.00 | 0.05 | 1.98 | 4.18 | 6.20 | 8.22 |
|---|---|---|---|---|---|---|---|
| published cd | 0.00700 | 0.00650 | 0.00650 | 0.00680 | 0.00760 | 0.00680 | 0.00800 |
| own cd | 0.00690 | 0.00642 | 0.00627 | 0.00647 | 0.00689 | 0.00768 | 0.00935 |
| own x_tr (upper) | 0.632 | 0.461 | 0.294 | 0.164 | 0.082 | 0.044 | 0.018 |

The lift slope is too steep by the missing displacement effect (§2.4), more at the lower Reynolds number where the boundary layer is thicker: exactly the trend expected. The zero-lift angle, which a symmetric section fixes by geometry, agrees within 0.05°. The drag is within 1–10 % at 5.97 × 10⁶ up to 4°; above 6° the predicted transition runs to the leading edge while the published drag falls (the report flags the 6.20° point as lower than at 4.18°, a laminar-bucket and transition effect of the free-transition model). At 2 × 10⁶ the estimate is 15–25 % high, because Michel's criterion trips the flow too early there. With transition forced at 0.05c the own code gives cd_min = 0.00821 against 0.00895 for the wraparound grit strip at 6.00 × 10⁶ (−8.3 %): the grit adds roughness drag a trip model does not.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment**, with optional practice on free tools.

**Safety first:** none beyond normal computer use; download tools only from the official pages named in the References.

**4.1 The experiment everyone does (no tools needed).**

1. The validation table above at 2.0 and 5.97 × 10⁶, then at your own Reynolds number (Problem 4). Predict the sign of each difference before running.
2. Sensitivity of the own code to its settings at 5.97 × 10⁶: panel count (80, 160, 320 panels: cd at 0° = 0.00596, 0.00622, 0.00633, cl at 4° changing only in the fourth figure), Squire–Young station (0.98–0.995: ±1.4 %), transition (free 0.00627 against tripped at 0.05c 0.00822). *Finding:* transition ≫ panels > Squire–Young station. These are the settings a header must carry.

**4.2 Optional: generating labelled tool polars.** No tool output ships with this lesson because the tools were not available when it was written. A learner who installs them can produce polars and run exactly the same comparison. Name each file `<tool>_<airfoil>_Re<Re>_<version>.csv`, for example `xfoil_naca0012_Re6e6_6.99.csv`, and give it the header of §2.2.

*XFOIL* (6.99; type each line, an empty line means press Enter):

```
NACA 0012
PPAR
N 160

OPER
VISC 6e6
MACH 0.15
VPAR
N 9

ITER 100
PACC
xfoil_naca0012_Re6e6_M015_N9.pol

ASEQ -4 16 0.5
PACC
QUIT
```

`PPAR` / `N 160` sets the panel count (Enter twice leaves the menu and repanels); `VISC` turns on the viscous solver at that Re; `VPAR` / `N 9` sets Ncrit (for a tripped case add `XTR 0.05 0.05` in VPAR); `ITER 100` raises the iteration limit; `PACC` opens the polar file (the second prompt, a dump file, can be left empty); `ASEQ` sweeps α. Then convert: `xfoil_polar_to_labelled_csv(open(".pol").read(), "xfoil_naca0012_Re6e6_6.99.csv", panels=160, iter_limit=100, produced_by="<name>, <date>")`, list any angles XFOIL reported as unconverged in a `# note:` line, and compare with `compare_polars(load_ladson("I", 5.97e6), read_tool_polar(...))`. Check the prompts against `xfoil_doc.txt` and `sessions.txt` for your version.

*AVL* (3.52). Geometry file `trainer.avl` for the reference trainer's wing (rectangular, b = 10.0 m, c = 1.2 m, NACA 2412):

```
Reference trainer wing (P10.21)
#Mach
0.0
#iYsym  iZsym  Zsym
0       0      0.0
#Sref   Cref   Bref
12.0    1.2    10.0
#Xref   Yref   Zref
0.3     0.0    0.0
#
SURFACE
Wing
#Nchord  Cspace  Nspan  Sspace
8        1.0     24     -2.0
YDUPLICATE
0.0
SECTION
#Xle  Yle  Zle  Chord  Ainc
0.0   0.0  0.0  1.2    0.0
NACA
2412
SECTION
#Xle  Yle  Zle  Chord  Ainc
0.0   5.0  0.0  1.2    0.0
NACA
2412
```

Commands: `avl trainer.avl`, then `OPER`, `A A 5` (set α = 5°), `X` (execute), and read CLtot, CDind and e from the total-forces listing; repeat for the α values you need, and record Nchord = 8, Nspan = 24 (per half, with YDUPLICATE), the spacing parameters and the reference values in the header. Compare with your `P10.11` lattice and the `P10.20` free-air values.

*OpenVSP* (3.55, VSPAERO). Add a Wing; set span 10 m (full), root and tip chord 1.2 m, sweep 0, airfoil NACA 4-series 2412; Analysis → VSPAERO; choose the vortex-lattice method; set Sref 12, cref 1.2, bref 10; α from −4 to 12 in 9 steps, Mach 0, Re 6 × 10⁶; three wake iterations; run, then write CL, CDi and CMy from the results into a labelled CSV (record method, tessellation counts and wake iterations). The OpenVSP wiki's VSPAERO tutorial is partly out of date, so treat these GUI steps as a guide and record what you actually set.

Licences: XFOIL and AVL are released under the GNU GPL; OpenVSP under the NASA Open Source Agreement 1.3. This curriculum links to them and does not bundle them.

Raw data format: labelled CSV polars as in §2.2; the validation tables `alpha_deg,cl_pub,cl_own,cd_pub,cd_own,xtr_own`; the sensitivity table `setting,value,cl_4deg,cd_0deg,cd_4deg`.

## 5. Validate and improve

| quantity | published (Ladson, 5.97 × 10⁶) | own code | difference | cause |
|---|---|---|---|---|
| cl_α [/deg] | 0.1069 | 0.1205 | +12.7 % | no displacement feedback (one-way coupling) |
| α_L0 | −0.008° | ≈ 0 | 0.009° | symmetric section; data scatter |
| cd_min | 0.00650 | 0.00627 | −3.5 % | transition model, fixed-H turbulent closure |
| cd at 8.2° | 0.00800 | 0.00935 | +17 % | predicted transition at the leading edge |
| tripped cd_min (grit) | 0.00895 (6.00 × 10⁶) | 0.00821 (trip at 0.05c) | −8.3 % | grit roughness drag not modelled |
| panel cl convergence | — | order 2.0 | — | verified |

The engineering answer: **your own codes reproduce the zero-lift angle and, at 6 × 10⁶, the minimum drag within a few per cent, but overpredict the lift slope by 11–16 % (worse at lower Re) because one-way coupling leaves out the boundary layer's displacement; a tool polar is usable as a cross-check only when its header records the tool, version and settings, and the wind-tunnel data remain the reference.** Improve next: displacement-thickness feedback (Extension 1), which addresses the largest difference. Not worth improving for this purpose: the panel count beyond 200, which changes cl in the fourth figure.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; `tool_validation.png`; the validation tables at two Reynolds numbers and your own; the settings-sensitivity table; and the validation statement (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7). If you ran any tool, add its labelled CSVs and their comparison; never present your own computation or a test fixture as tool output.

### Adversarial self-checks

- **What if the thickness goes to zero?** The panel lift slope tends to 2π per radian (the tests check NACA 0001 to 2 %), and the viscous correction shrinks with the boundary layer's effect on the shape.
- **What if Re goes up?** The boundary layer thins and the published slope rises (0.1036 → 0.1069 → 0.1087 per degree from 2 to 8.86 × 10⁶) toward the inviscid one; the own code's error falls from 16 % to 11 %.
- **Which assumption fails first?** One-way coupling, for lift; the transition model, for drag.
- **What measurement could falsify the code?** A published cl_α above the inviscid panel value would mean the panel code (or the data) is wrong; a drag below the laminar flat-plate value would mean the same.
- **Does the result survive a change of settings?** cl does; cd moves 6 % with panel count and 31 % with the transition setting. That is why the settings go into the header.
- **Which decision changes if the uncertainty doubles?** None of the qualitative conclusions; a drag estimate used for performance would need the tripped (conservative) value rather than the free-transition one.

## Check yourself

1. Why must a polar file record the tool, its version and its settings?
2. Name the five XFOIL settings this lesson's reader requires.
3. What tool output ships with this lesson, and what do the tests use instead?
4. Which four metrics does the lesson use to compare two polars?
5. By how much does your inviscid panel code overpredict the NACA 0012 lift slope against Ladson's data, and what physical effect is missing?
6. How close is the one-way drag estimate to the published minimum drag, and which setting dominates its uncertainty?
7. Under which licences are XFOIL, AVL and OpenVSP released?
