# Author brief: Phase 08 lessons

This brief is what every lesson author in Phase 08 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-08.md`, `references/sources.yaml`, the Phase 08 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/`, `solutions/_template/` and `projects/_template/` folders. **Phases 00–05 are the quality bar:** before writing, read `curriculum/phase-05/lesson-12-beam-deflection/` and `curriculum/phase-03/lesson-09-iterative-sparse-solvers/` (lesson, problems, starter, tests and the matching `solutions/` folder) end to end and match their depth, tone and structure. For the home lab, also read `curriculum/phase-05/lesson-18-home-lab-cantilever-deflection/lesson.md` and every file in `projects/p05-cantilever/`; that is the bar for the experiment, the bill of materials and the safety box.

## Where files go

```
curriculum/phase-08/lesson-NN-<slug>/lesson.md
curriculum/phase-08/lesson-NN-<slug>/problems.md
curriculum/phase-08/lesson-NN-<slug>/starter/main.py
curriculum/phase-08/lesson-NN-<slug>/tests/test_main.py
solutions/phase-08/lesson-NN-<slug>/hints.md
solutions/phase-08/lesson-NN-<slug>/method.md
solutions/phase-08/lesson-NN-<slug>/worked-solutions.md
solutions/phase-08/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
projects/p08-flow-pressure/bom.yaml, safety.md, procedure.md, report-template.md, data/README.md   <- lesson 19 only
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build | mode of part 4 |
|---|---|---|---|---|
| 01 | fluid-properties-hydrostatics | Fluid Properties & Hydrostatics | Standard atmosphere from scratch (ISA layers to 47 km by integrating the hydrostatic equation, checked against the closed forms and the ISA anchor values; viscosity by Sutherland's law; manometer and buoyancy) | computational experiment |
| 02 | fluid-kinematics | Fluid Kinematics | Streamline, pathline and streakline plotter (unsteady 2-D field; the three lines coincide for steady flow; material derivative, strain rate and rotation) | computational experiment |
| 03 | control-volume-mass-momentum | Control Volume: Mass & Momentum | Jet thrust from a momentum balance (Reynolds transport theorem; mass and momentum balances that close to round-off; pressure thrust; moving control volumes) | computational experiment |
| 04 | bernoulli-energy-equation | Bernoulli & Energy Equation | Pitot-static airspeed (Bernoulli along a streamline; IAS, CAS, EAS and TAS; static-port error; the energy equation with losses and shaft work) | computational experiment |
| 05 | dimensional-analysis-similarity | Dimensional Analysis & Similarity | Buckingham Pi tool (nullspace of the dimension matrix, by SymPy and by SVD) and scaling a wind-tunnel model (Re and Mach cannot both be matched) | computational experiment |
| 06 | navier-stokes-equations | Navier-Stokes Equations | Derivation and exact solutions (Couette, plane and pipe Poiseuille, Couette–Poiseuille, Stokes's first problem; each checked by substitution into the equations and by a finite-difference solve) | computational experiment |
| 07 | pipe-flow-losses | Pipe Flow & Losses | Propellant feed-line pressure drop (Darcy–Weisbach, laminar 64/Re, Colebrook solved by Newton against the Haaland explicit form, minor losses, a tank-to-injector budget) | computational experiment |
| 08 | potential-flow | Potential Flow | Superposition of elementary flows (uniform flow, source, sink, doublet, vortex; Rankine oval and lifting cylinder; zero normal flow on the body, Kutta–Joukowski by contour integration, d'Alembert's paradox) | computational experiment |
| 09 | vorticity-circulation | Vorticity & Circulation | Kelvin's theorem demo (circulation round a material loop advected through a vortex flow stays constant in inviscid flow and decays as Lamb–Oseen predicts with viscosity; Stokes's theorem; Biot–Savart; the starting vortex) | computational experiment |
| 10 | laminar-boundary-layers | Boundary Layers: Laminar | Blasius solution (shooting with RK4 and a root finder, checked against `scipy.integrate.solve_bvp`; δ99, δ*, θ, H, c_f) | computational experiment |
| 11 | integral-boundary-layer-methods | Integral Boundary Layer Methods | Thwaites method (momentum integral; Thwaites against Blasius on a flat plate and against the Hiemenz stagnation-point solution; Falkner–Skan wedges; separation prediction on a cylinder and a decelerating flow) | computational experiment |
| 12 | transition | Transition | e^N method intro (Orr–Sommerfeld for a plane Poiseuille flow by a Chebyshev collocation eigenvalue solver, checked against the published critical Reynolds number 5772.22; an e^N envelope for a flat plate from a supplied stability chart labelled as an approximate fit; N = 9) | computational experiment |
| 13 | turbulence-fundamentals | Turbulence Fundamentals | Turbulence statistics from data (Reynolds decomposition, moments, autocorrelation and integral scale, Taylor's hypothesis, spectra with the −5/3 range, Kolmogorov scales) | data-based investigation (authentic public dataset with source and license, or a seeded synthetic signal labelled synthetic everywhere) |
| 14 | turbulent-boundary-layers | Turbulent Boundary Layers | Log law and skin friction (law of the wall, Spalding's composite profile, Clauser fit of u_τ, power-law and 1/7th-law estimates, flat-plate c_f correlations, the mixed laminar–turbulent plate) | computational experiment (with a public DNS mean profile if the reference map verifies one, otherwise synthetic data labelled synthetic) |
| 15 | turbulence-modeling | Turbulence Modeling | Mixing-length and k-ε overview (RANS closure; Prandtl mixing length with van Driest damping integrated across a channel and compared with the log law and the friction law; the standard k-ε constants and the log-layer equilibrium check) | computational experiment |
| 16 | flow-separation-drag | Flow Separation & Drag | Drag of bluff bodies (pressure and friction drag; c_d(Re) of a sphere from the Schiller–Naumann and Morrison correlations, labelled correlations; drag crisis; Strouhal shedding; terminal velocity; streamlining) | computational experiment |
| 17 | intro-compressible-flow | Intro to Compressible Flow | Speed of sound and Mach number (a = √(γRT) from the ISA of lesson 01; isentropic total-to-static relations; where the Bernoulli airspeed of lesson 04 fails; compressible pitot below M = 1; the Rayleigh pitot formula above it) | computational experiment |
| 18 | incompressible-cfd-lid-driven-cavity | Incompressible CFD from Scratch | 2-D lid-driven cavity solver (finite differences, projection or stream function–vorticity; Re = 100 against the Ghia, Ghia & Shin 1982 centreline table, labelled as published data; grid convergence with observed order; mass conservation) | computational experiment |
| 19 | home-lab-flow-pressure-measurement | Home Lab: Flow & Pressure Measurement | Low-energy flow rig: a draining bottle (Torricelli) with a discharge-coefficient fit, and a water U-tube manometer read against a guarded fan or a hair dryer on cool, to estimate airspeed; Reynolds-number reasoning and calibrated pressure 🔬 | physical experiment, `projects/p08-flow-pressure` |

Lesson ids are `P08.NN`. Prerequisites: every lesson names the previous Phase 08 lessons it uses and the earlier lessons it builds on, in the metadata row. The ROADMAP's `Needs` column fixes the cross-phase ones: `P07.01` (systems, properties and state; ideal gas) for 01, `P01.11` (vector fields, line and surface integrals) for 02, `P02.11` (PDEs) for 06, `P02.12` (complex variables) for 08, `P03.05` (RK4) for 10, `P02.19` (random processes and power spectra) for 13, `P07.05` (isentropic processes) for 17, `P03.12` (finite volume methods) for 18, `P00.14` (measurement and uncertainty) for 19. Phase 07 is being written in parallel; name its lessons by id and title (`P07.01` Systems, Properties & State; `P07.05` Isentropic Processes & Efficiencies) and restate in two or three sentences anything you use from them, because the learner may arrive before Phase 07 is published. Other earlier lessons worth naming where used: `P00.03` (vectorized ISA table), `P00.05` (units), `P00.13` (home-experiment safety), `P00.15` (first experiment), `P01.08` (partial derivatives), `P01.12` (divergence and Stokes theorems), `P01.14` (rank and nullspace), `P01.16` (eigenvalues), `P01.17` (least squares), `P01.18` (SVD), `P02.02` (second-order ODEs), `P02.10` (FFT), `P02.18` (regression), `P03.02` (root finding), `P03.03` (interpolation), `P03.04` (quadrature), `P03.06` (adaptive and stiff integrators), `P03.08`–`P03.09` (direct and iterative linear solvers), `P03.10` (Newton for systems), `P03.11` (finite differences, stability, observed order), `P03.18` (verification and validation), `P04.01`–`P04.03` (Newton's laws, energy, momentum), `P04.07` (kinetic theory). The learner has finished Phases 00–05: Python, the scientific stack, calculus, linear algebra, ODEs, PDEs, numerical methods, physics and statics are known. What is new is fluid mechanics: continuum fluid properties, kinematics, conservation laws in control-volume and differential form, viscous and inviscid flows, boundary layers, turbulence, an introduction to compressibility, and a first CFD code.

## What Phase 08 must achieve

Fluid mechanics at the level of a first undergraduate fluids course and the fluids part of an aerospace unified-engineering sequence (MIT 2.06 and the Unified Engineering fluids notes; Çengel–Cimbala or White scope), taught through aircraft and rocket examples. The organising idea of the phase is *conservation of mass, momentum and energy, applied to a control volume or to a fluid particle, with a model for the stress*. Every lesson says which conservation law it uses, in which form (integral or differential), and what it assumes about the stress (hydrostatic, inviscid, Newtonian laminar, or a turbulence model). Every main result is checked by an independent route, for example:

- mass and momentum balances that close to round-off;
- the standard atmosphere checked against the ISA anchor values;
- Poiseuille and Couette exact solutions against a finite-difference solve;
- Blasius by shooting against `scipy.integrate.solve_bvp`;
- Thwaites against Blasius on a flat plate;
- potential-flow superposition satisfying the zero normal flow condition on the body;
- the Orr–Sommerfeld critical Reynolds number against the published value;
- the lid-driven cavity against the Ghia et al. 1982 benchmark, labelled as published data;
- grid convergence with observed order wherever a grid is used.

A learner who finishes can compute fluid properties and a standard atmosphere, describe a flow field and its lines, write mass and momentum balances for a control volume, measure airspeed with a pitot-static tube and say when Bernoulli fails, find dimensionless groups and scale a model test, solve the exact viscous flows, size a feed line, build potential flows and the lift on a cylinder, reason with vorticity and circulation, compute laminar boundary layers exactly and by an integral method, estimate transition, analyse a turbulent signal, use the log law and skin-friction correlations, explain what a turbulence model closes, estimate bluff-body drag, handle low-speed compressibility, write and verify a small incompressible Navier–Stokes solver, and measure a discharge coefficient and an airspeed at home and compare them with theory, with a tested Python module for each.

### Shared numbers (use these so lessons connect)

- **ISA sea level:** T₀ = 288.15 K, p₀ = 101 325 Pa, ρ₀ = 1.225 kg/m³, R = 287.05287 J/(kg·K), g₀ = 9.80665 m/s², γ = 1.4, lapse −6.5 K/km to 11 km, isothermal 216.65 K to 20 km, +1.0 K/km to 32 km, +2.8 K/km to 47 km (geopotential altitude). Anchor values for checking: p = 22 632.06 Pa at 11 km, 5 474.889 Pa at 20 km, 868.019 Pa at 32 km, 110.906 Pa at 47 km; ρ(11 km) = 0.36392 kg/m³. Sutherland: μ = 1.458 × 10⁻⁶ T^{3/2}/(T + 110.4) Pa·s (μ₀ = 1.7894 × 10⁻⁵ Pa·s at sea level).
- **The light aeroplane** (the same one the earlier phases used for its wing): cruise at 3000 m pressure altitude at a true airspeed of 60 m/s, mean wing chord 1.5 m, semi-span 5 m. Its pitot-static system, its wing boundary layer, its Reynolds number and its model in a wind tunnel recur in lessons 04, 05, 10, 11, 12, 14, 16 and 17.
- **Water at 20 °C:** ρ = 998.2 kg/m³, μ = 1.002 × 10⁻³ Pa·s. **Air at 20 °C, 101 325 Pa:** ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s.
- Fluid properties for propellants (lesson 07) are typical values, stated with their source and labelled typical.

Numbers must be realistic, and every number in `lesson.md` and `worked-solutions.md` must be reproduced by `solution.py` or by a probe script that calls it (keep probe scripts in your scratchpad, uniquely named, e.g. `probe_p08_07.py`; they are not committed).

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-08/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-08/lesson-NN-<slug>` runs against `solutions/phase-08/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish. Also run plain `python3 -m pytest curriculum/phase-08/lesson-NN-<slug>`: every failure must be a `NotImplementedError` (no import errors, no syntax errors, no tests that pass or fail for other reasons against the starter except tests of functions the starter supplies complete). Run `python3 scripts/check_lessons.py` (structure linter); it must report no ERROR for your lessons.

The starter has the same function signatures and docstrings as the solution, with `raise NotImplementedError` bodies, plus a working `main()` (and any helpers the lesson says are given) that produces the lesson's output once the learner completes the functions. Tests assert properties, not one magic number:

- conservation residuals (mass and momentum balances to round-off; ∇·u = 0 to discretisation order);
- exact solutions substituted back into the governing equations and boundary conditions;
- limiting cases (zero velocity, zero viscosity, infinite Reynolds number, an isothermal layer);
- agreement between two methods (shooting against `solve_bvp`; Colebrook by Newton against Haaland; SymPy against SVD nullspace);
- published anchor values (ISA, Blasius f″(0) = 0.33206, the critical Reynolds number 5772.22, the Ghia table), labelled as published;
- observed order of convergence from a grid sequence;
- dimensional sanity and invariance (results unchanged by rescaling units).

12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. **Each lesson's tests run in under 10 s** (check with `--durations=5`); the CFD lesson uses a small grid (for example 17 × 17 and 33 × 33 at Re = 100 in tests, with a finer sweep only in `main()` or the experiment). Anything random is seeded (`numpy.random.default_rng(seed)`) with statistical tolerances justified in a comment.

Code: Python 3.11+, NumPy, SciPy, Matplotlib and SymPy only (no CFD, fluids or atmosphere packages in tested code). Type hints, docstrings, `if __name__ == "__main__": main()`. Matplotlib uses the `Agg` backend (`matplotlib.use("Agg")` before importing pyplot) in any code that runs under tests; `main(out_dir)` saves figures to a path and closes every figure (`plt.close(fig)`), never `plt.show()`.

## Data honesty

- **Synthetic data is labelled synthetic everywhere it appears:** in the file name or header, in the docstring that generates it, in the lesson text, in figure titles, and in worked solutions.
- **Turbulence data** (lessons 13–15) is either an authentic public dataset whose source, license and measurement method are named (only if `references/phase-08.md` verifies one; the dataset itself must not be downloaded in tests) or a seeded synthetic signal or profile, labelled synthetic. A synthetic turbulence signal is generated with a stated model spectrum (for example a von Kármán spectrum) so the lesson can check that its analysis recovers the known integral scale and variance.
- **Published data** (the ISA anchor values, the Ghia et al. Re = 100 centreline table, the critical Reynolds number, drag correlations) is labelled "published data" or "correlation" with its citation; a correlation is never called a measurement.
- **Simulation is never presented as measurement.**

## Physical experiment (lesson 19)

- Lesson 19 owns `projects/p08-flow-pressure/` with `bom.yaml` (three tiers: `minimum_cost`, `recommended`, `enhanced`, as in the earlier kits; prices as ranges with the date checked; `reuse_from` naming the p00 and p04 kits it reuses, e.g. `p00-calibration`, `p00-first-experiment`, `p04-pendulum`, `p04-drag-drop`), `safety.md` (the full safety box in the format of `projects/p05-cantilever/safety.md`), `procedure.md`, `report-template.md` and `data/README.md` (raw-data formats), plus a seeded synthetic dataset under `data/synthetic/` labelled synthetic in its file names and headers.
- **Experiment 1, draining bottle (Torricelli).** A plastic bottle with a small round hole near the base, filled with room-temperature water (≤ 2 L) and drained into a basin or sink; the water level h(t) is read against a ruler taped to the bottle, timed with a phone video. The learner predicts the draining time from Torricelli with a discharge coefficient C_d, fits C_d from h(t) (the √h law is linear in t), and compares it with the textbook range for a sharp-edged orifice (labelled typical). Reynolds number of the jet, and why viscosity barely matters.
- **Experiment 2, U-tube manometer and airspeed.** A clear tube bent into a U, part-filled with water (food colouring optional), one leg connected to a pitot probe (a drinking straw or a small tube facing the flow) and the other open to still air or connected to a static port; held in the stream of a small guarded fan or of a hair dryer on its **cool** setting. Δh is read, Δp = ρ_w g Δh, V = √(2Δp/ρ_air). The learner states the resolution limit (a 1 mm reading at a few metres per second), computes the Reynolds number of the probe and of the fan, and decides whether an inclined manometer improves it. A second check of airspeed by an independent route where possible (the enhanced tier: a cheap anemometer; otherwise the fan's volume flow from its rating, labelled as a manufacturer's figure).
- The prediction, with its uncertainty, is written down and dated before measuring, in both experiments.
- **Hard limits, stated in `safety.md` and repeated in the lesson:** water at room temperature, no more than 2 L in total; fans must be guarded (finger-proof grille) and draw at most 30 W; a hair dryer is used on its cool setting only and is kept away from water (separate tables, the manometer's water sealed or below the dryer, never in the same basin); no compressed gas, no pressurised containers (no squeezing or pressurising bottles to drive flow), no mains appliance near spilled water; mop spills before switching anything on; stop criteria.
- The learner without a kit completes the lesson with the supplied synthetic dataset, labelled synthetic everywhere it appears.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`. Part 4 heading is exactly `## 4. Experiment at home, or investigate real data`; the practical-activity row names one of: physical experiment, computational experiment, data-based investigation, design-build.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phases 00–05. Plain language, short sentences, define each term on first use. SI units (knots and feet only where aviation practice demands, and say so).
- Every lesson states a concrete *engineering question* in the first paragraph (and in the metadata row, with the answer), derives in part 2, builds in part 3, runs an **experiment** in part 4 (physical in lesson 19; otherwise a computational experiment such as a parameter sweep, a sensitivity study, a method comparison or a grid-convergence study, or a data-based investigation), validates in part 5, and ends with `### Adversarial self-checks` and `## Check yourself` with 5–8 questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: the light aeroplane's pitot-static system, wing and wind-tunnel model; a turbojet or rocket nozzle's thrust; a liquid-propellant feed line from tank to injector; a starting vortex behind a wing; a sounding rocket's drag; a cabin air duct; a wind-tunnel contraction.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-08.md`; tag each required-free / supplementary-free / optional-paid. Where the map says GAP or partial GAP, write original instruction and say nothing about missing references. Never cite or require a paid resource; paid books only under an "Optional textbooks for deeper study" sentence or line, never on a line containing "required". Never paste text from a source; write your own explanation and link. Leishman's book is CC BY-NC-ND: link to it by chapter title and URL, never by chapter number, and never adapt its figures.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Honesty rules: never claim the course grants a degree or credential or is equivalent to one; never present simulation as physical measurement; say what a model assumes and where it fails; no invented time estimates.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders (and `projects/p08-flow-pressure/` for the lesson 19 author).
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
- Keep scratch scripts in your own scratchpad with unique names; write code and tests first, prose second.
