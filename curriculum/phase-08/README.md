# Phase 08: Fluid Mechanics

First undergraduate fluid mechanics, at the level of a first fluids course and the fluids part of an aerospace unified-engineering sequence (the scope of MIT 2.06 or of a Çengel–Cimbala or White course), taught through aircraft and rocket examples. The phase is organised around one idea: *conservation of mass, momentum and energy, applied to a control volume or to a fluid particle, with a model for the stress*. Every lesson says three things. It names the conservation law it uses. It says whether that law is written in integral form (a balance on a control volume) or in differential form (an equation at every point of the flow). And it says what it assumes about the stress: hydrostatic (a fluid at rest, pressure only), inviscid (a moving fluid with pressure only), Newtonian laminar (pressure plus a viscous stress proportional to the rate of strain), or a turbulence model (the mean stress plus a modelled Reynolds stress). Changing the stress model is what takes you from Bernoulli to Navier–Stokes to the log law, and each model is stated with the place where it fails.

Every main result is checked by an independent route:
- mass and momentum balances that close to round-off;
- the standard atmosphere against the ISA anchor values (published data);
- Couette and Poiseuille exact solutions against a finite-difference solve;
- Colebrook by Newton's method against the Haaland explicit form;
- potential-flow superposition satisfying zero normal flow on the body, and Kutta–Joukowski against pressure integration;
- Blasius by shooting against `scipy.integrate.solve_bvp`;
- Thwaites against Blasius on a flat plate and against the Hiemenz stagnation-point solution;
- the Orr–Sommerfeld critical Reynolds number against the published value 5772.22;
- a synthetic turbulence signal whose analysis must recover its known variance and integral scale;
- the lid-driven cavity against the Ghia, Ghia and Shin (1982) Re = 100 table, labelled as published data;
- grid convergence with an observed order wherever a grid is used;
- in lesson 19, Torricelli's law and Bernoulli against your own measurement.

By the end you can compute fluid properties and a standard atmosphere, and describe a flow field and its lines. You can write mass and momentum balances for a control volume, measure airspeed with a pitot-static tube, and say when Bernoulli fails. You can find dimensionless groups and scale a model test, and solve the exact viscous flows. You can size a propellant feed line, build potential flows and the lift on a cylinder, and reason with vorticity and circulation. You can compute laminar boundary layers exactly and by an integral method, estimate transition, analyse a turbulent signal, and use the log law and skin-friction correlations. You can explain what a turbulence model closes, estimate bluff-body drag, and handle low-speed compressibility. You can write and verify a small incompressible Navier–Stokes solver, and measure a discharge coefficient and an airspeed at home and compare them with theory. You get a tested Python module for each of these. The mathematics, numerical methods and physics are not new: Phases 01–04 supplied them. What is new is fluid mechanics: continuum fluid properties, kinematics, the conservation laws in control-volume and differential form, viscous and inviscid flow, boundary layers, turbulence, an introduction to compressibility, and a first CFD code.

One lesson has a physical experiment (marked 🔬): the flow-and-pressure home lab. A learner without the kit completes it with the synthetic dataset supplied, which is labelled synthetic wherever it appears. Lesson 13 is a data-based investigation, on a seeded synthetic turbulence signal generated from a stated model spectrum (or on an authentic public dataset only where the reference map verifies one, with its source and license named). Every other lesson's part 4 is a computational experiment, such as a parameter sweep, a sensitivity study, a comparison of two methods or a grid-convergence study. A calculation or a simulation is never presented as a measurement, and a correlation is never called one.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Fluid Properties & Hydrostatics](lesson-01-fluid-properties-hydrostatics/lesson.md) | Standard atmosphere from scratch (ISA layers to 47 km by integrating the hydrostatic equation, checked against the closed forms and the anchor values; Sutherland viscosity; manometer and buoyancy) | computational |
| 02 | [Fluid Kinematics](lesson-02-fluid-kinematics/lesson.md) | Streamline, pathline and streakline plotter (unsteady 2-D field; material derivative, strain rate and rotation) | computational |
| 03 | [Control Volume: Mass & Momentum](lesson-03-control-volume-mass-momentum/lesson.md) | Jet thrust from a momentum balance (Reynolds transport theorem; pressure thrust; moving control volumes) | computational |
| 04 | [Bernoulli & Energy Equation](lesson-04-bernoulli-energy-equation/lesson.md) | Pitot-static airspeed (IAS, CAS, EAS and TAS; static-port error; the energy equation with losses and shaft work) | computational |
| 05 | [Dimensional Analysis & Similarity](lesson-05-dimensional-analysis-similarity/lesson.md) | Buckingham Pi tool (nullspace by SymPy and by SVD) and scaling a wind-tunnel model | computational |
| 06 | [Navier-Stokes Equations](lesson-06-navier-stokes-equations/lesson.md) | Derivation and exact solutions (Couette, Poiseuille, Couette–Poiseuille, Stokes's first problem; substitution and finite-difference checks) | computational |
| 07 | [Pipe Flow & Losses](lesson-07-pipe-flow-losses/lesson.md) | Propellant feed-line pressure drop (Darcy–Weisbach, Colebrook by Newton against Haaland, minor losses, a tank-to-injector budget) | computational |
| 08 | [Potential Flow](lesson-08-potential-flow/lesson.md) | Superposition of elementary flows (Rankine oval, lifting cylinder, Kutta–Joukowski by contour integration, d'Alembert's paradox) | computational |
| 09 | [Vorticity & Circulation](lesson-09-vorticity-circulation/lesson.md) | Kelvin's theorem demo (a material loop in a vortex flow, inviscid and Lamb–Oseen; Stokes's theorem; Biot–Savart; the starting vortex) | computational |
| 10 | [Boundary Layers: Laminar](lesson-10-laminar-boundary-layers/lesson.md) | Blasius solution (RK4 shooting against `solve_bvp`; δ₉₉, δ*, θ, H, c_f) | computational |
| 11 | [Integral Boundary Layer Methods](lesson-11-integral-boundary-layer-methods/lesson.md) | Thwaites method (against Blasius and Hiemenz; Falkner–Skan wedges; separation on a cylinder and in a decelerating flow) | computational |
| 12 | [Transition](lesson-12-transition/lesson.md) | e^N method intro (Orr–Sommerfeld by Chebyshev collocation against Re = 5772.22; an e^N envelope from an approximate fit; N = 9) | computational |
| 13 | [Turbulence Fundamentals](lesson-13-turbulence-fundamentals/lesson.md) | Turbulence statistics from data (Reynolds decomposition, moments, autocorrelation and integral scale, Taylor's hypothesis, spectra, Kolmogorov scales) | data-based investigation, on a seeded synthetic signal labelled synthetic (or a verified public dataset) |
| 14 | [Turbulent Boundary Layers](lesson-14-turbulent-boundary-layers/lesson.md) | Log law and skin friction (law of the wall, Spalding's profile, Clauser fit, power laws, flat-plate correlations, the mixed plate) | computational, with a public or synthetic mean profile |
| 15 | [Turbulence Modeling](lesson-15-turbulence-modeling/lesson.md) | Mixing-length and k-ε overview (van Driest mixing length across a channel; k–ε constants and the log-layer equilibrium check) | computational |
| 16 | [Flow Separation & Drag](lesson-16-flow-separation-drag/lesson.md) | Drag of bluff bodies (sphere c_d correlations, drag crisis, Strouhal shedding, terminal velocity, streamlining) | computational |
| 17 | [Intro to Compressible Flow](lesson-17-intro-compressible-flow/lesson.md) | Speed of sound and Mach number (isentropic relations; where Bernoulli airspeed fails; compressible pitot below M = 1 and the Rayleigh pitot formula above it) | computational |
| 18 | [Incompressible CFD from Scratch](lesson-18-incompressible-cfd-lid-driven-cavity/lesson.md) | 2-D lid-driven cavity solver (Re = 100 against the Ghia et al. 1982 table, labelled published; grid convergence with observed order; mass conservation) | computational |
| 19 | [Home Lab: Flow & Pressure Measurement](lesson-19-home-lab-flow-pressure-measurement/lesson.md) 🔬 | Low-energy flow rig: a draining bottle with a discharge-coefficient fit, and a water U-tube manometer read against a guarded fan or a hair dryer on cool; Reynolds-number reasoning and calibrated pressure | physical experiment ([project](../../projects/p08-flow-pressure/)) |

Lesson ids are `P08.01` to `P08.19`. Each lesson's metadata table gives its exact prerequisites, and every lesson's metadata table reads:

| **Time** | unmeasured |
|---|---|

The Time row stays `unmeasured` until pilots measure it. No time estimate is given here for the same reason.

## Prerequisites

Phases 00–05, plus two lessons of Phase 07, which is being written in parallel. `P07.01` (Systems, Properties & State) supplies the ideal-gas law p = ρRT and the idea of a state fixed by two properties, used from lesson 01 on. `P07.05` (Isentropic Processes & Efficiencies) supplies the isentropic relation p/ρ^γ = constant for a reversible adiabatic ideal-gas process, used in lesson 17. Each lesson restates what it uses from them in two or three sentences, so you can start before Phase 07 is published. From the earlier phases:
- `P00.03` (a vectorised ISA table) and `P00.05` (units);
- `P00.13` (safety for home experiments), `P00.14` (measurement, calibration and uncertainty) and `P00.15` (the first predict–measure–compare experiment), for lesson 19;
- `P01.08` (partial derivatives), `P01.11` (vector fields, line and surface integrals) and `P01.12` (divergence and Stokes theorems), for kinematics, control volumes and circulation;
- `P01.14` (rank and nullspace) and `P01.18` (SVD), for the Buckingham Pi tool; `P01.16` (eigenvalues), for the Orr–Sommerfeld problem; `P01.17` (least squares), for the fits;
- `P02.02` (second-order ODEs), `P02.11` (PDEs), `P02.12` (complex variables, for potential flow), `P02.10` (FFT) and `P02.19` (random processes and power spectra, for turbulence);
- `P03.02` (root finding), `P03.03` (interpolation), `P03.04` (quadrature), `P03.05` (RK4, for Blasius), `P03.06` (adaptive and stiff integrators), `P03.08`–`P03.09` (direct and iterative linear solvers), `P03.10` (Newton for systems), `P03.11` (finite differences, stability, observed order), `P03.12` (finite volume methods, for the cavity) and `P03.18` (verification and validation);
- `P04.01`–`P04.03` (Newton's laws, energy, momentum) and `P04.07` (kinetic theory, for viscosity and the continuum).

Skill checks before starting:
- you can integrate a separable ODE with a boundary condition;
- you can find the nullspace of a 3 × 6 matrix by hand;
- you can solve a scalar nonlinear equation by Newton's method;
- you can write the real and imaginary parts of a complex function;
- you can propagate an uncertainty through a square root (the airspeed from a manometer goes as √Δh).

## The project kit

The kit lives in `projects/p08-flow-pressure/` and contains:
- `bom.yaml`: three cost tiers (`minimum_cost`, `recommended`, `enhanced`) with dated price ranges, reusing the kitchen scale and phone-video habits of the `p00-calibration` and `p00-first-experiment` kits and of the Phase 04 kits;
- `safety.md`: the full safety box;
- `procedure.md`, `report-template.md` and `data/README.md` (the raw-data formats);
- a seeded synthetic dataset under `data/synthetic/`, labelled synthetic in its file names and headers.

Read `safety.md` before buying anything and again before every session. The kit's `safety.md` is authoritative. The limits below summarise the hard rules that every version of this kit must keep, and nothing in the kit may relax them. The prediction, with its uncertainty, is written down and dated before measuring, in both experiments. Raw data are never edited, and the report states which numbers were measured, which were computed, which are manufacturer's figures and which are synthetic.

- **`p08-flow-pressure`** ([folder](../../projects/p08-flow-pressure/), lesson 19). Experiment 1 drains a plastic bottle of room-temperature water through a small round hole near its base into a basin or sink. The level is read against a taped-on ruler from a phone video. You predict the draining time from Torricelli's law with a discharge coefficient, fit the coefficient from h(t) (√h is linear in t), and compare it with the typical range for a sharp-edged orifice. Experiment 2 uses a clear U-tube part-filled with water and connected to a straw pitot probe. It is read in the stream of a small guarded fan or of a hair dryer on its cool setting, giving Δp = ρ_w gΔh and V = √(2Δp/ρ_air). You state the reading's resolution limit, decide whether an inclined manometer improves it, and compute the Reynolds numbers of the probe and of the fan. Hard limits:
  - water at room temperature, no more than 2 L in total;
  - fans guarded (finger-proof grille) and drawing at most 30 W;
  - the hair dryer used on its cool setting only and kept away from water (separate tables; the manometer's water sealed or below the dryer; never in the same basin);
  - no compressed gas and no pressurised containers (no squeezing or pressurising bottles to drive flow);
  - no mains appliance near spilled water, and spills mopped before anything is switched on;
  - stop at any spill near a plug or appliance, or at any appliance that sparks, smells or overheats.

Nothing in this phase involves heated water, combustion, propellants, pressure vessels, compressed gas, open rotors or flight. The propellant feed line of lesson 07 and the rocket nozzle of lesson 03 are computed, never built. **A report without a completed safety box is not graded.**

## Data honesty

Synthetic data are labelled synthetic everywhere they appear: in the file name or header, in the docstring that generates them, in the lesson text, in figure titles and in worked solutions. The turbulence signal of lesson 13 and any synthetic mean profile in lesson 14 are seeded and generated from a stated model, so the analysis can be checked against the known answer. Published values are labelled "published data" with their citation: the ISA anchor values, the Blasius constants, the Orr–Sommerfeld critical Reynolds number, and the Ghia, Ghia and Shin (1982) cavity table, which is itself a numerical computation, not a measurement and not an exact solution. Drag, friction and transition formulas are labelled "correlation". Propellant properties are typical values, stated with their source and labelled typical. A simulation is never presented as a measurement.

## Running the code

From the repo root, `pytest curriculum/phase-08/lesson-NN-<slug>` tests your `starter/main.py`. `AEFS_IMPL=solution pytest curriculum/phase-08/lesson-NN-<slug>` runs the same tests against the reference solution in `solutions/phase-08/lesson-NN-<slug>/solution.py`. Tests check properties rather than single numbers:
- conservation residuals (mass and momentum balances to round-off; ∇·u = 0 to discretisation order where it is not automatic);
- exact solutions substituted back into the governing equations and boundary conditions;
- limiting cases (zero velocity, zero viscosity, infinite Reynolds number, an isothermal layer);
- agreement between two methods (shooting against `solve_bvp`, Colebrook by Newton against Haaland, SymPy against SVD nullspace);
- published anchor values, labelled as published;
- observed order of convergence from a grid sequence;
- dimensional sanity and invariance (results unchanged by rescaling units).

The code uses NumPy, SciPy, Matplotlib and SymPy only; no CFD, fluids or atmosphere package is used in tested code. Anything random is seeded. Tests need no network and write files only to a temporary directory, and the cavity tests use small grids so that each lesson's tests stay fast. `python scripts/check_lessons.py` checks every lesson's structure.

## The exam and the mini-project

**Phase exam and mini-project:** [assessments/phase-08/](../../assessments/phase-08/). The exam has four parts:
- **Part A, concepts:** the continuum and the stress models, inviscid flow and circulation, boundary layers and transition, turbulence, similarity and compressibility.
- **Part B, analysis by hand** over lessons 01 to 17: hydrostatics and the ISA, kinematics, control volumes, pitot-static airspeed, Buckingham Pi, exact solutions, a feed line, potential flow, the starting vortex, Blasius, Thwaites, transition, turbulence statistics, the log law, turbulence modelling, drag and compressibility.
- **Part C, computation:** a tested Blasius solver by RK4 shooting, checked against `solve_bvp`, the published f''(0) and the momentum integral.
- **Part D, engineering interpretation:** its first question asks you to interpret your own `p08-flow-pressure` measurement, with a seeded synthetic dataset for learners without a kit. Its second reviews a lid-driven-cavity verification report against the Ghia et al. table (published data) and tests the kit's safety limits.

The exam has not been piloted, so no time limit is suggested. Choose one in advance and record it. Results are reported per competency: (1) reasoning, (2) analysis, (3) computation, (4) interpretation. The worked solutions and pass thresholds are in [solutions/phase-08/exam/](../../solutions/phase-08/exam/worked-solutions.md). The mini-project is the flow-and-pressure report from lesson 19, graded with [the rubric](../../assessments/phase-08/grading-rubric.md). Seeded unseen cases for the computational lessons are planned in [parameterized-cases/](../../assessments/phase-08/parameterized-cases/README.md). Outcomes are recorded as *practice completed*, *self-assessed mastery* or *externally reviewed mastery*. None of these is a credential, and nothing in this phase qualifies anyone to calibrate a real air-data system, size a real propellant feed system or approve a CFD code for design.

**References:** [references/phase-08.md](../../references/phase-08.md). Every required reading is free; each lesson names the exact sections it uses. The curriculum links to every source and writes its own explanations.
