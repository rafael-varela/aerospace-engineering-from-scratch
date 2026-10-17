# Phase 10: Aerodynamics

Incompressible aerodynamics at the level of an undergraduate aerodynamics course and the start of a graduate one: the scope of MIT 16.100 and the fluids part of Unified Engineering. It is taught through the design of real aircraft, drones and rotors, and organised around one idea: *a hierarchy of models, each checked against the one above and below it*.

| Family | The hierarchy |
|---|---|
| 2-D sections | thin-airfoil theory, then panel methods, then boundary layers |
| 3-D wings | horseshoe vortices, then lifting line, then vortex lattices |
| Rotors | momentum theory, then blade element, then blade element momentum (BEM) |
| Time | quasi-steady theory, then Theodorsen's frequency response and Wagner's step response |

Every lesson says where its model sits in that hierarchy, what it neglects, and checks its main result by an independent route:

- thin-airfoil theory against the panel method as the thickness goes to zero;
- the panel method against the exact Joukowski or Karman–Trefftz solution, with the observed order of convergence;
- lifting line against the elliptic-wing closed form;
- the vortex lattice against lifting line for unswept, high-aspect-ratio wings, and against its own lattice convergence;
- BEM against actuator-disk momentum theory;
- Theodorsen against its quasi-steady limit, and Wagner against Theodorsen through the Fourier transform;
- in lesson 23, all of them against your own measurement on a bench rig.

Every grid, panel or collocation study reports its observed order.

By the end you can:

- compute and transform aerodynamic coefficients, and generate airfoils;
- explain circulation and the Kutta condition, and solve thin-airfoil and panel problems;
- add a boundary layer and predict section drag;
- design and select airfoils, and quantify high-lift devices;
- compute induced drag, span loads and stability derivatives with lifting line and a vortex lattice, and optimise a planform;
- build a drag polar, and reason about low Reynolds number and stall;
- estimate the forces on bodies, and analyse propellers and helicopter rotors;
- compute unsteady lift and estimate ground effect;
- validate your own codes against free professional tools and authentic wind-tunnel data;
- measure lift and drag at home and explain why your numbers differ from the published ones.

You get a tested Python module for each of these.

The mathematics, numerics and fluid mechanics are not new: Phases 01–03 and the fluid mechanics of Phase 08 supply them. What is new is aerodynamics: the forces on lifting surfaces and bodies, and the models that predict them. The reference vehicles recur throughout:

- **The reference trainer:** 750 kg, a 12.0 m² rectangular wing of 10.0 m span and 1.2 m chord with a NACA 2412 section, cruise 55 m/s and approach 30 m/s at sea level, a 7.0 m fuselage of 1.1 m maximum diameter, and a 2.0 m² tail 4.5 m behind the wing's quarter chord.
- **The reference small UAV:** 1.5 m span, 0.20 m chord, 1.2 kg, cruising at 12 m/s on a 10 × 4.5 in two-blade propeller.
- **The reference rotor:** a 1100 kg light helicopter with a two-bladed rotor of 5.0 m radius and 0.25 m chord, NACA 0012 blades and a 210 m/s tip speed.

All of them fly in ISA sea-level air unless a lesson says otherwise.

One lesson has a physical experiment (marked 🔬): the airfoil bench rig. A learner without the kit completes it with the seeded synthetic dataset supplied, which is labelled synthetic wherever it appears. A calculation or a simulation is never presented as a measurement. Lessons 13 and 22 are data-based investigations of authentic NACA and NASA wind-tunnel data. Every other lesson's part 4 is a computational experiment: a convergence study, a method comparison, a parameter sweep, a sensitivity study, or a comparison with an authentic public dataset.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Forces, Moments & Coefficients](lesson-01-forces-moments-coefficients/lesson.md) | Coefficient calculator from balance data (balance axes to wind axes, moment transfer, dynamic pressure, cl, cd, cm; the effect of a moment-reference shift) | computational, on a synthetic balance dataset |
| 02 | [Airfoil Geometry](lesson-02-airfoil-geometry/lesson.md) | NACA 4- and 5-digit generator (thickness and camber lines, cosine spacing, closed trailing edge, area, thickness and camber checks) | computational |
| 03 | [The Kutta Condition & Circulation](lesson-03-kutta-condition-circulation/lesson.md) | Kutta–Joukowski lift (cylinder with circulation, Joukowski and Karman–Trefftz maps, the Kutta condition fixing Γ, exact surface Cp) | computational |
| 04 | [Thin Airfoil Theory](lesson-04-thin-airfoil-theory/lesson.md) | Thin-airfoil solver for any camber line (A₀, A₁, A₂ by quadrature; α_L0, cl, c_m,c/4, centre of pressure; a flap as a camber line) | computational |
| 05 | [Vortex Panel Method](lesson-05-vortex-panel-method/lesson.md) | 2-D linear-strength vortex panel method from scratch (Kutta condition, cl, Cp, cm) | computational: panel convergence with observed order, against the exact Karman–Trefftz solution |
| 06 | [Viscous-Inviscid Coupling](lesson-06-viscous-inviscid-coupling/lesson.md) | Panel method + Thwaites laminar, Michel transition, Head turbulent boundary layer, Squire–Young drag | computational, compared with authentic NACA 0012 data |
| 07 | [Airfoil Design & Selection](lesson-07-airfoil-design-selection/lesson.md) | Inverse design (target load to camber line by thin-airfoil inversion) and a selection matrix over candidate sections | computational |
| 08 | [High-Lift Devices](lesson-08-high-lift-devices/lesson.md) | Flap and slat effects (flap effectiveness, Δα_L0 and Δcm; a multi-element panel model; a Δcl_max estimate; the landing-speed payoff) | computational |
| 09 | [Finite Wings & Induced Drag](lesson-09-finite-wings-induced-drag/lesson.md) | Downwash and induced drag (horseshoe vortex, Biot–Savart, elliptic loading, span efficiency, the drag polar and best L/D) | computational |
| 10 | [Prandtl Lifting-Line Theory](lesson-10-prandtl-lifting-line/lesson.md) | Lifting-line solver (the monoplane equation by collocation; arbitrary taper and twist; δ and τ) | computational: collocation convergence and the elliptic closed form |
| 11 | [Vortex Lattice Method](lesson-11-vortex-lattice-method/lesson.md) | VLM from scratch (horseshoe vortices on a planform; CL, CDi in the Trefftz plane, spanwise load, Cm, neutral point) | computational: lattice convergence; VLM against lifting line for unswept high-AR wings |
| 12 | [Wing Planform Effects](lesson-12-wing-planform-effects/lesson.md) | Taper, sweep and twist optimizer (span efficiency, local cl and stall-onset station, washout to minimise induced drag under constraints) | computational |
| 13 | [Drag Build-Up](lesson-13-drag-build-up/lesson.md) | Component drag estimation (skin friction, form factors, interference, wetted areas, C_D0 of the reference trainer, the polar) | data-based investigation against authentic NACA 0012 drag data |
| 14 | [Low Reynolds Number Aerodynamics](lesson-14-low-reynolds-number-aerodynamics/lesson.md) | Drone and glider airfoils (laminar separation bubbles, the Re-dependence of cd and cl_max, Thwaites separation, transition) | computational |
| 15 | [Stall & Post-Stall](lesson-15-stall-post-stall/lesson.md) | Stall progression model (Kirchhoff/Helmholtz separation-point model, nonlinear lifting line, flat-plate post-stall limit, spanwise stall pattern) | computational, fitted to authentic NACA 0012 data |
| 16 | [Fuselage & Body Aerodynamics](lesson-16-fuselage-body-aerodynamics/lesson.md) | Slender-body theory (Munk moment, apparent mass, wetted area and fineness ratio, body drag, the body's contribution to Cm_α) | computational |
| 17 | [Propeller & Rotor Theory](lesson-17-propeller-rotor-theory/lesson.md) | Blade element momentum code (Prandtl tip loss, C_T, C_P, C_Q, efficiency against advance ratio) | computational: BEM against actuator-disk momentum theory |
| 18 | [Helicopter Aerodynamics](lesson-18-helicopter-aerodynamics/lesson.md) | Hover and forward-flight performance (momentum theory, figure of merit, blade-element hover, Glauert inflow, power against airspeed, autorotation) | computational |
| 19 | [Unsteady Aerodynamics](lesson-19-unsteady-aerodynamics/lesson.md) | Theodorsen and Wagner functions (C(k) from Bessel functions, indicial response, Duhamel integral, lift loops) | computational: Theodorsen against its quasi-steady limit, Wagner against Theodorsen via the Fourier transform |
| 20 | [Ground Effect & Interference](lesson-20-ground-effect-interference/lesson.md) | Wing-in-ground-effect estimate (image vortex system in the lattice, induced-drag reduction, lift change, wing–tail downwash) | computational |
| 21 | [Using XFOIL, AVL & OpenVSP](lesson-21-xfoil-avl-openvsp/lesson.md) | Validate your codes against the tools (optional practice with free tools; nothing in the tests imports or runs them) | computational: your codes against published NACA data and, where you have it, labelled tool output |
| 22 | [Aerodynamic Validation Project](lesson-22-aerodynamic-validation-project/lesson.md) | Compare your VLM to wind-tunnel data (code verification, solution verification, validation against authentic finite-wing and section data, with uncertainty) | data-based investigation |
| 23 | [Home Lab: Airfoil Lift & Drag on a Bench Rig](lesson-23-home-lab-airfoil-bench/lesson.md) 🔬 | Guarded fan with a straw-honeycomb straightener, foam or card airfoil on a kitchen-scale lift balance and a pivot-arm drag balance; CL and CD against α; theory against computation against measurement | physical experiment ([project](../../projects/p10-airfoil-bench/)) |

Lesson ids are `P10.01` to `P10.23`. Each lesson's metadata table gives its exact prerequisites, and every lesson's metadata table reads:

| **Time** | unmeasured |
|---|---|

The Time row stays `unmeasured` until pilots measure it. No time estimate is given here for the same reason.

## Prerequisites

Phases 00–05 and, in the roadmap's order, 06–09. In particular:

- `P00.13` (safety for home experiments), `P00.14` (measurement, calibration and uncertainty) and `P00.15` (the first predict-measure-compare experiment), for lessons 01 and 23;
- `P01.05` (Taylor series), `P01.09` (gradients and optimisation), `P01.11` (vector fields, line integrals and circulation), `P01.12` (Stokes' and Green's theorems) and `P01.13` (linear systems);
- `P02.08` (the Laplace transform) for unsteady lift, `P02.09` (Fourier series) for thin-airfoil theory and lifting line, `P02.12` (complex variables and conformal maps) for the Joukowski airfoil, and `P02.18` (regression and curve fitting);
- `P03.05` (ODE integrators), `P03.08` (direct linear solvers), `P03.11` (finite differences and observed order), `P03.14` (optimisation), `P03.16` (multidisciplinary design optimisation) and `P03.18` (verification and validation);
- `P04.04` (rotation) for balance axes;
- `P05.12` (beam deflection) where a wing bends.

The roadmap also names Phase 08 (fluid mechanics) lessons as prerequisites: `P08.03` (control-volume momentum), `P08.05` (dimensional analysis), `P08.08` (potential flow), `P08.11` (integral boundary layers and Thwaites' method), `P08.12` (transition), `P08.14` (turbulent boundary layers) and `P08.19` (the flow and pressure home lab, whose manometer method lesson 23 uses). **Phase 08 is not written yet.** Its lessons are therefore referred to by id only, with no link. Each Phase 10 lesson carries a short, self-contained recap (a few lines, with the formula used) of any Phase 08 result it needs, so the phase works before Phase 08 exists.

Skill checks before starting: you can resolve a force into components along two rotated axes; integrate a polynomial and a Fourier cosine coefficient by hand; solve a small linear system; differentiate a complex function and evaluate a line integral round a circle; state the momentum balance for a control volume; and fit a straight line with the standard error of its slope.

## Shared data

`curriculum/phase-10/data/` holds the authentic, public-domain wind-tunnel data that several lessons read (lessons 03, 06, 13, 14, 15, 21, 22 and 23 among them). It is read-only: never edit a value; apply any correction in code and say so. Its `README.md` gives each file's full citation, NTRS id and URL, licence, facility, Reynolds and Mach numbers, what was measured and how, and how it was transcribed and checked.

- `ladson_naca0012_M015.csv`: Ladson, NASA TM-4074 (1988), NACA 0012 at M = 0.15 in the Langley Low-Turbulence Pressure Tunnel. It holds Table I (free transition, Re = 2.0–18.9 × 10⁶) and Table VII (a wraparound No. 60 grit strip, Re = 2.0–12.1 × 10⁶): cl and cm from surface pressures, cd from a wake rake, wall-corrected by the report. There are 161 rows, transcribed from the rendered page images and checked against the printed l/d.
- `finite_wing_cla.csv`: measured lift-curve slopes of twelve wing planforms (fifteen rows) from NACA TN-1669, TN-1671 and TN-1422 (aspect ratio 1.34–9, sweep 0–60°, taper 0.25–1). The TN-1669 and TN-1671 slopes were digitised from the reports' figures, and their reading uncertainties are stated.

Synthetic datasets live with the lesson or project that generates them (for example `lesson-01-forces-moments-coefficients/data/` and `projects/p10-airfoil-bench/data/synthetic/`). Their file names start with `SYNTHETIC_` and their first line says which seeded function made them.

## The project kit

The kit lives in `projects/p10-airfoil-bench/` with these files:

- `bom.yaml`: three cost tiers, `minimum_cost`, `recommended` and `enhanced`, with dated price ranges, reusing items from `p00-calibration`, `p00-first-experiment`, `p04-drag-drop`, `p04-electronics-bench` and `p05-cantilever`;
- `safety.md`: the full safety box;
- `procedure.md`, `report-template.md` and `data/README.md`: the raw-data format;
- `data/synthetic/`: the seeded synthetic dataset for learners without the kit.

Read `safety.md` before buying anything and again before every session. The kit's `safety.md` is authoritative. The limits below summarise the hard rules that every version of this kit must keep, and nothing in the kit may relax them. The prediction, with its uncertainty and the predicted airspeed, is written down and dated before the first reading. Raw data are never edited. The report states which numbers were measured, which were computed, which are authentic published data and which are synthetic.

- **`p10-airfoil-bench`** ([folder](../../projects/p10-airfoil-bench/), lesson 23). A guarded desk or box fan with a straw-honeycomb flow straightener and a short settling length; a foam (XPS or EPS insulation board) or folded-card airfoil of known chord and span between end plates; a kitchen scale under a sting for lift, and a pivot arm (bell-crank) on a second scale, or a calibrated spring, for drag; the angle of attack set with a protractor template; the airspeed from a pitot tube and an inclined water manometer, or a vane anemometer. The dynamic pressure is tiny (q ≈ 10 Pa at 4 m/s), and that sets the uncertainty. You sweep α, compute CL and CD with uncertainties, and compare them with thin-airfoil theory, your own panel code (with and without a boundary-layer estimate) and Ladson's NACA data. Then you explain the Reynolds-number gap: Re ≈ 3 × 10⁴ on the bench against 2–9 × 10⁶ in the NACA tunnels means laminar separation, a lower lift-curve slope, much higher cd, and an earlier and softer stall. Hard limits:
  - the fan is guarded on both sides, with front and back grilles intact and finger-proof;
  - the fan is rated at most 50 W, and stands on a stable surface where it cannot tip or walk;
  - no propellers, drones or unguarded blades, ever;
  - no cutting tools near the running fan;
  - hot-wire foam cutting only with adult supervision and ventilation; otherwise a craft knife with the fan unplugged;
  - unplug the fan before adjusting the model, and before reaching anywhere near the fan;
  - the cable routed so it cannot be tripped over or pulled into the fan;
  - loose paper, hair, scarves and string kept away from the intake;
  - the model and balance secured so nothing can be blown into the fan;
  - stop if the fan rattles, the guard is damaged, the motor smells hot, or anything touches a grille.

Lessons 17 and 18 analyse propellers and rotors on the computer only: no propeller or rotor is run at home, at any size. **A report without a completed safety box is not graded.**

## XFOIL, AVL and OpenVSP

Lesson 21 gives the exact commands and settings for three free professional tools, for learners who choose to install them: XFOIL and AVL (Mark Drela, MIT; GNU GPL) and OpenVSP with VSPAERO (NASA Open Source Agreement). It treats them as cross-checks of your own codes, never as substitutes for them. The tools were not available where the lesson was written, so no tool output ships with the phase. A polar is accepted only with a header that records the tool, its version and its settings. Nothing in the tests imports or runs a tool, so every test runs without them installed.

## Running the code

From the repo root, `pytest curriculum/phase-10/lesson-NN-<slug>` tests your `starter/main.py`. `AEFS_IMPL=solution pytest curriculum/phase-10/lesson-NN-<slug>` runs the same tests against the reference solution in `solutions/phase-10/lesson-NN-<slug>/solution.py`. Tests check properties rather than single numbers:

- residuals of the boundary condition and the Kutta condition;
- conservation (momentum and energy in the actuator disk);
- symmetry (a symmetric airfoil at α = 0 gives cl = 0; a symmetric wing gives a symmetric span load) and linearity in α;
- limiting cases (thickness → 0, AR → ∞, k → 0, h/b → ∞);
- agreement between two methods, and a textbook closed form;
- the observed order of every convergence study, and dimensional sanity.

The code uses NumPy, SciPy and Matplotlib only; no aerodynamics library is used in tested code. Lessons that need code from an earlier lesson (a NACA generator, a panel solver, a horseshoe-vortex function) carry a compact copy as a given helper, so no lesson imports another. Anything random is seeded. Tests need no network and write files only to a temporary directory. `python scripts/check_lessons.py` checks every lesson's structure.

## The exam and the mini-project

**Phase exam and mini-project:** [assessments/phase-10/](../../assessments/phase-10/). The exam has four parts:

- **Part A, concepts:** coefficients, similarity and the hierarchy of models; circulation, the Kutta condition and induced drag; boundary layers, Reynolds number and stall; rotors, bodies, unsteady flow and ground effect.
- **Part B, analysis by hand:** over every lesson from 01 to 22.
- **Part C, computation:** a tested lifting-line solver, checked against the elliptic closed form and its observed orders of convergence.
- **Part D, engineering interpretation:** its first question asks you to interpret your own `p10-airfoil-bench` measurement, and a seeded synthetic dataset with a deliberate flaw is supplied for learners without a kit. Its second question reviews a synthetic aerodynamic data sheet for out-of-bounds inputs and 2-D-to-3-D confusions, and tests the bench rig's safety limits.

The exam has not been piloted, so no time limit is suggested; choose one in advance and record it. Results are reported per competency: (1) reasoning, (2) analysis, (3) computation, (4) interpretation. The worked solutions and pass thresholds are in [solutions/phase-10/exam/](../../solutions/phase-10/exam/worked-solutions.md). The mini-project is the bench-rig report from lesson 23, graded with [the rubric](../../assessments/phase-10/grading-rubric.md), which also gives the pass criteria for every lesson. Seeded unseen cases for the computational lessons are planned in [parameterized-cases/](../../assessments/phase-10/parameterized-cases/README.md). Outcomes are recorded as *practice completed*, *self-assessed mastery* or *externally reviewed mastery*. None of these is a credential, and nothing in this phase qualifies anyone to produce or approve the aerodynamic data, performance figures or handling assessments of a real aircraft.

**References:** [references/phase-10.md](../../references/phase-10.md). Every required reading is free. The core is:

- the MIT OpenCourseWare aerodynamics and fluids notes (`@mit-ocw-16-100-aerodynamics`, `@mit-ocw-16-unified-fluids`);
- Leishman's *Introduction to Aerospace Flight Vehicles* (`@leishman-intro-aerospace-flight-vehicles`, CC BY-NC-ND, linked by chapter title and never adapted);
- the NASA Glenn *Beginner's Guide to Aeronautics* (`@nasa-grc-beginners-guide-aeronautics`);
- NACA and NASA reports, among them the airfoil summary (`@ntrs-naca-tr-824-summary-of-airfoil-data`), Ladson's NACA 0012 data (`@ntrs-nasa-tm-4074-ladson-naca0012`), the panel-method primer (`@ntrs-nasa-tp-2995-panel-methods`), the trailing-edge high-lift summary (`@ntrs-naca-tr-938-trailing-edge-high-lift`), the stall-types note (`@ntrs-naca-tn-2502-airfoil-stall-types`), the Eppler 387 low-Re data (`@ntrs-nasa-tm-4062-eppler-387-low-re`), Stepniewski's rotary-wing aerodynamics (`@ntrs-nasa-cr-3082-rotary-wing-aerodynamics`), and the finite-wing reports (`@ntrs-naca-tn-1669-untapered-wings-aspect-ratio-sweep`, `@ntrs-naca-tn-1671-swept-wings-taper-ratio`, `@ntrs-naca-tn-1422-sivells-three-wings`).

The free tools' documentation (`@drela-xfoil`, `@drela-avl`, `@openvsp-docs`) and the MIT viscous-flow notes (`@mit-ocw-16-13-aerodynamics-viscous-fluids`) are supplementary. Paid textbooks appear only as optional further reading in the reference map, with what each adds. The curriculum links to every source and writes its own explanations.
