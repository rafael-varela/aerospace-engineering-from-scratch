# Author brief: Phase 10 lessons

This brief is what every lesson author in Phase 10 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-10.md`, `references/sources.yaml`, the Phase 10 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/`, `solutions/_template/` and `projects/_template/` folders. **Phases 00–05 are the quality bar:** before writing, read `curriculum/phase-05/lesson-12-beam-deflection/` and `curriculum/phase-03/lesson-16-multidisciplinary-design-optimization/` end to end (lesson, problems, starter, tests, and the matching `solutions/` folder) and match their depth, tone and structure. For the lesson with a physical experiment (lesson 23), also read `curriculum/phase-05/lesson-18-home-lab-cantilever-deflection/` and every file in `projects/p05-cantilever/`; that is the bar for the experiment, the bill of materials, the safety box and the synthetic dataset.

## Where files go

```
curriculum/phase-10/lesson-NN-<slug>/lesson.md
curriculum/phase-10/lesson-NN-<slug>/problems.md
curriculum/phase-10/lesson-NN-<slug>/starter/main.py
curriculum/phase-10/lesson-NN-<slug>/tests/test_main.py
curriculum/phase-10/lesson-NN-<slug>/data/...            <- optional: small CSVs the lesson reads (authentic data with provenance, or SYNTHETIC-labelled)
solutions/phase-10/lesson-NN-<slug>/hints.md
solutions/phase-10/lesson-NN-<slug>/method.md
solutions/phase-10/lesson-NN-<slug>/worked-solutions.md
solutions/phase-10/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
curriculum/phase-10/data/                                 <- shared authentic datasets (read-only for lesson authors)
projects/p10-airfoil-bench/bom.yaml, safety.md, procedure.md, report-template.md, data/README.md   <- lesson 23 only
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build | mode of part 4 |
|---|---|---|---|---|
| 01 | forces-moments-coefficients | Forces, Moments & Coefficients | Coefficient calculator from balance data (balance axes → wind axes, moment transfer, dynamic pressure, CL, CD, Cm; the effect of a moment-reference shift) | computational experiment on a synthetic balance dataset |
| 02 | airfoil-geometry | Airfoil Geometry | NACA 4- and 5-digit generator (thickness and camber lines, cosine spacing, closed trailing edge option, area, thickness and camber checks) | computational experiment |
| 03 | kutta-condition-circulation | The Kutta Condition & Circulation | Kutta–Joukowski lift (flow past a cylinder with circulation, Joukowski and Karman–Trefftz maps, the Kutta condition fixing Γ, exact surface Cp) | computational experiment |
| 04 | thin-airfoil-theory | Thin Airfoil Theory | Thin-airfoil solver for any camber line (Fourier coefficients A₀, A₁, A₂ by quadrature; α_L0, cl, cm_c/4, centre of pressure; flap as a camber line) | computational experiment |
| 05 | vortex-panel-method | Vortex Panel Method | 2-D linear-strength vortex panel method from scratch (Kutta condition, cl, Cp, cm) | computational experiment: panel convergence with observed order, against the exact Karman–Trefftz solution |
| 06 | viscous-inviscid-coupling | Viscous-Inviscid Coupling | Panel method + Thwaites laminar, Michel transition, Head turbulent boundary layer, Squire–Young drag (one-way coupling, with an optional displacement-thickness iteration) | computational experiment, compared with the NACA 0012 data of `curriculum/phase-10/data/` |
| 07 | airfoil-design-selection | Airfoil Design & Selection | Inverse design exercise (target Cp or target camber-line load → camber line by thin-airfoil inversion; a selection matrix over candidate sections for a mission) | computational experiment |
| 08 | high-lift-devices | High-Lift Devices | Flap and slat effects (thin-airfoil flap effectiveness, Δα_L0 and Δcm; a two-element vortex-panel model of a slat or slotted flap; Δcl_max estimate; landing-speed payoff) | computational experiment |
| 09 | finite-wings-induced-drag | Finite Wings & Induced Drag | Downwash and induced drag (horseshoe vortex, Biot–Savart, elliptic loading, span efficiency, CDi = CL²/(π e AR), the drag polar and best L/D) | computational experiment |
| 10 | prandtl-lifting-line | Prandtl Lifting-Line Theory | Lifting-line solver (Fourier/monoplane equation by collocation; arbitrary taper and twist; δ and τ factors) | computational experiment: collocation convergence and the elliptic closed form |
| 11 | vortex-lattice-method | Vortex Lattice Method | VLM from scratch (horseshoe vortices on a planform; CL, CDi by Trefftz plane, spanwise load, Cm, neutral point) | computational experiment: lattice convergence; VLM against lifting line for unswept high-AR wings |
| 12 | wing-planform-effects | Wing Planform Effects | Taper, sweep, twist optimizer (span efficiency, local cl distribution and stall-onset location, washout to minimise induced drag subject to constraints) | computational experiment |
| 13 | drag-build-up | Drag Build-Up | Component drag estimation (flat-plate skin friction laminar/turbulent/transitional, form factors, interference, wetted areas, CD0 of the reference aircraft; induced drag; polar) | data-based investigation: the 2-D skin-friction × form-factor estimate against the authentic NACA 0012 drag data |
| 14 | low-reynolds-number-aerodynamics | Low Reynolds Number Aerodynamics | Drone and glider airfoils (laminar separation bubbles, the Re-dependence of cd and cl_max, Thwaites separation criterion, transition, the scaling of Re for small UAVs) | computational experiment |
| 15 | stall-post-stall | Stall & Post-Stall | Stall progression model (Kirchhoff/Helmholtz separation-point model fitted to section data, nonlinear lifting line with a post-stall section curve, flat-plate post-stall limit, stall pattern along the span) | computational experiment fitted to the authentic NACA 0012 data |
| 16 | fuselage-body-aerodynamics | Fuselage & Body Aerodynamics | Slender body theory (Munk moment, normal force CNα = 2 per base area for a body ending in a full section, apparent mass, wetted area and fineness ratio, body drag, a body's contribution to Cmα) | computational experiment |
| 17 | propeller-rotor-theory | Propeller & Rotor Theory | Blade element momentum code (Glauert/Prandtl tip loss, CT, CP, CQ, efficiency vs advance ratio) | computational experiment: BEM against actuator-disk momentum theory |
| 18 | helicopter-aerodynamics | Helicopter Aerodynamics | Hover and forward-flight performance (momentum theory, figure of merit, blade-element hover with uniform inflow, Glauert forward-flight inflow, power vs airspeed, autorotation estimate) | computational experiment |
| 19 | unsteady-aerodynamics | Unsteady Aerodynamics | Theodorsen and Wagner functions (C(k) from Bessel functions, Jones/R.T. approximations, indicial response, Duhamel integral, pitching-plunging lift loops) | computational experiment: Theodorsen against its quasi-steady limit, Wagner against Theodorsen via Fourier transform |
| 20 | ground-effect-interference | Ground Effect & Interference | Wing-in-ground-effect estimate (image vortex system in the VLM/horseshoe model, Wieselsberger-type induced-drag reduction, lift change, wing–tail downwash interference) | computational experiment |
| 21 | xfoil-avl-openvsp | Using XFOIL, AVL & OpenVSP | Validate your codes against the tools (optional practice with free tools; nothing in the tests imports or runs them) | computational experiment: your codes against published NACA data and, where the learner has them, tool output stored as labelled CSVs |
| 22 | aerodynamic-validation-project | Aerodynamic Validation Project | Compare your VLM to wind-tunnel data (a validation report: code verification, solution verification, validation against authentic finite-wing and section data, with uncertainty) | data-based investigation |
| 23 | home-lab-airfoil-bench | Home Lab: Airfoil Lift & Drag on a Bench Rig | Guarded fan with a straw-honeycomb flow straightener, foam or card airfoil on a kitchen-scale lift balance and a pivot-arm drag balance; CL and CD vs α; theory vs computation vs measurement 🔬 | physical experiment, `projects/p10-airfoil-bench` |

Lesson ids are `P10.NN`. Prerequisites: every lesson needs the previous Phase 10 lessons it uses and the earlier lessons it builds on (name them in the metadata row): `P00.13` safety, `P00.14` measurement and uncertainty, `P00.15` first experiment, `P01.05` Taylor series, `P01.09` gradients and optimization, `P01.11` vector fields, line integrals and circulation, `P01.12` Stokes and Green's theorems, `P01.13` linear systems, `P02.08` Laplace transform, `P02.09` Fourier series, `P02.12` complex variables and conformal maps (the Joukowski airfoil), `P02.18` regression and curve fitting, `P03.05` ODE integrators, `P03.08` direct linear solvers, `P03.11` finite differences and observed order, `P03.14` optimization, `P03.16` multidisciplinary design optimization, `P03.18` verification and validation, `P04.04` rotation, `P05.12` beam deflection where a wing bends. The ROADMAP's Phase 10 table also names Phase 08 lessons (`P08.05` dimensional analysis, `P08.08` potential flow, `P08.11` integral boundary layers, `P08.12` transition, `P08.14` turbulent boundary layers, `P08.19` the flow and pressure home lab). **Phase 08 is not written yet.** Refer to those lessons by id only (for example "`P08.11` (Thwaites' method)"), never link to a file in `curriculum/phase-08/`, and give each Phase 10 lesson a short self-contained recap (a few lines, with the formula used) of any Phase 08 result it needs, so the lesson works before Phase 08 exists. The learner has finished Phases 00–05 and, by the roadmap order, 06–09: Python, the scientific stack, calculus, linear algebra, ODEs, numerical methods, physics, statics, and fluid mechanics. What is new is aerodynamics: the forces on lifting surfaces and bodies, and the models that predict them.

## What Phase 10 must achieve

Incompressible aerodynamics at the level of an undergraduate aerodynamics course and the start of a graduate one (MIT 16.100 and 16.110, Unified Engineering fluids; Anderson's *Fundamentals of Aerodynamics* or Drela's *Flight Vehicle Aerodynamics* scope), taught through the design of real aircraft, drones and rotors. The organising idea of the phase is *a hierarchy of models, each checked against the one above and below it*: thin-airfoil theory, panel methods and boundary layers in 2-D; horseshoe vortices, lifting line and vortex lattices in 3-D; momentum, blade element and BEM for rotors; quasi-steady, Theodorsen and Wagner in time. Every lesson says where its model sits in that hierarchy, what it neglects, and checks its main result by an independent route:

- thin-airfoil theory against the panel method as thickness → 0;
- the panel method against the exact Joukowski or Karman–Trefftz solution, with the observed order of convergence;
- lifting line against the elliptic-wing closed form (CL = 2π α AR/(AR + 2) for a flat elliptic wing with a 2π section, e = 1);
- the VLM against lifting line for unswept, high-aspect-ratio wings, and the VLM's lattice convergence;
- BEM against actuator-disk momentum theory;
- Theodorsen against its quasi-steady limit (k → 0) and Wagner against Theodorsen through the Fourier transform;
- every grid, panel or collocation study reports the observed order.

A learner who finishes can compute and transform aerodynamic coefficients; generate airfoils; explain circulation and the Kutta condition; solve thin-airfoil and panel problems; add a boundary layer and predict drag; design and select airfoils; quantify high-lift devices; compute induced drag, span loads and stability derivatives with lifting line and VLM; optimize planforms; build a drag polar; reason about low Reynolds number and stall; estimate body forces; analyse propellers and rotors; compute unsteady lift; estimate ground effect; validate their codes against free professional tools and authentic wind-tunnel data; and measure lift and drag at home on a bench rig and explain why the numbers differ from the published ones. Each comes with a tested Python module.

## Shared conventions (use these so the lessons agree)

- **Air:** ISA sea level unless the lesson says otherwise: ρ = 1.225 kg/m³, μ = 1.789 × 10⁻⁵ Pa·s, ν = 1.461 × 10⁻⁵ m²/s, a = 340.3 m/s, g = 9.81 m/s². Room air for the home lab: 20 °C, 101.3 kPa, ρ = 1.204 kg/m³, μ = 1.81 × 10⁻⁵ Pa·s.
- **Axes and signs:** α positive nose-up; lift perpendicular and drag parallel to the free stream; pitching moment positive nose-up, about the quarter chord unless stated; x along the chord from the leading edge, z up; 2-D section coefficients lower-case (cl, cd, cm), 3-D wing and aircraft coefficients upper-case (CL, CD, Cm). Circulation Γ positive clockwise for a lifting airfoil in a left-to-right stream, so that L′ = ρ V∞ Γ is positive (say this explicitly in lessons 03–05).
- **The reference aircraft** (use it wherever a lesson needs an aeroplane; every number derived from it must be reproduced by the code): a two-seat light trainer, mass 750 kg, wing area S = 12.0 m², span b = 10.0 m (AR = 8.33), rectangular unswept wing with chord 1.2 m, NACA 2412 section, cruise 55 m/s at sea level, approach 30 m/s, fuselage length 7.0 m and maximum diameter 1.1 m, horizontal tail area 2.0 m² at 4.5 m behind the wing quarter chord. Lessons may vary one parameter for a study (taper, sweep) and must say so.
- **The reference small UAV** (lessons 14, 17, 23): span 1.5 m, chord 0.20 m, mass 1.2 kg, cruise 12 m/s, 10 × 4.5 in (0.254 m × 0.114 m) two-blade propeller.
- **The reference rotor** (lesson 18): a light helicopter, gross mass 1100 kg, main rotor radius 5.0 m, 2 blades, chord 0.25 m, tip speed 210 m/s, NACA 0012 blades.
- Lessons that need code from an earlier lesson (a NACA generator, a panel solver, a horseshoe-vortex induced-velocity function) **copy a compact version into their own starter as a "given" helper** with a one-line docstring saying which lesson it comes from. No lesson imports another lesson's file.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-10/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-10/lesson-NN-<slug>` runs against `solutions/phase-10/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish. Also run plain `python3 -m pytest <your lesson folders>` against the starters: **every failure must be `NotImplementedError`** (a test that fails for any other reason against the starter, e.g. a broken given helper or an import error, is a bug). Then run `python3 scripts/check_lessons.py`; it must report no ERROR for your lessons.

The starter has the same function signatures, type hints and docstrings as the solution with `raise NotImplementedError` bodies for the functions the learner writes; helpers that are "given" are complete in both. Keep a working `main()` in the starter (it may raise `NotImplementedError` when it reaches unfinished code, but must import and parse cleanly). Tests assert properties: residuals of the boundary condition, the Kutta condition, conservation (momentum, energy in the actuator disk), symmetry (a symmetric airfoil at α = 0 gives cl = 0; a symmetric wing gives a symmetric span load), linearity in α, limiting cases (thickness → 0, AR → ∞, k → 0, h/b → ∞), agreement between two methods, a textbook closed form, the observed order of a convergence study, dimensional sanity; not one magic number. **12–30 tests per lesson; each lesson's tests run in under 10 s** (time them: `python3 -m pytest --durations=5`). No network access in tests. No files written outside `tmp_path`. Anything random is seeded (`numpy.random.default_rng(seed)`) with statistical tolerances justified in a comment.

Every number in `lesson.md` and `worked-solutions.md` must be reproduced by `solution.py` (its `main()` or a function it exposes) or by a probe script you run while writing; write the code and tests first, then the prose, and copy numbers from the code's output, never from memory.

Code: Python 3.11+, NumPy, SciPy, Matplotlib and SymPy only (no XFOIL, AVL, OpenVSP, aerosandbox or any aerodynamics library in tested code); type hints; docstrings; `if __name__ == "__main__": main()`. Use `matplotlib.use("Agg")` before importing pyplot, close every figure (`plt.close(fig)`), and have `main()` save figures to a path argument (default the current directory) rather than calling `plt.show()`. Tests that exercise plotting write only to `tmp_path`.

## Data rules

- **Authentic data** lives in `curriculum/phase-10/data/` (shared) or a lesson's own `data/` folder, as small CSVs whose first lines are `#` comments giving the source (full citation, NTRS id and URL), the license ("US Government work, public domain"), what was measured and how, the facility and Reynolds and Mach numbers, and how it was transcribed. A `README.md` beside the files repeats this. Lessons 13 and 22 (and 06, 15, 21 where they compare with data) use these files. Never alter a value; corrections (for example a wind-tunnel correction) are made in code and stated.
- **Synthetic data** is labelled synthetic everywhere it appears: file names start with `SYNTHETIC_`, the first comment line says "SYNTHETIC: generated by <function>(seed=…) in solution.py; not a measurement", and every table caption, plot title and sentence that uses it says synthetic. The generator is in the solution and the starter (given).
- **Tool output** (lesson 21): XFOIL, AVL and OpenVSP were not available in the authoring environment, so no tool output ships. Lesson 21 gives the exact commands and settings to produce CSVs, a reader that requires a header recording tool, version and settings, and validates the learner's own codes against the published NACA data instead. Nothing in the tests imports or runs the tools.
- Never present a simulation as a measurement. Never call a computed number "measured".

## Physical experiment (lesson 23)

- Lesson 23 owns `projects/p10-airfoil-bench/` with `bom.yaml` (three tiers: `minimum_cost`, `recommended`, `enhanced`, as in the earlier kits; prices as ranges with the date checked; `reuse_from` naming `p00-calibration`, `p00-first-experiment`, `p04-drag-drop`, `p04-electronics-bench` and `p05-cantilever` where items are reused), `safety.md` (the full safety box in the format of `projects/p05-cantilever/safety.md`), `procedure.md`, `report-template.md` and `data/README.md`, plus `data/synthetic/SYNTHETIC_seed<NNN>_*.csv` generated by the lesson's solution.
- Apparatus: a guarded desk or box fan with a straw-honeycomb flow straightener (drinking straws packed in a frame) and a short settling length; a foam (XPS/EPS insulation board) or folded-card airfoil section of known chord and span, with end plates; a simple balance: a kitchen scale under a sting or stand for lift (reading change with the fan on), and a pivot arm (a bell-crank) bearing on a second scale, or a calibrated spring, for drag. Angle of attack set with a protractor template. Airspeed estimated with the `P08.19` manometer method (a pitot tube and an inclined water manometer) or a cheap vane anemometer; the lesson explains how small the dynamic pressure is (q ≈ 10 Pa at 4 m/s) and what that does to the uncertainty. The learner sweeps α, computes CL and CD with uncertainties, and compares with thin-airfoil theory, their own panel code (with and without a boundary-layer estimate), and the NACA data at its much higher published Reynolds number, and explains the Reynolds-number gap (Re ≈ 3 × 10⁴ on the bench vs 2–9 × 10⁶ in the NACA tunnels: laminar separation, lower lift-curve slope, much higher cd, earlier and softer stall). Blockage and the fan's swirl and non-uniformity are discussed and, where possible, measured.
- **The prediction is written down and dated before measuring** (the predicted CL and CD at each planned α with their uncertainty, and the predicted airspeed).
- Hard limits, stated in `safety.md` and repeated in the lesson's safety box at the top of part 4: **the fan is guarded on both sides (front and back grilles intact, finger-proof), rated at most 50 W, and stands on a stable surface where it cannot tip or walk; no propellers, drones or unguarded blades, ever; no cutting tools near the running fan; hot-wire foam cutting only with adult supervision and ventilation, otherwise a craft knife with the fan unplugged; unplug the fan before adjusting the model** (and before reaching anywhere near the fan). Also: cable routed so it cannot be tripped over or pulled into the fan; loose paper, hair, scarves and string kept away from the intake; model and balance secured so nothing can be blown into the fan; stop if the fan rattles, the guard is damaged, the motor smells hot, or anything touches a grille.
- The learner without a kit completes the lesson with a supplied seeded synthetic dataset, labelled synthetic everywhere it appears.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `| **Time** | unmeasured |`.
- Audience: a motivated adult who has finished Phases 00–05. Plain language, short sentences, define each term on first use. SI units (knots and feet only where industry practice demands, and say so).
- Every lesson states a concrete *engineering question* in the first paragraph with the numerical answer the lesson reaches, derives in part 2, builds in part 3, runs an **experiment** in part 4 (physical in lesson 23; data-based in 13 and 22; otherwise a computational experiment: a parameter sweep, a sensitivity study, a method comparison, a convergence study, or a comparison with an authentic public dataset whose source and license are named), validates in part 5, and ends part 6 with `### Adversarial self-checks` and then `## Check yourself` with 5–8 questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: the reference trainer, gliders, small UAVs, airliner high-lift systems, propellers, helicopter rotors, gust and flutter, a wing-in-ground-effect craft. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-10.md` (written before the lessons; re-read it before you finalise citations). Where the map says GAP or Partial GAP, write original instruction and say nothing about missing references; cite only ids that exist in `references/sources.yaml` at the time you write. Never cite or require a paid resource: Anderson, Katz & Plotkin, Leishman's *Principles of Helicopter Aerodynamics*, Abbott & von Doenhoff, Hoerner, Raymer and Drela's *Flight Vehicle Aerodynamics* appear only under an "Optional textbooks for deeper study" line in the References row, with what they add. Never paste text from a source; write your own explanation and link. Leishman's *Introduction to Aerospace Flight Vehicles* is CC BY-NC-ND: link to it by chapter title and URL, never by chapter number, and never adapt its figures.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation, one unfamiliar scenario, one that uses the lesson's code, one judgment question about model limits. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Honesty rules: never claim the course grants a degree, credential or equivalence; never present simulation as physical measurement; say what a model assumes and where it fails; never invent time estimates.
- Length: `lesson.md` roughly 250–450 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders (and `projects/p10-airfoil-bench/` for the lesson 23 author). `curriculum/phase-10/data/` is read-only for lesson authors.
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Name any scratch or probe scripts uniquely (prefix them with your lesson numbers, e.g. `p10_l05_probe.py`) and keep them in the scratchpad, not in the repository.
- Do not add time estimates, grades, or badges.
