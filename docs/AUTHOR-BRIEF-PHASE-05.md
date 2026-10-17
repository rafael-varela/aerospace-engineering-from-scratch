# Author brief: Phase 05 lessons

This brief is what every lesson author in Phase 05 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-05.md`, `references/sources.yaml`, the Phase 05 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/`, `solutions/_template/` and `projects/_template/` folders. **Phases 00–04 are the quality bar:** before writing, read `curriculum/phase-04/lesson-04-rotation-angular-momentum/lesson.md` and `curriculum/phase-04/lesson-09-magnetism-induction/` (lesson, starter, tests, solution) end to end and match their depth, tone and structure. For a lesson with a physical experiment, also read `curriculum/phase-04/lesson-06-oscillations-waves/lesson.md` and every file in `projects/p04-pendulum/`; that is the bar for the experiment, the bill of materials and the safety box.

## Where files go

```
curriculum/phase-05/lesson-NN-<slug>/lesson.md
curriculum/phase-05/lesson-NN-<slug>/problems.md
curriculum/phase-05/lesson-NN-<slug>/starter/main.py
curriculum/phase-05/lesson-NN-<slug>/tests/test_main.py
solutions/phase-05/lesson-NN-<slug>/hints.md
solutions/phase-05/lesson-NN-<slug>/method.md
solutions/phase-05/lesson-NN-<slug>/worked-solutions.md
solutions/phase-05/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
projects/p05-<id>/bom.yaml, safety.md, procedure.md, report-template.md, data/README.md   <- physical-experiment lessons only
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build | mode of part 4 |
|---|---|---|---|---|
| 01 | forces-moments-free-body-diagrams | Forces, Moments & Free-Body Diagrams | FBD solver for a landing gear | computational experiment |
| 02 | equilibrium-rigid-bodies | Equilibrium of Rigid Bodies | Reaction solver (2-D and 3-D, statically determinate and a check for indeterminacy) | computational experiment |
| 03 | trusses-joints-sections | Trusses: Joints & Sections | Truss solver from scratch (method of joints as a linear system, method of sections as a check) | computational experiment |
| 04 | frames-machines-friction | Frames, Machines & Friction | Control linkage analysis (pushrod–bellcrank–horn, mechanical advantage, friction in pins and cables) | computational experiment |
| 05 | centroids-moments-of-inertia | Centroids & Moments of Inertia | Section property calculator (composite and polygon sections, parallel-axis, principal axes) | computational experiment |
| 06 | stress-strain | Stress & Strain | Stress-strain curve analyzer (modulus, 0.2 % offset yield, UTS, elongation, toughness) | computational experiment on a public or synthetic tensile dataset |
| 07 | axial-loading-thermal-stress | Axial Loading & Thermal Stress | Bolted joint preload (joint stiffness diagram, external load sharing, thermal mismatch) | computational experiment |
| 08 | torsion-of-shafts | Torsion of Shafts | Turbine shaft sizing (solid and hollow circular, thin-walled closed sections with Bredt–Batho) | computational experiment |
| 09 | shear-bending-moment-diagrams | Shear & Bending Moment Diagrams | V–M diagram generator (point, distributed, moment loads; singularity functions) | computational experiment |
| 10 | bending-stress-beams | Bending Stress in Beams | Wing spar stress (flexure formula, unsymmetric bending preview, spar cap sizing) | computational experiment |
| 11 | shear-stress-beams | Shear Stress in Beams | Shear flow in an I-beam (VQ/It, open thin-walled sections, shear centre of a channel) | computational experiment |
| 12 | beam-deflection | Beam Deflection | Deflection by integration and superposition (and a finite-difference check) | computational experiment |
| 13 | stress-transformation-mohrs-circle | Stress Transformation & Mohr's Circle | Mohr's circle tool (2-D and 3-D principal stresses via eigenvalues) | computational experiment |
| 14 | failure-theories | Failure Theories | von Mises, Tresca margin calculator (margins of safety with aerospace factors 1.0 yield / 1.5 ultimate, brittle Mohr–Coulomb) | computational experiment |
| 15 | column-buckling | Column Buckling | Euler and Johnson buckling (end conditions, effective length, tangent-modulus discussion) | computational experiment |
| 16 | energy-methods | Energy Methods | Castigliano deflection solver (strain energy in axial, bending and torsion; unit-load method; statically indeterminate redundant) | computational experiment |
| 17 | cad-technical-drawings | CAD & Technical Drawings I | Model and dimension a bracket; produce a drawing | practical activity in free CAD software (FreeCAD, open source) plus a computational check |
| 18 | home-lab-cantilever-deflection | Home Lab: Cantilever Beam Deflection | Load a clamped beam, measure deflection, compare with theory | physical experiment, `projects/p05-cantilever` |

Lesson ids are `P05.NN`. Prerequisites: every lesson needs the previous Phase 05 lessons it uses and the earlier lessons it builds on (name them: `P00.13` safety, `P00.14` measurement and uncertainty, `P00.15` first experiment, `P01.07` vectors, `P01.10` multiple integrals, `P01.13` linear systems, `P01.16` eigenvalues, `P02.02` second-order ODEs, `P02.14` tensors, `P03.08` direct linear solvers, `P03.10` nonlinear systems, `P03.13` finite elements, `P04.01` Newton's laws, `P04.04` rotation); name them in the metadata row. The learner has finished Phases 00–04: Python, the scientific stack, calculus, linear algebra, ODEs, numerical methods and introductory physics are known. What is new is engineering mechanics: equilibrium, internal forces, stress, strain, deflection and failure of structural members.

## What Phase 05 must achieve

Statics and mechanics of materials at the level of a first-year statics course and a second-year strength-of-materials course (MIT 2.001, 16.001; Hibbeler or Gere scope), taught through aerospace structures. The organising idea of the phase is *equilibrium, compatibility and constitutive law*: every lesson says which of the three it uses, draws the free-body diagram first, and checks its answer by an independent route (a second method, a limiting case, a finite-element or finite-difference solution from Phase 03, or a measurement). A learner who finishes can draw correct free-body diagrams, solve 2-D and 3-D equilibrium, analyse trusses, frames and linkages with friction, compute section properties, relate stress and strain, size members in tension, torsion and bending, compute shear flow, deflections and buckling loads, transform stress and apply failure criteria with aerospace margins of safety, use energy methods for deflection and redundants, model and dimension a part in free CAD software, and measure a beam deflection at home and compare it with theory, with a tested Python module for each.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-05/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-05/lesson-NN-<slug>` runs against `solutions/phase-05/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run `python scripts/check_lessons.py` (structure linter; it must report no ERROR for your lessons). The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (equilibrium residuals ΣF = 0 and ΣM = 0 to round-off, symmetry, superposition, limiting cases such as a vanishing load or an infinitely stiff member, agreement between two methods, a textbook closed form, dimensional sanity), not one magic number. 12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. Tests run in under 10 s per lesson.

Every lesson has an independent check that its main result must reproduce: method of joints against method of sections, a parallel-axis composite against polygon (Green's theorem) section properties, the flexure formula against a 1-D FEM or finite-difference beam from Phase 03, Mohr's circle against `numpy.linalg.eigh`, Castigliano against direct integration, Euler buckling against a finite-difference eigenvalue problem. Use NumPy/SciPy only; no structural or CAD libraries in tested code. Anything random is seeded (`numpy.random.default_rng(seed)`) with statistical tolerances justified in a comment.

Lesson 17 (CAD) is done in FreeCAD (free, open source, any platform) or any free CAD tool the learner prefers; the lesson gives tool-neutral instructions with FreeCAD menu names cited from its documentation. Its tested code checks the drawing, not the CAD file: a bracket described by parameters, its mass properties by polygon/extrusion formulas, a tolerance stack-up, and a check that a dimension scheme fully constrains the part without redundancy. Nothing in the tests imports a CAD package.

## Physical experiment (lesson 18)

- Lesson 18 owns `projects/p05-cantilever/` with `bom.yaml` (three tiers: `minimum_cost`, `recommended`, `enhanced`, as in the earlier kits; prices as ranges with the date checked; `reuse_from` naming the p00 and p04 kits), `safety.md` (the full safety box in the format of `projects/p04-pendulum/safety.md`), `procedure.md`, `report-template.md` and `data/README.md`.
- Apparatus: a ruler, wooden or aluminium strip, or a plastic ruler clamped to a table edge; known masses (coins, water bottles, a kitchen scale from `p00-calibration` to weigh them); deflection read with a second ruler, a phone camera against a grid, or a dial gauge at the enhanced tier. The prediction (from E and I, both with uncertainties) is written down and dated before loading. The learner measures load–deflection, fits E, and compares tip deflection along the beam with the cubic.
- Hard limits stated in `safety.md` and repeated in the lesson: total hanging mass ≤ 2 kg; no loading of glass, carbon-fibre rods or anything that can splinter; safety glasses when loading strips that can snap; clamp checked before each load; nothing hung over feet; stop at visible permanent set or cracking.
- The learner without a kit completes the lesson with a supplied synthetic dataset, labelled synthetic everywhere it appears.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phase 00. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs an **experiment** in part 4 (physical in lesson 18, CAD practice in 17; otherwise a computational experiment such as a parameter sweep, a sensitivity study, a method comparison or a comparison with an authentic public dataset whose source and license are named), validates in part 5, and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: a landing-gear strut and drag brace, an engine mount truss, a control-surface pushrod and bellcrank, wing-spar and stringer sections, bolted and riveted joints, a turbine or propeller shaft, a wing spar under a spanwise lift distribution, a fuselage frame, thin-walled shear flow, margins of safety with aerospace factors of safety (1.0 limit-to-yield, 1.5 ultimate), struts and longerons in compression. Materials: 2024-T3 and 7075-T6 aluminium, Ti-6Al-4V, 4130 steel, with properties stated with their source and labelled typical. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-05.md` (being written now; read it when you start, and re-read it before you finalise citations). Where the map says GAP or partial GAP, write original instruction and say nothing about missing references; cite only ids that exist in `references/sources.yaml` at the time you write. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link. Leishman's book is CC BY-NC-ND: link to it by chapter title and URL, never by chapter number, and never adapt its figures.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only; type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders (and `projects/p05-cantilever/` for the lesson 18 author).
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not call any `mcp__hearthbot__` tools, do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
