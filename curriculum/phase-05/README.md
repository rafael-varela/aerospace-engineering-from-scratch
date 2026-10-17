# Phase 05: Statics & Mechanics of Materials

First-year statics and second-year strength of materials, at the level of MIT 2.001 and 16.001 or a Hibbeler or Gere course, taught through aerospace structures and organised around one idea: *equilibrium, compatibility and constitutive law*. Every lesson draws the free-body diagram first and says which of the three it is using at each step. Equilibrium means the forces and moments on any part you cut free sum to zero. Compatibility means the deformed pieces still fit together. The constitutive law is the material's stress–strain relation. Every main result is then checked by an independent route: method of joints against method of sections, composite section properties against the polygon (Green's theorem) formulas, the flexure formula against a finite-difference or finite-element beam from Phase 03, Mohr's circle against `numpy.linalg.eigh`, Castigliano against direct integration, Euler buckling against a finite-difference eigenvalue problem, and in lesson 18 beam theory against your own measurement. By the end you can draw correct free-body diagrams and solve 2-D and 3-D equilibrium. You can analyse trusses, frames and control linkages with friction, and compute section properties. You can relate stress and strain and size members in tension, torsion and bending. You can compute shear flow, deflections and buckling loads, transform stress, and apply failure criteria with aerospace margins of safety. You can use energy methods for deflections and redundant forces, model and dimension a part in free CAD software, and measure a beam's deflection at home and compare it with theory. You get a tested Python module for each of these. The mathematics is not new: Phases 01–03 supplied it. What is new is engineering mechanics: internal forces, stress, strain, deflection and failure of structural members.

One lesson has a physical experiment (marked 🔬): the cantilever home lab. A learner without the kit completes it with the synthetic dataset supplied, which is labelled synthetic wherever it appears; a calculation or a simulation is never presented as a measurement. Lesson 17 is a practical activity in free CAD software with a computational check. Every other lesson's part 4 is a computational experiment, such as a parameter sweep, a sensitivity study, a comparison of two methods, or an analysis of a public or synthetic dataset.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Forces, Moments & Free-Body Diagrams](lesson-01-forces-moments-free-body-diagrams/lesson.md) | FBD solver for a landing gear | computational |
| 02 | [Equilibrium of Rigid Bodies](lesson-02-equilibrium-rigid-bodies/lesson.md) | Reaction solver (2-D and 3-D, statically determinate, with a check for indeterminacy) | computational |
| 03 | [Trusses: Joints & Sections](lesson-03-trusses-joints-sections/lesson.md) | Truss solver from scratch (method of joints as a linear system, method of sections as a check) | computational |
| 04 | [Frames, Machines & Friction](lesson-04-frames-machines-friction/lesson.md) | Control linkage analysis (pushrod–bellcrank–horn, mechanical advantage, friction in pins and cables) | computational |
| 05 | [Centroids & Moments of Inertia](lesson-05-centroids-moments-of-inertia/lesson.md) | Section property calculator (composite and polygon sections, parallel axis, principal axes) | computational |
| 06 | [Stress & Strain](lesson-06-stress-strain/lesson.md) | Stress–strain curve analyzer (modulus, 0.2% offset yield, UTS, elongation, toughness) | computational, on a public or synthetic tensile dataset |
| 07 | [Axial Loading & Thermal Stress](lesson-07-axial-loading-thermal-stress/lesson.md) | Bolted joint preload (joint stiffness diagram, external load sharing, thermal mismatch) | computational |
| 08 | [Torsion of Shafts](lesson-08-torsion-of-shafts/lesson.md) | Turbine shaft sizing (solid and hollow circular, thin-walled closed sections by Bredt–Batho) | computational |
| 09 | [Shear & Bending Moment Diagrams](lesson-09-shear-bending-moment-diagrams/lesson.md) | V–M diagram generator (point, distributed and moment loads; singularity functions) | computational |
| 10 | [Bending Stress in Beams](lesson-10-bending-stress-beams/lesson.md) | Wing spar stress (flexure formula, unsymmetric bending preview, spar cap sizing) | computational |
| 11 | [Shear Stress in Beams](lesson-11-shear-stress-beams/lesson.md) | Shear flow in an I-beam (VQ/It, open thin-walled sections, shear centre of a channel) | computational |
| 12 | [Beam Deflection](lesson-12-beam-deflection/lesson.md) | Deflection by integration and superposition, with a finite-difference check | computational |
| 13 | [Stress Transformation & Mohr's Circle](lesson-13-stress-transformation-mohrs-circle/lesson.md) | Mohr's circle tool (2-D and 3-D principal stresses via eigenvalues) | computational |
| 14 | [Failure Theories](lesson-14-failure-theories/lesson.md) | von Mises and Tresca margin calculator (aerospace factors 1.0 yield and 1.5 ultimate; brittle Mohr–Coulomb) | computational |
| 15 | [Column Buckling](lesson-15-column-buckling/lesson.md) | Euler and Johnson buckling (end conditions, effective length, tangent modulus) | computational |
| 16 | [Energy Methods](lesson-16-energy-methods/lesson.md) | Castigliano deflection solver (axial, bending and torsion energy; unit-load method; one redundant) | computational |
| 17 | [CAD & Technical Drawings I](lesson-17-cad-technical-drawings/lesson.md) | Model and dimension a bracket; produce a drawing | practical activity in free CAD software, plus a computational check |
| 18 | [Home Lab: Cantilever Beam Deflection](lesson-18-home-lab-cantilever-deflection/lesson.md) 🔬 | Load a clamped beam, measure deflection, compare with theory | physical experiment ([project](../../projects/p05-cantilever/)) |

Lesson ids are `P05.01` to `P05.18`. Each lesson's metadata table gives its exact prerequisites, and every lesson's metadata table reads:

| **Time** | unmeasured |
|---|---|

The Time row stays `unmeasured` until pilots measure it. No time estimate is given here for the same reason.

## Prerequisites

Phases 00–04, in particular: `P00.13` (safety for home experiments) and `P00.14` (measurement, calibration and uncertainty) for lesson 18, and `P00.15` (the first predict-measure-compare experiment), whose method lesson 18 repeats at a higher level; `P01.07` (vectors) for forces and moments as r × F; `P01.10` (multiple integrals) for centroids and moments of inertia; `P01.13` (linear systems) for the reaction and truss solvers; `P01.16` (eigenvalues) for principal axes, principal stresses and buckling; `P02.02` (second-order ODEs) for beam deflection and the column equation; `P02.14` (tensors) for stress transformation; `P03.08` (direct linear solvers) and `P03.10` (nonlinear systems) for the solvers and the friction iterations; `P03.13` (finite elements) for the independent beam check; `P04.01` (Newton's laws) for free-body diagrams; and `P04.04` (rotation) for moments and torque. Skill checks before starting: you can draw a free-body diagram that names the body exerting each force, solve a 3 × 3 linear system by hand, find the eigenvalues of a symmetric 2 × 2 matrix, integrate a polynomial twice with two boundary conditions, and propagate an uncertainty through a product of powers (the thickness of a beam enters its stiffness cubed).

## The project kit

The kit lives in `projects/p05-cantilever/` with a `bom.yaml` (three cost tiers, `minimum_cost`, `recommended` and `enhanced`, with dated price ranges, reusing the kitchen scale from `p00-calibration` and the phone-video habits of the Phase 04 kits), `safety.md` (the full safety box), `procedure.md`, `report-template.md` and `data/README.md` (the raw-data format). Read `safety.md` before buying anything and again before every session. The kit's `safety.md` is authoritative. The limits below summarise the hard rules that every version of this kit must keep, and nothing in the kit may relax them. The prediction, with its uncertainty, is written down and dated before the first load, raw data are never edited, and the report states which numbers were measured, which were computed, and which are synthetic.

- **`p05-cantilever`** ([folder](../../projects/p05-cantilever/), lesson 18). A ruler, a wooden or aluminium strip, or a plastic ruler clamped to a table edge; known masses (coins, water bottles, weighed on the kitchen scale); the deflection read with a second ruler, a phone camera against a grid, or a dial gauge at the enhanced tier. You predict the deflection from E and I, both with uncertainties, before loading; then you measure load against deflection, fit E, and compare the deflection along the beam with the cubic of beam theory. Hard limits:
  - total hanging mass ≤ 2 kg;
  - no loading of glass, carbon-fibre rods or anything that can splinter;
  - safety glasses when loading strips that can snap;
  - the clamp checked before each load;
  - nothing hung over feet;
  - stop at visible permanent set or cracking.

Nothing in this phase involves power tools, pressure vessels, heating, combustion, projectiles, flight or loading anything to failure. Lesson 06 analyses a tensile test supplied as a public or synthetic dataset; there is no destructive testing at home. **A report without a completed safety box is not graded.**

## Lesson 17 and FreeCAD

Lesson 17 is done in FreeCAD (free, open source, available for Windows, macOS and Linux) or any free CAD tool you prefer. The lesson's instructions are tool-neutral, with FreeCAD menu names given where they help. The cited manual, `@freecad-manual-van-havre` ("Modeling for product design" and "Generating 2D drawings"), was written for FreeCAD 0.16, so some menus and the drawing workflow look different in current releases; the lesson's own steps are written for the current release. Drawing conventions (projection, dimensioning, tolerancing) follow `@ford-engineering-graphics-and-design` ("Orthographic Projection", "Dimensioning Standards", "Conventional Tolerancing"), which uses third-angle projection. The tested code checks the *drawing*, not the CAD file: the bracket described by parameters, its mass properties from polygon and extrusion formulas, a tolerance stack-up, and a check that the dimension scheme fully constrains the part without redundancy. Nothing in the tests imports a CAD package, so the tests run without FreeCAD installed.

## Running the code

From the repo root, `pytest curriculum/phase-05/lesson-NN-<slug>` tests your `starter/main.py`. `AEFS_IMPL=solution pytest curriculum/phase-05/lesson-NN-<slug>` runs the same tests against the reference solution in `solutions/phase-05/lesson-NN-<slug>/solution.py`. Tests check properties rather than single numbers:
- equilibrium residuals ΣF = 0 and ΣM = 0 to round-off;
- agreement between two methods (joints and sections, composite and polygon, closed form and finite differences, Castigliano and integration, Mohr's circle and `eigh`);
- limiting cases (a vanishing load, an infinitely stiff member, a symmetric section, a pinned end with zero moment);
- symmetry, superposition and dimensional sanity;
- convergence with grid size where a numerical method is used.

The code uses NumPy, SciPy and Matplotlib only; no structural or CAD library is used in tested code. Anything random is seeded. Tests need no network and write files only to a temporary directory. `python scripts/check_lessons.py` checks every lesson's structure.

## The exam and the mini-project

**Phase exam and mini-project:** [assessments/phase-05/](../../assessments/phase-05/). The exam has four parts. Part A covers concepts: free-body diagrams and determinacy, stress and strain, failure and margins, beam theory and energy methods. Part B covers analysis by hand over every lesson from 01 to 17. Part C is computation: a tested finite-difference buckling solver checked against the Euler load. Part D is engineering interpretation. Its first question asks you to interpret your own `p05-cantilever` measurement; a seeded synthetic dataset is supplied for learners without a kit. Its second question reviews a synthetic stress check sheet for wrong factors and missing failure modes, and tests the cantilever kit's safety limits. The exam has not been piloted, so no time limit is suggested. Choose one in advance and record it. Results are reported per competency: (1) reasoning, (2) analysis, (3) computation, (4) interpretation. The worked solutions and pass thresholds are in [solutions/phase-05/exam/](../../solutions/phase-05/exam/worked-solutions.md). The mini-project is the cantilever report from lesson 18, graded with [the rubric](../../assessments/phase-05/grading-rubric.md). Seeded unseen cases for the computational lessons are planned in [parameterized-cases/](../../assessments/phase-05/parameterized-cases/README.md). Outcomes are recorded as *practice completed*, *self-assessed mastery* or *externally reviewed mastery*. None of these is a credential, and nothing in this phase qualifies anyone to analyse or approve real aircraft structure.

**References:** [references/phase-05.md](../../references/phase-05.md). Every required reading is free. Baker and Haynes's *Engineering Statics: Open and Interactive* (`@libretexts-baker-haynes-engineering-statics`) carries the statics lessons 01–05 and 09. Roylance's MIT 3.11 modules (`@mit-ocw-3-11-mechanics-of-materials`) and the MIT 16.001 lecture notes and labs (`@mit-ocw-16-001-unified-materials-structures`) carry most of the mechanics of materials, with `@nasa-rp-1228-fastener-design-manual` for bolted joints. The CAD lesson uses `@freecad-manual-van-havre` and `@ford-engineering-graphics-and-design`. The lessons add the further sources named in the reference map. The curriculum links to every source and writes its own explanations.
