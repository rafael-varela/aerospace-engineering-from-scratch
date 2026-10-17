# Author brief: Phase 04 lessons

This brief is what every lesson author in Phase 04 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-04.md`, `references/sources.yaml`, the Phase 04 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/`, `solutions/_template/` and `projects/_template/` folders. **Phases 00–03 are the quality bar:** before writing, read `curriculum/phase-03/lesson-05-ode-integrators-euler-rk4/lesson.md` and `curriculum/phase-03/lesson-10-nonlinear-systems/` (lesson, starter, tests, solution) end to end and match their depth, tone and structure. For a lesson with a physical experiment, also read `curriculum/phase-00/lesson-15-first-experiment/lesson.md` and every file in `projects/p00-first-experiment/`; that is the bar for the experiment, the bill of materials and the safety box.

## Where files go

```
curriculum/phase-04/lesson-NN-<slug>/lesson.md
curriculum/phase-04/lesson-NN-<slug>/problems.md
curriculum/phase-04/lesson-NN-<slug>/starter/main.py
curriculum/phase-04/lesson-NN-<slug>/tests/test_main.py
solutions/phase-04/lesson-NN-<slug>/hints.md
solutions/phase-04/lesson-NN-<slug>/method.md
solutions/phase-04/lesson-NN-<slug>/worked-solutions.md
solutions/phase-04/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
projects/p04-<id>/bom.yaml, safety.md, procedure.md, report-template.md, data/README.md   <- physical-experiment lessons only
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build | mode of part 4 |
|---|---|---|---|---|
| 01 | kinematics-newtons-laws | Kinematics & Newton's Laws | Projectile with drag | physical experiment, `projects/p04-drag-drop` |
| 02 | work-energy-power | Work, Energy & Power | Energy-height of an aircraft | computational experiment |
| 03 | momentum-collisions | Momentum & Collisions | Variable-mass system intuition | computational experiment |
| 04 | rotation-angular-momentum | Rotation & Angular Momentum | Spinning top simulation | computational experiment |
| 05 | gravitation-keplers-laws | Gravitation & Kepler's Laws | Two-body simulation | computational experiment |
| 06 | oscillations-waves | Oscillations & Waves | Wave propagation visualizer | physical experiment, `projects/p04-pendulum` |
| 07 | kinetic-theory-of-gases | Kinetic Theory of Gases | Molecular speed distributions | computational experiment |
| 08 | electrostatics-circuits | Electrostatics & Circuits | DC circuit solver | computational experiment (optional measurement with the `p04-electronics-bench` kit, clearly marked optional) |
| 09 | magnetism-induction | Magnetism & Induction | Magnetorquer field model | computational experiment |
| 10 | ac-circuits-electronics-basics | AC Circuits & Electronics Basics | RC filter design | computational experiment (optional measurement with the `p04-electronics-bench` kit, clearly marked optional) |
| 11 | electromagnetic-waves-antennas | Electromagnetic Waves & Antennas | Dipole radiation pattern | computational experiment |
| 12 | optics-sensors | Optics & Sensors | Star tracker optics budget | computational experiment |
| 13 | special-relativity-gps | Special Relativity for GPS | Clock corrections in GNSS | computational experiment |
| 14 | atomic-radiation-physics | Atomic & Radiation Physics | Radiation dose basics | computational experiment (no sources, no detectors at home) |
| 15 | electronics-bench-lab | Electronics Bench Lab | Breadboard, multimeter, microcontroller: measure voltage, current and a sensor | physical experiment, `projects/p04-electronics-bench` |

Lesson ids are `P04.NN`. Prerequisites: every lesson needs the previous Phase 04 lessons it uses and the earlier lessons it builds on (name them: `P00.13` safety, `P00.14` measurement and uncertainty, `P01.07` vectors, `P01.13` linear systems, `P02.02` second-order ODEs, `P02.04` systems of ODEs, `P02.07` phasors, `P02.16` distributions, `P03.05` RK4, `P03.07` symplectic integrators); name them in the metadata row. The learner has finished Phases 00–03: Python, the scientific stack, calculus, linear algebra, ODEs, transforms, probability and numerical methods are known. Phase 04 is the first *physics* phase, so what is new is the physics, not the mathematics: derive from first principles, then compute.

## What Phase 04 must achieve

Calculus-based introductory physics at the level of a two-semester university sequence (OpenStax University Physics volumes 1–3; MIT 8.01, 8.02 and parts of 8.03), taught through aerospace problems and with every result checked against measurement or an independent computation. The organising idea of the phase is *conservation and modelling*: every lesson states the conservation law or field equation that governs its system, builds the simplest model that captures the engineering question, says what the model assumes and where it fails, and compares it against data (a home experiment, a public dataset, or an independent numerical solution). A learner who finishes can set up and solve Newtonian mechanics problems (particles, rigid bodies, orbits, oscillators, waves), explain pressure and temperature from kinetic theory, analyse DC and AC circuits and simple electronics, model magnetic fields and induction, compute electromagnetic wave power, antenna patterns and link losses, size a camera or star-tracker optic, compute the relativistic clock corrections GNSS needs, estimate radiation dose and shielding, and build, measure and debug a low-voltage circuit on a breadboard with a multimeter and a microcontroller, with a tested Python module for each.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-04/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-04/lesson-NN-<slug>` runs against `solutions/phase-04/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run `python scripts/check_lessons.py` (structure linter; it must report no ERROR for your lessons). The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (conservation of energy, momentum and angular momentum to a stated tolerance; limiting cases such as zero drag, zero damping, small angles, non-relativistic speeds; symmetry; dimensional sanity; convergence of a simulation with step size; agreement with a closed form or with `scipy.integrate.solve_ivp` as the oracle), not one magic number. 12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. Tests run in under 10 s per lesson.

Every lesson has a closed-form or textbook result that its simulation must reproduce in a limit (free fall without drag, the vis-viva relation, the small-angle pendulum period, the RC time constant, the Biot–Savart field on the axis of a loop, the inverse-square intensity, the GPS 38 μs/day net offset, the half-life law), and at least one test checks that limit. Simulations use the integrators the learner built in Phase 03 only by calling `scipy.integrate.solve_ivp` (the lesson may remind them which method is which); a lesson may implement one short integrator by hand when the point is conservation (leapfrog in 05). Anything random is seeded (`numpy.random.default_rng(seed)`) with statistical tolerances justified in a comment.

## Physical experiments (lessons 01, 06, 15; optional in 08 and 10)

- Each physical-experiment lesson owns a folder `projects/p04-<id>/` with `bom.yaml` (three tiers: minimum cost, standard, extended; prices as ranges with the date checked; `reuse_from` naming earlier kits), `safety.md` (the full safety box in the format of `projects/p00-first-experiment/safety.md`, with hazards, energy limits, controls, supervision, location, regulations, stop criteria, disposal, not permitted), `procedure.md`, `report-template.md` and `data/README.md` describing the raw-data format. The lesson's metadata rows **Practical activity**, **Equipment** and **Safety** link to that folder, and part 4 puts a safety summary box first.
- `p04-drag-drop` (01): drop or toss light objects (ping-pong ball, paper cone, bean bag) and time with a phone camera; reuses the `p00-first-experiment` kit; the prediction (with and without quadratic drag) is written down before measuring.
- `p04-pendulum` (06): a string pendulum; period against length and against amplitude; phone video or phone accelerometer; string ≤ 2 m, bob ≤ 0.2 kg, no swinging near people or screens.
- `p04-electronics-bench` (15): breadboard, resistors, LEDs, a potentiometer, a thermistor or photoresistor, a multimeter and a 3.3 V or 5 V microcontroller board programmed over USB (any board with an ADC; the procedure must not depend on one vendor). Hard limits, stated in `safety.md` and repeated in the lesson: battery or USB power only, never mains; supply ≤ 9 V DC; no lithium-cell charging circuits; every LED through a resistor; current limited to what the supply can give (≤ 500 mA); no soldering required (if the learner solders, the standard safety for a soldering iron applies and is written down); no capacitors above 25 V or 1000 μF; discharge capacitors before touching.
- Every physical activity is optional in the sense that a learner without the kit completes the lesson with the supplied synthetic or public data; the lesson says so and labels every synthetic dataset as synthetic everywhere it appears. Nothing in this phase involves combustion, pressure vessels, projectiles launched by stored energy, flight, lasers above class 1, radioactive sources, or mains voltage.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phase 00. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs an **experiment** in part 4 (physical where the table says so; otherwise a computational experiment such as a parameter sweep, a sensitivity study or a comparison with an authentic public dataset whose source and license are named), validates in part 5, and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: the wind triangle and terminal speed, the energy height and specific excess power of an aircraft, the rocket equation from momentum conservation, reaction wheels and gyroscopic precession, orbital energy and Kepler's laws, the ISA atmosphere from kinetic theory, a satellite power bus and its wiring losses, a magnetorquer and the Earth's dipole field, an RC anti-aliasing filter for a sensor, a link budget and radiation pressure on a solar sail, the f-number and pixel scale of a star tracker, the GPS clock offset, dose behind aluminium shielding. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-04.md`. Where the map says GAP or partial GAP (antennas in 04-11; op-amps, logic and microcontrollers in 04-10; generic sensors in 04-12; the GPS worked example in 04-13; RTGs and single-event upsets in 04-14; breadboarding and multimeter practice in 04-15), write original instruction and say nothing about missing references; cite only ids that exist in `references/sources.yaml` at the time you write. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link. Leishman's book is CC BY-NC-ND: link to it by chapter title and URL, never by chapter number, and never adapt its figures.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only; type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders (and your assigned `projects/p04-*` folder).
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not call any `mcp__hearthbot__` tools, do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
