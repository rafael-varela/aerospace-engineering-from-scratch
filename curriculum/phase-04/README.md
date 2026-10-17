# Phase 04: Physics Foundations

Calculus-based introductory physics at the level of a two-semester university sequence, taught through aerospace problems and organised around one idea: *conservation and modelling*. Every lesson states the conservation law or field equation that governs its system, builds the simplest model that answers an engineering question, says what that model assumes and where it fails, and checks it against data: a home experiment, a public dataset, or an independent numerical solution. Every simulation must also reproduce a closed-form result in a limit: free fall without drag, the vis-viva relation, the small-angle pendulum period, the RC time constant, the Biot–Savart field on the axis of a loop, the inverse-square intensity, the 38 μs/day net GPS clock offset, and the half-life law. By the end you can set up and solve Newtonian mechanics problems (particles, rigid bodies, orbits, oscillators, waves), explain pressure and temperature from kinetic theory, analyse DC and AC circuits and simple electronics, model magnetic fields and induction, compute electromagnetic wave power, antenna patterns and link losses, size a camera or star-tracker optic, compute the relativistic clock corrections that GNSS needs, estimate radiation dose and shielding, and build, measure and debug a low-voltage circuit on a breadboard with a multimeter and a microcontroller. You get a tested Python module for each of these. The mathematics is not new: Phases 01–03 supplied it. The physics is new, so every result is derived from first principles before it is computed.

Three lessons have a physical experiment (marked 🔬), and two more have an optional measurement with the electronics kit. Every physical activity is optional in one sense: a learner without the kit completes the lesson with the synthetic or public data supplied. Every synthetic dataset is labelled synthetic wherever it appears, and a simulation is never presented as a measurement. Every other lesson's part 4 is a computational experiment, such as a parameter sweep, a sensitivity study, or a comparison with an independent solution or a public dataset.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Kinematics & Newton's Laws](lesson-01-kinematics-newtons-laws/lesson.md) 🔬 | Projectile with drag | physical experiment ([project](../../projects/p04-drag-drop/)) |
| 02 | [Work, Energy & Power](lesson-02-work-energy-power/lesson.md) | Energy-height of an aircraft | computational |
| 03 | [Momentum & Collisions](lesson-03-momentum-collisions/lesson.md) | Variable-mass system intuition | computational |
| 04 | [Rotation & Angular Momentum](lesson-04-rotation-angular-momentum/lesson.md) | Spinning top simulation | computational |
| 05 | [Gravitation & Kepler's Laws](lesson-05-gravitation-keplers-laws/lesson.md) | Two-body simulation | computational |
| 06 | [Oscillations & Waves](lesson-06-oscillations-waves/lesson.md) 🔬 | Wave propagation visualizer | physical experiment ([project](../../projects/p04-pendulum/)) |
| 07 | [Kinetic Theory of Gases](lesson-07-kinetic-theory-of-gases/lesson.md) | Molecular speed distributions | computational |
| 08 | [Electrostatics & Circuits](lesson-08-electrostatics-circuits/lesson.md) | DC circuit solver | computational (optional measurement with the [bench kit](../../projects/p04-electronics-bench/)) |
| 09 | [Magnetism & Induction](lesson-09-magnetism-induction/lesson.md) | Magnetorquer field model | computational |
| 10 | [AC Circuits & Electronics Basics](lesson-10-ac-circuits-electronics-basics/lesson.md) | RC filter design | computational (optional measurement with the [bench kit](../../projects/p04-electronics-bench/)) |
| 11 | [Electromagnetic Waves & Antennas](lesson-11-electromagnetic-waves-antennas/lesson.md) | Dipole radiation pattern | computational |
| 12 | [Optics & Sensors](lesson-12-optics-sensors/lesson.md) | Star tracker optics budget | computational |
| 13 | [Special Relativity for GPS](lesson-13-special-relativity-gps/lesson.md) | Clock corrections in GNSS | computational |
| 14 | [Atomic & Radiation Physics](lesson-14-atomic-radiation-physics/lesson.md) | Radiation dose basics | computational (no sources, no detectors at home) |
| 15 | [Electronics Bench Lab](lesson-15-electronics-bench-lab/lesson.md) 🔬 | Breadboard, multimeter, microcontroller: measure voltage, current and a sensor | physical experiment ([project](../../projects/p04-electronics-bench/)) |

Lesson ids are `P04.01` to `P04.15`. Each lesson's metadata table gives its exact prerequisites, and every lesson's **Time** row reads `unmeasured` until pilots measure it. No time estimate is given here for the same reason.

## Prerequisites

Phases 00–03, in particular: `P00.13` (safety for home experiments) and `P00.14` (measurement, calibration and uncertainty) for every experiment; `P00.15` (the first predict-measure-compare experiment), whose method lessons 01, 06 and 15 repeat at a higher level; `P01.07` (vectors) for forces, fields and orbits; `P01.13` (linear systems) for the circuit solver; `P02.02` (second-order ODEs) and `P02.04` (systems of ODEs) for oscillators and orbits; `P02.07` (phasors) for AC circuits; `P02.16` (distributions) for molecular speeds and radioactive decay; `P03.05` (RK4) and `P03.07` (symplectic integrators) for every simulation. Simulations call `scipy.integrate.solve_ivp`, and lesson 05 also writes a short leapfrog integrator by hand because there the point is conservation. Skill checks before starting: you can solve a 3 × 3 linear system by hand, write the general solution of a damped oscillator, convert between a complex phasor and an amplitude and phase, and propagate an uncertainty through a product of powers.

## The three project kits

Each kit lives in `projects/p04-<id>/` with a `bom.yaml` (three cost tiers with dated price ranges and reuse of earlier kits), `safety.md` (the full safety box), `procedure.md`, `report-template.md` and `data/README.md` (the raw-data format). Read `safety.md` before buying anything and again before every session. The kit's `safety.md` is authoritative. The limits below summarise the hard rules that every version of these kits must keep, and nothing in a kit may relax them. In every kit the prediction, with its uncertainty, is written down and dated before the first measurement, raw data are never edited, and the report states which numbers were measured, which were computed, and which are synthetic.

- **`p04-drag-drop`** ([folder](../../projects/p04-drag-drop/), lesson 01). Light objects (a ping-pong ball, a paper cone, a bean bag) are dropped or tossed gently by hand and timed with phone video. The kit reuses `p00-first-experiment`. Before measuring, you write down the prediction both without drag and with quadratic drag. Hard limits: light, soft objects only, at the mass and height limits stated in the kit's `safety.md`; hand release only, with no launcher of any kind (nothing that stores energy to launch a projectile is allowed anywhere in this phase); a clear drop zone; release from the floor or a stable step, never from a chair, ladder, window, balcony or stairwell.
- **`p04-pendulum`** ([folder](../../projects/p04-pendulum/), lesson 06). A string pendulum, used to measure period against length and against amplitude with phone video or the phone's accelerometer. Hard limits: string ≤ 2 m; bob ≤ 0.2 kg; no swinging near people or screens; a clear swing zone with stop criteria.
- **`p04-electronics-bench`** ([folder](../../projects/p04-electronics-bench/), lesson 15; optional in lessons 08 and 10). A breadboard, resistors, LEDs, a potentiometer, a thermistor or photoresistor, a multimeter, and a 3.3 V or 5 V microcontroller board with an ADC, programmed over USB. Any such board works; the procedure does not depend on one vendor. Hard limits:
  - battery or USB power only, never mains;
  - supply ≤ 9 V DC;
  - no lithium-cell charging circuits;
  - every LED through a resistor;
  - current limited to what the supply can give, ≤ 500 mA;
  - no soldering required (if you do solder, the standard soldering-iron safety applies and is written into your safety box);
  - no capacitors above 25 V or 1000 μF;
  - discharge capacitors before touching them.

Nothing in this phase involves combustion, pressure vessels, projectiles launched by stored energy, flight, lasers above class 1, radioactive sources or mains voltage. Lesson 14 computes dose and shielding from tabulated data only; there are no sources and no detectors at home. **A report without a completed safety box is not graded.**

## Running the code

From the repo root, `pytest curriculum/phase-04/lesson-NN-<slug>` tests your `starter/main.py`. `AEFS_IMPL=solution pytest curriculum/phase-04/lesson-NN-<slug>` runs the same tests against the reference solution in `solutions/phase-04/lesson-NN-<slug>/solution.py`. Tests check properties rather than single numbers:
- conservation of energy, momentum and angular momentum to a stated tolerance;
- limiting cases (zero drag, zero damping, small angles, non-relativistic speeds);
- symmetry and dimensional sanity;
- convergence with step size;
- agreement with a closed form or with `solve_ivp` as the oracle.

Anything random is seeded. Tests need no network and write files only to a temporary directory. `python scripts/check_lessons.py` checks every lesson's structure.

## The exam and the mini-project

**Phase exam and mini-project:** [assessments/phase-04/](../../assessments/phase-04/). The exam has four parts. Part A covers concepts. Part B covers analysis by hand over every lesson. Part C is computation: a tested pendulum-period function checked against the elliptic integral. Part D is engineering interpretation. Part D's first question asks you to interpret your own `p04-pendulum` or `p04-drag-drop` measurement; a synthetic dataset is supplied for learners without a kit. Its second question interprets a synthetic bench-lab table and tests the electronics kit's safety limits. The exam has not been piloted, so no time limit is suggested. Choose one in advance and record it. Results are reported per competency: (1) reasoning, (2) analysis, (3) computation, (4) interpretation. The worked solutions and pass thresholds are in [solutions/phase-04/exam/](../../solutions/phase-04/exam/worked-solutions.md). The mini-project is the pendulum report from lesson 06, graded with [the rubric](../../assessments/phase-04/grading-rubric.md). The drag-drop report or the bench-lab report may be substituted on the same terms. Seeded unseen cases for the computational lessons are planned in [parameterized-cases/](../../assessments/phase-04/parameterized-cases/README.md). Outcomes are recorded as *practice completed*, *self-assessed mastery* or *externally reviewed mastery*. None of these is a credential.

**References:** [references/phase-04.md](../../references/phase-04.md). Every required reading is free. The OpenStax *University Physics* volumes (`@openstax-university-physics-vol1`, `@openstax-university-physics-vol2`, `@openstax-university-physics-vol3`) carry most of the phase. The lessons on electronics, antennas, sensors, GPS, radiation effects and the bench lab add the sources named in the reference map. The curriculum links to every source and writes its own explanations.
