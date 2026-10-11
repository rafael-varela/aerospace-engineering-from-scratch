# Phase 00: Setup, Engineering Computing & First Measurements

The phase every path starts with. By the end you have a working Python engineering environment, you can compute, plot, test, and document like an engineer, you have read a classic NACA report, and you have done your first real experiment: a prediction written down first, a safe measurement with calibrated instruments, and an honest comparison.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Dev Environment for Engineers](lesson-01-dev-environment/lesson.md) | Python, uv, VS Code, Jupyter set up and verified | computational |
| 02 | [Git & Collaboration](lesson-02-git-and-collaboration/lesson.md) | Fork, branch, PR on this repo | design-build |
| 03 | [NumPy as an Engineering Calculator](lesson-03-numpy-engineering-calculator/lesson.md) | Vectorized ISA atmosphere table | computational |
| 04 | [Plotting That Engineers Trust](lesson-04-plotting-engineers-trust/lesson.md) | Publication-quality plot style module | computational |
| 05 | [Units & Dimensional Consistency](lesson-05-units-and-dimensions/lesson.md) | Unit-checked calculations with pint; Mars Climate Orbiter post-mortem | data-based |
| 06 | [SciPy Tour](lesson-06-scipy-tour/lesson.md) | Integrate, optimize, interpolate on a real dataset | computational |
| 07 | [Testing Engineering Code](lesson-07-testing-engineering-code/lesson.md) | pytest suite with tolerance-based assertions | computational |
| 08 | [Reading Data: CSV, HDF5, Telemetry Logs](lesson-08-reading-data/lesson.md) | Parser for a flight log | data-based |
| 09 | [Symbolic Math with SymPy](lesson-09-symbolic-math-sympy/lesson.md) | Derive and lambdify an equation of motion | computational |
| 10 | [Writing Technical Reports](lesson-10-writing-technical-reports/lesson.md) | Report template and a report linter | design-build |
| 11 | [Reading a Technical Paper](lesson-11-reading-a-technical-paper/lesson.md) | Annotated summary of NACA Report 824 | data-based |
| 12 | [Engineering Judgment & Sanity Checks](lesson-12-engineering-judgment/lesson.md) | Order-of-magnitude checker for an airliner | computational |
| 13 | [Safety & Standards for Home Experiments](lesson-13-safety-for-home-experiments/lesson.md) | Hazard checklist and stop-criteria template | design-build |
| 14 | [Measurement, Calibration & Uncertainty](lesson-14-measurement-calibration-uncertainty/lesson.md) 🔬 | Calibrate a scale and a ruler; report uncertainty | physical experiment ([project](../../projects/p00-calibration/)) |
| 15 | [First Experiment: Predict, Measure, Compare](lesson-15-first-experiment/lesson.md) 🔬 | Time a falling object with phone video; compare with prediction | physical experiment ([project](../../projects/p00-first-experiment/)) |

**Phase exam and mini-project:** [assessments/phase-00/](../../assessments/phase-00/). The free-fall report from lesson 15 is the mini-project.

**References:** [references/phase-00.md](../../references/phase-00.md). Every required reading is free.

**Running the code:** from the repo root, `pytest curriculum/phase-00/lesson-NN-<slug>` tests your `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-00/lesson-NN-<slug>` runs the same tests against the reference solution.
