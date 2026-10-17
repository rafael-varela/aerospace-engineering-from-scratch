# Author brief: Phase 02 lessons

This brief is what every lesson author in Phase 02 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-02.md`, `references/sources.yaml`, the Phase 02 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/` and `solutions/_template/` folders. **Phases 00 and 01 are the quality bar:** before writing, read `curriculum/phase-01/lesson-02-applications-of-derivatives/lesson.md` and `curriculum/phase-01/lesson-16-eigenvalues-eigenvectors/` (lesson, starter, tests, solution) end to end and match their depth, tone and structure.

## Where files go

```
curriculum/phase-02/lesson-NN-<slug>/lesson.md
curriculum/phase-02/lesson-NN-<slug>/problems.md
curriculum/phase-02/lesson-NN-<slug>/starter/main.py
curriculum/phase-02/lesson-NN-<slug>/tests/test_main.py
solutions/phase-02/lesson-NN-<slug>/hints.md
solutions/phase-02/lesson-NN-<slug>/method.md
solutions/phase-02/lesson-NN-<slug>/worked-solutions.md
solutions/phase-02/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build |
|---|---|---|---|
| 01 | first-order-odes | First-Order ODEs | Rocket mass depletion model |
| 02 | second-order-linear-odes | Second-Order Linear ODEs | Damped spring-mass solver |
| 03 | forced-response-resonance | Forced Response & Resonance | Resonance explorer |
| 04 | systems-and-phase-portraits | Systems of ODEs & Phase Portraits | Phase portrait plotter |
| 05 | matrix-exponential | Matrix Exponential & Linear Systems | State transition matrix |
| 06 | nonlinear-equilibria-linearization | Nonlinear ODEs, Equilibria & Linearization | Pendulum linearization |
| 07 | complex-numbers-phasors | Complex Numbers & Phasors | Phasor calculator |
| 08 | laplace-transform | Laplace Transform | Transfer-function solver for ODEs |
| 09 | fourier-series | Fourier Series | Square-wave synthesis |
| 10 | fourier-transform-fft | Fourier Transform & the FFT | FFT of vibration data |
| 11 | pdes-heat-wave-laplace | PDEs: Heat, Wave & Laplace Equations | Separation-of-variables solutions |
| 12 | complex-variables-conformal-maps | Complex Variables & Conformal Maps | Joukowski airfoil from a circle |
| 13 | calculus-of-variations | Calculus of Variations | Brachistochrone and minimum-time paths |
| 14 | tensor-notation | Tensor Notation & Index Gymnastics | Stress tensor transformations |
| 15 | probability-fundamentals | Probability Fundamentals | Monte Carlo dice-to-dispersion |
| 16 | random-variables-distributions | Random Variables & Distributions | Landing dispersion ellipse |
| 17 | statistics-hypothesis-testing | Statistics & Hypothesis Testing | Is this material batch in spec? |
| 18 | regression-curve-fitting | Regression & Curve Fitting | Thrust curve fit with confidence bands |
| 19 | random-processes-power-spectra | Random Processes & Power Spectra | Gust spectrum (Dryden) generator |
| 20 | bayesian-inference | Bayesian Inference | Update a failure-rate estimate |

Lesson ids are `P02.NN`. Prerequisites: every lesson needs the previous Phase 02 lessons it uses, the Phase 01 lessons it builds on (name them: `P01.04` integration techniques, `P01.05` series, `P01.12` vector theorems, `P01.15` rotations, `P01.16` eigenvalues, `P01.17` least squares), and `P00.09` (SymPy) where symbolic work appears; name them in the metadata row. The learner has finished Phases 00 and 01: Python, NumPy, SciPy, Matplotlib, SymPy and pytest are routine, and single- and multivariable calculus and linear algebra are known.

## What Phase 02 must achieve

Differential equations, transforms and probability at the level of a first-year-plus engineering sequence (MIT 18.03, OpenStax Calculus 2–3 and Introductory Statistics), taught through aerospace problems and *always* paired with computation: every analytic solution is checked against a numerical integration or a Monte Carlo run, and every numerical result is checked against a limiting case or a closed form. Derivations are shown in full in part 2 of each lesson. A learner who finishes can solve first- and second-order ODEs and linear systems by hand and by `scipy.integrate.solve_ivp`, use Laplace and Fourier methods to reason about response and spectra, set up separation-of-variables PDE solutions, map a circle to an airfoil, derive Euler–Lagrange equations, manipulate tensors in index notation, and model uncertainty with distributions, hypothesis tests, regression with confidence bands, random-process spectra and Bayesian updates, with a tested Python module for each.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-02/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-02/lesson-NN-<slug>` runs against `solutions/phase-02/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run `python scripts/check_lessons.py` (structure linter; it must report no ERROR for your lessons). The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (residuals, limiting cases, conservation, symmetry, convergence order, agreement between symbolic and numeric), not one magic number. 12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. Tests run in under 10 s per lesson.

Numerical lessons implement the core algorithm *from scratch* (Euler and RK4 steppers, a DFT before calling the FFT, a finite-difference PDE marcher, a Metropolis or grid posterior) and the tests compare against SciPy (`solve_ivp`, `scipy.fft`, `scipy.signal`, `scipy.stats`, `scipy.linalg.expm`) as the oracle; the lesson explains why we still call the library in practice. Anything random is seeded (`numpy.random.default_rng(seed)`) and tests use fixed seeds with statistical tolerances justified in a comment.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phase 00. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs a **computational experiment** or **data-based investigation** in part 4 (no physical experiments in this phase; synthetic data must be labelled as synthetic everywhere it appears), validates in part 5, and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: rocket mass depletion and burnout, landing-gear and panel vibration, flutter-style resonance, phugoid and short-period modes, attitude dynamics, the pendulum and the inertia swing test, transfer functions of an actuator, vibration spectra from an accelerometer, heat conduction in a re-entry tile and the vibrating string or membrane, the Joukowski airfoil, minimum-time climb paths, stress tensor rotation, landing dispersion, material acceptance testing, thrust-curve fits, Dryden gust spectra, failure-rate updates. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-02.md`. Where the map says GAP or lists a lesson under "Lessons without a mapped free source", write original instruction and say nothing about missing references; cite only ids that exist in `references/sources.yaml` at the time you write. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only (no statsmodels, no control library); type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders.
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not call any `mcp__hearthbot__` tools, do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
