# Author brief: Phase 03 lessons

This brief is what every lesson author in Phase 03 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-03.md`, `references/sources.yaml`, the Phase 03 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/` and `solutions/_template/` folders. **Phases 00–02 are the quality bar:** before writing, read `curriculum/phase-02/lesson-05-matrix-exponential/lesson.md` and `curriculum/phase-02/lesson-10-fourier-transform-fft/` (lesson, starter, tests, solution) end to end and match their depth, tone and structure.

## Where files go

```
curriculum/phase-03/lesson-NN-<slug>/lesson.md
curriculum/phase-03/lesson-NN-<slug>/problems.md
curriculum/phase-03/lesson-NN-<slug>/starter/main.py
curriculum/phase-03/lesson-NN-<slug>/tests/test_main.py
solutions/phase-03/lesson-NN-<slug>/hints.md
solutions/phase-03/lesson-NN-<slug>/method.md
solutions/phase-03/lesson-NN-<slug>/worked-solutions.md
solutions/phase-03/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build |
|---|---|---|---|
| 01 | floating-point-numerical-error | Floating Point & Numerical Error | Catastrophic cancellation demo |
| 02 | root-finding | Root Finding | Newton, secant, bisection; solve Kepler's equation |
| 03 | interpolation-splines | Interpolation & Splines | Atmospheric table interpolator |
| 04 | numerical-differentiation-integration | Numerical Differentiation & Integration | Simpson, Gauss quadrature |
| 05 | ode-integrators-euler-rk4 | ODE Integrators I: Euler to RK4 | RK4 from scratch |
| 06 | ode-integrators-adaptive-stiff | ODE Integrators II: Adaptive & Stiff | RK45 with step control, implicit Euler |
| 07 | symplectic-integrators | Symplectic Integrators | Energy-conserving orbit propagation |
| 08 | direct-linear-solvers | Direct Linear Solvers | LU, Cholesky from scratch |
| 09 | iterative-sparse-solvers | Iterative & Sparse Solvers | Conjugate gradient on a 2D Poisson problem |
| 10 | nonlinear-systems | Nonlinear Systems | Newton–Raphson for trim |
| 11 | finite-differences-pdes | Finite Differences for PDEs | 1D heat equation, stability analysis |
| 12 | finite-volume-methods | Finite Volume Methods | 1D advection with upwinding |
| 13 | finite-element-basics | Finite Element Basics | 1D bar FEM solver |
| 14 | unconstrained-optimization | Unconstrained Optimization | BFGS from scratch |
| 15 | constrained-optimization | Constrained Optimization | SQP for a sizing problem |
| 16 | multidisciplinary-design-optimization | Multidisciplinary Design Optimization | Wing aero-structural MDO (NumPy/SciPy; OpenMDAO as an optional extension) |
| 17 | uncertainty-quantification | Uncertainty Quantification | Monte Carlo and polynomial chaos |
| 18 | verification-and-validation | Verification & Validation | Method of manufactured solutions |

Lesson ids are `P03.NN`. Prerequisites: every lesson needs the previous Phase 03 lessons it uses and the Phase 01/02 lessons it builds on (name them: `P01.05` series, `P01.09` optimization, `P01.13` linear systems, `P01.16` eigenvalues, `P01.17` least squares, `P02.04` systems of ODEs, `P02.11` PDEs, `P02.16` distributions); name them in the metadata row. The learner has finished Phases 00–02: Python and the scientific stack are routine, and calculus, linear algebra, ODEs, transforms and probability are known.

## What Phase 03 must achieve

Numerical methods and scientific computing at the level of a first course in numerical analysis plus an introduction to computational engineering (MIT 18.330, 18.335J, 16.90, 2.086), taught through aerospace problems. The organising idea of the phase is *error*: every method is introduced with its error model (truncation, rounding, stability, convergence order), every implementation is checked against that model with a convergence study, and every result carries an error estimate. A learner who finishes can implement and analyse root finders, interpolants, quadrature, ODE integrators (explicit, adaptive, implicit, symplectic), direct and iterative linear solvers, Newton for nonlinear systems, finite-difference, finite-volume and finite-element discretisations, unconstrained and constrained optimizers, a small multidisciplinary optimisation, Monte Carlo and polynomial-chaos uncertainty propagation, and a verification study by the method of manufactured solutions, with a tested Python module for each and SciPy as the oracle.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-03/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-03/lesson-NN-<slug>` runs against `solutions/phase-03/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run `python scripts/check_lessons.py` (structure linter; it must report no ERROR for your lessons). The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (residuals, limiting cases, conservation, symmetry, convergence order, agreement between symbolic and numeric), not one magic number. 12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. Tests run in under 10 s per lesson.

Every lesson implements its core algorithm *from scratch* in NumPy (bisection/Newton/secant, Lagrange and cubic splines, Gauss–Legendre nodes, RK4 and an embedded RK45 with step control, implicit Euler with Newton, leapfrog/Störmer–Verlet, LU and Cholesky, conjugate gradient, Newton–Raphson, FTCS/BTCS/Crank–Nicolson, upwind/Lax–Wendroff finite volume, linear-element FEM assembly, BFGS with line search, an SQP or augmented-Lagrangian step, a Gauss–Seidel multidisciplinary analysis with a coupled derivative, Monte Carlo and a Hermite polynomial-chaos expansion, Richardson extrapolation) and the tests compare against SciPy (`scipy.optimize`, `scipy.interpolate`, `scipy.integrate`, `scipy.linalg`, `scipy.sparse.linalg`) as the oracle; the lesson explains why we still call the library in practice. Every lesson includes at least one *observed convergence order* test: run the method at two or three resolutions and assert the measured order against the theoretical one within a stated tolerance. Anything random is seeded (`numpy.random.default_rng(seed)`) with statistical tolerances justified in a comment.

Lesson 16 implements the MDO workflow in NumPy/SciPy only (a two-discipline wing aero-structural model with a fixed-point multidisciplinary analysis, finite-difference and coupled-adjoint derivatives, and `scipy.optimize.minimize` as the driver); OpenMDAO appears only in part 6 and the Extensions row as an optional, separately installable tool, and nothing in the tests imports it. Anything random is seeded (`numpy.random.default_rng(seed)`) and tests use fixed seeds with statistical tolerances justified in a comment.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phase 00. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs a **computational experiment** in part 4 (a convergence, stability or cost study; no physical experiments in this phase; synthetic data must be labelled as synthetic everywhere it appears), validates in part 5, and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: the ISA atmosphere table, Kepler's equation, orbit propagation and energy drift, re-entry heating, trim of an aircraft, a 2-D Poisson problem for potential flow or heat, a bar or truss FEM, a wing-box or fuselage sizing problem, a wing aero-structural trade (span vs weight), uncertainty in a range or Δv estimate, verification of a diffusion or advection solver. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-03.md`. Where the map says GAP (adaptive step-size theory in 03-06, symplectic integration in 03-07, polynomial chaos in 03-17, grid-convergence studies in 03-18), write original instruction and say nothing about missing references; cite only ids that exist in `references/sources.yaml` at the time you write. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only (no statsmodels, no control library, no OpenMDAO in tested code); type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders.
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not call any `mcp__hearthbot__` tools, do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
