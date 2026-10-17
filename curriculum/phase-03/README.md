# Phase 03: Numerical Methods & Scientific Computing

Numerical methods and scientific computing at the level of a first course in numerical analysis plus an introduction to computational engineering, taught through aerospace problems and organised around one idea: *error*. Every method is introduced with its error model (truncation, rounding, stability, convergence order), every implementation is written from scratch in NumPy and checked against that model with a convergence study, every result carries an error estimate, and SciPy is called afterwards as the oracle, with a reason given for why the library is still what you call in practice. By the end you can implement and analyse root finders, interpolants and splines, quadrature and finite differences, ODE integrators (explicit, adaptive, implicit and symplectic), direct and iterative linear solvers, Newton for nonlinear systems, finite-difference, finite-volume and finite-element discretisations, unconstrained and constrained optimisers, a small multidisciplinary optimisation with coupled derivatives, Monte Carlo and polynomial-chaos uncertainty propagation, and a verification study by the method of manufactured solutions, with a tested Python module for each. There are no physical experiments in this phase; every lesson's part 4 is a computational experiment (a convergence, stability or cost study), and any synthetic data is labelled as synthetic wherever it appears.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Floating Point & Numerical Error](lesson-01-floating-point-numerical-error/lesson.md) | Catastrophic cancellation demo | computational |
| 02 | [Root Finding](lesson-02-root-finding/lesson.md) | Newton, secant, bisection; solve Kepler's equation | computational |
| 03 | [Interpolation & Splines](lesson-03-interpolation-splines/lesson.md) | Atmospheric table interpolator | computational |
| 04 | [Numerical Differentiation & Integration](lesson-04-numerical-differentiation-integration/lesson.md) | Simpson, Gauss quadrature | computational |
| 05 | [ODE Integrators I: Euler to RK4](lesson-05-ode-integrators-euler-rk4/lesson.md) | RK4 from scratch | computational |
| 06 | [ODE Integrators II: Adaptive & Stiff](lesson-06-ode-integrators-adaptive-stiff/lesson.md) | RK45 with step control, implicit Euler | computational |
| 07 | [Symplectic Integrators](lesson-07-symplectic-integrators/lesson.md) | Energy-conserving orbit propagation | computational |
| 08 | [Direct Linear Solvers](lesson-08-direct-linear-solvers/lesson.md) | LU, Cholesky from scratch | computational |
| 09 | [Iterative & Sparse Solvers](lesson-09-iterative-sparse-solvers/lesson.md) | Conjugate gradient on a 2D Poisson problem | computational |
| 10 | [Nonlinear Systems](lesson-10-nonlinear-systems/lesson.md) | Newton–Raphson for trim | computational |
| 11 | [Finite Differences for PDEs](lesson-11-finite-differences-pdes/lesson.md) | 1D heat equation, stability analysis | computational |
| 12 | [Finite Volume Methods](lesson-12-finite-volume-methods/lesson.md) | 1D advection with upwinding | computational |
| 13 | [Finite Element Basics](lesson-13-finite-element-basics/lesson.md) | 1D bar FEM solver | computational |
| 14 | [Unconstrained Optimization](lesson-14-unconstrained-optimization/lesson.md) | BFGS from scratch | computational |
| 15 | [Constrained Optimization](lesson-15-constrained-optimization/lesson.md) | SQP for a sizing problem | computational |
| 16 | [Multidisciplinary Design Optimization](lesson-16-multidisciplinary-design-optimization/lesson.md) | Wing aero-structural MDO (NumPy/SciPy; OpenMDAO as an optional extension) | computational |
| 17 | [Uncertainty Quantification](lesson-17-uncertainty-quantification/lesson.md) | Monte Carlo and polynomial chaos | computational |
| 18 | [Verification & Validation](lesson-18-verification-and-validation/lesson.md) | Method of manufactured solutions | computational |

**Phase exam and mini-project:** [assessments/phase-03/](../../assessments/phase-03/). The verification study by the method of manufactured solutions from lesson 18, applied to your own lesson 11 or lesson 12 solver, is the mini-project.

**References:** [references/phase-03.md](../../references/phase-03.md). Every required reading is free.

**Running the code:** from the repo root, `pytest curriculum/phase-03/lesson-NN-<slug>` tests your `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-03/lesson-NN-<slug>` runs the same tests against the reference solution.
