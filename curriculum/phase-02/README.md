# Phase 02: Mathematics II: Differential Equations, Transforms & Probability

Differential equations, transforms and probability at the level of a first-year-plus engineering sequence, taught through aerospace problems and always paired with computation: every analytic solution is checked against a numerical integration or a Monte Carlo run, and every numerical result is checked against a limiting case or a closed form. By the end you can solve first- and second-order ODEs and linear systems by hand and with `scipy.integrate.solve_ivp`, use Laplace and Fourier methods to reason about response and spectra, set up separation-of-variables PDE solutions, map a circle to an airfoil, derive Euler–Lagrange equations, manipulate tensors in index notation, and model uncertainty with distributions, hypothesis tests, regression with confidence bands, random-process spectra and Bayesian updates, with a tested Python module for each. There are no physical experiments in this phase; every lesson's part 4 is a computational experiment or a data-based investigation, and any synthetic data is labelled as synthetic wherever it appears.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [First-Order ODEs](lesson-01-first-order-odes/lesson.md) | Rocket mass depletion model | computational |
| 02 | [Second-Order Linear ODEs](lesson-02-second-order-linear-odes/lesson.md) | Damped spring-mass solver | computational |
| 03 | [Forced Response & Resonance](lesson-03-forced-response-resonance/lesson.md) | Resonance explorer | computational |
| 04 | [Systems of ODEs & Phase Portraits](lesson-04-systems-and-phase-portraits/lesson.md) | Phase portrait plotter | computational |
| 05 | [Matrix Exponential & Linear Systems](lesson-05-matrix-exponential/lesson.md) | State transition matrix | computational |
| 06 | [Nonlinear ODEs, Equilibria & Linearization](lesson-06-nonlinear-equilibria-linearization/lesson.md) | Pendulum linearization | computational |
| 07 | [Complex Numbers & Phasors](lesson-07-complex-numbers-phasors/lesson.md) | Phasor calculator | computational |
| 08 | [Laplace Transform](lesson-08-laplace-transform/lesson.md) | Transfer-function solver for ODEs | computational |
| 09 | [Fourier Series](lesson-09-fourier-series/lesson.md) | Square-wave synthesis | computational |
| 10 | [Fourier Transform & the FFT](lesson-10-fourier-transform-fft/lesson.md) | FFT of vibration data | data-based |
| 11 | [PDEs: Heat, Wave & Laplace Equations](lesson-11-pdes-heat-wave-laplace/lesson.md) | Separation-of-variables solutions | computational |
| 12 | [Complex Variables & Conformal Maps](lesson-12-complex-variables-conformal-maps/lesson.md) | Joukowski airfoil from a circle | computational |
| 13 | [Calculus of Variations](lesson-13-calculus-of-variations/lesson.md) | Brachistochrone and minimum-time paths | computational |
| 14 | [Tensor Notation & Index Gymnastics](lesson-14-tensor-notation/lesson.md) | Stress tensor transformations | computational |
| 15 | [Probability Fundamentals](lesson-15-probability-fundamentals/lesson.md) | Monte Carlo dice-to-dispersion | computational |
| 16 | [Random Variables & Distributions](lesson-16-random-variables-distributions/lesson.md) | Landing dispersion ellipse | computational |
| 17 | [Statistics & Hypothesis Testing](lesson-17-statistics-hypothesis-testing/lesson.md) | Is this material batch in spec? | data-based |
| 18 | [Regression & Curve Fitting](lesson-18-regression-curve-fitting/lesson.md) | Thrust curve fit with confidence bands | data-based |
| 19 | [Random Processes & Power Spectra](lesson-19-random-processes-power-spectra/lesson.md) | Gust spectrum (Dryden) generator | computational |
| 20 | [Bayesian Inference](lesson-20-bayesian-inference/lesson.md) | Update a failure-rate estimate | computational |

**Phase exam and mini-project:** [assessments/phase-02/](../../assessments/phase-02/). The thrust-curve fit with confidence bands from lesson 18 is the mini-project.

**References:** [references/phase-02.md](../../references/phase-02.md). Every required reading is free.

**Running the code:** from the repo root, `pytest curriculum/phase-02/lesson-NN-<slug>` tests your `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-02/lesson-NN-<slug>` runs the same tests against the reference solution.
