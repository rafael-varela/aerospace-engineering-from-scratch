# Phase 01: Mathematics I: Calculus & Linear Algebra

Calculus and linear algebra at the level of a first-year engineering sequence, taught through aerospace problems and always paired with computation: every analytic result is checked numerically and every numerical result is checked against a limiting case. By the end you can differentiate and integrate by hand, expand in Taylor series and bound the error, work in vectors and multivariable calculus, and solve, factor, decompose and condition linear systems, with a tested Python module for each. There are no physical experiments in this phase; every lesson's part 4 is a computational experiment or a data-based investigation.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Limits, Continuity & Derivatives](lesson-01-limits-and-derivatives/lesson.md) | Numerical vs symbolic derivative comparison | computational |
| 02 | [Applications of Derivatives](lesson-02-applications-of-derivatives/lesson.md) | Max range glide from a drag polar | computational |
| 03 | [Integrals & the Fundamental Theorem](lesson-03-integrals-fundamental-theorem/lesson.md) | Area and centroid of an airfoil section | computational |
| 04 | [Techniques of Integration](lesson-04-techniques-of-integration/lesson.md) | Moment of inertia of a tapered beam | computational |
| 05 | [Sequences, Series & Taylor Expansions](lesson-05-series-and-taylor/lesson.md) | Small-angle approximation error explorer | computational |
| 06 | [Parametric Curves & Polar Coordinates](lesson-06-parametric-and-polar/lesson.md) | Plot a conic orbit from r(θ) | computational |
| 07 | [Vectors, Dot & Cross Products](lesson-07-vectors-dot-cross/lesson.md) | Vector library for 3D mechanics | computational |
| 08 | [Multivariable Functions & Partial Derivatives](lesson-08-partial-derivatives/lesson.md) | Sensitivity of lift to α and V | computational |
| 09 | [Gradients, Directional Derivatives & Optimization](lesson-09-gradients-and-optimization/lesson.md) | Minimum-drag speed by gradient search | computational |
| 10 | [Multiple Integrals](lesson-10-multiple-integrals/lesson.md) | Mass properties of a wing box | computational |
| 11 | [Vector Fields, Line & Surface Integrals](lesson-11-line-and-surface-integrals/lesson.md) | Circulation around a vortex | computational |
| 12 | [Divergence, Stokes & Green's Theorems](lesson-12-divergence-stokes-green/lesson.md) | Verify the divergence theorem numerically | computational |
| 13 | [Systems of Linear Equations](lesson-13-linear-systems/lesson.md) | Gaussian elimination from scratch | computational |
| 14 | [Vector Spaces, Basis & Rank](lesson-14-vector-spaces-basis-rank/lesson.md) | Detect a mechanism in a truss | computational |
| 15 | [Matrix Transformations & Rotation Matrices](lesson-15-rotation-matrices/lesson.md) | Direction cosine matrix tools | computational |
| 16 | [Eigenvalues & Eigenvectors](lesson-16-eigenvalues-eigenvectors/lesson.md) | Principal axes of an inertia tensor | computational |
| 17 | [Orthogonality, Least Squares & QR](lesson-17-least-squares-qr/lesson.md) | Fit a drag polar to wind-tunnel data | data-based |
| 18 | [SVD & Matrix Norms](lesson-18-svd-and-norms/lesson.md) | Conditioning of a structural stiffness matrix | computational |

**Phase exam and mini-project:** [assessments/phase-01/](../../assessments/phase-01/). The wind-tunnel drag-polar fit report from lesson 17 is the mini-project.

**References:** [references/phase-01.md](../../references/phase-01.md). Every required reading is free.

**Running the code:** from the repo root, `pytest curriculum/phase-01/lesson-NN-<slug>` tests your `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-01/lesson-NN-<slug>` runs the same tests against the reference solution.
