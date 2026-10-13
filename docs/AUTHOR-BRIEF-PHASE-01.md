# Author brief: Phase 01 lessons

This brief is what every lesson author in Phase 01 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-01.md`, `references/sources.yaml`, the Phase 01 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/` and `solutions/_template/` folders. **Phase 00 is the quality bar:** before writing, read `curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md` and `curriculum/phase-00/lesson-09-symbolic-math-sympy/` (lesson, starter, tests, solution) end to end and match their depth, tone and structure.

## Where files go

```
curriculum/phase-01/lesson-NN-<slug>/lesson.md
curriculum/phase-01/lesson-NN-<slug>/problems.md
curriculum/phase-01/lesson-NN-<slug>/starter/main.py
curriculum/phase-01/lesson-NN-<slug>/tests/test_main.py
solutions/phase-01/lesson-NN-<slug>/hints.md
solutions/phase-01/lesson-NN-<slug>/method.md
solutions/phase-01/lesson-NN-<slug>/worked-solutions.md
solutions/phase-01/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build |
|---|---|---|---|
| 01 | limits-and-derivatives | Limits, Continuity & Derivatives | Numerical vs symbolic derivative comparison |
| 02 | applications-of-derivatives | Applications of Derivatives | Max range glide from a drag polar |
| 03 | integrals-fundamental-theorem | Integrals & the Fundamental Theorem | Area and centroid of an airfoil section |
| 04 | techniques-of-integration | Techniques of Integration | Moment of inertia of a tapered beam |
| 05 | series-and-taylor | Sequences, Series & Taylor Expansions | Small-angle approximation error explorer |
| 06 | parametric-and-polar | Parametric Curves & Polar Coordinates | Plot a conic orbit from r(θ) |
| 07 | vectors-dot-cross | Vectors, Dot & Cross Products | Vector library for 3D mechanics |
| 08 | partial-derivatives | Multivariable Functions & Partial Derivatives | Sensitivity of lift to α and V |
| 09 | gradients-and-optimization | Gradients, Directional Derivatives & Optimization | Minimum-drag speed by gradient search |
| 10 | multiple-integrals | Multiple Integrals | Mass properties of a wing box |
| 11 | line-and-surface-integrals | Vector Fields, Line & Surface Integrals | Circulation around a vortex |
| 12 | divergence-stokes-green | Divergence, Stokes & Green's Theorems | Verify the divergence theorem numerically |
| 13 | linear-systems | Systems of Linear Equations | Gaussian elimination from scratch |
| 14 | vector-spaces-basis-rank | Vector Spaces, Basis & Rank | Detect a mechanism in a truss |
| 15 | rotation-matrices | Matrix Transformations & Rotation Matrices | Direction cosine matrix tools |
| 16 | eigenvalues-eigenvectors | Eigenvalues & Eigenvectors | Principal axes of an inertia tensor |
| 17 | least-squares-qr | Orthogonality, Least Squares & QR | Fit a drag polar to wind-tunnel data |
| 18 | svd-and-norms | SVD & Matrix Norms | Conditioning of a structural stiffness matrix |

Lesson ids are `P01.NN`. Prerequisites: every lesson needs the previous Phase 01 lessons it uses and `P00.03` (NumPy), `P00.09` (SymPy) where symbolic work appears; name them in the metadata row. The learner has finished Phase 00: they can run Python, NumPy, SciPy, Matplotlib, SymPy and pytest, and they know how to state assumptions, report uncertainty and check an order of magnitude.

## What Phase 01 must achieve

Calculus and linear algebra at the level of a first-year engineering sequence (OpenStax Calculus 1–3, MIT 18.06), taught through aerospace problems and *always* paired with computation: every analytic result is checked numerically and every numerical result is checked against a limiting case. Derivations are shown in full in part 2 of each lesson. A learner who finishes can differentiate and integrate by hand, expand in Taylor series and bound the error, work in vectors and multivariable calculus, and solve, factor, decompose and condition linear systems, with a tested Python module for each.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-01/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-01/lesson-NN-<slug>` runs against `solutions/phase-01/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run `python scripts/check_lessons.py` (structure linter; it must report no ERROR for your lessons). The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (residuals, limiting cases, conservation, symmetry, convergence order, agreement between symbolic and numeric), not one magic number. 12–30 tests per lesson. No network access in tests. No files written outside a `tmp_path`. Tests run in under 10 s per lesson.

Linear-algebra lessons (13–18) implement the algorithm *from scratch* in NumPy (elimination, Gram–Schmidt, power iteration or Jacobi rotation, etc.) and the tests compare against `numpy.linalg` as the oracle; the lesson explains why we still call the library in practice.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who passed Bridge B1 (algebra, trig) and Phase 00. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs a **computational experiment** or **data-based investigation** in part 4 (no physical experiments in this phase), validates in part 5, and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: drag polars and glide range, airfoil sections, tapered beams and wing boxes, conic orbits, inertia tensors and principal axes, direction cosine matrices, vortices and circulation, control volumes, trusses, wind-tunnel fits, stiffness matrices. Numbers must be realistic and every number in the text must be reproduced by the solution code.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-01.md`. Where the map says GAP (01-15 Euler angles and DCMs; 01-18 matrix norms), write original instruction and say nothing about missing references. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only; type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders.
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, or templates.
- Do not call any `mcp__hearthbot__` tools, do not run git commands, and do not write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
