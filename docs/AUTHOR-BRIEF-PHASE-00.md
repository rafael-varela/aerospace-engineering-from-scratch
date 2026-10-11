# Author brief: Phase 00 lessons

This brief is what every lesson author in Phase 00 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-00.md`, `references/sources.yaml`, the Phase 00 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/` and `solutions/_template/` folders. Mirror the template exactly.

## Where files go

```
curriculum/phase-00/lesson-NN-<slug>/lesson.md
curriculum/phase-00/lesson-NN-<slug>/problems.md
curriculum/phase-00/lesson-NN-<slug>/starter/main.py
curriculum/phase-00/lesson-NN-<slug>/tests/test_main.py
solutions/phase-00/lesson-NN-<slug>/hints.md
solutions/phase-00/lesson-NN-<slug>/method.md
solutions/phase-00/lesson-NN-<slug>/worked-solutions.md
solutions/phase-00/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
```

Slugs (fixed, use exactly these):

| # | slug | title |
|---|---|---|
| 01 | dev-environment | Dev Environment for Engineers |
| 02 | git-and-collaboration | Git & Collaboration |
| 03 | numpy-engineering-calculator | NumPy as an Engineering Calculator |
| 04 | plotting-engineers-trust | Plotting That Engineers Trust |
| 05 | units-and-dimensions | Units & Dimensional Consistency |
| 06 | scipy-tour | SciPy Tour |
| 07 | testing-engineering-code | Testing Engineering Code |
| 08 | reading-data | Reading Data: CSV, HDF5, Telemetry Logs |
| 09 | symbolic-math-sympy | Symbolic Math with SymPy |
| 10 | writing-technical-reports | Writing Technical Reports |
| 11 | reading-a-technical-paper | Reading a Technical Paper |
| 12 | engineering-judgment | Engineering Judgment & Sanity Checks |
| 13 | safety-for-home-experiments | Safety & Standards for Home Experiments |
| 14 | measurement-calibration-uncertainty | Measurement, Calibration & Uncertainty |
| 15 | first-experiment | First Experiment: Predict, Measure, Compare |

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-00/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-00/lesson-NN-<slug>` runs against `solutions/phase-00/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish. The starter has the same function signatures and docstrings with `raise NotImplementedError` bodies, plus a `main()` that produces the lesson's output when complete. Tests assert properties (residuals, limiting cases, conservation, monotonicity, convergence, units round-trips), not one magic number. No network access in tests. No files written outside a `tmp_path`. Tests must run in under 10 s.

Lessons whose deliverable is not numerical (01, 02, 10, 11, 13) still get a small, honest `starter/main.py` and tests: for example an environment-check script that reports versions (01), a script that parses `git log` output from a fixture string (02), a report-structure linter for a Markdown file (10), a reading-notes template validator (11), a hazard-checklist completeness checker (13). Keep those small and meaningful; do not pad.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult beginner who passed Bridge B1 (algebra, trig) and may know no Python. Plain language, short sentences, define each term on first use. SI units.
- Every lesson states a concrete *engineering question* and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- Aerospace examples throughout: the ISA atmosphere, a drag polar, a rocket's mass depletion, a flight log, a wind-tunnel balance, the Mars Climate Orbiter.
- The practical activity is labelled with its mode. For 14 and 15 (🔬 physical experiments) the lesson links to `projects/p00-calibration/` and `projects/p00-first-experiment/` respectively, which another author writes; describe the activity in the lesson, and put the safety box first in part 4.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-00.md`. Where the map says GAP, write original instruction and say nothing about missing references. Never cite or require a paid resource. Never paste text from a source; write your own explanation and link.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or argument, one unfamiliar scenario, one that uses the lesson's code, one judgment question. Prefer parameterized problems (the learner's own measurements or a seed).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy/pint only; type hints; docstrings; `if __name__ == "__main__": main()`. Matplotlib must use the `Agg` backend in code that runs under tests, and `main()` saves figures to a path rather than calling `plt.show()`.
- Honesty rules: never claim the course grants a credential; never present simulation as physical measurement; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 250–500 lines. Depth over padding.

## Do not

- Do not edit files outside your assigned lesson folders.
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, or templates.
- Do not call any `mcp__hearthbot__` tools or write to `/mnt/project-files`.
- Do not install packages beyond `requirements.txt`.
- Do not add time estimates, grades, or badges.
