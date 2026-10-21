# Contributing

Thank you for helping build a free aerospace engineering education. This file says how to contribute so that a lesson written by anyone meets the same standard.

## Ground rules

1. **Free and original.** Required readings link to free sources. Explanations, problems, and solutions are written by you, in your own words. Never paste protected textbook text, instructor answer keys, or test banks.
2. **Safe.** Every hands-on activity has a safety box before the procedure: hazards, controls, supervision, operating envelope, stop criteria. No at-home combustion, motor firing, pressure vessels, high-speed unguarded rotors, or autonomous outdoor flight, ever, as a requirement.
3. **Honest.** Say what a model assumes and where it fails. Never claim a simulation is a physical test. Never claim the course grants a credential.
4. **Checked.** `python scripts/check_paths.py`, `python scripts/check_references.py` and `python scripts/check_lessons.py` pass, and `AEFS_IMPL=solution pytest <lesson dir>` passes for any lesson you touch.

## Writing a lesson

1. Pick a lesson marked ⬚ Planned in [ROADMAP.md](ROADMAP.md) and open an issue titled `Lesson Pxx.yy: <title>` so nobody duplicates it.
2. Copy `curriculum/_template/lesson-template/` to `curriculum/phase-XX/lesson-YY-<slug>/` and `solutions/_template/` to `solutions/phase-XX/lesson-YY-<slug>/`.
3. Follow [LESSON_TEMPLATE.md](LESSON_TEMPLATE.md) exactly: the metadata block, then the six parts in order.
   Draw diagrams as SVG figures in the house style of [docs/FIGURES.md](docs/FIGURES.md), not as text art in code blocks.
4. Your "Build" code lives in `starter/` with tests in `tests/`. Tests check residuals, conservation laws, limiting cases, and convergence, not just one expected number.
5. Add every reference you cite to `references/sources.yaml` (or reuse an existing `id`) and to `references/phase-XX.md` with exact section numbers.
6. Mark the lesson 🚧 in ROADMAP.md in your PR; a reviewer marks it ✅ after the review checklist below passes.

## Review checklist

- [ ] Metadata block complete; learning outcomes are observable skills.
- [ ] Engineering problem is concrete and states success criteria.
- [ ] Derivation names assumptions, units, boundary conditions, and failure cases.
- [ ] Code runs from a clean environment; tests pass; numerical error is characterized.
- [ ] Practical activity labelled (physical / computational / data / design-build), with equipment tiers and safety box where physical.
- [ ] Validation compares at least two of analysis, computation, measurement, and discusses discrepancies.
- [ ] Adversarial self-checks included.
- [ ] Hints → method → full worked solution exist in `solutions/`, in that order of disclosure.
- [ ] Every reference is free to access, mapped to an exact section, and in `sources.yaml` with license and verification date.
- [ ] No paid resource is required anywhere.

## Writing a project (home experiment)

Copy `projects/_template/`. The BOM must have three tiers (minimum-cost, recommended, enhanced), vendor-neutral search terms, substitutions, a reuse map to earlier projects, and a dated price check. The procedure must demand a prediction before the measurement.

## Style

- Plain language. Short sentences. Define a term the first time it appears.
- SI units everywhere; show the imperial equivalent only when industry practice demands it (altitudes in feet, airspeed in knots) and say so.
- Code: Python 3.11+, NumPy, SciPy, Matplotlib; type hints; a docstring on every public function; no notebooks as the primary artifact (notebooks may accompany a lesson).
- Commit messages describe the change in one plain sentence.

## Code of conduct

See [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).
