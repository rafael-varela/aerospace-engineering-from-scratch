# Lesson Template

Every lesson in `curriculum/` follows this standard. Copy `curriculum/_template/lesson-template/` and `solutions/_template/` to start. The point of the standard is that a lesson written by anyone meets the same bar: observable outcomes, real engineering evidence, and honest validation.

## Folder structure

```
curriculum/phase-XX/lesson-YY-<slug>/
├── lesson.md          the lesson (metadata block + six parts)
├── problems.md        problem set (prompts only; no answers here)
├── starter/
│   └── main.py        scaffold the learner completes; runs end to end when done
├── tests/
│   └── test_main.py   public tests: residuals, limits, conservation, convergence
└── notebook/          optional companion notebook; never the primary artifact

solutions/phase-XX/lesson-YY-<slug>/
├── hints.md           one targeted hint per problem
├── method.md          method outline per problem (no final numbers)
└── worked-solutions.md full solutions: assumptions, units, derivation, checks, common mistakes

assessments/phase-XX/   exam.md, grading-rubric.md, parameterized-cases/   (one set per phase)
projects/<project-id>/  bom.yaml, safety.md, procedure.md, data/, report-template.md
```

Solutions are public and free. They live in a separate tree so the learner attempts first; the website hides them behind reveal controls, the repository never does.

## `lesson.md` format

```markdown
# Pxx.yy  <Lesson Title>

> <One-line motto: the idea that should stick.>

| | |
|---|---|
| **Learning outcomes** | 3–5 observable skills: "derive…", "implement and test…", "measure…", "decide whether…". Not "understand X". |
| **Prerequisites** | Lesson IDs (`P01.13`, `P03.05`) plus skill checks ("can solve a 3×3 linear system by hand"). |
| **Engineering question** | The prediction or decision the lesson makes. One sentence. |
| **Assumptions and limits** | The model's conditions of validity and its known failure cases. |
| **Analytical work** | What is derived and calculated by hand. |
| **Computational deliverable** | What the code does and what its tests prove (accuracy, convergence, conservation). |
| **Practical activity** | Mode: **physical experiment** / **computational experiment** / **data-based investigation** / **design-build**. One line on what is done. |
| **Equipment** | Required vs optional, substitutes, tier (none / Level 1 / Level 2 / Level 3), link to the project's `bom.yaml` if any. |
| **Safety** | "None beyond normal computer use" or a link to the project's `safety.md`. |
| **Evidence** | What the learner submits: data, plots, uncertainty, report. |
| **Assessment** | Link to rubric; which of the four competencies are assessed; pass thresholds. |
| **References** | Exact sections, each tagged required-free / supplementary-free / optional-paid. |
| **Extensions** | Optional ★ or specialist challenges. |
| **Time** | Omitted until measured in pilots. Write "unmeasured". |

## 1. The engineering problem

Two or three paragraphs. A real physical question, design trade, or observable phenomenon. State constraints and success criteria. Make it concrete: what goes wrong if you cannot do this.

## 2. Understand and derive

The physics and the governing equations, derived, with every assumption named where it is made. Units and dimensions checked. Boundary and initial conditions. Where the model stops being valid. One or two representative problems worked fully. Diagrams (SVG or ASCII) wherever a picture carries the idea.

## 3. Build the computational model

Step-by-step implementation in Python from the equations in part 2. Each step is a runnable code block with a short explanation. Then:

- **Verify:** test against an analytic case, a limiting case, and a conservation law; show the convergence or error plot.
- **Compare:** run the established tool or library for the same case and explain any difference.

## 4. Experiment at home, or investigate real data

When the physics allows a safe, affordable measurement: the apparatus (by tier), the calibration, the prediction *written down before measuring*, the procedure, the raw data format. When it does not: an authentic dataset (source, license, what was measured and how) or a defined computational experiment (sensitivity study, parameter sweep, method comparison). Label the mode. Physical activities link to `projects/<id>/safety.md` and put the safety box first.

## 5. Validate and improve

Put analysis, computation, and measurement side by side. Quantify uncertainty. Explain every discrepancy; "experimental error" is not an explanation. Say what you would change and, where feasible, change it and re-run.

## 6. Defend your work

What to submit: calculations, code, raw data, plots, interpretation, limitations, and a short technical report (template in `projects/_template/report-template.md`). Then the adversarial self-checks:

- What if the load doubles? What if velocity goes to zero? What are the limiting cases?
- Which assumption fails first in the real world?
- What measurement could falsify the model?
- Does the result survive a change in mesh, time step, solver, or input?
- Which decision would change if the measurement uncertainty were twice as large?

## Check yourself

5–8 quick questions on terminology, units, limiting cases, and intuition (auto-gradable; answers in `solutions/`).
```

## `problems.md` format

Numbered problems, graded easy → hard. Each states givens with units, what to find, and what form the answer takes (number with uncertainty, derivation, plot, short argument). Include at least one derivation, one unfamiliar scenario, one numerical problem that exercises the lesson's code, and one judgment question about model limits. Prefer parameterized problems (the learner's own data or a seed) so answers cannot simply be copied.

## Code standards

- Python 3.11+, NumPy, SciPy, Matplotlib, SymPy, pint. No proprietary dependencies.
- `starter/main.py` has the function signatures, docstrings, and `raise NotImplementedError` bodies; a completed version runs from the command line and produces the lesson's plots.
- `tests/test_main.py` asserts physics, not magic numbers: residuals below tolerance, limiting cases, conservation, and a convergence rate. Hidden validation cases for assessment live in `assessments/phase-XX/parameterized-cases/`.

## Assessment standard

Four competencies, assessed separately and reported separately:

1. Mathematical and physical reasoning
2. Analytical problem solving
3. Computational competence
4. Engineering interpretation

A lesson is passed when all four meet the rubric threshold. Passing the code tests alone never passes the lesson.

## Reference standard

Each reference entry: `id` from `references/sources.yaml`, exact chapter/section number and title, what to extract, access tier. Required readings are free. Paid books only under an "Optional textbooks for deeper study" heading, with the specific value they add and non-purchase options (library, used copy).
