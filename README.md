# Aerospace Engineering from Scratch

> **Learn aerospace engineering from first principles. Design it on your computer. Build it at home. Test it against reality.**

A free, open-source curriculum that takes a motivated self-learner from algebra to the capabilities expected of a graduate from a top aerospace engineering program: aerodynamics, structures, propulsion, flight dynamics and control, guidance and navigation, orbital mechanics, spacecraft systems, avionics, testing, systems engineering, and design.

You need a computer, Python, basic hand tools, a modest and incremental budget, and no university lab. Every required step can be completed for free.

**Status:** curriculum roadmap and lesson standard are drafted; lesson content is being written phase by phase. See [ROADMAP.md](ROADMAP.md) for every phase and lesson and their status.

## How it works

Every lesson follows the same six-part arc:

1. **The engineering problem.** A real physical question or design trade, with constraints and success criteria.
2. **Understand and derive.** The physics, the governing equations, the assumptions and where they fail, worked analytic problems.
3. **Build the computational model.** Implement the method yourself in Python, test it, characterize its error, then compare with an established tool.
4. **Experiment at home or investigate real data.** A safe, low-cost measurement when the physics allows it; an authentic dataset or a defined simulation study otherwise.
5. **Validate and improve.** Compare analysis, computation, and measurement; quantify uncertainty; explain discrepancies; iterate.
6. **Defend your work.** Calculations, code, raw data, plots, interpretation, limitations, and a concise technical report.

You write the panel method before you open XFOIL, and the orbit propagator before you touch GMAT. You load a real beam before you trust your beam theory.

**Completion means competence, not finished pages.** Each major lesson assesses mathematical and physical reasoning, analytical problem solving, computational competence, and engineering interpretation separately. A passing code test alone does not pass a lesson.

## Learning paths

All paths share the Core (bridge modules as needed, then Phases 00–09: computing, mathematics, numerical methods, physics, mechanics, thermodynamics, fluids, signals and control). Then:

| Path | Adds | Ends with |
|---|---|---|
| Aerospace Engineer (full) | Phases 10–25 | One design capstone and one build-and-test capstone |
| Aeronautical Engineer | Aerodynamics, structures, air-breathing propulsion, performance, flight dynamics, GNC, avionics, testing, aircraft design | Aircraft design plus UAV or wind-tunnel build |
| Astronautical Engineer | Compressible flow, structures, rocket propulsion, GNC, astrodynamics, spacecraft systems, avionics, testing, mission design | Mission design plus thrust-stand or attitude-demonstrator build |
| Propulsion, GNC & Autonomy, Structures & Materials | Specialist subsets | A design plus a build capstone in the specialty |

Details, prerequisites, and the home-laboratory levels are in [ROADMAP.md](ROADMAP.md).

## Repository layout

```
curriculum/phase-XX/lesson-YY/   lesson.md, problems.md, starter/, tests/
solutions/phase-XX/lesson-YY/    hints.md, worked-solutions.md   (public, separate tree)
assessments/phase-XX/            exam.md, grading-rubric.md, parameterized-cases/
projects/<project-id>/           bom.yaml, safety.md, procedure.md, data/, report-template.md
references/                      sources.yaml (manifest), phase-XX.md (lesson-to-section maps)
scripts/                         check_paths.py, check_references.py (run in CI)
docs/                            design documents and reviews
```

Solutions are public and free, in their own tree so you can attempt first. The website hides them behind "reveal hint / method / solution" controls; the repository never does.

## Principles

- **Free, always.** Every required reading is free to access. Paid textbooks appear only in clearly marked optional sections and are never prerequisites.
- **Safe, always.** Every hands-on activity opens with a safety box. No lesson or capstone requires at-home combustion, rocket motor firing, pressure vessels, or autonomous outdoor flight; those happen only at qualified supervised clubs, as optional extensions.
- **Accessible, always.** Three equipment tiers per experiment, household and phone-sensor substitutions where scientifically sound, and no mandatory 3D printer, machine shop, or proprietary software.
- **Honest, always.** This roadmap's scope is mapped to ABET aerospace criteria and to top university programs. Finishing it does not confer a degree, a credential, or a license, and we never say otherwise. Time estimates are published only once measured with real learners.

## Contributing

Read [CONTRIBUTING.md](CONTRIBUTING.md) and [LESSON_TEMPLATE.md](LESSON_TEMPLATE.md). `python scripts/check_paths.py` and `python scripts/check_references.py` must pass on every change; CI runs them.

## License

Code (everything under `curriculum/**/starter`, `curriculum/**/tests`, `scripts/`, and other `.py` files) is licensed under the [MIT License](LICENSE). Text, problems, figures, and lesson content are licensed under [CC BY-SA 4.0](LICENSE-CONTENT.md). Third-party references keep their own licenses, recorded per source in `references/sources.yaml`.

Inspired by the structure of [AI Engineering from Scratch](https://github.com/rohitg00/ai-engineering-from-scratch).
