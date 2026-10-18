# Website

The landing page plus one page for every phase and every lesson, generated from the Markdown in `curriculum/`. There is one lesson template: every lesson page is built by the same code from the same files (`lesson.md`, `problems.md`, `starter/`, `tests/`), so a new lesson or phase gets its pages on the next deploy with no UI work.

```
npm ci            # once
npm run build     # writes dist/
npm run serve     # build, then serve dist/ at http://localhost:8000
```

Vercel deploys `main` with the root directory set to `site`; `vercel.json` tells it to run the build and publish `dist/`. The build reads `../curriculum`, so the Vercel project must keep **Include files outside the root directory in the Build Step** turned on (the default). CI runs the same build on every pull request, so a lesson that breaks the page build fails CI.

## What is generated

| URL | From |
|---|---|
| `/` | `index.html`, with every phase that exists in `curriculum/` linked and the status line filled in |
| `/curriculum/` | every `curriculum/phase-NN/README.md` |
| `/curriculum/phase-NN/` | that phase's `README.md`; its lesson table becomes lesson cards |
| `/curriculum/phase-NN/lesson-MM-slug/` | that lesson's `lesson.md` and `problems.md` |

A lesson page shows the title, the one-line tagline and the metadata table as a lesson specification (learning outcomes, prerequisites and the engineering question up front, the rest folded, and the safety note always visible when it is not "None"), then the lesson body, the problem set, a link to the solutions, the test commands, and links to the starter code and tests on GitHub. The sidebar lists the phase's lessons and the right-hand column lists the sections.

## Conventions the build relies on

These hold for every lesson today; `scripts/check_lessons.py` enforces the metadata rows and the section headings:

- `lesson.md` starts with `# PNN.MM  Title`, then a `>` tagline, then the metadata table, then `## 1. …` sections.
- The phase `README.md` starts with `# Phase NN: Title`, then an intro paragraph, then the lesson table with the header `| # | Lesson | Build it | Mode |`.
- Maths is written as `$…$` inline and `$$ … $$` on its own lines (optionally followed by an equation number, `$$ … $$   (1)`), rendered with KaTeX at build time. The build prints a warning for any formula KaTeX cannot parse.
- A lesson id in backticks (`` `P10.03` ``) becomes a link to that lesson. A file path in backticks becomes a link to the file on GitHub when the file exists. Relative Markdown links to lessons and phases go to their pages; links to anything else go to the file on GitHub.

## Files

- `index.html`: the landing page (content, styles and the interactive airfoil), hand-written.
- `config.js`: the project name and repository URL. Change the name there once; the landing page fills it in on load and the build uses it for every generated page. The fallback text in `index.html` (search for `data-site-name`) is only for readers without JavaScript, so update it too when the name is final.
- `generator/`: the build (`build.mjs` finds phases and lessons and writes the pages, `markdown.mjs` renders Markdown, maths and code, `templates.mjs` holds the page templates).
- `assets/lesson.css`, `assets/lesson.js`: styles and the small script (copy buttons, current-section highlight) for the generated pages. They share the landing page's colour tokens and type.

Content is drawn from `ROADMAP.md`, `README.md` and `curriculum/`. Keep the page's claims in step with them: scope comparable to a degree, competence demonstrated rather than conferred, no credential, and no total time estimate until pilots have measured one.
