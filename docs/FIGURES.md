# Figures

Diagrams in lessons, problem sets, solutions and project kits are SVG files drawn in one house style. They are not text art in code blocks. One SVG renders on the website, on GitHub and in any Markdown viewer, scales cleanly on a phone, and follows the reader's light or dark mode.

## Where they go and how they are referenced

- Put the file in a `figures/` folder next to the Markdown that uses it, for example `curriculum/phase-10/lesson-01-forces-moments-coefficients/figures/body-wind-axes.svg`. Name it for what it shows, in lower-case kebab case.
- Reference it with an image line, followed on the next line by an italic caption that starts with `Figure N.` (the convention taught in `P00.10`). Number figures in order within each Markdown file.

```markdown
![Airfoil with body axes x along the chord and z normal to it; the stream arrives at angle α; lift is perpendicular to the stream and drag along it.](figures/body-wind-axes.svg)
*Figure 1. Body axes and wind axes. The balance measures N and A; lift and drag are the same force resolved perpendicular to and along the stream.*
```

- The alt text says what the figure shows, in a sentence or two, for a reader who cannot see it. The caption says what to notice. Neither repeats the other word for word.
- The site build fails if a Markdown file points at an image that does not exist.

## House style

Start every figure from [`figure-style.svg`](figure-style.svg). Copy its `<style>` and `<defs>` block unchanged and draw with its classes. Never write a colour into an element: the classes switch between light and dark palettes.

| Use | Class | Arrowhead |
|---|---|---|
| Bodies, structures, outlines | `line`, or `body` for a filled shape | `url(#ah-ink)` |
| Axes, coordinates, flow and streamlines | `axis`, `flow`; labels `t-blue` | `url(#ah-blue)` |
| Forces, moments, loads, the quantity the figure is about | `force`; labels `t-orange` | `url(#ah-orange)` |
| Dimensions, construction lines, leaders | `thin`, add `dash` for hidden or reference lines | `url(#ah-ink2)` |
| Ground, walls, clamps | `ground` (hatched) with a `line` edge | |
| A second filled region (fluid, a highlighted part) | `body-2` | |

- Text: labels in the default sans (`text`), symbols and variables in italic serif (`sym`), secondary notes in `note`, panel headings in `head`. Use real Unicode for Greek letters, primes, subscripts where they read well (`α`, `Γ`, `x′`, `V∞`), and keep the symbols exactly as the lesson text writes them.
- Size: `viewBox` 720 units wide (narrower is fine for a small figure), height as needed, with matching `width` and `height` attributes. Keep labels at 12.5 px or larger at that width.
- Give the root `role="img"` and a `<title>` with the same text as the alt text.
- Draw what the text uses: the axes, sign conventions, labels and numbers of the lesson. A figure that disagrees with the equations next to it is a bug. Exaggerate small angles or deflections when it helps, and say so in the caption.
- Qualitative sketches (a phase portrait, a stress–strain curve, a T–s diagram) are drawn smoothly from the shape the lesson describes, with labelled axes and no invented tick values. Plots of data or computed results belong in the lesson's own code output, not in a hand-drawn figure.

## What stays as text

Program output, tables of numbers, matrices, directory trees, terminal sessions and Markdown or code examples stay in code blocks. They are text, and a reader may want to copy them.
