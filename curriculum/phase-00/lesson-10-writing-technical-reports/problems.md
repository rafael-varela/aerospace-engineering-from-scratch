# P00.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer.

## 1. (Warm-up: significant figures from uncertainty)

A balance reports a mass of 0.247 83 kg and your calibration in `P00.14` style gives a standard uncertainty of 0.0012 kg. A stopwatch gives 12.4371 s with uncertainty 0.21 s. A pressure sensor gives 101 492 Pa with uncertainty 350 Pa. Write each quantity correctly, once in `value ± u unit` form and once in concise `value(u) unit` form. Answer form: six written quantities.

## 2. (Derivation: propagating uncertainty into a derived quantity)

The density of a machined aluminium cylinder is $\rho = m / (\pi r^2 h)$. Starting from the statement that for a product of powers $q = a^{\alpha} b^{\beta} c^{\gamma}$ with small, independent uncertainties the relative uncertainties combine as $(u_q/q)^2 = \sum (\text{exponent} \times u/x)^2$, derive the expression for $u_\rho / \rho$ in terms of $u_m/m$, $u_r/r$, and $u_h/h$. State where you use "small" and where you use "independent". Then evaluate it for $m = 0.2478 \pm 0.0012$ kg, $r = 15.02 \pm 0.02$ mm, $h = 130.1 \pm 0.1$ mm, and say which input dominates. Answer form: a derivation, a number with uncertainty in kg/m³, and one sentence naming the dominant term.

## 3. (Unfamiliar scenario: can the stated cause produce the stated effect?)

A report on a pendulum experiment measures a period 2.0 % longer than the small-angle prediction $T_0 = 2\pi\sqrt{L/g}$ and attributes it, in the Results section, to "air resistance". The release amplitude was 20°. The large-amplitude correction to the period is, to second order, $T \approx T_0\,(1 + \theta_0^2/16)$ with $\theta_0$ in radians. Compute the amplitude correction for 20°. Decide whether the report's attribution is supported, and state which section the sentence belongs in. Answer form: a percentage, a verdict with one sentence of reasoning, and the section name.

## 4. (Numerical, uses your code)

Write a 25-line Markdown report, from the template, about anything you measured this week. Make three deliberate mistakes: omit the Reproducibility heading, leave one figure without a caption, and write one bare quantity in the Results. Run `python starter/main.py yourfile.md` and paste the output. Then fix the three mistakes and paste the output again. Finally, add a *fourth* mistake that your linter does not catch, and say in one sentence why it slips through. Answer form: two terminal outputs and one sentence.

## 5. (Code extension: referenced figures)

Add a function `unreferenced_figures(text: str) -> list[int]` that returns the numbers of captioned figures that are never mentioned as "Figure n" (or "Fig. n") anywhere in the prose outside their own caption line. Write two tests first, one with a referenced figure and one without. Answer form: the function, the tests, and the pytest output. Do not modify the existing tests.

## 6. (Judgment: form versus truth)

Write, in at most 12 lines of Markdown, a report that passes every check in your linter and whose conclusion is nevertheless not supported by its own evidence. Then, in three sentences, say which human check would have caught it and whether it could be automated. Answer form: the report and the three sentences.

## 7. (Judgment: when the structure does not fit)

You are asked to write up a computational parameter sweep (no physical measurement): 200 runs of a drag-polar model with the Oswald efficiency factor varied from 0.6 to 0.95. Which of the nine sections change meaning, which become empty, and what would you put under "Apparatus, calibration, and safety"? Answer form: a short table of the nine sections with one line each.
