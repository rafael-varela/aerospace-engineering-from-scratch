# P00.10  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1.

**Assumptions.** Uncertainties are standard uncertainties; the rounding convention is one significant figure, or two when the leading digit is 1 (so that rounding does not throw away more than about 10 % of the uncertainty).
**Derivation.** Mass: $u = 0.0012$ kg; leading digit 1, keep two figures: 0.0012 kg; value to four decimals: 0.2478 kg. Time: $u = 0.21$ s → 0.2 s; value to one decimal: 12.4 s. Pressure: $u = 350$ Pa → 400 Pa (one figure) or 350 Pa (two figures are acceptable here; either is defensible, be consistent); value to the hundreds: 101 500 Pa.
**Result.**
- $0.2478 \pm 0.0012$ kg, or $0.2478(12)$ kg
- $12.4 \pm 0.2$ s, or $12.4(2)$ s
- $(1.015 \pm 0.004) \times 10^5$ Pa, or $1.015(4) \times 10^5$ Pa; equivalently $101.5 \pm 0.4$ kPa, $101.5(4)$ kPa
**Check.** Each written value has its last digit in the same decimal place as the uncertainty's last digit. In concise form, the parenthesised digits occupy exactly the positions of the last digits shown.
**Common mistakes.** Keeping all six digits of 0.247 83 "because the balance showed them"; writing `12.4371 ± 0.2 s`; rounding 101 492 to 101 000 (that is the thousands place, coarser than the uncertainty).

## 2.

**Assumptions.** Uncertainties are small relative to values (so a first-order Taylor expansion of $\ln\rho$ is accurate) and uncorrelated (so cross terms vanish). Cylinder is a right circular cylinder (shape error ignored).
**Derivation.** $\rho = m\,r^{-2}\,h^{-1}\,\pi^{-1}$. Take the logarithm: $\ln\rho = \ln m - 2\ln r - \ln h - \ln\pi$. Differentiate: $\delta\rho/\rho = \delta m/m - 2\,\delta r/r - \delta h/h$ (this is where "small" is used: the first-order term is kept). Square and take expectations; the cross terms are covariances, zero for independent inputs (this is where "independent" is used). Hence

$$\frac{u_\rho}{\rho} = \sqrt{\left(\frac{u_m}{m}\right)^2 + \left(2\,\frac{u_r}{r}\right)^2 + \left(\frac{u_h}{h}\right)^2}.$$

Numbers: $u_m/m = 0.0012/0.2478 = 4.84 \times 10^{-3}$; $2u_r/r = 2 \times 0.02/15.02 = 2.66 \times 10^{-3}$; $u_h/h = 0.1/130.1 = 0.77 \times 10^{-3}$. Sum of squares: $23.4 + 7.09 + 0.59 = 31.1 \times 10^{-6}$; root: $5.58 \times 10^{-3}$. Value: $r = 0.01502$ m, $h = 0.1301$ m, volume $= \pi (0.01502)^2 (0.1301) = 9.221 \times 10^{-5}$ m³, $\rho = 0.2478 / 9.221 \times 10^{-5} = 2687$ kg/m³. $u_\rho = 2687 \times 5.58 \times 10^{-3} = 15$ kg/m³.
**Result.** $\rho = 2687 \pm 15$ kg/m³ (or $2.687(15) \times 10^3$ kg/m³). The mass term dominates (75 % of the variance).
**Check.** Dimensions: kg / (m² · m) = kg/m³. Plausibility: aluminium alloys are 2650–2800 kg/m³; the value sits inside. Limiting case: set $u_m = 0$ and the relative uncertainty drops to $2.77 \times 10^{-3}$, confirming the mass term carried most of it.
**Common mistakes.** Forgetting the factor 2 on the radius term; adding relative uncertainties linearly (gives $8.3 \times 10^{-3}$, an overestimate for independent inputs); mixing mm and m in the volume.

## 3.

**Assumptions.** Simple pendulum with a dense bob; the second-order amplitude series is adequate at 20°.
**Derivation.** $\theta_0 = 20° = 0.349$ rad. $\theta_0^2/16 = 0.1218/16 = 7.6 \times 10^{-3}$, i.e. the period is 0.76 % longer than $T_0$ from amplitude alone. Air drag on a dense bob removes energy and reduces the amplitude swing by swing; its effect on the period is through the amplitude term above and is therefore smaller than 0.76 % and *shrinks* over the run. Drag cannot produce a 2 % period increase here.
**Result.** Amplitude correction: $+0.76\,\%$. Verdict: the attribution to air resistance is not supported; it cannot produce the magnitude. About 0.8 % of the 2 % is explained by finite amplitude; the remaining 1.2 % needs another cause (timing reaction bias, a length measured to the top of the bob rather than its centre of mass, or string stretch). The sentence belongs in the Discussion, not the Results, and should appear in a ranked list of causes with the size each can produce.
**Check.** Limiting case: $\theta_0 \to 0$ gives $T \to T_0$, as the small-angle formula requires. Order of magnitude: measuring $L$ to the top of a 4 cm bob on a 1 m string gives a 2 % length error, hence a 1 % period error, so that candidate alone is the right size.
**Common mistakes.** Squaring 20 (degrees) instead of 0.349 (radians); treating "air resistance" as an explanation without estimating its size; placing interpretation in the Results.

## 4.

**Assumptions.** The learner's report is their own; the outputs below illustrate the expected form.
**Derivation.** First run, typical output:

```
Linting report.md
 - missing heading: reproducibility
 - line 31: figure 'figures/trend.png' has no 'Figure <n>' caption
 - results, line 28: '2.3 kg' has no uncertainty
3 problem(s) found.
```

Second run: `OK: all required sections present, ...`. Fourth mistake examples that slip through: `12 kn` (unit not in the list), `9.8 m/s²` written in the Summary (not the Results section), a captioned figure never referred to in the prose, a Prediction section written after the data were taken (no code can detect this; only a dated commit can).
**Result.** Two outputs and one sentence, as above.
**Check.** The exit status should be 1 on the first run and 0 on the second (`echo $?` in a shell).
**Common mistakes.** Treating a passing linter as a passing report.

## 5.

**Assumptions.** References appear as "Figure n", "Fig. n", or "Fig n" in prose; caption lines are excluded.
**Derivation.**

```python
_REF_RE = re.compile(r"\bFig(?:ure|\.)?\s*(\d+)", re.IGNORECASE)

def unreferenced_figures(text: str) -> list[int]:
    figures = find_figures(text)
    caption_lines = set()
    lines = text.splitlines()
    for f in figures:
        for n in range(f.line, len(lines)):
            if lines[n].strip():
                caption_lines.add(n + 1)
                break
    referenced = {
        int(m.group(1))
        for i, line in enumerate(lines, start=1)
        if i not in caption_lines
        for m in _REF_RE.finditer(line)
    }
    return [f.number for f in figures if f.number is not None and f.number not in referenced]
```

Tests: a two-line report with `![a](a.png)`, `*Figure 1. x*`, and the prose "see Figure 1" returns `[]`; the same without the prose returns `[1]`.
**Result.** Function and two tests, pytest green.
**Check.** Limiting case: a report with no figures returns `[]`. A report whose only mention is inside the caption returns `[n]`, confirming captions are excluded.
**Common mistakes.** Counting the caption itself as a reference; matching "Figure" in a sentence like "a figure of merit" (the `\d+` requirement prevents this).

## 6.

**Assumptions.** None beyond the linter as written.
**Derivation.** Example (12 lines):

```markdown
# Spring constant: Technical report
## Summary
k = 20.1 ± 0.3 N/m, in agreement with prediction.
## 1. Objective and requirements
Predict k to 5 %.
## 2. Theory and assumptions
Hooke's law, F = k x.
## 3. Prediction (recorded before measurement)
k = 10.0 ± 0.5 N/m from the manufacturer's data sheet.
## 4. Apparatus, calibration, and safety
Ruler, known masses.
## 5. Results
k = 20.1 ± 0.3 N/m (Figure 1).
![F vs x](f.png)
*Figure 1. Force against extension.*
## 6. Discussion
Agreement is good.
## 7. Conclusions and redesign
The data sheet is confirmed.
## 8. Reproducibility
Commit abc1234.
```

Prediction $10.0 \pm 0.5$ and result $20.1 \pm 0.3$ differ by 10.1 N/m, about 17 combined standard uncertainties. The Conclusion is false; every check passes.
**Result.** The human check is: compare the Prediction and Results numbers against their combined uncertainty and check that the Discussion's verdict matches. It can be partly automated when both are in `value ± u unit` form with the same unit: extract both, compute $|a-b|/\sqrt{u_a^2+u_b^2}$, and flag a Discussion containing "agree" when that exceeds 2 or 3. Whether the Conclusion *follows* in general cannot be automated.
**Check.** The linter's own test fixture `GOOD_REPORT` passes for the same reason: nothing in it checks meaning.
**Common mistakes.** Breaking the form (omitting a heading) rather than the content; the problem asks for a report that passes.

## 7.

| Section | For a computational sweep |
|---|---|
| Summary | Unchanged: question, method, headline trend with the spread across runs, agreement with expectation, main limitation. |
| 1. Objective | Unchanged: what trend or threshold you set out to find and the criterion for "found". |
| 2. Theory | Unchanged, plus the model's known failure cases (the Oswald factor is an empirical fit; outside 0.6–0.95 the polar is not meaningful). |
| 3. Prediction | Still exists: the trend expected before running (L/D max should scale as $\sqrt{e}$), with a rough number. |
| 4. Apparatus, calibration, safety | Software environment (Python, NumPy versions), solver settings and tolerances, random seeds; "calibration" becomes a verification run against a hand-computed case; safety: none beyond computer use. |
| 5. Results | Unchanged in meaning; "uncertainty" becomes numerical tolerance and sensitivity to step size, which must still be stated. |
| 6. Discussion | Unchanged; ranked causes of any deviation from the expected trend; the follow-up run. |
| 7. Conclusions and redesign | Unchanged; redesign refers to the model or the sweep design. |
| 8. Reproducibility | Unchanged and more important: commit hash plus the exact command line reproduce everything. |

**Check.** No section becomes empty; the ones that change most are 4 and the meaning of "uncertainty" in 5.
**Common mistakes.** Deleting the Prediction because "it is only a computation"; writing "N/A" for Reproducibility when it is the one section a computation can satisfy perfectly.

## Check-yourself answers

1. Discussion. It interprets a result (assigns a cause); Results state only what was observed.
2. $1.000 \pm 0.0005$ m if you take half a division, commonly written $1.000 \pm 0.001$ m to include alignment at both ends; Type B (from the instrument's specification, not from repeated readings).
3. $u_g/g = \sqrt{0.002^2 + (2 \times 0.0025)^2} = \sqrt{4 + 25}\times 10^{-3} = 5.4 \times 10^{-3}$, i.e. 0.54 %. Improve the timing: its term (squared) is six times the length's.
4. $2.010 \pm 0.005$ s.
5. So that the comparison in the Discussion is a test of the model rather than a restatement of the result. An undated or edited prediction cannot be shown to have preceded the data.
6. Not a bug. The rule is scoped to Results because the Theory section legitimately quotes nominal values (a 300 K reference temperature) that are inputs to a model, not measurements.
7. Alt text is read aloud by screen readers and shown when the image cannot load; it describes the image. A caption labels and numbers the figure and explains axes, units, and error bars to a sighted reader; the linter looks for the caption on the line after the image.
8. `figure_numbering_issues`. For the third captioned figure it returns a message of the form "line N: expected Figure 3, found Figure 4".
