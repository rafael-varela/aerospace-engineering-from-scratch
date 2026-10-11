# P00.11  Reading a Technical Paper

> A paper is not read front to back. It is interrogated: what did they ask, what did they do, what did they see, and where does it stop being true?

| | |
|---|---|
| **Learning outcomes** | 1. Extract from a technical report its question, method, key figure, headline result, and limitations, and write them into a seven-field reading-notes file. 2. Decode a NACA four-digit designation into maximum camber, camber position, and thickness, and draw the section it describes. 3. Read a lift curve (C_l against α) and a drag polar (C_d against C_l) to obtain the lift-curve slope, the zero-lift angle, the minimum drag, and the maximum lift-to-drag ratio, with a reading uncertainty. 4. Implement and test the designation parser, the section geometry, the curve readers, and a reading-notes validator. 5. Decide whether a 1945 wind-tunnel dataset applies to a given modern problem, and say why. |
| **Prerequisites** | `P00.04` (reading and making plots with labelled axes), `P00.10` (report structure); skill check: can read a value off a graph between grid lines and state how precisely. |
| **Engineering question** | From NACA Report 824 alone, can you state what a NACA 2412 wing section looks like and what lift and drag it produces at a given angle of attack and Reynolds number, and where that data stops applying? |
| **Assumptions and limits** | The four-digit geometry formulas describe the NACA four-digit family exactly, and nothing else (the 5-digit and 6-series families use different camber lines). The lift and drag data in the report are two-dimensional (infinite wing), low-speed (Mach well below the drag rise), at Reynolds numbers of 3 to 9 million, in a low-turbulence tunnel; applying them to a finite wing, a small model, or a transonic airliner needs corrections you do not yet have. Reading values off a scanned plot carries an uncertainty of roughly one minor grid division. |
| **Analytical work** | Deriving the two-parabola camber line from four conditions; the thin-airfoil lift-curve slope of 2π per radian as a check value; the tangent construction for maximum L/D on a drag polar. |
| **Computational deliverable** | `starter/main.py`: `parse_naca4` (tested on 2412, 0012, 4415, 6409 and malformed input), the camber line and thickness distribution (tested for peak position, peak value, continuity, and derivative consistency), a lift-curve fit that recovers a known slope, a drag-polar reader that finds the analytic L/D optimum, and `validate_reading_notes`, which reports missing or empty fields. |
| **Practical activity** | **data-based investigation**: download NACA Report 824 from NTRS, locate the section characteristics plots for NACA 2412 and 0012, digitise about eight points from each curve by eye, fit them with your code, and write the reading notes. |
| **Equipment** | None beyond a computer and a PDF reader. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `reading-notes.md` that passes your validator; your digitised CSV; the fitted slope and zero-lift angle with a reading uncertainty; the plotted section for 2412 from your code; passing tests. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 1 (the camber-line derivation), 3 (the code), and 4 (applicability judgment). Pass: tests green, notes complete, fitted slope within reading uncertainty of the report's curve, and problem 6 argued with a computed Reynolds number. |
| **References** | `@ntrs-naca-tr-824-summary-of-airfoil-data`, NACA Report 824 "Summary of Airfoil Data": the summary and introduction (the first pages), the description of the four-digit series, and the section-characteristics plots for NACA 0012 and NACA 2412 (required-free; download the PDF from the NTRS citation page, the inline preview is unavailable). `@leishman-intro-aerospace-flight-vehicles` "Airfoil Shapes" (supplementary-free; a modern overview that cites the report). |
| **Extensions** | ★ Add the NACA five-digit camber line (230xx series) to the parser. ★ Digitise the moment-coefficient curve too and locate the aerodynamic centre. |
| **Time** | unmeasured |

## 1. The engineering problem

Every wing you will design in Phase 10 and Phase 24 starts from a wing *section*: the two-dimensional shape you would see if you sliced the wing with a plane parallel to the flight direction. Before you can choose one, you need to know how much lift and drag it makes at each angle to the oncoming air. Measuring that is a wind-tunnel programme. Between 1930 and 1945 the National Advisory Committee for Aeronautics (NACA, the predecessor of NASA) ran exactly such a programme on hundreds of systematically shaped sections, and in 1945 Abbott, von Doenhoff, and Stivers summarised it in Report 824. It is a public-domain document, it is still cited in design textbooks, and the four-digit sections it describes (the NACA 2412 on countless light aircraft, the symmetric 0012 on tails and helicopter blades) are still flying.

So the first "paper" you read in this course is a 1945 government technical report, for three reasons. Its data are real measurements, with the Reynolds number and the surface condition stated. Its figures are the kind you will read for the rest of your career: a lift curve and a drag polar side by side. And it is written for an engineer who needs numbers, not for a reader who wants to be told what to think, which means it rewards the reading method this lesson teaches and punishes skimming.

The success criterion is concrete. At the end you hand in a one-page reading-notes file, a drawing of the NACA 2412 produced by your own code from the four digits alone, and a lift-curve slope and zero-lift angle for the 2412 that you read off the report's plot, with an uncertainty, and that agree with what thin-airfoil theory says they should be. If you cannot do this, you will not be able to use any reference data in later phases without someone else telling you what it says.

## 2. Understand and derive

### 2.1 How to read a technical document

Read in three passes, and write something down after each.

**Pass 1, five minutes: what is this?** Read the title, the abstract or summary, the section headings, and the captions of every figure. Do not read the body. Write one sentence: the question the authors are answering. For Report 824 you should be able to write "a compilation of measured section lift, drag, and moment characteristics for the NACA airfoil families, with a description of how the families are defined" before reading any paragraph of the body.

**Pass 2, thirty minutes: how did they do it, and what did they see?** Read the method (for a wind-tunnel report: which tunnel, what model size, what Reynolds numbers, what surface condition), then study the figures properly. For each figure you care about: what is on each axis, with units; what each curve or symbol represents; what one grid division is worth. Then find the headline result, the one number or trend the whole document exists to deliver. Write the method in two sentences and the result in one.

**Pass 3, as long as it takes: where does it stop being true?** Every measurement was made under conditions. List them. Two-dimensional or three-dimensional? Which Mach number? Which Reynolds number? Smooth or rough? Tunnel-wall corrections applied or not? Then ask the question an engineer always asks: *does my problem sit inside those conditions?* Write the limitations in three sentences, and write down every question the document raised that it did not answer. Those questions are the most valuable thing you take away.

The reading-notes template for this course has seven fields, one per pass item:

```markdown
# Reading notes: <short title>

## Citation
Authors (year), title, report number, where you got it.

## Question
What were they trying to find out? One or two sentences.

## Method
How? Apparatus, conditions, what was measured. Two or three sentences.

## Key figure
Which figure matters most, and what is on its axes (with units).

## Result
The headline number or trend, with the uncertainty or the spread they report.

## Limitations
The conditions the result holds under, and what it does not cover.

## My questions
What the document raised and did not answer.
```

Your validator will check that each field exists and has something in it. It cannot check whether what you wrote is true. That is your job.

### 2.2 Abstract, summary, and the trap of reading only one

Report 824 opens with a summary, and NTRS shows an abstract on its citation page. They are not the same thing. An **abstract** is written to help you decide whether to read the document; it says what the document contains. A **summary** is written for someone who will not read the rest; it states the findings. A good reading note never cites either as the result. The result lives in the figures and tables, and your note points at the specific figure.

### 2.3 The NACA four-digit designation

Three words first. The **chord** is the straight line from the leading edge (the nose) to the trailing edge (the tail); its length *c* is the reference length, and every other length is quoted as a fraction of it. The **mean camber line** is the curve halfway between the upper and lower surfaces; its maximum distance from the chord is the **camber**. The **thickness** is the distance between the upper and lower surfaces measured perpendicular to the camber line; its maximum is the section's thickness.

A NACA four-digit designation `MPTT` is a recipe:

| Digit | Meaning | NACA 2412 | NACA 0012 | NACA 4415 |
|---|---|---|---|---|
| M | maximum camber, in hundredths of chord: m = M/100 | 0.02 | 0 | 0.04 |
| P | position of maximum camber from the leading edge, in tenths of chord: p = P/10 | 0.4 | 0 (unused) | 0.4 |
| TT | maximum thickness, in hundredths of chord: t = TT/100 | 0.12 | 0.12 | 0.15 |

So a NACA 2412 has 2 % camber at 40 % chord and is 12 % thick. A NACA 0012 has no camber (it is symmetric; the second digit has no meaning) and the same thickness. A NACA 4415 is more cambered and thicker. The parser you will write turns the string into these three numbers and refuses strings that are not four digits, that have zero thickness, or that claim camber with no position (`2012` is not a valid section).

### 2.4 Deriving the camber line

The NACA four-digit camber line is two parabolic arcs. Writing $y_c(x)$ for the camber line with $x$ and $y_c$ in chord fractions, the designers required:

1. $y_c(0) = 0$ at the leading edge;
2. $y_c(p) = m$: maximum camber $m$ at $x = p$;
3. $y_c'(p) = 0$: the maximum is a smooth peak;
4. $y_c(1) = 0$ at the trailing edge.

A parabola has three coefficients, so one parabola can satisfy three conditions. Use one parabola on $0 \le x < p$ for conditions 1–3 and another on $p \le x \le 1$ for conditions 2–4. For the front arc, write $y_c = a x^2 + b x$ (condition 1 kills the constant). Condition 3 gives $2ap + b = 0$; condition 2 gives $ap^2 + bp = m$. Substituting $b = -2ap$: $ap^2 - 2ap^2 = -ap^2 = m$, so $a = -m/p^2$ and $b = 2m/p$:

$$ y_c = \frac{m}{p^2}\left(2px - x^2\right), \qquad 0 \le x < p. \qquad (1) $$

For the rear arc, the same algebra with the trailing-edge condition (problem 2 asks you to do it) gives

$$ y_c = \frac{m}{(1-p)^2}\left[(1 - 2p) + 2px - x^2\right], \qquad p \le x \le 1. \qquad (2) $$

Check equation (2) at $x = p$: $(1 - 2p) + 2p^2 - p^2 = 1 - 2p + p^2 = (1-p)^2$, so $y_c = m$. Check at $x = 1$: $(1 - 2p) + 2p - 1 = 0$. The slopes are $y_c' = \frac{2m}{p^2}(p - x)$ in front and $\frac{2m}{(1-p)^2}(p - x)$ behind; both vanish at $x = p$, so the line is smooth there. Both are zero everywhere when $m = 0$: a symmetric section has a straight camber line, as it must.

### 2.5 The thickness distribution

The four-digit thickness is one polynomial for the whole family, scaled by $t$:

$$ y_t = 5t\left(0.2969\sqrt{x} - 0.1260\,x - 0.3516\,x^2 + 0.2843\,x^3 - 0.1015\,x^4\right). \qquad (3) $$

$y_t$ is the *half*-thickness: the surface sits $y_t$ above and below the camber line. The $\sqrt{x}$ term gives the rounded nose (infinite slope at $x = 0$). The maximum of the bracket is 0.1000 at $x \approx 0.30$, so the maximum thickness is $2 \times 5t \times 0.1 = t$, as the designation promises. At $x = 1$ the bracket is $0.0021$, not zero: the four-digit trailing edge has a small finite thickness, which is why real drawings often modify the last coefficient. The surface points are laid off perpendicular to the camber line: with $\theta = \arctan y_c'$, the upper surface is $(x - y_t\sin\theta,\ y_c + y_t\cos\theta)$ and the lower is $(x + y_t\sin\theta,\ y_c - y_t\cos\theta)$. For a symmetric section $\theta = 0$ and the surfaces are simply $\pm y_t$.

```
   y/c
    ^           NACA 2412 (schematic)
0.08|      .--------------.._
    |    .'                  ''--.._          upper surface
0.00|---'===========================''--.___  mean camber line (dashed in your plot)
    |    '.,                           __.-'  lower surface
-0.06|       ''-----------------------''
    +------------------------------------------> x/c
    0          0.3  0.4                      1.0
               max t   max camber
```

### 2.6 Reading a lift curve

The report's section-characteristics plots put the **lift coefficient** $C_l = L' / (\tfrac12 \rho V^2 c)$ (lift per unit span, divided by dynamic pressure and chord; dimensionless) against the **angle of attack** α in degrees, the angle between the chord line and the oncoming flow. Three features matter:

- **The linear region.** From a few degrees below zero to about 10°, $C_l$ rises in a straight line. Its slope is the **lift-curve slope** $a_0 = dC_l/d\alpha$. Thin-airfoil theory, which you will derive in Phase 10, predicts $a_0 = 2\pi$ per radian $= 0.1097$ per degree for *any* thin section; measured values are a little lower, about 0.10 to 0.11 per degree. This is your first check number: a slope you read as 0.05 or 0.2 per degree means you misread an axis.
- **The zero-lift angle** $\alpha_{L0}$, where the line crosses $C_l = 0$. For a symmetric section it is 0° by symmetry. Positive camber moves it negative: a NACA 2412 makes lift at zero angle of attack and has $\alpha_{L0} \approx -2°$. More camber, more negative.
- **Stall.** Above roughly 12° to 16° the flow separates from the upper surface, the curve bends over, reaches a maximum $C_{l,\max}$ (about 1.5 to 1.7 for these sections at these Reynolds numbers), then falls. The stall region depends strongly on Reynolds number and roughness; the linear region barely does.

The plots show several curves per section, labelled by **Reynolds number** $Re = \rho V c / \mu$ (the ratio of inertial to viscous effects; $\mu$ is the air's viscosity). "3, 6, and 9 million" are the values for a full-size wing at flight speeds. A 30 cm model glider wing at 8 m/s has $Re \approx 1.6 \times 10^5$, forty times smaller, and its stall and drag can be quite different (problem 6).

### 2.7 Reading a drag polar

The second panel plots the **drag coefficient** $C_d$ against $C_l$ (not against α). This is the **drag polar**. Read:

- **Minimum drag** $C_{d,\min}$, the lowest point, and the $C_l$ at which it occurs (near zero for a symmetric section, slightly positive for a cambered one). Values of 0.006 to 0.008 are typical for smooth four-digit sections at these Reynolds numbers. The "standard roughness" curve, with grit applied near the leading edge to mimic a dirty wing, roughly doubles it.
- **Maximum lift-to-drag ratio.** $L/D = C_l/C_d$ is the slope of a line from the origin to a point on the polar. The maximum is where a line from the origin is *tangent* to the polar. For a parabolic polar $C_d = C_{d0} + kC_l^2$ this happens at $C_l^* = \sqrt{C_{d0}/k}$ with $(L/D)_{\max} = 1/(2\sqrt{C_{d0}k})$; your test uses this analytic case to check the reader. Section L/D values are much higher than whole-aircraft values, because the wing-tip (induced) drag and the fuselage are not there.

Precision: on a scanned plot you can read $C_l$ to about ±0.02 and $C_d$ to about ±0.0005, one minor division. Record that as the uncertainty of every digitised point (`P00.10`), and propagate it into your fitted slope by refitting with the points perturbed.

### 2.8 Where the data stop being true

Two-dimensional: a real wing has tips, and its lift-curve slope and drag are both changed by the trailing vortices (Phase 10). Low speed: the tunnel Mach number was low, so none of this applies above about Mach 0.6, where compressibility changes everything (Phase 11). Reynolds number 3–9 million: small models and large airliners both fall outside. Smooth surface or one standard roughness: a wing with ice, bugs, or rivets is neither. Steady flow: a flapping or gusting condition is not covered. A reading note that omits these is incomplete, however good its numbers.

## 3. Build the computational model

Open `starter/main.py`. `REQUIRED_NOTE_FIELDS`, the `Naca4` dataclass, `EXAMPLE_NOTES`, and `main` are given.

### Step 1: Parse the designation

Accept `"2412"`, `"NACA 2412"`, and `"naca2412"`: strip an optional `naca` prefix (any case) and whitespace, then demand exactly four digits.

```python
_NACA4_RE = re.compile(r"^\s*(?:naca)?\s*(\d)(\d)(\d{2})\s*$", re.IGNORECASE)

def parse_naca4(designation: str) -> Naca4:
    m = _NACA4_RE.match(designation)
    if not m:
        raise ValueError(f"not a NACA four-digit designation: {designation!r}")
    d1, d2, d34 = int(m.group(1)), int(m.group(2)), int(m.group(3))
    if d34 == 0:
        raise ValueError("thickness digits '00' mean a zero-thickness section")
    if d1 > 0 and d2 == 0:
        raise ValueError("a cambered section needs a non-zero camber position")
    return Naca4(f"{d1}{d2}{d34:02d}", d1 / 100.0, d2 / 10.0, d34 / 100.0)
```

Raising `ValueError` on bad input is the right behaviour: a silent guess (`"2412 "` becoming `241`) is how a wing gets built to the wrong section.

### Step 2: Camber line and slope

Implement equations (1) and (2) with `np.where(x < p, front, back)`. Return zeros when `m == 0` before dividing by `p`, which is zero for a symmetric section. Write `camber_slope` the same way with the derivative expressions from §2.4.

### Step 3: Thickness and surfaces

Equation (3) is one line. For `airfoil_surfaces`, use cosine spacing, `x = 0.5 * (1 - cos(beta))` with `beta` from 0 to π, so points crowd toward the nose where the curvature is high; then lay the thickness off perpendicular to the camber line as in §2.5.

### Step 4: Lift-curve fit

Select the points inside `alpha_range` (default −6° to 8°, the linear region), fit a straight line with `np.polyfit(alpha, cl, 1)`, and return the slope and $\alpha_{L0} = -\text{intercept}/\text{slope}$. Raise `ValueError` with fewer than two points in range: a fit through one point is not a fit.

### Step 5: Drag polar summary

`cd_min` and its `cl` come from `np.argmin(cd)`. For L/D, mask out `cl <= 0` (the ratio is meaningless there), compute `cl/cd`, and take the argmax.

### Step 6: Notes validator

Find the headings (any level). For each required field, find the first heading containing it (case-insensitive). If none, report `missing field: <name>`. Otherwise collect the lines up to the next heading, join and strip them, and if what remains is empty or a placeholder (`...`, `TODO`, `tbd`), report `empty field: <name>`.

### Verify

`pytest curriculum/phase-00/lesson-11-reading-a-technical-paper`. The analytic cases: the camber line peaks at exactly $x = p$ with value $m$ and is continuous with zero slope there; `camber_slope` equals a central-difference derivative of `camber_line`; the thickness maximum is $t/2$ near 30 % chord and the trailing edge is thin but finite; the 0012 surfaces mirror. The limiting cases: `m = 0` gives a flat camber line; a fit with one point raises. The lift-curve fit must recover a planted slope of 0.11/deg and $\alpha_{L0} = -2°$ despite a planted stall; the polar reader must find $C_l^* = \sqrt{C_{d0}/k}$.

### Compare

Run `python starter/main.py`. It saves `naca_airfoils.png` with the 2412, 0012, and 4415 overlaid. Compare it with the airfoil drawings in the report: the 4415 should be visibly fatter and more arched than the 2412, and the 0012 should be the 2412 with the arch removed. Then compare your fitted slope from part 4 with the thin-airfoil value 0.1097/deg; it should be a few percent lower, not higher.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** The dataset is NACA Report 824 itself: measured section characteristics from the Langley two-dimensional low-turbulence pressure tunnel, 1940s, public domain as a US Government work. Download the PDF from the NTRS citation page (the inline preview does not work; the file is a large scan).

1. **Pass 1.** Read the summary, the table of contents, and flip through every page looking only at figure captions. Write the Question field of your notes.
2. **Pass 2.** Read the description of the four-digit series and the description of the tests. Write the Method field. Then find the section characteristics plots for NACA 0012 and NACA 2412 (they are in the plotted-data part of the report, in numerical order of designation). For the 2412 at the middle Reynolds number (about 6 million, smooth), digitise by eye about eight points of $C_l$ against α between −6° and +14° and about eight points of $C_d$ against $C_l$ across the polar. Record each with its reading uncertainty in a CSV: `alpha_deg, cl, u_cl` and `cl, cd, u_cd`.
3. Run `lift_curve_fit` on the first table and `drag_polar_summary` on the second. Then perturb every $C_l$ by $+u$ and refit, and by $-u$ and refit; half the spread in the slope is your reading uncertainty on it.
4. Repeat the lift-curve fit for the 0012 and confirm $\alpha_{L0} \approx 0$.
5. Write the Key figure, Result, Limitations, and My questions fields. Run `validate_reading_notes` on your file (add a few lines to `main` or a short script) until it reports nothing.

Record the page numbers you used. Someone checking your notes should land on the same plot.

## 5. Validate and improve

Put three things side by side: thin-airfoil theory, your reading of the report, and your code's drawing.

- **Slope.** Theory: 0.1097/deg. Your 2412 reading: expect roughly 0.10 to 0.11/deg with an uncertainty of a few thousandths. If your value is outside 0.09–0.115, the most likely causes, in order, are a misread axis scale, points taken into the stall region, and a skewed scan. Check each.
- **Zero-lift angle.** Your 2412 reading should be about −2°; your 0012 about 0°. A 0012 reading of −1° tells you the scan or your eye has a systematic offset; apply the same offset to the 2412 and say so.
- **Geometry.** Measure the maximum thickness and its position on your plotted 2412 (`np.max(yu - yl)` and the `x` where it occurs). You should get 0.120 at $x \approx 0.30$, the designation and equation (3) agreeing because one was built from the other. This confirms the code, not the physics.
- **Drag.** Your $C_{d,\min}$ for the smooth 2412 should be near 0.006–0.007; the roughness curve near double. If the smooth and rough readings are close, you read two Reynolds-number curves instead.

Then improve the one thing that limited you most. If reading precision dominated, digitise twice as many points and refit; the slope uncertainty should fall roughly as $1/\sqrt{n}$. If it was figure location, note the page number in your Key figure field so the next reader skips the search.

## 6. Defend your work

Submit: `reading-notes.md` (validator clean), both CSVs with uncertainties, the fitted slope and $\alpha_{L0}$ for 2412 and 0012 with reading uncertainties, `naca_airfoils.png`, the pytest output, and a paragraph on which limitation in §2.8 matters most for a light aircraft wing and which for a model glider.

### Adversarial self-checks

- What if the camber doubled (NACA 4412)? Which of slope, $\alpha_{L0}$, $C_{l,\max}$, and $C_{d,\min}$ would you expect to change, and in which direction? Check the report.
- Which assumption of the two-dimensional data fails first on a real wing: infinite span, low Mach, or the smooth surface? For which aircraft does the answer change?
- What measurement could falsify your reading? (Someone else digitising the same plot and getting a slope outside your stated uncertainty.)
- Does your fitted slope survive a change in `alpha_range` from (−6, 8) to (−4, 6)? If it moves by more than your uncertainty, your uncertainty was too small.
- Which decision would change if your reading uncertainty were twice as large? Could you still distinguish the 2412's slope from the thin-airfoil value?

## Check yourself

1. Decode NACA 6412: camber, camber position, thickness, as fractions of chord.
2. Why does the second digit of a NACA 00xx section carry no information?
3. The lift-curve slope of thin-airfoil theory is $2\pi$ per radian. Convert it to per degree.
4. A cambered section has $\alpha_{L0} = -2°$. Is its lift at α = 0 positive, negative, or zero?
5. Define the Reynolds number and state what "9 million" tells you about the test conditions.
6. On a drag polar, what geometric construction locates the maximum lift-to-drag ratio?
7. Why is the trailing-edge thickness of a four-digit section not exactly zero, and does that matter for your plot?
8. Your validator passes a notes file whose Result field says "see figure". Is that a bug in the validator or in the notes?
