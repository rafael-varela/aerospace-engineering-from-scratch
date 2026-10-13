# P00.04  Plotting That Engineers Trust

> A plot is a claim. Label it so a stranger can check the claim without asking you.

| | |
|---|---|
| **Learning outcomes** | 1. Build a figure with Matplotlib's explicit Figure/Axes API and label every axis with a quantity and a unit. 2. Choose and justify linear or logarithmic axes from the shape of the data. 3. Draw measured data with error bars and models as lines, in a fixed colour scheme. 4. Implement and test a small plotting module whose tests inspect the figure's objects, not its pixels. 5. Decide whether a given figure is fit to go into a technical report, and say why not. |
| **Prerequisites** | `P00.01`, `P00.03`; skill check: you can write a Python function with a docstring and run `pytest`; you can compute the ISA pressure array from P00.03. |
| **Engineering question** | Can a reviewer who has never spoken to you read the pressure at 15 km, the drag at CL = 1.0, and the uncertainty on a thrust sample off your figures alone, and get the same numbers you did? |
| **Assumptions and limits** | The house style targets static figures in printed or PDF reports; it is not a style for interactive dashboards. Colours are chosen to be distinguishable in print and for common colour-vision deficiencies, but the style check cannot see whether a line is hidden behind another. The ISA example model is valid 0–20 km only. |
| **Analytical work** | Why pressure needs a log axis (the exponential decay derivation); why a drag polar does not; the meaning of an error bar as ±1 standard deviation and how it propagates to the eye. |
| **Computational deliverable** | A `main.py` module with `axis_label`, `new_figure`, `plot_series`, `set_log_axes`, `save_vector`, `check_figure` and three example figures. Tests prove: labels carry units, log scale is set where required, error bars exist, colours come from the house palette, saved files are vector, and the style checker flags a bad figure. |
| **Practical activity** | **computational experiment**: re-plot the same ISA data four ways (linear/log, with/without units) and record what a reader can and cannot recover from each. |
| **Equipment** | None. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Your `main.py`, passing tests, the three SVG figures, and a one-page critique of a published plot of your choice. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 3 (computational competence) and 4 (engineering interpretation); pass: all public tests green and a critique that names at least three defects by their effect on the reader. |
| **References** | `@matplotlib-userguide` "Quick start guide" (required-free); "Axes and subplots" (required-free); "Text" (required-free); "Colors" (required-free); "Customizing Matplotlib with style sheets and rcParams" (required-free). |
| **Extensions** | ★ Write the house style as a `.mplstyle` file and load it with `plt.style.use`. ★ Add a `check_figure` rule that detects two series drawn in the same colour. |
| **Time** | unmeasured |

## 1. The engineering problem

A flight-test team sends you a plot of a motor's thrust against time. The vertical axis reads "Thrust". Is that newtons or pounds-force? The points have no error bars, so you cannot tell whether the wobble near the peak is real or noise. The curve is drawn in a default blue, the same blue as the model curve in the previous figure, and the two figures use different axis ranges, so you cannot compare them by eye. You write back to ask. The answer takes two days. The review slips.

This is the normal failure mode of engineering plots. A plot is a compressed report: it makes a claim about how one quantity depends on another, with some uncertainty. If the reader has to ask what the units are, the plot has failed at its only job. If the uncertainty is not shown, the reader cannot judge whether a difference between two curves means anything. If the colours change meaning from one figure to the next, the reader has to re-learn each figure.

The fix is not taste; it is a *house style*, a small set of rules applied by code so that you cannot forget them. In this lesson you build that code. The success criteria are concrete: a stranger can read the pressure at 15 km from your ISA figure to within the thickness of the line, the unit of every axis is printed on the axis, measured data and model curves are visibly different kinds of thing, every figure is saved as a vector file that stays sharp when a report is zoomed, and an automatic checker refuses a figure that breaks any of those rules.

## 2. Understand and derive

### 2.1 Anatomy of a figure

Matplotlib has two layers. The convenience layer (`plt.plot`, `plt.xlabel`) acts on "the current figure", which is whichever one you touched last. The object layer gives you a `Figure` (the page) that holds one or more `Axes` (a plotting area with an x axis and a y axis), and each `Axes` holds the drawn `Line2D`, text, and legend objects. Engineers use the object layer, because a function that receives an `Axes` can draw into it without caring where it lives, and a test can ask the `Axes` what it contains.

```
Figure ──────────────────────────────────────────────┐
│  Axes                                               │
│   ┌ title ─────────────────────────────────────┐    │
│   │ y label [unit]  ┆ Line2D (model, red)      │    │
│   │                 ┆ ErrorbarContainer (grey) │    │
│   │                 ┆ Legend                   │    │
│   └──────────────── x label [unit] ────────────┘    │
└─────────────────────────────────────────────────────┘
```

Everything a reader needs is on this diagram: the quantity and unit on each axis, a legend naming each series, and a visual distinction between a measured series and a modelled one.

### 2.2 When does an axis need a logarithm?

Take the ISA atmosphere from P00.03. In the isothermal layer above 11 km the pressure obeys the hydrostatic balance `dp/dh = -ρ g` with `ρ = p / (R T)` and constant `T`, so

```
dp/dh = -(g / (R T)) p      ⇒      p(h) = p11 · exp(-(h − h11) g / (R T)).
```

Pressure falls by the same *factor* for every kilometre climbed: roughly ×0.855 per km at 216.65 K. On a linear axis that curve is steep at the bottom and flat at the top, and the stratospheric values crowd into a strip a few pixels tall. Taking logarithms,

```
ln p = ln p11 − (g / (R T)) (h − h11),
```

a straight line. Where a quantity changes by factors rather than by increments, a logarithmic axis turns the curve into something the eye can read and compare, and makes a 5 % departure visible at 20 km as clearly as at sea level. The troposphere is not isothermal, so `ln p` is only approximately straight there; the curvature you see on the log plot is real information about the lapse rate, not a plotting artefact.

Rule of thumb: if the ratio of the largest to the smallest value you want the reader to see exceeds about 30, you need a log axis. The drag polar in part 3 has `CD` between 0.02 and 0.15: a ratio of about 7. A log axis would hide the parabola's shape and gain nothing. Pressure from sea level to 20 km spans a factor of 18.5; from sea level to 50 km it would span 1300. The ISA plot gets a log axis for pressure and a linear axis for altitude, because altitude is what you *set*, and a reader wants equal kilometres to look equal.

### 2.3 What an error bar says

An error bar drawn from `y − σ` to `y + σ` says: "if I repeated this measurement many times, about two thirds of the values would land inside this bar." It is a statement about the *measurement*, not about the model. When a model curve passes outside the bars of several neighbouring points, the model is in trouble; when it passes inside all of them, the data cannot distinguish that model from any other curve that also fits inside. The bars therefore set how much the plot can claim. A plot of measured data without error bars is making a claim of infinite precision, which is never true.

In this lesson the thrust-curve "data" are synthetic: generated in code from a known model plus Gaussian noise with a known standard deviation, so the bars are exactly right by construction. Real error bars come from P00.14, where you learn how to estimate `σ` from repeats and calibration.

### 2.4 Colour as meaning

Colour is cheap to change and expensive to read. The house style ties colour to *role*, not to the order in which you happened to plot things: measured data is always one blue, model curves always one red, published reference values always one grey, fitted curves always one teal, limits and envelopes always one brown. A reader who has seen one of your figures can read the rest without the legend. The checker in part 3 refuses any line drawn in a colour outside this palette.

### 2.5 The six defects that cost reviewers the most time

| Defect | What the reader cannot do | House rule that prevents it |
|---|---|---|
| No unit on an axis | Convert, compare with another source, or trust the number at all | `axis_label` always writes `[unit]`; `check_figure` refuses labels without one |
| Linear axis on a quantity spanning decades | Read any value in the small-number end | Log axis whenever max/min > 30 (part 2.2) |
| Measured points joined by lines | Tell measured from interpolated | Error-bar series are markers only (`fmt="o"`) |
| No error bars on measurements | Judge whether model and data disagree | `plot_series(..., yerr=σ)` for every measured series |
| Colours that change meaning between figures | Read the second figure without re-learning it | Role-based palette `HOUSE_COLORS` |
| Raster output | Zoom, select text, print cleanly | `save_vector` accepts `.svg` and `.pdf` only |

None of these rules is about beauty. Each removes one question a reviewer would otherwise have to ask you.

### 2.6 Worked example: what the reader can recover

Suppose the figure shows ISA pressure on a linear axis from 0 to 20 km with ticks every 20 kPa. At 15 km the true value is 12.04 kPa. The curve there sits at about a tenth of the axis height; a reader can estimate "around 10–15 kPa", a ±20 % reading. On a log axis with decades at 10 and 100 kPa, the same point lies 0.08 decades above the 10 kPa line, and a careful reader gets 12 kPa within ±5 %. The plot did not change the data; it changed how much of the data the reader can get back. That is the whole purpose of a plotting style.

## 3. Build the computational model

Open `starter/main.py`. Every function has a docstring that says what the tests expect. Work through them in order; `python main.py` runs end to end once all of them are implemented.

### Step 1: force the headless backend

Tests and scripts must never open a window. The first lines of the module select the `Agg` backend *before* `pyplot` is imported; importing `pyplot` first would lock in whatever backend the machine defaults to.

```python
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
```

### Step 2: labels with units

```python
def axis_label(quantity: str, unit: str) -> str:
    quantity = quantity.strip()
    unit = unit.strip() or "-"
    if not quantity:
        raise ValueError("quantity name must not be empty")
    return f"{quantity} [{unit}]"
```

`"Pressure [Pa]"`, `"Lift coefficient [-]"`. The dash for dimensionless quantities is deliberate: a reader seeing a label without brackets cannot tell "no unit" from "forgot the unit".

### Step 3: a figure factory

```python
def new_figure(xquantity, xunit, yquantity, yunit, *, title=None):
    apply_house_style()
    fig, ax = plt.subplots()
    ax.set_xlabel(axis_label(xquantity, xunit))
    ax.set_ylabel(axis_label(yquantity, yunit))
    if title:
        ax.set_title(title)
    return fig, ax
```

`apply_house_style` updates `plt.rcParams`: figure size, font size, a light grid, no top or right spines, `savefig.format = "svg"`. Everything drawn afterwards inherits those settings.

### Step 4: series by role

```python
def plot_series(ax, x, y, *, role, label, yerr=None, xerr=None):
    color = HOUSE_COLORS[role]            # KeyError for an unknown role is the right behaviour
    x, y = np.asarray(x, float), np.asarray(y, float)
    if x.shape != y.shape:
        raise ValueError(...)
    if yerr is not None or xerr is not None:
        return ax.errorbar(x, y, yerr=yerr, xerr=xerr, fmt="o", color=color, ecolor=color,
                           capsize=2.5, label=label)
    (line,) = ax.plot(x, y, color=color, label=label)
    return line
```

Two decisions are built in. Series with uncertainty are drawn as markers with bars, never joined by a line, because a line between measured points invents data that was not measured. Series without uncertainty are models and are drawn as lines, because a model is defined everywhere. `ax.errorbar` returns an `ErrorbarContainer`, which lands in `ax.containers`; `ax.plot` returns a list of `Line2D`, which land in `ax.get_lines()`. The tests use both lists to count what you drew.

### Step 5: scales, legend, vector output

```python
def set_log_axes(ax, *, x=False, y=False):
    if x: ax.set_xscale("log")
    if y: ax.set_yscale("log")

def save_vector(fig, path):
    path = Path(path)
    if path.suffix.lower() not in VECTOR_SUFFIXES:
        raise ValueError(...)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, format=path.suffix.lstrip(".").lower())
    return path
```

SVG and PDF store the lines as geometry, so a reviewer zooming into your report sees sharp lines and selectable text. PNG stores pixels; it blurs when scaled and triples the report's size for no gain. The house style refuses PNG; if a journal demands one, convert the vector file at the last step.

### Step 6: the checker

`check_figure(fig)` walks `fig.get_axes()` and returns a list of plain-language problems. It returns an empty list for a compliant figure. `main()` runs it on every figure before saving, so a style violation shows up in the terminal, not in a review two weeks later.

```python
def check_figure(fig):
    problems = []
    allowed = {matplotlib.colors.to_hex(c) for c in HOUSE_COLORS.values()}
    for i, ax in enumerate(fig.get_axes()):
        for name, text in (("x", ax.get_xlabel()), ("y", ax.get_ylabel())):
            if not text:
                problems.append(f"axes {i}: missing {name} label")
            elif not (text.rstrip().endswith("]") and "[" in text):
                problems.append(f"axes {i}: {name} label {text!r} has no [unit] tag")
        handles, labels = ax.get_legend_handles_labels()
        if len(labels) > 1 and ax.get_legend() is None:
            problems.append(f"axes {i}: {len(labels)} labelled series but no legend")
        for line in ax.get_lines():
            if matplotlib.colors.to_hex(line.get_color()) not in allowed:
                problems.append(f"axes {i}: line color {line.get_color()!r} is not a house color")
    return problems
```

Three details matter. Colours are compared after `to_hex`, because Matplotlib stores whatever you passed (`"#1f5fa8"`, `"tab:blue"`, an RGB tuple) and two spellings of one colour must compare equal. `get_legend_handles_labels()` returns only series that were given a `label=`; a single labelled series needs no legend. And the checker reads the same objects the tests read, so when the tests pass, the checker is known to be looking at the right things.

### Step 7: the example data and figures

`isa_pressure` is the two-layer ISA from P00.03 (lapse-rate troposphere to 11 km, isothermal above). `drag_polar` is `CD = CD0 + k CL²`. `synthetic_thrust_curve` samples `12 t e^{−t/0.35}` newtons at 25 times with noise of standard deviation `0.08 F + 0.15` N, seeded so the figure is reproducible. The three figure builders put those through the house style: pressure on a log axis with tabulated reference points and 0.5 % error bars; two drag polars (clean and gear-down) on linear axes; and the thrust samples with error bars under the model that generated them.

### Verify

Run `pytest tests/`. The tests are the specification:

- *Analytic case:* `axis_label("Pressure", "Pa") == "Pressure [Pa]"`; `isa_pressure` at the tropopause equals 22 632 Pa within 0.1 % and is continuous across it.
- *Limiting case:* `drag_polar` at `CL = 0` returns `CD0` and is symmetric in `CL`; an empty unit string becomes `[-]`.
- *Conservation:* the data a `Line2D` holds (`line.get_data()`) are exactly the arrays you passed; nothing is resampled or sorted.
- *Object model:* `ax.get_yscale() == "log"` on the ISA figure; exactly one `ErrorbarContainer` after one error-bar series; two `Line2D` after two line series; every line colour maps to a palette colour under `matplotlib.colors.to_hex`.
- *Negative case:* a deliberately bad figure (no y label, no unit, no legend, magenta lines) produces messages mentioning "unit", "y label", "legend" and "color".

### Compare

Build the ISA figure once with the convenience API, `plt.plot(h, p); plt.yscale("log"); plt.xlabel("Altitude [km]")`, and once with `new_figure` and `plot_series`. The pictures are identical. The difference is where the information lives: the first version cannot be checked by `check_figure` without finding "the current figure" through `plt.gcf()`, and a second plot anywhere in the program silently changes which figure that is. Matplotlib's own "Quick start guide" recommends the object-oriented interface for exactly this reason.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Make four versions of the ISA pressure figure: (a) linear axis, no units; (b) linear axis, units; (c) log axis, no units; (d) log axis, units, with the tabulated reference points and error bars. Save all four as SVG. Then hand each one, without comment, to another person (or to yourself after a day) with three questions: what is the pressure at 15 km, what is the unit, and is the model consistent with the reference points? Record the answers and how long each took.

Raw-data format: a small table with columns `figure, question, answer_given, true_value, seconds`. You will reuse this format in P00.10 when you write up the result.

Write down your prediction *before* collecting answers: which figure will give the most accurate reading at 15 km, and by what factor will its error beat version (a)? Part 2.5 gives you a basis for a number.

## 5. Validate and improve

Put the predicted and recorded reading errors side by side. Expect the linear figure's 15 km reading to be off by 20–50 % and the log figure's by under 10 %. If the log figure did *not* win, look at the tick labels: a log axis with only major ticks at 10⁴ and 10⁵ Pa is nearly as hard to read as a linear one, and you should add minor ticks or label them. If the reader gave the wrong unit on versions (b) or (d), the label font is too small or the brackets are missing; fix `apply_house_style` and re-run.

Quantify the uncertainty in your experiment: with one reader per figure, your reading-error estimate has no spread at all. Three readers give a crude standard deviation. Say so in the write-up rather than reporting a single number as if it were exact.

Then improve the module. The most common upgrade after this experiment is a `finish_axes` that also adds minor ticks on log axes and sets the x limits to the data range so the model line does not run past the last reference point.

## 6. Defend your work

Submit: `main.py`, the test output, the four experiment SVGs with their results table, and a one-page critique of a published figure (any free NASA technical report has candidates) naming each defect and the misreading it invites. Use `projects/_template/report-template.md`.

### Adversarial self-checks

- What if the altitude range doubles to 40 km? The pressure ratio becomes about 1300. Does your log axis still have readable ticks? Does the ISA model still apply? (It does not: the model stops at 20 km.)
- What if the thrust noise doubles? The error bars double. Does the model still pass through most of them? What claim could the plot make then?
- Which house rule fails first in the real world? Colour: a greyscale printer makes blue and red the same shade. Check your palette in greyscale and decide whether line style must also carry the role.
- What measurement could falsify the ISA model line on the plot? A reference point whose error bar does not touch the line. Do you have one? Should you?
- Does the result survive a change of backend or Matplotlib version? The tests inspect objects, so yes, unless a default colour name changes; that is why the palette is written in hex.
- Which decision changes if the error bars were twice as large? Whether the model and data "agree" at the peak of the thrust curve. Write down which points would newly agree.

## Check yourself

1. Why must `matplotlib.use("Agg")` appear before `import matplotlib.pyplot`?
2. A label reads "Velocity". Name two distinct ways a reader can misinterpret the plotted numbers.
3. Pressure from 0 to 20 km spans a factor of 18.5. Give the ratio that would be visible as one tick spacing on a log axis with decade ticks, and state whether a linear axis could show the 20 km value to ±5 %.
4. In the house style, measured data are markers with bars and models are lines. Give the physical reason for not joining measured points with a line.
5. `ax.get_lines()` reports two lines but `check_figure` says there is no legend. What did the author forget, and which Matplotlib call fixes it?
6. Why does the style refuse `.png` output for report figures, and when is converting to PNG still legitimate?
7. A drag polar is plotted with `CD` on a log axis. Name one thing the reader loses and one thing they might gain.
8. Two series are drawn in the same house colour because both are "models". Does `check_figure` catch it? What rule would?
