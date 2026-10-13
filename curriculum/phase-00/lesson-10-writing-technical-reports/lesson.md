# P00.10  Writing Technical Reports

> A result nobody can check is not a result. The report is how you let someone check.

| | |
|---|---|
| **Learning outcomes** | 1. Write a technical report in the nine-section structure of `projects/_template/report-template.md`, with a prediction recorded before the measurement. 2. Number and caption every figure, table, and equation so a reader can be sent to one by name. 3. Attach an uncertainty to every measured or derived number and state how it was obtained. 4. Implement and test a Markdown report linter that finds missing sections, uncaptioned figures, and unqualified numbers in the Results. 5. Decide, for a given report, whether its conclusion is supported by its own evidence. |
| **Prerequisites** | `P00.04` (saving labelled figures), `P00.05` (units), `P00.07` (pytest basics); skill check: can write a Python function that loops over the lines of a string. |
| **Engineering question** | Can a reader who was not in the room reproduce your measurement and reach the same conclusion from your report alone? |
| **Assumptions and limits** | The structure here suits a short experimental or computational report (2–10 pages). It is not a design review package, a journal paper, or a certification document, which add sections (requirements traceability, compliance matrices). The linter checks form, not truth: a report can pass every check and still be wrong. |
| **Analytical work** | Propagating one uncertainty into a derived quantity by hand; deciding significant figures from the uncertainty; checking that a conclusion follows from a prediction–measurement comparison. |
| **Computational deliverable** | `starter/main.py`: a linter that reads a Markdown report and lists missing required headings, figures without a numbered caption, non-sequential figure numbers, and numbers-with-units in the Results section that carry no ±. Tests prove each rule on fixture strings and on a file in `tmp_path`. |
| **Practical activity** | **data-based investigation**: take one measurement you already have (from any earlier lesson or a kitchen scale and a ruler) and write it up as a complete report that your linter passes. |
| **Equipment** | None beyond a computer. Optional: a kitchen scale and a ruler to generate a fresh dataset. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | The report (Markdown), the linter output showing zero problems, the passing test run, and one paragraph on which check you think the linter should also make but does not. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; competencies 3 (computational) and 4 (interpretation) are the focus; competency 1 through the uncertainty propagation in problem 2. Pass: linter tests green, report passes the linter, and the Discussion ranks at least two causes for any disagreement. |
| **References** | `@psu-style-for-students` Chapter 4 "Equations, Figures, and Tables" (required-free); Chapter 5 "Using Sources" (required-free); Chapter 6 "Writing Documents for Classes" (required-free). `@psu-writing-engineer` "Reports" / "Tutorial: Reports" (required-free). |
| **Extensions** | ★ Convert the report to PDF with `pandoc` and a LaTeX engine, and make the linter also check that every `Figure n` mentioned in prose exists. ★ Add a check that each equation is numbered `(n)` in order. |
| **Time** | unmeasured |

## 1. The engineering problem

In 1999 the Mars Climate Orbiter was lost because one team's software produced thruster impulse in pound-force seconds and another team's software read it as newton seconds (`P00.05`). The mishap report lists, among the contributing causes, that the anomaly had been noticed during the cruise and was not written up and escalated in a form that forced a decision. The numbers existed. The report that would have made someone act on them did not.

That is the job of a technical report: to move a result from your head into a form that another engineer can check, reproduce, and act on without you present. "Act on" is the point. A report is read by someone who has to decide something: accept the design, repeat the test, change the model, ground the aircraft. Every section of the structure you will learn exists to make that decision possible and to make it hard to fool yourself along the way.

Success criteria for this lesson. By the end you can take a measurement you have already made, write it up so that (a) every number has units and an uncertainty, (b) the prediction you made *before* measuring is on the page next to the result, (c) every figure and equation can be referred to by number, and (d) a program you wrote confirms the structure is complete. If you cannot do this, your later lab reports in `P00.14`, `P00.15`, and every project in the course will be rejected at review for reasons that have nothing to do with the physics.

## 2. Understand and derive

### 2.1 The nine sections and why each one exists

The repository's template, `projects/_template/report-template.md`, has a header line and nine sections. Each one answers a question a reviewer is going to ask. Learn them as questions, not as headings.

| Section | The reviewer's question | What fails if it is missing |
|---|---|---|
| Summary | What did you find, how sure are you, and should I read on? | A busy reader stops here and takes away nothing. |
| 1. Objective and requirements | What were you trying to predict or decide, and what counted as success *before* you started? | Success gets redefined after the fact to match whatever came out. |
| 2. Theory and assumptions | What model did you use, and where is each assumption used? | Nobody can tell which assumption to blame when the model disagrees with the data. |
| 3. Prediction | What number did you expect, with what uncertainty, computed how? | You cannot be surprised by your data, so you cannot learn from it. |
| 4. Apparatus, calibration, and safety | What did you measure with, how well did it read, and what did you do to stay safe? | The result cannot be reproduced and the uncertainty cannot be believed. |
| 5. Results | What did you observe? (Numbers, uncertainties, plots. No interpretation yet.) | Facts and opinions blur, and a reader cannot form their own view. |
| 6. Discussion | Do prediction and measurement agree beyond uncertainty? If not, what are the ranked causes, and what follow-up test did you run? | "Experimental error" gets written down as if it were an explanation. |
| 7. Conclusions and redesign | What does the evidence support, and what would you change? | The report describes activity instead of reaching a decision. |
| 8. Reproducibility | Where are the code (commit hash), raw data, calibration files, and BOM tier? | The work cannot be audited or extended. |

Two of these rows are the ones beginners skip, and they are the ones that matter most. The **Prediction** section is written and dated before the measurement, and it is never edited afterwards. If the prediction turns out wrong, you write that in the Discussion. A report whose prediction matches its result to four digits because the author adjusted it afterwards is worthless, and a reviewer can usually tell. The **Reproducibility** section is what makes the report a scientific document rather than an essay: a commit hash (from `P00.02`) pins the exact code; a path to the raw data lets someone redo the reduction.

### 2.2 Results and Discussion are different sections for a reason

The Results section states what was observed: "The mean period over 20 swings was 2.010 ± 0.005 s (Figure 1)." The Discussion interprets it: "This is 0.004 s longer than predicted, which is less than one standard deviation, so prediction and measurement agree." Keeping them apart lets a reader disagree with your interpretation while trusting your observations. Mixing them invites a particular failure: you notice the result is 2 % off, you write "which is probably because of air drag", and you never check whether air drag could produce 2 %. (For a pendulum bob it cannot; see problem 3.)

### 2.3 Every number has units and an uncertainty

A measured quantity is not a number. It is a number, a unit, and a statement of how far the true value might lie from it. Write it as

$$ L = 1.000 \pm 0.002\ \text{m} \qquad (1) $$

and say in the Apparatus section where the ± 0.002 m came from (half the smallest division of the tape measure, in this case, which is a Type B evaluation; `P00.14` goes deeper). Three rules keep this honest:

1. **The uncertainty decides the significant figures.** Round the uncertainty to one or two significant figures, then round the value to the same decimal place. `2.0103 ± 0.0051 s` becomes `2.010 ± 0.005 s`. Writing `2.0103 ± 0.005 s` claims a digit you do not have.
2. **Derived quantities carry a propagated uncertainty.** If $g = 4\pi^2 L / T^2$, then for small, independent uncertainties the relative uncertainties add in quadrature:

$$ \frac{u_g}{g} = \sqrt{\left(\frac{u_L}{L}\right)^2 + \left(2\,\frac{u_T}{T}\right)^2} \qquad (2) $$

   The factor 2 on the period term is the exponent of $T$. Equation (2) is a special case of the general law of propagation of uncertainty; it holds when the uncertainties are small compared with the values and uncorrelated. For the pendulum numbers above: $u_L/L = 0.002$, $2u_T/T = 2 \times 0.005/2.010 = 0.0050$, so $u_g/g = \sqrt{0.002^2 + 0.0050^2} = 0.0054$ and $g = 9.77 \pm 0.05\ \text{m/s}^2$. Notice which term dominated. That tells you what to improve next time: the timing, not the ruler.
3. **An uncertainty is not a guess at how wrong you are.** It is a width you can defend from the instrument, the scatter of repeated readings, or a calculation. "± 10 %" with no source is not an uncertainty; it is a shrug.

Two notations are acceptable in this course: `9.81 ± 0.05 m/s²` and the concise `9.81(5) m/s²`, where the digits in parentheses are the uncertainty in the last digits shown. Your linter accepts both.

### 2.4 Numbering figures, tables, and equations

A reader must be able to say "look at Figure 2" or "from equation (3)" and have you both land on the same thing. The rules are mechanical:

- Figures and tables are numbered separately, each from 1, in order of first appearance. Every one has a caption that begins with its label, "Figure 2." or "Table 1.", followed by a sentence that lets the figure stand on its own: what is plotted against what, with units, and what the error bars mean. The caption goes *below* a figure and *above* a table (the usual convention; the point is to pick one and be consistent).
- Every figure and table is mentioned in the prose by number before or near where it appears. A figure nobody refers to is decoration.
- Display equations that you refer to later are numbered `(1)`, `(2)`, … at the right margin, as in equations (1) and (2) above. Refer to them as "equation (2)" or "Eq. (2)". Equations you never refer to need no number.
- Axes carry quantity and unit, `Period T (s)`, and the plotted points carry error bars sized to the stated uncertainty (`P00.04`).

In Markdown a figure is an image line followed by its caption as the next paragraph:

```markdown
![Measured period against string length](figures/period_vs_length.png)
*Figure 1. Measured pendulum period against string length. Error bars are one standard deviation of 20 swings; the line is equation (1) with g = 9.81 m/s².*
```

The alt text inside the square brackets is for screen readers and for when the image fails to load. It is not the caption, and your linter will treat it as such.

### 2.5 Writing that a tired reviewer can follow

Short sentences. One idea per paragraph. The subject of a sentence is the thing that acts: "The balance read 0.3 N high at full scale", not "It was observed that a high reading was present." Past tense for what you did, present tense for what is true ("the model predicts"). Define a symbol the first time you use it and use the same symbol throughout. Put derivations longer than a few lines in an appendix and refer to them. Never write "clearly" or "obviously"; if it were, you would not need the word.

### 2.6 Where this structure stops being adequate

The nine sections suit an experiment or a computational study with one or a few questions. A design report adds requirements, trade studies, margins, and a verification matrix (Phase 23). A failure investigation adds a timeline and a fault tree. A paper for a journal compresses Apparatus and Reproducibility into "Method" and "Data availability" and adds a literature review. The habit that transfers to all of them is the same: a reviewer's questions, answered in a fixed order, with every number qualified.

## 3. Build the computational model

The "model" in this lesson is a program that reads a Markdown report and checks its structure. The lesson in doing this is not the checker itself; it is that the rules in part 2 are precise enough to be written as code. Open `starter/main.py`. The constants `REQUIRED_HEADINGS`, `_NUM`, `_UNIT`, and `_PM` are given.

### Step 1: Find the headings

A Markdown heading is a line that starts with one to six `#` characters, a space, then text. Headings inside fenced code blocks (between triple backticks) are code, not headings, so track whether you are inside a fence.

```python
_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")

def find_headings(text: str) -> list[tuple[int, int, str]]:
    out = []
    in_fence = False
    for i, raw in enumerate(text.splitlines(), start=1):
        line = raw.rstrip()
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            continue
        m = _HEADING_RE.match(line)
        if m:
            out.append((i, len(m.group(1)), m.group(2)))
    return out
```

Line numbers start at 1 because that is what editors show. The regular expression's groups are the run of `#` (whose length is the heading level) and the heading text.

### Step 2: Required headings by keyword

A real report writes "## 3. Prediction (recorded before measurement)", not "prediction". So match each required keyword as a case-insensitive substring of any heading. Join all heading texts into one lowercase string and test each keyword against it.

```python
def missing_headings(text: str) -> list[str]:
    present = " || ".join(h.lower() for _, _, h in find_headings(text))
    return [k for k in REQUIRED_HEADINGS if k not in present]
```

### Step 3: The body of one section

The uncertainty rule applies only to the Results section (the Theory section legitimately says "a 1.000 m string"). So you need the lines that belong to a section: from the line after its heading up to, but not including, the next heading of the same or a higher level (fewer `#`). Return `(line_number, text)` pairs so that later messages can say which line.

### Step 4: Figures and captions

An image is `![alt](path)`. Its caption is the first non-blank line after it, if that line starts with `Figure <n>` after any `*` or `_` emphasis markers. If the next non-blank line is anything else, the figure has no caption. Build a `Figure` dataclass for each image; `Figure.number` (already written for you) parses the number from the caption.

Then `figure_numbering_issues` checks that the captioned figures, in document order, are numbered 1, 2, 3, …: the i-th captioned figure must be "Figure i".

### Step 5: Numbers with units but no ±

This is the heart of the linter and the hardest regular expression you have written so far. Build it from the three given pieces. There are four shapes to recognise, and the order you list them in matters because the regex engine takes the first alternative that matches at a position:

| Shape | Example | Verdict |
|---|---|---|
| `NUM ± NUM UNIT` | `9.81 ± 0.05 m/s²` | accepted |
| `NUM UNIT ± NUM [UNIT]` | `2.50 kg ± 0.01 kg` | accepted |
| `NUM(digits) UNIT` | `21.5(5) °C` | accepted |
| `NUM UNIT` | `9.81 m/s²` | **flagged** |

Two guards stop false alarms. A number must not be glued to a preceding word character or dot, so `x2 kg` and `v1.2 m` do not start a match. A unit must not be followed by a word character, so the `1 s` inside "Table 1 shows" does not match (`s` is followed by `h`). Name the alternatives with `(?P<a>...)`, `(?P<b>...)`, `(?P<d>...)`, `(?P<c>...)` and return only the matches whose group `c` is non-empty.

```python
_END = r"(?![\w/²^°])"
_START = r"(?<![\w.])"
_QUANTITY_RE = re.compile(
    rf"{_START}(?:"
    rf"(?P<a>{_NUM}\s*{_PM}\s*{_NUM}\s*{_UNIT}{_END})"
    rf"|(?P<b>{_NUM}\s*{_UNIT}\s*{_PM}\s*{_NUM}(?:\s*{_UNIT})?{_END})"
    rf"|(?P<d>{_NUM}\(\d+\)\s*{_UNIT}{_END})"
    rf"|(?P<c>{_NUM}\s*{_UNIT}{_END})"
    rf")"
)
```

Read it once slowly. `rf"..."` is a raw f-string, so backslashes stay literal and `{_NUM}` is substituted. `(?<![\w.])` is a negative look-behind: "not preceded by a word character or a dot". `(?![\w/²^°])` is a negative look-ahead. Unit lists are ordered longest-first (`m/s²` before `m/s` before `m`) so the engine does not stop at `m` and leave `/s²` dangling.

### Step 6: Assemble the report and the file entry point

`lint_report` collects messages in the order a reader fixes them: missing headings, uncaptioned figures, numbering, then bare quantities in Results with their line numbers. `lint_file` reads a path and calls it. `main` is already written: it returns 0 on a clean report and 1 otherwise, so it can run as a pre-commit hook later.

### Verify

Run `pytest curriculum/phase-00/lesson-10-writing-technical-reports`. The tests are the "analytic cases" of this lesson: a known-good fixture report must produce zero problems; removing one heading must produce exactly that keyword; a figure with no caption, a swapped numbering, and a single bare `2.010 s` in Results must each be reported, with the line number; and the repository's own template must satisfy the heading rule. The limiting cases are an empty section (returns an empty list, not an error) and a heading inside a code fence (ignored).

### Compare

Run the finished linter on `projects/_template/report-template.md` (`python starter/main.py` with no arguments does this). It passes, which is a sanity check on the linter, not on the template: an empty skeleton has no figures and no numbers to get wrong. Then run it on any Markdown file that is *not* a report, such as this `lesson.md`. Read the list of "missing" headings it produces. That list is the difference between an explanation and a report.

## 4. Experiment at home, or investigate real data

**Mode: data-based investigation.** You will write one complete report and make it pass your linter.

Choose a measurement you already have. Options, in order of preference:

1. Data you took in `P00.08` (a flight log you parsed) or any plot from `P00.04` or `P00.06`.
2. A fresh ten-minute measurement: weigh five identical coins one at a time on a kitchen scale, then all five together. Predict the total from the singles, with uncertainty, *before* weighing the group.
3. Any dataset from earlier in your own life (a bicycle computer log, weather records) for which you can write down a prediction before analysing it.

Then write `report.md` from the template. Fill every section, however briefly. The Prediction section gets a date. The Results section contains at least one figure with a numbered caption and at least two numbers with units; every one carries a ± or the concise notation. The Discussion ranks causes if prediction and result disagree by more than the combined uncertainty, and says so explicitly if they agree. The Reproducibility section names the file with your raw numbers and the commit hash of the code that reduced them.

Run `python starter/main.py report.md`. Fix everything it lists. Then read the report aloud once; sentences you stumble on are sentences to shorten.

## 5. Validate and improve

Put the linter's view and a human's view side by side. The linter says the report is complete. Now check what it cannot:

- Does the Summary state the headline number with its uncertainty, and does that number appear identically in Results? A report that says `9.79 ± 0.06 m/s²` on page 1 and `9.8 m/s²` on page 3 has an inconsistency the linter missed because the summary is not the Results section.
- Does every figure that exists get mentioned in the prose? The linter checks captions, not references.
- Does the Discussion's ranked cause list include at least one cause you then *tested*? "The follow-up test and what it showed" is in the template for a reason.
- Does the stated uncertainty on a derived number follow from the stated uncertainties on the inputs by equation (2) or its equivalent? Recompute one by hand.

Each gap you find is a candidate rule. Pick the one you think is most valuable, write it as a test first (`P00.07`), then implement it. Record which you chose and why in your evidence paragraph. If you implemented "every figure is referenced in prose", rerun the linter on your report: did it pass?

A known limitation of the uncertainty check: it keys on a short list of units. A result written in knots or feet (which industry practice sometimes demands for altitude and airspeed) will slip past unflagged, and a quantity written in words ("about two seconds") is invisible to it. Say this in your report's limitations if you rely on the linter.

## 6. Defend your work

Submit: `report.md`, the figure files it references, the raw data file, the terminal output of `python starter/main.py report.md` showing no problems, the pytest output, and your evidence paragraph on the missing check.

### Adversarial self-checks

- What if a report had a hundred figures? Would `figure_numbering_issues` still run in a blink, and would a human still want one message per figure, or a single "figures renumbered from 37 onward"?
- Which rule in the linter fails first on a real report from a different field, say a chemistry lab? (Hint: the unit list.)
- What measurement, or rather what report, could falsify the claim that "a complete structure means a trustworthy report"? Write one that passes every check and is still wrong.
- Does the uncertainty check survive a change of notation, say `9.81 m/s² (±0.05)`? What should it do with that?
- Which decision in your Discussion would change if every uncertainty in your Results were doubled? If the answer is "none", your experiment was not sensitive enough to test the prediction; say so.

## Check yourself

1. In which section does the sentence "The 2 % discrepancy is most likely caused by the stopwatch reaction time" belong, and why not in Results?
2. A tape measure has 1 mm divisions. You read 1.000 m. State the length with an uncertainty and name the type of evaluation you used.
3. You compute $g = 4\pi^2 L/T^2$ with $u_L/L = 0.2\,\%$ and $u_T/T = 0.25\,\%$. What is $u_g/g$, and which measurement should you improve first?
4. Write `2.0103 ± 0.0051 s` with the correct number of significant figures.
5. Why must the Prediction section be dated and never edited after the measurement?
6. Your linter flags `300 K` in the Results but not `300 K` in the Theory section. Is that a bug?
7. What does the alt text of a Markdown image do, and why is it not a substitute for a caption?
8. A report has figures numbered 1, 2, 4. Which function reports this, and what exactly does it say?
