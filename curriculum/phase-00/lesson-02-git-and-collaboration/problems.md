# P00.02  Problem set

Attempt every problem before opening `solutions/`. Record the exact commands you ran and what Git printed.

## 1. (Warm-up) Read the status table

Given this `git status --porcelain` output, state for each file where it sits in the three-tree diagram and whether its change (or which of its changes) will be in the next commit if you run `git commit` right now:

```text
A  curriculum/phase-00/lesson-03-numpy-engineering-calculator/tests/test_extra.py
 M curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md
MM solutions/phase-00/lesson-03-numpy-engineering-calculator/solution.py
?? scratch_plot.png
UU projects/p00-first-experiment/bom.yaml
```

Then write the minimal sequence of commands that commits *only* the lesson text change. Answer form: a five-row table and a command list.

## 2. (Argument) Why small commits

A teammate commits once a week with the message "weekly update". Write a short argument (under 200 words) for small, single-purpose commits with one-sentence messages. Use the drag-coefficient scenario from part 1 of the lesson: show concretely what `git log` and `git diff` can answer under each habit when a number changes unexpectedly. State the strongest case for the teammate's habit and when it would be acceptable.

## 3. (Unfamiliar scenario) A conflict in a bill of materials

You merge a colleague's branch and Git stops with this in `bom.yaml`:

```yaml
- item: load cell, 5 kg, strain gauge type
<<<<<<< HEAD
  quantity: 2
  search_terms: "5kg load cell HX711"
=======
  quantity: 1
  search_terms: "5 kg load cell with HX711 amplifier board"
  note: one spare is not needed; the amplifier is the failure point
>>>>>>> feature/bom-review
```

Decide what the resolved block should contain and justify each line. Then state the commands that finish the merge, and what `git status --porcelain` shows for this file before and after you `git add` it. Finally: your colleague's commit timestamp is `2026-10-09T23:58:00-07:00` and yours is `2026-10-10T08:30:00+01:00`. Which change was made first, in UTC?

## 4. (Numerical, uses your code) Your own history

Produce the fixed-format log of this repository (or any public repository with more than 50 commits that you have cloned) and feed it to your `main.py` with `-`. Report: number of commits, the author table, how many subjects each lint rule flags, and the value of `is_newest_first`. Compare the author table with `git shortlog -sn` and explain every difference using the table in part 5 of the lesson. If `is_newest_first` is False, find the first adjacent pair that is out of order and explain how it came to be (hint: `git log --graph --oneline` around those commits).

## 5. (Numerical, extends your code) Committer versus author

Extend a copy of your parser to read the format `%H|%an|%aI|%cn|%cI|%s` (author name and date, committer name and date, subject). On the same log, count the commits where author and committer differ, and the commits where the committer date is more than one day after the author date. Explain what each count tells you about how the repository's maintainers merge contributions. Answer form: two counts, a two-sentence interpretation, and the lines you changed.

## 6. (Judgment) What goes in the repository

For each of the following, decide whether it belongs in this course's Git repository, in a separate archive referenced by hash, or nowhere, and give the one-sentence reason: (a) the 40 kB CSV of raw timings from `P00.15`; (b) the 300 MB phone video the timings were read from; (c) the PNG plot of the results; (d) the Python script that made the plot; (e) a `.env` file containing a weather-API key used by a later lesson; (f) `__pycache__/`; (g) the virtual environment folder `.venv/`. Then write the `.gitignore` lines that enforce your answers for (c), (f) and (g).
