# P00.02  Worked solutions

Each solution states: assumptions, the reasoning, the result, at least one independent check, and the common mistakes.

## 1. Read the status table

**Result.**

| line | staging area (col 1) | working tree (col 2) | in the next commit? |
|---|---|---|---|
| `A  …/tests/test_extra.py` | new file added | unchanged since | yes, the whole new file |
| ` M …/lesson.md` | nothing staged | modified | no |
| `MM …/solution.py` | modified (staged) | modified again | yes, but only the *first* edit; the later edit is not |
| `?? scratch_plot.png` | untracked | untracked | no, and it should be ignored, not added |
| `UU …/bom.yaml` | unmerged | unmerged | the commit is refused until this is resolved |

**Commands to commit only the lesson text.** First the conflict must be dealt with: either resolve `bom.yaml` and `git add` it (which puts it in the commit too), or `git merge --abort` to leave the merge. Assuming the merge is aborted (or was not in progress in your case): `git add curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md` then `git commit -m "Clarify the geopotential paragraph in P00.03"`. Because `test_extra.py` and the first `solution.py` edit are already staged, committing *only* the lesson text also requires unstaging them first: `git restore --staged <file>` for each. The minimal honest sequence is therefore: abort or resolve the merge; `git restore --staged` the two staged files; `git add lesson.md`; `git commit`.
**Check.** Run `git status --porcelain` after staging: only `M  …/lesson.md` should show in column 1. Run `parse_status_porcelain` on the original text: `staged` has two paths, `unstaged` two, `untracked` one, `conflicted` one; `solution.py` appears in both `staged` and `unstaged`.
**Common mistakes.** Reading ` M` as staged; thinking `git commit` with an unmerged file will succeed; forgetting that already-staged files go along for the ride.

## 2. Why small commits

**Argument (one acceptable version).** When the 10 km density changed from 0.4127 to 0.4135, the question is "which change did that?" Under weekly "weekly update" commits, `git log` offers 25 identical messages and `git diff` between two of them shows every line touched in a week across every file: finding the gas constant among them is a manual search. Under small, single-purpose commits, `git log --oneline -- main.py` shows "Use R = 287.05287 J/(kg K) in the ISA constants" and `git show` of that one commit is a two-line diff. `git bisect` can locate the change automatically among 300 small commits in nine steps; among 25 weekly blobs it finds the week, not the change. Small commits can also be reverted individually: undo the constant change without losing the unrelated test you added the same day. Each commit costs ten seconds.
**Strongest case for the teammate.** Solo exploratory work whose intermediate states will never be cited or reverted, such as a scratch notebook, gains nothing from fine-grained history, and a coarse commit as a backup is acceptable there. The moment a number from that work is reported, the habit must change.
**Common mistakes.** Arguing from rules rather than from a question you need answered; not conceding anything.

## 3. A conflict in a bill of materials

**Resolution.**

```yaml
- item: load cell, 5 kg, strain gauge type
  quantity: 1
  search_terms: "5 kg load cell with HX711 amplifier board"
  note: one spare is not needed; the amplifier is the failure point
```

Justification: the item line is common to both. Quantity 1 is adopted because the incoming branch gives a *reason* (the failure point is the amplifier, so a spare cell adds cost without reducing risk) and HEAD gives none; if you disagreed, you would answer the reason, not just pick your own number. The longer search term is kept because it names the amplifier that the note refers to. The note is kept because a BOM that records why is the one that survives the next review. All three fence lines are deleted.
**Finishing the merge.** `git add projects/p00-first-experiment/bom.yaml` then `git commit` (Git proposes a merge message; keep or edit it). Porcelain before: `UU projects/p00-first-experiment/bom.yaml`; after `git add`: `M  projects/p00-first-experiment/bom.yaml` (resolved and staged; nothing further in the working tree).
**Timestamps.** Colleague: 2026-10-09T23:58:00−07:00 = 2026-10-10T06:58:00Z. You: 2026-10-10T08:30:00+01:00 = 2026-10-10T07:30:00Z. The colleague's change was made first, by 32 minutes, even though your clock face (08:30) reads earlier than theirs (23:58) on a different calendar day.
**Check.** `datetime.fromisoformat` on both strings and subtract: `timedelta(minutes=32)`. The file must parse as YAML after resolution (`python -c "import yaml, sys; yaml.safe_load(open('bom.yaml'))"`).
**Common mistakes.** Deleting one side wholesale because it is "theirs"; leaving a fence line in the file (YAML will fail to parse); forgetting `git add` after editing, so Git still shows `UU`.

## 4. Your own history

**Assumptions.** Results depend on the repository and on the day you run it; the method is what is graded. For a repository with N commits, expect to report N, the author table, per-rule counts, and the ordering flag.
**Result (form).** "N = 312 commits; 4 authors; `trailing-period` 3, `lowercase-start` 7, `vague` 2, `too-long` 1, `empty` 0, `surrounding-whitespace` 0; `is_newest_first` = False." Then the comparison: `git shortlog -sn` may show, say, `Ada Learner 180` and `A. Learner 12` where your parser shows the same two lines; both agree, and the fix is a `.mailmap` entry, not a parser change. If `shortlog` were run with `--no-merges`, its total would be smaller than yours by the number of merge commits; `git log --merges --oneline | wc -l` gives that number and closes the gap.
**First inversion.** Loop `for a, b in zip(commits, commits[1:])` and print the first pair with `a.date < b.date`. Typically it is a merge commit followed by the branch's commits, which were authored days earlier than the merge; `git log --graph --oneline` shows the branch joining. This is correct behaviour: author dates record when work was done, and `git log` orders by parent links, not by date.
**Check.** `sum(counts.values()) == len(commits)` must hold. Re-run with `%cI` instead of `%aI`: committer dates are far more often monotonic.
**Common mistakes.** Calling the inversion a parser bug; comparing against `shortlog` run on a different branch or with different flags.

## 5. Committer versus author

**Changes.** `LOG_FORMAT = "%H|%an|%aI|%cn|%cI|%s"`; `parts = line.split("|", maxsplit=5)`; `Commit` gains `committer: str` and `cdate: datetime`; both dates parsed with `fromisoformat` and checked for a timezone.
**Counts (form).** "Author ≠ committer: 41 of 312 (13 %). Committer date more than one day after author date: 27 of 312." Your numbers will differ.
**Interpretation.** The first count is the number of commits that were rebased, cherry-picked or squash-merged by someone other than the author: it measures how much of the history was rewritten on its way into `main`. The second count is contributions that waited in review for more than a day before being integrated: a measure of review latency. A repository with both counts near zero merges with plain merge commits and quickly; one with a high first count and low second count rebases aggressively.
**Check.** Every commit where the committer date differs from the author date must have *some* reason; on a few, `git show --format=fuller <sha>` will show the rebase explicitly.
**Common mistakes.** Forgetting that `maxsplit` must grow with the field count, which silently merges the committer date into the subject; comparing naive and aware datetimes (Python raises `TypeError`, which is the right outcome).

## 6. What goes in the repository

| item | where | reason |
|---|---|---|
| (a) 40 kB CSV of raw timings | repository | small, text, primary evidence; it must travel with the analysis that cites it |
| (b) 300 MB phone video | separate archive, SHA-256 recorded in the repository | binary and huge; Git would store every version in full; the hash preserves provenance |
| (c) PNG plot | nowhere in Git (ignored) | regenerated by (d) from (a); committing it creates binary conflicts and stale figures |
| (d) plotting script | repository | the text that produces the evidence; without it the plot is unverifiable |
| (e) `.env` with an API key | never | history is permanent and public; if it was ever committed, rotate the key today |
| (f) `__pycache__/` | ignored | compiled cache, machine-specific, regenerated on every run |
| (g) `.venv/` | ignored | hundreds of megabytes of installed packages; `requirements.txt` is the reproducible description of it |

**`.gitignore` lines.** `*.png` (or, more selectively, `figures/*.png`), `__pycache__/`, `.venv/`.
**Check.** After adding the lines, `git status` must no longer list the plot, the cache or the environment as untracked; `git check-ignore -v scratch_plot.png` prints the rule that matched.
**Common mistakes.** Committing the plot "so reviewers can see it" (put it in the PR description or the report instead); committing the key and then deleting it in a later commit, which removes nothing from history.

## Check-yourself answers

1. `git add` (working tree → staging area); `git commit` (staging area → repository).
2. Modified in the working tree, nothing staged. No, it will not be in the next commit.
3. Hash, author and date cannot contain `|`; a subject can. Putting it last and splitting with `maxsplit=3` keeps any pipes in the subject intact.
4. The first (09:15−05:00 = 14:15 UTC) happened later than the second (13:40+02:00 = 11:40 UTC), by 155 minutes.
5. It is generated from committed code and data, so it is redundant and goes stale; it is binary, so Git cannot diff or merge it and the repository grows with every regeneration.
6. The author wrote the change (name and date from `user.name` at commit time); the committer created the commit object that is in the history now. They differ after a rebase, cherry-pick, `git am`, or a squash-merge performed by someone else.
7. Skipping would silently lose a commit, and a summary based on a partial log would be wrong without any sign of it. Raising with the line number makes the corruption visible and locatable.
