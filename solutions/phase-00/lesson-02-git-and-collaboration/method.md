# P00.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. For each line split off the two-letter code and the path. 2. Classify by the rules in §2.1 of the lesson (and the `conflicted` set for `UU`). 3. For each, say whether column 1 is non-blank (goes in the commit). 4. For the command list: the lesson text is unstaged, so stage only it and commit; but note the conflicted file blocks the commit until resolved or the merge is aborted.

## 2.
1. State the claim. 2. Scenario: the density changes; under weekly commits, `git diff` between two weekly snapshots shows hundreds of unrelated lines; under small commits, `git log --oneline -- main.py` shows one sentence naming the gas-constant change. 3. Mention bisecting and reverting one change without losing others. 4. Concede: a solo exploratory notebook with nothing cited may be committed coarsely.

## 3.
1. Keep the item line. 2. Choose quantity: the note gives a reason for 1 that the other side lacks; adopt it unless you have a counter-reason. 3. Choose the more descriptive search term. 4. Keep the note. 5. Remove the fence lines. 6. `git add bom.yaml`, `git commit`. 7. Porcelain before: `UU`; after add: `M ` (merge resolved, staged). 8. Convert both timestamps to UTC and subtract.

## 4.
1. `git log --pretty=format:'%H|%an|%aI|%s' > log.txt`. 2. `python main.py - < log.txt`. 3. Count lint codes with a `Counter` over all flagged codes. 4. `git shortlog -sn`; line up names; attribute each gap to one of the four causes. 5. If not newest-first, loop over adjacent pairs and print the first inversion; inspect with `git log --graph --oneline`.

## 5.
1. Change `LOG_FORMAT` and the `maxsplit`; add two fields to `Commit` (or a subclass). 2. Parse both dates as aware datetimes. 3. Count `c.author != c.committer`. 4. Count `(c.cdate - c.date) > timedelta(days=1)`. 5. Interpret: the first count is rebased/squashed commits; the second is contributions that waited in review.

## 6.
1. Apply the four questions from the hint to each item. 2. Small text data: commit. Large binary: archive and cite its SHA-256. Generated: ignore. Secret: never, and rotate it if it was ever committed. Environment: ignore. 3. Write three `.gitignore` patterns.
