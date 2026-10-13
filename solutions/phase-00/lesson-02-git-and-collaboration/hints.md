# P00.02  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Column 1 is the staging area, column 2 is the working tree. A letter in column 1 means "this will be committed"; a letter in column 2 means "there is a change that will not".
2. Ask what `git log --oneline` looks like after six months under each habit, and what `git bisect` (which finds the commit that introduced a change by halving the history) can do with 300 small commits versus 25 big ones.
3. The two sides agree on the item; they disagree on quantity and description, and one side adds a note explaining *why*. Resolution is a decision, not a mechanical merge. For the time question, convert both to UTC.
4. If the author tables differ, look first at names (same person, two spellings), then at merges (`git log --merges`), then at author versus committer.
5. `split("|", maxsplit=5)` now; and the one-day rule is `(c_date - a_date) > timedelta(days=1)` on aware datetimes.
6. Ask of each item: is it text written by a human? Can it be regenerated? Would leaking it hurt? Does it grow with every change?
