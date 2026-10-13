"""P00.02 Git & Collaboration: reference solution.

Small tools for reading what Git tells you: a parser for ``git log`` output in a fixed format, a
parser for ``git status --porcelain``, a commit-subject linter that encodes this repository's
commit-message rule, and a few summaries. Everything works on strings, so it runs on a fixture
without needing Git installed.

Run:  python solution.py            (summarises the built-in fixture)
      git log --pretty=format:'%H|%an|%aI|%s' | python solution.py -   (summarise a real log)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-02-git-and-collaboration
"""
from __future__ import annotations

import sys
from collections import Counter
from dataclasses import dataclass
from datetime import datetime

# The format string to give `git log --pretty=format:...`. Fields: full hash, author name,
# author date in strict ISO 8601, subject line. '|' separates fields; the subject comes last so
# that a '|' inside a subject is harmless.
LOG_FORMAT: str = "%H|%an|%aI|%s"
SUBJECT_MAX: int = 72

# A fixture: what `git log --pretty=format:'%H|%an|%aI|%s'` printed on a feature branch of this
# repository after a learner's first contribution. Newest commit first, as Git prints it.
FIXTURE_LOG: str = """\
9f1c2d3e4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d|Ada Learner|2026-10-11T14:02:33+02:00|Add ISA density check to lesson 03 tests
8e0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6f7a8b|Ada Learner|2026-10-11T13:40:10+02:00|Fix tropopause pressure tolerance in test_main.py
7d9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6f7a|Ben Reviewer|2026-10-10T09:15:00-05:00|Merge pull request #41 from ada/lesson-03-isa
6c8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e6f|Ada Learner|2026-10-09T21:05:47+02:00|wip
5b7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d5e|Ada Learner|2026-10-09T20:48:12+02:00|Implement vectorized pressure() | gradient and isothermal branches
4a6d7e8f9a0b1c2d3e4f5a6b7c8d9e0f1a2b3c4d|Ben Reviewer|2026-10-08T16:30:00-05:00|Document the geopotential conversion in lesson 03.
"""

FIXTURE_STATUS: str = """\
M  curriculum/phase-00/lesson-03-numpy-engineering-calculator/starter/main.py
 M curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md
MM solutions/phase-00/lesson-03-numpy-engineering-calculator/solution.py
A  curriculum/phase-00/lesson-03-numpy-engineering-calculator/tests/test_extra.py
?? scratch_plot.png
?? notes.txt
UU projects/p00-first-experiment/bom.yaml
"""


@dataclass(frozen=True)
class Commit:
    """One commit as reported by ``git log --pretty=format:LOG_FORMAT``."""

    sha: str
    author: str
    date: datetime      # timezone-aware
    subject: str


def parse_git_log(text: str) -> list[Commit]:
    """Parse the output of ``git log --pretty=format:LOG_FORMAT`` into Commit objects.

    Blank lines are skipped. Each remaining line must have at least four '|'-separated fields;
    the subject is everything after the third '|', so a subject may itself contain '|'.
    The date is parsed as timezone-aware ISO 8601. Raise ValueError for a malformed line.
    """
    commits: list[Commit] = []
    for lineno, raw in enumerate(text.splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        parts = line.split("|", maxsplit=3)
        if len(parts) != 4:
            raise ValueError(f"line {lineno}: expected 4 '|'-separated fields, got {len(parts)}: {raw!r}")
        sha, author, date_text, subject = parts
        if len(sha) != 40 or any(c not in "0123456789abcdef" for c in sha):
            raise ValueError(f"line {lineno}: {sha!r} is not a 40-character hex hash")
        try:
            date = datetime.fromisoformat(date_text)
        except ValueError as exc:
            raise ValueError(f"line {lineno}: bad ISO date {date_text!r}") from exc
        if date.tzinfo is None:
            raise ValueError(f"line {lineno}: date {date_text!r} has no timezone offset")
        commits.append(Commit(sha, author, date, subject))
    return commits


def commits_per_author(commits: list[Commit]) -> dict[str, int]:
    """Count commits by author name, most commits first (ties in first-seen order)."""
    counts = Counter(c.author for c in commits)
    return dict(counts.most_common())


def is_newest_first(commits: list[Commit]) -> bool:
    """True when author dates are non-increasing down the list, as ``git log`` prints them.

    Dates carry their own UTC offsets, so comparing them compares instants, not clock faces.
    """
    return all(a.date >= b.date for a, b in zip(commits, commits[1:]))


def lint_subject(subject: str) -> list[str]:
    """Check one commit subject against this repository's rule: one plain sentence.

    Returns a list of problem codes (empty when the subject is fine):
    ``"empty"``, ``"too-long"`` (over SUBJECT_MAX characters), ``"trailing-period"``,
    ``"surrounding-whitespace"``, ``"lowercase-start"``, ``"vague"`` (a bare word such as
    'wip', 'fix', 'update', 'changes', 'stuff').
    """
    problems: list[str] = []
    if subject.strip() == "":
        return ["empty"]
    if subject != subject.strip():
        problems.append("surrounding-whitespace")
    s = subject.strip()
    if len(s) > SUBJECT_MAX:
        problems.append("too-long")
    if s.endswith("."):
        problems.append("trailing-period")
    if s[0].isalpha() and not s[0].isupper():
        problems.append("lowercase-start")
    if s.lower().rstrip(".!") in {"wip", "fix", "fixes", "update", "updates", "changes", "stuff", "misc", "test", "tmp"}:
        problems.append("vague")
    return problems


def parse_status_porcelain(text: str) -> dict[str, list[str]]:
    """Classify the lines of ``git status --porcelain`` (version 1 format).

    Each line is ``XY path`` where X is the index (staged) state and Y the working-tree state.
    Returns a dict with keys ``"staged"``, ``"unstaged"``, ``"untracked"``, ``"conflicted"``;
    a file with both a staged and an unstaged change appears in both lists. Order is preserved.
    """
    out: dict[str, list[str]] = {"staged": [], "unstaged": [], "untracked": [], "conflicted": []}
    conflict_codes = {"DD", "AU", "UD", "UA", "DU", "AA", "UU"}
    for raw in text.splitlines():
        if len(raw) < 4:
            continue
        xy, path = raw[:2], raw[3:]
        if xy == "??":
            out["untracked"].append(path)
        elif xy in conflict_codes:
            out["conflicted"].append(path)
        else:
            if xy[0] not in " ?!":
                out["staged"].append(path)
            if xy[1] not in " ?!":
                out["unstaged"].append(path)
    return out


def summarize(commits: list[Commit]) -> str:
    """Render a short human-readable summary of a commit list."""
    lines = [f"{len(commits)} commits, newest first: {is_newest_first(commits)}"]
    for author, n in commits_per_author(commits).items():
        lines.append(f"  {author}: {n}")
    flagged = [(c, lint_subject(c.subject)) for c in commits]
    flagged = [(c, p) for c, p in flagged if p]
    lines.append(f"{len(flagged)} subject line(s) need attention:")
    for c, p in flagged:
        lines.append(f"  {c.sha[:7]}  {c.subject!r}: {', '.join(p)}")
    return "\n".join(lines)


def main() -> None:
    """Summarise the fixture log, or a real log read from stdin when called with '-'."""
    text = sys.stdin.read() if len(sys.argv) > 1 and sys.argv[1] == "-" else FIXTURE_LOG
    commits = parse_git_log(text)
    print(summarize(commits))
    print()
    status = parse_status_porcelain(FIXTURE_STATUS)
    print("fixture `git status --porcelain`:")
    for key, paths in status.items():
        print(f"  {key:<11}{len(paths)}  {paths}")


if __name__ == "__main__":
    main()
