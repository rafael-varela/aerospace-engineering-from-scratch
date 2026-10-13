"""P00.02 Git & Collaboration: starter code.

Complete every function marked NotImplementedError. The finished script holds small tools
for reading what Git tells you: a parser for ``git log`` output in a fixed format, a
parser for ``git status --porcelain``, a commit-subject linter that encodes this repository's
commit-message rule, and a few summaries. Everything works on strings, so it runs on a fixture
without needing Git installed.

Run:  python main.py            (summarises the built-in fixture)
      git log --pretty=format:'%H|%an|%aI|%s' | python main.py -   (summarise a real log)
Test: pytest curriculum/phase-00/lesson-02-git-and-collaboration   (from the repo root)
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
    raise NotImplementedError


def commits_per_author(commits: list[Commit]) -> dict[str, int]:
    """Count commits by author name, most commits first (ties in first-seen order)."""
    raise NotImplementedError


def is_newest_first(commits: list[Commit]) -> bool:
    """True when author dates are non-increasing down the list, as ``git log`` prints them.

    Dates carry their own UTC offsets, so comparing them compares instants, not clock faces.
    """
    raise NotImplementedError


def lint_subject(subject: str) -> list[str]:
    """Check one commit subject against this repository's rule: one plain sentence.

    Returns a list of problem codes (empty when the subject is fine):
    ``"empty"``, ``"too-long"`` (over SUBJECT_MAX characters), ``"trailing-period"``,
    ``"surrounding-whitespace"``, ``"lowercase-start"``, ``"vague"`` (a bare word such as
    'wip', 'fix', 'update', 'changes', 'stuff').
    """
    raise NotImplementedError


def parse_status_porcelain(text: str) -> dict[str, list[str]]:
    """Classify the lines of ``git status --porcelain`` (version 1 format).

    Each line is ``XY path`` where X is the index (staged) state and Y the working-tree state.
    Returns a dict with keys ``"staged"``, ``"unstaged"``, ``"untracked"``, ``"conflicted"``;
    a file with both a staged and an unstaged change appears in both lists. Order is preserved.
    """
    raise NotImplementedError


def summarize(commits: list[Commit]) -> str:
    """Render a short human-readable summary of a commit list."""
    raise NotImplementedError


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
