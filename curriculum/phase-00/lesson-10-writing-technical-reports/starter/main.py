"""P00.10 Writing Technical Reports: starter code.

A small linter for Markdown technical reports written to ``projects/_template/report-template.md``.
It checks three things a reviewer checks first: that every required section is present, that every
figure carries a numbered caption, and that every measured quantity in the Results section carries
an uncertainty.

Complete every function marked NotImplementedError. When finished, running this file lints a report.

Run:  python main.py [path/to/report.md]
Test: pytest curriculum/phase-00/lesson-10-writing-technical-reports
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

# Keywords that must appear (case-insensitively) in a heading of the report, in report order.
# They mirror the sections of projects/_template/report-template.md.
REQUIRED_HEADINGS: tuple[str, ...] = (
    "summary",
    "objective",
    "theory",
    "prediction",
    "apparatus",
    "results",
    "discussion",
    "conclusions",
    "reproducibility",
)

# Building blocks for your regular expressions (part 3 of the lesson explains each one).
_NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"          # 12, -3.5, 1.2e3
_UNIT = (
    r"(?:m/s\^?2|m/s²|m/s|km/h|rad/s|N·m|N\*m|kN|MN|MPa|kPa|GPa|Pa|kg|mg|g|mm|cm|km|m|"
    r"ms|min|s|Hz|kHz|°C|K|kJ|J|kW|W|mV|V|mA|A|%)"
)
_PM = r"(?:±|\+/-|\+-)"                                  # the three ways people write "plus or minus"


@dataclass(frozen=True)
class Figure:
    """One Markdown image found in the report.

    ``line`` is 1-based. ``caption`` is the caption line's text, or ``None`` when the next
    non-blank line after the image does not start with "Figure <n>". ``number`` is the parsed
    figure number, or ``None`` when there is no caption.
    """

    line: int
    path: str
    alt: str
    caption: str | None

    @property
    def number(self) -> int | None:
        if self.caption is None:
            return None
        m = re.match(r"^[*_]*\s*Figure\s+(\d+)", self.caption, re.IGNORECASE)
        return int(m.group(1)) if m else None


def find_headings(text: str) -> list[tuple[int, int, str]]:
    """Return ``(line, level, heading_text)`` for every ATX heading in ``text``.

    ``line`` is 1-based; ``level`` is the number of leading '#'. Headings inside fenced code
    blocks (``` ... ```) are ignored.
    """
    raise NotImplementedError


def missing_headings(text: str) -> list[str]:
    """Return the required heading keywords that no heading of ``text`` contains.

    Matching is case-insensitive and by substring, so "## 3. Prediction (recorded before
    measurement)" satisfies the keyword "prediction".
    """
    raise NotImplementedError


def section_lines(text: str, keyword: str) -> list[tuple[int, str]]:
    """Return ``(line, text)`` for the body lines of the first section whose heading contains ``keyword``.

    The body runs from the line after the heading to the line before the next heading of the
    same or higher level. Returns an empty list if no such heading exists.
    """
    raise NotImplementedError


def find_figures(text: str) -> list[Figure]:
    """Return every Markdown image in ``text`` with its caption, if any.

    A caption is the first non-blank line after the image when that line starts with
    "Figure <n>" (optionally wrapped in ``*`` or ``_`` emphasis).
    """
    raise NotImplementedError


def figure_numbering_issues(figures: list[Figure]) -> list[str]:
    """Return one message per captioned figure whose number is not 1, 2, 3, ... in document order."""
    raise NotImplementedError


def bare_quantities(text: str) -> list[str]:
    """Return every "number unit" token in ``text`` that carries no uncertainty.

    ``"9.81 ± 0.05 m/s²"``, ``"2.50 kg ± 0.01 kg"`` and ``"9.81(5) m/s²"`` are accepted;
    ``"9.81 m/s²"`` is returned. Numbers with no recognised unit ("Figure 3", "n = 20") are ignored.
    """
    raise NotImplementedError


def lint_report(text: str) -> list[str]:
    """Lint a Markdown report and return a list of problems. An empty list means the report passes.

    Checks: required headings, a numbered caption under every figure, sequential figure
    numbering, and an uncertainty on every quantity with units in the Results section.
    """
    raise NotImplementedError


def lint_file(path: str | Path) -> list[str]:
    """Read a Markdown file and return ``lint_report`` of its contents."""
    raise NotImplementedError


def main(argv: list[str] | None = None) -> int:
    """Lint the file named on the command line (default: the repo's report template) and print the result.

    Returns 0 when the report passes and 1 otherwise, so the linter can gate a commit hook.
    """
    args = sys.argv[1:] if argv is None else argv
    if args:
        target = Path(args[0])
    else:
        target = Path(__file__).resolve().parents[4] / "projects" / "_template" / "report-template.md"
    problems = lint_file(target)
    print(f"Linting {target}")
    if not problems:
        print("OK: all required sections present, every figure captioned, every result has an uncertainty.")
        return 0
    for p in problems:
        print(" -", p)
    print(f"{len(problems)} problem(s) found.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
