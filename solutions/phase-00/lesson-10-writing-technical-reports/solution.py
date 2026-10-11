"""P00.10 Writing Technical Reports: reference solution.

A small linter for Markdown technical reports written to ``projects/_template/report-template.md``.
It checks three things a reviewer checks first: that every required section is present, that every
figure carries a numbered caption, and that every measured quantity in the Results section carries
an uncertainty.

Run:  python solution.py [path/to/report.md]
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-10-writing-technical-reports
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

# A Markdown ATX heading: one to six '#', a space, then the heading text.
_HEADING_RE = re.compile(r"^(#{1,6})\s+(.*?)\s*#*\s*$")
# A Markdown image: ![alt text](path "optional title")
_IMAGE_RE = re.compile(r"^!\[([^\]]*)\]\(([^)\s]+)(?:\s+\"[^\"]*\")?\)\s*$")
# A caption line: optional emphasis markers, then "Figure 3", "Figure 3.", or "Figure 3:".
_CAPTION_RE = re.compile(r"^[*_]*\s*Figure\s+(\d+)\s*[.:]?", re.IGNORECASE)

# Numbers and units for the uncertainty check.
_NUM = r"[-+]?\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
_UNIT = (
    r"(?:m/s\^?2|m/s²|m/s|km/h|rad/s|N·m|N\*m|kN|MN|MPa|kPa|GPa|Pa|kg|mg|g|mm|cm|km|m|"
    r"ms|min|s|Hz|kHz|°C|K|kJ|J|kW|W|mV|V|mA|A|%)"
)
_PM = r"(?:±|\+/-|\+-)"
# A unit must not be immediately followed by another word character (so "1 s" inside "Table 1 shows" does
# not match), and a number must not be glued to a preceding word character ("x2 kg" is not a number).
_END = r"(?![\w/²^°])"
_START = r"(?<![\w.])"
_QUANTITY_RE = re.compile(
    rf"{_START}(?:"
    rf"(?P<a>{_NUM}\s*{_PM}\s*{_NUM}\s*{_UNIT}{_END})"          # 9.81 ± 0.05 m/s²
    rf"|(?P<b>{_NUM}\s*{_UNIT}\s*{_PM}\s*{_NUM}(?:\s*{_UNIT})?{_END})"  # 9.81 m/s² ± 0.05 m/s²
    rf"|(?P<d>{_NUM}\(\d+\)\s*{_UNIT}{_END})"                  # 9.81(5) m/s²  (concise notation)
    rf"|(?P<c>{_NUM}\s*{_UNIT}{_END})"                         # 9.81 m/s²   <- bare, flagged
    rf")"
)


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
        m = _CAPTION_RE.match(self.caption)
        return int(m.group(1)) if m else None


def find_headings(text: str) -> list[tuple[int, int, str]]:
    """Return ``(line, level, heading_text)`` for every ATX heading in ``text``.

    ``line`` is 1-based; ``level`` is the number of leading '#'. Headings inside fenced code
    blocks (``` ... ```) are ignored.
    """
    out: list[tuple[int, int, str]] = []
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


def missing_headings(text: str) -> list[str]:
    """Return the required heading keywords that no heading of ``text`` contains.

    Matching is case-insensitive and by substring, so "## 3. Prediction (recorded before
    measurement)" satisfies the keyword "prediction".
    """
    present = " || ".join(h.lower() for _, _, h in find_headings(text))
    return [k for k in REQUIRED_HEADINGS if k not in present]


def section_lines(text: str, keyword: str) -> list[tuple[int, str]]:
    """Return ``(line, text)`` for the body lines of the first section whose heading contains ``keyword``.

    The body runs from the line after the heading to the line before the next heading of the
    same or higher level. Returns an empty list if no such heading exists.
    """
    lines = text.splitlines()
    headings = find_headings(text)
    for idx, (start, level, title) in enumerate(headings):
        if keyword.lower() in title.lower():
            end = len(lines) + 1
            for later_start, later_level, _ in headings[idx + 1:]:
                if later_level <= level:
                    end = later_start
                    break
            return [(n, lines[n - 1]) for n in range(start + 1, end)]
    return []


def find_figures(text: str) -> list[Figure]:
    """Return every Markdown image in ``text`` with its caption, if any.

    A caption is the first non-blank line after the image when that line starts with
    "Figure <n>" (optionally wrapped in ``*`` or ``_`` emphasis).
    """
    lines = text.splitlines()
    figures: list[Figure] = []
    for i, raw in enumerate(lines, start=1):
        m = _IMAGE_RE.match(raw.strip())
        if not m:
            continue
        caption: str | None = None
        for nxt in lines[i:]:
            if nxt.strip():
                if _CAPTION_RE.match(nxt.strip()):
                    caption = nxt.strip()
                break
        figures.append(Figure(line=i, path=m.group(2), alt=m.group(1), caption=caption))
    return figures


def figure_numbering_issues(figures: list[Figure]) -> list[str]:
    """Return one message per captioned figure whose number is not 1, 2, 3, ... in document order."""
    issues: list[str] = []
    captioned = [f for f in figures if f.number is not None]
    for expected, fig in enumerate(captioned, start=1):
        if fig.number != expected:
            issues.append(f"line {fig.line}: expected Figure {expected}, found Figure {fig.number}")
    return issues


def bare_quantities(text: str) -> list[str]:
    """Return every "number unit" token in ``text`` that carries no uncertainty.

    ``"9.81 ± 0.05 m/s²"``, ``"2.50 kg ± 0.01 kg"`` and ``"9.81(5) m/s²"`` are accepted;
    ``"9.81 m/s²"`` is returned. Numbers with no recognised unit ("Figure 3", "n = 20") are ignored.
    """
    return [m.group("c") for m in _QUANTITY_RE.finditer(text) if m.group("c")]


def lint_report(text: str) -> list[str]:
    """Lint a Markdown report and return a list of problems. An empty list means the report passes.

    Checks: required headings, a numbered caption under every figure, sequential figure
    numbering, and an uncertainty on every quantity with units in the Results section.
    """
    problems: list[str] = []
    for key in missing_headings(text):
        problems.append(f"missing heading: {key}")
    figures = find_figures(text)
    for fig in figures:
        if fig.caption is None:
            problems.append(f"line {fig.line}: figure '{fig.path}' has no 'Figure <n>' caption")
    problems.extend("figure numbering: " + msg for msg in figure_numbering_issues(figures))
    for line_no, line in section_lines(text, "results"):
        for q in bare_quantities(line):
            problems.append(f"results, line {line_no}: '{q}' has no uncertainty")
    return problems


def lint_file(path: str | Path) -> list[str]:
    """Read a Markdown file and return ``lint_report`` of its contents."""
    return lint_report(Path(path).read_text(encoding="utf-8"))


def main(argv: list[str] | None = None) -> int:
    """Lint the file named on the command line (default: the repo's report template) and print the result.

    Returns 0 when the report passes and 1 otherwise, so the linter can gate a commit hook.
    """
    args = sys.argv[1:] if argv is None else argv
    if args:
        target = Path(args[0])
    else:
        target = Path(__file__).resolve().parents[3] / "projects" / "_template" / "report-template.md"
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
