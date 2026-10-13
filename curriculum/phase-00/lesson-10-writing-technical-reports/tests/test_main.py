"""Public tests for P00.10 Writing Technical Reports.

They check the linter's behaviour on small fixture reports, not one magic output string.

Run `pytest curriculum/phase-00/lesson-10-writing-technical-reports` against your starter/main.py, or
`AEFS_IMPL=solution pytest ...` against the reference solution.
"""
from pathlib import Path

import pytest

from conftest import load_impl

impl = load_impl(__file__)

GOOD_REPORT = """\
# Pendulum period: Technical report

## Summary (5 sentences)
We measured g and got 9.79 ± 0.06 m/s².

## 1. Objective and requirements
Predict the period to within 2 %.

## 2. Theory and assumptions
Small angles: T = 2π sqrt(L/g). The 1.000 m string is assumed inextensible.

## 3. Prediction (recorded before measurement)
T = 2.006 ± 0.004 s from L = 1.000 ± 0.002 m.

## 4. Apparatus, calibration, and safety
Tape measure, phone stopwatch.

## 5. Results
Mean period 2.010 ± 0.005 s over n = 20 swings (Figure 1). Length 1.000 m ± 0.002 m.
Air temperature 21.5(5) °C.

![Period against length](period_vs_length.png)
*Figure 1. Measured period against string length; error bars are one standard deviation.*

![Residuals](residuals.png)

Figure 2: Residuals of the fit in Figure 1.

## 6. Discussion
Agreement within 1 standard deviation.

## 7. Conclusions and redesign
Use a longer string.

## 8. Reproducibility
Commit abc1234; data in data/raw.csv.
"""


def _without(section_keyword: str) -> str:
    """GOOD_REPORT with every heading containing ``section_keyword`` demoted to plain text."""
    return "\n".join(
        line.lstrip("# ") if line.startswith("#") and section_keyword in line.lower() else line
        for line in GOOD_REPORT.splitlines()
    )


# --- headings ---------------------------------------------------------------------------------

def test_good_report_has_no_missing_headings():
    assert impl.missing_headings(GOOD_REPORT) == []


@pytest.mark.parametrize("keyword", ["prediction", "reproducibility", "apparatus"])
def test_missing_heading_is_reported_by_keyword(keyword):
    missing = impl.missing_headings(_without(keyword))
    assert missing == [keyword]


def test_heading_match_is_case_insensitive_and_ignores_numbering():
    text = "# Report\n## 3. PREDICTION (recorded before measurement)\nbody\n"
    assert "prediction" not in impl.missing_headings(text)


def test_headings_inside_code_fences_are_ignored():
    text = "# Report\n```\n## Summary\n```\n"
    heads = impl.find_headings(text)
    assert [h[2] for h in heads] == ["Report"]
    assert "summary" in impl.missing_headings(text)


def test_find_headings_reports_line_and_level():
    text = "intro\n# Title\n\n### Deep\n"
    assert impl.find_headings(text) == [(2, 1, "Title"), (4, 3, "Deep")]


def test_repo_template_satisfies_the_heading_rule():
    root = Path(__file__).resolve().parents[4]
    template = root / "projects" / "_template" / "report-template.md"
    assert impl.missing_headings(template.read_text(encoding="utf-8")) == []


# --- sections -----------------------------------------------------------------------------------

def test_section_lines_stop_at_next_heading_of_same_level():
    text = "## A\na1\na2\n### A.1\na3\n## B\nb1\n"
    got = impl.section_lines(text, "a")
    assert [t for _, t in got] == ["a1", "a2", "### A.1", "a3"]
    assert got[0][0] == 2  # 1-based line of 'a1'


def test_section_lines_missing_section_is_empty():
    assert impl.section_lines("## Only\nx\n", "results") == []


# --- figures ----------------------------------------------------------------------------------

def test_figures_found_with_and_without_captions():
    figs = impl.find_figures(GOOD_REPORT)
    assert [f.path for f in figs] == ["period_vs_length.png", "residuals.png"]
    assert [f.number for f in figs] == [1, 2]
    assert all(f.caption is not None for f in figs)


def test_uncaptioned_figure_is_reported():
    text = "## 5. Results\n\n![plot](plot.png)\n\nThe plot shows the trend.\n"
    figs = impl.find_figures(text)
    assert len(figs) == 1 and figs[0].caption is None and figs[0].number is None
    problems = impl.lint_report(text)
    assert any("plot.png" in p and "caption" in p for p in problems)


def test_alt_text_alone_is_not_a_caption():
    text = "![Figure 1. A plot](plot.png)\n\nSome prose.\n"
    assert impl.find_figures(text)[0].caption is None


def test_figure_numbering_must_be_sequential():
    text = (
        "![a](a.png)\n*Figure 2. Second first.*\n\n"
        "![b](b.png)\n*Figure 1. First second.*\n"
    )
    issues = impl.figure_numbering_issues(impl.find_figures(text))
    assert issues, "swapped figure numbers must be reported"
    assert "Figure 1" in issues[0] and "Figure 2" in issues[0] and "line 1" in issues[0]
    skipped = text.replace("Figure 2.", "Figure 1.").replace("Figure 1. First", "Figure 3. First")
    assert len(impl.figure_numbering_issues(impl.find_figures(skipped))) == 1
    assert impl.figure_numbering_issues(impl.find_figures(GOOD_REPORT)) == []


# --- uncertainties ------------------------------------------------------------------------------

@pytest.mark.parametrize(
    "snippet",
    [
        "9.81 ± 0.05 m/s²",
        "9.81 +/- 0.05 m/s^2",
        "2.50 kg ± 0.01 kg",
        "2.50 kg ± 0.01",
        "21.5(5) °C",
        "1.2e3 ± 50 Pa",
    ],
)
def test_quantities_with_uncertainty_are_accepted(snippet):
    assert impl.bare_quantities(f"The value was {snippet} today.") == []


@pytest.mark.parametrize(
    "snippet, expected",
    [
        ("9.81 m/s²", "9.81 m/s²"),
        ("a mass of 2.5 kg", "2.5 kg"),
        ("pressure 101325 Pa", "101325 Pa"),
        ("efficiency 45 %", "45 %"),
    ],
)
def test_bare_quantities_are_reported(snippet, expected):
    assert impl.bare_quantities(snippet) == [expected]


def test_numbers_without_units_are_ignored():
    text = "See Figure 3 and Table 1 shows n = 20 runs at Mach 0.8; Equation (2) applies; x = 0.25c."
    assert impl.bare_quantities(text) == []


def test_bare_quantity_in_results_is_flagged_but_theory_is_not():
    good = impl.lint_report(GOOD_REPORT)
    assert good == [], good  # the 1.000 m string in Theory is allowed; Results are all ± or (5)
    bad = GOOD_REPORT.replace("Mean period 2.010 ± 0.005 s", "Mean period 2.010 s")
    problems = impl.lint_report(bad)
    assert len(problems) == 1
    assert "2.010 s" in problems[0] and "results" in problems[0].lower()


# --- whole-file behaviour -----------------------------------------------------------------------

def test_lint_file_matches_lint_report(tmp_path):
    p = tmp_path / "report.md"
    bad = _without("discussion").replace("![Residuals](residuals.png)\n\nFigure 2:", "![Residuals](residuals.png)\n\nResiduals:")
    p.write_text(bad, encoding="utf-8")
    from_file = impl.lint_file(p)
    assert from_file == impl.lint_report(bad)
    assert any("discussion" in m for m in from_file)
    assert any("residuals.png" in m for m in from_file)


def test_main_exit_code(tmp_path, capsys):
    good = tmp_path / "good.md"
    good.write_text(GOOD_REPORT, encoding="utf-8")
    assert impl.main([str(good)]) == 0
    bad = tmp_path / "bad.md"
    bad.write_text(_without("summary"), encoding="utf-8")
    assert impl.main([str(bad)]) == 1
    out = capsys.readouterr().out
    assert "summary" in out
