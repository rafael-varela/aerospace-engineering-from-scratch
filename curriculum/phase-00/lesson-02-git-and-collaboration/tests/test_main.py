"""Public tests for P00.02. They run on fixture strings; Git itself is not required.

Run `pytest curriculum/phase-00/lesson-02-git-and-collaboration` against your starter/main.py,
or with `AEFS_IMPL=solution` against the reference solution.
"""
from datetime import datetime, timedelta, timezone

import pytest

from conftest import load_impl

impl = load_impl(__file__)


# --- parse_git_log -------------------------------------------------------------------------
def test_fixture_parses_to_one_commit_per_nonblank_line():
    commits = impl.parse_git_log(impl.FIXTURE_LOG)
    nonblank = [ln for ln in impl.FIXTURE_LOG.splitlines() if ln.strip()]
    assert len(commits) == len(nonblank) == 6
    assert all(len(c.sha) == 40 for c in commits)


def test_fields_have_the_right_types_and_timezones():
    c = impl.parse_git_log(impl.FIXTURE_LOG)[0]
    assert isinstance(c.date, datetime)
    assert c.date.tzinfo is not None, "dates must be timezone-aware"
    assert c.date.utcoffset().total_seconds() == 2 * 3600
    assert c.author == "Ada Learner"
    assert c.subject.startswith("Add ISA density check")


def test_subject_may_contain_the_separator():
    commits = impl.parse_git_log(impl.FIXTURE_LOG)
    with_pipe = [c for c in commits if "|" in c.subject]
    assert len(with_pipe) == 1
    assert with_pipe[0].subject == "Implement vectorized pressure() | gradient and isothermal branches"


def test_blank_lines_are_skipped_and_order_is_preserved():
    text = "\n\n" + impl.FIXTURE_LOG.replace("\n", "\n\n") + "\n"
    a = impl.parse_git_log(text)
    b = impl.parse_git_log(impl.FIXTURE_LOG)
    assert [c.sha for c in a] == [c.sha for c in b]


@pytest.mark.parametrize(
    "bad",
    [
        "9f1c2d3e4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d|Ada Learner|2026-10-11T14:02:33+02:00",  # 3 fields
        "notahash|Ada Learner|2026-10-11T14:02:33+02:00|Subject",                           # bad sha
        "9f1c2d3e4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d|Ada Learner|yesterday|Subject",             # bad date
        "9f1c2d3e4a5b6c7d8e9f0a1b2c3d4e5f6a7b8c9d|Ada Learner|2026-10-11T14:02:33|Subject",   # no tz
    ],
)
def test_malformed_lines_raise_value_error(bad):
    with pytest.raises(ValueError):
        impl.parse_git_log(bad)


def test_empty_log_gives_empty_list():
    assert impl.parse_git_log("") == []
    assert impl.commits_per_author([]) == {}
    assert impl.is_newest_first([])


# --- Summaries -----------------------------------------------------------------------------
def test_commits_per_author_counts_and_sorts():
    commits = impl.parse_git_log(impl.FIXTURE_LOG)
    counts = impl.commits_per_author(commits)
    assert counts == {"Ada Learner": 4, "Ben Reviewer": 2}
    assert list(counts) == ["Ada Learner", "Ben Reviewer"]
    assert sum(counts.values()) == len(commits)


def test_is_newest_first_compares_instants_not_clock_faces():
    commits = impl.parse_git_log(impl.FIXTURE_LOG)
    assert impl.is_newest_first(commits)
    assert not impl.is_newest_first(list(reversed(commits)))
    # 10:00 at UTC-5 is 15:00 UTC: later than 14:00 UTC even though the clock face is smaller.
    later = impl.Commit("a" * 40, "x", datetime(2026, 1, 1, 10, 0, tzinfo=timezone(timedelta(hours=-5))), "s")
    earlier = impl.Commit("b" * 40, "x", datetime(2026, 1, 1, 14, 0, tzinfo=timezone.utc), "s")
    assert later.date > earlier.date
    assert impl.is_newest_first([later, earlier])
    assert not impl.is_newest_first([earlier, later])


# --- lint_subject --------------------------------------------------------------------------
def test_good_subject_has_no_problems():
    assert impl.lint_subject("Add ISA density check to lesson 03 tests") == []


@pytest.mark.parametrize(
    "subject, code",
    [
        ("", "empty"),
        ("   ", "empty"),
        ("wip", "vague"),
        ("Fix", "vague"),
        ("Document the geopotential conversion in lesson 03.", "trailing-period"),
        ("add density test", "lowercase-start"),
        (" Add density test", "surrounding-whitespace"),
        ("A" * (impl.SUBJECT_MAX + 1), "too-long"),
    ],
)
def test_lint_subject_detects_each_rule(subject, code):
    assert code in impl.lint_subject(subject)


def test_lint_subject_length_boundary():
    assert "too-long" not in impl.lint_subject("A" * impl.SUBJECT_MAX)
    assert "too-long" in impl.lint_subject("A" * (impl.SUBJECT_MAX + 1))


def test_fixture_has_exactly_two_flagged_subjects():
    commits = impl.parse_git_log(impl.FIXTURE_LOG)
    flagged = {c.subject for c in commits if impl.lint_subject(c.subject)}
    assert flagged == {"wip", "Document the geopotential conversion in lesson 03."}


# --- parse_status_porcelain ------------------------------------------------------------------
def test_status_classification():
    st = impl.parse_status_porcelain(impl.FIXTURE_STATUS)
    assert st["untracked"] == ["scratch_plot.png", "notes.txt"]
    assert st["conflicted"] == ["projects/p00-first-experiment/bom.yaml"]
    assert "curriculum/phase-00/lesson-03-numpy-engineering-calculator/starter/main.py" in st["staged"]
    assert "curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md" in st["unstaged"]
    assert "curriculum/phase-00/lesson-03-numpy-engineering-calculator/lesson.md" not in st["staged"]
    assert "curriculum/phase-00/lesson-03-numpy-engineering-calculator/tests/test_extra.py" in st["staged"]


def test_status_file_with_both_changes_appears_in_both_lists():
    st = impl.parse_status_porcelain("MM a.py\n")
    assert st["staged"] == ["a.py"] and st["unstaged"] == ["a.py"]
    assert st["untracked"] == [] and st["conflicted"] == []


def test_status_every_line_lands_somewhere():
    st = impl.parse_status_porcelain(impl.FIXTURE_STATUS)
    paths = {p for lst in st.values() for p in lst}
    expected = {ln[3:] for ln in impl.FIXTURE_STATUS.splitlines() if ln.strip()}
    assert paths == expected
    assert impl.parse_status_porcelain("") == {"staged": [], "unstaged": [], "untracked": [], "conflicted": []}
