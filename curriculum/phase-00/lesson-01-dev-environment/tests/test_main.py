"""Public tests for P00.01. They assert behaviour, not one magic version number.

Run `pytest curriculum/phase-00/lesson-01-dev-environment` against your starter/main.py, or
with `AEFS_IMPL=solution` against the reference solution. The logic is tested with injected
fake version lookups; one test then runs the check against the real environment, because
the whole point of the lesson is that your environment passes it.
"""
import pytest

from conftest import load_impl

impl = load_impl(__file__)


# --- parse_version ------------------------------------------------------------------------
@pytest.mark.parametrize(
    "text, expected",
    [("1.26.4", (1, 26, 4)), ("3.13", (3, 13)), ("2.0.0rc1", (2, 0, 0)),
     ("0.23.post1", (0, 23)), (" 9.1.1", (9, 1, 1)), ("10", (10,))],
)
def test_parse_version(text, expected):
    assert impl.parse_version(text) == expected


def test_parse_version_rejects_garbage():
    with pytest.raises(ValueError):
        impl.parse_version("latest")


# --- meets_minimum: numeric, not lexicographic ------------------------------------------
def test_meets_minimum_equal_and_padding():
    assert impl.meets_minimum("1.26", "1.26.0")
    assert impl.meets_minimum("1.26.0", "1.26")


def test_meets_minimum_is_numeric_not_string_order():
    assert impl.meets_minimum("1.26.4", "1.9")          # "1.26.4" < "1.9" as strings: wrong
    assert not impl.meets_minimum("1.9", "1.26")
    assert impl.meets_minimum("10.0", "9.9")


def test_meets_minimum_is_a_total_order_on_examples():
    versions = ["0.9", "1.0", "1.0.1", "1.10", "2.0rc1", "2.0"]
    for i, lo in enumerate(versions):
        for hi in versions[i:]:
            assert impl.meets_minimum(hi, lo), (hi, lo)


# --- check_packages with an injected lookup --------------------------------------------
def test_check_packages_flags_each_case():
    fake = {"numpy": "1.26.4", "scipy": "1.9.3", "pint": None}
    reqs = {"numpy": "1.26", "scipy": "1.11", "pint": "0.23"}
    results = impl.check_packages(reqs, lookup=fake.get)
    by_name = {r.name: r for r in results}
    assert [r.name for r in results] == list(reqs)
    assert by_name["numpy"].ok and by_name["numpy"].installed == "1.26.4"
    assert not by_name["scipy"].ok, "1.9.3 is below 1.11"
    assert not by_name["pint"].ok and by_name["pint"].installed is None
    assert by_name["pint"].minimum == "0.23"


def test_check_python_with_injected_versions():
    assert impl.check_python((3, 11, 0), (3, 11)).ok
    assert impl.check_python((3, 13, 2), (3, 11)).ok
    assert not impl.check_python((3, 10, 14), (3, 11)).ok
    assert impl.check_python((3, 12, 1), (3, 11)).installed == "3.12.1"


# --- The real environment -----------------------------------------------------------------
def test_real_environment_passes():
    results = [impl.check_python(), *impl.check_packages()]
    failing = [r for r in results if not r.ok]
    assert not failing, f"environment check failed for: {failing}"


def test_installed_version_handles_missing_package():
    assert impl.installed_version("numpy") is not None
    assert impl.installed_version("definitely-not-a-real-package-aefs") is None


# --- Smoke test and report ----------------------------------------------------------------
def test_smoke_speed_of_sound():
    a0 = impl.speed_of_sound_smoke()
    assert isinstance(a0, float)
    assert a0 == pytest.approx(340.29, abs=0.05)


def test_format_report_mentions_every_component_and_status():
    results = [
        impl.CheckResult("numpy", "1.26.4", "1.26", True),
        impl.CheckResult("pint", None, "0.23", False),
    ]
    text = impl.format_report(results)
    assert "numpy" in text and "pint" in text
    assert "OK" in text and "FAIL" in text
    assert "missing" in text.lower()
    assert len(text.splitlines()) == 3  # header + one line per result
