"""Public tests for P00.13. They check the checker's behaviour on fixture strings, not magic numbers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. No files are read or written.
"""
import pytest

from conftest import load_impl

impl = load_impl(__file__)

COMPLETE = """\
| | |
|---|---|
| **Hazards** | falling mass; pinch between clamp and bench; eye (snapping rubber band). |
| **Energy limits** | mass 0.5 kg, drop height 1.2 m, 9 V battery, 0.3 A, 3 m/s. |
| **Controls** | safety glasses; 1 m exclusion circle; clamp to bench; mat. |
| **Supervision** | none needed for an adult learner. |
| **Location** | indoors on a stable bench. |
| **Regulations** | none apply. |
| **Stop criteria** | anyone inside the circle; any crack in the fixture; battery warm to the touch. |
| **Disposal** | battery to a collection point. |
| **Not permitted in this project** | mains power; drops onto people or pets. |
"""

MISSING_STOP = COMPLETE.replace(
    "| **Stop criteria** | anyone inside the circle; any crack in the fixture; battery warm to the touch. |\n", ""
)

PLACEHOLDER_DISPOSAL = COMPLETE.replace("battery to a collection point.", "…")

TEMPLATE_LIKE = """\
| | |
|---|---|
| **Hazards** | |
| **Energy limits** | TBD |
| **Controls** | ... |
"""

TOO_HOT = COMPLETE.replace("mass 0.5 kg, drop height 1.2 m, 9 V battery, 0.3 A, 3 m/s.",
                           "mass 3 kg, drop height 2.5 m, 230 V mains, 500 mA, 12 m/s, 45 °C.")


def test_parse_reads_every_bold_row():
    box = impl.parse_safety_box(COMPLETE)
    assert set(impl.REQUIRED_ROWS) <= set(box)
    assert box["Supervision"].startswith("none needed")
    # separator and blank header rows must not appear as labels
    assert all(set(k) - set("-: ") for k in box)


def test_parse_accepts_unbolded_labels_and_extra_whitespace():
    box = impl.parse_safety_box("|  Hazards   |   eye   |\n| Controls | glasses |")
    assert box == {"Hazards": "eye", "Controls": "glasses"}


def test_complete_box_passes():
    report = impl.check_safety_box(COMPLETE)
    assert report.missing == []
    assert report.exceeded == []
    assert report.ok
    assert "PASS" in report.format()


def test_missing_row_is_reported():
    report = impl.check_safety_box(MISSING_STOP)
    assert report.missing == ["Stop criteria"]
    assert not report.ok
    assert "Stop criteria" in report.format()


def test_placeholder_counts_as_missing():
    assert impl.missing_rows(impl.parse_safety_box(PLACEHOLDER_DISPOSAL)) == ["Disposal"]


def test_template_like_box_reports_all_rows_in_order():
    missing = impl.missing_rows(impl.parse_safety_box(TEMPLATE_LIKE))
    assert missing == list(impl.REQUIRED_ROWS)


def test_row_matching_is_case_insensitive():
    box = {"hazards": "x", "ENERGY LIMITS": "0.1 kg"}
    missing = impl.missing_rows(box, required=("Hazards", "Energy limits", "Controls"))
    assert missing == ["Controls"]


@pytest.mark.parametrize(
    "text, key, expected",
    [
        ("mass 0.5 kg", "mass_kg", 0.5),
        ("250 g steel ball", "mass_kg", 0.25),
        ("drop height 1.2 m", "drop_height_m", 1.2),
        ("fall of 80 cm", "drop_height_m", 0.8),
        ("a 9 V battery", "voltage_v", 9.0),
        ("500 mA fuse", "current_a", 0.5),
        ("up to 3 m/s", "speed_m_s", 3.0),
        ("36 km/h", "speed_m_s", 10.0),
        ("below 45 °C", "temperature_c", 45.0),
        ("1,5 kg (decimal comma)", "mass_kg", 1.5),
    ],
)
def test_energy_limit_units_convert_to_si_keys(text, key, expected):
    limits = impl.parse_energy_limits(text)
    assert limits[key] == pytest.approx(expected, rel=1e-6)


def test_speed_is_not_misread_as_height():
    limits = impl.parse_energy_limits("3 m/s maximum")
    assert "drop_height_m" not in limits
    assert limits["speed_m_s"] == pytest.approx(3.0)


def test_milliamps_are_not_misread_as_metres():
    limits = impl.parse_energy_limits("500 mA")
    assert "drop_height_m" not in limits and "mass_kg" not in limits


def test_maximum_is_kept_when_quantity_repeats():
    limits = impl.parse_energy_limits("0.2 kg ball or 0.9 kg bag")
    assert limits["mass_kg"] == pytest.approx(0.9)


def test_no_numbers_gives_empty_limits():
    assert impl.parse_energy_limits("no stored energy; hand power only") == {}


def test_values_at_the_envelope_pass_and_above_fail():
    env = impl.HOME_LAB_ENVELOPE
    at_limit = {k: v for k, v in env.items() if k not in ("energy_j", "pressure_kpa")}
    assert impl.check_envelope(at_limit) == []
    above = {k: v * 1.01 for k, v in at_limit.items()}
    exceeded = impl.check_envelope(above)
    assert set(above) <= set(exceeded)
    assert set(exceeded) <= set(above) | {"energy_j"}  # 2.02 kg at 2.02 m is just over 40 J


def test_stored_energy_is_derived_from_mass_and_height():
    # 2.0 kg at 2.0 m is 39.2 J: inside the 40 J limit. 2.0 kg at 2.1 m exceeds height and energy.
    assert impl.check_envelope({"mass_kg": 2.0, "drop_height_m": 2.0}) == []
    exceeded = impl.check_envelope({"mass_kg": 2.0, "drop_height_m": 2.1})
    assert set(exceeded) == {"drop_height_m", "energy_j"}
    # a light object from a legal height but heavy mass: only energy trips
    assert impl.check_envelope({"mass_kg": 1.9, "drop_height_m": 1.9}) == []
    assert "energy_j" not in impl.check_envelope({"mass_kg": 1.0, "drop_height_m": 1.0})


def test_over_envelope_box_fails_with_each_offender_named():
    report = impl.check_safety_box(TOO_HOT)
    assert report.missing == []
    assert {"mass_kg", "drop_height_m", "voltage_v", "speed_m_s", "energy_j"} <= set(report.exceeded)
    assert "current_a" not in report.exceeded and "temperature_c" not in report.exceeded
    assert not report.ok
    assert "FAIL" in report.format()


def test_custom_envelope_is_respected():
    strict = dict(impl.HOME_LAB_ENVELOPE, voltage_v=5.0)
    assert impl.check_envelope({"voltage_v": 9.0}, envelope=strict) == ["voltage_v"]
    assert impl.check_envelope({"voltage_v": 9.0}) == []


def test_builtin_example_passes():
    assert impl.check_safety_box(impl.EXAMPLE_SAFETY_BOX).ok
