"""P00.13 Safety & Standards for Home Experiments: starter code.

A completeness checker for the safety box that every physical (🔬) lesson in this course
carries. It reads the Markdown table from a ``safety.md``, reports the rows that are missing
or still hold placeholder text, extracts the numbers in the "Energy limits" row, and compares
them with the course's home-lab envelope.

Complete every function marked NotImplementedError. When finished, running this file prints
a report for the built-in example, or for the file named on the command line.

Run:  python main.py [path/to/safety.md]
Test: pytest curriculum/phase-00/lesson-13-safety-for-home-experiments
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

G0 = 9.80665
"""Standard gravitational acceleration, m/s^2, used for the stored-energy check."""

REQUIRED_ROWS: tuple[str, ...] = (
    "Hazards",
    "Energy limits",
    "Controls",
    "Supervision",
    "Location",
    "Regulations",
    "Stop criteria",
    "Disposal",
    "Not permitted in this project",
)
"""Row labels every safety box must carry, in the order of ``projects/_template/safety.md``."""

PLACEHOLDERS: frozenset[str] = frozenset({"", "…", "...", "tbd", "todo", "n/a", "na", "-", "?"})
"""Cell contents that count as 'not yet written'."""

HOME_LAB_ENVELOPE: dict[str, float] = {
    "mass_kg": 2.0,
    "drop_height_m": 2.0,
    "voltage_v": 12.0,
    "current_a": 2.0,
    "speed_m_s": 5.0,
    "temperature_c": 60.0,
    "pressure_kpa": 0.0,
    "energy_j": 40.0,
}
"""Upper limits for any home experiment in this course (inclusive).

``energy_j`` bounds the gravitational potential energy m*g0*h of a dropped or falling mass.
``pressure_kpa`` is zero gauge: pressure vessels are never permitted at home.
"""


@dataclass
class SafetyReport:
    """Result of checking one safety box."""

    missing: list[str] = field(default_factory=list)
    limits: dict[str, float] = field(default_factory=dict)
    exceeded: list[str] = field(default_factory=list)

    @property
    def ok(self) -> bool:
        """True when no row is missing and no limit exceeds the envelope."""
        return not self.missing and not self.exceeded

    def format(self) -> str:
        """Human-readable report, one finding per line."""
        raise NotImplementedError


def parse_safety_box(text: str) -> dict[str, str]:
    """Return ``{row label: cell text}`` for every two-column Markdown table row in ``text``.

    A row looks like ``| **Hazards** | falling mass, pinch |``. The bold markers are
    optional; the label is stripped of ``*`` and surrounding whitespace. Header and
    separator rows (``|---|---|``, ``| | |``) are ignored. The first occurrence of a label wins.
    """
    raise NotImplementedError


def missing_rows(box: dict[str, str], required: tuple[str, ...] = REQUIRED_ROWS) -> list[str]:
    """Return the required labels that are absent or still hold a placeholder.

    Matching is case-insensitive. A cell counts as a placeholder when its text, lower-cased
    and stripped of surrounding punctuation, is in ``PLACEHOLDERS``.
    """
    raise NotImplementedError


def parse_energy_limits(text: str) -> dict[str, float]:
    """Extract the largest value per quantity from the text of an "Energy limits" cell.

    Recognised units: kg, g (mass); m, cm, mm (drop height); V, mV, kV (voltage);
    A, mA (current); m/s, km/h (speed); °C (temperature); kPa, bar, psi (pressure); J (energy).
    Values are converted to the keys of ``HOME_LAB_ENVELOPE``. When a quantity appears
    several times, the maximum is kept, because the envelope is a ceiling.
    A decimal comma is accepted.

    Hint: match longer unit strings ("m/s", "mA", "mm") before shorter ones ("m", "A").
    """
    raise NotImplementedError


def check_envelope(limits: dict[str, float], envelope: dict[str, float] = HOME_LAB_ENVELOPE) -> list[str]:
    """Return the quantity keys whose value exceeds the envelope (a value equal to the limit passes).

    When both a mass and a drop height are given, the gravitational potential energy
    ``m * G0 * h`` is also compared with ``envelope["energy_j"]`` and reported as ``"energy_j"``.
    Unknown keys are ignored.
    """
    raise NotImplementedError


def check_safety_box(text: str, envelope: dict[str, float] = HOME_LAB_ENVELOPE) -> SafetyReport:
    """Run the complete check on the Markdown text of a safety box."""
    raise NotImplementedError


EXAMPLE_SAFETY_BOX = """\
# Safety box: Drop-test timing (example)

| | |
|---|---|
| **Hazards** | falling mass (foot, toe); tripping on the tape measure; eye (bouncing ball). |
| **Energy limits** | mass 0.20 kg, drop height 1.8 m, no electricity beyond a 5 V phone. |
| **Controls** | clear 1 m exclusion circle; safety glasses; soft landing mat; phone on a tripod. |
| **Supervision** | none needed (adult learner); children only with an adult present. |
| **Location** | indoors, hard floor with a mat, clear of furniture edges. |
| **Regulations** | none apply. |
| **Stop criteria** | anyone enters the exclusion circle; ball bounces toward a person; any damage to the floor. |
| **Disposal** | none (no consumables). |
| **Not permitted in this project** | drops from ladders or windows; metal or glass objects; drops onto people or pets. |
"""


def main() -> None:
    """Check the file named on the command line, or the built-in example."""
    if len(sys.argv) > 1:
        text = Path(sys.argv[1]).read_text(encoding="utf-8")
        print(f"checking {sys.argv[1]}")
    else:
        text = EXAMPLE_SAFETY_BOX
        print("checking the built-in example (pass a path to check your own safety.md)")
    report = check_safety_box(text)
    print(report.format())


if __name__ == "__main__":
    main()
