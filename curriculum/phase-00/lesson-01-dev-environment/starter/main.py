"""P00.01 Dev Environment for Engineers: starter code.

Complete every function marked NotImplementedError. The finished script is an environment
check. It reports the Python version and the installed version of each
package the course needs, compares each against a minimum, and runs one small numerical smoke
test (the ISA sea-level speed of sound) to prove that the numerical stack actually computes.

Run:  python main.py            (exit status 0 when every check passes, 1 otherwise)
Test: pytest curriculum/phase-00/lesson-01-dev-environment   (from the repo root)
"""
from __future__ import annotations

import importlib.metadata
import re
import sys
from collections.abc import Callable
from dataclasses import dataclass

import numpy as np

# Minimum versions the course relies on (mirrors requirements.txt at the repo root).
REQUIREMENTS: dict[str, str] = {
    "numpy": "1.26",
    "scipy": "1.11",
    "matplotlib": "3.8",
    "sympy": "1.12",
    "pint": "0.23",
    "pytest": "8.0",
}
PYTHON_MINIMUM: tuple[int, ...] = (3, 11)

# ISA sea-level constants used by the smoke test.
GAMMA = 1.4
R_AIR = 287.05287   # J/(kg K)
T0 = 288.15         # K


@dataclass(frozen=True)
class CheckResult:
    """One line of the environment report."""

    name: str
    installed: str | None   # version string, or None when the package is not importable
    minimum: str
    ok: bool


def parse_version(text: str) -> tuple[int, ...]:
    """Turn a version string into a tuple of integers for comparison.

    ``"1.26.4"`` -> ``(1, 26, 4)``; ``"2.0.0rc1"`` -> ``(2, 0, 0)`` (pre-release tags are
    dropped); ``"3.13"`` -> ``(3, 13)``. Raise ValueError when the string does not start with
    a number.
    """
    raise NotImplementedError


def meets_minimum(installed: str, minimum: str) -> bool:
    """True when ``installed`` is at least ``minimum`` under numeric, component-wise comparison.

    Shorter tuples are padded with zeros, so ``"1.26"`` and ``"1.26.0"`` compare equal, and
    ``"1.9"`` is correctly *below* ``"1.26"`` (string comparison would get that wrong).
    """
    raise NotImplementedError


def installed_version(name: str) -> str | None:
    """Return the installed distribution version of ``name``, or None if it is not installed."""
    raise NotImplementedError


def check_python(
    version_info: tuple[int, ...] = tuple(sys.version_info[:3]),
    minimum: tuple[int, ...] = PYTHON_MINIMUM,
) -> CheckResult:
    """Check the interpreter version. Both arguments default to the running interpreter."""
    raise NotImplementedError


def check_packages(
    requirements: dict[str, str] = REQUIREMENTS,
    lookup: Callable[[str], str | None] = installed_version,
) -> list[CheckResult]:
    """Check every package in ``requirements`` using ``lookup`` to find installed versions.

    ``lookup`` is injectable so the logic can be tested without touching the real environment.
    A package that is not installed gets ``installed=None`` and ``ok=False``.
    """
    raise NotImplementedError


def speed_of_sound_smoke() -> float:
    """Return the ISA sea-level speed of sound, sqrt(gamma R T0), computed with NumPy.

    A one-line proof that the numerical stack imports and computes. Expected: about 340.29 m/s.
    """
    raise NotImplementedError


def format_report(results: list[CheckResult]) -> str:
    """Render the check results as a fixed-width text table with an OK/FAIL column."""
    raise NotImplementedError


def main() -> int:
    """Print the report and the smoke test; return 0 when everything passes, else 1."""
    results = [check_python(), *check_packages()]
    print(format_report(results))
    a0 = speed_of_sound_smoke()
    smoke_ok = abs(a0 - 340.29) < 0.05
    print(f"\nsmoke test: ISA sea-level speed of sound = {a0:.2f} m/s "
          f"(expected 340.29)  {'OK' if smoke_ok else 'FAIL'}")
    all_ok = all(r.ok for r in results) and smoke_ok
    print("\nenvironment ready" if all_ok else "\nenvironment NOT ready: fix the FAIL lines above")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
