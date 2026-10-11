"""P00.11 Reading a Technical Paper: starter code.

Tools for reading NACA Report 824 ("Summary of Airfoil Data"): a parser for NACA four-digit
designations, the four-digit camber line and thickness distribution so you can draw the section
the designation describes, helpers for reading a lift curve and a drag polar, and a validator for
the reading-notes template used in this lesson.

Complete every function marked NotImplementedError. When finished, running this file prints the
parsed designations and saves ``naca_airfoils.png``.

Run:  python main.py
Test: pytest curriculum/phase-00/lesson-11-reading-a-technical-paper
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Headings every reading-notes file must contain, in order of appearance in the template.
REQUIRED_NOTE_FIELDS: tuple[str, ...] = (
    "Citation",
    "Question",
    "Method",
    "Key figure",
    "Result",
    "Limitations",
    "My questions",
)
# Body text that counts as "not written yet".
_PLACEHOLDERS = {"", "...", "…", "tbd", "todo", "n/a", "-"}


@dataclass(frozen=True)
class Naca4:
    """A NACA four-digit section, all lengths as fractions of chord.

    ``max_camber``: maximum camber ``m`` (first digit / 100).
    ``camber_position``: chordwise position of maximum camber ``p`` (second digit / 10).
    ``thickness``: maximum thickness ``t`` (last two digits / 100).
    """

    designation: str
    max_camber: float
    camber_position: float
    thickness: float

    @property
    def is_symmetric(self) -> bool:
        return self.max_camber == 0.0


def parse_naca4(designation: str) -> Naca4:
    """Parse "2412", "NACA 2412" or "naca2412" into its camber, camber position and thickness.

    Raises ``ValueError`` for anything that is not four digits, for zero thickness, and for a
    cambered section (first digit > 0) whose camber position digit is 0.
    """
    raise NotImplementedError


def camber_line(x: np.ndarray, m: float, p: float) -> np.ndarray:
    """Mean camber line ``y_c(x)`` of a NACA four-digit section, ``x`` and ``y_c`` in chord fractions.

    Two parabolic arcs meeting at ``x = p`` with value ``m`` and zero slope. A symmetric
    section (``m == 0``) returns zeros.
    """
    raise NotImplementedError


def camber_slope(x: np.ndarray, m: float, p: float) -> np.ndarray:
    """Slope ``dy_c/dx`` of the four-digit camber line."""
    raise NotImplementedError


def thickness_distribution(x: np.ndarray, t: float) -> np.ndarray:
    """Half-thickness ``y_t(x)`` of a NACA four-digit section (chord fractions).

    The NACA polynomial: ``y_t = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - 0.1015 x^4)``.
    Its maximum is ``t/2`` at about 30 % chord; the trailing edge has a small finite thickness.
    """
    raise NotImplementedError


def airfoil_surfaces(designation: str, n: int = 101) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Upper and lower surface coordinates ``(xu, yu, xl, yl)`` of a four-digit section.

    Points are cosine-spaced in ``x`` to resolve the nose. The thickness is laid off
    perpendicular to the camber line, as in the NACA construction.
    """
    raise NotImplementedError


def lift_curve_fit(alpha_deg: np.ndarray, cl: np.ndarray, alpha_range: tuple[float, float] = (-6.0, 8.0)) -> tuple[float, float]:
    """Fit the linear part of a lift curve; return ``(slope_per_degree, zero_lift_angle_deg)``.

    Only points with ``alpha_range[0] <= alpha <= alpha_range[1]`` are used (the pre-stall,
    linear region). Raises ``ValueError`` if fewer than two points fall in the range.
    """
    raise NotImplementedError


def drag_polar_summary(cl: np.ndarray, cd: np.ndarray) -> dict[str, float]:
    """Read the two numbers an engineer takes from a drag polar.

    Returns ``cd_min`` and the ``cl_at_cd_min`` where it occurs, and ``ld_max`` (maximum
    lift-to-drag ratio) with the ``cl_at_ld_max`` where it occurs. Points with ``cl <= 0`` are
    excluded from the L/D search.
    """
    raise NotImplementedError


def validate_reading_notes(text: str) -> list[str]:
    """Check a reading-notes Markdown file against ``REQUIRED_NOTE_FIELDS``.

    Each field must appear as a heading (any level, case-insensitive, substring match) and have a
    body that is not empty or a placeholder ("...", "TODO", "tbd"). Returns one message per
    problem; an empty list means the notes are complete.
    """
    raise NotImplementedError


EXAMPLE_NOTES = """\
# Reading notes: NACA Report 824

## Citation
Abbott, von Doenhoff, Stivers (1945), NACA Report 824, "Summary of Airfoil Data", NTRS 19930090976.

## Question
What are the measured lift and drag characteristics of the NACA airfoil families, and how do the
designations encode geometry?

## Method
Wind-tunnel tests in the Langley two-dimensional low-turbulence pressure tunnel at Reynolds numbers
of 3, 6 and 9 million; smooth and standard-roughness conditions.

## Key figure
The section characteristics plot for NACA 2412: left panel Cl against alpha (deg), right panel
Cd against Cl, three Reynolds numbers plus a roughness case.

## Result
Lift-curve slope near 0.1 per degree for all sections in the linear range; minimum drag rises with
thickness; roughness roughly doubles minimum drag.

## Limitations
Two-dimensional data at subsonic speeds; no sweep or finite-span effects; roughness is a standard
grain, not a real surface.

## My questions
Why does the drag bucket exist for the 6-series and not the 4-digit sections?
"""


def main(out_dir: str | Path = ".") -> None:
    """Print the three worked designations, save ``naca_airfoils.png``, and validate the example notes."""
    for d in ("2412", "0012", "4415"):
        s = parse_naca4(d)
        print(f"NACA {s.designation}: m = {s.max_camber:.2f}, p = {s.camber_position:.1f}, t = {s.thickness:.2f}"
              + ("  (symmetric)" if s.is_symmetric else ""))
    fig, ax = plt.subplots(figsize=(8, 3))
    for d in ("2412", "0012", "4415"):
        xu, yu, xl, yl = airfoil_surfaces(d, 201)
        ax.plot(np.r_[xu[::-1], xl], np.r_[yu[::-1], yl], label=f"NACA {d}")
    ax.set_aspect("equal")
    ax.set_xlabel("x / c")
    ax.set_ylabel("y / c")
    ax.grid(True, alpha=0.3)
    ax.legend()
    out = Path(out_dir) / "naca_airfoils.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {out}")
    problems = validate_reading_notes(EXAMPLE_NOTES)
    print("example notes:", "complete" if not problems else problems)


if __name__ == "__main__":
    main()
