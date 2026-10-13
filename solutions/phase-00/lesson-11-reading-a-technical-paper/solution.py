"""P00.11 Reading a Technical Paper: reference solution.

Tools for reading NACA Report 824 ("Summary of Airfoil Data"): a parser for NACA four-digit
designations, the four-digit camber line and thickness distribution so you can draw the section
the designation describes, helpers for reading a lift curve and a drag polar, and a validator for
the reading-notes template used in this lesson.

Run:  python solution.py            (prints parsed designations, saves naca_airfoils.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-00/lesson-11-reading-a-technical-paper
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
_NACA4_RE = re.compile(r"^\s*(?:naca)?\s*(\d)(\d)(\d{2})\s*$", re.IGNORECASE)


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
    m = _NACA4_RE.match(designation)
    if not m:
        raise ValueError(f"not a NACA four-digit designation: {designation!r}")
    d1, d2, d34 = int(m.group(1)), int(m.group(2)), int(m.group(3))
    if d34 == 0:
        raise ValueError(f"{designation!r}: thickness digits '00' mean a zero-thickness section")
    if d1 > 0 and d2 == 0:
        raise ValueError(f"{designation!r}: a cambered section needs a non-zero camber position")
    return Naca4(
        designation=f"{d1}{d2}{d34:02d}",
        max_camber=d1 / 100.0,
        camber_position=d2 / 10.0,
        thickness=d34 / 100.0,
    )


def camber_line(x: np.ndarray, m: float, p: float) -> np.ndarray:
    """Mean camber line ``y_c(x)`` of a NACA four-digit section, ``x`` and ``y_c`` in chord fractions.

    Two parabolic arcs meeting at ``x = p`` with value ``m`` and zero slope. A symmetric
    section (``m == 0``) returns zeros.
    """
    x = np.asarray(x, dtype=float)
    if m == 0.0:
        return np.zeros_like(x)
    front = m / p**2 * (2 * p * x - x**2)
    back = m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2)
    return np.where(x < p, front, back)


def camber_slope(x: np.ndarray, m: float, p: float) -> np.ndarray:
    """Slope ``dy_c/dx`` of the four-digit camber line."""
    x = np.asarray(x, dtype=float)
    if m == 0.0:
        return np.zeros_like(x)
    front = 2 * m / p**2 * (p - x)
    back = 2 * m / (1 - p) ** 2 * (p - x)
    return np.where(x < p, front, back)


def thickness_distribution(x: np.ndarray, t: float) -> np.ndarray:
    """Half-thickness ``y_t(x)`` of a NACA four-digit section (chord fractions).

    The NACA polynomial: ``y_t = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - 0.1015 x^4)``.
    Its maximum is ``t/2`` at about 30 % chord; the trailing edge has a small finite thickness.
    """
    x = np.asarray(x, dtype=float)
    return 5 * t * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x**2 + 0.2843 * x**3 - 0.1015 * x**4)


def airfoil_surfaces(designation: str, n: int = 101) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Upper and lower surface coordinates ``(xu, yu, xl, yl)`` of a four-digit section.

    Points are cosine-spaced in ``x`` to resolve the nose. The thickness is laid off
    perpendicular to the camber line, as in the NACA construction.
    """
    sec = parse_naca4(designation)
    beta = np.linspace(0.0, np.pi, n)
    x = 0.5 * (1 - np.cos(beta))
    yc = camber_line(x, sec.max_camber, sec.camber_position)
    yt = thickness_distribution(x, sec.thickness)
    theta = np.arctan(camber_slope(x, sec.max_camber, sec.camber_position))
    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)
    return xu, yu, xl, yl


def lift_curve_fit(alpha_deg: np.ndarray, cl: np.ndarray, alpha_range: tuple[float, float] = (-6.0, 8.0)) -> tuple[float, float]:
    """Fit the linear part of a lift curve; return ``(slope_per_degree, zero_lift_angle_deg)``.

    Only points with ``alpha_range[0] <= alpha <= alpha_range[1]`` are used (the pre-stall,
    linear region). Raises ``ValueError`` if fewer than two points fall in the range.
    """
    alpha_deg = np.asarray(alpha_deg, dtype=float)
    cl = np.asarray(cl, dtype=float)
    mask = (alpha_deg >= alpha_range[0]) & (alpha_deg <= alpha_range[1])
    if mask.sum() < 2:
        raise ValueError("need at least two points inside alpha_range")
    slope, intercept = np.polyfit(alpha_deg[mask], cl[mask], 1)
    return float(slope), float(-intercept / slope)


def drag_polar_summary(cl: np.ndarray, cd: np.ndarray) -> dict[str, float]:
    """Read the two numbers an engineer takes from a drag polar.

    Returns ``cd_min`` and the ``cl_at_cd_min`` where it occurs, and ``ld_max`` (maximum
    lift-to-drag ratio) with the ``cl_at_ld_max`` where it occurs. Points with ``cl <= 0`` are
    excluded from the L/D search.
    """
    cl = np.asarray(cl, dtype=float)
    cd = np.asarray(cd, dtype=float)
    i_min = int(np.argmin(cd))
    pos = cl > 0
    if not pos.any():
        raise ValueError("drag polar has no positive-lift points")
    ld = np.where(pos, cl / cd, -np.inf)
    i_ld = int(np.argmax(ld))
    return {
        "cd_min": float(cd[i_min]),
        "cl_at_cd_min": float(cl[i_min]),
        "ld_max": float(ld[i_ld]),
        "cl_at_ld_max": float(cl[i_ld]),
    }


def validate_reading_notes(text: str) -> list[str]:
    """Check a reading-notes Markdown file against ``REQUIRED_NOTE_FIELDS``.

    Each field must appear as a heading (any level, case-insensitive, substring match) and have a
    body that is not empty or a placeholder ("...", "TODO", "tbd"). Returns one message per
    problem; an empty list means the notes are complete.
    """
    lines = text.splitlines()
    heads: list[tuple[int, str]] = []
    for i, line in enumerate(lines):
        m = re.match(r"^#{1,6}\s+(.*?)\s*$", line)
        if m:
            heads.append((i, m.group(1)))
    problems: list[str] = []
    for field in REQUIRED_NOTE_FIELDS:
        hit = next(((i, h) for i, h in heads if field.lower() in h.lower()), None)
        if hit is None:
            problems.append(f"missing field: {field}")
            continue
        start = hit[0] + 1
        end = next((j for j, _ in heads if j > hit[0]), len(lines))
        body = " ".join(l.strip() for l in lines[start:end] if l.strip())
        if body.strip().lower().strip("*_ ") in _PLACEHOLDERS:
            problems.append(f"empty field: {field}")
    return problems


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
