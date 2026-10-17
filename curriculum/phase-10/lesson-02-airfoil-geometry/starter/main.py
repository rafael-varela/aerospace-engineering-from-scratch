"""P10.02 Airfoil Geometry: starter code.

Every 2-D aerodynamic code in this phase starts from a list of (x, y) points around an airfoil.
This module generates the NACA 4- and 5-digit families from their defining equations and checks
them the way you would check anyone else's coordinates:

* cosine (half-cosine) spacing, which clusters points at the leading and trailing edges;
* the NACA 4-digit thickness form, open (standard) or closed trailing edge;
* the 4-digit two-parabola camber line, and the 5-digit "230-type" cubic-then-linear camber line
  whose constants r and k1 are *derived* here from the two conditions that define them (the camber
  is greatest at x = p; the ideal lift coefficient is 0.15 L for designation L p 0), not copied;
* thickness laid off perpendicular to the camber line, points returned in Selig order
  (trailing edge -> upper surface -> leading edge -> lower surface -> trailing edge);
* checks: maximum thickness and its position, maximum camber and its position, the cross-section
  area (shoelace formula; 0.685 t c for the thickness form, from exact integration), the
  leading-edge radius 1.1019 t^2 c, and the trailing-edge thickness.

Coordinates are in units of the chord unless a chord is passed. x from the leading edge, y up.

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes airfoil_geometry.png)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import optimize

# NACA 4-digit thickness-form coefficients (open trailing edge); a4 = -0.1036 closes it
A_THICK = np.array([0.2969, -0.1260, -0.3516, 0.2843, -0.1015])
A4_CLOSED = -0.1036

# the (r, k1) pairs as commonly quoted for the 210-250 camber lines, checked (not used) in main()
QUOTED_5DIGIT = {0.05: (0.0580, 361.4), 0.10: (0.1260, 51.64), 0.15: (0.2025, 15.957), 0.20: (0.2900, 6.643),
                 0.25: (0.3910, 3.230)}

TRAINER_CHORD = 1.2     # m, the reference trainer's NACA 2412 wing
TRAINER_SPAN = 10.0     # m


# ---------------------------------------------------------------- spacing

def cosine_spacing(n: int) -> np.ndarray:
    """n + 1 chordwise stations x_i = (1 - cos(beta_i))/2, beta_i = i pi / n, from 0 to 1."""
    raise NotImplementedError


# ---------------------------------------------------------------- thickness and camber

def naca4_thickness(x: np.ndarray, t: float, closed_te: bool = False) -> np.ndarray:
    """Half-thickness y_t(x) of the NACA 4-digit thickness form (chord units):

        y_t = 5 t (0.2969 sqrt(x) - 0.1260 x - 0.3516 x^2 + 0.2843 x^3 - 0.1015 x^4),

    with -0.1036 in place of -0.1015 when closed_te (zero trailing-edge thickness)."""
    raise NotImplementedError


def naca4_camber(x: np.ndarray, m: float, p: float) -> tuple[np.ndarray, np.ndarray]:
    """NACA 4-digit camber line and its slope (chord units): two parabolas meeting at x = p with
    zero slope, height m there:

        y_c = m/p^2 (2 p x - x^2)              for x <= p,
        y_c = m/(1-p)^2 ((1 - 2p) + 2 p x - x^2)  for x > p.
    Returns (y_c, dy_c/dx). m = 0 gives a symmetric section."""
    raise NotImplementedError


def naca5_camber(x: np.ndarray, r: float, k1: float) -> tuple[np.ndarray, np.ndarray]:
    """NACA 5-digit (standard, non-reflexed) camber line and slope: a cubic up to x = r, then a
    straight line to the trailing edge:

        y_c = k1/6 (x^3 - 3 r x^2 + r^2 (3 - r) x)   for x <= r,
        y_c = k1 r^3 / 6 (1 - x)                      for x > r."""
    raise NotImplementedError


def naca5_r(p: float) -> float:
    """The constant r that puts the 5-digit camber line's maximum at x = p: the cubic's slope
    k1/6 (3x^2 - 6 r x + r^2 (3 - r)) vanishes at x = p, i.e. r^3 - 3 r^2 + 6 p r - 3 p^2 = 0
    (the root between p and 1)."""
    raise NotImplementedError


def naca5_ideal_cl(r: float, k1: float, n: int = 200) -> float:
    """Thin-airfoil ideal lift coefficient pi A1 of the 5-digit camber line (r, k1), with
    A1 = (2/pi) int_0^pi dy_c/dx cos(theta) d(theta) and x = (1 - cos theta)/2 (P10.04 derives this).
    Gauss-Legendre quadrature with n points on each side of theta_r = arccos(1 - 2r), where the
    camber slope has a kink, so each piece is smooth and the quadrature converges fast."""
    raise NotImplementedError


def naca5_constants(p: float, cl_design: float = 0.3) -> tuple[float, float]:
    """(r, k1) of the standard 5-digit camber line with maximum camber at p and thin-airfoil ideal
    lift coefficient cl_design (0.15 times the first digit; 0.3 for the 230-type series).
    r from ``naca5_r``; y_c is linear in k1, so k1 = cl_design / naca5_ideal_cl(r, 1)."""
    raise NotImplementedError


# ---------------------------------------------------------------- coordinates

def parse_naca(code: str) -> dict:
    """Decode a 4-digit ('2412') or standard 5-digit ('23012') designation into its parameters:
    4-digit: m = d1/100, p = d2/10, t = d34/100.
    5-digit (L P Q TT, Q = 0 only): cl_design = 0.15 L, p = P/20, t = TT/100."""
    code = code.strip()
    if len(code) == 4:
        return dict(kind=4, m=int(code[0]) / 100, p=int(code[1]) / 10, t=int(code[2:]) / 100)
    if len(code) == 5:
        if code[2] != "0":
            raise ValueError("only standard (non-reflexed, third digit 0) 5-digit sections are supported")
        return dict(kind=5, cl_design=0.15 * int(code[0]), p=int(code[1]) / 20, t=int(code[3:]) / 100)
    raise ValueError(f"not a NACA 4- or 5-digit code: {code}")


def naca_coordinates(code: str, n: int = 100, closed_te: bool = False,
                     chord: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    """Airfoil surface in Selig order, 2n + 1 points: trailing edge -> upper -> leading edge ->
    lower -> trailing edge, with n + 1 cosine-spaced stations per surface (the leading edge shared).

    Thickness is laid off perpendicular to the camber line:
        x_u = x - y_t sin(th), y_u = y_c + y_t cos(th),  x_l = x + y_t sin(th), y_l = y_c - y_t cos(th),
    th = arctan(dy_c/dx). Scaled by chord."""
    raise NotImplementedError


# ---------------------------------------------------------------- checks

def polygon_area(x: np.ndarray, y: np.ndarray) -> float:
    """Signed area of the closed polygon through the points (shoelace formula); positive for
    counter-clockwise order, which Selig order is. The closing edge last -> first is included."""
    raise NotImplementedError


def thickness_area_coefficient(closed_te: bool = False) -> float:
    """Exact area of the thickness form per (t c^2): int_0^1 2 y_t/t dx = 10 (2a0/3 + a1/2 + a2/3 +
    a3/4 + a4/5) = 0.68508 (open) or 0.68088 (closed)."""
    raise NotImplementedError


def max_thickness(t: float, closed_te: bool = False) -> tuple[float, float]:
    """(x of maximum thickness, maximum thickness 2 y_t) of the thickness form, by bounded
    minimisation of -y_t on (0.1, 0.6)."""
    raise NotImplementedError


def leading_edge_radius(t: float) -> float:
    """Leading-edge radius of the thickness form (chord units). Near x = 0, y_t ~ 5 t a0 sqrt(x), a
    parabola y^2 = 2 r x with r = (5 t a0)^2 / 2 = 1.1019 t^2."""
    raise NotImplementedError


def surface_thickness_camber(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Split Selig-ordered coordinates at the leading edge (minimum x), interpolate both surfaces
    to common stations, and return (xs, vertical thickness y_u - y_l, mean line (y_u + y_l)/2).
    This is how you check coordinates of unknown origin."""
    i_le = int(np.argmin(x))
    xu, yu = x[:i_le + 1][::-1], y[:i_le + 1][::-1]
    xl, yl = x[i_le:], y[i_le:]
    lo = max(xu.min(), xl.min())
    hi = min(xu.max(), xl.max())
    xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 801))) * (hi - lo) + lo
    Yu = np.interp(xs, xu, yu)
    Yl = np.interp(xs, xl, yl)
    return xs, Yu - Yl, 0.5 * (Yu + Yl)


def area_between(code: str, x_front: float, x_rear: float, chord: float = 1.0, n: int = 4000) -> float:
    """Cross-section area between two vertical planes x_front and x_rear (fractions of chord), for
    a spar box: int (y_u - y_l) dx on the generated surface, by the trapezoid rule on n points."""
    raise NotImplementedError


def area_convergence(code: str, ns: list[int], spacing: str = "cosine") -> np.ndarray:
    """Area error |polygon area - exact thickness-form area| for symmetric sections (camber zero)
    with n panels per surface on cosine or uniform spacing; used to measure observed order."""
    raise NotImplementedError


def observed_order(errors: np.ndarray, ratio: float = 2.0) -> np.ndarray:
    """Observed orders log(e_k/e_{k+1})/log(ratio) for grids refined by ``ratio`` (given)."""
    e = np.asarray(errors, float)
    return np.log(e[:-1] / e[1:]) / np.log(ratio)


def write_selig_dat(code: str, path: str | Path, n: int = 80, closed_te: bool = False) -> None:
    """Write coordinates as a Selig-format .dat file: a name line, then 'x y' rows (given)."""
    X, Y = naca_coordinates(code, n, closed_te)
    rows = [f"NACA {code}{' closed TE' if closed_te else ''}"] + [f"{a:.6f} {b:.6f}" for a, b in zip(X, Y)]
    Path(path).write_text("\n".join(rows) + "\n")


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save airfoil_geometry.png and naca2412.dat."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    print(f"thickness-form area coefficient: open TE {thickness_area_coefficient():.5f}, "
          f"closed TE {thickness_area_coefficient(True):.5f}  (area = coeff * t * c^2)")
    xm, tm = max_thickness(0.12)
    xmc, tmc = max_thickness(0.12, True)
    print(f"NACA 0012: max thickness {tm:.6f} c at x = {xm:.4f} c (closed TE: {tmc:.6f} at {xmc:.4f}); "
          f"t/c target 0.12 -> ratio {tm / 0.12:.5f}")
    print(f"TE thickness (open) = {2 * naca4_thickness(1.0, 0.12):.5f} c = {2 * naca4_thickness(1.0, 1.0):.4f} t; "
          f"LE radius = {leading_edge_radius(0.12):.6f} c = {leading_edge_radius(1.0):.4f} t^2")

    for p in (0.05, 0.10, 0.15, 0.20, 0.25):
        r, k1 = naca5_constants(p)
        xx = np.linspace(0, 1, 200001)
        yc, _ = naca5_camber(xx, r, k1)
        print(f"5-digit {int(round(p * 20)):d}x0 series (p = {p:.2f}): r = {r:.4f}, k1 = {k1:.3f}, "
              f"max camber {yc.max() * 100:.3f} % at x = {xx[np.argmax(yc)]:.4f}; quoted table (r, k1) = "
              f"{QUOTED_5DIGIT[p]} gives ideal cl {naca5_ideal_cl(*QUOTED_5DIGIT[p]):.4f}")

    c = TRAINER_CHORD
    X, Y = naca_coordinates("2412", n=400)
    A = polygon_area(X, Y)
    xs, th, cam = surface_thickness_camber(X, Y)
    print(f"\ntrainer NACA 2412, chord {c} m: area {A * c ** 2:.5f} m^2 ({A / 0.12:.4f} t c^2); "
          f"max vertical thickness {th.max() * c * 1000:.1f} mm at x = {xs[np.argmax(th)]:.3f} c; "
          f"max camber {cam.max() * c * 1000:.2f} mm at x = {xs[np.argmax(cam)]:.3f} c")
    Abox = area_between("2412", 0.15, 0.60, c)
    print(f"spar box 15-60 % chord: {Abox:.5f} m^2 per metre span; x {TRAINER_SPAN} m span = {Abox * TRAINER_SPAN:.4f} m^3 "
          f"= {Abox * TRAINER_SPAN * 1000:.0f} L gross (before ribs, skins, spar webs)")
    X5, Y5 = naca_coordinates("23012", n=400)
    xs5, th5, cam5 = surface_thickness_camber(X5, Y5)
    print(f"NACA 23012: max camber {cam5.max() * 100:.3f} % at x = {xs5[np.argmax(cam5)]:.3f}; area {polygon_area(X5, Y5):.5f} c^2")

    ns = [10, 20, 40, 80, 160, 320]
    e_cos = area_convergence("0012", ns, "cosine")
    e_uni = area_convergence("0012", ns, "uniform")
    print("\narea error vs panels per surface (NACA 0012):")
    print("     n     cosine       uniform")
    for n, a, b in zip(ns, e_cos, e_uni):
        print(f"  {n:4d}   {a:.3e}   {b:.3e}")
    print(f"  observed orders cosine:  {np.array2string(observed_order(e_cos), precision=2)}")
    print(f"  observed orders uniform: {np.array2string(observed_order(e_uni), precision=2)}")

    write_selig_dat("2412", out_dir / "naca2412.dat")
    fig, ax = plt.subplots(2, 1, figsize=(10, 6.5))
    for code, ls in (("0012", "-"), ("2412", "-"), ("23012", "--"), ("4415", ":")):
        Xc, Yc = naca_coordinates(code, 120)
        ax[0].plot(Xc, Yc, ls, label=f"NACA {code}")
    xc = cosine_spacing(120)
    ax[0].plot(xc, naca4_camber(xc, 0.02, 0.4)[0], "k-", lw=0.6, label="2412 camber line")
    ax[0].set_aspect("equal"); ax[0].set_xlabel("x/c"); ax[0].set_ylabel("y/c"); ax[0].legend(ncol=5, fontsize=8)
    ax[0].set_title("NACA sections (Selig order, cosine spacing)")
    h = 1.0 / np.array(ns)
    ax[1].loglog(h, e_cos, "o-", label="cosine spacing")
    ax[1].loglog(h, e_uni, "s-", label="uniform spacing")
    ax[1].loglog(h, e_cos[0] * (h / h[0]) ** 2, "k--", lw=0.8, label="slope 2")
    ax[1].loglog(h, e_uni[0] * (h / h[0]) ** 1.5, "k:", lw=0.8, label="slope 1.5")
    ax[1].set_xlabel("1 / panels per surface"); ax[1].set_ylabel("|area - 0.68508 t|")
    ax[1].set_title("polygon area convergence, NACA 0012"); ax[1].legend()
    fig.tight_layout()
    fig.savefig(out_dir / "airfoil_geometry.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'airfoil_geometry.png'} and {out_dir / 'naca2412.dat'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
