"""P05.09 Shear & Bending Moment Diagrams: starter code.

Cut a beam at x and keep the part to the left: the internal shear force V(x) and bending moment
M(x) at the cut are whatever keeps that part in equilibrium. Written with singularity (Macaulay)
functions, every point force, couple and linearly varying distributed load contributes one closed
form term, and one formula gives the whole diagram. This module builds:

* load objects (point force, point couple, linearly varying distributed load, and a general
  distributed load given as a function) and the reactions of statically determinate beams
  (simply supported with or without overhangs, cantilever), with a check for indeterminacy;
* V(x) and M(x) by singularity functions, the jumps at point loads and couples, and the extreme
  moments from the zeros of V;
* an independent check: the method of sections, summing the forces and moments on the left free
  body directly, with distributed loads integrated by scipy.integrate.quad;
* the lesson's aerospace cases: a wing spar under an elliptical lift distribution with inertia
  relief from structure, fuel and an engine, and a fuselage floor beam.

Sign convention: x along the beam from the left end; forces positive UP; distributed load w
positive UP [N/m]; applied couples positive COUNTERCLOCKWISE. Internal shear V(x) is the sum of
the upward forces to the left of the cut; internal moment M(x) is positive sagging (compression
on top). Then dV/dx = w and dM/dx = V. Units: N, m, N m.

Run:  python main.py [out_dir]     (prints the lesson's tables, writes vm_diagrams.png)
Test: pytest curriculum/phase-05/lesson-09-shear-bending-moment-diagrams
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
from scipy.integrate import quad

G0 = 9.80665     # m/s^2


# ---------------------------------------------------------------- loads

@dataclass(frozen=True)
class PointForce:
    """Concentrated force P [N] at x = a, positive up."""
    a: float
    P: float


@dataclass(frozen=True)
class PointMoment:
    """Concentrated couple M0 [N m] at x = a, positive counterclockwise."""
    a: float
    M0: float


@dataclass(frozen=True)
class DistributedLoad:
    """Load intensity varying linearly from w_a at x = a to w_b at x = b [N/m], positive up."""
    a: float
    b: float
    w_a: float
    w_b: float


@dataclass(frozen=True)
class DistributedFunction:
    """Load intensity w(x) [N/m] given as a function on a <= x <= b (positive up). Handled by the
    method of sections (quadrature); use :func:`piecewise_linear_load` to feed it to the
    singularity-function solution."""
    a: float
    b: float
    w: Callable[[float], float]


def macaulay(x: np.ndarray | float, a: float, n: int, left: bool = False) -> np.ndarray:
    """Singularity function <x - a>^n. For n >= 0 it is (x - a)^n where x >= a and zero before
    (so <x - a>^0 is the unit step, taken as 1 at x = a; with ``left=True`` it is taken as 0 at
    x = a, which gives limits from the left). For n = -1, -2 (the Dirac delta and the doublet that
    represent a point force and a point couple in the load function) it is zero everywhere except
    at x = a, where it is not a function; this returns zeros for them."""
    raise NotImplementedError


def load_resultant(load) -> tuple[float, float]:
    """(total force, moment about x = 0 counterclockwise) of one load."""
    raise NotImplementedError


def solve_reactions(supports: list[tuple[str, float]], loads: list) -> list:
    """Reactions of a statically determinate beam with transverse loads only.

    supports: list of (kind, x) with kind 'pin' or 'roller' (one vertical force each) or 'fixed'
    (a vertical force and a couple). Exactly two unknowns are allowed: two pins/rollers, or one
    fixed end; anything else is statically indeterminate (or a mechanism) for transverse loading
    and raises ValueError. Solves sum F = 0 and sum M about x = 0 = 0. Returns the reactions as
    PointForce / PointMoment objects."""
    raise NotImplementedError


def shear_moment(x: np.ndarray | float, loads: list, left: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """V(x) and M(x) from singularity functions; ``loads`` must include the reactions.

        point force P at a:          V += P <x-a>^0,                M += P <x-a>^1
        couple M0 at a (ccw):        V += 0,                        M -= M0 <x-a>^0
        linear load w_a..w_b on [a, b], slope k = (w_b - w_a)/(b - a):
            V += w_a <x-a>^1 + k/2 <x-a>^2 - w_b <x-b>^1 - k/2 <x-b>^2
            M += w_a/2 <x-a>^2 + k/6 <x-a>^3 - w_b/2 <x-b>^2 - k/6 <x-b>^3

    Values are right-continuous (at a point load, the value just to the right); ``left=True``
    gives the limits from the left instead."""
    raise NotImplementedError


def sections_check(x: float, loads: list) -> tuple[float, float]:
    """Independent V and M at one cut x by the method of sections: sum the upward forces to the
    left of the cut, and the moments about the cut of those forces and couples, integrating every
    distributed load over [a, min(x, b)] with scipy.integrate.quad. Point loads at exactly x count
    as to the left (right-continuous, like :func:`shear_moment`)."""
    raise NotImplementedError


def piecewise_linear_load(a: float, b: float, w: Callable[[float], float], n: int) -> list[DistributedLoad]:
    """Approximate a distributed load w(x) on [a, b] by n linear panels through its values at
    n + 1 equally spaced points (so the singularity solution can use it)."""
    raise NotImplementedError


def breakpoints(loads: list, L: float) -> np.ndarray:
    """Sorted x positions where the load description changes (ends, supports, load ends)."""
    pts = {0.0, L}
    for l in loads:
        if isinstance(l, (PointForce, PointMoment)):
            pts.add(l.a)
        else:
            pts.update((l.a, l.b))
    return np.array(sorted(p for p in pts if 0.0 <= p <= L))


def extreme_moments(loads: list, L: float) -> tuple[tuple[float, float], tuple[float, float]]:
    """Largest positive and most negative bending moment on [0, L] as ((x, M_max), (x, M_min)).
    Candidates: both one-sided values at every breakpoint (M jumps at a couple) and every zero of
    V inside each interval, where V is a polynomial of degree <= 2 (recovered exactly from three
    interior samples)."""
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's aerospace cases

@dataclass
class WingCase:
    """Half-wing of a light twin at a symmetric pull-up, root at the aircraft centre-line (x = 0),
    tip at x = semispan. Lift elliptical; structure, fuel and engine weigh n times their weight."""
    mass: float = 1900.0          # kg, aircraft
    n: float = 3.8                # limit load factor
    semispan: float = 6.7         # m
    wing_mass: float = 80.0       # kg per side, uniform over the semispan
    fuel_mass: float = 75.0       # kg per side, uniform on [fuel_a, fuel_b]
    fuel_a: float = 0.5
    fuel_b: float = 3.0
    engine_mass: float = 150.0    # kg per side (engine + nacelle), at engine_x
    engine_x: float = 1.9

    @property
    def w0(self) -> float:
        """Peak lift per unit span at the root [N/m]: half the lift n m g spread elliptically,
        integral of w0 sqrt(1 - (x/s)^2) over [0, s] = pi w0 s / 4."""
        return 0.5 * self.n * self.mass * G0 / (np.pi * self.semispan / 4.0)

    def lift(self, x: float) -> float:
        return self.w0 * np.sqrt(max(0.0, 1.0 - (x / self.semispan) ** 2))

    def inertia_loads(self) -> list:
        ng = self.n * G0
        s = self.semispan
        return [DistributedLoad(0.0, s, -self.wing_mass * ng / s, -self.wing_mass * ng / s),
                DistributedLoad(self.fuel_a, self.fuel_b, -self.fuel_mass * ng / (self.fuel_b - self.fuel_a),
                                -self.fuel_mass * ng / (self.fuel_b - self.fuel_a)),
                PointForce(self.engine_x, -self.engine_mass * ng)]

    def loads(self, panels: int = 400, relief: bool = True) -> list:
        """All applied loads for the singularity solution (lift as linear panels)."""
        lift = piecewise_linear_load(0.0, self.semispan, self.lift, panels)
        return lift + (self.inertia_loads() if relief else [])

    def exact_lift_root(self) -> tuple[float, float]:
        """Root shear and moment of the elliptical lift alone in closed form: the lift on the half
        wing, pi w0 s / 4, and its moment about the root, integral w0 sqrt(1-(x/s)^2) x dx = w0 s^2 / 3.
        As root reactions of the cantilever these are (-pi w0 s / 4, ...); this returns the
        magnitudes (V, M) the root carries."""
        return np.pi * self.w0 * self.semispan / 4.0, self.w0 * self.semispan ** 2 / 3.0


FLOOR_SPAN = 3.6
FLOOR_LOADS = [PointForce(0.45, -2.0e3), PointForce(1.25, -2.0e3), PointForce(2.35, -2.0e3),
               PointForce(3.15, -2.0e3), DistributedLoad(0.0, FLOOR_SPAN, -500.0, -500.0),
               PointMoment(2.0, 600.0)]


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save vm_diagrams.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)

    # ---- textbook checks
    L, w, P = 2.0, -1000.0, -3000.0
    ss = [DistributedLoad(0, L, w, w)]
    allss = ss + solve_reactions([("pin", 0.0), ("roller", L)], ss)
    (xm, Mm), _ = extreme_moments(allss, L)
    print(f"simply supported, UDL {-w:.0f} N/m over {L} m: M_max = {Mm:.1f} N m at x = {xm:.3f} m "
          f"(w L^2 / 8 = {-w * L ** 2 / 8:.1f})")
    pl = [PointForce(0.5, P)]
    allpl = pl + solve_reactions([("pin", 0.0), ("roller", L)], pl)
    (xm, Mm), _ = extreme_moments(allpl, L)
    print(f"simply supported, {-P:.0f} N at 0.5 m: M_max = {Mm:.1f} N m at x = {xm:.3f} m (P a b / L = {-P * 0.5 * 1.5 / L:.1f})")
    cl = [PointForce(L, P)]
    allcl = cl + solve_reactions([("fixed", 0.0)], cl)
    print(f"cantilever, {-P:.0f} N at the tip: root V = {shear_moment(0.0, allcl)[0]:.1f} N, "
          f"root M = {shear_moment(0.0, allcl)[1]:.1f} N m (-P L = {P * L:.1f})")

    # ---- fuselage floor beam
    R = solve_reactions([("pin", 0.0), ("roller", FLOOR_SPAN)], FLOOR_LOADS)
    allf = FLOOR_LOADS + R
    (xmax, Mmax), (xmin, Mmin) = extreme_moments(allf, FLOOR_SPAN)
    print(f"\nfloor beam, span {FLOOR_SPAN} m, four 2.0 kN seat-track loads, 0.5 kN/m floor, 0.6 kN m couple at 2.0 m:")
    print(f"  reactions {R[0].P:.1f} N and {R[1].P:.1f} N; M_max = {Mmax:.1f} N m at x = {xmax:.3f} m; "
          f"M just left/right of the couple: {float(shear_moment(2.0, allf, left=True)[1]):.1f} / {float(shear_moment(2.0, allf)[1]):.1f} N m")
    xs = np.array([0.3, 1.0, 2.0, 2.9])
    diff = max(abs(np.array(sections_check(x, allf)) - np.array([v[0] for v in shear_moment(np.array([x]), allf)])).max()
               for x in xs)
    print(f"  method of sections against singularity functions at 4 cuts: max |diff| = {diff:.1e}")

    # ---- the wing spar
    wc = WingCase()
    s = wc.semispan
    print(f"\nwing spar: {wc.mass:.0f} kg twin, n = {wc.n}, semispan {s} m, root lift intensity w0 = {wc.w0:.0f} N/m")
    V_ex, M_ex = wc.exact_lift_root()
    print(f"  elliptical lift alone, exact: root shear {V_ex:.0f} N, root moment {M_ex:.0f} N m, "
          f"centre of lift at {M_ex / V_ex:.4f} m = 4 s / 3 pi = {4 * s / (3 * np.pi):.4f} m")
    print("  panels   root V (N)   root M (N m)   rel. error in M (lift only)")
    errs = []
    for npan in (10, 20, 40, 80, 160, 320):
        lo = wc.loads(npan, relief=False)
        allw = lo + solve_reactions([("fixed", 0.0)], lo)
        V0, M0 = shear_moment(0.0, allw)
        errs.append(abs(float(M0) - M_ex) / M_ex)
        print(f"  {npan:6d}   {float(V0):9.1f}   {float(M0):11.1f}     {errs[-1]:.2e}")
    print(f"  observed orders {np.array2string(np.log2(np.array(errs[:-1]) / errs[1:]), precision=2)}")
    lo = wc.loads()
    allw = lo + solve_reactions([("fixed", 0.0)], lo)
    V0, M0 = (float(v) for v in shear_moment(0.0, allw))
    print(f"  with inertia relief (400 panels): root V = {V0:.0f} N, root M = {M0:.0f} N m "
          f"({M0 / M_ex * 100:.1f} % of lift alone)")
    for name, part in (("wing structure", wc.inertia_loads()[0:1]), ("fuel", wc.inertia_loads()[1:2]),
                       ("engine", wc.inertia_loads()[2:3])):
        F, Mo = np.sum([load_resultant(p) for p in part], axis=0)
        print(f"    relief from {name:15s}: {F:8.0f} N, root moment {Mo:8.0f} N m")
    Vl_, _ = shear_moment(wc.engine_x, allw, left=True)
    Vr_, _ = shear_moment(wc.engine_x, allw)
    print(f"  shear jump at the engine: {float(Vr_ - Vl_):.0f} N (engine n m g = {wc.n * wc.engine_mass * G0:.0f} N)")
    # sections check with the exact elliptic function
    exact_loads = [DistributedFunction(0.0, s, wc.lift)] + wc.inertia_loads()
    allx = exact_loads + solve_reactions([("fixed", 0.0)], exact_loads)
    print("  x (m)   V singularity   V sections   M singularity   M sections")
    for x in (0.0, 1.0, 1.9, 3.0, 5.0):
        Vs, Ms = shear_moment(np.array([x]), allw)
        Vc, Mc = sections_check(x, allx)
        print(f"  {x:4.1f}   {Vs[0]:12.1f}   {Vc:10.1f}   {Ms[0]:12.1f}   {Mc:10.1f}")

    # ---- experiment: engine position and fuel placement
    print("\nexperiment: root bending moment (N m) against engine position and fuel tank span")
    print("  engine x (m)   fuel 0.5-3.0 m   fuel 3.0-5.5 m   no fuel")
    for ex in (1.2, 1.9, 2.6, 3.3):
        row = f"  {ex:10.1f}  "
        for fa, fb, fm in ((0.5, 3.0, 75.0), (3.0, 5.5, 75.0), (0.5, 3.0, 0.0)):
            c = WingCase(engine_x=ex, fuel_a=fa, fuel_b=fb, fuel_mass=fm)
            l = c.loads()
            row += f"   {float(shear_moment(0.0, l + solve_reactions([('fixed', 0.0)], l))[1]):12.0f}"
        print(row)

    # ---- figure
    fig, ax = plt.subplots(2, 2, figsize=(12, 7), sharex="col")
    xf = np.linspace(0, FLOOR_SPAN, 2001)
    Vf, Mf = shear_moment(xf, allf)
    ax[0, 0].plot(xf, Vf / 1e3); ax[0, 0].set_ylabel("V (kN)"); ax[0, 0].set_title("fuselage floor beam")
    ax[1, 0].plot(xf, Mf / 1e3); ax[1, 0].set_ylabel("M (kN m)"); ax[1, 0].set_xlabel("x (m)")
    xw = np.linspace(0, s, 2001)
    Vw, Mw = shear_moment(xw, allw)
    lo_l = wc.loads(relief=False)
    Vl, Ml = shear_moment(xw, lo_l + solve_reactions([("fixed", 0.0)], lo_l))
    ax[0, 1].plot(xw, Vw / 1e3, label="lift + inertia relief"); ax[0, 1].plot(xw, Vl / 1e3, "--", label="lift only")
    ax[0, 1].set_title(f"wing spar at n = {wc.n}"); ax[0, 1].legend()
    ax[1, 1].plot(xw, Mw / 1e3); ax[1, 1].plot(xw, Ml / 1e3, "--"); ax[1, 1].set_xlabel("x from root (m)")
    for a in ax.ravel():
        a.axhline(0, color="k", lw=0.5)
    fig.tight_layout()
    fig.savefig(out_dir / "vm_diagrams.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'vm_diagrams.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
