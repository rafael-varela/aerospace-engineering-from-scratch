"""P08.08 Potential Flow: starter code.

How much circulation must the light aeroplane's wing carry in cruise, and what does the simplest
lifting flow, a cylinder with circulation (and its image under the Joukowski map, a flat plate),
say about where that lift comes from and why it costs no drag in an inviscid fluid? This module
builds plane potential flows by superposing elementary solutions of Laplace's equation, written with
the complex potential W(z) = phi + i psi of P02.12, so that dW/dz = u - i v:

* uniform flow U e^{-i alpha} z, source (m/2pi) log(z - z0), doublet kappa e^{i beta}/(2pi (z - z0)),
  and vortex (i Gamma/2pi) log(z - z0) with Gamma positive CLOCKWISE (the aerodynamics convention,
  so that a lifting body in flow from left to right has Gamma > 0 and lift rho U Gamma upward);
* the Rankine oval (uniform flow + source + sink) and the cylinder with and without circulation;
* line integrals round any closed contour: circulation (the real part of the contour integral of
  dW/dz dz; counterclockwise it equals minus the clockwise Gamma) and volume flux (its imaginary part);
* the force on a body two independent ways: integrating the Bernoulli pressure round its surface,
  and Blasius's theorem F_x - i F_y = (i rho/2) contour integral of (dW/dz)^2 dz on any contour
  enclosing it. Both give lift rho U Gamma (Kutta-Joukowski) and zero drag (d'Alembert's paradox);
* the Joukowski map z = zeta + a^2/zeta, which turns the circle |zeta| = a into a flat plate of
  chord 4a, with the Kutta condition fixing Gamma = 4 pi U a sin(alpha).

Contours are arrays of points z_k, counterclockwise, uniformly spaced in a periodic parameter;
integrals use the periodic trapezoid rule with dz/dtheta from the FFT, which is spectrally accurate
for smooth contours.

Units: SI (m, s, kg/m^3, Pa); strengths m and Gamma in m^2/s, kappa in m^3/s.

Run:  python main.py [out_dir]     (prints the lesson's numbers, writes potential_flow.png)
Test: pytest curriculum/phase-08/lesson-08-potential-flow
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

G0 = 9.80665
RHO_3000 = 0.90913          # kg/m^3, ISA density at 3000 m (P08.01)


@dataclass(frozen=True)
class Element:
    """One elementary plane potential flow: kind is 'uniform', 'source', 'doublet' or 'vortex';
    strength is U, m, kappa or Gamma; z0 its position; angle the flow direction (uniform) or the
    doublet axis (doublet), in radians."""
    kind: str
    strength: float
    z0: complex = 0j
    angle: float = 0.0


def uniform(U: float, alpha: float = 0.0) -> Element:
    """Uniform flow of speed U at angle alpha to the x axis: W = U e^{-i alpha} z."""
    return Element("uniform", U, 0j, alpha)


def source(m: float, z0: complex = 0j) -> Element:
    """Source of volume flow m per unit depth (a sink if m < 0): W = (m / 2 pi) log(z - z0)."""
    return Element("source", m, complex(z0))


def doublet(kappa: float, z0: complex = 0j, beta: float = 0.0) -> Element:
    """Doublet of strength kappa with axis angle beta: W = kappa e^{i beta} / (2 pi (z - z0)).
    It is the limit of a source at z0 - eps e^{i beta} and a sink at z0 + eps e^{i beta}, m 2 eps = kappa."""
    return Element("doublet", kappa, complex(z0), beta)


def vortex(Gamma: float, z0: complex = 0j) -> Element:
    """Point vortex of circulation Gamma, positive clockwise: W = (i Gamma / 2 pi) log(z - z0)."""
    return Element("vortex", Gamma, complex(z0))


def complex_potential(elements: list[Element], z: np.ndarray) -> np.ndarray:
    """W(z) = phi + i psi of the superposition (principal branch of the logarithms)."""
    raise NotImplementedError


def complex_velocity(elements: list[Element], z: np.ndarray) -> np.ndarray:
    """dW/dz = u - i v of the superposition, differentiated analytically."""
    raise NotImplementedError


def velocity(elements: list[Element], z: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Velocity components (u, v) at the points z."""
    raise NotImplementedError


def pressure_coefficient(elements: list[Element], z: np.ndarray, U: float) -> np.ndarray:
    """Bernoulli pressure coefficient Cp = (p - p_inf) / (rho U^2 / 2) = 1 - |V|^2 / U^2."""
    raise NotImplementedError


# ---------------------------------------------------------------- bodies

def cylinder(U: float, a: float, Gamma: float = 0.0, alpha: float = 0.0) -> list[Element]:
    """Flow past a circular cylinder of radius a at the origin, free stream U at angle alpha, with
    clockwise circulation Gamma: uniform flow + doublet kappa = 2 pi U a^2 (axis alpha) + vortex."""
    raise NotImplementedError


def cylinder_surface_cp(theta: np.ndarray, U: float, a: float, Gamma: float = 0.0) -> np.ndarray:
    """Closed-form surface pressure coefficient of the lifting cylinder (alpha = 0): the surface speed
    is u_theta = -2 U sin(theta) - Gamma/(2 pi a), so Cp = 1 - (2 sin(theta) + Gamma/(2 pi a U))^2."""
    raise NotImplementedError


def cylinder_stagnation_points(U: float, a: float, Gamma: float) -> np.ndarray:
    """Stagnation points of the lifting cylinder (alpha = 0). For Gamma <= 4 pi U a they lie on the
    surface at sin(theta) = -Gamma / (4 pi U a); above that, a single point lies in the flow below the
    cylinder on the negative y axis, at r = (Gamma + sqrt(Gamma^2 - (4 pi U a)^2)) / (4 pi U)."""
    raise NotImplementedError


def rankine_oval(U: float, m: float, a: float) -> list[Element]:
    """Uniform flow U + a source m at (-a, 0) + a sink -m at (+a, 0): a closed body (zero net source)."""
    raise NotImplementedError


def rankine_oval_size(U: float, m: float, a: float) -> tuple[float, float]:
    """Half-length and half-thickness of the Rankine oval. Stagnation points at x = +-l with
    l = a sqrt(1 + m / (pi U a)); half-thickness t from psi(0, t) = 0, i.e. t = (m / (pi U)) atan(a / t),
    solved by Brent's method."""
    raise NotImplementedError


def stream_function_oval(U: float, m: float, a: float, z: np.ndarray) -> np.ndarray:
    """Stream function of the Rankine oval with the angles measured so that psi = 0 on the dividing
    streamline: psi = U y + (m / 2 pi)(theta_1 - theta_2), theta_1, theta_2 the angles from the source
    and the sink, both in [0, 2 pi) above the axis (the principal branch would put a cut on the body)."""
    z = np.asarray(z, dtype=complex)
    th1 = np.angle(z + a)
    th2 = np.angle(z - a)
    return U * z.imag + m / (2 * np.pi) * (th1 - th2)


def oval_contour(U: float, m: float, a: float, n: int = 256) -> np.ndarray:
    """n points of the Rankine oval's surface, counterclockwise, uniformly spaced in the polar angle
    about its centre: for each angle the radius where psi = 0, found by Brent's method (the oval is
    star-shaped about its centre)."""
    raise NotImplementedError


# ---------------------------------------------------------------- contour integrals

def _dz_dtheta(contour: np.ndarray) -> np.ndarray:
    """Spectral derivative dz/dtheta of a closed contour sampled uniformly in theta on [0, 2 pi)."""
    n = contour.size
    k = np.fft.fftfreq(n, d=1.0 / n)
    if n % 2 == 0:
        k[n // 2] = 0.0
    return np.fft.ifft(1j * k * np.fft.fft(contour))


def contour_integral(values: np.ndarray, contour: np.ndarray) -> complex:
    """Contour integral of a function sampled at the contour points, sum f(z_k) z'(theta_k) dtheta
    (periodic trapezoid rule)."""
    raise NotImplementedError


def circulation_and_flux(elements: list[Element], contour: np.ndarray) -> tuple[float, float]:
    """Counterclockwise circulation (the integral of u dx + v dy) and outward volume flux (the
    integral of u dy - v dx) round the contour: the real and imaginary parts of the integral of dW/dz dz."""
    raise NotImplementedError


def force_by_pressure(elements: list[Element], contour: np.ndarray, rho: float, U: float,
                      p_inf: float = 0.0) -> tuple[float, float]:
    """Force (F_x, F_y) per unit span on the body whose surface is ``contour``, from integrating the
    Bernoulli pressure p = p_inf + rho (U^2 - |V|^2)/2 round it: F = -contour integral of p n ds,
    and with n ds = -i dz for a counterclockwise contour, F_x + i F_y = i contour integral of p dz."""
    raise NotImplementedError


def force_blasius(elements: list[Element], contour: np.ndarray, rho: float) -> tuple[float, float]:
    """Force per unit span from Blasius's theorem, F_x - i F_y = (i rho / 2) contour integral of
    (dW/dz)^2 dz, on any counterclockwise contour enclosing the body and no other singularity."""
    raise NotImplementedError


def circle_contour(radius: float, n: int = 128, center: complex = 0j) -> np.ndarray:
    """n points counterclockwise on a circle."""
    return center + radius * np.exp(2j * np.pi * np.arange(n) / n)


# ---------------------------------------------------------------- Joukowski and the Kutta condition

def kutta_circulation(U: float, a: float, alpha: float) -> float:
    """Clockwise circulation that puts the rear stagnation point of a cylinder of radius a (at angle
    of attack alpha) at zeta = a, the image of the plate's trailing edge: Gamma = 4 pi U a sin(alpha)."""
    raise NotImplementedError


def joukowski(zeta: np.ndarray, a: float) -> np.ndarray:
    """The Joukowski map z = zeta + a^2 / zeta: the circle |zeta| = a becomes the flat plate
    -2a <= x <= 2a, y = 0 (chord 4a); far away z ~ zeta, so the free stream is unchanged."""
    raise NotImplementedError


def zeta_from_z(z: np.ndarray, a: float) -> np.ndarray:
    """Inverse Joukowski map onto the exterior |zeta| >= a: zeta = (z + sqrt(z^2 - 4 a^2)) / 2, choosing
    the root with |zeta| >= a."""
    raise NotImplementedError


def plate_complex_velocity(z: np.ndarray, U: float, a: float, alpha: float, Gamma: float) -> np.ndarray:
    """dW/dz of the flow past the flat plate of chord 4a, from the cylinder flow in the zeta plane:
    dW/dz = (dW/dzeta) / (dz/dzeta) with dz/dzeta = 1 - a^2/zeta^2."""
    raise NotImplementedError


def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save potential_flow.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    mass, span, chord, V, rho = 750.0, 10.0, 1.5, 60.0, RHO_3000
    Wt = mass * G0
    Lp = Wt / span
    Gam = Lp / (rho * V)
    cl = Lp / (0.5 * rho * V ** 2 * chord)
    alpha = np.arcsin(cl / (2 * np.pi))
    print(f"light aeroplane: weight {Wt:.1f} N, span {span} m, chord {chord} m, cruise {V} m/s at 3000 m "
          f"(rho = {rho} kg/m^3)")
    print(f"  mean lift per span {Lp:.2f} N/m -> circulation Gamma = L'/(rho V) = {Gam:.3f} m^2/s, "
          f"c_l = {cl:.4f}")
    print(f"  flat plate of the same chord (Kutta): alpha = asin(c_l / 2 pi) = {np.degrees(alpha):.3f} deg")
    a = chord / 4
    print(f"  equivalent cylinder a = c/4 = {a} m: Gamma / (4 pi U a) = {Gam / (4 * np.pi * V * a):.5f}; "
          f"stagnation points at {np.degrees(np.angle(cylinder_stagnation_points(V, a, Gam))).round(3)} deg")
    els = cylinder(V, a, Gam)
    body = circle_contour(a, 128)
    Fp = force_by_pressure(els, body, rho, V)
    Fb = force_blasius(els, circle_contour(3 * a, 128), rho)
    print(f"  pressure integration on the surface: drag {Fp[0]:.2e} N/m, lift {Fp[1]:.4f} N/m")
    print(f"  Blasius on r = 3a:                    drag {Fb[0]:.2e} N/m, lift {Fb[1]:.4f} N/m; rho U Gamma = {rho * V * Gam:.4f}")
    Gk = kutta_circulation(V, a, alpha)
    big = circle_contour(5 * chord, 256)
    dWz = plate_complex_velocity(big, V, a, alpha, Gk)
    I = 0.5j * rho * contour_integral(dWz ** 2, big)
    Fx, Fy = I.real, -I.imag                               # plate (body) axes
    print(f"  Kutta circulation for that alpha {Gk:.3f} m^2/s; Blasius in the physical (plate) plane, body axes: "
          f"F_x = {Fx:.3f} N/m (leading-edge suction), F_y = {Fy:.3f} N/m")
    print(f"    wind axes: lift {-Fx * np.sin(alpha) + Fy * np.cos(alpha):.4f} N/m, "
          f"drag {Fx * np.cos(alpha) + Fy * np.sin(alpha):.1e} N/m")
    circ, flux = circulation_and_flux(els, circle_contour(2 * a))
    print(f"  circulation round r = 2a (counterclockwise) {circ:.4f} m^2/s = -Gamma; flux {flux:.1e} m^2/s")

    # Rankine oval
    Uo, mo, ao = 1.0, 2 * np.pi, 1.0
    l, t = rankine_oval_size(Uo, mo, ao)
    oval = rankine_oval(Uo, mo, ao)
    cont = oval_contour(Uo, mo, ao, 256)
    Fo = force_by_pressure(oval, cont, 1.0, Uo)
    dzt = _dz_dtheta(cont)
    u_n = np.abs((complex_velocity(oval, cont) * dzt).imag) / np.abs(dzt)       # u dy - v dx per ds
    print(f"\nRankine oval m/(pi U a) = {mo / (np.pi * Uo * ao):.1f}: half-length {l:.5f} a, half-thickness {t:.5f} a, "
          f"fineness {l / t:.4f}")
    print(f"  max |normal velocity| on the traced surface {u_n.max():.1e} U; force {Fo[0]:.1e}, {Fo[1]:.1e} "
          f"(rho U^2 a units): d'Alembert")
    print(f"  max speed on the oval {np.abs(complex_velocity(oval, cont)).max():.5f} U, min Cp "
          f"{pressure_coefficient(oval, cont, Uo).min():.5f}")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    xg, yg = np.meshgrid(np.linspace(-3, 3, 301), np.linspace(-2, 2, 201))
    zg = xg + 1j * yg
    psi = stream_function_oval(Uo, mo, ao, zg)
    ax[0].contour(xg, yg, psi, levels=np.linspace(-2, 2, 21), colors="C0", linewidths=0.7)
    ax[0].plot(cont.real, cont.imag, "k-", lw=1.5)
    ax[0].plot([-ao, ao], [0, 0], "ro", ms=4)
    ax[0].set_aspect("equal")
    ax[0].set_title("Rankine oval: streamlines")
    xg2, yg2 = np.meshgrid(np.linspace(-3 * a, 3 * a, 301), np.linspace(-2 * a, 2 * a, 201))
    zc = xg2 + 1j * yg2
    with np.errstate(divide="ignore", invalid="ignore"):
        psi_c = complex_potential(cylinder(V, a, 4 * np.pi * V * a * 0.5), zc).imag
    psi_c = np.where(np.abs(zc) < a, np.nan, psi_c)
    ax[1].contour(xg2, yg2, psi_c, levels=40, colors="C0", linewidths=0.7)
    th = np.linspace(0, 2 * np.pi, 200)
    ax[1].plot(a * np.cos(th), a * np.sin(th), "k-", lw=1.5)
    ax[1].set_aspect("equal")
    ax[1].set_title("lifting cylinder, Gamma = 2 pi U a")
    thp = np.linspace(0, 2 * np.pi, 361)
    for G_, lab in ((0.0, "Gamma = 0"), (Gam, f"Gamma = {Gam:.2f} m^2/s (cruise)"), (2 * np.pi * V * a, "Gamma = 2 pi U a")):
        ax[2].plot(np.degrees(thp), cylinder_surface_cp(thp, V, a, G_), label=lab)
    ax[2].invert_yaxis()
    ax[2].set_xlabel("theta (deg)")
    ax[2].set_ylabel("Cp")
    ax[2].set_title("surface pressure, lifting cylinder")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(out_dir / "potential_flow.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'potential_flow.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
