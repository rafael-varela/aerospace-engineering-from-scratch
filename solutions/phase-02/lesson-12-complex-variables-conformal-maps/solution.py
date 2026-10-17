"""P02.12 Complex Variables & Conformal Maps: reference solution.

A Joukowski airfoil from a circle: analytic functions checked by the Cauchy-Riemann equations, the
complex potential of flow past a cylinder with circulation, the Kutta condition, the Joukowski map
z = zeta + c^2 / zeta, surface velocity and pressure on the airfoil, the Kutta-Joukowski lift, and
contour integrals checked against the residue theorem.

Conventions:
    * z = x + i y is the physical (airfoil) plane, zeta the circle plane. Lengths in metres.
    * Complex potential F = phi + i psi; complex velocity w = dF/dz = u - i v [m/s].
    * Flow past a cylinder of radius R centred at zeta_0, free stream U at angle alpha, circulation
      Gamma (positive clockwise, so that positive Gamma gives positive lift at positive alpha):
          F(zeta) = U [ (zeta - zeta_0) e^{-i alpha} + R^2 e^{i alpha} / (zeta - zeta_0) ]
                    + i Gamma / (2 pi) * log(zeta - zeta_0).
    * The Joukowski circle passes through zeta = c (the trailing edge maps to z = 2c) and is centred at
      zeta_0 = x0 + i y0 with x0 < 0 (thickness) and y0 > 0 (camber). R = |c - zeta_0|.
    * Kutta condition: the rear stagnation point sits at zeta = c, which fixes
          Gamma = 4 pi U R sin(alpha + beta),  beta = atan2(y0, c - x0).

Run:  python solution.py          (prints the airfoil's geometry, lift and checks, writes joukowski_airfoil.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-12-complex-variables-conformal-maps
"""
from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

import numpy as np


# GIVEN
def lesson_airfoil() -> dict[str, float]:
    """Return the lesson's airfoil and flow: c = 0.25 m, centre (-0.025, 0.025) m, U = 50 m/s, alpha = 4 deg, rho = 1.225 kg/m^3.

    The circle radius is R = |c - zeta_0| = 0.2761 m and the chord about 1.03 m: a general-aviation-size
    section with roughly 13 % thickness and 5 % camber. Plausible, not a catalogue section.
    """
    return {"c": 0.25, "x0": -0.025, "y0": 0.025, "U": 50.0, "alpha": np.deg2rad(4.0), "rho": 1.225}


# ---------------------------------------------------------------- analytic functions

def cauchy_riemann_residual(f: Callable[[np.ndarray], np.ndarray], z: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """Return |df/dx + i df/dy| at z by central differences: zero (to O(h^2)) exactly where f is analytic.

    With f = u + i v, df/dx + i df/dy = (u_x - v_y) + i (v_x + u_y): both Cauchy-Riemann equations at once.
    """
    z = np.asarray(z, dtype=complex)
    fx = (f(z + h) - f(z - h)) / (2 * h)
    fy = (f(z + 1j * h) - f(z - 1j * h)) / (2 * h)
    return np.abs(fx + 1j * fy)


def complex_derivative(f: Callable[[np.ndarray], np.ndarray], z: np.ndarray, h: float = 1e-6) -> np.ndarray:
    """Return df/dz by the central difference along the real direction, (f(z + h) - f(z - h)) / 2h.

    For an analytic f the derivative is the same in every direction, so this is df/dz; for a
    non-analytic f it is merely df/dx.
    """
    z = np.asarray(z, dtype=complex)
    return (f(z + h) - f(z - h)) / (2 * h)


def laplacian_of_real_part(f: Callable[[np.ndarray], np.ndarray], z: np.ndarray, h: float = 1e-4) -> np.ndarray:
    """Return u_xx + u_yy for u = Re f at z by the five-point stencil: zero where f is analytic (u is harmonic)."""
    z = np.asarray(z, dtype=complex)
    u = lambda w: np.real(f(w))
    return (u(z + h) + u(z - h) + u(z + 1j * h) + u(z - 1j * h) - 4 * u(z)) / h ** 2


def conformality_defect(f: Callable[[np.ndarray], np.ndarray], z: complex, theta1: float, theta2: float, h: float = 1e-6) -> float:
    """Return |angle between the images of two short segments from z| - |theta2 - theta1| (radians).

    Segments leave z in directions theta1 and theta2; their images leave f(z) in directions
    arg(f(z + h e^{i theta}) - f(z)). An analytic map with f'(z) != 0 preserves the angle, so the
    defect is O(h).
    """
    z = complex(z)
    d1 = f(z + h * np.exp(1j * theta1)) - f(z)
    d2 = f(z + h * np.exp(1j * theta2)) - f(z)
    img = np.angle(d2 / d1)
    orig = np.angle(np.exp(1j * (theta2 - theta1)))
    return float(abs(img - orig))


# ---------------------------------------------------------------- flow past a cylinder

def cylinder_potential(zeta: np.ndarray, U: float, alpha: float, R: float, Gamma: float, zeta0: complex = 0.0) -> np.ndarray:
    """Return F(zeta) = U[(zeta - zeta0) e^{-i alpha} + R^2 e^{i alpha} / (zeta - zeta0)] + i Gamma/(2 pi) log(zeta - zeta0).

    Gamma > 0 is a clockwise circulation (the sign that lifts at positive alpha).
    """
    s = np.asarray(zeta, dtype=complex) - zeta0
    return U * (s * np.exp(-1j * alpha) + R ** 2 * np.exp(1j * alpha) / s) + 1j * Gamma / (2 * np.pi) * np.log(s)


def cylinder_velocity(zeta: np.ndarray, U: float, alpha: float, R: float, Gamma: float, zeta0: complex = 0.0) -> np.ndarray:
    """Return w = dF/dzeta = U[e^{-i alpha} - R^2 e^{i alpha} / (zeta - zeta0)^2] + i Gamma / (2 pi (zeta - zeta0)), i.e. u - i v."""
    s = np.asarray(zeta, dtype=complex) - zeta0
    return U * (np.exp(-1j * alpha) - R ** 2 * np.exp(1j * alpha) / s ** 2) + 1j * Gamma / (2 * np.pi * s)


def kutta_circulation(U: float, R: float, alpha: float, beta: float) -> float:
    """Return Gamma = 4 pi U R sin(alpha + beta): the circulation that puts the rear stagnation point at angle -beta.

    beta is the angle below the centre at which the trailing-edge point sits on the circle,
    beta = atan2(y0, c - x0) for the Joukowski circle.
    """
    return 4.0 * np.pi * U * R * np.sin(alpha + beta)


def stagnation_points(U: float, alpha: float, R: float, Gamma: float, zeta0: complex = 0.0) -> np.ndarray:
    """Return the two zeros of the complex velocity (the stagnation points) in the zeta plane.

    Solving w = 0: s^2 e^{-i alpha} + i Gamma s / (2 pi U) - R^2 e^{i alpha} = 0, a quadratic in s = zeta - zeta0.
    For Gamma > 4 pi U R the roots leave the circle (one inside, one outside).
    """
    a = np.exp(-1j * alpha)
    b = 1j * Gamma / (2 * np.pi * U)
    cc = -R ** 2 * np.exp(1j * alpha)
    disc = np.sqrt(b * b - 4 * a * cc + 0j)
    return np.array([(-b + disc) / (2 * a), (-b - disc) / (2 * a)]) + zeta0


# ---------------------------------------------------------------- the Joukowski map

def joukowski(zeta: np.ndarray, c: float) -> np.ndarray:
    """Return z = zeta + c^2 / zeta. Maps the circle |zeta| = c to the slit [-2c, 2c]; critical points zeta = +-c."""
    zeta = np.asarray(zeta, dtype=complex)
    return zeta + c ** 2 / zeta


def joukowski_derivative(zeta: np.ndarray, c: float) -> np.ndarray:
    """Return dz/dzeta = 1 - c^2 / zeta^2, which vanishes at zeta = +-c (the map is not conformal there)."""
    zeta = np.asarray(zeta, dtype=complex)
    return 1.0 - c ** 2 / zeta ** 2


def airfoil_circle(c: float, x0: float, y0: float) -> tuple[complex, float]:
    """Return (zeta0, R): the centre and radius of the circle through zeta = c centred at x0 + i y0.

    Raise ``ValueError`` unless x0 < 0 <= y0 (the circle must enclose zeta = -c for a closed airfoil).
    """
    if not (x0 < 0.0 and y0 >= 0.0):
        raise ValueError("need x0 < 0 and y0 >= 0 for a Joukowski airfoil")
    zeta0 = complex(x0, y0)
    return zeta0, abs(c - zeta0)


def airfoil_coordinates(c: float, x0: float, y0: float, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
    """Return (zeta, z): n points on the circle (starting and ending at the trailing edge) and their images.

    The circle is traversed counter-clockwise from zeta = c; the trailing edge is z = 2c exactly.
    """
    zeta0, R = airfoil_circle(c, x0, y0)
    beta = np.arctan2(y0, c - x0)
    theta = -beta + np.linspace(0.0, 2 * np.pi, n)
    zeta = zeta0 + R * np.exp(1j * theta)
    zeta[0] = zeta[-1] = c
    return zeta, joukowski(zeta, c)


def airfoil_geometry(z: np.ndarray) -> dict[str, float]:
    """Return chord, maximum thickness and maximum camber (all in metres) of a closed airfoil polygon z.

    Chord: x_max - x_min. The upper and lower surfaces are split at the leading edge (minimum x) and
    interpolated onto a common x grid; thickness is the largest gap, camber the largest mean-line height.
    """
    x, y = np.real(z), np.imag(z)
    i_le = int(np.argmin(x))
    first, second = (x[: i_le + 1], y[: i_le + 1]), (x[i_le:], y[i_le:])
    # Traversed counter-clockwise from the trailing edge: the first arc is the upper surface.
    xu, yu = first
    xl, yl = second
    order_u, order_l = np.argsort(xu), np.argsort(xl)
    xg = np.linspace(x.min(), x.max(), 2001)[1:-1]
    yu_g = np.interp(xg, xu[order_u], yu[order_u])
    yl_g = np.interp(xg, xl[order_l], yl[order_l])
    return {"chord": float(x.max() - x.min()), "thickness": float(np.max(yu_g - yl_g)),
            "camber": float(np.max(0.5 * (yu_g + yl_g)) - 0.0)}


# ---------------------------------------------------------------- the airfoil flow

def airfoil_surface_velocity(zeta: np.ndarray, c: float, x0: float, y0: float, U: float, alpha: float, Gamma: float) -> np.ndarray:
    """Return the complex velocity w_z = w_zeta / (dz/dzeta) at circle points zeta, mapped to the airfoil.

    At the trailing edge dz/dzeta = 0; with the Kutta circulation w_zeta = 0 there too, and the limit is
    finite: it is evaluated by l'Hopital as dw_zeta/dzeta / d^2z/dzeta^2. Without the Kutta circulation
    the velocity there is infinite and this function returns inf.
    """
    zeta0, R = airfoil_circle(c, x0, y0)
    zeta = np.asarray(zeta, dtype=complex)
    w_zeta = cylinder_velocity(zeta, U, alpha, R, Gamma, zeta0)
    dz = joukowski_derivative(zeta, c)
    out = np.empty_like(w_zeta)
    at_te = np.abs(zeta - c) < 1e-12 * max(c, 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        out[~at_te] = w_zeta[~at_te] / dz[~at_te]
    if np.any(at_te):
        s = c - zeta0
        dw = U * 2 * R ** 2 * np.exp(1j * alpha) / s ** 3 - 1j * Gamma / (2 * np.pi * s ** 2)   # dw_zeta/dzeta
        d2z = 2 * c ** 2 / c ** 3                                                             # d^2 z / dzeta^2 at zeta = c
        if abs(w_zeta[at_te][0]) < 1e-9 * U:
            out[at_te] = dw / d2z
        else:
            out[at_te] = np.inf
    return out


def pressure_coefficient(w: np.ndarray, U: float) -> np.ndarray:
    """Return C_p = 1 - |w|^2 / U^2 (Bernoulli, steady incompressible flow)."""
    return 1.0 - np.abs(np.asarray(w, dtype=complex)) ** 2 / U ** 2


def lift_per_span(rho: float, U: float, Gamma: float) -> float:
    """Return the Kutta-Joukowski lift per unit span, L' = rho U Gamma [N/m]."""
    return rho * U * Gamma


def lift_coefficient(Gamma: float, U: float, chord: float) -> float:
    """Return C_l = L' / (1/2 rho U^2 chord) = 2 Gamma / (U chord)."""
    if chord <= 0.0 or U <= 0.0:
        raise ValueError("chord and U must be positive")
    return 2.0 * Gamma / (U * chord)


# ---------------------------------------------------------------- contour integrals

def circle_path(center: complex, radius: float, n: int = 2000) -> np.ndarray:
    """Return n + 1 points on the circle |z - center| = radius, counter-clockwise, closed (last = first exactly)."""
    t = np.linspace(0.0, 2 * np.pi, n + 1)
    path = center + radius * np.exp(1j * t)
    path[-1] = path[0]
    return path


def contour_integral(f: Callable[[np.ndarray], np.ndarray], path: np.ndarray) -> complex:
    """Return the integral of f along the path.

    Closed path (last point equals the first): the points are taken as samples z(t_k) of a smooth periodic
    curve at equal parameter spacing; z'(t) is obtained by spectral (FFT) differentiation and the integral
    is the trapezoidal sum of f(z) z' dt, which converges faster than any power of the spacing when f is
    analytic near the path: 2000 points give machine precision. Open path: the chord trapezoidal rule,
    sum (f(z_k) + f(z_{k+1})) (z_{k+1} - z_k) / 2, accurate to second order.
    """
    path = np.asarray(path, dtype=complex)
    if path.size > 2 and path[0] == path[-1]:
        z = path[:-1]
        n = z.size
        k = np.fft.fftfreq(n, 1.0 / n)                     # integer wavenumbers
        dz_dt = np.fft.ifft(1j * k * np.fft.fft(z))        # spectral derivative with respect to t in [0, 2 pi)
        return complex(np.sum(f(z) * dz_dt) * (2 * np.pi / n))
    fz = f(path)
    return complex(np.sum(0.5 * (fz[1:] + fz[:-1]) * np.diff(path)))


def circulation_and_flux(w: Callable[[np.ndarray], np.ndarray], path: np.ndarray) -> tuple[float, float]:
    """Return (Gamma, Q) from the contour integral of the complex velocity: int w dz = Gamma + i Q.

    Gamma is the circulation around the path (counter-clockwise positive in this integral; the lesson's
    clockwise-positive Gamma appears with a minus sign) and Q the net outflow through it.
    """
    val = contour_integral(w, path)
    return float(np.real(val)), float(np.imag(val))


# ---------------------------------------------------------------- main

# GIVEN
def main(out_dir: str | Path = ".") -> None:
    """Build the airfoil, check analyticity and the residue theorem, apply the Kutta condition, plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    p = lesson_airfoil()
    c, x0, y0, U, alpha, rho = p["c"], p["x0"], p["y0"], p["U"], p["alpha"], p["rho"]
    zeta0, R = airfoil_circle(c, x0, y0)
    beta = np.arctan2(y0, c - x0)
    print(f"circle: centre {zeta0:.4f} m, R = {R:.4f} m, beta = {np.degrees(beta):.3f} deg; map constant c = {c} m")

    # Analyticity checks.
    pts = np.array([0.7 + 0.2j, -1.1 + 0.5j, 2.0 - 1.5j])
    for name, f in (("z^2", lambda z: z ** 2), ("exp(z)", np.exp), ("1/z", lambda z: 1 / z), ("conj(z)", np.conj), ("|z|^2", lambda z: np.abs(z) ** 2)):
        print(f"Cauchy-Riemann residual of {name:8s}: {np.max(cauchy_riemann_residual(f, pts)):.2e};"
              f" Laplacian of Re: {np.max(np.abs(laplacian_of_real_part(f, pts))):.2e}")
    print(f"conformality defect of the Joukowski map at zeta = 0.4 + 0.3i: {conformality_defect(lambda z: joukowski(z, c), 0.4 + 0.3j, 0.0, 1.2):.2e} rad;"
          f" at the critical point zeta = c: {conformality_defect(lambda z: joukowski(z, c), c, 0.0, 1.2):.3f} rad (angles double there)")

    # Residue theorem.
    path = circle_path(0.0, 1.0)
    print(f"contour integral of 1/z around |z| = 1: {contour_integral(lambda z: 1 / z, path):.12f} (2 pi i = {2 * np.pi:.12f} i)")
    print(f"contour integral of 1/(z - 3) around |z| = 1: {abs(contour_integral(lambda z: 1 / (z - 3), path)):.2e} (pole outside)")
    print(f"contour integral of z^2 around |z| = 1: {abs(contour_integral(lambda z: z ** 2, path)):.2e} (entire function)")

    # Kutta condition and lift.
    Gamma = kutta_circulation(U, R, alpha, beta)
    zeta, z = airfoil_coordinates(c, x0, y0)
    geom = airfoil_geometry(z)
    print(f"\nairfoil: chord {geom['chord']:.4f} m, thickness {geom['thickness'] / geom['chord'] * 100:.1f} % chord,"
          f" camber {geom['camber'] / geom['chord'] * 100:.1f} % chord; trailing edge at z = {z[0]:.3f}")
    print(f"Kutta circulation at alpha = {np.degrees(alpha):.1f} deg: Gamma = {Gamma:.3f} m^2/s")
    sp = stagnation_points(U, alpha, R, Gamma, zeta0)
    print(f"stagnation points in the zeta plane: {sp[0]:.4f}, {sp[1]:.4f}; trailing edge zeta = {c}")
    print(f"|w_zeta| at the trailing edge: {abs(cylinder_velocity(np.array([c]), U, alpha, R, Gamma, zeta0)[0]):.2e} m/s")
    w = airfoil_surface_velocity(zeta, c, x0, y0, U, alpha, Gamma)
    cp = pressure_coefficient(w, U)
    print(f"surface speed at the trailing edge (finite by l'Hopital): {abs(w[0]):.3f} m/s; max surface speed {np.max(np.abs(w)):.2f} m/s,"
          f" min C_p {cp.min():.3f}, max C_p {cp.max():.4f}")
    w_nokutta = airfoil_surface_velocity(zeta, c, x0, y0, U, alpha, 0.0)
    print(f"without circulation: speed at the trailing edge {abs(w_nokutta[0])} and near it {abs(w_nokutta[1]):.1f} m/s")
    L = lift_per_span(rho, U, Gamma)
    Cl = lift_coefficient(Gamma, U, geom["chord"])
    print(f"lift per span rho U Gamma = {L:.1f} N/m; C_l = {Cl:.4f}; thin-airfoil 2 pi (alpha + beta) = {2 * np.pi * (alpha + beta):.4f}")
    # Lift from the pressure distribution: closed-contour integral of -C_p dz, normal force per unit chord.
    dz = np.diff(z)
    cp_mid = 0.5 * (cp[1:] + cp[:-1])
    force = 1j * 0.5 * rho * U ** 2 * np.sum(cp_mid * dz)      # F_x + i F_y = -closed-int p n ds, n ds = -i dz on a ccw contour
    lift_from_cp = np.imag(force) * np.cos(alpha) - np.real(force) * np.sin(alpha)   # component normal to the free stream
    print(f"lift from integrating C_p around the surface: {lift_from_cp:.1f} N/m (Kutta-Joukowski {L:.1f} N/m)")
    # Circulation from a contour integral of the mapped velocity.
    big = joukowski(circle_path(zeta0, 3.0 * R, 4000), c)
    # velocity in the z plane along the big contour: w_z = w_zeta / (dz/dzeta) at the corresponding zeta points
    zeta_big = circle_path(zeta0, 3.0 * R, 4000)
    w_big = cylinder_velocity(zeta_big, U, alpha, R, Gamma, zeta0) / joukowski_derivative(zeta_big, c)
    # The velocity is known at the mapped points, so integrate w z'(t) dt with the same spectral derivative.
    zb = big[:-1]
    kk = np.fft.fftfreq(zb.size, 1.0 / zb.size)
    val = np.sum(w_big[:-1] * np.fft.ifft(1j * kk * np.fft.fft(zb))) * (2 * np.pi / zb.size)
    print(f"contour integral of w dz around the airfoil (ccw): {val.real:.3f} + {val.imag:.1e} i  (expected -Gamma = {-Gamma:.3f}, flux 0)")
    for a_deg in (0.0, 2.0, 4.0, 6.0, 8.0):
        G = kutta_circulation(U, R, np.deg2rad(a_deg), beta)
        print(f"  alpha = {a_deg:3.0f} deg: C_l = {lift_coefficient(G, U, geom['chord']):.4f}")

    # Figure: streamlines (psi contours) in the z plane, and C_p.
    r = np.linspace(R, 4 * R, 160)
    th = np.linspace(0.0, 2 * np.pi, 361)
    RR, TH = np.meshgrid(r, th, indexing="ij")
    ZETA = zeta0 + RR * np.exp(1j * TH)
    PSI = np.imag(cylinder_potential(ZETA, U, alpha, R, Gamma, zeta0))
    Z = joukowski(ZETA, c)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    ax = axes[0]
    ax.contour(np.real(Z), np.imag(Z), PSI, levels=np.linspace(-1.5 * U * R, 1.5 * U * R, 31), linewidths=0.6)
    ax.fill(np.real(z), np.imag(z), color="0.3")
    ax.set_aspect("equal"); ax.set_xlim(-2.5 * c, 3.5 * c); ax.set_ylim(-2 * c, 2 * c)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_title(f"streamlines, alpha = {np.degrees(alpha):.0f} deg, Kutta condition")
    ax = axes[1]
    ax.plot(np.real(z), cp, lw=1.0)
    ax.invert_yaxis(); ax.set_xlabel("x [m]"); ax.set_ylabel("$C_p$ (inverted)")
    ax.set_title(f"surface pressure, $C_l$ = {Cl:.3f}"); ax.grid(True, lw=0.4)
    fig.suptitle("P02.12  Joukowski airfoil from a circle")
    fig.tight_layout()
    fig.savefig(out_dir / "joukowski_airfoil.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'joukowski_airfoil.png'}")


if __name__ == "__main__":
    main()
