"""P10.03 The Kutta Condition & Circulation: starter code.

Lift in a steady, inviscid, incompressible 2-D flow is L' = rho V Gamma (Kutta-Joukowski). This
module builds the exact potential flows that show it, and checks the theorem four independent ways:

* the cylinder with circulation: complex potential, surface velocity and Cp, stagnation points;
* conformal maps from the circle plane (zeta) to the airfoil plane (z): Joukowski and
  Karman-Trefftz (a finite trailing-edge angle tau), with their derivatives and inverse;
* the Kutta condition: the trailing-edge point of the circle is a stagnation point, which fixes
  Gamma = 4 pi V a sin(alpha + beta);
* exact airfoil surface velocity and Cp;
* lift from (1) Kutta-Joukowski, (2) integrating the surface pressure, (3) Blasius' theorem on a
  contour away from the body, (4) a momentum balance on a control surface; circulation by the line
  integral of velocity (`P01.11`) around any contour that encloses the airfoil.

Sign convention (Phase 10): Gamma is positive CLOCKWISE for a lifting airfoil in a left-to-right
stream, so that L' = rho V Gamma > 0. In complex form w(zeta) = V (zeta' e^{-i alpha} +
a^2 e^{i alpha}/zeta') + i Gamma/(2 pi) ln zeta', zeta' = zeta - mu; a counter-clockwise line
integral of velocity returns -Gamma.

Units: SI. Lengths may be in chords; V in m/s; alpha in radians unless a name says _deg.

Run:  python main.py [out_dir]   (prints the lesson's numbers, writes kutta_circulation.png)
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy import optimize

RHO = 1.225
TRAINER = dict(mass=750.0, S=12.0, b=10.0, c=1.2, V=55.0)
G0 = 9.81


# ---------------------------------------------------------------- the cylinder with circulation

def circle_dwdz(zeta: np.ndarray, V: float, alpha: float, a: float, Gamma: float,
                mu: complex = 0.0) -> np.ndarray:
    """Complex velocity u - i v = dw/dzeta of a uniform stream V at angle alpha past a circle of
    radius a centred at mu, with clockwise circulation Gamma:

        dw/dzeta = V (e^{-i alpha} - a^2 e^{i alpha} / zeta'^2) + i Gamma / (2 pi zeta'),  zeta' = zeta - mu."""
    raise NotImplementedError


def circle_surface_speed(theta: np.ndarray, V: float, alpha: float, a: float, Gamma: float) -> np.ndarray:
    """Tangential (counter-clockwise positive) velocity on the circle at polar angle theta:
    u_theta = -2 V sin(theta - alpha) - Gamma / (2 pi a)."""
    raise NotImplementedError


def circle_surface_cp(theta: np.ndarray, V: float, alpha: float, a: float, Gamma: float) -> np.ndarray:
    """Exact surface pressure coefficient on the circle, Cp = 1 - (u_theta / V)^2."""
    raise NotImplementedError


def stagnation_angles(V: float, alpha: float, a: float, Gamma: float) -> np.ndarray:
    """Polar angles of the two surface stagnation points, from u_theta = 0:
    sin(theta - alpha) = -Gamma / (4 pi a V). Returns an empty array when |Gamma| > 4 pi a V (the
    single stagnation point has left the body)."""
    raise NotImplementedError


def kutta_circulation(V: float, a: float, alpha: float, beta: float) -> float:
    """Circulation that puts a stagnation point at the circle angle -beta (the image of the sharp
    trailing edge): u_theta(-beta) = 0 gives Gamma = 4 pi V a sin(alpha + beta)."""
    raise NotImplementedError


# ---------------------------------------------------------------- conformal maps

def circle_geometry(b: float, mu: complex) -> tuple[float, float]:
    """Radius a and angle beta of the circle centred at mu that passes through the critical point
    zeta = b: b - mu = a e^{-i beta}. beta > 0 when the centre is above the real axis (camber)."""
    raise NotImplementedError


def karman_trefftz(zeta: np.ndarray, b: float, tau: float = 0.0) -> np.ndarray:
    """Karman-Trefftz map with trailing-edge angle tau (radians), n = 2 - tau/pi:

        (z - n b)/(z + n b) = ((zeta - b)/(zeta + b))^n,   i.e.  z = n b (1 + r^n)/(1 - r^n),
        r = (zeta - b)/(zeta + b).
    n = 2 (tau = 0) is the Joukowski map z = zeta + b^2/zeta. The principal power is used; r is
    never on the negative real axis outside the circle, so the map is analytic there."""
    raise NotImplementedError


def karman_trefftz_dzdzeta(zeta: np.ndarray, b: float, tau: float = 0.0) -> np.ndarray:
    """dz/dzeta of the Karman-Trefftz map: 4 n^2 b^2 r^{n-1} / ((1 - r^n)^2 (zeta + b)^2).
    It vanishes at zeta = b (the trailing edge), which is why the velocity there is finite only if
    the circle flow has a stagnation point there."""
    raise NotImplementedError


def karman_trefftz_inverse(z: np.ndarray, b: float, tau: float, mu: complex) -> np.ndarray:
    """zeta(z) outside the circle: invert r^n = (z - n b)/(z + n b) with the principal root, then
    polish with Newton steps on karman_trefftz(zeta) = z. Valid for contours that keep away from the
    trailing edge (checked: the result lies outside the circle)."""
    n = 2.0 - tau / np.pi
    z = np.asarray(z, complex)
    r = ((z - n * b) / (z + n * b)) ** (1.0 / n)
    zeta = b * (1 + r) / (1 - r)
    for _ in range(3):
        zeta = zeta - (karman_trefftz(zeta, b, tau) - z) / karman_trefftz_dzdzeta(zeta, b, tau)
    a, _ = circle_geometry(b, mu)
    if np.any(np.abs(zeta - mu) < a * (1 - 1e-9)):
        raise ValueError("contour too close to the airfoil for the principal inverse")
    return zeta


# ---------------------------------------------------------------- the exact airfoil flow

def airfoil_flow(b: float, mu: complex, tau: float, alpha: float, V: float = 1.0, n: int = 256,
                 Gamma: float | None = None) -> dict:
    """Exact flow on the Karman-Trefftz airfoil from the circle (centre mu, through zeta = b).

    Circle angles theta_k = -beta + 2 pi (k + 1/2)/n (half a step from the trailing edge, so the
    0/0 limit there is never evaluated), counter-clockwise. Gamma defaults to the Kutta value.
    Returns z (surface points), dz/dtheta, the complex velocity u - i v on the airfoil, speed q,
    cp = 1 - (q/V)^2, theta, Gamma, a, beta and the chord (TE to the farthest surface point)."""
    raise NotImplementedError


def pressure_force_coefficients(z_dtheta: np.ndarray, cp: np.ndarray, alpha: float, chord: float,
                                dtheta: float) -> tuple[float, float]:
    """(cl, cd) by integrating surface pressure around a counter-clockwise surface parameterised by
    an angle with uniform step dtheta: the force on the body is -oint p n ds with n ds = (dy, -dx),
    so c_x = -oint cp dy / c and c_y = oint cp dx / c; then rotate to wind axes:
    cl = c_y cos(alpha) - c_x sin(alpha), cd = c_x cos(alpha) + c_y sin(alpha).
    The trapezoid rule on a smooth periodic integrand converges spectrally."""
    raise NotImplementedError


def blasius_force(w_on_contour: np.ndarray, dz: np.ndarray, rho: float, alpha: float) -> tuple[float, float]:
    """Blasius' theorem X - i Y = (i rho / 2) oint (dw/dz)^2 dz on a counter-clockwise contour, given
    the complex velocity there and the contour increments dz. Returns (lift, drag) per unit span."""
    raise NotImplementedError


def momentum_force(w_on_contour: np.ndarray, dz: np.ndarray, rho: float, V: float,
                   alpha: float) -> tuple[float, float]:
    """Force on the body from a momentum balance on a counter-clockwise control surface:
    F = -oint (p n + rho u (u . n)) ds, with p - p_inf = rho (V^2 - |u|^2)/2 (Bernoulli) and
    n ds = (dy, -dx). Returns (lift, drag) per unit span."""
    raise NotImplementedError


def circulation_line_integral(w_on_contour: np.ndarray, dz: np.ndarray) -> float:
    """Clockwise circulation from a counter-clockwise contour: oint u . dl = Re oint (dw/dz) dz is
    the counter-clockwise circulation, so Gamma = -Re oint (dw/dz) dz."""
    raise NotImplementedError


def contour_in_z(kind: str, size: float, center: complex = 0.0, n: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """Counter-clockwise closed contour in the airfoil plane (given): 'circle' of radius size,
    'ellipse' with semi-axes size and 0.4 size (trapezoid rule on a periodic parameter, n points),
    or 'rectangle' with half-sides size and 0.5 size (n/4 Gauss-Legendre points on each straight
    edge). Returns nodes z_k and weighted increments dz_k, so that oint f dz ~ sum f(z_k) dz_k."""
    if kind in ("circle", "ellipse"):
        t = 2 * np.pi * (np.arange(n) + 0.5) / n
        ry = size if kind == "circle" else 0.4 * size
        z = size * np.cos(t) + 1j * ry * np.sin(t)
        dzdt = -size * np.sin(t) + 1j * ry * np.cos(t)
        return z + center, dzdt * (2 * np.pi / n)
    if kind == "rectangle":
        hx, hy = size, 0.5 * size
        corners = [hx - 1j * hy, hx + 1j * hy, -hx + 1j * hy, -hx - 1j * hy, hx - 1j * hy]
        xg, wg = np.polynomial.legendre.leggauss(max(n // 4, 4))
        zs, dzs = [], []
        for p0, p1 in zip(corners[:-1], corners[1:]):
            zs.append(0.5 * (p0 + p1) + 0.5 * (p1 - p0) * xg)
            dzs.append(0.5 * (p1 - p0) * wg)
        return np.concatenate(zs) + center, np.concatenate(dzs)
    raise ValueError(kind)


def airfoil_velocity_at(z: np.ndarray, b: float, mu: complex, tau: float, alpha: float, V: float,
                        Gamma: float) -> np.ndarray:
    """Complex velocity u - i v at airfoil-plane points z (outside the airfoil): map back to the
    circle plane and divide by dz/dzeta."""
    zeta = karman_trefftz_inverse(z, b, tau, mu)
    a, _ = circle_geometry(b, mu)
    return circle_dwdz(zeta, V, alpha, a, Gamma, mu) / karman_trefftz_dzdzeta(zeta, b, tau)


# ---------------------------------------------------------------- shape fitting

def section_shape(b: float, mu: complex, tau: float, n: int = 2048) -> dict:
    """Chord, chord-line inclination phi (nose-up positive, so alpha_chord = alpha + phi), maximum
    thickness and camber (fractions of chord, measured normal to the chord line) of the
    Karman-Trefftz section (given)."""
    a, beta = circle_geometry(b, mu)
    theta = -beta + 2 * np.pi * (np.arange(n) + 0.5) / n
    z = karman_trefftz(mu + a * np.exp(1j * theta), b, tau)
    z_te = karman_trefftz(np.array([b + 0j]), b, tau)[0]
    i_le = int(np.argmax(np.abs(z - z_te)))
    z_le = z[i_le]
    chord = abs(z_te - z_le)
    phi = np.angle(z_te - z_le)            # rotation of the chord line from +x
    # chord frame: LE at 0, TE at 1, nose-up chord inclination removed
    zc = (z - z_le) * np.exp(-1j * phi) / chord
    up, lo = zc[:i_le + 1], zc[i_le:]      # theta increasing: TE -> upper -> LE -> lower -> TE
    xu, yu = up.real[::-1], up.imag[::-1]
    xl, yl = lo.real, lo.imag
    xs = 0.5 * (1 - np.cos(np.linspace(0, np.pi, 801)))
    ou, ol = np.argsort(xu), np.argsort(xl)
    Yu = np.interp(xs, xu[ou], yu[ou])
    Yl = np.interp(xs, xl[ol], yl[ol])
    th, cam = Yu - Yl, 0.5 * (Yu + Yl)
    return dict(chord=chord, phi=-phi, t=float(th.max()), x_t=float(xs[np.argmax(th)]), m=float(cam.max()),
                x_m=float(xs[np.argmax(cam)]), xs=xs, yu=Yu, yl=Yl)


def fit_karman_trefftz(t_target: float, m_target: float, tau: float, b: float = 1.0) -> complex:
    """Circle centre mu = -eps_x + i eps_y (b = 1) whose Karman-Trefftz section has the target
    maximum thickness and camber (fractions of chord), by fsolve from a thin-section starting guess
    (given). The trailing-edge angle alone gives some thickness, so small t with large tau has no
    solution; for symmetric sections the camber equation is satisfied by eps_y = 0."""
    tfun = lambda ex, ey: section_shape(b, complex(-ex, ey), tau, 1024)
    ex = optimize.brentq(lambda e: tfun(e, 0.0)["t"] - t_target, -0.2 * b, 0.5 * b, xtol=1e-13)
    if m_target == 0.0:
        return complex(-ex, 0.0)

    def res(p):
        sh = tfun(p[0], p[1])
        return [sh["t"] - t_target, sh["m"] - m_target]
    sol = optimize.fsolve(res, [ex, 2.0 * m_target * b], xtol=1e-12)
    return complex(-sol[0], sol[1])


# ---------------------------------------------------------------- main

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's numbers and save kutta_circulation.png."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    t = TRAINER
    W = t["mass"] * G0
    Lp = W / t["b"]
    Gam = Lp / (RHO * t["V"])
    cl_need = 2 * Gam / (t["V"] * t["c"])
    print(f"trainer cruise: lift per metre span W/b = {Lp:.2f} N/m -> Gamma = L'/(rho V) = {Gam:.3f} m^2/s, "
          f"cl = 2 Gamma/(V c) = {cl_need:.4f} (uniform-loading estimate; finite wings in P10.09-P10.10)")

    # cylinder
    a, V = 1.0, 1.0
    for G in (0.0, 2 * np.pi, 4 * np.pi, 6 * np.pi):
        st = stagnation_angles(V, 0.0, a, G)
        th = 2 * np.pi * (np.arange(512) + 0.5) / 512
        z = a * np.exp(1j * th)
        cl, cd = pressure_force_coefficients(1j * z, circle_surface_cp(th, V, 0.0, a, G), 0.0, 2 * a, 2 * np.pi / 512)
        print(f"cylinder Gamma/(4 pi a V) = {G / (4 * np.pi):.2f}: stagnation at "
              f"{np.array2string(np.degrees(st), precision=1)} deg; cl(pressure) = {cl:.6f}, "
              f"rho V Gamma/(q 2a) = {2 * G / (2 * a):.6f}, cd = {cd:.1e}")

    # trainer-like Karman-Trefftz section: 12 % thick, 2 % camber, 16 deg trailing-edge angle (as NACA 2412)
    tau = np.radians(16.0)
    mu = fit_karman_trefftz(0.12, 0.02, tau)
    sh = section_shape(1.0, mu, tau)
    a, beta = circle_geometry(1.0, mu)
    print(f"\nKarman-Trefftz section (tau = 16 deg): mu = {mu.real:.5f} + {mu.imag:.5f} i (b = 1), a = {a:.5f}, "
          f"beta = {np.degrees(beta):.4f} deg; chord {sh['chord']:.4f}, t = {sh['t']:.4f} at {sh['x_t']:.3f}, "
          f"m = {sh['m']:.4f} at {sh['x_m']:.3f}, chord-line inclination {np.degrees(sh['phi']):.4f} deg")
    c = sh["chord"]
    # cl = 8 pi a sin(alpha + beta)/c, alpha from the x-axis; chord-line alpha = alpha_x + phi
    alpha_x = np.arcsin(cl_need * c / (8 * np.pi * a)) - beta
    alpha_c = alpha_x + sh["phi"]
    aL0_c = -beta + sh["phi"]
    slope = 8 * np.pi * a * np.cos(alpha_x + beta) / c
    print(f"zero-lift angle (chord line) {np.degrees(aL0_c):.3f} deg; cruise cl {cl_need:.4f} at alpha = "
          f"{np.degrees(alpha_c):.3f} deg (chord line); lift slope {slope:.4f} /rad = {slope / (2 * np.pi):.4f} x 2 pi")
    f = airfoil_flow(1.0, mu, tau, alpha_x, 1.0, 512)
    cl_p, cd_p = pressure_force_coefficients(f["dzdtheta"], f["cp"], alpha_x, c, 2 * np.pi / 512)
    print(f"exact flow, n = 512: cl from Kutta-Joukowski {2 * f['Gamma'] / c:.10f}, from surface pressure "
          f"{cl_p:.10f}, cd from pressure {cd_p:.1e}; min cp {f['cp'].min():.4f}")
    for n in (16, 32, 64, 128):
        fn = airfoil_flow(1.0, mu, tau, alpha_x, 1.0, n)
        cln, _ = pressure_force_coefficients(fn["dzdtheta"], fn["cp"], alpha_x, c, 2 * np.pi / n)
        print(f"   pressure-integration error with n = {n:4d}: {abs(cln - 2 * f['Gamma'] / c):.2e}")
    for kind, size in (("circle", 1.5 * c), ("ellipse", 0.8 * c), ("rectangle", 0.7 * c)):
        zc, dz = contour_in_z(kind, size, center=0.05 * c)
        wc = airfoil_velocity_at(zc, 1.0, mu, tau, alpha_x, 1.0, f["Gamma"])
        G = circulation_line_integral(wc, dz)
        Lb, Db = blasius_force(wc, dz, 1.0, alpha_x)
        Lm, Dm = momentum_force(wc, dz, 1.0, 1.0, alpha_x)
        print(f"  {kind:9s} contour: Gamma = {G:.8f} (Kutta {f['Gamma']:.8f}); Blasius L' = {Lb:.8f}, D' = {Db:.1e}; "
              f"momentum L' = {Lm:.8f}, D' = {Dm:.1e}  (rho V Gamma = {f['Gamma']:.8f})")
    # the Kutta condition in numbers: velocity near the trailing edge with and without it
    fK = airfoil_flow(1.0, mu, tau, alpha_x, 1.0, 4096)
    f0 = airfoil_flow(1.0, mu, tau, alpha_x, 1.0, 4096, Gamma=0.0)
    print(f"speed at the surface point nearest the trailing edge (n = 4096): with Kutta {fK['q'][0]:.4f} V, "
          f"with Gamma = 0 {f0['q'][0]:.1f} V")
    # thickness effect on the lift slope (symmetric KT, tau = 16 deg)
    print("lift slope / 2 pi for symmetric Joukowski sections (tau = 0), against 1 + 0.77 t/c:")
    for tt in (0.04, 0.08, 0.12, 0.16, 0.20):
        mus = fit_karman_trefftz(tt, 0.0, 0.0)
        ss = section_shape(1.0, mus, 0.0)
        aa, _ = circle_geometry(1.0, mus)
        print(f"   t/c = {tt:.2f}: {8 * np.pi * aa / ss['chord'] / (2 * np.pi):.4f}  (1 + 0.77 t/c = {1 + 0.77 * tt:.4f})")

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    th = np.linspace(0, 2 * np.pi, 361)
    for G in (0.0, 2 * np.pi, 4 * np.pi):
        ax[0].plot(np.degrees(th), circle_surface_cp(th, 1.0, 0.0, 1.0, G), label=f"Gamma/(4 pi a V) = {G / (4 * np.pi):.1f}")
    ax[0].invert_yaxis(); ax[0].set_xlabel("theta (deg)"); ax[0].set_ylabel("Cp"); ax[0].legend()
    ax[0].set_title("cylinder with circulation")
    zeta_g = mu + (a * 1.0001 + np.linspace(0, 3, 120)[:, None]) * np.exp(1j * np.linspace(0, 2 * np.pi, 361)[None, :])
    zz = karman_trefftz(zeta_g, 1.0, tau)
    psi = np.imag(1.0 * ((zeta_g - mu) * np.exp(-1j * alpha_x) + a ** 2 * np.exp(1j * alpha_x) / (zeta_g - mu))
                  + 1j * f["Gamma"] / (2 * np.pi) * np.log(zeta_g - mu))
    ax[1].contour(zz.real, zz.imag, psi, levels=30, linewidths=0.6)
    ax[1].fill(f["z"].real, f["z"].imag, color="0.6")
    ax[1].set_aspect("equal"); ax[1].set_xlim(-3, 3); ax[1].set_ylim(-1.5, 1.5)
    ax[1].set_title(f"streamlines, Kutta condition, alpha = {np.degrees(alpha_c):.2f} deg")
    xs_c = ((f["z"] - f["z"][np.argmax(np.abs(f["z"] - f["z_te"]))]) * np.exp(1j * sh["phi"])).real / c
    ax[2].plot(xs_c, f["cp"], ".", ms=2)
    ax[2].invert_yaxis(); ax[2].set_xlabel("x/c"); ax[2].set_ylabel("Cp")
    ax[2].set_title("exact Karman-Trefftz surface Cp (trainer cruise cl)")
    fig.tight_layout()
    fig.savefig(out_dir / "kutta_circulation.png", dpi=120)
    plt.close(fig)
    print(f"\nwrote {out_dir / 'kutta_circulation.png'}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else ".")
