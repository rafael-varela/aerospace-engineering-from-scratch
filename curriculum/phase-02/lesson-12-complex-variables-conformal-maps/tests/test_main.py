"""Public tests for P02.12 Complex Variables & Conformal Maps. They assert properties: Cauchy-Riemann residuals,
harmonic real parts, the residue theorem, conformality away from critical points, the Kutta condition, and
agreement between the Kutta-Joukowski lift and the integrated surface pressure, with SymPy as the oracle
for derivatives.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
import sympy as sp

from conftest import load_impl

impl = load_impl(__file__)

P = impl.lesson_airfoil()
RNG = np.random.default_rng(212)
PTS = np.array([0.7 + 0.2j, -1.1 + 0.5j, 2.0 - 1.5j, 0.3 - 0.9j])


# ---------------------------------------------------------------- analytic functions

def test_cauchy_riemann_residual_vanishes_for_analytic_functions_only():
    for f in (lambda z: z ** 2, np.exp, lambda z: 1 / z, np.sin, lambda z: z ** 3 - 2 * z + 1j):
        assert np.max(impl.cauchy_riemann_residual(f, PTS)) < 1e-7
        assert np.max(np.abs(impl.laplacian_of_real_part(f, PTS))) < 1e-5          # Re f is harmonic
    assert np.min(impl.cauchy_riemann_residual(np.conj, PTS)) == pytest.approx(2.0, rel=1e-6)   # d/dx conj = 1, d/dy conj = -i
    assert np.min(impl.cauchy_riemann_residual(lambda z: np.abs(z) ** 2, PTS)) > 0.5
    assert np.min(np.abs(impl.laplacian_of_real_part(lambda z: np.abs(z) ** 2, PTS))) == pytest.approx(4.0, rel=1e-4)
    # Im f is harmonic too (the conjugate harmonic): check via -i f.
    assert np.max(np.abs(impl.laplacian_of_real_part(lambda z: -1j * np.exp(z), PTS))) < 1e-5


def test_complex_derivative_matches_sympy():
    z = sp.symbols("z")
    for expr in (z ** 3 - 2 * z, sp.exp(z) / z, sp.log(z) * z):
        f = sp.lambdify(z, expr, "numpy")
        df = sp.lambdify(z, sp.diff(expr, z), "numpy")
        npt.assert_allclose(impl.complex_derivative(f, PTS), df(PTS), rtol=1e-6)
    # The cylinder velocity is dF/dzeta exactly.
    U, a, R, G = 50.0, 0.1, 0.3, 20.0
    F = lambda s: impl.cylinder_potential(s, U, a, R, G)
    npt.assert_allclose(impl.complex_derivative(F, PTS), impl.cylinder_velocity(PTS, U, a, R, G), rtol=1e-6)
    zs = sp.symbols("zeta")
    Fs = U * (zs * sp.exp(-sp.I * a) + R ** 2 * sp.exp(sp.I * a) / zs) + sp.I * G / (2 * sp.pi) * sp.log(zs)
    ws = sp.lambdify(zs, sp.diff(Fs, zs), "numpy")
    npt.assert_allclose(ws(PTS), impl.cylinder_velocity(PTS, U, a, R, G), rtol=1e-12)


def test_joukowski_map_is_conformal_except_at_its_critical_points():
    c = P["c"]
    f = lambda s: impl.joukowski(s, c)
    for s0 in (0.4 + 0.3j, -0.5 + 0.1j, 1.0 + 1.0j, 0.0 + 0.6j):
        for th1, th2 in ((0.0, 1.2), (0.3, -0.9), (2.0, 2.5)):
            assert impl.conformality_defect(f, s0, th1, th2) < 1e-6
    # At zeta = c the derivative vanishes and angles double: a right angle becomes a straight line.
    assert impl.conformality_defect(f, c, 0.0, np.pi / 2) == pytest.approx(np.pi / 2, abs=1e-4)
    npt.assert_allclose(impl.joukowski_derivative(np.array([c, -c]), c), 0.0, atol=1e-15)
    npt.assert_allclose(impl.joukowski_derivative(PTS, c), impl.complex_derivative(f, PTS), rtol=1e-7)
    # A non-analytic map is not conformal.
    assert impl.conformality_defect(np.conj, 0.4 + 0.3j, 0.0, 1.2) == pytest.approx(2.4, abs=1e-6)


def test_joukowski_maps_the_circle_of_radius_c_to_a_slit_and_is_symmetric():
    c = 0.25
    th = np.linspace(0, 2 * np.pi, 201)
    z = impl.joukowski(c * np.exp(1j * th), c)
    npt.assert_allclose(np.imag(z), 0.0, atol=1e-15)
    assert z.real.min() == pytest.approx(-2 * c) and z.real.max() == pytest.approx(2 * c)
    # Large zeta: z ~ zeta (identity at infinity), so the far field is undisturbed.
    big = 1e3 * np.exp(1j * th[:50])
    npt.assert_allclose(impl.joukowski(big, c), big, rtol=1e-6)
    # zeta and c^2 / zeta have the same image.
    npt.assert_allclose(impl.joukowski(PTS, c), impl.joukowski(c ** 2 / PTS, c), rtol=1e-12)


# ---------------------------------------------------------------- the airfoil

def test_airfoil_is_closed_with_a_sharp_trailing_edge_and_plausible_geometry():
    c, x0, y0 = P["c"], P["x0"], P["y0"]
    zeta0, R = impl.airfoil_circle(c, x0, y0)
    assert zeta0 == complex(x0, y0) and R == pytest.approx(abs(c - zeta0))
    assert abs(impl.joukowski(np.array([c]), c)[0] - 2 * c) < 1e-15                    # the trailing edge
    zeta, z = impl.airfoil_coordinates(c, x0, y0, 600)
    assert z[0] == z[-1] == 2 * c
    npt.assert_allclose(np.abs(zeta - zeta0), R, rtol=1e-12)
    g = impl.airfoil_geometry(z)
    assert 4 * c < g["chord"] < 4.2 * c                                                  # a little over 4c
    assert 0.08 < g["thickness"] / g["chord"] < 0.16
    assert 0.02 < g["camber"] / g["chord"] < 0.07
    # No camber when y0 = 0: the airfoil is symmetric about the x axis.
    _, zs = impl.airfoil_coordinates(c, x0, 0.0, 600)
    assert impl.airfoil_geometry(zs)["camber"] == pytest.approx(0.0, abs=1e-9)
    npt.assert_allclose(np.sort(np.imag(zs)), np.sort(-np.imag(zs)), atol=1e-12)
    # Thickness grows with -x0.
    _, zt = impl.airfoil_coordinates(c, 2 * x0, y0, 600)
    assert impl.airfoil_geometry(zt)["thickness"] > 1.5 * g["thickness"]
    with pytest.raises(ValueError):
        impl.airfoil_circle(c, 0.01, 0.0)


def test_streamfunction_is_constant_on_the_cylinder_and_far_field_is_uniform():
    U, alpha, R, G = P["U"], P["alpha"], 0.3, 15.0
    zeta0 = -0.02 + 0.03j
    th = np.linspace(0, 2 * np.pi, 400)
    on = zeta0 + R * np.exp(1j * th)
    psi = np.imag(impl.cylinder_potential(on, U, alpha, R, G, zeta0))
    assert np.ptp(psi) < 1e-10 * U * R                                                   # the circle is a streamline
    # Normal velocity zero on the surface: Re(w * n) = 0 with n = e^{i th}.
    w = impl.cylinder_velocity(on, U, alpha, R, G, zeta0)
    npt.assert_allclose(np.real(w * np.exp(1j * th)), 0.0, atol=1e-10 * U)
    far = zeta0 + 1e4 * R * np.exp(1j * th[:10])
    npt.assert_allclose(impl.cylinder_velocity(far, U, alpha, R, G, zeta0), U * np.exp(-1j * alpha), rtol=1e-4)
    # Zero circulation, zero alpha: velocity doubles at the top of the cylinder.
    assert abs(impl.cylinder_velocity(np.array([zeta0 + 1j * R]), U, 0.0, R, 0.0, zeta0)[0]) == pytest.approx(2 * U)


def test_kutta_condition_puts_a_stagnation_point_at_the_trailing_edge():
    c, x0, y0, U, alpha = P["c"], P["x0"], P["y0"], P["U"], P["alpha"]
    zeta0, R = impl.airfoil_circle(c, x0, y0)
    beta = np.arctan2(y0, c - x0)
    for a in (-0.05, 0.0, alpha, 0.15):
        G = impl.kutta_circulation(U, R, a, beta)
        assert abs(impl.cylinder_velocity(np.array([c]), U, a, R, G, zeta0)[0]) < 1e-9 * U
        sp_ = impl.stagnation_points(U, a, R, G, zeta0)
        assert np.min(np.abs(sp_ - c)) < 1e-9
        npt.assert_allclose(np.abs(sp_ - zeta0), R, rtol=1e-9)                           # both on the circle
    # Without the Kutta circulation the trailing-edge speed on the airfoil is infinite; with it, finite.
    zeta, _ = impl.airfoil_coordinates(c, x0, y0, 400)
    G = impl.kutta_circulation(U, R, alpha, beta)
    w = impl.airfoil_surface_velocity(zeta, c, x0, y0, U, alpha, G)
    assert np.all(np.isfinite(w)) and abs(w[0]) < 2 * U
    assert abs(w[0]) == pytest.approx(abs(w[1]), rel=0.05)                               # continuous into the edge
    w0 = impl.airfoil_surface_velocity(zeta, c, x0, y0, U, alpha, 0.0)
    assert not np.isfinite(w0[0]) and abs(w0[1]) > 5 * U
    # Symmetric airfoil at zero incidence: no circulation, no lift.
    assert impl.kutta_circulation(U, R, 0.0, 0.0) == 0.0


def test_pressure_and_lift_are_consistent_with_kutta_joukowski():
    c, x0, y0, U, alpha, rho = P["c"], P["x0"], P["y0"], P["U"], P["alpha"], P["rho"]
    zeta0, R = impl.airfoil_circle(c, x0, y0)
    beta = np.arctan2(y0, c - x0)
    G = impl.kutta_circulation(U, R, alpha, beta)
    zeta, z = impl.airfoil_coordinates(c, x0, y0, 2000)
    w = impl.airfoil_surface_velocity(zeta, c, x0, y0, U, alpha, G)
    cp = impl.pressure_coefficient(w, U)
    assert cp.max() == pytest.approx(1.0, abs=2e-3)                                      # a stagnation point exists
    assert cp.min() < -1.0                                                               # suction peak on the upper surface
    npt.assert_allclose(impl.pressure_coefficient(np.array([0.0, U, 2 * U]), U), [1.0, 0.0, -3.0])
    # Integrate the surface pressure: F = i (1/2) rho U^2 closed-int C_p dz; lift = component normal to U.
    cp_mid = 0.5 * (cp[1:] + cp[:-1])
    force = 1j * 0.5 * rho * U ** 2 * np.sum(cp_mid * np.diff(z))
    lift = np.imag(force) * np.cos(alpha) - np.real(force) * np.sin(alpha)
    drag = np.real(force) * np.cos(alpha) + np.imag(force) * np.sin(alpha)
    L = impl.lift_per_span(rho, U, G)
    assert lift == pytest.approx(L, rel=2e-3)
    assert abs(drag) < 2e-3 * L                                                          # d'Alembert's paradox
    Cl = impl.lift_coefficient(G, U, impl.airfoil_geometry(z)["chord"])
    assert Cl == pytest.approx(L / (0.5 * rho * U ** 2 * impl.airfoil_geometry(z)["chord"]))
    # Lift slope close to 2 pi per radian (a little above, for thickness).
    G2 = impl.kutta_circulation(U, R, alpha + 0.01, beta)
    slope = (impl.lift_coefficient(G2, U, impl.airfoil_geometry(z)["chord"]) - Cl) / 0.01
    assert 2 * np.pi < slope < 1.2 * 2 * np.pi
    with pytest.raises(ValueError):
        impl.lift_coefficient(G, U, 0.0)


# ---------------------------------------------------------------- contour integrals

def test_contour_integrals_obey_the_residue_theorem():
    path = impl.circle_path(0.0, 1.0)
    assert path[0] == path[-1]
    assert impl.contour_integral(lambda z: 1 / z, path) == pytest.approx(2j * np.pi, abs=1e-10)
    assert abs(impl.contour_integral(lambda z: 1 / (z - 3.0), path)) < 1e-12               # pole outside
    assert abs(impl.contour_integral(lambda z: z ** 2 + np.exp(z), path)) < 1e-12           # analytic inside
    assert abs(impl.contour_integral(lambda z: 1 / z ** 2, path)) < 1e-12                   # residue zero
    # Sum of residues: 1/((z - 0.3)(z + 0.4)) has residues 1/0.7 and -1/0.7 inside, sum zero; 3/(z - 0.5i) has residue 3.
    assert abs(impl.contour_integral(lambda z: 1 / ((z - 0.3) * (z + 0.4)), path)) < 1e-10
    assert impl.contour_integral(lambda z: 3 / (z - 0.5j), path) == pytest.approx(6j * np.pi, abs=1e-10)
    # Path independence: a different circle around the same pole.
    assert impl.contour_integral(lambda z: 1 / z, impl.circle_path(0.2 + 0.1j, 0.5)) == pytest.approx(2j * np.pi, abs=1e-10)


def test_circulation_and_flux_from_the_velocity_contour_integral():
    U, alpha, R, G = P["U"], P["alpha"], 0.3, 20.0
    zeta0 = 0.1 - 0.05j
    w = lambda s: impl.cylinder_velocity(s, U, alpha, R, G, zeta0)
    for rad in (1.5 * R, 4.0 * R):
        gam, q = impl.circulation_and_flux(w, impl.circle_path(zeta0, rad))
        assert gam == pytest.approx(-G, rel=1e-9)                                        # clockwise Gamma, ccw path
        assert abs(q) < 1e-9 * U * R                                                     # no net outflow (no source)
    # A source adds flux and no circulation.
    src = lambda s: w(s) + 5.0 / (2 * np.pi * (s - zeta0))
    gam, q = impl.circulation_and_flux(src, impl.circle_path(zeta0, 2 * R))
    assert gam == pytest.approx(-G, rel=1e-9) and q == pytest.approx(5.0, rel=1e-9)
    # A contour not enclosing the cylinder sees nothing.
    gam, q = impl.circulation_and_flux(w, impl.circle_path(zeta0 + 3 * R, 0.5 * R))
    assert abs(gam) < 1e-9 * U * R and abs(q) < 1e-9 * U * R


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "joukowski_airfoil.png").stat().st_size > 1000
