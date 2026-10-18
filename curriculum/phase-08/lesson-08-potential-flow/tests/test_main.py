"""Public tests for P08.08 Potential Flow. They assert the physics of plane potential flow, not magic
numbers: every element's analytic velocity is the gradient of its potential and satisfies the
Cauchy-Riemann equations; potentials are harmonic; a source has flux m and no circulation, a
vortex the reverse; a doublet is the limit of a source-sink pair; superposition is linear; the
cylinder and the traced Rankine oval have zero normal velocity on their surfaces; the stagnation
points are where the velocity vanishes; the force from integrating the surface pressure equals the
force from Blasius's theorem on any enclosing contour, with lift rho U Gamma (Kutta-Joukowski) and
zero drag (d'Alembert); and the Kutta condition on a Joukowski flat plate gives c_l = 2 pi sin(alpha).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

U, A = 3.0, 0.5
PTS = np.array([1.3 + 0.7j, -0.9 + 1.1j, 0.4 - 1.6j, -2.0 - 0.3j])


def all_kinds():
    return [impl.uniform(U, 0.3), impl.source(2.0, 0.1 + 0.2j), impl.doublet(1.5, -0.2j, 0.7),
            impl.vortex(4.0, 0.3 - 0.1j)]


# ---------------------------------------------------------------- elements

def test_velocity_is_the_gradient_of_the_potential_and_satisfies_cauchy_riemann():
    h = 1e-6
    for e in all_kinds():
        W = lambda z: impl.complex_potential([e], z)
        dW = impl.complex_velocity([e], PTS)
        # complex derivative along x and along y must agree (analyticity)
        npt.assert_allclose((W(PTS + h) - W(PTS - h)) / (2 * h), dW, rtol=1e-7, atol=1e-9)
        npt.assert_allclose((W(PTS + 1j * h) - W(PTS - 1j * h)) / (2j * h), dW, rtol=1e-7, atol=1e-9)
        u, v = impl.velocity([e], PTS)
        phi = lambda z: W(z).real
        psi = lambda z: W(z).imag
        npt.assert_allclose(u, (phi(PTS + h) - phi(PTS - h)) / (2 * h), rtol=1e-6, atol=1e-8)     # u = phi_x
        npt.assert_allclose(u, (psi(PTS + 1j * h) - psi(PTS - 1j * h)) / (2 * h), rtol=1e-6, atol=1e-8)  # u = psi_y
        npt.assert_allclose(v, -(psi(PTS + h) - psi(PTS - h)) / (2 * h), rtol=1e-6, atol=1e-8)    # v = -psi_x


def test_potentials_are_harmonic():
    h = 1e-3
    for e in all_kinds():
        W = lambda z: impl.complex_potential([e], z)
        lap = (W(PTS + h) + W(PTS - h) + W(PTS + 1j * h) + W(PTS - 1j * h) - 4 * W(PTS)) / h ** 2
        assert np.max(np.abs(lap)) < 1e-4                  # Laplace's equation for phi and psi


def test_source_flux_and_vortex_circulation():
    c = impl.circle_contour(1.0, 64, 0.2 + 0.1j)
    circ, flux = impl.circulation_and_flux([impl.source(2.5, 0.3)], c)
    assert flux == pytest.approx(2.5, rel=1e-12) and abs(circ) < 1e-12
    circ, flux = impl.circulation_and_flux([impl.vortex(4.0, 0.3)], c)
    assert circ == pytest.approx(-4.0, rel=1e-12) and abs(flux) < 1e-12   # clockwise Gamma: ccw integral is -Gamma
    # a contour that does not enclose the singularity sees neither
    far = impl.circle_contour(0.5, 64, 3.0 + 0j)
    for e in (impl.source(2.5, 0.3), impl.vortex(4.0, 0.3)):
        circ, flux = impl.circulation_and_flux([e], far)
        assert abs(circ) < 1e-12 and abs(flux) < 1e-12
    # the clockwise sense: on the right of a positive vortex the flow goes down
    u, v = impl.velocity([impl.vortex(1.0)], np.array([1.0 + 0j]))
    assert abs(u[0]) < 1e-15 and v[0] == pytest.approx(-1 / (2 * np.pi))


def test_doublet_is_the_limit_of_a_source_sink_pair():
    kappa, beta = 1.2, 0.4
    z = np.array([1.5 + 0.5j, -1.0 + 2.0j])
    for eps, tol in ((1e-2, 2e-4), (1e-3, 2e-6)):
        m = kappa / (2 * eps)
        d = eps * np.exp(1j * beta)
        pair = [impl.source(m, -d), impl.source(-m, d)]
        npt.assert_allclose(impl.complex_velocity(pair, z), impl.complex_velocity([impl.doublet(kappa, 0j, beta)], z),
                            rtol=tol)


def test_superposition_is_linear():
    els = all_kinds()
    total = impl.complex_velocity(els, PTS)
    parts = sum(impl.complex_velocity([e], PTS) for e in els)
    npt.assert_allclose(total, parts, rtol=1e-14)
    assert np.all(impl.complex_velocity([], PTS) == 0)


# ---------------------------------------------------------------- cylinder

@pytest.mark.parametrize("Gamma,alpha", [(0.0, 0.0), (5.0, 0.0), (-3.0, 0.25)])
def test_no_flow_through_the_cylinder_and_uniform_flow_far_away(Gamma, alpha):
    th = np.linspace(0, 2 * np.pi, 37)
    zs = A * np.exp(1j * th)
    u, v = impl.velocity(impl.cylinder(U, A, Gamma, alpha), zs)
    u_r = u * np.cos(th) + v * np.sin(th)
    assert np.max(np.abs(u_r)) < 1e-12 * U
    far = impl.complex_velocity(impl.cylinder(U, A, Gamma, alpha), np.array([1e6 + 3e5j]))
    assert far[0] == pytest.approx(U * np.exp(-1j * alpha), rel=1e-5)


def test_surface_pressure_closed_form_and_stagnation_points():
    Gamma = 2 * np.pi * U * A * 0.8
    th = np.linspace(0, 2 * np.pi, 73)
    cp = impl.pressure_coefficient(impl.cylinder(U, A, Gamma), A * np.exp(1j * th), U)
    npt.assert_allclose(cp, impl.cylinder_surface_cp(th, U, A, Gamma), atol=1e-12)
    cp0 = impl.cylinder_surface_cp(th, U, A, 0.0)
    assert cp0.min() == pytest.approx(-3.0) and cp0.max() == pytest.approx(1.0)     # 2U at the shoulders
    for G in (0.0, Gamma, 4 * np.pi * U * A, 6 * np.pi * U * A):
        zs = impl.cylinder_stagnation_points(U, A, G)
        assert np.max(np.abs(impl.complex_velocity(impl.cylinder(U, A, G), zs))) < 1e-10 * U
        assert np.all(zs.imag <= 1e-12)                    # clockwise circulation moves them down
    assert len(impl.cylinder_stagnation_points(U, A, 6 * np.pi * U * A)) == 1
    assert abs(impl.cylinder_stagnation_points(U, A, 6 * np.pi * U * A)[0]) > A   # off the body
    assert impl.cylinder_stagnation_points(U, A, 4 * np.pi * U * A) == pytest.approx(-1j * A)  # merged at the bottom


def test_kutta_joukowski_and_dalembert_by_two_routes():
    rho = 1.2
    for Gamma in (0.0, 1.7, 2 * np.pi * U * A):
        els = impl.cylinder(U, A, Gamma)
        Fx, Fy = impl.force_by_pressure(els, impl.circle_contour(A, 96), rho, U, p_inf=1e5)
        assert abs(Fx) < 1e-9 * rho * U ** 2 * A                         # no drag
        assert Fy == pytest.approx(rho * U * Gamma, abs=1e-9 * rho * U ** 2 * A)   # lift rho U Gamma, up
        for contour in (impl.circle_contour(1.5 * A, 96), impl.circle_contour(4 * A, 96, 0.3 - 0.2j)):
            Bx, By = impl.force_blasius(els, contour, rho)
            assert Bx == pytest.approx(Fx, abs=1e-9) and By == pytest.approx(Fy, abs=1e-9)


def test_cp_does_not_depend_on_the_speed_for_fixed_gamma_over_ua():
    th = np.linspace(0, 2 * np.pi, 19)
    g = 0.6
    cp1 = impl.cylinder_surface_cp(th, 1.0, A, g * 1.0 * A)
    cp2 = impl.cylinder_surface_cp(th, 50.0, A, g * 50.0 * A)
    npt.assert_allclose(cp1, cp2, atol=1e-13)


# ---------------------------------------------------------------- Rankine oval

def test_rankine_oval_dimensions_and_stagnation_points():
    Uo, m, a = 2.0, 5.0, 0.8
    l, t = impl.rankine_oval_size(Uo, m, a)
    els = impl.rankine_oval(Uo, m, a)
    assert np.max(np.abs(impl.complex_velocity(els, np.array([l, -l])))) < 1e-12 * Uo
    assert abs(impl.stream_function_oval(Uo, m, a, np.array([1j * t]))[0]) < 1e-12
    # a thin oval (weak source) tends to the two points +-a; a strong one grows like a cylinder
    l_small, t_small = impl.rankine_oval_size(Uo, 1e-6, a)
    assert l_small == pytest.approx(a, rel=1e-6) and t_small < 1e-6
    assert l > a and t < l


def test_traced_oval_has_no_normal_flow_and_no_force():
    Uo, m, a = 1.0, 2 * np.pi, 1.0
    cont = impl.oval_contour(Uo, m, a, 128)
    els = impl.rankine_oval(Uo, m, a)
    assert np.max(np.abs(impl.stream_function_oval(Uo, m, a, cont))) < 1e-10
    # normal velocity from the geometric normal of the traced points (FFT tangent): u dy - v dx = 0
    zc = np.fft.ifft(1j * np.fft.fftfreq(cont.size, 1 / cont.size) * np.fft.fft(cont))
    u_n = (impl.complex_velocity(els, cont) * zc).imag / np.abs(zc)
    assert np.max(np.abs(u_n)) < 1e-8 * Uo
    Fx, Fy = impl.force_by_pressure(els, cont, 1.0, Uo)
    assert abs(Fx) < 1e-9 and abs(Fy) < 1e-9                             # d'Alembert
    # net source strength zero: no flux out of the body
    assert abs(impl.circulation_and_flux(els, impl.circle_contour(5.0, 128))[1]) < 1e-12


# ---------------------------------------------------------------- Joukowski and the Kutta condition

def test_joukowski_maps_the_circle_to_a_plate_and_inverts():
    a = 0.375
    zc = impl.joukowski(impl.circle_contour(a, 64), a)
    assert np.max(np.abs(zc.imag)) < 1e-14 and zc.real.max() == pytest.approx(2 * a) and zc.real.min() == pytest.approx(-2 * a)
    z = np.array([3.0 + 0.1j, -1.0 - 2.0j, 0.1 + 0.05j, 5j])
    zeta = impl.zeta_from_z(z, a)
    assert np.all(np.abs(zeta) >= a)
    npt.assert_allclose(impl.joukowski(zeta, a), z, rtol=1e-13)


def test_kutta_condition_and_thin_plate_lift():
    Uo, a, rho = 60.0, 0.375, 0.9
    for alpha in np.radians([1.0, 4.0, 8.0]):
        G = impl.kutta_circulation(Uo, a, alpha)
        # rear stagnation point exactly at the image of the trailing edge
        assert abs(impl.complex_velocity(impl.cylinder(Uo, a, G, alpha), np.array([a + 0j]))[0]) < 1e-10 * Uo
        # in the plate plane the trailing-edge velocity stays finite; the leading edge's does not
        te = impl.plate_complex_velocity(np.array([2 * a + 1e-6j, 2 * a + 1e-8j]), Uo, a, alpha, G)
        le = impl.plate_complex_velocity(np.array([-2 * a + 1e-6j, -2 * a + 1e-8j]), Uo, a, alpha, G)
        assert abs(te[1]) < 2 * Uo
        assert abs(le[1]) > 5 * abs(le[0])                # ~ d^(-1/2): 100 times closer, about 10 times faster
        # Blasius on a big circle round the plate: force normal to the free stream, magnitude rho U Gamma
        c = impl.circle_contour(6 * a, 256)
        I = 0.5j * rho * impl.contour_integral(impl.plate_complex_velocity(c, Uo, a, alpha, G) ** 2, c)
        Fx, Fy = I.real, -I.imag
        lift = -Fx * np.sin(alpha) + Fy * np.cos(alpha)
        drag = Fx * np.cos(alpha) + Fy * np.sin(alpha)
        assert lift == pytest.approx(rho * Uo * G, rel=1e-9) and abs(drag) < 1e-9 * lift
        assert lift / (0.5 * rho * Uo ** 2 * 4 * a) == pytest.approx(2 * np.pi * np.sin(alpha), rel=1e-9)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "potential_flow.png").stat().st_size > 10_000
