"""Public tests for P10.03 The Kutta Condition & Circulation. They assert the physics of exact
potential flow: the circle is a streamline and the far field is the free stream; stagnation points
have zero speed and leave the body when Gamma > 4 pi a V; the Karman-Trefftz map reduces to Joukowski
at tau = 0, has the derivative it claims and inverts; the Kutta circulation makes the trailing-edge
point of the circle a stagnation point; and lift is rho V Gamma by four independent routes (surface
pressure, Blasius, momentum balance, and the circulation line integral around any enclosing contour),
with zero drag (d'Alembert). A contour that does not enclose the airfoil has zero circulation.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

B = 1.0
MU = complex(-0.08, 0.06)       # a cambered Joukowski/Karman-Trefftz circle centre
TAU = np.radians(10.0)


# ---------------------------------------------------------------- the circle

def test_circle_is_a_streamline_and_far_field_is_the_stream():
    V, al, a, G = 3.0, 0.2, 1.3, 5.0
    th = np.linspace(0, 2 * np.pi, 37)
    zeta = MU + a * np.exp(1j * th)
    w = impl.circle_dwdz(zeta, V, al, a, G, MU)
    ur = np.real(w * np.exp(1j * th))                        # radial velocity Re((u - iv) e^{i theta})
    npt.assert_allclose(ur, 0.0, atol=1e-12)
    ut = -np.imag(w * np.exp(1j * th))
    npt.assert_allclose(ut, impl.circle_surface_speed(th, V, al, a, G), atol=1e-12)
    far = impl.circle_dwdz(np.array([1e7 + 3e6j]), V, al, a, G, MU)[0]
    assert far == pytest.approx(V * np.exp(-1j * al), abs=1e-6)


def test_stagnation_points():
    V, a = 2.0, 1.0
    for G in (0.0, 3.0, 4 * np.pi * a * V * 0.9):
        st = impl.stagnation_angles(V, 0.1, a, G)
        assert st.size == 2
        npt.assert_allclose(impl.circle_surface_speed(st, V, 0.1, a, G), 0.0, atol=1e-12)
    assert impl.stagnation_angles(V, 0.0, a, 4 * np.pi * a * V * 1.01).size == 0


def test_cylinder_pressure_lift_is_rho_V_Gamma_and_drag_zero():
    a, V, n = 1.0, 1.0, 256
    th = 2 * np.pi * (np.arange(n) + 0.5) / n
    z = a * np.exp(1j * th)
    for al, G in ((0.0, 2.0), (0.3, -1.5), (-0.2, 9.0)):
        cl, cd = impl.pressure_force_coefficients(1j * z, impl.circle_surface_cp(th, V, al, a, G), al, 2 * a,
                                                  2 * np.pi / n)
        assert cl == pytest.approx(2 * G / (V * 2 * a), abs=1e-12)            # rho V Gamma / (q 2a)
        assert abs(cd) < 1e-12


# ---------------------------------------------------------------- the maps

def test_karman_trefftz_with_zero_angle_is_joukowski():
    zeta = np.array([2.0 + 1.0j, -1.5 + 0.3j, 0.2 - 2.0j])
    npt.assert_allclose(impl.karman_trefftz(zeta, B, 0.0), zeta + B ** 2 / zeta, atol=1e-12)


def test_map_derivative_and_inverse():
    zeta = MU + 1.4 * np.exp(1j * np.linspace(0.3, 6.0, 7))
    h = 1e-6
    fd = (impl.karman_trefftz(zeta + h, B, TAU) - impl.karman_trefftz(zeta - h, B, TAU)) / (2 * h)
    npt.assert_allclose(impl.karman_trefftz_dzdzeta(zeta, B, TAU), fd, rtol=1e-7)
    z = impl.karman_trefftz(zeta, B, TAU)
    npt.assert_allclose(impl.karman_trefftz_inverse(z, B, TAU, MU), zeta, atol=1e-10)


def test_circle_geometry_passes_through_the_critical_point():
    a, beta = impl.circle_geometry(B, MU)
    assert abs(MU + a * np.exp(-1j * beta) - B) < 1e-14
    assert beta > 0                                                    # centre above the axis: camber


# ---------------------------------------------------------------- the Kutta condition

def test_kutta_circulation_makes_the_trailing_edge_a_stagnation_point():
    a, beta = impl.circle_geometry(B, MU)
    for al in (-0.1, 0.0, 0.15):
        G = impl.kutta_circulation(1.0, a, al, beta)
        assert impl.circle_surface_speed(-beta, 1.0, al, a, G) == pytest.approx(0.0, abs=1e-12)
    assert impl.kutta_circulation(1.0, a, -beta, beta) == pytest.approx(0.0, abs=1e-15)   # zero-lift angle


def test_kutta_condition_keeps_the_trailing_edge_velocity_finite():
    al = np.radians(5.0)
    fk = impl.airfoil_flow(B, MU, 0.0, al, 1.0, 4096)
    f0 = impl.airfoil_flow(B, MU, 0.0, al, 1.0, 4096, Gamma=0.0)
    # Joukowski cusp: with Kutta the two surface points next to the TE have equal, finite speed
    assert fk["q"][0] == pytest.approx(fk["q"][-1], rel=1e-2) and fk["q"][0] < 2.0
    assert f0["q"][0] > 20 * fk["q"][0]                               # without it, flow turns the edge


# ---------------------------------------------------------------- lift four ways

@pytest.mark.parametrize("al_deg", [-3.0, 2.0, 8.0])
def test_surface_pressure_lift_equals_kutta_joukowski(al_deg):
    al = np.radians(al_deg)
    f = impl.airfoil_flow(B, MU, TAU, al, 1.0, 1024)
    cl, cd = impl.pressure_force_coefficients(f["dzdtheta"], f["cp"], al, f["chord"], 2 * np.pi / 1024)
    assert cl == pytest.approx(2 * f["Gamma"] / f["chord"], abs=1e-6)
    assert abs(cd) < 1e-6
    assert 2 * f["Gamma"] / f["chord"] == pytest.approx(8 * np.pi * f["a"] * np.sin(al + f["beta"]) / f["chord"])


def test_symmetric_section_at_zero_incidence():
    f = impl.airfoil_flow(B, complex(-0.1, 0.0), TAU, 0.0, 1.0, 512)
    assert abs(f["Gamma"]) < 1e-15
    npt.assert_allclose(f["cp"], f["cp"][::-1], atol=1e-12)            # upper and lower mirror images


@pytest.mark.parametrize("kind,size", [("circle", 6.0), ("ellipse", 3.2), ("rectangle", 2.8)])
def test_contour_integrals_give_Gamma_and_rho_V_Gamma(kind, size):
    al, rho, V = np.radians(4.0), 1.2, 1.0
    a, beta = impl.circle_geometry(B, MU)
    G = impl.kutta_circulation(V, a, al, beta)
    zc, dz = impl.contour_in_z(kind, size, center=0.2 + 0.1j, n=1200)
    w = impl.airfoil_velocity_at(zc, B, MU, TAU, al, V, G)
    assert impl.circulation_line_integral(w, dz) == pytest.approx(G, rel=1e-9)
    Lb, Db = impl.blasius_force(w, dz, rho, al)
    Lm, Dm = impl.momentum_force(w, dz, rho, V, al)
    assert Lb == pytest.approx(rho * V * G, rel=1e-9) and abs(Db) < 1e-9
    assert Lm == pytest.approx(rho * V * G, rel=1e-9) and abs(Dm) < 1e-9


def test_contour_not_enclosing_the_airfoil_has_no_circulation_or_force():
    al = np.radians(4.0)
    a, beta = impl.circle_geometry(B, MU)
    G = impl.kutta_circulation(1.0, a, al, beta)
    zc, dz = impl.contour_in_z("circle", 1.0, center=6.0 + 2.0j, n=600)
    w = impl.airfoil_velocity_at(zc, B, MU, TAU, al, 1.0, G)
    assert abs(impl.circulation_line_integral(w, dz)) < 1e-10
    assert abs(impl.momentum_force(w, dz, 1.0, 1.0, al)[0]) < 1e-10


def test_thickness_raises_the_lift_slope_above_2_pi():
    slopes = []
    for t in (0.04, 0.12, 0.20):
        mu = impl.fit_karman_trefftz(t, 0.0, 0.0)
        sh = impl.section_shape(B, mu, 0.0)
        a, _ = impl.circle_geometry(B, mu)
        slopes.append(8 * np.pi * a / sh["chord"] / (2 * np.pi))
        assert slopes[-1] == pytest.approx(1 + 0.77 * t, abs=2e-3)
    assert np.all(np.diff(slopes) > 0)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "kutta_circulation.png").stat().st_size > 1000
