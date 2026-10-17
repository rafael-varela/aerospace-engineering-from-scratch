"""Public tests for P05.08 Torsion of Shafts. They assert physics, not magic numbers: torque times
speed is power; tau = T r / J grows linearly in r and integrates back to T (equilibrium of the
section); sizing formulas return shafts that are exactly at the allowable, and the hollow formula
tends to the solid one; stepped twists add and the indeterminate shaft obeys equilibrium and
compatibility; Bredt-Batho has constant shear flow whose moment is T, reproduces the exact thin
tube as t -> 0, and is independent of the polygon's starting vertex and orientation; the open
section is (r/t)^2 times softer. Independent check: Prandtl's stress function by finite
differences converges at second order to the exact series for a solid rectangle, tends to b t^3/3
for a thin strip, and approaches Bredt-Batho for a thin square tube.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

G_TI = 113.8e3 / 2.684


# ---------------------------------------------------------------- circular shafts

def test_power_speed_torque():
    T = impl.shaft_torque(1.4e6, 20_900)                       # N mm
    omega = 20_900 * 2 * np.pi / 60
    assert T / 1e3 * omega == pytest.approx(1.4e6)            # P = T omega
    assert impl.shaft_torque(1.4e6, 2 * 20_900) == pytest.approx(T / 2)


def test_polar_moment_and_stress_integrate_back_to_torque():
    d_o, d_i, T = 38.0, 30.0, 8e5
    J = impl.polar_moment(d_o, d_i)
    assert J == pytest.approx(quad(lambda r: 2 * np.pi * r ** 3, d_i / 2, d_o / 2)[0], rel=1e-12)
    # equilibrium of the cut section: integral of tau r dA = T
    Tint = quad(lambda r: impl.torsion_shear_stress(T, r, J) * r * 2 * np.pi * r, d_i / 2, d_o / 2)[0]
    assert Tint == pytest.approx(T, rel=1e-12)
    assert impl.torsion_shear_stress(T, 0.0, J) == 0.0
    assert impl.torsion_shear_stress(T, d_o / 2, J) == pytest.approx(2 * impl.torsion_shear_stress(T, d_o / 4, J))
    assert impl.polar_moment(d_o) == pytest.approx(np.pi * d_o ** 4 / 32)


def test_twist_is_linear_and_stepped_twists_add():
    T, L, d = 8e5, 600.0, 35.0
    J = impl.polar_moment(d)
    phi = impl.twist_angle(T, L, G_TI, J)
    assert impl.twist_angle(2 * T, L, G_TI, J) == pytest.approx(2 * phi)
    assert impl.stepped_shaft_twist([T, T], [200.0, 400.0], [G_TI, G_TI], [J, J]) == pytest.approx(phi)
    # kinematics: surface shear strain gamma = r phi / L equals tau / G (Hooke's law in shear)
    assert (d / 2) * phi / L == pytest.approx(impl.torsion_shear_stress(T, d / 2, J) / G_TI)


def test_indeterminate_shaft_equilibrium_and_compatibility():
    T0 = 8e5
    GJ1, GJ2 = G_TI * impl.polar_moment(24.0), G_TI * impl.polar_moment(30.0, 24.0)
    phi, T1, T2 = impl.fixed_fixed_shaft(T0, 200.0, GJ1, 400.0, GJ2)
    assert T1 - T2 == pytest.approx(T0)                            # junction equilibrium
    assert T1 * 200.0 / GJ1 + T2 * 400.0 / GJ2 == pytest.approx(0.0, abs=1e-12)   # zero total twist
    assert phi == pytest.approx(T1 * 200.0 / GJ1)
    # infinitely stiff segment 2 takes all the torque
    _, T1r, T2r = impl.fixed_fixed_shaft(T0, 200.0, GJ1, 400.0, 1e12 * GJ2)
    assert T1r == pytest.approx(0.0, abs=1e-6 * T0) and T2r == pytest.approx(-T0, rel=1e-6)


def test_sizing_hits_the_allowable_and_hollow_tends_to_solid():
    T, tau = 1.2e6, 570.0
    d = impl.size_solid_shaft(T, tau)
    assert impl.torsion_shear_stress(T, d / 2, impl.polar_moment(d)) == pytest.approx(tau)
    assert impl.size_hollow_shaft(T, tau, 0.0) == pytest.approx(d)
    for ratio in (0.5, 0.8, 0.9):
        do = impl.size_hollow_shaft(T, tau, ratio)
        assert do > d
        assert impl.torsion_shear_stress(T, do / 2, impl.polar_moment(do, ratio * do)) == pytest.approx(tau)
        # same stress, less material
        assert impl.shaft_mass_per_length(do, ratio * do, 4430.0) < impl.shaft_mass_per_length(d, 0.0, 4430.0)
    dt = impl.size_for_twist(T, 600.0, G_TI, np.radians(5.0), 0.8)
    assert impl.twist_angle(T, 600.0, G_TI, impl.polar_moment(dt, 0.8 * dt)) == pytest.approx(np.radians(5.0))


def test_shear_allowables_and_material_table():
    for name, m in impl.MATERIALS.items():
        Fsy, Fsu = impl.shear_allowables(name)
        assert Fsy == pytest.approx(m["Fty"] / np.sqrt(3)) and Fsu == pytest.approx(0.6 * m["Ftu"])
        nu = m["E"] / (2 * m["G"]) - 1
        assert 0.25 < nu < 0.36, name                              # G consistent with isotropy


# ---------------------------------------------------------------- thin-walled sections

def test_bredt_batho_constant_shear_flow_and_moment_equals_torque():
    T = 10e6
    v = impl.WING_BOX
    t = impl.WING_BOX_T
    r = impl.bredt_batho(T, v, t, 27.5e3)
    npt.assert_allclose(r["tau"] * t, r["q"])                      # q = tau t on every wall
    # moment of the shear flow about an arbitrary point = T
    p = np.array([137.0, -55.0])
    M = 0.0
    for k in range(len(v)):
        a, b = v[k], v[(k + 1) % len(v)]
        F = r["q"] * (b - a)                                       # force along the wall
        arm = a - p
        M += arm[0] * F[1] - arm[1] * F[0]
    assert abs(M) == pytest.approx(T, rel=1e-12)
    assert r["A_m"] == pytest.approx(400 * 120)
    assert r["twist_rate"] == pytest.approx(T / (27.5e3 * r["J"]))


def test_bredt_batho_is_independent_of_vertex_order_and_start():
    v, t = impl.WING_BOX, impl.WING_BOX_T
    a = impl.bredt_batho(1e6, v, t, 27.5e3)
    b = impl.bredt_batho(1e6, np.roll(v, 2, axis=0), np.roll(t, 2), 27.5e3)
    vr = v[::-1]
    tr = np.roll(t[::-1], -1)                                      # walls follow their vertices
    c = impl.bredt_batho(1e6, vr, tr, 27.5e3)
    for k in ("J", "q", "A_m", "twist_rate"):
        assert a[k] == pytest.approx(b[k]) and a[k] == pytest.approx(c[k])
    npt.assert_allclose(np.sort(a["tau"]), np.sort(c["tau"]))


def test_bredt_batho_thin_circular_tube_tends_to_exact():
    rm = 20.0
    ratios = []
    for t in (4.0, 2.0, 1.0, 0.5):
        exact = impl.polar_moment(2 * (rm + t / 2), 2 * (rm - t / 2))
        n = 4000
        ang = np.linspace(0, 2 * np.pi, n, endpoint=False)
        poly = impl.bredt_batho(1.0, rm * np.column_stack([np.cos(ang), np.sin(ang)]), t, 1.0)["J"]
        assert poly == pytest.approx(impl.thin_tube_J(rm, t), rel=1e-5)   # polygon -> circle
        ratios.append(impl.thin_tube_J(rm, t) / exact)
    # exact ratio is 1 / (1 + (t / 2 r_m)^2): the error falls as t^2
    assert np.all(np.diff(ratios) > 0) and ratios[-1] == pytest.approx(1 / (1 + (0.5 / 40) ** 2), rel=1e-12)


def test_open_section_is_much_softer_than_closed():
    rm, t = 20.0, 1.0
    J_closed = impl.thin_tube_J(rm, t)
    J_open = impl.open_section_J([2 * np.pi * rm], [t])            # the slit tube
    assert J_closed / J_open == pytest.approx(3 * (rm / t) ** 2)
    r = impl.bredt_batho(1e6, impl.WING_BOX, impl.WING_BOX_T, 27.5e3)
    assert r["J"] / impl.open_section_J(r["lengths"], impl.WING_BOX_T) > 1000


# ---------------------------------------------------------------- independent check: Prandtl FD

def test_rectangle_series_limits():
    # square: J = 0.1406 a^4 (classical coefficient); thin strip: b t^3 / 3 (1 - 0.630 t / b)
    assert impl.rectangle_torsion_constant(1.0, 1.0) == pytest.approx(0.1406, abs=1e-4)
    a, b = 20.0, 1.0
    assert impl.rectangle_torsion_constant(a, b) == pytest.approx(a * b ** 3 / 3 * (1 - 0.6302 * b / a), rel=1e-4)
    assert impl.rectangle_torsion_constant(1.0, 2.0) == impl.rectangle_torsion_constant(2.0, 1.0)


def test_prandtl_fd_converges_at_second_order_to_series():
    exact = impl.rectangle_torsion_constant(2.0, 1.0)
    errs = [abs(impl.prandtl_fd(2.0, 1.0, h)[0] - exact) for h in (0.1, 0.05, 0.025)]
    orders = np.log2(np.array(errs[:-1]) / np.array(errs[1:]))
    assert np.all(orders > 1.9) and np.all(orders < 2.1)
    J, phi = impl.prandtl_fd(2.0, 1.0, 0.05)
    assert np.all(phi >= 0) and phi[0, :].max() == 0 and phi[:, -1].max() == 0   # zero on the boundary
    npt.assert_allclose(phi, phi[::-1, ::-1], atol=1e-12)                       # symmetric


def test_prandtl_fd_hollow_square_approaches_bredt_as_wall_thins():
    ratios = []
    for t in (0.2, 0.1, 0.05):
        J, phi = impl.prandtl_fd(1.0, 1.0, t / 8, hole=(1 - 2 * t, 1 - 2 * t))
        jb = impl.bredt_batho(1.0, np.array([[0, 0], [1 - t, 0], [1 - t, 1 - t], [0, 1 - t]]), t, 1.0)["J"]
        ratios.append(jb / J)
        # the hole carries one constant value: the stress function is flat across it
        c = phi.shape[0] // 2
        assert phi[c, c] == pytest.approx(phi[c, c + 2]) and phi[c, c] == phi.max()
    assert np.all(np.diff(ratios) > 0)                   # Bredt-Batho improves as the wall thins
    assert 0.96 < ratios[-1] < 1.0                       # within 4 % at t / a = 0.05 (corners add stiffness)
    # removing material softens: the tube is less stiff than the solid square
    assert impl.prandtl_fd(1.0, 1.0, 0.025, hole=(0.8, 0.8))[0] < impl.prandtl_fd(1.0, 1.0, 0.025)[0]


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "torsion.png").stat().st_size > 10_000
