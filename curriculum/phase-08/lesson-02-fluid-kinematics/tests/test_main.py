"""Public tests for P08.02 Fluid Kinematics. They assert kinematics, not magic numbers: RK4
pathlines and streaklines of an oscillating stream reproduce the closed forms and converge at
fourth order; in a steady flow streamlines, pathlines and streaklines coincide, and in an unsteady
one they do not; a decaying Taylor-Green vortex is unsteady yet its lines coincide; the
velocity-gradient tensor splits into strain and rotation with the right divergence and vorticity
for stagnation flow, shear and rigid rotation; the material derivative equals the acceleration of
a computed pathline and vanishes for a quantity carried by the flow; an incompressible flow keeps
the area of a material loop.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

U, V, OMEGA = 10.0, 0.875, 4 * np.pi           # a 10 m/s stream swinging +/- 5 deg at 2 Hz


# ---------------------------------------------------------------- the oscillating stream: closed forms

def test_pathline_matches_closed_form_and_converges_at_fourth_order():
    vel = impl.oscillating_flow(U, V, OMEGA)
    t, x, y = impl.pathline(vel, 0.1, -0.2, 0.05, 0.6, 200)
    xe, ye = impl.oscillating_pathline_exact(U, V, OMEGA, 0.1, -0.2, 0.05, t)
    npt.assert_allclose(x, xe, atol=1e-12)
    npt.assert_allclose(y, ye, atol=1e-9)
    errs = [abs(impl.pathline(vel, 0.0, 0.0, 0.0, 0.3, n)[2][-1] - impl.oscillating_pathline_exact(U, V, OMEGA, 0, 0, 0, 0.3)[1])
            for n in (10, 20, 40)]
    orders = impl.observed_order(errs)
    assert np.all(orders > 3.8) and np.all(orders < 4.2)


def test_closed_form_pathline_satisfies_the_ode():
    t = np.linspace(0.0, 1.0, 11)
    h = 1e-6
    x1, y1 = impl.oscillating_pathline_exact(U, V, OMEGA, 0.0, 0.0, 0.0, t + h)
    x0, y0 = impl.oscillating_pathline_exact(U, V, OMEGA, 0.0, 0.0, 0.0, t - h)
    npt.assert_allclose((x1 - x0) / (2 * h), U, rtol=1e-8)
    npt.assert_allclose((y1 - y0) / (2 * h), V * np.cos(OMEGA * t), atol=1e-6)


def test_streakline_matches_closed_form():
    vel = impl.oscillating_flow(U, V, OMEGA)
    tau, x, y = impl.streakline(vel, 0.0, 0.0, 0.0, 1.0, 100)
    assert tau.shape == x.shape == y.shape == (101,)
    assert x[-1] == 0.0 and y[-1] == 0.0                                  # newest particle at the source
    xe, ye = impl.oscillating_streakline_exact(U, V, OMEGA, 0.0, 0.0, 1.0, tau)
    npt.assert_allclose(x, xe, atol=1e-10)
    npt.assert_allclose(y, ye, atol=1e-8)


def test_streamline_of_the_oscillating_stream_is_straight_at_each_instant():
    vel = impl.oscillating_flow(U, V, OMEGA)
    for t in (0.0, 0.1, 0.2):
        x, y = impl.streamline(vel, 0.0, 0.0, t, 2.0, 50)
        npt.assert_allclose(y, V * np.cos(OMEGA * t) / U * x, atol=1e-12)
        assert np.hypot(x[-1], y[-1]) == pytest.approx(2.0, rel=1e-12)    # arc length 2 m


def test_unsteady_lines_differ():
    vel = impl.oscillating_flow(U, V, OMEGA)
    # (not t = 0.5 s: after a whole period this pathline and streakline happen to trace the same curve)
    tau, xs, ys = impl.streakline(vel, 0.0, 0.0, 0.0, 0.3, 100)
    xl, yl = impl.streamline(vel, 0.0, 0.0, 0.3, 3.0, 300)
    assert np.max(impl.distance_to_curve(xs, ys, xl, yl)) > 0.05             # tens of mm apart
    _, xp, yp = impl.pathline(vel, 0.0, 0.0, 0.0, 0.3, 100)
    assert np.max(impl.distance_to_curve(xp, yp, xs, ys)) > 0.05


def test_smoke_error_vanishes_as_the_reduced_frequency_goes_to_zero():
    e = [max(impl.smoke_error(U, V, OMEGA, L, t) for t in np.linspace(0, 0.5, 41)) for L in (0.01, 0.02, 0.04)]
    k = OMEGA * np.array([0.01, 0.02, 0.04]) / U
    npt.assert_allclose(np.array(e) / k, 0.5, rtol=0.02)                     # error ~ k/2 for small k


# ---------------------------------------------------------------- steady flow: the lines coincide

def test_steady_stagnation_flow_lines_coincide_on_a_hyperbola():
    vel = impl.stagnation_flow(1.5)
    _, xs, ys = impl.streakline(vel, 0.2, 1.0, 0.0, 1.5, 150)
    _, xp, yp = impl.pathline(vel, 0.2, 1.0, 0.0, 1.5, 150)
    xl, yl = impl.streamline(vel, 0.2, 1.0, 3.7, 2.0, 400)                   # any frozen time
    for x, y in ((xs, ys), (xp, yp), (xl, yl)):
        npt.assert_allclose(x * y, 0.2, rtol=1e-8)
    npt.assert_allclose(np.sort(xs), np.sort(xp), rtol=1e-9)                # same particles, same places


def test_taylor_green_unsteady_but_lines_coincide():
    vel = impl.taylor_green(nu=0.05)
    psi = lambda x, y: np.sin(x) * np.sin(y)                               # streamfunction shape
    _, xp, yp = impl.pathline(vel, 1.0, 0.7, 0.0, 3.0, 300)
    xl, yl = impl.streamline(vel, 1.0, 0.7, 2.0, 1.5, 300)
    npt.assert_allclose(psi(xp, yp), psi(1.0, 0.7), atol=1e-9)
    npt.assert_allclose(psi(xl, yl), psi(1.0, 0.7), atol=1e-9)


def test_streamline_depends_on_direction_only():
    slow, fast = impl.stagnation_flow(0.5), impl.stagnation_flow(5.0)
    x1, y1 = impl.streamline(slow, 0.3, 1.0, 0.0, 1.0, 100)
    x2, y2 = impl.streamline(fast, 0.3, 1.0, 0.0, 1.0, 100)
    npt.assert_allclose(x1, x2, rtol=1e-12)
    npt.assert_allclose(y1, y2, rtol=1e-12)


def test_pathline_converges_against_an_independent_integrator():
    vel = impl.pulsating_stagnation_flow(1.0, 0.5, 3.0)
    sol = solve_ivp(lambda t, z: np.array(vel(z[0], z[1], t), dtype=float).ravel(), (0.0, 2.0), [0.3, 1.0],
                    rtol=1e-12, atol=1e-14)
    errs = []
    for n in (20, 40, 80):
        _, x, y = impl.pathline(vel, 0.3, 1.0, 0.0, 2.0, n)
        errs.append(np.hypot(x[-1] - sol.y[0, -1], y[-1] - sol.y[1, -1]))
    assert errs[-1] < 1e-7
    assert np.all(impl.observed_order(errs) > 3.7)


# ---------------------------------------------------------------- strain, rotation, divergence

def test_stagnation_flow_is_pure_strain():
    d = impl.decompose_gradient(impl.velocity_gradient(impl.stagnation_flow(2.0), 0.4, -0.3, 0.0))
    npt.assert_allclose(d["S"], [[2.0, 0.0], [0.0, -2.0]], atol=1e-8)
    npt.assert_allclose(d["W"], 0.0, atol=1e-8)
    assert abs(d["divergence"]) < 1e-8 and abs(d["vorticity"]) < 1e-8


def test_shear_is_half_strain_half_rotation():
    g = 3.0
    d = impl.decompose_gradient(impl.velocity_gradient(impl.shear_flow(g), 0.1, 0.2, 0.0))
    npt.assert_allclose(d["S"], [[0.0, g / 2], [g / 2, 0.0]], atol=1e-8)
    npt.assert_allclose(d["W"], [[0.0, g / 2], [-g / 2, 0.0]], atol=1e-8)
    assert d["vorticity"] == pytest.approx(-g, abs=1e-8)
    npt.assert_allclose(d["principal_strain_rates"], [-g / 2, g / 2], atol=1e-8)   # extension along 45 deg


def test_rigid_rotation_has_no_strain_and_vorticity_twice_omega():
    d = impl.decompose_gradient(impl.velocity_gradient(impl.solid_body_rotation(0.7), 1.2, -0.5, 0.0))
    npt.assert_allclose(d["S"], 0.0, atol=1e-8)
    assert d["vorticity"] == pytest.approx(1.4, abs=1e-8)


def test_example_fields_are_divergence_free():
    fields = [impl.oscillating_flow(U, V, OMEGA), impl.taylor_green(0.1), impl.pulsating_stagnation_flow(1.0, 0.5, 3.0)]
    for fld in fields:
        for (x, y, t) in ((0.3, 0.4, 0.1), (1.1, -0.6, 0.7)):
            assert abs(impl.decompose_gradient(impl.velocity_gradient(fld, x, y, t))["divergence"]) < 1e-8


# ---------------------------------------------------------------- the material derivative

def test_particle_acceleration_equals_pathline_second_derivative():
    vel = impl.pulsating_stagnation_flow(1.0, 0.5, 3.0)
    t, x, y = impl.pathline(vel, 0.3, 1.0, 0.0, 1.0, 1000)
    i, dt = 500, t[1] - t[0]
    a_path = np.array([(x[i + 1] - 2 * x[i] + x[i - 1]) / dt ** 2, (y[i + 1] - 2 * y[i] + y[i - 1]) / dt ** 2])
    npt.assert_allclose(impl.particle_acceleration(vel, x[i], y[i], t[i]), a_path, rtol=1e-5)


def test_steady_flow_has_convective_acceleration_only():
    a = 1.5
    acc = impl.particle_acceleration(impl.stagnation_flow(a), 2.0, 1.0, 0.0)
    npt.assert_allclose(acc, [a ** 2 * 2.0, a ** 2 * 1.0], rtol=1e-7)
    acc = impl.particle_acceleration(impl.oscillating_flow(U, V, OMEGA), 5.0, 1.0, 0.1)   # uniform: local only
    npt.assert_allclose(acc, [0.0, -V * OMEGA * np.sin(OMEGA * 0.1)], atol=1e-6)


def test_material_derivative_of_a_conserved_quantity_is_zero():
    # x y is constant along the particles of stagnation flow, but not at a fixed point of a moving field
    vel = impl.stagnation_flow(2.0)
    assert abs(impl.material_derivative(lambda x, y, t: x * y, vel, 0.7, 1.3, 0.0)) < 1e-7
    assert impl.material_derivative(lambda x, y, t: x, vel, 0.7, 1.3, 0.0) == pytest.approx(1.4, rel=1e-8)


def test_incompressible_flow_keeps_the_area_of_a_material_loop():
    vel = impl.taylor_green(0.02)
    th = np.linspace(0.0, 2 * np.pi, 400, endpoint=False)
    x, y = 1.2 + 0.3 * np.cos(th), 0.9 + 0.3 * np.sin(th)
    A0 = impl.polygon_area(x, y)
    for i in range(100):                                                    # advect to t = 2
        x, y = impl.rk4_step(vel, x, y, 0.02 * i, 0.02)
    assert impl.polygon_area(x, y) == pytest.approx(A0, rel=1e-3)          # loop resolution limits it


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "flow_lines.png").stat().st_size > 10_000
