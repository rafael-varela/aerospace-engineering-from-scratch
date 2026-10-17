"""Public tests for P04.09 Magnetism & Induction. They assert physics, not magic numbers: the
Lorentz force is perpendicular to v and B; the discretised Biot-Savart sum converges at second order
to the closed form on the axis of a loop and agrees with a ``scipy.integrate.quad`` oracle off the
axis; a long straight polyline reproduces Ampere's mu0 I / (2 pi d); the finite solenoid tends to
mu0 n I and equals the integral of loop fields; the Earth's dipole is divergence-free with the right
1/r^3 and 2:1 pole-to-equator ratio; the torque has no component along B; Faraday's emf is -dPhi/dt;
the RL closed form matches ``solve_ivp``; the B-dot detumble decays exponentially with I / (k B^2).

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad, solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

MU0 = 4e-7 * np.pi


# ---------------------------------------------------------------- Lorentz force

def test_lorentz_force_is_perpendicular_to_v_and_B_and_vanishes_for_parallel_motion():
    rng = np.random.default_rng(409)
    for _ in range(5):
        q = rng.uniform(-2, 2)
        v, B = rng.normal(size=3), rng.normal(size=3)
        F = impl.lorentz_force(q, v, np.zeros(3), B)
        assert abs(np.dot(F, v)) < 1e-12 * np.linalg.norm(F) * np.linalg.norm(v) + 1e-15
        assert abs(np.dot(F, B)) < 1e-12 * np.linalg.norm(F) * np.linalg.norm(B) + 1e-15
        sin_th = np.linalg.norm(np.cross(v, B)) / (np.linalg.norm(v) * np.linalg.norm(B))
        assert np.linalg.norm(F) == pytest.approx(abs(q) * np.linalg.norm(v) * np.linalg.norm(B) * sin_th)
    # parallel motion: no magnetic force; the electric term adds linearly
    npt.assert_allclose(impl.lorentz_force(1.0, [0, 0, 3.0], np.zeros(3), [0, 0, 2.0]), 0.0, atol=1e-15)
    npt.assert_allclose(impl.lorentz_force(2.0, [1.0, 0, 0], [0, 0, 5.0], [0, 1.0, 0]), [0, 0, 2 * (5.0 + 1.0)])


def test_gyroradius_and_cyclotron_frequency_are_consistent():
    m, v, q, B = impl.M_PROTON, 1e6, impl.Q_PROTON, 3e-5
    r = impl.gyroradius(m, v, q, B)
    w = impl.cyclotron_frequency(q, B, m)
    assert r * w == pytest.approx(v)                             # v = omega r
    assert impl.gyroradius(m, 2 * v, q, B) == pytest.approx(2 * r)
    assert impl.gyroradius(m, v, q, 2 * B) == pytest.approx(r / 2)
    assert impl.cyclotron_frequency(-q, B, m) == pytest.approx(w)  # sign of q does not change the rate
    assert 100.0 < r < 1000.0                                      # a few hundred metres in LEO


# ---------------------------------------------------------------- Biot-Savart

def test_loop_axis_closed_form_limits():
    R, I = 0.05, 2.0
    assert impl.loop_field_axis(0.0, R, I) == pytest.approx(MU0 * I / (2 * R))
    z = 40 * R
    m = I * np.pi * R ** 2
    assert impl.loop_field_axis(z, R, I) == pytest.approx(MU0 * m / (2 * np.pi * z ** 3), rel=2e-3)   # dipole far field
    zz = np.array([-0.1, 0.0, 0.1])
    b = impl.loop_field_axis(zz, R, I)
    assert b.shape == (3,) and b[0] == pytest.approx(b[2])        # even in z


def test_segmented_loop_converges_to_the_axis_closed_form_at_second_order():
    R, I = 0.04, 0.12
    z = 0.5 * R
    exact = impl.loop_field_axis(z, R, I)
    errs = []
    for n in (16, 32, 64, 128):
        B = impl.loop_field(np.array([0.0, 0.0, z]), R, I, n)
        assert abs(B[0]) < 1e-12 * exact and abs(B[1]) < 1e-12 * exact   # purely axial on the axis
        errs.append(abs(B[2] - exact) / exact)
    orders = np.diff(np.log(errs)) / np.log(0.5)
    assert np.all(orders > 1.8) and np.all(orders < 2.2)
    assert errs[-1] < 1e-3


def test_segmented_loop_matches_scipy_quad_off_axis():
    R, I = 0.04, 0.12
    rho, z = 0.5 * R, 0.3 * R
    # oracle built here, independently of the implementation: exact Biot-Savart integrand over phi'
    d3 = lambda p: (rho ** 2 + R ** 2 - 2 * rho * R * np.cos(p) + z ** 2) ** 1.5
    bz = MU0 * I / (4 * np.pi) * quad(lambda p: R * (R - rho * np.cos(p)) / d3(p), 0, 2 * np.pi)[0]
    br = MU0 * I / (4 * np.pi) * quad(lambda p: R * z * np.cos(p) / d3(p), 0, 2 * np.pi)[0]
    B = impl.loop_field(np.array([rho, 0.0, z]), R, I, 720)
    assert B[0] == pytest.approx(br, rel=1e-4)
    assert B[2] == pytest.approx(bz, rel=1e-4)
    assert abs(B[1]) < 1e-12 * abs(bz)                              # B_phi = 0 by symmetry
    q = impl.loop_field_quad(rho, z, R, I)
    assert q[0] == pytest.approx(br, rel=1e-9) and q[1] == pytest.approx(bz, rel=1e-9)


def test_loop_field_symmetry_and_batch_shape():
    R, I = 0.03, 1.0
    pts = np.array([[0.01, 0.0, 0.02], [0.01, 0.0, -0.02], [0.0, 0.01, 0.02]])
    B = impl.loop_field(pts, R, I, 180)
    assert B.shape == (3, 3)
    assert B[0, 2] == pytest.approx(B[1, 2])                       # B_z even in z
    assert B[0, 0] == pytest.approx(-B[1, 0])                      # B_rho odd in z
    assert B[2, 1] == pytest.approx(B[0, 0]) and B[2, 2] == pytest.approx(B[0, 2])   # rotation about z
    # reversing the current reverses the field; doubling it doubles the field
    npt.assert_allclose(impl.loop_field(pts[0], R, -I, 180), -B[0], rtol=1e-12)
    npt.assert_allclose(impl.loop_field(pts[0], R, 2 * I, 180), 2 * B[0], rtol=1e-12)


def test_long_straight_polyline_reproduces_amperes_law():
    d, I = 0.05, 1.5
    pts = np.column_stack([np.zeros(8001), np.zeros(8001), np.linspace(-20.0, 20.0, 8001)])   # 5 mm segments at 5 cm
    B = impl.wire_field(pts, I, np.array([d, 0.0, 0.0]))
    assert np.linalg.norm(B) == pytest.approx(impl.straight_wire_field(d, I), rel=1e-3)
    assert B[1] > 0 and abs(B[0]) < 1e-12 * B[1] and abs(B[2]) < 1e-12 * B[1]   # azimuthal, right-hand rule
    assert impl.straight_wire_field(d, I) == pytest.approx(MU0 * I / (2 * np.pi * d))


def test_solenoid_limits_and_consistency_with_loop_integral():
    R, I, n_pm, L = 0.005, 1.0, 2000.0, 1.0
    assert impl.solenoid_field_axis(0.0, R, I, n_pm, L) == pytest.approx(MU0 * n_pm * I, rel=1e-4)       # long: mu0 n I
    assert impl.solenoid_field_axis(0.5 * L, R, I, n_pm, L) == pytest.approx(0.5 * MU0 * n_pm * I, rel=1e-4)  # end: half
    # short solenoid: integral of loop fields over the windings (independent quad oracle)
    L = 0.02
    z0 = 0.004
    ref = quad(lambda zp: n_pm * impl.loop_field_axis(z0 - zp, R, I), -L / 2, L / 2)[0]
    assert impl.solenoid_field_axis(z0, R, I, n_pm, L) == pytest.approx(ref, rel=1e-8)
    # inductance scales as N^2 and as A / L
    assert impl.solenoid_inductance(200, 1e-4, 0.1) == pytest.approx(4 * impl.solenoid_inductance(100, 1e-4, 0.1))
    assert impl.solenoid_inductance(100, 2e-4, 0.2) == pytest.approx(impl.solenoid_inductance(100, 1e-4, 0.1))


# ---------------------------------------------------------------- the Earth's dipole

def test_dipole_field_has_2_to_1_pole_to_equator_ratio_and_falls_as_r_cubed():
    Re = impl.R_EARTH
    Beq = impl.earth_dipole_field(np.array([Re, 0.0, 0.0]))
    Bpole = impl.earth_dipole_field(np.array([0.0, 0.0, Re]))
    assert np.linalg.norm(Beq) == pytest.approx(impl.B0_EQUATOR, rel=1e-12)
    assert np.linalg.norm(Bpole) == pytest.approx(2 * impl.B0_EQUATOR, rel=1e-12)
    assert Beq[2] > 0 and Bpole[2] < 0                              # north at the equator, downward at the north pole
    B2 = impl.earth_dipole_field(np.array([2 * Re, 0.0, 0.0]))
    assert np.linalg.norm(B2) == pytest.approx(impl.B0_EQUATOR / 8)
    # magnitude at 500 km is between 20 and 50 uT everywhere
    _, Bo = impl.field_along_orbit(500e3, 98.0)
    mag = np.linalg.norm(Bo, axis=1)
    assert mag.min() > 20e-6 and mag.max() < 50e-6


def test_dipole_field_is_divergence_free_and_curl_free():
    m = np.array([0.3, -0.5, 0.8]) * 1e20
    r0 = np.array([5e6, -3e6, 4e6])
    h = 1.0
    def B(r):
        return impl.dipole_field(r, m)
    div = 0.0
    curl = np.zeros(3)
    grad = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = h
        grad[:, j] = (B(r0 + e) - B(r0 - e)) / (2 * h)
    div = np.trace(grad)
    curl = np.array([grad[2, 1] - grad[1, 2], grad[0, 2] - grad[2, 0], grad[1, 0] - grad[0, 1]])
    scale = np.linalg.norm(B(r0)) / np.linalg.norm(r0)
    assert abs(div) < 1e-6 * scale
    assert np.linalg.norm(curl) < 1e-6 * scale
    # batch evaluation
    assert impl.dipole_field(np.vstack([r0, 2 * r0]), m).shape == (2, 3)


def test_field_along_orbit_equatorial_is_constant_and_polar_has_2_to_1_ratio():
    u, B0 = impl.field_along_orbit(500e3, 0.0)
    mag0 = np.linalg.norm(B0, axis=1)
    assert np.ptp(mag0) < 1e-12 * mag0.mean()
    assert mag0[0] == pytest.approx(impl.B0_EQUATOR * (impl.R_EARTH / (impl.R_EARTH + 500e3)) ** 3)
    u, B90 = impl.field_along_orbit(500e3, 90.0)
    mag90 = np.linalg.norm(B90, axis=1)
    assert mag90.max() / mag90.min() == pytest.approx(2.0, rel=1e-6)
    pos = impl.orbit_positions(500e3, 63.4, 100)
    npt.assert_allclose(np.linalg.norm(pos, axis=1), impl.R_EARTH + 500e3, rtol=1e-12)


# ---------------------------------------------------------------- torque

def test_torque_has_no_component_along_B_and_peaks_at_m_B():
    rng = np.random.default_rng(1)
    for _ in range(5):
        m, B = rng.normal(size=3), rng.normal(size=3) * 1e-5
        tau = impl.torque_on_dipole(m, B)
        assert abs(np.dot(tau, B)) < 1e-12 * np.linalg.norm(tau) * np.linalg.norm(B)
        assert abs(np.dot(tau, m)) < 1e-12 * np.linalg.norm(tau) * np.linalg.norm(m)
        assert np.linalg.norm(tau) <= np.linalg.norm(m) * np.linalg.norm(B) * (1 + 1e-12)
    npt.assert_allclose(impl.torque_on_dipole([0, 0, 1.0], [0, 0, 3e-5]), 0.0, atol=1e-20)
    npt.assert_allclose(impl.torque_on_dipole([1.0, 0, 0], [0, 3e-5, 0]), [0, 0, 3e-5])
    assert impl.magnetic_moment(200, 0.1, 0.0064) == pytest.approx(0.128)


def test_torquer_coil_resistance_power_and_moment_scale_correctly():
    R1 = impl.coil_resistance(100, 0.08, 0.25e-3)
    assert impl.coil_resistance(200, 0.08, 0.25e-3) == pytest.approx(2 * R1)
    assert impl.coil_resistance(100, 0.08, 0.5e-3) == pytest.approx(R1 / 4)
    assert R1 == pytest.approx(impl.RHO_COPPER * 100 * 4 * 0.08 / (np.pi * (0.25e-3) ** 2 / 4))
    tq = impl.Torquer(N=200, side=0.08, wire_diameter=0.25e-3, I=0.12)
    assert tq.moment == pytest.approx(200 * 0.12 * 0.08 ** 2)
    assert tq.power == pytest.approx(0.12 ** 2 * tq.resistance)
    assert 0.1 < tq.power < 1.0 and 0.05 < tq.moment < 0.5          # a CubeSat-class coil


# ---------------------------------------------------------------- Faraday, Lenz, RL

def test_faraday_emf_is_minus_dphi_dt_and_generator_peak_is_NABomega():
    t = np.linspace(0.0, 0.1, 20001)
    N, A, B, w = 50, 2e-3, 0.2, 2 * np.pi * 60
    flux = N * A * B * np.cos(w * t)
    emf = impl.faraday_emf(t, flux)
    npt.assert_allclose(emf, impl.generator_emf(t, N, A, B, w), atol=2e-3 * N * A * B * w)
    assert emf.max() == pytest.approx(N * A * B * w, rel=1e-3)
    # Lenz: flux increasing linearly gives a constant negative emf
    npt.assert_allclose(impl.faraday_emf(t, 3.0 * t)[1:-1], -3.0, atol=1e-9)
    assert impl.motor_back_emf(0.02, 300.0) == pytest.approx(6.0)


def test_rl_closed_form_matches_solve_ivp_and_time_constant():
    V, R, L = 5.0, 20.0, 4e-3
    tau = L / R
    t = np.linspace(0.0, 5 * tau, 201)
    i_cf = impl.rl_current(t, V, R, L)
    assert i_cf[0] == pytest.approx(0.0, abs=1e-15)
    assert impl.rl_current(tau, V, R, L) == pytest.approx((V / R) * (1 - np.exp(-1)))
    assert np.gradient(i_cf, t)[0] == pytest.approx(V / L, rel=2e-2)       # initial slope V / L
    sol = solve_ivp(lambda tt, y: [(V - R * y[0]) / L], (0, 5 * tau), [0.0], t_eval=t, rtol=1e-10, atol=1e-13)
    npt.assert_allclose(i_cf, sol.y[0], atol=1e-9)
    ts, i_s = impl.rl_simulate(5 * tau, lambda tt: V, R, L, n=201)
    npt.assert_allclose(i_s, i_cf, atol=1e-8)
    # with initial current the decay to V/R has the same tau
    assert impl.rl_current(tau, V, R, L, i0=1.0) == pytest.approx(V / R + (1.0 - V / R) * np.exp(-1))


# ---------------------------------------------------------------- detumble

def test_bdot_detumble_decays_exponentially_when_unsaturated():
    I_sc, w0, B, k = 0.04, 0.05, 3e-5, 5e4
    T_d = impl.detumble_time_constant(I_sc, k, B)
    assert T_d == pytest.approx(I_sc / (k * B ** 2))
    m_max = 10.0                                                   # never saturates
    t, w = impl.bdot_detumble(I_sc, w0, B, k, m_max, 3 * T_d, n=301)
    npt.assert_allclose(w, w0 * np.exp(-t / T_d), rtol=1e-6, atol=1e-12)
    assert np.all(np.diff(w) < 0)


def test_bdot_detumble_saturated_phase_is_linear_and_the_whole_run_is_monotone():
    I_sc, w0, B, k, m_max = 0.04, np.radians(10.0), 3e-5, 5e4, 0.15
    w_s = impl.saturation_rate(m_max, k, B)
    assert w_s == pytest.approx(m_max / (k * B))
    assert w_s < w0                                                # starts saturated
    t, w = impl.bdot_detumble(I_sc, w0, B, k, m_max, 4000.0, n=4001)
    slope = m_max * B / I_sc
    t_sat = (w0 - w_s) / slope
    early = t < 0.8 * t_sat
    npt.assert_allclose(w[early], w0 - slope * t[early], rtol=1e-5)
    assert np.all(np.diff(w) <= 1e-12)
    assert w[-1] > 0.0 and w[-1] < 0.1 * w0
    assert t_sat > impl.momentum_dump_time(I_sc, w0 - w_s, m_max * B) * (1 - 1e-9)
    assert impl.momentum_dump_time(I_sc, w0, m_max * B) == pytest.approx(I_sc * w0 / (m_max * B))


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "magnetorquer.png").exists()
