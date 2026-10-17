"""Public tests for P06.18 Multibody Simulation. They assert mechanics, not magic numbers: the gas energy
is the integral of the gas force; the orifice force is odd, quadratic and always dissipative; the tyre
cannot pull; the reduced mass matrix is the kinetic energy of the axle kinematics; the generalized
forces are minus the gradient of the potential when nothing dissipates (virtual work); the energy put
in equals the energy stored plus the energy dissipated at every instant, and nothing is dissipated when
the dampers are removed; constrained Newton-Euler reproduces Lagrange, and its multipliers equal the
joint and guide forces computed independently from the Lagrange accelerations; an initial constraint
violation drifts linearly without stabilisation and is pulled back by Baumgarte's method; and the drop
metrics behave sensibly with sink rate.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

G0 = 9.80665
SINK = 3.05
T_END = 0.25                      # covers touchdown to beyond peak stroke (about 0.17 s)


@pytest.fixture(scope="module")
def gear():
    return impl.Gear()


@pytest.fixture(scope="module")
def reduced(gear):
    return impl.simulate_reduced(gear, SINK, t_end=T_END, n_out=501)


@pytest.fixture(scope="module")
def full(gear):
    return impl.simulate_full(gear, SINK, t_end=T_END, n_out=501)


# ---------------------------------------------------------------- force elements

def test_gas_force_and_energy(gear):
    assert impl.gas_force(0.0, gear) == pytest.approx(gear.A_a * (gear.p0 - 101_325.0), rel=1e-12)
    s = np.linspace(-0.01, 0.29, 31)
    assert np.all(np.diff(impl.gas_force(s, gear)) > 0)                     # stiffens as it compresses
    assert impl.gas_energy(0.0, gear) == pytest.approx(0.0, abs=1e-12)
    for s1 in (0.05, 0.2, 0.29):
        num = quad(lambda x: float(impl.gas_force(x, gear)), 0.0, s1, epsabs=1e-10, epsrel=1e-12)[0]
        assert impl.gas_energy(s1, gear) == pytest.approx(num, rel=1e-10)
    # isentropic gas (n = 1.4) is stiffer than near-isothermal (n = 1.1) at the same stroke
    g14 = impl.Gear(n_poly=1.4)
    assert impl.gas_force(0.25, g14) > impl.gas_force(0.25, gear)


def test_orifice_force_is_odd_quadratic_and_dissipative(gear):
    v = np.array([-2.0, -0.5, 0.5, 1.0, 2.0])
    F = impl.orifice_force(v, gear)
    npt.assert_allclose(F[:2], -impl.orifice_force(-v[:2], gear), rtol=1e-14)
    assert impl.orifice_force(2.0, gear) == pytest.approx(4.0 * impl.orifice_force(1.0, gear), rel=1e-14)
    assert np.all(F * v > 0)                                                # power always into the oil
    k = gear.rho * gear.A_h ** 3 / (2 * (gear.Cd * gear.A_o) ** 2)
    assert impl.orifice_force(1.0, gear) == pytest.approx(k, rel=1e-14)


def test_stop_and_tyre_are_one_sided(gear):
    assert impl.stop_force(0.01, gear) == 0.0 and impl.stop_force(-0.001, gear) == pytest.approx(1.0e4)
    assert impl.tyre_force(gear.R + 0.01, -1.0, gear) == (0.0, 0.0)          # in the air
    Ft, Fe = impl.tyre_force(gear.R - 0.02, -1.0, gear)                      # compressing
    assert Fe == pytest.approx(gear.k_t * 0.02) and Ft == pytest.approx(Fe + gear.c_t * 1.0)
    Ft, Fe = impl.tyre_force(gear.R - 1e-4, 10.0, gear)                      # extending fast: cannot pull
    assert Ft == 0.0 and Fe > 0.0


# ---------------------------------------------------------------- Lagrange model

def test_mass_matrix_is_the_kinetic_energy_of_the_kinematics(gear):
    Mr = impl.mass_matrix_reduced(gear)
    npt.assert_allclose(Mr, Mr.T, atol=0)
    assert np.all(np.linalg.eigvalsh(Mr) > 0)
    # T = M y1'^2/2 + m |r2'|^2/2 with r2 from axle_position, differentiated numerically
    y1, s, v = 1.2, 0.1, np.array([-2.0, 1.5])
    h = 1e-6
    x_a, y_a = impl.axle_position(y1 + h * v[0], s + h * v[1], gear)
    x_b, y_b = impl.axle_position(y1 - h * v[0], s - h * v[1], gear)
    r2dot = np.array([x_a - x_b, y_a - y_b]) / (2 * h)
    T = 0.5 * gear.M * v[0] ** 2 + 0.5 * gear.m * r2dot @ r2dot
    assert 0.5 * v @ Mr @ v == pytest.approx(T, rel=1e-8)


def test_generalized_forces_are_minus_the_potential_gradient(gear):
    # with no damping and at rest, Q = -dV/dq: virtual work checked against the potential
    g = impl.Gear(c_t=0.0)
    lift = g.lift_ratio * g.W

    def V(y1, s):
        _, y2 = impl.axle_position(y1, s, g)
        d = max(g.R - y2, 0.0)
        return (g.M * G0 * y1 + g.m * G0 * y2 - lift * y1 + impl.gas_energy(s, g)
                + 0.5 * g.k_t * d * d + 0.5 * g.k_stop * max(-s, 0.0) ** 2)

    for y1, s in ((1.27, 0.05), (1.20, 0.15), (1.30, -0.002)):
        Q = impl.generalized_forces(y1, s, 0.0, 0.0, g)
        h = 1e-7
        grad = [(V(y1 + h, s) - V(y1 - h, s)) / (2 * h), (V(y1, s + h) - V(y1, s - h)) / (2 * h)]
        npt.assert_allclose(Q, -np.array(grad), rtol=1e-6, atol=1e-3)


def test_hanging_stroke(gear):
    s0 = impl.initial_stroke(gear)
    assert -0.01 < s0 < 0.0                                               # just onto the stop
    assert impl.gas_force(s0, gear) - impl.stop_force(s0, gear) == pytest.approx(-gear.m * G0 * np.cos(gear.beta))
    s_ff = impl.initial_stroke(impl.Gear(lift_ratio=0.0))
    assert impl.gas_force(s_ff, gear) == pytest.approx(impl.stop_force(s_ff, gear))   # free fall: no strut load


def test_energy_input_equals_stored_plus_dissipated(gear, reduced):
    en = impl.energy_terms(reduced, gear, SINK)
    scale = en["input"].max()
    assert scale > 0.9 * 0.5 * (gear.M + gear.m) * SINK ** 2              # most of the KE is taken in
    # DOP853 at rtol 1e-10: the budget closes to 1e-7 of the energy absorbed
    assert np.max(np.abs(en["residual"])) < 1e-7 * scale
    assert np.all(np.diff(en["dissipated"]) >= -1e-9)                      # dissipation never runs backwards
    assert np.all(en["stored"] >= -1e-6)


def test_no_dampers_no_dissipation():
    g = impl.Gear(A_o=1.0, c_t=0.0)          # a 1 m^2 orifice: hydraulic force ~ 1e-8 of the default
    res = impl.simulate_reduced(g, SINK, t_end=0.15, n_out=301)
    en = impl.energy_terms(res, g, SINK)
    assert res["D"][-1] < 1e-6 * en["input"].max()
    assert np.max(np.abs(en["input"] - en["stored"])) < 1e-6 * en["input"].max()


# ---------------------------------------------------------------- Newton-Euler against Lagrange

def test_newton_euler_reproduces_lagrange(reduced, full):
    npt.assert_allclose(full["q"][1], reduced["y1"], atol=1e-9)
    npt.assert_allclose(full["s"], reduced["s"], atol=1e-9)
    npt.assert_allclose(full["Ft"], reduced["Ft"], atol=1e-4 * reduced["Ft"].max())
    npt.assert_allclose(full["D"], reduced["D"], atol=1e-6 * reduced["D"].max())
    assert np.max(np.abs(full["phi"])) < 1e-12                             # consistent start stays consistent


def test_multipliers_are_the_joint_and_guide_forces(gear, reduced, full):
    """lambda_4 (bearing side load) and lambda_1 (guide) from Newton-Euler, against the forces needed to
    give the axle the accelerations of the independent Lagrange model."""
    e = np.array([np.sin(gear.beta), np.cos(gear.beta)])
    n = np.array([np.cos(gear.beta), -np.sin(gear.beta)])
    for k in (50, 150, 300):
        z = np.array([reduced["y1"][k], reduced["s"][k], reduced["y1dot"][k], reduced["sdot"][k], 0.0])
        a1, a_s = impl.rhs_reduced(0.0, z, gear)[2:4]
        a2 = np.array([a_s * e[0], a1 + a_s * e[1]])                     # axle acceleration
        F2 = np.array([0.0, -gear.m * G0 + reduced["Ft"][k]])             # applied on the axle (strut force is along e)
        assert full["lam"][3][k] == pytest.approx((gear.m * a2 - F2) @ n, rel=1e-6, abs=1e-3)
        # the only external horizontal force on the two bodies is the guide reaction
        assert full["lam"][0][k] == pytest.approx(gear.m * a2[0], rel=1e-6, abs=1e-3)
    # at peak tyre force the side load is close to Ft sin(beta), less the axle's inertia
    k = int(np.argmax(full["Ft"]))
    assert 0.8 < full["side_load"][k] / (full["Ft"][k] * np.sin(gear.beta)) < 1.0


def test_baumgarte_pulls_back_an_initial_constraint_violation(gear):
    dz = np.zeros(13)
    dz[3], dz[9] = 1.0e-3, 0.05
    c = np.cos(gear.beta)
    free = impl.simulate_full(gear, SINK, t_end=0.1, n_out=11, dz0=dz)
    # without stabilisation Phi'' = 0 is enforced exactly, so Phi = Phi0 + Phi0' t
    npt.assert_allclose(free["phi"][3], 1.0e-3 * c + 0.05 * c * free["t"], rtol=1e-6)
    stab = impl.simulate_full(gear, SINK, t_end=0.2, n_out=21, alpha=50.0, beta_b=50.0, dz0=dz)
    assert abs(stab["phi"][3][-1]) < 0.01 * abs(stab["phi"][3][0])


# ---------------------------------------------------------------- the drop

def test_drop_metrics_and_sink_rate(gear, reduced):
    m = impl.drop_metrics(reduced, gear)
    assert not m["bottomed"] and 0.0 < m["s_max"] < gear.stroke_max
    assert 0.7 < m["efficiency"] < 1.0                                      # an oleo is efficient; a spring would be 0.5
    assert m["n_gear"] == pytest.approx(m["Ft_max"] / gear.W)
    slow = impl.drop_metrics(impl.simulate_reduced(gear, 2.0, t_end=T_END, n_out=501), gear)
    assert slow["n_gear"] < m["n_gear"] and slow["s_max"] < m["s_max"]
    # energy bound: the strut force times the stroke at its peak efficiency must absorb most of the KE
    assert m["Fs_max"] * (m["s_max"] - reduced["s"][0]) * m["efficiency"] > 0.5 * 0.5 * (gear.M + gear.m) * SINK ** 2


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "landing_gear.png").stat().st_size > 10_000


def test_without_lift_the_gear_absorbs_more_energy(gear, reduced):
    # lift = weight is the drop-test convention; with no lift, gravity keeps working through the stroke
    no_lift = impl.Gear(lift_ratio=0.0)
    r0 = impl.simulate_reduced(no_lift, SINK, t_end=T_END, n_out=501)
    assert impl.drop_metrics(r0, no_lift)["n_gear"] > impl.drop_metrics(reduced, gear)["n_gear"]
    assert np.max(np.abs(impl.energy_terms(r0, no_lift, SINK)["residual"])) < 1e-7 * impl.energy_terms(r0, no_lift, SINK)["input"].max()
