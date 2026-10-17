"""Public tests for P04.03 Momentum & Collisions. They assert physics, not magic numbers: momentum is
conserved in every collision and kinetic energy only in the elastic one, the energy loss formula, the
impulse-momentum theorem, the pellet rocket's convergence to Tsiolkovsky at order 1/n with the derived
constant, the continuous rocket against the closed-form burnout state, conservation of the momentum of
rocket plus exhaust in free space, the gravity-loss law, and the staging results.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-04/parameterized-cases/.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

G = impl.G0
U = 300.0 * G
M0, MF = 10_000.0, 2_500.0


# ---------------------------------------------------------------- collisions and impulse

@pytest.mark.parametrize("m1,v1,m2,v2", [(1.0, 2.0, 3.0, 0.0), (5.0, -1.0, 2.0, 4.0), (12e3, 0.1, 420e3, 0.0)])
def test_collisions_conserve_momentum_and_elastic_ones_conserve_energy(m1, v1, m2, v2):
    for e in (0.0, 0.3, 0.7, 1.0):
        a, b = impl.collide_1d(m1, v1, m2, v2, e)
        assert m1 * a + m2 * b == pytest.approx(m1 * v1 + m2 * v2, rel=1e-12, abs=1e-12)
        assert b - a == pytest.approx(-e * (v2 - v1), abs=1e-12)
        K0 = 0.5 * m1 * v1 ** 2 + 0.5 * m2 * v2 ** 2
        K1 = 0.5 * m1 * a ** 2 + 0.5 * m2 * b ** 2
        assert K0 - K1 == pytest.approx(impl.collision_energy_loss(m1, v1, m2, v2, e), rel=1e-9, abs=1e-9)
        assert K1 <= K0 + 1e-9
    a, b = impl.collide_1d(m1, v1, m2, v2, 0.0)
    assert a == pytest.approx(b)                                          # perfectly inelastic: they move together


def test_collision_limits_and_validation():
    a, b = impl.collide_1d(2.0, 3.0, 2.0, 0.0, 1.0)                     # equal masses, elastic: velocities swap
    assert (a, b) == (pytest.approx(0.0, abs=1e-15), pytest.approx(3.0))
    a, b = impl.collide_1d(1.0, 1.0, 1e12, 0.0, 1.0)                    # ball on a wall: bounces back
    assert a == pytest.approx(-1.0, rel=1e-9) and b == pytest.approx(0.0, abs=1e-9)
    # in the centre-of-mass frame the e = 0 loss is all the kinetic energy
    m1, m2, v1 = 3.0, 1.0, 2.0
    vcm = m1 * v1 / (m1 + m2)
    loss = impl.collision_energy_loss(m1, v1 - vcm, m2, -vcm, 0.0)
    assert loss == pytest.approx(0.5 * m1 * (v1 - vcm) ** 2 + 0.5 * m2 * vcm ** 2)
    with pytest.raises(ValueError):
        impl.collide_1d(1.0, 1.0, 1.0, 0.0, 1.5)


def test_impulse_momentum_theorem_and_centre_of_mass():
    t = np.linspace(0.0, 0.01, 1001)
    F = 2000.0 * np.sin(np.pi * t / 0.01)                                  # a half-sine kick
    J = impl.impulse(t, F)
    assert J == pytest.approx(2000.0 * 0.01 * 2 / np.pi, rel=1e-5)
    m = 0.45                                                               # the kick on a 0.45 kg ball at rest
    assert J / m == pytest.approx(28.29, rel=1e-3)
    assert impl.centre_of_mass([1.0, 3.0], [0.0, 4.0]) == pytest.approx(3.0)
    npt.assert_allclose(impl.centre_of_mass([1.0, 1.0, 2.0], [[0, 0], [2, 0], [1, 2]]), [1.0, 1.0])


# ---------------------------------------------------------------- the pellet rocket

def test_one_pellet_is_a_single_momentum_balance():
    v = impl.pellet_rocket_velocities(M0, MF, U, 1)
    assert v.shape == (2,) and v[0] == 0.0
    # (m0 - mp) v + mp (v - u) = 0
    assert v[1] == pytest.approx(U * (M0 - MF) / M0)
    v = impl.pellet_rocket_velocities(M0, MF, U, 50)
    assert v.shape == (51,)
    assert np.all(np.diff(v) > 0) and np.all(np.diff(np.diff(v)) > 0)      # each throw gains more as the rocket lightens


def test_pellets_converge_to_tsiolkovsky_at_first_order_with_the_derived_constant():
    exact = impl.tsiolkovsky_dv(U, M0, MF)
    assert exact == pytest.approx(U * math.log(4.0))
    ns = [10, 100, 1000, 10000]
    errs = impl.pellet_convergence(M0, MF, U, ns)
    assert np.all(errs > 0)                                               # always slower than the continuous rocket
    slope = np.polyfit(np.log(ns), np.log(errs), 1)[0]
    assert slope == pytest.approx(-1.0, abs=0.02)
    C = U * (M0 - MF) * (1 / MF - 1 / M0) / 2
    assert errs[-1] * ns[-1] == pytest.approx(C, rel=1e-3)


def test_mass_ratio_inverts_tsiolkovsky():
    for dv in (1000.0, 4000.0, 9400.0):
        R = impl.mass_ratio_for_dv(dv, U)
        assert impl.tsiolkovsky_dv(U, R, 1.0) == pytest.approx(dv, rel=1e-12)
    assert impl.tsiolkovsky_dv(U, 5.0, 5.0) == 0.0
    # doubling the exhaust speed takes the square root of the mass ratio
    assert impl.mass_ratio_for_dv(9400.0, 2 * U) == pytest.approx(math.sqrt(impl.mass_ratio_for_dv(9400.0, U)))


# ---------------------------------------------------------------- the continuous rocket

def test_free_space_rocket_reaches_tsiolkovsky_and_conserves_total_momentum():
    fl = impl.simulate_rocket(M0, MF, 50.0, U, g=0.0)
    assert fl.v[-1] == pytest.approx(impl.tsiolkovsky_dv(U, M0, MF), rel=1e-9)
    P = fl.m * fl.v + fl.p_exhaust
    assert np.max(np.abs(P)) < 1e-6 * M0 * U                               # zero initial momentum stays zero
    assert fl.p_exhaust[-1] < 0                                           # the exhaust carries momentum backward
    npt.assert_allclose(fl.m, M0 - 50.0 * fl.t, rtol=1e-12)
    # v(m) = u ln(m0/m) all along the burn, not just at the end
    npt.assert_allclose(fl.v, U * np.log(M0 / fl.m), rtol=1e-8, atol=1e-6)


def test_vertical_ascent_matches_the_closed_form_burnout_state():
    for mdot in (40.0, 50.0, 120.0):                                       # T/W0 = 1.2, 1.5, 3.6
        fl = impl.simulate_rocket(M0, MF, mdot, U)
        vb, zb = impl.burnout_closed_form(M0, MF, mdot, U)
        assert fl.v[-1] == pytest.approx(vb, rel=1e-9)
        assert fl.z[-1] == pytest.approx(zb, rel=1e-9)
        assert fl.burn_time == pytest.approx((M0 - MF) / mdot)
        ideal, grav, drag = impl.velocity_losses(fl, U)
        assert grav == pytest.approx(G * fl.burn_time)
        assert abs(drag) < 1e-5
    # g -> 0 limit of the closed form
    vb, zb = impl.burnout_closed_form(M0, MF, 50.0, U, g=0.0)
    assert vb == pytest.approx(impl.tsiolkovsky_dv(U, M0, MF))


def test_gravity_loss_falls_with_thrust_to_weight_and_drag_costs_velocity():
    vb = []
    for tw in (1.2, 2.0, 5.0):
        mdot = tw * M0 * G / U
        fl = impl.simulate_rocket(M0, MF, mdot, U, n_samples=5)
        vb.append(fl.v[-1])
        assert impl.velocity_losses(fl, U)[1] == pytest.approx(G * (M0 - MF) / mdot)
    assert vb[0] < vb[1] < vb[2] < impl.tsiolkovsky_dv(U, M0, MF)
    k = 0.5 * 1.225 * 0.3 * math.pi * 0.25
    fd = impl.simulate_rocket(M0, MF, 50.0, U, k_drag=k)
    fn = impl.simulate_rocket(M0, MF, 50.0, U)
    loss = impl.velocity_losses(fd, U)[2]
    assert 0 < loss < 0.1 * impl.tsiolkovsky_dv(U, M0, MF)
    assert fd.v[-1] == pytest.approx(fn.v[-1] - loss, rel=1e-9)
    # constant density (H = inf) costs far more than an atmosphere that thins with height
    fc = impl.simulate_rocket(M0, MF, 50.0, U, k_drag=k, scale_height=math.inf)
    assert impl.velocity_losses(fc, U)[2] > 5 * loss


def test_rocket_that_cannot_lift_off_is_rejected():
    with pytest.raises(ValueError):
        impl.simulate_rocket(M0, MF, 0.9 * M0 * G / U, U)
    fl = impl.simulate_rocket(M0, MF, 0.9 * M0 * G / U, U, g=0.0)           # fine in free space
    assert fl.v[-1] > 0


def test_rhs_thrust_is_momentum_flux():
    y = np.array([0.0, 100.0, 5000.0, 0.0])
    d = impl.rocket_rhs(0.0, y, 20.0, 3000.0, 0.0, 0.0)
    assert d[1] == pytest.approx(20.0 * 3000.0 / 5000.0)                  # a = mdot u / m
    assert d[2] == -20.0
    assert d[3] == pytest.approx(20.0 * (100.0 - 3000.0))
    # momentum balance of the rhs: d(m v)/dt + dP_ex/dt = external force (zero here)
    assert y[2] * d[1] + y[1] * d[2] + d[3] == pytest.approx(0.0, abs=1e-9)
    d = impl.rocket_rhs(0.0, y, 20.0, 3000.0, G, 0.0)
    assert y[2] * d[1] + y[1] * d[2] + d[3] == pytest.approx(-y[2] * G)      # ... and -m g with gravity


# ---------------------------------------------------------------- staging

def test_stack_dv_single_stage_and_two_stage_identities():
    total, parts = impl.stack_dv(500.0, [(20_000.0, 2_000.0, U)])
    assert total == pytest.approx(U * math.log(22_500.0 / 2_500.0))
    assert len(parts) == 1
    st = [(15_000.0, 1_500.0, U), (5_000.0, 500.0, 1.2 * U)]
    total, parts = impl.stack_dv(300.0, st)
    m01 = 300.0 + 22_000.0
    assert parts[0] == pytest.approx(U * math.log(m01 / (m01 - 15_000.0)))
    assert parts[1] == pytest.approx(1.2 * U * math.log(5_800.0 / 800.0))
    assert total == pytest.approx(sum(parts))
    assert impl.structure_mass(9_000.0, 0.1) == pytest.approx(1_000.0)


def test_staging_beats_a_single_stage_and_the_single_stage_limit():
    payload, mp, eps = 500.0, 20_000.0, 0.10
    one = impl.stack_dv(payload, [(mp, impl.structure_mass(mp, eps), U)])[0]
    f, two = impl.best_two_stage_split(payload, mp, eps, U)
    assert two > one
    assert one < U * math.log(1 / eps)                                     # a single stage can never beat u ln(1/eps)
    assert two > U * math.log(1 / eps)                                     # ... but two stages can
    assert 0.5 < f < 0.99                                                  # the first stage carries most of the propellant
    # the optimum is a maximum
    def dv_at(ff):
        return impl.stack_dv(payload, [(ff * mp, impl.structure_mass(ff * mp, eps), U), ((1 - ff) * mp, impl.structure_mass((1 - ff) * mp, eps), U)])[0]
    assert dv_at(f) >= dv_at(f - 0.03) and dv_at(f) >= dv_at(f + 0.03)
    # with equal u and eps the optimum has equal stage mass ratios
    m1, m2 = f * mp, (1 - f) * mp
    s1, s2 = impl.structure_mass(m1, eps), impl.structure_mass(m2, eps)
    R1 = (payload + m1 + m2 + s1 + s2) / (payload + m2 + s1 + s2)
    R2 = (payload + m2 + s2) / (payload + s2)
    assert R1 == pytest.approx(R2, rel=1e-3)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "rocket_momentum.png").exists()
