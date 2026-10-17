"""Public tests for P04.08 Electrostatics & Circuits. They assert Coulomb's law and superposition, E = -grad V,
the dipole far field, resistance scaling, series/parallel/divider identities, Kirchhoff's current law and
the power balance at every nodal solution, the Wheatstone bridge against a hand-written two-loop solution,
linearity, singular-circuit detection, meter loading in its limits, the RC transient against solve_ivp,
and the half-energy rule of capacitor charging.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. The only file written goes to tmp_path.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import solve_ivp

from conftest import load_impl

impl = load_impl(__file__)

K = impl.K_COULOMB


# ---------------------------------------------------------------- electrostatics

def test_coulomb_inverse_square_signs_and_third_law():
    r1, r2 = np.zeros(3), np.array([0.0, 0.0, 2.0])
    F = impl.coulomb_force(1e-6, 1e-6, r1, r2)
    assert F[2] > 0.0                                                      # like charges repel: pushed away from r1
    assert np.linalg.norm(F) == pytest.approx(K * 1e-12 / 4.0)
    assert np.linalg.norm(impl.coulomb_force(1e-6, 1e-6, r1, 2 * r2)) == pytest.approx(np.linalg.norm(F) / 4)
    assert impl.coulomb_force(1e-6, -1e-6, r1, r2)[2] < 0.0                # unlike charges attract
    npt.assert_allclose(impl.coulomb_force(1e-6, -3e-6, r2, r1), -impl.coulomb_force(-3e-6, 1e-6, r1, r2))


def test_field_is_minus_gradient_of_potential_and_superposes():
    rng = np.random.default_rng(8)
    q = rng.uniform(-1e-9, 1e-9, 4)
    pos = rng.uniform(-0.1, 0.1, (4, 3))
    r = np.array([0.25, -0.15, 0.3])
    E = impl.electric_field(q, pos, r)
    npt.assert_allclose(E, impl.field_from_potential(q, pos, r), rtol=1e-6)
    # superposition: the field of all four is the sum of the fields of each
    E_sum = sum(impl.electric_field(q[i:i + 1], pos[i:i + 1], r) for i in range(4))
    npt.assert_allclose(E, E_sum, rtol=1e-12)
    # a single charge: k q / r^2 radially, k q / r for the potential
    assert impl.electric_potential([2e-9], [[0.0, 0.0, 0.0]], [0.0, 0.5, 0.0]) == pytest.approx(K * 2e-9 / 0.5)
    npt.assert_allclose(impl.electric_field([2e-9], [[0.0, 0.0, 0.0]], [0.0, 0.5, 0.0]), [0.0, K * 2e-9 / 0.25, 0.0])


def test_dipole_far_field_falls_as_inverse_cube():
    q = np.array([1e-9, -1e-9])
    pos = np.array([[0.0, 0.0, 0.5e-3], [0.0, 0.0, -0.5e-3]])
    p = impl.dipole_moment(q, pos)
    assert p[2] == pytest.approx(1e-12)
    z1, z2 = 0.05, 0.10
    E1 = impl.electric_field(q, pos, [0.0, 0.0, z1])[2]
    E2 = impl.electric_field(q, pos, [0.0, 0.0, z2])[2]
    assert E1 / E2 == pytest.approx(8.0, rel=1e-3)                        # 1/z^3, not 1/z^2
    assert E1 == pytest.approx(2 * K * p[2] / z1 ** 3, rel=1e-3)
    assert impl.electric_potential(q, pos, [0.1, 0.0, 0.0]) == pytest.approx(0.0, abs=1e-12)   # equatorial plane


# ---------------------------------------------------------------- resistance

def test_resistance_scaling_and_combinations():
    assert impl.resistance(impl.RHO_COPPER, 2.0, 1e-6) == pytest.approx(2 * impl.resistance(impl.RHO_COPPER, 1.0, 1e-6))
    assert impl.wire_resistance(1.0, 2e-3) == pytest.approx(impl.wire_resistance(1.0, 1e-3) / 4)       # ~ 1/d^2
    assert impl.wire_resistance(10.0, 1.291e-3) == pytest.approx(0.128, abs=0.003)                   # AWG 16, ~13 mohm/m
    assert impl.series([1.0, 2.0, 3.0]) == 6.0
    assert impl.parallel([2.0, 2.0]) == pytest.approx(1.0)
    assert impl.parallel([10.0, 1e12]) == pytest.approx(10.0, rel=1e-9)                              # an open branch
    assert impl.voltage_divider(12.0, 10e3, 4.7e3) == pytest.approx(12.0 * 4.7 / 14.7)
    assert impl.power_dissipated(5.0, 0.02) == pytest.approx(0.1)


# ---------------------------------------------------------------- the solver

def _divider():
    return impl.Circuit().add(impl.VoltageSource(1, 0, 12.0)).add(impl.Resistor(1, 2, 10e3)).add(impl.Resistor(2, 0, 4.7e3))


def test_nodal_solver_reproduces_the_divider_and_ohms_law():
    c = _divider()
    sol = impl.solve_dc(c)
    assert sol.node_voltages[0] == 0.0
    assert sol.node_voltages[1] == pytest.approx(12.0)
    assert sol.node_voltages[2] == pytest.approx(impl.voltage_divider(12.0, 10e3, 4.7e3))
    I = impl.branch_current(c.resistors[0], sol.node_voltages)
    assert I == pytest.approx(12.0 / 14.7e3)
    assert sol.source_currents[0] == pytest.approx(I)                      # the source delivers the same current
    A, z = impl.build_mna(c)
    npt.assert_allclose(A, A.T)                                            # MNA of a resistor network is symmetric
    assert A.shape == (3, 3) and z[-1] == 12.0


def test_kcl_and_power_balance_hold_on_a_random_network():
    rng = np.random.default_rng(88)
    c = impl.Circuit().add(impl.VoltageSource(1, 0, 9.0)).add(impl.CurrentSource(0, 3, 2e-3))
    nodes = 5
    for _ in range(14):
        a, b = rng.integers(0, nodes, 2)
        if a != b:
            c.add(impl.Resistor(int(a), int(b), float(rng.uniform(100.0, 5000.0))))
    sol = impl.solve_dc(c)
    npt.assert_allclose(impl.kcl_residuals(c, sol), 0.0, atol=1e-12)
    delivered, dissipated = impl.power_balance(c, sol)
    assert delivered == pytest.approx(dissipated, rel=1e-10)
    assert dissipated > 0.0


def test_wheatstone_bridge_matches_a_hand_written_two_loop_solution_and_balances():
    Vs, R1, R2, R3, Rx, Rg = 5.0, 1000.0, 2000.0, 1000.0, 2200.0, 500.0
    c, ig = impl.wheatstone_bridge(Vs, R1, R2, R3, Rx, Rg)
    # Kirchhoff by hand, loop currents i1 (Vs-R1-R2), i2 (R1-Rg-R3), i3 (R2-Rx-Rg): three KVL equations
    M = np.array([[R1 + R2, -R1, -R2],
                  [-R1, R1 + Rg + R3, -Rg],
                  [-R2, -Rg, R2 + Rx + Rg]])
    i1, i2, i3 = np.linalg.solve(M, [Vs, 0.0, 0.0])
    # with these loop orientations i2 crosses the galvanometer from 3 to 2 and i3 from 2 to 3
    assert ig == pytest.approx(i3 - i2, rel=1e-10)
    V = impl.solve_dc(c).node_voltages
    assert ig == pytest.approx((V[2] - V[3]) / Rg, rel=1e-12)
    assert abs(ig) == pytest.approx(56.18e-6, rel=2e-3)
    # balanced when R1/R2 = R3/Rx: no galvanometer current, whatever Rg
    for Rgal in (10.0, 500.0, 1e6):
        assert impl.wheatstone_bridge(Vs, R1, R2, R3, 2000.0, Rgal)[1] == pytest.approx(0.0, abs=1e-15)
    # the sign flips when Rx crosses the balance value
    assert impl.wheatstone_bridge(Vs, R1, R2, R3, 1800.0, Rg)[1] * ig < 0.0


def test_solver_is_linear_in_the_sources():
    c1 = _divider()
    c2 = impl.Circuit().add(impl.VoltageSource(1, 0, 24.0)).add(impl.Resistor(1, 2, 10e3)).add(impl.Resistor(2, 0, 4.7e3))
    npt.assert_allclose(impl.solve_dc(c2).node_voltages, 2 * impl.solve_dc(c1).node_voltages)
    # superposition of a voltage and a current source
    ca = impl.Circuit().add(impl.VoltageSource(1, 0, 5.0)).add(impl.Resistor(1, 2, 1e3)).add(impl.Resistor(2, 0, 2e3))
    cb = impl.Circuit().add(impl.VoltageSource(1, 0, 0.0)).add(impl.Resistor(1, 2, 1e3)).add(impl.Resistor(2, 0, 2e3)).add(impl.CurrentSource(0, 2, 1e-3))
    cab = impl.Circuit().add(impl.VoltageSource(1, 0, 5.0)).add(impl.Resistor(1, 2, 1e3)).add(impl.Resistor(2, 0, 2e3)).add(impl.CurrentSource(0, 2, 1e-3))
    npt.assert_allclose(impl.solve_dc(cab).node_voltages, impl.solve_dc(ca).node_voltages + impl.solve_dc(cb).node_voltages, atol=1e-12)
    # a current source into a resistor to ground: V = I R
    cs = impl.Circuit().add(impl.CurrentSource(0, 1, 2e-3)).add(impl.Resistor(1, 0, 1500.0))
    assert impl.solve_dc(cs).node_voltages[1] == pytest.approx(3.0)


def test_singular_circuits_are_rejected():
    floating = impl.Circuit().add(impl.VoltageSource(1, 0, 5.0)).add(impl.Resistor(1, 0, 100.0)).add(impl.Resistor(2, 3, 100.0))   # nodes 2-3 touch nothing else
    with pytest.raises(ValueError):
        impl.solve_dc(floating)
    two_sources = impl.Circuit().add(impl.VoltageSource(1, 0, 5.0)).add(impl.VoltageSource(1, 0, 3.0)).add(impl.Resistor(1, 0, 10.0))
    with pytest.raises(ValueError):
        impl.solve_dc(two_sources)


# ---------------------------------------------------------------- meters

def test_voltmeter_loading_limits():
    c = impl.Circuit().add(impl.VoltageSource(1, 0, 12.0)).add(impl.Resistor(1, 2, 100e3)).add(impl.Resistor(2, 0, 100e3))
    assert impl.voltmeter_reading(c, 2, 0, np.inf) == pytest.approx(6.0)
    v10m = impl.voltmeter_reading(c, 2, 0, 10e6)
    assert v10m < 6.0 and v10m == pytest.approx(12.0 * impl.parallel([100e3, 10e6]) / (100e3 + impl.parallel([100e3, 10e6])))
    assert impl.voltmeter_reading(c, 2, 0, 100e3) == pytest.approx(4.0)   # a meter equal to R2 halves R2
    # the loading error scales as R_th / R_meter for a large meter resistance
    err1 = 6.0 - impl.voltmeter_reading(c, 2, 0, 1e7)
    err2 = 6.0 - impl.voltmeter_reading(c, 2, 0, 2e7)
    assert err1 / err2 == pytest.approx(2.0, rel=0.01)
    assert impl.voltmeter_reading(c, 1, 0, 1e6) == pytest.approx(12.0)   # across an ideal source: no loading


def test_ammeter_burden():
    c = impl.Circuit().add(impl.VoltageSource(1, 0, 5.0)).add(impl.Resistor(1, 0, 220.0))
    assert impl.ammeter_reading(c, 0, 0.0) == pytest.approx(5.0 / 220.0)
    assert impl.ammeter_reading(c, 0, 5.0) == pytest.approx(5.0 / 225.0)
    assert impl.ammeter_reading(c, 0, 5.0) < impl.ammeter_reading(c, 0, 0.0)
    # inserting the meter does not change the original circuit object
    assert c.resistors[0].R == 220.0 and len(c.components) == 2


# ---------------------------------------------------------------- RC

def test_rc_closed_form_matches_solve_ivp_and_the_time_constant():
    R, C, Vs = 4.7e3, 47e-6, 3.3
    tau = impl.rc_time_constant(R, C)
    assert tau == pytest.approx(R * C)
    t = np.linspace(0.0, 6 * tau, 301)
    sol = solve_ivp(impl.rc_rhs, (0.0, 6 * tau), [0.0], args=(Vs, R, C), t_eval=t, method="DOP853", rtol=1e-11, atol=1e-13)
    npt.assert_allclose(impl.rc_capacitor_voltage(t, Vs, R, C), sol.y[0], atol=1e-9)
    assert impl.rc_capacitor_voltage(tau, Vs, R, C) / Vs == pytest.approx(1 - np.exp(-1))
    assert impl.rc_capacitor_voltage(0.0, Vs, R, C) == 0.0
    assert impl.rc_capacitor_voltage(100 * tau, Vs, R, C) == pytest.approx(Vs)
    # discharge from V0 with Vs = 0, and the current is (Vs - V_C)/R
    assert impl.rc_capacitor_voltage(tau, 0.0, R, C, V0=2.0) == pytest.approx(2.0 * np.exp(-1))
    npt.assert_allclose(impl.rc_current(t, Vs, R, C), (Vs - impl.rc_capacitor_voltage(t, Vs, R, C)) / R)


def test_charging_a_capacitor_wastes_half_the_energy_whatever_the_resistor():
    C, Vs = 100e-6, 5.0
    for R in (100.0, 10e3, 1e6):
        tau = R * C
        drawn, stored, lost = impl.rc_energy_accounting(Vs, R, C, 25 * tau)
        assert drawn == pytest.approx(C * Vs ** 2, rel=1e-6)
        assert stored == pytest.approx(0.5 * C * Vs ** 2, rel=1e-6)
        assert lost == pytest.approx(0.5 * C * Vs ** 2, rel=1e-6)
        assert drawn == pytest.approx(stored + lost, rel=1e-6)


# ---------------------------------------------------------------- the bus

def test_bus_circuit_drop_and_wire_sizing():
    loads = [impl.load_resistance_for_power(P, 28.0) for P in (120.0, 60.0, 40.0)]
    assert sum(28.0 ** 2 / RL for RL in loads) == pytest.approx(220.0)
    sol0 = impl.solve_dc(impl.bus_circuit(28.0, 0.0 + 1e-9, loads))
    assert sol0.node_voltages[2] - sol0.node_voltages[3] == pytest.approx(28.0, rel=1e-6)   # no wire, no drop
    Rw = impl.wire_resistance(8.0, 1.291e-3)
    sol = impl.solve_dc(impl.bus_circuit(28.0, Rw, loads))
    V_load = sol.node_voltages[2] - sol.node_voltages[3]
    I = sol.source_currents[0]
    assert I > 0.0 and V_load < 28.0
    assert V_load == pytest.approx(28.0 - 2 * I * Rw)                                    # KVL around the harness
    assert V_load == pytest.approx(28.0 * impl.parallel(loads) / (impl.parallel(loads) + 2 * Rw))
    d = impl.wire_diameter_for_drop(220.0 / 28.0, 8.0, 28.0, 0.02)
    assert 2 * impl.wire_resistance(8.0, d) * 220.0 / 28.0 == pytest.approx(0.02 * 28.0)
    assert impl.wire_diameter_for_drop(220.0 / 28.0, 8.0, 28.0, 0.02, impl.RHO_ALUMINIUM) > d


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "dc_circuits.png").exists()
