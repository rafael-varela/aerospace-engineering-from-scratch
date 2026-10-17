"""Public tests for P05.07 Axial Loading & Thermal Stress. They assert physics, not magic numbers:
the bar finite-element solver satisfies equilibrium (reactions balance the loads) and reproduces
the closed forms (stepped bar, fixed-fixed bar with load and temperature, constrained thermal
stress) to round-off and the tapered bar at second order; the pressure-cone member stiffness
matches a scipy.integrate.quad integration of the cone and the single-frustum closed form; the
joint diagram obeys equilibrium (F_b - F_m = P) and compatibility (bolt stretch = member
unloading), agrees with an independent two-element FE model of the joint with a misfit, has the
right limits (rigid members, rigid bolt, separation); the thermal change of preload vanishes for
matched materials and is linear in dT.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.integrate import quad

from conftest import load_impl

impl = load_impl(__file__)

E_AL, E_TI, E_ST = 71.7e3, 113.8e3, 205e3
A_AL, A_TI = 23.6e-6, 8.6e-6


# ---------------------------------------------------------------- bars

def test_axial_deformation_and_stepped_bar():
    assert impl.axial_deformation(15e3, 1200.0, 100.0, 73.1e3) == pytest.approx(15e3 * 1200 / (73.1e3 * 100))
    # two equal halves in series give the same as one bar
    d = impl.stepped_bar_deformation([5e3, 5e3], [300.0, 300.0], [80.0, 80.0], [E_ST, E_ST])
    assert d == pytest.approx(impl.axial_deformation(5e3, 600.0, 80.0, E_ST))
    # linear in load, and a compressive segment shortens
    assert impl.stepped_bar_deformation([1e3, -2e3], [100.0, 100.0], [10.0, 10.0], [E_AL, E_AL]) < 0


def test_fe_reproduces_stepped_bar_and_satisfies_equilibrium():
    # free end loaded bar with three segments, intermediate loads: N_i known by sections
    x = np.array([0.0, 200.0, 450.0, 600.0])
    EA = np.array([E_ST * 80, E_AL * 120, E_TI * 60])
    loads = np.array([0.0, -3e3, 5e3, 8e3])
    u, N, R = impl.bar_fe_solve(x, None, EA, loads, [0])
    npt.assert_allclose(N, [10e3, 13e3, 8e3], rtol=1e-12)          # method of sections
    assert R[0] == pytest.approx(-10e3) and R.sum() + loads.sum() == pytest.approx(0.0, abs=1e-9)
    assert u[-1] == pytest.approx(impl.stepped_bar_deformation(N, np.diff(x), EA / [E_ST, E_AL, E_TI],
                                                               [E_ST, E_AL, E_TI]), rel=1e-12)


def test_fixed_fixed_bar_matches_fe_with_load_and_temperature():
    for P, dT in ((10e3, 0.0), (0.0, 40.0), (10e3, 40.0), (-7e3, -60.0)):
        args = (P, 300.0, E_ST * 80, 12.2e-6, 400.0, 73.1e3 * 150, 23.2e-6)
        u, N1, N2 = impl.fixed_fixed_two_segment(*args, dT=dT)
        uf, Nf, Rf = impl.bar_fe_solve([0.0, 300.0, 700.0], None, [E_ST * 80, 73.1e3 * 150], [0.0, P, 0.0], [0, 2],
                                       eps0=[12.2e-6 * dT, 23.2e-6 * dT])
        assert uf[1] == pytest.approx(u, abs=1e-12) and Nf == pytest.approx([N1, N2], abs=1e-8)
        assert N2 - N1 + P == pytest.approx(0.0, abs=1e-8)            # junction equilibrium
        assert Rf.sum() + P == pytest.approx(0.0, abs=1e-8)            # overall equilibrium
        # compatibility: the two elongations (mechanical + thermal) cancel
        e1 = N1 * 300 / (E_ST * 80) + 12.2e-6 * dT * 300
        e2 = N2 * 400 / (73.1e3 * 150) + 23.2e-6 * dT * 400
        assert e1 + e2 == pytest.approx(0.0, abs=1e-12)


def test_superposition_of_load_and_temperature():
    base = (300.0, E_ST * 80, 12.2e-6, 400.0, 73.1e3 * 150, 23.2e-6)
    a = np.array(impl.fixed_fixed_two_segment(10e3, *base, dT=0.0))
    b = np.array(impl.fixed_fixed_two_segment(0.0, *base, dT=40.0))
    c = np.array(impl.fixed_fixed_two_segment(10e3, *base, dT=40.0))
    npt.assert_allclose(a + b, c, rtol=1e-12)


def test_constrained_bar_thermal_stress_and_its_fe_limit():
    E, al, dT = 73.1e3, 23.2e-6, 50.0
    s = impl.constrained_thermal_stress(E, al, dT)
    assert s == pytest.approx(-E * al * dT) and s < 0
    # uniform bar between walls, any number of elements: no displacement, stress -E alpha dT
    u, N, R = impl.bar_fe_solve(np.linspace(0, 1000.0, 6), None, E * 100.0, np.zeros(6), [0, 5], eps0=al * dT)
    npt.assert_allclose(u, 0.0, atol=1e-12)
    npt.assert_allclose(N / 100.0, s, rtol=1e-12)
    # the stress does not depend on length or area
    u2, N2, _ = impl.bar_fe_solve(np.linspace(0, 37.0, 3), None, E * 7.0, np.zeros(3), [0, 2], eps0=al * dT)
    npt.assert_allclose(N2 / 7.0, s, rtol=1e-12)
    # a free bar: no stress, free expansion alpha dT L
    u3, N3, _ = impl.bar_fe_solve(np.linspace(0, 1000.0, 6), None, E * 100.0, np.zeros(6), [0], eps0=al * dT)
    npt.assert_allclose(N3, 0.0, atol=1e-9)
    assert u3[-1] == pytest.approx(al * dT * 1000.0)


def test_tapered_bar_fe_converges_at_second_order():
    P, L, A0, A1, E = 20e3, 500.0, 200.0, 50.0, 73.1e3
    exact = impl.tapered_bar_deformation(P, L, A0, A1, E)
    assert exact == pytest.approx(quad(lambda x: P / (E * (A0 + (A1 - A0) * x / L)), 0, L)[0], rel=1e-10)
    errs = []
    for m in (8, 16, 32, 64):
        x = np.linspace(0, L, m + 1)
        xm = 0.5 * (x[1:] + x[:-1])
        u, _, _ = impl.bar_fe_solve(x, None, E * (A0 + (A1 - A0) * xm / L), np.r_[np.zeros(m), P], [0])
        errs.append(abs(u[-1] - exact) / exact)
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.9) and np.all(orders < 2.1)
    # uniform limit
    assert impl.tapered_bar_deformation(P, L, A0, A0, E) == pytest.approx(P * L / (E * A0))


# ---------------------------------------------------------------- bolt and members

def test_tensile_stress_area_and_bolt_stiffness():
    At = impl.tensile_stress_area(8.0, 1.25)
    assert 36.0 < At < 37.0                                      # M8 x 1.25: 36.6 mm^2
    assert At < np.pi * 8.0 ** 2 / 4
    kb_all_shank = impl.bolt_stiffness(8.0, 1.25, 20.0, 0.0, E_TI)
    assert kb_all_shank == pytest.approx(E_TI * np.pi * 16.0 / 20.0)
    kb = impl.bolt_stiffness(8.0, 1.25, 14.0, 6.0, E_TI)
    assert impl.bolt_stiffness(8.0, 1.25, 0.0, 20.0, E_TI) < kb < kb_all_shank   # threads are softer
    assert impl.bolt_stiffness(8.0, 1.25, 28.0, 12.0, E_TI) == pytest.approx(kb / 2)  # twice the grip


def test_frustum_stiffness_matches_quadrature_of_the_cone():
    E, d, Dw, t = E_AL, 8.0, 12.0, 10.0
    ta = np.tan(np.radians(30.0))
    comp = quad(lambda s: 1.0 / (E * np.pi * ((Dw / 2 + s * ta) ** 2 - (d / 2) ** 2)), 0, t)[0]
    assert impl.frustum_stiffness(E, d, Dw, t) == pytest.approx(1.0 / comp, rel=1e-10)
    # thin-plate limit: an annulus of area pi (Dw^2 - d^2)/4 and thickness t
    tt = 1e-4
    assert impl.frustum_stiffness(E, d, Dw, tt) == pytest.approx(E * np.pi * (Dw ** 2 - d ** 2) / 4 / tt, rel=1e-4)


def test_member_stiffness_symmetric_and_layered():
    d, Dw = 8.0, 12.0
    km = impl.member_stiffness(d, Dw, [(10.0, E_AL), (10.0, E_AL)])
    assert km == pytest.approx(impl.frustum_stiffness(E_AL, d, Dw, 10.0) / 2, rel=1e-12)
    # splitting a layer changes nothing; the order of layers on one side does
    assert impl.member_stiffness(d, Dw, [(5.0, E_AL)] * 4) == pytest.approx(km, rel=1e-12)
    # layered stack against an independent quadrature over the double cone
    layers = [(4.0, E_ST), (9.0, E_AL), (7.0, E_TI)]
    l = 20.0
    ta = np.tan(np.radians(30.0))
    def E_at(x):
        z = 0.0
        for t, E in layers:
            if x <= z + t:
                return E
            z += t
    integrand = lambda x: 1.0 / (E_at(x) * np.pi * (((Dw + 2 * min(x, l - x) * ta) / 2) ** 2 - (d / 2) ** 2))
    comp = quad(integrand, 0, l, points=[4.0, 10.0, 13.0], epsabs=1e-16, epsrel=1e-12, limit=200)[0]
    assert impl.member_stiffness(d, Dw, layers) == pytest.approx(1.0 / comp, rel=1e-8)
    # stiffer material, stiffer joint; proportional to E for a single material
    assert impl.member_stiffness(d, Dw, [(20.0, 2 * E_AL)]) == pytest.approx(2 * km, rel=1e-12)


# ---------------------------------------------------------------- the joint diagram

def test_joint_loads_obey_equilibrium_and_compatibility():
    kb, km, Fi = 257e3, 577e3, 20e3
    C = impl.joint_constant(kb, km)
    P = np.linspace(0, 0.99 * impl.separation_load(Fi, kb, km), 50)
    Fb, Fm = impl.joint_loads(Fi, P, kb, km)
    npt.assert_allclose(Fb - Fm, P, atol=1e-8)                             # equilibrium
    npt.assert_allclose((Fb - Fi) / kb, (Fi - Fm) / km, rtol=1e-12, atol=1e-15)   # compatibility
    npt.assert_allclose(np.gradient(Fb, P), C, rtol=1e-10)
    assert 0 < C < 1


def test_joint_separation_and_limits():
    kb, km, Fi = 257e3, 577e3, 20e3
    Psep = impl.separation_load(Fi, kb, km)
    Fb, Fm = impl.joint_loads(Fi, np.array([Psep * 0.999, Psep * 1.001, 2 * Psep]), kb, km)
    assert Fm[0] > 0 and Fm[1] == 0 and Fm[2] == 0
    npt.assert_allclose(Fb[1:], [Psep * 1.001, 2 * Psep])                  # bolt carries everything
    assert Fb[0] == pytest.approx(Psep * 0.999, rel=2e-3)                  # continuous at separation
    # rigid members: the bolt feels nothing until separation; rigid bolt: the bolt takes it all
    Fb_r, _ = impl.joint_loads(Fi, 5e3, kb, 1e12)
    assert float(Fb_r) == pytest.approx(Fi, rel=1e-6)
    assert impl.joint_constant(1e12, km) == pytest.approx(1.0, rel=1e-6)
    # no external load: both carry the preload
    Fb0, Fm0 = impl.joint_loads(Fi, 0.0, kb, km)
    assert float(Fb0) == Fi and float(Fm0) == Fi


def test_joint_diagram_against_two_element_fe_model_with_misfit():
    J = impl.LESSON_JOINT
    kb, km, l = J.k_b, J.k_m, J.grip
    misfit = J.F_i * (1 / kb + 1 / km)          # the nut turn that produces the preload
    for P in (0.0, 5e3, 12e3, 0.95 * impl.separation_load(J.F_i, kb, km)):
        u, N, R = impl.bar_fe_solve(np.array([0.0, l]), [[0, 1], [0, 1]], [kb * l, km * l], [0.0, P], [0],
                                    eps0=[-misfit / l, 0.0])
        Fb, Fm = impl.joint_loads(J.F_i, P, kb, km)
        assert N[0] == pytest.approx(float(Fb), rel=1e-10) and -N[1] == pytest.approx(float(Fm), rel=1e-9, abs=1e-6)
        assert R[0] + P == pytest.approx(0.0, abs=1e-8)


def test_thermal_preload_change():
    kb, km = 257e3, 577e3
    layers_al = [(10.0, A_AL), (10.0, A_AL)]
    # matched materials: no change
    assert impl.thermal_preload_change(-75.0, A_AL, layers_al, kb, km) == pytest.approx(0.0, abs=1e-12)
    dF_cold = impl.thermal_preload_change(-75.0, A_TI, layers_al, kb, km)
    dF_hot = impl.thermal_preload_change(50.0, A_TI, layers_al, kb, km)
    assert dF_cold < 0 < dF_hot                                  # aluminium shrinks more than titanium
    assert dF_hot / dF_cold == pytest.approx(50.0 / -75.0)       # linear in dT
    # independent check: two parallel FE elements with thermal initial strains
    l = 20.0
    u, N, _ = impl.bar_fe_solve(np.array([0.0, l]), [[0, 1], [0, 1]], [kb * l, km * l], [0.0, 0.0], [0],
                                eps0=[A_TI * -75.0, A_AL * -75.0])
    assert N[0] == pytest.approx(dF_cold, rel=1e-10) and N[1] == pytest.approx(-dF_cold, rel=1e-10)


def test_lesson_joint_is_in_a_realistic_range():
    J = impl.LESSON_JOINT
    C = impl.joint_constant(J.k_b, J.k_m)
    assert 0.15 < C < 0.5                                        # members stiffer than the bolt
    assert 0.4 < J.F_i / J.A_t / impl.MATERIALS["Ti-6Al-4V"]["Fty"] < 0.8
    assert J.preload_at(impl.DT_COLD) < J.F_i < J.preload_at(impl.DT_HOT)
    assert impl.preload_from_torque(0.2 * J.F_i * J.d, 0.2, J.d) == pytest.approx(J.F_i)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "bolted_joint.png").stat().st_size > 10_000
