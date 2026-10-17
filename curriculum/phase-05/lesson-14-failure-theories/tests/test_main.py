"""Public tests for P05.14 Failure Theories. They assert properties, not magic numbers: principal
stresses agree with numpy.linalg.eigh and with Mohr's circle; von Mises from components, from principal
stresses and from the distortion energy agree; both criteria vanish for hydrostatic stress, agree in
uniaxial tension and differ by exactly 2/sqrt(3) in pure shear; the Tresca/von Mises ratio stays in
[1, 2/sqrt(3)]; both are invariant under rotation and homogeneous of degree one (checked by root finding);
the tube's section stresses integrate back to the applied N, M and T (equilibrium); Coulomb-Mohr
reproduces its three closed forms and collapses to Tresca when S_t = S_c.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)


def random_tensor(rng, scale=100e6):
    a = rng.normal(scale=scale, size=(3, 3))
    return 0.5 * (a + a.T)


def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q = q @ np.diag(np.sign(np.diag(r)))
    if np.linalg.det(q) < 0:
        q[:, 0] *= -1
    return q


# ---------------------------------------------------------------- principal stresses

def test_principal_stresses_match_eigh_sorted_and_keep_the_trace():
    rng = np.random.default_rng(1401)
    for _ in range(10):
        s = random_tensor(rng)
        p = impl.principal_stresses(s)
        npt.assert_allclose(p, np.linalg.eigh(s)[0][::-1], rtol=0, atol=1e-6)   # 1e-6 Pa on 1e8 Pa
        assert p[0] >= p[1] >= p[2]
        assert p.sum() == pytest.approx(np.trace(s), abs=1e-5)


def test_mohr_circle_plane_principal_stresses_agree_with_eigh():
    rng = np.random.default_rng(1402)
    for _ in range(10):
        sx, sy, txy = rng.normal(scale=200e6, size=3)
        a, b = impl.principal_stresses_plane(sx, sy, txy)
        ev = np.linalg.eigvalsh(np.array([[sx, txy], [txy, sy]]))
        npt.assert_allclose([a, b], ev[::-1], rtol=1e-12, atol=1e-4)
        assert a + b == pytest.approx(sx + sy)                 # first invariant
        assert a * b == pytest.approx(sx * sy - txy ** 2, rel=1e-9, abs=1e3)   # second invariant


# ---------------------------------------------------------------- von Mises and Tresca

def test_von_mises_three_routes_agree():
    rng = np.random.default_rng(1403)
    for _ in range(10):
        s = random_tensor(rng)
        vm = impl.von_mises(s)
        assert vm == pytest.approx(impl.von_mises_principal(impl.principal_stresses(s)), rel=1e-10)
        dev = s - np.trace(s) / 3.0 * np.eye(3)
        J2 = 0.5 * np.sum(dev * dev)
        assert vm == pytest.approx(np.sqrt(3.0 * J2), rel=1e-10)


def test_hydrostatic_stress_causes_no_yield():
    for p in (-300e6, 50e6, 1e9):
        s = p * np.eye(3)
        assert impl.von_mises(s) == pytest.approx(0.0, abs=1e-3)
        assert impl.tresca(s) == pytest.approx(0.0, abs=1e-3)
        # adding a hydrostatic part to any state changes neither criterion
    rng = np.random.default_rng(1404)
    s = random_tensor(rng)
    assert impl.von_mises(s + 2e8 * np.eye(3)) == pytest.approx(impl.von_mises(s), rel=1e-9)
    assert impl.tresca(s + 2e8 * np.eye(3)) == pytest.approx(impl.tresca(s), rel=1e-9)


def test_uniaxial_and_pure_shear_limits():
    sig = 250e6
    s = impl.plane_stress_tensor(sig, 0.0, 0.0)
    assert impl.von_mises(s) == pytest.approx(sig)
    assert impl.tresca(s) == pytest.approx(sig)
    tau = 100e6
    s = impl.plane_stress_tensor(0.0, 0.0, tau)
    assert impl.von_mises(s) == pytest.approx(np.sqrt(3.0) * tau)
    assert impl.tresca(s) == pytest.approx(2.0 * tau)
    assert impl.max_shear_stress(s) == pytest.approx(tau)


def test_plane_stress_tresca_counts_the_zero_principal_stress():
    # equal biaxial tension: in-plane principal stresses equal, but sigma_3 = 0 out of plane
    s = impl.plane_stress_tensor(200e6, 200e6, 0.0)
    assert impl.tresca(s) == pytest.approx(200e6)
    assert impl.von_mises(s) == pytest.approx(200e6)


def test_tresca_over_von_mises_lies_between_one_and_two_over_root_three():
    rng = np.random.default_rng(1405)
    ratios = []
    for _ in range(200):
        s = random_tensor(rng)
        ratios.append(impl.tresca(s) / impl.von_mises(s))
    ratios = np.array(ratios)
    assert np.all(ratios >= 1.0 - 1e-12)
    assert np.all(ratios <= 2.0 / np.sqrt(3.0) + 1e-12)


def test_criteria_are_invariant_under_rotation_of_axes():
    rng = np.random.default_rng(1406)
    for _ in range(5):
        s = random_tensor(rng)
        Q = random_rotation(rng)
        sr = Q @ s @ Q.T
        assert impl.von_mises(sr) == pytest.approx(impl.von_mises(s), rel=1e-10)
        assert impl.tresca(sr) == pytest.approx(impl.tresca(s), rel=1e-10)


def test_distortion_energy_is_the_von_mises_energy():
    rng = np.random.default_rng(1407)
    E, nu = 205e9, 0.29
    for _ in range(5):
        s = random_tensor(rng)
        u = impl.strain_energy_density(s, E, nu)
        assert u > 0.0
        ud = impl.distortion_energy_density(s, E, nu)
        assert ud == pytest.approx((1.0 + nu) / (3.0 * E) * impl.von_mises(s) ** 2, rel=1e-9)
    # pure hydrostatic stress stores only volumetric energy
    s = 1e8 * np.eye(3)
    assert impl.distortion_energy_density(s, E, nu) == pytest.approx(0.0, abs=1e-6 * impl.strain_energy_density(s, E, nu))


def test_load_factor_by_root_finding_equals_the_homogeneous_closed_form():
    rng = np.random.default_rng(1408)
    for crit, f in (("von_mises", impl.von_mises), ("tresca", impl.tresca)):
        s = random_tensor(rng, scale=10e6)
        lam = impl.load_factor_to_yield(s, 435e6, crit)
        assert lam == pytest.approx(435e6 / f(s), rel=1e-9)
        assert f(3.0 * s) == pytest.approx(3.0 * f(s), rel=1e-12)
        assert f(-s) == pytest.approx(f(s), rel=1e-12)


# ---------------------------------------------------------------- margins

def test_margin_of_safety_definition():
    assert impl.margin_of_safety(450e6, 300e6, 1.5) == pytest.approx(0.0)
    assert impl.margin_of_safety(400e6, 200e6, 1.0) == pytest.approx(1.0)
    assert impl.margin_of_safety(400e6, -200e6, 1.0) == pytest.approx(1.0)     # sign of the stress irrelevant
    assert impl.margin_of_safety(400e6, 0.0) == np.inf
    # doubling the load more than halves the margin + 1
    assert impl.margin_of_safety(400e6, 400e6) + 1 == pytest.approx(0.5 * (impl.margin_of_safety(400e6, 200e6) + 1))


def test_margins_use_yield_at_limit_and_ultimate_at_one_and_a_half():
    mat = impl.Material("test", 70e9, 0.33, 300e6, 450e6, 2700.0)
    s = impl.plane_stress_tensor(200e6, 0.0, 0.0)
    m = impl.margins(s, mat)
    assert m["MS_yield_vm"] == pytest.approx(300 / 200 - 1)
    assert m["MS_ult_vm"] == pytest.approx(450 / (1.5 * 200) - 1)
    assert m["MS_yield_tresca"] == pytest.approx(m["MS_yield_vm"])      # uniaxial: criteria agree
    s = impl.plane_stress_tensor(0.0, 0.0, 100e6)
    m = impl.margins(s, mat)
    assert m["MS_yield_tresca"] < m["MS_yield_vm"]                      # Tresca is the conservative one
    for mt in impl.MATERIALS.values():
        assert 0.2 < mt.nu < 0.4 and mt.F_tu > mt.F_ty > 0 and mt.E > 50e9


# ---------------------------------------------------------------- the tube

def test_tube_section_properties_and_thin_wall_limit():
    s = impl.tube_section(0.05, 0.003)
    assert s["J"] == pytest.approx(2 * s["I"])
    assert s["A"] == pytest.approx(np.pi / 4 * (0.05 ** 2 - 0.044 ** 2))
    # thin wall: I -> pi r^3 t
    r, t = 0.1, 1e-5
    s = impl.tube_section(2 * r, t)
    assert s["I"] == pytest.approx(np.pi * (r - t / 2) ** 3 * t, rel=1e-6)


def test_tube_stresses_integrate_back_to_the_applied_loads():
    M, T, P = 1.2e3, 0.8e3, -20e3
    N, Mb, Tq = impl.tube_resultants(M, T, P, 0.05, 0.003)
    # midpoint rule in r and theta; the integrands are polynomial in r, so the error is O(dr^2)
    assert N == pytest.approx(P, rel=1e-4)
    assert Mb == pytest.approx(M, rel=1e-4)
    assert Tq == pytest.approx(T, rel=1e-4)


def test_critical_point_is_the_compression_side_for_compressive_axial_load():
    M, T, P = 1.2e3, 0.8e3, -20e3
    th, vm = impl.critical_point(M, T, P, 0.05, 0.003)
    assert np.cos(th) == pytest.approx(-1.0, abs=1e-6)
    s = impl.tube_section(0.05, 0.003)
    sig = P / s["A"] - M * s["c"] / s["I"]
    tau = T * s["c"] / s["J"]
    assert vm == pytest.approx(np.sqrt(sig ** 2 + 3 * tau ** 2), rel=1e-9)
    # reversing the axial load moves it to the tension side
    th2, _ = impl.critical_point(M, T, -P, 0.05, 0.003)
    assert np.cos(th2) == pytest.approx(1.0, abs=1e-6)


# ---------------------------------------------------------------- brittle criteria

def test_coulomb_mohr_closed_forms():
    St, Sc = 214e6, 752e6
    s = 50e6
    assert impl.coulomb_mohr_factor(impl.plane_stress_tensor(s, 0, 0), St, Sc) == pytest.approx(St / s)
    assert impl.coulomb_mohr_factor(impl.plane_stress_tensor(-s, 0, 0), St, Sc) == pytest.approx(Sc / s)
    assert impl.coulomb_mohr_factor(impl.plane_stress_tensor(0, 0, s), St, Sc) == pytest.approx(St * Sc / ((St + Sc) * s))
    # tension cut-off: equal biaxial tension fails at S_t
    assert impl.coulomb_mohr_factor(impl.plane_stress_tensor(s, s, 0), St, Sc) == pytest.approx(St / s)
    # in plane stress (one principal stress is zero, so sigma_1 >= 0 >= sigma_3) never more optimistic
    # than Rankine; under triaxial compression it is, by design (confinement strengthens a brittle solid)
    rng = np.random.default_rng(1409)
    for _ in range(50):
        st = impl.plane_stress_tensor(*rng.normal(scale=100e6, size=3))
        assert impl.coulomb_mohr_factor(st, St, Sc) <= impl.rankine_factor(st, St, Sc) * (1 + 1e-12)


def test_coulomb_mohr_with_equal_strengths_is_tresca():
    rng = np.random.default_rng(1410)
    S = 300e6
    for _ in range(20):
        sx, sy, txy = rng.normal(scale=100e6, size=3)
        st = impl.plane_stress_tensor(sx, sy, txy)
        assert impl.coulomb_mohr_factor(st, S, S) == pytest.approx(S / impl.tresca(st), rel=1e-10)


def test_envelopes_lie_on_their_criteria():
    F = 435e6
    sig, tau = impl.sigma_tau_envelope("von_mises", F, n=73)
    npt.assert_allclose(sig ** 2 + 3 * tau ** 2, F ** 2, rtol=1e-10)
    sig, tau = impl.sigma_tau_envelope("tresca", F, n=73)
    npt.assert_allclose(sig ** 2 + 4 * tau ** 2, F ** 2, rtol=1e-10)
    sig, tau = impl.sigma_tau_envelope("coulomb_mohr", 214e6, n=73, S_c=752e6)
    assert np.max(sig) == pytest.approx(214e6, rel=1e-9) and np.min(sig) == pytest.approx(-752e6, rel=1e-9)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "failure_theories.png").stat().st_size > 10_000
