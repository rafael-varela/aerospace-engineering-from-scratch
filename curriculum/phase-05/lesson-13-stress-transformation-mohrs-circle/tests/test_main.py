"""Public tests for P05.13 Stress Transformation & Mohr's Circle. They assert mechanics, not magic
numbers: the 2-D transformation is a rotation of the tensor (it matches R sigma R^T and preserves the
invariants); the rotated stress points lie on Mohr's circle; the principal planes carry no shear and
the principal stresses are the extreme normal stresses; the closed-form 3-D cubic agrees with the
independent numerical solver numpy.linalg.eigh on random tensors, including repeated roots;
principal directions agree with eigh's eigenvectors; the traction on any plane falls inside the
3-D Mohr region; the strain rosette inverts the strain transformation; limiting cases (pure shear,
hydrostatic, uniaxial) give their textbook answers.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

MPa = 1e6


def random_tensor(rng, scale=100e6):
    a = rng.normal(scale=scale, size=(3, 3))
    return 0.5 * (a + a.T)


# ---------------------------------------------------------------- 2-D transformation

def test_2d_transformation_equals_tensor_rotation_and_preserves_invariants():
    rng = np.random.default_rng(1301)
    for _ in range(10):
        sx, sy, txy = rng.normal(scale=80e6, size=3)
        th = rng.uniform(-np.pi, np.pi)
        a, b, c = impl.transform_stress_2d(sx, sy, txy, th)
        R = np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])
        S = R @ np.array([[sx, txy], [txy, sy]]) @ R.T                 # built here, independently
        npt.assert_allclose([a, b, c], [S[0, 0], S[1, 1], S[0, 1]], atol=1e-6)
        assert a + b == pytest.approx(sx + sy)                          # first invariant
        assert a * b - c ** 2 == pytest.approx(sx * sy - txy ** 2, rel=1e-9, abs=1e-3)   # second invariant


def test_rotated_states_lie_on_mohrs_circle_and_180_degrees_is_identity():
    sx, sy, txy = 37.5 * MPa, 75 * MPa, 30 * MPa
    C, R = impl.mohr_circle_2d(sx, sy, txy)
    th = np.linspace(0, np.pi, 37)
    a, b, c = impl.transform_stress_2d(sx, sy, txy, th)
    npt.assert_allclose(np.hypot(a - C, c), R, rtol=1e-12)
    npt.assert_allclose(np.hypot(b - C, c), R, rtol=1e-12)
    npt.assert_allclose(impl.transform_stress_2d(sx, sy, txy, np.pi), (sx, sy, txy), atol=1e-6)
    npt.assert_allclose(impl.transform_stress_2d(sx, sy, txy, np.pi / 2), (sy, sx, -txy), atol=1e-6)


def test_principal_planes_carry_no_shear_and_hold_the_extreme_normal_stresses():
    rng = np.random.default_rng(1302)
    for _ in range(10):
        sx, sy, txy = rng.normal(scale=80e6, size=3)
        s1, s2, tp = impl.principal_stresses_2d(sx, sy, txy)
        a, b, c = impl.transform_stress_2d(sx, sy, txy, tp)
        assert abs(c) < 1e-6 * max(abs(s1), abs(s2))                    # no shear on principal planes
        assert a == pytest.approx(s1) and b == pytest.approx(s2)        # the angle points at s1
        all_sn = impl.transform_stress_2d(sx, sy, txy, np.linspace(0, np.pi, 3601))[0]
        assert s2 - 1e-3 <= all_sn.min() and all_sn.max() <= s1 + 1e-3   # extremes
        ev = np.linalg.eigvalsh(np.array([[sx, txy], [txy, sy]]))      # independent check
        npt.assert_allclose([s2, s1], ev, rtol=1e-12, atol=1e-4)


def test_max_inplane_shear_plane_is_45_degrees_from_principal():
    sx, sy, txy = 77.5 * MPa, 75 * MPa, 8 * MPa
    R, ts = impl.max_inplane_shear(sx, sy, txy)
    _, _, tp = impl.principal_stresses_2d(sx, sy, txy)
    assert ts == pytest.approx(tp - np.pi / 4)
    a, b, c = impl.transform_stress_2d(sx, sy, txy, ts)
    assert abs(c) == pytest.approx(R)
    assert a == pytest.approx((sx + sy) / 2) and b == pytest.approx((sx + sy) / 2)


def test_limiting_cases_pure_shear_uniaxial_hydrostatic():
    s1, s2, tp = impl.principal_stresses_2d(0.0, 0.0, 50 * MPa)          # pure shear (torsion)
    assert s1 == pytest.approx(50 * MPa) and s2 == pytest.approx(-50 * MPa)
    assert tp == pytest.approx(np.pi / 4)
    s1, s2, tp = impl.principal_stresses_2d(100 * MPa, 0.0, 0.0)          # uniaxial
    assert (s1, s2, tp) == pytest.approx((100 * MPa, 0.0, 0.0))
    assert impl.max_inplane_shear(100 * MPa, 0.0, 0.0)[0] == pytest.approx(50 * MPa)
    C, R = impl.mohr_circle_2d(60 * MPa, 60 * MPa, 0.0)                   # equal biaxial: a point
    assert R == 0.0 and C == pytest.approx(60 * MPa)
    # but the absolute max shear of equal biaxial tension is not zero: sigma_z = 0 is a principal stress
    assert impl.plane_stress_abs_max_shear(60 * MPa, 60 * MPa, 0.0) == pytest.approx(30 * MPa)


def test_absolute_max_shear_exceeds_inplane_for_same_sign_principals():
    sx, sy, txy = impl.fuselage_skin_stress(60e3, 2.0, 1.6e-3, 0.0, 30e6)
    s1, s2, _ = impl.principal_stresses_2d(sx, sy, txy)
    assert s1 > s2 > 0
    assert impl.plane_stress_abs_max_shear(sx, sy, txy) == pytest.approx(s1 / 2)
    assert impl.plane_stress_abs_max_shear(sx, sy, txy) > impl.max_inplane_shear(sx, sy, txy)[0]
    # opposite signs: in-plane is the absolute maximum
    assert impl.plane_stress_abs_max_shear(50e6, -30e6, 0.0) == pytest.approx(impl.max_inplane_shear(50e6, -30e6, 0.0)[0])


def test_fuselage_membrane_stresses():
    sx, sy, txy = impl.fuselage_skin_stress(60e3, 2.0, 1.6e-3)
    assert sy == pytest.approx(2 * sx)                                     # hoop = 2 x axial
    assert sy == pytest.approx(60e3 * 2.0 / 1.6e-3) and txy == 0.0


# ---------------------------------------------------------------- 3-D

def test_cubic_principal_stresses_match_eigh_on_random_tensors():
    rng = np.random.default_rng(1303)
    for _ in range(200):
        S = random_tensor(rng)
        p = impl.principal_stresses_3d(S)
        ev = np.linalg.eigvalsh(S)[::-1]
        scale = np.max(np.abs(ev))
        npt.assert_allclose(p, ev, atol=1e-9 * scale)
        assert p[0] >= p[1] >= p[2]


def test_repeated_and_degenerate_roots():
    # Near a double root cos(3a) is near +-1, where arccos has infinite slope: the trigonometric
    # formula loses about half the digits, an error of order sqrt(machine eps) ~ 1e-8 of the stress
    # scale. The tolerance below (1e-7 of the scale) allows for that and no more.
    npt.assert_allclose(impl.principal_stresses_3d(-5e6 * np.eye(3)), [-5e6] * 3)          # hydrostatic
    S = impl.stress_tensor(100e6, 0.0, 0.0, 0.0, 0.0, 0.0)                                  # uniaxial: double root 0
    npt.assert_allclose(impl.principal_stresses_3d(S), [100e6, 0.0, 0.0], atol=1e-7 * 100e6)
    S = impl.stress_tensor(50e6, 50e6, -20e6, 0.0, 0.0, 0.0)                                # double root on top
    npt.assert_allclose(impl.principal_stresses_3d(S), [50e6, 50e6, -20e6], atol=1e-7 * 50e6)
    # a rotated double root is still found
    rng = np.random.default_rng(1304)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    npt.assert_allclose(impl.principal_stresses_3d(Q @ S @ Q.T), [50e6, 50e6, -20e6], atol=1e-7 * 50e6)


def test_invariants_are_unchanged_by_rotation_and_give_the_characteristic_polynomial():
    rng = np.random.default_rng(1305)
    S = random_tensor(rng)
    Q, _ = np.linalg.qr(rng.normal(size=(3, 3)))
    npt.assert_allclose(impl.invariants(impl.transform_tensor(S, Q)), impl.invariants(S), rtol=1e-9)
    I1, I2, I3 = impl.invariants(S)
    for s in np.linalg.eigvalsh(S):
        assert abs(s ** 3 - I1 * s ** 2 + I2 * s - I3) < 1e-9 * np.max(np.abs(S)) ** 3
    # rotation_z reproduces the 2-D formulas inside the 3-D tensor
    S2 = impl.stress_tensor(37.5e6, 75e6, 0.0, 30e6, 0.0, 0.0)
    Sr = impl.transform_tensor(S2, impl.rotation_z(0.4))
    npt.assert_allclose([Sr[0, 0], Sr[1, 1], Sr[0, 1]], impl.transform_stress_2d(37.5e6, 75e6, 30e6, 0.4), atol=1e-6)


def test_principal_directions_match_eigh_eigenvectors_and_are_orthogonal():
    rng = np.random.default_rng(1306)
    for _ in range(20):
        S = random_tensor(rng)
        p = impl.principal_stresses_3d(S)
        ev, V = np.linalg.eigh(S)
        N = np.array([impl.principal_direction(S, s) for s in p])
        for k in range(3):
            assert abs(N[k] @ V[:, 2 - k]) == pytest.approx(1.0, abs=1e-8)     # same axis, either sign
            npt.assert_allclose(S @ N[k], p[k] * N[k], atol=1e-6 * np.max(np.abs(p)))
        npt.assert_allclose(np.abs(N @ N.T), np.eye(3), atol=1e-8)


def test_traction_on_principal_planes_and_mohr_region():
    rng = np.random.default_rng(1307)
    S = impl.stress_tensor(180e6, -40e6, 90e6, 60e6, 30e6, -25e6)
    p = impl.principal_stresses_3d(S)
    ev, V = np.linalg.eigh(S)
    sn, tn = impl.traction_on_plane(S, V[:, -1])
    assert sn == pytest.approx(p[0]) and tn < 1e-6 * p[0]
    # every plane's (sigma_n, tau_n) lies in the region bounded by the three circles
    n = rng.normal(size=(500, 3))
    tau_max = 0.0
    for v in n:
        sn, tn = impl.traction_on_plane(S, v)
        assert impl.in_mohr_region(sn, tn, p)
        tau_max = max(tau_max, tn)
    assert tau_max <= impl.absolute_max_shear(p) * (1 + 1e-12)
    # the 45-degree plane between s1 and s3 carries the absolute max shear
    sn, tn = impl.traction_on_plane(S, V[:, -1] + V[:, 0])
    assert tn == pytest.approx(impl.absolute_max_shear(p), rel=1e-9)
    assert sn == pytest.approx((p[0] + p[2]) / 2, rel=1e-9)
    assert not impl.in_mohr_region(p[0] + 1e6, 0.0, p)                     # outside: impossible traction


# ---------------------------------------------------------------- rosette

def test_rosette_inverts_the_strain_transformation():
    rng = np.random.default_rng(1308)
    for angles in ((0.0, 45.0, 90.0), (0.0, 60.0, 120.0)):
        ex, ey, g = rng.normal(scale=1e-3, size=3)
        t = np.radians(angles)
        eps = ex * np.cos(t) ** 2 + ey * np.sin(t) ** 2 + g * np.sin(t) * np.cos(t)
        npt.assert_allclose(impl.rosette_strains(eps, angles), (ex, ey, g), rtol=1e-10, atol=1e-15)


def test_plane_stress_hooke_round_trip_and_uniaxial_limit():
    E, nu = 73.1e9, 0.33
    sx, sy, txy = 90e6, 8e6, 25e6
    G = E / (2 * (1 + nu))
    ex, ey, g = (sx - nu * sy) / E, (sy - nu * sx) / E, txy / G       # compliance, independently
    npt.assert_allclose(impl.plane_stress_from_strain(ex, ey, g, E, nu), (sx, sy, txy), rtol=1e-12)
    # uniaxial stress: ey = -nu ex gives sy = 0
    assert impl.plane_stress_from_strain(1e-3, -nu * 1e-3, 0.0, E, nu)[1] == pytest.approx(0.0, abs=1e-3)
    # principal strain and principal stress directions coincide for an isotropic material
    ex, ey, g = impl.rosette_strains(np.array([1200e-6, 900e-6, -300e-6]))
    tp_strain = 0.5 * np.arctan2(g, ex - ey)
    sx, sy, txy = impl.plane_stress_from_strain(ex, ey, g, E, nu)
    assert impl.principal_stresses_2d(sx, sy, txy)[2] == pytest.approx(tp_strain)


def test_main_writes_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "mohr_circle.png").stat().st_size > 10_000
