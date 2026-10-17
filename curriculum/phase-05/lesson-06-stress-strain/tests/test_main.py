"""Public tests for P05.06 Stress & Strain. They assert physics and properties, not magic numbers:
engineering and true measures obey constant volume; the 3-D Hooke's law is invertible, gives the
Poisson contraction, G = E / 2(1 + nu) in pure shear and the bulk modulus under pressure; the
analyzer recovers the 0.2 % offset yield of a Ramberg-Osgood curve exactly (the independent
check: the yield is known by construction), the modulus of a straight line with a toe, the
ultimate point at Considere's criterion (computed independently here by root finding on the true
curve), and the toughness of an elastic-perfectly-plastic curve in closed form; results scale
correctly with stress; the offset construction is immune to a linear machine compliance; seeded
noise moves the results by statistically bounded amounts.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest
from scipy.optimize import brentq

from conftest import load_impl

impl = load_impl(__file__)

E_AL, SY_AL, N_RO = 73.1e3, 345.0, 9.0
COUPON = dict(E=E_AL, sigma_y=SY_AL, n=N_RO, eps_fracture=0.18, drop=0.10)


# ---------------------------------------------------------------- stress and strain measures

def test_engineering_measures_are_linear_and_dimensionally_consistent():
    F = np.array([0.0, 10e3, 20e3])
    s = impl.engineering_stress(F, 75.0)                     # N / mm^2 = MPa
    npt.assert_allclose(s, F / 75.0)
    npt.assert_allclose(impl.engineering_stress(2 * F, 75.0), 2 * s)
    npt.assert_allclose(impl.engineering_strain(np.array([0.0, 0.5]), 50.0), [0.0, 0.01])


def test_true_stress_strain_follow_constant_volume():
    rng = np.random.default_rng(6)
    eps = rng.uniform(0, 0.3, 50)
    sig = rng.uniform(1, 600, 50)
    st, et = impl.true_stress_strain(sig, eps)
    npt.assert_allclose(et, np.log1p(eps))
    # constant volume: A = A0 L0 / L, so F / A = sigma (1 + eps) = sigma exp(eps_t)
    npt.assert_allclose(st, sig * np.exp(et))
    assert np.all(st >= sig) and np.all(et <= eps)
    # small-strain limit: true and engineering agree to first order
    st0, et0 = impl.true_stress_strain(np.array([300.0]), np.array([1e-5]))
    assert st0[0] == pytest.approx(300.0, rel=2e-5) and et0[0] == pytest.approx(1e-5, rel=1e-5)


# ---------------------------------------------------------------- Hooke's law

def test_uniaxial_hooke_gives_poisson_contraction():
    E, nu, s = 73.1e3, 0.33, 200.0
    e = impl.hooke_strain(np.diag([s, 0.0, 0.0]), E, nu)
    assert e[0, 0] == pytest.approx(s / E)
    assert e[1, 1] == pytest.approx(-nu * s / E) and e[2, 2] == pytest.approx(-nu * s / E)
    npt.assert_allclose(e - np.diag(np.diag(e)), 0.0, atol=1e-18)


def test_hooke_stress_inverts_hooke_strain_and_is_symmetric():
    rng = np.random.default_rng(66)
    for E, nu in ((73.1e3, 0.33), (113.8e3, 0.342), (205e3, 0.29)):
        for _ in range(5):
            a = rng.normal(scale=100.0, size=(3, 3))
            S = 0.5 * (a + a.T)
            e = impl.hooke_strain(S, E, nu)
            npt.assert_allclose(e, e.T, atol=1e-18)
            npt.assert_allclose(impl.hooke_stress(e, E, nu), S, atol=1e-9)


def test_pure_shear_and_hydrostatic_pressure_give_G_and_K():
    E, nu, tau, p = 71.7e3, 0.33, 50.0, 100.0
    S = np.zeros((3, 3)); S[0, 1] = S[1, 0] = tau
    e = impl.hooke_strain(S, E, nu)
    gamma = 2 * e[0, 1]                                      # engineering shear strain
    assert tau / gamma == pytest.approx(impl.shear_modulus(E, nu))
    assert impl.shear_modulus(E, nu) == pytest.approx(E / (2 * (1 + nu)))
    e = impl.hooke_strain(-p * np.eye(3), E, nu)
    assert -p / np.trace(e) == pytest.approx(impl.bulk_modulus(E, nu))
    # incompressible limit: no volume change under any stress
    a = np.random.default_rng(1).normal(size=(3, 3))
    assert np.trace(impl.hooke_strain(a + a.T, E, 0.5)) == pytest.approx(0.0, abs=1e-15)


def test_material_table_is_physically_consistent():
    for name, m in impl.MATERIALS.items():
        assert 60e3 < m["E"] < 220e3, name
        assert 0.25 < m["nu"] < 0.36, name
        assert m["Fty"] < m["Ftu"], name                    # yield below ultimate
        assert 0.0 < m["elongation"] < 0.4, name
        assert impl.bulk_modulus(m["E"], m["nu"]) > impl.shear_modulus(m["E"], m["nu"]) > 0


# ---------------------------------------------------------------- the analyzer: independent checks

def test_offset_yield_of_a_ramberg_osgood_curve_is_its_sigma_y():
    for E, sy, n in ((E_AL, SY_AL, 9.0), (113.8e3, 880.0, 25.0), (205e3, 435.0, 6.0)):
        s = np.linspace(0, 1.3 * sy, 20001)
        e = impl.ramberg_osgood_strain(s, E, sy, n)
        assert impl.ramberg_osgood_strain(sy, E, sy, n) == pytest.approx(sy / E + 0.002)
        ey, syf = impl.offset_yield(e, s, E)
        assert syf == pytest.approx(sy, rel=1e-6)
        assert ey == pytest.approx(sy / E + 0.002, rel=1e-6)


def test_fit_modulus_recovers_slope_and_toe_of_a_line():
    eps0, E = 4e-4, 71.7e3
    e = np.linspace(eps0, 0.01, 400)
    s = E * (e - eps0)
    s = np.minimum(s, 500.0)                                   # flat top so there is a maximum
    Ef, e0f = impl.fit_modulus(e, s)
    assert Ef == pytest.approx(E, rel=1e-9) and e0f == pytest.approx(eps0, rel=1e-6)


def test_ultimate_point_is_at_considere_criterion():
    E, sy, n = E_AL, SY_AL, N_RO
    # independent: d(sigma_t)/d(eps_t) = sigma_t on the true Ramberg-Osgood curve, by root finding
    g = lambda s: 1.0 / (1.0 / E + 0.002 * n * s ** (n - 1) / sy ** n) - s
    s_c = brentq(g, sy, 3 * sy)
    et = s_c / E + 0.002 * (s_c / sy) ** n
    eps_u, sig_u = np.expm1(et), s_c * np.exp(-et)
    a, b = impl.considere_point(E, sy, n)
    assert a == pytest.approx(eps_u, rel=1e-8) and b == pytest.approx(sig_u, rel=1e-10)
    x, y = impl.synthetic_tensile_test(**COUPON)
    eu, su = impl.ultimate_strength(x, y)
    assert su == pytest.approx(sig_u, rel=1e-5)
    assert eu == pytest.approx(eps_u, abs=2 * 0.18 / 1999)      # within one strain step


def test_clean_synthetic_coupon_reduces_to_its_inputs():
    x, y = impl.synthetic_tensile_test(**COUPON)
    r = impl.analyze_tensile(x, y)
    assert r.E == pytest.approx(E_AL, rel=0.01)                 # RO curvature in the window: < 1 %
    assert r.eps0 == pytest.approx(0.0, abs=2e-5)
    # engineering vs true stress differ by (1 + eps) ~ 0.7 % at yield: the reduction is within 1 %
    assert r.sigma_y == pytest.approx(SY_AL, rel=0.01) and r.sigma_y < SY_AL
    sig_f = 0.9 * r.sigma_u
    assert y[-1] == pytest.approx(sig_f, rel=1e-6)          # sampled maximum: 1e-6
    assert r.elongation == pytest.approx(0.18 - sig_f / r.E - r.eps0, rel=1e-6)
    assert r.resilience == pytest.approx(r.sigma_y ** 2 / (2 * r.E))
    assert r.sigma_y < r.sigma_u and r.eps_y < r.eps_u < 0.18


def test_toughness_of_elastic_perfectly_plastic_curve():
    E, sy, ef = 205e3, 435.0, 0.25
    ey = sy / E
    e = np.concatenate([np.linspace(0, ey, 200), np.linspace(ey, ef, 2000)[1:]])
    s = np.minimum(E * e, sy)
    assert impl.toughness(e, s) == pytest.approx(sy * (ef - ey / 2), rel=1e-12)   # exact for piecewise linear
    assert impl.toughness(e[e <= ey], s[e <= ey]) == pytest.approx(impl.resilience(sy, E), rel=1e-12)


def test_toughness_trapezoid_converges_at_second_order():
    E, sy, n = E_AL, SY_AL, N_RO
    s_ref = np.linspace(0, 450.0, 200001)
    ref = impl.toughness(impl.ramberg_osgood_strain(s_ref, E, sy, n), s_ref)
    errs = []
    for m in (101, 201, 401):
        s = np.linspace(0, 450.0, m)
        errs.append(abs(impl.toughness(impl.ramberg_osgood_strain(s, E, sy, n), s) - ref))
    orders = np.log(np.array(errs[:-1]) / np.array(errs[1:])) / np.log(2)
    assert np.all(orders > 1.8) and np.all(orders < 2.2)


def test_results_scale_with_stress_and_strain_units():
    x, y = impl.synthetic_tensile_test(**COUPON)
    r1 = impl.analyze_tensile(x, y)
    r2 = impl.analyze_tensile(x, 2.5 * y)                      # e.g. a stronger material, same shape
    assert r2.E == pytest.approx(2.5 * r1.E) and r2.sigma_y == pytest.approx(2.5 * r1.sigma_y)
    assert r2.sigma_u == pytest.approx(2.5 * r1.sigma_u) and r2.toughness == pytest.approx(2.5 * r1.toughness)
    assert r2.elongation == pytest.approx(r1.elongation) and r2.eps_y == pytest.approx(r1.eps_y)


def test_vanishing_offset_and_monotone_order_of_points():
    s = np.linspace(0, 450.0, 20001)
    e = impl.ramberg_osgood_strain(s, E_AL, SY_AL, N_RO)
    ys = [impl.offset_yield(e, s, E_AL, offset=o)[1] for o in (0.0005, 0.001, 0.002, 0.005)]
    assert np.all(np.diff(ys) > 0)                             # larger offset, higher proof stress
    # offset o gives sigma_y (o / 0.002)^(1/n) on a Ramberg-Osgood curve
    assert ys[0] == pytest.approx(SY_AL * 0.25 ** (1 / N_RO), rel=1e-6)


# ---------------------------------------------------------------- measurement imperfections

def test_machine_compliance_lowers_apparent_modulus_but_not_offset_yield():
    c = 5e-6
    x0, y0 = impl.synthetic_tensile_test(**COUPON)
    x, y = impl.synthetic_tensile_test(**COUPON, compliance=c)
    r0, r = impl.analyze_tensile(x0, y0), impl.analyze_tensile(x, y)
    assert r.E == pytest.approx(1.0 / (1.0 / r0.E + c), rel=1e-5)   # springs in series (fit on a curve: 1e-5)
    assert r.sigma_y == pytest.approx(r0.sigma_y, rel=1e-6)          # offset line moves with the data
    assert r.elongation == pytest.approx(r0.elongation, rel=1e-6)
    rc = impl.analyze_tensile(impl.remove_compliance(x, y, c), y)
    assert rc.E == pytest.approx(r0.E, rel=1e-9)
    # but the handbook modulus with crosshead strain under-reads the yield
    assert impl.offset_yield(x, y, E_AL, r.eps0)[1] < 0.9 * r0.sigma_y


def test_seeded_noise_moves_results_by_bounded_amounts():
    clean = impl.analyze_tensile(*impl.synthetic_tensile_test(**COUPON))
    Es, Ys, Us = [], [], []
    for seed in range(40):
        r = impl.analyze_tensile(*impl.synthetic_tensile_test(**COUPON, noise_stress=1.0, noise_strain=2e-5,
                                                              seed=seed), window=21)
        Es.append(r.E); Ys.append(r.sigma_y); Us.append(r.sigma_u)
    # single-run spread of E is about 0.6 GPa and of Fty about 1.2 MPa (lesson table, 200 seeds);
    # the mean of 40 runs is then within 0.1 GPa and 0.2 MPa at 1 sigma; allow about 4 sigma.
    assert abs(np.mean(Es) - clean.E) < 0.4e3
    assert abs(np.mean(Ys) - clean.sigma_y) < 1.0
    # a 21-point average shrinks the 1 MPa noise to ~0.2 MPa; the max of a flat top is biased up a little
    assert 0.0 < np.mean(Us) - clean.sigma_u < 1.0
    # same seed, same data
    a = impl.synthetic_tensile_test(**COUPON, noise_stress=1.0, seed=3)
    b = impl.synthetic_tensile_test(**COUPON, noise_stress=1.0, seed=3)
    npt.assert_array_equal(a[1], b[1])


def test_generator_rejects_fracture_before_necking():
    with pytest.raises(ValueError):
        impl.synthetic_tensile_test(E=E_AL, sigma_y=SY_AL, n=N_RO, eps_fracture=0.05)


def test_main_writes_figure(tmp_path, capsys):
    impl.main(tmp_path)
    assert (tmp_path / "stress_strain.png").stat().st_size > 10_000
    assert "SYNTHETIC" in capsys.readouterr().out
