"""Public tests for P05.18 Home Lab: Cantilever Beam Deflection. They assert properties: the cubic
satisfies its boundary conditions, is continuous with continuous slope at the load and straight beyond
it, and agrees with a second-order finite-difference solution of E I w'' = M; the moment from the cubic
balances the load (equilibrium); the self-weight formula gives qL^4/8EI; first-order uncertainty
propagation agrees with a seeded Monte Carlo; the safe load puts the clamp stress exactly at the
allowable; the weighted fit recovers a known slope and intercept from noise-free data and has a sensible
chi^2 on seeded noise; E from the slope inverts the forward model; the profile fit recovers a planted
clamp rotation; the elastica tends to linear theory for small loads and falls below it for large ones;
the SYNTHETIC files say SYNTHETIC and read back.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution.
"""
import numpy as np
import numpy.testing as npt
import pytest

from conftest import load_impl

impl = load_impl(__file__)

b, h, a, L = 0.025, 0.003, 0.290, 0.300
E = 69e9
EI = E * b * h ** 3 / 12
G = 9.81


# ---------------------------------------------------------------- theory

def test_second_moment_and_its_uncertainty():
    I, uI = impl.rect_second_moment(b, h, 0.0, 0.0)
    assert I == pytest.approx(b * h ** 3 / 12) and uI == 0.0
    I, uI = impl.rect_second_moment(b, h, 0.0, 0.01 * h)
    assert uI / I == pytest.approx(0.03)                             # thickness enters cubed
    I, uI = impl.rect_second_moment(b, h, 0.01 * b, 0.0)
    assert uI / I == pytest.approx(0.01)


def test_cubic_boundary_conditions_continuity_and_straight_tail():
    P = 5.0
    x = np.linspace(0, L, 3001)
    w = impl.deflection_point_load(x, P, a, EI)
    assert w[0] == 0.0
    # zero slope at the clamp: the one-sided difference is w'' h / 2, about 3e-4 of the slope at the load
    assert (w[1] - w[0]) / (x[1] - x[0]) < 1e-3 * P * a ** 2 / (2 * EI)
    assert impl.deflection_point_load(a, P, a, EI) == pytest.approx(P * a ** 3 / (3 * EI))
    eps = 1e-7
    left = (impl.deflection_point_load(a, P, a, EI) - impl.deflection_point_load(a - eps, P, a, EI)) / eps
    right = (impl.deflection_point_load(a + eps, P, a, EI) - impl.deflection_point_load(a, P, a, EI)) / eps
    assert left == pytest.approx(right, rel=1e-5)                                    # continuous slope
    tail = x[x > a]
    npt.assert_allclose(np.diff(impl.deflection_point_load(tail, P, a, EI), 2), 0.0, atol=1e-15)   # straight
    assert impl.deflection_point_load(L, P, a, EI) == pytest.approx(P * a ** 2 * (3 * L - a) / (6 * EI))


def test_cubic_curvature_balances_the_load_moment():
    # equilibrium of the piece beyond a cut at x: E I w''(x) = P (a - x) for x < a, 0 beyond
    P = 5.0
    x = np.linspace(0.01, L - 0.01, 50)
    hh = 1e-4
    wpp = (impl.deflection_point_load(x + hh, P, a, EI) - 2 * impl.deflection_point_load(x, P, a, EI)
           + impl.deflection_point_load(x - hh, P, a, EI)) / hh ** 2
    M = np.where(x < a, P * (a - x), 0.0)
    mask = np.abs(x - a) > 2 * hh
    npt.assert_allclose(EI * wpp[mask], M[mask], atol=1e-6 * P * a)


def test_finite_differences_converge_to_the_cubic():
    P = 6.867
    errs = []
    for n in (60, 120, 240):         # node spacing divides a = 0.29 of L = 0.3 only at n = 30 k; the kink is O(h^2) anyway
        x, w = impl.fd_cantilever(P, a, L, EI, n)
        errs.append(np.max(np.abs(w - impl.deflection_point_load(x, P, a, EI))) / (P * a ** 3 / (3 * EI)))
    assert errs[-1] < 1e-4
    assert errs[0] / errs[-1] > 10                                   # about 16 for second order


def test_self_weight_closed_form_and_linearity():
    q = 2.0
    assert impl.deflection_self_weight(L, q, L, EI) == pytest.approx(q * L ** 4 / (8 * EI))
    assert impl.deflection_self_weight(0.0, q, L, EI) == 0.0
    # superposition: the change in reading from adding P does not depend on the self-weight
    x = np.linspace(0, L, 7)
    tot = impl.deflection_self_weight(x, q, L, EI) + impl.deflection_point_load(x, 3.0, a, EI)
    npt.assert_allclose(tot - impl.deflection_self_weight(x, q, L, EI), impl.deflection_point_load(x, 3.0, a, EI))


def test_safe_load_puts_the_clamp_stress_at_the_allowable():
    Pmax = impl.max_safe_load(72.5e6, a, b, h)
    assert impl.clamp_stress(Pmax, a, b, h) == pytest.approx(72.5e6)
    # sigma = M c / I with M = P a, c = h/2
    assert impl.clamp_stress(1.0, a, b, h) == pytest.approx(a * (h / 2) / (b * h ** 3 / 12))


# ---------------------------------------------------------------- uncertainty

def test_gum_propagation_matches_monte_carlo():
    args = (4.0, a, E, 2e9, b, 0.05e-3, h, 0.02e-3, 1e-3)
    d, u = impl.predict_deflection(*args)
    dm, um = impl.predict_deflection_mc(*args, n=200_000, seed=7)
    assert d == pytest.approx(4.0 * a ** 3 / (3 * EI))
    # MC sd of a 2e5 sample is good to ~0.2 %; the first-order formula misses curvature terms of order
    # (u/x)^2 ~ 0.1 %, so 2 % is a generous band
    assert u == pytest.approx(um, rel=0.02)
    assert dm == pytest.approx(d, rel=0.005)
    d0, u0 = impl.predict_deflection(4.0, a, E, 0.0, b, 0.0, h, 0.0, 0.0)
    assert u0 == 0.0
    # reading at a fixed mark beyond the load: same propagation, a enters as a^2 (3x - a)
    d1, u1 = impl.predict_deflection(4.0, a, E, 0.0, b, 0.0, h, 0.0, 1e-3, x=L)
    dm1, um1 = impl.predict_deflection_mc(4.0, a, E, 0.0, b, 0.0, h, 0.0, 1e-3, x=L, n=200_000, seed=8)
    assert u1 == pytest.approx(um1, rel=0.02)


# ---------------------------------------------------------------- fits

def test_weighted_line_fit_recovers_exact_data():
    P = np.linspace(0, 8, 9)
    d = 2.1e-3 * P + 0.2e-3
    f = impl.fit_line(P, d, 0.4e-3)
    assert f["slope"] == pytest.approx(2.1e-3) and f["intercept"] == pytest.approx(0.2e-3, abs=1e-12)
    assert f["chi2"] == pytest.approx(0.0, abs=1e-12) and f["dof"] == 7


def test_line_fit_chi2_and_slope_uncertainty_on_seeded_noise():
    rng = np.random.default_rng(1801)
    P = np.linspace(0, 8, 9)
    slopes, chis = [], []
    for _ in range(400):
        d = 2.1e-3 * P + rng.normal(0, 0.4e-3, P.size)
        f = impl.fit_line(P, d, 0.4e-3)
        slopes.append(f["slope"]); chis.append(f["chi2"])
    # over 400 repeats the scatter of the slope matches the quoted u_slope to about 5 % (sd of a sd
    # estimate from 400 samples is 3.5 %), and the mean chi^2 is the dof to within 4 %
    assert np.std(slopes) == pytest.approx(f["u_slope"], rel=0.1)
    assert np.mean(chis) == pytest.approx(7.0, rel=0.1)


def test_modulus_from_slope_inverts_the_forward_model():
    for x in (None, L, 0.2):
        slope = float(impl.deflection_point_load(a if x is None else x, 1.0, a, EI))
        Ef, uE = impl.modulus_from_slope(slope, 0.0, a, 0.0, b, 0.0, h, 0.0, x)
        assert Ef == pytest.approx(E, rel=1e-9) and uE == pytest.approx(0.0, abs=1e-6 * E)
    slope = float(impl.deflection_point_load(a, 1.0, a, EI))
    _, uE = impl.modulus_from_slope(slope, 0.0, a, 1e-3, b, 0.0, h, 0.0)
    assert uE / E == pytest.approx(3e-3 / a, rel=1e-3)               # reading at the load: g = a^3/3
    _, uE = impl.modulus_from_slope(slope * 1.0, 0.03 * slope, a, 0.0, b, 0.0, h, 0.0)
    assert uE / E == pytest.approx(0.03, rel=1e-9)


def test_profile_fit_recovers_a_planted_clamp_rotation():
    P, theta0 = 6.867, 4e-3
    x = np.linspace(0.025, L, 12)
    w = impl.deflection_point_load(x, P, a, EI) + theta0 * x
    f = impl.fit_profile(x, w, P, a, 0.3e-3)
    assert f["EI"] == pytest.approx(EI, rel=1e-9)
    assert f["theta0"] == pytest.approx(theta0, rel=1e-9)
    assert f["chi2"] == pytest.approx(0.0, abs=1e-9)


def test_en_number():
    assert impl.en_number(70.0, 1.0, 69.0, 1.0) == pytest.approx(1 / (2 * np.sqrt(2)))
    assert impl.en_number(70.0, 1.0, 70.0, 5.0) == 0.0


# ---------------------------------------------------------------- large deflection

def test_elastica_tends_to_linear_theory_and_falls_below_it():
    small = impl.elastica_deflection(0.5, a, EI)
    assert small == pytest.approx(0.5 * a ** 3 / (3 * EI), rel=1e-3)
    big_P = 50.0
    big = impl.elastica_deflection(big_P, a, EI)
    assert big < big_P * a ** 3 / (3 * EI)
    assert big < a                                                  # cannot exceed the strip's length


# ---------------------------------------------------------------- SYNTHETIC data

def test_synthetic_dataset_is_labelled_seeded_and_reads_back(tmp_path):
    d1 = impl.make_synthetic_cantilever_data(5)
    d2 = impl.make_synthetic_cantilever_data(5)
    assert d1["loads"] == d2["loads"]
    assert max(r[0] for r in d1["loads"]) <= 2000.0                 # never above the 2 kg kit limit
    paths = impl.write_synthetic_csv(tmp_path, 5)
    for p in paths:
        assert "SYNTHETIC" in p.name
        assert "SYNTHETIC" in p.read_text().splitlines()[0]
        arr = impl.read_csv(p)
        assert arr.shape[1] == 3 and arr.shape[0] >= 7
    loads = impl.read_csv(paths[0])
    npt.assert_allclose(loads, np.array(d1["loads"]), atol=0.051)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "cantilever_lab.png").stat().st_size > 10_000
