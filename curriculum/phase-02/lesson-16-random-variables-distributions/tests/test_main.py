"""Public tests for P02.16 Random Variables & Distributions. They assert properties, with scipy.stats as the oracle.

Run `pytest <lesson dir>` against your starter/main.py, or `AEFS_IMPL=solution pytest <lesson dir>`
against the reference solution. Hidden validation cases live in assessments/phase-02/parameterized-cases/.

Statistical tolerances: every simulation uses a fixed seed, so results are deterministic; tolerances are
nevertheless set at four standard errors (two-sided miss probability 6e-5 for a different valid generator
stream) and each one is justified in a comment where it appears.
"""
import math

import numpy as np
import numpy.testing as npt
import pytest
from scipy import stats

from conftest import load_impl

impl = load_impl(__file__)


# ---------------------------------------------------------------- Simpson and discrete moments

def test_simpson_is_exact_for_cubics_and_rejects_odd_panels():
    assert impl.simpson_integral(lambda x: x ** 3 - 2 * x + 1, -1.0, 2.0, n=2) == pytest.approx(15 / 4 - 3 + 3, rel=1e-14)
    assert impl.simpson_integral(np.sin, 0.0, math.pi, n=200) == pytest.approx(2.0, abs=1e-8)
    with pytest.raises(ValueError):
        impl.simpson_integral(np.sin, 0.0, 1.0, n=3)


@pytest.mark.parametrize("n,p", [(3, 0.98), (4, 0.9), (10, 0.3)])
def test_binomial_pmf_moments_and_scipy(n, p):
    k = np.arange(0, n + 1)
    pmf = impl.binomial_pmf(k, n, p)
    npt.assert_allclose(pmf, stats.binom.pmf(k, n, p), rtol=1e-12)
    assert pmf.sum() == pytest.approx(1.0)
    assert impl.discrete_mean(k, pmf) == pytest.approx(n * p)
    assert impl.discrete_variance(k, pmf) == pytest.approx(n * p * (1 - p))
    assert impl.binomial_pmf(n + 1, n, p) == 0.0 and impl.binomial_pmf(-1, n, p) == 0.0
    with pytest.raises(ValueError):
        impl.discrete_mean(k, pmf * 0.9)                                           # not a distribution


@pytest.mark.parametrize("lam", [0.5, 2.4, 30.0])
def test_poisson_pmf_moments_and_scipy(lam):
    k = np.arange(0, int(lam + 12 * math.sqrt(lam) + 20))
    pmf = impl.poisson_pmf(k, lam)
    npt.assert_allclose(pmf, stats.poisson.pmf(k, lam), rtol=1e-10)
    assert pmf.sum() == pytest.approx(1.0, abs=1e-10)
    assert impl.discrete_mean(k, pmf / pmf.sum()) == pytest.approx(lam, rel=1e-8)
    assert impl.discrete_variance(k, pmf / pmf.sum()) == pytest.approx(lam, rel=1e-6)
    assert impl.poisson_pmf(500, lam) >= 0.0                                        # logs: no overflow
    assert impl.poisson_pmf(0, 0.0) == 1.0


def test_binomial_approaches_poisson_as_one_over_n():
    gaps = np.array([impl.binomial_poisson_gap(n, 2.4) for n in (10, 100, 1000)])
    assert np.all(np.diff(gaps) < 0.0)
    assert gaps[1] < gaps[0] / 5 and gaps[2] < gaps[1] / 5                         # roughly 1/n


# ---------------------------------------------------------------- continuous distributions

def test_uniform_pdf_cdf_and_moments():
    a, b = -0.5, 0.5
    pdf = lambda x: impl.uniform_pdf(x, a, b)
    assert impl.simpson_integral(pdf, -1.0, 1.0, n=4000) == pytest.approx(1.0, abs=1e-3)   # step edges limit Simpson
    assert impl.continuous_mean(pdf, a, b) == pytest.approx(0.0, abs=1e-12)
    assert impl.continuous_variance(pdf, a, b) == pytest.approx(1 / 12, rel=1e-9)
    x = np.linspace(-1, 1, 11)
    npt.assert_allclose(impl.uniform_cdf(x, a, b), stats.uniform.cdf(x, a, b - a), atol=1e-14)
    assert impl.uniform_cdf(np.array(b), a, b) == 1.0 and impl.uniform_cdf(np.array(a - 1), a, b) == 0.0
    with pytest.raises(ValueError):
        impl.uniform_pdf(x, 1.0, 0.0)


def test_exponential_moments_cdf_and_memorylessness():
    lam = 1 / 5000.0
    pdf = lambda x: impl.exponential_pdf(x, lam)
    assert impl.simpson_integral(pdf, 0.0, 2e5, n=20000) == pytest.approx(1.0, abs=1e-8)
    assert impl.continuous_mean(pdf, 0.0, 2e5, n=20000) == pytest.approx(1 / lam, rel=1e-6)
    assert impl.continuous_variance(pdf, 0.0, 2e5, n=20000) == pytest.approx(1 / lam ** 2, rel=1e-5)
    t = np.array([-1.0, 0.0, 1000.0, 10000.0])
    npt.assert_allclose(impl.exponential_cdf(t, lam), stats.expon.cdf(t, scale=1 / lam), atol=1e-14)
    assert impl.exponential_pdf(np.array(-5.0), lam) == 0.0
    assert impl.memoryless_gap(lam, 3000.0, 1000.0) == pytest.approx(0.0, abs=1e-12)
    # A uniform lifetime on [0, 1] is NOT memoryless: P(T > 0.75 | T > 0.5) = 0.5 but P(T > 0.25) = 0.75.
    surv_u = lambda u: 1 - float(impl.uniform_cdf(np.array(u), 0.0, 1.0))
    assert surv_u(0.75) / surv_u(0.5) - surv_u(0.25) == pytest.approx(-0.25)


def test_normal_pdf_cdf_quantile_and_empirical_rule():
    x = np.linspace(-6, 6, 241)
    npt.assert_allclose(impl.normal_pdf(x, 1.0, 2.0), stats.norm.pdf(x, 1.0, 2.0), rtol=1e-12)
    npt.assert_allclose(impl.normal_cdf(x, 1.0, 2.0), stats.norm.cdf(x, 1.0, 2.0), atol=1e-12)
    # The cdf is the integral of the pdf.
    for xe in (-1.0, 0.3, 2.0):
        integral = impl.simpson_integral(lambda u: impl.normal_pdf(u), -12.0, xe, n=4000)
        assert integral == pytest.approx(float(impl.normal_cdf(np.array(xe))), abs=1e-9)
    # Symmetry and the empirical rule.
    assert float(impl.normal_cdf(np.array(0.0))) == pytest.approx(0.5)
    for z, frac in ((1, 0.6827), (2, 0.9545), (3, 0.9973)):
        assert 2 * float(impl.normal_cdf(np.array(float(z)))) - 1 == pytest.approx(frac, abs=5e-5)
    for p in (0.001, 0.025, 0.5, 0.84, 0.975, 0.995, 0.9999):
        q = impl.normal_quantile(p)
        assert q == pytest.approx(stats.norm.ppf(p), abs=1e-9)
        assert float(impl.normal_cdf(np.array(q))) == pytest.approx(p, abs=1e-11)   # quantile inverts the cdf
    lo, hi = impl.normal_interval(0.0, 0.02, 0.95)
    assert lo == pytest.approx(-1.959964 * 0.02, rel=1e-6) and hi == -lo
    npt.assert_allclose(impl.z_score(np.array([3.0, 5.0]), 3.0, 2.0), [0.0, 1.0])
    with pytest.raises(ValueError):
        impl.normal_quantile(1.0)
    with pytest.raises(ValueError):
        impl.normal_pdf(x, 0.0, -1.0)


# ---------------------------------------------------------------- samples and the CLT

def test_sample_statistics_match_numpy_and_the_empirical_rule_on_gaussian_data():
    rng = np.random.default_rng(1616)
    x = rng.normal(2.0, 3.0, size=50_000)
    mu, sd = impl.sample_mean_sd(x)
    assert mu == pytest.approx(x.mean()) and sd == pytest.approx(x.std(ddof=1))
    fr = impl.empirical_rule_fractions(x, 2.0, 3.0)
    # 4 SE tolerance with SE = sqrt(p(1-p)/N): 0.0021, 0.0009, 0.0002.
    for got, want, tol in zip(fr, (0.6827, 0.9545, 0.9973), (0.0085, 0.0040, 0.0010)):
        assert got == pytest.approx(want, abs=tol)
    with pytest.raises(ValueError):
        impl.sample_mean_sd(np.array([1.0]))


def test_clt_standardised_sums_approach_the_normal():
    z12 = impl.clt_sum_of_uniforms(12, 200_000, seed=16)
    assert z12.shape == (200_000,)
    # Mean within 4 SE (SE = 1/sqrt(N) = 0.0022); variance within 4 * sqrt(2/N) = 0.013 of 1.
    assert abs(z12.mean()) < 0.009
    assert z12.var(ddof=1) == pytest.approx(1.0, abs=0.013)
    frac1 = np.mean(np.abs(z12) <= 1.0)
    assert abs(frac1 - 0.6827) < 0.01                                               # 12 uniforms: within 1 % of Gaussian
    gaps = [impl.clt_distance_to_normal(n, 100_000, seed=16) for n in (1, 2, 5, 12)]
    assert np.all(np.diff(gaps) < 0.0)                                              # monotone approach
    assert gaps[0] > 0.04 and gaps[-1] < 0.01
    # n = 1 is the uniform itself: its standardised range is +/- sqrt(3), so nothing beyond 2 sigma.
    z1 = impl.clt_sum_of_uniforms(1, 1000, seed=0)
    assert np.all(np.abs(z1) <= math.sqrt(3.0) + 1e-12)
    with pytest.raises(ValueError):
        impl.clt_sum_of_uniforms(0, 10)


# ---------------------------------------------------------------- covariance and the ellipse

def test_sample_covariance_matches_numpy_and_propagation_matches_simulation():
    rng = np.random.default_rng(7)
    pts = rng.normal(size=(500, 3))
    npt.assert_allclose(impl.sample_covariance(pts), np.cov(pts, rowvar=False), rtol=1e-12)
    J, Sigma, names = impl.landing_sensitivities()
    assert J.shape == (2, 4) and Sigma.shape == (4, 4) and len(names) == 4
    cov = impl.propagate_covariance(J, Sigma)
    npt.assert_allclose(cov, cov.T)
    assert np.all(np.linalg.eigvalsh(cov) > 0.0)
    # Trace = sum of squared 1-sigma contributions (conservation of variance).
    contrib = J * np.sqrt(np.diag(Sigma))
    assert np.trace(cov) == pytest.approx(np.sum(contrib ** 2))
    landings = impl.simulate_landings(20_000, J, Sigma, seed=2016)
    assert landings.shape == (20_000, 2)
    # Sample covariance within 3 % of J Sigma J^T: SE of a variance estimate ~ sqrt(2/N) = 1 %.
    npt.assert_allclose(impl.sample_covariance(landings), cov, rtol=0.03)
    # Identity sensitivities return Sigma itself; a diagonal J scales variances by its squares.
    npt.assert_allclose(impl.propagate_covariance(np.eye(4), Sigma), Sigma)
    npt.assert_allclose(impl.propagate_covariance(np.diag([2.0, 3.0]), np.eye(2)), np.diag([4.0, 9.0]))
    with pytest.raises(ValueError):
        impl.propagate_covariance(J, np.eye(3))


def test_eig2x2_matches_eigh_and_handles_diagonal_input():
    rng = np.random.default_rng(22)
    for _ in range(10):
        B = rng.normal(size=(2, 2))
        C = B @ B.T + 0.1 * np.eye(2)
        lam, V = impl.eig2x2_symmetric(C)
        lam_np, V_np = np.linalg.eigh(C)
        npt.assert_allclose(lam, lam_np, rtol=1e-12)
        npt.assert_allclose(V.T @ V, np.eye(2), atol=1e-12)
        npt.assert_allclose(np.abs(np.sum(V * V_np, axis=0)), 1.0, atol=1e-10)      # same axes up to sign
        npt.assert_allclose(V @ np.diag(lam) @ V.T, C, atol=1e-10)
        assert lam.sum() == pytest.approx(np.trace(C))
    lam, V = impl.eig2x2_symmetric(np.diag([4.0, 1.0]))
    npt.assert_allclose(lam, [1.0, 4.0])
    npt.assert_allclose(np.abs(V), [[0, 1], [1, 0]])
    with pytest.raises(ValueError):
        impl.eig2x2_symmetric(np.array([[1.0, 2.0], [0.0, 1.0]]))


def test_ellipse_containment_matches_chi_square_and_inverts():
    for k in (0.5, 1.0, 2.0, 2.448, 3.035):
        p = impl.ellipse_containment_probability(k)
        assert p == pytest.approx(stats.chi2.cdf(k * k, 2), rel=1e-12)
        assert impl.k_for_probability(p) == pytest.approx(k, rel=1e-10)
    assert impl.ellipse_containment_probability(1.0) == pytest.approx(0.3935, abs=5e-5)
    assert impl.ellipse_containment_probability(2.0) == pytest.approx(0.8647, abs=5e-5)
    assert impl.k_for_probability(0.0) == 0.0
    with pytest.raises(ValueError):
        impl.k_for_probability(1.0)


def test_dispersion_ellipse_geometry_orientation_and_simulated_fill():
    J, Sigma, _ = impl.landing_sensitivities()
    cov = impl.propagate_covariance(J, Sigma)
    e = impl.dispersion_ellipse(cov, 0.99)
    lam = np.linalg.eigvalsh(cov)
    assert e["a"] >= e["b"] > 0.0
    assert e["a"] == pytest.approx(e["k"] * math.sqrt(lam[1])) and e["b"] == pytest.approx(e["k"] * math.sqrt(lam[0]))
    assert e["area"] == pytest.approx(math.pi * e["a"] * e["b"])
    assert -90.0 < e["angle_deg"] <= 90.0
    # Orientation: a covariance built by rotating diag(9, 1) by 30 degrees has its major axis at 30 degrees.
    th = math.radians(30.0)
    R = np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    e30 = impl.dispersion_ellipse(R @ np.diag([9.0, 1.0]) @ R.T, 0.5)
    assert e30["angle_deg"] == pytest.approx(30.0, abs=1e-9)
    assert e30["a"] / e30["b"] == pytest.approx(3.0)
    # A circular covariance gives a circle.
    ec = impl.dispersion_ellipse(4.0 * np.eye(2), 0.3935)
    assert ec["a"] == pytest.approx(ec["b"]) and ec["a"] == pytest.approx(2.0 * impl.k_for_probability(0.3935))
    # Simulated fill: fraction inside k-ellipse within 4 SE of 1 - exp(-k^2/2); N = 20000.
    landings = impl.simulate_landings(20_000, J, Sigma, seed=2016)
    for prob in (0.3935, 0.8647, 0.95, 0.99):
        k = impl.k_for_probability(prob)
        frac = impl.fraction_inside_ellipse(landings, np.zeros(2), cov, k)
        se = math.sqrt(prob * (1 - prob) / 20_000)
        assert abs(frac - prob) < 4 * se
    # Mahalanobis distance of a point on a principal axis is its sigma multiple.
    lam2, V = impl.eig2x2_symmetric(cov)
    pt = (2.0 * math.sqrt(lam2[1]) * V[:, 1])[None, :]
    assert impl.mahalanobis_squared(pt, np.zeros(2), cov)[0] == pytest.approx(4.0)
    with pytest.raises(ValueError):
        impl.dispersion_ellipse(np.array([[1.0, 2.0], [2.0, 1.0]]), 0.5)           # not positive definite


def test_circle_containment_reduces_to_closed_form_and_is_monotone():
    sigma = 469.6
    circ = sigma ** 2 * np.eye(2)
    for R in (200.0, 552.88, 1000.0, 2000.0):
        # 1e-5: the integrand has a sqrt(R^2 - u^2) end-point singularity, which limits Simpson to about 1e-6.
        assert impl.circle_containment_probability(circ, R) == pytest.approx(1 - math.exp(-R ** 2 / (2 * sigma ** 2)), abs=1e-5)
    assert impl.circle_containment_probability(circ, 0.0) == 0.0
    J, Sigma, _ = impl.landing_sensitivities()
    cov = impl.propagate_covariance(J, Sigma)
    ps = [impl.circle_containment_probability(cov, R) for R in np.linspace(50, 3000, 30)]
    assert np.all(np.diff(ps) > 0.0) and ps[-1] > 0.999
    # Against simulation, 4 SE.
    landings = impl.simulate_landings(20_000, J, Sigma, seed=2016)
    r = np.sqrt((landings ** 2).sum(axis=1))
    for R in (354.0, 800.0, 1255.0):
        p = impl.circle_containment_probability(cov, R)
        assert abs(np.mean(r <= R) - p) < 4 * math.sqrt(p * (1 - p) / 20_000)


def test_radius_for_probability_gives_cep_and_beats_the_ellipse_in_area():
    sigma = 100.0
    assert impl.radius_for_probability(sigma ** 2 * np.eye(2), 0.5) == pytest.approx(sigma * math.sqrt(2 * math.log(2)), rel=1e-4)
    assert impl.radius_for_probability(sigma ** 2 * np.eye(2), 0.99) == pytest.approx(sigma * math.sqrt(-2 * math.log(0.01)), rel=1e-4)
    J, Sigma, _ = impl.landing_sensitivities()
    cov = impl.propagate_covariance(J, Sigma)
    cep = impl.radius_for_probability(cov, 0.5)
    r99 = impl.radius_for_probability(cov, 0.99)
    assert impl.circle_containment_probability(cov, cep) == pytest.approx(0.5, abs=1e-5)
    assert cep < r99
    e99 = impl.dispersion_ellipse(cov, 0.99)
    assert e99["b"] < r99 < e99["a"]                                                # the circle sits between the axes
    assert math.pi * r99 ** 2 > e99["area"]                                         # and wastes area
    with pytest.raises(ValueError):
        impl.radius_for_probability(cov, 1.0)


def test_main_writes_the_figure(tmp_path):
    impl.main(tmp_path)
    assert (tmp_path / "landing_dispersion.png").stat().st_size > 1000
