"""P02.16 Random Variables & Distributions: reference solution.

Discrete and continuous random variables, their pmf/pdf/cdf, expectation and variance, the
binomial, Poisson, uniform, exponential and normal distributions written from scratch (the one
borrowed primitive is the error function), the central limit theorem shown by simulation, and
the propagation of input uncertainties through a linear sensitivity matrix to a two-dimensional
landing dispersion whose constant-probability contours are ellipses. ``scipy.stats`` is the oracle
in the tests, never the implementation.

All landing data in this module are SYNTHETIC: they are generated from the stated sensitivity
matrix and input standard deviations, not measured.

Conventions:
    * Angles returned by ``dispersion_ellipse`` are in degrees, measured from the +x (down-range)
      axis towards +y (cross-range).
    * Every random function takes ``seed`` and builds ``numpy.random.default_rng(seed)`` itself.
    * Covariance matrices are 2x2 symmetric positive definite unless a function says otherwise.

Run:  python solution.py          (prints the tables, writes landing_dispersion.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-16-random-variables-distributions
"""
from __future__ import annotations

import math
from collections.abc import Callable
from pathlib import Path

import numpy as np
from scipy.special import erf

# ---------------------------------------------------------------- numerical integration (from scratch)

def simpson_integral(f: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int = 2000) -> float:
    """Return the composite Simpson's-rule estimate of the integral of ``f`` over [a, b] with ``n`` panels.

    ``n`` must be even and at least 2. Exact for cubics; error O(h^4). Raise ``ValueError`` otherwise.
    """
    if n < 2 or n % 2:
        raise ValueError("n must be an even integer >= 2")
    x = np.linspace(a, b, n + 1)
    y = np.asarray(f(x), dtype=float)
    h = (b - a) / n
    return float(h / 3.0 * (y[0] + y[-1] + 4.0 * y[1:-1:2].sum() + 2.0 * y[2:-1:2].sum()))


# ---------------------------------------------------------------- discrete random variables

def discrete_mean(values: np.ndarray, probs: np.ndarray) -> float:
    """Return E[X] = sum x P(x). Raise ``ValueError`` unless the probabilities sum to 1 and are non-negative."""
    values = np.asarray(values, dtype=float)
    probs = np.asarray(probs, dtype=float)
    if values.shape != probs.shape or np.any(probs < 0.0) or abs(probs.sum() - 1.0) > 1e-9:
        raise ValueError("probs must be a probability distribution matching values")
    return float(np.dot(values, probs))


def discrete_variance(values: np.ndarray, probs: np.ndarray) -> float:
    """Return Var X = E[(X - mu)^2] = sum (x - mu)^2 P(x)."""
    values = np.asarray(values, dtype=float)
    mu = discrete_mean(values, probs)
    return float(np.dot((values - mu) ** 2, np.asarray(probs, dtype=float)))


def binomial_pmf(k: np.ndarray | int, n: int, p: float) -> np.ndarray | float:
    """Return P(X = k) = C(n, k) p^k (1-p)^(n-k) for the number of successes in n Bernoulli trials.

    Vectorised over ``k``; zero outside 0..n. Uses ``math.comb`` (exact integers).
    """
    if not 0.0 <= p <= 1.0 or n < 0:
        raise ValueError("need 0 <= p <= 1 and n >= 0")
    k_arr = np.atleast_1d(np.asarray(k, dtype=int))
    out = np.zeros(k_arr.shape, dtype=float)
    for i, kk in enumerate(k_arr.ravel()):
        if 0 <= kk <= n:
            out.ravel()[i] = math.comb(n, int(kk)) * p ** kk * (1.0 - p) ** (n - kk)
    return out if np.ndim(k) else float(out[0])


def poisson_pmf(k: np.ndarray | int, lam: float) -> np.ndarray | float:
    """Return P(X = k) = lam^k exp(-lam) / k!, evaluated through logarithms so large k and lam do not overflow."""
    if lam < 0.0:
        raise ValueError("lam must be non-negative")
    k_arr = np.atleast_1d(np.asarray(k, dtype=int))
    out = np.zeros(k_arr.shape, dtype=float)
    for i, kk in enumerate(k_arr.ravel()):
        if kk >= 0:
            if lam == 0.0:
                out.ravel()[i] = 1.0 if kk == 0 else 0.0
            else:
                out.ravel()[i] = math.exp(kk * math.log(lam) - lam - math.lgamma(kk + 1))
    return out if np.ndim(k) else float(out[0])


def binomial_poisson_gap(n: int, lam: float) -> float:
    """Return max_k |Binomial(n, lam/n)(k) - Poisson(lam)(k)|: how far the binomial is from its Poisson limit."""
    if n < 1 or lam <= 0.0 or lam > n:
        raise ValueError("need n >= 1 and 0 < lam <= n")
    k = np.arange(0, n + 1)
    return float(np.max(np.abs(binomial_pmf(k, n, lam / n) - poisson_pmf(k, lam))))


# ---------------------------------------------------------------- continuous random variables

def uniform_pdf(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Return the density 1/(b-a) on [a, b], zero outside."""
    if b <= a:
        raise ValueError("need b > a")
    x = np.asarray(x, dtype=float)
    return np.where((x >= a) & (x <= b), 1.0 / (b - a), 0.0)


def uniform_cdf(x: np.ndarray, a: float, b: float) -> np.ndarray:
    """Return P(X <= x) for the uniform distribution on [a, b]."""
    if b <= a:
        raise ValueError("need b > a")
    x = np.asarray(x, dtype=float)
    return np.clip((x - a) / (b - a), 0.0, 1.0)


def exponential_pdf(x: np.ndarray, lam: float) -> np.ndarray:
    """Return lam exp(-lam x) for x >= 0 (zero for x < 0); ``lam`` is the rate, 1/lam the mean."""
    if lam <= 0.0:
        raise ValueError("rate must be positive")
    x = np.asarray(x, dtype=float)
    return np.where(x >= 0.0, lam * np.exp(-lam * np.clip(x, 0.0, None)), 0.0)


def exponential_cdf(x: np.ndarray, lam: float) -> np.ndarray:
    """Return P(X <= x) = 1 - exp(-lam x) for x >= 0."""
    if lam <= 0.0:
        raise ValueError("rate must be positive")
    x = np.asarray(x, dtype=float)
    return np.where(x >= 0.0, 1.0 - np.exp(-lam * np.clip(x, 0.0, None)), 0.0)


def memoryless_gap(lam: float, s: float, t: float) -> float:
    """Return P(T > s + t | T > s) - P(T > t) for the exponential distribution; zero if it is memoryless."""
    surv = lambda u: 1.0 - float(exponential_cdf(np.array(u), lam))
    return surv(s + t) / surv(s) - surv(t)


def normal_pdf(x: np.ndarray, mu: float = 0.0, sigma: float = 1.0) -> np.ndarray:
    """Return the Gaussian density exp(-(x-mu)^2 / (2 sigma^2)) / (sigma sqrt(2 pi))."""
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    x = np.asarray(x, dtype=float)
    z = (x - mu) / sigma
    return np.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi))


def normal_cdf(x: np.ndarray, mu: float = 0.0, sigma: float = 1.0) -> np.ndarray:
    """Return P(X <= x) = (1 + erf(z / sqrt 2)) / 2 with z = (x - mu) / sigma.

    The error function is the one primitive borrowed from ``scipy.special``; the test checks that
    this equals the Simpson integral of ``normal_pdf``.
    """
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    x = np.asarray(x, dtype=float)
    return 0.5 * (1.0 + erf((x - mu) / (sigma * math.sqrt(2.0))))


def normal_quantile(p: float, tol: float = 1e-12) -> float:
    """Return z with Phi(z) = p for the standard normal, by bisection on ``normal_cdf``.

    Raise ``ValueError`` unless 0 < p < 1. Bisection on [-40, 40] converges to ``tol`` in about 50 steps.
    """
    if not 0.0 < p < 1.0:
        raise ValueError("p must lie strictly inside (0, 1)")
    lo, hi = -40.0, 40.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if float(normal_cdf(np.array(mid))) < p:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def z_score(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """Return (x - mu) / sigma: how many standard deviations x lies from the mean."""
    if sigma <= 0.0:
        raise ValueError("sigma must be positive")
    return (np.asarray(x, dtype=float) - mu) / sigma


def normal_interval(mu: float, sigma: float, prob: float) -> tuple[float, float]:
    """Return the symmetric interval (mu - z sigma, mu + z sigma) that holds probability ``prob``."""
    if not 0.0 < prob < 1.0:
        raise ValueError("prob must lie strictly inside (0, 1)")
    z = normal_quantile(0.5 + 0.5 * prob)
    return mu - z * sigma, mu + z * sigma


def continuous_mean(pdf: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int = 4000) -> float:
    """Return E[X] = integral of x f(x) over [a, b] by Simpson's rule (the pdf must be negligible outside)."""
    return simpson_integral(lambda x: x * pdf(x), a, b, n)


def continuous_variance(pdf: Callable[[np.ndarray], np.ndarray], a: float, b: float, n: int = 4000) -> float:
    """Return Var X = integral of (x - mu)^2 f(x) over [a, b] by Simpson's rule."""
    mu = continuous_mean(pdf, a, b, n)
    return simpson_integral(lambda x: (x - mu) ** 2 * pdf(x), a, b, n)


# ---------------------------------------------------------------- samples and the central limit theorem

def sample_mean_sd(x: np.ndarray) -> tuple[float, float]:
    """Return (mean, standard deviation with the n - 1 divisor) of a sample, by the two-pass formula."""
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        raise ValueError("need at least two samples")
    mu = float(x.mean())
    return mu, math.sqrt(float(np.sum((x - mu) ** 2)) / (x.size - 1))


def empirical_rule_fractions(x: np.ndarray, mu: float, sigma: float) -> np.ndarray:
    """Return the fractions of ``x`` within 1, 2 and 3 sigma of mu (a Gaussian gives 0.6827, 0.9545, 0.9973)."""
    z = np.abs(z_score(x, mu, sigma))
    return np.array([np.mean(z <= 1.0), np.mean(z <= 2.0), np.mean(z <= 3.0)])


def clt_sum_of_uniforms(n_terms: int, n_trials: int, seed: int = 0) -> np.ndarray:
    """Return standardised sums (S - n/2) / sqrt(n/12) of ``n_terms`` independent U(0, 1) draws, ``n_trials`` times.

    Each uniform has mean 1/2 and variance 1/12; the standardised sum tends to N(0, 1) as n grows.
    """
    if n_terms < 1 or n_trials < 1:
        raise ValueError("need n_terms >= 1 and n_trials >= 1")
    rng = np.random.default_rng(seed)
    s = rng.random((n_trials, n_terms)).sum(axis=1)
    return (s - 0.5 * n_terms) / math.sqrt(n_terms / 12.0)


def clt_distance_to_normal(n_terms: int, n_trials: int = 200_000, seed: int = 0) -> float:
    """Return the largest gap between the empirical cdf of ``clt_sum_of_uniforms`` and the standard normal cdf."""
    z = np.sort(clt_sum_of_uniforms(n_terms, n_trials, seed))
    ecdf = np.arange(1, z.size + 1) / z.size
    return float(np.max(np.abs(ecdf - normal_cdf(z))))


# ---------------------------------------------------------------- two-dimensional dispersion

def sample_covariance(points: np.ndarray) -> np.ndarray:
    """Return the sample covariance matrix (n - 1 divisor) of rows ``points`` (n, d), from the definition."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[0] < 2:
        raise ValueError("points must be (n >= 2, d)")
    d = points - points.mean(axis=0)
    return d.T @ d / (points.shape[0] - 1)


def propagate_covariance(J: np.ndarray, Sigma: np.ndarray) -> np.ndarray:
    """Return J Sigma J^T: the covariance of y = J x when x has covariance Sigma (linear propagation)."""
    J = np.asarray(J, dtype=float)
    Sigma = np.asarray(Sigma, dtype=float)
    if Sigma.shape != (J.shape[1], J.shape[1]):
        raise ValueError("Sigma must be square with the same size as J has columns")
    return J @ Sigma @ J.T


def landing_sensitivities() -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return (J, Sigma, names) of the SYNTHETIC landing-dispersion model.

    Landing offset (down-range x, cross-range y) [m] = J @ (input errors). Inputs and their 1-sigma
    values: wind-speed error 3.0 m/s; drag-coefficient fractional error 0.05; entry flight-path angle
    error 0.3 deg; heading error 0.5 deg. J holds the metres of landing offset per unit of each input.
    These are illustrative sensitivities for a parachute-recovered sounding-rocket payload, not data.
    """
    J = np.array([[150.0, -1200.0, -400.0, 0.0],
                  [40.0, 0.0, 0.0, 250.0]])
    sig = np.array([3.0, 0.05, 0.3, 0.5])
    names = ["wind speed [m/s]", "drag coefficient [-]", "flight-path angle [deg]", "heading [deg]"]
    return J, np.diag(sig ** 2), names


def simulate_landings(n_trials: int, J: np.ndarray, Sigma: np.ndarray, seed: int = 0) -> np.ndarray:
    """Return SYNTHETIC landing offsets (n_trials, 2) [m]: Gaussian input errors with covariance ``Sigma``, mapped by ``J``.

    Sigma may be non-diagonal; the inputs are drawn with ``rng.multivariate_normal``.
    """
    if n_trials < 1:
        raise ValueError("need n_trials >= 1")
    rng = np.random.default_rng(seed)
    J = np.asarray(J, dtype=float)
    x = rng.multivariate_normal(np.zeros(J.shape[1]), np.asarray(Sigma, dtype=float), size=n_trials)
    return x @ J.T


def eig2x2_symmetric(C: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (eigenvalues ascending, eigenvectors as columns) of a symmetric 2x2 matrix in closed form.

    With m = (a + d)/2 and r = sqrt(((a - d)/2)^2 + b^2): lambda = m -/+ r; the eigenvector of lambda_2
    is (b, lambda_2 - a) unless b = 0. Raise ``ValueError`` for a non-symmetric input.
    """
    C = np.asarray(C, dtype=float)
    if C.shape != (2, 2) or abs(C[0, 1] - C[1, 0]) > 1e-12 * max(1.0, np.abs(C).max()):
        raise ValueError("C must be symmetric 2x2")
    a, b, d = C[0, 0], C[0, 1], C[1, 1]
    m = 0.5 * (a + d)
    r = math.sqrt((0.5 * (a - d)) ** 2 + b * b)
    lam = np.array([m - r, m + r])
    if abs(b) < 1e-15 * max(1.0, abs(a), abs(d)):
        V = np.eye(2) if a <= d else np.array([[0.0, 1.0], [1.0, 0.0]])
    else:
        v2 = np.array([b, lam[1] - a])
        v2 /= np.linalg.norm(v2)
        v1 = np.array([-v2[1], v2[0]])
        V = np.column_stack([v1, v2])
    return lam, V


def ellipse_containment_probability(k: float) -> float:
    """Return P(inside the k-sigma ellipse) = 1 - exp(-k^2 / 2) for a two-dimensional Gaussian."""
    if k < 0.0:
        raise ValueError("k must be non-negative")
    return 1.0 - math.exp(-0.5 * k * k)


def k_for_probability(prob: float) -> float:
    """Return the ellipse scale k = sqrt(-2 ln(1 - prob)) that contains probability ``prob`` in two dimensions."""
    if not 0.0 <= prob < 1.0:
        raise ValueError("prob must lie in [0, 1)")
    return math.sqrt(-2.0 * math.log(1.0 - prob))


def dispersion_ellipse(cov: np.ndarray, prob: float) -> dict[str, float]:
    """Return the constant-probability ellipse of a 2-D Gaussian with covariance ``cov`` containing ``prob``.

    Keys: 'k' (scale), 'a' and 'b' (semi-major and semi-minor axes, metres, a >= b), 'angle_deg'
    (orientation of the major axis from +x towards +y, in (-90, 90]), 'area' (pi a b).
    """
    lam, V = eig2x2_symmetric(cov)
    if lam[0] <= 0.0:
        raise ValueError("covariance must be positive definite")
    k = k_for_probability(prob)
    a, b = k * math.sqrt(lam[1]), k * math.sqrt(lam[0])
    vx, vy = V[:, 1]
    angle = math.degrees(math.atan2(vy, vx))
    if angle <= -90.0:
        angle += 180.0
    elif angle > 90.0:
        angle -= 180.0
    return {"k": k, "a": a, "b": b, "angle_deg": angle, "area": math.pi * a * b}


def mahalanobis_squared(points: np.ndarray, mean: np.ndarray, cov: np.ndarray) -> np.ndarray:
    """Return (p - mean)^T cov^{-1} (p - mean) for every row of ``points``."""
    d = np.asarray(points, dtype=float) - np.asarray(mean, dtype=float)
    Ci = np.linalg.inv(np.asarray(cov, dtype=float))
    return np.einsum("ij,jk,ik->i", d, Ci, d)


def fraction_inside_ellipse(points: np.ndarray, mean: np.ndarray, cov: np.ndarray, k: float) -> float:
    """Return the fraction of ``points`` with Mahalanobis distance <= k (inside the k-sigma ellipse)."""
    return float(np.mean(mahalanobis_squared(points, mean, cov) <= k * k))


def circle_containment_probability(cov: np.ndarray, radius: float, n: int = 2000) -> float:
    """Return P(x^2 + y^2 <= radius^2) for a zero-mean 2-D Gaussian with covariance ``cov``, by quadrature.

    In principal axes the components are independent with variances lambda_1, lambda_2, so
    P = integral_{-R}^{R} f_1(u) [Phi(sqrt(R^2 - u^2)/s_2) - Phi(-sqrt(R^2 - u^2)/s_2)] du, done with
    Simpson's rule in ``n`` panels. For a circular covariance it must equal 1 - exp(-R^2 / (2 sigma^2)).
    """
    if radius < 0.0:
        raise ValueError("radius must be non-negative")
    if radius == 0.0:
        return 0.0
    lam, _ = eig2x2_symmetric(cov)
    s1, s2 = math.sqrt(lam[1]), math.sqrt(lam[0])

    def integrand(u: np.ndarray) -> np.ndarray:
        h = np.sqrt(np.clip(radius ** 2 - u ** 2, 0.0, None))
        return normal_pdf(u, 0.0, s1) * (normal_cdf(h, 0.0, s2) - normal_cdf(-h, 0.0, s2))

    return simpson_integral(integrand, -radius, radius, n)


def radius_for_probability(cov: np.ndarray, prob: float, tol: float = 1e-6) -> float:
    """Return the radius of the circle (centred on the mean) that contains probability ``prob``, by bisection.

    With prob = 0.5 this is the circular error probable (CEP). For a circular covariance sigma^2 I it
    equals sigma sqrt(-2 ln(1 - prob)), i.e. 1.1774 sigma for the CEP.
    """
    if not 0.0 < prob < 1.0:
        raise ValueError("prob must lie strictly inside (0, 1)")
    lam, _ = eig2x2_symmetric(cov)
    lo, hi = 0.0, 10.0 * math.sqrt(lam[1])
    while hi - lo > tol * hi:
        mid = 0.5 * (lo + hi)
        if circle_containment_probability(cov, mid) < prob:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the lesson's tables and write landing_dispersion.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=4, suppress=True)

    print("== discrete random variables ==")
    k = np.arange(0, 4)
    pmf = binomial_pmf(k, 3, 0.98)
    print(f"working IMUs out of 3 (p = 0.98): pmf {pmf}, mean {discrete_mean(k, pmf):.4f} (np = 2.94), "
          f"variance {discrete_variance(k, pmf):.4f} (np(1-p) = {3 * 0.98 * 0.02:.4f})")
    lam = 2.4
    k = np.arange(0, 15)
    pp = poisson_pmf(k, lam)
    print(f"debris conjunction alerts per month, Poisson(2.4): P(0) = {pp[0]:.4f}, P(>= 5) = {1 - pp[:5].sum():.4f}, "
          f"mean {discrete_mean(k, pp / pp.sum()):.4f}, variance {discrete_variance(k, pp / pp.sum()):.4f}")
    for n in (10, 100, 1000):
        print(f"  binomial(n = {n}, p = 2.4/n) vs Poisson(2.4): max pmf gap {binomial_poisson_gap(n, lam):.5f}")

    print("\n== continuous random variables ==")
    mu_u = continuous_mean(lambda x: uniform_pdf(x, -0.5, 0.5), -0.5, 0.5)
    var_u = continuous_variance(lambda x: uniform_pdf(x, -0.5, 0.5), -0.5, 0.5)
    print(f"quantisation error U(-0.5, 0.5) LSB: mean {mu_u:.2e}, variance {var_u:.6f} (1/12 = {1 / 12:.6f})")
    rate = 1.0 / 5000.0
    print(f"exponential lifetime, MTBF 5000 h: P(fails before 1000 h) = {float(exponential_cdf(np.array(1000.0), rate)):.4f}, "
          f"P(survives 10000 h) = {1 - float(exponential_cdf(np.array(10000.0), rate)):.4f}, "
          f"memoryless gap {memoryless_gap(rate, 3000.0, 1000.0):.1e}, "
          f"mean by Simpson {continuous_mean(lambda x: exponential_pdf(x, rate), 0.0, 2e5):.1f} h")
    print(f"normal: P(|z| <= 1, 2, 3) = {[round(2 * float(normal_cdf(np.array(z))) - 1, 4) for z in (1, 2, 3)]}; "
          f"z for 95 % = {normal_quantile(0.975):.4f}, for 99 % = {normal_quantile(0.995):.4f}")
    lo, hi = normal_interval(0.0, 0.02, 0.95)
    print(f"accelerometer noise sigma = 0.02 m/s^2: 95 % of samples within ({lo:.4f}, {hi:.4f}) m/s^2")

    print("\n== central limit theorem (SYNTHETIC uniform noise) ==")
    for n in (1, 2, 5, 12, 30):
        z = clt_sum_of_uniforms(n, 200_000, seed=16)
        fr = empirical_rule_fractions(z, 0.0, 1.0)
        print(f"  n = {n:2d}: sd of standardised sum {z.std(ddof=1):.4f}, within 1/2/3 sigma {fr}, "
              f"max cdf gap to N(0,1) {clt_distance_to_normal(n, 200_000, seed=16):.4f}")

    print("\n== landing dispersion (SYNTHETIC sensitivity model) ==")
    J, Sigma, names = landing_sensitivities()
    cov = propagate_covariance(J, Sigma)
    sx, sy = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
    rho = cov[0, 1] / (sx * sy)
    print("input 1-sigma:", dict(zip(names, np.sqrt(np.diag(Sigma)))))
    for i, name in enumerate(names):
        print(f"  {name:26s} contributes {J[0, i] * math.sqrt(Sigma[i, i]):8.1f} m down-range, "
              f"{J[1, i] * math.sqrt(Sigma[i, i]):7.1f} m cross-range (1 sigma)")
    print(f"covariance [m^2]:\n{cov}")
    print(f"sigma_x = {sx:.1f} m, sigma_y = {sy:.1f} m, correlation {rho:.3f}")
    landings = simulate_landings(20_000, J, Sigma, seed=2016)
    cov_hat = sample_covariance(landings)
    print(f"sample covariance from 20000 synthetic landings:\n{cov_hat}")
    print(f"  relative error {np.abs(cov_hat - cov).max() / np.abs(cov).max():.4f}")
    for prob in (0.3935, 0.8647, 0.95, 0.99):
        e = dispersion_ellipse(cov, prob)
        inside = fraction_inside_ellipse(landings, np.zeros(2), cov, e["k"])
        print(f"  {prob:.4f} ellipse: k = {e['k']:.3f}, a = {e['a']:.0f} m, b = {e['b']:.0f} m, "
              f"angle {e['angle_deg']:.1f} deg, area {e['area'] / 1e6:.3f} km^2; simulated fraction inside {inside:.4f}")
    cep = radius_for_probability(cov, 0.5)
    r99 = radius_for_probability(cov, 0.99)
    e99 = dispersion_ellipse(cov, 0.99)
    print(f"CEP (50 % circle) = {cep:.0f} m; 99 % circle radius = {r99:.0f} m, area {math.pi * r99 ** 2 / 1e6:.3f} km^2 "
          f"vs 99 % ellipse area {e99['area'] / 1e6:.3f} km^2 (ratio {math.pi * r99 ** 2 / e99['area']:.2f})")
    r = np.sqrt((landings ** 2).sum(axis=1))
    print(f"simulated fraction within CEP {np.mean(r <= cep):.4f}, within 99 % circle {np.mean(r <= r99):.4f}")
    circ = np.array([[sx ** 2, 0.0], [0.0, sx ** 2]])
    print(f"check, circular sigma = {sx:.1f} m: CEP by quadrature {radius_for_probability(circ, 0.5):.2f} m, "
          f"closed form 1.1774 sigma = {sx * math.sqrt(2 * math.log(2)):.2f} m")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    x = np.linspace(-4, 4, 400)
    for n in (1, 2, 12):
        z = clt_sum_of_uniforms(n, 200_000, seed=16)
        axes[0].hist(z, bins=80, range=(-4, 4), density=True, histtype="step", label=f"n = {n}")
    axes[0].plot(x, normal_pdf(x), "k--", label="N(0, 1)")
    axes[0].set_xlabel("standardised sum of n uniforms")
    axes[0].set_ylabel("density")
    axes[0].legend()
    axes[0].set_title("central limit theorem (synthetic)")
    xs = np.linspace(0, 20000, 400)
    axes[1].plot(xs, exponential_pdf(xs, rate) * 5000, label="exponential pdf x MTBF")
    axes[1].plot(xs, exponential_cdf(xs, rate), label="exponential cdf")
    axes[1].plot(xs, normal_cdf(xs, 5000, 1500), label="normal cdf (mean 5000, sd 1500)")
    axes[1].set_xlabel("time [h]")
    axes[1].legend(fontsize=8)
    axes[1].set_title("pdf and cdf")
    axes[2].plot(landings[:3000, 0], landings[:3000, 1], ".", ms=2, alpha=0.4, label="3000 synthetic landings")
    th = np.linspace(0, 2 * np.pi, 200)
    lam, V = eig2x2_symmetric(cov)
    for prob, col in ((0.3935, "C1"), (0.8647, "C2"), (0.99, "C3")):
        e = dispersion_ellipse(cov, prob)
        pts = V @ np.vstack([e["b"] * np.cos(th), e["a"] * np.sin(th)])
        axes[2].plot(pts[0], pts[1], col, label=f"{prob:.0%} ellipse (k = {e['k']:.2f})")
    axes[2].plot(r99 * np.cos(th), r99 * np.sin(th), "k:", label="99 % circle")
    axes[2].set_aspect("equal")
    axes[2].set_xlabel("down-range offset [m]")
    axes[2].set_ylabel("cross-range offset [m]")
    axes[2].legend(fontsize=7, loc="upper left")
    axes[2].set_title("landing dispersion (synthetic)")
    fig.tight_layout()
    fig.savefig(out_dir / "landing_dispersion.png", dpi=150)
    print(f"saved {out_dir / 'landing_dispersion.png'}")


if __name__ == "__main__":
    main()
