"""P02.20 Bayesian Inference: reference solution.

Update a failure-rate estimate for a spacecraft reaction wheel from heritage data and new fleet
hours: Bayes' theorem with discrete and continuous priors, the Gamma-Poisson and Beta-binomial
conjugate updates derived in the lesson, a grid posterior and a Metropolis sampler written from
scratch and checked against ``scipy.stats`` as the oracle, credible intervals (equal-tailed and
highest-density), posterior predictive reliability, the frequentist chi-square bound for
comparison, test planning (hours or firings needed to demonstrate a requirement) and a
calibration experiment on synthetic fleets.

Conventions:
    * A constant failure rate lambda [1/h]: failures in exposure T hours are Poisson with mean lambda T.
    * Gamma(alpha, beta) is the *rate* parametrisation: pdf beta^alpha lambda^(alpha-1) e^(-beta lambda) / Gamma(alpha),
      mean alpha/beta, variance alpha/beta^2. alpha counts "prior failures", beta counts "prior hours".
    * Beta(a, b) for a per-demand failure probability p: mean a/(a+b).
    * Every number used as data in this lesson is a lesson parameter or a seeded synthetic draw,
      not a measurement.

Run:  python solution.py          (prints the update, writes failure_rate_posterior.png)
Test: AEFS_IMPL=solution pytest curriculum/phase-02/lesson-20-bayesian-inference
"""
from __future__ import annotations

from collections.abc import Callable
from math import lgamma
from pathlib import Path

import numpy as np
from scipy.special import gammainc, gammaincinv
from scipy.stats import chi2

HOURS_PER_YEAR = 8766.0

# Reaction-wheel reliability problem (lesson parameters, not measurements).
HERITAGE = {"failures": 2, "hours": 1.2e6}          # predecessor wheel: 2 failures in 1.2 million unit-hours
FLEET = {"failures": 1, "hours": 6.0e4}             # new wheel: 10 units x 6000 h, one failure
MISSION_HOURS = 3.0 * HOURS_PER_YEAR                # 3-year mission: 26 298 h
RELIABILITY_REQUIRED = 0.95                         # P(no failure in the mission) >= 0.95
PYRO = {"failures": 0, "firings": 40}               # separation nuts: 40 test firings, no failure


# ---------------------------------------------------------------- Bayes' theorem, discrete

def bayes_discrete(prior: np.ndarray, likelihood: np.ndarray) -> tuple[np.ndarray, float]:
    """Return (posterior, evidence) for hypotheses with the given prior probabilities and likelihoods of the data.

    posterior_i = prior_i * likelihood_i / evidence, evidence = sum_i prior_i * likelihood_i.
    Raise ``ValueError`` if the prior does not sum to 1, has negative entries, or the evidence is zero.
    """
    prior = np.asarray(prior, dtype=float)
    like = np.asarray(likelihood, dtype=float)
    if prior.shape != like.shape or np.any(prior < 0.0) or np.any(like < 0.0):
        raise ValueError("prior and likelihood must be non-negative arrays of the same shape")
    if not np.isclose(prior.sum(), 1.0, rtol=0.0, atol=1e-9):
        raise ValueError("prior must sum to 1")
    joint = prior * like
    evidence = float(joint.sum())
    if evidence <= 0.0:
        raise ValueError("the data are impossible under every hypothesis")
    return joint / evidence, evidence


def odds_update(prior_odds: float, bayes_factor: float) -> float:
    """Return the posterior odds = prior odds x Bayes factor (likelihood ratio of the data under the two hypotheses)."""
    if prior_odds < 0.0 or bayes_factor < 0.0:
        raise ValueError("odds and Bayes factors are non-negative")
    return prior_odds * bayes_factor


# ---------------------------------------------------------------- likelihoods and densities (log-space, from scratch)

def log_binomial_likelihood(k: int, n: int, p: np.ndarray) -> np.ndarray:
    """Return log P(k failures in n demands | p) = log C(n, k) + k log p + (n - k) log(1 - p), elementwise in p."""
    if not 0 <= k <= n:
        raise ValueError("need 0 <= k <= n")
    p = np.asarray(p, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        t_fail = np.where(k == 0, 0.0, k * np.log(p))                 # 0 log 0 = 0 at p = 0
        t_ok = np.where(n - k == 0, 0.0, (n - k) * np.log1p(-p))      # and at p = 1
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1) + t_fail + t_ok


def log_poisson_likelihood(k: int, T: float, lam: np.ndarray) -> np.ndarray:
    """Return log P(k failures in T hours | lambda) = k log(lambda T) - lambda T - log k!, elementwise in lambda."""
    if k < 0 or T <= 0.0:
        raise ValueError("need k >= 0 and T > 0")
    lam = np.asarray(lam, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = np.where(lam > 0.0, k * np.log(lam * T), 0.0 if k == 0 else -np.inf)
    return ll - lam * T - lgamma(k + 1)


def log_gamma_pdf(lam: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """Return log of the Gamma(alpha, beta) density in the rate parametrisation, elementwise; -inf for lambda <= 0.

    For beta = 0 (an improper prior such as Jeffreys' alpha = 1/2 or the flat alpha = 1) the
    normalising constant is dropped and (alpha - 1) log lambda is returned.
    """
    if alpha <= 0.0 or beta < 0.0:
        raise ValueError("alpha must be positive and beta non-negative")
    lam = np.asarray(lam, dtype=float)
    with np.errstate(divide="ignore"):
        ll = np.where(lam > 0.0, (alpha - 1.0) * np.log(np.where(lam > 0.0, lam, 1.0)) - beta * lam, -np.inf)
    if alpha == 1.0:
        ll = np.where(lam >= 0.0, -beta * lam, -np.inf)                       # lambda^0 = 1 at the origin
    if beta > 0.0:
        ll = ll + alpha * np.log(beta) - lgamma(alpha)
    return ll


def log_beta_pdf(p: np.ndarray, a: float, b: float) -> np.ndarray:
    """Return log of the Beta(a, b) density, elementwise; -inf outside [0, 1]."""
    if a <= 0.0 or b <= 0.0:
        raise ValueError("a and b must be positive")
    p = np.asarray(p, dtype=float)
    inside = (p >= 0.0) & (p <= 1.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ll = (a - 1.0) * np.log(p) + (b - 1.0) * np.log1p(-p)
    ll = ll + lgamma(a + b) - lgamma(a) - lgamma(b)
    return np.where(inside, ll, -np.inf)


# ---------------------------------------------------------------- conjugate updates

def gamma_poisson_update(alpha: float, beta: float, k: int, T: float) -> tuple[float, float]:
    """Return the posterior (alpha + k, beta + T) of a Gamma(alpha, beta) prior after k failures in T hours."""
    if alpha <= 0.0 or beta < 0.0 or k < 0 or T < 0.0:
        raise ValueError("bad update parameters")
    return alpha + k, beta + T


def beta_binomial_update(a: float, b: float, k: int, n: int) -> tuple[float, float]:
    """Return the posterior (a + k, b + n - k) of a Beta(a, b) prior after k failures in n demands."""
    if a <= 0.0 or b <= 0.0 or not 0 <= k <= n:
        raise ValueError("bad update parameters")
    return a + k, b + n - k


def sequential_update(alpha: float, beta: float, events: list[tuple[int, float]]) -> tuple[float, float]:
    """Apply ``gamma_poisson_update`` once per (k, T) record in ``events``, in order, and return the final (alpha, beta)."""
    for k, T in events:
        alpha, beta = gamma_poisson_update(alpha, beta, k, T)
    return alpha, beta


# ---------------------------------------------------------------- grid posterior (from scratch)

def grid_posterior(grid: np.ndarray, log_prior: np.ndarray, log_likelihood: np.ndarray) -> np.ndarray:
    """Return the posterior density on ``grid`` from log prior and log likelihood values, normalised by the trapezoid rule.

    Works in log space (subtracting the maximum before exponentiating) so that likelihoods like
    lambda^k e^(-lambda T) with T ~ 10^6 h do not underflow. Raise ``ValueError`` if nothing is finite.
    """
    grid = np.asarray(grid, dtype=float)
    lp = np.asarray(log_prior, dtype=float) + np.asarray(log_likelihood, dtype=float)
    if grid.ndim != 1 or lp.shape != grid.shape or np.any(np.diff(grid) <= 0.0):
        raise ValueError("grid must be a strictly increasing 1-D array matching the log arrays")
    finite = np.isfinite(lp)
    if not np.any(finite):
        raise ValueError("the posterior is zero everywhere on the grid")
    dens = np.zeros_like(grid)
    dens[finite] = np.exp(lp[finite] - lp[finite].max())
    area = np.trapezoid(dens, grid)
    if area <= 0.0:
        raise ValueError("the posterior has zero area on the grid")
    return dens / area


def grid_mean(grid: np.ndarray, pdf: np.ndarray) -> float:
    """Return the mean of a density tabulated on ``grid`` (trapezoid rule)."""
    return float(np.trapezoid(grid * pdf, grid))


def grid_cdf(grid: np.ndarray, pdf: np.ndarray) -> np.ndarray:
    """Return the cumulative distribution on ``grid`` (cumulative trapezoid), from 0 to 1."""
    grid = np.asarray(grid, dtype=float)
    pdf = np.asarray(pdf, dtype=float)
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pdf[1:] + pdf[:-1]) * np.diff(grid))])
    return cdf / cdf[-1]


def grid_quantile(grid: np.ndarray, pdf: np.ndarray, q: float) -> float:
    """Return the value x with CDF(x) = q by linear interpolation of the tabulated CDF."""
    if not 0.0 <= q <= 1.0:
        raise ValueError("q must lie in [0, 1]")
    return float(np.interp(q, grid_cdf(grid, pdf), grid))


def credible_interval(grid: np.ndarray, pdf: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    """Return the equal-tailed credible interval (lower, upper) that holds probability ``level``."""
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie in (0, 1)")
    tail = 0.5 * (1.0 - level)
    return grid_quantile(grid, pdf, tail), grid_quantile(grid, pdf, 1.0 - tail)


def hpd_interval(grid: np.ndarray, pdf: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    """Return the highest-posterior-density interval: the shortest interval holding probability ``level``.

    Found by lowering a density threshold until the mass above it reaches ``level``; assumes a
    unimodal density on the grid. Raise ``ValueError`` for a bad level.
    """
    if not 0.0 < level < 1.0:
        raise ValueError("level must lie in (0, 1)")
    grid = np.asarray(grid, dtype=float)
    pdf = np.asarray(pdf, dtype=float)
    cdf = grid_cdf(grid, pdf)
    for thr in np.sort(pdf)[::-1]:
        inside = pdf >= thr
        lo, hi = np.flatnonzero(inside)[[0, -1]]
        if cdf[hi] - cdf[lo] >= level:
            return float(grid[lo]), float(grid[hi])
    return float(grid[0]), float(grid[-1])


# ---------------------------------------------------------------- Metropolis sampler (from scratch)

def metropolis(log_target: Callable[[float], float], x0: float, n_samples: int, step: float,
               rng: np.random.Generator, burn_in: int = 1000) -> tuple[np.ndarray, float]:
    """Return (samples, acceptance rate) from a random-walk Metropolis chain targeting exp(log_target).

    Proposal x' = x + step * N(0, 1); accept with probability min(1, exp(log_target(x') - log_target(x))).
    ``burn_in`` initial draws are discarded. A proposal with log_target = -inf is always rejected,
    which is how a support boundary (lambda > 0) is enforced.
    """
    if n_samples <= 0 or step <= 0.0 or burn_in < 0:
        raise ValueError("n_samples and step must be positive, burn_in non-negative")
    x = float(x0)
    lx = float(log_target(x))
    if not np.isfinite(lx):
        raise ValueError("x0 must have finite log density")
    out = np.empty(n_samples)
    accepted = 0
    total = n_samples + burn_in
    innov = rng.normal(0.0, step, size=total)
    uni = np.log(rng.uniform(size=total))
    for i in range(total):
        xp = x + innov[i]
        lxp = float(log_target(xp))
        if np.isfinite(lxp) and uni[i] < lxp - lx:
            x, lx = xp, lxp
            if i >= burn_in:
                accepted += 1
        if i >= burn_in:
            out[i - burn_in] = x
    return out, accepted / n_samples


def autocorrelation_time(samples: np.ndarray, max_lag: int = 200) -> float:
    """Return the integrated autocorrelation time 1 + 2 sum_k rho_k (summing until rho first drops below 0).

    Effective sample size = n / autocorrelation_time.
    """
    x = np.asarray(samples, dtype=float)
    x = x - x.mean()
    var = float(np.dot(x, x) / x.size)
    if var == 0.0:
        return 1.0
    tau = 1.0
    for k in range(1, min(max_lag, x.size - 1) + 1):
        rho = float(np.dot(x[:-k], x[k:]) / x.size) / var
        if rho <= 0.0:
            break
        tau += 2.0 * rho
    return tau


# ---------------------------------------------------------------- what the posterior means for the mission

def prob_rate_below(alpha: float, beta: float, lam: float) -> float:
    """Return P(lambda <= lam) under Gamma(alpha, beta): the regularised lower incomplete gamma P(alpha, beta lam)."""
    if alpha <= 0.0 or beta <= 0.0 or lam < 0.0:
        raise ValueError("alpha, beta must be positive and lam non-negative")
    return float(gammainc(alpha, beta * lam))


def posterior_predictive_reliability(alpha: float, beta: float, t: float) -> float:
    """Return P(no failure in the next t hours) averaged over Gamma(alpha, beta): (beta / (beta + t))^alpha.

    E[e^(-lambda t)] is the Gamma moment-generating function at -t. It is always *above*
    e^(-mean t) because e^(-lambda t) is convex (Jensen's inequality): the chance that the rate is
    low helps more than the chance that it is high hurts, for the probability of zero failures.
    """
    if alpha <= 0.0 or beta <= 0.0 or t < 0.0:
        raise ValueError("alpha, beta must be positive and t non-negative")
    return float((beta / (beta + t)) ** alpha)


def hours_to_demonstrate(alpha: float, beta: float, lam_req: float, conf: float) -> float:
    """Return the additional failure-free hours T such that P(lambda <= lam_req | data + T) = conf.

    P(alpha, (beta + T) lam_req) = conf gives (beta + T) lam_req = P^{-1}(alpha, conf); T = that / lam_req - beta,
    clipped at 0 when the requirement is already demonstrated.
    """
    if not 0.0 < conf < 1.0 or lam_req <= 0.0:
        raise ValueError("conf must lie in (0, 1) and lam_req be positive")
    x = float(gammaincinv(alpha, conf))
    return max(0.0, x / lam_req - beta)


def frequentist_upper_bound(k: int, T: float, conf: float) -> float:
    """Return the one-sided upper confidence bound on lambda: chi2_{2k+2}(conf) / (2T).

    It equals the Bayesian upper credible bound from the flat prior Gamma(1, 0), whose posterior is
    Gamma(k + 1, T) = chi2 with 2k + 2 degrees of freedom scaled by 1/(2T).
    """
    if k < 0 or T <= 0.0 or not 0.0 < conf < 1.0:
        raise ValueError("bad inputs")
    return float(chi2.ppf(conf, 2 * k + 2) / (2.0 * T))


def zero_failure_tests_required(p_max: float, conf: float) -> int:
    """Return the smallest n with P(p <= p_max | 0 failures in n, flat prior) >= conf, i.e. 1 - (1 - p_max)^(n+1) >= conf.

    With the flat Beta(1, 1) prior the posterior after n successes is Beta(1, n + 1), whose CDF at
    p_max is 1 - (1 - p_max)^(n + 1).
    """
    if not 0.0 < p_max < 1.0 or not 0.0 < conf < 1.0:
        raise ValueError("p_max and conf must lie in (0, 1)")
    n = np.log(1.0 - conf) / np.log(1.0 - p_max) - 1.0
    return int(np.ceil(n - 1e-12))


def beta_cdf_grid(a: float, b: float, p_max: float, n_grid: int = 20001) -> float:
    """Return P(p <= p_max) under Beta(a, b) by integrating the tabulated density (a from-scratch incomplete beta).

    The density can be singular at p = 0 or p = 1 (a < 1 or b < 1, as for Jeffreys' Beta(1/2, 1/2)), so the
    grid is placed in u with p = sin^2(pi u / 2): the Jacobian dp/du = (pi/2) sin(pi u) vanishes at both
    ends and tames any singularity weaker than 1/p, after which the trapezoid rule converges.
    """
    if not 0.0 <= p_max <= 1.0:
        raise ValueError("p_max must lie in [0, 1]")
    u = np.linspace(0.0, 1.0, n_grid)
    p = np.sin(0.5 * np.pi * u) ** 2
    with np.errstate(divide="ignore", invalid="ignore"):
        log_dens = log_beta_pdf(p, a, b) + np.log(0.5 * np.pi * np.sin(np.pi * u))
    log_dens[~np.isfinite(log_dens)] = -np.inf          # the two end points (0 * inf): a single zero costs O(1/n_grid)
    dens_u = grid_posterior(u, log_dens, np.zeros_like(u))
    u_max = 2.0 / np.pi * np.arcsin(np.sqrt(p_max))
    return float(np.interp(u_max, u, grid_cdf(u, dens_u)))


# ---------------------------------------------------------------- synthetic fleets and calibration

def simulate_fleet_failures(lam: float, T: float, rng: np.random.Generator) -> int:
    """Return a synthetic failure count: one Poisson draw with mean lambda T."""
    if lam < 0.0 or T < 0.0:
        raise ValueError("lam and T must be non-negative")
    return int(rng.poisson(lam * T))


def credible_interval_coverage(alpha0: float, beta0: float, T: float, n_trials: int, level: float,
                               seed: int) -> float:
    """Return the fraction of synthetic fleets whose ``level`` equal-tailed posterior interval contains the true rate.

    Each trial draws lambda from the prior Gamma(alpha0, beta0), a Poisson count in T hours, updates
    to Gamma(alpha0 + k, beta0 + T) and checks whether lambda lies inside the interval (computed with
    ``scipy``'s gamma quantiles). Averaged over draws from the prior, Bayesian intervals are calibrated:
    the fraction tends to ``level``.
    """
    from scipy.stats import gamma as gamma_dist

    if n_trials <= 0 or not 0.0 < level < 1.0:
        raise ValueError("bad inputs")
    rng = np.random.default_rng(seed)
    tail = 0.5 * (1.0 - level)
    hits = 0
    for _ in range(n_trials):
        lam = rng.gamma(alpha0, 1.0 / beta0)
        k = simulate_fleet_failures(lam, T, rng)
        a, b = gamma_poisson_update(alpha0, beta0, k, T)
        lo, hi = gamma_dist.ppf([tail, 1.0 - tail], a, scale=1.0 / b)
        hits += int(lo <= lam <= hi)
    return hits / n_trials


def main(out_dir: str | Path = ".") -> None:
    """Run the reaction-wheel update, the pyro demonstration, the prior-sensitivity table and the calibration experiment."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.stats import gamma as gamma_dist

    out_dir = Path(out_dir)
    M = 1e6                                               # report rates per million hours
    lam_req = -np.log(RELIABILITY_REQUIRED) / MISSION_HOURS
    print(f"mission {MISSION_HOURS:.0f} h, reliability {RELIABILITY_REQUIRED}: required lambda <= {lam_req * M:.3f} per 10^6 h")

    # 1. Discrete update: is the lot good or bad?
    p_hyp = np.array([0.01, 0.10])
    prior = np.array([0.9, 0.1])
    like = np.exp(log_binomial_likelihood(1, 20, p_hyp))
    post, ev = bayes_discrete(prior, like)
    print(f"\nlot test, 1 failure in 20: likelihoods {like.round(4)}, evidence {ev:.4f}, posterior P(good) = {post[0]:.4f}")
    print(f"   prior odds good:bad 9:1, Bayes factor {like[0] / like[1]:.3f}, posterior odds {odds_update(9.0, like[0] / like[1]):.3f}")

    # 2. Reaction wheel: three priors, one data set.
    kh, Th = HERITAGE["failures"], HERITAGE["hours"]
    kf, Tf = FLEET["failures"], FLEET["hours"]
    priors = {
        "Jeffreys (alpha=1/2, beta=0)": (0.5, 0.0),
        "flat (alpha=1, beta=0)": (1.0, 0.0),
        "heritage (Jeffreys + 2 in 1.2e6 h)": gamma_poisson_update(0.5, 0.0, kh, Th),
    }
    print(f"\nnew wheel: {kf} failure in {Tf:.0f} h; maximum-likelihood rate {kf / Tf * M:.2f} per 10^6 h")
    print(f"{'prior':38s} {'post alpha':>10s} {'post beta/1e6':>13s} {'mean':>7s} {'90% ETI':>16s} {'P(req)':>7s} {'R(3 yr)':>8s} {'hours to 90%':>13s}")
    results = {}
    for name, (a0, b0) in priors.items():
        a, b = gamma_poisson_update(a0, b0, kf, Tf)
        lo, hi = gamma_dist.ppf([0.05, 0.95], a, scale=1.0 / b)
        results[name] = (a, b)
        print(f"{name:38s} {a:10.1f} {b / M:13.3f} {a / b * M:7.3f} [{lo * M:6.3f}, {hi * M:6.3f}] "
              f"{prob_rate_below(a, b, lam_req):7.3f} {posterior_predictive_reliability(a, b, MISSION_HOURS):8.4f} "
              f"{hours_to_demonstrate(a, b, lam_req, 0.9) / M:13.3f}e6")
    a, b = results["heritage (Jeffreys + 2 in 1.2e6 h)"]
    print(f"heritage posterior: mean {a / b * M:.3f}, SD {np.sqrt(a) / b * M:.3f} per 10^6 h;"
          f" e^(-mean t) = {np.exp(-a / b * MISSION_HOURS):.4f} vs predictive {posterior_predictive_reliability(a, b, MISSION_HOURS):.4f}")
    print(f"frequentist 95% upper bound from the fleet data alone: {frequentist_upper_bound(kf, Tf, 0.95) * M:.3f} per 10^6 h;"
          f" flat-prior 95% credible upper bound: {gamma_dist.ppf(0.95, kf + 1, scale=1.0 / Tf) * M:.3f}")

    # 3. Grid posterior and Metropolis against the conjugate closed form (heritage prior).
    grid = np.linspace(0.0, 15.0 / M, 6001)
    a0, b0 = priors["heritage (Jeffreys + 2 in 1.2e6 h)"]
    pdf = grid_posterior(grid, log_gamma_pdf(grid, a0, b0), log_poisson_likelihood(kf, Tf, grid))
    eti = credible_interval(grid, pdf, 0.9)
    hpd = hpd_interval(grid, pdf, 0.9)
    print(f"\ngrid posterior: mean {grid_mean(grid, pdf) * M:.4f} (closed form {a / b * M:.4f}),"
          f" 90% ETI [{eti[0] * M:.3f}, {eti[1] * M:.3f}], 90% HPD [{hpd[0] * M:.3f}, {hpd[1] * M:.3f}]"
          f" (width {(hpd[1] - hpd[0]) * M:.3f} vs {(eti[1] - eti[0]) * M:.3f})")
    log_target = lambda x: float(log_gamma_pdf(x, a0, b0) + log_poisson_likelihood(kf, Tf, x)) if x > 0.0 else -np.inf
    samples, acc = metropolis(log_target, a / b, 40_000, 1.5 * np.sqrt(a) / b, np.random.default_rng(20))
    tau_int = autocorrelation_time(samples)
    print(f"Metropolis: 40000 draws, acceptance {acc:.2f}, autocorrelation time {tau_int:.1f}, effective samples {40000 / tau_int:.0f};"
          f" mean {samples.mean() * M:.4f}, SD {samples.std() * M:.4f} (closed form {np.sqrt(a) / b * M:.4f}),"
          f" P(req) {np.mean(samples <= lam_req):.3f} (closed form {prob_rate_below(a, b, lam_req):.3f})")

    # 4. Pyro separation nuts: zero-failure demonstration.
    n, kp = PYRO["firings"], PYRO["failures"]
    ap, bp = beta_binomial_update(1.0, 1.0, kp, n)
    print(f"\npyro nuts: {kp} failures in {n} firings. Flat prior -> Beta({ap:.0f}, {bp:.0f}); P(p <= 0.05) = {beta_cdf_grid(ap, bp, 0.05):.4f}"
          f" (closed form {1 - 0.95 ** bp:.4f}); firings needed for 95%: {zero_failure_tests_required(0.05, 0.95)}")
    aj, bj = beta_binomial_update(0.5, 0.5, kp, n)
    print(f"   Jeffreys prior -> Beta({aj}, {bj}); P(p <= 0.05) = {beta_cdf_grid(aj, bj, 0.05):.4f}")

    # 5. Calibration experiment on synthetic fleets.
    cov = credible_interval_coverage(a0, b0, Tf, 4000, 0.9, seed=5)
    print(f"\ncalibration: 90% intervals from the heritage prior cover the true rate in {cov:.3f} of 4000 synthetic fleets")

    # 6. Sequential updating is order-independent.
    events = [(0, 2.0e4), (1, 1.5e4), (0, 2.5e4)]
    print("sequential update:", sequential_update(a0, b0, events), "reversed:", sequential_update(a0, b0, events[::-1]))

    # Figure.
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    x = grid * M
    axs[0, 0].plot(x, np.exp(log_gamma_pdf(grid, a0, b0)) / M, "k--", label="prior (heritage)")
    like_curve = np.exp(log_poisson_likelihood(kf, Tf, grid))
    axs[0, 0].plot(x, like_curve / np.trapezoid(like_curve, x), ":", label="likelihood (scaled)")
    axs[0, 0].plot(x, pdf / M, label="posterior (grid)")
    axs[0, 0].hist(samples * M, bins=80, density=True, alpha=0.3, label="posterior (Metropolis)")
    axs[0, 0].axvline(lam_req * M, color="red", lw=0.8, label="requirement")
    axs[0, 0].set_xlabel("failure rate [per 10^6 h]")
    axs[0, 0].set_ylabel("density")
    axs[0, 0].set_title("reaction-wheel rate: prior x likelihood = posterior")
    axs[0, 0].legend(fontsize=8)
    for name, (aa, bb) in results.items():
        axs[0, 1].plot(x, gammainc(aa, bb * grid), label=name.split(" (")[0])
    axs[0, 1].axvline(lam_req * M, color="red", lw=0.8)
    axs[0, 1].axhline(0.9, color="gray", lw=0.5)
    axs[0, 1].set_xlabel("failure rate [per 10^6 h]")
    axs[0, 1].set_ylabel("posterior CDF")
    axs[0, 1].set_title("prior sensitivity")
    axs[0, 1].legend(fontsize=8)
    ns = np.arange(0, 121)
    axs[1, 0].plot(ns, 1.0 - 0.95 ** (ns + 1), label="flat prior")
    axs[1, 0].plot(ns, [beta_cdf_grid(0.5, 0.5 + nn, 0.05, 4001) for nn in ns], "--", label="Jeffreys prior")
    axs[1, 0].axhline(0.95, color="gray", lw=0.5)
    axs[1, 0].axvline(zero_failure_tests_required(0.05, 0.95), color="red", lw=0.8)
    axs[1, 0].set_xlabel("failure-free firings n")
    axs[1, 0].set_ylabel("P(p <= 0.05 | 0 failures in n)")
    axs[1, 0].set_title("zero-failure demonstration")
    axs[1, 0].legend(fontsize=8)
    T_add = np.linspace(0.0, 6.0e6, 200)
    for name, (aa, bb) in results.items():
        axs[1, 1].plot(T_add / M, [prob_rate_below(aa, bb + t, lam_req) for t in T_add], label=name.split(" (")[0])
    axs[1, 1].axhline(0.9, color="gray", lw=0.5)
    axs[1, 1].set_xlabel("additional failure-free hours [10^6 h]")
    axs[1, 1].set_ylabel("P(lambda <= requirement)")
    axs[1, 1].set_title("test planning")
    axs[1, 1].legend(fontsize=8)
    for ax in axs.ravel():
        ax.grid(True, lw=0.3)
    fig.suptitle("Bayesian update of a failure rate (lesson parameters and synthetic data)")
    fig.tight_layout()
    fig.savefig(out_dir / "failure_rate_posterior.png", dpi=150)
    print(f"\nsaved {out_dir / 'failure_rate_posterior.png'}")


if __name__ == "__main__":
    main()
