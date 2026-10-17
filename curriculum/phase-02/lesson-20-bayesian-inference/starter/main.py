"""P02.20 Bayesian Inference: starter code.

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

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the update, writes failure_rate_posterior.png)
Test: pytest curriculum/phase-02/lesson-20-bayesian-inference
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
    raise NotImplementedError


def odds_update(prior_odds: float, bayes_factor: float) -> float:
    """Return the posterior odds = prior odds x Bayes factor (likelihood ratio of the data under the two hypotheses)."""
    raise NotImplementedError


# ---------------------------------------------------------------- likelihoods and densities (log-space, from scratch)

def log_binomial_likelihood(k: int, n: int, p: np.ndarray) -> np.ndarray:
    """Return log P(k failures in n demands | p) = log C(n, k) + k log p + (n - k) log(1 - p), elementwise in p."""
    raise NotImplementedError


def log_poisson_likelihood(k: int, T: float, lam: np.ndarray) -> np.ndarray:
    """Return log P(k failures in T hours | lambda) = k log(lambda T) - lambda T - log k!, elementwise in lambda."""
    raise NotImplementedError


def log_gamma_pdf(lam: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """Return log of the Gamma(alpha, beta) density in the rate parametrisation, elementwise; -inf for lambda <= 0.

    For beta = 0 (an improper prior such as Jeffreys' alpha = 1/2 or the flat alpha = 1) the
    normalising constant is dropped and (alpha - 1) log lambda is returned.
    """
    raise NotImplementedError


def log_beta_pdf(p: np.ndarray, a: float, b: float) -> np.ndarray:
    """Return log of the Beta(a, b) density, elementwise; -inf outside [0, 1]."""
    raise NotImplementedError


# ---------------------------------------------------------------- conjugate updates

def gamma_poisson_update(alpha: float, beta: float, k: int, T: float) -> tuple[float, float]:
    """Return the posterior (alpha + k, beta + T) of a Gamma(alpha, beta) prior after k failures in T hours."""
    raise NotImplementedError


def beta_binomial_update(a: float, b: float, k: int, n: int) -> tuple[float, float]:
    """Return the posterior (a + k, b + n - k) of a Beta(a, b) prior after k failures in n demands."""
    raise NotImplementedError


def sequential_update(alpha: float, beta: float, events: list[tuple[int, float]]) -> tuple[float, float]:
    """Apply ``gamma_poisson_update`` once per (k, T) record in ``events``, in order, and return the final (alpha, beta)."""
    raise NotImplementedError


# ---------------------------------------------------------------- grid posterior (from scratch)

def grid_posterior(grid: np.ndarray, log_prior: np.ndarray, log_likelihood: np.ndarray) -> np.ndarray:
    """Return the posterior density on ``grid`` from log prior and log likelihood values, normalised by the trapezoid rule.

    Works in log space (subtracting the maximum before exponentiating) so that likelihoods like
    lambda^k e^(-lambda T) with T ~ 10^6 h do not underflow. Raise ``ValueError`` if nothing is finite.
    """
    raise NotImplementedError


def grid_mean(grid: np.ndarray, pdf: np.ndarray) -> float:
    """Return the mean of a density tabulated on ``grid`` (trapezoid rule)."""
    raise NotImplementedError


def grid_cdf(grid: np.ndarray, pdf: np.ndarray) -> np.ndarray:
    """Return the cumulative distribution on ``grid`` (cumulative trapezoid), from 0 to 1."""
    raise NotImplementedError


def grid_quantile(grid: np.ndarray, pdf: np.ndarray, q: float) -> float:
    """Return the value x with CDF(x) = q by linear interpolation of the tabulated CDF."""
    raise NotImplementedError


def credible_interval(grid: np.ndarray, pdf: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    """Return the equal-tailed credible interval (lower, upper) that holds probability ``level``."""
    raise NotImplementedError


def hpd_interval(grid: np.ndarray, pdf: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    """Return the highest-posterior-density interval: the shortest interval holding probability ``level``.

    Found by lowering a density threshold until the mass above it reaches ``level``; assumes a
    unimodal density on the grid. Raise ``ValueError`` for a bad level.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- Metropolis sampler (from scratch)

def metropolis(log_target: Callable[[float], float], x0: float, n_samples: int, step: float,
               rng: np.random.Generator, burn_in: int = 1000) -> tuple[np.ndarray, float]:
    """Return (samples, acceptance rate) from a random-walk Metropolis chain targeting exp(log_target).

    Proposal x' = x + step * N(0, 1); accept with probability min(1, exp(log_target(x') - log_target(x))).
    ``burn_in`` initial draws are discarded. A proposal with log_target = -inf is always rejected,
    which is how a support boundary (lambda > 0) is enforced.
    """
    raise NotImplementedError


def autocorrelation_time(samples: np.ndarray, max_lag: int = 200) -> float:
    """Return the integrated autocorrelation time 1 + 2 sum_k rho_k (summing until rho first drops below 0).

    Effective sample size = n / autocorrelation_time.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- what the posterior means for the mission

def prob_rate_below(alpha: float, beta: float, lam: float) -> float:
    """Return P(lambda <= lam) under Gamma(alpha, beta): the regularised lower incomplete gamma P(alpha, beta lam)."""
    raise NotImplementedError


def posterior_predictive_reliability(alpha: float, beta: float, t: float) -> float:
    """Return P(no failure in the next t hours) averaged over Gamma(alpha, beta): (beta / (beta + t))^alpha.

    E[e^(-lambda t)] is the Gamma moment-generating function at -t. It is always *above*
    e^(-mean t) because e^(-lambda t) is convex (Jensen's inequality): the chance that the rate is
    low helps more than the chance that it is high hurts, for the probability of zero failures.
    """
    raise NotImplementedError


def hours_to_demonstrate(alpha: float, beta: float, lam_req: float, conf: float) -> float:
    """Return the additional failure-free hours T such that P(lambda <= lam_req | data + T) = conf.

    P(alpha, (beta + T) lam_req) = conf gives (beta + T) lam_req = P^{-1}(alpha, conf); T = that / lam_req - beta,
    clipped at 0 when the requirement is already demonstrated.
    """
    raise NotImplementedError


def frequentist_upper_bound(k: int, T: float, conf: float) -> float:
    """Return the one-sided upper confidence bound on lambda: chi2_{2k+2}(conf) / (2T).

    It equals the Bayesian upper credible bound from the flat prior Gamma(1, 0), whose posterior is
    Gamma(k + 1, T) = chi2 with 2k + 2 degrees of freedom scaled by 1/(2T).
    """
    raise NotImplementedError


def zero_failure_tests_required(p_max: float, conf: float) -> int:
    """Return the smallest n with P(p <= p_max | 0 failures in n, flat prior) >= conf, i.e. 1 - (1 - p_max)^(n+1) >= conf.

    With the flat Beta(1, 1) prior the posterior after n successes is Beta(1, n + 1), whose CDF at
    p_max is 1 - (1 - p_max)^(n + 1).
    """
    raise NotImplementedError


def beta_cdf_grid(a: float, b: float, p_max: float, n_grid: int = 20001) -> float:
    """Return P(p <= p_max) under Beta(a, b) by integrating the tabulated density (a from-scratch incomplete beta).

    The density can be singular at p = 0 or p = 1 (a < 1 or b < 1, as for Jeffreys' Beta(1/2, 1/2)), so the
    grid is placed in u with p = sin^2(pi u / 2): the Jacobian dp/du = (pi/2) sin(pi u) vanishes at both
    ends and tames any singularity weaker than 1/p, after which the trapezoid rule converges.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- synthetic fleets and calibration

def simulate_fleet_failures(lam: float, T: float, rng: np.random.Generator) -> int:
    """Return a synthetic failure count: one Poisson draw with mean lambda T."""
    raise NotImplementedError


def credible_interval_coverage(alpha0: float, beta0: float, T: float, n_trials: int, level: float,
                               seed: int) -> float:
    """Return the fraction of synthetic fleets whose ``level`` equal-tailed posterior interval contains the true rate.

    Each trial draws lambda from the prior Gamma(alpha0, beta0), a Poisson count in T hours, updates
    to Gamma(alpha0 + k, beta0 + T) and checks whether lambda lies inside the interval (computed with
    ``scipy``'s gamma quantiles). Averaged over draws from the prior, Bayesian intervals are calibrated:
    the fraction tends to ``level``.
    """
    raise NotImplementedError


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
