"""P02.17 Statistics & Hypothesis Testing: starter code.

Is this material batch in spec? Sample statistics, confidence intervals, one- and two-sample
t-tests, a proportion test, Type I and Type II errors, power by Monte Carlo and by the noncentral
t distribution, sample-size planning, a percentile bootstrap, a permutation test and one-sided
tolerance bounds (the "B-basis"-style allowable). Test statistics, p-value logic, bootstrap,
permutation and the power Monte Carlo are written from scratch; the Student-t and noncentral-t
quantiles come from ``scipy.stats`` and ``scipy.stats.ttest_1samp`` / ``ttest_ind`` are the oracles
in the tests.

All coupon data in this module are SYNTHETIC: they are drawn from stated normal distributions with
stated seeds, not measured on any real material.

Conventions:
    * ``alternative`` is one of 'two-sided', 'greater' (H1: mean > mu0) or 'less' (H1: mean < mu0).
    * Sample variance uses the n - 1 divisor.
    * Every random function takes ``seed`` and builds ``numpy.random.default_rng(seed)`` itself.
    * Strengths are in MPa.

Complete every function marked NotImplementedError.

Run:  python main.py          (prints the acceptance tables, writes acceptance_testing.png)
Test: pytest curriculum/phase-02/lesson-17-statistics-hypothesis-testing
"""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
from scipy import stats

ALTERNATIVES = ("two-sided", "greater", "less")


def _check_alternative(alternative: str) -> None:
    if alternative not in ALTERNATIVES:
        raise ValueError(f"alternative must be one of {ALTERNATIVES}")


# ---------------------------------------------------------------- sample statistics

def sample_mean(x: np.ndarray) -> float:
    """Return the arithmetic mean of a sample."""
    raise NotImplementedError


def sample_variance(x: np.ndarray) -> float:
    """Return the sample variance with the n - 1 divisor, by the two-pass formula sum (x - mean)^2 / (n - 1)."""
    raise NotImplementedError


def standard_error(x: np.ndarray) -> float:
    """Return s / sqrt(n): the estimated standard deviation of the sample mean."""
    raise NotImplementedError


def synthetic_batch(n: int, mean: float, sd: float, seed: int) -> np.ndarray:
    """Return ``n`` SYNTHETIC coupon strengths [MPa] drawn from N(mean, sd^2) with the given seed."""
    raise NotImplementedError


# ---------------------------------------------------------------- confidence intervals

def t_critical(conf: float, df: float) -> float:
    """Return the two-sided critical value t* with P(|T_df| <= t*) = conf (from scipy.stats.t)."""
    raise NotImplementedError


def confidence_interval_mean(x: np.ndarray, conf: float = 0.95) -> tuple[float, float]:
    """Return the two-sided ``conf`` confidence interval for the population mean: mean -/+ t* s / sqrt(n)."""
    raise NotImplementedError


# ---------------------------------------------------------------- hypothesis tests

def p_value_from_t(t: float, df: float, alternative: str) -> float:
    """Return the p-value of a t statistic under H0 for the given alternative.

    two-sided: 2 P(T >= |t|); greater: P(T >= t); less: P(T <= t).
    """
    raise NotImplementedError


def one_sample_t_test(x: np.ndarray, mu0: float, alternative: str = "two-sided") -> dict[str, float]:
    """Test H0: mean = mu0 against the alternative. Return {'t', 'df', 'p', 'mean', 'se'}.

    t = (mean - mu0) / (s / sqrt n) with n - 1 degrees of freedom. Assumes the observations are
    independent draws from a normal population (or n large enough for the CLT).
    """
    raise NotImplementedError


def welch_t_test(x: np.ndarray, y: np.ndarray, alternative: str = "two-sided") -> dict[str, float]:
    """Two-sample Welch t-test of H0: mean_x = mean_y, without assuming equal variances.

    t = (mean_x - mean_y) / sqrt(s_x^2/n_x + s_y^2/n_y); the Welch-Satterthwaite degrees of freedom are
    (s_x^2/n_x + s_y^2/n_y)^2 / [(s_x^2/n_x)^2/(n_x - 1) + (s_y^2/n_y)^2/(n_y - 1)].
    Return {'t', 'df', 'p', 'diff', 'se'} where diff = mean_x - mean_y.
    """
    raise NotImplementedError


def proportion_z_test(k: int, n: int, p0: float, alternative: str = "two-sided") -> dict[str, float]:
    """Normal-approximation test of H0: proportion = p0 from k successes in n trials.

    z = (k/n - p0) / sqrt(p0 (1 - p0) / n); valid when n p0 and n (1 - p0) are both at least about 10.
    Return {'z', 'p', 'p_hat'}. The tests compare it with the exact binomial test.
    """
    raise NotImplementedError


def decide(p_value: float, alpha: float = 0.05) -> str:
    """Return 'reject H0' if p < alpha, else 'fail to reject H0'. Raise ``ValueError`` for alpha outside (0, 1)."""
    raise NotImplementedError


# ---------------------------------------------------------------- errors, power, sample size

def type_i_error_rate_mc(mu0: float, sigma: float, n: int, alpha: float = 0.05, n_trials: int = 4000,
                         seed: int = 0, alternative: str = "two-sided") -> float:
    """Return the fraction of SYNTHETIC samples drawn under H0 (mean = mu0) that the t-test rejects.

    Should equal ``alpha`` to within the Monte Carlo standard error sqrt(alpha (1 - alpha) / n_trials).
    """
    raise NotImplementedError


def power_monte_carlo(delta: float, sigma: float, n: int, alpha: float = 0.05, n_trials: int = 4000,
                      seed: int = 0, alternative: str = "greater", mu0: float = 0.0) -> float:
    """Return the fraction of SYNTHETIC samples with true mean mu0 + delta that the t-test of H0: mean = mu0 rejects.

    Each trial draws n values from N(mu0 + delta, sigma^2) and runs ``one_sample_t_test``. With delta = 0
    this is the Type I error rate; otherwise it is the power 1 - beta.
    """
    raise NotImplementedError


def power_analytic(delta: float, sigma: float, n: int, alpha: float = 0.05, alternative: str = "greater") -> float:
    """Return the power of the one-sample t-test for a true shift ``delta`` from the noncentral t distribution.

    Under H1 the statistic follows a noncentral t with df = n - 1 and noncentrality delta sqrt(n) / sigma.
    """
    raise NotImplementedError


def required_sample_size(delta: float, sigma: float, alpha: float = 0.05, power: float = 0.9,
                         alternative: str = "greater", n_max: int = 10_000) -> int:
    """Return the smallest n >= 2 whose analytic power to detect ``delta`` reaches ``power``."""
    raise NotImplementedError


# ---------------------------------------------------------------- resampling

def bootstrap_ci_mean(x: np.ndarray, n_boot: int = 5000, conf: float = 0.95, seed: int = 0) -> tuple[float, float]:
    """Return the percentile-bootstrap confidence interval for the mean.

    Resample ``x`` with replacement ``n_boot`` times, compute each mean, and return the
    (1 - conf)/2 and (1 + conf)/2 quantiles. Makes no normality assumption.
    """
    raise NotImplementedError


def permutation_test(x: np.ndarray, y: np.ndarray, n_perm: int = 5000, seed: int = 0) -> float:
    """Return the two-sided permutation p-value for H0: x and y come from the same distribution.

    Statistic: difference of means. Pool the data, shuffle ``n_perm`` times, split into groups of the
    original sizes and count how often |difference| is at least the observed one; p = (count + 1) / (n_perm + 1).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- tolerance bounds

def one_sided_tolerance_factor(n: int, p: float = 0.90, conf: float = 0.95) -> float:
    """Return k such that mean - k s is, with confidence ``conf``, below the ``p`` fraction of a normal population.

    k = t_nc(conf; df = n - 1, nc = z_p sqrt n) / sqrt n with t_nc the noncentral-t quantile and z_p the
    standard-normal quantile of p. With p = 0.90 and conf = 0.95 this is the "B-basis" style factor.
    """
    raise NotImplementedError


def lower_tolerance_bound(x: np.ndarray, p: float = 0.90, conf: float = 0.95) -> float:
    """Return mean - k s: the one-sided lower tolerance bound covering fraction ``p`` with confidence ``conf``."""
    raise NotImplementedError


def tolerance_bound_coverage_mc(n: int, p: float = 0.90, conf: float = 0.95, n_trials: int = 4000,
                                seed: int = 0) -> float:
    """Return the fraction of SYNTHETIC N(0, 1) samples whose lower tolerance bound lies below the true p-quantile.

    Should equal ``conf`` within the Monte Carlo standard error: that is what "confidence" means here.
    """
    raise NotImplementedError


# ---------------------------------------------------------------- the lesson's numbers

def main(out_dir: str | Path = ".") -> None:
    """Print the acceptance-test tables and write acceptance_testing.png."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(out_dir)
    np.set_printoptions(precision=2, suppress=True)
    spec_min, claimed, sigma_hist, alpha = 440.0, 470.0, 12.0, 0.05

    print("== the batch (SYNTHETIC coupon data) ==")
    batch = synthetic_batch(8, 452.0, 12.0, seed=17)
    print(f"8 coupons [MPa]: {batch}")
    m, s, se = sample_mean(batch), math.sqrt(sample_variance(batch)), standard_error(batch)
    print(f"mean {m:.2f} MPa, sample sd {s:.2f} MPa, standard error {se:.2f} MPa")
    lo, hi = confidence_interval_mean(batch, 0.95)
    print(f"95 % CI for the batch mean: ({lo:.2f}, {hi:.2f}) MPa, half-width {0.5 * (hi - lo):.2f} MPa (t* = {t_critical(0.95, 7):.3f})")
    blo, bhi = bootstrap_ci_mean(batch, 5000, 0.95, seed=17)
    print(f"percentile bootstrap 95 % CI: ({blo:.2f}, {bhi:.2f}) MPa")

    print("\n== question 1: is the batch mean above the 440 MPa requirement? (H0: mu = 440, H1: mu > 440) ==")
    r1 = one_sample_t_test(batch, spec_min, "greater")
    print(f"t = {r1['t']:.3f}, df = {r1['df']:.0f}, p = {r1['p']:.5f} -> {decide(r1['p'], alpha)}")
    print(f"scipy.stats.ttest_1samp: t = {stats.ttest_1samp(batch, spec_min, alternative='greater').statistic:.3f}, "
          f"p = {stats.ttest_1samp(batch, spec_min, alternative='greater').pvalue:.5f}")

    print("\n== question 2: does the batch match the supplier's claimed 470 MPa? (two-sided) ==")
    r2 = one_sample_t_test(batch, claimed, "two-sided")
    print(f"t = {r2['t']:.3f}, p = {r2['p']:.5f} -> {decide(r2['p'], alpha)}")
    reference = synthetic_batch(10, 470.0, 12.0, seed=1701)
    r3 = welch_t_test(batch, reference, "two-sided")
    print(f"against 10 reference coupons (mean {sample_mean(reference):.2f}): Welch t = {r3['t']:.3f}, df = {r3['df']:.2f}, p = {r3['p']:.5f}")
    print(f"scipy ttest_ind(equal_var=False): p = {stats.ttest_ind(batch, reference, equal_var=False).pvalue:.5f}; "
          f"permutation test p = {permutation_test(batch, reference, 5000, seed=17):.4f}")

    print("\n== question 3: rivet defect rate (H0: p = 0.02, H1: p > 0.02) ==")
    r4 = proportion_z_test(5, 120, 0.02, "greater")
    print(f"5 defective of 120: p_hat = {r4['p_hat']:.4f}, z = {r4['z']:.3f}, p = {r4['p']:.4f}; "
          f"exact binomial p = {stats.binomtest(5, 120, 0.02, alternative='greater').pvalue:.4f} -> {decide(r4['p'], alpha)}")

    print("\n== Type I error, power and sample size (SYNTHETIC) ==")
    a_hat = type_i_error_rate_mc(spec_min, sigma_hist, 8, alpha, 20000, seed=170, alternative="greater")
    print(f"Type I error rate by simulation (20000 batches at exactly 440): {a_hat:.4f} (alpha = {alpha}, SE {math.sqrt(alpha * (1 - alpha) / 20000):.4f})")
    for delta in (5.0, 10.0, 15.0, 20.0):
        pa = power_analytic(delta, sigma_hist, 8, alpha, "greater")
        pm = power_monte_carlo(delta, sigma_hist, 8, alpha, 20000, seed=171, alternative="greater")
        print(f"  true margin of {delta:4.1f} MPa above 440, n = 8: power analytic {pa:.4f}, Monte Carlo {pm:.4f}")
    for target in (0.8, 0.9, 0.95):
        print(f"  coupons needed to demonstrate a 10 MPa margin with power {target:.2f}: n = {required_sample_size(10.0, sigma_hist, alpha, target, 'greater')}")

    print("\n== tolerance bound (B-basis style: 90 % of the population with 95 % confidence) ==")
    k8 = one_sided_tolerance_factor(8, 0.90, 0.95)
    print(f"k for n = 8: {k8:.3f}; for n = 30: {one_sided_tolerance_factor(30, 0.90, 0.95):.3f}; "
          f"for n = 1000: {one_sided_tolerance_factor(1000, 0.90, 0.95):.3f} (z_0.90 = {stats.norm.ppf(0.9):.3f})")
    print(f"lower tolerance bound of the batch: {lower_tolerance_bound(batch, 0.90, 0.95):.1f} MPa "
          f"(mean - z_0.90 s = {m - stats.norm.ppf(0.9) * s:.1f} MPa would ignore the uncertainty in mean and s)")
    print(f"coverage check by simulation, n = 8: {tolerance_bound_coverage_mc(8, 0.90, 0.95, 20000, seed=172):.4f} (target 0.95)")

    # Figure: sampling distribution of t under H0 and H1, power curve, CI coverage.
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.0))
    tt = np.linspace(-4, 8, 400)
    axes[0].plot(tt, stats.t.pdf(tt, 7), label="t under H0 (df = 7)")
    nc = 10.0 * math.sqrt(8) / sigma_hist
    axes[0].plot(tt, stats.nct.pdf(tt, 7, nc), label="t under H1 (true mean 450)")
    tc = stats.t.ppf(1 - alpha, 7)
    axes[0].axvline(tc, color="k", ls="--", label=f"critical t = {tc:.2f}")
    axes[0].axvline(r1["t"], color="C3", label=f"observed t = {r1['t']:.2f}")
    axes[0].set_xlabel("t statistic")
    axes[0].set_ylabel("density")
    axes[0].legend(fontsize=7)
    axes[0].set_title("one-sided t-test, n = 8")
    deltas = np.linspace(0, 25, 51)
    for n in (4, 8, 16, 32):
        axes[1].plot(deltas, [power_analytic(d, sigma_hist, n, alpha, "greater") for d in deltas], label=f"n = {n}")
    axes[1].plot([5, 10, 15, 20], [power_monte_carlo(d, sigma_hist, 8, alpha, 20000, seed=171, alternative="greater") for d in (5.0, 10.0, 15.0, 20.0)],
                 "ko", label="Monte Carlo, n = 8")
    axes[1].axhline(alpha, color="gray", ls=":")
    axes[1].set_xlabel("true margin of the batch mean above 440 MPa [MPa]")
    axes[1].set_ylabel("power (probability of demonstrating compliance)")
    axes[1].legend(fontsize=8)
    axes[1].set_title("power curve, sigma = 12 MPa")
    rng = np.random.default_rng(173)
    for i in range(40):
        sample = rng.normal(452.0, 12.0, size=8)
        lo_i, hi_i = confidence_interval_mean(sample, 0.95)
        axes[2].plot([lo_i, hi_i], [i, i], color="C0" if lo_i <= 452.0 <= hi_i else "C3")
    axes[2].axvline(452.0, color="k", ls="--", label="true mean (synthetic)")
    axes[2].set_xlabel("95 % CI for the mean [MPa]")
    axes[2].set_ylabel("synthetic batch")
    axes[2].legend(fontsize=8)
    axes[2].set_title("40 batches of 8: about 2 intervals miss")
    fig.tight_layout()
    fig.savefig(out_dir / "acceptance_testing.png", dpi=150)
    print(f"saved {out_dir / 'acceptance_testing.png'}")


if __name__ == "__main__":
    main()
