"""P03.17 Uncertainty Quantification: starter code.

Propagation of input uncertainty through the Breguet range equation by three routes written from
scratch: linear (first-order Taylor, GUM-style) propagation, Monte Carlo with its 1/sqrt(N) error
estimate, bootstrap intervals and antithetic variance reduction, and a Hermite polynomial-chaos
expansion (PCE) built by Gauss-Hermite quadrature (Golub-Welsch nodes from scratch) with moments and
Sobol sensitivity indices read off the coefficients. ``scipy.stats`` and ``numpy.polynomial.hermite_e``
are the oracles.

Conventions:
    * Uncertain inputs are independent Gaussians written as X_i = mu_i (1 + s_i xi_i) with xi_i standard
      normal; the vector xi has shape (N, d) with d = 3 (c_T, L/D, W_fuel).
    * Probabilists' Hermite polynomials He_n(x): He_0 = 1, He_1 = x, He_{n+1} = x He_n - n He_{n-1},
      orthogonal under the standard normal density with E[He_m He_n] = n! delta_mn.
    * A PCE of total degree p in d variables uses the multi-indices alpha with |alpha| <= p; its basis
      functions are products He_alpha(xi) = prod_i He_{alpha_i}(xi_i) with norm prod_i alpha_i!.
    * Every probability distribution in this lesson is an *assumed* input, not a measurement.

Given: the dataclass, ``range_model``, ``total_degree_multi_index`` and ``main``; everything else is yours.

Run:  python main.py          (prints the MC, linear and PCE tables, writes range_uq.png)
Test: pytest curriculum/phase-03/lesson-17-uncertainty-quantification
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

G0 = 9.80665


# ---------------------------------------------------------------- the model (GIVEN)


@dataclass(frozen=True)
class RangeInputs:
    """Means and relative standard deviations of the uncertain inputs of the Breguet range, plus the
    fixed quantities. The numbers describe the light twin of P03.16 at its optimised wing; the
    uncertainty levels are assumed engineering judgements, not measured statistics."""

    V: float = 80.0                # cruise speed [m/s], fixed
    W_zf: float = 29126.0          # zero-fuel weight [N], fixed
    c_T_mean: float = 1.4e-5       # thrust-specific fuel consumption [kg/(N s)]
    c_T_rel_std: float = 0.03      # 3 % (engine-to-engine and deterioration)
    LD_mean: float = 17.9          # cruise lift-to-drag ratio
    LD_rel_std: float = 0.05       # 5 % (drag prediction)
    W_fuel_mean: float = 4512.0    # usable fuel loaded [N]
    W_fuel_rel_std: float = 0.02   # 2 % (fuel gauging and density)
    R_required: float = 1.35e6     # the mission range that must be met [m]

    @property
    def means(self) -> np.ndarray:
        return np.array([self.c_T_mean, self.LD_mean, self.W_fuel_mean])

    @property
    def rel_stds(self) -> np.ndarray:
        return np.array([self.c_T_rel_std, self.LD_rel_std, self.W_fuel_rel_std])


def range_model(xi: np.ndarray, p: RangeInputs) -> np.ndarray:
    """Return the Breguet range R [m] for standard-normal inputs xi of shape (N, 3) or (3,).

    R = (V / (g c_T)) (L/D) ln((W_zf + W_fuel) / W_zf) with c_T = c_T_mean (1 + s1 xi1), L/D = LD_mean (1 + s2 xi2),
    W_fuel = W_fuel_mean (1 + s3 xi3). Nonlinear in c_T (1/x) and in W_fuel (logarithm), linear in L/D.
    """
    xi = np.atleast_2d(np.asarray(xi, dtype=float))
    X = p.means * (1.0 + p.rel_stds * xi)
    c_T, LD, W_fuel = X[:, 0], X[:, 1], X[:, 2]
    R = (p.V / (G0 * c_T)) * LD * np.log((p.W_zf + W_fuel) / p.W_zf)
    return R if R.size > 1 else R[0]


# ---------------------------------------------------------------- sampling from scratch


def standard_normals_box_muller(n: int, seed: int) -> np.ndarray:
    """Return n standard normal samples from uniforms by the Box-Muller transform.

    With U1, U2 uniform on (0, 1): Z1 = sqrt(-2 ln U1) cos(2 pi U2), Z2 = sqrt(-2 ln U1) sin(2 pi U2) are
    independent standard normals. Uses ``numpy.random.default_rng(seed).random`` for the uniforms only.
    """
    raise NotImplementedError


def inverse_cdf_sample(u: np.ndarray, inverse_cdf: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Return samples X = F^{-1}(U) by the inversion method: if U is uniform on (0, 1) then X has CDF F.
    Raise ``ValueError`` if any u is outside (0, 1)."""
    raise NotImplementedError


# ---------------------------------------------------------------- Monte Carlo


def monte_carlo(f: Callable[[np.ndarray], np.ndarray], d: int, n: int, seed: int) -> dict[str, float]:
    """Return {'mean', 'std', 'se_mean', 'samples'} from n evaluations of f at seeded standard normals.

    The standard error of the mean is s / sqrt(n) with s the sample standard deviation (ddof = 1): the
    Monte Carlo error falls as 1/sqrt(n) whatever the dimension d.
    """
    raise NotImplementedError


def mc_probability(samples: np.ndarray, threshold: float) -> tuple[float, float]:
    """Return (p, se) with p the fraction of samples below ``threshold`` and se = sqrt(p (1 - p) / n)
    its binomial standard error. A tail probability needs about 10 / p samples to be known to 30 %."""
    raise NotImplementedError


def bootstrap_ci(samples: np.ndarray, statistic: Callable[[np.ndarray], float], n_boot: int, seed: int,
                 alpha: float = 0.05) -> tuple[float, float]:
    """Return the percentile bootstrap (1 - alpha) confidence interval of ``statistic``.

    Resample the data with replacement n_boot times, compute the statistic each time, and take the
    alpha/2 and 1 - alpha/2 quantiles. Needs no formula for the statistic's standard error, which is why
    it is used for medians, quantiles and tail probabilities.
    """
    raise NotImplementedError


def antithetic_monte_carlo(f: Callable[[np.ndarray], np.ndarray], d: int, n: int, seed: int) -> dict[str, float]:
    """Return {'mean', 'se_mean'} from n evaluations in antithetic pairs (xi, -xi).

    The pair average (f(xi) + f(-xi)) / 2 removes every odd component of f, so for a nearly linear
    model the variance of the mean estimate falls far below s^2 / n. The standard error is computed
    from the n/2 independent pair averages.
    """
    raise NotImplementedError


def mc_error_study(f: Callable[[np.ndarray], np.ndarray], d: int, ns: np.ndarray, n_rep: int, seed: int,
                   exact_mean: float) -> np.ndarray:
    """Return the root-mean-square error of the Monte Carlo mean over ``n_rep`` independent replications
    for each sample size in ``ns``: the observed order in N should be -1/2."""
    raise NotImplementedError


# ---------------------------------------------------------------- linear (GUM) propagation


def linear_propagation(f: Callable[[np.ndarray], np.ndarray], d: int, h: float = 1e-4) -> dict[str, float | np.ndarray]:
    """Return the first-order (GUM) propagation of standard-normal inputs through f at xi = 0.

    Gradient g_i = df/dxi_i by central differences; u^2 = sum g_i^2 (independent inputs of unit
    variance); the mean is f(0) at first order, and the second-order mean correction is
    (1/2) sum d^2f/dxi_i^2 (the trace of the Hessian; the mean shifts when f is curved). Returns
    {'mean1', 'mean2', 'std', 'grad', 'contributions'} with contributions g_i^2 / u^2.
    f is evaluated on single-row arrays of shape (1, d), as every other function here evaluates it on (N, d).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- Hermite polynomials and quadrature


def hermite_e(n: int, x: np.ndarray) -> np.ndarray:
    """Return the probabilists' Hermite polynomial He_n(x) by the three-term recurrence
    He_{k+1} = x He_k - k He_{k-1} with He_0 = 1, He_1 = x. Raise ``ValueError`` for n < 0."""
    raise NotImplementedError


def gauss_hermite_e(n: int) -> tuple[np.ndarray, np.ndarray]:
    """Return n Gauss-Hermite nodes and weights for the *standard normal* weight (weights sum to 1).

    Golub-Welsch: the nodes are the eigenvalues of the symmetric tridiagonal Jacobi matrix of the
    recurrence (zero diagonal, off-diagonal sqrt(k), k = 1..n-1), and the weights are the squared first
    components of the normalised eigenvectors. The rule integrates polynomials of degree <= 2n - 1
    exactly against the normal density, so E[xi^2k] = (2k - 1)!! is reproduced for 2k <= 2n - 1.
    """
    raise NotImplementedError


def total_degree_multi_index(d: int, p: int) -> np.ndarray:
    """Return all multi-indices alpha in N^d with |alpha| <= p, ordered by total degree then
    lexicographically: the basis of a total-degree PCE. Count: (d + p)! / (d! p!)."""
    idx = [a for a in itertools.product(range(p + 1), repeat=d) if sum(a) <= p]
    idx.sort(key=lambda a: (sum(a), a))
    return np.array(idx, dtype=int)


def pce_basis(index: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Return the matrix Psi[n, k] = He_alpha_k(xi_n) = prod_i He_{alpha_k,i}(xi_n,i) for the multi-indices
    ``index`` (K x d) at the points ``xi`` (N x d)."""
    raise NotImplementedError


def pce_norms(index: np.ndarray) -> np.ndarray:
    """Return E[He_alpha^2] = prod_i alpha_i! for each multi-index."""
    raise NotImplementedError


def pce_projection(f: Callable[[np.ndarray], np.ndarray], d: int, p: int, n_quad: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (coefficients, index) of the total-degree-p Hermite PCE of f by non-intrusive spectral
    projection: c_alpha = E[f He_alpha] / E[He_alpha^2], with the expectation computed by a tensor-product
    Gauss-Hermite rule of n_quad points per dimension (n_quad^d model evaluations). Exact for polynomial f
    of degree <= p when n_quad >= p + 1 (the integrand f He_alpha has degree <= 2p <= 2 n_quad - 1)."""
    raise NotImplementedError


def pce_regression(f: Callable[[np.ndarray], np.ndarray], d: int, p: int, n_samples: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (coefficients, index) by point collocation: least-squares fit of the basis to f at
    n_samples seeded standard-normal points (needs n_samples >= 2 x number of basis terms to be
    well posed). Agrees with projection to the sampling error; it is the route when the model can only
    be evaluated at points you do not choose."""
    raise NotImplementedError


def pce_evaluate(coeffs: np.ndarray, index: np.ndarray, xi: np.ndarray) -> np.ndarray:
    """Return the PCE surrogate sum_alpha c_alpha He_alpha(xi) at the points xi (N x d)."""
    raise NotImplementedError


def pce_moments(coeffs: np.ndarray, index: np.ndarray) -> tuple[float, float]:
    """Return (mean, variance): mean = c_0 (the constant term), variance = sum_{alpha != 0} c_alpha^2 E[He_alpha^2],
    by orthogonality. No sampling involved."""
    raise NotImplementedError


def sobol_indices(coeffs: np.ndarray, index: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return (first-order S_i, total S_Ti) Sobol indices from the PCE coefficients.

    S_i sums the variance of the terms that depend on xi_i alone; S_Ti sums every term that involves xi_i.
    The first-order indices sum to 1 exactly when the model is additive (no interaction terms).
    """
    raise NotImplementedError


# ---------------------------------------------------------------- main


def main(out_dir: str | Path = ".") -> None:
    """Propagate the input uncertainty three ways, compare, and plot."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy import stats

    out_dir = Path(out_dir)
    p = RangeInputs()
    f = lambda xi: range_model(xi, p)
    d = 3
    R0 = float(range_model(np.zeros(3), p))
    print(f"ASSUMED input distributions: c_T {p.c_T_mean:.2e} +- {100 * p.c_T_rel_std:.0f} %, L/D {p.LD_mean} +- {100 * p.LD_rel_std:.0f} %,"
          f" W_fuel {p.W_fuel_mean:.0f} N +- {100 * p.W_fuel_rel_std:.0f} %; nominal range {R0 / 1e3:.1f} km, required {p.R_required / 1e3:.0f} km")

    lin = linear_propagation(f, d)
    print(f"\nlinear (GUM) propagation: mean {lin['mean1'] / 1e3:.2f} km (second order {lin['mean2'] / 1e3:.2f}), std {lin['std'] / 1e3:.2f} km,"
          f" variance shares {np.round(lin['contributions'], 3)}")

    print("\nMonte Carlo (seed 17):")
    for n in (100, 1000, 10_000, 100_000):
        mc = monte_carlo(f, d, n, seed=17)
        pr, se = mc_probability(mc["samples"], p.R_required)
        print(f"  N = {n:7d}: mean {mc['mean'] / 1e3:.2f} +- {mc['se_mean'] / 1e3:.2f} km, std {mc['std'] / 1e3:.2f} km,"
              f" P(R < required) = {pr:.4f} +- {se:.4f}")
    mc = monte_carlo(f, d, 10_000, seed=17)
    lo, hi = bootstrap_ci(mc["samples"], lambda y: float(np.mean(y < p.R_required)), n_boot=2000, seed=1)
    print(f"  bootstrap 95 % interval of the shortfall probability at N = 10^4: [{lo:.4f}, {hi:.4f}]")
    ant = antithetic_monte_carlo(f, d, 10_000, seed=17)
    print(f"  antithetic N = 10^4: mean {ant['mean'] / 1e3:.2f} +- {ant['se_mean'] / 1e3:.3f} km (plain +- {mc['se_mean'] / 1e3:.2f})")
    ns = np.array([100, 400, 1600, 6400])
    c2, _ = pce_projection(f, d, 6, 8)
    mean_ref, var_ref = pce_moments(c2, total_degree_multi_index(d, 6))
    rms = mc_error_study(f, d, ns, n_rep=200, seed=3, exact_mean=mean_ref)
    slope = np.polyfit(np.log(ns), np.log(rms), 1)[0]
    print(f"  RMS error of the MC mean over 200 replications: {np.round(rms / 1e3, 2)} km for N = {ns.tolist()}; slope {slope:.3f} (theory -0.5)")

    print("\npolynomial chaos (Hermite, tensor Gauss-Hermite with p + 1 points per dimension):")
    for deg in (1, 2, 3, 4):
        c, idx = pce_projection(f, d, deg, deg + 1)
        m, v = pce_moments(c, idx)
        print(f"  degree {deg}: {idx.shape[0]:3d} terms, {(deg + 1) ** d:4d} model runs, mean {m / 1e3:.3f} km, std {math.sqrt(v) / 1e3:.3f} km")
    c, idx = pce_projection(f, d, 4, 5)
    S, ST = sobol_indices(c, idx)
    print(f"  Sobol first-order (c_T, L/D, W_fuel): {np.round(S, 4)}, total: {np.round(ST, 4)}, sum of first-order {S.sum():.4f}")
    rng = np.random.default_rng(5)
    xi_big = rng.standard_normal((200_000, d))
    R_sur = pce_evaluate(c, idx, xi_big)
    pr_sur = float(np.mean(R_sur < p.R_required))
    print(f"  shortfall probability from 2 x 10^5 samples of the degree-4 surrogate: {pr_sur:.4f}"
          f" (normal approximation with the linear std: {stats.norm.cdf((p.R_required - lin['mean1']) / lin['std']):.4f})")
    c_reg, _ = pce_regression(f, d, 2, 400, seed=11)
    c_proj, _ = pce_projection(f, d, 2, 3)
    print(f"  regression vs projection, degree 2: max |difference| / c_0 = {np.max(np.abs(c_reg - c_proj)) / c_proj[0]:.2e}")

    # Figure.
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    ax = axes[0]
    ax.hist(mc["samples"] / 1e3, bins=60, density=True, alpha=0.5, label="Monte Carlo, N = 10$^4$")
    xs = np.linspace(mc["samples"].min(), mc["samples"].max(), 300)
    ax.plot(xs / 1e3, stats.norm.pdf(xs, lin["mean1"], lin["std"]) * 1e3, "k--", label="linear (GUM) normal")
    ax.hist(R_sur / 1e3, bins=120, density=True, histtype="step", color="C3", label="PCE surrogate, 2 x 10$^5$")
    ax.axvline(p.R_required / 1e3, color="C2", ls=":", label="required")
    ax.set_xlabel("range [km]"); ax.set_ylabel("density [1/km]"); ax.set_title("range distribution (assumed inputs)"); ax.legend(fontsize=7)
    ax = axes[1]
    ax.loglog(ns, rms / 1e3, "o-", label="RMS error of MC mean")
    ax.loglog(ns, rms[0] / 1e3 * np.sqrt(ns[0] / ns), "k--", label="N$^{-1/2}$")
    ax.set_xlabel("N"); ax.set_ylabel("error [km]"); ax.set_title("Monte Carlo convergence"); ax.legend(); ax.grid(lw=0.3, which="both")
    ax = axes[2]
    degs = np.sum(idx, axis=1)
    ax.semilogy(np.arange(len(c)), np.abs(c) / abs(c[0]) + 1e-17, "o", ms=4)
    for dd in range(1, 5):
        ax.axvline(np.argmax(degs == dd) - 0.5, color="0.7", lw=0.5)
    ax.set_xlabel("coefficient index (grouped by total degree 0..4)"); ax.set_ylabel("|c| / c$_0$"); ax.set_title("PCE coefficient decay"); ax.grid(lw=0.3)
    fig.suptitle("P03.17  Uncertainty in the Breguet range: MC, linear and polynomial chaos")
    fig.tight_layout()
    fig.savefig(out_dir / "range_uq.png", dpi=150)
    plt.close(fig)
    print(f"\nsaved {out_dir / 'range_uq.png'}")


if __name__ == "__main__":
    main()
