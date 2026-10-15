# P02.16  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Four thrusters: the pmf has five entries, C(4, k) pᵏ(1 − p)⁴⁻ᵏ. The quarter in (b) has three times the rate, and "no alert in a quarter" is also (P(no alert in a month))³: check that the two agree.
2. For (a) expand E[(aX + b − aμ − b)²]. For (d), the change to principal axes makes the two components independent standard normals *after* dividing by √λᵢ; the density is then a function of ρ² = z₁² + z₂² only and the area element is ρ dρ dθ.
3. Expected count = flux × area × days; a count of independent events at constant rate is Poisson; the time to the first event is exponential with rate = flux × area per day. For (d), a z-score is (observed − λ)/√λ.
4. Build `Sigma = np.diag([sw**2, 0.05**2, 0.3**2, 0.5**2])` and reuse `J` from `landing_sensitivities()`. The 99 % circle is `radius_for_probability(cov, 0.99)`; the area ratio is πR²/(πab). For the rotation, think about which input the off-diagonal covariance comes from.
5. Fit log(gap) against log(n) with `np.polyfit`. For the exponential sum, standardise by (S − n)/√n. The exponential is skewed; skewness decays as 1/√n and dominates the cdf gap, whereas the uniform's skewness is zero from the start.
6. (a) P(x > 1200) = 1 − Φ(1200/469.6). The t₄ distribution has a finite variance but tails that fall like a power law, so matched-σ tails beyond 2.5σ are *heavier*. (b) Σ₁₄ = Σ₄₁ = ρ σ_w σ_ψ; the extra covariance adds to Cov_xy through J₁₁J₂₄ and J₂₁J₂₄. (c) How many simulated landings fall outside a 99.9 % ellipse out of 20 000, and how many *real* flights would be needed to check that fraction?
