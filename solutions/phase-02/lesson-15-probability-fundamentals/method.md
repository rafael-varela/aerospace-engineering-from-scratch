# P02.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. |S| = 6³. 2. P(sum = T): count ordered triples with that total (list them, or read it off `dice_sum_pmf(3)`). 3. P(≥ 1 six) = 1 − (5/6)³. 4. P(all different) = (6 × 5 × 4)/6³. 5. Independence: count triples that are all different *and* contain exactly one six (3 positions × 5 × 4 ordered non-six faces); compare that fraction with the product of the two probabilities. 6. `enumerate_dice(3)` with three predicates.

## 2.
(a) Write B as the disjoint union of (B and Aᵢ); axiom (iii) gives P(B) = Σ P(B and Aᵢ); the multiplication rule turns each term into P(B | Aᵢ)P(Aᵢ). Bayes: P(Aᵢ | B) = P(Aᵢ and B)/P(B), expand the numerator with the multiplication rule and the denominator with total probability. (b) One specific sequence with h heads has probability pʰ(1 − p)ⁿ⁻ʰ by independence; the sequences with exactly h heads are mutually exclusive and there are C(n, h) of them; add. (c) Binomial theorem with a = p, b = 1 − p gives (p + 1 − p)ⁿ = 1. (d) k = 1: 1 − P(0 successes). k = n: only the j = n term survives.

## 3.
(a) p_go = p_w p_l. (b) 1 − (1 − p_go)³. (c) Leaves: wind-only (1 − p_w) p_l, lightning-only p_w (1 − p_l), both (1 − p_w)(1 − p_l); they sum to 1 − p_go. P(wind involved | scrub) = (wind-only + both)/(1 − p_go). (d) Positive correlation raises P(both bad) above (1 − p_w)(1 − p_l) and so, since the marginals are fixed, raises P(both good) too.

## 4.
(a) `k_of_n_reliability(1, 1, p)`, `(3, 2, p)`, `(5, 3, p)`; nines = −log₁₀(1 − P). (b) Sampler `rng.random((n, 5)) < p`, event `row sum >= 3`; `monte_carlo_probability(..., 20_000, seed=s)`; error/SE. (c) `monte_carlo_error_study`, `np.polyfit(log N, log rms, 1)`; the slope's own uncertainty from 20 repeats is about 0.1. (d) 3p² − 2p³ − p = p(2p − 1)(p − 1) < 0 for 0 < p < ½; check with p = 0.4.

## 5.
(a) σ = √n d; max = n d. (b) `radius_for_probability_exact(n, 12.0, 0.9)` and `(…, 0.99)`; Gaussian σ√(−2 ln 0.1), σ√(−2 ln 0.01). (c) `random_walk_landings(n, 12.0, 20_000, seed=s)`; `probability_within_radius_mc(landings, R90)`; compare with `probability_within_radius_exact(n, 12.0, R90)` in SE units. (d) Scatter plot with two circles.

## 6.
(a) One flip decides all 20 layers: offsets ±n d only, spread n d; compare with √n d; real gust correlation lengths sit between. (b) Expected hits N p; SE formula gives zero at p̂ = 0 and is a poor approximation for N p < 10; increase N, or use importance sampling / a variance-reduction method, or compute the exact tail. (c) Posterior odds = prior odds × (sensitivity/false-alarm); evaluate which factor must change and by how much to reach odds 9:1.
