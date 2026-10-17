# P02.15  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. The sample space has 6³ outcomes. "At least one six" is easiest through its complement. "All different" is a permutation count: 6 × 5 × 4. For independence, count the outcomes that are both all-different *and* contain a six (choose the six's position, then two different non-six faces in order).
2. Every rule is a rewriting of P(A and B) = P(A | B) P(B). For the binomial formula, first find the probability of *one specific* sequence with h heads, then count how many such sequences there are. For the sum, set a = p and b = 1 − p in the binomial theorem.
3. A "go" needs both conditions, so multiply; "at least one go in three days" is the complement of three scrubs; the scrub tree has three leaves and the question in (c) is a conditional probability over two of them. For (d), compare P(both bad) under independence with what positive correlation does to it.
4. `k_of_n_reliability(n, k, p)` for the exact values; "nines" is −log₁₀(1 − P). For the Monte Carlo, the sampler draws `(n_trials, 5)` uniforms and the event is "at least three below p". For (d), solve 3p² − 2p³ < p.
5. The grid positions are d(2h − n); the largest offset is n d along one axis. `radius_for_probability_exact` returns a grid radius, so the exact probability at it is ≥ P, not = P; compare the simulation with the *exact grid probability at that radius*, not with 0.90.
6. Perfectly correlated layers reduce 20 coin flips to one. For rare events, N p is the expected number of hits; think about what √(p̂(1 − p̂)/N) is when p̂ = 0. For the supplier, write the posterior odds as (prior odds) × (sensitivity / false-alarm rate) and see which factor can realistically move.
