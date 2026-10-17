# P02.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result (exact, or simulated with its standard error), at least one independent check (limiting case, complement, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py`.

## 1. Three dice by counting ($s = 15$: T = 7 + 3 = 10)

**Assumptions.** Fair dice, independent rolls, so the 6³ = 216 ordered triples are equally likely and P(A) = |A|/216.

**Derivation.** P(sum = 10): the ordered triples with total 10 number 27 (the compositions of 10 into three parts between 1 and 6: six permutations each of {1,3,6}, {1,4,5}, {2,3,5}; three each of {2,2,6}, {2,4,4}, {3,3,4}; 18 + 9 = 27), so P = 27/216 = **0.1250**. P(at least one six) = 1 − (5/6)³ = 91/216 = **0.4213**. P(all different) = 6 × 5 × 4/216 = 120/216 = **0.5556**.

Independence: triples that are all different *and* contain a six: the six takes one of 3 positions, the other two faces are an ordered pair of distinct non-six faces, 5 × 4 = 20, so 60 triples and P(both) = 60/216 = 0.2778. The product is 0.5556 × 0.4213 = 0.2341 ≠ 0.2778, so the events are **not independent**: knowing all faces are different makes a six more likely (a triple with a repeated face has fewer distinct faces available to be a six).

**Check.** `enumerate_dice(3)` with the four predicates gives 0.125, 0.421296, 0.555556 and 0.277778; `dice_sum_pmf(3)` has 27/216 at sum 10. Dimensionless throughout.

**Common mistakes.** Counting unordered sets {1,3,6} once instead of six times; using (5/6)³ as P(at least one six) instead of its complement; concluding independence because the events "feel unrelated".

## 2. The rules from the axioms

**Derivation.** (a) B = ∪ᵢ (B and Aᵢ) because the Aᵢ cover S; the pieces are mutually exclusive because the Aᵢ are. Axiom (iii), applied m − 1 times, gives P(B) = Σᵢ P(B and Aᵢ). The definition of conditional probability, rearranged to P(B and Aᵢ) = P(B | Aᵢ) P(Aᵢ) (valid when P(Aᵢ) > 0; terms with P(Aᵢ) = 0 contribute nothing either way), gives the law of total probability. For Bayes: P(Aᵢ | B) = P(Aᵢ and B)/P(B) by definition; the numerator is P(B | Aᵢ)P(Aᵢ) by the multiplication rule and the denominator is the sum just derived.

(b) Fix one sequence with h heads in n flips. By independence (multiplication rule applied n − 1 times) its probability is pʰ(1 − p)ⁿ⁻ʰ, regardless of *which* positions hold the heads. The event "exactly h heads" is the union of all such sequences; they are mutually exclusive (different sequences are different outcomes), and there are C(n, h) of them (choose the h positions). Axiom (iii): P(h) = C(n, h) pʰ(1 − p)ⁿ⁻ʰ.

(c) The binomial theorem states (a + b)ⁿ = Σₕ C(n, h) aʰ bⁿ⁻ʰ. With a = p, b = 1 − p the right side is Σₕ P(h) and the left side is 1ⁿ = 1.

(d) k = 1: Σⱼ₌₁ⁿ P(j) = 1 − P(0) = 1 − C(n, 0) p⁰(1 − p)ⁿ = 1 − (1 − p)ⁿ, the parallel formula. k = n: only j = n survives, C(n, n) pⁿ(1 − p)⁰ = pⁿ, the series formula.

**Check.** For n = 3, p = 0.98: 1 − 0.02³ = 0.999992 and 0.98³ = 0.941192, matching `parallel_reliability` and `series_reliability` in `main()`.

**Common mistakes.** Forgetting that the law of total probability needs a *partition* (both exhaustive and mutually exclusive); claiming the sequences in (b) are independent (they are mutually exclusive, the opposite).

## 3. Launch-window weather ($s = 15$: p_w = 0.775, p_l = 0.90)

**Assumptions.** The two conditions are independent on a day, days are independent, the probabilities are constant over the three days.

**Derivation.** (a) p_go = 0.775 × 0.90 = **0.6975**. (b) P(at least one go in three days) = 1 − (1 − 0.6975)³ = 1 − 0.3025³ = 1 − 0.02768 = **0.9723**. (c) The scrub tree: wind only, 0.225 × 0.90 = 0.2025; lightning only, 0.775 × 0.10 = 0.0775; both, 0.225 × 0.10 = 0.0225; total 0.3025 = 1 − p_go ✓. P(wind involved | scrub) = (0.2025 + 0.0225)/0.3025 = **0.7438**.

(d) With fixed marginals P(wind bad) = 0.225 and P(lightning) = 0.10, positive correlation raises P(both bad) above 0.0225. Since P(wind bad) = P(both bad) + P(wind only), P(wind only) must *fall* by the same amount, and so must P(lightning only); P(any scrub) = both + wind-only + lightning-only therefore falls by that amount, and P(go) *rises*. The bad weather bunches up on the same days, leaving more clean days. (Negative correlation, bad conditions avoiding each other, would do the opposite.)

**Check.** The three leaves plus p_go sum to 1 (0.6975 + 0.3025). Limiting case: p_l = 1 gives p_go = p_w and P(wind | scrub) = 1.

**Common mistakes.** Adding p_w and p_l; computing P(wind only | scrub) = 0.669 and forgetting the "both" leaf; asserting that correlation always makes things worse.

## 4. A voting flight computer ($s = 15$: p = 0.9375)

**Result.** (a) Single: 0.9375 (1.20 nines). 2-of-3: 3p²(1 − p) + p³ = **0.98877** (1.95 nines). 3-of-5: **0.99778** (2.65 nines). Each added pair of units buys about 0.7 of a nine at this p.

(b) `monte_carlo_probability` with a `(20 000, 5)` uniform draw and the event "at least three entries below p", seed 15: p̂ = **0.99785 ± 0.00033**, error +0.00007 = 0.21 SE.

(c) `monte_carlo_error_study([200, 800, 3200, 12800], n_repeats=20, seed=15)` gives RMS errors 0.02593, 0.01670, 0.00784, 0.00384 and a log-log slope of **−0.468**. Each RMS value is estimated from 20 repeats and carries about 1/√40 = 16 % relative scatter, which propagates to roughly ±0.1 in a four-point slope, so −0.47 is consistent with −0.5. The predicted RMS values √(p(1 − p)/N) are 0.0264, 0.0132, 0.0066, 0.0033; the measured ones scatter around them within the expected 16 %.

(d) 2-of-3 beats a single unit when 3p² − 2p³ > p, i.e. p(2p − 1)(p − 1) < 0, which for 0 < p < 1 holds exactly when p > ½. For p < ½ the vote is worse: at p = 0.4, `k_of_n_reliability(3, 2, 0.4)` = 0.352 < 0.4 (and at p = 0.6 it is 0.648 > 0.6).

**Check.** At p = ½ the vote returns exactly ½ (the test asserts it). `scipy.stats.binom.sf(2, 5, 0.9375)` = 0.99778.

**Common mistakes.** Using the number of *units* as the exponent in the nines calculation; reporting the Monte Carlo error without the standard error; fitting the slope with a linear rather than log-log fit.

## 5. A different sounding ($s = 15$: n = 25 layers, d = 12 m)

**Assumptions.** Independent fair coins per layer and axis; the model is synthetic.

**Result.** (a) σ = √25 × 12 = **60.0 m** per axis; the largest possible offset is n d = 300 m along one axis (all 25 layers the same way), 424 m on the diagonal. (b) `radius_for_probability_exact(25, 12, 0.9)` = **132.5 m** (Gaussian 60√(−2 ln 0.1) = 128.8 m, 2.9 % lower); for 99 %, **180.4 m** (Gaussian 182.1 m, 0.9 % higher). The grid probabilities at those radii are 0.9132 and 0.9909, not exactly 0.90 and 0.99, because probability accumulates in jumps as the circle takes in whole rings of grid points. (c) 20 000 simulated descents with seed 15: P(r ≤ 132.5 m) = **0.9162 ± 0.0020**, against the exact grid value 0.9132: a 1.5 SE difference, consistent. Against 0.90 it would look 8 SE high, which is why the comparison must be made against the grid value. (d) The plot shows a square-ish lattice of landing points (spacing 24 m) with the two circles; about 1 in 12 points lies outside the inner circle and about 1 in 110 outside the outer one.

**Check.** The simulated per-axis standard deviation is within 1 % of 60 m. The Gaussian and exact 99 % radii differ in the opposite sense from the 90 % radii: at 25 layers the discreteness error has the same size as the tail-lightness error and the two trade places between the two quantiles.

**Common mistakes.** Comparing the simulation to the nominal probability instead of the grid probability; treating n d as the "typical" offset instead of the maximum; forgetting that an odd n has no landing point at the origin.

## 6. Where the model breaks

(a) If all 20 layers push the same way, one coin per axis decides everything: the offset is ±n d = ±300 m on each axis, so the landing is at one of four corners 424 m from the aim point, with spread n d = 300 m against √n d = 67 m for independent layers, a factor √20 = 4.5. The real atmosphere lies between: gusts have a correlation length of a few hundred metres, so the 50 m layers come in correlated groups of perhaps 4–8, and the effective number of independent layers is 20 divided by the group size. With 4 independent groups of 5 the spread would be √4 × 5 × 15 = 150 m. Lesson 19 replaces the coins by a gust spectrum that carries this correlation length.

(b) Expected hits N p = 10⁵ × 10⁻⁵ = 1. The outcome is 0, 1, 2 or 3 hits, so p̂ is 0, 1, 2 or 3 × 10⁻⁵, and √(p̂(1 − p̂)/N) is either 0 (meaningless) or about equal to p̂ itself: a 100 % relative error. The rule of thumb is N p ≥ 10, so N = 10⁶ at least for a 30 % relative error and 10⁸ for 1 %. Better: compute the tail exactly if the model allows (for the coin-flip walk it does, with (10)), or use importance sampling that biases the coins outward and reweights the hits.

(c) To the supplier: "0.269 is the probability that a flagged part is actually cracked, given that only 2 % of parts are cracked to begin with; 73 % of the parts your inspection flags are good, because 5 % false alarms on 98 % good parts (4.9 % of all parts) outnumber 90 % detections on 2 % cracked ones (1.8 %). The posterior odds are the prior odds (0.02/0.98 = 1:49) times the likelihood ratio 0.90/0.05 = 18, i.e. 18:49 ≈ 0.37, or 27 %. To exceed 0.9 (odds 9:1) you would need the prior odds to reach 1:2 (a 33 % defect rate, which would mean a far worse batch, not a better test: `inspection_tree(0.3334, 0.9, 0.05)` gives 0.9000) or the false-alarm rate to fall to 0.90/(9 × 49) = 0.0020 (a 25-fold improvement; `inspection_tree(0.02, 0.9, 0.002)` gives 0.9018). Raising the sensitivity cannot do it alone: even at 1.0 the posterior is 0.02/(0.02 + 0.049) = 0.29." Scrapping all flagged parts throws away three good parts for each cracked one; re-inspecting flagged parts with an independent second method is the usual answer.

## Check-yourself answers

1. (i) P(A) ≥ 0; (ii) P(S) = 1; (iii) P(A or B) = P(A) + P(B) for mutually exclusive A, B. A and not-A are mutually exclusive with union S, so P(A) + P(not A) = 1.
2. P(A)P(B) = 0.15 = P(A and B): independent. Not mutually exclusive (P(A and B) ≠ 0). P(A or B) = 0.3 + 0.5 − 0.15 = 0.65.
3. Mutually exclusive means P(A and B) = 0, but independence requires P(A and B) = P(A)P(B) > 0 when both probabilities are positive; the two cannot hold at once.
4. P(B) = P(B | A₁)P(A₁) + P(B | A₂)P(A₂); P(flag) = 0.90 × 0.02 + 0.05 × 0.98 = 0.018 + 0.049 = 0.067.
5. About 1 %: posterior = 0.99 × 10⁻⁴/(0.99 × 10⁻⁴ + 0.01 × 0.9999) = 0.0098. The false alarms (1 % of the 99.99 % good parts) are a hundred times more numerous than the true detections.
6. SE = √(0.25/N): N = 0.25/10⁻⁶ = 250 000 for 0.001 and 2.5 × 10⁷ for 0.0001 (a hundred times more for ten times the precision).
7. Offsets d(2h − n) for h = 0, …, n, i.e. −n d, −(n − 2)d, …, n d in steps of 2d (n + 1 values). The largest, n d, has probability 1/2ⁿ (9.5 × 10⁻⁷ for n = 20).
8. Above: 0.6897 against 0.6753. The binomial endpoint has bounded support (|x| ≤ n d) and lighter tails than the Gaussian, so more of its probability lies inside any central region; the coarse grid adds a further, radius-dependent step.
