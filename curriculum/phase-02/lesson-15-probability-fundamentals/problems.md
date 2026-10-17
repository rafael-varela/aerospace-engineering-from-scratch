# P02.15  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. Give every probability to four significant figures, state whether it is exact or simulated, and give the standard error of every simulated value.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: three dice by counting)

Three fair dice are rolled. Let T = 7 + (s mod 6), a target sum between 7 and 12.

(a) Write the size of the sample space. (b) By counting, find P(sum = T), P(at least one six) and P(all three faces different). Use the complement rule where it helps. (c) Are "all three faces different" and "at least one six" independent? Decide by computing P(both) and comparing with the product. (d) Confirm every number with `enumerate_dice(3)` and `event_probability`.

Answer form: four exact fractions and decimals, one independence verdict with the two numbers that decide it.

## 2. (Derivation: the rules from the axioms)

(a) Starting from the definition P(A | B) = P(A and B)/P(B), derive the law of total probability for a partition {A₁, …, Aₘ} and then Bayes' rule, stating where you use each axiom.
(b) Derive P(h heads in n flips of a coin with P(head) = p) = C(n, h) pʰ(1 − p)ⁿ⁻ʰ, naming the rule used at each step (independence, mutual exclusivity, counting).
(c) Prove that these probabilities sum to 1 over h = 0, …, n.
(d) Show that the probability of at least k of n components working, Σⱼ₌ₖⁿ C(n, j) pʲ(1 − p)ⁿ⁻ʲ, equals the parallel formula 1 − (1 − p)ⁿ when k = 1 and the series formula pⁿ when k = n.

Answer form: four short derivations.

## 3. (Unfamiliar scenario: launch-window weather)

A launch needs two independent weather conditions on the same day: winds aloft within limits, probability p_w = 0.70 + s/200, and no lightning within 10 nautical miles, probability p_l = 0.90. Days are independent of each other.

(a) What is the probability of a "go" on a given day? (b) What is the probability of launching on at least one of three consecutive days? (c) Day 1 is scrubbed. Draw the tree of the three ways a day can be scrubbed (wind only, lightning only, both) and compute the probability that the wind was at least part of the reason. (d) The two conditions are in reality correlated: thunderstorms bring strong winds aloft. Without any numbers, say whether the true probability of a go on a given day is higher or lower than (a) when the two bad conditions tend to occur together, and why.

Answer form: three probabilities, a tree diagram, one argued direction.

## 4. (Numerical, uses your code: a voting flight computer)

Each of several identical flight computers returns a correct answer with probability p = 0.90 + s/400 on a given cycle, independently.

(a) With `k_of_n_reliability`, compute the probability of a correct majority for a single computer, a 2-of-3 vote and a 3-of-5 vote. Express each improvement as "number of nines" (−log₁₀ of the failure probability).
(b) Estimate the 3-of-5 value by Monte Carlo with `monte_carlo_probability`, N = 20 000 trials and seed s. Report the estimate, its standard error, and the error in units of the standard error.
(c) Run `monte_carlo_error_study([200, 800, 3200, 12800], n_repeats=20, seed=s)` for the dice-sum problem and fit the log-log slope. State the slope and whether it is consistent with −½ given that each RMS value is itself estimated from only 20 repeats.
(d) For which range of p does the 2-of-3 vote do *worse* than a single computer? Derive the condition and verify at one value of p inside the range.

Answer form: three exact probabilities with their nines, one simulated value with standard error, one slope with a judgement, one inequality with a check.

## 5. (Dispersion, uses your code: a different sounding)

A different payload descends through n = 10 + s layers, each pushing it d = 12 m east or west and north or south with equal probability.

(a) Give the per-axis standard deviation √n d and the largest possible offset. (b) With `radius_for_probability_exact`, find the smallest radius that contains 90 % and 99 % of landings. Compare with the Gaussian values σ√(−2 ln(1 − P)). (c) Simulate 20 000 descents with seed s and estimate the probability of landing within the 90 % radius from (b); report it with its standard error and say whether it is consistent with 0.90 (remember that the exact grid probability at that radius is at least 0.90, not exactly 0.90). (d) Plot the landing scatter with the two circles.

Answer form: two numbers, two radii with their Gaussian comparisons, one simulated probability with standard error, one plot.

## 6. (Judgment: where the model breaks)

(a) Suppose the 20 wind layers are perfectly correlated instead of independent (every layer pushes the same way). What is the per-axis spread of the landing offset then, and what are the only possible landing points? Compare with the independent case and explain which real feature of the atmosphere lies between the two.
(b) The Monte Carlo standard error formula √(p̂(1 − p̂)/N) is to be used to estimate a probability of about 10⁻⁵ (a landing beyond a 1 km exclusion zone). With N = 10⁵ trials, how many hits do you expect, and why is the formula then useless? What would you change?
(c) The inspection posterior P(cracked | flagged) was 0.269 for a 2 % prior. A supplier argues that because the inspection "is 90 % accurate" the flagged parts should all be scrapped. In three sentences, explain to the supplier what the number 0.269 means and which of the three inputs (prior, sensitivity, false-alarm rate) would have to change, and by how much, to make the posterior exceed 0.9.

Answer form: three short paragraphs, each with at least one number.
