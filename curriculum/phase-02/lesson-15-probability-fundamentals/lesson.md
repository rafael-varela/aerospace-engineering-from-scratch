# P02.15  Probability Fundamentals

> A probability is a fraction of a sample space. Count it exactly when you can; simulate it when you cannot; and know the error of the simulation before you trust it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the addition rule, the multiplication rule, the law of total probability and Bayes' rule from the definition of conditional probability, and apply them to a series/parallel reliability chain and to an inspection tree. 2. Implement and test, from scratch, exact probability by enumeration and binomial coefficients, k-of-n reliability, Bayes' rule through a tree, and a seeded Monte Carlo estimator with its standard error, against `scipy.special.comb` and `scipy.stats.binom` as oracles. 3. Measure how the Monte Carlo error falls with the number of trials and show that it obeys the 1/√N law. 4. Decide, for a parachute payload drifting through random wind layers, what recovery radius captures 90 % and 99 % of landings, and say when the exact count and the continuous approximation disagree by more than the Monte Carlo error. |
| **Prerequisites** | `P00.03` (NumPy), `P00.06` (SciPy), `P00.07` (pytest), `P00.14` (reporting uncertainty), `P01.05` (series, for the binomial expansion). Skill check: you can expand (a + b)³, compute 5!/(3! 2!) by hand, and write a `for` loop that counts how often a condition holds. |
| **Engineering question** | A parachute payload released at 1000 m is pushed 15 m east or west and 15 m north or south by each of 20 wind layers, at random; how large a recovery circle contains 90 % of landings, how sure is a 20 000-run simulation of that number, and when is the bell-curve shortcut wrong? |
| **Assumptions and limits** | Equally likely outcomes where "classical" probabilities are used (fair dice, fair coins); independence between layers, components and trials where the multiplication rule is used; the coin-flip wind model is *synthetic* and deliberately crude (fixed drift per layer, no correlation between layers, no gust spectrum: Lesson 19 does that properly). Reliability numbers are illustrative inputs, not data on any real vehicle. Monte Carlo standard errors assume the trials are independent and N large enough that p̂ is roughly normal, which fails for rare events (p̂ N < 10). |
| **Analytical work** | Sample spaces and counting (permutations, combinations, the binomial coefficient); the three axioms; complement and addition rules; conditional probability, independence, the multiplication rule; the law of total probability and Bayes' rule derived from it; series, parallel and k-of-n reliability; the binomial distribution of coin-flip sums derived by counting; the standard error of a Monte Carlo proportion and the 1/√N law; the exact landing distribution of the coin-flip walk as a product of binomials. |
| **Computational deliverable** | `starter/main.py`: `n_choose_k`, `permutations_count`, `enumerate_dice`, `event_probability`, `dice_sum_pmf`, `complement`, `addition_rule`, `conditional_probability`, `are_independent`, `total_probability`, `bayes_posterior`, `inspection_tree`, `series_reliability`, `parallel_reliability`, `binomial_probability`, `k_of_n_reliability`, `proportion_standard_error`, `monte_carlo_probability`, `monte_carlo_dice_sum`, `monte_carlo_error_study`, `random_walk_landings`, `walk_endpoint_pmf`, `probability_within_radius_exact`, `probability_within_radius_mc`, `gaussian_radius_probability`, `radius_for_probability_exact`. Tests prove: binomial coefficients match SciPy and obey Pascal's rule; every pmf sums to one; enumeration reproduces hand counts; the rules reduce to their limiting cases (disjoint events, independent events, uninformative tests); k-of-n reduces to series and parallel at k = n and k = 1 and is monotone in p; Monte Carlo estimates lie within four standard errors of the exact value and their error falls as N^(−1/2); the coin-flip walk has the binomial endpoint distribution with variance n d² and the simulated recovery probability agrees with exact counting. |
| **Practical activity** | **computational experiment**: a Monte Carlo convergence study (error against N over five decades of N, twenty repeats each) and a sweep of the coin-flip wind model over layer counts to find where the continuous approximation becomes acceptable. All data are synthetic. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `probability_fundamentals.png`, the hand derivation of Bayes' rule and of the inspection posterior, the convergence table with its fitted slope, and a short note recommending a recovery radius with its uncertainty. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivations in part 2 are reproduced without the lesson open; all public tests pass; the hand values for the inspection tree and the 2-of-3 IMU match the code to four figures; the recovery-radius note states the radius, the method, the Monte Carlo standard error and the size of the Gaussian-approximation error. |
| **References** | `@openstax-introductory-statistics-2e` §3.1 "Terminology", §3.2 "Independent and Mutually Exclusive Events", §3.3 "Two Basic Rules of Probability", §3.5 "Tree and Venn Diagrams" (required-free); §3.4 "Contingency Tables" (supplementary-free). `@openstax-precalculus-2e` §11.5 "Counting Principles", §11.6 "Binomial Theorem", §11.7 "Probability" (required-free). `@mit-ocw-18-05-introduction-to-probability-and-statistics` Reading 1b "Counting and Sets", Reading 2 "Probability: Terminology and Examples", Reading 3 "Conditional Probability, Independence and Bayes' Theorem" (supplementary-free; Bayes' theorem is stated and named there, whereas Introductory Statistics reaches it only through tree diagrams in §3.5). |
| **Extensions** | ★ Replace the fair coin in each wind layer by a biased one (a prevailing wind) and derive the landing offset's mean and the recovery circle's new centre. ★ Let the drift per layer be drawn from a uniform distribution instead of ±d and show, by simulation, that the landing scatter is nearly unchanged when the variance per layer is matched (a preview of Lesson 16's central limit theorem). ★ Implement the inspection tree for a two-stage inspection (a second, independent test on flagged parts) and show that the posterior after both tests is the product of two likelihood ratios times the prior odds. |
| **Time** | unmeasured |

## 1. The engineering problem

A sounding-rocket payload comes down under a parachute from 1000 m. On the way it crosses 20 layers of air, 50 m thick, and in each one the wind pushes it sideways. Nobody knows the wind in each layer in advance; what is known, from past soundings, is that a layer's push is about 15 m and that its direction is as likely east as west, as likely north as south, independently of the layers above and below. The recovery team needs to know how big a circle around the aim point they must search, and the range-safety officer needs to know how likely the payload is to land outside the cleared zone. Both questions have the same form: *what fraction of all the ways the winds could blow puts the payload inside a given circle?* That fraction is a **probability**.

This lesson builds the vocabulary and the four rules that every later statistical lesson rests on, and it builds them twice: once by **counting** (writing out every equally likely outcome, or using the binomial coefficient when there are too many to write) and once by **simulation** (rolling the dice on a computer many times and counting how often the event occurs). The two must agree, and the simulation must come with an error bar. The same two tools answer the reliability questions of part 2 (what is the probability that a five-event flight sequence completes, or that a 2-of-3 inertial-measurement-unit vote returns a correct attitude) and the inspection question (if an eddy-current test flags a fastener hole, how likely is the hole really cracked).

Success criteria, written down first: enumeration reproduces the hand counts for two dice exactly; the k-of-n reliability formula reduces to the series and parallel formulas at k = n and k = 1; the Monte Carlo estimate of P(sum = 7) with N = 10⁶ trials lands within four standard errors (±0.0015) of 1/6; the root-mean-square Monte Carlo error over N = 100…25 600 falls with a log-log slope of −0.5 ± 0.1; and for the 20-layer wind model the simulated probability of landing within 100 m agrees with the exact count to within its standard error of 0.0033.

## 2. Understand and derive

### 2.1 Sample space, events, and the three axioms

An **experiment** is anything with an uncertain outcome: one roll of a die, one descent of the payload. The **sample space** S is the set of every possible outcome, and an **event** A is a subset of S. For two dice S has 36 ordered pairs (1,1), (1,2), …, (6,6); the event "the sum is 7" is {(1,6), (2,5), (3,4), (4,3), (5,2), (6,1)}.

A **probability** P assigns a number to every event such that (i) P(A) ≥ 0, (ii) P(S) = 1, and (iii) P(A or B) = P(A) + P(B) whenever A and B are **mutually exclusive** (they share no outcome). Those three statements are the axioms; everything else in the lesson is derived from them. When the outcomes are **equally likely** (fair dice, fair coins), the axioms force the classical formula

**P(A) = (number of outcomes in A) / (number of outcomes in S)**,   (1)

so P(sum = 7) = 6/36 = 0.166667, P(sum ≥ 10) = 6/36 and P(double) = 6/36. The code's `enumerate_dice` writes S out and `event_probability` counts; the first test checks that it reproduces these fractions exactly.

Two consequences follow directly. The **complement** rule: since A and "not A" are mutually exclusive and together make S, P(not A) = 1 − P(A). It is often the fastest route: P(at least one six in four rolls) = 1 − P(no six in four rolls) = 1 − (5/6)⁴ = 0.517747, which the 1296-outcome enumeration confirms. And the general **addition rule** (`@openstax-introductory-statistics-2e` §3.3): if A and B overlap, adding P(A) and P(B) counts the overlap twice, so

**P(A or B) = P(A) + P(B) − P(A and B)**,   (2)

which reduces to axiom (iii) when the overlap is empty. "Sum = 7" and "double" never happen together (7 is odd), so P(sum = 7 or double) = 6/36 + 6/36 − 0 = 1/3.

### 2.2 Counting: permutations, combinations, the binomial coefficient

For large sample spaces you count instead of listing (`@openstax-precalculus-2e` §11.5). The number of ordered selections of k items from n is P(n, k) = n!/(n − k)!; the number of unordered selections, the **combination** or **binomial coefficient**, is

**C(n, k) = n! / (k! (n − k)!)**.   (3)

Three dice out of ten can be chosen in C(10, 3) = 120 ways, arranged in P(10, 3) = 720. The code computes (3) by the multiplicative form C(n, k) = ∏ᵢ₌₁ᵏ (n − k + i)/i, in integers, so that C(52, 5) = 2 598 960 is exact and C(1000, 500) does not overflow a float. The name comes from the **binomial theorem** (§11.6): (a + b)ⁿ = Σₖ C(n, k) aᵏ bⁿ⁻ᵏ, which is the key to 2.5.

### 2.3 Conditional probability and independence

The probability of A *given that* B has happened is the fraction of B's outcomes that are also in A:

**P(A | B) = P(A and B) / P(B)**,   P(B) > 0.   (4)

Rearranged, (4) is the **multiplication rule**, P(A and B) = P(A | B) P(B). Two events are **independent** when knowing B tells you nothing about A, P(A | B) = P(A), which by (4) is the same as

**P(A and B) = P(A) P(B)**.   (5)

For independent components (5) multiplies: a flight sequence of five events (ignition, liftoff clear of the tower, stage separation, fairing release, payload release) with reliabilities 0.995, 0.998, 0.990, 0.997, 0.999 succeeds with probability 0.995 × 0.998 × 0.990 × 0.997 × 0.999 = 0.97915. This is the **series** formula: every link must hold. For a **parallel** arrangement, where any one of the components suffices, use the complement: all fail with probability ∏(1 − pᵢ), so P(at least one works) = 1 − ∏(1 − pᵢ). Three inertial measurement units of reliability 0.98 each give 1 − 0.02³ = 0.999992 if one working unit were enough.

Do not confuse "mutually exclusive" with "independent": mutually exclusive events with positive probability are *never* independent, because P(A and B) = 0 ≠ P(A)P(B) (`@openstax-introductory-statistics-2e` §3.2).

### 2.4 Total probability and Bayes' rule, through a tree

Let A₁, …, Aₘ be a **partition** of S (mutually exclusive, exhaustive). Any event B is the union of the pieces (B and Aᵢ), so by the axioms and the multiplication rule

**P(B) = Σᵢ P(B | Aᵢ) P(Aᵢ)**   (law of total probability).   (6)

Apply (4) in the other direction and substitute (6) for the denominator:

**P(Aᵢ | B) = P(B | Aᵢ) P(Aᵢ) / Σⱼ P(B | Aⱼ) P(Aⱼ)**   (Bayes' rule).   (7)

Introductory Statistics reaches (7) through tree diagrams (§3.5) without naming it; the 18.05 Reading 3 states it as above. The tree is the clearest way to compute it. An eddy-current inspection of fastener holes is applied to a batch in which 2 % of holes are cracked (the **prior**). The method detects a crack with probability 0.90 (its **sensitivity**) and flags an uncracked hole with probability 0.05 (its **false-alarm rate**):

```
                 ┌── flag   0.90  → cracked & flagged   0.02 × 0.90 = 0.0180
   cracked 0.02 ─┤
                 └── pass   0.10  → cracked & passed    0.02 × 0.10 = 0.0020
 ●
                 ┌── flag   0.05  → good & flagged      0.98 × 0.05 = 0.0490
   good   0.98 ──┤
                 └── pass   0.95  → good & passed       0.98 × 0.95 = 0.9310
```

The four leaves sum to 1 (a check the code makes). P(flag) = 0.0180 + 0.0490 = 0.0670 by (6), and the **posterior** probability that a flagged hole is cracked is 0.0180/0.0670 = 0.26866 by (7). Nearly three of four flagged holes are good: the false alarms on the 98 % of good holes outnumber the true detections on the 2 % of cracked ones. The posterior for a hole that passed is 0.0020/0.9330 = 0.00214, nine times lower than the prior. Both numbers drive real decisions: whether to scrap flagged parts or re-inspect them, and whether a passed part can be trusted. Lesson 20 turns this single update into a full Bayesian inference.

### 2.5 From coin flips to the binomial distribution

Flip a fair coin n times. The number of heads h has the sample space {0, …, n}, but its outcomes are *not* equally likely. Each particular sequence of n flips has probability (½)ⁿ by independence, and the number of sequences with exactly h heads is C(n, h) by 2.2, so by the addition rule over those sequences

**P(h heads) = C(n, h) / 2ⁿ**,   and for a biased coin   **P(h) = C(n, h) pʰ (1 − p)ⁿ⁻ʰ**.   (8)

The probabilities sum to 1 by the binomial theorem with a = p, b = 1 − p. Equation (8) is the **binomial distribution**; Lesson 16 studies it as a random variable with a mean and a variance, but its derivation is pure counting and belongs here.

(8) answers the voting question. A 2-of-3 IMU vote succeeds when at least two units work: P = C(3,2)p²(1 − p) + C(3,3)p³. At p = 0.98 that is 3 × 0.9604 × 0.02 + 0.941192 = 0.998816. Compare with a single unit (0.98) and with "all three must work" (0.941192, the series case): voting buys two extra nines of reliability from three mediocre units, which is why flight computers vote. In general the probability that at least k of n identical independent components work is Σⱼ₌ₖⁿ C(n, j) pʲ(1 − p)ⁿ⁻ʲ; k = 1 is the parallel formula and k = n the series one, which the tests confirm.

### 2.6 Monte Carlo estimation and its error

When S is too large to enumerate, or the event is defined by a simulation rather than a formula, **Monte Carlo** estimates P(A) by generating N independent outcomes and taking p̂ = (number in A)/N. The count is binomial with parameters N and p, and its variance is Np(1 − p) (Lesson 16 derives this; here take it as the fact that the error of a count of N fair-coin flips is about √N/2). Dividing by N,

**standard error of p̂: SE = √(p̂(1 − p̂)/N)**.   (9)

The error falls as 1/√N: ten times more precision costs a hundred times more trials. That is the price of Monte Carlo and the reason to count exactly whenever you can. The code seeds every generator (`numpy.random.default_rng(seed)`) so that a run is reproducible; a result that changes when the seed changes by more than about 2 SE is a bug, not randomness.

Running `monte_carlo_dice_sum(2, 7, N, seed=15)`:

| N | p̂ | SE from (9) | p̂ − 1/6 |
|---|---|---|---|
| 100 | 0.2100 | 0.0407 | +0.0433 |
| 10 000 | 0.1614 | 0.0037 | −0.0053 |
| 1 000 000 | 0.16735 | 0.00037 | +0.00068 |

Every error is within 1.5 SE of zero. The 1/√N law is tested properly in part 4.

### 2.7 Dice to dispersion: the coin-flip wind model, exactly

Now the payload. In each of n layers the east–west push is +d or −d with probability ½, so after n layers the east offset is x = d(2h − n), where h is the number of "east" layers, binomial by (8): P(x = d(2h − n)) = C(n, h)/2ⁿ. The north offset y has the same distribution and, because the two coins in each layer are independent, P(x, y) = P(x) P(y). The probability of landing inside a circle of radius R is then an exact finite sum over the (n + 1)² grid points:

**P(r ≤ R) = Σ_{x² + y² ≤ R²} P(x) P(y)**.   (10)

For n = 20 and d = 15 m that is 441 terms; `probability_within_radius_exact` evaluates them and gives P(r ≤ 100 m) = 0.6897. A seeded simulation of 20 000 descents gives 0.6875 ± 0.0033, a difference of 0.0022, two-thirds of one standard error. The spread of each axis is √n d = 67.08 m (the variance of a sum of independent ±d steps is n d², a fact Lesson 16 derives; the simulation returns 66.9 and 67.8 m).

Lesson 16 will show that for large n the walk's endpoint is nearly Gaussian and that the probability of landing within R is then 1 − exp(−R²/(2σ²)) with σ² = n d². At n = 20 that formula gives 0.6708 for R = 100 m, 2.7 % below the exact count, and the gap is *larger* than the Monte Carlo error: the grid is coarse, with landing points only every 30 m in each direction, so the circle of radius 100 m captures whole rows of grid points at once. Part 4 measures how the gap closes as n grows.

### 2.8 Where the model stops being valid

- **Equally likely outcomes.** Formula (1) is wrong for a loaded die or a prevailing wind; then the probabilities of the outcomes are inputs, not counts (extension 1).
- **Independence.** Real wind layers are correlated (a gust extends over hundreds of metres); real redundant components share power supplies and software. Correlation makes series systems *better* than (5) predicts and parallel systems *worse*; a common-cause failure can wipe out all three IMUs at once, and no vote helps.
- **Fixed drift per layer.** A real layer pushes by a continuous, random amount. The *shape* of the dispersion barely changes (extension 2) but the grid structure, and the exactness of (10), disappear.
- **Monte Carlo standard error.** (9) assumes N p̂ and N(1 − p̂) are both at least about 10; for a rare event with p = 10⁻⁴, N = 10⁴ trials give one hit or none and (9) is meaningless. Rare-event estimation needs importance sampling, beyond this lesson.

## 3. Build the computational model

Open `starter/main.py`. Every function is independent of the others except where noted, so build and test them in the order below.

### Step 1: counting and enumeration

```python
def n_choose_k(n, k):
    if n < 0: raise ValueError("n must be non-negative")
    if k < 0 or k > n: return 0
    k = min(k, n - k)
    result = 1
    for i in range(1, k + 1):
        result = result * (n - k + i) // i        # exact: each partial product is itself a binomial coefficient
    return result

def enumerate_dice(n_dice, sides=6):
    return np.array(list(itertools.product(range(1, sides + 1), repeat=n_dice)), dtype=int)

def event_probability(outcomes, predicate):
    mask = np.asarray(predicate(outcomes), dtype=bool)
    return float(np.count_nonzero(mask)) / outcomes.shape[0]
```

`dice_sum_pmf` counts the ways to reach each total by convolving the single-die count vector (six ones) with itself n − 1 times, then divides by 6ⁿ: the convolution is the addition rule over all the ways the last die can complete a total.

### Step 2: the rules, with their domain checks

`complement`, `addition_rule`, `conditional_probability` and `are_independent` are one-liners from (2), (4) and (5); each raises `ValueError` for an input outside [0, 1] or for P(A and B) > P(B), because a probability routine that silently accepts 1.3 will produce a plausible-looking wrong answer downstream. `total_probability(priors, likelihoods)` is the dot product in (6) after checking that the priors form a distribution; `bayes_posterior` divides by it, and `inspection_tree(prior, sensitivity, false_alarm)` returns the four leaves and the two posteriors as a dictionary.

### Step 3: reliability

```python
def binomial_probability(k, n, p):
    if k < 0 or k > n: return 0.0
    return n_choose_k(n, k) * p ** k * (1.0 - p) ** (n - k)        # (8)

def k_of_n_reliability(n, k, p):
    return float(sum(binomial_probability(j, n, p) for j in range(k, n + 1)))
```

`series_reliability` is `np.prod(ps)`; `parallel_reliability` is `1 - np.prod(1 - ps)`.

### Step 4: the Monte Carlo engine

```python
def monte_carlo_probability(sampler, event, n_trials, seed=0):
    rng = np.random.default_rng(seed)
    outcomes = sampler(rng, n_trials)                 # any array whose first axis indexes trials
    hits = np.asarray(event(outcomes), dtype=bool)
    p_hat = float(np.count_nonzero(hits)) / n_trials
    return p_hat, proportion_standard_error(p_hat, n_trials)
```

Separating the *sampler* (how an outcome is generated) from the *event* (what we are counting) lets the same engine estimate any probability. `monte_carlo_dice_sum` is the two-line specialisation; `monte_carlo_error_study` repeats it `n_repeats` times for each N with different seeds and returns the root-mean-square error against the exact pmf.

### Step 5: the coin-flip wind model

```python
def random_walk_landings(n_layers, drift, n_trials, seed=0):
    rng = np.random.default_rng(seed)
    signs = rng.choice([-1.0, 1.0], size=(n_trials, n_layers, 2))   # one coin per layer per axis
    return drift * signs.sum(axis=1)                                  # (n_trials, 2) offsets in metres

def walk_endpoint_pmf(n_layers, drift):
    h = np.arange(n_layers + 1)
    probs = np.array([n_choose_k(n_layers, int(k)) for k in h], float) / 2.0 ** n_layers
    return drift * (2 * h - n_layers), probs
```

`probability_within_radius_exact` forms the (n + 1) × (n + 1) grid of products P(x)P(y) and sums where x² + y² ≤ R², equation (10); use a relative tolerance of 10⁻¹² on the inequality so that a grid point exactly on the circle is counted. `probability_within_radius_mc` counts simulated landings inside the circle and returns (9). `radius_for_probability_exact` sorts the grid points by radius, accumulates their probabilities and returns the first radius at which the cumulative sum reaches the target. `gaussian_radius_probability` is the Lesson 16 formula, included only for comparison.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-15-probability-fundamentals`. The tests:

1. **Counting.** C(n, k) equals `scipy.special.comb(exact=True)` for 0 ≤ k ≤ n ≤ 30, obeys Pascal's rule C(n, k) = C(n−1, k−1) + C(n−1, k), and sums to 2ⁿ over k; the two-dice sample space has 36 outcomes and P(sum = 7) = 1/6 by enumeration; the dice-sum pmf sums to 1 and matches enumeration for three dice.
2. **The rules.** The addition rule reduces to P(A) + P(B) for disjoint events; conditional probability of an event given itself is 1; independence holds for the two dice faces and, surprisingly, for "sum = 7" against "first die = 1" (P(both) = 1/36 = 1/6 × 1/6, a classic trap), but fails for "sum = 6" against "first die = 1" (1/36 ≠ 5/36 × 1/6); Bayes' posterior equals the prior when the test is uninformative (sensitivity = false-alarm rate), approaches 1 as the false-alarm rate → 0, and the tree leaves sum to 1.
3. **Reliability.** Series ≤ every component, parallel ≥ every component; k-of-n at k = 1 equals parallel and at k = n equals series; it is monotone in p and in k; it reproduces the 2-of-3 value 0.998816 at p = 0.98 against `scipy.stats.binom.sf`.
4. **Monte Carlo.** With seed 15, N = 10⁵ trials estimate 1/6 within four standard errors (justified in the test's comment: a 4σ miss has probability 6 × 10⁻⁵ even for an unseeded run); the RMS error over N = 100…25 600 has log-log slope −0.5 ± 0.15; `proportion_standard_error` vanishes at p̂ = 0 and 1 and peaks at ½.
5. **Dispersion.** The walk endpoint pmf sums to 1, is symmetric, and matches `scipy.stats.binom`; simulated landings have mean within 3 SE of zero and variance within 5 % of n d²; the exact and simulated recovery probabilities agree within four standard errors; P(within R) is non-decreasing in R, 0 for R < d when n is odd, and 1 for R ≥ √2 n d; the Gaussian formula converges to the exact count as n grows.
6. `main()` writes `probability_fundamentals.png` to the directory it is given.

### Compare

`scipy.special.comb` and `scipy.stats.binom.pmf` / `.sf` are the oracles for counting; our integer loop matches `comb(exact=True)` to the last digit, and `binom.sf(k − 1, n, p)` reproduces `k_of_n_reliability` to 10⁻¹². We still wrote our own because the library hides the two ideas the lesson is about (that a combination counts sequences and that a reliability is a sum over an event) and because `scipy.stats` routines use logarithms and gamma functions internally that are faster and safer for n = 10⁶ but not more *correct* for n = 20. In practice you call the library and keep a counting check for one small case; that habit is what the tests encode.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All data are synthetic, generated by the coin-flip model with stated seeds.

1. **The 1/√N law.** Run `monte_carlo_error_study([100, 400, 1600, 6400, 25600], n_repeats=20, seed=15)`. Before running, predict the RMS error at each N from (9) with p = 1/6: 0.0373, 0.0186, 0.0093, 0.0047, 0.0023. Measured: 0.0372, 0.0212, 0.0111, 0.0038, 0.0019; fitted log-log slope −0.551. Record both and the slope.
2. **Grid against bell.** For n = 4, 10, 20, 50 layers compute the exact probability of landing within 1.5σ (σ = √n d) and the Gaussian formula's value 1 − e^{−1.125} = 0.6753. Predict first whether the exact value is above or below, and whether the gap closes. Measured: 0.7656, 0.7383, 0.6897, 0.6764 against 0.6753: always above, and the gap falls from 0.090 to 0.001.
3. **Recovery radius.** For n = 20, d = 15 m, find the radius containing 50 %, 90 % and 99 % of landings exactly (`radius_for_probability_exact`) and by the Gaussian formula σ√(−2 ln(1 − P)). Measured: 84.9, 150.0, 201.2 m exact against 79.0, 144.0, 203.6 m Gaussian.
4. **Seed sensitivity.** Repeat the 20 000-descent simulation with seeds 1–5 and record P(r ≤ 100 m); every value should lie within about 2 × 0.0033 of 0.6897.

Raw data format: three CSV tables (`N, rms_measured, rms_predicted`; `n_layers, p_exact, p_gauss`; `prob, R_exact, R_gauss`) and one list of five seeded estimates.

What you should find: the error law holds in slope, with the measured points scattering around the prediction by the amount expected from only 20 repeats; the exact count sits *above* the Gaussian value because the binomial endpoint has lighter tails than a Gaussian for small n (no outcome can exceed n d, whereas a Gaussian extends to infinity), and the gap at n = 20 is already below the 1 % level for the 99 % radius, where the decision is made.

## 5. Validate and improve

| quantity | exact counting (§2) | simulation / library | difference | cause |
|---|---|---|---|---|
| P(sum = 7), two dice | 6/36 = 0.166667 | 0.16735 ± 0.00037 (N = 10⁶) | 1.8 SE | sampling |
| P(≥ 1 six in 4 rolls) | 1 − (5/6)⁴ = 0.517747 | 0.517747 (enumeration of 1296) | 0 | both exact |
| 2-of-3 vote, p = 0.98 | 0.998816 | 0.998816 (`binom.sf`) | < 10⁻¹² | both exact |
| P(cracked ∣ flagged) | 0.26866 | 0.26866 (`bayes_posterior`) | 0 | both exact |
| P(r ≤ 100 m), n = 20 | 0.6897 | 0.6875 ± 0.0033 (N = 20 000) | 0.7 SE | sampling |
| same, Gaussian formula | 0.6897 | 0.6708 | 0.019 = 5.7 SE | discrete grid, light tails |
| per-axis σ | 67.08 m | 66.9, 67.8 m | 0.4 % | sampling (SE of an SD ≈ σ/√(2N) = 0.3 m) |

The only differences larger than the standard error come from the Gaussian approximation, and part 4 shows they vanish as n grows. That is the right kind of discrepancy: a model simplification whose size we can compute, not an unexplained residual.

What the agreement does not show: that real winds behave like fair coins. The model has no prevailing wind (which would shift the centre of the scatter), no correlation between layers (which would widen it: perfectly correlated layers give a spread of n d = 300 m instead of √n d = 67 m), and no variation in the size of the push. Each of those is a one-line change to `random_walk_landings`; the exact formula (10) survives only the first. A real dispersion analysis (Lesson 16 onward) replaces the coins by measured wind statistics and the fixed drift by the payload's descent dynamics, and then only the Monte Carlo route is available, which is why its error law matters.

Improvement worth making now: add a biased coin (P(east) = 0.6) and check that the mean offset becomes n d(2p − 1) = 60 m east and the variance 4n d² p(1 − p) = 4320 m² (σ = 65.7 m), both of which (8) gives directly. Improvement worth *not* making: refining the grid by using 200 layers of 1.5 m each; the exact sum (10) still works but the Gaussian formula is then within 0.1 %, and you would be doing counting to confirm what the central limit theorem already guarantees.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (6) and (7) and the inspection tree with its four leaves; `probability_fundamentals.png`; the four tables from part 4; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that recommends a recovery radius for the 90 % and 99 % cases, states how it was computed and how large the simulation's standard error is, and names the assumption you trust least.

### Adversarial self-checks

- **What if the number of layers doubles?** The per-axis spread grows by √2 (to 94.9 m for n = 40), not by 2; the recovery radius for a given probability grows by the same factor. If instead the drift per layer doubles, the spread doubles. Check both with `walk_endpoint_pmf`.
- **What if the per-component reliability goes to 1? To 0?** Series, parallel and k-of-n all go to 1 and to 0 respectively; between, 2-of-3 beats a single unit only when p > ½ (solve 3p² − 2p³ > p). Below p = ½ voting makes things *worse*: three unreliable units outvote the one that is right.
- **Which assumption fails first in the real world?** Independence between wind layers. A gust that spans several layers makes the pushes correlated, and the spread then lies between √n d and n d. The second to fail is the fixed push size.
- **What measurement could falsify the model?** The landing scatter of twenty real descents under similar soundings: if its per-axis standard deviation is far from √n d, or its centre is far from the aim point, the coins are biased or correlated. (Twenty descents give an SD with about 16 % uncertainty, so "far" means more than a factor 1.3.)
- **Does the result survive a change in seed, N, or method?** The exact count does not depend on any of them. The simulation moves by about one SE per seed change and by 1/√N with N; if it moves by more, the generator is being re-seeded inside a loop or the trials are not independent.
- **Which decision would change if the uncertainty doubled?** With SE = 0.0066 instead of 0.0033 on P(r ≤ 100 m), the 95 % interval would span 0.674–0.701 and the Gaussian value 0.6708 would fall inside it: you could no longer tell the two models apart from this simulation, and would have to either run four times more trials or trust the exact count. The recovery-radius recommendation (150 m for 90 %) would not change, because it comes from the exact sum.

## Check yourself

1. State the three axioms of probability and derive the complement rule from them in one line.
2. Two events have P(A) = 0.3, P(B) = 0.5 and P(A and B) = 0.15. Are they independent? Are they mutually exclusive? What is P(A or B)?
3. Why can two mutually exclusive events with positive probability never be independent?
4. Write the law of total probability for a partition {A₁, A₂} and use it to compute P(flag) in the inspection tree with prior 0.02, sensitivity 0.90 and false-alarm rate 0.05.
5. A test with sensitivity 0.99 and false-alarm rate 0.01 is applied where the prior defect rate is 10⁻⁴. Roughly what fraction of flagged parts are actually defective, and why is it so low?
6. How many Monte Carlo trials are needed to estimate a probability near 0.5 to a standard error of 0.001? Of 0.0001?
7. After n fair-coin layers of drift d, what are the possible landing offsets along one axis, and what is the probability of the largest one?
8. For the 20-layer model, is the exact probability of landing within 1.5σ above or below the Gaussian value, and what feature of the binomial distribution explains the sign?
