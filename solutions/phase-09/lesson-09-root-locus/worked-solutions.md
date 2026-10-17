# P09.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. Your altitude loop (s = 15: T = 1.1 s)

**Assumptions.** G = 50/(s(1.1s + 1)(0.5s + 1)), proportional control, K > 0.

(a) Poles **0, −0.9091, −2 s⁻¹**. Segments **[−0.909, 0]** and **(−∞, −2]**. Asymptotes **60°, 180°, 300°** from **σ_a = (0 − 0.9091 − 2)/3 = −0.9697**.

(b) N D′ − N′D ∝ 1.65s² + 3.2s + 1 = 0: s = −0.3916 (on [−0.909, 0]) or −1.548 (rejected). **Breakaway −0.3916, K = 0.003585 rad/m.**

(c) den = 0.55s³ + 1.6s² + s, so the characteristic polynomial is 0.55s³ + 1.6s² + s + 50K. Routh's s¹ entry (1.6 − 0.55·50K)/1.6 vanishes at **K = 1.6/(0.55 × 50) = 0.05818 rad/m**; the auxiliary 1.6s² + 50K = 0 gives **ω = √(1/0.55) = 1.3484 rad/s**. `jw_crossings` returns (0.058182, 1.348400).

(d) **K = 0.009814 rad/m**, pair **−0.3125 ± 0.5413j** (ω_n = 0.625 rad/s), third pole **−2.284**. Predicted 16.3 % and 4/(0.5 × 0.625) = 12.8 s; simulated **15.6 % and 13.3 s**. The pair dominates (the third pole is 7.3 times further left), so the prediction is close.

**Check.** The root sum −0.3125·2 − 2.284 = −2.909 = −1.6/0.55.

**Common mistakes.** Putting the centroid at the mean of the poles' *magnitudes*; keeping the breakaway candidate on the wrong segment.

## 2. The rules

(a) 1 + KG(s) = 0 ⇔ G(s) = −1/K. For real K > 0, −1/K is a negative real number: angle 180° (mod 360°), magnitude 1/K. Conversely, if ∠G(s) = 180°, G(s) is a negative real number −g, and K = 1/g > 0 makes s a root. So the angle decides membership; the magnitude gives **K = 1/\|G(s)\|**.

(b) For real s, a real pole or zero r contributes ∠(s − r) = 0 if r < s and 180° if r > s; a complex pair contributes θ + (−θ) = 0. The total angle is 180° × (number of real poles and zeros to the right of s), up to sign. It is an odd multiple of 180° exactly when that number is odd.

(c) For large s, D(s)/N(s) = s^{n−m}(1 − Σp/s + …)/(1 − Σz/s + …) = s^{n−m}(1 − (Σp − Σz)/s + O(s⁻²)) = (s − σ_a)^{n−m}(1 + O(s⁻²)), with σ_a = (Σp − Σz)/(n − m). On the locus D/N = −K, so (s − σ_a)^{n−m} ≈ −K: s − σ_a = K^{1/(n−m)} e^{j(2k+1)π/(n−m)}, rays from σ_a at **(2k + 1)180°/(n − m)**.

(d) On the real axis, K(σ) = −D(σ)/N(σ). Two branches meet where K(σ) has a local extremum: **dK/dσ = −(D′N − DN′)/N² = 0 ⇔ N D′ − N′D = 0**. Root sum: D + KN = a_n sⁿ + a_{n−1}s^{n−1} + …; if m ≤ n − 2, KN does not reach s^{n−1}, so −a_{n−1}/a_n, the sum of the roots, does not change with K.

**Common mistakes.** Forgetting that candidates from N D′ − N′D = 0 must also have K > 0.

## 3. An unstable airframe

**Assumptions.** L = K(s + 1)/((s² − 4)(0.1s + 1)); poles +2, −2, −10; zero −1.

(a) Segments **[−10, −2]** and **[−1, +2]**. Two asymptotes at **±90°** from σ_a = (2 − 2 − 10 − (−1))/2 = **−4.5**. Breakaway on [−10, −2] at **−5.248**. The branch from +2 runs left along the axis to the zero at −1. The branches from −2 and −10 meet at −5.248 and leave along the ±90° asymptotes.

(b) Characteristic polynomial 0.1s³ + s² + (K − 0.4)s + (K − 4). Routh's first column is 0.1, 1, (K − 0.4) − 0.1(K − 4) = 0.9K, K − 4. All positive ⇔ **K > 4**. The crossing is at **s = 0, K = 4**: the branch from +2 passes through the origin (`jw_crossings` returns (4.0, 0.0)). The roots confirm it: one unstable root at K = 2 (+0.809), a root at 0 at K = 4, all stable at 4.5, 10, 50, 200.

(c) **K = 9.742**, pair **−4.671 ± 8.090j**, third pole **−0.658**. The pair is *not* dominant: the real pole at −0.658 is seven times slower and only partly cancelled by the zero at −1. The step response will be a slow creep set by that pole.

(d) The airframe diverges on its own (the pole at +2 has a 0.35 s time-to-double); a weak controller cannot catch it. Stability needs a minimum gain. Losing gain (a weak actuator, a sensor fault, saturation) is as dangerous as having too much. `P09.10` shows the same fact as a counterclockwise encirclement in the Nyquist plot.

**Check.** At K = 4 the constant term of the characteristic polynomial vanishes: s = 0 is a root.

**Common mistakes.** Assuming "more gain is less stable" for every plant; missing the crossing at ω = 0.

## 4. The climb-rate feedback (s = 15: T_d = 2.0 s)

(a)

| T_d | K | pair (ω_n) | third pole | predicted OS, t_s | simulated OS, t_s |
|---|---|---|---|---|---|
| 2.0 s | 0.02851 rad/m | −1.284 ± 2.224j (2.568) | −0.433 | 16.3 %, 3.12 s | 5.4 %, 4.40 s |
| 1.0 s | 0.04000 | −1.000 ± 1.732j (2.000) | −1.000 | 16.3 %, 4.00 s | 16.3 %, 4.04 s |
| 3.0 s | 0.02056 | −1.361 ± 2.358j (2.723) | −0.277 | 16.3 %, 2.94 s | 0.6 %, 8.03 s |

(b) T_d = 1 s: the zero at −1 cancels the open-loop pole at −1, the closed loop is exactly second order, and the prediction is exact (the 0.04 s difference is the 2 % band crossing on a 0.01 s grid). T_d = 2 s: the third pole at −0.433 sits near the zero at −0.5. The residue of a closed-loop pole near a zero is small, so it adds a slow tail that suppresses the overshoot (5.4 %) and stretches the settling (4.4 s). T_d = 3 s: the pole at −0.277 is even slower and its tail, though small, lasts long enough to set the 2 % settling at 8.0 s.

(c) **T_d = 1 s** settles fastest in the simulation (4.04 s), although T_d = 3 s has the largest ω_n. The slowest pole, not the fastest pair, sets the settling time.

**Common mistakes.** Reporting the predicted metrics without simulating; ignoring the zero's own effect (a zero adds overshoot, a near-cancelled pole adds a tail).

## 5. Routh's special cases

(a) s⁴ + 2s³ + 3s² + 4s + 5: first column **1, 2, 1, −6, 5**: two sign changes, **2 RHP roots** (0.288 ± 1.416j). Regular.

(b) s³ + 3s + 2: rows [1, 3], [0, 2]. **Lone zero in the first column**: replace by ε; s¹: 3 − 2/ε → −∞; s⁰: 2. Column 1, ε, −∞, 2: **2 RHP roots** (0.298 ± 1.807j).

(c) s⁵ + s⁴ + 4s³ + 4s² + 3s + 3: the s³ row is **all zeros**; auxiliary from s⁴: s⁴ + 4s² + 3, derivative 4s³ + 8s; continue: column **1, 1, 4, 2, 2, 3**: no sign change, **0 RHP roots**, but the auxiliary (s² + 1)(s² + 3) puts **four roots on the imaginary axis** (±j, ±j√3). Marginally stable, not stable.

(d) 2s⁴ + s³ + 2s² + s + 1: the s² row starts with (1·2 − 2·1)/1 = 0, rest nonzero: **ε case**; s¹: (ε − 1)/ε → −1/ε; s⁰: 1. Column 2, 1, ε, −∞, 1: **2 RHP roots** (0.226 ± 0.875j).

(e) s³ + 3s² + 3s + 1 = (s + 1)³: column **1, 3, 8/3, 1**: **0 RHP roots**. Regular despite the triple root.

All five agree with `routh_rhp_count` and `np.roots`.

**Common mistakes.** Reading "no sign changes" in (c) as "stable"; for ε, taking ε < 0 (the count is the same either way only if done consistently).

## 6. Is ζ = 0.5 enough?

(a) With K = 0.01037 rad/m fixed:

| change | ζ of the pair | critical gain |
|---|---|---|
| τ = 0.4 s | 0.536 | 0.0700 rad/m |
| τ = 0.5 s | 0.500 | 0.0600 |
| τ = 0.7 s | 0.438 | 0.0486 |
| V = 40 m/s | 0.595 | 0.0750 |
| V = 50 m/s | 0.500 | 0.0600 |
| V = 60 m/s | 0.429 | 0.0500 |

A 20 % change in airspeed moves ζ about as much as a 40 % change in the inner-loop lag. Airspeed is known well and can be scheduled; the inner-loop lag is the less certain quantity, so both deserve margin.

(b) The climb-rate term K T_d ḣ passes sensor noise straight to the pitch command, amplified by K T_d = 0.04 rad/(m/s) for T_d = 1 s; in turbulence the vertical-speed signal is noisy and a filter (with its own lag) is needed. The "exact" second-order result rests on the zero cancelling the flight-path pole at −1/T exactly; T varies with airspeed and weight, so the cancellation is never exact in flight.

(c) **No.** 10 % overshoot needs ζ ≥ 0.591, and 2 % settling within 5 s needs ζω_n ≥ 4/5 = 0.8 s⁻¹. On the P-control locus the complex pair leaves the axis at the breakaway −0.4226 and moves right as K grows; the largest decay rate ζω_n it ever has is **0.4226 s⁻¹**, so the settling time can never be shorter than **9.5 s**. The requirement needs a zero (rate feedback or lead) to pull the branches left.

**Common mistakes.** Trying gains until the numbers look right; the locus shows at once that no gain works.

## Check-yourself answers

1. ∠G(s) = 180° (mod 360°) and K\|G(s)\| = 1. The angle condition decides membership; the magnitude condition gives the gain.
2. [−1, 0] and (−∞, −2].
3. 60°, 180°, 300° from σ_a = (0 − 1 − 2)/3 = −1.
4. At −1 + 1/√3 = −0.4226; the other root, −1.5774, lies between −2 and −1, which is not on the locus (it would need K < 0).
5. Column 1, 3, (6 − K)/3, K: stable for **0 < K < 6**; at K = 6 the auxiliary 3s² + 6 gives **ω = √2 rad/s**.
6. The coefficient of s^{n−1} in D + KN comes only from D when deg N ≤ n − 2, and minus that coefficient over the leading one is the sum of the roots.
7. A slow real closed-loop pole near a zero (T_d = 3 or 4 s here) adds a slow tail that removes the overshoot and dominates the settling.
