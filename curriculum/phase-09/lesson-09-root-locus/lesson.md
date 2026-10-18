# P09.09  Root Locus

> As the gain rises from zero, the closed-loop poles leave the open-loop poles and travel to the zeros or to infinity along paths fixed by one rule: the loop gain must equal −1. Draw the paths and you can read the gain that gives the damping you want.

| | |
|---|---|
| **Learning outcomes** | 1. State the angle condition ∠G(s) = 180° (mod 360°) and the magnitude condition K\|G(s)\| = 1, and verify that every computed closed-loop root satisfies both. 2. Derive and apply the construction rules: real-axis segments, the number of branches, asymptote angles and centroid, breakaway points from dK/ds = 0, and the conservation of the root sum. 3. Build the Routh array (including the zero-in-first-column and row-of-zeros cases), count right-half-plane roots from it, and find the gain and frequency at which the locus crosses the imaginary axis, checked two ways. 4. Choose the altitude-hold gain that places the dominant pair at ζ = 0.5, judge from the third pole whether the dominant-pair prediction can be trusted, and verify it by simulating the step response. |
| **Prerequisites** | `P09.03` (transfer functions, feedback), `P09.04` (poles and time response, ζ and ω_n, the effect of a zero), `P09.08` (closed-loop characteristic polynomial, Routh for a cubic), `P02.07` (complex numbers in polar form), `P01.16` (polynomial roots as eigenvalues). Skill check: you can add the angles of the vectors from three points to a fourth point in the complex plane. |
| **Engineering question** | An altitude-hold loop on a small aircraft or UAV commands pitch from altitude error through a closed inner pitch loop (τ = 0.5 s) and a flight-path lag (T = 1.0 s) at V = 50 m/s. What proportional gain (rad of pitch per metre of altitude error) gives a dominant closed-loop pair with ζ = 0.5, how close is that gain to the one that makes the loop oscillate, and does adding climb-rate feedback (a zero) help? Answer: K = 0.01037 rad/m (0.594°/m), 17 % of the 0.0600 rad/m at which the locus crosses the imaginary axis at 1.414 rad/s; with climb-rate feedback the loop never crosses, but whether the ζ = 0.5 pair still dominates depends on where the zero sits. |
| **Assumptions and limits** | A linear, low-order model: the closed inner pitch loop is approximated by a first-order lag, the flight-path response by another, altitude by V∫γ dt, with constant airspeed. Gains are positive (the K < 0 "complementary" locus is not drawn). The ζ–ω_n formulas for overshoot and settling assume two dominant poles and no nearby zero; the lesson tests that assumption by simulation. Pitch limits, sensor noise on the climb-rate signal and the inner loop's own dynamics beyond one lag are ignored. The numbers are representative, not one vehicle's data. |
| **Analytical work** | The characteristic equation 1 + KG(s) = 0 and its two conditions; real-axis segments by counting poles and zeros to the right; n − m asymptotes at (2k + 1)180°/(n − m) from the centroid (Σp − Σz)/(n − m); breakaway from d/ds[D/N] = 0; constancy of the root sum when n − m ≥ 2; the Routh array and its special cases; the imaginary-axis crossing from Im[D(jω)N(−jω)] = 0 and from the auxiliary polynomial; the ζ line s = ω_n(−ζ + j√(1 − ζ²)). |
| **Computational deliverable** | `starter/main.py`: `char_poly`, `locus_roots`, `angle_condition_error`, `gain_at_point`, `real_axis_segments`, `asymptotes`, `breakaway_points`, `routh_array`, `routh_rhp_count`, `jw_crossings`, `gain_for_damping` (`poly_add`, `tf_eval`, `altitude_plant`, `match_roots`, `tf_to_companion`, `step_response`, `overshoot_and_settling`, `pd_design_table` and `main` are given). Tests prove: every locus point of three plants satisfies the angle condition to 10⁻⁶ degrees and K\|G\| = 1 to 10⁻⁶; points off the locus fail it; branches start at the open-loop poles, match `np.roots`, and never jump; the root sum is constant; far branches follow the asymptotes; the real-axis rule agrees with the angle condition at 141 points; the breakaway is a double root and the off-locus candidate is rejected; the Routh count equals the `np.roots` count for six random quartics and three special cases; the jω crossing gain and frequency agree with Routh and with the auxiliary polynomial; rate feedback removes the crossing; the damping-target gain lands a root on the ζ line and equals the closed form 0.04 rad/m when a zero cancels a pole; the step response matches `scipy.signal.step`. |
| **Practical activity** | **computational experiment**: sweep the climb-rate time constant T_d of a controller K(T_d s + 1), find the gain for ζ = 0.5 at each, and compare the predicted overshoot and settling with simulated step responses; and test the sensitivity of the design to the inner-loop lag and the airspeed. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `root_locus.png`, a hand sketch of the locus with segments, asymptotes, breakaway and jω crossing marked before running the code, the Routh table at the critical gain, the T_d experiment table, and a half-page recommendation of a gain (and T_d if used) with its predicted and simulated step metrics. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the hand sketch has the right segments, asymptotes and crossing; all public tests pass; the recommendation explains with numbers when the dominant-pair prediction fails. |
| **References** | `@mit-ocw-16-06-principles-of-automatic-control` Lecture 8 "The Routh criterion" (The Routh Stability Criterion; The Routh Array), Lecture 11 "The root locus method", Lecture 12 "Root locus rules" and Lecture 13 "Root locus rules, lead compensation" (PDF "Root Locus Rules (cont-d)", with a Routh array for the jω crossing) (required-free): the angle and magnitude conditions and the construction rules (real-axis segments, asymptotes and centroid, breakaway, jω crossing by Routh); Lecture 15 "Zero degree root locus" (required-free) for the K < 0 extension. `@mit-ocw-16-30-feedback-control-systems` Topic 2 "Basic root locus: analysis and examples": Root Locus Basics; Root Locus Analysis; Root Locus Examples; Dynamic Compensation; Classic Root Locus Approach (required-free): the same rules from an aerospace course, and picking K for a damping-ratio target. No verified free source treats an outer altitude loop around an inner pitch loop by root locus: the altitude-hold loop is original instruction in this lesson. |
| **Extensions** | ★ Draw the complementary locus (K < 0) and explain what it would mean physically for the altitude loop. ★ Add a second-order inner pitch loop (ω_n = 3 rad/s, ζ = 0.7) instead of the first-order lag and find how much the critical gain changes. ★ Place a lead compensator (zero at −0.7, pole at −5) and design for ζ = 0.5 with the fastest ω_n you can get; compare with `P09.12`. |
| **Time** | unmeasured |

## 1. The engineering problem

An autopilot holds altitude by commanding pitch: below the target, nose up; above it, nose down. The obvious law is proportional, θ_c = K(h_c − h). The questions are how large K should be and what happens as it changes. Too small and the aircraft drifts back to altitude over tens of seconds; too large and it porpoises, overshooting the target, and at some gain the oscillation stops decaying. Between them lies a gain that gives a well-damped return. The designer wants that gain, its distance from the oscillating one, and a picture of how the closed-loop behaviour changes in between.

The **root locus** is that picture. The closed-loop poles are the roots of the characteristic polynomial D(s) + KN(s), and as K runs from 0 to ∞ they move continuously in the complex plane. Evans's insight was that their paths obey simple rules that can be sketched by hand from the open-loop poles and zeros, and that every point on them satisfies one angle condition. The sketch shows at a glance whether a design is possible, and the angle and magnitude conditions give the exact gain. This lesson builds the locus by computing polynomial roots, checks every computed point against the two conditions, and derives the construction rules, including the Routh test for where the locus crosses into the right half plane, as independent checks on the computation.

Success criteria, written first: every computed locus point satisfies ∠G = 180° within 10⁻⁶ degrees and K\|G\| = 1 within 10⁻⁶; the imaginary-axis crossing from the Im/Re conditions agrees with the Routh array (the right-half-plane count jumps from 0 to 2 there) and with its auxiliary polynomial to 10⁻⁸; the gain for ζ = 0.5 places a closed-loop root on the ζ line to 10⁻⁶.

## 2. Understand and derive

### 2.1 The plant

Model the closed inner pitch loop as θ/θ_c = 1/(τs + 1), the flight path as γ/θ = 1/(Ts + 1) (the aircraft's path follows its attitude with a lag set by how quickly lift builds), and altitude as ḣ = Vγ. Then

**G(s) = h/θ_c = V / (s(Ts + 1)(τs + 1)),  V = 50 m/s, T = 1.0 s, τ = 0.5 s,**   (1)

with poles at 0, −1 and −2 s⁻¹, and no finite zeros. The loop gain is KG with K in rad/m; with proportional control the characteristic polynomial is

**D(s) + KN(s) = 0.5s³ + 1.5s² + s + 50K.**   (2)

### 2.2 The two conditions

A point s is a closed-loop pole for some K > 0 if and only if 1 + KG(s) = 0, i.e. KG(s) = −1. A complex number equals −1 when its angle is 180° (mod 360°) and its magnitude is 1:

**∠G(s) = Σ∠(s − z_i) − Σ∠(s − p_j) = 180° (mod 360°),  K = 1/\|G(s)\|.**   (3)

The angle condition alone (16.06 Lecture 11) decides whether s is on the locus; the magnitude condition then gives the gain that puts a pole there. Every point computed in this lesson is checked against both: of 1800 points (600 gains, three branches), the worst angle error is 1.4 × 10⁻¹³ degrees and the worst \|K\|G\| − 1\| is 3.2 × 10⁻¹².

### 2.3 Construction rules

All follow from (3) and from the polynomial D + KN (16.06 Lecture 12; 16.30 Topic 2, "Root Locus Analysis").

1. **Branches.** There are n = deg D branches. At K = 0 they start at the open-loop poles; as K → ∞, m of them end at the zeros (D/K + N → N) and n − m go to infinity.
2. **Real axis.** For a real s, each complex-conjugate pair contributes zero net angle, and each real pole or zero to the left contributes 0°, to the right 180°. So **a real point is on the locus when the number of real poles and zeros to its right is odd.** For (1): (−∞, −2] and [−1, 0].
3. **Asymptotes.** For large \|s\|, G ≈ s^{−(n−m)} and (3) gives n − m directions at **(2k + 1)180°/(n − m)**, here 60°, 180° and 300°. They radiate from the **centroid** σ_a = (Σp − Σz)/(n − m) = (0 − 1 − 2)/3 = **−1**. Its origin is rule 4.
4. **Root sum.** If n − m ≥ 2, the coefficient of s^{n−1} in D + KN does not depend on K, so the sum of the closed-loop poles is constant: here −a₂/a₃ = −3 for every K. The code confirms it to 10⁻⁹ over the whole sweep.
5. **Breakaway.** On a real segment between two poles, K = −D(s)/N(s) rises from 0, reaches a maximum and the two branches leave the axis. At the maximum dK/ds = 0, that is **N D′ − N′D = 0**. For (1) this is 1.5s² + 3s + 1 = 0, roots −0.4226 and −1.5774. Only −0.4226 = −1 + 1/√3 lies on a segment (K = 0.003849 rad/m there, a double root); −1.5774 would need K < 0 and is rejected.
6. **Imaginary-axis crossing.** Where a branch crosses, s = jω and D(jω) + KN(jω) = 0 with K real. Equivalently Im[D(jω)·conj(N(jω))] = 0, a real polynomial in ω, then K = −D(jω)/N(jω). For (2): Im part −0.5ω³ + ω = 0 gives ω = √2, and the real part −1.5ω² + 50K = 0 gives **K = 0.0600 rad/m at ω = 1.4142 rad/s**.

### 2.4 Routh: the independent check on the crossing

The **Routh array** (16.06 Lecture 8) decides how many roots of a real polynomial lie in the right half plane without computing them. Write the even-indexed coefficients in the first row and the odd-indexed ones in the second; each further entry is R[i, j] = (R[i−1, 0]R[i−2, j+1] − R[i−2, 0]R[i−1, j+1])/R[i−1, 0]. The number of sign changes down the first column equals the number of right-half-plane roots. For (2):

| row | first column |
|---|---|
| s³ | 0.5 |
| s² | 1.5 |
| s¹ | (1.5 · 1 − 0.5 · 50K)/1.5 |
| s⁰ | 50K |

The s¹ entry vanishes at K = 1.5/(0.5 · 50) = 0.06: below it the column is all positive (stable), above it two sign changes (two unstable roots). At K = 0.06 the s² row gives the **auxiliary polynomial** 1.5s² + 3 = 0, whose roots ±j√2 are the crossing. `main()` prints the first column at 0.5, 0.999, 1.001 and 2 times the critical gain: [0.5, 1.5, 0.5, 1.5], [0.5, 1.5, 0.001, 2.997], [0.5, 1.5, −0.001, 3.003] and [0.5, 1.5, −1, 6], with 0, 0, 2 and 2 right-half-plane roots, the same as `np.roots`.

Two special cases need care. A **zero in the first column** with the rest of the row nonzero: replace it by a small ε > 0 and continue; the sign changes are counted as ε → 0⁺ (s⁴ + s³ + 2s² + 2s + 3 has two unstable roots). A **whole row of zeros** means the polynomial has roots placed symmetrically about the origin (on the imaginary axis, or ±a real pairs). Replace the row by the coefficients of the derivative of the auxiliary polynomial from the row above. For (s² − 4)(s + 1) = s³ + s² − 4s − 4 the s¹ row vanishes; the auxiliary s² − 4 has derivative 2s, the column becomes 1, 1, 2, −4, and the single sign change is the root at s = +2.

### 2.5 Choosing the gain for a damping ratio

For a pair of poles −ζω_n ± jω_n√(1 − ζ²), the ratio ζ fixes the angle of the pole from the negative real axis, arccos ζ = 60° for ζ = 0.5 (`P09.04`). Searching along that ray for the angle condition (3) gives the point where the locus crosses it, and the magnitude condition gives the gain. For (1): **s = −0.3333 ± 0.5774j (ω_n = 0.667 rad/s), K = 0.01037 rad/m = 0.594 deg/m**, with the third pole at −2.333 (the root sum −3 again). The gain is 17 % of the critical 0.06, a gain margin of 5.8 (`P09.10` gives the frequency-domain view).

The **dominant-pair prediction** is overshoot e^{−πζ/√(1−ζ²)} = 16.3 % and 2 % settling time 4/(ζω_n) = 12.0 s. It assumes the third pole is far to the left and there is no zero nearby. Here the third pole is 7 times further left than the pair's real part, and the simulated step gives **15.5 % and 12.5 s**: close.

### 2.6 A zero: climb-rate feedback

Add climb-rate feedback, θ_c = K(h_c − h − T_d ḣ), equivalent to a controller K(T_d s + 1) with a zero at −1/T_d. The pole excess falls to 2, the asymptotes become ±90° from the centroid (−3 + 1/T_d)/2, and for (1) the characteristic polynomial 0.5s³ + 1.5s² + (1 + 50KT_d)s + 50K has a positive Routh column for every K > 0 when T_d > 1/3 s: **the locus never crosses the imaginary axis**. The zero bends the branches left. But a zero also creates a slow closed-loop pole between itself and the origin, and that pole, not the ζ = 0.5 pair, may dominate. That is the experiment of part 4.

### 2.7 Where the model stops

The first-order inner loop is the weakest piece: a real pitch loop has its own second-order dynamics and actuator lag, which add phase and lower the critical gain (extension 2). Pitch is limited (±10 to 15° for passenger comfort and stall margin), so a large altitude error saturates the command and the linear analysis no longer applies; altitude-hold modes typically switch to a climb-rate mode for large errors for this reason. Airspeed enters (1) linearly: at 60 m/s the critical gain falls to 0.050 rad/m and the ζ of a fixed-gain design to 0.43 (Problem 6), so a fixed gain must be scheduled with speed.

## 3. Build the computational model

Open `starter/main.py`. The plant, the root matcher, the step simulator and `main()` are given.

### Step 1: the locus

```python
def char_poly(num, den, K):
    return poly_add(den, K * np.asarray(num, dtype=float))

def locus_roots(num, den, gains):
    out, prev = np.empty((len(gains), len(den) - 1), complex), None
    for i, K in enumerate(gains):
        r = np.roots(char_poly(num, den, K)).astype(complex)
        r = np.sort_complex(r) if prev is None else match_roots(prev, r)
        out[i], prev = r, r
    return out
```

`match_roots` (given) reorders each new set of roots to be closest to the previous one, so each column is a continuous branch; without it `np.roots` returns roots in an arbitrary order and a plotted branch would jump between paths.

### Step 2: the conditions and the rules

`angle_condition_error` returns ∠G(s) − 180°, wrapped to (−180°, 180°]. `gain_at_point` is 1/\|G(s)\|. `real_axis_segments` sorts the real poles and zeros in descending order and pairs them (first with second, third with fourth, …, the last with −∞ if the count is odd). `asymptotes` returns the centroid and the angles. `breakaway_points` solves N D′ − N′D = 0 with `np.polyder` and keeps the real roots where −D/N > 0.

### Step 3: Routh and the crossing

`routh_array` fills the table row by row, handling the two special cases of 2.4. `routh_rhp_count` counts sign changes in its first column. `jw_crossings` writes D(jω) and N(jω) as real polynomials in ω (the coefficient of s^k contributes to the real or imaginary part according to j^k), solves Im[D conj(N)] = 0, and keeps the roots with ω ≥ 0 at which K = −D/N is real and positive.

### Step 4: the damping target

`gain_for_damping` scans ω_n along the ζ ray, finds a sign change of the angle error that is not a ±180° wrap (a jump of less than 90°), refines it by bisection and returns (K, s).

Running `main()` prints the plant, the segments and asymptotes, the breakaway, the crossing and the Routh table of 2.4, the locus check of 2.2, then the experiment table of part 4, and writes `root_locus.png` (the locus with poles, asymptotes, the ζ = 0.5 line and the design point; the simulated step responses for T_d = 0, 1, 2, 4 s).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-09-root-locus`. The 24 tests check:

1. **The conditions** (independent check 1) at every one of 300 gains for three plants: the altitude P loop, the altitude loop with T_d = 2 s, and a plant with complex poles. The angle error is below 10⁻⁶ degrees, K\|G\| = 1 to 10⁻⁶, and `gain_at_point` returns K. Four points off the locus fail by more than 1°.
2. **Branches**: start at the open-loop poles, equal the `np.roots` set at every 50th gain, no step larger than 0.15 between neighbouring gains.
3. **Root sum** constant at −3; **asymptotes** at −1 with 60°, 180°, 300°, and the roots at K = 10⁶ lie in those directions within 0.5°; two asymptotes at ±90° from −1.25 with the zero.
4. **Real axis**: the segment rule and the angle condition agree at 141 points for all three plants.
5. **Breakaway** at −1 + 1/√3, a double root; −1.577 rejected.
6. **Routh**: equal to the `np.roots` count for six random quartics; the ε case (2 unstable roots), the row-of-zeros cases (s² − 4)(s + 1) (1 unstable root) and (s² + 1)(s + 1) (none), and the auxiliary-derivative entry 2.
7. **Crossing** (independent check 2): K_c = a₂a₁/(a₃V) and ω = √2 from `jw_crossings`; Routh count 0 just below K_c and 2 just above; ω from the auxiliary polynomial; a root within 10⁻⁶ of jω at K_c. With T_d = 2 s there is no crossing and Routh reports 0 unstable roots up to K = 1000.
8. **Damping target** for ζ = 0.3, 0.5, 0.7: a root on the ray, stable; the closed form K = 0.04 rad/m and ω_n = 2 rad/s for T_d = 1 s.
9. **Step response** against `scipy.signal.step` to 10⁻⁶ and the exact 16.3 % overshoot of the cancelled design.

### Compare

There is no root-locus routine in `scipy.signal` (and the lesson avoids the python-control package by design). The independent checks are the two conditions, which are the definition of the locus, and Routh, which never computes a root. `np.roots` computes the eigenvalues of the companion matrix and is accurate to about 10⁻¹⁴ relative for well-separated roots; near the breakaway (a double root) its accuracy degrades to about the square root of that, yet the angle condition still holds to 10⁻¹³ degrees, because evaluating G at a slightly wrong double root still gives a value very close to −1/K.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Sketch the locus by hand first: poles, segments, asymptotes, breakaway, crossing. Then write each prediction down before running the code.

1. **Climb-rate feedback sweep.** For T_d = 0, 1, 2 and 4 s, use `pd_design_table` to find the gain for ζ = 0.5, the pair, the third pole, the predicted overshoot and 2 % settling time from the pair, and the simulated values. Predict which designs the dominant-pair formulas will describe well, from where the third pole sits relative to the pair and to the zero at −1/T_d.
2. **Sensitivity to the model.** Keep K = 0.01037 rad/m and change the inner-loop lag τ to 0.4 and 0.7 s, then the airspeed V to 40 and 60 m/s. Record ζ of the complex pair and the critical gain. Predict the direction of each change from the locus.

Raw data format: CSV tables `Td_s, K_rad_per_m, pair_re, pair_im, third_pole, pred_os_pct, pred_ts_s, sim_os_pct, sim_ts_s`; `parameter, value, zeta_pair, K_crit_rad_per_m`, each with a header line naming the plant.

What you should find:

| T_d | K for ζ = 0.5 | pair (ω_n) | third pole | predicted OS, t_s | simulated OS, t_s |
|---|---|---|---|---|---|
| 0 | 0.01037 rad/m | −0.333 ± 0.577j (0.667) | −2.333 | 16.3 %, 12.00 s | 15.5 %, 12.49 s |
| 1 s | 0.04000 | −1.000 ± 1.732j (2.000) | −1.000 | 16.3 %, 4.00 s | 16.3 %, 4.04 s |
| 2 s | 0.02851 | −1.284 ± 2.224j (2.568) | −0.433 | 16.3 %, 3.12 s | 5.4 %, 4.40 s |
| 4 s | 0.01597 | −1.398 ± 2.421j (2.796) | −0.204 | 16.3 %, 2.86 s | 0.0 %, 11.18 s |

With T_d = 1 s the zero at −1 cancels the pole at −1 exactly. The closed loop is then a true second-order system, s² + 2s + 100K, and the prediction is exact. With T_d = 2 and 4 s the third pole lies close to the zero (−0.433 near −0.5, −0.204 near −0.25) and well to the right of the pair. It is slow and only partly cancelled, so it dominates the tail: the overshoot vanishes but the settling time grows to 11 s. A ζ = 0.5 pair is not a ζ = 0.5 response. For the sensitivity study, ζ = 0.536, 0.500 and 0.438 for τ = 0.4, 0.5 and 0.7 s, with critical gains 0.070, 0.060 and 0.0486 rad/m; and ζ = 0.595, 0.500 and 0.429 for V = 40, 50 and 60 m/s, with critical gains 0.075, 0.060 and 0.050 rad/m. A slower inner loop or a faster aircraft pushes the design toward the crossing.

## 5. Validate and improve

| quantity | analysis (§2) | code | independent check | difference |
|---|---|---|---|---|
| centroid, angles | −1; 60°, 180°, 300° | same | roots at K = 10⁶ point along them | < 0.5° |
| breakaway | −1 + 1/√3 = −0.4226 | −0.4226 | double root at K = 0.003849 | < 10⁻⁴ |
| crossing | K = 0.06, ω = √2 | 0.060000, 1.414214 | Routh sign change; auxiliary root | < 10⁻⁸ |
| locus points | ∠G = 180°, K\|G\| = 1 | 1800 points | — | 1.4 × 10⁻¹³ °, 3.2 × 10⁻¹² |
| ζ = 0.5 gain, P | — | 0.01037 rad/m | root on the ray | < 10⁻⁶ |
| ζ = 0.5 gain, T_d = 1 s | 100K = ω_n² = 4: 0.04 | 0.040000 | closed form | < 10⁻⁸ |
| overshoot, P | 16.3 % (pair) | — | simulated 15.5 % | third pole at −2.33 |
| overshoot, T_d = 2 s | 16.3 % (pair) | — | simulated 5.4 % | slow pole −0.43 near zero −0.5 |

The locus computation is exact to rounding. What fails is an approximation layered on top of it, the dominant-pair formulas, and the simulation shows when. These are statements about the linear model (1); a flight test would add pitch limits, a richer inner loop and turbulence. For this model the recommendation is **K = 0.0104 rad/m (0.59 deg/m) with no rate feedback for a gentle, well-predicted response (15.5 % overshoot, 12.5 s), or K = 0.040 rad/m with T_d = 1 s for a four-second response (16.3 %, 4.0 s) if a clean climb-rate signal is available**, the latter relying on the inner-loop lag being where the model puts it.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand sketch made before running the code; the derivations of (3), the asymptote centroid, the breakaway condition and the Routh table at the critical gain; `root_locus.png`; the two experiment tables; and a half-page recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) of a gain and, if used, T_d, with predicted and simulated metrics and the parameter that most threatens the design.

### Adversarial self-checks

- **What if the airspeed rises by 20 %?** V multiplies the loop gain, so the locus is unchanged and only the gain labels along it change: the same K now sits 20 % further along it. ζ falls from 0.50 to 0.43 and the critical gain drops to 0.050 rad/m.
- **What if K → 0 or K → ∞?** At K → 0 the poles sit at 0, −1, −2 and the altitude error never decays (the pole at 0). At K → ∞ two poles head off along ±60° into the right half plane; the loop becomes unstable at K = 0.06 rad/m.
- **Which assumption fails first?** The first-order inner pitch loop. A second-order inner loop adds phase lag near 1 rad/s and moves the crossing to a lower gain; then pitch saturation for large altitude errors.
- **What measurement could falsify the model?** Raise the gain in simulation or flight until the altitude oscillation stops decaying; its period should be 2π/1.414 = 4.4 s at 0.06 rad/m. A different period or gain means the inner loop or the path lag is not as modelled.
- **Does the result survive a change of solver?** The angle condition, Routh and `np.roots` agree; the step responses match `scipy.signal.step` to 10⁻⁶.
- **Which decision changes if τ is uncertain by ±40 %?** At τ = 0.7 s the fixed P gain gives ζ = 0.44 and a critical gain of 0.049 rad/m, still a margin of 4.7. The T_d = 1 s design relies on cancelling the pole at −1/T; if T is uncertain, the cancellation is imperfect, and a slow residual pole appears as in the T_d = 2 s row.

## Check yourself

1. State the angle and magnitude conditions. Which one decides whether a point is on the locus?
2. For G = 1/(s(s + 1)(s + 2)), which parts of the real axis are on the locus for K > 0?
3. Give the asymptote angles and the centroid for that G.
4. Where is the breakaway point, and why is the other root of N D′ − N′D = 0 rejected?
5. Build the Routh array for s³ + 3s² + 2s + K. For what K is the loop stable, and at what frequency does it oscillate at the limit?
6. Why does the sum of the closed-loop poles stay constant when the pole excess is at least two?
7. A design places a pair at ζ = 0.5 but the step response has no overshoot. Give a reason found in this lesson.
