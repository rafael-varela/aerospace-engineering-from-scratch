# P01.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; symbolic results are exact, numerical ones are given to four significant figures and come from `solution.py`.

## 1. Your own glider

**Assumptions.** Steady glide, parabolic polar, sea-level ISA density, small angles for V* and the range.

**Derivation** (s = 15: m = 375 kg, S = 11.5 m², AR = 18.75, C_D0 = 0.015). W = 375 × 9.80665 = 3677.5 N. (a) K = 1/(π × 0.9 × 18.75) = 0.018863. (b) C_L* = √(0.015/0.018863) = 0.8917; (L/D)max = 1/(2√(0.015 × 0.018863)) = 29.72. (c) V* = √(2 × 3677.5/(1.225 × 11.5 × 0.8917)) = 24.20 m/s; γ = arctan(1/29.72) = 1.93°. (d) Sink = 24.20/29.72 = 0.814 m/s (small angle); exact 0.813 m/s. (e) Range = 1000 × 29.72 = 29.72 km.

**Result.** K = 0.01886, C_L* = 0.8917, (L/D)max = 29.72, V* = 24.20 m/s (24.19 exact), γ = 1.93°, sink 0.813 m/s, range 29.72 km.

**Check.** `best_glide_cl_newton(0.015, 0.018863)` = 0.891657; `airspeed(0.8917, 3677.5, 1.225, 11.5)` = 24.20; `sink_rate(...)` = 0.813; `glide_range(1000.0, ...)` = 29 720 m. Dimensional: K is dimensionless, V in m/s, range in m. Limiting case: compared with the lesson's glider (AR 20, same C_D0), the lower aspect ratio costs 3 % of L/D, consistent with (L/D)max ∝ √AR.

**Common mistakes.** Using mass instead of weight in V; forgetting π in K; quoting the range in metres as if it were kilometres.

## 2. Minimum sink and the sensitivity's error term

**Derivation.** (a) w = V sin γ ≈ V tan γ = V C_D/C_L with V = √(2W/(ρS)) C_L^(−1/2) (small angle, L = W), so w = √(2W/(ρS)) · C_D C_L^(−3/2). The prefactor is fixed, so minimise g(C_L) = (C_D0 + K C_L²) C_L^(−3/2). Product rule: g′ = 2K C_L · C_L^(−3/2) + (C_D0 + K C_L²)(−3/2) C_L^(−5/2) = C_L^(−5/2)[2K C_L² − (3/2)C_D0 − (3/2)K C_L²] = C_L^(−5/2)[(1/2)K C_L² − (3/2)C_D0]. Zero when K C_L² = 3C_D0, so C_L,ms = √(3C_D0/K). It is a minimum because g → ∞ as C_L → 0 (C_L^(−3/2) blows up) and as C_L → ∞ (K C_L^(1/2) grows), and this is the only critical point. Equivalently, the bracket (1/2)K C_L² − (3/2)C_D0 is increasing in C_L and changes sign from − to +, so g′ does too: the first-derivative test for a minimum.

(b) At C_L*: K C_L*² = C_D0, C_D = 2C_D0. At C_L,ms = √3 C_L*: K C_L,ms² = 3C_D0, C_D = 4C_D0. Ratio w*/w_ms = [2C_D0 C_L*^(−3/2)]/[4C_D0 (√3 C_L*)^(−3/2)] = (1/2) · 3^(3/4) = 1.1398. For the glider the code gives 0.7835/0.6872 = 1.1401; the 0.03 % difference is the cos γ factors dropped in the small-angle derivation.

(c) f(C_D0) = ½ K^(−1/2) C_D0^(−1/2). Power rule: f′ = −¼ K^(−1/2) C_D0^(−3/2) = −f/(2C_D0). Again: f″ = (3/8) K^(−1/2) C_D0^(−5/2) = (3/4) f/C_D0². Taylor to second order: f(C_D0 + Δ) = f + f′Δ + ½ f″Δ² + …, so the linear approximation's error is −½ f″ Δ² = −(3/8) f (Δ/C_D0)². For the glider, f = 30.70, δ = 0.1: −(3/8)(30.70)(0.01) = −0.115. The code reports −0.106; the remaining 0.009 is the third-order term, −(1/6)f‴Δ³ with f‴ = −(15/8)f/C_D0³, which gives +0.0096. ✓

**Check.** Units: f′ is in (L/D) per unit C_D0, both dimensionless. The numerical second derivative (central difference with 0.1 % steps) is 102 333 against (3/4)(30.70)/0.015² = 102 333. ✓

**Common mistakes.** Minimising C_D/C_L instead of C_D/C_L^(3/2) (gives best glide again); dropping the ½ in the second-order term; comparing w*/w_ms with 3^(3/4) = 2.28 and forgetting the factor ½ from C_D = 2C_D0 against 4C_D0.

## 3. A jet's best-range lift coefficient

**Derivation.** (a) h(C_L) = C_L^(1/2)/(C_D0 + K C_L²). Quotient rule numerator: ½ C_L^(−1/2)(C_D0 + K C_L²) − C_L^(1/2) · 2K C_L. Multiply by 2C_L^(1/2): C_D0 + K C_L² − 4K C_L² = C_D0 − 3K C_L² = 0, so C_L,jet = √(C_D0/(3K)) = C_L*/√3 = 0.5317 for the glider's polar (C_L* = 0.9210). It is a maximum because h(0) = 0 and h → 0 as C_L → ∞ (l'Hôpital: C_L^(1/2)/(K C_L²) → 0), with one critical point.

(b) At C_L,jet: C_D = C_D0 + C_D0/3 = (4/3)C_D0; L/D = (C_L*/√3)/((4/3)C_D0) = (√3/4)(C_L*/C_D0) = (√3/4) · 2(L/D)max = (√3/2)(L/D)max = 0.8660 (L/D)max. For the glider: 0.8660 × 30.70 = 26.59. (The same value as at the minimum-sink point, which is no accident: C_L,ms = √3 C_L* and C_L,jet = C_L*/√3 are mirror images about C_L* on a log scale, and L/D is symmetric in log C_L about its maximum.)

(c) V ∝ C_L^(−1/2), so V_jet/V* = (C_L*/C_L,jet)^(1/2) = 3^(1/4) = 1.316. For the glider's weight and area: V_jet = 31.68 m/s against V* = 24.07 m/s. The ratio depends only on the exponent of C_L in the range expression, not on C_D0, K, W or ρ.

(d) With h′(C_L) written out and `newton(hprime, hprime_derivative, 0.5)`, the root 0.53174 is found in 3 iterations; from 0.3 in 4. (A central difference of h′ as the derivative also works, at the cost of one more iteration.)

**Check.** C_L,jet × C_L,ms = C_L*², so the three optima are in geometric progression with ratio √3. ✓ Units: C_L dimensionless, speeds in m/s.

**Common mistakes.** Maximising C_L^(1/2)/C_D by setting dC_D/dC_L = 0; forgetting that the jet-range optimum is a *faster* condition than best glide (lower C_L); concluding that L/D = 0.866 (L/D)max is "13 % worse" without noting that the range expression is different and this is the optimum for it.

## 4. Wind, and Newton's basin

**Result** (the glider at sea level; V* = 24.07 m/s, min-sink speed 18.29 m/s).

| V_w (m/s) | V_best (m/s) | C_L | ground L/D |
|---|---|---|---|
| −20 | 21.27 | 1.180 | 57.79 |
| −15 | 21.70 | 1.134 | 50.82 |
| −10 | 22.26 | 1.077 | 43.95 |
| −5 | 23.02 | 1.007 | 37.22 |
| 0 | 24.07 | 0.921 | 30.70 |
| 5 | 25.59 | 0.815 | 24.52 |
| 10 | 27.81 | 0.690 | 18.87 |
| 15 | 31.03 | 0.554 | 14.01 |
| 20 | 35.42 | 0.425 | 10.18 |

(a) The speed increase is 1.52, 3.74, 6.96, 11.35 m/s for 5, 10, 15, 20 m/s of headwind: 0.30, 0.37, 0.46, 0.57 of the wind. A straight-line fit through the four headwind points has slope 0.65 and intercept −2.3 m/s, which is not a useful rule; "add 0.4 × headwind" is within 1 m/s for 5–15 m/s, while "add half the headwind" is within 1 m/s only near 15 m/s: it is 1.0–1.3 m/s too fast at 5–10 m/s and 1.3 m/s too slow at 20 m/s. The increase is convex in V_w because the polar steepens as V³.

(b) At −20 m/s the best speed is 21.27 m/s, still 3 m/s above minimum sink, and the C_L is 1.18 against 1.60. The tangent from a point (V_w, 0) with V_w < 0 must touch the polar where the polar's slope equals the slope of the line to that point, which is negative (sink increases away from the point); left of the minimum-sink speed the polar slopes the other way, so no tangent from any axis point can touch there. Flying slower than minimum sink always loses both time and distance.

(c) Starts 0.2–1.3 converge to C_L* = 0.9210 (4–8 iterations); 1.4 converges to −0.9210 (the inverted-flight minimum, a genuine root of (4)); 0.1 and every start from 1.5 upward raise `RuntimeError`. The basin's upper edge sits just below the inflection point √(3C_D0/K) = 1.595, where (L/D)″ = 0 and the Newton step becomes infinite; a start at 1.4 or 1.5 is thrown far across the root by the small denominator. The lower edge is near 0.1, where f″ → 0 again as C_L → 0.

(d) From 0.8: |error| = 0.121, 0.0193, 5.9 × 10⁻⁴, 5.6 × 10⁻⁷, 5.1 × 10⁻¹³, then 0 to machine precision. Correct digits 1, 2, 3, 6, 12, 16: doubling once the iterate is close.

**Check.** At V_w = 0 the table reproduces (L/D)max = 30.70 and V* exactly. The headwind entries all have C_L < C_L* and the tailwind entries C_L > C_L*, as the tangent construction requires.

**Common mistakes.** Fitting one straight line through tailwind and headwind points together; treating the −0.921 root as a convergence failure (it is a success of the method on the wrong basin); reporting the basin edge as exactly the inflection point (Newton fails a little before it).

## 5. Can we reach the field?

**Assumptions.** Uniform headwind during the descent; the pilot flies the headwind-optimum airspeed for the *forecast* wind; sea-level density throughout (a 2 % error on speeds over 1500 m, nil on ratios).

**Derivation.** (a) R = d/h, quotient rule: dR/dt = (ḋ h − d ḣ)/h² = [−(V_b − V_w) h + d w]/h² = (w/h)[d/h − (V_b − V_w)/w]. Since w/h > 0, dR/dt < 0 exactly when (V_b − V_w)/w, the ground glide ratio, exceeds the required ratio d/h.

(b) d/h = 26 000/1500 = 17.33. With V_w = 7: V_b = 26.37 m/s, w = 0.873 m/s, G = 22.18, dR/dt = −2.82 × 10⁻³ s⁻¹. V_w = 10: V_b = 27.81, w = 0.944, G = 18.87, dR/dt = −9.66 × 10⁻⁴ s⁻¹. V_w = 13: V_b = 29.60, w = 1.048, G = 15.84, dR/dt = +1.04 × 10⁻³ s⁻¹. The field is reached for 7 and 10 m/s; at 13 m/s the required ratio grows with time and the glider lands short.

(c) Height needed d/G: 1172 m (7 m/s), 1378 m (10 m/s), 1641 m (13 m/s). Time d/(V_b − V_w): 22.4, 24.3, 26.1 minutes.

(d) The margin at the forecast wind is 122 m (8 %), and the pessimistic wind needs 141 m more than is available. A rule that makes the decision: *commit to a field only if the pessimistic wind still reaches it*. By that rule the upwind field is rejected and the downwind field (12 km needed; the tailwind ground range from 1500 m at −10 ± 3 m/s is at least 1500 × 37 = 55 km) is chosen with a margin of a factor four. The derivative dR/dt is what the pilot would *see*: with the upwind field the ratio of distance-to-go to height would be falling only slowly, and any gust would flip its sign.

**Check.** At the forecast wind the time to the field, 24.3 min, times the sink 0.944 m/s gives 1376 m, equal to the height needed (1378 m) to rounding. ✓

**Common mistakes.** Using the still-air (L/D)max × h = 46 km as the range; flying the still-air best-glide speed in the headwind (G = 17.95 at 10 m/s: the margin shrinks to 1 %); deciding on the forecast case instead of the pessimistic one.

## 6. Where the model stops

(a) **Parabolic polar fails near stall.** The minimum-sink C_L of 1.6 exceeds C_Lmax = 1.4 by 14 %; the wing cannot produce it. The real minimum sink is at or just below C_Lmax, with a sink rate higher than the model's 0.687 m/s because C_D rises steeply there. Observable: stall warning or buffet before the computed minimum-sink airspeed is reached, and a vertical-speed reading that worsens as the airspeed is reduced below about 19 m/s.

(b) **Steady glide fails.** The excess kinetic energy ½(40² − 24²)/g = 52 m of height is recovered only if the pilot zooms gently; a rapid deceleration with the speed brakes loses it to drag. During the 10–20 s transition the sink rate is not the polar's. Observable: airspeed decaying while the altimeter barely moves, then settling.

(c) **Uniform wind fails.** The optimum airspeed for 10 m/s is 27.8 m/s with G = 18.9; for 4 m/s it is 25.3 m/s with G = 25.7. Flying 27.8 m/s in the 4 m/s wind near the ground gives G = 25.2, 2 % below the optimum, and flying 25.2 m/s up high in 10 m/s gives 18.5, 2 % below; the loss is small because optima are flat, but the *range estimate* made at 1000 m with 10 m/s is pessimistic by about 20 % over the lower half of the descent. Observable: groundspeed (GPS) increasing as the glider descends at constant airspeed.

(d) **Fixed weight.** Weight does not appear in C_L*, (L/D)max or the range, so the still-air range is unchanged, but every speed scales by √(420/350) = 1.095: best glide 26.4 m/s instead of 24.1, sink 0.858 m/s instead of 0.783. In a headwind the heavier glider does *better* (its polar is shifted to higher speeds, closer to the wind point). Observable: the measured sink at the old best-glide airspeed is higher than the polar says, because that airspeed is now below the optimum.

**Common mistakes.** Treating the minimum-sink prediction as safe because the arithmetic is right; forgetting that kinetic energy at cruise is worth real height; assuming ballast shortens the glide.

## 7. Sinking air

**Derivation.** (a) G(V) = V/(w(V) + u). G′ = [(w + u) − V w′]/(w + u)² = 0 ⟹ w + u = V w′(V). This is (11) with V_w = 0 and w → w + u. Geometrically: shift the polar down by u (the air itself sinks) and draw the tangent from the origin; the touching point is the best speed.

(b) Copy `headwind_best_speed`, write `f = lambda v: (w(v) + u) - v * w1(v)` and keep `df`. Results for the glider:

| u (m/s) | V_best (m/s) | C_L | total sink (m/s) | ground L/D |
|---|---|---|---|---|
| 0 | 24.07 | 0.921 | 0.784 | 30.70 |
| 0.5 | 27.61 | 0.700 | 1.434 | 19.26 |
| 1.0 | 30.63 | 0.569 | 2.116 | 14.48 |
| 2.0 | 35.57 | 0.422 | 3.531 | 10.08 |

(c) 1 m/s of sinking air cuts the range from 1000 m from 30.7 km to 14.5 km, a loss of 16.2 km; a 10 m/s headwind cuts it to 18.9 km, a loss of 11.8 km. Sinking air of only 1 m/s is worse than a 10 m/s headwind, because it adds to a sink rate that is itself only 0.78 m/s.

(d) The condition w + u = V w′ must hold at the optimum; adding u > 0 raises the left side, so V w′(V) must be larger. Above best glide w′ = 3aV² − b/V² is positive and increasing, so V w′ grows with V, and the equality is restored only at a higher V: the glider must fly faster to spend less *time* in the sinking air, which the derivative condition encodes without any appeal to intuition.

**Check.** u → 0 reproduces 24.07 m/s and 30.70. ✓ For u = 1 the C_L of 0.569 is close to the jet-range optimum C_L*/√3 = 0.532 of Problem 3; both are "fly faster than best glide" conditions.

**Common mistakes.** Adding u to the airspeed instead of to the sink; keeping V_w = 10 in the function while asking about still horizontal wind; using the still-air speed and concluding the range is 1000/(0.784 + 1) × 24.07 = 13.5 km (close, but 7 % short of the optimum 14.5 km, which is the whole point).

## Check-yourself answers

1. tan γ = D/L = C_D/C_L; weight cancels out of the ratio. A heavier glider flies the same angle at a higher speed (V ∝ √W), so it sinks faster but not more steeply.
2. Equal: K C_L*² = C_D0, so C_D = 2C_D0 at best glide.
3. If f′(c) = 0 and f″(c) < 0 then f has a local maximum at c. At C_L*, f″ = −K C_L*/(2C_D0²) < 0, so (L/D) is maximal there.
4. The limit is of the form ∞/∞; l'Hôpital gives 1/(2K C_L) → 0. With f(0) = 0 as well, the positive function f vanishes at both ends of (0, ∞) and its single critical point must be the maximum.
5. Lower: V ∝ C_L^(−1/2), so the minimum-sink speed is 3^(−1/4) = 0.760 times the best-glide speed (18.3 against 24.1 m/s for the glider).
6. (L/D)max(C_D0 + Δ) ≈ (L/D)max − (L/D)max Δ/(2C_D0); the error is ½f″Δ² = (3/8)(L/D)max (Δ/C_D0)², so 0.36 % of (L/D)max for a 10 % change in C_D0 (0.115 against the exact 0.106).
7. The chain rule, dV*/dt = (dV*/dρ)(dρ/dh)(dh/dt). In a descent dh/dt < 0 and dρ/dh < 0 while dV*/dρ < 0, so dV*/dt < 0: the best-glide speed falls as the air thickens.
8. Quadratic convergence (the error is squared each step, digits double). Starting at C_L = 2 lies beyond the inflection point of L/D at √(3C_D0/K) = 1.6, where (L/D)″ = 0; the first Newton step is huge and the iteration diverges, so the code raises `RuntimeError`.
