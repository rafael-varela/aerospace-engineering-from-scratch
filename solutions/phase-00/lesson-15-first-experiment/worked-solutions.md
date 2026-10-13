# P00.15  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty, at least one independent check, and the common mistakes. g = 9.807 m/s² unless given; standard uncertainties unless marked k = 2.

## 1. One prediction

**Assumptions.** Release from rest, no drag, inputs independent.

**Derivation.** t = √(2 × 1.500 / 9.81) = 0.5530 s. Relative terms: u_h/h = 0.010/1.5 = 0.00667; u_g/g = 0.03/9.81 = 0.00306. u_t/t = ½ √(0.00667² + 0.00306²) = ½ × 0.00734 = 0.00367; u_t = 0.00203 s.

**Result.** t_pred = 0.5530 ± 0.0020 s, or 0.553 ± 0.004 s (k = 2). The height dominates: its squared term (4.4×10⁻⁵) is almost five times the g term (9.4×10⁻⁶). Equal contributions need u_h/h = 0.00306, that is u_h = 4.6 mm: a careful tape reading to the ball's bottom, not a heroic one.

**Check.** `predict_fall_time(1.5, 0.010, 9.81, 0.03)` returns (0.55300, 0.00203). Dimensions: √(m / (m s⁻²)) = s.

**Common mistakes.** Forgetting the ½ from the square root (gives u_t = 4 ms). Quoting t to five decimals with a 2 ms uncertainty.

## 2. A hand that pushes

**(a)** ½ g t² + v₀ t − h = 0 gives t = (−v₀ ± √(v₀² + 2gh)) / g; the positive root is t = (−v₀ + √(v₀² + 2gh)) / g. At v₀ = 0 it reduces to √(2h/g) ✓.

**(b)** dt/dv₀ = (−1 + v₀/√(v₀² + 2gh)) / g. At v₀ = 0 the second term vanishes: dt/dv₀ = −1/g = −0.102 s per m/s, with no h in it. For v₀ = 0.1 m/s the fall is shorter by 10.2 ms. Exact check at h = 1.8 m: t(0.1) = (−0.1 + √(0.01 + 35.31)) / 9.81 = 0.5957 s against t₀ = 0.6058 s: 10.1 ms ✓.

**(c)** Δt_push = −v₀/g: constant in h. Δt_drag = t₀ · gh/(6v_t²) = √(2/g) · g/(6v_t²) · h^{3/2}: grows as h^{3/2}, tripling the height multiplies it by 5.2. Each is below a single frame at 30 fps (33 ms), but with five heights and six drops each you have thirty timings: fit the residuals t_meas − t_pred to a constant plus a term in h^{3/2}. A push shows as a negative constant with no h-dependence; drag as a positive term rising steeply toward 1.8 m. The two have opposite signs and opposite shapes, which is why they can be separated even when neither is visible in one drop.

**Common mistakes.** Keeping the negative root. Treating the push as random and expecting it to average out.

## 3. The same experiment on Mars

**Assumptions.** g_Mars = 3.71 m/s², ρ = 0.020 kg/m³, C_D = 0.5, same ball and camera.

**(a)** t = √(2h/3.71): 0.569 s at 0.6 m, 0.985 s at 1.8 m. Falls take 1.63 times longer (√(9.81/3.71)).

**(b)** A = π(0.067)²/4 = 3.53×10⁻³ m². v_t = √(2 × 0.058 × 3.71 / (0.020 × 0.5 × 3.53×10⁻³)) = 110 m/s (Earth: 23 m/s). Correction gh/(6v_t²) = 3.71 × 1.8 / (6 × 12 200) = 9.1×10⁻⁵, about 0.09 ms, sixty times smaller than Earth's 0.55 % (3.3 ms). Thin air makes free fall a far better model on Mars.

**(c)** u_frame = 13.6 ms either way. Relative to t at 1.8 m: 1.4 % on Mars versus 2.2 % on Earth. Sharper on Mars, because the same clock resolves a longer interval. In general a slower fall is an easier timing experiment; that is also why tall drops are better than short ones, up to the envelope.

**(d)** h = ½ × 9.81 × 1² = 4.9 m. The `P00.13` envelope caps drop height at 2 m (and stored energy at 40 J). Not permitted at home; a 4.9 m drop means a stairwell or a window, both on the "not permitted" list. The engineer's goal is reached instead by a faster camera.

**Check.** v_t scaling: v_t ∝ √(g/ρ); √((3.71/0.020)/(9.81/1.2)) = √(185.5/8.175) = 4.76, and 110/23.2 = 4.76 ✓.

**Common mistakes.** Using Earth's g with Mars's density. Forgetting the envelope when "making the timing easier".

## 4. Numerical, with the code

Synthetic, seed 0, six drops per height, u_h = 5 mm; your numbers will differ in detail.

**(a)** At 30 fps: u_meas ≈ 15 ms at every height (the 13.6 ms frame term dominates; the Type A part is 3 to 7 ms). E_n = 0.01, 0.04, 0.02, 0.20, 0.34: all far below 2, as they must be for data generated from the model. The plot shows measured points sitting on the √h curve with k = 2 bars of ±30 ms.

**(b)** At 120 fps: u_meas falls to about 3.7 ms. The frame term dropped by four (13.6 → 3.4 ms); the Type A term is similar in size, so it now matters. E_n = 0.3, 0.5, 0.3, 0.1, 0.7: still consistent, and the test is four times sharper. The same model survives a stricter test, which is the result worth reporting.

**(c)** Subtracting 0.15 s from every timing: E_n = 9.6, 10.2, 10.2, 9.8, 9.4, all with t_meas < t_pred. A constant bias of 150 ms over a combined uncertainty of 15 ms gives E_n ≈ 10 everywhere. More drops shrink the Type A term only; the frame term and the bias are untouched, and even with the Type A term at zero E_n would be 150/13.6 = 11. A systematic error has to be removed, not averaged.

**(d)** `fit_g` on the 30 fps synthetic means: g = 9.77 ± 0.06 m/s². |9.77 − 9.81| = 0.04 < u_g: consistent. The uncertainty from the fit is honest here because the frame quantisation enters the synthetic data at random phases; with real data read with a fixed frame rule, add the frame term as a systematic before claiming agreement at this level.

**Common mistakes.** Passing frame counts instead of seconds to `compare`. Reporting E_n without u_meas beside it.

## 5. Drag by design: the ping-pong ball

**Assumptions.** C_D = 0.5 (a sphere at Reynolds numbers of 10⁴ to 10⁵), ρ = 1.2 kg/m³, no spin.

**(a)** A = π(0.040)²/4 = 1.257×10⁻³ m². v_t = √(2 × 0.0027 × 9.807 / (1.2 × 0.5 × 1.257×10⁻³)) = 8.38 m/s.

**(b)**

| h (m) | t₀ (s) | first-order Δt (ms) | `fall_time_with_drag` Δt (ms) |
|---|---|---|---|
| 0.6 | 0.3498 | 4.88 | 4.90 |
| 1.2 | 0.4947 | 13.81 | 13.92 |
| 1.8 | 0.6059 | 25.38 | 25.67 |

At 1.8 m the first-order formula is 1.1 % low (gh/v_t² = 0.25 there, so the next term is small but not zero). Good enough for design; use the closed form for the comparison.

**(c)** E_n = Δt / √(u_meas² + u_pred²) with u_pred ≈ 1.5 ms:

| h (m) | 30 fps (u_meas 15 ms) | 120 fps (u_meas 4.5 ms) |
|---|---|---|
| 0.6 | 0.33 | 1.0 |
| 1.2 | 0.92 | 2.9 |
| 1.8 | 1.70 | 5.4 |

At 30 fps the drag-free model is never rejected (E_n < 2 everywhere): the ping-pong ball looks like free fall to a 30 fps camera. At 120 fps drag is detected at 1.2 m and decisively at 1.8 m. The experiment's power to test the assumption is set by the frame rate, and you knew that before dropping anything.

**(d)** At 1.8 m the tennis ball arrives 3.3 ms late and the ping-pong ball 25.7 ms late relative to free fall; the gap is 22.3 ms, or 2.7 frames at 120 fps. Design: rest both balls on one ruler at the 1.8 m mark, release by tilting the ruler sharply away (same instant, same v₀ ≈ 0 for both), camera level with the floor at 120 fps with both impact points in frame. Read the two impact frames; predict the tennis ball first by 2 to 3 frames. A single video then shows an assumption failing for one object and holding for another, which is the most convincing kind of validation.

**Check.** Ratio of corrections ping-pong/tennis at 1.8 m should be (v_t,tennis/v_t,pp)² = (23.2/8.38)² = 7.7; 25.7/3.3 = 7.8 ✓.

**Common mistakes.** Using the diameter in A = πd²/4 as if it were the radius (four times too large an area). Expecting the ping-pong ball to "float" visibly from 1.8 m: 26 ms is not visible to the eye, only to the camera.

## 6. What E_n does and does not tell you

**(a)** Learned: nothing in the model is wrong by more than roughly 2 × 50 ms = 100 ms, which is 15 to 30 % of the fall time. Not learned: whether drag, release speed, the height reference, or g itself are right at the few-percent level, because all of those produce effects of 3 to 25 ms that disappear inside a 50 ms uncertainty. E_n < 2 with a blunt instrument is agreement with almost anything. The single change: replace the stopwatch with frame counting, which cuts u_meas by a factor of three at 30 fps and twelve at 120 fps.

**(b)** t_meas longer at every height, discrepancy growing with h, crossing 2 at 1.8 m. Candidates, ranked: (1) drag, which produces exactly a positive shift rising as h^{3/2}, and whose predicted size for a tennis ball (3.3 ms at 1.8 m) is too small for E_n = 2.6 at u_meas ≈ 15 ms, so either the ball is lighter or larger than a tennis ball or this is not the whole story; (2) the tall marks measured too low (the tape sagged or was read at an angle), which makes t_pred too short at the tall heights only; a random-by-height pattern, but with the worst error plausibly at the top; (3) a late impact-frame rule, which is constant and does not fit the growth. First follow-up: re-measure the 1.5 and 1.8 m marks with the tape held taut and the ball's bottom at the mark, because it is five minutes' work and distinguishes (1) from (2) immediately.

**(c)** Shorter measured times, discrepancy largest at the smallest height and shrinking with h: the heights were recorded too large by a constant amount, so the ball actually fell less than h. The classic cause is measuring to the centre (or top) of the ball instead of its bottom; a 34 mm error is 5.6 % of 0.6 m but 1.9 % of 1.8 m, and in time 2.9 % versus 0.9 %. Fix: re-measure every height to the bottom of the ball at the release position, correct `data/raw` only by recording the corrected heights in a new file with a note, and re-run.

**Common mistakes.** Choosing the cause you understand best rather than the one whose signature matches the pattern.

## 7. Where is h = 0?

**(a)** t_wrong = √(2h_c/g) uses the centre height h_c; t_right = √(2(h_c − r)/g) uses the distance the bottom falls. Ratio: √(h_c / (h_c − r)) = (1 − r/h_c)^{−½} ≈ 1 + r/(2h_c) for small r/h_c.

**(b)** r = 33.5 mm. At 0.6 m: ratio 1.0291, error 2.9 % of t_right = 0.3397 s: 9.9 ms. At 1.8 m: ratio 1.0094, error 0.9 % of 0.6002 s: 5.7 ms. Prediction too long in both cases, more so at the short height.

**(c)** E_n with u_pred ≈ 1.5 ms: at 30 fps (u_meas 15 ms) 0.66 at 0.6 m and 0.38 at 1.8 m: the mistake is invisible, buried in the frame quantisation. At 120 fps (u_meas 4.5 ms) 2.1 at 0.6 m and 1.2 at 1.8 m: caught at the shortest height, borderline elsewhere. Signature: measured shorter than predicted, discrepancy decreasing with height, the mirror image of drag. If you see E_n highest at your *lowest* height, measure your ball's radius before anything else.

**Check.** r/(2h_c) at 0.6 m = 0.028, against the exact 0.029 ✓.

**Common mistakes.** Correcting by r in the time rather than the height. Concluding from E_n < 2 at 30 fps that the heights were right.

## Check-yourself answers

1. t = √(2h/g); assumptions: g constant over the drop, release from rest (v₀ = 0), no air drag (plus a rigid floor at a known h = 0).
2. √2 ≈ 1.41.
3. u_t/t = ½ √((u_h/h)² + (u_g/g)²).
4. 0.41/60 s = 6.8 ms.
5. E_n = |t_pred − t_meas| / √(u_pred² + u_meas²); consistent at k = 2 when E_n < 2.
6. No: a bias is the same on every drop, so the mean is wrong by 0.2 s however many drops are averaged. E_n is large at every height with t_meas longer than t_pred by roughly the same amount.
7. Drag lengthens the fall; Δt ≈ t₀ gh/(6v_t²); the ping-pong ball, whose terminal speed (8.4 m/s) is a third of the tennis ball's (23 m/s), shows it about eight times more strongly.
8. The frame term: 13.6 ms at 30 fps against u_pred ≈ 1.6 ms, nearly nine times larger. The comparison is limited by the camera, not the model.
