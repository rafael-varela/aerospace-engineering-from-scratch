# P00.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. t = √(2h/g). Relative uncertainties u_h/h and u_g/g; u_t/t = ½ √(sum of squares).
2. Multiply by t; round u_t to two significant figures and t to the same place.
3. The larger squared term dominates. For equal contributions set u_h/h = u_g/g and solve for u_h.

## 2.
1. (a) ½g t² + v₀ t − h = 0; quadratic formula; keep the root that is positive for v₀ = 0.
2. (b) Differentiate t(v₀) with respect to v₀; the derivative of √(v₀² + 2gh) is v₀/√(…), which vanishes at v₀ = 0, leaving −1/g. Multiply by 0.1 m/s.
3. (c) Δt_push = −v₀/g: no h. Δt_drag = √(2h/g) · gh/(6v_t²) ∝ h^{3/2}. Over five heights one is flat and one rises steeply; fit a constant and a power law to the residuals t_meas − t_pred and see which fits. Averaging across heights brings each shift below the single-drop resolution.

## 3.
1. (a) Same formula with g = 3.71.
2. (b) v_t = √(2mg/(ρC_D A)); then gh/(6v_t²). Compare with Earth's v_t ≈ 23 m/s and 0.55 % correction.
3. (c) u_frame is fixed in seconds; divide by t at each planet; longer t means smaller relative uncertainty.
4. (d) h = ½ g (1 s)²; compare with the 2 m height row of the envelope.

## 4.
1. (a) `compare(H, 0.005, timings, u_single=frame_timing(1, fps)[1])`; tabulate; `plot_comparison`.
2. (b) Regenerate or re-film at 120 fps; the frame term falls by four; the Type A term is similar; E_n values stay below 2 if the data are honest.
3. (c) Subtract 0.15 s from every timing before `compare`; all E_n ≈ 0.15 / u_meas, all with the same sign; the standard deviation is unchanged so more drops cannot fix a mean that is wrong.
4. (d) `fit_g(H, comp.t_meas)`; compare |g_fit − 9.81| with u_g (and with the frame quantisation, which is systematic).

## 5.
1. (a) A = πd²/4; v_t = √(2mg/(ρC_D A)).
2. (b) t₀ at each height; first-order Δt = t₀ gh/(6v_t²); exact Δt = `fall_time_with_drag(h, v_t) − t₀`; compare at 1.8 m.
3. (c) E_n = Δt / √(u_meas² + u_pred²) with u_pred ≈ 1.5 ms, for each height and u_meas.
4. (d) Gap = Δt(ping-pong) − Δt(tennis) at 1.8 m; multiply by 120 for frames. Describe the frame: both balls released by one ruler tilt, camera level with the floor, impact frames read for each ball.

## 6.
1. (a) State what E_n < 2 with a large denominator means; the improvement is to shrink u_meas (video instead of stopwatch).
2. (b) Candidates: drag (rises with h^{3/2}), far marks mis-measured (random by height), impact frame rule late (constant). Only one rises with h. First test: compute the ping-pong/tennis drag numbers or re-measure the 1.8 m mark.
3. (c) A constant length error in h is a larger fraction at small h; t_meas shorter means the true drop is shorter than the recorded h. Re-measure to the bottom of the ball.

## 7.
1. (a) t_wrong = √(2h_c/g), t_right = √(2(h_c − r)/g); ratio = √(h_c/(h_c − r)).
2. (b) (ratio − 1) × t_right at each height, in ms.
3. (c) Divide by √(u_meas² + u_pred²). Signature: discrepancy shrinks with h, prediction longer than measurement.
