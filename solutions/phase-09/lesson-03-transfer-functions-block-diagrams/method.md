# P09.03  Method

The route to each answer, without the final numbers. `worked-solutions.md` has the numbers for s = 15.

## 1. Three connections

Series: 2(s + a)/((s + 1)(s² + 2s + 5)); expand the denominator. Parallel: [2(s² + 2s + 5) + (s + a)(s + 1)]/((s + 1)(s² + 2s + 5)). Loop: 2(s + a)/[(s + 1)(s² + 2s + 5) + 2(s + a)]. Poles: roots of the denominators (the series and parallel share them); zeros: roots of the numerators; DC gain: constant terms. Check with the code.

## 2. Derivations

(a) E = R − HY, Y = GE ⇒ Y(1 + GH) = GR. Positive feedback: E = R + HY ⇒ T = G/(1 − GH). (b) T = n_G d_H/(d_G d_H + n_G n_H): the numerator contains d_H; dividing by n_G n_H, T → d_H/n_H = 1/H when n_G grows (G = K n/d, K → ∞). (c) A = [[0, 1], [−a₀, −a₁]], B = (0, 1)ᵀ; solve (sI − A)v = B: v = (1, s)/(s² + a₁s + a₀); C = (b₀ − b₂a₀, b₁ − b₂a₁), D = b₂; Cv + D, put over the common denominator. (d) Pick-off moved before G: the branch must still carry G·x, so insert G in the branch. Summing junction moved after G: the added signal must still pass through G, so insert G in the side input.

## 3. Reaction wheel

Inner: wheel × body = 1/((τ_w s + 1)Js); feedback K_d: 1/((τ_w s + 1)Js + K_d). Times 1/s, times K_p, unity feedback: K_p/(Jτ_w s³ + Js² + K_d s + K_p). Poles by `poles`; DC gain 1. `simulate_diagram` with blocks wheel, body, angle and `connect = lambda y, r: {"wheel": Kp*(r - y["angle"]) - Kd*y["body"], "body": y["wheel"], "angle": y["body"]}`; compare with `simulate_tf`. Metrics from the simulated angle (peak, last exit from the ±2 % band) and the wheel output (peak torque). Relate the settling time to the slowest pair's real part.

## 4. Pitch damper

Loop K_q; `reduce_pitch_loop(dict(p, K_q=Kq, K_theta=Kt))["closed"]`; `poles(...).real.max()`; `simulate_tf` on linspace(0, 8, 1601); overshoot = max − 1; settling = time after the last sample outside ±0.02. Compare the K_q of least settling with the K_q of most negative max Re; the residues (how strongly each mode is excited) and the zeros decide the first, the poles alone the second. Confirm with `simulate_diagram(pitch_blocks(q), pitch_connect(q), ...)`.

## 5. Judgment

(a) `series(ctrl, plant)` then `minreal`. Unreduced: blocks `ctrl_lag` = 1/(s + 2) and `plant` = 1/(s − 1); wiring u = r − 3 y_lag, plant input u (+ d). Sample the plant output at 5, 10, 20 s with and without d. The unstable mode is unobservable in the reduced TF, not absent; any disturbance or initial condition excites e^{t}. (b) `np.roots(np.poly(-np.ones(n)))`; error against −1. (c) 25 deg / 2.295. (d) For example: actuator saturation and rate limits; sampling and computation delay; sensor noise and bias; changes of the aerodynamic coefficients with speed.

## 6. Altitude hold

(a) h/θ_c = T(s) · V(1/T_θ2)/(s(s + 1/T_θ2)); loop L = K_h T hθ; h/h_c = L/(1 + L): order 5 + 2 = 7, DC gain 1 (integrator in the loop). The factor (s + 1/T_θ2) appears in n_T and in d_hθ, so it divides both terms of d_T d_hθ + K_h n_T n_hθ: a closed-loop pole for every K_h, cancelled in h/h_c. (b) `poles(feedback(series(gain(Kh), T, hth))).real.max()` for each K_h; bisection between a stable and an unstable gain; the crossing frequency is the imaginary part of the pole that reaches the axis. (c) Extend the pitch blocks with `gamma_lag` = V(1/T_θ2)/(s + 1/T_θ2) driven by θ and `alt` = 1/s driven by the lag; θ_c = K_h(h_c − h) inside the wiring. Compare with `simulate_tf` of the reduced h/h_c for h_c = 10 m; metrics from the simulated altitude; initial pitch command K_h × 10 m. (d) Compare with the pitch authority of a small UAV and the comfort of a 10 m step; add a pitch-command limit and a climb-rate limit.
