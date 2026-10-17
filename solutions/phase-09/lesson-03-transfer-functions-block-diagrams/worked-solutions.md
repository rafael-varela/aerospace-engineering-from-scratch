# P09.03  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (its `main()` or the functions named) or by a probe script that calls them; the "independent check" line says how it was verified by a different route. The pitch model is representative, not data for a specific aircraft.

## 1. Three connections (a = 3)

**Result.** (s + 1)(s² + 2s + 5) = s³ + 3s² + 7s + 5.

| connection | transfer function | poles | zeros | DC gain |
|---|---|---|---|---|
| series G₁G₂ | (2s + 6)/(s³ + 3s² + 7s + 5) | −1, −1 ± 2j | −3 | 6/5 = **1.2** |
| parallel G₁ + G₂ | (3s² + 8s + 13)/(s³ + 3s² + 7s + 5) | −1, −1 ± 2j | −1.3333 ± 1.5986j | 13/5 = **2.6** |
| loop G₁G₂/(1 + G₁G₂) | (2s + 6)/(s³ + 3s² + 9s + 11) | −1.6258, −0.6871 ± 2.5087j | −3 | 6/11 = **0.5455** |

The code returns the same polynomials and roots.

**Independent check.** Parallel DC gain = G₁(0) + G₂(0) = 2 + 0.6 = 2.6; loop DC gain = 1.2/(1 + 1.2) = 0.5455. Both from the block values alone.

**Common mistakes.** Adding numerators for the parallel connection without cross-multiplying; forgetting that the loop's poles are not the open-loop poles.

## 2. Derivations

**(a)** Negative feedback: E = R − HY, Y = GE ⇒ Y = GR − GHY ⇒ **T = G/(1 + GH)**. Positive feedback: E = R + HY ⇒ **T = G/(1 − GH)**.
**(b)** With G = n_G/d_G and H = n_H/d_H, T = n_G d_H/(d_G d_H + n_G n_H). The numerator is n_G d_H, so the zeros of T are the zeros of G and the roots of d_H, the poles of H (unless a cancellation occurs). For G = K n/d, T = K n d_H/(d d_H + K n n_H) → n d_H/(n n_H) = d_H/n_H = 1/H as K → ∞ (wherever n n_H ≠ 0).
**(c)** A = [[0, 1], [−a₀, −a₁]], B = (0, 1)ᵀ. (sI − A) = [[s, −1], [a₀, s + a₁]], determinant s² + a₁s + a₀ = d(s), and (sI − A)⁻¹B = (1, s)ᵀ/d. With D = b₂ and C = (b₀ − b₂a₀, b₁ − b₂a₁): C(sI − A)⁻¹B + D = [(b₀ − b₂a₀) + (b₁ − b₂a₁)s + b₂(s² + a₁s + a₀)]/d = (b₂s² + b₁s + b₀)/d.
**(d)** A pick-off point after G carries Gx; moved before G it carries x, so a block G must be inserted in the branch. A summing junction before G adds a signal v that then passes through G; moved after G, the side input must pass through a block G first so that the sum is still G(x + v).

## 3. Reaction-wheel attitude loop

**Reduction.** Inner: 1/((0.1s + 1)·10s) with K_d = 8 feedback: 1/(s² + 10s + 8) after normalising (Jτ_w = 1). Outer: θ/θ_c = **2/(s³ + 10s² + 8s + 2)**.
**Poles** −9.1495 and −0.4252 ± 0.1943j; **DC gain 1** (the body is a double integrator, so the loop has an integrator in its forward path).
**Simulation.** The unreduced diagram (three blocks, wiring τ_c = K_p(θ_c − θ) − K_d ω) and `simulate_tf` of the reduced TF agree to **2.2 × 10⁻¹⁶ rad** over 40 s.
**Metrics** (1 rad step): peak 1.0010 rad at 16.3 s (**overshoot 0.10 %**), **2 % settling 10.39 s**, **peak wheel torque 1.705 N m** per rad of command. The settling time is set by the slow, well-damped pair: its real part −0.425 s⁻¹ gives 4/0.425 ≈ 9.4 s, close to the measured 10.4 s.

**Independent check.** Final-value: θ(∞) = T(0) = 2/2 = 1.

**Common mistakes.** Writing the body as 1/(Js²) and then also integrating it (double-counting the integrator); feeding back θ instead of ω through K_d.

## 4. Tuning the pitch damper (K_θ = 2.5)

| K_q | 0.0 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| max Re (s⁻¹) | −0.388 | −0.933 | −1.017 | −1.002 | −0.986 | −0.970 | −0.954 | −0.938 | −0.921 | −0.904 | −0.887 |
| overshoot (%) | 37.75 | 22.21 | 9.51 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 |
| 2 % settling (s) | > 8 | 4.01 | 3.04 | 2.97 | 3.04 | 3.15 | 3.26 | 3.37 | 3.48 | 3.60 | 3.73 |

**(b)** The smallest settling time is at **K_q = 0.3** (2.97 s); the fastest slowest pole is at **K_q = 0.2** (−1.017 s⁻¹). They differ because settling depends on how strongly each mode is excited (its residue), not only on where it sits: at K_q = 0.2 the complex pair still produces 9.5 % overshoot that lingers outside the band, while at 0.3 the response approaches from below on the slow real pole alone. Above 0.3 the slow real pole moves back toward the axis and every metric worsens: too much damper makes the aircraft sluggish.
**(c)** For K_q = 0.3 the unreduced diagram and the reduced TF agree to **1.2 × 10⁻¹⁵ rad**.

**Common mistakes.** Choosing the gain by the poles alone; reporting "0 % overshoot" as "well damped" without looking at the slow tail.

## 5. Judgment

**(a)** series(ctrl, plant) = (s − 1)/((s + 2)(s − 1)) = (s − 1)/(s² + s − 2); `minreal` gives **1/(s + 2)**, which predicts a plant output of 0.5 for a unit step.

| plant output at | 5 s | 10 s | 20 s |
|---|---|---|---|
| unreduced, as drawn | 0.49998 | 0.50000 | 0.50000 |
| unreduced, 10⁻⁶ disturbance at the plant input | 0.50012 | 0.522 | **485.7** |

Exactly as drawn, from rest and in near-exact arithmetic, the unstable mode is never excited and the output looks like the reduced prediction. A disturbance of one part in a million at the plant input grows as e^{t}: 10⁻⁶ e^{20} ≈ 485. The cancellation removed the mode from the transfer function from r to y, not from the system: the loop is **internally unstable**. Never cancel an unstable pole with a controller zero.
**(b)** Largest root error: n = 4: **2.2 × 10⁻⁴**; n = 10: **0.051**; n = 16: **0.23**. A root of multiplicity n moves by about ε^{1/n} for a coefficient perturbation ε ≈ 10⁻¹⁶; for many identical or clustered blocks keep the product factored (multiply poles lists, not polynomials) and only expand at the end, or work in state space.
**(c)** 25 deg / 2.295 = **10.9 deg** of attitude step; beyond it the elevator saturates and the linear prediction is wrong (the response is slower and can overshoot more when the saturation ends).
**(d)** Any two of: elevator deflection and rate limits; the sampling and computation delay of the flight computer; gyro noise, bias and quantisation; variation of the short-period coefficients with speed, altitude and centre of gravity; the phugoid mode; structural modes picked up by the gyro.

## 6. Altitude hold (V = 30 m/s, 1/T_θ2 = 1.5 s⁻¹)

**(a)** h/θ = 45/(s(s + 1.5)), so **h/h_c = K_h T(s)·45/(s(s + 1.5)) / [1 + K_h T(s)·45/(s(s + 1.5))]**: order 5 + 2 = **7**, **DC gain 1** (the altitude integrator). T(s) = n_T/d_T contains the factor (s + 1.5) in n_T, and h/θ contains it in its denominator. The characteristic polynomial d_T s(s + 1.5) + 45K_h n_T then contains (s + 1.5) in both terms: **s = −1.5 is a closed-loop pole for every K_h**. It cancels in h/h_c (a command cannot excite it), but it is a real mode of the aircraft: the slow flight-path response to pitch, which no altitude gain can move.

**(b)**

| K_h (rad/m) | 0.001 | 0.005 | 0.01 | 0.02 | 0.05 |
|---|---|---|---|---|---|
| max Re of closed-loop poles (s⁻¹) | −0.031 | −0.186 | −0.515 | −0.495 | −0.434 |

The loop becomes unstable at **K_h = 0.2276 rad/m**, where a pair crosses the imaginary axis at **±2.745j rad/s**. The fastest slowest pole in a fine sweep is at K_h ≈ 0.0084 rad/m (−0.518 s⁻¹).

**(c)** For K_h = 0.01 rad/m and a 10 m step: the unreduced simulation (pitch blocks + γ-lag + altitude integrator, θ_c = K_h(h_c − h) in the wiring) agrees with the reduced TF to **6.2 × 10⁻¹⁴ m**. Peak **10.011 m** at 13.4 s (0.11 % overshoot), **2 % settling 8.58 s**, initial pitch command K_h × 10 = 0.1 rad = **5.73 deg**. Closed-loop poles −100.90, −17.24, −2.01 ± 7.34j, −1.5, −0.515 ± 0.237j.

**(d)** 0.57 deg per metre gives 5.7 deg of pitch for a 10 m step: comfortable for a small UAV at 30 m/s (a climb rate of about Vθ ≈ 3 m/s at the start). A 100 m step would command 57 deg, which is absurd: add a pitch-command limit (for example ±10 deg) and a climb-rate limit, after which the loop is nonlinear and the transfer-function prediction holds only for small altitude errors. The design also ignores the speed loss in a climb (the phugoid and the throttle), which a real altitude hold must manage.

**Common mistakes.** Using ḣ = Vθ (ignoring the γ lag, which hides the −1.5 pole); forgetting the altitude integrator and finding a non-unit DC gain.

## Check-yourself answers (lesson.md)

1. T = G/(1 + GH); closed-loop poles are the roots of 1 + G(s)H(s) = 0, i.e. of d_G d_H + n_G n_H.
2. As a closed-loop zero at s = −100 (the poles of H become zeros of the closed loop).
3. 5/(s² + 2s + 5): poles −1 ± 2j; DC gain 1.
4. T → 1/H: the output follows the inverse of the sensor's transfer function, so sensor errors pass straight to the output.
5. Polynomial multiplication keeps every factor; cancelling is a separate step (`minreal`). It is dangerous when the cancelled pole is unstable (or poorly damped): the mode remains inside the loop and any disturbance excites it.
6. A = [[0, 1], [−5, −3]], B = (0, 1)ᵀ, C = (1, 2), D = 0.
7. It checks the algebra of the reduction (both are the same system if the reduction is right). It does not check the integration error (both share RK4; compare with the exact `scipy.signal.step` for that) and it does not check the model or the drawing.
