# P09.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
K = 10 + 2s. (a) GM = 120/K, ω_pc = 1/√(0.05 × 0.01). (b) `margins(*servo_loop(K))`. (c) `nyquist_count(*servo_loop(K))`; M_s from `np.max(1/np.abs(1 + tf_eval(num, den, 1j*w)))` on a fine log grid. (d) Compare 20 log₁₀ GM with 6 and PM with 45.

## 2.
(a) Argument principle on F = 1 + L around the clockwise D-contour; F's zeros = closed-loop poles, F's poles = open-loop poles; N_cw(F about 0) = Z − P. (b) L ≈ (c/ε^k)e^{−jkφ} as φ goes −90° → +90°: angle goes from +k·90° to −k·90°. (c) L(jω) = K/(jω(1 + jτ_mω)(1 + jτ_fω)); the denominator is jω(1 − τ_mτ_fω²) − (τ_m + τ_f)ω²; it is real (L real and negative) when 1 − τ_mτ_fω² = 0; then \|L\| = K/((τ_m + τ_f)ω²) = Kτ_mτ_f/(τ_m + τ_f). Compare with Routh for τ_mτ_f s³ + (τ_m + τ_f)s² + s + K. (d) e^{−jωT} has unit magnitude and phase −ωT; at each gain crossover the phase available is PM_i; the first to be used up sets DM = min PM_i/ω_gc,i.

## 3.
`margins(*unstable_airframe_loop(20.0))`; for f in (0.5, 0.9, 1.1, 1.5): `nyquist_count(num, den, delay=f*DM)` and `closed_loop_rhp_roots(num, den, f*DM, pade_order=8)`. For (d): L(0) = −3 and the closed-loop polynomial's constant term 3K − 20.

## 4.
Loop over K: `margins`, M_s on a log grid. For (c) bisect on K for GM_dB = 6 (K = 120/10^{0.3}) and for PM = 45° separately; the smaller K binds.

## 5.
`margins(*satellite_loop(kp=kp))`. (b) d/dω[arctan(T_dω) − arctan(τω)] = T_d/(1 + T_d²ω²) − τ/(1 + τ²ω²) = 0 ⇒ ω² = 1/(T_dτ); evaluate. (c) PM_with_delay = PM − ω_gc × 0.5 s (in degrees); keep K_p with ≥ 30°; confirm with `phase_margin(..., delay=0.5)`.

## 6.
(a) `servo_loop(20, tau_f=0.02)`; K_crit = (τ_m + τ_f)/(τ_mτ_f). (b) `margins(*servo_loop(20), delay=0.02)`. (c) M_s is the inverse of the smallest distance from L(jω) to −1 over all ω. (d) Think of saturation, rate limits, backlash, friction.
