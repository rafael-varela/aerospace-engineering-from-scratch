# P09.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
c = 2D₀/V₀. (a) `nl, dl = open_loop(*p_controller(Kp), *cruise_plant())`; `error_constants`; `steady_state_error(nl, dl, 'step', 2.0)`; drag loss 200/(c + K_p), or `final_value` on 200·`disturbance_to_output`. (b) `critical_gain(*cruise_plant())`. (c) `pi_controller(Kp, Ki)`; `system_type`; `error_constants`; `steady_state_error(..., 'ramp', 0.5)`; `np.roots(closed_loop(nl, dl)[1])`.

## 2.
(a) V = G(C E + D), E = R − V ⇒ E(1 + CG) = R − GD ⇒ E = SR − GSD; T + S = (L + 1)/(1 + L). (b) sE = s·A/(s^{q+1}(1 + K₀N₀/(sⁿD₀))) = A sⁿD₀/(s^q(sⁿD₀ + K₀N₀)); take s → 0 for each pair (n, q). (c) Final value of −G S F/s·(−1)… keep track of the sign: e = r − v, v = GS d, d = −F/s; lim s·GS·F/s = G(0)F/(1 + C(0)G(0)) with G(0) = 1/c. (d) Unity feedback on L = 1/s² (poles at ±j) or on any L with closed-loop poles in the right half plane.

## 3.
L = (K_d s + K_p)/(J s²): count integrators; K_a = K_p/J; closed loop J s² + K_d s + K_p. Parabola error A/K_a. Disturbance: θ/τ_d = G/(1 + CG) = 1/(J s² + K_d s + K_p) ⇒ θ_ss = τ_d/K_p. PID: L = (K_d s² + K_p s + K_i)/(J s³); closed loop J s³ + K_d s² + K_p s + K_i; Routh for the cubic (K_d K_p > J K_i); θ/τ_d = s/(J s³ + …) ⇒ 0. Simulate `simulate_tf(nd, dd, t, 1e-5*np.ones_like(t))` over 3000 s.

## 4.
(a) Loop over K_p; `np.roots(closed_loop(...)[1])`, ζ = −Re p/\|p\|. (b) Loop over K_i. (c) `fvt_vs_simulation(*pi_controller(400, 40), 'ramp', dt=dt, amplitude=0.5)`; compare the difference with A dt/2. (d) Solve 1/(1 + K/c) = 0.0065 for K; compare with `critical_gain`; `final_value` raises, `simulate_tf` diverges.

## 5.
F = m g sin(2°); open loop F/c and m/c. `disturbance_to_output(*ctrl, *cruise_plant())`; simulate with d = −F for 300 s at dt = 0.05; min, argmin, and the last time \|v − v_final\| > 0.1. Steady thrust = F.

## 6.
(a) `sensitivity(nl, dl)`, `tf_eval(ns, ds, 1j*w)` at the six frequencies. (b) m·a plus c·Δv at the end of the ramp. (c) Separate low-frequency results (they depend on c and the controller's DC behaviour) from dynamic ones (which depend on the lags and on the dynamics the model leaves out).
