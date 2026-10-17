# P09.20  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a scalar MPC by hand)

A scalar discrete plant x_{k+1} = a x_k + b u_k with a = 1 + s/100, b = 0.1, weights q = 1, r = 0.1 and terminal weight P = q.

(a) For horizon N = 1, write J(u₀) = q x₀² + r u₀² + P(a x₀ + b u₀)², minimise it by hand and give K_mpc. (b) Do the same for N = 2 by the Riccati recursion P₁ = q + a²P₂ − (abP₂)²/(r + b²P₂) with P₂ = P, and give K₀. (c) Check (a) and (b) with `unconstrained_mpc_gain` and `finite_horizon_lqr`. (d) Iterate the recursion to convergence (`dare`) and give the infinite-horizon gain; check that `unconstrained_mpc_gain` with P = P_∞ returns it for N = 1.

Answer form: two gains with derivations, two agreeing pairs, one gain with a check.

## 2. (Derivation: why MPC is LQR when nothing is active)

(a) Derive the prediction matrices Sx and Su of X = Sx x₀ + Su U. (b) Show that J = ½UᵀHU + x₀ᵀFᵀU + x₀ᵀYx₀ with H = 2(SuᵀQ̄Su + R̄), F = 2SuᵀQ̄Sx, and that H is positive definite when R > 0. (c) Show by dynamic programming that the first block of the unconstrained minimiser −H⁻¹Fx₀ is −K₀x₀, K₀ being the first gain of the backward Riccati recursion started from P_N = P. (d) Show that if P solves the DARE, the recursion is stationary (P_k = P for all k) and therefore K_mpc equals the infinite-horizon gain for every N. (e) Write the KKT conditions of min ½zᵀHz + fᵀz s.t. Gz ≤ w and explain why they certify the global optimum.

Answer form: five derivations.

## 3. (Unfamiliar scenario: a reaction-wheel slew with a rate limit)

A small satellite (J = 10 kg·m², the `P09.16` model) must slew by (10 + s)° with a wheel torque limit of 0.02 N·m and a body-rate limit of 0.5°/s (a star tracker loses lock above it). Discretise J θ̈ = u by ZOH at 1 s.

(a) Using `prediction_matrices`, `condensed_cost`, `mpc_constraints` (input bound, a loose rate-of-change bound, and the state row c = (0, 1) for the body rate) and `solve_qp`, write a receding-horizon loop with N = 40, Q = diag(1/θ_max², 1/ω_max²), R = 1/u_max² (θ_max = 5°, ω_max = 0.5°/s, u_max = 0.02 N·m) and the DARE terminal cost. (b) Report the peak torque, the peak body rate and the time to come within 0.1° of the target. (c) Compare with the bang-coast-bang lower bound on the slew time: accelerate at 0.002 rad/s² to 0.5°/s, coast, decelerate. (d) How much longer is the MPC slew, and why?

Answer form: a loop (code), three numbers, a lower bound with its arithmetic, one sentence.

## 4. (Numerical, uses your code: your altitude change)

The autopilot of the lesson must climb (20 + s) m.

(a) With N = 30 and the DARE terminal cost, give the peak climb rate, the peak elevator, the peak elevator rate, the overshoot and the time to come within 0.4 m. (b) Compare the time with the lower bound Δh/(2.5 m/s). (c) Give the largest number of QP iterations in any step and the number of active constraints in the first QP. (d) Repeat (a) with the unconstrained LQR clipped to the elevator limits and report its peak climb rate.

Answer form: five numbers, a bound, two counts, one number.

## 5. (Judgment: MPC in a flight computer)

(a) The plant in the lesson equals the prediction model. List three ways a real aircraft differs and say which one threatens the climb-rate constraint most. (b) An updraft adds 1 m/s to the climb rate for 3 s while the aircraft is already climbing at the limit. What happens to the QP, and what is the standard remedy? (c) The QP takes up to 71 iterations at N = 50. Estimate whether it fits in a 10 Hz loop on a flight computer that does 10⁸ floating-point operations per second, counting each iteration as a few solves of size up to N. (d) Would you certify a flight-control MPC by its closed-loop stability proof, by simulation, or by both? What does the proof need that the lesson does not provide?

Answer form: four short paragraphs with numbers where asked.

## 6. (Hard: horizon, feasibility and the terminal cost)

(a) For the 20 m climb with P = DARE, find the shortest horizon N for which the receding-horizon loop stays feasible over 25 s. (b) Explain, with the elevator-rate limit, the climb-rate limit and the short-period dynamics, why short horizons become infeasible a few samples after the start. (c) With P = Q, short horizons remain feasible but do not reach the target in 25 s. Explain the difference between the two failures. (d) For N = 5, 10, 20 report the spectral radius of the unconstrained closed loop with P = Q and with P = DARE, and say what it predicts for (c).

Answer form: one horizon, two explanations, a six-entry table.
