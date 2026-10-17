# P09.20  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (the functions named in each problem) or checked with `scipy.linalg.solve_discrete_are`; the "independent check" line says how. Units: SI; elevator in degrees where stated.

## 1. Scalar MPC

**Assumptions.** a = 1 + 15/100 = **1.15** (unstable), b = 0.1, q = 1, r = 0.1, P = q = 1.

**(a)** dJ/du₀ = 2ru₀ + 2Pb(ax₀ + bu₀) = 0 ⇒ **K_mpc = abP/(r + b²P) = 0.115/0.11 = 1.04545**.

**(b)** P₁ = 1 + 1.3225 − 0.115²/0.11 = **2.20227**; **K₀ = abP₁/(r + b²P₁) = 0.253261/0.122023 = 2.07553**.

**(c)** `unconstrained_mpc_gain` gives 1.04545 (N = 1) and 2.07553 (N = 2); `finite_horizon_lqr(…, N = 2)` gives [K₀, K₁] = [2.07553, 1.04545]: the last gain of the two-step problem is the one-step gain, as dynamic programming says.

**(d)** `dare` gives **P_∞ = 5.91548**, **K_∞ = 4.27433**, and `unconstrained_mpc_gain(…, P_∞, 1)` returns 4.27433. Closed-loop pole a − bK_∞ = 0.72257.

**Common mistakes.** Using P₂ instead of P₁ in K₀; forgetting the terminal term in J for N = 1 (then K = 0).

## 2. Derivations

**(a)** x_k = Aᵏx₀ + Σ_{j=0}^{k−1}A^{k−1−j}Bu_j: block row k of **Sx is Aᵏ**, block (k, j) of **Su is A^{k−1−j}B for j < k and 0 otherwise** (block lower-triangular: causality).

**(b)** J = x₀ᵀQx₀ + XᵀQ̄X + UᵀR̄U with X = Sx x₀ + Su U: XᵀQ̄X = UᵀSuᵀQ̄SuU + 2x₀ᵀSxᵀQ̄SuU + x₀ᵀSxᵀQ̄Sxx₀. So **J = ½UᵀHU + x₀ᵀFᵀU + x₀ᵀYx₀** with H = 2(SuᵀQ̄Su + R̄), F = 2SuᵀQ̄Sx, Y = Q + SxᵀQ̄Sx. UᵀHU ≥ 2UᵀR̄U > 0 for U ≠ 0: H ≻ 0.

**(c)** Define V_k(x) = min over u_k…u_{N−1} of the remaining cost, V_N = xᵀPx. If V_{k+1} = xᵀP_{k+1}x, then V_k(x) = min_u[xᵀQx + uᵀRu + (Ax + Bu)ᵀP_{k+1}(Ax + Bu)], minimised by **u = −K_kx, K_k = (R + BᵀP_{k+1}B)⁻¹BᵀP_{k+1}A**, with V_k = xᵀP_kx and P_k = Q + AᵀP_{k+1}(A − BK_k). Both the batch problem (4) and dynamic programming minimise the same strictly convex J, whose minimiser is unique, so the first move of U* is u₀ = −K₀x₀: **K_mpc = K₀**.

**(d)** If P_{k+1} = P_DARE then P_k = Q + AᵀPA − AᵀPB(R + BᵀPB)⁻¹BᵀPA = P_DARE, so by induction every P_k equals P_DARE and every K_k = K_∞; hence **K_mpc = K_∞ for all N**.

**(e)** HU + f + Gᵀλ = 0, GU ≤ w, λ ≥ 0, λ_i(GU − w)_i = 0. For any feasible V, J(V) − J(U) = ½(V − U)ᵀH(V − U) + (HU + f)ᵀ(V − U) ≥ −λᵀG(V − U) = −λᵀ(GV − w) + λᵀ(GU − w) = −λᵀ(GV − w) ≥ 0. So U is the **global** minimum, unique because H ≻ 0.

## 3. Reaction-wheel slew

**Assumptions.** J = 10 kg·m², Δθ = 10 + 15 = **25°**, T = 1 s, N = 40, torque ≤ 0.02 N·m, |ω| ≤ 0.5°/s, Bryson weights, DARE terminal cost; the model is the plant.

**(a)** The loop is as in `method.md`; it solved every QP (at most 64 iterations).

**(b)** Peak torque **0.0200 N·m**, peak body rate **0.500°/s** (the limit, reached and held), within 0.1° of the target after **79 s**.

**(c)** a_max = 0.02/10 = 0.002 rad/s²; t_a = 0.008727/0.002 = **4.363 s**; each of acceleration and deceleration covers ½·0.002·4.363² = 0.01904 rad = **1.091°**; coast (25 − 2.182)/0.5 = 45.64 s; **T = 54.36 s**.

**(d)** The MPC slew takes 79 s, **24.6 s (45 %) longer** than bang-coast-bang. The quadratic cost does not ask for minimum time: near the target the constraints are inactive and the closed loop approaches with the exponential tail of its LQR poles, and the deceleration starts earlier than a time-optimal law would start it. A minimum-time objective (or a larger weight on θ) shortens it.

## 4. Your climb (s = 15: 35 m)

**(a)** Peak climb rate **2.5000 m/s** (above the limit by 4.4 × 10⁻¹⁰ m/s, the QP tolerance), peak elevator **10.00°**, peak elevator rate **30.00°/s**, overshoot **6.1 mm**, within 0.4 m after **15.4 s**.

**(b)** Lower bound 35/2.5 = **14 s**; the MPC needs 1.4 s more, the time to start and stop the climb through the short-period dynamics and the elevator-rate limit (the same 1.4 s as for 20 m: 9.4 against 8 s).

**(c)** At most **78** QP iterations in any step; the first QP ends with **29** active constraints after 59 iterations.

**(d)** Clipped LQR peaks at **12.43 m/s** of climb, five times the limit.

## 5. Judgment

**(a)** Updrafts and gusts (they change the climb rate directly, and the model has no disturbance input), speed and thrust changes (the constant-speed assumption), errors in L_α and M_δ, sensor noise and delay. The updraft threatens the climb-rate constraint most: it acts on exactly the constrained quantity without passing through the elevator.

**(b)** If the aircraft is already climbing at 2.5 m/s and an updraft adds 1 m/s, the measured state violates the constraint at k = 0 or forces a violation at k = 1 whatever the elevator does: the QP is infeasible and `solve_qp` raises. The standard remedy is to soften the state constraint with a slack variable s ≥ 0, |climb| ≤ 2.5 + s, heavily penalised in the cost (Extension 1), so the QP stays feasible and returns the least-violating input; hard limits are kept only for physical actuator bounds.

**(c)** Estimate for N = 50 with our dense implementation: per iteration, forming N̄H⁻¹N̄ᵀ with up to 50 active rows costs about 2 × 50³ = 2.5 × 10⁵ flops, its solve about 4 × 10⁴, and the products with the 300 × 50 matrix G about 3 × 10⁴: roughly 3 × 10⁵ flops. Times 71 iterations ≈ 2 × 10⁷ flops, against 10⁷ available in 0.1 s at 10⁸ flop/s: **too slow by about a factor 2** (an estimate). Incremental factorisation updates (O(n²) per iteration), warm starts and a shorter horizon (N = 20 needs at most 33 iterations) bring it well inside the budget.

**(d)** Both. The classical stability proof for constrained MPC needs a terminal cost (we have it), a terminal constraint set that is invariant under the terminal LQR law and inside the constraints (we do not construct it), and an exact model; without the terminal set, recursive feasibility is not guaranteed (Problem 6 shows it failing). Simulation over the flight envelope, with disturbances and model errors, and flight test then check what the proof assumes.

## 6. Horizon study

**(a)** With P = DARE the loop is infeasible at N = 3 (step 8), N = 4 (step 7) and N = 5 (step 7); **N = 6** is the shortest horizon feasible over 25 s.

**(b)** Early in the climb the controller pitches up hard (the DARE cost-to-go promises an aggressive LQR afterwards). To stop at 2.5 m/s the flight-path angle must be levelled off, which takes the elevator from up to down through the short-period response, and the elevator can move only 3° per sample. If the horizon (0.3–0.5 s) is shorter than that manoeuvre, the plan chosen at one step leaves the next step with a state from which every input sequence over the horizon overshoots the climb-rate limit: the QP is infeasible.

**(c)** With P = Q there is no promise beyond the horizon, the unconstrained gain is weak, the controller climbs slowly and never approaches the limit, so it stays feasible, but it reaches only −8.21 m (N = 3) or −1.49 m (N = 5) of error after 25 s. One failure is too timid (performance), the other too bold (feasibility).

**(d)**

| N | spectral radius, P = Q | spectral radius, P = DARE |
|---|---|---|
| 5 | 0.9886 | 0.8905 |
| 10 | 0.9504 | 0.8905 |
| 20 | 0.9006 | 0.8905 |

With P = Q the slowest mode of the unconstrained loop decays by only 1.1 % per sample at N = 5 (a time constant of about 9 s), which predicts the slow creep in (c); with the DARE cost the radius is that of K_∞ for every N.

## Check-yourself answers (lesson.md)

1. X = Sx x₀ + Su U with Sx = [A; …; A^N] and (Su)_{k,j} = A^{k−1−j}B for j < k; block lower-triangular because x_k cannot depend on inputs applied at or after step k.
2. H = 2(SuᵀQ̄Su + R̄), F = 2SuᵀQ̄Sx; H ≻ 0 because R̄ ≻ 0 adds a positive definite term.
3. u₀ = −K_mpc x₀ with K_mpc the first m rows of H⁻¹F; it equals K₀ of the finite-horizon Riccati recursion started at P_N = P.
4. The DARE solution is the fixed point of the Riccati recursion, so every stage gain equals K_∞; the terminal cost is exactly the cost of continuing with LQR.
5. Rate: ±(u_k − u_{k−1}) ≤ Δu_max with u_{−1} from the last applied input (depends on the previous input). State: ±C̄Su U ≤ x_max ∓ C̄Sx x₀ (depends on the current state).
6. HU + f + Gᵀλ = 0, GU ≤ w, λ ≥ 0, λ_i(GU − w)_i = 0; the problem is convex with affine constraints, so KKT points are global minima (unique since H ≻ 0).
7. Clipping acts on the current command only and cannot see where the state is going; the climb rate responds to the elevator through the dynamics, so by the time it reaches the limit the clipped command can no longer stop it.
8. That a feasible QP at one step implies a feasible QP at the next. At N = 3 with the DARE terminal cost the QP became infeasible 0.8 s into the climb.
