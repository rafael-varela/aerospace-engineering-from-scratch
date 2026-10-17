# P09.20  Method

The route to each answer, without the final numbers. Pair with `hints.md`; `worked-solutions.md` has the numbers for s = 15.

## 1. Scalar MPC

(a) dJ/du₀ = 2ru₀ + 2Pb(ax₀ + bu₀) = 0 ⇒ u₀ = −abP/(r + b²P)·x₀, so K_mpc = abP/(r + b²P) with P = q. (b) P₁ = q + a²P₂ − (abP₂)²/(r + b²P₂), K₀ = abP₁/(r + b²P₁). (c) 1 × 1 arrays into `unconstrained_mpc_gain(A, B, Q, R, Q, N)` and `finite_horizon_lqr(A, B, Q, R, Q, N)[0][0]`. (d) `dare`, `lqr_gain`, then `unconstrained_mpc_gain(..., P_inf, 1)`; the closed-loop pole a − bK_∞.

## 2. Derivations

(a) x_k = Aᵏx₀ + Σ_{j=0}^{k−1}A^{k−1−j}Bu_j; stack. (b) Σ_{k=1}^{N} x_kᵀQ̄_kx_k = (Sx x₀ + Su U)ᵀQ̄(Sx x₀ + Su U); expand; add UᵀR̄U and x₀ᵀQx₀; H = 2(SuᵀQ̄Su + R̄) is positive definite because UᵀR̄U > 0 for U ≠ 0. (c) Dynamic programming: V_N = xᵀPx; V_k(x) = min_u[xᵀQx + uᵀRu + V_{k+1}(Ax + Bu)] gives (5) and u_k = −K_kx_k; the minimiser of a strictly convex function is unique, so the batch first move equals −K₀x₀. (d) P_N = P_DARE ⇒ P_{N−1} = P_DARE by the DARE ⇒ induction; all K_k = K_∞. (e) Stationarity, primal and dual feasibility, complementarity; for convex objective and affine constraints, KKT points are global minimisers (and here unique since H > 0).

## 3. Reaction-wheel slew

(a) Build A_d, B_d by `c2d_zoh`, Q and R by Bryson, P by `dare`, H, F once and (Sx, Su) once; at each step `mpc_constraints(Ad, Bd, 40, x, u_prev, 0.02, 1.0, c=(0, 1), c_max=0.5°/s, pred=...)`, `solve_qp`, apply U[0]. (b) Track max |u|, max |ω| and the last time |θ − θ_target| > 0.1°. (c) a_max = u_max/J; t_a = ω_max/a_max; angle during acceleration ½a_max t_a²; coast the rest at ω_max; T = 2t_a + (Δθ − 2·½a_max t_a²)/ω_max. (d) Compare and relate to the quadratic cost's exponential final approach.

## 4. Your climb

`simulate_mpc` from (−(20 + s), 0, 0, 0) for 40 s; peaks from x[:, 1]·25, u and diff(u)/0.1; overshoot max(0, max Δh); `settling_time`; the bound Δh/2.5; `iters.max()`; for the first QP count λ > 10⁻⁹. Then `simulate_lqr(..., clip=altitude_limits())`.

## 5. Judgment

(a) Gusts and updrafts, speed and thrust changes, model errors in L_α and M_δ, sensor noise and delay; the updraft acts directly on the climb rate. (b) Hard constraint violated at the current state ⇒ QP infeasible; soften with slack variables, or drop/relax the state constraint for a few steps. (c) Estimate flops per iteration (forming and solving the active-set system, matrix-vector products with G) times the iteration count, against 10⁷ flops available per 0.1 s; consider incremental factorisations. (d) Both: the proof needs a terminal set (recursive feasibility) and exact model assumptions; simulation and flight test check the rest.

## 6. Horizon study

(a) Increase N until `simulate_mpc(x0, 250, N, "dare")` raises no error. (b) Count the samples needed to bring γ from its limit to zero with the elevator moving at most 3° per sample through the short-period dynamics; if the horizon is shorter, the controller cannot see the stop. (c) P = Q: feasible but weak (timid) gain; P = DARE: aggressive but short-sighted. (d) Spectral radii of A − BK_mpc.
