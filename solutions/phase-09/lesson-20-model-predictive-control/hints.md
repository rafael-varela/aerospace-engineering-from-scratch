# P09.20  Hints

One hint per problem. Read only the one you need.

1. J(u₀) is a quadratic in one variable: set dJ/du₀ = 0. For N = 2, run the recursion once backwards from P₂ = q to get P₁, then use P₁ in the gain formula.

2. (a) Unroll x_k = Ax_{k−1} + Bu_{k−1}. (b) Substitute X into Σ xᵀQx and collect the terms quadratic, linear and constant in U; note the k = 0 term x₀ᵀQx₀ does not involve U. (c) Both methods minimise the same strictly convex function; dynamic programming returns the optimal first move as a feedback of x₀. (d) If P_{k+1} = P solves the DARE, what is P_k? (e) A convex problem with linear constraints: KKT is necessary and sufficient.

3. ZOH of J θ̈ = u at T = 1 s: A_d = [[1, 1], [0, 1]], B_d = [[0.5/J], [1/J]]. The state constraint row is c = (0, 1) with limit 0.5°/s in rad/s. Pass a large rate-of-change bound (for example 1 N·m) so it never binds. For the bound, the angle covered while accelerating is ½a t_a².

4. `simulate_mpc(np.array([-(20 + s), 0, 0, 0]), 400, N=30)`; settling with `settling_time`; for the first QP build `condensed_cost`, `mpc_constraints` and count multipliers above 10⁻⁹.

5. Think about what the model cannot predict (gusts, speed, unmodelled dynamics). Ask whether a hard constraint can be met if the state is already outside it. Count flops of forming and solving N̄H⁻¹N̄ᵀ for up to N active rows.

6. Loop N upward from 3, catching `ValueError` from `simulate_mpc`. For the spectral radius use `unconstrained_mpc_gain` with P = Q and with P = `dare(...)`.
