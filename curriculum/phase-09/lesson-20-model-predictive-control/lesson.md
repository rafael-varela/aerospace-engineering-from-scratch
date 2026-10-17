# P09.20  Model Predictive Control ★

> LQR finds the best linear law and then ignores the limits. Model predictive control solves a small optimisation problem at every sample, with the limits written in, and applies only the first move. When no limit is active, it is LQR.

| | |
|---|---|
| **Learning outcomes** | 1. Build the prediction matrices X = Sx x₀ + Su U of a discrete linear model and write the finite-horizon quadratic cost as a condensed quadratic program (QP) in the input sequence. 2. Prove and verify that the unconstrained MPC law equals the first gain of the finite-horizon LQR from the backward Riccati recursion and, with the terminal cost P from the discrete algebraic Riccati equation (DARE), the infinite-horizon LQR gain for every horizon. 3. Write input, input-rate and state constraints as linear inequalities, solve the QP from scratch with a dual active-set (Goldfarb–Idnani) method, and certify the solution by its KKT conditions. 4. Run the receding horizon on an altitude-hold autopilot with elevator, elevator-rate and climb-rate limits, show that no constraint is violated, and explain the effect of the horizon and the terminal cost on performance and feasibility. |
| **Prerequisites** | `P09.16` (LQR and the Riccati equation), `P09.17` (ZOH discretisation), `P09.13` (state space), `P09.11` (pitch loop, actuator limits), `P03.15` (constrained optimisation: Lagrangian, KKT conditions, active sets, QP), `P03.14` (unconstrained optimisation). Skill check: you can write the KKT conditions of min ½zᵀHz + fᵀz subject to Gz ≤ w and solve a 2-variable equality-constrained QP by hand. |
| **Engineering question** | A small UAV's altitude-hold autopilot commands the elevator directly. It must climb 20 m without exceeding ±10° of elevator, 30°/s of elevator rate or 2.5 m/s of climb rate. An LQR tuned for a crisp response asks for 77° of elevator and 9.8 m/s of climb; clipping its command keeps the elevator legal but still climbs at 8.2 m/s. Can an MPC respect all three limits at every sample, how fast does it reach the new altitude, how long a horizon does it need, and how do we know its optimisation is solved correctly? |
| **Assumptions and limits** | A linear, time-invariant longitudinal model at constant speed (25 m/s, an autothrottle holds speed), with the representative short-period derivatives of `P09.17` and altitude as the integral of V γ; ZOH at 10 Hz; the plant equals the prediction model (no disturbance, no model error); full state measured. Constraints are hard; the QP is convex (H positive definite) and solved to a tolerance of 10⁻¹⁰ relative to the constraint scale. Stability and recursive feasibility of constrained MPC need a terminal set as well as a terminal cost, which this lesson does not construct; it demonstrates the loss of feasibility at short horizons instead. |
| **Analytical work** | The prediction matrices; the condensed cost H, F, Y; the unconstrained minimiser and its equality with the Riccati recursion (dynamic programming); the stationarity of the recursion at the DARE solution; the constraint matrices for amplitude, rate and state limits; the KKT conditions and the Goldfarb–Idnani step (primal direction, dual direction, full and partial steps); the lower bound on climb time Δh/ḣ_max. |
| **Computational deliverable** | `starter/main.py`: `finite_horizon_lqr`, `dare`, `lqr_gain`, `prediction_matrices`, `condensed_cost`, `unconstrained_mpc_gain`, `mpc_constraints`, `solve_qp`, `kkt_residuals`, `mpc_step` (`expm`, `c2d_zoh`, `altitude_model`, `altitude_limits`, `altitude_weights`, `simulate_mpc`, `simulate_lqr`, `settling_time` and `main` are given). Tests prove: the ZOH model and its two integrators; the DARE by iteration matches `scipy.linalg.solve_discrete_are` with residual below 10⁻⁸ and a stable closed loop; the Riccati recursion converges to the DARE; the prediction matrices reproduce a direct simulation and are causal; the condensed cost equals the stage-wise sum to 10⁻¹⁰ and H > 0; **the unconstrained MPC gain equals the finite-horizon Riccati gain for N = 1, 3, 10, 40 and, with P from the DARE, the infinite-horizon gain for N = 1, 5, 25**; the QP solver reproduces a hand solution with its multipliers, satisfies the KKT conditions on five random problems and beats every feasible perturbation, returns the unconstrained minimum when nothing is active, and detects infeasibility; **the receding-horizon loop never violates the elevator, elevator-rate or climb-rate limit** and uses them; it reaches the target without overshoot; LQR would break the limits; with loose limits MPC is the LQR loop; a three-step horizon loses feasibility and says so. |
| **Practical activity** | **computational experiment**: unconstrained LQR, clipped LQR and MPC on the 20 m climb; a horizon sweep N = 3 to 50 with and without the DARE terminal cost (performance, feasibility, QP iterations); a KKT optimality check of the first QP. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `mpc_altitude.png`, the derivations of the condensed cost and of the MPC–LQR equivalence, the horizon table, and a half-page note choosing a horizon and terminal cost for the autopilot with its constraint record, climb time and worst-case QP iterations. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the condensed QP and the equivalence with the Riccati recursion are derived; all public tests pass; the constraint record of the closed loop is shown; the note explains the loss of feasibility at short horizons and states what a flight implementation would add. |
| **References** | `@mit-ocw-16-323-principles-of-optimal-control` Lecture 16 "Model predictive control" (PDF heading "Model Predictive Control") (https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec16/) (required-free; graduate): stacking the predictions over a horizon into a standard QP, the full-state assumption, stability with constraints, and a helicopter example (§2.2, §2.4, §2.6). Lecture 3 "Dynamic programming: principle of optimality, dynamic programming, discrete LQR" (.../lec3/) (required-free): the backward Riccati recursion that the unconstrained MPC must match (§2.3). `@rawlings-mayne-diehl-mpc-2nd-edition` (https://sites.engineering.ucsb.edu/~jbraw/mpc/) (supplementary-free; free PDF from the authors but all rights reserved: link only, nothing copied or adapted): Ch. 1 "Getting Started with Model Predictive Control", 1.3 Introductory MPC Regulator (1.3.1–1.3.6, the linear quadratic problem, its dynamic-programming solution and the infinite-horizon LQ problem); Ch. 2 "Model Predictive Control—Regulation", 2.2 Model Predictive Control, 2.4 Stability, 2.6 Is a Terminal Constraint Set Xf Necessary?; Ch. 8 "Numerical Optimal Control", 8.7 Newton-Type Optimization with Inequalities, 8.8 Structure in Discrete Time Optimal Control: the receding horizon, terminal cost and terminal set at graduate depth. `@ocw-18-335j` Week 11 "Lecture 29: Lagrange Dual Problems" (supplementary-free): the Lagrangian, KKT conditions and duality behind the active-set QP (as used in `P03.15`). `@scipy-userguide` API reference "Linear algebra (scipy.linalg)": Matrix Equation Solvers (`solve_discrete_are`) (supplementary-free): the oracle of the DARE tests. Partial GAP (original instruction): no verified free source at this level walks through a from-scratch active-set QP solver for the constrained problem; the Goldfarb–Idnani dual active-set solver (§2.5) and the altitude-hold MPC example are written here. |
| **Extensions** | ★ Soften the climb-rate constraint with a slack variable and a large linear and quadratic penalty, and show that the updraft of Problem 5 no longer makes the QP infeasible. ★ Warm-start the active-set solver with the previous solution shifted by one step and measure the reduction in iterations. ★ Compute the maximal output-admissible set for the LQR terminal law and use it as a terminal constraint; check recursive feasibility at N = 3. |
| **Time** | unmeasured |

## 1. The engineering problem

Every autopilot lives inside limits. The elevator moves ±10° and no faster than its servo allows; the aircraft must not climb or descend faster than the airspace, the passengers or the engine allow; a quadrotor must not tilt beyond the angle at which its thrust can hold altitude. Linear controllers (PID, `P09.11`; LQR, `P09.16`) do not know these limits. They are tuned so that, for the commands expected, the limits are rarely reached, and when they are reached, the command is clipped and anti-windup logic patches the consequences. That works when limits are an exception. When the limits *are* the problem, as in a large altitude change where the best possible climb is the fastest legal one, clipping is the wrong tool: it can enforce actuator limits but not state limits, because it does not know where the state is going.

**Model predictive control** (MPC) looks ahead. At every sample it uses the model to predict the next N states as a function of the next N inputs, chooses the input sequence that minimises a quadratic cost subject to every limit over the whole horizon, applies the first input, and repeats at the next sample from the new measured state (the **receding horizon**). The limits are part of the problem, so the controller slows the climb *before* the climb-rate limit is reached and stops it in time.

The price is an optimisation problem solved in real time. For a linear model, quadratic cost and linear limits it is a **quadratic program** (QP), convex and solvable exactly in a finite number of steps. This lesson builds the whole chain from scratch for the altitude hold of a small UAV and checks it against the theory it must reduce to: when no limit is active, MPC must *be* LQR.

Success criteria, written first: the unconstrained MPC gain equals the Riccati recursion's gain to 10⁻⁹; with the DARE terminal cost it equals the infinite-horizon gain for every N; the QP solutions satisfy the KKT conditions to 10⁻⁹; in closed loop the elevator, its rate and the climb rate never exceed their limits by more than the solver tolerance (10⁻⁹); the 20 m climb completes without overshoot beyond a few millimetres.

## 2. Understand and derive

### 2.1 The model and the limits

The longitudinal model at constant speed V = 25 m/s, with state x = (Δh, γ, θ, q): altitude error Δh = h − h_ref, flight-path angle γ, pitch angle θ, pitch rate q, and input u = δ_e (elevator, positive nose-up):

Δḣ = Vγ,  γ̇ = L_α(θ − γ),  θ̇ = q,  q̇ = M_α(θ − γ) + M_q q + M_δ δ_e,   (1)

with α = θ − γ and the representative values L_α = 1.0 s⁻¹, M_α = −9 s⁻², M_q = −2.5 s⁻¹, M_δ = 12 s⁻². It has two eigenvalues at 0 (altitude integrates climb rate; and a trimmed climb at γ = θ needs no elevator in this constant-speed model, because the autothrottle supplies the energy) and the short-period pair −1.75 ± 2.905j s⁻¹. It is discretised exactly by ZOH at T = 0.1 s (`P09.17`): x_{k+1} = Ax_k + Bu_k.

The limits: |δ_e| ≤ 10°; |δ_{e,k} − δ_{e,k−1}| ≤ 3° per sample (30°/s); |Vγ_k| ≤ 2.5 m/s. The weights come from Bryson-like scales (2 m, 0.1 rad, 0.2 rad, 0.5 rad/s; 0.17 rad of elevator).

### 2.2 Prediction and the condensed cost

Over a horizon of N samples, x_k = Aᵏx₀ + Σ_{j<k} A^{k−1−j}Bu_j. Stacking X = (x₁, …, x_N) and U = (u₀, …, u_{N−1}):

**X = Sx x₀ + Su U**, Sx = [A; A²; …; A^N], (Su)_{k,j} = A^{k−1−j}B for j < k, 0 otherwise.   (2)

Su is block lower-triangular: an input cannot affect earlier states (causality). The cost

J = Σ_{k=0}^{N−1}(x_kᵀQx_k + u_kᵀRu_k) + x_NᵀPx_N   (3)

becomes, with Q̄ = blkdiag(Q, …, Q, P) and R̄ = blkdiag(R, …, R),

**J = ½UᵀHU + x₀ᵀFᵀU + x₀ᵀYx₀,  H = 2(SuᵀQ̄Su + R̄),  F = 2SuᵀQ̄Sx.**   (4)

H is positive definite whenever R > 0 (R̄ alone already is). Without constraints, ∇J = HU + Fx₀ = 0 gives U* = −H⁻¹Fx₀, and the controller applies the first block: **u₀ = −K_mpc x₀, K_mpc = (H⁻¹F)_{first m rows}**. The unconstrained MPC is a linear state feedback. (16.323 Lecture 16 stacks the predictions into the same standard QP.)

### 2.3 Why the unconstrained MPC is LQR

The same minimisation can be done backwards in time by dynamic programming. The cost-to-go from step k with N − k steps left is quadratic, V_k(x) = xᵀP_kx, with P_N = P, and minimising stage k gives

K_k = (R + BᵀP_{k+1}B)⁻¹BᵀP_{k+1}A,  P_k = Q + AᵀP_{k+1}(A − BK_k).   (5)

This is the **backward Riccati recursion** (the discrete counterpart of `P09.16`; 16.323 Lecture 3 derives it by dynamic programming). Since dynamic programming and the batch least-squares (4) minimise the same convex function, their first moves agree: **K_mpc = K₀**. The test checks it for N = 1, 3, 10 and 40 to 10⁻⁹.

If the terminal weight P is the fixed point of (5), the **DARE** P = Q + AᵀPA − AᵀPB(R + BᵀPB)⁻¹BᵀPA, then every P_k equals P and every K_k equals the infinite-horizon LQR gain K_∞, so **K_mpc = K_∞ for every horizon N**, even N = 1. The terminal cost is the exact cost of continuing with LQR after the horizon, so the finite horizon loses nothing. `dare` solves it by iterating (5) from P = Q to convergence; it matches `scipy.linalg.solve_discrete_are` to 10⁻⁸. For the altitude model, K_∞ = [0.0673, 2.1851, 0.9248, 0.3195] with closed-loop eigenvalue magnitudes 0.746 (twice) and 0.891 (twice); with P = Q instead, N = 5 gives a much weaker gain [0.0026, 0.5221, 0.1985, 0.1848] and N = 20 approaches K_∞ ([0.0558, 1.9418, 0.8355, 0.3049]).

### 2.4 Constraints as linear inequalities

For a single input, write the limits over the horizon as **GU ≤ w**:

- amplitude: U ≤ u_max·1 and −U ≤ u_max·1;
- rate: DU ≤ Δu_max·1 + e and −DU ≤ Δu_max·1 − e, where (DU)_k = u_k − u_{k−1} and e = (u_{−1}, 0, …, 0) carries the input applied at the previous sample;
- state: with c = (0, V, 0, 0) picking the climb rate, C̄ = blkdiag(c, …, c): C̄Su U ≤ ḣ_max − C̄Sx x₀ and −C̄Su U ≤ ḣ_max + C̄Sx x₀.

The first two depend only on the previous input; the third depends on the current state. For N = 30 there are 180 inequalities in 30 unknowns.

### 2.5 Solving the QP: KKT conditions and the dual active set

The QP min ½UᵀHU + fᵀU s.t. GU ≤ w (f = Fx₀) is convex, so a point U with multipliers λ is optimal **if and only if** the KKT conditions hold (`P03.15`):

HU + f + Gᵀλ = 0 (stationarity), GU ≤ w (primal feasibility), λ ≥ 0 (dual feasibility), λ_i(GU − w)_i = 0 (complementarity).   (6)

`kkt_residuals` evaluates all four; they are the solver's proof of optimality.

`solve_qp` is the **dual active-set method of Goldfarb and Idnani** (1983). It starts at the unconstrained minimum −H⁻¹f, which satisfies everything in (6) except primal feasibility, and keeps (6) except primal feasibility true at every step. At each outer step it picks the most violated constraint p and moves along

z = −H⁻¹(g_p + N̄ᵀr),  r = −(N̄H⁻¹N̄ᵀ)⁻¹N̄H⁻¹g_p,   (7)

where N̄ holds the rows of G in the current active set. Along z the active constraints stay at equality (N̄z = 0), λ_p grows by the step t and the active multipliers change by tr. The step is the smaller of the **full step** that makes constraint p active and the **partial step** at which an active multiplier reaches zero (that constraint is dropped and the step repeated). The objective rises monotonically, so no active set repeats, and the method ends in finitely many steps, either with no violated constraint (optimal) or with no possible step (infeasible: no U satisfies all constraints). For the first QP of the 20 m climb (N = 30, DARE terminal cost) it takes 41 iterations, ends with 27 active constraints, and has KKT residuals of 1.0 × 10⁻¹¹ (stationarity), 4.3 × 10⁻¹⁵ (primal), 0 (dual) and 3.0 × 10⁻¹² (complementarity).

Why a dual method rather than projected gradient? Projected gradient handles box constraints on U cheaply, but the climb-rate limit is a general linear inequality in U, and projected gradient converges only to a tolerance; the active-set method ends at an exact vertex of the active set, which is what "never violated" needs.

### 2.6 The receding horizon, feasibility and the horizon length

Applying u₀ and re-solving from the measured x₁ turns the open-loop plan into feedback. Two properties that LQR has for free are not automatic:

- **Recursive feasibility.** The QP at the next step may have no solution. With hard limits, a short horizon can let the state reach a point from which the climb can no longer be stopped within the elevator-rate limit (part 4: N = 3, 4 and 5 fail 0.7–0.8 s into the climb).
- **Stability.** With P from the DARE and a horizon long enough that the constraints are inactive at its end, the closed loop inherits LQR's stability; a rigorous guarantee also needs a terminal constraint set (Extension 3; Rawlings, Mayne and Diehl §2.4 and §2.6).

### 2.7 Where the model stops being valid

- **Model equals plant.** Here the prediction is exact. A real aircraft has unmodelled dynamics, speed changes and gusts; the climb-rate limit can then be exceeded between predictions or the QP can become infeasible (Problem 5). Flight implementations soften state constraints with slack variables (Extension 1).
- **Constant speed.** The autothrottle assumption hides the energy cost of climbing; at the limit of engine power, the achievable climb rate itself depends on speed.
- **Computation.** The QP must finish within the 0.1 s sample; iteration counts (41 to 78 here) bound the work but not the timing jitter (`P09.17`).
- **Hard limits only on the samples.** The constraints hold at t = kT; between samples the continuous climb rate can exceed the limit slightly.

## 3. Build the computational model

Open `starter/main.py`. The matrix exponential, ZOH, the model, the limits and weights, the two simulations and `main` are given.

### Step 1: Riccati

`finite_horizon_lqr` runs (5) backwards from P_N, storing K_k in a list indexed from 0. `dare` iterates it from P = Q until the relative change is below 10⁻¹² (raise `RuntimeError` otherwise). `lqr_gain` is one `np.linalg.solve`.

### Step 2: prediction and cost

```python
def prediction_matrices(Ad, Bd, N):
    n, m = Bd.shape
    powers = [np.eye(n)]
    for _ in range(N):
        powers.append(Ad @ powers[-1])
    Sx = np.vstack(powers[1:])
    Su = np.zeros((N * n, N * m))
    for k in range(1, N + 1):
        for j in range(k):
            Su[(k - 1) * n:k * n, j * m:(j + 1) * m] = powers[k - 1 - j] @ Bd
    return Sx, Su
```

`condensed_cost` builds Q̄ (the last block P) and R̄ = I_N ⊗ R and returns (H, F, Y) of (4); symmetrise H. `unconstrained_mpc_gain` is the first m rows of `np.linalg.solve(H, F)`.

### Step 3: constraints and QP

`mpc_constraints` stacks the six blocks of 2.4 in the order +u, −u, +Δu, −Δu, +c, −c (D = I − shifted identity). `solve_qp` implements (7) with an explicit H⁻¹ (fine for N ≤ 50), the full and partial steps, and the drop; it raises `ValueError` when no step is possible. `kkt_residuals` returns the four residuals of (6). `mpc_step` assembles H, F (cached between calls) and G, w, solves, and returns the first move.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-20-model-predictive-control`. The tests:

1. **Model**: A_d = `expm`(A·0.1); two zero eigenvalues; a trimmed climb needs no elevator.
2. **Riccati**: `dare` against `scipy.linalg.solve_discrete_are` on random (2,1), (4,1), (4,2) systems, residual and closed-loop stability; the recursion from P_N = 0 converges to the DARE in 400 steps; the end-of-horizon gains differ.
3. **Condensed MPC**: prediction equals simulation; Su is causal; the condensed cost equals the stage-wise sum for five random (x₀, U); H > 0.
4. **The independent checks**: K_mpc = K₀ (Riccati) for N = 1, 3, 10, 40; with the DARE terminal cost K_mpc = K_∞ for N = 1, 5, 25, and `lqr_gain` with SciPy's DARE gives the same K_∞.
5. **QP**: a hand case (solution (0, 1), multiplier 2 on the sum constraint); KKT residuals below 10⁻⁹ and no better feasible perturbation on five random 12-variable, 30-constraint problems; zero iterations when nothing is violated; infeasibility detected.
6. **Closed loop**: over 15 s at N = 20 the elevator, its rate and the climb rate never exceed their limits by more than 10⁻⁹ and each limit is reached; the altitude converges within 5 cm without overshooting 5 cm; unconstrained LQR exceeds the elevator and climb limits three times over and clipped LQR still doubles the climb limit; with loose limits MPC reproduces the LQR trajectory to 10⁻⁹; N = 3 raises `ValueError`.

### Compare

There is no QP solver among the permitted oracles, so the QP is certified by its own KKT conditions, which for a convex problem is a proof of optimality, and by the random-perturbation test. Production MPC uses specialised solvers (active-set with warm starts, interior-point, or first-order methods such as ADMM) and often the uncondensed (sparse) formulation, which keeps the states as variables and scales linearly in N instead of cubically; the condensed form here is the clearest for teaching and fast enough for N ≤ 50. `scipy.linalg.solve_discrete_are` uses a QZ decomposition of the symplectic pencil; our fixed-point iteration is slower (it converges at the rate of the closed-loop spectral radius squared) but agrees to 10⁻⁸.

Running `main()` prints (abridged):

```
infinite-horizon LQR K = [0.0673, 2.1851, 0.9248, 0.3195]; closed-loop |z| = [0.746 0.746 0.891 0.891]
N =  5: |K_mpc(P=Q) - K_0 Riccati| = 1.4e-16; |K_mpc(P=DARE) - K_inf| = 4.4e-12; K_mpc(P=Q) = [0.0026, 0.5221, 0.1985, 0.1848]
N = 20: |K_mpc(P=Q) - K_0 Riccati| = 1.3e-15; |K_mpc(P=DARE) - K_inf| = 8.1e-12; K_mpc(P=Q) = [0.0558, 1.9418, 0.8355, 0.3049]
20 m climb, unconstrained LQR: peak elevator 77.1 deg, peak rate 771 deg/s, peak climb 9.84 m/s, settling (0.4 m) 4.2 s
   LQR clipped to the actuator: peak climb 8.21 m/s (limit 2.5), overshoot 0.02 m, settling 5.2 s
   MPC N = 30 (P = DARE): peak elevator 10.00 deg, peak rate 30.00 deg/s, peak climb 2.5000 m/s, overshoot 0.006 m, settling 9.4 s, QP iterations max 53
first QP: 180 constraints, 27 active, 41 iterations; KKT residuals stationarity 1.0e-11, primal 4.3e-15, dual 0.0e+00, complementarity 3.0e-12
```

and writes `mpc_altitude.png` (altitude, climb rate and elevator for unconstrained LQR, clipped LQR and MPC, with the limits drawn).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **Three controllers.** Run the 20 m climb with unconstrained LQR, LQR clipped to the elevator and elevator-rate limits, and MPC (N = 30, DARE terminal cost). Record peak elevator, peak elevator rate, peak climb rate, overshoot and the time to stay within 0.4 m. Predict: which controllers violate which limits, and the shortest possible climb time 20 m / 2.5 m/s = 8 s plus the time to start and stop the climb.
2. **Horizon and terminal cost.** For N = 3, 5, 10, 20 (and 30, 50 with the DARE cost), with P = Q and with P = DARE, record feasibility, overshoot, settling time, the error at 25 s and the largest number of QP iterations. Also record the spectral radius of the *unconstrained* closed loop A − BK_mpc for P = Q.
3. **Optimality check.** For the first QP record the number of constraints, active constraints, iterations and the four KKT residuals.

Raw data format: CSV with headers `controller, peak_elev_deg, peak_rate_deg_s, peak_climb_m_s, overshoot_m, settle_s`; `N, terminal, feasible, infeasible_step, overshoot_m, settle_s, error_25s_m, max_qp_iter, rho_unconstrained`; `n_constraints, n_active, iterations, stationarity, primal, dual, complementarity`.

What you should find. Unconstrained LQR is fast (4.2 s) and impossible: 77.1° of elevator, 771°/s of elevator rate, 9.84 m/s of climb. Clipping the elevator and its rate keeps the actuator legal, but the climb still reaches 8.21 m/s: clipping cannot see a state limit. MPC holds all three: peak elevator 10.00°, peak rate 30.00°/s, peak climb 2.5000 m/s (above the limit by at most 4 × 10⁻¹⁰ m/s, the solver tolerance), overshoot 6 mm, settled in 9.4 s against the 8 s bound for a climb at the limit. The horizon table:

| N | P = Q | P = DARE |
|---|---|---|
| 3 | feasible; never settles; error at 25 s −8.21 m | infeasible at step 8 (0.8 s) |
| 5 | feasible; never settles; error at 25 s −1.49 m | infeasible at step 7 (0.7 s) |
| 10 | settles in 12.4 s, no overshoot | settles in 9.6 s, overshoot 6 mm, ≤ 23 iterations |
| 20 | settles in 9.6 s, overshoot 8 mm | settles in 9.4 s, ≤ 33 iterations |
| 30 / 50 | — | 9.4 s, ≤ 53 / ≤ 71 iterations |

The two failures at short horizons are different. With P = Q a short horizon makes a *timid* controller: the unconstrained gain is weak (spectral radius 0.9957 at N = 3, 0.9886 at N = 5, against 0.9504 at N = 10 and 0.9006 at N = 20), so the aircraft creeps toward the target and is still 8 m short after 25 s. With P = DARE the terminal cost promises the aggressive LQR behaviour after the horizon, the controller starts the climb hard, and three to five samples are too few to see that the climb must be stopped within the elevator-rate limit: within 0.8 s no input sequence satisfies the climb-rate limit over the horizon and the QP reports infeasibility. From N = 6 on (Problem 6) the DARE controller is feasible throughout. The QP iterations grow roughly in proportion to N.

## 5. Validate and improve

| quantity | theory (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| K_mpc vs Riccati K₀, P = Q | equal | N = 5: 1.4 × 10⁻¹⁶, N = 20: 1.3 × 10⁻¹⁵ | (5) | — | rounding |
| K_mpc vs K_∞, P = DARE | equal for every N | N = 5: 4.4 × 10⁻¹², N = 20: 8.1 × 10⁻¹² | `solve_discrete_are` | 10⁻⁸ rel | DARE tolerance |
| first QP | KKT (6) | 41 iterations, 27 active | residuals ≤ 1.0 × 10⁻¹¹ | — | rounding |
| climb rate | ≤ 2.5 m/s | max 2.5000000004 m/s | — | 4 × 10⁻¹⁰ | QP tolerance |
| elevator, rate | ≤ 10°, ≤ 30°/s | 10.00°, 30.00°/s | — | rounding | — |
| climb time | ≥ 20/2.5 = 8 s (+ transients) | 9.4 s | — | 1.4 s | starting and stopping the climb |
| loose limits | MPC = LQR | trajectories equal | `simulate_lqr` | < 10⁻⁹ | — |

The equivalences hold to rounding, the solver is verified by its KKT conditions, and the constraint record is clean to the solver tolerance. What remains is the model question: the plant here is the model. The improvements that matter in flight are a disturbance model (an integrator on the climb-rate error, or an observer of a constant updraft), softened state constraints (Extension 1), warm starts to bound the iterations (Extension 2), and a terminal set for guaranteed feasibility (Extension 3). What not to improve: the QP accuracy. It is already at 10⁻¹¹.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (4), the equivalence K_mpc = K₀ and its DARE corollary; `mpc_altitude.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) choosing the horizon and terminal cost for the autopilot, with the constraint record, the climb time, the worst-case iterations and the first thing you would add before flight.

### Adversarial self-checks

- **What if the altitude change is larger?** The climb time grows by Δh/2.5 m/s and the structure of the solution does not change: for 35 m the MPC settles in 15.4 s against a 14 s bound (Problem 4), again about 1.4 s for starting and stopping.
- **What if the horizon goes to 1?** With P = DARE and no active constraint MPC is still exactly LQR (test), but with the limits it is infeasible within the first second; with P = Q it barely moves the aircraft.
- **Which assumption fails first?** Model equals plant. An updraft or a model error moves the climb rate the prediction did not expect, and hard state constraints turn that into infeasibility.
- **What measurement could falsify the design?** A flight log in which the climb rate exceeds 2.5 m/s while the elevator is not at a limit would show that the model's climb-rate prediction is wrong (L_α or V off), not that the optimiser failed.
- **Does the result survive a change of solver?** The KKT certificate does not depend on the solver; any correct QP solver must return the same U (H > 0 makes the solution unique).
- **Which decision changes if the climb-rate limit is tightened to 1.5 m/s?** The climb takes longer (about 20/1.5 = 13 s plus transients), and the horizon must be long enough to stop a slower climb, which is easier; the minimum feasible N should not increase. Check it.

## Check yourself

1. Write X = Sx x₀ + Su U and say why Su is block lower-triangular.
2. Write H and F of the condensed cost, and say why H is positive definite.
3. What is the unconstrained MPC law, and which LQR gain does it equal?
4. Why does a terminal cost from the DARE make the MPC gain independent of the horizon when no constraint is active?
5. How are an input-rate limit and a state limit written as GU ≤ w? Which of the two depends on the current state?
6. State the KKT conditions of the QP and say why they certify the global optimum here.
7. Why can clipping an LQR command never enforce a climb-rate limit?
8. What is recursive feasibility, and what happened to it at N = 3 in this lesson?

## References

The sources, with the exact lectures and sections to read, are in the **References** row at the top of this lesson: 16.323 Lectures 16 and 3 for the condensed QP and the Riccati recursion; Rawlings, Mayne and Diehl Chapters 1, 2 and 8 as graduate-level further reading (link only); 18.335J Lecture 29 for the KKT conditions; SciPy for the DARE oracle. The dual active-set QP solver and the altitude-hold study are original instruction (Partial GAP in `references/phase-09.md`).
