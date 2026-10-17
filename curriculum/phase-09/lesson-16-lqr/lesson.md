# P09.16  LQR

> Pole placement asks you where the poles should go. LQR asks what you are willing to pay, in error and in effort, and then puts the poles where that price is lowest, with margins you get for free.

| | |
|---|---|
| **Learning outcomes** | 1. Pose the infinite-horizon linear-quadratic regulator, derive the continuous algebraic Riccati equation (CARE) from the costate equations, and derive the optimal gain K = R⁻¹BᵀP and the optimal cost x₀ᵀPx₀. 2. Solve the CARE from scratch two independent ways, from the stable invariant subspace of the Hamiltonian matrix and by Newton–Kleinman iteration on Lyapunov equations, and verify both against `scipy.linalg.solve_continuous_are` and the ARE residual. 3. Choose Q and R by Bryson's rule, predict peak torque and settling time, and check the design by simulation with and without the actuator limit. 4. Derive the Kalman return-difference inequality and the guaranteed margins of single-input LQR (gain margin [½, ∞), phase margin ≥ 60°), and decide which real-world effects those guarantees do not cover. |
| **Prerequisites** | `P09.15` (state feedback, pole placement, observers), `P09.14` (controllability; Gramians by Lyapunov solves), `P09.13` (state-space models), `P09.10` (Nyquist, gain and phase margins), `P02.05` (matrix exponential), `P01.16` (eigenvalues and eigenvectors), `P03.14` (unconstrained optimization: Newton's method). Skill check: you can solve a 2 × 2 Lyapunov equation AᵀX + XA + Q = 0 by hand and read a phase margin off a Nyquist plot. |
| **Engineering question** | A 10 kg·m² small satellite must remove a 5° single-axis attitude error with a reaction wheel limited to 0.02 N·m of torque and 0.2 N·m·s of momentum. Which state-feedback gains minimise a sensible quadratic cost, how fast does the error settle, does the wheel stay inside its limits, what happens to a 20° error, and what stability margins does the design guarantee? |
| **Assumptions and limits** | Rigid body about one principal axis, decoupled from the other two (small angles, low rates, no gyroscopic coupling); the wheel is an ideal torque source with a hard torque limit and no friction or motor dynamics; the full state (angle and rate) is measured exactly and instantly; no disturbance torques (solar pressure, gravity gradient, magnetic) and no flexible appendages. LQR theory needs (A, B) stabilisable and (A, Q^{1/2}) detectable; the margins hold at the plant input, for one input, with state feedback, and are lost when an observer, a delay or unmodelled dynamics sit inside the loop. |
| **Analytical work** | The quadratic cost and its Bellman (dynamic-programming) argument; the costate equations and the Hamiltonian matrix; the invariant-subspace solution P = X₂X₁⁻¹; the eigenvalue symmetry of H; Newton–Kleinman as Newton's method on the Riccati map, with its monotone, quadratic convergence; the closed-form LQR for the double integrator, k₁ = √(q₁/r), k₂ = √(q₂/r + 2Jk₁); the cost of an arbitrary gain by a Lyapunov equation; Bryson's rule; the return-difference identity and the margins it implies. |
| **Computational deliverable** | `starter/main.py`: `lyapunov_solve`, `hamiltonian_matrix`, `care_hamiltonian`, `care_newton_kleinman`, `are_residual`, `lqr`, `bryson_weights`, `cost_of_gain`, `loop_transfer` (`satellite_model`, `wheel_limits`, `expm`, `bass_stabilizing_gain`, `phase_margin`, `gain_margin_interval`, `simulate`, `quadratic_cost`, `weight_sweep` and `main` are given). Tests prove: the Lyapunov solver matches `scipy.linalg.solve_continuous_lyapunov` with residual below 10⁻¹⁰; the Hamiltonian CARE solution matches `solve_continuous_are` to 10⁻⁸ on random 2- to 6-state systems with ARE residual below 10⁻¹⁰ and a stable closed loop; Newton–Kleinman agrees with both, decreases monotonically and converges quadratically; H has a spectrum symmetric about the imaginary axis whose stable half are the closed-loop poles; an unstabilisable system is refused; the double-integrator gains equal the closed form; the optimal cost equals x₀ᵀPx₀, equals the simulated integral, and is beaten by no perturbed gain; the unsaturated simulation equals the matrix exponential; the saturated loop respects the torque limit and converges; |1 + L(jω)| ≥ 1, PM ≥ 60° and stability for gains from 0.51 to 1000 on random single-input systems. |
| **Practical activity** | **computational experiment**: a sweep of the control weight R = ρR_Bryson from ρ = 0.1 to 10 recording peak torque, settling time, crossover and phase margin; the 20° slew with and without the torque limit; a method comparison (Hamiltonian, Newton–Kleinman, SciPy) on the satellite and on random systems. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `lqr_satellite.png`, the hand derivations of the CARE from the costate equations and of the double-integrator gains, the weight-sweep table, and a half-page design note stating the chosen weights, gains, predicted settling time, peak torque and momentum, and the margins with what they do not cover. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the CARE and the double-integrator closed form are derived; all public tests pass; the weight sweep is explained in terms of the error-effort trade; the note states which assumption behind the guaranteed margins fails first on a real spacecraft. |
| **References** | `@mit-ocw-16-30-feedback-control-systems` Topic 18 "Deterministic linear quadratic regulator (LQR)" (https://ocw.mit.edu/courses/16-30-feedback-control-systems-fall-2010/resources/mit16_30f10_lec18/): Linear Quadratic Regulator (LQR); Constrained Optimization; LQR Optimization; LQR Stability Margins (required-free): the costate equations, the Hamiltonian matrix and the algebraic Riccati equation, K = R⁻¹BᵀP, and the guaranteed gain margin (½, ∞) and phase margin of at least 60°. Topic 12 "Pole placement approach" (.../mit16_30f10_lec12/): Linear Quadratic Regulator; LQR Weight Matrix Selection (Bryson's Rules); the symmetric root locus; Example: B747 (required-free): choosing Q and R and what the weights do to the poles. `@mit-ocw-16-323-principles-of-optimal-control` Lecture 4 "HJB equation: ... continuous LQR" (PDF heading "HJB Equation") (https://ocw.mit.edu/courses/16-323-principles-of-optimal-control-spring-2008/resources/lec4/) (supplementary-free; graduate): the Bellman/HJB route to the differential Riccati equation and its steady state (§2.2); Lecture 3 "Dynamic programming: ... discrete LQR" (supplementary-free): the discrete counterpart used in `P09.20`. `@mit-ocw-18-06-linear-algebra` Lecture 25 "Symmetric matrices and positive definiteness" (supplementary-free): Q ⪰ 0, R ≻ 0, P ≻ 0. `@scipy-userguide` API reference "Linear algebra (scipy.linalg)": Matrix Equation Solvers (`solve_continuous_are`, `solve_continuous_lyapunov`) (supplementary-free): the test oracles only. Partial GAP (original instruction): the loaded notes form the Hamiltonian matrix but do not derive P = X₂X₁⁻¹ from its stable invariant subspace step by step, and no verified free source covers Newton–Kleinman iteration (with Bass's stabilising starting gain); both solvers (§2.3) are written here. |
| **Extensions** | ★ Add the wheel speed as a third state with a momentum-dumping torque from magnetorquers and design a three-state LQR that keeps the wheel away from saturation. ★ Replace the measured rate by a Luenberger observer (`P09.15`) and measure how much of the 77° phase margin survives (Doyle's 1978 counterexample shows it can be almost all lost). ★ Implement the matrix-sign-function or Schur method for the CARE and compare its accuracy with the eigenvector method on a badly scaled 6-state system. |
| **Time** | unmeasured |

## 1. The engineering problem

A small satellite points its camera, antenna or solar panels by turning its whole body. The usual actuator is a **reaction wheel**: a flywheel driven by a motor. When the motor applies torque τ to the wheel, the wheel applies −τ to the body, so the body turns one way and the wheel spins up the other way; the total angular momentum stays constant. Wheels are limited in two ways. The motor has a **torque limit** (here 0.02 N·m, typical of a wheel for a 50–100 kg spacecraft), and the wheel has a **momentum capacity** (here 0.2 N·m·s): spin it faster and it saturates, after which it can push the body no further in that direction.

`P09.15` placed the closed-loop poles of an attitude loop wherever we chose. That leaves the real question open: *where should they go?* Fast poles remove an error quickly but demand large torque, which a 0.02 N·m wheel cannot supply; slow poles are gentle but take minutes. The **linear-quadratic regulator** (LQR) turns that trade into a cost we can write down,

J = ∫₀^∞ (xᵀQx + uᵀRu) dt,

with Q weighting state errors and R weighting control effort, and finds the one linear feedback u = −Kx that minimises J. The engineer's choices move from pole locations, which have no units of "cost", to weights, which can be read off physical limits (Bryson's rule, 2.5).

Success criteria, written first. For a 5° error the design must keep the wheel torque at or below 0.02 N·m and the momentum change below 0.2 N·m·s, and settle to 2 % in under a minute. The Riccati solution must agree three ways (Hamiltonian eigenvectors, Newton–Kleinman, `scipy.linalg.solve_continuous_are`) to 10⁻⁸, with an ARE residual below 10⁻¹⁰. The cost x₀ᵀPx₀ must equal the simulated integral of the cost to 10⁻⁴ relative. And the loop must show the guaranteed margins: |1 + L(jω)| ≥ 1 at every frequency.

## 2. Understand and derive

### 2.1 The plant

About one principal axis, with small angles and rates so that coupling to the other axes is negligible, Newton's second law for rotation (`P04`-level mechanics) gives J θ̈ = u, where J is the moment of inertia [kg·m²] and u the wheel torque on the body [N·m]. With x = (θ, ω):

ẋ = A x + B u,  A = [[0, 1], [0, 0]],  B = [[0], [1/J]].   (1)

This is a **double integrator**: both eigenvalues at 0, controllable (rank [B, AB] = 2, `P09.14`), unstable without feedback (a constant disturbance torque makes θ grow as t²). The wheel's momentum changes by h_w(t) = −J(ω(t) − ω(0)) because the total is conserved, so the peak body rate sets the peak wheel momentum.

### 2.2 The cost and the Riccati equation

Take Q ⪰ 0 (symmetric, positive semidefinite) and R ≻ 0 (positive definite). The **value function** V(x₀) = min J starting from x₀ must be quadratic for a linear plant and a quadratic cost; write V(x) = xᵀPx with P symmetric. Over a short time dt the cost incurred is (xᵀQx + uᵀRu)dt and the state moves to x + (Ax + Bu)dt, so the principle of optimality (Bellman; 16.323 Lecture 4 carries this out as the HJB equation) gives

0 = min_u [ xᵀQx + uᵀRu + 2xᵀP(Ax + Bu) ].   (2)

The bracket is quadratic in u; setting its gradient 2Ru + 2BᵀPx to zero gives the optimal control

**u = −Kx,  K = R⁻¹BᵀP.**   (3)

Substituting back, (2) must hold for every x, so the symmetric matrix in the middle must vanish:

**AᵀP + PA − PBR⁻¹BᵀP + Q = 0.**   (4)

This is the **continuous algebraic Riccati equation** (CARE; 16.30 Topic 18, "LQR Optimization"). It is quadratic in P and has many solutions; the one we want is the **stabilising** solution, the unique symmetric P for which A − BK is Hurwitz (all eigenvalues in the left half plane). It exists when (A, B) is stabilisable and (A, Q^{1/2}) is detectable (every unstable mode is either controllable, and penalised by Q so the regulator cares about it). Then P ⪰ 0, and the **optimal cost is J* = x₀ᵀPx₀**.

Units check for the satellite: x is (rad, rad/s), u in N·m. Q₁₁ has units 1/rad², R has 1/(N·m)², so J is dimensionless times seconds, and P₁₁ has units s/rad². K = R⁻¹BᵀP then has (N·m)²·(1/(kg·m²))·s/rad² = N·m/rad, as a gain from angle to torque must. Doing this check once catches most errors in building Q and R.

### 2.3 Two ways to solve the CARE

**The Hamiltonian matrix.** Pontryagin's principle (or the same calculation with a costate λ = 2Px, Problem 2) gives the pair of linear equations ẋ = Ax − BR⁻¹Bᵀ(λ/2), d(λ/2)/dt = −Qx − Aᵀ(λ/2), i.e. d/dt(x, λ/2) = H(x, λ/2) with

**H = [[A, −BR⁻¹Bᵀ], [−Q, −Aᵀ]].**   (5)

(16.30 Topic 18 forms this Hamiltonian matrix from the costate equations.)

H is a **Hamiltonian matrix**: with J₂ = [[0, I], [−I, 0]], J₂HJ₂⁻¹ = −Hᵀ, so H and −Hᵀ have the same eigenvalues, and the spectrum of H is symmetric about the imaginary axis (λ and −λ come in pairs). If no eigenvalue lies on the axis there are exactly n stable ones. Stack their eigenvectors as the columns of [X₁; X₂]. Along the optimal trajectory, (x, λ/2) stays in that stable invariant subspace, and λ/2 = Px means X₂ = PX₁:

**P = X₂ X₁⁻¹.**   (6)

The stable eigenvalues of H are precisely the closed-loop poles eig(A − BK), which the tests check. For the satellite, H has eigenvalues ±0.1159, ±0.1977 rad/s, and the closed-loop poles are −0.1159 and −0.1977 rad/s. This is the method behind `care_hamiltonian`: one 2n × 2n eigenproblem, then one n × n inverse. (Production codes use an ordered Schur form instead of eigenvectors because eigenvectors can be ill-conditioned when eigenvalues nearly coincide; Extension 3.)

**Newton–Kleinman.** Write the CARE as F(P) = 0 and apply Newton's method (`P03.14`). Linearising F about the current iterate turns each Newton step into a **Lyapunov equation**:

(A − BK_k)ᵀP_k + P_k(A − BK_k) + Q + K_kᵀRK_k = 0,  K_{k+1} = R⁻¹BᵀP_k.   (7)

Read (7) as "P_k is the cost matrix of the gain K_k" (2.4): each step evaluates the current gain exactly and then improves it, which is policy iteration. Kleinman's theorem (1968): if K₀ is stabilising, every K_k is stabilising, P₁ ⪰ P₂ ⪰ … ⪰ P, and the convergence is quadratic near the solution. The catch is the start: the double integrator is not stable, so K₀ = 0 is not allowed. `bass_stabilizing_gain` supplies one by Bass's method: solve (A + βI)W + W(A + βI)ᵀ = 2BBᵀ and take K₀ = BᵀW⁻¹, which moves every pole to the left of −β. For the satellite with β = 2 it gives K₀ = [80, 40] and poles at −2 ± 2j: absurdly aggressive (80 N·m per radian) but stabilising, which is all Newton needs. From there the iteration takes 14 steps; the relative error falls 32, 13.6, 6.6, 3.5, 1.9, 1.0, 0.48, 0.17, 0.036, 1.9 × 10⁻³, 5.6 × 10⁻⁶, 4.6 × 10⁻¹¹, 1.7 × 10⁻¹⁵: slow while far away, then each error roughly 1.5 times the square of the last.

Each Lyapunov equation AᵀX + XA + Q = 0 is linear in X. `lyapunov_solve` writes it as (I ⊗ Aᵀ + Aᵀ ⊗ I) vec X = −vec Q, an n² × n² linear system (⊗ is the Kronecker product), which is fine for the small systems here (n ≤ 10) and is how the method is taught; libraries use the Bartels–Stewart algorithm, O(n³) instead of O(n⁶).

### 2.4 The cost of any gain

For any stabilising K, not just the optimal one, u = −Kx gives ẋ = (A − BK)x and J = x₀ᵀX x₀ with

(A − BK)ᵀX + X(A − BK) + Q + KᵀRK = 0.   (8)

Proof: d/dt(xᵀXx) = xᵀ[(A − BK)ᵀX + X(A − BK)]x = −xᵀ(Q + KᵀRK)x; integrate from 0 to ∞ with x(∞) = 0. `cost_of_gain` implements (8); the tests use it to show that no random perturbation of the LQR gain lowers the cost. For the 5° error, J(K) = 13.6845, the trapezoidal integral of the simulated trajectory is 13.6847 (a difference of 10⁻⁵ relative, from the 0.05 s quadrature), and J(0.5K) = J(2K) = 14.4817: halving or doubling the gain costs 6 % more. That J(kK) = J(K/k) is an accident of this plant and initial state (Problem 6).

### 2.5 The satellite by hand, and Bryson's rule

For the double integrator with Q = diag(q₁, q₂) and R = r, write P = [[p₁, p₂], [p₂, p₃]] and b = 1/J. The (1,1), (1,2) and (2,2) entries of (4) give q₁ − b²p₂²/r = 0, p₁ − b²p₂p₃/r = 0 and 2p₂ − b²p₃²/r + q₂ = 0, so

**k₁ = √(q₁/r),  k₂ = √(q₂/r + 2J k₁).**   (9)

The angle gain depends only on the ratio of the angle weight to the torque weight; the rate gain grows with inertia, as it should (a heavier body needs more damping torque per unit rate).

**Bryson's rule** (16.30 Topic 12, "LQR Weight Matrix Selection (Bryson's Rules)") sets each weight to one over the square of the largest acceptable value: Q = diag(1/θ_max², 1/ω_max²), R = 1/u_max². It makes each term of the cost dimensionless and of order one when its variable reaches its limit, which is a sensible first guess and nothing more. With θ_max = 5°, ω_max = 0.5°/s and u_max = 0.02 N·m: Q = diag(131.31, 13131.2), R = 2500, and (9) gives **k₁ = 0.2292 N·m/rad, k₂ = 3.1363 N·m·s/rad**, the same to every printed digit as both numerical CARE solutions. The closed-loop poles −0.1159 and −0.1977 rad/s are real (ζ = 1.000): the rule happens to give a critically damped loop here. For a 5° error the initial torque is k₁ × 5° = 0.0200 N·m: exactly the wheel limit, because Bryson's rule with θ_max equal to the initial error makes √(q₁/r) θ₀ = u_max. The loop settles to 2 % in 41.2 s, and the body rate peaks at 0.272°/s, so the wheel's momentum changes by 0.0475 N·m·s, a quarter of its capacity.

### 2.6 Guaranteed margins

Break the loop at the plant input. The loop transfer function is L(s) = K(sI − A)⁻¹B. Multiply (4) by −1, add and subtract sP, and pre- and post-multiply by Bᵀ(−sI − Aᵀ)⁻¹ and (sI − A)⁻¹B. After some algebra (Problem 2d) this gives the **return-difference identity**

(1 + L(−jω))ᵀ R (1 + L(jω)) = R + G(−jω)ᵀ Q G(jω),  G = (jωI − A)⁻¹B.   (10)

The right side is at least R because Q ⪰ 0, so for a single input **|1 + L(jω)| ≥ 1 at every frequency**: the Nyquist curve of L never enters the unit disk centred at −1. Two consequences follow from geometry. That disk covers the segment (−2, 0) of the real axis, so the curve can cross the negative real axis only at points −x with x ≥ 2. Scaling L by k moves such a crossing to −kx, and the encirclement count of −1 (`P09.10`) can change only when a crossing passes through −1, at k = 1/x ≤ ½. So every k > ½ leaves the closed loop stable: **gain margin [½, ∞)**. And at the gain crossover |L| = 1 the point L lies on the unit circle about 0 and outside the unit circle about −1; the two circles meet at angles ±120°, so **phase margin ≥ 60°** (the margins 16.30 Topic 18 lists under "LQR Stability Margins"). For the satellite: PM = 77.2° at ω_c = 0.3216 rad/s, min |1 + L| = 1.0000 (reached at high frequency, where L → 0), and a scan of k from 0.01 to 10⁴ found no instability (the double integrator with these gains is stable for every k > 0, which is more than the guarantee).

What the guarantee does not say: it holds for one input, with the full state fed back, at the plant input. A pure delay τ adds phase lag ωτ; the **delay margin** PM/ω_c = 1.347 rad / 0.3216 rad/s = 4.2 s is large for this slow loop. An observer between the sensors and the gain can remove the margins entirely (Doyle, 1978; Extension 2), and an unmodelled flexible mode near crossover can do the same (`P09.18`).

### 2.7 Where the model stops being valid

- **Torque saturation.** LQR is linear; u = −Kx does not know the wheel limit. For a 20° error the linear law asks 0.0800 N·m, four times the limit (2.5, part 4).
- **Momentum saturation.** Not modelled as a state here. For the 20° slew the wheel momentum change reaches 0.1763 N·m·s with the torque limit, 88 % of capacity; any stored momentum from earlier manoeuvres or disturbances could saturate it (Extension 1).
- **Coupling.** A rotation about one principal axis stays decoupled; simultaneous rates about two axes couple through the gyroscopic term ω × Jω, and a three-axis slew needs the full rigid-body model (Phase 06).
- **Sensors.** The rate comes from gyros and the angle from a star tracker, both noisy and delayed; the regulator assumes perfect knowledge.
- **Flexibility.** Solar arrays and antennas have modes at fractions of a hertz; a rigid-body design is safe only if its crossover is well below them.

## 3. Build the computational model

Open `starter/main.py`. The plant, the limits, a matrix exponential, Bass's stabilising gain, the margin and simulation helpers and `main` are given.

### Step 1: Lyapunov equations

```python
def lyapunov_solve(A, Q):
    n = A.shape[0]
    Mk = np.kron(np.eye(n), A.T) + np.kron(A.T, np.eye(n))     # acts on vec(X), column-major
    if np.linalg.cond(Mk) > 1e14:
        raise ValueError("eigenvalues of A sum to zero")
    X = np.linalg.solve(Mk, -Q.reshape(-1, order="F")).reshape((n, n), order="F")
    return 0.5 * (X + X.T)                                       # symmetric Q gives symmetric X
```

The `order="F"` matters: vec stacks columns, and NumPy's default reshape stacks rows.

### Step 2: the Hamiltonian solution

`hamiltonian_matrix` builds (5) with `np.linalg.solve(R, B.T)` rather than an explicit inverse. `care_hamiltonian` takes `np.linalg.eig(H)`, refuses if any eigenvalue has |Re λ| below 10⁻⁹ times the largest |λ| (no stabilising solution), selects the n columns with Re λ < 0, splits them into X₁ and X₂, and returns the real part of X₂X₁⁻¹, symmetrised. The eigenvectors are complex when the poles are, but X₂X₁⁻¹ is real up to rounding because complex eigenvectors come in conjugate pairs.

### Step 3: Newton–Kleinman

`care_newton_kleinman` starts from K₀ (zero if A is Hurwitz, Bass otherwise), checks that A − BK₀ is Hurwitz, and loops (7) until ‖P_k − P_{k−1}‖_F ≤ 10⁻¹² ‖P_k‖_F, returning P and the iteration count.

### Step 4: the regulator and its checks

`are_residual` returns ‖AᵀP + PA − PBR⁻¹BᵀP + Q‖_F divided by the larger of ‖Q‖_F and ‖AᵀP + PA‖_F, so that it measures cancellation relative to the size of the terms. `lqr` calls either solver and returns (K, P, eig(A − BK)). `bryson_weights` is two `np.diag` calls. `cost_of_gain` solves (8) and returns `inf` for a destabilising gain. `loop_transfer` evaluates K(jωI − A)⁻¹B with a linear solve at each frequency.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-16-lqr`. The tests:

1. **Lyapunov**: residual below 10⁻¹⁰‖Q‖ and agreement with `scipy.linalg.solve_continuous_lyapunov(A.T, -Q)` for n = 2, 3, 5 (SciPy solves AX + XAᴴ = Q, hence the transposes and the sign); X positive semidefinite; a singular case raises.
2. **CARE, the independent check**: Hamiltonian P against `solve_continuous_are` to 10⁻⁸ on random (2,1), (3,1), (4,2), (6,3) systems, residual below 10⁻¹⁰, P positive definite, A − BK Hurwitz; Newton–Kleinman against both on three more systems.
3. **Newton–Kleinman theory**: on the satellite, every P_k ⪰ P_{k+1} ⪰ P and, once the error is below 0.1, e_{k+1} ≤ 3e_k²; K = 0 is refused as a start.
4. **Hamiltonian structure**: the spectrum equals its own negative; its stable half are the closed-loop poles; an unstabilisable plant is refused.
5. **Satellite**: both solvers reproduce (9); Bryson's rule; J* = x₀ᵀPx₀ = the simulated integral to 10⁻⁴; 20 random perturbed gains all cost more; −K costs ∞.
6. **Simulation**: unsaturated RK4 equals e^{(A−BK)t}x₀ to 10⁻¹⁰; the saturated loop never exceeds 0.02 N·m and converges from 20°; a zero initial state produces no motion.
7. **Margins**: on four random single-input systems, |1 + L| ≥ 1 − 10⁻⁹ on 3000 frequencies, PM ≥ 60° and stability for every k on a grid from 0.51 to 1000; the low- and high-frequency asymptotes of L; the weight sweep ordering.

### Compare

`scipy.linalg.solve_continuous_are` uses a QZ (generalised Schur) decomposition of the Hamiltonian pencil, which avoids both the explicit R⁻¹ and the eigenvectors. On the satellite it agrees with `care_hamiltonian` to 1.9 × 10⁻¹⁵ relative and with Newton–Kleinman to 10⁻¹⁶: the satellite is well conditioned, and all three methods are at rounding level. The differences appear on hard problems: eigenvectors become inaccurate when two closed-loop poles nearly coincide (here they are 0.08 rad/s apart, ample), and Newton–Kleinman needs a stabilising start and n² × n² solves in our Kronecker form. Its practical virtue is refinement: one or two Newton steps from any other method's answer remove the residual, which is how some production codes polish a Schur solution.

Running `main()` prints the design and its checks:

```
Bryson: theta_max = 5 deg, omega_max = 0.5 deg/s, u_max = 0.02 N m -> Q = diag(131.31, 13131.2), R = 2500
P (Hamiltonian) = [[1796.94697, 5729.57795], [5729.57795, 78406.6026]]
Newton-Kleinman: 14 iterations from the Bass gain, max |P_n - P_h| / max|P_h| = 8.1e-17
ARE residuals: Hamiltonian 1.4e-15, Newton 4.2e-16
K = [0.2292 N m/rad, 3.1363 N m s/rad]; closed-loop poles [-0.1159+0.j -0.1977+0.j] rad/s
5 deg error: peak torque 0.0200 N m, cost x0'Px0 = 13.6845, simulated cost 13.6847
   2 % settling 41.2 s, peak body rate 0.272 deg/s, peak wheel momentum change 0.0475 N m s
loop at the plant input: phase margin 77.2 deg at 0.3216 rad/s; min |1 + L| = 1.0000; delay margin 4.2 s
```

and writes `lqr_satellite.png` (attitude and torque for 5° and 20°, with and without the limit; the Nyquist plot of 1 + L outside the unit circle).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running.

1. **The error–effort trade.** Scale the control weight, R = ρR_Bryson, for ρ = 0.1, 0.3, 1, 3, 10, and record the slowest closed-loop pole, the 2 % settling time of a 5° error, the peak torque, the crossover and the phase margin (`weight_sweep`). Predict from (9): k₁ ∝ ρ^{−1/2}, so the peak torque k₁θ₀ scales as ρ^{−1/2}; the crossover should fall as ρ rises.
2. **Saturation.** Command a 20° error with and without the 0.02 N·m limit (`simulate(..., u_max=0.02)`). Predict the time the wheel spends at the limit and whether the limited response is slower or overshoots.
3. **Method comparison.** For 20 random systems of each size n = 2, 4, 6, 8 (seeded), record the relative difference between `care_hamiltonian`, `care_newton_kleinman` and `solve_continuous_are`, the ARE residual of each, and the Newton iteration count. Predict which method degrades first as n grows.
4. **Inertia.** Repeat the design for J = 5, 10, 20, 40 kg·m² with the same weights. Predict from (9) which gain changes and how the settling time moves.

Raw data format: CSV tables with headers `rho, slowest_pole, settle_s, peak_torque_Nm, pm_deg, wc_rad_s`; `case, peak_torque_Nm, time_at_limit_s, settle_s, peak_rate_deg_s, peak_momentum_Nms`; `n, seed, diff_ham_scipy, diff_nk_scipy, res_ham, res_nk, nk_iters`; `J, k1, k2, poles, settle_s`. Each header line names the fixed parameters.

What you should find. The sweep gives peak torques of 0.0632, 0.0365, 0.0200, 0.0115 and 0.0063 N·m: the predicted ρ^{−1/2} law exactly, since the peak is k₁θ₀ at t = 0. The settling time hardly moves (40.2, 40.8, 41.2, 42.0, 46.6 s) while the crossover falls from 0.8234 to 0.1487 rad/s and the phase margin from 83.9° to 70.9°. That is a lesson in itself: below ρ = 1 the extra torque buys a faster *initial* response but not faster settling, because the slowest pole stays near −0.10 rad/s. Formula (9) explains it: as r → 0 the characteristic polynomial s² + (k₂/J)s + k₁/J has a fast root near −k₂/J and a slow root near −k₁/k₂ → −√(q₁/q₂) = −ω_max/θ_max = −0.1 rad/s (at ρ = 10⁻⁴ the slow pole is −0.100001 rad/s). The cost asks the rate to stay below ω_max, so the error cannot close faster than ω_max/θ_max however much torque is offered. The weight to change for faster settling is ω_max, not u_max. For the 20° error the limited loop sits at 0.02 N·m for 7.6 s and settles in 44.6 s instead of 41.2 s, without overshoot (the critically damped design leaves room), and the wheel momentum change rises to 0.1763 N·m·s, 88 % of capacity: the torque limit is harmless, the momentum limit is the real constraint for large slews. The method comparison (seed 16, 20 systems per size, one input, R = 1) is the surprise. The Hamiltonian solution stays within 5 × 10⁻¹¹ of SciPy for every size, with residuals below 2 × 10⁻¹⁰. Newton–Kleinman matches it for n = 2, 4 and 6 (median 8, 15 and 22 iterations; one n = 4 case needs 58 and lands 2.6 × 10⁻⁹ from SciPy). At n = 8 it fails in 19 of 20 cases: Bass's gains reach 10⁷–10⁸, the start is not stabilising once rounded (3 cases), or an early iterate loses stability (16 cases). Before the solver checked every iterate, those runs converged to an *indefinite* P with a residual of 10⁻¹⁵ that differed from SciPy's by about 100 %: a non-stabilising root of the same equation. **A small residual does not certify the stabilising solution; check eig(A − BK) as well.** A gentler start fixes it: Bass with β = 0.5 + max Re λ(A) instead of 1 + ‖A‖₁ reaches SciPy's answer to 10⁻⁸ in all 20 n = 8 cases. Newton is an excellent refiner and a fragile starter. Finally, inertia: with the same weights k₁ = 0.2292 N·m/rad for every J (as (9) says), while k₂ is 2.7467, 3.1363, 3.7973 and 4.8567 N·m·s/rad for J = 5, 10, 20 and 40 kg·m²; the poles go from real (−0.1026, −0.4467 rad/s at 5 kg·m²) to complex (−0.0607 ± 0.0452j rad/s at 40 kg·m²) and the 2 % settling time from 40.7 s to 49.9 s.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| k₁ (N·m/rad) | √(q₁/r) = 0.2292 | 0.2292 (both solvers) | SciPy P: same | < 10⁻¹⁴ rel | rounding |
| k₂ (N·m·s/rad) | √(q₂/r + 2Jk₁) = 3.1363 | 3.1363 | SciPy P: same | < 10⁻¹⁴ rel | rounding |
| P, Hamiltonian vs SciPy | — | — | 1.9 × 10⁻¹⁵ rel | — | QZ vs eigenvectors |
| P, Newton vs Hamiltonian | — | 8.1 × 10⁻¹⁷ rel | — | — | rounding |
| ARE residual | 0 | 1.4 × 10⁻¹⁵ / 4.2 × 10⁻¹⁶ | — | — | rounding |
| optimal cost, 5° | x₀ᵀPx₀ = 13.6845 | Lyapunov (8): 13.6845 | simulated integral 13.6847 | 1.5 × 10⁻⁵ rel | trapezoid, 0.05 s |
| peak torque, 5° | k₁θ₀ = 0.0200 N·m | 0.0200 N·m | — | 0 | u(0) is the peak |
| min |1 + L| | ≥ 1 | 1.0000 | — | — | Kalman |
| phase margin | ≥ 60° | 77.2° at 0.3216 rad/s | — | — | — |
| closed-loop poles | stable half of H | −0.1159, −0.1977 rad/s | eig(H): ±0.1159, ±0.1977 | rounding | — |

The numbers agree to rounding: the Riccati machinery is right. The validation question that remains is about the model. The rigid single-axis double integrator is the textbook idealisation of a slew; the first real effects to add are the wheel's momentum (a state that integrates −u), disturbance torques that make the wheel drift toward saturation over an orbit, and the observer that turns gyro and star-tracker data into the state. Each makes LQR more useful (more states to trade against each other through Q) but erodes the guaranteed margins, which is why a flight design is checked by the robustness analysis of `P09.18` rather than by the theorem.

What to improve: add the wheel momentum state (Extension 1); weight θ and ω by the pointing requirement rather than by Bryson's rule; check a disturbance torque of 10⁻⁵ N·m. What not to improve: the CARE solver. All three agree to 10⁻¹⁵ here.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (3)–(4) by the Bellman argument and of (6) from the Hamiltonian; the double-integrator closed form (9); `lqr_satellite.png`; the four experiment tables; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving your chosen weights, the gains, predicted settling time, peak torque and momentum for 5° and 20°, the guaranteed and measured margins, and the assumption you would test first.

### Adversarial self-checks

- **What if the inertia doubles?** k₁ is unchanged (it depends on q₁/r only), k₂ grows, the poles slow down: for J = 20 kg·m² they become −0.0949 ± 0.0495j rad/s (Problem 4) and the 5° response settles in 42.6 s. The peak torque is still k₁θ₀ = 0.02 N·m.
- **What if the initial error is 20°?** The linear law asks 0.08 N·m; saturation lengthens settling only by 3.4 s, but the momentum change approaches capacity. A designer who checks only the torque misses the real limit.
- **Which assumption fails first?** Perfect state knowledge. The gyro rate is noisy and the star tracker is delayed; with an observer the guaranteed margins are gone, and the 4.2 s delay margin is the more honest robustness figure for latency.
- **What measurement could falsify the model?** A slew on the spacecraft (or on an air-bearing table) whose rate history departs from e^{(A−BK)t}x₀ by more than the gyro noise: a steady drift means a disturbance torque, an oscillation near the crossover means a flexible mode.
- **Does the result survive a change of solver?** Hamiltonian, Newton–Kleinman and QZ agree to 10⁻¹⁵; the closed form (9) agrees to every digit.
- **Which decision changes if the torque limit were half as large?** At 0.01 N·m, Bryson's rule makes R four times larger (ρ = 4), k₁ halves, the peak torque for 5° is 0.0100 N·m, and the settling time rises only from 41.2 s to 42.6 s with a phase margin of 73.1°: the weaker wheel would be acceptable for this manoeuvre.

## Check yourself

1. Write the LQR cost, the CARE and the optimal gain, and give the units of K for the satellite.
2. Why must Q be positive semidefinite and R positive definite, and what goes wrong if R = 0?
3. What is the Hamiltonian matrix of the LQR problem, why are its eigenvalues symmetric about the imaginary axis, and which half are the closed-loop poles?
4. Why does Newton–Kleinman need a stabilising initial gain, and what does each iteration solve?
5. State Bryson's rule. For θ_max = 2°, ω_max = 0.2°/s and u_max = 0.01 N·m, write Q and R.
6. For the double integrator Jθ̈ = u, give the LQR gains in closed form. Which gain depends on J?
7. State the guaranteed gain and phase margins of single-input LQR and the inequality they come from.
8. Name two conditions under which the guaranteed margins no longer hold.

## References

The sources, with the exact topics and sections to read, are in the **References** row at the top of this lesson: 16.30 Topics 18 and 12 for the Riccati equation, Bryson's rule and the guaranteed margins; 16.323 Lecture 4 for the dynamic-programming (HJB) derivation; SciPy for the test oracles. The eigenvector derivation of the CARE solution and Newton–Kleinman iteration are original instruction (Partial GAP in `references/phase-09.md`).
