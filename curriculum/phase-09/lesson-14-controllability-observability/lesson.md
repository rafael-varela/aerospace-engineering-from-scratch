# P09.14  Controllability & Observability

> Before you design a controller, ask whether the actuators can move every state. Before you design an estimator, ask whether the sensors can see every state. Two rank tests answer yes or no; the Gramians say how hard.

| | |
|---|---|
| **Learning outcomes** | 1. Define controllability and observability, derive the Kalman rank tests from Cayley–Hamilton and the Popov–Belevitch–Hautus (PBH) eigenvalue test from left eigenvectors, and show that observability of (A, C) is controllability of (Aᵀ, Cᵀ). 2. Implement and test both rank tests, the PBH test (which names the lost mode), the unobservable subspace, a Lyapunov solver and the controllability and observability Gramians (infinite-horizon by Lyapunov, finite-horizon by Van Loan's block exponential), and verify that the rank tests and PBH agree and that each Gramian equals quadrature of its defining integral. 3. Apply them to a reaction-wheel satellite: decide which wheel failures leave the attitude controllable and which star-tracker and gyro combinations leave it observable, and quantify "how controllable" by the smallest Gramian eigenvalue. 4. Compute a minimum-energy rest-to-rest slew from the Gramian and check its torque against a wheel limit. |
| **Prerequisites** | `P09.13` (state-space models, similarity transforms), `P01.14` (rank, null space, SVD), `P01.16` (eigenvalues and eigenvectors), `P02.05` (matrix exponential, Van Loan-type block matrices), `P03.14` (least-norm solutions, used for minimum energy). Skill check: you can find the rank and null space of a 3 × 4 matrix with `np.linalg.svd`. |
| **Engineering question** | A small satellite holds attitude with reaction wheels and measures it with a star tracker and gyros. With three orthogonal wheels, does losing one wheel lose an axis? With a four-wheel pyramid, does it? If the gyros fail, or if the star tracker measures only two axes, can the attitude still be estimated? And how fast can a 10° slew be done with 10 mN m wheels? Answer: one orthogonal-wheel failure leaves rank 4 of 6 (the rigid-body mode at s = 0 about that axis is uncontrollable); the pyramid stays controllable after any single failure (smallest torque-authority eigenvalue 1/3 of a nominal wheel) but not after two; a full star tracker keeps everything observable even with a failed gyro, while a two-axis tracker leaves the boresight angle unobservable whatever gyros are added; the minimum-energy 10° slew about a 3 kg m² axis needs 6JΔθ/T² of peak torque, so at least 17.7 s. |
| **Assumptions and limits** | Linear, time-invariant models. The satellite is a rigid body linearised about an inertial attitude hold: three decoupled double integrators J_i θ̈_i = τ_i, with no gyroscopic coupling (small rates, zero stored momentum), no gravity-gradient or orbital terms, no wheel momentum limits, and ideal sensors. The inertia diag(4, 5, 3) kg m² and the 10 mN m wheel limit are representative small-satellite numbers, not a specific spacecraft. Rank tests are yes/no and need a numerical tolerance; they say nothing about actuator limits, which the Gramian begins to address. |
| **Analytical work** | The controllability and observability matrices from the series of e^{At}B and Cayley–Hamilton; PBH from left eigenvectors; duality; invariance under x = Tz; the Lyapunov equation for the infinite Gramian and the differential Lyapunov identity for the finite one; Van Loan's block exponential; the double-integrator Gramian in closed form; the minimum-energy input and its energy; the rest-to-rest slew u(t) = (6JΔθ/T²)(1 − 2t/T) and its 4/3 advantage over bang-bang. |
| **Computational deliverable** | `starter/main.py`: `ctrb`, `obsv`, `is_controllable`, `is_observable`, `pbh_uncontrollable_modes`, `pbh_unobservable_modes`, `unobservable_subspace`, `lyap`, `ctrb_gramian`, `obsv_gramian`, `ctrb_gramian_finite`, `min_energy_input`, `satellite_model` (`expm`, `numerical_rank`, `satellite_inertia`, `wheel_axes`, `short_period_model`, `main` are given). Tests prove: the Kalman structure and duality; **PBH and the rank test agree on six random systems** (three controllable, three built uncontrollable), and PBH names the right eigenvalue; ranks survive x = Tz; the wheel-failure and sensor cases give the verdicts above, with the boresight angle as the unobservable direction; `lyap` equals `scipy.linalg.solve_continuous_lyapunov` with zero residual; **the infinite Gramian and the finite Gramian equal Simpson quadrature of their integrals**, and the finite one satisfies AW + WAᵀ + BBᵀ = e^{AT}BBᵀe^{AᵀT}; the double-integrator closed form; a singular Gramian when a wheel axis is lost; the minimum-energy input reaches its target in `solve_ivp` with the predicted energy, matches the closed-form slew, beats bang-bang by 4/3 and scales as 1/T³. |
| **Practical activity** | **computational experiment**: a failure-case study (wheel and sensor configurations: rank, PBH, unobservable directions, smallest Gramian eigenvalue) and a slew-time study (peak torque and energy against slew time, against the wheel limit). |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `controllability.png`, the derivations of the rank tests, PBH and the minimum-energy input, the failure-case and slew tables, and a half-page note on which wheel and sensor configuration you would fly and why. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the rank tests and PBH are derived; all public tests pass; the note identifies the uncontrollable or unobservable direction in each failed case and uses a Gramian number, not just a rank, to compare the surviving configurations. |
| **References** | `@mit-ocw-18-06-linear-algebra` Lecture 10 "The four fundamental subspaces", Lecture 29 "Singular value decomposition", Lecture 21 "Eigenvalues and eigenvectors" (required-free: rank, null space and the SVD behind every test here). `@scipy-userguide` "Linear Algebra (scipy.linalg)" (supplementary-free: the oracle routines used in the tests). The Phase 09 reference map (`references/phase-09.md`) did not exist when this lesson was written, so no section of a control text is mapped yet: the PID, compensation, state-space and observer material here is original instruction (GAP), to be linked to exact sections when the map is published. |
| **Extensions** | ★ Add the stored wheel momentum h along one axis (gyroscopic coupling Jω̇ + ω × h = τ) and find whether a single remaining wheel can then control two axes. ★ Compute the Kalman decomposition (controllable/observable subspaces) of the failed configurations by SVD. ★ Balanced realisation: transform so that the controllability and observability Gramians are equal and diagonal, and truncate the flexible-array model of `P09.13` Problem 4. |
| **Time** | unmeasured |

## 1. The engineering problem

Spacecraft lose actuators and sensors. A reaction wheel's bearing wears out; a star tracker is blinded by the Sun or the Moon; a gyro drifts beyond use. Several missions have had to continue with fewer wheels or sensors than they launched with, and some could not continue at all. The first question after a failure is structural, not numerical: *can* the remaining actuators still move every attitude state, and *can* the remaining sensors still determine every state? If the answer is no, no controller or estimator can fix it, however clever. If the answer is yes, the second question is quantitative: how much torque, energy or time does it now cost?

The first question has an exact answer: a pair of rank tests on matrices built from A, B and C. Rank tests have a weakness. They say "rank 6" or "rank 4" but not *which* motion was lost, and a nearly-failed wheel still passes. The PBH test names the lost mode by its eigenvalue, and the null space of the observability matrix names the unseen state. The **Gramian**, an integral of how the input spreads into the state, measures how hard each direction is to reach. Its smallest eigenvalue drops toward zero as the system approaches uncontrollability, and its inverse gives the minimum energy needed for a manoeuvre.

Success criteria, written first: the rank test and PBH give the same verdict on every test system; the wheel and sensor verdicts match physical reasoning (an axis without a wheel cannot be turned; an angle without a sensor cannot be found); each Gramian equals quadrature of its integral to 10⁻⁵ or better; the minimum-energy slew reaches its target with the energy the Gramian predicts.

## 2. Understand and derive

### 2.1 Controllability and the Kalman rank test

(A, B) is **controllable** if for every initial state x₀, target x_f and time T > 0 there is an input that takes x₀ to x_f in time T. From x(T) = e^{AT}x₀ + ∫₀ᵀ e^{A(T−τ)}Bu(τ)dτ (`P02.05`), every reachable change lies in the span of the vectors e^{At}B for 0 ≤ t ≤ T. Expanding e^{At} = Σ Aᵏtᵏ/k!, and using Cayley–Hamilton (`P09.13`) to write every Aᵏ with k ≥ n as a combination of I, A, …, A^{n−1}, that span equals the column space of

**W_c = [B, AB, A²B, …, A^{n−1}B]**   (n × nm).   (1)

So (A, B) is controllable **if and only if rank W_c = n**. In different coordinates x = Tz, W_c becomes T⁻¹W_c: the rank does not change.

### 2.2 Observability and duality

(A, C) is **observable** if the output y(t) = Ce^{At}x₀ over any interval determines x₀ uniquely. A state x₀ is invisible exactly when Ce^{At}x₀ = 0 for all t, i.e. CAᵏx₀ = 0 for every k, i.e. x₀ is in the null space of

**W_o = [C; CA; CA²; …; CA^{n−1}]**,   (2)

so (A, C) is observable **if and only if rank W_o = n**. Since W_o = W_c(Aᵀ, Cᵀ)ᵀ, observability of (A, C) is controllability of the **dual** (Aᵀ, Cᵀ). One algorithm serves both (`P09.15` designs observers by this route). The null space of W_o (`unobservable_subspace`) is the set of states the sensors can never distinguish from zero.

### 2.3 The PBH test: which mode is lost

**(A, B) is uncontrollable if and only if there is an eigenvalue λ of A with rank [λI − A, B] < n.** If such a λ exists, there is a row vector w ≠ 0 with w(λI − A) = 0 and wB = 0: a left eigenvector that the input cannot touch. Then w x(t) = e^{λt} w x₀ whatever u does. The mode with eigenvalue λ evolves on its own. Conversely, if the system is uncontrollable, the left null space of W_c is invariant under Aᵀ and therefore contains a left eigenvector. `pbh_uncontrollable_modes` tests each distinct eigenvalue and returns the bad ones. For observability, rank [λI − A; C] < n (the dual). PBH and the rank test are logically equivalent, but numerically they are different computations, which is why the tests check that they agree.

### 2.4 The reaction-wheel satellite

Linearised about an inertial hold with small rates (gyroscopic terms dropped), each axis is J_iθ̈_i = τ_i. With x = (θ_x, θ_y, θ_z, ω_x, ω_y, ω_z) and wheel torques u along the spin-axis columns of W (3 × n_w):

A = [[0, I], [0, 0]],  B = [[0], [J⁻¹W]].   (3)

A has the eigenvalue 0 six times. W_c = [B, AB, 0, …] = [[0, J⁻¹W], [J⁻¹W, 0]] has rank 2 rank W. The satellite is **controllable if and only if the working wheel axes span three dimensions**:

- **three orthogonal wheels**: rank W = 3. With one failed, rank W = 2, rank W_c = 4: the angle and rate about the failed axis form a double integrator nobody can push, and PBH reports λ = 0 as uncontrollable.
- **four-wheel pyramid**: spin axes at β = 54.74° from +z (cos²β = 1/3), azimuths 0°, 90°, 180°, 270°. Nominally WWᵀ = (4/3)I, isotropic. With any one wheel failed, the remaining three axes still span 3D (rank 6); with two failed, rank W = 2 and rank W_c = 4 whichever two they are.

The eigenvalues of WWᵀ measure torque authority per direction: 1, 1, 1 for the orthogonal set; 4/3 three times for the pyramid; and **1/3, 4/3, 4/3** for the pyramid with one wheel failed. Every direction can still be pushed, but the weakest only with a quarter of the nominal authority.

### 2.5 Sensors

A star tracker measures attitude angles; gyros measure rates. With C selecting the measured states:

- a three-axis tracker alone, or with one gyro failed: W_o has rank 6. Rates are derivatives of angles, so the observer can recover them (`P09.15`).
- **gyros only**: rank 3. The angles are integrals of the rates with unknown constants. The unobservable subspace is all three angle directions, and PBH reports λ = 0.
- **a two-axis tracker** measuring only θ_x and θ_y (a single tracker determines rotation about its boresight poorly; modelled here as not at all) **plus three gyros**: rank 5. The unobservable direction is (0, 0, 1, 0, 0, 0), the boresight angle θ_z. Gyros give its rate but never its value.

### 2.6 Gramians

The **controllability Gramian** over a horizon T is

**W(T) = ∫₀ᵀ e^{At}BBᵀe^{Aᵀt} dt**,   (4)

symmetric and positive semidefinite, and positive definite exactly when (A, B) is controllable. Differentiating e^{At}BBᵀe^{Aᵀt} and integrating from 0 to T gives the **differential Lyapunov identity** AW(T) + W(T)Aᵀ + BBᵀ = e^{AT}BBᵀe^{AᵀT}. When A is stable, the right side vanishes as T → ∞ and the infinite-horizon Gramian solves the **Lyapunov equation**

**AW + WAᵀ + BBᵀ = 0.**   (5)

`lyap` solves (5) by vectorisation: with column-stacking vec, vec(AX) = (I ⊗ A)vec X and vec(XAᵀ) = (A ⊗ I)vec X, so (I ⊗ A + A ⊗ I)vec W = −vec(BBᵀ), an n² × n² linear solve. For the stable short-period model of `P09.11` (input elevator) it gives W = [[2.16216, 2.59459], [2.59459, 22.57297]]. The satellite is not stable (eigenvalues at 0), so its infinite Gramian diverges and only W(T) exists. `ctrb_gramian_finite` uses **Van Loan's** block exponential: with H = [[−A, BBᵀ], [0, Aᵀ]], e^{HT} = [[F₁₁, F₁₂], [0, F₂₂]] and W(T) = F₂₂ᵀF₁₂. For a double integrator with B = (0, 1/J)ᵀ it reproduces the closed form W(T) = [[T³/3, T²/2], [T²/2, T]]/J². The **observability Gramian** is the dual, AᵀM + MA + CᵀC = 0.

### 2.7 Minimum-energy steering

Among all inputs that take x₀ to x_f in time T, the one with the least ∫₀ᵀ|u|²dt is

**u(t) = Bᵀe^{Aᵀ(T−t)} W(T)⁻¹ d,  d = x_f − e^{AT}x₀,  energy = dᵀW(T)⁻¹d.**   (6)

(It reaches x_f: substitute into the convolution integral and the Gramian appears. It is minimal: any other admissible u differs by a v with ∫e^{A(T−t)}Bv dt = 0, which is orthogonal to (6) in the L² inner product, so the energy only adds.) The inverse Gramian is literally the price of reaching a direction. Directions with small eigenvalues of W are expensive. For a rest-to-rest slew Δθ about one axis with inertia J, (6) gives

u(t) = (6JΔθ/T²)(1 − 2t/T),  energy = 12J²Δθ²/T³:   (7)

a linear torque ramp from +6JΔθ/T² to −6JΔθ/T². A bang-bang slew (±4JΔθ/T² for half the time each) needs energy 16J²Δθ²/T³, **4/3 times as much**, but a lower peak torque. The peak of (7) must stay below the wheel limit: T ≥ √(6JΔθ/τ_max).

### 2.8 Where the model stops being valid

- **Gyroscopic coupling.** With stored wheel momentum or body rates, the axes couple through ω × (Jω + h). Coupling can restore controllability that the decoupled model denies (extension 1), and it can also complicate it.
- **Momentum and torque limits.** Controllability says nothing about saturation. A pyramid with a failed wheel is controllable, but the weakest direction has a quarter of the nominal authority and the wheels saturate sooner.
- **Numerical rank.** A wheel degraded to 10⁻⁶ of its torque still passes the rank test. The Gramian's smallest eigenvalue and the condition number of W_c tell the truth.
- **Boresight observability.** Real single trackers do measure roll about the boresight, typically several times less accurately than the cross-boresight angles. The model's "unobservable" is "poorly observable" in practice; the observability Gramian quantifies it.

## 3. Build the computational model

Open `starter/main.py`. `expm`, `numerical_rank` (SVD with a relative tolerance), `satellite_inertia`, `wheel_axes` and `short_period_model` are given.

### Step 1: rank tests and PBH

`ctrb` stacks B, AB, A²B, … horizontally; `obsv` stacks C, CA, … vertically. `is_controllable` and `is_observable` compare `numerical_rank` with n. `pbh_uncontrollable_modes` merges repeated eigenvalues, scales [λI − A, B] by the larger of ‖A‖ and ‖B‖, and counts singular values above the tolerance. `pbh_unobservable_modes` calls it on (Aᵀ, Cᵀ). `unobservable_subspace` returns the right singular vectors of W_o beyond its rank. `satellite_model` is (3).

### Step 2: Gramians and minimum energy

```python
def lyap(A, Q):
    n = A.shape[0]
    K = np.kron(np.eye(n), A) + np.kron(A, np.eye(n))
    X = np.linalg.solve(K, -Q.reshape(-1, order="F")).reshape((n, n), order="F")
    return 0.5 * (X + X.T)
```

`ctrb_gramian` checks stability, then calls `lyap(A, B @ B.T)`. `obsv_gramian` uses the dual. `ctrb_gramian_finite` builds H and exponentiates it. `min_energy_input` returns a function u(t) per (6) and the energy, raising if W(T) is singular.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-14-controllability-observability`. The tests:

1. **Structure**: the third and fourth blocks of W_c are A²B and A³B; W_o = W_c(Aᵀ, Cᵀ)ᵀ.
2. **Independent check, part 1**: on six random 5-state systems (three generic, three built as T·blockdiag·T⁻¹ with an untouched block), the rank test and PBH give the same verdict, and it is the right one; PBH names −3 in a 2 × 2 example; ranks survive x = Tz, and W_c transforms as T⁻¹W_c.
3. **Satellite**: orthogonal wheels lose rank 2 with one failure, PBH says λ = 0; the pyramid survives each single failure and fails each tested double one; three sensor sets give the verdicts of §2.5, and the unobservable direction of the two-axis tracker is ±e₃.
4. **Independent check, part 2**: `lyap` against `scipy.linalg.solve_continuous_lyapunov` for n = 2, 4, 6, with residual below 10⁻⁹; the short-period Gramian equals Simpson quadrature of (4) to 10 s (12 s horizon) within 10⁻⁷ relative; the pyramid's W(8 s) equals quadrature within 10⁻⁵ and satisfies the differential Lyapunov identity to 10⁻¹⁰; the double-integrator closed form; the infinite Gramian refuses the unstable satellite.
5. **Minimum energy**: with one pyramid wheel failed, (6) takes the satellite to a 5°, −3°, 8° attitude at rest in 30 s in `solve_ivp` to 10⁻⁸, and the quadrature of |u|² equals dᵀW⁻¹d; the single-axis slew equals (7), bang-bang costs 4/3 as much, and halving T costs eight times the energy.

### Compare

The oracles are `scipy.linalg.solve_continuous_lyapunov` (a Bartels–Stewart Schur-form solver, O(n³), against our O(n⁶) Kronecker solve) and quadrature of the defining integrals with `scipy.linalg.expm`. Running `main()` prints the failure table:

```
  orthogonal failed ()     : rank ctrb = 6/6, PBH uncontrollable eigenvalues [], ... lambda_min W(20 s) = 1.989e-01
  orthogonal failed (2,)   : rank ctrb = 4/6, PBH uncontrollable eigenvalues [0j], ... lambda_min W(20 s) = 0.000e+00
  pyramid    failed ()     : rank ctrb = 6/6, PBH uncontrollable eigenvalues [], ... lambda_min W(20 s) = 2.652e-01
  pyramid    failed (0,)   : rank ctrb = 6/6, PBH uncontrollable eigenvalues [], ... lambda_min W(20 s) = 1.191e-01
  pyramid    failed (3,)   : rank ctrb = 6/6, PBH uncontrollable eigenvalues [], ... lambda_min W(20 s) = 8.082e-02
  pyramid    failed (0, 2) : rank ctrb = 4/6, PBH uncontrollable eigenvalues [0j], ... lambda_min W(20 s) = 2.102e-45
```

and the short-period Gramian agrees with a trapezoid quadrature at 1 ms steps to 6 × 10⁻⁵ (the trapezoid's own error; the test uses Simpson). The Kronecker solve is fine for the n ≤ 20 of this phase. For a large structural model, the Schur method is the one to use.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Failure cases.** For orthogonal and pyramid wheel sets with no failure, each single failure and two double failures, record the rank of W_c, the PBH eigenvalues, the eigenvalues of WWᵀ and the smallest eigenvalue of W(20 s). Predict which configurations are controllable and order the survivors by the Gramian number.
2. **Sensor cases.** For the five sensor sets of `main()`, record the rank of W_o, the PBH eigenvalues and the unobservable directions. Predict them from physics before running.
3. **Slew time.** For a 10° rest-to-rest slew about z (J = 3 kg m²), compute the minimum-energy peak torque and energy for T = 10, 20 and 40 s, and the bang-bang energy. Predict the scaling with T and the shortest slew within 10 mN m.

Raw data format: CSV tables `config, failed, rank_ctrb, pbh_modes, eig_WWT, lambda_min_W20`; `sensors, rank_obsv, pbh_modes, unobservable_directions`; `T_s, peak_torque_mNm, energy, bang_bang_energy`, each with a header stating the model and that the values are computed, not measured.

What you should find: rank 6 for the orthogonal set and every single-failure pyramid; rank 4 with PBH λ = 0 for an orthogonal failure and for two pyramid failures. The smallest Gramian eigenvalue is 0.199 for the orthogonal set, 0.265 for the full pyramid, 0.119 with pyramid wheel 0 failed and 0.081 with wheel 3 failed, against 0 (to rounding: 2 × 10⁻⁴⁵) for the uncontrollable cases. The two single failures differ because the inertia is not isotropic. Sensors: a three-axis tracker gives rank 6 with or without the x gyro; a two-axis tracker plus gyros gives rank 5 with the boresight angle unobservable; gyros alone give rank 3 (all angles); a two-axis tracker alone gives rank 4 (θ_z and ω_z). Slew: peak torque 31.4, 7.85 and 1.96 mN m and energy 3.29 × 10⁻³, 4.11 × 10⁻⁴ and 5.14 × 10⁻⁵ N² m² s for T = 10, 20, 40 s (the 1/T² and 1/T³ laws), bang-bang 4/3 of that, and a shortest minimum-energy slew of **17.7 s** within the 10 mN m limit.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference |
|---|---|---|---|---|
| rank, orthogonal with one failure | 2 × rank W = 4 | 4 | PBH: λ = 0 | — |
| pyramid WWᵀ, one failure | smallest eigenvalue 1/3 | 0.3333 | — | — |
| two-axis tracker unobservable direction | boresight angle e₃ | (0, 0, −1, 0, 0, 0) | W_o N = 0 (test) | sign |
| short-period W | (5) | [[2.16216, 2.59459], [2.59459, 22.57297]] | quadrature: 6 × 10⁻⁵ (trapezoid), < 10⁻⁷ rel (Simpson, test) | quadrature error |
| W(T), double integrator | [[T³/3, T²/2], [T²/2, T]]/J² | same to 10⁻¹² | — | — |
| slew peak, T = 20 s | 6JΔθ/T² = 7.854 mN m | 7.854 mN m | — | — |
| slew energy, T = 20 s | 12J²Δθ²/T³ = 4.112 × 10⁻⁴ | 4.1123 × 10⁻⁴ | quadrature of |u|² (test) | — |
| shortest slew in 10 mN m | √(6JΔθ/τ_max) = 17.72 s | 17.72 s | — | — |

The verdicts agree three ways (rank test, PBH, physical reasoning), and the Gramians agree with direct quadrature. The engineering lesson is in the Gramian column. Controllability is a cliff (rank 6 or not), but the smallest Gramian eigenvalue slides: 0.265 for the full pyramid, 0.119 or 0.081 with one wheel gone, 0 with two. It is the number to watch as a wheel degrades. These are model results. In flight, wheel friction, momentum build-up and sensor noise all add, and the decoupled model ignores the gyroscopic terms.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), PBH, (5), (6) and (7); `controllability.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a wheel arrangement and a sensor suite for a small satellite that must survive any single failure, with the rank verdicts, the Gramian numbers and the slew time you would accept.

### Adversarial self-checks

- **What if a wheel is degraded rather than failed (10⁻⁶ of its torque)?** For the orthogonal set with the z wheel at 10⁻⁶ of its torque, the rank test still says 6 with a tolerance of 10⁻⁹. The smallest eigenvalue of W(20 s) falls from 0.199 to 5.5 × 10⁻¹³ (torque squared: a factor of about 10¹²), and the energy to reach that direction rises by the same factor. Which number would you put on the operator's display?
- **What if the inertia were isotropic?** All single-failure pyramid cases give the same smallest Gramian eigenvalue (0.1036 with J = 4I, computed for wheels 0, 1 and 3). The difference in the lesson (0.119 against 0.081) is the inertia, not the geometry.
- **Which assumption fails first in the real world?** Zero stored momentum. A momentum-biased spacecraft has gyroscopic coupling, and its controllability with a missing wheel must be re-tested on the coupled model (extension 1).
- **What measurement could falsify a verdict?** Command a slow sinusoidal torque on the remaining wheels and log the attitude. An axis that the model says is uncontrollable must show no response beyond disturbances.
- **Does the result survive a change of coordinates or solver?** Ranks are invariant under x = Tz (tested). The Gramian changes as T⁻¹WT⁻ᵀ under x = Tz, so "smallest eigenvalue" depends on the units of the states. Compare configurations in the same coordinates.
- **Which decision changes if the wheel limit is halved?** The shortest minimum-energy slew grows by √2 to 25.1 s; a bang-bang slew, whose peak is 2/3 of the minimum-energy peak, reaches the same angle in 20.5 s at the cost of 4/3 the energy.

## Check yourself

1. Define controllability and write the rank test.
2. What does the PBH test tell you that the rank test does not?
3. Why is observability of (A, C) the same question as controllability of (Aᵀ, Cᵀ)?
4. Three orthogonal reaction wheels, one fails. What is rank W_c, and which state pair is uncontrollable?
5. Why can gyros alone never determine attitude, however accurate they are?
6. What equation does the infinite-horizon controllability Gramian satisfy, and why does it not exist for the free rigid-body satellite?
7. A 10° slew about an axis with J = 3 kg m² must be done in 20 s with minimum energy. What is the peak torque?
8. A wheel degrades to one millionth of its torque. What does the rank test say, and what number tells you the truth?
