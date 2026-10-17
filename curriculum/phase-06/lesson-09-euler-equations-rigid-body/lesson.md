# P06.09  3D Rigid Body Dynamics: Euler's Equations

> A free body's rates wander; its angular momentum in space does not move. Two conserved quantities fix the path, and an elliptic function gives the timetable.

| | |
|---|---|
| **Learning outcomes** | 1. Find the principal moments and axes of an inertia tensor with products of inertia as a symmetric eigenproblem, and write Euler's equations in principal and non-principal axes. 2. Implement and test a tumbling-body simulator (Euler's equations with quaternion attitude) that conserves the kinetic energy, \|H\| and the inertial angular-momentum vector to 10⁻⁹, and holds the rate vector on Poinsot's invariable plane. 3. Derive the energy and momentum ellipsoids, classify the polhodes by the sign of H² − 2TI₂, and locate the separatrix through the intermediate axis. 4. Derive and implement the exact torque-free solution in Jacobi elliptic functions, check it against the simulation, and predict the flip interval 2K(m)/λ of the Dzhanibekov effect. 5. Explain why the flip interval grows only logarithmically as the disturbance shrinks, and what that means for spinning a spacecraft about its intermediate axis. |
| **Prerequisites** | `P04.04` (inertia tensor, Euler's equations, the linearised intermediate-axis instability), `P06.08` (quaternion attitude), `P06.07` (DCMs), `P03.05` (Runge–Kutta methods), `P01.16` (eigenvalues of symmetric matrices), `P05.05` (moments of inertia), `P06.04` (the transport theorem). Skill check: you can say why a real symmetric matrix has orthogonal eigenvectors, and you have met sn, cn and dn or are willing to treat them as "sine-like functions that scipy evaluates". |
| **Engineering question** | A T-handle spun at 4 rev/s about its intermediate axis with a 1 % wobble flips end over end. How often, and how does the interval change if the wobble is a thousand times smaller? And what does the same mathematics say about a 500 kg box satellite spun at 0.5 rad/s about its intermediate axis, or a tumbling piece of debris? |
| **Assumptions and limits** | Rigid body; no external torque except where stated; distinct principal moments for the elliptic solution (a symmetric body is `P06.10`); principal axes known exactly. The rigid model conserves energy exactly; a real body with flexible parts or liquid dissipates it, which changes the long-term answer completely (`P06.10`). Near the separatrix the motion is extremely sensitive to the initial state, so a predicted flip time is only as good as the measured disturbance. |
| **Analytical work** | Principal axes as eigenvectors; Euler's equations in any axes and their covariance; the power balance; conservation of T, \|H\| and H_i; the two ellipsoids and the polhode classification; the separatrix planes; the invariable plane ω · H = 2T; the Jacobi elliptic solution (amplitudes, rate λ, parameter m), its period 4K/λ and flip interval 2K/λ; the small-amplitude limit (the `P04.04` wobble frequency) and the near-separatrix limit (2/σ) ln(4/k′). |
| **Computational deliverable** | `starter/main.py`: `principal_axes`, `kinetic_energy`, `angular_momentum_body`, `angular_momentum_inertial`, `euler_rhs`, `simulate`, `invariant_drift`, `invariable_plane_residual`, `polhode_type`, `elliptic_parameters`, `elliptic_solution`, `small_wobble_frequency`, `intermediate_growth_rate` (`quat_multiply`, `dcm_from_quat`, `cross3`, `Tumble`, `flip_times` and `main` are given). Tests prove: principal axes form a proper eigenbasis; Euler's equations in component form, for a sphere, in rotated axes, and with dT/dt = ω·τ; T, \|H\| and H_i conserved to 10⁻⁹; the invariable plane; dH_i/dt equal to an inertial torque; non-principal axes give the same motion; the polhode classification including an exact separatrix; the elliptic solution equal to the simulation to 10⁻⁸ for five initial states on both sides of the separatrix, for any axis order, satisfying Euler's equations; the flip interval equal to the simulated flips to 10⁻⁴; the period tending to the wobble frequency; 2 ln 10/σ per decade of disturbance near the separatrix. |
| **Practical activity** | **computational experiment**: the Dzhanibekov T-handle simulated and predicted in closed form; the flip interval against the disturbance over five decades; polhodes on both sides of the separatrix; a tolerance study of the conserved quantities. |
| **Equipment** | None. A computer with the environment from `P00.01`. (A book with a rubber band around it, tossed spinning about each of its three axes, shows the effect; nothing is measured with it.) |
| **Safety** | None beyond normal computer use. If you toss a book, toss it gently over a bed or sofa. |
| **Evidence** | Completed `main.py`, passing tests, `tumbling_body.png`, the derivations of the conservation laws, the polhode classification and the elliptic solution's period, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the conservation laws and the polhode classification are derived; the elliptic solution is verified by substitution into Euler's equations; all public tests pass; the simulated flip interval is compared with 2K/λ with numbers. |
| **References** | `@mit-ocw-16-07-dynamics` L26 "3D rigid body dynamics: the inertia tensor" (Angular Momentum; Principal Axes of Inertia); L27 "3D rigid body dynamics: kinetic energy, instability, equations of motion" (Governing Equations for Rotational Motion of a Three-Dimensional Body); L28 "3D rigid body dynamics: equations of motion; Euler's equations" (Body-Fixed Axis; Stability of Free Motion about a Principal Axis) (required-free). `@libretexts-cline-variational-principles` Ch. 13 "Rigid-body Rotation": 13.4 Inertia Tensor; 13.6 Principal Axis System; 13.11 Angular Momentum and Angular Velocity Vectors; 13.12 Kinetic Energy of Rotating Rigid Body; 13.17 Euler's equations of motion for rigid-body rotation; 13.21 Torque-free rotation of an asymmetric rigid rotor; 13.22 Stability of torque-free rotation of an asymmetric body (required-free). `@mit-ocw-2-003sc-engineering-dynamics` Lecture 10: Equations of Motion, Torque, Angular Momentum of Rigid Bodies; Lecture 13: Four Classes of Problems With Rotational Motion (supplementary-free). `@mit-ocw-16-61-aerospace-dynamics` lecture 12 "Rigid Body Dynamics" (supplementary-free). `@openstax-university-physics-vol1` 11.2 Angular Momentum; 11.3 Conservation of Angular Momentum (required-free, the fixed-axis background). `@libretexts-baker-haynes-engineering-statics` 10.7 Products of Inertia (required-free). `@mit-ocw-18-06-linear-algebra` Lecture 25: Symmetric matrices and positive definiteness (required-free). `@ocw-16-90` Unit 1 "1.9 Runge-Kutta Methods"; `@scipy-userguide` solve_ivp (required-free). |
| **Extensions** | ★ Draw the herpolhode: the path of the rate vector's tip on the invariable plane in inertial axes, and show that it is not closed in general. ★ Add a small constant body torque about the intermediate axis and find how it shifts the flip times. ★ Compute the attitude in closed form too (Jacobi's solution for the Euler angles involves an elliptic integral of the third kind) and compare with the quaternion integration. |
| **Time** | unmeasured |

## 1. The engineering problem

Most of what a spacecraft attitude engineer does assumes that the body is turning slowly about a well-behaved axis. But satellites tumble: a spent upper stage after separation, a satellite that has lost control, a piece of debris that a servicing spacecraft must approach and grab. To plan an approach, or to understand why a stage sprayed its propellant in an unexpected direction, one needs to know how a free rigid body turns when it is not spinning about a principal axis. `P04.04` derived Euler's equations, showed that torque-free motion conserves the inertial angular momentum and the energy, and linearised the instability of spin about the intermediate axis. This lesson takes the next step: the full nonlinear motion, its geometry, and its exact solution.

The memorable example is the **Dzhanibekov effect**, usually named after the cosmonaut who reported seeing a wing nut, spun off its thread in weightlessness, flip over every few seconds. Spin a T-handle (or a tennis racket, or a book) about its intermediate axis and it does the same: it spins steadily, flips 180° end over end, spins steadily the other way, flips back, again and again. No torque acts. The angular momentum is fixed in space the whole time. The flip is exact rigid-body dynamics, and its timetable follows from two conservation laws and an elliptic function. This lesson's T-handle has illustrative principal moments 2.2, 6.0 and 7.4 × 10⁻⁶ kg·m², is spun at 4 rev/s (25.133 rad/s) about axis 2 with a 1 % disturbance on axis 1, and flips every **0.8125 s**, by simulation and by closed form alike.

Success criteria, written first: T, |H| and every inertial component of H conserved to 10⁻⁹ of |H|; the simulated body rates equal to the elliptic solution to 10⁻⁸ for generic states; the simulated flip interval equal to 2K/λ to 10⁻⁴; the closed form's two limits (small wobble, near separatrix) reproduced.

## 2. Understand and derive

### 2.1 Principal axes

The inertia tensor about the centre of mass, I = Σmᵢ(|rᵢ|²E − rᵢrᵢᵀ) (`P04.04`, 16.07 L26, Cline 13.4), is symmetric and positive definite. Its off-diagonal entries, the **products of inertia** (Baker & Haynes 10.7), are not zero in arbitrary body axes: a satellite with an instrument boom off to one side, or a body whose reference axes were chosen for manufacturing convenience. Because I is real and symmetric, it has three real positive eigenvalues and three orthogonal eigenvectors (18.06 Lecture 25):

**I = R diag(I₁, I₂, I₃) Rᵀ,**   (1)

with the columns of R the **principal axes** (made a proper rotation by flipping one column if det R = −1) and I₁ ≤ I₂ ≤ I₃ the **principal moments** (Cline 13.6). In principal axes L = Iω has components (I₁ω₁, I₂ω₂, I₃ω₃). The lesson's code takes the box satellite of `P04.04`, I = diag(135.4, 208.3, 260.4) kg·m², describes it in misaligned axes (full tensor with products as large as 54 kg·m²) and recovers the three moments to 10⁻¹³ with `numpy.linalg.eigh`.

### 2.2 Euler's equations in any axes

From dH/dt = τ in the inertial frame and the transport theorem (`P06.04`; 16.07 L27–L28; Cline 13.17):

**I ω̇ + ω × (I ω) = τ**   (2)

in any body-fixed axes, with I constant in those axes. It is **covariant**: in axes rotated by a constant R, I′ = RIRᵀ, ω′ = Rω, τ′ = Rτ satisfy the same equation (a test). In principal axes,

I₁ω̇₁ = (I₂ − I₃)ω₂ω₃ + τ₁,  I₂ω̇₂ = (I₃ − I₁)ω₃ω₁ + τ₂,  I₃ω̇₃ = (I₁ − I₂)ω₁ω₂ + τ₃.   (3)

Dot (2) with ω: ω · Iω̇ = ω · τ because ω · (ω × Iω) = 0, and ω · Iω̇ = d(½ω · Iω)/dt, so **dT/dt = ω · τ**, the power balance. The attitude is carried by a quaternion with q̇ = ½ q ⊗ (0, ω) (`P06.08`), and the inertial angular momentum is H_i = Cᵀ(Iω).

### 2.3 Two ellipsoids and the polhodes

With τ = 0 two scalars are conserved (Cline 13.11, 13.12):

**2T = I₁ω₁² + I₂ω₂² + I₃ω₃²,   H² = I₁²ω₁² + I₂²ω₂² + I₃²ω₃².**   (4)

In body-rate space each is an ellipsoid; ω lies on both, so it moves along their intersection, the **polhode**. Multiply the first by I₂ and subtract it from the second:

**H² − 2T I₂ = I₁(I₁ − I₂)ω₁² + I₃(I₃ − I₂)ω₃².**   (5)

The first coefficient is negative, the second positive. If H² > 2TI₂, ω₃ can never be zero: the polhode circles the **major axis** (axis 3). If H² < 2TI₂, ω₁ can never be zero: it circles the **minor axis**. If H² = 2TI₂ exactly, (5) is a pair of planes through the ω₂ axis, I₁(I₂ − I₁)ω₁² = I₃(I₃ − I₂)ω₃², and the polhodes are the two ellipses where these planes cut the energy ellipsoid: the **separatrix**, which passes through pure spin about the intermediate axis. A small disturbance of intermediate-axis spin puts the state just on one side of it, on a polhode that runs close to the separatrix all the way round to the opposite pole −ω₂: that is the flip.

**Poinsot's invariable plane.** In the inertial frame, ω_i · H_i = ω · Iω = 2T is constant and H_i is fixed, so the tip of ω_i stays in a fixed plane perpendicular to H at distance 2T/|H| from the origin. The body moves as if its energy ellipsoid rolled without slipping on that plane. The test measures the plane's constancy to 10⁻⁹.

### 2.4 The exact solution

Equations (3) with τ = 0 can be solved in closed form (Cline 13.21 treats the torque-free asymmetric rotor). Order the axes so that I₁ < I₂ < I₃, and take H² > 2TI₂ (around the major axis). With the shorthand a = 2TI₃ − H² ≥ 0 and b = H² − 2TI₁ ≥ 0, try ω₁ = A₁ cn τ, ω₂ = A₂ sn τ, ω₃ = A₃ dn τ with τ = λt + τ₀, where sn, cn and dn are the **Jacobi elliptic functions** of parameter m. They satisfy sn² + cn² = 1, dn² + m sn² = 1, and

d sn/dτ = cn dn,   d cn/dτ = −sn dn,   d dn/dτ = −m sn cn.   (6)

Substituting into (3) and matching (Problem 2):

**A₁² = a / (I₁(I₃ − I₁)),  A₂² = a / (I₂(I₃ − I₂)),  A₃² = b / (I₃(I₃ − I₁)),**
**λ² = (I₃ − I₂) b / (I₁I₂I₃),  m = (I₂ − I₁) a / ((I₃ − I₂) b).**   (7)

For H² < 2TI₂ the roles of axes 1 and 3 swap (dn on axis 1, cn on axis 3, with the matching changes in A, λ and m). The signs of the three amplitudes and the starting phase τ₀ come from ω(0) and from (3) itself (the code tries the few sign choices and keeps the one whose derivative matches (3)). With m = 0 the functions become sin, cos and 1, and with m = 1 they become tanh, sech and sech: the separatrix, where the period is infinite.

**Period and flip interval.** sn and cn have period 4K(m) in τ, with K(m) the complete elliptic integral of the first kind; dn has period 2K. So the body rates repeat with

**period = 4K(m)/λ,   flip interval = 2K(m)/λ,**   (8)

the second being the time between sign reversals of ω₂, i.e. between flips of a body spun near its intermediate axis. For the T-handle: m = 0.99986383, λ = 14.3681 rad/s, K = 5.8373, so the period is **1.6251 s** and the flip interval **0.8125 s**.

**Two limits.** *Small wobble about an extreme axis:* spin ω₃ with a small ω₂. Then a ≈ 0, m → 0, K → π/2, and 4K/λ → 2π/λ with λ² → (I₃ − I₂)(I₃ − I₁)ω₃²/(I₁I₂), exactly `P04.04`'s wobble frequency. For the T-handle spun about axis 3 with a 1 % disturbance the period is 0.33663 s against 2π/(wobble frequency) = 0.33664 s. *Near the separatrix:* with k′² = 1 − m → 0, K ≈ ln(4/k′), and λ tends to the intermediate-axis growth rate σ = |ω₂|√((I₃ − I₂)(I₂ − I₁)/(I₁I₃)) of `P04.04`. Since k′ is proportional to the disturbance ε,

**flip interval ≈ (2/σ) ln(4/k′):** each tenfold reduction of the disturbance adds only 2 ln 10/σ.   (9)

A near-perfect intermediate-axis spin still flips; it merely waits a little longer. Computing 1 − m directly as (I₃ − I₁)(H² − 2TI₂)/((I₃ − I₂)b), rather than as 1 minus (7), keeps it accurate down to 10⁻¹²; SciPy's `ellipkm1` then evaluates K from 1 − m without cancellation.

### 2.5 Where the model fails

- **Rigidity.** A real body dissipates energy (flexing, sloshing); at fixed |H| that lowers T and moves the polhode towards the major axis, so a tumble about the minor axis or near the separatrix does not persist. `P06.10` builds that.
- **No torque.** Gravity gradient, drag and magnetic torques act on real satellites at 10⁻⁶–10⁻³ N·m; over many periods they change H.
- **Exact initial state.** Near the separatrix the flip timing depends on the logarithm of an initial disturbance that is never known well. The *existence* and *shape* of the flip are robust; its *phase* is not.
- **Distinct moments.** If two moments are equal, (7) degenerates (the motion is the simpler coning of `P06.07` and `P06.10`); if they are nearly equal, the intermediate axis is barely distinguishable.

## 3. Build the computational model

Open `starter/main.py`. The quaternion helpers from `P06.08`, `cross3`, the result type `Tumble`, `flip_times` and `main` are given.

### Step 1: inertia, energy, momentum

`principal_axes` is `np.linalg.eigh` with the determinant fix. `kinetic_energy`, `angular_momentum_body` and `angular_momentum_inertial` accept single states or stacks.

### Step 2: the simulator

```python
def euler_rhs(I, w, torque=None):
    tau = np.zeros(3) if torque is None else np.asarray(torque, dtype=float)
    return np.linalg.solve(I, tau - cross3(w, I @ w))
```

`simulate` integrates [q, ω] with DOP853 (rtol 10⁻¹¹) and a torque function of (t, q, ω) in body axes; the quaternion is renormalised on output. `invariant_drift` reports the three conservation residuals, `invariable_plane_residual` Poinsot's.

### Step 3: the closed form

`polhode_type` is the sign of (5). `elliptic_parameters` sorts the axes with a proper signed permutation, computes 2T, H², a, b and H² − 2TI₂ term by term from (4)–(5) so that nothing cancels, builds (7), finds the signs and τ₀ (`scipy.special.ellipkinc` gives τ₀ from the starting amplitude φ₀ = atan2(sn, cn)), and returns the period and flip interval (8). `elliptic_solution` evaluates `scipy.special.ellipj` and maps back. `small_wobble_frequency` and `intermediate_growth_rate` are the `P04.04` limits.

Running `main()` prints the principal axes of the misaligned satellite, the T-handle's closed form and simulation, the flip-interval table of Experiment 1 and the stable-axis check, and writes `tumbling_body.png`: the T-handle's body rates with the elliptic ω₂ on top; polhodes on both sides of the separatrix; the inertial H drift.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-09-euler-equations-rigid-body`. The tests:

1. **Principal axes**: eigenvalues recovered, RᵀIR diagonal, R a proper rotation.
2. **Euler's equations**: the component form (3); no gyroscopic term for a sphere; covariance under a rotation of axes; dT/dt = ω · τ.
3. **Known states**: energy, H_b and H_i for single and stacked states.
4. **Conservation**: T, |H| and H_i to 10⁻⁹ over 20 s while ω₁ changes by more than 0.1 rad/s; |q| = 1; the invariable plane.
5. **Torque**: H_i(t) = H_i(0) + τ_i t for an inertially fixed torque.
6. **Non-principal axes**: the same motion, rotated.
7. **Polhodes**: major, minor and separatrix (including a state exactly on the separatrix planes); the closed form refuses the separatrix.
8. **Elliptic solution**: equal to the simulation to 10⁻⁸ for five initial states with mixed signs on both sides; invariant under relabelling the axes; satisfying (3) by finite differences; conserving T and |H| to 10⁻¹².
9. **Flips**: the simulated T-handle flips at intervals equal to 2K/λ to 10⁻⁴.
10. **Limits**: the period tends to 2π over the wobble frequency; near the separatrix λ → σ and each decade of disturbance adds 2 ln 10/σ.
11. `main` writes the figure.

### Compare

The independent checks are the conservation laws and the closed form. Over three periods (4.88 s, 6001 samples) the T-handle simulation conserves T to **3.5 × 10⁻¹¹**, |H| to **1.8 × 10⁻¹¹** and the inertial H vector to **5.9 × 10⁻¹¹** of |H|; the rate vector stays on the invariable plane to **1.1 × 10⁻¹⁵**. Its flips come at 0.4063, 1.2188, 2.0313, 2.8439, 3.6564 and 4.4689 s, a mean interval of **0.81253 s** against 2K/λ = **0.81253 s**, and the body rates differ from the elliptic solution by at most **1.1 × 10⁻¹⁰** of |ω₀|. With the disturbance on axis 3 instead (the major-axis side of the separatrix), the agreement is **5.7 × 10⁻¹¹**. These checks catch the classic mistakes: a sign error in one of (3) breaks energy conservation at once; using I instead of I⁻¹ (or forgetting `solve`) in non-principal axes breaks the covariance test; and computing 1 − m by subtraction loses digits near the separatrix: at ε = 10⁻⁶ the flip interval moves by 1.4 × 10⁻⁷ relative, and once 1 − m falls below about 10⁻¹⁶ it rounds to zero and K becomes infinite.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements. Write each prediction down before running it.

1. **The flip interval against the disturbance.** Spin the T-handle at 4 rev/s about axis 2 with a disturbance ε × spin on axis 1, for ε = 10⁻¹, 10⁻², 10⁻³, 10⁻⁴ and 10⁻⁶. Predict from (9) how much each decade adds (2 ln 10/σ with σ = 14.367 s⁻¹: 0.3205 s). *Finding:*

   | ε | 1 − m | 2K/λ (s) | (2/λ) ln(4/k′) (s) |
   |---|---|---|---|
   | 10⁻¹ | 1.344 × 10⁻² | 0.4908 | 0.4896 |
   | 10⁻² | 1.362 × 10⁻⁴ | 0.8125 | 0.8125 |
   | 10⁻³ | 1.362 × 10⁻⁶ | 1.1331 | 1.1331 |
   | 10⁻⁴ | 1.362 × 10⁻⁸ | 1.4536 | 1.4536 |
   | 10⁻⁶ | 1.362 × 10⁻¹² | 2.0947 | 2.0947 |

   A disturbance a million times smaller than 1 % only multiplies the flip interval by 2.6. The logarithmic estimate is exact to four figures from ε = 10⁻² down.
2. **Both sides of the separatrix.** Move the 1 % disturbance from axis 1 to axis 3. Predict the polhode type of each. *Finding:* axis 1 gives a 'minor' polhode (period 1.6251 s), axis 3 a 'major' one (period 1.5952 s); both flip, because both run close to the separatrix, and both match the simulation to 10⁻¹⁰. The figure's polhodes show the two families and the separatrix between them.
3. **Stable spin.** Spin about axis 3 with 1 % on axis 2. Predict the period from `P04.04`'s wobble frequency. *Finding:* 0.33663 s by (8) against 0.33664 s from the linear wobble: the linear theory is exact to the amplitude correction.
4. **Solver tolerance.** Repeat a T-handle run (Problem 4) at rtol 10⁻⁶ and 10⁻¹¹. Predict which output degrades: the conserved quantities or the flip times. *Finding (Problem 4's case):* the invariants drift by 10⁻⁵ at rtol 10⁻⁶ and 10⁻¹¹ at 10⁻¹¹; the flip times agree to 10⁻⁴ s. Near the separatrix, a solver error acts like a change in the disturbance, which moves the flip times only logarithmically.

Raw data format: CSV tables `eps, one_minus_m, flip_closed_s, flip_log_s` (1); `disturbed_axis, polhode_type, period_s, max_dev` (2); `case, period_closed_s, period_linear_s` (3); `rtol, T_drift, H_drift, flip_times` (4); plus `tumbling_body.png`.

## 5. Validate and improve

| quantity | closed form / conservation law | simulation | difference | cause |
|---|---|---|---|---|
| T, \|H\|, H_i over 4.88 s | constant | 3.5, 1.8, 5.9 × 10⁻¹¹ | — | solver tolerance |
| invariable plane | ω_i · Ĥ = 2T/\|H\| | 1.1 × 10⁻¹⁵ | — | rounding |
| T-handle flip interval (s) | 2K/λ = 0.81253 | 0.81253 (mean of 5) | < 10⁻⁴ relative | sampling of the zero crossings |
| T-handle body rates | Jacobi solution | — | 1.1 × 10⁻¹⁰ of \|ω₀\| | solver tolerance |
| major-axis period (s) | 0.33663 (8) | 0.33664 (`P04.04` linear) | 3 × 10⁻⁵ | amplitude correction |
| flip interval at ε = 10⁻⁶ (s) | 2.0947 | (9): 2.0947 | < 10⁻⁴ | — |

The simulation and the closed form agree to the solver tolerance wherever both apply, and the closed form agrees with both of `P04.04`'s linear results in their limits. The engineering answer: **the T-handle at 4 rev/s with a 1 % wobble flips every 0.81 s; a wobble a thousand times smaller stretches that only to 1.13 s. For `P04.04`'s box satellite at 0.5 rad/s about its intermediate axis (σ = 0.1641 s⁻¹), a 10⁻⁶ rad/s disturbance gives the first reversal after 86.1 s and a flip every 172.2 s (Problem 5): intermediate-axis spin cannot be held passively, and making the initial spin more precise buys almost nothing.**

What to improve: the rigid-body assumption. A body that dissipates energy drifts off these closed polhodes; `P06.10` adds the dissipation and shows where it goes.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)'s power balance, (4), (5), the substitution of (7) into (3), and the two limits of (8); `tumbling_body.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the spin doubles?** λ and σ double, so every time in (8) halves: the T-handle flips every 0.406 s; the number of flips per revolution is unchanged (the motion is self-similar in time).
- **What if I₂ → I₁ (a nearly symmetric body)?** σ → 0 and the instability becomes very slow; at I₂ = I₁ the body is symmetric and the closed form becomes the coning of `P06.10`. The elliptic code refuses equal moments.
- **Which assumption fails first in the real world?** Rigidity, for anything with fuel, flexible appendages or loose parts: energy dissipation changes which axis the body ends up spinning about.
- **What measurement could falsify the model?** Film a spinning T-handle or book in slow motion, measure the spin and the time between flips at two disturbance levels: (9) predicts the difference, independent of the exact disturbance, as 2 ln(ε₁/ε₂)/σ.
- **Does the result survive a change of solver?** Yes: rtol 10⁻⁶ and 10⁻¹¹ give the same flip times to 10⁻⁴ s (Experiment 4), and the closed form needs no solver.
- **Which decision changes if the uncertainty doubles?** If the moments are uncertain by 10 %, σ (and so every flip time) changes by up to about 10 %, but the conclusion that intermediate-axis spin is unstable does not change unless the uncertainty can swap the order of the moments.

## Check yourself

1. Why are the principal moments real and the principal axes orthogonal?
2. Which two quantities are conserved in torque-free motion, and which vector is conserved in which frame?
3. Write the condition that separates polhodes around the major axis from those around the minor axis.
4. What is the invariable plane, and how far is it from the origin?
5. In the elliptic solution, which axis carries the dn function, and why can its rate never change sign?
6. What is the period of the body-rate motion in terms of K(m) and λ, and what happens to it on the separatrix?
7. By how much does the flip interval of a body spun near its intermediate axis grow when the disturbance is made ten times smaller?
8. What does the small-amplitude limit of the elliptic period reduce to?
