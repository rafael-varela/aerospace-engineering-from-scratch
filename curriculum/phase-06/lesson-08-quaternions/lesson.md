# P06.08  Quaternions

> Four numbers, one constraint, no singularity: a quaternion is a rotation axis and half its angle, and it multiplies the way rotations compose.

| | |
|---|---|
| **Learning outcomes** | 1. Implement the quaternion algebra (Hamilton product, conjugate, norm, inverse) from its rules and prove that the norm is multiplicative. 2. Derive the DCM of a unit quaternion from q* v q, convert quaternion ↔ DCM (Shepperd's method) ↔ 3-2-1 Euler angles, and show that q and −q are the same attitude. 3. Derive the kinematic equation q̇ = ½ q ⊗ (0, ω), integrate it with `solve_ivp`, with RK4 and with the exact exponential step, and check every result against the P06.07 DCM route and a closed form. 4. Measure the norm drift of RK4 (one order above its attitude error), and explain what renormalisation fixes and what it does not. 5. Interpolate star-tracker attitude samples with SLERP, show its second-order convergence and its shortest-arc rule, and decide when normalised linear interpolation is good enough. |
| **Prerequisites** | `P06.07` (DCMs, Euler sequences, Ċ = −[ω×]C, coning), `P01.15` (Rodrigues, passive rotations), `P01.07` (dot and cross products), `P03.05` (Runge–Kutta), `P03.06` (adaptive integrators), `P04.04` (a quaternion used without derivation). Skill check: you can expand (a × b) × c and you know that sin²x + cos²x = 1 is all trigonometry ever needs here. |
| **Engineering question** | A spacecraft's star tracker reports its attitude as a quaternion a few times a second while the spacecraft cones at 2 rad/s. How should the attitude computer propagate the quaternion between reports, how often must it renormalise it, and how far apart may the reports be if the attitude between them is interpolated to 10⁻³ rad? |
| **Assumptions and limits** | Rigid body; exact body rates; a non-rotating reference frame; passive convention throughout (q describes the frame b relative to n, as C_b/n does); the star-tracker samples are exact (real trackers have noise of order 10⁻⁵–10⁻⁴ rad, which sets a floor the interpolation does not need to beat). Quaternions have no singular attitude; their only defect is the norm constraint, which integration breaks slowly. |
| **Analytical work** | The Hamilton rules; the product in scalar–vector form; the norm and the conjugate of a product; the unit quaternion of an axis and angle; v_b = q*(0, v_n)q and the DCM formula; the composition order; Shepperd's method; the 3-2-1 quaternion and its inverse; q̇ = ½ q ⊗ (0, ω) and q · q̇ = 0; the exact step for a constant rate; the RK4 norm defect; SLERP and its constant rate. |
| **Computational deliverable** | `starter/main.py`: `quat_multiply`, `quat_conjugate`, `quat_norm`, `quat_normalize`, `quat_inverse`, `quat_from_axis_angle`, `axis_angle_from_quat`, `transform`, `dcm_from_quat`, `quat_from_dcm`, `quat_from_euler_321`, `euler_321_from_quat`, `quat_rate`, `integrate_quat`, `integrate_quat_rk4`, `quat_exp_step`, `quat_error_angle`, `slerp`, `nlerp`, `interpolate_samples` (the DCM helpers carried over from `P06.07`, the coning closed form and `main` are given). Tests prove: the Hamilton rules; associativity and a multiplicative norm; conjugate and inverse; axis–angle both ways with q and −q equal; the quaternion DCM equal to Rodrigues and to SciPy; Shepperd in every branch including 180° turns; transformation and composition order; 3-2-1 both ways including gimbal lock; q · q̇ = 0 and agreement with Ċ = −[ω×]C; the integrated quaternion equal to the DCM route and to coning to 10⁻⁹; the exact step to 10⁻¹³; RK4 at order 4 with a norm defect at order 5; SLERP's endpoints, constant rate, great arc and shortest arc; star-tracker interpolation at second order. |
| **Practical activity** | **computational experiment**: quaternion against DCM propagation of a coning spacecraft (accuracy, drift, cost); an RK4 step-size study with and without renormalisation; star-tracker interpolation by SLERP and NLERP against the sample interval; the double cover in practice. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `quaternions.png`, the derivations of the DCM formula, the kinematic equation and SLERP, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: C(q) and q̇ are derived (not quoted); all public tests pass; the quaternion route is compared with the DCM route with numbers; the note gives a sample interval and a renormalisation policy with their reasons. |
| **References** | `@ntrs-nasa-tm-74839-euler-angles-quaternions` 2.2 Transformation Matrices Using the Hamilton Quaternion; 2.3 Euler Angle and Quaternion Relationships; Appendix A Relationships for the Three-Axis Euler Angle Rotation Sequences (required-free: the quaternion–matrix–Euler relationships for all twelve sequences, a check table for your code; **its matrices map the rotated frame into the original frame, the transpose of this lesson's C_b/n**, so compare with Cᵀ and check its quaternion sign convention before trusting a match). `@mit-ocw-18-06-linear-algebra` Lecture 17: Orthogonal matrices and Gram-Schmidt (required-free: the DCM side of every check). Context: `@nasa-basics-of-space-flight` Chapter 11 "Onboard Systems", Attitude and Articulation Control Subsystems (AACS), the "Celestial Reference" and "Inertial Reference" paragraphs (supplementary-free, qualitative: star trackers and inertial references). |
| **Extensions** | ★ Average N noisy star-tracker quaternions by the eigenvector method (the largest eigenvector of Σ qqᵀ) and show it is immune to sign flips. ★ Write the quaternion kinematic equation as q̇ = ½ Ω(ω) q with a 4 × 4 matrix Ω, show Ω is skew-symmetric, and derive the exact step as exp(½ΩΔt). ★ Implement the multiplicative error quaternion δq = q̂* ⊗ q and its small-angle form δq ≈ (1, δθ/2), the starting point of the attitude Kalman filter in Phase 18 (`P18.10`). |
| **Time** | unmeasured |

## 1. The engineering problem

`P06.07` ended with a choice. Three Euler angles are compact and readable, but singular somewhere; the DCM is never singular, but it carries nine numbers for three degrees of freedom and drifts away from being a rotation unless repaired. Spacecraft attitude systems almost all take a third way, the **unit quaternion**: four numbers, one constraint (|q| = 1), no singular attitude, and a composition rule that is a single product. Star trackers, the cameras that recognise star patterns and report the spacecraft's orientation, output quaternions; attitude filters estimate quaternions; ground software and telemetry carry them.

This lesson builds a quaternion library from its definition and tests every piece against the `P06.07` DCM route, which is carried into this lesson as a set of given helpers so that the two routes are independent implementations of the same physics. The test motion is the coning spacecraft of `P06.07`: spin φ̇ = 2.0 rad/s, precession ψ̇ = 0.4 rad/s, half-cone β = 20°, whose attitude is known in closed form. Integrating its body rates by quaternion and by DCM must give the same attitude to the solver tolerance (it does: 1.8 × 10⁻¹⁰ rad over 30 s). Three engineering questions follow. How does the quaternion's norm drift under fixed-step integration, and how often must it be repaired? How far apart may star-tracker samples be if the attitude between them is interpolated? And what goes wrong when software forgets that q and −q are the same attitude?

Success criteria, written first: the quaternion DCM equal to Rodrigues, to the 3-2-1 DCM and to SciPy to 10⁻¹⁴; the integrated quaternion equal to the DCM route and to the closed form to 10⁻⁹; RK4 at fourth order; SLERP interpolation converging at second order.

## 2. Understand and derive

### 2.1 The algebra

A quaternion is q = q₀ + q₁**i** + q₂**j** + q₃**k**, four real numbers, with Hamilton's rules

**i² = j² = k² = ijk = −1,**   (1)

from which ij = k = −ji, jk = i = −kj, ki = j = −ik (multiply ijk = −1 by i on the left, and so on). Write q = (q₀, **v**) with a scalar part and a vector part. Expanding a product term by term with (1) gives the scalar–vector form

**p ⊗ q = (p₀q₀ − p_v · q_v,  p₀ q_v + q₀ p_v + p_v × q_v).**   (2)

The cross product makes it **non-commutative**; it is associative (a test). The **conjugate** is q* = (q₀, −**v**), and from (2) q ⊗ q* = (|q|², 0) with |q|² = q₀² + |**v**|². Two facts follow and do the work of the lesson: (p ⊗ q)* = q* ⊗ p* (check with (2): the cross product changes sign when the order is swapped), and therefore

**|p ⊗ q|² = (p⊗q)(p⊗q)* = p ⊗ q ⊗ q* ⊗ p* = |p|²|q|²:**   (3)

the norm is multiplicative, so the product of two unit quaternions is a unit quaternion. The inverse is q⁻¹ = q*/|q|²; for a unit quaternion it is just the conjugate.

### 2.2 Rotations

A **unit quaternion** can always be written q = (cos(a/2), sin(a/2) **n**) with **n** a unit vector: this is the attitude of the frame b obtained from n by a turn a about **n** (Euler's theorem, `P01.15`). Embed a vector as the pure quaternion (0, **v**). The map

**(0, v_b) = q* ⊗ (0, v_n) ⊗ q**   (4)

keeps the scalar part zero and preserves length by (3). Expanding it with (2) (Problem 2) gives v_b = C(q) v_n with

**C_b/n(q) = (q₀² − |v|²) E + 2 v vᵀ − 2 q₀ [v×].**   (5)

Substitute q₀ = cos(a/2), **v** = sin(a/2)**n**: q₀² − |v|² = cos a, 2vvᵀ = (1 − cos a)**nn**ᵀ, 2q₀[v×] = sin a [**n**×]. That is exactly the passive Rodrigues matrix of `P01.15`. Equation (5) is quadratic in q, so **q and −q give the same C**: every attitude has two quaternions, a **double cover**. The angle a and the turn a − 2π about the same axis are the same attitude with opposite quaternions. Equation (5) also shows what a drifting norm does: if |q| = 1 + ε, C(q) = (1 + ε)² × a rotation, an orthogonality error ‖CᵀC − E‖ of about 4√3 ε ≈ 6.9 ε.

**Composition.** Turn n into b by q_b/n, then b into c by q_c/b (a turn defined about the axes of b). From (4) twice, v_c = q_c/b* ⊗ q_b/n* ⊗ v_n ⊗ q_b/n ⊗ q_c/b, so

**q_c/n = q_b/n ⊗ q_c/b,**   (6)

the newer turn on the right, mirroring C_c/n = C_c/b C_b/n with the newer matrix on the left. The 3-2-1 attitude is therefore q = q₃(ψ) ⊗ q₂(θ) ⊗ q₁(φ), with qₖ(a) the turn a about axis k.

**Example: the lesson's attitude.** (ψ, θ, φ) = (30°, 10°, −5°) gives q = (0.960350, −0.064509, 0.072859, 0.261261): one turn of 32.3767° about (−0.2314, 0.2613, 0.9371), the same axis and angle `P01.15` found from the DCM. Gravity (0, 0, 9.80665) m/s² in NED becomes (−1.7029, −0.8417, 9.6209) m/s² in body axes by (4), the same as by the DCM.

### 2.3 Back from a DCM, and to and from Euler angles

From (5): tr C = 3q₀² − |v|² = 4q₀² − 1, C₁₁ = q₀² + q₁² − q₂² − q₃², and so on; the off-diagonal differences give C₂₃ − C₃₂ = 4q₀q₁, C₃₁ − C₁₃ = 4q₀q₂, C₁₂ − C₂₁ = 4q₀q₃, and the sums give C₁₂ + C₂₁ = 4q₁q₂, C₁₃ + C₃₁ = 4q₁q₃, C₂₃ + C₃₂ = 4q₂q₃. The naive inversion q₀ = ½√(1 + tr C), qᵢ = (difference)/(4q₀) fails near a = 180°, where q₀ → 0. **Shepperd's method** computes all four of 4q₀² = 1 + tr C, 4q₁² = 1 + 2C₁₁ − tr C, 4q₂² = 1 + 2C₂₂ − tr C, 4q₃² = 1 + 2C₃₃ − tr C, takes the square root of the largest (at least 1/4 of the sum, 1), and gets the other three components by dividing a sum or difference by it. No small number is ever square-rooted or divided by. The sign is then chosen with q₀ ≥ 0.

The 3-2-1 angles follow from the DCM entries written in q (TM-74839 2.3 tabulates all twelve sequences, in its transposed convention):

ψ = atan2(2(q₀q₃ + q₁q₂), 1 − 2(q₂² + q₃²)),  θ = arcsin(2(q₀q₂ − q₁q₃)),  φ = atan2(2(q₀q₁ + q₂q₃), 1 − 2(q₁² + q₂²)).   (7)

Gimbal lock lives here, in the angles, not in the quaternion: at θ = ±90° the code puts the rotation into ψ, as `P06.07` did.

### 2.4 Kinematics

Let q(t) be the attitude and ω the body rate. Over a short time Δt the body turns by the small angle |ω|Δt about ω/|ω|, a turn defined in body axes, so by (6) q(t + Δt) = q(t) ⊗ (1, ωΔt/2) to first order. Hence

**q̇ = ½ q ⊗ (0, ω).**   (8)

It is linear in q, has no singularity, and uses four states. From (2), q · q̇ = ½ (−q₀ ω·v + q₀ v·ω + v·(v × ω)) = 0: the exact solution keeps |q| = 1. Differentiating (5) along (8) gives Ċ = −[ω×]C, the `P06.07` equation (a test checks it by finite differences). For a rate held constant over a step the exact solution is

**q(t + Δt) = q(t) ⊗ (cos(|ω|Δt/2), sin(|ω|Δt/2) ω/|ω|),**   (9)

the exponential step, free of truncation error and of norm drift (to rounding).

**RK4 and the norm.** For a constant ω, an RK4 step multiplies q by the series of exp(½(0, ω)Δt) cut after five terms. That is RK4's amplification factor for an oscillator at frequency |ω|/2, whose squared magnitude is 1 − x⁶/72 + O(x⁸) with x = |ω|Δt/2 (`P06.07` Problem 2e has the same calculation at x = |ω|Δt): the norm falls by about x⁶/144 per step, so over a fixed time it drifts at order Δt⁵, one order above the attitude error. Because the quaternion turns at half the body rate, its truncation error at a given Δt is about 2⁴ = 16 times smaller than the DCM's: at Δt = 0.0125 s on the coning case, 2.92 × 10⁻⁸ rad against `P06.07`'s 4.52 × 10⁻⁷ rad. Renormalising, q ← q/|q|, is one division; like the DCM's SVD it restores the constraint but not the attitude.

### 2.5 Interpolation: SLERP

Unit quaternions lie on the unit sphere in four dimensions, and a constant-rate turn about a fixed axis is a great circle on it. Between q₀ and q₁ with cos Ω = q₀ · q₁, the great-circle path at constant speed is

**slerp(q₀, q₁, u) = [sin((1 − u)Ω) q₀ + sin(uΩ) q₁] / sin Ω,  0 ≤ u ≤ 1.**   (10)

Check: at u = 0 and 1 it gives the endpoints; its norm is 1 for all u (expand |…|² with q₀·q₁ = cos Ω); and the angle between successive points is proportional to the step in u, so the attitude turns at a constant rate about a fixed axis, the "straight line" between two attitudes. Because q₁ and −q₁ are the same attitude but sit on opposite sides of the sphere, (10) must first replace q₁ by −q₁ if q₀ · q₁ < 0, or it takes the long way round: a 170° turn interpolated halfway gives 85° along the short arc but 95° the other way without the check. **NLERP**, the normalised straight line ((1 − u)q₀ + uq₁)/|…|, follows the same arc but not at constant speed; for small Ω the difference is negligible.

Interpolating samples of a smoothly turning body with (10) treats each interval as a constant-rate turn; the error is that of linear interpolation, of order Δt² in the sample interval.

### 2.6 Where the model stops being valid

- **Exact samples and rates.** Real star trackers have noise and occasional outages; real gyros have bias. A filter, not interpolation, combines them.
- **Sign continuity.** A tracker, a filter and a ground tool may each choose the sign of q differently (q₀ ≥ 0, or continuity with the last sample). Any arithmetic on quaternions from different sources (averages, differences, interpolation) must align signs first.
- **Convention.** Scalar-first or scalar-last, active or passive, Hamilton or the "JPL" product order: published formulas differ in all three. TM-74839 itself uses the transposed matrix. Every match in this lesson is checked numerically, not by eye.

## 3. Build the computational model

Open `starter/main.py`. The `P06.07` helpers (`rot`, `skew`, `dcm_from_euler_321`, `dcm_from_axis_angle`, `misalignment_angle`, `integrate_dcm`), the coning closed form (`coning_rate`, `coning_dcm`) and `main` are given.

### Step 1: the algebra

`quat_multiply` is (2) written out in components; `quat_conjugate`, `quat_norm`, `quat_normalize`, `quat_inverse` follow 2.1. `quat_from_axis_angle` and `axis_angle_from_quat` are 2.2, with the sign chosen so that q₀ ≥ 0 and the angle from 2 atan2(|v|, q₀).

### Step 2: rotations

`transform(q, v)` is (4); `dcm_from_quat` is (5) without normalising (so that a drifting norm shows); `quat_from_dcm` is Shepperd's method; `quat_from_euler_321` composes by (6); `euler_321_from_quat` is (7) with the gimbal-lock branch.

### Step 3: kinematics and integration

```python
def quat_rate(q, w):
    return 0.5 * quat_multiply(q, np.array([0.0, w[0], w[1], w[2]]))
```

`integrate_quat` passes it to `solve_ivp` (DOP853, rtol 10⁻¹⁰) and does not renormalise, so the drift is visible; `integrate_quat_rk4` has an optional renormalisation; `quat_exp_step` is (9). `quat_error_angle(qa, qb)` is the angle of qb* ⊗ qa, 2 atan2(|vector|, |scalar|), which ignores the sign of either quaternion.

### Step 4: interpolation

`slerp` is (10) with the shortest-arc check; `nlerp` the normalised linear form (with an option to skip the check, to show what it prevents); `interpolate_samples` applies either sample interval by sample interval.

Running `main()` prints the lesson's attitude in every form, the coning comparison, the RK4 study, the exact step, the interpolation study and the double-cover example, and writes `quaternions.png`.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-08-quaternions`. The tests:

1. **Algebra**: Hamilton's rules, non-commutativity, associativity, multiplicative norm, identity, conjugate of a product, inverse.
2. **Axis–angle**: round trip for 20 turns, for q and −q.
3. **DCM**: equal to Rodrigues (the `P06.07` route) and to SciPy's `Rotation` (scalar-last, active, transposed) to 10⁻¹⁴; a rotation; the same for −q; a scaled q scales C by |q|².
4. **Shepperd**: round trip for 30 random attitudes and for turns of π − 10⁻³ and π about each axis.
5. **Transform** equal to C v, length preserved, the 30° yaw sign check of `P01.15`.
6. **Composition**: C(q_b/n ⊗ q_c/b) = C(q_c/b) C(q_b/n).
7. **Euler**: the 3-2-1 quaternion's DCM equal to the 3-2-1 DCM; angles back; gimbal lock rebuilds the DCM.
8. **Error angle**: symmetric, sign-blind, equal to the DCM misalignment.
9. **Kinematics**: q · q̇ = 0; a finite difference of C(q) along q̇ equals −[ω×]C.
10. **Integration**: the quaternion route against the DCM route and against the closed form to 10⁻⁹; the exact step after 500 steps to 10⁻¹³; RK4 order 3.7–4.3 for the attitude and 4.5–5.5 for the norm; renormalisation holds |q| = 1 and leaves the attitude error unchanged.
11. **SLERP**: endpoints, unit norm, equal steps, one great arc; NLERP on the same arc but not at constant rate; the shortest arc for q and −q, the long arc without the check.
12. **Star tracker**: SLERP interpolation of coning samples converges at order 1.8–2.3.
13. `main` writes the figure.

### Compare

The independent route is the `P06.07` DCM. Integrating the coning spacecraft's body rates for 30 s by (8) and by Ċ = −[ω×]C, both with `solve_ivp` at rtol 10⁻¹⁰, the two attitudes differ by at most **1.80 × 10⁻¹⁰ rad**; the quaternion differs from the closed form by **7.21 × 10⁻¹⁰ rad**, and its norm drifts by **2.85 × 10⁻¹⁰**. SciPy's `Rotation` (used only in the tests) agrees with (5) to 10⁻¹⁴ once its scalar-last order and active convention are accounted for. The exact step (9) carries a 30°/s rotation for 20 s in 2000 steps to **6.1 × 10⁻¹⁵ rad** with a norm error of **8.7 × 10⁻¹⁴**. The usual mistakes show up at once: writing q̇ = ½(0, ω) ⊗ q (the reference-frame form) fails the finite-difference test; a scalar-last array passed where scalar-first is expected fails the SciPy comparison; and forgetting the factor ½ in (8) doubles every rotation.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements. Write each prediction down before running it.

1. **RK4 on the quaternion.** Integrate the coning spacecraft for 30 s with RK4 at Δt = 0.1, 0.05, 0.025 and 0.0125 s, without and with renormalisation. Predict the order of the attitude error and of the norm drift, and the ratio to `P06.07`'s DCM errors. *Finding:* attitude error 1.191 × 10⁻⁴, 7.472 × 10⁻⁶, 4.674 × 10⁻⁷, 2.922 × 10⁻⁸ rad (orders 3.99, 4.00, 4.00), about 16 times below the DCM's at each step; norm drift 5.98 × 10⁻⁶, 1.87 × 10⁻⁷, 5.85 × 10⁻⁹, 1.83 × 10⁻¹⁰ (order 5.0), which makes the raw C(q) non-orthogonal by 4.1 × 10⁻⁵ … 1.3 × 10⁻⁹; with renormalisation every step the attitude error is identical to four figures.
2. **Star-tracker interpolation.** Sample the coning attitude at intervals of 1, 0.5, 0.25 and 0.125 s, interpolate at 100 Hz with SLERP and with NLERP, and record the largest error over 30 s. Predict the order and the interval that meets 10⁻³ rad. *Finding:*

   | interval (s) | SLERP (rad) | NLERP (rad) |
   |---|---|---|
   | 1.0 | 3.953 × 10⁻² | 6.426 × 10⁻² |
   | 0.5 | 8.847 × 10⁻³ | 9.560 × 10⁻³ |
   | 0.25 | 2.152 × 10⁻³ | 2.154 × 10⁻³ |
   | 0.125 | 5.347 × 10⁻⁴ | 5.348 × 10⁻⁴ |

   SLERP converges at order 2.16, 2.04, 2.01. Meeting 10⁻³ rad needs samples about every 0.17 s (between the last two rows; interpolate on log–log axes); the coning body turns at |ω| = 2.38 rad/s, so that is about 0.40 rad between samples. NLERP is as good as SLERP once the turn per interval is small (0.6 rad or less here), because both follow the same arc and differ only in speed along it.
3. **The double cover.** Take q_a = identity and q_b = a 170° yaw. Interpolate halfway with SLERP from q_a to q_b and to −q_b, and with NLERP to −q_b without the shortest-arc check. *Finding:* 85.0° both ways with SLERP; 95.0° the other way round without the check, a 190° turn instead of a 170° one. Then average ten noisy samples whose signs are random (Problem 3): the naive mean is almost zero.

Raw data format: CSV tables `dt_s, err_raw_rad, norm_drift, orth_raw, err_renorm_rad` (1); `interval_s, slerp_max_rad, nlerp_max_rad` (2); `case, angle_deg` (3); plus `quaternions.png`.

## 5. Validate and improve

| quantity | independent route | quaternion code | difference | cause |
|---|---|---|---|---|
| C(q) for random turns | Rodrigues; SciPy | `dcm_from_quat` | < 10⁻¹⁴ | rounding |
| lesson attitude axis–angle | `P01.15`: 32.38° about (−0.231, 0.261, 0.937) | 32.3767° about (−0.2314, 0.2613, 0.9371) | — | — |
| coning, 30 s, solve_ivp | DCM route | `integrate_quat` | 1.80 × 10⁻¹⁰ rad | solver tolerance |
| coning, 30 s, solve_ivp | closed form | `integrate_quat` | 7.21 × 10⁻¹⁰ rad | solver tolerance |
| coning, RK4 Δt = 0.0125 s | closed form; DCM RK4 4.52 × 10⁻⁷ rad | 2.92 × 10⁻⁸ rad | factor 15.5 | quaternion turns at ω/2 |
| constant rate, exact step | closed form | `quat_exp_step` | 6.1 × 10⁻¹⁵ rad | rounding |
| SLERP order | 2 (linear interpolation) | 2.16, 2.04, 2.01 | — | asymptotic range |

The quaternion reproduces the DCM route and the closed forms to the solver tolerance. The engineering answer: **propagate the quaternion between star-tracker reports with the exact step (9) for each gyro sample (or RK4 if the rate varies within a sample), renormalise every step (one division costs nothing and keeps C(q) orthogonal to rounding), and interpolate between reports with SLERP: for this 2 rad/s coning motion, reports every 0.17 s keep the interpolated attitude within 10⁻³ rad; a filter that propagates with the gyros, rather than interpolating, removes the requirement altogether.**

What to improve: real star-tracker noise (10⁻⁵–10⁻⁴ rad) is far below the interpolation errors at these intervals, so the sample interval, not the tracker, limits the interpolated attitude; propagating with gyro rates (the attitude filter of `P18.10`) is the next step.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (5), (6), (8) and the properties of (10); `quaternions.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the spin rate doubles?** The turn per star-tracker interval doubles and the SLERP error grows about four times; the interval must halve to keep 10⁻³ rad. The RK4 error at fixed Δt grows about 2⁴ = 16 times.
- **What if the rotation per step goes to zero?** `slerp` divides by sin Ω; the code falls back to the normalised linear form when q₀ · q₁ > 1 − 10⁻¹², and `quat_exp_step` uses the first-order form when |ω|Δt < 10⁻¹². Which tests exercise those branches?
- **Which assumption fails first in the real world?** Exact samples: a real tracker's noise of 10⁻⁴ rad would make it pointless to interpolate below that level; and sign conventions between tools.
- **What measurement could falsify the code?** Rotate a phone through 180° about one axis on a table while logging its gyro; integrate with `quat_exp_step` and compare the final axis and angle with the known turn. A sign error in (8) gives the right angle about the wrong axis direction after a two-axis sequence.
- **Does the result survive a change of method?** solve_ivp, RK4 and the exact step agree to their tolerances; the DCM route agrees to 10⁻¹⁰. Changing the sign convention (q₀ ≥ 0 or continuity) changes the numbers stored but not the attitude.
- **Which decision changes if the uncertainty doubles?** If the attitude budget were 2 × 10⁻³ rad instead of 10⁻³, the sample interval could grow by about √2, to roughly 0.24 s; the renormalisation policy would not change.

## Check yourself

1. Compute i ⊗ j and j ⊗ i. What does the difference tell you about composing rotations?
2. Write the unit quaternion of a 90° turn about axis 3. What is its conjugate, and what attitude does it describe?
3. Why do q and −q describe the same attitude? Name one computation where this matters.
4. In which order are the quaternions multiplied for "first turn by q_b/n, then by q_c/b about the new axes"?
5. Show in one line that q̇ = ½ q ⊗ (0, ω) keeps |q| constant.
6. Why is RK4's error on the quaternion about 16 times smaller than on the DCM at the same step?
7. What does renormalising a quaternion fix, and what does it not fix?
8. When does NLERP differ noticeably from SLERP, and why does SLERP need a shortest-arc check?
