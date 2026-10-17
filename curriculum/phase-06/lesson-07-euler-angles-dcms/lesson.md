# P06.07  3D Rotations: Euler Angles & DCMs

> Three angles are a good way to talk about an attitude and a bad way to integrate one. Carry the matrix; read the angles off it when a human needs them.

| | |
|---|---|
| **Learning outcomes** | 1. Build the direction cosine matrix (DCM) of any of the twelve Euler-angle sequences and extract the angles of any sequence from a DCM, including the gimbal-lock case. 2. Derive the matrix B that turns Euler-angle rates into body rates, show that \|det B\| = \|cos a₂\| (or \|sin a₂\| for a symmetric sequence), and explain gimbal lock as the singularity of B. 3. Derive the DCM kinematic equation Ċ = −[ω×]C from the orthogonality of C and integrate it with `solve_ivp` and with fixed-step RK4, with SVD re-orthonormalisation. 4. Verify the integration against two closed-form attitude histories (a fixed-axis rotation and steady coning) and against SciPy's rotation library, and measure the convergence order. 5. Decide, for a given manoeuvre and computer rate, whether Euler angles can be integrated or whether the DCM must be. |
| **Prerequisites** | `P01.15` (elementary rotations, the passive 3-2-1 DCM, Rodrigues, the small-angle strapdown step), `P01.07` (cross products), `P01.16` (eigenvalues), `P03.05` (Runge–Kutta methods), `P03.06` (adaptive integrators), `P06.04` (rotating frames, the transport theorem), `P06.05` (angular velocity in the plane). Skill check: you can multiply three 3 × 3 matrices and write [v×] for a given v. |
| **Engineering question** | An aerobatic aeroplane flies a stall turn: it pulls up to within δ of the vertical, yaws through 180°, and comes down. Its attitude computer integrates gyro rates at 100 Hz. For which δ can it integrate yaw, pitch and roll directly, and for which must it integrate the DCM? And how often must a spinning spacecraft's DCM be re-orthonormalised to stay a rotation? |
| **Assumptions and limits** | Rigid body, so one matrix describes the whole vehicle; the reference frame (north–east–down, or inertial) is treated as non-rotating; body rates are known exactly (gyro errors are a later topic); every matrix is passive (`P01.15`): v_b = C v_n. The Euler-angle rate equations are singular at a₂ = ±90° (asymmetric sequences) or a₂ = 0, 180° (symmetric); the DCM has no singularity anywhere. |
| **Analytical work** | The DCM of a general sequence as a product of three elementary rotations; the extraction formulas for asymmetric and symmetric sequences with a parity sign; the gimbal-lock case; the rate matrix B and its determinant; Ċ = −[ω×]C from CCᵀ = E; d(CᵀC)/dt = 0; the closed forms for fixed-axis rotation and for 3-1-3 coning; the RK4 truncation order; why re-orthonormalisation does not remove attitude error. |
| **Computational deliverable** | `starter/main.py`: `dcm_from_euler`, `euler_from_dcm`, `euler_rate_matrix`, `euler_angle_rates`, `dcm_rate`, `orthogonality_error`, `orthonormalize_svd`, `misalignment_angle`, `integrate_dcm`, `integrate_dcm_rk4`, `integrate_euler_angles`, `integrate_euler_angles_rk4`, `single_axis_dcm`, `coning_dcm`, `coning_body_rate`, `hammerhead` (`rot`, `skew`, `dcm_from_axis_angle`, `is_symmetric`, `piecewise_rate`, `integrate_segments_rk4` and `main` are given). Tests prove: every DCM is a rotation; all twelve sequences agree with SciPy to 10⁻¹⁴; the 3-2-1 matrix entry by entry; angles round-trip and rebuild any random attitude; gimbal lock still rebuilds the matrix; B against the textbook 3-2-1 equations, against a finite difference of the DCM, and its determinant; the integrated DCM against both closed forms to 10⁻⁸; RK4 at order 4.0 ± 0.3; SVD gives the nearest rotation and fixes orthogonality without fixing the attitude; Euler-angle integration fails near vertical while the DCM does not. |
| **Practical activity** | **computational experiment**: a stall turn flown ever closer to the vertical, integrated by Euler angles and by DCM, adaptively and at a fixed 100 Hz; a step-size study of RK4 on a coning spacecraft with and without re-orthonormalisation; the gimbal-lock map of all twelve sequences for a rocket on the pad. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `attitude_kinematics.png`, the derivations of B, its determinant and Ċ = −[ω×]C, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: B and the DCM kinematic equation are derived (not quoted); all public tests pass; the integrated DCM is compared with a closed form with numbers; the note states the δ below which Euler-angle integration at 100 Hz fails and why. |
| **References** | `@ntrs-nasa-tm-74839-euler-angles-quaternions` 2.1 Euler Angle Transformation Matrices; Appendix A Relationships for the Three-Axis Euler Angle Rotation Sequences (required-free: all twelve sequences and their extraction formulas, a check table for your code; **its matrices map the rotated frame into the original frame, the transpose of this lesson's C**, so transpose before comparing). `@mit-ocw-16-07-dynamics` L3 "Vectors, matrices and coordinate transformations"; L25 "3D rigid body kinematics" (Rotation About a Fixed Point; Angular Velocity About a Fixed Point) (required-free). `@libretexts-cline-variational-principles` 13.13 Euler Angles; 13.14 Angular Velocity (required-free: Euler-angle rates from ω). `@mit-ocw-18-06-linear-algebra` Lecture 17: Orthogonal matrices and Gram-Schmidt (required-free: re-orthonormalisation). `@mit-ocw-16-61-aerospace-dynamics` lecture 11 "Kinematics of Rigid Bodies" (supplementary-free). Context: `@nasa-basics-of-space-flight` Chapter 11 "Onboard Systems", Attitude and Articulation Control Subsystems (AACS) (supplementary-free, qualitative). |
| **Extensions** | ★ Replace RK4 by the exact exponential step C ← exp(−[ω×]Δt) C (Rodrigues) for a piecewise-constant rate and show the error drops to rounding; then make the rate vary within the step and measure the coning error that remains. ★ Find, for each of the twelve sequences, the attitudes a vertical rocket passes through during a gravity turn, and pick the sequence that stays furthest from its singularity. ★ Implement P01.15's first-order re-orthonormalisation and compare its cost and residual with the SVD. |
| **Time** | unmeasured |

## 1. The engineering problem

Every guidance, navigation and control computer carries an attitude: which way the vehicle points. Gyros measure the body's angular velocity; the computer integrates it into an orientation; the autopilot, the navigator and the display all read that orientation. `P01.15` built the tool for one fixed attitude, the 3-2-1 direction cosine matrix (DCM), and integrated a constant rotation with small-angle steps. This lesson finishes the job. It generalises to all twelve Euler-angle sequences (spacecraft and orbit work use 3-1-3; some launch vehicles use 2-1-3 or 1-2-3), it writes the equations that propagate an attitude under rates that change in time, and it shows where three angles stop working.

The test case is a stall turn (a hammerhead). An aerobatic aeroplane pulls up at 0.5 rad/s until its nose is δ short of vertical, kicks the rudder and yaws through 180° at 0.5 rad/s, and pulls through to come down the way it went up. Nothing about the motion is violent: the body rates never exceed 0.5 rad/s. Yet as δ shrinks, the yaw-pitch-roll angles that describe it change faster and faster: at δ = 0.01°, yaw must change at 2292 rad/s at the top. A computer that integrates the three angles at 100 Hz loses 1.1° of attitude in that one manoeuvre; one that integrates the DCM at the same rate loses 10⁻¹¹ rad. A second case, a spacecraft coning about its angular momentum (spin 2 rad/s, precession 0.4 rad/s, half-cone 20°), tests the integrator: the DCM's nine numbers are not forced to stay a rotation, and they drift unless something puts them back.

Success criteria, written first: every DCM orthonormal to 10⁻¹⁴; all twelve sequences equal to an independent library; the integrated DCM equal to the closed-form attitude to 10⁻⁸ for both test motions; RK4 converging at fourth order; re-orthonormalisation holding CᵀC − E at rounding level.

## 2. Understand and derive

### 2.1 Any sequence of three turns

`P01.15` defined the passive elementary rotations R₁(a), R₂(a), R₃(a): the components of a fixed vector in a frame turned by a about axis 1, 2 or 3. Any attitude is reached by three successive turns about body axes, the second and third about axes that the previous turns have moved. Write a sequence as three axis digits in the order of the turns, s₁s₂s₃, with angles a₁, a₂, a₃. Each turn is passive from the previous intermediate frame, so the components chain right to left:

**C = C_b/n = R_{s₃}(a₃) R_{s₂}(a₂) R_{s₁}(a₁).**   (1)

Two neighbouring axes may not be the same (a turn about 3 followed by another about 3 is one turn). That leaves 3 × 2 × 2 = **twelve sequences**: six **asymmetric** ones that use all three axes (1-2-3, 1-3-2, 2-1-3, 2-3-1, 3-1-2, 3-2-1) and six **symmetric** ones whose first and last axes coincide (1-2-1, 1-3-1, 2-1-2, 2-3-2, 3-1-3, 3-2-3). Aircraft use 3-2-1 (yaw ψ, pitch θ, roll φ). Orbit work uses 3-1-3: right ascension of the ascending node Ω, inclination i, argument of periapsis ω (`P01.15` extension). A spinning spacecraft is often described by 3-1-3 as precession, nutation and spin (16.07 L29). Henderson's NASA report lists all twelve matrices with their extraction formulas (TM-74839 Appendix A); it uses the transposed convention, mapping rotated-frame components into the original frame, so compare its tables with Cᵀ.

Every product of rotations is a rotation (CᵀC = E, det C = +1, `P01.15` and 18.06 Lecture 17), so (1) gives a valid DCM for any angles.

**Example: the lesson's aircraft attitude** (ψ, θ, φ) = (30°, 10°, −5°) has the 3-2-1 DCM

```
C = [[ 0.85287  0.49240 -0.17365]
     [-0.51120  0.85516 -0.08583]
     [ 0.10623  0.16197  0.98106]]
```

The same physical attitude is (146.7402°, 11.1690°, −116.3025°) in 3-1-3, (−9.3749°, 6.0982°, 30.9382°) in 1-2-3 and (6.1802°, −9.3214°, 29.9334°) in 2-1-3. One attitude, twelve sets of angles: the angles are coordinates, not the attitude.

### 2.2 Getting the angles back for any sequence

The extraction is easiest on the transpose. Since Rₙ(a)ᵀ is the active rotation Aₙ(a), Cᵀ = A_{s₁}(a₁) A_{s₂}(a₂) A_{s₃}(a₃). Number the axes from 0 (i = s₁ − 1, j = s₂ − 1) and let k be s₃ − 1 for an asymmetric sequence, or the unused axis for a symmetric one. Let ε = +1 if (i, j, k) is a cyclic order of (0, 1, 2) and −1 otherwise. Multiplying out the general product (do it once for 3-2-1 and once for 3-1-3; the tests check all twelve) gives, with M = Cᵀ:

| | middle angle | first angle | last angle |
|---|---|---|---|
| asymmetric | sin a₂ = ε M_ik | a₁ = atan2(−ε M_jk, M_kk) | a₃ = atan2(−ε M_ij, M_ii) |
| symmetric | cos a₂ = M_ii | a₁ = atan2(M_ji, −ε M_ki) | a₃ = atan2(M_ij, ε M_ik) |

For 3-2-1 (i = 2, j = 1, k = 0, ε = −1) the first row reduces to `P01.15`'s θ = −arcsin C₁₃, ψ = atan2(C₁₂, C₁₁), φ = atan2(C₂₃, C₃₃). The middle angle lies in [−90°, 90°] for an asymmetric sequence and [0°, 180°] for a symmetric one.

**Gimbal lock in the extraction.** When cos a₂ = 0 (asymmetric) or sin a₂ = 0 (symmetric), all four entries that fix a₁ and a₃ vanish, and only their sum or difference is defined: the first and last axes have become parallel. The code then sets a₃ = 0 and recovers a₁ from A_{s₁}(a₁) = M A_{s₂}(a₂)ᵀ. The rebuilt matrix is exact; only the split between a₁ and a₃ is a convention.

### 2.3 Euler-angle rates and the singularity

Each angle rate is a rotation rate about an axis of an intermediate frame: ȧ₁ about the first axis (seen before the second and third turns), ȧ₂ about the second (seen before the third), ȧ₃ about the body axis. Carrying each into body components with the turns that follow it (16.07 L25 "Angular Velocity About a Fixed Point"; Cline 13.14):

**ω_b = R_{s₃}(a₃) R_{s₂}(a₂) e_{s₁} ȧ₁ + R_{s₃}(a₃) e_{s₂} ȧ₂ + e_{s₃} ȧ₃ = B(a) ȧ.**   (2)

For 3-2-1, with (ȧ₁, ȧ₂, ȧ₃) = (ψ̇, θ̇, φ̇) and ω_b = (p, q, r):

p = φ̇ − ψ̇ sin θ,  q = θ̇ cos φ + ψ̇ sin φ cos θ,  r = −θ̇ sin φ + ψ̇ cos φ cos θ.   (3)

The gyros measure (p, q, r); to propagate the angles the computer needs the inverse,

ψ̇ = (q sin φ + r cos φ)/cos θ,  θ̇ = q cos φ − r sin φ,  φ̇ = p + (q sin φ + r cos φ) tan θ.   (4)

**The determinant.** The columns of B are three unit vectors: the first axis carried by two turns, the second carried by one, and the third. Its determinant is their triple product. The second and third columns are perpendicular (adjacent axes), so det B depends only on how the first column leans out of their plane. For an asymmetric sequence the first axis is perpendicular to the second, and after the turn a₂ it makes the angle 90° − a₂ with the third: **|det B| = |cos a₂|**. For a symmetric sequence the first and third axes coincide before the middle turn, so **|det B| = |sin a₂|**. For 3-2-1, det B = −cos θ: at the lesson's attitude, −0.984808. When det B = 0 the three rotation axes are coplanar and (4) divides by zero: **gimbal lock** is the singularity of B. It is not in the motion; it is in the angles. Near it the angle rates are huge even for small body rates: at θ = 89.9° with φ = 30°, a yaw rate r = 0.1 rad/s needs ψ̇ = 0.1 cos 30°/cos 89.9° = 49.62 rad/s.

### 2.4 The DCM kinematic equation

Differentiate CCᵀ = E: Ċ Cᵀ + C Ċᵀ = 0, so Ω = −ĊCᵀ is skew-symmetric, Ω = [ω×] for some vector ω. Then

**Ċ = −[ω×] C,**   (5)

and the transport theorem of `P06.04` identifies ω: for a vector fixed in n, its body components v_b = C v_n change at v̇_b = Ċv_n = −ω × v_b, which is how a fixed vector looks from a frame turning at ω. So ω is the angular velocity of b relative to n in body components, exactly what a strapdown gyro triad measures. Equation (5) is linear in C, has no singularity anywhere, and costs nine multiply-adds per entry; its price is nine numbers for three degrees of freedom. The exact solution keeps CᵀC = E (differentiate: ĊᵀC + CᵀĊ = Cᵀ([ω×] − [ω×])C = 0, a test), but a numerical integration does not.

**Two closed forms to check it against.**

- **Fixed axis.** If ω = ω n with n constant in body axes, n is also fixed in n's frame (Cᵀn is constant), and C(t) = C(n, ωt) C₀ with C(n, a) the passive Rodrigues matrix of `P01.15`.
- **Steady coning.** A 3-1-3 motion with constant precession rate ψ̇, constant half-cone angle β and constant spin φ̇ has C(t) = R₃(φ̇t) R₁(β) R₃(ψ̇t), and from (2) the body rates

ω_b = (ψ̇ sin β sin φ, ψ̇ sin β cos φ, φ̇ + ψ̇ cos β),  φ = φ̇t.   (6)

The transverse body rate rotates in the body at φ̇ while the body's axis 3 sweeps a cone of half-angle β around the reference axis 3. This is the torque-free motion of a symmetric spacecraft (`P06.10`) and the motion that makes strapdown integration hard ("coning"): the rate vector never stays put in the body.

### 2.5 Integration error, drift and repair

Classical RK4 on (5) has a global error of order Δt⁴. Its update is a polynomial in [ω×]Δt, not the exponential, so CᵀC − E grows each step. Re-orthonormalising replaces C by the nearest rotation. In the Frobenius norm that is the polar factor: write C = U S Vᵀ (singular value decomposition) and take

**C ← U diag(1, 1, det(UVᵀ)) Vᵀ,**   (7)

the exact version of `P01.15`'s first-order step C(3E − CᵀC)/2, and of Gram–Schmidt (18.06 Lecture 17), with the advantage that it treats the three rows alike. It removes the drift of CᵀC, but not the attitude error: the nearest rotation to a slightly wrong matrix is a slightly wrong rotation. The tests check both halves of that sentence.

### 2.6 Where the model stops being valid

- **Rigid body.** One matrix per vehicle. A large wing flexes by degrees.
- **Exact rates.** Real gyros have bias, scale-factor error and noise; their effect usually exceeds the integration error measured here by orders of magnitude.
- **Non-rotating reference.** North–east–down turns with the Earth at 7.29 × 10⁻⁵ rad/s; `P06.04` adds it.
- **Piecewise-constant rates.** The stall turn's rates switch instantly; a real aeroplane's rates change smoothly. A fixed-step integrator that steps across a switch loses an order of accuracy, which is why `integrate_segments_rk4` restarts at each switch.

## 3. Build the computational model

Open `starter/main.py`. `rot`, `skew`, `dcm_from_axis_angle` (all from `P01.15`), `is_symmetric`, `piecewise_rate`, `integrate_segments_rk4` and `main` are given.

### Step 1: DCMs and angles

`dcm_from_euler(angles, seq)` is (1) with `rot`. `euler_from_dcm(C, seq)` is the table of 2.2, with the parity sign computed from (i, j, k), `acos`/`asin` arguments clipped to [−1, 1] against rounding, and the gimbal-lock branch.

### Step 2: rates

`euler_rate_matrix(angles, seq)` builds the three columns of (2). `euler_angle_rates` solves B ȧ = ω with `np.linalg.solve` and raises `ValueError` when |det B| is below a tolerance.

### Step 3: the kinematic equation and its tools

```python
def dcm_rate(C, w):
    return -skew(w) @ C

def orthonormalize_svd(C):
    U, _, Vt = np.linalg.svd(C)
    return U @ np.diag([1.0, 1.0, np.sign(np.linalg.det(U @ Vt))]) @ Vt
```

`misalignment_angle(Ca, Cb)` is the rotation angle of Ca Cbᵀ, computed with atan2 as in `P01.15`. `integrate_dcm` passes the nine entries to `solve_ivp` (DOP853, rtol 10⁻¹⁰); `integrate_dcm_rk4` is a fixed-step RK4 with an optional call to (7) every k steps. `integrate_euler_angles` and `integrate_euler_angles_rk4` integrate (4) the same two ways.

### Step 4: the closed forms and the manoeuvre

`single_axis_dcm` and `coning_dcm`, `coning_body_rate` are 2.4. `hammerhead(q, r, theta_top)` returns the three constant-rate segments (pitch, yaw 180°, pitch) and the exact final attitude, the product of three single-axis rotations.

Running `main()` prints the lesson's DCM in four sequences, the rate matrix and the gimbal-lock rates, both closed-form checks, the RK4 study and the stall-turn study, and writes `attitude_kinematics.png`: the 3-2-1 angles through a stall turn 1° from vertical, the Euler-angle integration's error against the DCM along it, and the RK4 error and drift against step size.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-07-euler-angles-dcms`. The tests:

1. **Rotation matrices**: 60 random DCMs over all twelve sequences have CᵀC = E and det C = 1 to 10⁻¹⁴.
2. **3-2-1 entries**: equal to the written-out matrix of `P01.15`.
3. **Oracle**: every sequence equal to `scipy.spatial.transform.Rotation.from_euler` (intrinsic, upper-case axes) transposed, to 10⁻¹⁴. SciPy is used only in the tests.
4. **Bad sequences** refused.
5. **Round trip** for 120 random angle sets; **extraction** of 20 random attitudes in all twelve sequences rebuilds them to 10⁻¹², with the middle angle in range.
6. **Gimbal lock** in four sequences rebuilds the matrix with a₃ = 0.
7. **Rate matrix**: equal to (3) and its inverse (4); |det B| as in 2.3; `ValueError` at lock; equal to a central finite difference of the DCM to 10⁻⁸; d(CᵀC)/dt = 0.
8. **Closed forms**: `solve_ivp` against the fixed-axis rotation (10⁻⁹) and coning (10⁻⁸); the half-cone angle constant; (6) equal to B(ψ̇, 0, φ̇).
9. **RK4**: observed order 3.7–4.3; drift above 10⁻⁴ at Δt = 0.1 s without repair, below 10⁻¹³ with it, and the same attitude error either way.
10. **SVD**: the result is a rotation, a rotation is left alone, and no nearby rotation is closer.
11. **Misalignment**: zero, symmetric, the right angle, π handled.
12. **Euler angles**: integrated 3-1-3 angles equal coning away from lock; in the stall turn at 0.1° from vertical their error at 100 Hz is over 1000 times that at 10°, while the DCM stays below 10⁻⁹.
13. `main` writes the figure.

### Compare

The independent checks are the two closed forms and the SciPy oracle. A 30°/s rotation about (1, 2, 2)/3 for 20 s integrated by `solve_ivp` stays within **4.3 × 10⁻¹¹ rad** of the Rodrigues closed form, with CᵀC − E below **3.9 × 10⁻¹⁰**. The coning spacecraft (ψ̇ = 0.4, φ̇ = 2.0 rad/s, β = 20°) over 30 s stays within **7.6 × 10⁻¹⁰ rad** of its 3-1-3 closed form, with orthogonality error **9.9 × 10⁻¹⁰**. The stall turn's final attitude matches the product of its three rotations to **2–4 × 10⁻¹¹ rad** for every δ. These checks catch the usual mistakes immediately: Ċ = +[ω×]C (the active sign) puts the coning attitude 1.53 rad (88°) away from the closed form after only 1 s; building C in the wrong order (R₁R₂R₃ for 3-2-1) misses the attitude (0.5, 0.3, −0.2) rad by 0.186 rad and fails the SciPy oracle; and a parity sign that is right for one sequence is wrong for others, which is why the tests run all twelve.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model, not measurements. Write each prediction down before running it.

1. **The stall turn.** Fly the hammerhead with q = r = 0.5 rad/s for δ = 10°, 1°, 0.1° and 0.01°. For each, predict the largest yaw rate ψ̇ on the path (from (4), about r/sin δ), then integrate by DCM and by Euler angles, adaptively (rtol 10⁻¹⁰) and with RK4 at 100 Hz, and compare the final attitudes with the closed form. *Finding:*

   | δ | largest \|ψ̇\| (rad/s) | adaptive DCM (rad) | adaptive Euler (rad) | adaptive Euler rhs calls | RK4 100 Hz DCM (rad) | RK4 100 Hz Euler (rad) |
   |---|---|---|---|---|---|---|
   | 10° | 2.9 | 2.9 × 10⁻¹¹ | 1.7 × 10⁻¹¹ | 2654 | 1.6 × 10⁻¹¹ | 4.0 × 10⁻¹⁰ |
   | 1° | 28.1 | 2.1 × 10⁻¹¹ | 8.8 × 10⁻¹² | 3146 | 1.6 × 10⁻¹¹ | 3.6 × 10⁻⁷ |
   | 0.1° | 229.2 | 3.8 × 10⁻¹¹ | 7.0 × 10⁻¹² | 3515 | 1.6 × 10⁻¹¹ | 1.0 × 10⁻³ |
   | 0.01° | 2291.7 | 3.6 × 10⁻¹¹ | 6.2 × 10⁻¹² | 3956 | 1.6 × 10⁻¹¹ | 2.0 × 10⁻² |

   An adaptive integrator survives the singularity by taking ever smaller steps (its cost grows); a flight computer at a fixed rate cannot, and loses 2.0 × 10⁻² rad (1.1°) in one manoeuvre at δ = 0.01°. The DCM's error does not depend on δ at all. At δ = 1° the 3-2-1 yaw read from the DCM changes by 9.72° in 6.2 ms while the aeroplane turns 0.179°, and over the turn yaw goes from 0° to 180°: a correct description that no rate-limited display or filter can follow.
2. **RK4 and re-orthonormalisation.** Integrate the coning spacecraft for 30 s with RK4 at Δt = 0.1, 0.05, 0.025 and 0.0125 s, with and without (7) after every step. Predict the order of the attitude error and what (7) does to it. *Finding:* misalignment 1.812 × 10⁻³, 1.150 × 10⁻⁴, 7.217 × 10⁻⁶, 4.515 × 10⁻⁷ rad (orders 3.98, 3.99, 4.00); orthogonality error 1.04 × 10⁻³, 3.26 × 10⁻⁵, 1.02 × 10⁻⁶, 3.19 × 10⁻⁸ without repair (order 5) and 6–13 × 10⁻¹⁶ with it; the misalignment with repair is the same to four figures. Re-orthonormalisation keeps the matrix a rotation; only a smaller step or a better step reduces the attitude error.
3. **A rocket on the pad.** A launch vehicle standing vertically with its body x axis up has the 3-2-1 angles (0°, 90°, 0°): it is at gimbal lock before it moves. Use `euler_from_dcm` and `euler_rate_matrix` to tabulate |det B| for all twelve sequences at that attitude and predict which are singular. *Finding:* 1-2-3, 2-1-2, 2-3-2 and 3-2-1 are singular (|det B| = 0); the other eight have |det B| = 1. A launch vehicle's guidance either uses one of those eight or, more often, carries a quaternion (`P06.08`).

Raw data format: CSV tables `delta_deg, max_psidot, err_dcm_adapt, err_euler_adapt, nfev_euler, err_dcm_rk4, err_euler_rk4` (1); `dt_s, misalign_raw, orth_raw, misalign_renorm, orth_renorm` (2); `sequence, a1_deg, a2_deg, a3_deg, abs_det_B` (3); plus `attitude_kinematics.png`.

## 5. Validate and improve

| quantity | closed form or oracle | code | difference | cause |
|---|---|---|---|---|
| 12 sequences vs SciPy | SciPy `Rotation` | `dcm_from_euler` | < 10⁻¹⁴ | rounding |
| fixed-axis rotation, 20 s | Rodrigues | `integrate_dcm` | 4.3 × 10⁻¹¹ rad | solver tolerance |
| coning, 30 s | 3-1-3 closed form | `integrate_dcm` | 7.6 × 10⁻¹⁰ rad | solver tolerance |
| coning, RK4 Δt = 0.0125 s | 3-1-3 closed form | `integrate_dcm_rk4` | 4.5 × 10⁻⁷ rad, order 4.00 | RK4 truncation |
| stall turn, δ = 0.01°, 100 Hz | product of three rotations | DCM RK4: 1.6 × 10⁻¹¹ rad; Euler RK4: 2.0 × 10⁻² rad | 9 orders of magnitude | singular B |
| det B, 3-2-1 at (30°, 10°, −5°) | −cos 10° = −0.984808 | −0.984808 | 0 | — |

Every closed form is reproduced to the integrator's tolerance and the RK4 order is the theoretical one. The engineering answer: **at 100 Hz, Euler angles integrate the stall turn to 4 × 10⁻¹⁰ rad when it stays 10° from vertical and to 3.6 × 10⁻⁷ rad at 1°, but lose 10⁻³ rad at 0.1° and 2 × 10⁻² rad at 0.01°; the DCM's error is 1.6 × 10⁻¹¹ rad regardless.** An aerobatic or launch vehicle must integrate the DCM (or a quaternion) and compute angles only for display. For the coning spacecraft, RK4 at 0.0125 s drifts from orthogonality by 3.2 × 10⁻⁸ in 30 s; renormalising every 50 steps holds the largest drift to 6.7 × 10⁻¹⁰ and every 200 steps to 2.7 × 10⁻⁹, at negligible cost.

What to improve: the step itself. For piecewise-constant rates the exact exponential step (★ extension) is free of truncation error; for rates that vary inside a step, coning compensation is the subject of the inertial-navigation lessons. The quaternion (`P06.08`) carries the same information in four numbers with one constraint instead of six.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), det B for both kinds of sequence, (5) and (6); `attitude_kinematics.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the rates double?** The stall turn takes half the time, but ψ̇ at the top doubles too. At δ = 1° and 100 Hz the Euler-angle RK4 error rises from 3.6 × 10⁻⁷ to 5.6 × 10⁻⁶ rad (a factor of 16, as (ωΔt)⁴ with half as many steps predicts), and the DCM error from 1.6 × 10⁻¹¹ to 2.6 × 10⁻¹⁰ rad: both scale the same way, but from very different starting points.
- **What if δ → 0?** The angles stay defined everywhere except at δ = 0 itself, but their rates diverge as 1/sin δ. Exactly at δ = 0 the Euler-angle integration divides by zero; the DCM does not notice.
- **Which assumption fails first in the real world?** Exact gyro rates. A 10⁻³ °/s gyro bias adds 0.03° in 30 s, far more than any integration error here.
- **What measurement could falsify the code?** Rotate a phone through a known sequence of turns on a table (90° about each axis) while logging its gyro, integrate with `integrate_dcm`, and compare the final DCM with the known product of rotations: a wrong sign in (5) or a wrong multiplication order shows up as a wrong final orientation, not as a small error.
- **Does the result survive a change of solver?** The DCM result survives adaptive or fixed-step integration, any δ and any sequence used to read it; the Euler-angle result depends on all of them. That asymmetry is the lesson.
- **Which decision changes if the uncertainty doubles?** If the attitude budget were 10⁻³ rad instead of 2 × 10⁻³ rad, Euler-angle RK4 at 100 Hz would fail already at δ = 0.1° rather than only below it; the DCM decision would not change.

## Check yourself

1. How many distinct Euler-angle sequences are there, and why are sequences like 3-3-1 excluded?
2. In (1), which elementary rotation is applied first, and which matrix appears rightmost?
3. For which values of the middle angle is a 3-2-1 sequence singular? A 3-1-3 sequence?
4. Write det B for an asymmetric sequence and explain in one sentence what it measures geometrically.
5. Derive in two lines why Ċ Cᵀ must be skew-symmetric.
6. What does SVD re-orthonormalisation fix, and what does it not fix?
7. A strapdown integrator at fixed Δt halves its step. By how much should the attitude error fall with RK4, and by how much should the drift from orthogonality fall?
8. The yaw angle of an aeroplane read from its DCM jumps by 180° in a fraction of a second while the aeroplane barely turns. Is the DCM wrong? Explain.
