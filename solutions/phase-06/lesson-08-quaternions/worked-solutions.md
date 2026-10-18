# P06.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` (and the probe scripts that call it) and are given to the figures shown. All "simulated" values are outputs of the model, not measurements.

## 1. Two turns in two orders

**Assumptions.** Passive convention; q_b/n composed with the newer turn on the right (6).

(a) q₃ = (c, 0, 0, c) and q₂ = (c, 0, c, 0) with c = cos 45° = sin 45°, so c² = ½. By (2), yaw then pitch, p = q₃, q = q₂: the scalar part is p₀q₀ − p_v·q_v = c² − 0 = c²; the vector part is p₀q_v + q₀p_v + p_v × q_v = (0, c², 0) + (0, 0, c²) + (0, 0, c) × (0, c, 0) = (0, c², c²) + (−c², 0, 0) = (−c², c², c²). So **q = (0.5, −0.5, 0.5, 0.5)**. Pitch then yaw: q₂ ⊗ q₃ has p_v × q_v = (0, c, 0) × (0, 0, c) = (c², 0, 0), so **q = (0.5, 0.5, 0.5, 0.5)**. The code agrees.

(b) Yaw then pitch: **120° about (−1, 1, 1)/√3**; pitch then yaw: **120° about (1, 1, 1)/√3**. Same angle, different axes: the two orders give different attitudes.

(c) Row 1 of C: yaw then pitch, nose **(0, 0, −1)**, straight up; pitch then yaw, nose **(0, 1, 0)**, pointing east. (Pitching first puts the nose straight up and body axis 3 horizontal, pointing north; the yaw about that axis then swings the nose down to the horizon, towards east.)

(d) North in body axes: yaw then pitch **(0, −1, 0)**; pitch then yaw **(0, 0, 1)**.

**Check.** Each final attitude is a 120° turn about a body diagonal, the attitude that cyclically permutes the axes: consistent with the two 90° turns about perpendicular axes.

**Common mistakes.** Multiplying in matrix order (q₂ ⊗ q₃ for yaw then pitch); taking a column of C instead of a row for the nose.

## 2. The quaternion as a rotation

(a) pq = (p₀ + p_v)(q₀ + q_v) = p₀q₀ + p₀q_v + q₀p_v + p_v q_v. For pure quaternions a = a₁i + a₂j + a₃k, b likewise, the nine products with i² = −1, ij = k = −ji etc. give a b = −(a·b) + a × b. Hence **p ⊗ q = (p₀q₀ − p_v·q_v, p₀q_v + q₀p_v + p_v × q_v)**.

(b) Conjugating (2): (p₀q₀ − p_v·q_v, −p₀q_v − q₀p_v − p_v × q_v). And q* ⊗ p* = (q₀p₀ − q_v·p_v, −q₀p_v − p₀q_v + q_v × p_v); since q_v × p_v = −p_v × q_v, they are equal. Then |pq|² = (pq)(pq)* = p q q* p* = p |q|² p* = |p|²|q|² (the real number |q|² commutes with everything).

(c) With q = (q₀, w) and the vector v: first q* ⊗ (0, v) = (w·v, q₀v − w × v). Then times q: the scalar part is (w·v)q₀ − (q₀v − w × v)·w = 0; the vector part is (w·v)w + q₀(q₀v − w × v) + (q₀v − w × v) × w = (w·v)w + q₀²v − 2q₀ w × v − (w × v) × w. With (w × v) × w = v|w|² − w(w·v): **v_b = (q₀² − |w|²)v + 2(w·v)w − 2q₀ w × v**, which is (5). With q₀ = cos(a/2), w = sin(a/2)n: q₀² − |w|² = cos a, 2wwᵀ = (1 − cos a)nnᵀ, 2q₀[w×] = sin a[n×]: **C = cos a E + (1 − cos a)nnᵀ − sin a[n×]**, the passive Rodrigues matrix.

(d) A turn δ about the body axes during Δt composes on the right: q(t + Δt) = q(t) ⊗ (cos(|ω|Δt/2), sin(|ω|Δt/2) ω/|ω|) = q(t) ⊗ (1, ωΔt/2) + O(Δt²). So [q(t + Δt) − q(t)]/Δt → **½ q ⊗ (0, ω)**. With (2), q ⊗ (0, ω) = (−w·ω, q₀ω + w × ω), and q · that = −q₀ w·ω + w·(q₀ω) + w·(w × ω) = **0**.

(e) |slerp|² sin²Ω = sin²((1 − u)Ω) + sin²(uΩ) + 2 sin((1 − u)Ω) sin(uΩ) cos Ω. Using 2 sin A sin B = cos(A − B) − cos(A + B) with A + B = Ω and the identity, this equals sin²Ω, so **|slerp| = 1**. And q₀ · slerp = [sin((1 − u)Ω) + sin(uΩ) cos Ω]/sin Ω = cos(uΩ) (expand sin((1 − u)Ω)), so the quaternion has moved through the angle uΩ on the sphere and the attitude through **2uΩ**.

**Common mistakes.** Treating the vector part of a pure-quaternion product as the cross product only (forgetting the −a·b scalar part); confusing the quaternion angle Ω with the attitude angle 2Ω.

## 3. Averaging star-tracker samples (s = 15)

**Assumptions.** Ten independent errors of 10⁻³ rad per axis (standard deviation); signs random. Numbers from `numpy.random.default_rng(15)`.

(a) The signs come out (−, +, +, −, +, +, +, −, −, −): five of each. The naive mean has norm **2.36 × 10⁻⁴** (almost zero) and its normalised value is **3.141 rad (180°)** from q_t: meaningless.

(b) After aligning signs with the first sample: error **2.58 × 10⁻⁴ rad**; the mean error of a single sample is **1.66 × 10⁻³ rad**.

(c) Eigenvector method: **2.58 × 10⁻⁴ rad**, the same to three figures. Σ qqᵀ is unchanged when any q is replaced by −q (qqᵀ = (−q)(−q)ᵀ), so the signs cannot matter.

(d) Independent errors partly cancel, so the error of the mean falls as 1/√N: **about 0.32 for N = 10**. The observed ratio, 2.58 × 10⁻⁴/1.66 × 10⁻³ = 0.16, is one random draw of the averaged error compared with a ten-sample mean of single errors, so it scatters around that value; repeating with other seeds shows the spread.

**Check.** With all signs equal (no flips), the naive average equals the aligned one.

**Common mistakes.** Averaging Euler angles (fine for tiny spreads, wrong near ±180° yaw); trusting a sign convention across tools.

## 4. Quaternion against DCM (s = 15: ψ̇ = 0.35, φ̇ = 3.0 rad/s, β = 17.5°)

(a) Over 60 s the two routes differ by at most **1.19 × 10⁻¹⁰ rad**; the norm drifts by at most **3.8 × 10⁻¹⁰**.

(b)

| Δt (s) | error (rad) | | \|q\| − 1 | |
|---|---|---|
| 0.04 | 3.265 × 10⁻⁵ | 9.22 × 10⁻⁷ |
| 0.02 | 2.043 × 10⁻⁶ | 2.88 × 10⁻⁸ |
| 0.01 | 1.277 × 10⁻⁷ | 9.00 × 10⁻¹⁰ |

Orders **4.00, 4.00** (error) and **5.00, 5.00** (norm).

(c) The `P06.07` DCM errors at the same steps were 5.090 × 10⁻⁴, 3.197 × 10⁻⁵ and 2.000 × 10⁻⁶ rad: ratios **15.6, 15.6, 15.7**. The quaternion's kinematics run at ω/2 (8), so its RK4 truncation error at a given step is that of a DCM at half the rate: (½)⁴ = 1/16.

(d) Integrating 1 s at Δt = 0.02 s: |q| − 1 reaches **4.8 × 10⁻¹⁰**, and C(q) is off orthogonal by **3.3 × 10⁻⁹** (≈ 6.9 × the norm error, as 2.2 predicts) before each renormalisation. Acceptable for almost any purpose; renormalising every step costs one division and removes the question.

**Common mistakes.** Comparing the raw (unnormalised) q with the DCM and attributing the 10⁻⁹ difference to the integration.

## 5. Interpolating between star-tracker reports (s = 15: h = 0.25 s)

(a) h = 0.25 s: SLERP **2.152 × 10⁻³ rad**, NLERP **2.154 × 10⁻³ rad**; h = 0.125 s: SLERP **5.347 × 10⁻⁴**, NLERP **5.348 × 10⁻⁴ rad**.

(b) Linear 3-2-1 angles: **1.133 × 10⁻² rad** at 0.25 s and **2.847 × 10⁻³ rad** at 0.125 s, five times worse than SLERP.

(c) Orders: SLERP **2.01**, NLERP **2.01**, angles **1.99**. All are linear interpolations of something, so all are second order; the angles have the largest constant because their rates vary more than the quaternion's along a coning motion (and they would break completely near pitch ±90°).

(d) **SLERP** (or NLERP when the turn per interval is guaranteed small): it is exact for a constant-rate turn, has no singularity, and costs a few operations; linear angle interpolation is five times less accurate here and singular at gimbal lock.

**Common mistakes.** Interpolating angles without unwrapping (a jump from +179° to −179° is interpolated through 0°).

## 6. Judgment

(a) Reading (x, y, z, w) as (w, x, y, z) puts q₃ into the scalar slot: the result is a different, usually large rotation: a small true turn (≈ (1, ε₁, ε₂, ε₃) stored as (ε₁, ε₂, ε₃, 1)) is read as nearly a 180° turn about axis 3. The SciPy comparison test, which explicitly reorders to scalar-last, catches the swap.

(b) The unit quaternions are the unit sphere in four dimensions, and dividing by the norm gives the nearest point on it, exactly. A 3 × 3 matrix has six orthonormality conditions, and the nearest rotation needs the polar decomposition (SVD), not a scale factor.

(c) Probably not. If the software keeps q₀ ≥ 0, the stored sign flips whenever the attitude passes a 180° turn from the reference; if it keeps continuity, a jump would indicate a fault. A jump from +0.7 to −0.7 with the vector part also changing sign is the same attitude; check the error angle between successive samples, which should be small.

(d) The direction of the matrix (TM-74839 maps rotated-frame components into the original frame, the transpose of C_b/n) and the quaternion's sign and product convention (which side the conjugate sits on). Check numerically: build random attitudes with `quat_from_euler_321`, evaluate the report's formula, and compare with `dcm_from_quat` and its transpose, and with ±q.

(e) When the turn per interval is small: here at 0.25 s (0.6 rad) SLERP and NLERP agree to 10⁻³ of the error; at 1 s (2.4 rad) NLERP is 60 % worse. NLERP is also adequate wherever the speed along the path does not matter (e.g. a display).

(f) Euler angles are readable (a pilot or an operator understands yaw, pitch, roll) and three numbers are the minimum; the quaternion costs a fourth number with a constraint to maintain, and the double cover must be handled in every comparison, average or interpolation.

## Check-yourself answers

1. i ⊗ j = k and j ⊗ i = −k: the product is not commutative. (Here k and −k happen to be the same attitude, a 180° turn about axis 3, because the factors are 180° turns; for general turns the two orders give different attitudes, as Problem 1 shows.)
2. (cos 45°, 0, 0, sin 45°) = (0.7071, 0, 0, 0.7071); its conjugate (0.7071, 0, 0, −0.7071) is the inverse turn, −90° about axis 3. The first describes a frame yawed 90° from the reference.
3. Because C(q) is quadratic in q. It matters in averaging, interpolation (shortest arc), differences between successive samples, and comparisons between tools.
4. q_c/n = q_b/n ⊗ q_c/b: the newer turn on the right.
5. d|q|²/dt = 2 q · q̇ = q · (q ⊗ (0, ω)) = −q₀ w·ω + q₀ w·ω + w·(w × ω) = 0.
6. The quaternion turns at half the body rate, and RK4's error scales with the fourth power of the rate times the step: (½)⁴ = 1/16.
7. It restores |q| = 1, so C(q) is a rotation again; it does not remove the attitude error the integration made.
8. When the turn per interval is large (around a radian or more): NLERP's speed along the arc is not uniform. SLERP needs the check because q and −q are the same attitude but opposite points on the sphere; without it the interpolation can take the long way round.
