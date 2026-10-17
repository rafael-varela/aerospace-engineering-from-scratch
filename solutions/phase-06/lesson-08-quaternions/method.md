# P06.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. q₃(90°) = (cos 45°, 0, 0, sin 45°), q₂(90°) = (cos 45°, 0, sin 45°, 0). 2. Yaw then pitch: q₃ ⊗ q₂ by (2); pitch then yaw: q₂ ⊗ q₃. 3. `axis_angle_from_quat` on each. 4. Row 1 of `dcm_from_quat(q)`. 5. `transform(q, [1, 0, 0])`.

## 2.
(a) pq = p₀q₀ + p₀q_v + q₀p_v + p_v q_v, and for pure quaternions a b = −a·b + a × b (from i² = −1, ij = k, ji = −k). (b) Conjugate each term of (2); then |pq|² = pq(pq)* = p(qq*)p* = |q|² pp*. (c) Compute q* (0, v) first, then multiply by q; collect terms in v, (q_v·v)q_v and q_v × v; write as a matrix; substitute the half-angle forms. (d) q(t + Δt) = q(t) ⊗ δq, δq = (cos(|ω|Δt/2), sin(|ω|Δt/2) ω/|ω|) ≈ (1, ωΔt/2); subtract and divide by Δt. q · (q ⊗ (0, ω)) by (2). (e) |slerp|² with q₀ · q₁ = cos Ω and the identity sin²((1 − u)Ω) + sin²(uΩ) + 2 sin((1 − u)Ω) sin(uΩ) cos Ω = sin²Ω; the dot product of q₀ with slerp is cos(uΩ).

## 3.
(a) `meas.mean(axis=0)`, its norm, `quat_error_angle(quat_normalize(mean), q_t)`. (b) Flip samples with `m @ meas[0] < 0`, average, normalise. (c) `np.linalg.eigh(sum(np.outer(m, m) for m in meas))`, last column. (d) Independent errors: standard error falls as 1/√N.

## 4.
(a) Both integrations with `n_samples=601`; `misalignment_angle(dcm_from_quat(quat_normalize(Q[m])), Cs[m])`; `np.linalg.norm(Q, axis=1)`. (b) `integrate_quat_rk4(q0, w, 60, dt)`, `quat_error_angle` against `quat_from_dcm(coning_dcm(t[-1], …))`, `quat_norm`. (c) Divide the `P06.07` errors by these. (d) Integrate 1 s at Δt = 0.02 s; norm drift and `dcm_from_quat` orthogonality.

## 5.
(a) `ts = np.arange(0, 30 + 1e-9, h)`; samples from `quat_from_dcm(coning_dcm(t, …))`; `interpolate_samples` both ways; max `quat_error_angle` against the truth on `np.linspace(0, 30, 3001)`. (b) Unwrapped angles, `np.interp` per angle, back to q. (c) log₂ of the ratio between h and h/2. (d) Accuracy, cost, singularity.

## 6.
(a) The element order permutes the axis and turns the angle into a different one; the SciPy comparison test. (b) The unit sphere: dividing by the norm is the nearest point; the 3 × 3 orthogonality needs six conditions. (c) Sign convention q₀ ≥ 0 or a crossing through q₀ = 0. (d) Matrix direction (the report maps rotated to original frame) and the quaternion's sign or handedness; compare with `dcm_from_quat` on random angles. (e) Small turns per interval. (f) Human readability; four numbers and a norm constraint, the double cover.
