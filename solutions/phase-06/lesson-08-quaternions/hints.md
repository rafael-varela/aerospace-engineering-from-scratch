# P06.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. q_b/n = q₃(90°) ⊗ q₂(90°) for yaw-then-pitch, with the newer turn on the right; row 1 of C(q) is the nose in reference components.
2. Write p = p₀ + p_v and expand p_v q_v with i² = −1 and ij = k; the product of two pure quaternions is (−a·b, a × b). For (c), use the identity a × (b × c) = b(a·c) − c(a·b).
3. q and −q average to zero; any operation that is quadratic in q (qqᵀ) cannot see the sign.
4. `integrate_quat` and the given `integrate_dcm`; the quaternion turns at half the body rate, so its RK4 error scales with (ωΔt/2)⁴.
5. `interpolate_samples(ts, qs, t_fine, "slerp")`; for the angles, `euler_321_from_quat`, `np.unwrap(…, axis=0)`, `np.interp`, then `quat_from_euler_321`.
6. Element order; what the nearest unit quaternion is; the double cover; matrix direction and quaternion sign conventions; the size of the turn between samples; readability and the constraint.
