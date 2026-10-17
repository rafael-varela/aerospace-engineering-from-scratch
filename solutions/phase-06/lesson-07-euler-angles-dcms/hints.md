# P06.07  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. `dcm_from_euler([Ω, i, ω], "313")`; v_p = C v_i, so the inertial components of a perifocal axis are a row of C.
2. Each column of B is a unit vector; the angle between the first and last columns decides the triple product. For (c), differentiate CCᵀ and name the skew-symmetric matrix you get.
3. Row i of a passive DCM is body axis i written in the reference frame. Loop over `SEQUENCES` with `euler_from_dcm` and `euler_rate_matrix`; then repeat for `dcm_from_euler([0, θ, 0], "321")`.
4. `coning_body_rate`, `integrate_dcm`, `integrate_dcm_rk4` with and without `renormalize_every=1`; the order is log₂ of the error ratio for each halving.
5. `hammerhead(0.5, 0.5, radians(90 − δ))` and `integrate_segments_rk4(…, 1/f, "euler")` and `"dcm"`; compare with the returned closed-form attitude.
6. What each representation is good at; the parity sign; what the nearest rotation can and cannot know; the convention of the source; a mechanical gimbal against a strapdown computer; disturbances and the cost of being wrong.
