# P09.18  Hints

One hint per problem. Read only the one you need.

1. The disk bounds are formulas (4) and (5) of the lesson. A first-order lag peaks at ω = 0; a second-order resonance with ζ < 1/√2 peaks at ω_n√(1 − 2ζ²) with height 1/(2ζ√(1 − ζ²)). Convert a gain ratio to dB with 20 log₁₀.

2. (a) Solve y = Gu + d, u = C(r − y − n) for y first. (b) |1 + L| = 1/|S|. Draw the disk of radius 1/Ms about −1 and the unit circle about 0; the phase margin is the angle, seen from the origin, between −1 and the nearest intersection point. (c) Factor L₀ out of 1 + L₀ + L₀WΔ. For necessity, an all-pass c(a − s)/(a + s) has constant magnitude |c| and a phase you can tune with a. (d) Write (jωI − H(γ))v = 0 in two halves and eliminate one of them.

3. Split the integral into [0, ω₁], [ω₁, ω₂] and [ω₂, ∞) and bound ln|S| by ln ε, ln Ms and 0 on the three pieces.

4. `lead(Kc, 20, 80)`, `closed_loop_tfs`, `tf_hinf_norm` on S, T and KS; `tf_to_ss` then both norm functions; `robust_stability_index`; `margins` on L₀(1 + Δ_true) with `true_servo_relative_error` on a grid starting at 10⁻⁴ rad/s.

5. Read the test as "for every Δ in the ball". Copy `true_servo_relative_error` with 120 rad/s and call `weight_covers`. For the real gain, scale the plant numerator by k and test stability over a range of k.

6. Bisect K_c on `robust_stability_index(...) < 1` and on `is_closed_loop_stable(..., *true_servo_plant(4))`. For (d), compare Δ_true(jω*)/W(jω*) with the destabilising value −1/(WT(jω*)), both in magnitude and in phase.
