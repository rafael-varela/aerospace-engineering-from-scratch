# P09.19  Hints

One hint per problem. Read only the one you need.

1. Multiply out AᵀP + PA for a symmetric P = [[p₁₁, p₁₂], [p₁₂, p₂₂]]: the (1,1) entry contains only p₁₂, the (2,2) entry p₁₂ and p₂₂, the (1,2) entry all three. Solve in that order.

2. (a) Differentiate e^{Aᵀt}Qe^{At} and integrate from 0 to ∞. (b) For A > δ the output is k·A sin θ until sin θ = δ/A and kδ after; integrate y sin θ over a quarter period and use symmetry. (c) The output is a square wave shifted by the angle arcsin(h/A). (d) Solve the pendulum equation for u with θ̈ replaced by v.

3. The phase of K/(s(τs + 1)) is −90° − arctan(ωτ). For an ideal relay, try simulating and watch the switching interval shrink. Each period of the cycle fires each thruster once.

4. The critical gain is 1/|G(j50)| = 1/0.24. For the region of attraction, reuse the lesson's code with your Q and seed 19 for the random starts, sampling uniformly in the ellipse (radius ∝ √uniform).

5. Ask what the first-harmonic approximation throws away and when that is large. Compare the shape of an ellipse with the shape of the true region near the saturated corners. Look at the peak torque printed by `main()`. Linearise the mismatched closed loop at θ = 0 and check the signs of its coefficients.

6. Pass a `dict(p, m=1.1 * p["m"])` as the plant parameters while the controller uses the nominal ones. Wrap each controller in `min(6, max(-6, u))` and bisect on the initial angle.
