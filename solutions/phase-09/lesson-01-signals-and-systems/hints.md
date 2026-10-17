# P09.01  Hints

One hint per problem. Read only the one you need, then `method.md` if the hint is not enough.

1. ∫₀^∞ e^{−4t} dt = 1/4. Expand (3 cos θ + 1)² and average each term over whole periods. For the even part of e^{−t}u(t), write the two halves t > 0 and t < 0 separately.

2. (a) Triangle inequality for sufficiency; for necessity pick x[m] = sign(h[N − m]) and look at y[N]. (b) Compare S(x₁ + x₂) with S(x₁) + S(x₂); then compute S(x₁) − S(x₂). (c) Feed δ[n] and δ[n − 1] and compare the outputs; for causality look at which input index y[1] reads. (d) Substitute: if g(t) = x(t − b), what is g(at)?

3. A median of three numbers is not a weighted sum of them: try x₁ = (0, 0, 1) and x₂ = (0, 1, 0). A deadband is the identity for large inputs and zero for small ones, so test homogeneity with a scale factor that moves a sample across the threshold.

4. Angle = bias × time. For (b) predict from the step response 1 − aⁿ: by how much does its peak still grow between n and 4n? Then look at which of the four test inputs actually trips the 1.5 threshold.

5. For each tester, ask what the test never tries: amplitudes, times after the record ends, inputs that excite a particular frequency, look-ahead at an index that is not perturbed.

6. (a) Check that h satisfies the recursion with h[−1] = h[−2] = 0 using the identity sin((n + 1)θ) = 2 cos θ sin(nθ) − sin((n − 1)θ). (b) |h[n]| ≤ rⁿ/|sin θ|; for r = 1 the terms do not decay. (c) Resonance: the input frequency equals the frequency of the undamped mode. (d) Print the peak ratio for each of the four inputs separately.
