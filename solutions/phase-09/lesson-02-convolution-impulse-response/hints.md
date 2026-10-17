# P09.02  Hints

One hint per problem. Read only the one you need, then `method.md` if the hint is not enough.

1. Write h reversed under x and slide it one place at a time; there are 3 + 3 − 1 positions. For a circular length L, fold y: add y[n + L] onto y[n].

2. (a) Each strip of the staircase is a short pulse of area x(τ_k)Δ; what does an LTI system do with a short pulse? (b) Combine the exponentials in the integrand into e^{−bt} e^{(b−a)τ}. Laplace: 1/((s + a)(s + b)) by partial fractions. (c) Integrate t e^{−t/τ} by parts. (d) Compare the rectangle sum with the trapezoid rule term by term: the difference is the two end-point half weights.

3. The input is a ramp r(t) = 2t for t ≤ 10 s. The response of a first-order lag to a ramp of slope c is c[t − τ(1 − e^{−t/τ})]. After the ramp stops, the remaining lag decays like e^{−(t−10)/τ}.

4. The FFT length L = next power of two ≥ N + M − 1 does not change while N + M − 1 stays below the same power of two. For the circular error, the wrap-around region is the first M − 1 samples.

5. (a) The fastest rate of a first-order response is at the start of a step: (size of the step)/τ. At the reversal of a doublet the step is twice the amplitude. (b) Evaluate both expressions with dt = 0.02 s and τ = 0.05 s. (c) ∫₀^T h/∫₀^∞ h = 1 − e^{−T/τ}. (d) Think of what in the actuator depends on temperature.

6. (a) Split ∫₀^{t_k} x(τ)h(t_k − τ)dτ into the strips where x is constant and integrate h analytically on each. (c) The formula makes no approximation when its assumption holds. (d) What does the formula assume about x inside a strip, and what if x is a smooth curve instead?
