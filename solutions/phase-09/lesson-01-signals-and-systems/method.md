# P09.01  Method

The route to each answer, without the final numbers. `worked-solutions.md` has the numbers for s = 15.

## 1. Energy, power, even and odd

(a) E = (s/10)² ∫₀^∞ e^{−4t} dt. (b) (3 cos θ + 1)² = 9 cos² θ + 6 cos θ + 1; the means over whole periods are 9/2, 0 and 1. (c) For t > 0: x_e = x_o = ½e^{−t}; for t < 0: x_e = ½e^{t}, x_o = −½e^{t}; each energy is 2 × ¼ ∫₀^∞ e^{−2t} dt. (d) `energy(x, dt)` on t = arange(0, 20, 1e-4); `average_power` over 10 s (a whole number of periods for integer s); `even_part`/`odd_part` on a symmetric grid such as linspace(−20, 20, 400001). The left Riemann sum of a decreasing function overestimates by a relative amount ≈ (rate) × dt.

## 2. Four proofs

(a) |y[n]| ≤ Σ_k |h[k]||x[n − k]| ≤ M Σ|h|. Necessity: x[m] = sign(h[N − m]) for 0 ≤ m ≤ N gives y[N] = Σ_{k=0}^{N} |h[k]|. (b) S(x₁ + x₂) = x₁ + x₂ − b ≠ x₁ + x₂ − 2b; S(x₁) − S(x₂) = x₁ − x₂. (c) S(δ) = δ (the sample at n = 0 is kept); S(δ[n − 1]) = 0 (x[1] is never read at an integer output index); a delayed input does not give a delayed output. y[1] = x[2]: the output at n = 1 needs the input at n = 2. (d) g(t) = x(t − b), g(at) = x(at − b); h(t) = x(at), h(t − b/a) = x(at − b); h(t − b) = x(at − ab).

## 3. Star tracker

Median: compare median(x₁ + x₂) with median(x₁) + median(x₂) for the hint's pair; time-invariant (the same rule at every n); causal (uses n, n − 1, n − 2); has memory; bounded (the median lies between the smallest and largest input). Deadband: homogeneity fails (scale a sample across d); time-invariant; causal; memoryless; |y| ≤ |x|. Code: write each as a function of the array (pad with two zeros for the median, `np.where` for the deadband) and call `classify`. A transfer function needs LTI: the median is not linear at any amplitude (it is a rank filter); the deadband is linear only for inputs that never cross d, i.e. for |x| ≥ d always (identity) or |x| < d always (zero).

## 4. Bias drift and record length

(a) `accumulator(1/100)(np.full(N, b))[-1]` with N = 100 × T; by hand b × T. (b) Loop n over 10 · 2^k, k = 0..7; `is_bibo_stable(first_order_recursion(a, 1 − a), n=n)`; first True. Divide by 10 s. (c) Print, for the failing n, the ratio peak(4n)/peak(n) for each input of `bounded_test_inputs`: the step's ratio is 1 + aⁿ + a²ⁿ + a³ⁿ; the random-sign input's output wanders slowly and its peak keeps growing until the record is several time constants long.

## 5. Judgment

One fooling system per tester: linearity (a saturation, or any block that is linear only below an amplitude larger than the test amplitude); time invariance (a gain switched at a time beyond the record); causality (a block that looks ahead only from some index onward that the test never perturbs, or one whose look-ahead is longer than the gap between perturbation points if the test perturbed sparsely); memory (a block whose dependence on the past is only triggered by large inputs); stability (a slow pole near 1, or a resonance at a frequency none of the inputs contains). Fixes: amplitude sweeps up to the operating limit, records longer than the slowest time constant and the longest schedule, broadband and swept-frequency inputs, and analysis of the code itself. The paragraph should say that tests can falsify but not certify, and that LTI blocks are certified by their poles and coefficients.

## 6. Resonator

(a) Substitute h into the recursion using the identity in the hint; h[0] = 1, h[1] = 2r cos θ. (b) Σ|h| ≤ (1/|sin θ|) Σ rⁿ < ∞ for r < 1; for r = 1, |sin((n + 1)θ)| does not tend to zero, so the sum diverges. `abs_sum(impulse_response(...))` over 5000 samples. (c) With r = 1 the characteristic roots are e^{±jθ} on the unit circle; an input at the same frequency gives a particular solution of the form n sin(θn + φ) (resonance); a step gives a constant particular solution 1/(2 − 2 cos θ) plus a bounded oscillation. (d) Run `is_bibo_stable` and, for each of the four inputs, the ratio peak(8000)/peak(2000).
