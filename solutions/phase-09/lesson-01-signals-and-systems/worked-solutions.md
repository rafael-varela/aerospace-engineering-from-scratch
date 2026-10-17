# P09.01  Worked solutions

Worked for **s = 15**. Every number below is produced by `solution.py` (its `main()` or the functions named) or by a probe script that calls them; the "independent check" line says how it was verified by a different route.

## 1. Energy, power, even and odd

**Assumptions.** Continuous signals sampled at dt = 10⁻⁴ s; energy by the left Riemann sum.

**Derivation and result.** (a) E = (1.5)² ∫₀^∞ e^{−4t} dt = 2.25/4 = **0.5625** (units of x² × s). (b) (3 cos θ + 1)² = 9 cos² θ + 6 cos θ + 1, whose averages over whole periods are 4.5, 0 and 1: **P = 5.5**. The cross term is a sinusoid and averages to zero over whole periods; the power of a sum is the sum of the powers when the components are orthogonal (here: different frequencies, 0 and 15 Hz). (c) x_e(t) = ½e^{−|t|}, x_o(t) = ½ sgn(t) e^{−|t|}; each has energy 2 × ¼ × ½ = **0.25**, together the 0.5 of x. (d) Code: 0.56261 for (a), 5.5 for (b), and 0.25008 and 0.24998 for the two parts of (c) (their sum 0.50005). The integrand of (a) is 2.25e^{−4t}, and its left Riemann sum is 2.25 dt/(1 − e^{−4dt}) ≈ 0.5625(1 + 2dt): a relative error of 2dt = 2 × 10⁻⁴, first order in dt, exactly what the code shows.

**Independent check.** The even and odd energies add to the total, as orthogonality requires (0.25008 + 0.24998 = 0.50005 against E(x) = 0.50005 on the same grid).

**Common mistakes.** Forgetting to square the amplitude; averaging over a non-integer number of periods; putting the factor ½ on only one side of t = 0.

## 2. Four proofs

**(a)** Sufficiency: if |x[m]| ≤ M for all m, |y[n]| = |Σ_k h[k]x[n − k]| ≤ Σ_k |h[k]||x[n − k]| ≤ M Σ_k |h[k]|, a bound valid for every n. Necessity: if Σ|h| = ∞, take any N and the bounded input x[m] = sign(h[N − m]) for 0 ≤ m ≤ N (zero elsewhere, |x| ≤ 1). Then y[N] = Σ_{k=0}^{N} h[k] sign(h[k]) = Σ_{k=0}^{N} |h[k]|, which exceeds any bound as N → ∞: no single bound M′ works for all bounded inputs, so the system is not BIBO stable.
**(b)** S(x₁ + x₂) = x₁ + x₂ − b, but S(x₁) + S(x₂) = x₁ + x₂ − 2b: they differ by b ≠ 0 (and S(0) = −b ≠ 0). But S(x₁) − S(x₂) = x₁ − x₂: the response to a *change* of input is linear in that change, which is what "incrementally linear" means, and why a bias can be handled by linearising around an operating point.
**(c)** Time variance: S(δ)[n] = δ[2n] = δ[n] (only n = 0 survives), but S(δ[n − 1])[n] = δ[2n − 1] = 0 for every integer n. A delayed input gives zero output, not a delayed output. Non-causality: y[1] = x[2], so the output at n = 1 depends on the input at n = 2.
**(d)** Let g(t) = x(t − b); then g(at) = x(at − b). Let h(t) = x(at); then h(t − b/a) = x(a(t − b/a)) = x(at − b). And h(t − b) = x(at − ab) ≠ x(at − b) unless a = 1 or b = 0. The test of 2.3 checks the same three cases numerically.

**Common mistakes.** Proving necessity with "an unbounded h gives an unbounded output": h can be bounded (the accumulator's h = u is bounded) and still not summable. Calling the bias correction linear because "it is a straight line".

## 3. Star tracker

**Hand verdicts.**

| block | linear | TI | causal | memoryless | BIBO |
|---|---|---|---|---|---|
| median of 3 | no: median(0,1,1) = 1 but median(0,0,1) + median(0,1,0) = 0 | yes: same rule at every n | yes: uses n, n−1, n−2 | no | yes: the median lies between the min and max input |
| deadband d | no: x = 0.4 gives 0 but 2x = 0.8 gives 0.8 ≠ 2·0 | yes | yes | yes | yes: \|y\| ≤ \|x\| |

**Code.** `classify` on the two callables returns, for the median, `linear False, time_invariant True, causal True, memoryless False, bibo_stable True`, and for the deadband with d = 0.5, `False, True, True, True, True`: identical to the hand table.

**(c)** Neither block has a transfer function in general. The median filter is a rank filter and is not linear at any amplitude (it is linear only on monotone input sequences, where it acts as a one-sample delay). The deadband is the identity when |x| ≥ d at every sample and zero when |x| < d at every sample; a frequency-domain analysis that treats it as a unit gain is valid only for errors well above 0.5 arcsec. Near the threshold it behaves like a gain that depends on amplitude (a describing function, `P09.19`), and it can cause the slow limit cycles typical of attitude controllers with deadbands.

**Common mistakes.** Testing linearity of the deadband only with large inputs (it then looks like a wire); calling the median filter "a kind of moving average".

## 4. Bias drift and the record-length trap

**(a)** b = 0.002 × 15 = **0.03 deg/s**. After 10 min: 0.03 × 600 = **18.00 deg**; after 1 h: **108.00 deg**. `accumulator(0.01)` on a constant 0.03 at 100 Hz gives 18.00 and 108.00 (to 10⁻⁹).

**(b)** a = 1 − 1/150, time constant 150 samples.

| n | 10 | 20 | 40 | 80 | 160 | 320 | 640 | 1280 |
|---|---|---|---|---|---|---|---|---|
| `is_bibo_stable` | False | False | False | False | False | False | **True** | True |

The smallest correct n in the sweep is **640 = 4.3 time constants**.

**(c)** The tester compares the output peak over 4n samples with the peak over the first n. For the step input the ratio is (1 − a^{4n})/(1 − aⁿ) = 1 + aⁿ + a²ⁿ + a³ⁿ, which falls below 1.5 once aⁿ < 0.39, about one time constant: the step alone gives a false alarm only up to n ≈ 160 (ratio 1.501 at n = 160, 1.133 at 320). The false alarm at n = 320 comes from the random-sign input: the filtered noise wanders on the scale of the time constant, so its peak keeps growing (ratio 1.865) until the record is several time constants long (ratio 1.000 at n = 640).

**Independent check.** Σ|h| = (1 − a)Σaⁿ = 1: the filter is stable whatever the test says.

**Common mistakes.** Reading the verdict at n = 320 as evidence of instability; expecting the threshold from the step response alone.

## 5. Judgment

- **Linearity**: a ±1 saturation passes `is_linear` at amplitude 0.01; fix: sweep the amplitude up to the operating limit (or beyond, to know the margin).
- **Time invariance**: y[n] = x[n] for n < 10 000 and 2x[n] after (a gain schedule that switches at 100 s at 100 Hz) passes on a 256-sample record; fix: records covering every mode and schedule of the software, or inspection of the code for any dependence on n or time.
- **Causality**: a block that reads ahead only when the input exceeds a threshold (an outlier rejector using the next sample) passes with Gaussian inputs of small amplitude; fix: perturb with large values at many indices.
- **Memory**: a block that stores the maximum input seen so far and uses it only once it exceeds a limit (a latched fault flag) looks memoryless with small inputs; fix: perturbations of every size.
- **Stability**: the low-pass with a = 0.999 is called unstable on a 50-sample record (false alarm), and a resonator with a pole on the unit circle at a frequency absent from the test inputs can be called stable (missed); fix: records many time constants long and broadband or swept inputs, and, for LTI blocks, the pole test.

Reply to the reviewer: "`classify` can prove that a property fails, by exhibiting a counterexample; it cannot prove that it holds, because it tries finitely many inputs at finitely many amplitudes over a finite record. For a filter whose code is a linear difference equation, linearity and time invariance follow from the code itself, and stability follows from its poles (all strictly inside the unit circle); those are the facts to certify. The numerical tests are useful regression checks that would catch an implementation error, provided their amplitudes and record lengths cover the flight envelope and the slowest time constant."

## 6. A marginally stable resonator

**(a)** With h[−1] = h[−2] = 0: h[0] = 1 = sin θ/sin θ; h[1] = 2r cos θ = r sin 2θ/sin θ. For n ≥ 2, substituting h[n] = rⁿ sin((n + 1)θ)/sin θ into 2r cos θ h[n − 1] − r² h[n − 2] gives (rⁿ/sin θ)[2 cos θ sin(nθ) − sin((n − 1)θ)] = rⁿ sin((n + 1)θ)/sin θ by the identity in the hint. The code's impulse response matches the closed form to 1.3 × 10⁻¹⁴ for r = 0.95, θ = 2π/50.

**(b)** |h[n]| ≤ rⁿ/|sin θ|, so Σ|h| ≤ 1/((1 − r)|sin θ|) < ∞ for r < 1. For r = 1 the terms |sin((n + 1)θ)|/|sin θ| do not tend to zero (for θ not a multiple of π they exceed a fixed positive number infinitely often), so the sum diverges. For r = 0.95, θ = 2π/50: **Σ|h| = 101.12** over 5000 samples (the bound gives 1/(0.05 × 0.1253) = 159.6).

**(c)** For r = 1 the characteristic equation z² − 2 cos θ z + 1 = 0 has roots e^{±jθ} on the unit circle. The input cos(θn) has the same frequency, so the particular solution is of the form n·A sin(θn + φ): the output grows linearly (resonance). A step has frequency 0 ≠ θ; its particular solution is the constant 1/(2 − 2 cos θ), and the homogeneous part is a bounded oscillation. For θ = 2π/50 that constant is 63.4 and the step response oscillates between 0 and 126.8; the code's peak is 127 at both 2000 and 8000 samples.

**(d)**

| (r, θ) | `is_bibo_stable` | peak ratio (8000 vs 2000 samples): step, ±1, random signs, sine of period 50 |
|---|---|---|
| (0.95, 2π/50) | True | 1.00, 1.00, 1.00, 1.00 |
| (1, 2π/50) | **False** | 1.00, 1.00, 1.00, **3.99** |
| (1, 2π/37) | **False** | 1.00, 1.00, **1.94**, 1.00 |

For θ = 2π/50 the slow sine of `bounded_test_inputs` happens to sit exactly at the resonance and its output grows as n (ratio 4 for 4 times the record). For θ = 2π/37 no test input is at resonance, and only the random-sign input, which has power at every frequency, reveals the growth (a random walk at the resonant frequency grows like √n: ratio about 2 for 4 times the record). The step and the alternating input never reveal it. A test-input set must be broadband, or swept in frequency, to find a marginal mode; a step test alone would have certified this resonator as stable.

**Common mistakes.** Concluding from a bounded step response that the resonator is stable; forgetting that a finite sum over 5000 samples is not Σ|h| for r = 1 (it grows with the number of samples).

## Check-yourself answers (lesson.md)

1. E = 1/(2 × 3) = 1/6 (energy signal; its power is zero). P = 25/2 = 12.5 (power signal; its energy is infinite).
2. x_e = ½e^{−|t|}, x_o = ½ sgn(t)e^{−|t|}. Their product is odd, so ∫x_e x_o dt = 0, and E(x) = ∫(x_e + x_o)² = E(x_e) + E(x_o) + 2∫x_e x_o = E(x_e) + E(x_o).
3. x(3t − 6) = x(3(t − 2)): after scaling by 3, shift right (delay) by 2.
4. Not linear: S(0) = 2 ≠ 0 (it is affine). Time-invariant: the rule is the same at every n, and a delayed input gives a delayed output.
5. Its output at index n is the input at index 2n, which for n ≥ 1 lies after n. Run in real time at the output rate, the block would need inputs that have not yet arrived; reading a stored file hides this because the whole file exists.
6. Σ_n |h[n]| < ∞. The accumulator has h[n] = u[n] = 1 for n ≥ 0, so Σ|h| = ∞ (and the step input gives a ramp).
7. No: Σ|h| = 1, the pole 0.999 is inside the unit circle. The record is shorter than the time constant, so the step response is still rising; use a record several time constants long (at least about 4000 samples here) or decide from the pole.
8. 0.005 × 1200 = 6 deg. The integrator (accumulator) is not BIBO stable: a bounded constant input produces an output that grows linearly in time.
