# P02.07  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. r = √(x² + y²), θ = atan2(y, x). 2. z⁵ has modulus r⁵ and argument 5θ; reduce modulo 2π into (−π, π]; `z ** 5` for the check. 3. `nth_roots(z, 3)`: modulus r^{1/3}, arguments (θ + 2πk)/3; `np.prod`, `np.sum`. 4. The roots are r^{1/3} e^{iθ/3} times the three cube roots of unity, which sum to zero.

## 2.
(a) Substitute z = iθ in Σ z^k/k!, separate even and odd k using i² = −1, recognise the two series. (b) Expand both sides of e^{i(a+b)} = e^{ia}e^{ib} with Euler's formula and compare real and imaginary parts. (c) Write F₀ cos ωt = Re(F₀ e^{iωt}); substitute x = Re(X e^{iωt}); since m, c, k are real constants, m Re(w) = Re(m w) etc., and the derivatives bring down iω and (iω)²; collect the coefficient of e^{iωt}; equality of real parts for all t gives equality of the coefficients (evaluate at two instants a quarter period apart). (d) With u = r², D(u) = (1 − u)² + 4ζ²u; dD/du = −2(1 − u) + 4ζ² = 0 gives u = 1 − 2ζ²; D at that point is 4ζ²(1 − ζ²); |H|_peak = 1/(2ζ√(1 − ζ²)).

## 3.
(a) Z(t) = 1.2 + 0.9 e^{i2πΔf t}. (b) |Z(t)| by the cosine rule with δ = 2πΔf t; the envelope repeats every 1/Δf. (c) `t = np.linspace(0, 2, 200001)`, build x, loop over windows, `fit_sinusoid(t[w], x[w], 2*np.pi*40)`; compare with `combined_amplitude(1.2, 0.9, 2*np.pi*df*(t0 + 0.05))`. (d) The amplitude sweeps through its full range once per beat period; the synchrophaser sets Δf = 0 and freezes Z at the chosen δ.

## 4.
(a) `add_sinusoids([(1.2, 0), (0.9, np.radians(180 - s))])`; displacement = A/ω². (b) `transmissibility(40/30, 0.08)` times A. (c) `mount_frequency_for_transmissibility(0.2, 40, 0.08)`; k = m ω_n², c = 2ζ m ω_n; sag = g/ω_n². (d) `transmissibility(80/f_new, zeta)`, `transmissibility(20/f_new, zeta)`; for the ODE, rhs returns (v, (c ẏ + k y − c v − k x)/m) with y = Y cos ωt, ẏ = −Yω sin ωt; `solve_ivp` to 3 s at rtol 1e-10; `fit_sinusoid` on t > 2 s; ratio A_fit/Y.

## 5.
(a) `combined_amplitude(1.2, 0.9, pi)` and `combined_amplitude(1.0, 1.0, pi + pi/3)`. (b) Sweep δ from 0 to 360° in 1° steps; location 1 minimum at δ = 180°, location 2 at δ + 60° = 180°; the sum's minimum by `np.argmin`. (c) Blade-passing (80 Hz) content and broadband noise; the first adds a second phasor problem with its own best δ, the second does not shift the 40 Hz minimum but makes it harder to find. (d) Compare the sum at the compromise angle with the sum at each location's own optimum.

## 6.
(a) Geometric series with ratio q = e^{2πi/N} ≠ 1: sum = (1 − q^N)/(1 − q) = 0 since q^N = 1. (b) With δ = 2πh/N the formula gives sin(πh)/sin(πh/N): the numerator is always zero; the denominator is zero only when h/N is an integer, where the sum is N by direct count. `add_sinusoids([(1.0, 2*pi*h*k/N) for k in range(N)])`. (c) `add_sinusoids([(1.05, 0), (1, 2*pi/3), (1, 4*pi/3)])`; subtract the balanced sum (zero) to see that only the excess 0.05 at phase 0 survives. (d) As in the hint; take the modulus at the end.
