# P02.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. T = N/f_s; Δf = 1/T; f_Nyq = f_s/2; bins = N/2 + 1 (N even). 2. a cos(2πf_k t + φ) = ½a e^{iφ}e^{2πikn/N} + c.c.; orthogonality gives |X[k]| = |X[N − k]| = Na/2; `amplitude_spectrum` doubles and divides by N. 3. Rectangular main lobe is 2Δf wide (first zeros at ±Δf), so separation needs Δf ≤ 0.1 Hz → T ≥ 1/0.1; the Hann main lobe is 4Δf, so T doubles.

## 2.
(a) Set r = e^{2πi(m−n)/N}; if m = n the sum is N; otherwise r ≠ 1, rᴺ = 1, sum = 0. Multiply (1) by e^{2πikm/N}, sum over k, swap sums, apply the identity. (b) Σ_n |x[n]|² = Σ_n x[n] x̄[n]; insert (2) for both factors; the n-sum of the exponentials gives Nδ_{kl}; the result is (1/N)Σ|X[k]|². (c) Σ_{n=0}^{N−1} e^{2πi δ n/N} = (1 − e^{2πiδ})/(1 − e^{2πiδ/N}); magnitudes via |1 − e^{iθ}| = 2|sin(θ/2)|; for small argument sin(πδ/N) ≈ πδ/N. (d) w = −i; X[k] = Σ_n x[n] w^{kn} with only n = 0 and n = 2 nonzero: X[k] = 1 − (−i)^{2k} = 1 − (−1)^k.

## 3.
(a) f = rpm/60 for 1P; × blades for blade passing. (b) For each, r = f mod f_s, alias = min(r, f_s − r); compare with f_s/2. (c) An alias is indistinguishable from a real tone at the same frequency; the anti-aliasing low-pass filter before the ADC removes content above f_s/2. (d) f_s ≥ 2 × 1.2 × f_max.

## 4.
(a) Δf = 1/T; nearest bin = round(f₀/Δf)Δf; δ = |f₀ − nearest|/Δf. (b) Factor sin(πδ)/(πδ); run `amplitude_spectrum(x, 1/1024)` on 2 cos(2πf₀t) and take the maximum. (c) Same with `window="hann"`. (d) The periodic extension changed; the window-dependent bias must be reported with the number.

## 5.
(a) Geometric series: Σ(1 − α)αⁿ = 1; H(f) = (1 − α)Σ(α e^{−2πif/f_s})ⁿ. (b) `y = convolve_fft(a, h)[:a.size]`; analyse `y[256:]` and `a[256:]` with the Hann window; read the three peaks. (c) Ratio after/before against |H(f₀)|; against 1/√(1 + (f/f_c)²), which the discrete filter approximates when f_c ≪ f_s. (d) `np.max(np.abs(convolve_fft(a, h) − np.convolve(a, h)))`; time both with `time.perf_counter`.

## 6.
(a) f_s ≥ 2 × f_max × margin; filter corner between f_max and f_s/2 with enough roll-off. (b) T from Δf = 0.5 Hz and the window's main-lobe width; bias from the worst-case δ for the chosen window. (c) Drift × T against Δf: the tone must stay within a bin during the record. (d) Zero-padding evaluates the same windowed spectrum at more frequencies; the main-lobe width in hertz is unchanged, so separation is not improved.
