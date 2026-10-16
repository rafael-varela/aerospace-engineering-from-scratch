# P02.10  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything follows from T = N Δt: Δf = 1/T, f_Nyquist = f_s/2, bins = N/2 + 1. For (b), write the cosine as two complex exponentials and use orthogonality; for (c), two tones are separable when they are at least one main-lobe width apart.
2. For (a) the geometric series Σ rᵏ = (1 − rᴺ)/(1 − r) with r = e^{2πi(m−n)/N} has rᴺ = 1. For (b) substitute (2) into Σ|x[n]|² and swap the order of summation. For (c) the sum is again geometric with ratio e^{2πiδ/N}; use |1 − e^{iθ}| = 2|sin(θ/2)|. For (d) write out the four terms of (1) with w = e^{−2πi/4} = −i.
3. Blade-passing frequency = (rpm/60) × blades. Anything above f_s/2 folds: apply (6) and check it lands in [0, f_s/2]. The hardware fix is in front of the sampler, not after it.
4. δ is the distance from f₀ to the *nearest* bin, in units of Δf; it changes with record length because the bin grid changes. The rectangular factor is sin(πδ)/(πδ) times the true amplitude.
5. Σ(1 − α)αⁿ is a geometric series; H(f) is the DFT-style sum Σ h[n] e^{−2πifn/f_s} of the same series. The first 256 output samples contain the filter's start-up transient (the time constant is 1/(1 − α) ≈ 4 samples, so 256 is generous); windowing a transient would leak.
6. Each requirement maps to one number from part 2: Nyquist to f_s, resolution to T, amplitude bias to the window, drift to the product of T and the drift rate. Zero-padding interpolates; it does not add information.
