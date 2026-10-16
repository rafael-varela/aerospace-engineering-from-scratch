# P02.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. L = (40 + s) mm; τ = L²/α; Fo = 1200/τ. 2. One term: T_back ≈ T_s + (T₀ − T_s)(4/π)e^{−(π²/4)Fo}; then `tile_temperature(np.array([L]), 1200, L, alpha, 300, 1300)`. The n = 1 term enters with (−1)¹ at x = L. 3. Solve T_back(t) = 450 K: invert the one-term formula for Fo, or `brentq` on `tile_temperature` in t.

## 2.
(a) Flux in minus flux out = rate of stored energy; Fourier's law for q; divide by dx and let dx → 0. (b) θ = XΘ; X″ + λ²X = 0; X(0) = 0 ⇒ X = sin λx; X′(L) = λcos λL = 0 ⇒ λL = (2n+1)π/2. With λ = nπ/L the slope at L is λcos nπ ≠ 0. (c) b_n = (2/L)∫₀ᴸ θ₀ sin λ_n x dx; evaluate; use cos λ_n L = 0 and λ_n L = (2n+1)π/2. (d) ∂²/∂x² gives −k², ∂²/∂y² gives +k² for sinh and −k² for sin; Laplace needs the sum zero; the heat equation needs X″/X = Θ′/(αΘ) negative for decay.

## 3.
(a) X′(0) = 0 kills the sine; X′(L) = 0 gives sin kL = 0, k = nπ/L including n = 0 (uniform pressure, a constant). f_n = nc/2L. (b) Insert L and c. (c) X(0) = 0, X′(L) = 0 (or vice versa): quarter-wave modes, f = (2n − 1)c/4L. (d) `string_mode_frequency(n, L, T, mu)` depends on T and μ only through √(T/μ) = c; closed–closed and fixed–fixed share f_n = nc/2L; the closed–open set is not n c/2L.

## 4.
(a) dt = rΔx²/α; n_steps = ⌈600/dt⌉; `heat_ftcs(..., left=("dirichlet", Ts), right=("neumann", 0.0))`; compare with `tile_temperature`. (b) Copy the interior update and the ghost-point update into a loop without the check; record `np.max(np.abs(u))`. Growth per step ≈ |1 − 4r|. (c) At κΔx = π the mode is e^{iπj} = (−1)^j. (d) n_steps ∝ 1/Δx² at fixed r.

## 5.
(a) `plucked_string_coefficients(L, a, h, 8) / b[0]`; look for sin(nπa/L) = 0. (b) `dalembert(x, t, c, L, f)` against `wave_series(x, t, L, c, b)`; the difference peaks at the kinks. (c) `wave_leapfrog` with dt = 0.8Δx/c, n_steps = round(0.01/dt); `string_energy` before and after. (d) Repeat with 121 points; the energy error is a kink-cell effect, so it scales with Δx.

## 6.
(a) Fixed hot face: heat enters at full rate from t = 0 (real surface heats up over tens of seconds) → model high. Insulated back: real structure absorbs heat → model high. Constant cold α: real hot α larger → model low. Magnitudes from the lesson's α sensitivity. (b) `tile_temperature` and `tile_thickness_for_limit` with 2α. (c) Two thermocouples at different depths: the ratio of excess temperatures is a function of Fo only, so a single α must fit both histories. (d) Exact for the model; the model's inputs are uncertain; state the uncertainty band.
