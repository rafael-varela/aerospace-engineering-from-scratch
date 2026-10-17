# P04.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. L = 0.5 + s/50. 2. `pendulum_period_small(L)`, `pendulum_period_exact(L, radians(10))`, `(…, radians(20))`; excess = T/T₀ − 1. 3. u = √2 × 0.5/(√3 × 30); per period u/N for each N. 4. N = u/(0.1 × T₀ × excess₁₀); round up; say how much the amplitude falls in that many swings for a lightly damped bob.

## 2.
(a) dE/dt = ẋ(mẍ + kx) = 0; at a turning point ẋ = 0, x = A. (b) x(t) = Ae^{−ζω₀t}cos(ω_d t − φ); peaks at t and t + 2π/ω_d; ratio = e^{ζω₀ 2π/ω_d}; δ = 2πζ/√(1 − ζ²); square and solve for ζ. (c) X² ∝ 1/((ω₀² − ω²)² + 4ζ²ω₀²ω²); at ω = ω₀(1 + ε): ω₀² − ω² ≈ −2εω₀², so X² ∝ 1/(4ω₀⁴(ε² + ζ²)); half power where ε² = ζ², ε = ±ζ. (d) Energy gives θ̇ = √(2g/L)√(cos θ − cos θ₀); dt = dθ/θ̇; integrate 0 → θ₀ and multiply by 4. Write cos θ − cos θ₀ = 2(k² − sin²(θ/2)); substitute sin(θ/2) = k sin φ, dθ = 2k cos φ dφ/√(1 − k² sin² φ). Expand (1 − k² sin² φ)^{−1/2} ≈ 1 + ½k² sin² φ; ∫ sin² φ = π/4; k² ≈ θ₀²/4. (e) F_y = T_s(sin α(x + dx) − sin α(x)), sin α ≈ tan α = u_x for small slopes; F_y = T_s u_xx dx = μ dx u_tt.

## 3.
(a) ζ, ω₀ = 2πf₀; Q = 1/2ζ; `resonance_peak_frequency`; `steady_state_amplitude` at the peak over the value at ω = 0. (b) Fine grid of ω around ω₀; keep where the ratio > 10; convert Hz to rpm (× 60); compare with `half_power_bandwidth`. (c) e^{−ζω₀t} = 0.01 → t = ln 100/(ζω₀); cycles = t f₀. (d) Ratio at 2ω₀ and ω₀/2; rule: keep the wheel's speed range, and its harmonics, outside the ×10 band with margin.

## 4.
(a) `string_wave_speed`, `standing_wave_frequencies`. (b) `wave_equation_fd(c, L, 200, 2L/c, C, pulse)`; compare every time level with `dalembert`; `wave_energy` spread = (max − min)/mean. (c) `wave_equation_fd(c, L, nx, 1.3/f₃, 0.5, mode 3)`; error at the last level against `standing_wave`; slope of log error against log Δx. (d) Same run at C > 1; max |u|; the highest grid mode's amplification factor exceeds 1 per step and the number of steps per round trip is about 2 × nx/C.

## 5.
(a) For each row: `period_from_timings(row, 10, u_single=0.41/30)`; `g_from_pendulum(L, T, 0.002, u_T)`. (b) Multiply each g by `amplitude_correction(0.15)**2`; weighted mean Σg/u² / Σ1/u², uncertainty 1/√Σ(1/u²). (c) `fit_power_law(lengths, T)`; g = 4π²/A²; times the correction. (d) E_n = |ḡ − 9.81|/u_ḡ (the generating value has no uncertainty); repeat with the uncorrected mean.

## 6.
Identify in each case the term that breaks: drag force comparable to weight (light bob), moving support (curtain rail), amplitude-dependent damping (nonlinear joint friction), numerical dispersion (Courant < 1 with sharp fronts), aerodynamic coupling (damping falling toward flutter). Give the observable and the remedy.
