# P06.13  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. I = bt³/12; `cantilever_stiffness`; m_beam = ρbtL; `effective_mass`. 2. `natural_frequency` with m and with m_tip. 3. z_static = g/ω_n²; `bracket_root_stress`; n_yield = σ_y/(σ_1g × 1/(2ζ)).

## 2.
(a) T = ½ẇ_tip² ∫₀ᴸ ρA φ(x)² dx with φ = (3x²L − x³)/(2L³); substitute u = x/L; the integral is 33/140 L. (b) z = Z e^{iωt}, y = Y e^{iωt}: (k − mω² + icω)Z = mω²Y; X = Y + Z; take moduli; at r = √2 compare numerator and denominator. (c) Steady part Z = −A/(2iζω_n²) at r = 1; homogeneous part with the opposite initial values; for small ζ, ω_d ≈ ω_n, combine. 99 %: e^{−ζω_n t} = 0.01. (d) Half cycle: x = ±d + (x_start ∓ d) cos ω_n t; next turning point 2(±d) − x_start; energy per cycle: 4FX for friction, πcωX² for viscous.

## 3.
(a) Write |Z|ω_n²/|A| as a function of r; solve |gain − 1| = 0.05 numerically for each ζ (a fine grid or `brentq`); f_n = 2000 Hz/r_max; phase arctan(2ζr/(1 − r²)). (b) Solve |r²/D − 1| ≤ 0.05 for all r above r_min. (c) Compare gain flatness and phase linearity.

## 4.
(a) `synthetic_tap_test(s, ...)`, `find_peaks_refined`, `log_decrement` for each n, `zeta_from_log_decrement`; mean interval of 10 peaks. (b) Steady formulas of 2.4 and (9) with your ζ. (c) `simulate_base_sine(t, 2g, ω_n, ω_n, ζ)` for 1 s; max |z|; ratio. (d) Stress ∝ 1/ζ: evaluate at ζ ± 0.0004.

## 5.
(a) Add −2ζω_n ẋ to the right-hand side; keep the event logic. (b) `np.polyfit` of peaks against cycle number (linear) and of log peaks (exponential); residuals. (c) `np.polyfit(x_n, x_{n+1}, 1)`; ζ from δ = −ln q; derive by two half cycles, each multiplying the distance from the current centre by √q. (d) The first turning point with |x| ≤ d.

## 6.
(a) Required input = 2 g × 30/50; discuss the uncertainty in f_n and ζ. (b) Friction (12) versus microslip and joint nonlinearity. (c) Geometry and loading effects not in a tip-loaded beam. (d) Compare with the continuous-beam frequency equation with a tip mass.
