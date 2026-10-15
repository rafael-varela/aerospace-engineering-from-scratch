# P02.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. I_pivot = I_cg + m d²; L_eff = I_pivot/(m d); ω₀ = √(g/L_eff); T₀ = 2π/ω₀. 2. `swing_test_error(np.radians(15), case)` returns T_exact, the inferred I_cg and the relative error. 3. Lever = I_pivot/I_cg; compare the error with lever × θ₀²/8 (θ₀ in radians). 4. `np.degrees(max_amplitude_for_inertia_error(0.01, case))`.

## 2.
(a) y = (θ, θ̇); f = (θ̇, −2ζω₀θ̇ − ω₀² sin θ); equilibria where θ̇ = 0 and sin θ = 0; J = [[0, 1], [−ω₀² cos θ, −2ζω₀]]. (b) At θ = π, cos θ = −1; det(λI − J) = λ² + 2ζω₀λ − ω₀²; quadratic formula; the root with the + sign is positive since √(1 + ζ²) > ζ; doubling time ln 2/λ₊. (c) Energy at release ω₀²(1 − cos θ₀); θ̇ = ω₀√(2(cos θ − cos θ₀)); T/4 = ∫₀^θ₀ dθ/θ̇; half-angle identity; substitution sin(θ/2) = k sin φ gives dθ = 2k cos φ dφ/cos(θ/2) and √(cos θ − cos θ₀) = √2 k cos φ. (d) Insert k² = sin²(θ₀/2) ≈ θ₀²/4 in K's series and keep the first correction.

## 3.
(a) Equilibria at sin 2θ = 0: θ = 0, π/2, π. J = [[0, 1], [−3n²k cos 2θ, 0]]; det J = 3n²k cos 2θ: positive (center) at 0 and π, negative (saddle) at π/2. The stable one has the axis of least inertia (x, since k > 0 means I_x > I_z... check the sign convention: θ = 0 is the orientation in which the gravity-gradient torque restores) along the local vertical. (b) ω_lib = n√(3k); T_lib = 2π/ω_lib; ratio 1/√(3k); only the inertia *ratio* enters. (c) φ = 2θ: φ̈ + 3n²k sin φ = 0, a pendulum with ω₀² = 3n²k; amplitude φ₀ = 2θ₀; `period_exact(2*theta0, omega_lib) / (2*pi/omega_lib)`. (d) k < 0 means I_z > I_x; det J changes sign at every equilibrium: 0 and π become saddles, π/2 becomes the center; the satellite prefers the other axis vertical.

## 4.
(a) `integrate_rk4(lambda t, y: pendulum_rhs(t, y, w0, zeta), [radians(20), 0], 12*T0, 1e-3)`; peaks via `np.where((Y[1:-1,0] > Y[:-2,0]) & (Y[1:-1,0] >= Y[2:,0]))[0] + 1`; δ and ζ formulas as given. (b) `period_from_trajectory(t, Y[:, 0])`; compare with T₀, with `period_exact(radians(20), w0)` and with T₀/√(1 − ζ²); note the amplitude decays during the record so the mean period sits between the 20° value and T₀. (c) `inertia_from_period(T_measured, m, d)`; compare the size of (T/T₀ − 1) from amplitude (≈ θ_mean²/16) and from damping (≈ ζ²/2). (d) Measure T(20°) and T(10°); with T = T₀(1 + θ₀²/16), two equations give T₀ = (T₁ − T₂)/((θ₁² − θ₂²)/16) ... solve the pair for T₀ by elimination; then `inertia_from_period(T0_est, m, d)`.

## 5.
(a) ∂I/∂T = 2 m g d T/(4π²) = 2(I_pivot)/T; ∂I/∂d = m g T²/(4π²) − 2 m d = I_pivot/d − 2 m d; multiply each by its uncertainty and divide by I_cg. (b) `swing_test_error` at the three amplitudes. (c) Compare the three relative errors; the amplitude error is a bias (always positive), the others are random; a bias is removed by correction, not by averaging. (d) Quadratic damping, pivot compliance or out-of-plane swing; the measurement is the log decrement at several amplitudes, or the period at two suspension points.

## 6.
(a) J = [[0, 1], [−1, μ]]; tr = μ, det = 1, disc = μ² − 4: spiral for μ < 2, node for μ > 2, unstable in both cases. (b) `limit_cycle_amplitude(mu, [0.5, 0])` and `(mu, [3, 0])`. (c) W = μA² ∫₀^{2π} (1 − A² cos² t) sin² t dt = μA²(π − A²π/4) = μπA²(1 − A²/4). (d) A linearisation describes an O(|u|) neighbourhood of one point; the cycle's amplitude is set by a balance over a closed orbit far from the equilibrium, where the cubic term is not small.
