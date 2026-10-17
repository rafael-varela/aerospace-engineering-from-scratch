# P06.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. M₀ = [[(m₁ + m₂)l₁², m₂l₁l₂], [m₂l₁l₂, m₂l₂²]], K₀ = diag((m₁ + m₂)gl₁, m₂gl₂). 2. a = m₁l₁l₂, b = −(m₁ + m₂)g(l₁ + l₂), c = (m₁ + m₂)g²; ω² from the quadratic; periods 2π/ω. 3. Shape from the first row: θ₂/θ₁ = (K₁₁ − ω²M₁₁)/(ω²M₁₂). 4. Single pendulum 2π√((l₁ + l₂)/g); hook estimate as given; percentage differences. 5. `normal_modes(DoublePendulum(m1, m2, l1, l2))`.

## 2.
(a) Differentiate (1); square and add; use cos θ₁cos θ₂ + sin θ₁sin θ₂ = cos Δ. Form ∂L/∂θ̇ᵢ, differentiate in time, subtract ∂L/∂θᵢ; collect θ̈ terms into M, the rest into f. (b) Multiply each Lagrange equation by q̇ᵢ, sum, recognise dL/dt with ∂L/∂t = 0; then Euler's theorem for the quadratic form T. (c) Substitute (6) into (7); dot the m₂ equation with n₂ (eliminates T₂); add the two equations and dot with n₁ (eliminates T₂ and T₁); multiply by l₂ and l₁ respectively; compare with (5) using nᵢ · nⱼ = cos(θᵢ − θⱼ) and nᵢ · eⱼ = sin(θᵢ − θⱼ) with the right sign. (d) Expand the determinant, divide by m₂l₁l₂, take the limits.

## 3.
(a) x = r sin θ, y = −r cos θ; T = ½m(ṙ² + r²θ̇²); V = mgy + ½k(r − l₀)²; `euler_lagrange` for [r, θ]. (b) Set rates and θ to zero: k(r₀ − l₀) = mg. Linearise each equation about (r₀, 0). (c) Solve k r₀ = 4mg with r₀ = l₀ + mg/k for k; integrate with `solve_ivp` (DOP853, rtol 10⁻¹⁰); compute the two energies. (d) Look for the time at which the swing energy peaks and what the bounce energy does then.

## 4.
(a) Loop over 1000 states from `default_rng(s)`; compare. (b) Solve Φc = q(0) with Φ the shape matrix of `normal_modes`; simulate 120 s; FFT (a Hann window); find the two peaks; compare with Problem 1. If the fast peak is off, release a pure mode 2 at amplitudes 0.1°, 1° and c₂ and measure its frequency. (c) Call `newton_euler` along the run; divide T₁ by (m₁ + m₂)g and T₂ by m₂g; compare with a mode-1 release of the same θ₁. (d) Two runs at different tolerances; energy error relative to the swing energy E − E(0, 0, 0, 0).

## 5.
(a) r₂ = (l₁ sin θ₁ + l₂ sin θ₂, −l₁ cos θ₁ − l₂ cos θ₂); differentiate with respect to each θ; dot with (−D, 0). (b) With all rates zero, (2) becomes ∂V/∂θᵢ = Qᵢ; solve for tan θᵢ. (c) Mass 2: T₂e₂ balances m₂g and D. Mass 1: T₁e₁ balances m₁g and T₂e₂. (d) Numbers; T₁ = √(((m₁ + m₂)g)² + D²); check that f(q, 0) + Q = 0 at the angles.

## 6.
(a) Required δ₀ = Δ e^{−λt} with Δ ≈ 10⁻⁴ rad; compare with double-precision round-off (≈ 10⁻¹⁶). (b) Conserved quantities and statistical properties. (c) `lyapunov_estimate` for the slung load at 10° over increasing times; regular motion gives an estimate that falls like (ln t)/t. (d) The hook mode depends on the hook mass and the cable tensions and on what the model holds fixed; list what in a real aircraft changes those.
