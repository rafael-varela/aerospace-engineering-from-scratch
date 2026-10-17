# P10.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. AR; CL_α = 2π/(1 + 2/AR); CL at 4°; α_i = CL/(πAR); CDi = CL α_i. 2. `elliptic_wing(10, 100/AR)`, `solve_monoplane(w, 4°, 10)`, `wing_coefficients`. 3. Span doubles at fixed S: AR × 4.

## 2.
(a) dΓ/dη = (dΓ/dφ)/(dη/dφ), substitute (3), use Glauert's integral for each n; put w into (2), divide by 2bV∞, multiply by sin θ. (b) L = ρV∞(b/2)∫Γ sin θ dθ; D = ρ∫Γw dy, products of sines. (c) μ = μ0 sin θ; try A1 alone. (d) L′ = −ρV∞∫Γ y dy with the sign convention; only n = 2 survives.

## 3.
1. `rectangular_wing(10, 1.2, alpha_L0=lambda y: np.where(abs(y) < 3, -0.609 * 10°, 0))`; solve at α = 0 with three N. 2. 2πτδ × 0.6; and × CL_α/(2π) with CL_α from the clean wing. 3. `span_loading` at y = 2.9 and 3.1 m (N = 160).

## 4.
1. `tapered_wing(10, 12, s/40)`; `convergence_study(w, 5°, [10, 20, 40, 80])`. 2. `wing_coefficients` at N = 80; compare λ with 0.35. 3. Linear twist ∝ |y| has a kink at the root.

## 5.
1. `tapered_wing(10, 100/9, 0.4, twist=linear_twist(10, −2°), a0=k 2π)`; `lift_slope_and_tau`; convert to per degree. 2. Interpolate k for 0.085–0.086. 3. Section slope at the test Reynolds number; fuselage or tunnel corrections; lifting-line's missing chordwise effect; reading uncertainty.

## 6.
(a) AR 2, sweep, chordwise loading. (b) The lattice resolves the tips. (c) The integral equation smooths: Γ stays continuous, only its derivative jumps.
