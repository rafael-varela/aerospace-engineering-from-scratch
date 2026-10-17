# P08.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Ω = 2πN/60, U = ΩR. 2. τ = μU/h (`wall_shear_stress(h, U, 0, mu)`). 3. Torque = τ × (2πRL) × R; power = torque × Ω. 4. Re = ρUh/μ, Ta = Re² h/R; compare with 1700. 5. The power heats the oil; think about what a warmer oil does to μ and hence to the torque and the film.

## 2.
(a) Continuity with ∂u/∂x = 0 → ∂v/∂y = 0 → v = 0 (no-slip). y-momentum → ∂p/∂y = 0. x-momentum → μu″ = dp/dx; since the left side depends only on y and the right only on x, both equal a constant −G. Integrate twice, apply u(0) = 0, u(h) = U. Integrate u over the gap for q; differentiate for τ. Subtract the two stresses. (b) μ(1/r)(r w′)′ = −G; integrate once (r w′ = −Gr²/(2μ) + C₁; C₁ = 0 for finite w′ on the axis), again with w(R) = 0; integrate 2πr w dr for Q; τ_w = μ|w′(R)|; V = Q/(πR²); form f and simplify with Re. (c) u_t = U F′ η_t with η_t = −η/(2t); u_yy = U F″/(4νt); substitute and cancel; integrate F″/F′ = −2η twice; fix constants with F(0) = 1, F(∞) = 0. The similarity follows because the problem contains no length scale.

## 3.
1. Momentum equations with gx = g sin θ, gy = −g cos θ. 2. y-momentum: ∂p/∂y = −ρg cos θ → p = p_a + ρg cos θ (H − y). 3. x-momentum: μu″ = −ρg sin θ; u(0) = 0, u′(H) = 0. 4. q = ∫u dy → H³ formula; invert for H. 5. Surface speed u(H); Re = ρq/μ. 6. `ns_residuals(u, 0, p, rho, mu, gx=g sin θ, gy=-g cos θ)`.

## 4.
1. `actuator_leakage(D, h, L, dp, U, MU_HYD, RHO_HYD)` for U = 0, ±0.1; read Q, Re_gap, force, force_pressure, force_couette, t99. 2. Set q = Uh/2 + Gh³/(12μ) = 0 and solve for U. 3. Ratios ((h ± 1 µm)/h)³. 4. `crank_nicolson` across the gap from zeros to `startup_time(h, NU_HYD)` with source G/ρ; compare u at mid-gap with the steady value.

## 5.
(a) Differentiate; compare the two sides; depth from e^{−ky} = 0.01. (b) Exact initial data, `u_lo=lambda t: np.cos(om*t)`, `u_hi` the exact value at y = 10/k; errors at t = 3 periods; `observed_order`. (c) Same with zero initial data and u_hi = 0; loop over the number of periods. (d) k = √(πf/ν); compare 1/k with the clearance.

## 6.
(a) ΔT = Δp/(ρc_p); then the sign of dμ/dT. (b) Compare 2.5 with the other factors. (c) Re from Problem 4 or the lesson; 0.06 Re h against the land length. (d) Put the factors side by side: temperature (through μ), eccentricity, tolerance on h, entrance length.
