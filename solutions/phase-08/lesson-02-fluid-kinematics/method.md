# P08.02  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) G = [[U₀/L, 0], [0, −U₀/L]]; tr G, G₂₁ − G₁₂, S = (G + Gᵀ)/2, eigenvalues of S. (b) Du/Dt = (u·∇)u: u ∂u/∂x and v ∂v/∂y; `particle_acceleration(vel, 1.0, 0.2, 0.0)`. (c) Separate variables: dx/(1 + x/L) = U₀ dt; t = (L/U₀) ln 2. Differentiate y(1 + x/L) along the path with the chain rule; steady flow, so pathlines are streamlines.

## 2.
(a) Integrate dx/dt = U, dy/dt = V cos ωt from t₀; the streamline at frozen t has slope V cos ωt/U; the streakline is the pathline with t₀ = τ evaluated at t, then eliminate τ = t − (x − x_s)/U; slope dy/dx at x = x_s. (b) Uniqueness of solutions of an autonomous ODE with a smooth right-hand side. (c) d/dt f(x(t), y(t), t) = f_t + f_x ẋ + f_y ẏ with ẋ = u, ẏ = v. (d) Edges δx, δy of a rectangle: d(δx)/dt = (∂u/∂x)δx, d(δy)/dt = (∂v/∂y)δy, so d(δxδy)/dt = (∇·u)δxδy; the rotation of the two edges is ∂v/∂x and −∂u/∂y; their average is ω/2.

## 3.
(a) Integrate in time: x = (W/Ω) sin Ωt, y = (W/Ω)(1 − cos Ωt); streamline at time t is a straight line along (cos Ωt, sin Ωt). (b) Particle released at τ: x = (W/Ω)(sin Ωt − sin Ωτ), y = (W/Ω)(cos Ωτ − cos Ωt); compute (x − x_c)² + (y − y_c)². (c) `pathline`, `streakline`, `streamline`; distances of the points from the circles. (d) Compare the streakline tangent at the source with the wind direction at t.

## 4.
(a) V = U tan 5°, ω = 2πf, wavelength 2πU/ω, amplitude V/ω. (b) Loop over L; for each, max of `smoke_error` over 121 phases in one period. (c) `brentq` on worst-error(L) − 0.05; Taylor-expand sin ω(t − x/U) to second order. (d) Compare k for 0.3 m with your 5 % limit.

## 5.
(a) Build the loop, 200 RK4 steps; `polygon_area`, perimeter as the sum of segment lengths, `np.ptp(x)`/0.2; compare with e². (b) Same with u = ax, v = 0; `decompose_gradient(velocity_gradient(...))`. (c) Area preservation is a property of each point following the exact flow map; the error is in resolving the loop with straight segments.

## 6.
For each method, decide: instant or exposure; released from a point or everywhere; steady or unsteady flow. Then name the line.
