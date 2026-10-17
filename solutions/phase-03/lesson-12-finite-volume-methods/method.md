# P03.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Δx = L/n; Δt = νΔx/a; t = 2 m/a; steps = t/Δt. 2. D = ½aΔx(1 − ν); σ² = σ₀² + 2Dt; repeat with ν = 0.4. 3. 2Dt ≤ 0.21σ₀² gives Δx ≤ 0.21σ₀²/(a(1 − ν)t), then n = L/Δx.

## 2.
(a) ∫∫(u_t + (au)_x) dx dt = 0 over the cell and the step; the first term gives Δx(ū^{m+1} − ū^m), the second Δt(F_{i+1/2} − F_{i−1/2}) with F the time average of au at the interface; divide by Δx. Only F is approximated. (b) u^{m+1} = u − aΔt u_x + ½a²Δt² u_xx; centred differences; group as −(Δt/Δx)[F_{i+1/2} − F_{i−1/2}]. (c) Expand ū_{i−1} = ū − Δx u_x + ½Δx² u_xx and ū^{m+1} = ū + Δt u_t + ½Δt² u_tt; substitute u_tt = a²u_xx; solve for u_t + au_x. At ν = 1 the update is ū_i^{m+1} = ū_{i−1}^m. (d) g = 1 − ν + νe^{−iθ}; |g|² = (1 − ν)² + ν² + 2ν(1 − ν)cos θ = 1 − 2ν(1 − ν)(1 − cos θ); the bracket is non-negative iff ν(1 − ν) ≥ 0.

## 3.
(a) q_t = au_t = −a(au)_x = −aq_x. Along dx/dt = a, q is constant: u = q/a falls as a rises; the slug's length grows in proportion to a (front faster than tail). Under the advective form u is constant along characteristics and only the length changes. (b) Intensive → advective form; extensive (mass per metre) → conservation form; mass per metre = ρA × mixture excess, and a ∝ 1/A at fixed mass flow, so the two forms describe the same flow in different variables. (c) Sum (1) over cells: the fluxes telescope for any F; Σ(aū) does not telescope because a differs between cells. The total mass of excess species is conserved; the integral of the ratio is not. (d) A periodic wrap needs a(0) = a(L); otherwise the flux into cell 0 is not the flux out of cell n − 1. The real line needs an inflow value and an outflow extrapolation.

## 4.
(a) Build u0 from `cell_averages(square_pulse(…))`; i_s = int(x_sensor/dx); `sensor_history` for 0.2/dt steps; `arrival_time` at 0.1 and 0.5; max and min of the trace. (b) Same with n = 480. (c) Compare peaks with 0.9. (d) The 10 % point moves ahead by about the smearing width √(2Dt); Lax–Wendroff's phase speed ratio is below 1 for resolved modes, so its front lags.

## 5.
(a) Three calls to `convergence_study`; `observed_order`; predicted ratio (1 − ν₁)/(1 − ν₂). (b) `advect` at Δt = Δx/a; copy `fv_step` without the check; seed the sawtooth; record max|u| at 100, 200, 300 steps; compare with 1.1ᵏ × 10⁻⁶. (c) Only the modes with |g| > 1 grow: the shortest ones, which the Gaussian barely contains. (d) Variable speed makes ν vary across the grid; nonlinear fluxes change the wave speed; a Runge–Kutta integrator's stability region is smaller than the single-step CFL bound.

## 6.
(a) Δx = 2D_phys/(a(1 − ν)); n = L/Δx; steps = 0.1 a/(νΔx). (b) Boundedness, arrival, peak against the sensor table; Lax–Wendroff fails boundedness on every grid. (c) The limiter zeroes the slope at the extremum of a smooth peak every step, so a smooth maximum is eroded; a flat-topped slug has no extremum to erode. (d) Clipping adds mass each step; a limiter changes the flux between cells so that mass moves rather than appears.
