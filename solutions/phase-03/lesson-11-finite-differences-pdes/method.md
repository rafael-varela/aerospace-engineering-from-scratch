# P03.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. α = k/(ρc_p); τ = L²/α; count = 10 s/τ. 2. Δx = L/n; Δt_max = Δx²/(2α); steps = 10 s/Δt_max; doubling n quarters Δt_max and quadruples the steps. 3. r = αΔt/Δx²; steps = 10 s/Δt; Crank–Nicolson for accuracy, after checking g(π) = (1 − 2r)/(1 + 2r) against the initial condition's grid-scale content.

## 2.
(a) Write u(x ± Δx) to fourth order and u(t + Δt) to second order; subtract the scheme from the PDE; use u_tt = α²u_xxxx; collect the Δx² and Δt terms into αΔx²(1/12 − r/2)u_xxxx. For BTCS the time expansion about t^{m+1} flips the sign of the Δt term, giving (1/12 + r/2). (b) Substitute gᵐe^{ijθ}; divide by gᵐe^{ijθ}; solve for g. FTCS: |1 − 4rs| ≤ 1 for s ∈ [0, 1] iff r ≤ ½. BTCS: denominator ≥ 1. Crank–Nicolson: |1 − z/2| ≤ |1 + z/2| for z ≥ 0. (c) Series: e^{−z} = 1 − z + z²/2 − …; 1 − z misses the +z²/2 (too small); 1/(1 + z) = 1 − z + z² − … (too large); (1 − z/2)/(1 + z/2) = 1 − z + z²/2 − z³/4 + … (matches to z², first error z³/12 relative to e^{−z}).

## 3.
(a) α = k/(ρc_p), τ = L²/α; compare 3600 s with τ. (b) Replace u_{n+1} by u_{n−1} in the FTCS stencil at j = n; write the update matrix rows; apply Gershgorin to the boundary row (centre 1 − 2r, radius 2r). (c) Copy `step_ftcs`, add the boundary row, loop with Δt = rΔx²/α; read the values at the mid-plane node and the last node for three grids. (d) Evaluate the series with 2000 terms at the two points and compare.

## 4.
(a) Steady mid-plane = (T_hot + T_cold)/2 > T_limit; call `time_to_limit_exact`. (b) Three calls to `time_to_limit`; subtract the exact time. (c) `cold_face_flux_exact` at the exact time; steady flux k(T_hot − T_cold)/L; divide. (d) First order: the error halves; run BTCS at Δt = 0.05 s.

## 5.
(a) Two calls to `temporal_convergence`; `observed_order` on each. (b) Δt* = (t_f/5)·e_CN/e_BTCS; steps ratio = (t_f/5)/Δt*; pick the nearest Δt that divides t_f and run. (c) Two calls to `spatial_convergence`; ratio of the last entries; orders from `observed_order`; the Δx² term of (4) vanishes at r = 1/6 leaving Δx⁴. (d) Constant α, uniform grid, u_tt = α²u_xxxx; temperature-dependent k breaks the first and the third.

## 6.
(a) r from Δt and Δx; g(π); N = ln(0.01)/ln|g(π)|; compare with 5 steps. (b) `march` with "cn" and "btcs" for 5 steps of 1 s; max |u − series| and max |second difference|. (c) Two BTCS steps of 0.5 s then four Crank–Nicolson steps of 1 s; same two measures. (d) Mean k over 300–700 K; α scales with k; the time to a threshold scales as 1/α; compare the shift with the numerical errors of Problems 4–6.
