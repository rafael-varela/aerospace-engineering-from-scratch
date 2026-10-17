# P10.01  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. W = m g; q = ρV²/2; C_L = W/(qS); Re = Vc/ν; M = V/a. 2. Model: V_model = Re_full ν_room / c_model; M_model = V_model / a_room. 3. Compare M_model with the low-speed limit (about 0.3). 4. Choose which parameter governs the physics of interest (boundary layer vs compressibility).

## 2.
(a) Unit vectors of the stream (cos α, sin α) and of lift (−sin α, cos α) in chord axes; project (A, N). Show L² + D² = N² + A² by expanding. Inverse: rotate by −α. (b) Moment about P of F = (A, N) at O: (z_O − z_P) A − (x_O − x_P) N; rearrange. Check: N > 0 at x_O > x_P gives a nose-down contribution about P. (c) M_ref = M_ac + (x_ref − x_ac) N (no vertical offset); divide by qSc; differentiate with respect to c_n. c_n is the force the lever arm multiplies; c_l = c_n cos α − c_a sin α differs from it at second order in α. (d) ∂D/∂α = N cos α − A sin α = L; divide by qS; multiply by u_α in radians.

## 3.
1. α_b = α − 2°. 2. L, D from N_b, A_b rotated by α_b. 3. Moment transfer in balance axes from the balance centre to the quarter chord (Δx = −0.050 m, Δz = +0.030 m). 4. Divide by qS and qSc. 5. Repeat with α instead of α_b, and with Δz = 0; compare.

## 4.
1. `d = synthetic_balance_data(s)`; `r = reduce_balance(...)` with the MODEL geometry. 2. `np.polyfit(np.radians(d["alpha_deg"]), r["cl"], 1, cov=True)` gives a₀ and the intercept; α_L0 = −intercept/a₀. 3. `np.polyfit(r["cn"], r["cm"], 1, cov=True)`: x_ac/c = 0.25 − slope, c_m,ac = intercept; uncertainties from the square roots of the covariance diagonal. 4. Loop over seeds 1–20; standard deviation of x_ac/c (ddof = 1) vs mean fit uncertainty. 5. Number of σ = |x_ac − truth|/u.

## 5.
1. For each angle build the input vector [N, A, M, α, q] and the standard uncertainties from `NOISE`. 2. `propagate_uncertainty(point_coefficients-wrapper, x, u)`; shares = contributions/total. 3. `monte_carlo_uncertainty` for one angle. 4. Sensitivity to α: sqrt(α contribution)/u_α; the α uncertainty that equalises the shares is sqrt(A contribution)/sensitivity. 5. With α perfect, u(c_d) is the root-sum-square of the rest.

## 6.
For each quantity, decide whether it is mostly inviscid (shape: α_L0, x_ac, c_m,ac to first order, a₀ to first order) or viscous (c_d, c_l,max, and the small viscous corrections to a₀ and c_m). Rank by sensitivity to a factor of 6 in Re. Propose a measurement that probes the boundary-layer state (a transition trip, a Reynolds sweep, surface flow visualisation).
