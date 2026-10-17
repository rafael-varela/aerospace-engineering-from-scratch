# P03.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. Sizing an advection run (s = 15: a = 11.0 m/s)

**Assumptions.** Constant speed; upwind scheme; the modified-equation diffusivity applies (the slug is resolved: σ₀ = 2.5 cells).

**Derivation.** (a) Δx = 3/150 = **2.0 cm**; Δt = νΔx/a = 0.8 × 0.02/11 = **1.455 ms**; t = 2/11 = **0.1818 s**; steps = 0.1818/0.001455 = **125**.

(b) D = ½ × 11 × 0.02 × 0.2 = **0.0220 m²/s**; σ² = 0.05² + 2 × 0.022 × 0.1818 = 0.0025 + 0.0080, σ = **0.1025 m** (the slug doubles in width). At ν = 0.4: D = **0.0660 m²/s**, σ = **0.163 m**.

(c) Need 2Dt ≤ (1.1² − 1)σ₀² = 0.21 × 0.0025 = 5.25 × 10⁻⁴ m²; D ≤ 1.44 × 10⁻³ m²/s; Δx ≤ 2D/(a(1 − ν)) = 1.31 mm; **n ≥ 2290 cells** (and 125 × 2290/150 ≈ 1900 steps). A first-order scheme needs a grid about fifteen times finer than the one that merely resolves the slug before its own smearing becomes a 10 % effect.

**Check.** At ν = 1 the formula gives D = 0 and σ = σ₀, the exact shift ✓. The σ² growth at ν = 0.4 is three times that at ν = 0.8, as (1 − ν) predicts ✓.

**Common mistakes.** Using σ instead of σ² in the spreading law; forgetting that the step count also grows when the grid is refined at fixed ν.

## 2. The integral law, the fluxes, and the modified equation

(a) ∫_{t^m}^{t^{m+1}} ∫_{x_{i−1/2}}^{x_{i+1/2}} [u_t + (au)_x] dx dt = 0. The first term integrates in t to ∫[u(x, t^{m+1}) − u(x, t^m)] dx = Δx(ū_i^{m+1} − ū_i^m) by the definition of the cell average. The second integrates in x to ∫[au(x_{i+1/2}, t) − au(x_{i−1/2}, t)] dt = Δt(F_{i+1/2} − F_{i−1/2}) with F_{i+1/2} = (1/Δt)∫ au(x_{i+1/2}, t) dt. Dividing by Δx gives (1). The cell averages and the identity are exact; **only F**, the time-averaged interface flux, is approximated, because the scheme knows u only as cell averages at t^m.

(b) u^{m+1} = u + Δt u_t + ½Δt² u_tt + O(Δt³) = u − aΔt u_x + ½a²Δt² u_xx + …. With u_x ≈ (ū_{i+1} − ū_{i−1})/(2Δx) and u_xx ≈ (ū_{i+1} − 2ū_i + ū_{i−1})/Δx²: ū_i^{m+1} = ū_i − ½ν(ū_{i+1} − ū_{i−1}) + ½ν²(ū_{i+1} − 2ū_i + ū_{i−1}). Write the right side as ū_i − (Δt/Δx)[F_{i+1/2} − F_{i−1/2}] with F_{i+1/2} = ½a(ū_i + ū_{i+1}) − ½aν(ū_{i+1} − ū_i): the first bracket gives the central term, the second the ½ν² D₂ term ✓.

(c) Upwind: ū_i^{m+1} = ū_i − ν(ū_i − ū_{i−1}). Expand: left side u + Δt u_t + ½Δt² u_tt; right side u − ν(Δx u_x − ½Δx² u_xx + O(Δx³)). Divide by Δt (ν/Δt = a/Δx): u_t + ½Δt u_tt = −au_x + ½aΔx u_xx + O(Δx²). With u_tt = a²u_xx + O(Δx): u_t + au_x = ½aΔx u_xx − ½a²Δt u_xx = **½aΔx(1 − ν) u_xx**, which is (2). At ν = 1 the update reads ū_i^{m+1} = ū_{i−1}^m: the cell averages move one cell per step, and since the exact solution moves exactly aΔt = Δx, the averages are exact for ever.

(d) g = 1 − ν + νe^{−iθ}. |g|² = (1 − ν + ν cos θ)² + ν² sin²θ = (1 − ν)² + 2ν(1 − ν) cos θ + ν²(cos²θ + sin²θ) = 1 − 2ν(1 − ν)(1 − cos θ). Since 1 − cos θ ∈ [0, 2], |g|² ≤ 1 for all θ iff ν(1 − ν) ≥ 0, i.e. **0 ≤ ν ≤ 1**; and for ν(1 − ν) < 0 the mode θ = π has |g|² = 1 + 4|ν(1 − ν)| > 1.

**Common mistakes.** In (c), dropping the ½Δt² u_tt term and concluding D = ½aΔx (missing the −ν); in (d), forgetting that stability must hold for every θ.

## 3. A line whose speed changes (a₁ = 10 m/s, a₂ = 20 m/s)

(a) From u_t + (au)_x = 0 with a = a(x): q = au gives q_t = au_t = −a(au)_x = −aq_x, so **q_t + aq_x = 0**: q is constant along the characteristics dx/dt = a(x). A slug entering at speed a₁ with peak u₁ leaves at speed a₂ with q unchanged, so u₂ = u₁a₁/a₂ = **0.5u₁**: the peak halves. The front of the slug reaches the fast section before the tail, so the slug's length grows by a₂/a₁ = **2**. Under the advective form u itself is constant along characteristics: the peak is unchanged and the length still doubles.

(b) A mixture ratio is **intensive**: a parcel of fluid keeps its composition when it speeds up, so the advective form is right for it. The mass of excess oxidiser per metre of line is **extensive**: when the parcel stretches to twice its length the same mass is spread over twice the distance, so its density per metre halves, which is the conservation form. The link is mass flow: at fixed mass flow ṁ = ρAa, a ∝ 1/A, and (excess mass per metre) = ρA × (excess fraction) ∝ (fraction)/a, which is exactly q/a with q = a × (mass per metre)... the two are the same physics in different variables, and the product au of the conservation form is the flux of excess mass, constant along characteristics in a steady flow.

(c) Σ_i(ū_i^{m+1} − ū_i^m)Δx = −Δt Σ_i(F_{i+1/2} − F_{i−1/2}) = 0 on a periodic ring, whatever F is, so **Σū_iΔx is conserved**. Σ(aū)_iΔx is not: multiplying (1) by a_i before summing leaves a_i F_{i+1/2} − a_i F_{i−1/2}, which does not telescope when a_i ≠ a_{i+1}. In the narrowing line the physically conserved quantity is the excess mass, i.e. Σū_iΔx *when u is mass per metre*; if u is the mixture ratio, nothing of the form Σu is conserved and nothing should be.

(d) The periodic ring joins cell n − 1 to cell 0, so it requires a(L⁻) = a(0⁺); with a(x) rising from 10 to 20 m/s the flux leaving the last cell (20ū) is not the flux entering the first (10ū), and the "conservation" check will fail by exactly that mismatch each step. The real problem has an **inflow** condition (a prescribed u at x = 0, zero after the slug has entered) and an **outflow** condition at x = L (extrapolate ū from the last cell; information leaves and none returns).

**Check.** Steady state of the conservation form: (au)_x = 0, au = const, u ∝ 1/a: the slug halves in height exactly as the characteristic argument says ✓.

**Common mistakes.** Treating "conservation form" as always right; it is right for the conserved variable and wrong for a ratio.

## 4. The sensor reading (s = 15: slug width 0.35 m)

Exact front arrival at the sensor cell centre: n = 120 → sensor cell 90 at x = 1.810 m, (1.810 − 0.6 − 0.175)/12 = **86.25 ms**; n = 480 → cell 360 at 1.8025 m, **85.62 ms**.

(a) n = 120, ν = ½:

| scheme | 10 % [ms] | 50 % [ms] | peak | min |
|---|---|---|---|---|
| upwind | 75.90 | 86.13 | 0.891 | 0.000 |
| Lax–Friedrichs | 69.49 | 86.42 | 0.650 | 0.000 |
| Lax–Wendroff | 81.82 | 87.72 | 1.210 | −0.205 |
| minmod | 80.88 | 86.12 | 0.993 | 0.000 |

(b) n = 480:

| scheme | 10 % [ms] | 50 % [ms] | peak | min |
|---|---|---|---|---|
| upwind | 80.38 | 85.63 | 0.999 | 0.000 |
| Lax–Friedrichs | 76.64 | 85.39 | 0.935 | 0.000 |
| Lax–Wendroff | 83.91 | 86.29 | 1.233 | −0.236 |
| minmod | 83.49 | 85.62 | 1.000 | 0.000 |

Converged: the 50 % times (within 1 ms of exact for all), the upwind and minmod peaks (0.999, 1.000). Not converged: the 10 % times of the first-order schemes (still 5–9 ms early). Worse: the Lax–Wendroff overshoot and undershoot grew (1.210 → 1.233, −0.205 → −0.236).

(c) Threshold 0.9: on 120 cells **Lax–Wendroff and minmod** trip (upwind at 0.891 just misses, Lax–Friedrichs is far below); on 480 cells all but Lax–Friedrichs trip. The exact solution (peak 1.0) trips. The upwind model on 120 cells gives the wrong decision by 0.009.

(d) The modified equation says upwind spreads the front like diffusion with D = 0.06 m²/s (120 cells); in 0.086 s the smearing width √(2Dt) = 0.10 m, which the leading edge reaches 0.10/12 = 8.6 ms early while the symmetric 50 % point barely moves (observed: 10 % at 75.9 ms, 10.3 ms early; 50 % at 86.1 ms, 0.1 ms early). Lax–Wendroff has no diffusion but a phase speed below the true one for every resolved mode (0.93 at θ = π/4), so the sharp content of the front travels late and the 50 % crossing is 1.5 ms late, with the lagging modes appearing as the trailing wiggles.

**Common mistakes.** Comparing arrival times with the front time at x_sensor rather than at the sensor cell's centre (a 0.8 ms offset on 120 cells); reading the peak at a fixed time rather than from the time trace.

## 5. The Courant number as a knob (Gaussian slug, a = 12 m/s)

(a)

| ν | L1 errors (n = 240, 480, 960, 1920) | orders | e/e(0.9) at n = 1920 | predicted (1 − ν)/(1 − 0.9) |
|---|---|---|---|---|
| 0.25 | 7.71 × 10⁻², 4.49 × 10⁻², 2.46 × 10⁻², 1.29 × 10⁻² | 0.78, 0.87, 0.93 | 7.16 | 7.5 |
| 0.5 | 5.68 × 10⁻², 3.18 × 10⁻², 1.69 × 10⁻², 8.77 × 10⁻³ | 0.84, 0.91, 0.95 | 4.86 | 5.0 |
| 0.9 | 1.37 × 10⁻², 7.07 × 10⁻³, 3.58 × 10⁻³, 1.81 × 10⁻³ | 0.95, 0.98, 0.99 | 1 | 1 |

The ratios approach the (1 − ν) prediction from below as the grid is refined (5.6 → 7.2 for ν = 0.25), because on coarse grids the error saturates when the slug is heavily diffused.

(b) ν = 1: L1 error **3.9 × 10⁻¹⁶**, the exact shift. ν = 1.05 with the 10⁻⁶ sawtooth seed: max|u| = **1.039, 191, 2.62 × 10⁶** after 100, 200, 300 steps; |g(π)|ᵏ × 10⁻⁶ = 1.1ᵏ × 10⁻⁶ = 0.014, 190, 2.62 × 10⁶. After 100 steps the growing sawtooth (0.014) is still hidden under the slug (peak ≈ 1); by 200 steps it dominates exactly as predicted.

(c) Instability is a statement about the modes with |g| > 1, which for upwind at ν = 1.05 are the shortest ones (|g|² = 1 + 2 × 0.05 × 1.05 × (1 − cos θ), largest at θ = π). A smooth Gaussian contains almost none of them; only rounding error and the seed do, and 10⁻⁶ × 1.1ᵏ takes about 130 steps to reach 10⁻⁰. The run is unstable from the first step; it merely takes a while to show.

(d) From this lesson: (i) with a variable speed ν varies across the grid and the global step is set by the fastest cell, so most cells run below 1; (ii) the nonlinear fluxes of real problems change the local wave speed during the run, and a margin prevents a transient from crossing the limit. From `P03.05`: a multi-stage Runge–Kutta time integrator has a stability region of its own, and the combined limit for RK2 or RK3 with these spatial fluxes is below ν = 1.

**Common mistakes.** Declaring the ν = 1.05 run "stable" after 100 steps; expecting the ratio e(0.25)/e(0.9) to be 7.5 on the coarsest grid.

## 6. Choosing the scheme for the controller model

(a) ½aΔx(1 − ν) < 10⁻³ m²/s. At ν = ½: Δx < 2 × 10⁻³/(12 × 0.5) = **0.333 mm**, **7200 cells** on 2.4 m and 0.1 × 12/(0.5 × 0.333 × 10⁻³) = **7200 steps** for 0.1 s. At ν = 0.9: Δx < 1.667 mm, **1440 cells**, **800 steps**. Either way, a first-order scheme needs 10–60 times the grid of a limited second-order one to make its own smearing smaller than the physics.

(b) On 120 cells (lesson table, 0.3 m slug): bounded in [0, 1]: upwind, Lax–Friedrichs, minmod pass; Lax–Wendroff fails (−0.20). Arrival within 2 ms (50 % crossing against the sensor-centre time 86.25 ms + 0.83 ms): upwind (88.4 ms, +1.3 ms) and minmod (88.3) pass, Lax–Friedrichs (90.4) and Lax–Wendroff (89.9) fail. Peak within 5 %: only minmod (0.98) passes; upwind 0.83, Lax–Friedrichs 0.57, Lax–Wendroff 1.20 fail. Ranking: **minmod** (3/3), upwind (2/3), Lax–Wendroff (0/3), Lax–Friedrichs (0/3 but bounded). The boundedness requirement eliminates Lax–Wendroff on every grid, because its overshoot (1.20 → 1.23 on refinement) does not shrink.

(c) At the maximum of a smooth Gaussian the two one-sided differences have opposite signs, so minmod sets the slope to zero there every step: the reconstruction is locally first-order upwind and the peak erodes (0.875 on 120 cells, converging at order 1 at the peak while the sides converge at order 2). A square slug has a flat top: no extremum to erode, and the peak reads 0.981. For the controller the question is the plateau value of a valve-generated slug, which is flat-topped, so the square-slug number is the relevant one, and the 12 % Gaussian clipping is a warning for *smooth* signals (a slow mixture drift), not for this event.

(d) Clipping ū_i ← max(ū_i, 0) after each step adds mass wherever an undershoot occurred: Σū_iΔx is no longer conserved, and the added mass (0.20 × several cells × Δx per event) accumulates downstream as a spurious excess that the controller will see as a second, longer slug. A limiter achieves boundedness by changing the *flux* between cells, so that the undershoot never forms and every gram still goes from one cell to the next: it is the conservative, consistent form of "do not let it go negative". The Lax–Wendroff scheme with a minmod limiter on its antidiffusive part is, in fact, essentially the minmod scheme of this lesson.

**Common mistakes.** Judging schemes by their order alone; comparing the Gaussian-slug peak with the controller requirement when the physical slug is flat-topped.

## Check-yourself answers

1. ū_i^{m+1} = ū_i^m − (Δt/Δx)(F_{i+1/2} − F_{i−1/2}). Summing over all cells on a periodic ring, every F_{i+1/2} appears once with each sign and cancels, so Σū_iΔx is unchanged whatever F is.
2. ν = aΔt/Δx; |ν| ≤ 1. At ν = 2 the exact solution moves two cells per step while the stencil reaches one: the numerical domain of dependence cannot contain the physical one, so no explicit scheme on that stencil can be convergent, and the amplification factors confirm |g| > 1.
3. a > 0: F = aū_i; a < 0: F = aū_{i+1}. Update (a > 0): ū_i^{m+1} = (1 − ν)ū_i + νū_{i−1}, weights in [0, 1] summing to 1: the new value lies between the old ones, so no new extrema and TV non-increasing (monotone, TVD).
4. u_t + au_x = ½aΔx(1 − ν)u_xx. The variance grows by 2Dt = aΔx(1 − ν)t; zero at ν = 1.
5. Its error term is dispersive (u_xxx), not diffusive: short waves travel at the wrong speed and separate from the front as oscillations; the scheme is not monotone. Refinement moves the wiggles closer to the front but does not shrink their height (the sensor overshoot grows from 1.20 to 1.23 between 120 and 480 cells).
6. It limits the reconstructed slope to the smaller one-sided difference and to zero at extrema, so the scheme never creates new extrema (TVD) while keeping second order where the solution is smooth and monotone. At jumps and extrema it is first order, which is acceptable because no scheme does better than order ½–⅔ in L1 at a jump anyway and boundedness is worth more there.
7. At a jump the front is smeared over a width ∝ √(Dt) ∝ √Δx, so the L1 error (height × width) is O(Δx^{1/2}); on smooth data the error is the O(Δx) diffusion term acting on a bounded u_xx.
8. The grid. The modified equation gives D = ½ × 12 × 0.02 × 0.5 = 0.06 m²/s, sixty times the line's eddy diffusivity (10⁻³), so the smearing that lowers the peak to 0.83 is numerical; the exact solution's peak is 1.0.
