# P03.18  Worked solutions

Each solution states: assumptions, units, the derivation, the result with its uncertainty, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. Every data set in this lesson is synthetic and is labelled so.

## 1. Grid convergence by hand ($s = 15$)

f₁ = 0.02015, f₂ = 0.02060, f₃ = 0.02255 (fine, medium, coarse), r = 2.

(a) d₂₁ = 4.5 × 10⁻⁴, d₃₂ = 1.95 × 10⁻³, ratio 4.333: **p = ln 4.333/ln 2 = 2.115**. f_ex = f₁ − d₂₁/(2^{2.115} − 1) = 0.02015 − 4.5 × 10⁻⁴/3.333 = **0.020015**. ε₂₁ = 4.5 × 10⁻⁴/0.02015 = 0.02233; **GCI_fine = 1.25 × 0.02233/3.333 = 0.837 %**.

(b) With p = 2: f_ex = 0.02015 − 4.5 × 10⁻⁴/3 = **0.020000** (the sequence was built as 0.0200 + C h², so this is exact), **GCI_fine = 1.25 × 0.02233/3 = 0.931 %**. ε₃₂ = 1.95 × 10⁻³/0.0206 = 0.09466; ratio = (0.09466/0.02233)/4 = **1.06**: within 6 % of 1, the grids are in the asymptotic range (the observed 2.115 against the formal 2 says the same; the coarse grid carries a small higher-order term).

(c) Factor = (0.837 %/0.1 %)^{1/2.115} = 8.37^{0.473} = **2.73** (2.89 with p = 2): a grid about 2.7× finer in each direction, 20× the cells in 3D.

**Check.** The code's `richardson(0.02015, 0.0206, 0.02255)` returns p = 2.1155, f_extrap = 0.020015, gci_fine = 0.008375 ✓, and with `p_formal=2.0` f_extrap = 0.0200, gci 0.009305, ratio 1.0597 ✓.

**Common mistakes.** Using the absolute difference instead of ε = d/f₁ in the GCI; forgetting that f₃ is the coarse grid in the NPARC convention; reporting a GCI without the order it was computed with.

## 2. Truncation error, exact discrete solution, Richardson, planted bug

(a) About t_{n+½}: u^{n+1} = u + (Δt/2)u_t + (Δt²/8)u_tt + (Δt³/48)u_ttt + …, uⁿ the same with −Δt/2. Difference: Δt u_t + (Δt³/24)u_ttt. Average θu^{n+1} + (1 − θ)uⁿ = u + (θ − ½)Δt u_t + (Δt²/8)u_tt + (θ − ½)(Δt³/48)·2u_ttt + …; the same for any quantity g(t) in place of u. Space: D₂u = u_xx + (Δx²/12)u_xxxx + O(Δx⁴). Put these into (1) with g = αu_xx − βu + q and subtract the PDE u_t = g at (x_i, t_{n+½}):

τ = (Δt²/24)u_ttt − (θ − ½)Δt g_t − (Δt²/8)g_tt − (Δx²/12)α u_xxxx(θ-averaged) + …,

so **τ = −(θ − ½)Δt g_t + O(Δt²) + O(Δx²)**: the first-order term vanishes only for θ = ½ ✓ (2).

(b) u = 2t x(L − x). Space: u_{i−1} − 2u_i + u_{i+1} = 2t[(x − Δx)(L − x + Δx) − 2x(L − x) + (x + Δx)(L − x − Δx)] = 2t(−2Δx²), so D₂u = −4t exactly. Time: (u^{n+1} − uⁿ)/Δt = 2x(L − x) exactly. Right side of (1): the quantity g = αD₂u + q = −4αt + 2x(L − x) + 4αt = 2x(L − x) is independent of t, so its θ-average is itself. Both sides equal 2x(L − x) identically, for every Δx, Δt and θ; the ends x = 0, L give u = 0 ✓. (For general θ this is still exact; for a source quadratic in t it would not be.)

(c) f₂ − f₁ = C(rh)ᵖ − Chᵖ = Chᵖ(rᵖ − 1); f₃ − f₂ = Chᵖrᵖ(rᵖ − 1). Ratio (f₃ − f₂)/(f₂ − f₁) = rᵖ, so **p = ln((f₃ − f₂)/(f₂ − f₁))/ln r**. From the first, Chᵖ = (f₂ − f₁)/(rᵖ − 1), and f = f₁ − Chᵖ = **f₁ + (f₁ − f₂)/(rᵖ − 1)** ✓ (5). GCI₃₂/GCI₂₁ = |ε₃₂|/|ε₂₁| = |(f₃ − f₂)/f₂|·|f₁/(f₂ − f₁)| = rᵖ f₁/f₂, so the asymptotic ratio GCI₃₂/(rᵖGCI₂₁) = f₁/f₂ ≈ 1 whatever the data: with the observed p it is not a test.

(d) The correct CN source is ½(q^{n+1} + qⁿ) = q(t_{n+½}) + (Δt²/8)q_tt + …; the buggy one is qⁿ = q(t_{n+½}) − (Δt/2)q_t + (Δt²/8)q_tt − …. The difference, (Δt/2)q_t, is a first-order truncation term. With Δt = cΔx the global error is C₁Δx² + C₂(c/2)Δx·(something of order q_t), and the observed order ln(e_coarse/e_fine)/ln 2 tends to 1 from above as Δx → 0 — exactly the 1.14, 1.08, 1.04, 1.02 of the lesson ✓.

**Common mistakes.** Expanding about t_n instead of t_{n+½} (the Crank–Nicolson second order then appears only after cancellation); in (b), forgetting to check the source's time dependence.

## 3. Verifying an ODE integrator by manufactured solution ($s = 15$: T(0) = 365 K)

(a) T_m = 250 + 115e^{−t/500}; T_m' = −0.23e^{−t/500} K/s. **q_m(t) = −0.23e^{−t/500} + 2 × 10⁻¹¹[(250 + 115e^{−t/500})⁴ − 250⁴]**. At t = 0: k(365⁴ − 250⁴) = 2 × 10⁻¹¹ × (1.7748 × 10¹⁰ − 3.906 × 10⁹) = 0.2769 K/s, so **q_m(0) = −0.230 + 0.277 = 0.0469 K/s**: the manufactured probe needs a small heater because the chosen decay (0.23 K/s) is slower than radiation alone (0.28 K/s).

(b) Heun to 1000 s:

| n | 10 | 20 | 40 | 80 | 160 |
|---|---|---|---|---|---|
| error [K] | 0.2777 | 0.0663 | 0.01619 | 0.00400 | 0.000994 |
| order | | 2.067 | 2.034 | 2.017 | 2.008 |

Second order, converging to 2 from above ✓.

(c) With k₂ evaluated at t: errors 0.651, 0.231, 0.0940, 0.0419, 0.0197 K, orders 1.49, 1.30, 1.17, 1.09 → **1**. The second stage uses q_m(t) where q_m(t + h) was needed, an error of hq_m' in k₂ and ½h²q_m' in the step: O(h²) per step, O(h) over 1/h steps, so the global order drops to one. The orders approach 1 slowly because the O(h²) part of the error is still large at coarse h (the 1.49 at n = 20); this is why an order study needs several refinements, not two.

(d) The source is additive; the code path for the nonlinear term −k(T⁴ − T_a⁴), the stage structure and the step update are the same whether q_m is present or not. A wrong k, a wrong power, a wrong stage weight or a wrong time argument all change the manufactured error's order or level. The only bug it could miss is one that depends on q being identically zero (e.g. a branch `if q is None`), which an explicit test of the unmodified equation against a converged reference covers.

**Check.** Dimensions: k [K⁻³ s⁻¹] × T⁴ [K⁴] = K/s ✓. Limiting case: at t → ∞, T_m → T_a and q_m → 0 ✓.

**Common mistakes.** Building q_m from the unmodified equation's residual with the wrong sign; concluding from the first order 1.49 that the bug "did not break the order".

## 4. Your own manufactured solution ($s = 15$: amplitude 35 K, β = 1.5 × 10⁻³ s⁻¹)

(a) At x = 0.1 m, t = 50 s (L = 0.5, τ = 200, α = 9.75 × 10⁻⁵): term 1: 35 sin(0.2π)(1 − e^{−0.25}) = 35 × 0.5878 × 0.2212 = 4.551 K; term 2: 15 × 0.1 × 0.4/0.25 × sin(π/2) = 2.400 K; **u_m = 6.951 K** ✓ (code 6.9506). q_m = u_t − αu_xx + βu: u_t = 35 × 0.5878 × e^{−0.25}/200 + 15 × 0.16 × (2π/200)cos(π/2) = 0.08011 + 0 = 0.08011; u_xx = −35(π/L)² sin(0.2π)(0.2212) − 2 × 15/L² × 1 = −179.7 − 120 = −299.7 K/m², −αu_xx = 0.02922; βu = 0.01043; **q_m = 0.1198 K/s** ✓ (code 0.11975).

(b) Crank–Nicolson, Δt = 300/nx:

| nx | 25 | 50 | 100 | 200 | 400 |
|---|---|---|---|---|---|
| error [K] | 1.379 × 10⁻² | 3.456 × 10⁻³ | 8.646 × 10⁻⁴ | 2.162 × 10⁻⁴ | 5.405 × 10⁻⁵ |
| order | | 1.997 | 1.999 | 2.000 | 2.000 |

`bug="source_time"`: errors 0.441, 0.217, 0.108, 0.0535 K, orders 1.023, 1.012, 1.006: first order, caught. (The bug's error is 5× larger than in the lesson's case because this u_m has a faster source variation, the sin(2πt/τ) term.)

(c) Heater with β = 1.5 × 10⁻³ s⁻¹, temperature at x = 0.2 m, t = 600 s: 30.6835 (25), 30.5727 (50), 30.5457 (100) K; **p = 2.037**, **f_ex = 30.5370 K**, **GCI_fine = 0.0356 %**; the 800-cell run gives 30.5369 K, so the extrapolation is within 10⁻⁴ K and the actual fine-grid error (0.0288 %) is inside the GCI ✓. Refinement for GCI 0.01 %: (0.0356/0.01)^{1/2.037} = **1.87**. The loss term lowers the heater temperature from 38.2 K (β = 0) to 30.5 K: 20 %, far more than any numerical error.

**Common mistakes.** Passing β to `manufactured_source` but not to `heat_theta` (the "manufactured" solution then is not a solution and the order collapses: a useful accident, since it shows what a model-form mismatch looks like); using the maximum over nodes as the functional.

## 5. Validation with your own data set ($s = 15$; SYNTHETIC data)

(a) Lossless model, 200 cells, Δt = 1 s, u_num = 0.012 K, u_val = 0.300 K: **3 of 24 points within 2u_val, reduced χ² = 62.3, max |E|/u_val = 16.9**, E at 600 s = (3.36, 5.07, 4.17, 2.21) K at x = 0.1, 0.2, 0.3, 0.4 m. All positive, growing with time, largest at the heater: a missing loss term.

(b) Hold-out fit on t ≤ 500 s: **β̂_ho = 7.60 × 10⁻⁴ s⁻¹** (10 model runs). Hold-out row E/u_val at 600 s: (2.15, 2.56, 0.71, −0.17), E = (0.65, 0.77, 0.21, −0.05) K: two of four slightly outside 2u_val, both positive, so the five-row fit slightly underestimates the loss and the model form is right but the parameter is not yet pinned down (the early rows carry little information on β because βu·t is small there). Full fit: **β̂ = 8.00 × 10⁻⁴ s⁻¹**, χ²'' from β̂ ± 2 × 10⁻⁵ gives **σ_β = 0.23 × 10⁻⁴ s⁻¹**; the truth 8.20 × 10⁻⁴ is **0.87σ** away ✓; reduced χ² = 1.08 at β̂ (0.99 at β_true), 22 of 24 within 2u_val.

(c) h = β̂ρcA/P = 8.00 × 10⁻⁴ × 2700 × 900 × 10⁻⁴/0.04 = **4.86 W/m²K** (truth 5.0; ±0.14 from σ_β). Heater temperature at 1200 s: **35.37 K with β̂, 35.25 K with β_true, 41.32 K with β = 0**. The design report gets 35.4 K ± 0.4 K, where 0.4 K is the propagated σ_β (the 1200 s prediction changes by about 0.17 K per 10⁻⁵ s⁻¹ of β), plus the sentence: "the loss coefficient was calibrated on the same synthetic data set (24 points, σ = 0.3 K, seed 15) that the model is compared with; the prediction at 1200 s extrapolates beyond the 600 s of data, and a second run at a different heater power is required before the calibrated value can be called validated." Numerics (0.012 K) do not appear in the uncertainty because they are 30× smaller than the calibration term.

**Check.** The lossless prediction is 6 K (17 %) high at 1200 s: the model-form error grows with time, so extrapolating the uncalibrated model would have been the costly mistake, not the grid.

**Common mistakes.** Quoting β̂ without σ_β; calling the calibrated model validated; forgetting that u_num and u_data add in quadrature (0.012 K changes u_val by 0.08 %).

## 6. A vendor says "validated"

(a) The plot shows a *single* comparison, on a *different* geometry, against a *correlation* (itself ±10–15 % on most duct data), on *one* grid. In NASA-STD-7009 terms it is weak validation evidence for a different intended use and no verification evidence at all. Ask first for: (1) a code-verification report with an observed order of accuracy on an exact or manufactured solution for the discretisation used (expect 2, accept ≥ 1.8); (2) a grid-convergence study on a geometry with a bend giving p, the extrapolated value and the GCI of the heat-transfer coefficient (expect a GCI below 3 % on the production grid); (3) the iterative-convergence criterion and a demonstration that tightening it by 10× changes the coefficient by less than 1 %.

(b) Requirements: (1) "The solver shall exhibit an observed spatial order of accuracy of 2.0 ± 0.2 for the energy equation on the manufactured solution of the V&V plan over three grids refined by r ≥ 1.3, verified by analysis of the discrete L2 error." (2) "The heat-transfer coefficient in the bend region, averaged over the wall patch defined in drawing X, shall have a grid convergence index not exceeding 3 % on the production grid, verified by a three-grid Richardson study with F_s = 1.25." (3) "The simulated coefficient shall agree with the channel rig measurement at Re = 20 000 within a comparison error not exceeding the validation uncertainty √(u_exp² + u_num²), verified by test and analysis, with u_exp documented by the rig's calibration report."

(c) f₁ = 390, f₂ = 395, f₃ = 412 W/m²K, r = 1.5: d₂₁ = 5, d₃₂ = 17, p = ln 3.4/ln 1.5 = **3.02** (above the formal order 2: the coarse grid is not yet asymptotic, so treat the result with caution). With observed p: f_ex = 390 − 5/(1.5^{3.02} − 1) = 390 − 5/2.40 = 387.9; GCI_fine = 1.25 × (5/390)/2.40 = **0.67 %**. With the formal p = 2: f_ex = 390 − 5/1.25 = 386.0, GCI_fine = 1.25 × 0.01282/1.25 = **1.28 %**, asymptotic ratio (17/395)/(5/390)/2.25 = 1.49 (not asymptotic). Either way the fine grid is adequate for ±10 % and, numerically, for ±2 %; but with a non-asymptotic sequence the honest numerical uncertainty is the conservative 1.3 %, which eats two thirds of a ±2 % budget before any model error is counted.

(d) "The simulated bend heat-transfer coefficient agrees with the single rig measurement to 4 %, within the 8 % uncertainty of that measurement combined with a 1.3 % numerical uncertainty (u_val = 8.1 %); this demonstrates consistency with the one available referent at the ±8 % level and does not demonstrate accuracy to the ±10 % the margin requires with any confidence beyond that, because a model error of up to 12 % is compatible with the observed agreement. A second test point at a different Reynolds number with a thermocouple uncertainty of 3 % or better is needed before the fan decision can rest on the simulation alone."

**Common mistakes.** Treating "agrees within 4 %" as "accurate to 4 %"; computing a GCI from a non-asymptotic sequence and quoting the smaller of the two numbers.

## Check-yourself answers

1. Verification: the implementation solves the chosen equations correctly (mathematics; exact solutions, orders, backward error). Validation: the chosen equations represent the real system for the intended use (physics; comparison with referent data and its uncertainty). The lossless heat model converges at second order with 10⁻⁴ K extrapolation error (verified) and misses the thermocouple data by 4.6 K (not validated).
2. Crank–Nicolson: O(Δx²) + O(Δt²), both shrink by 4 per halving → order 2. Backward Euler: O(Δx²) + O(Δt), the Δt term dominates → order 1.
3. For the polynomial, u_xxxx = 0 kills the spatial truncation term and the linear time dependence kills the temporal ones, so (1) is satisfied exactly; the sine solution has u_xxxx ≠ 0 and e^{−λt} has all time derivatives, so its errors are O(Δx²) + O(Δt²).
4. q_m = A sin(πx/L) e^{−t/τ}[−1/τ + α(π/L)² + β].
5. f_h = f + Chᵖ on h, rh, r²h: f₂ − f₁ = Chᵖ(rᵖ − 1), f₃ − f₂ = Chᵖrᵖ(rᵖ − 1); ratio rᵖ gives p; f = f₁ − (f₂ − f₁)/(rᵖ − 1).
6. Order ≈ 0: the error has a floor, the code converges to the solution of a different problem (wrong constant, wrong BC node, wrong sign in a small term). Order 1 instead of 2: a first-order term has crept in (a boundary difference, a source at the wrong time level, a one-sided stencil).
7. Tested: that the three values are consistent with a single power law and give a 0.03 % band. Not tested: that the grids are asymptotic, because the ratio equals f₁/f₂ by construction when p is observed; recompute with the formal order, or use a fourth grid.
8. Calibration adjusts parameters to minimise the disagreement with a data set; validation compares predictions with data *not used* for fitting, with stated uncertainties. The distinguishing test is a hold-out or a second condition (different heater power): a physical β̂ is unchanged, a fudge factor moves.
