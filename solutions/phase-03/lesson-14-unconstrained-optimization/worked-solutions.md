# P03.14  Worked solutions

Worked for seed **s = 15** (κ = 25 in problem 1, a = 2.5 in problem 2, M = 0.715 in problem 3, W = 570 kN in problem 4, x₀ = (−1.3, 1) in problem 5). Every number below is reproduced by `solution.py` with the stated inputs; the "Check yourself" line at the end of each problem gives the values to compare against for a different s.

## 1. Steepest descent on a quadratic by hand

(a) Q = diag(1, 25), x₀ = (25, 1): g₀ = (25, 25), gᵀg = 1250, gᵀQg = 625 + 15 625 = 16 250, α₀ = 1250/16 250 = 1/13 = 0.07692 (the general formula 2/(1 + κ)). x₁ = (25, 1) − (25/13)(1, 1) = (300/13, −12/13) = (23.077, −0.9231) = (24/26) × (25, −1). The iterates stay on the line pair through (κ, ±1), shrinking by (κ − 1)/(κ + 1) = 0.9231 and flipping the sign of x₂ each step; the code's first four iterates are (25, 1), (23.077, −0.923), (21.302, 0.852), (19.663, −0.787).

(b) f is homogeneous of degree 2 and even in x₂, so f(x₁) = (24/26)² f(x₀): the ratio is **0.85207** per step, exactly the bound (1) for κ = 25. To reduce f by 10⁻¹⁰ needs k = ln(10⁻¹⁰)/ln(0.85207) = **143.8 → 144 iterations**. Check after 10 steps: 0.85207¹⁰ = 0.2017.

(c) `steepest_descent` with `exact_quadratic=Q`: the measured ratio is 0.852071 at every step (first 0.8520710, worst 0.8520710), f₁₀/f₀ = 0.20172, and the run stops after **333 iterations** at ‖g‖ ≤ 10⁻¹⁰, which is f ≈ 10⁻²⁰: ln(10⁻²⁰/325)/ln(0.85207) = 323, consistent with the extra 180-odd iterations beyond (b). `bfgs` with the exact line search finishes in **2 iterations** (n = 2) with H = diag(1, 0.04) = Q⁻¹ to rounding.

Check yourself: ratio ((κ − 1)/(κ + 1))², iterations ln(10⁻¹⁰)/ln(ratio) (≈ 144 for κ = 25, 93 for κ = 20, 219 for κ = 31); BFGS always 2 iterations, H = diag(1, 1/κ).

## 2. The BFGS update

(a) ∇f(x_{k+1}) = ∇f(x_k) + ∇²f(x_k)(x_{k+1} − x_k) + O(|s|²), so y ≈ Bs with B the Hessian; asking the inverse approximation to reproduce it exactly gives H_{k+1}y = s. For a quadratic the relation y = Qs is exact. The secant equation is n scalar equations; a symmetric n × n matrix has n(n + 1)/2 free entries, so for n ≥ 2 infinitely many H satisfy it. BFGS picks the one closest to H_k in the weighted norm ‖W^{1/2}(H − H_k)W^{1/2}‖_F with W = B (`@ocw-18-335j` Lecture 31).

(b) Right-multiplying (3) by y: (I − ρysᵀ)y = y − ρy(sᵀy) = y − y = 0, so the first term vanishes and ρs(sᵀy) = s. Symmetry: VH_kVᵀ is symmetric when H_k is, and ssᵀ is symmetric. Positive definiteness: for z ≠ 0, zᵀH_{k+1}z = (Vᵀz)ᵀH_k(Vᵀz) + ρ(sᵀz)² with Vᵀ = I − ρysᵀ. Both terms are ≥ 0 (H_k ≻ 0, ρ > 0). If the first vanishes then Vᵀz = 0, i.e. z = ρy(sᵀz); if also sᵀz = 0 then z = 0, contradiction; so the second term is strictly positive.

(c) φ'(α) = ∇f(x + αp)ᵀp and φ'(0) = ∇f(x)ᵀp < 0. The curvature condition gives φ'(α) ≥ −c₂|φ'(0)| = c₂φ'(0) > φ'(0) since 0 < c₂ < 1 and φ'(0) < 0. Then yᵀs = (∇f(x + αp) − ∇f(x))ᵀ(αp) = α(φ'(α) − φ'(0)) > 0.

(d) s = (1, 0), y = (2.5, 1), yᵀs = 2.5, ρ = 0.4. V = I − ρsyᵀ = I − 0.4[[2.5, 1], [0, 0]] = [[0, −0.4], [0, 1]]. VVᵀ = [[0.16, −0.4], [−0.4, 1]]; ρssᵀ = [[0.4, 0], [0, 0]]. **H₊ = [[0.56, −0.4], [−0.4, 1]]**. Check: H₊y = (0.56 × 2.5 − 0.4, −0.4 × 2.5 + 1) = (1, 0) = s. Trace 1.56, determinant 0.56 − 0.16 = 0.40, eigenvalues 0.3235 and 1.2365: positive. The one-variable update gives H = s/y = 1/2.5 = 0.40 = det H₊; H₊'s (1, 1) entry 0.56 is larger because the second component of y says the gradient also changed in the direction not travelled, and the update distributes the curvature accordingly.

Check yourself: with a = 1 + s/10, H₊ = [[1/a² + 1/a, −1/a], [−1/a, 1]], det = 1/a, H₊y = (1, 0).

## 3. Step size and scaling

(a) At x = (0.715, 1.0), f = 2.45614 kg/km and the Richardson reference gradient is (−1.6047, −1.0019) kg/km per unit of (M, h/10 km). Errors in ∂f/∂M (second column ∂f/∂(h/10 km)):

| h | central | forward |
|---|---|---|
| 10⁻¹ | 1.2 × 10⁻¹, 4.6 × 10⁻³ | — |
| 10⁻² | 2.4 × 10⁻³, 4.2 × 10⁻⁵ | 1.0 × 10⁻¹, 2.7 × 10⁻² |
| 10⁻³ | 2.5 × 10⁻⁵, 4.2 × 10⁻⁷ | — |
| 10⁻⁴ | 2.5 × 10⁻⁷, 4.2 × 10⁻⁹ | 1.1 × 10⁻³, 2.7 × 10⁻⁴ |
| 10⁻⁵ | 2.5 × 10⁻⁹, 5.7 × 10⁻¹¹ | — |
| 10⁻⁶ | **2.6 × 10⁻¹⁰**, 4.5 × 10⁻¹⁰ | 1.1 × 10⁻⁵, 2.7 × 10⁻⁶ |
| 10⁻⁷ | 7.0 × 10⁻¹⁰, 3.3 × 10⁻⁹ | — |
| 10⁻⁸ | 7.0 × 10⁻¹⁰, 1.9 × 10⁻⁸ | **1.6 × 10⁻⁷**, 6.3 × 10⁻⁸ |
| 10⁻⁹ | 2.0 × 10⁻⁷, 2.0 × 10⁻⁷ | — |
| 10⁻¹⁰ | 6.4 × 10⁻⁷, 8.6 × 10⁻⁷ | — |

Central: slope exactly +2 (a factor 100 per decade) down to h = 10⁻⁵, best at h ≈ 10⁻⁶ (2.6 × 10⁻¹⁰), then rising roughly as 1/h. Forward: slope +1 down to 10⁻⁶, best near 10⁻⁸ (1.6 × 10⁻⁷), three decades worse than the central best. The turning points are where C h² (C ≈ 24 for ∂/∂M) meets ε|f|/h ≈ 5 × 10⁻¹⁶/h: h ≈ (5 × 10⁻¹⁶/48)^{1/3} ≈ 2 × 10⁻⁶ for central, and where C'h ≈ 10 meets the same rounding term, h ≈ 7 × 10⁻⁹, for forward.

(b) In (M, h_m) the gradient is (−1.6047, −1.0019 × 10⁻⁴): the altitude component is 10⁴ times smaller and the (2, 2) Hessian entry 10⁸ times smaller. Condition numbers: **1.4 × 10⁹** in (M, h_m) against **19.8** in (M, h/10 km) at this point (9.9 at the optimum).

(c) Scaled: converged in **9 iterations**, f = 2.252861 at (0.7760, 1.1955). Unscaled: **not converged in 500 iterations** by the gradient test ‖g‖ ≤ 10⁻⁷, although it sits at (0.7760, 11 954.6 m) with f = 2.252861 to 13 digits. The unscaled gradient's altitude component is of order 10⁻⁴ × (relative error), so ‖g‖ is dominated by the Mach component and the stopping value 10⁻⁷ asks the Mach gradient to be smaller than the finite-difference gradient's own noise (10⁻⁹ to 10⁻⁷ relative to a function of order 2); the same test in scaled variables is reachable because the finite-difference noise is spread over two components of comparable size and the BFGS matrix is well-conditioned.

(d) Replace ‖g‖ ≤ tol by ‖Dg‖ ≤ tol with D = diag(1, 10⁴) (typical sizes of the variables), which is identical to rescaling; or stop on |f_k − f_{k−1}| ≤ 10⁻¹⁰|f_k|, which would have stopped the unscaled run at the same point. The cost: a relative-decrease test can stop on a plateau or during a short step that is not near a stationary point, and it certifies nothing about the gradient; the scaled-gradient test keeps the guarantee and should be preferred, which is why the design variables are scaled to order one from the start.

Check yourself: for s ≠ 15 the gradient values change in the third digit but the table's slopes, the best h (10⁻⁶ central, 10⁻⁸ forward) and the two condition numbers (within 10 %) do not.

## 4. Cruise-climb

(a) Predictions from §2.7: M* = 0.7760 and C_L* = 0.585 unchanged; h* = 11 955 + 6341.6 × ln(600/570) = 11 955 + 325 = **12 280 m**; fuel per km = 2.2529 × 570/600 = **2.1403 kg/km**.

(b) `optimal_cruise(dict(ac, W=570e3))`: M = 0.7760, h = 12 280 m, C_L = 0.585, L/D = 16.14, **2.1402 kg/km** in 8 iterations. Table:

| W [kN] | M* | h* [m] | C_L | fuel [kg/km] |
|---|---|---|---|---|
| 600 | 0.7760 | 11 955 | 0.585 | 2.2529 |
| 570 | 0.7760 | 12 280 | 0.585 | 2.1402 |
| 540 | 0.7760 | 12 623 | 0.585 | 2.0276 |
| 500 | 0.7760 | 13 111 | 0.585 | 1.8774 |

h* against ln W is a straight line of slope −6341 m (from 600 to 500 kN: 1156 m over ln 1.2 = 0.1823), the ISA scale height RT/g₀ = 6341.6 m; the local slope is −10.8 m/kN between 600 and 570 kN against the tangent value −10.6 m/kN at 600 kN.

(c) At 570 kN: M 0.78 at 11 km burns 2.2034 kg/km, **+2.95 %**; continuing at the 600 kN optimum point (0.7760, 11 955 m) burns 2.1445 kg/km, **+0.20 %**. Staying 325 m low costs a fifth of a per cent; the flight level, not the Mach, is what cruise-climb buys.

(d) Hessian at the new optimum (h = 10⁻³): [[81.59, 2.86], [2.86, 8.37]] kg/km per unit²; 1 % half-widths ΔM = √(2 × 0.01 × 2.1402/81.59) = **0.0229**, Δ(h/10 km) = √(2 × 0.01 × 2.1402/8.37) = 0.0715 → **715 m**. The altitude is the likelier binding constraint: 12 280 m (FL 403) is near the certified ceiling of a typical narrow-body and the optimum keeps climbing as fuel burns (13 111 m at 500 kN); the Mach window of ±0.023 sits comfortably below a typical M_MO of 0.82.

Check yourself: h* = 11 955 + 6341.6 ln(600/W) m with W = 600 − 2s kN; fuel = 2.2529 W/600; the Hessian scales with W so the half-widths do not change with s.

## 5. The three methods on Rosenbrock

(a) From x₀ = (−1.3, 1), f₀ = 52.9:

| method | iterations | f evaluations | ∇f evaluations | final f |
|---|---|---|---|---|
| steepest descent (tol 10⁻⁶) | 14 546 | 158 802 | 44 342 | 6.3 × 10⁻¹³ |
| damped Newton (tol 10⁻¹⁰) | 23 | 54 | 24 | 0 |
| BFGS (tol 10⁻¹⁰) | 25 | 63 | 76 | 3.3 × 10⁻²⁴ |

SciPy's BFGS with the analytic gradient takes 35 iterations and 43 function evaluations from the same start. Newton's 23 Hessian evaluations are not in the table; each costs a 2 × 2 solve.

(b) Last five error ratios: steepest descent 0.9992, 0.9992, 0.9992, 0.9992, 0.9992 (order estimates 1.00 throughout): **linear**, with a factor so close to 1 that 14 500 steps are needed. Newton 0.381, 0.247, 0.062, 0.0061, ~0 with errors 0.143, 0.0545, 0.0135, 8 × 10⁻⁴, 5 × 10⁻⁶, 10⁻¹⁰, 0 (orders 1.45, 1.99, 1.84, 2.07): **quadratic** once the full step is taken, the ratio being proportional to the error itself. BFGS 0.267, 0.036, 0.0098, 0.021, 0.012 (orders 2.51, 1.39, 0.83, 1.14): **superlinear**, the ratios tending to zero but not monotonically and without a fixed power. The order estimator divides log(e_{k+1}/e_k) by log(e_k/e_{k−1}); at the end both errors are within a few digits of the rounding floor 10⁻¹⁶, so a last error of 10⁻¹⁷ against 10⁻¹⁰ gives an order that depends on rounding alone (Newton's last estimate 1.63 in `main()` is of this kind).

(c) The exact inverse Hessian at (1, 1) is [[0.5, 1], [1, 2.005]]; the final BFGS H differs from it by up to **5.2 × 10⁻³** (3.4 × 10⁻⁵ from (−1.2, 1)). BFGS needs only that −Hg be a good direction and that the secant equations hold along the last few steps s_k; those steps span the valley direction well and the cross-valley direction poorly, and the curvature across the valley is never fully learned because the iterates never travel that way at the end. The Dennis–Moré condition for superlinear convergence requires H to be accurate only along the step directions, which is exactly what the secant updates enforce.

(d) Pure Newton (α = 1, no check) from (−1.3, 1) converges in **6 iterations**, but its path passes through f = **2197**, forty times the starting value 52.9: the first Newton step overshoots far up the far wall of the valley, where the Hessian happens to be positive definite and the next step recovers. It converges here by the structure of this particular function; on a function with an indefinite Hessian along that excursion it would have climbed away. The line search is what makes Newton's method a descent method; the positive-definiteness check is what makes its direction a descent direction.

Check yourself: Newton 20–25 iterations, BFGS 25–40, steepest descent > 10⁴ for any start x₀ = (−1 − s/50, 1); Newton's last ratios ∝ errors, BFGS's last ratios < 0.1 on average; pure Newton 6–8 iterations with a large excursion in f for every seed.

## 6. What the optimiser cannot tell you

(a) Established: the gradient vanishes to 5 × 10⁻¹⁰ at (0.7760, 11 955 m); the Hessian there is positive definite (eigenvalues 8.7 and 86.0), so the point is a strict local minimum; the minimum value 2.2529 kg/km and its curvature (hence the 1 % window) are known. Not established: that it is the global minimum (the single-basin evidence comes from three other starting points, not from the optimiser); that the point is reachable (ceiling, buffet margin at C_L = 0.585, M_MO, the discrete flight levels ATC offers); that the model is right (drag rise, TSFC law, no wind); how the answer moves with weight or temperature (a separate analysis, §2.7).

(b) Hessian estimate: ΔM = 0.80 − 0.7760 = 0.0240, Δf = ½ × 85.88 × 0.0240² = 0.0246 kg/km, **1.09 %**. Direct evaluation at (0.80, 11 955 m): 2.2813 kg/km, **1.26 %**; re-optimising the altitude for M = 0.80 (11 818 m) gives 2.2805, 1.22 %. Quote the direct value: the quadratic model misses the quartic growth of the wave-drag term and under-estimates the penalty by a seventh already at the edge of the 1 % window. For 3.1 % more speed the airline pays 1.2–1.3 % more fuel per km, which is why a cost index moves the cruise Mach up.

(c) With 0.5 % noise in f the "rounding" term in §2.6 becomes δ|f|/h with δ|f| ≈ 0.011 kg/km instead of 5 × 10⁻¹⁶: against the truncation term 24h² the best central step is h ≈ (0.011/48)^{1/3} ≈ 0.06 in the scaled variables, i.e. ΔM = 0.06 and Δh = 600 m, and the gradient error at that step is about 0.4 kg/km per unit, a quarter of the gradient at the starting point and larger than the gradient over the whole 1 % window. BFGS steps near the optimum would then have gradient changes y dominated by noise, yᵀs would be negative about half the time and the updates skipped, and the method would stall at the noise floor. Use Nelder–Mead (`@scipy-userguide` "Nelder-Mead Simplex algorithm") or fit a smooth quadratic surrogate to many evaluations and minimise that; or reduce the noise by averaging before differencing.

(d) Minimise f(M, h) = fuel per km subject to h ≤ h_ceiling(W) (thrust and buffet limits, a function of weight), M ≤ M_MO, C_L(M, h) ≤ C_L,buffet(M), optionally h ∈ {flight levels} (which makes the problem mixed-integer: solve the continuous problem and compare the two adjacent levels). SciPy: `minimize(method='SLSQP')` or `'trust-constr'` with inequality constraints, or `'L-BFGS-B'` if only bounds remain. Check in the output: which constraints are active, that the Lagrange multipliers of the active ones are non-negative (the objective would improve by crossing them), that the reduced gradient is zero to tolerance, and that the solution is the same from two starting points. For the 570 kN aircraft of problem 4 the ceiling constraint will be active and the optimum will sit on it with M re-optimised slightly higher (the one-dimensional optimum at fixed h).

Check yourself: the Hessian penalty at M = 0.80 is ½H₁₁(0.80 − M*)²/f* ≈ 1.1 % and the direct one 1.2–1.3 % for every seed (the aircraft is the given one); the noise-floor step h ≈ (δ|f|/2C)^{1/3} ≈ 0.06.
