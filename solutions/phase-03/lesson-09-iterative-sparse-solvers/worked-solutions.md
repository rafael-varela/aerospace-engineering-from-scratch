# P03.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Storage and the flop budget ($s = 15$: n = 250, N = 62 500)

(a) nnz = 5n² − 4n = **311 500**. Dense: 8N² = **3.1 × 10¹⁰ bytes (31 GB)**. CSR: 12 × 311 500 + 4 × 62 501 = **4.0 × 10⁶ bytes (4 MB)**, a factor 7 800 less.

(b) Matvec: 2 nnz = **6.2 × 10⁵ flops**. CG iteration: 2 nnz + 10N ≈ **1.25 × 10⁶ flops** (the vector operations cost as much as the matvec when there are only five nonzeros per row).

(c) κ ≈ (2 × 251/π)² = 2.55 × 10⁴, √κ = 159.8. CG: m = ln(2 × 10⁶)/ln((√κ + 1)/(√κ − 1)) = **1 159 iterations**, 1.45 × 10⁹ flops. Jacobi: m = ln(10⁻⁶)/ln cos(π/251) = **176 400 sweeps**, each 2 nnz + 3N = 8.1 × 10⁵ flops: **1.4 × 10¹¹ flops**, a hundred times CG.

(d) Dense Cholesky N³/3 = 8.1 × 10¹³ flops, **22.6 hours** at 10⁹ flops/s (and 31 GB to hold it); banded Cholesky N n² = 3.9 × 10⁹, **3.9 s**; CG **1.45 s**; Jacobi **143 s**. At this size banded-direct and CG are comparable; CG's advantage grows as n, and in 3-D the band becomes n² wide and the direct cost n⁷.

**Check.** The CG count ≈ (√κ/2) ln(2/ε) = 79.9 × 14.5 = 1 159, the same. Jacobi ≈ 2(n + 1)² ln(10⁶)/π² = 176 378: the small-angle form of ln cos is accurate to 3 × 10⁻⁵.

**Common mistakes.** Forgetting the vector operations in the CG cost (a factor 2 here); using κ instead of √κ in the CG bound; counting the CSR row pointer as one per nonzero.

## 2. Eigenvalues, Jacobi, and the CG recurrence

(a) With θ = pπ/(n + 1) and v_i = sin iθ: (Tv)_i = 2 sin iθ − sin(i − 1)θ − sin(i + 1)θ = 2 sin iθ − 2 sin iθ cos θ = (2 − 2cos θ) v_i, and the boundary rows work because v₀ = sin 0 = 0 and v_{n+1} = sin pπ = 0. So **λ_p(T) = 2 − 2cos(pπ/(n + 1)) = 4 sin²(pπ/(2(n + 1)))**. For a Kronecker sum, (I ⊗ T + T ⊗ I)(v_q ⊗ v_p) = v_q ⊗ Tv_p + Tv_q ⊗ v_p = (λ_p + λ_q)(v_q ⊗ v_p), giving (2). Extremes: p = q = 1 and p = q = n, with sin(nπ/(2(n + 1))) = cos(π/(2(n + 1))), so **κ = cos²θ₁/sin²θ₁ = cot²(π/(2(n + 1))) ≈ (2(n + 1)/π)²** for large n, the error being the O(θ²) terms of cot θ ≈ 1/θ − θ/3 (0.3 % at n = 24).

(b) D = (4k/h²)I, so G = I − D⁻¹A = I − (h²/4k)A has the same eigenvectors as A with eigenvalues 1 − [sin²(pθ/2) + sin²(qθ/2)] = **(cos pθ + cos qθ)/2**, θ = π/(n + 1). The largest modulus is at p = q = 1 (and, with opposite sign, p = q = n): **ρ_J = cos(π/(n + 1))**. Then ln ρ_J = ln cos θ = −θ²/2 − θ⁴/12 − …, and the number of sweeps for ‖e‖ ≤ ε‖e₀‖ is m = ln ε/ln ρ_J ≈ **2(n + 1)² ln(1/ε)/π²**: quadratic in n.

(c) φ(x + αp) = φ(x) + α pᵀ(Ax − b) + ½α² pᵀAp is a convex quadratic in α (pᵀAp > 0); setting the derivative −pᵀr + α pᵀAp to zero gives **α = pᵀr/pᵀAp**. The new residual is r_new = b − A(x + αp) = r − αAp and pᵀr_new = pᵀr − α pᵀAp = 0. For the direction update: p_mᵀAp_{m−1} = r_mᵀAp_{m−1} + β p_{m−1}ᵀAp_{m−1}. From the residual recurrence, Ap_{m−1} = (r_{m−1} − r_m)/α_{m−1}, so r_mᵀAp_{m−1} = (r_mᵀr_{m−1} − r_mᵀr_m)/α_{m−1} = −r_mᵀr_m/α_{m−1} using r_mᵀr_{m−1} = 0. Also α_{m−1} p_{m−1}ᵀAp_{m−1} = p_{m−1}ᵀr_{m−1} = r_{m−1}ᵀr_{m−1} (since p_{m−1} = r_{m−1} + β'p_{m−2} and p_{m−2}ᵀr_{m−1} = 0). Hence p_mᵀAp_{m−1} = (−r_mᵀr_m + β r_{m−1}ᵀr_{m−1})/α_{m−1}, which vanishes exactly for **β = r_mᵀr_m/r_{m−1}ᵀr_{m−1}**.

**Check.** For n = 6 the code's `poisson_eigen_bounds` and `jacobi_spectral_radius` agree with `eigvalsh` to 10⁻¹⁰ (test 4); at n = 24, ρ_J = 0.99211 is the measured residual ratio.

**Common mistakes.** Taking p = n for λ_min; forgetting that ρ_J belongs to the *smoothest* mode, so Jacobi is a good smoother and a bad solver; in (c), using r_mᵀr_{m−1} = 0 without noting that it comes from pᵀr_new = 0 and p_{m−1} = r_{m−1} + β'p_{m−2}.

## 3. A heated fin with an insulated tip ($s = 15$: q = 2.5 × 10⁶ W/m³)

(a) k/h² = 200/0.025² = 3.2 × 10⁵ W/(m³·K). Rows 1–3: (k/h²)(−T_{i−1} + 2T_i − T_{i+1}) = q with T₀ = 300 K moved right; row 4 with the ghost node T₅ = T₃: (k/h²)(2T₄ − 2T₃) = q, halved to (k/h²)(T₄ − T₃) = q/2 so that a₃₄ = a₄₃.

**A = 3.2 × 10⁵ × [[2, −1, 0, 0], [−1, 2, −1, 0], [0, −1, 2, −1], [0, 0, −1, 1]] W/(m³·K),  b = (9.85 × 10⁷, 2.5 × 10⁶, 2.5 × 10⁶, 1.25 × 10⁶) W/m³.**

A is symmetric, diagonally dominant in every row (|2| ≥ 1 + 1, |1| ≥ 1) and strictly so in row 1, and irreducible (every node is connected to the next); an irreducibly diagonally dominant matrix is nonsingular, and with positive diagonal and real eigenvalues that are ≥ 0 by Gershgorin, all eigenvalues are positive. Alternatively the Cholesky pivots are 2, 2 − 1/2 = 3/2, 2 − 2/3 = 4/3, 1 − 3/4 = 1/4 (times k/h²): all positive, so **A is SPD**. Its eigenvalues are 3.86 × 10⁴, 3.20 × 10⁵, 7.51 × 10⁵, 1.13 × 10⁶, matching (k/h²)(2 − 2cos((2p − 1)π/9)).

(b) From x = 0, r₀ = b, ‖r₀‖ = 9.857 × 10⁷. Jacobi sweep 1: x = r₀/diag(A) = (153.91, 3.91, 3.91, 3.91) K, ‖r‖ = **5.059 × 10⁷**. Sweep 2: x = (155.86, 82.81, 7.81, 7.81) K, ‖r‖ = **3.667 × 10⁷**. The root information moves one node per sweep: after two sweeps node 3 has barely heard of the wall. Steepest descent from 0: α = r₀ᵀr₀/r₀ᵀAr₀ = 1.605 × 10⁻⁶, x₁ = (158.08, 4.01, 4.01, 2.01) K, ‖r₁‖ = **5.189 × 10⁷** (ratio 0.526): one SD step is about as good as one Jacobi sweep here.

(c) G = I − D⁻¹A has ρ(G) = **0.92388** (= cos(π/8), not cos(π/9): the tip row with its diagonal 1 changes the Jacobi spectrum although A's eigenvalues follow the (2p − 1)π/9 law). Predicted sweeps for 10⁻⁶: Jacobi ln(10⁻⁶)/ln 0.92388 = **174**, Gauss–Seidel 87. Measured: **163 and 81** (the right-hand side is dominated by the wall term in row 1, which is not the slowest mode). Steepest descent took 188.

(d) Exact: **T = (327.34, 346.88, 358.59, 362.50) K**, and the continuous solution 300 + (q/k)(Lx − x²/2) at x = 0.025, 0.05, 0.075, 0.1 m gives the same four numbers to all digits. The discrete solution is exact because the second difference of a quadratic is exact and the ghost-node condition T₅ = T₃ is exact for a parabola symmetric about its vertex at the tip. CG must finish in at most 4 steps (N = 4 orthogonal residuals): it does, with residuals 9.86 × 10⁷, 5.19 × 10⁷, 3.65 × 10⁷, 2.98 × 10⁷, 8.5 × 10⁻⁸ and an error of 3 × 10⁻¹³ K.

**Check.** Tip rise q L²/(2k) = 2.5 × 10⁶ × 0.01/400 = 62.5 K: 362.5 K.

**Common mistakes.** Not halving the tip row (the matrix is then unsymmetric and CG is invalid); using ρ_J = 1 − λ_min/(2k/h²) when D is not a multiple of I; expecting Jacobi's count to match the prediction exactly for a right-hand side that is not the slowest mode.

## 4. The insert and the tolerance ($s = 15$: insert conductivity kt/15)

(a) `eigsh`: λ_min = 53.77, λ_max = 7.317 × 10⁴, **κ₂ = 1.36 × 10³** against 440.7 for the uniform plate (3.1 times). Maximum rise **61.09 K** (361.1 K at the centre) against 17.29 K: the insert conducts fifteen times worse, so the heat absorbed on it must cross it to reach the frame.

(b) To 10⁻⁸: **CG 113** iterations (error 6.7 × 10⁻⁶ K), **Jacobi-PCG 59** (2.0 × 10⁻⁶ K), **IC(0)-PCG 29** (1.2 × 10⁻⁵ K); the uniform plate needs 62 with plain CG. Diagonal scaling removes the contrast's share of κ; IC(0) halves the count again, as it does on the uniform plate.

(c)

| tol | iterations | ‖r‖/‖b‖ reached | max error (K) | ‖e‖₂ (K) | bound ‖r‖₂/λ_min (K) |
|---|---|---|---|---|---|
| 10⁻⁴ | 81 | 8.9 × 10⁻⁵ | 0.19 | 1.02 | 53.7 |
| 10⁻⁶ | 100 | 4.7 × 10⁻⁷ | 6.8 × 10⁻⁴ | 4.4 × 10⁻³ | 0.28 |
| 10⁻⁸ | 113 | 4.2 × 10⁻⁹ | 6.7 × 10⁻⁶ | 3.8 × 10⁻⁵ | 2.5 × 10⁻³ |
| 10⁻¹⁰ | 125 | 5.5 × 10⁻¹¹ | 8.9 × 10⁻⁸ | 6.3 × 10⁻⁷ | 3.3 × 10⁻⁵ |

Each decade of tolerance costs 12–19 iterations and buys two decades of error; the bound holds with a margin of about 50–70, because the remaining error is not in the λ_min mode.

(d) With a discretisation error near 0.03 K, **10⁻⁶** is the right tolerance for the temperature (solver error 7 × 10⁻⁴ K, 40 times below the discretisation error; 10⁻⁴ leaves 0.19 K, six times above it). The heat flux into the frame is kt (T₁ − T_wall)/h, which multiplies the error at the first grid line by kt/h = 56 W/(m·K): the rough part of a 10⁻⁶ solution error is harmless, but a safe rule is one more decade, **10⁻⁸**, which costs 13 iterations.

**Common mistakes.** Reading ‖r‖/‖b‖ as the error; using the uniform plate's κ for the insert; solving to 10⁻¹² "to be safe" at 25 extra iterations for nothing.

## 5. Successive over-relaxation ($s = 15$: n = 23)

(a) ρ_J = cos(π/24) = **0.99144**; ω_opt = 2/(1 + sin(π/24)) = **1.76909**; ρ_SOR = ω_opt − 1 = **0.76909**.

(b) To 10⁻⁶: **Jacobi 1 603** (predicted ln(10⁻⁶)/ln 0.99144 = 1 608), **Gauss–Seidel 803** (predicted 804), **SOR 71** (predicted 53), **CG 12**. SOR with the optimal ω changes the exponent: ρ_SOR = ω_opt − 1 ≈ 1 − 2π/(n + 1), so its count is O(n), the same as CG's, though with a constant six times larger here.

(c) ω = 1.0: 803; 1.5: 262; 1.7: 126; 1.769: 71; 1.9: 148; 1.95: 302. The optimum is sharp on the high side (ρ_SOR = ω − 1 for ω > ω_opt, so 1.9 costs twice as much) and gentler below it.

(d) At ω_opt the two eigenvalues of the SOR iteration matrix that come from ρ_J coalesce into a double eigenvalue ω_opt − 1 with a single eigenvector (a 2 × 2 Jordan block), so the error component in that mode decays like m ρ^m, not ρ^m. The observed ratio of successive residuals is ρ(1 + 1/m) → ρ from above (measured 0.785 over the last twenty sweeps against 0.769), and the count that solves m ρ^m = 10⁻⁶ is 71 rather than 53: 34 % more. The spectral radius gives the asymptotic rate, not the count for a finite tolerance.

**Common mistakes.** Testing the ω that minimises *ρ* and concluding the theory is wrong when the count is not ln ε/ln ρ; using SOR on a matrix that is not consistently ordered with the Laplacian's ω_opt.

## 6. Choosing the solver

(a) 4 096 unknowns in 2-D: `spsolve` (sparse-direct, fill O(n² log n)) or CG, both in milliseconds; sparse-direct is simplest.
(b) 500 right-hand sides: sparse Cholesky once, 500 pairs of triangular solves; an iterative method would redo 60–150 iterations per case.
(c) 3-D SPD with 2 × 10⁶ elements: PCG with IC(0) or algebraic multigrid; a direct factorisation's fill in 3-D is O(N^{4/3}) in storage and O(N²) in flops and does not fit.
(d) Nonsymmetric, 10⁵: GMRES (restarted) or BiCGSTAB with an ILU preconditioner; CG is invalid because pᵀAp can vanish.
(e) 50 × 50 SPD in a Newton loop: dense Cholesky, 4 × 10⁴ flops; any iterative method is slower and adds a tolerance to worry about.
(f) 10⁴:1 contrast: κ is multiplied by about the contrast, plain CG stalls; Jacobi preconditioning removes the diagonal's share, IC or multigrid does the rest.
(g) Same SPD system every step: factor once if it fits (then each step is two substitutions); otherwise PCG started from the previous solution, which cuts the count to a handful because the initial error is small and smooth.
