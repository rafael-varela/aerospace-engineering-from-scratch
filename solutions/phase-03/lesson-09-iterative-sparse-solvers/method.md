# P03.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. n = 100 + 10s, N = n², nnz = 5n² − 4n. 2. Dense 8N² bytes; CSR 12 nnz + 4(N + 1) bytes. 3. Matvec 2 nnz flops; CG iteration ≈ 2 nnz + 10N. 4. κ ≈ (2(n + 1)/π)²; CG iterations m = ln(2/ε)/ln((√κ + 1)/(√κ − 1)) ≈ (√κ/2) ln(2/ε); Jacobi m = ln ε/ln cos(π/(n + 1)) ≈ 2(n + 1)² ln(1/ε)/π². 5. Multiply by the per-iteration cost; compare with N³/3 and N n² at 10⁹ flops/s.

## 2.
(a) (Tv)_i = 2 sin iθ − sin(i − 1)θ − sin(i + 1)θ = (2 − 2cos θ) sin iθ with θ = pπ/(n + 1); boundary terms vanish because sin 0 = sin((n + 1)θ) = 0. (I ⊗ T + T ⊗ I)(v_p ⊗ v_q) = (λ_p + λ_q)(v_p ⊗ v_q). λ_min = 8 sin²(π/(2(n + 1)))·k/h², λ_max = 8 cos²(π/(2(n + 1)))·k/h²; ratio cot². cot θ ≈ 1/θ for small θ. (b) G = I − (h²/4k)A has eigenvalues 1 − [sin² + sin²] = (cos pθ + cos qθ)/2; largest at p = q = 1. ln cos θ = −θ²/2 + O(θ⁴); m = ln ε/ln ρ. (c) φ(x + αp) is a quadratic in α with derivative pᵀA(x + αp) − pᵀb = −pᵀr + αpᵀAp. r_new = r − αAp, so pᵀr_new = pᵀr − αpᵀAp = 0. For (c) second part: p_mᵀAp_{m−1} = r_mᵀAp_{m−1} + βp_{m−1}ᵀAp_{m−1}; with Ap_{m−1} = (r_{m−1} − r_m)/α_{m−1} and r_mᵀr_{m−1} = 0, the first term is −r_mᵀr_m/α_{m−1}; with α_{m−1} = p_{m−1}ᵀr_{m−1}/p_{m−1}ᵀAp_{m−1} and p_{m−1}ᵀr_{m−1} = r_{m−1}ᵀr_{m−1} (by induction), the sum is zero exactly when β = r_mᵀr_m/r_{m−1}ᵀr_{m−1}.

## 3.
(a) Rows 1–3: (k/h²)(−T_{i−1} + 2T_i − T_{i+1}) = q, with T₀ = 300 moved to the right; row 4: (k/h²)(−T₃ + T₄) = q/2. A = (k/h²) tridiag(−1, [2, 2, 2, 1], −1). Irreducible (tridiagonal with nonzero off-diagonals), diagonally dominant in every row and strictly in row 1 → nonsingular and, being symmetric with positive diagonal, positive definite; or do the four Cholesky pivots: 2, 3/2, 4/3, 1/4 (times k/h²), all positive. (b) Jacobi: x ← x + r/diag(A), recompute r each sweep; SD: α = rᵀr/rᵀAr with r = b. (c) G = I − D⁻¹A, 4 × 4, `np.linalg.eigvals`; sweeps = ln(10⁻⁶)/ln ρ and ln(10⁻⁶)/ln ρ²; compare with `jacobi` and `gauss_seidel` on `sp.csr_matrix(A)`. (d) `np.linalg.solve`; T(x) at x = h, …, 4h. Second differences of a quadratic are exact and the ghost-node condition T₅ = T₃ is exact for a symmetric quadratic about the tip. CG: at most 4 steps.

## 4.
(a) `kf = plate_conductivity(...)`, `A = poisson_matrix(32, h, kf)`, `b = poisson_rhs(32, h, np.full((32, 32), 2000.), 300., kf)`; `eigsh` for λ_min, λ_max; `spsolve` for the reference. (b) `conjugate_gradient(A, b, tol=1e-8, M=...)`, `hist.size - 1`. (c) Loop over tolerances; `np.linalg.norm(x - x_ref)` against `hist[-1] / lmin`. (d) Compare each tolerance's error with 0.03 K and with 0.03 K × (something like 1/h in relative terms) for the gradient.

## 5.
(a) ρ_J = cos(π/(n + 1)); ω_opt = 2/(1 + sin(π/(n + 1))); ρ_SOR = ω_opt − 1. (b) `jacobi`, `gauss_seidel`, `sor`, `conjugate_gradient` with tol = 1e-6; predictions ln(10⁻⁶)/ln ρ. (c) Loop over ω. (d) Jordan block of size 2 at ω_opt: error ∝ m ρ^m; the apparent contraction ρ(1 + 1/m) exceeds ρ and the number of sweeps for ε solves m ρ^m ≈ ε.

## 6.
(a) 4 096 unknowns in 2-D: `spsolve` (sparse-direct) or CG, either is instant. (b) Factor once (sparse Cholesky) and substitute 500 times; iterative would repeat the work. (c) 3-D, millions: PCG with IC or algebraic multigrid; direct fill is prohibitive. (d) Nonsymmetric: GMRES or BiCGSTAB with ILU. (e) Dense Cholesky: 50³/3 flops, overhead dominates. (f) High contrast → large κ: PCG with Jacobi at least, better IC or multigrid; plain CG would stall. (g) CG with the previous step's solution as x₀, or factor once if the matrix is fixed and small enough.
