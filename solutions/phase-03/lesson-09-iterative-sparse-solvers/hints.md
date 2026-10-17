# P03.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. nnz = 5n² − 4n; a matvec is 2 nnz flops; iterations are ln(2/ε)/ln((√κ + 1)/(√κ − 1)) for CG and ln ε/ln ρ_J for Jacobi; the totals are iterations × per-iteration cost. The dense answer is in hours, the others in seconds.
2. (a) Use sin((i ± 1)θ) = sin iθ cos θ ± cos iθ sin θ; Kronecker-sum eigenvalues add; cot² = (1 − sin²)/sin². (b) D = (4k/h²)I, so G = I − h²A/(4k) and its eigenvalues follow from (2); ln cos θ ≈ −θ²/2. (c) dφ(x + αp)/dα = pᵀ(A(x + αp) − b) = 0; substitute Ap_{m−1} from the residual recurrence into p_mᵀAp_{m−1}.
3. Row 4 before halving is (k/h²)(2T₄ − 2T₃) = q; after halving, (k/h²)(T₄ − T₃) = q/2. Diagonal dominance is strict in row 1 only and weak elsewhere, but the matrix is irreducible, which suffices. The continuous solution is a quadratic and second differences of a quadratic are exact.
4. `eigsh(A, k=1, which="SA", return_eigenvectors=False)`. The bound comes from ‖e‖₂ = ‖A⁻¹r‖₂ ≤ ‖r‖₂/λ_min. For the flux you differentiate the temperature, which multiplies the rough error by 1/h.
5. Copy `gauss_seidel`, blend the old and new value with ω. ρ_SOR = ω_opt − 1 is the spectral radius, but the iteration matrix at ω_opt is not diagonalisable, so the error behaves like m ρ^m rather than ρ^m: the observed ratio approaches ρ from above and the count is larger than ln ε/ln ρ.
6. Match: symmetric or not, size, fill of a direct factorisation (2-D against 3-D), number of right-hand sides, contrast in coefficients (κ), and whether a good initial guess exists.
