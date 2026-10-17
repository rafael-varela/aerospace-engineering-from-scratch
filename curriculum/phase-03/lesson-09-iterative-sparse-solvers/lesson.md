# P03.09  Iterative & Sparse Solvers

> When the matrix has a million rows and five nonzeros per row, elimination fills it in and dense storage does not even fit. Store what is there, multiply by it, and let the solution converge: Jacobi in n² iterations, conjugate gradients in n, preconditioned CG in fewer still.

| | |
|---|---|
| **Learning outcomes** | 1. Assemble the five-point finite-difference matrix of a 2-D heat-conduction (Poisson) problem in sparse (CSR) storage, state its eigenvalues, its condition number κ ≈ (2n/π)² and why dense elimination is the wrong tool. 2. Derive the Jacobi, Gauss–Seidel and steepest-descent iterations and their convergence rates, and the conjugate-gradient method from the minimisation of the energy norm over a Krylov subspace, with its bound in terms of √κ. 3. Implement and test, from scratch, all four methods plus Jacobi and incomplete-Cholesky preconditioning, against `scipy.sparse.linalg` as the oracle, and show that the iteration counts scale as n², n² and n. 4. Choose a solver and a preconditioner for a given matrix from its size, symmetry, conditioning and the number of right-hand sides, and state how the residual tolerance relates to the error in the answer. |
| **Prerequisites** | `P03.08` (LU, Cholesky, condition number, backward error), `P01.16` (eigenvalues of symmetric matrices, Rayleigh quotient), `P01.18` (norms, κ), `P03.04` (second differences). Skill check: you can write the second-difference approximation of T'' and say what a positive definite matrix is. |
| **Engineering question** | A 0.2 m square aluminium skin panel, 2 mm thick, edges held at 300 K by its frame, absorbs 2 kW/m² of solar flux: how hot does its centre get, and how do we solve for it on a 64 × 64 grid (4 096 unknowns), then on a 1 000 × 1 000 grid (10⁶ unknowns), when the matrix cannot be stored densely and must not be eliminated? |
| **Assumptions and limits** | Steady in-plane conduction in a thin plate with uniform absorbed flux, no radiation or convection from the faces, isotropic conductivity (piecewise constant where an insert is modelled), Dirichlet edges; five-point stencil with harmonic-mean face conductivities. All matrices are symmetric positive definite, which is what CG needs; nonsymmetric systems (GMRES, BiCGSTAB) are discussed, not built. Pure-Python loops in Gauss–Seidel and the triangular solves, so timings are indicative only. |
| **Analytical work** | Eigenvalues of the 1-D and 2-D second-difference matrices by separation of variables; κ ≈ (2(n+1)/π)²; the Jacobi iteration matrix and ρ_J = cos(πh); Gauss–Seidel ρ = ρ_J²; steepest descent and the A-norm; CG as the minimiser over the Krylov subspace, its three-term recurrence, finite termination and the √κ bound; why preconditioning works. |
| **Computational deliverable** | `starter/main.py`: `poisson_matrix`, `poisson_rhs`, `jacobi`, `gauss_seidel`, `steepest_descent`, `conjugate_gradient`, `jacobi_preconditioner`, `incomplete_cholesky`, `sparse_triangular_solve`, `ic_preconditioner`, `solve_plate`, `discretisation_errors`, `iterations_vs_n` (the grid, the manufactured solution, the eigenvalue formulas, the CG bound, the A-norm, the conductivity field, `convergence_order` and `scaling_exponent` are given). Tests prove: the matrix equals the Kronecker form I ⊗ T + T ⊗ I with 5n² − 4n nonzeros and is symmetric; the variable-k matrix uses harmonic means and is SPD; a sourceless plate sits at the wall temperature; the eigenvalue formula and ρ_J match `eigvalsh`; Jacobi converges at ρ_J to within 2 × 10⁻⁴ and in the predicted number of iterations ± 5 %; Gauss–Seidel takes 0.42–0.58 of Jacobi's count; every method stops at once from the exact solution; the A-norm error of steepest descent falls monotonically; CG matches `spsolve` and `scipy.sparse.linalg.cg`, terminates in N steps on a 9 × 9 SPD matrix, has a monotone A-norm error below the √κ bound, and needs less than a tenth of steepest descent's iterations; the Jacobi preconditioner reproduces plain CG when k is constant and cuts the count by more than 20 % with a 100:1 insert; IC(0) equals exact Cholesky on a tridiagonal matrix, keeps no fill, refuses an indefinite matrix, and cuts the count further; the sparse triangular solves match `spsolve_triangular`; the plate's centre rise matches the series solution 0.07367 qL²/(kt) to 0.3 % and is symmetric; the discretisation order is 2 ± 0.15; the CG iteration count scales as n^1.0 and Jacobi's as n^1.9. |
| **Practical activity** | **computational experiment**: a solver race at n = 24 (iterations, time and error of Jacobi, Gauss–Seidel, steepest descent, CG and `scipy.sparse.linalg.cg`), a scaling study of iterations against n, a grid-convergence study with a manufactured solution, and a preconditioning study on a plate with a low-conductivity insert. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `iterative_solvers.png`, the hand derivation of the eigenvalues, of ρ_J and of the CG coefficients, the four experiment tables, and a half-page note stating the centre temperature of the panel with its discretisation and solver errors separated. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the eigenvalues and ρ_J are derived; all public tests pass; the observed iteration exponents are explained from κ; the note separates the h² error from the residual-controlled solver error with the numbers of the lesson. |
| **References** | `@ocw-18-335j` Week 8 "Lecture 21: Preconditioning Techniques. The Conjugate-Gradient Method", "Lecture 22: Convergence of Conjugate Gradient" (required-free: CG, its bound, preconditioning); Week 7 "Lecture 18: Krylov Methods and the Arnoldi Algorithm" (required-free: what a Krylov subspace is); Week 7 "Lecture 20: The GMRES Algorithm and Convergence of GMRES and Arnoldi", Week 8 "Lecture 23: Biconjugate Gradient Algorithms" (supplementary-free: the nonsymmetric case); `@scipy-userguide` "Sparse Arrays (scipy.sparse)" (required-free: CSR storage, `scipy.sparse.linalg`); "Sparse eigenvalue problems with ARPACK" (supplementary-free: `eigsh` for the extreme eigenvalues used in part 4). |
| **Extensions** | ★ Implement SOR with ω_opt = 2/(1 + sin πh) and verify ρ_SOR = ω_opt − 1 (Problem 5 starts it). ★ Implement GMRES(m) and solve the convection–diffusion matrix of `P03.11`, which is not symmetric. ★ Write a two-grid cycle (Jacobi smoothing, coarse-grid correction) and show that its iteration count is independent of n. |
| **Time** | unmeasured |

## 1. The engineering problem

`P03.08` ended with a warning: Gaussian elimination costs 2n³/3 flops and n² storage, and a Laplace problem on a 1 000 × 1 000 grid has n = 10⁶. Dense storage is 8 terabytes; elimination is 7 × 10¹⁷ flops. Even a banded Cholesky, which keeps the band of width 1 000, costs 10⁶ × 10⁶ = 10¹² flops and fills the band completely. Yet the matrix itself has only five nonzeros per row, 5 × 10⁶ numbers in all, and multiplying a vector by it costs 10⁷ flops. **Iterative methods** solve Ax = b using nothing but that multiplication: they never modify A, never fill it in, and stop when the residual is small enough for the purpose. This lesson builds the four classical ones and shows, by theory and by measurement, that the conjugate-gradient method with a preconditioner is the one that scales.

The vehicle is a square aluminium skin panel, L = 0.2 m, thickness t = 2 mm, conductivity k = 170 W/(m·K), whose edges are clamped to a frame at 300 K and whose face absorbs q = 2 kW/m² of solar flux. Steady in-plane conduction gives −(kt) ∇²T = q, a Poisson equation, whose five-point discretisation on an n × n interior grid is a sparse SPD system of n² unknowns. The series solution for a square gives a centre rise of 0.07367 qL²/(kt) = **17.334 K**; the code reaches 17.186, 17.295 and 17.324 K at n = 16, 32, 64, converging as h². A second version of the panel carries a square insert of lower conductivity (a composite patch, or a hole filled with potting), which makes the matrix ill-conditioned in a way that plain CG feels and a preconditioner removes.

Success criteria, written first: every solver reaches the `spsolve` solution to the requested residual tolerance (error below 2 × 10⁻⁶ K at tol = 10⁻⁶ relative); Jacobi's observed contraction equals cos(πh) to four digits and Gauss–Seidel's iteration count is half of Jacobi's; CG's A-norm error lies under the √κ bound at every step; the iteration counts scale as n^2 (Jacobi, steepest descent) and n^1 (CG) within ± 0.25 in the exponent; the discretisation error scales as h^2 ± 0.15; and the panel centre rise agrees with the series solution to 0.3 % at n = 32.

## 2. Understand and derive

### 2.1 The five-point matrix and why it is sparse

Divide the square into (n + 1)² cells of side h = L/(n + 1) and number the n² interior points row by row, p = i + nj. Central second differences in x and y give, for constant conductivity k (here k stands for the conductance kt of the plate),

**(k/h²)(4T_p − T_{p−1} − T_{p+1} − T_{p−n} − T_{p+n}) = q_p,**   (1)

with a neighbour on the boundary replaced by T_wall and moved to the right-hand side: that is `poisson_rhs`. Row p has at most five nonzeros; the matrix has 5n² − 4n of them in n⁴ entries, 0.84 % at n = 24, 0.48 % at n = 32, and 5 × 10⁻⁴ % at n = 1 000. `scipy.sparse` stores it in **CSR** form (`@scipy-userguide`, "Sparse Arrays"): one array of the nonzero values, one of their column indices, and one of n² + 1 row pointers, 12 bytes per nonzero, 60 kB at n = 32 against 8.4 MB dense. A matrix–vector product costs 2 nnz ≈ 10n² flops, and that is the only operation every method in this lesson performs with A.

The matrix is the Kronecker sum A = (k/h²)(I ⊗ T + T ⊗ I) with T = tridiag(−1, 2, −1), and this is what makes its spectrum known. The 1-D matrix T has eigenvectors sin(pπ i h/L) and eigenvalues 2 − 2cos(pπ/(n + 1)) = 4 sin²(pπ/(2(n + 1))), p = 1, …, n (check: T applied to the sine gives 2 sin(piθ) − sin((i−1)pθ) − sin((i+1)pθ) = (2 − 2cos pθ) sin(ipθ)). Eigenvalues of a Kronecker sum add, so

**λ_{p,q} = (4k/h²)[sin²(pπ/(2(n + 1))) + sin²(qπ/(2(n + 1)))],  κ₂ = λ_max/λ_min = cot²(π/(2(n + 1))) ≈ (2(n + 1)/π)²**   (2)

(`poisson_eigen_bounds`). At n = 24 the extreme eigenvalues are 492.8 and 1.245 × 10⁵ (formula and `eigsh` agree to five digits), κ = 252.6 against the approximation 253.3; at n = 32, κ = 440.7. The matrix is SPD (all λ > 0) and κ grows as n²: the finer the grid, the harder the system, and every convergence rate below is a function of κ.

With a variable conductivity the face between cells of conductivities k_a and k_b gets the harmonic mean 2k_a k_b/(k_a + k_b), which is the series resistance of two half-cells and keeps the matrix symmetric. The insert of part 4 with k ratio 100 raises κ to 8.0 × 10³.

### 2.2 Stationary iterations: Jacobi and Gauss–Seidel

Split A = D − L − U (diagonal, strictly lower, strictly upper). **Jacobi** solves each equation for its own unknown using the old values of the others: x^{(m+1)} = x^{(m)} + D⁻¹ r^{(m)} with r = b − Ax, one matvec and one vector division per sweep. The error obeys e^{(m+1)} = G e^{(m)} with G = I − D⁻¹A, so convergence is governed by the spectral radius ρ(G). For (1), D = (4k/h²) I and G has the same eigenvectors as A with eigenvalues 1 − λ_{p,q}h²/(4k) = [cos(pπh/L) + cos(qπh/L)]/2, so

**ρ_J = cos(π/(n + 1)) ≈ 1 − π²/(2(n + 1)²).**   (3)

The error falls by ρ_J per sweep; m sweeps for a reduction ε need m = ln ε/ln ρ_J ≈ (2(n + 1)²/π²) ln(1/ε): **n² sweeps**, each costing n² flops, n⁴ in all, the same as a banded direct solve. At n = 24, ρ_J = 0.99211 and the prediction for ε = 10⁻⁶ is 1 745 sweeps; the code took 1 740 and its measured residual ratio over the last fifty sweeps was 0.99211. The smooth error modes (small p, q) are the ones that survive: ρ_J belongs to p = q = 1, the half-sine over the whole plate. Jacobi kills the rough error in a few sweeps and then crawls.

**Gauss–Seidel** uses each new value as soon as it exists: sweep the rows in order, x_i ← (b_i − Σ_{j≠i} a_ij x_j)/a_ii with the x_j already updated for j < i. In splitting form x^{(m+1)} = (D − L)⁻¹(U x^{(m)} + b), a forward substitution per sweep. For matrices like (1) (consistently ordered), ρ_GS = ρ_J²: **half the sweeps** of Jacobi, 873 predicted and 872 measured at n = 24. The gain is a constant factor, not a change of exponent; the n² remains. SOR, which over-relaxes the Gauss–Seidel update by ω ≈ 2/(1 + sin πh), does change the exponent to n (Extension 1, Problem 5), but it needs the right ω, which for a general matrix nobody knows.

### 2.3 Minimisation: steepest descent and the energy norm

For SPD A, solving Ax = b is the same as minimising φ(x) = ½xᵀAx − bᵀx, because ∇φ = Ax − b = −r. The error e = x − x* satisfies φ(x) − φ(x*) = ½ eᵀAe = ½‖e‖²_A: the **A-norm** (energy norm) of the error is what a minimiser of φ reduces. **Steepest descent** moves along −∇φ = r and chooses the step that minimises φ on that line: α = rᵀr/rᵀAr, then x ← x + αr, r ← r − αAr. One matvec per step, no matrix modification. Its A-norm error falls by at least (κ − 1)/(κ + 1) per step, which for κ = (2n/π)² is again 1 − O(1/n²): **n² iterations**, the same exponent as Jacobi (measured exponent 2.03 against Jacobi's 1.89 over n = 8–32, with 166, 696 and 2 760 iterations). The reason is the zigzag: each new direction is orthogonal to the previous one, so the method undoes part of what it did two steps before. The residual norm itself is not even monotone; only ‖e‖_A is, which is why the test watches the A-norm.

### 2.4 Conjugate gradients: the right directions

The cure is to search in directions that do not fight each other. After m steps, steepest descent and CG have both only ever formed the vectors b, Ab, A²b, …: the **Krylov subspace** K_m = span{r₀, Ar₀, …, A^{m−1}r₀} (`@ocw-18-335j` Lecture 18). CG is the method whose m-th iterate is *the* minimiser of φ (equivalently of ‖e‖_A) over x₀ + K_m. Three facts, derived in Lecture 21 and in Problem 2:

1. **A-orthogonal directions.** Choose directions p₀, …, p_{m−1} with p_iᵀAp_j = 0 (i ≠ j) spanning K_m. Then minimising over K_m decouples into m one-dimensional minimisations, α_j = p_jᵀr_j/p_jᵀAp_j, and the minimiser is reached by m successive line searches. Each residual is orthogonal to all previous directions: r_mᵀp_j = 0 for j < m, hence r_mᵀr_j = 0 too.
2. **The two-term recurrence.** Take p_m = r_m + β_m p_{m−1} with β_m = r_mᵀr_m/r_{m−1}ᵀr_{m−1}. This p_m is A-orthogonal to *all* earlier directions, not just the last one, because Ar_m lies in K_{m+2} and the residuals are mutually orthogonal. So CG needs one matvec, two inner products and three vector updates per step, and stores four vectors, however large m becomes:

   **α = rᵀr/pᵀAp,  x ← x + αp,  r ← r − αAp,  β = r_newᵀr_new/rᵀr,  p ← r_new + βp.**   (4)

3. **Finite termination and the κ bound.** Since the residuals are orthogonal, at most N of them are nonzero: CG is a direct method in exact arithmetic, finishing in N steps (the test checks it on a 9 × 9 matrix). The useful statement is the other one: because x_m minimises ‖e‖_A over the Krylov subspace, ‖e_m‖_A = min over polynomials P of degree m with P(0) = 1 of ‖P(A)e₀‖_A ≤ max_λ |P(λ)| ‖e₀‖_A, and the Chebyshev polynomial scaled to [λ_min, λ_max] gives (`@ocw-18-335j` Lecture 22)

   **‖e_m‖_A/‖e₀‖_A ≤ 2 [(√κ − 1)/(√κ + 1)]^m.**   (5)

For (1), √κ ≈ 2n/π, so the factor per step is about 1 − π/n and a reduction ε needs **m ≈ (n/π) ln(2/ε) ∝ n** steps (`cg_error_bound`): one power of n better than Jacobi, Gauss–Seidel and steepest descent. The code took 5, 9, 19, 44 iterations at n = 8, 16, 32, 64 for 10⁻⁶ (exponent 1.05), against Jacobi's 222, 802, 3 034. At n = 24, CG needed 13 iterations and 0.2 ms where Jacobi needed 1 740 and 18 ms and steepest descent 1 554. The bound (5) is pessimistic when the eigenvalues are clustered or the right-hand side is smooth: after 40 iterations at n = 24 the bound allows 1.3 × 10⁻², the actual A-norm error is 3.3 × 10⁻¹⁵, round-off. What (5) buys is the exponent, not the constant.

### 2.5 Preconditioning

Everything above depends on κ. **Preconditioning** replaces Ax = b by M⁻¹Ax = M⁻¹b with M ≈ A cheap to apply; CG on the preconditioned system (kept symmetric by working with M^{−1/2}AM^{−1/2}, which the algorithm never forms) is (4) with z = M⁻¹r replacing r in the inner products and in the direction update (`@ocw-18-335j` Lecture 21). The iteration count then follows κ(M⁻¹A). Two preconditioners are built:

- **Jacobi**, M = D. For constant k it is a scalar and changes nothing (the test checks that the residual histories coincide); for the plate with a 100:1 insert the diagonal varies by a factor 100 across the patch boundary and M = D brings the count from 180 to 58 at tolerance 10⁻⁸.
- **Incomplete Cholesky IC(0)**: run the Cholesky recursion of `P03.08` but keep only the entries where A has nonzeros, discarding all fill. The result L Lᵀ ≈ A costs the same storage as A and two sparse triangular solves per application. On the plate it brings the count to 29–30 whether or not there is an insert; the factorisation itself costs about as much as a dozen CG iterations. IC(0) can break down (a non-positive pivot) even on an SPD matrix, which is why `incomplete_cholesky` raises and practice uses shifted or modified variants.

The library offers the same algorithm as `scipy.sparse.linalg.cg(A, b, M=...)`, and the extreme eigenvalues that the κ column of part 4 uses come from `scipy.sparse.linalg.eigsh`, ARPACK's Lanczos iteration, which is itself a Krylov method (`@scipy-userguide`, "Sparse eigenvalue problems with ARPACK"). For nonsymmetric matrices, such as the convection–diffusion problems of `P03.11`, CG's short recurrence is lost; GMRES keeps the whole Krylov basis and minimises the residual, BiCGSTAB keeps a short recurrence at the price of a non-monotone residual (`@ocw-18-335j` Lectures 20 and 23).

### 2.6 What the residual tolerance means

All solvers stop at ‖r‖ ≤ tol·‖b‖. By `P03.08`, ‖e‖/‖x‖ ≤ κ(A)·‖r‖/(‖A‖‖x‖), so a residual tolerance of 10⁻⁶ on a matrix with κ = 253 can leave a relative error up to about 2.5 × 10⁻⁴; the measured errors at n = 24 are 1.6 × 10⁻⁶ K for Jacobi and Gauss–Seidel, 1.1 × 10⁻⁶ K for steepest descent and 1.1 × 10⁻⁹ K for CG against `spsolve`, far below the bound, because the error that remains is a smooth mode whose A-norm is small. The discretisation error at n = 24 is of order 10⁻³ K (part 4: 1.6 × 10⁻², 4.5 × 10⁻³, 1.2 × 10⁻³, 3.1 × 10⁻⁴ for the manufactured solution at n = 8–64). **The solver tolerance should be chosen below the discretisation error and no further**; solving the discrete system to 10⁻¹² when the discrete system is wrong at 10⁻³ is wasted work, which is why engineering codes default to 10⁻⁶ or 10⁻⁸.

### 2.7 Where the methods stop being valid

CG requires A to be SPD: on an indefinite or nonsymmetric matrix the inner product pᵀAp can vanish or go negative and the method fails without warning. Jacobi requires ρ(G) < 1, which diagonal dominance guarantees but a general matrix does not. Gauss–Seidel converges for any SPD matrix; Jacobi does not. All iteration counts above assume the smooth, low-frequency error of the Laplacian; a problem with a strongly varying coefficient or anisotropy can push κ far above (2n/π)² and the counts with it, which is when a preconditioner stops being an optimisation and becomes the method. The harmonic-mean stencil is second order only where k is smooth; at a jump it is first order locally, though the global error of part 4 still fell as h².

## 3. Build the computational model

Open `starter/main.py`. `grid`, `manufactured_solution`, `poisson_eigen_bounds`, `jacobi_spectral_radius`, `cg_error_bound`, `a_norm`, `plate_conductivity`, `convergence_order`, `scaling_exponent` and `main` are given.

### Step 1: the matrix

`poisson_matrix(n, h, k)` loops over the interior points, computes the four face conductivities (harmonic means from the (n + 2) × (n + 2) field, or the constant), appends −k_face/h² for each interior neighbour and their sum/h² on the diagonal, and builds `sp.csr_matrix((vals, (rows, cols)))`. `poisson_rhs` adds k_face T_wall/h² for every neighbour that falls outside the grid. Check the Kronecker identity of the test before going further.

### Step 2: stationary methods

```python
def jacobi(A, b, x0=None, tol=1e-8, maxiter=20000):
    d = A.diagonal(); x = np.zeros(A.shape[0]) if x0 is None else np.array(x0, float)
    hist = []
    for _ in range(maxiter):
        r = b - A @ x; hist.append(np.linalg.norm(r))
        if hist[-1] <= tol * np.linalg.norm(b): break
        x = x + r / d
    return x, np.array(hist)
```

`gauss_seidel` walks the CSR arrays `indptr`, `indices`, `data` row by row and updates `x[i]` in place (a vectorised version is impossible: the sweep is sequential). Record the residual before each sweep, as above, so that `hist.size - 1` is the number of sweeps taken and a converged initial guess records exactly one entry.

### Step 3: Krylov methods

`steepest_descent` and `conjugate_gradient` are (4) with the residual updated recursively (`r -= alpha * A @ p`), one matvec per iteration; the preconditioned version applies `M(r)` once per iteration and uses `r @ z` in α and β. Keep `hist[0] = ‖r₀‖`.

### Step 4: preconditioners

`jacobi_preconditioner` returns a closure dividing by the diagonal. `incomplete_cholesky` runs the recursion (1) of `P03.08` on a dense work copy restricted to the pattern `A != 0`, raises `ValueError` on a non-positive pivot, and returns the lower factor as CSR. `sparse_triangular_solve` is forward or back substitution over CSR rows; `ic_preconditioner` composes the two solves.

### Step 5: the plate and the studies

`solve_plate` assembles and calls CG to 10⁻¹², returning the grid, T as (n, n) and the count. `discretisation_errors` solves the manufactured problem (whose exact solution is (1 + x/L) sin(πx/L) sin(πy/L), chosen *not* to be an eigenvector of A, since CG would then converge in one step) at each n and returns the max-norm errors; `iterations_vs_n` returns iteration counts for a given method.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-09-iterative-sparse-solvers`. The tests:

1. **Grid** spacing L/(n + 1) and layout X[j, i] = x_i.
2. **Matrix** equals (k/h²)(I ⊗ T + T ⊗ I) to 10⁻¹², 5n² − 4n nonzeros, symmetric; variable k: symmetric, harmonic-mean face value 3.2 between k = 2 and 8, positive eigenvalues, constant field equals constant k.
3. **Right-hand side**: no source gives T = T_wall to 10⁻¹².
4. **Eigenvalues** (2) and ρ_J = cos(π/7) against `eigvalsh` at n = 6.
5. **Jacobi** at n = 10: converged, error < 10⁻⁶, measured ratio = ρ_J ± 2 × 10⁻⁴, count within 5 % of ln(ε‖b‖/‖r₀‖)/ln ρ_J.
6. **Gauss–Seidel** at n = 8: converged, 0.42–0.58 of Jacobi's sweeps.
7. **Initial guess**: all four methods record one residual and take no step from the exact solution.
8. **Steepest descent**: converged, A-norm error strictly decreasing over 20 steps.
9. **CG**: equals `spsolve` and `scipy.sparse.linalg.cg` at n = 12 in fewer than 40 steps; finishes in 9 steps on a 9 × 9 SPD matrix; A-norm error monotone and under (5) for 10 steps at n = 12; `cg_error_bound(1, 5) = 0`; at least ten times fewer iterations than steepest descent.
10. **Jacobi preconditioner**: same histories as plain CG for constant k; **IC(0)** equals `scipy.linalg.cholesky` on a tridiagonal matrix, lower triangular, no fill, refuses [[1, 2], [2, 1]]; **triangular solves** equal `spsolve_triangular`; with a 100:1 insert at n = 12, Jacobi-PCG takes under 80 % of CG's count and IC-PCG under 85 % of that and under 60 % of CG's, all with error below 10⁻⁵ of the temperature rise.
11. **Plate** at n = 32: shape, centre rise = 0.07367 qL²/(kt) to 0.3 %, symmetric about both diagonals, 40–120 CG iterations; the conductivity field marks the insert.
12. **Discretisation order** 2 ± 0.15 for n = 6, 12, 24; `convergence_order` refuses non-positive errors.
13. **Scaling**: CG exponent 0.8–1.25, Jacobi 1.7–2.2 over n = 8, 16, 32; Jacobi more than 50 times CG at n = 32.
14. `main` writes the figure.

### Compare

`scipy.sparse.linalg.cg` is the same algorithm (4) in compiled vector operations with the matvec dispatched to the CSR kernel; it returned the same solution as ours at n = 24 to 10⁻⁹ K and took 0.6 ms against our 0.2 ms (both dominated by overhead at this size: 13 iterations of a 576-row matvec). `spsolve` is SuperLU, a sparse-direct solver with fill-reducing ordering: at n = 24–64 it is the fastest of all, and it is what one should call for a 2-D problem of this size, because the fill of a 2-D nested-dissection ordering is only O(n² log n). The iterative methods win in 3-D, where sparse-direct fill grows as n⁴ for n³ unknowns, for matrix-free problems where A is never formed, and when a good initial guess (last time step) is available. We write them so that "CG converges in O(√κ) iterations and a preconditioner is a change of κ" is something we have measured, and because the preconditioner is where all the engineering judgment sits.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **The race.** Solve the manufactured problem at n = 24 to 10⁻⁶ with Jacobi, Gauss–Seidel, steepest descent, CG and `scipy.sparse.linalg.cg`; record iterations, time and max error against `spsolve`. Predict: 1 745, 873, O(10³), O(10¹) iterations from (3), ρ_J² and (5).
2. **Scaling.** Iterations for 10⁻⁶ at n = 8, 16, 32 (64 for CG); fit exponents. Predict 2, 2, 1.
3. **Grid convergence.** Max error against the manufactured solution at n = 8, 16, 32, 64, CG to 10⁻¹²; predict order 2.
4. **Preconditioning.** The plate with a 0.08 m square insert of conductivity k/1, k/10, k/100: κ from `eigsh`, iterations of CG, Jacobi-PCG, IC-PCG to 10⁻⁸, max temperature rise. Predict: κ and the plain-CG count grow with the ratio; the preconditioned counts do not.

Raw data format: CSV tables (`method, iterations, time_ms, max_err_K`; `n, it_jacobi, it_sd, it_cg`; `n, h, max_err`; `k_ratio, kappa, it_cg, it_jpcg, it_icpcg, max_rise_K`), each with a header naming the fixed parameters.

What you should find: 1 740, 872, 1 554, 13 and 13 iterations (18, 1 086, 22, 0.2 and 0.6 ms; the Gauss–Seidel time is the pure-Python row loop), errors 1.6 × 10⁻⁶, 1.5 × 10⁻⁶, 1.1 × 10⁻⁶, 1.1 × 10⁻⁹, 1.1 × 10⁻⁹ K; exponents 1.89, 2.03, 1.05; errors 1.57 × 10⁻², 4.47 × 10⁻³, 1.19 × 10⁻³, 3.07 × 10⁻⁴, order 1.99; κ = 4.4 × 10², 9.8 × 10², 8.0 × 10³ with CG 62, 102, 180, Jacobi-PCG 62, 59, 58, IC-PCG 30, 29, 29, and rises 17.29, 45.54, 325.25 K (the last is unphysical for a real panel, which would radiate; it is the linear model's answer).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| nnz at n = 24 | 5n² − 4n = 2 784 | 2 784 | 0 | exact |
| λ_min, λ_max at n = 24 | 492.8, 1.245 × 10⁵ | 492.8, 1.245 × 10⁵ (`eigsh`) | < 10⁻⁴ | exact formula |
| κ at n = 24 | (2·25/π)² = 253.3 | 252.6 | 0.3 % | small-angle approximation in κ |
| ρ_J at n = 24 | cos(π/25) = 0.99211 | 0.99211 (last 50 sweeps) | < 10⁻⁵ | exact |
| Jacobi sweeps for 10⁻⁶ | 1 745 | 1 740 | 0.3 % | initial residual has rough components that die faster |
| Gauss–Seidel sweeps | 873 | 872 | 0.1 % | ρ_GS = ρ_J² |
| CG A-norm error after 40 steps | ≤ 1.3 × 10⁻² (bound) | 3.3 × 10⁻¹⁵ | — | bound (5) is a worst case |
| CG exponent | 1 | 1.05 | 0.05 | √κ ∝ n; small-n constants |
| Jacobi / SD exponents | 2 / 2 | 1.89 / 2.03 | ≤ 0.11 | pre-asymptotic at n = 8 |
| discretisation order | 2 | 1.99 | 0.01 | asymptotic |
| plate centre rise (K) | 17.334 (series) | 17.186, 17.295, 17.324 at n = 16, 32, 64 | 0.9, 0.2, 0.06 % | h² |
| Jacobi-PCG for constant k | = CG | = CG (62 = 62) | 0 | M = scalar |

Everything in the table is a check against an exact statement except the last two rows of part 4 (the preconditioned counts), which are measurements; they are also the rows that answer the engineering question of scale. The plate's centre temperature is **317.3 K** at n = 64, with a discretisation error of about 0.01 K (Richardson from the three grids: the differences 0.109 and 0.029 K fall by 3.8, consistent with h²) and a solver error of 10⁻¹⁰ K. The result is limited by the model (no radiation, uniform flux) long before the arithmetic.

Improvements worth making: a red–black ordering for Gauss–Seidel, which vectorises the sweep and halves its cost without changing ρ; a shifted IC (add αD before factoring) that never breaks down; `scipy.sparse.linalg.cg` with our preconditioner as a `LinearOperator`, which keeps our M and the library's speed. Improvement worth *not* making: solving to 10⁻¹² "to be safe": the discretisation error is 10⁻³ and the extra 25 iterations buy nothing.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (2), (3), (4) and (5); `iterative_solvers.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the panel centre temperature, its discretisation error from the three grids, the solver tolerance used and why it is sufficient, and which solver you would use at 10⁶ unknowns.

### Adversarial self-checks

- **What if the flux doubles?** The model is linear: the rise doubles to 34.7 K, and no iteration count changes (counts depend on A, not b, apart from the rough content of r₀). What changes is the validity of ignoring radiation: at 350 K an aluminium face radiates εσT⁴ ≈ 85 W/m² for ε = 0.1 (bare) and 400 W/m² for an anodised ε = 0.5: 4–20 % of the absorbed flux, which the model omits.
- **What if the tolerance is loosened?** On the plate with a 15:1 insert at n = 32 (Problem 4), CG at 10⁻⁴ stops after 81 iterations with a maximum error of 0.19 K, at 10⁻⁶ after 100 with 6.8 × 10⁻⁴ K, at 10⁻⁸ after 113 with 6.7 × 10⁻⁶ K. The discretisation error of a 32-grid is about 0.03 K (17.295 against 17.324 K at n = 64), so 10⁻⁴ is too loose, 10⁻⁶ is enough for the temperature, and the gradient (the heat flux into the frame), which amplifies the rough part of the error by 1/h, wants one more decade.
- **Which assumption fails first in the real world?** Uniform conductivity: rivets, doublers and a composite patch change k by factors of 10–100 and κ with it, which is what part 4 models; then the Dirichlet edge, since a real frame has a finite conductance.
- **What measurement could falsify the model?** A thermocouple at the centre: 317.3 K predicted; a reading of 312 K would point to radiation or a warmer edge than 300 K, not to the solver.
- **Does the result survive a change of solver?** Jacobi, Gauss–Seidel, steepest descent, CG, PCG, `scipy.sparse.linalg.cg` and `spsolve` agree to the residual tolerance; CG and `spsolve` to 10⁻⁹ K. A disagreement beyond κ·tol would mean a bug in the matrix, not in the solver.
- **Which decision changes if the uncertainty doubles?** If k is known to 10 % instead of 5 %, the centre rise is known to 10 %, and solving to 10⁻⁶ rather than 10⁻³ changes nothing. The decision that does depend on the numerics is the choice of solver at 10⁶ unknowns: CG with IC at a few hundred iterations against a sparse-direct solve with O(n² log n) fill, and in 3-D only the iterative choice survives.

## Check yourself

1. How many nonzeros does the five-point matrix have for an n × n grid, and how many flops does one matrix–vector product cost?
2. Write the eigenvalues of the 2-D five-point matrix and the approximation κ ≈ (2n/π)².
3. What is the Jacobi iteration in splitting form, and why is its spectral radius cos(πh) for the Laplacian?
4. Why does Gauss–Seidel take half the sweeps of Jacobi, and why does neither change the n² scaling?
5. Define the A-norm of the error and explain what steepest descent and CG each minimise.
6. State the CG error bound and derive the O(n) iteration count from it.
7. What does a preconditioner change, and why does Jacobi preconditioning do nothing for constant k?
8. A CG solve reports ‖r‖/‖b‖ = 10⁻⁶ on a matrix with κ = 10⁴. What is the worst-case relative error of the solution, and how should the tolerance be related to the discretisation error?
