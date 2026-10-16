# P03.08  Direct Linear Solvers

> Gaussian elimination is a factorisation: A = LU. Do it once, keep the factors, and every new load case costs n² instead of n³. The residual tells you the solver worked; the condition number tells you whether the answer means anything.

| | |
|---|---|
| **Learning outcomes** | 1. Derive LU factorisation from Gaussian elimination, explain why partial pivoting is needed and what the growth factor bounds, and count the cost of factorisation, substitution, inversion and the tridiagonal (Thomas) algorithm. 2. Derive the Cholesky factorisation and show that its pivots are positive exactly when the matrix is positive definite, so that it is both the cheapest solver and the test for a stiffness or conductance matrix. 3. Implement and test, from scratch, LU with pivoting, Cholesky, the triangular substitutions, the Thomas algorithm, determinants, inverses, condition numbers, backward error and iterative refinement, against `scipy.linalg` as the oracle. 4. Distinguish backward from forward error, predict the forward error from κ(A)·ε, and decide when a computed solution of an ill-conditioned system can be trusted. |
| **Prerequisites** | `P01.13` (Gaussian elimination by hand, inverses, determinants), `P01.16` (eigenvalues, symmetric matrices), `P01.17` (least squares, QR), `P01.18` (norms, SVD, condition number), `P03.01` (floating point, machine epsilon). Skill check: you can solve a 3×3 system by elimination and say what a positive definite matrix is. |
| **Engineering question** | A cantilevered wing-spar truss (2 m span, 0.3 m depth, four bays, aluminium bars of 2 cm²) carries 4 kN at the tip: how far does the tip deflect, which member is most stressed, how many flops does the answer cost, how much less does the *next* load case cost, and how many digits of the answer are real? |
| **Assumptions and limits** | Pin-jointed bars, linear elasticity, small displacements, loads at joints only (the truss model, `P05` will do it properly); a 1-D rod with constant conductivity and a given volumetric source; exact arithmetic in the derivations and IEEE double in the code. Dense storage throughout (the sparse-direct and iterative families are `P03.09`); partial pivoting only (complete pivoting and rook pivoting are mentioned, not built). |
| **Analytical work** | LU from elimination; existence and the role of pivoting; flop counts 2n³/3, n², 2n³, n³/3, 8n; the growth factor and Wilkinson's matrix; Cholesky column by column and its pivot test; the bar stiffness matrix and the assembly rule; statics of the root bay by hand; the backward-error definition and the bound ‖δx‖/‖x‖ ≤ κ(A)·(backward error); the exact solution of the heated rod. |
| **Computational deliverable** | `starter/main.py`: `forward_substitution`, `back_substitution`, `lu_decompose`, `lu_solve`, `solve`, `determinant_from_lu`, `inverse_from_lu`, `growth_factor`, `cholesky`, `cholesky_solve`, `thomas_solve`, `condition_number_1`, `backward_error`, `iterative_refinement`, `truss_stiffness`, `apply_supports`, `member_forces`, `solve_truss`, `heated_rod_errors` (the truss geometry, the rod discretisation, `hilbert`, `wilkinson_matrix`, `norm1`, `convergence_order` are given). Tests prove: the substitutions match `solve_triangular`; A[perm] = LU with unit L and multipliers ≤ 1, factors equal to `scipy.linalg.lu`, solutions to 10⁻¹⁰ and backward error below 10⁻¹⁵; unpivoted elimination fails where pivoted succeeds and singular matrices raise; determinants, inverses and κ₁ match NumPy; the growth factor is exactly 2ⁿ⁻¹ on Wilkinson's matrix; Cholesky reconstructs, matches SciPy, and refuses unsymmetric, indefinite and singular matrices; Thomas matches a dense solve; on Hilbert matrices the backward error stays at round-off while the forward error grows with slope 0.94 in κ and stays below κε; refinement reduces the error; the truss stiffness has exactly three zero modes, is SPD once supported, satisfies equilibrium and the hand statics of the root bay, matches `scipy.linalg.solve`, and superposes; the heated rod converges at order 2.00. |
| **Practical activity** | **computational experiment**: a cost study (our LU against `scipy.linalg.lu_factor` at n = 100–400, Thomas against dense LU at n = 600), a conditioning study on Hilbert matrices, a convergence study on the rod, and a load-case study on the truss. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `direct_solvers.png`, the hand derivation of the Cholesky pivot test and of the root-bay member forces, the Hilbert table, and a half-page note stating how many digits of the truss deflection are trustworthy and why. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: LU and Cholesky are derived with their costs; all public tests pass; the root-bay forces by statics agree with the code; the note distinguishes backward from forward error with the numbers of the lesson. |
| **References** | `@ocw-18-335j` Week 5 "Lecture 13: LU Factorization and Partial Pivoting", "Lecture 14: Cholesky Factorization and other Specialized Solvers" (required-free: pivoting, growth factor, cost); Week 3 "Lecture 7: The SVD, its Applications, and Condition Numbers" (required-free: κ and backward error); `@ocw-18-085` "Lecture 3: Solving a Linear System", "Lecture 7: Positive Definite Day" (required-free: elimination as LU, the four tests for positive definiteness); `@ocw-2-092` "Physical explanation of Gauss elimination" (SES 13) (supplementary-free: elimination as static condensation of a structure); `@scipy-userguide` "Linear Algebra (scipy.linalg)" (required-free); `@ocw-18-335j` Week 9 "Lecture 24: Sparse-Direct Solvers" (supplementary-free: what changes when A is sparse). |
| **Extensions** | ★ Implement banded LU for a matrix of bandwidth b in O(n b²) and verify on the rod with a five-point stencil. ★ Implement LDLᵀ for symmetric indefinite matrices with Bunch–Kaufman pivoting and apply it to the KKT matrix of `P03.15`. ★ Replace the truss root supports by elastic springs and show how κ(K) grows as the spring stiffness goes to zero. |
| **Time** | unmeasured |

## 1. The engineering problem

Almost every computation in this phase ends in a linear system. The implicit Euler step of `P03.06` solved one per Newton iteration; the finite-element, finite-difference and optimisation lessons ahead solve one per time step, per load case or per design iteration. `P01.13` taught elimination by hand. This lesson turns it into an algorithm with three properties an engineer must know: its **cost** (2n³/3 flops to factor, 2n² to solve, so factor once and reuse), its **stability** (partial pivoting keeps the multipliers bounded; without it a pivot of 10⁻¹⁷ destroys the answer), and its **accuracy** (the residual is always tiny; the error is tiny only when the matrix is well conditioned).

The vehicle is a wing-spar truss: ten nodes, seventeen aluminium bars, a 2 m cantilever 0.3 m deep, 4 kN at the tip. Its stiffness matrix is symmetric and, once the root is held, positive definite, so Cholesky applies at half the cost of LU and without any pivoting, and the factorisation serves every load case that follows. A second vehicle is a heated rod, whose finite-difference matrix is tridiagonal and is solved in O(n) by the Thomas algorithm, with the exact solution available to measure the discretisation order. A third is the Hilbert matrix, which is small, symmetric, positive definite and so ill-conditioned that at n = 10 the computed solution has five correct digits out of sixteen while its residual is at round-off.

Success criteria, written first: every factorisation reconstructs its matrix to 10⁻¹³ relative and matches `scipy.linalg` to 10⁻¹²; the truss solution satisfies equilibrium to 10⁻¹⁰ and the hand statics of the root bay to 10⁻⁹; the rod shows observed order 2.00 ± 0.1; the Hilbert forward error follows κε with slope within 0.75–1.15 while the backward error stays below 10⁻¹⁵.

## 2. Understand and derive

### 2.1 Elimination is a factorisation

Gaussian elimination on A subtracts multiples l_ik = a_ik/a_kk of row k from the rows below, for k = 1, …, n − 1, leaving an upper-triangular U. Each elimination step is a left-multiplication by a unit lower-triangular matrix, and their product's inverse is unit lower triangular with the multipliers in it: **A = LU** (`@ocw-18-085` Lecture 3, `@ocw-18-335j` Lecture 13). The step k update of the trailing block is a rank-1 update, a_ij ← a_ij − l_ik a_kj for i, j > k, which costs 2(n − k)² flops; summed, **2n³/3**. The solve is two triangular substitutions, Ly = b then Ux = y, each n² flops: **2n²** per right-hand side. Hence the central rule of this lesson: with m load cases the cost is 2n³/3 + 2mn², not m·2n³/3, and *never* form A⁻¹ to multiply by it (A⁻¹ costs 2n³, three times the factorisation, and A⁻¹b costs the same 2n² as the substitutions while being less accurate).

**Pivoting.** The multiplier l_ik = a_ik/a_kk is undefined when the pivot a_kk is zero ([[0, 1], [1, 0]] has no LU factorisation) and catastrophic when it is tiny: for A = [[10⁻¹⁷, 1], [1, 1]], b = (1, 2), the multiplier is 10¹⁷, the update 1 − 10¹⁷ rounds to −10¹⁷, and the computed solution is (0, 1) instead of (1, 1). **Partial pivoting** swaps the row with the largest |a_ik| in column k into the pivot position before eliminating, so every multiplier satisfies |l_ik| ≤ 1 and the factorisation is **PA = LU** with P a permutation (stored as the index array `perm`: row i of PA is row `perm[i]` of A). The test checks both failures and the cure.

**Growth.** With |l_ik| ≤ 1, each update can at most double the largest entry, so after n − 1 steps the entries of U are bounded by 2ⁿ⁻¹ max|A|; the ratio ρ = max|U|/max|A| is the **growth factor**. Wilkinson's matrix (1 on the diagonal, −1 below, 1 in the last column) reaches the bound exactly (ρ = 32 at n = 6, 5.4 × 10⁸ at n = 30), yet for random matrices ρ is about 1 (0.96 in the run), and in fifty years of practice no naturally occurring matrix has been found to grow badly. The backward error of LU with partial pivoting is bounded by about ρ n ε, which is why the Hilbert residuals in part 4 are at 10⁻¹⁷ regardless of how ill-conditioned the matrix is.

### 2.2 Cholesky: the pivots are the test

If A is symmetric, elimination preserves the symmetry of the trailing block, and A = LDLᵀ with D the pivots. If moreover every pivot is positive, write L̃ = LD^{1/2} and **A = L̃L̃ᵀ**: the Cholesky factorisation, computed column by column from

**L_jj = √(A_jj − Σ_{k<j} L_jk²),  L_ij = (A_ij − Σ_{k<j} L_ik L_jk)/L_jj  (i > j),**   (1)

at n³/3 flops, half of LU, with no pivoting and half the storage. The quantity under the square root is the j-th pivot. By `@ocw-18-085` Lecture 7, A is positive definite (xᵀAx > 0 for all x ≠ 0) if and only if all its pivots are positive, so **a non-positive pivot in (1) is a proof that A is not positive definite**, and `cholesky` raises exactly there. The test demonstrates it three ways: [[1, 2], [2, 1]] (eigenvalues 3 and −1) fails at column 1; the singular [[1, 1], [1, 1]] fails with a zero pivot; the full truss stiffness matrix, before supports are applied, fails at column 19 because its three rigid-body modes make it singular; and for a nearly singular SPD matrix the last pivot equals the small eigenvalue to four digits. Stiffness, conductance, mass, covariance and normal-equation matrices are SPD by construction, so Cholesky is both the solver of choice and a free sanity check that the assembly is right.

### 2.3 Tridiagonal systems in O(n)

Central differences on a 1-D problem give a matrix with three nonzero diagonals. Elimination then touches one entry per row: the **Thomas algorithm** sweeps forward, b_i ← b_i − (a_i/b_{i−1}) c_{i−1}, d_i ← d_i − (a_i/b_{i−1}) d_{i−1}, and back, x_i = (d_i − c_i x_{i+1})/b_i, in **8n** flops, against 2n³/3 for a dense solve. It is LU without pivoting, so it is safe only when the matrix is diagonally dominant or SPD, which finite-difference matrices are. For the heated rod with n = 600 interior points the run shows 0.80 ms against 225 ms for our dense LU: a factor 280, growing as n².

### 2.4 The truss: assembling a stiffness matrix

A bar between nodes i and j of length L, with direction cosines (c, s), stiffness EA/L and end displacements (u_i, v_i, u_j, v_j), has strain (Δ of the projection of the end displacements onto the bar)/L and force N = (EA/L)[−c, −s, c, s]·u_e. Its contribution to the global stiffness matrix is the 4×4 block

**k_e = (EA/L) [[k, −k], [−k, k]],  k = [c, s]ᵀ[c, s],**   (2)

added into the rows and columns of its four degrees of freedom. The global K is symmetric, and for an unsupported planar truss it has exactly three zero eigenvalues (two translations, one rotation): the test counts them. Fixing the two root nodes (degrees of freedom 0, 1, 10, 11) deletes four rows and columns and leaves a 16×16 SPD matrix with eigenvalues from 6.1 × 10⁴ to 1.4 × 10⁸ N/m, κ₂ = 2 227: a well-conditioned problem, as structural stiffness matrices of this size are. The reactions are R = Ku − f on the fixed degrees of freedom, and they must balance the applied loads (`P01.13` never let you check that; the code does).

The root bay can be solved by statics, which is the independent check. The tip load F = 4 kN at 2 m produces a moment 8 kN·m at the root, carried by the chords as a couple over the depth d = 0.3 m: the upper chord is in tension at **F·2/d = 26.67 kN** (133 MPa in 2 cm²: well below the yield of 2024-T3, about 40 % of it). The lower chord's root member is in compression at F·1.5/d = **20.00 kN** (the moment about the upper node at x = 0.5 m, where the diagonal takes over), and the root diagonal carries the shear, F/sin θ with tan θ = 0.3/0.5: **7.77 kN** in compression. The code reproduces all three to 10⁻⁹ relative, and the tip deflects **20.32 mm**, with 10.63 mm under four 1 kN loads on the upper chord. Problem 4 varies the depth.

### 2.5 The heated rod: a convergence check for free

A rod of length L = 0.5 m, conductivity k = 15 W/m/K, with a volumetric source q₀ sin(πx/L), q₀ = 2 × 10⁵ W/m³, and both ends at T_w = 300 K, obeys −kT'' = q₀ sin(πx/L), whose exact solution is T = T_w + (q₀L²/kπ²) sin(πx/L): a peak rise of 337.7 K. Central differences on n interior points with spacing h = L/(n + 1) give the tridiagonal system (k/h²)(−T_{i−1} + 2T_i − T_{i+1}) = q_i with the wall temperatures moved to the right-hand side at both ends. The truncation error of the second difference is h²T''''/12, so the error should fall as h²: the run gives 2.282, 0.629, 0.165, 0.042, 0.011 K for n = 10, 20, 40, 80, 160, observed order **2.00**. This is the lesson's convergence test, and its purpose is as much to check `thomas_solve` as the discretisation: a wrong solver would not converge at all.

### 2.6 Backward error, forward error, and what the residual can tell you

Let x̂ be the computed solution. Its **residual** r = b − Ax̂ is cheap to compute and defines the **normwise relative backward error** η = ‖r‖/(‖A‖‖x̂‖ + ‖b‖): the smallest relative perturbation of (A, b) for which x̂ is the exact solution. A **backward stable** solver returns η = O(ε) always, and LU with partial pivoting does (η ≤ 3 × 10⁻¹⁷ for every Hilbert matrix in part 4). The **forward error** ‖x̂ − x‖/‖x‖ is what you want and cannot compute without the answer; the two are linked by the condition number κ(A) = ‖A‖‖A⁻¹‖ (`@ocw-18-335j` Lecture 7, `P01.18`):

**‖x̂ − x‖/‖x‖ ≤ κ(A) · η · (1 + ‖b‖/(‖A‖‖x̂‖)) ≈ κ(A) · η.**   (3)

The derivation is three lines: A(x̂ − x) = −r, so ‖x̂ − x‖ ≤ ‖A⁻¹‖‖r‖ = κ‖r‖/‖A‖, and divide by ‖x‖ ≈ ‖x̂‖. The consequence is the lesson's slogan: **a small residual proves the solver worked; it says nothing about the answer unless κ is small.** For the Hilbert matrices, κ₁ grows from 2.8 × 10⁴ at n = 4 to 3.5 × 10¹³ at n = 10, and the forward error from 1.3 × 10⁻¹³ to 5.2 × 10⁻⁵, a slope of 0.94 in log κ against the theoretical 1, always below the bound κε (6 × 10⁻¹² to 8 × 10⁻³). At n = 12 the matrix is singular to working precision and `lu_decompose` refuses it. For the truss, κ₂ = 2 227 and the deflection is good to 12 digits; for a Hilbert-like normal-equations matrix in Problem 5 it is good to none.

**Iterative refinement** uses the factors once more: compute r in extended precision, solve A δ = r with the stored L and U (2n² flops), and set x̂ ← x̂ + δ. With a double-length residual each pass gains about log₁₀(1/(κε)) digits until the working-precision limit; the Hilbert errors fall by a factor 1.5–2.5 after three passes here, not more, because NumPy's `longdouble` is only 80-bit on this platform and the residual itself is barely more accurate than the solution. The point of the technique is less the gain than the cost: n² per pass against n³ to re-solve more carefully.

### 2.7 Where the methods stop being valid

- **Dense storage.** A 10⁴-node finite-element model has a 2 × 10⁴ square K with 4 × 10⁸ entries (3 GB) almost all zero; dense LU costs 5 × 10¹² flops. Banded, sparse-direct (`@ocw-18-335j` Lecture 24) and iterative solvers (`P03.09`) exist for that reason, and the fill-in of sparse LU is their central problem.
- **Partial pivoting** can fail in principle (ρ = 2ⁿ⁻¹) and complete pivoting is the cure at 2× the cost; in practice nobody pays it.
- **Thomas without pivoting** is unsafe for a tridiagonal matrix that is not diagonally dominant (a convection-dominated discretisation with large cell Péclet number, `P03.12`).
- **Cholesky on a nearly singular SPD matrix** (a structure with a very soft mode, extension 3) gives tiny pivots and a large κ; the factorisation succeeds and the answer is garbage. The pivot test is a test of definiteness, not of conditioning.
- **The truss model** itself: pin joints, no bending, no buckling of the compressed members (the lower chord's 20 kN compression in a solid 2 cm² bar 0.5 m long exceeds that bar's Euler load of about 9 kN; a real spar uses a tube or a deeper section, and `P05` computes it), small deflections (20 mm on 2 m is fine).

## 3. Build the computational model

Open `starter/main.py`. `spar_truss`, `heated_rod_system`, `heated_rod_exact`, `hilbert`, `wilkinson_matrix`, `norm1`, `convergence_order` and `main` are given.

### Step 1: substitutions and LU

```python
def lu_decompose(A, pivot=True):
    A = np.array(A, dtype=float); n = len(A); perm = np.arange(n)
    for k in range(n):
        if pivot:
            p = k + np.argmax(np.abs(A[k:, k]))
            A[[k, p]] = A[[p, k]]; perm[[k, p]] = perm[[p, k]]
        A[k + 1:, k] /= A[k, k]                                  # multipliers, stored in place
        A[k + 1:, k + 1:] -= np.outer(A[k + 1:, k], A[k, k + 1:])   # rank-1 update
    return np.tril(A, -1) + np.eye(n), np.triu(A), perm
```

Add the two `ValueError` checks (zero pivot without pivoting; pivot below 10⁻¹⁴ max|A| with it). `forward_substitution` and `back_substitution` are the loops of 2.1; `lu_solve` permutes b then calls both; `solve` handles several right-hand sides. `determinant_from_lu` multiplies the diagonal of U by the sign of `perm` (count its cycles: sign = (−1)^{n − cycles}); `inverse_from_lu` solves for each column of I; `growth_factor` factors and takes the ratio of maxima.

### Step 2: Cholesky and Thomas

`cholesky` is (1) with the symmetry check and the pivot test; `cholesky_solve` is two substitutions with L and Lᵀ. `thomas_solve` is the two sweeps of 2.3 on copies of the four arrays.

### Step 3: errors

`condition_number_1` uses `norm1` and `inverse_from_lu`; `backward_error` is η; `iterative_refinement` loops r (in `np.longdouble`), δ = `lu_solve(... r)`, x += δ, and records ‖r‖₁.

### Step 4: the truss

`truss_stiffness` loops over elements and adds (2) into K; `apply_supports` returns the sub-matrix and the free indices; `member_forces` evaluates N_e; `solve_truss` assembles, reduces, factors with `cholesky` once, solves each column of the loads, and returns displacements and reactions Ku − f.

### Step 5: the rod study

`heated_rod_errors` builds the system for each n, solves with `thomas_solve`, and returns the max-norm error against `heated_rod_exact`.

### Verify

Run `python -m pytest curriculum/phase-03/lesson-08-direct-linear-solvers`. The tests:

1. **Substitutions** equal `solve_triangular`; unit-diagonal mode; a zero diagonal raises.
2. **LU** on n = 1, 3, 6, 12: A[perm] = LU to 10⁻¹³, unit L, |multipliers| ≤ 1, same L and U as `scipy.linalg.lu`, solution to 10⁻¹⁰, backward error below 10⁻¹⁵.
3. **Pivoting**: [[0, 1], [1, 0]] raises without pivoting and solves with it; the 10⁻¹⁷ pivot gives an error above 0.1 unpivoted and 10⁻¹² pivoted; singular and non-square inputs raise; several right-hand sides.
4. **Determinant, inverse, κ₁** against NumPy; a row swap flips the determinant's sign.
5. **Growth** exactly 2ⁿ⁻¹ on Wilkinson for n = 4, 6, 10; below 10 on a random 50×50.
6. **Cholesky** on SPD matrices of size 1, 4, 9: reconstruction, lower-triangular with positive diagonal, equal to SciPy; refuses unsymmetric, indefinite and singular; the last pivot equals the small eigenvalue.
7. **Thomas** equals a dense solve; a zero pivot raises.
8. **Hilbert**: η < 10⁻¹⁵ and forward error < κε for n = 4–10, slope of log(error) against log κ in (0.75, 1.15); refinement does not worsen the error and reduces the residual.
9. **Truss**: 10 nodes, 17 members; K symmetric with exactly three zero modes; 16 free degrees of freedom; SPD when supported and Cholesky refuses the unsupported K; a single bar gives EA/L on the x degrees of freedom; equilibrium of reactions to 10⁻¹⁰; root chord and diagonal forces by statics to 10⁻⁹; equal to `scipy.linalg.solve`; superposition of load cases.
10. **Rod**: errors decrease, order 2.00 ± 0.1, error below 0.2 K at n = 40, the exact peak, and Cholesky on the SPD tridiagonal matrix agrees with Thomas.
11. `main` writes the figure.

### Compare

`scipy.linalg.lu_factor`/`lu_solve`, `cholesky`/`cho_solve` and `solve` call LAPACK: the same algorithms, blocked for cache and vectorised, in compiled code. Our pure-NumPy LU takes 2.3, 9.3 and 62 ms at n = 100, 200, 400 against 0.32, 1.03 and 2.81 ms for `lu_factor`; the observed time exponent of ours is 2.37 because Python overhead per step hides the n³ until n is larger, while LAPACK's blocked algorithm runs near the machine's peak. The factors agree to 10⁻¹² and the solutions to 10⁻¹⁰, because both pivot identically. We call the library in practice for speed, for the banded, packed and sparse variants (`scipy.linalg.solve_banded`, `scipy.sparse.linalg.splu`), and for the condition-number estimator that comes free with `solve(..., check_finite)` and `lu_factor`'s error flags; we write our own once so that "LU with partial pivoting is backward stable but the forward error is κ times larger" is something we have measured. `scipy.linalg.solve(A, b, assume_a="pos")` is Cholesky and should be used whenever the matrix is known to be SPD: it is twice as fast and raises on an indefinite matrix, which, as 2.2 argues, is a feature.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Cost.** Time `lu_decompose` and `scipy.linalg.lu_factor` at n = 100, 200, 400 and fit the exponent of time against n for each. Predict 3 for both; expect the pure-Python one to be lower at these sizes. Time `thomas_solve` against `solve` on the rod at n = 600 and predict the ratio from 8n against 2n³/3.
2. **Conditioning.** For Hilbert n = 4, 6, 8, 10 with x = ones, record κ₁, the backward error, the forward error and κε; fit the slope of forward error against κ. Predict: backward error at round-off throughout, forward error growing with slope 1, staying below κε.
3. **Convergence.** Rod errors for n = 10 … 160; predict order 2.
4. **Load cases.** Solve the truss for the tip load, the distributed load and their sum from one Cholesky factorisation; verify superposition and count flops (n³/3 once, 2n² per case).

Raw data format: CSV tables (`n, t_ours_ms, t_scipy_ms`; `n, kappa, backward, forward, kappa_eps`; `n, h, max_err_K`; `case, tip_deflection_mm, max_member_force_kN`), each with a header naming the fixed parameters.

What you should find: exponents 2.37 (ours) and about 1.5 (SciPy, whose cost at these sizes is dominated by overhead); Thomas 0.80 ms against 225 ms dense (predicted ratio n²/12 = 30 000 in flops; the measured 280 reflects the per-row Python loop in Thomas against vectorised rank-1 updates in LU); the Hilbert table of 2.6 with slope 0.94; rod order 2.00; tip deflections 20.32 mm, 10.63 mm and 30.95 mm, the last the sum of the first two to 10⁻¹⁰, with maximum member forces 26.67 kN and 16.67 kN for the first two cases (the distributed load has a root moment of 1 kN × (0.5 + 1.0 + 1.5 + 2.0) m = 5 kN·m, so 5/0.3 = 16.67 kN).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| A[perm] − LU, random 6×6 | 0 | 2.2 × 10⁻¹⁶ | — | round-off |
| growth factor, Wilkinson n = 6 / 30 | 32 / 5.369 × 10⁸ | 32 / 5.369 × 10⁸ | 0 | exact |
| zero modes of unsupported K | 3 | 3 | 0 | exact |
| κ₂ of supported K | — | 2 227 | — | well conditioned |
| root upper chord (kN) | 26.67 | 26.67 | 10⁻⁹ | statics is exact for a determinate bay |
| root lower chord (kN) | −20.00 | −20.00 | 10⁻⁹ | same |
| root diagonal (kN) | −7.77 | −7.77 | 10⁻⁹ | same |
| vertical reactions (kN) | 4.00 | 4.000 | 10⁻¹⁰ | equilibrium |
| K − LLᵀ (relative) | 0 | 2.0 × 10⁻¹⁶ | — | round-off |
| u vs `scipy.linalg.solve` (m) | 0 | 8.7 × 10⁻¹⁸ | — | identical algorithm |
| rod order | 2 | 2.00 | 0 | asymptotic |
| rod peak rise (K) | 337.74 | 337.74 − 0.011 at n = 160 | h² | truncation |
| Hilbert n = 10: η / forward error | ≤ ε / ≤ κε = 7.9 × 10⁻³ | 3.1 × 10⁻¹⁷ / 5.2 × 10⁻⁵ | — | bound (3) holds with margin 150 |

Everything numerical here is checked against an exact statement: a factorisation identity, a count of zero modes, statics, an exact solution, a theorem. The one row that is a *measurement* rather than a check is κ₂ = 2 227 for the truss, and it is the row that answers the engineering question "how many digits are real": with η ≈ 10⁻¹⁶ and κ ≈ 2 × 10³, about 12 of the 16 digits of the 20.32 mm deflection are determined by the arithmetic. The remaining question is whether the *model* deserves them: E for an aluminium alloy is known to about 2 %, bar areas to 1 %, and pin joints do not exist, so the deflection is engineering-accurate to perhaps 5 %, and the solver's 12 digits are a luxury that costs nothing. The Hilbert rows make the opposite point: with κ = 3.5 × 10¹³ the solver, still backward stable, returns five digits, and no amount of care in the elimination improves that; only a better-conditioned formulation does (Problem 5 shows QR doing exactly that for a least-squares problem).

Improvements worth making: banded storage for the rod (extension 1), which turns the 225 ms into the 0.8 ms of Thomas without a special-purpose routine; `assume_a="pos"` wherever a matrix is SPD; refinement with a true double-double residual. Improvement worth *not* making: computing A⁻¹ "to have it available". It costs three factorisations, loses accuracy, and nobody needs the inverse itself.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of the LU cost, the Cholesky pivot test, (2) and (3), and the statics of the root bay; `direct_solvers.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the truss tip deflection with the number of digits you trust, separating arithmetic from model uncertainty, and the cost in flops of adding ten more load cases.

### Adversarial self-checks

- **What if the load doubles?** Everything doubles: the model is linear, and `solve_truss` with 2f returns 2u to round-off (the superposition test). Nothing in the arithmetic changes; the lower chord's 40 kN compression is what changes the *model's* validity (buckling).
- **What if the depth goes to zero?** The chord forces F·2/d → ∞ and the truss becomes a mechanism; κ(K) grows and the three zero modes become four. Problem 4 measures κ against d.
- **Which assumption fails first in the real world?** Pin joints: a riveted or welded truss carries bending at the joints, which changes member stresses by 5–20 % (secondary stresses) but the deflection by much less; then buckling of the compressed chord.
- **What measurement could falsify the model?** A dial gauge at the tip: 20.3 mm predicted; 25 mm measured would indicate joint slip or a lower E; strain gauges on the root chords should read 133 and −100 MPa (± 5 % for secondary bending).
- **Does the result survive a change of solver?** LU, Cholesky, `scipy.linalg.solve` and `assume_a="pos"` agree to 10⁻¹⁷ m. Thomas and Cholesky agree on the rod to 10⁻¹⁰. The Hilbert result does *not* survive a change of formulation, which is the point: it is the conditioning, not the solver.
- **Which decision changes if the uncertainty doubles?** If E is known to 4 % instead of 2 %, the deflection is known to 4 %; the solver's contribution (10⁻¹²) is irrelevant either way, and the decision to spend n³/3 rather than 2n³/3 flops changes nothing. If κ doubled (a softer support) the forward error would still be 10⁻¹². Only at κ ≈ 10¹² (Hilbert n = 9) does the arithmetic start to matter, and the decision then is to reformulate, not to refine.

## Check yourself

1. Write the cost in flops of LU factorisation, one solve from the factors, forming the inverse, and Cholesky, for an n×n matrix.
2. Why is partial pivoting needed? Give a 2×2 matrix with no LU factorisation and one where elimination without pivoting loses all accuracy.
3. What is the growth factor, what is its bound under partial pivoting, and which matrix attains it?
4. State the Cholesky pivot test: what does a non-positive pivot prove?
5. Define the normwise relative backward error and write the inequality linking it to the forward error.
6. A solver returns a residual of 10⁻¹⁶ for a system with κ = 10¹⁰. How many digits of the solution are reliable?
7. Why is the Thomas algorithm O(n), and when is it unsafe?
8. For the spar truss, how is the upper chord force at the root obtained without a computer?
