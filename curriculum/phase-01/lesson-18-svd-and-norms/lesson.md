# P01.18  SVD & Matrix Norms

> Every matrix stretches some direction most and some direction least. The ratio of the two is how much it can lie to you.

| | |
|---|---|
| **Learning outcomes** | 1. Define vector and matrix norms, prove the consistency and submultiplicative inequalities, and show that the induced 2-norm is the largest singular value and the Frobenius norm the root-sum-square of all of them. 2. Derive the singular value decomposition from the eigen-decomposition of AᵀA, and the error-amplification bound ‖δu‖/‖u‖ ≤ cond(K) ‖δf‖/‖f‖ with the load patterns that attain it. 3. Implement and test, from scratch, the four norms, a one-sided Jacobi SVD, the condition number, the pseudoinverse, the null space and the best rank-k approximation, against `numpy.linalg.svd`, `norm`, `cond` and `pinv`. 4. Measure how the condition number of a bar's stiffness matrix grows with element count and with a damaged element, decide whether a 1 % load error matters, and detect a mechanism from the singular values. |
| **Prerequisites** | `P01.13` (solving K u = f), `P01.14` (rank, null space), `P01.16` (eigenvalues of symmetric matrices, power iteration, Jacobi rotations), `P01.17` (least squares, orthonormal columns), `P00.03`, `P00.07`. Skill check: you can find the eigenvalues of a 2×2 symmetric matrix and state what the null space of a matrix is. |
| **Engineering question** | A 1 m aluminium bar, modelled as ten axial elements and loaded with 1 kN known to 1 %, deflects 0.143 mm: how wrong can the computed displacements be, which load pattern makes them that wrong, how much worse does a damaged element make it, and how does the matrix tell you the structure has become a mechanism? |
| **Assumptions and limits** | Linear elastic elements, small displacements, axial loads only, node 0 clamped; the condition number bounds the *worst case* over all load perturbations, which for a bar is an alternating-sign pattern no real load applies, so the bound is pessimistic by two orders of magnitude for typical loads (part 5). The one-sided Jacobi SVD needs m ≥ n and is slower than the library for large matrices; the pseudoinverse's rank decision depends on a threshold you must choose. |
| **Analytical work** | Norm axioms; ‖Ax‖ ≤ ‖A‖‖x‖; ‖A‖₂ = σ₁ and ‖A‖_F² = Σσ_i² from A = UΣVᵀ; the SVD from AᵀA; cond(K) and the perturbation bound with its attaining vectors; exact eigenvalues of the uniform chain and cond ≈ 16n²/π²; the one-sided Jacobi rotation angle. Hand calculation of all four norms and the condition number of a 2×2 matrix. |
| **Computational deliverable** | `starter/main.py`: `vector_norm`, `matrix_norm_1`, `matrix_norm_inf`, `frobenius_norm`, `matrix_norm_2_power`, `svd_jacobi`, `condition_number`, `pseudoinverse`, `low_rank_approximation`, `null_space`, `element_stiffness`, `chain_stiffness_matrix`, `uniform_chain_condition_analytic`, `tip_deflection_analytic`, `error_amplification`, `worst_case_loads`, `mechanism_modes`. Tests prove: norms match NumPy and satisfy the axioms, consistency and submultiplicativity, with the induced norms attained; the SVD reconstructs A with orthonormal factors and matches `numpy.linalg.svd`, including rank-deficient inputs; cond(K) equals the exact eigenvalue ratio and grows as n²; the perturbation bound holds for random loads and is attained for the worst-case pair; a failed element produces exactly one null-space mode of the right shape; the pseudoinverse matches `pinv` in the tall, wide, invertible and singular cases; low-rank approximations satisfy Eckart–Young. |
| **Practical activity** | **computational experiment**: sweep the element count and the stiffness of one damaged element, measure cond₂(K) and the actual error amplification for random and worst-case loads, and compare with the analytic bound. |
| **Equipment** | None |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `stiffness_conditioning.png`, the amplification table (random versus worst case), and a short note on whether the 1 % load uncertainty matters for this bar. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4; pass when the norm inequalities and the perturbation bound are derived with the attaining vectors identified, all public tests pass, and the note distinguishes the bound from the typical case with numbers. |
| **References** | `@mit-ocw-18-06-linear-algebra` Lecture 29 "Singular value decomposition", Lecture 33 "Left and right inverses; pseudoinverse" (required-free); Lecture 25 "Symmetric matrices and positive definiteness" (prerequisite, for AᵀA). |
| **Extensions** | ★ Build a bar whose elements differ in stiffness by 10⁴ *without* a soft element in series (alternate stiff and soft elements in parallel pairs) and show that diagonal scaling then removes most of the condition number, unlike the damaged bar. ★ Implement the Golub–Kahan bidiagonalisation step and compare its cost with one-sided Jacobi at n = 200. ★ Apply `low_rank_approximation` to a matrix of pressure-distribution snapshots and show that two modes capture 99.9 % of the energy (problem 6 starts this). |
| **Time** | unmeasured |

## 1. The engineering problem

A finite-element model of a structure ends in a linear system K **u** = **f**: K is the **stiffness matrix** (N/m), **f** the vector of nodal loads (N) and **u** the nodal displacements (m). Lesson 13 solved such systems. This lesson asks a question that no solver answers on its own: *how much should you trust the solution?* Loads are never known exactly. A 1 kN test load is applied with a jack whose calibration is good to 1 %; a flight load comes from an aerodynamic model good to 5 %. If a 1 % error in **f** can become a 100 % error in **u**, the model is useless for that load case, however accurate the solver.

The quantity that answers the question is the **condition number** of K, the ratio of its largest to its smallest singular value. It is a property of the matrix, not of the solver or the data, and it bounds the amplification of relative error from **f** to **u**. The same number says whether a least-squares fit is trustworthy (Lesson 17 used it), whether a control system is close to uncontrollable, and whether a structure is close to being a mechanism: as an element loses stiffness, the smallest singular value goes to zero and the condition number to infinity, and the corresponding singular vector *is* the mechanism's shape.

The concrete problem: an aluminium bar (E = 70 GPa, 10 mm × 10 mm, 1 m long) clamped at one end and modelled as 10 axial elements of stiffness EA/L_e = 7 × 10⁷ N/m each. We want the norms and condition number of its 10×10 stiffness matrix, the worst-case and typical amplification of a 1 % load error, how both change when element 5 is damaged to 0.1 %, 0.01 %, … of its stiffness, and what the matrix looks like when the element has failed entirely. Success criteria: our from-scratch SVD agrees with `numpy.linalg.svd` to 10⁻¹⁰ and reconstructs K to rounding; the condition number equals the exact eigenvalue formula for the uniform chain to 10⁻⁹; the perturbation bound is attained to 10⁻⁸ by the predicted load pair; a failed element gives exactly one null-space mode, and its shape is the one a structural engineer would draw.

## 2. Understand and derive

### 2.1 Vector norms

A **norm** is a function ‖·‖ that assigns a length to a vector and satisfies three axioms: ‖**x**‖ ≥ 0 with equality only for **x** = 0; ‖c**x**‖ = |c| ‖**x**‖; and the triangle inequality ‖**x** + **y**‖ ≤ ‖**x**‖ + ‖**y**‖. Three are in daily use:

$$
\|\mathbf x\|_1 = \sum_i |x_i|, \qquad \|\mathbf x\|_2 = \Big(\sum_i x_i^2\Big)^{1/2}, \qquad \|\mathbf x\|_\infty = \max_i |x_i|.
$$

For a load vector, ‖**f**‖₁ is the total absolute load, ‖**f**‖∞ the largest single nodal load, ‖**f**‖₂ the usual Euclidean length. They are never far apart: ‖**x**‖∞ ≤ ‖**x**‖₂ ≤ ‖**x**‖₁ ≤ √n ‖**x**‖₂ ≤ n ‖**x**‖∞. The tests check the axioms and this ordering on random vectors. Which norm you use changes the numbers by at most a factor √n, never the conclusions.

### 2.2 Matrix norms: how much can A stretch a vector?

The natural size of a matrix is its largest stretch factor. For a vector norm ‖·‖_p, the **induced** (operator) norm is

$$
\|A\|_p = \max_{\mathbf x \ne 0} \frac{\|A\mathbf x\|_p}{\|\mathbf x\|_p}.
$$

Two consequences follow straight from the definition and are the two inequalities the whole lesson rests on. **Consistency:** ‖A**x**‖ ≤ ‖A‖ ‖**x**‖ for every **x** (the maximum is at least the value at **x**). **Submultiplicativity:** ‖AB‖ ≤ ‖A‖ ‖B‖, because ‖AB**x**‖ ≤ ‖A‖ ‖B**x**‖ ≤ ‖A‖ ‖B‖ ‖**x**‖.

For p = 1 and p = ∞ the maximum can be read off the entries. ‖A**x**‖₁ = Σ_i |Σ_j a_ij x_j| ≤ Σ_j |x_j| Σ_i |a_ij| ≤ (max_j Σ_i |a_ij|) ‖**x**‖₁, and the bound is attained by **x** = **e**_j for the column with the largest absolute sum. So

$$
\|A\|_1 = \max_j \sum_i |a_{ij}| \ \text{(largest column sum)}, \qquad \|A\|_\infty = \max_i \sum_j |a_{ij}| \ \text{(largest row sum)},
$$

the second attained by the sign vector of the largest row. Transposing swaps the two. The tests check attainment, not just the inequality.

The **Frobenius norm** ‖A‖_F = (Σ_ij a_ij²)^{1/2} treats the matrix as a long vector. It is not induced by any vector norm, but it is submultiplicative and easy, and it dominates the 2-norm: ‖A‖₂ ≤ ‖A‖_F ≤ √n ‖A‖₂.

The 2-norm has no entry formula. ‖A**x**‖₂² = **x**ᵀAᵀA**x**, and the maximum of this Rayleigh quotient over unit **x** is the largest eigenvalue of the symmetric matrix AᵀA (Lesson 16), attained at its eigenvector. Hence ‖A‖₂ = √λ_max(AᵀA), and power iteration on AᵀA computes it; that is `matrix_norm_2_power`. The square roots of all the eigenvalues of AᵀA have a name.

### 2.3 The singular value decomposition

AᵀA is symmetric positive semi-definite, so it has an orthonormal eigenbasis **v**₁, …, **v**_n with eigenvalues λ_i ≥ 0 (`@mit-ocw-18-06-linear-algebra` Lecture 29). Define the **singular values** σ_i = √λ_i, ordered σ₁ ≥ σ₂ ≥ … ≥ 0, and for each σ_i > 0 the unit vector **u**_i = A**v**_i/σ_i. These **u**_i are orthonormal: **u**_iᵀ**u**_j = **v**_iᵀAᵀA**v**_j/(σ_iσ_j) = λ_j **v**_iᵀ**v**_j/(σ_iσ_j) = δ_ij. Collecting them,

$$
A = U\,\Sigma\,V^{\mathsf T} = \sum_{i=1}^{r} \sigma_i\,\mathbf u_i \mathbf v_i^{\mathsf T},
$$

with U (m×n) and V (n×n) having orthonormal columns and Σ = diag(σ_i). Geometrically: A maps the unit sphere to an ellipsoid whose semi-axes are σ_i **u**_i; **v**_i is the input direction that lands on the i-th axis. The rank r of A is the number of non-zero singular values; the **v**_i with σ_i = 0 span the null space, the **u**_i with σ_i > 0 span the column space.

Everything in 2.2 now reads off Σ: ‖A‖₂ = σ₁ (the largest semi-axis); ‖A‖_F² = tr(AᵀA) = Σλ_i = Σσ_i² (the Frobenius norm is unchanged by orthogonal U and V); and the smallest stretch of a full-rank square matrix is σ_n. For a symmetric positive definite K, such as a stiffness matrix, AᵀA = K² has eigenvalues λ_i², so the singular values of K are its eigenvalues and U = V: the SVD is the eigen-decomposition. The tests check this on the bar.

Computing the SVD *through* AᵀA squares the condition number (Lesson 17's lesson), so real SVD codes avoid it. **One-sided Jacobi** (Hestenes' method) rotates pairs of *columns* of A directly: for columns p and q with α = ‖**w**_p‖², β = ‖**w**_q‖², γ = **w**_pᵀ**w**_q, the rotation by the angle with t = tan θ = sign(ζ)/(|ζ| + √(1 + ζ²)), ζ = (β − α)/(2γ), makes the two columns orthogonal. Sweep over all pairs until every pair is orthogonal to rounding; then the column norms are the σ_i, the normalised columns are the **u**_i, and the accumulated rotations are V. It is the same rotation as Lesson 16's Jacobi applied implicitly to AᵀA, without ever forming it. Each rotation touches two columns of length m, so a sweep of n(n − 1)/2 rotations costs about 3mn² flops; convergence is quadratic after the first few sweeps, and for the 10×10 stiffness matrix seven sweeps rotate something and the eighth finds every pair already orthogonal and stops.

### 2.4 The condition number and the perturbation bound

Solve K**u** = **f**, and let the load carry an error δ**f**, so that K(**u** + δ**u**) = **f** + δ**f**. Then δ**u** = K⁻¹δ**f** and, by consistency, ‖δ**u**‖₂ ≤ ‖K⁻¹‖₂ ‖δ**f**‖₂. Also ‖**f**‖₂ = ‖K**u**‖₂ ≤ ‖K‖₂ ‖**u**‖₂, so 1/‖**u**‖₂ ≤ ‖K‖₂/‖**f**‖₂. Multiplying,

$$
\frac{\|\delta\mathbf u\|_2}{\|\mathbf u\|_2} \le \underbrace{\|K\|_2\,\|K^{-1}\|_2}_{\kappa_2(K)}\;\frac{\|\delta\mathbf f\|_2}{\|\mathbf f\|_2}, \qquad \kappa_2(K) = \frac{\sigma_1}{\sigma_n},
$$

since ‖K⁻¹‖₂ is the largest singular value of K⁻¹, which is 1/σ_n. κ ≥ 1 always (σ₁ ≥ σ_n), and κ = 1 exactly for an orthogonal matrix, which stretches nothing. The bound is *attained*: take **f** = **u**₁ (the direction K stretches most, so **u** = **v**₁/σ₁ is as short as a unit load can make it) and δ**f** along **u**_n (the direction K⁻¹ stretches most, so δ**u** = δ**f**-size × **v**_n/σ_n is as long as it can be). The reverse pairing, **f** along **u**_n and δ**f** along **u**₁, gives the *smallest* possible amplification, 1/κ. Both are tested.

In floating point the data error is at least the rounding of the entries, ε_mach ≈ 2 × 10⁻¹⁶, so a solve loses about log₁₀ κ digits in the worst case: κ = 175 costs two digits at most, κ = 10¹⁶ costs all of them. Which κ you quote matters: κ₁ and κ∞ use the entry formulas and differ from κ₂ by at most a factor n; `numpy.linalg.cond` defaults to κ₂.

### 2.5 The stiffness matrix of a bar and its exact condition number

```
   clamp ║──k₁──●──k₂──●──k₃──●── … ──k_n──●──▶ f_n
         ║      u₁     u₂     u₃          u_n
```

Element i joins nodes i − 1 and i with stiffness k_i = EA_i/L_i (N/m: Pa × m² / m). Equilibrium of node i (free body: the two springs pulling on it and the applied load f_i) gives k_i(u_i − u_{i−1}) − k_{i+1}(u_{i+1} − u_i) = f_i, with u₀ = 0 at the clamp and no k_{n+1} beyond the last node. In matrix form K**u** = **f** with K tridiagonal: K_ii = k_i + k_{i+1}, K_{i,i+1} = −k_{i+1}. K is symmetric and, with the clamp, positive definite. Two sanity checks the tests use: a rigid shift **u** = (1, 1, …, 1) loads only the clamp, K**1** = (k₁, 0, …, 0); and a tip load F gives tip displacement F Σ 1/k_i (springs in series), which for the uniform bar is FL/(EA) = 1000 × 1/(70 × 10⁹ × 10⁻⁴) = 0.1429 mm.

For n identical elements the eigenvalues are known in closed form, λ_j = 2k[1 − cos((2j − 1)π/(2n + 1))], j = 1..n. (Check n = 1: λ = 2k[1 − cos(π/3)] = k ✓.) Hence

$$
\kappa_2(K) = \frac{1 - \cos\frac{(2n-1)\pi}{2n+1}}{1 - \cos\frac{\pi}{2n+1}} \approx \frac{2}{\tfrac12(\pi/(2n+1))^2} = \frac{4(2n+1)^2}{\pi^2} \approx \frac{16 n^2}{\pi^2}.
$$

The condition number of a finite-element stiffness matrix grows as the *square* of the number of elements along a load path (for 3D solids, as h⁻² in the mesh size): a refined mesh is a worse-conditioned system, which is why large FE solvers precondition. For n = 10: exact 175.09, asymptotic 162.1; the run reproduces 175.087 from the SVD, from `numpy.linalg.cond` and from the formula.

A damaged element with stiffness rk, r ≪ 1, makes the structure nearly a mechanism: the nodes beyond it can move almost freely, σ_n drops to about rk times a geometric factor, and κ grows as 1/r (part 5 measures the exponent as 1.00). At r = 0, σ_n = 0 exactly, K is singular, and the null vector **v**_n is the mechanism: zero for the nodes before the break, equal for all nodes after it, (0, 0, 0, 0, 1, 1, 1, 1, 1, 1)/√6 for a break in element 5. That is `mechanism_modes`, and it is the structural meaning of Lesson 14's null space.

### 2.6 The pseudoinverse and the best low-rank approximation

Inverting the SVD where it can be inverted gives the **pseudoinverse** A⁺ = V Σ⁺ Uᵀ, with Σ⁺ holding 1/σ_i for σ_i > 0 and 0 otherwise (`@mit-ocw-18-06-linear-algebra` Lecture 33). For a tall full-rank A it equals (AᵀA)⁻¹Aᵀ and returns the least-squares solution (Lesson 17); for a wide full-rank A it equals Aᵀ(AAᵀ)⁻¹ and returns the solution of A**x** = **b** of smallest norm (any other solution differs by a null-space vector, which is orthogonal and only adds length); for an invertible square A it is A⁻¹; and for a singular A it inverts on the column space and ignores the null space, which is the only sensible thing to do with a mechanism: a load along the mechanism has no equilibrium solution, and A⁺ silently drops it. The threshold that decides "σ_i > 0" is a modelling decision, not a numerical detail.

Truncating the sum at k terms gives A_k = Σ_{i≤k} σ_i**u**_i**v**_iᵀ, and the **Eckart–Young theorem** says this is the closest rank-k matrix to A in both the 2-norm (error σ_{k+1}) and the Frobenius norm (error √(Σ_{i>k} σ_i²)). This is the mathematics of data compression and of proper orthogonal decomposition in fluid dynamics: a matrix of pressure snapshots whose singular values fall off fast is a flow with few modes. The test checks both error formulas and that a competing rank-k matrix does worse.

### 2.7 Where the model stops being valid

The bound of 2.4 is a worst case over *all* δ**f**; for the bar the worst-case pattern alternates in sign node by node (part 5 prints it), which no physical load does, so the bound overstates typical errors by two orders of magnitude. The linear model fails for large displacements and for buckling, where the stiffness itself depends on the load. One-sided Jacobi is fine up to a few hundred columns and then loses to the bidiagonalisation methods in the library by a factor of several. The pseudoinverse threshold `rcond` is a guess about which singular values are "really" zero; on a damaged-but-intact element (r = 10⁻⁵) the default threshold keeps it, which is right, but a threshold of 10⁻⁴ would declare the bar a mechanism, which is wrong.

## 3. Build the computational model

Open `starter/main.py`. The bar constants are given.

### Step 1: norms

`vector_norm(x, p)` implements the three formulas and raises `ValueError` for any other p. `matrix_norm_1` and `matrix_norm_inf` are a column-sum and a row-sum with `np.abs`; `frobenius_norm` is `np.sqrt(np.sum(A * A))`. `matrix_norm_2_power` runs power iteration on AᵀA (ones as the start, stop when the Rayleigh quotient changes by less than `tol` relative) and returns the square root.

### Step 2: the SVD

```python
for sweep in range(max_sweeps):
    rotated = False
    for p in range(n - 1):
        for q in range(p + 1, n):
            alpha, beta, gamma = W[:, p] @ W[:, p], W[:, q] @ W[:, q], W[:, p] @ W[:, q]
            if abs(gamma) <= tol * np.sqrt(alpha * beta): continue
            rotated = True
            zeta = (beta - alpha) / (2 * gamma)
            t = np.sign(zeta) / (abs(zeta) + np.sqrt(1 + zeta**2))
            c = 1 / np.sqrt(1 + t**2); s = c * t
            W[:, p], W[:, q] = c * W[:, p] - s * W[:, q], s * W[:, p] + c * W[:, q]   # same rotation on V
    if not rotated: break
```

Afterwards σ_i = ‖**w**_i‖, sorted descending with the columns of W and V permuted alike; **u**_i = **w**_i/σ_i. For a zero singular value, set σ_i = 0 and complete the column of U by Gram–Schmidt of a coordinate vector against the columns found (the test on a rank-2 matrix checks UᵀU = I). Raise `ValueError` for m < n (transpose first) and `RuntimeError` if the sweeps run out.

### Step 3: condition number, pseudoinverse, null space, low rank

`condition_number` is σ₁/σ_n (inf if σ_n = 0). `pseudoinverse(A, rcond)` transposes a wide input, inverts the singular values above `rcond * σ₁`, and forms V Σ⁺ Uᵀ. `null_space` collects the right singular vectors with σ ≤ rtol σ₁ and, for a wide matrix, completes the orthogonal complement of the row space by Gram–Schmidt. `low_rank_approximation(A, k)` returns `(U[:, :k] * s[:k]) @ Vt[:k]`.

### Step 4: the bar

`element_stiffness` is EA/L with validation. `chain_stiffness_matrix(k)` assembles the tridiagonal K from a vector of element stiffnesses, allowing zeros (a failed element) and rejecting negatives. `uniform_chain_condition_analytic(n)` evaluates the eigenvalue ratio. `error_amplification(K, f, df)` solves both systems with `np.linalg.solve` (Lesson 13's tool) and returns the actual relative error and the bound κ‖δ**f**‖/‖**f**‖. `worst_case_loads(K)` returns **u**₁ and **u**_n from the SVD. `mechanism_modes(K)` is `null_space(K, rtol)`.

Running `main()` prints:

```
norms: |K|_1 = 2.800e+08, |K|_inf = 2.800e+08, |K|_F = 5.191e+08, |K|_2 (power) = 2.738e+08
cond_2(K) = 175.087;  numpy.linalg.cond = 175.087;  analytic = 175.087;  16 n^2 / pi^2 = 162.1
tip load 1 kN: tip displacement 0.1429 mm (analytic F L / E A = 0.1429 mm)
1 % random load error: relative displacement error 1.36 % (bound 175.1 %)
1 % worst-case load error: relative displacement error 175.1 % (bound 175.1 %)
```

then the damaged-element table, the mechanism mode and the n-sweep.

### Verify

- **Analytic case:** the 3-element hand-assembled K; the uniform-bar κ equals the eigenvalue formula to 10⁻⁹; the tip deflection equals FL/EA to 10⁻¹²; the strain is uniform; a diagonal matrix's SVD is its sorted diagonal; an orthogonal matrix has 2-norm 1 and κ = 1.
- **Limiting cases:** rank-deficient inputs give exactly two non-zero singular values, κ = ∞ and a completed orthonormal U; a failed element gives one null vector of the predicted shape; an intact bar has an empty null space; scaling K by 10⁻⁹ leaves κ unchanged.
- **Conservation / identities:** A = UΣVᵀ to 10⁻¹², UᵀU = VᵀV = I; AᵀA**v**_i = σ_i²**v**_i; ‖A‖_F² = Σσ²; the induced norms are attained by the predicted vectors; A A⁺ A = A for a singular matrix; the minimum-norm solution is shorter than any other solution.
- **Scaling / bounds:** κ grows with n at order 1.9–2.05; with the damaged-element ratio at order −1 (0.95–1.05); the perturbation bound holds for ten random perturbations, is attained for the worst-case pair, and the reverse pair gives 1/κ; Eckart–Young errors match σ_{k+1} and the tail root-sum-square.

### Compare

`numpy.linalg.svd` calls LAPACK's `gesdd`: Householder bidiagonalisation followed by a divide-and-conquer eigen-solver, about 4mn² + 8n³ flops and backward stable. One-sided Jacobi agrees with it to 10⁻¹⁰ on every test matrix here and is slightly *more* accurate for the smallest singular values of a badly scaled matrix (it never mixes columns of very different size), but it is several times slower at n = 100 and much slower at n = 1000. `numpy.linalg.cond`, `pinv` and `matrix_rank` are all thin wrappers around the SVD, which is why the SVD is the last decomposition in this phase: it is the one the others are judged by. We wrote it to know what σ_n and **v**_n mean for a structure, and to be able to choose the rank threshold ourselves instead of accepting a default.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Map the conditioning of the bar.

1. For n = 2, 4, 8, …, 128 identical elements compute κ₂(K) with `condition_number`, with `numpy.linalg.cond` and with the analytic formula; record all three. Predict the slope on log-log axes before running.
2. For n = 10, damage element 5 to r = 10⁻¹, 10⁻², …, 10⁻⁵ of its stiffness; record κ₂ and the smallest singular value. Predict the exponent of κ against r.
3. For the intact bar and the r = 10⁻³ bar, apply a tip load of 1 kN with (i) ten random 1 % perturbations and (ii) the worst-case pair from `worst_case_loads` scaled to 1 %; record the actual relative displacement errors and the bound.
4. Set r = 0 and record the null-space dimension and mode. Plot κ against n and against r; save as `stiffness_conditioning.png`.

Raw data format: CSV with columns `n, cond_svd, cond_numpy, cond_analytic` and `ratio, cond, sigma_min`.

## 5. Validate and improve

The run gives, for the element-count sweep:

| n | κ₂ from the Jacobi SVD | exact eigenvalue ratio | 16n²/π² |
|---|---|---|---|
| 2 | 6.9 | 6.9 | 6.5 |
| 8 | 113.5 | 113.5 | 103.8 |
| 32 | 1708.7 | 1708.7 | 1660.0 |
| 128 | 26 765.0 | 26 765.0 | 26 560.7 |

The SVD and the exact ratio agree to every printed digit; the asymptotic formula is 6 % low at n = 2 and 0.8 % low at n = 128, as the (2n + 1)² against 4n² correction predicts. For the damaged element (n = 10, element 5): κ = 175, 357, 2364, 2.25 × 10⁴, 2.24 × 10⁵, 2.24 × 10⁶ for r = 1, 10⁻¹, …, 10⁻⁵. From r = 10⁻² on, κ ≈ 22.4/r: the exponent is exactly −1, and the constant is σ₁/σ_n with σ₁ = 2.61 × 10⁸ = 3.73k (barely changed by the damage) and σ_n ≈ 0.166 rk ≈ rk/6 (the soft element in series with the six nodes beyond it, which move almost as a rigid block). At r = 0, σ_n = 0 and `mechanism_modes` returns the single vector (0, 0, 0, 0, 0.408, 0.408, 0.408, 0.408, 0.408, 0.408): nodes 1–4 fixed, nodes 5–10 translating together, which is the mechanism an engineer would sketch.

The amplification table is the point of the lesson:

| bar | load error | actual relative error in **u** | bound κ ‖δ**f**‖/‖**f**‖ |
|---|---|---|---|
| intact | 1 %, random direction | 1.36 % | 175 % |
| intact | 1 %, worst-case pair | 175.1 % | 175.1 % |
| r = 10⁻³ | 1 %, random direction | 1.73 % | 2.25 × 10⁴ % |
| r = 10⁻³ | 1 %, worst-case pair | 2.25 × 10⁴ % | 2.25 × 10⁴ % |

Explain both columns. The bound is attained, to eight digits, by the pair the derivation predicts: that is a *validation of the theorem*, not of the design. The random case sits two orders of magnitude below it for the intact bar and four orders below it for the damaged bar (1.73 % against 2.25 × 10⁴ %), because the worst-case load pattern, **u**₁ = (−0.13, 0.25, −0.34, 0.41, −0.44, 0.43, −0.38, 0.30, −0.19, 0.07), alternates in sign from node to node: it is the highest "mode" of the chain, and a random perturbation has only a 1/√10 projection on it while a physical load (a tip force, a distributed load) has almost none. For a tip load, the displacement is dominated by the *soft* directions, which are also where δ**f** does its damage, so the actual error stays near the 1 % of the input even when κ is 10⁴. Conclusion for the note: the 1 % load uncertainty gives about 1.4 % displacement uncertainty for this bar and this load; κ = 175 is not a concern; κ = 2 × 10⁴ for the damaged bar is a warning about *rounding* (four of sixteen digits) and about any load with a component along the high mode, not a prediction of 10⁴ % errors.

Improvement, feasible now, with a result worth seeing: scale the damaged bar by D = √diag(K) and compute κ(D⁻¹KD⁻¹). Diagonal scaling is the standard first preconditioner, and for a matrix whose rows merely have different *sizes* it removes most of the condition number. Here it does not: for r = 10⁻³ it takes κ from 2.25 × 10⁴ to 2.01 × 10⁴, an 11 % reduction, because the ill-conditioning is a near-mechanism, a physical soft mode of the structure, not an artefact of entry scaling. A preconditioner can undo bad scaling; it cannot undo a nearly free degree of freedom. Report both κ values and say which kind of ill-conditioning this is.

## 6. Defend your work

Submit `main.py`, the test log, `stiffness_conditioning.png` with its CSVs, the amplification table, and a one-page note (template in `projects/_template/report-template.md`) that states κ for the intact and damaged bars, the fitted exponents in n and in r, the actual-versus-bound comparison with the reason for the gap, and the answer to "does the 1 % load uncertainty matter?".

### Adversarial self-checks

- What if the load doubles? Nothing changes: κ is a property of K, the relative error bound is homogeneous in **f**, and the random-case error stays at 1.4 %. What if the *stiffness* doubles everywhere? Also nothing: κ is invariant to scaling K.
- What if the element count goes to one? κ = 1: a single spring cannot amplify anything. What if it goes to 10⁴? κ ≈ 1.6 × 10⁷, and a double-precision solve still has nine good digits; the mesh is not the problem until κ reaches 10¹².
- Which assumption fails first in the real world? Linearity: a bar damaged to 0.1 % stiffness is a crack, and a crack opens and closes with the sign of the load, so K depends on **u** and no single κ describes it.
- What measurement could falsify the model? Static test with strain gauges on each element: the measured strain pattern must be uniform for the intact bar and must concentrate 99.9 % of the elongation in the damaged element for r = 10⁻³; a different pattern means the damage model is wrong, not the arithmetic.
- Does the result survive a change of solver? Jacobi SVD, `numpy.linalg.cond` and the eigenvalue formula agree to 10⁻⁹. Which decision changes if the load uncertainty doubles? None for the intact bar (2.7 % displacement error); for a load with a large high-mode component, the test rig would need to be redesigned to apply smoother loads, since no solver change can help.

## Check yourself

1. State the three norm axioms and give ‖(3, −4)‖₁, ‖·‖₂ and ‖·‖∞.
2. Define the induced matrix norm and prove ‖AB‖ ≤ ‖A‖ ‖B‖ in two lines.
3. Why is ‖A‖₂ equal to the largest singular value, and ‖A‖_F to the root-sum-square of all of them?
4. Write the perturbation bound for K**u** = **f**, name the two vectors that attain it, and say what pairing gives the smallest amplification.
5. Compute κ₂ of [[1, 1], [1, 1.0001]] to two figures and say how many digits a solve can lose.
6. How does the condition number of the uniform bar scale with the element count, and why does a refined mesh make the solve harder?
7. What is the null space of a stiffness matrix with a failed element, physically, and which part of the SVD exposes it?
8. In what sense is A_k = Σ_{i≤k} σ_i**u**_i**v**_iᵀ the best rank-k approximation, and what is its 2-norm error?
