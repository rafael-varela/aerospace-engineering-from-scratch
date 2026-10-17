# P03.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. The flop budget ($s = 15$: n = 1 500)

(a) LU: 2n³/3 = **2.25 × 10⁹**; substitutions: 2n² = **4.5 × 10⁶**; Cholesky: n³/3 = **1.13 × 10⁹**; inverse: 2n³ = **6.75 × 10⁹**; Thomas: 8n = **1.2 × 10⁴** flops.

(b) 2.25 × 10⁹/10⁹ = **2.25 s** for the factorisation; in the same time **500** load cases (= n/3) can be solved from the factors.

(c) Inverse plan: 6.75 × 10⁹ + m × 4.5 × 10⁶; factor-and-substitute: 2.25 × 10⁹ + m × 4.5 × 10⁶. At m = 10: 6.80 × 10⁹ against 2.30 × 10⁹ (3.0×); at m = 1 000: 1.13 × 10¹⁰ against 6.75 × 10⁹ (1.7×). The per-case cost is the same; only the set-up differs, and the inverse's set-up is three times larger. The second reason: multiplying by a computed A⁻¹ is not backward stable (the residual of A⁻¹b can be κ times larger than that of two substitutions), so the inverse plan is both slower and less accurate.

**Check.** The crossover where the inverse plan ever wins does not exist: its per-case cost is not lower.

**Common mistakes.** Counting multiply and add as one flop (the counts above use two); believing that A⁻¹b is "one matrix–vector product and therefore cheap" without counting what produced A⁻¹.

## 2. Pivots, definiteness and the error bound

(a) If A = LU with L unit lower triangular, the leading k×k block of A is the product of the leading blocks, A_k = L_k U_k, because the entries of A_k involve only the first k rows of L and the first k columns of U, and L_k, U_k are those blocks. Then det A_k = det L_k det U_k = Π_{i≤k} u_ii, and dividing consecutive determinants, **u_kk = det A_k/det A_{k−1}** (det A₀ = 1). For [[0, 1], [1, 0]]: det A₁ = 0, so u₁₁ = 0, and u₂₂ = det A₂/det A₁ = −1/0 is undefined: no factorisation. For a symmetric matrix with all leading minors positive, every ratio is positive, so every pivot is positive, and A = LDLᵀ with D > 0: that is Sylvester's criterion meeting elimination.

(b) After j − 1 steps of symmetric elimination, A = [[A₁₁, A₁₂], [A₁₂ᵀ, A₂₂]] with A₁₁ the (j − 1)×(j − 1) block already factored (A₁₁ = L₁₁L₁₁ᵀ, so it is SPD), and the quantity under the square root in (1) is the (1, 1) entry of the Schur complement S = A₂₂ − A₁₂ᵀA₁₁⁻¹A₁₂ (the sum Σ_{k<j} L_jk² is exactly A₁₂ᵀA₁₁⁻¹A₁₂ at that entry). Take x = (−A₁₁⁻¹A₁₂e₁, e₁): xᵀAx = e₁ᵀSe₁ = s₁₁. If s₁₁ ≤ 0 then xᵀAx ≤ 0 for a nonzero x, so **A is not positive definite**. Conversely, if A is SPD then S is SPD and every pivot is positive, so the recursion never fails: the Cholesky pivot test is exact.

(c) From Ax̂ = b − r and Ax = b, A(x̂ − x) = −r, so ‖x̂ − x‖ ≤ ‖A⁻¹‖‖r‖ = κ(A)‖r‖/‖A‖, and dividing by ‖x‖,

**‖x̂ − x‖/‖x‖ ≤ κ(A) · ‖r‖/(‖A‖‖x‖).**

For H₄ with x = (1, 1, 1, 1): ‖H₄‖₁ = 1 + 1/2 + 1/3 + 1/4 = 2.083, ‖x‖₁ = 4, so a residual ‖r‖₁ = 10⁻¹⁶ allows a forward error up to 2.84 × 10⁴ × 10⁻¹⁶/(2.083 × 4) = **3.4 × 10⁻¹³** relative; the code found 1.3 × 10⁻¹³, within the bound (the bound in the lesson, κε = 6.3 × 10⁻¹², uses the normalised backward error and is looser).

**Common mistakes.** Applying (a) with pivoting (then the minors are those of PA); reading (b) as "a positive pivot proves positive definiteness" (all pivots must be positive); forgetting ‖A‖ in the bound.

## 3. Satellite thermal network ($s = 15$: G₁₂ = 3.5 W/K)

**Assumptions.** Steady state, conduction only, lumped nodes, constant conductances.

**Derivation.** (a) Node balances Σ_j G_ij(T_i − T_j) = Q_i with T_r = 250 K known:

G = [[4.5, −3.5, 0, −1.0], [−3.5, 7.0, −3.0, 0], [0, −3.0, 4.5, −1.5], [−1.0, 0, −1.5, 6.5]] W/K,  rhs = (20, 10 + 0.5 × 250, 0, 5 + 4 × 250) = **(20, 135, 0, 1 005) W**.

G is symmetric (each conductance appears in both its rows). Row sums of off-diagonal magnitudes: 4.5, 6.5, 4.5, 2.5, against diagonals 4.5, 7.0, 4.5, 6.5: every row is diagonally dominant and rows 2 and 4 (the ones touching the radiator) strictly so. Without the radiator, every row would sum to zero, G·(1, 1, 1, 1)ᵀ = 0, and a uniform temperature shift would be a zero mode (the thermal analogue of a rigid-body mode); fixing one temperature removes that mode and leaves G positive definite (an irreducibly diagonally dominant symmetric matrix with positive diagonal and at least one strictly dominant row is SPD).

(b) Cholesky by hand (W/K)^{1/2}: L₁₁ = √4.5 = 2.1213; L₂₁ = −3.5/2.1213 = −1.6499; L₃₁ = 0; L₄₁ = −1/2.1213 = −0.4714. L₂₂ = √(7 − 1.6499²) = √4.2778 = 2.0683; L₃₂ = −3/2.0683 = −1.4505; L₄₂ = (0 − (−0.4714)(−1.6499))/2.0683 = −0.3761. L₃₃ = √(4.5 − 1.4505²) = √2.3961 = 1.5479; L₄₃ = (−1.5 − (−0.3761)(−1.4505))/1.5479 = −1.3214. L₄₄ = √(6.5 − 0.4714² − 0.3761² − 1.3214²) = √4.3903 = 2.0953. Forward substitution Ly = rhs: y = (9.428, 72.79, 68.21, 537.8); back substitution Lᵀx = y: **T = (268.72, 266.44, 263.19, 256.69) K**.

(c) Heat to the radiator: 0.5 × (266.44 − 250) + 4 × (256.69 − 250) = 8.22 + 26.78 = **35.0 W** = ΣQ ✓. Node 1 is hottest although node 4 also dissipates and node 2 dissipates 10 W: node 1 has no path to the radiator except through nodes 2 and 4, and its 20 W must cross 3.5 W/K and 1 W/K of conductance; node 4 sits on the 4 W/K link to the radiator and stays within 7 K of it.

(d) G is a stiffness matrix: symmetric, assembled from element contributions g[[1, −1], [−1, 1]], singular until a "support" (the fixed temperature) is applied, then SPD; the energy xᵀGx/2 is the dissipation. For 5 000 nodes, sparse Cholesky (the matrix has a few nonzeros per row) or conjugate gradients (`P03.09`); dense Cholesky at n³/3 = 4 × 10¹⁰ flops would still run in seconds but wastes 200 MB on zeros.

**Check.** `cholesky_solve(cholesky(G), rhs)` returns the same four temperatures to 10⁻¹²; the eigenvalues of G are 0.64, 4.40, 6.66, 10.80 W/K, all positive.

**Common mistakes.** Putting the radiator conductances only on the diagonal and forgetting their G_ir T_r term on the right; writing G_ij = +G_ij.

## 4. The depth of the spar ($s = 15$: d = 0.35 m)

(a) Predictions: upper chord F·2/d = 4 × 2/0.35 = **22.86 kN** (tension); diagonal −F√(0.25 + 0.1225)/0.35 = −4 × 0.6103/0.35 = **−6.98 kN**. Code: N₁ = +22.857 kN, N₂ = −6.975 kN: agreement to 10⁻⁹ relative.

(b) Tip deflection **−15.35 mm**; most stressed member: the root upper chord at 22.86 kN, **114.3 MPa**; κ₂ = **1 643**.

(c) Ratio to the d = 0.3 m case: 15.35/20.32 = **0.755**; pure bending predicts (0.3/0.35)² = **0.735**. The truss deflection has a shear part, carried by the diagonals and verticals, whose stiffness does not grow with d² (the diagonals get *longer* and more steeply inclined as d grows), so the deflection falls a little more slowly than 1/d².

(d) κ₂ for d = 0.1, 0.2, 0.3, 0.4, 0.5, 0.6 m: 3.8 × 10⁴, 5.2 × 10³, 2.2 × 10³, 1.3 × 10³, 8.5 × 10², 6.4 × 10² (from `np.linalg.eigvalsh` of the reduced K; your values should agree to the figures shown). Conditioning worsens as the truss gets shallower: the softest mode (bending of a shallow truss) gets softer while the stiffest (axial stretching of a chord) does not change. At d = 0 all nodes are collinear, every member is horizontal, the structure cannot resist any vertical load, and the three rigid-body modes become a continuum of mechanisms: the reduced K is singular.

**Check.** Equilibrium holds at every d (reactions sum to 4 kN); the root chord force times d equals 8 kN·m at every d.

**Common mistakes.** Reading the most stressed member off the deflection plot instead of the force vector; comparing κ₂ between runs with different units of K.

## 5. Normal equations against QR ($s = 15$: p = 5)

(a) Max coefficient error: **Cholesky on VᵀV 1.3 × 10⁻⁸**; **QR 4.3 × 10⁻¹³**; **lstsq 1.1 × 10⁻¹³**. The three solve the same least-squares problem; the first loses five more digits.

(b) κ₂(V) = **1.62 × 10⁴**, κ₂(VᵀV) = **2.63 × 10⁸** = κ(V)² ✓. The bounds: κ(VᵀV)ε = 5.8 × 10⁻⁸ (Cholesky's error 1.3 × 10⁻⁸ sits just below it); κ(V)ε = 3.6 × 10⁻¹² (QR's 4.3 × 10⁻¹³ sits below that). Each method is backward stable *for the problem it solves*; the normal-equations route solves a problem whose conditioning is the square.

(c) Degree 8: κ(V) = 7.1 × 10⁶, Cholesky error 1.0 × 10⁻³, QR 8.6 × 10⁻¹¹. Degree 10: κ(V) = 4.2 × 10⁸, Cholesky error 1.0 (no digits), QR 2.8 × 10⁻⁹. Degree 12: κ(V) = 2.6 × 10¹⁰, κ(VᵀV) = 6.6 × 10²⁰ > 1/ε: **`cholesky` refuses VᵀV** with a negative pivot (−5.8 × 10⁻¹³ at column 12), because the computed VᵀV is no longer positive definite in floating point: its rounding errors exceed its smallest eigenvalue. QR still returns seven correct digits.

(d) **Never form AᵀA to solve a least-squares problem unless κ(A) is small**: it squares the condition number, and with κ(A) > 10⁸ the normal equations are numerically singular in double precision. Use QR (`P01.17`) or the SVD.

**Common mistakes.** Blaming Cholesky for the loss (it is backward stable; the matrix handed to it is the problem); using `np.vander` with the default decreasing powers and mismatching the coefficients.

## 6. Choosing the solver

(a) **Cholesky**, once (16³/3 ≈ 1 400 flops), then 200 × 2 × 16² ≈ 10⁵ flops of substitutions: the factorisation is 1 % of the work.
(b) **Thomas**: 8 × 10⁶ flops and three arrays of 10⁶; a dense matrix would be 8 TB.
(c) **Sparse-direct with nested dissection, or preconditioned CG** (`P03.09`): 10⁶ unknowns with five nonzeros per row; dense LU would be 7 × 10¹⁷ flops.
(d) **Dense LU** (`np.linalg.solve`): 6×6 is 144 flops; any cleverness costs more than it saves, and the Jacobian is not symmetric.
(e) **Sparse-direct Cholesky** (SPD, one load case, fill-in controlled by reordering), or PCG if memory for the factor is short; the model is far beyond dense storage (10¹⁰ entries).
(f) **QR** of the 10⁴ × 8 matrix (Householder, 2 × 10⁴ × 64 flops): the normal equations would square a condition number that may already be large for correlated regressors.

## Check-yourself answers

1. LU 2n³/3; one solve 2n²; inverse 2n³; Cholesky n³/3.
2. The multiplier a_ik/a_kk needs a_kk ≠ 0 and small |a_kk| amplifies rounding. [[0, 1], [1, 0]] has no LU; [[10⁻¹⁷, 1], [1, 1]] with b = (1, 2) gives (0, 1) instead of (1, 1) without pivoting.
3. ρ = max|U|/max|A|; bounded by 2ⁿ⁻¹ under partial pivoting; Wilkinson's matrix (1 diagonal, −1 below, 1 last column) attains it.
4. A non-positive quantity under the square root at column j proves A is not positive definite (an explicit x with xᵀAx ≤ 0 exists); for SPD A the recursion never fails.
5. η = ‖b − Ax̂‖/(‖A‖‖x̂‖ + ‖b‖); ‖x̂ − x‖/‖x‖ ≲ κ(A) η.
6. About 6: the forward error can be up to κη ≈ 10¹⁰ × 10⁻¹⁶ = 10⁻⁶.
7. Each elimination step touches one sub-diagonal entry and each back-substitution one super-diagonal entry: 8n flops; it does not pivot, so it is unsafe when the matrix is not diagonally dominant (or SPD), e.g. a convection-dominated stencil.
8. Moment about the root, 4 kN × 2 m = 8 kN·m, divided by the depth 0.3 m: 26.67 kN tension, by statics alone (the root bay is determinate).
