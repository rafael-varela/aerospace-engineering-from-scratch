# P03.08  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: the flop budget)

A finite-element model has n = 100·s unknowns and m load cases.

(a) Give the flop counts of one LU factorisation, one pair of substitutions, one Cholesky factorisation, forming the inverse, and one Thomas solve (if the matrix were tridiagonal).
(b) At 10⁹ flops per second, how long does the LU factorisation take? How many load cases can be solved from the factors in the same time?
(c) A colleague computes K⁻¹ once and multiplies each load vector by it. Compare the cost of that plan with factor-and-substitute for m = 10 and m = 1 000, and name the second reason (not cost) why it is the worse plan.

Answer form: five counts, two times, one count, one sentence.

## 2. (Derivation: pivots, definiteness and the error bound)

(a) Show that if A = LU with L unit lower triangular (no pivoting), the k-th pivot is u_kk = det(A_k)/det(A_{k−1}), where A_k is the leading k×k block of A. Hence show that A = [[0, 1], [1, 0]] has no such factorisation and that a symmetric matrix with all leading minors positive has positive pivots.
(b) Show that if the Cholesky recursion (1) reaches a non-positive quantity under the square root at column j, then A is not positive definite. (Hint: the pivot is the (j, j) entry of the Schur complement, and the Schur complement of an SPD matrix is SPD; or construct an x with xᵀAx ≤ 0.)
(c) Derive the bound ‖x̂ − x‖/‖x‖ ≤ κ(A) ‖r‖/(‖A‖‖x‖) from A(x̂ − x) = −r. Then show that for the Hilbert matrix of order 4 (κ₁ = 2.84 × 10⁴) a residual of 10⁻¹⁶ allows a forward error of up to about 3 × 10⁻¹², and compare with the 1.3 × 10⁻¹³ the code found.

Answer form: three short proofs, two numbers.

## 3. (Unfamiliar scenario: a satellite thermal network)

Four equipment nodes of a small satellite exchange heat by conduction with each other and with a radiator held at T_r = 250 K. Conductances in W/K: G₁₂ = 2 + s/10, G₂₃ = 3, G₃₄ = 1.5, G₁₄ = 1, G₂r = 0.5, G₄r = 4. Dissipations: Q = (20, 10, 0, 5) W at nodes 1–4. Steady state: for each node, Σ_j G_ij (T_i − T_j) = Q_i.

(a) Write the 4×4 conductance matrix G and the right-hand side (the radiator terms move to the right). Show that G is symmetric, that every diagonal entry is at least the sum of the off-diagonal magnitudes in its row, and explain why fixing the radiator temperature makes G positive definite when the unfixed network alone would be singular.
(b) Factor G by Cholesky by hand (four columns) and solve for the four temperatures.
(c) Check: the total heat reaching the radiator equals ΣQ. Which node is hottest, and why is it not the one with the largest dissipation?
(d) Which structural quantity of this lesson does G correspond to, and which solver would you use for a 5 000-node version?

Answer form: a matrix, a vector, a lower-triangular factor, four temperatures with units, two sentences.

## 4. (Numerical, uses your code: the depth of the spar)

Set the truss depth to d = 0.2 + s/100 m (`spar_truss(depth=d)`), E = 70 GPa, A = 2 cm², 4 kN at the lower tip node.

(a) Predict the root upper-chord force F·2/d and the root diagonal force F·√(0.5² + d²)/d before running; then run `solve_truss` and `member_forces` and compare.
(b) Report the tip deflection, the most stressed member and its stress, and κ₂ of the reduced stiffness matrix (`np.linalg.eigvalsh`).
(c) Compare the tip deflection with the d = 0.3 m value (20.32 mm) and with the (0.3/d)² scaling a pure bending model would predict. Explain the difference in one sentence.
(d) Sweep d from 0.1 to 0.6 m in 0.1 m steps and tabulate κ₂. Which way does conditioning go as the truss gets shallower, and what happens to the three zero modes of the unsupported K at d = 0?

Answer form: two predictions and two code values, three numbers, one sentence, one table.

## 5. (Numerical, uses your code: normal equations against QR)

Fit a polynomial of degree p = 4 + ⌊s/8⌋ to 50 points x in [0.2, 1] with y = Σ c_k x^k for random coefficients (`np.random.default_rng(s).standard_normal(p + 1)`), using the Vandermonde matrix V = `np.vander(x, p + 1, increasing=True)`.

(a) Solve the normal equations (VᵀV) c = Vᵀy with `cholesky` and `cholesky_solve`; solve the same problem with `numpy.linalg.qr` and `back_substitution` (`P01.17`); compare both with `numpy.linalg.lstsq` and with the true coefficients.
(b) Report κ₂(V) and κ₂(VᵀV) and verify κ(VᵀV) = κ(V)². Compare the error of each method with κ(VᵀV)·ε and κ(V)·ε.
(c) Increase the degree to 8, 10 and 12 and repeat. At which degree does `cholesky` refuse VᵀV, and what has happened to the matrix?
(d) State the rule this problem illustrates about forming AᵀA.

Answer form: a table of errors by method and degree, two condition numbers, one rule.

## 6. (Judgment: choosing the solver)

For each case say which solver you would use (dense LU, Cholesky, banded or Thomas, a sparse-direct solver, or an iterative method of `P03.09`) and why, in one sentence each: (a) the 16×16 truss of this lesson, 200 load cases; (b) a 1-D heat-conduction problem with 10⁶ cells; (c) a 2-D potential-flow Laplace problem on a 1 000 × 1 000 grid; (d) the 6×6 Jacobian systems inside a Newton iteration for aircraft trim, 20 iterations; (e) a 10⁵-degree-of-freedom wing finite-element model with one load case; (f) a least-squares fit with 10⁴ observations and 8 parameters.

Answer form: six labelled lines.
