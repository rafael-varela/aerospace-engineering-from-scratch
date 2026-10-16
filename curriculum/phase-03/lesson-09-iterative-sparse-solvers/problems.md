# P03.09  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: storage and the flop budget)

A plate problem is discretised on an n × n interior grid with n = 100 + 10·s, N = n² unknowns, five-point stencil.

(a) Give the number of nonzeros of A, the dense storage in bytes (8 per entry) and the CSR storage (8 bytes per value, 4 per column index, 4 per row pointer).
(b) Give the cost of one matrix–vector product and of one CG iteration (matvec plus about ten vector operations of length N).
(c) Using κ ≈ (2(n + 1)/π)², predict the number of CG iterations for a 10⁻⁶ reduction from the bound (5) and the number of Jacobi sweeps from (3). Give the total flops of each.
(d) Compare with dense Cholesky (N³/3) and banded Cholesky (about N n²) at 10⁹ flops per second.

Answer form: four storage numbers, two per-iteration costs, two iteration counts, four times.

## 2. (Derivation: eigenvalues, Jacobi, and the CG recurrence)

(a) Show that v_p = (sin(pπ i/(n + 1)))_{i=1..n} is an eigenvector of T = tridiag(−1, 2, −1) with eigenvalue 2 − 2cos(pπ/(n + 1)), and hence derive (2) for the 2-D matrix (k/h²)(I ⊗ T + T ⊗ I). Show that κ = cot²(π/(2(n + 1))) and expand for large n.
(b) Derive the Jacobi iteration matrix G = I − D⁻¹A for the 2-D Laplacian and show that its eigenvalues are (cos(pπ/(n + 1)) + cos(qπ/(n + 1)))/2, so that ρ_J = cos(π/(n + 1)). Expand ln ρ_J to show that the number of sweeps for a reduction ε is about 2(n + 1)² ln(1/ε)/π².
(c) For SPD A and φ(x) = ½xᵀAx − bᵀx, show that the minimiser of φ along x + αp is α = pᵀr/pᵀAp with r = b − Ax, and that the new residual satisfies pᵀr_new = 0. Then show that if p_m = r_m + βp_{m−1} with β = r_mᵀr_m/r_{m−1}ᵀr_{m−1}, and r_mᵀr_{m−1} = 0, then p_mᵀAp_{m−1} = 0. (Hint: Ap_{m−1} = (r_{m−1} − r_m)/α_{m−1}.)

Answer form: three derivations; the final formulas boxed.

## 3. (Unfamiliar scenario: a heated fin with an insulated tip)

A straight aluminium fin (k = 200 W/(m·K)) of length 0.1 m is heated uniformly by eddy currents at q = (10 + s) × 10⁵ W/m³; its root is at 300 K and its tip is insulated. Discretise with four unknowns at x = h, 2h, 3h, 4h (h = 0.025 m, node 4 at the tip), using a ghost node T₅ = T₃ for the insulated end and halving the tip equation so that the matrix stays symmetric.

(a) Write the 4 × 4 matrix A and the right-hand side b (the root temperature enters the first equation). Show that A is SPD without computing eigenvalues (irreducible diagonal dominance, or a Cholesky by hand).
(b) Perform two Jacobi sweeps from x = 0 and one steepest-descent step from x = 0; give the residual norm after each.
(c) Compute ρ(G) for Jacobi numerically (or from the eigenvalues of A, which are (k/h²)(2 − 2cos((2p − 1)π/9)), p = 1, …, 4; note that D is not a multiple of I here, so ρ_J is *not* simply 1 − λ_min h²/(2k)). Predict the Jacobi and Gauss–Seidel sweeps for a 10⁻⁶ reduction and verify by running your code.
(d) Solve exactly and compare with the continuous solution T(x) = 300 + (q/k)(Lx − x²/2). Why is the discrete solution exact here? In how many iterations must CG finish, and does it?

Answer form: a matrix, a vector, two residual norms after Jacobi and one after steepest descent, ρ_J, two sweep counts (predicted and measured), four temperatures with units, two sentences.

## 4. (Numerical, uses your code: the insert and the tolerance)

Take the plate of the lesson (n = 32, L = 0.2 m, kt = 0.34 W/K, q = 2 kW/m², edges at 300 K) with a 0.08 m square insert centred on the plate whose conductivity is kt/s (`plate_conductivity(32, 0.2, kt, kt/s, (0.06, 0.14, 0.06, 0.14))`).

(a) Compute κ₂ of A with `scipy.sparse.linalg.eigsh` (`which="SA"` and `"LA"`) and compare with the uniform plate's 440.7. Report the maximum temperature rise and compare with the uniform plate's 17.29 K.
(b) Solve to 10⁻⁸ with CG, Jacobi-PCG and IC(0)-PCG; report iteration counts and the maximum error against `spsolve`.
(c) Solve with plain CG at tolerances 10⁻⁴, 10⁻⁶, 10⁻⁸, 10⁻¹⁰; tabulate iterations, ‖r‖/‖b‖ reached, and the maximum and 2-norm error against `spsolve`. Compare the 2-norm error with the bound ‖r‖/λ_min.
(d) The discretisation error of the uniform plate at n = 32 is about 0.03 K. Which tolerance would you use for the temperature, and which for the heat flux into the frame? Justify in two sentences.

Answer form: κ and a rise, three counts with errors, a four-row table, two sentences.

## 5. (Numerical, extends your code: successive over-relaxation)

Write `sor(A, b, omega, tol, maxiter)`: the Gauss–Seidel sweep with x_i ← (1 − ω)x_i + ω·(Gauss–Seidel value). For the manufactured problem at n = 8 + s:

(a) Compute ρ_J, ω_opt = 2/(1 + sin(π/(n + 1))) and the predicted ρ_SOR = ω_opt − 1.
(b) Run Jacobi, Gauss–Seidel, SOR(ω_opt) and CG to 10⁻⁶; report the sweep counts and compare with the predictions ln(10⁻⁶)/ln ρ for each.
(c) Run SOR for ω = 1.0, 1.5, 1.7, ω_opt, 1.9, 1.95 and plot or tabulate the counts. Is the optimum sharp on both sides?
(d) The measured SOR contraction over the last sweeps is a little above ω_opt − 1 and the count exceeds the prediction by 30 %. Explain (hint: the SOR iteration matrix at ω_opt has a repeated eigenvalue with a single eigenvector).

Answer form: three numbers, four counts against four predictions, a six-row table, one paragraph.

## 6. (Judgment: choosing the solver)

For each case say which solver you would use (dense LU or Cholesky, banded, sparse-direct, Jacobi/Gauss–Seidel, CG, preconditioned CG, GMRES) and name the preconditioner if any, in one sentence each: (a) the 64 × 64 plate of this lesson, one right-hand side; (b) the same plate, 500 flux distributions; (c) a 3-D thermal model of a satellite structure with 2 × 10⁶ tetrahedra, SPD; (d) the convection–diffusion matrix of a boundary-layer solver, nonsymmetric, 10⁵ unknowns; (e) a 50 × 50 SPD matrix inside a Newton loop; (f) a thermal model with a 10⁴:1 contrast in conductivity between a heat pipe and foam insulation; (g) a transient simulation that solves the same SPD system at every time step with a slowly changing right-hand side.

Answer form: seven labelled lines.
