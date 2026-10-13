# P01.13  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are floating-point, with the residual where you have one.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: elimination by hand)

Solve by hand, showing the augmented array after each column is cleared:

2x + y − z = s − 7,  −3x − y + 2z = −11,  −2x + y + 2z = −3.

(a) First without row swaps. List the three pivots and their product.
(b) Then with partial pivoting. List the pivots, the row order, and their product with its sign.
(c) Write the L and U of (a) and multiply them back to A. Check your x with `gaussian_elimination` and your L, U with `lu_decompose` (whose factors will differ from yours by the pivoting; say how).

Answer form: three numbers; two pivot lists; two 3 × 3 matrices and a residual.

## 2. (Derivation: the cost and the factorisation)

(a) Derive the operation count of forward elimination, Σ_k 2(n − k)(n − k + 1), and show it is 2n³/3 + O(n²). Add back substitution. Then derive the cost of Gauss–Jordan inversion and show it is about three solves.
(b) Show that the elimination step "row i ← row i − ℓ row k" is left-multiplication by E = I − ℓ **e**_i **e**_kᵀ, that E⁻¹ = I + ℓ **e**_i **e**_kᵀ, and that the product of all the E⁻¹ in reverse order is unit lower triangular with the multipliers in place. Conclude A = LU.
(c) Prove that det A equals the product of the pivots times (−1)^{number of row swaps}. Which properties of the determinant do you need?

Answer form: three derivations with each step named.

## 3. (Unfamiliar scenario: a spacecraft thermal network)

An electronics box (node 1) dissipates 40 W. It conducts to a mounting panel (node 2) through G₁₂ = 2 W/K; the panel conducts to a radiator (node 3) through G₂₃ = 1.5 W/K; the radiator sheds heat to deep space, linearised as G₃ₛ = 0.5 W/K to a reference temperature of 0 (all temperatures below are rises above that reference). Heat balance at each node: Σ G (T_other − T_node) + Q_node = 0.

(a) Write the three balances as G **T** = **Q** and show G is symmetric. What is the sign pattern of its rows and why?
(b) Solve by hand (elimination) and with `gaussian_elimination`. Verify that all 40 W leave through the radiator.
(c) Double G₂₃. Which temperatures change and which do not? Explain physically.
(d) Compute G⁻¹ with `inverse_gauss_jordan` and interpret its columns: what does G⁻¹[i, j] mean in watts and kelvin? Why are all its entries positive?

Answer form: a matrix; three temperatures in K; two sentences; a 3 × 3 matrix with an interpretation.

## 4. (Numerical, uses your code: cost, pivoting, conditioning)

(a) Time `gaussian_elimination` and `np.linalg.solve` for n = 50, 100, 200, 400, 800 (best of three). Fit the exponent of each. Compute your solver's speed in flop/s from `elimination_flops`.
(b) For ε = 10⁻⁸, 10⁻¹⁰, 10⁻¹², 10⁻¹⁴ solve [[ε, 1], [1, 1]] x = (1, 2) with `pivot=True` and `pivot=False`, and tabulate the error in x₁ against the prediction 10⁻¹⁶/ε. Add one, two and three steps of `iterative_refinement` with `pivot=False`.
(c) For the Hilbert matrices n = 4, 6, 8, 10, tabulate κ(H), the residual of `gaussian_elimination(H, H @ ones)` and its maximum error, before and after three refinement steps. State in one sentence what refinement can and cannot do.

Answer form: three tables; two fitted exponents; one sentence.

## 5. (The truss with your load)

Add to the lesson's two-bay truss a second load of s kN downward at node 4 (the top middle joint), keeping the 15 kN at node 1.

(a) Compute the reactions by global equilibrium by hand before running anything.
(b) Run `solve_truss` and list the member forces. Which member carries the largest compression? Which members are zero-force, and why (argue from the joints, not the numbers)?
(c) Now remove the s kN load and instead tilt the 15 kN load at node 1 by adding a horizontal component of 3 kN in +x. Which forces change and which do not? Which reaction takes the horizontal load?
(d) Replace the diagonal 1-5 by a member 3-5. The count m + r = 2j is unchanged. Run `solve_truss` and explain, in terms of the right bay's geometry, what the error message means.

Answer form: three numbers; a table of nine forces with two identified members; two sentences; one paragraph.

## 6. (Judgment: when to distrust a solve)

For each case, say whether the solver's answer is trustworthy, give the measurable sign you would look at, and name the cause:

(a) A 12 × 12 truss solve returns a residual of 10⁻¹⁶ and member forces of ±10⁵ kN for a 15 kN load.
(b) A 1000 × 1000 stiffness matrix solve returns a residual of 10⁻³ relative.
(c) `np.linalg.cond(A)` is 10¹⁴ and the forces look reasonable.
(d) A thermal network solve returns one negative temperature rise with all heat sources positive.

Answer form: four short paragraphs, each with a number or a measurable quantity.

## 7. (Hard, uses your code: a tridiagonal system in O(n))

A 1 km space tether is modelled as n = 1000 segments of 1 m, each a spring of stiffness k = 10⁵ N/m (EA/L for a 1 m segment of Kevlar line), with a mass m = 0.05 kg lumped at each node and gravity g = 9.80665 m/s² pulling along the tether. The top node is fixed; the displacement u_i of node i satisfies k(2u_i − u_{i−1} − u_{i+1}) = m g for interior nodes and k(u_n − u_{n−1}) = m g at the free end.

(a) Write the stiffness matrix K: it is tridiagonal with 2k on the diagonal (k in the last row) and −k off it. Derive by hand the exact tip displacement and the tension in the top segment.
(b) Implement the Thomas algorithm (elimination that only touches the three diagonals: forward sweep, back sweep) and solve the system. Count its operations and compare with `elimination_flops(1000)`.
(c) Time the Thomas algorithm, `gaussian_elimination` and `np.linalg.solve` on this system and compare all three solutions with (a).
(d) Why does partial pivoting never trigger a swap for this matrix? (Hint: compare the diagonal with the sum of the off-diagonal magnitudes in each row.)

Answer form: a derivation and two numbers; code and an operation count; three times and three errors; one argument.
