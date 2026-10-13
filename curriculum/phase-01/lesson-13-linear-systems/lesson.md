# P01.13  Systems of Linear Equations

> Elimination is the one algorithm every structural, thermal and flow solver runs at its heart. Write it once by hand so that you know what the library is doing, and why it sometimes cannot.

| | |
|---|---|
| **Learning outcomes** | 1. Set up the equilibrium of a planar truss as one linear system A x = b and solve a 3 × 3 system by hand with elimination and back substitution. 2. Implement Gaussian elimination with partial pivoting, LU factorisation, a determinant and a Gauss–Jordan inverse from scratch in NumPy, and test each against `numpy.linalg`. 3. Derive the 2n³/3 operation count and measure it. 4. Explain, with a numerical example, why pivoting is necessary, why a tiny residual does not guarantee a small error, and why iterative refinement fixes one of those problems but not the other. 5. Decide from the equilibrium matrix alone whether a truss is a mechanism. |
| **Prerequisites** | `P01.07` (vectors, unit vectors, dot product), `P00.03` (NumPy slicing and fancy indexing), `P00.07` (pytest). Bridge B1: you can solve two equations in two unknowns by substitution. Skill check: you can resolve a force along a unit vector and write Σ F = 0 at a pin joint. |
| **Engineering question** | A 1.5 t engine hangs from the middle joint of a two-bay planar truss on a test stand; what force does each of the nine members carry, which are in compression, and how do we know the numbers are right? |
| **Assumptions and limits** | Pin-jointed truss: members carry axial force only, loads act at joints, deformations are small enough that equilibrium is written on the undeformed geometry, and the structure is statically determinate (as many unknowns as equations and a non-singular equilibrium matrix). Elimination in floating point: the answer is exact for the problem *as stored*; a relative error of order (condition number) × 10⁻¹⁶ is the best any solver can do, and a singular or nearly singular matrix means the *question* is ill-posed, not the solver. |
| **Analytical work** | Method of joints written as a matrix; hand elimination of a 3 × 3 system with pivoting; the operation count; elimination as a product of elementary matrices, hence A = LU; the determinant as the product of pivots; global equilibrium as a check on the reactions. |
| **Computational deliverable** | `starter/main.py`: `back_substitution`, `forward_substitution`, `forward_elimination`, `gaussian_elimination`, `lu_decompose`, `lu_solve`, `determinant`, `inverse_gauss_jordan`, `residual_norm`, `iterative_refinement`, `elimination_flops`, `truss_equilibrium_matrix`, `truss_load_vector`, `solve_truss`. Tests prove: agreement with `numpy.linalg.solve`, `det` and `inv` to 10⁻⁹ on random systems; P A = L U with unit lower L and bounded multipliers; residuals below 10⁻¹²; correct handling of zero and tiny pivots; singular and malformed input raise; the truss satisfies joint and global equilibrium, is linear in the load, detects zero-force members and raises for a mechanism. |
| **Practical activity** | **computational experiment**: timing and error of the hand-written solver against LAPACK as n grows, the effect of pivoting on a tiny-pivot system, and the Hilbert matrix as a case where the residual is tiny and the answer wrong. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `truss_forces.png`, the hand solution of the 3 × 3 system with its pivot sequence, the timing table with the fitted exponent, the Hilbert table, and a short note on when to trust a solve. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: the truss system is set up with the sign convention stated and the reactions checked by global equilibrium; all public tests pass; the operation count is derived and the measured exponent reported; the judgment problem names a measurable sign of an untrustworthy solve. |
| **References** | `@mit-ocw-18-06-linear-algebra` Lecture 1: The geometry of linear equations (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 2: Elimination with matrices (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 3: Multiplication and inverse matrices (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 4: Factorization into A = LU (required-free); `@openstax-precalculus-2e` 9.5 Matrices and Matrix Operations, 9.6 Solving Systems with Gaussian Elimination, 9.7 Solving Systems with Inverses (supplementary-free, a gentler text backup). |
| **Extensions** | ★ Implement the Thomas algorithm for tridiagonal systems (problem 7) and show O(n) cost; ★ store L and U in one array in place, as LAPACK does; ★ add a 3D truss (three equations per joint) and solve a tetrahedral space-frame. |
| **Time** | unmeasured |

## 1. The engineering problem

A test stand holds a 1.5 t engine from the middle joint of the bottom chord of a two-bay truss, 4 m long and 1.5 m tall, pinned at one end and on a roller at the other. Before the stand is used, the engineer must know the force in every member to check it against the allowable stress, and must know which members are in compression, because a slender bar in compression can buckle long before it yields. With nine members and three support reactions there are twelve unknowns. At each of the six pin joints two equations say that the forces balance: twelve equations. The problem is a system of twelve linear equations in twelve unknowns, and the method of joints that Phase 06 teaches by hand, joint by joint, is nothing but a way of solving that system with pencil.

This lesson writes the system as a matrix equation **A x = b** and solves it by **Gaussian elimination**, the algorithm that every finite-element, finite-volume and circuit solver runs at its core. We write it ourselves, in NumPy, not because `numpy.linalg.solve` is inadequate (it calls LAPACK, which is the right choice in practice, and the lesson ends by measuring how much faster it is) but because an engineer who has never watched a pivot go to zero does not know what a singular stiffness matrix is telling them. The lesson's second question is therefore: when the solver returns a number, what has it actually guaranteed, and what has it not?

Success criteria, written before we start: the hand-written solver agrees with `numpy.linalg.solve` to 10⁻⁹ relative on random systems up to n = 60 and leaves a residual ‖b − Ax‖/‖b‖ below 10⁻¹²; the LU factors reconstruct P A to 10⁻¹²; the truss solution satisfies every joint equation and the three global equilibrium equations, and the reactions come out as 7.5 kN each by symmetry; removing one diagonal makes the matrix singular and the code says so rather than returning garbage; and the operation count scales as n³, measured.

## 2. Understand and derive

### 2.1 From joints to a matrix

```
   3 ──────── 4 ──────── 5        y = 1.5 m        members: 0-1, 1-2 (bottom chord)
   │  ╲       │  ╲       │                                  3-4, 4-5 (top chord)
   │    ╲     │    ╲     │                                  0-3, 1-4, 2-5 (verticals)
   │      ╲   │      ╲   │                                  0-4, 1-5 (diagonals)
   0 ──────── 1 ──────── 2        y = 0           supports: pin at 0 (R0x, R0y), roller at 2 (R2y)
   ▲△         │          ○        x = 0, 2, 4 m   load:     15 kN down at node 1
            15 kN
```

**Assumption 1 (pin-jointed truss):** every member is a two-force member, so the force it exerts on a joint is along the member. Let F_m be the force in member m, **tension positive**, and **e** the unit vector from node i to node k along the member. A member in tension pulls node i toward k: force +F_m **e** on node i and −F_m **e** on node k (Newton's third law). **Assumption 2 (loads at joints, small deflections):** external loads **P**_i act at joints and the geometry does not change under load.

Equilibrium at node i, in each of x and y:

Σ_{members at i} F_m (±**e**_m) + **R**_i + **P**_i = **0**.   (1)

Collect the m member forces and r reactions in one unknown vector **x** = (F₁, …, F_m, R₁, …, R_r). Then (1) at all joints is

**A x = −P**,   (2)

where **A** has 2j rows (two per joint) and m + r columns. The column of member (i, k) holds +**e** in the two rows of node i and −**e** in the two rows of node k; the column of a reaction holds a single 1 in the row of its node and direction. The right-hand side is minus the load vector. For the two-bay truss, 2j = 12 and m + r = 9 + 3 = 12: **A** is square. That count, m + r = 2j, is necessary for a determinate truss but not sufficient; §2.6 and `P01.14` return to it.

This is the **column picture** of `@mit-ocw-18-06-linear-algebra` Lecture 1: the load vector −**P** must be written as a combination of the twelve columns of **A**, and the coefficients are the forces. The **row picture** is twelve planes in twelve dimensions meeting at a point. Lecture 1 is worth watching for the 3 × 3 case before you believe the 12 × 12 one.

### 2.2 Elimination by hand, with a pivot

Take the 3 × 3 system (problem 1 uses it with your seed):

2x + y − z = 8,  −3x − y + 2z = −11,  −2x + y + 2z = −3.

Write it as an augmented array [A | b] and eliminate column by column (`@mit-ocw-18-06-linear-algebra` Lecture 2). The entry used to clear the column below it is the **pivot**. Without swapping rows the first pivot is 2: subtract (−3/2) × row 1 from row 2 and (−2/2) × row 1 from row 3:

```
[ 2   1  -1 |  8 ]        [ 2    1    -1  |  8 ]        [ 2   1   -1  |  8 ]
[-3  -1   2 | -11]  --->  [ 0   0.5   0.5 |  1 ]  --->  [ 0  0.5  0.5 |  1 ]
[-2   1   2 | -3 ]        [ 0    2     1  |  5 ]        [ 0   0   -1  |  1 ]
```

The pivots are 2, 0.5, −1. **Back substitution** from the bottom: z = −1, then 0.5y + 0.5(−1) = 1 gives y = 3, then 2x + 3 + 1 = 8 gives x = 2. Check in the original third equation: −4 + 3 − 2 = −3. ✓

With **partial pivoting**, at each column we first swap into the pivot row whichever remaining row has the largest |entry| in that column. Column 1: |−3| is largest, so rows 1 and 2 swap and the pivots become −3, 5/3, 1/5; the answer is the same (2, 3, −1), the arithmetic is different. The product of the pivots is 2 × 0.5 × (−1) = −1 without pivoting and −3 × 5/3 × 1/5 = −1 with pivoting; the row order after pivoting is (2, 3, 1), a cyclic permutation of even sign, so both give det A = −1 (§2.4).

**Why pivot at all?** Consider

10⁻¹³ x + y = 1,  x + y = 2,

whose solution is x = y = 1 to thirteen digits. Without pivoting the multiplier is 10¹³; row 2 becomes (1 − 10¹³) y = 2 − 10¹³, so y = 1 − 10⁻¹³ correctly, but then x = (1 − y)/10⁻¹³ divides a number known only to 16 digits (1 − y ≈ 10⁻¹³ carries an absolute error of 10⁻¹⁶, i.e. 0.1 % relative) by 10⁻¹³: x comes out as 1.0003. Swap the rows first and the multiplier is 10⁻¹³, nothing is amplified, and x = 1.000000000000. The rule "largest pivot first" keeps every multiplier at most 1 in magnitude, which bounds the growth of rounding errors. A zero pivot with a non-zero entry below it (the matrix [[0, 1], [1, 0]]) is the extreme case: without a swap the algorithm divides by zero on a perfectly good system.

### 2.3 The cost: 2n³/3

Step k of elimination clears column k below the pivot: for each of the n − k rows below, one multiplier and then n − k + 1 multiply-subtract pairs across the row and the right-hand side. That is about 2(n − k)(n − k + 1) floating-point operations. Summed over k:

Σ_{k=1}^{n} 2(n − k)(n − k + 1) = 2 Σ_{j=0}^{n−1} j(j + 1) ≈ 2n³/3.   (3)

Back substitution adds n² (one dot product of average length n/2 per row). So the cost grows as the *cube* of the size: doubling n costs eight times as much. A 1000 × 1000 system is 7 × 10⁸ operations, a tenth of a second on one core; a 10⁶ × 10⁶ finite-element model would be 7 × 10¹⁷, never solved this way. Real solvers exploit sparsity (problem 7 shows a tridiagonal system solved in O(n)), and that is where the n³ intuition decides what is feasible.

### 2.4 Elimination is a factorisation: A = LU

Each elimination step "subtract ℓ × row k from row i" is multiplication on the left by an elementary matrix E_{ik} = I − ℓ e_i e_kᵀ (`@mit-ocw-18-06-linear-algebra` Lecture 2–4). Doing all the steps gives E A = U with U upper triangular. The inverse of each E is I + ℓ e_i e_kᵀ, and the product of the inverses, in the reverse order, is a unit lower triangular matrix L whose entry (i, k) is exactly the multiplier ℓ_{ik} used at that step. Hence

**A = L U**,   (4)

and with row swaps, **P A = L U** for a permutation matrix P. For the 3 × 3 example without pivoting, L has entries −3/2 and −1 below the diagonal in column 1 and 4 in column 2 (the multiplier of 2/0.5), and you can multiply L U back to A as a check.

Why this matters: elimination is the expensive part (n³); back substitution is cheap (n²). Once L and U are stored, a new right-hand side costs only two triangular solves: **L y = P b**, then **U x = y**. A structure analysed for 500 load cases is factored once. The determinant comes free as well: det A = det(Pᵀ) det L det U = (±1) × 1 × Π u_{kk}: **the product of the pivots, times −1 for each row swap**. A zero pivot that cannot be fixed by swapping means det A = 0 and the system has either no solution or infinitely many: a mechanism, in the truss.

The **inverse** is elimination applied to [A | I] all the way to [I | A⁻¹] (Gauss–Jordan, `@mit-ocw-18-06-linear-algebra` Lecture 3): three times the cost of a solve, and you should almost never form it; "x = A⁻¹ b" in a textbook means "solve A x = b" in code.

### 2.5 What a solve guarantees, and what it does not

Floating-point elimination with partial pivoting is **backward stable**: the computed x̂ is the exact solution of a nearby problem, (A + δA) x̂ = b with ‖δA‖/‖A‖ of order 10⁻¹⁶ (times a modest growth factor). Consequently the **residual** r = b − A x̂ is always tiny. But the **error** x̂ − x is not: perturbing A by δA perturbs x by up to κ(A) ‖δA‖/‖A‖ relative, where κ(A) = ‖A‖ ‖A⁻¹‖ is the **condition number** (`P01.18` defines it properly through singular values). The Hilbert matrix H_{ij} = 1/(i + j + 1) is the standard villain: for n = 10, κ ≈ 10¹³, and solving H x = H·**1** returns x with errors of 10⁻⁴ although the residual is 10⁻¹⁶. Nothing is wrong with the solver; the problem itself does not determine its answer to better than four digits in sixteen-digit arithmetic.

**Iterative refinement** (compute r = b − A x̂, solve A d = r, set x̂ ← x̂ + d, repeat) fixes a *poor solve* of a well-conditioned problem, for example an unpivoted elimination: two steps bring the 10⁻¹³-pivot example to full accuracy. It does not fix conditioning: for the Hilbert matrix, the refined answer is no better, because d is computed from a residual that is already at rounding level. Knowing which of the two situations you are in is the whole art, and the measurable sign is the condition number, not the residual.

### 2.6 Where the model stops being valid

- **Determinacy by count is not enough.** m + r = 2j makes A square; it does not make it invertible. Replace the diagonal 1-5 by a second member 3-5 (same count) and A is singular: the right bay is a four-bar mechanism and the left bay has a redundant member. Elimination hits an exact zero pivot at step 8 and the code raises. `P01.14` turns this into rank, null space and the counting formula.
- **Indeterminate trusses.** Add a member and A has more columns than rows: the forces depend on the stiffness of the members, not on equilibrium alone, and the system to solve is K u = f with the stiffness matrix (Phase 07). Elimination is still the solver; the matrix is different.
- **Large deflections and buckling.** Equilibrium on the undeformed geometry fails when deflections are comparable to member length, and a compression member buckles at the Euler load regardless of what this linear analysis says its force is. The −25 kN diagonal in problem 5 is the member to check first.
- **Conditioning.** A nearly collinear pair of members at a joint, or a very long slender bay, gives a nearly singular A: the forces are correct in exact arithmetic but hypersensitive to geometry; a 1 mm survey error in the joint positions can move a member force by tens of percent. Problem 6 asks for the sign.

## 3. Build the computational model

Open `starter/main.py`. Every routine copies its input (`np.array(..., copy=True)`) so that elimination never destroys the caller's matrix; `_as_system` (given) does the copying and the shape checks.

### Step 1: back substitution

```python
def back_substitution(U, c):
    U, c = _as_system(U, c)
    n = len(c)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        if U[i, i] == 0.0:
            raise ValueError(f"zero pivot on the diagonal at row {i}")
        x[i] = (c[i] - U[i, i + 1:] @ x[i + 1:]) / U[i, i]
    return x
```

Last row first. `U[i, i + 1:] @ x[i + 1:]` is the dot product of the already-known unknowns with their coefficients, one NumPy call instead of an inner Python loop: that is what keeps the n² step cheap. `forward_substitution` is the mirror image, first row first, with `L[i, :i] @ y[:i]`.

### Step 2: forward elimination with partial pivoting

```python
def forward_elimination(A, b, pivot=True):
    U, c = _as_system(A, b)
    n = len(c)
    perm = np.arange(n)
    scale = max(float(np.max(np.abs(U))), np.finfo(float).tiny)
    for k in range(n):
        if pivot:
            p = k + int(np.argmax(np.abs(U[k:, k])))
            if p != k:
                U[[k, p]] = U[[p, k]]
                c[[k, p]] = c[[p, k]]
                perm[[k, p]] = perm[[p, k]]
        if abs(U[k, k]) < SINGULAR_TOL * scale:
            raise ValueError("matrix is singular to working precision")
        for i in range(k + 1, n):
            factor = U[i, k] / U[k, k]
            U[i, k:] -= factor * U[k, k:]
            c[i] -= factor * c[k]
    return U, c, perm
```

`np.argmax(np.abs(U[k:, k]))` finds the largest remaining entry in column k; `U[[k, p]] = U[[p, k]]` swaps two rows with fancy indexing (the right-hand side is evaluated first, so no temporary is needed). The singularity test is *relative* to the matrix's largest entry: a pivot of 10⁻¹⁴ in a matrix of entries of order 1 is zero to working precision, while the same pivot in a matrix scaled by 10⁻¹⁴ is perfectly good. `U[i, k:] -= factor * U[k, k:]` updates the whole remaining row in one vectorised operation; only the entries from column k onward change, and those before it are already zero. `gaussian_elimination` is three lines: eliminate, then back-substitute.

### Step 3: LU

`lu_decompose` is the same loop storing each `factor` in `L[i, k]` instead of applying it to `c`. One subtlety: when rows k and p are swapped *after* some multipliers have been stored in columns 0..k−1 of L, those multipliers must be swapped too (`L[[k, p], :k] = L[[p, k], :k]`), because they belong to the rows, not to the positions. The permutation matrix is `np.eye(n)[perm]`: row i of P picks out row perm[i] of A. `lu_solve` is `back_substitution(U, forward_substitution(L, P @ b))`. `determinant` multiplies the diagonal of U by det P (±1). `inverse_gauss_jordan` eliminates above *and* below each pivot on [A | I] and normalises each pivot row to 1.

### Step 4: the truss

```python
def truss_equilibrium_matrix(nodes, members, supports):
    j = nodes.shape[0]
    m, r = len(members), len(supports)
    A = np.zeros((2 * j, m + r))
    for col, (i, k) in enumerate(members):
        d = nodes[k] - nodes[i]
        e = d / np.linalg.norm(d)
        A[2 * i:2 * i + 2, col] += e
        A[2 * k:2 * k + 2, col] -= e
    for col, (node, dof) in enumerate(supports):
        A[2 * node + dof, m + col] = 1.0
    return A
```

Equation (2) assembled column by column: row 2i is the x equation at node i, row 2i + 1 the y equation. `truss_load_vector` returns −**P**; `solve_truss` checks the matrix is square, eliminates, and splits **x** into forces and reactions. The sign convention is in the docstring and in the tests: a positive F_m is tension.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-13-linear-systems` from the repo root. The oracle throughout is `numpy.linalg`; the assertions are properties. In order of how much they would hurt if wrong:

1. **Triangular solves** agree with `np.linalg.solve` to 10⁻¹² and refuse a zero diagonal; inputs are not modified.
2. **Elimination** yields an upper-triangular U with a valid permutation and the same solution; `gaussian_elimination` matches NumPy to 10⁻⁹ for n up to 60 with residuals below 10⁻¹²; identity, scaling and superposition hold.
3. **Pivoting** solves [[0, 1], [1, 0]] and the 10⁻¹³-pivot system correctly, while `pivot=False` divides by zero on the first and is wrong by more than 10⁻⁶ on the second.
4. **Singular input raises**: an exactly singular 2 × 2, a rank-2 5 × 5 built from two outer products, a non-square A, a mismatched b.
5. **LU**: P A = L U to 10⁻¹²; L unit lower with |ℓ| ≤ 1; U upper; P a permutation; `lu_solve` equals elimination for three right-hand sides. **Determinant**: matches NumPy; det(AB) = det A det B; det(2A) = 2⁴ det A; a row swap flips the sign; a singular matrix gives exactly 0. **Inverse**: matches `np.linalg.inv` and A A⁻¹ = A⁻¹ A = I to 10⁻¹¹.
6. **Refinement and conditioning**: two refinement steps repair the unpivoted solve to 10⁻¹²; Hilbert n = 10 has residual below 10⁻¹² yet error above 10⁻⁷, before and after refinement.
7. **Cost**: `elimination_flops(800)/elimination_flops(400)` is between 7.5 and 8.5.
8. **Truss**: A is 12 × 12; every member column sums to zero in x and in y and has norm √2 (Newton's third law built in); the reactions are (0, 7.5, 7.5) kN by global force and moment balance; the joint residual is below 10⁻¹²; the diagonal at node 0 carries −7.5/0.6 = −12.5 kN; forces double when the load doubles; node 3's two members are zero-force; a missing or misplaced member raises.

### Compare

`main()` prints the comparison. On random systems the hand-written solver and `numpy.linalg.solve` agree to 10⁻¹⁴ for n = 50, P A − L U is at 10⁻¹⁵, and the determinant and inverse are at 10⁻¹⁵ relative. The two solvers are doing the *same arithmetic* (LAPACK's `dgesv` is elimination with partial pivoting), so agreement to rounding is expected and a discrepancy of 10⁻⁸ would mean a bug, not a difference of method. Where they differ is cost: 10 ms against 0.7 ms at n = 100, 160 ms against 2.4 ms at n = 400, on the author's machine. The Python loop over k is the outer loop only (the inner work is vectorised), so our code is within a factor of 50–70 of LAPACK, which is blocked, cache-aware and multithreaded; a pure-Python triple loop would be a thousand times slower again. Use the library. Keep the hand-written one for understanding its error messages.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No measurement is involved; we are characterising the solver.

**Experiment A: cost.** Prediction: time ∝ n³, so doubling n multiplies the time by 8. Time `gaussian_elimination` and `np.linalg.solve` on random systems with n = 50, 100, 200, 400, 800 (use `time.perf_counter` and the best of three runs). Fit the exponent of each with `np.polyfit(np.log(n), np.log(t), 1)`. Expect close to 3 for yours at large n and a smaller exponent for NumPy at small n, where fixed overheads dominate. Also record `elimination_flops(n)` divided by the time: that is your solver's speed in flop/s, to compare with the ~10¹⁰ flop/s a single core achieves with LAPACK.

**Experiment B: pivoting and refinement.** For ε = 10⁻⁸, 10⁻¹⁰, 10⁻¹², 10⁻¹⁴, solve [[ε, 1], [1, 1]] x = (1, 2) with and without pivoting and record the error in x₁. Prediction: the unpivoted error grows as 10⁻¹⁶/ε, the pivoted one stays at 10⁻¹⁶. Then apply one, two and three steps of `iterative_refinement` with `pivot=False` and record the error after each.

**Experiment C: conditioning.** For n = 4, 6, 8, 10, solve H x = H·**1** and record κ(H) (`np.linalg.cond`), the residual and the maximum error. Prediction: error ≈ κ × 10⁻¹⁶, residual ≈ 10⁻¹⁶ throughout. Then do the same for the truss matrix after moving node 4 to (2, 0.015) m, so the diagonals become nearly collinear with the chords: record κ and the change in member forces for a 1 mm further move. Prediction: κ grows from 7.4 to several hundred, the chord forces grow to about 100 times the load, and 1 mm changes them by several percent.

Raw data format: three plain-text tables, (n, t_ours, t_numpy), (ε, err_nopivot, err_pivot, err_refined_1, err_refined_2), (n, κ, residual, error).

## 5. Validate and improve

Put the sources side by side for the truss with 15 kN at node 1:

| quantity | hand (global equilibrium, joint 0) | code | `numpy.linalg.solve` |
|---|---|---|---|
| R0x, R0y, R2y | 0, 7.5, 7.5 kN | 0.000, 7.500, 7.500 | same to 1.8 × 10⁻¹² N |
| diagonal 0-4 | −7.5/0.6 = −12.5 kN (compression) | −12.500 | same |
| bottom chord 0-1 | +0.8 × 12.5 = +10.0 kN (tension) | +10.000 | same |
| members 0-3, 3-4, 1-2 | 0 (unloaded two-member joints) | 0.000 | same |
| relative residual | — | 5.2 × 10⁻¹⁷ | — |
| det A | — | −0.960000 | −0.960000 |

Every difference must be explained, and here each has a cause you can name:

- **1.8 × 10⁻¹² N between ours and NumPy** on forces of 10⁴ N is 10⁻¹⁶ relative: both solvers round differently in the last bit because LAPACK blocks the operations in a different order. Not a discrepancy.
- **The reactions are symmetric although the truss is not.** Global equilibrium (three equations: ΣF_x, ΣF_y, ΣM) fixes the three reactions without any knowledge of the members: that is why they are the first thing to check by hand, and why a wrong reaction means the *matrix* is wrong (a sign in a reaction column or in **b**), never the elimination.
- **Member 1-2 is zero while 0-1 carries 10 kN.** Node 2 has the roller reaction 7.5 kN up, the vertical 2-5 and the chord 1-2. The vertical must take all of the 7.5 kN (compression), and then the x equation at node 2 forces F_{12} = 0. The load path goes 1 → 5 (diagonal) → 2 (vertical), not along the chord. The asymmetric diagonal layout makes the asymmetric force pattern; the plot shows it.
- **det A = −0.96** has no units (the entries are direction cosines and ones) and its magnitude by itself says little; a crude conditioning check compares it with the product of the column norms (√2⁹ × 1³ ≈ 22.6), and `np.linalg.cond(A)` = 7.4 says the truss is well posed. Moving node 4 down to y = 15 mm, so that the diagonals are nearly collinear with the chords, gives det −0.0077, κ = 606, a bottom-chord force of 1000 kN for the same 15 kN load, and a further 1 mm move changes that force by 6 %: the geometry, not the solver, is the problem.

What is still unverified when everything agrees: the model matches itself. The real stand has welded or bolted joints that carry moment, members with self-weight, and an engine that is not a point load; the truss idealisation is good to a few percent for slender members and is exactly how preliminary sizing is done, but the final check is a strain gauge on the diagonal (Phase 07).

Improvements worth making: return the condition number estimate alongside the solution (`np.linalg.cond` is itself an SVD, `P01.18`); store L and U in place in one array as LAPACK does; let `solve_truss` accept a list of load cases and factor once (`lu_decompose` then `lu_solve` per case). Improvements worth not making: forming A⁻¹ to "solve faster for many loads": it costs three solves to build and is less accurate than LU for each use.

## 6. Defend your work

Submit: `main.py` with all tests passing; `truss_forces.png`; the hand solution of problem 1 with the pivot sequence; the three tables of part 4 with fitted exponents; and a one-page note (`projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the truss assumptions, the sign convention, and your rule for when to trust a solve.

### Adversarial self-checks

- **What if the load doubles?** Every force and reaction doubles exactly (the system is linear; a test asserts it). What if the load is moved to node 4? The reactions do not change (same x position) but the load path does; re-solve and see which members switch sign.
- **What if the truss height goes to zero?** Nodes 3, 4, 5 fall onto 0, 1, 2; the diagonals become collinear with the chords; det A → 0 and the forces → ∞ before that. A truss needs depth; the forces scale as (span/height).
- **Which assumption fails first in the real world?** Pin joints: welded joints carry secondary bending moments of a few percent of the axial stress; next, the point load: an engine mount spreads over a bracket; last, small deflections, which only fail near buckling.
- **What measurement could falsify the model?** A strain gauge on member 0-4 reading −12.5 kN/(E A) within 5 %; a reading 30 % off means joint fixity or a load not at the joint.
- **Does the result survive a change of solver, of node numbering?** It must: NumPy agrees to rounding, and renumbering the nodes (swap the labels of nodes 1 and 4, say) permutes rows and columns of A without changing any member force, to 10⁻¹² N. Elimination *without* pivoting, on the other hand, fails on this matrix at step 1: the y equation of node 0 has no entry for member 1-2, so the second pivot is an exact zero that only a row swap can fix. A zero pivot is not the same thing as a singular matrix, and the code must tell them apart.
- **Which decision changes if the uncertainty doubled?** With member forces known to ±2 % instead of ±1 %, the margin on the compression diagonal against buckling must grow by the same amount; if the stand was sized at a safety factor of 1.5 on a buckling load known to ±10 %, the force uncertainty is not the driver, the buckling model is.

## Check yourself

1. In the equilibrium matrix, what does a column with entries (+0.8, +0.6) in rows 0–1 and (−0.8, −0.6) in rows 8–9 represent?
2. Why does each member column sum to zero in both x and y, and what physical law is that?
3. Elimination without pivoting on [[0, 1], [1, 0]] fails. On [[10⁻¹³, 1], [1, 1]] it "succeeds". What goes wrong in the second case?
4. How many operations does elimination take for n = 2000 compared with n = 1000, and roughly how long at 10⁹ flop/s?
5. State what L, U and P are in P A = L U, and where the multipliers end up.
6. A solve returns a residual of 10⁻¹⁶ and the condition number is 10¹². How many digits of the answer can you trust?
7. The truss has m + r = 2j. Name one way it can still be a mechanism, and what the code does in that case.
8. You need the forces for 500 load cases on the same truss. What do you factor, how many times, and what does each additional case cost?
