# P01.14  Vector Spaces, Basis & Rank

> Counting unknowns against equations tells you whether a structure *might* stand up. Rank tells you whether it *does*.

| | |
|---|---|
| **Learning outcomes** | 1. Define subspace, span, linear independence, basis and dimension, and compute the four fundamental subspaces of a small matrix by hand from its reduced row echelon form. 2. Implement reduced row echelon form, rank, null space, left null space, column and row space bases and the general solution of A x = b from scratch in NumPy, and test them against `numpy.linalg.matrix_rank` and an SVD null space. 3. Derive Maxwell's counting rule m + r − 2j = s − k for a planar truss from rank–nullity, and explain why the count can be zero for a structure that is still a mechanism. 4. Detect a mechanism in a truss from the left null space of its equilibrium matrix, display the motion, and verify that it stretches no member. 5. Decide what a numerical rank means when the matrix is nearly singular, and say which tool answers the question rank cannot. |
| **Prerequisites** | `P01.13` (elimination, the truss equilibrium matrix, `truss_equilibrium_matrix` is reused unchanged), `P01.07` (vectors, dot product), `P00.03` (NumPy), `P00.07` (pytest). Skill check: you can row-reduce a 3 × 4 matrix by hand and read a solution off the result. |
| **Engineering question** | Four two-bay trusses have the same six joints and the same three supports; two of them have exactly as many unknowns as equations. Which of them can move without any member stretching, which carry a locked-in force with no load at all, and how do we tell from the equilibrium matrix alone? |
| **Assumptions and limits** | Pin-jointed planar truss, small displacements: a "mechanism" here is a first-order (infinitesimal) mechanism, a displacement with zero member elongation to first order; some such structures stiffen at second order (cable nets, tensegrity) and are usable when pre-stressed. Rank is computed in floating point with a relative tolerance; a truss whose members are nearly collinear has full rank at every reasonable tolerance yet is practically a mechanism, and the rref pivots are a poor measure of that (the SVD of `P01.18` is the right one). |
| **Analytical work** | rref of a 3 × 4 matrix by hand with its four subspaces; proof that the special solutions are independent and span the null space; rank–nullity for A and Aᵀ; Maxwell's rule from rank–nullity; the mechanism and self-stress of the Maxwell-trap truss found by inspection and confirmed by the code. |
| **Computational deliverable** | `starter/main.py`: `rref`, `rank`, `nullspace`, `left_nullspace`, `column_space_basis`, `row_space_basis`, `is_independent`, `four_subspace_dimensions`, `is_consistent`, `solve_general`, `truss_mechanisms`, `truss_self_stress`, `member_elongations`, `classify_truss`. Tests prove: rref of a textbook matrix; pivot columns are unit vectors and the form is idempotent; rank equals `numpy.linalg.matrix_rank` on matrices of known rank and is invariant under transpose, scaling and Gram products; the null space is annihilated, has dimension n − r and spans the same subspace as the SVD null space; rank–nullity both ways; the general solution reproduces b for any null-space combination and inconsistency raises; the four trusses classify as determinate, mechanism, indeterminate and both; the mechanism mode stretches no member and moves no support; the self-stress state lives in the doubled bay with the bay's cos/sin force ratio. |
| **Practical activity** | **computational experiment**: a tolerance study of numerical rank on a truss driven toward a mechanism, and a search for the single member that repairs a mechanism. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `truss_mechanism.png`, the hand rref with its four subspaces, the tolerance table, and a short note on what rank can and cannot tell a structural engineer. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4. Pass: Maxwell's rule derived from rank–nullity with the role of each subspace named; all public tests pass; the Maxwell-trap mechanism drawn and explained; the judgment problem answered with a measurable sign. |
| **References** | `@mit-ocw-18-06-linear-algebra` Lecture 5: Transposes, permutations, spaces R^n (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 6: Column space and nullspace (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 7: Solving Ax = 0: pivot variables, special solutions (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 8: Solving Ax = b: row reduced form R (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 9: Independence, basis, and dimension (required-free); `@mit-ocw-18-06-linear-algebra` Lecture 10: The four fundamental subspaces (required-free). |
| **Extensions** | ★ Extend the truss functions to three dimensions (three rows per joint) and classify a tetrahedral space frame; ★ compute the smallest singular value of the equilibrium matrix with `np.linalg.svd` for the near-mechanism sweep and compare with the rref pivots; ★ implement the parity-space fault detector of problem 7 for a 5-gauge balance and show it isolates any single faulty gauge. |
| **Time** | unmeasured |

## 1. The engineering problem

In `P01.13` the two-bay truss had 9 members and 3 reactions against 12 joint equations, and elimination solved it. A designer who adds or removes a member will first check the count: m + r against 2j. Maxwell's rule says that if m + r < 2j the truss must be a mechanism, and if m + r > 2j it must be statically indeterminate (the forces depend on stiffness, not on equilibrium alone). What the rule does *not* say is that m + r = 2j makes the truss rigid. Put both diagonals in the left bay and none in the right: the count is still 12 = 12, and the right bay folds like a parallelogram while the left bay carries a self-balancing force that no load produces. This "Maxwell trap" is a classic design error, and the count cannot see it.

The equilibrium matrix **A** of `P01.13` can. The question "can the truss move without any member stretching?" is the question "is there a joint displacement **d** ≠ 0 with **Aᵀ d = 0**?", and the question "can the members carry force with no load?" is "is there **x** ≠ 0 with **A x = 0**?". Both are questions about the **null spaces** of **A** and **Aᵀ**, and their dimensions are fixed by one number, the **rank** of **A**, through the rank–nullity theorem. Maxwell's count turns out to be the rank–nullity theorem with the rank left out. This lesson supplies the vocabulary (subspace, span, independence, basis, dimension, the four fundamental subspaces of `@mit-ocw-18-06-linear-algebra` Lecture 10), the algorithm (reduced row echelon form), and the engineering use (mechanism detection), then asks what the numbers mean when the matrix is only nearly singular.

Success criteria, written before we start: `rank` agrees with `numpy.linalg.matrix_rank` on random matrices of every rank from 0 to full, for tall, wide and square shapes; the computed null space is annihilated by **A** to 10⁻¹⁰ and spans the same subspace as the SVD null space (projectors equal to 10⁻¹⁰); the four trusses classify as determinate (rank 12), mechanism (rank 11, k = 1), indeterminate (rank 12, s = 1) and Maxwell trap (rank 11, k = 1, s = 1), with the count m + r − 2j equal to s − k in every case; the mechanism mode stretches no member to 10⁻¹⁰ and the self-stress state balances every joint with no load.

## 2. Understand and derive

### 2.1 Subspaces, span, independence, basis, dimension

A **subspace** of Rⁿ is a set of vectors closed under addition and scaling: it contains **0**, and with **u** and **v** it contains every a**u** + b**v** (`@mit-ocw-18-06-linear-algebra` Lecture 5). Lines and planes through the origin are subspaces; a line that misses the origin is not. The **span** of a set of vectors is the subspace of all their linear combinations. Vectors are **linearly independent** if the only combination giving **0** is the trivial one: c₁**v**₁ + ⋯ + c_k**v**_k = **0** forces every c_i = 0. A **basis** of a subspace is an independent set that spans it; every basis of a given subspace has the same number of vectors, its **dimension** (Lecture 9). In R³, three vectors are a basis exactly when the 3 × 3 matrix with those columns is invertible.

In the truss, a joint displacement vector **d** ∈ R¹² and a force vector **x** ∈ R^{m+r} live in different spaces, and **A** maps forces to net joint loads. Each column of **A** is the pattern of joint loads that one unit member force (or one unit reaction) produces. The question "can this load be carried?" is "is −**P** in the span of the columns?".

### 2.2 The four fundamental subspaces

For an m × n matrix **A** (`@mit-ocw-18-06-linear-algebra` Lectures 6 and 10):

| subspace | lives in | definition | dimension | truss meaning (A is 2j × (m+r)) |
|---|---|---|---|---|
| column space C(A) | R^m | span of the columns = {A x} | r | loads the members can balance |
| null space N(A) | R^n | {x : A x = 0} | n − r | self-stress: forces with no load |
| row space C(Aᵀ) | R^n | span of the rows | r | force patterns that *do* produce load |
| left null space N(Aᵀ) | R^m | {y : Aᵀ y = 0} | m − r | mechanisms: motions no member resists |

The number r is the **rank**, the dimension of the column space, and the dimension counts are the **rank–nullity theorem**: dim N(A) = n − r and dim N(Aᵀ) = m − r. Two orthogonality facts complete the picture and are checked by the tests: every vector in N(A) is orthogonal to every row of **A** (that is what A x = 0 says row by row), so N(A) ⟂ C(Aᵀ); likewise N(Aᵀ) ⟂ C(A).

### 2.3 Reduced row echelon form: the algorithm that exposes all four

Elimination (`P01.13`) taken all the way, with each pivot scaled to 1 and its column cleared above as well as below, gives the **reduced row echelon form** R = rref(A) (`@mit-ocw-18-06-linear-algebra` Lectures 7–8). Row operations do not change the row space or the null space, so R has the same rank and null space as **A**, and in R everything is visible.

**Worked example** (Strang's, Lecture 7):

```
A = [ 1  2  2  2 ]        rref       R = [ 1  2  0 -2 ]     pivot columns 1 and 3
    [ 2  4  6  8 ]        ---->          [ 0  0  1  2 ]     free columns 2 and 4
    [ 3  6  8 10 ]                       [ 0  0  0  0 ]     rank r = 2
```

Steps: row 2 ← row 2 − 2 row 1 = (0, 0, 2, 4); row 3 ← row 3 − 3 row 1 = (0, 0, 2, 4); column 2 has no pivot (all zeros below row 1), move on; pivot 2 in column 3 of row 2, scale to (0, 0, 1, 2); clear it from row 1 (row 1 − 2 row 2 = (1, 2, 0, −2)) and row 3 (zero row). Pivot columns {1, 3}; free columns {2, 4}; rank 2.

- **Null space.** Set one free variable to 1, the others to 0, and read the pivot variables from R: free x₂ = 1 gives x₁ = −2, x₃ = 0: **s**₁ = (−2, 1, 0, 0); free x₄ = 1 gives x₁ = 2, x₃ = −2: **s**₂ = (2, 0, −2, 1). These **special solutions** are independent (look at the free rows: they form an identity block) and every null vector is a combination of them (its free components fix the combination, and the pivot components then agree): a basis, dimension n − r = 2. Check: A **s**₁ = (−2 + 2, −4 + 4, −6 + 6) = 0. ✓
- **Column space.** The pivot columns of **A** itself, (1, 2, 3) and (2, 6, 8), are a basis (the columns of R are *not*; row operations change the column space). Dimension 2.
- **Row space.** The two non-zero rows of R. Dimension 2.
- **Left null space.** Solve Aᵀ y = 0, or track the row operations: here row 3 − row 2 − row 1 = 0, so **y** = (−1, −1, 1) up to scale. Dimension m − r = 1.
- **Solving A x = b.** rref of the augmented [A | b]: if a pivot appears in the last column, the system is inconsistent (b ∉ C(A)); otherwise a particular solution has the free variables 0 and the pivot variables read from the last column, and every solution is **x**_p + N(A) (Lecture 8).

### 2.4 The truss, and Maxwell's rule as rank–nullity

**A** is 2j × (m + r). Apply rank–nullity twice:

dim N(A) = (m + r) − rank =: s, the number of independent **states of self-stress**;
dim N(Aᵀ) = 2j − rank =: k, the number of independent **mechanisms**.

Subtract: **m + r − 2j = s − k**.   (1)

That is Maxwell's rule, with its missing half restored: the count equals the *difference* between self-stress states and mechanisms, not either one. A count of zero says s = k, which is satisfied by s = k = 0 (determinate) and equally by s = k = 1 (the trap). A negative count forces k ≥ 1; a positive count forces s ≥ 1; nothing else follows from counting.

**Why Aᵀ d = 0 means "no member stretches".** The column of member (i, k) is +**e** at node i and −**e** at node k. Its dot product with a displacement vector **d** is **e** · (**d**_i − **d**_k) = −(elongation of the member to first order). A reaction column dotted with **d** is the displacement of that restrained degree of freedom. So Aᵀ **d** = 0 says every member has zero first-order elongation *and* every support stays put, which is the definition of an (infinitesimal) mechanism. This is the **principle of virtual work** in matrix form: the equilibrium matrix and the compatibility matrix are transposes of each other.

**Worked example: the Maxwell trap by inspection.** Members: chords 0-1, 1-2, 3-4, 4-5; verticals 0-3, 1-4, 2-5; diagonals 0-4 and 1-3, both in the left bay. Count: 9 + 3 = 12 = 2j. The right bay 1-2-5-4 is a quadrilateral with no diagonal, so it can shear. The code finds the mode **d** = (0, 0; 0, −1; 0, 0; 0.75, 0; 0.75, −1; 0.75, 0) (one (x, y) pair per node): joints 1 and 4 drop by 1, the whole top chord shifts right by 0.75, the supported joints 0 and 2 stay put. Check each member with **e** · (**d**_k − **d**_i): chord 0-1, (1, 0) · (0, −1) = 0; chord 1-2, (1, 0) · (0, 1) = 0; vertical 1-4, (0, 1) · (0.75, 0) = 0; top chords, (1, 0) · (0, ∓1) = 0; verticals 0-3 and 2-5, (0, 1) · (0.75, 0) = 0; diagonal 0-4 with **e** = (0.8, 0.6), (0.8, 0.6) · (0.75, −1) = 0.6 − 0.6 = 0; diagonal 1-3 with **e** = (−0.8, 0.6), (−0.8, 0.6) · (0.75, 1) = −0.6 + 0.6 = 0. Every elongation is zero to first order; the motion is real (a vertical displacement of a horizontal bar's end changes its length only at second order). The self-stress state is the left bay's two diagonals in tension (+1 each) with its two chords at −0.8 and two verticals at −0.6 in compression: at node 0 the diagonal pulls (0.8, 0.6), the chord pushes back −0.8 in x and the vertical −0.6 in y, so the joint balances with no reaction; the same happens at the other three corners, and the right bay carries nothing.

### 2.5 Where the model stops being valid

- **First order only.** A mechanism of **Aᵀ** is a motion with zero elongation *to first order*. The Maxwell trap's right bay really folds (a finite mechanism). A square of four bars with a diagonal cable, or a cable net, has first-order mechanisms that stiffen at second order and are perfectly usable once pre-stressed (tensegrity, sails, tents). Rank says "check this"; it does not by itself say "unsafe".
- **Numerical rank is a decision, not a fact.** In floating point no pivot is ever exactly zero unless the geometry is. `rank` treats a pivot below `tol × max|A|` as zero. Move node 4 of the determinate truss from y = 1.5 m toward the bottom chord: the rank stays 12 at every tolerance down to 10⁻¹⁰ even at y = 1.5 µm, where the smallest singular value is 3.6 × 10⁻⁷ and the forces for any load are 10⁶ times the load. rref pivots are a poor proxy for near-singularity because scaling a row to 1 hides how small it was; the honest measure is the smallest singular value (`P01.18`), and part 4 shows the two side by side.
- **Geometry errors.** A real truss is built to tolerances; a mechanism that exists only for one exact geometry (members exactly collinear) becomes, after a millimetre of error, a nearly singular but technically rigid structure with enormous forces. The design rule is to stay far from rank deficiency, not merely off it.
- **Beyond trusses.** The same rank test answers "can three reaction wheels produce torque about every axis?" (problem 3) and "which sensor failed?" (problem 7). Rank is a property of a matrix, not of a structure.

## 3. Build the computational model

Open `starter/main.py`. `truss_equilibrium_matrix`, `truss_variant` and `maxwell_count` are given; everything else is yours. All bases are returned as matrices whose *columns* are the basis vectors, so that `A @ N` is a one-line check.

### Step 1: rref with a relative tolerance

```python
def rref(A, tol=DEFAULT_TOL):
    R = np.array(A, dtype=float, copy=True)
    m, n = R.shape
    scale = float(np.max(np.abs(R)))
    pivots, row = [], 0
    for col in range(n):
        if row >= m:
            break
        p = row + int(np.argmax(np.abs(R[row:, col])))
        if abs(R[p, col]) <= tol * scale:
            R[row:, col] = 0.0
            continue
        if p != row:
            R[[row, p]] = R[[p, row]]
        R[row] /= R[row, col]
        for i in range(m):
            if i != row and R[i, col] != 0.0:
                R[i] -= R[i, col] * R[row]
        pivots.append(col)
        row += 1
    R[np.abs(R) <= tol * scale] = 0.0
    return R, pivots
```

Differences from `P01.13`: the loop is over *columns* and the pivot row index advances only when a pivot is found, so a column without a pivot is skipped and recorded as free; the pivot row is scaled to 1 and cleared *above* as well as below (Gauss–Jordan); and the zero test is `tol * scale`, relative to the largest entry of the original matrix. The final sweep removes rounding dust so that the zero rows are exactly zero. `rank` is `len(pivots)`.

### Step 2: the null space from the special solutions

```python
def nullspace(A, tol=DEFAULT_TOL):
    R, pivots = rref(A, tol)
    n = R.shape[1]
    free = [c for c in range(n) if c not in pivots]
    N = np.zeros((n, len(free)))
    for j, f in enumerate(free):
        N[f, j] = 1.0
        for i, p in enumerate(pivots):
            N[p, j] = -R[i, f]
    return N
```

Exactly §2.3: column j of N is the special solution for free variable f, with a 1 in row f and −R[i, f] in the row of the i-th pivot variable (from the i-th reduced equation x_p + R[i, f] x_f = 0). With no free columns, N has shape (n, 0), which `A @ N` and `np.linalg.matrix_rank` both accept. `left_nullspace(A)` is `nullspace(A.T)`.

### Step 3: the other bases, consistency and the general solution

`column_space_basis` returns `A[:, pivots]` (the pivot columns of **A**, not of R). `row_space_basis` returns the non-zero rows of R, transposed to columns. `is_independent(V)` is `rank(V) == V.shape[1]`. `four_subspace_dimensions` packages r, n − r and m − r. `is_consistent(A, b)` compares `rank([A | b])` with `rank(A)`. `solve_general` row-reduces the augmented matrix, raises if the last column is a pivot column, otherwise sets the pivot variables from the last column and returns the particular solution with the null-space basis.

### Step 4: the truss

```python
def truss_mechanisms(A, tol=DEFAULT_TOL):
    return left_nullspace(A, tol)

def truss_self_stress(A, tol=DEFAULT_TOL):
    return nullspace(A, tol)

def classify_truss(nodes, members, supports, tol=DEFAULT_TOL):
    A = truss_equilibrium_matrix(nodes, members, supports)
    j, m, r = len(nodes), len(members), len(supports)
    rk = rank(A, tol)
    k, s = 2 * j - rk, m + r - rk
    ...
```

Two lines each, because the work was done in §2.4: a mechanism is a left null vector, a self-stress state is a null vector, and the counts are rank–nullity. `member_elongations(nodes, members, d)` computes **e** · (**d**_k − **d**_i) for each member directly from the geometry, independently of **A**, so that the test "a mechanism stretches nothing" does not merely re-check Aᵀ d = 0 with the same matrix.

### Verify

Run `python -m pytest curriculum/phase-01/lesson-14-vector-spaces-basis-rank` from the repo root. The oracles are `np.linalg.matrix_rank` and `np.linalg.svd`. In order of how much they would hurt if wrong:

1. **rref by hand.** Strang's 3 × 4 matrix gives the R of §2.3 with pivots [0, 2] and leaves the input untouched; the identity is its own rref; the zero matrix has no pivots; pivot columns of R are unit vectors; rows beyond the rank are zero; rref(R) = R.
2. **Rank** equals `matrix_rank` on seven (m, n, r) shapes including rank 0 and full rank; rank is invariant under transpose, scaling by 10^{±6}, forming AᵀA and AAᵀ, and repeating columns, and rises by one when a generic column is added.
3. **Tolerance.** Two columns that differ from dependence by 10⁻¹³ count as dependent at tol = 10⁻¹⁰ and as independent at 10⁻¹⁶.
4. **Null space.** A N = 0 to 10⁻¹⁰, N has n − r independent columns, the projector N N⁺ equals the SVD projector; the free rows of N are an identity block; the left null space annihilates from the left.
5. **Bases.** The column-space basis consists of actual columns of A, is independent and spans every column; the row-space basis is orthogonal to the null space; `is_independent` on hand cases.
6. **Rank–nullity** both ways on three shapes.
7. **A x = b.** A consistent b gives a particular solution plus a 5-dimensional null space, every combination of which reproduces b; a generic perturbation of b is inconsistent and raises; a square invertible system gives the unique solution with an empty N.
8. **The trusses.** The four variants give (rank, k, s) = (12, 0, 0), (11, 1, 0), (12, 0, 1), (11, 1, 1) and the right verdicts, rank agreeing with NumPy; Maxwell's count equals s − k; the mechanism mode of the two defective trusses is annihilated by Aᵀ, stretches no member (by `member_elongations`), moves no support, and moves joint 1; the determinate truss has no mechanism; the self-stress state of the two doubled trusses balances every joint, needs no reactions, is confined to the left bay, and has both diagonals of one sign with the chords at −0.8 times the diagonal (the bay's cosine); elongations are linear in **d**, vanish for a rigid translation and equal −A_membersᵀ **d**.

### Compare

`np.linalg.matrix_rank` does not row-reduce: it computes the singular values and counts those above a threshold (by default max(m, n) × σ_max × ε). On every test matrix it agrees with our rref rank, because those matrices are either exactly rank-deficient or comfortably full rank. The two methods disagree on *nearly* singular matrices, and that disagreement is the content of part 4. The SVD null space (`vt[r:].T` from `np.linalg.svd`) is orthonormal; ours is not (special solutions have an identity block in the free rows and whatever the reduced rows dictate elsewhere), so the comparison is between the subspaces they span, via the projectors N N⁺ and N_svd N_svdᵀ: they agree to 4 × 10⁻¹⁵ on the Maxwell trap. In practice you would use `scipy.linalg.null_space` (SVD-based, orthonormal, with a sensible tolerance) and keep rref for small matrices you want to read by eye and for the special solutions, which are the ones with a structural meaning (one free member force at a time).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No measurement is involved; we are probing what numerical rank means.

**Experiment A: a truss driven toward a mechanism.** Take the determinate truss and move node 4 from y = 1.5 m to 0.15, 0.015, …, 1.5 × 10⁻⁶ m, so that the diagonals become collinear with the chords. For each geometry record: `rank(A)` at tol = 10⁻¹⁰, 10⁻⁶ and 10⁻³; `np.linalg.matrix_rank(A)`; the smallest singular value σ_min from `np.linalg.svd(A, compute_uv=False)`; and the condition number σ_max/σ_min. Prediction before running: σ_min falls in proportion to y₄ (the diagonal's vertical component), the condition number rises as 1/y₄, and the rref rank stays at 12 until the tolerance exceeds roughly y₄/1.5, i.e. the drop happens at different y₄ for each tolerance. Then compute the member forces for the 15 kN load at each geometry and watch them grow as 1/y₄.

**Experiment B: repair a mechanism.** Remove one member from the determinate truss (use member index s mod 9 with your seed) and classify. Then try adding back every possible member between two unconnected joints, one at a time, and classify each. Prediction: only members that connect the moving part of the mechanism to the fixed part restore rank 12; members inside the already-rigid part add a self-stress state and leave the mechanism.

Raw data format: two plain-text tables, (y₄, rank at three tolerances, numpy rank, σ_min, κ, max member force) and (candidate member, rank, k, s, verdict).

## 5. Validate and improve

Put the sources side by side for the four trusses:

| truss | m | m + r − 2j (hand) | rank (ours) | rank (NumPy) | k, s (ours) | k, s (inspection, §2.4) |
|---|---|---|---|---|---|---|
| determinate | 9 | 0 | 12 | 12 | 0, 0 | rigid, no redundancy |
| missing diagonal | 8 | −1 | 11 | 11 | 1, 0 | right bay folds |
| extra diagonal | 10 | +1 | 12 | 12 | 0, 1 | left bay doubled |
| Maxwell trap | 9 | 0 | 11 | 11 | 1, 1 | right bay folds *and* left bay doubled |

Every entry must be explained, and here each has a cause you can name:

- **Rank agrees with NumPy in every row** because each matrix is either exactly deficient (a member column is an exact combination of others: the geometry is in whole numbers and 0.6/0.8) or full rank with σ_min = 0.29. Agreement here says nothing about nearly singular cases.
- **The count is right in every row and wrong as a verdict in one.** Rows 1 and 4 have the same count and opposite verdicts; the rank distinguishes them. Equation (1) holds in all four rows: 0 = 0 − 0, −1 = 0 − 1, +1 = 1 − 0, 0 = 1 − 1.
- **The mechanism mode has zero elongation to 10⁻¹⁷**, not 10⁻¹⁰: the special solution is built from an rref whose entries are exact decimals (0.75, 1) for this geometry, so the check is exact up to a few rounding operations. On a truss with irrational direction cosines expect 10⁻¹⁵.
- **The self-stress forces are 1, 1, −0.8, −0.8, −0.6, −0.6**: the bay's geometry (2 m by 1.5 m, diagonal 2.5 m) gives cos = 0.8, sin = 0.6, and a unit tension in each diagonal is balanced at every corner by a compression of 0.8 in the chord and 0.6 in the vertical. That this state exists at all means the left bay could be "wound up" by shortening one diagonal in manufacture, without any external load: the structure would carry locked-in stress forever, invisible to a load test.

From Experiment A (author's numbers): at y₄ = 1.5 × 10⁻⁴ m the rref rank is 12 at tol 10⁻¹⁰ and 10⁻⁶ but 11 at tol 10⁻³; `matrix_rank` says 12; σ_min = 3.6 × 10⁻⁵ and κ = 6 × 10⁴; the chord force for 15 kN is 10⁵ kN. Every tool reports "full rank" while the structure is useless. The lesson: numerical rank answers "is there an exact dependency within tolerance?"; conditioning answers "how close is the nearest dependency?"; a designer needs the second.

What is still unverified when everything agrees: a first-order mechanism found here may or may not be a finite one, and a structure flagged "indeterminate" may be the better design (redundancy is a safety feature in aircraft structures; what rank tells you is that you now need stiffness, not just equilibrium, to find the forces). Rank is the beginning of the structural analysis, not its verdict.

Improvements worth making: return σ_min alongside the rank in `classify_truss`; make `nullspace` orthonormalise its basis with Gram–Schmidt (`P01.17`) so that mode amplitudes are comparable. Improvements worth not making: lowering the tolerance to "catch" near-mechanisms; a tolerance is a decision about exact dependency and cannot be tuned into a conditioning measure.

## 6. Defend your work

Submit: `main.py` with all tests passing; `truss_mechanism.png`; your hand rref of §2.3 with all four subspaces; the two tables of part 4; and a one-page note (`projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating what the rank test assumes, what it found for the four trusses, and what it cannot find.

### Adversarial self-checks

- **What if the loads double?** Nothing in this lesson changes: rank, null spaces and mechanisms are properties of **A**, which depends on geometry and connectivity only. What if the truss is scaled by 2? The direction cosines are unchanged, so **A** is unchanged.
- **What if a support is removed?** r drops by one: the count drops by one and, unless a self-stress state is also lost, k rises by one: a rigid-body motion appears (the roller's removal lets the whole truss rotate about the pin). The left null space now contains that rigid-body mode; check that it stretches nothing.
- **Which assumption fails first in the real world?** First-order elongation: a mechanism with finite motion changes its own geometry and either locks (second-order stiffening) or collapses; the linear analysis cannot say which. Next, pin joints: welded joints make every planar truss a frame with no mechanism in the truss sense but possibly a very flexible one.
- **What measurement could falsify the model?** Push the Maxwell-trap truss at joint 1 with a small lateral force: the model predicts a finite motion at zero stiffness; a real model built with slightly stiff joints shows a small but finite stiffness, and its magnitude measures how far the "mechanism" is from reality. Conversely, cutting the diagonal 1-3 in a pre-stressed left bay should release a measurable strain in the chord: the self-stress state made visible.
- **Does the result survive a change of tolerance, node numbering, member orientation?** Renumbering permutes rows and columns and changes no dimension; reversing a member's (i, k) negates its column and changes no subspace; the tolerance changes nothing for the four exact cases and everything for the near-mechanism sweep, which is the point of Experiment A.
- **Which decision changes if the geometric uncertainty doubled?** For a truss near rank deficiency, member forces scale as 1/(distance to the mechanism), so doubling the build tolerance doubles the force uncertainty; the design decision "add a member" is made when that uncertainty exceeds the margin, and the rank test is what names the member to add.

## Check yourself

1. State the dimensions of the four fundamental subspaces of a 5 × 8 matrix of rank 3.
2. Why are the special solutions of A x = 0 automatically independent?
3. The pivot columns of R = rref(A) are unit vectors. Why is the column space of A *not* spanned by the columns of R?
4. Maxwell's count for a planar truss is m + r − 2j. Write what it equals in terms of the null spaces of **A** and **Aᵀ**, and give a truss where it is zero but the truss moves.
5. What does Aᵀ **d** = 0 say, member by member and support by support?
6. The self-stress state of the Maxwell trap has forces (1, 1, −0.8, −0.8, −0.6, −0.6). Where do 0.8 and 0.6 come from?
7. `rank(A)` returns 12 for a truss whose σ_min is 10⁻⁷. Is the truss safe to build? What number would you quote instead?
8. Three reaction wheels have spin axes (1, 0, 0), (0, 1, 0) and (1, 1, 0)/√2. Which torque direction cannot be produced, and which subspace of which matrix contains it?
