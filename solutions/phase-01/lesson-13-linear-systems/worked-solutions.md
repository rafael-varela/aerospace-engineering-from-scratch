# P01.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15. Every number below is reproduced by `solution.py`; timings are from the author's machine and will differ on yours.

## 1. Elimination by hand

**Derivation.** s = 15, so the right-hand side is (8, −11, −3).

(a) No swaps. Column 1, pivot 2, multipliers −3/2 (row 2) and −1 (row 3):

```
[ 2   1  -1 |  8 ]      [ 2    1    -1  |  8 ]      [ 2   1   -1  |  8 ]
[-3  -1   2 | -11]  ->  [ 0   0.5   0.5 |  1 ]  ->  [ 0  0.5  0.5 |  1 ]
[-2   1   2 | -3 ]      [ 0    2     1  |  5 ]      [ 0   0   -1  |  1 ]
```

Column 2, pivot 0.5, multiplier 2/0.5 = 4. Pivots 2, 0.5, −1; product −1. Back substitution: z = −1, y = (1 − 0.5(−1))/0.5 = 3, x = (8 − 3 + (−1))/2 = 2.

(b) Partial pivoting. Column 1: |−3| is the largest, swap rows 1 and 2. Pivot −3; multipliers 2/(−3) = −2/3 and (−2)/(−3) = 2/3:

```
[-3  -1   2 | -11 ]      [-3   -1     2   | -11  ]      [-3   -1     2   | -11  ]
[ 2   1  -1 |  8  ]  ->  [ 0   1/3   1/3  |  2/3 ]  ->  [ 0   5/3   2/3  | 13/3 ]
[-2   1   2 | -3  ]      [ 0   5/3   2/3  | 13/3 ]      [ 0   1/3   1/3  |  2/3 ]
```

Column 2: 5/3 > 1/3, swap rows 2 and 3; multiplier (1/3)/(5/3) = 1/5; last row becomes [0 0 1/5 | −1/5]. Pivots −3, 5/3, 1/5, product −1. Row order (2, 3, 1): two swaps, sign +1, so det A = −1. Back substitution gives the same (2, 3, −1).

(c) From (a): L = [[1, 0, 0], [−3/2, 1, 0], [−1, 4, 1]], U = [[2, 1, −1], [0, 0.5, 0.5], [0, 0, −1]]. Multiply: row 2 of LU = −1.5 × (2, 1, −1) + (0, 0.5, 0.5) = (−3, −1, 2) ✓; row 3 = −1 × (2, 1, −1) + 4 × (0, 0.5, 0.5) + (0, 0, −1) = (−2, 1, 2) ✓. `lu_decompose` returns P = [[0, 1, 0], [0, 0, 1], [1, 0, 0]], L = [[1, 0, 0], [2/3, 1, 0], [−2/3, 1/5, 1]], U = [[−3, −1, 2], [0, 5/3, 2/3], [0, 0, 1/5]]: the pivoted factors of (b), with P A = L U to 10⁻¹⁶. `gaussian_elimination` returns (2, 3, −1) with residual 6 × 10⁻¹⁷; `determinant` returns −1.000.

**Result.** x = (2, 3, −1); pivots (2, 0.5, −1) and (−3, 5/3, 1/5); det = −1.

**Check.** Third original equation: −2(2) + 3 + 2(−1) = −3. ✓ Both pivot products are −1, as the determinant must be independent of the elimination route.

**Common mistakes.** Forgetting to apply the multiplier to b; swapping rows of U but not of the stored multipliers when pivoting late; reporting det = +1 after an odd number of swaps.

## 2. The cost and the factorisation

**Derivation.** (a) At step k (1-based) there are n − k rows below the pivot. Each needs one division for the multiplier and then, for each of the n − k entries to the right of the pivot plus the right-hand side, one multiply and one subtract: 2(n − k + 1) flops, plus the division, ≈ 2(n − k)(n − k + 1) per step. With j = n − k, Σ_{j=0}^{n−1} 2j(j + 1) = 2 Σ j² + 2 Σ j = (n−1)n(2n−1)/3 + (n−1)n = 2n(n−1)(n+1)/3 = 2n³/3 − 2n/3: the leading term is 2n³/3 (the divisions and the last column add O(n²)). Back substitution: row i needs n − i multiplies and subtracts and one division, Σ ≈ n². Gauss–Jordan on the n × 2n array [A | I]: each of n pivot steps updates all n − 1 other rows over about 2n columns, ≈ 2n³ flops: three solves' worth, plus the n² multiply A⁻¹ b you then need, which a solve never does.

(b) E = I − ℓ **e**_i **e**_kᵀ has −ℓ at position (i, k) and ones on the diagonal; E A replaces row i by row i − ℓ row k. (I − ℓ **e**_i **e**_kᵀ)(I + ℓ **e**_i **e**_kᵀ) = I + ℓ**e**_i**e**_kᵀ − ℓ**e**_i**e**_kᵀ − ℓ²**e**_i(**e**_kᵀ**e**_i)**e**_kᵀ = I since **e**_kᵀ **e**_i = 0 for i ≠ k. Elimination is E_N ⋯ E₁ A = U, so A = E₁⁻¹ ⋯ E_N⁻¹ U. For two inverses, (I + ℓ**e**_i**e**_kᵀ)(I + ℓ′**e**_{i′}**e**_{k′}ᵀ) = I + ℓ**e**_i**e**_kᵀ + ℓ′**e**_{i′}**e**_{k′}ᵀ + ℓℓ′**e**_i(**e**_kᵀ**e**_{i′})**e**_{k′}ᵀ, and the cross term survives only if k = i′, i.e. if the later step eliminates an entry in a row that served as the pivot row of the earlier step. In elimination order that never happens: a pivot row is finished once it has been used and is never modified again, so every later step has i′ > k′ ≥ k > … and i′ ≠ k. The product of all the inverses is therefore I plus each ℓ_{ik} in its own position (i, k), a unit lower triangular matrix L. Hence A = L U.

(c) det(E A) = det A for E = I − ℓ**e**_i**e**_kᵀ (det E = 1, triangular with unit diagonal); det(P_swap A) = −det A; det U = Π u_kk (triangular). So det A = (−1)^{swaps} Π pivots. Needed: multiplicativity det(XY) = det X det Y, the triangular rule, and the swap rule.

**Check.** For n = 1000, 2n³/3 = 6.7 × 10⁸ matches `elimination_flops(1000)` = 6.68 × 10⁸. ✓

**Common mistakes.** Counting multiply-subtract as one flop (gives n³/3; both conventions exist, say which); concluding "A⁻¹ is only three times a solve, so it is fine" (it is also less accurate per use).

## 3. A spacecraft thermal network

**Assumptions.** Linear conduction; radiation to space linearised as a conductance to a fixed reference; steady state.

**Derivation.** (a) Node 1: 2(T₂ − T₁) + 40 = 0; node 2: 2(T₁ − T₂) + 1.5(T₃ − T₂) = 0; node 3: 1.5(T₂ − T₃) − 0.5 T₃ = 0. As G **T** = **Q**:

G = [[2, −2, 0], [−2, 3.5, −1.5], [0, −1.5, 2.0]] W/K, **Q** = (40, 0, 0) W. Symmetric because conduction between i and j is the same conductance seen from either end. Each diagonal is the sum of the conductances touching that node (positive), each off-diagonal is minus the shared conductance; rows 1 and 2 sum to zero (no path to the reference), row 3 sums to +0.5 (the radiator's link to space).

(b) Elimination: multiplier −1 for row 2: row 2 ← row 2 + row 1 = [0, 1.5, −1.5 | 40]; multiplier −1 for row 3 relative to the new row 2: row 3 ← row 3 + row 2 = [0, 0, 0.5 | 40]. T₃ = 80 K; T₂ = (40 + 1.5 × 80)/1.5 = 106.67 K; T₁ = (40 + 2 × 106.67)/2 = 126.67 K. Code: (126.667, 106.667, 80.000) K. Radiator: 0.5 × 80 = 40 W. ✓

(c) G₂₃ = 3 W/K: (113.33, 93.33, 80.00) K. T₃ is unchanged: in steady state all 40 W must leave through the only path to space, so T₃ = Q₁/G₃ₛ regardless of the interior conductances. T₁ and T₂ drop by 13.3 K because the drop across G₂₃ halves (from 26.7 K to 13.3 K).

(d) G⁻¹ = [[3.167, 2.667, 2.0], [2.667, 2.667, 2.0], [2.0, 2.0, 2.0]] K/W. Column j is the temperature rise of every node per watt injected at node j: 1 W at the radiator raises every node by 2 K (it all goes to space through 0.5 W/K); 1 W at the box raises the box by 3.167 K = 2 + 1/1.5 + 1/2, the series resistance to space. All entries positive because heat injected anywhere can only raise temperatures: G is diagonally dominant with non-positive off-diagonals (an M-matrix), and such matrices have non-negative inverses.

**Result.** T = (126.7, 106.7, 80.0) K above reference; T₃ unchanged at 80 K when G₂₃ doubles; G⁻¹ as above.

**Check.** T₁ − T₂ = 20 K = 40 W / 2 W/K; T₂ − T₃ = 26.67 K = 40/1.5. ✓ Dimensions: W/(W/K) = K. ✓

**Common mistakes.** Writing +G off the diagonal; forgetting G₃ₛ on the radiator's diagonal (then G is singular: with no path to the reference, temperatures are defined only up to a constant).

## 4. Cost, pivoting, conditioning

**Result.** (a) Author's machine, best of three:

| n | ours [s] | numpy [s] | flops | ours [flop/s] |
|---|---|---|---|---|
| 50 | 3.5 × 10⁻³ | 3.8 × 10⁻⁵ | 8.6 × 10⁴ | 2.5 × 10⁷ |
| 100 | 1.3 × 10⁻² | 2.2 × 10⁻⁴ | 6.8 × 10⁵ | 5.1 × 10⁷ |
| 200 | 3.7 × 10⁻² | 4.3 × 10⁻⁴ | 5.4 × 10⁶ | 1.4 × 10⁸ |
| 400 | 1.5 × 10⁻¹ | 1.8 × 10⁻³ | 4.3 × 10⁷ | 2.9 × 10⁸ |
| 800 | 7.4 × 10⁻¹ | 9.0 × 10⁻³ | 3.4 × 10⁸ | 4.6 × 10⁸ |

Fitted exponents: ours 1.9 over the whole range, 2.2 over n ≥ 200 and still rising toward 3 (the Python-level loop over k costs a fixed amount per step that dominates at small n); NumPy 2.0 overall, 2.3 at the top end. Our speed climbs from 2.5 × 10⁷ to 4.6 × 10⁸ flop/s as the vectorised inner work grows; LAPACK is at 4 × 10¹⁰ flop/s at n = 800, a factor of 80.

(b)

| ε | pivot | no pivot | prediction 10⁻¹⁶/ε | refined ×1, ×2, ×3 |
|---|---|---|---|---|
| 10⁻⁸ | 0 | 5.0 × 10⁻⁹ | 10⁻⁸ | 0, 0, 0 |
| 10⁻¹⁰ | 0 | 8.3 × 10⁻⁸ | 10⁻⁶ | 0, 0, 0 |
| 10⁻¹² | 0 | 2.2 × 10⁻⁵ | 10⁻⁴ | 0, 0, 0 |
| 10⁻¹⁴ | 0 | 8.0 × 10⁻⁴ | 10⁻² | 0, 0, 0 |

The prediction is an upper bound (the actual rounding of 1 − y is a fraction of an ulp); the trend 1/ε is right. One refinement step already restores the exact answer: the unpivoted factorisation is a poor but *consistent* approximation of A, and the residual, computed with the true A, carries the correction.

(c)

| n | κ(H) | residual | error | after refinement |
|---|---|---|---|---|
| 4 | 1.6 × 10⁴ | 0 | 1.5 × 10⁻¹³ | 1.5 × 10⁻¹³ |
| 6 | 1.5 × 10⁷ | 7.8 × 10⁻¹⁷ | 3.2 × 10⁻¹⁰ | 3.1 × 10⁻¹⁰ |
| 8 | 1.5 × 10¹⁰ | 2.7 × 10⁻¹⁷ | 7.0 × 10⁻⁸ | 2.9 × 10⁻⁸ |
| 10 | 1.6 × 10¹³ | 1.5 × 10⁻¹⁶ | 1.6 × 10⁻⁴ | 5.9 × 10⁻⁴ |

Error ≈ κ × 10⁻¹⁶ to within a factor of 10, residual at rounding throughout. Refinement can repair a solve whose factorisation is inaccurate; it cannot supply digits that the problem itself does not determine, and for n = 10 it even wanders (the residual it works from is pure rounding noise).

**Check.** The n = 400 flop count 4.3 × 10⁷ at 2.9 × 10⁸ flop/s predicts 0.15 s. ✓

**Common mistakes.** Timing a single run (the first call includes import and cache warm-up); fitting the exponent through the overhead-dominated small-n points and reporting "elimination is O(n²)".

## 5. The truss with your load

**Assumptions.** Lesson truss; 15 kN at node 1 and s = 15 kN at node 4, both downward.

**Derivation.** (a) ΣF_x: R0x = 0. ΣM about node 0: R2y × 4 = 15 × 2 + 15 × 2 ⇒ R2y = 15 kN. ΣF_y: R0y = 30 − 15 = 15 kN.

(b) `solve_truss` (kN, tension positive): 0-1: +20.0; 1-2: 0; 3-4: 0; 4-5: −20.0; 0-3: 0; 1-4: 0; 2-5: −15.0; 0-4: −25.0; 1-5: +25.0; reactions (0, 15, 15). Largest compression: diagonal 0-4 at −25 kN (it carries the whole left reaction: 15/0.6 = 25). Zero-force members: 0-3 and 3-4 because node 3 is an unloaded joint with only two non-collinear members; 1-2 because at node 2 the vertical 2-5 must balance the 15 kN roller reaction and then nothing is left for the horizontal chord; 1-4 because at node 4 the two chord members 3-4 (zero) and 4-5 are collinear, the diagonal 0-4 and vertical 1-4 and the 15 kN load are the only other forces, and once 0-4's vertical component (+15 kN on node 4, since it pushes) cancels the load, 1-4 must be zero.

(c) Load (3, −15) kN at node 1: forces (13.0, 0, 0, −10.0, 0, +7.5, −7.5, −12.5, +12.5) kN, reactions (−3, 7.5, 7.5) kN. Only member 0-1 changes (from 10 to 13 kN): the horizontal load travels along the bottom chord to the pin, which takes it as R0x = −3 kN; nothing else sees it, because the vertical loads and the moment balance are unchanged.

(d) With member 3-5 in place of 1-5: the right bay 1-2-5-4 is a quadrilateral with no diagonal, a four-bar mechanism that can shear; the left bay has three verticals/diagonals (0-3, 0-4, 3-4, and now 3-5 spanning both bays) and is redundant. The count is still 12 = 12 but the columns are dependent: elimination reaches an exact zero pivot at step 8 and raises `ValueError: matrix is singular to working precision`. The message means the equilibrium equations do not determine the forces: there is a displacement mode with no member resisting it.

**Result.** Reactions (0, 15, 15) kN; largest compression 0-4 at −25 kN; four zero-force members; pin takes the 3 kN.

**Check.** Node 0 by hand with R0y = 15: F_{04} = −15/0.6 = −25, F_{01} = 0.8 × 25 = 20. ✓

**Common mistakes.** Calling a member "zero-force" from the numbers alone without the joint argument (the numbers are exactly 0 here, but 10⁻¹⁷ kN would also print as 0); forgetting that the pin, not the roller, carries horizontal load.

## 6. When to distrust a solve

(a) Not trustworthy as a design number. Forces 10⁴ times the load with a tiny residual is the signature of near-singularity: check `np.linalg.cond(A)` (expect ≥ 10⁴) or compare the determinant with the product of column norms. Cause: geometry close to a mechanism, e.g. nearly collinear members at a joint; the forces are real in the model but hypersensitive to millimetres.

(b) Not trustworthy. A relative residual of 10⁻³ from a direct solver means the solver did not do what you think: the matrix was modified after factorisation, or a mixed-precision or iterative solver stopped early. Sign: the residual itself; a direct solve should give ≈ 10⁻¹⁶ × κ at worst.

(c) Partly trustworthy: κ = 10¹⁴ with 16-digit arithmetic leaves about two significant digits. Sign: perturb the input by 10⁻¹⁰ relative and re-solve; if the forces move in the second digit, that is the real uncertainty. Cause: ill-conditioning (a near mechanism, or wildly different scales in the unknowns, which can be fixed by scaling).

(d) Not trustworthy. With all sources positive and a connected network, every temperature rise must be ≥ 0 (maximum principle for diffusion). A negative value means a sign error in G (an off-diagonal written +G, or a missing diagonal term); sign: a row of G that does not sum to the conductance to the reference.

**Common mistakes.** Treating the residual as the error; trusting a plausible-looking answer from a κ = 10¹⁴ system.

## 7. A tridiagonal system in O(n)

**Assumptions.** Linear springs, lumped masses, uniform gravity, static.

**Derivation.** (a) Number segments 1..n from the top. Segment i carries the weight of nodes i..n: T_i = (n − i + 1) m g; its extension is T_i/k. Tip displacement u_n = Σ_i T_i/k = m g n(n+1)/(2k) = 0.05 × 9.80665 × 1000 × 1001/(2 × 10⁵) = 2.4541 m over 1 km (0.25 % strain). Top tension T₁ = n m g = 490.33 N.

(b) Thomas algorithm: forward sweep c′₁ = c₁/b₁, d′₁ = d₁/b₁, then for i = 2..n, c′_i = c_i/(b_i − a_i c′_{i−1}), d′_i = (d_i − a_i d′_{i−1})/(b_i − a_i c′_{i−1}); back sweep x_n = d′_n, x_i = d′_i − c′_i x_{i+1}. About 8n = 8000 flops against `elimination_flops(1000)` = 6.7 × 10⁸: 80 000 times fewer.

(c) Thomas 0.8 ms (a pure Python loop), `np.linalg.solve` on the dense 1000 × 1000 matrix 10 ms, `gaussian_elimination` 1.4 s. All three give u_n = 2.4541141625 m, within 1.1 × 10⁻¹², 0.7 × 10⁻¹², 1.1 × 10⁻¹² of the exact value; top tension k u₁ = 490.3325 N. ✓

(d) The matrix is diagonally dominant: row 1 has 2k against one off-diagonal k; interior rows have 2k against k + k (equality, weak dominance). After step i of elimination the pivot is b′_i = k(i + 1)/i > k, strictly larger than the single off-diagonal entry −k below it, so `argmax` over the column always picks the current row and no swap ever happens. Elimination without pivoting is stable for diagonally dominant (and for symmetric positive definite) matrices, which is why banded solvers can skip it and keep the band.

**Result.** u_tip = 2.454 m, T_top = 490.3 N; 8n flops; speed-up 10³–10⁴ over dense elimination.

**Check.** The continuum limit: a uniform bar of length L, line density μ = m/1 m and stiffness EA = k × 1 m extends by μ g L²/(2EA) = 0.05 × 9.80665 × 10⁶/(2 × 10⁵) = 2.4517 m; the discrete sum has the extra factor (n+1)/n = 1.001. ✓

**Common mistakes.** Forgetting the free-end row has k, not 2k, on the diagonal (then the tip is held by a phantom spring); indexing c′_n (there is no super-diagonal in the last row).

## Check-yourself answers

1. A member from node 0 to node 4 with unit vector (0.8, 0.6): it pulls node 0 toward node 4 and node 4 toward node 0 when in tension.
2. The member exerts equal and opposite forces on its two end joints (Newton's third law), so the sum over the structure of its contributions is zero in each direction; an internal force cannot produce a net force on the whole truss.
3. The multiplier 10¹³ amplifies the rounding error of 1 − y, which is known only to about 10⁻¹⁶ absolute, into an error of 10⁻³ in x₁; no division by zero, just a loss of ten digits.
4. Eight times as many: 2(2000)³/3 = 5.3 × 10⁹ against 6.7 × 10⁸; about 5 s against 0.7 s at 10⁹ flop/s.
5. P is the permutation matrix recording the row swaps, L is unit lower triangular holding the multipliers ℓ_{ik} in position (i, k), U is the upper triangular result of elimination; P A = L U.
6. About 16 − 12 = 4 digits: the error is of order κ × 10⁻¹⁶ relative, whatever the residual says.
7. Members arranged so that one bay has no diagonal and another has two (same count, dependent columns); elimination meets an exact zero pivot and the code raises `ValueError` instead of returning numbers.
8. Factor A once into P, L, U (2n³/3 operations); each additional load case costs two triangular solves, about 2n² operations.
