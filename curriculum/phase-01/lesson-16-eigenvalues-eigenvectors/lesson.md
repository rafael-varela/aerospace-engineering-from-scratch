# P01.16  Eigenvalues & Eigenvectors

> An eigenvector is a direction the matrix does not turn. Find those directions and the matrix becomes three numbers.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the characteristic polynomial of a 3×3 matrix from its three invariants, and prove that a symmetric matrix has real eigenvalues and orthogonal eigenvectors. 2. Implement and test, from scratch, the closed-form symmetric 3×3 eigenvalues, power iteration, shifted inverse iteration and cyclic Jacobi rotations, against `numpy.linalg.eigh` as the oracle. 3. Find the principal axes and principal moments of a spacecraft inertia tensor, check them by a geometric calculation that uses no linear algebra, and decide about which axis the spacecraft can spin stably. 4. Measure the convergence of power iteration and Jacobi sweeps and predict the iteration count from the eigenvalue ratio. |
| **Prerequisites** | `P01.10` (moments of inertia as integrals), `P01.13` (determinants, solving linear systems), `P01.14` (null space: an eigenvector spans the null space of A − λI), `P01.15` (rotation matrices, C I Cᵀ), `P00.03`, `P00.07`. Skill check: you can expand a 3×3 determinant and solve a quadratic. |
| **Engineering question** | A 122 kg satellite carries its antenna and propellant tank off the symmetry planes: in which directions are its principal axes, how far are they tilted from the body axes the designers drew, and can it be spun about its z axis without wobbling? |
| **Assumptions and limits** | Rigid body; point-mass model (each component lumped at its centre of mass, so the tensor omits each component's own inertia about its own centre); torque-free motion for the spin-stability rule; the "with dissipation" rule assumes some internal energy loss, however small. The rule says nothing about *how fast* an unstable spin diverges, and nothing at all when two principal moments are equal. Power iteration fails when two eigenvalues have equal magnitude; Jacobi applies only to symmetric matrices. |
| **Analytical work** | Characteristic polynomial from the invariants; spectral theorem for symmetric matrices; the Rayleigh quotient as the moment of inertia about an axis; triangle inequalities for principal moments; convergence ratio of power iteration; the angle that zeroes one off-diagonal entry in a Jacobi rotation; hand solution of the 2×2 case. |
| **Computational deliverable** | `starter/main.py`: `center_of_mass`, `inertia_tensor`, `moment_about_axis`, `rotate_tensor`, `characteristic_coefficients`, `symmetric_eigenvalues_3x3`, `power_iteration`, `shifted_inverse_iteration`, `off_diagonal_norm`, `jacobi_eigen`, `rayleigh_quotient`, `principal_axes`, `triangle_inequalities_hold`, `spin_axis_stability`, `matrix_function_symmetric`, `linear_system_decay_rates`. Tests prove: the tensor obeys the parallel-axis theorem, the characteristic coefficients match `np.poly` and the invariants, every algorithm matches `numpy.linalg` to 10⁻¹⁰, eigenvectors are orthonormal and reconstruct the matrix, principal moments are rotation-invariant and satisfy the triangle inequalities, the moment about each principal axis computed from the geometry equals its eigenvalue, power iteration's cost scales as log(tol)/log(λ₂/λ₁), and Jacobi's off-diagonal norm shrinks faster each sweep. |
| **Practical activity** | **computational experiment**: measure the convergence of power iteration against the eigenvalue ratio and of Jacobi against sweep count, on the satellite tensor and on a 10×10 stiffness-pattern matrix. |
| **Equipment** | None |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `eigen_convergence.png`, the principal-axis table with the geometric check, and a short note on the spin-axis decision. |
| **Assessment** | `assessments/phase-01/grading-rubric.md`; competencies 1–4; pass when the spectral-theorem proof is complete, all public tests pass, the geometric check agrees with the eigenvalues to six figures, and the spin-axis note names the rule and its assumption. |
| **References** | `@mit-ocw-18-06-linear-algebra` Lecture 21 "Eigenvalues and eigenvectors", Lecture 22 "Diagonalization and powers of A", Lecture 25 "Symmetric matrices and positive definiteness" (required-free); Lecture 23 "Differential equations and exp(At)" (supplementary-free, for the stability preview); Lectures 18–20 on determinants (prerequisite). |
| **Extensions** | ★ Replace the point masses by solid blocks and add each block's own inertia; show the principal axes move by less than 1°. ★ Implement the QR algorithm (unshifted) for a symmetric matrix and compare its sweep count with Jacobi. ★ Integrate Euler's equations for the satellite spinning about each principal axis with a 1° perturbation and watch the intermediate axis tumble. |
| **Time** | unmeasured |

## 1. The engineering problem

A spacecraft's **inertia tensor** is the 3×3 matrix that relates its angular velocity to its angular momentum, **H** = I **ω**. Its diagonal entries are the familiar moments of inertia about the body axes; its off-diagonal entries, the **products of inertia**, are non-zero whenever mass sits off the symmetry planes, which on a real satellite with a boom antenna, a propellant tank and asymmetric payloads is always. With products of inertia present, spinning the satellite about its body z axis does *not* produce an angular momentum along z: **H** points somewhere else, and the satellite wobbles (the attitude-control term is nutation). Every attitude-control design starts by finding the three directions for which **H** is parallel to **ω**. Those are the **principal axes**, and the three moments about them are the **principal moments**.

Mathematically, "**H** parallel to **ω**" is I **ω** = λ **ω**: an **eigenvalue problem**. The eigenvectors of the inertia tensor are the principal axes, the eigenvalues are the principal moments, and the rotation matrix whose rows are the eigenvectors (Lesson 15) is the DCM from body to principal axes. The same mathematics decides whether a vibration mode is stable, whether a control loop diverges, what the natural frequencies of a wing are, and how fast a Markov process forgets its start; this lesson uses the inertia tensor because every number in it can be checked by a geometric calculation that involves no matrices at all.

The concrete problem: a 122 kg satellite modelled as eight point masses (bus corners, two solar arrays, an antenna boom tip, a propellant tank). Find its principal axes and moments, measure how far the major axis is tilted from body z, and decide whether a spin about body z is a sensible stabilisation choice. Success criteria: our from-scratch eigen-solvers agree with `numpy.linalg.eigh` to 10⁻¹⁰; the moment of inertia about each computed principal axis, evaluated directly as Σ m r⊥², equals the eigenvalue to six figures; the DCM to principal axes is a proper rotation; and the power-iteration count is predicted from the eigenvalue ratio to within 15 %.

## 2. Understand and derive

### 2.1 The inertia tensor of point masses

For point masses m_k at positions **r**_k = (x_k, y_k, z_k) from a reference point, the angular momentum about that point for angular velocity **ω** is **H** = Σ m_k **r**_k × (**ω** × **r**_k). The vector identity **a** × (**b** × **c**) = **b**(**a**·**c**) − **c**(**a**·**b**) turns this into **H** = Σ m_k [(**r**_k·**r**_k) **ω** − **r**_k (**r**_k·**ω**)] = I **ω**, with

$$
I = \sum_k m_k \left[(\mathbf{r}_k\cdot\mathbf{r}_k)\,I_3 - \mathbf{r}_k\mathbf{r}_k^{\mathsf T}\right]
= \sum_k m_k \begin{pmatrix} y_k^2 + z_k^2 & -x_k y_k & -x_k z_k \\ -x_k y_k & x_k^2 + z_k^2 & -y_k z_k \\ -x_k z_k & -y_k z_k & x_k^2 + y_k^2 \end{pmatrix}.
$$

Units: kg m². The matrix is **symmetric** by construction. Its diagonal entries are moments of inertia, each the sum of m times a squared distance from an axis, so they are positive. The off-diagonals are minus the products of inertia; the sign convention matters when you compare with a textbook that tabulates I_xy = Σ m x y without the minus sign.

The tensor depends on the reference point. For a shift **d** from the centre of mass, the **parallel-axis theorem** gives I_O = I_cm + M[(**d**·**d**) I₃ − **dd**ᵀ]; the test checks it and that the diagonal is smallest about the centre of mass. The reference point for attitude dynamics is the centre of mass, and that is the default in the code.

```
                 z (body)
                 ▲          ○ antenna tip (4 kg)      Bus corners: 4 x 22 kg at (±0.30, ±0.20, ±0.15) m
                 │        ╱                           Arrays:      2 x 7 kg at (0, ±1.60, 0.10) m
   array ○───────┼───────○ array                      Antenna tip: 4 kg at (0.90, 0.25, -0.55) m
         (y -)   │  ■ tank (16 kg)       (y +)         Tank:        16 kg at (-0.20, 0.10, 0.28) m
                 │
                 ●──────────▶ x (body, along the antenna)
```

### 2.2 The moment about any axis is a quadratic form

Project the angular momentum onto the spin axis **n** (a unit vector): the moment of inertia about that axis is I_n = **n**ᵀ I **n**. Check against the geometry: **n**ᵀ[(**r**·**r**) I₃ − **rr**ᵀ]**n** = |**r**|² − (**r**·**n**)² = |**r**⊥|², the squared perpendicular distance from the axis, so **n**ᵀ I **n** = Σ m_k |**r**_⊥,k|², the elementary definition. The expression **x**ᵀA**x**/**x**ᵀ**x** is called the **Rayleigh quotient**; for the inertia tensor it has this direct physical meaning, and `moment_about_axis` computes it from distances alone, giving us an oracle that owes nothing to linear algebra.

### 2.3 Eigenvalues: the characteristic polynomial

I **v** = λ **v** means (I − λ I₃) **v** = **0** has a non-zero solution, which by Lesson 14 means I − λ I₃ is singular: det(I − λ I₃) = 0 (`@mit-ocw-18-06-linear-algebra` Lecture 21). Expanding the determinant of a 3×3 matrix A gives a cubic whose coefficients are the three **invariants**:

$$
\det(\lambda I_3 - A) = \lambda^3 - I_1\lambda^2 + I_2\lambda - I_3, \qquad
I_1 = \operatorname{tr}A, \quad I_2 = \sum \text{principal } 2\times2 \text{ minors}, \quad I_3 = \det A.
$$

They are invariants because the eigenvalues do not change when the body is described in another frame (A → C A Cᵀ has the same characteristic polynomial, since det(λI − CACᵀ) = det(C(λI − A)Cᵀ) = det(λI − A)). Two checks fall out immediately and the tests use both: the sum of the eigenvalues is the trace, and their product is the determinant. For the satellite, I₁ = 111.973 kg m², I₂ = 3824.83 kg² m⁴, I₃ = 36 885.0 kg³ m⁶, and the three roots are 16.259, 43.200 and 52.515 kg m² (sum 111.973 ✓).

For a symmetric 3×3 matrix the cubic has a closed form. Shift by the mean eigenvalue q = tr A/3 so the quadratic term vanishes, set B = A − qI, p = tr(B²)/6, and the three roots are

$$
\lambda_k = q + 2\sqrt{p}\,\cos\!\left(\frac{\phi + 2\pi k}{3}\right), \quad k = 0, 1, 2, \qquad \cos 3\phi = \frac{\det B}{2p^{3/2}}.
$$

This is `symmetric_eigenvalues_3x3`, and it is the one place in the lesson where we get eigenvalues without iterating. It does not generalise: there is no closed form for n ≥ 5 (that is Abel's theorem), which is *why* every eigenvalue routine for larger matrices is iterative.

### 2.4 The spectral theorem: why symmetric matrices are easy

Two facts make the inertia tensor, and every stiffness or covariance matrix after it, tractable (`@mit-ocw-18-06-linear-algebra` Lecture 25).

**Real eigenvalues.** Let A **v** = λ **v** with A real symmetric and allow **v** complex for the moment. Take the conjugate transpose of both sides against **v**: **v***A**v** = λ **v***v**. The left side equals its own conjugate (because A = Aᵀ and A is real), so it is real; **v***v** = |**v**|² is real and positive; hence λ is real.

**Orthogonal eigenvectors.** If A **v**₁ = λ₁**v**₁ and A **v**₂ = λ₂**v**₂ with λ₁ ≠ λ₂, then **v**₂ᵀA**v**₁ = λ₁ **v**₂ᵀ**v**₁ and, transposing, **v**₂ᵀA**v**₁ = (A**v**₂)ᵀ**v**₁ = λ₂ **v**₂ᵀ**v**₁. Subtract: (λ₁ − λ₂) **v**₂ᵀ**v**₁ = 0, so **v**₁ ⊥ **v**₂. (With a repeated eigenvalue, any orthonormal basis of its eigenspace will do.)

So a symmetric A has an orthonormal eigenvector basis, A = V Λ Vᵀ with Vᵀ V = I₃. For the inertia tensor that means the three principal axes are mutually perpendicular: they form a frame, and if we order them so that det V = +1 the matrix C_p/b = Vᵀ is a rotation matrix in the sense of Lesson 15, with I_principal = C I Cᵀ = Λ diagonal. The rows of C are the principal axes in body coordinates.

The same decomposition gives functions of the matrix for free: A² = V Λ² Vᵀ, A⁻¹ = V Λ⁻¹ Vᵀ, √A = V √Λ Vᵀ, and exp(A) = V exp(Λ) Vᵀ (`@mit-ocw-18-06-linear-algebra` Lecture 22). The last one is the key to linear differential equations: d**u**/dt = A**u** has solution **u**(t) = exp(At) **u**(0), a sum of terms e^{λt} along the eigenvectors, so the system decays if every eigenvalue has a negative real part (Lecture 23). We test it on a damped oscillator; Phase 02 builds on it.

### 2.5 Two things the eigenvalues tell you about a body

**Triangle inequalities.** From the definition, I_xx + I_yy = Σ m (x² + y² + 2z²) ≥ Σ m (x² + y²) = I_zz, with equality only if every mass has z = 0, i.e. the body is flat. Since this holds in every frame, the principal moments satisfy I₁ + I₂ ≥ I₃ and its permutations. A "measured" tensor that violates them cannot belong to a real mass distribution: it is a sign error or a bad fit, and `triangle_inequalities_hold` is the check.

**Spin stability.** Euler's equations for torque-free motion, linearised about a pure spin at rate Ω about principal axis 3, give for the small transverse rates ω₁, ω₂:

$$
\ddot\omega_1 = -\frac{(I_3 - I_2)(I_3 - I_1)}{I_1 I_2}\,\Omega^2\,\omega_1 .
$$

The coefficient is positive, giving bounded oscillation, when I₃ is the largest or the smallest moment, and negative, giving exponential growth, when I₃ is the intermediate one. That is the **intermediate axis theorem**, and you can see it by flipping a book about each of its three axes. With internal energy dissipation (sloshing propellant, flexing arrays) only the *major* axis survives: at fixed **H**, kinetic energy H²/(2I) is smallest about the largest moment, and a dissipating body drifts to the lowest-energy state. The first US satellite, Explorer 1 (1958), was spun about its long, minor axis and began tumbling within hours because its whip antennas dissipated energy. `spin_axis_stability` encodes both rules.

### 2.6 Power iteration

Write a starting vector in the eigenvector basis, **x**₀ = Σ c_k **v**_k. Then A^N **x**₀ = Σ c_k λ_k^N **v**_k = λ₁^N [c₁**v**₁ + Σ_{k>1} c_k (λ_k/λ₁)^N **v**_k]. If |λ₁| > |λ_k| for all other k, every term but the first dies as (λ_k/λ₁)^N, and the normalised iterate converges to **v**₁. The error in the vector falls by the factor r = |λ₂/λ₁| per step; the Rayleigh quotient, being quadratic in the error for a symmetric matrix, converges as r^{2N}. The number of steps to reach a tolerance ε is about log ε / log r. For the satellite r = 43.200/52.515 = 0.8226 and ε = 10⁻¹² gives 142; the run needs 130 (the start vector happens to have a smaller component along **v**₂). When two eigenvalues have the same magnitude, r = 1 and the iteration never settles: `test_power_iteration_handles_negative_dominant_eigenvalue_and_raises_when_stuck` forces exactly that with diag(3, −3, 1).

**Shifted inverse iteration** applies the same idea to (A − σI)⁻¹, whose dominant eigenvalue is 1/(λ_j − σ) for the λ_j nearest the shift σ. Each step is a linear solve (Lesson 13), never an explicit inverse. With σ = 0 it finds the smallest moment: ratio 16.259/43.200 = 0.376, predicted 28 steps, measured 29.

### 2.7 Jacobi rotations

Power iteration gives one eigenpair at a time. For a symmetric matrix, Jacobi's method finds all of them by applying plane rotations that zero one off-diagonal pair at a time. With J the identity except for a 2×2 rotation block in rows/columns p and q (c = cos θ, s = sin θ at (p,p), (q,q), (p,q), (q,p) as +s, −s), the update A ← JᵀAJ keeps A symmetric, keeps its eigenvalues (it is a similarity by an orthogonal matrix), and changes only rows and columns p and q. Writing out the new (p, q) entry and setting it to zero gives, with t = tan θ and τ = (a_qq − a_pp)/(2a_pq),

$$
t^2 + 2\tau t - 1 = 0 \quad\Longrightarrow\quad t = \frac{\operatorname{sign}\tau}{|\tau| + \sqrt{1 + \tau^2}}, \qquad c = \frac{1}{\sqrt{1 + t^2}}, \quad s = tc.
$$

The smaller root is chosen so that |θ| ≤ 45°, which keeps the method stable. Zeroing one entry disturbs the others, but the sum of squares of all off-diagonal entries falls by exactly 2a_pq² at each rotation (the Frobenius norm is preserved by the similarity and the diagonal gains what the off-diagonal loses), so the **off-diagonal norm** decreases monotonically and, once small, roughly squares each sweep. Accumulating V ← VJ gives the eigenvectors. On the 10×10 matrix in `main()` the off-diagonal norm goes 1.669, 0.392, 0.0406, 1.7 × 10⁻⁴, 1.4 × 10⁻⁹, 7 × 10⁻²⁰: reduction factors 0.23, 0.10, 0.004, 8 × 10⁻⁶, …, which is what quadratic convergence looks like.

### 2.8 Where the model stops being valid

The point-mass model omits each component's inertia about its own centre, which for a compact bus is a few per cent of the total (extension 1 quantifies it). The spin-stability rule is a *linear* rule: it tells you the intermediate axis is unstable, not how long the tumble takes to develop, and it says nothing when two moments are equal (a symmetric top), where the full non-linear equations are needed. Power iteration assumes a dominant eigenvalue; Jacobi assumes symmetry. Neither handles the non-symmetric state matrices of flight dynamics, whose eigenvalues are complex: that is the QR algorithm's job, and it is why `linear_system_decay_rates` calls `numpy.linalg.eigvals` instead of reinventing it.

## 3. Build the computational model

Open `starter/main.py`. `satellite_point_masses()` is given; everything else is yours.

### Step 1: the tensor and its geometric oracle

```python
def inertia_tensor(masses, positions, about=None):
    origin = center_of_mass(masses, positions) if about is None else np.asarray(about, float)
    d = positions - origin
    I = np.zeros((3, 3))
    for mk, rk in zip(masses, d):
        I += mk * (np.dot(rk, rk) * np.eye(3) - np.outer(rk, rk))
    return I
```

Validate the inputs (positive masses, matching shapes) and raise `ValueError` otherwise. `moment_about_axis` must *not* use the tensor: compute Σ m (|**r**|² − (**r**·**n**)²) directly, so that it can serve as an independent check. `rotate_tensor(I, C)` is `C @ I @ C.T`.

### Step 2: invariants and the closed form

`characteristic_coefficients` returns `[1, -I1, I2, -I3]` with the determinant by cofactor expansion (write it out; it is six products). `symmetric_eigenvalues_3x3` implements 2.3; clip the cosine argument to [−1, 1] before `arccos`, because rounding can push it to 1.0000000000000002, and return the roots sorted. The test compares both against `np.poly` and `np.linalg.eigvalsh`.

### Step 3: power iteration and its shifted inverse

```python
def power_iteration(A, x0=None, tol=1e-12, max_iter=10000):
    x = np.ones(A.shape[0]) if x0 is None else np.array(x0, float)
    x /= np.linalg.norm(x)
    for k in range(1, max_iter + 1):
        y = A @ x
        y /= np.linalg.norm(y)
        if y @ x < 0: y = -y          # a negative eigenvalue flips the sign every step
        delta = np.linalg.norm(y - x); x = y
        if delta < tol: return x @ A @ x, x, k
    raise RuntimeError("no convergence")
```

Return the iteration count: the test uses it to confirm the log ε / log r law. `shifted_inverse_iteration` replaces `A @ x` by `np.linalg.solve(A - shift * I, x)` and must raise `RuntimeError` if the shift is exactly an eigenvalue (NumPy raises `LinAlgError` for the singular solve; catch it).

### Step 4: Jacobi

Loop over sweeps; within a sweep loop over all pairs p < q; skip a pair whose entry is already zero; compute τ, t, c, s as in 2.7; build J; update `A = J.T @ A @ J` and `V = V @ J`; set the (p, q) and (q, p) entries to exactly zero (they are zero by construction; this removes rounding dust). After each sweep compute `off_diagonal_norm(A)`, append it to `history` if one was passed, and stop when it is below `tol` times the Frobenius norm of the input. Return the diagonal sorted ascending with the columns of V permuted to match. Two subtleties the tests catch: a diagonal input must return after one sweep with its eigenvalues *sorted*, and a non-symmetric input must raise `ValueError`.

### Step 5: principal axes and the two rules

`principal_axes(I)` calls `jacobi_eigen`, flips the sign of the third eigenvector if det V < 0 so the triad is right-handed, and returns `(moments, V.T)`. `triangle_inequalities_hold` sorts the moments and checks a + b ≥ c (with a relative tolerance so that a flat plate, where equality holds, passes). `spin_axis_stability(moments, k)` returns `{"rigid": is_max or is_min, "with_dissipation": is_max}` and raises `ValueError` when the axis is not isolated.

### Step 6: functions of a matrix and the stability preview

`matrix_function_symmetric(A, f)` is `V @ np.diag(f(lam)) @ V.T` with the Jacobi output. `linear_system_decay_rates(A)` returns the sorted real parts of `np.linalg.eigvals(A)`.

Running `main()` prints, for the satellite:

```
inertia tensor about the CM [kg m^2]:
 [[44.188  -0.5715  2.8881]
  [-0.5715 16.2721  0.1804]
  [ 2.8881  0.1804 51.5133]]
principal moments: [16.2588 43.1996 52.515 ]
principal axis 0: moment 16.259, nearest body axis y, tilted  1.27 deg
principal axis 1: moment 43.200, nearest body axis x, tilted 19.17 deg
principal axis 2: moment 52.515, nearest body axis z, tilted 19.13 deg
```

followed by the geometric check (identical to six figures), the triangle test, the three spin verdicts, the iteration counts and the Jacobi history.

### Verify

- **Analytic case:** the characteristic coefficients equal the trace, the sum of principal minors and the determinant, and every eigenvalue from `numpy` is a root; the closed-form roots match `eigvalsh` to 10⁻¹⁰; a diagonal matrix returns its diagonal; a multiple of the identity returns a triple root.
- **Limiting cases:** a negative dominant eigenvalue (`diag(−5, 2, 1)`) is found with the right sign; equal-magnitude eigenvalues (`diag(3, −3, 1)`) make power iteration raise; a shift equal to an eigenvalue makes inverse iteration raise; the flat plate (1, 1, 2) sits exactly on the triangle boundary.
- **Conservation:** the trace, determinant and Frobenius norm of A equal the sum, product and root-sum-square of the Jacobi eigenvalues; V Λ Vᵀ reconstructs A; the principal moments are unchanged when the tensor is expressed in a rotated frame.
- **Convergence:** the extra iterations between tol = 10⁻⁶ and 10⁻¹² equal log(10⁻⁶)/log(λ₂/λ₁) within 15 %; the Jacobi off-diagonal norm decreases monotonically with a reduction factor that shrinks by more than ten between the first and last sweep.

### Compare

`numpy.linalg.eigh` (symmetric) and `eigvals` (general) are the oracles. They use LAPACK: a Householder reduction to tridiagonal form followed by the implicit QR algorithm, which costs about 4n³/3 flops and is backward stable. Jacobi costs 3 to 5 times more for large n but is simpler to write, embarrassingly parallel, and slightly more accurate for tiny eigenvalues; power iteration costs n² per step and is the right tool when you want one eigenpair of a huge sparse matrix (it is the core of PageRank and of the Lanczos method). We wrote our own to know what "the eigenvector" means operationally and to be able to recognise when the library's answer cannot be trusted: for instance, eigenvectors of nearly equal eigenvalues are ill-defined individually, and the library returns *some* orthonormal pair without warning.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Map the convergence of the two iterative solvers.

1. Power iteration on the satellite tensor from **x**₀ = (1, 1, 1)/√3: record the relative error of the Rayleigh quotient against λ_max after each of 60 iterations. Before running, predict the slope on a semilog plot from r = λ₂/λ₁.
2. Jacobi on the satellite tensor and on the 10×10 matrix K = tridiag(−1, 2, −1) + diag(0, 1/9, …, 1) with `tol = 1e-13`, recording the off-diagonal norm after each sweep. Predict how many sweeps the 10×10 needs.
3. Change the satellite by moving the antenna tip to (0.90, 0.25, −0.85) m and repeat step 1. Predict whether the iteration count rises or falls from the new eigenvalue ratio *before* running.
4. Plot both panels and save as `eigen_convergence.png`.

Raw data format: CSV with columns `k, rel_error_power` and `sweep, off_norm`.

## 5. Validate and improve

The run gives, for power iteration on the satellite tensor:

| iteration k | relative error in λ_max | r^{2k}, r = 0.8226 |
|---|---|---|
| 1 | 5.65 × 10⁻² | 6.77 × 10⁻¹ |
| 5 | 5.34 × 10⁻³ | 1.42 × 10⁻¹ |
| 10 | 7.77 × 10⁻⁴ | 2.01 × 10⁻² |
| 20 | 1.57 × 10⁻⁵ | 4.05 × 10⁻⁴ |
| 40 | 6.37 × 10⁻⁹ | 1.64 × 10⁻⁷ |
| 60 | 2.58 × 10⁻¹² | 6.66 × 10⁻¹¹ |

The *slope* is right: from k = 20 to k = 40 the error falls by a factor 2470 and r⁴⁰ = 2470. The *level* is 25 times below the bound because the bound assumes the starting vector has equal components along every eigenvector; (1, 1, 1)/√3 happens to have a small component along **v**₂. The prediction of the iteration count (142 against 130 measured) inherits the same offset: the law gives the slope exactly and the intercept only to within the starting vector's geometry. State this in your note; a convergence "law" that predicts the slope and not the level is still a law, and it is the slope that tells you what a tighter tolerance will cost.

For Jacobi, the 3×3 tensor is diagonal to rounding after two sweeps (off-diagonal norm 0.111 then exactly 0); the 10×10 needs six, with the factors listed in 2.7. The tolerance 10⁻¹³ is not arbitrary: each rotation leaves rounding of order ε|A| in the touched rows, so the off-diagonal norm cannot fall below about n ε |A|, and a tolerance of 10⁻¹⁵ never terminates for n = 10. That failure mode is worth seeing once; make it happen and then set the tolerance from the floor rather than from hope.

Improvement, feasible now: start power iteration from a random vector instead of (1, 1, 1)/√3 and show that the level moves while the slope does not; then implement the Rayleigh-quotient shift in inverse iteration (update σ to the current Rayleigh quotient each step) and show the convergence become cubic, typically four steps to 10⁻¹² from any reasonable start.

## 6. Defend your work

Submit `main.py`, the test log, `eigen_convergence.png` with its CSV, the table of principal moments with the geometric check, and a one-page note (template in `projects/_template/report-template.md`) that states the three principal moments, the tilt of each principal axis from its nearest body axis, the spin-axis recommendation with the rule it rests on, and the convergence-law comparison.

### Adversarial self-checks

- What if every mass doubles? Every moment doubles, every principal axis stays put, every eigenvalue ratio is unchanged, so the iteration counts are unchanged. What if every position doubles instead?
- What if the antenna tip moves onto the x axis and the tank onto the z axis? The products of inertia I_xy and I_yz vanish but I_xz does not: one principal axis is exactly body y and the other two are tilted in the x–z plane. Check with the code.
- Which assumption fails first in the real world? The rigid body: propellant slosh in a half-full tank changes the tensor with time and dissipates energy, which is precisely the case where only major-axis spin is safe.
- What measurement could falsify the model? A spin-balance test: spin the satellite about body z on a bearing and measure the transverse angular momentum (a wobble torque on the bearing). Its magnitude predicts the product-of-inertia column I_xz, I_yz directly.
- Does the result survive a change of solver? Jacobi, the closed form and `eigh` agree to 10⁻¹⁰; power iteration agrees for the largest moment. Which decision changes if the two largest moments were within 1 % of each other? The major axis becomes ill-defined (its eigenvector swings with tiny changes in the tensor), and a spin-stabilised design would be replaced by three-axis control.

## Check yourself

1. Define eigenvalue and eigenvector in one sentence each, and say what they mean physically for an inertia tensor.
2. What are the three invariants of a 3×3 matrix, and why are they called invariants?
3. Prove in three lines that a real symmetric matrix has real eigenvalues.
4. Why is the Rayleigh quotient **n**ᵀI**n** the moment of inertia about **n**? What are its minimum and maximum over all unit **n**?
5. Give the convergence ratio of power iteration and the number of steps to reach 10⁻⁹ when λ₂/λ₁ = 0.9.
6. State the intermediate-axis theorem and the additional rule that applies when energy is dissipated.
7. The satellite's principal moments are 16.26, 43.20 and 52.52 kg m². Verify the triangle inequalities and compute the trace and determinant of the body-axis tensor without the tensor.
8. Why does `linear_system_decay_rates` call `numpy.linalg.eigvals` instead of `jacobi_eigen`?
