# P02.05  Matrix Exponential & Linear Systems

> Every linear system with constant coefficients has one solution operator, e^{At}. Build it three ways from scratch, prove what it must satisfy, and use it to fly a spacecraft to a docking port with two burns computed by a 2×2 linear solve.

| | |
|---|---|
| **Learning outcomes** | 1. Define the matrix exponential by its series, derive the properties of the state transition matrix Φ(t) = e^{At} (Φ(0) = I, Φ(t)Φ(s) = Φ(t + s), Φ⁻¹(t) = Φ(−t), dΦ/dt = AΦ, det Φ = e^{tr A·t}) and the condition under which e^{A+B} = e^A e^B. 2. Derive the forced solution x(t) = Φ(t)x₀ + ∫Φ(t − s)Bu(s) ds, the augmented-matrix form of the zero-order-hold discretisation and the stability criterion Re λ < 0. 3. Implement and test the exponential by Taylor series, by scaling and squaring and by eigen-decomposition, with `scipy.linalg.expm` as oracle and the Clohessy–Wiltshire closed-form transition matrix as an independent one, and explain when each method fails. 4. Use Φ to solve a two-impulse rendezvous as a linear system and judge the model's validity. |
| **Prerequisites** | `P02.04` (eigen-solution of planar systems, classification), `P02.01` (RK4, `solve_ivp`), `P01.16` (eigen-decomposition A = VΛV⁻¹, defective matrices), `P01.13` (matrix products, inverses, block matrices), `P00.07` (Taylor series). Skill check: you can diagonalise a 2×2 matrix and write e^x as a series. |
| **Engineering question** | A chaser spacecraft is 50 m below and 1 km behind a station in a 400 km circular orbit, at rest relative to it. Where does it drift in half an orbit if nothing is done? What two impulsive burns bring it to the station, at rest, in a quarter, a half and three quarters of an orbit, and how much Δv does each plan cost? And is a 10-minute, 0.2 mm/s² thrust arc computed by holding the thrust constant over 10 s steps exactly right, or only approximately? |
| **Assumptions and limits** | Clohessy–Wiltshire (Hill) relative motion: a circular target orbit, a chaser within a few kilometres (distances ≪ orbit radius, so gravity is linearised), two-body point-mass gravity with no J₂, drag or third bodies, impulsive burns (instantaneous velocity changes) and, for the thrust arc, a constant acceleration held exactly over each step. In-plane motion only. The matrix is defective (eigenvalues 0, 0, ±in), which is a feature of the physics, not a numerical accident. |
| **Analytical work** | The exponential series and its convergence; the group, inverse, derivative and determinant properties; the commutation condition; variation of parameters for the forced system; the augmented-matrix identity for the zero-order hold; the closed-form Clohessy–Wiltshire transition matrix by direct integration; the two-impulse rendezvous as Φ_rv Δv₁ = −(Φ_rr r₀ + Φ_rv v₀). |
| **Computational deliverable** | `starter/main.py`: `expm_taylor`, `expm_scaling_squaring`, `expm_eig`, `state_transition`, `propagate`, `is_stable`, `forced_response_zoh`, `discretize`, `modal_decomposition`, `cw_stm_closed_form`, `rendezvous_delta_v` (`iss_orbit`, `cw_matrix`, `cw_rhs` are given). Tests prove: the Taylor series is exact for nilpotent matrices, converges for small norm and fails for large; scaling and squaring matches `scipy.linalg.expm` to 10⁻⁹ on random 2×2 to 6×6 matrices; the eigen method matches for diagonalisable matrices and raises for defective ones; e^{A+B} = e^A e^B exactly when A and B commute; Φ satisfies all five properties and dΦ/dt = AΦ to 10⁻⁷; `propagate` agrees with `solve_ivp` and with `P02.04`'s closed form; the zero-order-hold response equals a tight `solve_ivp` integration of the same piecewise-constant input to 10⁻⁹; the Clohessy–Wiltshire closed form matches our exponential and SciPy's to 10⁻⁸ at five times including one and two half-orbits, has unit determinant, and shows the −12π x₀ per-orbit drift; the two-impulse rendezvous arrives at the origin at rest for four transfer times and raises for a whole orbit. |
| **Practical activity** | **computational experiment**: map the Δv cost of the two-impulse rendezvous against transfer time from 0.1 to 1.9 orbits and locate the singularity at one orbit; measure where the plain Taylor series breaks down as ‖At‖ grows and how many squarings scaling-and-squaring uses; compare zero-order-hold steps of 1, 10 and 60 s against `solve_ivp` for the thrust arc. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `state_transition.png`, the hand derivations of the five Φ properties and of the Clohessy–Wiltshire first column, the Δv-versus-transfer-time table, and a half-page note recommending a transfer time with its cost and risk. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the group and derivative properties and the augmented-matrix identity are derived; all public tests pass; the Taylor breakdown and the singular transfer time are explained with numbers; the note states a Δv and a reason. |
| **References** | `@mit-ocw-18-03-differential-equations` Lecture 29 "Matrix Exponentials", Lecture 28 "Matrix Methods for Inhomogeneous Systems", Lecture 30 "Decoupling Linear Systems with Constant Coefficients" (required-free); notes 36 "Normal modes and the matrix exponential" (supplementary-free). `@mit-ocw-18-06-linear-algebra` Lecture 23 "Differential equations and exp(At)" (required-free): e^{At} = S e^{Λt} S⁻¹ and stability iff Re λ < 0. |
| **Extensions** | ★ Implement the Padé approximant of order (6, 6) used inside library `expm` routines and compare its term count with the Taylor series at the same accuracy. ★ Add the out-of-plane Clohessy–Wiltshire equation (a decoupled harmonic oscillator) and plan a three-dimensional rendezvous. ★ Replace the impulsive burns by finite thrust arcs via `forced_response_zoh` and quantify the gravity-loss-like penalty. |
| **Time** | unmeasured |

## 1. The engineering problem

`P02.04` solved d**x**/dt = A**x** by finding eigenvalues and eigenvectors and assembling **x**(t) = V diag(e^{λt}) V⁻¹ **x**₀. The matrix in the middle depends only on A and t, not on **x**₀. It is the **state transition matrix** Φ(t): the single linear map that carries every initial state forward by time t. For a scalar equation ẋ = ax it is e^{at}; for a matrix it is, by definition, **e^{At}**, and the whole of linear systems theory, from flight control to orbit determination, is built on it. This lesson constructs it from scratch three different ways, proves what it must satisfy, and tests each construction against the others.

The vehicle for the lesson is spacecraft rendezvous. Relative to a target in a circular orbit the chaser's motion, for separations small compared with the orbit radius, is governed by the **Clohessy–Wiltshire equations**, a 4×4 linear system whose matrix has a *repeated zero eigenvalue with only one eigenvector*. That defect is physical (a chaser slightly below the target orbits faster and drifts ahead without bound), and it breaks the eigen-decomposition method of `P02.04`. The series methods do not care, and the closed-form transition matrix, which can be derived by hand, is an independent oracle against which `scipy.linalg.expm` and our own code are both checked. With Φ in hand, the two-impulse rendezvous is a linear solve: pick the arrival time, partition Φ(T) into position and velocity blocks, solve a 2×2 system for the first burn.

Success criteria, written first: our scaling-and-squaring exponential matches `scipy.linalg.expm` to 10⁻⁹ relative on random matrices up to 6×6 and the Clohessy–Wiltshire closed form to 10⁻⁸ over a full orbit (where ‖At‖ is in the thousands); all five properties of Φ hold to 10⁻¹⁰; the zero-order-hold response of a thrust arc matches `solve_ivp` to 10⁻⁹ m; and the planned rendezvous arrives at the origin to better than 10⁻⁸ m in the model.

## 2. Understand and derive

### 2.1 Definition and the first properties

Define, for any square matrix A,

**e^{A} = I + A + A²/2! + A³/3! + ⋯ = Σ_{k≥0} A^k/k!**.   (1)

The series converges for every A because ‖A^k/k!‖ ≤ ‖A‖^k/k! and the scalar series Σ‖A‖^k/k! = e^{‖A‖} converges (18.03 Lecture 29). Replace A by At to get Φ(t) = e^{At}. Differentiating (1) term by term (legitimate for a power series), dΦ/dt = A + A²t + A³t²/2! + ⋯ = A e^{At}, so **dΦ/dt = AΦ, Φ(0) = I**, and **x**(t) = Φ(t)**x**₀ satisfies ẋ = A**x** with the right initial condition: Φ is the solution operator, and by uniqueness (`P02.01`) the only one.

Three more properties follow:

- **Group property** Φ(t)Φ(s) = Φ(t + s). Both sides, as functions of t, satisfy dX/dt = AX with X(0) = Φ(s); uniqueness makes them equal. Consequence: Φ(t)Φ(−t) = Φ(0) = I, so **Φ⁻¹(t) = Φ(−t)** and e^{At} is always invertible, which is why trajectories of a linear system never cross or merge.
- **Determinant** (Liouville): d(det Φ)/dt = tr(A) det Φ, so **det Φ(t) = e^{tr(A) t}**. For the Clohessy–Wiltshire matrix tr A = 0 and Φ preserves area in state space; for a damped system tr A < 0 and volumes shrink.
- **Commutation.** e^{A+B} = e^A e^B **if and only if AB = BA** (for the "if": expand both sides and compare, which is the binomial theorem and needs commutation; the "only if" direction is shown by counterexample in Problem 2). With A = [[0, 1], [0, 0]] and B = [[0, 0], [1, 0]]: e^A e^B = [[2, 1], [1, 1]] but e^{A+B} = [[cosh 1, sinh 1], [sinh 1, cosh 1]] = [[1.543, 1.175], [1.175, 1.543]]. This is why e^{At} cannot be built one matrix element at a time.

### 2.2 Three ways to compute it

**Taylor series** (`expm_taylor`). Sum (1) to N terms with the recurrence T_k = T_{k−1}A/k. Exact for nilpotent matrices (A^m = 0 for some m: the series terminates), excellent when ‖A‖ ≲ 1, and hopeless when ‖A‖ is large: for the short-period matrix of `P02.04` times 30 s, ‖A‖₁ = 144 and the 30-term sum is wrong by 10²² while e^{At} itself has entries of order e^{−30}. The terms grow like 144^k/k! until k ≈ 144 before they shrink, and the floating-point cancellation among terms of size 10⁶⁰ destroys everything.

**Scaling and squaring** (`expm_scaling_squaring`). Use e^{A} = (e^{A/2^s})^{2^s}: choose s so that ‖A‖/2^s ≤ ½, sum the series for the small matrix (20 terms give 10⁻¹⁶), then square the result s times. For ‖A‖ = 144, s = 9 squarings. This is the idea inside every library `expm` (SciPy replaces the Taylor series by a Padé rational approximant, extension 1, but the scaling and squaring are the same). Its weakness is that squaring amplifies rounding error by about 2^s · cond, which is why the test asks for 10⁻⁹ rather than 10⁻¹⁵.

**Eigen-decomposition** (`expm_eig`, 18.06 Lecture 23). If A = VΛV⁻¹ then A^k = VΛ^kV⁻¹ and (1) gives **e^{A} = V e^{Λ} V⁻¹** with e^{Λ} = diag(e^{λ_k}): the normal-mode form, the matrix statement of `P02.04`'s solution, and the one that *explains* the behaviour (each mode grows or decays as e^{λt}). It fails when V is singular or nearly so, that is when A is **defective**: a repeated eigenvalue with too few eigenvectors. Then e^{At} contains terms t e^{λt} (`P02.04` eq. 2) that no diagonal e^{Λ} can produce. `expm_eig` checks cond(V) and raises; the Clohessy–Wiltshire matrix is the lesson's live example (cond(V) = ∞).

### 2.3 Stability

From the eigen form, **x**(t) = Σ_k c_k e^{λ_k t}**v**_k. Every solution decays to zero if and only if **Re λ_k < 0 for all k** (asymptotic stability). If some Re λ = 0 and the rest are negative, solutions stay bounded (marginal stability) unless the zero-real-part eigenvalue is defective, in which case t e^{λt} grows: the Clohessy–Wiltshire drift is exactly this. If any Re λ > 0, some solutions grow. `is_stable` tests the strict condition; for a 2×2 it reduces to τ < 0 and Δ > 0 (`P02.04`), and for a damped oscillator to c > 0.

### 2.4 Forced systems: variation of parameters and the zero-order hold

For ẋ = A**x** + B**u**(t), try **x** = Φ(t)**c**(t) (Lecture 28). Then Φ**ċ** = B**u**, so **ċ** = Φ(−t)B**u** and

**x**(t) = Φ(t)**x**₀ + ∫₀ᵗ Φ(t − s) B **u**(s) ds.   (2)

The integral is a superposition of impulse responses: the input at time s is carried forward by Φ(t − s). When **u** is **held constant** over an interval of length h (a zero-order hold, which is what a digital controller or a thruster firing at constant acceleration does), (2) over one interval gives

**x**_{k+1} = A_d **x**_k + B_d **u**_k,  A_d = e^{Ah},  B_d = ∫₀ʰ e^{As} ds · B.   (3)

The neat way to get both at once is the **augmented matrix**: with M = [[A, B], [0, 0]] (size n + m),

e^{Mh} = [[A_d, B_d], [0, I]],   (4)

because the series for e^{Mh} has upper-right block Σ_{k≥1} A^{k−1}B h^k/k! = (∫₀ʰ e^{As} ds) B. One exponential, no quadrature, and (3) is *exact* for the piecewise-constant input, not an approximation; the test checks this against `solve_ivp` at 10⁻⁹. For small h, A_d ≈ I + Ah and B_d ≈ Bh (the explicit Euler step is the first-order truncation of the exact map); for the Clohessy–Wiltshire thrust arc at h = 10 s the position rows of B_d are ½h² = 50 m per m/s², which Euler's Bh misses entirely.

### 2.5 The Clohessy–Wiltshire equations

Place the origin at the target, x radially outward, y along the velocity. In this frame rotating at the mean motion n = √(μ/a³), the chaser's equations, linearised in the small separation (`P02.06` will show the linearisation; here we take the result), are

**ẍ = 3n²x + 2nẏ,  ÿ = −2nẋ**.   (5)

The 2nẏ and −2nẋ terms are Coriolis; 3n²x is the tidal term (gravity gradient plus centrifugal). As a first-order system in (x, y, ẋ, ẏ):

A_CW = [[0, 0, 1, 0], [0, 0, 0, 1], [3n², 0, 0, 2n], [0, 0, −2n, 0]],  eigenvalues 0, 0, ±in, tr A = 0.

The y equation integrates once: ẏ = −2nx + (ẏ₀ + 2nx₀). Substituting into the x equation gives ẍ + n²x = 2n(ẏ₀ + 2nx₀): a harmonic oscillator at the orbital frequency driven by a constant, whose solution is x(t) = (4 − 3 cos nt)x₀ + (sin nt/n)ẋ₀ + (2(1 − cos nt)/n)ẏ₀. Integrating ẏ then gives y(t) = 6(sin nt − nt)x₀ + y₀ − (2(1 − cos nt)/n)ẋ₀ + ((4 sin nt − 3nt)/n)ẏ₀, and differentiating both gives the velocity rows. These four rows are the **closed-form transition matrix** `cw_stm_closed_form`, derived with nothing but one-variable calculus, and they contain the physics: the **secular terms** −6nt x₀ and −3nt ẏ₀ in y are the t e^{0·t} of the defective zero eigenvalue. A chaser 50 m below the target (x₀ = −50) is on a faster orbit and drifts *ahead* by −6n·(2π/n)·x₀ = −12π x₀ = **1885 m per orbit**; the test checks this coefficient and that a chaser at the same radius but ahead (x₀ = 0, y₀ = 500) stays put for ever.

The example in `main()`: n = 1.1314 × 10⁻³ rad/s (period 92.56 min) for a = 6778 km. From (−50, −1000) m at rest, after half an orbit the chaser is at (−350.0, −57.5) m moving at 0.679 m/s along-track; after one orbit it is back at x = −50 but at y = +885 m, 1885 m ahead of where it started. Our exponential agrees with the closed form to 5.5 × 10⁻¹⁰ and SciPy's to 1.3 × 10⁻¹¹ at half an orbit; det Φ = 1.000000000000; Φ(P/2)Φ(P/2) − Φ(P) = 4.6 × 10⁻¹³.

### 2.6 Two-impulse rendezvous as a linear solve

Partition Φ(T) into 2×2 blocks, Φ = [[Φ_rr, Φ_rv], [Φ_vr, Φ_vv]]. An impulse Δ**v**₁ at t = 0 changes the velocity to **v**₀ + Δ**v**₁; the position at T is **r**(T) = Φ_rr**r**₀ + Φ_rv(**v**₀ + Δ**v**₁). Demanding **r**(T) = 0:

**Φ_rv Δv₁ = −(Φ_rr r₀ + Φ_rv v₀)**,   (6)

a 2×2 linear system (`P01.13`). The arrival velocity Φ_vr**r**₀ + Φ_vv(**v**₀ + Δ**v**₁) is cancelled by Δ**v**₂. The system is singular when det Φ_rv = 0; from the closed form det Φ_rv = [8(1 − cos nT) − 3nT sin nT]/n², which vanishes at nT = 2π: after exactly one orbit every chaser returns to the same radial position whatever its velocity, so no first burn can fix the position there. `rendezvous_delta_v` raises in that case. For the example: a quarter-orbit transfer costs 1.494 m/s (0.620 + 0.423 then 0.676 + 0.310 by component), half an orbit 0.518 m/s, three quarters 0.205 m/s. Slower is cheaper, down to the singularity, and the ascending branch beyond one orbit is checked in the tests at 1.3 orbits. The arrival position error in the model is 10⁻¹³ m: the solve is exact, and everything that goes wrong in a real rendezvous comes from the model (2.7).

### 2.7 Where the model stops being valid

- **Linearisation.** (5) drops terms of order (separation/a)². At 1 km and a = 6778 km the relative error is about 10⁻⁴ per orbit, a few decimetres over a 1 km transfer: fine for planning, not for the final approach, which uses relative navigation sensors and repeated corrections.
- **Circular target orbit.** For eccentricity e the equations acquire periodic coefficients (the Tschauner–Hempel equations); at e = 0.001 the error is similar to the linearisation error.
- **Perturbations.** J₂ affects chaser and target almost equally and mostly cancels in the relative motion; differential drag (different ballistic coefficients) does not, and at 400 km can be millimetres per second squared, comparable to the thrust arc in `main()`.
- **Impulsive burns.** A real burn lasts seconds to minutes; the finite-thrust correction is the subject of extension 3 and is small when the burn is short compared with the orbit.
- **The matrix is given.** In practice n is known to about 10⁻⁷ relative from tracking, which over a half-orbit transfer is a few centimetres: negligible against the linearisation error.

## 3. Build the computational model

Open `starter/main.py`. `iss_orbit`, `cw_matrix` and `cw_rhs` are given.

### Step 1: three exponentials

```python
def expm_taylor(A, n_terms=30):
    term, total = np.eye(len(A)), np.eye(len(A))
    for k in range(1, n_terms):
        term = term @ A / k                 # A^k / k! from A^{k-1} / (k-1)!
        total = total + term
    return total

def expm_scaling_squaring(A, n_terms=20):
    norm = np.linalg.norm(A, 1)
    s = max(0, int(np.ceil(np.log2(norm / 0.5)))) if norm > 0.5 else 0
    E = expm_taylor(A / 2 ** s, n_terms)
    for _ in range(s):
        E = E @ E
    return E
```

`expm_eig` is `V @ np.diag(np.exp(lam)) @ np.linalg.inv(V)` with a `cond(V)` check that raises `ValueError`; return the real part (the imaginary part is rounding noise when A is real). `state_transition(A, t)` is `expm_scaling_squaring(A * t)`; `propagate` applies it at each time; `is_stable` tests `eigvals(A).real < 0`.

### Step 2: forced response

`discretize(A, B, dt)` builds the augmented matrix (4) of size n + m, exponentiates it once, and slices A_d and B_d out. `forced_response_zoh` loops (3) over the columns of **u**. Accept **u** either as shape (m, N) or (N,) for a single input.

### Step 3: Clohessy–Wiltshire

`cw_stm_closed_form(n, t)` types in the four rows of 2.5. `rendezvous_delta_v` partitions it, checks det Φ_rv against 10⁻¹²‖Φ_rv‖², solves (6) with `np.linalg.solve`, and returns (Δ**v**₁, Δ**v**₂ = −arrival velocity). `modal_decomposition` returns (λ, V, V⁻¹) from `np.linalg.eig` for the diagonalisable systems of `P02.04`.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-05-matrix-exponential`. The tests:

1. **The series.** Exact for a nilpotent matrix in two terms; equal to the scalar exponential; converges to `expm` at ‖A‖ ≈ 0.3 and is off by more than 1 at ‖A‖ = 144 with 30 terms.
2. **Scaling and squaring** against `scipy.linalg.expm` to 10⁻⁹ on random 2×2, 3×3 (‖A‖ ≈ 5), 4×4 (‖A‖ ≈ 40) and 6×6 matrices.
3. **Eigen method** equal to `expm` for random and complex-eigenvalue matrices, real-valued, and raising for Clohessy–Wiltshire and for [[1, 1], [0, 1]].
4. **Commutation**: e^{A+B} = e^A e^B for A and B = 2A + ½I, and not for the nilpotent pair.
5. **Φ properties** on four matrices: Φ(0) = I, Φ(0.7)Φ(1.9) = Φ(2.6) in both orders, Φ(−t)Φ(t) = I, det Φ = e^{tr A·t}, and dΦ/dt = AΦ by central differences to 10⁻⁷.
6. **Propagation**: equal to `solve_ivp` (rtol 10⁻¹¹) and to the `P02.04` eigen closed form on the short-period mode; stability verdicts for a spiral, a saddle, an undamped oscillator and CW.
7. **Discretisation**: A_d = e^{Ah}; B_d equals the quadrature of ∫e^{As}ds B to 10⁻⁷; the small-h limits; the zero-order-hold response of a random piecewise-constant force equals hold-by-hold `solve_ivp` to 10⁻⁹; a constant force reaches the static deflection F/k.
8. **Clohessy–Wiltshire**: eigenvalues 0, 0, ±in; closed form = our Φ = SciPy's at 10 s, 600 s, P/2, P and 5000 s; det = 1; closed form = `solve_ivp` over one orbit; the −12π x₀ drift; a co-orbital chaser stays put; rendezvous arrives at rest for T = 0.25, 0.5, 0.75 and 1.3 orbits, confirmed by `solve_ivp`; T = P raises; the thrust arc by zero-order hold equals `solve_ivp` to 10⁻⁷ m.

### Compare

`scipy.linalg.expm` is the production tool: Padé(13) with scaling and squaring and a careful choice of s from norm estimates, accurate to near machine precision for well-scaled matrices. Our scaling-and-squaring differs from it by 10⁻¹⁰ on the Clohessy–Wiltshire half-orbit and by 10⁻⁹ relative on random 6×6 matrices, the difference being our cruder Taylor core and fixed threshold ½. The hand closed form is better than both (it is exact, the test tolerance is set by the oracle), which is the usual situation in orbital mechanics: closed-form transition matrices exist for the problems that matter most and are used in flight software, with `expm` reserved for the matrices that have none. `solve_ivp` reproduces everything to 10⁻⁷ with rtol 10⁻¹² but needs a new integration for every initial condition, while Φ(T), computed once, answers for all of them: that is why the rendezvous is a linear solve and not an optimisation. The eigen method is the one to use for *understanding* (which mode is slow, which is unstable) and the one that fails on exactly the matrices where understanding is most needed.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Δv against transfer time.** For the chaser at (−50, −1000) m at rest, compute total Δv for T from 0.1 to 1.9 orbits in steps of 0.05, skipping T = 1 orbit. Predict the shape: large for short transfers (the chaser must be thrown at the target), falling toward one orbit, a singular spike at one orbit, and a second, cheaper window beyond it. Record the T of the cheapest transfer in the sweep.
2. **Taylor breakdown.** For A = A_sp (the `P02.04` short-period matrix) × t with t = 1, 3, 10, 30, 100 s, compute the error of `expm_taylor` with 20, 30 and 60 terms against `expm`, and the number of squarings `expm_scaling_squaring` uses. Predict the t beyond which 30 terms fail (where ‖At‖₁ passes roughly 15, i.e. where the largest term ‖At‖^k/k! exceeds 10¹⁶ times the size of the answer).
3. **Hold length.** For the 10-minute thrust arc, compare `forced_response_zoh` with holds of 1, 10, 60 and 600 s against `solve_ivp` with a constant acceleration. Predict: all equal to rounding, because the acceleration *is* constant; then make the acceleration ramp linearly from 0 to 0.4 mm/s², sample it at the midpoint of each hold, and predict that the error now scales as h² (and as h if you sample at the left end instead).
4. **The co-orbital trap.** Start a chaser at (0, −1000) m with a small radial velocity (1 cm/s) and propagate one orbit; predict where it ends up from the closed form (x returns, y shifts by −2(1 − cos 2π)/n · ẋ₀ = 0).

Raw data format: CSV tables (`T_orbits, dv1_x, dv1_y, dv2_x, dv2_y, total_m_per_s`; `t_s, norm_At, err_20, err_30, err_60, squarings`; `hold_s, end_error_m`), each with a header naming the fixed parameters.

What you should find: the Δv curve falls from 1.49 m/s at a quarter orbit through 0.52 at a half and 0.21 at three quarters, spikes at one orbit, and the cheapest transfer in the sweep is 1.75 orbits at 0.13 m/s (1.25 orbits costs 0.24, 1.5 orbits 0.39); plain Taylor with 30 terms is fine at t = 1 s (error 10⁻¹⁵), marginal at 3 s (4 × 10⁻⁸ with ‖At‖₁ = 14.4), and useless at 10 s (error 10⁸ against true entries of 10⁻⁴) and 30 s (10²²), while scaling and squaring uses 4, 5, 7, 9 and 10 squarings for t = 1, 3, 10, 30, 100 s; the constant-thrust arc is exact for every hold length (end error below 10⁻⁹ m) and the ramped one gives end errors of 2.3 × 10⁻⁵, 2.3 × 10⁻³ and 8.4 × 10⁻² m for 1, 10 and 60 s midpoint holds (the h² law), against 0.05, 0.5 and 3.0 m for left-end holds.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| n (rad/s) | √(3.986 × 10¹⁴ / 6.778 × 10⁶)³... = 1.1314 × 10⁻³ | 1.131401 × 10⁻³ | 0 | exact |
| drift per orbit (m) | −12π x₀ = 1885 | 1885.0 | 0 | exact |
| Φ(P/2): ours vs closed form | 0 | 5.5 × 10⁻¹⁰ | — | 10 squarings × rounding |
| Φ(P/2): SciPy vs closed form | 0 | 1.3 × 10⁻¹¹ | — | Padé core |
| det Φ(P/2) | 1 | 1.000000000000 | < 10⁻¹² | tr A = 0 |
| propagate vs `solve_ivp`, one orbit (m) | — | 9.5 × 10⁻⁸ | — | integrator rtol 10⁻¹² |
| Δv, 0.25/0.5/0.75 orbit (m/s) | — | 1.494 / 0.518 / 0.205 | — | 2×2 solve, exact |
| arrival position error (m) | 0 | ≤ 1.1 × 10⁻¹³ | — | rounding |
| thrust arc, ZOH vs `solve_ivp` (m) | 0 (ZOH is exact) | 4.8 × 10⁻¹⁰ | — | integrator |
| Taylor, 30 terms, ‖A‖ = 144 | fails | error 1.2 × 10²² | — | cancellation |

The numerics are exact to the limits of floating point, so the validation question is about the *model*: how far is a Clohessy–Wiltshire plan from what a real chaser would do? The honest answer is "a few decimetres per kilometre per orbit from linearisation, plus whatever differential drag and burn execution errors contribute", and operational rendezvous handles this not by a better transition matrix but by re-planning: the same linear solve is repeated every few minutes from the measured relative state, each time with a shorter transfer and a smaller residual. Φ makes that cheap enough to do in flight software at every step. The second validation is of the exponential itself, and here the lesson is that three independent constructions agreeing to 10⁻⁹ (series, Padé in SciPy, hand calculus) is strong evidence, while any one of them alone is not: the Taylor series "converges" to garbage at ‖A‖ = 144 without raising an error.

Improvements worth making: the Padé core (extension 1), which halves the squarings for the same accuracy; the out-of-plane equation (extension 2); finite burns (extension 3). Improvement worth *not* making: a smaller threshold than ½ in the scaling step: it adds squarings and rounding without gaining accuracy.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of the five Φ properties, the commutation counterexample, (2), (4) and the first column of the Clohessy–Wiltshire matrix; `state_transition.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a transfer time for the example chaser, with the Δv, the sensitivity to a 1 % error in n and the operational risk of a plan near the singular time.

### Adversarial self-checks

- **What if the target orbit were at 800 km instead of 400 km?** n falls by 8 %, the period rises to 100.9 min, and the Δv for a given fraction of an orbit falls in proportion to n (the velocities in Φ_rv scale with n). Check with `iss_orbit` replaced.
- **What if the chaser starts 10 km behind instead of 1 km?** Everything scales linearly (the model is linear), so Δv is ten times larger; but the linearisation error is now 100 times larger in absolute terms (it is quadratic in separation), several metres per orbit. Which part of the answer do you trust?
- **Which assumption fails first in the real world?** Differential drag at 400 km for dissimilar vehicles, then burn execution error (typically 1 % of Δv, i.e. 5 mm/s for the half-orbit plan, which over the 46 min coast is a 14 m miss, against 4 m for a 10⁻⁶ m/s² differential drag and 0.15 m for linearisation); linearisation only at larger separations.
- **What measurement could falsify the model?** Relative GPS or lidar range during a coasting arc: a residual growing as t² where the model predicts the closed-form path indicates an unmodelled constant acceleration (drag or venting); a residual at twice the orbital frequency indicates eccentricity.
- **Does the result survive a change of solver?** Series, scaling and squaring, SciPy's Padé, the hand closed form and `solve_ivp` agree to 10⁻⁷ m. The eigen method does not survive at all on this matrix, and says so.
- **Which decision changes if the uncertainty doubles?** If burn error doubles to 2 %, the half-orbit plan misses by about 29 m and the three-quarter plan by 17 m (less Δv, so less absolute velocity error, despite the longer coast), while a drag-type error grows as T² and favours the shorter transfer. Either way a mid-course correction becomes necessary, and its cost (a few cm/s) is comparable to the Δv difference between the plans: state which plan you would fly and why.

## Check yourself

1. Write the series that defines e^{A} and explain why it converges for every square matrix.
2. State the five properties of Φ(t) = e^{At} and say which one implies that trajectories of a linear system cannot cross.
3. Why is e^{A+B} not equal to e^A e^B in general? Give the condition under which it is.
4. Which of the three computation methods fails on the Clohessy–Wiltshire matrix, and why is that failure a statement about the physics?
5. Write the zero-order-hold update x_{k+1} = A_d x_k + B_d u_k and the augmented matrix whose exponential yields A_d and B_d.
6. A chaser is 100 m below the target at rest. How far along-track has it drifted after one orbit, and in which direction?
7. Why is the two-impulse rendezvous a linear solve, and at what transfer time does it have no solution?
8. State the condition on the eigenvalues for every solution of ẋ = Ax to decay, and what happens when an eigenvalue with zero real part is defective.
