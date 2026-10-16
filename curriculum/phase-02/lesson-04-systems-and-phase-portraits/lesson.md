# P02.04  Systems of ODEs & Phase Portraits

> Two numbers, the trace and the determinant, decide whether a linear system spirals home, runs away along a line, or circles for ever. Draw the plane and you can see the answer before you solve anything.

| | |
|---|---|
| **Learning outcomes** | 1. Write a planar linear system d**x**/dt = A**x**, derive its characteristic polynomial λ² − τλ + Δ from the trace and determinant, and classify the fixed point (saddle, node, spiral, center, degenerate cases) from (τ, Δ) alone. 2. Solve the system in closed form from eigenvalues and eigenvectors, including the complex and repeated-root cases, and read natural frequency, damping ratio, period and time-to-half from the eigenvalues. 3. Implement and test a phase-portrait module: trace–determinant classification, hand-computed eigenstructure, closed-form solutions, RK4 trajectories, vector field, nullclines and the conserved quadratic of a center, checked against `numpy.linalg`, `scipy.linalg.expm` and `solve_ivp`. 4. Decide whether an aircraft's short-period and phugoid modes are acceptable and what happens to the short period as the centre of gravity moves aft. |
| **Prerequisites** | `P02.02` (second-order equations as 2×2 systems, ζ and ω_n), `P02.01` (RK4), `P01.16` (eigenvalues and eigenvectors of a 2×2 matrix, A = VΛV⁻¹), `P01.13` (2×2 linear solves, determinants), `P00.03`, `P00.07`. Skill check: you can find the eigenvalues and eigenvectors of a 2×2 matrix by hand and sketch e^{σt} cos ωt. |
| **Engineering question** | A jet transport in cruise has the short-period matrix A = [[−0.8, 1], [−4.0, −1.2]] s⁻¹ for angle of attack and pitch rate, and a phugoid with V₀ = 230 m/s and L/D = 17. How fast does a 5° gust-induced angle-of-attack disturbance die away, what are the two modes' periods and damping ratios, and at what aft shift of static stability (M_α) does the aircraft become a saddle that diverges? |
| **Assumptions and limits** | Small perturbations about trimmed, straight, level flight, so the equations are linear with constant coefficients; the short-period approximation holds speed constant (valid because the mode is over the speed response); Lanchester's phugoid holds angle of attack constant and neglects pitching dynamics and compressibility (it gives the period well and the damping only roughly); no control inputs; no stall (the linear lift slope fails beyond about 12° angle of attack). A 2×2 model has at most one mode pair: the real longitudinal dynamics are 4×4 and couple the two modes weakly. |
| **Analytical work** | Trace and determinant as the coefficients of the characteristic polynomial; the trace–determinant plane and its regions; eigenvectors by hand as columns of the adjugate; real solutions from a complex eigenpair; the repeated-root solution e^{λt}(I + (A − λI)t) via Cayley–Hamilton; Lanchester's phugoid from energy; the conserved quadratic form of a traceless system; nullclines and the direction of crossing. |
| **Computational deliverable** | `starter/main.py`: `trace_determinant`, `eigen_structure`, `classify`, `modal_parameters`, `solve_linear`, `rk4_trajectory`, `vector_field`, `nullcline_slopes`, `center_invariant`, `phase_portrait` (`short_period_matrix`, `phugoid_matrix`, `second_order_matrix` are given). Tests prove: the hand eigenstructure matches `numpy.linalg.eig` and satisfies A**v** = λ**v** for random matrices; all eleven fixed-point types are classified correctly and consistently with the eigenvalues on 40 random matrices; the closed form equals `expm(At)`**x**₀ and satisfies the equation in every regime including the degenerate and phugoid cases; RK4 converges at order 4 and `solve_ivp` agrees to 10⁻⁸; saddle trajectories follow the eigen-directions and every trajectory ends parallel to the unstable one; the center's quadratic is conserved to 10⁻⁹ and equals −2E/m for the undamped oscillator; the spiral's zero crossings and decay match the eigenvalues; nullclines are crossed vertically and horizontally. |
| **Practical activity** | **computational experiment**: walk the short-period matrix across the trace–determinant plane by varying M_α from −4 to +3 s⁻² and record the fixed-point type, the eigenvalues and the time to half or double; sweep the phugoid over speed and L/D. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `phase_portraits.png`, the hand classification and eigenvectors for the three short-period cases, the M_α sweep table with the two transition values you predicted, and a short note on whether the modes meet a handling-qualities-style requirement you state. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the trace–determinant derivation and the complex-pair real form are complete; all public tests pass; the hand eigenvalues of the three short-period matrices match the code to three figures; the note states the two M_α transitions and a judgment with numbers. |
| **References** | `@mit-ocw-18-03-differential-equations` Lecture 24 "Introduction to First-order Systems of ODEs", Lecture 25 "Homogeneous Linear Systems with Constant Coefficients", Lecture 26 "Continuation: Repeated Real Eigenvalues", Lecture 27 "Sketching Solutions of 2x2 Homogeneous Linear System with Constant Coefficients" (required-free); notes 32 "Linear systems and matrices", 33 "Eigenvalues, eigenvectors", 34 "Complex or repeated eigenvalues", 35 "Qualitative behavior of linear systems; phase plane" (supplementary-free). No OpenStax coverage: Calculus Volume 3 stops at single second-order equations, so `P02.02` is the bridge. |
| **Extensions** | ★ Build the 4×4 longitudinal matrix by combining the two approximations with the coupling terms set to small values, and show the four eigenvalues split into the two pairs found here. ★ Add a pitch-rate feedback q → elevator (M_q → M_q + k) and plot the short-period eigenvalues as k varies (a root locus, previewed). ★ Draw the trace–determinant plane itself and mark every case in `main()` on it. |
| **Time** | unmeasured |

## 1. The engineering problem

When an aircraft in cruise meets a gust, two things happen on two time scales. Within a few seconds the nose pitches, the angle of attack overshoots and settles: the **short-period mode**. Over the next minute or two the aircraft trades speed for height and back, climbing and sinking in a long, slow oscillation: the **phugoid**. Each mode is described, to a good approximation, by *two* coupled first-order equations, angle of attack and pitch rate for the short period, speed and flight-path angle for the phugoid. Whether each mode dies away, how fast, and with how much overshoot are the first questions of flight dynamics, and they decide whether the aircraft can be flown by hand.

Mathematically both are **planar linear systems**, d**x**/dt = A**x** with A a constant 2×2 matrix. This lesson is about such systems in general: how to classify the behaviour near the fixed point at the origin from the trace and determinant of A without solving anything, how to solve exactly when you need numbers, and how to draw the **phase portrait**, the picture of all trajectories at once. The same picture describes a damped oscillator, a predator–prey balance, a chemical reactor and a spacecraft's attitude; the aircraft modes are the example because their matrices are realistic and their eigenvalues have names a pilot would recognise.

The concrete problem: for the short-period matrix A = [[−0.8, 1], [−4.0, −1.2]] (units s⁻¹, 1, s⁻², s⁻¹), classify the mode, find its frequency, damping ratio, period and time to half amplitude, and follow a 5° angle-of-attack disturbance until it is within 0.5°. Do the same for the phugoid at 230 m/s and L/D = 17. Then move the static-stability derivative M_α from −4 toward positive values (the centre of gravity moving aft) and find where the stable spiral becomes a node and where it becomes a saddle. Success criteria, written first: the hand eigenstructure matches `numpy.linalg.eig` to 10⁻⁹; the closed-form solution matches `scipy.linalg.expm` to 10⁻⁸ for every type of fixed point; RK4 converges at order 4; and the transition values of M_α predicted from τ² = 4Δ and Δ = 0 agree with the sweep.

## 2. Understand and derive

### 2.1 From one second-order equation to a system, and back

`P02.02` wrote m ẍ + c ẋ + k x = 0 as d/dt (x, ẋ) = [[0, 1], [−k/m, −c/m]] (x, ẋ): a **first-order system** with state vector **x** = (x₁, x₂). The reverse also works for any 2×2 system: differentiate ẋ₁ = a x₁ + b x₂ once more and eliminate x₂ to get ẍ₁ − (a + d) ẋ₁ + (ad − bc) x₁ = 0 (18.03 Lecture 24). The coefficients are the **trace** τ = a + d and the **determinant** Δ = ad − bc of A, so every planar linear system *is* a second-order equation with ω_n² = Δ and 2ζω_n = −τ, and everything learned about ζ and ω_n carries over.

### 2.2 The eigenvalue solution and the characteristic polynomial

Try **x** = **v** e^{λt} (Lecture 25). Then λ**v** e^{λt} = A**v** e^{λt}, so A**v** = λ**v**: λ is an **eigenvalue** and **v** an **eigenvector** (`P01.16`). Non-trivial **v** needs det(A − λI) = 0:

**λ² − τλ + Δ = 0,   λ = τ/2 ± √(τ²/4 − Δ)**.   (1)

The eigenvalues sum to τ and multiply to Δ. For an eigenvalue λ the eigenvector is any non-zero column of the adjugate of (A − λI), that is (b, λ − a) or (λ − d, c); `eigen_structure` picks the longer of the two and normalises it, so the code never calls `numpy.linalg.eig` except to check itself. For the short period: τ = −2.0, Δ = (−0.8)(−1.2) − (1)(−4.0) = 4.96, λ = −1 ± √(1 − 4.96) = **−1 ± 1.99i s⁻¹**.

With two independent eigenvectors the general solution is **x** = c₁**v**₁e^{λ₁t} + c₂**v**₂e^{λ₂t}, the constants fixed by **x**(0): **c** = V⁻¹**x**₀ with V = [**v**₁ **v**₂]. Equivalently **x**(t) = V diag(e^{λt}) V⁻¹ **x**₀, which is `solve_linear`, and `P02.05` will name the matrix in the middle e^{At}.

### 2.3 Complex eigenvalues: the real solution

When τ² < 4Δ, λ = σ ± iω with σ = τ/2 and ω = √(Δ − τ²/4), and the eigenvectors are a conjugate pair **v**, **v̄** (notes 34). The complex solution **v**e^{(σ+iω)t} has real and imaginary parts

e^{σt}[Re **v** cos ωt − Im **v** sin ωt]  and  e^{σt}[Re **v** sin ωt + Im **v** cos ωt],

each a real solution (the equation is real, so real and imaginary parts solve it separately), and together they span the solution space. The trajectory is an ellipse in the (Re **v**, Im **v**) frame traversed at angular rate ω, scaled by e^{σt}: a **spiral** inward for σ < 0, outward for σ > 0, a closed **center** for σ = 0. The engineer's names: **ω_n = √Δ**, **ζ = −τ/(2√Δ)**, damped frequency ω, **period 2π/ω**, **time to half amplitude ln 2/|σ|** (or time to double if σ > 0). For the short period: ω_n = 2.227 rad/s, ζ = 0.449, period 3.16 s, half-amplitude time 0.69 s. A 5° disturbance peaks at a pitch rate of −5.15°/s after 0.56 s and is inside 0.5° after 2.0 s: a mode a pilot does not have to fly, which is what "acceptable short period" means.

The computer does not need the real form (it exponentiates complex numbers and takes the real part at the end), but you do, to read ω and σ off a matrix, and the test checks that the zero crossings of α are spaced by half the period and that the amplitude falls by e^{σT} per period.

### 2.4 Repeated eigenvalues: the degenerate node and the star

When τ² = 4Δ there is one eigenvalue λ = τ/2. If A = λI every vector is an eigenvector and trajectories are straight rays: a **star**. Otherwise there is a single eigenvector line and the second solution is not another exponential. Cayley–Hamilton says A satisfies its own characteristic equation, (A − λI)² = 0, so the series for e^{At} terminates after two terms (Lecture 26):

**x**(t) = e^{λt}[I + (A − λI)t] **x**₀.   (2)

The factor t is the same one that appeared in the critically damped oscillator (`P02.02` eq. 7) and at resonance (`P02.03` eq. 5): wherever a root is repeated, a t appears. The trajectories of this **degenerate node** all become tangent to the single eigenvector line. `solve_linear` switches to (2) when the eigenvalues coincide to 10⁻⁹ relative, which also handles the star (A − λI = 0).

### 2.5 The trace–determinant plane

Collect the cases (Lecture 27, notes 35). For an **isolated** fixed point (Δ ≠ 0):

```
                 Δ
                 ▲          τ² = 4Δ  (parabola)
   stable        │        unstable
   spirals   ╲   │   ╱    spirals             center on the Δ axis (τ = 0, Δ > 0)
              ╲  │  ╱
   stable      ╲ │ ╱      unstable            degenerate nodes / stars ON the parabola
   nodes        ╲│╱       nodes
  ───────────────┼─────────────────▶ τ
            saddles (Δ < 0)                   non-isolated fixed points ON the τ axis (Δ = 0)
```

- **Δ < 0: saddle.** Two real eigenvalues of opposite sign. Trajectories approach along the stable eigenvector, leave along the unstable one, and every other trajectory ends up parallel to the unstable direction. Always unstable.
- **Δ > 0, τ² > 4Δ: node.** Two real eigenvalues of the same sign: stable if τ < 0, unstable if τ > 0. Trajectories are tangent to the *slow* eigenvector (smaller |λ|) at the origin.
- **Δ > 0, τ² < 4Δ: spiral** (τ ≠ 0) or **center** (τ = 0).
- **Δ > 0, τ² = 4Δ:** degenerate node, or star if A = λI.
- **Δ = 0:** a line of fixed points (one eigenvalue is zero); trajectories are straight lines onto it.

The determinant decides saddle or not; the trace decides stable or not; the parabola decides oscillatory or not. That is the whole classification, and `classify` is twenty lines that implement it. The aircraft cases from `main()`:

```
case                                                   tau   Delta  type                     eigenvalues
short period (jet transport, cruise)                -2.000   4.960  stable spiral            [-1.-1.99j -1.+1.99j]
phugoid (V0 = 230 m/s, L/D = 17)                    -0.005   0.004  stable spiral            [-0.0025-0.0603j -0.0025+0.0603j]
short period, CG too far aft (M_alpha = +2)         -2.000  -1.040  saddle                   [-2.4283  0.4283]
short period with tau forced to 0 (M_q = +0.8): center  0.000   3.360  center               [0.-1.833j 0.+1.833j]
undamped oscillator m x'' + k x = 0 (center)         0.000   4.000  center                   [0.-2.j 0.+2.j]
critically damped oscillator (degenerate node)      -4.000   4.000  stable degenerate node   [-2. -2.]
overdamped oscillator (stable node)                 -5.000   4.000  stable node              [-4. -1.]
```

### 2.6 The aircraft matrices

**Short period.** With speed frozen, the perturbation equations of pitch are α̇ = (Z_α/V₀) α + q and q̇ = M_α α + M_q q. Z_α/V₀ < 0: more angle of attack means more lift, which curves the flight path up and so *reduces* α (the aircraft "weathercocks" into the relative wind). M_α < 0 is **static stability**: a nose-up disturbance produces a nose-down moment, and its magnitude depends on how far the centre of gravity is ahead of the neutral point. M_q < 0 is **pitch damping** from the tail moving through the air. Δ = (Z_α/V₀)M_q − M_α and τ = Z_α/V₀ + M_q. As the centre of gravity moves aft, M_α rises toward zero and beyond. With Z_α/V₀ = −0.8 and M_q = −1.2 fixed: Δ = 0.96 − M_α, τ = −2. The spiral becomes a node when τ² = 4Δ, i.e. 4 = 4(0.96 − M_α), **M_α = −0.04 s⁻²**; the node becomes a saddle when Δ = 0, **M_α = +0.96 s⁻²**. Beyond that the aircraft diverges along the unstable eigenvector; at M_α = +2 the eigenvalues are −2.43 and +0.428 s⁻¹, the unstable direction is (0.631, 0.775) in (α, q), and the disturbance doubles every 1.62 s. Part 4 walks the whole path.

**Phugoid** (Lanchester, 1908). Hold α, hence the lift coefficient, constant. A speed excess u above trim raises lift by the factor (1 + u/V₀)² ≈ 1 + 2u/V₀, and the excess lift mg·2u/V₀ curves the path upward: θ̇ = (2g/V₀²) u, i.e. d(θ)/dt = (2g/V₀)(u/V₀). Climbing at angle θ, gravity decelerates the aircraft: u̇ = −gθ, and the drag excess at higher speed, 2D u/V₀ = (2mg/(L/D)) u/V₀, removes energy: u̇ = −gθ − (2g/(L/D))(u/V₀). In the state (u/V₀, θ):

A_ph = [[−2g/(V₀ L/D), −g/V₀], [2g/V₀, 0]],  Δ = 2g²/V₀²,  τ = −2g/(V₀ L/D),

so **ω_n = √2 g/V₀** and **ζ = 1/(√2 L/D)**, independent of mass, wing area and altitude. For V₀ = 230 m/s and L/D = 17: ω_n = 0.0603 rad/s, period **104 s**, ζ = 0.0416, half-amplitude time 276 s. The period formula is remarkably good for real aircraft (it predicts 60 s for a 130 m/s trainer); the damping is only roughly right, because compressibility and the pitching moment from speed change (M_u) add damping the model omits. The phugoid is slow enough that a pilot, or an autopilot, damps it without noticing, which is why a ζ of 0.04 is acceptable here and would not be for the short period.

### 2.7 Nullclines and the conserved quadratic

Two more tools for reading a portrait. The **x₁-nullcline** is the line a x₁ + b x₂ = 0 where ẋ₁ = 0, which trajectories cross *vertically*; the **x₂-nullcline** c x₁ + d x₂ = 0 is crossed *horizontally*. They meet only at the origin when Δ ≠ 0 and tell you the sense of rotation at a glance. For a **center** (τ = 0) something is conserved: with A = [[a, b], [c, −a]], the quadratic

**Q = c x₁² − 2a x₁x₂ − b x₂²**   satisfies   dQ/dt = −2b(a + d) x₂² = 0,   (3)

(expand each term, use ẋ = A**x** and a + d = 0), so the orbits are the level curves of Q, ellipses when Δ = −a² − bc > 0. For the undamped oscillator A = [[0, 1], [−k/m, 0]] this gives Q = −(k/m)x² − v² = −(2/m)E: the conserved quadratic is the energy. For a damped system the same Q decays, and the test checks that it falls to 1 % along the short period.

### 2.8 Where the model stops being valid

- **Linearity.** The aerodynamic derivatives are the slopes at trim; beyond about 12° angle of attack the lift curve bends and M_α itself changes. A saddle's exponential divergence is a statement about the first few seconds only; what happens next is nonlinear (`P02.06`).
- **Decoupled approximations.** The real longitudinal system is 4×4; the two 2×2 models are good because the modes differ in frequency by a factor of 40 here. For an aircraft with a slow short period (a sailplane) they interact.
- **Constant coefficients.** The derivatives depend on speed, altitude and configuration: the matrix for approach is not the matrix for cruise.
- **No controls.** M_q includes only the natural pitch damping; a stability augmentation system adds to it (extension 2), which is how an aircraft with M_α > 0 (relaxed static stability) is flown at all.

## 3. Build the computational model

Open `starter/main.py`. The three matrix builders are given; everything else is yours.

### Step 1: trace, determinant, eigenstructure by hand

```python
def eigen_structure(A):
    tau, delta = trace_determinant(A)
    disc = tau * tau / 4.0 - delta
    sq = np.sqrt(complex(disc))
    lam = np.array([tau / 2.0 - sq, tau / 2.0 + sq])                      # (1)
    a, b, c, d = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
    vecs = []
    for l in lam:
        v1, v2 = np.array([b, l - a]), np.array([l - d, c])              # columns of adj(A - l I)
        v = v1 if np.linalg.norm(v1) >= np.linalg.norm(v2) else v2
        if np.linalg.norm(v) == 0.0:                                     # A = l I: any vector will do
            v = np.array([1.0, 0.0]) if not vecs else np.array([0.0, 1.0])
        vecs.append(v / np.linalg.norm(v))
    ...                                                                   # return real arrays when the eigenvalues are real
```

Use complex arithmetic for the square root so that a negative discriminant is not an error, and return real arrays when the eigenvalues are real so that the saddle's eigen-directions can be plotted. The test compares with `numpy.linalg.eig` on random matrices and checks A**v** = λ**v** directly.

### Step 2: classification

`classify(A)` follows 2.5 in order: Δ < 0 → saddle; |Δ| ≈ 0 → non-isolated; then the sign of τ (negative, positive, zero → stable, unstable, center); then the discriminant τ² − 4Δ (zero → star or degenerate node, positive → node, negative → spiral or center). Use tolerances scaled by the matrix entries, not absolute ones. `modal_parameters` returns ω_n = √Δ, ζ = −τ/(2√Δ), ω_d, period and ln 2/|τ/2|, and raises for real eigenvalues.

### Step 3: solutions

```python
def solve_linear(A, x0, t):
    lam, V = eigen_structure(A)
    if abs(lam[0] - lam[1]) > 1e-9 * max(1.0, np.abs(lam).max()):
        coef = np.linalg.solve(V.astype(complex), x0.astype(complex))     # c = V^{-1} x0
        E = np.exp(np.outer(t, lam))                                      # e^{lambda_k t}
        return np.real((E * coef) @ V.T.astype(complex))                  # sum_k c_k e^{lambda_k t} v_k
    l = float(np.real(lam[0]))
    N = A - l * np.eye(2)
    return np.exp(l * t)[:, None] * (x0[None, :] + t[:, None] * (N @ x0)[None, :])   # (2)
```

`rk4_trajectory` is the stepper of `P02.01` with f(**x**) = A**x**; `phase_portrait` calls `solve_linear` from each start; `vector_field` returns (U, V) on a mesh for `plt.streamplot`; `nullcline_slopes` returns −a/b and −c/d (∞ for a vertical line); `center_invariant` is (3).

### Verify

Run `python -m pytest curriculum/phase-02/lesson-04-systems-and-phase-portraits`. The tests:

1. **Analytic structure.** τ and Δ are the trace and determinant; both eigenvalues are roots of λ² − τλ + Δ; the hand eigenstructure matches `numpy.linalg.eig` and A**v** = λ**v** holds to 10⁻¹⁰ on seven matrices including random ones; the degenerate node returns one eigen-direction twice and the star returns any orthogonal pair.
2. **Classification.** Eleven named cases, one per region of the plane; on 40 random matrices the name agrees with the eigenvalues (complex ↔ spiral/center, opposite signs ↔ saddle, sign of the real part ↔ stable/unstable).
3. **Solutions.** For seven matrices (spiral, saddle, node, center, star, degenerate node, phugoid) the closed form equals `expm(At)`**x**₀ to 10⁻⁸ at 61 times and satisfies d**x**/dt = A**x** by central differences; RK4 orders are 4.0 ± 0.3; `solve_ivp` agrees to 10⁻⁸.
4. **Geometry.** Starts on a saddle's eigenvectors stay on them; any other start ends parallel to the unstable direction (cos > 1 − 10⁻⁶ at t = 12 s); the center's Q is conserved to 10⁻⁹ and equals −2E/m for the oscillator; along the damped spiral Q falls below 1 %; the spiral's zero crossings are spaced by half the period and its amplitude falls by e^{−ζω_nT} per period to 10⁻³.
5. **Flight dynamics.** Lanchester's ω_n = √2 g/V₀ and ζ = 1/(√2 L/D) hold for three (V₀, L/D) pairs; `second_order_matrix` reproduces `P02.02`'s roots; `modal_parameters` raises for a saddle.

### Compare

`numpy.linalg.eig` (LAPACK's Hessenberg–QR algorithm) and `scipy.linalg.expm` (Padé approximation with scaling and squaring, `P02.05`) are the oracles; our hand eigenstructure agrees with the first to 10⁻⁹ and our closed form with the second to 10⁻⁸. For a 2×2 matrix the hand formulas are exact and cheaper, and they are what you should use when you want to *see* the dependence on a parameter (2.6 would be impossible with a numerical eigensolver as a black box). For 4×4 and beyond there is no closed form and `eig` is the tool; `expm` matters when you want the solution at many times without stepping. `solve_ivp` agrees with both to 10⁻⁸ and is the only one of the three that will still work when `P02.06` makes the right-hand side nonlinear. We wrote our own classification because `eig` returns eigenvalues, not a verdict: the verdict is the engineer's, and the trace–determinant plane is how it is reached.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Walk across the plane.** With Z_α/V₀ = −0.8 and M_q = −1.2 fixed, vary M_α from −4 to +3 in steps of 0.25 s⁻² and record (τ, Δ), the classification, the eigenvalues and the time to half (or double) amplitude. Predict the two transition values from τ² = 4Δ and Δ = 0 before running (2.6 gives −0.04 and +0.96).
2. **Portraits at the transitions.** Draw the phase portrait at M_α = −2 (spiral), 0 (node, eigenvalues −1.2 and −0.8), +0.5 (node, −1.73 and −0.27) and +2 (saddle) from eight starts on a circle. Predict which eigenvector the node trajectories hug (the slow one).
3. **Phugoid sweep.** Period and ζ for V₀ = 60, 130, 230, 250 m/s and L/D = 10, 17, 20. Predict the scaling (period ∝ V₀, ζ ∝ 1/(L/D)).
4. **RK4 cost.** Steps needed for 10⁻⁶ accuracy over 8 s of short period against over 300 s of phugoid; predict from the periods.

Raw data format: CSV tables (`M_alpha, tau, Delta, type, re_lambda1, im_lambda1, re_lambda2, im_lambda2, t_half_or_double`; `V0, LD, period_s, zeta`), each with a header naming the fixed parameters.

What you should find: the type changes exactly where predicted; the doubling time of the saddle is 35 s at M_α = +1 (barely unstable, λ = +0.020 s⁻¹), 2.9 s at +1.5, 1.62 s at +2 and 0.93 s at +3, with the unstable direction tilting from 39° to 57° in the (α, q) plane; node trajectories bend onto the slow eigenvector; the phugoid period runs from 27 s at 60 m/s to 113 s at 250 m/s independent of L/D, and ζ from 0.071 at L/D = 10 to 0.035 at L/D = 20 independent of speed.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| short period τ, Δ | −2.0, 4.96 | −2.000, 4.960 | 0 | exact |
| eigenvalues (s⁻¹) | −1 ± 1.99i | −1 ± 1.98997i | 0 | exact |
| ω_n, ζ | 2.227, 0.449 | 2.2271, 0.4490 | 0 | exact |
| period, t_half (s) | 3.16, 0.69 | 3.157, 0.693 | 0 | exact |
| phugoid ω_n (rad/s) | √2 × 9.81/230 = 0.0603 | 0.0603 | 0 | exact |
| phugoid period (s) | 104 | 104.26 | 0 | exact |
| phugoid ζ | 1/(17√2) = 0.0416 | 0.0416 | 0 | exact |
| M_α transitions (s⁻²) | −0.04, +0.96 | −0.04, +0.96 (sweep) | 0 | exact |
| saddle eigenvalues at M_α = +2 | −1 ± √(1 + 1.04) = −2.428, 0.428 | −2.4283, 0.4283 | 0 | exact |
| RK4 (800 steps, 8 s) vs closed form | ∝ h⁴ | 2.5 × 10⁻¹⁰ | — | order 4 |
| `solve_ivp` vs closed form | — | 2.7 × 10⁻¹² | — | rtol 10⁻¹⁰ |

Every number is exact or explained, which is the easy part. The hard part is that the *inputs* are approximate: Z_α/V₀, M_α and M_q are aerodynamic derivatives estimated from wind-tunnel data or computed from geometry, typically to 10–20 %. A 20 % error in M_α moves Δ from 4.96 to 4.16 or 5.76, ζ from 0.449 to 0.49 or 0.42, the period by ±8 %: the *character* of the short period is robust. Near the saddle transition it is not: at M_α = +0.96 ± 0.2 the aircraft is either a slow node or a saddle doubling every 10 s, so an aircraft designed close to neutral stability needs its derivatives known far better than 20 %, or a stability augmentation system that makes M_q large enough that the margin does not matter. That is the engineering reading of the trace–determinant plane: distance from the boundaries is margin.

What validates these models is flight test: a pulse on the elevator excites the short period, the recorded α(t) and q(t) are fitted to e^{σt} cos(ωt − φ), and σ and ω are compared with the eigenvalues. For the phugoid a small speed change and hands-off flight give the period directly. The decoupled approximations typically agree with flight data to 10 % in frequency and worse in damping, which is where the 4×4 model (extension 1) earns its complexity.

Improvements worth making: the 4×4 longitudinal matrix; the pitch-rate feedback root locus (extension 2), which is the cheapest way to see how a saddle at M_α = +2 can be made a spiral; and a trace–determinant plot with all seven `main()` cases marked. Improvement worth *not* making: a finer RK4 step; the closed form is exact.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (1), the real form of the complex solution, (2) and (3); the hand eigenvalues and eigenvectors of the three short-period matrices; `phase_portraits.png`; the M_α sweep table with your predicted transitions; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating a requirement for each mode (for instance "short-period ζ between 0.3 and 2, phugoid stable") and whether this aircraft meets it, with the assumption you trust least.

### Adversarial self-checks

- **What if the pitch damping doubles (M_q = −2.4)?** τ = −3.2, Δ = 1.92 + 4 = 5.92, ζ = 0.66: still a spiral, more heavily damped. What if it goes to zero? τ = −0.8, ζ = 0.19, period 3.1 s with 11 overshoots to settle. Check both with the code.
- **What if M_α → −∞ (very forward CG)?** Δ → ∞, ω_n → ∞ with τ fixed: a fast, lightly damped short period (ζ → 0) that is unpleasant to fly for the opposite reason. The spiral region has a comfortable middle.
- **Which assumption fails first in the real world?** Constant speed in the short-period model fails if the mode is slow (sailplanes); linearity fails first for large gusts (α beyond the linear lift range).
- **What measurement could falsify the model?** A flight-test pulse response whose α(t) shows two frequencies (coupling with the phugoid or a structural mode) or a decay rate that depends on the pulse size (nonlinearity).
- **Does the result survive a change of solver?** Hand eigenvalues, `eig`, `expm`, RK4 and `solve_ivp` agree to 10⁻⁸. Does the *classification* survive a 20 % change in M_α? At M_α = −4, yes; at M_α = +0.9, no. Say which design you would accept without an augmentation system.
- **Which decision changes if the uncertainty doubles?** With derivatives known to 40 %, the cruise short period is still a well-damped spiral (ζ between 0.37 and 0.56); but a design at M_α = −0.5 (a light node/spiral boundary) could be anywhere from ζ = 0.75 to a saddle, and would be rejected without augmentation.

## Check yourself

1. Write the characteristic polynomial of a 2×2 matrix in terms of its trace and determinant, and state what the sign of Δ tells you.
2. Classify A = [[1, 2], [−3, −4]] from (τ, Δ) without computing eigenvalues.
3. The eigenvalues of a system are −0.5 ± 3i s⁻¹. Give ω_n, ζ, the period and the time to half amplitude.
4. Why does a repeated eigenvalue produce a factor t in the solution, and what is the solution formula?
5. Which eigenvector do the trajectories of a stable node hug as they approach the origin, and why?
6. What quantity is conserved along the trajectories of a center, and what is it for the undamped oscillator?
7. For the short-period model with Z_α/V₀ = −0.8 and M_q = −1.2, at what M_α does the aircraft become a saddle, and what does that mean physically?
8. Lanchester's phugoid has ω_n = √2 g/V₀. What is the period at 100 m/s, and on what does the damping ratio depend?
