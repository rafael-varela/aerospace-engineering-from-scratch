# P02.06  Nonlinear ODEs, Equilibria & Linearization

> Near a point where nothing moves, every smooth system is linear. Find those points, linearise there, and let the Jacobian's eigenvalues tell you which way the world falls.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the equilibria and the Jacobian of a two-state nonlinear system, classify each equilibrium from the trace and determinant, and state when that verdict carries over to the nonlinear system and when it does not. 2. Implement and test, from scratch, a fourth-order Runge–Kutta integrator for the damped nonlinear pendulum, verified against the closed-form linear response, against `scipy.integrate.solve_ivp`, and against energy conservation. 3. Derive the exact large-amplitude period of a pendulum as an elliptic integral and its series T₀(1 + θ₀²/16 + …), and measure the period from a simulated trajectory. 4. Decide how large a release amplitude an aircraft inertia swing test may use before the small-angle formula corrupts the inferred moment of inertia by more than a stated tolerance. |
| **Prerequisites** | `P02.01` (autonomous first-order equations, the phase line), `P02.02` (damped oscillator, damping ratio), `P02.04` (phase portraits, trace–determinant classification), `P01.05` (Taylor series and the remainder), `P01.08` (partial derivatives, the Jacobian), `P01.16` (eigenvalues of a 2×2 matrix), `P00.09` (SymPy). Skill check: you can write the Taylor series of sin θ to the cubic term, find the eigenvalues of a 2×2 matrix from its trace and determinant, and write a second-order ODE as a first-order system. |
| **Engineering question** | A 4.2 kg fixed-wing UAV is hung from a knife edge 0.35 m above its centre of mass and swung to measure its pitch moment of inertia; how large a release amplitude may the technician use before the small-angle formula I = m g d T²/(4π²) − m d² is wrong by more than 1 %? |
| **Assumptions and limits** | Rigid body on a frictionless pivot with linear viscous damping of ratio ζ (the test uses ζ = 0.01); no aerodynamic damping proportional to speed squared, no motion of the pivot, no swing out of the plane; standard gravity. The linearisation describes the motion only near an equilibrium; the exact period formula holds for the undamped pendulum released from rest and for amplitudes below 180°. The classification by the Jacobian is reliable for a saddle, node or spiral and *silent* for a center: energy conservation, not linear algebra, is what makes the undamped pendulum's center real. |
| **Analytical work** | The compound pendulum reduced to a simple one (effective length); the first-order system and its two equilibria; Taylor expansion of the vector field, the Jacobian, eigenvalues at both equilibria; the trace–determinant classification and the Hartman–Grobman idea; the energy integral, the separatrix, the exact period as a complete elliptic integral and its series; the error of the small-angle inertia formula and its amplification by the lever ratio I_pivot/I_cg; the Van der Pol oscillator as a system no linearisation can describe. |
| **Computational deliverable** | `starter/main.py`: `effective_length`, `natural_frequency`, `inertia_from_period`, `pendulum_rhs`, `linear_pendulum_rhs`, `pendulum_equilibria`, `jacobian_numeric`, `pendulum_jacobian`, `pendulum_jacobian_sympy`, `classify_equilibrium`, `pendulum_eigenvalues`, `rk4_step`, `integrate_rk4`, `linear_response`, `energy`, `separatrix_energy`, `is_rotating`, `period_small_angle`, `period_exact`, `period_series`, `period_from_trajectory`, `max_amplitude_for_period_error`, `max_amplitude_for_inertia_error`, `van_der_pol_rhs`, `limit_cycle_amplitude`, `swing_test_error`. Tests prove: the analytic Jacobian matches central differences and SymPy; the classification follows the trace–determinant plane; RK4 converges at fourth order against the closed form and agrees with `solve_ivp` to 10⁻⁷ rad; the undamped energy drifts by less than 10⁻⁹ and the damped energy never rises; the exact period reduces to T₀, matches the series to O(θ₀⁶), is 18.034 % longer at 90°, and is reproduced by zero crossings of a simulated trajectory; the separatrix separates swinging from rotating; the inertia error grows as θ₀² with the predicted lever factor; the Van der Pol amplitude is the same from inside and outside the cycle. |
| **Practical activity** | **computational experiment**: a sweep of the swing-test error against release amplitude, a convergence study of RK4 against the closed form, and a test of how far the linear verdict for each equilibrium survives finite perturbations. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pendulum_linearization.png`, the hand derivation of the Jacobian and the period series, the amplitude-error table, and a one-paragraph recommendation for the swing-test procedure with its numbers. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the Jacobian and both classifications are derived by hand and agree with the code; all public tests pass; the period series is derived to the θ₀² term; the swing-test recommendation gives an amplitude limit with the tolerance it comes from. |
| **References** | `@mit-ocw-18-03-differential-equations` L31 "Non-linear Autonomous Systems", L32 "Limit Cycles", notes 37 "Nonlinear systems", 38 "Linearization near equilibria; the nonlinear pendulum" (required-free); L33 "Relation Between Non-linear Systems and First-order ODEs", notes 39 "Limitations of the linear: limit cycles and chaos" (supplementary-free). |
| **Extensions** | ★ Add damping proportional to θ̇|θ̇| (aerodynamic) and show that the log decrement is no longer constant. ★ Linearise the gravity-gradient pitch libration of a satellite, θ̈ + 3n²(I_x − I_z)/I_y sin θ cos θ = 0, about all its equilibria and find the libration period of a 90-minute orbit. ★ Drive the damped pendulum sinusoidally at the pivot and find the forcing amplitude above which the motion stops repeating (the route to chaos in notes 39). |
| **Time** | unmeasured |

## 1. The engineering problem

Before a new aircraft flies, its moments of inertia are measured, because every prediction of its short-period and roll response depends on them and no drawing-office estimate is trusted to better than 10 %. The classic method is the **swing test**: hang the airframe from a knife edge, displace it, release it, and time the swings. The airframe is a **compound pendulum**, a rigid body of mass m with its centre of mass a distance d below the pivot and moment of inertia I_pivot about the pivot. For small swings its period is T₀ = 2π√(I_pivot/(m g d)), so I_pivot = m g d T₀²/(4π²) and the parallel-axis theorem (P01.10) gives the quantity that matters, I_cg = I_pivot − m d².

The formula has a hidden assumption: that sin θ may be replaced by θ. P01.05 showed that the period is then within 1 % of the truth only for amplitudes below 22.8°. But the quantity reported is not the period; it is I_cg, obtained by squaring the period and subtracting a large number. For our UAV, I_pivot = 0.7245 kg m² and m d² = 0.5145 kg m², so I_cg = 0.21 kg m² is the small difference of two large numbers, and a 1 % error in the period becomes a 2 % error in I_pivot and a 6.9 % error in I_cg. The question "how large may the amplitude be?" therefore has a sharper answer than P01.05 gave, and it needs the full nonlinear equation.

That equation, θ̈ + 2ζω₀θ̇ + ω₀² sin θ = 0, cannot be solved by the methods of P02.02 because sin θ is not linear. This lesson is about what *can* be done with a nonlinear system: find the points where it stands still, linearise there, read the local behaviour from the eigenvalues of a matrix, and integrate numerically for everything else. The pendulum is the standard first example because every step can be checked against something exact (energy, an elliptic integral) and because the same picture, a stable equilibrium, an unstable one and a separatrix between swinging and tumbling, is the picture of an aircraft's pitch dynamics, a satellite's libration and a rocket's attitude.

Success criteria, written down first: the hand Jacobian agrees with central differences to 10⁻⁶ and with SymPy to 10⁻¹²; RK4 shows fourth-order convergence (error ratio 16 per halving) and agrees with `solve_ivp` to 10⁻⁷ rad; the undamped energy is conserved to 10⁻⁹; the period measured from a simulated 60° swing agrees with the elliptic integral to 10⁻⁵; and the amplitude limit for a 1 % inertia error is reported with the lever factor that explains why it is smaller than the 22.8° of P01.05.

## 2. Understand and derive

### 2.1 From the compound pendulum to the simple one

```
          ▲ pivot (knife edge)
          │╲
          │ ╲  d = 0.35 m
          │  ╲
          │ θ ╲
          │    ● centre of mass (m = 4.2 kg, I_cg = 0.21 kg m²)
          ▼
         m g
```

**Assumption 1 (rigid body, planar swing, frictionless pivot).** The moment of gravity about the pivot is −m g d sin θ, so Newton's second law for rotation gives I_pivot θ̈ = −m g d sin θ − c θ̇, with c a viscous damping coefficient [N m s]. Divide by I_pivot:

**θ̈ + 2ζω₀ θ̇ + ω₀² sin θ = 0**,   ω₀² = m g d / I_pivot = g / L_eff,   L_eff = I_pivot/(m d).   (1)

L_eff is the length of the simple pendulum with the same small-angle period: for the UAV, L_eff = 0.7245/(4.2 × 0.35) = 0.4929 m, ω₀ = √(9.80665/0.4929) = 4.4607 rad/s, T₀ = 2π/ω₀ = 1.4086 s. Units: m g d is N m = kg m² s⁻², over kg m² gives s⁻². Good. The damping ratio ζ = c/(2 I_pivot ω₀) is defined as in P02.02 so that the linearised equation is the familiar damped oscillator.

### 2.2 The first-order system and its equilibria

Write the state as **y** = (θ, θ̇). Then (1) is the **autonomous system** (no explicit t on the right)

**ẏ** = **f**(**y**) = (θ̇, −2ζω₀ θ̇ − ω₀² sin θ).   (2)

An **equilibrium** is a state where **f** = **0**: θ̇ = 0 and sin θ = 0, so θ = 0, ±π, ±2π, … Physically distinct are only two: the pendulum hanging down, (0, 0), and balanced upside down, (π, 0). Every other equilibrium is one of these two seen after a whole number of turns. An equilibrium is where the system can stay forever; whether it *will* stay there when nudged is the stability question.

### 2.3 Linearisation: the Jacobian

Near an equilibrium **y*** write **y** = **y*** + **u** with **u** small and expand each component of **f** in a Taylor series (P01.08) to first order:

**f**(**y*** + **u**) = **f**(**y***) + J **u** + O(|**u**|²),   J_ij = ∂f_i/∂y_j evaluated at **y***.   (3)

The first term is zero by definition of an equilibrium, so **u̇** ≈ J **u**: a *linear* system with constant coefficients, exactly what P02.04 and P02.05 solve. The matrix J is the **Jacobian**. For (2),

$$
J(\theta) = \begin{pmatrix} 0 & 1 \\ -\omega_0^2\cos\theta & -2\zeta\omega_0 \end{pmatrix}.
$$

At the hanging equilibrium cos 0 = 1 and J is the matrix of the damped linear oscillator θ̈ + 2ζω₀θ̇ + ω₀²θ = 0, which is what "replace sin θ by θ" means. At the inverted equilibrium cos π = −1 and the sign of the stiffness flips: **u̇** describes θ̈ + 2ζω₀θ̇ − ω₀²θ = 0, a negative spring. For the UAV with ζ = 0.01, J(0) = [[0, 1], [−19.898, −0.0892]] and J(π) = [[0, 1], [+19.898, −0.0892]]; the code checks both against a central-difference Jacobian (agreement 3 × 10⁻¹²) and against SymPy's `Matrix.jacobian`.

### 2.4 Classifying the equilibria

The eigenvalues of a 2×2 matrix are λ = ½[tr J ± √(tr² J − 4 det J)], so the trace and determinant decide everything (P02.04; `@mit-ocw-18-03-differential-equations` notes 38):

| det J | tr J | tr² − 4 det | eigenvalues | type |
|---|---|---|---|---|
| < 0 | any | > 0 | real, opposite signs | **saddle**: unstable |
| > 0 | < 0 | < 0 | complex, negative real part | **stable spiral** |
| > 0 | < 0 | > 0 | real, both negative | **stable node** |
| > 0 | = 0 | < 0 | pure imaginary | **center**: neutrally stable (linear) |
| > 0 | > 0 | | positive real part | **unstable spiral / node** |

Hanging: det J = ω₀² cos 0 = +19.898, tr J = −2ζω₀ = −0.0892, so a **stable spiral** with eigenvalues −ζω₀ ± iω₀√(1 − ζ²) = −0.0446 ± 4.4604i s⁻¹. Without damping tr J = 0 and the linear verdict is a **center**. Inverted: det J(θ) = ω₀² cos θ in general, so at θ = π it is −ω₀² = −19.898 < 0 for every ζ: a **saddle**, eigenvalues +4.416 and −4.506 s⁻¹ with ζ = 0.01, ±ω₀ = ±4.4607 s⁻¹ without damping. A small tilt from the inverted position grows as e^{4.42 t}: it doubles every ln 2/4.42 = 0.157 s. That number is as real as the period: it is how fast a balanced rocket falls over without control.

### 2.5 When does the linear verdict hold for the nonlinear system?

The **Hartman–Grobman theorem** (stated, not proved, in notes 38) says: if no eigenvalue of J has zero real part, the nonlinear phase portrait near the equilibrium is a continuous distortion of the linear one. Saddles, nodes and spirals are **hyperbolic** and their verdicts transfer. The **center** is the exception: its eigenvalues are pure imaginary, and the neglected O(|**u**|²) terms can turn closed linear orbits into a slow spiral either way. Example: ẋ = −y − x(x² + y²), ẏ = x − y(x² + y²) has J = [[0, −1], [1, 0]] at the origin, a center, but r² = x² + y² obeys d(r²)/dt = −2r⁴ < 0, so every orbit spirals in. The linearisation cannot tell you this.

For the undamped pendulum the center *is* real, but the reason is not the Jacobian; it is the **energy integral** of 2.6. The lesson's `classify_equilibrium` returns "center" when the trace is zero and the lesson's rule is: a center from a Jacobian is a question, not an answer. Find a conserved quantity or integrate.

### 2.6 Energy, the separatrix, and the exact period

Multiply (1) with ζ = 0 by θ̇ and integrate: d/dt[½θ̇² + ω₀²(1 − cos θ)] = 0. The quantity

**E = ½ θ̇² + ω₀² (1 − cos θ)**   [rad²/s²]   (4)

(the mechanical energy per unit I_pivot) is conserved, so each orbit in the phase plane is a level curve of E. Near (0, 0), E ≈ ½θ̇² + ½ω₀²θ² gives ellipses: the center. The curve through the saddle (π, 0) has E = 2ω₀², and it is the **separatrix**: inside it the pendulum swings back and forth, outside it (E > 2ω₀²) it goes over the top and rotates. `is_rotating` is this inequality and the test confirms it by integration: released at θ = 0 with θ̇ = 0.98 × 2ω₀ the pendulum turns back, at 1.02 × 2ω₀ it passes θ = π.

For a swing released from rest at amplitude θ₀, E = ω₀²(1 − cos θ₀) and (4) gives θ̇ = ω₀√(2(cos θ − cos θ₀)). The quarter period is the time from θ = 0 to θ = θ₀:

$$
\frac{T}{4} = \frac{1}{\omega_0}\int_0^{\theta_0}\frac{d\theta}{\sqrt{2(\cos\theta - \cos\theta_0)}} .
$$

With cos θ = 1 − 2 sin²(θ/2) and the substitution sin(θ/2) = k sin φ, k = sin(θ₀/2), the integral becomes

$$
T = \frac{4}{\omega_0}\int_0^{\pi/2}\frac{d\phi}{\sqrt{1 - k^2\sin^2\phi}} = \frac{4}{\omega_0}\,K(k),
\qquad K(k) = \frac{\pi}{2}\left(1 + \frac{k^2}{4} + \frac{9k^4}{64} + \cdots\right),   (5)
$$

K being the **complete elliptic integral of the first kind** (`scipy.special.ellipk`, which takes m = k²). Expanding k² = sin²(θ₀/2) = θ₀²/4 − θ₀⁴/48 + … and collecting:

**T = T₀ (1 + θ₀²/16 + 11θ₀⁴/3072 + …)**.   (6)

The period grows with amplitude because the restoring torque ω₀² sin θ is weaker than ω₀²θ: the pendulum spends longer near the turning points. At 90°, T/T₀ = 2K(1/√2)/π = 1.18034; at 179°, 3.90; as θ₀ → π the period diverges logarithmically, because the orbit creeps past the saddle. The code evaluates (5) and (6) and the test checks that they agree to O(θ₀⁶) and that (5) reduces to T₀ as θ₀ → 0.

### 2.7 The swing-test error

A technician who releases at θ₀ measures T = T₀(1 + θ₀²/16 + …) and inserts it in the small-angle formula. The inferred inertia is

I_cg,inferred = I_pivot (T/T₀)² − m d²,   so   **(I_inferred − I_cg)/I_cg = (I_pivot/I_cg)[(T/T₀)² − 1] ≈ (I_pivot/I_cg) θ₀²/8**.   (7)

The factor I_pivot/I_cg is the **lever ratio**; for the UAV it is 0.7245/0.21 = 3.45, so the inertia error is about 6.9 times the period error. The code's table, exact period in the small-angle formula:

| release amplitude | T/T₀ | inferred I_cg [kg m²] | error |
|---|---|---|---|
| 5° | 1.00048 | 0.21069 | +0.33 % |
| 10° | 1.00191 | 0.21277 | +1.32 % |
| 20° | 1.00767 | 0.22116 | +5.31 % |
| 30° | 1.01741 | 0.23544 | +12.1 % |
| 45° | 1.03997 | 0.26908 | +28.1 % |

Equation (7) at 20°: 3.45 × 0.349²/8 = 5.25 %, against the exact 5.31 %; the difference is the θ₀⁴ term. Bisection on the exact expression gives the answer to the engineering question: **release below 8.7° for a 1 % inertia error**, 6.2° for 0.5 %, 12.3° for 2 %. P01.05's 22.8° is the right number for the period and the wrong number for the inertia, and the reason is the subtraction in I_cg = I_pivot − m d²: hang the airframe closer to its centre of mass (smaller d) and the lever ratio falls, but so does the restoring torque and the period gets long and hard to time. Real swing tests use d of about half a metre and amplitudes of 5–10°, and now you know why.

### 2.8 A limit cycle, for contrast

Not every closed orbit is a center's. The **Van der Pol oscillator** ẍ − μ(1 − x²)ẋ + x = 0 has one equilibrium, the origin, with J = [[0, 1], [−1, μ]]: tr J = μ > 0, det J = 1, an unstable spiral for 0 < μ < 2. Every trajectory leaves the origin, yet none runs away: for |x| > 1 the damping term changes sign and removes energy. The result is a **limit cycle**, a single isolated closed orbit that attracts every trajectory from inside and outside (`@mit-ocw-18-03-differential-equations` L32). Its amplitude is 2 + O(μ²): the code finds 2.0086 for μ = 1 starting from (0.1, 0) and the same 2.0086 from (4, 0). No linearisation anywhere predicts that number; it is a property of the whole phase plane. Limit cycles are the mathematics of wing flutter at its stable amplitude, of a badly tuned autopilot hunting, and of every self-excited oscillation.

### 2.9 Where the model stops being valid

- **Linear damping.** Air drag on an airframe swinging at 0.2 m/s is proportional to speed squared; the log decrement then falls with amplitude instead of staying constant (extension). The period is affected only at O(ζ²), which is why the swing test ignores damping for the inertia and uses it only to check the amplitude is decaying slowly.
- **Rigid pivot.** A knife edge that rolls or a cable that stretches adds a degree of freedom and a second period.
- **Planar swing.** Any out-of-plane component makes the motion a spherical pendulum, whose period at the same amplitude is slightly different.
- **Linearisation.** Valid within the basin of the equilibrium only. For the inverted pendulum the basin is empty; for the hanging one it is the whole region inside the separatrix, but the *quantitative* linear solution is good only while θ₀²/16 is negligible.

## 3. Build the computational model

Open `starter/main.py`. `uav_swing_test()` holds the nominal case; every function takes its parameters explicitly so the tests can vary them.

### Step 1: the pendulum and its linearisation

```python
def pendulum_rhs(t, y, omega0, zeta=0.0):
    theta, theta_dot = y
    return np.array([theta_dot, -2.0 * zeta * omega0 * theta_dot - omega0**2 * np.sin(theta)])   # (2)

def pendulum_jacobian(theta, omega0, zeta=0.0):
    return np.array([[0.0, 1.0], [-omega0**2 * np.cos(theta), -2.0 * zeta * omega0]])            # (3)

def jacobian_numeric(f, x, eps=1e-6):
    x = np.asarray(x, float); n = x.size
    f0 = np.asarray(f(x), float); J = np.zeros((f0.size, n))
    for j in range(n):
        dx = np.zeros(n); dx[j] = eps
        J[:, j] = (f(x + dx) - f(x - dx)) / (2.0 * eps)
    return J
```

The right-hand sides use the SciPy signature `f(t, y, ...)` so that the same function feeds `solve_ivp`. `jacobian_numeric` takes a function of the state alone; wrap the right-hand side in a lambda. With eps = 10⁻⁶ the truncation error is O(10⁻¹²) and the rounding error O(10⁻¹⁶/10⁻⁶) = 10⁻¹⁰: that is why the test asks for 10⁻⁶ agreement, not 10⁻¹². `pendulum_jacobian_sympy` builds the same matrix with `sympy.Matrix.jacobian` and evaluates it; use it once to check your hand derivative, then keep the hand version.

### Step 2: classification

```python
def classify_equilibrium(J, tol=1e-12):
    tr = J[0, 0] + J[1, 1]; det = J[0, 0] * J[1, 1] - J[0, 1] * J[1, 0]
    scale = max(1.0, np.abs(J).max() ** 2)
    if abs(det) <= tol * scale:   return "degenerate"
    if det < 0.0:                 return "saddle"
    disc = tr * tr - 4.0 * det
    if abs(tr) <= tol * np.sqrt(scale): return "center"
    if abs(disc) <= tol * scale:  return "degenerate"
    if disc < 0.0:                return "stable spiral" if tr < 0.0 else "unstable spiral"
    return "stable node" if tr < 0.0 else "unstable node"
```

Order matters: a negative determinant is a saddle whatever the trace; a zero trace with positive determinant is a center; the discriminant decides spiral against node. The tolerances are relative to the size of J so that the same code works for ω₀ = 4 and for the satellite libration at ω₀ = 0.0016.

### Step 3: RK4 from scratch

```python
def rk4_step(f, t, y, h):
    k1 = f(t, y)
    k2 = f(t + 0.5 * h, y + 0.5 * h * k1)
    k3 = f(t + 0.5 * h, y + 0.5 * h * k2)
    k4 = f(t + h, y + h * k3)
    return y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
```

`integrate_rk4(f, y0, t_end, h)` takes ceil((t_end − t₀)/h) steps and shortens the last one to land exactly on t_end, returning `(t, Y)` with one row per step. RK4 samples the slope four times per step and combines them with Simpson's weights 1, 2, 2, 1; the local error is O(h⁵) and the global error O(h⁴), so halving h should divide the error by 16. On the linearised UAV pendulum over three periods the code finds errors 2.69 × 10⁻⁶, 1.68 × 10⁻⁷, 1.05 × 10⁻⁸, 6.54 × 10⁻¹⁰ rad for h = 20, 10, 5, 2.5 ms: ratios 16.0, 16.0, 16.0.

### Step 4: energy and period

`energy` is (4); `separatrix_energy` is 2ω₀²; `period_exact` is 4K(k)/ω₀ with `scipy.special.ellipk(k**2)`; `period_series` is (6); `period_from_trajectory` finds the upward zero crossings of θ(t), locates each by linear interpolation between the samples that bracket it, and returns the mean spacing. Linear interpolation on a 1 ms grid gives the crossing time to about 10⁻⁷ s, which is why the test asks for 10⁻⁵ relative agreement with (5).

### Step 5: the swing-test study

`swing_test_error(theta0, case)` returns T₀, the exact period, the inertia the small-angle formula infers from it, and the relative error (7). `max_amplitude_for_period_error(tol)` and `max_amplitude_for_inertia_error(tol, case)` bisect the monotone error curves; `main()` prints the table of 2.7.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-06-nonlinear-equilibria-linearization`. The tests:

1. **Analytic cases.** The right-hand side vanishes at both equilibria and differs from the linear model by ω₀²θ³/6; the analytic Jacobian matches central differences to 10⁻⁶ and SymPy to 10⁻¹²; the eigenvalues at (0, 0) are −ζω₀ ± iω₀√(1 − ζ²) and at (π, 0) are ±ω₀ for ζ = 0.
2. **Classification.** Center, stable spiral, stable node at the hanging equilibrium for ζ = 0, 0.05, 1.5; saddle at the inverted one for every ζ; unstable spiral and node for negative damping; degenerate for a singular J.
3. **Convergence and agreement.** RK4 error ratios above 2^3.7 per halving against `linear_response`; agreement with `solve_ivp` at rtol = 10⁻¹⁰ to 10⁻⁷ rad over four periods at 60°.
4. **Conservation.** Undamped energy drift below 10⁻⁹ over five periods at 120°; with damping the energy at successive turning points never increases.
5. **Limiting cases.** `period_exact` → T₀ as θ₀ → 0; the series agrees to O(θ₀⁶); T/T₀ = 1.18034 at 90°; period ∝ √L_eff; measured and exact period agree to 10⁻⁵ at 60°.
6. **The decision.** The inertia error at 10° is four times that at 5° and (2 I_pivot/I_cg) times the period error; 22.8° for a 1 % period error reproduces P01.05; the inertia limit is smaller.
7. **Limit cycle.** The Van der Pol origin is an unstable spiral and the amplitude from (0.1, 0) equals that from (4, 0) to 10⁻³, within 1 % of 2.

### Compare

`solve_ivp` with its default RK45 (Dormand–Prince) is the oracle. It is also Runge–Kutta, but it takes two estimates of different order at each step, uses their difference to estimate the local error and chooses the next step size to keep that error at the requested tolerance. Our fixed-step RK4 at h = 2 ms and `solve_ivp` at rtol = 10⁻¹⁰ agree to 5.8 × 10⁻¹⁰ rad over five periods at 20°, and both conserve energy to 2.5 × 10⁻¹¹. The library wins on cost: it takes a few hundred steps where we take 3500, because it enlarges the step where θ is slowly varying. We wrote our own to know what a step *is*, to see the order of convergence with our own eyes, and because in Phase 03 the same four lines become the core of every simulator you will write; in practice you call `solve_ivp` and check its tolerance by halving it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three studies with `main.py` as the instrument.

1. **Amplitude sweep.** For release amplitudes 5°, 10°, 15°, 20°, 30°, 45° compute T/T₀ and the inertia error. **Prediction to write down first:** the error should go as θ₀² with coefficient I_pivot/(8 I_cg) = 0.431 rad⁻²; check the ratio error(20°)/error(10°) = 4.
2. **Lever sweep.** Keep I_cg and m and move the pivot: d = 0.20, 0.35, 0.50, 0.80 m. For each, the lever ratio, T₀, and the amplitude limit for a 1 % inertia error. Predict first whether a longer d helps.
3. **Linear verdict against finite perturbation.** Start the damped (ζ = 0.01) pendulum at 1°, 10°, 60°, 170° and at 179.9° and 180.1° with θ̇ = 0. Record whether it settles to (0, 0), to (2π, 0) or stays near (π, 0), and the time to fall below 1° of the equilibrium it reaches. The linear verdict (stable spiral at 0, saddle at π) says which way each goes but not how long it takes.

Raw data format: three CSV tables, `amp_deg, T_ratio, I_error`; `d, lever, T0, amp_limit_deg`; `theta0_deg, final_equilibrium, settling_time_s`, each with a one-line header naming the fixed parameters. All data in this part are synthetic (generated by your own integrator) and must be labelled as such in your report.

What you should find: the quadratic law holds to 2 % up to 20° and the coefficient matches; the amplitude limit *falls* as the pivot moves away from the centre of mass, because the lever ratio I_pivot/I_cg = 1 + m d²/I_cg rises (1.8, 3.45, 6.0, 13.8 for the four pivots, with limits 12.1°, 8.7°, 6.6°, 4.4°), so a long suspension makes the period easier to time and the amplitude harder to tolerate; the 179.9° and 180.1° starts go to (0, 0) and (2π, 0) respectively and reach 1° of them at the same time, 109.8 s, after lingering near the saddle for about 1 s (six doubling times).

## 5. Validate and improve

| quantity (UAV, ζ = 0 unless stated) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| L_eff, ω₀, T₀ | 0.4929 m, 4.4607 s⁻¹, 1.4086 s | 0.49286, 4.46067, 1.40858 | 0 | exact formulas |
| J(0), ζ = 0.01 | [[0, 1], [−19.898, −0.0892]] | same; numeric to 3 × 10⁻¹² | rounding | central-difference step |
| eigenvalues at (0, 0), ζ = 0.01 | −0.0446 ± 4.4604i | −0.04461 ± 4.46044i | 0 | — |
| doubling time at (π, 0) | ln 2/ω₀ = 0.155 s | 0.1554 s (0.1570 with ζ = 0.01) | 1 % | damping slows the fall |
| T/T₀ at 20° | 1 + 0.349²/16 = 1.00761 | 1.00767 (elliptic) | 6 × 10⁻⁵ | θ₀⁴ term, 11θ₀⁴/3072 = 5.3 × 10⁻⁵ |
| inertia error at 20° | 3.45 × θ₀²/8 = 5.25 % | 5.31 % | 0.06 % | same term, doubled by the square |
| period from trajectory, 20° | 1.41938 s (elliptic) | 1.41938 s (zero crossings) | < 10⁻⁵ | interpolation on a 2 ms grid |
| Van der Pol amplitude, μ = 1 | 2 (leading order) | 2.0086 | 0.4 % | O(μ²) correction |

Every discrepancy has a named cause and each is far below what a real swing test can resolve: timing ten swings with a stopwatch gives T to about 0.1 %, which is 0.2 % in I_pivot and 0.7 % in I_cg through the lever, and d is rarely known to better than 2 mm, another 0.8 %. The amplitude error at 10° (1.3 %) is therefore comparable to the other errors and at 20° it dominates them. The recommendation: release below 9°, time ten swings, and correct the period with (6) anyway, since the amplitude is easy to read off a protractor and the correction costs nothing.

What the agreement does *not* show: that the real airframe swings as a rigid body on a frictionless knife edge with linear damping. Measuring the log decrement at several amplitudes is the check: a constant decrement means linear damping; one that falls with amplitude means quadratic (aerodynamic) damping, and then the period correction for damping, though still O(ζ²), should be rechecked.

Improvements worth making: fit the amplitude correction to the measured period by two releases at different amplitudes and solve for T₀ (problem 4); replace the fixed-step RK4 with step doubling and an error estimate; add the pivot's own compliance as a second degree of freedom. Improvement worth not making: trusting the linear verdict at an equilibrium whose Jacobian has a zero eigenvalue.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of J at both equilibria with the classification; the derivation of (6) to the θ₀² term; `pendulum_linearization.png`; the three CSV tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that gives the recommended release amplitude, the tolerance it comes from, the lever ratio, and the assumption you trust least.

### Adversarial self-checks

- **What if the mass doubles?** With the same geometry and I_cg ∝ m (same shape, denser), L_eff, ω₀, T₀ and every error are unchanged: the period of a pendulum does not depend on its mass. If instead only m doubles and I_cg stays, L_eff falls, the lever ratio rises and the amplitude limit tightens. Check both with the code.
- **What if the amplitude goes to zero? To 180°?** T → T₀ and the linear model becomes exact; T → ∞ logarithmically and the swing lingers at the saddle. No release amplitude gives a *shorter* period than T₀: the small-angle inertia is always an overestimate.
- **Which assumption fails first in the real world?** Linear damping: a wing swinging through air is quadratically damped. Second: the rigid knife edge.
- **What measurement could falsify the model?** The period at two amplitudes: (6) predicts T(20°)/T(10°) = 1.00575. A measured ratio of 1.000 means the amplitudes were not what you thought; one of 1.02 means something other than gravity is restoring the airframe (a cable, a stop).
- **Does the result survive a change in step or solver?** RK4 at h = 2 ms and `solve_ivp` at rtol = 10⁻¹⁰ agree to 10⁻⁹ rad; halving h changes the measured period by less than 10⁻⁸ s. The classification does not depend on the solver at all.
- **Which decision changes if the uncertainty doubles?** If the period were known only to 0.2 %, the inertia uncertainty from timing alone would be 1.4 %, and the amplitude limit for "1 % total" would have to tighten from 8.7° to about 5°, where the swing is hard to see; the better decision is then to apply the correction (6) rather than to restrict the amplitude.

## Check yourself

1. Define equilibrium and Jacobian for the system **ẏ** = **f**(**y**) in one sentence each, and write the Jacobian of the pendulum at a general θ.
2. From the trace and determinant alone, classify J = [[0, 1], [−9, −0.6]] and J = [[0, 1], [9, −0.6]].
3. Why is the linear verdict "center" not a verdict about the nonlinear system, and what replaces it for the undamped pendulum?
4. What is the energy of the separatrix in units of ω₀², and what does a state above it do?
5. The period at 90° is 18 % longer than T₀. Give the first-order series estimate at 90° and say why it is poor there.
6. A swing test releases at 15°. Using I_pivot/I_cg = 3.45, estimate the error in I_cg from (7).
7. RK4 at step h gives an error of 10⁻⁶; what error do you expect at h/2 and at 2h?
8. Name one physical system whose steady oscillation is a limit cycle rather than a center, and say what a linearisation at its equilibrium shows.
