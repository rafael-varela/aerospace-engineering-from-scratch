# P06.11  Lagrangian Mechanics

> Write two energies in the coordinates you like, differentiate, and the equations of motion fall out with every constraint force already gone.

| | |
|---|---|
| **Learning outcomes** | 1. Choose generalized coordinates for a constrained system, count its degrees of freedom, and write its kinetic and potential energies in those coordinates. 2. Derive the equations of motion of a double pendulum from Lagrange's equations by hand and in SymPy, put them in the form M(q) q̈ = f(q, q̇), and show that the energy function is conserved when the Lagrangian does not depend on time. 3. Derive the same accelerations independently by Newton–Euler (free-body diagrams with the link tensions as unknowns) and show numerically that the two agree to round-off at random states. 4. Linearise about the hanging equilibrium, solve K φ = ω² M φ for the normal modes in closed form, and confirm by simulation that a small motion started in a mode stays in it at its frequency. 5. Measure the sensitivity to initial conditions at large amplitude with a finite-time Lyapunov exponent, and decide over what horizon a simulation of a chaotic motion can be trusted. |
| **Prerequisites** | `P02.13` (calculus of variations, the Euler–Lagrange equation, Hamilton's principle), `P02.06` (linearisation about an equilibrium), `P01.16` (eigenvalues of symmetric matrices), `P02.02` (second-order linear ODEs), `P03.06` (adaptive integrators and `solve_ivp`), `P04.04` (rotation, angular momentum), `P06.01` (polar and path coordinates), `P06.02` (particle kinetics, free-body diagrams), `P00.09` (SymPy). Skill check: you can differentiate l sin θ(t) twice with respect to t, and solve a 2 × 2 linear system by Cramer's rule. |
| **Engineering question** | A helicopter hovers with an 80 kg hook block on a 12 m cable and a 1500 kg container slung 4 m below the hook (illustrative numbers). What are the two swing periods the pilot must not excite, what do the cable loads do in a 10° swing, at what amplitude does the linear model stop being good enough, and how long can a simulation of a violently swinging double pendulum be believed? The answers: 7.99 s and 0.79 s; the cable tension swings between 0.985 and 1.031 times the static weight; the period error passes 1 % near 22°; and, at a Lyapunov exponent near 1.4 s⁻¹, each hundredfold improvement in the initial state buys only about 2.7 s of trustworthy prediction. |
| **Assumptions and limits** | Planar motion; point masses on massless, inextensible links (the cable stays taut: tension > 0); a fixed suspension point (a perfectly station-keeping helicopter); no aerodynamic forces (no rotor downwash, no drag on the load); no friction at the joints. The chaos study uses the textbook equal-mass, equal-length pendulum, not the slung load. A real slung load also swings sideways, yaws, bounces on the cable's elasticity and is pushed by downwash; all are left out. |
| **Analytical work** | Generalized coordinates and degrees of freedom; Lagrange's equations from Hamilton's principle and why constraint forces do not appear; T and V of the double pendulum; M(q), the velocity-squared terms and gravity; the energy function and its conservation; the Newton–Euler free-body diagrams with tensions; linearisation, the frequency equation m₁l₁l₂ω⁴ − (m₁ + m₂)g(l₁ + l₂)ω² + (m₁ + m₂)g² = 0, its limits and mode shapes; the hook-mode estimate; the prediction horizon from a Lyapunov exponent. |
| **Computational deliverable** | `starter/main.py`: `euler_lagrange`, `double_pendulum_lagrangian`, `derive_equations`, `derive_accelerations`, `mass_matrix`, `generalized_forces`, `rhs`, `newton_euler`, `simulate`, `energy`, `linearized_matrices`, `normal_modes`, `normal_mode_frequencies_closed_form`, `lyapunov_estimate` (`DoublePendulum`, `slung_load`, `lab_pendulum`, `lagrange_accelerations`, `positions`, `dominant_frequency` and `main` are given). Tests prove: Lagrange's equations of a spring–mass and a simple pendulum; the SymPy M(q) and f(q, q̇) equal the hand forms; m₂ → 0 leaves a simple pendulum; M(q) is symmetric positive definite; Lagrange equals Newton–Euler at 150 random states for three parameter sets; tensions in the static and rigid-rotation limits; energy conserved to 10⁻⁸ and equal to an independent Cartesian formula; time reversibility; the modes solve the frequency equation, are M-orthogonal and tend to two independent pendulums; the SymPy Jacobian at equilibrium is −M⁻¹K; small motions oscillate in their mode at its frequency and a 20° swing is slower by the amplitude correction; a 120° release is chaotic and a 5° release is not. |
| **Practical activity** | **computational experiment**: an amplitude sweep of the slung load's frequency and cable tension; a mass-ratio sweep of the two modes; the Lyapunov exponent against release amplitude; and the prediction horizon of a chaotic run against solver tolerance and against the size of the initial error. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `double_pendulum.png`, the hand derivation of M(q) and f(q, q̇) by Lagrange and of the accelerations by Newton–Euler, the frequency equation and its two limits, the four experiment tables with your predictions, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the equations of motion are derived by both routes and shown to agree; all public tests pass; the slung load's two periods are given with their mode shapes and the amplitude at which the linear period is 1 % wrong; the note states a prediction horizon and what sets it. |
| **References** | `@mit-ocw-16-07-dynamics` L20 "Energy methods: Lagrange's equations" (https://ocw.mit.edu/courses/16-07-dynamics-fall-2009/resources/mit16_07f09_lec20/): Derivation of Lagrange's Equations in Cartesian Coordinates; Simple Pendulum by Lagrange's Equations; Double Pendulum by Lagrange's Equations; Derivation of Lagrange's Equation for General Coordinate Systems (required-free): the same double pendulum derived the same way; the Newton–Euler check is this lesson's own. `@mit-ocw-2-003sc-engineering-dynamics` Lecture 15: Introduction to Lagrange With Examples; Lecture 16: Kinematic Approach to Finding Generalized Forces; Recitation 9 Notes: Generalized Forces with Double Pendulum Example (required-free): generalized coordinates and generalized forces (Problem 5). `@libretexts-cline-variational-principles` Ch. 6 "Lagrangian Dynamics": 6.3 Lagrange Equations from d'Alembert's Principle; 6.4 Lagrange equations from Hamilton's Principle; 6.5 Constrained Systems; 6.13 The Lagrangian versus the Newtonian approach to classical mechanics; Ch. 4 "Nonlinear Systems and Chaos": 4.6 Differentiation Between Ordered and Chaotic Motion (required-free). `@mit-ocw-16-61-aerospace-dynamics` lecture 7 "Lagrange's Equations"; lecture 9 "Virtual Work and the Derivation of Lagrange's Equations" (supplementary-free). `@sympy-tutorial` "Introductory Tutorial": "SymPy Features" > "Calculus", "Solvers", "Matrices" (required-free): the symbolic derivation and `lambdify`. `@mit-ocw-18-03-differential-equations` notes 38 "Linearization near equilibria; the nonlinear pendulum" (required-free). `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (required-free). |
| **Extensions** | ★ Let the suspension point accelerate horizontally (a helicopter starting forward flight): add x_p(t) to both positions, derive the new terms in SymPy, and find the cable angle in a steady 1 m/s² acceleration. ★ Replace the point load by a rigid body (a container of length 6 m swinging about its sling point) and repeat the modes. ★ Add a Rayleigh dissipation function D = ½ c (θ̇₁² + θ̇₂²) and show the energy decays at the rate −2D. |
| **Time** | unmeasured |

## 1. The engineering problem

Helicopters carry loads that will not fit inside them: containers, bridge sections, firefighting buckets, other aircraft. The load hangs from a cable on a hook, and the hook, the cable and the load form a pendulum that the helicopter drags along. A slung load that swings is a hazard: it moves the helicopter's centre of mass, it can couple with the pilot's corrections (the pilot pushes the stick, the load swings back, the pilot corrects again a half-period late), and in the worst case it strikes the ground or the aircraft. The first thing a flight-test engineer wants to know is the load's natural periods, because those are the frequencies the pilot and the flight-control system must not feed. The second is what the cable carries as the load swings, because the hook and cable are certified to a limit load.

With a hook block of mass m₁ on a cable of length l₁, and the load's centre of mass m₂ a further l₂ below on its sling, the system is a **double pendulum**. Writing its equations by free-body diagrams is possible but tedious: there are two unknown link tensions, the accelerations are those of points on two rotating links, and one sign error ruins everything. **Lagrange's equations**, which you met as a consequence of Hamilton's principle in `P02.13` (Cline 6.4), give the same equations from two scalar energies, with no tensions in sight. This lesson derives the double pendulum both ways, symbolically and by hand, and uses each to check the other. Then it asks what the equations say: two normal modes at small amplitude, a slow change of period with amplitude, and, at large amplitude, chaos.

The concrete numbers (illustrative, not from a particular aircraft): m₁ = 80 kg, l₁ = 12 m, m₂ = 1500 kg, l₂ = 4 m, a fixed hover point. Success criteria, written first: the SymPy equations equal the hand-derived M(q) and f(q, q̇) to 10⁻¹²; the Lagrange accelerations equal the Newton–Euler ones to 10⁻¹⁰ at random states; the energy is conserved to 10⁻⁸ in a 10 s large-amplitude run; the numerical modes solve the closed-form frequency equation; a 10⁻³ rad motion started in a mode keeps its shape to 10⁻⁴ and oscillates at the mode's frequency.

## 2. Understand and derive

### 2.1 Generalized coordinates

A configuration of two point masses in a plane needs four numbers (x₁, y₁, x₂, y₂). The links impose two **holonomic constraints**, relations among the coordinates alone:

x₁² + y₁² = l₁²,  (x₂ − x₁)² + (y₂ − y₁)² = l₂².

Four coordinates minus two constraints leave **two degrees of freedom**. Any two numbers that fix the configuration and can be varied independently are **generalized coordinates**; the natural pair is the two link angles θ₁, θ₂ from the downward vertical:

x₁ = l₁ sin θ₁,  y₁ = −l₁ cos θ₁,  x₂ = x₁ + l₂ sin θ₂,  y₂ = y₁ − l₂ cos θ₂.   (1)

```
          O  pivot (hover point)
          |\
          | \  l1, theta1 from the vertical
          |  \
          |   o m1  (hook block)
          |   |\
          |   | \  l2, theta2
          |   |  \
          |   |   O m2  (load)
          y up, x right; theta positive counter-clockwise
```

These coordinates satisfy the constraints automatically. That is the whole trick: the constraint forces (the link tensions) do no work on any motion the coordinates allow, because such motions keep the link lengths fixed, so they will not appear in equations written in θ₁, θ₂.

### 2.2 Lagrange's equations

`P02.13` §2.5 stated **Hamilton's principle**: the motion between two configurations at two times makes the action S = ∫ L dt stationary, with the **Lagrangian** L = T − V. The Euler–Lagrange equation of the calculus of variations, applied once for each generalized coordinate, gives

**d/dt (∂L/∂q̇ᵢ) − ∂L/∂qᵢ = Qᵢ,  i = 1 … n,**   (2)

with Qᵢ = 0 when every force is either a constraint force or derived from the potential V. (When there are other forces, such as drag or a pilot's push, Qᵢ = Σₖ Fₖ · ∂rₖ/∂qᵢ is the **generalized force**: the work the forces do per unit change of qᵢ. It has units of force for a length coordinate and of moment for an angle.) A second route to (2), which explains why the constraint forces vanish, is d'Alembert's (Cline 6.3; 16.07 L20 derives it in Cartesian coordinates first, then for general coordinates): write Newton's law for each particle, dot it with the virtual displacement ∂rₖ/∂qᵢ δqᵢ, and sum. Constraint forces are perpendicular to every displacement the coordinates allow, so their dot products are zero, and the inertia terms rearrange into the left side of (2).

Two checks you can do by hand (and `euler_lagrange` passes in the tests): a spring–mass, L = ½mẋ² − ½kx², gives mẍ + kx = 0; a pendulum, L = ½ml²θ̇² + mgl cos θ, gives ml²θ̈ + mgl sin θ = 0.

### 2.3 The double pendulum by Lagrange

Differentiating (1), the velocities are ẋ₁ = l₁θ̇₁ cos θ₁, ẏ₁ = l₁θ̇₁ sin θ₁ and similarly for the second link. Squaring and adding, with Δ = θ₁ − θ₂ and cos θ₁ cos θ₂ + sin θ₁ sin θ₂ = cos Δ:

T = ½(m₁ + m₂) l₁² θ̇₁² + ½ m₂ l₂² θ̇₂² + m₂ l₁ l₂ θ̇₁ θ̇₂ cos Δ,
V = −(m₁ + m₂) g l₁ cos θ₁ − m₂ g l₂ cos θ₂.   (3)

The kinetic energy is a quadratic form in the rates, T = ½ q̇ᵀ M(q) q̇, with the **mass matrix**

**M(q) = [[(m₁ + m₂) l₁², m₂ l₁ l₂ cos Δ], [m₂ l₁ l₂ cos Δ, m₂ l₂²]].**   (4)

M is symmetric and positive definite for every configuration (its determinant is m₂ l₁² l₂² (m₁ + m₂ sin²Δ) > 0), and its off-diagonal term is the coupling between the links. Now apply (2). For θ₁: ∂L/∂θ̇₁ = (m₁ + m₂) l₁² θ̇₁ + m₂ l₁ l₂ θ̇₂ cos Δ; its time derivative contains θ̈₁, θ̈₂ and a term −m₂l₁l₂θ̇₂ sin Δ (θ̇₁ − θ̇₂); ∂L/∂θ₁ = −m₂ l₁ l₂ θ̇₁ θ̇₂ sin Δ − (m₁ + m₂) g l₁ sin θ₁. The θ̇₁θ̇₂ terms cancel, and the same happens for θ₂. The result is

**M(q) q̈ = f(q, q̇),  f₁ = −m₂ l₁ l₂ θ̇₂² sin Δ − (m₁ + m₂) g l₁ sin θ₁,  f₂ = m₂ l₁ l₂ θ̇₁² sin Δ − m₂ g l₂ sin θ₂.**   (5)

Units check: every entry of M q̈ is kg·m²·s⁻², a moment, and so is every term of f. The θ̇² terms are the **centripetal** coupling: the lower link feels the upper one's rotation and vice versa. `derive_equations` does this in SymPy with `linear_eq_to_matrix`, which splits each equation into its q̈ part (the rows of M) and the rest (f); the tests compare the result with (4) and (5) at random states. Solving the 2 × 2 system by Cramer's rule gives θ̈₁ and θ̈₂ explicitly; SymPy prints the familiar long fraction with the denominator l₁(m₁ + m₂ − m₂ cos²Δ).

### 2.4 The same accelerations by Newton–Euler

The independent route uses only free-body diagrams. Let eᵢ = (sin θᵢ, −cos θᵢ) be the unit vector along link i (pointing down the link) and nᵢ = (cos θᵢ, sin θᵢ) the perpendicular in the direction of increasing θᵢ. A point at distance l on a link turning at θ̇ with angular acceleration θ̈ has the acceleration of `P06.01`'s polar coordinates, l θ̈ n − l θ̇² e, so

a₁ = l₁(θ̈₁ n₁ − θ̇₁² e₁),  a₂ = a₁ + l₂(θ̈₂ n₂ − θ̇₂² e₂).   (6)

The upper link pulls m₁ towards the pivot with tension T₁; the lower link pulls m₁ down towards m₂ and m₂ up towards m₁ with tension T₂. Newton's second law for each mass:

m₂ a₂ = −T₂ e₂ + m₂ g,  m₁ a₁ = −T₁ e₁ + T₂ e₂ + m₁ g,  g = (0, −g).   (7)

That is four scalar equations, linear in the four unknowns (θ̈₁, θ̈₂, T₁, T₂). `newton_euler` assembles them as a 4 × 4 system and solves it. This route never mentions energy, a Lagrangian or M(q), so agreement with (5) is a genuine check of both. It also returns the tensions, which Lagrange's method eliminated: the price of not seeing constraint forces is not knowing them. Two limits fix the signs: hanging at rest, T₂ = m₂g and T₁ = (m₁ + m₂)g; hanging straight and turning together at ω, T₂ = m₂(g + (l₁ + l₂)ω²) and T₁ = T₂ + m₁(g + l₁ω²), the centripetal loads on top of the weights.

### 2.5 Energy

The Lagrangian (3) does not contain t explicitly. Multiply each equation (2) (with Q = 0) by q̇ᵢ and sum: Σ q̇ᵢ d/dt(∂L/∂q̇ᵢ) − Σ q̇ᵢ ∂L/∂qᵢ = 0. Since dL/dt = Σ(∂L/∂qᵢ q̇ᵢ + ∂L/∂q̇ᵢ q̈ᵢ), this is d/dt(Σ q̇ᵢ ∂L/∂q̇ᵢ − L) = 0. The conserved quantity, the **energy function** h = Σ q̇ᵢ ∂L/∂q̇ᵢ − L, equals 2T − (T − V) = T + V because T is a quadratic form in the rates (Euler's theorem on homogeneous functions). So

**E = T + V is constant along every motion,**   (8)

whatever the amplitude. It is the first diagnostic of any simulation of (5). With a moving suspension point, L would depend on t and E would not be conserved (Extension 1).

### 2.6 Small angles: two normal modes

At the hanging equilibrium q = 0, q̇ = 0. Keeping first-order terms in (5), cos Δ → 1, sin θ → θ and the θ̇² terms drop out (they are second order):

**M₀ q̈ + K₀ q = 0,  M₀ = M(0),  K₀ = diag((m₁ + m₂) g l₁, m₂ g l₂).**   (9)

Trying q = φ cos ωt gives the generalized eigenproblem K₀ φ = ω² M₀ φ (lesson `P06.14` treats it in general). Its characteristic equation det(K₀ − ω²M₀) = 0, after dividing by m₂ l₁ l₂, is

**m₁ l₁ l₂ ω⁴ − (m₁ + m₂) g (l₁ + l₂) ω² + (m₁ + m₂) g² = 0.**   (10)

Two limits check it. If m₂ → 0, (10) becomes m₁(l₁ω² − g)(l₂ω² − g) = 0: two independent pendulums, √(g/l₁) and √(g/l₂). If m₁ → 0 with m₂ fixed, the ω⁴ term vanishes and ω² = g/(l₁ + l₂): the hook has no mass, the two links line up, and the system is a single pendulum of length l₁ + l₂.

For the slung load, M₀ = [[227 520, 72 000], [72 000, 24 000]] kg·m² and K₀ = diag(185 934, 58 840) N·m/rad. The roots are

- **mode 1: ω₁ = 0.78667 rad/s, period 7.987 s, shape θ₂/θ₁ = +1.0129.** The cable and sling swing together, almost as one straight pendulum; a single 16 m pendulum would have a period of 8.026 s.
- **mode 2: ω₂ = 7.99635 rad/s, period 0.786 s, shape θ₂/θ₁ = −3.1196.** The load barely moves (θ₂ ≈ −3.12 θ₁ with l₂ = l₁/3 puts m₂ almost straight below the pivot) while the light hook jiggles sideways between two taut lines. Treating the load as fixed, the hook's lateral stiffness is T₁/l₁ + T₂/l₂ = (m₁ + m₂)g/l₁ + m₂g/l₂, which gives ω ≈ 7.881 rad/s, 1.4 % below the exact value.

The two modes are **M₀-orthogonal**, φ₁ᵀ M₀ φ₂ = 0 (a test), so any small motion is a sum of the two, each at its own frequency. The engineering reading: a pilot's slow corrections couple into the 8 s mode, and anything near 1.3 Hz (a rotor or airframe vibration, or a flight-control loop) can feed the hook mode.

### 2.7 Large amplitude and chaos

At larger amplitude the sin θ in (5) is smaller than θ and the swing slows. For mode 1, which is almost a rigid pendulum, the simple-pendulum correction from `P04.06`, ω ≈ ω₁(1 − a²/16) for an amplitude a in radians, is a good estimate (part 4 tests it).

The double pendulum has two degrees of freedom and, for a general motion, only one conserved quantity, the energy. One integral is not enough to confine a two-degree-of-freedom motion to a regular, predictable path, and at large energy the motion is **chaotic** (Cline 4.6 shows how to tell ordered from chaotic motion): two motions that start a distance δ₀ apart separate, on average, as δ(t) ≈ δ₀ e^{λt}, with λ > 0 the largest **Lyapunov exponent**. Equation (5) is perfectly deterministic; what fails is our ability to know the initial state well enough. If an answer is useless once the error reaches Δ, a simulation is good for

**t_h ≈ (1/λ) ln(Δ/δ₀),**   (11)

the **prediction horizon**. Its message is logarithmic: a hundredfold better initial state adds only ln(100)/λ to the horizon. At small amplitude the motion is close to the sum of two modes, nearby motions drift apart only linearly in time (through the small amplitude dependence of the frequencies), and λ = 0.

`lyapunov_estimate` measures λ by **Benettin's method**: integrate a reference motion and a neighbour d₀ = 10⁻⁸ away for τ = 0.5 s; record ln(d/d₀); pull the neighbour back to distance d₀ along the new separation (so that it never leaves the linear regime); repeat, and average. Rates are scaled by √(l₁/g) so that the four state components have the same units.

### 2.8 Where the model stops being valid

- **A taut cable.** The links are rigid only while they pull. If T₁ or T₂ falls to zero (a violent swing, a sudden climb), the cable goes slack and the load falls freely until it snatches the cable tight: a large impulsive load, outside this model. `newton_euler` returns the tensions, so you can check that they stay positive.
- **A fixed suspension point.** A real helicopter moves; its motion is an input that drives both modes. The coupled helicopter–load dynamics is the actual flight-test problem.
- **No aerodynamics.** Rotor downwash pushes on the load, and a container in forward flight has drag and can fly, yaw and spin on its sling. These are often the dominant effects.
- **A planar, point-mass, inextensible model.** Real loads swing in two directions and rotate; cables stretch (a vertical bounce mode); the sling is a set of legs, not one link.

## 3. Build the computational model

Open `starter/main.py`. `DoublePendulum`, the two parameter sets, `lagrange_accelerations` (a `lambdify` of your SymPy result), `positions`, `dominant_frequency` and `main` are given.

### Step 1: Lagrange's equations in SymPy

`euler_lagrange(L, qs, t)` returns d/dt(∂L/∂q̇) − ∂L/∂q for each q; SymPy differentiates with respect to `sp.diff(q, t)` directly. `double_pendulum_lagrangian()` builds (1), differentiates the positions with respect to t and returns T − V. `derive_equations()` applies `euler_lagrange`, replaces θ̈ᵢ, θ̇ᵢ and θᵢ by plain symbols, and splits the two equations with `sp.linear_eq_to_matrix`:

```python
eqs = [sp.expand(e.subs({th1: a1, th2: a2})) for e in eqs]      # after replacing th'' and th'
M, f = sp.linear_eq_to_matrix(eqs, [d1, d2])                    # eqs = M [d1, d2] - f
tidy = lambda e: sp.factor(sp.trigsimp(e))                       # entry by entry: fast
return M.applyfunc(tidy), f.applyfunc(tidy), (a1, a2, w1, w2, m1, m2, l1, l2, g)
```

Simplify entry by entry: `trigsimp` on the whole equations took about 5 s, on the five small entries well under a second. Decorate both `derive_*` functions with `functools.lru_cache` so the tests derive once. `derive_accelerations()` is `M.LUsolve(f)`.

### Step 2: the hand form and the integrator

`mass_matrix` and `generalized_forces` are (4) and (5). `rhs(t, y, p)` returns [θ̇₁, θ̇₂, θ̈₁, θ̈₂] by Cramer's rule on the 2 × 2 system. It is about as fast as the `lambdify`d SymPy result (within 20 % in a timing test), and it keeps the integrator independent of the symbolic step; it is what `simulate` integrates (DOP853, rtol = atol = 10⁻¹⁰ by default). `energy` is (3).

### Step 3: Newton–Euler

`newton_euler(y, p)` builds the 4 × 4 system (7) with (6) substituted, unknowns u = [θ̈₁, θ̈₂, T₁, T₂], and solves it with `np.linalg.solve`. Rows 0–1 are mass 2 (x and y), rows 2–3 mass 1.

### Step 4: modes and chaos

`linearized_matrices` is (9); `normal_modes` calls `scipy.linalg.eigh(K0, M0)` and scales each shape to θ₁ = 1; `normal_mode_frequencies_closed_form` solves (10), with the small root written as 2c/(−b + √(b² − 4ac)) to avoid cancellation. `lyapunov_estimate` is 2.7: integrate the reference and the neighbour together as an eight-state system for each interval.

Running `main()` prints

```
Lagrange (SymPy): M(q) = [[l1**2*(m1 + m2), l1*l2*m2*cos(a1 - a2)], [l1*l2*m2*cos(a1 - a2), l2**2*m2]]
                  f(q, q') = [-l1*(g*m1*sin(a1) + g*m2*sin(a1) + l2*m2*w2**2*sin(a1 - a2)), l2*m2*(-g*sin(a2) + l1*w1**2*sin(a1 - a2))]

slung load: hook 80 kg on 12 m, load 1500 kg on 4 m
  mode 1: w = 0.78667 rad/s (closed form 0.78667), period 7.987 s, shape th2/th1 = +1.0129
  mode 2: w = 7.99635 rad/s (closed form 7.99635), period 0.786 s, shape th2/th1 = -3.1196
  single pendulum of length l1 + l2: period 8.026 s

Lagrange vs Newton-Euler at 200 random states: max relative difference 1.5e-14
10 deg mode-1 swing over 60 s: energy drift 2.3e-13 of the swing energy; cable tension 0.9847-1.0305 x static weight
  simulated frequency 0.78517 rad/s vs linear 0.78667
lab pendulum released at th1 = th2 =   5.0 deg: finite-time Lyapunov exponent 0.012 1/s
lab pendulum released at th1 = th2 =  30.0 deg: finite-time Lyapunov exponent 0.038 1/s
lab pendulum released at th1 = th2 =  90.0 deg: finite-time Lyapunov exponent 0.850 1/s
lab pendulum released at th1 = th2 = 120.0 deg: finite-time Lyapunov exponent 1.396 1/s
```

and writes `double_pendulum.png`: the path of the lab pendulum's lower mass after a 120° release; the separation of two releases 10⁻⁸ rad apart on a log scale (a straight rise, then saturation at the size of the motion); and the slung load's two angles in a 10° mode-1 swing.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-11-lagrangian-mechanics`. The tests:

1. **Lagrange's equations**: spring–mass and simple pendulum symbolically; the double pendulum's SymPy M and f equal (4) and (5) at 20 random states for two parameter sets; with m₂ = 0 the first equation is a simple pendulum; M(q) symmetric positive definite.
2. **Newton–Euler**: the Lagrange accelerations, symbolic and hand, equal (7)'s at 50 random states for each of three parameter sets, to 10⁻¹⁰; the static and rigid-rotation tensions.
3. **Integration**: energy constant to 10⁻⁸ over 10 s from a large-amplitude start, and equal to ½Σm|v|² + Σmgy computed from the Cartesian positions with fourth-order differences; a run forwards and back returns to its start to 10⁻⁷.
4. **Modes**: eigenproblem roots equal (10) for three parameter sets; det(K₀ − ω²M₀) = 0; M₀-orthogonality; in-phase and out-of-phase shapes; the SymPy Jacobian at q = 0 equals −M₀⁻¹K₀; a vanishing m₂ gives √(g/l₁) and √(g/l₂); a 10⁻³ rad start in either mode oscillates at its ω to 0.2 % and keeps its shape to 10⁻⁴; a 20° swing is slower by (1 − a²/16) to 0.2 %.
5. **Chaos**: λ > 0.5 s⁻¹ for the 120° release, λ < 0.1 s⁻¹ for 5°, and the ratio above 10.
6. `main` writes the figure.

### Compare

The established tool here is SymPy's own mechanics package, `sympy.physics.mechanics`, which builds the Lagrangian from reference frames, points and particles and forms Lagrange's equations itself (`LagrangesMethod`). It is not used in `main.py`, but it is the natural cross-check: built for the same double pendulum and `lambdify`d, its accelerations differ from the lesson's by at most 1.8 × 10⁻¹⁴ (relative) over 400 random states of the slung load and the lab pendulum. The Newton–Euler route agrees to 1.5 × 10⁻¹⁴. Three derivations, two of them sharing nothing but the geometry, give the same numbers to round-off: the equations are right.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Every result below is a simulation of the model, not a measurement of a helicopter. Write each prediction down before you run it.

1. **Amplitude.** Start the slung load in its mode-1 shape at 5°, 10°, 20°, 30° and 45° and simulate 20 periods. Measure the frequency (`dominant_frequency`), the largest and smallest cable tension T₁ as a multiple of the static weight W = (m₁ + m₂)g = 15 494 N, and the largest departure from the mode shape, max|θ₂ − 1.0129 θ₁|/max|θ₂|. Predict the frequency shift from (1 − a²/16). *Finding:*

   | amplitude | ω (rad/s) | shift | (1 − a²/16) estimate | T₁/W range | shape departure |
   |---|---|---|---|---|---|
   | 5° | 0.78627 | −0.051 % | −0.048 % | 0.9962–1.0076 | 7.2 × 10⁻⁵ |
   | 10° | 0.78505 | −0.206 % | −0.190 % | 0.9847–1.0305 | 2.9 × 10⁻⁴ |
   | 20° | 0.78020 | −0.822 % | −0.762 % | 0.9393–1.1210 | 1.2 × 10⁻³ |
   | 30° | 0.77252 | −1.798 % | −1.713 % | 0.8651–1.2688 | 3.0 × 10⁻³ |
   | 45° | 0.75683 | −3.793 % | −3.855 % | 0.7043–1.5878 | 9.5 × 10⁻³ |

   The period error passes 1 % between 20° and 30° (near 22° by interpolation). The cable load is the sharper lesson: at 30° the hook sees 1.27 W at the bottom of each swing, a dynamic increment the static weight hides.
2. **Mass ratio.** Keep m₂ = 1500 kg, l₁ = 12 m, l₂ = 4 m, and vary the hook mass m₁ = 20, 80, 300, 1500, 6000 kg. Predict which period changes most. *Finding:* mode-1 periods 8.016, 7.987, 7.895, 7.594, 7.251 s; mode-2 periods 0.399, 0.786, 1.442, 2.597, 3.440 s; mode-2 shapes θ₂/θ₁ = −3.03, −3.12, −3.45, −5.16, −11.3. The slow mode hardly cares about the hook (it tends to the 16 m pendulum's 8.026 s as m₁ → 0); the hook mode scales roughly as √m₁ when the hook is light, the stiffness-on-a-small-mass picture of 2.6.
3. **Chaos against amplitude.** For the lab pendulum (1 kg, 1 m, equal), release from rest at θ₁ = θ₂ = 5°, 30°, 60°, 90°, 120°, 150° and compute λ over 30 s. *Finding:* 0.012, 0.038, 0.074, 0.850, 1.396, 1.354 s⁻¹. The small values are not zero only because a finite-time estimate of a regular motion decays like (ln t)/t; the jump between 60° and 90° is the onset of chaos for this family of releases.
4. **Solver tolerance and the horizon.** Release the lab pendulum at 120° and integrate 20 s with rtol = atol = 10⁻⁶, 10⁻⁸, 10⁻¹⁰, 10⁻¹². For each, record the largest energy error and the time at which θ₂ first differs from a 10⁻¹³ run by 0.1 rad. *Finding:* energy errors 2.1 × 10⁻⁴, 1.1 × 10⁻⁵, 1.3 × 10⁻⁷, 6.4 × 10⁻¹⁰ J (of E = 14.71 J); horizons 10.70, 12.70, 15.49, 18.07 s. A run that conserves energy to ten digits still loses the actual trajectory after 18 s: energy conservation is necessary, not sufficient.
5. **Initial error and the horizon.** With rtol = 10⁻¹³, perturb θ₁(0) by 10⁻⁴, 10⁻⁶, 10⁻⁸, 10⁻¹⁰ rad. *Finding:* horizons 6.93, 10.41, 12.46, 15.11 s, about 2.7 s gained per factor of 100, i.e. an effective growth rate of ln(100)/2.7 ≈ 1.7 s⁻¹, consistent in size with experiment 3's λ (the horizon also depends on the threshold and on where in the motion the error grows fastest).

Raw data format: CSV tables `amplitude_deg, omega_radps, shift_pct, estimate_pct, T1_min_W, T1_max_W, shape_departure` (1); `m1_kg, T1_s, T2_s, shape1, shape2` (2); `release_deg, lambda_per_s` (3); `rtol, max_dE_J, horizon_s` (4); `delta0_rad, horizon_s` (5); plus `double_pendulum.png`.

## 5. Validate and improve

| quantity | analysis (§2) | computation | difference | cause |
|---|---|---|---|---|
| equations of motion | (4), (5) by hand | SymPy `derive_equations` | 0 (symbolic) | — |
| accelerations at random states | Newton–Euler (7) | Lagrange (5) | 1.5 × 10⁻¹⁴ relative | round-off |
| same | `sympy.physics.mechanics` | Lagrange (5) | 1.8 × 10⁻¹⁴ relative | round-off |
| mode frequencies (rad/s) | (10): 0.78667, 7.99635 | `eigh`: 0.78667, 7.99635 | < 10⁻¹⁰ | — |
| hook mode (rad/s) | lateral-stiffness estimate 7.881 | 7.996 | 1.4 % | load treated as fixed |
| 10° swing frequency (rad/s) | 0.78667 (1 − a²/16) = 0.78517 | 0.78505 (20 periods) | 0.015 % | estimate treats mode 1 as a rigid pendulum |
| energy, 10° swing, 60 s | constant (8) | drift 2.3 × 10⁻¹³ of the swing energy | — | solver tolerance |
| energy, 120° chaotic release, 20 s | constant (8) | 1.3 × 10⁻⁷ J at rtol 10⁻¹⁰ | 10⁻⁸ relative | solver tolerance |

The independent check required of this lesson, Lagrange against Newton–Euler for the same system, holds to round-off at every state tested, and the energy is conserved to the solver tolerance; the linear modes are exactly what the small-amplitude simulation does. The engineering answer: **the slung load swings in a 7.99 s pendulum mode and a 0.79 s hook mode; in a 10° swing the cable load varies by −1.5 % to +3.1 % about the static weight, and by +27 % at 30°; the linear period is 1 % wrong near 22°.** For a chaotic motion, the model's equations are exact but its predictions are not: each factor of 100 in initial precision or solver tolerance buys only 2–3.5 s (2.7 s on average).

What to improve first is the physics the model leaves out, not the numerics: the helicopter's motion as an input (Extension 1), then aerodynamic forces on the load, then the cable's elasticity. What not to improve: the solver tolerance beyond 10⁻¹⁰ for the small-amplitude questions, which are already exact to ten digits.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (3)–(5), (7), (8) and (10); `double_pendulum.png`; the five experiment tables with your predictions; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the two periods and mode shapes, the cable-load range in a 10° and a 30° swing, the amplitude at which the linear model is 1 % wrong, and a prediction horizon for the chaotic release.

### Adversarial self-checks

- **What if the load mass doubles?** Mode 1 hardly moves (it is a pendulum; mass cancels), but the hook mode stiffens: the cable tensions double while the hook mass does not. Run it before you believe it.
- **What if the amplitude goes to zero?** The motion is exactly the sum of the two modes of (9); λ → 0; the cable tension → its static value. All three are tests.
- **Which assumption fails first in the real world?** For the slung load, the fixed hover point and the absence of downwash: both drive the modes far more than the nonlinearity at 10°. For the lab pendulum, joint friction, which removes energy and with it the chaos.
- **What measurement could falsify the model?** A video of a scaled slung load (a weight on a string hung from a second weight) released in the mode-1 shape: the measured period and shape should match (10) and the mode-1 ratio; a second frequency in the record means the release was not a pure mode.
- **Does the result survive a change of method?** Lagrange by hand, Lagrange in SymPy, `sympy.physics.mechanics` and Newton–Euler agree to 10⁻¹⁴. The chaotic trajectory does not survive a change of tolerance beyond 10–18 s, and should not be expected to; its statistics (λ, the energy) do.
- **Which decision changes if the uncertainty doubles?** If the hook-block mass were known only to ±20 %, the hook-mode period would move by about ±10 % (it scales roughly as √m₁); a flight-control notch filter placed at 1.27 Hz would then need a wider band, while the 8 s mode, set by the 16 m length, stays put.

## Check yourself

1. How many degrees of freedom does a planar double pendulum with rigid links have, and why do the link tensions not appear in Lagrange's equations?
2. Write the Lagrangian of a simple pendulum and derive its equation of motion.
3. What is the physical meaning of the off-diagonal term m₂l₁l₂ cos Δ in M(q), and when does it vanish?
4. Why is T + V conserved for the double pendulum? Name one change to the system that would break it.
5. What does the frequency equation (10) reduce to when m₂ → 0, and when m₁ → 0? Interpret each.
6. In the slung load's second mode, which mass moves most, and why is its frequency set mainly by the hook mass?
7. What is a Lyapunov exponent, and how does it limit the prediction horizon of a simulation?
8. A simulation conserves energy to 10⁻¹⁰. Does that mean its trajectory is right? Answer for a small and for a chaotic motion.
