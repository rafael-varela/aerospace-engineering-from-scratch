# P08.02  Fluid Kinematics

> A streamline is a snapshot of directions, a pathline is one particle's history, a streakline is what smoke shows. In a steady flow they are the same curve; the moment the flow changes in time, they part company.

| | |
|---|---|
| **Learning outcomes** | 1. Describe a flow field in the Eulerian (fields at fixed points) and Lagrangian (following particles) views, and define streamlines, pathlines and streaklines precisely. 2. Derive the closed-form pathline, streamline and streakline of an oscillating stream, and prove that the three coincide in a steady flow. 3. Implement and test a line plotter for any unsteady 2-D field by RK4 integration, checked against the closed forms, against an independent ODE solver and against steady-flow coincidence. 4. Compute the velocity-gradient tensor, split it into strain rate and rotation, and evaluate the divergence, the vorticity and the material derivative, verified against the acceleration of a computed pathline. 5. Decide when a smoke or dye photograph can be read as a streamline picture, from the reduced frequency ωL/U. |
| **Prerequisites** | `P01.11` (vector fields, line integrals, the velocity field as a vector field), `P01.08` (partial derivatives and the chain rule), `P01.12` (divergence theorem), `P03.05` (RK4 and observed order), `P02.02` (systems of ODEs), `P08.01` (the continuum). Skill check: you can integrate dx/dt = f(x, t) with RK4 and write the Jacobian matrix of a 2-D vector field. |
| **Engineering question** | A wind tunnel's gust generator swings a 10 m/s stream by ±5° at 2 Hz. Can the smoke lines photographed over a model be read as streamlines? The answer: only where the reduced frequency k = ωL/U is small. The smoke departs from the instantaneous streamline by about k/2 of the flow deflection: below 5 % only over the first 80 mm of the model (k ≤ 0.1), 19 % over a 0.3 m model, 85 % over the light aeroplane's full 1.5 m chord. |
| **Assumptions and limits** | Kinematics only: the velocity field is given, and no forces or conservation laws are imposed beyond what the example fields satisfy (all are divergence-free). Smoke and dye particles are ideal tracers: no inertia, no diffusion, no buoyancy. 2-D fields; RK4 with fixed steps; velocity gradients by central differences of step 10⁻⁵. Real smoke diffuses and has finite inertia, which blurs streaklines over long distances. |
| **Analytical work** | The three line definitions as ODEs; their closed forms for the oscillating stream and the sine-wave streakline; the uniqueness argument for steady flow; the material derivative from the chain rule; the decomposition G = S + W; the divergence as rate of area change and the vorticity as twice the rotation rate; the small-k expansion of the smoke error. |
| **Computational deliverable** | `starter/main.py`: `pathline`, `streakline`, `streamline`, `oscillating_pathline_exact`, `oscillating_streakline_exact`, `velocity_gradient`, `decompose_gradient`, `material_derivative`, `particle_acceleration`, `smoke_error` (the example fields, `rk4_step`, `polygon_area`, `distance_to_curve`, `observed_order` and `main` are given). Tests prove: RK4 pathlines and streaklines match the closed forms to 10⁻⁸ m and converge at order 4.0 ± 0.2; the closed forms satisfy their ODE; streamlines of the oscillating stream are straight with the instantaneous slope; unsteady lines differ by more than 50 mm; the smoke error tends to k/2; in steady stagnation flow all three lines lie on xy = const to 10⁻⁸; a decaying Taylor–Green vortex is unsteady but its lines coincide; streamlines depend only on direction; pathlines agree with `solve_ivp` at fourth order; stagnation flow is pure strain, shear is half strain and half rotation, rigid rotation has vorticity 2Ω; the example fields are divergence-free; the material derivative equals a pathline's second derivative to 10⁻⁵ and vanishes for a conserved quantity; a material loop keeps its area to 10⁻³. |
| **Practical activity** | **computational experiment**: a reduced-frequency sweep of the smoke-line error for a gust tunnel, a convergence study of the streakline release count and the pathline step, and a strain study of a material loop in stagnation and compressing flows. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `flow_lines.png`, the derivations of the three lines for the oscillating stream and of the k/2 error, the three experiment tables, and a half-page note on when your smoke pictures are streamline pictures. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the three lines are derived for the oscillating stream and their coincidence in steady flow is proved; all public tests pass; the k/2 law is derived and checked numerically; the note gives the largest model length for 5 % error at your tunnel's conditions. |
| **References** | `@mit-ocw-16-01-unified-engineering-fluids` Fluids Lecture 8 notes "Streamlines, Pathlines, Streaklines" (1. Streamlines; 2. Pathlines; 3. Streaklines) and Lecture 10 notes "Substantial Derivative", https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f08_new_fall/ (required-free: the definitions and the material derivative at exactly this level). `@mit-ocw-2-25-advanced-fluid-mechanics` "The Continuum Viewpoint and the Equations of Motion" class notes: Overview of Lagrangian and Eulerian Descriptions; Particle Kinematics; Geometric Interpretation of Fluid Kinematics In Steady Shear Flow (required-free: strain and rotation). `@libretexts-smyth-all-things-flow` Ch. 5 "Fluid Kinematics": 5.1 Lagrangian and Eulerian descriptions; 5.3 Rotation and strain- the relative motion of two nearby particles (supplementary-free). `@nasa-glenn-beginners-guide-aeronautics` "Definition of Streamlines" (supplementary-free). `@leishman-intro-aerospace-flight-vehicles` chapter "Types of Fluid Flows": Steady Versus Unsteady Flows (supplementary-free; link only). |
| **Extensions** | ★ Add tracer inertia: integrate a small particle's equation of motion with Stokes drag (`P08.05`) and find the particle size at which smoke lags the streakline by 1 mm over the model. ★ Compute streaklines for the Kármán vortex street of `P08.16` from a model field and compare with a cylinder-wake photograph. ★ Replace fixed-step RK4 by `solve_ivp` with dense output and compare the cost of a 500-release streakline. |
| **Time** | unmeasured |

## 1. The engineering problem

Flow visualisation is the oldest diagnostic in aerodynamics. A smoke wand upstream of a model in a tunnel, dye from a port in a water tunnel, tufts of wool on a wing in flight: each turns an invisible velocity field into a picture. The pictures are only useful if you know which curves they show. A smoke wand releases tracer at a fixed point, so a photograph shows the present position of everything that has passed that point: a **streakline**. What the aerodynamicist usually wants is the **streamline**, the curve tangent to the velocity everywhere at one instant, because streamlines are what potential-flow and CFD results are drawn as and what the pressure follows (Bernoulli, `P08.04`). A long exposure of one tracer shows a third curve, its **pathline**.

In a steady flow the three are the same curve, so the distinction never matters in a steady tunnel. It does matter in a gust tunnel, behind a pitching model, in a vortex wake, and in any flow that changes during the time a tracer takes to cross the region of interest. The question this lesson answers is quantitative: for a gust generator that swings a 10 m/s stream by ±5° at 2 Hz, how far do the smoke lines over a model depart from the streamlines, and over what length of model can they be trusted?

Success criteria, written first: the RK4 pathlines and streaklines reproduce the closed forms of the oscillating stream to 10⁻⁸ m and converge at fourth order; all three lines coincide (to 10⁻⁸) in a steady stagnation-point flow; the velocity-gradient decomposition gives the exact strain and rotation of three standard flows; the material derivative equals the second derivative of a computed pathline; and the smoke error follows k/2 for small reduced frequency.

## 2. Understand and derive

### 2.1 Two ways to describe a flow

The **Eulerian** description gives the fields at fixed points: u(x, y, t), v(x, y, t), p(x, y, t). The **Lagrangian** description follows fluid particles: x_P(t), y_P(t) for each particle P. The bridge is that a particle's velocity is the field's value where it is (16.01 Lecture 8; 2.25 "Overview of Lagrangian and Eulerian Descriptions"):

**dx_P/dt = u(x_P, y_P, t),  dy_P/dt = v(x_P, y_P, t).**   (1)

### 2.2 The three lines

- **Streamline** at time t₀: a curve everywhere tangent to the velocity at that instant, dy/dx = v/u at fixed t₀. In arc-length form, dx/ds = u/|u|, dy/ds = v/|u|: the speed does not matter, only the direction.
- **Pathline**: the solution of (1) for one particle over time.
- **Streakline** at time t from a source at (x_s, y_s): the set of positions at time t of all particles that passed the source at earlier times τ. Each point is the pathline of the particle released at τ, evaluated at t.

**Steady flow.** If u does not depend on t, (1) is an autonomous ODE. Through each point it has exactly one solution curve (uniqueness, for a smooth u), and that curve is independent of when the particle started. The pathline through a point is then the streamline through it, and every particle released at a source travels the same curve: **all three coincide** (Problem 2). They also coincide in some unsteady flows: if only the speed changes in time and never the direction, as in the decaying Taylor–Green vortex u = sin x cos y e^{−2νt}, v = −cos x sin y e^{−2νt}, the streamlines are fixed curves and particles follow them. The tests check both cases.

### 2.3 The oscillating stream, worked fully

Take u = U, v = V cos ωt: uniform in space, unsteady in time, and divergence-free. It is the simplest model of a gust tunnel whose flow direction swings by ±arctan(V/U).

**Streamline** at time t through the origin: slope v/u = (V/U) cos ωt, a straight line. **Pathline** of the particle at (x₀, y₀) at t₀: integrate (1),

**x = x₀ + U(t − t₀),  y = y₀ + (V/ω)(sin ωt − sin ωt₀),**   (2)

a sine wave of amplitude V/ω. **Streakline** at time t from a source at the origin: the particle released at τ is at x = U(t − τ), y = (V/ω)(sin ωt − sin ωτ). Eliminating τ = t − x/U:

**y = (V/ω)[sin ωt − sin ω(t − x/U)],**   (3)

a sine wave of wavelength 2πU/ω, frozen into the smoke and carried downstream at U. Its slope at the source is (V/U) cos ωt, the slope of the present streamline: smoke leaves along the current flow direction and then keeps the direction it was given. For U = 10 m/s, V = U tan 5° = 0.8749 m/s and ω = 4π rad/s (2 Hz), the wavelength is 5.000 m and the pathline amplitude is 69.62 mm.

**How far does the smoke depart from the streamline?** Expand (3) for small x: sin ω(t − x/U) = sin ωt − (ωx/U) cos ωt − ½(ωx/U)² sin ωt + …, so

y_streak − y_stream = (V/ω)·½(ωx/U)² sin ωt + O(x³).

Divide by the streamline's lateral excursion over a length L, VL/U, and take the worst phase at x = L:

**error ≈ k/2,  k = ωL/U (the reduced frequency).**   (4)

The smoke is a faithful streamline picture over a length L only if k ≪ 1. For 5 % at 2 Hz and 10 m/s, k ≤ 0.1 gives **L ≤ 0.080 m**.

### 2.4 The velocity gradient: strain and rotation

Near a point, u(x + δx) ≈ u(x) + G δx with the **velocity-gradient tensor** G_ij = ∂u_i/∂x_j. Split it into symmetric and antisymmetric parts (2.25 "Particle Kinematics"; Smyth 5.3):

**G = S + W,  S = (G + Gᵀ)/2 (strain rate),  W = (G − Gᵀ)/2 (rotation).**   (5)

Following the edges of a small rectangle (Problem 2): its area changes at the rate **tr G = ∇·u** per unit area (the **divergence**; zero for an incompressible flow), and its edges turn at the average angular velocity **ω/2**, where **ω = ∂v/∂x − ∂u/∂y** is the **vorticity**. The eigenvalues of S are the principal strain rates: the rates of stretching along two perpendicular directions. Three examples settle the intuition:

| flow | G | divergence | vorticity | principal strain rates |
|---|---|---|---|---|
| stagnation u = ax, v = −ay | [[a, 0], [0, −a]] | 0 | 0 | ±a |
| shear u = γy | [[0, γ], [0, 0]] | 0 | −γ | ±γ/2 (at 45°) |
| rigid rotation u = −Ωy, v = Ωx | [[0, −Ω], [Ω, 0]] | 0 | 2Ω | 0 |

Simple shear, the flow next to every wall, is half strain and half rotation. That split is why a boundary layer both deforms and spins its fluid elements, and why vorticity is born at walls (`P08.09`).

### 2.5 The material derivative

The rate of change of a quantity f seen by a particle is, by the chain rule along (1),

**Df/Dt = ∂f/∂t + u ∂f/∂x + v ∂f/∂y = ∂f/∂t + (u·∇)f**   (6)

(16.01 Lecture 10, "Substantial Derivative"). The first term is the **local** change at a fixed point; the second is the **convective** change from moving to where f is different. Applied to the velocity itself, Du/Dt is the particle's acceleration, the left-hand side of Newton's second law for a fluid element in every lesson from `P08.04` on. A steady flow has no local acceleration, but its particles still accelerate wherever the speed or direction changes along a streamline: in the stagnation flow with a = 1 s⁻¹, the particle at (2, 1) m accelerates at (a²x, a²y) = (2, 1) m/s².

### 2.6 Where the model stops

Real tracers have inertia and diffuse. Smoke particles a micrometre across follow air flows faithfully at these speeds, but heavier particles lag in strongly curved flow, and smoke diffuses into a widening band over metres. Dye in water diffuses slowly but is denser than water if not matched. Far from the source a streakline is a record of the flow's history *and* of the tracer's physics.

## 3. Build the computational model

Open `starter/main.py`. The example fields, `rk4_step`, `polygon_area`, `distance_to_curve`, `observed_order` and `main` are given. Every field is a function `vel(x, y, t) -> (u, v)` that accepts arrays, so one RK4 step can move a thousand particles at once.

### Step 1: pathlines and streamlines

`pathline(vel, x0, y0, t0, t1, n)` is n calls of `rk4_step`. `streamline(vel, x0, y0, t, length, n)` is the same with the velocity replaced by its unit vector at the frozen time t, `unit = lambda x, y, s: (u/|u|, v/|u|)`, stepping in arc length; it raises at a stagnation point, where the direction is undefined.

### Step 2: streaklines

`streakline(vel, xs, ys, t_start, t_now, n)` keeps an array of released particles, advances them all by one RK4 step, then appends a new particle at the source, n times:

```python
x, y = np.array([xs]), np.array([ys])
for i in range(n):
    x, y = rk4_step(vel, x, y, t_start + i * dt, dt)
    x, y = np.append(x, xs), np.append(y, ys)
```

The oldest particle is first; the newest is at the source.

### Step 3: closed forms, gradients and derivatives

`oscillating_pathline_exact` and `oscillating_streakline_exact` are (2) and (3). `velocity_gradient` evaluates G by central differences of step h = 10⁻⁵ (exact for linear fields up to round-off). `decompose_gradient` returns S, W, the divergence, the vorticity and the eigenvalues of S. `material_derivative(f, vel, x, y, t)` is (6) by central differences; `particle_acceleration` applies it to each velocity component. `smoke_error` evaluates the closed forms on a grid of x and returns the largest gap divided by VL/U.

Running `main()` prints

```
gust tunnel: U = 10.0 m/s, flow angle +/- 5.0 deg (V = 0.8749 m/s), 2.0 Hz (omega = 12.566 rad/s)
  streakline wavelength 2 pi U/omega = 5.000 m, pathline amplitude V/omega = 69.62 mm
  streakline from RK4 (200 releases) vs closed form: max error 3.77e-10 m
  pathline from RK4 (200 steps) vs closed form: max error 3.77e-10 m
  pathline to t = 0.3 s, RK4 errors [2.88e-07 1.80e-08 1.12e-09 7.01e-11], orders [4. 4. 4.]
  smoke-line error / streamline excursion over a model of length L (worst phase):
     L (m)   k = omega L / U   error
     0.05        0.063        0.031
     0.15        0.188        0.094
     0.30        0.377        0.188
     1.50        1.885        0.853
     5.00        6.283        1.000

steady stagnation flow a = 1/s: x y along streakline 5.6e-13 spread, pathline 5.6e-13, streamline 4.1e-12 (all on x y = 0.2)

velocity gradients (central differences):
  stagnation a = 1/s     div +0.000  vorticity +0.000  principal strain rates [-1.  1.]
  shear gamma = 1/s      div +0.000  vorticity -1.000  principal strain rates [-0.5  0.5]
  rotation Omega = 1/s   div +0.000  vorticity +2.000  principal strain rates [0. 0.]
  particle acceleration in steady stagnation flow at (2, 1): [2. 1.] m/s^2 (a^2 x, a^2 y)
```

and writes `flow_lines.png`: the three different lines of the oscillating stream, and the coinciding lines of the stagnation flow over its hyperbolic streamlines.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-02-fluid-kinematics`. The 19 tests:

1. **Closed forms**: RK4 pathline and streakline against (2) and (3) (10⁻⁹ m); pathline order 3.8–4.2; (2) satisfies (1) by central differences; streamlines straight with slope (V/U) cos ωt and arc length exact.
2. **Unsteadiness**: streakline, pathline and streamline at t = 0.3 s more than 50 mm apart; smoke error / k = 0.50 ± 2 % for small k.
3. **Steadiness**: in stagnation flow all three lines lie on xy = 0.2 to 10⁻⁸ and the streakline particles sit at the pathline's points; the Taylor–Green pathline and streamline stay on one streamfunction contour; streamlines of fields differing only in speed are identical.
4. **Independent integrator**: pathline in a pulsating, sheared stagnation flow against `scipy.integrate.solve_ivp` (rtol 10⁻¹²), error below 10⁻⁷ and order above 3.7.
5. **Gradients**: stagnation pure strain; shear S and W with ±γ/2; rigid rotation S = 0 and vorticity 2Ω; three example fields divergence-free.
6. **Material derivative**: equals the pathline's second difference to 10⁻⁵; (a²x, a²y) in steady stagnation flow and purely local in the uniform oscillating stream; D(xy)/Dt = 0 in stagnation flow.
7. **Area**: a 400-point loop advected for 2 s in the Taylor–Green field keeps its area to 10⁻³.
8. `main` writes the figure.

### Compare

RK4 (`rk4_step`) and SciPy's adaptive Runge–Kutta (`solve_ivp`, an independent implementation with error control) agree on a pathline through a non-uniform unsteady field to 2.7 × 10⁻⁸ m at 80 steps, with the observed order 4.09, 4.05, 4.03 as the step halves. For the oscillating stream the closed forms are exact: the RK4 pathline error is 3.8 × 10⁻¹⁰ m over 1 s at 200 steps, and the streakline error falls from 1.6 × 10⁻⁶ to 3.8 × 10⁻¹⁰ m as the release count goes from 25 to 200 (fourth order, since each released particle is an RK4 pathline). The decomposition was checked against hand values for three flows, and the material derivative against a finite-difference acceleration of a computed pathline: (1.194265, 1.111814) m/s² against (1.194264, 1.111813) m/s², a relative difference of 9 × 10⁻⁷, the O(dt²) error of the second difference. At that point the local part ∂u/∂t is −0.391 m/s²: both parts of (6) are present.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Reduced frequency.** For the 10 m/s, ±5° gust tunnel, compute the worst-phase smoke error over model lengths 0.05–5 m at 2 Hz, and over a 0.3 m model at 0.5, 1, 2, 4 and 8 Hz. Predict from (4). *Finding:* the table above; and at 0.3 m: errors 0.047, 0.094, 0.188, 0.371, 0.707 for k = 0.094, 0.188, 0.377, 0.754, 1.51. The k/2 law holds to 1 % up to k ≈ 0.4 and then saturates (the error cannot exceed the full excursion). A root find gives 5 % at L = 0.0796 m, k = 0.1000.
2. **Release count and step.** Streakline error against (3) at t = 1 s for 25, 50, 100, 200 releases; pathline error at t = 0.3 s for 10–80 steps. Predict the orders. *Finding:* 1.55 × 10⁻⁶, 9.6 × 10⁻⁸, 6.0 × 10⁻⁹, 3.8 × 10⁻¹⁰ m (order 4.0); pathline 2.88 × 10⁻⁷ … 7.0 × 10⁻¹¹ m (orders 4.00). A streakline costs n²/2 particle-steps; at 200 releases that is 20 000, still trivial.
3. **Stretching a material loop.** A circle of radius 0.1 m at (0.3, 1.0) m in stagnation flow (a = 1 s⁻¹) for 2 s, and in the compressing-expanding field u = ax, v = 0. Predict area and perimeter. *Finding:* stagnation flow: area ratio 1.000000000003, perimeter ratio 4.71, x-extent stretched by 7.389 = e²; u = ax: area ratio 7.389 = e^{(∇·u)t}. Incompressible means area-preserving, not shape-preserving: strain stretches material lines exponentially, which is how mixing works (Problem 5).

Raw data format: CSV tables `f_Hz, L_m, k, smoke_error` (1); `n, streak_err_m, path_err_m` (2); `field, t_s, area_ratio, perimeter_ratio` (3); plus `flow_lines.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| pathline of oscillating stream | (2) | RK4, 200 steps | 3.8 × 10⁻¹⁰ m | RK4 truncation, order 4 |
| streakline | (3) | 200 releases | 3.8 × 10⁻¹⁰ m | as above |
| pathline in pulsating flow | `solve_ivp` (rtol 10⁻¹²) | RK4, 80 steps | 2.7 × 10⁻⁸ m | order 4.03 |
| steady lines coincide | xy = 0.2 exactly | spread 5.6 × 10⁻¹³ (streak, path), 4.1 × 10⁻¹² (stream) | < 10⁻¹¹ | round-off |
| smoke error at small k | k/2 | 0.49999 k at k = 0.013 | 10⁻⁵ | O(k²) term |
| Du/Dt | pathline second difference | `particle_acceleration` | 9 × 10⁻⁷ relative | O(dt²) of the difference |
| loop area, Taylor–Green, t = 2 s | conserved | ratio 0.999998 | 2 × 10⁻⁶ | 400-point polygon |

The engineering answer: **in the 10 m/s gust tunnel at 2 Hz, smoke lines can be read as streamlines to 5 % only over the first 80 mm downstream of the wand (k = ωL/U ≤ 0.1); over a 0.3 m model they are off by 19 % of the flow deflection and over a 1.5 m chord by 85 %.** The fix is not a better smoke generator: measure the instantaneous velocity field (particle image velocimetry) and compute streamlines from it, or compare the photographs with computed streaklines rather than streamlines.

What to improve: tracer inertia and diffusion (Extension 1), which add a second error to the one computed here. What not to improve: the RK4 step; the kinematic error is ten orders of magnitude larger than the numerical one.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), (4), (5) and (6); `flow_lines.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the largest model length for which your tunnel's smoke pictures are streamline pictures, and the measurement you would make instead.

### Adversarial self-checks

- **What if the frequency goes to zero?** k → 0, the flow is quasi-steady, and the three lines coincide over any length: (4) gives zero error and the code's smoke error falls linearly with k.
- **What if the velocity is zero somewhere?** A streamline has no direction at a stagnation point; `streamline` raises rather than returning a meaningless curve. Pathlines and streaklines approach it but take infinite time to reach it.
- **Which assumption fails first in the real world?** Tracer fidelity: smoke diffuses and heavy particles lag, so the measured streakline differs from the computed one by more than the RK4 error long before any kinematic subtlety.
- **What measurement could falsify the analysis?** Phase-locked photographs of smoke at several gust phases: the smoke's lateral position at x = L should follow (3), a sine of wavelength 2πU/ω; a different wavelength would mean the tunnel's convection speed is not U (a boundary-layer or blockage effect).
- **Does the result survive a change of integrator?** RK4 and `solve_ivp` agree to 3 × 10⁻⁸ m; the closed forms to 4 × 10⁻¹⁰ m.
- **Which decision changes if the uncertainty doubles?** If the gust amplitude is uncertain by a factor of 2, the smoke error *fraction* does not change (it depends on k only), but the absolute gap does: a decision based on the absolute lateral offset (say, whether smoke crosses a probe) would change.

## Check yourself

1. Define a streamline, a pathline and a streakline.
2. When do the three coincide?
3. Which of them does smoke from a fixed wand show?
4. Write the material derivative and say what its two terms mean.
5. Can a particle accelerate in a steady flow? Give an example.
6. What does the divergence of the velocity measure?
7. Define the strain-rate tensor, the rotation tensor and the vorticity.
8. How much of a simple shear flow's velocity gradient is strain and how much is rotation?
