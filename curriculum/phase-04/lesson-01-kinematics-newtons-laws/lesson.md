# P04.01  Kinematics & Newton's Laws 🔬

> Draw the forces, write F = ma, and the trajectory follows. The whole of flight dynamics is this sentence with more forces.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the constant-acceleration kinematics of a projectile from Newton's second law, state the three assumptions behind it, and write the time of flight, range and apex height in closed form. 2. Draw the free-body diagram of a body in air, derive the quadratic-drag equation of motion, its terminal speed, and the closed-form solutions for purely vertical and purely horizontal motion, and expand the vertical one to first order in the drag. 3. Implement the two-dimensional equations of motion with drag and integrate them with `scipy.integrate.solve_ivp` to a terminal event, verifying the three closed-form limits and the convergence of the result with the solver tolerance. 4. Compose velocities for relative motion (the wind triangle), solve the inverse problem (which heading makes good a track) and prove that a round trip in wind always takes longer than in still air. 5. Design, predict, run and reduce a drop experiment that measures the terminal speed, hence the drag coefficient, of a ping-pong ball with a phone camera, with a prediction and an uncertainty budget written before the first drop. |
| **Prerequisites** | `P00.13` (safety box), `P00.14` (uncertainty: Type A and B, propagation, E_n), `P00.15` (the drop timing method and kit, reused here), `P01.07` (vectors: components, addition, magnitude), `P02.04` (a second-order ODE as a first-order system), `P03.05` (what `solve_ivp` does and which method is which). Skill check: you can resolve a vector into components with a sine and a cosine, and you can write x″ = f(x, x′) as a first-order system. |
| **Engineering question** | A ping-pong ball (2.7 g, 40 mm) dropped from 2.0 m takes 4.7 % longer to land than free fall predicts, and tossed at 6 m/s and 45° it lands at 73 % of the vacuum range. Can a phone camera at 30 or 60 fps measure that drag, and recover the ball's drag coefficient to ±20 %, and what frame rate and how many drops does that take? |
| **Assumptions and limits** | A point mass with a constant drag coefficient: F_D = ½ρC_D A v² opposite to the velocity, with C_D = 0.5 for a smooth sphere at Reynolds numbers 10³–10⁵ (the ball's 10⁴ at 5 m/s qualifies). No spin (a tossed ball with spin feels a Magnus force comparable to drag, which is why the toss is a prediction in this lesson and the drop is the measurement). Still air, constant ρ = 1.20 kg/m³ and g = 9.807 m/s², release from rest. The model fails when the ball spins, when the Reynolds number crosses 2 × 10⁵ (the drag crisis; not reached here), when the speed approaches the speed of sound (compressibility), and in a draught. The frame-quantisation uncertainty 0.41/fps limits what the camera can see: at 30 fps the 2 m drag excess (30 ms) is only 2.2 frame uncertainties. |
| **Analytical work** | Projectile kinematics from F = ma; the free-body diagram with drag; terminal speed; v(t) = v_t tanh(gt/v_t) and s(t) = (v_t²/g) ln cosh(gt/v_t) by separation of variables; the first-order drag correction t ≈ t₀(1 + gh/6v_t²); the horizontal decay x = ln(1 + kv₀t)/k; the wind triangle and its inverse; the round-trip inequality; the hand prediction of every number in the experiment table. |
| **Computational deliverable** | `starter/main.py`: `drag_parameter`, `terminal_speed`, `projectile_rhs`, `vacuum_state`, `vacuum_time_of_flight`, `vacuum_range`, `vertical_drop_state`, `vertical_drop_time`, `horizontal_decay_state`, `simulate_projectile`, `drag_range`, `optimal_angle`, `wind_triangle`, `heading_for_track`, `predict_fall_times`, `fit_terminal_speed`, `drag_coefficient_from_terminal_speed` (`frame_timing`, `reduce_timings`, `normalized_error`, the SYNTHETIC data generator and `main` are given). Tests prove: k scales as area over mass; drag equals weight at v_t and the acceleration there is zero; the drag force opposes and is parallel to the velocity; the vacuum closed forms are mutually consistent and symmetric in complementary angles; the vertical and horizontal closed forms satisfy their ODEs; the drop time inverts the distance law to 10⁻¹⁰, keeps its digits at v_t = 10⁶ m/s and matches the first-order expansion; the simulation reproduces the vacuum formulas to 10⁻⁸ without drag, the vertical closed form to 10⁻⁹ for a drop, the horizontal decay at early times, and an independent DOP853 solution to 10⁻⁶ m; drag shortens range, time, apex and impact speed, steepens the descent, and moves the apex before mid-flight; the range is monotone in k and converges with tolerance; the optimal angle is 45° without drag and lower with drag; the wind triangle is vector addition with the right symmetries and `heading_for_track` inverts it; the fit recovers the generating terminal speed from SYNTHETIC frame-quantised drops within 2.5 standard uncertainties and C_D within 25 %. |
| **Practical activity** | **physical experiment**: drop a ping-pong ball and a dense object from five heights, time the falls by frame counting, and fit the ball's terminal speed; toss the ball and compare the landing point with the drag prediction. See `projects/p04-drag-drop/`. A learner without the kit completes the lesson with the SYNTHETIC dataset that `main()` generates, labelled as such throughout. |
| **Equipment** | Level 1, reusing the `p00-first-experiment` kit: a ping-pong ball, a dense soft object (bean bag), a tape measure, a phone camera (30 fps minimum; 60 or 120 fps recommended), a phone stand, a soft mat. See `projects/p04-drag-drop/bom.yaml`. |
| **Safety** | See `projects/p04-drag-drop/safety.md`; summary box in part 4. |
| **Evidence** | `data/prediction.md` dated before the first drop; `data/raw/drops.csv` and `data/raw/tosses.csv`; the fitted v_t ± u and C_D with the fit's residuals; `projectile_drag.png`; the E_n table against the vacuum and the drag model; the derivations; the report. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: the drag equation and its two closed forms are derived with assumptions named; all public tests pass; a prediction with uncertainty exists dated before the measurement; the fitted C_D carries an uncertainty and every E_n > 2 has a named, tested cause. |
| **References** | `@openstax-university-physics-vol1` Ch. 3 "Motion Along a Straight Line": 3.4 Motion with Constant Acceleration, 3.5 Free Fall (required-free: the kinematic equations); Ch. 4 "Motion in Two and Three Dimensions": 4.3 Projectile Motion, 4.5 Relative Motion in One and Two Dimensions (required-free: the vacuum trajectory, velocity addition); Ch. 5 "Newton's Laws of Motion": 5.3 Newton's Second Law, 5.5 Newton's Third Law, 5.7 Drawing Free-Body Diagrams (required-free: the law and the diagram discipline); Ch. 6 "Applications of Newton's Laws": 6.4 Drag Force and Terminal Speed (required-free: the quadratic drag law and v_t). `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (supplementary-free: events and dense output). |
| **Extensions** | ★ Add a Magnus force a = S (ω × v) with a spin parameter and show that backspin lengthens the toss. ★ Replace the constant C_D by a tabulated C_D(Re) for a sphere and see whether the fit improves. ★ Build the break-beam timer of the `p00-first-experiment` enhanced tier and measure the 2 m drop to 1 ms. |
| **Time** | unmeasured |

## 1. The engineering problem

Every trajectory you will compute in this course, a sounding rocket, a re-entry capsule, a satellite, a glider, a rotor blade element, is Newton's second law integrated in time. The method never changes: draw the body alone, draw every force on it, add them, divide by the mass, and you have the acceleration; integrate twice and you have the path. What changes is the list of forces. This lesson has the shortest list that is still real, weight and air drag, and its subject is the point at which the "no drag" textbook model stops being good enough and the camera on your phone can tell.

The two objects are a ping-pong ball and a bean bag. Dropped from 2.0 m in a vacuum both would land after 0.639 s; in air the bean bag (terminal speed 18 m/s) lands after 0.645 s and the ping-pong ball (terminal speed 8.4 m/s) after 0.669 s, 30 ms or 4.7 % late. Tossed at 6 m/s and 45° the vacuum range is 3.67 m; the bean bag reaches 92 % of it and the ping-pong ball 73 %, 2.67 m, and the best launch angle for the ball is no longer 45° but 42.6°. These are the predictions. The experiment asks whether a frame-counted timing with a Type B uncertainty of 0.41/fps (13.6 ms at 30 fps, 6.8 ms at 60 fps, 3.4 ms at 120 fps) can resolve the 5–30 ms drag excess, and whether a least-squares fit of the terminal speed to five heights recovers the ball's drag coefficient. The answer, worked out in part 4 before any drop, is that 30 fps is marginal (the vacuum model is rejected at 2 m with E_n = 2.2 only if the Type A scatter is small), 60 fps gives a clean rejection above 1.2 m and C_D to about ±20 % with ten drops per height, and 120 fps gives ±10 %.

Success criteria, written first: the simulation reproduces every closed form in its limit to better than 10⁻⁸; the predictions above are in `data/prediction.md` with uncertainties before the first drop; the fitted v_t has an uncertainty and lies within 2 u of the nominal 8.4 m/s or the discrepancy has a named cause; the toss lands within the predicted range band or the reason (spin) is argued from the data.

## 2. Understand and derive

### 2.1 Newton's second law and the free-body diagram

Newton's second law for a particle is **ΣF = m a**: the vector sum of all forces on the body equals its mass times its acceleration (OpenStax 5.3). Two disciplines make it usable. First, the **free-body diagram** (5.7): draw the body alone, and draw only forces that *act on it*, each one exerted by something you can name. Second, the third law (5.5): every force in the diagram has a partner acting on the other body, so the drag that the air exerts on the ball is matched by the ball pushing the air, which is where the energy goes.

For a ball in flight there are two forces. The weight **W = m g**, downward, exerted by the Earth. The drag **F_D = ½ ρ C_D A v²**, exerted by the air, directed against the velocity, with ρ the air density, A = πd²/4 the frontal area and C_D a dimensionless coefficient (6.4). Nothing else: no "force of motion", no force from the hand after release. Units: ρ [kg/m³] × A [m²] × v² [m²/s²] = kg·m/s² = N.

```
          z
          ^              F_D = ½ ρ C_D A v²   (opposite to v)
          |                   ↖
          |                    ●  ball, mass m
          |                    | 
          |                    v  W = m g
          +----------------------> x
```

Divide by m and write the acceleration in components with **k = ρC_D A/(2m)** (units 1/m):

**ẍ = −k |v| ẋ,  z̈ = −g − k |v| ż,  |v| = √(ẋ² + ż²).**   (1)

The two equations are coupled through |v|: the horizontal drag depends on the vertical speed and vice versa. That coupling is why the general case has no closed form and why three special cases, which decouple, carry the verification.

For the ping-pong ball: A = 1.257 × 10⁻³ m², k = 1.20 × 0.5 × 1.257 × 10⁻³/(2 × 0.0027) = **0.1396 m⁻¹**. Its reciprocal, 7.2 m, is the distance scale of the drag: over a few metres of flight the drag has a large effect, over a few centimetres none. For the bean bag (100 g, 80 mm, C_D taken as 1.0 for a lumpy shape, which is a guess to be checked) k = 0.0302 m⁻¹, 33 m.

### 2.2 No drag: constant acceleration

With k = 0, (1) is ẍ = 0, z̈ = −g, and integrates twice (3.4, 4.3) from the origin with speed v₀ at angle θ:

x = v₀ cos θ · t,  z = v₀ sin θ · t − ½ g t²,  ẋ = v₀ cos θ,  ż = v₀ sin θ − g t.

Three assumptions were used and each has a name: no drag; constant g (true to 10⁻⁶ over 2 m); a point mass (no rotation affects the path). Setting z = 0 gives the **time of flight T = 2v₀ sin θ/g**, the **range R = v₀² sin 2θ/g** and the **apex H = (v₀ sin θ)²/2g**. The range is maximal at 45° and symmetric about it (30° and 60° give the same R); the path is a parabola symmetric about its apex; the impact speed equals the launch speed. For v₀ = 6 m/s at 45°: T = 0.865 s, R = 3.671 m, H = 0.918 m. For a drop from rest: t₀ = √(2h/g), 0.6387 s at 2 m.

### 2.3 Terminal speed

Drop the ball from rest. Only the vertical equation matters, and with s the distance fallen (positive down, ṡ = v):

**dv/dt = g − k v².**   (2)

The acceleration vanishes when kv² = g, at the **terminal speed v_t = √(g/k) = √(2mg/(ρC_D A))** (6.4). A body released from rest approaches it from below and never exceeds it; a body thrown down faster than v_t slows to it. Ping-pong ball: v_t = √(9.807/0.1396) = **8.38 m/s**; at that speed the drag ½ × 1.2 × 0.5 × 1.257 × 10⁻³ × 8.38² = 0.0265 N equals the weight 0.0027 × 9.807 = 0.0265 N, as it must. Bean bag: 18.0 m/s. A tennis ball (58 g, 67 mm): 23.2 m/s. The impact speed from 2 m in vacuum is √(2gh) = 6.26 m/s, 75 % of the ping-pong ball's v_t: the ball is deep in the drag regime. For the bean bag it is 35 % and for the tennis ball 27 %.

### 2.4 The vertical closed form

Write (2) as dv/dt = g(1 − v²/v_t²) and separate: ∫dv/(1 − v²/v_t²) = ∫g dt, so v_t artanh(v/v_t) = gt and

**v(t) = v_t tanh(gt/v_t),  s(t) = (v_t²/g) ln cosh(gt/v_t).**   (3)

Check it against (2): tanh′ = 1 − tanh², so dv/dt = g(1 − v²/v_t²). Check the limits: for gt ≪ v_t, tanh u ≈ u and ln cosh u ≈ u²/2, so v ≈ gt and s ≈ ½gt²; for gt ≫ v_t, v → v_t and s ≈ v_t t − (v_t²/g) ln 2, a straight line at terminal speed that started "late" by (v_t/g) ln 2. Inverting s = h for the fall time:

**t(h) = (v_t/g) arccosh(exp(gh/v_t²)).**   (4)

At 2 m: gh/v_t² = 0.2793, exp = 1.3221, arccosh = 0.7826, v_t/g = 0.8546 s, **t = 0.6688 s**, against 0.6387 s in vacuum: +30.1 ms. The ball hits at 5.48 m/s (65 % of v_t) instead of 6.26 m/s. One more expansion is worth having. With ln cosh u ≈ u²/2 − u⁴/12 and a little algebra (Problem 2),

**t ≈ t₀ (1 + gh/(6v_t²)),  t₀ = √(2h/g),**   (5)

which gives 0.6684 s at 2 m (0.06 % below the exact 0.6688 s, because gh/v_t² = 0.28 is not that small) and says in one line what the experiment is about: the drag excess grows linearly with h and inversely with v_t². The table that the code prints, and that `data/prediction.md` must contain:

| h (m) | t₀ vacuum (s) | t ping-pong (s) | t bean bag (s) | excess ping-pong (ms) | (%) |
|---|---|---|---|---|---|
| 0.6 | 0.3498 | 0.3547 | 0.3509 | 4.9 | 1.4 |
| 0.9 | 0.4284 | 0.4375 | 0.4304 | 9.0 | 2.1 |
| 1.2 | 0.4947 | 0.5086 | 0.4977 | 13.9 | 2.8 |
| 1.5 | 0.5531 | 0.5726 | 0.5573 | 19.5 | 3.5 |
| 1.8 | 0.6059 | 0.6316 | 0.6114 | 25.7 | 4.2 |
| 2.0 | 0.6387 | 0.6688 | 0.6451 | 30.1 | 4.7 |

A numerical detail: for weak drag gh/v_t² is tiny and arccosh(exp a) loses its digits to cancellation; the code evaluates (4) as (v_t/g) log1p(x + √(x(x + 2))) with x = expm1(a), and a test checks the vacuum limit at v_t = 10⁶ m/s to 10⁻⁹.

### 2.5 The horizontal closed form

Horizontal motion with drag and no gravity (a puck on ice in air, or the x-motion of a projectile before gravity has built up any vertical speed): dv/dt = −kv², so 1/v = 1/v₀ + kt and

**v(t) = v₀/(1 + kv₀t),  x(t) = ln(1 + kv₀t)/k.**   (6)

The speed halves after a distance ln 2/k, 5.0 m for the ping-pong ball, whatever v₀ is: quadratic drag has a length scale, not a time scale. The distance grows without bound but only logarithmically; (6) is the second oracle for the simulation (a horizontal launch from a great height, checked over the first 50 ms before gravity matters).

### 2.6 The coupled case, and what drag does to a toss

With both gravity and drag, (1) must be integrated numerically, but two qualitative results follow from it without computing. First, ẋ = v₀ cos θ − ∫k|v|ẋ dt decreases monotonically (the integrand is positive while ẋ > 0), so the horizontal speed at impact is below the launch value and the descent is steeper than the ascent: the path is not a parabola and its apex comes before mid-flight. Second, the kinetic energy decreases throughout (drag does negative work: the subject of `P04.02`), so the impact speed is below the launch speed. For the 6 m/s toss the code gives: range 2.670 m (73 %), time of flight 0.789 s (less than the vacuum 0.865 s: a lower apex, 0.765 m against 0.918 m, is reached sooner), apex at 0.380 s of 0.789, ẋ from 4.24 m/s at launch to 3.37 at the apex and 2.72 at impact, impact speed 4.51 m/s (44 % of the kinetic energy gone), impact angle 52.8° below the horizontal against 45° at launch. And the best angle moves down to 42.6° (range 2.679 m, 0.3 % more than at 45°: the optimum is flat); at 15 m/s it would be 37.1°, because the drag matters more at higher speed.

### 2.7 Relative motion: the wind triangle

Velocities add as vectors (4.5): the velocity of an aircraft over the ground is its velocity relative to the air plus the velocity of the air relative to the ground, **v_ground = v_air + v_wind**. In aviation v_air has magnitude TAS (true airspeed) and direction the *heading*; v_ground has magnitude the ground speed and direction the *track*; the wind is quoted by the direction it blows *from*. With east and north components, v_air = TAS (sin ψ, cos ψ) for heading ψ measured clockwise from north, and a wind of speed W from direction φ is W(−sin φ, −cos φ). An aircraft at 60 m/s on heading 090° in 15 m/s from 030°: v_air = (60, 0), v_wind = (−7.5, −12.99), v_ground = (52.5, −12.99), ground speed 54.1 m/s, track 103.9°: it drifts 13.9° to the right. The inverse problem, which heading makes good a wanted track, is the **wind correction angle** WCA with sin WCA = (W/TAS) sin(φ − χ) for track χ (the crosswind component must be cancelled by the across-track component of the airspeed), heading = χ + WCA, and ground speed TAS cos WCA − W cos(φ − χ): here heading 077.5° at 51.1 m/s. It has no solution when the crosswind exceeds the airspeed. And a round trip into and back from a wind along the track always takes longer than in still air, because d/(V − W) + d/(V + W) = 2dV/(V² − W²) > 2d/V (Problem 5): 59.3 min against 55.6 min for 100 km each way at 60 m/s in 15 m/s of wind.

### 2.8 Where the model fails

- **Spin.** A ball thrown by hand spins; the Magnus force on a spinning sphere is of the order of the drag at a few hundred rpm and bends the path up (backspin) or down (topspin). It is why the toss is a prediction to be compared with, not the measurement that fixes C_D; the drop, released from rest with two fingers, spins negligibly.
- **C_D is not a constant.** For a smooth sphere C_D ≈ 0.47–0.5 over Re = vd/ν from about 10³ to 2 × 10⁵; the ping-pong ball at 2–8 m/s is at Re = 0.5–2.2 × 10⁴, inside that range. Above 2 × 10⁵ the boundary layer turns turbulent and C_D drops to 0.1 (the drag crisis, which the dimples of a golf ball exploit); nothing in this lesson gets there.
- **Air motion.** A draught of 0.5 m/s is 6 % of v_t and shifts the landing point of a toss by centimetres. Close the window.
- **Density.** ρ = 1.20 kg/m³ is 20 °C at sea level; at 1500 m altitude it is about 1.06 and the drag parameter is 12 % smaller. Use your own ρ from p and T.
- **The point-mass assumption** is exact for the translation of a rigid body: the centre of mass obeys ΣF = ma whatever the body does about it (the subject of `P04.04`).

## 3. Build the computational model

Open `starter/main.py`. `frame_timing`, `reduce_timings` and `normalized_error` are the functions of `P00.15`; `make_synthetic_drop_data` (SYNTHETIC: frame-quantised times generated from the drag model, for testing the pipeline and for learners without the kit) and `main` are given.

### Step 1: the drag parameter and the right-hand side

```python
def drag_parameter(mass, diameter, cd=0.5, rho=RHO_AIR):
    return rho * cd * (math.pi * diameter ** 2 / 4) / (2 * mass)

def projectile_rhs(t, y, k, g=G0):
    vx, vz = y[2], y[3]
    speed = math.hypot(vx, vz)
    return np.array([vx, vz, -k * speed * vx, -g - k * speed * vz])
```

`terminal_speed` returns √(g/k), and `inf` when k = 0.

### Step 2: the closed forms

`vacuum_state(t, v0, angle_deg)` returns [x, z, vx, vz] for scalar or array t (use `np.stack(..., axis=-1)`); `vacuum_time_of_flight` and `vacuum_range` are one line each. `vertical_drop_state(t, v_t)` is (3), with the `v_t = inf` branch and an overflow-safe ln cosh (for u > 20 use u − ln 2 + log1p(e^{−2u})); `vertical_drop_time(h, v_t)` is (4) in the `expm1`/`log1p` form, with the branch a > 300 where arccosh(eᵃ) = a + ln 2; `horizontal_decay_state(t, v0, k)` is (6).

### Step 3: the simulation

```python
def simulate_projectile(v0, angle_deg, k, g=G0, z0=0.0, rtol=1e-10, n_samples=400, max_time=60.0):
    a = math.radians(angle_deg)
    y0 = np.array([0.0, z0, v0 * math.cos(a), v0 * math.sin(a)])
    def hit_ground(t, y, *args): return y[1]
    hit_ground.terminal, hit_ground.direction = True, -1
    sol = solve_ivp(projectile_rhs, (0.0, max_time), y0, args=(k, g), rtol=rtol, atol=1e-12,
                    events=hit_ground, dense_output=True)
    ...
```

The event is armed only for a downward crossing (`direction = -1`) so a launch from z₀ = 0 does not stop at t = 0. Sample the dense output at `n_samples` times up to the event, read the impact state from `sol.y_events`, and return the `Trajectory` dataclass with time of flight, range, impact speed and apex height. `drag_range` calls it at rtol 10⁻⁸ (the range converges to 10⁻⁷ m there, which is all a range study needs); `optimal_angle` maximises the range with `minimize_scalar(..., bounds=(5, 85), method="bounded")`.

### Step 4: relative motion

`wind_triangle` adds the two vectors in east/north components and returns the ground speed and the track from `atan2(east, north)` reduced to [0, 360); `heading_for_track` is the closed form of 2.7 and raises `ValueError` when |sin WCA| > 1 or the ground speed is not positive.

### Step 5: the experiment pipeline

`predict_fall_times(heights, k)` returns the vacuum and drag times. `fit_terminal_speed(heights, times, u_times)` fits (4) to the mean times by weighted least squares, with `scipy.optimize.least_squares` on the parameter w = 1/v_t² (so that "no drag" is the finite value w = 0 and the bound w ≥ 0 is natural), and returns v_t with the standard uncertainty from the Jacobian at the solution scaled by the reduced chi-square; `drag_coefficient_from_terminal_speed` inverts v_t² = 2mg/(ρC_D A).

Running `main()` prints the terminal speeds, the prediction table of 2.4, the solve_ivp drop against (4) (agreement 7 × 10⁻¹⁴ s), the toss numbers of 2.6 with the range against angle table, the convergence of the range with tolerance (2.67022430, 2.67023276, 2.67023286, 2.67023286 m at rtol 10⁻⁴ to 10⁻¹⁰), the wind triangle of 2.7, and the SYNTHETIC drop dataset through the pipeline, then writes `projectile_drag.png` (the toss with and without drag for both objects; range against angle; fall time against height with the SYNTHETIC points).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-01-kinematics-newtons-laws`. The tests:

1. **Drag parameter**: proportional to C_D, to d², inversely to m.
2. **Terminal speed**: kv_t² = g; between 7 and 10 m/s for the ping-pong ball; infinite for k = 0; the right-hand side gives zero acceleration at v_t.
3. **Right-hand side**: free fall at k = 0; the drag acceleration opposes and is parallel to v with magnitude k|v|².
4. **Vacuum closed forms**: z(T) = 0, x(T) = R, mirrored impact velocity, constant ẋ, R(45°) = v₀²/g, R(30°) = R(60°).
5. **Vertical closed form**: ṡ = v and v̇ = g(1 − v²/v_t²) by central differences; v → v_t; the vacuum branch.
6. **Drop time**: inverts s(t) = h to 10⁻¹⁰ at four heights; always longer than √(2h/g); the vacuum limit at v_t = ∞ and at 10⁶ m/s to 10⁻⁹; the first-order expansion (5) at v_t = 40 m/s to 10⁻⁵; a 10⁶ m fall takes 10⁶/v_t s without overflow.
7. **Horizontal closed form**: satisfies its ODE; x(2 s) exact; the k = 0 branch.
8. **Simulation, no drag**: time of flight, range and impact speed to 10⁻⁸, apex to 10⁻⁴, the sampled path against `vacuum_state` to 10⁻⁷ m.
9. **Simulation, drop**: time to 10⁻⁹ against (4); path and speed against (3) to 10⁻⁸; impact below v_t.
10. **Simulation, horizontal launch from 1000 m**: x(t) against (6) over the first 50 ms to 0.2 %; the whole path against an independent `solve_ivp(DOP853, rtol=1e-12)` call to 10⁻⁶ m.
11. **What drag does**: shorter range, time, apex, impact speed; steeper descent; apex before mid-flight; range monotone in k over five values.
12. **Convergence**: the range at rtol 10⁻⁷ is within 10⁻⁶ m of the rtol 10⁻¹⁰ value and closer than the 10⁻⁴ one.
13. **Optimal angle**: 45.00 ± 0.05° and v₀²/g without drag; between 35° and 45° with drag and a local maximum; lower still at 4k.
14. **Wind triangle**: no wind, headwind, tailwind, pure crosswind, mirror symmetry; `heading_for_track` inverted by `wind_triangle` at three cases to 10⁻⁸°; the straight-headwind case; `ValueError` for an impossible track.
15. **Pipeline on SYNTHETIC data**: predictions ordered with a growing excess; the fit recovers v_t within 2.5 u (u < 0.6 m/s) and C_D within 25 % for three seeds at 120 fps with 8 drops; the fit on exact times is exact to 10⁻⁶; a dense object fits a v_t more than twice the ball's; E_n against the drag model below 2.5 at every height for four seeds at 60 fps, and the vacuum model rejected at 1.8 m in at least three of them.
16. `main` writes the figure.

### Compare

`scipy.integrate.solve_ivp` is both the tool and, in a second configuration (DOP853 at rtol 10⁻¹², no events), the oracle; the two agree to 10⁻⁶ m on a 14 s trajectory, and the closed forms agree with RK45 to 10⁻¹³ s on the drop. The closed forms are the real comparison: (3), (4) and (6) are exact, cost nothing, and bracket the coupled problem from both sides. What they cannot give is the toss, where the coupling through |v| is the whole effect (the range drops 27 % while a naive "drag along x only" estimate from (6) with the vacuum flight time, ln(1 + 0.1396 × 4.24 × 0.865)/0.1396 = 2.96 m, overestimates it by 11 %). The lesson's habit: integrate numerically, but never without a closed-form limit to land on.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box and the full procedure are in `projects/p04-drag-drop/` (`bom.yaml`, `safety.md`, `procedure.md`, `report-template.md`, `data/README.md`). The kit is the `p00-first-experiment` kit plus a ping-pong ball. A learner without a camera completes the lesson with the SYNTHETIC dataset that `main()` generates: it is generated from the drag model and labelled SYNTHETIC in the code, the printout and the figure, and must be labelled so in anything you write.

> **Safety (summary; the authoritative box is `projects/p04-drag-drop/safety.md`).** Hazards: a falling or tossed object (eye, face, phone); a bounce toward something fragile; a trip over the tape measure; reaching a release height. Energy limits: object mass ≤ 0.1 kg, drop height ≤ 2.0 m, toss speed by hand only, ≤ 8 m/s, never toward a person, a window or a screen; phone battery only. Controls: soft mat at the impact point; a 1 m exclusion circle for drops and a 4 m × 2 m clear lane for tosses; safety glasses; phone on a stand 1.5 m to the side of the lane; tape fixed flat and removed after marking. Supervision: none needed for an adult; children with an adult. Location: indoors, hard floor, clear of furniture edges and glass, door closed, windows closed (draughts). Stop criteria: anyone in the lane or circle; a bounce toward a person or glass; a draught you can feel; any floor damage; experimenter tired or rushed. Not permitted: chairs or ladders as release points; metal, glass or heavy objects; any launcher, elastic, spring or compressed air (hand tosses only); drops or tosses outdoors in wind, or near pets.

The setup for the drops is that of `P00.15`, with the ping-pong ball and the bean bag released by the same method (two fingers opened from below) at the same five heights, 0.6 to 1.8 m, at least eight drops each. The toss lane:

```
   release                            landing zone (mat, 1 m long)
   ○ → v0 ≈ 6 m/s, 45°                     ▓▓▓▓▓▓▓▓
   |                                       |
   |<------------- 3.5 m lane ------------>|      phone on a stand 1.5 m to the side,
   |                                       |      lens at 0.5 m height, whole lane in frame
   ═══════════════════════════════════════════  floor, tape marks every 0.5 m
```

Procedure in outline (details in `procedure.md`):

1. **Predict first.** Run `main()` with your own ball's mass (kitchen scale, 0.1 g) and diameter (callipers or a ruler across the ball, 0.5 mm), your air density from today's pressure and temperature, and C_D = 0.5 ± 0.1. Write into `data/prediction.md`, dated: the table of 2.4 for your heights for both objects with uncertainties (propagate u_h = 5 mm, u_g = 0.03 m/s², u_CD = 0.1 through `vertical_drop_time` by finite differences, as `P00.14` taught), the expected E_n of each height against the vacuum model at your frame rate, the toss range band for 5–7 m/s at 40–50°, and the frame rate and number of drops you chose and why (the design study of Problem 4).
2. **Calibrate** the camera's frame rate against a stopwatch display (10 s, count frames) as in `P00.15`, and write your frame rule for release and impact.
3. **Drop.** Alternate the two objects at each height; eight drops each; record `height_m, object, drop, release_frame, impact_frame, fps_nominal, note` in `data/raw/drops.csv`. Never edit it.
4. **Toss.** Ten hand tosses of the ping-pong ball from a marked release point at about 0.5 m height, aiming at 45° and 6 m/s; from the video read the release frame, the frame and position of the apex (against the 0.5 m floor marks), the landing frame and position. Record in `data/raw/tosses.csv`. The launch speed and angle come from the first three frames after release (displacement over 2/fps); their uncertainty is large and is part of the result.
5. **Reduce.** For each height and object: `reduce_timings` with the frame term; E_n against vacuum and against the drag model with the predicted uncertainty; `fit_terminal_speed` on the ping-pong means with the combined uncertainties; C_D from the fit with its uncertainty (u_CD/C_D = 2u_v/v_t); the bean bag's fit as a control (expect a v_t too large to be determined: the fit should return a large value with a large uncertainty, and that *is* the result). For each toss: `simulate_projectile` at the measured launch speed and angle, k from your fitted v_t, and E_n on the range with the uncertainty from the launch speed propagated through the simulation.

What you should find (from the model; your numbers will differ): the ping-pong ball's times above the vacuum curve by 5–26 ms, growing with h, and the bean bag's within a frame of it; the vacuum model rejected (E_n > 2) from 1.2 m up at 60 fps (2.05 at 1.2 m, 3.8 at 1.8 m, frame term alone) and never at 30 fps; v_t = 8.4 ± 0.8 m/s at 60 fps with ten drops per height (±0.4 at 120 fps), hence C_D = 0.50 ± 0.10 (±0.05); the toss landing 2.4–3.0 m out at 6 m/s, systematically beyond the prediction if the ball has backspin and short of it with topspin.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| k, ping-pong ball (m⁻¹) | 0.1396 | 0.1396 | 0 | — |
| v_t (m/s) | 8.38 | 8.38 | 0 | — |
| drop from 2 m, vacuum (s) | 0.6387 | 0.6387 | 0 | — |
| drop from 2 m, drag, (4) (s) | 0.6688 | 0.668758 (solve_ivp) | 7 × 10⁻¹⁴ | solver tolerance 10⁻¹² |
| drop from 2 m, first order (5) (s) | 0.6684 | — | 0.06 % of exact | gh/v_t² = 0.28 is not small |
| toss, vacuum range (m) | 3.671 | 3.671 | 0 | — |
| toss, drag range (m) | — | 2.670 | — | no closed form |
| drag along x only, (6) with vacuum T | 2.96 | 2.670 | +11 % | the coupling through \|v\| |
| best angle with drag | < 45° | 42.6° | — | as argued |
| range convergence (rtol 10⁻⁷ vs 10⁻¹⁰) | — | 1 × 10⁻⁷ m | — | — |
| ground speed, 60 m/s in 15 m/s from 030° | 54.1 m/s, track 103.9° | 54.08, 103.9 | 0 | — |
| SYNTHETIC fit at 60 fps, 5 drops | v_t = 8.38 | 8.00 ± 0.71 | 0.5 u | quantisation at a random phase |
| SYNTHETIC fit at 30 fps, 5 drops | 8.38 | 8.0 ± 1.4 | 0.3 u | frame term doubled |

The physics rows all close; the two SYNTHETIC rows show what the uncertainty budget predicts, that the fit is limited by the frame term and that 30 fps gives ±17 % on v_t, ±35 % on C_D, with five drops. With your own data the table has three more rows and they are the lesson: E_n of each object against each model, the fitted v_t against the nominal, and the toss range against the prediction. The diagnosis patterns of `P00.15` §5 apply unchanged (a constant offset at every height is the release; an error that shrinks with h is the height reference), with two new ones: the bean bag's times above the vacuum curve by more than a frame mean your C_D guess of 1.0 for it was low (it does not matter: the bag is the control for the release and the frame rule, not a measurement of its own drag); and a toss that lands consistently long or short with the measured v₀ and θ means spin, which you can confirm by marking the ball with a pen and watching it rotate in the video.

What to improve, in order of effect: the frame rate (each doubling halves the frame term and roughly halves u_v); the number of drops only after the frame term is below the Type A scatter; the launch speed measurement of the toss (the three-frame displacement at 60 fps is 20 cm with a 1 cm reading error, 5 %: a higher frame rate helps here too). What not to improve: the solver tolerance (10⁻⁷ m against a 1 cm landing-mark reading), or C_D beyond the ±10 % that the sphere literature itself quotes for this Reynolds range.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4), (5) and (6) and of the wind correction angle; `data/prediction.md` (dated before the first drop), `data/raw/drops.csv`, `data/raw/tosses.csv`; `projectile_drag.png` with your data replacing the SYNTHETIC points; the E_n table for both objects against both models; the fitted v_t ± u and C_D ± u; and a report following `projects/p04-drag-drop/report-template.md`.

### Adversarial self-checks

- **What if the mass doubles at the same diameter?** k halves, v_t grows by √2 to 11.9 m/s, the 2 m drag excess falls from 30 ms to about 15 ms (by (5), roughly as 1/v_t²), and at 30 fps the vacuum model can no longer be rejected at any height: the experiment's sensitivity is set by k, not by h alone.
- **What if the launch speed goes to zero?** The toss becomes a drop and the simulation must reproduce (4): it does, to 10⁻¹³ s. What if k → 0? Every drag result must reduce to the vacuum formula, and the tests hold the code to that at v_t = 10⁶ m/s.
- **Which assumption fails first in the real world?** For the toss, no-spin: a hand toss has spin, and the Magnus force is of the order of the drag. For the drop, the constant C_D: the ball crosses Re = 0.5–2.2 × 10⁴ during a fall, over which C_D for a sphere varies by about ±5 %. Still air and constant ρ fail only with an open window.
- **What measurement could falsify the drag model?** A fall-time excess that does not grow with height (it must, by (5)); a bean bag that falls as late as the ping-pong ball; or a fitted v_t that depends on which heights are included.
- **Does the result survive a change of solver or tolerance?** RK45 at rtol 10⁻⁴ to 10⁻¹⁰ and DOP853 at 10⁻¹² agree on the range to 10⁻⁵ m and on the drop time to 10⁻¹³ s; a result that depended on the tolerance at the 1 cm level would be a bug.
- **Which decision changes if the measurement uncertainty doubles?** At 60 fps with ten drops, u_v = 0.8 m/s gives C_D ± 20 %, and the experiment's conclusion ("C_D consistent with 0.5") stands; doubled to 1.6 m/s it would be ±40 %, consistent with anything from 0.3 to 0.8, and the right decision would be to switch to 120 fps rather than take more drops, because the frame term does not average down.

## Check yourself

1. Write Newton's second law for a ball in air, draw its free-body diagram, and name the body that exerts each force.
2. Define the drag parameter k, give its units, and say what 1/k means physically.
3. What is the terminal speed of a body with k = 0.14 m⁻¹? What fraction of it does a vacuum drop from 2 m reach?
4. State the vertical closed form v(t) and s(t) and show that they reduce to free fall for gt ≪ v_t.
5. By how much, to first order, does drag lengthen a fall from height h? Does the correction grow or shrink with v_t?
6. Why is the descent of a projectile with drag steeper than its ascent, and why does the best launch angle fall below 45°?
7. An aircraft flies at 50 m/s TAS heading 000° in 20 m/s of wind from 090°. What are the ground speed and track, and what heading makes good a track of 000°?
8. At 30 fps, can a drop from 1.2 m distinguish the ping-pong ball's drag from free fall? Give the numbers that decide it.
