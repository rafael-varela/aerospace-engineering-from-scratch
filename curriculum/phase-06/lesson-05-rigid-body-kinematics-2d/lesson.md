# P06.05  Rigid Body Kinematics in 2D

> A rigid body has one angular velocity. Every point on it moves as if the body were turning about one point, and the relative-velocity equation finds that motion from any two of its points.

| | |
|---|---|
| **Learning outcomes** | 1. State and use the rigid-body relations v_B = v_A + ω × r_B/A and a_B = a_A + α × r_B/A − ω² r_B/A, and find instant centres with Kennedy's theorem. 2. Write the loop-closure equations of a four-bar linkage and solve them for position by Newton's method and, independently, by intersecting two circles. 3. Derive the velocity and acceleration of the coupler and the rocker as linear systems with the loop's Jacobian, and implement and test them. 4. Check the analytic rates against finite differences of the position solution (second-order convergence) and against instant centres. 5. Classify a linkage by Grashof's criterion, and decide whether a control-surface linkage's gearing, symmetry and transmission angle are acceptable. |
| **Prerequisites** | `P06.01` (relative motion, turning unit vectors), `P04.04` (angular velocity), `P01.07` (cross products), `P01.13` (2 × 2 linear systems), `P03.02` (root finding), `P03.10` (Newton's method for nonlinear systems), `P03.04` (finite differences), `P05.02` (equilibrium of rigid bodies, virtual work). Skill check: you can solve a 2 × 2 linear system by hand and write a two-variable Newton step. |
| **Engineering question** | A servo arm (20 mm) drives a control-surface horn (25 mm) through a 120 mm pushrod; the servo shaft is 120 mm from the hinge line (an illustrative small-UAV installation). Over ±45° of servo travel, how far does the surface move, how linear and how symmetric is the gearing, how close does the linkage come to toggle, and what servo torque does a 0.3 N·m hinge moment need? The answers: −34.19° to +34.68° (0.49° asymmetric), gearing 0.800 at neutral falling to 0.669 at the ends (best straight line 0.780, worst departure 0.83°), transmission angle never below 52.9°, and 0.240 N·m at neutral. |
| **Assumptions and limits** | Rigid links, ideal pin joints with no clearance (backlash) and no friction, planar motion. The servo angle is prescribed (kinematics only; no inertia or compliance). The linkage stays on one assembly branch. The torque relation uses virtual work for a massless, frictionless linkage. |
| **Analytical work** | Rigid-body velocity and acceleration relations; instant centres and Kennedy's theorem; the loop-closure equations; Newton's method and its Jacobian, whose determinant r₃r₄ sin(θ₃ − θ₄) vanishes at toggle; the circle-intersection closed form and its two branches; velocity and acceleration as linear systems; velocity ratios from instant centres; the transmission angle by the law of cosines; Grashof's criterion; mechanical advantage by virtual work. |
| **Computational deliverable** | `starter/main.py`: `grashof`, `loop_residual`, `loop_jacobian`, `solve_position_newton`, `position_circle_intersection`, `joint_positions`, `velocity_analysis`, `acceleration_analysis`, `coupler_point`, `instant_centres`, `velocity_ratio_ic`, `transmission_angle` (`_cross2`, `_intersect`, `sweep_positions`, `central_rates`, `observed_order` and `main` are given). Tests prove: six Grashof classes; Newton closes the loop to 10⁻¹³ with link lengths preserved, agrees with the closed form on both branches and converges quadratically; the Jacobian against finite differences and singular at toggle; assembly failures are reported; the velocity solution satisfies v_B = v_A + ω₃ × r_BA and matches instant centres to 10⁻¹⁰, with Kennedy's collinearities; a parallelogram's coupler translates; rates are linear in the crank rate; the transmission angle against the law of cosines; finite differences of positions converge to the analytic ω and α at order 2.0 ± 0.1; a coupler point's velocity and acceleration match finite differences. |
| **Practical activity** | **computational experiment**: a ±45° servo sweep of the control linkage (deflection, gearing, transmission angle, servo torque); Newton's convergence; finite differences against analytic rates; and a horn-length design sweep across Grashof classes. |
| **Equipment** | None. A computer with the environment from `P00.01`. A hobby servo, pushrod and horn are a good way to see the linkage move, but nothing is measured with them. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `four_bar.png`, the derivations of the loop closure, the velocity and acceleration systems and the instant-centre ratio, the experiment tables, and a half-page note on the control linkage's acceptability. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the velocity system is derived from v_B = v_A + ω × r (not quoted); all public tests pass; the finite-difference convergence orders and the instant-centre agreement are reported; the note judges gearing, symmetry and transmission angle with numbers. |
| **References** | `@mit-ocw-16-07-dynamics` L21 "2D rigid body dynamics": Kinematics of Two-Dimensional Rigid Body Motion; Instantaneous Center of Rotation (required-free: v_P = v_O + ω × r and the instantaneous centre). `@openstax-university-physics-vol1` 10.3 Relating Angular and Translational Quantities (required-free: v = ωr and the tangential and centripetal accelerations of a point on a rotating body). `@scipy-userguide` "Optimization (scipy.optimize)" > "Root finding" > "Sets of equations" and `@ocw-16-90` Unit 1 "1.7.4 Newton-Raphson Implement Implicit Methods on Nonlinear Problems" (supplementary-free: Newton's method for systems). |
| **Extensions** | ★ Add pin clearance: let each joint move within a 0.05 mm circle and find the worst-case surface deflection error (the dead band) over the servo travel. ★ Find the servo neutral angle and horn length that make the gearing most linear over ±45° (minimise the worst departure from a straight line). ★ Trace a coupler-point curve of a crank-rocker and find a point whose path is nearly straight over part of the cycle (a Watt or Chebyshev straight-line linkage). |
| **Time** | unmeasured |

## 1. The engineering problem

Aircraft are full of linkages. Every control surface on a light aircraft or a UAV is driven by one: a control stick or a servo arm turns, a pushrod or a cable carries the motion, and a horn on the surface turns it about its hinge. Landing gear retracts through linkages; flaps follow tracks driven by linkages; engine cowl flaps, doors and thrust reversers open by linkages. Most of them are four-bars: four rigid links joined by four pins, one of them the fixed structure.

What a designer needs from a control linkage is simple to state and not obvious to get: the surface must move as far as the specification asks, its motion should be close to proportional to the servo's (so that the flight controller's gains mean the same thing everywhere), it should move the same amount up and down, and the linkage must never approach toggle, the position where the pushrod lines up with the horn, the mechanical advantage collapses and a small error in the servo makes the surface jump or lock. The servo's torque must also be enough for the hinge moment the airflow puts on the surface.

All of that is kinematics of rigid bodies in a plane. This lesson builds the tools: the relative-velocity and relative-acceleration equations, instant centres, and the loop-closure equations of a four-bar solved by Newton's method. The independent check is finite differences: the analytic angular velocities and accelerations must equal the time derivatives of the position solution, to second order in the step. Success criteria, written first: Newton's solution closes the loop to 10⁻¹³ and agrees with the circle-intersection closed form; the velocity analysis agrees with instant centres to 10⁻¹⁰; finite differences converge to the analytic ω₄ and α₄ at second order; the servo sweep gives the deflection, gearing and transmission angle with which to judge the design.

## 2. Understand and derive

### 2.1 Velocity and acceleration of a rigid body in a plane

A rigid body in a plane has one angular velocity ω (a vector along e_z) and one angular acceleration α. For two points A and B on it, the vector r_B/A = r_B − r_A has a fixed length and only turns, at ω. Its derivative is therefore ω × r_B/A (the derivative of a turning vector, as in `P06.01`), and

**v_B = v_A + ω × r_B/A,**   (1)

**a_B = a_A + α × r_B/A − ω² r_B/A.**   (2)

The last term is the centripetal acceleration of B about A; the middle one is the tangential acceleration (OpenStax 10.3 has the special case of a fixed axis: v = ωr, a_t = αr, a_c = ω²r). In two dimensions ω × r = ω(−r_y, r_x): a quarter turn of r scaled by ω (16.07 L21, "Kinematics of Two-Dimensional Rigid Body Motion").

**Instant centre.** If ω ≠ 0 there is one point of the body's plane with zero velocity at this instant, the **instant centre** I: from (1), v_P = ω × r_P/I for every point P (16.07 L21, "Instantaneous Center of Rotation"). Every point's velocity is perpendicular to the line from I, and proportional to its distance. If you know the directions of the velocities of two points, I is where the perpendiculars to them meet.

For several bodies, the instant centre I_jk of bodies j and k is the point where they have the same velocity. **Kennedy's theorem**: for any three bodies, the three instant centres I_jk, I_kl, I_jl lie on one straight line. (Proof sketch: if I_jl were off the line through I_jk and I_kl, its velocities as a point of j and as a point of l, each perpendicular to its own line to the other centres, could not be equal.)

### 2.2 The four-bar and its loop closure

Number the links: 1 the ground (from the servo shaft O₂ to the hinge O₄, length r₁, along +x), 2 the crank (servo arm, O₂ to A, length r₂, angle θ₂), 3 the coupler (pushrod, A to B, r₃, θ₃), 4 the rocker (horn, O₄ to B, r₄, θ₄). Angles are counterclockwise from +x.

```
                 A ●───────── r3 (pushrod) ─────────● B
                  /  θ3                               \
         r2 (arm)/                                     \ r4 (horn)
                / θ2                                θ4  \
     O2 ●──────────────────── r1 (ground) ──────────────● O4
```

Going round the loop from O₂ to A to B to O₄ and back must close:

**r₂ e^{iθ₂} + r₃ e^{iθ₃} − r₄ e^{iθ₄} − r₁ = 0,**   (3)

two real equations (real and imaginary parts) in the two unknowns θ₃, θ₄ for a given θ₂.

**Newton's method** (`P03.10`). With f(θ₃, θ₄) the residual of (3), the Jacobian is

**J = [[−r₃ sin θ₃, r₄ sin θ₄], [r₃ cos θ₃, −r₄ cos θ₄]],  det J = r₃ r₄ sin(θ₃ − θ₄),**   (4)

and the iteration is [θ₃, θ₄] ← [θ₃, θ₄] − J⁻¹f. It converges quadratically from a good guess. The determinant vanishes when the coupler and the rocker are in line (θ₃ = θ₄ or θ₃ = θ₄ + π): a **toggle** position, where Newton fails and, physically, the crank cannot drive the rocker.

**An independent closed form: two circles.** B lies on the circle of radius r₃ about A and on the circle of radius r₄ about O₄. With d the distance from A to O₄, the common chord's foot is a = (r₃² − r₄² + d²)/(2d) along A→O₄ and its half-length h = √(r₃² − a²). The two intersections, to the left and to the right of A→O₄, are the linkage's two **assembly branches** (open and crossed); if r₃² < a² the circles do not meet and the linkage cannot be assembled at that θ₂. Newton converges to whichever branch its guess is near; sweeping θ₂ and starting each solve from the last answer (continuation) keeps it on one branch.

**Grashof's criterion.** Let s and l be the shortest and longest links and p, q the other two. If s + l < p + q, the shortest link can turn fully relative to its neighbours: with the shortest as the crank, a **crank-rocker** (the servo could turn all the way round while the horn rocks); as the ground, a **double-crank**; as the coupler, a **double-rocker**. If s + l = p + q, a **change-point** linkage that can fold flat (a parallelogram is one); if s + l > p + q, no link turns fully (a **triple-rocker**). The lesson's servo linkage, 20 + 120 < 25 + 120, is a crank-rocker.

### 2.3 Velocity and acceleration of the four-bar

Differentiate (3) in time: each term r e^{iθ} gives r i θ̇ e^{iθ}, which is exactly ω × r of (1) written in complex form. With ω₂ = θ̇₂ known:

**J [ω₃, ω₄]ᵀ = [r₂ ω₂ sin θ₂, −r₂ ω₂ cos θ₂]ᵀ:**   (5)

the relative-velocity equation v_B = v_A + ω₃ × r_B/A = ω₄ × r_B/O₄, a *linear* system with the same Jacobian as Newton's method. Differentiate again: each term gives r(iα − ω²)e^{iθ}, which is (2):

**J [α₃, α₄]ᵀ = −[known terms in α₂, ω₂², ω₃², ω₄²].**   (6)

The rates are linear in ω₂ at fixed position (double the servo speed and every rate doubles), and both systems fail at toggle.

**Velocity ratios from instant centres.** The four pins are four of the six instant centres: I₁₂ = O₂, I₂₃ = A, I₃₄ = B, I₁₄ = O₄. Kennedy's theorem places I₁₃ (the coupler's instant centre relative to the ground) at the intersection of lines O₂A and O₄B, and I₂₄ at the intersection of the ground line O₂O₄ with the coupler line AB. Since I₂₄ is a point of both 2 and 4 with one velocity, ω₂(I₂₄ − O₂) = ω₄(I₂₄ − O₄) and

**ω₄/ω₂ = x₂₄/(x₂₄ − r₁).**   (7)

The gearing of the linkage can be read from a drawing: it is the ratio of I₂₄'s distances from the two ground pivots.

**Transmission angle.** The acute angle μ between the coupler and the rocker. The force the pushrod can exert is along the pushrod; only its component perpendicular to the horn (proportional to sin μ) turns the horn. By the law of cosines across the diagonal A–O₄, cos γ = (r₃² + r₄² − r₁² − r₂² + 2r₁r₂ cos θ₂)/(2r₃r₄), with μ = min(γ, π − γ). A rule of thumb keeps μ above about 40°.

**Torque.** For an ideal massless linkage, virtual work gives T₂ω₂ = M₄ω₄: the servo torque needed for a hinge moment M₄ is **T₂ = M₄ (ω₄/ω₂)**. A gearing below 1 trades surface travel for servo torque.

### 2.4 The control linkage by hand at neutral

At neutral the servo arm is vertical, θ₂ = 90°, A = (0, 20) mm, O₄ = (120, 0) mm. The circle intersection gives θ₃ = 2.388° and θ₄ = 90.239°: the horn is almost vertical and the pushrod almost horizontal, so μ = 87.85°. With both arms nearly perpendicular to the pushrod, the pushrod moves almost in pure translation: v_A = ω₂r₂ ≈ v_B = ω₄r₄, so ω₄/ω₂ ≈ r₂/r₄ = 0.800 (the code gives 0.7999). As the servo turns away from neutral, both arms tilt and their effective lever arms (r cos of the tilt) shrink at different rates: the gearing falls and the motion becomes nonlinear.

### 2.5 Where the model stops

- **Backlash and compliance.** Real pins have clearance and real pushrods bend; a few hundredths of a millimetre at a 25 mm horn is a dead band of about 0.1°, which flight controllers notice as limit cycles (Extension 1).
- **Toggle.** Near det J = 0 small errors in link length or servo angle produce large changes of θ₄; the model is exact but the real linkage is not.
- **Dynamics.** The servo's own torque limit, the surface's inertia and the aerodynamic hinge moment change with deflection; this lesson is kinematics.

## 3. Build the computational model

Open `starter/main.py`. `_cross2`, `_intersect`, `sweep_positions`, `central_rates`, `observed_order` and `main` are given.

### Step 1: classification and position

`grashof(r1, r2, r3, r4)` returns one of the six names. `loop_residual` and `loop_jacobian` are (3) and (4). `solve_position_newton(th2, links, guess)` iterates until |f| < 10⁻¹³ r₁ and returns the iteration count; it raises if it does not converge. `position_circle_intersection(th2, links, branch)` is the closed form, raising `ValueError` when the circles do not meet. `joint_positions` returns O₂, A, B, O₄.

### Step 2: velocity and acceleration

`velocity_analysis` solves (5) and `acceleration_analysis` solves (6):

```python
rhs = np.array([r2 * w2 * np.sin(th2), -r2 * w2 * np.cos(th2)])
w3, w4 = np.linalg.solve(loop_jacobian(th3, th4, links), rhs)
```

`coupler_point` applies (1) and (2) from A to a point P on the coupler.

### Step 3: instant centres and the transmission angle

`instant_centres` returns the six centres (I₁₃ and I₂₄ by line intersection); `velocity_ratio_ic` returns ω₃/ω₂ and ω₄/ω₂ from them alone; `transmission_angle` folds |θ₄ − θ₃| into [0, π/2].

Running `main()` prints the neutral position, Newton's error sequence, the ±45° sweep (deflection, gearing, transmission angle, torque), the instant-centre check, the finite-difference table and the horn-length sweep, and writes `four_bar.png` (the linkage at −45°, 0 and +45°; horn rotation against servo travel with its best straight line and the transmission angle; the finite-difference error against step).

### Verify

Run `python -m pytest curriculum/phase-06/lesson-05-rigid-body-kinematics-2d`. The tests:

1. **Grashof**: crank-rocker, double-crank, double-rocker, rocker-crank, change-point, triple-rocker.
2. **Position**: Newton closes the loop to 10⁻¹³ at four crank angles in at most 8 iterations and preserves all link lengths; it agrees with the closed form on both branches; the branches differ.
3. **Newton's convergence** is quadratic: e_{k+1}/e_k² stays bounded as the error falls below 10⁻⁶.
4. **Jacobian** against finite differences; singular when coupler and rocker are in line.
5. **Assembly**: impossible positions raise; a triple-rocker sweeps over its allowed range.
6. **Velocity**: v_B from the coupler equals v_B from the rocker to 10⁻¹²; instant-centre ratios equal the loop-closure ratios to 10⁻¹⁰; I₁₃ and I₂₄ lie on their Kennedy lines; a coupler point moves about I₁₃.
7. **Parallelogram**: the coupler translates (ω₃ = 0, α₃ = 0) and ω₄ = ω₂.
8. **Linearity** of the rates in ω₂; a small crank step moves the rocker by ω₄/ω₂ times as much.
9. **Transmission angle** against the law of cosines.
10. **Finite differences** of positions along a crank motion with angular acceleration: ω₃, ω₄, α₃, α₄ converge at order 2.0 ± 0.1.
11. **Coupler point**: velocity and acceleration from (1) and (2) match finite differences of its position.
12. `main` writes the figure.

### Compare

Two independent checks. **Instant centres** at servo +30°: ω₃/ω₂ = 0.019009 and ω₄/ω₂ = 0.762392 from the loop closure; 0.019009 and 0.762392 from (7) and the I₁₃ construction, with I₁₃ = (516.1, −893.8) mm and I₂₄ = (−385.0, 0) mm. **Finite differences** of the Newton position solution, with the servo at θ₂ = 67.30°, ω₂ = 2.150 rad/s and α₂ = 3.0 rad/s² (analytic ω₄ = 1.65802151 rad/s, α₄ = 2.98704383 rad/s²):

| h (s) | \|ω₄ FD − ω₄\| (rad/s) | \|α₄ FD − α₄\| (rad/s²) |
|---|---|---|
| 0.0200 | 1.087 × 10⁻⁴ | 5.247 × 10⁻⁴ |
| 0.0100 | 2.717 × 10⁻⁵ | 1.311 × 10⁻⁴ |
| 0.0050 | 6.791 × 10⁻⁶ | 3.278 × 10⁻⁵ |
| 0.0025 | 1.698 × 10⁻⁶ | 8.195 × 10⁻⁶ |

Observed orders 2.00, 2.00, 2.00 (velocity) and 2.001, 2.000, 2.000 (acceleration). A sign error in one term of (6), such as the ω₄² r₄ term, would leave the velocity column converging and the acceleration column stuck at a constant error: the two columns test different equations.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are kinematic calculations for an illustrative linkage, not measurements. Write each prediction down before running it.

1. **Newton's convergence.** Start Newton at neutral from a deliberately poor guess (11.5°, 68.8°). Predict the number of correct digits after each iteration. *Finding:* errors 4.1 × 10⁻¹, 3.7 × 10⁻², 2.3 × 10⁻⁴, 3.9 × 10⁻⁸, 2.5 × 10⁻¹⁶ rad: the digits double each step, then round-off.
2. **The servo sweep.** ±45° of servo travel about neutral. Predict the surface travel from the gearing r₂/r₄ and whether it is symmetric. *Finding:* −34.19° to +34.68° (asymmetry 0.49°), against ±36° from 0.8 × 45°; gearing 0.7999 at neutral, falling to 0.669 at −45°; best straight line 0.7797 °/° with a worst departure of 0.83°; transmission angle 87.85° at neutral and 52.90° at −45° (its minimum); servo torque for a 0.3 N·m hinge moment 0.240 N·m at neutral, 0.201 and 0.210 N·m at −45° and +45°. The servo works hardest at neutral, where the gearing is highest.
3. **Two checks.** Instant centres and finite differences as in Compare.
4. **Horn length.** Repeat the sweep with horns of 15, 20, 25, 30 and 40 mm (servo arm 20 mm). Predict which one is linear. *Finding:*

| horn (mm) | Grashof class | deflection (°) | asymmetry (°) | minimum transmission angle (°) |
|---|---|---|---|---|
| 15 | rocker-crank | −67.81 … +76.91 | 9.10 | 7.51 |
| 20 | change-point | −45.00 … +45.00 | 0.00 | 45.00 |
| 25 | crank-rocker | −34.19 … +34.68 | 0.49 | 52.90 |
| 30 | crank-rocker | −27.48 … +28.64 | 1.16 | 57.26 |
| 40 | crank-rocker | −19.40 … +21.58 | 2.18 | 61.41 |

The 20 mm horn makes a parallelogram: perfectly linear (gearing exactly 1) and symmetric, but a change-point linkage that can fold into a crossed configuration. A 15 mm horn (shorter than the servo arm) is a rocker-crank, multiplies the motion by 1.33 and comes within 7.5° of toggle: unacceptable. Longer horns give less travel, more servo torque margin and better transmission, at the price of growing asymmetry.

Raw data format: CSV tables `iteration, error_rad` (1); `servo_deg, horn_deg, gearing, mu_deg, torque_Nm` (2); `h_s, err_w4, err_al4` (3); `horn_mm, class, defl_min_deg, defl_max_deg, asym_deg, mu_min_deg` (4); plus `four_bar.png`.

## 5. Validate and improve

| quantity | analytic (§2) | code | difference | cause |
|---|---|---|---|---|
| neutral horn angle (°) | 90.2388 (circle intersection) | 90.2388 (Newton, 4 iterations) | < 10⁻¹² | — |
| neutral gearing | ≈ r₂/r₄ = 0.8000 | 0.7999 | 1.3 × 10⁻⁴ | the arms are not exactly perpendicular to the pushrod |
| ω₄/ω₂ at +30° | 0.762392 (instant centres, (7)) | 0.762392 (loop closure, (5)) | < 10⁻¹⁰ | — |
| ω₄ at 67.30° (rad/s) | 1.65802151 (5) | FD, h = 2.5 ms: within 1.7 × 10⁻⁶ | O(h²) | finite differences, order 2.00 |
| α₄ at 67.30° (rad/s²) | 2.98704383 (6) | FD, h = 2.5 ms: within 8.2 × 10⁻⁶ | O(h²) | finite differences, order 2.00 |
| parallelogram (20 mm horn) | linear, gearing 1 | −45.00° … +45.00° | 0 | — |

The engineering answer: **over ±45° of servo travel the 20/120/25 mm linkage moves the surface −34.19° to +34.68°; the gearing is 0.800 at neutral and 0.669 at the ends, so a straight line is good to 0.83°; the transmission angle never falls below 52.9°, comfortably above 40°; and the servo needs 0.240 N·m for a 0.3 N·m hinge moment, most at neutral.** It is acceptable for a UAV control surface. If linearity mattered more than torque, a 20 mm horn (a parallelogram) would be exactly linear; if torque margin mattered more, a 30 mm horn.

What to improve: backlash (Extension 1), which sets the real dead band; the hinge moment's variation with deflection; and the neutral angle that minimises nonlinearity (Extension 2). What not to improve: the position solution, which closes to 10⁻¹³.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(7); `four_bar.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) judging the control linkage on travel, linearity, symmetry, transmission angle and servo torque.

### Adversarial self-checks

- **What if the servo travel doubles to ±90°?** The gearing keeps falling and the transmission angle reaches about 33° at −90° (a probe sweep to ±120° gives its minimum there): below the 40° rule. The linkage is a crank-rocker, so the servo could turn fully, but the horn would reverse at its extremes.
- **What if the horn shrinks toward the servo arm's length?** At r₄ = r₂ the linkage becomes a parallelogram (exactly linear); below it, the horn can turn fully and the transmission angle collapses (7.5° at 15 mm).
- **Which assumption fails first in the real world?** Ideal pins: clearance in four joints gives a dead band, and a rubbery pushrod adds compliance under the hinge moment.
- **What measurement could falsify the model?** Measure surface deflection against servo angle with a protractor at 5° steps: the curve must follow the sweep within the measuring error plus the dead band. A systematic asymmetry larger than 0.49° would point to a wrong neutral angle or link length.
- **Does the result survive a change of method?** Newton and the closed form agree to 10⁻¹²; instant centres and the loop closure agree to 10⁻¹⁰; finite differences converge to the analytic rates at order 2.
- **Which decision changes if the uncertainty doubles?** If the pushrod length were uncertain by ±2 mm instead of ±1 mm, the neutral horn angle would shift by up to about ±4.6° instead of ±2.3° (dθ₄ ≈ dr₃/r₄ at neutral): the trim would have to be set on the aircraft, not from the drawing.

## Check yourself

1. Write the velocity of point B of a rigid body in terms of the velocity of A and the angular velocity. Why does the relative term have no component along r_B/A?
2. What is an instant centre, and how do you find it from the velocity directions of two points?
3. State Kennedy's theorem.
4. Why does Newton's method for the four-bar fail at toggle, and what is happening physically?
5. Classify a linkage with ground 100, crank 30, coupler 110, rocker 60 (mm) by Grashof's criterion.
6. Why are the velocity and acceleration problems linear, when the position problem is not?
7. A parallelogram linkage has ω₂ = 3 rad/s. What are ω₃ and ω₄?
8. A linkage has a gearing ω₄/ω₂ = 0.6. What servo torque does a 0.5 N·m hinge moment need, ignoring friction?
