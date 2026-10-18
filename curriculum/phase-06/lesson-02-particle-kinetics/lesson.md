# P06.02  Particle Kinetics

> Write F = m a along the path and across it. Across it, the lift has to hold up the weight *and* bend the path, and that sum is the load factor.

| | |
|---|---|
| **Learning outcomes** | 1. Draw the free-body diagram of an aircraft as a particle and write Newton's second law in path (tangential–normal) coordinates, m V̇ = T − D − W sin γ and m V γ̇ = L − W cos γ. 2. Derive the load factor of a pull-up, n = 1 + V²/(gR) at the bottom and cos γ + V²/(gR) around a loop, of a level turn (n = 1/cos φ), and the seat and wing loads that follow. 3. Derive the altitude lost and the time taken in a constant-n, constant-speed dive recovery in closed form, and integrate the point-mass flight path with `solve_ivp` to check it. 4. Recover the load factor from integrated positions alone, and show that energy is conserved when thrust equals drag. 5. Decide how much a finite g-onset rate costs in a dive recovery and which limit (pilot, structure, lift) governs. |
| **Prerequisites** | `P06.01` (path coordinates, a_n = V²/ρ), `P04.01` (Newton's laws, free-body diagrams), `P04.02` (work and energy), `P05.01` (forces and free-body diagrams), `P02.02` (second-order ODEs), `P03.06` (adaptive integrators and events in `solve_ivp`), `P03.04` (numerical differentiation). Skill check: you can integrate dh/dγ = f(γ) by substitution and you can stop `solve_ivp` with an event. |
| **Engineering question** | A fighter (12 000 kg, illustrative) pulls up from a 60° dive at 250 m/s. How much altitude does it lose reaching level flight at 7 g, how much more if the g comes on at a finite rate, what does the 100 kg pilot press on the seat and what must the wing lift, and what radius gives 9 g at the bottom of a pull-up? The answers: 510.1 m in 4.33 s with the g applied at once; 608.8 m at a 6 g/s onset and 1053 m at 1 g/s; 6.86 kN on the seat and 824 kN of lift at 7 g (8.83 kN and 1059 kN at 9 g); and R = 797 m. |
| **Assumptions and limits** | The aircraft is a particle in a vertical plane (no rotation dynamics: the attitude follows the commanded load factor instantly, apart from the stated onset ramp). Uniform gravity, flat Earth, still air. Either constant speed (throttle and speed brakes hold V) or thrust equal to drag. Lift is whatever the load factor needs: the wing's maximum lift and the structure's limit are checked separately, not modelled. The pilot is a rigid mass; human g tolerance (which depends on g-onset and duration) is outside the model. |
| **Analytical work** | The path-coordinate equations of motion from the free-body diagram; n = L/W and its relation to an accelerometer; n at the bottom and top of a loop; the level-turn bank, radius and rate; the closed-form altitude loss (V²/g) ln((n − cos γ₀)/(n − 1)) and time of a recovery; its small-angle limit; a delay estimate of the cost of g-onset; the energy integral with thrust equal to drag. |
| **Computational deliverable** | `starter/main.py`: `loop_load_factor`, `radius_for_load_factor`, `level_turn`, `apparent_weight`, `pullout_altitude_loss`, `pullout_time`, `flight_path_rhs`, `simulate_pullout`, `kinematic_load_factor` (`ramp_load_factor`, `recovery_altitude_loss` and `main` are given). Tests prove: n = 1 + V²/(gR) and its inverse; bottom minus top is 2 g at constant speed; the level turn satisfies both force balances; seat load n m g; the closed forms equal quadrature and their limits; straight and level is an equilibrium; `solve_ivp` reproduces the closed-form altitude loss and time to 10⁻⁸; the load factor recovered from positions alone matches the command to 2 × 10⁻³; a loop flown with n = cos γ + V²/(gR) is a circle to 10⁻⁵ m; V²/2 + gh is conserved to 10⁻⁸ with thrust equal to drag; slower g-onset always costs altitude. |
| **Practical activity** | **computational experiment**: the dive-recovery altitude loss by closed form and by simulation, a sweep of g-onset rate, constant speed against thrust-equals-drag, a constant-speed loop, and a sensitivity study of the recovery to speed and load factor. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pullout.png`, the free-body diagram and the derivations of the path equations, the loop load factor and the closed-form recovery, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the free-body diagram is drawn in path coordinates and the normal equation derived from it; all public tests pass; the closed-form and simulated altitude losses are compared with numbers; the note states the effect of g-onset and the limit that governs. |
| **References** | `@mit-ocw-16-07-dynamics` L6 "Intrinsic coordinates": Equations of Motion in Intrinsic Coordinates (required-free: ΣF_t = m V̇ and ΣF_n = m V²/ρ). `@openstax-university-physics-vol1` Ch. 6 "Applications of Newton's Laws": 6.1 Solving Problems with Newton's Laws; 6.3 Centripetal Force, subsection Banked Curves (required-free: the banked curve is the coordinated level turn with the road replaced by lift). `@scipy-userguide` "Integration (scipy.integrate)" > "Ordinary differential equations (solve_ivp)" (required-free: events and dense output). |
| **Extensions** | ★ Add a lift limit: L ≤ ½ρV²S C_L,max with illustrative S = 27.9 m² and C_L,max = 1.4, and find the speed below which the 7 g recovery becomes lift-limited (the corner speed). ★ Replace the onset ramp by a first-order lag n̈ + … (a pilot-plus-airframe response with a time constant) and fit its time constant to the ramp results. ★ Add drag D = ½ρV²S C_D(n) with an induced-drag term and redo the energy audit for the thrust-equals-drag case. |
| **Time** | unmeasured |

## 1. The engineering problem

Every manoeuvre that bends an aircraft's flight path does it with lift. In straight and level flight the lift balances the weight; in a pull-up or a turn the wing must make more lift than the weight, and the ratio, the **load factor** n = L/W, is the number that sizes the wing, limits the manoeuvre and presses the pilot into the seat. A fighter's structure is typically designed for around 9 g; a transport's for 2.5 g. A pilot can tolerate high g only if it comes on slowly enough for the body to respond, which is why the rate of g-onset matters as much as its peak.

The question here is a dive recovery, the manoeuvre that automatic ground-collision avoidance systems have to perform at the last safe moment. A fighter is diving at 60° and 250 m/s. It pulls to 7 g. How much height does it need? The particle model is enough to answer: Newton's second law written along and across the flight path, with the lift chosen to give the commanded load factor. It gives a closed form for an idealised recovery, and a simulation for the real one with a finite g-onset rate. The same equations answer the classic pull-up question: what radius puts 9 g on the pilot at the bottom of a loop?

Success criteria, written first: the closed-form altitude loss and time match `solve_ivp` to its tolerance; the load factor recovered from the simulated positions alone matches the commanded one; a loop flown with the derived n(γ) closes into a circle; with thrust equal to drag the specific energy V²/2 + gh is constant; and the g-onset sweep is monotonic, with the instant-g closed form as its limit.

## 2. Understand and derive

### 2.1 The free-body diagram in path coordinates

Treat the aircraft as a particle of mass m moving in a vertical plane at speed V, with flight-path angle γ (positive climbing). Four forces act: the weight W = mg down; the lift L perpendicular to the velocity, in the plane; the thrust T and the drag D along the velocity (assume the thrust line is close to the flight path).

![Free-body diagram of the aircraft as a particle of mass m on a flight path climbing at angle γ above the horizontal: the velocity V, thrust T and drag D lie along the path, the lift L is perpendicular to it, and the weight W = m g is vertical, at angle γ from the path normal.](figures/path-coordinates-fbd.svg)
*Figure 1. Path coordinates: T and D act along the velocity, L across it, and only the weight has components on both axes, −W sin γ along the path and −W cos γ across it, as in (1).*

Use path coordinates (`P06.01`): e_t along the velocity and e_n perpendicular to it, toward the centre of curvature of a pull-up. The acceleration is a = V̇ e_t + V γ̇ e_n, because the velocity direction turns at γ̇ and V γ̇ = V²/ρ is the normal acceleration (16.07 L6, "Equations of Motion in Intrinsic Coordinates"). The weight has components −W sin γ along e_t and −W cos γ along e_n. Newton's second law along each axis:

**m V̇ = T − D − W sin γ,   m V γ̇ = L − W cos γ.**   (1)

Two more equations carry the position: ẋ = V cos γ and ḣ = V sin γ. Units check: every term of (1) is a force.

### 2.2 Load factor

Divide the normal equation by W and write n = L/W:

**γ̇ = g (n − cos γ)/V,  or  n = cos γ + V γ̇/g = cos γ + V²/(g ρ).**   (2)

The load factor is what an accelerometer along the aircraft's normal axis reads, in g: an accelerometer measures the specific force f = a − g_vec, and its component along e_n is V γ̇ + g cos γ = n g. In level flight n = 1; in a free-fall parabola n = 0.

**The bottom of a pull-up.** At γ = 0 on a path of radius R: **n = 1 + V²/(gR)**, and **R = V²/(g(n − 1))**. At 250 m/s, 9 g needs R = 250²/(9.80665 × 8) = **796.7 m**. The pilot (100 kg with gear, illustrative) presses on the seat with n m g = **8.83 kN**, nine times his weight; the wing must lift 9 × 117.68 kN = **1059 kN**, 530 kN per wing.

**Around a constant-speed loop.** On a circle of radius R at constant V, ρ = R everywhere and (2) gives n = cos γ + V²/(gR): largest at the bottom, smallest at the top (γ = π), and the difference is always exactly 2 g. At 250 m/s: R = 1000 m gives 7.37 g at the bottom and 5.37 g at the top; R = 4000 m gives 2.59 g and 0.59 g. A real loop loses speed going up, so the top is gentler still.

**A sustained level turn.** Bank the lift by φ. Vertical balance L cos φ = W gives **n = 1/cos φ**; the horizontal component provides the centripetal force, L sin φ = m V²/R, so **R = V²/(g√(n² − 1))** and the turn rate is ω = V/R = g√(n² − 1)/V (OpenStax 6.3, Banked Curves: the same balance with the road's normal force in place of the lift). At 9 g and 250 m/s: φ = 83.62°, R = 712.5 m, ω = 20.10 °/s. The level turn is tighter than the pull-up at the same n (712 m against 797 m) because in the pull-up one g of the lift is spent holding up the weight, while in the turn the weight is balanced by the vertical component and √(n² − 1) g is left for turning.

### 2.3 The dive recovery in closed form

Hold the speed constant (the throttle and speed brake absorb the weight component, so V̇ = 0) and hold n constant from the start. The normal equation (2) and ḣ = V sin γ give, dividing one by the other,

dh/dγ = V² sin γ / (g (n − cos γ)).

The numerator is the derivative of n − cos γ, so the integral is a logarithm. From γ₀ < 0 to level flight:

**Δh = (V²/g) ln((n − cos γ₀)/(n − 1)).**   (3)

The time follows from dt = V dγ/(g(n − cos γ)); with the half-angle substitution u = tan(γ/2),

**t = (V/g) [F(0) − F(γ₀)],  F(γ) = (2/√(n² − 1)) atan(√((n + 1)/(n − 1)) tan(γ/2)).**   (4)

For γ₀ = −60°, V = 250 m/s and n = 7: **Δh = 510.13 m, t = 4.327 s.** For a 30° dive, 140.7 m.

Checks. Dimensions: V²/g is a length and V/g a time. Small angle: for |γ₀| ≪ 1, n − cos γ₀ ≈ (n − 1) + γ₀²/2, so Δh ≈ (V²/g) γ₀²/(2(n − 1)) = R γ₀²/2 with R = V²/(g(n − 1)), the drop of a circular arc of radius R through an angle γ₀: right, because near level flight the path is almost a circle. Doubling V quadruples Δh. As n → 1 the loss is infinite: with lift equal to weight the dive never ends.

### 2.4 g-onset and the delay estimate

A real pull does not jump from 1 g to 7 g. If the load factor ramps at ṅ g/s, it takes (n − 1)/ṅ seconds to reach n, and during the ramp the aircraft, still diving steeply, has an average extra load factor of only about half the final increment. A rough estimate treats the ramp as a delay of half its length, during which the aircraft descends at V sin|γ₀|:

extra Δh ≈ V sin|γ₀| (n − 1)/(2ṅ).   (5)

At 6 g/s this is 250 × 0.866 × 0.5 = 108.3 m. The simulation (part 4) gives 98.7 m: the estimate is the right size and slightly high, because the flight path starts turning during the ramp.

### 2.5 Energy with thrust equal to drag

If instead T = D, the tangential equation (1) is V̇ = −g sin γ, and multiplying by V: d(V²/2)/dt = −g V sin γ = −g ḣ. So **V²/2 + g h is constant**: the lift does no work (it is perpendicular to the velocity), and the aircraft trades height for speed in the dive. That is an independent check on the integrator, and it makes the recovery longer, because the aircraft speeds up while it is still pointing down.

### 2.6 Where the model stops

- **Rotation is ignored.** A real aircraft must first pitch to the angle of attack that gives n; that takes a fraction of a second and is part of what the onset ramp stands for.
- **Lift is unlimited here.** At 250 m/s near 3 km (ρ ≈ 0.909 kg/m³) and an illustrative wing area of 27.9 m², 9 g needs C_L ≈ 1.34, close to a fighter wing's maximum. Below a certain speed (the corner speed) the wing, not the pilot or the structure, sets the limit (Extension 1).
- **Speed control is idealised.** Constant speed in a 60° dive requires a large speed brake; thrust equal to drag is the other extreme.
- **The human.** g tolerance depends on onset rate, duration and the anti-g suit; the model only reports n(t).

## 3. Build the computational model

Open `starter/main.py`. `ramp_load_factor`, `recovery_altitude_loss` and `main` are given.

### Step 1: the free-body results

`loop_load_factor(V, R, gamma)` is (2) on a circle; `radius_for_load_factor(V, n, gamma)` its inverse (inf when n = cos γ); `level_turn(V, n)` returns (φ, R, ω) and refuses n < 1; `apparent_weight(m, n)` is n m g.

### Step 2: the closed-form recovery

`pullout_altitude_loss(V, n, gamma0)` is (3) and `pullout_time` is (4). Both refuse n ≤ 1.

### Step 3: the equations of motion

`flight_path_rhs(t, y, n_cmd, speed_hold, excess_thrust)` returns the derivatives of y = [x, h, V, γ] from (1) and (2). With `speed_hold` V̇ = 0; otherwise V̇ = g((T − D)/W − sin γ).

```python
x, h, V, gam = y
n = n_cmd(t, y)
dV = 0.0 if speed_hold else G0 * (excess_thrust - np.sin(gam))
return [V * np.cos(gam), V * np.sin(gam), dV, G0 * (n - np.cos(gam)) / V]
```

### Step 4: integrate to level flight

`simulate_pullout(V0, gamma0, n_cmd, ...)` calls `solve_ivp` with dense output and an event on γ = 0 rising (terminal unless `stop_level=False`). rtol 10⁻¹⁰, atol 10⁻⁹ by default.

### Step 5: the load factor from positions alone

`kinematic_load_factor(t, x, h)` differentiates the positions twice with `np.gradient` (second order), finds γ = atan2(ḣ, ẋ), forms the specific force f = a + g e_h and returns f·e_n/g with e_n = (−sin γ, cos γ). It never sees the commanded n: this is how you would check a flight test from a tracking record.

Running `main()` prints the pull-up and turn numbers, the closed-form and simulated recovery, the load factor recovered from positions, the g-onset table, the thrust-equals-drag case with its energy drift, and the constant-speed loop; it writes `pullout.png` (recovery paths at 1 g/s, 4 g/s and instant g; commanded and recovered n(t); altitude loss against onset rate with the closed form).

### Verify

Run `python -m pytest curriculum/phase-06/lesson-02-particle-kinetics`. The tests:

1. **Pull-up**: n = 1 + V²/(gR) and R(n) are inverses; doubling V and quadrupling R leaves n unchanged.
2. **Loop**: top n = V²/(gR) − 1; bottom minus top = 2 g; an infinite radius gives n = cos γ.
3. **Level turn**: L cos φ = W, L sin φ = mV²/R, ω = V/R; n = 1 is straight; n < 1 refused.
4. **Seat and wing loads** scale with n and vanish at n = 0.
5. **Closed forms** against `scipy.integrate.quad` to 10⁻¹⁰; zero dive angle loses nothing; the small-angle limit Rγ₀²/2; more g, less loss; Δh ∝ V².
6. **Equations of motion**: straight and level is an equilibrium; an unloaded climb curves down at g cos γ/V and loses speed at g sin γ.
7. **Simulation**: stops at level flight; altitude loss and time equal (3) and (4) to 10⁻⁸.
8. **Positions alone**: n recovered from x(t), h(t) equals the commanded ramp to 2 × 10⁻³ (away from the ramp's kink and the ends, where `np.gradient` is less accurate).
9. **Loop**: n = cos γ + V²/(gR) flies a circle of radius R to 10⁻⁵ m and turns the velocity through 2π.
10. **Energy**: V²/2 + gh constant to 10⁻⁸ with thrust equal to drag; excess thrust adds energy at g V (T − D)/W.
11. **g-onset**: altitude loss falls monotonically as the onset rate rises, toward the closed form; thrust equal to drag needs more height than constant speed.
12. `main` writes the figure.

### Compare

Two independent routes meet. The closed form (3) and (4) gives 510.130 m in 4.3272 s; `solve_ivp` (RK45, rtol 10⁻¹⁰) gives the same to 2.3 × 10⁻⁸ m and 1.9 × 10⁻¹⁰ s. The load factor recovered from the integrated positions alone, on 2001 points, averages 6.999999 with a largest error of 3.6 × 10⁻⁴ (at the ends, where `np.gradient` is one-sided). With thrust equal to drag, the specific energy drifts by 1.0 × 10⁻¹² relative. A constant-speed loop at R = 1500 m flown with the law n = cos γ + V²/(gR) stays on its circle to 3.7 × 10⁻⁸ m, with n between 3.249 at the top and 5.249 at the bottom. Each check would catch a different mistake: a sign error in cos γ breaks the loop, a wrong V in γ̇ breaks the closed form, and a drag term in the wrong equation breaks the energy.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are simulations of the particle model with illustrative numbers, not flight measurements. Write each prediction down before running it.

1. **Closed form against simulation.** 60° dive, 250 m/s, 7 g applied at once. Predict the altitude loss from (3). *Finding:* 510.13 m in 4.327 s, both routes.
2. **g-onset.** The same recovery with the load factor ramped from 1 g at 1, 2, 4, 6 and 10 g/s. Predict the extra loss at 6 g/s with (5). *Finding:*

| onset (g/s) | altitude loss (m) | extra over instant g (m) |
|---|---|---|
| 1 | 1053.1 | 543.0 |
| 2 | 797.5 | 287.3 |
| 4 | 657.2 | 147.0 |
| 6 | 608.8 | 98.7 |
| 10 | 569.7 | 59.5 |
| instant | 510.1 | 0 |

The extra loss is close to inversely proportional to the onset rate, as (5) says (108.3 m predicted at 6 g/s). A slow pull at 1 g/s more than doubles the height needed.
3. **Speed control.** Repeat at 6 g/s with thrust equal to drag instead of constant speed. Predict whether the loss grows and the exit speed. *Finding:* 669.6 m (61 m more), leaving the dive at 275.0 m/s; energy drift 10⁻¹².
4. **Pull-up and loop.** Tabulate n at the bottom and top of constant-speed loops of 800, 1000, 2000 and 4000 m at 250 m/s, and fly the 1500 m loop with the law (2). *Finding:* 8.97/6.97, 7.37/5.37, 4.19/2.19, 2.59/0.59 g; the 1500 m loop closes to 4 × 10⁻⁸ m.
5. **Sensitivity.** Recompute the 6 g/s recovery for V = 225, 250 and 275 m/s and for 6 g instead of 7 g. Predict which matters more, a 10 % speed error or one g less. *Finding:* 501.9, 608.8 and 725.9 m (−18 % and +19 %, close to the V² scaling of the instant-g part); 6 g at 250 m/s needs 688.8 m (+13 %).

Raw data format: CSV tables `onset_gps, dh_m, extra_m` (2); `mode, dh_m, V_exit_mps, energy_drift` (3); `R_m, n_bottom, n_top` (4); `V_mps, n, dh_m` (5); plus `pullout.png`.

## 5. Validate and improve

| quantity | closed form (§2) | simulation | difference | cause |
|---|---|---|---|---|
| altitude loss, 60°, 7 g, instant (m) | 510.130 (3) | 510.130 | 2.3 × 10⁻⁸ | solver tolerance |
| recovery time (s) | 4.3272 (4) | 4.3272 | 1.9 × 10⁻¹⁰ | solver tolerance |
| load factor from positions | 7 (commanded) | 6.999999 mean, 3.6 × 10⁻⁴ max error | 3.6 × 10⁻⁴ | `np.gradient` at the ends |
| loop radius, n = cos γ + V²/gR (m) | 1500 | 1500 ± 3.7 × 10⁻⁸ | 4 × 10⁻⁸ | solver tolerance |
| specific energy, T = D | constant | drift 1.0 × 10⁻¹² | 10⁻¹² | solver tolerance |
| extra loss at 6 g/s (m) | 108.3 (5) | 98.7 | 10 % | (5) ignores turning during the ramp |

Every exact result is reproduced to the solver's tolerance; the only real discrepancy is the delay estimate (5), which is an approximation and is 10 % high. The engineering answer: **from a 60° dive at 250 m/s a 7 g recovery needs 510 m if the g arrives at once, 609 m at a realistic 6 g/s onset and over 1 km at 1 g/s; the pilot presses 6.86 kN on the seat and the wing lifts 824 kN at 7 g; a 9 g pull-up at 250 m/s has a radius of 797 m.** The onset rate and the speed are the levers; the structure's 9 g limit is not reached in this recovery.

What to improve: a lift limit (C_L,max) and drag, which decide whether 7 g is available at all at lower speed; and a pitch-dynamics model in place of the ramp. What not to improve: the solver tolerance.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the free-body diagram and the derivations of (1)–(5); `pullout.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the speed doubles?** Every length in (3) scales as V²: the instant-g loss quadruples to 2041 m. The delay term (5) only doubles, so the onset rate matters relatively less at high speed.
- **What if n → 1?** (3) diverges: with lift equal to weight the path never levels. Near 1 g the recovery is hypersensitive to n.
- **Which assumption fails first in the real world?** Unlimited lift at low speed, then the idealised speed control. At 150 m/s the same 7 g would need C_L ≈ 2.9 with the illustrative wing: not available, so the aircraft cannot pull 7 g there at all.
- **What measurement could falsify the model?** A flight-test record of n(t) from the aircraft's accelerometer and of x, h from tracking: the n recovered from tracking must match the accelerometer, and the altitude loss must match the simulation driven by the recorded n(t).
- **Does the result survive a change of solver?** The closed form does not depend on any solver; the simulation agrees with it to 10⁻⁸ m, and the energy check holds to 10⁻¹².
- **Which decision changes if the uncertainty doubles?** A ground-collision avoidance trigger set at the predicted 609 m plus a 10 % speed margin (about 120 m) would need about 240 m of margin if the speed were known only to 20 %: the trigger would fire earlier and more often.

## Check yourself

1. Write Newton's second law along and across the flight path for an aircraft in the vertical plane.
2. Define the load factor. What does an accelerometer along the aircraft's normal axis read in a free-fall parabola?
3. What radius gives 5 g at the bottom of a pull-up at 200 m/s?
4. Why is the load factor at the bottom of a constant-speed loop exactly 2 g more than at the top?
5. At the same load factor and speed, which is tighter, a level turn or a pull-up, and why?
6. How does the altitude lost in a dive recovery scale with speed?
7. Why does a slower g-onset cost altitude, and roughly how much?
8. With thrust equal to drag, which quantity is conserved, and why does the lift not change it?
