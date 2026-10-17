# P06.18  Multibody Simulation

> A multibody model is bodies, joints and force elements, and two honest ways to write its equations. If both ways agree and the energy books balance at every instant, you may start believing the loads.

| | |
|---|---|
| **Learning outcomes** | 1. Model a landing-gear drop test as a planar two-body system: airframe and unsprung mass joined by an oleo-pneumatic strut (polytropic gas spring, quadratic orifice damping, top-out stop), on a compression-only tyre, with lift. 2. Derive its equations of motion twice: by Lagrange's equations in two reduced coordinates, and by constrained Newton–Euler in six coordinates with four constraint equations and Lagrange multipliers. 3. Implement and test both, show that they agree, that the multipliers are the physical joint and guide forces, and that Baumgarte stabilisation removes constraint drift. 4. Close the energy books: show that the energy put in equals the energy stored plus the energy dissipated at every instant, and use the model to compute a drop test's load factor, stroke, strut efficiency and bearing side load. 5. Decide from a parameter study which strut parameter most changes the landing load. |
| **Prerequisites** | `P06.11` (Lagrange's equations, generalised forces by virtual work), `P06.05` and `P06.06` (planar rigid-body kinematics and kinetics), `P06.03` (work–energy for impacts), `P06.13` (SDOF dynamics), `P04.07` (ideal gas, polytropic compression), `P03.06` (adaptive and stiff integrators), `P01.16` (linear systems; a block matrix solve). Skill check: you can derive the generalised force for a coordinate from the virtual work of a set of forces. |
| **Engineering question** | One main gear of a light business jet carries 1540 kg (1500 kg sprung, 40 kg unsprung) and touches down at 3.05 m/s (10 ft/s, a typical limit sink rate) with lift equal to weight. Its strut has a 100 mm piston, a 22.6 mm orifice, 1.5 MPa precharge and 300 mm of stroke, raked 8° (illustrative values). Does it absorb the impact without bottoming, what is the gear load factor, how efficient is the strut, and what side load do its bearings carry? The answers: no bottoming (229 mm of stroke), peak tyre force 32.2 kN (load factor 2.13), strut efficiency 0.91, bearing side load 4.3 kN at the peak. |
| **Assumptions and limits** | Planar motion; a drop-test rig: the airframe is guided to move vertically without pitching; rigid strut with a frictionless sliding joint (no bearing friction, which real struts have and TR-1154 includes); gas compressed polytropically with one exponent; orifice flow quasi-steady with a constant discharge coefficient and the same orifice in both directions (no recoil valve); linear tyre with light damping and no tyre mass; frictionless ground (no wheel spin-up, so no drag load); constant lift. The numbers are illustrative, of the size of a light business jet's main gear; they are not data for any aircraft. |
| **Analytical work** | Kinematics of the axle in reduced coordinates; the kinetic energy and the constant mass matrix; generalised forces by virtual work, including an internal force along a sliding joint; the gas force and its stored energy in closed form; the orifice force; the hanging stroke before touchdown; the constraint equations, their Jacobian and the acceleration right-hand side γ; the augmented (saddle-point) system; Baumgarte's stabilised constraint equation and why explicit Runge–Kutta preserves linear constraints; the energy balance; strut efficiency. |
| **Computational deliverable** | `starter/main.py`: the force elements (`gas_force`, `gas_energy`, `orifice_force`, `strut_force`, `tyre_force`), the Lagrange model (`axle_position`, `mass_matrix_reduced`, `generalized_forces`, `rhs_reduced`, `dissipated_power`, `energy_terms`) and the constrained Newton–Euler model (`constraints`, `constraint_jacobian`, `constraint_gamma`, `stroke_full`, `applied_forces_full`, `accelerations_full`); `Gear`, `stop_force`, `initial_stroke`, `initial_state`, both simulators, `drop_metrics`, `full_state_from_reduced` and `main` are given. Tests (14) prove: the gas energy is the integral of the gas force; the orifice force is odd, quadratic and dissipative; the tyre cannot pull; the mass matrix is the kinetic energy of the kinematics; the generalised forces are minus the potential gradient (virtual work); the hanging stroke; the energy balance closes to 10⁻⁷ and nothing is dissipated without dampers; Newton–Euler reproduces Lagrange to 10⁻⁹ m; the multipliers equal the joint and guide forces computed from the Lagrange accelerations; an initial violation drifts linearly without stabilisation and decays with it; the load factor and stroke rise with sink rate; and a drop without lift loads the gear harder. |
| **Practical activity** | **computational experiment**: a drop-test parameter study of sink rate, orifice size, discharge coefficient, polytropic exponent, lift and rake, with the energy residual reported for every run; and a constraint-stabilisation study. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `landing_gear.png`, the derivation of both sets of equations, the energy budget at peak stroke, the parameter-study tables, and a half-page note on which strut parameter you would change to cut the load factor by 10 % and what it costs. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: both formulations derived and implemented; all public tests pass; the energy residual reported for every run; the note gives a load factor with its stroke and efficiency and names the assumption most likely to make it wrong. |
| **References** | `@ntrs-naca-tr-1154-landing-gear-behavior` Milwitzky and Cook, "Analysis of Landing-Gear Behavior" (NACA-TR-1154, 1953): Mechanics of Landing Gear (Dynamics of System; Forces in Shock Strut: Hydraulic force, Pneumatic force, Internal friction force; Forces on Tire); Equations of Motion; Parameter Studies (Effect of Orifice Discharge Coefficient; Effect of Air-Compression Process); Evaluation of Analysis by Comparison of Calculated Results with Experimental Data (supplementary-free): the two-mass model and its force laws, which this lesson follows. `@mit-ocw-16-07-dynamics` L20 "Energy methods: Lagrange's equations" (required-free): generalised coordinates and forces. `@mit-ocw-16-61-aerospace-dynamics` lecture 6 "Numerical Solution of Nonlinear Differential Equations" (supplementary-free). `@ocw-16-90` "1.7 Stiffness and Implicit Methods" (required-free): why the stop and tyre make the system stiff. `@scipy-userguide` "Ordinary differential equations (solve_ivp)" (required-free). Constraint-force formulations (Lagrange multipliers, the augmented system, Baumgarte stabilisation) are original instruction. |
| **Extensions** | ★ Add Coulomb friction in the strut bearings, proportional to the side load (TR-1154 "Internal friction force"), and show how rake now changes the load factor. ★ Let the airframe pitch (free θ₁ with the gear ahead of the centre of mass): the constraint on the sliding joint becomes nonlinear and drift appears without stabilisation. ★ Add wheel spin-up: a rotating wheel with friction at the contact patch (`P06.06`) produces a drag load and a fore-and-aft strut bending moment. ★ Add a recoil (rebound) valve: a smaller orifice for extension, and show that it stops the airframe leaving the ground after the bounce. |
| **Time** | unmeasured |

## 1. The engineering problem

A landing gear is an energy absorber. At touchdown the aircraft's vertical velocity, the sink rate, must be brought to zero over the strut's stroke and the tyre's deflection, and the force that does it is what the gear, the wing or fuselage attachment and everything in the aircraft feel. Certification sets a limit sink rate (3.05 m/s, 10 ft/s, is typical for transport aircraft at their design landing weight) and proves the gear by **drop tests**: the gear on a heavy carriage, dropped onto a platform from a height that gives that sink rate, with wing lift simulated by reducing the effective weight. The quantities recorded are the peak vertical ground force, often quoted as a **load factor** (force / weight), the peak stroke (the strut must not bottom), and the strut's efficiency.

The heart of most gears is the **oleo-pneumatic strut**: a cylinder with gas above and oil below, and a piston that forces oil through an orifice as it slides in. The gas is a nonlinear spring; the orifice is a damper whose force goes with the square of the stroke rate. Between them they make a nearly constant force over the stroke, which is what an energy absorber should do: for a given energy and stroke, a constant force is the smallest possible peak force.

This lesson builds a two-body model of a drop test and writes its equations two ways. It follows the classical two-mass analysis of Milwitzky and Cook (NACA-TR-1154), and adds what modern multibody codes do: the joint as a constraint with a force you can compute. Success criteria, written first: the two formulations agree to 10⁻⁹ m; the constraint forces equal those computed independently; the energy put in equals the energy stored plus dissipated to 10⁻⁷ of the total at every instant; and the strut does not bottom at the limit sink rate.

## 2. Understand and derive

### 2.1 The bodies, the joints and the force elements

The model has two bodies and a ground:

```
            ║ guide (x1 = 0, th1 = 0)
        ┌───╨───┐  body 1: airframe / drop carriage, mass M, at height y1; lift L = W acts on it
        │   ●   │  attachment point = body 1's reference point
        └───┬───┘
            ╲  strut axis e, raked beta from vertical      gas + orifice + stop: force F along e
             ╲ ← sliding joint: body 2 moves along e only
              ╲     stroke s (s > 0 compresses)
               ◎  body 2: piston, axle, wheel; mass m; axle at (x2, y2)
        ═══════════════  ground: tyre force F_t vertical, compression only; no friction
```

- **Body 1**, the airframe or drop carriage, mass M = 1500 kg. In a drop rig it is guided vertically and cannot pitch.
- **Body 2**, the "unsprung mass": the piston, axle, brake and wheel, m = 40 kg. A sliding (prismatic) joint lets it move only along the strut axis, which is raked β = 8° from the vertical.
- **The strut** acts along the axis between the two bodies with a net force F (positive pushing them apart), made of three parts (TR-1154, "Forces in Shock Strut").
  - *Pneumatic*: gas at p₀ = 1.5 MPa in a volume V₀ = 3.14 L at full extension, compressed polytropically (pVⁿ constant, n = 1.1 here) by the piston of area A_a = 7.854 × 10⁻³ m²:
    **F_gas(s) = A_a[p₀(V₀/(V₀ − A_a s))ⁿ − p_atm],**
    10.99 kN at full extension and 53.3 kN at the 0.30 m end of stroke.
  - *Hydraulic*: oil of density ρ pushed by an area A_h through an orifice of area A_o with discharge coefficient C_d. Bernoulli across the orifice, with the oil flow rate A_h ṡ through an effective area C_d A_o, gives a pressure drop ρ(A_h ṡ)²/(2(C_d A_o)²) acting on A_h:
    **F_hyd = ρ A_h³ ṡ|ṡ| / (2(C_d A_o)²),**
    written with ṡ|ṡ| so that it always opposes the motion.
  - *Top-out stop*: a stiff spring k_stop (10⁷ N/m) that resists extension beyond s = 0, so the gas precharge holds the strut fully extended in the air.
- **The tyre**: a spring k_t = 800 kN/m with a light damper c_t = 500 N·s/m, acting vertically at the axle when the axle is lower than the tyre radius R. It can push but not pull (TR-1154, "Forces on Tire", uses a measured force–deflection curve instead).
- **Lift** L = W = (M + m)g on body 1, the drop-test convention, so the aircraft's weight is balanced and only the kinetic energy of the sink has to be absorbed.

### 2.2 Formulation 1: Lagrange in reduced coordinates

The guide and the sliding joint leave two degrees of freedom. Choose the airframe height y₁ and the stroke s. The axle sits a length ℓ₀ − s below the attachment along the axis:

x₂ = −(ℓ₀ − s) sin β,  y₂ = y₁ − (ℓ₀ − s) cos β,

so ẋ₂ = ṡ sin β and ẏ₂ = ẏ₁ + ṡ cos β. The kinetic energy is

T = ½Mẏ₁² + ½m[(ṡ sin β)² + (ẏ₁ + ṡ cos β)²] = ½ q̇ᵀ M_r q̇,  **M_r = [[M + m, m cos β], [m cos β, m]].**

M_r is constant, so Lagrange's equations (`P06.11`; 16.07 L20) reduce to **M_r q̈ = Q**. The generalised forces come from virtual work. A virtual δy₁ moves both bodies up by δy₁; a virtual δs moves the axle by δs along the axis (δy₂ = cos β δs) and compresses the strut by δs, against its force F. With F_t the tyre force (vertical, at the axle):

**Q_y₁ = −(M + m)g + L + F_t,  Q_s = (F_t − mg) cos β − F.**

The internal strut force does not appear in Q_y₁: moving the whole gear up does not change the strut's length. The sliding joint's normal force and the guide's reaction do no virtual work and vanish from the equations. That is the strength of the reduced formulation, and its weakness: the forces the joints carry are not computed.

**Check by potential.** When nothing dissipates (c_t = 0, ṡ = 0), Q must be −∂V/∂q with V = Mgy₁ + mgy₂ − Ly₁ + U_gas(s) + ½k_t d² + ½k_stop(−s)²₊. The test does this numerically: it is the virtual-work derivation checked by an independent route.

### 2.3 The gas spring's stored energy

The work done compressing the gas from s = 0 to s is

U_gas(s) = ∫₀ˢ F_gas ds′ = **p₀V₀/(n − 1) [(V₀/(V₀ − A_a s))^{n−1} − 1] − p_atm A_a s**

(substitute V = V₀ − A_a s′ and integrate pVⁿ = p₀V₀ⁿ). The energy test checks it against quadrature. As n → 1 (isothermal) it tends to p₀V₀ ln(V₀/V) − p_atm A_a s. TR-1154's parameter study "Effect of Air-Compression Process" varies n between the isothermal and adiabatic values; real drop impacts sit in between because they are fast but the oil and gas mix.

### 2.4 Before touchdown

In the air the strut sits on its top-out stop. With lift equal to weight the gear descends at constant speed and the strut holds the unsprung weight; with no lift the gear falls freely and the strut carries nothing. In general the gear accelerates at −g(1 − L/W) with s̈ = 0, and the s-equation gives F = −(L/W) m g cos β. So the hanging stroke s₀ solves F_gas(s₀) − k_stop(−s₀) = −(L/W) mg cos β: **s₀ = −1.134 mm** for the lesson's gear. The initial state is then: tyre just touching (y₂ = R), s = s₀, ẏ₁ = −3.05 m/s, ṡ = 0.

### 2.5 Formulation 2: constrained Newton–Euler in full coordinates

Now give each body its three planar coordinates, q = (x₁, y₁, θ₁, x₂, y₂, θ₂), and write each joint as a **constraint equation** Φ(q) = 0:

Φ₁ = x₁ (guide),  Φ₂ = θ₁ (guide),  Φ₃ = θ₂ − θ₁ (sliding joint: no relative rotation),  Φ₄ = (r₂ − r₁)·n(θ₁) (sliding joint: no motion across the axis),

where e(θ₁) and n(θ₁) are the axis direction and its normal, fixed in body 1: e₀ = (sin β, cos β), n₀ = (cos β, −sin β), rotated by θ₁. Six coordinates and four constraints leave the same two degrees of freedom.

Newton–Euler for each body, with the joints replaced by unknown **constraint forces**:

**M q̈ = Q + Φ_qᵀ λ,**

with M = diag(M, M, I₁, m, m, I₂), Q the applied forces (gravity, lift, tyre, and the strut force: +F e on body 1, −F e on body 2, plus a θ₁ term −F (r₂ − r₁)·e′ from the virtual work), Φ_q the 4 × 6 Jacobian of the constraints and λ the four **Lagrange multipliers**. The term Φ_qᵀλ is the set of joint forces: it acts only in directions the constraints forbid, so it does no work on any allowed motion. Each multiplier is a physical force: λ₁ is the guide's horizontal reaction, λ₄ the force across the sliding joint, the **side load on the strut bearings**.

Six equations, ten unknowns (q̈ and λ). The other four come from differentiating the constraints twice: Φ̇ = Φ_q q̇ = 0 and

**Φ_q q̈ = γ,  γ = −(Φ_q q̇)_q q̇.**

For the linear constraints Φ₁–Φ₃, γ = 0. For Φ₄ = d·n with d = r₂ − r₁ and ṅ = θ̇₁ n′, ṅ′ = −θ̇₁ n: Φ̈₄ = d̈·n + θ̈₁ d·n′ + 2θ̇₁ ḋ·n′ − θ̇₁² d·n, so γ₄ = −2θ̇₁ ḋ·n′ + θ̇₁² d·n. Together, one linear solve per evaluation:

```
 [ M     −Φ_qᵀ ] [ q̈ ]   [ Q ]
 [ Φ_q    0    ] [ λ ] = [ γ ]
```

This is the **augmented** or saddle-point form used by general multibody codes: no reduced coordinates have to be found by hand, at the price of a larger system and a new problem, drift.

### 2.6 Drift and Baumgarte stabilisation

The integrator enforces Φ̈ = 0, not Φ = 0. Any error in Φ or Φ̇ (from round-off, the integrator, or an inconsistent start) is not corrected: Φ(t) = Φ(0) + Φ̇(0)t grows linearly. Baumgarte's method replaces Φ̈ = 0 with a damped equation,

**Φ̈ + 2αΦ̇ + β_B²Φ = 0,**

that is, γ → γ − 2αΦ_q q̇ − β_B²Φ in the augmented system, so violations decay like a critically damped oscillator when α = β_B. Choosing the gains is a compromise: too small and drift is slowly removed; too large and the equations become stiff.

There is a subtlety. In this drop rig, θ₁ stays exactly zero, so along the motion all four constraints are *linear* in q. Explicit Runge–Kutta methods preserve linear invariants exactly (each stage is a linear combination of slopes that satisfy Φ_q q̇ = 0), so from a consistent start the drift stays at round-off, 1.9 × 10⁻¹⁶ here. Drift becomes a practical problem with nonlinear constraints, such as a pitching airframe (Extension 2). The lesson therefore tests stabilisation by starting with a deliberate violation.

### 2.7 Energy accounting

Energy is the independent check of both formulations. Count, from touchdown to time t:

- **input** = the kinetic energy lost, ½(M + m)v₀² − T(t), plus the work done by gravity and lift;
- **stored** = the gas energy U_gas(s) − U_gas(s₀), the tyre's ½k_t d², the stop's ½k_stop(−s)²₊ (relative to touchdown);
- **dissipated** = ∫(F_hyd ṡ + (F_t − k_t d)ḋ) dt, integrated as an extra state of the ODE.

The balance input = stored + dissipated must hold at every instant. With lift = weight, gravity and lift cancel on body 1, and the input at the end of the compression stroke is essentially the sink energy: 7073 J of the 7163 J. The other 90 J is the work of the lift's excess over the airframe's own weight: lift = (M + m)g acts on body 1 only, so it does −mg cos β Δs on the stroke.

### 2.8 Load factor and efficiency

The **gear load factor** is the peak vertical ground force divided by the weight, n = F_t,max/W. The **strut efficiency** is the work the strut force does over the compression stroke divided by the peak strut force times the stroke:

η = ∫F ds / (F_max Δs).

A linear spring has η = 0.5 (a triangle under its force–stroke curve); a perfectly constant force has η = 1. An oleo strut's gas force rises with stroke while its orifice force falls as the stroke rate falls, and their sum can be nearly flat: efficiencies of 0.8–0.9 are typical.

### 2.9 Where the model fails

- **Friction in the strut bearings** (TR-1154 "Internal friction force"), proportional to the side load: it adds to the load factor, and it couples the rake to the vertical response.
- **Wheel spin-up.** At touchdown the stationary wheel slides on the runway until friction spins it up (`P06.06`); the drag load bends the strut backwards and springs it forwards. It is excluded here (frictionless ground).
- **Tyre.** A real tyre's force–deflection curve is nonlinear and it can bottom (TR-1154 "Impact With Tire Bottoming").
- **Orifice.** Real struts use metering pins (orifice area varying with stroke) to flatten the force, and recoil valves to damp the rebound. Without one, this model's airframe leaves the ground after the bounce.
- **Rigid airframe.** The airframe is a mass, not a flexible structure; wing bending modes respond to the landing impulse and can raise local loads.
- **A guided airframe.** A real aircraft pitches; the drop rig does not.

## 3. Build the computational model

Open `starter/main.py`. `Gear` (the parameters), `stop_force`, `initial_stroke`, `initial_state`, both simulators, `drop_metrics`, `_rot`, `_axis`, `full_state_from_reduced` and `main` are given.

### Step 1: force elements

`gas_force`, `gas_energy`, `orifice_force`, `strut_force` (gas + orifice − stop) and `tyre_force` (returning the total force, clipped at zero, and its elastic part) are 2.1 and 2.3.

### Step 2: the Lagrange model

`axle_position` and `mass_matrix_reduced` are 2.2. `generalized_forces` evaluates the tyre force at the axle height and rate, then Q. `rhs_reduced` solves M_r q̈ = Q and appends the dissipated power from `dissipated_power`. `energy_terms` builds the three budget terms of 2.7 from a run.

### Step 3: the constrained Newton–Euler model

`constraints`, `constraint_jacobian` and `constraint_gamma` are 2.5 (use `_axis`, which returns e, n and their θ₁-derivatives). `stroke_full` gives s = ℓ₀ + (r₂ − r₁)·e and its rate. `applied_forces_full` assembles Q. `accelerations_full` builds the 10 × 10 augmented matrix and solves it with `np.linalg.solve`, with the Baumgarte terms of 2.6.

```python
A = np.zeros((10, 10))
A[:6, :6], A[:6, 6:], A[6:, :6] = Mq, -J.T, J
sol = np.linalg.solve(A, np.concatenate([Q, gamma - 2 * alpha * (J @ qd) - beta_b ** 2 * Phi]))
qdd, lam = sol[:6], sol[6:]
```

Both simulators use `solve_ivp` with DOP853 at rtol 10⁻¹⁰ and a maximum step of 2 ms. The stop (10⁷ N/m on 40 kg, about 80 Hz) and the tyre (about 22 Hz on the unsprung mass) make the problem mildly stiff, and the contact switches make the right-hand side non-smooth; an explicit high-order method with tight tolerances handles both here at modest cost (`P03.06`; 16.90 "1.7 Stiffness and Implicit Methods").

Running `main()` prints

```
gear: M = 1500 kg, m = 40 kg, rake 8 deg, W = 15102 N, lift = 1.0 W, sink 3.05 m/s, KE = 7163 J
strut: gas force at full extension 10985 N, at 0.3 m 53335 N; hanging stroke s0 = -1.1337 mm
static stroke on the ground (no lift): 85.8 mm

Lagrange (reduced) drop: peak tyre force 32191 N, gear load factor 2.132, peak stroke 229.4 mm (bottomed: False), peak airframe deceleration 2.149 g, strut efficiency 0.907
  at peak stroke (t = 0.1658 s): input 7073.4 J = stored 4578.8 J + dissipated 2494.6 J; residual -1.34e-09 J
  max |residual| over the run 5.83e-05 J (8.2e-09 of the largest input)
  at t = 0.50 s: stored 1.0 J, dissipated 4045.1 J, stroke -1.5 mm

Newton-Euler (full, 6 coordinates, 4 constraints) against Lagrange: max |y1 diff| 4.23e-12 m, max |s diff| 4.37e-11 m, peak tyre force 32191 N
  constraint drift max |Phi| from a consistent start: 1.90e-16 (explicit Runge-Kutta preserves the constraints, linear along this motion, to round-off)
  perturbed start, Baumgarte alpha = beta = 0/s: |Phi_4| at 0, 0.1, 0.2 s = 0.990, 5.942, 10.893 mm
  perturbed start, Baumgarte alpha = beta = 50/s: |Phi_4| at 0, 0.1, 0.2 s = 0.990, 0.073, 0.001 mm
  strut bearing side load at peak tyre force: 4309 N (tyre force x sin(beta) = 4480 N); guide reaction -24 N
```

then the sink-rate table of part 4, and writes `landing_gear.png`: forces against time (tyre by both formulations, strut, bearing side load), the strut force against stroke with the gas force alone, and the energy budget.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-18-multibody-simulation`. The 14 tests: the gas energy against quadrature to 10⁻¹⁰ and its monotone stiffening; the orifice force odd, quadratic and dissipative; the stop and tyre one-sided; the mass matrix against the kinetic energy of `axle_position` by numerical differentiation; the generalised forces against −∇V at three states (virtual work, independently); the hanging stroke for lift = weight and for free fall; the energy balance to 10⁻⁷ of the energy absorbed, with monotone dissipation; zero dissipation and stored = input without dampers; Newton–Euler against Lagrange (y₁ and s to 10⁻⁹ m, tyre force and dissipated energy); the multipliers λ₄ and λ₁ against the forces the Lagrange accelerations require; linear drift without stabilisation and decay below 1 % in 0.2 s with it; the drop metrics; and a heavier load without lift.

### Compare

The two formulations agree to 4.2 × 10⁻¹² m in y₁ and 4.4 × 10⁻¹¹ m in s over 0.5 s: they share the force elements but not the equations, the coordinates or the way the joints are handled. The multipliers agree with the independent route of the test: the bearing side load at peak tyre force is **4309 N**, less than F_t sin β = 4480 N by the inertia and weight of the unsprung mass across the axis; the guide reaction, −24 N, is m ẍ₂, the only horizontal acceleration in the system. For comparison with an established analysis, the force laws and the two-mass structure are those of NACA-TR-1154, whose "Evaluation of Analysis" section compares calculated and measured drop tests and found good agreement for normal impacts; this lesson's numbers are not that report's gear, so the comparison is of method, not of values.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Every run reports its energy residual; a run whose budget does not close is not a result. Write each prediction before running.

1. **Sink rate.** 1.5–4.0 m/s, lift = weight. Predict how n scales (if the force were constant, n ∝ v²/stroke). *Finding:*

   | sink [m/s] | 1.50 | 2.50 | 3.05 | 4.00 |
   |---|---|---|---|---|
   | gear load factor | 1.304 | 1.765 | 2.132 | 2.986 |
   | peak stroke [mm] | 96.3 | 189.1 | 229.4 | 281.1 |
   | efficiency | 0.882 | 0.901 | 0.907 | 0.876 |
   | max energy residual [J] | 3.1 × 10⁻⁷ | 1.3 × 10⁻⁷ | 1.4 × 10⁻⁷ | 2.2 × 10⁻⁷ |

   The stroke grows less than v² (the gas spring stiffens) and the load factor more than linearly; at 4 m/s the strut uses 94 % of its stroke and its efficiency falls as the gas force climbs.

2. **Orifice size** (the designer's main knob). A_o from 3.0 to 6.0 cm² at 3.05 m/s. Predict whether the optimum is at the smallest or the largest orifice. *Finding:* n = 2.622, 2.345, 2.132, 2.205, 2.360 and stroke 204, 219, 229, 244, 253 mm for 3.0, 3.5, 4.0, 5.0, 6.0 cm² (19.5–27.6 mm diameter); efficiencies 0.848, 0.876, 0.907, 0.796, 0.705. Too small an orifice makes a large damping spike early in the stroke; too large leaves the gas spring to stop the aircraft, with a late peak. The 4.0 cm² orifice is near the minimum load factor, and the minimum coincides with the highest efficiency.

3. **Discharge coefficient** (TR-1154 "Effect of Orifice Discharge Coefficient"). C_d = 0.6, 0.7, 0.8, 0.9. *Finding:* n = 2.380, 2.132, 2.121, 2.231; strokes 217, 229, 239, 246 mm. C_d enters only as C_d A_o, so it acts exactly like the orifice area; an uncertain C_d of ±0.1 moves the load factor by up to 12 %, which is why it is measured.

4. **Polytropic exponent** (TR-1154 "Effect of Air-Compression Process"). n = 1.0 (isothermal), 1.1, 1.2, 1.4 (adiabatic). *Finding:* load factors 2.127, 2.132, 2.137, 2.273; strokes 235, 229, 224, 213 mm. A stiffer gas law shortens the stroke and raises the load only at the adiabatic end, where the gas force near the end of stroke takes over from the orifice.

5. **Lift and rake.** Lift = 0 against lift = weight; rake 0°, 8°, 15°. *Finding:* without lift the gear must absorb 12.6 kJ instead of 7.1 kJ; the load factor rises to 3.652 and the stroke reaches 300.8 mm, past the 300 mm available: **the strut bottoms**. Rake barely changes the vertical response (n = 2.100, 2.132, 2.213) but sets the bearing side load at the peak: 0, 4309, 8323 N. The side load is what bearing friction acts on, so in a real strut rake does change the load factor (Extension 1).

6. **Baumgarte gains.** Start with the axle 1 mm and 0.05 m/s off the strut axis; α = β_B = 0, 10, 50, 200 s⁻¹. *Finding:* |Φ₄| at 0.2 s = 10.89, 1.742, 0.00094 and 5 × 10⁻¹¹ mm; the peak stroke changes by under 0.5 mm in every case. Without stabilisation the violation grows exactly linearly (Φ₀ + Φ̇₀t, 10.89 mm at 0.2 s); gains of 50 s⁻¹ remove it within the first impact.

Raw data format: CSV tables `sink_mps, n_gear, s_max_mm, efficiency, max_residual_J` (1); `A_o_cm2, d_mm, n_gear, s_max_mm, efficiency` (2); `Cd, n_gear, s_max_mm, efficiency` (3); `n_poly, n_gear, s_max_mm, efficiency` (4); `lift_ratio, rake_deg, n_gear, s_max_mm, side_load_N` (5); `alpha, phi4_0p05_mm, phi4_0p2_mm` (6); plus `landing_gear.png`.

## 5. Validate and improve

Three independent checks support the result.

| check | quantity | result | tolerance |
|---|---|---|---|
| Lagrange against constrained Newton–Euler | max \|Δy₁\|, max \|Δs\| over 0.5 s | 4.2 × 10⁻¹² m, 4.4 × 10⁻¹¹ m | integrator, rtol 10⁻¹⁰ |
| energy: input = stored + dissipated | at peak stroke | 7073.4 = 4578.8 + 2494.6 J, residual 1.3 × 10⁻⁹ J | — |
| energy, whole run | max \|residual\| / largest input | 8.2 × 10⁻⁹ | test: 10⁻⁷ |
| multipliers against the Lagrange accelerations | λ₄ (side load), λ₁ (guide) at three instants | agree to 10⁻⁶ (test) | — |
| constraint drift, consistent start | max \|Φ\| | 1.9 × 10⁻¹⁶ | round-off |
| constraint drift, 1 mm violation, no stabilisation | Φ₄ at 0.1 s, 0.2 s | 5.942, 10.893 mm = Φ₀ + Φ̇₀t | exact |
| constraint drift, Baumgarte 50 s⁻¹ | Φ₄ at 0.2 s | 0.00094 mm | < 1 % of Φ₀ |

At the peak stroke the energy budget says where the sink energy went: 4021 J into the gas, 564 J into the tyre (37.6 mm deflection), 2495 J into heat in the oil and the tyre, and almost nothing left as kinetic energy. The gas then gives most of its energy back: by 0.5 s, 4045 J has been dissipated and the airframe has rebounded; with no recoil valve the tyre leaves the ground. That is the model's most visible flaw, and a real strut fixes it with a separate orifice for extension.

The engineering answer: **at 3.05 m/s with lift = weight the strut uses 229 of its 300 mm, the peak ground force is 32.2 kN (load factor 2.13) at 34 ms after touchdown, the strut efficiency is 0.91 and the bearings carry 4.3 kN of side load.** The orifice size and its discharge coefficient move the load factor most (±10–20 % across plausible values); the gas exponent matters only at the adiabatic end; rake matters for the bearings. A drop *without* lift would bottom the strut: the 10 ft/s sink-rate requirement is only meaningful together with its lift assumption.

What to improve: bearing friction (Extension 1), a measured tyre curve, a metering pin and a recoil valve; and, before trusting any of these numbers for a real aircraft, a comparison with drop-test data as TR-1154 did.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of M_r, Q, the gas energy, Φ_q and γ; `landing_gear.png`; the six experiment tables, each with its energy residual; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the load factor, the stroke, the efficiency and the side load, the parameter you would change to cut the load factor by 10 %, and what that change costs in stroke.

### Adversarial self-checks

- **What if the sink rate doubles?** The energy quadruples. From experiment 1, 4 m/s already uses 94 % of the stroke; 6.1 m/s takes it to 339 mm, past the end, with a load factor of 6.3, and bottoming is not modelled (the stroke runs on with only the gas force). The `bottomed` flag is the model telling you it is outside its range.
- **What if the orifice goes to zero or to infinity?** Zero: an infinitely stiff damper, the strut locks and the tyre alone absorbs the energy, with a very high load factor. Infinite: only the gas spring acts, the efficiency falls towards that of a stiffening spring, and the airframe bounces. Experiment 2 shows both trends.
- **Which assumption fails first?** The frictionless sliding joint: real bearings carry the 4.3 kN side load with a friction coefficient of perhaps 0.05–0.1, adding a few hundred newtons of axial force. Then the tyre law, then the absent spin-up drag.
- **What measurement could falsify the model?** A drop test with the ground force, stroke and stroke rate recorded (TR-1154's "Evaluation of Analysis"). The model predicts the time of the peak force (34 ms) and its shape; a measured peak much later than predicted would mean the orifice force is overestimated (C_d too low or flow not quasi-steady).
- **Does the result survive a change of formulation or integrator?** Lagrange and Newton–Euler agree to 10⁻¹¹ m; rtol from 10⁻⁸ to 10⁻¹² changes the load factor by less than 10⁻⁹; the energy residual stays below 10⁻⁷.
- **Which decision changes if the uncertainty doubles?** If C_d is known to ±0.2 instead of ±0.1, the load factor spans 2.1–2.7 (2.715 at C_d = 0.5) instead of 2.1–2.4 (experiment 3), and a gear sized to a 2.5 limit has no margin at the low C_d end. The decision changes from "trust the analysis" to "measure C_d in a flow test before the drop test".

## Check yourself

1. Name the three force elements of an oleo-pneumatic strut and say how each depends on stroke or stroke rate.
2. Why does the strut's internal force not appear in the generalised force for the airframe height y₁?
3. What does a Lagrange multiplier represent physically in the constrained Newton–Euler equations? Which multiplier is the strut's bearing side load?
4. Why does constraint drift stay at round-off in this drop rig, and when would it not?
5. Write Baumgarte's stabilised constraint equation and say what the gains do.
6. What is strut efficiency, and what values would a linear spring and a constant force give?
7. With lift equal to weight, what energy must the gear absorb at a sink rate v? What changes without lift?
8. Which strut parameter most changes the load factor in this lesson's study, and why does the discharge coefficient act exactly like it?
