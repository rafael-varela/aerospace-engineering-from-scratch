# P04.04  Rotation & Angular Momentum

> A spinning top does not fall. It turns sideways, at a rate you can predict to four decimals from one vector equation.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the inertia tensor of simple bodies, move it to another point with the parallel-axis theorem in tensor form, and compute rotational kinetic energy and angular momentum from it. 2. Derive Euler's equations I ω̇ + ω × (Iω) = τ from dL/dt = τ in a rotating frame, carry the attitude with a unit quaternion, and integrate both with `scipy.integrate.solve_ivp`. 3. Show that torque-free motion conserves the inertial angular momentum and the kinetic energy, derive the body-frame nutation rate of a symmetric body and the exponential growth rate of spin about the intermediate axis, and confirm both in simulation. 4. Explain how a reaction wheel turns a spacecraft without external torque, and size a rest-to-rest slew in closed form. 5. Derive the gyroscopic precession rate of a top, both the fast-top estimate mgl/(I₃ω₃) and the exact steady-precession root, the minimum spin for steady precession and the nutation turning point, and verify them against a full simulation that conserves the top's three invariants. |
| **Prerequisites** | `P04.01` (Newton's laws), `P04.02` (work and energy), `P04.03` (momentum, systems of particles), `P01.07` (dot and cross products), `P01.15` (rotation matrices), `P01.16` (eigenvalues), `P02.06` (linearisation about an equilibrium), `P03.06` (adaptive integrators). Skill check: you can compute a × b for two 3-vectors by hand, and you can say why the eigenvalues of a symmetric 3 × 3 matrix are real. |
| **Engineering question** | A 100 g disk spinning at 50 rev/s on a 40 mm stem is tilted 30° from the vertical. How fast does it precess, how slowly may it spin before steady precession becomes impossible, and how far does it nod if released without being started on its cone? And the same physics on a spacecraft: about which axes can a 500 kg box spin stably, and what wheel torque and wheel momentum turn it 30° in a minute? |
| **Assumptions and limits** | Rigid bodies; principal axes known; a frictionless pivot for the top (real tops slow down and their pivot slides); uniform gravity; for the spacecraft, no external torque (no gravity gradient, drag, solar pressure or magnetic torque) and a perfectly balanced wheel with an ideal motor; small perturbations in the stability analysis. The equations of motion themselves are exact for a rigid body. |
| **Analytical work** | Inertia tensor of a cylinder and a box; the parallel-axis theorem I = I_cm + m(|d|²E − ddᵀ); rotational kinetic energy and angular momentum; Euler's equations from the transport theorem; quaternion kinematics; conservation of L and T in torque-free motion; the symmetric body's nutation rate (I₃ − I₁)ω₃/I₁; the intermediate-axis growth rate; the reaction-wheel gyrostat equation; bang-bang slew; the gyroscope Ω = τ/(L sin θ); the steady-precession quadratic and its minimum spin; the nutation turning point from the invariants. |
| **Computational deliverable** | `starter/main.py`: `inertia_solid_cylinder`, `inertia_box`, `parallel_axis`, `rotational_energy`, `angular_momentum_inertial`, `euler_rates`, `quaternion_rate`, `simulate_rigid_body`, `symmetric_nutation_rate`, `intermediate_axis_growth_rate`, `principal_axis_wobble_frequency`, `bang_bang_slew`, `gravity_torque_body`, `top_initial_state`, `steady_precession_rate`, `fast_top_precession`, `top_invariants`, `nutation_turning_tilt`, `axis_tilt_and_azimuth`, `measured_precession_rate` (`cross3`, the quaternion helpers, `Rotation`, `Top` and `main` are given). Tests prove: the inertia formulas and their rod, disk and cube limits; the parallel-axis theorem against the perpendicular distance; Euler's equations for a sphere, from rest and with a wheel; the quaternion rate keeps \|q\| = 1; torque-free motion conserves L to 10⁻⁸ and T to 10⁻⁹; stable spin about the extreme axes wobbles at the predicted frequency; the intermediate axis grows at the predicted rate to 2 %; the symmetric nutation rate to 10⁻⁶; the simulated slew against the closed form with body plus wheel momentum zero; the gravity torque; the precession root, its fast-top and 90° limits and its refusal of too slow a spin; the simulated top precessing at the root with three invariants held; the nutation turning point to 10⁻⁴. |
| **Practical activity** | **computational experiment**: a spinning-top simulation (steady precession against spin, the release without precession, the minimum spin), the intermediate-axis instability of a box satellite, and a reaction-wheel slew with a momentum audit. |
| **Equipment** | None. A computer with the environment from `P00.01`. A toy top or a bicycle wheel is a welcome way to see the effects, but nothing is measured with it. |
| **Safety** | None beyond normal computer use. If you try a spinning bicycle wheel by hand, hold it by the axle, keep fingers away from the spokes and keep it away from your face. |
| **Evidence** | Completed `main.py`, passing tests, `spinning_top.png`, the derivations of Euler's equations, the intermediate-axis growth rate, the steady-precession quadratic and the turning point, the experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: Euler's equations are derived (not quoted); all public tests pass; the simulated precession is compared with both the fast-top and the exact formula with numbers; the conservation checks are reported. |
| **References** | `@openstax-university-physics-vol1` Ch. 10 "Fixed-Axis Rotation": 10.1 Rotational Variables, 10.3 Relating Angular and Translational Quantities, 10.4 Moment of Inertia and Rotational Kinetic Energy, 10.5 Calculating Moments of Inertia, 10.6 Torque, 10.7 Newton's Second Law for Rotation; Ch. 11 "Angular Momentum": 11.2 Angular Momentum, 11.3 Conservation of Angular Momentum, 11.4 Precession of a Gyroscope (required-free: the fixed-axis foundations and the fast-gyroscope precession formula, which this lesson extends to three dimensions). |
| **Extensions** | ★ Add a passive nutation damper (a small mass on a spring inside the body) and show that a body spinning about its minor axis ends up spinning about its major axis. ★ Replace the frictionless pivot by a pivot with a small friction torque against the spin and find when the top falls. ★ Three orthogonal reaction wheels: slew about an arbitrary axis with a quaternion-error feedback law and compare its wheel momentum with the bang-bang slew. |
| **Time** | unmeasured |

## 1. The engineering problem

Every spacecraft must point: its antenna at Earth, its camera at a target, its solar arrays at the Sun. In space there is nothing to push against, so a spacecraft turns by exchanging angular momentum between its body and spinning wheels inside it, or keeps its orientation by spinning like a top. Both rely on one law, dL/dt = τ, applied to a body whose inertia changes direction as it turns. That law in three dimensions holds three surprises that a fixed-axis treatment never shows: the spin axis of a heavy top moves sideways instead of falling; a body spun about the wrong axis tumbles by itself; and a free body's rates change all the time while its angular momentum in space stays fixed.

This lesson writes the general equations, Euler's equations for the body rates and quaternion kinematics for the attitude, and checks them against closed forms in three settings. A **spinning top**: a 100 g disk of 30 mm radius on a stem 40 mm from the pivot, spinning at 50 rev/s and tilted 30°. The gyroscope formula predicts a precession rate of 2.7747 rad/s; the exact steady-precession equation gives 2.8667 rad/s, 3.3 % faster; the full simulation precesses at 2.86665 rad/s. Below 17.6 rev/s no steady precession exists at this tilt. Released without being pushed onto its cone, the top nods 2.2642° down, exactly as energy and angular-momentum conservation predict. A **box satellite** of 500 kg (1.0 × 1.5 × 2.0 m): spun at 0.5 rad/s about its largest or smallest axis it holds steady; about the middle axis a disturbance of 10⁻⁶ rad/s grows by e every 6.1 s and the satellite turns over after 82 s. And a **reaction wheel** that turns the same satellite 30° in 60 s with a 78.8 mN·m torque, ending at rest with the body and the wheel's angular momenta cancelling to 10⁻¹⁵.

Success criteria, written first: the inertial angular momentum and the energy conserved to the solver tolerance in every torque-free run; the simulated precession equal to the steady root to five figures; the intermediate-axis growth equal to the predicted rate; the slew equal to its closed form.

## 2. Understand and derive

### 2.1 From fixed-axis rotation to the inertia tensor

For rotation about a fixed axis, OpenStax Ch. 10 builds the angular variables θ, ω = dθ/dt, α = dω/dt (10.1), relates them to the speed v = rω of each particle (10.3), and defines the **moment of inertia** I = Σmᵢrᵢ² with rᵢ the distance to the axis, so that the kinetic energy is ½Iω² (10.4) and Newton's second law becomes τ = Iα (10.7). In three dimensions the angular velocity is a vector **ω**, each particle moves at **v** = **ω** × **r**, and the angular momentum is **L** = Σmᵢ **r**ᵢ × (**ω** × **r**ᵢ) (11.2). Expanding the double cross product gives a linear map,

**L = I ω,  I = Σmᵢ (|rᵢ|² E − rᵢ rᵢᵀ),**   (1)

the **inertia tensor**, a symmetric 3 × 3 matrix with E the identity. Its diagonal entries are the moments of inertia about the coordinate axes; the off-diagonal entries, the products of inertia, make L point away from ω. Because I is symmetric it has three real eigenvalues and orthogonal eigenvectors: the **principal moments** I₁, I₂, I₃ and **principal axes**. In principal axes I is diagonal and L = (I₁ω₁, I₂ω₂, I₃ω₃). The kinetic energy is

**T = ½ ω · I ω.**   (2)

For a uniform solid cylinder of mass m, radius r and length ℓ (10.5 integrates the simpler cases): I₃ = ½mr² about the axis and I₁ = I₂ = m(3r² + ℓ²)/12 about a diameter through the centre; ℓ → 0 is a disk (I₁ + I₂ = I₃, the perpendicular-axis theorem), r → 0 a thin rod (mℓ²/12). For a box with sides a, b, c: I = (m/12) diag(b² + c², a² + c², a² + b²). The 500 kg box satellite of 1.0 × 1.5 × 2.0 m has I = diag(260.4, 208.3, 135.4) kg·m².

**The parallel-axis theorem** in tensor form follows from (1) with rᵢ = d + ρᵢ, Σmᵢρᵢ = 0:

**I_P = I_cm + m (|d|² E − d dᵀ).**   (3)

Along any axis this adds m times the squared perpendicular distance; off the diagonal it creates products of inertia unless d lies along a principal axis. The lesson's top is a disk of 100 g, radius 30 mm and thickness 4 mm whose centre is 40 mm up its axis from the pivot: about its centre I₁ = 2.263 × 10⁻⁵ kg·m², and about the pivot I₁ = 2.263 × 10⁻⁵ + 0.1 × 0.04² = **1.8263 × 10⁻⁴ kg·m²**, eight times more; I₃ = ½mr² = **4.5 × 10⁻⁵ kg·m²** is unchanged.

### 2.2 Euler's equations

Newton's law for rotation about a fixed point or the centre of mass is dL/dt = τ with both vectors in an inertial frame (11.2; the torque τ = r × F is 10.6). The trouble is that I is constant only in **body axes**, which rotate with ω. For any vector **A**, the rate seen in the inertial frame equals the rate of its body components plus the rotation of the axes:

(dA/dt)_inertial = (dA/dt)_body + ω × A.   (4)

Applied to L = Iω with I constant in body axes:

**I ω̇ + ω × (I ω) = τ**   (Euler's equations).   (5)

In principal axes, component by component: I₁ω̇₁ = (I₂ − I₃)ω₂ω₃ + τ₁ and its cyclic permutations. The term ω × (Iω) is the **gyroscopic term**: it couples the three axes even with no torque, and it vanishes only when L ∥ ω, i.e. for spin about a principal axis or for a sphere.

**Attitude.** The body rates say how fast the body turns, not where it points. The orientation is a rotation matrix R (body to inertial; `P01.15`), carried here by a unit **quaternion** q = [cos(φ/2), sin(φ/2) n] for a rotation by φ about the unit axis n, which has no singularities and four numbers instead of nine. Its kinematics are

**q̇ = ½ q ⊗ [0, ω],**   (6)

with ⊗ the quaternion product and ω in body axes. Equation (6) keeps |q| = 1 exactly (q · q̇ = 0, a test); numerically it drifts at the level of the solver tolerance, and the code renormalises on output.

### 2.3 Torque-free motion

With τ = 0, the inertial angular momentum L_i = R I ω is constant (11.3), and dotting (5) with ω gives d(½ω·Iω)/dt = 0: **the kinetic energy is constant too**. In body axes the two conservation laws say that ω lies on an ellipsoid of constant energy and on an ellipsoid of constant |L|², and their intersection is the curve the tip of ω follows. The body rates change while L in space is fixed: the body tumbles around a fixed vector.

**Symmetric body (I₁ = I₂).** Euler's third equation gives ω̇₃ = 0, and the first two become ω̇₁ = −λω₂, ω̇₂ = λω₁ with

**λ = (I₃ − I₁) ω₃ / I₁:**   (7)

the transverse rate rotates in the body at the constant rate λ, in the sense of the spin for an oblate body (I₃ > I₁) and against it for a prolate one. For I₁ = 800, I₃ = 1200 kg·m² at ω₃ = 1 rad/s, λ = 0.5000 rad/s, and the simulation's transverse rate turns at 0.5000 rad/s.

**Spin about a principal axis: stable or not.** Let the body spin at ω₂ about axis 2 with small rates ε₁, ε₃ on the other axes. Linearising Euler's equations (`P02.06`):

I₁ ε̇₁ = (I₂ − I₃) ω₂ ε₃,  I₃ ε̇₃ = (I₁ − I₂) ω₂ ε₁  ⇒  ε̈₁ = (I₂ − I₃)(I₁ − I₂) ω₂² ε₁ / (I₁ I₃).   (8)

If axis 2 has the largest or smallest moment the two brackets have opposite signs and ε oscillates at **√((Iₐ − I_j)(Iₐ − I_k)/(I_j I_k)) |ω|** about axis a: 0.2402 rad/s about axis 1 and 0.2049 rad/s about axis 3 for the box at 0.5 rad/s. If axis 2 is the **intermediate** axis the brackets have the same sign and ε grows as e^{σt} with

**σ = |ω₂| √((I₁ − I₂)(I₂ − I₃)/(I₁ I₃)).**   (9)

For the box at 0.5 rad/s, σ = **0.1641 s⁻¹**. Starting from 10⁻⁶ rad/s, the disturbance reaches the size of the spin after about ln(0.5/10⁻⁶)/σ = 80 s; the simulation's spin first reverses sign at 82.4 s, then flips back and forth periodically: the tennis-racket or Dzhanibekov effect. Every flip conserves L and T exactly; nothing external happens.

An honest caveat that the rigid model hides: a real spacecraft is not rigid and dissipates energy (fuel slosh, flexing antennas). At fixed |L|, T = L²/(2I) is lowest for spin about the **major** axis, so dissipation drives any spinner there. A spinner about its minor axis is stable in the rigid model but not in practice (Problem 3).

### 2.4 Reaction wheels

A wheel inside the body stores angular momentum h (body axes). The total, R(Iω + h), is conserved if no external torque acts. With the motor torque on the wheels ḣ, Euler's equation for the body becomes

**I ω̇ + ω × (I ω + h) = τ − ḣ:**   (10)

the body feels the motor's reaction. For a rest-to-rest turn by an angle Θ in a time T about a principal axis, the simplest profile is **bang-bang**: accelerate with a constant torque for T/2, decelerate for T/2. Each half covers ½(τ/I)(T/2)², so

**τ = 4IΘ/T²,  ω_max = 2Θ/T,  h_max = I ω_max.**   (11)

For the box turned 30° about axis 3 in 60 s: τ = **78.78 mN·m**, ω_max = **1.000 °/s**, h_max = **2.3635 N·m·s**. A wheel of I_w = 0.02 kg·m² must spin at 118.2 rad/s = 1128 rpm relative to the body to hold that. Wheels have two limits: torque (motor size) and momentum (top speed); a constant external torque, however small, fills a wheel up over time and must be dumped with thrusters or magnetic torquers.

### 2.5 The gyroscope and the heavy top

**Fast-top estimate (OpenStax 11.4).** A top spinning with L₃ = I₃ω₃ about its axis, tilted θ from the vertical, feels the gravity torque τ = mgl sin θ, horizontal and perpendicular to the axis. If L₃ dominates every other angular momentum, dL/dt = τ turns L sideways without changing its size: the axis moves on a cone at the rate Ω with Ω L₃ sin θ = mgl sin θ, so

**Ω_fast = m g l / (I₃ ω₃),**   (12)

independent of θ. For the lesson's top, mgl = 0.03923 N·m and L₃ = 4.5 × 10⁻⁵ × 314.16 = 1.4137 × 10⁻² kg·m²/s, so Ω_fast = **2.7747 rad/s** (one cone every 2.26 s).

**Exact steady precession.** The fast-top estimate ignores the angular momentum of the precession itself. Put the top on a cone at constant θ, precessing at Ω about the vertical, with body spin component ω₃. In body axes (axis 1 horizontal, 2 in the vertical plane of the axis) the vertical is (0, sin θ, cos θ) and the body rate is ω = (0, Ω sin θ, ω₃). For this motion to be steady, the change of L is the precession itself, dL/dt = Ω ẑ × L, and it must equal the gravity torque. Working out ẑ × Iω and keeping the horizontal component:

Ω sin θ (I₃ω₃ − I₁ Ω cos θ) = m g l sin θ  ⇒  **I₁ cos θ Ω² − I₃ ω₃ Ω + m g l = 0.**   (13)

The smaller root is the slow precession that tops show (the larger, about I₃ω₃/(I₁ cos θ), is a fast precession that needs a hard push to start); it tends to (12) as ω₃ → ∞, and at θ = 90° (13) is linear with root (12) exactly. For the lesson's top at 30° the slow root is **2.8667 rad/s**, 3.3 % above the estimate. The discriminant must be positive:

**(I₃ ω₃)² ≥ 4 I₁ m g l cos θ,**   (14)

so below **ω₃ = 110.7 rad/s = 17.6 rev/s** at 30° no steady precession exists, and as friction slows a real top it starts to wobble and falls.

**Invariants and nutation.** The heavy top has three conserved quantities: the energy T + mgl cos θ (gravity does work only through the height of the centre of mass); the vertical component of L (gravity's torque is horizontal); and the spin component L₃ = I₃ω₃ (the torque has no component along the symmetry axis). Release the top at θ₀ with no precession. As it starts to fall, conservation of L_z forces it to precess, and the energy it gains by dropping goes into that precession. The axis nods between θ₀ and a lower turning point where θ̇ = 0; at that point the two conservation laws give (Problem 2)

**cos θ₀ − cos θ = 2 I₁ m g l sin²θ / L₃².**   (15)

For a fast top, θ − θ₀ ≈ 2I₁mgl sin θ₀ / L₃² = 2.054°; the exact root of (15) is **2.2642°**. The nodding is called **nutation**; for a fast top its frequency is near L₃/I₁ = 77 rad/s, much faster than the precession, which is why a fast top looks as if it precesses smoothly. In the simulation the nod repeats at 72.5 rad/s at 50 rev/s and 152.4 rad/s at 100 rev/s (against L₃/I₁ = 154.8): the fast-top estimate improves as the spin grows.

### 2.6 Where the model fails

- **Rigid body.** Real spacecraft flex and slosh; energy dissipation turns minor-axis spin into a tumble (2.3), which no rigid model shows.
- **Frictionless pivot.** A real top's point slides and rubs; the spin decays, (14) is eventually violated and the top falls. A rolling contact can even make a top rise (the "sleeping top"); that needs a friction model.
- **No external torque on the spacecraft.** Gravity gradient, magnetic, aerodynamic and solar-pressure torques of 10⁻⁶–10⁻³ N·m act on real satellites and accumulate in the wheels.
- **Ideal wheel.** Real wheels have friction, imbalance (vibration) and a maximum speed.

## 3. Build the computational model

Open `starter/main.py`. The helpers `cross3`, `quat_multiply`, `quat_to_matrix`, `quat_from_axis_angle`, the dataclasses `Rotation` and `Top`, and `main` are given.

### Step 1: inertia, energy and angular momentum

`inertia_solid_cylinder`, `inertia_box` and `parallel_axis` are 2.1 and (3), each returning a 3 × 3 array. `rotational_energy(I, w)` must accept one rate or a stack: `0.5 * np.sum(w * (w @ I.T), axis=-1)`. `angular_momentum_inertial(q, w, I, h_wheel)` returns R(q)(Iω + h) for one state or a stack.

### Step 2: Euler's equations and the attitude

```python
def euler_rates(w, I, torque, h_wheel=None, h_wheel_dot=None):
    h = np.zeros(3) if h_wheel is None else np.asarray(h_wheel, dtype=float)
    hd = np.zeros(3) if h_wheel_dot is None else np.asarray(h_wheel_dot, dtype=float)
    return np.linalg.solve(I, torque - cross3(w, I @ w + h) - hd)
```

`quaternion_rate` is (6). `simulate_rigid_body` integrates the ten states [q, ω, h] with DOP853 (rtol 10⁻¹¹ by default, atol 10⁻¹²), taking an external torque and a wheel torque as functions of (t, q, ω); it returns a `Rotation` with the quaternion renormalised. (`cross3` is a written-out cross product; `np.cross` gives the same numbers but is several times slower on single vectors, which matters when the right-hand side is called tens of thousands of times.)

### Step 3: closed forms

`symmetric_nutation_rate` is (7); `intermediate_axis_growth_rate` is (9); `principal_axis_wobble_frequency` is the oscillation frequency of 2.3 and raises `ValueError` for the intermediate axis; `bang_bang_slew` returns the three numbers of (11).

### Step 4: the top

`gravity_torque_body(q, m, l)` is (l e₃) × (m g_body) with g_body = Rᵀ(0, 0, −g). `top_initial_state(theta, w3, Omega)` returns the quaternion of a tilt θ about the inertial x axis and ω = (0, Ω sin θ, ω₃). `steady_precession_rate` is the slow root of (13), written as 2c/(−b + √(b² − 4ac)) to avoid cancellation when the spin is fast, with a `ValueError` when (14) fails; `fast_top_precession` is (12); `top_invariants` returns the energy, L_z and L₃ along a run; `nutation_turning_tilt` solves (15) with `brentq`; `axis_tilt_and_azimuth` and `measured_precession_rate` read the symmetry axis R e₃ and fit a line to its unwrapped azimuth.

Running `main()` prints the top's numbers, the steady and released runs with their invariants, a sweep of spin rates, the box satellite about each axis, the symmetric spinner, and the slew with its momentum audit, and writes `spinning_top.png` (the path of the axis tip of the released top; precession against spin rate with the two formulas and the simulated points; the body rates of the box spun about its intermediate axis).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-04-rotation-angular-momentum`. The tests:

1. **Inertia**: cylinder and box formulas; rod, disk (perpendicular axes) and cube limits.
2. **Parallel axis**: symmetric; moment about a line equal to I_cm + m r⊥²; products of inertia appear; zero offset changes nothing; positive definite.
3. **Energy and momentum** of known states, single and stacked, rotated, with a wheel.
4. **Euler's equations**: a sphere has no gyroscopic term; from rest ω̇ = τ/I; the three component formulas; the wheel torque acts on the body with the opposite sign.
5. **Quaternion rate**: q · q̇ = 0; the rate of a pure spin.
6. **Torque-free conservation**: L_i constant to 10⁻⁸, T to 10⁻⁹, |q| = 1, while the body rates change.
7. **Stable axes**: spin about axes 1 and 3 stays within 10 times the disturbance and wobbles at the predicted frequency (one FFT bin); the intermediate axis is refused.
8. **Intermediate axis**: the disturbance grows at σ of (9) to 2 %.
9. **Symmetric nutation**: λ = 0.5 rad/s to 10⁻⁶; ω₃ and the transverse magnitude constant; a prolate body turns the other way.
10. **Slew**: closed form; final angle to 10⁻⁶; peak rate; body plus wheel momentum zero to 10⁻¹².
11. **Gravity torque**: zero upright, mgl sin θ when tilted, no component along the axis or the vertical.
12. **Precession root**: solves (13); above the fast-top value below 90°; equal to it for a fast spin and at 90°; too slow a spin is refused.
13. **Steady top**: tilt constant to 10⁻⁶ rad, precession equal to the root to 10⁻⁵, three invariants held to 10⁻⁷.
14. **Released top**: starts at the upper turning point, reaches the turning point of (15) to 10⁻⁴; the fast-top estimate is within 15 % and low.
15. `main` writes the figure.

### Compare

The independent checks are the conservation laws and the closed forms (7), (9), (11), (13) and (15). The steady top precesses at 2.86665 rad/s, equal to the root of (13) to 10⁻¹⁰ relative; its tilt changes by 2.5 × 10⁻⁹ degrees in 1.5 s; its energy drifts by 4.4 × 10⁻¹³ and L_z by 2.8 × 10⁻¹¹ relative. The released top's nod, 2.2642°, matches (15) to every printed digit. These checks catch the common mistakes at once: a sign error in ω × Iω reverses the sense of the symmetric body's nutation (test 9), a gravity torque expressed in inertial instead of body axes breaks energy and L_z conservation at the 10⁻³ level within 1.5 s, and the fast-top formula alone would leave an unexplained 3 % gap.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results here are simulations of the model, not measurements of a real top. Write each prediction down before running it.

1. **Precession against spin.** For 20, 35, 50 and 100 rev/s at 30°, compute the fast-top and steady-root precession rates by hand, then simulate a steady start for 1 s. Predict which formula the simulation follows and where they part.
2. **Release without precession.** At 50 rev/s and 30°, release the top with Ω = 0. Predict the tilt swing from (15) and the fast-top estimate, then simulate 1.5 s and plot the path of the axis tip (it traces small loops, cusps, because the axis starts at rest).
3. **Minimum spin.** Compute (14) at tilts of 10°, 30°, 60° and 80°. Predict how it varies with the tilt.
4. **Box satellite.** Spin the 500 kg box at 0.5 rad/s about each principal axis with a 10⁻⁴ rad/s disturbance for 120 s. Predict which axis is unstable and how fast the disturbance grows; then fit the growth rate of a 10⁻⁶ disturbance over its small-amplitude stage only.
5. **Slew audit.** Turn the box 30° about axis 3 in 60 s with a bang-bang wheel torque. Predict the torque, peak rate and peak wheel momentum, then check the final angle and that the total angular momentum stays zero.

Raw data format: CSV tables `spin_revps, fast_radps, steady_radps, simulated_radps` (1); `t_s, tilt_deg, azimuth_deg` (2); `tilt_deg, min_spin_revps` (3); `axis, max_offaxis_radps, L_drift, T_drift` and `t_s, w1, w2, w3` (4); `t_s, angle_deg, body_rate, wheel_momentum, total_momentum` (5); plus `spinning_top.png`.

What you should find (simulation results): (1) fast-top 6.9368, 3.9639, 2.7747, 1.3874 rad/s; steady root 9.4173, 4.2530, 2.8667, 1.3983 rad/s; simulated equal to the root to four decimals at every spin, so the fast-top estimate is 26 % low at 20 rev/s and 0.8 % low at 100 rev/s; (2) swing 2.2642° against (15)'s 2.2642° and the estimate's 2.054°; mean precession 2.8652 rad/s; (3) 18.8, 17.6, 13.4 and 7.9 rev/s, falling to zero at 90°; (4) axes 1 and 3: off-axis rates stay below 1.8 × 10⁻⁴ rad/s; axis 2: they reach 0.40 rad/s, the spin reverses at 82.4 s, and the fitted growth rate is 0.1641 s⁻¹ against (9)'s 0.1641 s⁻¹; L and T held to 4 × 10⁻¹¹ in every case; (5) 78.78 mN·m, 1.000 °/s, 2.3635 N·m·s; final angle 30.0000°, final rate 1.2 × 10⁻¹² rad/s, total angular momentum below 2 × 10⁻¹⁵ N·m·s.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| I₁ of the top about the pivot (kg·m²) | 1.8263 × 10⁻⁴ | 1.8263 × 10⁻⁴ | 0 | — |
| fast-top precession (rad/s) | 2.7747 | 2.7747 | 0 | formula |
| steady precession, 50 rev/s (rad/s) | 2.86665 (root of 13) | 2.86665 (simulated) | 10⁻¹⁰ relative | solver tolerance |
| nutation swing, released (°) | 2.2642 (15); 2.054 fast-top | 2.2642 | 0; 10 % | fast-top estimate drops terms |
| nutation frequency (rad/s) | ≈ 77.4 (L₃/I₁) | 72.5 | 6 % | fast-top estimate |
| minimum spin at 30° (rev/s) | 17.6 | 17.6 | 0 | formula |
| intermediate-axis growth (s⁻¹) | 0.1641 | 0.1641 (fitted) | 10⁻⁵ relative | fit window |
| symmetric nutation (rad/s) | 0.5000 | 0.5000 | 10⁻¹² | — |
| slew final angle (°) | 30 | 30.0000 | 7 × 10⁻¹¹ relative | solver |
| total angular momentum, slew (N·m·s) | 0 | 1.8 × 10⁻¹⁵ | rounding | — |

Every quantity with an exact closed form matches to the solver tolerance; the only real differences are those of the fast-top approximations, which the lesson has quantified. The model's weaknesses are physical: the rigid-body assumption (no dissipation, so the minor-axis spinner looks stable), the frictionless pivot, and the absence of environmental torques on the satellite. The first improvement would be a damper (★ extension), the second a pivot friction torque.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (3), (5), (7), (9), (11), (13), (14) and (15); `spinning_top.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the precession rate by both formulas and the simulation, the minimum spin, the nod on release, the stable axes of the box and the slew's torque and wheel momentum.

### Adversarial self-checks

- **What if the spin halves?** The fast-top precession doubles (12), but the exact root rises faster (at 25 rev/s, 6.49 rad/s against 5.55 in the estimate) and the nod grows four times (it scales as 1/L₃²). Halve it again, to 12.5 rev/s, and (14) fails: no steady precession at 30°.
- **What if the top were a long rod instead of a disk?** I₃ becomes tiny; L₃ shrinks for the same spin rate and the fast-top assumption fails. A gyroscope is a disk for a reason: a large I₃/I₁ makes the spin dominate.
- **Which assumption fails first in the real world?** For the top, the frictionless pivot: the spin decays within seconds to tens of seconds. For a spinning spacecraft, rigidity: any dissipation drives it to its major axis.
- **What measurement could falsify the model?** Timing a real top's precession at a known spin (a strobe or a phone's slow-motion video) and comparing with (13); a result well outside the fast-top-to-exact range at a known tilt would point at friction or at a wrong inertia.
- **Does the result survive a change of solver?** The conservation laws are the test: at rtol 10⁻⁹, E and L_z stay within 10⁻⁹ relative in a 1.5 s run; a solver that drifted would show it there first. The closed forms (13) and (15) do not depend on any solver.
- **Which decision changes if the uncertainty doubles?** If the satellite's inertias were uncertain by 10 %, the box's ordering (260, 208, 135) would survive and the stability conclusions would not change; for a nearly symmetric body (two inertias within a few per cent) the intermediate axis is ambiguous, and spin about either near-equal axis must be treated as marginal.

## Check yourself

1. Write the inertia tensor of a uniform solid cylinder about its centre. Which entry changes when you move the reference point along the axis?
2. Why does L not point along ω for a general body? When does it?
3. Derive Euler's equations from dL/dt = τ in a rotating frame.
4. A free body's rates change continuously. Which two quantities stay constant, and in which frame?
5. Which principal axis is unstable for spin, and why does a real spacecraft prefer its major axis?
6. A reaction wheel spins up. What happens to the spacecraft, and why can the wheel not cancel a constant external torque forever?
7. State the fast-top precession rate. Why is it independent of the tilt, and what does the exact formula add?
8. Below what spin can a top not precess steadily, and what happens to a real top as friction brings it there?
