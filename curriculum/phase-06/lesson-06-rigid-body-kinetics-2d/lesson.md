# P06.06  Rigid Body Kinetics in 2D

> A sliding tyre turns exactly as much kinetic energy into heat as it gives the wheel, whatever the friction coefficient. Kinetics tells you how long it takes; momentum tells you how much it costs.

| | |
|---|---|
| **Learning outcomes** | 1. Write the plane equations of motion of a rigid body, ΣF = m a_G and ΣM_G = I_G α, with a rolling constraint or a Coulomb slip law, and decide from them whether a body rolls or slides. 2. Derive the touchdown spin-up of a landing-gear wheel in closed form by impulse–momentum: the friction impulse, the spin-up time, the final speed, the slip distance and the friction energy, for a constant and for a ramped normal load. 3. Prove that the friction energy equals ½ M_eff v_s0², independent of the friction coefficient and the load history, and check it two ways: as ∫F v_s dt and as the loss of kinetic energy. 4. Implement and test a hybrid sliding/rolling/locked simulation with `solve_ivp` events that reproduces the closed forms and closes its energy audit to 10⁻⁹. 5. Decide what brake torque a wheel can take without skidding, and explain why a locked tyre stops an aircraft more slowly than a rolling one. |
| **Prerequisites** | `P06.05` (plane rigid-body kinematics, the rolling constraint v = ωR), `P06.03` (work–energy and impulse–momentum), `P04.04` (rotation, moment of inertia, τ = Iα), `P05.05` (moments of inertia), `P04.02` (work and energy), `P04.03` (momentum), `P03.06` (adaptive integrators and event detection). Skill check: you can draw the free-body diagram of a wheel on the ground and say in which direction friction acts on it when the axle moves faster than the rim. |
| **Engineering question** | A 60 t airliner touches down at 70 m/s. One main wheel (120 kg rotating, rolling radius 0.56 m) is not turning, and its normal load rises to 120 kN in 0.1 s. How long does the tyre slide, how far, how large is the drag load on the gear, how much energy goes into the tyre as heat, and how much brake torque can the wheel take afterwards without skidding? |
| **Assumptions and limits** | Planar motion on a flat runway; the aircraft's share carried by the wheel is a single translating mass M; the wheel is rigid with a fixed rolling radius; Coulomb friction with constant kinetic and static coefficients; the normal load is a given function of time (the strut dynamics are `P06.18`); aerodynamic drag, lift changes and rolling resistance during the 0.1 s spin-up are neglected. Real tyres deform, and their friction depends on slip ratio, speed, temperature and the runway surface: the model gives the right energy and impulse but only an estimate of the time history. |
| **Analytical work** | Free-body diagrams of the wheel and the axle; ΣF = m a_G, ΣM_G = I_G α; the incline: rolling acceleration, the friction needed to roll, the sliding case; the effective mass 1/M_eff = 1/M + R²/I; impulse–momentum for the spin-up; the spin-up time and slip distance for constant and ramped loads; the friction energy two ways and its limits (heavy aircraft, pre-spin); the friction needed to roll under braking and the largest brake torque without skid. |
| **Computational deliverable** | `starter/main.py`: `incline_motion`, `min_friction_to_roll`, `effective_mass`, `spin_up_closed_form`, `kinetic_energy_loss`, `rolling_friction_for_brake`, `max_brake_torque`, `simulate_touchdown`, `energy_audit` (`Wheel`, `ramp_load`, `Touchdown` and `main` are given). Tests prove: the incline equations and the hoop–disk–sphere order; the sliding threshold; the impulse balance of the spin-up; the friction energy by impulse–momentum equal to the kinetic-energy loss for 20 random wheels; its independence of μ_k and of the load ramp; the heavy-aircraft and pre-spin limits; continuity of the two ramp branches; the simulation equal to the closed forms to 10⁻⁸ for three load histories; the energy audit closed to 10⁻⁹; free rolling with zero friction; braking below the limit rolling at the predicted deceleration; braking above it skidding and locking; an overspun wheel pushing the axle forward. |
| **Practical activity** | **computational experiment**: sweeps of the friction coefficient, the load-ramp time and the wheel pre-spin; a brake-torque study around the skid limit; a solver-tolerance study of the energy audit. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `wheel_spinup.png`, the derivations of the effective mass, the spin-up closed form and the friction energy by both routes, the five experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the friction energy is derived by both routes (not quoted); all public tests pass; the simulated spin-up time and friction energy are compared with the closed forms with numbers; the note states the brake-torque limit and the assumption that most limits the answer. |
| **References** | `@mit-ocw-16-07-dynamics` L21 "2D rigid body dynamics" (sections Conservation of Angular Momentum about the Center of Mass; Radius of Gyration; Summary: Governing Equations), L22 "2D rigid body dynamics: work and energy", L23 "2D rigid body dynamics: impulse and momentum" (required-free: ΣF = m a_G and ΣM_G = I_G α, and the work and impulse forms this lesson uses as its two routes). `@openstax-university-physics-vol1` 11.1 Rolling Motion (Rolling Motion without Slipping; Rolling Motion with Slipping; Rolling Down an Inclined Plane with Slipping) (required-free: the incline warm-up and the end of slip at v = ωR). `@libretexts-baker-haynes-engineering-statics` 9.1 Dry Friction (9.1.1 Coulomb Friction); 10.8 Mass Moment of Inertia (required-free). `@libretexts-cline-variational-principles` 13.24 The Rolling Wheel (supplementary-free). `@mit-ocw-2-003sc-engineering-dynamics` Recitation 4 Notes: Torque and Angular Momentum, Pendulum with Torsional Spring, Rolling Pipe on Moving Truck (supplementary-free: a worked rolling problem with an accelerating support). |
| **Extensions** | ★ Replace the Coulomb law by a slip-ratio friction curve μ(κ), κ = v_s/V, rising to a peak near κ ≈ 0.1 and falling towards the locked value, and show that the friction energy is unchanged while the time history is not. ★ Model an anti-skid controller that releases the brake when the slip ratio exceeds a threshold and compare its stopping distance with the rolling and locked limits. ★ Add a gear leg as a spring in the fore–aft direction and compute the spring-back oscillation the spin-up drag load excites. |
| **Time** | unmeasured |

## 1. The engineering problem

When an airliner lands, its main wheels hit the runway at 60–80 m/s without turning. For a tenth of a second the tyres slide. The friction force that drags them backwards is the same force that spins them up: it acts at the contact patch, below the axle, so it both decelerates the axle and turns the wheel. You can see the result on any runway as black rubber marks at the touchdown zone, and as the puff of smoke at every landing. The same force reaches the airframe as a **spin-up drag load**: a sharp backward pull on the landing gear, followed by a forward spring-back when the wheel reaches rolling speed and the load suddenly drops. Gear legs and their attachments are sized for it.

This lesson asks four questions about one illustrative main-gear wheel of a 60 t single-aisle airliner. The wheel, tyre and brake rotors weigh 120 kg; the rolling radius is 0.56 m; the radius of gyration is 0.65 R, so I = 15.90 kg·m². The wheel carries one quarter of the aircraft, an equivalent translating mass M = 15 000 kg. It touches down at 70 m/s with no rotation; the normal load rises linearly to 120 kN in 0.1 s as the strut compresses; the tyre–runway kinetic friction coefficient is 0.55 and the static coefficient 0.8 (both *illustrative*; dry concrete under a rolling or skidding tyre spans roughly 0.5–0.9). How long does the tyre slide, and how far? What is the peak drag load? How much energy goes into the tyre as heat? And once the wheel rolls, how hard can the pilot brake before it skids?

The answers come from plane rigid-body kinetics: Newton's second law for the centre of mass, its rotational partner about the centre of mass, and a law for the force at the contact. They are checked two ways. A closed form from impulse–momentum gives the final state and the friction energy without integrating anything. A hybrid time integration with event detection gives the time history. Success criteria, written first: the two routes to the friction energy agree to the solver tolerance; the simulated spin-up time equals the closed form to 10⁻⁸; the energy audit (kinetic energy + friction heat + brake energy) closes to 10⁻⁹ at every sample.

## 2. Understand and derive

### 2.1 Plane kinetics of a rigid body

`P06.05` described the motion of a rigid body in a plane: the velocity of its centre G and one angular velocity ω. Kinetics adds the forces. For a rigid body in plane motion, two vector laws hold. The centre of mass moves as if all the mass were there and all the forces acted there, and the angular momentum about G changes at the rate of the moments about G:

**ΣF = m a_G,   ΣM_G = I_G α.**   (1)

The second is `P04.04`'s τ = Iα applied about the centre of mass (16.07 L21, "Summary: Governing Equations"). It holds about G even when G accelerates; about any other point it picks up an extra term, which is why the safe habit is to take moments about G. Equation (1) gives three scalar equations (two forces, one moment) for a body in a plane.

Two kinds of contact close the system. In **rolling** (no slip) the contact point has zero velocity, so v_G = ωR and a_G = αR (`P06.05`): one kinematic equation, and the friction force f is an unknown to be found. It is acceptable only if |f| ≤ μ_s N. In **sliding** the contact point moves at the slip velocity v_s = v_G − ωR ≠ 0, and the friction force is the Coulomb law (Baker & Haynes 9.1.1; OpenStax 11.1, "Rolling Motion with Slipping"):

**f = μ_k N, directed against the slip velocity.**   (2)

The rule for solving any problem of this kind: assume rolling, solve for the friction it needs, and check it against μ_s N. If the check fails, the body slides; replace the constraint by (2) and solve again.

### 2.2 Warm-up: rolling or sliding down an incline

A round body of mass m, radius R and moment of inertia I = k² m R² about its centre is released on an incline of angle θ. With x down the slope, f the friction up the slope and N the normal force, (1) gives

m a = m g sin θ − f,   N = m g cos θ,   I α = f R.

If it rolls, a = αR, and eliminating f:

**a = g sin θ / (1 + I/(mR²)),   f = m g sin θ · I/(I + mR²).**   (3)

A hoop (k² = 1) accelerates at g sin θ/2, a disk (½) at 2g sin θ/3, a solid sphere (2/5) at 5g sin θ/7: the more of the mass sits far from the axis, the more of gravity's work goes into spin and the slower the body goes. On a 10° slope the code gives 0.8515, 1.1353 and 1.2164 m/s² for the three, against g sin θ = 1.7029 m/s² without friction. Rolling needs f ≤ μ_s N, which with (3) is

**μ_s ≥ tan θ / (1 + mR²/I).**   (4)

For the disk, μ_s ≥ tan θ/3. Below that the body slides (OpenStax 11.1 works this case under "Rolling Down an Inclined Plane with Slipping"): a = g(sin θ − μ_k cos θ) and α = μ_k m g cos θ R/I, and the contact point slips down the slope because the centre outruns the rim.

### 2.3 The touchdown: free-body diagrams

```
                 axle (share of aircraft, mass M)
                    │  V →
               ┌────●────┐          wheel: I about the axle, radius R
               │    │ ω  │          ω positive when rolling forward (clockwise in this view)
               │    ↻    │
               └────┬────┘
        ←── F  ═════╧═════  runway
                    ↑ N
                contact: slips forward at v_s = V − ωR
```

Take the axle and the aircraft share it carries as one body of mass M moving at V, and the rotating parts as a wheel of moment of inertia I about the axle. At the contact the runway pushes up with N and backwards with F. While the tyre slides forward (v_s > 0), (2) gives F = μ_k N. The axle's equation and the wheel's moment equation about its centre are

**M dV/dt = −F,   I dω/dt = F R − T_b,**   (5)

with T_b a brake torque (zero for now). The friction both slows the axle and spins the wheel. The slip velocity changes at

d v_s/dt = dV/dt − R dω/dt = −F (1/M + R²/I) = −F / M_eff,   **1/M_eff = 1/M + R²/I.**   (6)

M_eff, the **effective mass**, is the mass the slip "feels": a friction force F destroys slip as if it acted on a single mass M_eff. Because M is huge here, M_eff = 50.529 kg is almost exactly I/R² = 50.70 kg: the aircraft barely notices, and it is the wheel's rotational inertia, referred to the rim, that has to be brought up to speed.

### 2.4 The spin-up by impulse–momentum

The tyre stops sliding when v_s = 0. Integrate (5) over the sliding phase with J = ∫F dt, the friction impulse (`P06.03`; 16.07 L23):

M (V₀ − V_f) = J,   I (ω_f − ω₀) = J R,   V_f = ω_f R.

Three equations, three unknowns. With v_s0 = V₀ − ω₀R the initial slip,

**J = M_eff v_s0,   V_f = V₀ − J/M,   ω_f = ω₀ + J R/I.**   (7)

Equation (7) needs neither μ_k nor N(t): the final state of the spin-up is fixed by momentum alone. For the lesson's wheel, J = 3537.0 N·s, V_f = 69.7642 m/s (the aircraft loses 0.24 m/s), and ω_f = 124.579 rad/s (1190 rpm).

**Spin-up time.** The time is set by how fast friction delivers that impulse, J = μ_k ∫₀^{t_s} N dt. For a constant load, t_s = J/(μ_k N). With N ramping linearly to N_max in t_r, ∫N dt = N_max t²/(2t_r) during the ramp and N_max(t − t_r/2) after it, so

**t_s = √(2 J t_r / (μ_k N_max))  if J ≤ ½ μ_k N_max t_r,  otherwise  t_s = J/(μ_k N_max) + t_r/2.**   (8)

For the lesson's wheel: constant 120 kN load, t_s = 3537.0/66 000 = **53.59 ms**; with the 0.1 s ramp, J > ½ μ_k N_max t_r = 3300 N·s, so the second branch applies: t_s = 53.59 + 50 = **103.59 ms**. The two branches meet when J = ½ μ_k N_max t_r, where both give t_s = t_r (a test).

**Slip distance.** The tread slides over the runway a distance s = ∫v_s dt. From (6), v_s falls at F/M_eff. For a constant load it falls linearly, so s = ½ v_s0 t_s = **1.876 m**. During a ramp it falls as v_s0 − c t² with c = μ_k N_max/(2 M_eff t_r); integrating each piece gives **4.831 m** for the lesson's ramp. That is the length of a black mark on the runway.

**Peak drag load.** The drag on the axle is F = μ_k N, so it peaks at μ_k N at the moment spin-up ends, here 0.55 × 120 kN = **66.0 kN**, more than half the vertical load, applied for 0.1 s and then removed in an instant. That sudden release is what makes a gear leg spring forward.

### 2.5 The friction energy, two ways

**Route 1: heat at the contact.** The friction force does work on the tyre–runway pair at the rate F v_s (force times sliding speed). The energy dissipated is E_f = ∫F v_s dt. With F dt = −M_eff dv_s from (6),

**E_f = ∫ F v_s dt = −M_eff ∫_{v_s0}^{0} v_s dv_s = ½ M_eff v_s0².**   (9)

Again neither μ_k nor N(t) appears. A higher friction coefficient makes a larger force for a shorter time over a shorter slip distance; the energy is the same.

**Route 2: kinetic energy.** Energy is conserved overall, so the heat must equal the kinetic energy lost:

E_f = ½ M V₀² + ½ I ω₀² − (½ M V_f² + ½ I ω_f²),   with V_f, ω_f from (7).   (10)

Substituting (7) and simplifying gives (9) again (Problem 2). The code evaluates (10) directly, without using (9), and the test checks the two against each other for 20 random wheels. For the lesson's wheel both give **E_f = 123.797 kJ**.

**Two limits that check it.**

- **Heavy aircraft (M → ∞).** M_eff → I/R², V_f → V₀, and (9) becomes ½ I (V₀/R)² = ½ I ω_f²: **the heat equals the kinetic energy the wheel gains**. The aircraft loses twice the wheel's final kinetic energy: half goes into spin, half into heat. For the lesson's wheel the aircraft loses 247.18 kJ, the wheel gains 123.38 kJ, and 123.80 kJ becomes heat. It is the same split as charging a capacitor through a resistor, and for the same reason: a fixed "source" (the aircraft's speed) brings a store up to its own level through a dissipative link.
- **Pre-spin (ω₀ = f V₀/R).** v_s0 = (1 − f)V₀, so E_f falls as (1 − f)²: half the touchdown speed as pre-spin removes three quarters of the heat (123.797 → 30.949 kJ), 90 % removes 99 %. Pre-rotating wheels to save tyres has been proposed many times; the energy is easy to save, the difficulty is matching the speed reliably.

### 2.6 Rolling and braking

Once rolling, V = ωR. With no brake, (5) with F = 0 is consistent: the wheel rolls without needing friction, and V stays constant (rolling resistance and drag are neglected). With a brake torque T_b, eliminate dV/dt = R dω/dt from (5):

**F = T_b M R / (I + M R²),   max T_b without skid = μ_s N (I + M R²)/(M R).**   (11)

For M ≫ I/R² the limit is μ_s N R, the obvious result, increased slightly because the brake must also decelerate the wheel's own spin. For the lesson's wheel at 120 kN and μ_s = 0.8: **53.94 kN·m**. A rolling, braked tyre at 90 % of that limit pulls the aircraft share back at F/M = 5.760 m/s². Above the limit the tyre skids: the friction falls from the static value it was supplying to μ_k N, which no longer balances the brake, and the wheel decelerates at (T_b − μ_k N R)/I until it **locks**. A locked tyre drags at only μ_k N, 4.40 m/s² here, and grinds one spot of the tread. That is why anti-skid systems release the brake as soon as the wheel slows faster than the aircraft: the best braking is at the edge of static friction, not beyond it.

### 2.7 Where the model fails

- **Tyre deformation.** A real tyre is not rigid: the contact patch deforms, part of the tread grips and part slips, and the friction depends on the **slip ratio** κ = v_s/V. Peak friction typically comes at κ ≈ 0.1–0.2, not at full sliding. The impulse (7) and the energy (9) still hold, because they need only momentum and the final rolling state; the time history changes.
- **Constant μ.** At 70 m/s, the tread surface heats by hundreds of kelvin in milliseconds, and the rubber's friction changes with temperature. A wet runway can cut μ by a factor of two or more, doubling t_s and the slip distance but not the energy.
- **Given normal load.** N(t) is really set by the strut and the sink rate (`P06.18`). A tyre that bounces leaves the ground with a partly spun-up wheel; the model treats this only through ω₀.
- **One translating mass.** The gear leg is flexible fore and aft, so the axle does not move with the aircraft. The spin-up drag load excites the leg (the ★ extension).

## 3. Build the computational model

Open `starter/main.py`. The dataclass `Wheel` (the lesson's illustrative wheel), `ramp_load`, the result dataclass `Touchdown` and `main` are given.

### Step 1: the incline

`incline_motion(theta, m, I, R, mu_s, mu_k)` follows 2.1: compute N and the friction needed to roll from (3); if it is within μ_s N, return the rolling accelerations; otherwise the sliding ones. `min_friction_to_roll` is (4).

### Step 2: the closed forms

`effective_mass` is (6). `spin_up_closed_form` computes J from (7), then the time and slip distance from the right branch of (8), and returns the final state, the energy (9) and the peak drag μ_k N(t_s). It raises `ValueError` if v_s0 ≤ 0 (the closed form assumes forward slip). `kinetic_energy_loss` is (10), written independently of (9). `rolling_friction_for_brake` and `max_brake_torque` are (11).

### Step 3: the hybrid simulation

`simulate_touchdown` integrates the states [x, V, ω, E_f, E_b, s] with DOP853. Two of them are bookkeeping: E_f grows at |F v_s| and E_b (the brake's energy) at T_b ω. The motion has three phases, each with its own right-hand side, and each ends with a terminal `solve_ivp` event:

```python
if phase == "sliding":                                   # F = mu_k N sign(v_s)
    rhs:    [V, -F/M, (F*R - T_b)/I, |F v_s|, T_b*omega, |v_s|]
    events: v_s = 0 (-> rolling),  omega = 0 under braking (-> locked)
elif phase == "rolling":                                 # F from (11)
    rhs:    [V, -F/M, -F/(M R), 0, T_b*omega, 0]
    events: mu_s N - |F| = 0 (-> sliding: a skid)
else:                                                    # locked, omega = 0, F = mu_k N
    rhs:    [V, -F/M, 0, F*V, 0, V]
    events: T_b - mu_k N R = 0 (-> sliding)
```

At each event the loop restarts `solve_ivp` from the event state in the new phase. Two details matter. On entering rolling, set V = ωR exactly, so that the constraint holds from the start of the phase and not only to the event tolerance. On leaving rolling, nudge ω by 10⁻¹² so the slip has the skid's sign and the new phase's slip event does not fire at once. Pass `max_step` small enough (here t_end/200) that the integrator cannot step over the switch-on of the brake. `energy_audit` returns the residual KE(0) − KE(t) − E_f − E_b, which must be zero at every sample.

Running `main()` prints the incline cases, the closed forms for a constant and a ramped load, the simulation and its differences from the closed form, the sweeps of Experiments 1 and 3, and two braking runs, and writes `wheel_spinup.png`: axle speed and rim speed meeting at spin-up; the drag load following μ_k N and dropping to zero; and the friction energy by both routes on one axis.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-06-rigid-body-kinetics-2d`. The tests:

1. **Incline**: (3) for a hoop, a disk and a sphere, with the returned friction satisfying both equations of (1); the order hoop < disk < sphere; sliding just below (4), with the sliding accelerations; no spin without friction.
2. **Effective mass**: the limits I/R² and M.
3. **Spin-up closed form**: rolling at the end; the linear and angular impulses both equal J; t_s = J/(μ_k N); s = ½ v_s0 t_s.
4. **Two routes**: (9) equal to (10) for 20 random wheels, loads, ramps and pre-spins.
5. **Independence and limits**: the same energy for nine combinations of μ_k, N and t_r; t_s halves when N doubles; the heavy-aircraft limit; the (1 − f)² law; the overspun case refused.
6. **Ramp branches**: equal at J = ½ μ_k N_max t_r; tending to the constant load as t_r → 0.
7. **Simulation**: t_s, V_f, ω_f, E_f and s equal to the closed forms to 10⁻⁸ for a constant load and for ramps that end before and after spin-up.
8. **Energy audit**: residual below 10⁻⁹; friction heat never decreasing; final heat equal to the kinetic-energy loss.
9. **Rolling and braking**: free rolling with zero slip and zero friction; braking at 0.9 × the limit rolling at the deceleration of (11); braking at 1.3 × the limit skidding, locking at the predicted time, and then dragging at μ_k N/M; an overspun wheel spinning down and pushing the axle forward.
10. `main` writes the figure.

### Compare

The independent check is the closed form against the integration, and the two routes to the energy against each other. For the lesson's ramped touchdown, the simulation's spin-up time is 103.5916 ms against 103.59 ms from (8), a difference of 2.2 × 10⁻¹¹ s; its friction energy ∫F v_s dt is 123.7966 kJ against (9), a difference of 2.3 × 10⁻⁵ J; the kinetic-energy loss is the same 123.7966 kJ; and the largest energy-audit residual over the run is 6.2 × 10⁻¹³ of the initial kinetic energy. These checks catch the classic mistakes at once: a friction torque with the wrong sign makes the wheel spin backwards and v_s never reaches zero; taking moments about the contact point without the extra term gives the wrong α and a spin-up time off by the ratio M/M_eff; and an event that is not terminal lets the integrator run on with Coulomb friction across v_s = 0, where the friction reverses every step and chatters.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are outputs of the model with illustrative parameters, not measurements of a real tyre. Write each prediction down before running it.

1. **Friction coefficient.** μ_k = 0.3, 0.55 and 0.8 with the 0.1 s ramp. Predict how t_s, the slip distance and E_f change. *Finding:* t_s 148.25, 103.59, 85.84 ms; slip 6.642, 4.831, 4.006 m; E_f 123.797 kJ in all three. On a wet runway the tyre slides further and longer, but heats no more.
2. **Load ramp.** t_r = 0, 0.05, 0.1 and 0.2 s at μ_k = 0.55. Predict which branch of (8) each falls in and what happens to the peak drag. *Finding:* t_s 53.59, 78.59, 103.59, 146.41 ms; slip 1.876, 3.490, 4.831, 6.833 m; peak drag 66.0 kN for the first three and 48.32 kN for t_r = 0.2 s, where spin-up ends while the load is still rising (the first branch of (8)). A softer strut lowers the spin-up drag load.
3. **Pre-spin.** ω₀ = 0, 0.5 and 0.9 of V₀/R. Predict E_f from (9). *Finding:* 123.797, 30.949 and 1.238 kJ, exactly (1 − f)²; t_s 103.59, 73.21 and 32.74 ms (less than proportional, because the load is still ramping).
4. **Braking at the limit.** After spin-up, apply a constant brake torque at 0.3 s of 0.9 and 1.2 times (11)'s 53.94 kN·m. Predict which run keeps rolling and the two decelerations. *Finding:* at 0.9 × the wheel rolls and the aircraft share decelerates at 5.760 m/s²; at 1.2 × the tyre skids at 0.3 s and locks 71.3 ms later (the prediction ω I/(T_b − μ_k N R) gives 71.3 ms too), after which it drags at 4.40 m/s². Held from 70 m/s, these decelerations would stop the share in 425 m and 557 m (a rough comparison only: lift, reverse thrust and the changing load are left out).
5. **Solver tolerance.** Repeat the ramped touchdown with rtol = atol = 10⁻⁴, 10⁻⁶, 10⁻⁸ and 10⁻¹⁰. Predict which quantity is most sensitive. *Finding:* the spin-up time error is −5.2 × 10⁻⁸, 2.5 × 10⁻⁸, 1.2 × 10⁻⁹ and −2.2 × 10⁻¹¹ s; the friction-energy error 4.2 × 10⁻², 2.2 × 10⁻², 5.5 × 10⁻² and 2.3 × 10⁻⁵ J; the audit residual 1.1 × 10⁻⁹, 3.4 × 10⁻⁹, 1.5 × 10⁻⁹ and 6.2 × 10⁻¹³. Even a loose tolerance gets the physics right to parts per million here, because the right-hand side is smooth inside each phase and the events locate the switches; `max_step` does more for accuracy at loose tolerances than rtol does.

Raw data format: CSV tables `mu_k, t_spin_ms, slip_m, E_f_kJ` (1); `t_ramp_s, t_spin_ms, slip_m, peak_drag_kN` (2); `prespin_fraction, E_f_kJ, t_spin_ms` (3); `brake_fraction, phases, decel_mps2, E_brake_kJ` (4); `rtol, t_spin_err_s, E_err_J, audit_residual` (5); plus `wheel_spinup.png`.

## 5. Validate and improve

| quantity | closed form (§2) | simulation | difference | cause |
|---|---|---|---|---|
| spin-up time, 0.1 s ramp (ms) | 103.59 (8) | 103.5916 | 2.2 × 10⁻¹¹ s | event location |
| final speed (m/s) | 69.7642 (7) | 69.7642 | < 10⁻⁹ relative | — |
| final wheel rate (rad/s) | 124.579 (7) | 124.579 | < 10⁻⁸ relative | — |
| slip distance (m) | 4.831 | 4.8315 | < 10⁻⁸ relative | — |
| friction energy, ∫F v_s dt (kJ) | 123.797 (9) | 123.7966 | 2.3 × 10⁻⁵ J | solver tolerance |
| kinetic-energy loss (kJ) | 123.797 (10) | 123.7966 | 2.3 × 10⁻⁵ J | solver tolerance |
| energy audit residual | 0 | 6.2 × 10⁻¹³ relative | — | rounding |
| lock time after a 1.2 × skid (ms) | 71.33 | 71.33 | < 10⁻⁹ s | — |

Every quantity the closed forms predict, the simulation reproduces to the solver tolerance, and the energy is the same by both routes. The engineering answer: **the tyre slides for 104 ms and 4.8 m, the drag load peaks at 66 kN (0.55 of the vertical load), and 124 kJ goes into the tyre and runway as heat whatever the friction coefficient; the wheel then takes up to 54 kN·m of brake torque at 120 kN before it skids.**

Two things to improve, in order. First the friction law: a slip-ratio curve (★ extension) would change the time history and the peak drag load, which matter for the gear leg, but not the energy. Second the load history: coupling the wheel to a strut model (`P06.18`) would replace the assumed ramp. What not to improve: the solver tolerance, which Experiment 5 shows is already far below the uncertainty in μ.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4), (6), (7), (8), (9) and (11), and the reduction of (10) to (9); `wheel_spinup.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the touchdown speed doubles?** J doubles and E_f quadruples (v_s0²): 495.2 kJ instead of 123.8 kJ. With the 0.1 s ramp, t_s grows only from 103.6 to 157.2 ms (the ramp's 50 ms does not scale), but the slip distance nearly triples, from 4.83 to 13.96 m. Tyre heat rises with the square of the landing speed, one reason fast landings wear tyres.
- **What if μ_k goes to zero?** The tyre never spins up: t_s → ∞, the slip distance → ∞, and E_f stays ½ M_eff v_s0², delivered ever more slowly. The limit is singular in time and regular in energy.
- **What if the wheel had no inertia?** I → 0 gives M_eff → 0: no impulse, no heat, instant spin-up. A light wheel is cheap to spin up; the brake rotors, the heaviest rotating part, dominate I.
- **Which assumption fails first in the real world?** The rigid tyre with constant μ: real spin-up is shaped by the tyre's slip–friction curve. The energy result survives; the drag-load history does not.
- **What measurement could falsify the model?** A wheel-speed record at touchdown (anti-skid systems measure it): ω(t) should rise as ∫N dt; a spin-up much faster than (8) would mean μ is higher than assumed, a slower one that the load builds more slowly.
- **Does the result survive a change of solver or tolerance?** Yes: Experiment 5. The closed forms do not depend on any solver.
- **Which decision changes if the uncertainty doubles?** If μ_s were known only to ±0.2 instead of ±0.1, the safe brake torque would have to be set from the lower value, 40.5 kN·m at μ_s = 0.6 instead of 47.2 kN·m at μ_s = 0.7: the anti-skid system exists precisely because μ_s is not known in advance.

## Check yourself

1. Write the two plane equations of motion of a rigid body. About which point is the moment equation simplest, and why?
2. A disk and a hoop of the same mass and radius roll down the same incline. Which arrives first, and what is the ratio of their accelerations?
3. Define the effective mass M_eff of the touchdown problem. What is it when the aircraft is very heavy?
4. Why does the final rolling speed after spin-up not depend on the friction coefficient?
5. Give the friction energy dissipated during spin-up in terms of M_eff and the initial slip velocity. What fraction of the aircraft's lost kinetic energy is it when M → ∞?
6. A wheel is pre-spun to 70 % of the touchdown speed. By what factor does the friction energy fall?
7. Why does a locked tyre decelerate an aircraft less than a braked rolling one?
8. In the hybrid simulation, why must the sliding phase end with a terminal event at v_s = 0 instead of letting Coulomb friction run on?
