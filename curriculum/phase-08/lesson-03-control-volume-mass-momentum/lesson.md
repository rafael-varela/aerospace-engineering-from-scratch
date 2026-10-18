# P08.03  Control Volume: Mass & Momentum

> Draw a box, count what crosses its faces, and Newton's second law becomes bookkeeping: the force on whatever is inside equals the momentum that leaves minus the momentum that arrives, plus the pressure on the faces.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the integral mass and linear-momentum balances for a fixed or uniformly moving control volume from the Reynolds transport theorem, and state the stress model they use (pressure only on the faces the fluid crosses). 2. Derive the jet-engine thrust equation F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e from a control volume of arbitrary size, and the rocket thrust from one moving with the vehicle. 3. Implement and test a control-volume tool (a list of faces with uniform properties) and show that its mass balance closes to round-off and its momentum balance returns the formula's thrust whatever the box size and in any inertial frame. 4. Compute momentum-flux correction factors for non-uniform profiles, integrate the rocket equation with gravity loss, and analyse a moving vane. 5. Decide where to draw a control volume so that the unknown force appears once and every flux is known. |
| **Prerequisites** | `P04.03` (momentum and impulse; the rocket equation by the particle route), `P04.01` (Newton's laws, free-body diagrams), `P01.12` (the divergence theorem, surface integrals), `P03.04` (Gauss–Legendre quadrature), `P03.05` (RK4), `P08.01` (pressure; ISA at 11 km), `P08.02` (Lagrangian and Eulerian views). Skill check: you can write the net force on a body from a free-body diagram and integrate ∫u² 2πr dr for a polynomial profile. |
| **Engineering question** | How much thrust does an illustrative small turbojet produce on the test stand and in cruise at 11 km and 236 m/s, how much of it is pressure thrust, and does the answer depend on where the control volume is drawn? The answers: 12 240 N on the stand; 6030 N in cruise (gross momentum thrust 6630 N, ram drag 2360 N, pressure thrust 1760 N, 29 % of the net); and boxes of 0.5 to 100 m² all return 6030.277 N, within 5 × 10⁻¹⁰ N of the formula, with mass residuals of 7 × 10⁻¹⁶ kg/s. |
| **Assumptions and limits** | Steady flow in the chosen frame (fixed or moving at constant velocity); uniform properties on each face, or a known profile with a correction factor; only pressure acts on the faces the fluid crosses (viscous normal stresses neglected, faces drawn where the flow is not shearing strongly); the exhaust is an ideal gas with R of air for continuity at the exit; gravity on the fluid inside negligible for horizontal balances. Fails when the flow inside is unsteady (startup, surge), when the faces cut a mixing layer, and when the exit plane is not where the exhaust is uniform. Engine and vehicle numbers are illustrative, not data of a real engine. |
| **Analytical work** | The Reynolds transport theorem and the two balances; the box around a turbojet with its side face; the thrust equation and its independence of box size; Galilean invariance as a consequence of mass conservation; momentum-flux correction factors of the laminar and power-law profiles in closed form; the rocket thrust and Tsiolkovsky's equation with gravity loss; the moving vane and its optimum speed; the thrust reverser. |
| **Computational deliverable** | `starter/main.py`: `mass_outflow`, `momentum_outflow_x`, `pressure_force_x`, `solid_force_x`, `exit_area`, `turbojet_thrust`, `turbojet_control_volume`, `profile_fluxes`, `beta_laminar`, `beta_power_law`, `rocket_thrust`, `rocket_burn`, `tsiolkovsky_delta_v`, `vane_force`, `thrust_reverser_fraction` (`Face`, `observed_order`, `engine_case`, the `STATIC` and `CRUISE` cases and `main` are given). Tests prove: a uniform pressure exerts no net force on a closed surface; a straight duct carries no net flux or force; boxes of 0.3 to 1000 m² balance mass to 10⁻¹³ and return the formula's thrust to 10⁻¹²; continuity at the exit; the thrust is unchanged in a moving frame and changes when the side face is left out; the breakdown and its limits; dF/dp_e = A_e; momentum-flux factors 1, 4/3 and the power law to 10⁻⁵, never below 1 for seeded random profiles; vacuum minus sea-level rocket thrust is p_aA_e; the RK4 burn matches Tsiolkovsky with gravity loss at order 4.0 ± 0.2 and reduces to c ln(m₀/m_f) without gravity; vane forces and the V/3 optimum; reverser fractions. |
| **Practical activity** | **computational experiment**: control-volume size and frame independence, thrust against flight speed with ram drag and propulsive efficiency, the effect of nozzle expansion on cruise thrust, rocket thrust against altitude, and the gravity-loss trade of burn time. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. No home experiment with jets, rockets or pressurised gas: those are outside the course's safety rules (`P00.13`). |
| **Evidence** | Completed `main.py`, passing tests, `jet_thrust.png`, the derivations of the thrust equation from the box and of Galilean invariance, the five experiment tables, and a half-page note giving the cruise thrust breakdown and what the box-size check proves. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the mass and momentum balances are derived from the transport theorem and applied with every face accounted for; all public tests pass; the box-size and frame independence are shown and explained; the note states the thrust breakdown and the stress assumption on the faces. |
| **References** | `@mit-ocw-2-25-advanced-fluid-mechanics` Sonin, "Fundamental Laws of Motion for Particles, Material Volumes, and Control Volumes": 2 The transformation to control volumes (2.4 Reynolds' material-volume to control-volume transformation); 3 Basic laws for control volumes (3.1 Mass conservation; 3.2 Linear momentum theorem), and "On Choosing and Using Control Volumes", https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/resources/mit2_25f13_fundam_law-son/ (required-free: the transport theorem and the control-volume laws, carefully). `@leishman-intro-aerospace-flight-vehicles` chapters "Equations of Fluid Motion" (Finite Control Volume Approach; Reynolds Transport Theorem), "Momentum Equations" (Momentum Equation – Integral Approach) and "Applications of the Conservation Laws" (Force on a Pipe Bend; Jet Engine Performance) (required-free; link only). `@libretexts-bar-meir-fluid-mechanics` Ch. 5: 5.4 Reynolds Transport Theorem; Ch. 6 "Momentum Conservation for Control Volume": 6.1 Momentum Governing Equation; 6.2 Momentum Equation Application (supplementary-free). `@mit-ocw-16-01-unified-engineering-fluids` F6 "Control Volume, Mass Conservation", F7 "Momentum Conservation", F9 "Momentum Integral Applications" (supplementary-free). `@mit-ocw-1-060-engineering-mechanics-ii` L9 "Conservation Laws for Finite Control Volumes, Reynolds Transport Theorem", L10 and L11 (supplementary-free). `@nasa-glenn-beginners-guide-aeronautics` "Conservation of Mass" (supplementary-free). |
| **Extensions** | ★ Add a bypass stream to the box (a turbofan with two exits) and show that the thrust is the sum of the two streams' contributions; compute the propulsive efficiency for bypass ratios 0–10. ★ Write the angular-momentum balance for the box and compute the torque on a turbine stage from the swirl of the flow. ★ Make the box unsteady: an engine spooling up, with the mass stored inside rising; show what the steady formula misses. |
| **Time** | unmeasured |

## 1. The engineering problem

A turbojet swallows air at the flight speed, burns a little fuel in it, and throws it out of the nozzle faster than it came in. The thrust is what the engine mount must carry. You cannot measure the pressure on every blade, vane and duct wall inside the engine, and you do not need to: draw a box around the engine, account for what crosses its faces, and Newton's second law gives the force on everything inside in one line. That is the **control-volume** method, and it sizes engines (`P12`), nozzles, wind-tunnel balances, pipe bends in a propellant feed (`P08.07`), thrust reversers and rotors.

The method is powerful because it ignores the inside, and dangerous for the same reason: forget one face, or put a face where you do not know the pressure, and the answer is wrong without any sign of it. This lesson builds the method as a tool that must pass two hard checks. First, a box around an engine can be drawn any size: the thrust must not change, and mass must balance to round-off. Second, the thrust is a force, so it must be the same for an observer moving past the engine at constant speed (Galilean invariance), and that turns out to be true only because mass is conserved.

Success criteria, written first: for the illustrative cruise engine, boxes from 0.5 to 100 m² give mass residuals at round-off and the same thrust as F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e to 10⁻¹² relative; the thrust is unchanged when every x-velocity is shifted by 50 m/s; the momentum-flux factors of the laminar and 1/7th-law profiles are reproduced by quadrature to 10⁻⁵; and an RK4 integration of a rocket burn reproduces Tsiolkovsky's equation with gravity loss at fourth order.

## 2. Understand and derive

### 2.1 From a system to a control volume

Newton's second law applies to a **system**: a fixed collection of matter. In a flow the matter moves through the region we care about, so we use a **control volume** (CV), a region of space chosen by the analyst, bounded by a control surface (CS). The **Reynolds transport theorem** connects them (Sonin, 2.4; Leishman, "Equations of Fluid Motion", Reynolds Transport Theorem). For any extensive property B of the system with intensive value b = dB/dm,

**dB_sys/dt = d/dt ∫_CV ρb dV + ∮_CS ρb (u_rel·n) dA,**   (1)

with u_rel the fluid velocity relative to the moving surface and n the outward normal. The rate of change of the system's B equals the rate of change of the amount inside the CV plus the net outflow through its surface.

**Mass** (B = m, b = 1, dm_sys/dt = 0), steady:

**Σ_faces ρ u_n A = 0.**   (2)

**Linear momentum** (B = mu, b = u, dB_sys/dt = ΣF), steady, x-component:

**ΣF_x = Σ_faces ρ u_n u_x A.**   (3)

Here u_n is the normal velocity relative to the face (positive out) and u_x is the fluid's x-velocity in the inertial frame of the balance. ΣF_x is the pressure on the control surface, −Σ p n_x A, plus the forces through any solid the surface cuts (a pylon, a load cell, a pipe flange). The **stress model** is the simplest one: on faces the fluid crosses we keep only pressure, and we draw those faces where viscous stresses are negligible. A face is uniform or carries a known profile; `Face(area, nx, rho, un, ux, p)` represents one.

### 2.2 The box around a turbojet

Work in the engine's frame: air approaches at V₀ along +x (the aircraft flies towards −x). Draw a cylinder of cross-section A_box around the engine, far enough upstream that the air is undisturbed, and with its exit plane at the nozzle exit.

![A control-volume box of cross-section A_box around a turbojet in the engine frame. Air at V₀ and p_a enters through the inlet face (n_x = −1); the side face (n_x = 0) carries ṁ_side with u_x = V₀; the exit plane at the nozzle exit carries the jet (V_e, p_e, ρ_e, area A_e) and ambient flow over area A_box − A_e. The cut pylon carries force F and the fuel line brings in f ṁ_a with no x-velocity.](figures/turbojet-control-volume.svg)
*Figure 1. Control volume around a turbojet, in the engine frame. The side face carries mass although its normal has no x-component; leave it out and the mass balance fails and the thrust depends on the box size.*

**Mass** (2), face by face: inlet −ρ₀V₀A_box; jet ṁ_e = (1 + f)ṁ_a; outer exit ρ₀V₀(A_box − A_e); fuel −fṁ_a. The side face must carry the rest:

**ṁ_side = ρ₀V₀A_e − ṁ_a.**   (4)

It is the difference between the free-stream air that would have crossed the jet's area and the air the engine actually swallowed. At the lesson's cruise point ρ₀V₀A_e = 8.71 kg/s and ṁ_a = 10 kg/s, so ṁ_side = −1.295 kg/s: air is drawn *in* through the side of the box towards the intake.

**Momentum** (3) with the pylon force F and the pressures: inlet −ρ₀V₀²A_box + p_aA_box; jet ṁ_eV_e + p_eA_e; outer exit ρ₀V₀²(A_box − A_e) + p_a(A_box − A_e); side ṁ_sideV₀ (its normal has no x-component, so no pressure term); fuel zero x-momentum. Adding, and using (4),

**F = ṁ_eV_e − ṁ_aV₀ + (p_e − p_a)A_e**   (5)

(Leishman, "Applications of the Conservation Laws", Jet Engine Performance). A_box has cancelled: the box can be any size. The three terms are the **gross momentum thrust**, the **ram drag** (the momentum the engine must give the swallowed air just to bring it along with the aircraft) and the **pressure thrust** (non-zero when the nozzle does not expand the exhaust to ambient pressure). The exit area follows from continuity, A_e = ṁ_e/(ρ_eV_e) with ρ_e = p_e/(RT_e).

**The side face is not optional.** Leave it out and (2) fails by ṁ_side; then (3) gives a box-size-dependent answer. The cleanest way to see this is to change frame.

### 2.3 Galilean invariance needs mass conservation

An observer moving at −W along x sees every x-velocity increased by W; the normal velocities relative to the faces do not change. The right-hand side of (3) becomes Σρu_n(u_x + W)A = Σρu_nu_xA + WΣρu_nA. The thrust is a force and must not depend on the observer, so

**the momentum balance is frame-independent if and only if Σρu_nA = 0.**   (6)

The test `test_thrust_is_galilean_invariant_because_mass_balances` shifts the frame by 50 m/s: the full box returns the same force to 2 × 10⁻¹¹ N; the box without its side face (mass residual 1.295 kg/s) changes by 64.8 N = 1.295 × 50. In the ground frame the engine moves and air is at rest far ahead; in the engine frame the air moves; (6) is why both give (5).

### 2.4 Non-uniform faces

Real exit profiles are not uniform. With a profile u(r) on a circular face of radius R, the mass flow is ∫ρu 2πr dr and the momentum flow ∫ρu² 2πr dr. Writing the momentum flow as βṁV̄, with V̄ the mean velocity, defines the **momentum-flux correction factor**

**β = (1/A)∫(u/V̄)² dA ≥ 1**   (7)

(≥ 1 by the Cauchy–Schwarz inequality). For the laminar pipe profile u = u_max(1 − r²/R²), β = 4/3; for the power law u = u_max(1 − r/R)^{1/n},

**β = (n + 1)(2n + 1)²/(4n²(n + 2)),**   (8)

1.0204 for the turbulent 1/7th law. Using ṁV̄ for a turbulent jet underestimates the gross thrust by 2 %; for a laminar one by 25 % (a factor 3/4).

### 2.5 A control volume that moves: rockets and vanes

**Rocket.** Take a CV moving with the vehicle (an inertial frame momentarily, `P04.03`). No air enters; propellant leaves at V_e relative to the nozzle. With ambient pressure p_a on the whole surface except the exit,

**F = ṁV_e + (p_e − p_a)A_e.**   (9)

The thrust grows as p_a falls; the vacuum thrust exceeds the sea-level thrust by p_aA_e. Newton's law for the vehicle, m dv/dt = F − mg with m = m₀ − ṁt, integrates for constant F and vertical flight to

**Δv = c ln(m₀/m_f) − g t_b,  c = F/ṁ,**   (10)

Tsiolkovsky's equation with the **gravity loss** g t_b.

**Moving vane.** A jet of area A and speed V strikes a vane moving downstream at u and is turned through θ. In a CV moving with the vane the flow is steady, arrives at V − u, and (inviscid, at constant pressure) leaves at V − u in the turned direction. Mass flow through the vane's CV is ρA(V − u), so

**F_x = ρA(V − u)²(1 − cos θ),  power = F_x u,**   (11)

which peaks at u = V/3 for a single vane. A **thrust reverser** that turns the exhaust through θ at the same speed delivers −cos θ of the gross thrust forward: 0.707 at 135°.

### 2.6 Where the model stops

- **Faces in shear.** A face drawn across a mixing layer or a boundary layer carries viscous stresses the model drops; draw it where the flow is uniform.
- **Unsteadiness.** (2) and (3) above drop d/dt ∫ρ dV and d/dt ∫ρu dV: an engine spooling up stores mass and momentum (Extension 3).
- **The exit state.** (5) needs ṁ_e, V_e, p_e and A_e consistent with each other (continuity) and with the nozzle; the formula does not know whether the nozzle can produce them.

## 3. Build the computational model

Open `starter/main.py`. `Face`, `observed_order`, `engine_case`, the `STATIC` and `CRUISE` dictionaries and `main` are given.

### Step 1: face sums

`mass_outflow`, `momentum_outflow_x` and `pressure_force_x` are one-line sums over faces of ρu_nA, ρu_nu_xA and −pn_xA. `solid_force_x` is the force the cut solids must carry, (3) rearranged: momentum outflow minus pressure force. Check that a closed surface at uniform pressure gives zero.

### Step 2: the engine

`exit_area` is continuity; `turbojet_thrust` is (5) with its breakdown. `turbojet_control_volume` builds the five faces of 2.2, with the side face's flow from (4):

```python
m_side = rho_a*V0*A_box + f*mdot_a - (1 + f)*mdot_a - rho_a*V0*(A_box - Ae)
faces = [Face(A_box, -1, rho_a, -V0, V0, pa), Face(Ae, 1, rho_e, Ve, Ve, pe),
         Face(A_box - Ae, 1, rho_a, V0, V0, pa), Face(1, 0, 1, m_side, V0, pa), Face(1, 0, 1, -f*mdot_a, 0, pa)]
```

(The side and fuel faces use unit area and density so that u_n carries the mass flow.)

### Step 3: profiles, rockets, vanes

`profile_fluxes` maps Gauss–Legendre points (`np.polynomial.legendre.leggauss`) to [0, R] and sums ρu 2πr w and ρu² 2πr w; `beta_laminar` and `beta_power_law` are 4/3 and (8). `rocket_thrust` is (9); `rocket_burn` integrates dv/dt = F/(m₀ − ṁt) − g, dh/dt = v by RK4; `tsiolkovsky_delta_v` is (10); `vane_force` and `thrust_reverser_fraction` are 2.5.

Running `main()` prints

```
illustrative turbojet (numbers chosen for the lesson, not a specific engine)
  static, sea level        A_e = 0.0771 m^2  gross  12240.0 N  ram drag     0.0 N  pressure     0.0 N  net F =  12240.0 N
  cruise, 11 km, 236 m/s   A_e = 0.1014 m^2  gross   6630.0 N  ram drag  2360.0 N  pressure  1760.3 N  net F =   6030.3 N

control-volume box around the cruise engine, any size:
   A_box (m^2)   mass residual (kg/s)   F from box (N)   F - formula (N)
         0.5     -7.22e-16            6030.277      +1.8e-12
         2.0     -7.22e-16            6030.277      +0.0e+00
        10.0     -7.22e-16            6030.277      +2.7e-11
       100.0     -7.22e-16            6030.277      +5.2e-10

momentum-flux correction factors (quadrature vs closed form):
  laminar parabola: beta = 1.333333 (4/3 = 1.333333)
  power law 1/6: beta = 1.026908 (closed form 1.026910; ...)
  power law 1/7: beta = 1.020407 (closed form 1.020408; ...)
  power law 1/10: beta = 1.010624 (closed form 1.010625; ...)

sounding rocket (illustrative): m0 = 100 kg, 60 kg propellant, 8 s burn, c = 2000 m/s
  thrust 15000 N, initial thrust/weight 15.30
  burnout speed: Tsiolkovsky - g t_b = 1754.128 m/s (gravity loss 78.45 m/s); RK4 n = 80: 1754.128264 m/s
  RK4 errors [2.03e-03 1.28e-04 8.02e-06 5.02e-07], orders [3.99 4.   4.  ]
  burnout altitude 5912.4 m

vane turning a water jet (A = 1 cm^2, V = 20 m/s) by 180 deg: power peaks at u = 6.70 m/s (V/3 = 6.67), 236.6 W
thrust reverser turning the jet by 135 deg: 0.707 of gross thrust reversed
```

and writes `jet_thrust.png`: net thrust and ram drag against flight speed, and the rocket's speed during the burn by RK4 against (10).

### Verify

Run `python -m pytest curriculum/phase-08/lesson-03-control-volume-mass-momentum`. The 22 tests:

1. **Face sums**: a closed surface at uniform pressure has no net force; a straight duct no net mass, momentum or force.
2. **The box** (four sizes, 0.3–1000 m²): mass residual below 10⁻¹³ of the inlet flow and the formula's thrust to 10⁻¹²; continuity at the exit to 10⁻¹⁴; a box smaller than the jet is refused.
3. **Frames**: a 50 m/s shift leaves the force unchanged to 10⁻¹²; without the side face the mass residual exceeds 1 kg/s and the shift changes the force by more than 10 N.
4. **Breakdown and limits**: F = gross − ram + pressure; positive pressure thrust when underexpanded; static and perfectly expanded: F = (1 + f)ṁV_e; no change in velocity or pressure: no thrust; dF/dp_e = A_e.
5. **Profiles**: uniform β = 1 and mass flow exact; laminar mean u_max/2 and β = 4/3 to 10⁻¹²; power law 1/6, 1/7, 1/9 to 10⁻⁵ and the flat-profile limit; twenty seeded random positive profiles all have β ≥ 1.
6. **Rockets**: vacuum minus sea-level thrust = p_aA_e; RK4 burn against (10), order 3.8–4.2, climbing throughout; without gravity Δv = c ln(m₀/m_f) to 10⁻⁹.
7. **Vane**: 2ρAV² at 180°, zero at 0°, (ρAV², −ρAV²) at 90°, and the power peak at V/3; reverser fractions at 90°, 135°, 180°.
8. `main` writes the figure.

### Compare

Two independent routes agree for the thrust: the formula (5) and the face-by-face momentum balance of a box that knows nothing of the formula. Their difference, 0 to 5 × 10⁻¹⁰ N on 6030 N, grows with the box size only because the inlet and outer-exit fluxes (each about 2 × 10⁵ N for a 10 m² box, 2 × 10⁶ N at 100 m²) cancel in floating point: a relative round-off of 10⁻¹⁶ on the largest term. For the rocket, RK4 and (10) agree to 5 × 10⁻⁷ m/s at 80 steps, converging at order 4.00. The momentum-flux factors from 200-point Gauss–Legendre quadrature match the closed forms to 1.5 × 10⁻⁶; quadrature converges only algebraically for the power law (errors 1.3 × 10⁻⁴, 3.0 × 10⁻⁵, 6.7 × 10⁻⁶, 1.5 × 10⁻⁶ at 25, 50, 100, 200 points for 1/7), because the profile's slope is infinite at the wall, and to round-off for the laminar polynomial.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Box size and frame.** The table above for box size; then shift the frame by 50 m/s with and without the side face. *Finding:* full box: force change 1.8 × 10⁻¹¹ N; without the side face (residual 1.2954 kg/s): 64.77 N, exactly 1.2954 × 50. A control volume that does not conserve mass gives an observer-dependent "force".
2. **Ram drag and flight speed.** The static engine's ṁ_a = 20 kg/s and V_e = 600 m/s held fixed (perfectly expanded), V₀ = 0–300 m/s. Predict the fraction of static thrust left at 236 m/s and the propulsive efficiency η_p = 2/(1 + V_e/V₀). *Finding:* F = 12 240, 11 040, 10 240, 8240, 7520, 6240 N at V₀ = 0, 60, 100, 200, 236, 300 m/s (61 % left at 236 m/s); η_p = 0.18, 0.29, 0.50, 0.56, 0.67 at 60–300 m/s. Thrust falls as efficiency rises: the reason fans and propellers, which accelerate a lot of air a little, are efficient and jets, which accelerate a little air a lot, are not.
3. **Nozzle expansion.** At the cruise point hold the gas state at the nozzle (the reference 40 kPa, 900 K, 650 m/s) and let a longer nozzle expand it isentropically (γ = 1.33, c_p = 1150 J/(kg·K), typical hot-gas values) to p_e = 30, 22.632, 18, 15 kPa; recompute V_e, A_e and F. Predict where F peaks. *Finding:* V_e = 650, 751.7, 833.8, 891.4, 932.5 m/s; A_e = 0.1014, 0.1088, 0.1212, 0.1347, 0.1477 m²; F = 6030.3, 6109.3, **6145.1**, 6108.0, 6024.8 N. The peak is at p_e = p_a: pressure thrust is zero there, but the momentum thrust gained by expanding is larger than the pressure thrust given up. Over-expanding (p_e < p_a) loses again.
4. **Rocket thrust with altitude.** The rocket engine of Problem 1 (ṁ = 12 kg/s, V_e = 2800 m/s, p_e = 40 kPa, A_e = 0.065 m²) at 0, 3000, 7185 and 11 000 m. *Finding:* 29 614, 31 643, 33 600, 34 729 N, against 36 200 N in vacuum; the pressure thrust changes sign at 7185 m, where p_a = p_e.
5. **Burn time and gravity loss.** The 100 kg sounding rocket with 60 kg of propellant and c = 2000 m/s, burning in 4, 8, 16, 32 and 60 s. *Finding:* burnout speeds 1793.4, 1754.1, 1675.7, 1518.8, 1244.2 m/s (RK4 and (10) agree to 10⁻⁹); initial thrust-to-weight 30.6 down to 2.0. A short burn wins in vacuum-like conditions; in the real atmosphere drag (`P08.16`) pushes the other way.

Raw data format: CSV tables `A_box_m2, mass_residual, F_box_N, F_formula_N` and `frame_shift, with_side_face, dF_N` (1); `V0_ms, F_N, eta_p` (2); `pe_Pa, Ve_ms, Ae_m2, F_N` (3); `h_m, pa_Pa, F_N` (4); `t_b_s, F_over_W0, v_bo_RK4, v_bo_closed, h_bo_m` (5); plus `jet_thrust.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| cruise thrust | (5): 6030.277 N | box, 0.5–100 m² | ≤ 5 × 10⁻¹⁰ N | round-off on large cancelling fluxes |
| mass balance | (2): 0 | −7.2 × 10⁻¹⁶ kg/s | round-off | — |
| frame shift 50 m/s | (6): no change | 1.8 × 10⁻¹¹ N | round-off | — |
| side-face omission | W × ṁ_side = 64.77 N | 64.77 N | 0 | (6) |
| β laminar, 1/7th | 4/3, 1.020408 | 1.333333, 1.020407 | < 1.5 × 10⁻⁶ | quadrature at the wall singularity |
| rocket burnout speed | (10): 1754.128 m/s | RK4, 80 steps | 5 × 10⁻⁷ m/s | RK4, order 4 |
| vane optimum | V/3 = 6.67 m/s | 6.70 m/s (grid of 0.1 m/s) | grid | — |

The engineering answer: **the illustrative turbojet produces 12 240 N on the stand and 6030 N in cruise at 11 km and 236 m/s, of which 1760 N (29 %) is pressure thrust from its underexpanded nozzle and 2360 N is lost to ram drag; any control volume, of any size, in any inertial frame, gives the same number, because mass balances.** A nozzle that expanded to ambient pressure would add 1.9 % (Experiment 3).

What to improve: the exit state. The thrust is only as good as ṁ_e, V_e and p_e; a real analysis takes them from a nozzle model (`P08.17`, `P12`). What not to improve: the box. Its job is to prove the bookkeeping, and it does so to round-off.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (4), (5) and (6) from (1); `jet_thrust.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the cruise thrust breakdown, what the box-size and frame checks prove, and which face assumption you would worry about in a real installation.

### Adversarial self-checks

- **What if the flight speed goes to zero?** Ram drag vanishes and the static thrust is (1 + f)ṁV_e + (p_e − p_a)A_e; with the lesson's static engine, 12 240 N. What if V₀ → V_e? The net thrust tends to fṁV_e plus any pressure thrust: an engine cannot push an aircraft faster than its own jet.
- **What if a face is forgotten?** The mass residual is non-zero and the "thrust" depends on the box size and the observer's frame: Experiment 1 shows a 64.8 N change for a 50 m/s shift. Always print the mass residual.
- **Which assumption fails first in the real world?** The exit plane: a real exhaust is not uniform (β ≈ 1.02 for a turbulent jet, Problem 5) and its pressure varies across it; the formula's inputs carry the uncertainty, not the method.
- **What measurement could falsify the analysis?** A test stand's load cell against (5) computed from measured mass flow, exit pressure and a traversed exit velocity profile; a systematic gap of 2 % would point at β, a gap that varies with power setting at p_e.
- **Does the result survive a change of control volume?** Yes, by construction and by test: four box sizes and two frames.
- **Which decision changes if the uncertainty doubles?** If the exit velocity is known to ±2 % instead of ±1 %, the gross thrust is known to ±2 % (130 N in cruise) but the net thrust to ±2.2 % (because the ram drag is subtracted exactly): still within a typical engine-sizing margin, but no longer good enough to resolve the 1.9 % gain of a better nozzle.

## Check yourself

1. What is the difference between a control volume and a system, and which one do Newton's laws apply to?
2. State the Reynolds transport theorem in words.
3. Write the thrust equation of a jet engine and name its three terms.
4. Why does the thrust from a box around the engine not depend on the box's size?
5. What is the momentum-flux correction factor, and what are its values for laminar and turbulent pipe profiles?
6. How does a control volume moving with a rocket give the rocket equation, and where does the gravity loss come from?
7. A jet strikes a single moving vane. At what vane speed is the power greatest, and why?
8. Why does a jet engine's net thrust fall as flight speed rises, even with the same exhaust speed?
