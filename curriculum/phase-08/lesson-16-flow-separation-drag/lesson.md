# P08.16  Flow Separation & Drag

> A bluff body's drag is mostly the pressure it fails to recover behind it. Keep the flow attached, by shape or by a turbulent boundary layer, and the drag falls by an order of magnitude.

| | |
|---|---|
| **Learning outcomes** | 1. Split drag into pressure (form) drag and friction drag, compute pressure drag by integrating a surface pressure distribution, and show why attached potential flow gives none (d'Alembert's paradox) while a separated wake with a low base pressure gives most of a cylinder's drag. 2. Use the Stokes law and the Schiller–Naumann and Morrison correlations for c_d(Re) of a sphere, and White's for a cylinder, labelled as correlations, and explain the drag crisis. 3. Predict vortex-shedding frequencies from Strouhal-number correlations and check them against a structure's natural frequency. 4. Implement and test a terminal-velocity solver (root finding on the force balance), checked against Stokes's closed form, a time integration of the equation of motion and Reynolds similarity. 5. Decide whether to fair a strut, and estimate what the fairing saves. |
| **Prerequisites** | `P08.08` (potential flow round a cylinder, c_p = 1 − 4 sin²θ, d'Alembert's paradox), `P08.10`–`P08.11` (laminar boundary layers, adverse pressure gradients and separation), `P08.12` (transition), `P08.14` (turbulent skin friction, C_f = 0.074 Re^(−1/5)), `P08.05` (Reynolds number and similarity), `P08.01` (the ISA and Sutherland's law), `P03.02` (root finding), `P03.04` (quadrature), `P03.06` (ODE integrators), `P04.01` (Newton's second law). Skill check: you can write the drag force ½ρV²c_dA and say what A is for a sphere and for a cylinder. |
| **Engineering question** | One metre of 25 mm round tube is exposed as a strut on the light aeroplane at 60 m/s TAS at 3000 m. What drag and power does it cost, at what frequency does it shed vortices, and how much does a 4:1 fairing of the same thickness save? The answers: Re_d = 80 500, c_d ≈ 1.0 from White's correlation (published values near 1.2 make it up to 20 % more), so 41–49 N and 2.5–2.9 kW of propulsive power per metre; shedding at about 509 Hz (St ≈ 0.21); a fairing of thickness ratio 0.25 has c_d ≈ 0.081 on its thickness, 3.3 N and 0.2 kW, twelve times less: the drag of a 2 mm round wire. |
| **Assumptions and limits** | Steady, incompressible flow (M = 0.18 here) round smooth bodies in a uniform, low-turbulence stream; two-dimensional cylinders (no end effects; a short strut with free ends has less drag); correlations are fits to experiments with scatter of order ±10 % and stated Re ranges, outside which they must not be used; the drag crisis depends strongly on surface roughness and free-stream turbulence, which the smooth-sphere fits do not include. Terminal velocity: a rigid sphere falling steadily in still fluid (no wall effects, no added mass or history force in the transient, no spin). The fairing estimate is a classical form-factor correlation for attached flow. |
| **Analytical work** | Pressure drag as ∮(p − p∞) n dS and c_d = ½∫c_p cos θ dθ; d'Alembert's paradox from the symmetry of the potential-flow c_p; a closed-form separated-wake model c_d = sin θ_s − (4/3)sin³θ_s − c_pb sin θ_s; Stokes drag 24/Re and terminal velocity (ρ_p − ρ_f)gd²/(18μ); the force balance for terminal velocity at any Re and its Archimedes-number similarity; the shedding frequency f = St V/d; the fairing's friction-dominated drag. |
| **Computational deliverable** | `starter/main.py`: `cp_potential_cylinder`, `cp_separated_model`, `pressure_drag_coefficient`, `pressure_lift_coefficient`, `cd_sphere_stokes`, `cd_sphere_schiller_naumann`, `cd_sphere_morrison`, `cd_cylinder_white`, `drag_crisis`, `strouhal_cylinder`, `shedding_frequency`, `stokes_terminal_velocity`, `terminal_velocity`, `simulate_fall`, `fairing_cd`, `strut_drag` (`isa_troposphere`, `reynolds`, `flat_plate_cf_turbulent`, the fluid constants and `main` are given). Tests (12, under 2 s) prove: attached potential flow has zero drag and side force; the drag integral is exact for Fourier modes and reproduces the separated model's closed form; the sphere correlations tend to 24/Re, sit in the Newton range and (Morrison only) show the drag crisis between Re = 2 × 10⁵ and 6 × 10⁵; St is undefined below Re = 47, continuous, and near 0.21 when subcritical; f scales as V/d; terminal velocity balances weight and drag to 10⁻¹⁰, equals Stokes's law with 24/Re, and is approached monotonically by the time integration; terminal Re depends only on the Archimedes number; a fairing saves between 5 and 25 times. |
| **Practical activity** | **computational experiment**: a pressure-drag study with a separated-wake model, a comparison of three sphere correlations across nine decades of Re, terminal velocities of four objects by root finding and by time integration, and a strut-fairing trade study. |
| **Equipment** | None. A computer with the environment from `P00.01`. (The ping-pong ball case connects to the drop experiment of `projects/p04-drag-drop`.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `drag.png`, the derivations of the pressure-drag integral, the separated-wake closed form and Stokes's terminal velocity, the four experiment tables, and a half-page recommendation on fairing the strut. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: pressure drag is derived and d'Alembert's paradox explained; all public tests pass; every drag coefficient used is labelled as a correlation with its range; the recommendation states the drag saving with its uncertainty and checks the shedding frequency. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapter "Bluff Body Flows": Drag Coefficients; Flow Patterns Around a Circular Cylinder; Drag of Circular Cylinders & Spheres; Roughness Effects on Cylinders & Spheres; Why Do Golf Balls Have Dimples?; Vortex-Induced Vibration; Streamlining (required-free; https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/bluff-body-flows/); chapter "Time-Dependent Flows": Strouhal Number (supplementary-free). `@nasa-glenn-beginners-guide-aeronautics` "Drag of a Sphere"; "Shape Effects on Drag" (required-free). `@mit-ocw-1-060-engineering-mechanics-ii` L20 "Drag Forces on Immersed Bodies" (supplementary-free). `@mit-ocw-2-25-advanced-fluid-mechanics` "Boundary Layers, Separation, and Drag" (supplementary-free; its textbook readings are paid and not needed). |
| **Extensions** | ★ Replace the constant base pressure by a measured-shape model (a linear recovery after separation) and fit θ_s and c_pb to the c_d of White's correlation and to St. ★ Add a rough-sphere crisis by shifting Morrison's crisis term to lower Re, and recompute the payload of Problem 4. ★ Integrate a ping-pong ball's fall with added mass (½ρ_f V_ball) and compare with the `projects/p04-drag-drop` data from Phase 04. |
| **Time** | unmeasured |

## 1. The engineering problem

Drag costs fuel, range and climb rate, and on a light aeroplane a surprising share of it comes from small things: wheels, struts, antennas, the pitot mast, steps and hinges. Each is a **bluff body**: the flow round it cannot follow its rear surface, separates, and leaves a wide wake at low pressure. That wake, not the friction on the surface, is where most of the drag comes from. The same physics sets the terminal velocity of anything falling through air or water, from a raindrop to a sounding-rocket payload whose parachute failed, and the frequency at which a mast or a cable sheds vortices and may start to vibrate.

Potential flow (`P08.08`) predicts zero drag for any body: d'Alembert's paradox. Boundary-layer theory (`P08.10`–`P08.11`) explains why: the boundary layer cannot climb the adverse pressure gradient behind the body's shoulder and separates. Where it separates, how low the pressure in the wake is, and whether the boundary layer is laminar or turbulent when it separates, depend on the Reynolds number and the surface. There is no closed-form theory for any of that at engineering Reynolds numbers; engineers use **correlations**: curve fits to experiments, with stated ranges and scatter.

This lesson builds the tools: pressure drag from a surface pressure, the c_d(Re) correlations for spheres and cylinders with the drag crisis, vortex shedding, terminal velocity, and the effect of streamlining. It applies them to a round strut on the light aeroplane. Success criteria, written first: zero drag from attached potential flow to round-off; correlations reproduce the Stokes limit; terminal velocities balance weight and drag to 10⁻¹⁰ and agree with a time integration; every drag coefficient is labelled as a correlation and used inside its range.

## 2. Understand and derive

### 2.1 Pressure drag and friction drag

The force of the fluid on a body is the integral over its surface of the pressure, acting along the inward normal, and of the wall shear stress, acting along the surface. Their components along the stream are the **pressure drag** (or form drag) and the **friction drag**. A flat plate along the stream has only friction drag; the same plate across the stream has only pressure drag. A streamlined body is one shaped so that the pressure drag is small and friction dominates.

For a circular cylinder of diameter d, measure θ from the front stagnation point. A surface element at angle θ has outward normal with stream-wise component −cos θ, so the pressure pushes it downstream with (p − p∞) cos θ (d/2) dθ per unit span. Dividing by the dynamic pressure q = ½ρV² and the frontal area per unit span d:

**c_d,p = ½ ∫₀^{2π} c_p cos θ dθ,  c_p = (p − p∞)/q.**   (1)

A side force follows from −sin θ in place of cos θ. For a periodic, smooth c_p the trapezoid rule on n equal intervals is exact for every Fourier mode below n/2, which makes (1) easy to evaluate accurately.

### 2.2 d'Alembert's paradox

Attached potential flow round a cylinder has c_p = 1 − 4 sin²θ (`P08.08`). It is symmetric front to back (θ → π − θ), while cos θ is antisymmetric, so (1) is zero: **no drag**. Physically, the pressure recovers fully on the rear half and pushes the body forward exactly as hard as the front half pushes it back. A real flow cannot recover that pressure: the slow fluid near the wall in the boundary layer cannot climb from c_p = −3 at the shoulder back to +1 at the rear, and separates.

### 2.3 A separated-wake model

A simple model keeps the potential-flow c_p ahead of a separation angle θ_s on each side and a constant base pressure c_pb in the wake behind. Integrating (1) piecewise:

front: ½∫_{−θs}^{θs}(1 − 4 sin²θ) cos θ dθ = sin θ_s − (4/3) sin³θ_s;  wake: ½ c_pb ∫_{θs}^{2π−θs} cos θ dθ = −c_pb sin θ_s.

**c_d,p = sin θ_s − (4/3) sin³θ_s − c_pb sin θ_s.**   (2)

With laminar separation near θ_s = 82° and c_pb = −1.0 to −1.2 (typical subcritical values), (2) gives 0.69–0.88. The wake term −c_pb sin θ_s is most of it: the low base pressure *is* the drag. When the boundary layer becomes turbulent before separating, it carries more momentum near the wall, separates later (θ_s ≈ 110–120°) and the base pressure rises (c_pb ≈ −0.2 to −0.3): (2) then gives 0.12–0.17. That is the **drag crisis**.

### 2.4 Drag correlations

Write the drag of a sphere as D = c_d ½ρV² (πd²/4) and Re = ρVd/μ.

- **Stokes's law** (Re ≪ 1): c_d = 24/Re, D = 3πμVd. Exact for creeping flow; one third of it is pressure drag and two thirds friction.
- **Schiller–Naumann** (a correlation, Re ≤ 1000): c_d = (24/Re)(1 + 0.15 Re^0.687), and the Newton-regime constant 0.44 above 1000. It has no drag crisis.
- **Morrison** (2013, a correlation for smooth spheres up to Re = 10⁶):
  c_d = 24/Re + 2.6(Re/5)/(1 + (Re/5)^1.52) + 0.411(Re/2.63 × 10⁵)^(−7.94)/(1 + (Re/2.63 × 10⁵)^(−8)) + 0.25(Re/10⁶)/(1 + Re/10⁶).
  The third term switches the drag down through the crisis; its minimum, 0.090, falls at Re = 4.4 × 10⁵.
- **Cylinder, White's fit** (a correlation, 1 < Re < 2 × 10⁵): c_d = 1 + 10 Re^(−2/3), per unit span on area d. It gives 1.0 near Re = 10⁵; many published curves show about 1.2 there, so treat it as ±20 %.

These are fits, with scatter of order ±10 % against the experiments behind them. Near Re = 1 the two sphere fits differ by about 12 %, near 1000 by about 10 %. `@nasa-glenn-beginners-guide-aeronautics` "Drag of a Sphere" and `@leishman-intro-aerospace-flight-vehicles` "Drag of Circular Cylinders & Spheres" show the measured curves these fits follow; look at them rather than at any single number.

### 2.5 Vortex shedding

Above Re ≈ 47 the wake of a cylinder is unsteady: vortices are shed alternately from each side at a frequency f, and the **Strouhal number** St = fd/V is nearly constant over a wide range. Roshko's correlations: St = 0.212(1 − 21.2/Re) for 50 ≤ Re ≤ 150 and St = 0.212(1 − 12.7/Re) above 300 (fitted to 2000; St stays at 0.20–0.21 up to about 2 × 10⁵), with an irregular range in between. The alternating vortices push the body sideways at f; if f approaches a natural frequency of the structure, it can lock in and vibrate (vortex-induced vibration, `@leishman-intro-aerospace-flight-vehicles` "Vortex-Induced Vibration").

### 2.6 Terminal velocity

A sphere of density ρ_p falling at constant speed in a fluid of density ρ_f and viscosity μ has weight minus buoyancy equal to drag:

(ρ_p − ρ_f) g (πd³/6) = ½ ρ_f V² c_d(Re) (πd²/4),  Re = ρ_f V d/μ.   (3)

With c_d = 24/Re, (3) is linear in V: **V_Stokes = (ρ_p − ρ_f) g d²/(18μ)**. In general c_d depends on V, so (3) is a scalar root-finding problem (`P03.02`): the residual is negative at small V and positive at large V for a monotone drag, and `brentq` brackets it between 10⁻⁶ and 10³ times the Stokes speed. Multiply (3) by ρ_f d²/μ² and it becomes **c_d Re² = (4/3) Ar**, with the Archimedes number Ar = (ρ_p − ρ_f)ρ_f g d³/μ²: the terminal Re depends only on Ar. That is Reynolds similarity for falling objects, and the tests use it.

The time integration of m dV/dt = (m − m_f)g − D from rest (`P03.06`, `P04.01`) must approach the root of (3); its initial acceleration is g(1 − ρ_f/ρ_p).

### 2.7 Streamlining

A fairing of thickness t and chord c round a strut keeps the flow attached; its drag is mostly friction on both surfaces, increased by a form factor for the thickness. A classical correlation for 2-D sections is c_d,c = 2C_f(1 + 2t/c + 60(t/c)⁴) on the chord, with C_f the turbulent flat-plate coefficient (`P08.14`). On the thickness, c_d,t = c_d,c · c/t. Thinner sections have less form drag but more wetted area per unit thickness; the optimum is near t/c = 0.25 (Problem 5 of the problem set explores it).

## 3. Build the computational model

Open `starter/main.py`. `isa_troposphere`, `reynolds`, `flat_plate_cf_turbulent`, the fluid constants and `main` are given.

### Step 1: pressure drag

`cp_potential_cylinder`, `cp_separated_model` (fold θ into [0, π] and switch at θ_s with `numpy.where`), and `pressure_drag_coefficient` (1) on a uniform periodic grid: `0.5 * np.sum(cp * np.cos(theta)) * 2 * np.pi / n`.

### Step 2: correlations and shedding

The four c_d functions are one line each. `drag_crisis` minimises Morrison's c_d over log₁₀Re with `minimize_scalar(..., method="bounded")`. `strouhal_cylinder` returns NaN below 47 and interpolates across 150–300.

### Step 3: terminal velocity

```python
W = (rho_p - rho_f) * G0 * np.pi * d**3 / 6
f = lambda V: 0.5 * rho_f * V**2 * cd(rho_f * V * d / mu) * np.pi * d**2 / 4 - W
V = brentq(f, 1e-6 * Vs, 1e3 * Vs)
```

`simulate_fall` integrates the equation of motion with `solve_ivp` (RK45, rtol 10⁻⁹).

### Step 4: the strut

`fairing_cd` is 2.7; `strut_drag` evaluates both struts at the flight condition, with the ISA density and Sutherland viscosity at 3000 m (`P08.01`).

### Verify

`pytest curriculum/phase-08/lesson-16-flow-separation-drag` runs 12 tests in under 2 s. Three of them check analysis against computation: the drag of attached potential flow is below 10⁻¹³; the separated model matches (2) to 10⁻⁴ at three angles (the discontinuity at θ_s limits the trapezoid rule to first order, hence 400 000 points); the time integration of a steel ball in water reaches the root of (3) to 10⁻⁴ and starts with acceleration g(1 − ρ_f/ρ_p).

### Compare

`main()` prints:

```
ISA 3000 m: T = 268.65 K, p = 70109 Pa, rho = 0.9091 kg/m^3, mu = 1.6937e-05 Pa s
  Re_d = 80514, c_d = 1.005 (White correlation), drag 41.1 N, power 2.47 kW
  shedding: St = 0.212 (Roshko correlation), f = 509 Hz
  faired, t/c = 0.25 (chord 100 mm, Re_c = 322056): c_d (on thickness) = 0.0813, drag 3.32 N, power 199 W
  fairing cuts the drag by a factor 12.4; same drag as a round wire of 2.02 mm
  attached potential flow: c_d = -1.1e-16 (d'Alembert)
  separation at 82 deg, base c_p -1.2: c_d = 0.884
  separation at 120 deg, base c_p -0.2: c_d = 0.173
  Morrison drag-crisis minimum: c_d = 0.0902 at Re = 4.44e+05; c_d at 1e5 is 4.7 times larger
  ping-pong ball, 40 mm, 2.7 g, air: V_t = 9.344 m/s (Re = 2.49e+04, c_d = 0.395); Schiller-Naumann 8.852 m/s
```

and writes `drag.png`: c_d(Re) of the three sphere correlations and the cylinder fit on log axes (titled as correlations, not data), and the separated and attached c_p(θ) of the cylinder. The separated model's 0.88 is below the cylinder correlation's 1.0 because it leaves out friction (a few per cent at this Re) and, more, because a real wake's pressure is not constant and the real separation point moves; it explains where the drag comes from, not its value to 10 %.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Where the drag comes from.** Evaluate (2) for θ_s = 80–120° and c_pb from −1.2 to −0.2. Predict which matters more, θ_s or c_pb. *Finding:* at θ_s = 80°: 0.70 (c_pb −1.0) and 0.89 (−1.2); moving separation to 85° with c_pb −1.2 changes c_d by only 2 % (0.87); a supercritical wake, 110° and −0.3, gives 0.12, and 120° and −0.2 gives 0.17. The base pressure dominates; the separation angle matters through it.
2. **Three sphere correlations.** c_d at Re = 0.1 to 10⁶. Predict where they disagree most. *Finding:* within 3.1 % of 24/Re at Re = 0.1; Schiller–Naumann is 12 % above Morrison at Re = 1, 5 % above at 10–100 and 9 % below at 1000; above 10⁵ only Morrison has the crisis, falling from 0.426 at 10⁵ to 0.093 at 4 × 10⁵ (minimum 0.090 at 4.4 × 10⁵), then rising to 0.130 at 10⁶.
3. **Terminal velocities.** Four objects, by root finding (Morrison and Schiller–Naumann) and by time integration. Predict which obeys Stokes's law. *Finding:*

   | object | V_t Morrison | Re | V_t Schiller–Naumann | Stokes law | time integration |
   |---|---|---|---|---|---|
   | ping-pong ball, 40 mm, 2.7 g, air 20 °C | 9.34 m/s | 2.5 × 10⁴ | 8.85 m/s | 3822 m/s | 9.344 m/s |
   | glass bead 0.1 mm in water | 8.01 mm/s | 0.80 | 7.29 mm/s | 8.17 mm/s | 8.01 mm/s |
   | steel ball 10 mm in water | 1.513 m/s | 1.5 × 10⁴ | 1.428 m/s | 373 m/s | 1.513 m/s |
   | water drop 1 mm in air at 3000 m | 4.48 m/s | 240 | 4.36 m/s | 32 m/s | 4.475 m/s |

   Only the bead is near Stokes (2 % off at Re = 0.8). The correlations differ by 3–10 %, which is the honest uncertainty of any terminal velocity from a correlation.
4. **The strut trade.** Round strut against fairings of t/c = 0.1–0.5 at 60 m/s, and the round strut at 30, 45 and 60 m/s. Predict the best thickness ratio and how the round strut's drag scales with V. *Finding:* the fairing's c_d on thickness is 0.118, 0.094, 0.081, 0.089, 0.155 for t/c = 0.1, 0.15, 0.25, 0.33, 0.5 (savings 8.5, 10.7, 12.4, 11.3, 6.5 times): best near 0.25. The round strut's drag is 10.3, 23.2, 41.1 N at 30, 45, 60 m/s, almost exactly ∝ V² (c_d changes from 1.009 to 1.005), and it sheds at 254, 382, 509 Hz.

Raw data format: CSV `theta_s_deg, cp_base, cd` (1), `Re, cd_stokes, cd_sn, cd_morrison` (2), `object, d_m, rho_p, V_morrison, Re, V_sn, V_stokes, V_ode` (3), `t_over_c, V, cd_round, cd_faired, ratio, f_shed_Hz` (4), plus `drag.png`.

## 5. Validate and improve

| quantity | exact or published | code | difference | cause |
|---|---|---|---|---|
| potential-flow c_d | 0 (d'Alembert, exact) | 1.1 × 10⁻¹⁶ | round-off | — |
| separated model, 82°, −1.2 | (2): 0.8839 | 0.884 (2048 points) | < 10⁻³ | discontinuous c_p |
| Stokes limit, 10 µm glass bead in water | Stokes's law: 81.66 µm/s | Morrison 81.66 µm/s | 4 × 10⁻⁶ relative, Re = 8 × 10⁻⁴ | the correlation's other terms |
| terminal velocity, steel ball | root of (3): 1.513 m/s | time integration 1.513 m/s | < 10⁻⁴ | ODE tolerance |
| strut c_d | White fit 1.005; many published curves near 1.2 | — | up to 20 % | correlation scatter, aspect ratio, turbulence |

The code is verified against exact results; its physics is only as good as the correlations, and their spread (3–10 % between sphere fits, up to 20 % for the cylinder) is the uncertainty of the engineering answer. For the strut: **41–49 N and 2.5–2.9 kW per metre bare, 3.3 N and 0.2 kW faired**: the saving, roughly 38–46 N and 2.3–2.7 kW per metre, is far larger than any of the uncertainties, so the decision to fair is robust. The shedding frequency, about 500 Hz at cruise and 250 Hz at 30 m/s, must be checked against the strut's bending frequencies: a bare strut passes through every frequency between zero and 509 Hz during the take-off run.

What to improve: an end-effect correction for a strut of finite length, a measured c_d for the actual tube (surface, joints), and a fairing correlation checked against a wind-tunnel test. What not to improve: the numerics.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2) and Stokes's terminal velocity; `drag.png`; the four experiment tables; and a half-page recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on fairing the strut, with the drag saving as a range and the shedding check.

### Adversarial self-checks

- **What if the speed doubles?** The round strut's drag rises fourfold to about 165 N (c_d barely changes until the drag crisis near Re = 2–5 × 10⁵, which a 25 mm tube reaches only above about 150 m/s at this altitude); the fairing's rises slightly less than fourfold, because C_f falls with Re.
- **What if the velocity goes to zero?** Below Re ≈ 47 there is no shedding (St is undefined); below Re ≈ 1 the drag becomes Stokes drag, linear in V, and viscosity, not the wake, sets it.
- **Which assumption fails first?** The smooth surface and the quiet stream: roughness or free-stream turbulence trigger the drag crisis at lower Re (golf-ball dimples are designed for that), and the smooth-sphere correlations then overpredict drag by a factor of up to four.
- **What could falsify the model?** A measured c_d for the tube outside 0.9–1.3 at Re = 8 × 10⁴, or a measured shedding frequency not near 0.2V/d.
- **Does the result survive a change in method?** Root finding and time integration agree to 10⁻⁴; two independent correlations agree within 10 %; the decision (fair the strut) survives a factor of two on any of them.
- **Which decision changes if the uncertainty doubles?** None for the strut: even 40 % uncertainty leaves a saving of at least 25 N per metre. A terminal-velocity estimate used to size a recovery system would need its own margin: a ±10 % correlation spread is ±10 % on the landing speed and ±20 % on the impact energy.

## Check yourself

1. What is the difference between pressure drag and friction drag? Give a shape dominated by each.
2. Why does attached potential flow round a cylinder produce no drag, and what does a real flow do differently?
3. In the separated-wake model, which term carries most of the drag?
4. Write Stokes's drag law and the terminal velocity it gives. When is it valid?
5. What is the drag crisis, at roughly what Reynolds number does it happen for a smooth sphere, and why does a turbulent boundary layer reduce drag?
6. A 6 mm antenna wire on an aeroplane at 50 m/s sheds vortices at what frequency (St ≈ 0.2)?
7. Why is a correlation not a measurement, and what two things must you state when you use one?
8. Why does a fairing with t/c = 0.25 have less drag than both a much thinner and a much thicker one of the same thickness?
