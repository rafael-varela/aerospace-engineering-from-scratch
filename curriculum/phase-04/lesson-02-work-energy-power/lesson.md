# P04.02  Work, Energy & Power

> Altitude and speed are two bank accounts holding the same currency. The engine pays in at the rate (T − D)V; the pilot moves money between the accounts for free.

| | |
|---|---|
| **Learning outcomes** | 1. Compute the work of a force along a path as a line integral, show numerically and by argument that the work of gravity is path-independent and the work of drag is not, and state the work–energy theorem with non-conservative forces. 2. Derive the energy height h_E = h + V²/2g and the specific excess power P_s = (T − D)V/W = dh_E/dt of an aircraft from the work–energy theorem, and use them for rate of climb, zoom climb, ceilings, turn performance and the acceleration factor. 3. Derive the speeds of minimum drag and minimum power, (L/D)max and the glide angle, and Lanchester's phugoid period √2πV/g together with its correction for the density gradient. 4. Integrate the longitudinal point-mass equations with `scipy.integrate.solve_ivp`, verify that the energy height is conserved to 10⁻⁸ m without thrust and drag and changes at exactly P_s with them, and reproduce the phugoid and glide closed forms. 5. Find the minimum-time climb schedule by the energy method (Rutowski) and judge, with numbers, how much it gains over simpler schedules for a subsonic aircraft. |
| **Prerequisites** | `P04.01` (Newton's second law, drag, free-body diagrams), `P01.07` (dot product: F·dr), `P02.04` (first-order systems), `P03.05` (what `solve_ivp` does), `P03.10` (the regional-jet trim: drag polar, ISA density). Skill check: you can write the drag of an aircraft as qS(C_D0 + kC_L²) and find the minimum of a one-variable function. |
| **Engineering question** | A 5 t jet trainer flies at sea level at 150 m/s. How fast can it climb, how high can it zoom, what are its absolute and service ceilings, and what speed schedule takes it to an energy height of 10 km in the least time, and how much does that schedule gain over flying a constant speed? |
| **Assumptions and limits** | Point mass with constant mass (fuel burn ignored), flat non-rotating Earth, constant g; parabolic drag polar C_D = C_D0 + kC_L² with constant coefficients (no compressibility: the speed is capped at Mach 0.75 and above it the model is invalid); thrust T₀(ρ/ρ₀) along the flight path and independent of speed (a crude turbojet model, the first assumption to check against an engine deck); ISA atmosphere 0–20 km. The energy method assumes that speed and altitude can be exchanged instantly along a line of constant h_E, which a real zoom or dive takes tens of seconds to do. |
| **Analytical work** | Work as a line integral; work of gravity and of drag; the work–energy theorem with non-conservative forces; h_E and P_s; V_md, V_mp, (L/D)max = 1/(2√(kC_D0)); the glide ratio; the steady-climb rate and the acceleration factor; Lanchester's phugoid from energy conservation at constant C_L and its density-gradient correction ω² = 2g²/V² + gβ; the zoom altitude; the minimum-time condition ∂P_s/∂V = 0 along h_E = const. |
| **Computational deliverable** | `starter/main.py`: `work_along_path`, `kinetic_energy`, `energy_height`, `zoom_altitude`, `drag`, `thrust`, `specific_excess_power`, `power_required_available`, `min_drag_speed`, `min_power_speed`, `max_lift_to_drag`, `best_climb`, `ceiling`, `point_mass_rhs`, `simulate_flight`, `phugoid_period_lanchester`, `density_scale_gradient`, `phugoid_period_with_gradient`, `oscillation_period`, `rutowski_path`, `time_to_energy`, `constant_speed_climb_ps`, `climb_rate_with_acceleration` (`isa`, `speed_of_sound`, `Aircraft`, `stall_speed`, `climb_time_accelerate_then_constant` and `main` are given). Tests prove: the ISA anchors and continuity; gravity's work is path-independent and zero round a loop while drag's is −D × arc length; the work–energy theorem on a spring; the drag and power minima at the closed-form speeds and D_min = W/(L/D)max; induced drag grows as n²; P_s as defined and vectorised; the best-climb search finds a maximum and the ceilings have P_s = 0 and 0.508 m/s and move the right way with mass and thrust; h_E constant to 10⁻⁶ m in a thrust-free, drag-free flight while h swings 50 m; dh_E/dt = (T − D)V/W along a powered trajectory to 10⁻⁴ m/s; the phugoid period equal to Lanchester's to 0.2 % at constant density and to the gradient formula to 0.5 % in the ISA; a steady glide at (L/D)max to 10⁻⁸; the Rutowski path on its energy curves and maximising P_s; the energy-method climb no slower than three alternatives. |
| **Practical activity** | **computational experiment**: a P_s map over speed and altitude with the minimum-time path, a phugoid study (constant and ISA density, small and large amplitude), a glide study, a comparison of climb schedules, and a sensitivity study of the ceiling to mass and thrust. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `energy_height.png`, the derivations of P_s, V_md, V_mp, (L/D)max and the phugoid period with its correction, the five experiment tables, and a half-page note answering the engineering question with numbers and stating which assumption limits each answer. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: P_s and the phugoid period are derived from energy; all public tests pass; the climb-schedule comparison is reported with times and the small gain explained; the note names the assumption behind each number. |
| **References** | `@openstax-university-physics-vol1` Ch. 7 "Work and Kinetic Energy": 7.1 Work, 7.2 Kinetic Energy, 7.3 Work-Energy Theorem, 7.4 Power (required-free: work as a line integral, the theorem, P = F·v); Ch. 8 "Potential Energy and Conservation of Energy": 8.1 Potential Energy of a System, 8.2 Conservative and Non-Conservative Forces, 8.3 Conservation of Energy, 8.4 Potential Energy Diagrams and Stability (required-free: path independence, energy bookkeeping with losses, turning points). `@leishman-intro-aerospace-flight-vehicles` chapter "Airplane Equations of Motion", https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/airplane-equations-of-motion/ (supplementary-free: the point-mass equations along and normal to the flight path). |
| **Extensions** | ★ Let the mass fall with fuel burn (specific fuel consumption 2.5 × 10⁻⁵ kg/(N·s)) and find the energy-optimal climb that minimises fuel instead of time (P_s divided by fuel flow). ★ Add a compressible drag rise above Mach 0.75 and a supersonic thrust model and reproduce the classic dive-through-Mach-1 minimum-time path. ★ Fly the Rutowski schedule in `simulate_flight` with a load-factor controller and measure how far the real path deviates from the instantaneous-exchange idealisation. |
| **Time** | unmeasured |

## 1. The engineering problem

A pilot at sea level with the throttle open has two ways to spend the engine's work: climb, or accelerate. A fighter pilot sees a third: trade speed for height in a zoom, or height for speed in a dive, almost for free. Which of these is "free" and which costs fuel is the subject of energy bookkeeping, and it is the most practical piece of mechanics in aircraft performance. It replaces a tangle of forces along and across the flight path by one scalar, the **energy height** h_E = h + V²/2g, the height the aircraft would reach if it converted all its kinetic energy into altitude, and one rate, the **specific excess power** P_s = (T − D)V/W, the rate at which the engine, net of drag, adds to it. Climb rate, ceilings, acceleration, sustained turns and the minimum-time path to a given altitude and speed all follow from these two quantities.

The vehicle is a twin-engine jet trainer: 5000 kg, wing area 20 m², C_D0 = 0.022, k = 0.060, sea-level thrust 16 kN falling with density, capped at Mach 0.75. At sea level and 150 m/s its energy height is 1147 m, its kinetic energy 56.3 MJ, the drag 6.59 kN against 16.0 kN of thrust, and P_s = 28.8 m/s: it could climb at 28.8 m/s at constant speed, or accelerate at P_s g/V = 1.9 m/s² in level flight, or any mixture. An ideal zoom to 80 m/s reaches 821 m. The best climb rate at sea level is 28.9 m/s at 143 m/s; it falls to 7.4 m/s at 9 km; the absolute ceiling is 12 828 m and the service ceiling (100 ft/min, 0.508 m/s) 12 572 m. The minimum-time climb to h_E = 10 km takes 577 s along a schedule from 153 to 167 m/s; flying a constant 150 m/s takes 585 s, and a constant 120 m/s takes 714 s. That a careful optimisation gains only 1 % over a sensible constant speed, and a poor choice costs 24 %, is an engineering answer too.

Success criteria, written first: every closed form (V_md, V_mp, (L/D)max, the glide ratio, the phugoid period) reproduced by the simulation or by a numerical minimum to the stated tolerance; h_E conserved to 10⁻⁶ m in conservative flight; dh_E/dt equal to P_s along a powered trajectory; the climb comparison with times and an explanation of the size of the gain.

## 2. Understand and derive

### 2.1 Work along a path

The **work** done by a force F on a body that moves along a path C is the line integral **W = ∫_C F · dr** (OpenStax 7.1): only the component of the force along the displacement does work. Units: N·m = J. Two forces in this lesson behave very differently.

The weight F = (0, −mg) is constant, so W = −mg ∫dz = −mg(z₁ − z₀) whatever the path: gravity is **conservative** (8.2), its work depends only on the endpoints, it is zero round any closed loop, and it can be written as minus the change of a **potential energy** U = mgz (8.1). The code checks this on a straight and a wavy path to the same point 500 m higher: both give −24.5166 MJ for a 5000 kg mass, equal to −mgΔz to all printed digits.

Drag is the opposite. It always points against the motion, F_D = −D v̂, so F_D · dr = −D |dr|: the work is minus the drag times the *arc length*, always negative, and larger on a longer path. A constant 5 kN drag does −5.59 MJ along the 1118 m straight path and −10.82 MJ along the wavy one. Drag is **non-conservative** (8.2): no potential energy can be defined for it, and its work is a loss that no later motion recovers.

### 2.2 The work–energy theorem with losses

Newton's second law dotted with the velocity gives m v · dv/dt = F · v, i.e. d(½mv²)/dt = F · v: the rate of change of the **kinetic energy** K = ½mv² (7.2) equals the **power** P = F · v (7.4) of the net force. Integrated along the path, ΔK = W_net (7.3). Splitting the forces into the conservative weight and the rest,

**Δ(K + U) = W_nc**,   (1)

the mechanical energy changes only by the work of the non-conservative forces (8.3): here thrust (positive) and drag (negative). This is the whole of the bookkeeping; the rest is choosing good variables.

### 2.3 Energy height and specific excess power

For an aircraft with thrust T along the flight path and drag D against it, lift L perpendicular to the path, and weight W = mg, the lift does no work (it is perpendicular to v), so (1) in rate form is d(½mV² + mgh)/dt = (T − D)V. Divide by the weight:

**h_E = h + V²/(2g),   dh_E/dt = P_s ≡ (T − D)V/W.**   (2)

h_E is the mechanical energy per unit weight, in metres; P_s is the excess power per unit weight, in m/s. At sea level and 150 m/s: q = ½ × 1.225 × 150² = 13 781 Pa, C_L = W/(qS) = 49 033/(13 781 × 20) = 0.1779, C_D = 0.022 + 0.06 × 0.1779² = 0.02390, D = qSC_D = 6.59 kN, T = 16.0 kN, P_s = (16.0 − 6.59) × 150/49.03 = **28.8 m/s**. Three uses follow at once.

- **Steady climb** at constant V: the whole of P_s goes into h, so the rate of climb is P_s. The best climb speed maximises P_s at each altitude.
- **Level acceleration**: the whole of P_s goes into V²/2g, so dV/dt = gP_s/V = 1.88 m/s² at 150 m/s.
- **Climb while the speed changes** with altitude (dV/dh ≠ 0): h_Ė = ḣ(1 + (V/g)dV/dh), so **ḣ = P_s/(1 + (V/g)dV/dh)**. The bracket is the **acceleration factor**: accelerating 5 m/s per km at 150 m/s costs 7 % of the climb rate (18.6 instead of 20 m/s).

A **zoom climb** is the case P_s ≈ 0 over a short time: h_E is conserved and h_end = h_E − V_end²/2g. From 150 m/s at sea level to 80 m/s: 1147 − 326 = **821 m**. This is the potential-energy diagram of 8.4 in disguise: along h_E = const, the "total energy" line, the kinetic energy is the gap between the line and the altitude, and the turning point of the zoom is where the gap closes to the minimum flying speed.

### 2.4 Drag, power and the special speeds

With C_L = W/(qS) in level flight, D(V) = qSC_D0 + kW²/(qS): a zero-lift term growing as V² and an induced term falling as 1/V². The minimum is where they are equal, C_D0 = kC_L²:

**V_md = √(2W/(ρS)) (k/C_D0)^{1/4},   D_min = W/(L/D)max,   (L/D)max = 1/(2√(kC_D0)).**   (3)

The power required P_r = DV has its minimum where the zero-lift drag is a third of the induced drag, C_D0 = kC_L²/3, at **V_mp = V_md/3^{1/4}** (Problem 2). For the trainer at sea level: V_md = 81.3 m/s, V_mp = 61.8 m/s, stall 53.5 m/s (C_Lmax = 1.4), (L/D)max = **13.76**, D_min = 3.56 kN. The tests find both minima numerically on a 0.005 m/s grid and recover the closed forms to 0.02 m/s.

**Gliding.** With T = 0 and a steady glide at flight-path angle γ < 0, L = W cos γ and D = −W sin γ, so tan(−γ) = D/L: the best glide angle is arctan(1/(L/D)max) = **4.16°**, at C_L = √(C_D0/k) = 0.606, and the aircraft covers 13.76 m forward per metre of height. Energy says the same thing: the glider's P_s = −DV/W is its sink rate.

### 2.5 The phugoid: energy exchange at constant lift coefficient

Hold C_L constant, remove thrust and drag, and perturb the speed. Then h_E is exactly conserved and the aircraft trades speed for height in a slow oscillation, the **phugoid**. Linearise about level flight at V₀: energy gives V δV = −g δh; the lift L = ½ρV²SC_L changes by δL/L = 2δV/V + δρ/ρ; and the vertical equation is m ḧ = δL. With constant density, δL/W = 2δV/V₀ = −2gδh/V₀², so ḧ = −(2g²/V₀²)δh:

**ω² = 2g²/V₀²,  T = √2 πV₀/g**   (Lanchester, 1908).   (4)

At 150 m/s, T = 67.96 s. In the real atmosphere the density falls with height, δρ/ρ = −βδh with **β = −(1/ρ)dρ/dh** (1.03 × 10⁻⁴ m⁻¹ at 3 km, a scale height of 9.7 km), which adds a restoring term:

**ω² = 2g²/V₀² + gβ.**   (5)

At 150 m/s and 3 km the second term is 12 % of the first (1.01 × 10⁻³ against 8.55 × 10⁻³ s⁻²) and the period falls to **64.27 s**. The phugoid is an energy-conserving oscillation in its ideal form; drag damps it slowly and thrust changes its equilibrium, which is why a pilot hardly notices it and an autopilot must not ignore it.

### 2.6 Ceilings, turns and the minimum-time climb

**Ceilings.** The maximum P_s over speed falls with altitude because the thrust lapses with density while the minimum drag does not. The **absolute ceiling** is where max P_s = 0; the **service ceiling** where it is 0.508 m/s (100 ft/min). Both are root-finding problems in h on a function that is itself a maximum over V, solved by `brentq` around a bounded `minimize_scalar`.

**Turns.** In a level turn at bank angle φ the lift is nW with load factor n = 1/cos φ, the induced drag grows as n², and P_s at n > 1 is smaller: the aircraft bleeds energy in a hard turn. The **sustained** load factor is the n at which P_s = 0 (Problem 5).

**Minimum-time climb.** To reach a final energy height in least time, minimise ∫dh_E/P_s. The integrand depends only on where you are on each h_E level, so the minimum is reached by choosing, independently on each level, the point of maximum P_s: **∂P_s/∂V = 0 along h = h_E − V²/2g** (Rutowski, 1954). The idealisation is that the aircraft can move along an h_E line instantly (zoom or dive for free); the result is a schedule V*(h_E). For a supersonic fighter this schedule includes a dive to break through the transonic drag rise and is dramatically better than any steady climb; for a subsonic aircraft with no drag rise in its envelope it lies close to the best-climb schedule, and the gain is small.

### 2.7 Where the model fails

- **Compressibility** above about Mach 0.75: the drag polar changes and the thrust model fails. The code caps V at 0.75a(h) (221 m/s at 11 km); every number near the cap is suspect.
- **Thrust independent of speed**: a turbojet's thrust varies with speed and temperature; T₀(ρ/ρ₀) is a placeholder for an engine deck.
- **Constant mass**: a 10-minute climb burns a few per cent of the mass; the sensitivity study (part 4, item 5) shows the ceiling moves by about 600 m per 10 % of mass.
- **Instantaneous exchange**: a real zoom from 150 to 80 m/s takes tens of seconds and loses energy to drag on the way; the 821 m is an upper bound.
- **Point mass**: no pitching dynamics, so the phugoid here has no coupling with the short-period mode (acceptable: the short-period pitch mode of an aircraft like this lasts seconds, the phugoid a minute).

## 3. Build the computational model

Open `starter/main.py`. `isa`, `speed_of_sound`, `Aircraft`, `stall_speed`, `climb_time_accelerate_then_constant` and `main` are given.

### Step 1: work along a path

```python
def work_along_path(force, path):
    dr = np.diff(path, axis=0)
    return sum(0.5 * np.dot(force(path[i], dr[i]) + force(path[i + 1], dr[i]), dr[i]) for i in range(len(dr)))
```

The force receives the segment direction so that drag can be written `lambda r, dr: -D * dr / np.linalg.norm(dr)`. The trapezoid rule is exact for a force linear in position (the spring test).

### Step 2: energies, forces and power

`kinetic_energy`, `energy_height` and `zoom_altitude` are one line each. `drag(ac, rho, V, n)` builds q, C_L = nW/(qS) and qS(C_D0 + kC_L²); `thrust(ac, rho)` is T₀(ρ/ρ₀)^n; `specific_excess_power(ac, h, V, n)` broadcasts h and V together and loops over the elements (the ISA call is scalar), returning an array of their common shape. `power_required_available` returns DV and TV; `min_drag_speed`, `min_power_speed`, `max_lift_to_drag` are (3).

### Step 3: best climb and ceilings

`best_climb(ac, h)` maximises P_s over [V_stall(h), M_max a(h)] with `minimize_scalar(method="bounded")`; `ceiling(ac, ps_target)` runs `brentq` on 0–20 km, raising `ValueError` if the target is out of reach at sea level or still exceeded at 20 km.

### Step 4: the equations of motion

```python
def point_mass_rhs(t, y, ac, CL, thrust_on=True, drag_on=True, const_rho=None):
    x, h, V, gam = y
    rho = isa(h)[2] if const_rho is None else const_rho
    q = 0.5 * rho * V * V
    L, D, T = q * ac.S * CL, q * ac.S * (ac.CD0 + ac.k * CL ** 2), thrust(ac, rho)
    return [V * cos(gam), V * sin(gam), (T - D) / ac.mass - G0 * sin(gam), (L - ac.W * cos(gam)) / (ac.mass * V)]
```

(with the switches applied). These are Newton's second law along the flight path and normal to it (Leishman, "Airplane Equations of Motion"); C_L is held constant, which is a pilot holding a fixed angle of attack. `simulate_flight` integrates it with DOP853 at rtol 10⁻¹⁰. `oscillation_period` interpolates upward mean-crossings; `phugoid_period_lanchester`, `density_scale_gradient` (central difference of the ISA density) and `phugoid_period_with_gradient` are (4) and (5).

### Step 5: the energy method

`rutowski_path(ac, hE_grid)` maximises P_s(h_E − V²/2g, V) over V on each level, between the stall speed and the Mach cap at the altitude reached (both found by a short fixed-point iteration on h(V)) and with h ≥ 0. `time_to_energy` integrates 1/P_s over h_E by the trapezoid rule and refuses a non-positive P_s; `constant_speed_climb_ps` evaluates P_s along h = h_E − V²/2g at fixed V.

Running `main()` prints the sea-level numbers of part 1, the best climb table (28.91 m/s at 143.2 m/s at sea level, 20.45 at 3 km, 13.40 at 6 km, 7.39 at 9 km, 3.77 at 11 km, 1.66 at 12 km), the ceilings, the two work comparisons, the phugoid and glide studies, the climb schedules and the ceiling sensitivity, and writes `energy_height.png` (P_s contours over V and h with the energy-height lines, the stall and Mach limits and the minimum-time path; the phugoid altitude and the energy-height drift in micrometres).

### Verify

Run `python -m pytest curriculum/phase-04/lesson-02-work-energy-power`. The tests:

1. **ISA**: 288.15 K, 101 325 Pa and 1.225 kg/m³ at sea level; continuous pressure and 216.65 K at 11 km; 22 632 Pa there; a = 340.3 m/s; `ValueError` above 20 km.
2. **Work**: gravity −mgΔz on a straight and a wavy path to 10⁻¹⁰, zero round a loop; drag −D × 2π round a unit circle to 10⁻⁴ and more negative on the longer path; the spring's work equals −Δ(cx²/2).
3. **Energies**: K, h_E, the zoom formula and its trivial case; K/W equals the speed part of h_E.
4. **Special speeds**: the numerical minima of D and DV within 0.02 m/s of (3); D_min = W/(L/D)max; C_D0 = kC_L² at V_md; induced drag × n²; the thrust lapse.
5. **P_s**: the definition to 10⁻¹²; vectorisation to a (2, 3) grid; less at altitude and in a turn.
6. **Climb and ceilings**: P_s falls and V* rises with altitude; V* is a local maximum; P_s at the ceilings equals 0 and 0.508 m/s to 10⁻³; heavier lowers and more thrust raises the service ceiling; an aircraft with 1 kN cannot climb.
7. **Conservation**: h_E constant to 10⁻⁶ m over 150 s of thrust-free, drag-free flight while h swings more than 50 m.
8. **Power balance**: dh_E/dt = (T − D)V/W to 10⁻⁴ m/s along a powered 60 s trajectory; the trajectory against an independent RK45 call to 10⁻⁴ m.
9. **Phugoid**: Lanchester's period to 0.2 % at constant density; the ISA period within 0.5 % of (5) and shorter than Lanchester's; β between 0.9 and 1.2 × 10⁻⁴ m⁻¹; (5) reduces to (4) at β = 0; the period estimator on a sine and its `ValueError`.
10. **Glide**: a steady glide at constant density covers (L/D)max per metre of height to 10⁻⁸ at constant speed, loses energy monotonically, and other C_L glide worse.
11. **Energy method**: the path lies on its energy levels to 10⁻¹⁰, maximises P_s within ±5 m/s, P_s falls with h_E; the time integral on a constant P_s; `ValueError` on a negative P_s; the energy-method time no more than 1 % above any of three accelerate-then-constant-speed climbs and below the constant 150 m/s climb.
12. `main` writes the figure.

### Compare

There is no library for aircraft performance in the course's stack, so the comparisons are internal and independent: the closed forms (3)–(5) against numerical minima and against the ODE solution; the energy rate against the ODE; and two solvers on the same right-hand side. The phugoid is the instructive one. At constant density the simulated period is 67.96 s, Lanchester's to four figures; with the ISA density it is 64.33 s, 5.3 % shorter, and formula (5) predicts 64.27 s. A discrepancy that looked like an error of the textbook formula was a physical effect the formula had assumed away, and a one-line extension of the same energy argument accounts for it to 0.1 %. With a 10 % initial speed perturbation (466 m altitude swing) the constant-density period is 67.99 s: the phugoid is remarkably linear.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **The P_s map.** Contour P_s over 50–240 m/s and 0–13 km, with the stall line, the Mach-0.75 line and the h_E lines. Predict: the zero contour encloses the flight envelope and closes at the ceiling near 12.8 km; the best-climb speed (the ridge) rises from about 140 m/s at sea level to about 170 m/s at the ceiling.
2. **Phugoid.** Start 1 % fast from level trim at 150 m/s, 3 km, without thrust and drag, at constant and at ISA density; then 10 % fast. Measure the period, the altitude swing and the drift of h_E. Predict 68.0 s, 64.3 s, no change with amplitude, and an h_E drift at round-off level.
3. **Glide.** At (L/D)max from 3 km, at constant density over 60 s and in the ISA over 300 s. Predict the ratio 13.76 exactly in the first case; in the second, decide before running whether descending into denser air at constant C_L lengthens or shortens the glide.
4. **Climb schedules.** From (0 m, 150 m/s) to h_E = 10 km: the energy-method schedule, and "change speed at sea level, then constant speed" at 120, 150 and 170 m/s. Predict the ranking and the size of the energy method's advantage.
5. **Sensitivity.** Service ceiling with ±10 % mass and ±10 % thrust. Predict which matters more.

Raw data format: CSV tables `V_mps, h_m, Ps_mps` (1); `case, period_s, swing_m, hE_drift_m` (2); `case, distance_over_height, V_end_mps` (3); `schedule, time_s` (4); `mass_factor, thrust_factor, service_ceiling_m` (5); plus `energy_height.png`.

What you should find: (1) the zero contour closing at 12.8 km at 172 m/s, below the Mach cap there; (2) 67.96, 64.33 and 67.99 s, swings of 41 m and 466 m, h_E drift 9 × 10⁻⁹ m over 300 s; (3) 13.7620 at constant density; 14.33 in the ISA, *longer*, because at constant C_L the speed falls from 94.3 to 85.8 m/s as the air thickens and the released kinetic energy (78 m of energy height) buys 13.8 × 78 ≈ 1.1 km of extra distance; (4) 577 s for the energy method, 585 s at 150 and at 170 m/s, 714 s at 120 m/s; (5) 12 572 m baseline, 11 967 m and 13 240 m for mass × 1.1 and × 0.9, 13 188 m and 11 890 m for thrust × 1.1 and × 0.9: about 600–700 m per 10 % either way, thrust slightly stronger, which is what a thrust-to-weight ratio argument predicts (both enter P_s through T/W).

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| P_s at sea level, 150 m/s (m/s) | 28.8 | 28.80 | 0 | — |
| zoom 150 → 80 m/s (m) | 821 | 821 | 0 | — |
| V_md, V_mp at sea level (m/s) | 81.3, 61.8 | numerical minima within 0.02 | 0 | grid |
| (L/D)max, glide angle | 13.76, 4.16° | 13.7620 (simulated, constant ρ) | 10⁻⁸ | — |
| glide in the ISA | 13.76 | 14.33 | +4 % | kinetic energy released as ρ rises |
| phugoid, constant ρ (s) | 67.96 | 67.96 | 0 | — |
| phugoid, ISA (s) | 64.27 by (5) | 64.33 | 0.1 % | higher-order terms in β |
| h_E drift, conservative flight (m) | 0 | 9 × 10⁻⁹ over 300 s | — | solver tolerance |
| dh_E/dt − P_s along a powered flight (m/s) | 0 | 5.5 × 10⁻⁶ | — | finite-difference derivative |
| work of gravity, two paths (MJ) | −24.5166 | −24.5166, −24.5166 | 0 | — |
| min-time climb vs constant 150 m/s | small gain | 577 vs 585 s | 1 % | subsonic: the schedule is near best-climb |

Every closed form is met; the two rows that differ (the ISA glide and the ISA phugoid) differ because the textbook formula assumed constant density, and in both cases the energy argument extended by one term explains the difference. That is the pattern of this phase: a discrepancy between a simple model and a more complete computation is a measurement of the neglected physics.

What to improve, if the question were real: the thrust model first (an engine deck changes the ceiling by more than any other single input), then fuel burn, then compressibility near the cap. What not to improve: the solver tolerance, or the Rutowski schedule for a subsonic trainer, where the gain over a constant 150 m/s climb (8 s in 585) is smaller than the effect of a 2 % thrust error.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (3), V_mp, the glide ratio and (4)–(5); `energy_height.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the climb rate, the zoom, the ceilings, the minimum-time schedule and its gain, each with the assumption that limits it.

### Adversarial self-checks

- **What if the mass increases by 10 %?** The service ceiling drops by 605 m to 11 967 m; the minimum-time climb to h_E = 10.5 km takes 746 s instead of 643 s (Problem 4); the zoom altitude does not change at all (h_E does not depend on mass), which is a good test of whether you understand what h_E is.
- **What if the speed goes to zero?** P_s = (T − D)V/W → 0 because the power TV does: a jet at rest has thrust but no power, and the model also leaves its domain below the stall speed. Limiting cases of the formulas: k → 0 gives V_md → 0 and no induced drag; C_D0 → 0 gives (L/D)max → ∞.
- **Which assumption fails first in the real world?** Thrust independent of speed: a real turbojet's thrust at 150 m/s can differ from its static value by 10–20 %, larger than any other modelling error here, and part 4 shows that ±10 % of thrust moves the ceiling by 600–700 m.
- **What measurement could falsify the model?** A flight-test sawtooth climb: measured rates of climb at several speeds at one altitude give the P_s curve directly; a measured best-climb speed far from 143 m/s, or a measured peak far from 28.9 m/s, falsifies the drag polar or the thrust.
- **Does the result survive a change of solver?** DOP853 and RK45 agree to 10⁻⁴ m over 60 s; the conservation check holds to 10⁻⁸ m. A different integration grid for the climb time (81 against 161 points) changes it by less than 1 s.
- **Which decision changes if the uncertainty doubles?** If thrust were known only to ±20 %, the ceiling would be uncertain by ±1.3 km and the 1 % advantage of the energy-method schedule would be invisible inside it: fly the constant-speed climb and spend the effort on the engine data.

## Check yourself

1. Define the energy height and the specific excess power, with units, and state the relation between them.
2. Why does lift do no work on an aircraft, and why does drag always do negative work?
3. The work of gravity from A to B is the same along every path; the work of drag is not. Which property of each force decides this?
4. An aircraft at 1000 m and 100 m/s has what energy height? How high can it zoom to 50 m/s?
5. At what condition on the drag coefficients is the drag a minimum, and what is the best glide ratio of an aircraft with C_D0 = 0.02 and k = 0.05?
6. What is Lanchester's phugoid period at 100 m/s, and does the density gradient of the atmosphere make the real period longer or shorter?
7. In a 60° bank level turn, what is the load factor, and why does P_s fall?
8. Why is the minimum-time climb schedule of a subsonic trainer close to its best-climb speed, while a supersonic fighter's schedule includes a dive?
