# P05.04  Frames, Machines & Friction

> Take the machine apart, one free body per part; then let friction take its cut, in the pins and in the cables.

| | |
|---|---|
| **Learning outcomes** | 1. Analyse a frame or machine by dismembering it into two-force and multi-force members with equal and opposite pin forces, and solve a control linkage (stick, pushrod, bellcrank, pushrod, control horn) member by member. 2. Compute the linkage's kinematics from loop closure (two circle intersections) and its gearing by differencing the kinematics, and verify the frictionless stick force independently by virtual work. 3. State Coulomb's law of dry friction, the friction angle and the friction circle of a pin, and solve the linkage with pin friction by fixed-point iteration, finding the driving force, the back-driven force and the mechanical efficiency. 4. Derive the capstan law T₂/T₁ = e^{μθ} and the tension ratio of a cable over a pulley turning on a pin, verify the first by numerical integration and the second by its limits, and compare the friction of a pushrod and a cable control run. 5. Judge where friction, toggle positions and gearing nonlinearity matter in a flight-control system. |
| **Prerequisites** | `P05.01` (two-force members), `P05.02` (rigid-body equilibrium), `P03.10` (nonlinear systems and fixed-point iteration), `P02.01` (first-order ODEs, for the capstan check), `P01.07` (vectors). Skill check: you can find the intersection points of two circles and solve a scalar equation with `scipy.optimize.brentq`. |
| **Engineering question** | An elevator is driven from the control stick through a pushrod, a bellcrank, a second pushrod and a control horn. For an aerodynamic hinge moment that grows with deflection, what stick force must the pilot apply over the travel, how does the gearing change with stick angle, and how much do friction in the pins and, for comparison, in a cable run add to the force and to the dead band between pushing and releasing? |
| **Assumptions and limits** | Rigid planar links; pushrods pinned at both ends (rod-end friction neglected, so they stay two-force members); quasi-static motion (no inertia); Coulomb friction with one coefficient for sticking and sliding; a linear, illustrative hinge-moment model H = −kδ; small pin clearances. Cable stretch and pretension changes are ignored. |
| **Analytical work** | The method of frames and machines; loop-closure kinematics; moment equilibrium of the stick, bellcrank and elevator about their pivots; virtual work for a frictionless machine; Coulomb friction, the friction angle and the friction circle; mechanical efficiency; the capstan law; the pulley-on-a-pin tension ratio and its limits; the toggle (over-centre) mechanism. |
| **Computational deliverable** | `starter/main.py`: `circle_intersections`, `solve_kinematics`, `gearing`, `hinge_moment`, `friction_radius`, `linkage_forces`, `stick_force_virtual_work`, `friction_power_balance`, `mechanical_advantage`, `efficiency`, `capstan_ratio`, `capstan_numeric`, `pulley_ratio`, `cable_run_ratio` (`cross2`, `unit`, the `Linkage` dataclass with its geometry methods, and `main` are given). Tests prove: circle intersections lie on both circles; pushrods keep their lengths over the travel and neutral maps to neutral; the small-angle gearing equals the product of arm ratios; every free body is in force and moment equilibrium with and without friction; the chain equals virtual work; the aerodynamic resultant reproduces the hinge moment; zero hinge moment needs no force; the friction band brackets the ideal force and shrinks with μ; the power balance closes; the capstan law matches an ODE; the pulley ratio has the right limits; a cable run multiplies its elements. |
| **Practical activity** | **computational experiment**: a sweep of stick angle over the full travel (elevator angle, gearing, stick force with and without friction), a sweep of the pin friction coefficient, and a comparison of a pushrod run with cable runs on plain-bushing and ball-bearing pulleys. |
| **Equipment** | None. A computer with the environment from `P00.01`. A bicycle brake cable or a rope over a broom handle is a welcome way to feel the capstan effect, but nothing is measured with it. |
| **Safety** | None beyond normal computer use. If you try the rope-and-broom-handle demonstration, use light loads and keep fingers clear of the rope where it wraps. |
| **Evidence** | Completed `main.py`, passing tests, `control_linkage.png`, free-body diagrams of the stick, bellcrank and elevator, the hand calculation of the small-angle gearing and stick force, the derivation of the capstan law, the three experiment tables, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-05/grading-rubric.md`; competencies 1–4. Pass: each member's free-body diagram is drawn with equal and opposite pin forces; the frictionless stick force is computed both member by member and by virtual work; all public tests pass; the friction results are interpreted as forces the pilot feels. |
| **References** | `@libretexts-baker-haynes-engineering-statics` Ch. 6 "Equilibrium of Structures": 6.2 Interactions between members, 6.6 Frames and Machines (required-free: dismembering a frame); Ch. 9 "Friction" on engineeringstatics.org (Chapter_09-dry-friction): 9.1 Dry Friction (9.1.1 Coulomb Friction, 9.1.2 Friction Angle and Friction Resultant), 9.5 Flexible Belts, 9.6 Journal Bearings (required-free: Coulomb friction, the capstan law and pin friction). The control-linkage geometry, its loads and the pulley-on-a-pin ratio are original to this lesson. |
| **Extensions** | ★ Add friction at the four rod-end pins: each pushrod then carries a small moment and is no longer a two-force member; re-solve and compare the efficiency. ★ Replace the linear hinge moment by one with a tab or an aerodynamic balance and design the gearing so the stick force per degree is constant over the travel. ★ Model a closed cable loop with pretension T₀: show that the friction in a pulley depends on the sum of the two cable tensions, so pretension costs friction. |
| **Time** | unmeasured |

## 1. The engineering problem

A pilot moves the elevator through a machine. In a typical light aircraft with pushrods, the control stick is a lever pivoted on the cockpit floor; a pushrod runs from its lower end to a **bellcrank**, a pivoted plate with two arms that turns the motion through an angle; a second pushrod runs from the bellcrank to the **control horn**, an arm fixed to the elevator near its hinge. Everything the pilot feels is the aerodynamic hinge moment on the elevator, reflected back through this chain, plus whatever friction the chain adds.

Two questions drive the design. **Gearing**: how many degrees of elevator per degree of stick, and is that ratio constant over the travel? A ratio that changes makes the aircraft feel different at large deflections. **Friction**: certification rules and pilots both dislike a control that needs a "breakout" force before it moves, or that does not return to centre. A machine with friction needs more force to drive than the air needs to hold, and less force to hold than the air needs to push it back: the difference is a dead band the pilot feels.

This lesson analyses the linkage as a **machine**, one free-body diagram per part, solves its kinematics by loop closure, checks the frictionless force by **virtual work**, adds **pin friction** with the friction circle, and compares a pushrod run with a **cable run** whose friction follows the capstan law. With an illustrative hinge moment of 36.4 N·m at 25° of elevator, the frictionless stick force at full travel is **84.2 N**; pin friction with μ = 0.2 on 6.35 mm bolts raises it to **87.3 N** when pushing and lowers it to **81.2 N** when the air pushes back, an efficiency of **96.4 %**. A cable run through four plain-bushing pulleys and three fairleads widens that band from 6.1 N to **18.8 N**.

Success criteria, written first: every member in equilibrium; the member-by-member stick force equal to the virtual-work force over the whole travel; the friction solution closing an energy balance; the capstan law equal to an ODE solution.

## 2. Understand and derive

### 2.1 Frames and machines

A **truss** has only two-force members. A **frame** or **machine** has at least one **multi-force member**: a part loaded at three or more points, or loaded between its pins, which carries bending as well as axial force (Baker & Haynes 6.2, 6.6). A frame is meant to stay put; a machine transmits forces and motion. The method is the same for both:

1. Draw the free-body diagram of the whole, if it helps to find external reactions.
2. **Dismember**: draw each part separately. Where two parts share a pin, the pin force on one is equal and opposite to the pin force on the other (Newton's third law).
3. Identify the **two-force members** first: each has a single unknown, along the line of its pins.
4. Write equilibrium for each multi-force member, choosing moment centres that eliminate unknowns.

The elevator linkage, in side view (x aft, y up), has three multi-force members, each on a fixed pivot: the **stick** (pivot S₀, grip 0.55 m above it, pushrod pin 0.08 m below it), the **bellcrank** (pivot C₀ at x = 1.20 m, arm 1 of 0.08 m pointing down, arm 2 of 0.10 m pointing aft), and the **elevator** (hinge E₀ at (1.38, 0.60) m, horn of 0.08 m pointing forward). Two pushrods, 1.20 m and 0.60 m long, are two-force members. At neutral (stick vertical, φ = 0) the first pushrod is horizontal and the second vertical, each perpendicular to the arms it joins.

```
   grip
    |
    |                                        E0 o==== elevator (chord aft)
    |                                    H o---'
    o S0 (stick pivot)                     |  horn 0.08 m
    |                                      |  pushrod 2 (0.60 m)
    o Ps ------------ pushrod 1 ------ o B1|
                                      |    |
                               C0 o---+----o B2   bellcrank: arm 1 down, arm 2 aft
```

### 2.2 Kinematics by loop closure

A pushrod keeps its length. So for a stick angle φ, the bellcrank pin B₁ lies both on the circle of radius b₁ about C₀ and on the circle of radius L₁ about the stick's pushrod pin P_s(φ): it is one of the two **intersections of two circles**. Of the two, the one that is a continuous deformation of the neutral assembly is chosen (the other is the linkage assembled "crossed"). The bellcrank angle ψ follows, then the bellcrank's second pin B₂(ψ), and the elevator angle δ from a second circle intersection: the horn pin on the circle of radius h about E₀ and the circle of radius L₂ about B₂.

The **gearing** dδ/dφ is the velocity ratio of the machine. For small angles, with each pushrod perpendicular to its arms, a stick rotation dφ moves the pushrod by a_s dφ, which turns the bellcrank by a_s dφ/b₁, which moves the second pushrod by (b₂/b₁) a_s dφ, which turns the elevator by (b₂/b₁)(a_s/h) dφ:

**dδ/dφ ≈ −(a_s/b₁)(b₂/h) = −(0.08/0.08)(0.10/0.08) = −1.25**   (1)

(negative: pushing the stick forward, φ > 0, moves the elevator trailing edge down, δ < 0). Here the stick and the bellcrank form a parallelogram (equal 0.08 m arms, a pushrod as long as the pivot spacing), so ψ = φ exactly at every angle; the nonlinearity comes from the bellcrank–horn stage, whose arms are unequal. Over the travel the code finds dδ/dφ = −1.2500 at neutral, −1.2737 at φ = +20° and −1.3261 at −20°, and δ = −25.19° at φ = +20° and +25.44° at −20°: the gearing grows by 2–6 % at the ends, slightly asymmetrically.

### 2.3 Statics member by member

Take each pushrod force P₁, P₂ as a tension. For the **elevator**, moments about its hinge: the aerodynamic hinge moment H(δ) plus the moment of the pushrod-2 force on the horn must vanish. For the **bellcrank**, moments about C₀: the two pushrod forces. For the **stick**, moments about S₀: the pushrod-1 force and the pilot's grip force F (taken perpendicular to the stick, positive in the direction of increasing φ). In each equation the pivot reaction drops out:

H + P₂ c_e = 0,  P₁ c_b1 + P₂ c_b2 = 0,  F c_g + P₁ c_s = 0,   (2)

where each c is the moment about the pivot of a unit force along the corresponding pushrod (or, for c_g, a unit grip force). Solve from the elevator back to the stick. The pivot reactions then follow from each member's force balance, for example R_stick = −(F t_g + P₁ u₁).

**The hinge moment.** The aerodynamic moment on the elevator about its hinge is H = C_h q S_e c_e, with q the dynamic pressure, S_e and c_e the elevator's area and chord, and C_h a hinge-moment coefficient that grows with deflection. The lesson uses an illustrative linear model, H = −kδ with k = −(∂C_h/∂δ) q S_e c_e: for ∂C_h/∂δ = −0.15 per radian, q = ½ × 1.225 × 50² = 1 531 Pa (50 m/s at sea level), S_e = 1.2 m² and c_e = 0.30 m, k = **82.7 N·m/rad**. At δ = −25.19° the hinge moment is **36.4 N·m**, restoring.

At φ = 20° the chain gives P₂ = −497.0 N, P₁ = −616.0 N (both pushrods in compression when pushing the stick forward) and **F = 84.19 N**. The small-angle estimate F ≈ |H| × 1.25/0.55 with the hinge moment at 25° gives 82.0 N, 2.6 % low, because it ignores the growth of the gearing.

### 2.4 An independent check: virtual work

For a frictionless machine in equilibrium, the work done by all external forces in any small motion compatible with the constraints is zero. The only such motion here has one degree of freedom: the stick turns by δφ, the grip moves L_g δφ along its path, and the elevator turns by δδ = (dδ/dφ) δφ. So F L_g δφ + H δδ = 0, or

**F = −H (dδ/dφ) / L_g.**   (3)

Equation (3) needs no pushrod force, no lever arm and no reaction: only the hinge moment and the gearing, which `gearing` measures by differencing the kinematics. It gives 84.1920 N at φ = 20°, equal to the chain value to 10⁻⁷. The **mechanical advantage**, hinge moment resisted per unit of moment applied at the stick, is 1/|dδ/dφ|: 0.800 at neutral, 0.785 at 20°. Virtual work is often the fastest way to the force of a machine with one degree of freedom; the member-by-member analysis is still needed for the pushrod and pivot loads that size the parts.

### 2.5 Dry friction

**Coulomb's law** (Baker & Haynes 9.1.1): two surfaces pressed together by a normal force N resist sliding with a tangential force up to μ_s N, and when sliding, with μ_k N, opposite to the relative motion. This lesson uses one coefficient μ for both. The resultant of the normal and friction forces leans away from the normal by at most the **friction angle** φ_f = arctan μ (9.1.2). Example: an aircraft parked with all wheels braked on a slope of angle α stays put while tan α ≤ μ; tyre-on-concrete friction coefficients of 0.5 and more hold it on any ramp it will meet.

**The friction circle of a pin.** A pin of radius r in a bearing carries a reaction R. When the joint turns, the contact point creeps up the bearing until the reaction leans by the friction angle; its line then misses the pin centre by **r_f = r sin φ_f = r μ/√(1 + μ²)**, the radius of the friction circle (9.6). The reaction therefore has a moment **r_f |R|** about the pin centre, always opposing the rotation of the joint. For a 6.35 mm (1/4 in) bolt, r = 3.175 mm, and μ = 0.2 (an illustrative value for an unlubricated bushing), r_f = **0.6227 mm**; the small-μ approximation μr gives 0.635 mm.

### 2.6 The linkage with pin friction

Add the friction moments M_f,s, M_f,b, M_f,e at the three pivots to (2). Each opposes the rotation of its member: if the stick moves toward +φ (direction +1), the stick turns positively, the bellcrank as dψ/dφ and the elevator as dδ/dφ, and each friction moment has the opposite sign. Its size, r_f |R|, depends on the pivot reaction, which depends on the forces we are solving for. The system is **nonlinear**, but only weakly (r_f is a fraction of a millimetre against arms of centimetres), so **fixed-point iteration** (`P03.10`) converges fast: solve (2) without friction, compute the reactions and friction moments, re-solve with them, repeat until the friction moments stop changing (10⁻¹² relative). It takes 7 iterations here.

Results at φ = 20°, μ = 0.2: **F = 87.31 N driving** (the pilot pushes the elevator further against the air) and **81.19 N back-driven** (the pilot holds the stick while the air pushes it back): a **friction band of 6.12 N**. The pivot reactions are 715.5 N (stick), 800.8 N (bellcrank) and 968.6 N (elevator hinge, which also carries the aerodynamic load). The **efficiency** in driving, F_ideal/F_actual, is **0.9643**.

A second check closes the energy balance: power in at the grip equals power into the air plus power lost in the pins, F L_g φ̇ + H δ̇ = Σ r_f |R_i| |θ̇_i|. With the converged reactions the balance closes to 3 × 10⁻⁸ N·m (per unit φ̇).

### 2.7 Cables: the capstan law

Many aircraft use cables instead of pushrods: lighter over long runs, easy to route around corners with pulleys and **fairleads** (guides that the cable slides through). A cable sliding over a fixed curved surface loses tension to friction. Take a short element subtending dθ: the tensions T and T + dT pull at its ends, the surface pushes it with a normal force T dθ, and friction μT dθ opposes the sliding. Tangential equilibrium gives dT = μ T dθ, so

**T₂/T₁ = e^{μθ}**   (4)

for a wrap angle θ (9.5). The exponential is ruthless: μ = 0.2 over 180° gives a ratio of 1.874; over two full turns, 12.3. That is why a sailor can hold a boat with a few turns on a post, and why a control cable must only graze a fairlead (3° at μ = 0.2 costs 1.05 %) and must turn corners on **pulleys**, which roll instead of slide.

### 2.8 A pulley on a pin

A pulley of radius R turns on a pin of radius r_b. The cable does not slide on the pulley; the friction is in the pin. The cable tensions T_t (tight side) and T_s (slack side) both pull on the pulley, their resultant has magnitude √(T_t² + T_s² − 2T_tT_s cos θ) for a wrap angle θ, and moment equilibrium about the centre with the friction-circle moment gives

(T_t − T_s) R = r_f √(T_t² + T_s² − 2 T_t T_s cos θ).   (5)

With k = T_t/T_s and ρ = r_f/R, squaring gives a quadratic, (1 − ρ²)k² − 2(1 − ρ² cos θ)k + (1 − ρ²) = 0, whose root above 1 is the tension ratio. Two limits check it: at θ = 180° the two tensions are parallel and k = **(1 + ρ)/(1 − ρ)** exactly; for small friction, the resultant is about 2T sin(θ/2) and k ≈ 1 + 2ρ sin(θ/2). A 38.1 mm pulley on a 5 mm-radius pin: with a plain bushing (μ_b = 0.15), k = 1.0279 at 90°; with a ball bearing (an equivalent μ_b of 0.005), 1.00093. If a pulley **seizes**, the cable slides over it and (4) applies instead: 180° at μ = 0.3 gives 2.566.

### 2.9 Where the model stops being valid

- **Stick-slip.** Real friction is higher at rest than in motion; a control can jump when it breaks free. The single-μ model gives the average band, not the breakout spike.
- **Rod ends and bearings.** Friction at the pushrod ends adds moments to the pushrods, which then bend slightly; sealed bearings add a torque that does not scale with load.
- **Cable systems.** A real cable loop is pretensioned; pulley friction then depends on the sum of both cable tensions, and temperature changes the pretension (aluminium airframe, steel cable).
- **Hinge moment.** The linear model is illustrative; real hinge moments depend on angle of attack, tabs, balance and Mach number.
- **Toggle positions.** If a pushrod becomes collinear with an arm, the gearing goes to zero and the force ratio to infinity; designs keep well away from these positions except where an over-centre lock is wanted (Problem 3).

## 3. Build the computational model

Open `starter/main.py`. `cross2`, `unit`, the `Linkage` dataclass (its geometry, the neutral pushrod lengths set in `__post_init__`, and the methods `stick_points`, `bellcrank_points`, `horn_point`) and `main` are given.

### Step 1: kinematics

`circle_intersections(c0, r0, c1, r1)` returns the two points or raises `ValueError` when the circles do not meet: with d = |c₁ − c₀|, a = (r₀² − r₁² + d²)/(2d), h = √(r₀² − a²), the points are c₀ + a ê ± h n̂. `solve_kinematics(lk, φ)` applies it twice, picking the solution nearer the neutral position, and returns (ψ, δ) wrapped to (−π, π]. `gearing` differences it with a step of 10⁻⁶ rad.

### Step 2: statics with friction

```python
for it in range(1, max_iter + 1):
    P2 = -(Hm + Mf[2]) / c_e
    P1 = -(P2 * c_b2 + Mf[1]) / c_b1
    F = -(P1 * c_s + Mf[0]) / c_g
    R_hinge = -(F_aero - P2 * u2)
    R_bell = -(-P1 * u1 + P2 * u2)
    R_stick = -(F * tg + P1 * u1)
    new = -sgn * rf * np.array([norm(R_stick), norm(R_bell), norm(R_hinge)])
    if np.max(np.abs(new - Mf)) <= tol * max(1.0, abs(Hm)):
        break
    Mf = new
```

Before the loop, `linkage_forces` computes the configuration, the unit vectors along the pushrods, the moment coefficients c of (2), the aerodynamic resultant (a force H/e_cp normal to the chord at e_cp = 0.075 m aft of the hinge, so that it reproduces H) and the rotation signs. It returns the forces, reactions, friction moments and the iteration count. `stick_force_virtual_work` is (3); `friction_power_balance` the energy balance of §2.6; `mechanical_advantage` and `efficiency` follow from them.

### Step 3: cables

`capstan_ratio` is (4); `capstan_numeric` integrates dT/dθ = μT with `solve_ivp` (rtol 10⁻¹¹); `pulley_ratio` is the root of the quadratic of §2.8, raising `ValueError` when ρ ≥ 1 (a friction circle larger than the pulley: it cannot turn); `cable_run_ratio` multiplies the ratios of a list of fairleads and pulleys.

Running `main()` prints the hinge stiffness, the gearing, the force table over the travel by both routes, the friction results, the cable comparison and the sweep, and writes `control_linkage.png` (the linkage drawn at φ = 15°; δ and the gearing against φ; the stick force frictionless, driving and back-driven).

### Verify

Run `python -m pytest curriculum/phase-05/lesson-04-frames-machines-friction`. The tests check:

1. **Circle intersections** lie on both circles; disjoint circles are refused.
2. **Kinematics**: neutral maps to (0, 0); both pushrods keep their lengths and the horn its radius over the travel.
3. **Gearing**: dψ/dφ = a_s/b₁ and dδ/dφ = −(a_s/b₁)(b₂/h) at neutral; halving b₂ halves the gearing and doubles the mechanical advantage.
4. **Equilibrium of every free body** (forces and moments of stick, bellcrank and elevator, with and without friction, both directions).
5. **Chain against virtual work** at nine angles over ±20°.
6. **Aerodynamic resultant** reproduces the hinge moment.
7. **Limits**: zero hinge moment needs no force even with friction; the force is linear in k.
8. **Friction band**: back-driven < ideal < driving; the band shrinks with μ and halves when μ halves (small μ); efficiency 1 at μ → 0 and falling with μ.
9. **Power balance** closes with friction in both directions.
10. **Friction radius**: zero, μr and r limits.
11. **Capstan** against `solve_ivp`; e^{μ·2π} = (e^{μπ})².
12. **Pulley**: (1 + ρ)/(1 − ρ) at 180°; 1 + 2ρ sin(θ/2) for small friction; no friction gives 1; the root satisfies (5); an oversized friction circle is refused.
13. **Cable run** is the product of its elements.
14. `main` writes the figure.

### Compare

The member-by-member chain and virtual work agree to 10⁻⁷ over the whole travel (the residual difference is the central-difference error of the gearing). The friction solution passes two independent tests: each free body is in equilibrium with the friction moments, and the energy balance closes. The capstan law matches the ODE to 10⁻⁸. A sign error in a friction moment (friction helping the motion) is the classic mistake; it makes the driving force smaller than the ideal one, which the band test catches at once.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results are computations with the model; the hinge-moment model and friction coefficients are illustrative, not measured. Write each prediction down before running it.

1. **Travel sweep.** For φ from −20° to +20°, predict from (1) the elevator angle at the ends and whether the gearing stays at −1.25; then compute δ, dδ/dφ and the stick force with and without friction (μ = 0.2).
2. **Friction coefficient.** At φ = 20°, sweep μ = 0.1, 0.2, 0.4, 0.8. Predict how the band and the efficiency scale. Then try 12.7 mm (1/2 in) pins at μ = 0.2.
3. **Pushrods against cables.** Compare the pushrod linkage's band with a cable run of four pulleys (wraps 90°, 45°, 30°, 90°; R = 38.1 mm; pin radius 5 mm) and three fairleads (3° each, μ = 0.2), with plain-bushing pulleys (μ_b = 0.15) and with ball-bearing pulleys (equivalent μ_b = 0.005), assuming for simplicity that the cable carries the same ideal force path so the driving force is F_ideal × k and the back-driven force F_ideal/k. Then seize the 90° pulley (μ = 0.3 sliding).

Raw data format: CSV tables `phi_deg, delta_deg, gearing, F_ideal_N, F_drive_N, F_back_N` (1); `mu, r_pin_mm, F_drive_N, F_back_N, efficiency` (2); `system, ratio, F_drive_N, F_back_N, band_N` (3); plus `control_linkage.png`.

What you should find (model results): (1) δ runs from +25.44° at φ = −20° to −25.19° at +20°; the gearing from −1.3261 through −1.2500 at neutral to −1.2737; the frictionless force rises from 41.37 N at 10° to 84.19 N at 20° (and −41.61, −88.52 N on the pull side), slightly more than linearly because the gearing grows. (2) Band 3.10, 6.12, 11.59 and 19.53 N; efficiency 0.9817, 0.9643, 0.9336 and 0.8909. The band is nearly proportional to μ at small μ and less than proportional at large μ, because the friction radius saturates at r. Doubling the pin diameter to 12.7 mm doubles the friction moments: efficiency 0.9300. (3) Pushrods: band 6.12 N. Cable run with plain pulleys: ratio 1.1179, driving 94.12 N, back-driven 75.31 N, **band 18.81 N**, three times the pushrods'; with ball-bearing pulleys: ratio 1.0347, band 5.75 N, comparable to the pushrods. A seized 90° pulley alone makes the run's ratio 1.656: 139.4 N to drive and 50.8 N to hold, a control that feels broken. Pulley bearings matter more than anything else in a cable system.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| gearing at neutral | −1.25 (1) | −1.2500 | 10⁻⁹ | finite difference |
| stick force at φ = 20° (N) | 82.0 (small-angle) | 84.19 (chain) = 84.19 (virtual work) | 2.6 %; 10⁻⁷ | gearing growth; difference step |
| friction radius (mm) | 0.635 (μr) | 0.6227 (exact) | 2 % | small-μ approximation |
| capstan, μ = 0.2, 180° | e^{0.2π} = 1.874456 | 1.874456 (ODE) | 10⁻⁸ | solver tolerance |
| pulley at 180° | (1 + ρ)/(1 − ρ) | equal | 10⁻¹² | — |
| power balance with friction (N·m) | 0 | 3 × 10⁻⁸ | iteration tolerance | — |

The model agrees with every closed form and with its own energy balance. Its weaknesses are the inputs: the hinge-moment model and the friction coefficients are illustrative. The most useful improvement would be measured values: a hinge-moment coefficient from wind-tunnel or flight data (Phase 10) and friction from a simple test, hanging weights on the stick of a real system and recording when it moves each way, which gives the band directly.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the free-body diagrams of the stick, bellcrank and elevator with equal and opposite pin forces; the hand derivations of (1), (3) and (4); `control_linkage.png`; the three experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question: the stick force over the travel, the gearing nonlinearity, and the friction band of the pushrod and cable alternatives.

### Adversarial self-checks

- **What if the airspeed doubles?** The hinge moment, the stick force, the pushrod loads, the reactions and hence the friction moments all scale with q, four times; the efficiency does not change (friction is proportional to load), but the band in newtons quadruples.
- **What if the hinge moment goes to zero?** No force is needed and, in this model, no friction appears. A real system still has a small no-load friction (seals, bearing preload) that the model omits.
- **What if a pushrod lines up with an arm?** The gearing goes to zero, the force ratio to infinity, and a small error in geometry makes a large error in force: a toggle. Over-centre locks use this on purpose.
- **Which assumption fails first in the real world?** The single friction coefficient: breakout (static) friction is higher than running friction, so the pilot feels a step when the control starts to move.
- **What measurement could falsify the model?** Hang known weights on the grip with a known hinge moment applied (a spring on the elevator) and record the forces at which the stick starts to move each way; the mean should match the ideal force and the difference the predicted band.
- **Does the result survive a change of method?** Yes: chain and virtual work agree to 10⁻⁷, and the friction solution passes both the member equilibrium and the energy balance.

## Check yourself

1. What makes a member a multi-force member, and why does it matter for the analysis?
2. Why are the pin forces on two connected parts equal and opposite, and how many unknowns does a pin between two parts introduce?
3. Write the small-angle gearing of a stick–pushrod–bellcrank–pushrod–horn linkage in terms of its four arm lengths.
4. State the principle of virtual work for a frictionless machine with one degree of freedom, and use it to write the stick force.
5. What is the friction circle of a pin, and what is its radius for μ = 0.2 and r = 5 mm?
6. Why is the force needed to push a control against the air larger than the force needed to hold it, and what is the difference called here?
7. Derive the capstan law. What tension ratio does μ = 0.2 give over a full turn?
8. Why do control cables go round corners on pulleys rather than through fairleads, and what happens if a pulley seizes?
