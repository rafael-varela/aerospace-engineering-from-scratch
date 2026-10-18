# P08.09  Vorticity & Circulation

> Circulation round a loop of fluid particles cannot change unless viscosity, baroclinic density gradients or a non-conservative force act on it. That one theorem explains the starting vortex, the bound vortex and the trailing vortices of every wing.

| | |
|---|---|
| **Learning outcomes** | 1. Define vorticity and circulation, prove Stokes's theorem's link between them for a plane flow, and prove Kelvin's circulation theorem with its three conditions. 2. Derive the Lamb–Oseen vortex as an exact solution of the vorticity diffusion equation, with its core growth r_c² = r_c0² + 4νt, and use it to show what viscosity does to the circulation of a material loop. 3. Implement and test circulation as a line integral and as an area integral, a material loop advected exactly and by RK4, and the Biot–Savart law for straight filaments and a horseshoe vortex, with every result checked by an independent route. 4. Estimate the light aeroplane's starting and trailing vortex strength, the downwash at its tailplane, and decide when a laminar vortex model is (and is not) a sound basis for a decision. |
| **Prerequisites** | `P08.02` (vorticity, material derivative, pathlines), `P08.08` (circulation, point vortex, Kutta–Joukowski), `P08.06` (Navier–Stokes; the diffusion of momentum), `P01.11` and `P01.12` (line and surface integrals, Stokes's theorem), `P03.05` (RK4), `P03.04` (Gauss–Legendre and the periodic trapezoid rule). Skill check: you can compute ∮ u·dl round a circle in a given velocity field, and you can write one RK4 step. |
| **Engineering question** | When the light aeroplane (750 kg, 10 m span) starts to fly and reaches 60 m/s at 3000 m, how strong is the starting vortex it leaves behind, how quickly would viscosity alone spread and weaken it, and what downwash does the wing's vortex system induce at a tailplane 4.5 m behind? The answers: Γ₀ = 17.2 m²/s (horseshoe model of an elliptic loading; 13.5 m²/s as a 2-D mean); laminar diffusion would take 5.4 h to halve the circulation within 1 m (real wakes decay in minutes, by turbulence, so the laminar model fails for this decision); and the tail sees a downwash of 1.62 m/s, a downwash angle of 1.55°, more than the far-wake 1.33°. |
| **Assumptions and limits** | Incompressible flow; Kelvin's theorem for an inviscid, barotropic fluid with conservative body forces; the Lamb–Oseen vortex is axisymmetric, laminar and isolated; the wing's vortex system is a single horseshoe of span b′ = πb/4 (an elliptic loading lumped onto two trailing legs), straight, undeformed and with no core (points on a filament are excluded). Fails for turbulent wake vortices (which decay orders of magnitude faster than laminar diffusion predicts), for rolled-up and deforming wakes, and close to a filament, where a real vortex has a finite core. |
| **Analytical work** | ω = ∇×u, Γ = ∮u·dl and Stokes's theorem; Kelvin's theorem from the Euler equation; the viscous term ν∮∇²u·dl; the vorticity equation in plane flow and the Lamb–Oseen solution, its swirl, its circulation within r and the radius of peak swirl; the exact turning angle of a particle in a diffusing vortex (with the exponential integral); the Biot–Savart law for a straight segment, its infinite and semi-infinite limits, and the horseshoe downwash formulas; the starting vortex's downwash in a point-vortex model. |
| **Computational deliverable** | `starter/main.py`: `lamb_oseen_utheta`, `lamb_oseen_vorticity`, `lamb_oseen_circulation`, `vortex_velocity`, `circulation_polygon`, `disk_area_integral`, `turning_angle`, `material_loop_exact`, `advect_rk4`, `material_circulation_by_area`, `biot_savart_segment`, `horseshoe_velocity`, `point_vortices_velocity` (`core_radius2`, `circle_loop`, `observed_order` and `main` are given). Tests prove: the vorticity is the curl of the velocity and integrates to Γ, the field is divergence-free; the axis and far-field limits and the peak-swirl radius; the vorticity diffusion equation is satisfied; Stokes's theorem on three disks; the disk quadrature is exact for polynomials; circulation in irrotational flow depends only on what is enclosed (circles and a square); the turning angle matches quadrature and RK4 converges at order 4; Kelvin's theorem holds for a loop stretched more than fourfold, with a drift that falls as Δs²; viscosity lowers the circulation exactly as Lamb–Oseen predicts and line and area integrals agree on the deformed loop; Biot–Savart reproduces the infinite and semi-infinite lines, reverses with the filament and closes Stokes's theorem in 3-D; the horseshoe formulas, upwash ahead and outboard; a starting vortex makes the total circulation zero. |
| **Practical activity** | **computational experiment**: Kelvin's theorem demonstrated on a material loop advected through a vortex with and without viscosity, an RK4 convergence study against the exact particle paths, the laminar spreading of the starting vortex, and the horseshoe downwash along the centreline. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vorticity.png`, the proofs of Kelvin's theorem and of the Lamb–Oseen solution, the Biot–Savart derivations, the four experiment tables, and a half-page note on the starting vortex, its decay and the tail downwash. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: Kelvin's theorem is proved with its conditions; the Lamb–Oseen solution is verified by substitution; all public tests pass; the inviscid drift is shown to be discretisation error by refinement; the note states Γ₀, the tail downwash angle and why the laminar decay time cannot be used. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapter "Potential Flow Theory": Understanding Circulation and Vorticity, and chapter "Lifting Line & Finite Wing Theory": Helmholtz's Vortex Theorems; Theory of Vortex Lines (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/lifting-line-theory/) (required-free; link only). `@mit-ocw-2-25-advanced-fluid-mechanics` topic "Vorticity and Circulation" (https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/vorticity-and-circulation/) and the class note "Viscous Decay of a Free Vortex" (required-free: a viscous vortex solved by a similarity variable). `@mit-ocw-16-01-unified-engineering-fluids` 16.01–16.02 note F11 "Vorticity and Circulation" and 16.03–16.04 note F6 "3D Vortex Filaments, Lifting Line Theory" (required-free). `@libretexts-smyth-all-things-flow` 4.3 "Vorticity and Circulation"; 7.2 "Vortex dynamics in a homogeneous, inviscid fluid" (7.2.1 The vorticity equation; 7.2.2 Vortex filaments); 7.3 "Viscous effects" (7.3.1 The vorticity equation for a viscous fluid) (https://eng.libretexts.org/Bookshelves/Civil_Engineering/All_Things_Flow_-_Fluid_Mechanics_for_the_Natural_Sciences_(Smyth)/07%3A_Vortices/7.02%3A_Vortex_dynamics_in_a_homogeneous%2C_inviscid_fluid) (supplementary-free). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics* (paid; library), and Kundu & Cohen's *Fluid Mechanics* chapter on vorticity dynamics, which MIT 2.25 assigns (paid; library). |
| **Extensions** | ★ Replace the frozen vortex with two co-rotating Lamb–Oseen vortices and advect a material loop through their merger with RK4 and point insertion. ★ Build a lifting line: discretise the wing into many horseshoes, solve for their strengths, and compare the downwash with the single horseshoe. ★ Give each filament a finite core (a Lamb–Oseen-like cut-off in Biot–Savart) and compute the velocity a light aeroplane meets when it flies through an airliner's trailing vortex (Problem 3). |
| **Time** | unmeasured |

## 1. The engineering problem

`P08.08` showed that a wing's lift is ρVΓ per unit span, set by its circulation Γ. It left two questions open. Where does the circulation come from, if the air started at rest with none? And where does it go at the wing tips, where the lift must fall to zero? The answer to both is Kelvin's circulation theorem, which says that in an inviscid fluid the circulation round a loop of particles cannot change. A wing that acquires bound circulation +Γ must leave −Γ in the air: the **starting vortex** shed from the trailing edge as it starts to move. A finite wing's bound vortex cannot simply end at the tips (Helmholtz's theorems, Leishman "Lifting Line & Finite Wing Theory"), so it turns downstream as two **trailing vortices** that join the starting vortex in a closed loop. The trailing vortices induce a downwash at the tailplane, which sets the tail's angle of attack and so the aeroplane's trim; behind large aircraft they are a hazard to the aircraft that follow.

For the light aeroplane (750 kg, span 10 m, chord 1.5 m) at 60 m/s and 3000 m (ρ = 0.90913 kg/m³), the questions are: how strong is the starting vortex, how long would viscosity alone take to spread it out, and what downwash do its vortices induce at a tailplane 4.5 m behind the wing? Success criteria, written first: Stokes's theorem holds numerically between a line integral and an independent area integral; Kelvin's theorem holds for a strongly stretched material loop to a discretisation error that is shown to fall as Δs²; with viscosity the material loop's circulation follows the Lamb–Oseen prediction by two routes; RK4 converges at fourth order against exact particle paths; the Biot–Savart code reproduces the closed-form line-vortex limits.

## 2. Understand and derive

### 2.1 Vorticity, circulation and Stokes's theorem

**Vorticity** is the curl of the velocity, ω = ∇×u; in plane flow it has one component, ω = ∂v/∂x − ∂u/∂y, equal to twice the local angular velocity of a fluid particle (`P08.02`). **Circulation** round a closed curve C is Γ = ∮_C u·dl, counterclockwise positive in this lesson (the mathematical sense; `P08.08` counted clockwise as positive so that lift came out positive). **Stokes's theorem** (`P01.12`) connects them:

**Γ = ∮_C u·dl = ∬_A ω dA.**   (1)

Circulation counts the vorticity enclosed (Leishman, "Potential Flow Theory", Understanding Circulation and Vorticity; Smyth 4.3). A point vortex has ω = 0 everywhere except at its centre, so the circulation of any loop that encloses it is Γ, of any loop that does not is zero, whatever the loop's shape: the tests check a circle, a larger circle and a square.

### 2.2 Kelvin's theorem

Follow a closed loop C(t) made of fluid particles. Its circulation changes for two reasons: the velocity of each particle changes, and the loop's elements stretch, d(dl)/dt = du. So

dΓ/dt = ∮ (Du/Dt)·dl + ∮ u·du.

The second integral is ∮ d(|u|²/2) = 0. In an inviscid fluid, Euler's equation gives Du/Dt = −∇p/ρ + g. If the fluid is **barotropic** (p a function of ρ only, so ∇p/ρ = ∇P for a pressure function P) and the body force is **conservative** (g = −∇Φ), the first integral is −∮ d(P + Φ) = 0. Hence

**dΓ/dt = 0**   (Kelvin's circulation theorem)   (2)

(MIT 2.25, "Vorticity and Circulation"; Smyth 7.2.) Viscosity adds ν∮∇²u·dl = −ν∮(∇×ω)·dl, which is non-zero wherever the vorticity varies along the loop. A baroclinic fluid and a non-conservative force also generate circulation.

**The starting vortex.** Before the take-off run, a loop enclosing the wing and a large region of air round it has Γ = 0. Away from the boundary layers the air is inviscid, so the loop's circulation stays zero as it moves with the air. When the wing carries bound circulation (counterclockwise −Γ₀ in this lesson's sign, the aerodynamic clockwise Γ₀), the loop must also enclose +Γ₀ somewhere: in the vortex shed from the trailing edge as the Kutta condition (`P08.08`) is established. The code checks it with two point vortices, the wing's −17.168 m²/s and a starting vortex +17.168 m²/s 50 m behind: the circulation round a loop enclosing both is 0.00, round the wing alone −17.168 m²/s. While the starting vortex is near, it induces a downwash at the wing, Γ/(2πd) (−0.0546 m/s at 50 m), which reduces the wing's effective angle of attack. The lift therefore builds up over the first few chords of travel (Problem 6).

### 2.3 The Lamb–Oseen vortex

Take the curl of the Navier–Stokes equations (`P08.06`). In plane incompressible flow the pressure drops out and the vorticity obeys

**Dω/Dt = ν∇²ω.**   (3)

Vorticity is carried by the flow and diffuses like heat. For an axisymmetric vortex with purely azimuthal velocity u_θ(r, t), the convective term (u_θ/r)∂ω/∂θ vanishes and (3) is the diffusion equation in polar form, ω_t = ν(ω_rr + ω_r/r). Its solution starting from a concentrated vortex is a spreading Gaussian:

**ω = (Γ/(π r_c²)) e^{−r²/r_c²},   r_c² = r_c0² + 4νt,**   (4)

**u_θ = (Γ/(2πr))(1 − e^{−r²/r_c²}),   Γ(r) = Γ(1 − e^{−r²/r_c²}).**   (5)

This is the **Lamb–Oseen vortex**, an exact solution of the Navier–Stokes equations (verified by substitution in Problem 2 and in `test_lamb_oseen_satisfies_the_vorticity_diffusion_equation`). Near the axis it rotates as a solid body, u_θ ≈ (ω(0)/2) r; far outside its core it is the point vortex Γ/(2πr); the swirl peaks at r = 1.1209 r_c; the core grows as √(νt). With ν = 0 (r_c fixed) the same profile is a steady solution of Euler's equations, a vortex with a finite core, and that is the field used for the inviscid Kelvin demonstration. MIT 2.25's class note "Viscous Decay of a Free Vortex" solves a closely related problem by the same similarity variable r/√(νt).

### 2.4 A material loop in a vortex

Every particle of a Lamb–Oseen vortex moves on a circle. So a material loop can be advected **exactly**: a particle at radius r turns through the angle ∫₀ᵗ u_θ(r, s)/r ds. With ν = 0 that is Ω(r)t. With ν > 0 the integral of e^{−r²/r_c²(s)} can be done by substituting w = r²/r_c²(s), which gives, with E₁ the exponential integral,

θ(r, t) = (Γ/(2πr²))[t − (r²/4ν)(g(w₀) − g(w_t))],   g(w) = E₁(w) − e^{−w}/w,   w₀ = r²/r_c0²,  w_t = r²/r_c²(t).   (6)

(The first version of this formula in the code had the two g terms the wrong way round; the test against `scipy.integrate.quad` found it.) Because the flow only turns each circle rigidly, the material region keeps, at every radius, the same arc length it started with, and ω depends on r only. The circulation of the material loop at time t is therefore, by (1), **the integral of the current vorticity ω(r, t) over the initial disk**. That is an area integral entirely independent of the deformed loop: `material_circulation_by_area`. Without viscosity ω does not change and the circulation is constant (Kelvin); with viscosity it falls as the vorticity diffuses out of the loop.

### 2.5 The Biot–Savart law and the horseshoe vortex

In 3-D, the velocity induced by a vortex filament of strength Γ is dV = (Γ/4π) dl × r/|r|³ (the same law as the magnetic field of a current; Leishman, "Lifting Line & Finite Wing Theory", Theory of Vortex Lines). For a straight segment from A to B, with r₀ = B − A, r₁ = P − A, r₂ = P − B,

**V = (Γ/4π) (r₁ × r₂)/|r₁ × r₂|² [r₀·(r₁/|r₁| − r₂/|r₂|)].**   (7)

For an infinite line it gives Γ/(2πd); for a semi-infinite one, at a point abeam its end, Γ/(4πd). A wing's simplest vortex system is a **horseshoe**: a bound vortex along the span (its vector along +y, so that L = ρV × Γ points up) whose ends turn downstream into two trailing legs. For an elliptic lift distribution the legs are b′ = πb/4 apart, and the centre-span circulation is Γ₀ = L/(ρVb′). At the centre of the bound vortex the legs induce **−Γ₀/(πb′)**; far downstream, where they look infinite, **−2Γ₀/(πb′)**. The downwash far behind is twice the downwash at the wing.

**The aeroplane.** b′ = 7.854 m and **Γ₀ = 7355.0/(0.90913 × 60 × 7.854) = 17.168 m²/s** (the mean 2-D value, L/(ρVb), is 13.484 m²/s). The horseshoe gives a downwash of 0.6958 m/s at the wing (0.664°), **1.619 m/s at the tail 4.5 m behind (1.546°)** and 1.392 m/s far downstream (1.329°). The tail sees more than the far wake because the bound vortex, seen from behind, adds its own downwash there. The elliptic-wing estimate of the far downwash angle, 2C_L/(πAR) with C_L = 0.2996 and AR = 6.667, is 1.639°; the single horseshoe gives 1.329°. The 19 % difference is the cost of lumping all the trailing vorticity into two lines (Extension 2).

### 2.6 Where the models stop

- **Turbulence.** A real trailing or starting vortex is turbulent. Its effective viscosity is many orders of magnitude above ν, and it decays by atmospheric turbulence, by instability of the vortex pair and near the ground. Laminar Lamb–Oseen numbers (Experiment 3) are lower bounds on spreading, not predictions.
- **Roll-up.** The trailing vorticity leaves the whole span as a sheet and rolls up into two cores over a few spans; the horseshoe assumes it is concentrated from the start and does not deform.
- **Cores.** Biot–Savart's straight filament has infinite speed on its own line. The code returns zero there (no self-induction), and any point closer to a filament than a core radius is outside the model.
- **Kelvin's conditions.** Inside the boundary layers and the wake sheet viscosity acts; Kelvin applies to loops that stay in the inviscid outer flow, which is why the loop round "the wing" must be taken outside its boundary layer.

## 3. Build the computational model

Open `starter/main.py`. `core_radius2`, `circle_loop`, `observed_order` and `main` are given.

### Step 1: the vortex

`lamb_oseen_utheta`, `lamb_oseen_vorticity` and `lamb_oseen_circulation` are (4)–(5). Use `-np.expm1(-r**2/rc2)` for 1 − e^{−x}, which keeps full precision near the axis, and the series Γr/(2π r_c²) for r² < 10⁻⁸ r_c² to avoid 0/0. `vortex_velocity(xy, Gamma, rc2, center)` turns u_θ into Cartesian components (N, 2).

### Step 2: two routes to the circulation

`circulation_polygon(loop, velocity)` is the midpoint rule on the polygon's chords, second order in the segment length:

```python
nxt = np.roll(loop, -1, axis=0)
mid = 0.5 * (loop + nxt)
return float(np.sum(velocity(mid) * (nxt - loop)))
```

`disk_area_integral(f, center, radius, nr, nt)` uses Gauss–Legendre points in r (with the Jacobian r) and the periodic trapezoid rule in θ: spectrally accurate for smooth integrands, and exact for low-degree polynomials (a test).

### Step 3: the material loop

`turning_angle` is Ω(r)t for ν = 0 and (6) for ν > 0 (`scipy.special.exp1`); `material_loop_exact` rotates each point of the initial loop by its angle. `advect_rk4` integrates dx/dt = u(x, t) with classical RK4, the velocity evaluated with the core size at each stage time. `material_circulation_by_area` is the area integral of ω(r, t) over the initial disk.

### Step 4: filaments

`biot_savart_segment(P, A, B, Gamma)` is (7), vectorised over the points P, returning zero where |r₁ × r₂| vanishes. `horseshoe_velocity` adds three segments, the legs ending at a large `length` that stands for infinity. `point_vortices_velocity` superposes plane point vortices.

Running `main()` prints (abridged)

```
  mean 2-D circulation L/(rho V b) = 13.484 m^2/s; horseshoe (b' = pi b/4 = 7.8540 m): Gamma0 = 17.168 m^2/s
  Kelvin, plane model: circulation round both vortices 0.00e+00 m^2/s; round the wing alone -17.1682 m^2/s
  starting vortex 50 m behind induces -0.0546 m/s at the wing (2-D, Gamma/(2 pi d))
air at 3000 m: mu = 1.6937e-05 Pa s, nu = 1.8630e-05 m^2/s
  t =      1 s: core radius of peak swirl 0.97 cm, peak swirl 201.99 m/s, circulation within 1 m 100.0000 % of Gamma0
  t =     60 s: core radius of peak swirl 7.50 cm, peak swirl 26.08 m/s, circulation within 1 m 100.0000 % of Gamma0
  t =    600 s: core radius of peak swirl 23.70 cm, peak swirl 8.25 m/s, circulation within 1 m 100.0000 % of Gamma0
  laminar diffusion takes 19360 s (5.4 h) to halve the circulation within 1 m
  horseshoe downwash on the centreline      4.5 m behind the bound vortex: w = -1.6193 m/s, angle 1.546 deg
Kelvin's theorem, material loop: circle radius 0.6 centred at (0.9, 0); Gamma = 1, rc0 = 1
  Stokes at t = 0: line integral 0.1538153197, area integral 0.1538153786
  RK4 against the exact particle paths over one turnover: errors [4.98e-05 3.17e-06 1.99e-07 1.25e-08], orders [3.98 3.99 3.99]
  RK4-advected loop (400 steps, 2 turnovers): circulation 0.15381505, change -2.7e-07
```

and writes `vorticity.png`: the material loop at 0, 1, 2 and 4 turnovers, its circulation against time with and without viscosity by both routes, and the horseshoe downwash along the centreline.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-09-vorticity-circulation`. The tests:

1. **Vortex**: finite-difference curl of the velocity equals (4); divergence below 10⁻⁶; ∫ω dA = Γ by quadrature; zero swirl on the axis, the point vortex at r = 20, a maximum at 1.1209 r_c; (4) satisfies the polar diffusion equation by finite differences in r and t.
2. **Stokes**: line and area integrals agree to 10⁻⁵ on three disks (the midpoint rule's O((2π/n)²) error on 2000 chords) and match (5) on a concentric circle; the disk quadrature integrates polynomials to 10⁻¹³.
3. **Irrotational flow**: the circulation of a point vortex is the same on two circles and a square and zero on a loop that misses it.
4. **Material loops**: (6) matches `quad` to 10⁻¹⁰; RK4 converges at order 4 ± 0.3; particles keep their radius; the inviscid circulation after 120 time units changes by under 2 × 10⁻⁵ of itself while the loop stretches more than fourfold, the change falls by 4 when the points double, and an RK4-advected loop gives the same; the viscous concentric loop follows (5) and the off-centre loop's line and area integrals agree to 10⁻⁵, both below 80 % of the initial value.
5. **Biot–Savart**: the infinite and semi-infinite lines, reversal of the filament, zero on its own line; the circulation of the induced field round the filament equals Γ (Stokes in 3-D) and is zero round a loop that misses it; the horseshoe's −Γ/(πb′), −2Γ/(πb′), upwash ahead and outboard.
6. **Starting vortex**: total circulation zero, the wing's alone −Γ, and the induced downwash fading with distance.
7. **main** writes the figure.

### Compare

The independent routes are built into every result. Stokes's theorem is tested by comparing a polygon line integral with a Gauss–Legendre area integral that never touches the loop; Kelvin's theorem by an exact particle map and by RK4, which agree to 10⁻⁷; the viscous decay by the closed form (5) and by two numerical integrals on a loop deformed out of all recognition. `scipy.integrate.quad` is the established tool that checked (6), and it found a sign error in the first version. For the 3-D filaments, the infinite-line formula and the line integral of the induced field round the filament are the checks.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Kelvin's theorem.** A circle of radius 0.6 centred at (0.9, 0) in the vortex Γ = 1, r_c0 = 1 (4000 points), advected for 0, 1, 2 and 4 turnovers (the turnover time at r = 0.9 is 57.6), first with ν = 0, then with ν = 0.01. Predict both circulation histories and the loop's length. *Finding:* inviscid, the line integral reads 0.15381532, 0.15381525, 0.15381505, 0.15381423 while the area integral stays at 0.15381538 and the loop's length grows 3.77 → 10.1 → 19.2 → 37.7; viscous, 0.15381532 → 0.08185026 → 0.05406037 → 0.03202916 (area integrals 0.08185032, 0.05406042, 0.03202920), with the loop growing only to 6.42 because the spreading vortex shears it less. The inviscid drift (7 × 10⁻⁶ of Γ at four turnovers) falls by a factor of 4.00 per doubling of the points (Problem 4): discretisation, not physics.
2. **RK4 against the exact paths.** 20 particles of the loop over one turnover with 25–200 steps. Predict the order. *Finding:* errors 4.98 × 10⁻⁵, 3.17 × 10⁻⁶, 1.99 × 10⁻⁷, 1.25 × 10⁻⁸, orders 3.98, 3.99, 3.99. With 400 steps over two turnovers the RK4-advected loop's circulation differs from the initial value by 2.7 × 10⁻⁷, the same as the exact map's.
3. **The starting vortex, laminar.** Γ₀ = 17.168 m²/s from a zero core in air at 3000 m (ν = 1.863 × 10⁻⁵ m²/s) at t = 1, 10, 60, 600 s. Predict the core size at one minute. *Finding:* the radius of peak swirl is 0.97, 3.06, 7.50 and 23.7 cm and the peak swirl 202, 63.9, 26.1 and 8.25 m/s; the circulation within 1 m stays at 100.0000 % for ten minutes, and halving it takes 5.4 h. Real starting and trailing vortices are gone in minutes: the laminar model is the wrong tool for any decision about wake decay.
4. **Downwash along the centreline.** The horseshoe from 0 to 10 km behind the bound vortex. Predict where the downwash is largest. *Finding:* −0.6958 m/s at the wing, −1.6193 at 4.5 m, −1.3937 at 50 m and −1.3916 m/s far behind (the formulas give −0.6958 and −1.3916). The downwash overshoots its far value just behind the wing, where the bound vortex still contributes, which is where the tail sits.

Raw data format: CSV tables `t, nu, gamma_line, gamma_area, loop_length` (1); `nt, max_position_error` (2); `t_s, r_peak_cm, u_peak_m_s, frac_circ_1m` (3); `x_m, w_m_s, eps_deg` (4); plus `vorticity.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| Stokes at t = 0 | line = area | 0.1538153197 vs 0.1538153786 | 6 × 10⁻⁸ | midpoint rule, O(Δs²) |
| Kelvin, ν = 0, four turnovers | constant | drift 1.1 × 10⁻⁶ (7.1 × 10⁻⁶ relative) | → 0 as Δs² | discretisation |
| viscous loop | line = area integral of ω(t) | agree to 6 × 10⁻⁸ | — | — |
| RK4 particle paths | order 4 | 3.98, 3.99, 3.99 | — | — |
| downwash at the bound vortex | −Γ₀/(πb′) = −0.6958 m/s | −0.6958 | < 10⁻⁴ | — |
| far downwash | −2Γ₀/(πb′) = −1.3916 m/s | −1.3916 (at 10 km) | < 10⁻⁴ | finite leg length |
| total circulation with starting vortex | 0 | 0.00 | — | — |

Every model result agrees with its independent check. The engineering answer: **the starting vortex has the wing's circulation, 17.2 m²/s at mid-span (13.5 m²/s as a 2-D mean); the tail 4.5 m behind sees a downwash of 1.62 m/s, a 1.55° downwash angle, which the tail's incidence must allow for.** The laminar decay time of 5.4 hours is a correct answer to the wrong model and must not be used: real vortices decay by turbulence. The tail downwash is probably within about 20 % (the horseshoe against the elliptic-wing estimate, 1.33° against 1.64° far downstream). A lifting-line model (Extension 2) or a flight measurement of the tail's angle of attack would settle it.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the proofs of (2) and of (4)–(5) and the derivation of the horseshoe formulas; `vorticity.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the starting-vortex strength, the tail downwash angle and its uncertainty, and why the laminar decay time is not a basis for a decision.

### Adversarial self-checks

- **What if the aeroplane flies at half the speed?** At the same weight Γ₀ doubles to 34.3 m²/s; the downwash velocity doubles and the downwash angle w/V quadruples, to about 6.2° at the tail. Slow flight, with flaps, is where the tail's setting is most demanding.
- **What if viscosity goes to zero?** The Lamb–Oseen core stops growing and Kelvin's theorem holds exactly; the starting vortex would persist forever. That does not happen, which shows that diffusion and turbulence do the work.
- **Which assumption fails first in the real world?** Laminar flow in the vortex cores; then the straight, undeformed horseshoe (the wake rolls up and descends at Γ₀/(2πb′)).
- **What measurement could falsify the model?** A vane or pitot rake at the tail measuring the flow angle at several lift coefficients: the downwash angle should be proportional to C_L. A different slope measures how good the horseshoe's b′ is; a nonlinearity means separation or a rolled-up wake.
- **Does the result survive a change of method?** Exact particle maps and RK4 give the same circulation to 10⁻⁷; line and area integrals agree on every loop; the horseshoe matches its closed forms.
- **Which decision changes if the uncertainty doubles?** If the tail downwash were uncertain by 40 % instead of 20 %, ±0.6° on the tail's angle of attack, the tail incidence could not be set from analysis alone and an adjustable trim surface or a flight test would become necessary.

## Check yourself

1. Define vorticity and circulation, and state the theorem that connects them.
2. State Kelvin's circulation theorem and its three conditions.
3. Why must a wing that starts to fly leave a starting vortex behind?
4. What is the vorticity of a point vortex away from its centre, and the circulation round any loop that encloses it?
5. How does the core of a Lamb–Oseen vortex grow, and what happens to the circulation far outside the core?
6. State the Biot–Savart result for a straight filament and its infinite-line limit.
7. Why is a horseshoe vortex's downwash far behind the wing twice its value at the wing?
8. In the material-loop experiment the inviscid circulation changed slightly. How do you know it was not a failure of Kelvin's theorem?
