# P10.17  Propeller & Rotor Theory

> Momentum theory tells you the least power any propeller can need; blade element momentum theory tells you how much more yours will need, and why.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the actuator-disk results from conservation of mass, momentum and energy: the induced velocity, the ideal power, the static power T^{3/2}/√(2ρA) and the ideal efficiency η_i = 2/(1 + √(1 + 8C_T/(πJ²))). 2. Derive the blade element momentum (BEM) equations for an annulus, with axial and swirl induced velocities and Prandtl's tip-loss factor, and reduce them to a single residual in the inflow angle. 3. Implement and test a BEM code two ways (relaxed fixed-point iteration and a bracketing solve) and show that its ideal limit reproduces momentum theory exactly. 4. Compute C_T, C_P, C_Q and η against advance ratio for a small-UAV propeller and split its losses into induced, swirl, tip and profile parts. 5. Decide which propeller pitch suits a cruise condition and say what the section model's uncertainty does to the answer. |
| **Prerequisites** | `P10.04` (lift-curve slope and zero-lift angle of a section), `P10.10` (trailing vortices and induced velocity on a lifting line), `P03.02` (bisection and Brent's method), `P03.11` (observed order of convergence), `P00.14`. From Phase 08 (recapped in §2.1): `P08.03` (control-volume momentum balance) and `P08.05` (dimensionless groups). Skill check: you can write the momentum and energy balance for a control volume with one inlet and one outlet. |
| **Engineering question** | At what speed must the reference UAV's 10 × 4.5 in propeller turn to give the 1.20 N cruise thrust at 12 m/s, how efficient is it there, and how far is that below the momentum-theory ideal? The answer, for the stated blade geometry and section model: 6203 rpm (J = 0.457), η = 0.676 against an ideal 0.940, drawing 21.4 W of shaft power; profile drag of the low-Reynolds-number sections costs 0.24 of efficiency, tip loss 0.013 and swirl 0.003. A 6 in pitch would do better in cruise (η = 0.753). |
| **Assumptions and limits** | Steady, incompressible, axisymmetric inflow along the axis (no incidence of the propeller disk, no fuselage blockage); annuli independent of each other (no radial flow); 2-D section data apply at each radius; Prandtl's tip-loss factor models the finite blade number; the stated section model replaces measured section data at Re ≈ 6 × 10⁴; tip Mach number 0.24, so compressibility is neglected. Momentum theory fails in the vortex-ring state and when the wake would move upstream (descent through the rotor's own wake); BEM inherits that. |
| **Analytical work** | The actuator-disk balance and its four results; the propeller coefficients and their relation C_P = 2π C_Q, η = J C_T/C_P; the annulus balance with tip loss; the inflow-angle residual R(φ); the energy identity dP = dT (V + u_a) for an ideal annulus and the proof that uniform inflow gives the Froude efficiency. |
| **Computational deliverable** | `starter/main.py`: `actuator_disk`, `ideal_efficiency`, `static_ideal_power`, `prandtl_tip_loss`, `element_residual`, `bem`, `bem_fixed_point`, `design_uniform_inflow_blade`, `rpm_for_thrust` (`Section`, `Blade`, `reference_propeller`, `annuli`, `uav_cruise_thrust`, `sweep_advance_ratio`, `observed_order` and `main` are given). Tests prove: momentum and energy are conserved in the actuator disk; the ideal-efficiency formula agrees with the disk; the BEM residual vanishes at the solution and each annulus satisfies its momentum balance; a blade designed for uniform inflow reproduces momentum theory exactly in the ideal limit; no blade beats momentum theory; each loss lowers the efficiency; coefficients collapse at fixed J; the two solvers agree; the element study converges at about second order. |
| **Practical activity** | **computational experiment**: BEM against actuator-disk momentum theory (the independent check), a loss budget at the cruise point, a sweep in advance ratio, a pitch comparison and convergence studies in the number of elements and the relaxation factor. No propeller is run. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. This lesson involves no running propeller: do not test propellers or rotors at home. A spinning 10 in propeller at 6000 rpm has a tip speed of 80 m/s and cuts skin and eyes; the course never requires an unguarded rotor. |
| **Evidence** | Completed `main.py`, passing tests, `propeller_bem.png`, the derivations of §2, the five experiment tables of §4 and a half-page note recommending a propeller pitch for the UAV's cruise, with the uncertainty that matters most. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: η_i is derived from the disk balance; all public tests pass; the uniform-inflow check reproduces momentum theory to 10⁻⁶; the loss budget is explained term by term; the note states a pitch and the assumption that most affects it. |
| **References** | `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics` (Stepniewski, *Rotary-Wing Aerodynamics* Vol. I, https://ntrs.nasa.gov/citations/19790013868) Chapter II "Momentum Theory": 3. Actuator Disc; 3.1 Induced Velocity and Thrust in Axial Translation; 3.3 Ideal Power in Climb and Hovering; Chapter III "Blade Element Theory": 2.1 Basic Considerations of Thrust and Torque Predictions; 2.2 Combined Blade-Element and Momentum Theory; 2.3 Nondimensional Coefficients; 2.5 Tip Losses (required-free). `@leishman-intro-aerospace-flight-vehicles` "Aircraft Propellers" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/propellers/): Propeller Fundamentals; Momentum Theory of the Single Propeller; Propeller Performance Coefficients; General Propeller Performance (required-free, link only). `@nasa-glenn-beginners-guide-aeronautics` "Propeller Analysis" (https://www.grc.nasa.gov/www/k-12/airplane/propanl.html) (supplementary-free: simple momentum theory). `@openvsp-docs` "Definitions of Aerodynamic Coefficients" (https://openvsp.org/wiki/doku.php?id=coeffdefn) (supplementary-free: C_T, C_Q, C_P, J, η conventions). `@ntrs-naca-tr-640-full-scale-propellers` (supplementary-free: measured full-scale propeller charts for Extension 2). `@scipy-userguide` "Optimization (scipy.optimize)" > "Root finding" > "Scalar functions" (supplementary-free: `brentq`). Optional textbooks for deeper study: J. G. Leishman, *Principles of Helicopter Aerodynamics* (BEMT with non-uniform inflow and its validation); Bertin & Cummings, *Aerodynamics for Engineers* (worked propeller examples). |
| **Extensions** | ★ Add a hub-loss factor (the same form as the tip loss, with r − r_hub) and see how much thrust the inner blade loses. ★ Digitise one curve of the free full-scale propeller charts in NACA TR-640, model that blade, and compare C_T(J) and η(J); explain the difference. ★ Design the minimum-induced-loss blade for the cruise point (constant u_a with tip loss, Betz's condition) and compare its η with the reference blade at the same thrust. |
| **Time** | unmeasured |

## 1. The engineering problem

The reference UAV (span 1.5 m, chord 0.20 m, 1.2 kg) cruises at 12 m/s on an electric motor and a 10 × 4.5 in two-blade propeller: diameter 0.254 m, nominal geometric pitch 0.114 m, the distance the blade would advance in one turn if it screwed through the air like a bolt through a nut. The battery holds a fixed amount of energy. Every watt the propeller wastes is endurance lost, and the propeller is usually the least efficient link between the battery and the thrust: a typical small motor converts 80 % of its electrical power, and a small propeller only 60–75 % of the shaft power into thrust power. To pick a propeller, or to know how fast it must turn, you need its thrust and torque as functions of speed and rotational speed.

Two models answer at two levels. **Momentum theory** treats the propeller as an ideal disk that adds momentum to the air passing through it. It ignores everything about the blades and so gives a bound: the *least* power any propeller of that diameter can absorb to make that thrust. **Blade element momentum theory** (BEM) cuts the disk into rings, lets each ring obey momentum theory, and makes the blade element in it produce lift and drag from its section data. It predicts what a particular blade does, and it must collapse to momentum theory when the blade is ideal. That collapse is this lesson's independent check.

Success criteria, written first: the actuator disk conserves momentum and energy to round-off; the BEM residual is below 10⁻⁹ at every element; a blade designed for uniform inflow, run with no drag, no swirl and no tip loss, reproduces the momentum-theory thrust, inflow and efficiency to 10⁻⁶; no blade beats the momentum-theory efficiency at its own thrust; the thrust converges with the number of elements at about second order.

## 2. Understand and derive

### 2.1 Recap: the momentum balance (`P08.03`)

For a steady flow through a control volume, the net force on the fluid equals the momentum flux out minus the momentum flux in; the power put into the fluid equals the kinetic-energy flux out minus in (no heat, no change of height). Coefficients make results independent of size (`P08.05`): a propeller with diameter D turning at n revolutions per second is described by the **advance ratio** J = V/(nD), the distance flown per revolution in diameters, and

C_T = T/(ρ n² D⁴),  C_Q = Q/(ρ n² D⁵),  C_P = P/(ρ n³ D⁵) = 2π C_Q,  η = TV/P = J C_T/C_P.   (1)

(Rotor engineers use C_T = T/(ρA(ΩR)²) instead; lesson `P10.18` does. The conventions differ only in the reference quantities.)

### 2.2 The actuator disk

Replace the propeller by a thin disk of area A that raises the pressure of the air passing through it without changing its speed across the disk. Far ahead, the air moves at V relative to the disk; at the disk it moves at V + v (v is the **induced velocity**); far behind, in the contracted **slipstream**, it moves at V + w.

![An actuator disk of area A, with a pressure jump Δp, inside a stream tube. Far ahead the air moves at V, at the disk at V + v and far behind, in the contracted slipstream, at V + w. Below, the axial speed along the tube rises smoothly from V through V + v at the disk to V + w far behind.](figures/actuator-disk.svg)
*Figure 1. The actuator disk. The speed is continuous through the disk and the pressure is not; the slipstream contracts because the air keeps accelerating behind the disk, and w = 2v (equation (2)). Drawn for heavy loading, v = V/2, so the contraction is visible.*

*Mass:* ṁ = ρA(V + v). *Momentum:* T = ṁ w. *Energy:* the power absorbed equals the rise of kinetic-energy flux, P = ½ṁ((V + w)² − V²) = ½ṁ w(2V + w). *Work at the disk:* the disk pushes with force T on air passing through it at V + v, so it does work at the rate P = T(V + v). Equating the two expressions for P with T = ṁ w gives

**w = 2v**: half the final velocity increase happens ahead of the disk.   (2)

Then T = 2ρA(V + v)v, a quadratic for v:

**v = −V/2 + √(V²/4 + T/(2ρA)),  P_ideal = T(V + v),  η_i = TV/P = V/(V + v).**   (3)

At V = 0 (static thrust, hover), v_h = √(T/(2ρA)) and

**P_ideal = T^{3/2}/√(2ρA).**   (4)

Doubling the thrust costs 2^{3/2} = 2.83 times the power; doubling the disk area saves a factor √2. That is why helicopters have big rotors and why a small drone's hover is expensive.

**In propeller coefficients.** Write a = v/V. Then T = 2ρ(πD²/4)V²a(1 + a) and with V = JnD, C_T = (π/2)J²a(1 + a), so a = [−1 + √(1 + 8C_T/(πJ²))]/2 and

**η_i = 1/(1 + a) = 2/(1 + √(1 + 8C_T/(πJ²))).**   (5)

Light loading (C_T/J² small) gives η_i → 1; heavy loading gives η_i → 0.

**The UAV at cruise.** The airframe's polar is assumed (C_D0 = 0.035, e = 0.8, AR = 7.5): C_L = 0.445 and the thrust required is **1.204 N**. For the 10 in disk, (3) gives v = 0.760 m/s, a slipstream of 13.52 m/s, an ideal power of 15.36 W and **η_i = 0.940**. That is the bound.

### 2.3 Blade element momentum theory

Cut the disk into annuli of radius r and width dr. At radius r the blade section sees the axial velocity V + u_a and the tangential velocity Ωr − u_t, where u_a is the axial induced velocity and u_t the **swirl**, the rotation the blade leaves in the air (it turns with the blade, so the blade sees less tangential speed). The **inflow angle** φ and the angle of attack α follow from the velocity triangle:

![Velocity triangle of a blade element at radius r, seen from the blade tip, with the flight direction up and the blade moving left. The tangential component Ωr − u_t and the axial component V + u_a add to the relative velocity W, which meets the plane of rotation at the inflow angle φ. The chord line is at the blade pitch β to the plane of rotation, so the angle of attack is α = β − φ. Lift dL is perpendicular to W and drag dD lies along it.](figures/blade-element-velocity-triangle.svg)
*Figure 2. The velocity triangle of (6). Resolving dL and dD along the axis and the plane of rotation gives (7): lift contributes cos φ to thrust and sin φ to torque. Angles are drawn larger than on a cruising propeller.*

tan φ = (V + u_a)/(Ωr − u_t),  α = β − φ,  W² = (V + u_a)² + (Ωr − u_t)².   (6)

*Blade element.* B blades of chord c produce lift and drag per unit span ½ρW²c·(cl, cd), perpendicular and parallel to W. Resolved along the axis and the rotation,

dT/dr = B c ½ρW² (cl cos φ − cd sin φ),  dQ/dr = B c ½ρW² (cl sin φ + cd cos φ) r.   (7)

*Momentum on the annulus.* Applying (2) to the annulus of area 2πr dr, and the angular-momentum balance (the far wake swirls at 2u_t),

dT/dr = 4πrρF(V + u_a)u_a,  dQ/dr = 4πr²ρF(V + u_a)u_t,   (8)

where F is **Prandtl's tip-loss factor**,

F = (2/π) arccos[exp(−B(R − r)/(2r sin φ))].   (9)

A real propeller has B discrete blades, not a disk: near the tip, air leaks round the blade from the high-pressure to the low-pressure side, the load falls to zero at the tip, and the momentum the annulus receives is less than the disk value. F → 0 at the tip and F → 1 inboard or as B → ∞. It is a correction from Prandtl's model of the helical wake as a stack of sheets, not a law; it is very good for lightly loaded propellers.

**One residual per element.** Equate (7) and (8) and write σ′ = Bc/(2πr) (the local solidity), C_n = cl cos φ − cd sin φ, C_t = cl sin φ + cd cos φ:

k = σ′C_n/(4F sin²φ),  u_a = Vk/(1 − k);  k′ = σ′C_t/(4F sin φ cos φ),  u_t = Ωr k′/(1 + k′).   (10)

Substituting into the velocity triangle (V + u_a) cos φ = (Ωr − u_t) sin φ gives a single equation in φ:

**R(φ) = sin φ (1 − k) − (V/(Ωr)) cos φ (1 + k′) = 0.**   (11)

R has no division by V, so the static case (V = 0, root at k = 1) needs no special code. For a thrusting element R < 0 as φ → 0 (k grows like 1/sin φ) and R > 0 at φ = π/2 (C_n = −cd < 0 there), so a root is bracketed and bisection cannot fail. This formulation follows the published "guaranteed convergence" idea of solving each element in φ alone (it is written here from the equations above).

**The classical iteration** instead guesses (u_a, u_t), computes φ, α and the loads, solves (8) for new values and relaxes: u ← u + ω(u_new − u). It is simple and usually works; §4 shows where it fails.

### 2.4 The ideal limit: BEM must reproduce momentum theory

Set cd = 0, F = 1 and u_t = 0. Then (7) gives dQ/dr ⋅ Ω = dT/dr ⋅ Ωr tan φ = dT/dr (V + u_a): each annulus absorbs exactly **dP = dT (V + u_a)**, the disk result (3) ring by ring. If u_a is the same on every ring, the whole blade has P = T(V + u_a) and η = V/(V + u_a), the Froude efficiency at that thrust. If u_a varies, the same thrust costs more power (at fixed total momentum flux, kinetic energy is least when the velocity increase is uniform), so η < η_i. Hence:

- a blade *designed* for uniform inflow must reproduce (3) exactly in the ideal limit;
- no blade can beat (5).

`design_uniform_inflow_blade` builds such a blade: choose u_a from (3) on the annulus area π(R² − r_hub²), give every ring the thrust (8) with that u_a, and size the chord at a design lift coefficient (c = (dT/dr)/(B ½ρW² cl_d cos φ)), with pitch β = φ + α(cl_d).

### 2.5 The blade and the section (stated assumptions)

The reference propeller's real geometry is not published, so the lesson **assumes** a smooth one: chord c/R = 0.045 + 0.55x(1 − x)^{1.3} with x = r/R (largest, 0.16R, near mid-radius; 0.045R at the tip), constant geometric pitch 0.114 m so β(r) = arctan(p/(2πr)) (8.1° at the tip), hub at 0.15R. The sections are thin and cambered, at Re ≈ 6 × 10⁴ (`P10.14`), where 2-D data scatter widely; the **section model** is cl = cl_max tanh(a₀(α − α₀)/cl_max) with a₀ = 5.7 /rad, α₀ = −4°, cl_max = 1.1, and cd = 0.020 + 0.025(cl − 0.3)² plus a post-stall term 1.8 sin²(excess angle). These are modelling choices, not data. Experiment 5 shows that the drag level matters most.

### 2.6 Where the model stops

- **Annulus independence.** Each ring is solved alone; real wakes have radial flow, and heavily loaded rotors (hover) have a strongly contracting, non-uniform wake.
- **Section data.** At Re ≈ 6 × 10⁴ the section behaviour depends on laminar separation bubbles (`P10.14`); a ±50 % change in cd0 moves η by about ±0.09 (Experiment 5).
- **Inflow.** An aircraft propeller sees the fuselage's blockage and, in a climb or turn, an inclined flow; BEM here assumes axial, uniform inflow.
- **Momentum breakdown.** Windmilling at large J and the vortex-ring state in descent violate the slipstream picture.

## 3. Build the computational model

Open `starter/main.py`. `Section`, `Blade`, `reference_propeller`, `annuli`, `uav_cruise_thrust`, `sweep_advance_ratio`, `observed_order` and `main` are given.

### Step 1: momentum theory

`actuator_disk(T, V, A)` returns v, P, η, ṁ and the wake speed from (2)–(3); `ideal_efficiency(CT, J)` is (5); `static_ideal_power(T, A)` is (4).

### Step 2: tip loss and the element residual

`prandtl_tip_loss(r, R, B, phi)` is (9) (guard sin φ away from zero and clip the exponential into [0, 1]). `element_residual(phi, r, c, beta, blade, V, Omega, ...)` evaluates α = β − φ, the section coefficients, C_n, C_t, F, k, k′ and R(φ) of (10)–(11) for all elements at once.

### Step 3: the robust solver

`bem(blade, V, n_rps, ...)` brackets every element's φ in (10⁻⁶, π/2] and bisects all of them together:

```python
lo, hi = np.full_like(r, 1e-6), np.full_like(r, np.pi / 2)
R_lo = element_residual(lo, ...)[0]
while np.max(hi - lo) > tol:
    mid = 0.5 * (lo + hi)
    R_mid = element_residual(mid, ...)[0]
    same = np.sign(R_mid) == np.sign(R_lo)
    lo, R_lo, hi = np.where(same, mid, lo), np.where(same, R_mid, R_lo), np.where(same, hi, mid)
```

Then u_t = Ωr k′/(1 + k′), u_a = (Ωr − u_t) tan φ − V, and the loads (7) summed over the annuli give T, Q, P and the coefficients (1). The switches `tip_loss`, `swirl`, `drag` turn each loss off for the ideal-limit test and the loss budget.

### Step 4: the classical iteration

`bem_fixed_point(...)` relaxes (u_a, u_t) as described in §2.3, solving (8) for u_a = −V/2 + √(V²/4 + (dT/dr)/(4πrρF)). Keep it: it is the textbook method, and the comparison of the two solvers is part of the experiment.

### Step 5: the uniform-inflow blade and the operating point

`design_uniform_inflow_blade(T, V, n_rps, R, ...)` is §2.4 with an exactly linear section cl = 2πα and no drag. `rpm_for_thrust(blade, T_req, V)` finds n with `brentq` on bem(...)["T"] − T_req.

Running `main()` prints (abridged):

```
UAV cruise: C_L = 0.445, thrust required = 1.204 N (assumed CD0 = 0.035, e = 0.8)
actuator disk at 12 m/s: v = 0.760 m/s, wake 13.520 m/s, ideal power 15.36 W, ideal efficiency 0.9404
BEM cruise: n = 103.38 rev/s = 6203 rpm, J = 0.4570, CT = 0.02209, CP = 0.01494, power 21.37 W, eta = 0.6759 (momentum ideal at same CT, J: 0.9404); residual 8.0e-13 in 41 it
  section Reynolds number at 0.75 R: 6.19e+04;  tip Mach 0.242
loss budget at the cruise point (same rpm, same V):
  ideal: no drag, no swirl, no tip loss    T = 1.324 N  P = 17.06 W  eta = 0.9314  (momentum ideal at this CT: 0.9352)
  + swirl                                  T = 1.292 N  P = 16.70 W  eta = 0.9286  (momentum ideal at this CT: 0.9366)
  + tip loss                               T = 1.231 N  P = 16.14 W  eta = 0.9156  (momentum ideal at this CT: 0.9392)
  + profile drag (full BEM)                T = 1.204 N  P = 21.37 W  eta = 0.6759  (momentum ideal at this CT: 0.9404)
uniform-inflow blade (n = 100.0 rev/s), ideal BEM: T = 1.20395 N (target 1.20395), u_a spread 4.80e-11 m/s, eta = 0.939223; actuator disk on the annulus: eta = 0.939223, v = 0.77652 m/s; BEM u_a = 0.77652 m/s
static at 6000 rpm: T = 3.266 N, P = 24.59 W, CT = 0.0641, CP = 0.0190, figure of merit 0.681, ideal power 16.75 W
sweep at 6203 rpm: eta_max = 0.6842 at J = 0.425; zero thrust at J = 0.5943 (geometric pitch/D = 0.4488)
```

and writes `propeller_bem.png`: C_T and C_P against J, η against J with the momentum-theory bound at the same C_T, and the radial distributions of thrust, axial inflow and F at cruise.

### Verify

Run `python -m pytest curriculum/phase-10/lesson-17-propeller-rotor-theory`. The tests check: momentum and energy conservation of the disk for four cases including static; the static power law (2^{3/2}, 2^{−1/2}); (5) against the disk and its light-loading limit; F between 0 and 1, zero at the tip and → 1 as B → ∞; R(φ) < 10⁻⁹ at the solution, the velocity triangle and the annulus momentum balance; static thrust ∝ n² with fixed coefficients; collapse of C_T, C_Q at fixed J and C_P = 2πC_Q; the uniform-inflow blade against the disk (thrust, inflow, efficiency, power); the ideal-annulus energy identity dP = dT(V + u_a); η ≤ η_i at three J with and without losses; each loss lowering η; C_T falling with J through zero to windmilling; the operating-point root; the fixed-point and bracketing solvers agreeing to 10⁻⁷ (including V = 0); the element study converging with order 1.5–2.5.

### Compare

The ideal blade reproduces momentum theory: thrust 1.20395 N on target, u_a = 0.77652 m/s on every ring (spread 5 × 10⁻¹¹ m/s) equal to the disk value on the annulus area, and η = 0.939223 in both. The real blade with the same switches off gives 0.9314 against 0.9352 for the disk at its thrust: its loading is not uniform, and (§2.4) any non-uniformity costs efficiency. With every loss on, BEM's 0.676 sits 0.264 below the bound.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No propeller is run (see Safety). Write each prediction down before running it.

1. **BEM against momentum theory, the loss budget.** At the cruise rpm and speed, switch on swirl, tip loss and profile drag one at a time. Predict which loss is largest. *Finding:* η = 0.931 (no losses but the induced one), 0.929 with swirl, 0.916 with tip loss, 0.676 with profile drag. The lossless blade sits 0.004 below the disk at its own thrust (non-uniform inflow) and the disk itself loses 0.060 to the induced velocity; swirl costs 0.003, tip loss 0.013 and profile drag 0.240. At J = 0.457 the sections work at α between −3.3° and −1.1° (cl ≈ 0.27 at 0.75R), where the assumed section's lift-to-drag ratio is only about 13: the blade is operating far below its best cl.

2. **Sweep in advance ratio.** At 6203 rpm, J = 0.05–0.60. Predict where η peaks and where thrust vanishes. *Finding:* C_T falls from 0.0591 (J = 0.1) to 0.0155 (J = 0.5) and through zero at J = 0.594; η rises to 0.684 at J = 0.425 and falls to 0.633 at J = 0.5. The momentum bound at the same C_T rises monotonically, from 0.399 to 0.963: the gap grows at high J because the blade's thrust shrinks while its profile drag does not. Zero thrust comes at J = 0.594, not at the geometric pitch ratio 0.449, because a cambered section still lifts at zero geometric angle: with α₀ = −4°, the aerodynamic pitch at 0.75R is 0.622 D.

3. **Static thrust.** At 6000 and 8000 rpm. Predict the scaling. *Finding:* 3.266 N and 5.807 N (ratio 1.778 = (8/6)²), 24.59 W and 58.29 W (ratio 2.37 = (8/6)³), C_T = 0.0641, C_P = 0.0190, figure of merit (ideal over actual power, `P10.18`) 0.681 at both speeds.

4. **Element and relaxation convergence.** Thrust error against n_el = 1280 for n_el = 10–160; the fixed-point iteration for n_el = 20–640 and relaxation 0.05–0.5. Predict the order and which elements make the iteration fail. *Finding:* errors 2.09 × 10⁻³, 6.51 × 10⁻⁴, 1.87 × 10⁻⁴, 5.04 × 10⁻⁵, 1.30 × 10⁻⁵ N, orders 1.68, 1.80, 1.89, 1.95, tending to 2 (midpoint rule on tip-clustered annuli). The fixed point converges at relaxation 0.3 up to n_el = 80, needs 0.2 at 160, 0.1 at 320, and 0.05 at 640; at the cruise point it never needs more than about 420 iterations when it converges. Its failures come from the elements nearest the tip, where F → 0 (F = 0.009 in the last of 160 elements) and the update (8) becomes very sensitive to φ. The bracketing solver converges for every resolution with no tuning (about 40 bisections).

5. **Section drag.** Repeat the cruise point with cd0 = 0.010, 0.020, 0.030. Predict the change in η. *Finding:* 6189, 6203, 6217 rpm and η = 0.779, 0.676, 0.597 (18.56, 21.37, 24.22 W). The rpm barely changes, the power by ±13 %: cd0 is the input to know best.

Raw data format: CSV tables `label, T_N, P_W, eta, eta_ideal` (1); `J, CT, CP, eta, eta_ideal` (2); `rpm, T_N, P_W, CT, CP, FM` (3); `n_el, err_N` and `n_el, relax, converged, iterations` (4); `cd0, rpm, eta, P_W` (5).

## 5. Validate and improve

| quantity | momentum theory | BEM | difference | cause |
|---|---|---|---|---|
| uniform-inflow blade, ideal limit: η | 0.939223 | 0.939223 | < 10⁻⁶ | none: the identity of §2.4 |
| uniform-inflow blade: u_a | 0.77652 m/s | 0.77652 m/s, spread 5 × 10⁻¹¹ | < 10⁻⁶ | none |
| reference blade, ideal limit: η at its C_T | 0.9352 | 0.9314 | −0.004 | non-uniform inflow, hub cut-out |
| reference blade, full BEM: η at cruise | 0.9404 | 0.6759 | −0.264 | profile drag 0.240, tip 0.013, swirl 0.003 |
| static: power for 3.266 N | 16.75 W | 24.59 W | FM 0.681 | profile drag, tip loss, non-uniform inflow |

Momentum theory checks BEM where they should agree, and bounds it where they should not. The engineering answer: **6203 rpm (J = 0.457) for the 1.20 N cruise thrust, with η = 0.676 and 21.4 W of shaft power; the ideal is 0.940 and 15.4 W.** What the comparison cannot validate is the section model: BEM against momentum theory tests the code, not the cd. The next step is measured data. For full-scale propellers, NACA TR-640 gives C_T, C_P and η charts (Extension 2); for small UAV propellers, wind-tunnel databases of measured C_T(J) and C_P(J) exist and are the right comparison for the reference propeller. Experiment 5 says where to spend the effort: cd0 changes the power by ±13 %; the element count and the solver change it by 10⁻⁵.

A design improvement, computed in Problem 4: a 10 × 6 in propeller (pitch 0.152 m) needs 5158 rpm at cruise, works at J = 0.550 near its own best efficiency and gives **η = 0.753, 19.18 W**, 10 % less power. Its price is static thrust: 2.76 N at its cruise rpm against 3.49 N for the 4.5 in pitch at its own cruise rpm, so its take-off and climb are weaker for the same motor.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(5), (10)–(11) and the energy identity of §2.4; `propeller_bem.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a pitch for the UAV's cruise and stating the uncertainty that matters most.

### Adversarial self-checks

- **What if the thrust doubles?** At 12 m/s, (3) gives a larger v and a lower ideal efficiency (heavier loading). BEM needs a higher rpm, which lowers J and raises the sections' cl towards their best lift-to-drag ratio, so the profile loss may shrink while the induced loss grows. `rpm_for_thrust(blade, 2.4, 12.0)` gives 7301 rpm (J = 0.388) and η = 0.678 against an ideal 0.893: the two effects almost cancel, and η stays below (5).
- **What if V → 0?** J → 0 and η → 0 by definition (no useful work), though the thrust is largest; the right figure of merit is then the ideal power over the actual power (0.681 here), as in `P10.18`.
- **Which assumption fails first in the real world?** The section data at Re ≈ 6 × 10⁴ (laminar separation bubbles, `P10.14`), then the axial-inflow assumption behind a fuselage or in a climb.
- **What measurement could falsify the model?** Measured thrust and torque on a test stand in a wind tunnel at several J (performed by a facility, not at home): if η(J) peaks at a different J or the zero-thrust J differs from 0.594, the blade geometry or α₀ is wrong; if C_P is off at all J, the cd level is.
- **Does the result survive a change of method or grid?** The two solvers agree to 10⁻¹² N where both converge; 60 elements are within 10⁻⁴ N of the converged thrust.
- **Which decision changes if the uncertainty doubles?** The pitch choice: if cd0 were uncertain by a factor of 2 rather than 1.5, the absolute η would be unknown between about 0.6 and 0.8, but both pitches move together, so the ranking (6 in better in cruise) survives; the motor sizing for 21 W would not.

## Check yourself

1. Define the advance ratio J and the coefficients C_T, C_P and η of a propeller. Which two are related by a factor 2π?
2. In the actuator-disk model, what fraction of the final slipstream velocity increase has occurred at the disk, and which conservation laws together give that result?
3. Why can no propeller beat the Froude efficiency at its own thrust, and why does uniform inflow achieve it?
4. What physical effect does Prandtl's tip-loss factor represent, and what are its values at the tip and for infinitely many blades?
5. Static ideal power scales as T^{3/2}/√(2ρA). By what factor does it change if the thrust doubles? If the diameter doubles?
6. Why does the zero-thrust advance ratio of the reference propeller (0.594) exceed its geometric pitch ratio (0.449)?
7. Why can the relaxed fixed-point iteration fail when the elements are crowded near the tip, while bisection of R(φ) cannot?
