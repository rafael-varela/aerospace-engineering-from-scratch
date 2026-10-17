# P08.06  Navier-Stokes Equations

> Newton's second law for a fluid particle, with a stress that is linear in the rate of strain. Most flows have no exact solution; the few that do are the test cases every viscous solver must pass.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the incompressible Navier–Stokes equations from the Cauchy momentum equation and the Newtonian stress law, and state the no-slip and no-penetration conditions. 2. Reduce them by stated assumptions to the exact Couette, plane Poiseuille, Couette–Poiseuille and Hagen–Poiseuille solutions and to Stokes's first problem, and derive flow rates, wall stresses and diffusion times. 3. Implement and test the closed forms, a SymPy substitution check of the full equations, a tridiagonal finite-difference solver for steady channel and pipe flow and a Crank–Nicolson solver for the unsteady problems, with observed order of convergence. 4. Decide whether the exact laminar solution describes a real narrow-gap flow (an actuator clearance, a bearing film) by checking its Reynolds number, its diffusion time and the assumptions it makes. |
| **Prerequisites** | `P08.03` (mass and momentum balances for a control volume), `P08.02` (material derivative, strain rate and rotation), `P08.05` (Reynolds number), `P02.11` (the diffusion equation and similarity solutions), `P01.08` (partial derivatives), `P03.11` (finite differences, truncation error, observed order), `P03.08` (tridiagonal systems). Skill check: you can integrate u″ = constant with two boundary conditions, and you can write the central difference for u″. |
| **Engineering question** | A flight-control hydraulic actuator's 25 mm piston slides at up to 0.1 m/s in its bore with a 10 µm radial clearance across a 12 mm land holding 20.7 MPa: how much fluid leaks past it, what viscous force does the piston feel, and is a steady laminar exact solution the right model? The answers: 0.057 L/min with the piston at rest (0.060 L/min moving with the leak), a shear force of 8.0–8.1 N (0.08 % of the 10.2 kN actuator force), a gap Reynolds number below 1 and a start-up time of 3.4 µs: yes, as long as the clearance (cubed) and the piston's eccentricity are known. |
| **Assumptions and limits** | Incompressible Newtonian fluid with constant density and viscosity; laminar flow; walls rigid and smooth; flows fully developed (no change along the flow); the narrow annular gap treated as a plane slot (error of order (h/R)², checked against the exact annulus); concentric piston. Fails for turbulent or transitional flow, for large viscous heating (viscosity falls as the fluid warms), for an eccentric piston (up to 2.5 times the leakage) and for non-Newtonian fluids. Hydraulic-fluid properties are typical values, not data-sheet values. |
| **Analytical work** | Cauchy's equation, the Newtonian stress and the reduction to ρDu/Dt = −∇p + μ∇²u + ρg; the exact parallel-flow solutions and their flow rates and wall stresses; the momentum balance τ(h) − τ(0) = −Gh; the Hagen–Poiseuille flow and f = 64/Re; the similarity solution erfc(y/(2√(νt))); the start-up series of plane Poiseuille flow and its 99 % time; the truncation error of the central difference and why it vanishes for a parabola. |
| **Computational deliverable** | `starter/main.py`: `couette_poiseuille`, `plane_flow_rate`, `wall_shear_stress`, `pipe_poiseuille`, `pipe_flow_rate`, `annulus_flow_rate`, `stokes_first_problem`, `startup_channel_series`, `startup_time`, `ns_residuals`, `pipe_axial_residual`, `fd_steady_channel`, `fd_steady_pipe`, `crank_nicolson` (`_tridiag_solve`, `observed_order`, `actuator_leakage` and `main` are given). Tests prove: each exact solution leaves a zero SymPy residual in the full Navier–Stokes equations (Couette–Poiseuille, Taylor–Green including its nonlinear terms, hydrostatics, Stokes's first problem, pipe flow) while a wrong profile does not; flow rates equal the integrals of the profiles; the wall stresses close the slab's momentum balance; superposition and zero limits; the annulus formula matches its own profile and tends to the slot; the steady finite-difference solvers reproduce the parabolas to round-off and converge at order 2.00 ± 0.05 on a manufactured problem; Crank–Nicolson converges at order 2 on Stokes's first problem and on the start-up flow, and relaxes to Poiseuille; the 99 % start-up time; the actuator's h³ scaling and force balance. |
| **Practical activity** | **computational experiment**: a convergence study of the steady and unsteady finite-difference solvers against the exact solutions (including a manufactured solution that is not a polynomial), a sweep of clearance and piston speed for the actuator, and a check of the plane-slot approximation against the exact annulus. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `navier_stokes.png`, the derivations of the Couette–Poiseuille profile, the pipe flow and Stokes's similarity solution, the four experiment tables, and a half-page note on the actuator's leakage and friction with the assumption that most limits the answer. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the Navier–Stokes equations are derived with every assumption named; the exact solutions are reduced from them, not quoted; all public tests pass; observed orders are reported for the manufactured and the unsteady problems with the reason the parabolas are reproduced exactly; the note states the leakage, the gap Reynolds number and the start-up time, and ranks the assumptions. |
| **References** | `@mit-ocw-2-25-advanced-fluid-mechanics` topic "Equations of Viscous Flow", class notes "Equation of Motion for Viscous Flow", "Couette & Poiseuille Flows" and "Criteria for Locally Fully Developed Viscous Flow" (https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/equations-of-viscous-flow/) (required-free: the equations, the parallel-flow solutions, flow rates and wall stresses); topic "More Complex Viscous-Dominated Flows", "Flow Inside a Cylinder Which is Suddenly Rotated" (supplementary-free: an unsteady viscous diffusion problem). `@leishman-intro-aerospace-flight-vehicles` chapter "Momentum Equations": Momentum Equation – Differential Form; Navier-Stokes Equations, and chapter "Viscous-Dominated Flows": Stokes Flow Between Parallel Plates; Stokes' First and Second Problems; Taylor-Couette Flow (required-free; link only). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 8 "Introduction to Viscous Flows": Viscous Stress; Navier-Stokes Equation; Examples of Exact Solutions to Navier-Stokes Equation (https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/08%3A_Introduction_to_Viscous_Flows) (required-free). `@libretexts-bar-meir-fluid-mechanics` Ch. 8 "Differential Analysis": 8.5 Derivations of the Momentum Equation; 8.6 Boundary Conditions and Driving Forces; 8.7 Examples for Differential Equation (Navier-Stokes) (supplementary-free). `@libretexts-smyth-all-things-flow` Ch. 6 "Fluid Dynamics": 6.3 Momentum conservation; 6.8 Summary: the equations of motion; 6.9 Boundary conditions (supplementary-free). Optional textbooks for deeper study: F. M. White, *Viscous Fluid Flow* (paid; a library copy is enough), the standard collection of exact solutions. |
| **Extensions** | ★ Add Stokes's second problem (an oscillating wall) to the tests and measure the start-up transient (Problem 5). ★ Solve the eccentric annulus numerically on a polar grid and reproduce the factor 1 + 1.5ε² for the leakage. ★ Add a temperature equation with viscous dissipation and a viscosity that falls with temperature, and find how much the 20.7 MPa drop raises the leakage. |
| **Time** | unmeasured |

## 1. The engineering problem

A flight-control actuator is a piston in a cylinder, pushed by hydraulic fluid at 20.7 MPa (3000 psi, the pressure most aircraft hydraulic systems use; aviation practice states it in psi). The piston must slide freely, so it cannot be a perfect fit: there is a radial clearance of a few micrometres between piston and bore, and fluid leaks through it from the high-pressure side to the low-pressure side. Leakage wastes pump power and heats the fluid, and in a servo valve it sets how stiffly the actuator holds its position. The same narrow gap also shears the fluid as the piston moves, which is a friction force on the piston.

The piston here is 25 mm in diameter, with a 12 mm land (the sealing length) and a radial clearance h = 10 µm, and it moves at up to 0.1 m/s. The fluid is a synthetic hydrocarbon hydraulic fluid of the MIL-PRF-83282 type, with typical properties near 40 °C of ρ = 845 kg/m³ and ν = 14 mm²/s (μ = 0.0118 Pa·s). The questions: how much fluid leaks, what force acts on the piston, and can a steady laminar flow with a closed-form solution be trusted to answer?

The tool is the Navier–Stokes equations, the momentum balance of `P08.03` written for a fluid particle with a model for the stress. This lesson derives them, reduces them to the handful of flows that have exact solutions, and checks each solution twice: by substituting it back into the full equations with SymPy, and by solving the reduced equations numerically with finite differences. Success criteria, written first: every exact solution leaves a zero residual; the finite-difference solutions converge at second order (or are exact, with the reason stated); the flow rates and wall stresses satisfy the slab's momentum balance to round-off; the actuator's leakage, force, Reynolds number and start-up time are computed with the assumptions that limit them.

## 2. Understand and derive

### 2.1 Newton's law for a fluid particle

Follow a small particle of fluid of density ρ. Its acceleration is the material derivative of `P08.02`, Du/Dt = ∂u/∂t + (u·∇)u: the change at a fixed point plus the change from moving to where the velocity is different. The forces on it are body forces ρg per unit volume and surface forces from its neighbours. The surface force on a face with unit normal n is σ·n, with σ the **stress tensor** (nine components, symmetric by the balance of moments, so six independent). Summing the surface forces over the faces of the particle and dividing by its volume gives the divergence of σ (the divergence theorem of `P01.12` applied to a vanishing volume). So

**ρ Du/Dt = ∇·σ + ρg**   (Cauchy's equation)   (1)

This is exact for any continuum. It has more unknowns than equations until we say how σ depends on the motion. Together with mass conservation, which for an incompressible fluid is **∇·u = 0** (`P08.03` in differential form), it is the momentum equation in differential form (Bar-Meir 8.5, "Derivations of the Momentum Equation"; Smyth 6.3).

### 2.2 The Newtonian fluid

A fluid at rest carries only pressure: σ = −pI. In motion, the extra (viscous) stress must depend on how the particle deforms, not on how it translates or rotates rigidly, because a rigid motion does not deform anything. `P08.02` split the velocity gradient into the rate-of-strain tensor S = ½(∇u + ∇uᵀ) and the rotation. A **Newtonian fluid** is one whose viscous stress is linear in S and the same in every direction:

**σ = −pI + 2μS + λ(∇·u)I.**   (2)

μ is the dynamic viscosity (Pa·s) and λ a second coefficient that multiplies the volume change rate ∇·u and so vanishes for incompressible flow. Air, water, fuels and hydraulic fluids are Newtonian to a very good approximation; blood, paint and polymer solutions are not. In 2-D the shear component is τ_xy = μ(∂u/∂y + ∂v/∂x), which for a flow u(y) is μ du/dy: Newton's law of viscosity from `P08.01`.

Substitute (2) into (1) with ∇·u = 0 and constant μ. ∇·(2μS) = μ(∇²u + ∇(∇·u)) = μ∇²u, so

**ρ(∂u/∂t + u·∇u) = −∇p + μ∇²u + ρg,   ∇·u = 0.**   (3)

These are the **incompressible Navier–Stokes equations**: four equations (three momentum components and continuity) for four unknowns (u, v, w, p). Check the units: every term is N/m³; μ∇²u is Pa·s × (m/s)/m² = Pa/m ✓. The left side is inertia, which is nonlinear in u; the right side is pressure, viscous and body forces. Their ratio, inertia over viscous, is the Reynolds number ρUL/μ of `P08.05`.

**Boundary conditions.** At a solid wall the fluid does not pass through it (no penetration) and does not slip along it: **u = u_wall** (no slip). No slip is an observation, not a theorem; it holds for liquids and for gases except in very rarefied flow. At a free surface the shear stress is zero and the pressure is the outside pressure (Problem 3).

The nonlinear term u·∇u is why (3) has so few exact solutions: in general it couples every component to every other. It vanishes identically in one important family, **parallel flows**, where u = (u(y, t), 0, 0): then u·∇u = u ∂u/∂x = 0. Every exact solution in this lesson is of that kind, plus one (the Taylor–Green vortex, used only as a test) where the nonlinear term is balanced exactly by the pressure gradient.

### 2.3 Couette–Poiseuille flow

Take two parallel plates a distance h apart, the lower one at rest and the upper moving at U along x, with a pressure gradient dp/dx = −G between them (G > 0 when pressure falls in +x). Assume the flow is **steady** and **fully developed**: nothing changes along x except the pressure (MIT 2.25, "Criteria for Locally Fully Developed Viscous Flow", gives the conditions: a gap much smaller than the length, and a small Reynolds number based on the gap times h/L). Then:

- Continuity: ∂u/∂x + ∂v/∂y = 0 with ∂u/∂x = 0 gives ∂v/∂y = 0; v = 0 on both walls, so **v = 0 everywhere**.
- y-momentum: all velocity terms vanish, so ∂p/∂y = ρg_y. Pressure across the gap is hydrostatic and does not drive the flow.
- x-momentum: ρ(u u_x + v u_y) = 0, so **μ d²u/dy² = dp/dx = −G**. The left side depends only on y and the right only on x, so both are constant: the pressure falls linearly along the gap.

Integrating twice with u(0) = 0 and u(h) = U:

**u(y) = U y/h + G y (h − y)/(2μ).**   (4)

The first term is **Couette flow** (a linear profile dragged by the wall) and the second is **plane Poiseuille flow** (a parabola pushed by the pressure gradient). The equation is linear, so the two simply add. The flow rate per unit width and the shear stresses τ = μ du/dy at the walls are

**q = Uh/2 + Gh³/(12μ),   τ(0) = μU/h + Gh/2,   τ(h) = μU/h − Gh/2.**   (5)

The two stresses differ by **τ(h) − τ(0) = −Gh**, and this is a control-volume check (`P08.03`): a slab of fluid of length dx feels a net pressure force Gh dx forward, and the walls hold it back with τ(0)dx on its lower face and −τ(h)dx on its upper face. The test `test_wall_shear_closes_the_momentum_balance` checks it to round-off. The force per unit area the fluid exerts on the lower wall is +τ(0) in x and on the upper wall −τ(h).

```
   y = h  ═══════════════════════  wall moving at U  →
          ─────────→                Couette: linear
          ──────────────→           Poiseuille: parabola
          ─────────→                (sum: equation 4)
   y = 0  ═══════════════════════  wall at rest
          p high  ─────────────→  p low      G = −dp/dx > 0
```

**The actuator.** The clearance is an annulus of mean circumference πD and height h = 10 µm, with h/R = 8 × 10⁻⁴: unrolled, it is a plane slot of width πD. With G = Δp/L = 20.7 MPa/12 mm = 1.725 × 10⁹ Pa/m and the piston at rest, (5) gives Q = πD Gh³/(12μ) = **9.544 × 10⁻⁷ m³/s = 0.0573 L/min**. With the piston moving at 0.1 m/s in the direction of the leak, the Couette part adds πD Uh/2 = 3.93 × 10⁻⁸ m³/s, for **0.0596 L/min**. The mean gap speed is 1.27 m/s and the gap Reynolds number ρ V h/μ = **0.90**: the flow is laminar by a factor of a thousand.

The shear force on the piston is −τ(h) × πDL. The pressure-driven part, (Gh/2)πDL = **8.129 N**, is exactly half of the pressure force on the fluid in the gap, Δp πDh = 16.258 N; the bore carries the other half. The Couette part, −(μU/h)πDL = −0.112 N at 0.1 m/s, is small, so the total is **8.017 N**: 0.08 % of the 10.16 kN that 20.7 MPa exerts on the piston face. Leakage, not friction, is the engineering issue. And it goes as **h³**: a clearance of 5, 10 and 20 µm leaks 0.0072, 0.0573 and 0.458 L/min.

### 2.4 Pipe flow and the annulus

In a round pipe of radius R the same assumptions (steady, fully developed, w = w(r) only) reduce the axial momentum equation, in cylindrical coordinates, to **μ (1/r) d/dr(r dw/dr) = dp/dz = −G**. Integrating once, r w′ = −Gr²/(2μ) + C₁; a finite shear on the axis needs C₁ = 0. Again with w(R) = 0:

**w(r) = G(R² − r²)/(4μ),   Q = πR⁴G/(8μ),   τ_w = GR/2.**   (6)

This is **Hagen–Poiseuille flow**: the flow rate goes as R⁴, so halving a line's radius at the same pressure gradient cuts its flow by sixteen. With the mean speed V = Q/(πR²) = GR²/(8μ) and the Darcy friction factor f = 8τ_w/(ρV²), (6) gives **f = 64/Re**, Re = ρVD/μ (Problem 2): the laminar friction factor that `P08.07` uses for feed lines.

The concentric annulus Ri < r < Ro has the exact solution w = (G/4μ)[Ro² − r² + (Ro² − Ri²) ln(r/Ro)/ln(Ro/Ri)] and

Q = (πG/8μ)[Ro⁴ − Ri⁴ − (Ro² − Ri²)²/ln(Ro/Ri)].   (7)

For a thin gap this is a small difference of large terms; the code factors it to avoid losing digits (a first version typed directly from (7) put the annulus-to-slot ratio 3.5 × 10⁻⁷ away from 1 at h/R = 8 × 10⁻⁴, where the true difference is about 10⁻⁸). The ratio of (7) to the slot formula π(Ri + Ro)h³G/(12μ) is 1.000001 for a 0.1 mm gap and 1.000099 for a 1 mm gap on the 25 mm piston: the error of unrolling the annulus is of order (h/R)², below 10⁻⁷ for the actuator. Eccentricity is a far bigger effect: a fully eccentric piston leaks 2.5 times as much (Problem 6).

### 2.5 Unsteady flows: diffusion of momentum

Drop the steady assumption but keep parallel flow. The x-momentum equation becomes

**∂u/∂t = ν ∂²u/∂y² + G/ρ,**   (8)

the diffusion equation of `P02.11` with ν = μ/ρ as the diffusivity of momentum, plus a source. Its solutions answer "how long until the flow is steady?".

**Stokes's first problem.** A plate at y = 0 is started impulsively at speed U under fluid at rest; G = 0. The problem has no length scale (only U, ν and t), so the only length available is √(νt), and the solution must have the form u/U = F(η) with η = y/(2√(νt)). Substituting into (8) gives F″ + 2ηF′ = 0 with F(0) = 1 and F(∞) = 0:

**u = U erfc(y/(2√(νt))).**   (9)

The layer of moving fluid is 3.64√(νt) thick (where u = 0.01U) and grows as √t (Leishman, "Viscous-Dominated Flows", Stokes' First and Second Problems; Liburdy Ch. 8). In the hydraulic fluid, after 1 ms it is 431 µm thick: 43 clearances.

**Start-up of Poiseuille flow.** Fluid at rest between fixed plates; at t = 0 the pressure gradient G is switched on. Write u = u_s(y) − w(y, t) with u_s the steady parabola. Then w satisfies the plain diffusion equation with w = 0 at both walls and w(y, 0) = u_s(y). Separating variables, each Fourier sine mode sin(nπy/h) decays as exp(−n²π²νt/h²), and the sine series of the parabola has coefficients 8h²G/(2μ(nπ)³) for odd n:

**u = u_s(y) − (G/2μ) Σ_{n odd} (8h²/(nπ)³) sin(nπy/h) exp(−n²π²νt/h²).**   (10)

The slowest mode, n = 1, sets the time to steady state. At mid-gap the series becomes u/u_s = 1 − (32/π³)e^{−π²νt/h²} + …, so the flow reaches 99 % of its steady centreline speed at

**t₉₉ = (h²/(π²ν)) ln(32/(π³ × 0.01)) = 0.470 h²/ν.**   (11)

For the actuator, h²/ν = 7.14 µs and **t₉₉ = 3.36 µs**. A piston that changes speed over milliseconds sees a gap flow that is steady at every instant: quasi-steady. The test checks (11) against the full series.

### 2.6 Where the exact solutions stop

- **Turbulence.** (4)–(6) are exact for any Reynolds number, but they are stable only at low ones. Pipe flow typically becomes turbulent near Re ≈ 2300 (`P08.07`); a Taylor–Couette film between a rotating shaft and its bearing becomes unstable to Taylor vortices when the Taylor number Re²(h/R) exceeds about 1700 (the flow is described in Leishman, "Viscous-Dominated Flows", Taylor-Couette Flow; Problem 1 applies the criterion). The actuator, at Re = 0.9, is far from either.
- **Entrance length.** Fully developed flow takes about 0.06 Re h to establish in a laminar slot: 0.5 µm here, against a 12 mm land.
- **Viscous heating.** The leaking fluid drops 20.7 MPa doing no work; if the energy stays in it, it warms by Δp/(ρc_p) ≈ 13 K (Problem 6). Viscosity falls as fluids warm, so the constant-μ answer underestimates the leakage.
- **Geometry.** Eccentricity (up to ×2.5), taper of the bore, and the bore's elastic expansion under 20 MPa all change h, and h enters cubed.
- **The fluid.** Typical properties, not a data sheet; the fluid's viscosity at the actual operating temperature is the input that matters.

## 3. Build the computational model

Open `starter/main.py`. `_tridiag_solve` (a wrapper of `scipy.linalg.solve_banded`), `observed_order`, `actuator_leakage` and `main` are given.

### Step 1: the closed forms

`couette_poiseuille`, `plane_flow_rate`, `wall_shear_stress`, `pipe_poiseuille` and `pipe_flow_rate` are (4)–(6). `annulus_flow_rate` is (7) written with δ = (Ro − Ri)/Ro and k = Ri/Ro as δ(1 + k)[1 + k² − δ(1 + k)/ln(1/k)], with ln(1/k) = −`log1p`(−δ), so the thin-gap cancellation happens inside a factor rather than between terms of size Ro⁴. `stokes_first_problem` is (9) with `scipy.special.erfc`; `startup_channel_series` is (10) summed over 200 odd terms; `startup_time` is (11) for any fraction.

### Step 2: substitution with SymPy

`ns_residuals(u, v, p, rho, mu, gx, gy)` takes SymPy expressions in the module's symbols X, Y, T and returns the three residuals of (3) in 2-D (continuity, x- and y-momentum), simplified. Every term is there, including the nonlinear ones, so the function is a test of any claimed solution:

```python
d = sp.diff
cont = d(u, X) + d(v, Y)
xm = rho*(d(u, T) + u*d(u, X) + v*d(u, Y)) + d(p, X) - mu*(d(u, X, 2) + d(u, Y, 2)) - rho*gx
```

`pipe_axial_residual(w, p, mu)` does the same for the axisymmetric axial equation in R_SYM and Z_SYM.

### Step 3: steady finite differences

`fd_steady_channel(n, h, U, f, mu)` solves μu″ = −f(y) on n intervals with u(0) = 0, u(h) = U, the central difference (u_{i−1} − 2u_i + u_{i+1})/Δy² at the interior nodes and the boundary values moved to the right-hand side: a tridiagonal system of n − 1 unknowns (`P03.08`). `f` is a constant G or any function f(y), so the same solver handles a manufactured problem. `fd_steady_pipe(n, R, G, mu)` uses face radii, [r_{i+½}(w_{i+1} − w_i) − r_{i−½}(w_i − w_{i−1})]/(r_i Δr²), and on the axis the limit (1/r)(rw′)′ → 2w″ with the mirror node w₋₁ = w₁, i.e. 4(w₁ − w₀)/Δr².

The truncation error of the central difference is (Δy²/12)u'''' + O(Δy⁴). For the parabolas (4) and (6), u'''' = 0, so **the discrete solution is exact at the nodes on any grid**: the test checks n = 2, 5 and 32 to 10⁻¹² (the pipe's face-radius stencil is exact for a parabola too). That is a correct result, but a useless convergence test. So the second steady test uses a **manufactured solution**, u = sin(πy/h), with the forcing f = μ(π/h)² sin(πy/h) that makes it exact, and measures the observed order.

### Step 4: Crank–Nicolson for the unsteady problems

`crank_nicolson(y, u0, t0, t1, nt, nu, source, u_lo, u_hi)` marches (8) by averaging the central difference at the old and new time levels:

(u_i^{k+1} − u_i^k)/Δt = (ν/2)(δ²u^{k+1} + δ²u^k)_i + s,

a tridiagonal solve per step, second order in Δt and Δy and unconditionally stable. The wall values may depend on time (callables). For Stokes's first problem the test starts from the exact profile at t₀ = 0.01 (in units with ν = 1) rather than from the impulsive start, whose discontinuity at t = 0 would spoil the order; the far boundary at y = 2 takes the exact value. The start-up problem starts from rest with source G/ρ; its initial data are continuous, and second order survives.

Running `main()` prints

```
hydraulic fluid (typical): rho = 845 kg/m^3, mu = 0.01183 Pa s, nu = 14.0 mm^2/s
actuator: D = 25 mm, clearance h = 10 um, land 12 mm, dp = 20.7 MPa (G = 1.725e+09 Pa/m), piston speed 0.1 m/s
  leakage, piston at rest:   Q = 9.5436e-07 m^3/s = 0.0573 L/min
  leakage, piston moving:    Q = 9.9363e-07 m^3/s = 0.0596 L/min (Couette part 3.927e-08 m^3/s)
  mean gap speed 1.265 m/s, gap Reynolds number 0.904
  shear force on piston 8.017 N = pressure-driven 8.129 N + wall-driven -0.1115 N
  pressure force on the land annulus dp * pi D h = 16.258 N
  start-up to 99 %: 3.356 us (h^2/nu = 7.143 us)
  concentric annulus Q = 9.54743e-07 m^3/s, plane slot pi (Ri+Ro) h^3 G/(12 mu) = 9.54743e-07: |ratio - 1| < 1e-7 is True
  annulus / slot for a 0.1 mm gap on a 25 mm piston: 1.000001
  annulus / slot for a 1.0 mm gap on a 25 mm piston: 1.000099
  clearance    5 um: leakage 0.0072 L/min
  clearance   10 um: leakage 0.0573 L/min
  clearance   20 um: leakage 0.4581 L/min

finite differences, steady Couette-Poiseuille (max |u_FD - u_exact| / max u):
  n =   4: 1.19e-16  (exact: u'''' = 0)
  n =   8: 2.37e-16  (exact: u'''' = 0)
  n =  16: 2.37e-16  (exact: u'''' = 0)
  manufactured u = sin(pi y): errors [8.265e-03 2.059e-03 5.142e-04 1.285e-04], orders [2.005 2.001 2.   ]
  Stokes's first problem, Crank-Nicolson (dt ~ dy): errors [3.441e-03 8.621e-04 2.146e-04 5.372e-05], orders [1.997 2.006 1.998]
  start-up Poiseuille from rest: errors [7.401e-05 1.857e-05 4.645e-06 1.162e-06], orders [1.995 1.999 2.   ]
```

and writes `navier_stokes.png`: the clearance profiles for the piston at rest and moving both ways with finite-difference points; Stokes's first problem at three times; the log–log convergence of the three finite-difference studies with a slope-2 line.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-06-navier-stokes-equations`. The tests:

1. **Substitution**: the Couette–Poiseuille profile, the Taylor–Green vortex (an unsteady 2-D solution with all nonlinear terms active), the hydrostatic state and Stokes's first problem leave zero residuals; a cubic profile and Taylor–Green without its pressure do not; the pipe profile satisfies the axial equation and a doubled pressure gradient does not.
2. **Boundary conditions and integrals**: no slip on both walls; erfc's wall value, decay and self-similarity, and its 99 % thickness 3.64√(νt); the pipe and channel flow rates equal quadratures of the profiles; the pipe's mean speed is half its centreline speed.
3. **Momentum balance**: τ(h) − τ(0) = −Gh to 10⁻¹²; the stresses equal μ du/dy by differencing the profile; pure Couette flow has one stress.
4. **Superposition and limits**: Couette plus Poiseuille equals Couette–Poiseuille; no driving gives no flow.
5. **Annulus**: (7) equals the quadrature of its own profile and tends to the slot with error below (h/R)².
6. **Finite differences**: the parabolas to round-off on coarse grids; order 2.00 ± 0.05 on the manufactured problem; the pipe's flow rate by the trapezoid rule; Crank–Nicolson at order 2 ± 0.1 on Stokes's first problem and above 1.8 on the start-up flow, which also relaxes to the parabola; the 99 % time against the series.
7. **Actuator**: Q = πDh³Δp/(12μL), the h³ law, Re < 10, the two walls sharing the pressure force, the force decomposition, t₉₉ below 0.1 ms.
8. **main** writes the figure.

### Compare

The established tool for the unsteady problem is the method of lines: discretise in space with the same central difference and hand the ODE system to `scipy.integrate.solve_ivp` with its BDF integrator at a tight tolerance. On Stokes's first problem (80 intervals, t from 0.01 to 0.05) that gives a maximum error of 2.29 × 10⁻⁴ against (9); `crank_nicolson` with 40 time steps gives 2.15 × 10⁻⁴. With 4000 steps it agrees with the BDF solution to 1.6 × 10⁻⁹, and with 40 steps to 1.5 × 10⁻⁵. The two codes solve the same semi-discrete system, so they converge to each other as Δt → 0, and the remaining error against the exact solution is the spatial one: on this grid, Δy² dominates and the time step chosen (Δt ∝ Δy) is about right. SymPy plays the role of the established tool for the substitution checks; it is independent of anything written in this lesson.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Steady convergence.** `fd_steady_channel` on the actuator gap for n = 4, 8, 16, and on the manufactured problem for n = 10–80. Predict the error on the parabola and the order on the sine. *Finding:* 1–2 × 10⁻¹⁶ relative on every grid (round-off); manufactured errors 8.265 × 10⁻³, 2.059 × 10⁻³, 5.142 × 10⁻⁴, 1.285 × 10⁻⁴, orders 2.005, 2.001, 2.000. A convergence study on the exact parabola would have "shown" infinite order; the manufactured problem is the honest test.
2. **Unsteady convergence.** Crank–Nicolson with Δt ∝ Δy for Stokes's first problem (n = 20–160) and the start-up flow. Predict both orders. *Finding:* the table above, orders 2.00 for both. The Stokes run starts from the exact profile at t₀ rather than from the impulsive start, whose jump at the wall is not smooth initial data and would test the start-up, not the scheme.
3. **Clearance and piston speed.** Leakage and force for h = 5, 10, 20 µm at rest, and for U = −1 to +1 m/s at h = 10 µm. Predict the leakage ratio for doubling h and the speed at which the leak would stop. *Finding:* 0.0072, 0.0573, 0.4581 L/min (×8 per doubling); at U = −1.0, −0.3, −0.1, 0, 0.1, 0.3, 1.0 m/s the leakage is 0.0337, 0.0502, 0.0549, 0.0573, 0.0596, 0.0643, 0.0808 L/min and the force 9.244, 8.463, 8.240, 8.129, 8.017, 7.794, 7.014 N; the leak would stop only at U = −Gh²/(6μ) = −2.43 m/s, an unreachable piston speed. The piston's motion modulates the leak by tens of per cent; the clearance changes it by factors of eight.
4. **Slot against annulus.** (7) against the slot for gaps of 10 µm, 0.1 mm and 1 mm on a 25 mm piston. Predict how the ratio scales. *Finding:* within 10⁻⁷, 1.000001, 1.000099: the ratio departs from 1 as (h/R)², so the plane-slot model is not where the uncertainty lies.

Raw data format: CSV tables `n, max_err_parabola, max_err_manufactured` (1); `n, err_stokes, err_startup` (2); `h_um, U_m_s, Q_L_min, force_N, Re_gap` (3); `gap_mm, ratio_annulus_slot` (4); plus `navier_stokes.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| Couette–Poiseuille residual | 0 | SymPy: (0, 0, 0) | 0 | exact solution |
| channel profile, FD n = 4 | (4) | FD | 1.2 × 10⁻¹⁶ (relative) | u'''' = 0: no truncation error |
| manufactured sine, FD | order 2 | 2.005, 2.001, 2.000 | — | (Δy²/12)u'''' |
| Stokes's first problem, CN | order 2 | 1.997, 2.006, 1.998 | — | Δt² + Δy² |
| start-up flow, CN | (10), order 2 | 1.995, 1.999, 2.000 | — | — |
| slab momentum balance | τ(h) − τ(0) = −Gh | to 10⁻¹² | round-off | — |
| annulus vs slot (10 µm) | 1 + O((h/R)²) | < 10⁻⁷ | — | thin-gap limit |
| t₉₉ | 0.470h²/ν = 3.356 µs | CN across the gap: 0.98998 of steady at t₉₉ (Problem 4, 7 µm) | 2 × 10⁻⁵ | one-term formula |

Every computed quantity agrees with its independent check to the accuracy the method predicts. The engineering answer: **the 10 µm actuator clearance leaks 0.057 L/min at rest and 0.060 L/min with the piston moving at 0.1 m/s; the piston feels about 8 N of viscous force, negligible against its 10.2 kN; the gap flow is laminar (Re ≈ 0.9) and steady within 3.4 µs, so the exact solution is the right model.** What the exact solution cannot give is the clearance itself: ±1 µm changes the leakage by −27 %/+33 %, eccentricity up to ×2.5, and heating of the leaking fluid lowers μ. The arithmetic is exact; the geometry and the temperature are not.

What to improve: measure the leakage at operating temperature and infer an effective clearance from (5); then add the eccentricity and the thermal effect (Extensions 2 and 3). What not to improve: the grid, or the plane-slot approximation.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4)–(5), (6) and (9); `navier_stokes.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the actuator's leakage and friction force, the gap Reynolds number and start-up time that justify the model, and the assumptions ranked by their effect.

### Adversarial self-checks

- **What if the pressure doubles?** Leakage and the pressure-driven force double (the problem is linear in G); Re doubles to about 1.8, still laminar. The heating doubles too, to about 26 K, and that is where the linear answer starts to fail.
- **What if the viscosity goes to zero?** Q ∝ 1/μ grows without bound, Re grows as 1/μ², and the laminar exact solution is eventually unstable: a gap leaking water instead of hydraulic fluid (μ smaller by a factor of 12) would need its Reynolds number checked before (4) is trusted.
- **Which assumption fails first in the real world?** The concentric, uniform clearance: pistons are pushed to one side by side loads, and h³ makes small changes large. Then the constant viscosity.
- **What measurement could falsify the model?** Leakage measured against Δp at constant temperature: (5) predicts a straight line through the origin with slope πDh³/(12μL). Curvature means heating or a deforming bore; a slope off by a factor of 2 means eccentricity or a different clearance.
- **Does the result survive a change of method?** The closed form, SymPy, the steady finite-difference solver and Crank–Nicolson relaxing to steady state all give the same profile; the method of lines with `solve_ivp` agrees with Crank–Nicolson to 10⁻⁹ when the time step is small.
- **Which decision changes if the uncertainty doubles?** If the clearance is known to ±2 µm instead of ±1 µm on 10 µm, the leakage spans 0.029 to 0.099 L/min, a factor of 3.4: a pump sized for the nominal leak would need a margin that large, or a measured clearance.

## Check yourself

1. Write the incompressible Navier–Stokes momentum equation and name each term and its units.
2. What does "Newtonian fluid" assume about the stress, and what does incompressibility do to the viscous term?
3. Why does the central-difference solver reproduce plane Poiseuille flow exactly on any grid?
4. How do the leakage through a slot and the flow through a pipe scale with the gap h and the radius R at a fixed pressure gradient?
5. What is the time for viscous diffusion across a gap of height h, and what does it tell you about whether a gap flow is quasi-steady?
6. In Couette–Poiseuille flow the two wall shear stresses differ by Gh. Where does that difference come from?
7. Why is Stokes's first-problem profile self-similar, and how thick is the layer that moves?
8. Name two assumptions of the actuator calculation that a real actuator violates, and say which way each moves the leakage.
