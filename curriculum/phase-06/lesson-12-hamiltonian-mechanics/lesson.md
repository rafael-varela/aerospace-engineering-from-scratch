# P06.12  Hamiltonian Mechanics ★

> Trade velocities for momenta and every motion becomes a flow in phase space that carries volume like an incompressible fluid: shapes stretch, volumes never change.

This is an optional (★) lesson. Nothing later in the phase depends on it, but lessons on orbit determination and uncertainty propagation in later phases use its language.

| | |
|---|---|
| **Learning outcomes** | 1. Derive a Hamiltonian from a Lagrangian by the Legendre transform, write Hamilton's equations, and use Poisson brackets to show which quantities are conserved (energy, a cyclic momentum, the Laplace–Runge–Lenz vector). 2. Draw and read the pendulum's phase portrait (libration, rotation, separatrix), compute the exact period with the complete elliptic integral, and check it independently as the energy derivative of the enclosed phase area, T = dA/dE. 3. State Liouville's theorem, demonstrate it by carrying the boundary of a blob of states along the exact flow, and show that linear damping shrinks the area exactly as e^{−ct}. 4. Compute the area factor (Jacobian determinant) of the explicit Euler, RK4 and Störmer–Verlet steps and predict what each does to a blob and to the energy. 5. Treat the Kepler problem in phase space: effective potential, radial action and period, the symplectic flow Jacobian, and the along-track growth of an orbit-determination error. |
| **Prerequisites** | `P06.11` (Lagrange's equations), `P03.07` (symplectic integrators, the modified energy, the Kepler oracle), `P02.04` (phase portraits), `P04.05` (Kepler's laws), `P04.06` (the pendulum's amplitude-dependent period), `P01.10` (area integrals), `P01.08` (partial derivatives), `P00.09` (SymPy). Skill check: you can compute the determinant of a 2 × 2 Jacobian and evaluate ∮ p dq around an ellipse. |
| **Engineering question** | Orbit determination places a satellite in a 7000 km, e = 0.1 orbit to within about 10 m in position and 1 cm/s in velocity (illustrative). How does that cloud of possible states evolve: does its phase-space volume grow, or only its shape; how large is the along-track uncertainty after 1 and after 20 orbits; and which integrators carry such a cloud without creating or destroying phase-space volume? The answers: the volume is conserved (det Σ constant to the 10⁻⁵ accuracy of a finite-difference Jacobian); the along-track standard deviation grows from 10 m to 335 m in one orbit and 6.7 km in 20, linearly in time; Störmer–Verlet preserves area to round-off at any step, RK4 loses a little (1.4 × 10⁻⁵ in the pendulum test), explicit Euler gains 25 % in one second. |
| **Assumptions and limits** | Conservative, time-independent Hamiltonians (damping enters only as a deliberate counter-example); one-degree-of-freedom pictures for the pendulum, the planar two-body problem for orbits (no J₂, drag or third bodies); linearised (first-order) propagation of small orbit uncertainties; a blob's area measured as a polygon through its boundary points, accurate while the points stay dense along the stretched boundary. |
| **Analytical work** | The Legendre transform and Hamilton's equations; Poisson brackets and conservation laws; the pendulum's energy levels, separatrix and period by the elliptic integral K; the theorem T = dA/dE; Liouville's theorem from the divergence of the phase velocity; area contraction e^{−ct} under linear damping; the step Jacobians of Euler and Verlet; the Kepler Hamiltonian in polar coordinates, V_eff, the radial action J_r = 2π(μ/√(−2E) − p_θ) and dJ_r/dE = T; the along-track drift 3π δa √((1 + e)/(1 − e)) per orbit at periapsis. |
| **Computational deliverable** | `starter/main.py`: `legendre_transform`, `hamilton_equations`, `poisson_bracket`, `pendulum_hamiltonian`, `pendulum_flow`, `separatrix_momentum`, `pendulum_period`, `phase_area`, `flow_blob`, `euler_step`, `rk4_step`, `verlet_step`, `step_jacobian_det`, `kepler_hamiltonian`, `kepler_polar_hamiltonian`, `effective_potential`, `kepler_flow`, `lrl_vector`, `radial_action`, `flow_jacobian`, `along_track_drift` (`polygon_area`, `blob_boundary`, `iterate`, `symplectic_residual`, `periapsis_state`, `energy_error_history` and `main` are given). Tests prove: the Legendre transforms of the pendulum and of polar Kepler; Hamilton's equations equal Lagrange's; canonical brackets, {h, H} = 0 and {A_x, H} = 0; the separatrix lies on H = ω₀²; T(E) by K equals dA/dE at six energies (libration and rotation) to 10⁻⁶ and a simulated period to 10⁻⁸; the exact flow keeps a blob's area to 2 × 10⁻⁵ and damping shrinks it as e^{−ct}; Verlet's step determinant is 1, Euler's is 1 + h²ω₀² cos q, RK4's slightly below 1, with the matching blob areas; Kepler's energy and LRL vector conserved, |A| = μe; V_eff minimum at the circular orbit; J_r and dJ_r/dE = T; the flow Jacobian symplectic; the along-track drift to 10⁻³; Verlet's energy error bounded, RK4's growing. |
| **Practical activity** | **computational experiment**: the pendulum period from small amplitude to the separatrix by two routes; the area of a blob under three integrators and three step sizes; the shearing (filamentation) of a blob at constant area; the along-track drift against the size of the semi-major-axis error; the energy behaviour of Verlet and RK4 against step size; and the propagation of an orbit-determination covariance with the flow Jacobian. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `phase_space.png`, the derivations of Hamilton's equations for the pendulum and polar Kepler, of T = dA/dE and of Liouville's theorem, the experiment tables with your predictions, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-06/grading-rubric.md`; competencies 1–4. Pass: the Legendre transform and Liouville's theorem are derived (not quoted); all public tests pass; the period is checked by two routes and the area by the exact flow and three integrators with numbers; the note states the along-track growth and why volume is conserved while the cloud stretches. |
| **References** | `@libretexts-cline-variational-principles` Ch. 8 "Hamiltonian Mechanics": 8.2 Legendre Transformation between Lagrangian and Hamiltonian mechanics; 8.3 Hamilton's Equations of Motion; 8.4 Hamiltonian in Different Coordinate Systems; 8.5 Applications of Hamiltonian Dynamics; Ch. 3: 3.4 Geometrical Representations of Dynamical Motion (configuration, state and phase space; the plane-pendulum phase portrait); Ch. 15 "Advanced Hamiltonian Mechanics": 15.3 Canonical Transformations in Hamiltonian Mechanics; 15.7 Symplectic Representation; Ch. 11 "Conservative two-body Central Forces" (required-free). `@langtangen-linge-fd-computing-pdes` Ch. 1: 1.6 Energy Considerations; 1.7 The Euler-Cromer Method (1.7.4 The Störmer-Verlet Algorithm) (required-free). `@mit-ocw-18-03-differential-equations` notes 35 "Qualitative behavior of linear systems; phase plane"; notes 38 "Linearization near equilibria; the nonlinear pendulum" (supplementary-free). `@openstax-university-physics-vol1` 13.3 Gravitational Potential Energy and Total Energy; 13.5 Kepler's Laws of Planetary Motion (supplementary-free). |
| **Extensions** | ★ Action–angle variables for the pendulum: compute the angle variable from the time along an orbit and show that the blob of part 4 is a straight shear in (J, φ). ★ Propagate a 4-D cloud of 1000 sampled orbit states with `solve_ivp` and compare its sample covariance with M Σ₀ Mᵀ after 20 orbits; find the error size at which the linear prediction fails. ★ Add J₂ to the Kepler Hamiltonian and show with Poisson brackets which invariants survive. |
| **Time** | unmeasured |

## 1. The engineering problem

No spacecraft's state is known exactly. Orbit determination, from radar, GNSS or optical tracking, ends with a best estimate and a cloud of states that fit the data almost as well: a covariance. Mission design then asks how that cloud evolves. When will the uncertainty in where the satellite is along its orbit exceed the beam width of a ground antenna, or the separation used in a conjunction screening? Does the cloud grow without bound in every direction, or only in some? And can the numerical integrator used to propagate it be trusted to keep the cloud's size right over many orbits?

The natural language for these questions is **phase space**: the space of positions *and* momenta, in which a state is a point and a motion is a curve. Lagrange's equations of `P06.11` are second-order equations in configuration space; the **Hamiltonian** form rewrites them as first-order equations in phase space, q̇ = ∂H/∂p, ṗ = −∂H/∂q, with the energy function H as the generator of the motion. That form has a geometric consequence with practical weight: **Liouville's theorem**, which says that the flow carries any region of phase space to a region of the same volume. A cloud of states may stretch into a long thin filament, but it does not inflate or deflate. Integrators that respect this, the symplectic methods of `P03.07`, keep long-term statistics right; integrators that do not, slowly invent or destroy volume.

This lesson builds the Hamiltonian machinery in SymPy, tests it on the pendulum, whose phase portrait is the textbook picture of phase space (Cline 3.4), measures Liouville's theorem directly by carrying the boundary of a blob of states, and then turns to the Kepler problem in phase space. Success criteria, written first: the Legendre transforms reproduce the known Hamiltonians symbolically; the pendulum period by the elliptic integral equals dA/dE to 10⁻⁶; the exact flow keeps a blob's area to the polygon's resolution; damping shrinks it as e^{−ct} to the same accuracy; the Kepler flow Jacobian satisfies MᵀJM = J to the accuracy of its finite differences; and the along-track drift agrees with first-order theory to 10⁻³.

## 2. Understand and derive

### 2.1 From Lagrange to Hamilton

For a Lagrangian L(q, q̇) define the **canonical momenta**

**pᵢ = ∂L/∂q̇ᵢ.**   (1)

For a particle in Cartesian coordinates this is mẋ, the ordinary momentum; for an angle it is an angular momentum; in general it is whatever (1) gives. The **Legendre transform** (Cline 8.2) replaces the rates q̇ by the momenta p as the independent variables:

**H(q, p) = Σ pᵢ q̇ᵢ − L,  with q̇ expressed through p.**   (2)

This is the energy function h of `P06.11` §2.5, now written as a function of q and p; when T is quadratic in the rates, H = T + V. Its differential is dH = Σ(q̇ᵢ dpᵢ + pᵢ dq̇ᵢ) − Σ(∂L/∂qᵢ dqᵢ + ∂L/∂q̇ᵢ dq̇ᵢ). The dq̇ terms cancel by (1), and Lagrange's equation says ∂L/∂qᵢ = ṗᵢ, so dH = Σ(q̇ᵢ dpᵢ − ṗᵢ dqᵢ). Reading off the coefficients:

**q̇ᵢ = ∂H/∂pᵢ,  ṗᵢ = −∂H/∂qᵢ:**   (3)

**Hamilton's equations** (Cline 8.3), 2n first-order equations instead of n second-order ones.

**The pendulum.** L = ½ml²θ̇² + mgl cos θ gives p = ml²θ̇ and H = p²/(2ml²) − mgl cos θ; (3) gives θ̇ = p/(ml²), ṗ = −mgl sin θ, which is ml²θ̈ = −mgl sin θ. The code works per unit ml², with q = θ and p = θ̇, so H = p²/2 − ω₀² cos q with ω₀² = g/l.

**Kepler in polar coordinates.** For unit mass, L = ½(ṙ² + r²θ̇²) + μ/r gives p_r = ṙ and p_θ = r²θ̇ (the specific angular momentum), and

**H = p_r²/2 + p_θ²/(2r²) − μ/r.**   (4)

The angle θ does not appear: it is **cyclic**, so ṗ_θ = −∂H/∂θ = 0 and p_θ = h is constant (Kepler's second law, `P04.05`). The radial motion is then a one-degree-of-freedom system with Hamiltonian p_r²/2 + V_eff(r), with the **effective potential**

**V_eff(r) = h²/(2r²) − μ/r.**   (5)

Its minimum, at r = h²/μ with V_eff = −μ²/(2h²), is the circular orbit; a bound orbit oscillates between the two roots of V_eff(r) = E, periapsis and apoapsis (Cline 8.4 and Ch. 11 treat the coordinates and the central-force problem).

### 2.2 Poisson brackets and conservation laws

For any function f(q, p) along a motion, df/dt = Σ(∂f/∂qᵢ q̇ᵢ + ∂f/∂pᵢ ṗᵢ) = {f, H}, with the **Poisson bracket**

**{f, g} = Σᵢ (∂f/∂qᵢ ∂g/∂pᵢ − ∂f/∂pᵢ ∂g/∂qᵢ).**   (6)

So **f is conserved exactly when {f, H} = 0**. The canonical pairs satisfy {qᵢ, pⱼ} = δᵢⱼ. Three conservation laws checked symbolically in the tests: {H, H} = 0 (energy), {p_θ, H} = 0 in (4), and, for the Cartesian Kepler Hamiltonian, {h, H} = 0 for h = xp_y − yp_x and {A_x, H} = 0 for the **Laplace–Runge–Lenz vector**

**A = p × h − μ q/|q|,**   (7)

which points from the focus to periapsis and has magnitude μe. The LRL vector is the reason Kepler orbits close: its constancy fixes the direction of periapsis, so the ellipse does not precess. Add any perturbation (J₂, a third body, relativity) and {A, H} ≠ 0: the orbit's line of apsides turns.

### 2.3 The pendulum's phase portrait and its period

Each motion of the pendulum stays on a level curve of H (Cline 3.4):

- **libration**, −ω₀² < H < ω₀²: closed curves around (0, 0), the back-and-forth swing;
- **rotation**, H > ω₀²: open curves over the top, p never zero;
- the **separatrix**, H = ω₀²: p = ±2ω₀ cos(q/2), the two curves that join the upright equilibria (±π, 0) and take infinite time to traverse.

**The exact period.** With H = E on a libration curve, p = q̇ = √(2(E + ω₀² cos q)). Write E + ω₀² cos q = 2ω₀²(k² − sin²(q/2)) with k² = (E + ω₀²)/(2ω₀²) (k = sin(a/2) for an amplitude a), and substitute sin(q/2) = k sin φ. The quarter-period integral becomes ∫₀^{π/2} dφ/(ω₀√(1 − k² sin²φ)), the **complete elliptic integral of the first kind** K(k²)/ω₀ (`scipy.special.ellipk` takes the parameter m = k²):

**T_lib = 4 K(k²)/ω₀,  T_rot = 2 K(1/k²)/(k ω₀)** (one turn, k > 1).   (8)

As k → 0, K → π/2 and T → 2π/ω₀; as k → 1, K → ∞ logarithmically, and a pendulum released near the top takes ever longer to fall away.

**An independent route: T = dA/dE.** Let A(E) be the phase area enclosed by the orbit of energy E, A = ∮ p dq. Differentiating under the integral sign with ∂p/∂E = 1/p (from H = E) and p = dq/dt,

**dA/dE = ∮ dq/p = ∮ dt = T.**   (9)

This is a theorem about any one-degree-of-freedom Hamiltonian system, and it connects time (a dynamical quantity) to area (a geometric one). The code computes A(E) by quadrature of √(2(E + ω₀² cos q)) between the turning points and differentiates numerically; (8) uses a special function. They share nothing, and they agree to 10⁻⁶ at every energy tested, including rotations, where A is the area under one turn of the curve. For the harmonic oscillator, A = 2π(E − E_min)/ω₀ (an ellipse with semi-axes √(2ΔE)/ω₀ and √(2ΔE)), so dA/dE = 2π/ω₀: the period does not depend on the amplitude, as it must.

### 2.4 Liouville's theorem

Think of Hamilton's equations (3) as a velocity field v = (∂H/∂p, −∂H/∂q) in the (q, p) plane. The rate at which the flow of a field changes the area of a small patch is its divergence (`P01.12`):

**div v = ∂²H/∂q∂p − ∂²H/∂p∂q = 0.**   (10)

The mixed partials are equal, so the divergence vanishes identically: the phase flow of any Hamiltonian system is **incompressible**, and it carries every region to a region of the same area (in n degrees of freedom, the same 2n-dimensional volume). That is **Liouville's theorem**. It says nothing about shape. In the pendulum, an orbit near the separatrix takes longer than one near the bottom, so a blob of states that straddles several energies is sheared: its outer edge lags its inner edge, and over time the blob is drawn into a thin spiral filament whose area is still exactly the initial area.

**Damping breaks it in a predictable way.** With ṗ = −ω₀² sin q − cp, div v = −c everywhere, and the area obeys dA/dt = −cA:

**A(t) = A₀ e^{−ct},**   (11)

exactly, for any blob and any amplitude. A damped system contracts phase space, which is why it settles to an attractor; a Hamiltonian one cannot.

**Measuring an area.** `flow_blob` carries 2000–4000 points on the boundary of an ellipse of initial states along the exact flow with one `solve_ivp` call on the stacked state, and `polygon_area` applies the shoelace formula to the result. Two errors limit it: the solver's (10⁻¹⁰) and the polygon's chords, which cut corners on the curved, stretched boundary. The second dominates and grows as the boundary lengthens (part 4), so every area comparison is made against the exact flow carried with the same points.

### 2.5 One-step maps as area maps

`P03.07` showed that symplectic maps conserve a modified energy. In one degree of freedom, being symplectic is exactly preserving area: det J = 1 for the step's Jacobian J. For q̇ = p, ṗ = a(q) with step h:

- **explicit Euler**, q ← q + hp, p ← p + h a(q): J = [[1, h], [h a′(q), 1]], **det J = 1 − h² a′(q) = 1 + h²ω₀² cos q** for the pendulum. Near the bottom every step multiplies area by 1 + h²ω₀² > 1: Euler inflates phase space, and its energy spirals out.
- **Störmer–Verlet**: a half kick, a drift and a half kick, each with Jacobian of unit determinant ([[1, 0], [·, 1]] or [[1, h], [0, 1]]), so **det J = 1 exactly**, at any step size (Langtangen & Linge 1.7.4).
- **RK4**: det J is a polynomial in h that is not 1; for the harmonic oscillator it is |R(iωh)|² = 1 − (ωh)⁶/72 + …, slightly below 1, so RK4 slowly **contracts** phase space and its energy drifts down (`P03.07`).

`step_jacobian_det` computes det J by central differences; at the blob's centre with h = 0.05 it returns 1.024674 for Euler, 0.999999723 for RK4 and 1.000000000 for Verlet.

### 2.6 Kepler in phase space

**Radial action and period.** The radial motion in V_eff is a one-degree-of-freedom system, so (9) applies: the **radial action** J_r = ∮ p_r dr = 2∫ √(2(E − V_eff)) dr between periapsis and apoapsis has energy derivative equal to the radial period, which for Kepler is the orbital period. The integral has a closed form,

**J_r = 2π(μ/√(−2E) − p_θ),  dJ_r/dE = 2πμ(−2E)^{−3/2} = 2π√(a³/μ),**   (12)

with E = −μ/(2a): Kepler's third law from a phase area. For the 7000 km, e = 0.1 orbit, `radial_action` by quadrature gives 1.663633 × 10⁹ m²/s, equal to (12) to all digits printed, and its numerical energy derivative is 5828.517 s, the period.

**The flow Jacobian is symplectic.** For a 2-degree-of-freedom flow, Liouville's theorem is a consequence of a stronger property: the Jacobian M = ∂y(t)/∂y(0) of the flow satisfies **MᵀJM = J** with J = [[0, I], [−I, 0]] (Cline 15.7), and det M = 1 follows. `flow_jacobian` builds M by central differences of DOP853 runs. In SI units its entries span many orders of magnitude (metres per metre, metres per m/s of order the period), so the test and the printout use **canonical units** (length a, time 1/n, μ = 1); there, after one orbit, det M = 1.00000005 and max|MᵀJM − J| = 7.9 × 10⁻⁸, the finite-difference accuracy.

**An error that grows.** Two orbits that start together at periapsis with semi-major axes a and a + δa have periods differing by δT = (3/2)(δa/a)T. Each orbit the second arrives δT later and is a distance v_p δT behind. With v_p T = 2πa√((1 + e)/(1 − e)),

**along-track separation ≈ 3π δa n √((1 + e)/(1 − e))  after n orbits.**   (13)

A 10 m error in a, which a GNSS-tracked satellite might have, puts it 104.2 m behind after one orbit and 2084 m after 20. This is Liouville at work: the cloud of states shears along the orbit (along-track grows without bound, linearly), while the radial and energy spread stay fixed, so the volume is unchanged. It is the reason along-track error dominates every long orbit prediction.

### 2.7 Where the picture stops being valid

- **Non-conservative forces.** Atmospheric drag, thrust and outgassing make the flow non-Hamiltonian; phase volume then contracts (drag) or is changed by control, and (11)-like corrections apply.
- **Time-dependent and perturbed Hamiltonians.** With J₂ or third bodies, H is still Hamiltonian (Liouville holds), but the LRL vector and the closed form (12) do not.
- **Linear propagation.** Equation (13) and the covariance mapping Σ = MΣ₀Mᵀ are first order. When the along-track spread becomes a sizable fraction of the orbit, the cloud curves along the ellipse and a Gaussian description fails (part 4, experiment 4 shows the first-order law failing at 1.7 % for δa = 10 km).
- **The polygon.** A blob's boundary stretches; with a fixed number of points the polygon area error grows (experiment 3). Adaptive point insertion would fix it.

## 3. Build the computational model

Open `starter/main.py`. `polygon_area`, `blob_boundary`, `iterate`, `symplectic_residual`, `periapsis_state`, `energy_error_history` and `main` are given.

### Step 1: the symbolic machinery

`legendre_transform(L, qs, qdots, ps)` forms the momenta (1) with `sp.diff`, solves pᵢ = ∂L/∂q̇ᵢ for the q̇ᵢ with `sp.solve(..., dict=True)`, and substitutes into Σpq̇ − L. `hamilton_equations` is (3) and `poisson_bracket` is (6), both one line of `sp.diff`. Use plain symbols for q and q̇ here, not functions of t: H is a function on phase space, not along a motion.

### Step 2: the pendulum

`pendulum_hamiltonian`, `pendulum_flow` (which accepts a stack [q₁ … q_N, p₁ … p_N], so a whole blob moves in one `solve_ivp` call), `separatrix_momentum` and `pendulum_period` (8). `phase_area` integrates p dq with `scipy.integrate.quad` between the turning points ±arccos(−E/ω₀²), doubled, or over (−π, π) for a rotation. `flow_blob` is one DOP853 call with rtol 10⁻¹⁰.

### Step 3: maps and their Jacobians

`euler_step`, `rk4_step`, `verlet_step` take an acceleration function and arrays q, p, so they move every point of a blob at once. `step_jacobian_det` perturbs q and p by ±10⁻⁶ and forms the 2 × 2 Jacobian.

### Step 4: Kepler

`kepler_hamiltonian`, `kepler_polar_hamiltonian` (4), `effective_potential` (5), `kepler_flow`, `lrl_vector` (7) for one state or a stack, `radial_action` (quadrature between the roots of Er² + μr − p_θ²/2 = 0), `flow_jacobian` (central differences, perturbation 10⁻⁷ of each component's scale, rtol 10⁻¹²) and `along_track_drift` (13).

Running `main()` prints

```
pendulum: H = -g*l*m*cos(theta) + p_theta**2/(2*l**2*m)   Hamilton: ([p_theta/(l**2*m)], [-g*l*m*sin(theta)])
Kepler (polar): H = -mu/r + p_r**2/2 + p_theta**2/(2*r**2)   {p_theta, H} = 0

pendulum w0 = 3.14159 rad/s (small-angle period 2.000 s)
   amplitude   T elliptic (s)   T = dA/dE (s)
      10.0     2.003814         2.003814
      60.0     2.146364         2.146364
     120.0     2.745761         2.745761
     170.0     4.878725         4.878725

Liouville: blob area 0.8882640; after 6 s exact flow 0.8882630 (change -1.2e-06, polygon error); with damping c = 0.3/s 0.1468290 (A0 exp(-c t) = 0.1468291)
   Euler  h = 0.05: area after 1 s 1.1098860 (change +2.50e-01); one-step det at the blob centre 1.024674011
   RK4    h = 0.05: area after 6 s 0.8882510 (change -1.46e-05); one-step det at the blob centre 0.999999723
   Verlet h = 0.05: area after 6 s 0.8882630 (change -1.20e-06); one-step det at the blob centre 1.000000000

Kepler a = 7000 km, e = 0.1: period 5828.52 s; over 5 orbits (DOP853, rtol 1e-12):
   energy 8.5e-12, angular momentum 5.0e-12, LRL direction 1.2e-10 rad; |A|/mu = 0.100000
   radial action J_r = 1.663633e+09 m^2/s (closed form 1.663633e+09); dJ_r/dE = 5828.517 s vs period 5828.517 s
   flow Jacobian over one orbit (canonical units): det 1.0000000496, max|M^T J M - J| 7.9e-08
   da = 10 m: after  1 orbits along-track    104.2 m, total    104.6 m (first order    104.2 m)
   da = 10 m: after  5 orbits along-track    521.0 m, total    521.1 m (first order    521.0 m)
   da = 10 m: after 20 orbits along-track   2083.9 m, total   2083.9 m (first order   2083.9 m)
   RK4    200 steps/orbit: largest |energy error| in orbits 1-10 7.65e-08, in orbits 91-100 7.32e-07, at the end -7.28e-07
   Verlet 200 steps/orbit: largest |energy error| in orbits 1-10 1.06e-04, in orbits 91-100 1.06e-04, at the end +8.96e-07
```

and writes `phase_space.png`: the pendulum's phase portrait with the separatrix and the blob at 0 s and 6 s; V_eff with the orbit's energy; and the energy error of RK4 and Verlet over 100 orbits.

### Verify

Run `python -m pytest curriculum/phase-06/lesson-12-hamiltonian-mechanics`. The tests:

1. **Symbolic**: the pendulum and polar-Kepler Legendre transforms; Hamilton's equations equal Lagrange's; {q, p} = 1, {p, q} = −1; {p_θ, H} = 0 and {p_r, H} = h²/r³ − μ/r²; {h, H} = 0 and {A_x, H} = 0 for Cartesian Kepler.
2. **Pendulum**: H = ω₀² on the separatrix; T by (8) equals dA/dE to 10⁻⁶ at six energies; the small-amplitude limit 2π/ω₀; the period of a simulated 120° swing (from `solve_ivp` events) equals (8) to 10⁻⁸; A = 2πΔE/ω₀ for small swings; energies below the minimum are refused.
3. **Liouville**: the blob's area kept by the exact flow to 2 × 10⁻⁵ while its points move by more than 0.5 rad; damping gives A₀e^{−ct} to the same accuracy.
4. **Maps**: det J = 1 for Verlet and 1 + h²ω₀² cos q for Euler at random points, slightly below 1 for RK4; after 6 s Verlet's blob area equals the exact flow's to 10⁻⁶, RK4's is smaller by more than five times that, and Euler's is 15 % larger after 1 s.
5. **Kepler**: E = −μ/2a, energy conserved to 10⁻¹⁰, |A| = μe and A along +x; polar and Cartesian H agree; V_eff minimum; J_r = (12) to 10⁻⁶ and dJ_r/dE = T to 10⁻⁵; the flow Jacobian has det 1 and MᵀJM = J to 10⁻⁵ while having entries above 5; the along-track drift (13) to 10⁻³; Verlet's energy error has the same peak in the last 5 orbits as in the first 5, RK4's more than four times larger.
6. `main` writes the figure.

### Compare

The independent checks are internal, by design: the period by a special function (8) against the period by an area derivative (9) against the period of a `solve_ivp` run; the area of the exact flow against Liouville's theorem and against the e^{−ct} law with damping; the radial action by quadrature against its closed form (12); MᵀJM = J for a flow Jacobian built from nothing but finite differences of an integrator. Every pair agrees to the accuracy of the cheaper of the two. The integrator comparison repeats `P03.07` from the area side: Verlet's energy error swings by 1.06 × 10⁻⁴ at 200 steps per orbit and does not grow (1.06 × 10⁻⁴ in orbits 1–10 and in orbits 91–100), RK4's is smaller but grows tenfold from orbit 10 to orbit 100, the same −7.3 × 10⁻⁹ per orbit that `P03.07` measured.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** All results here are outputs of the model. Write each prediction before running it.

1. **Period to the separatrix.** For amplitudes 1°, 10°, 30°, 60°, 90°, 120°, 150°, 170°, 179°, compute T by (8), by dA/dE, and by the two-term series 2π/ω₀ (1 + a²/16) of `P04.06`. *Finding:* 2.000038, 2.003814, 2.034818, 2.146364, 2.360681, 2.745761, 3.524407, 4.878725, 7.80213 s; dA/dE agrees to six figures throughout (7.802135 s at 179°, where quadrature near the turning point starts to cost a digit); the series is within 0.03 % to 30°, 2.2 % low at 90°, and useless beyond (3.22 s against 7.80 s at 179°). In rotation (H = 1.5, 3, 10 ω₀²) one turn takes 1.2853, 0.8346, 0.4481 s, and dA/dE agrees to 10⁻⁹ relative.
2. **Area under three integrators.** The blob of 4000 points centred on (0, 1.2ω₀) with semi-axes 0.3 rad and 0.3ω₀; step sizes 0.025, 0.05, 0.1 s; Euler over 1 s, RK4 and Verlet over 6 s, each compared with the exact flow. Predict the sign of each change. *Finding:* Euler +13.4 %, +25.0 %, +39.1 %; RK4 −4.2 × 10⁻⁷, −1.35 × 10⁻⁵, −4.15 × 10⁻⁴ (roughly h⁵ over a fixed time, consistent with a per-step defect of order h⁶); Verlet −7.5 × 10⁻⁹, −3.1 × 10⁻⁸, −1.3 × 10⁻⁷ (the polygon's own error, which changes slightly with where the points end up). Euler's one-step determinant at the blob centre is 1.006169, 1.024674, 1.098696.
3. **Shear at constant area.** Carry the same blob for 2, 6, 12, 24 s with the exact flow and record its perimeter (in q and p/ω₀) and area. *Finding:* perimeter 1.885, 2.980, 6.334, 12.378, 24.310; area 0.888264, 0.8882639, 0.888263, 0.8882597, 0.8882463. The perimeter grows linearly with time (the period varies across the blob), the area stays put; its apparent loss of 2 × 10⁻⁵ at 24 s is the polygon's chords, which grows with the perimeter. A harmonic oscillator would rotate the blob rigidly (all periods equal) with no stretching at all.
4. **Along-track drift.** Two orbits starting at periapsis with δa = 1, 10, 100, 1000, 10 000 m; along-track separation after 20 orbits against (13). *Finding:* 208.39, 2083.9, 20 838.9, 208 347.7, 2 049 063 m against 208.39, 2083.9, 20 839.0, 208 389.7, 2 083 897 m: relative errors −2.5 × 10⁻⁷, −4.1 × 10⁻⁷, −5.2 × 10⁻⁶, −2.0 × 10⁻⁴, −1.7 × 10⁻². The relative error of the first-order law grows faster than δa itself: by δa = 10 km the two satellites are 2000 km apart, a sizeable arc of the ellipse, and a straight-line, first-order picture no longer fits.
5. **Energy against step.** Verlet and RK4 at 100, 200, 400 steps per orbit for 40 orbits. *Finding:* Verlet swings 4.25 × 10⁻⁴, 1.06 × 10⁻⁴, 2.64 × 10⁻⁵ (order 2, bounded); RK4 drifts −2.33 × 10⁻⁷, −7.28 × 10⁻⁹, −2.28 × 10⁻¹⁰ per orbit (a factor of 32 per halving of the step), secular.
6. **An orbit-determination cloud.** Take a diagonal covariance Σ₀ with σ = 10 m in each position component and 0.01 m/s in each velocity component at periapsis, map it with the flow Jacobian in canonical units, Σ(t) = M Σ₀ Mᵀ, after 1, 5 and 20 orbits. *Finding:* radial σ stays 10.0 m; along-track σ 335 m, 1672 m, 6689 m (linear in time; the position uncertainty contributes more than the velocity uncertainty, 257.5 m against 213.7 m after one orbit, almost all of the former from the radial 10 m, and the two add in quadrature); det Σ/det Σ₀ = 1.0000001, 1.0000010, 1.0000074, conserved to the accuracy of the finite-difference Jacobian. A 1 cm/s velocity error matters almost as much as 10 m of radial position: through the energy, each sets δa, and δa sets the drift.

Raw data format: CSV tables `amplitude_deg, T_elliptic_s, T_dAdE_s, T_series_s` (1); `h_s, euler_change, rk4_change, verlet_change, euler_det` (2); `t_s, perimeter, area` (3); `da_m, along_sim_m, along_first_order_m, rel_error` (4); `steps_per_orbit, verlet_swing, rk4_drift_per_orbit` (5); `orbits, sigma_radial_m, sigma_along_m, det_ratio` (6); plus `phase_space.png`.

## 5. Validate and improve

| quantity | analysis (§2) | computation | difference | cause |
|---|---|---|---|---|
| pendulum period, 120° (s) | (8): 2.745761 | dA/dE 2.745761; `solve_ivp` events (test) | < 10⁻⁶; < 10⁻⁸ | quadrature; solver |
| blob area, exact flow, 6 s | A₀ (Liouville) = 0.8882640 | 0.8882630 | 1.2 × 10⁻⁶ | polygon chords |
| blob area, damped c = 0.3/s, 6 s | A₀e^{−1.8} = 0.1468291 | 0.1468290 | 7 × 10⁻⁷ | polygon chords |
| one-step area factor, Euler, h = 0.05 | 1 + h²ω₀² cos 0 = 1.024674 | 1.024674011 | < 10⁻⁹ | — |
| one-step area factor, Verlet | 1 | 1.000000000 | < 10⁻⁹ | — |
| radial action (m²/s) | (12): 1.663633 × 10⁹ | quadrature 1.663633 × 10⁹ | < 10⁻⁶ | quadrature |
| dJ_r/dE (s) | T = 5828.517 | 5828.517 | < 10⁻⁶ | — |
| flow Jacobian, one orbit | MᵀJM = J, det 1 | residual 7.9 × 10⁻⁸, det 1.00000005 | 10⁻⁷ | finite differences |
| along-track, δa = 10 m, 20 orbits (m) | (13): 2083.9 | 2083.9 | 4 × 10⁻⁷ | second order in δa/a |

Every Hamiltonian result is confirmed by an independent route to the accuracy of the method, and the two required checks of this lesson hold: Liouville's area preservation (exact flow and Verlet) against the area change of a non-symplectic method, and bounded symplectic energy error against RK4 drift, with the Kepler invariants (energy, angular momentum, the LRL vector) conserved to 10⁻¹⁰. The engineering answer: **a 10 m, 1 cm/s orbit-determination cloud keeps its phase-space volume but shears along the orbit, its along-track σ growing linearly to 335 m after one orbit and 6.7 km after 20; a symplectic integrator carries it without changing its volume, explicit Euler inflates it at a rate that depends on the step, and RK4 shrinks it slowly.** For a covariance propagated over weeks, the linearisation (experiment 4) fails before the integrator does.

What to improve: adaptive insertion of boundary points when the polygon stretches (experiment 3), and a sampled (Monte Carlo) cloud to test the linear covariance mapping (Extension 2). What not to improve: the integration tolerance, which is already far below the polygon error.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (8), (9), (10), (11) and (13); `phase_space.png`; the six experiment tables with your predictions; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the damping goes to zero?** Equation (11) gives A = A₀; the damped and undamped tests are the same test. What if the step goes to zero? Euler's determinant tends to 1 and its area error to zero, but only linearly in h over a fixed time (experiment 2: +13 % at h = 0.025 s).
- **Which assumption fails first in the real world?** For low orbits, drag: it makes the flow non-Hamiltonian, shrinks the semi-major axis and makes along-track error grow quadratically in time, not linearly. For the covariance, the linearisation, once the along-track spread is a measurable fraction of the orbit.
- **What measurement could falsify the model?** Repeated orbit determinations of the same satellite: the along-track difference between a prediction and the next fit should grow linearly with time at the rate (13) predicts from the first fit's δa. A quadratic growth means an unmodelled force (drag); a periodic error means a missing perturbation (J₂).
- **Does the result survive a change of method?** The pendulum period survives three methods; the blob area survives a change from 2000 to 8000 points at the 10⁻⁶ level for the exact flow and Verlet; the RK4 area loss changes by the predicted power of h. A result that depended on the polygon's resolution would be a polygon artefact, not physics.
- **Which decision changes if the uncertainty doubles?** Doubling σ_v to 2 cm/s doubles the velocity part of the along-track growth rate (it is linear in δa), and the total by a factor of 1.49 (the position part is unchanged and the two add in quadrature): a conjunction screening threshold reached after 20 orbits would be reached after about 13.4, and the operator would need a fresh orbit determination half as often again. Doubling both σ_r and σ_v doubles the rate and halves the interval.

## Check yourself

1. Define the canonical momentum of a coordinate q. What is it for the angle of a pendulum, and for the polar angle of an orbit?
2. Write Hamilton's equations for H = p²/(2m) + ½kq² and show they give the spring–mass equation.
3. What does {f, H} = 0 mean physically? Give two Kepler examples.
4. Sketch the pendulum's phase portrait and label libration, rotation and the separatrix. Which motion has an infinite period?
5. State the relation between the period of a one-degree-of-freedom orbit and the phase area it encloses, and use it to show that a harmonic oscillator's period does not depend on its amplitude.
6. State Liouville's theorem. Why does a blob of pendulum states stretch into a filament without changing its area?
7. What is the area factor per step of explicit Euler for the pendulum, and of Störmer–Verlet? What does each do to the energy over many steps?
8. Why does along-track uncertainty grow linearly in time while radial uncertainty stays bounded for a Kepler orbit?
