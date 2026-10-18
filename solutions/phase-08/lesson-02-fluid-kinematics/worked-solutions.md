# P08.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown.

## 1. The centreline of a wind-tunnel contraction (s = 15: L = 3.5 m)

**Assumptions.** The model field is 2-D, steady, and valid near the axis only.

(a) G = [[∂u/∂x, ∂u/∂y], [∂v/∂x, ∂v/∂y]] = [[U₀/L, 0], [0, −U₀/L]] = **[[1.4286, 0], [0, −1.4286]] s⁻¹**. Divergence **0** (incompressible), vorticity **0** (irrotational), S = G (pure strain, no rotation), principal strain rates **±1.4286 s⁻¹**: stretching along the axis, squeezing across it. That is what a contraction does to a fluid element, and why it straightens turbulent eddies.

(b) Du/Dt = ∂u/∂t + u ∂u/∂x + v ∂u/∂y = 0 + U₀(1 + x/L) U₀/L = **9.184 m/s²**; Dv/Dt = v ∂v/∂y = (−U₀y/L)(−U₀/L) = **0.408 m/s²**. Both are convective: the flow is steady, yet each particle accelerates because it moves to where the velocity is larger. `particle_acceleration` gives (9.1837, 0.4082) m/s².

(c) dx/dt = U₀(1 + x/L) gives 1 + x/L = exp(U₀t/L), so the axis particle reaches x = L at **t = (L/U₀) ln 2 = 0.4852 s** (`pathline` with 200 steps ends at x = 3.500000 m). Along a pathline, d/dt[y(1 + x/L)] = v(1 + x/L) + y u/L = −U₀(y/L)(1 + x/L) + (y/L)U₀(1 + x/L) = 0. Because the flow is steady, pathlines are streamlines: the streamlines are the hyperbolas **y(1 + x/L) = const**, the stagnation-point flow of the lesson centred at x = −L.

**Common mistakes.** Calling the acceleration zero because the flow is steady; forgetting the cross-stream component.

## 2. Lines and derivatives

(a) **Pathline:** x = x₀ + U(t − t₀), y = y₀ + (V/ω)(sin ωt − sin ωt₀). **Streamline** at time t: dy/dx = v/u = (V/U) cos ωt, a straight line y − y₀ = (V cos ωt/U)(x − x₀). **Streakline** at time t: the particle released at τ is at x = x_s + U(t − τ), y = y_s + (V/ω)(sin ωt − sin ωτ). With τ = t − (x − x_s)/U:

**y − y_s = (V/ω)[sin ωt − sin ω(t − (x − x_s)/U)],**

a sine wave in x of wavelength **2πU/ω**, frozen in the smoke and carried downstream at U. Its slope at the source is dy/dx = (V/ω)(ω/U) cos ωt = **(V/U) cos ωt**, the slope of the streamline through the source at that instant: the smoke leaves the wand along the current flow direction and then keeps the direction it was given.

(b) In a steady flow, dx/dt = u(x) is autonomous: its solution curves do not depend on the starting time, and with a smooth u there is exactly one curve through each point (uniqueness). A pathline is such a curve; a streamline at any frozen instant is the same curve parametrised by arc length; every particle released at the source travels along it, so the streakline lies on it too. Uniqueness is the property needed.

(c) Along a pathline (x(t), y(t)), d/dt f(x(t), y(t), t) = ∂f/∂t + (∂f/∂x)ẋ + (∂f/∂y)ẏ = **∂f/∂t + u ∂f/∂x + v ∂f/∂y = Df/Dt**.

(d) A rectangle with edges δx along x and δy along y: the far end of the δx edge moves faster than the near end by (∂u/∂x)δx in x and (∂v/∂x)δx in y. So d(δx)/dt = (∂u/∂x)δx, and the edge turns anticlockwise at ∂v/∂x. Similarly d(δy)/dt = (∂v/∂y)δy and the δy edge turns at −∂u/∂y. Then **d(δA)/dt = (∂u/∂x + ∂v/∂y)δA = (∇·u)δA**, and the average rotation of the two edges is **(∂v/∂x − ∂u/∂y)/2 = ω/2**. The difference of the two turning rates is the shear strain rate, which is the off-diagonal part of S.

**Check.** The oscillating-stream formulas are what `oscillating_pathline_exact` and `oscillating_streakline_exact` implement; RK4 reproduces them to 4 × 10⁻¹⁰ m (lesson part 3).

**Common mistakes.** Writing the streakline as the pathline from the source (that is one particle, not the locus of all of them); dropping the minus sign in the rotation of the δy edge.

## 3. A smoke plume in a veering wind (s = 15: T = 750 s)

**Assumptions.** Uniform wind, smoke as passive tracers, no diffusion.

(a) Ω = 2π/750 = 8.378 × 10⁻³ s⁻¹, R = W/Ω = **596.8 m**. **Pathline:** x = R sin Ωt, y = R(1 − cos Ωt): a **circle of radius R centred at (0, R)**. **Streamline** at any time t: the wind is the same everywhere, so it is a **straight line through the chimney along (cos Ωt, sin Ωt)**.

(b) The particle released at τ is at x = R(sin Ωt − sin Ωτ), y = R(cos Ωτ − cos Ωt). Then (x − R sin Ωt)² + (y + R cos Ωt)² = R²(sin²Ωτ + cos²Ωτ) = R²: a **circle of radius R centred at (R sin Ωt, −R cos Ωt)**, which at t = T/4 is **(596.8 m, 0)**. The smoke occupies the quarter of it from the chimney to the oldest particle at (596.8, 596.8) m, where it meets the pathline (the oldest smoke particle is the one whose pathline (a) traces).

(c) RK4 with 300 steps: the pathline points lie on their circle and the streakline points on theirs to **1.6 × 10⁻¹⁰ m**; the streamline at T/4 points along +y to 3 × 10⁻¹³ m in 1000 m.

(d) Near the chimney the drawing is right: the streakline's tangent at the source is (cos Ωt, sin Ωt), the present wind. Far from it the smoke curves back towards the directions the wind had when it was released; the photograph shows the wind's history, not its present field. Here the present wind is uniform, but the smoke is a circle.

**Common mistakes.** Expecting a curved plume to mean a curved wind field; computing the streakline as one particle's path.

## 4. When does smoke show streamlines? (s = 15: U = 18 m/s, f = 2.5 Hz)

(a) V = 18 tan 5° = **1.5748 m/s**, ω = **15.708 rad/s**, wavelength 2πU/ω = **7.20 m**, pathline amplitude V/ω = **100.3 mm**.

(b) Worst-phase error over one gust period (121 phases):

| L (m) | k = ωL/U | smoke error |
|---|---|---|
| 0.1 | 0.087 | 0.0436 |
| 0.2 | 0.175 | 0.0872 |
| 0.3 | 0.262 | 0.1306 |
| 0.5 | 0.436 | 0.2170 |
| 1.0 | 0.873 | 0.4271 |

(c) The error reaches 5 % at **L = 0.1146 m, k = 0.1000**. Expanding the streakline for small ξ = x − x_s: sin ω(t − ξ/U) = sin ωt − (ωξ/U) cos ωt − ½(ωξ/U)² sin ωt + …, so y_streak − y_stream = (V/ω)·½(ωξ/U)² sin ωt + O(ξ³). The worst phase (sin ωt = ±1) at ξ = L, divided by the streamline's excursion VL/U, gives **error ≈ k/2**, and k/2 = 0.05 gives L = 0.1U/ω = 0.1146 m: the estimate is exact to three figures here.

(d) No: a 0.3 m model has k = 0.26 and a 13 % error, so the smoke is a record of the gust's history over the model, not its instantaneous streamlines. Use phase-locked particle image velocimetry, which measures the instantaneous velocity field, and compute streamlines from it; or keep the smoke and compare it with computed streaklines, not streamlines.

**Common mistakes.** Using one phase of the gust (the error depends on phase and vanishes at some); confusing the reduced frequency k = ωL/U with the Strouhal number based on frequency in hertz.

## 5. Stretching a material loop

(a) In the stagnation flow at t = 2 s: area ratio **1.0000000000028**, perimeter ratio **4.708**, x-extent stretch **7.389**, equal to e^{at} = e² = 7.389 to 10⁻⁹. The loop becomes an ellipse with semi-axes 0.1e² and 0.1e⁻²: same area, much longer perimeter.

(b) With u = ax, v = 0 the area ratio is **7.389 = e²**: tr G = **1.000 s⁻¹** from `decompose_gradient`, and d(δA)/dt = (∇·u)δA integrates to e^{at}. This field is not divergence-free; a gas being compressed or expanded behaves like this.

(c) Area preservation is exact for the flow map of a divergence-free field, and RK4 tracks each point to about 10⁻¹² here; the shoelace area of 800 points along an ellipse whose curvature stays mild is accurate to the same order. As the ellipse thins further, the ends become sharply curved and the straight segments between points cut the corners: the representation of the loop, not the flow, limits the accuracy. Inserting points where segments grow (adaptive redistribution) cures it.

**Common mistakes.** Expecting the perimeter to be conserved; confusing incompressible (area-preserving) with rigid (shape-preserving).

## 6. Which line does a flow-visualisation method show?

(a) Smoke wand, steady tunnel: a streakline, which equals the streamline in steady flow.
(b) With the gust generator running: a streakline, which differs from the streamlines by about k/2 of the deflection (Problem 4).
(c) Long exposure of bubbles: pathline segments (each streak is one bubble's path during the exposure).
(d) PIV: the instantaneous velocity field (short pathline segments much shorter than any flow scale), from which streamlines at that instant can be computed.
(e) Wool tufts: the local flow direction at the surface at an instant, the near-wall streamline direction (with the tuft's own lag and flutter).
(f) Dye from a nose port while the model pitches: a streakline in an unsteady flow; it shows where fluid that passed the port has gone, which near separation is very different from the streamlines.

## Check-yourself answers

1. A streamline is everywhere tangent to the velocity at one instant; a pathline is the track of one particle over time; a streakline is the locus at one instant of all particles that passed a fixed point.
2. In a steady flow (and in flows such as the decaying Taylor–Green vortex whose velocity direction never changes at any point).
3. The streakline is what smoke or dye from a fixed source shows.
4. Df/Dt = ∂f/∂t + u·∇f: the rate of change of f following a fluid particle; the first term is the local change at a fixed point, the second the change from moving to where f is different.
5. Yes: in a steady flow the convective acceleration (u·∇)u is nonzero wherever the speed or direction changes along a streamline, as in a contraction or round a bend.
6. The divergence ∇·u is the rate of change of volume (area in 2-D) per unit volume of a fluid element; it is zero for an incompressible flow.
7. The strain-rate tensor is the symmetric part of the velocity gradient (stretching and shearing); the rotation tensor is the antisymmetric part; the vorticity is twice the local angular velocity, ω = ∂v/∂x − ∂u/∂y in 2-D.
8. Simple shear u = γy has vorticity −γ and strain rates ±γ/2 along the 45° diagonals: half its velocity gradient deforms the element and half rotates it.
