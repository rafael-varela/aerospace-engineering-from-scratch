# P08.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact solutions of the Navier–Stokes equations, which are approximations (a narrow gap treated as a plane slot, a typical fluid property), and which are numerical.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

Fluid properties (typical values, not data-sheet values): hydraulic fluid ρ = 845 kg/m³, ν = 14 mm²/s (μ = 0.01183 Pa·s); bearing oil ρ = 870 kg/m³, μ = 0.03 Pa·s.

## 1. (Warm-up: a turbopump journal bearing)

A shaft of diameter 30 mm turns at N = 6000 + 200s rpm in an oil-lubricated journal bearing 50 mm long with a radial clearance of 25 µm. Treat the concentric film as plane Couette flow (the gap is 1/600 of the radius).

(a) Find the surface speed U, the shear stress in the film, the friction torque and the power dissipated.
(b) Compute the gap Reynolds number ρUh/μ and the Taylor number Ta = Re²(h/R). Laminar Couette flow between a rotating inner and a fixed outer cylinder becomes unstable to Taylor vortices at Ta ≈ 1700. Is the Couette solution valid here?
(c) Where does the dissipated power go, and why does it matter for the viscosity you assumed?

Answer form: four numbers with units; two numbers and a sentence; two sentences.

## 2. (Derivation: three exact solutions)

(a) Starting from the incompressible Navier–Stokes equations in 2-D, state the assumptions (steady, fully developed, no body force along x, walls at y = 0 and y = h) that reduce them to μ u″ = dp/dx. Show from continuity and the no-slip condition that v = 0, and from the y-momentum equation that p depends on x only. Integrate to obtain the Couette–Poiseuille profile u(y) = Uy/h + Gy(h − y)/(2μ), its flow rate per unit width and its two wall shear stresses. Check that τ(h) − τ(0) = −Gh is the momentum balance of a slab of fluid.
(b) Write the axial momentum equation for steady, fully developed flow in a round pipe in cylindrical coordinates, solve it with a finite velocity on the axis and no slip at r = R, and derive Q = πR⁴G/(8μ). Define the Darcy friction factor f = 8τ_w/(ρV²) and show f = 64/Re with Re = ρVD/μ.
(c) For Stokes's first problem (u_t = ν u_yy, u(0, t) = U, u(∞, t) = 0, u(y, 0) = 0) look for u = U F(η) with η = y/(2√(νt)). Show that F″ + 2ηF′ = 0 and that F = erfc η. Why can the solution depend on y and t only through η?

Answer form: three derivations with the assumptions listed where they are used.

## 3. (Unfamiliar scenario: a de-icing fluid film)

Before take-off, a wing panel inclined at θ = 10 + s/2 degrees to the horizontal carries a film of de-icing fluid (ρ = 1050 kg/m³, μ = 0.01 Pa·s, typical) draining under gravity at q = 1.0 × 10⁻⁵ m² /s per metre of span. The free surface carries no shear stress (the air above is still).

(a) Set up the problem in coordinates along (x) and normal to (y) the panel. Which body-force components appear in the two momentum equations? What is the boundary condition at the free surface y = H?
(b) Solve for u(y) and the pressure p(y), and express the film thickness H in terms of q.
(c) Find H, the surface speed and the film Reynolds number ρq/μ for your θ. Is the film laminar?
(d) Use `ns_residuals` (with the gravity components as gx and gy) to confirm that your u and p are an exact solution, and check that the free surface is shear-free.

Answer form: two equations and a boundary condition; a profile, a pressure and a formula; three numbers; a SymPy output.

## 4. (Numerical, uses your code: your own actuator)

Take the actuator of the lesson (D = 25 mm, land length 12 mm, Δp = 20.7 MPa) with a radial clearance h = 4 + s/5 µm.

(a) With `actuator_leakage`, find the leakage (in L/min), the gap Reynolds number, the shear force on the piston and the start-up time for the piston at rest and moving at ±0.1 m/s. Which part of the force changes with the piston speed?
(b) At what piston speed (sign and size) would the Couette flow cancel the leak exactly? Is that speed reachable in service?
(c) The clearance is held to ±1 µm in manufacture. By what factors does the leakage change at the two tolerance limits? What does this say about specifying a clearance?
(d) Check the start-up time with `crank_nicolson` (40 intervals across the gap, 400 time steps, from rest): what fraction of the steady centreline speed is reached at t₉₉?

Answer form: a table of four quantities for three speeds and a sentence; a speed and a sentence; two factors and a sentence; a fraction.

## 5. (Numerical, extends your code: Stokes's second problem)

A plate at y = 0 oscillates in its own plane with u(0, t) = U cos ωt under fluid at rest far away. The periodic exact solution is u = U e^{−ky} cos(ωt − ky) with k = √(ω/(2ν)).

(a) Substitute the solution into u_t = ν u_yy (by hand or with `ns_residuals`) and confirm it is exact. How deep does the oscillation reach before its amplitude falls to 1 %?
(b) With ν = 1, ω = 2π and a domain 0 ≤ y ≤ 10/k, start `crank_nicolson` from the exact solution at t = 0 and march three periods with u_hi set to the exact value. Show that the error converges at second order when n = 40, 80, 160 intervals and 2n time steps are used.
(c) Now start from rest (u = 0) with u_hi = 0. How many periods does it take for the difference from the periodic solution to fall below 10⁻³ on n = 200 intervals with 400 steps per period?
(d) A hydraulic actuator's spool dithers at 10 Hz and 100 Hz. How far into the hydraulic fluid does the oscillating shear layer reach? Compare with a 10 µm clearance: is the gap flow quasi-steady?

Answer form: a residual and a depth; three errors and two orders; a number of periods; two depths and a sentence.

## 6. (Judgment: when does the exact solution stop describing the actuator?)

The lesson's leakage of 0.057 L/min comes from the plane Couette–Poiseuille solution with a constant viscosity, a concentric piston and rigid walls.

(a) The fluid that leaks drops 20.7 MPa with no work done. If all that energy stays in the leaking fluid (c_p ≈ 1900 J/(kg·K), typical), by how much does it warm? Which way does that move the leakage, and why?
(b) A fully eccentric piston (touching the bore on one side) leaks (1 + 1.5ε²) = 2.5 times as much as a concentric one at ε = 1, a standard lubrication result. What does that do to your confidence in the leakage figure?
(c) The developing length of a laminar slot flow is about 0.06 Re h. Is the fully developed assumption justified?
(d) Rank (a)–(c) and the uncertainty in h from Problem 4(c) by their effect on the leakage, and say which you would measure first.

Answer form: a temperature rise and a direction; a factor and a sentence; a length and a verdict; a ranked list.
