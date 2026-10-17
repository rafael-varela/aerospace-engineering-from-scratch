# P01.12  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical, with the resolution used.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: divergence, curl and flux by hand)

Let **v** = (s x, −y, (1 − s) z) and **F** = (y², z², x²), both in SI units with coordinates in metres.

(a) Compute ∇ · **v**, ∇ × **v**, ∇ · **F** and ∇ × **F** by hand and check them with `divergence_sym` and `curl_sym`.
(b) Without integrating anything, state the outward flux of **v** and of **F** through the unit cube [0, 1]³. Check with `box_flux` at n = 3.
(c) For **G** = (s x², 0, 0), compute the flux through the unit cube by hand on the two x faces and the volume integral of ∇ · **G**. Run `divergence_theorem_check` at n = 2 and explain the size of the residual.

Answer form: four expressions; two numbers with a one-line reason each; two numbers and one sentence.

## 2. (Derivation: Green for a rectangle, and the continuity equation)

(a) Derive Green's theorem ∮(P dx + Q dy) = ∬(∂Q/∂x − ∂P/∂y) dA for the rectangle [a, b] × [c, d] from the fundamental theorem of calculus, handling the P and Q terms separately and stating where the counter-clockwise orientation enters.
(b) Starting from "the mass inside a fixed control volume changes only because mass crosses its surface", derive ∂ρ/∂t + ∇ · (ρ**v**) = 0. Name the step that uses the divergence theorem and the step that uses "for every control volume".
(c) Show that for steady flow the mass flow rate ρ v A is the same at every cross-section of a stream tube (a tube whose side walls are everywhere parallel to **v**). Which faces of the control volume carry flux and why?

Answer form: two derivations with named steps; one short argument.

## 3. (Unfamiliar scenario: weighing the Earth with a flux)

Inside a sphere of uniform density ρ the gravitational acceleration is **g** = −(4/3)πGρ **r** (G = 6.674 × 10⁻¹¹ m³ kg⁻¹ s⁻²).

(a) Compute ∇ · **g** by hand and show that the flux of **g** through any closed surface inside the body is −4πG × (mass enclosed). This is Gauss's law for gravity.
(b) Take ρ = 5515 kg/m³ and a 1 km cube at a depth of 3000 km (any position inside the body). Compute the flux of **g** through it with `box_flux` and `uniform_sphere_gravity`, divide by −4πG and compare with the mass of rock in the cube.
(c) The measured surface gravity is 9.82 m/s² and R = 6371 km. Using **g** at the surface, infer Earth's mean density and compare with the density of crustal rock, 2700 kg/m³. What does the comparison tell you about the interior?
(d) In the uniform model, what is g at a depth of 3000 km? Real seismic models give about 10 m/s² there, slightly *more* than at the surface. Which assumption of the uniform model is responsible for the difference?

Answer form: a derivation; three numbers with units; two short paragraphs.

## 4. (Numerical, uses your code: residuals and where they come from)

(a) Run `divergence_theorem_convergence` for **F** = (y eˣ, z sin y, x cos z) on the contraction box ((−0.2, 0.2), (−0.2, 0.2), (1, 2)) with n = 2, 4, 8, 16, 32, 64. Tabulate, fit the slope, and compare the n = 64 volume integral with `scipy.integrate.tplquad`.
(b) Repeat for **F** = (x²y, yz², z sin x) and for **F** = (x³, 0, 0). One of them has a residual at rounding level for every n. Prove, for the x component alone, that the discrete face-flux difference equals the discrete volume sum of ∂F_x/∂x exactly whenever F_x is at most quadratic in x, and explain why the cubic field breaks it.
(c) Run `greens_theorem_check` for the Rankine vortex (Γ = 300 m²/s, a = 1.5 m) on the 6 m square with n = 25, 50, 100, 200, 400, 800. Tabulate the errors of both sides against Γ, fit both slopes and explain why one of them is noisy.

Answer form: three tables with slopes; a short proof; one paragraph.

## 5. (Stokes on two surfaces and one closed surface)

Use **F** = (y eˣ, z sin y, x cos z) with its SymPy curl, R = 1.3 m.

(a) Run `stokes_theorem_check` for n = 10, 20, 40, 80. Tabulate rim, cap and disc.
(b) The cap and the disc share the rim. Form the closed surface "cap plus disc with its normal reversed" and show that the flux of ∇ × **F** through it, cap − disc, tends to zero at second order. Which identity makes it exactly zero in the continuum?
(c) By hand, for the rigid rotation (−ωy, ωx, 0), compute the flux of the curl through the upper hemisphere of radius R and show it equals 2πωR². (You need ∬_hemisphere z dA = πR³.)
(d) Suppose you had parametrised the hemisphere with s and t swapped. What would `stokes_theorem_check` have returned for the cap, and how would you detect the mistake without knowing the answer?

Answer form: a table; a short table with an order; a derivation; two sentences.

## 6. (Judgment: when "flux in equals flux out" fails)

A control volume is drawn around the test section of a wind tunnel and the six face fluxes are measured with a pitot rake. For each situation, say whether the net *volume* flux should be zero, whether the net *mass* flux should be zero, which term of ∂ρ/∂t + ∇ · (ρ**v**) = 0 is responsible, and give one measurable sign:

(a) Steady run at 40 m/s.
(b) The first two seconds after the fan is switched on.
(c) Steady run at Mach 0.8 in a transonic tunnel.
(d) A steady run at 40 m/s with a 2 % leak through a badly sealed access hatch on one wall.

Answer form: four short paragraphs, each with a number or a measurable quantity.

## 7. (Hard, uses your code: volume and centroid from a surface alone)

CAD systems compute the volume of a solid from its closed surface mesh without ever integrating over the interior.

(a) Show from the divergence theorem that V = (1/3)∯ **r** · **n** dA and that the x coordinate of the centroid is x_c = (1/V)∯ (x²/2) n_x dA.
(b) Parametrise the ellipsoid (x/a)² + (y/b)² + (z/c)² = 1 with a = 2, b = 1, c = 0.5 m as **r**(s, t) = (a sin s cos t, b sin s sin t, c cos s) and compute V with `surface_flux` for n = 10, 20, 40, 80. Compare with (4/3)πabc and report the convergence order.
(c) Shift the ellipsoid to be centred at (1, 0, 0) m and compute x_c from the flux formula with n = 80.
(d) The surface of a real wing is a triangulated mesh of flat facets, not a smooth parametrisation. For a flat triangle with vertices **p**₁, **p**₂, **p**₃, write the exact contribution to (1/3)∯ **r** · **n** dA in terms of a scalar triple product, and say what happens to the computed volume if one triangle's vertex order is reversed.

Answer form: two derivations; a table with an order; one number; a formula and one sentence.
