# P08.08  Potential Flow

> Outside the thin layer where viscosity acts, the flow past a body is the gradient of a potential that obeys Laplace's equation. Its solutions add, so a few singularities build a body, and one vortex gives it lift ρUΓ and no drag at all.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the velocity potential and stream function of irrotational incompressible plane flow, the complex potential W(z) with dW/dz = u − iv, and the elementary flows (uniform stream, source, sink, doublet, vortex). 2. Superpose them into the Rankine oval and the cylinder with and without circulation, and derive the surface pressure, the stagnation points and the Kutta–Joukowski lift ρUΓ, by pressure integration and by Blasius's theorem. 3. Implement and test the superposition, contour integrals for circulation, flux and force, and a traced body surface, checking zero normal flow from the geometry, lift = ρUΓ and zero drag (d'Alembert) by two independent routes. 4. Use the Joukowski map and the Kutta condition to obtain c_l = 2π sin α for a flat plate, apply it to the light aeroplane's wing, and say which of its predictions survive in a real, viscous flow. |
| **Prerequisites** | `P08.02` (stream function, vorticity, irrotational flow), `P08.04` (Bernoulli), `P08.03` (momentum), `P02.12` (complex variables: analytic functions, Laurent series, residues, conformal maps), `P01.11` (line integrals), `P03.04` (quadrature: the periodic trapezoid rule), `P02.10` (the FFT, for spectral differentiation). Skill check: you can write the Cauchy–Riemann equations and evaluate ∮ dz/z round the origin. |
| **Engineering question** | The light aeroplane (750 kg, 10 m span, 1.5 m chord) cruises at 60 m/s at 3000 m: what circulation must its wing carry, at what angle of attack does a flat plate of its chord carry it, where does that put the stagnation points, and how much drag does the potential-flow model predict? The answers: Γ = 13.48 m²/s for the mean 735.5 N/m of lift (c_l = 0.300); α = 2.73° for the flat plate (Kutta condition); the stagnation point of the equivalent cylinder moves to −2.73°, exactly the plate's angle; and the drag is zero to 10⁻¹³ N/m by both routes, d'Alembert's paradox. |
| **Assumptions and limits** | Steady, incompressible, inviscid, irrotational plane flow (no boundary layer, no separation, no wake); bodies are closed streamlines; the Kutta condition is imposed, not derived. The model is accurate for the pressure and lift on streamlined bodies in attached flow at high Reynolds number; it predicts zero drag for every closed body, fails completely behind bluff bodies (where the real flow separates), and gives an infinite suction at a sharp leading edge. The 2-D strip treats the wing's lift as uniform along the span (no tip vortices; those are `P08.09`). |
| **Analytical work** | φ and ψ from irrotationality and continuity; Laplace's equation; the complex potential and velocity; the five elementary flows; the cylinder W = U(z + a²/z) plus a vortex; its surface speed, Cp = 1 − (2 sin θ + Γ/(2πaU))², stagnation points and lift by hand; Blasius's theorem and Kutta–Joukowski from a residue; the Rankine oval's length and thickness; the Joukowski map and the Kutta circulation 4πUa sin α. |
| **Computational deliverable** | `starter/main.py`: `complex_potential`, `complex_velocity`, `velocity`, `pressure_coefficient`, `cylinder`, `cylinder_surface_cp`, `cylinder_stagnation_points`, `rankine_oval`, `rankine_oval_size`, `oval_contour`, `contour_integral`, `circulation_and_flux`, `force_by_pressure`, `force_blasius`, `kutta_circulation`, `joukowski`, `zeta_from_z`, `plate_complex_velocity` (the element constructors, `circle_contour`, `stream_function_oval`, the spectral derivative `_dz_dtheta` and `main` are given). Tests prove: each element's velocity is the derivative of its potential (Cauchy–Riemann both ways) and its potentials are harmonic; source flux m and vortex circulation −Γ counterclockwise, zero for a contour that misses them; the doublet is the limit of a source–sink pair; superposition is linear; zero radial velocity on the cylinder with and without circulation and incidence; the closed-form Cp; velocity zero at the computed stagnation points; lift ρUΓ and drag 0 by pressure integration and by Blasius on three contours; the Rankine oval's stagnation points, thickness, zero normal flow on the traced surface from its own geometry, zero force and zero net flux; the Joukowski map and its inverse; the Kutta condition at the trailing edge, a finite trailing-edge and infinite leading-edge speed, and c_l = 2π sin α. |
| **Practical activity** | **computational experiment**: a circulation sweep for the lifting cylinder (stagnation points, lift, drag), a contour-independence and resolution study of the two force routes, the Rankine oval family, and c_l(α) with the leading-edge suction for the flat plate. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `potential_flow.png`, the hand derivations of the cylinder's surface pressure and lift and of Blasius's theorem, the four experiment tables, and a half-page note on the wing's circulation and what the potential-flow answer leaves out. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the elementary flows and the cylinder are derived; lift ρUΓ is obtained by hand and checked two ways in code; all public tests pass; the note states Γ, α and the zero drag, and explains which predictions survive viscosity. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapter "Potential Flow Theory": Definition of a Potential Flow; Velocity Potential & Stream Function; Elementary Flows; Superposition of Elementary Flows; Complex Variable Methods in Potential Flow (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/potential-flows/) (required-free; link only). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 6 "Potential Flows": Basic Flows; Superposition; DOUBLET; Uniform Flow over a Cylinder; Rankine Half-Body; Enclosed Bodies; Blasius Theorem and Lift Force for an Arbitrary Body (https://eng.libretexts.org/Bookshelves/Civil_Engineering/Intermediate_Fluid_Mechanics_(Liburdy)/06%3A_Potential_Flows) (required-free). `@mit-ocw-16-01-unified-engineering-fluids` 16.01–16.02 notes F12 "Stream Function and Velocity Potential", F15 "Sources, Sinks and Doublets", F16 "Vortex, Lifting Flow" (required-free: the aerospace treatment, short notes). `@libretexts-bar-meir-fluid-mechanics` Ch. 10 "Inviscid Flow or Potential Flow": 10.3.1 Flow Around a Circular Cylinder; 10.3.1.1 Adding Circulation to a Cylinder; 10.4.1 Complex Potential and Complex Velocity (supplementary-free). `@mit-ocw-2-25-advanced-fluid-mechanics` topic "Potential Flow Theory", problem "The Magnus Effect" with solution (https://ocw.mit.edu/courses/2-25-advanced-fluid-mechanics-fall-2013/pages/potential-flow-theory/) (supplementary-free). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics* (paid; library), the text the Unified Engineering notes follow. |
| **Extensions** | ★ Build a source panel method: N straight panels of constant source strength on the cylinder, solve for zero normal flow at the panel midpoints, and compare Cp with 1 − 4 sin²θ as N grows. ★ Offset the Joukowski circle's centre to make a cambered, thick aerofoil; impose the Kutta condition at its trailing edge and plot its Cp. ★ Use the method of images to put the lifting cylinder near a wall (ground effect) and find how its lift changes with height. |
| **Time** | unmeasured |

## 1. The engineering problem

A wing makes lift by turning the air downward, and the momentum balance of `P08.03` can say how much once the flow is known. But how is the flow round a wing to be found? The Navier–Stokes equations of `P08.06` have no exact solution for it. At the light aeroplane's Reynolds number, about 5 × 10⁶ on its chord, viscosity matters only in a boundary layer that is a small fraction of the chord thick, and outside it the air behaves as if it had no viscosity at all. If the flow also starts from rest, so that its particles carry no rotation, the velocity outside the boundary layer is the gradient of a potential, and the whole problem reduces to Laplace's equation: linear, with solutions that can be added.

That is the idea of this lesson. A handful of elementary solutions (a uniform stream, sources and sinks, doublets and vortices) superpose into the flow past an oval, a cylinder, and, through a conformal map, a flat-plate wing. One quantity, the **circulation** Γ round the body, sets the lift: the Kutta–Joukowski theorem L′ = ρUΓ. And the same theory says the drag is exactly zero, which is wrong for every real body and tells us precisely what viscosity contributes.

The light aeroplane (750 kg, span 10 m, chord 1.5 m) cruises at 60 m/s at 3000 m, where ρ = 0.90913 kg/m³ (the ISA of `P08.01`). Treat the lift as spread uniformly along the span, 735.5 N per metre. The questions: what Γ must the wing carry, at what angle does a flat plate carry it, where are the stagnation points, and what drag does the model give? Success criteria, written first: every elementary flow is checked to satisfy Laplace's equation and the Cauchy–Riemann equations; bodies have zero normal flow, checked from their geometry rather than from ψ; the force is computed two independent ways (surface pressure and Blasius's contour integral) and agrees with ρUΓ to round-off.

## 2. Understand and derive

### 2.1 Potential and stream function

In plane incompressible flow, continuity ∂u/∂x + ∂v/∂y = 0 is satisfied identically by a **stream function** ψ with u = ∂ψ/∂y, v = −∂ψ/∂x (`P08.02`). Lines of constant ψ are streamlines, and the flow between two streamlines is the difference of their ψ values. If the flow is also **irrotational**, ω = ∂v/∂x − ∂u/∂y = 0, there is a **velocity potential** φ with u = ∂φ/∂x, v = ∂φ/∂y. Substituting each into the other condition:

**∇²φ = 0,  ∇²ψ = 0.**   (1)

Both are harmonic, and u = φ_x = ψ_y, v = φ_y = −ψ_x are the **Cauchy–Riemann equations** of `P02.12`. So W = φ + iψ is an analytic function of z = x + iy, the **complex potential**, and its derivative is the **complex velocity**

**dW/dz = u − iv.**   (2)

Every analytic function is a possible flow (Leishman, "Potential Flow Theory", Complex Variable Methods in Potential Flow; Bar-Meir 10.4.1). Laplace's equation and the wall condition (no flow through the body: the body is a streamline) are linear, so solutions **superpose**. Bernoulli's equation of `P08.04`, p + ½ρ|V|² = p∞ + ½ρU², then gives the pressure, and it is *not* linear: the pressure must always be computed from the total velocity, never added.

Why should a real flow be irrotational? Kelvin's theorem (`P08.09`) says that in an inviscid fluid the circulation round a loop of particles never changes. Air that starts at rest, or comes from a uniform stream far upstream, has no vorticity and acquires none outside the boundary layer.

### 2.2 The elementary flows

| flow | W(z) | dW/dz | what it is |
|---|---|---|---|
| uniform | U e^{−iα} z | U e^{−iα} | speed U at angle α |
| source at z₀ | (m/2π) log(z − z₀) | m/(2π(z − z₀)) | radial outflow, volume flux m per unit depth (m < 0: sink) |
| doublet at z₀ | κ e^{iβ}/(2π(z − z₀)) | −κ e^{iβ}/(2π(z − z₀)²) | a source and sink merged, axis β |
| vortex at z₀ | (iΓ/2π) log(z − z₀) | iΓ/(2π(z − z₀)) | circular flow, circulation Γ clockwise |

(Leishman, "Elementary Flows"; Liburdy Ch. 6, Basic Flows and DOUBLET.) Two of them carry a number that any closed contour round them measures. The contour integral ∮(dW/dz) dz = ∮(u dx + v dy) + i∮(u dy − v dx) has as its real part the counterclockwise **circulation** and as its imaginary part the outward **volume flux**. By the residue theorem, a source gives flux m and no circulation; a vortex gives circulation and no flux. With the aerodynamic convention used here (Γ positive **clockwise**, so that lifting flow from left to right has Γ > 0), the counterclockwise integral of a vortex is −Γ. The doublet is the limit of a source m at z₀ − ε e^{iβ} and a sink −m at z₀ + ε e^{iβ} as ε → 0 with 2mε = κ fixed: the tests check the limit to O(ε²).

### 2.3 The cylinder

Uniform flow plus a doublet of strength κ = 2πUa² gives

**W = U(z + a²/z).**   (3)

On |z| = a, with z = ae^{iθ}, a²/z = ae^{−iθ} and W = 2Ua cos θ is real, so ψ = 0: the circle is a streamline, and (3) is the flow past a cylinder of radius a (Liburdy, Uniform Flow over a Cylinder). Adding a vortex (iΓ/2π) log z adds (Γ/2π) ln r to ψ, which is constant on the circle, so the circle stays a streamline for any Γ. On the surface the velocity is purely tangential:

**u_θ = −2U sin θ − Γ/(2πa),   Cp = 1 − (2 sin θ + Γ/(2πaU))².**   (4)

Without circulation, Cp = 1 − 4 sin²θ: stagnation (Cp = 1) at the front and back, a suction peak Cp = −3 (speed 2U) at the shoulders, and a pattern symmetric fore and aft and top and bottom. A symmetric pressure gives no force at all. With clockwise Γ the flow speeds up over the top and slows underneath, the **stagnation points move down to sin θ_s = −Γ/(4πUa)**, and they merge at the bottom when Γ = 4πUa. Above that, a single stagnation point leaves the body and sits in the flow below it.

**The lift by hand.** The force per unit span is F = −∮ p n ds with p = p∞ + ½ρU²Cp. Expanding (4), Cp = 1 − 4 sin²θ − (4Γ/(2πaU)) sin θ − (Γ/(2πaU))². Against n = (cos θ, sin θ), every term integrates to zero over 0..2π except the one in sin θ against sin θ:

L′ = −∫₀^{2π} ½ρU² (−(4Γ/(2πaU)) sin θ) sin θ a dθ = ½ρU² a (4Γ/(2πaU)) π = **ρUΓ,  D = 0.**   (5)

### 2.4 Blasius's theorem and Kutta–Joukowski

With n ds = −i dz on a counterclockwise contour, F_x + iF_y = i∮ p dz. Using Bernoulli, conjugating, and using the fact that dW is real along the body (a streamline), the pressure integral becomes **Blasius's theorem**:

**F_x − iF_y = (iρ/2) ∮ (dW/dz)² dz,**   (6)

and because (dW/dz)² is analytic outside the body, the contour can be moved to any curve enclosing it (Liburdy, Blasius Theorem and Lift Force for an Arbitrary Body). Far away, the velocity of any closed body (no net source) has the Laurent series dW/dz = U + iΓ/(2πz) + O(z⁻²). Squaring, the 1/z coefficient is iUΓ/π, and the residue theorem gives ∮ = −2UΓ, so

**F_x = 0,  F_y = ρUΓ**   (Kutta–Joukowski)   (7)

for **any** body shape, not just the cylinder: the lift depends only on the circulation, and the drag is zero. That is **d'Alembert's paradox**: a steady potential flow exerts no drag on any closed body. Real bodies have drag because of viscosity (skin friction, and separation behind bluff bodies, `P08.16`), not because the inviscid argument is wrong.

The code computes (5) and (6) separately. `force_by_pressure` integrates Bernoulli's pressure round the body surface; `force_blasius` integrates (dW/dz)² round any enclosing contour. Both use the periodic trapezoid rule, with dz/dθ by FFT differentiation of the contour points. For a smooth closed contour the periodic trapezoid rule converges faster than any power of the spacing (`P03.04`). On a circle, Cp is a trigonometric polynomial of degree 2, and even four points integrate it exactly.

### 2.5 The Rankine oval

A source m at (−a, 0) and an equal sink at (+a, 0) in a uniform stream U make a closed body, because the net source strength is zero: everything the source emits, the sink swallows. Its stagnation points, where U + (m/2π)(1/(x + a) − 1/(x − a)) = 0, are at **x = ±l, l = a√(1 + m/(πUa))**. The body is the dividing streamline ψ = 0 with ψ = Uy + (m/2π)(θ₁ − θ₂), and its half-thickness t solves **t = (m/(πU)) arctan(a/t)** (Leishman, "Superposition of Elementary Flows"; Liburdy, Enclosed Bodies). The angles θ₁, θ₂ from the source and sink must be taken in a branch that keeps ψ continuous on the body, which the principal branch of `np.angle` does once the sign conventions are right; `stream_function_oval` is given. For m/(πUa) = 2: l = 1.732a, t = 1.307a, fineness ratio 1.33.

### 2.6 A wing: the Joukowski map and the Kutta condition

The map **z = ζ + a²/ζ** sends the circle |ζ| = a to the segment −2a ≤ x ≤ 2a: a flat plate of chord c = 4a. Far away z ≈ ζ, so the free stream is the same in both planes. A flow past the circle in the ζ-plane (free stream at angle α, circulation Γ) maps to a flow past the plate. dW/dz = (dW/dζ)/(dz/dζ) with dz/dζ = 1 − a²/ζ², which vanishes at ζ = ±a, the images of the plate's edges. The velocity there is infinite unless dW/dζ also vanishes there. Choosing Γ so that the rear stagnation point sits at ζ = a makes the flow leave the **trailing edge** smoothly with finite speed. That choice is the **Kutta condition**:

**Γ = 4πUa sin α = πcU sin α,   c_l = L′/(½ρU²c) = 2π sin α.**   (8)

Viscosity enforces the Kutta condition in a real flow: a flow that went round a sharp trailing edge would separate there at once, and the starting vortex of `P08.09` is what carries off the difference. At the **leading edge** the speed stays infinite (the plate has no thickness), and the pressure there produces a finite **leading-edge suction** force along the plate. It exactly cancels the drag component of the normal force: in plate (body) axes the force is tilted, in wind axes it is pure lift (Experiment 4).

**The wing.** L′ = 7355.0/10 = 735.50 N/m and **Γ = L′/(ρV) = 13.484 m²/s**; c_l = L′/(½ρV²c) = **0.2996**; a flat plate needs **α = asin(c_l/2π) = 2.733°**. The equivalent cylinder of radius a = c/4 = 0.375 m with that Γ has Γ/(4πUa) = 0.0477 and stagnation points at **−2.733° and −177.27°**: the rear one sits at θ = −α. That is the Kutta condition seen in the cylinder plane with the stream along x. Real thin aerofoils in attached flow have a lift slope close to 2π per radian, shifted by camber, so the circulation is right; the angle depends on the section's camber.

## 3. Build the computational model

Open `starter/main.py`. The `Element` dataclass and its constructors `uniform`, `source`, `doublet`, `vortex`, and `circle_contour`, `stream_function_oval`, `_dz_dtheta` and `main` are given.

### Step 1: superposition

`complex_potential(elements, z)` and `complex_velocity(elements, z)` loop over the elements and add each one's W or dW/dz from the table in 2.2, differentiated by hand (not numerically: the tests compare them with finite differences of W). `velocity` returns (u, v) = (Re, −Im) of dW/dz; `pressure_coefficient` is 1 − |dW/dz|²/U².

### Step 2: bodies

`cylinder(U, a, Gamma, alpha)` returns [uniform(U, α), doublet(2πUa², 0, α), vortex(Γ)]; `cylinder_surface_cp` is (4); `cylinder_stagnation_points` returns the two surface points, or for Γ > 4πUa the single point below the body at r = (Γ + √(Γ² − (4πUa)²))/(4πU). `rankine_oval` and `rankine_oval_size` are 2.5 (the thickness by `brentq`). `oval_contour(U, m, a, n)` traces the surface: for n equally spaced polar angles about the centre, it finds the radius where `stream_function_oval` changes sign, by `brentq` along the ray (the oval is star-shaped about its centre).

### Step 3: contour integrals

`contour_integral(values, contour)` is Σ f(z_k) z′(θ_k) Δθ with z′ from `_dz_dtheta`:

```python
def contour_integral(values, contour):
    return complex(np.sum(values * _dz_dtheta(contour)) * 2 * np.pi / contour.size)
```

`circulation_and_flux` applies it to dW/dz; `force_by_pressure` to the Bernoulli pressure (F_x + iF_y = i∮p dz); `force_blasius` to (dW/dz)² with (6).

### Step 4: the plate

`kutta_circulation` and `joukowski` are (8) and the map. `zeta_from_z` inverts the map, choosing the root of ζ² − zζ + a² = 0 with |ζ| ≥ a (the exterior of the circle). `plate_complex_velocity` divides the cylinder's dW/dζ by dz/dζ.

Running `main()` prints

```
light aeroplane: weight 7355.0 N, span 10.0 m, chord 1.5 m, cruise 60.0 m/s at 3000 m (rho = 0.90913 kg/m^3)
  mean lift per span 735.50 N/m -> circulation Gamma = L'/(rho V) = 13.484 m^2/s, c_l = 0.2996
  flat plate of the same chord (Kutta): alpha = asin(c_l / 2 pi) = 2.733 deg
  equivalent cylinder a = c/4 = 0.375 m: Gamma / (4 pi U a) = 0.04769; stagnation points at [  -2.733 -177.267] deg
  pressure integration on the surface: drag -5.58e-15 N/m, lift 735.4988 N/m
  Blasius on r = 3a:                    drag -1.32e-13 N/m, lift 735.4988 N/m; rho U Gamma = 735.4987
  Kutta circulation for that alpha 13.484 m^2/s; Blasius in the physical (plate) plane, body axes: F_x = -35.075 N/m (leading-edge suction), F_y = 734.662 N/m
    wind axes: lift 735.4987 N/m, drag 1.4e-13 N/m
  circulation round r = 2a (counterclockwise) -13.4836 m^2/s = -Gamma; flux -3.5e-15 m^2/s

Rankine oval m/(pi U a) = 2.0: half-length 1.73205 a, half-thickness 1.30654 a, fineness 1.3257
  max |normal velocity| on the traced surface 2.1e-13 U; force -2.2e-17, 0.0e+00 (rho U^2 a units): d'Alembert
  max speed on the oval 1.73881 U, min Cp -2.02346
```

and writes `potential_flow.png`: the Rankine oval's streamlines and traced surface, the streamlines of a lifting cylinder, and the surface Cp for Γ = 0, the cruise Γ and Γ = 2πUa.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-08-potential-flow`. The tests:

1. **Elements**: dW/dz equals the derivative of W taken along x and along y (analyticity); u = φ_x = ψ_y, v = −ψ_x; the five-point Laplacian of W vanishes.
2. **Integrals**: a source's flux m and zero circulation; a vortex's counterclockwise circulation −Γ and zero flux; both zero on a contour that misses them; the clockwise sense of a positive vortex.
3. **Doublet** as the limit of a source–sink pair, error falling as ε².
4. **Superposition**: total velocity = sum of parts; no elements, no flow.
5. **Cylinder**: zero radial velocity on the surface for three (Γ, α) pairs, uniform flow far away; the closed-form Cp equals Bernoulli's; extremes −3 and 1; zero velocity at every computed stagnation point, all on the lower half, merged at −ia for Γ = 4πUa, one point off the body above that.
6. **Forces**: pressure integration and Blasius on two other contours (one off-centre) agree, give lift ρUΓ and no drag, for three Γ; Cp depends only on Γ/(Ua).
7. **Rankine oval**: zero velocity at ±l; ψ(0, t) = 0; thin and thick limits; the traced points have ψ = 0 and, from their own FFT-computed tangents, zero normal velocity to 10⁻⁸; zero force; zero net flux.
8. **Plate**: the map sends the circle to [−2a, 2a] and inverts; the Kutta circulation puts the stagnation point at ζ = a; finite speed at the trailing edge and a speed growing as d^{−1/2} at the leading edge; lift ρUΓ and zero drag in wind axes from Blasius in the plate plane; c_l = 2π sin α at 1°, 4° and 8°.
9. **main** writes the figure.

### Compare

The two force routes are independent: one integrates a pressure on the body with Bernoulli, the other a squared velocity on a contour that can be anywhere outside the body, and neither knows the answer ρUΓ. They agree with it to 10⁻¹² for the cruise wing. The third route is the residue calculation (7), done by hand. The established tool for bodies of arbitrary shape is the panel method (Extension 1), which replaces the analytic singularities by a numerical distribution on the surface. The tests run it in miniature for the oval: the surface is found numerically (by root-finding on ψ) and its zero normal flow is then checked from its own geometry. That is the check a panel code makes at its control points.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Circulation sweep.** The cruise-speed cylinder (U = 60 m/s, a = 0.375 m) with Γ/(4πUa) = 0, 0.25, 0.5, 1, 1.5. Predict the stagnation angles and the lift. *Finding:* stagnation at (0°, 180°), (−14.48°, −165.52°), (−30°, −150°), (−90°, −90°) merged, and for 1.5 a single point off the body at r = 2.618a on the negative y axis; lift by pressure integration 0, 3855.757, 7711.514, 15 423.03, 23 134.54 N/m, equal to ρUΓ to all printed figures; drag below 6 × 10⁻¹³ N/m in every case. Lift is linear in Γ, whatever the stagnation points do.
2. **Contour independence and resolution.** Blasius on circles of radius 1.01a, 2a, 10a and an off-centre circle of 3a; pressure integration on 4–32 surface points. Predict whether the result depends on the contour or on n. *Finding:* lift errors of 10⁻¹² to 2 × 10⁻¹² N/m (round-off) on every contour; pressure integration exact to 10⁻¹² N/m with as few as four points, because Cp on a circle is a trigonometric polynomial of degree 2 and the periodic trapezoid rule integrates it exactly. For the oval, the force cancels by symmetry even on a 16-point contour (2 × 10⁻¹⁶), so the stronger oval test is the zero normal flow on the traced surface (2 × 10⁻¹³ U on 256 points).
3. **The oval family.** m/(πUa) from 0.1 to 100. Predict the fineness ratio and the peak surface speed. *Finding:* l/t falls from 7.34 to 1.007 and the peak speed rises from 1.29U to 1.99U (the cylinder's 2U) (Problem 5 has the table); the force is zero to 10⁻¹⁵ for every one.
4. **The flat plate.** c_l and the leading-edge suction for α = 1°, 2.733°, 5° and 10°. Predict c_l. *Finding:* c_l = 0.10966, 0.29959, 0.54762, 1.09106, equal to 2π sin α to five figures; the leading-edge suction force is 4.70, 35.07, 117.2 and 465.1 N/m, which is L′ sin α in every case: exactly what is needed to turn the plate-normal force into a force perpendicular to the stream. Real plates cannot sustain that suction at a sharp edge; the flow separates there at a few degrees, which is why aerofoils have rounded noses.

Raw data format: CSV tables `Gamma_ratio, theta_s1_deg, theta_s2_deg, r_s_over_a, lift_pressure, lift_rhoUGamma, drag` (1); `contour, radius_over_a, n, lift_error, drag` (2); `m_ratio, l_a, t_a, fineness, vmax_U, Fx, Fy` (3); `alpha_deg, cl, cl_2pi_sin, LE_suction_N_m` (4); plus `potential_flow.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| lift of the cruise cylinder | ρUΓ = 735.49875 N/m | pressure and Blasius 735.49875 | 2 × 10⁻¹³ N/m | round-off (the print rounds the last digit differently) |
| drag | 0 (d'Alembert) | −5.6 × 10⁻¹⁵, −1.3 × 10⁻¹³ N/m | round-off | — |
| stagnation point | −α = −2.733° | −2.733° | 0 | Kutta condition in the cylinder plane |
| circulation round 2a | −Γ (counterclockwise) | −13.4836 m²/s | — | residue theorem |
| plate c_l at 2.733° | 2π sin α = 0.29959 | 0.29959 (Blasius, plate plane) | < 10⁻⁵ | — |
| oval stagnation and thickness | l = √3 a, t from (2.5) | 1.73205a, 1.30654a; ψ = 0 there | < 10⁻¹² | — |
| oval normal flow | 0 | 2.1 × 10⁻¹³ U | — | FFT tangent of a smooth curve |

The model reproduces itself exactly, as it should: these are the checks that the code solves the potential-flow problem it claims to. The engineering answer: **the wing must carry Γ = 13.5 m²/s (c_l = 0.30); a flat plate of its chord does so at 2.73° with the rear stagnation point at the trailing edge; and the predicted drag is zero.** The validation against reality is the part potential flow cannot do by itself. The lift and the pressure outside the boundary layer survive in attached flow; the drag does not (it is all viscous: skin friction from `P08.10` and `P08.14`, and induced drag from the tip vortices of `P08.09`, which a 2-D strip leaves out); and the leading-edge suction peak of a sharp plate is not physical.

What to improve: thickness and camber (Extension 2), the finite span (`P08.09`, and the lifting-line theory that follows), and the boundary layer for drag. What not to improve: the resolution of the contour integrals, which are already exact to round-off.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (4), (5) and (6)–(7); `potential_flow.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the wing's circulation, its flat-plate angle of attack, the zero drag and what it leaves out, and which predictions you would use in design.

### Adversarial self-checks

- **What if the speed doubles at the same lift?** Γ = L′/(ρV) halves to 6.74 m²/s and c_l falls by four, to 0.075, so α falls to 0.68°; the drag stays zero in the model. In reality the skin-friction drag would roughly quadruple (it scales with V²).
- **What if Γ goes to zero?** The flow is symmetric, Cp = 1 − 4 sin²θ, and there is no force at all, by symmetry and by (7).
- **Which assumption fails first in the real world?** Attached flow: behind a cylinder the boundary layer separates and the measured drag is large, so the cylinder results are a mathematical device for the plate, not a model of a real cylinder. For the wing, the 2-D strip: the real lift falls to zero at the tips, which changes the local Γ and adds induced drag.
- **What measurement could falsify the model?** A pressure-tapped aerofoil in a wind tunnel: the measured Cp outside the boundary layer should follow the potential-flow prediction for the same lift; a plateau of constant Cp near the trailing edge means separation, and the model has stopped applying there.
- **Does the result survive a change of method?** Lift by surface pressure, by Blasius on four different contours and by the residue theorem agree to 10⁻¹²; c_l from Blasius in the physical plane equals 2π sin α.
- **Which decision changes if the uncertainty doubles?** If the lift slope were uncertain by 10 % instead of 5 %, the cruise angle of attack would be known to ±0.27° instead of ±0.14°. That matters for setting the wing's incidence on the fuselage (and so the cabin's attitude in cruise), not for the circulation the wing must carry, which is fixed by the weight.

## Check yourself

1. What is the complex potential, and how are u and v obtained from it?
2. Why can elementary flows be superposed, and why can pressures not be?
3. Which elementary flows make a cylinder, and what doublet strength does it need?
4. Where are the stagnation points of a lifting cylinder, and what happens when Γ exceeds 4πUa?
5. State the Kutta–Joukowski theorem with its sign convention.
6. State d'Alembert's paradox and say what potential flow leaves out.
7. What does the Kutta condition fix for a flat plate, and what lift coefficient results?
8. Why is a body surface a line of constant ψ, and how did the tests check zero normal flow without using ψ?
