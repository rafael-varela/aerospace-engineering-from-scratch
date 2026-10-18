# P10.03  The Kutta Condition & Circulation

> Potential flow allows any circulation round an airfoil. The sharp trailing edge chooses one, and that circulation is the lift: L′ = ρVΓ.

| | |
|---|---|
| **Learning outcomes** | 1. Write the complex potential of a uniform stream past a circular cylinder with circulation and derive its surface velocity, pressure coefficient and stagnation points. 2. Map the circle to Joukowski and Karman–Trefftz airfoils, and derive the Kutta condition Γ = 4πVa sin(α + β) from the requirement that the trailing-edge velocity be finite. 3. Implement and test the exact airfoil flow and verify L′ = ρVΓ four independent ways: surface-pressure integration, Blasius' theorem, a momentum balance on a control surface, and the circulation line integral round contours of any shape. 4. Explain why exact potential flow over-predicts the measured lift slope, and what the Kutta condition stands in for. |
| **Prerequisites** | `P08.08` (potential flow, superposition of elementary flows; recapped in §2.1), `P02.12` (complex variables and conformal maps, the Joukowski map), `P01.11` (line integrals and circulation), `P01.12` (Stokes' theorem: the same circulation on every enclosing contour), `P03.04` (the trapezoid rule and its accuracy on periodic functions). Skill check: you can differentiate a complex function and evaluate a line integral round a circle. |
| **Engineering question** | What circulation must each metre of the reference trainer's wing carry at cruise, and at what angle of attack does an exact potential-flow section with the NACA 2412's thickness (12 %), camber (2 %) and trailing-edge angle (16°) produce it? Answers: Γ = 10.92 m²/s (c_l = 0.331, uniform loading assumed); α = 0.38° from the chord line, with a zero-lift angle of −2.38° and a lift slope of 1.094 × 2π. The lift from four independent routes agrees to 10⁻⁷ and the drag is zero, as it must be. |
| **Assumptions and limits** | Steady, two-dimensional, inviscid, incompressible, irrotational flow outside the body; a sharp trailing edge; the Kutta condition stands in for the viscous shedding that sets the circulation. Real sections have boundary layers: the measured NACA 0012 lift slope (Ladson, Problem 6) is 0.94–0.97 × 2π, 12–17 % below the exact 1.094 × 2π. No drag is predicted (d'Alembert). The Karman–Trefftz family is not the NACA family: same thickness and camber, different distribution. |
| **Analytical work** | The complex potential and velocity of the cylinder with circulation; u_θ = −2V sin(θ − α) − Γ/(2πa); stagnation angles; the cylinder's lift by direct pressure integration; the Joukowski and Karman–Trefftz maps and their derivatives; the Kutta condition; c_l = 8πa sin(α + β)/c and its thin-section limit 2π sin α; Blasius' theorem and the momentum balance. |
| **Computational deliverable** | `starter/main.py`: `circle_dwdz`, `circle_surface_speed`, `circle_surface_cp`, `stagnation_angles`, `kutta_circulation`, `circle_geometry`, `karman_trefftz`, `karman_trefftz_dzdzeta`, `airfoil_flow`, `pressure_force_coefficients`, `blasius_force`, `momentum_force`, `circulation_line_integral` (`karman_trefftz_inverse`, `contour_in_z`, `airfoil_velocity_at`, `section_shape`, `fit_karman_trefftz` and `main` are given). Tests prove: the circle is a streamline and the far field is the stream; stagnation points have zero speed and leave the body at Γ > 4πaV; the cylinder's pressure lift is ρVΓ with zero drag; Karman–Trefftz with τ = 0 is Joukowski; the map derivative and inverse are right; the Kutta circulation makes the trailing-edge point a stagnation point and keeps the edge velocity finite; pressure lift equals Kutta–Joukowski at three angles; Blasius, momentum and the line integral give Γ and ρVΓ on a circle, an ellipse and a rectangle; a contour that misses the body gives zero; thickness raises the slope as 1 + 0.77 t/c. |
| **Practical activity** | **computational experiment**: lift four ways on three contours; the trailing-edge speed with and without the Kutta condition; the convergence of pressure integration for cusped and finite-angle trailing edges; the thickness effect on lift slope; the exact slope against authentic NACA 0012 data. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `kutta_circulation.png`, the derivations of §2.2–2.5, the experiment tables, and a half-page note on why the exact lift slope differs from the measured one. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the Kutta condition derived from the finite-velocity requirement; all public tests pass; lift shown equal to ρVΓ by at least three routes with the agreement quantified; the note names the viscous effect behind the measured slope. |
| **References** | `@mit-ocw-16-100-aerodynamics` Lecture 13 "Kutta Condition", Lecture 14 "Kutta Condition Solution", study material "Lift Generation and Streamline Curvature" (https://ocw.mit.edu/courses/16-100-aerodynamics-fall-2005/resources/16100lectre13_cg/, …/16100lectre14_cg/, …/4liftgen/) (required-free). `@mit-ocw-16-01-unified-engineering-fluids` 16.01 Fluids Lecture 16 "Vortex, Lifting Flow" (Vortex; Lifting flow about circular cylinder) and spring Lecture 1 "Formation of Lifting Flow" (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f16_fall/, …/f01_sp/) (required-free). `@leishman-intro-aerospace-flight-vehicles` "Classic Airfoil Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/classic-airfoil-theory/): Kutta Condition; Starting Vortex & Kelvin's Theorem; Conformal Transformations & Joukowsky Airfoils; and "Potential Flow Theory" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/potential-flows/): Understanding Circulation and Vorticity; Complex Variable Methods in Potential Flow (required-free, link only). `@nasa-glenn-beginners-guide-aeronautics` "Lift of a Rotating Cylinder", "Conformal Mapping" (https://www.grc.nasa.gov/www/k-12/airplane/cyl.html, …/map.html) (supplementary-free). `@ntrs-naca-tr-452-arbitrary-wing-sections` (supplementary-free, optional depth: mapping an arbitrary section to a circle). Data: `@ntrs-nasa-tm-4074-ladson-naca0012` via `curriculum/phase-10/data/` (Problem 6). Optional textbooks for deeper study: J. D. Anderson, *Fundamentals of Aerodynamics* (paid; library), the standard treatment of the Kutta–Joukowski theorem and the cylinder with circulation. |
| **Extensions** | ★ Compute the pitching moment about the quarter chord from Blasius' moment theorem and compare with thin-airfoil theory (`P10.04`). ★ Use the Theodorsen–Garrick idea (TR-452) to map a NACA 2412 to a near-circle iteratively and compare its exact Cp with `P10.05`. ★ Animate the starting vortex: a point vortex −Γ shed at the trailing edge and convected downstream, and the bound circulation that Kelvin's theorem requires. |
| **Time** | unmeasured |

## 1. The engineering problem

Why does a wing lift? The flow-turning and pressure-difference explanations are both true and both incomplete: they describe the lifting flow without saying why the flow chooses that pattern. The complete answer for a steady, two-dimensional, high-Reynolds-number flow fits in one line. The lift per unit span is **L′ = ρV∞Γ**, where Γ is the circulation, the line integral of velocity round the airfoil; and the circulation is set by the **Kutta condition**: the flow must leave the sharp trailing edge smoothly instead of whipping round it.

The reference trainer at cruise (750 kg, 55 m/s, span 10 m) supports 735.75 N per metre of span, so each metre of its wing must carry Γ = L′/(ρV) = 10.92 m²/s, a section lift coefficient of 0.331 on the 1.2 m chord (ignoring, until `P10.09`, that the loading of a finite wing is not uniform). This lesson builds the exact flows that show how an airfoil produces that circulation, using the conformal maps of `P02.12`, and checks the Kutta–Joukowski theorem until there is no room left to doubt it.

Success criteria, written first: on the exact flow, the pressure-integrated lift equals ρVΓ to the accuracy of the quadrature and the pressure drag is zero; Blasius' theorem, a momentum balance on a control surface and the circulation line integral give the same Γ on contours of three different shapes; a contour not enclosing the body gives zero; the trailing-edge velocity is finite only with the Kutta circulation.

## 2. Understand and derive

### 2.1 Recap: potential flow and complex potentials (`P08.08`)

In steady, incompressible, irrotational 2-D flow the velocity is the gradient of a potential φ that satisfies Laplace's equation, and the stream function ψ does too. Together they form the **complex potential** w(z) = φ + iψ, an analytic function of z = x + iy, whose derivative is the **complex velocity** dw/dz = u − iv. Elementary flows superpose: a uniform stream at angle α, w = Ve^{−iα}z; a doublet, w = μ/z; a point vortex, w = (iΓ/2π) ln z, which with this sign turns **clockwise** for Γ > 0. Bernoulli gives the pressure, Cp = 1 − |u|²/V².

### 2.2 The cylinder with circulation

A stream plus a doublet of the right strength makes the circle |ζ′| = a a streamline; adding a clockwise vortex keeps it one. With the circle centred at μ and ζ′ = ζ − μ:

**w(ζ) = V(ζ′e^{−iα} + a²e^{iα}/ζ′) + (iΓ/2π) ln ζ′.**   (1)

On the surface ζ′ = ae^{iθ}, (dw/dζ)e^{iθ} = u_r − iu_θ gives u_r = 0 and

**u_θ = −2V sin(θ − α) − Γ/(2πa),  Cp = 1 − (u_θ/V)².**   (2)

The stagnation points satisfy sin(θ − α) = −Γ/(4πaV): with no circulation they are fore and aft; clockwise circulation slides both down the back, they merge at the bottom when Γ = 4πaV, and beyond that a single stagnation point leaves the body. The flow is faster over the top, so the pressure is lower there. Integrating (2) round the circle (Problem 2d) gives **L′ = ρVΓ and D′ = 0** exactly. `main()` confirms it numerically: at Γ/(4πaV) = 0.5 the pressure-integrated c_l on the diameter is 6.283185 = 4πg, drag 3 × 10⁻¹⁶.

Equation (1) holds for every Γ. Potential theory cannot pick one.

### 2.3 Mapping the circle to an airfoil

An analytic map z = f(ζ) carries the circle flow to the flow past the image of the circle, because an analytic function of an analytic function is analytic. The **Joukowski map** z = ζ + b²/ζ sends a circle that passes through the critical point ζ = b (where dz/dζ = 0) and encloses ζ = −b to an airfoil with a cusped trailing edge at z = 2b. Moving the centre left (−μ_x) makes it thicker; moving it up (+μ_y) cambers it. The **Karman–Trefftz map** generalises it to a finite trailing-edge angle τ, with n = 2 − τ/π:

**(z − nb)/(z + nb) = ((ζ − b)/(ζ + b))ⁿ,  z = nb(1 + rⁿ)/(1 − rⁿ),  r = (ζ − b)/(ζ + b),**   (3)

**dz/dζ = 4n²b² r^{n−1} / ((1 − rⁿ)²(ζ + b)²).**   (4)

For n = 2 this is Joukowski (a test checks it). Outside the circle r never lies on the negative real axis, so the principal power is analytic there.

The airfoil-plane velocity is the circle-plane velocity divided by the map's derivative:

**u − iv = (dw/dζ)/(dz/dζ).**   (5)

![Left: the ζ-plane, with a circle of radius a centred at μ, left of and above the origin, passing through ζ = b and enclosing ζ = −b; the radius to ζ = b lies at angle β below the horizontal. Right: its image under the Karman–Trefftz map z = f(ζ), a cambered airfoil whose trailing edge, of angle τ, is the image of ζ = b at z = nb.](figures/karman-trefftz-map.svg)
*Figure 1. The map (3) carries the circle to the airfoil and the critical point ζ = b to the trailing edge z = nb. Moving μ left thickens the section, moving it up cambers it. Computed from (3) for b = 1, μ = −0.10 + 0.14i, τ = 10°.*

### 2.4 The Kutta condition

At ζ = b, (4) is zero: the map folds the smooth circle into a sharp edge. By (5) the airfoil velocity there is infinite (the flow whips round the trailing edge) unless dw/dζ also vanishes at ζ = b. The point ζ = b is on the circle at θ = −β, where b − μ = ae^{−iβ}. Setting u_θ(−β) = 0 in (2):

**Γ = 4πVa sin(α + β),  c_l = ρVΓ/(½ρV²c) = 8πa sin(α + β)/c.**   (6)

This is the **Kutta condition**: of all the flows (1), the one that leaves the sharp edge smoothly. It gives the zero-lift angle α_L0 = −β (set by camber) and a lift slope 8πa/c (set by thickness; a thin symmetric section has c = 4a and slope 2π). For a cusped edge the finite trailing-edge velocity is the same on both surfaces; for τ > 0 the trailing edge is a stagnation point.

**Why it is the right choice.** A real airfoil starting from rest sheds a vortex from its trailing edge, by viscous action in the boundary layer there, until the flow leaves the edge smoothly. Kelvin's theorem then requires an equal and opposite circulation round the airfoil. The Kutta condition is the end state of that viscous process written as a rule for the inviscid flow (`@leishman-intro-aerospace-flight-vehicles`, "Classic Airfoil Theory": Starting Vortex & Kelvin's Theorem).

**Worked example.** Fitting a Karman–Trefftz section to the 2412's t = 0.12, m = 0.02 and τ = 16° gives (b = 1) μ = −0.04850 + 0.04387i, a = 1.04942, β = 2.3957°, chord 3.8330, chord line inclined 0.0172° nose-up. Then α_L0 = −2.3957° + 0.0172° = −2.378° from the chord line, the slope is 8πa cos(α + β)/c = 6.873 /rad = 1.094 × 2π, and c_l = 0.3309 needs sin(α_x + β) = 0.3309 × 3.8330/(8π × 1.04942), α_x = 0.361°, α = 0.378° from the chord line. (The section's thickness peaks at 38 % chord and its camber at 50 %, not at 30 % and 40 % as on the 2412: same numbers, different shape.)

### 2.5 Four ways to measure lift

1. **Surface pressure.** The force is −∮p n ds; with n ds = (dy, −dx) on a counter-clockwise surface, c_x = −∮Cp dy/c and c_y = ∮Cp dx/c, then rotate by α. Parameterised by the circle angle, the trapezoid rule converges spectrally for a cusp and algebraically for a finite angle (Problem 5).
2. **Blasius' theorem.** X − iY = (iρ/2)∮(dw/dz)² dz on any contour enclosing the body: the residue of (dw/dz)² at infinity gives −iρVΓe^{−iα}.
3. **Momentum balance.** On any control surface C, F = −∮(p n + ρu(u·n)) ds: the force on the body equals the momentum flux and pressure on C.
4. **Circulation.** Γ = −Re∮(dw/dz) dz counter-clockwise (the minus for clockwise-positive Γ). In irrotational flow, by Stokes' theorem (`P01.12`), every contour enclosing the body gives the same value, and one that does not enclose it gives zero.

Routes 2–4 use contours away from the body, built in the z-plane (a circle, an ellipse and a rectangle; the given `karman_trefftz_inverse` maps their points back to ζ).

### 2.6 Where the model stops

- **Viscosity.** The boundary layer thickens towards the trailing edge and decambers the section: measured lift slopes are below the exact ones (Problem 6: 0.938 × 2π at Re = 2.0 × 10⁶ and 0.974 × 2π at 5.97 × 10⁶ for the NACA 0012, against 1.094 × 2π). At high angles it separates and the Kutta condition loses its meaning (stall, `P10.15`).
- **Drag.** Zero, always. Pressure drag and skin friction need the boundary layer (`P10.06`).
- **Three dimensions.** A finite wing sheds trailing vortices and its loading is not uniform (`P10.09`).
- **Shape family.** Karman–Trefftz sections are not NACA sections; `P10.05` handles any shape, and uses this exact solution as its benchmark.

## 3. Build the computational model

Open `starter/main.py`. The inverse map, the contours, the off-body velocity, the shape measurement and the fit are given.

### Step 1: the circle flow

`circle_dwdz(zeta, V, alpha, a, Gamma, mu)` is the derivative of (1). `circle_surface_speed` and `circle_surface_cp` are (2); `stagnation_angles` solves u_θ = 0 and returns an empty array if |Γ| > 4πaV.

### Step 2: maps and the Kutta condition

`circle_geometry(b, mu)` returns a = |b − μ| and β = −arg(b − μ). `karman_trefftz` and `karman_trefftz_dzdzeta` are (3) and (4). `kutta_circulation` is (6).

### Step 3: the exact airfoil flow

`airfoil_flow(b, mu, tau, alpha, V, n, Gamma=None)` places n points at circle angles −β + 2π(k + ½)/n (never exactly on the edge, where (5) is 0/0), maps them, and returns z, dz/dθ = (dz/dζ) i a e^{iθ}, the complex velocity (5), the speed, Cp, Γ (Kutta by default) and the chord.

### Step 4: forces and circulation

```python
def pressure_force_coefficients(z_dtheta, cp, alpha, chord, dtheta):
    dx, dy = np.real(z_dtheta) * dtheta, np.imag(z_dtheta) * dtheta
    cx, cy = -np.sum(cp * dy) / chord, np.sum(cp * dx) / chord
    return cy * np.cos(alpha) - cx * np.sin(alpha), cx * np.cos(alpha) + cy * np.sin(alpha)
```

`blasius_force`, `momentum_force` and `circulation_line_integral` take the complex velocity at contour nodes and the weighted increments dz from `contour_in_z`.

Running `main()` prints (abridged)

```
trainer cruise: lift per metre span W/b = 735.75 N/m -> Gamma = L'/(rho V) = 10.920 m^2/s, cl = 2 Gamma/(V c) = 0.3309
Karman-Trefftz section (tau = 16 deg): mu = -0.04850 + 0.04387 i (b = 1), a = 1.04942, beta = 2.3957 deg; chord 3.8330, t = 0.1200 at 0.383, m = 0.0200 at 0.504
zero-lift angle (chord line) -2.378 deg; cruise cl 0.3309 at alpha = 0.378 deg (chord line); lift slope 6.8730 /rad = 1.0939 x 2 pi
exact flow, n = 512: cl from Kutta-Joukowski 0.3309158374, from surface pressure 0.3309159267, cd from pressure -1.0e-06; min cp -0.5562
  circle    contour: Gamma = 0.63420348 (Kutta 0.63420348); Blasius L' = 0.63420348, D' = 5.6e-16; momentum L' = 0.63420348, D' = 3.7e-15
  rectangle contour: Gamma = 0.63420348 (Kutta 0.63420348); Blasius L' = 0.63420348, D' = -1.9e-15; momentum L' = 0.63420348, D' = -2.2e-15
speed at the surface point nearest the trailing edge (n = 4096): with Kutta 0.5203 V, with Gamma = 0 32.1 V
```

and writes `kutta_circulation.png`: cylinder Cp for three circulations, streamlines of the Kutta flow round the section, and its exact surface Cp.

### Verify

`python -m pytest curriculum/phase-10/lesson-03-kutta-condition-circulation` runs 18 tests (listed in the metadata). The tolerances are set by the method: 10⁻¹² for closed forms, 10⁻⁶ for surface pressure with 1024 points on a 10° trailing edge, 10⁻⁹ for contour integrals with Gauss or periodic trapezoid nodes.

### Compare

The independent routes are the point of this lesson: four estimates of the lift from different integrals, one of them over the body and three over contours of different shapes, all agreeing.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment**, with one comparison against authentic data. Predict before you run.

1. **Lift four ways.** Predict whether Blasius on a rectangle will match Kutta–Joukowski as closely as on a circle. *Finding:* Γ, Blasius L′ and momentum L′ all equal 0.63420348 (V = 1, ρ = 1) on the circle (radius 1.5c), the ellipse and the rectangle, with drags of 10⁻¹⁵; a circle of radius 1 placed beside the body gives Γ = 0 (a test). Surface pressure with n = 512 gives c_l = 0.3309159 against 0.3309158 and c_d = −1.0 × 10⁻⁶.
2. **What the Kutta condition does.** Predict the speed next to the trailing edge with Γ = 0. *Finding:* 32.1V at the nearest of 4096 points, growing without bound as the point approaches the edge; with the Kutta circulation, 0.52V there, falling (very slowly, as distance^{0.05}) to zero at a 16° edge.
3. **Convergence of pressure integration.** Predict which trailing edge converges faster. *Finding:* (Problem 5) the cusped Joukowski section reaches round-off by n = 128; the 16° section's error falls 2.2 × 10⁻³ → 5.4 × 10⁻⁷ from n = 16 to 256, irregularly, limited by the branch point of (4) at the edge.
4. **Thickness and lift slope.** Predict the slope of a 20 % Joukowski section. *Finding:* 1.0308, 1.0616, 1.0926, 1.1236, 1.1549 × 2π for t/c = 0.04–0.20, against 1 + 0.77 t/c = 1.0308 … 1.1540: thickness adds lift slope almost linearly.
5. **Against measurement.** Predict whether the exact slope will be above or below the authentic NACA 0012 data. *Finding:* (Problem 6) Ladson, Table I, free transition: 0.938 × 2π at Re = 2.0 × 10⁶ and 0.974 × 2π at 5.97 × 10⁶, against 1.094 × 2π for a 12 % symmetric Karman–Trefftz section. The difference shrinks with Re: it is the boundary layer.

Raw data format: CSV `contour, Gamma, L_blasius, D_blasius, L_momentum, D_momentum` (1); `n, err_tau0, err_tau16` (3); `t_over_c, slope_over_2pi` (4); `Re, slope_measured_over_2pi` (5), citing the source file.

## 5. Validate and improve

| quantity | exact | numerical route | difference |
|---|---|---|---|
| cylinder lift, g = 0.5 | ρVΓ: c_l = 6.283185 | surface pressure, 512 points | < 10⁻⁶; drag 3 × 10⁻¹⁶ |
| section c_l at cruise | 0.3309158 (Kutta–Joukowski) | surface pressure, n = 512 | 9 × 10⁻⁸ |
| Γ on three contours | 0.63420348 | line integral | < 10⁻⁹ |
| L′ on three contours | 0.63420348 | Blasius; momentum | < 10⁻⁹; drag < 4 × 10⁻¹⁵ |
| lift slope, 12 % Joukowski | 1 + 0.77t/c = 1.0924 | 1.0926 | 2 × 10⁻⁴ (the rule is a fit) |
| lift slope vs NACA 0012 data | 1.094 × 2π | measured 0.938–0.974 × 2π | 12–17 %: viscosity |

The engineering answer: **each metre of the trainer's wing carries Γ = 10.92 m²/s at cruise; an exact section with the 2412's thickness, camber and trailing-edge angle produces it at α = 0.38° (zero-lift angle −2.38°, lift slope 1.094 × 2π), and every route to its lift agrees to 10⁻⁷.** The real section will need more angle than that: its measured slope is lower (Problem 6), and the finite wing needs more again (`P10.09`). The mathematics is not the uncertainty; the boundary layer is.

What to improve: add the boundary layer's displacement effect (`P10.06` does, and finds a 2 % loss of c_l on the 2412). What not to improve: the number of points, which already makes the quadrature error negligible.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (6) and the cylinder's lift; `kutta_circulation.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on why the exact lift slope is above the measured one.

### Adversarial self-checks

- **What if V → 0?** Γ ∝ V and L′ ∝ V²: no flow, no lift. Kutta–Joukowski is a statement about a steady flow at speed.
- **What if α + β = 0?** Γ = 0 and the section is at zero lift with both stagnation points at the edges of the circle's diameter through b; c_l changes sign with α + β.
- **Which assumption fails first?** The thin boundary layer at the trailing edge: at high α it separates, the flow no longer leaves the edge smoothly, and (6) over-predicts lift badly (stall, `P10.15`).
- **What measurement could falsify the theorem?** A wake survey in a 2-D tunnel gives the circulation from the velocity round a contour; a balance gives the lift. Their ratio should be ρV within the uncertainty. A momentum-integral check on PIV data does exactly this.
- **Does the result survive a change of contour or method?** Yes, by construction and by test: three contours, three integrals, one surface integral.
- **Which decision changes if the lift slope is 15 % lower?** The trainer's cruise angle rises by about 15 % of (α − α_L0), roughly 0.4°, which matters for the fuselage attitude and the cabin floor angle, and the design must carry that margin.

## Check yourself

1. Define the circulation Γ and state the sign convention used in this phase.
2. State the Kutta–Joukowski theorem. What drag does it predict?
3. Why does potential flow need the Kutta condition, and what does the condition require?
4. Write the Kutta circulation for a Joukowski or Karman–Trefftz section and say which point it makes a stagnation point.
5. Why would the velocity at a sharp trailing edge be infinite without the right circulation?
6. What property of an irrotational flow makes the circulation the same on every contour round the airfoil? What is it on a contour that does not enclose the airfoil?
7. How does thickness change the lift slope of a Joukowski section? Give the slope ratio to 2π at 4 %, 12 % and 20 % thickness.
