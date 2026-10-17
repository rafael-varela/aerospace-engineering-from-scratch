# P10.05  Vortex Panel Method

> Cover the airfoil with straight panels carrying a vortex sheet whose strength varies linearly, make every panel a streamline, add the Kutta condition, and a linear system gives the flow round any shape. Then prove it converges.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the velocity induced by a straight panel with a linearly varying vortex sheet, and assemble the (N + 1) × (N + 1) system of flow tangency at the panel midpoints plus the Kutta condition γ₀ + γ_N = 0. 2. Implement and test the method from scratch, with c_l from the circulation and from the surface pressure, the pressure drag and c_m,c/4. 3. Verify it against the exact Karman–Trefftz solution of `P10.03` with observed orders of convergence for c_l and Cp, and against thin-airfoil theory (`P10.04`) as the thickness goes to zero. 4. Decide how many panels a given accuracy needs, and where the method's own answer is not to be trusted. |
| **Prerequisites** | `P01.13` (linear systems), `P03.08` (direct solvers), `P10.03` (Kutta condition, exact Karman–Trefftz flow), `P10.04` (thin-airfoil theory), `P10.02` (NACA coordinates in Selig order), `P03.11` (observed order, Richardson extrapolation), `P01.11` (line integrals). Skill check: you can integrate 1/(t² + η²) and t/(t² + η²), and rotate a vector between two frames. |
| **Engineering question** | What are c_l and c_m,c/4 of the trainer's NACA 2412 section at α = 2° by an inviscid panel method, how many panels does 0.1 % accuracy in c_l need, and how far from thin-airfoil theory is the cruise attitude? Answers: c_l = 0.5009 (Richardson-extrapolated), c_m,c/4 = −0.0583; 80 panels give 0.08 % and 160 give 0.014 %; the cruise c_l of 0.331 comes at α = 0.59° against thin-airfoil theory's 0.94°, because thickness raises the lift slope to 1.099 × 2π. Against the exact Karman–Trefftz flow, c_l converges at an observed order of 1.7–2.2 and Cp (5–95 % of chord) at order 1.9–2.0. |
| **Assumptions and limits** | Steady, 2-D, inviscid, incompressible potential flow; a closed polygon with straight panels; a linear vortex distribution on each; flow tangency at panel midpoints; the fluid inside the body at rest, so the surface speed equals the sheet strength; the Kutta condition at a sharp trailing edge. No drag (the computed pressure drag measures discretisation error); no boundary layer, so no viscous loss of lift and no stall; the corner flow at a finite-angle trailing edge is not resolved by any panel. |
| **Analytical work** | The four influence integrals I₁–I₄ of a linear vortex panel; the velocity in panel and global axes; the self-influence on the outer side; counting unknowns and equations; the Kutta condition; c_l from Γ and from Cp; the circle with a "Kutta node" (Problem 3). |
| **Computational deliverable** | `starter/main.py`: `panel_geometry`, `linear_vortex_influence`, `node_influence`, `vortex_panel`, `panel_forces`, `kt_convergence` (the NACA generator from `P10.02`, the exact Karman–Trefftz helpers from `P10.03`, `thin_airfoil_naca4` from `P10.04`, `observed_order`, `richardson` and `main` are given). Tests prove: unit outward normals; the influence equals a direct quadrature of point vortices to 10⁻¹⁰ and a point vortex from far away; flow tangency at every midpoint and the Kutta condition to 10⁻¹⁰; zero lift and symmetric Cp for a 0012 at α = 0; superposition of two free streams to 10⁻¹⁰; c_l from Cp within 0.5 % of c_l from Γ and pressure drag falling with N; against Karman–Trefftz, c_l converging at mean order 1.6–2.4 and Cp monotonically, below 10⁻³ at N = 320; thin-section lift approaching 2πα and thin-airfoil c_l and c_m from above. |
| **Practical activity** | **computational experiment**: panel convergence with observed order against the exact Karman–Trefftz solution; Richardson extrapolation on the NACA 2412; the thin-airfoil limit as thickness → 0; the trailing-edge corner. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `vortex_panel.png`, the derivation of the influence integrals and the system, the convergence tables with observed orders, and a half-page note on how many panels you would use and why. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: influence integrals derived; all public tests pass; observed orders for c_l and Cp reported against the exact solution, with the trailing-edge behaviour explained; the thin-airfoil limit shown. |
| **References** | `@mit-ocw-16-unified-fluids` 16.01 Fluids (fall) Lecture 17 "Flowfield Prediction" (Flowfield prediction; Source Sheets) and Lecture 18 "Lift Prediction, Vortex Sheets" (Prediction of Lift; Vortex Sheets) (https://ocw.mit.edu/courses/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/resources/f17_fall/, …/f18_fall/) (required-free). `@ntrs-nasa-tp-2995-panel-methods` Erickson, *Panel Methods: An Introduction*: 2. What Panel Methods Can and Cannot Do; 3. Common Roots of Panel Methods; 5. Modeling (https://ntrs.nasa.gov/citations/19910009745) (required-free). `@leishman-intro-aerospace-flight-vehicles` "Potential Flow Theory": General Method of Singularities; Superposition of Elementary Flows (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/potential-flows/) (supplementary-free, link only). `@drela-xfoil` XFOIL User Primer "Inviscid Formulation" (https://web.mit.edu/drela/Public/web/xfoil/xfoil_doc.txt) (supplementary-free: a production linear-vorticity panel method). `@ntrs-naca-tr-452-arbitrary-wing-sections` (supplementary-free: the conformal-mapping alternative). Optional textbooks for deeper study: J. Katz and A. Plotkin, *Low-Speed Aerodynamics* (paid; library), the standard step-by-step derivation of 2-D and 3-D panel methods including this one. |
| **Extensions** | ★ Rewrite the method in XFOIL's stream-function form (ψ = constant at the nodes) and compare its convergence. ★ Add a constant-strength source panel distribution and model the boundary layer's displacement by transpiration (`P10.06`). ★ Make a two-element version (a slat or slotted flap, `P10.08`) with a Kutta condition on each element. |
| **Time** | unmeasured |

## 1. The engineering problem

Thin-airfoil theory (`P10.04`) ignores thickness; the conformal maps of `P10.03` give exact answers only for special shapes. A designer needs the flow round *any* section: a NACA 2412, a laminar-flow section, a slotted flap. The tool that made this routine, from the 1960s on, is the **panel method**: cover the surface with panels carrying singularities (sources, doublets or vortices), make the surface a streamline at one point per panel, and solve a linear system. Every production airfoil code in use (XFOIL's inviscid core among them) is a refinement of this idea.

A panel method is a numerical approximation, so its answers have an error, and you must know how big. This lesson builds a linear-strength vortex panel method from scratch, runs it on the exact Karman–Trefftz section of `P10.03` to measure its order of convergence, and on thinner and thinner NACA sections to recover thin-airfoil theory. Then it answers the trainer's question: what are its section's c_l and c_m, and how many panels are enough?

Success criteria, written first: flow tangency and the Kutta condition hold to round-off; a symmetric section at zero angle has zero lift; c_l from the circulation and from the pressure agree and the pressure drag tends to zero; against the exact solution, c_l converges at an observed order near 2 and Cp away from the trailing edge converges too; as the thickness goes to zero, the panel c_l tends to thin-airfoil theory.

## 2. Understand and derive

### 2.1 The model

The airfoil surface, listed in Selig order (trailing edge → upper surface → leading edge → lower surface → trailing edge, counter-clockwise), is cut into N straight panels between N + 1 nodes. On panel j a **vortex sheet** of strength γ(s) varies linearly from γ_j at its start to γ_{j+1} at its end. γ is positive clockwise (the Phase 10 convention, so that Γ > 0 and L′ = ρVΓ > 0 for lift). The node values γ₀ … γ_N are the unknowns. The trailing edge appears twice (nodes 0 and N), so there are N + 1 unknowns even for a closed edge.

Two facts make the method work. The sheet makes the velocity jump across the surface by γ; and if the flow inside the body is at rest, the velocity just outside equals the jump. So **the surface speed is |γ|**, and Cp = 1 − (γ/V)².

### 2.2 Velocity induced by one panel

In panel coordinates (ξ along the panel from its start, η to its left), a clockwise element γ ds at (s, 0) induces du = γη/(2πD), dv = −γ(ξ − s)/(2πD) at (ξ, η), with D = (ξ − s)² + η². With r₁, r₂ the distances to the panel's ends and β the angle the panel subtends at the point,

**I₁ = ∫₀^L η/D ds = β,  I₂ = ∫₀^L (ξ − s)/D ds = ln(r₁/r₂),**
**I₃ = ∫₀^L sη/D ds = η ln(r₂/r₁) + ξβ,  I₄ = ∫₀^L s(ξ − s)/D ds = −L + ηβ + ξ ln(r₁/r₂).**   (1)

With γ(s) = g_a(1 − s/L) + g_b s/L,

**u = [g_a(I₁ − I₃/L) + g_b I₃/L]/2π,  v = −[g_a(I₂ − I₄/L) + g_b I₄/L]/2π,**   (2)

rotated back to global axes. These are the panel's **influence coefficients** for its two end nodes. On the panel itself (η → 0) β jumps from +π to −π: the outside of a counter-clockwise surface is to the right of the panel (η < 0), so the self-influence uses β = −π. That is the velocity jump of §2.1 appearing in the algebra. A test checks (2) against direct numerical quadrature of point vortices along the panel; far away the panel looks like a point vortex of strength L(g_a + g_b)/2.

**Worked example.** A panel from (0, 0) to (1, 0) with uniform γ = 1 induces at (0.5, 0.5) (where it subtends a right angle) u = β/2π = 0.25, v = 0; at (2, 0), u = 0, v = −ln 2/2π = −0.1103. A point vortex of strength 1 at the midpoint would give 0.318 and −0.106: close to a panel, its length matters.

### 2.3 The linear system and the Kutta condition

At each panel midpoint the normal velocity must vanish:

**Σ_k A_ik γ_k = −V(cos α n_x,i + sin α n_y,i),  i = 0 … N − 1,**   (3)

where A_ik collects the normal components of (2) from the two panels that share node k. That is N equations for N + 1 unknowns. The missing one is the **Kutta condition**:

**γ₀ + γ_N = 0.**   (4)

γ₀ and γ_N are minus the surface velocities along the node direction at the two trailing-edge nodes, which point forward on the upper surface and aft on the lower; equal and opposite γ means equal speeds leaving the edge on both sides, the flow leaving the edge smoothly. Without (4) the system has a one-parameter family of solutions, the circulation, exactly as in `P10.03`; Problem 3 shows that on a circle the same row produces Γ = 4πVa sin α, the Kutta circulation of a "trailing edge" placed at the rear point.

### 2.4 Forces

The circulation is Γ = Σ L_j(γ_j + γ_{j+1})/2 (exact for linear γ), so **c_l,Γ = 2Γ/(Vc)**. Independently, the pressure on each panel, −Cp n L per unit q, summed and rotated to wind axes gives c_l,p and c_d,p; the moment about the quarter chord gives c_m,c/4. In exact inviscid flow c_d = 0, so c_d,p measures the discretisation error. The chord runs from the trailing-edge node to the node farthest from it.

### 2.5 The benchmarks

**Karman–Trefftz.** Put the nodes on the exact surface of a Karman–Trefftz section at equal steps of the circle angle (which clusters them at both edges). The exact c_l, c_m and Cp at each midpoint's circle angle are known (`P10.03`, copied here as given helpers). The test section has b = 1, μ = −0.08 + 0.06i, τ = 10°: exact c_l = 0.868141 and c_m,c/4 = −0.096619 at α = 4°.

**Thin-airfoil theory.** For NACA sections with 2 % camber at 40 % and thickness t → 0, the panel c_l and c_m must tend to `P10.04`'s c_l = 0.66644 and c_m,c/4 = −0.05312 at α = 4°.

### 2.6 Where the model stops

- **The trailing-edge corner.** At a finite-angle edge the exact speed falls to zero over a region smaller than any panel, as distance to the power 2π/(2π − τ) − 1 (0.03 for τ = 10°). The panel next to the edge always has an O(1) Cp error; it does not spoil c_l, but it does stall the c_m error at about 1.5 × 10⁻⁵ (§4).
- **Viscosity.** The panel method's lift slope (1.099 × 2π for the 2412) is the inviscid one; the real section's is lower (`P10.03`, Problem 6). `P10.06` adds the boundary layer.
- **Geometry.** A polygon is an O(h²) approximation of the surface; with cosine spacing that is the method's order too.

## 3. Build the computational model

Open `starter/main.py`. The NACA generator, the exact Karman–Trefftz helpers, `thin_airfoil_naca4`, `observed_order`, `richardson` and `main` are given.

### Step 1: geometry

`panel_geometry(x, y)` returns midpoints, lengths, unit tangents along the node order and outward normals n = (t_y, −t_x).

### Step 2: influence

`linear_vortex_influence(px, py, geo, outside_self)` evaluates (1) and (2) for every point and every panel at once with NumPy broadcasting, returning the (a) and (b) parts in global axes:

```python
xi = rx * tx + ry * ty;  eta = -rx * ty + ry * tx
beta = np.arctan2(eta, xi - L) - np.arctan2(eta, xi)
lnr = 0.5 * np.log((xi**2 + eta**2) / ((xi - L)**2 + eta**2))
I3, I4 = -eta * lnr + xi * beta, -L + eta * beta + xi * lnr
ua, ub = (beta - I3 / L) / (2*np.pi), I3 / L / (2*np.pi)
va, vb = -(lnr - I4 / L) / (2*np.pi), -I4 / L / (2*np.pi)
```

With `outside_self=True` the diagonal uses η = −0 and β = −π. `node_influence` adds the (a) part of panel k and the (b) part of panel k − 1 into the column of node k.

### Step 3: solve

`vortex_panel(x, y, alpha)` builds the normal-velocity rows (3), the Kutta row (4), calls `np.linalg.solve`, and returns γ, the surface velocity −(γ_j + γ_{j+1})/2 at each midpoint, Cp, Γ, the residual and the coefficients from `panel_forces`.

### Step 4: the convergence study

`kt_convergence(b, mu, tau, alpha, ns)` runs the method on the exact surface for each N and returns the errors in c_l (both ways), c_m and Cp (on 5–95 % of the chord and over the whole surface).

Running `main()` prints

```
Karman-Trefftz test section (b = 1, mu = (-0.08+0.06j), tau = 10 deg), alpha = 4 deg: exact cl = 0.868141, cm_c/4 = -0.096619
     N    |cl_G - cl|   |cl_p - cl|   |cm - cm_ex|  max|dCp| 5-95%   max|dCp| all
    20    1.053e-02    1.352e-02    3.817e-04    1.522e-02      2.635e-01
    40    2.922e-03    4.229e-03    2.494e-04    8.144e-03      5.025e-01
    80    8.574e-04    1.182e-03    1.917e-04    4.281e-03      5.641e-01
   160    1.997e-04    2.790e-04    1.403e-05    1.309e-03      5.462e-01
   320    6.161e-05    8.098e-05    1.758e-05    3.287e-04      5.102e-01
   640    1.365e-05    1.838e-05    1.462e-05    8.582e-05      4.729e-01
  observed order, cl from Gamma: [1.85 1.77 2.1  1.7  2.17]
  observed order, cl from Cp   : [1.68 1.84 2.08 1.78 2.14]
  observed order, Cp 5-95 %    : [0.9  0.93 1.71 1.99 1.94]
NACA 2412 (closed TE), alpha = 2 deg: cl_Gamma by N = 40: 0.49921, 80: 0.50051, 160: 0.50084, 320: 0.50088, 640: 0.50090
  observed order from successive differences: [1.97 2.93 1.09]; Richardson (p = 2) cl = 0.50091
  N = 160: lift slope 6.9037 /rad = 1.0988 x 2 pi; alpha_L0 = -2.154 deg; cruise cl 0.3309 at alpha = 0.592 deg (thin airfoil 0.940 deg)
```

(abridged) and writes `vortex_panel.png`: panel Cp against the exact Karman–Trefftz Cp, the convergence plot with a slope-2 line, and the thin-airfoil limit.

### Verify

`python -m pytest curriculum/phase-10/lesson-05-vortex-panel-method` runs 13 tests (listed in the metadata).

### Compare

Two exact references (Karman–Trefftz for thick, cambered, finite-angle sections; thin-airfoil theory for the limit t → 0), and a third route within the method (c_l from Γ against c_l from Cp). XFOIL's inviscid mode (optional, `P10.21`) solves the same problem with a stream-function formulation.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment: panel convergence with observed order, against the exact Karman–Trefftz solution.** Predict before you run.

1. **Convergence against the exact solution.** Predict the order for c_l. *Finding:* the table above. c_l from Γ: 1.05 × 10⁻² (N = 20) to 1.37 × 10⁻⁵ (N = 640), orders 1.85, 1.77, 2.10, 1.70, 2.17; from Cp: orders 1.68–2.14. Cp on 5–95 % of the chord: 1.5 × 10⁻² to 8.6 × 10⁻⁵, orders rising from 0.9 to 1.99 and 1.94 once the leading edge is resolved. Over the whole surface the Cp error stays near 0.5 (the trailing-edge panel, §2.6), and the c_m error stalls at about 1.5 × 10⁻⁵ from N = 160.
2. **How many panels for the trainer?** NACA 2412 at α = 2°. Predict N for 0.1 %. *Finding:* c_l = 0.49921, 0.50051, 0.50084, 0.50088, 0.50090 for N = 40–640; Richardson (p = 2) 0.50091. Errors 0.339 %, 0.080 %, 0.014 %, 0.005 %, 0.001 %: **N = 80 meets 0.1 %**, N = 160 gives 0.014 %. c_m,c/4 = −0.0584, −0.0583, −0.0583, −0.0583, −0.0583. (The successive-difference orders, 1.97, 2.93, 1.09, wander because cosine node positions shift relative to the curvature change at p; the extrapolated value is safe to 10⁻⁵.)
3. **Thin-airfoil limit.** 2 % camber at 40 %, α = 4°, N = 320. Predict the lift excess at 12 % thickness. *Finding:* c_l = 0.74156, 0.70305, 0.68413, 0.67477, 0.67011 for t/c = 0.12, 0.06, 0.03, 0.015, 0.0075, against 0.66644: ratios 1.1127, 1.0549, 1.0265, 1.0125, 1.0055. c_m,c/4 = −0.0612, −0.0569, −0.0550, −0.0541, −0.0535 against −0.0531. The excess per unit thickness, (c_l − c_l,thin)/t, falls from 0.63 to 0.49: nearly, not exactly, linear in t (a thin round nose is not a thin-airfoil nose).
4. **The trainer's attitude.** *Finding:* lift slope 6.904 /rad = 1.099 × 2π, α_L0 = −2.154°: the cruise c_l of 0.331 at α = 0.59°, against 0.94° from thin-airfoil theory.

Raw data format: CSV `N, err_clG, err_clp, err_cm, err_cp_window, err_cp_all` (1); `N, cl, cm` (2); `t_over_c, cl, cl_thin, cm, cm_thin` (3).

## 5. Validate and improve

| check | reference | panel method | difference |
|---|---|---|---|
| flow tangency, Kutta | 0 | residual < 10⁻¹⁰ | round-off |
| symmetric section, α = 0 | c_l = 0 | < 10⁻¹² | round-off |
| c_l, Karman–Trefftz, N = 320 | 0.868141 (exact) | error 6.2 × 10⁻⁵ | order ≈ 2 |
| Cp, 5–95 % chord, N = 640 | exact | 8.6 × 10⁻⁵ | order ≈ 2 |
| NACA 2412 c_l, N = 160 | Richardson 0.50091 | 0.50084 | 0.014 % |
| thin limit, t = 0.75 % | thin airfoil 0.66644 | 0.67011 | +0.55 %, ∝ t |
| circle with a rear node | 4πVa sin α | N = 160: 1.09516 vs 1.09523 | 7 × 10⁻⁵ |

The engineering answer: **the trainer's NACA 2412 has c_l = 0.5009 and c_m,c/4 = −0.0583 at α = 2° in inviscid flow; 80 panels give 0.1 % accuracy and 160 give 0.014 %; the cruise c_l of 0.331 is reached at α = 0.59°, 0.35° below thin-airfoil theory, because 12 % thickness raises the lift slope by 10 %.** The numerical error is now far below the modelling error: the boundary layer (`P10.06`) and the finite wing (`P10.09`) each change the answer by percent, not hundredths of a percent.

What to improve: the physics (the boundary layer), not the panel count. What not to trust: the Cp of the last panel at the trailing edge.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1)–(4); `vortex_panel.png`; the four experiment tables with observed orders; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on how many panels you would use for the trainer and why.

### Adversarial self-checks

- **What if α changes?** The system matrix does not; only the right-hand side does, and the solution is exactly cos α times the α = 0 solution plus sin α times the α = 90° one (a test checks it to 10⁻¹⁰). Lift is linear in sin and cos of α: no stall, ever.
- **What if the thickness → 0?** The panels on the two surfaces approach each other and the system becomes ill-conditioned; at t/c = 0.75 % it still gives 0.55 % above thin-airfoil theory with 320 panels.
- **Which assumption fails first?** The inviscid, attached flow: at the trailing edge the boundary layer is thick enough to change the circulation by several per cent.
- **What could falsify the code?** A symmetric section with lift at α = 0; a c_l that changes with the starting node of the list; an observed order far from 2 on the exact benchmark.
- **Does the result survive a change of N?** Yes, with a known error and order; Richardson extrapolation gives the N → ∞ value.
- **Which decision changes if the panel error were 1 %?** None for the trainer, whose c_l is uncertain by more than that from viscosity; but a stability derivative or a hinge moment computed from small Cp differences would need the finer grid.

## Check yourself

1. How many unknowns and how many equations does a linear-vortex panel method with N panels have, and what closes the system?
2. What does the Kutta condition γ₀ + γ_N = 0 express physically?
3. Name two ways to get c_l from a panel solution. What should their difference do as N grows?
4. What observed order of convergence do you expect for c_l against an exact solution, and what did the Karman–Trefftz study show?
5. Why does the Cp error over the whole surface stop falling, while the error away from the trailing edge keeps falling?
6. What happens to the panel c_l as the thickness goes to zero, and what theory does it approach?
7. What exact solution is used to verify the panel method, and how are the panel nodes placed on it?
