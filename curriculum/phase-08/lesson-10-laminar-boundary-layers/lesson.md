# P08.10  Boundary Layers: Laminar

> Viscosity acts only in a thin layer next to the wall; inside it, one ordinary differential equation describes every station of a flat plate, and its wall slope f″(0) = 0.332 is the skin friction.

| | |
|---|---|
| **Learning outcomes** | 1. Derive Prandtl's boundary-layer equations from the Navier–Stokes equations of `P08.06` by an order-of-magnitude argument, and state the scaling δ/x ~ Re_x^(−1/2) and why the pressure is constant across the layer. 2. Reduce them for a flat plate to the Blasius equation f‴ + ½ff″ = 0 with a similarity variable, and derive the displacement and momentum thicknesses, the shape factor and the momentum integral dθ/dx = c_f/2. 3. Implement and test the Blasius solution by RK4 shooting with a secant root finder, check it against `scipy.integrate.solve_bvp` and against a marching finite-difference solution of the boundary-layer PDEs that never assumes similarity, and measure the convergence of all three. 4. Compute δ99, δ*, θ, H, c_f and friction drag for a real wing, and decide where the laminar flat-plate model stops being valid. |
| **Prerequisites** | `P08.06` (Navier–Stokes equations, no slip, Stokes's first problem), `P08.05` (Reynolds number), `P08.01` (ISA and Sutherland's law), `P08.03` (control-volume momentum balance), `P03.05` (RK4), `P03.02` (root finding: secant method), `P02.02` (second-order ODEs and boundary conditions), `P03.11` (finite differences, observed order), `P03.08` (banded solvers). Skill check: you can turn a third-order ODE into a first-order system and integrate it with RK4. |
| **Engineering question** | If the boundary layer on the light aeroplane's wing stayed laminar from the leading edge to the 1.5 m trailing edge at 3000 m and 60 m/s, how thick would it be and how much friction drag would the wing's 30 m² of wetted surface carry? Answer: δ99 = 3.35 mm at the trailing edge (δ* = 1.17 mm, θ = 0.45 mm), C_F = 6.04 × 10⁻⁴ and 29.7 N of friction drag; this is the floor that transition (`P08.12`) and a turbulent layer (`P08.14`, 166 N if turbulent from the nose) raise. |
| **Assumptions and limits** | Steady, two-dimensional, incompressible flow of a Newtonian fluid; laminar; layer thin compared with the distance from the leading edge (Re_x ≫ 1, so the approximation fails within about 1 mm of the wing's leading edge); zero pressure gradient (a flat plate at zero incidence; an aerofoil's pressure gradient is `P08.11`); smooth wall; no transition (it happens, `P08.12`). |
| **Analytical work** | Order-of-magnitude reduction of Navier–Stokes to the boundary-layer equations; the similarity transformation; δ*, θ, H; the momentum integral equation and θ = 2f″(0)x/√Re_x; local and average skin friction (C_F = 2c_f(L)); Töpfer's scaling, which gives f″(0) from one integration (Problem 2). |
| **Computational deliverable** | `starter/main.py`: `blasius_rhs`, `rk4_blasius`, `shoot_blasius`, `blasius_bvp`, `blasius_integrals`, `plate_boundary_layer`, `friction_drag`, `march_boundary_layer`, `wall_cf_from_profile` (`isa_troposphere`, `observed_order`, `_cumtrapz` and `main` are given). Tests prove: RK4 is fourth order; the shooting solution meets all three boundary conditions, satisfies the Blasius equation when substituted back, and agrees with `solve_bvp` to 10⁻⁷ and with the published f″(0) = 0.332057 to 2 × 10⁻⁶; δ99, δ*, H and θ = 2f″(0) are right; the layer scales as √x and is unchanged by a change of units; the drag equals the integrated wall shear; and the marched PDE conserves mass, converges at first order in x and collapses onto the Blasius profile. |
| **Practical activity** | **computational experiment**: RK4 step convergence, the effect of truncating the domain at η_max, the `solve_bvp` tolerance, a grid study of the marching solver and its similarity collapse, and a speed and altitude sweep of the wing's laminar friction drag. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `blasius.png`, the derivation of the Blasius equation and of θ = 2f″(0)x/√Re_x, the five experiment tables, and a half-page note stating the wing's laminar boundary-layer thickness and friction drag with the assumptions that limit them. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the boundary-layer equations and the Blasius equation are derived with every assumption named; all public tests pass; the RK4 order and the domain-truncation error are measured and separated; the note gives δ99, C_F and the drag at the trailing edge and says where the laminar model fails first. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` "Boundary Layer Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/): Laminar & Turbulent Boundary Layers; Development of a Boundary Layer; Boundary Layer Equations; Laminar Boundary Layer – Blasius Solution (required-free; link only, CC BY-NC-ND). `@mit-ocw-16-100-aerodynamics` Lecture 29 "Laminar Boundary Layer Order of Magnitude Analysis" and Lecture 30 "Solutions of the Laminar Boundary Layer Equations" (Flat Plate (Laminar): Blasius Solution) (required-free). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 9 "Boundary Layer Flows": Physical Attributes of a Boundary Layer; Boundary Layer Equations; Laminar Boundary Layer Solution (required-free). `@mit-ocw-2-25-advanced-fluid-mechanics` McKinley, "The Generalized Boundary Layer Equations" (supplementary-free: the same reduction at graduate depth). `@nasa-glenn-beginners-guide-aeronautics` "Boundary Layer" (https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/boundary-layer/) (supplementary-free: the picture before the mathematics). Optional textbooks for deeper study: Schlichting and Gersten, *Boundary-Layer Theory*, and White, *Viscous Fluid Flow* (paid; a library copy is enough; see `references/phase-08.md`). |
| **Extensions** | ★ Use Töpfer's transformation (Problem 2b) to solve Blasius with one integration and compare its cost with shooting. ★ Add a second-order (Crank–Nicolson in x) version of `march_boundary_layer` and show its observed order is 2. ★ Add a constant suction velocity v_w < 0 at the wall to the marching solver and show the layer stops growing (the asymptotic suction profile u = U(1 − e^(v_w y/ν))). |
| **Time** | unmeasured |

## 1. The engineering problem

Air does not slip over a wing. At the surface it is at rest relative to the skin (the no-slip condition of `P08.06`); a few millimetres away it moves at the full flight speed. All of that shearing happens in a thin **boundary layer**, and the shear stress at the bottom of it, τ_w = μ(∂u/∂y) at the wall, summed over the wetted surface, is the **skin-friction drag**, a large share of the drag of a clean aeroplane at cruise. The layer also decides whether the flow stays attached (`P08.11`), whether it stays smooth or becomes turbulent (`P08.12`), and how thick a region a CFD mesh must resolve near the wall (`P08.14`, `P08.18`).

Prandtl's insight of 1904 was that at large Reynolds number the viscous term matters only in this thin layer, so the Navier–Stokes equations simplify there, and the outer flow can be treated as inviscid (`P08.08`). For the simplest case, a flat plate at zero incidence, the simplified equations have a **similarity solution**: every station's velocity profile is the same curve once the height is scaled by √(νx/U). Blasius found it in 1908. Its one unknown number, the wall slope f″(0), sets the skin friction everywhere.

This lesson computes it three ways (RK4 shooting, a collocation solver, and a marching finite-difference solution of the full boundary-layer PDEs) and applies it to the light aeroplane's wing at cruise: 3000 m, 60 m/s, chord 1.5 m, semi-span 5 m. Success criteria, written first: f″(0) from shooting and from `solve_bvp` agree to 10⁻⁷ and match the published 0.332057; RK4 converges at order 4; θ√Re_x/x equals 2f″(0) (the momentum integral); the marched PDE, started from a uniform profile, collapses onto Blasius and converges to its skin friction at the order its scheme predicts; and the dimensional answers are invariant to the choice of units.

## 2. Understand and derive

### 2.1 The boundary-layer approximation

Take steady, two-dimensional, incompressible flow along a wall (x along it, y normal to it), with the outer velocity U(x). Navier–Stokes (`P08.06`) in two dimensions:

u u_x + v u_y = −p_x/ρ + ν(u_xx + u_yy),
u v_x + v v_y = −p_y/ρ + ν(v_xx + v_yy),
u_x + v_y = 0.

Suppose the layer has thickness δ(x) ≪ x. Then x-derivatives scale as 1/x and y-derivatives as 1/δ. Continuity gives v ~ Uδ/x: the normal velocity is small. In the x-momentum equation, u_xx ~ U/x² is negligible against u_yy ~ U/δ². Inertia u u_x ~ U²/x and viscosity ν u_yy ~ νU/δ² must be the same size inside the layer, or viscosity would not matter there at all. Equating them:

**δ/x ~ (νx/U)^(1/2)/x = Re_x^(−1/2),  Re_x = Ux/ν.**   (1)

The layer is thin exactly when Re_x is large. Every term of the y-momentum equation is smaller by δ/x than its x-momentum partner, so to leading order ∂p/∂y = 0: **the pressure is impressed on the layer by the outer flow**, which satisfies Bernoulli (`P08.04`), −p_x/ρ = U dU/dx. Prandtl's boundary-layer equations are

**u u_x + v u_y = U dU/dx + ν u_yy,  u_x + v_y = 0,**   (2)

with u = v = 0 at y = 0 (no slip, impermeable wall) and u → U(x) at the edge. The conservation laws are mass and x-momentum in differential form; the stress model is Newtonian and laminar (τ = μ ∂u/∂y). The equations are parabolic in x: information travels downstream only, which is why they can be **marched** (Step 6 does it). For a flat plate U is constant and the pressure term vanishes.

### 2.2 Similarity: the Blasius equation

A flat plate has no length scale in x, so the only height scale available at station x is √(νx/U). Try a profile that depends on y only through

η = y √(U/(νx)),  ψ = √(νUx) f(η),

where ψ is the stream function (u = ψ_y, v = −ψ_x, so continuity holds automatically). Then u = U f′(η) and v = ½√(νU/x)(η f′ − f). Substituting into (2) with dU/dx = 0 (every term turns out to carry the factor U²/x; Problem 2 has the details):

**f‴ + ½ f f″ = 0,  f(0) = f′(0) = 0,  f′(∞) = 1.**   (3)

A partial differential equation in x and y has become an ordinary differential equation in η: one curve f′(η) is the velocity profile at every station. The wall shear stress is τ_w = μU f″(0)√(U/(νx)), so the **local skin-friction coefficient** is

**c_f = τ_w/(½ρU²) = 2f″(0)/√Re_x.**   (4)

### 2.3 Thicknesses and the momentum integral

Three measures of thickness are used, each in units of x/√Re_x = √(νx/U):

- **δ99**, where u = 0.99U: a visual thickness, sensitive to the 0.99.
- **Displacement thickness** δ* = ∫₀^∞(1 − u/U)dy = lim(η − f) as η → ∞: the distance the wall would have to move outward for an inviscid flow to carry the same mass flux. The outer flow sees a body thickened by δ*, and the fluid just outside the layer moves outward at v_e = U dδ*/dx = ½(η f′ − f)U/√Re_x.
- **Momentum thickness** θ = ∫₀^∞(u/U)(1 − u/U)dy: the same idea for the momentum flux; the momentum the wall has removed from the flow is ρU²θ per unit span.

Their ratio is the **shape factor** H = δ*/θ. Integrating (2) across the layer for a flat plate (Problem 2c, and the control-volume momentum balance of `P08.03` applied to a slab between 0 and x) gives the **von Kármán momentum integral equation**

**dθ/dx = c_f/2,**   (5)

so the drag of one side of a plate of length L is ρU²θ(L) per unit span: the momentum deficit carried away in the wake. With θ ∝ √x and (4), (5) forces θ√Re_x/x = 2f″(0) exactly, a check the code must pass (the test asserts it to 10⁻⁴). Averaging (4) over 0 ≤ x ≤ L, ∫x^(−1/2)dx = 2√L, gives the **friction coefficient of the whole plate**

**C_F = (1/L)∫₀^L c_f dx = 4f″(0)/√Re_L = 2c_f(L).**   (6)

### 2.4 Solving (3): shooting

(3) is a boundary-value problem: two conditions at the wall and one at infinity. Write it as the first-order system Y = (f, f′, f″), Y′ = (f′, f″, −½ff″). From the wall we know f and f′ but not f″(0) = s. **Shooting** guesses s, integrates outward with RK4 (`P03.05`) to a finite η_max standing in for infinity, and reads the miss r(s) = f′(η_max; s) − 1. A root finder (`P03.02`) adjusts s until r = 0. The secant method from s = 0.2 (r = −0.287) and 0.5 (r = +0.314) reaches |r| < 10⁻¹² in five updates (s = 0.34327, 0.33123, 0.332062, 0.3320573391, 0.3320573372). The result depends on two numerical choices: the RK4 step (error ∝ h⁴) and η_max (too small an η_max forces the profile to reach U too soon and steepens it at the wall). Part 4 separates them.

There is a neater way. (3) is unchanged by the scaling f(η) → aF(aη), so integrating once with F″(0) = 1 and rescaling gives f″(0) = [F′(∞)]^(−3/2) with no iteration at all (Töpfer's transformation; Problem 2b). We still build shooting, because it works for problems that have no such symmetry (`P08.11`'s Falkner–Skan wedges with a continuation, and the stability problem of `P08.12`).

### 2.5 The wing at cruise

At 3000 m the ISA of `P08.01` gives T = 268.65 K, p = 70 109 Pa, ρ = 0.9091 kg/m³, μ = 1.6937 × 10⁻⁵ Pa·s, ν = 1.8630 × 10⁻⁵ m²/s, and at 60 m/s the dynamic pressure is q = 1636.4 Pa. With the constants of Step 5 (δ99 = 4.910, δ* = 1.72079, θ = 0.66411 in units of x/√Re_x, H = 2.5911, f″(0) = 0.332057):

| x (m) | Re_x | δ99 (mm) | δ* (mm) | θ (mm) | c_f | τ_w (Pa) |
|---|---|---|---|---|---|---|
| 0.010 | 3.221 × 10⁴ | 0.274 | 0.096 | 0.037 | 3.70 × 10⁻³ | 6.056 |
| 0.100 | 3.221 × 10⁵ | 0.865 | 0.303 | 0.117 | 1.170 × 10⁻³ | 1.915 |
| 0.300 | 9.662 × 10⁵ | 1.499 | 0.525 | 0.203 | 6.76 × 10⁻⁴ | 1.106 |
| 0.750 | 2.415 × 10⁶ | 2.369 | 0.830 | 0.320 | 4.27 × 10⁻⁴ | 0.699 |
| 1.500 | 4.831 × 10⁶ | 3.351 | 1.174 | 0.453 | 3.02 × 10⁻⁴ | 0.494 |

The layer is 3.35 mm thick at the trailing edge of a 1.5 m chord: 0.22 % of the chord. The wall shear falls by a factor of 12 from 1 cm to 1.5 m. Over both surfaces of the 15 m² wing (30 m² wetted), (6) gives **C_F = 6.04 × 10⁻⁴ and D = q C_F S_wet = 29.7 N**. For scale, the same wing turbulent from the nose carries 166 N (`P08.14`).

### 2.6 Where the model stops being valid

- **Near the leading edge** δ99/x = 4.91/√Re_x is not small: it reaches 0.1 at Re_x = 2.4 × 10³, 0.75 mm from the wing's leading edge (Problem 6). The singular c_f ∝ x^(−1/2) there is integrable, so the drag is unaffected.
- **Pressure gradient.** A real aerofoil accelerates the flow to a suction peak and then decelerates it. Blasius is the zero-gradient member of a family (Falkner–Skan, `P08.11`); in decelerating flow a laminar layer separates after only about a 12 % fall in U.
- **Transition.** The Blasius profile becomes unstable to small waves well before the trailing edge, and the layer turns turbulent somewhere past that point (`P08.12`); then (3) no longer applies.
- **Roughness, insects, rain** trip transition early. Compressibility is negligible at M ≈ 0.18; heat transfer and curvature are ignored.

## 3. Build the computational model

Open `starter/main.py`. `isa_troposphere`, `observed_order`, `_cumtrapz` and `main` are given.

### Step 1: the right-hand side

```python
def blasius_rhs(eta, Y):
    f, fp, fpp = Y
    return np.array([fp, fpp, -0.5 * f * fpp])
```

Unpacking `Y` along its first axis makes the same function work for one state (shape (3,)) and for `solve_bvp`, which passes a (3, m) array.

### Step 2: RK4

`rk4_blasius(fpp0, eta_max, n)` takes n equal steps h = η_max/n from Y(0) = (0, 0, fpp0) with the classical four stages, storing every state: it returns η (n + 1 points) and Y (n + 1 × 3).

### Step 3: shooting

`shoot_blasius(eta_max, n, guesses)` defines r(s) = f′(η_max; s) − 1 and runs the secant update s ← s − r(s)(s − s_old)/(r(s) − r(s_old)) until |r| < 10⁻¹², then returns s and the converged profile. Raise `RuntimeError` if it fails to converge, rather than returning a wrong number quietly.

### Step 4: the independent solver

`blasius_bvp(eta_max, tol)` passes `blasius_rhs`, the residuals (f(0), f′(0), f′(η_max) − 1) and the guess f′ = 1 − e^(−η) on 101 points to `scipy.integrate.solve_bvp`, which solves the whole boundary-value problem at once by collocation, adapting its mesh to the tolerance. It shares nothing with Steps 2–3 except the equation.

### Step 5: thicknesses and the plate

`blasius_integrals(eta, Y)` computes δ99 by linear interpolation where f′ crosses 0.99, δ* and θ by the trapezoid rule, H, f″(0), and the edge value of ½(ηf′ − f). `plate_boundary_layer(x, U, nu, rho, consts)` multiplies by x/√Re_x and applies (4); `friction_drag` applies (6).

### Step 6: marching the PDEs, the third route

`march_boundary_layer(U, nu, L, y_max, ny, nx)` solves (2) on a grid without assuming similarity. It starts at x₀ = 10⁻⁴L with u = U everywhere off the wall (nothing like Blasius) and steps downstream on nx steps spaced geometrically (the layer grows like √x, so steps should grow with x). At each step the unknown profile u^(j+1) satisfies

u^j_i (u^(j+1)_i − u^j_i)/Δx + v^j_i (u^(j+1)_(i+1) − u^(j+1)_(i−1))/(2Δy) = ν (u^(j+1)_(i+1) − 2u^(j+1)_i + u^(j+1)_(i−1))/Δy²,

backward Euler in x with the convecting velocities lagged, central differences in y: a tridiagonal system, solved with `scipy.linalg.solve_banded` (`P03.08`). Then continuity gives v^(j+1) by integrating −u_x upward with the trapezoid rule. `wall_cf_from_profile` takes the wall slope from the second-order one-sided difference (−3u₀ + 4u₁ − u₂)/(2Δy).

```python
lower = -vj / (2 * dy) - nu / dy**2;  diag = uj / dx + 2 * nu / dy**2;  upper = vj / (2 * dy) - nu / dy**2
rhs = uj * uj / dx;  rhs[-1] -= upper[-1] * U          # u = U at the top of the grid
```

Running `main()` prints

```
Blasius f''(0): shooting (RK4, n = 1000, eta_max = 10) 0.33205734; solve_bvp 0.33205734; published 0.332057
  difference shooting - bvp -1.7e-11
  delta99 4.9100, delta* 1.72079, theta 0.66411, H 2.5911, v_edge 0.8604 (bvp: theta 0.66410, H 2.5912)
RK4 convergence of f''(0) (eta_max = 10, reference n = 4000):
  n =   25: |error| = 1.479e-06
  n =   50: |error| = 1.202e-07
  n =  100: |error| = 8.407e-09
  n =  200: |error| = 5.535e-10
  observed orders [3.62 3.84 3.92]
domain truncation: f''(0) for eta_max = 4, 5, 6, 8, 10, 12 (n = 100 per unit eta):
  eta_max  4.0: 0.35274625
  eta_max  5.0: 0.33615234
  eta_max  6.0: 0.33256593
  eta_max  8.0: 0.33205919
  eta_max 10.0: 0.33205734
  eta_max 12.0: 0.33205734

wing at 3000 m: T 268.65 K, p 70109 Pa, rho 0.9091 kg/m^3, mu 1.6937e-05 Pa s, nu 1.8630e-05 m^2/s
     x (m)     Re_x      delta99 (mm)  delta* (mm)  theta (mm)   c_f        tau_w (Pa)
    0.010   3.221e+04      0.274        0.096       0.037   3.7006e-03  6.056
    ...
    1.500   4.831e+06      3.351        1.174       0.453   3.0216e-04  0.494
all-laminar wing (two sides, 15 m^2 planform): Re_c 4.831e+06, C_F 6.0431e-04, q 1636.4 Pa, friction drag 29.7 N

marching the boundary-layer PDEs from a uniform profile (U = 10 m/s, nu = 1.4607e-5, L = 1 m, ny = 300, nx = 1600):
  x = 0.010 m: c_f sqrt(Re_x) = 0.6674 (Blasius 0.6641)
  x = 0.100 m: c_f sqrt(Re_x) = 0.6623 (Blasius 0.6641)
  x = 0.501 m: c_f sqrt(Re_x) = 0.6629 (Blasius 0.6641)
  x = 1.000 m: c_f sqrt(Re_x) = 0.6631 (Blasius 0.6641)
```

and writes `blasius.png`: the profile f′(η) from both solvers with the marched PDE's profiles at three stations on top, and the RK4 convergence plot with a slope-4 line.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-10-laminar-boundary-layers`. The tests:

1. **The ODE and RK4**: `blasius_rhs` is the system and works on (3, m) arrays; RK4 starts at the wall, has the right shape, converges at order 3.5–4.5 on a fixed initial-value problem; f″(0) = 0 gives the trivial solution f ≡ 0.
2. **Shooting**: f(0) = f′(0) = 0 and f′(η_max) = 1 to 10⁻¹⁰; f″(0) within 2 × 10⁻⁶ of the published 0.332057 and within 10⁻⁷ of `solve_bvp`; the profiles of the two methods agree to 2 × 10⁻⁵; the finite-difference f‴ + ½ff″ is below 10⁻⁵ everywhere; f′ rises monotonically and f″ falls monotonically to below 10⁻⁶ of its wall value; f″(0) approaches the published value from above as η_max grows.
3. **Thicknesses**: δ99 = 4.91, δ* = 1.7208, H = 2.591, θ = 2f″(0) to 10⁻⁴, v_edge = δ*/2, and the `solve_bvp` profile gives the same θ and H.
4. **Dimensions**: δ ∝ √x and c_f ∝ x^(−1/2); c_f√Re_x = 2f″(0); τ_w = ½ρU²c_f; the same plate in millimetres gives the same Re_x and c_f and 1000× the thicknesses; the drag equals the numerically integrated wall shear to 10⁻⁶ and C_F = 2c_f(L).
5. **The marched PDE**: shape and boundary values; collapse onto f′(η) to 6 × 10⁻³ at x = 0.3 and 1 m on a moderate grid, with c_f√Re_x within 0.6 % of 0.664; observed order 0.8–1.3 in nx; and continuity: the growth of the mass-flux deficit equals the outflow v at the top of the grid to 10⁻⁸, with v > 0 there.
6. **main** writes the figure.

### Compare

`solve_bvp` and shooting agree to 1.7 × 10⁻¹¹ at η_max = 10: both are far more accurate than the six figures the published value is quoted to. `solve_bvp` needed 318 mesh points at tol = 10⁻⁸ and gave f″(0) within 3.1 × 10⁻⁶ even at tol = 10⁻³ (101 points). The marching solver is a different kind of check: it solves the PDE (2), not the ODE (3), and starts from the wrong profile, so its agreement (0.2 % in c_f at nx = 1600, converging at first order towards Blasius) tests the similarity assumption itself, not just the ODE solver. That is the only one of the three that would still work with a pressure gradient that has no similarity form (`P08.11` uses it).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **RK4 step.** f″(0) with n = 25, 50, 100, 200 steps on η_max = 10, against n = 4000. Predict the order. *Finding:* errors 1.48 × 10⁻⁶, 1.20 × 10⁻⁷, 8.41 × 10⁻⁹, 5.54 × 10⁻¹⁰; orders 3.62, 3.84, 3.92, approaching 4 from below as the step resolves the wall region. Against a `solve_bvp` reference at the same η_max = 9 (Problem 4) the orders are 3.93, 3.97, 4.03.
2. **Domain truncation.** f″(0) for η_max = 4–12 with 100 steps per unit η. Predict the sign of the error. *Finding:* 0.35275, 0.33615, 0.33257, 0.332059 and 0.3320573 at η_max = 4, 5, 6, 8, 10: errors of +2.1 × 10⁻², +4.1 × 10⁻³, +5.1 × 10⁻⁴, +1.9 × 10⁻⁶ and +1.2 × 10⁻⁹, always positive. A domain that ends inside the layer (δ99 = 4.91) is wrong by 6 % at η_max = 4; six figures need η_max ≥ 9.
3. **Collocation tolerance.** `blasius_bvp` at tol = 10⁻³, 10⁻⁵, 10⁻⁸. *Finding:* 101, 101 and 318 mesh points; f″(0) differs from the shooting value by 3.1 × 10⁻⁶, 5.8 × 10⁻⁷ and 1.7 × 10⁻¹¹. The tolerance controls the residual of the collocation equations, not directly the error in f″(0), which came out 300, 17 and 600 times smaller than the tolerance.
4. **The marched PDE.** U = 10 m/s, sea-level ν, L = 1 m, y_max = 3δ99(L), ny = 300, nx = 100–1600. Predict the order in nx. *Finding:* c_f√Re_x at x = L is 0.64830, 0.65673, 0.66047, 0.66222, 0.66306 for nx = 100, 200, 400, 800, 1600 (errors 1.58 × 10⁻², 7.38 × 10⁻³, 3.65 × 10⁻³, 1.90 × 10⁻³, 1.05 × 10⁻³; orders 1.10, 1.02, 0.94, 0.85); Richardson with p = 1 on the two finest gives 0.66391, 0.03 % from 2f″(0). The profile at x = 0.25 m and 1 m differs from f′(η) by at most 1.5 × 10⁻³ and 1.1 × 10⁻³. The PDE has forgotten its uniform starting profile.
5. **Speed and altitude.** Laminar friction drag of the wing (two sides, 30 m²) at 40, 60 and 80 m/s at sea level and at 3000 m. Predict how D scales with U. *Finding:* sea level 19.27, 35.40, 54.50 N; 3000 m 16.15, 29.67, 45.68 N; trailing-edge δ99 at 3000 m 4.10, 3.35, 2.90 mm. D ∝ U^1.5 (q ∝ U², C_F ∝ U^(−1/2)) and ∝ √(ρμ) (q ∝ ρ, C_F ∝ √(μ/ρ)): the thinner, less viscous air at altitude lowers the laminar drag at the same true airspeed by 16 %.

Raw data format: CSV tables `n, fpp0, err` (1); `eta_max, fpp0, err` (2); `tol, nodes, fpp0, diff` (3); `nx, ny, cf_sqrtRe, err` (4); `h_m, U_mps, Re_c, CF, D_N, delta99_TE_mm` (5); plus `blasius.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| f″(0) | published 0.332057 | shooting 0.33205734; `solve_bvp` 0.33205734 | 3 × 10⁻⁷ (published rounding); 1.7 × 10⁻¹¹ between methods | — |
| θ√Re_x/x | 2f″(0) = 0.664115 (5) | 0.66411 (trapezoid) | < 10⁻⁵ | trapezoid on h = 0.01 |
| v_e√Re_x/U | δ*/2 = 0.8604 | 0.8604 | < 10⁻⁵ | — |
| RK4 order | 4 | 3.62 → 3.92; 3.93 → 4.03 | — | pre-asymptotic at large steps |
| C_F | 2c_f(L) (6) | equal; integrated τ_w to 10⁻⁶ | — | — |
| c_f√Re_x from the PDE | 0.664115 | 0.66306 (nx = 1600); 0.66391 extrapolated | 0.16 %; 0.03 % | first-order x-steps |

Every number is reproduced by an independent route with the discrepancy its error source predicts. **The engineering answer: if the wing's boundary layer stayed laminar to the trailing edge, it would be 3.35 mm thick there (δ* 1.17 mm, θ 0.45 mm) and the wing would carry C_F = 6.04 × 10⁻⁴, 29.7 N of friction drag at 60 m/s and 3000 m.** The uncertainty in that number is not numerical (the three solvers agree to better than 0.2 %); it is the model. The laminar flat plate is a lower bound: the pressure gradient, transition and the surface condition all move the real drag, and `P08.12` and `P08.14` show that transition alone roughly triples it.

What to improve: the first-order x-step of the marching solver (Extension 2), if it is to be used for pressure-gradient flows; the domain truncation is already negligible at η_max = 10. What not to improve: more RK4 steps (10⁻⁹ is far below any physical uncertainty).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3) and of θ = 2f″(0)x/√Re_x; `blasius.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the wing's laminar boundary-layer thickness and friction drag and the assumption that most limits them.

### Adversarial self-checks

- **What if the speed doubles?** δ falls by √2 and C_F by √2, q rises 4×: the drag rises by 2^1.5 = 2.8, and Re_c doubles to 9.7 × 10⁶, which makes a laminar layer over the whole chord even less plausible.
- **What if viscosity goes to zero?** δ and c_f go to zero as √ν, but Re_x → ∞ and the layer becomes turbulent long before; the inviscid limit of `P08.08` (d'Alembert: no drag) is never reached by a real fluid.
- **Which assumption fails first in the real world?** Laminar flow. On a flat plate in quiet air the Blasius profile is predicted to become unstable at Re_x of order 10⁵ (`P08.12`: about 28 mm from the wing's leading edge at the Orr–Sommerfeld critical point Re_δ* ≈ 520, 41 mm by the e^N chart's onset) and the layer to become turbulent after roughly half the chord; rivets, insects or a gust trip it sooner.
- **What measurement could falsify the model?** A boundary-layer rake or a hot-wire traverse at x = 0.75 m: the predicted laminar layer is 2.4 mm thick with u/U = 0.5 at η = 1.54 (0.75 mm above the skin). A profile several times thicker and much fuller near the wall means the layer is turbulent there.
- **Does the result survive a change of method?** Shooting, collocation and marching agree to 0.2 % or better; f″(0) changes by 10⁻⁹ between η_max = 10 and 12; the result is unchanged when the units change (test 4).
- **Which decision changes if the uncertainty doubles?** None of the numerical uncertainty matters. If the transition point is uncertain by a factor of two, the drag estimate moves by tens of newtons, which is why transition (not the Blasius constants) is the question for a laminar-flow wing.

## Check yourself

1. Define δ99, δ*, θ and H, and give their Blasius values in units of x/√Re_x.
2. How does the laminar layer's thickness depend on x and on U? What happens to δ at the trailing edge if the speed is quadrupled?
3. Why does c_f fall along the plate, and why is the plate's average C_F exactly twice c_f at the trailing edge?
4. What is f″(0) physically, and why does a shooting method have to search for it?
5. Why does truncating the domain at too small an η_max make f″(0) too large?
6. List the assumptions of Prandtl's boundary-layer equations. Which fails near the leading edge?
7. Why is there an outward velocity v at the edge of a boundary layer on a flat plate, and how is it related to δ*?
8. Which conservation laws does this lesson use, in which form, and what does it assume about the stress?
