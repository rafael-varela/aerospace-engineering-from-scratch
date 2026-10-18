# P08.11  Integral Boundary Layer Methods

> Give up the profile, keep the momentum balance: one integral equation and one empirical closure turn a boundary layer in any pressure gradient into a single quadrature, and tell you where it will separate.

| | |
|---|---|
| **Learning outcomes** | 1. Derive von Kármán's momentum integral equation with a pressure gradient, and from it Thwaites's method, θ²Ue⁶ = θ₀²Ue₀⁶ + 0.45ν∫Ue⁵dx, with the closure functions l(λ), H(λ) and the separation criterion λ = −0.09. 2. Implement and test Thwaites's method for any edge-velocity distribution, including a stagnation point, and the exact Falkner–Skan wedge-flow solutions with `scipy.integrate.solve_bvp`. 3. Check Thwaites against Blasius on a flat plate, against Hiemenz at a stagnation point and across the Falkner–Skan family, and check its separation predictions against a marching finite-difference solution of the boundary-layer PDEs and Howarth's published separation point. 4. Decide how much deceleration a laminar boundary layer can take, and what that means for where an aerofoil's suction peak goes. |
| **Prerequisites** | `P08.10` (boundary-layer equations, Blasius, δ*, θ, H, the momentum integral on a flat plate, the marching solver), `P08.08` (potential flow round a cylinder), `P08.04` (Bernoulli: the outer flow sets the pressure), `P03.04` (quadrature), `P02.02` (boundary-value problems), `P03.02` (root finding). Skill check: you can evaluate ∫₀ˣ(1 − t)⁵dt in closed form and with the trapezoid rule. |
| **Engineering question** | Behind the suction peak of the light aeroplane's wing the flow decelerates. How much deceleration can a laminar boundary layer take before it separates, where does a model wing section's laminar layer separate, and where does the layer leave a round strut? Answers: from a fresh start only about a 12 % fall in edge velocity (Thwaites 12.3 %, a marching PDE 12.0 %, published 0.1199 L); on a model section with its peak (1.25U) at 0.3c the layer separates at x/c = 0.51 after a 7.1 % fall; on a cylinder in its potential-flow pressure field at 103.1° (Thwaites) or 104.6° (marched PDE) from the front stagnation point. |
| **Assumptions and limits** | Everything in `P08.10` (steady, two-dimensional, incompressible, laminar, thin layer) except that the edge velocity Ue(x) now varies and is prescribed by the outer flow; one-parameter closure: every profile is assumed to belong to a family labelled by λ alone; the closure fits are valid for −0.1 ≤ λ ≤ 0.1; the method stops at separation, where the boundary-layer approximation itself fails; the outer flow is not corrected for the layer (no viscous–inviscid interaction). |
| **Analytical work** | The momentum integral equation; Thwaites's linearisation and quadrature; closed forms for the flat plate (θ = 0.671x/√Re_x), wedge flows (λ = 0.45m/(5m + 1)), the stagnation point (θ² = 0.075ν/a) and Howarth's flow (x_s/L = 1 − 2.2^(−1/6) = 0.1231); the cylinder's separation angle is independent of U, R and ν. |
| **Computational deliverable** | `starter/main.py`: `thwaites_l_H`, `thwaites`, `separation_point`, `wedge_thwaites`, `falkner_skan`, `march_boundary_layer`, `march_separation` (`isa_nu`, `cumtrapz`, `observed_order`, `cylinder_edge_velocity`, `model_wing_edge_velocity` and `main` are given). Tests prove: the closure fits are continuous and give l = 0 at λ = −0.09; Thwaites reproduces its closed forms (flat plate, wedge flows, stagnation point, Howarth) and satisfies its own ODE to 10⁻³ and the momentum integral to 6 %; it is within 1 % of Blasius on θ; the cylinder's separation angle is independent of U, R, ν and units; Falkner–Skan reduces to Blasius and Hiemenz and loses its wall shear near β = −0.1988; the closure tracks the exact wedge flows to 8 %; the marched PDE matches Blasius's shear and finds Howarth's separation within 0.5 % of the published value. |
| **Practical activity** | **computational experiment**: Thwaites against the exact Falkner–Skan family; the separation point of Howarth's flow and of a cylinder by two methods and on refined grids; and the position of a model wing's suction peak against where its laminar layer separates. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `thwaites.png`, the derivations of the momentum integral equation and of Thwaites's quadrature, the closed forms for wedge flows and Howarth's flow, the experiment tables, and a half-page note on how much deceleration a laminar layer tolerates and what it means for the wing. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the momentum integral and Thwaites's quadrature are derived; all public tests pass; Thwaites's error against Blasius, Hiemenz and Falkner–Skan is reported with its sign; separation of Howarth's flow is computed by Thwaites and by the marched PDE and compared with the published value; the note says why a laminar layer cannot follow a steep pressure recovery. |
| **References** | `@mit-ocw-16-100-aerodynamics` Lecture 22 "Integral Boundary Layer Equations" and Lecture 23 "Correlation Methods for Integral Boundary Layers" (required-free: the integral equation, λ and Thwaites's quadrature); Lecture 24 "Method of Assumed Profiles" (supplementary-free: the Kármán–Pohlhausen alternative); Lecture 30 "Solutions of the Laminar Boundary Layer Equations" (Falkner-Skan Solutions) (supplementary-free). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 5 "Integral Boundary Layer Relationships": General Integral Boundary Layer Equations; Solution Method (required-free). `@leishman-intro-aerospace-flight-vehicles` "Boundary Layer Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/): von Kármán Momentum Integral Equation; Pressure Gradients & Flow Separation (required-free; link only). Optional textbooks for deeper study: White, *Viscous Fluid Flow*, and Schlichting and Gersten, *Boundary-Layer Theory* (paid; library copy; see `references/phase-08.md`). |
| **Extensions** | ★ Add Head's entrainment method for a turbulent layer and continue the wing's layer past transition (it is what `P10.06` couples to a panel method). ★ Make `march_boundary_layer` second order in x (Crank–Nicolson) and show that its Howarth separation converges at order 2 (up to the singularity). ★ Couple Thwaites to the outer flow: displace the cylinder by δ* and recompute Ue by superposition (`P08.08`); does the separation angle move? |
| **Time** | unmeasured |

## 1. The engineering problem

A wing makes lift by accelerating the flow over its upper surface to a **suction peak** and then decelerating it back to roughly free-stream speed at the trailing edge. By Bernoulli (`P08.04`) the deceleration is a rising pressure, an **adverse pressure gradient**. The outer, inviscid flow can climb that pressure hill using its kinetic energy, but the fluid near the wall has lost most of its kinetic energy to friction. Push it hard enough and the near-wall fluid stops, then reverses: the boundary layer **separates**, lift collapses and drag rises. Whether and where this happens is a design question for every aerofoil, diffuser, wind-tunnel contraction and engine intake.

The Blasius solution of `P08.10` cannot answer it, because a flat plate has no pressure gradient, and the marching PDE solver of `P08.10` can, but it is a full numerical calculation for each new Ue(x). In 1949 Thwaites found a way to reduce the problem to a single integral of Ue⁵. His method is still how many aerofoil codes start the boundary layer, and it is what `P10.06` couples to a panel method. This lesson derives it, builds it, and checks it against three independent routes: the exact Falkner–Skan similarity solutions (by `solve_bvp`), the marching PDE solver, and the published separation point of Howarth's linearly decelerating flow.

Success criteria, written first: Thwaites within 1–2 % of Blasius on a flat plate; θ constant at a stagnation point and within 7 % of Hiemenz's; λ constant along every wedge flow and equal to 0.45m/(5m + 1); Howarth's separation in closed form and within 4 % of the published 0.1199L, the marched PDE within 0.5 %; a cylinder's separation angle independent of speed, size and viscosity.

## 2. Understand and derive

### 2.1 The momentum integral equation with a pressure gradient

Start from Prandtl's equations (`P08.10`, equation (2)), now with U(x) written Ue(x), the edge velocity imposed by the outer flow:

u u_x + v u_y = Ue Ue′ + ν u_yy,  u_x + v_y = 0.

Integrate from the wall to a height h above the layer. Continuity turns the v u_y term into terms in u (Problem 2a), and the wall shear appears from ν u_yy. With δ* = ∫(1 − u/Ue)dy and θ = ∫(u/Ue)(1 − u/Ue)dy:

**d(Ue²θ)/dx + Ue Ue′ δ* = τ_w/ρ,  or  dθ/dx + (2 + H)(θ/Ue)dUe/dx = c_f/2.**   (1)

This is the x-momentum balance of a control volume one layer thick (`P08.03`): the rate at which the momentum deficit grows equals the wall's friction plus the pressure gradient's push on the deficit. With Ue′ = 0 it is `P08.10`'s dθ/dx = c_f/2. It is exact, but it has three unknowns (θ, H, c_f) for one equation: it needs a **closure**.

### 2.2 Thwaites's closure

Multiply (1) by 2Ueθ/ν and define the **pressure-gradient parameter** and the **shear parameter**

**λ = θ² Ue′/ν,  l = τ_w θ/(μ Ue).**

(1) becomes Ue d(θ²)/dx = 2[l − (2 + H)λ] ≡ F(λ). Thwaites's two observations, from exact and experimental laminar profiles: (i) l and H are, to a good approximation, functions of λ alone (a **one-parameter family** of profiles); (ii) F(λ) is close to the straight line **F ≈ 0.45 − 6λ**. Then Ue d(θ²)/dx + 6θ²Ue′ = 0.45ν, which is d(θ²Ue⁶)/dx = 0.45νUe⁵, and

**θ²(x) = [θ₀²Ue₀⁶ + 0.45ν ∫_(x₀)^x Ue⁵ dx′] / Ue⁶(x).**   (2)

Everything else follows: λ from θ, then l(λ) and H(λ) from the closure fits, then c_f = 2νl/(Ueθ) and δ* = Hθ. This lesson uses Cebeci and Bradshaw's fits to the exact solutions:

| range | l(λ) | H(λ) |
|---|---|---|
| 0 ≤ λ ≤ 0.1 | 0.22 + 1.57λ − 1.8λ² | 2.61 − 3.75λ + 5.24λ² |
| −0.1 ≤ λ < 0 | 0.22 + 1.402λ + 0.018λ/(λ + 0.107) | 2.088 + 0.0731/(λ + 0.14) |

l falls to zero at **λ = −0.09: laminar separation**. λ is the ratio of the pressure-gradient force to the viscous force across the layer: positive when the flow accelerates (thin, full profiles, high shear), negative when it decelerates.

### 2.3 Closed forms

- **Flat plate** (Ue = U): θ² = 0.45νx/U, so θ√Re_x/x = √0.45 = **0.6708**, against Blasius's 0.6641 (+1.0 %); λ = 0, l = 0.22, H = 2.61, c_f√Re_x = 2l/0.6708 = **0.6559**, against 0.6641 (−1.2 %).
- **Wedge flows** Ue = Cx^m (the flow over a wedge of half-angle βπ/2, β = 2m/(m + 1)): ∫x^(5m) = x^(5m+1)/(5m + 1), so θ² = 0.45νx/((5m + 1)Ue) and **λ = 0.45m/(5m + 1)**, the same at every station. These flows have exact similarity solutions (2.4), which makes them the natural test of the closure.
- **Stagnation point** (m = 1, Ue = ax): **θ² = 0.075ν/a**, constant: the inflow towards the wall exactly balances viscous thickening. (2) is 0/0 at x = 0; take the limit.
- **Howarth's flow** Ue = U₀(1 − x/L): λ = −0.075[(1 − x/L)^(−6) − 1], which reaches −0.09 at **x_s/L = 1 − 2.2^(−1/6) = 0.1231**, when Ue has fallen 12.3 %. Exact numerical solutions put it at 0.1199 (published).
- **Cylinder** (Ue = 2U sin(x/R), potential flow): with x = Rφ, U, R and ν cancel from λ (Problem 2d). The separation angle is universal: 103.1° by (2).

### 2.4 The exact check: Falkner–Skan

For Ue = Cx^m the boundary-layer equations have similarity solutions u = Ue f′(η), η = y√((m + 1)Ue/(2νx)), with

**f‴ + f f″ + β(1 − f′²) = 0,  f(0) = f′(0) = 0,  f′(∞) = 1,  β = 2m/(m + 1).**   (3)

β = 0 is Blasius (in this scaling f″(0) = 0.332057√2 = 0.4696), β = 1 is **Hiemenz's stagnation-point flow** (f″(0) = 1.23259), and the wall shear vanishes at β = −0.1988, the separating profile. In Thwaites's variables, λ_exact = θ_η²·2m/(m + 1) and l_exact = f″(0)θ_η. The code solves (3) with `solve_bvp`, using a continuation in β from 0 for β < 0 so it stays on the attached branch.

### 2.5 Where the method stops being valid

- **At separation.** λ = −0.09 is where the method stops, not where the flow stops: beyond it the layer leaves the wall, the pressure is no longer imposed by the outer flow, and neither (1) nor the marching PDE applies (Problem 5c).
- **One parameter is not enough** in strongly accelerated or decelerated flow: errors grow from 1 % on a flat plate to 6 % in θ at a stagnation point (2.4 and part 4). Two-parameter methods (Head's, or the lag-entrainment methods for turbulent layers) do better.
- **Laminar only.** The closure belongs to laminar profiles; after transition (`P08.12`) a turbulent closure is needed (`P08.14`, `P10.06`).
- **The outer flow is fixed.** A real separation changes the outer flow, which is why a real cylinder's layer separates far ahead of the 103–105° that the potential-flow pressure distribution gives.

## 3. Build the computational model

Open `starter/main.py`. `isa_nu`, `cumtrapz`, `observed_order`, `cylinder_edge_velocity`, `model_wing_edge_velocity` and `main` are given.

### Step 1: the closure

`thwaites_l_H(lam)` evaluates the table of 2.2 with `np.where`, after clipping λ to [−0.1, 0.1].

### Step 2: the quadrature

`thwaites(x, Ue, dUedx, nu, theta0)` computes `cumtrapz(Ue**5, x)`, then (2), with the stagnation-point limit 0.075ν/Ue′(0) when Ue[0] = 0; then λ, (l, H), δ*, and c_f (infinite where Ue = 0). `separation_point(x, lam)` interpolates linearly where λ first falls through −0.09, and returns nan if it never does.

```python
I = cumtrapz(Ue**5, x)
th2 = (theta0**2 * Ue[0]**6 + 0.45 * nu * I) / Ue**6
if Ue[0] == 0.0:
    th2[0] = 0.075 * nu / dU[0]
```

The trapezoid rule is poor on the first few cells next to a stagnation point (on the first cell it gives h⁶/2 for ∫x⁵dx instead of h⁶/6), which shows as a spike in λ near x = 0 in `thwaites.png`; it is gone by the twentieth station (the test checks θ from there on).

### Step 3: wedge flows, both ways

`wedge_thwaites(m)` is the closed form of 2.3 in Falkner–Skan's scaling, θ_η = √(0.45(m + 1)/(2(5m + 1))). `falkner_skan(beta)` solves (3) with `solve_bvp` (guess f′ = 1 − e^(−η)), continuing from β = 0 in six steps for β < 0, and returns f″(0), δ*, θ, H and Thwaites's λ and l.

### Step 4: the marching PDE with a pressure gradient

`march_boundary_layer(x, Ue, dUedx, nu, y_max, ny)` is `P08.10`'s marcher with the source Ue Ue′ added to the right-hand side and Ue(x) at the top boundary. It records the wall shear and stops at the first station where it is not positive. `march_separation(x, tau)` fits a straight line to τ² over the last ten stations and returns its zero: near separation τ_w ∝ √(x_s − x) (Goldstein's singularity), so τ² is linear.

Running `main()` prints

```
flat plate (U = 60.0 m/s, nu = 1.8630e-05 m^2/s, x = 1.5 m, Re_x = 4.831e+06):
  theta sqrt(Re_x)/x: Thwaites 0.6708, Blasius 0.6641 (+1.01 %)
  c_f sqrt(Re_x): Thwaites 0.6559, Blasius 0.6641 (-1.23 %)
  H: Thwaites 2.610, Blasius 2.591

wedge flows Ue = C x^m: exact (Falkner-Skan, solve_bvp) against Thwaites
     m       beta     lam exact  lam Thw    l exact   l Thw     H exact   H Thw    theta_eta err
  -0.0904  -0.1988    -0.0681   -0.0742    0.0041   0.0751     3.970    3.200    +4.40 %
  -0.0600  -0.1277    -0.0360   -0.0386    0.1414   0.1558     2.904    2.809    +3.55 %
  -0.0300  -0.0619    -0.0152   -0.0159    0.1896   0.1946     2.701    2.677    +2.18 %
   0.0000   0.0000     0.0000    0.0000    0.2205   0.2200     2.591    2.610    +1.01 %
   0.1000   0.1818     0.0310    0.0300    0.2764   0.2655     2.422    2.502    -1.59 %
   0.3333   0.5000     0.0613    0.0562    0.3249   0.3026     2.297    2.416    -4.24 %
   1.0000   1.0000     0.0855    0.0750    0.3603   0.3276     2.216    2.358    -6.32 %
Hiemenz stagnation point: f''(0) = 1.23259, theta sqrt(a/nu) = 0.2923 exact, 0.2739 Thwaites; delta* sqrt(a/nu) = 0.6479, H = 2.216
Falkner-Skan at beta = -0.19883 (m = -0.0904): f''(0) = 0.00236, H = 4.009, lambda = -0.0681

Howarth flow Ue = U (1 - x/L), L = 1.0 m: separation x_s/L Thwaites 0.1231, marched PDE 0.1202, published 0.1199
  Thwaites: the layer separates after Ue has fallen by 12.3 %
cylinder (potential-flow Ue): separation at 103.1 deg (Thwaites), 104.6 deg (marched PDE), from the front stagnation point

model wing upper surface (U = 60.0 m/s at 3000 m, chord 1.5 m, peak 1.25 U, trailing edge 0.95 U):
  peak at 0.1 c: theta at peak 0.106 mm; laminar separation at x/c = 0.506, 609 mm behind the peak, after Ue fell 10.8 %
  peak at 0.3 c: theta at peak 0.183 mm; laminar separation at x/c = 0.508, 312 mm behind the peak, after Ue fell 7.1 %
  peak at 0.5 c: theta at peak 0.236 mm; laminar separation at x/c = 0.540, 60 mm behind the peak, after Ue fell 1.9 %
  peak at 0.7 c: theta at peak 0.280 mm; laminar separation at x/c = 0.700, 0 mm behind the peak, after Ue fell 0.0 %
```

and writes `thwaites.png`: l(λ) and H(λ) with the exact Falkner–Skan points, and λ(x) for Howarth's flow and the cylinder down to the separation line.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-11-integral-boundary-layer-methods`. The tests:

1. **Closure**: l = 0.22 at λ = 0 from both sides; H continuous to 0.01; |l(−0.09)| < 2 × 10⁻³ with 3.4 < H < 3.8; l increasing and H decreasing in λ.
2. **Quadrature**: flat plate θ√Re_x/x = √0.45 to 10⁻¹² and within 1.1 % of Blasius, c_f within 1.5 %; Thwaites's ODE Ue(θ²)′ + 6θ²Ue′ = 0.45ν to 10⁻³ and the momentum integral (1) to 6 % away from separation; λ constant along three wedge flows; θ = √(0.075ν/a) at a stagnation point; Howarth's λ(x) in closed form and its separation within 4 % of the published value; `separation_point` interpolates and returns nan without separation; the cylinder's angle independent of (U, R, ν) to 10⁻⁶; unit invariance.
3. **Falkner–Skan**: β = 0 gives 0.332057√2 and H = 2.591; Hiemenz f″(0) = 1.23259 (published) and Thwaites's θ 6 % low; f″(0) < 0.03 and H > 3.5 at β = −0.198 with λ between −0.075 and −0.06; the closure within 8 % of l and 7 % of H, and θ_η within 7 %, for four β.
4. **Marching PDE**: a flat plate runs to the end with c_f√Re_x within 1 % of 0.664; Howarth's flow stops, with τ_w falling monotonically, and `march_separation` lands within 0.5 % of 0.1199; Thwaites within 3 % of it; the τ² extrapolation is exact for a square root.
5. `main` writes the figure.

### Compare

Three routes, three kinds of agreement. Against the exact Falkner–Skan solutions (the "established tool" here is `solve_bvp`), Thwaites's θ is within 1–2 % for mild gradients and 4–6 % at the ends of the family, always erring the same way: too thick in decelerating flow, too thin in accelerating flow. Against the marching PDE, which uses no closure at all, Thwaites puts Howarth's separation 2.4 % late (0.1231 against 0.1202) and the cylinder's 1.5° early (103.1° against 104.6°). The marched value itself converges towards the published 0.1199 as the grid is refined (part 4). The closure's separation λ (−0.09) is not the exact wedge value (−0.068); it works for the separation *point* of decelerating flows like Howarth's, which is what matters for design (Problem 4c).

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **The closure across the family.** Exact (Falkner–Skan) against Thwaites for m = −0.0904 to 1. Predict where Thwaites is best. *Finding:* the table above; θ error +4.4 % at the separating profile, +1.0 % at the flat plate, −6.3 % at the stagnation point; H error up to 0.8 near separation (3.20 against 3.97) because the fit was not made to reproduce the separating profile.
2. **Howarth's flow on refined grids.** The marching PDE on linear x-grids with nx = 1000, 2000, 4000, 8000 (U = 10 m/s, sea level). Predict whether it converges to 0.1199. *Finding:* x_s/L = 0.12047, 0.12036, 0.12024, 0.12009, heading for the published 0.1199; Thwaites's 0.1231 is independent of the grid (its quadrature converges at once: 0.123141 at 501 to 4001 points). Geometric x-grids, good for a flat plate, are worse here (0.12099 at nx = 8000) because the steps near separation are long.
3. **The cylinder.** Thwaites and the marched PDE, (nx, ny) from (2000, 150) to (16000, 600). *Finding:* marched 104.84°, 104.69°, 104.61°, 104.54°; Thwaites 103.11° for every grid from 751 to 6001 points and for (U, R, ν) = (20 m/s, 0.05 m, 1.5 × 10⁻⁵), (5, 0.01, 10⁻⁶) and (60, 1.0, 1.86 × 10⁻⁵).
4. **A developed layer meets a diffuser.** Problem 3's duct (U₀ = 10 m/s, L_d = 3 m). Predict the effect of 0.35 m of straight duct upstream. *Finding:* without it the layer separates 0.475 m into the diffuser after a 13.7 % fall in Ue; with it, after 0.158 m and a 5.0 % fall (marched PDE: 0.153 m).
5. **Where to put the suction peak.** The given model upper-surface distribution (an accelerating wedge flow, m = 0.1, to a peak of 1.25U at x_p, then a linear fall to 0.95U at the trailing edge), wing at 3000 m and 60 m/s, x_p = 0.1c to 0.7c. Predict whether a later peak delays separation. *Finding:* separation at x/c = 0.506, 0.508, 0.540, 0.700, i.e. 609, 312, 60 and 0 mm behind the peak, after Ue fell 10.8, 7.1, 1.9 and 0.0 %. A later peak keeps the layer attached further back, but the layer arrives at the peak thicker (θ = 0.106 → 0.280 mm) and meets a steeper recovery, so the margin behind the peak shrinks to nothing: λ just behind a 0.7c peak is −0.168. A laminar layer cannot make the recovery of any of these sections; the wing relies on transition (`P08.12`) ahead of separation.

Raw data format: CSV tables `m, beta, theta_exact, theta_thw, cf_exact, cf_thw, H_exact, H_thw` (1); `nx, xs_over_L` (2); `nx, ny, angle_deg` (3); `l0_m, lambda0, xs_m, Ue_fall_pct` (4); `x_peak_over_c, theta_peak_mm, xs_over_c, Ue_fall_pct` (5); plus `thwaites.png`.

## 5. Validate and improve

| quantity | analysis (§2) | code | difference | cause |
|---|---|---|---|---|
| flat-plate θ√Re_x/x | Blasius 0.6641 (published) | Thwaites 0.6708 | +1.0 % | the closure |
| stagnation θ√(a/ν) | Hiemenz 0.2923 (`solve_bvp`) | Thwaites 0.2739 | −6.3 % | one-parameter closure |
| Howarth x_s/L | published 0.1199; closed form 0.12314 | Thwaites 0.12314; marched 0.1202 → 0.1201 | +2.7 %; +0.2 % | closure; first-order x-steps |
| cylinder separation | universal angle (2.3) | Thwaites 103.1°; marched 104.6° → 104.5° | 1.4° | closure |
| wedge λ | 0.45m/(5m + 1) | constant along x to 0.2 % | — | trapezoid start-up |

The engineering answer: **a laminar boundary layer starting fresh separates after its edge velocity has fallen by only about 12 % (Thwaites 12.3 %, marched PDE 12.0 %, published 0.1199); one that arrives already developed takes even less (5 % in the diffuser of experiment 4); and behind a model wing's suction peak (1.25U at 0.3c) it separates at x/c = 0.51, after a 7 % fall.** On a round strut in its potential-flow pressure field the laminar layer would separate at about 103–105°; a real strut's wake moves that point forward. A laminar layer cannot follow the pressure recovery of a wing; the turbulent layer of `P08.14`, which carries much more momentum near the wall, can.

What to improve: a two-parameter closure (or a turbulent one after transition) if the method is to be used beyond λ ≈ −0.05; viscous–inviscid coupling if the separation point itself is the question. What not to improve: the quadrature grid (Thwaites converged at 501 points).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1) and (2) and the closed forms of 2.3; `thwaites.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating how much deceleration a laminar layer tolerates, where the model wing's layer separates, and the assumption that most limits the answer.

### Adversarial self-checks

- **What if the speed doubles?** For Howarth's flow and the cylinder, nothing: x_s/L and the angle are independent of U and ν (2.3). θ shrinks as 1/√U, but λ does not change, so the separation point does not move. Speed matters through transition, not through laminar separation.
- **What if the deceleration is gentler?** For Ue = U₀/(1 + x/L_d) separation always comes after a 13.7 % fall from a fresh start, whatever L_d: a laminar layer's tolerance is a fraction of Ue, not a gradient.
- **Which assumption fails first in the real world?** That the layer stays laminar up to the computed separation point. At the wing's Reynolds numbers transition (`P08.12`) usually comes first or immediately after, and laminar separation becomes a short bubble rather than a stall.
- **What measurement could falsify the model?** Surface pressure taps on a model wing: a laminar separation bubble shows as a pressure plateau starting near the predicted x_s; oil-flow visualisation shows the separation line.
- **Does the result survive a change of method or grid?** Thwaites, the marched PDE and the published value agree to 3 % on Howarth's flow; Thwaites does not change with the grid, and the marched value converges towards the published one.
- **Which decision changes if the uncertainty doubles?** A 3 % uncertainty in x_s changes nothing in practice; a factor-of-two uncertainty in the transition point would decide whether the wing has a laminar bubble at all, which is why `P08.12` comes next.

## Check yourself

1. Define λ and say what its sign means physically.
2. Write von Kármán's momentum integral equation with a pressure gradient and say what each term represents.
3. What two assumptions make Thwaites's method work, and why is it called a one-parameter method?
4. By how much does Thwaites miss Blasius on a flat plate, in θ and in c_f?
5. Why is θ constant near a two-dimensional stagnation point?
6. By how much can Ue fall before a laminar layer separates, starting fresh? Why does a layer that arrives already developed tolerate less?
7. Why can the marching PDE solver not continue past separation?
8. Thwaites and the marched PDE put separation on a cylinder at 103–105°. Why does a real subcritical cylinder's layer separate earlier?
