# P08.15  Turbulence Modeling

> Averaging the Navier–Stokes equations leaves one term nobody knows: the Reynolds stress. A turbulence model is a guess for it, calibrated on a few simple flows. Know which flows, and you know where to trust it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Reynolds-averaged momentum equation for fully developed channel flow, show that the total shear stress is exactly linear across the channel, and state what the closure problem is. 2. Close it with Prandtl's mixing length and van Driest's damping, reduce the model to a quadratic for the velocity gradient at each point, and integrate it across the channel to the velocity profile, the bulk velocity and the skin friction. 3. Implement and test the model against its exact limits (laminar Poiseuille flow, the viscous sublayer, the log layer), show the bulk velocity converges at second order, and compare the friction with the log-law friction law and with Dean's correlation. 4. State the standard k–ε model and its constants, derive its log-layer equilibrium solution, and show by substitution that the ε equation closes only if κ² = σ_ε √C_μ (C_ε2 − C_ε1). 5. Decide how far a two-constant model can be trusted for a duct's pressure drop and how to place the first grid cell of a CFD mesh. |
| **Prerequisites** | `P08.13` (Reynolds decomposition, Reynolds stresses), `P08.14` (wall units, the law of the wall, the log law, skin-friction correlations), `P08.06` (Navier–Stokes, plane Poiseuille flow), `P08.03` (control-volume momentum balance), `P03.02` (root finding), `P03.04` (quadrature), `P03.11` (observed order), `P01.17` (least squares). Skill check: you can convert a wall distance in metres to y⁺ given u_τ and ν. |
| **Engineering question** | A flat cabin-air duct with a 50 mm gap carries air at 20 °C at a bulk speed of 5 m/s. What wall shear stress and pressure gradient does a mixing-length model predict, how close is it to Dean's correlation, and how thin must the first cell of a CFD mesh be to resolve the wall? The answers: Re_b = 16 630, Re_τ = 488; τ_w = 0.104 Pa and dp/dx = −4.15 Pa/m, 7.3 % above Dean's 3.87 Pa/m with Prandtl's outer mixing length and 7.8 % below it with Nikuradse's; the first cell for y⁺ = 1 is 51 µm thick (1.5 mm for wall functions at y⁺ = 30). |
| **Assumptions and limits** | Incompressible, statistically steady, fully developed flow between infinite parallel plates (no side walls, no entrance length); smooth walls; Newtonian fluid; a scalar eddy viscosity aligned with the mean strain (the Boussinesq hypothesis); the mixing length is algebraic, so the turbulence has no memory of upstream history; van Driest's damping constant A⁺ = 26 and κ = 0.41 are typical calibrated values, not laws. The k–ε analysis is restricted to the log layer, where production equals dissipation. The model fails where history matters: separation, strong pressure gradients, curvature, rotation, transition. |
| **Analytical work** | Reynolds averaging of the x-momentum equation; the linear total-stress profile τ⁺ = 1 − y⁺/Re_τ from a control-volume balance; the mixing-length closure and the quadratic for du⁺/dy⁺ with its cancellation-free root; the laminar (κ → 0) and fully turbulent limits; the log law from l = κy; the log-law friction law; the k–ε equations, their log-layer solution k = u_τ²/√C_μ, ε = u_τ³/(κy), ν_t = κu_τy; the constraint on the constants. |
| **Computational deliverable** | `starter/main.py`: `total_stress`, `mixing_length`, `velocity_gradient`, `mixing_length_profile`, `laminar_profile`, `bulk_velocity`, `friction`, `fit_log_law`, `log_law_friction`, `re_tau_for_re_b`, `duct_flow`, `keps_log_layer`, `keps_residuals`, `kappa_from_keps`, `c_mu_from_stress_ratio` (`channel_grid`, `cumtrapz0`, `log_law`, `dean_cf`, the constants and `main` are given). Tests (13, under 2 s) prove: the total stress is the linear momentum balance; the gradient solves the quadratic to 10⁻¹² with both limits; the damping and outer mixing lengths have their limits; κ → 0 reproduces Poiseuille flow with c_f Re_b = 12; u⁺ = y⁺ in the sublayer and the slope is 1/κ in the log layer, where ν_t = κy⁺√τ⁺; viscous plus turbulent stress equals the total; the log-law fit recovers κ and A⁺ sets B; the bulk velocity converges at order 2; the friction lies within 20 % of Dean's correlation and the two outer closures bracket it; the results depend only on Re_b; the k–ε log layer has P = ε and −⟨u′v′⟩/k = √C_μ; κ implied by the constants is 0.433; the ε residual converges to zero at order 2 only for that κ. |
| **Practical activity** | **computational experiment**: a Reynolds-number sweep of the mixing-length channel against two friction laws, a sensitivity study of the friction to κ, A⁺ and the outer mixing length, and the k–ε log-layer consistency check with standard and modified constants. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `turbulence_models.png`, the derivations of τ⁺ = 1 − y⁺/Re_τ, of the velocity-gradient root and of the k–ε constraint, the sweep and sensitivity tables, and a half-page note on the duct's pressure drop with its model uncertainty. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the total-stress profile and the k–ε constraint are derived; all public tests pass; the note quotes the duct's pressure gradient with a band set by the sensitivity study, not a single number, and says what the model cannot predict. |
| **References** | `@mit-ocw-2-29-numerical-fluid-mechanics` Lecture 25: Reynolds-averaged Navier-Stokes (RANS); Turbulence closures: Eddy viscosity and diffusivity, Mixing-length Models, k-ε Models (required-free). `@leishman-intro-aerospace-flight-vehicles` chapter "Turbulent Flows": Eddy Viscosity Models; Turbulence Models; Numerical Solution of the RANS Equations (required-free; https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/). `@mit-ocw-2-27-turbulent-flow-transport` "Concepts in Turbulence" 3.8 Prandtl's (momentum-transfer) mixing length hypothesis; "Introduction to Turbulence Models" 8.1–8.6, including 8.4 the calibration of the constants and 8.5 the equation for ε (supplementary-free). |
| **Extensions** | ★ Download a Lee & Moser channel DNS mean profile (`@lee-moser-2015-channel-dns-data`, Re_τ = 550 or 1000; follow their citation request; do not redistribute) and plot the mixing-length profile against it: where does the model miss, and by how much in U_b⁺? ★ Replace the algebraic mixing length by a one-dimensional k–ε solver with wall functions at y⁺ = 30 and compare its c_f with this lesson's. ★ Add a Couette (moving-wall) case, where τ⁺ = 1 everywhere, and show the log layer then extends to the centre. |
| **Time** | unmeasured |

## 1. The engineering problem

Most flows on an aircraft are turbulent: the boundary layer over most of the wing, the air in the cabin ducts, the propellant in the feed lines, the exhaust of every engine. The Navier–Stokes equations describe them exactly, but resolving every eddy is out of reach for engineering work: a direct numerical simulation of a channel at a friction Reynolds number of 5200 used billions of grid points. Engineers solve instead for the **mean** flow, and that costs one unknown more than there are equations. A **turbulence model** supplies the missing relation. Every RANS (Reynolds-averaged Navier–Stokes) CFD result you will ever see rests on one.

This lesson builds the simplest useful model, Prandtl's mixing length with van Driest's correction at the wall, for the flow in which turbulence models are calibrated: fully developed flow between two parallel plates. That geometry is a fair model of a flat cabin-air duct. It answers a real sizing question (what pressure gradient must the fan supply per metre of duct?) and a CFD question (how thin must the first grid cell be?). It then looks at the k–ε model, the workhorse of industrial CFD, and shows by a direct substitution that its constants are tied together by the log law: the model was calibrated to reproduce it, imperfectly.

Success criteria, written first: the model reproduces Poiseuille flow exactly when the mixing length is zero; it has a viscous sublayer u⁺ = y⁺ and a log layer of slope 1/κ; its bulk velocity converges at second order with the grid; its friction coefficient lies within 20 % of Dean's correlation across Re_τ = 180–5200; and the k–ε log-layer residual vanishes only for the κ that its constants imply.

## 2. Understand and derive

### 2.1 Reynolds averaging and the closure problem

Split each velocity component into a mean and a fluctuation, u = U + u′, with ⟨u′⟩ = 0 (`P08.13`). Averaging the incompressible Navier–Stokes equations leaves them unchanged except for one new term: the average of the product of fluctuations, ρ⟨u′v′⟩, which acts on the mean flow like a stress, the **Reynolds shear stress**. For flow between plates 2h apart, fully developed (nothing changes with x except the pressure) and statistically steady, V = 0 and the mean x-momentum equation is

0 = −dP/dx + d/dy [μ dU/dy − ρ⟨u′v′⟩].   (1)

That is one equation for two unknowns, U(y) and ⟨u′v′⟩(y): the **closure problem**. Conservation laws cannot close it; a model must.

### 2.2 The total stress is exactly linear

Integrate (1) from the wall, where the viscous stress is τ_w and the Reynolds stress is zero (no slip):

τ(y) ≡ μ dU/dy − ρ⟨u′v′⟩ = τ_w + y dP/dx.

By symmetry τ = 0 on the centreline y = h, so dP/dx = −τ_w/h: **the pressure gradient balances the wall shear on both walls**, which is also a control-volume momentum balance on a slab of the duct (`P08.03`). Hence

**τ⁺ ≡ τ/τ_w = 1 − y⁺/Re_τ,**   (2)

in wall units (`P08.14`): u_τ = √(τ_w/ρ), y⁺ = y u_τ/ν, u⁺ = U/u_τ, Re_τ = h u_τ/ν. Equation (2) holds for any closure: it is exact. The model only decides how τ⁺ is shared between the viscous part du⁺/dy⁺ and the turbulent part −⟨u′v′⟩⁺.

### 2.3 The mixing length

Prandtl's idea (`@mit-ocw-2-27-turbulent-flow-transport` "Concepts in Turbulence" 3.8): an eddy moves fluid a distance l (the mixing length) across the mean gradient, carrying a velocity difference l dU/dy, so the Reynolds stress scales as

−⟨u′v′⟩ = ν_t dU/dy,  ν_t = l² |dU/dy|.   (3)

ν_t is the **eddy viscosity** (the Boussinesq hypothesis: turbulence acts like an extra viscosity). Near a wall the largest eddies are limited by the distance to it, so l = κy, with κ ≈ 0.41 the von Kármán constant. Right at the wall the fluctuations are damped by viscosity; **van Driest** multiplied l by D = 1 − exp(−y⁺/A⁺) with A⁺ = 26, so l⁺ ≈ κy⁺²/A⁺ for y⁺ ≪ A⁺. Far from the wall l = κy keeps growing; Nikuradse's pipe-flow fit l/h = 0.14 − 0.08(1 − y/h)² − 0.06(1 − y/h)⁴ caps it at 0.14h and tends to 0.4y near the wall. Both outer forms are in the code.

### 2.4 One quadratic per point

Substitute (3) in (2), in wall units, with g = du⁺/dy⁺ ≥ 0:

(1 + l⁺² g) g = τ⁺  ⇒  l⁺² g² + g − τ⁺ = 0.

The positive root is (−1 + √(1 + 4l⁺²τ⁺))/(2l⁺²), which suffers cancellation when l⁺ is small. Multiply top and bottom by (1 + √…):

**g = 2τ⁺ / (1 + √(1 + 4 l⁺² τ⁺)).**   (4)

Two limits check it. With l⁺ = 0, g = τ⁺: laminar flow, and integrating gives u⁺ = y⁺ − y⁺²/(2Re_τ), plane Poiseuille flow (`P08.06`), with U_b⁺ = Re_τ/3 and c_f Re_b = 12. With l⁺ large, g → √τ⁺/l⁺; with l⁺ = κy⁺ and τ⁺ ≈ 1 near the wall, g = 1/(κy⁺), and integrating gives **u⁺ = (1/κ) ln y⁺ + B**: the log law (`P08.14`) comes out of the model, with B set by the damping. Then u⁺(y⁺) = ∫₀^{y⁺} g dy⁺ by the trapezoid rule on a grid clustered near the wall, and the bulk velocity, Reynolds number and friction coefficient follow from the definitions:

U_b⁺ = (1/Re_τ)∫₀^{Re_τ} u⁺ dy⁺,  Re_b = U_b(2h)/ν = 2U_b⁺Re_τ,  c_f = τ_w/(½ρU_b²) = 2/U_b⁺².   (5)

### 2.5 Friction laws

**The log-law friction law.** If the log law held from wall to centreline, U_b⁺ = (1/Re_τ)∫₀^{Re_τ}[(1/κ)ln y⁺ + B] dy⁺ = (1/κ) ln Re_τ + B − 1/κ. With (5) that is an implicit c_f(Re_b); it ignores the viscous layer (which lowers U_b⁺) and the wake (which raises it), and the two nearly cancel.

**Dean's correlation.** Dean (1978) fitted channel experiments with c_f = 0.073 Re_b^(−1/4) for roughly 6 × 10³ < Re_b < 6 × 10⁵. It is a **correlation**, a fit with a few per cent of scatter, not a law and not a measurement of your duct.

### 2.6 The k–ε model and its log layer

The standard k–ε model (Launder and Spalding, 1974; `@mit-ocw-2-29-numerical-fluid-mechanics` Lecture 25 "k-ε Models"; `@mit-ocw-2-27-turbulent-flow-transport` 8.2–8.5) carries two transport equations, for the turbulent kinetic energy k and its dissipation rate ε, and sets ν_t = C_μ k²/ε. For fully developed flow they read

0 = d/dy[(ν_t/σ_k) dk/dy] + P − ε,   0 = d/dy[(ν_t/σ_ε) dε/dy] + (C_ε1 P − C_ε2 ε) ε/k,   P = ν_t (dU/dy)²,   (6)

with the published standard constants **C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_k = 1.0, σ_ε = 1.3**.

**The log-layer solution.** In the log layer the stress is constant, ν_t dU/dy = u_τ², and production balances dissipation, P = ε. With dU/dy = u_τ/(κy) these give ν_t = κu_τy, ε = u_τ³/(κy), and from ν_t = C_μk²/ε, **k = u_τ²/√C_μ**: k is constant, so the k equation is satisfied exactly. The ratio −⟨u′v′⟩/k = √C_μ ≈ 0.3 is what measurements show in the log layer, and that is where C_μ = 0.09 comes from. Note that ν_t = κu_τy is the mixing-length eddy viscosity: the two closures agree in the log layer.

**The constraint.** Substitute into the ε equation. The diffusion term is d/dy[(κu_τy/σ_ε)(−u_τ³/(κy²))] = u_τ⁴/(σ_ε y²); the source is (C_ε1 − C_ε2)ε²/k = (C_ε1 − C_ε2)√C_μ u_τ⁴/(κ²y²). They cancel only if

**κ² = σ_ε √C_μ (C_ε2 − C_ε1).**   (7)

With the standard constants, κ = √(1.3 × 0.3 × 0.48) = **0.4327**, not 0.41. The k–ε model's own log law has a slightly smaller slope than the measured one; with κ = 0.41 the ε equation is out of balance by −2.6 % of its sink term. Reproducing κ = 0.41 would need C_ε2 = 1.871 or σ_ε = 1.167, which would spoil the decay of grid turbulence that fixed C_ε2. The constants are a compromise between flows, which is the nature of a calibrated model.

**Wall functions.** The standard k–ε model is not valid in the viscous sublayer. A CFD code using it puts the first cell centre in the log layer (30 < y⁺ < 0.1–0.2 Re_τ) and imposes the log-law values there (`@leishman-intro-aerospace-flight-vehicles` "Turbulent Flows": Turbulence Models). A low-Reynolds-number model instead needs y⁺ ≈ 1 at the first cell.

### 2.7 Where the model fails

Algebraic mixing length has no memory: l depends on where you are, not on what happened upstream. It cannot represent a separated flow, a boundary layer recovering after a step, or the effect of streamline curvature or rotation. The Boussinesq hypothesis also assumes the Reynolds stresses are aligned with the mean strain, which fails in secondary flows (the corner vortices in a square duct are produced by stress anisotropy and are invisible to any eddy-viscosity model). The friction it predicts in a plane channel is only as good as the calibration flows it was fitted to, which this lesson measures in part 4.

## 3. Build the computational model

Open `starter/main.py`. The stretched grid `channel_grid`, the cumulative trapezoid `cumtrapz0`, `log_law`, `dean_cf`, the constants and `main` are given.

### Step 1: the closure

`total_stress` is (2). `mixing_length` returns l⁺ with van Driest damping and the chosen outer form. `velocity_gradient` is (4), vectorised.

### Step 2: the profile and its integrals

`mixing_length_profile(re_tau)` builds the grid, evaluates τ⁺, l⁺ and g, and integrates:

```python
yp = channel_grid(re_tau, n, beta)
g = velocity_gradient(total_stress(yp, re_tau), mixing_length(yp, re_tau, kappa, A_plus, outer))
u_plus = cumtrapz0(g, yp)
```

It also returns ν_t⁺ = l⁺²g and the Reynolds stress ν_t⁺g. `bulk_velocity` and `friction` are (5); `laminar_profile` is the κ → 0 limit; `fit_log_law` is a least-squares fit of u⁺ against ln y⁺ for 30 ≤ y⁺ ≤ 0.2Re_τ (`P01.17`); `log_law_friction` is 2.5.

### Step 3: a duct in SI units

`re_tau_for_re_b` inverts the model with `brentq` (`P03.02`): Re_b grows monotonically with Re_τ. `duct_flow(V, gap)` computes Re_b, finds Re_τ, and converts back: u_τ = Re_τν/h, τ_w = ρu_τ², dP/dx = −τ_w/h, and the first-cell heights ν/u_τ (y⁺ = 1) and 30ν/u_τ.

### Step 4: k–ε in the log layer

`keps_log_layer` returns the profiles of 2.6. `keps_residuals` substitutes them into (6), with the diffusion terms by a conservative finite difference (fluxes at mid-points, then their difference), and returns both residuals normalised by their sink terms. `kappa_from_keps` is (7).

### Verify

`pytest curriculum/phase-08/lesson-15-turbulence-modeling` runs 13 tests in under 2 s, listed in the metadata. Two are worth reading: the κ → 0 limit, which recovers Poiseuille flow to 10⁻⁸ and c_f Re_b = 12, and the ε-equation test, which shows the residual falling at order 2 with the grid when κ = 0.4327 (it is truncation error) and tending to the constant −0.0255 when κ = 0.41 (it is a model inconsistency). On a geometric grid the conservative difference happens to be exact for ε ∝ 1/y, and the residual is round-off; the test uses a uniform grid to see the order.

### Compare

`main()` prints the duct, a sweep, and the k–ε check:

```
  Re_b = 16630, Re_tau = 488.4, u_tau = 0.2937 m/s, tau_w = 0.1038 Pa, dp/dx = -4.154 Pa/m
  c_f = 0.00690 (model) vs 0.00643 (Dean correlation): +7.3 %
  with Nikuradse's outer mixing length: Re_tau = 452.6, tau_w = 0.0892 Pa, dp/dx = -3.567 Pa/m, c_f = 0.00593 (-7.8 %)
  Dean correlation: tau_w = 0.0967 Pa, dp/dx = -3.870 Pa/m
  first-cell height for y+ = 1: 51.2 um; for y+ = 30 (wall functions): 1.54 mm
  Re_tau     Re_b      c_f model   c_f Dean   diff %   kappa_fit  B_fit   U_b+   u+_cl
     180       5112   0.00992    0.00863   +14.9     0.3445   2.870  14.200  16.148
     550      19075   0.00665    0.00621    +7.1     0.3906   4.470  17.341  19.075
    1000      37783   0.00560    0.00524    +7.0     0.4037   4.896  18.892  20.576
    2000      82554   0.00470    0.00431    +9.0     0.4112   5.163  20.638  22.293
    5200     239232   0.00378    0.00330   +14.5     0.4155   5.349  23.003  24.640
  kappa implied by the standard constants: sqrt(1.3 sqrt(0.09) (1.92 - 1.44)) = 0.4327
```

and writes `turbulence_models.png`: the Re_τ = 1000 profile on a log axis with the log law and u⁺ = y⁺ (labelled as a model, not data), and c_f against Re_b for the model, the log-law friction law and Dean's correlation. The profile follows u⁺ = y⁺ to within 1 % for y⁺ < 3, crosses into the log region near y⁺ = 30 (u⁺ = 13.11 against the log law's 13.30), and the viscous and turbulent stresses are equal at y⁺ = 10.5, the middle of the buffer layer; the Reynolds stress peaks at 0.90 τ_w at y⁺ = 58. At Re_τ = 180 the fitted κ (0.34) and B (2.9) are meaningless: there is almost no log layer to fit, which is a warning about fitting a log law to low-Reynolds-number data.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Reynolds-number sweep.** Re_τ = 180, 550, 1000, 2000, 5200 (table in §3). Predict whether the model's error against Dean grows or shrinks with Re. *Finding:* it is smallest (7 %) in the middle of the range and grows at both ends (15 % at 180 and 5200). At low Re the log layer is short and the damping dominates; at high Re Prandtl's l = κy keeps growing to the centreline, the profile has no wake, U_b⁺ is too low and c_f too high. The log-law friction law tracks Dean within 0.7 % at Re_τ = 500 and 2 % at 1000, drifting to 5 % at 2000.
2. **Sensitivity to the constants** at Re_τ = 1000. Predict which matters most: κ, A⁺ or the outer form. *Finding:*

   | change | c_f | against Dean | fitted B |
   |---|---|---|---|
   | baseline (κ 0.41, A⁺ 26, Prandtl outer) | 0.00560 | +7.0 % | 4.90 |
   | κ = 0.38 | 0.00512 | −1.1 % | 4.69 |
   | κ = 0.42 | 0.00577 | +9.8 % | 4.96 |
   | A⁺ = 22 | 0.00611 | +15.5 % | 4.39 |
   | A⁺ = 30 | 0.00518 | −0.1 % | 5.31 |
   | Nikuradse outer | 0.00477 | −7.0 % | (wake; fit not meaningful) |

   A⁺ and the outer closure matter most, each about ±8 %: the damping sets B (the intercept of the log law) and the outer length sets the wake. A ±0.01 change in κ moves c_f by about 3 %. This is the model's real uncertainty, and it is larger than Dean's scatter.
3. **The duct.** With both outer closures, the 50 mm duct at 5 m/s needs dP/dx = −4.15 or −3.57 Pa/m against Dean's −3.87 Pa/m: **3.9 ± 0.3 Pa per metre**, about 39 Pa over 10 m of duct, before bends, inlets and outlets (`P08.07`). Laminar flow at the same Re_b would need ten times less, which is why ducts are not sized with laminar formulas.
4. **k–ε consistency.** Compute κ from (7) for the standard constants and for one modified set. Predict the sign of the ε residual when κ = 0.41 is imposed. *Finding:* κ = 0.4327; with κ = 0.41 the normalised residual is −0.0255 everywhere in the log layer (diffusion falls short of the source imbalance); C_ε2 = 1.871 or σ_ε = 1.167 would restore κ = 0.41.

Raw data format: CSV `re_tau, re_b, cf_model, cf_dean, kappa_fit, B_fit, Ub_plus, ucl_plus` (1), `case, cf, cf_over_dean, B_fit` (2), `closure, re_tau, tau_w_Pa, dpdx_Pa_per_m, y1_um` (3), plus `turbulence_models.png`.

## 5. Validate and improve

| quantity | exact or published | model | difference | cause |
|---|---|---|---|---|
| κ → 0 limit, c_f Re_b | 12 (Poiseuille, exact) | 12.000 | < 10⁻⁵ | quadrature only |
| sublayer, y⁺ < 3 | u⁺ = y⁺ (exact) | within 0.4 % | — | the damped l⁺² term |
| log-layer slope, 200 < y⁺ < 1000 | 1/κ (the closure) | within 1 % | — | τ⁺ < 1 away from the wall |
| U_b⁺ grid convergence | — | order 2.0 | 8 × 10⁻⁵ at n = 2001 | trapezoid rule |
| c_f, duct | Dean 0.00643 (correlation) | 0.00690 / 0.00593 | +7.3 % / −7.8 % | outer mixing length; A⁺ |
| k–ε log-layer κ | 0.41 (log-law fits) | 0.4327 from the constants | +5.5 % | constants calibrated on several flows |

The model is verified (its limits are exact and its numerics converge) and partly validated: against a correlation, not against a measurement, and the correlation carries a few per cent of its own scatter. A comparison with the Lee & Moser DNS profiles (Extension 1) is the next validation step; their mean profiles are free to download with a citation request but no license, so nothing from them is shipped here.

What to improve: the outer region. A mixing length capped in the core (Nikuradse) and a damping constant tuned to the DNS would bring c_f within a few per cent over the sweep, but every such tuning is calibration on the same flow it is then tested on. What not to improve: the grid (U_b⁺ is converged to 10⁻⁵ at n = 2001).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2), (4) and (7); `turbulence_models.png`; the sweep and sensitivity tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the duct's pressure gradient as a band with its causes, and the first-cell height you would specify for a CFD mesh of the duct with a low-Re model and with wall functions.

### Adversarial self-checks

- **What if the velocity goes to zero?** Re_b → 0 and the flow becomes laminar; the mixing-length model does not know about transition, and at Re_τ ≈ 60 it still predicts a turbulent-looking profile. Below Re_b ≈ 2000–3000 use c_f = 12/Re_b.
- **What if the duct is twice as wide at the same speed?** Re_b doubles; c_f falls by about 2^(−1/4) = 16 %, and the pressure gradient per metre falls by more than half (τ_w/h with h doubled).
- **Which assumption fails first in a real duct?** Full development: a duct needs tens of hydraulic diameters to develop, and cabin ducts have bends every few metres. Then the side walls (a 50 mm gap in a 150 mm wide duct is not an infinite channel) and the secondary corner flows no eddy-viscosity model predicts.
- **What could falsify the model?** A measured pressure drop in a long straight duct more than about 15 % from the model's band, or a measured velocity profile without a log region at Re_τ ≈ 500.
- **Does the result survive a change of grid or method?** U_b⁺ changes by 8 × 10⁻⁵ between n = 1001 and 2001; the log-law friction law, an independent route, gives c_f within 1 % of Dean at this Re.
- **Which decision changes if the uncertainty doubles?** A fan sized with a 30 % margin on a ±8 % pressure-drop band is safe; at ±16 % the margin on the duct friction alone would be half used, and you would measure the pressure drop of a prototype duct before freezing the fan.

## Check yourself

1. What new term appears when the Navier–Stokes equations are averaged, and why does it make the system unclosed?
2. Why is the total shear stress across a fully developed channel linear in y whatever the turbulence model?
3. Write the mixing-length eddy viscosity and say why l = κy near a wall.
4. What does van Driest damping change, and which constant of the log law does it set?
5. Why is g = 2τ/(1 + √(1 + 4l²τ)) better than the textbook root of the quadratic?
6. In the log layer, what are k, ε and ν_t in terms of u_τ, κ and y in the k–ε model?
7. With C_μ = 0.09, C_ε1 = 1.44, C_ε2 = 1.92, σ_ε = 1.3, what κ does the k–ε model imply?
8. Where should the first grid cell be for a low-Reynolds-number model, and where for wall functions?
