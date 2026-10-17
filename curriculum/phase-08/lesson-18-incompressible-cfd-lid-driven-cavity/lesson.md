# P08.18  Incompressible CFD from Scratch

> A CFD code earns trust in three steps: it solves its own equations (verification), it converges at the order it was designed for (grid convergence), and it agrees with a trusted answer to a standard problem (validation). The lid-driven cavity is where every incompressible code takes that test.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the stream function–vorticity form of the 2-D incompressible Navier–Stokes equations, explain why it removes the pressure and satisfies continuity identically, and derive Thom's wall-vorticity condition from a Taylor series. 2. Implement, from scratch with NumPy and SciPy, a second-order finite-difference solver for the steady lid-driven cavity by two independent routes (Newton's method on the coupled sparse system, and explicit pseudo-time marching) and show they reach the same discrete solution. 3. Verify each operator against an exact or manufactured solution, show the discrete divergence vanishes to round-off, and run a three-grid convergence study with observed order, Richardson extrapolation and a grid-convergence index. 4. Validate the Re = 100 centreline velocities and primary vortex against the published data of Ghia, Ghia & Shin (1982), and explain why the difference stops shrinking on fine grids. 5. Decide what grid a stated accuracy needs, and what the solver cannot be trusted for. |
| **Prerequisites** | `P08.06` (the Navier–Stokes equations), `P08.09` (vorticity, circulation, the stream function), `P08.02` (stream function and streamlines), `P03.12` (finite-volume methods: conservation on a grid), `P03.11` (finite differences, stability limits, observed order), `P03.08`–`P03.09` (sparse direct and iterative solvers), `P03.10` (Newton's method for systems), `P03.18` (verification and validation), `P02.11` (PDEs), `P01.12` (divergence and Stokes's theorem). Skill check: you can write the five-point Laplacian and assemble it as a sparse matrix with `scipy.sparse.kron`. |
| **Engineering question** | Before a home-grown Navier–Stokes code is used on a cabin-air or wheel-well problem, can it reproduce the Re = 100 lid-driven cavity: what grid puts the centreline velocities within 1 % of the lid speed of the published data, at what order does it converge, and does it conserve mass? The answers: a 33 × 33 grid is within 0.9 % (largest differences 0.0089 in u and 0.0073 in v); the centre velocity converges at observed order 2.07 on 33/65/129 grids with a grid-convergence index of 0.28 %; the primary-vortex stream function extrapolates to −0.10352 against Ghia's −0.103423 (0.09 %); the discrete divergence is below 10⁻¹⁴. |
| **Assumptions and limits** | Two-dimensional, steady, incompressible, Newtonian, constant viscosity; a square cavity with a uniform lid speed (which is discontinuous at the two top corners, where the vorticity is singular); second-order central differences on a uniform collocated grid; Thom's first-order wall formula. Central convection differences lose robustness when the cell Reynolds number U h Re exceeds about 2: at Re = 1000 a 33 × 33 grid has no converged steady solution. The flow at Re = 100 is steady; above roughly Re = 8000 the real cavity flow is not, and a steady solver is then the wrong tool. |
| **Analytical work** | Curl of the momentum equation to the vorticity transport equation; ψ from continuity; Poisson equation ∇²ψ = −ω; Thom's formula; the commuting-difference proof of discrete continuity; the Jacobian of the discrete system for Newton; the stability limits of forward-Euler marching; observed order and Richardson extrapolation from three grids; the grid-convergence index. |
| **Computational deliverable** | `starter/main.py`: `laplacian_matrix`, `solve_poisson`, `velocities`, `wall_vorticity`, `divergence`, `convection`, `cavity_residual`, `solve_cavity_newton`, `march_cavity`, `compare_with_ghia`, `vortex_centre` (`grid`, `observed_order`, `difference_operators`, `centreline_profiles`, `centre_values`, the Ghia constants and `main` are given). Tests (15, under 4 s) prove: the Laplacian is symmetric negative definite with the textbook extreme eigenvalues; the Poisson solver converges at order 2.0 to a manufactured solution and is exact for a harmonic quadratic; velocities are exact on quadratics and divergence-free to round-off; Thom's formula is exact for the lid term and first order at a wall; the convection operator is exact on quadratics; Newton converges in at most 6 iterations and satisfies the discrete equations to 10⁻⁹; marching and Newton agree to 10⁻⁷; the net flux through the centreline falls as h²; the Stokes limit is mirror-symmetric and the Re = 100 vortex moves downstream; Reynolds similarity holds; the observed order on 17/33/65 grids is 1.7–2.3; the centreline profiles approach the published data. |
| **Practical activity** | **computational experiment**: a four-grid convergence study at Re = 100 (17 to 129), a comparison of the two solution routes in cost, a validation against published centreline data, and a Reynolds-number sweep by continuation. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `cavity.png`, the derivations of the ψ–ω equations and of Thom's formula, the grid-convergence table with observed orders, extrapolated values and GCI, the Ghia comparison table, and a half-page verification-and-validation note. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the ψ–ω equations and Thom's formula are derived; all public tests pass; the convergence study reports observed orders with the grids used and separates discretisation error from the published data's own error; the note states what grid meets a 1 % criterion and what the solver must not be used for. |
| **References** | `@barba-cfd-python-12-steps` "Step 10: 2D Poisson Equation" and "Step 11: Cavity Flow with Navier–Stokes" (Discretized equations) (required-free: the same cavity solved by a pressure-Poisson route, to compare with this lesson's stream-function route). `@mit-ocw-2-29-numerical-fluid-mechanics` Lecture 19: Streamfunction-Vorticity Methods (scheme and boundary conditions); Projection Methods; Fractional Step Methods (required-free); Lecture 17 "Numerical Methods for the Navier-Stokes Equations": Choice of Variable Arrangement on the Grid; Calculation of the Pressure (supplementary-free). `@nasa-nparc-cfd-vv-tutorial` "Examining Spatial (Grid) Convergence" (required-free: observed order, Richardson extrapolation, GCI). The benchmark data are cited in §5. |
| **Extensions** | ★ Write the same solver in primitive variables with a projection (fractional-step) method on a staggered grid, and check that it reproduces this lesson's centreline profiles to the discretisation error. ★ Replace the uniform lid by the regularised lid u = 16x²(1 − x)² and show that the corner vorticity then converges instead of growing like 1/h. ★ Add second-order upwinding of the convection term and find the Reynolds number at which a 33 × 33 grid again has a steady solution; measure what the upwinding costs in accuracy at Re = 100. |
| **Time** | unmeasured |

## 1. The engineering problem

Computational fluid dynamics (CFD) is how most flow problems on an aircraft are now first answered: the airflow in a cabin, the recirculation in an open wheel well or weapons bay, the cooling flow round an avionics rack. Every one of those answers comes from a code that discretises the Navier–Stokes equations (`P08.06`) and solves a very large algebraic system. A CFD result is only as good as the evidence that the code solves the equations it claims to solve, on a grid fine enough for the question. That evidence is not a pretty picture. It is a set of numbers: residuals, convergence rates, conservation errors, and agreement with trusted answers to standard problems.

The **lid-driven cavity** is the standard problem for incompressible codes. A square box of side L is full of fluid; its top wall (the lid) slides sideways at speed U; the other three walls are still. The lid drags the fluid round in a large vortex, with small counter-rotating vortices in the bottom corners. The geometry is trivial, so all the difficulty is in the flow: a recirculating region, thin shear layers along the walls, and two corners where the moving lid meets a stationary wall and the vorticity is infinite. It is the incompressible analogue of the closed recirculating flows in a cavity under a moving stream, and its Reynolds number Re = UL/ν is the only parameter.

In 1982 Ghia, Ghia and Shin published high-quality solutions at Reynolds numbers from 100 to 10 000, with tables of the velocity along the two centrelines. Those tables became the benchmark that thousands of codes have been checked against. This lesson builds a cavity solver from scratch and takes it through the whole test: operator checks, two independent solution routes, mass conservation, a grid-convergence study, and the comparison with the published data. Success criteria, written first: residuals of the discrete equations below 10⁻⁹; the two solution routes agree to 10⁻⁷; the discrete divergence at round-off; an observed order between 1.7 and 2.3 for a second-order scheme; and centreline velocities within 1 % of the lid speed of Ghia's values on the coarsest grid that achieves it.

## 2. Understand and derive

### 2.1 Equations and assumptions

The flow is two-dimensional, steady and incompressible, with constant density ρ and viscosity μ (a Newtonian fluid). Scale lengths by L, velocities by U and pressure by ρU². Mass and momentum conservation in differential form (`P08.06`) are

∂u/∂x + ∂v/∂y = 0,   (1)

u ∂u/∂x + v ∂u/∂y = −∂p/∂x + (1/Re) ∇²u,   u ∂v/∂x + v ∂v/∂y = −∂p/∂y + (1/Re) ∇²v,   (2)

with Re = UL/ν. The stress model is Newtonian and laminar: the viscous terms are (1/Re)∇² of the velocity because the flow is incompressible.

**The stream function.** Any velocity field of the form u = ∂ψ/∂y, v = −∂ψ/∂x satisfies (1) identically, because the mixed derivatives are equal (`P08.02`). The lines ψ = constant are streamlines, and the volume flow (per unit depth) between two streamlines is the difference of their ψ values. The walls are streamlines, so ψ = 0 on all four walls.

**The vorticity.** Define ω = ∂v/∂x − ∂u/∂y (`P08.09`). Substituting the stream function gives the kinematic link

**∇²ψ = −ω.**   (3)

Take ∂/∂x of the v-equation minus ∂/∂y of the u-equation. The pressure terms cancel (∂²p/∂x∂y − ∂²p/∂y∂x = 0), and with (1) the convective terms collapse to the transport of ω:

**u ∂ω/∂x + v ∂ω/∂y = (1/Re) ∇²ω.**   (4)

Equations (3) and (4) are two scalar equations for two unknowns, ψ and ω. The pressure has gone, and continuity is automatic. That is the whole appeal of the formulation in two dimensions. The price is the boundary condition: the walls give two conditions on ψ (ψ = 0, and ∂ψ/∂n = wall tangential velocity) and none on ω.

### 2.2 Thom's wall-vorticity formula

Take the lid, y = 1, where ψ = 0 and u = ∂ψ/∂y = U. Along the wall ψ is constant, so ∂²ψ/∂x² = 0 and (3) gives ω = −∂²ψ/∂y² there. Expand ψ one grid step h into the fluid:

ψ(x, 1 − h) = ψ(x, 1) − h ∂ψ/∂y + (h²/2) ∂²ψ/∂y² − … = 0 − hU − (h²/2) ω_lid + O(h³).

Solve for the wall value:

**ω_lid = −2ψ₁/h² − 2U/h,  and on a fixed wall ω_w = −2ψ₁/h²,**   (5)

where ψ₁ is the value at the first interior node. The neglected term is O(h³)/h² = O(h): Thom's formula is first-order accurate at the wall. Experience and the tests below show the interior solution still converges at second order, because the error is confined to a layer one cell thick. At the top corners the lid's speed jumps from U to 0, so ω at the corner nodes is undefined: the corners are left out of every stencil, and the vorticity next to them grows like 1/h as the grid is refined (Problem 6).

### 2.3 The discrete equations

Use N × N nodes with spacing h = 1/(N − 1) and central differences:

- the five-point Laplacian (ψ_{i−1,j} + ψ_{i+1,j} + ψ_{i,j−1} + ψ_{i,j+1} − 4ψ_{ij})/h², error (h²/12)(ψ_xxxx + ψ_yyyy);
- u_ij = (ψ_{i,j+1} − ψ_{i,j−1})/(2h), v_ij = −(ψ_{i+1,j} − ψ_{i−1,j})/(2h), and the same centred differences for ω_x and ω_y.

**Discrete continuity.** With D_x and D_y the centred difference operators, the discrete divergence is D_x(D_yψ) − D_y(D_xψ). The two operators act on different indices, so they commute, and the divergence is zero to round-off at every interior node, for any ψ that vanishes on the walls. Mass is conserved exactly in this discrete sense. A weaker, integral check is not exact: the trapezoid rule for the net flow through the vertical centreline, ∫u dy, which must be zero, gives (ψ_{N−2} − ψ_1)/2 + h U/2, an O(h²) error (0.0088, 0.0019, 0.00045 on 17, 33 and 65 grids).

**The cell Reynolds number.** Central differencing of (4) keeps the discrete operator well behaved while U h Re ≲ 2. At Re = 100 the cell Reynolds number is 6.25, 3.1, 1.6 and 0.8 on the 17, 33, 65 and 129 grids. The coarse grids are outside the comfortable range at Re = 100 and still converge; at Re = 1000 a 33 × 33 grid (cell Reynolds number 31) has no converged steady solution (Problem 4).

### 2.4 Two routes to the steady solution

**Newton's method (`P03.10`).** Collect the interior values into z = [ψ, ω], 2(N − 2)² unknowns. The residuals of (3) and (4) at the interior nodes, with the wall ω from (5), are

F₁ = Lψ + ω,  F₂ = (D_yψ)(D_xω) − (D_xψ)(D_yω) − Lω/Re,

where L, D_x and D_y are sparse matrices on the full grid and the wall ω is a linear function of ψ (a matrix W and the lid vector). F is quadratic in z, so its Jacobian is exact and easy to assemble from the same matrices: ∂F₂/∂ψ contains diag(D_xω)D_y − diag(D_yω)D_x plus the wall terms, and ∂F₂/∂ω contains diag(D_yψ)D_x − diag(D_xψ)D_y − L/Re. Start from the Stokes solution (Re → 0, convection dropped, one linear solve) and iterate z ← z − J⁻¹F with a sparse direct solve (`P03.08`). Newton converges quadratically once close: four iterations at Re = 100 on every grid.

**Pseudo-time marching.** Add ∂ω/∂t to (4) and march forward in time from rest until nothing changes: solve (3) for ψ, set the wall ω by (5), and update the interior ω by forward Euler. Forward Euler with central differences is stable only if dt ≤ h²Re/4 (diffusion) and dt ≤ 2/(Re U²) (the convection–diffusion limit for central differences, `P03.11`). The code takes 0.8 of the smaller. The steady state of the march satisfies exactly the same discrete equations as Newton's solution, so the two must agree to the stopping tolerance. That is the independent check of the Newton Jacobian.

### 2.5 Grid convergence

A second-order scheme has error ≈ C h² on fine enough grids. From values f₁, f₂, f₃ on grids refined by r = 2 (`P03.11`, `@nasa-nparc-cfd-vv-tutorial` "Examining Spatial (Grid) Convergence"):

**p = ln[(f₁ − f₂)/(f₂ − f₃)]/ln r,  f_∞ ≈ f₃ + (f₃ − f₂)/(r^p − 1),  GCI = 1.25 |f₃ − f₂|/(|f₃|(r^p − 1)).**   (6)

The observed order p checks that the grids are in the asymptotic range; Richardson extrapolation estimates the grid-independent value; the grid-convergence index (GCI) is a banded estimate of the fine grid's relative error with a safety factor 1.25. The quantities used should be point values that need no interpolation: ψ and u at the centre node (0.5, 0.5), which is a node on every odd grid, and the minimum of ψ located by a local quadratic fit.

## 3. Build the computational model

Open `starter/main.py`. The grid, `observed_order`, `difference_operators` (the sparse D_x, D_y, L, the selection matrices and Thom's W and lid vector), `centreline_profiles`, `centre_values`, the published data and `main` are given.

### Step 1: the Poisson solver

`laplacian_matrix(N)` builds the interior five-point matrix as kron(T, I) + kron(I, T) with T = tridiag(1, −2, 1)/h². `solve_poisson(f, boundary)` moves known edge values to the right-hand side and calls `spsolve`.

### Step 2: the field operators

`velocities`, `wall_vorticity` (5), `divergence`, `convection` and `cavity_residual` are each a few lines of array slicing, for example

```python
u[1:-1, 1:-1] = (psi[1:-1, 2:] - psi[1:-1, :-2]) / (2 * h)
w[1:-1, -1] = -2.0 * psi[1:-1, -2] / h**2 - 2.0 * U / h      # Thom, lid
```

### Step 3: Newton

`solve_cavity_newton(N, Re)` uses the operators of `difference_operators` to write F and J as in 2.4:

```python
dpsi = P @ (sp.diags(Dx @ omf) @ Dy - sp.diags(Dy @ omf) @ Dx) @ S \
     + P @ (sp.diags(Dy @ psf) @ Dx - sp.diags(Dx @ psf) @ Dy) @ W
dom  = P @ (sp.diags(Dy @ psf) @ Dx - sp.diags(Dx @ psf) @ Dy) @ S
J = sp.bmat([[A11, I], [dpsi - P @ L @ W / Re, dom - P @ L @ S / Re]], format="csc")
```

Stop when the largest residual is below 10⁻¹⁰ or the step is at round-off. An optional `initial` solution allows continuation in Re.

### Step 4: marching, and the post-processing

`march_cavity` is 2.4 with one `splu` factorisation reused at every step. `compare_with_ghia` interpolates the centreline profiles to the published stations with a cubic spline; `vortex_centre` fits a quadratic surface to the 3 × 3 nodes round the smallest ψ.

### Verify

`pytest curriculum/phase-08/lesson-18-incompressible-cfd-lid-driven-cavity` runs 15 tests on grids no larger than 65 × 65 in under 4 s. The operator tests come first: the manufactured solution φ = sin πx sin 2πy + eˣ cos y converges at order 2.00 (the eˣ cos y part puts non-zero values on the edges); x² − y² + 3xy is reproduced to 10⁻¹²; Thom's error on ψ = x²(1 − x)²y²(1 − y)² falls 0.0151, 0.0077, 0.0039, 0.0019 on 17–129 grids, order 1.0; the divergence of the velocity of a random ψ with zero edges is round-off. Then the solvers and the physics, as listed in the metadata.

### Compare

Running `main()` prints (times on one laptop core; yours will differ):

```
    N   unknowns  Newton its  time (s)   psi_min     x_c     y_c    max|u-Ghia|  max|v-Ghia|  max|div|
   17      450       4       0.023   -0.091060  0.6371  0.7591    0.04180      0.05025    2.2e-16
   33     1922       4       0.051   -0.100427  0.6208  0.7432    0.00889      0.00728    8.9e-16
   65     7938       4       0.252   -0.102768  0.6170  0.7387    0.00234      0.00591    3.6e-15
  129    32258       4       2.017   -0.103337  0.6161  0.7377    0.00431      0.00842    7.1e-15
N = 17: marching (2463 steps of dt = 0.0160, 0.36 s) against Newton (4 iterations, 0.023 s): max|psi difference| = 5.2e-11
N = 129 residuals: kinematics 1.7e-12, vorticity transport 5.5e-12
Stokes limit (Re = 0.001, N = 65): psi_min -0.100021 at (0.5000, 0.7648); symmetry error max|psi(x) - psi(1-x)| = 5.6e-07
```

and writes `cavity.png`: the streamlines on the 129 × 129 grid with the two bottom corner vortices, and the two centreline profiles on all four grids over the published points.

The two routes agree to 5 × 10⁻¹¹, far below the discretisation error: the Newton Jacobian is right. Newton is faster here by an order of magnitude, and its cost grows with the sparse factorisation (0.25 s at 65, 2.0 s at 129); marching needs dt ∝ h² once diffusion limits it (8041 steps and 5.3 s at 65 × 65), so its step count grows like N² and each step is a Poisson solve on N² unknowns: the cost grows at least like N⁴. The `@barba-cfd-python-12-steps` "Step 11: Cavity Flow with Navier–Stokes" notebook solves the same cavity in primitive variables with a pressure-Poisson equation and explicit time stepping; comparing its centreline profile at the same grid with this lesson's is a check that two different formulations converge to the same flow, which is the point of Extension 1.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Grid convergence at Re = 100.** Grids 17, 33, 65 and 129. Predict the observed order of the centre values for a second-order scheme with a first-order wall condition. *Finding:*

   | quantity | 33 | 65 | 129 | order p | extrapolated | GCI (129) |
   |---|---|---|---|---|---|---|
   | u(0.5, 0.5) | −0.20099 | −0.20719 | −0.20867 | 2.07 | −0.20913 | 0.28 % |
   | ψ(0.5, 0.5) | −0.062840 | −0.065646 | −0.066326 | 2.04 | −0.066544 | 0.41 % |
   | ψ_min | −0.100427 | −0.102768 | −0.103337 | 2.04 | −0.103519 | 0.22 % |

   On 17/33/65 the orders are 1.88 (ψ) and 2.02 (u). The first-order wall formula does not spoil the interior order.
2. **Validation against the published data.** The largest differences from Ghia's centreline tables, for each grid (table in §3). Predict whether the 129 × 129 grid, the same grid as Ghia's, matches them best. *Finding:* no. The differences fall from 0.042 to 0.0089 to 0.0023 (u) between 17 and 65, then rise to 0.0043 at 129; Richardson-extrapolated profiles differ from the table by up to 0.0051 (u) and 0.0093 (v). Part 5 explains why.
3. **Cost of the two routes.** Time Newton and marching on 17, 33 and 65 grids. Predict how the marching step count scales. *Finding:* 2463 steps at 17 and 33 (both limited by the convection–diffusion bound dt ≤ 2/Re, so dt = 0.016), 8041 at 65 (now limited by diffusion, dt = 0.0049 ∝ h²); Newton takes four iterations on every grid.
4. **Reynolds-number sweep by continuation.** Re = 100, 400, 1000 on the 65 × 65 grid, each started from the previous solution (Problem 4). Predict where the vortex centre moves. *Finding:* Newton from the Stokes start fails at Re = 1000 on the 65 grid but converges in 5–6 iterations by continuation; the centre moves from (0.617, 0.739) to (0.558, 0.610) to (0.535, 0.571), towards the cavity centre, and ψ_min is −0.1028, −0.1091, −0.1053. At Re = 1000 the 129 grid gives −0.1155: the 65 grid is no longer converged, which the study at one Re could not have told you.

Raw data format: CSV `N, unknowns, iterations, seconds, psi_min, x_c, y_c, max_du, max_dv, max_div` (1, 2), `N, steps, dt, seconds` (3), `Re, N, iterations, psi_min, x_c, y_c` (4), plus `cavity.png`.

## 5. Validate and improve

The benchmark is published data from Ghia, Ghia & Shin (1982), J. Comput. Phys. 48, 387–411, `@ghia-ghia-shin-1982-cavity`: Table I (u on the vertical centreline) and Table II (v on the horizontal centreline) for Re = 100, and the primary vortex in Table V (ψ = −0.103423 at (0.6172, 0.7344)). The article is paywalled; nothing in this lesson needs it, because the Re = 100 values are typed into the module as constants with that citation. They are a computed solution, not a measurement: Ghia et al. solved the same equations on a 129 × 129 grid with their own discretisation.

| quantity | this code | published (Ghia et al. 1982) | difference | cause |
|---|---|---|---|---|
| ψ_min, 129 × 129 | −0.103337 | −0.103423 | 8.6 × 10⁻⁵ (0.08 %) | both are discretisation errors of 129 grids |
| ψ_min, Richardson (33/65/129) | −0.103519 | −0.103423 | −9.6 × 10⁻⁵ | the table's own error |
| vortex centre, 129 × 129 | (0.6161, 0.7377) | (0.6172, 0.7344) | (−0.001, +0.003) | Ghia report a node of their grid (spacing 0.0078) |
| max |u − u_G|, 33 / 65 / 129 | 0.0089 / 0.0023 / 0.0043 | Table I | see text | |
| max |v − v_G|, 33 / 65 / 129 | 0.0073 / 0.0059 / 0.0084 | Table II | see text | |
| discrete divergence | < 10⁻¹⁴ | — | — | commuting differences |

The comparison is honest only with the convergence study beside it. This code's own error, from the GCI, is about 0.3 % on the 129 grid and falls as h². The difference from the published tables falls with h only until it reaches about 0.004–0.008 of the lid speed, and then stays there, with a smooth, systematic shape along the centrelines (positive near the lid, negative in the core). A transcription error would appear as a single outlier; a coding error would not converge at order 2 to a value that differs from the table by the same amount as the extrapolate does. The simplest explanation is that the benchmark table carries a discretisation error of about half a per cent of the lid speed. That is a general lesson: validation data have uncertainty too, and a code converging past them is not a code that has gone wrong. For the engineering question, the 33 × 33 grid already meets a 1 % criterion against the published table, and a 65 × 65 grid is needed for the code's own discretisation error (GCI) to be below about 1 %.

**What to improve.** The wall formula (first order) and the corner singularity limit accuracy near the lid; the near-lid differences in the Ghia comparison are the largest at every grid. A regularised lid (Extension 2) or a higher-order wall formula improves the near-wall error. For higher Re, the central convection scheme must be replaced or the grid refined until U h Re ≲ 2. **What not to improve:** the Newton tolerance (10⁻¹⁰ is far below the discretisation error) or the grid beyond 129 at Re = 100.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), (4) and (5); `cavity.png`; the grid-convergence table with orders, extrapolated values and GCI; the Ghia comparison table; and a half-page verification-and-validation note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that states the grid needed for 1 % agreement and the error budget of the result.

### Adversarial self-checks

- **What if Re goes to zero?** Convection vanishes and the problem is linear (Stokes flow); the solution must be mirror-symmetric about x = 0.5. At Re = 10⁻⁶ the code gives a symmetry error of 6 × 10⁻¹⁰ and a vortex centre at (0.5000, 0.7648) on the 65 grid; the test checks it.
- **What if the lid speed doubles?** In dimensionless form only Re = UL/ν matters: doubling U at fixed L and ν doubles Re. Doubling U while halving 1/ν must give exactly twice the ψ and ω; the test of Reynolds similarity checks it to 10⁻¹⁰.
- **Which assumption fails first?** The steady, two-dimensional flow: real cavities are three-dimensional (end-wall effects, Taylor–Görtler-like vortices along the lid) and become unsteady at high Re. Numerically, central convection fails first as Re grows on a fixed grid.
- **What could falsify the code?** An observed order far from 2 on refined grids; a divergence above round-off; marching and Newton disagreeing by more than the stopping tolerance; a Stokes solution that is not symmetric. Each is a test.
- **Does the result survive a change of solver or grid?** Two solvers agree to 5 × 10⁻¹¹; three grids give a consistent order and an extrapolate.
- **Which decision changes if the uncertainty doubles?** If the code's error were 0.6 % rather than 0.3 % on the 129 grid, it would be the same size as the differences from the published table, and you could no longer say which of the two is closer to the exact solution. You would refine once more before claiming anything finer than 1 %.

## Check yourself

1. Why does the stream function satisfy continuity automatically, and why is ψ constant on a wall?
2. Which term of the momentum equations disappears when you take the curl, and why does that help?
3. Write Thom's formula for the wall vorticity under a lid moving at U. What is its order of accuracy?
4. Why is the discrete divergence of (D_yψ, −D_xψ) zero to round-off on a uniform grid?
5. Three grids give 0.20099, 0.20719, 0.20867 (magnitudes). What is the observed order, and what does it tell you?
6. What is the cell Reynolds number, and why does it matter for central differences of the convection term?
7. Why can the difference between a converging code and a published benchmark stop decreasing on fine grids, and how do you tell that apart from a bug?
8. Name two reasons the lid-driven cavity is a good benchmark and one reason it is a poor model of a real wheel well.
