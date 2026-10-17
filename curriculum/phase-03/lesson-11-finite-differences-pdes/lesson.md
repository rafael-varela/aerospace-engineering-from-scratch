# P03.11  Finite Differences for PDEs

> A finite-difference scheme is a recipe for turning a partial differential equation into arithmetic; whether the arithmetic converges to the equation is decided by two numbers, the truncation error and the amplification factor, and you must compute both before you trust a single plot.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the FTCS, BTCS and Crank–Nicolson schemes for the one-dimensional heat equation from Taylor expansions, state the truncation error of each, and show that FTCS at r = 1/6 is fourth order in Δx. 2. Derive the von Neumann amplification factor of each scheme, prove the FTCS stability limit r ≤ ½, prove that the implicit schemes are stable for every r, and show that matrix (eigenvalue) stability analysis gives the same answer for Dirichlet boundaries. 3. Implement and test the three schemes and a Thomas tridiagonal solver from scratch, with the exact series, the exact semi-discrete solution, `scipy.linalg.solve_banded`, `scipy.linalg.expm` and `scipy.integrate.solve_ivp` as oracles, and measure the observed convergence orders in Δx and Δt. 4. Decide which scheme and time step to use for a transient heating problem by cost (steps × work per step) at a stated accuracy, and recognise the ringing that Crank–Nicolson produces on a step initial condition. |
| **Prerequisites** | `P02.11` (the heat equation, separation of variables, a first FTCS marcher and the Fourier number), `P03.05` (local and global error, order of accuracy, stability of a time stepper), `P03.06` (implicit methods for stiff systems; the semi-discrete heat equation is the standard stiff system), `P03.08` (banded elimination), `P01.16` (eigenvalues of a symmetric tridiagonal matrix), `P01.05` (Taylor series). Skill check: you can write the Taylor expansion of u(x + Δx, t) to fourth order and you know that |1 − z| ≤ 1 for z in [0, 2]. |
| **Engineering question** | A 10 mm titanium bulkhead separates a hot-gas duct from an avionics bay. At t = 0 the duct-side face jumps to 700 K while the bay-side face is held at 300 K by a coolant loop. How long until the mid-plane, where a bonded bracket must stay below 450 K, reaches that limit, and what heat load reaches the coolant at that moment? And which of three finite-difference schemes answers this in the fewest arithmetic operations at an accuracy of 0.01 s? |
| **Assumptions and limits** | One-dimensional conduction (the plate is wide compared with its thickness); constant k, ρ, c_p (Ti-6Al-4V's conductivity rises by roughly half between 300 and 700 K, which the problem set addresses); both faces at prescribed temperatures (Dirichlet), which idealises a hot-gas film and a coolant film as perfect; no contact resistance at the bracket; no radiation. The schemes are for a uniform grid; FTCS is stable only for r = αΔt/Δx² ≤ ½; Crank–Nicolson is stable for every r but at large r it does not damp grid-scale content, so a step initial condition rings for hundreds of steps. The series solution is exact only for the constant-property Dirichlet problem. |
| **Analytical work** | Taylor-series derivation of the three stencils and their truncation errors; the r = 1/6 cancellation; von Neumann analysis of all three schemes; the eigenvalues of the discrete Laplacian with Dirichlet boundaries and the equivalence of matrix and Fourier stability; the exact semi-discrete solution; the Fourier coefficients of the step problem and the cold-face heat flux; the cost comparison. |
| **Computational deliverable** | `starter/main.py`: `thermal_diffusivity`, `heat_exact_mode`, `heat_exact_step`, `cold_face_flux_exact`, `semi_discrete_exact`, `second_difference`, `solve_tridiagonal`, `step_ftcs`, `step_btcs`, `step_crank_nicolson`, `march`, `amplification_factor`, `max_stable_dt`, `update_matrix`, `spectral_radius`, `observed_order`, `spatial_convergence`, `temporal_convergence`, `time_to_limit`, `time_to_limit_exact` (`titanium_bulkhead`, `method_of_lines_oracle` and `main` are given). Tests prove: both exact solutions satisfy the PDE to a relative residual of 10⁻³ and their boundary, initial and steady limits; the second difference is exact for quadratics and converges at order 2.00 ± 0.05; the Thomas solver matches `solve_banded` to 10⁻¹²; each step holds the Dirichlet data, is the identity at r = 0, and the implicit steps satisfy their defining equations to 10⁻⁹; the update-matrix eigenvalues equal the von Neumann factors at θ_j = jπ/(n + 1) to 10⁻⁹ and FTCS is stable exactly when r ≤ ½; FTCS at r = 0.6 grows past 10⁶ in 200 steps and `march` refuses it; BTCS keeps a step monotone at r = 100 while Crank–Nicolson rings by tens of kelvin; the semi-discrete solution equals `expm` of the discrete operator and `solve_ivp`'s BDF result; the observed spatial order is 2.0 ± 0.1 for all three schemes at fixed r and 4.0 ± 0.1 for FTCS at r = 1/6; the observed temporal order is 1.0 ± 0.1 for BTCS and 2.0 ± 0.1 for Crank–Nicolson; and the time to the mid-plane limit agrees across schemes and with the series bisection. |
| **Practical activity** | **computational experiment**: a stability sweep across the r = ½ boundary, a two-axis convergence study (Δx at fixed r, Δt at fixed grid), a cost-at-fixed-accuracy comparison of the three schemes on the bulkhead, and a study of Crank–Nicolson ringing against the time step. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `heat_fd_schemes.png`, the hand derivations of the three truncation errors and three amplification factors, the stability and convergence tables, and a half-page note recommending a scheme and time step for the bulkhead with the cost and the accuracy it buys. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the truncation error and amplification factor of each scheme are derived; all public tests pass; the observed orders in both tables are explained from the derivations; the note states steps, operations and an error for each scheme and picks one with a reason. |
| **References** | `@ocw-16-90` Unit 2 "Numerical Methods for PDEs": "2.2 Partial Differential Equations", "2.3 Introduction to Finite Difference Methods", "2.4 Analysis of Finite Difference Methods", "2.7 Eigenvalue Stability of Finite Difference Methods" (required-free: the model PDEs, consistency, truncation error, von Neumann and eigenvalue stability); "2.6 Upwinding and the CFL Condition" (supplementary-free: the convection analogue, taken up in `P03.12`). `@ocw-18-085` "Lecture 23: Laplace's Equation", "Lecture 25: Fast Poisson Solver (part 1)" (supplementary-free: the eigenvalues of the discrete second-difference matrix used in §2.5). |
| **Extensions** | ★ Implement the DuFort–Frankel scheme, show it is unconditionally stable but inconsistent unless Δt/Δx → 0, and watch it converge to the wrong equation. ★ Replace the hot-face Dirichlet condition by a convective one, −k T_x = h(T_gas − T), with a ghost node, and redo the stability analysis. ★ Make k(T) linear in temperature and march the nonlinear problem with Crank–Nicolson and one Newton iteration per step (`P03.10`). |
| **Time** | unmeasured |

## 1. The engineering problem

A hot-gas duct passes through a bay full of avionics. The bulkhead between them is a 10 mm plate of Ti-6Al-4V (k = 7 W/(m·K), ρ = 4430 kg/m³, c_p = 560 J/(kg·K)). When the duct comes alive the duct-side face goes to 700 K almost at once; the bay-side face is bonded to a cold plate and stays at 300 K. Half-way through the thickness a composite bracket is bonded, and its adhesive is rated to 450 K. The steady state is a straight line from 700 K to 300 K, which puts the mid-plane at 500 K, so the bracket *will* exceed its limit; the question is **when**, because the duct runs for a known time, and what heat flux the coolant must carry at that moment. The answer should be good to 0.01 s, far better than the physics deserves, because we want to see the numerical error separately from the modelling error.

The problem is the one-dimensional heat equation of `P02.11`, T_t = αT_xx, with α = k/(ρc_p) = 2.822 × 10⁻⁶ m²/s and a diffusion time L²/α = 35.4 s. It has a series solution, so every number a scheme produces can be checked. That is the point: this lesson is not about the bulkhead but about **how to discretise a PDE and know that the discretisation is right**. The method is the one `@ocw-16-90` Unit 2 lays out: write the stencil, expand in Taylor series to find what equation it actually solves (consistency and truncation error), ask whether errors grow or decay from step to step (stability), and only then run it. A scheme can be consistent and useless (FTCS with too large a step), stable and wrong (Crank–Nicolson ringing on a step), or stable and consistent and simply expensive (FTCS with its tiny step). Each of those appears in this lesson with numbers attached.

Success criteria, written down first: all three schemes reproduce the series to better than 0.01 K at 5 s when their steps are small; the observed orders match the derived ones to 0.1; the spectral radius of the update matrix equals the largest von Neumann factor to 10⁻⁹; the time to the 450 K limit agrees between the series and the schemes to 5 ms where the derivation says it should; and the cost comparison names the number of steps and the arithmetic per step for each scheme at that accuracy.

## 2. Understand and derive

### 2.1 The grid, the stencil and the diffusion number

Put n + 1 nodes x_j = jΔx, j = 0 … n, Δx = L/n, across the plate and times t^m = mΔt. Write u_j^m ≈ T(x_j, t^m). Nodes 0 and n carry the prescribed face temperatures and never change; the unknowns are the n − 1 interior values. The **centred second difference**

D₂u_j = (u_{j−1} − 2u_j + u_{j+1})/Δx²   (1)

is the discrete T_xx. Expanding u_{j±1} = u_j ± Δx u_x + ½Δx² u_xx ± ⅙Δx³ u_xxx + (1/24)Δx⁴ u_xxxx ± …, the odd terms cancel and D₂u_j = u_xx + (Δx²/12) u_xxxx + O(Δx⁴): the second difference is **second-order accurate** and exact for any polynomial of degree three or less. `second_difference` implements (1), and the test measures the order 2.00 on sin 3x + x³.

For the time direction there are three natural choices, all built from D₂ and the dimensionless **diffusion number**

**r = αΔt/Δx²**,   (2)

the ratio of the time step to the time α⁻¹Δx² diffusion needs to cross one cell.

### 2.2 Three schemes and what equation each one solves

**FTCS** (forward time, centred space): evaluate D₂ at the old time.

u_j^{m+1} = u_j^m + r(u_{j−1}^m − 2u_j^m + u_{j+1}^m).   (3)

Explicit: each new value is an arithmetic combination of three old ones, no solve. Expand everything about (x_j, t^m): the left side is u + Δt u_t + ½Δt² u_tt + …, the right is u + αΔt(u_xx + (Δx²/12)u_xxxx + …). Divide by Δt and use u_t = αu_xx: the scheme satisfies u_t − αu_xx = −½Δt u_tt + (αΔx²/12) u_xxxx + O(Δt², Δx⁴). For a solution of the heat equation u_tt = α² u_xxxx, so the truncation error is

τ_FTCS = αΔx² (1/12 − r/2) u_xxxx + O(Δt², Δx⁴).   (4)

It is first order in Δt and second in Δx, and the two leading terms **cancel at r = 1/6**, leaving a fourth-order scheme. That is not a trick to rely on in general (it needs constant α and a uniform grid) but it is a sharp test of the derivation: the code's observed order jumps from 2 to 4 when r is set to 1/6, and nothing else changes.

**BTCS** (backward time; implicit Euler in time): evaluate D₂ at the new time.

u_j^{m+1} − r(u_{j−1}^{m+1} − 2u_j^{m+1} + u_{j+1}^{m+1}) = u_j^m.   (5)

The unknowns at the new level are coupled: (5) is a tridiagonal linear system with diagonal 1 + 2r and off-diagonals −r, solved once per step. The same expansion (now about t^{m+1}) gives τ_BTCS = αΔx²(1/12 + r/2) u_xxxx + …: first order in Δt, no cancellation.

**Crank–Nicolson**: average the two.

u_j^{m+1} − (r/2)(u_{j−1} − 2u_j + u_{j+1})^{m+1} = u_j^m + (r/2)(u_{j−1} − 2u_j + u_{j+1})^m.   (6)

This is the trapezoidal rule of `P03.05` applied to the semi-discrete system; expanding about the half-time t^{m+1/2} the odd terms in Δt cancel and τ_CN = (αΔx²/12) u_xxxx − (Δt²/12)(u_ttt − …) + …: **second order in both Δt and Δx**. The price is the same tridiagonal solve as BTCS, with diagonal 1 + r and off-diagonals −r/2 (`@ocw-16-90` 2.3, 2.4 give the general procedure for deriving these truncation errors).

### 2.3 Stability by von Neumann's method

Consistency says the scheme approaches the PDE as Δx, Δt → 0; it says nothing about whether errors already present grow. For a linear scheme with constant coefficients, substitute a single Fourier mode u_j^m = g^m e^{ijθ}, with θ = kΔx the phase advance per cell (θ runs from 0 to π; θ = π is the grid-scale sawtooth). The second difference of e^{ijθ} is (e^{−iθ} − 2 + e^{iθ})e^{ijθ}/Δx² = −(4 sin²(θ/2)/Δx²) e^{ijθ}. Writing s = sin²(θ/2) ∈ [0, 1]:

- FTCS: **g = 1 − 4rs**. |g| ≤ 1 for all s requires 4r ≤ 2, i.e. **r ≤ ½**. At r = 0.6 the sawtooth has g = −1.4 and grows by 1.4²⁰⁰ ≈ 10²⁹ in 200 steps, which the test confirms (the code reports > 10⁶; it is far more).
- BTCS: **g = 1/(1 + 4rs)**. 0 < g ≤ 1 for every r: **unconditionally stable**, and the sawtooth is damped almost completely (g(π) = 1/(1 + 4r), 0.0009 at r = 282).
- Crank–Nicolson: **g = (1 − 2rs)/(1 + 2rs)**. |g| ≤ 1 for every r because |1 − z| ≤ |1 + z| for z ≥ 0: unconditionally stable. But g(π) = (1 − 2r)/(1 + 2r) → −1 as r → ∞: at large r the sawtooth is **not damped**, it flips sign each step. Compare with the true factor e^{−αk²Δt} = e^{−4rs} for the same mode, which is essentially zero. Crank–Nicolson is stable but misrepresents the fast modes, and any initial condition with grid-scale content (a step) shows it.

The exact amplification per step is e^{−4rs}. FTCS's 1 − 4rs is below it (FTCS **over-damps** every mode, slightly), BTCS's 1/(1 + 4rs) is above it (BTCS **under-damps**, so a transient lingers and the bulkhead reaches its limit *late*), and Crank–Nicolson's factor matches to third order in rs. The temporal-convergence test checks the signs of these errors as well as their orders.

### 2.4 Stability by eigenvalues, and why the two agree here

Von Neumann's method assumes a periodic or infinite line. Our problem has Dirichlet ends. The honest analysis (`@ocw-16-90` 2.7) writes the update on the interior unknowns as a matrix: with A the (n − 1) × (n − 1) tridiagonal matrix of (1, −2, 1), FTCS is u^{m+1} = (I + rA)u^m, BTCS is (I − rA)u^{m+1} = u^m and Crank–Nicolson is (I − rA/2)u^{m+1} = (I + rA/2)u^m. The scheme is stable iff the **spectral radius** (largest |eigenvalue|) of the update matrix G is at most 1. The eigenvectors of A are the discrete sines sin(jθ_l) with θ_l = lπ/n, l = 1 … n − 1 (`@ocw-18-085` Lecture 23 derives this for the Laplacian), with eigenvalues −4 sin²(θ_l/2). Every G is a rational function of A, so its eigenvalues are g(θ_l): exactly the von Neumann factors, sampled at the n − 1 phases the grid can hold. For this problem the two analyses agree to rounding (the test asks for 10⁻⁹), and the table in `main()` shows ρ(G) = 1.038 for FTCS at r = 0.51 and 0.998 at r = 0.50: the stability boundary is sharp. The agreement fails for other boundary conditions and for variable coefficients, where the matrix analysis is the one to trust.

### 2.5 The exact solutions used as references

Three references, each isolating one error:

1. **The exact PDE solution** for one mode, sin(kx) e^{−αk²t} with k = mπ/L (`heat_exact_mode`). Comparing a scheme with it measures the total error.
2. **The exact semi-discrete solution** (`semi_discrete_exact`): the second difference turns sin(kx_j) into −λ_h sin(kx_j) with λ_h = (4/Δx²) sin²(kΔx/2), so the method-of-lines ODE system has the exact solution sin(kx_j) e^{−αλ_h t}. Comparing a scheme with it measures the **time** error alone; comparing it with reference 1 measures the **space** error alone, which is e^{−αk²t}(1 − e^{−α(λ_h − k²)t}) ≈ αk²t (kΔx)²/12 · e^{−αk²t}: 8.2 × 10⁻⁵ relative at n = 100 for mode 1 at t = 1/(αk²). The test checks this formula to 10⁻⁶. The same solution is `expm(αA t/Δx²)` applied to the initial data (`P02.05`), which is how the test ties it to SciPy.
3. **The step problem** (`heat_exact_step`), by the separation of variables of `P02.11`: T = T_hot + (T_cold − T_hot)x/L + Σ b_n sin(nπx/L) e^{−α(nπ/L)²t}, with b_n the sine coefficients of the initial deviation from the steady line, b_n = (2/(nπ))[(T₀ − T_hot)(1 − (−1)ⁿ) + (T_cold − T_hot)(−1)ⁿ]. For T₀ = T_cold = 300 K and T_hot = 700 K this is b_n = −800/(nπ). The cold-face flux is q = −k T_x(L) = k(T_hot − T_cold)/L − k Σ b_n (nπ/L)(−1)ⁿ e^{−α(nπ/L)²t} (`cold_face_flux_exact`): zero at t = 0 (the cold face has not yet felt anything), 143.0 kW/m² at 5 s, 280 kW/m² at steady state.

### 2.6 The cost of each scheme

FTCS costs about 5 flops per node per step; the Thomas algorithm costs about 8 per node per step (3 for elimination, 2 for back substitution, plus the right-hand side), so an implicit step is less than twice an explicit one. The decisive factor is the number of steps. On n = 100 cells, Δx = 0.1 mm, the FTCS limit is Δt ≤ Δx²/(2α) = 1.77 ms: reaching t ≈ 6 s takes 3665 steps at r = 0.45. BTCS and Crank–Nicolson can take Δt = 0.1 s (r = 28) and need 60 steps. Crank–Nicolson at Δt = 0.1 s puts the mid-plane at 450 K at **5.846 s** against the series' **5.845 s**; BTCS at the same step says 5.927 s, 82 ms late, because of the under-damping of §2.3; it needs Δt = 0.01 s (586 steps) to get within 8 ms. FTCS at 3665 steps is about 1 ms early. The cost ranking at 0.01 s accuracy is therefore Crank–Nicolson (59 steps × 8 flops/node), then FTCS (3665 × 5), then BTCS (586 × 8), and a problem with 10 000 cells instead of 100 would make FTCS 10⁴ times more expensive while leaving the implicit counts unchanged: the stability limit, not the accuracy, sets its cost.

### 2.7 Where the model stops being valid

- **Constant properties.** k of Ti-6Al-4V rises from about 7 to about 10–11 W/(m·K) between 300 and 700 K; α rises with it and the limit arrives sooner than the constant-α model says. Problem 6 asks by how much.
- **Dirichlet faces.** A real hot-gas film has a finite heat-transfer coefficient; the face temperature rises over seconds, not instantly, and the first-second behaviour of the solution (the Gibbs-like steep front) is an artefact of the idealised jump. The convective condition is extension 2.
- **One dimension.** Valid away from the plate edges and the bracket, which conducts heat sideways.
- **Uniform grid, linear PDE.** Von Neumann analysis needs both. The matrix analysis survives non-uniform grids; nothing here survives a nonlinear k(T) without a Newton step (extension 3).

## 3. Build the computational model

Open `starter/main.py`. `titanium_bulkhead`, `method_of_lines_oracle` and `main` are given; everything else is yours.

### Step 1: stencil and tridiagonal solver

```python
def second_difference(u, dx):
    return (u[:-2] - 2.0 * u[1:-1] + u[2:]) / dx ** 2                     # (1), interior nodes

def solve_tridiagonal(lower, diag, upper, rhs):
    b, d = diag.copy(), rhs.copy()
    for i in range(1, n):                                                 # forward elimination
        w = lower[i - 1] / b[i - 1]
        b[i] -= w * upper[i - 1]
        d[i] -= w * d[i - 1]
    x[-1] = d[-1] / b[-1]
    for i in range(n - 2, -1, -1):                                        # back substitution
        x[i] = (d[i] - upper[i] * x[i + 1]) / b[i]
```

No pivoting: the BTCS and Crank–Nicolson matrices are strictly diagonally dominant, so the pivots never vanish. Raise `ValueError` on a zero pivot anyway.

### Step 2: the three steps

`step_ftcs` is (3) on `u[1:-1]` with the ends copied. `step_btcs` builds the right-hand side `u[1:-1]` and adds `r * u[0]` to its first entry and `r * u[-1]` to its last (the known Dirichlet values move across the equals sign), then calls the Thomas solver with diagonal 1 + 2r and off-diagonals −r. `step_crank_nicolson` does the same with diagonal 1 + r, off-diagonals −r/2, and the right-hand side (1 − r)u_j + (r/2)(u_{j−1} + u_{j+1}) plus (r/2) times the boundary values. `march` dispatches, computes r, and refuses FTCS with r > ½: an unstable run is an error, not a result.

### Step 3: analysis functions

`amplification_factor(scheme, r, theta)` types in the three formulas of §2.3. `update_matrix` builds A with `np.eye(n, k=±1)` and returns I + rA, (I − rA)⁻¹ or (I − rA/2)⁻¹(I + rA/2); `spectral_radius` takes `np.linalg.eigvals`. `max_stable_dt` returns Δx²/(2α) for FTCS and `np.inf` otherwise.

### Step 4: convergence drivers and the engineering answer

`spatial_convergence(scheme, ns, L, alpha, t_final, r)` runs mode 1 on each grid with Δt = rΔx²/α, rounds t_final to a whole number of steps, and returns the max-norm error against `heat_exact_mode` *at the time actually reached*. `temporal_convergence` fixes the grid, requires each Δt to divide t_final, and measures against `semi_discrete_exact`. `observed_order` is log(e_i/e_{i+1})/log 2. `time_to_limit` marches the step problem from a uniform 300 K field with the faces set, watches the probe node, and interpolates linearly between the two steps that bracket T_limit; `time_to_limit_exact` bisects the series.

Running `main()` prints

```
Ti-6Al-4V bulkhead: alpha = 2.822e-06 m^2/s, diffusion time L^2/alpha = 35.4 s
n = 100 cells, dx = 0.10 mm: FTCS needs dt <= 1.772 ms; implicit schemes have no limit

exact (series) time for the mid-plane to reach 450 K: 5.845 s
  ftcs  dt =    1.595 ms  r =    0.450: t = 5.844 s (-1.1 ms) in 3665 steps
  btcs  dt =  100.000 ms  r =   28.217: t = 5.927 s (+81.5 ms) in 60 steps
  cn    dt =  100.000 ms  r =   28.217: t = 5.846 s (+0.2 ms) in 59 steps
  btcs  dt =   10.000 ms  r =    2.822: t = 5.854 s (+8.3 ms) in 586 steps
  cn    dt =   10.000 ms  r =    2.822: t = 5.846 s (+0.2 ms) in 585 steps
```

followed by the profile errors at 5 s, the spectral-radius table, the two convergence tables and the figure.

### Verify

- **Analytic cases:** both exact solutions satisfy T_t = αT_xx to a relative residual of 10⁻³ (central differences in t, (1) in x); the step series is 700 K and 300 K at the faces, 300 K inside at t = 0 (to 0.5 K with 20 000 terms), the straight line at t = 10⁴ s, and monotone in between; the flux is zero at t → 0 and k ΔT/L at steady state.
- **Limiting cases:** r = 0 is the identity for all three steps; all three agree with u + rD₂u to O(r²); the implicit steps satisfy (5) and (6) to 10⁻⁹; `update_matrix` satisfies its defining identities.
- **Stability:** eigenvalues of G equal g(θ_l) to 10⁻⁹ for four values of r; FTCS is stable iff r ≤ ½; the sawtooth grows past 10⁶ in 200 FTCS steps at r = 0.6 and `march` raises; BTCS keeps a 400 K step monotone at r = 100 while Crank–Nicolson's second difference exceeds 20 K.
- **Convergence (observed order):** space, fixed r: 2.0 ± 0.1 for FTCS (r = 0.4), BTCS (0.4) and Crank–Nicolson (0.5); 4.0 ± 0.1 for FTCS at r = 1/6. Time, fixed grid: 1.0 ± 0.1 for BTCS, 2.0 ± 0.1 for Crank–Nicolson, with the under-damping of BTCS and over-damping of FTCS checked by sign.
- **Agreement with the oracles:** Thomas = `solve_banded` to 10⁻¹²; semi-discrete exact = `expm` to 10⁻¹⁰ and = `solve_ivp` BDF to 10⁻⁷; the time to the limit from FTCS and Crank–Nicolson within 0.5 % of the series bisection and from BTCS within 3 %.

### Compare

`scipy.integrate.solve_ivp(method="BDF")` on the semi-discrete system is the library way to do what BTCS and Crank–Nicolson do by hand: an implicit multistep method with adaptive order and step, a Newton solve (trivial here, the system is linear) and error control to a requested tolerance. It returns the semi-discrete solution to 10⁻⁷ K and differs from the series by the 0.002 K spatial error of the grid, which it cannot see. What it does not give you is r, g(θ) or the knowledge that the mid-plane answer from BTCS at Δt = 0.1 s is 82 ms late for a reason you can derive. In practice the heat-conduction solver in a thermal code is Crank–Nicolson or BDF2 with a few fully implicit steps at the start, exactly the combination §4 arrives at, and `solve_banded` or a sparse LU does the tridiagonal solve. Call them; but set the step from r and the modes you need to resolve, not from a tolerance you have not derived.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before the run.

1. **The stability boundary.** On n = 40 cells, start from mode 1 plus 1 % of the sawtooth (−1)ʲ and run 200 FTCS steps at r = 0.40, 0.49, 0.50, 0.51, 0.55. Predict from g(π) = 1 − 4r which runs stay bounded and by what factor the sawtooth has grown or decayed after 200 steps (|1 − 4r|²⁰⁰). Record `r, max|u|, predicted factor`.
2. **Two-axis convergence.** Run `spatial_convergence` for the three schemes at r = 0.4 (and FTCS at r = 1/6) on n = 10, 20, 40, 80, and `temporal_convergence` for BTCS and Crank–Nicolson on n = 200 at Δt = t_f/5 … t_f/40 with t_f = L²/(απ²) = 3.59 s. Predict orders 2, 2, 2, 4 and 1, 2 before you run.
3. **Cost at fixed accuracy.** For the bulkhead, find for each scheme the largest Δt that puts the mid-plane time within 10 ms of the series value, and tabulate `scheme, dt, steps, flops ≈ steps × nodes × (5 or 8)`. Predict from §2.6 that Crank–Nicolson wins by a factor of about 40 over FTCS and about 10 over BTCS.
4. **Crank–Nicolson ringing.** March the step problem with Crank–Nicolson at Δt = 1, 0.5, 0.1, 0.05 s to t = 5 s and record the maximum error and the maximum |second difference| of the profile. Then repeat Δt = 0.5 s with two BTCS half-steps first and Crank–Nicolson after. Predict: the ringing amplitude scales with |g(π)|^{steps}, so it persists for hundreds of steps at r ≈ 140 and is removed by two fully implicit steps that multiply the sawtooth by (1/(1 + 4r))² ≈ 3 × 10⁻⁶.

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `heat_fd_schemes.png`.

What you should find: the r = 0.50 run holds its sawtooth at amplitude 1 (g = −1 exactly) and r = 0.51 grows it by 1.04²⁰⁰ ≈ 2500; the observed orders are 2.01, 2.00, 2.00 (FTCS 0.4), 1.98, 2.00, 2.00 (BTCS), 1.97, 1.99, 2.00 (Crank–Nicolson) and 4.01, 4.00, 4.00 (FTCS 1/6), then 0.94, 0.97, 0.99 (BTCS in time) and 2.00, 2.00, 2.00 (Crank–Nicolson in time); Crank–Nicolson at Δt = 0.5 s (r = 141) is wrong by 211 K at 5 s with a second difference of 325 K, while BTCS at the same step is wrong by 6.2 K with no ringing and Crank–Nicolson at 0.05 s by 0.010 K; two BTCS half-steps before Crank–Nicolson at Δt = 1 s bring the 5 s error from 330 K to 1.2 K.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| α, L²/α | 7/(4430 × 560) = 2.822 × 10⁻⁶ m²/s, 35.4 s | 2.822 × 10⁻⁶, 35.4 | 0 | definitions |
| FTCS limit, n = 100 | Δx²/(2α) = 1.772 ms | 1.772 ms | 0 | — |
| ρ(G) vs max|g|, r = 0.51 | equal | 1.037987 both | < 10⁻⁹ | §2.4, exact eigenvectors |
| spatial order, FTCS r = 0.4 / 1/6 | 2 / 4 | 2.00 / 4.00 | 0.01 | τ in (4) |
| temporal order, BTCS / CN | 1 / 2 | 0.94–0.99 / 2.00 | pre-asymptotic at Δt = t_f/5 | O(Δt²) term of BTCS visible at the coarsest step |
| spatial error, mode 1, n = 100 | αk²t(kΔx)²/12 · e⁻¹ ≈ 3.0 × 10⁻⁵ | 3.0 × 10⁻⁵ | 10⁻⁶ rel. | — |
| time to 450 K, CN Δt = 0.1 s | 5.845 s (series) | 5.846 s | 0.2 ms | O(Δt²) + linear interpolation |
| time to 450 K, BTCS Δt = 0.1 s | later than series (under-damping) | 5.927 s | +82 ms | O(Δt) |
| cold-face flux at 5 s | 142.99 kW/m² (series) | 143.02 (FTCS), 143.01 (CN 0.05 s) | 0.03 | one-sided difference at the face |
| profile at 5 s, CN Δt = 0.5 s | stable, |g(π)| = 0.993 | error 211 K | — | ringing, not instability |
| method of lines (BDF) vs series | spatial error only | 0.002 K | — | the grid's (kΔx)²/12 |

Every discrepancy has a derived cause and none is "numerical noise". The one that matters for the engineer is the last-but-two row: a stable scheme produced a profile wrong by 211 K, with no warning, because stability (|g| ≤ 1) and accuracy (g ≈ e^{−4rs}) are different requirements and only the first is a theorem. The remedy is cheap (a few BTCS steps at the start, or a smaller first step) and the lesson is general: for any problem with sharp initial data, look at g at θ = π before choosing Crank–Nicolson.

What agreement with the series does *not* validate: the bulkhead model. Real conductivity rises with temperature, the hot face is heated through a gas film, and the bracket conducts sideways. Those change the answer by seconds, not milliseconds, and no refinement of Δx or Δt can see them. The numerical study tells you the model's prediction to 1 ms; the model tells you the bracket's fate to perhaps ±2 s. State both.

Improvements worth making: Rannacher start-up (two BTCS half-steps) before Crank–Nicolson, which keeps second-order accuracy and removes the ringing; a convective hot-face condition (extension 2); k(T) with a Newton step (extension 3). Improvement worth *not* making: a finer grid to cure Crank–Nicolson ringing. At fixed Δt, halving Δx quadruples r and pushes g(π) closer to −1: the ringing gets worse.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivations of (4), its BTCS and Crank–Nicolson counterparts, and the three amplification factors with the stability conclusions; `heat_fd_schemes.png`; the four tables of §4; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a scheme and time step for the bulkhead with the steps, the arithmetic and the error each alternative would have cost.

### Adversarial self-checks

- **What if the plate were aluminium (α = 6.9 × 10⁻⁵ m²/s, 24 times larger)?** The diffusion time falls to 1.4 s, the FTCS step limit at Δx = 0.1 mm falls to 72 µs, and the implicit schemes' step is set by the time scale you want to resolve, so their advantage grows by the same factor of 24. Which number in your note changes, and which does not?
- **What if the grid is halved at fixed Δt?** r quadruples. FTCS becomes unstable if it was anywhere near r = ½; BTCS is unaffected; Crank–Nicolson rings more. Check with `amplification_factor("cn", 4 r, π)`.
- **Which assumption fails first in the real world?** The Dirichlet hot face. A gas-side heat-transfer coefficient of 500 W/(m²·K) gives a Biot number hL/k ≈ 0.7: the face temperature lags the gas by hundreds of kelvin in the first seconds, and the mid-plane reaches 450 K later than the model says by a margin larger than every numerical error in this lesson combined.
- **What measurement could falsify the model?** A thermocouple at the cold face: the model predicts the flux history of `cold_face_flux_exact` (zero for the first fraction of a second, 143 kW/m² at 5 s). A measured flux rising earlier means the plate is thinner or more conductive than assumed; a plateau below 280 kW/m² means the hot face never reached 700 K.
- **Does the result survive a change in grid, step or scheme?** The mid-plane time agrees to 1 ms across FTCS at r = 0.45, Crank–Nicolson at 0.1 and 0.01 s and the series; BTCS needs Δt ≤ 0.01 s to join them. That spread is the numerical uncertainty, and it is 1000 times smaller than the modelling uncertainty above.
- **Which decision changes if the uncertainty doubles?** If the adhesive limit were known only to ±20 K, the time window would be ±0.5 s, and the choice between Crank–Nicolson at 0.1 s and BTCS at 0.1 s (82 ms apart) would be immaterial; you would spend the effort on the gas-film model instead.

## Check yourself

1. Define the diffusion number r and state the stability limit of FTCS. For aluminium (α = 6.9 × 10⁻⁵ m²/s) on a 0.1 mm grid, what is the largest stable FTCS step?
2. Write the von Neumann amplification factors of FTCS, BTCS and Crank–Nicolson and the exact factor they approximate. Which one is below the exact value, and which above?
3. Why is Crank–Nicolson stable for every r and yet unreliable for a step initial condition at large r? Which value of g tells you?
4. What is special about r = 1/6 for FTCS, and why does it not make FTCS a fourth-order method in general?
5. How does the matrix (eigenvalue) stability analysis relate to the von Neumann analysis for Dirichlet boundaries, and when do they disagree?
6. Two references are used in the convergence tests: the exact PDE solution and the exact semi-discrete solution. Which error does each comparison isolate?
7. A BTCS run at Δt = 0.1 s puts the mid-plane at 450 K 82 ms late. Is the sign of that error predictable? Roughly how late would Δt = 0.05 s be?
8. The implicit step costs about 1.6 times the explicit one per node. Under what condition is the implicit scheme nevertheless cheaper, and by how much for the bulkhead on 100 cells?
