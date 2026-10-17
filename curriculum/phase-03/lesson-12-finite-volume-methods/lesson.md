# P03.12  Finite Volume Methods

> A finite-volume scheme never loses a gram: whatever leaves one cell enters the next. What it can lose is the *shape* of what it carries, and the flux formula decides whether the loss is smearing, wiggles, or neither.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the finite-volume update from the integral form of a conservation law, explain why it conserves the total exactly for any flux formula, and derive the upwind, Lax–Friedrichs, Lax–Wendroff and minmod-limited fluxes for linear advection. 2. Derive the von Neumann amplification factors of the linear schemes, the CFL condition |ν| ≤ 1, and the modified equation of the upwind scheme with its numerical diffusivity a Δx(1 − ν)/2. 3. Implement and test the four schemes from scratch in flux form, with the exact translated solution, `scipy.integrate.quad`, `scipy.linalg.circulant` and `scipy.integrate.solve_ivp` as oracles, and measure the observed convergence order on smooth data (1 and 2) and on a discontinuity (½ and about ⅔). 4. Decide which scheme a feed-line model should use when the quantity carried must stay bounded, and state the grid a first-order scheme would need before its numerical diffusion falls below the physical one. |
| **Prerequisites** | `P03.11` (finite differences, von Neumann analysis, observed order), `P03.05` (explicit time stepping and its stability), `P03.04` (Gauss–Legendre quadrature, used for cell averages), `P02.11` (the wave and advection equations, d'Alembert's translating solution), `P01.12` (the divergence theorem: the integral form of a conservation law), `P02.10` (Fourier modes on a periodic grid). Skill check: you can integrate a conservation law over a box in x and t and you know that e^{iθ} has unit modulus. |
| **Engineering question** | A slug of off-mixture propellant, 0.3 m long, is created by a valve transient in a feed line where the flow moves at 12 m/s. A mixture-ratio sensor 1.2 m downstream feeds the engine controller. When does the slug reach the sensor, how sharp is its front, and what peak does the sensor read? And when the model replaces the real line, which numerical scheme reports an arrival time, front steepness and peak that belong to the propellant rather than to the grid? |
| **Assumptions and limits** | One-dimensional plug flow at constant speed (no velocity profile across the pipe, no turbulent mixing: the real line has an eddy diffusivity of order 10⁻³ m²/s which this lesson omits and §5 returns to); the mixture fraction is a passive scalar, carried but not reacting; a periodic model domain of 2.4 m stands in for a long line (the slug never wraps during the runs). The exact solution is a pure translation, so every deformation is numerical. Explicit schemes on a three-cell stencil are stable only for |ν| = |a|Δt/Δx ≤ 1; the first-order schemes smear a front over a width that grows like √t; Lax–Wendroff produces undershoots (negative mixture fraction) at a sharp front; the minmod scheme is bounded but drops to first order at extrema. |
| **Analytical work** | The integral form of u_t + (au)_x = 0 over a cell and a time step; the update as an exact identity and the numerical flux as its only approximation; the four flux formulas; von Neumann factors and |g|² for the three linear schemes; the CFL condition; the modified equation of upwind by Taylor expansion and its diffusivity; the variance growth 2Dt of a Gaussian; the phase lag of Lax–Wendroff; the L1 convergence orders expected on smooth and discontinuous data. |
| **Computational deliverable** | `starter/main.py`: `cell_centres`, `cell_averages`, `gaussian_pulse`, `square_pulse`, `exact_advection`, `courant_number`, `minmod`, `numerical_flux`, `fv_step`, `advect`, `sensor_history`, `arrival_time`, `amplification_factor`, `phase_speed_ratio`, `numerical_diffusivity`, `pulse_variance`, `total_variation`, `l1_error`, `observed_order`, `convergence_study` (`feed_line`, `semi_discrete_upwind_oracle`, `main` are given). Tests prove: cell averages match `quad` and are exact for cubics; the exact solution translates and is periodic; each flux reduces to its textbook form, is consistent (uniform state → au), and Lax–Wendroff and Lax–Friedrichs collapse to upwind at ν = 1; every scheme conserves Σu_iΔx to 10⁻¹³, mirrors under a → −a, refuses |ν| > 1, and shifts the data by exactly one cell at ν = 1; upwind and minmod are monotone and total-variation-diminishing on a square slug while Lax–Wendroff overshoots by more than 10 % both ways; the von Neumann factors equal the eigenvalues of the circulant update matrix to 10⁻⁸ and exceed 1 past the CFL limit; the upwind scheme spreads a Gaussian's variance by exactly 2Dt to 0.2 % for three Courant numbers; the gap to `solve_ivp`'s semi-discrete solution is O(ν); observed L1 orders are 1 (from below) for upwind, 2.00 ± 0.05 for Lax–Wendroff and between 1.5 and 2 for minmod on the Gaussian, and 0.50 ± 0.05 for upwind on the square slug; the sensor reads the front's 50 % crossing within 4 ms for every scheme, but the peak is 0.83 (upwind), 0.57 (Lax–Friedrichs), 1.20 with a −0.20 undershoot (Lax–Wendroff) and 0.98 (minmod) against a true 1.00. |
| **Practical activity** | **computational experiment**: a Courant-number sweep for accuracy and stability (including ν = 1 exactness and ν > 1 growth), a convergence study on smooth and discontinuous data, a measured-against-predicted test of the modified equation, and a sensor-reading study on four grids. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `fv_advection.png`, the hand derivations of the integral form, the four fluxes, the three amplification factors and the upwind modified equation, the convergence and sensor tables, and a half-page note recommending a scheme and grid for the feed-line model with the error it leaves in the sensor prediction. |
| **Assessment** | `assessments/phase-03/grading-rubric.md`; competencies 1–4. Pass: the integral-form derivation and the modified equation are complete; all public tests pass; the observed orders on smooth and discontinuous data are explained; the note states a peak-reading error and a front-arrival error for the chosen scheme and grid. |
| **References** | `@ocw-16-90` Unit 2 "2.5 Introduction to Finite Volume Methods", "2.6 Upwinding and the CFL Condition" (required-free: conservation form, cell averages, fluxes, the upwind flux and the CFL condition); "2.3 Introduction to Finite Difference Methods", "2.4 Analysis of Finite Difference Methods" (supplementary-free: the Taylor-expansion machinery behind the modified equation). `@ocw-16-920j` topic "Finite Difference and Finite Volume Discretizations" (supplementary-free: the course page lists topics rather than note titles; the lecture notes for that topic cover the same conservation-form discretisation). |
| **Extensions** | ★ Replace the periodic boundary by an inflow condition u(0, t) = 0 and a first-order outflow extrapolation, and check that the total in the domain falls by exactly the outflow flux integrated in time. ★ Apply the minmod scheme to inviscid Burgers' equation u_t + (u²/2)_x = 0 with the Godunov flux and watch a shock form from the Gaussian slug; verify the Rankine–Hugoniot shock speed. ★ Add a physical diffusion term with an explicit central stencil and derive the combined stability condition ν + 2DΔt/Δx² ≤ 1. |
| **Time** | unmeasured |

## 1. The engineering problem

A valve in a liquid-propellant feed line closes a few milliseconds late and a slug of propellant with the wrong mixture ratio enters the line. The slug is about 0.3 m long and the flow carries it at 12 m/s toward the injector. A mixture-ratio sensor 1.2 m downstream watches for exactly this kind of event: if the slug's front arrives sharp and its peak is above a threshold the controller trims the other propellant; if it arrives smeared and low, nothing happens. The question for the line's designer is when the slug arrives, how steep its front is and how high its peak reads. The question for the engineer *modelling* the line is harder: if the simulation says the front arrives 10 ms early and the peak is 17 % low, is that the propellant or the grid?

The physics is the simplest conservation law there is, linear advection u_t + a u_x = 0, whose exact solution is the initial profile translated unchanged at speed a (`P02.11`). Every change of shape a scheme produces is therefore error, and can be measured. Yet this toy equation is where every CFD method is first tested, because the two ways a scheme can fail on it (smearing a front, or decorating it with wiggles) are the two ways schemes fail on the Euler equations of a nozzle or an inlet, and the cures (upwinding, flux limiting, the CFL condition) are the same. `@ocw-16-90` 2.5 introduces the method on exactly this equation.

Success criteria, written first: every scheme conserves the total slug mass to rounding; the exact Courant number ν = 1 reproduces the translation exactly; the measured spreading of a Gaussian under the upwind scheme matches the modified equation's 2Dt to 1 %; the observed L1 orders are 1 and 2 on smooth data and ½ for the first-order scheme at a jump; and the sensor study names the error each scheme puts into the arrival time and the peak, so that the designer's question can be separated from the modeller's.

## 2. Understand and derive

### 2.1 Conservation form and the integral law

Write the advection equation as u_t + f(u)_x = 0 with flux f(u) = au: the rate of change of u in any interval equals the flux in minus the flux out. Divide [0, L) into n cells of width Δx = L/n, cell i covering [x_{i−1/2}, x_{i+1/2}] with centre x_i, and define the **cell average** ū_i(t) = (1/Δx)∫ u dx over cell i. Integrate the PDE over the cell and over a time step from t^m to t^{m+1} (`P01.12` in one dimension):

**ū_i^{m+1} = ū_i^m − (Δt/Δx)(F_{i+1/2} − F_{i−1/2}),   F_{i+1/2} = (1/Δt)∫_{t^m}^{t^{m+1}} f(u(x_{i+1/2}, t)) dt.**   (1)

This is **exact**: no approximation has been made. F_{i+1/2} is the time-averaged flux through the interface between cells i and i + 1. The whole of finite-volume numerics is the approximation of that one number from the cell averages on either side, and whatever formula is used, (1) adds F_{i+1/2} to cell i + 1 and subtracts it from cell i, so **Σ_i ū_i Δx is conserved exactly** (to rounding) on a periodic domain. The test checks this to 10⁻¹³ for all four schemes; it is a property of the form, not of the accuracy.

Because the unknowns are averages, the initial data and the reference solution must be averaged too (`cell_averages`, by four-point Gauss–Legendre quadrature from `P03.04`). Sampling at cell centres differs from the average by Δx² u''/24, enough to spoil a second-order convergence study. The test checks the averages against `scipy.integrate.quad`.

### 2.2 Four numerical fluxes

Let a > 0 and ν = aΔt/Δx (the **Courant number**: the fraction of a cell the solution moves per step). The interface flux must be built from ū_i and ū_{i+1}.

**Upwind.** The information at the interface comes from where the flow comes from: F_{i+1/2} = a ū_i. Substituting in (1), ū_i^{m+1} = (1 − ν)ū_i + ν ū_{i−1}: for 0 ≤ ν ≤ 1 a convex combination of two old values, so the new value lies between them. The scheme can create no new maximum or minimum (**monotone**), and the total variation TV = Σ|ū_{i+1} − ū_i| never increases (**TVD**). At ν = 1 it is ū_i^{m+1} = ū_{i−1}: the data shift by exactly one cell, which is the exact solution.

**Lax–Friedrichs.** Average the two sides and add the largest diffusion a monotone scheme can carry: F = ½a(ū_i + ū_{i+1}) − (Δx/2Δt)(ū_{i+1} − ū_i). The update is ū_i^{m+1} = ½(ū_{i−1} + ū_{i+1}) − ½ν(ū_{i+1} − ū_{i−1}): monotone for |ν| ≤ 1, twice as diffusive as upwind at ν = ½, and it does not need to know the sign of a, which is why it survives in systems where the wave speeds have both signs.

**Lax–Wendroff.** Expand in time: u(t + Δt) = u + Δt u_t + ½Δt² u_tt + …, and use the PDE twice, u_t = −au_x, u_tt = a²u_xx. Replacing u_x and u_xx by centred differences gives ū_i^{m+1} = ū_i − ½ν(ū_{i+1} − ū_{i−1}) + ½ν²(ū_{i+1} − 2ū_i + ū_{i−1}), which is (1) with F = ½a(ū_i + ū_{i+1}) − ½aν(ū_{i+1} − ū_i). Second order in both Δx and Δt. Its diffusion ½ν² D₂ is exactly the amount that cancels the first-order error; it is *less* than the ½ν D₂ of upwind, so the scheme is not monotone and a sharp front acquires wiggles. At ν = 1 it too reduces to the exact shift.

**Minmod-limited upwind (a TVD scheme).** Reconstruct u inside each cell as a straight line ū_i + σ_i(x − x_i)/Δx with a **limited slope** σ_i = minmod(ū_i − ū_{i−1}, ū_{i+1} − ū_i): the smaller of the two one-sided differences when they agree in sign, zero when they do not (at an extremum). The upwind flux of this reconstruction, averaged over the step, is F = a[ū_i + ½(1 − ν)σ_i]. Where the solution is smooth and monotone the slope is close to the true one and the scheme is second order; at a jump or an extremum the limiter cuts the slope and the scheme falls back to first-order upwind, which keeps it bounded and TVD. The price is that even a smooth peak is clipped slightly, so the observed order on a Gaussian is between 1 and 2 (1.88 at the finest grid in `main()`), and the scheme is nonlinear: no single amplification factor describes it.

### 2.3 Stability: von Neumann analysis and the CFL condition

Substitute ū_j^m = gᵐ e^{ijθ} (`P03.11` 2.3). With c = cos θ, s = sin θ:

- upwind: g = 1 − ν(1 − e^{−iθ}), **|g|² = 1 − 2ν(1 − ν)(1 − c)**: ≤ 1 iff 0 ≤ ν ≤ 1.
- Lax–Friedrichs: g = c − iνs, **|g|² = c² + ν²s²**: ≤ 1 iff ν² ≤ 1.
- Lax–Wendroff: g = 1 − iνs − ν²(1 − c), **|g|² = 1 − ν²(1 − ν²)(1 − c)²**: ≤ 1 iff ν² ≤ 1.

All three require **|ν| ≤ 1**: in one step the solution may not travel further than one cell, because the stencil only reaches one cell. This is the **CFL condition** (`@ocw-16-90` 2.6): the numerical domain of dependence must contain the physical one. It is necessary for every explicit scheme on this stencil and here it is also sufficient. The exact factor is e^{−iνθ}, of unit modulus: the true solution loses no amplitude and shifts every mode's phase by νθ. The numerical factors lose amplitude (|g| < 1: **dissipation**, worst for Lax–Friedrichs) and shift the phase by the wrong amount (**dispersion**): `phase_speed_ratio` returns −arg g/(νθ), which for Lax–Wendroff is below 1 for every resolved mode (0.93 at θ = π/4, ν = ½), so its short waves lag the front and the wiggles trail behind it. For a periodic grid the update is a **circulant matrix** whose eigenvalues are exactly g(2πk/n), k = 0 … n − 1; the test builds the matrix column by column from `fv_step`, checks it against `scipy.linalg.circulant`, and matches its eigenvalues to the formulas to 10⁻⁸.

### 2.4 The modified equation of the upwind scheme

Stability says errors do not grow; it does not say what the scheme *does* to the solution. Expand the upwind update in Taylor series about (x_i, t^m), exactly as for FTCS in `P03.11` 2.2: ū + Δt u_t + ½Δt² u_tt = ū − ν(Δx u_x − ½Δx² u_xx) + O(Δx³). Divide by Δt, use u_tt = a²u_xx (valid to the order kept), and collect:

**u_t + a u_x = D u_xx + O(Δx²),  D = ½aΔx(1 − ν).**   (2)

The upwind scheme solves an **advection–diffusion** equation with a numerical diffusivity D that vanishes only at ν = 1. This is the single most useful fact about first-order schemes, because it is quantitative: a Gaussian slug of variance σ² becomes one of variance **σ² + 2Dt** (the diffusion equation's spreading law, `P02.11`). For the feed line on 120 cells (Δx = 2 cm) at ν = ½, D = 0.06 m²/s, and in the 0.1 s to the sensor the variance should grow by 0.0120 m²: the run gives 0.0120 m² (σ from 0.100 m to 0.148 m). The same calculation for Lax–Wendroff gives no u_xx term; its leading error is −(aΔx²/6)(1 − ν²)u_xxx, a **dispersive** term, which is why it does not smear but wiggles. And it explains the sensor readings: numerical diffusion carries some of the slug ahead of the true front (early 10 % arrival) and lowers the peak; the 50 % crossing is almost unaffected because diffusion is symmetric about the front.

### 2.5 What convergence order to expect

For smooth data the orders are the formal ones: 1 for upwind and Lax–Friedrichs, 2 for Lax–Wendroff, between 1 and 2 for the limited scheme. The first-order schemes approach order 1 **from below** (0.84, 0.91, 0.95 on 240 → 1920 cells for upwind) because on coarse grids the numerical diffusion has already flattened the slug and the error saturates; one must refine until 2Dt ≪ σ² before the asymptotic rate shows. For a **discontinuity** the picture changes: a first-order scheme smears a jump over a width proportional to √(Dt) ∝ √Δx, so the L1 error is O(Δx^{1/2}), **order ½**, and Lax–Wendroff's wiggles, which do not shrink in height as the grid is refined, hold it to about 0.6. The TVD scheme reaches about ⅔. The run gives 0.50, 0.60 and 0.65. No scheme on a fixed stencil achieves first order in L1 at a jump; that is a theorem about the problem, not a defect of the code, and it is why a *front* in CFD is always under-resolved and why the modeller must say so.

### 2.6 Where the model stops being valid

- **Plug flow.** A real pipe flow has a velocity profile and turbulent mixing. The eddy diffusivity of a turbulent liquid line is of order 10⁻³ m²/s; the upwind scheme's 0.06 m²/s on 120 cells is sixty times larger, so the simulated smearing is all numerical. §5 computes the grid at which that reverses.
- **Passive scalar.** A mixture-ratio excursion that reaches the injector changes the combustion, which changes the pressure, which changes a. The one-way coupling holds only upstream of the injector.
- **Periodic domain.** A convenience for testing (every quantity is conserved exactly and the exact solution is a shift). A real line has an inlet and an outlet; extension 1.
- **Linear advection.** The conservation form (1) is exact for any flux f(u); the four flux formulas and the amplification factors are for f = au only. Extension 2 takes the limited scheme to Burgers' equation, where the same code with a different flux captures a shock.

## 3. Build the computational model

Open `starter/main.py`. `feed_line`, `semi_discrete_upwind_oracle` and `main` are given.

### Step 1: grid and data

`cell_centres(L, n)` returns ((i + ½)Δx, Δx). `cell_averages(f, L, n, n_gauss)` evaluates f at the Gauss–Legendre points of each cell (`np.polynomial.legendre.leggauss`, mapped to the cell) and takes the weighted sum divided by 2. `exact_advection(f, x, a, t, L)` is f((x − at) mod L).

### Step 2: fluxes and the update

```python
def numerical_flux(u, a, nu, scheme):          # a > 0; mirror the array for a < 0
    u_r = np.roll(u, -1)                        # u_{i+1}
    if scheme == "upwind":
        return a * u
    if scheme == "lax_friedrichs":
        return 0.5 * a * (u + u_r) - 0.5 * (a / nu) * (u_r - u)
    if scheme == "lax_wendroff":
        return 0.5 * a * (u + u_r) - 0.5 * a * nu * (u_r - u)
    slope = minmod(u - np.roll(u, 1), u_r - u)  # minmod: limited slope
    return a * (u + 0.5 * (1.0 - nu) * slope)

def fv_step(u, a, dt, dx, scheme):
    F = numerical_flux(u, a, a * dt / dx, scheme)
    return u - (dt / dx) * (F - np.roll(F, 1))  # (1): F_{i+1/2} - F_{i-1/2}
```

For a < 0, reverse the array, compute with −a, reverse the fluxes back and change their sign (the test checks that reversing the flow mirrors the solution). `fv_step` raises `ValueError` for |ν| > 1. `advect` loops; `sensor_history` records one cell after each step; `arrival_time` interpolates the first threshold crossing.

### Step 3: analysis

`amplification_factor` types in the three complex formulas of §2.3 and raises for `"minmod"`. `phase_speed_ratio` is −arg g/(νθ). `numerical_diffusivity` is (2). `pulse_variance` computes the variance of a periodic distribution about its circular mean (so the slug can sit anywhere on the ring). `total_variation`, `l1_error` and `observed_order` are one-liners. `convergence_study` averages the data, marches at fixed ν, averages the exact solution at the time reached, and returns the L1 errors.

Running `main()` prints

```
feed line: a = 12.0 m/s, L = 2.4 m, n = 120 cells (dx = 2.0 cm), nu = 0.5, dt = 0.833 ms, 120 steps to t = 0.1 s (slug travels 1.20 m)

Gaussian slug after 120 steps: L1 error, peak (exact 0.993), min, total variation (initial 1.987), mass change
  upwind          L1 = 0.0945, peak = 0.672, min = +0.0000, TV = 1.343, sum(u) dx change = +2.4e-19
  lax_friedrichs  L1 = 0.1769, peak = 0.466, min = +0.0000, TV = 0.931, sum(u) dx change = +4.5e-17
  lax_wendroff    L1 = 0.0219, peak = 0.978, min = -0.0076, TV = 1.971, sum(u) dx change = -3.5e-17
  minmod          L1 = 0.0198, peak = 0.875, min = +0.0000, TV = 1.751, sum(u) dx change = -4.1e-17

square slug after 120 steps: L1 error, peak (exact 1.000), min, total variation (initial 2.000), mass change
  upwind          L1 = 0.1556, peak = 0.826, min = +0.0000, TV = 1.652, sum(u) dx change = -3.1e-17
  lax_friedrichs  L1 = 0.2673, peak = 0.570, min = +0.0000, TV = 1.140, sum(u) dx change = +5.2e-17
  lax_wendroff    L1 = 0.1111, peak = 1.200, min = -0.2044, TV = 3.104, sum(u) dx change = +4.7e-17
  minmod          L1 = 0.0643, peak = 0.981, min = +0.0000, TV = 1.962, sum(u) dx change = -5.1e-17
```

followed by the modified-equation check, the sensor table, the stability table, the convergence tables and the figure.

### Verify

- **Analytic cases:** cell averages equal `quad` to 10⁻⁶ and are exact for cubics; the exact solution returns to itself after one lap; every flux equals its textbook form; a uniform state has flux au for every scheme; at ν = 1 every scheme shifts the data by one cell to 10⁻¹².
- **Conservation and symmetry:** Σū_iΔx constant to 10⁻¹³ for all schemes; a → −a mirrors the solution to 10⁻¹².
- **Stability:** |g| ≤ 1 for |ν| ≤ 1 and > 1 at ν = 1.2 for the three linear schemes; `fv_step` refuses |ν| > 1; the update matrix is circulant and its eigenvalues are g(2πk/n) to 10⁻⁸; the closed forms of |g|² hold to 10⁻¹².
- **Monotonicity:** upwind and minmod stay within [0, 1] on the square slug and do not increase TV; Lax–Wendroff goes below −0.1 and above 1.1 and raises TV by more than 0.5.
- **Modified equation:** variance growth equals 2Dt to 0.2 % at ν = 0.25, 0.5 and 0.8; the gap between fully discrete upwind and `solve_ivp`'s semi-discrete solution falls at order 1.0 ± 0.1 in ν.
- **Convergence (observed order):** Gaussian: upwind 0.8–1.05 and rising, Lax–Wendroff 2.00 ± 0.05, minmod 1.5–2; square: upwind 0.50 ± 0.05, minmod above 0.6.
- **The sensor:** the 50 % crossing within 4 ms of the exact front for all schemes; the peak and undershoot of each scheme in the ranges of §2.4.

### Compare

There is no single SciPy routine for a finite-volume advection solve; what SciPy offers is the parts, and the lesson uses them as oracles: `quad` for the averages, `circulant` and `eigvals` for the stability analysis, `solve_ivp` for the semi-discrete (method-of-lines) system that the fully discrete upwind scheme approaches as ν → 0 (the test measures the O(ν) gap). In production CFD codes the time integration *is* done this way (a spatial flux evaluation handed to an explicit Runge–Kutta integrator, `P03.05`), and the limited flux of §2.2 is a close relative of what those codes do at every cell face of a nozzle or wing mesh; `scipy.integrate.solve_ivp` with the flux function as right-hand side is a perfectly good way to run the semi-discrete scheme with error control. What writing the fully discrete steps bought is the Courant number as an explicit knob: the exactness at ν = 1, the instability at ν > 1 and the diffusivity ½aΔx(1 − ν) are statements about the fully discrete scheme that an adaptive integrator hides.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before the run.

1. **Courant-number sweep.** On 240 cells, run the upwind scheme on the Gaussian slug to t = 0.1 s at ν = 0.25, 0.5, 0.9, 1.0 and (with the CFL check bypassed in a copy of `fv_step`, seeding 10⁻⁶ of the sawtooth (−1)ⁱ) 1.05. Predict from (2) that the L1 error scales with 1 − ν, that ν = 1 is exact, and that at ν = 1.05 the sawtooth grows by |g(π)| = 1.1 per step: 10⁻⁶ × 1.1³⁰⁰ ≈ 2.6 × 10⁶ after 300 steps. Record `nu, L1 error, max|u|`.
2. **Convergence on two slugs.** Run `convergence_study` for the four schemes on the Gaussian and for upwind, Lax–Wendroff and minmod on the square slug, n = 240 … 1920, ν = ½. Predict orders 1, 1, 2, 1.5–2 and ½, 0.6, ⅔.
3. **The modified equation.** On 240 cells, measure the variance growth of the Gaussian under upwind at ν = 0.25, 0.5, 0.8 and compare with 2Dt. Then do the same for Lax–Wendroff and predict that the growth is near zero (it has no u_xx term) while the TV rises.
4. **The sensor on four grids.** For the square slug, record the 10 % and 50 % arrival times, peak and minimum at the sensor for n = 60, 120, 240, 480 and all four schemes. Predict: the 50 % time converges to the exact front time (87.5 ms) for every scheme, the upwind peak rises toward 1 as 1 − O(√Δx), the Lax–Wendroff overshoot does **not** shrink with refinement, and the minmod peak is within 2 % of 1 on every grid.

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `fv_advection.png`.

What you should find: upwind L1 errors at n = 1920 in the ratio 7.2 : 4.9 : 1 for ν = 0.25, 0.5, 0.9 against the prediction 7.5 : 5 : 1 from 1 − ν (the ratio approaches the prediction as the grid is refined); the orders 0.84 → 0.95 (upwind), 1.99 → 2.00 (Lax–Wendroff), 1.64 → 1.88 (minmod) on the Gaussian and 0.50, 0.60, 0.65 on the square slug; a variance growth of 0.0120 m² against 0.0120 predicted at ν = ½; and a Lax–Wendroff sensor peak of 1.20–1.23 with a −0.20 to −0.24 undershoot on every grid, while the upwind peak reads 0.83 on 120 cells.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| mass Σū_iΔx, any scheme | conserved exactly | change ≤ 5 × 10⁻¹⁷ | rounding | conservation form (1) |
| ν = 1, any scheme | exact shift | L1 error 4 × 10⁻¹⁶ | rounding | g = e^{−iθ} |
| max|g|, ν = 1.1 | 1.2 (upwind), 1.1 (LF), 1.42 (LW) | 1.2000, 1.1000, 1.4200 | 0 | closed forms |
| upwind D, 120 cells, ν = ½ | 0.0600 m²/s; 2Dt = 0.0120 m² | 0.01200 m² measured | < 10⁻⁵ | modified equation (2) |
| Gaussian order: upwind / LW | 1 / 2 | 0.95 / 2.00 at n = 1920 | pre-asymptotic | 2Dt not yet ≪ σ² |
| square order: upwind | ½ | 0.50 | 0 | front width ∝ √Δx |
| LW phase speed, θ = π/4, ν = ½ | < 1 (lagging) | 0.93 | — | dispersive error term |
| sensor 50 % crossing | 87.5 ms (exact front) | 88.4 / 90.4 / 89.9 / 88.3 ms | ≤ 3 ms | sensor cell centre at 1.81 m is 0.8 ms downstream of 1.80 m; diffusion is symmetric |
| sensor peak | 1.000 | 0.83 / 0.57 / 1.20 / 0.98 | up to −43 % | numerical diffusion; LW overshoot |

The engineering question has two answers. For the designer: the slug front reaches the sensor at 87.5 ms and the peak is 1.0 by definition of the slug, because the exact solution is a translation; any physical smearing comes from turbulent mixing (D ≈ 10⁻³ m²/s, which would widen σ from 0.100 to 0.101 m in 0.1 s, a 1 % effect). For the modeller: on a 120-cell grid the upwind scheme reports the peak 17 % low and the 10 % arrival 10 ms early, errors sixty times larger than the physics it omits, because its numerical diffusivity 0.06 m²/s dwarfs the eddy diffusivity. To make the numerical diffusion smaller than the physical one at ν = ½ requires ½aΔx(1 − ν) < 10⁻³, Δx < 0.33 mm, **7200 cells** and 7200 steps for 0.1 s; the limited scheme gets the peak to 2 % on 120 cells. That is the whole case for second-order TVD schemes in one table, and it is why first-order upwind survives only as the fallback inside limiters.

The Lax–Wendroff row is the other lesson: a second-order scheme whose observed order on the smooth slug is a perfect 2.00 reports a mixture fraction of −0.20 at the sensor, which no propellant can have, and the undershoot does not go away with refinement. A controller fed that signal would trim the wrong way. Order of accuracy is a statement about smooth solutions; boundedness is a separate property and for this application a mandatory one.

Improvements worth making: an inlet/outlet boundary (extension 1); a physical diffusion term so that the model has a smearing of its own to compare with the numerical one (extension 3); a smoother limiter (van Leer or superbee) to recover the clipped peak of the minmod scheme. Improvement worth *not* making: running Lax–Wendroff at ν = 1 "because it is exact there". It is, for constant a on a uniform grid, and for nothing else; a real line has a varying speed and the exactness evaporates while the wiggles remain.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), the four fluxes, the three |g|² formulas and (2); `fv_advection.png`; the four tables of §4; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a scheme and grid for the feed-line model, with the peak and arrival errors it leaves and the physical diffusivity it should be compared with.

### Adversarial self-checks

- **What if the slug is half as long (0.15 m)?** The exact answer is unchanged in arrival time and peak; the numerical peak error roughly doubles for the first-order schemes (the front width ∝ √(Dt) is unchanged while the slug is shorter), and the Lax–Wendroff overshoot is unchanged. Check with `square_pulse(x, 0.6, 0.15)`.
- **What if a → 0?** ν → 0 at fixed Δt and the schemes become the identity; D → 0 too but so does the distance travelled, and the slug sits still. Nothing is learned and nothing breaks. At a → −a everything mirrors (tested).
- **Which assumption fails first in the real world?** Plug flow. A turbulent pipe profile is flat to within 10–20 %, but the near-wall fluid lags, so a real slug develops a tail that no scalar advection model has. Second: the passive-scalar assumption at the injector.
- **What measurement could falsify the model?** The sensor's own trace: a 10 %–90 % rise time at the sensor of a few milliseconds is consistent with turbulent mixing; a rise time of 15 ms (the upwind result on 120 cells) is not, and if the measured trace were that slow the slug itself would be longer than assumed.
- **Does the result survive a change of grid, Courant number or scheme?** The 50 % arrival does (within 3 ms on every grid and scheme); the peak and the 10 % arrival do not for the first-order schemes until n ≳ 500, and the Lax–Wendroff undershoot never does. Report the quantity that survives and flag the ones that do not.
- **Which decision changes if the uncertainty doubles?** If the controller's threshold were 0.5 instead of 0.8, every scheme's peak would clear it and the scheme choice would not matter for the trip decision, only for its timing (±3 ms). At a threshold of 0.8 the upwind model on 120 cells says the controller never trips; the exact solution says it does. That is a decision reversed by numerical diffusion.

## Check yourself

1. Write the finite-volume update (1) and explain, without reference to any flux formula, why it conserves Σū_iΔx.
2. Define the Courant number and state the CFL condition for the schemes of this lesson. Why can an explicit three-cell stencil never be stable at ν = 2?
3. Give the upwind flux for a > 0 and for a < 0, and write the update as a convex combination. What property follows for 0 ≤ ν ≤ 1?
4. State the modified equation of the upwind scheme. By how much does the variance of a Gaussian grow in time t, and for which ν is the growth zero?
5. Why does Lax–Wendroff overshoot at a jump although it is second order and stable, and does refining the grid cure it?
6. What is the role of the minmod limiter, where does it reduce the scheme to first order, and why is that acceptable?
7. Why is the observed L1 order of the upwind scheme ½ at a discontinuity and 1 on smooth data?
8. On the 120-cell grid the upwind scheme reads the slug peak as 0.83. Is that a property of the propellant or of the grid, and what one number from §2.4 lets you say so?
