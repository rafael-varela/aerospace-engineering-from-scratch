# P02.11  PDEs: Heat, Wave & Laplace Equations

> Separate the variables and a partial differential equation becomes a family of ordinary ones; the boundary conditions pick the family members, and the initial condition pays each one its Fourier coefficient.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the one-dimensional heat equation from an energy balance, the wave equation from Newton's law on a string element, and Laplace's equation as the steady limit, naming each assumption. 2. Solve each by separation of variables for given boundary conditions, including the mixed (fixed face, insulated face) case of a re-entry tile, and state the Fourier number that governs the heat solution. 3. Implement and test, from scratch, an explicit finite-difference marcher for the heat equation (FTCS), a leapfrog marcher for the wave equation, d'Alembert's travelling-wave solution with the odd periodic extension, and Jacobi relaxation for Laplace's equation, and prove each agrees with the series. 4. Decide, from the series, how thick a tile must be to keep its back face below a limit, and say which assumption of the model fails first. |
| **Prerequisites** | `P02.09` (Fourier sine series and orthogonality), `P02.02` (second-order linear ODEs: the separated equations), `P02.01` (first-order ODEs: the time factor), `P01.12` (the divergence theorem, used in the energy balance), `P01.04` (integration by parts for the coefficients), `P00.03`, `P00.07`. Skill check: you can integrate x sin(nπx/L) by parts and solve X″ + λX = 0 with X(0) = 0. |
| **Engineering question** | A 50 mm silica tile whose outer face is held at 1300 K for 1200 s: does its insulated back face stay below the 450 K limit of the aluminium beneath it, and if not, how thick must the tile be? Alongside it: at what frequencies does a 1.2 m tensioned sensor cable ring when plucked, and what is the steady temperature field in a plate heated along one edge? |
| **Assumptions and limits** | Heat: constant conductivity, density and specific heat (real silica tiles have k rising three-fold between 300 and 1300 K); one-dimensional conduction (edges far away); hot face at a fixed temperature (in reality a radiative balance); back face perfectly insulated (conservative: the structure is a heat sink). Wave: small slopes, uniform tension and mass, no bending stiffness, no damping (a real cable's modes decay and are slightly inharmonic). Laplace: steady state, constant conductivity, edge temperatures prescribed. The explicit heat marcher is stable only for αΔt/Δx² ≤ ½; leapfrog only for cΔt/Δx ≤ 1; Jacobi converges slowly (sweeps grow as the square of the grid count). |
| **Analytical work** | Energy-balance derivation of T_t = αT_xx; string-element derivation of u_tt = c²u_xx; separation of variables with Dirichlet–Dirichlet, Dirichlet–Neumann and rectangular boundary conditions; the eigenvalues λ_n = (2n+1)π/2L of the mixed problem and the coefficients 4/((2n+1)π); the Fourier number and the one-term approximation; d'Alembert's solution and the odd periodic extension; the plucked-string coefficients; the sinh series for the rectangle; von Neumann stability of FTCS; the CFL condition. |
| **Computational deliverable** | `starter/main.py`: `sine_coefficients`, `heat_series`, `tile_temperature`, `tile_thickness_for_limit`, `heat_ftcs`, `string_mode_frequency`, `plucked_string_coefficients`, `wave_series`, `odd_periodic_extension`, `dalembert`, `wave_leapfrog`, `string_energy`, `laplace_rectangle`, `laplacian_residual`, `laplace_jacobi`. Tests prove: the heat series satisfies the PDE (residual below 10⁻⁹), its boundary and initial conditions, the maximum principle and the Fourier-number scaling; FTCS matches the series and converges at second order in Δx; FTCS refuses r > ½; the wave series equals d'Alembert to 0.4 % of the pluck height at every tested time; leapfrog tracks d'Alembert, conserves energy to 0.5 % and is exact at Courant number 1; the Laplace series is harmonic, obeys the mean-value property, and Jacobi converges to it at second order with a sweep count that grows as the grid squared. |
| **Practical activity** | **computational experiment**: a grid-refinement study for FTCS and Jacobi, a Courant-number sweep for leapfrog, and a tile-thickness design sweep. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `pde_solutions.png`, the three comparison tables (series against FTCS, d'Alembert and leapfrog, Jacobi), the tile-thickness result, and a short note on which tile assumption you would remove first. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation of the mixed-boundary tile series is complete with the coefficients; all public tests pass; the hand one-term estimate of the back-face temperature and the required thickness agree with the code to 2 %; the note names an assumption and a measurable sign of its failure. |
| **References** | `@mit-ocw-18-303-linear-pdes` Lecture 1 "Overview of linear PDEs and analogies with matrix algebra", Lecture 2 "Poisson's equation and eigenfunctions in 1d: Fourier sine series", Lecture 8 "Music and wave equations, Separation of variables, in time and space", Lecture 10 "General Dirichlet and Neumann boundary conditions", Lecture 22 "Von Neumann analysis and the heat equation", Lecture 25 "Traveling waves: D'Alembert's solution" (required-free; the lecture summaries page with handout PDFs); Lecture 23 "Algebraic properties of wave equations and unitary time evolution, Conservation of energy in a stretched string" (supplementary-free); `@mit-ocw-18-03-differential-equations` notes 20 "Fourier series", 21 "Operations on Fourier series" (prerequisite); `@dawkins-differential-equations-notes` Chapter 9, sections "The Heat Equation", "Separation of Variables", "Solving the Heat Equation", "Laplace's Equation", "Vibrating String" (supplementary-free: a gentler walk-through of the same examples). |
| **Extensions** | ★ Replace the fixed hot-face temperature by a radiative balance εσT⁴ = q_in and march it with FTCS; compare the back-face history. ★ Add bending stiffness to the cable (u_tt = c²u_xx − (EI/μ)u_xxxx) and show the modes become inharmonic. ★ Solve Laplace's equation on the rectangle with Gauss–Seidel and SOR and measure how the sweep count falls. |
| **Time** | unmeasured |

## 1. The engineering problem

A re-entry vehicle's windward surface is covered with rigid silica-fibre tiles. Their job is not to *survive* 1300 K (silica does that easily) but to *delay* the heat: behind each tile is an aluminium skin that loses strength above about 450 K, and the heating pulse lasts about 20 minutes. The tile is a one-dimensional heat-conduction problem: a slab of thickness L, outer face held at the surface temperature T_s, back face (bonded to the structure) treated as insulated, starting uniform at T₀. The question is the back-face temperature after 1200 s and, if it exceeds 450 K, the thickness that would hold it there. The answer must come with a time scale, because the tile works by being slow, and the time scale is set by one number, the thermal diffusivity α = k/(ρc_p).

Two more problems ride along because the same method solves them. A 1.2 m sensor cable under 200 N tension is plucked 5 mm at one point and released: at which frequencies does it ring, and what does its shape look like 2 ms later? And a 0.3 m × 0.2 m equipment plate has one edge heated 60 K above the other three: what is the steady temperature in the middle? The three equations, **heat**, **wave** and **Laplace**, are the three classical linear PDEs; everything in structural dynamics, acoustics, heat transfer and potential flow is one of them with more dimensions and worse boundaries.

Success criteria, written down first: each series solution satisfies its PDE to 10⁻⁹ in a finite-difference residual and its boundary conditions exactly; the tile back-face temperature from a one-term hand calculation agrees with the full series to 0.5 %; the finite-difference marcher for the heat equation agrees with the series to 0.2 K on a 1000 K range and its error falls four-fold when the grid is halved; the cable's series and d'Alembert's formula agree to 0.4 % of the pluck height; the leapfrog marcher conserves the cable's energy to 0.5 %; and Jacobi relaxation agrees with the plate series to 0.03 K.

## 2. Understand and derive

### 2.1 The heat equation from an energy balance

Take a slab element between x and x + dx of unit cross-section. Fourier's law (**assumption: conduction only, isotropic, k constant**) says the heat flux is q = −k ∂T/∂x [W/m²]. The energy stored in the element is ρ c_p T dx, and its rate of change equals the flux in minus the flux out:

ρ c_p (∂T/∂t) dx = q(x) − q(x + dx) = −(∂q/∂x) dx = k (∂²T/∂x²) dx,

**∂T/∂t = α ∂²T/∂x², α = k/(ρ c_p) [m²/s].**   (1)

(In three dimensions the same balance over a volume, with the divergence theorem of `P01.12`, gives T_t = α∇²T.) Units: K/s = (m²/s)(K/m²). ✓ For the tile, k = 0.05 W/(m K), ρ = 144 kg/m³, c_p = 700 J/(kg K): α = 4.96 × 10⁻⁷ m²/s. The only length in the problem is L and the only time scale that can be built from α and L is **L²/α = 5040 s** for L = 50 mm. That number already says the tile will be well on its way to equilibrium after 1200 s: the ratio **Fo = αt/L² = 0.238** is the **Fourier number**, and every heat-conduction answer depends on it alone.

### 2.2 Separation of variables: the Dirichlet slab

Look for solutions of (1) of the form T(x, t) = X(x) Θ(t) (`@mit-ocw-18-303-linear-pdes` Lecture 2). Substituting and dividing by XΘ: Θ′/(αΘ) = X″/X. The left side depends only on t, the right only on x, so both equal a constant, −λ². Two ODEs (`P02.01`, `P02.02`): Θ′ = −αλ²Θ, so Θ = e^{−αλ²t}, and X″ + λ²X = 0, so X = A sin λx + B cos λx.

With both faces held at zero (T(0) = T(L) = 0), B = 0 and sin λL = 0: **λ_n = nπ/L**, n = 1, 2, …. The boundary conditions have chosen a discrete set of **eigenfunctions** sin(nπx/L), each with its own decay rate α(nπ/L)². Because (1) is linear, any sum of them is a solution:

**T(x, t) = Σ_n b_n sin(nπx/L) e^{−α(nπ/L)² t},**   (2)

and at t = 0 this must equal the initial temperature f(x). That is a Fourier sine series (`P02.09`), so by orthogonality **b_n = (2/L)∫₀ᴸ f(x) sin(nπx/L) dx**. For f = 1, b_n = 4/(nπ) for odd n and 0 for even n. The n-th mode decays n² times faster than the first: after a Fourier number of 0.2, the second surviving mode (n = 3) is down by e^{−9π²·0.2} ≈ 2 × 10⁻⁸ relative to the first. This is why heat problems become one-term problems so quickly, and why a tile's history is dominated by a single exponential.

### 2.3 The tile: a fixed face and an insulated face

Set θ = T − T_s. Then θ obeys (1) with θ(0, t) = 0 at the hot face, ∂θ/∂x(L, t) = 0 at the insulated back (**Neumann** condition: zero flux), and θ(x, 0) = T₀ − T_s. The separated X must satisfy X(0) = 0, so X = sin λx, and X′(L) = λ cos λL = 0, so

**λ_n = (2n + 1)π/(2L), n = 0, 1, 2, …**   (3)

(`@mit-ocw-18-303-linear-pdes` Lecture 10). The eigenfunctions are quarter-waves, not half-waves: the lowest has a node at the hot face and an antinode at the back. The coefficients of the constant θ₀ = T₀ − T_s are

b_n = (2/L)∫₀ᴸ θ₀ sin λ_n x dx = (2θ₀/(Lλ_n))(1 − cos λ_n L) = 2θ₀/(Lλ_n) = **4θ₀/((2n + 1)π)**,

since cos λ_n L = cos((2n+1)π/2) = 0. So

**T(x, t) = T_s + (T₀ − T_s) Σ_{n≥0} [4/((2n+1)π)] sin λ_n x · e^{−αλ_n² t}.**   (4)

At the back face sin λ_n L = (−1)ⁿ. With Fo = 0.238, the first exponent is αλ₀²t = (π²/4)Fo = 0.587, and the one-term back-face temperature is T_s + (T₀ − T_s)(4/π)e^{−0.587} = 1300 − 1000 × 1.2732 × 0.5557 = **592.4 K**; the full series gives **594.6 K** (the n = 1 term, with its minus sign and e^{−9×0.587}, adds 2.2 K). The back face is 145 K over the limit. The 50 mm tile is too thin for this pulse.

**How thick?** Invert the one-term formula: (450 − 1300)/(300 − 1300) = 0.85 = (4/π)e^{−(π²/4)Fo} gives Fo = 0.1638, so L = √(αt/Fo) = √(4.96 × 10⁻⁷ × 1200/0.1638) = **60.3 mm**. The code's root-finder on the full series gives **61.4 mm**; the 1.1 mm difference is the second term, which at the smaller Fourier number is less negligible (it is −0.7 % of θ₀ there). The required thickness scales as √t: four times the pulse, twice the tile, which the test checks exactly.

Dimensional reading: nothing in (4) depends on k, ρ or c_p separately, only on α, and nothing on L or t separately, only on Fo. Two tiles with the same Fo have the same temperature profile in x/L. That is the second test.

### 2.4 The wave equation from a string element

A string (or cable) of mass per length μ under tension T, displaced by u(x, t) with small slope (**assumption: |u_x| ≪ 1, so the tension is uniform and the length unchanged**). The transverse force on the element dx is the difference of the vertical components of tension at its ends, T u_x(x + dx) − T u_x(x) = T u_xx dx, and Newton's law gives μ dx u_tt = T u_xx dx:

**u_tt = c² u_xx, c = √(T/μ) [m/s].**   (5)

Units: (m/s)² (1/m) = m/s². ✓ For the cable, c = √(200/0.010) = **141.42 m/s**. Separation of variables with u = 0 at both ends gives X_n = sin(nπx/L) again, but now the time equation is Θ″ + (nπc/L)²Θ = 0, an undamped oscillator (`P02.02`) with angular frequency ω_n = nπc/L:

**u(x, t) = Σ_n [b_n cos ω_n t + (d_n/ω_n) sin ω_n t] sin(nπx/L), f_n = n c/(2L).**   (6)

b_n are the sine coefficients of the initial shape, d_n those of the initial velocity. The modes are **harmonic** (integer multiples of f₁ = c/2L = **58.93 Hz**; f₂ = 117.85 Hz, f₃ = 176.78 Hz), which is why a plucked string has a pitch. A cable with bending stiffness is slightly inharmonic (extension), and a real cable's modes decay; (5) has no damping.

**The plucked cable.** Initial shape a triangle of height h = 5 mm at a = 0.3 m, released from rest (d_n = 0). Integrating by parts (`P01.04`),

**b_n = 2hL²/(n²π² a(L − a)) sin(nπa/L):**   (7)

3.821, 1.351, 0.425, 0.000 mm for n = 1 … 4. The fourth coefficient vanishes because sin(4π × 0.3/1.2) = sin π = 0: plucking at a quarter of the length cannot excite the mode with a node there. The coefficients fall as 1/n², so 60 terms represent the kink to 0.045 mm (1 % of h); `sine_coefficients` by quadrature reproduces (7) to 7 × 10⁻¹¹ m.

### 2.5 d'Alembert: the same solution as travelling waves

Any function of x − ct or x + ct satisfies (5) (check: u_tt = c²F″, u_xx = F″). For a string released from rest in shape F on the whole line, u = ½[F(x − ct) + F(x + ct)]: the shape splits into two half-height copies travelling in opposite directions (`@mit-ocw-18-303-linear-pdes` Lecture 25). On a finite string with fixed ends, the trick is to extend f as an **odd, 2L-periodic** function F: then F(x − ct) + F(x + ct) vanishes at x = 0 (odd) and at x = L (odd about L, by periodicity), so the ends stay fixed by themselves, and each wave that reaches an end comes back inverted. The series (6) with d_n = 0 is the Fourier expansion of exactly this expression (cos ω_n t sin k_n x = ½[sin k_n(x + ct) + sin k_n(x − ct)]), so the two must agree. The run gives a maximum difference of 0.045 mm at t = 0 (the truncated series at the kink), 0.020 mm at 1 ms, 0.018 mm at 2 ms: the series' Gibbs error at the kinks and nothing else. After one round trip 2L/c = **16.97 ms** the string is back in its initial shape, and after half of that it is inverted and mirrored, which the tests check to 10⁻¹² because d'Alembert's formula is exact.

### 2.6 Laplace's equation: the steady state

Let t → ∞ in the two-dimensional heat equation: T_t → 0 and **T_xx + T_yy = 0**, Laplace's equation (`@mit-ocw-18-303-linear-pdes` Lecture 1). Its solutions, **harmonic functions**, have no local maxima or minima inside the domain (the **maximum principle**: heat cannot pile up in steady state) and the value at a point is the average over any circle around it (the **mean-value property**). Both are tests.

On the rectangle 0 < x < a, 0 < y < b with T = 0 on three edges and T = f(x) on y = b, separate T = X(x)Y(y): X″/X = −Y″/Y = −k². The x-conditions give X_n = sin(nπx/a); the y-equation Y″ = k²Y gives sinh and cosh, and Y(0) = 0 selects sinh:

**T(x, y) = Σ_n b_n sin(nπx/a) sinh(nπy/a)/sinh(nπb/a),**   (8)

with b_n the sine coefficients of f. Note how the separation constant has changed sign relative to the heat and wave cases: one direction oscillates, the other grows, which is why Laplace problems need boundary data all the way round and heat problems need it only at the start. For the plate (a = 0.3 m, b = 0.2 m) with a parabolic edge profile f = 60 × 4x(a − x)/a² K, b_n = 32 × 60/(nπ)³ for odd n: 61.92, 2.293, 0.495 K for n = 1, 3, 5. The centre temperature is **19.25 K** above the cold edges (one term alone: 19.35 K).

### 2.7 Finite differences: the independent check

Each series is checked by a method that knows nothing about eigenfunctions.

**FTCS for heat.** Replace T_t by a forward difference and T_xx by the centred second difference (`@mit-ocw-18-303-linear-pdes` Lecture 22):

**T_j^{m+1} = T_j^m + r (T_{j+1}^m − 2T_j^m + T_{j−1}^m), r = αΔt/Δx².**   (9)

Substituting a Fourier mode T_j^m = g^m e^{iκjΔx} gives the amplification factor g = 1 − 4r sin²(κΔx/2), and |g| ≤ 1 for every κ only if **r ≤ ½** (von Neumann stability). At r = 0.51 the shortest wave grows by 4 % per step and the solution reaches 10¹⁶ K in 1000 steps; `heat_ftcs` raises `ValueError` rather than produce that. The insulated face uses a ghost point T_{J+1} = T_{J−1}, which keeps the scheme second order there. With 51 points, r = 0.4 and 1489 steps, the FTCS tile profile at 1200 s differs from (4) by at most 0.179 K; at 600 s the error is 0.515, 0.129, 0.032 K for 26, 51, 101 points: ratio 4.0 per halving, second order in Δx, as the centred difference promises.

**Leapfrog for the wave.** Centred differences in both t and x:

**u_j^{m+1} = 2u_j^m − u_j^{m−1} + s² (u_{j+1}^m − 2u_j^m + u_{j−1}^m), s = cΔt/Δx.**   (10)

The **Courant number** s must satisfy s ≤ 1 (the **CFL condition**: the numerical domain of dependence must contain the physical one; a wave must not cross more than one cell per step). At s = 1 exactly the scheme reproduces d'Alembert's solution on the grid *exactly* (the test demands 10⁻¹²), a happy accident of the one-dimensional wave equation. At s = 0.8 with 241 points, the error against d'Alembert is 0.019 mm at 2 ms and 0.030 mm at 10 ms; the energy E = ½∫(μu_t² + Tu_x²)dx is 11.06 mJ at both times against 11.05 mJ initially (the exact triangle energy is ½Th²(1/a + 1/(L − a)) = 11.11 mJ; the 0.5 % deficit is the central difference smoothing the kink; `@mit-ocw-18-303-linear-pdes` Lecture 23 derives the conservation law).

**Jacobi for Laplace.** Replace each interior value by the average of its four neighbours (weighted by 1/Δx² and 1/Δy²) and repeat until nothing changes. This is the five-point Laplacian set to zero, solved by the slowest possible iteration; on a 61 × 41 grid it needs 7238 sweeps for 10⁻¹⁰, and the sweep count grows as the square of the grid count (526, 1962, 7238 for the three grids in part 4). The converged grid differs from (8) by 0.054, 0.014, 0.0035 K: second order again. Phase 03 replaces Jacobi by SOR, multigrid and direct sparse solves; the point here is that an iteration with no Fourier content in it lands on the sinh series.

### 2.8 Where the models stop being valid

- **Tile: constant properties.** Silica's conductivity roughly triples from 300 to 1300 K (radiation through the fibres); a constant-α model with the cold value under-predicts heat soak, with the hot value over-predicts it. The remedy is a temperature-dependent k in the FTCS marcher, which the series cannot do; that is why we wrote the marcher.
- **Tile: fixed hot-face temperature.** The surface is really in radiative balance with the boundary layer, and its temperature follows the heating profile. A step is the conservative bound.
- **Tile: insulated back.** The structure conducts heat away; the real back-face temperature is lower. Conservative.
- **Cable: no stiffness, no damping.** Modes are harmonic and ring forever in (5); real cables have inharmonicity of order (n²π²EI/TL²) and decay in tens of cycles.
- **Plate: steady state, prescribed edges.** Real edges have thermal resistance; the result is the steady limit only.

## 3. Build the computational model

Open `starter/main.py`. `reentry_tile`, `sensor_cable` and `main` are given.

### Step 1: sine coefficients and the heat series

```python
def sine_coefficients(f, L, n_terms, n_quad=8001):
    x = np.linspace(0.0, L, n_quad)
    n = np.arange(1, n_terms + 1)
    integrand = f(x)[None, :] * np.sin(np.outer(n, x) * np.pi / L)
    return (2.0 / L) * np.trapezoid(integrand, x, axis=1)

def heat_series(x, t, L, alpha, b):
    n = np.arange(1, b.size + 1); lam = n * np.pi / L
    return np.sin(np.outer(x, lam)) @ (b * np.exp(-alpha * lam**2 * t))   # (2)
```

Raise `ValueError` for t < 0: the heat equation run backwards amplifies every mode by e^{+αλ²|t|} and is ill-posed. `tile_temperature` is (4) with λ_n from (3); `tile_thickness_for_limit` wraps `scipy.optimize.brentq` on T(L, t; L) − T_limit, after checking the bracket.

### Step 2: FTCS

Compute r, refuse r > ½, then loop: interior update from (9); each boundary either reset to its Dirichlet value or updated with the ghost-point formula `u[0] + 2 r (u[1] − u[0])`. The tests check the insulated–insulated slab conserves ∫T dx, which is only true if the ghost point is done right.

### Step 3: the wave equation three ways

`string_mode_frequency` is n c/(2L). `plucked_string_coefficients` is (7). `wave_series` is (6); allow `b_vel=None`. `odd_periodic_extension(f, L)` returns a function: reduce x modulo 2L, use f on [0, L] and −f(2L − x) on (L, 2L). `dalembert` averages F(x − ct) and F(x + ct). `wave_leapfrog` needs two levels to start: take u¹ = u⁰ + Δt v⁰ + ½s² D₂u⁰ (a Taylor step), then iterate (10); return the velocity as a centred difference across the final level so that `string_energy` can be evaluated at one instant.

### Step 4: Laplace

`laplace_rectangle` is (8); write the sinh ratio as e^{k(y−b)}(1 − e^{−2ky})/(1 − e^{−2kb}) so that n = 40 on a 0.3 m plate does not overflow. `laplacian_residual` is the five-point stencil on the interior. `laplace_jacobi` builds the grid, sets the top edge, and sweeps with the weighted average until the largest change is below `tol` times the largest boundary value; count the sweeps and raise `RuntimeError` at `max_iter`.

Running `main()` prints

```
TILE: alpha = 4.960e-07 m^2/s, L = 50 mm, diffusion time L^2/alpha = 5040 s, Fourier number at 1200 s = 0.2381
back-face temperature after 1200 s: series 594.6 K, one-term 592.4 K, limit 450 K
thickness that holds the back face at 450 K after 1200 s: 61.4 mm
FTCS with 51 points, r = 0.400, 1489 steps: max |FD - series| = 0.179 K, back face 594.4 K
CABLE: c = sqrt(T/mu) = 141.42 m/s, f1 = 58.93 Hz, f2 = 117.85 Hz, round trip 2L/c = 16.97 ms
t = 2 ms: max |series(60 terms) - d'Alembert| = 0.0179 mm
leapfrog s = 0.80, 71 steps to 2.01 ms: max |FD - d'Alembert| = 0.0188 mm, energy 11.0604 mJ (initial 11.0494 mJ)
Jacobi on 61 x 41: 6205 sweeps; max |Jacobi - series| = 0.0289 K; centre temperature series 19.2511 K
```

(the Jacobi sweep count in `main()` uses tol = 10⁻⁹; part 4 uses 10⁻¹⁰) and saves `pde_solutions.png`.

### Verify

- **Analytic cases:** the sine coefficients of 1 are 4/(nπ) for odd n; of sin(3πx/L) they are the unit vector e₃; the closed-form plucked-string coefficients equal the quadrature ones to 10⁻⁹; a single sinh term equals its formula to 10⁻¹²; at Courant number 1 leapfrog equals d'Alembert to 10⁻¹².
- **PDE residuals:** heat series to 10⁻⁹ K/s; wave series to 10⁻⁴ relative; Laplace series below 10⁻³ of T/Δx² on a fine grid; the Jacobi grid below 10⁻⁶ of the same.
- **Boundary and initial conditions:** hot face at T_s to 10⁻⁹; insulated-face slope below 10⁻³ of the mean gradient; fixed cable ends to 10⁻¹²; plate edges to 10⁻¹².
- **Limiting cases and conservation:** t → ∞ gives T_s everywhere; the Fourier-number scaling (same x/L and αt/L² give the same T) to 10⁻¹⁰; the insulated–insulated slab conserves heat content to 10⁻⁶; leapfrog conserves energy to 0.5 %; the string returns to its shape after 2L/c and is inverted after L/c; the maximum principle and the mean-value property for the plate.
- **Convergence and stability:** FTCS error ratio > 3 per halving; Jacobi error ratio > 3 per halving and sweep ratio > 2.5; FTCS refuses r > ½; leapfrog refuses s > 1; Jacobi raises at `max_iter`.
- **Design:** the required tile thickness reproduces T_limit to 10⁻⁶ K, exceeds 50 mm, falls for a higher limit, and doubles for four times the pulse.

### Compare

`scipy` has no separation-of-variables solver, because there is nothing to solve: the method is the derivation, and the computation is a sum. What practice uses instead is the finite-difference (or finite-element) marcher with real material data, because it handles temperature-dependent properties, radiative faces and multidimensional geometry that no series can. The series is the *validation* of the marcher: before trusting an FTCS code with k(T), run it with constant k against (4) and demand the 0.18 K agreement and the second-order convergence shown here. The same is true for leapfrog against d'Alembert and for any Laplace solver against (8). An engineer who has both a closed form and a grid code, and knows how far apart they are allowed to be, has a validated tool; one who has only the grid code has a plot.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Four studies; write the prediction before each run.

1. **FTCS refinement.** Tile at 600 s with 26, 51, 101, 201 points, r = 0.25: tabulate the maximum error against (4). Predict the ratio between successive rows (4) and the step count (grows as 4× per halving, since Δt ∝ Δx²). Record `nx, n_steps, max_err_K`.
2. **Stability edge.** Same tile, 51 points, r = 0.25, 0.45, 0.5, 0.51 (remove the guard temporarily for the last one, or write your own loop): record the maximum |T| after 100, 500, 1000 steps. Predict from g = 1 − 4r at κΔx = π which run blows up and how fast (1.04 per step at r = 0.51).
3. **Courant sweep.** Cable with 241 points, s = 0.5, 0.8, 1.0, 1.01: error against d'Alembert at 10 ms and the energy drift. Predict that s = 1 is exact and s = 1.01 is refused.
4. **Tile design.** Required thickness for T_limit = 400, 450, 500, 550 K at 1200 s, and for 450 K at 600, 1200, 2400 s. Predict the √t scaling first.

Raw data format: four CSV tables with a one-line header naming the fixed parameters; `pde_solutions.png`.

What you should find: FTCS errors 0.515, 0.129, 0.032 K and then about 0.008 K, ratio 4; the r = 0.51 run reaches 10⁸ K by step 500 and 10¹⁶ K by step 1000 while r = 0.5 stays bounded (but with a visible 2Δx zigzag, since g = −1 there: stable, not accurate); leapfrog at s = 1 agrees with d'Alembert to rounding and at s = 0.5 is worse than at 0.8 (more steps, more dispersion); the thickness for 450 K goes 43.4, 61.4, 86.9 mm for 600, 1200, 2400 s, the √2 ratios of the Fourier scaling.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| tile back face, 1200 s | 592.4 K (one term) | 594.6 K (series), 594.4 K (FTCS) | 0.4 %, 0.03 % | second term; Δx² error |
| required thickness | 60.3 mm (one term) | 61.4 mm | 1.8 % | second term at Fo = 0.16 |
| FTCS error vs grid | ∝ Δx² | 0.515 → 0.129 → 0.032 K | ratio 4.0 | centred second difference |
| cable f₁, f₂ | 58.93, 117.85 Hz | same | 0 | exact |
| b₁ … b₄ | 3.821, 1.351, 0.425, 0 mm | same, quadrature to 7 × 10⁻¹¹ | 0 | exact |
| series vs d'Alembert, 2 ms | 0 | 0.018 mm | 0.4 % of h | 60-term Gibbs error at the kinks |
| leapfrog energy | 11.11 mJ (exact triangle) | 11.05–11.06 mJ | 0.5 % | central difference across the kink |
| plate centre | 19.35 K (one term) | 19.25 K | 0.5 % | n = 3 term, −0.10 K |
| Jacobi vs series | 0 | 0.054 → 0.014 → 0.0035 K | ratio 4 | five-point stencil, Δx² |

Every difference is named, every one shrinks under refinement at the predicted rate, and all are far below the uncertainty in the inputs: α of a silica tile is not known to better than 20 % (it depends on temperature, density batch and orientation), which moves the 50 mm tile's back-face temperature between 510 and 672 K and the required thickness between 55 and 67 mm (±10 %). The design answer is therefore "60 mm plus a margin set by the α uncertainty", not "61.4 mm".

What the agreement does not show: that a real tile follows (1) with constant α. The strongest effect left out is k(T); the second is the radiative surface. Both are handled by the FTCS marcher with a few more lines (extension 1), which is the improvement worth making. Improvement worth *not* making: adding terms to the series beyond n = 3 for the back-face temperature at Fo > 0.2; they change nothing at the fourth figure.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (3), (4) and (7); `pde_solutions.png`; the four tables of part 4; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the back-face temperature, the required thickness with its uncertainty, and the assumption you would remove first.

### Adversarial self-checks

- **What if the pulse is twice as long?** Fo doubles to 0.476; the back face of the 50 mm tile reaches about 1300 − 1000 × 1.273 × e^{−1.175} = 907 K; the required thickness grows by √2 to 86.9 mm. Check both with the code.
- **What if α goes to zero? To infinity?** α → 0: the tile never warms (Fo → 0, T = T₀ everywhere behind the face). α → ∞: the tile is instantly at T_s. Both are limits of (4).
- **Which assumption fails first in the real world?** Constant k: by 1000 K the tile conducts three times better than the cold value, so a cold-value α under-predicts the soak. Then the fixed hot-face step, which is conservative.
- **What measurement could falsify the model?** Two thermocouples, one at mid-depth and one at the back face. Fit α from the back-face history alone; (4) then predicts the mid-depth history with no free parameter (the ratio of the two excess temperatures is a function of Fo only: 0.86, 0.75, 0.71 at 300, 600, 1200 s for this tile). If the mid-depth reading misses that prediction, the properties are not constant.
- **Does the result survive a change of grid, step or solver?** FTCS at 51 and 101 points, and the series, agree to 0.2 K; leapfrog at s = 0.5, 0.8, 1.0 and d'Alembert agree to 0.03 mm; Jacobi on three grids and the sinh series agree to 0.05 K. State the smallest grid you would accept for each.
- **Which decision would change if the uncertainty doubled?** If α were uncertain by ±40 %, the required thickness would be uncertain by ±20 % (12 mm), and the 50 mm tile could not be ruled in or out without a test; you would ask for a measured α at temperature before cutting tiles.

## Check yourself

1. Write the heat equation and give the units of α. What single dimensionless group governs the slab solution?
2. For a slab with both faces at zero, what are the eigenfunctions and decay rates? How much faster does mode 3 decay than mode 1?
3. For a slab with one face fixed and one insulated, what are the eigenvalues λ_n, and why are they quarter-waves rather than half-waves?
4. Give the wave speed of a string and the frequency of its n-th mode. Why does plucking at L/4 produce no fourth harmonic?
5. State d'Alembert's solution for a string released from rest, and explain the role of the odd periodic extension.
6. Why does the separation constant change sign between the heat and Laplace problems, and what does that imply about where boundary data are needed?
7. State the stability conditions for FTCS and for leapfrog, and say what each means physically.
8. The 50 mm tile's back face reaches 595 K after 1200 s. Without the code, estimate the back-face temperature of a 100 mm tile after 4800 s.
