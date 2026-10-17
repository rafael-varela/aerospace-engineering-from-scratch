# P03.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and given to the figures shown.

## 1. The numbers that size a run (s = 15: L = 4.50 mm)

**Assumptions.** Constant properties; a uniform grid of n cells; FTCS limited by r ≤ ½.

**Derivation.** (a) α = 167/(2700 × 896) = **6.903 × 10⁻⁵ m²/s**; τ = L²/α = (4.5 × 10⁻³)²/6.903 × 10⁻⁵ = **0.293 s**; the 10 s run is **34 diffusion times** long (the skin is essentially at steady state after the first second).

(b) n = 50: Δx = 0.090 mm, Δt_max = Δx²/(2α) = **58.7 µs**, steps = 10/58.7 × 10⁻⁶ = **170 448** (round up). n = 100: Δt_max = **14.7 µs**, **681 789** steps. Steps ∝ n² (Δt_max ∝ Δx²).

(c) r = αΔt/Δx² = 6.903 × 10⁻⁵ × 10⁻³/(9 × 10⁻⁵)² = **8.52**; **10 000 steps**. Crank–Nicolson, for its second-order time error; but at r = 8.5 its grid-scale factor is g(π) = (1 − 17)/(1 + 17) = −0.89, so if the initial condition has a step (a sudden face temperature) the first steps will ring, and two BTCS half-steps at the start are the check-and-fix.

**Check.** Steps for n = 100 = 4 × steps for n = 50 ✓ (681 789/170 448 = 4.00).

**Common mistakes.** Using Δx = L/(n + 1); forgetting the factor 2 in Δx²/(2α); quoting r without stating Δx.

## 2. Truncation error and amplification factors

(a) With u and its derivatives at (x_j, t^m): u_{j±1}^m = u ± Δx u_x + ½Δx²u_xx ± ⅙Δx³u_xxx + (1/24)Δx⁴u_xxxx + O(Δx⁵), so D₂u = u_xx + (Δx²/12)u_xxxx + O(Δx⁴). And u_j^{m+1} = u + Δt u_t + ½Δt²u_tt + O(Δt³). FTCS divided by Δt reads (u_t + ½Δt u_tt) − α(u_xx + (Δx²/12)u_xxxx) = O(Δt², Δx⁴). The exact solution has u_t = αu_xx and u_tt = αu_xxt = α²u_xxxx, so the residual is τ = ½Δt α²u_xxxx − (αΔx²/12)u_xxxx = αΔx²(r/2 − 1/12)u_xxxx, which is (4) up to the sign convention; it vanishes at **r = 1/6**. For BTCS, expand about t^{m+1}: u_j^m = u − Δt u_t + ½Δt²u_tt − …, and the Δt term enters with the opposite sign: τ = −αΔx²(r/2 + 1/12)u_xxxx. Both terms have the same sign for r > 0, so no cancellation is possible.

(b) Substitute u_j^m = gᵐe^{ijθ}. FTCS: g = 1 + r(e^{−iθ} − 2 + e^{iθ}) = 1 − 4r sin²(θ/2) = 1 − 4rs. Stability needs −1 ≤ 1 − 4rs for s = 1 (the worst case), i.e. r ≤ ½; the other side, 1 − 4rs ≤ 1, is automatic. BTCS: g(1 + 4rs) = 1, g = 1/(1 + 4rs) ∈ (0, 1]. Crank–Nicolson: g(1 + 2rs) = 1 − 2rs, g = (1 − 2rs)/(1 + 2rs); with z = 2rs ≥ 0, |1 − z| ≤ 1 + z always, so |g| ≤ 1, with g → −1 as z → ∞.

(c) e^{−z} = 1 − z + z²/2 − z³/6 + …. FTCS: 1 − z, below e^{−z} by z²/2: every mode decays **too fast**. BTCS: 1/(1 + z) = 1 − z + z² − z³ + …, above e^{−z} by z²/2: every mode decays **too slowly**, so a transient lingers and a threshold is reached **late** (the +82 ms of the lesson). Crank–Nicolson: (1 − z/2)/(1 + z/2) = 1 − z + z²/2 − z³/4 + …, agreeing to z² with an error −z³/12 relative to e^{−z}: a slight over-damping for small z, but for large z the factor tends to −1 while e^{−z} → 0, which is the ringing.

**Check.** At r = 1/6 and small θ, 1 − 4rs ≈ 1 − (2/3)s versus e^{−(2/3)s}: the Δx² term of the exponent matches because the next term of sin² and of the exponential cancel, which is the fourth-order statement in Fourier form.

**Common mistakes.** Expanding BTCS about t^m and getting the FTCS sign; forgetting that stability is required for *all* θ, so the worst case s = 1 decides.

## 3. An insulated face by a ghost node

**Assumptions.** Constant foam properties; outer face perfectly insulated (a real tank has convection and radiation to ambient, which would keep the outer face far warmer); inner face instantly at 20 K.

(a) α = 0.03/(40 × 1500) = **5.0 × 10⁻⁷ m²/s**; τ = L²/α = (0.025)²/5 × 10⁻⁷ = **1250 s**. One hour is 2.9 diffusion times: long enough that the foam is nearly at its steady (uniform 20 K) state.

(b) The ghost node u_{n+1} = u_{n−1} enforces (u_{n+1} − u_{n−1})/(2Δx) = 0. The FTCS stencil at j = n becomes u_n^{m+1} = u_n + r(u_{n−1} − 2u_n + u_{n−1}) = **u_n + 2r(u_{n−1} − u_n)**. The update matrix on nodes 1 … n has interior rows (r, 1 − 2r, r) and a last row (2r, 1 − 2r). Gershgorin: eigenvalues lie in discs centred at 1 − 2r with radii 2r for every row, i.e. in [1 − 4r, 1]. The lower end is ≥ −1 iff **r ≤ ½**: the Neumann boundary does not change the limit. (The discrete-cosine eigenvectors of the Neumann Laplacian give the exact eigenvalues 1 − 4r sin²(θ_l/2) with θ_l = (l − ½)π/n, confirming it.) A run at r = 0.55 on 25 cells reaches 10²⁵⁶ by one hour.

(c) FTCS at r = 0.4 to 3600 s:

| n | T(L/2) [K] | T(L) [K] |
|---|---|---|
| 25 | 20.20 | 20.29 |
| 50 | 20.21 | 20.29 |
| 100 | 20.21 | 20.29 |

Converged to 0.01 K between n = 50 and 100.

(d) Series with λ_m = (2m + 1)π/(2L): T(L/2) = **20.21 K**, T(L) = **20.29 K**: agreement to the printed precision. The leading term alone, (300 − 20)(4/π) sin(π/4) e^{−απ²t/(4L²)} = 356.4 × 0.7071 × e^{−7.106} = 0.21 K at the mid-plane, is all that survives after 2.9 diffusion times.

**Common mistakes.** Writing the boundary row with r instead of 2r (that imposes a different, wrong, boundary condition and is still stable); forgetting that the ghost-node form is second-order accurate only because the node sits exactly on the boundary.

## 4. The bulkhead with a hotter duct (s = 15: T_hot = 750 K)

**Assumptions.** As in the lesson; constant α = 2.822 × 10⁻⁶ m²/s; n = 100.

(a) Steady mid-plane (750 + 300)/2 = **525 K** > 450 K, so the limit is reached. `time_to_limit_exact` → **4.812 s**.

(b)

| scheme | Δt | time [s] | error [ms] | steps |
|---|---|---|---|---|
| FTCS, r = 0.45 | 1.595 ms | 4.811 | −1.0 | 3017 |
| Crank–Nicolson | 0.1 s | 4.812 | −0.1 | 49 |
| BTCS | 0.1 s | 4.879 | +67.0 | 49 |

(c) q(4.812 s) = **153.0 kW/m²**; steady q = k(T_hot − T_cold)/L = 7 × 450/0.01 = **315.0 kW/m²**; **48.6 %** of the steady coolant load has arrived when the bracket reaches its limit.

(d) Prediction: BTCS is first order in Δt, so halving the step halves the error to about +33 ms. Run: **4.846 s, +33.5 ms, 97 steps** ✓.

**Check.** The hotter duct reaches the limit sooner than the lesson's 700 K case (4.81 s against 5.85 s) because the mid-plane has further to go in steady state (525 K against 500 K) and the same time constant.

**Common mistakes.** Reading the limit time at the end of the step that crosses rather than interpolating (an error of up to Δt); comparing BTCS with Crank–Nicolson rather than with the series.

## 5. Buying accuracy with order (L = 10 mm, α = 2.822 × 10⁻⁶ m²/s, t_f = 3.591 s)

(a) n = 200, Δt = 0.718, 0.359, 0.180, 0.0898 s:

| scheme | errors | observed orders |
|---|---|---|
| BTCS | 3.40 × 10⁻², 1.77 × 10⁻², 9.01 × 10⁻³, 4.55 × 10⁻³ | 0.94, 0.97, 0.99 |
| Crank–Nicolson | 1.23 × 10⁻³, 3.07 × 10⁻⁴, 7.67 × 10⁻⁵, 1.92 × 10⁻⁵ | 2.00, 2.00, 2.00 |

(b) Δt* ≈ 0.718 × (1.23 × 10⁻³/3.40 × 10⁻²) = **0.026 s**, about **28 times** more steps than Crank–Nicolson at t_f/5. Run with 138 steps (Δt = 0.02602 s): error **1.33 × 10⁻³**, within 8 % of the target (the first-order law is not quite asymptotic at Δt = t_f/5, order 0.94).

(c) FTCS at n = 80: r = 0.4 gives 6.62 × 10⁻⁵, r = 1/6 gives 1.62 × 10⁻⁹: a factor of **40 900**. Orders 2.01, 2.00, 2.00 at r = 0.4 and **4.01, 4.00, 4.00** at r = 1/6. From (4), the coefficient (1/12 − r/2) of the Δx² term vanishes at r = 1/6, and the next term is O(Δx⁴).

(d) The cancellation needs (i) constant α, (ii) a uniform grid, and (iii) u_tt = α²u_xxxx, i.e. the solution actually satisfies the constant-coefficient heat equation. Temperature-dependent conductivity breaks (i) and (iii) at once: α varies across the plate and u_t = (k(T)T_x)_x/(ρc_p) no longer gives u_tt = α²u_xxxx. The grid can be kept uniform, so (ii) survives longest.

**Common mistakes.** Choosing Δt values that do not divide t_f and comparing at different final times; expecting the r = 1/6 miracle to survive a Dirichlet step initial condition (it does in exact arithmetic, but the fourth-order term is multiplied by u_xxxx, which is enormous at a step).

## 6. A colleague's time step

(a) r = αΔt/Δx² = 2.822 × 10⁻⁶ × 1/(10⁻⁴)² = **282**; g(π) = (1 − 564)/(1 + 564) = **−0.9965**; N = ln 0.01/ln 0.9965 = **1299 steps** to decay the sawtooth to 1 %. The first 5 s take 5 steps. The grid-scale content of the step initial condition is therefore essentially undamped throughout the run, flipping sign every step.

(b) At t = 5 s:

| scheme, Δt = 1 s | max error vs series [K] | max |second difference| [K] |
|---|---|---|
| Crank–Nicolson | 329.8 | 423.9 |
| BTCS | 13.7 | 0.1 |

Crank–Nicolson is stable (nothing grows) and useless (a 330 K error on a 400 K range); BTCS is smooth and 14 K off through first-order under-damping.

(c) Two BTCS steps of 0.5 s, then four Crank–Nicolson steps of 1 s: error **1.24 K**, second difference **0.08 K**. The two implicit half-steps multiply the sawtooth by (1/(1 + 2 × 141))² = 1.2 × 10⁻⁵ and the rest of the run keeps second-order accuracy. (Four Crank–Nicolson steps of 1 s still carry a 1.2 K first-order-looking error from the start-up and the large step; at Δt = 0.1 s the same remedy gives 0.01 K.)

(d) Mean k over 300–700 K ≈ (7 + 10.5)/2 = 8.75 W/(m·K), so α is 1.25 times larger. With fixed boundary and initial temperatures the solution depends on x and αt only, so every time scales as 1/α: the lesson's 5.845 s becomes **4.68 s, 1.17 s earlier** (for s = 15, 4.81 s becomes 3.85 s, 0.96 s earlier). Every numerical error in this problem set is below 0.1 s once a sane step is used, and even the colleague's BTCS run is 14 K off, a fraction of a second in time. The modelling error is ten times the worst numerical one. The colleague should put the effort into k(T) and the hot-face boundary condition, and run Crank–Nicolson at Δt ≈ 0.1 s with a BTCS start.

**Common mistakes.** Calling the Crank–Nicolson result "unstable" (it is not: |g| ≤ 1 and nothing grows; it is inaccurate); believing a finer grid will fix it (r grows, g(π) → −1 faster).

## Check-yourself answers

1. r = αΔt/Δx²; FTCS requires r ≤ ½. For α = 6.9 × 10⁻⁵ m²/s and Δx = 0.1 mm, Δt_max = Δx²/(2α) = 72 µs.
2. FTCS 1 − 4rs, BTCS 1/(1 + 4rs), Crank–Nicolson (1 − 2rs)/(1 + 2rs), with s = sin²(θ/2); the exact factor is e^{−4rs}. FTCS is below it (over-damps), BTCS above it (under-damps).
3. Because |1 − z| ≤ |1 + z| for z ≥ 0, so |g| ≤ 1 for every r; but g(π) = (1 − 2r)/(1 + 2r) → −1, so the grid-scale content of a step is not damped, only flipped in sign each step. g at θ = π is the number to look at.
4. The Δt and Δx² terms of the truncation error, αΔx²(1/12 − r/2)u_xxxx, cancel, leaving O(Δx⁴). It requires constant α, a uniform grid and u_tt = α²u_xxxx, so it is a property of this one linear problem, not of the scheme.
5. For Dirichlet ends the eigenvectors of the update matrix are the discrete sines sin(jθ_l), θ_l = lπ/n, and its eigenvalues are the von Neumann factors g(θ_l): the two analyses agree exactly. They disagree for boundary conditions or coefficients that break the Fourier-mode assumption; the matrix analysis is then the valid one.
6. Against the exact PDE solution: the total (space + time) error. Against the exact semi-discrete solution: the time error alone. The difference between the two references is the spatial error of the grid.
7. Yes: BTCS's factor 1/(1 + z) exceeds e^{−z}, so every mode decays too slowly and a threshold is reached late. The error is first order in Δt, so Δt = 0.05 s gives about 41 ms late (the code gives 40.9 ms).
8. When the stability-limited explicit step is far smaller than the step the accuracy requires, i.e. when r_accuracy ≫ ½. For the bulkhead on 100 cells: 3665 FTCS steps × 5 flops against 59 Crank–Nicolson steps × 8 flops, a factor of about 40.
