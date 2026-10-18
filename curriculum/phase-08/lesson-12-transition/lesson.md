# P08.12  Transition

> A laminar layer does not turn turbulent all at once: small waves start to grow at a critical Reynolds number, and transition comes where they have grown by a factor of about e⁹.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Orr–Sommerfeld equation for small wavy disturbances of a parallel shear flow, state Rayleigh's inflection-point criterion, and explain why viscosity destabilises plane Poiseuille and Blasius flows. 2. Implement and test a Chebyshev collocation solver for the Orr–Sommerfeld eigenvalue problem, show spectral convergence, and reproduce the published critical point Re_c = 5772.22, α_c = 1.02056 and Orszag's eigenvalue at Re = 10 000. 3. Apply the e^N method to a flat plate and to wedge flows with a supplied approximate-fit stability chart, and predict where a wing's boundary layer becomes unstable and where it turns turbulent. 4. Judge how the transition point depends on free-stream turbulence, the pressure gradient and the shape factor, and when the e^N method should not be used. |
| **Prerequisites** | `P08.10` (Blasius layer, θ, H, Re_θ = 0.664√Re_x), `P08.11` (Thwaites's method, wedge flows, H(λ)), `P08.06` (Navier–Stokes, plane Poiseuille flow), `P01.16` (eigenvalues; generalised eigenproblems), `P03.08` (dense linear algebra), `P03.02` (root finding: Brent's method), `P02.02` (boundary-value problems), `P02.11` (linear PDEs and normal modes). Skill check: you can compute the eigenvalues of a 3 × 3 matrix with `scipy.linalg.eig` and say what a complex eigenvalue means for a solution e^(λt). |
| **Engineering question** | How far back can the boundary layer on the light aeroplane's wing stay laminar at cruise (3000 m, 60 m/s, chord 1.5 m), treating its upper surface as a flat plate, and how sensitive is that to the turbulence of the air? Answer: on the e^N chart, small waves start to grow 41 mm from the leading edge (Re_x = 1.33 × 10⁵; the Orr–Sommerfeld critical point of the Blasius profile itself is nearer, at about 28 mm); the e⁹ method puts transition at Re_x = 2.79 × 10⁶, 0.865 m or 58 % of the chord, in quiet air; at a turbulence intensity of 0.3 % (N_crit = 5.5) it moves forward to 28 % of the chord. The eigenvalue solver behind the method reproduces the published critical Reynolds number 5772.22 to 3 × 10⁻⁷. |
| **Assumptions and limits** | Linear stability (infinitesimal disturbances); two-dimensional waves (Squire's theorem makes them the first to go unstable in a parallel flow); parallel base flow (exact for a channel, approximate for a growing boundary layer); incompressible. The e^N part uses an **approximate fit** (Drela and Giles 1987) to Orr–Sommerfeld results for Falkner–Skan profiles as a supplied chart, and Mack's **correlation** for N_crit; it covers only transition by slowly amplified Tollmien–Schlichting waves, not bypass, roughness, crossflow or attachment-line transition. |
| **Analytical work** | Linearisation and the normal-mode ansatz; the Orr–Sommerfeld equation and its boundary conditions; Rayleigh's criterion; the meaning of α, c_r, c_i and the neutral curve; the e^N envelope N = (dN/dRe_θ)(Re_θ − Re_θ0); the amplitude ratio e⁹ ≈ 8100; the critical point of a channel in physical units (Problem 3). |
| **Computational deliverable** | `starter/main.py`: `cheb`, `os_matrices`, `os_eigenvalues`, `least_stable`, `neutral_reynolds`, `critical_point`, `en_envelope`, `transition_reynolds`, `flat_plate_transition` (`isa_nu`, `drela_giles_dn_dretheta`, `drela_giles_retheta0`, `thwaites_wedge`, `mack_ncrit`, `observed_decay` and `main` are given). Tests prove: the Chebyshev matrix differentiates polynomials exactly and converges spectrally; the solver reproduces Orszag's published eigenvalue to 10⁻⁷ and converges by more than 20× per ten points; Poiseuille flow is stable below Re_c and unstable above it; viscous decay scales as 1/Re at low Re; the neutral point has zero growth; the critical point matches the published value to 10⁻⁵; the envelope is zero before the critical Re_θ and linear after it; the flat-plate transition point is self-consistent and unit-invariant; a higher N_crit and a favourable gradient delay transition. |
| **Practical activity** | **computational experiment**: spectral convergence of the eigenvalue solver, the neutral curve and the critical point, the trustworthy part of the spectrum, and the sensitivity of the wing's predicted transition point to N_crit (free-stream turbulence), the pressure gradient and the shape factor. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `transition.png`, the derivation of the Orr–Sommerfeld equation and Rayleigh's criterion, the experiment tables, and a half-page note stating the wing's predicted instability and transition points, their sensitivity, and the transition routes the prediction ignores. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: the Orr–Sommerfeld equation is derived with its assumptions; all public tests pass; spectral convergence and the critical point are demonstrated; the note gives x₀ and x_tr for the wing, says that the chart is an approximate fit, and names at least two transition routes it does not model. |
| **References** | `@ntrs-malik-1987-stability-theory-lfc` (https://ntrs.nasa.gov/citations/19900003197): NEED FOR TRANSITION PREDICTION METHODS; PREDICTION METHODS BASED UPON BOUNDARY-LAYER STABILITY THEORY (Orr-Sommerfeld equation; e^N method; N-factor envelopes) (supplementary-free: the method as used in laminar-flow design). `@leishman-intro-aerospace-flight-vehicles` "Boundary Layer Flows" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/boundary-layers/): Laminar & Turbulent Boundary Layers; Surface Roughness; Laminar Separation Bubbles (required-free; link only). `@libretexts-liburdy-intermediate-fluid-mechanics` Ch. 10 "Introduction to Turbulence Effects": Initial Laminar Starting Length (required-free). `@mit-ocw-2-27-turbulent-flow-transport` "Concepts in Turbulence": 3.1 Comments on laminar flow, its stability, and the transition to turbulent flow (supplementary-free). Optional textbooks for deeper study: Schlichting and Gersten, *Boundary-Layer Theory*, and White, *Viscous Fluid Flow* (paid; library copy; see `references/phase-08.md`). |
| **Extensions** | ★ Solve Orr–Sommerfeld for the Blasius profile itself (map Chebyshev points onto [0, η_max]) and find its critical Re_δ*; compare it with the chart's Re_θ0 × H. ★ Compute spatial growth rates (complex α at real frequency) with Gaster's transformation α_i ≈ −ω_i/c_g and integrate N along the plate yourself. ★ Add Squire's transformation and verify numerically that an oblique wave at angle φ is neutral at Re/cos φ. |
| **Time** | unmeasured |

## 1. The engineering problem

`P08.10` found that a laminar layer over the whole 1.5 m chord of the light aeroplane's wing would carry 29.7 N of friction drag at cruise; `P08.14` will find 166 N for a layer that is turbulent from the leading edge. The real wing lies between, and where it lies depends on a single position: the **transition point**, where the laminar layer becomes turbulent. Laminar-flow aerofoils, sailplanes and modern airliner nacelles are designed around it, and a wind-tunnel test is only useful if its transition point is where the aircraft's would be.

Transition is not predicted by solving for turbulence. It starts with **instability**: above a critical Reynolds number, small wavy disturbances in the layer (Tollmien–Schlichting waves) grow instead of decaying. They grow slowly, by a factor that can be computed from linear theory, and turbulence appears where they have become large. The **e^N method**, used in aircraft design for decades, says: compute the amplification e^N of the most dangerous wave as it travels downstream, and predict transition where N reaches an empirical N_crit, about 9 in quiet flight.

This lesson builds the stability calculation for the flow where it is cleanest to check, plane Poiseuille flow between two plates, whose critical Reynolds number 5772.22 is a published benchmark (Orszag 1971), and then applies the e^N idea to the wing using a supplied chart of growth rates for boundary-layer profiles. Success criteria, written first: the eigenvalue solver converges spectrally to Orszag's published eigenvalue at α = 1, Re = 10⁴ and finds the critical point to five figures; the flow is stable below Re_c and unstable above; the e^N envelope is consistent with its chart; the predicted transition point responds in the right direction to turbulence and pressure gradient.

## 2. Understand and derive

### 2.1 Small disturbances of a parallel flow

Take a steady parallel flow U(y) between walls at y = ±1 (lengths scaled by the half-height h, velocities by the centreline speed, Re = U_c h/ν). Add a small disturbance and linearise the Navier–Stokes equations (`P08.06`) about U: products of disturbances are dropped. Eliminating the pressure by cross-differentiation and seeking **normal modes**,

ψ′ = φ(y) e^(iα(x − ct)),

a wave of wavenumber α (wavelength 2π/α) travelling at phase speed c_r and growing at the rate αc_i (c = c_r + ic_i), gives (Problem 2a)

**(U − c)(φ″ − α²φ) − U″φ = (φ⁗ − 2α²φ″ + α⁴φ)/(iαRe),  φ = φ′ = 0 at y = ±1.**   (1)

This is the **Orr–Sommerfeld equation**. For given α and Re it has nontrivial solutions only for particular c: an **eigenvalue problem**. If any eigenvalue has c_i > 0 the flow is linearly unstable. The curve c_i(α, Re) = 0 is the **neutral curve**, and its lowest point is the **critical Reynolds number**. Squire's theorem (Extension 3) shows that two-dimensional waves become unstable first, so (1) is enough to find Re_c.

Two physical points:

- **Viscosity destabilises.** Without the right-hand side (Re → ∞, Rayleigh's equation) an unstable mode needs an inflection point in U (Problem 2b). Poiseuille flow (U = 1 − y², U″ = −2) and the Blasius profile have none, so they are inviscidly stable; their Tollmien–Schlichting waves owe their growth to the viscous terms, which is why they exist only in a window of Re and α.
- **Channel flow transition is subcritical.** Real channel flows become turbulent at much lower Re than 5772 because large disturbances bypass the linear stage (Problem 3). Boundary layers in quiet air are the case where the slow linear stage dominates, which is why e^N works there.

### 2.2 Chebyshev collocation

The eigenfunctions are smooth, so represent φ by its values at the **Chebyshev–Gauss–Lobatto points** y_j = cos(πj/N), j = 0…N, which cluster near the walls where the viscous layers of φ are thin. The **differentiation matrix** D maps values to the derivative of the interpolating polynomial (Trefethen, *Spectral Methods in MATLAB*): D_ij = (c_i/c_j)(−1)^(i+j)/(y_i − y_j) for i ≠ j, with c₀ = c_N = 2 and c_j = 1 otherwise, and the diagonal chosen so that each row sums to zero. Its error falls faster than any power of 1/N for a smooth function (**spectral accuracy**).

Boundary conditions: φ = 0 at the walls is imposed by keeping only the N − 1 interior points. φ′ = 0 needs care in the fourth derivative; write φ = (1 − y²)g, which vanishes at the walls, and build the fourth-derivative matrix on g:

D4 = [diag(1 − y²) D⁴ − 8 diag(y) D³ − 12 D²] diag(0, 1/(1 − y_j²), 0),

which builds in φ′ = 0. Then (1) is the **generalised matrix eigenproblem**

**A φ = c B φ,  A = diag(U)(D2 − α²) − diag(U″) − (D4 − 2α²D2 + α⁴)/(iαRe),  B = D2 − α²,**   (2)

solved by `scipy.linalg.eig(A, B)`. With N = 50–60 the least-stable eigenvalue is good to eight digits.

### 2.3 The e^N envelope for a boundary layer

On a growing boundary layer the local Reynolds number Re_θ = Ueθ/ν rises downstream. A wave of fixed frequency is amplified while it is inside the neutral curve; its amplitude grows as A = A₀e^N with N = ∫(spatial growth rate)dx. The **envelope** over all frequencies, N(x), is what matters. For the Falkner–Skan family the envelope is almost linear in Re_θ beyond a critical value, with a slope that depends only on H. This lesson uses Drela and Giles's (1987) **approximate fit** to those Orr–Sommerfeld results, supplied in `starter/main.py` as a chart in formula form:

dN/dRe_θ = 0.01√([2.4H − 3.7 + 2.5 tanh(1.5H − 4.65)]² + 0.25),
log₁₀ Re_θ0 = (1.415/(H − 1) − 0.489) tanh(20/(H − 1) − 12.9) + 3.295/(H − 1) + 0.44,

**N = (dN/dRe_θ)(Re_θ − Re_θ0) for Re_θ > Re_θ0, and 0 before.**   (3)

Transition is predicted where N = N_crit. Mack's **correlation** relates N_crit to the free-stream turbulence intensity Tu (as a fraction): **N_crit = −8.43 − 2.4 ln Tu**, which gives 9.0 at Tu = 0.07 %, typical of quiet flight and the best low-turbulence tunnels.

### 2.4 The wing

For a flat plate (Blasius, H = 2.591, Re_θ = 0.664√Re_x) the chart gives Re_θ0 = 242 and dN/dRe_θ = 0.01039, so

- **instability** (the chart's onset) at Re_x0 = (242/0.664)² = 1.33 × 10⁵: x₀ = Re_x0 ν/U = **41.2 mm** from the leading edge at 3000 m and 60 m/s (ν = 1.8630 × 10⁻⁵ m²/s);
- **transition** at Re_θ = 242 + 9/0.01039 = 1108, Re_x = 2.79 × 10⁶: **x_tr = 0.865 m, 58 % of the 1.5 m chord**.

The 41 mm is the onset of the fitted envelope, not the linear-stability critical point. The Orr–Sommerfeld critical point of the Blasius profile itself is Re_δ* ≈ 520 (Extension 1 asks you to compute it), which with δ* = 1.7208x/√Re_x is Re_x = (520/1.7208)² = 9.13 × 10⁴, **x ≈ 28 mm**. The chart's Re_θ0 × H = 627 differs from 520 because the chart is an approximate fit to the amplification envelope, not a stability calculation. Either way the waves start to grow within a few centimetres of the leading edge, and the transition estimate below uses the chart. The waves are unstable over 95 % of the laminar run, and need about 0.8 m to grow by e⁹ = 8100. Everything therefore depends on N_crit and on the growth rate, which depends on H: part 4 measures both sensitivities.

### 2.5 Where the method stops being valid

- **Bypass routes**: high free-stream turbulence, wakes, roughness, crossflow on swept wings, attachment-line contamination: Problem 6.
- **Non-parallel effects**: a boundary layer grows, so the parallel-flow form (1) is an approximation for it, poorest at low Reynolds number where the layer grows fastest relative to its thickness.
- **The chart is a fit**: H enters exponentially, so a 0.7 % error in H moves x_tr by 17 % (part 4).
- **N_crit is calibrated, not derived**: it encodes the disturbance environment and the breakdown process, neither of which linear theory describes.

## 3. Build the computational model

Open `starter/main.py`. `isa_nu`, the two chart functions `drela_giles_dn_dretheta` and `drela_giles_retheta0` (labelled approximate fit), `thwaites_wedge`, `mack_ncrit`, `observed_decay` and `main` are given.

### Step 1: the differentiation matrix

```python
def cheb(N):
    x = np.cos(np.pi * np.arange(N + 1) / N)
    c = np.ones(N + 1); c[0] = c[-1] = 2.0; c *= (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(c, 1.0 / c) / (X - X.T + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x
```

Test it before anything else: `D @ x**k` must equal `k * x**(k-1)` to round-off for k ≤ N.

### Step 2: the Orr–Sommerfeld matrices

`os_matrices(N, alpha, Re)` builds D², D³, D⁴ from `cheb`, applies the (1 − y²) trick of 2.2, keeps the interior block, and returns A and B of (2) as complex arrays. `os_eigenvalues` calls `scipy.linalg.eig(A, B, right=False)`, drops infinite eigenvalues and sorts by decreasing c_i; `least_stable` returns the first.

### Step 3: neutral curve and critical point

`neutral_reynolds(alpha, N)` finds the root of c_i(Re) = 0 on the lower branch with `scipy.optimize.brentq` on [4000, 10 000] (raising `ValueError` if there is no sign change). `critical_point(N)` minimises the neutral Re over α with `minimize_scalar(..., method="bounded")` on [0.95, 1.08]: a root finder inside an optimiser, each evaluation an eigenvalue solve.

### Step 4: the e^N envelope

`en_envelope(Re_theta, H)` is (3). `transition_reynolds(H, k, N_crit)` inverts it for a similar layer with Re_θ = k√Re_x; `flat_plate_transition(U, nu, N_crit)` uses the Blasius values and converts to metres.

Running `main()` prints

```
Orr-Sommerfeld, plane Poiseuille flow, alpha = 1, Re = 10 000 (published c = 0.23752649 + 0.00373967i):
  N =  20: c = 0.23879009 + 0.00630869i, |c - published| = 2.86e-03
  N =  30: c = 0.23763053 + 0.00380254i, |c - published| = 1.22e-04
  N =  40: c = 0.23752655 + 0.00374080i, |c - published| = 1.13e-06
  N =  50: c = 0.23752649 + 0.00373967i, |c - published| = 3.80e-09
  N =  60: c = 0.23752649 + 0.00373967i, |c - published| = 1.25e-09
  N =  80: c = 0.23752649 + 0.00373967i, |c - published| = 2.06e-09
critical point (N = 60): Re_c = 5772.22, alpha_c = 1.02054; published 5772.22, 1.02056
  relative error in Re_c 3.1e-07
  N = 30: neutral Re at alpha = 1.02056: 5783.566
  N = 40: neutral Re at alpha = 1.02056: 5772.217
  N = 50: neutral Re at alpha = 1.02056: 5772.222

flat plate, Blasius H = 2.591 (approximate-fit stability chart): dN/dRe_theta = 0.01039, Re_theta0 = 242.0
  instability: Re_x0 = 1.327e+05; transition (N = 9): Re_theta = 1108, Re_x = 2.785e+06
  wing at 3000 m, 60.0 m/s (nu = 1.8630e-05): instability at x = 41.2 mm, transition at x = 0.865 m = 58 % of the 1.5 m chord
  (the chart's onset; the Orr-Sommerfeld critical point of the Blasius profile, Re_delta* = 520, is at Re_x = 9.13e+04, x = 28.4 mm)
  sensitivity to N_crit (Mack: Tu -> N_crit):
    Tu = 0.07 %: N_crit =  9.00, Re_x_tr = 2.79e+06, x_tr = 0.866 m (58 % chord)
    Tu = 0.10 %: N_crit =  8.15, Re_x_tr = 2.39e+06, x_tr = 0.742 m (49 % chord)
    Tu = 0.30 %: N_crit =  5.51, Re_x_tr = 1.35e+06, x_tr = 0.420 m (28 % chord)
    Tu = 1.00 %: N_crit =  2.62, Re_x_tr = 5.54e+05, x_tr = 0.172 m (11 % chord)
    Tu = 2.00 %: N_crit =  0.96, Re_x_tr = 2.53e+05, x_tr = 0.079 m (5 % chord)
  pressure gradient (wedge flows Ue ~ x^m, Thwaites H; N_crit = 9):
    m = -0.05: H = 2.753, Re_theta0 = 108, Re_x_tr = 6.26e+05
    m = -0.02: H = 2.650, Re_theta0 = 157, Re_x_tr = 1.46e+06
    m = +0.00: H = 2.610, Re_theta0 = 206, Re_x_tr = 2.27e+06
    m = +0.05: H = 2.544, Re_theta0 = 394, Re_x_tr = 5.77e+06
    m = +0.10: H = 2.502, Re_theta0 = 647, Re_x_tr = 1.2e+07
    m = +0.20: H = 2.452, Re_theta0 = 1126, Re_x_tr = 3.13e+07
```

and writes `transition.png`: the error against Orszag's eigenvalue on a log scale against N (a straight line falling until round-off: spectral convergence) and the N(x) envelope along the wing for the flat plate and two wedge flows, with N = 9 marked.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-12-transition`. The tests:

1. **Chebyshev**: endpoints ±1; D·1 = 0; exact derivatives of x, x³, x⁷, x^N; for e^(sin x) the error falls by more than 50× from N = 8 to 12 and again to 16.
2. **Orr–Sommerfeld**: (N − 1)² complex A, real B; Orszag's eigenvalue to 10⁻⁷ with c_i > 0; errors falling by more than 20× per 10 points from N = 20 to 40; eigenvalues sorted with 0 < c_r < 1 for the ten least stable; stable at Re = 5000 for α = 0.8, 1.0, 1.2, unstable at 6500 and stable at 5600 for α_c; αc_iRe the same at Re = 1 and 0.5 (viscous decay ∝ 1/Re); zero growth at the neutral point and `ValueError` with no root in the bracket; the critical point within 10⁻⁵ of 5772.22 and 10⁻³ of 1.02056, with higher neutral Re on both sides.
3. **e^N**: N = 0 up to Re_θ0 and slope dN/dRe_θ after it; the flat-plate transition is self-consistent (N = 9 at Re_θ,tr), lies between Re_x = 10⁶ and 5 × 10⁶, and is unchanged by a change of units; x_tr increases with N_crit; Mack gives 9.0 at Tu = 0.07 %; m = −0.03, 0, +0.05 give increasing Re_x,tr.
4. `main` writes the figure.

### Compare

There is no stability solver in SciPy to compare against, so the comparison is with published numbers: Orszag's eigenvalue (eight digits at N ≥ 50) and critical point (Re_c = 5772.22 to 3 × 10⁻⁷, α_c 1.02054 against 1.02056). A second route inside the code is the convergence in N itself: the neutral Re at α = 1.02056 is 5783.566, 5772.217, 5772.222 at N = 30, 40, 50, so the answer is a property of the equation, not of the grid. The e^N chart cannot be checked against an exact solution here; it is an approximate fit, and the lesson treats its output as an estimate whose sensitivity is part of the answer.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Spectral convergence.** |c − c_Orszag| for N = 20–80 at α = 1, Re = 10⁴. Predict the shape of the error curve on a log scale. *Finding:* 2.86 × 10⁻³, 1.22 × 10⁻⁴, 1.13 × 10⁻⁶, 3.8 × 10⁻⁹, 1.25 × 10⁻⁹, 2.1 × 10⁻⁹ for N = 20, 30, 40, 50, 60, 80: a factor of about 23, 110 and 300 per 10 points, then a floor near 10⁻⁹ set by the published value's eight decimals and round-off.
2. **Neutral curve.** Lower-branch neutral Re at α = 0.90, 0.95, 1.00, 1.05, 1.08 (N = 60). Predict where the minimum lies. *Finding:* 6965, 6208, 5815, 5890, 6478; minimum Re_c = 5772.22 at α = 1.0205. At N = 30 the neutral Re at α_c is 5783.6, already within 0.2 %; at N = 40, 5772.217.
3. **Which eigenvalues to trust.** The whole spectrum at α = 1, Re = 10⁴ for N = 60 and N = 100. Predict how many agree. *Finding:* one growing mode; only 7 of the 59 eigenvalues agree to 10⁻⁶ between the two N, all with c_i ≥ −0.183. Strongly damped modes and near-coincident pairs are resolution-dependent; only the least-stable few are physical at this N.
4. **N_crit (free-stream turbulence).** x_tr for Tu = 0.07, 0.1, 0.3, 1, 2 % on the wing. Predict the trend. *Finding:* N_crit = 9.00, 8.15, 5.51, 2.62, 0.96; x_tr = 0.866, 0.742, 0.420, 0.172, 0.079 m (58, 49, 28, 11, 5 % of the chord). A tunnel with 0.3 % turbulence halves the laminar run.
5. **Pressure gradient and shape factor.** Wedge flows m = −0.05 to +0.2 (Thwaites H), and H alone at 2.55–2.65 with N_crit = 9. *Finding:* Re_x,tr = 6.3 × 10⁵, 1.46 × 10⁶, 2.27 × 10⁶, 5.8 × 10⁶, 1.2 × 10⁷, 3.1 × 10⁷ for m = −0.05, −0.02, 0, 0.05, 0.1, 0.2; with Blasius's θ and N_crit = 9, x_tr = 1.361, 0.865, 0.720, 0.514 m for H = 2.55, 2.591, 2.61, 2.65. Changing H from Blasius's 2.591 to Thwaites's 2.61 moves the transition point 17 % forward: the chart is steep in H, so H must come from a good boundary-layer method.

Raw data format: CSV tables `N, c_real, c_imag, err` (1); `alpha, Re_neutral, N` (2); `N, c_real, c_imag, matched` (3); `Tu, N_crit, Re_x_tr, x_tr_m` (4); `m_or_H, Re_theta0, Re_x_tr, x_tr_m` (5); plus `transition.png`.

## 5. Validate and improve

| quantity | published / analysis | code | difference | cause |
|---|---|---|---|---|
| c at α = 1, Re = 10⁴ | 0.23752649 + 0.00373967i (Orszag) | same, N ≥ 50 | 4 × 10⁻⁹ | published to 8 decimals |
| Re_c, α_c | 5772.22, 1.02056 (published) | 5772.22, 1.02054 (N = 60) | 3 × 10⁻⁷; 2 × 10⁻⁵ | flat minimum: Re barely changes with α near α_c |
| low-Re decay | αc_iRe constant | equal at Re = 1 and 0.5 to 10⁻³ | — | inertia negligible |
| flat-plate Re_x,tr | chart (3): 2.79 × 10⁶ | `flat_plate_transition` | — | consistency only |
| N_crit at Tu = 0.07 % | Mack: 9.0 | 9.00 | — | correlation |

The stability solver is validated against published data to the digits published. The e^N half is not validated here, only made self-consistent: its inputs are an approximate fit and a correlation. **The engineering answer: on the light aeroplane's wing at cruise, treated as a flat plate, Tollmien–Schlichting waves start to grow a few centimetres from the leading edge (41 mm on the e^N chart, about 28 mm at the Orr–Sommerfeld critical point of the Blasius profile) and the e⁹ method puts transition at about 0.87 m, 58 % of the chord, in quiet air; 0.3 % free-stream turbulence moves it to 28 %, and a mild adverse gradient (m = −0.05) to about 13 %.** The answer's uncertainty is dominated by N_crit, the pressure gradient and H, not by the eigenvalue solver.

What to improve: compute the growth rates for the actual profiles instead of using the chart (Extension 1), and use the wing's real pressure distribution with `P08.11`. What not to improve: N beyond 60 for the least-stable mode.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1) and of Rayleigh's criterion; `transition.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) stating the wing's predicted instability and transition points with their sensitivity to N_crit and to the pressure gradient, and the transition routes the prediction ignores.

### Adversarial self-checks

- **What if the speed doubles?** Re_x doubles at each station, so x₀ and x_tr halve: transition moves to 29 % of the chord. A faster aeroplane has proportionally less laminar flow at the same N_crit.
- **What if the flow is very viscous (Re → 0)?** Every mode decays at a rate ∝ 1/Re (the test checks αc_iRe is constant): there is no instability below Re_c, and a laminar layer stays laminar.
- **Which assumption fails first in the real world?** That disturbances enter only as small Tollmien–Schlichting waves. Insects on the leading edge or a propeller wake bypass the linear stage; for the light aeroplane, the wing root behind the propeller is probably turbulent from near the leading edge.
- **What measurement could falsify the model?** Infrared thermography of the wing in smooth air: a transition line near 0.9 m would support it; one at 0.2 m on a clean wing would point to bypass or to an adverse gradient the flat-plate model ignores.
- **Does the result survive a change in solver or grid?** The critical point is the same to 10⁻⁶ for N = 40–60. The transition point does not survive a 0.7 % change in H (17 %), which is the honest statement of its uncertainty.
- **Which decision changes if the uncertainty doubles?** A laminar-flow design that needs transition behind 50 % chord fails at N_crit = 8 (49 %); doubling the uncertainty in N_crit turns "probably laminar to mid-chord" into "do not count on it".

## Check yourself

1. What is a normal mode, and what do α, c_r and c_i mean? When does the mode grow?
2. What is the critical Reynolds number of plane Poiseuille flow, on what velocity and length is it based, and what does it mean?
3. Why is Chebyshev collocation well suited to the Orr–Sommerfeld equation?
4. What is N in the e^N method, and why is N ≈ 9 used for flight?
5. On the wing, where does the layer become unstable and where does it become turbulent, and why are the two so far apart?
6. How does free-stream turbulence change N_crit and the transition point?
7. Why does a favourable pressure gradient delay transition?
8. Name one route to transition that the e^N method of this lesson does not model.
