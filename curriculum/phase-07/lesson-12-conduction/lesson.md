# P07.12  Heat Transfer: Conduction

> Every node is an energy balance. Get the balances right and the heat that leaves through the boundaries equals the heat made inside, to the last digit; get the closed forms right and you know how close the temperatures are.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the heat-conduction equation for a sheet with sources and face losses, and solve the 1D wall, the fin, the rectangle with one heated edge (separation of variables) and the semi-infinite solid (erfc) in closed form. 2. Build a vertex-centred finite-volume model of 2D steady conduction (Laplace/Poisson) and of 1D/2D transients with the θ method (explicit, Crank–Nicolson, implicit), with fixed-temperature, flux, convective and adiabatic boundaries. 3. Verify the solver against each closed form, show second-order convergence, close the global energy balance to round-off, and find the explicit stability limit. 4. Answer two aerospace questions with it: the hot spot of a conduction-cooled avionics card and the firing time of an uncooled rocket test-chamber wall. 5. Decide when the semi-infinite and lumped approximations may be used. |
| **Prerequisites** | `P03.11` (finite differences for PDEs, truncation error, stability), `P02.11` (the heat equation, separation of variables), `P01.13` and `P03.08` (linear systems; sparse direct solvers), `P03.02` (`brentq`), `P07.02` (energy balances), `P07.04` (entropy generation). Skill check: you can write the five-point Laplacian and explain why the forward-Euler heat equation needs α Δt/Δx² ≤ ½. |
| **Engineering question** | (1) A 160 × 100 mm conduction-cooled avionics card with four 35 µm copper planes is clamped along its long edges to a 70 °C chassis wall and carries a 10 W processor (30 × 30 mm) plus 5 W spread over the board (typical, illustrative). How hot is the board under the processor, and how thick an aluminium heat frame keeps it at or below 100 °C? (2) How long can a 25 mm copper heat-sink chamber wall fire (gas side h = 3000 W/(m²·K) at 3000 K, illustrative) before its hot face reaches 800 K? The answers: 129.2 °C without a frame, 0.33 mm of 6061 for 100 °C (1.5 mm gives 80.8 °C); 4.34 s, where the semi-infinite estimate (5.65 s) and the lumped one (5.87 s) are both 30 % too long, and an Inconel X-750 wall lasts 0.18 s. |
| **Assumptions and limits** | Fourier's law with constant, isotropic conductivity (properties at 20 °C from Lienhard Tables A.1 and A.2); the card is a thin sheet with uniform temperature through its thickness, faces adiabatic (sealed box), edges ideally clamped (no contact resistance), parallel layers adding their k·t; the chamber wall is a 1D slab, back face adiabatic, h and gas temperature constant, no radiation. Loads, h, gas temperature, wall temperature and the 800 K limit are typical, illustrative values. Fails: where k varies strongly with temperature (copper falls from 398 to 366 W/(m·K) by 600 °C), at contact interfaces (wedge-lock resistance can double the card's ΔT), at the component itself (junction-to-board resistance adds to the board temperature), and where the wall melts, ablates or radiates. |
| **Analytical work** | The heat equation for a sheet, div(k t ∇T) + q″ − h_f(T − T∞) = ρ c t ∂T/∂t; series resistances of a layered wall; the fin equation and its closed form; separation of variables for the rectangle; the erfc solutions for a surface step and for convection (Lienhard eq. 5.53) and the absorbed energy; the finite-volume row as an energy balance; global conservation from symmetry; truncation error and the explicit stability limit; Biot number and lumped capacity. |
| **Computational deliverable** | `starter/main.py`: closed forms (`wall_heat_flux`, `fin_closed_form`, `fin_heat_rate`, `rectangle_fourier`, `semi_infinite_step`, `semi_infinite_step_energy`, `semi_infinite_convection`), the finite-volume systems (`conduction_system_2d`, `conduction_system_1d`), `solve_steady`, `heat_balance`, `march`, `time_to_surface_limit`, `semi_infinite_time_to_limit` and `board_sheet_conductance`; `MATERIALS`, `diffusivity`, `observed_order`, `card_heat_load`, `solve_card` and `main` are given. Tests (17) prove the wall and Fourier closed forms' own boundary conditions, the semi-infinite limits, the rectangle against its Fourier series at second order, order 2.00 on a smooth problem, the exact parabola under uniform generation, the global energy balance to round-off with every boundary type, a layered wall to 10⁻¹², the 2D sheet against the fin closed form at second order, the card's linearity in 1/(k t), erfc agreement at second order for implicit and Crank–Nicolson, exact transient energy conservation and the absorbed-energy formula, the explicit stability limit, the steady limit of the transient, and the semi-infinite, thin-wall lumped and low-conductivity limits of the heat-sink wall. |
| **Practical activity** | **computational experiment**: a convergence study against the Fourier series, the fin and the erfc solutions; a heat-frame thickness sweep for the avionics card; a material and wall-thickness study for the heat-sink chamber, comparing finite-wall, semi-infinite and lumped answers; an explicit-scheme stability test. |
| **Equipment** | None (a computer with the course's Python environment). (Lesson `P07.15` measures a cooling curve at home and fits lumped and distributed models built on this lesson.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | `conduction.png`; the verification tables (Fourier, fin, wall, erfc, energy balance) with observed orders; the heat-frame sweep; the heat-sink table; a half-page note answering both engineering questions with the approximation each answer depends on. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the finite-volume row derived as an energy balance and the global balance shown to hold exactly; all public tests pass; second-order convergence demonstrated against a closed form; the heat-sink time given with the Biot and penetration-depth argument for why the semi-infinite and lumped answers fail. |
| **References** | `@lienhard-heat-transfer-textbook` Ch. 2: 2.1 The heat conduction equation ("Fourier's law"); 2.3 Thermal resistance and the electrical analogy; Ch. 4: 4.2 General solution of the heat conduction equation ("Separation of variables: A general solution of multidimensional problems", Example 4.1); 4.5 Fin design ("Analysis of a one-dimensional fin"; "Fin performance and design"); Ch. 5: 5.6 Transient heat conduction to a semi-infinite region (eq. 5.53 and Table 5.3); 5.7 Steady multidimensional heat conduction; Appendix A, Tables A.1 Properties of metallic solids and A.2 Properties of nonmetallic solids (required-free; https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf). `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 16 "Conductive Heat Transfer": 16.2–16.4; Ch. 18: 18.1 Temperature Distributions in the Presence of Heat Sources; 18.2 Heat Transfer From a Fin (required-free, read online; https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node115.html). `@mit-ocw-16-050-thermal-energy` Part 3 "Introduction to Engineering Heat Transfer": 2.0–2.2, 4.0, 5.0 (supplementary-free; https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/10_part3/). `@openstax-university-physics-vol2` 1.6 Mechanisms of Heat Transfer (Conduction) (supplementary-free). |
| **Extensions** | ★ Add a temperature-dependent conductivity (copper from Lienhard Table A.1, 0–600 °C) to the 1D transient by re-assembling K each step, and find how much it changes the heat-sink firing time. ★ Add a contact conductance per unit edge length at the wedge locks (the `('h', G, T)` edge type) and find the conductance at which the contact adds as much as the board. ★ Run the card's power-up transient with `march` and find the time to reach 95 % of the steady hot-spot rise; compare with ρ c t L²/(k t) for L the clamp spacing. |
| **Time** | unmeasured |

## 1. The engineering problem

Two conduction problems sit at opposite ends of an aerospace vehicle. In the avionics bay, a flight computer's circuit cards live in a sealed box with no air blown across them: each card is clamped along two edges (wedge locks) to a chassis wall that a cold plate or the airframe keeps cool, and every watt its processor dissipates must spread sideways through the board to those edges. The board is mostly glass-epoxy (FR4), a poor conductor, with thin copper planes inside it. Is that enough, or does the card need an aluminium heat frame bonded to it, and how thick?

In a rocket test cell, short firings often use an uncooled "heat-sink" chamber: a thick copper wall that simply soaks up heat until its hot face nears its limit, then the test stops. How long may the engine fire? A quick answer treats the wall as infinitely thick (the erfc solution) or as a single lump at one temperature; both are wrong in an instructive way.

The tools are the heat-conduction equation and its classic closed forms (Lienhard, chapters 2, 4 and 5), and a finite-volume solver built from `P03.11` that handles any shape, source and boundary. The closed forms check the solver; the solver answers the questions the closed forms cannot.

Success criteria, written first: the solver matches the Fourier series of a rectangle with one heated edge, the 1D layered wall (to round-off), the fin closed form and the erfc solutions; it converges at second order; heat out equals heat in to round-off in every steady case and the transient energy budget closes exactly; the explicit scheme is stable below Fourier number ½ and not above.

## 2. Understand and derive

### 2.1 The conduction equation for a sheet

Fourier's law says heat flows down the temperature gradient, q″ = −k ∇T (W/m²), with conductivity k (W/(m·K)). Take a thin sheet of thickness t whose temperature does not vary through its thickness. A patch of it of area dA stores ρ c t dA dT/dt; it gains heat by conduction across its edges, ∇·(k t ∇T) dA; it generates q″ dA (W/m² of sheet); and it may lose h_f (T − T∞) dA from its faces (h_f the total of both faces). The energy balance is

**ρ c t ∂T/∂t = ∇·(k t ∇T) + q″ − h_f (T − T∞)**.   (1)

The product k t (W/K) is the **sheet conductance**: layers side by side conduct in parallel, so their k t simply add. For the card: FR4 alone, 0.29 × 1.6 mm = 4.6 × 10⁻⁴ W/K; four 35 µm copper planes add 4 × 35 × 10⁻⁶ × 398 = 0.0557 W/K, 120 times more, for k t = 0.05614 W/K (an effective in-plane k of 35.1 W/(m·K)). The copper planes, 9 % of the thickness, carry 99 % of the heat. A 6061 frame adds 167 t_frame.

Steady with no face loss, (1) is Poisson's equation (Laplace's with no source). Without sources it is linear and every closed form below is a solution of it.

### 2.2 Closed forms that check the solver

**The layered wall.** In 1D steady conduction the flux is the same through every layer, so thermal resistances per unit area add in series:

q″ = (T_hot − T_cold)/(1/h_hot + Σ L_i/k_i + 1/h_cold).   (2)

For a cooled chamber wall (1 mm Inconel X-750 on 4 mm copper; gas h = 3000 W/(m²·K) at 3000 K; coolant h = 2 × 10⁴ at 300 K; illustrative), q″ = 5.630 MW/m², hot face 1123.4 K, interface 638.1 K, cold face 581.5 K. The gas-side film dominates; the thin Inconel layer, with 1/34 of copper's conductivity, takes 485 K.

**The fin.** A sheet with face losses and a fixed base temperature, uniform across its width, obeys θ″ = m² θ with θ = T − T∞ and m² = h_f/(k t). With an adiabatic tip,

θ = θ_b cosh(m(L − x))/cosh(mL),  Q = k t W m θ_b tanh(mL).   (3)

For a 1 mm 6061 sheet, 40 mm long and 20 mm wide, h_f = 50 W/(m²·K) on its faces in total, base 40 K above the air: mL = 0.6921, Q = 1.3855 W, fin efficiency 0.866.

**The rectangle with one heated edge.** On 0 < x < W, 0 < y < H with T = 0 on three edges and T = T₁ on the top, separation of variables (Lienhard §4.2; `P02.11`) gives

T = (4T₁/π) Σ_{n odd} sin(nπx/W) sinh(nπy/W)/(n sinh(nπH/W)).   (4)

The corners where T jumps from 0 to T₁ make the series converge slowly there, and any numerical method too; the checks use y ≤ 0.75 H. By symmetry the centre of a square is exactly T₁/4 (four rotated copies add to T₁ everywhere). With T = sin(πx/W) on the top edge only the first term survives, an exact smooth solution for convergence tests.

**The semi-infinite solid** (Lienhard §5.6). Initially at T_i, its surface held at T_s from t = 0:

(T − T_s)/(T_i − T_s) = erf(x/(2√(αt))),  absorbed energy 2k(T_s − T_i)√(t/(πα)),   (5)

with α = k/(ρc) the **diffusivity**. With convection to a fluid at T∞ instead (Lienhard eq. 5.53), Θ = (T − T∞)/(T_i − T∞) = erf(ζ/2) + exp(βζ + β²) erfc(ζ/2 + β), ζ = x/√(αt), β = h√(αt)/k. Heat has penetrated a distance of order √(αt): the solution describes a finite wall only while √(αt) is well below its thickness.

### 2.3 The finite-volume model

Put nodes on a uniform grid, Δx × Δy, including the boundaries. Give each node a **control volume**: the rectangle halfway to its neighbours (half a cell on an edge, a quarter at a corner). Between neighbours i and j heat flows at G_ij(T_i − T_j) with conductance G = k t (face length)/(spacing). The energy balance of node i is

C_i dT_i/dt = −Σ_j G_ij (T_i − T_j) + q″ V_i − h_f V_i (T_i − T∞) + (edge inflows),   (6)

with V_i the control-volume area and C_i = ρ c t V_i. Written for all nodes, **C dT/dt = f − K T**, with K symmetric. On a uniform interior grid the row is k t times the five-point Laplacian times V_i: the standard central difference, second-order accurate (truncation error (h²/12) × fourth derivatives). On an adiabatic edge the half-volume row is the second-order mirror condition, with no ghost nodes. Flux edges add q′ × (edge length) to f; convective edges add a conductance to the diagonal; fixed-temperature nodes have their rows replaced by T = T_b.

**Global conservation, exactly.** The conduction part of K has zero row sums and is symmetric, so Σ_i (K T)_i is just the face and edge-convection terms: the conduction terms cancel in pairs. Summing all rows, the heat drawn out through the fixed-temperature nodes (minus the residuals of their own, unenforced, energy balances) plus the face losses equals the generation plus the edge inflows, for any grid, to round-off. That is the global energy balance, and it is a property of the discretisation, not an accident of the solution. (A finite-difference scheme written point by point, with one-sided differences at the edges, would conserve energy only to its truncation error.)

**Steady**: solve K T = f on the free nodes with the fixed ones moved to the right-hand side (`scipy.sparse.linalg.spsolve`). **Transient**: the θ method,

(C/Δt + θK) Tⁿ⁺¹ = (C/Δt − (1 − θ)K) Tⁿ + f,   (7)

θ = 0 explicit (forward Euler), ½ Crank–Nicolson, 1 implicit Euler. The explicit scheme is stable only if Δt ≤ C_i/K_ii for every node: in 1D, **Fo = αΔt/Δx² ≤ ½** (`P03.11`). The implicit scheme is unconditionally stable and never oscillates; Crank–Nicolson is second order in time but does not damp the shortest wavelengths, so after a sudden change it can ring when Fo is large. Summing (7) over the free nodes gives the transient energy budget, again exact: stored energy change = Δt × (inflows at T^θ).

### 2.4 The two questions, by hand first

**The card.** If all 15 W were spread evenly, the card would be a 1D strip between the clamps: T_max − T_wall = q″H²/(8 k t) = (15/0.016)(0.1²)/(8 × 0.05614) = **20.9 K**. Concentrating 10 W on a 30 mm square at the centre makes it far worse; that needs the 2D model, which gives **59 K** (129.2 °C). Because (1) is linear in T with these boundaries, the rise above the wall scales as 1/(k t): doubling k t halves it (tested).

**The heat-sink wall.** Two quick estimates. Lumped (the whole wall at one temperature, valid if the Biot number Bi = hL/k ≪ 1): t = ρcL/h ln((T_g − T_i)/(T_g − T_lim)) = 5.87 s for 25 mm of copper. Semi-infinite (5.53), valid while √(αt) ≪ L: 5.65 s. For copper Bi = 0.188, not small, so the hot face runs ahead of the mean (the lumped answer is long), and √(αt) at 5 s is 24 mm, the whole wall (the semi-infinite answer is long, because the adiabatic back face reflects heat back). The finite-volume wall says **4.34 s**. For Inconel X-750 (Bi = 6.5) the hot face reaches 800 K in 0.18 s while the heat has penetrated only 0.75 mm: the semi-infinite answer (0.180 s) is right, and the lumped one (6.42 s) is absurd. Low-conductivity metals make poor heat sinks because their **effusivity** √(kρc) is low: 6605 J/(m²·K·s^½) against copper's 36 993.

### 2.5 What generates entropy, and where the model stops

Conduction is heat transfer across a finite temperature difference, so it always generates entropy: per unit volume k|∇T|²/T². In the card, 15 W falling from the processor footprint to a 343 K wall generates about 15 × (1/343.15 − 1/402.35) ≈ 6.4 mW/K at least; the steeper the gradients, the more.

The models stop where k(T) varies (copper loses 8 % of its k by 600 °C), at interfaces (wedge-lock contact resistance, solder and thermal pads under the processor), at the component (the junction runs hotter than the board under it), and in the chamber where radiation, melting or chemical attack set in.

## 3. Build the computational model

Open `starter/main.py`. `MATERIALS` (Lienhard Tables A.1 and A.2 at 20 °C), `diffusivity`, `observed_order`, `card_heat_load`, `solve_card` and `main` are given.

### Step 1: the closed forms

`wall_heat_flux` is (2), returning the interface temperatures by walking down the resistances. `fin_closed_form`, `fin_heat_rate` are (3); `rectangle_fourier` is (4), with the ratio of sinh written as exponentials so that n up to 800 does not overflow; `semi_infinite_step`, `semi_infinite_step_energy` and `semi_infinite_convection` are (5) and eq. 5.53 (use `scipy.special.erfcx(u) = exp(u²) erfc(u)` so exp(βζ + β²) erfc(ζ/2 + β) = exp(−ζ²/4) erfcx(ζ/2 + β) never overflows).

### Step 2: the finite-volume systems

`conduction_system_2d(W, H, nx, ny, kt, q_gen, bcs, h_face, T_face, rho_c_t)` builds K, f, C (6) by looping over x-faces and y-faces and calling a `link(a, b, G)` helper that adds G to both diagonals and −G to both off-diagonals, then the face term, then flux and convective edges, then the fixed-temperature edges (which win at corners). It records the per-node generation and edge inflows so `heat_balance` can audit them. `conduction_system_1d` is the same in 1D, with k per interval (a layered wall) and ρc per node.

```python
for j in range(ny):                                  # x-direction faces
    for i in range(nx - 1):
        G = float(kfun(x[i] + dx / 2, y[j])) * wy[j] / dx
        link(idx(i, j), idx(i + 1, j), G)
```

### Step 3: solving, auditing, marching

`solve_steady` partitions K into free and fixed blocks. `heat_balance` returns generation, edge inflows, face losses, the heat out through fixed nodes (minus the sum of their residuals) and the imbalance. `march` is (7) with one sparse LU factorisation (`splu`) reused for every step, the energy budget accumulated at T^θ, and an optional `stop(t, T)` callback. `time_to_surface_limit` marches a convective slab with Crank–Nicolson until the hot face reaches the limit and interpolates within the last step; `semi_infinite_time_to_limit` solves eq. 5.53 at x = 0 for the same time with `brentq`. `board_sheet_conductance` adds the k t of copper planes, FR4 and an optional frame.

Running `main()` prints

```
1D cooled wall (1 mm Inconel X-750 on 4 mm copper; gas h = 3000, 3000 K; coolant h = 20000, 300 K; illustrative):
  closed form q'' = 5.6298 MW/m^2, hot face 1123.40 K, interface 638.07 K, cold face 581.49 K
  finite volumes: hot face 1123.40 K, interface 638.07 K, cold face 581.49 K; max |diff| 4.9e-12 K
fin: ... m L = 0.6921
  max |FV - closed form| for n = 11, 21, 41, 81: 2.65e-03, 6.63e-04, 1.66e-04, 4.14e-05 K; orders [2. 2. 2.]
  heat rate: closed form 1.3855 W, FV 1.3855 W (energy imbalance 3.8e-12 W); efficiency 0.866
unit square, T = 1 on the top edge, 0 on the others: FV against the Fourier series (y <= 0.75)
  n =  21: max error 1.681e-03; centre T = 0.250000 (series 0.250000)
  n =  41: max error 4.367e-04; ...
  n =  81: max error 1.117e-04; ...
  smooth edge sin(pi x / W): errors 4.19e-04, 1.05e-04, 2.63e-05, 6.58e-06; observed orders [1.996 1.999 1.999]
```

followed by the card, transient and heat-sink tables of part 4, and writes `conduction.png`: the card's temperature map with no heat frame, and the copper heat-sink wall's temperature profiles at 0.5, 1, 2 and 4 s against the semi-infinite solution.

### Verify

Run `python -m pytest curriculum/phase-07/lesson-12-conduction`. The tests:

1. **Closed forms**: the wall's series resistances and equal flux through each; the Fourier series' zero edges and centre value ¼; erf and eq. 5.53 at their limits (h → ∞ gives the fixed-surface solution, h → 0 leaves the solid alone).
2. **Steady finite volumes**: the square against the Fourier series below 10⁻³ with an error ratio of 4 per halving; order 2.0 ± 0.1 on the smooth problem; the exact parabola under uniform generation (quadratic, so no truncation error); the global balance to 10⁻⁹ on the card and with flux, convection and face losses all present; the layered wall to 10⁻¹²; the 2D sheet against the fin at order > 1.9, with the base supplying what the faces lose; heat-frame and 1/(k t) scaling of the card.
3. **Transients**: erfc agreement below 0.25 K with order > 1.8 for implicit and Crank–Nicolson (Δt ∝ Δx²); stored energy equal to heat in to 10⁻¹⁰, and to the absorbed-energy formula within 0.2 %; eq. 5.53 within 0.3 K; explicit stable at Fo = 0.45 and blowing up at 0.55; implicit stable at Fo = 5; the 2D transient settling on the steady solution.
4. **Heat-sink limits**: a 100 mm wall matches the semi-infinite time within 1 %; a 10 mm wall is quicker; a 2 mm wall (Bi = 0.015) is within 3 % of the lumped time and slightly quicker; Inconel is more than ten times quicker than copper.
5. **main** writes the figure.

### Compare

Lienhard's book (§5.7) solves 2D steady problems by flux plots and shape factors, not numerically; MIT 16.Unified Ch. 18 gives the fin. The closed forms of 2.2 are the comparison, and every one is reproduced at the order the method predicts. A commercial thermal code would give the same card temperatures for the same model; its value is in what this model leaves out (contacts, the component stack), not in the numerics.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Convergence.** Fourier rectangle (n = 21, 41, 81: max error 1.681 × 10⁻³, 4.367 × 10⁻⁴, 1.117 × 10⁻⁴ for y ≤ 0.75), smooth edge (orders 1.996, 1.999, 1.999), fin (orders 2.00), erfc with Δt ∝ Δx² (implicit 3.511, 0.881, 0.221, 0.0553 K for n = 26–201, orders 2.00, 1.99, 2.00; Crank–Nicolson 0.636, 0.163, 0.0410, 0.0103 K, orders 1.96, 1.99, 2.00). Predict the order. *Finding:* 2 everywhere. The absorbed energy converges much faster with the implicit scheme (relative error 1.6 × 10⁻⁸ at n = 201) than with Crank–Nicolson (6.8 × 10⁻⁵): CN's undamped short-wavelength modes, excited by the sudden surface step, linger. The convective erfc case (h = 3000, n = 101, Δt = 1 ms, t = 1 s) agrees to 0.0815 K and its stored energy equals the heat in, 7.631421 × 10⁶ J/m², to all printed digits.
2. **Stability.** Explicit scheme, 400 steps: at Fo = 0.45 the largest temperature is 1000 K (the surface value; bounded); at Fo = 0.55 it is 1.085 × 10³¹ K. Predict before running. *Finding:* the ½ limit is sharp.
3. **The card.** Grid n = 41, 81, 161: hot spot 128.91, 129.20, 129.11 °C, with 15.000000 W out through the clamps of 15.000000 W generated (imbalance at most 1.8 × 10⁻¹¹ W). Heat frame 0, 0.5, 1.0, 1.5, 2.0, 3.0 mm of 6061: 129.20, 93.80, 84.89, 80.84, 78.52, 75.97 °C. Predict the thickness for 100 °C. *Finding:* **0.327 mm**; diminishing returns beyond about 1.5 mm, because the rise is proportional to 1/(k t) and k t is already dominated by the frame. A 1 mm frame halves the rise; the last 10 K would need a different path (a thermal pad to the cover, a second clamp edge).
4. **The heat-sink wall.** Gas h = 3000 W/(m²·K) at 3000 K, initial 300 K, limit 800 K, 25 mm walls:

| material | finite wall (s) | semi-infinite (s) | lumped (s) | √(αt) (mm) | Bi |
|---|---|---|---|---|---|
| copper | 4.341 | 5.654 | 5.87 | 22.42 | 0.188 |
| Inconel X-750 | 0.182 | 0.180 | 6.42 | 0.75 | 6.47 |
| stainless 304 | 0.184 | 0.182 | 5.46 | 0.89 | 5.44 |

Copper thickness 5, 10, 20, 25, 40, 60 mm: 1.116, 2.112, 3.729, 4.341, 5.393, 5.646 s. Predict which estimate fails for which wall. *Finding:* thin walls behave as lumps, thick walls as semi-infinite bodies, and a 25 mm copper wall is neither: past about 40 mm, extra copper buys little time, because heat never reaches it in time. Every run's stored energy equals its heat in to 10⁻⁵ J/m² or better.

Raw data format: CSV tables `n, err_fourier, err_smooth, err_fin` (1); `Fo, max_T_after_400_steps` (2); `n, T_hot_C, Q_out_W, imbalance_W` and `t_frame_mm, kt_W_per_K, T_hot_C` (3); `material, thickness_mm, t_fd_s, t_semi_inf_s, t_lumped_s` (4); plus `conduction.png`.

## 5. Validate and improve

| quantity | closed form | finite volumes | difference | cause |
|---|---|---|---|---|
| layered wall, hot face | 1123.40 K | 1123.40 K | 4.9 × 10⁻¹² K | exact (piecewise linear) |
| fin heat rate | 1.3855 W | 1.3855 W | order 2 in T | truncation |
| square, centre | 0.250000 | 0.250000 | — | symmetry |
| smooth rectangle | sinh solution | error 6.6 × 10⁻⁶ at n = 81 | order 2.00 | truncation |
| erfc at 0.5 s | (5) | 0.055 K at n = 201 | order 2.00 | truncation |
| card energy balance | 15 W in | 15.000000 W out | 10⁻¹¹ W | round-off |
| copper heat-sink time | semi-infinite 5.65 s | 4.34 s | −23 % | finite wall, not numerics |

The engineering answers: **the card's board reaches 129 °C under the 10 W processor with copper planes alone; 0.33 mm of 6061 brings it to 100 °C, and a 1.5 mm frame to 81 °C. The 25 mm copper heat-sink wall reaches 800 K after 4.3 s: both quick estimates promise 5.7–5.9 s, about 30 % too long; an Inconel or stainless wall of any thickness lasts under 0.2 s.**

What to improve: the contact resistance at the clamps and under the processor, and copper's k(T) (Extension 1) for the wall; the gas-side h, which sets the firing time almost in inverse proportion (Problem 6). What not to improve: the grid; n = 81 is converged to 0.1 K.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (6) and of the global balance; `conduction.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering both engineering questions and saying which approximation each depends on.

### Adversarial self-checks

- **What if the processor power doubles?** The hot spot goes from 129.2 to 181.4 °C (the processor's share of the rise doubles; the 5 W spread does not); doubling both loads doubles the whole rise, to 188.4 °C. A frame becomes mandatory.
- **What if Δt goes beyond Fo = ½?** Explicit marching explodes (10³¹ K in 400 steps); implicit stays bounded at Fo = 5; Crank–Nicolson stays bounded but can ring after a sudden step.
- **Which assumption fails first?** For the card, ideal clamps and no component stack; for the wall, constant h and constant k.
- **What measurement could falsify the model?** Thermocouples on the card at the processor, mid-span and at the clamp: the model predicts the shape (parabola-like across the clamps, a peak under the chip); a large step at the clamp would expose contact resistance. For the wall, a back-face thermocouple: semi-infinite predicts no rise for the first second, the finite wall shows it starting earlier.
- **Does the result survive a change of grid or method?** The card hot spot changes by 0.3 K from n = 41 to 161; the heat-sink time is converged in both Δx and Δt; three independent closed forms are reproduced.
- **Which decision changes if the uncertainty doubles?** If h is ±50 % instead of ±25 %, the copper firing time spans 2.35–10.3 s instead of a narrow band: the test plan would need a margin or a calibration firing.

## Check yourself

1. State Fourier's law and the units of k, α and k t.
2. Why do layers in series add their resistances L/k, and layers side by side add their k t?
3. What makes the finite-volume global energy balance close to round-off rather than to truncation error?
4. Why does the Fourier-series check exclude the corners of the heated edge?
5. What is the explicit stability limit in 1D, and why does the implicit scheme have none?
6. What do the Biot number and √(αt)/L each tell you about a transient wall?
7. Why is copper a good heat-sink wall and Inconel a poor one, even though Inconel has the higher melting point?
8. A grid refinement study shows the error falling by 4 each time h halves. What is the observed order, and what would an order of 1 suggest?
