# P08.07  Pipe Flow & Losses

> Every metre of line and every fitting costs pressure, and in a pressure-fed rocket the tank has to supply all of it. Friction goes as V², the velocity as D⁻²: the line diameter is the strongest lever in the budget.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the Darcy–Weisbach friction loss from the energy equation of `P08.04` and the laminar friction factor 64/Re from the Hagen–Poiseuille solution of `P08.06`, and explain what the Colebrook correlation adds for turbulent flow. 2. Implement and test Newton's method for Colebrook on x = 1/√f, prove that its iterates climb monotonically to the root after the first step, and check it against Haaland, the smooth and fully rough limits, Blasius and the laminar limit. 3. Build a tank-to-injector pressure budget with friction, minor losses, an injector orifice and the hydrostatic head under acceleration, for a line of several segments. 4. Size a propellant line for a loss limit, and decide which uncertainty (roughness, viscosity, fittings, injector) controls the tank pressure. |
| **Prerequisites** | `P08.04` (Bernoulli and the energy equation with losses), `P08.06` (Hagen–Poiseuille flow, f = 64/Re), `P08.05` (Reynolds number), `P03.02` (root finding: Newton and Brent), `P03.10` (Newton's convergence), `P00.05` (units). Skill check: you can write Newton's iteration for g(x) = 0 and say when it converges quadratically. |
| **Engineering question** | A small pressure-fed engine burns 1.6 kg/s of RP-1 at a chamber pressure of 2.0 MPa; the fuel goes through 2.5 m of 20 mm drawn tubing with a tank outlet, three bends, a ball valve, a filter and a check valve to an injector sized for a 0.4 MPa drop. What pressure must the fuel tank hold under 3 g of acceleration with 1 m of propellant above the injector, and what is the smallest line that keeps the line loss below 10 % of chamber pressure? The answers: 2.514 MPa, of which the line costs 137 kPa (6.9 % of p_c, with the fittings taking 95 kPa and the tube only 42 kPa); the smallest line is 18.3 mm. |
| **Assumptions and limits** | Incompressible, steady, isothermal liquid flow; fully developed friction along each segment (entrance effects lumped into K); the Colebrook correlation for turbulent flow (a fit to measurements, accurate to a few per cent), 64/Re for laminar flow, and an interpolation in the transition range 2300 < Re < 4000 where f is genuinely uncertain; loss coefficients and roughnesses are typical handbook values; the injector is an orifice with a constant discharge coefficient. No cavitation, no line dynamics (water hammer, pogo), no heat transfer; propellant properties are typical. |
| **Analytical work** | The energy equation with a head loss; Darcy–Weisbach; f = 64/Re and its equivalence to Hagen–Poiseuille; the structure of Colebrook (smooth and rough limits); Newton's iteration on 1/√f, its monotone convergence and quadratic rate; the D⁻⁵ law at fixed mass flow; the orifice equation; the pressure budget and the sign of the acceleration head. |
| **Computational deliverable** | `starter/main.py`: `reynolds`, `friction_laminar`, `friction_haaland`, `colebrook_residual`, `friction_colebrook`, `friction_blasius`, `friction_factor`, `darcy_weisbach_dp`, `minor_loss_dp`, `segment_loss`, `injector_dp`, `injector_area_for_dp`, `feed_line_budget`, `line_diameter_for_loss` (`fuel_line`, the typical property and loss-coefficient tables and `main` are given). Tests prove: f·Re = 64 from the integrated Poiseuille profile; Darcy–Weisbach equals Hagen–Poiseuille in laminar flow whatever the roughness; Newton solves Colebrook to a residual below 10⁻¹² in at most six iterations everywhere on the Moody chart and quadratically from a poor start; Haaland agrees within 2 %; the smooth limit is Prandtl's law and the fully rough limit is reached; Blasius within 3 %; f falls with Re and rises with roughness; the regimes join continuously; losses scale with L, V² and D^{−4.5…−5}; the injector inverts; the budget adds up, conserves mass through two segments, reduces to hydrostatics at zero flow and is unchanged in mm–g–ms units; the sizing root meets its target. |
| **Practical activity** | **computational experiment**: the Moody chart rebuilt from Colebrook and checked against Haaland across 4000 < Re < 10⁸; a line-diameter sweep for the fuel feed; the split of the budget between tube, fittings, injector and head; the effect of cleaning up the fittings. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. (Propellants and pressurised tanks are never handled in this course; the feed line is a calculation.) |
| **Evidence** | Completed `main.py`, passing tests, `feed_line.png`, the derivations of 64/Re, of Newton's convergence for Colebrook and of the D⁻⁵ law, the four experiment tables, and a half-page note giving the tank pressure, the chosen line size and the term that dominates the budget. |
| **Assessment** | `assessments/phase-08/grading-rubric.md`; competencies 1–4. Pass: Darcy–Weisbach and 64/Re are derived; Newton's convergence is argued from the shape of g; all public tests pass; the budget is presented term by term with units; the note names the dominant loss and justifies the line size with its margin. |
| **References** | `@leishman-intro-aerospace-flight-vehicles` chapter "Internal Flows": Entrance Length; Laminar Flow in Pipes; Turbulent Flows in Pipes; Colebrook-White Equation; Moody Chart; Pressure Drops Over Bends & Corners (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/internal-flows/) (required-free; link only). `@mit-ocw-1-060-engineering-mechanics-ii` L16 "Pipe Friction Losses, Moody Diagram, Minor Losses" (https://ocw.mit.edu/courses/1-060-engineering-mechanics-ii-spring-2006/resources/lecture16/) and L17 "Pipe Flow Analysis" (required-free: Darcy–Weisbach, the Moody diagram, minor losses); L15 "Velocity Distribution in Laminar Flows, Turbulence, The Reynolds Number" (supplementary-free). `@libretexts-ahmari-kabir-fluid-mechanics-lab-manual` 1.3 "Experiment #3: Energy Loss in Pipe Fittings" and 1.4 "Experiment #4: Energy Loss in Pipes" (supplementary-free: how K and f are measured). `@openstax-university-physics-vol1` 14.7 "Viscosity and Turbulence": Laminar Flow Confined to Tubes: Poiseuille's Law; Flow and Resistance as Causes of Pressure Drops (https://openstax.org/books/university-physics-volume-1/pages/14-7-viscosity-and-turbulence) (supplementary-free: a gentle first pass). |
| **Extensions** | ★ Add a cavitation check: the static pressure at every point of the line (including the vena contracta of the injector) against the propellant's vapour pressure. ★ Model a turbopump instead of a pressurised tank: the pump head required is the budget less the tank pressure, and the power is ṁ Δp/(ρη). ★ Add the line's inertance L/(A) and compute the time for the flow to start when a valve opens (a first look at feed-system dynamics). |
| **Time** | unmeasured |

## 1. The engineering problem

A pressure-fed rocket engine has no pump. Its propellant tanks are pressurised with gas, and the tank pressure alone must push the propellant through the feed lines, valves and filters, across the injector and into a combustion chamber that is itself at high pressure. Whatever the line costs in pressure, the tank must supply: and tank pressure is tank mass, because the wall thickness of a pressure vessel is proportional to its pressure. A feed line that wastes 1 MPa forces a heavier tank on the whole stage. Getting the budget right, and knowing which term in it is uncertain, is a central part of a feed-system design.

The engine here burns RP-1 (a rocket-grade kerosene; typical properties near 20 °C: ρ = 810 kg/m³, μ = 2.0 × 10⁻³ Pa·s) at 1.6 kg/s with a chamber pressure p_c = 2.0 MPa. Its fuel line is 2.5 m of drawn stainless tubing, 20 mm bore, roughness 1.5 µm (typical), with a sharp-edged tank outlet, three 90° bends, a ball valve, a filter, a check valve and the exit into the injector manifold. The injector's orifices are sized for a 0.4 MPa drop (20 % of p_c, a common choice that isolates the feed system from chamber pressure oscillations). In flight the vehicle accelerates at 3 g with 1 m of fuel above the injector.

Two questions: what tank pressure is needed, and how small can the line be before its loss exceeds 10 % of p_c? Success criteria, written first: the friction factor is computed by Newton's method to round-off and agrees with an explicit formula to its stated accuracy; the laminar limit reproduces 64/Re exactly; the budget closes term by term and conserves mass through a change of diameter; and the line size is found by a root-finder whose answer meets the limit to 10⁻⁸.

## 2. Understand and derive

### 2.1 Head loss in a pipe

The energy equation of `P08.04` between two sections of a pipe carrying an incompressible fluid, with no shaft work, is

p₁/ρ + V₁²/2 + gz₁ = p₂/ρ + V₂²/2 + gz₂ + g h_L,   (1)

where h_L, the **head loss**, is the mechanical energy per unit weight turned into heat by viscous dissipation. In a straight pipe of constant diameter, V₁ = V₂, and a momentum balance on the fluid in the pipe (`P08.03`) shows that the pressure drop is balanced by the wall shear stress: (p₁ − p₂)(πD²/4) = τ_w(πD)L. Dimensional analysis (`P08.05`) says that τ_w can depend only on ρ, V, D, μ and the wall roughness height ε, so τ_w/(ρV²) is a function of Re = ρVD/μ and ε/D. Define the **Darcy friction factor** f = 8τ_w/(ρV²); then

**Δp = f (L/D) ρV²/2**   (Darcy–Weisbach)   (2)

All the physics is now in f(Re, ε/D) (Leishman, "Internal Flows", Turbulent Flows in Pipes; MIT 1.060 L16).

### 2.2 Laminar flow: exact

For laminar flow, `P08.06` gave the exact Hagen–Poiseuille solution, τ_w = GR/2 and V = GR²/(8μ), so

**f = 64/Re**   (laminar, exact)   (3)

Substituting (3) into (2) gives Δp = 128μLQ/(πD⁴): the density cancels and so does the roughness. A laminar flow does not feel small wall roughness. Laminar flow in a pipe is typically observed up to Re ≈ 2300; above that, disturbances grow into turbulence, and by Re ≈ 4000 the flow is turbulent in practice. Between the two, f depends on how disturbed the inlet is, and no formula is reliable (Leishman, "Internal Flows", Laminar Flow in Pipes).

### 2.3 Turbulent flow: a correlation

Turbulent pipe flow has no exact solution. Its friction factor was measured by Nikuradse (sand-roughened pipes) and for commercial pipes, and the measurements are fitted by the **Colebrook equation**:

**1/√f = −2 log₁₀(ε/(3.7D) + 2.51/(Re √f)).**   (4)

Moody plotted (4) as the chart that bears his name. It is a correlation: accurate to a few per cent against the data it was fitted to, not a law. Its two limits are physical:

- **Smooth pipe** (ε = 0): 1/√f = 2 log₁₀(Re√f) − 0.8, Prandtl's universal law for smooth pipes, which comes from the logarithmic velocity profile near a wall (`P08.14`). f falls slowly as Re rises.
- **Fully rough** (Re → ∞): 1/√f = −2 log₁₀(ε/(3.7D)), von Kármán's rough-pipe law. The roughness elements stick out through the thin viscous layer at the wall, the drag on them is form drag, and f stops depending on Re.

Between them is the **transitionally rough** regime, where both terms in the logarithm matter. Because the two terms add inside one logarithm, f there can exceed both limiting values at the same Re (Problem 1). Typical roughnesses: drawn tubing 1.5 µm, commercial steel 45 µm, galvanised sheet metal 0.15 mm.

### 2.4 Solving Colebrook

(4) is implicit in f. Write x = 1/√f, a = ε/(3.7D) and b = 2.51/Re. Then

g(x) = x + 2 log₁₀(a + bx) = 0,   g′(x) = 1 + 2b/((a + bx) ln 10),   g″(x) = −2b²/((a + bx)² ln 10).

So g is **increasing and concave** wherever it is defined (a + bx > 0, which includes all x > 0), with exactly one root. Newton's iteration x_{k+1} = x_k − g(x_k)/g′(x_k) (`P03.02`) moves along the tangent, and the tangent of a concave function lies above it. So the tangent's zero has g ≤ 0: after one step every iterate is at or to the left of the root, and from there the iterates increase monotonically to it (g < 0 makes each step positive, and the next tangent zero is again at or left of the root). **From any start between 0 and the root, Newton converges monotonically**, and quadratically once close (`P03.10`). A start to the right of the root is thrown to the left by the first step; if it is far to the right, that step can leave the domain a + bx > 0 and the logarithm fails. That is one more reason to start from a good estimate. The fixed-point iteration x_{k+1} = −2 log₁₀(a + bx_k) also converges, but only linearly, at a rate of about 0.1 per step (Problem 5).

A good start costs nothing: **Haaland's explicit formula**

1/√f = −1.8 log₁₀[(ε/(3.7D))^{1.11} + 6.9/Re]   (5)

is within about 1.5 % of Colebrook over the whole turbulent chart. Starting Newton from (5) takes at most four iterations to round-off (the code measured 4 everywhere on 4000 < Re < 10⁸, 0 ≤ ε/D ≤ 0.05). For smooth pipes below Re ≈ 10⁵, **Blasius's** f = 0.3164 Re^{−1/4} is a third, simpler fit.

Every one of these is a function of Re and ε/D only. Check: the same line in millimetres, grams and milliseconds gives the same Re and the same f, with Δp in g/(mm·ms²) = 10⁶ Pa. The test `test_results_do_not_depend_on_the_units` makes that check.

### 2.5 Minor losses, the injector and the head

Fittings (bends, valves, tees, entrances, exits, filters) dissipate energy in separated flow and secondary flows. Their losses are written

**Δp = K ρV²/2,**   (6)

with a loss coefficient K measured for each fitting type (Leishman, "Internal Flows", Pressure Drops Over Bends & Corners; Ahmari & Kabir, Experiment #3). Typical values, used here and labelled typical because they vary with make and size: sharp-edged tank outlet 0.5, rounded outlet 0.05, 90° bend with r/D = 1.5 about 0.3, ball valve fully open 0.05, swing check valve 2.0, filter element 1.5, exit into a large manifold 1.0 (the jet's dynamic pressure is all lost). The fuel line's K adds to **5.95**.

The **injector** is a set of small orifices. An orifice turns pressure into velocity and loses it in the chamber: Δp_inj = (ṁ/(C_d A))²/(2ρ), with a discharge coefficient C_d (typically 0.6–0.8 for a plain orifice; 0.7 here). Sizing the area for 0.4 MPa at 1.6 kg/s gives A = 89.8 mm².

The propellant column of height H above the injector, in a vehicle accelerating axially at a (in g's), adds a hydrostatic pressure ρ a g₀ H that pushes the propellant toward the engine: it **reduces** the tank pressure needed. So the budget, with the tank's ullage gas at rest and the exit loss included in K, is

**p_tank = p_c + Δp_inj + Σ(f L/D + K) ρV²/2 − ρ a g₀ H.**   (7)

### 2.6 The budget, worked

At 1.6 kg/s in a 20 mm line: V = ṁ/(ρA) = **6.288 m/s**, Re = **50 930** (turbulent), ε/D = 7.5 × 10⁻⁵. Colebrook by Newton: **f = 0.02108** in 4 iterations (residual 9 × 10⁻¹⁶); Haaland 0.02084, Blasius 0.02106. The dynamic pressure ρV²/2 = 16.01 kPa. Then:

| term | value (kPa) |
|---|---|
| chamber pressure | 2000.0 |
| injector, C_d = 0.7, A = 89.79 mm² | 400.0 |
| tube friction, f L/D = 2.64 | 42.18 |
| fittings, K = 5.95 | 95.27 |
| head, 1 m at 3 g | −23.83 |
| **tank pressure** | **2513.6** |

The line costs 137.45 kPa, **6.87 %** of the chamber pressure, and the fittings cost more than twice the tube. The tube is short (L/D = 125) and the fittings many, so K = 5.95 beats fL/D = 2.64. The injector, chosen deliberately, is the largest loss of all.

**The D⁻⁵ law.** At fixed ṁ, V ∝ D⁻², so ρV²/2 ∝ D⁻⁴, and (2) gives Δp ∝ f D⁻⁵ (the minor losses go as D⁻⁴). With f constant (fully rough) the exponent is −5; with Blasius's f ∝ Re^{−1/4} ∝ D^{1/4} it is −4.75 (Problem 2). Halving the line multiplies its loss by about 27 to 32. Diameter is the lever.

### 2.7 Where the model stops

- **Transition.** Between Re = 2300 and 4000 the code interpolates; the real f may be anywhere between the laminar and turbulent values. A line designed to run there has an unpredictable loss (Problem 5).
- **Correlation accuracy.** Colebrook is a fit with scatter of a few per cent; the loss coefficients K are typical values and can be off by ±50 % for a particular fitting. Problem 6 shows that K and the injector's C_d dominate the uncertainty in the tank pressure.
- **Cavitation.** If the static pressure anywhere (at a valve seat, in the vena contracta of an orifice) falls to the vapour pressure, the liquid boils and the flow chokes. With a tank at 2.5 MPa that is far away here, but it governs pump inlets and cryogenic lines (Extension 1).
- **Unsteady flow.** Valve openings, closures and combustion oscillations drive pressure waves along the line (water hammer, pogo); the steady budget says nothing about them.
- **Temperature.** Viscosity varies by tens of per cent over a vehicle's temperature range; luckily f depends only weakly on Re.

## 3. Build the computational model

Open `starter/main.py`. `fuel_line`, the tables `K_TYPICAL`, the typical properties and `main` are given.

### Step 1: the friction factor

`reynolds` and `friction_laminar` are one line each. `friction_haaland` is (5); `colebrook_residual` is (4) written as a residual; `friction_blasius` is the fit.

`friction_colebrook(Re, rel_rough, x0=None, tol, maxiter)` is the Newton iteration of 2.4 on x = 1/√f, from Haaland unless x0 is given, stopping when the step is below tol × x and returning (f, the array of iterates). Returning the iterates lets the tests measure the convergence rate.

```python
a, b = rel_rough / 3.7, 2.51 / Re
x = 1.0 / np.sqrt(friction_haaland(Re, rel_rough)) if x0 is None else float(x0)
for _ in range(maxiter):
    g = x + 2.0 * np.log10(a + b * x)
    step = g / (1.0 + 2.0 * b / ((a + b * x) * np.log(10.0)))
    x -= step
    if abs(step) <= tol * abs(x):
        break
```

`friction_factor(Re, rel_rough)` chooses: 64/Re up to 2300, Colebrook from 4000, and a linear interpolation in Re between the two values in between (labelled in the docstring as an engineering patch).

### Step 2: losses and the budget

`darcy_weisbach_dp` and `minor_loss_dp` are (2) and (6). `segment_loss(mdot, rho, mu, D, L, roughness, K)` computes V from the mass flow, Re, f and both losses for one segment. `injector_dp` and `injector_area_for_dp` are the orifice equation and its inverse. `feed_line_budget` applies (7) over a list of segments (each with its own D, L, roughness and K, so a reducer or a change of tube is two segments carrying the same ṁ) and returns every term.

### Step 3: sizing

`line_diameter_for_loss` finds the D at which friction plus minor loss equals the limit, with `scipy.optimize.brentq` (`P03.02`): the loss falls monotonically with D, so a bracket of 1 mm to 0.5 m always contains one root.

Running `main()` prints

```
RP-1 (typical): rho = 810.0 kg/m^3, mu = 0.002 Pa s; mdot = 1.6 kg/s, chamber 2.0 MPa
line D = 20 mm, L = 2.5 m, roughness 1.5 um, total K = 5.95
  V = 6.288 m/s, Re = 50930, f (Colebrook) = 0.02108, Haaland 0.02084, Blasius 0.02106
  Newton from Haaland: 4 iterations, residual -8.9e-16
  dynamic pressure rho V^2/2 = 16.01 kPa
budget (kPa): chamber 2000.0 + injector 400.0 (Cd 0.7, A = 89.79 mm^2) + friction 42.18 + minor 95.27 - head 23.83 (1.0 m at 3.0 g) = tank 2513.6 kPa
  line loss 137.45 kPa = 6.87 % of chamber pressure
  smallest line for a loss of 10 % of p_c: D = 18.31 mm

diameter sweep: D (mm), V (m/s), Re, f, friction (kPa), minor (kPa), tank (MPa)
     10   25.150    101859  0.01870    1197.79    1524.27  5.0982
     12   17.466     84883  0.01921     494.54     735.08  3.6058
     16    9.824     63662  0.02019     123.32     232.59  2.7321
     20    6.288     50930  0.02108      42.18      95.27  2.5136
     25    4.024     40744  0.02207      14.47      39.02  2.4297
     32    2.456     31831  0.02329       4.45      14.54  2.3952
  log-log slope of line loss against D (16-32 mm): -4.23 (V^2 alone gives -4; f(Re) adds the rest)

Haaland vs Colebrook over 4e3 < Re < 1e8, 0 <= eps/D <= 0.05: max |difference| 1.37 %; Newton iterations from Haaland <= 4
```

and writes `feed_line.png`: the Moody chart rebuilt from 64/Re, Colebrook and Blasius, and the tank pressure against line diameter.

### Verify

Run `python -m pytest curriculum/phase-08/lesson-07-pipe-flow-losses`. The tests:

1. **Laminar**: f·Re = 64 from the Poiseuille profile integrated numerically and its wall stress; Darcy–Weisbach with `friction_factor` equals 128μLQ/(πD⁴) even with ε/D = 10⁻³.
2. **Colebrook by Newton**: residual below 10⁻¹² and at most six iterations over 25 Reynolds numbers and six roughnesses; e_{k+1}/e_k² bounded from x₀ = 1 (quadratic); Haaland within 2 %; the smooth limit is Prandtl's law to 10⁻¹²; the fully rough limit at Re = 10¹²; Blasius within 3 % for 5000 < Re < 10⁵; f decreasing with Re and increasing with roughness.
3. **Regimes**: 64/Re below 2300; continuity at 2300 and 4000; a ValueError for Re ≤ 0.
4. **Losses**: friction ∝ L; minor ∝ V²; turbulent friction grows by a factor between 3 and 4 when ṁ doubles (f falls with Re); the injector area inverts.
5. **Budget**: the terms add up to 10⁻¹⁴; ρVA = ṁ in each of two segments; the head is ρ a g₀ H; zero flow leaves only hydrostatics; mm–g–ms units give the same Re, f and V and the same Δp.
6. **Sizing**: the root meets the target to 10⁻⁸; a 5 % smaller line exceeds it; the friction exponent between −4.5 and −5 for a doubling of D.
7. **main** writes the figure.

### Compare

Haaland is the established explicit tool, and it agrees with our Newton solution to 1.37 % in the worst case (at Re ≈ 1.5 × 10⁵, ε/D = 10⁻⁴) and to 1.1 % for the fuel line. The difference is a property of Haaland's fit, not of the solver: Newton's residual is 10⁻¹⁵, so our f is Colebrook's. The second check is physical: the laminar limit is exact, and the smooth and fully rough limits are the log-law forms that Colebrook was built from. A Colebrook solver that passes all three limits and matches an independent fit within its stated accuracy is a trustworthy Colebrook solver; whether Colebrook is right for this tube is a separate question that only a flow test answers.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **The Moody chart and the explicit fit.** Colebrook by Newton against Haaland on 60 Reynolds numbers × 6 roughnesses. Predict the worst difference and the iteration count. *Finding:* worst 1.37 %, at Re = 1.47 × 10⁵ and ε/D = 10⁻⁴; never more than 4 Newton iterations from the Haaland start. Haaland alone is good enough for a budget; Newton costs four more evaluations and removes the doubt.
2. **Line diameter.** The fuel line from 10 to 32 mm (the table above). Predict the tank pressure at 10 mm. *Finding:* 5.10 MPa at 10 mm against 2.40 MPa at 32 mm; the line loss falls as D^{−4.23} between 16 and 32 mm (the minor losses go as D⁻⁴; the friction alone as D^{−4.79}); the 10 % limit is crossed at 18.3 mm. Below 15.6 mm the line loss exceeds the injector's 400 kPa and becomes the largest term after the chamber pressure.
3. **What the head buys.** Tank pressure at 1, 3 and 5 g. Predict the slope. *Finding:* 2529.5, 2513.6 and 2497.7 kPa: the 1 m column is worth 7.94 kPa per g. Small, but it is a free contribution that grows during the burn as acceleration rises (and shrinks as the column drains).
4. **Cleaning up the fittings.** Replace the sharp-edged outlet with a rounded one (K 0.5 → 0.05) and delete the check valve (K 2.0), so K falls from 5.95 to 3.50; separately, make the tube perfectly smooth. Predict which saves more. *Finding:* the fittings change saves 39.2 kPa (minor loss 95.3 → 56.0 kPa); a perfectly smooth tube saves 0.54 kPa. At ε/D = 7.5 × 10⁻⁵ and Re = 5 × 10⁴ the drawn tube is already hydraulically smooth. Effort spent on fittings pays; effort spent polishing a drawn tube does not.

Raw data format: CSV tables `Re, eps_D, f_colebrook, f_haaland, newton_iterations` (1); `D_mm, V, Re, f, dp_friction_kPa, dp_minor_kPa, p_tank_MPa` (2); `accel_g, p_tank_kPa` (3); `case, K, dp_friction_kPa, dp_minor_kPa, p_tank_kPa` (4); plus `feed_line.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| laminar f·Re | 64 (exact) | 64 (from the integrated profile) | < 10⁻¹⁰ | exact |
| Colebrook residual | 0 | 9 × 10⁻¹⁶ (4 iterations) | round-off | quadratic convergence |
| Haaland against Colebrook | ≈ 1.5 % (its stated accuracy) | ≤ 1.37 % over the chart; 1.1 % for the line | — | Haaland's fit |
| Blasius against Colebrook (Re = 50 930) | ≈ 3 % | 0.1 % | — | Blasius's fit, smooth tube |
| smooth limit | Prandtl's law | to 10⁻¹² | — | algebraic identity |
| budget sum | (7) | 2513.6 kPa, terms add to 10⁻¹⁴ | — | — |
| sizing root | loss = 200 kPa | 200 kPa to 10⁻⁸ at D = 18.31 mm | — | Brent tolerance |

Every check passes at the accuracy the method predicts, but these are checks of the arithmetic and of the correlation's internal consistency, not of the correlation against this line. The engineering answer: **the fuel tank must hold 2.51 MPa to deliver 1.6 kg/s at a 2.0 MPa chamber pressure; the 20 mm line costs 137 kPa (6.9 % of p_c), of which the fittings take 95 kPa; the smallest line that keeps the loss under 10 % is 18.3 mm.** The uncertainty sits in the fittings' K values and the injector's C_d, each worth about ±40–50 kPa (Problem 6), not in the friction factor.

What to improve: a cold-flow (water) test of the assembled line and injector, which measures K and C_d directly; then a rounded tank outlet and fewer fittings (Experiment 4). What not to improve: the friction-factor correlation, or the tube finish.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (3), of Newton's convergence for (4) and of the D⁻⁵ law; `feed_line.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the tank pressure term by term, the chosen line size and its margin, and the term whose uncertainty matters most.

### Adversarial self-checks

- **What if the flow doubles?** Every turbulent loss roughly quadruples (the injector exactly, the friction a little less because f falls with Re), so the line loss goes from 137 kPa to 528 kPa and the injector from 0.4 to 1.6 MPa: the feed system must be sized at the maximum flow, not the nominal.
- **What if the velocity goes to zero?** All losses vanish and p_tank = p_c − ρ a g₀ H, the hydrostatic limit the tests check.
- **Which assumption fails first in the real world?** The loss coefficients: a K of 2.0 for "a check valve" can be 1 or 4 for a particular valve. Then the transition patch, if a smaller line runs at low flow.
- **What measurement could falsify the model?** Pressure taps at the tank outlet, before and after the filter and at the injector manifold during a water-flow test: each Δp should scale as ṁ² (turbulent). A Δp that scales as ṁ means laminar flow (a blocked filter's pores), one that jumps means cavitation.
- **Does the result survive a change of method?** Newton, the fixed-point iteration, Haaland and Blasius agree to within their stated accuracies; the budget is unchanged by a change of units.
- **Which decision changes if the uncertainty doubles?** If K is uncertain by ±100 % instead of ±50 %, the line loss could reach 232 kPa, above the 10 % limit: the 20 mm line would then need either a measured K or the next tube size up.

## Check yourself

1. Define the Darcy friction factor and write the Darcy–Weisbach equation.
2. Why is f = 64/Re exact, while the Colebrook equation is a correlation?
3. What does ε/D measure, and when does f stop depending on Re?
4. Why does Newton's method on x = 1/√f converge in a few iterations from the Haaland start?
5. At fixed mass flow, by roughly what factor does the friction loss change if the line diameter is reduced by 20 %?
6. Write the tank-to-injector pressure budget and explain the sign of the hydrostatic term.
7. Why is the friction factor uncertain between Re = 2300 and 4000?
8. In the lesson's line the fittings cost more than the tube. Why, and what does that suggest about the design?
