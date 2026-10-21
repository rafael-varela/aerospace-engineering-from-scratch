# P07.06  Power Cycles: Brayton

> A gas turbine is two isentropes and two isobars. Its efficiency is set by the pressure ratio, its work by how far apart the two isobars are at the hot end, and Carnot sets a ceiling it can approach only by doing no work at all.

| | |
|---|---|
| **Learning outcomes** | 1. Build the state table of the ideal Brayton cycle (isentropic compression, constant-pressure heat addition, isentropic expansion, constant-pressure heat rejection) for a constant-cp gas and for NASA variable-cp air. 2. Derive η = 1 − π^(−(γ−1)/γ), the net work w_net = cpT₁[τ(1 − 1/x) − (x − 1)] and the pressure ratio of maximum work (T₃/T₁)^(γ/(2(γ−1))), and show that maximum work occurs where T₂ = T₄ for any ideal gas. 3. Show that the net work is the area enclosed by the cycle in both the p–v and T–s planes, and compute it by a numerical contour integral that converges at second order. 4. Explain, with entropic mean temperatures, why the Brayton efficiency is below the Carnot efficiency between the same extreme temperatures. 5. Decide what the ideal cycle can and cannot say about a real engine. |
| **Prerequisites** | `P07.05` (isentropes with variable cp; the `Gas` model), `P07.04` (entropy, the Carnot efficiency, the Clausius integral), `P07.02` (steady-flow energy equation), `P03.04` (trapezoid rule and its error), `P03.14` (one-dimensional optimisation). Skill check: you can compute T₂ₛ for air at a given pressure ratio with constant and with variable cp. |
| **Engineering question** | An illustrative gas-turbine core takes sea-level air (288.15 K, 101 325 Pa) through a pressure ratio of 12 to a turbine inlet temperature of 1400 K (typical, illustrative): at best, with every component ideal, what thermal efficiency and specific work can it deliver, how close is that to Carnot, and what pressure ratio would give the most work? The answers: η = 0.508 and 415.7 kJ/kg for constant cp (0.483 and 448.6 kJ/kg with NASA air), against a Carnot bound of 0.794; maximum work at π = 15.9 (19.4 with NASA air), where T₂ = T₄ = √(T₁T₃) for constant cp. |
| **Assumptions and limits** | Air-standard cycle: the working fluid is air throughout, combustion is replaced by heat addition and the exhaust and intake by heat rejection to close the loop, so the fuel's mass and the products' properties are ignored. Every component is ideal: isentropic compressor and turbine, no pressure losses, kinetic energy changes neglected. Ideal gas with constant cp (γ = 1.4) or the NASA dry-air model. Real engines lose 10–20 points of efficiency to component losses (`P07.07`), and combustion itself destroys about a third of the fuel's exergy (`P07.08`); the ideal cycle is a bound and a guide to trends, not a performance prediction. |
| **Analytical work** | The four processes and the state table; first law for each component and for the cycle; η = 1 − π^(−m); w_net(π, τ) and its maximum; dw_net/d ln π = R(T₄ − T₂) for any ideal gas; the enclosed-area identity ∮p dv = ∮T ds = w_net; entropic mean temperatures and the Carnot bound; the zero-work pressure ratio. |
| **Computational deliverable** | `starter/main.py`: `specific_entropy`, `brayton_states`, `brayton_performance`, `efficiency_closed_form`, `net_work_closed_form`, `carnot_efficiency`, `optimum_pressure_ratio_closed_form`, `optimum_pressure_ratio` (numerical, `scipy.optimize.minimize_scalar`), `cycle_path`, `enclosed_area` (the NASA air model, `Gas`, `constant_cp_gas`, the temperature inversions, `isentropic_temperature`, `observed_order` and `main` are given). Tests (17) prove: the state table has two isentropes and two isobars and satisfies pv = RT; the first law closes for the cycle; efficiency and work match the closed forms; zero work at π = 1 and at π_max, where η reaches Carnot; the enclosed area in p–v and T–s equals w_net for both gas models and converges at order 2; isentropic legs add nothing to ∮T ds; η < η_Carnot; the numerical optimum matches the closed form with T₂ = √(T₁T₃); and T₂ = T₄ at maximum work for NASA air too. |
| **Practical activity** | **computational experiment**: a sweep of pressure ratio and turbine inlet temperature for efficiency and work; a convergence study of the contour integrals; constant against variable cp; the numerical optimum against the closed form. |
| **Equipment** | None. A computer with the course's Python environment. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `brayton_cycle.png` (p–v and T–s diagrams, work and efficiency against π), the experiment tables, the derivations of η, w_net, the optimum and the Carnot bound, and a half-page note answering the engineering question. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: η = 1 − π^(−m) and the optimum π derived; all public tests pass; the state table at π = 12 and T₃ = 1400 K reproduced by hand with constant cp; the note states the ideal efficiency and work, the Carnot gap and its cause, and what the ideal cycle omits. |
| **References** | `@mit-16-unified-thermodynamics-propulsion-notes` Ch. 3 "The First Law Applied to Engineering Cycles": 3.1 Some Properties of Engineering Cycles; Work and Efficiency; 3.3 The Carnot Cycle; 3.7 Brayton Cycle (3.7.1 Work and Efficiency; 3.7.2 Gas Turbine Technology and Thermodynamics), https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node28.html; and Ch. 6 "Applications of the Second Law": 6.3 Representation of Thermodynamic Processes in T-s coordinates; 6.4 Brayton Cycle in T-s Coordinates (6.4.1 Net work per unit mass flow in a Brayton cycle), https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node47.html (required-free; read online, link only: the cycle, its efficiency, the pressure ratio for maximum work and the T–s picture). `@mit-ocw-16-050-thermal-energy` Part 1A 1.A.3 "Brayton Cycles (or Joule Cycles): The Power Cycle for a Gas Turbine Jet Engine" and 1.A.4 "Gas Turbine Technology and Thermodynamics", https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/03_part1a_to_cms/; Part 1C 1.C.4 "Brayton Cycle in T-s coordinates"; Part 2A 2.A.3 "Brayton Cycle", https://ocw.mit.edu/courses/16-050-thermal-energy-fall-2002/resources/07_part2a/ (supplementary-free; the openly licensed earlier version of the same notes). `@openstax-university-physics-vol2` 4.2 "Heat Engines"; 4.5 "The Carnot Cycle" (supplementary-free; efficiency and the Carnot bound from the physics side). `@nasa-tp-2002-211556-glenn-coefficients` "Empirical Equations for Fitting Thermodynamic Functions", https://ntrs.nasa.gov/citations/20020085330 (required-free; the variable-cp air model). |
| **Extensions** | ★ Add a ram-compression inlet (flight Mach number M₀) and a nozzle to make the ideal turbojet; show that the ideal ramjet (π_c = 1) has η = 1 − 1/(1 + (γ − 1)M₀²/2). ★ Replace heat addition at constant pressure by heat addition at constant volume (the Humphrey cycle of pulse-detonation concepts) and compare the efficiency at the same compressor pressure ratio. ★ Put the T–s area integral on the variable-cp cycle with Simpson's rule and show the order rises to 4. |
| **Time** | unmeasured |

## 1. The engineering problem

Every jet engine, turboprop, turboshaft, APU and land-based gas turbine runs on the same cycle: compress air, heat it at roughly constant pressure by burning fuel in it, expand it through a turbine that drives the compressor, and either take the surplus as shaft power or let it accelerate the exhaust into a jet. The cycle is named after George Brayton (and in Europe after James Joule). Before any component is designed, the cycle has to answer two questions: how much of the heat released becomes work (the **thermal efficiency**), and how much work each kilogram of air delivers (the **specific work**, which sets the engine's size for a given power, because a turboprop's power is its specific work times its air mass flow).

The two answers pull in different directions. Raising the pressure ratio raises efficiency without limit (in the ideal cycle), but the work per kilogram rises, peaks and falls back to zero. Raising the turbine inlet temperature raises the work and the pressure ratio at which it peaks, but turbine materials and cooling set the limit. The ideal cycle makes these trade-offs exact, and it gives the bound against which `P07.07` measures a real engine with component losses and `P07.08` measures where the remaining work goes.

The working example is an illustrative core: sea-level air, pressure ratio 12, turbine inlet temperature 1400 K. Lessons 07 and 08 use the same numbers. Success criteria, written first: the state table has two isentropes (s₂ = s₁, s₄ = s₃) and two isobars; q_in − q_out = w_net to round-off; η and w_net match the closed forms to 10⁻¹⁰ for constant cp; the areas enclosed in the p–v and T–s planes equal w_net and converge at second order; η < 1 − T₁/T₃ for every π; and the numerical work maximum agrees with (T₃/T₁)^(γ/(2(γ−1))).

## 2. Understand and derive

### 2.1 The air-standard cycle

![Air enters the compressor C at station 1 and leaves at 2 into the combustor CC, where the fuel adds heat q_in; it enters the turbine T at 3 and leaves as exhaust at 4, rejecting q_out to the atmosphere. A shaft carries turbine work w_t to the compressor, which absorbs w_c.](figures/brayton-components.svg)
*Figure 1. The air-standard Brayton cycle as components and stations. The combustor is modelled as heat q_in added at constant pressure, and rejecting q_out to the atmosphere stands in for the exhaust and the fresh intake.*

The **air-standard** idealisation (assumptions named here and used throughout): the working fluid is air at every station; the combustor is a heat exchanger that adds q_in at constant pressure; the exhaust and the fresh intake are replaced by a fourth process that rejects q_out at constant pressure and returns the air to state 1; the compressor and turbine are adiabatic and reversible (isentropic); kinetic energy changes are neglected. The four states:

1. compressor inlet, (T₁, p₁);
2. after isentropic compression to p₂ = πp₁: s°(T₂) = s°(T₁) + R ln π (`P07.05`);
3. after constant-pressure heating to the turbine inlet temperature T₃;
4. after isentropic expansion to p₁: s°(T₄) = s°(T₃) − R ln π.

The steady-flow energy equation applied to each component (`P07.02`) gives, per kilogram of air,

w_c = h₂ − h₁,  q_in = h₃ − h₂,  w_t = h₃ − h₄,  q_out = h₄ − h₁,   (1)

and the cycle's first law: **w_net = w_t − w_c = q_in − q_out**. Thermal efficiency η = w_net/q_in; the **back-work ratio** w_c/w_t is the fraction of the turbine's work the compressor takes back.

### 2.2 Efficiency for constant cp

With constant cp the temperature ratios across both isentropes are the same, T₂/T₁ = T₃/T₄ = π^m with m = (γ − 1)/γ. Then

q_out/q_in = (T₄ − T₁)/(T₃ − T₂) = T₁(T₄/T₁ − 1)/(T₂(T₃/T₂ − 1)) = T₁/T₂,

because T₄/T₁ = T₃/T₂ (rearrange the isentropic ratios). So

**η = 1 − T₁/T₂ = 1 − π^(−(γ−1)/γ).**   (2)

The efficiency depends only on the pressure ratio, not on T₃. At π = 12, η = 0.5083. A monatomic gas (γ = 5/3, m = 0.4) does better at the same π: 0.6299.

**The state table** at π = 12, T₃ = 1400 K, γ = 1.4 (R = 287.05 J/(kg·K), cp = 1004.7 J/(kg·K)):

| state | T (K) | p (kPa) | v (m³/kg) |
|---|---|---|---|
| 1 | 288.15 | 101.33 | 0.81630 |
| 2 | 586.08 | 1215.90 | 0.13836 |
| 3 | 1400.00 | 1215.90 | 0.33051 |
| 4 | 688.32 | 101.33 | 1.94995 |

w_c = 299.32 kJ/kg, w_t = 714.99 kJ/kg, **w_net = 415.68 kJ/kg**, q_in = 817.71 kJ/kg, q_out = 402.03 kJ/kg, **η = 0.5083**, back-work ratio 0.4186. The compressor eats 42 % of what the turbine produces, which is why compressor and turbine efficiencies matter so much in `P07.07`: a small loss on each of two large numbers is a large loss on their difference.

### 2.3 Specific work and its maximum

Write τ = T₃/T₁ and x = π^m. Then T₂ = xT₁ and T₄ = T₃/x, and

**w_net = cp T₁ [τ(1 − 1/x) − (x − 1)].**   (3)

It is zero at x = 1 (π = 1: no compression, no work) and again at x = τ (π = τ^(1/m): the compressor alone heats the air to T₃ and no heat can be added). Between them it peaks where dw_net/dx = cpT₁(τ/x² − 1) = 0:

**x = √τ,  π_opt = (T₃/T₁)^(γ/(2(γ−1))),  T₂ = T₄ = √(T₁T₃).**   (4)

For T₃ = 1400 K, π_opt = 15.90, T₂ = T₄ = 635.15 K and w_net = 419.80 kJ/kg (η = 0.5463 there): 1 % more work than at π = 12. The work curve is flat near its maximum, so an engine designer can move the pressure ratio for efficiency at little cost in size. The zero-work pressure ratio is τ^3.5 = 252.8.

**A general result.** The optimum T₂ = T₄ does not need constant cp. Along the compressor isentrope, dh₂ = v₂ dp₂ = RT₂ d ln π (the Gibbs equation with ds = 0); along the turbine isentrope, the exit pressure is fixed at p₁ but the inlet pressure rises with π, and dh₄ = −RT₄ d ln π. With h₁ and h₃ fixed,

**dw_net/d ln π = R(T₄ − T₂)**,   (5)

for any ideal gas. The work rises while the turbine exit is hotter than the compressor exit and peaks when they are equal. With NASA air the optimum is π = 19.35, w_net = 460.32 kJ/kg, and T₂ = T₄ = 661.08 K; the tests check T₂ = T₄ at the numerical optimum for both gases. Equation (5) will return in `P07.07`: the condition T₄ > T₂ is exactly when a recuperator can help.

### 2.4 Work as an enclosed area

**In the p–v plane.** For a steady-flow device the shaft work out is w = −∫v dp (from the Gibbs equation with ds = 0, dh = v dp); the cycle's net work is therefore −∮v dp. Since ∮d(pv) = 0 around a closed path, −∮v dp = ∮p dv: **the net work is the area enclosed by the cycle in the p–v plane**, traversed clockwise. Only the isentropes contribute to ∮v dp (dp = 0 on the isobars), while all four legs contribute to ∮p dv; the two integrals are equal anyway.

**In the T–s plane.** Every process is internally reversible, so δq = T ds on each leg, and ∮T ds = q_in − q_out = w_net. The isentropes are vertical lines and contribute nothing; **the net work is the area between the two isobars**. (`P07.05` showed that T ds is not heat in an irreversible process; the T–s area is the net work only for internally reversible cycles.)

![The ideal Brayton cycle in two planes. In the T–s plane the isentropes 1–2 and 3–4 are vertical and the p₂ and p₁ isobars rise with entropy and diverge; heat q_in is added along 2–3 and q_out rejected along 4–1. In the p–v plane the isobars are horizontal and the isentropes curve; in both planes the area enclosed, traversed clockwise, is w_net.](figures/brayton-t-s-p-v.svg)
*Figure 2. The ideal cycle in the T–s and p–v planes, drawn to scale from the state table (π = 12, T₃ = 1400 K, constant cp). The isentropes are vertical in T–s; in both planes the area enclosed, traversed clockwise, is the net work.*

The tests check both identities numerically: `cycle_path` lays the closed loop out as 4n + 1 points (the isentropes in geometric pressure steps with T from the isentrope; the isobars in equal temperature steps), and `enclosed_area(x, y)` sums ½(yₖ + yₖ₊₁)(xₖ₊₁ − xₖ) around the loop. This trapezoid contour integral is exact for straight segments and has an O(h²) error for curved ones (`P03.04`), so its error must fall by 4 each time n doubles.

### 2.5 Why Brayton is below Carnot

Carnot's efficiency between the coldest and hottest temperatures in the cycle, T₁ and T₃, is 1 − T₁/T₃ = 0.7942 (`P07.04`). The Brayton cycle takes in heat over a range of temperatures, T₂ to T₃, and rejects it over T₄ to T₁. Define the **entropic mean temperature** of a heat transfer as T̄ = q/Δs: the constant temperature at which the same heat would carry the same entropy. The heat added and rejected carry the same Δs (the isentropes connect the two isobars), so

**η = 1 − q_out/q_in = 1 − T̄_out/T̄_in.**   (6)

For the example, Δs = 874.83 J/(kg·K), T̄_in = 934.71 K and T̄_out = 459.56 K: 1 − 459.56/934.71 = 0.5083. Because T̄_in < T₃ and T̄_out > T₁, η is strictly less than Carnot. The gap closes only as π → τ^(1/m), where T₂ → T₃ and T₄ → T₁: η → 1 − T₁/T₃, but the cycle shrinks to a line and w_net → 0. A heat engine that matches Carnot does nothing useful: the trade between efficiency and work is real, and it is why real engines are designed for neither the efficiency limit nor the work maximum alone.

### 2.6 Constant cp against NASA air

With NASA air, the same cycle has T₂ = 580.54 K (cp rises, so the isentrope runs cooler; `P07.05`) and T₄ = 750.95 K (γ is about 1.33 in the hot expansion, so the gas cools less per unit of pressure ratio). The work is larger (448.64 against 415.68 kJ/kg, +8 %) because the turbine expands a gas with higher cp, but the heat input is larger still (928.87 against 817.71 kJ/kg), and the efficiency falls from 0.5083 to 0.4830. The efficiency now depends on T₃ too: at π = 12 it is 0.4909, 0.4830 and 0.4760 for T₃ = 1100, 1400 and 1700 K. Equation (2) is a cold-air-standard result; the trend it states (efficiency rises with π) survives, its numbers do not, to about 2.5 points.

### 2.7 Where the model stops

- **Component losses.** A real compressor (η_c ≈ 0.85) and turbine (η_t ≈ 0.88) cut the example's efficiency from 0.48 to about 0.34 and move the optimum pressure ratio down (`P07.07`).
- **Pressure losses** in the combustor, ducts and exhaust; **cooling air** taken from the compressor to keep 1400 K turbine blades alive.
- **The working fluid.** Combustion adds 2–3 % fuel mass and changes the composition; the products' cp is higher than air's (`P07.08`).
- **Combustion.** Treating it as heat transfer hides the largest single irreversibility in the engine: burning fuel at 1400 K, not at the flame's own temperature or reversibly (`P07.08`).
- **The atmosphere is not a heat exchanger.** The "heat rejection" process is the exhaust mixing with the air; the ideal cycle's q_out is a bookkeeping device that closes the loop.

## 3. Build the computational model

Open `starter/main.py`. The NASA dry-air model, `Gas`, `constant_cp_gas`, `temperature_from_s0`, `temperature_from_h`, `isentropic_temperature` (all from `P07.05`), `observed_order` and `main` are given.

### Step 1: the state table and performance

`specific_entropy(gas, T, p)` is s°(T) − R ln(p/p°). `brayton_states(gas, T1, p1, pi, T3)` computes T₂ and T₄ with `isentropic_temperature` (π and 1/π) and returns arrays T, p, v = RT/p, h and s for the four states. `brayton_performance` evaluates (1), η, the back-work ratio and the Carnot efficiency from the state table's enthalpies.

### Step 2: closed forms and the optimum

`efficiency_closed_form` is (2); `net_work_closed_form` is (3); `carnot_efficiency(T_cold, T_hot)` is 1 − T_cold/T_hot; `optimum_pressure_ratio_closed_form` is (4). `optimum_pressure_ratio(gas, T1, p1, T3)` maximises `brayton_performance(...)["w_net"]` over ln π with `scipy.optimize.minimize_scalar(method="bounded")`, between π = 1.01 and the π at which T₂ₛ = T₃ (from s°(T₃) − s°(T₁) = R ln π_max).

### Step 3: the closed path and the contour integral

```python
u = np.linspace(0.0, 1.0, n + 1)[:-1]              # n points per leg, end point excluded
pp = p[0] * (p[1] / p[0]) ** u                     # isentrope 1-2: geometric in p
TT = temperature_from_s0(gas, s[0] + gas.R * np.log(pp / P_STD))
...
def enclosed_area(x, y):
    return float(np.sum(0.5 * (y[1:] + y[:-1]) * np.diff(x)))
```

The isobars are straight lines in T at constant p (curved in the T–s plane, since s rises as ln T). Close the loop by repeating the first point.

Running `main()` prints, for each gas model, the state table, the performance, the enclosed areas and the optimum:

```
  constant cp, gamma = 1.4
  w_c = 299.32, w_t = 714.99, w_net = 415.68, q_in = 817.71, q_out = 402.03 kJ/kg; eta = 0.5083, back-work ratio = 0.4186, Carnot = 0.7942
  enclosed area: p-v 415.6792 kJ/kg, T-s 415.6777 kJ/kg (w_net 415.6776); first law q_in - q_out - w_net = -1.16e-10 J/kg
  optimum pi for work = 15.900: w_net = 419.80 kJ/kg, eta = 0.5463, T2 = 635.15 K, T4 = 635.15 K

  dry air, NASA 9-coefficient
  w_c = 298.44, w_t = 747.08, w_net = 448.64, q_in = 928.87, q_out = 480.23 kJ/kg; eta = 0.4830, back-work ratio = 0.3995, Carnot = 0.7942
  enclosed area: p-v 448.6434 kJ/kg, T-s 448.6415 kJ/kg (w_net 448.6414); first law q_in - q_out - w_net = -5.82e-11 J/kg
  optimum pi for work = 19.352: w_net = 460.32 kJ/kg, eta = 0.5455, T2 = 661.08 K, T4 = 661.08 K
```

then the convergence study and the sweep of part 4, and writes `brayton_cycle.png`. At 5 kg/s of air (turboprop or APU scale, illustrative) the ideal NASA-air cycle at π = 12 would deliver 2.243 MW.

### Verify

`python -m pytest curriculum/phase-07/lesson-06-brayton-cycle` runs 17 tests:

1. **State table** (both gases): s₂ = s₁, s₄ = s₃, p on the two isobars, pv = RT, h = h(T), the temperature ordering.
2. **First law** (both gases): q_in − q_out = w_t − w_c = w_net to 10⁻¹²; η and the back-work ratio from their definitions.
3. **Closed forms**: η and w_net for constant cp against (2) and (3) at four π and three T₃; η = 1 − T₁/T₂; γ = 5/3 beats γ = 1.4.
4. **Limits**: w_net = 0 at π = 1 and at π_max, where η reaches Carnot.
5. **Enclosed areas** (both gases): ∮p dv and ∮T ds equal w_net to 10⁻⁵ with n = 400; the contour integral is signed (unit square ±1, circle 4π).
6. **Convergence**: observed order 1.8–2.2 for both planes; the isentropic legs are constant in s and the isobar constant in p.
7. **Carnot** (both gases): 0 < η < 1 − T₁/T₃ at four π and two T₃.
8. **Optimum**: the numerical maximum matches (4) to 10⁻⁵ with T₂ = √(T₁T₃); for NASA air too, T₂ = T₄ at the maximum and the work falls on both sides.
9. **main** writes the figure.

### Compare

Two independent routes to the net work: differences of enthalpies at four states, and the area of a 1601-point closed curve in either of two planes. They agree to 10⁻⁵ at n = 400, and the remaining difference falls by exactly 4 per halving of the step (part 4). The numerical optimiser, which knows nothing of (4), finds π = 15.900 against the closed form 15.900, and for NASA air finds a point where T₂ = T₄ to 10⁻⁵, the condition (5) predicts.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **Pressure ratio and turbine inlet temperature.** Sweep π = 4–40 at T₃ = 1100, 1400 and 1700 K with NASA air. *Predict* whether efficiency depends on T₃ and where the work peaks. *Finding:*

   | T₃ (K) | π | w_net (kJ/kg) | η (NASA) | η (2), γ = 1.4 | Carnot |
   |---|---|---|---|---|---|
   | 1100 | 4 | 228.74 | 0.3123 | 0.3270 | 0.7380 |
   | 1100 | 12 | 282.05 | 0.4909 | 0.5083 | 0.7380 |
   | 1100 | 40 | 213.83 | 0.6336 | 0.6514 | 0.7380 |
   | 1400 | 4 | 332.40 | 0.3059 | 0.3270 | 0.7942 |
   | 1400 | 12 | 448.64 | 0.4830 | 0.5083 | 0.7942 |
   | 1400 | 40 | 433.33 | 0.6264 | 0.6514 | 0.7942 |
   | 1700 | 4 | 436.14 | 0.3005 | 0.3270 | 0.8305 |
   | 1700 | 12 | 615.72 | 0.4760 | 0.5083 | 0.8305 |
   | 1700 | 40 | 654.55 | 0.6196 | 0.6514 | 0.8305 |

   The optimum pressure ratio for work is 10.43, 15.90 and 22.33 by (4), and 11.79, 19.35 and 29.19 with NASA air: hotter turbines want higher pressure ratios, and variable cp pushes the optimum up by 13–31 %. With variable cp the efficiency falls slightly as T₃ rises (the hotter gas has lower γ), the opposite of what a real engine shows: in a real engine hotter turbines raise efficiency because they reduce the relative weight of component losses (`P07.07`).

2. **Contour-integral convergence.** n = 25 to 400 points per process, NASA air, π = 12, T₃ = 1400 K. *Predict* the order. *Finding:* |∮p dv − w_net| = 500.2, 125.0, 31.25, 7.813 and 1.954 J/kg; |∮T ds − w_net| = 26.46, 6.615, 1.655, 0.4147 and 0.1046 J/kg; observed orders 2.00 and 2.00 (1.99 at the last step). The T–s integral is 19 times more accurate at every n because its integrand varies less along the isobars, and its isentropic legs contribute exactly zero.

3. **Constant against variable cp.** Compare the two state tables of part 3. *Predict* the sign of the change in T₂, T₄, w_net and η. *Finding:* T₂ falls 5.5 K, T₄ rises 62.6 K, w_net rises 8 %, η falls 2.5 points; the back-work ratio falls from 0.419 to 0.400.

4. **The optimum two ways.** Run `optimum_pressure_ratio` for both gases at T₃ = 1100, 1400 and 1700 K and compare with (4) and with the condition T₂ = T₄. *Predict* which agrees with which. *Finding:* the constant-cp optimum matches (4) to 10⁻⁵; the NASA-air optimum does not match (4) but does satisfy T₂ = T₄, as (5) says it must.

Raw data format: CSV tables `T3_K, pi, w_net_kJkg, eta_nasa, eta_cf, eta_carnot` (1); `n, err_pv_Jkg, err_Ts_Jkg` (2); `gas, state, T_K, p_kPa, v, h_kJkg, s` (3); `gas, T3_K, pi_opt_numeric, pi_opt_cf, T2_K, T4_K` (4); plus `brayton_cycle.png`.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| η, constant cp, π = 12 | (2): 0.5083 | 0.5083 | < 10⁻¹¹ | identity |
| w_net, constant cp | (3): 415.68 kJ/kg | 415.68 kJ/kg | < 10⁻¹⁰ | identity |
| first law of the cycle | q_in − q_out = w_net | residual 1.2 × 10⁻¹⁰ J/kg | round-off | — |
| ∮p dv (n = 400) | w_net = 448.6414 kJ/kg (NASA) | 448.6434 kJ/kg | 1.95 J/kg | trapezoid, order 2 |
| ∮T ds (n = 400) | w_net | 448.6415 kJ/kg | 0.10 J/kg | trapezoid, order 2 |
| optimum π, constant cp | (4): 15.900 | 15.900 (minimiser) | < 10⁻⁵ | optimiser tolerance |
| T₂ − T₄ at the NASA optimum | 0 by (5) | 661.08 − 661.08 K | < 10⁻⁵ relative | optimiser tolerance |
| η against Carnot | (6): 1 − 459.56/934.71 | 0.5083 < 0.7942 | — | entropic mean temperatures |

The engineering answer: **an ideal core at π = 12 and 1400 K converts 50.8 % of its heat into work with constant cp (48.3 % with NASA air) and delivers 416 (449) kJ/kg; Carnot between the same temperatures allows 79.4 %. The gap is not a loss in any component: it is the price of taking heat in between 586 and 1400 K and rejecting it between 688 and 288 K. A pressure ratio of 15.9 (19.4) would give 1 % (2.6 %) more work.**

What to improve: the cycle's weakest assumptions for a real engine are ideal components and the air standard. `P07.07` adds the component losses, which matter far more than variable cp; `P07.08` replaces the heat addition by real combustion.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (2)–(6); `brayton_cycle.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) answering the engineering question, with the Carnot gap explained.

### Adversarial self-checks

- **What if π → 1?** No compression: η → 0 and w_net → 0. The cycle collapses onto one isobar.
- **What if π → π_max = τ^(1/m)?** T₂ → T₃, η → Carnot, w_net → 0 (test 4). An engine at the Carnot limit has no output.
- **What if T₃ doubles?** For constant cp, η is unchanged (2) and the work more than doubles by (3); the optimum π rises as τ^1.75, by a factor of 3.4.
- **Which assumption fails first in the real world?** Isentropic components: with η_c = 0.85 and η_t = 0.88 the back-work ratio of 0.42 amplifies their losses (`P07.07`).
- **What measurement could falsify the model?** None directly: no real engine runs an ideal cycle. What can be measured is a real engine's state points; plotting them on the T–s diagram shows entropy rising across the compressor and turbine, which the ideal cycle forbids, and the area between the measured isobars then overstates the work.
- **Does the result survive a change of method?** Enthalpy differences and two contour integrals agree to the trapezoid error; the optimiser agrees with calculus.
- **Which decision changes if T₃ is 100 K lower than planned?** The work at π = 12 falls by 12.4 % (NASA air: 448.64 → 393.02 kJ/kg), so the engine needs 14 % more air for the same power and grows accordingly; the ideal efficiency barely moves (0.4830 → 0.4856).

## Check yourself

1. Name the four processes of the ideal Brayton cycle and the component that performs each.
2. Write the ideal Brayton efficiency for constant cp. Which design parameter does it depend on, and which does it not?
3. What is the back-work ratio, and why does a large back-work ratio make component efficiencies important?
4. At the pressure ratio of maximum work, what is the relation between the compressor exit and turbine exit temperatures?
5. Why is the net work of the ideal cycle equal to the area enclosed in the T–s diagram? Would this hold with an irreversible turbine?
6. Explain in one sentence, using mean temperatures, why the Brayton efficiency is below the Carnot efficiency between the same extreme temperatures.
7. Can the ideal Brayton cycle reach the Carnot efficiency? At what cost?
