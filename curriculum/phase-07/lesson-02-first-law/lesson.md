# P07.02  The First Law

> Energy is never lost, only moved across a boundary or carried across it by mass. Write the balance for the box you drew, and the tank's temperature, pressure and stranded gas all follow from one ODE.

| | |
|---|---|
| **Learning outcomes** | 1. Write the first law for a closed system (ΔU = Q − W) and for a control volume with mass flow (dE/dt = Q̇ − Ẇ + Σṁ_in h_in − Σṁ_out h_out), and say which form a given aerospace problem needs. 2. Compute closed-system boundary work for isothermal, adiabatic and polytropic processes of an ideal gas. 3. Derive the choked orifice mass flow from the isentropic nozzle relations and the closed-form adiabatic and isothermal blowdown of a tank. 4. Implement and test a tank-blowdown model whose energy balance closes to round-off, and check it against p ∝ ρ^γ, the closed forms and the isothermal limit. 5. Decide how much pressurant a blowdown strands in the tank and how wall heat transfer changes it. |
| **Prerequisites** | `P02.01` (first-order ODEs: separation of variables, exponential decay), `P07.01` (ideal gas, cv, cp, γ, h(T)), `P03.05`/`P03.06` (ODE integrators, `scipy.integrate.solve_ivp`, events), `P04.02` (work and energy). Skill check: you can solve dy/dt = −k yⁿ by separation of variables and integrate an ODE system with `solve_ivp`. |
| **Engineering question** | A 50 L nitrogen pressurant tank (a composite-overwrapped pressure vessel, COPV) at 30 MPa and 300 K feeds a propellant-tank ullage held at 0.5 MPa through a 2 mm orifice with Cd = 0.8 (typical, illustrative values). How long until the tank falls to 3 MPa, the lowest pressure the regulator downstream accepts; how cold does the gas get; and how much nitrogen is stranded? The answers: adiabatic, 189.6 s, 155.4 K and 3.25 kg (19.3 % of the 16.85 kg load); isothermal, 224.2 s, 300 K and 1.68 kg (10.0 %). A real tank lies between, set by the wall's heat-transfer conductance hA: at hA = 10 W/K, 205.7 s, 183.2 K and 2.76 kg. |
| **Assumptions and limits** | Ideal gas with constant cp and cv (γ = 1.4 for N₂, 5/3 for He); a uniform state inside the tank (well mixed, negligible kinetic energy); an isentropic orifice with a discharge coefficient, choked while the tank pressure exceeds 1.893 times the back pressure (N₂); constant back pressure; constant wall temperature (the wall's own heat capacity is not modelled). Fails: the ideal gas at tens of MPa (lesson 03 shows how to measure Z); the uniform state in a fast blowdown (stratification); the constant wall temperature when the liner itself cools; and near the end of a long adiabatic blowdown, where the model's gas temperature reaches nitrogen's condensation temperature. |
| **Analytical work** | The closed-system first law and boundary work; the control-volume energy balance specialised to a rigid tank with one exit; the choked mass-flow function from the isentropic nozzle relations; the closed-form adiabatic blowdown p/p₀ = (1 + (γ − 1)t/(2τ))^(−2γ/(γ−1)) and the isothermal exponential; the remaining-mass fractions; the recovery pressure; the thermal time constant m cv/(hA). |
| **Computational deliverable** | `starter/main.py`: `critical_pressure_ratio`, `choked_flow_function`, `mass_flow` (choked and subsonic), `blowdown_rhs` (the ODE right-hand side for adiabatic, isothermal and wall heat transfer), `energy_residual`, `blowdown_time_constant`, the closed forms (`adiabatic_choked_closed_form`, `isothermal_choked_closed_form`, `adiabatic_time_to_pressure`), `remaining_mass_fraction`, `polytropic_work`, `recovery_pressure`; `simulate_blowdown` (the `solve_ivp` driver) and `main` are given. Tests (17) prove the critical ratios, continuity of the mass flow at choking and its scaling, the energy balance closing below 10⁻¹² in all three heat modes, the isothermal heat R T₀ Δm, p ∝ ρ^γ in the adiabatic tank, Δu = ∫p dv for the gas left behind, agreement with both closed forms to 10⁻⁷ for N₂ and He, the stop event, wall heat transfer bounded by the two limits, the stranded mass and recovery pressure, and the approach to the back pressure. |
| **Practical activity** | **computational experiment**: a sweep of wall heat-transfer conductance between the adiabatic and isothermal limits, nitrogen against helium, orifice size, and a check of where the model's own assumptions break. No pressure vessel is built, filled or vented at home in this course. |
| **Equipment** | None (a computer with the course's Python environment). |
| **Safety** | None beyond normal computer use. This lesson is computational only: **pressure vessels are never built, filled or vented at home in this course.** A 30 MPa tank stores enough energy to kill; tank work belongs in a qualified test facility with remote operation. |
| **Evidence** | `tank_blowdown.png`; the printed verification tables; the experiment tables; a half-page note answering the engineering question with the adiabatic, isothermal and finite-hA results and the assumption most likely to be wrong. |
| **Assessment** | `assessments/phase-07/grading-rubric.md`; competencies 1–4. Pass: the control-volume energy balance and the closed-form adiabatic blowdown derived; all public tests pass; the closed-form time to 3 MPa reproduced by hand to three figures; the note states the stranded mass in both limits and explains why the adiabatic tank strands more. |
| **References** | `@libretexts-yan-engineering-thermodynamics` Ch. 4 "The First Law of Thermodynamics for Closed Systems": 4.4 Work (4.3.1 Boundary work; 4.3.2 Polytropic process and its boundary work); 4.5 The first law of thermodynamics for closed systems; Ch. 5 "The First Law of Thermodynamics for a Control Volume": 5.3 Mass and energy conservation equations in a control volume (5.2.1 Steady flow and transient flow; 5.2.2 Mass conservation equation; 5.2.3 Energy conservation equation) (required-free; https://eng.libretexts.org/Bookshelves/Mechanical_Engineering/Introduction_to_Engineering_Thermodynamics_(Yan)). `@mit-16-unified-thermodynamics-propulsion-notes` 2.1 First Law of Thermodynamics; 2.3.3 Transient filling of a tank (the mirror image of this lesson's blowdown); 2.5 Control volume form of the conservation laws; 14.3.1 Quasi-one-dimensional compressible flow in a variable area duct (the choked mass flow) (required-free, read online, https://web.mit.edu/16.unified/www/FALL/thermodynamics/notes/node17.html). `@openstax-university-physics-vol2` 3.3 First Law of Thermodynamics; 3.6 Adiabatic Processes for an Ideal Gas (required-free). `@lienhard-heat-transfer-textbook` Appendix A, Table A.5 "Thermophysical properties of saturated vapors (p ≠ 1 atm)", nitrogen (required-free; https://ahtt.mit.edu/): the saturation data for the condensation check. |
| **Extensions** | ★ Give the tank wall a heat capacity (a second ODE for T_wall, with an outer convective coefficient) and find how much of the isothermal limit a carbon-fibre COPV of realistic mass can provide. ★ Replace the ideal gas by a compressibility-factor model Z(p, T) fitted to data you find in a verified source, and see how the stranded mass changes. ★ Model the propellant tank as well: the ullage volume grows as propellant is expelled at a given flow rate, and the pressurant flow keeps the ullage pressure up; find the orifice size that holds the ullage at 0.5 MPa for a 100 s burn. |
| **Time** | unmeasured |

## 1. The engineering problem

A liquid rocket stage pushes its propellants out of their tanks with gas. In a pressure-fed stage, and in the tank-pressurisation system of many pump-fed ones, that gas comes from small, very high pressure bottles: composite-overwrapped pressure vessels (COPVs) holding helium or nitrogen at tens of megapascals. As the propellant drains, gas flows from the bottle through a regulator or a flow-limiting orifice into the tank ullage (the gas space above the liquid) and keeps its pressure up. Sizing the bottle is a first-law problem: how much gas can be drawn out before the bottle's pressure falls below what the regulator needs, and what happens to the temperature of the gas left behind?

The question is subtle because the bottle is a **control volume** with mass leaving it, and because the gas inside does work on the gas leaving. If no heat reaches it, the remaining gas expands and cools, sometimes by more than 100 K, which matters for seals, for the composite's strength and for the mass that can still be drawn out. If the wall supplies heat fast enough, the gas stays near the wall temperature. Real bottles are in between. This lesson writes the first law for that tank, integrates it, and checks it three ways.

Success criteria, written first: the energy balance U(t) + ∫ṁh dt − ∫Q̇ dt = U(0) closes to 10⁻¹²; the adiabatic tank state stays on p ∝ ρ^γ; the numerical blowdown agrees with the closed-form choked solution to 10⁻⁷; the isothermal case reproduces its exponential and its heat R T₀ Δm; and any finite wall conductance gives an answer between the two limits.

**No pressure vessel is built, filled or vented at home in this course.** Everything here is computed.

## 2. Understand and derive

### 2.1 The first law for a closed system

For a closed system (a fixed mass), energy conservation says the change in the system's energy equals the heat added minus the work done by the system:

**ΔU = Q − W**  (kinetic and potential energy of the whole system neglected).

Sign convention: Q > 0 into the system, W > 0 done by the system. Q and W are not properties: they depend on the path (`P07.01`, 2.2). U is a property: over a cycle ΔU = 0. Per unit mass, Δu = q − w.

**Boundary work.** A gas whose boundary moves by dV does work δW = p dV on its surroundings when the process is quasi-static. Per unit mass, w = ∫p dv. For an ideal gas along a **polytropic** path p vⁿ = constant:

**w = (p₂v₂ − p₁v₁)/(1 − n)  for n ≠ 1,  w = p₁v₁ ln(v₂/v₁)  for n = 1 (isothermal).**

n = γ is the reversible adiabatic (isentropic) process: then q = 0 and w = −Δu = cv(T₁ − T₂). The code's `polytropic_work` computes both, and the test checks that n → 1 recovers the logarithm. A worked example: 1 kg of N₂ at 1 MPa and 300 K (v₁ = 0.089041 m³/kg) expands to three times its volume.

| path | p₂ (kPa) | T₂ (K) | w (kJ/kg) | Δu (kJ/kg) | q = Δu + w (kJ/kg) |
|---|---|---|---|---|---|
| isothermal, n = 1 | 333.3 | 300.0 | 97.82 | 0 | 97.82 |
| polytropic, n = 1.2 | 267.6 | 240.8 | 87.82 | −43.91 | 43.91 |
| adiabatic, n = γ = 1.4 | 214.8 | 193.3 | 79.16 | −79.16 | 0 |

Same end volume, three different amounts of work and heat: path functions.

### 2.2 The first law for a control volume

For a control volume (CV) with mass crossing its boundary, two things change. Mass carries its own energy u + V²/2 + gz across, and the surroundings do **flow work** p v on each kilogram pushed in or out. Grouping u + p v = h, the energy balance is

**dE_cv/dt = Q̇ − Ẇ_shaft + Σ_in ṁ (h + V²/2 + g z) − Σ_out ṁ (h + V²/2 + g z).**

Mass conservation is dm_cv/dt = Σ_in ṁ − Σ_out ṁ. Enthalpy appears, not internal energy, because of the flow work: this is why h is the useful property for flowing systems (`P07.01`). In steady flow (dE_cv/dt = 0) with one inlet and one exit, no shaft work and small kinetic energy, the balance reduces to q = h₂ − h₁: the combustor equation `P07.01` used.

### 2.3 The tank as a control volume

Draw the boundary around the inside of the tank. **Assumptions**, named where they enter: the tank is rigid (no boundary work, no shaft work); its contents are uniform (one T and p everywhere: well mixed); the gas inside is nearly at rest (E_cv = U = m u); the gas leaving has the tank's state, so it carries the tank's enthalpy h = cp T (its later acceleration in the orifice converts h into kinetic energy, but h + V²/2 is the tank's h); no mass enters. Then

**dm/dt = −ṁ,  d(m u)/dt = −ṁ h + Q̇.**

Expand the second: m du/dt + u dm/dt = −ṁ h + Q̇, so

**m cv dT/dt = −ṁ (h − u) + Q̇ = −ṁ R T + Q̇.**

The gas left behind cools at the rate ṁ R T/(m cv) unless heat arrives. That term, h − u = p v = R T, is the flow work the remaining gas does pushing the departing gas out. Three heat cases close the model:

- **adiabatic**: Q̇ = 0 (a fast blowdown, or a well-insulated tank);
- **isothermal**: Q̇ = ṁ R T exactly cancels the cooling, so T = T₀ (a slow blowdown through a conducting wall);
- **finite wall conductance**: Q̇ = hA (T_wall − T), with hA in W/K lumping the inner convective coefficient and the wetted area.

The p of an ideal gas in the rigid tank is p = m R T/V = (γ − 1) U/V for constant cv: the internal energy per unit volume.

### 2.4 The choked orifice

The orifice is a short nozzle. Take the tank as its **stagnation state** (p₀, T₀: the state where the gas is at rest) and assume the flow from the tank to the throat is isentropic. Energy gives h₀ = h + V²/2, i.e. T₀/T = 1 + (γ − 1)M²/2 with M = V/a the Mach number and a = √(γ R T); the isentropic relation gives p/p₀ = (T/T₀)^(γ/(γ−1)). The mass flow per unit area ρV has a maximum where M = 1 at the throat; once the back pressure is low enough that the throat reaches M = 1, lowering it further changes nothing upstream: the flow is **choked**. With T* = 2T₀/(γ + 1), p* = p₀(2/(γ + 1))^(γ/(γ−1)), ρ* = p*/(R T*) and V* = a* = √(γ R T*),

**ṁ = Cd A ρ* a* = Cd A Γ p₀ √(γ/(R T₀)),  Γ = (2/(γ + 1))^((γ+1)/(2(γ−1))).**

For γ = 1.4, Γ = 0.5787 and the flow is choked while p_back/p₀ < **0.5283** (for helium, Γ = 0.5625 and the ratio is 0.4871). Below choking the isentropic subsonic formula, ṁ = Cd A p₀ √(2γ/((γ − 1)R T₀) [r^(2/γ) − r^((γ+1)/γ)]) with r = p_back/p₀, joins the choked one continuously at the critical ratio (the tests check it). The discharge coefficient Cd < 1 lumps the vena contracta and friction; 0.8 is a typical, illustrative value for a sharp-edged orifice. Note what choked flow means for the tank: ṁ ∝ p₀/√T₀, so a cold tank pushes out more mass per second at the same pressure.

### 2.5 The closed-form blowdowns

While the flow is choked the ODEs have closed-form solutions. **Adiabatic.** With Q̇ = 0 the remaining gas has expanded reversibly and adiabatically (each kilogram still inside has only done p dv work on the gas ahead of it), so p/p₀ = (ρ/ρ₀)^γ and T/T₀ = (ρ/ρ₀)^(γ−1). Write r = ρ/ρ₀ and substitute into V dρ/dt = −ṁ:

V ρ₀ dr/dt = −Cd A Γ √(γ/R) p₀ r^γ /(√T₀ r^((γ−1)/2)) = −Cd A Γ √(γ R T₀) ρ₀ r^((γ+1)/2).

Define the **blowdown time constant** **τ = V/(Cd A Γ a₀)**, a₀ = √(γ R T₀); then dr/dt = −r^((γ+1)/2)/τ. Separating variables (`P02.01`): r^(−(γ−1)/2) = 1 + (γ − 1)t/(2τ), so

**p/p₀ = (1 + (γ − 1) t/(2τ))^(−2γ/(γ−1)),  T/T₀ = (p/p₀)^((γ−1)/γ),  m/m₀ = (p/p₀)^(1/γ).**

**Isothermal.** With T = T₀, V dρ/dt = −Cd A Γ √(γ R T₀) ρ, so **p/p₀ = m/m₀ = exp(−t/τ)**: the same τ. In both cases τ = m₀/ṁ₀, the time the initial flow rate would take to empty the tank.

For the case of the engineering question (N₂, V = 50 L, d = 2 mm, Cd = 0.8, T₀ = 300 K): a₀ = 353.07 m/s, m₀ = p₀V/(R T₀) = **16.846 kg**, ṁ₀ = **173.02 g/s** and **τ = 97.37 s**. The flow stays choked until the tank falls to 0.5/0.5283 = **0.946 MPa**. Time to 3 MPa:

- adiabatic: t = (2τ/(γ − 1))[(p/p₀)^(−(γ−1)/(2γ)) − 1] = 486.9 × (10^0.142 86 − 1) = **189.62 s**, with T = 300 × 0.1^0.2857 = **155.38 K** and m = 0.1^(1/1.4) m₀ = **0.1931 m₀ = 3.2525 kg**;
- isothermal: t = τ ln 10 = **224.20 s**, with m = **0.1000 m₀ = 1.6846 kg** and heat in Q = R T₀ (m₀ − m) = **1350.0 kJ**.

The adiabatic tank reaches 3 MPa **sooner** and with **more gas left**: its pressure falls both because mass leaves and because the gas cools, so it reaches the regulator limit with 1.93 times the mass the isothermal tank still has. That gas is stranded. Close the valve and let it warm back to 300 K at constant volume, and the pressure recovers to p = 3 × 300/155.38 = **5.79 MPa**: the "pressure recovery" familiar from tank tests.

### 2.6 Why the energy balance closes to round-off

`simulate_blowdown` integrates four states: m, U, H_out = ∫ṁh dt and Q_in = ∫Q̇ dt. The combination U + H_out − Q_in has time derivative (−ṁh + Q̇) + ṁh − Q̇ = 0 identically: it is a **linear invariant** of the ODE system. Every Runge–Kutta method, including `solve_ivp`'s DOP853, is a linear combination of right-hand-side evaluations, so it preserves linear invariants exactly, step after step, whatever its truncation error. The residual |U + H_out − Q_in − U₀|/U₀ is therefore a test of the code (a sign error or a wrong h would break it), not of the step size, and it comes out at 2.5 × 10⁻¹⁶. The truncation error shows up elsewhere: in the agreement with the closed form (1.2 × 10⁻¹⁰ with rtol = 10⁻¹⁰).

### 2.7 Wall heat transfer and its time scale

A tank wall of conductance hA warms the gas with a thermal time constant **τ_th = m cv/(hA)**. At the start m₀ cv = p₀V/((γ − 1)T₀) = 12 500 J/K for this tank, so τ_th = 12 500/hA s. Compared with the blowdown time τ = 97 s: hA = 1 W/K gives τ_th = 12 500 s (essentially adiabatic); hA = 100 W/K gives 125 s (comparable, so neither limit holds). As the tank empties, m and with it τ_th fall, so heat transfer matters more and more late in the blowdown. A real COPV's hA depends on the free-convection coefficient inside (lesson 13) and the wall area; this lesson treats it as an input to sweep.

### 2.8 Where the model fails

- **Real gas.** At 30 MPa nitrogen is far from ideal; the mass a 50 L bottle holds at 30 MPa and 300 K is not exactly p V/(R T). Lesson 03 shows how a table measures the departure (Z); the trends here survive, the numbers shift.
- **Condensation.** Run the adiabatic blowdown all the way down to the 0.5 MPa back pressure and the model's gas reaches **93.1 K**. Lienhard's Table A.5 gives nitrogen's saturation pressure as 0.3605 MPa at 90 K and 0.7783 MPa at 100 K, which interpolates (ln p linear in 1/T) to a saturation temperature of about 94.0 K at 0.5 MPa. The model's end state is at the condensation line: a single-phase ideal gas is no longer the right model there. With helium (normal boiling point near 4 K) this cannot happen.
- **Uniform state.** A fast blowdown stratifies the gas (cold near the outlet, warmer near the wall); one temperature is an average.
- **Constant wall temperature.** A thin liner cools along with the gas; the wall's heat capacity sets how much heat it can give (Extension 1).
- **Constant Cd and back pressure.** Both vary in practice; the regulator itself is a throttle (lesson 04 counts its entropy generation).

## 3. Build the computational model

Open `starter/main.py`. The gas data (M and γ for N₂ and He), the `Tank` and `Orifice` dataclasses, the case of the engineering question, `gas_properties`, `initial_mass`, the integration driver `simulate_blowdown` and `main` are given.

### Step 1: the orifice

`critical_pressure_ratio(γ)` and `choked_flow_function(γ)` are one line each. `mass_flow(p0, T0, p_back, Cd, A, R, γ)` returns 0 if p₀ ≤ p_back, the choked formula if p_back/p₀ is at or below the critical ratio, and the subsonic formula otherwise (clip the bracket at zero to protect the square root from round-off).

### Step 2: the right-hand side

```python
def blowdown_rhs(t, y, tank, orifice, p_back, mode, hA=0.0, T_wall=None):
    R, g, cv, cp = gas_properties(tank.gas)
    m, U = y[0], y[1]
    T = U / (m * cv)
    p = m * R * T / tank.V
    mdot = mass_flow(p, T, p_back, orifice.Cd, orifice.A, R, g)
    Q = ...                  # 0, mdot * R * T, or hA * (T_wall - T)
    h = cp * T
    return [-mdot, -mdot * h + Q, mdot * h, Q]
```

The internal-energy datum is u = cv T (zero at 0 K), consistent with h = cp T; any other datum would work if used for both.

### Step 3: the given driver

`simulate_blowdown` calls `solve_ivp(blowdown_rhs, ..., method="DOP853", rtol=1e-10)` with a terminal event at p = p_stop (p = (γ − 1)U/V, so the event needs no extra call), and returns the states at the integrator's steps and on an even grid. Read it: it is the pattern of every transient balance in this phase.

### Step 4: the checks and closed forms

`energy_residual(res)` evaluates max |U + H_out − Q_in − U₀|/U₀ on `res["steps_y"]`. `blowdown_time_constant`, `adiabatic_choked_closed_form`, `isothermal_choked_closed_form` and `adiabatic_time_to_pressure` are 2.5. `remaining_mass_fraction`, `polytropic_work` and `recovery_pressure` are one or two lines each.

Running `main()` prints

```
COMPUTATIONAL MODEL ONLY: no pressure vessel is built, filled or vented at home in this course.
N2 tank V = 50 L, p0 = 30.0 MPa, T0 = 300 K: m0 = 16.8462 kg; orifice d = 2 mm, Cd = 0.8; back pressure 0.5 MPa
critical pressure ratio 0.5283: choked while p > 0.946 MPa; initial mdot = 173.02 g/s; tau = V/(Cd A Gamma a0) = 97.37 s
adiabatic: reaches 3 MPa at t = 189.62 s, T = 155.38 K, m = 3.2525 kg (19.31 % of m0 left), heat in 0.00 kJ
  energy balance residual max |U + int(mdot h) - int(Q) - U0| / U0 = 2.5e-16
isothermal: reaches 3 MPa at t = 224.20 s, T = 300.00 K, m = 1.6846 kg (10.00 % of m0 left), heat in 1350.00 kJ
  energy balance residual max |U + int(mdot h) - int(Q) - U0| / U0 = 5.0e-16
adiabatic closed form: t to 3 MPa = 189.62 s; max |p_num/p_closed - 1| = 1.2e-10
  p / rho^gamma varies by 9.6e-11 (isentropic remaining gas)
  remaining gas, closed system: u0 - u = 107.306 kJ/kg, boundary work int p dv = 107.306 kJ/kg
isothermal closed form: t to 3 MPa = tau ln(p0/p) = 224.20 s; max |p_num/p_closed - 1| = 1.7e-10; heat = R T0 (m0 - m) = 1350.00 kJ
stranded mass at 3 MPa: adiabatic 0.1931 m0, isothermal 0.1000 m0
after the valve closes the cold adiabatic gas warms back to 300 K: p rises from 3.0 to 5.79 MPa
```

then the wall-conductance table and the helium case (Section 4), and writes `tank_blowdown.png`: tank pressure against time (ODE and closed form, adiabatic and isothermal, with the unchoking pressure marked); gas temperature for adiabatic and two wall conductances; and the first-law residual against time.

### Verify

`python -m pytest curriculum/phase-07/lesson-02-first-law` runs 17 tests:

1. **Orifice**: critical ratios 0.528 282 (γ = 1.4) and 0.487 139 (5/3); mass flow continuous across the critical ratio to 10⁻⁶; choked flow independent of back pressure and ∝ p₀/√T₀; subsonic flow falling monotonically to zero; no reverse flow.
2. **Conservation**: the energy residual below 10⁻¹² for adiabatic, isothermal and wall heat transfer; mass falling monotonically; the enthalpy carried out bounded by U₀ plus the heat in.
3. **Isothermal limit**: T constant to 10⁻⁹ and Q_in = R T₀ (m₀ − m).
4. **Adiabatic tank**: p/ρ^γ constant to 10⁻⁸; the closed-system first law for the gas left behind, cv(T₀ − T) = ∫p dv, to 10⁻⁸; the polytropic work tending to the isothermal logarithm as n → 1.
5. **Closed forms**: p, T and m against the adiabatic closed form to 10⁻⁷ for N₂ and He up to 5 % above the unchoking pressure; the isothermal exponential to 10⁻⁷; τ = m₀/ṁ₀.
6. **Events and limits**: the stop event lands on p_stop at the closed-form time; wall heat transfer lies between the limits, monotone in hA, with hA = 0 reproducing the adiabatic case and very large hA the isothermal one; the remaining-mass fraction and recovery pressure; the long blowdown ending at the back pressure; helium's shorter τ.

### Compare

The independent routes agree to the integrator's tolerance: the ODE and the closed-form adiabatic blowdown to 1.2 × 10⁻¹⁰ in pressure; the ODE and the isothermal exponential to 1.7 × 10⁻¹⁰; the tank's p/ρ^γ constant to 9.6 × 10⁻¹¹; and the closed-system view of the gas left behind (107.306 kJ/kg of internal energy lost, 107.306 kJ/kg of boundary work done) matches the control-volume view exactly. The energy balance, by construction a linear invariant, closes to 2.5 × 10⁻¹⁶. Each of these checks would catch a different bug: the closed form a wrong mass-flow constant, p ∝ ρ^γ a wrong energy equation, the residual a wrong sign.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** No pressure vessel is built, filled or vented at home in this course; every result here is computed. Write each prediction down before running it.

1. **Between the limits.** Sweep the wall conductance hA = 0, 1, 3, 10, 30 and 100 W/K for the case of the engineering question. *Predict* how the time to 3 MPa, the final temperature and the stranded mass vary. *Finding:*

   | hA (W/K) | t to 3 MPa (s) | T at 3 MPa (K) | m left (kg) | energy residual |
   |---|---|---|---|---|
   | 0 | 189.62 | 155.38 | 3.2525 | 2.5 × 10⁻¹⁶ |
   | 1 | 191.27 | 157.92 | 3.2002 | 2.5 × 10⁻¹⁶ |
   | 3 | 194.57 | 163.19 | 3.0970 | 3.7 × 10⁻¹⁶ |
   | 10 | 205.71 | 183.16 | 2.7592 | 2.5 × 10⁻¹⁶ |
   | 30 | 225.90 | 234.70 | 2.1533 | 3.7 × 10⁻¹⁶ |
   | 100 | 232.20 | 282.61 | 1.7883 | 3.7 × 10⁻¹⁶ |

   The temperature and the stranded mass move monotonically from the adiabatic to the isothermal limit. The time does not: at hA = 30 and 100 W/K it exceeds the isothermal 224.20 s, and at 300 W/K it is 228.46 s. The gas first cools (to a minimum of 250.5 K at 89.6 s for hA = 100 W/K), and then the wall warms it again while the pressure is falling; that late re-warming holds the pressure up and stretches the time. A design rule taken from either limit alone would miss this.
2. **Helium against nitrogen.** Repeat the adiabatic case with helium in the same tank. *Predict* which gas empties faster and which gets colder. *Finding:* helium's m₀ is 2.4070 kg (one seventh), τ = 34.70 s (its speed of sound is 2.89 times nitrogen's), it reaches 3 MPa at 60.89 s at 119.43 K (larger γ, more cooling), strands 25.1 % of its mass and recovers to 7.54 MPa on warming. With hA = 10 and 100 W/K: 65.13 s at 132.38 K and 84.20 s at 250.11 K.
3. **Orifice size.** Repeat with d = 1.5, 2.0, 2.5 and 3.0 mm. *Predict* how τ scales (∝ 1/d²) and whether the stranded mass at 3 MPa changes in the adiabatic and isothermal limits (it should not: the remaining-mass fractions depend only on p/p₀ and γ), and how it changes with hA fixed (smaller orifice, longer blowdown, more time for heat to arrive). *Finding:* τ = 173.10, 97.37, 62.32 and 43.27 s; adiabatic times to 3 MPa 337.10, 189.62, 121.36 and 84.28 s and isothermal 398.57, 224.20, 143.49 and 99.64 s, with the stranded mass fixed at 3.2525 kg and 1.6846 kg; at hA = 10 W/K the gas ends at 205.79, 183.16, 172.65 and 167.14 K and strands 2.4559, 2.7592, 2.9271 and 3.0238 kg. A small orifice buys warmth and usable gas at the price of flow rate.
4. **Where the model breaks.** Run the adiabatic case to the back pressure (t_end = 3000 s). *Predict* the final temperature. *Finding:* the flow unchokes at 310.82 s at 111.75 K with 1.427 kg left; the tank reaches 0.5 MPa at 93.13 K, at nitrogen's condensation temperature for that pressure (about 94.0 K from Lienhard Table A.5). State what the model can and cannot say after that point.

Raw data format: CSV tables `hA_WK, t_stop_s, T_end_K, m_end_kg, residual` (1); `gas, hA_WK, t_stop_s, T_end_K, m_end_kg` (2); `d_mm, tau_s, t_stop_s, T_end_K, m_end_kg` (3); `t_s, p_Pa, T_K, m_kg` (4); plus `tank_blowdown.png`. Every file is labelled as computed (`COMPUTED` in its header line); none is a measurement.

## 5. Validate and improve

| quantity | independent value | ODE | difference | cause |
|---|---|---|---|---|
| energy balance | U + ∫ṁh − ∫Q̇ = U₀ | residual 2.5 × 10⁻¹⁶ (adiabatic), 5.0 × 10⁻¹⁶ (isothermal) | round-off | linear invariant, exact for RK |
| p(t), adiabatic, choked | closed form 2.5 | max relative 1.2 × 10⁻¹⁰ | — | rtol = 10⁻¹⁰ |
| time to 3 MPa, adiabatic | 189.62 s (closed form) | event at 189.62 s | < 10⁻⁷ | — |
| p(t), isothermal | exp(−t/τ) | 1.7 × 10⁻¹⁰ | — | rtol |
| isothermal heat | R T₀ (m₀ − m) = 1350.00 kJ | 1350.00 kJ | — | — |
| adiabatic tank state | p ∝ ρ^γ | variation 9.6 × 10⁻¹¹ | — | rtol |
| gas left behind | Δu = ∫p dv = 107.306 kJ/kg | 107.306 kJ/kg | — | closed system inside the CV |
| finite hA | between the limits | monotone in T and m | — | — |

The engineering answer: **the 50 L, 30 MPa nitrogen tank reaches the 3 MPa regulator limit after 189.6 s if no heat reaches the gas, which by then is at 155.4 K, leaving 3.25 kg (19.3 %) stranded; with an isothermal wall it takes 224.2 s and strands 1.68 kg (10.0 %). A real tank lies between, set by hA: about 2.8 kg stranded at 10 W/K.** The model's numbers are as good as its ideal-gas assumption at 30 MPa and its single uniform temperature; the first thing to improve is the equation of state, the second the wall's heat capacity (Extension 1).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of the tank ODE, the choked mass flow and the closed-form adiabatic blowdown; `tank_blowdown.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) answering the engineering question.

### Adversarial self-checks

- **What if the orifice area doubles?** τ halves; the closed-form times to any pressure halve; the adiabatic and isothermal end states at 3 MPa do not change (they depend on p/p₀ and γ only); with a finite hA the gas has less time to warm, so it ends colder and strands more.
- **What if the back pressure goes to zero?** The flow stays choked all the way down; the closed forms then hold for the whole blowdown; the adiabatic temperature keeps falling, and the condensation limit arrives sooner.
- **What are the limiting cases?** hA = 0 is adiabatic, hA → ∞ isothermal; γ → 1 makes the two limits coincide (no cooling by expansion); a zero-area orifice gives no change; p_back ≥ p gives no flow.
- **Which assumption fails first in the real world?** The ideal gas at 30 MPa, then the uniform temperature in a fast blowdown, then the constant wall temperature.
- **What measurement could falsify the model?** A tank test with gas-temperature and pressure transducers: an adiabatic-looking blowdown whose p and T do not satisfy p ∝ T^(γ/(γ−1)), or a stranded mass that does not match m = pV/(RT) at the measured end state, would point at the equation of state or at stratification.
- **Does the result survive a change of integrator or tolerance?** The energy residual is round-off at any tolerance (linear invariant); the agreement with the closed form scales with rtol (10⁻¹⁰ here). Switching to RK45 changes neither conclusion.

## Check yourself

1. Write the first law for a closed system and say what sign convention you use for Q and W.
2. Why does enthalpy, not internal energy, appear in the energy carried by a flowing stream?
3. Write the energy equation of a rigid tank with one exit and heat transfer Q̇. What makes the gas inside cool?
4. What is the critical pressure ratio for γ = 1.4, and what does "choked" mean for the tank upstream?
5. Why does an adiabatic blowdown reach a given pressure sooner than an isothermal one, and strand more mass?
6. Give the isothermal blowdown law p(t) and the meaning of τ.
7. Why does `energy_residual` come out at round-off whatever the step size?
8. What happens to the tank pressure after the valve closes on a cold, adiabatically blown-down tank?
