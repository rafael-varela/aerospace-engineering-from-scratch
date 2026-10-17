# P04.08  Electrostatics & Circuits

> Charge is conserved and the electric field is conservative. At a node, the first fact is Kirchhoff's current law; around a loop, the second is his voltage law. Write both for every node, solve one linear system, and you have every voltage and current in any DC network.

| | |
|---|---|
| **Learning outcomes** | 1. Compute forces, fields and potentials of point-charge distributions by superposition, and verify E = −∇V numerically. 2. Compute the resistance of a conductor from its resistivity and geometry, and apply Ohm's law, the power law and the series, parallel and divider rules. 3. Derive Kirchhoff's current and voltage laws from charge conservation and the conservative field, formulate them as modified nodal analysis, and implement and test a DC circuit solver that handles resistors, voltage sources and current sources, with KCL residuals and the power balance as checks. 4. Predict and correct the loading error of a voltmeter and the burden of an ammeter. 5. Derive the RC charging transient, verify it against `scipy.integrate.solve_ivp`, and show that charging a capacitor from a fixed source always dissipates half the energy drawn, whatever the resistor. |
| **Prerequisites** | `P01.13` (systems of linear equations), `P03.08` (direct linear solvers, conditioning), `P02.01` (first-order linear ODEs), `P01.07` (vectors), `P04.02` (work, energy, power), `P04.05` (inverse-square fields and potentials by analogy). Skill check: you can write and solve a 3 × 3 linear system with `numpy.linalg.solve` and solve y′ = (a − y)/τ by hand. |
| **Engineering question** | On a 28 V satellite power bus feeding 220 W of loads through an 8 m harness, what voltage reaches the loads and how much power the wiring wastes for a given wire size, what wire diameter keeps the drop under 2 %, and how far can you trust a multimeter reading of that circuit? |
| **Assumptions and limits** | Electrostatics: point charges in vacuum. Circuits: steady DC, lumped ideal elements (linear resistors at a fixed temperature, ideal sources, ideal wires except where a wire is modelled as a resistor), node 0 as ground. Copper's resistivity rises by about 0.4 % per kelvin, so a warm harness has more resistance than at 20 °C; real loads on a regulated bus draw constant power rather than behaving as fixed resistors, which makes the drop slightly worse than the resistive model predicts (Problem 4); the RC model ignores the source's internal resistance and the capacitor's leakage and series resistance. |
| **Analytical work** | Coulomb's law, field, potential and their relation; the dipole far field; R = ρL/A; Ohm's law and P = VI; series, parallel and the divider; KCL and KVL and their matrix form (modified nodal analysis); the Wheatstone bridge and its balance; meter loading; the RC transient and its energy accounting; wire sizing for a drop. |
| **Computational deliverable** | `starter/main.py`: `coulomb_force`, `electric_field`, `electric_potential`, `field_from_potential`, `dipole_moment`; `resistance`, `wire_resistance`, `series`, `parallel`, `voltage_divider`; the `Circuit` model with `build_mna`, `solve_dc`, `branch_current`, `kcl_residuals`, `power_balance`, `wheatstone_bridge`, `voltmeter_reading`, `ammeter_reading`; the RC transient (`rc_rhs`, `rc_capacitor_voltage`, `rc_current`, `rc_energy_accounting`); the bus (`bus_circuit`, `load_resistance_for_power`, `wire_diameter_for_drop`). Tests (15) check Coulomb's law and Newton's third law, E = −∇V and superposition, the dipole's 1/z³ far field, resistance scaling and the combination rules, the solver on a divider, KCL and the power balance on a random network, the bridge against a hand-written loop-current solution, linearity and superposition, rejection of singular circuits, voltmeter loading and ammeter burden in their limits, the RC transient against `solve_ivp`, the half-energy rule, and the bus drop and wire sizing. |
| **Practical activity** | **computational experiment**: a wire-size and load sweep of the 28 V bus model, a meter-loading study across meter resistances, and an RC energy study across resistors. *Optional* (not required to pass): measure a divider and an RC time constant on a breadboard with the `projects/p04-electronics-bench` kit of `P04.15`, under its safety box. |
| **Equipment** | None required. Optional: the `projects/p04-electronics-bench` kit (breadboard, resistors, a ≤ 25 V, ≤ 1000 μF electrolytic capacitor, a multimeter, a 9 V battery or USB 5 V supply); see that project's `bom.yaml`. |
| **Safety** | None beyond normal computer use for the required work. The optional bench measurement follows `projects/p04-electronics-bench/safety.md` (summary box in part 4). |
| **Evidence** | `dc_circuits.png`; the bus table (wire diameter, load voltage, current, harness loss, drop); the meter-loading table; the RC energy table; a half-page note answering the engineering question. If you did the optional measurement: the raw readings and their comparison with the model, labelled as measurement. |
| **Assessment** | `assessments/phase-04/grading-rubric.md`; competencies 1–4. Pass: KCL/KVL written as the MNA system for the bridge by hand; all public tests pass; the 28 V bus drop at one wire size and the 2 % wire diameter reproduced by hand to three figures; the note states the wire diameter, the harness loss and the expected multimeter reading with its loading error. |
| **References** | `@openstax-university-physics-vol2` Ch. 5 "Electric Charges and Fields": 5.1 "Electric Charge", 5.3 "Coulomb's Law", 5.4 "Electric Field"; Ch. 7 "Electric Potential": 7.2 "Electric Potential and Potential Difference", 7.5 "Equipotential Surfaces and Conductors"; Ch. 9 "Current and Resistance": 9.1 "Electrical Current", 9.3 "Resistivity and Resistance", 9.4 "Ohm's Law", 9.5 "Electrical Energy and Power"; Ch. 10 "Direct-Current Circuits": 10.1 "Electromotive Force", 10.2 "Resistors in Series and Parallel", 10.3 "Kirchhoff's Rules", 10.4 "Electrical Measuring Instruments", 10.5 "RC Circuits" (all required-free; https://openstax.org/books/university-physics-volume-2/pages/10-3-kirchhoffs-rules). |
| **Extensions** | ★ Add temperature to the harness: R(T) = R₂₀(1 + 0.0039(T − 20 °C)) and the wire's own I²R heating with a radiative heat balance; iterate the circuit and the temperature together. ★ Extend the solver with capacitors (a backward-Euler companion model: a conductance C/Δt in parallel with a current source) and simulate the bus's inrush transient. ★ Model a battery as an EMF with internal resistance and find the load resistance that maximises delivered power. |
| **Time** | unmeasured |

## 1. The engineering problem

Every spacecraft has a power bus: a regulated voltage (28 V is a common standard; 50 V and 100 V for larger vehicles) distributed through a wiring harness to the loads, the transmitter, the computers, the heaters, the reaction wheels. The harness is a set of long copper wires, and every metre has resistance. Current through resistance drops voltage, so the loads see less than 28 V, and dissipates power, which the solar arrays had to generate and the thermal system has to remove. Thicker wire drops less and wastes less but weighs more, and on a spacecraft every kilogram costs. The trade between drop, loss and mass is decided by Ohm's law and Kirchhoff's rules, and checking it on the real hardware is done with a multimeter that itself changes the circuit it measures.

This lesson starts where electricity starts: charge, Coulomb's law, the field and the potential. The potential is what makes circuits computable. Because the electrostatic field is conservative, a potential is defined at every node of a circuit, and the sum of voltage changes around any loop is zero; because charge is conserved, the currents into any node sum to zero. Those two statements, applied to every node at once, are a linear system, and solving it is the same Gaussian elimination you wrote in `P01.13`. You build a general DC solver, test it against hand analysis and conservation of power, use it to size the bus wiring, and use it again to model the multimeter.

Success criteria, written first: the solver satisfies KCL to 10⁻¹² A and balances delivered and dissipated power to 10⁻¹⁰ on any network; it reproduces the hand-solved two-loop bridge to 10⁻¹⁰; it rejects a circuit with a floating node; the RC transient agrees with `solve_ivp` to 10⁻⁹ V; the 2 %-drop wire diameter is reproduced by hand.

## 2. Understand and derive

### 2.1 Charge and Coulomb's law

Electric charge comes in two signs, is quantised in units of e = 1.602 × 10⁻¹⁹ C, and is **conserved**: it moves but is never created or destroyed. Two point charges q₁ and q₂ a distance r apart exert forces on each other along the line joining them,

**F = k q₁ q₂ / r²**,  k = 1/(4πε₀) = 8.988 × 10⁹ N m²/C²,

repulsive for like signs, attractive for unlike, equal and opposite (OpenStax 5.3). The structure is gravity's (`P04.05`) with charge for mass and a sign. The magnitudes are not alike: an electron and a proton at the Bohr radius (5.29 × 10⁻¹¹ m) attract with 8.2 × 10⁻⁸ N, about 10³⁹ times their gravitational attraction.

### 2.2 Field and potential

The **electric field** E at a point is the force per unit positive test charge: for a point charge, E = kq r̂/r². For many charges the fields add (**superposition**): E(r) = Σ kqᵢ(r − rᵢ)/|r − rᵢ|³. The electrostatic force is conservative, so it has a potential energy, and the **electric potential** (potential energy per unit charge) is

**V(r) = Σ k qᵢ / |r − rᵢ|**,  with **E = −∇V**.

Potential is measured in volts (J/C); field in V/m. The potential is a scalar, so it adds without vector bookkeeping, and the field follows by differentiation. `field_from_potential` differentiates `electric_potential` numerically and the tests compare it with `electric_field` at random points: they agree to 10⁻⁶. A pair of equal and opposite charges ±q separated by d is a **dipole**, with moment p = qd pointing from − to +. Far away its field falls as **1/r³**, not 1/r², because the two charges' fields nearly cancel. On the axis it is 2kp/z³, measured from the dipole's centre. With q = 1 nC and d = 1 mm (p = 10⁻¹² C m), at 50 mm from the positive charge (49.5 mm from the centre) the code gives 148.2 V/m, and 2kp/z³ with z = 49.5 mm gives the same; with z = 50 mm it gives 143.8 V/m, a 3 % error from measuring z from the wrong point. Its potential is zero everywhere on the plane bisecting the pair (an **equipotential**).

### 2.3 Current, resistance and Ohm's law

A **current** I is charge per unit time through a cross-section (1 A = 1 C/s). In a conductor an applied field drives the free electrons, which collide constantly with the lattice and drift at a steady average speed. The result, for most metals at fixed temperature, is that the current density is proportional to the field: **Ohm's law**. For a uniform wire of length L and cross-section A it becomes

**V = I R**,  **R = ρ L / A**,

with ρ the resistivity (copper 1.68 × 10⁻⁸ Ω m at 20 °C; aluminium 2.65 × 10⁻⁸ Ω m). R grows with length and falls with area, so as 1/d² for a round wire. Ten metres of 1.291 mm copper (AWG 16, a common harness gauge) is 0.1283 Ω; the same length at 0.5 mm is 0.856 Ω. The power delivered to an element with voltage V across it and current I through it is **P = VI**, and for a resistor **P = I²R = V²/R** (OpenStax 9.5).

Two combination rules follow from the laws below. **Series** (same current through each): R = ΣRᵢ. **Parallel** (same voltage across each): 1/R = Σ1/Rᵢ. Two resistors in series across V_s form a **voltage divider**, with V_s R₂/(R₁ + R₂) across R₂: 12 V across 10 kΩ and 4.7 kΩ gives 3.8367 V.

### 2.4 Kirchhoff's rules

**Current law (KCL).** Charge is conserved and does not accumulate at a junction of wires in steady state, so the currents flowing out of any node sum to zero.

**Voltage law (KVL).** The electrostatic field is conservative, so its line integral around any closed path is zero: the potential differences around any loop sum to zero. Equivalently, every node has a single, well-defined potential.

Series–parallel reduction handles many circuits, but not all. The **Wheatstone bridge** (OpenStax 10.4 uses it as a measuring instrument) has a resistor across the middle that is neither in series nor in parallel with anything. With a 5 V source, arms of 1 kΩ/2 kΩ and 1 kΩ/2.2 kΩ, and a 500 Ω galvanometer across the middle, you must write Kirchhoff's rules. Two ways:

- **Loop currents.** Assign a circulating current to each of the three independent loops and write KVL around each: three equations, three unknowns (Problem 2 does it; the test does it too).
- **Node voltages.** Assign a potential to each node (ground = 0) and write KCL at each: one equation per unknown node. With a voltage source the source fixes one node directly.

The node method generalises mechanically to any circuit, and it is what circuit simulators use.

### 2.5 Modified nodal analysis

Number the nodes 0 (ground), 1, …, N − 1. A resistor R between nodes a and b carries current (V_a − V_b)/R from a to b, which leaves node a and enters node b. KCL at node a, summing over its resistors, is linear in the node voltages: Σ (V_a − V_b)/R = (current injected into a by sources). In matrix form, **G V = i**, where the conductance matrix G gets +1/R on the diagonals at a and b and −1/R in the two off-diagonal positions, for every resistor. That is "stamping": each component adds its own small pattern to the matrix.

A current source I from node a to b adds −I to the right-hand side at a and +I at b. A voltage source cannot be written as a conductance, so **modified** nodal analysis adds its current as an extra unknown and its constraint V_plus − V_minus = V as an extra equation:

```
[ G   B ] [ V ]   [ i ]
[ Bᵀ  0 ] [ I ] = [ e ]
```

B has a +1 at the source's positive node and a −1 at its negative node. For a network of resistors and sources the matrix is symmetric; it is non-singular when every node has a path to ground and no loop is made of voltage sources alone. **A floating node** (one connected to nothing that leads to ground) makes G singular: its potential is undefined, and the solver must say so instead of returning nonsense. `solve_dc` checks the rank and raises `ValueError`.

Two independent checks apply to every solution. KCL residuals: the net current out of every node, recomputed from the solved voltages, must be zero. Power balance: the power delivered by the sources must equal the power dissipated in the resistors (energy conservation). For the bridge: the galvanometer carries 56.18 μA from the 2.2 kΩ side (node 3, at 3.39888 V) to the 2 kΩ side (node 2, at 3.37079 V); KCL residuals are 10⁻¹⁸ A; the source delivers 16.1517 mW and the resistors dissipate 16.1517 mW. With R_x = 2 kΩ the bridge balances (R₁/R₂ = R₃/R_x) and the galvanometer current is exactly zero, whatever its resistance: that is how a bridge measures an unknown resistance, and how a strain-gauge bridge on a wing spar works.

### 2.6 Meters change what they measure

A **voltmeter** is connected in parallel with what it measures and draws a small current through its input resistance R_m (10 MΩ for a typical digital multimeter). The circuit it is attached to looks, from the meter's terminals, like an EMF V_th behind a resistance R_th (Thévenin's theorem), so the meter reads

V_read = V_th R_m/(R_m + R_th),

low by a fraction R_th/(R_m + R_th). Across the middle of a divider of two 100 kΩ resistors on 12 V (true value 6.0000 V, R_th = 50 kΩ): a 10 MΩ meter reads 5.9701 V (0.5 % low), a 1 MΩ meter 5.7143 V, and a 100 kΩ meter 4.0000 V. The rule: the meter's resistance must be much larger than the circuit's resistance at the point measured, or you must correct for it. An **ammeter** is connected in series and adds its own resistance (the "burden"). 5 V across 220 Ω is 22.727 mA; with a 5 Ω ammeter in the loop it is 22.222 mA, and the meter drops 0.111 V. OpenStax 10.4 describes both instruments.

### 2.7 The RC transient

A capacitor C stores charge Q = CV on its plates and energy ½CV². Connect an uncharged one through a resistor R to a source V_s at t = 0. KVL: V_s = IR + Q/C, and I = dQ/dt, so

**dV_C/dt = (V_s − V_C)/(RC)**,

a first-order linear ODE (`P02.01`) with solution

**V_C(t) = V_s(1 − e^{−t/τ})**,  **I(t) = (V_s/R) e^{−t/τ}**,  **τ = RC**.

After one time constant the capacitor is at 1 − e⁻¹ = 63.21 % of V_s; after 5τ, 99.33 %. With R = 10 kΩ and C = 100 μF, τ = 1.00 s, and the closed form agrees with DOP853 to 3.3 × 10⁻¹¹ V. Discharging (V_s = 0, starting from V₀) is the same equation: V_C = V₀e^{−t/τ}.

**Where the energy goes.** The source delivers ∫V_s I dt = V_s Q_final = CV_s². The capacitor stores ½CV_s². The resistor dissipates ∫I²R dt = (V_s²/R)∫e^{−2t/τ}dt = ½CV_s². **Exactly half the energy is lost in the resistor, independent of R.** A small resistor gives a short, intense current pulse; a large one gives a long, weak one; the integral is the same. The code integrates the power numerically to 20τ: 2.5000 mJ from the source, 1.2500 mJ stored and 1.2500 mJ dissipated for 5 V and 100 μF, and the same with ten times the resistance. On a spacecraft, the inrush current when a capacitor bank is switched onto the bus can trip protection; limiting it with a resistor lowers the peak, but not the energy lost.

### 2.8 The 28 V bus

Model the harness as a supply wire and a return wire, each of resistance R_w, with the loads in parallel between them. The loads are 120, 60 and 40 W at 28 V (R = V²/P: 6.53, 13.07 and 19.60 Ω). With an 8 m harness, the solver gives:

| wire diameter (mm) | R_w each (Ω) | V at loads (V) | current (A) | load power (W) | harness loss (W) | drop |
|---|---|---|---|---|---|---|
| 0.500 | 0.6845 | 20.229 | 5.676 | 114.83 | 44.11 | 27.75 % |
| 0.800 | 0.2674 | 24.347 | 6.832 | 166.33 | 24.96 | 13.05 % |
| 1.291 | 0.1027 | 26.474 | 7.429 | 196.68 | 11.33 | 5.45 % |
| 2.053 | 0.0406 | 27.376 | 7.682 | 210.31 | 4.79 | 2.23 % |

A drop limit of 2 % at the nominal current I = 220/28 = 7.86 A requires 2R_w ≤ 0.02 × 28/7.86 Ω, and with R_w = ρL/(πd²/4), **d = 2.19 mm** in copper (2.75 mm in aluminium, which is lighter per unit of conductance). Fixed-resistance loads slightly flatter the result (their current falls as the voltage falls); real regulated loads draw constant power and draw *more* current as the voltage falls (Problem 4).

### 2.9 Where the models fail

- **Temperature.** Copper's resistance rises by 0.39 % per kelvin; a harness at 60 °C has 16 % more resistance than at 20 °C.
- **Constant-power loads.** DC–DC converters hold their output power, so their input behaves like a *negative* incremental resistance; with long wiring that can make the bus voltage unstable, a problem the linear solver cannot show.
- **Non-ohmic elements.** Diodes, LEDs and transistors are not resistors; `P04.10` introduces them.
- **Time dependence.** The DC solver ignores inductance and capacitance; switching transients need the time-domain model (the ★ extension) or `P04.10`'s AC analysis.
- **Charge in vacuum.** Spacecraft surfaces charge up in the plasma environment, to kilovolts in geostationary orbit during eclipse, and discharges damage electronics. Coulomb's law and the potential of a sphere estimate the charge involved (Problem 3).

## 3. Build the computational model

Open `starter/main.py`. The constants and the component dataclasses (`Resistor`, `VoltageSource`, `CurrentSource`, `Circuit`, `DCSolution`) are given; the functions are yours.

### Step 1: electrostatics

`coulomb_force`, `electric_field`, `electric_potential`, `field_from_potential` (central differences), `dipole_moment`. Write the field as an explicit loop over charges first; then, if you want, vectorise it.

### Step 2: resistance and the combination rules

`resistance`, `wire_resistance`, `power_dissipated`, `series`, `parallel`, `voltage_divider`: one line each.

### Step 3: the solver

```python
def build_mna(circuit):
    N, vs = circuit.n_nodes, circuit.voltage_sources
    A = np.zeros((N - 1 + len(vs), N - 1 + len(vs)))
    z = np.zeros(N - 1 + len(vs))
    for r in circuit.resistors:                       # stamp each conductance
        g = 1.0 / r.R
        for a, b, sgn in ((r.n1, r.n1, 1), (r.n2, r.n2, 1), (r.n1, r.n2, -1), (r.n2, r.n1, -1)):
            if a and b:
                A[a - 1, b - 1] += sgn * g
    ...                                               # current sources into z; voltage sources as extra rows
```

Ground has no row or column (its potential is fixed at zero). `solve_dc` checks the rank, solves, prepends V₀ = 0, and reports each voltage source's current as the current it delivers (out of its positive terminal into the circuit). `branch_current`, `kcl_residuals` and `power_balance` are the checks of 2.5. `wheatstone_bridge`, `voltmeter_reading` (add a resistor across the two nodes) and `ammeter_reading` (add the meter's resistance in series with the chosen resistor, without modifying the original circuit) reuse the solver.

### Step 4: the RC transient and the bus

`rc_rhs` for `solve_ivp`, the closed forms, and `rc_energy_accounting` (trapezoidal integrals of V_sI and I²R). `bus_circuit` builds the four-node bus of 2.8; `load_resistance_for_power` and `wire_diameter_for_drop` are 2.8's formulas.

Running `main()` prints the dipole and Coulomb checks, the wire resistances and the divider, the bridge with its KCL and power checks, the meter-loading readings, the RC checks and energy accounting, and the bus table, and writes `dc_circuits.png` (the RC charging curve against `solve_ivp`, the voltmeter reading against its resistance, and the bus drop against wire diameter).

### Verify

`python -m pytest curriculum/phase-04/lesson-08-electrostatics-circuits` runs 15 tests. Coulomb: repulsion, attraction, 1/r², third law. Fields: E = −∇V at a random point to 10⁻⁶, superposition, a single charge's field and potential. Dipole: 1/z³ (ratio 8 for doubled distance), 2kp/z³, zero potential on the bisecting plane. Resistance: L, 1/d², AWG 16 at about 13 mΩ/m, the combination rules and the divider. Solver: the divider and Ohm's law, a symmetric MNA matrix; KCL to 10⁻¹² and the power balance to 10⁻¹⁰ on a random 5-node network with a voltage and a current source; the bridge against an independent three-loop KVL solution, its balance at any galvanometer resistance and the sign change across balance; linearity and superposition; floating node and parallel ideal sources rejected. Meters: no loading at infinite resistance, the exact loaded value, the 1/R_m scaling of the error, no loading across an ideal source; the ammeter burden, leaving the original circuit unchanged. RC: closed form against `solve_ivp`, 1 − e⁻¹ at τ, discharge, I = (V_s − V_C)/R. Energy: drawn CV², stored and dissipated ½CV² each, at three resistances. Bus: no wire, no drop; KVL around the harness; the series-parallel formula; the sizing formula and aluminium needing a larger diameter.

### Compare

The bridge test solves the same circuit by an independent method (loop currents and KVL, written by hand in the test) and the two agree to 10⁻¹⁰. The power balance compares energy in with energy out, which a sign error in a stamp would break. For the RC transient the oracle is `solve_ivp`. A professional circuit simulator (SPICE and its descendants) uses exactly modified nodal analysis with the same stamps, adds nonlinear elements through Newton's method (`P03.10`) and time dependence through implicit integrators (`P03.06`); for resistive DC circuits it would return the same numbers.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Wire size.** Sweep the bus wire diameter from 0.4 to 2.5 mm and record load voltage, current, load power, harness loss and drop. Then compute the harness copper mass (density 8960 kg/m³, supply plus return). *Predict:* drop ∝ 1/d² at small drop; the 2 % point at 2.19 mm; the mass ∝ d².
2. **Meter loading.** For the 100 kΩ/100 kΩ divider, sweep R_m from 10 kΩ to 100 MΩ. *Predict:* the error falls as R_th/R_m once R_m ≫ R_th, and is 0.5 % at 10 MΩ. Repeat for a 10 kΩ/10 kΩ divider and predict a tenfold smaller error.
3. **RC energy.** For C = 100 μF and V_s = 5 V, compute the energy dissipated for R = 100 Ω, 10 kΩ and 1 MΩ, integrating to 25τ. *Predict:* 1.25 mJ in every case.
4. **Linearity.** Double every source in the bridge and check that every node voltage doubles; then add a current source and verify superposition.

Raw data format: CSV files `d_mm, R_wire_ohm, V_load_V, I_A, P_load_W, P_loss_W, drop_pct, mass_kg` (1), `R_th_ohm, R_m_ohm, V_read_V, error_pct` (2), `R_ohm, E_source_J, E_stored_J, E_dissipated_J` (3); plus `dc_circuits.png`.

**Optional measurement (not required to pass).** If you have the `projects/p04-electronics-bench` kit from `P04.15`, build two circuits on the breadboard: a 100 kΩ/100 kΩ divider from a 9 V battery, measured with your multimeter, and an RC circuit (e.g. 10 kΩ and 100 μF, τ = 1 s) whose capacitor voltage you read every second after connecting it. Predict both readings first with this lesson's code, including your meter's input resistance (from its manual) for the divider. This is a measurement, so it carries uncertainty (resistor tolerance, typically ±1 % or ±5 %; capacitor tolerance, often ±20 %); compare with E_n, not by eye. Label the results as measured and keep them separate from the computed tables.

> **Safety (summary, for the optional measurement only; the authoritative box is `projects/p04-electronics-bench/safety.md`).** Battery or USB power only, never mains; supply ≤ 9 V DC; current limited to what the supply can give (≤ 500 mA); no lithium-cell charging circuits; every LED through a resistor; no capacitors above 25 V or 1000 μF; respect electrolytic capacitor polarity; discharge capacitors through a resistor before touching or rewiring; no soldering is needed (if you solder, the soldering-iron rules in that box apply). Stop if anything gets warm, smells, or a capacitor bulges.

## 5. Validate and improve

| quantity | theory (§2) | code | difference | cause |
|---|---|---|---|---|
| divider 12 V, 10 k / 4.7 k | 3.8367 V | 3.8367 V (solver) | 0 | — |
| bridge galvanometer current | loop-current KVL solution | 56.18 μA, same to 10⁻¹⁰ | 0 | independent methods agree |
| bridge KCL residual | 0 | 10⁻¹⁸ A | rounding | — |
| bridge power | delivered = dissipated | 16.1517 = 16.1517 mW | < 10⁻¹⁰ | — |
| balanced bridge | 0 | 0 | 0 | — |
| 10 MΩ voltmeter on 100 k/100 k | 6 × 10/(10 + 0.05) = 5.9701 V | 5.9701 V | 0 | — |
| RC at τ | 63.21 % | 0.6321 V_s | 0 | closed form |
| RC closed form vs DOP853 | 0 | 3.3 × 10⁻¹¹ V | — | solver tolerance |
| RC energy dissipated | ½CV² = 1.25 mJ | 1.2500 mJ at R and 10R | < 10⁻⁶ | trapezoid at 20 001 points |
| bus 2 % wire | 2.19 mm (formula at the nominal 7.86 A) | solver drop 1.96 % at 2.19 mm, loss 4.23 W | 0.04 % | resistive loads draw less than nominal current at reduced voltage |

Every computed result matches its independent check. The model's real uncertainties are in its inputs, not its arithmetic: copper at 60 °C instead of 20 °C raises R by 16 %; the loads are constant-power, which raises the current at reduced voltage; connectors and crimps add tens of milliohms each. For a flight design, add a margin to the drop budget or size at the hot-case resistance, and model the loads as constant-power (Problem 4 iterates it). What not to improve: the solver itself, which is exact for this model to rounding.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the MNA system for the bridge written by hand (matrix and right-hand side) and solved; the RC energy derivation; `dc_circuits.png`; the experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6, 7) answering the engineering question: the wire diameter for a 2 % drop, its harness loss and copper mass, the drop at 60 °C, and what a 10 MΩ multimeter would read at the load terminals and why the reading can be trusted there.

### Adversarial self-checks

- **What if the load power doubles?** The current doubles, the drop doubles (and more, with constant-power loads), the harness loss quadruples (I²R). The 2 %-drop diameter grows by √2.
- **What if the wire resistance goes to zero?** The load sees the full 28 V and the loss vanishes; the solver must return exactly the source voltage, which a test checks.
- **Which assumption fails first in the real world?** The fixed-temperature resistance (a 40 K rise is +16 %), then the constant-resistance load model.
- **What measurement could falsify the model?** A measured drop at a known current that differs from 2IR_w by more than the resistance's uncertainty: look for a bad crimp or connector (an extra series resistance) before doubting Ohm's law.
- **Does the result survive a change in method?** Series–parallel reduction by hand, MNA, and (for the bridge) loop currents give the same numbers; the power balance holds.
- **Which decision changes if the resistivity uncertainty doubles?** The wire diameter needed for 2 % moves by half the relative change in ρ (d ∝ √ρ); a design sized with no margin might then fail its drop budget, which is why harnesses are sized at the hot case.

## Check yourself

1. Write Coulomb's law. How does the force change when the distance halves? When one charge changes sign?
2. What is the relation between electric field and potential? In what units is each measured?
3. How does a round wire's resistance change when its diameter doubles? When its length doubles?
4. State Kirchhoff's two rules and the conservation law behind each.
5. Why does a floating node make the nodal-analysis matrix singular?
6. A 10 MΩ voltmeter measures the midpoint of a 1 MΩ/1 MΩ divider on 3.3 V. What does it read, and what is the true value?
7. What is the time constant of a 4.7 kΩ resistor and a 47 μF capacitor, and what fraction of the final voltage is reached after one time constant?
8. Charging a capacitor to V through a resistor R: what fraction of the energy drawn from the source is lost, and how does it depend on R?
