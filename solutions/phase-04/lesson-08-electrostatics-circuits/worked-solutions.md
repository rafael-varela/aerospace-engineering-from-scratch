# P04.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numbers come from `solution.py` and probe calculations that call it.

## 1. A sensor lead and a divider ($s = 15$: R₂ = 25 kΩ)

**Assumptions.** Copper at 20 °C; ideal source; resistors at nominal value.

(a) A = π(0.644 × 10⁻³)²/4 = 3.257 × 10⁻⁷ m²; R = 1.68 × 10⁻⁸ × 2.0/A = **0.1032 Ω**, **51.6 mΩ/m**; drop 0.5 × 0.1032 = **0.0516 V**; power 0.5² × 0.1032 = **25.8 mW**.

(b) V_out = 12 × 25/35 = **8.5714 V**; I = 12/35 kΩ = **0.3429 mA**; P = 12 × I = **4.114 mW**; R_th = 10 ∥ 25 kΩ = **7.143 kΩ**.

**Check.** Power by parts: I²R₁ + I²R₂ = (0.3429 mA)² × 35 kΩ = 4.114 mW. A 50 mV drop on a 2 m lead at 0.5 A is 1 % of a 5 V supply: enough to matter for an analogue sensor reference.

**Common mistakes.** Using the diameter as the radius (a factor 4 in R); mixing kΩ and Ω in the current.

## 2. Derivations

(a) Conductances g₁ = 1/R₁, …, g_g = 1/R_g. Unknowns (V₁, V₂, V₃, I), where I is the current flowing from node 1 into the source's positive terminal; it leaves node 1, so it enters node 1's KCL row (currents leaving = 0) with coefficient +1:

```
[ g1+g3     -g1          -g3          1 ] [V1]   [ 0  ]
[ -g1    g1+g2+gg        -gg          0 ] [V2] = [ 0  ]
[ -g3       -gg       g3+gx+gg        0 ] [V3]   [ 0  ]
[  1         0            0           0 ] [I ]   [ Vs ]
```

Symmetric: each resistor stamps a symmetric 2 × 2 block, and the source contributes the same +1 to row 1/column 4 and row 4/column 1. (The code reports the source current as −I, the current delivered to the circuit.)

(b) KCL at node 2: g₁(V₂ − V₁) + g₂V₂ + g_g(V₂ − V₃) = 0; at node 3: g₃(V₃ − V₁) + g_xV₃ + g_g(V₃ − V₂) = 0. If I_g = g_g(V₂ − V₃) = 0 then V₂ = V₁ R₂/(R₁ + R₂) and V₃ = V₁ R_x/(R₃ + R_x) (each node is then an unloaded divider), and V₂ = V₃ requires R₂/(R₁ + R₂) = R_x/(R₃ + R_x), i.e. **R₁/R₂ = R₃/R_x**. Conversely, if the ratios are equal, V₂ = V₃ = the common divider value satisfies both node equations with the g_g terms zero, and since the solution is unique (the matrix is non-singular), I_g = 0 for every R_g.

(c) Loop a: V_s → R₁ → R₂ → ground; loop b: R₁ → R_g → R₃ (the upper triangle); loop c: R₂ → R_x → R_g (the lower triangle). KVL:

(R₁ + R₂)i_a − R₁i_b − R₂i_c = V_s
−R₁i_a + (R₁ + R_g + R₃)i_b − R_g i_c = 0
−R₂i_a − R_g i_b + (R₂ + R_x + R_g)i_c = 0

R_g is shared by loops b and c, so it carries **i_c − i_b** (from node 2 to node 3 with these orientations). For the lesson's values this is −56.18 μA, matching the nodal solution to 10⁻¹⁰ (a test).

(d) With integral: I = (V_s/R)e^{−t/RC}, ∫₀^∞ I²R dt = (V_s²/R)(RC/2) = **½CV_s²**. Without: the charge that passes through the source is Q = CV_s, all at the source's fixed voltage V_s, so the source does work QV_s = CV_s². The capacitor ends with ½CV_s². Energy conservation puts the difference, **½CV_s²**, in the resistor, for any R (indeed for any resistive path, linear or not).

**Common mistakes.** In (a), adding the source as a conductance; in (b), proving only one direction; in (d), assuming the source's energy is ½CV_s² too.

## 3. A charged spacecraft ($s = 15$: V = −2.5 kV)

**Assumptions.** Isolated conducting sphere in vacuum (the plasma's screening is ignored), potential zero at infinity.

(a) C = 4πε₀R = **111.3 pF**; Q = CV = **−2.78 × 10⁻⁷ C**; **1.74 × 10¹² excess electrons**.

(b) Outside the sphere it acts as a point charge at its centre. Surface: |E| = k|Q|/R² = |V|/R = **2500 V/m**. At 2 m: **625 V/m** (a quarter), V = kQ/2 = **−1250 V** (a half). `electric_field([Q], [[0,0,0]], [2,0,0])` returns 625.0 V/m and `electric_potential` −1250.0 V.

(c) ½QV = **0.35 mJ**. Small, but an arc delivers it in nanoseconds: peak currents of amperes and fast voltage transients couple into nearby cables and can upset or destroy a sensitive input that is rated for microjoules, and arcs on solar arrays can trigger sustained arcs fed by the array itself. Spacecraft design for charging makes all surfaces conductive and bonded to a common ground so differential charging cannot build up.

**Check.** |E| = |V|/R at the surface of a sphere: 2500 V / 1 m. Units: C = F × V.

**Common mistakes.** Using 1/(4πε₀) as the capacitance; forgetting the sign of Q.

## 4. Sizing the bus harness ($s = 15$: loads 115, 60, 40 W)

(a) P = 215 W, I = 7.679 A. For 3 % (0.84 V over supply and return): copper **d = 1.769 mm**, mass **0.352 kg**; aluminium **d = 2.221 mm**, mass **0.167 kg**.

(b) At the copper diameter the solver gives V_load = **27.184 V**, drop **2.91 %**, current **7.455 A**, load power **202.7 W**, harness loss **6.08 W**. The drop is below 3 % because fixed-resistance loads draw less current at the lower voltage (7.455 A instead of 7.679 A), and the sizing formula used the nominal current.

(c) Constant-power loads: iteration converges to V_load = **27.133 V**, drop **3.10 %**, current **7.924 A**. Constant-power loads draw more current as the voltage falls, so they are the **conservative** model for sizing: the wire sized at the nominal current misses its 3 % budget by 0.1 %.

(d) Aluminium is half the mass for the same drop (its conductance per unit mass is about twice copper's), at a 26 % larger diameter. Before committing: connectors and crimps rated for aluminium (oxide layers and creep under clamping), galvanic corrosion at copper–aluminium junctions, the hot-case resistance (+16 % at 60 °C), and the constant-power case from (c).

**Common mistakes.** Using one wire's length instead of supply plus return; sizing at the load voltage instead of the bus voltage; forgetting the mass of both conductors.

## 5. Trusting the multimeter

(a) True midpoint **1.6500 V**. 10 MΩ meter: **1.5714 V** (4.8 % low). 1 MΩ meter: **1.1000 V** (33 % low).

(b) The divider seen from the midpoint is a source V_th = 1.65 V behind R_th = 0.5 MΩ. The meter forms a divider with it: V_read = V_th R_m/(R_m + R_th), so **V_true = V_read(1 + R_th/R_m)** = 1.5714 × 1.05 = **1.6500 V**.

(c) 0 Ω: **50.000 mA**; 1 Ω: **49.505 mA**; 10 Ω: **45.455 mA**.

(d) To first order the error is R_th/R_m for a voltmeter and R_a/R_circuit for an ammeter: for 0.1 % accuracy the meter's resistance must be at least 1000 times the circuit's Thévenin resistance (voltmeter), or at most a thousandth of the loop resistance (ammeter). A 10 MΩ meter is therefore fine for circuits below about 10 kΩ and needs correcting above about 100 kΩ.

**Common mistakes.** Using the full divider resistance (2 MΩ) instead of the Thévenin resistance (0.5 MΩ) at the midpoint.

## 6. Inrush into a capacitor bank

(a) τ = 10 Ω × 1000 μF = **10 ms**; peak current **2.8 A**; 99 % at τ ln 100 = **46.1 ms**.

(b) Stored ½CV² = ½ × 10⁻³ × 28² = **0.392 J**; dissipated **0.392 J**. `rc_energy_accounting(28, 10, 1e-3, 0.25)` returns 0.7840 J from the source, 0.3920 J stored and 0.3920 J dissipated.

(c) Maximum difference between the closed form and `solve_ivp` (DOP853, rtol 10⁻¹¹) over 5τ: **1.8 × 10⁻¹⁰ V**.

(d) With 0.1 Ω: peak **280 A**, τ = **0.1 ms**. The energy lost is the same 0.392 J (the half-energy rule), but it is delivered in a 0.1 ms pulse of hundreds of amperes, which pulls the bus voltage down for every other load, trips current protection, and stresses switches, connectors and the capacitors themselves. A limiter does not save energy; it controls the peak.

**Common mistakes.** Thinking the limiter "wastes" energy that would otherwise be saved; quoting 5τ as "fully charged" when the specification is 99 % (4.6τ).

## 7. Judgment

(a) **Fixed resistance.** The harness heats in sunlight; copper's resistance rises 0.39 % per kelvin, so at 60 °C the drop at the 2.19 mm size grows from 1.96 % to 2.26 %. Size the harness at the hot-case temperature.

(b) **Infinite meter resistance.** With R_m comparable to R_th the reading is about half the true value; use a meter or data-acquisition input with a much higher input resistance (or a buffer amplifier), or correct with V_true = V_read(1 + R_th/R_m).

(c) **Ideal wires.** A corroded crimp adds an unmodelled series resistance, perhaps tens of milliohms or more, which shows as an excess drop at a known current and as local heating. Measure the drop across each connector (a four-wire measurement) and replace it.

(d) **Constant-resistance loads.** Constant-power loads draw more current as the voltage sags, which sags it further; with a long thin harness the bus can collapse or oscillate, which a linear DC solve cannot show. Stagger the switch-on, use thicker wire, or regulate closer to the loads.

(e) **A reference node.** Without ground every node voltage is defined only up to a constant: the MNA matrix is singular. A solver that returns numbers anyway is hiding the error; `solve_dc` checks the rank and raises.

## Check-yourself answers

1. F = kq₁q₂/r². Halving the distance quadruples the force; changing one sign turns repulsion into attraction (same magnitude).
2. E = −∇V. V in volts (J/C); E in volts per metre (equivalently N/C).
3. Doubling the diameter divides R by 4 (area ∝ d²); doubling the length doubles R.
4. KCL: currents out of a node sum to zero (conservation of charge). KVL: voltages around a loop sum to zero (the electrostatic field is conservative, so the potential is single-valued).
5. A floating node's potential can be shifted by any constant without changing any current, so the equations have a one-parameter family of solutions: the matrix has a null vector and is singular.
6. The true value is 1.65 V; R_th = 0.5 MΩ, so the meter reads 1.65 × 10/10.5 = 1.571 V.
7. τ = 4.7 kΩ × 47 μF = 0.221 s; 1 − e⁻¹ = 63.2 % after one time constant.
8. Half, independent of R.
