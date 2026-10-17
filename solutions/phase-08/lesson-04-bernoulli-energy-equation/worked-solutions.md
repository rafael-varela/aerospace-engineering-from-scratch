# P08.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Propellant data are typical values.

## 1. A pitot tube on a water manometer (s = 15: Δh = 45 mm)

**Assumptions.** Steady, incompressible (M ≈ 0.08), inviscid along the stagnation streamline; perfect static holes; air at 1.204 kg/m³.

(a) Δp = (998.2 − 1.204) × 9.80665 × 0.045 = **439.97 Pa**; V = √(2 × 439.97/1.204) = **27.03 m/s** (52.6 kt).

(b) V = √(2Δp/ρ) ∝ Δh^{1/2}, so ln V = ½ ln Δh + const and **δV/V = ½ δ(Δh)/Δh**. With δ(Δh) = 1 mm: δV/V = 0.5/45 = **1.1 %**, δV = **±0.30 m/s**. At lower speeds it grows fast: Δh ∝ V², so the same 1 mm is 4.4 % at half the speed. That is why an inclined manometer or an electronic transducer is used at low speed (`P08.19`).

(c) The flow is steady and the density change is 0.2 % (M²/4 at M = 0.08), so the least well satisfied condition is inviscid flow very near the probe's nose, where a thin boundary layer forms. Along the stagnation streamline itself the losses are negligible for a probe Reynolds number of a few thousand (a 3 mm tube at 27 m/s: Re ≈ 5400); viscous corrections to pitot readings matter only below Re of order 100.

**Common mistakes.** Using ρ_water in the airspeed formula; forgetting the factor ½ in the uncertainty.

## 2. Bernoulli, the pitot tube and the energy equation

(a) An element of cross-section δA and length δs along a streamline. Forces along s: pressure, p δA − (p + (∂p/∂s)δs)δA, and the weight component −ρg δA δs (∂z/∂s). **No shear forces: inviscid.** Mass ρδAδs times the acceleration along s, which for a **steady** flow is V ∂V/∂s (the local part ∂V/∂t is zero). Dividing by δAδs: **∂p/∂s + ρV ∂V/∂s + ρg ∂z/∂s = 0**. With **ρ constant** this is ∂/∂s(p + ½ρV² + ρgz) = 0: **p + ½ρV² + ρgz is constant along the streamline**. **No shaft work** is implicit: nothing between the two points adds energy (a fan would add a pressure jump not in this force balance).

(b) From the free stream (p, V) to the stagnation point on the probe's nose (p₀, 0), at the same height: **p₀ − p = ½ρV²**; the static holes on the side sense p. For the compressible law, with M² = ρV²/(γp) and x = 0.2M²: (1 + x)^{3.5} − 1 = 3.5x + (3.5 × 2.5/2)x² + … = 0.7M² + 0.175M⁴ + …. So q_c = p(0.7M² + 0.175M⁴) = ½ρV²(1 + 0.25M² + …) using 0.7pM² = ½ρV²: **q_c = ½ρV²(1 + M²/4 + …)**.

(c) **EAS = TAS√(ρ/ρ₀)**, so ½ρ₀EAS² = ½ρTAS² = q. An incompressible indicator solves ½ρ₀V_ind² = (p₀ − p) = q, so V_ind = EAS. Lift, drag and stall depend on q (and on the angle of attack), so a wing stalls at the same EAS at every altitude: pilots fly the approach and the take-off by IAS ≈ CAS ≈ EAS and never need to compute TAS for it.

(d) Per unit weight, steady flow between sections 1 and 2: p₁/(ρg) + α₁V₁²/(2g) + z₁ + h_pump = p₂/(ρg) + α₂V₂²/(2g) + z₂ + h_turbine + h_L. With no shaft work (h_pump = h_turbine = 0), no heat transfer and no loss (h_L = 0), and uniform profiles (α = 1), this is Bernoulli. With a loss, **h_L = H₁ − H₂**, the drop in total head H = p/(ρg) + V²/(2g) + z: friction turns mechanical energy into heat that the equation does not recover.

**Check.** The lesson's code shows q_c/q − 1 = 0.00836 at M = 0.1826 against M²/4 = 0.00834.

**Common mistakes.** Applying Bernoulli across a fan or pump; dropping the steadiness condition.

## 3. A rocket turbopump (s = 15: ṁ = 27.5 kg/s)

**Assumptions.** Incompressible liquid, steady, uniform profiles (α = 1), losses folded into the efficiency.

(a) Q = 27.5/1141 = **0.02410 m³/s**. h_shaft = (50 − 3) × 10⁵/(1141 × 9.80665) + (20² − 3²)/(2 × 9.80665) + 0.5 = 420.04 + 19.94 + 0.5 = **440.48 m**. The pressure rise dominates (95 %).

(b) P = ρgQh/η = **169.7 kW**. The estimate ṁΔp/(ρη) = **161.8 kW** misses the kinetic-energy and height terms, 4.6 %.

(c) The 30 % (about 51 kW) is dissipated in the pump by friction, leakage and turbulence and appears as heat in the liquid: with a typical specific heat of 1700 J/(kg·K) for liquid oxygen, the temperature rise is about 51 000/(27.5 × 1700) = **1.1 K**. For a cryogenic liquid close to its boiling point any warming raises its vapour pressure and brings it closer to cavitation at the next low-pressure point.

(d) The turbine burns propellant (in a gas-generator cycle that propellant leaves with little thrust; in a staged-combustion cycle the turbine's pressure drop limits the chamber pressure), so every kilowatt of pump power is paid in specific impulse.

**Common mistakes.** Using gauge and absolute pressures inconsistently; forgetting to divide by the efficiency.

## 4. Your own pitot-static system (s = 15: h = 3500 m, TAS = 65 m/s)

(a) M = **0.1990**; q = **1823.6 Pa**; q_c = **1841.7 Pa**; TAS 65.00 m/s (**126.3 kt**); EAS **54.56 m/s (106.1 kt)**; CAS **54.66 m/s (106.2 kt)**; IAS with perfect ports = CAS.

(b) With c_p = −0.03 the static port senses 54.7 Pa less than the true static pressure: IAS = **55.46 m/s (107.8 kt)**, IAS − CAS = **+1.56 kt**, and the altimeter reads **6.5 m high**.

(c) For 100 kt IAS with that port: TAS = **60.31 m/s (117.2 kt)** by root finding. Treating 100 kt as CAS: **61.19 m/s (118.9 kt)**, 0.88 m/s (1.7 kt) too fast. The rule TAS = IAS/√σ with σ = 0.7047: **61.28 m/s (119.1 kt)**: it also ignores the compressibility difference between CAS and EAS, which is small here.

(d) For navigation, a 1.7 kt error in TAS is 1.4 % of the time between fixes, absorbed by wind uncertainty; it matters for fuel planning on long legs only in sum. For the stall, the stall speed is published as IAS for that aircraft's own position error, so the pilot is protected only if the published figure was measured with the same port; a 1.6 kt error near the stall is a real margin. Position error is calibrated in flight test for exactly this reason.

**Common mistakes.** Taking IAS as TAS at altitude; applying the port error with the wrong sign (a port in fast air reads low pressure, so IAS and altitude read high).

## 5. Where Bernoulli's airspeed fails, and a contraction (s = 15: A_in = 4.5 m², CR = 9)

(a) The ratio q_c/q depends only on M, so does the error of V = √(2q_c/ρ). It reaches **1 % at M = 0.2824** (TAS 96.1 m/s at sea level, 92.8 m/s at 3000 m) and **5 % at M = 0.6279** (213.7 m/s at sea level, 206.3 m/s at 3000 m). The light aeroplane at M = 0.18 is safe with Bernoulli to 0.4 %; a turboprop at M = 0.5 is not.

(b)

| n | max \|Euler − Bernoulli\| (Pa) | observed order |
|---|---|---|
| 25 | 9.86 | — |
| 50 | 2.54 | 1.95 |
| 100 | 0.641 | 1.99 |
| 200 | 0.161 | 2.00 |

The pressure drop across the contraction is **535.11 Pa**.

(c) Continuity V_in = V_t/CR; Bernoulli Δp = ½ρ(V_t² − V_in²) = **q_test(1 − 1/CR²)**: 541.8 × (1 − 1/81) = 535.11 Pa ✓. `venturi_flow_rate(535.11, 4.5, 0.5, 1.204)` returns **15.000 m³/s**, the true flow. A wall ring at each end of the contraction measures the speed without putting a probe in the test section, where it would disturb the model and its wake; it is calibrated once against a pitot in the empty tunnel, which also absorbs the small losses in the contraction.

**Common mistakes.** Thinking the Bernoulli error depends on altitude; inverting the contraction formula with CR instead of CR².

## 6. Where does Bernoulli fail?

(a) Yes: steady, inviscid outside the boundary layer, incompressible at M = 0.18 to 0.4 %, no work: the pitot formula.
(b) No: the propeller does work on the air (a pressure jump). Use Bernoulli separately upstream and downstream, joined by the momentum balance of the disc (Problem 3 of `P08.03`).
(c) No: viscosity matters inside the boundary layer; total pressure falls towards the wall. Use the boundary-layer equations (`P08.10`); note that the static pressure across a thin boundary layer is still the outer value, which is why a static port works.
(d) No: the separated wake is turbulent and dissipative, so total pressure is lost. Use the momentum balance with the measured wake profile (wake rake) to get the drag.
(e) No: at M = 0.8 the density change matters (Bernoulli's error is about 8 %, from (a)'s trend). Use the compressible relation (`P08.17`).
(f) Borderline: the reduced frequency is ωc/U = (2π/0.1 s)(1.5 m)/(60 m/s) ≈ 1.6, so the local acceleration ∂V/∂t is not negligible. Use the unsteady Bernoulli equation (with ∂φ/∂t) or an unsteady aerodynamic model; a pitot tube itself, being small, still reads the instantaneous local q.

**Static port.** The static pressure on the fuselage varies with position: above free stream at the nose (c_p → 1 at the stagnation point), below it where the flow speeds up round the shoulders, recovering towards zero along the cylindrical part. A port on the nose would read nearly the total pressure; a port on the side of the fuselage at a station where c_p ≈ 0 reads the free-stream static pressure, and placing ports on both sides and connecting them cancels the effect of sideslip.

## Check-yourself answers

1. Steady flow, inviscid (no shear stress along the streamline), incompressible, no shaft work or heat addition between the two points; and the two points on the same streamline (or anywhere in an irrotational flow).
2. Total (stagnation) pressure is the pressure a flow reaches when brought to rest without loss, p₀ = p + ½ρV² in incompressible flow; dynamic pressure is q = ½ρV².
3. TAS is the speed relative to the air; EAS is the sea-level speed with the same dynamic pressure, TAS√σ; CAS is what a perfect instrument calibrated to the standard (compressible) law reads; IAS is what the dial shows, including position and instrument errors.
4. Because a wing's lift depends on q (and α), so the stall happens at a fixed EAS ≈ CAS at low speed whatever the altitude.
5. It senses a pressure p + c_p q; with c_p < 0 it reads low, so the altimeter reads high (by c_p q/(ρg)) and the airspeed indicator, which sees q_c − c_p q, reads high.
6. Above M ≈ 0.3 the error of the incompressible formula exceeds about 1 % (the q_c/q ratio is 1 + M²/4 + …).
7. The head loss h_L, the drop in total head p/(ρg) + V²/(2g) + z between two sections, is mechanical energy turned into heat by friction; a pump or fan adds head h_pump, a turbine removes h_turbine.
8. It throws away the exit kinetic energy and loses energy in the screens, corners and diffuser; the fan supplies both. A larger diffuser exit lowers the discarded kinetic energy, so it lowers the power, which is why tunnels have long diffusers.
