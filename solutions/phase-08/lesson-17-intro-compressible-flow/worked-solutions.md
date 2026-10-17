# P08.17  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Air is a calorically perfect gas, γ = 1.4, R = 287.05287 J/(kg·K).

## 1. An airliner in cruise (s = 15: V = 245 m/s)

**Assumptions.** ISA at a geopotential altitude of 10 668 m (the tropopause is at 11 000 m, so this is still the gradient layer); isentropic stagnation in the pitot.

(a) **T = 218.808 K, p = 23 842.3 Pa, ρ = 0.37960 kg/m³, a = 296.535 m/s, M = 0.8262.**

(b) T₀ = T(1 + 0.2M²) = **248.68 K**; p₀ = p × 1.5650 = **37 314.1 Pa**; q_c = **13 471.8 Pa**, against ½ρV² = **11 392.6 Pa**: the impact pressure is 18 % larger than the dynamic pressure.

(c) √(2q_c/ρ) = **266.42 m/s**, too high by **8.74 %** (M²/8 alone predicts 8.5 %; the higher terms add the rest). The density at the pitot's stagnation point is **1.377** times the static density: a 38 % compression. An air-data computer that ignored compressibility would overstate the airliner's speed by 21 m/s.

**Independent check.** `subsonic_pitot_airspeed(13471.8, 23842.3, 218.808)` returns 245.0 m/s and M = 0.8262.

**Common mistakes.** Using the sea-level speed of sound; confusing q_c with ½ρV².

## 2. Derivations

(a) Steady adiabatic flow with no work: h + V²/2 = h₀. With h = c_pT, c_p = γR/(γ − 1) and a² = γRT: T₀/T = 1 + V²/(2c_pT) = **1 + (γ − 1)M²/2**. That needs only **adiabatic** flow (no heat added, no work), losses allowed. If the deceleration is also reversible (isentropic), p/ρ^γ is constant and, with p = ρRT, p ∝ T^{γ/(γ−1)} and ρ ∝ T^{1/(γ−1)}: **p₀/p = (T₀/T)^{γ/(γ−1)}, ρ₀/ρ = (T₀/T)^{1/(γ−1)}**. Behind a shock T₀ is unchanged but p₀ is lower.

(b) With x = (γ − 1)M²/2 and n = γ/(γ − 1), (1 + x)ⁿ = 1 + nx + n(n − 1)x²/2 + n(n − 1)(n − 2)x³/6 + …. Here nx = γM²/2, n(n − 1)x²/2 = γM⁴/8, n(n − 1)(n − 2)x³/6 = γ(2 − γ)M⁶/48. So p₀/p − 1 = (γM²/2)(1 + M²/4 + (2 − γ)M⁴/24 + …), and with ½ρV² = ½γpM²: **p₀ − p = ½ρV²(1 + M²/4 + (2 − γ)M⁴/24 + …)**. The Bernoulli speed √(2q_c/ρ) = V(1 + M²/4 + …)^{1/2} ≈ **V(1 + M²/8)**.

(c) A thin steady control volume across the shock (no area change, no wall forces): ρ₁u₁ = ρ₂u₂; p₁ + ρ₁u₁² = p₂ + ρ₂u₂²; c_pT₁ + u₁²/2 = c_pT₂ + u₂²/2. With ρu² = γpM², the momentum balance reads p₁(1 + γM₁²) = p₂(1 + γM₂²). Using mass and energy to express M₂ in terms of M₁, M₂² = (1 + (γ − 1)M₁²/2)/(γM₁² − (γ − 1)/2), and substituting: **p₂/p₁ = 1 + 2γ(M₁² − 1)/(γ + 1)**. Behind the shock the flow is subsonic and the probe brings it to rest isentropically, so p₀₂/p₁ = (p₂/p₁)(1 + (γ − 1)M₂²/2)^{γ/(γ−1)}. Substituting M₂ and simplifying gives **Rayleigh's formula** p₀₂/p₁ = [(γ + 1)²M₁²/(4γM₁² − 2(γ − 1))]^{γ/(γ−1)} (1 − γ + 2γM₁²)/(γ + 1). The test `test_rayleigh_pitot_is_shock_then_isentropic_and_continuous_at_mach_one` checks the identity numerically.

**Common mistakes.** Applying the isentropic p₀/p across the shock; using c_v instead of c_p in the energy equation.

## 3. A high-speed wind tunnel

**Assumptions.** Isentropic expansion from the chamber; perfect gas.

(a) T = 290/1.128 = **257.09 K**; p = 200 000/1.5243 = **131.20 kPa**; ρ = **1.7779 kg/m³**; V = 0.8 × 321.43 = **257.15 m/s**; ρV = **457.2 kg/(m²·s)**.

(b) p₀ − p = 68.80 kPa. With the static density: √(2 × 68 796/1.7779) = **278.2 m/s, +8.2 %**. With the chamber density ρ₀ = 2.4025 kg/m³: **239.3 m/s, −6.9 %**.

(c) Bernoulli would be exact with ρ* = 2(p₀ − p)/V² = **2.081 kg/m³**, between the two. The density varies along the streamline by 26 %, and the compressible relation is the exact integral of dp/ρ over that variation. No single constant density reproduces it at all Mach numbers. At M = 0.8 the tunnel must be calibrated with the compressible formula.

**Common mistakes.** Using the chamber temperature to compute the test-section speed of sound.

## 4. Pitot readings on a supersonic rocket (s = 15: h = 10 km)

**Assumptions.** A normal shock stands ahead of the probe's opening; isentropic compression behind it; ISA at 10 km: T = 223.15 K, p = 26 436.2 Pa, a = 299.46 m/s.

(a)–(c)

| M | p₀₂/p₁ | reading (kPa) | `mach_from_pitot` | subsonic inversion | Bernoulli V (m/s) vs true | p₀₂/p₀₁ |
|---|---|---|---|---|---|---|
| 1.5 | 3.4133 | 90.23 | 1.500000 | 1.4494 (−3.4 %) | 556.0 vs 449.2 (+24 %) | 0.9298 |
| 2.0 | 5.6404 | 149.11 | 2.000000 | 1.7879 (−10.6 %) | 771.0 vs 598.9 (+29 %) | 0.7209 |
| 2.5 | 8.5261 | 225.40 | 2.500000 | 2.0551 (−17.8 %) | 981.9 vs 748.7 (+31 %) | 0.4990 |
| 3.0 | 12.0610 | 318.85 | 3.000000 | 2.2769 (−24.1 %) | 1190.4 vs 898.4 (+32.5 %) | 0.3283 |

The subsonic inversion assumes the probe recovers the free-stream total pressure p₀₁; behind the shock it recovers only p₀₂ = (p₀₂/p₀₁)p₀₁, so it reads low, and the loss grows fast with M (two thirds of the total pressure is lost at M = 3). Incompressible Bernoulli over-reads by 24–32.5 %, because it ignores compressibility entirely.

(d) T₀ = **624.8 K** at M = 3 (from 223 K static): the probe tip runs at about 350 °C. It must be metallic (or ceramic) and the pressure transducer sits behind a tube, away from the tip, which in turn adds a pneumatic lag to the reading.

**Common mistakes.** Using the isentropic p₀/p above M = 1; forgetting that the Rayleigh formula already includes the isentropic part behind the shock.

## 5. Calibrated and equivalent airspeed

(a) a_SL = 340.294 m/s, ρ_SL = 1.225 kg/m³.

| condition | M | q_c (Pa) | CAS | EAS |
|---|---|---|---|---|
| light aeroplane, 3000 m, 60 m/s | 0.1826 | 1650.1 | 51.754 m/s (100.60 kt) | 51.689 m/s (100.47 kt) |
| airliner, FL350, 245 m/s | 0.8262 | 13 471.8 | 144.994 m/s (281.85 kt) | 136.383 m/s (265.11 kt) |

(b) CAS − EAS = **0.066 m/s** for the light aeroplane and **8.61 m/s (16.7 kt)** for the airliner. The CAS formula assumes sea-level pressure; where p is far below p_SL, the same q_c corresponds to a higher Mach number and a larger compressibility effect, which CAS does not remove and EAS does. At 3000 m and M = 0.18 the effect is a tenth of a per cent.

(c) √(2q_c/ρ_SL) = 51.904 m/s (light aeroplane, 0.3 % above CAS) and 148.306 m/s (airliner, 2.3 % above CAS). EAS belongs to aerodynamic loads and stall, because ½ρV² = ½ρ_SL EAS²: lift and loads depend on the dynamic pressure. CAS is what the instrument shows and is close to EAS at low altitude and speed, which is why limiting speeds are published as CAS for light aircraft. At high altitude the Mach number, not the dynamic pressure, often sets the limit.

**Common mistakes.** Using the local p in the CAS formula; quoting knots without saying so.

## 6. Is the flow compressible?

(a)

| M | ρ₀/ρ − 1 | Bernoulli airspeed error |
|---|---|---|
| 0.1 | 0.50 % | 0.125 % |
| 0.3 | 4.56 % | 1.129 % |
| 0.5 | 12.97 % | 3.154 % |

"Incompressible" means the density change is small compared with the accuracy needed. At M = 0.3 the density changes by 4.6 % at stagnation but the airspeed (a square root, with partial cancellation) is wrong by only 1.1 %. The density criterion is the stricter one.

(b) (i) Light aeroplane at M = 0.18: error 0.42 %, within the 1 % wanted; the incompressible formula would do, though instruments use the compressible one anyway. (ii) Cabin duct at 10 m/s (M ≈ 0.03): incompressible, the error is 10⁻⁴. (iii) Airliner at M = 0.83: error 8.7 %; compressible relations are essential. (iv) Tunnel at M = 0.8: compressible; errors of 7–8 % either way.

(c) (1) A probe misaligned with the flow (a large yaw or angle of attack) does not bring the flow to rest at its opening and reads below p₀. (2) At a very low probe Reynolds number (a tiny probe in slow flow) viscous effects make the stagnation pressure differ from the isentropic value. Other situations: a probe inside a separated wake or boundary layer, where the "free stream" it sees is not the one wanted; unsteady flow faster than the tube's response.

**Common mistakes.** Treating M < 0.3 as a law; it is a choice of tolerance.

## Check-yourself answers

1. a = √(γRT): 340.29 m/s at 288.15 K, 328.58 m/s at 3000 m, 295.07 m/s at 216.65 K. It depends only on temperature for a perfect gas.
2. T₀/T = 1 + (γ − 1)M²/2 (adiabatic); p₀/p = (T₀/T)^{γ/(γ−1)} and ρ₀/ρ = (T₀/T)^{1/(γ−1)} (isentropic).
3. At M = 0.18 the error is 0.42 %; it reaches 1 % at M = 0.282 (93 m/s at 3000 m) and 5 % at M = 0.628. It grows as M²/8.
4. Because p₀ − p = ½ρV²(1 + M²/4 + …) → ½ρV² as M → 0, which is Bernoulli's result for constant density.
5. A normal shock ahead of the probe raises the static pressure and loses total pressure; the probe then compresses subsonic air isentropically. Rayleigh's formula combines both.
6. At M = 1 both give p₀/p = 1.8929: the shock has zero strength there, so the two formulas join.
7. Mass, momentum and energy are conserved (ρu, p + ρu², h + u²/2 equal on both sides), but entropy rises, so the total pressure falls (p₀₂/p₀₁ = 0.7209 at M = 2).
8. The total temperature: T₀ = T + V²/(2c_p), 1.8 K above static for the light aeroplane, 173 K for the rocket at M = 2 at 11 km.
