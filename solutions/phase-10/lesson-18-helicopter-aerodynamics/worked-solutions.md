# P10.18  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results come from `solution.py` and a probe script that calls it, given to the figures shown. κ, C_d0, a, θ_tw, f and K are the lesson's stated assumptions; nothing here is a measurement.

## 1. Hover high and heavy (s = 15: 1150 kg)

**Assumptions.** Reference rotor, κ = 1.15, C_d0 = 0.010, ISA densities.

| | sea level (ρ = 1.225) | 2000 m (ρ = 1.0066) |
|---|---|---|
| v_h | 7.657 m/s | 8.447 m/s |
| ideal power | 86.38 kW | 95.29 kW |
| hover power (5) | 134.79 kW | 138.72 kW |
| FM | 0.641 | 0.687 |
| C_T/σ | 0.0835 | 0.1017 |

(b) The induced power grows as 1/√ρ (T^{3/2}/√(2ρA)): +10.3 %. The profile power, ρA(ΩR)³σC_d0/8, falls in proportion to ρ: −18 %. The total rises by only 2.9 % here because the profile share is large; the figure of merit rises because the ideal power grew.

(c) The mean blade lift coefficient is 6C_T/σ: 0.501 at sea level and 0.610 at 2000 m. The blades work closer to stall, so the margin for manoeuvre (and for the retreating blade in forward flight) shrinks; at some altitude the blades, not the engine, limit the helicopter.

**Independent check.** v_h² ratio = 1.225/1.0066 = 1.217 = (8.447/7.657)² ✓.

**Common mistakes.** Scaling the whole power by 1/√ρ; forgetting that C_T itself changes with ρ.

## 2. Derivations

(a) A blade element at radius y = rR has tangential speed Ωy and inflow λΩR, so its angle of attack is θ − φ with φ ≈ λΩR/(Ωy) = λ/r. Its lift per unit span is ½ρ(ΩrR)²c a(θ − λ/r). For B blades, dT = B½ρΩ²R²r²c a(θ − λ/r)R dr. Dividing by ρπR²(ΩR)² and writing σ = Bc/(πR): **dC_T = (σa/2)(θr² − λr) dr**. With θ = θ₀ + θ_tw r and constant λ: ∫₀¹(θ₀r² + θ_tw r³ − λr) dr = θ₀/3 + θ_tw/4 − λ/2, which is (3).

(b) The induced power of an element is its thrust times the induced velocity: dC_Pi = λ dC_T. With λ uniform, C_Pi = λC_T, and momentum theory gives λ = √(C_T/2), so **C_Pi = C_T^{3/2}/√2**, the ideal power (1) in coefficient form.

(c) (r + μ sin ψ)³ = r³ + 3r²μ sin ψ + 3rμ² sin²ψ + μ³ sin³ψ. Averaging over ψ, the odd powers of sin ψ vanish and sin²ψ averages to ½: r³ + (3/2)rμ². Integrating over 0 ≤ r ≤ 1: 1/4 + (3/4)μ² = (1/4)(1 + 3μ²), so **K = 3**. Counting the radial flow replaces u_T³ by (u_T² + u_R²)^{3/2} ≥ |u_T|³: the drag acts on the whole in-plane velocity, which is larger, so more power is absorbed. Expanding (r² + 2rμ sin ψ + μ²)^{3/2} to second order in μ and averaging over ψ gives r³ + (9/4)rμ², so K → 4 × (9/8) = 4.5 as μ → 0; the disk integral gives 4.51 at μ = 0.05, rising to 4.84 at μ = 0.4.

(d) At μ = 0, (6) reads λ = C_T/(2λ), so λ = √(C_T/2) = λ_h. For μ ≫ λ, √(μ² + λ²) ≈ μ and λ_i → C_T/(2μ). Multiplying by ΩR, with C_T = T/(ρA(ΩR)²) and μ ≈ V/(ΩR): **v = λ_iΩR = T/(2ρAV)**, the induced velocity of a wing of span 2R whose lift is carried by the air passing through the disk's area.

**Common mistakes.** Using φ = λ (not λ/r); averaging sin³ψ to a non-zero value; forgetting the θ_tw/4 term.

## 3. A quieter rotor

**Assumptions.** Same mass, chord, twist, κ, C_d0, f; only ΩR changes.

| ΩR (m/s) | P hover (kW) | induced | profile | FM | C_T/σ | 6C_T/σ | θ_0.75 | P_min forward (kW) |
|---|---|---|---|---|---|---|---|---|
| 210 | 128.38 | 92.93 | 35.45 | 0.629 | 0.0799 | 0.479 | 7.86° | 74.86 |
| 195 | 121.32 | 92.93 | 28.39 | 0.666 | 0.0927 | 0.556 | 8.86° | 67.65 |
| 180 | 115.26 | 92.93 | 22.33 | 0.701 | 0.1088 | 0.653 | 10.10° | 61.44 |

(b) The induced power is T^{3/2}/√(2ρA)·κ: it depends on thrust and disk area, not on how fast the blades turn. The profile power is σC_d0/8 × ρA(ΩR)³: it falls as the cube of the tip speed (180/210 cubed = 0.63).

(c) The blades must make the same thrust at lower speed, so they work at a higher lift coefficient (mean 0.653 at 180 m/s against 0.479): the stall margin shrinks. In forward flight it is worse: μ at 60 m/s is **0.286, 0.308, 0.333** for 210, 195, 180 m/s; the retreating tip moves at (1 − μ)ΩR = 150, 135, 120 m/s and must carry its share of the lift at even higher cl, so retreating-blade stall limits the speed and the gust and manoeuvre margins. The minimum tip speed is set by stall (and by the rotor's stored kinetic energy for autorotation entry), not by the power.

(d) Listed in the table: 74.86, 67.65, 61.44 kW (the minimum-power speed stays near 23.7 m/s).

**Independent check.** The induced power is identical to 10⁻¹² kW in all three rows, as (b) predicts.

**Common mistakes.** Concluding that the slowest rotor is best; this model has no stall, so it cannot see the penalty.

## 4. The power curve of your helicopter (s = 15: f = 1.10 m²)

(a) Minimum power **75.64 kW at 23.15 m/s**; maximum range at **36.91 m/s**, P = 93.30 kW. Least autorotative descent **7.01 m/s** (energy estimate) at 23.15 m/s; best glide at 36.9 m/s with descent 8.65 m/s, **glide ratio 4.27**.

(b) P_av = 0.9 × 160 = 144 kW. Rate of climb (P_av − P)/W is largest at V_mp: **6.33 m/s at 23.15 m/s**. Check: `forward_flight_power(23.15, rot, climb=6.334)` returns 144.04 kW, the available power ✓ (the 0.04 kW is the grid step in V and the small difference between T and W). Maximum level speed where P = P_av: **50.19 m/s (98 kn)**.

(c) The fuel per metre is proportional to P/V, the slope of the line from the origin to a point on the curve; that slope is least where the line just touches the curve.

**Independent check.** At V_max the parasite power ½ρV³f = 85.2 kW is the largest part, consistent with the dominance of f at high speed.

**Common mistakes.** Using the hover excess power for the best climb; reading V_max where P = 160 kW (forgetting the transmission).

## 5. The profile-power factor

(a) With radial flow (400-point grid): K = **4.507, 4.529, 4.601, 4.707, 4.842** at μ = 0.05, 0.1, 0.2, 0.3, 0.4. Without radial flow (lesson's 200-point values): 3.002 at 0.1, 3.015 at 0.2, 3.034 at 0.3, 3.060 at 0.4.

(b) Inside the reverse-flow circle (r < −μ sin ψ on the retreating side) u_T < 0. The analytic K = 3 integrates u_T³, which is negative there; the physical power is ∝ |u_T|³, positive. The difference grows with the size of the circle, roughly as μ⁴, hence 3.002 → 3.060.

(c)

| K | V_mp (m/s) | P_min (kW) | V_mr (m/s) | P at V_mr (kW) | P at 60 m/s (kW) |
|---|---|---|---|---|---|
| 3 | 24.00 | 74.11 | 38.16 | 91.27 | 188.13 |
| 4.65 | 23.66 | 74.86 | 37.88 | 92.49 | 192.71 |
| 5 | 23.60 | 75.02 | 37.82 | 92.74 | 193.69 |

K matters at high speed, where μ² is large (2.4 % at 60 m/s between 3 and 4.65), and hardly at the bucket (1 %).

**Common mistakes.** Treating K as a fitted constant without knowing it is μ-dependent.

## 6. Where does the model fail?

(a)

| input | hover P (kW) | FM | least V_d (m/s) |
|---|---|---|---|
| κ = 1.10 | 124.34 | 0.650 | 6.82 |
| κ = 1.15 | 128.38 | 0.629 | 6.94 |
| κ = 1.25 | 136.46 | 0.592 | 7.17 |
| C_d0 = 0.008 | 121.29 | 0.666 | 6.24 |
| C_d0 = 0.010 | 128.38 | 0.629 | 6.94 |
| C_d0 = 0.012 | 135.47 | 0.596 | 7.63 |

For hover both matter about equally over these ranges; for the autorotative descent rate (and the bucket) C_d0 matters most, because the profile power is about half the power at V_mp. A flight test should first measure power against speed near the bucket, which pins C_d0, and then hover power, which then pins κ.

(b)

| V (m/s) | μ | reverse flow r/R < | retreating tip (m/s) | advancing tip Mach |
|---|---|---|---|---|
| 40 | 0.190 | 0.190 | 170.2 | 0.734 |
| 60 | 0.280 | 0.280 | 151.2 | 0.790 |
| 70 | 0.321 | 0.321 | 142.6 | 0.815 |

At 70 m/s the model leaves out retreating-blade stall (the retreating side must produce lift at a much lower speed, hence a high cl) and advancing-tip compressibility (drag rise near Mach 0.8 on a 12 % thick section), both of which raise the power and limit speed.

(c) Momentum theory assumes a slipstream that leaves the disk in one direction. In a slow vertical descent the air moving up relative to the rotor meets the downward induced flow: the rotor descends into its own wake, the tip vortices pile up into a ring (the vortex-ring state), and the flow is unsteady and recirculating. Neither Glauert's equation nor the actuator disk applies; empirical descent data are needed.

(d) (1) The tail rotor, transmission and accessories still absorb power from the main rotor in autorotation. (2) The descent changes the inflow (the air comes up through the disk), and the rotor works in a less efficient state than the level-flight power curve assumes.

**Common mistakes.** Reporting the energy estimate as the descent rate; using the model at μ > 0.3 without comment.

## Check-yourself answers

1. C_T = T/(ρA(ΩR)²), C_P = P/(ρA(ΩR)³), λ = (V sin α + v)/(ΩR), μ = V cos α/(ΩR), σ = Bc/(πR). In hover λ = √(C_T/2).
2. FM = ideal power / actual power. The profile power is fixed while the induced power grows as T^{3/2}, so at a higher thrust the ideal part is a larger share of the total.
3. The local dynamic pressure grows as r², so ∫θr² dr = θ/3: the thrust is weighted to the outer blade, and the average of r² over the span is 1/3 (θ/2 would be the lift of a wing with uniform dynamic pressure, which a rotor is not).
4. λ = μ tan α + C_T/(2√(μ² + λ²)); λ_i → λ_h = √(C_T/2) as μ → 0, and λ_i → C_T/(2μ) (v = T/(2ρAV)) when μ ≫ λ.
5. Induced (largest in hover and at low speed, below about 20 m/s here), profile (roughly constant, rising slowly with μ²; dominant near the bucket), parasite (∝ V³, dominant above about 40 m/s).
6. Averaging (r + μ sin ψ)³ over the disk gives 1 + 3μ² (times the hover value); including the drag of the radial flow, the in-plane speed is larger and the integral gives about 4.5–4.8 over μ = 0.05–0.4, represented by the constant 4.65.
7. Draw the tangent from the origin to the power curve; the tangent point minimises P/V, the energy per metre. Beyond V_mp the power rises slowly at first, so the extra speed costs proportionally less power than it gains in distance, until the tangency.
8. It leaves out the tail-rotor and transmission power the rotor still drives, and it uses level-flight inflow rather than the descending-flow state.
