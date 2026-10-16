# P02.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown.

## 1. Avionics shelf ($s = 15$: m = 3.5 kg)

**Assumptions.** Linear isolators and damper, one degree of freedom, SI units.

**Derivation.** (a) ω_n = √(8000/3.5) = 47.81 rad/s = 7.609 Hz. √(km) = √(28 000) = 167.3 kg/s, ζ = 40/(2 × 167.3) = 0.1195. ω_d = 47.81√(1 − 0.01429) = 47.47 rad/s. Roots −ζω_n ± iω_d = **−5.714 ± 47.47i s⁻¹**; underdamped. (b) With x₀ = 5 mm, v₀ = 0: x(t) = e^{−5.714t}[5 cos 47.47t + (5.714 × 5/47.47) sin 47.47t] mm = e^{−5.714t}[5 cos 47.47t + 0.602 sin 47.47t] mm. x(0.1) = **−0.243 mm**, x(0.5) = **+0.0149 mm**. (c) τ = 1/(ζω_n) = 1/5.714 = 0.175 s; δ = 2π × 0.1195/√(1 − 0.01429) = 0.756 (peaks shrink by e^{0.756} = 2.13 per cycle). (d) c_crit = 2√(km) = **334.7 N s/m**, 8.4 times the fitted damper.

**Check.** `free_response([0.1, 0.5], 3.5, 40, 8000, 0.005, 0)` → (−2.430 × 10⁻⁴, 1.488 × 10⁻⁵) m ✓; `characteristic_roots` → −5.7143 ± 47.466i ✓. Dimensions: c/(2m) = (kg/s)/kg = s⁻¹ ✓. Limiting case: c → 0 gives ω_d → ω_n and the roots ±iω_n.

**Common mistakes.** Dropping the sine term for release from rest (v₀ = 0 does not make it zero); using ω_n instead of ω_d inside the trigonometric functions; quoting the roots in Hz.

## 2. The three regimes and the Wronskian

**Derivation.** (a) x = e^{rt}: (mr² + cr + k)e^{rt} = 0 ⟹ mr² + cr + k = 0. With c = 2ζω_n m and k = ω_n²m: r² + 2ζω_n r + ω_n² = 0, r = −ζω_n ± √(ζ²ω_n² − ω_n²) = ω_n(−ζ ± √(ζ² − 1)).

(b) For ζ < 1, √(ζ² − 1) = i√(1 − ζ²), so r = −σ ± iω_d with σ = ζω_n, ω_d = ω_n√(1 − ζ²). e^{(−σ ± iω_d)t} = e^{−σt}(cos ω_d t ± i sin ω_d t). The real and imaginary parts, e^{−σt} cos ω_d t and e^{−σt} sin ω_d t, are real solutions (a real linear equation maps real and imaginary parts separately), and they are independent (their Wronskian at 0 is ω_d ≠ 0). x = e^{−σt}(A cos ω_d t + B sin ω_d t): x(0) = A = x₀; ẋ(0) = −σA + ω_d B = v₀ ⟹ B = (v₀ + σx₀)/ω_d.

(c) ζ = 1: x = u e^{−ω_n t}, ẋ = (u̇ − ω_n u)e^{−ω_n t}, ẍ = (ü − 2ω_n u̇ + ω_n²u)e^{−ω_n t}. Then ẍ + 2ω_n ẋ + ω_n² x = [ü − 2ω_n u̇ + ω_n²u + 2ω_n u̇ − 2ω_n²u + ω_n²u]e^{−ω_n t} = ü e^{−ω_n t} = 0 ⟹ ü = 0 ⟹ u = A + Bt. The second solution is t e^{−ω_n t}; with x(0), ẋ(0): A = x₀, B = v₀ + ω_n x₀.

(d) W = x₁ẋ₂ − ẋ₁x₂. Ẇ = ẋ₁ẋ₂ + x₁ẍ₂ − ẍ₁x₂ − ẋ₁ẋ₂ = x₁ẍ₂ − ẍ₁x₂. Substitute ẍ_i = −(c/m)ẋ_i − (k/m)x_i: Ẇ = x₁[−(c/m)ẋ₂ − (k/m)x₂] − [−(c/m)ẋ₁ − (k/m)x₁]x₂ = −(c/m)(x₁ẋ₂ − ẋ₁x₂) = −(c/m)W. Hence W = W(0)e^{−ct/m}.

(e) Write x₀ cos ω_d t + B sin ω_d t = R cos(ω_d t − φ). Peaks of e^{−σt}R cos(ω_d t − φ): their spacing is governed by the cosine (the envelope only tilts them very slightly; to first order in σ/ω_d the shift is the same for every peak), so successive positive peaks are T_d = 2π/ω_d apart and their ratio is e^{σT_d}. δ = σT_d = 2πζω_n/(ω_n√(1 − ζ²)) = 2πζ/√(1 − ζ²). Square: δ²(1 − ζ²) = 4π²ζ² ⟹ ζ² = δ²/(4π² + δ²) ⟹ ζ = δ/√(4π² + δ²).

**Check.** (e) at small ζ: δ ≈ 2πζ, the usual approximation; at ζ = 0.2, exact δ = 1.283 against 2πζ = 1.257, a 2 % difference the code resolves.

**Common mistakes.** In (c), substituting into the equation with c and k instead of the ζ-ω_n form and losing track of the cancellation; in (e), using the small-ζ formula ζ = δ/2π for ζ > 0.1.

## 3. The RLC stage ($s = 15$: R = 15 Ω)

**Derivation.** (a) L ↔ m, R ↔ c, 1/C ↔ k. ω_n = 1/√(LC) = 1/√(10⁻² × 4.7 × 10⁻⁶) = 4613 rad/s = **734.1 Hz**. ζ = R/(2√(L/C)) = 15/(2 × 46.13) = **0.1626**. (b) Underdamped. R_crit = 2√(L/C) = **92.25 Ω**. (c) At 15 Ω the charge rings at ω_d = 4551 rad/s (724 Hz) inside an envelope with time constant 2L/R = 1.33 ms, taking several cycles to die; at 92 Ω it decays smoothly in about 1 ms with no reversal. A filter designer wants neither the long ringing (it smears transients into the next sample) nor the sluggishness of heavy damping; the usual compromise is ζ ≈ 0.7, which here is R ≈ 65 Ω. (d) `natural_frequency(0.01, 1/4.7e-6)` = 4612.7; `damping_ratio(0.01, 15, 1/4.7e-6)` = 0.1626; `classify_damping` → 'underdamped'; `damper_for_ratio(0.01, 1/4.7e-6, 1.0)` = 92.25 ✓.

**Check.** Units: 1/√(H F) = 1/√(s²) = s⁻¹ ✓; √(L/C) = √(H/F) = Ω ✓, so ζ is dimensionless.

**Common mistakes.** Mapping C (not 1/C) onto k; forgetting the factor 2 in R_crit.

## 4. Your landing ($s = 15$: v₀ = 3.5 m/s)

**Result.** (a)

| ζ | stroke (mm) | peak force (kN) | bounce | 2 % settle (s) |
|---|---|---|---|---|
| 0.2 | 299.0 | 47.7 | yes | 2.241 |
| 0.3 | 263.9 | 45.5 | no | 1.435 |
| 0.4 | 235.6 | 45.5 | no | 1.152 |
| 0.5 | 212.4 | 48.5 | no | 0.908 |
| 0.7 | 176.7 | 65.7 | no | 0.681 |
| 1.0 | 140.2 | 93.9 | no | 0.616 |

(b) Minimum peak force 45.2 kN at ζ = 0.347; bounce threshold ζ = 0.276. (c) Recommend **ζ ≈ 0.4**: 0.12 above the bounce threshold (margin for a heavier aircraft), peak force 45.5 kN within 1 % of the minimum, stroke 236 mm, settled in 1.15 s. (d) `solve_ivp` against the closed form at ζ = 0.4: maximum difference 4.6 × 10⁻¹² m.

**Check.** The whole table scales sensibly from the lesson's v₀ = 3 m/s table: strokes up by roughly the ratio of velocities, peak forces up by 3.5/3 = 1.17 (the damper spike cv₀ exactly so).

**Common mistakes.** Measuring the settling band about zero instead of about x_s; recommending critical damping because it "settles fastest" without looking at the force column.

## 5. What the strut model cannot say

(a) The linear damper force at t = 0 is c v₀ with the spring still slack: 56 kN at ζ = 0.7, 80 kN at ζ = 1. An orifice damper gives F_d = c₂ẋ|ẋ|, zero at any instant when the velocity is zero: it is zero at t = 0 *and* at maximum stroke, and peaks in between. The real force history has no touchdown spike, so the "peak force minimum at ζ = 0.36" is an artefact of the linear model; with an orifice the designer tunes the orifice area (and a metering pin that varies it with stroke) to make the force nearly flat over the stroke, which is the most efficient way to absorb a given energy in a given travel.

(b) Negative compression means the strut is at its fully extended stop and the tyre has unloaded: the runway can push on the tyre but not pull. The aircraft is airborne again with the upward velocity it had at that instant, and the only thing left in contact is the tyre, a stiffer and almost undamped spring, which is why the rebound is quick and why bounced landings often lead to a second, harder touchdown.

(c) Stroke: remove the linear spring first; the gas spring stiffens strongly near full compression and bounds the stroke. Load: remove the linear damper first; it sets the shape of the force history and the spurious spike.

(d) The force–velocity plot from the drop test: a linear damper gives a straight line through the origin (after subtracting the spring force), an orifice gives a parabola. Equivalently, F(0) ≈ 0 at the instant of contact.

**Common mistakes.** Treating the bounce flag as a prediction of what happens after the bounce; thinking "more damping is always safer".

## 6. The first peak and the least load

**Derivation.** (a) ẋ = e^{−σt}[(−σx₀ + bω_d) cos ω_d t − (σb + x₀ω_d) sin ω_d t]. ẋ = 0 ⟹ tan(ω_d t*) = (−σx₀ + bω_d)/(σb + x₀ω_d); if the arctan is negative, add π/ω_d (the first zero is in the first half-period). With ω_n = 11.180, x₀ = −0.0785 m, v₀ = 3 m/s: ζ = 0.3: σ = 3.354, ω_d = 10.665, b = 0.2566, t* = **0.1465 s**; ζ = 0.7: σ = 7.826, ω_d = 7.984, b = 0.2988, t* = **0.1318 s**. `max_stroke` returns 0.14654 and 0.13179 s ✓ (strokes 235.1 and 157.1 mm).

(b) Newton for the mass: m Ẍ = mg − kX − cẊ = mg − F, so F = m(g − Ẍ). F is largest where Ẍ is most negative, i.e. where the mass decelerates hardest. For ζ = 0.3 (fine 10⁻⁴ s grid): F_max = 40.2 kN at t = 0.089 s, *before* the stroke peak at 0.1465 s, where F is only 35.3 kN; F(0) = 24.1 kN. For ζ = 1.0: F_max = 80.5 kN at t = 0, the damper spike; at the stroke peak (0.126 s) F is 18.7 kN.

(c) Two terms compete: the touchdown spike F(0) = cv₀ = 2ζ√(km)v₀, proportional to ζ and v₀, and the mid-stroke force, dominated by k x_max, which falls with ζ. For small ζ the mid-stroke term is the maximum and falls as ζ rises; for large ζ the spike is the maximum and rises linearly; the crossing is the minimum, at ζ = 0.362 for v₀ = 3 m/s (39.8 kN). Halving v₀ halves the spike exactly but reduces x_max by less than half (the static deflection x_s does not scale with v₀), so the spike line crosses the mid-stroke curve later: the optimum moves **up**. Code: v₀ = 1.5 m/s → ζ_opt = 0.472 (24.0 kN); v₀ = 6 m/s → ζ_opt = 0.313 (72.2 kN) ✓.

**Check.** At ζ = 0 the force is kX and its maximum is at maximum stroke (t* = 0.166 s), consistent with (b) since Ẍ = g − kX/m is then most negative where X is largest.

**Common mistakes.** Assuming the force peak and the stroke peak coincide; computing Ẍ by differentiating the sampled velocity with too coarse a grid (a 1 ms grid is fine; a 10 ms one smears the t = 0 spike).

## Check-yourself answers

1. ẍ + 2ζω_n ẋ + ω_n² x = 0 with ω_n = √(k/m), ζ = c/(2√(km)). c is in N s/m = kg/s; so is √(km).
2. r = ω_n(−0.6 ± 0.8i); ω_d = 0.8 ω_n.
3. The double root gives only one exponential; substituting x = u e^{−ω_n t} leaves ü = 0, so u = A + Bt and t e^{−ω_n t} is the second solution.
4. W(t) = W(0)e^{−ct/m}. An exponential is never zero, so W cannot vanish later if it is non-zero at t = 0.
5. δ = ln(12.0/4.3)/3 = 0.342; ζ = 0.342/√(4π² + 0.117) = 0.0544.
6. d/dt(½mẋ² + ½kx²) = −cẋ². With c = 0 the mechanical energy is constant: the motion is a pure exchange between kinetic and spring energy.
7. The gear arrives with velocity, and a viscous damper applies cv₀ at the first instant with the spring still slack; that spike grows with ζ while the mid-stroke spring force falls with ζ. The minimum is near ζ = 0.36 for the lesson's strut at 3 m/s.
8. Downward: noise raises the measured height of each peak by a few tenths of σ, which is a larger fraction of the small late peaks than of the large early ones, so the fitted decay is too slow. Raising the acceptance threshold to about 20σ (keeping only peaks for which the offset is under 2 %) fixed it.
