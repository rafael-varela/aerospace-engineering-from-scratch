# P02.07  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. All simulated records are synthetic.

## 1. Polar form, a power and three roots ($s = 15$: z = 4.5 + 4i)

**Derivation and result.** (a) |z| = √(4.5² + 4²) = √36.25 = **6.0208**; arg z = atan2(4, 4.5) = **0.72664 rad = 41.63°**. (b) z⁵: modulus 6.0208⁵ = **7911.7**, argument 5 × 0.72664 = 3.6332 rad, minus 2π gives **−2.6500 rad (−151.8°)**; z⁵ = 7911.7 (cos(−2.6500) + i sin(−2.6500)) = **−6974.7 − 3734.8i**. `z ** 5` gives (−6974.71875 − 3734.75i): agreement to the rounding of the hand values. (c) Common modulus 6.0208^{1/3} = **1.8192**; arguments (41.63° + 360°k)/3 = **13.88°, 133.88°, 253.88° (= −106.12°)**; roots 1.7661 + 0.4363i, −1.2609 + 1.3113i, −0.5052 − 1.7477i. Product: (4.5 + 4i) to 10⁻¹⁵ ✓; sum: 3 × 10⁻¹⁶ ✓. (d) The three roots are r^{1/3}e^{iθ/3} times 1, e^{2πi/3}, e^{4πi/3}, and the three cube roots of unity sum to zero (geometric series with ratio e^{2πi/3}, numerator 1 − e^{2πi} = 0).

**Check.** The product of the roots of w³ − z = 0 is +z (Vieta: product = −(−z) for a cubic with leading coefficient 1), consistent with (c). Dimensions: none; the modulus 7911.7 = 36.25^{2.5} ✓.

**Common mistakes.** Reporting 3.633 rad for arg z⁵ without reducing it; taking only the principal cube root; using atan(y/x) instead of atan2, which fails in the left half-plane (not here, but the habit matters).

## 2. Euler, the addition formula, and the complex gain

**Derivation.** (a) e^{iθ} = Σ (iθ)^k/k!. Powers of i cycle 1, i, −1, −i. Even k = 2j: (iθ)^{2j} = (−1)^j θ^{2j}; odd k = 2j + 1: (iθ)^{2j+1} = i(−1)^j θ^{2j+1}. So e^{iθ} = Σ_j (−1)^j θ^{2j}/(2j)! + i Σ_j (−1)^j θ^{2j+1}/(2j+1)! = cos θ + i sin θ, the two bracketed series being the Maclaurin series of cos and sin (P01.05), both absolutely convergent so the rearrangement is legitimate.

(b) e^{i(a+b)} = cos(a + b) + i sin(a + b). Also e^{ia}e^{ib} = (cos a + i sin a)(cos b + i sin b) = (cos a cos b − sin a sin b) + i(sin a cos b + cos a sin b). Equating real and imaginary parts: cos(a + b) = cos a cos b − sin a sin b, sin(a + b) = sin a cos b + cos a sin b.

(c) F₀ cos ωt = Re(F₀e^{iωt}). Let x = Re(X e^{iωt}). Because m is a real constant, m ẍ = m Re(d²/dt² X e^{iωt}) = Re(m(iω)² X e^{iωt}) (differentiation commutes with Re for a real variable t, and a real constant moves inside Re). Likewise c ẋ = Re(c iω X e^{iωt}) and k x = Re(k X e^{iωt}). Linearity of Re lets us add: Re[(−mω² + icω + k)X e^{iωt}] = Re[F₀ e^{iωt}] for all t. Write the left coefficient as P and evaluate at ωt = 0 and ωt = −π/2: Re P = Re F₀ = F₀ and Re(−iP) = Im P = Im F₀ = 0, so P = F₀ and X = F₀/(k − mω² + icω). Linearity was used to add the three terms and to move constants inside Re; realness of m, c, k was used in "m Re(w) = Re(m w)", which is false for complex m.

(d) D(u) = (1 − u)² + 4ζ²u with u = r². dD/du = −2(1 − u) + 4ζ² = 0 ⟹ u = 1 − 2ζ², a minimum since d²D/du² = 2 > 0, valid when u > 0, i.e. ζ < 1/√2. D_min = (2ζ²)² + 4ζ²(1 − 2ζ²) = 4ζ²(1 − ζ²), so |H|_peak = 1/(2ζ√(1 − ζ²)). For ζ = 0.08: r_peak = √(1 − 0.0128) = **0.99358**, |H|_peak = 1/(0.16 × 0.99679) = **6.2701** (against 6.25 at r = 1).

**Check.** (a) at θ = π: the cosine series gives −1 and the sine series 0; `exp_series(1j*np.pi)` returns −1 + 3 × 10⁻¹⁶i. (d) ζ → 0 gives r_peak → 1 and |H| → ∞; ζ = 1/√2 gives r_peak = 0, the flattest response (no peak).

**Common mistakes.** In (c), substituting x = X e^{iωt} without Re and then wondering why x is complex: the answer is still right because the ODE is linear, but the justification is the one above. In (d), maximising |H| instead of minimising D and losing a factor in the chain rule.

## 3. Two unsynchronised propellers ($s = 15$: f₂ = 40.25 Hz, Δf = 0.25 Hz)

**Assumptions.** Two pure sinusoids, no structural filtering, synthetic record.

**Derivation.** (a) 0.9 cos(2πf₂t) = Re(0.9 e^{i2π(f₂−f₁)t} e^{i2πf₁t}), so x = Re(Z(t) e^{iω₁t}) with **Z(t) = 1.2 + 0.9 e^{i2πΔf t}**: a 40 Hz phasor whose second arm turns once every 1/Δf seconds. (b) |Z(t)| = √(1.2² + 0.9² + 2 × 1.2 × 0.9 cos(2πΔf t)) by the cosine rule; it runs from 2.1 to 0.3 and back with period 1/Δf = **4.0 s**. (c) Fits at 40 Hz over 0.1 s windows (synthetic record, 100 kHz): start 0.0 s: **2.1003** (envelope at 0.05 s: 2.0984); 0.5 s: **1.9114** (1.9112); 1.0 s: **1.4400** (1.4424); 1.5 s: **0.7767** (0.7794). The fit sees the mean phasor over the window, so it agrees with the envelope at the window centre to about 0.3 %; the fitted *phase* also drifts (1.9°, 20.9°, 38.4°, 48.5°), which is the arm turning. (d) The cabin vibration swells and fades every 4 s between 0.3 and 2.1 m/s²: that is the throb, and it is more noticeable than a steady 2.1 would be because the ear and the body sense change. The synchrophaser makes Δf = 0 and holds the arm at a fixed angle, so the amplitude is constant, and if the angle is 180° it is constant at the minimum.

**Check.** The envelope formula is the lesson's (3) with δ = 2πΔf t; at t = 0 it gives 2.1, the maximum, consistent with both cosines starting at their peaks. Limiting case Δf → 0: the envelope is constant and the fit over any window gives the same phasor.

**Common mistakes.** Fitting at the *mean* frequency 40.125 Hz, which gives a phasor that is nearly constant but misses the physics of the problem; using a window longer than a fraction of the beat period, which averages the envelope.

## 4. Redesign the tray mount ($s = 15$: δ = 165°)

**Assumptions.** Linear mount, one frequency, steady state; m = 2.5 kg, ζ = 0.08, g = 9.80665 m/s².

**Result.** (a) `add_sinusoids([(1.2, 0), (0.9, radians(165))])` gives **0.4045 m/s²** at **35.16°**; floor displacement 0.4045/(2π × 40)² = **6.40 µm**. (b) T(4/3, 0.08) = **1.2678**; tray acceleration **0.513 m/s²**. (c) f_n = **15.84 Hz** (r = 2.525); ω_n = 99.5 rad/s; k = 2.5 × 99.5² = **24 761 N/m**; c = 2 × 0.08 × 2.5 × 99.5 = **39.8 N s/m**; sag = g/ω_n² = **0.99 mm**. (d) At 80 Hz: r = 5.05, T = **0.0524** (old mount 0.177); at 20 Hz: r = 1.263, T = **1.625** (old mount at 20 Hz: 1.778): the new mount still amplifies an idle-speed excitation, by 63 %. ODE check with base amplitude Y = 6.40 µm: fitted tray amplitude after 2 s = 1.2807 µm, ratio **0.2000** against T = 0.2000 (the transient at 2 s is e^{−ζω_n t} = 1.2 × 10⁻⁷).

**Check.** The new mount's T at 40 Hz is 0.2 by construction; the independent integration reproduces it to 4 × 10⁻¹⁰. Dimensions: k/m = s⁻² ✓; c/(2√(km)) = 39.8/(2√(24761 × 2.5)) = 0.080 ✓.

**Common mistakes.** Using acceleration where displacement is needed in the ODE (the ratio is the same, so the error is invisible until the units are checked); forgetting the c ẏ term on the right-hand side, which changes T from (6) to |H| and gives 0.159 instead of 0.200.

## 5. Can one synchrophaser angle quiet the whole cabin?

**Result.** (a) Tray, δ = 180°: √(1.44 + 0.81 − 2.16) = **0.300 m/s²**. Two rows aft, δ + 60° = 240°: √(1 + 1 + 2 cos 240°) = **1.000 m/s²**. (b) The tray alone is minimised at δ = **180°** (0.300); the aft location alone at δ + 60° = 180°, i.e. δ = **120°** (0.000). The sum of the two amplitudes, swept in 1° steps, is least at δ ≈ **135°**, where the tray sees 0.850 and the aft seat 0.518, sum 1.37, against 1.30 at δ = 180° and 1.50... recompute: at δ = 180° the sum is 0.300 + 1.000 = 1.300; at δ = 120°, 1.044 + 0 = 1.044; at δ = 135°, 0.850 + 0.518 = 1.368. The minimum of the *sum* is therefore at δ = 120°, 1.044, where the aft seat is silent and the tray is at 1.04 m/s². (With other weightings, say the square of the amplitudes, the compromise moves; the point is that it is a weighting choice.) (c) A real spectrum contains the 80 Hz blade-passing line, which has its own phasors and its own best angle (the blade-passing phase difference is 2δ, so its minimum is at δ = 90° or 270°, not 180°), and broadband aerodynamic noise, which has no phase and does not move the 40 Hz minimum but masks it in a measurement. The first changes the recommendation; the second changes only how carefully it must be measured. (d) A fixed δ is a compromise between locations that cannot all be at their minima at once; a schedule that depends on flight phase (who sits where, which line dominates) is what production synchrophasers do, and the numbers here show why: 180° is best for the tray (0.30) and poor aft (1.00), 120° is the reverse.

**Check.** At δ = 120° the aft phasors are 1 at 0° and 1 at 180°: exact cancellation, as the arithmetic shows.

**Common mistakes.** Assuming the location with equal amplitudes "obviously" wins because it can reach zero; the tray's 0.3 residual at its own optimum is set by the amplitude mismatch and no angle removes it.

## 6. Why a balanced propeller has no 1/rev force

**Derivation.** (a) S = Σ_{k=0}^{N−1} q^k with q = e^{2πi/N}. For N ≥ 2, q ≠ 1 and the geometric sum is S = (1 − q^N)/(1 − q) = (1 − e^{2πi})/(1 − q) = 0. The N blade phasors are the N-th roots of unity and close into a polygon.

(b) With δ = 2πh/N, |S_h| = |sin(πh)/sin(πh/N)|. The numerator sin(πh) = 0 for every integer h. The denominator sin(πh/N) is non-zero unless h/N is an integer, so |S_h| = 0 whenever N does not divide h. When N divides h every phasor has phase 2πhk/N ≡ 0 and the sum is N by direct count (the formula is 0/0 there). `add_sinusoids` confirms (|sum| to 10⁻¹⁵ where zero):

| h | N = 2 | N = 3 | N = 4 |
|---|---|---|---|
| 1 | 0 | 0 | 0 |
| 2 | 2 | 0 | 0 |
| 3 | 0 | 3 | 0 |
| 4 | 2 | 0 | 4 |
| 5 | 0 | 0 | 0 |
| 6 | 2 | 3 | 0 |

A two-blade propeller transmits every even harmonic, a three-blade every third, a four-blade every fourth: more blades push the first transmitted harmonic higher and quieter.

(c) `add_sinusoids([(1.05, 0), (1.0, 2π/3), (1.0, 4π/3)])` = **0.0500 F at 0°**. The balanced set sums to zero, so what remains is exactly the excess (1.05 − 1)F of the heavy blade at its own phase: the residual 1/rev force is the imbalance of one blade, at that blade's phase, which is what a dynamic balancer measures and corrects.

(d) Σ_{k=0}^{N−1} e^{ikδ} = (1 − e^{iNδ})/(1 − e^{iδ}). Factor: 1 − e^{iNδ} = e^{iNδ/2}(e^{−iNδ/2} − e^{iNδ/2}) = −2i e^{iNδ/2} sin(Nδ/2), and 1 − e^{iδ} = −2i e^{iδ/2} sin(δ/2). The ratio is e^{i(N−1)δ/2} sin(Nδ/2)/sin(δ/2), whose modulus is |sin(Nδ/2)/sin(δ/2)| since |e^{iα}| = 1.

**Check.** The formula at N = 3, δ = 50° gives 2.2856, and `add_sinusoids([(1, 0), (1, 50°), (1, 100°)])` gives 2.2856. Limiting case δ → 0: sin(Nδ/2)/sin(δ/2) → N, all phasors aligned.

**Common mistakes.** Treating the 0/0 case by the formula instead of by direct count; forgetting that the blade-passing harmonic h = N carries the *full* N F, which is why the blade-passing line, not the 1/rev, dominates a balanced propeller's signature.

## Check-yourself answers

1. 1 + i = √2 e^{iπ/4}; (1 + i)⁸ = (√2)⁸ e^{i2π} = 16.
2. e^{iθ} = cos θ + i sin θ; the Maclaurin series of e^z, cos θ and sin θ.
3. 1, e^{2πi/3} = −½ + (√3/2)i, e^{4πi/3} = −½ − (√3/2)i; sum 0, because it is a geometric series with ratio e^{2πi/3} whose numerator 1 − e^{2πi} vanishes (or: the three arrows form a closed triangle).
4. Z = 3 e^{−iπ/6}; derivative iωZ = 30 e^{i(π/2 − π/6)} = 30 e^{iπ/3} (amplitude 30, leads by 90°); integral Z/(iω) = 0.3 e^{−i(π/6 + π/2)} = 0.3 e^{−i2π/3}.
5. √(1.2² + 0.9²) = 1.500. For the sum to equal 1.2: 1.44 = 1.44 + 0.81 + 2.16 cos δ gives cos δ = −0.375, δ = 112.0°.
6. Because for real constants a, a Re(w) = Re(aw) and Re commutes with sums and with d/dt, so the real part of a complex solution is a real solution. ½mẋ² is a product: Re(w)² ≠ Re(w²), so the energy must be formed from the real signal (or with the conjugate: |w|²/2 for the mean).
7. |H| = 1/(2 × 0.05) = 10, lag 90°. T(2, 0.1) = √(1 + 0.04)/√(9 + 0.04) = 0.356 < 1.
8. 950 rpm is 15.8 Hz, which is the new mount's own natural frequency: at idle the mount would resonate with magnification about 1/(2ζ) ≈ 6 and transmissibility about 6.3.
