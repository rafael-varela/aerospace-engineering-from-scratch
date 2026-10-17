# P02.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. All simulated records are synthetic.

## 1. A one-sided pulse train ($s = 15$: d = 0.35)

**Derivation.** (a) a₀ = (2/T)∫₀^{dT} 1 dt = **2d**. a_n = (2/T)∫₀^{dT} cos(nω₀t) dt = (2/T)[sin(nω₀t)/(nω₀)]₀^{dT} = **sin(2πnd)/(πn)**. b_n = (2/T)[−cos(nω₀t)/(nω₀)]₀^{dT} = **(1 − cos(2πnd))/(πn)**. Amplitude: a_n² + b_n² = [sin²x + 1 − 2cos x + cos²x]/(πn)² = 2(1 − cos x)/(πn)² with x = 2πnd; since 1 − cos x = 2 sin²(x/2), the amplitude is **2|sin(πnd)|/(πn)**.

**Result.** (b) `fourier_coefficients` with 8192 samples: a₀ = 0.69995 (formula 0.7000); amplitudes of harmonics 1–6: **0.5672, 0.2575, 0.0331, 0.1514, 0.0901, 0.0327**, against the formula 0.5672, 0.2575, 0.0332, 0.1514, 0.0900, 0.0328 (the trapezoid rule at a jump is good to about 10⁻⁴ here; the individual a_n and b_n agree to the same level). (c) Mean square d = 0.35. Parseval with six harmonics: **0.3331** (95.2 % of the power); with 200 harmonics: **0.3495** (99.8 %): the 1/n decay of a jump makes the tail long. (d) Harmonic n vanishes when sin(πnd) = 0, i.e. nd an integer: for d = 0.35 = 7/20 that first happens at n = 20. The second harmonic vanishes when 2d is an integer, d = 0.5: the symmetric square wave, which is why the lesson's wave has only odd harmonics.

**Check.** d = 0.5 reduces the amplitude formula to 2|sin(πn/2)|/(πn) = 2/(πn) for odd n: the lesson's 4/(πn) for a wave of peak-to-peak 2, halved for a wave of peak-to-peak 1 ✓. The mean a₀/2 = 0.35 = d ✓.

**Common mistakes.** Integrating over the whole period with the pulse height treated as ±1; forgetting that a one-sided train has a mean (a₀ ≠ 0) and even harmonics.

## 2. Orthogonality, the square wave, and Parseval

**Derivation.** (a) sin(mω₀t) sin(nω₀t) = ½[cos((m − n)ω₀t) − cos((m + n)ω₀t)]. Over a full period ∫₀^T cos(kω₀t) dt = 0 for any non-zero integer k and = T for k = 0. Since m + n ≥ 2 the second term always integrates to zero; the first gives T/2 when m = n and 0 otherwise. For sin × cos: sin(mω₀t) cos(nω₀t) = ½[sin((m + n)ω₀t) + sin((m − n)ω₀t)], and every sine of a non-zero integer multiple integrates to zero over a period, while m = n gives ½ sin(2mω₀t), also zero.

(b) Multiply (1) by sin(nω₀t) and integrate over (0, T). The constant and every cosine term vanish by (a); every sine term with m ≠ n vanishes; the m = n term gives b_n T/2. So b_n = (2/T)∫ f sin(nω₀t) dt. Likewise for a_n with cos; and integrating (1) alone gives a₀ T/2 = ∫ f dt.

(c) b_n = (2A/T)[∫₀^{T/2} sin(nω₀t) dt − ∫_{T/2}^T sin(nω₀t) dt]. First integral: (1 − cos nπ)/(nω₀); second: (cos nπ − cos 2nπ)/(nω₀) = (cos nπ − 1)/(nω₀). Difference: 2(1 − cos nπ)/(nω₀) = 2(1 − (−1)ⁿ)/(nω₀). With 2/T × 1/ω₀ = 1/π: b_n = (2A/(nπ))(1 − (−1)ⁿ) = 4A/(πn) for odd n, 0 for even n.

(d) f² = [a₀/2 + Σ(a_n cos + b_n sin)]²; integrate over a period: cross terms vanish by (a), the squares give (a₀/2)²T + Σ(a_n² + b_n²)T/2. Divide by T: (1/T)∫f² = a₀²/4 + ½Σ(a_n² + b_n²). For the square wave the left side is A² and the right is ½Σ_odd 16A²/(π²n²) = (8A²/π²)Σ_odd 1/n², so Σ_odd 1/n² = π²/8. Then Σ_all 1/n² = Σ_odd + Σ_even = π²/8 + ¼Σ_all, so Σ_all = (π²/8)/(3/4) = π²/6.

**Check.** `parseval_power(*square_wave_coefficients(4001))` = 0.9998 → 1 ✓; (8/π²)(1 + 1/9 + 1/25) = 0.933, the three-term value, matches the code's 0.933 for N = 5.

**Common mistakes.** Forgetting the k = 0 case in (a) when m = n; in (c), getting the second integral's sign wrong and finding b_n = 0 for all n.

## 3. Ripple on a rectified 400 Hz bus

**Assumptions.** Ideal rectification, no load; V₀ = 1 in the numbers.

**Derivation.** (a) |sin(2πf₀t)| repeats every half cycle of the sine: period **T_r = 1/(2f₀) = 1.25 ms**, so ω₀ = 2π/T_r = 4πf₀. v(−t) = |sin(−2πf₀t)| = v(t): even, so b_n = 0. With θ = 2πf₀t (θ runs 0 → π over one period), a_n = (2/T_r)∫₀^{T_r} sin(2πf₀t) cos(2nθ) dt = (2/π)∫₀^π sin θ cos 2nθ dθ. Product to sum: sin θ cos 2nθ = ½[sin((1 + 2n)θ) + sin((1 − 2n)θ)]; ∫₀^π sin(kθ) dθ = (1 − cos kπ)/k = 2/k for odd k. So a_n = (1/π)[2/(1 + 2n) + 2/(1 − 2n)] = (2/π)(2/(1 − 4n²)) = **−4V₀/(π(4n² − 1))**, and a₀ = (2/π)∫₀^π sin θ dθ = **4V₀/π**. (b) `fourier_coefficients` (T = 1.25 ms, 4096 samples): a₀ = 1.27324 (4/π = 1.27324); a₁…a₆ = **−0.42441, −0.08488, −0.03638, −0.02021, −0.01286, −0.00890**, against the formula −0.42441, −0.08488, −0.03638, −0.02021, −0.01286, −0.00890; largest |b_n| = 1.2 × 10⁻¹⁶. (c) The DC level is a₀/2 = 2V₀/π = 0.6366 V₀. The lowest ripple harmonic is n = 1 at **800 Hz**, amplitude 4/(3π) = 0.4244 V₀, i.e. **66.7 % of DC**. Rms ripple from the coefficients (six terms): √(½Σa_n²) = 0.3077 V₀, **48.3 % of DC**; from Parseval with the exact mean square ½: √(½ − 4/π²) = 0.3078 V₀, **48.3 %** ✓. (d) Required attenuation at 800 Hz: 0.01/0.667 = 0.015, so √(1 + (ωτ)²) = 66.7, ωτ = 66.7, τ = 66.7/(2π × 800) = **13.3 ms**. The higher harmonics start smaller (a₂/a₁ = 0.2) *and* are attenuated more (the filter falls as 1/ω), so if the first harmonic meets the target the others are below 0.2 × ½ = 0.1 of it; designing for the lowest harmonic is conservative by construction.

**Check.** Limiting cases: a_n → −V₀/(πn²) for large n, the 1/n² decay of a kink (|sin| has corners, not jumps) ✓. Σ a_n² from (a) equals ½ − 4/π² = 0.0947 by Parseval; the six-term sum is 0.0946 ✓.

**Common mistakes.** Using the period of sin (2.5 ms) instead of |sin|, which makes every odd harmonic of that period vanish and hides the structure; quoting the ripple as 42 % of V₀ instead of as a fraction of the DC level.

## 4. Choosing a pulse period ($s = 15$: T = 1.15 s)

**Result.** (a) Harmonics 3, 5, 7 are at **2.609, 4.348, 6.087 Hz**; the fifth is nearest the 5 Hz mode (r₅ = 0.870) but no longer on it. With 51 harmonics the dominant harmonic is **1** (amplitude 1.313 static deflections); harmonic 5 gives 1.042, harmonic 3 0.583; the peak of the summed response is **2.147 static deflections** against 14.33 at T = 1.00 s: changing the period by 15 % cuts the response by a factor 6.7. (b) Sweep (peak in static deflections): T = 0.90: 2.453; 1.00: **14.332**; 1.10: 2.494; 1.15: 2.147; 1.20: **2.060** (the least in 0.90–1.30; the 1.2 s pattern puts the 5th harmonic at 4.17 Hz and the 7th at 5.83 Hz, equally far from the mode). (c) `solve_ivp` to 40 s, last five seconds: peak 2.1468 against 2.1471 from the series; maximum pointwise difference **3.0 × 10⁻³** static deflections (0.14 %). (d) A mode in [4.8, 5.2] Hz is hit by harmonic n when T = n/f, i.e. T ∈ [n/5.2, n/4.8]: 5th harmonic **[0.962, 1.042] s**; 3rd **[0.577, 0.625] s**; 7th **[1.346, 1.458] s**. T = 1.15 s is outside all three, with the nearest band edge (1.042 s) 0.11 s away; so is 1.20 s.

**Check.** The response at T = 1.15 s is dominated by the fundamental with |H(0.174)| = 1.03: essentially the quasi-static square wave (amplitude 1) plus a modest fifth-harmonic contribution: 1.3 + 1.0 ≈ 2.1 ✓ in order of magnitude.

**Common mistakes.** Reading the peak from `harmonic_amplitudes` alone (the harmonics add with different phases, so the peak of the sum is less than the sum of the amplitudes: 2.147 < 1.313 + 1.042 + 0.583 + …); forgetting that a ± 0.2 Hz band in f is a *wider* band in T for higher harmonics.

## 5. What the series does and does not tell you

(a) The 9 % overshoot belongs to the *truncated* series S_N, which overshoots at a point T/(2(N + 1)) from the edge that moves to the edge as N grows; the full series converges to the square wave at every point except the jump itself, where it gives the midpoint. The torque is the square wave; nothing physical overshoots. Adding 9 % to the rating would be adding an artefact of a plotting choice. (The real valve, with its 40 ms rise, undershoots the corner if anything.)

(b) Ranking by effect on the 12.7: the damping ratio (response ∝ 1/ζ; ζ of a deployed array is known to a factor of about two, so **± 100 %**); the mode frequency (a 2 % shift, 0.1 Hz, takes r₅ from 1 to 0.98 and |H| from 50 to 21, **−60 %**); the valve rise time (lowers b₅ by a few per cent, see (c)); the Fourier coefficients themselves (exact). Measure ζ first, then f_mode, from a ring-down of the boom.

(c) A ramp of rise time t_r smooths the corner; harmonics with half-period shorter than about t_r, n ≳ T/(2t_r) ≈ 1/(2 × 0.04) = **12**, fall as 1/n² instead of 1/n (the extension's exact knee is n ≈ T/(πt_r) ≈ 8). The fifth harmonic is below the knee and its coefficient changes by a few per cent (a 40 ms linear ramp multiplies b₅ by sin(5π × 0.04)/(5π × 0.04) = 0.935): the answer to the engineering question (the 5th harmonic resonates; avoid T = 1.0 s) is unchanged, the number 12.7 becomes about 11.9.

(d) The resonant response grows as 1 − e^{−ζω_n t} with time constant 1/(ζω_n) = 1/(0.01 × 31.4) = **3.2 s**. In 20 s it reaches 1 − e^{−6.3} = 99.8 % of the steady value: yes, the boom reaches the table's amplitude long before the pattern changes. Only a pattern changed every second or two would keep it below, and then the response is a transient problem (P02.08), not a steady-state one.

**Check.** The 2 % shift in (b): |H(0.98, 0.01)| = 1/√((1 − 0.9604)² + (0.0196)²) = 1/√(0.00157 + 0.00038) = 22.6 ✓ (quoted 21, grid value).

**Common mistakes.** Treating the Gibbs overshoot as physics; assuming the exact Fourier coefficients are the uncertain part of the calculation.

## 6. The Gibbs constant and the convergence rate

**Derivation.** (a) With b_k = 4/(πk) for odd k, S_N = (4/π)Σ_{k odd ≤ N} sin(kω₀t)/k, and S_N′ = (4ω₀/π)Σ_{k odd ≤ N} cos(kω₀t). Let x = ω₀t and q = e^{2ix}. Σ_{k=1,3,…,N} e^{ikx} = e^{ix}(1 + q + … + q^{(N−1)/2}) = e^{ix}(1 − q^{(N+1)/2})/(1 − q) = e^{ix}(1 − e^{i(N+1)x})/(1 − e^{2ix}). Factor: numerator = −2i e^{i(N+1)x/2} sin((N + 1)x/2), denominator = −2i e^{ix} sin x; the ratio is e^{i(N+1)x/2} sin((N + 1)x/2)/sin x, whose real part is cos((N + 1)x/2) sin((N + 1)x/2)/sin x = sin((N + 1)x)/(2 sin x). Hence **S_N′ = (2ω₀/π) sin((N + 1)ω₀t)/sin(ω₀t)**.

(b) S_N′ > 0 just after t = 0 and first vanishes when (N + 1)ω₀t = π: **t* = T/(2(N + 1))**. `gibbs_overshoot`'s grid puts the maximum at t = 0.05003, 0.016684, 0.0050003 for N = 9, 29, 99 (T = 1), against T/(2(N + 1)) = 0.05, 0.016667, 0.005 ✓ (the grid resolution is 1/(800N)).

(c) S_N(t*) = ∫₀^{t*} S_N′ dt = (2ω₀/π)∫₀^{t*} sin((N + 1)ω₀t)/sin(ω₀t) dt. Substitute u = (N + 1)ω₀t, du = (N + 1)ω₀ dt, and for large N use sin(ω₀t) ≈ ω₀t = u/(N + 1) over the whole range (ω₀t* = π/(N + 1) is small): S_N(t*) ≈ (2ω₀/π)∫₀^π [sin u/(u/(N + 1))] du/((N + 1)ω₀) = (2/π)∫₀^π (sin u/u) du = (2/π) Si(π). `scipy.special.sici(np.pi)[0]` = 1.851937, so **(2/π) Si(π) = 1.178980**: the overshoot is 0.17898 = 8.949 % of the jump 2. The code's N = 99 value, 1.17901, differs by 3 × 10⁻⁵, the O(1/N²) remainder of the small-angle approximation.

(d) Parseval's tail: rms² = ½Σ_{odd n>N} 16/(π²n²) ≈ (8/π²)∫_N^∞ dn/(2n²) = 4/(π²N) (the odd integers are half of all integers, hence the ½ inside the integral), so **rms ≈ 2/(π√N)**: 0.1389 at N = 21 and 0.0707 at N = 81, against `truncation_error_rms` **0.1357 and 0.0703** (2.3 % and 0.6 % off, the integral approximation of a sum). The Gibbs ears have fixed height 0.18 but width of order T/(N + 1), so their contribution to the mean-square error is O(0.18² × 2/N) ≈ 0.06/N, the same order as the tail 0.4/N but smaller, and in any case *already included* in it: Parseval counts all of f − S_N, ears and all. The overshoot is invisible in the mean square because it is narrow, not because it is small; that is the difference between convergence in the mean and uniform convergence.

**Check.** Dimensions: S_N′ is in units of f per second, (2ω₀/π) × (dimensionless) ✓. Limit N → ∞ at fixed t > 0: t* → 0 and the ear squeezes against the jump, consistent with pointwise convergence at every fixed t.

**Common mistakes.** Summing the geometric series over all k instead of odd k only (ratio e^{ix} instead of e^{2ix}); in (c), integrating to T/(2N) instead of T/(2(N + 1)) and getting Si(πN/(N + 1)), which is fine in the limit but confusing at N = 9.

## Check-yourself answers

1. a_n = (2/T)∫₀^T f cos(nω₀t) dt, b_n = (2/T)∫₀^T f sin(nω₀t) dt, a₀ = (2/T)∫₀^T f dt. The factor is 2/T because ∫₀^T cos²(nω₀t) dt = T/2, not T: the basis functions have "length squared" T/2.
2. b₁ = 4/π = 1.273, b₃ = 4/(3π) = 0.424. The even b_n contain the factor 1 − cos nπ = 0 for even n; equivalently, the square wave is antisymmetric about T/4 and even harmonics are not.
3. As 1/n² (one integration by parts more than a jump); rms error ~ N^{−3/2}.
4. (1/T)∫f² = a₀²/4 + ½Σ(a_n² + b_n²). For the unit square wave: 1 = ½Σ_odd 16/(π²n²), so 1 + 1/9 + 1/25 + … = π²/8 = 1.2337.
5. 8.95 % of the jump height ((2/π)Si(π) − 1)/2); it does not decrease as N → ∞, it only moves closer to the jump and narrows.
6. c_n = (a_n − i b_n)/2 for n ≥ 1, c₀ = a₀/2. Because a_n and b_n are real, c_{−n} = (a_n + i b_n)/2 = c̄_n; the pair c_n e^{inω₀t} + c_{−n} e^{−inω₀t} is then real.
7. The 7th harmonic (7 × 1 Hz). Forcing amplitude 4/(7π) = 0.1819; response 0.1819 × 1/(2 × 0.02) = 4.55 static deflections; the response lags that harmonic by 90°.
8. Avoid T = 0.6, 1.0 and 1.4 s (harmonics 3, 5, 7 on 5 Hz); 1.2 s is safe (harmonics 5 and 7 at 4.17 and 5.83 Hz), as is 0.8 s (3.75 and 6.25 Hz).
