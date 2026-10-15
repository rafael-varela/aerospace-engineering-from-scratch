# P02.19  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are from `solution.py` and given to the figures shown. Every record is synthetic.

## 1. Dryden numbers by hand

**Givens** ($s = 15$). σ_w = 1.75 m/s, L = 350 m, V = 50 m/s.

**Derivation.** (a) b = V/L = 50/350 = **0.1429 rad/s** = **0.0227 Hz**. Φ_w(0) = σ_w²L/(πV) = 3.0625 × 350/(π × 50) = **6.824 (m/s)² per rad/s**.

(b) With g(x) = (1 + 3x²)/(1 + x²)², the quotient rule gives g′(x) = [6x(1 + x²) − 4x(1 + 3x²)]/(1 + x²)³ = 2x(1 − 3x²)/(1 + x²)³. It vanishes at x = 0 (a minimum: g′ > 0 just above 0) and at **x = 1/√3**. There g = (1 + 1)/(4/3)² = 2 × 9/16 = **9/8**. The peak is at ω = b/√3 = **0.0825 rad/s** (0.0131 Hz) with Φ_w = 9/8 × 6.824 = **7.677**.

(c) Integral time scale L/(2V) = **3.5 s**; R_w = 0 first at τ = 2L/V = **14.0 s**.

(d) `dryden_psd(0.0825, 1.75, 350, 50)` = 7.6767, and at ±10⁻⁴ rad/s either side 7.67674 (lower in the sixth figure, so it is a maximum); `dryden_autocorrelation(14.0, ...)` = 0 to rounding.

**Check.** Units: σ²L/V = (m/s)² m/(m/s) = (m/s)² s = (m/s)² per rad/s ✓. Limit: as L → ∞ the spectrum piles up at ω → 0 and Φ_w(0) → ∞ while σ² stays fixed, consistent with a process that is nearly constant in time.

**Common mistakes.** Reporting b in Hz as if it were rad/s; taking the peak of Φ_w at ω = 0 (true for u, not for w); forgetting that the first zero of R_w is at 2L/V, not L/V.

## 2. Normalisation, transform pair and scatter

**Derivation.** (a) ω = bx gives ∫₀^∞ Φ_w dω = (σ²/π)∫₀^∞ (1 + 3x²)/(1 + x²)² dx. With x = tan θ, dx = sec²θ dθ and 1 + x² = sec²θ: ∫₀^∞ dx/(1 + x²)² = ∫₀^{π/2} cos²θ dθ = π/4 and ∫₀^∞ x²dx/(1 + x²)² = ∫₀^{π/2} sin²θ dθ = π/4. Sum: (σ²/π)(π/4 + 3π/4) = σ² ∎.

(b) F{e^{−a|τ|}} = ∫_{−∞}^{∞} e^{−a|τ|}e^{−iωτ}dτ = 2a/(a² + ω²). Differentiating under the integral with respect to a: F{−|τ|e^{−a|τ|}} = 2/(a² + ω²) − 4a²/(a² + ω²)² = 2(ω² − a²)/(a² + ω²)², so F{|τ|e^{−a|τ|}} = 2(a² − ω²)/(a² + ω²)². Then S_w = σ²[2a/(a² + ω²) − (a/2)·2(a² − ω²)/(a² + ω²)²] = σ²a[2(a² + ω²) − (a² − ω²)]/(a² + ω²)² = σ²a(a² + 3ω²)/(a² + ω²)². Dividing by π and writing a = V/L, x = ω/a: σ²(L/(πV))(1 + 3x²)/(1 + x²)² ∎.

(c) ∫₀^∞(1 − aτ/2)²e^{−2aτ}dτ = ∫₀^∞(1 − aτ + a²τ²/4)e^{−2aτ}dτ = 1/(2a) − a·1/(2a)² + (a²/4)·2/(2a)³ = 1/(2a) − 1/(4a) + 1/(16a) = 5/(16a), so ∫₀^∞R_w² dτ = 5σ⁴/(16a) = 5σ⁴L/(16V). For a zero-mean Gaussian process, Cov(x²(t), x²(t′)) = 2R²(t − t′), so Var[(1/T)∫₀^T x²dt] = (1/T²)∫∫2R²(t − t′)dt dt′ ≈ (2/T)∫_{−∞}^{∞}R²dτ when T ≫ L/V (the double integral over the square becomes T times the integral over the lag). That is (2/T)·2·5σ⁴L/(16V) = 5σ⁴L/(4VT), relative SD √(5L/(4VT)) ∎. Code check: `quad` of ρ² for the $s = 15$ parameters gives 2.1875 = 5/(16 × 0.142857) ✓.

(d) ω²Φ_w = σ²(L/(πV))ω²(1 + 3x²)/(1 + x²)² → σ²(L/(πV))·3x²ω²/x⁴ = 3σ²V/(πL) as x → ∞: a constant, so ∫ω²Φ_w dω diverges linearly. The derivative of a Dryden gust has infinite variance: the model's gust velocity is continuous but nowhere differentiable, like Brownian motion, and any quantity that needs ẇ (a crossing rate, a gust gradient) is undefined in it.

**Common mistakes.** Dropping the factor 2 that comes from the two-sided integral in (c); using 1/(N − k) estimates to argue the bias is zero and the variance follows; in (b) forgetting that the transform of an even function is real so no imaginary part can survive.

## 3. Gyro noise integrates to an angle random walk ($s = 15$: G₀ = 0.035² = 1.225 × 10⁻³ (deg/s)²/Hz)

**Assumptions.** White, Gaussian, zero-mean rate noise; a perfect integrator; no bias.

**Derivation.** (a) Φ_θ(ω) = |H|²Φ_n = Φ_n/ω², and ∫₀^∞ dω/ω² diverges at ω → 0: a stationary θ does not exist. But θ(t) = ∫₀^t n dt′ at fixed t is a well-defined random variable; splitting [0, t] into m independent pieces of length t/m, each with variance S₀t/m (the area of the δ-correlation times the piece length), the variances add: Var θ(t) = S₀t. Converting: in this lesson's convention the one-sided density per hertz is G = 2πΦ = 2S, so S₀ = G₀/2 = 6.125 × 10⁻⁴ (deg/s)²·s.

(b) σ_θ(600 s) = √(6.125 × 10⁻⁴ × 600) = **0.606°**. It grows as √t: the **angle random walk**.

(c) With dt = 0.01 s, `white_noise` has variance 100; scaled by √S₀ and summed times dt, θ(600) over 200 seeds (1000–1199) has SD **0.590°**, 2.7 % below 0.606°. An SD estimated from n = 200 samples has relative scatter ≈ 1/√(2n) = 5 %, so this is a 0.5σ difference.

(d) `welch_psd` of one record with 4096-point segments, averaged over 1–40 Hz: **1.232 × 10⁻³** against G₀ = 1.225 × 10⁻³ (+0.6 %; with K = 28 segments and 3200 bins the band mean scatters by ~0.4 %).

**Check.** Dimensions: (deg/s)²/Hz × s = (deg/s)² s² = deg² ✓. Limiting case: a 10× quieter gyro (G₀/100) gives a 10× smaller σ_θ for the same t, and the same σ_θ needs 100× the time: the usual √t trade.

**Common mistakes.** Using G₀ instead of G₀/2 (a factor √2 in σ_θ); summing without the dt (a factor 1/dt in variance); expecting `welch_psd` of the *integrated* record to be flat.

## 4. A heavier aircraft ($s = 15$: m = 6500 kg)

**Result.** (a) τ = 2 × 6500/(1.19 × 100 × 30 × 5.5) = **0.662 s**; b = 100/300 = 0.333 rad/s, c = 1.510 rad/s. From (11): σ_v = **1.295 m/s**, σ_a = **1.144 m/s²**, σ_n = **0.1166 g**. The sharp-edged value σ_w/(gτ) = 0.231 g is twice the spectral RMS; the heavier aircraft at r = c/b = 4.5 follows less of the gust (σ_v/σ_w = 0.86 against 0.95 for the lesson aircraft) and is less bumpy per unit speed, but it flies faster, so σ_n is only 14 % below the lesson aircraft's 0.135 g.

(b) Seed 15, 2000 s: sample σ_w = 1.491 m/s (model 1.5; variance scatter √(5 × 300/(4 × 100 × 2000)) = 4.3 %, RMS 2.2 %: a 0.3σ draw), σ_v = 1.290 m/s (model 1.295), σ_n = **0.1153 g** (model 0.1166, −1.1 %; the acceleration decorrelates in ~τ, so its RMS scatters by ~1 %).

(c) P(|Δn| > 0.25 g) = erfc(0.25 × 9.807/(1.144√2)) = **3.21 %** from the Gaussian model; the record gives **3.13 %**.

(d) Welch/model band ratios: **1.07, 0.96, 1.00** in 0.01–0.1, 0.1–1, 1–10 Hz, each within the segment scatter (7 %, 5 %, 1 %). The load-factor spectrum peaks at **0.128 Hz**. Between b and c, |H_a|² rises as ω² while Φ_w falls as ω⁻², so their product would be flat; it is bent down by the roll-off of |H_a|² above c and by the flattening of Φ_w below b, and the maximum lands near the geometric mean √(bc)/(2π) = 0.113 Hz (the exact maximum of the rational product is a little higher because the Φ_w tail has the factor 3).

**Check.** Scaling: with r ≫ 1, σ_a² ≈ 1.5σ_w²bc; for the lesson aircraft 1.5 × 2.25 × 0.2 × 2.857 = 1.929 against the exact 1.763 (r = 14); for this aircraft r = 4.5 and the approximation is 30 % high, which is why the exact (11) is used. Trace check: the sample variance of Δn should equal the integral of its Welch spectrum: Σ P Δf = 0.01328 g² against 0.1153² = 0.01330 g² ✓ (the 0.1 % missing is the energy below the first Welch bin, removed with each segment's mean).

**Common mistakes.** Using m in the time constant without the factor 2; quoting the exceedance fraction without noticing it is a 2.2σ tail that moves 16 % for every 3 % of σ_a; comparing the Welch estimate with Φ_w in rad/s units instead of G in Hz units (a factor 2π).

## 5. What the model can and cannot say

(a) σ_w is set by the integral of Φ_w, which converges as ∫dω/ω² beyond b: the tail above 10b contributes 3σ²V/(πL × 10b) ÷ σ² = 3/(10π) = 9.5 % of the variance, so σ_w is tail-insensitive to within a few per cent. σ_v is even less sensitive (the aircraft filters the tail). σ_a is tail-sensitive: over half of σ_a² (55 % by `quad`) comes from above c, where Dryden and von Kármán differ by a growing factor, so a −5/3 tail would lower σ_a by perhaps 10 %. The crossing rate of v is a ratio σ_a/σ_v and inherits σ_a's sensitivity. The exceedance fraction at 0.5 g is the most sensitive of all: a 3.7σ tail, where 10 % in σ_a is a factor of 3 in probability.

(b) The variance of successive 5-minute blocks (or the Welch level of each block in the 0.1–1 Hz decade): (13) says blocks of 300 s scatter by 14 % in variance, so blocks that differ by more than ~40 % (3σ) are from different processes. The trustworthy record for σ_w is the longest run of blocks that agree; over mixed terrain that may be 10–15 minutes, giving σ_w to 5–7 % at best, not the 2 % of 260 minutes.

(c) Rice's formula needs σ_ȧ² = ∫ω²|H_a|²Φ_w dω, whose integrand tends to 3σ_w²V/(πLτ²), a constant: the integral diverges and the formula gives an infinite crossing rate for the load factor. The missing physics is the roll-off of the lift response above a few hertz: unsteady aerodynamics (the lift takes time to build as the wing crosses a gust, a lag of order c/V with c the chord) and the structure's filtering. With either added, |H_a|² falls as ω⁻² or faster beyond its corner and σ_ȧ is finite.

(d) With r ≫ 1, σ_a² ≈ (3/2)σ_w²bc and both b and c are proportional to V, so σ_n ∝ V and the criterion is met at V = 60 × 0.10/0.1354 = **44 m/s**, a 26 % speed reduction. The exact (11) at 44.3 m/s gives 0.0999 g. Before trusting it to 10 %: the quasi-steady, plunge-only model. At the lower speed the aircraft flies at a higher C_L, closer to stall where C_Lα is not constant, and the short-period pitch mode, which the model omits, changes the gust response by tens of per cent in either direction. The number 44 m/s is a first estimate, not a limit.

**Common mistakes.** Treating the ω⁻² tail as physics; using the crossing rate of v as if it were the bump rate; applying the V-scaling without stating that it holds only for r ≫ 1 (r = L/(Vτ) = Lρ S C_Lα/(2m) = 14.3, the same at every speed because τ ∝ 1/V).

## 6. Designing a flight-test record ($s = 15$)

**Result.** (a) T = 5L/(4V × 0.10²) = 1500/(240 × 0.01) = **625 s** (10.4 min).

(b) Resolution 0.005 Hz needs 200 s segments = **10 000 samples** at 50 Hz. In 625 s (31 250 samples) with 50 % overlap, K = floor((31 250 − 5000)/5000) = **5** segments; per-bin relative SD ≈ 1/√5 = **45 %**.

(c) Seed 15: variance ratio to σ_w² **1.007** (predicted scatter 10 %); the estimate at 0.01 Hz is **0.48** of the model (a single bin at 45 % scatter: a 1.2σ draw); the band average over 0.01–0.1 Hz (18 bins) is **0.95**; the per-bin SD/mean between 0.05 and 0.5 Hz is **0.37**, consistent with the 0.45 predicted for five Hann segments with overlap (overlap makes the effective K a little larger than 5).

(d) The variance requirement fixes T; the resolution requirement fixes the segment length; K is their ratio, and 625 s divided by 200 s cannot give more than five segments however the record is cut. The customer can have a variance to 10 % and *either* a spectrum at 0.005 Hz resolution with 45 % per-bin error *or* a spectrum smoothed to 0.05 Hz (ten bins averaged) with about 15 % error, but not a fine spectrum with small error from a 10-minute record. The honest report: σ_w² = (sample value) ± 10 %; the Welch spectrum at 0.005 Hz with ±45 % error bars on each bin, drawn on log axes so the bars are visible; and a band-averaged spectrum (three bins per decade, ~15 %) for comparison with the Dryden form. To halve the per-bin error, fly four times longer and check stationarity block by block.

**Check.** Σ P Δf over the Welch estimate is 0.90 of the record variance. The missing 10 % is real: each 200 s segment has its mean removed, which discards the energy below the first bin at 0.005 Hz, and for Dryden turbulence with b = 0.032 Hz that band holds about a tenth of σ_w². Short segments cost low-frequency energy as well as resolution; the 8192-point segments of the lesson (first bin 0.006 Hz) pay the same price, which is why the lesson quotes the variance from the record, not from the spectrum. Limit: with no overlap K = 3 and the per-bin error would be 58 %; the overlap is free information and is why 50 % is the default.

**Common mistakes.** Counting K as T/T_seg = 3.1 and forgetting the overlap, or as 6.25 and forgetting that partial segments are dropped; quoting the per-bin scatter as the scatter of a band average; believing a longer FFT gives a better spectrum (it gives a finer one).

## Check-yourself answers

1. A process is WSS when its mean is constant and E[x(t)x(t + τ)] depends only on τ. From the definition alone: R is even, and |R(τ)| ≤ R(0) (Cauchy–Schwarz).
2. Φ(ω) = (1/π)∫_{−∞}^{∞}R(τ)e^{−iωτ}dτ = (2/π)∫₀^∞R(τ)cos ωτ dτ for ω ≥ 0, and R(τ) = ∫₀^∞Φ(ω)cos ωτ dω. The check: ∫₀^∞Φ dω must equal R(0) = σ².
3. Φ_w = σ_w²(L/(πV))(1 + 3x²)/(1 + x²)², x = Lω/V. Φ_w(0) = σ_w²L/(πV); as ω → ∞, Φ_w → 3σ_w²V/(πLω²).
4. The discrete spike σ_e² at lag 0 stands for a continuous δ-function of area σ_e²Δt, and S₀ = 1 needs that area to be 1. With variance 1 the generated gust has variance Δt σ_w²: for Δt = 0.02 s, 2 % of the intended value, i.e. an RMS 7 times too small.
5. Each periodogram bin is the squared modulus of a complex Gaussian: an exponential variable with SD equal to its mean, whatever N; a longer record adds bins rather than averaging them. Welch trades frequency resolution (segments of length N/K have bins K times wider) for a factor 1/√K in scatter.
6. Most of the gust energy is below the break frequency 0.032 Hz, where the aircraft follows the gust (v ≈ w_g) and feels no acceleration; only the ω⁻² tail above 1/τ = 0.45 Hz produces acceleration, and it holds a small fraction of σ_w². The sharp-edged gust puts all of its energy into an instant step.
7. The Dryden spectrum falls only as ω⁻², so ∫ω²Φ_w dω diverges: w_g has no finite derivative variance, nor does a = dv/dt (its spectrum tends to a constant). v has a finite derivative variance because |H_v|² adds a further ω⁻² roll-off.
8. T = 5L/(4V·0.05²) = 1500/(240 × 0.0025) = 2500 s (42 min). T ∝ 1/V: flying faster sweeps through more independent gusts per second.
