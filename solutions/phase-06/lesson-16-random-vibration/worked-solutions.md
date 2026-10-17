# P06.16  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and a probe script that calls it, given to the figures shown. Every time history is SYNTHETIC. The workmanship spectrum is NASA-STD-7001C Table 2; the Problem 1 spectrum is illustrative.

## 1. A qualification-style spectrum (s = 15: G_p = 0.175 g²/Hz)

**Assumptions.** Straight lines on log-log axes between breakpoints; one-sided PSD in g²/Hz.

(a) +6 dB/oct is n = 6.02/3.01 = 2, so G(20) = G_p(20/50)² = **0.16 G_p = 0.028 g²/Hz**; −6 dB/oct is n = −2, so G(2000) = G_p(800/2000)² = **0.16 G_p = 0.028 g²/Hz**. `segment_slopes_db_per_octave` returns +6.02, 0, −6.02 dB/oct.

(b) Rising segment: G₁f₁[(f₂/f₁)³ − 1]/3 = G_p(50/3)(1 − 0.4³) = 15.6 G_p = **2.730 g²**. Plateau: 750 G_p = **131.25 g²**. Falling: G₁f₁[(f₂/f₁)⁻¹ − 1]/(−1) = 800 G_p (1 − 0.4) = 480 G_p = **84.00 g²**. Total 1245.6 G_p = **217.98 g²**, **Grms = 14.764**. `psd_area` and `grms` agree to round-off. (With G_p = 0.16 the same spectrum gives 14.12 Grms.)

(c) The trapezoid rule on the four breakpoints gives **16.003 Grms, 8.4 % high**. On linear axes a power law with n = 2 or n = −2 is convex, so the chord between its end points lies above it; the error is largest on the long falling segment (the chord from 0.175 to 0.028 g²/Hz over 1200 Hz has area 122 g², the curve 84 g²).

**Check.** Dimensions: g²/Hz × Hz = g². The plateau alone is 131 g² = 11.5 Grms: the slopes add 29 % to the overall level here (14.76 against 11.46 Grms).

**Common mistakes.** Using n = 6 for +6 dB/oct; using the n ≠ −1 formula at −3 dB/oct (division by zero); integrating on the breakpoints with the trapezoid rule.

## 2. Derivations

(a) With G = G₁(f/f₁)ⁿ, ∫_{f₁}^{f₂} G df = G₁f₁⁻ⁿ [fⁿ⁺¹/(n + 1)]_{f₁}^{f₂} = **G₁f₁[(f₂/f₁)^{n+1} − 1]/(n + 1)** for n ≠ −1, and G₁f₁ [ln f]_{f₁}^{f₂} = **G₁f₁ ln(f₂/f₁)** for n = −1, which is a slope of −10 log₁₀2 = **−3.01 dB/oct**. (The first formula tends to the second as n → −1, by l'Hôpital.)

(b) In I₂ put r = 1/u, dr = −du/u². D(1/u) = (1 − 1/u²)² + 4ζ²/u² = [(u² − 1)² + 4ζ²u²]/u⁴ = D(u)/u⁴. So r²dr/D(r) = (1/u²)(−du/u²)(u⁴/D(u)) = −du/D(u), and the limits ∞ → 0 swap back: **I₂ = I₀**. Then, with f = rfₙ,

∫₀^∞ |H|² df = fₙ(I₀ + 4ζ²I₂) = fₙ(1 + 4ζ²) π/(4ζ) = **(π/2) fₙ Q (1 + 4ζ²)**, Q = 1/(2ζ).

For a white input of level G, Grms² = G × that, exactly. Miles' equation follows from three approximations: (i) G is flat near fₙ, so it can be taken out of the integral as G(fₙ); (ii) the integral can run from 0 to ∞ although the spectrum is band-limited; (iii) 1 + 4ζ² ≈ 1. Hence **Grms ≈ √((π/2) fₙ Q G(fₙ))**.

(c) Parseval's theorem for the transfer function F(ω) = 1/(ωₙ² − ω² + 2iζωₙω) and its impulse response h(t) = e^{−σt} sin(ω_d t)/ω_d (σ = ζωₙ, ω_d = ωₙ√(1 − ζ²)) gives ∫_{−∞}^{∞}|F|² dω = 2π∫₀^∞ h² dt. With sin² = (1 − cos 2ω_d t)/2,

∫₀^∞ e^{−2σt} sin²(ω_d t) dt = ½[1/(2σ) − σ/(2(σ² + ω_d²))] = ω_d²/(4σ(σ² + ω_d²)) = ω_d²/(4σωₙ²),

so ∫₀^∞ h² dt = 1/(4ζωₙ³) and ∫_{−∞}^{∞}|F|²dω = π/(2ζωₙ³). In r = ω/ωₙ, |F|² = 1/(ωₙ⁴D(r)) and dω = ωₙdr, so ∫_{−∞}^{∞} dr/D = π/(2ζ); the integrand is even, so **I₀ = π/(4ζ)**.

(d) |Z/Ÿ|² = |F|². With G in g²/Hz, z_rms² = g₀²G∫₀^∞|F|² df = g₀²G (1/(2π)) (1/2) π/(2ζωₙ³) = g₀²G/(8ζωₙ³) = **g₀²G/(64π³ζfₙ³)**. And g₀²Grms_Miles²/ωₙ⁴ = g₀²(π/2)fₙG/(2ζ)/(16π⁴fₙ⁴) = g₀²G/(64π³ζfₙ³): **the same**. (`relative_displacement_rms` agrees with it to 10⁻³ in the tests.)

**Common mistakes.** Forgetting the factor 2π between ∫dω and ∫df; integrating |H|² over the whole axis and forgetting the one-sided PSD covers only f ≥ 0.

## 3. A star-tracker bracket (s = 15: fₙ = 1050 Hz)

**Assumptions.** Linear SDOF, ζ = 0.03 (Q = 16.7), rigid base, the lesson's NASA-STD-7001C Table 2 workmanship spectrum.

(a) G(1050) = 0.04 × 500/1050 = **0.01905 g²/Hz** (on the −3 dB/oct segment). Miles: **22.88 Grms**; full integral: **23.36 Grms**; Miles is **2.0 % low**. On the falling segment the spectrum just below fₙ is higher than G(fₙ) and above it lower, which nearly cancels in the resonance, but the quasi-static region (|H|² ≈ 1 below the resonance, across the energetic plateau) adds area Miles does not have.

(b) **7.4 %** of the mean square comes from below 0.7fₙ. Miles discards it completely: it only knows G(fₙ).

(c) 3σ: **70.1 g**, bracket force 3 × 23.36 × 1.5 × 9.80665 = **1031 N**. ν₀⁺ = 1016 /s, so N = 1.22 × 10⁵ cycles in 120 s and Davenport's factor is **4.96**: expected largest peak **115.8 g**.

(d) z_rms = **5.25 µm**, 3σ **15.8 µm**, against a 0.10 mm clearance: a margin of more than 6 at 3σ and about 3.8 even at the 4.96σ expected peak (26 µm). **Not a concern**; a stiff bracket moves little, as 2.4 predicts (z falls as fₙ^(−3/2)).

**Check.** Miles' 22.88 = √(π/2 × 1050 × 16.67 × 0.01905): ✓.

**Common mistakes.** Reading G at the plateau value 0.04 instead of the slope value; using 3σ for the peak.

## 4. Identify the box from a SYNTHETIC test (s = 15)

(a) RMS / full integral: **0.9417** at T = 4 s (u = 3.64 %, **−1.60σ**), **1.0114** at 16 s (u = 1.82 %, **+0.62σ**), **1.0111** at 64 s (u = 0.91 %, **+1.22σ**). All within 2σ: consistent. The 4 s record is 5.8 % low and would look like a model error to anyone who had not computed u.

(b) From the 64 s SYNTHETIC records (`nperseg = 8192`, Δf = 1 Hz): least squares of |H| gives **fₙ = 150.000 ± 0.003 Hz, ζ = 0.05028 ± 0.00001** (covariance). The half-power points of |H|², interpolated between bins, are 141.79 and 157.12 Hz around the peak bin at 149 Hz: **ζ ≈ 15.33/(2 × 149) = 0.0514**. For reference, the exact half-power band of the model |H|² gives 0.05025 (the numerator 1 + (2ζr)² shifts it slightly from ζ).

(c) The least-squares fit, but not its error bar. ζ is recovered 0.6 % high, which is 25 "standard deviations" of the covariance estimate: the covariance accounts only for the scatter of the Welch estimate, not for its bias (the Hann window smears the 15 Hz-wide resonance over neighbouring 1 Hz bins). The half-power method is limited by the 1 Hz resolution: the peak bin is 1 Hz from the true peak and each crossing is interpolated across a 1 Hz step, a few per cent of a 15 Hz bandwidth. On a real test, longer segments (finer Δf) and a check that ζ does not change with `nperseg` are the honest way to bound the bias.

**Common mistakes.** Not discarding the start-up transient; fitting |H|² with a model of |H|; quoting the covariance error as the total uncertainty.

## 5. The test came in low

(a) (7.2/9.66)² = 0.556, i.e. **−2.55 dB** on the PSD. The tolerance is ±1.5 dB, which alone could bring the response to 8.13 Grms, not 7.2: **outside tolerance**, so the input is not the whole story (and the control accelerometer says it was within tolerance anyway).

(b) Solving the full integral: **Q = 5.52** at fₙ = 150 Hz (Miles' scaling gives 5.56), or **fₙ = 86 Hz** at Q = 10. A 43 % drop in frequency would be obvious in a sine survey and would mean a much softer mount than designed; a Q of 5.5 instead of 10 is entirely plausible.

(c) Ranked: (1) **amplitude-dependent damping**: bolted joints and cable bundles dissipate more at high level, so Q at 9 Grms is lower than any low-level estimate; test by comparing a low-level random or sine survey with the full-level transfer function. (2) **The accelerometer location**: on a box face rather than at the mass centre, or near a node of a higher mode, it under-reads the mode's motion; test by moving it or adding a second. (3) **The model's fₙ**: check the measured peak frequency; if it matches 150 Hz, this cause is out. (4) The input at the low edge of tolerance, which contributes but cannot explain all of it.

(d) **Do not lower the design load to the test value alone.** The test is one unit at one level on a rigid fixture; flight units differ, the flight interface differs, and a Q that fell with amplitude may not fall the same way in another unit. A reasonable decision: keep Q = 10 for the design load (or use the test-derived Q with a stated uncertainty factor), and use the test to justify that the margin exists. Lowering the load is defensible only with several units and a measured Q at full level.

**Common mistakes.** Treating ±1.5 dB as ±1.5 dB on the RMS (it is on the PSD; ±0.75 dB on the RMS); explaining the gap as "test variability" without computing what would explain it.

## 6. How narrow is the band? (s = 15)

(a) 5426 positive maxima in 32 s of SYNTHETIC response (170 per second). Fractions above kσ, measured against Rayleigh exp(−k²/2):

| level | measured | Rayleigh |
|---|---|---|
| 1σ | 0.541 | 0.607 |
| 2σ | 0.123 | 0.135 |
| 3σ | 0.0076 | 0.0111 |

The measured fractions are lower at every level.

(b) α = m₂/√(m₀m₄) = **0.778**, and the rate of maxima √(m₄/m₂) = **194 per second** against ν₀⁺ = 151 per second. A narrow-band process has α = 1 and one maximum per upcrossing. The *absolute* acceleration of the box contains, above resonance, the input passed through with gain 2ζfₙ/f; that broadband content dominates m₄ and adds small extra maxima, many below 1σ, so the fraction of maxima above any level falls. (The relative displacement, which falls as 1/f² above resonance, is much closer to narrow-band.)

(c) Rayleigh with ν₀⁺T = 9061 cycles in 60 s: about **101 cycles above 3σ**. One exceedance is expected at k = √(2 ln 9061) = **4.27σ**; Davenport's expected largest peak, 4.40σ, is that level plus the correction 0.5772/4.27 for the mean of the extreme-value distribution.

**Common mistakes.** Counting maxima of both signs as positive peaks; applying the Rayleigh law to a broadband signal without checking α.

## Check-yourself answers

1. g²/Hz (or (m/s²)²/Hz); the area is the mean square of the acceleration, g², and its square root is Grms.
2. n = −2; G falls to a quarter over one octave.
3. On linear axes a log-log straight segment with n ≠ 1 is a convex curve below its chord (for n > 1 or n < 0), so the trapezoid's chord overestimates the area; for the lesson's −6 dB/oct example the overall error is 8 %.
4. Grms = √((π/2) fₙ Q G(fₙ)). It assumes a linear SDOF with light viscous damping (the 4ζ² term dropped), an input PSD flat around fₙ, and a response dominated by the resonance (nothing from the quasi-static region or from beyond the band edges).
5. High on a rising slope (G(fₙ) overstates the input below fₙ) and low on a falling slope (the energetic plateau below fₙ adds quasi-static response that Miles omits); the errors grow with damping.
6. 0.27 %; in 60 s a 150 Hz response makes about 9000 cycles, so dozens exceed 3σ and the expected largest peak is about 4.4σ.
7. Its statistical scatter: for an SDOF response the RMS of a record of length T scatters by about 1/(2√(ζωₙT)), 1.3 % for 32 s of the lesson's box; 2 % is then 1.5σ, not a disagreement.
8. Acceleration × √2 (Grms ∝ √fₙ); displacement × 2^(−3/2) ≈ 0.354.
