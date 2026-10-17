# P08.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result with its uncertainty, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and probe scripts that call it, given to the figures shown. **Every signal analysed here is synthetic** (a seeded Gaussian random process with the von Kármán model spectrum), not a measurement.

## 1. Statistics of one short record (synthetic, seed 15)

**Assumptions.** Stationary record, Taylor's hypothesis (built into the synthetic signal), the model's known σ = 0.3 m/s and L = 0.05 m as the truth to recover.

(a) Mean **10.0000 m/s**; variance **0.08980 m²/s²**; turbulence intensity **3.00 %**; skewness **0.009**; flatness **3.053**. The band variance for 2¹⁶ samples is 0.08858 m²/s², so the record is **1.37 % high**; the expected relative standard deviation is √(2 × 0.866 × 0.005/3.277) = **5.1 %**, so the difference is **0.27 standard deviations**: entirely ordinary.

(b) Integral time scale **4.49 ms**, so **L = UT = 44.9 mm**. The spectral plateau with 2¹² -point segments gives **37.3 mm**. The spectral estimate is the biased one: a 2¹²-sample segment lasts 0.2 s, its bins are 4.9 Hz apart, and the first three bins (4.9–14.6 Hz) already sit on the roll-off of the model spectrum (whose corner, 1.339kL = 1, is at 23.8 Hz), while removing each segment's mean empties the lowest bin's neighbourhood. The autocorrelation estimate is unbiased to within its large scatter (Problem 4).

(c) Slope **−1.671** (Kolmogorov −5/3 = −1.667); ε = **0.3860 m²/s³** against the model's **0.3855** (+0.1 %). The inertial range has hundreds of independent spectral estimates, so ε is the best-determined statistic of the record.

**Check.** Skewness and flatness are 0 and 3 within their scatter: the synthetic signal is Gaussian by construction.

**Common mistakes.** Comparing the variance with σ² = 0.09 instead of the band variance (the 1.3 % above the Nyquist frequency is not in any sampled record); quoting L to three figures from one short record.

## 2. Derivations

(a) Averaging ∂u_i/∂t + ∂(u_i u_j)/∂x_j = −(1/ρ)∂p/∂x_i + ν∇²u_i with u_i = Ū_i + u′_i (⟨u′_i⟩ = 0, averages commute with derivatives): ⟨u_i u_j⟩ = Ū_iŪ_j + ⟨u′_i u′_j⟩, so
**∂Ū_i/∂t + ∂(Ū_iŪ_j)/∂x_j = −(1/ρ)∂P̄/∂x_i + ∂/∂x_j[ν∂Ū_i/∂x_j − ⟨u′_i u′_j⟩].**
Averaging continuity gives ∂Ū_j/∂x_j = 0, and subtracting it from the instantaneous one gives ∂u′_j/∂x_j = 0. The mean flow sees an extra stress **−ρ⟨u′_i u′_j⟩**. Closure problem: these six new unknowns have no equation of their own that does not introduce higher moments; a model must supply them (`P08.15`).

(b) σ̂² = (1/T)∫₀ᵀu′²dt. For a Gaussian process, Cov(u′²(t), u′²(t′)) = 2R(t − t′)², so Var(σ̂²) = (1/T²)∫∫2R²dt dt′ ≈ (2/T)∫_{−∞}^{∞}R(τ)²dτ = (2σ⁴/T)∫ρ²dτ for T ≫ T_int. The relative standard error is √(2∫ρ²dτ/T) = ε_r, so **T = 2∫ρ²dτ/ε_r²**. For ρ = e^(−|τ|/T_int), ∫e^(−2|τ|/T_int)dτ = T_int and **T = 2T_int/ε_r²**. For the von Kármán model ∫ρ²dτ = 0.866T_int (computed in a probe from ½∫S²df/σ⁴), so the exponential formula is 15 % conservative.

(c) ν^a ε^b has units m^(2a+2b)s^(−a−3b). Length: 2a + 2b = 1, a + 3b = 0 ⇒ a = 3/4, b = −1/4: **η = (ν³/ε)^(1/4)**. Time: 2a + 2b = 0, a + 3b = 1 ⇒ **τ_η = (ν/ε)^(1/2)**. Velocity: 2a + 2b = 1, a + 3b = 1 ⇒ **v_η = (νε)^(1/4)**. Then η v_η/ν = ν^(3/4+1/4)ε^(−1/4+1/4)/ν = **1**.

(d) The one-sided spectrum is S(f) = 4∫₀^∞R(τ)cos(2πfτ)dτ (Wiener–Khinchin), so S(0) = 4∫₀^∞R dτ = **4σ²T_int**. For the model, ∫₀^∞E11dk = (2σ²L/π)(1/(1.339L))∫₀^∞(1 + s²)^(−5/6)ds = (2σ²/(1.339π)) × 2.1033 = σ² (because 2.1033/1.339 = π/2). The wavenumber analogue L = πE11(0)/(2σ²) gives π(2σ²L/π)/(2σ²) = **L**.

**Common mistakes.** Forgetting the factor 2 in the Gaussian fourth moment; mixing one-sided and two-sided spectra (a factor of 2 in S(0)).

## 3. The light aeroplane in rough air (s = 15: σ = 1.0 m/s, L = 350 m)

**Assumptions.** Frozen turbulence convected past the aircraft at 60 m/s (the aircraft moves through it), the von Kármán model with an inertial range continuing down to η, ν at 3000 m.

(a) Tu = σ/U = **1.67 %**; T_int = L/U = **5.83 s**; energy-containing frequency U/L = **0.17 Hz**.

(b) ε = **2.04 × 10⁻³ m²/s³**; η = **1.33 mm**, τ_η = **0.096 s**; the smallest eddies pass the probe at U/(2πη) = **7.2 kHz**. Resolving them needs a sampling rate above twice that, in practice **20 kHz or more** with an anti-aliasing filter, and a probe smaller than about a millimetre. (Gust-load work needs only the large eddies, a few hertz.)

(c) With the exponential formula, ±10 %: **1167 s (70 km)**; ±5 %: **4667 s (280 km)**. With the von Kármán ρ (factor 0.866): **1010 s (61 km)** and **4041 s (242 km)**. Atmospheric turbulence is rarely stationary over 60–250 km, which is why airborne gust statistics carry wide error bars.

(d) Safer in flight: the probe moves at 60 m/s through gusts of σ = 1 m/s, **u′/U = 1.7 %**, against 3 % in the tunnel; Taylor's hypothesis needs u′/U ≪ 1.

**Common mistakes.** Using the air's speed relative to the ground; forgetting that the record length scales with L/U, so a fast aircraft needs fewer seconds but the same number of kilometres per integral scale.

## 4. Is my error bar right? (seeds 15–54)

(a)–(b) For 40 records of 2¹⁶ samples: mean variance 0.08757 m²/s² (band variance 0.08858); relative scatter **4.64 %** against the predicted **5.14 %**. The standard deviation of 40 samples is itself uncertain by about 1/√78 = 11 %, so 4.64 against 5.14 (−10 %) **agrees**.

(c) Autocorrelation L: **51.2 ± 9.8 mm**; spectral plateau (nperseg 2¹²): **40.9 ± 4.4 mm**; ε **0.3894 ± 0.0061 m²/s³** (model 0.3855). The spectral-plateau estimator is **biased about 18 % low** at this segment length (Problem 1b); the autocorrelation estimator is unbiased within its standard error (9.8/√40 = 1.5 mm) but has twice the scatter, because the integral to the first zero crossing of a noisy ρ is itself noisy.

(d) With 10 records of 2²⁰ samples: relative scatter **1.36 %** against the predicted **1.29 %**; the ratio of scatters is 4.64/1.36 = 3.4 against √16 = 4, consistent given that a standard deviation from 10 samples is uncertain by about 24 %. At this length the two L estimates are 52.8 ± 2.2 mm and 48.8 ± 2.2 mm (nperseg 2¹⁴): longer segments remove most of the spectral bias.

**Common mistakes.** Using ddof = 0 for a 40-sample standard deviation (a 1 % effect, harmless here, but know which you used); judging a single record against the truth and concluding the method is "wrong".

## 5. What the instrument does to the spectrum (synthetic, seed 2026)

(a) At U = 10 m/s and 20 kHz, ℓ/U is 2 samples for 1 mm and 6 for 3 mm.

| wire | variance ratio | slope (300–3000 Hz) | ε (m²/s³) | error in ε | effect on η |
|---|---|---|---|---|---|
| none | 1 | −1.658 | 0.3849 | — | — |
| 1 mm | 0.989 | −1.759 | 0.3408 | −11.5 % | +3.1 % |
| 3 mm | 0.959 | −3.270 | 0.1391 | −63.9 % | +29 % |

The variance hardly changes (the small eddies hold little energy), but the inertial range is cut off from the top: ε, inferred from the spectral level in the band, is under-estimated, and η = (ν³/ε)^(1/4) over-estimated. With 3 mm the fit band is no longer an inertial range at all (slope −3.3).

(b) Keeping every fourth sample (5 kHz, Nyquist 2.5 kHz) without a filter: slope **−1.461**, ε **0.483 m²/s³ (+25 %)**. At 2 kHz the decimated spectrum is **1.64 × 10⁻⁶ m²/s²/Hz**, against 1.12 × 10⁻⁶ at the full rate and 1.11 × 10⁻⁶ from the model: **+46 %**. Energy between 2.5 and 10 kHz has been folded back (aliased) onto the band below 2.5 kHz, flattening the spectrum near the Nyquist frequency and inflating ε.

(c) Reject it. The 1 mm wire spans 3.3η and already costs 11 % in ε; at 5 kHz the Nyquist frequency (2.5 kHz) is below the Kolmogorov frequency (5.2 kHz) and, without a filter, aliasing adds 25 % to ε. Use a wire no longer than about η (0.3 mm), sample at 20 kHz or more and low-pass filter below the Nyquist frequency.

**Common mistakes.** Filtering *after* decimating (the aliased energy is already in the band); reading the 3 mm slope as "steeper turbulence".

## 6. What a synthetic signal cannot tell you

(a) (1) **Intermittency**: real velocity increments over small separations have flatness far above 3 (rare, violent events), which a Gaussian signal cannot produce; for gust loads the tails of the distribution, not the variance, set the fatigue and the peak load. (2) **The skewness of the velocity derivative** (negative in real turbulence), the signature of the energy cascade's direction (vortex stretching transfers energy to small scales); a Gaussian signal has zero skewness at every scale and no cascade, so it says nothing about mixing rates in a combustor or the fine-scale structure that a model must reproduce.

(b) Plan: use the database's cutout service or Python client (outside this repository; large queries need a token) to extract u along one line of the forced isotropic turbulence field at one time; treat the line as a frozen record (x plays the role of Ut); run `moments`, `autocorrelation`, `welch_psd` and `fit_inertial_range` on it; compare E11(k) with the −5/3 law and with the synthetic record, and compute the velocity-increment flatness and derivative skewness that the synthetic record cannot show. Credit: cite the database as its Citing page requests (ODC-By requires attribution). It is **not a measurement**: it is direct numerical simulation data, a solution of the Navier–Stokes equations, and must be labelled as such.

(c) No: at 20 % intensity the eddies are convected at speeds that differ from the mean by ±20 % and evolve while passing the probe, so the time-to-space conversion x = Ut is no longer accurate, and the integral scale inferred from it is uncertain by a comparable amount.

## Check-yourself answers

1. u = Ū + u′ with Ū the time (or ensemble) average and ⟨u′⟩ = 0; averaging the momentum equation leaves the Reynolds stress −ρ⟨u′_i u′_j⟩ as an extra unknown stress on the mean flow.
2. T_int = ∫₀^∞ρ(τ)dτ, the time over which the signal remembers itself; L = UT_int. Integrating only to the first zero crossing avoids summing the noisy tail of the estimated ρ, whose random part does not average out.
3. Turbulence is assumed to be convected past the probe unchanged ("frozen") at the mean speed, so x = Ut; it is valid when u′/U ≪ 1 and the eddies evolve slowly compared with their passage time (3 % in the tunnel, 1.7 % in flight).
4. Between the energy-containing scales and the Kolmogorov scale, E11(k) = C1ε^(2/3)k^(−5/3) depends only on ε; fitting the spectral level in that range gives ε (0.385 m²/s³ for the lesson's record).
5. η = (ν³/ε)^(1/4), τ_η = (ν/ε)^(1/2), v_η = (νε)^(1/4): 0.307 mm, 6.2 ms, 0.049 m/s in the tunnel. L/η scales as Re_L^(3/4): 163 against 998^(3/4) ≈ 178; a turbulent flow's scale range grows rapidly with Reynolds number, which is why direct simulation is expensive.
6. It is built as a sum of independent Gaussian Fourier modes, so it is a Gaussian process: skewness 0, flatness 3. It has the right spectrum and the right two-point statistics, but no intermittency, no phase coupling and no cascade.
7. T = 2T_int/ε_r² = 2 × 5 ms/0.05² = 4 s (3.5 s with the von Kármán ρ).
8. A wire longer than η averages away the small eddies and biases ε low (−11.5 % for 1 mm); a sampling rate below twice the highest frequency present aliases energy into the measured band (+25 % in ε at 5 kHz).
