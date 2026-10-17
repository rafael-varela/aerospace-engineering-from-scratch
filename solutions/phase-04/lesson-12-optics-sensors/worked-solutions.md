# P04.12  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15; numerical results are from `solution.py` and the probe scripts that call it, given to the figures shown. Monte Carlo results are SYNTHETIC: computed on simulated images, not measured.

## 1. A fibre and a window

**Assumptions.** Step-index fibre, meridional rays, geometric optics; normal incidence on the window, no coating, multiple reflections ignored.

(a) n₂ = 1.455. **θ_c = arcsin(1.455/1.48) = 79.45°**.

(b) A ray entering the end face at θ_a refracts to θ_core with sin θ_a = n₁ sin θ_core, and meets the wall at 90° − θ_core. Trapping needs 90° − θ_core ≥ θ_c, so at the limit sin θ_a = n₁ sin(90° − θ_c) = n₁ cos θ_c = n₁√(1 − n₂²/n₁²) = **√(n₁² − n₂²)**. **NA = 0.2709, θ_a = 15.72°**. Check: `refraction_angle(1.0, 1.48, 15.72°)` gives θ_core = 10.55°, and 90° − 10.55° = 79.45° = θ_c.

(c) **R = (0.76/2.76)² = 7.58 % per surface**; two surfaces pass **(1 − R)² = 0.854**, a 14.6 % loss, nearly twice that of glass because sapphire's index is higher. An anti-reflection coating is worth adding.

**Common mistakes.** Using the critical angle itself as the acceptance angle (it is measured from the wall's normal, not the axis); forgetting the second surface.

## 2. Four derivations

(a) A lens focused at infinity images a star at angle α to x = f tan α. The detector's edge, at x = np/2, sees **FOV/2 = arctan(np/(2f))**. Near the axis dx = f dα, so one pixel subtends **p/f** radians. Since arctan u < u for u > 0, 2 arctan(np/2f) < np/f; with arctan u = u − u³/3 + …, the two agree to relative order (np/2f)²/3: for the baseline (np/2f = 0.113) the arctan FOV is 0.4 % smaller.

(b) Photon i is recorded at x_i = x₀ + e_i + r_i: e_i the scatter about the true centre (variance σ²), r_i the rounding to its pixel's centre. If the photons' positions are spread over many grid positions relative to the pixel boundaries (true when σ is not much smaller than a pixel, or when x₀ is random within the pixel), r_i is uniform on [−½, ½] with variance 1/12 and independent of e_i. The centroid is the mean of N independent x_i, so **Var = (σ² + 1/12)/N**. When σ ≪ 1 the r_i are all equal (every photon in one pixel), not independent: that is pixel locking.

(c) Left divider with the gauge at the bottom: V_ex R(1 + x)/(R(1 + x) + R) = V_ex(1 + x)/(2 + x); right divider: V_ex/2. **V_out = V_ex[(1 + x)/(2 + x) − ½] = V_ex x/(2(2 + x))**, and for x ≪ 2, **≈ V_ex x/4 = V_ex GF ε/4**. Full bridge: one divider holds R(1 − x) over R(1 + x), giving (1 + x)/2; the other R(1 + x) over R(1 − x), giving (1 − x)/2: **V_out = V_ex x exactly**, linear and four times the quarter bridge.

(d) Let k = γ/(γ − 1) and ε = (γ − 1)M²/2. (1 + ε)^k − 1 = kε + k(k − 1)ε²/2 + …; kε = γM²/2, and k(k − 1)ε²/2 = (γ/(γ − 1))(1/(γ − 1))(γ − 1)²M⁴/8 = γM⁴/8. So q_c = ½γpM²(1 + M²/4 + …), and since ρV² = ρa²M² = γpM², **q_c = ½ρV²(1 + M²/4 + …)**. Check: the exact ratio at M = 0.3 is 1.0227 against 1 + 0.0225; at 0.6, 1.093 against 1.09. Using Bernoulli at M = 0.6 overestimates V by √1.093 = 1.046, the 4.6 % of the lesson.

**Common mistakes.** In (b), forgetting the 1/12 (0.0429 instead of 0.0447 pixel for the baseline); in (c), putting both active gauges in the same arm, which cancels their signals.

## 3. A CubeSat battery thermistor

**Assumptions.** Beta model; ideal 12-bit ADC with 3.3 V reference; divider supply equal to the ADC reference (ratiometric).

(a) **R(20 °C) = 12 535 Ω**. dV/dT = V_s R_top/(R_top + R_T)² × dR_T/dT. For fixed R_T, g(R_top) = R_top/(R_top + R_T)² has g′ = (R_T − R_top)/(R_top + R_T)³ = 0 at **R_top = R_T = 12 535 Ω**, where V = 1.650 V, mid-scale.

(b) At 20 °C:

| R_top (Ω) | dV/dT (mV/K) | K per count |
|---|---|---|
| 2 000 | −18.0 | 0.0448 |
| 5 000 | −30.9 | 0.0261 |
| **12 535** | **−37.9** | **0.0212** |
| 20 000 | −35.9 | 0.0224 |
| 50 000 | −24.3 | 0.0331 |

The optimum is flat: a standard 12 kΩ resistor loses almost nothing.

(c) With 12 535 Ω:

| T (°C) | |dV/dT| (mV/K) | K per count |
|---|---|---|
| −20 | 19.3 | 0.0417 |
| 0 | 34.6 | 0.0233 |
| 40 | 27.8 | 0.0290 |
| 60 | 16.2 | **0.0497** |

Coarsest at +60 °C, where the thermistor's resistance is small compared with R_top; still 0.05 K, far finer than the thermistor's typical ±0.2–1 K interchangeability tolerance. The ADC is not the limit: the thermistor's calibration is.

(d) P = V²/R_T = 1.65²/12 535 = **0.217 mW**. Self-heating error = P/δ; for it to equal one count (0.0212 K), **δ = 10.2 mW/K**. Small bead thermistors have dissipation constants of order 1 mW/K in still air, and less in vacuum where there is no convection, so the self-heating error would be several counts to a few tenths of a kelvin: not negligible at one-count level. Fixes: power the divider only while sampling, or use a larger R_top at a small cost in sensitivity.

**Common mistakes.** Optimising the voltage instead of its slope; quoting resolution as accuracy.

## 4. Your own star tracker (s = 15)

**Assumptions.** As in the lesson: the assumed star-count fit and photometric zero point, σ = 1 pixel defocus, Poisson and read noise only.

(a) f = 23.5 mm, **D = 16.79 mm**; **scale 30.28″/px; FOV 17.10°, 0.0884 sr; 17.56 stars** expected to magnitude 5.5.

(b) **r_Airy = 1.22 × 0.6 μm × 1.4 = 1.025 μm = 0.297 pixel**. Yes: a focused spot of a third of a pixel would pixel-lock (the lesson's σ = 0.2 pixel case has 0.145 pixel bias); defocus to σ ≈ 1 pixel.

(c) **N_e = 279 electrons**. SYNTHETIC, seed 15, 2000 images:

| window | rms (px) | per star (″) | cross-boresight (″) |
|---|---|---|---|
| 3 × 3 | 0.157 | 4.76 | 1.14 |
| **5 × 5** | **0.1055** | **3.20** | **0.763** |
| 7 × 7 | 0.166 | 5.02 | 1.20 |

(d) **0.76″ against the baseline's 0.75″: essentially the same.** The finer pixels (30.3 against 45.4″, a factor 0.67) were cancelled by fewer photons: (16.79/20.83)² × ½ × 10^{0.2} = 0.515 of the baseline's electrons (smaller aperture, half the exposure, half a magnitude brighter limit), so the centroid in pixels grew from 0.068 to 0.106. The star counts are nearly equal (17.6 against 17.1): the wider field and the brighter limit cancel. To beat the baseline, this design should lengthen the exposure.

**Check.** Pixel scale by hand: 3.45 μm/23.5 mm = 146.8 μrad = 30.28″.

**Common mistakes.** Using the 3 × 3 window because it gives the "tightest" spot; forgetting that D changes with f at fixed f-number.

## 5. A wing-spar strain channel

**Assumptions.** Ideal bridge, GF = 2, 3.3 V excitation; amplifier gain exact; ADC 0–3.3 V, 12 bits, zero strain at mid-scale.

(a) x = 2 × 1500 × 10⁻⁶ = 0.003. **Quarter bridge: 2.4713 mV exact, 2.4750 mV linear, linearity error 0.15 %. Full bridge: 9.900 mV.**

(b) Gains for a 1.6 V swing: **full bridge 161.6, quarter bridge 647.4** (four times more, with four times the amplified offset and noise).

(c) The full-bridge channel's sensitivity is 1.6 V/1500 με; the 3.3 V ADC range spans 3.3/(1.6/1500) = 3094 με, so **0.755 με per count**. `bits_needed(3000e-6, 1e-6)` = **12** (2¹² = 4096 > 3000).

(d) Larger than one count in practice: (1) temperature: the gauge's resistance and the spar's thermal expansion change with temperature, tens of με per 10 K for an uncompensated gauge; the full bridge cancels it because all four gauges see the same temperature. (2) Amplifier offset and its drift: at a gain of 160 a 50 μV offset is 8 mV at the ADC, about 10 counts; use an instrumentation amplifier with low drift and zero the channel on the ground. (3) Gauge-factor and bonding uncertainty (about 1 %, 15 με at full scale): calibrate the spar with known loads. Also: lead-wire resistance (use a three- or four-wire connection) and electrical noise and aliasing (shielded twisted pair, the anti-aliasing filter of `P04.10`).

**Common mistakes.** Taking 3.3 V/4096 as the strain resolution without dividing by the gain and bridge sensitivity; claiming 0.76 με accuracy.

## 6. Reviewing a star-tracker proposal

Results below for the lesson's baseline tracker at magnitude 6 (543 electrons) and 10 e⁻ read noise, SYNTHETIC, seed 12, 2000 images.

**Focus sharply: wrong.** A focused spot (σ = 0.055 pixel) gives a high peak but pixel-locks: noise-free biases up to 0.37 pixel (17″) that more photons do not reduce. With a 9 × 9 window and 10 e⁻ the rms is **0.517 pixel = 23.5″ per star**, already above the 10″ requirement for one star. Change: defocus to σ ≈ 1 pixel.

**9 × 9 window: wrong.** With 10 e⁻ read noise each of the 81 pixels adds noise weighted by its distance from the centre; even after defocus the rms is **0.477 pixel (21.6″)**, against **0.153 pixel (6.95″)** with a 5 × 5 window. Change: a 5 × 5 window, or a thresholded window.

**One brightest star: wrong, for a different reason.** A single star fixes the direction of one line of sight, which constrains two axes; rotation about that line (roll) is unknown, so it cannot give a three-axis attitude at all. Its centroid is good (a magnitude-3 star, 8603 electrons, 5 × 5, 10 e⁻: 0.034 pixel = 1.5″), but there is no averaging and no protection against a misidentified star or a hot pixel. Change: identify and use all stars above threshold (about 17 to magnitude 6); with the defocused 5 × 5 centroid the cross-boresight error is 6.95″/√17.1 = **1.7″**, inside 10″ with margin for calibration errors, and roll is determined.

**First check with real images.** The read noise and the photon count: take night-sky frames through the actual lens and detector, measure the frame-to-frame scatter of the centroid of a magnitude-6 star with the chosen window (ground images include atmospheric seeing, so compare with the scatter of bright stars for scale). Reject the design if the measured per-star scatter with the defocused 5 × 5 window exceeds about 0.3 pixel (14″), which would leave less than a factor of two of margin on 10″ after averaging over the stars actually detected.

## Check yourself: answers

1. Above **41.81°** (arcsin(1/1.5)); **4 %** per surface at normal incidence.
2. At **d_i = f**, the focal plane, because a star is at infinity; the detector placed there records every star as a spot.
3. **p/f = 220 μrad = 45.4″ per pixel; FOV = 2 arctan(5.632/50) = 12.85°.**
4. A focused star (0.16 pixel) falls almost entirely on one pixel, so its centroid **locks** to the pixel centre with a systematic bias of up to 0.37 pixel (17″) that more photons do not reduce; defocusing to σ ≈ 1 pixel spreads it over several pixels and removes the bias.
5. **√(1 + 1/12)/√543 = 0.0447 pixel**; more than σ/√N = 0.0429 because each photon's position is also rounded to its pixel, adding 1/12 pixel² of variance per photon.
6. Every pixel in the window adds read noise whether or not it holds starlight, and the centroid weights it by its distance from the centre, so outer pixels add large errors and little signal.
7. About **V_ex GF ε/4 = 1.65 mV** (1.6484 mV exact). A full bridge gives **four times** the output (6.6 mV), is exactly linear, and cancels temperature effects common to all four gauges.
8. **32 g/4096 = 7.81 mg** per count; 1 mg over 32 g needs 32 000 levels, **15 bits**.
