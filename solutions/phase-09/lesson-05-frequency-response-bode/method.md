# P09.05  Method

The route to each answer, without the final numbers. `worked-solutions.md` has the numbers for s = 15.

## 1. Servo loop

(a) K_B = K, type 1, real pole corner 1/τ = 10 rad/s (`bode_form`). (b) Straight line: 20 log₁₀ K − 20 log₁₀ ω below 10 rad/s, −40 dB/decade above; phase −90° below 1 rad/s, −180° above 100 rad/s, linear in log ω between. (c) `bode(G, [1, 10, 100])`; the differences are the real-pole corrections (−0.043, −3.01, −0.043 dB; 5.7°, 0°, 5.7°). (d) ω_c² (1 + ω_c²τ²) = K²; `gain_crossover`; the straight line on the −40 slope crosses at √(K/τ).

## 2. Derivations

(a) y = Re{Ae^{jωt}[G(jω) − ∫_t^∞ h(τ)e^{−jωτ}dτ]}; bound the tail by ∫_t^∞|h| ∝ e^{−σ_min t}. (b) |G₁G₂| = |G₁||G₂|, log of a product; arg of a product. (c) 20 log₁₀(1/√2), 20 log₁₀(1/√1.01), the decade rule's line versus −arctan(ω/ω_b) at ω = ω_b/10 and 10ω_b. (d) d/du² of (1 − u²)² + 4ζ²u² = 0 ⇒ u² = 1 − 2ζ²; substitute; at u = 1, |G| = 1/(2ζ) while the line is 0 dB.

## 3. Notch

(a) 20 log₁₀ M_r with ζ = 0.02 (`resonance`). (b) `bode(notch, 2π f)` at each f. (c) mode × notch with `np.polymul`; `bode` on `linspace(2π·1, 2π·40, 200001)`; `argmax`. (d) The notch's phase lag at 1 Hz is taken from the pitch loop's margin; a mode shifted 10 % sits at the edge of the notch, where the attenuation is far smaller than at its centre: compute it.

## 4. Sine-test settling

(a) `freq_response` at ω. (b) `sine_test(G, w, settle=…, n_periods=3)` minus the exact values. (c) The leftover transient after t_settle is about e^{−1.068 t_settle} of its initial size. (d) e^{−1.068 t} small enough for 0.01° (use the table to calibrate), plus a few periods at each frequency (2π/0.1 s per period at the low end); add over 20 log-spaced frequencies.

## 5. Judgment

(a) Shape, slopes, rough crossover: yes; margins and lag at a frequency: no. (b) `bode` of (s + 2)/… and (−s + 2)/…; same magnitude, the RHP version lags by 2 arctan(ω/2) more. (c) Backlash/dead zone and static friction (stiction) reduce the response at small amplitude; repeat at ±0.5 and ±2 deg and see whether the deficit shrinks with amplitude. (d) No steady state exists for an unstable open loop: the transient grows; test with a stabilising loop closed and divide out the controller.

## 6. Time delay

(a) e^{−jωT} = cos ωT − j sin ωT. Magnitude unchanged ⇒ ω_c unchanged. (b) `gain_crossover(L)`; subtract ω_c T in degrees. (c) −2 arctan(ωT/2) against −ωT (degrees); `phase_from_factors` of (1 − sT/2)/(1 + sT/2) at ω_c. (d) T = (180° + ∠L(jω_c)) in radians / ω_c.
