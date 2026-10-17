# P04.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check, and the common mistakes. Seeded problems are worked for $s = 15$; numbers come from `solution.py` and probe calculations that call it. Results from `make_synthetic_pendulum_data` are **SYNTHETIC** and labelled as such.

## 1. Planning a pendulum measurement ($s = 15$: L = 0.80 m)

**Assumptions.** Point bob, rigid support, g = 9.80665 m/s², frame timing as in `P00.15`.

(a) T₀ = 2π√(0.8/9.80665) = **1.7946 s**. Exact periods: 10° → excess **0.191 %**; 20° → **1.8083 s**, excess **0.767 %** (the 10° period is 1.7980 s).

(b) One timing: u = √2 × 0.5/(√3 × 30) = **13.6 ms**. Per period:

| N swings | 1 | 5 | 10 | 20 | 50 |
|---|---|---|---|---|---|
| u_T (ms) | 13.6 | 2.72 | 1.36 | 0.68 | 0.27 |
| u_T/T₀ | 0.76 % | 0.15 % | 0.076 % | 0.038 % | 0.015 % |

(c) The 10° correction is 0.191 % of 1.7946 s = 3.42 ms; one tenth is 0.342 ms; N = 13.6/0.342 = 39.8, so **N = 40 swings**. It is possible but costly: over 40 swings (72 s) a light bob loses a noticeable fraction of its amplitude, so the "10°" run is really an average over perhaps 10° down to 7°, and the predicted correction must use the mean amplitude. Five repeats of 10 swings achieve a Type A reduction too, but not of the frame term, which is the same for each run.

**Check.** The excess at 10° from the series θ₀²/16 = (0.1745)²/16 = 0.190 %, matching the exact 0.191 %.

**Common mistakes.** Dividing the frame term by √N instead of N (it is one timing of N periods, not N timings); expecting repeats to shrink the frame term.

## 2. Derivations

(a) dE/dt = mẋẍ + kxẋ = ẋ(mẍ + kx) = 0 by the equation of motion. At a turning point ẋ = 0 and |x| = A, so **E = ½kA²**.

(b) For ζ < 1, x = Ae^{−ζω₀t}cos(ω_d t − φ). Its maxima are separated by T_d = 2π/ω_d (the cosine factor repeats; the envelope only rescales). So x_n/x_{n+1} = e^{ζω₀T_d} = **exp(2πζ/√(1 − ζ²))**, and δ = 2πζ/√(1 − ζ²). Squaring: δ²(1 − ζ²) = 4π²ζ², so **ζ = δ/√(4π² + δ²)**. Check: the code's simulated peaks at ζ = 0.05 give δ = 0.3146 and ζ = 0.0500.

(c) X² = f₀²/((ω₀² − ω²)² + 4ζ²ω₀²ω²). With ω = ω₀(1 + ε), ω₀² − ω² = −ω₀²(2ε + ε²) ≈ −2εω₀², and 4ζ²ω₀²ω² ≈ 4ζ²ω₀⁴ to first order. So X² ≈ f₀²/(4ω₀⁴(ε² + ζ²)), maximal at ε = 0 and half of it at **ε = ±ζ**. Bandwidth **2ζω₀ = ω₀/Q**. Check: the code's grid search at ζ = 0.05 gives the band within 2 % of 2ζω₀ (a test).

(d) Energy (mL² about the pivot): ½mL²θ̇² − mgL cos θ = −mgL cos θ₀, so θ̇ = √(2g/L)√(cos θ − cos θ₀). The time from θ₀ to 0 is a quarter period: **T = 4√(L/g) ∫₀^{θ₀} dθ/√(2(cos θ − cos θ₀))**. With cos θ = 1 − 2 sin²(θ/2): 2(cos θ − cos θ₀) = 4(k² − sin²(θ/2)), k = sin(θ₀/2). Substituting sin(θ/2) = k sin φ gives ½cos(θ/2) dθ = k cos φ dφ and k² − sin²(θ/2) = k² cos² φ, so the integrand becomes dφ/√(1 − k² sin² φ) and the limits 0 → π/2: **T = 4√(L/g)K(k²)**. To first order (1 − k² sin² φ)^{−1/2} ≈ 1 + ½k² sin² φ; ∫₀^{π/2} = π/2 + (π/8)k². So T ≈ 2π√(L/g)(1 + k²/4) and with k² ≈ θ₀²/4, **T ≈ T₀(1 + θ₀²/16)**.

(e) An element dx of string with tension T_s at slope angle α at each end: the vertical force is T_s(sin α(x + dx) − sin α(x)). **Small slopes:** sin α ≈ tan α = u_x, and the horizontal components cancel (so the tension is uniform). Then F_y = T_s(u_x(x + dx) − u_x(x)) = T_s u_xx dx, and Newton: μ dx u_tt = T_s u_xx dx, so **u_tt = (T_s/μ) u_xx**. Units: √(N/(kg/m)) = √(kg m s⁻² m/kg) = m/s.

**Common mistakes.** In (b), measuring the decrement between a peak and the next trough (half a period: that gives half of δ); in (d), forgetting the factor 4 or using L instead of L² in the moment of inertia; in (e), keeping sin α instead of the slope when writing the force.

## 3. Solar-array boom and reaction wheel ($s = 15$: ζ = 0.010)

**Assumptions.** Single-degree-of-freedom boom, viscous damping, harmonic forcing from the wheel imbalance at the spin frequency, linear response.

(a) **Q = 50**. Peak at f₀√(1 − 2ζ²) = **0.79992 Hz**; amplification **50.0** times static.

(b) Response above 10 × static between **0.7597 Hz and 0.8382 Hz**, i.e. wheel speeds **45.6 to 50.3 rpm** (a band of 0.079 Hz, about −5.0 %/+4.8 % around f₀). The half-power bandwidth 2ζf₀ = **0.016 Hz** is five times narrower: the ×10 band is much wider than the "resonance width" because amplification falls off slowly in the skirts.

(c) Envelope e^{−ζω₀t} = 0.01 gives t = ln 100/(ζω₀) = **91.6 s**, **73.3 cycles** (log decrement 0.0628 per cycle). A boom this lightly damped rings for a minute and a half after every slew.

(d) At 2f₀: **0.333**, attenuated (above resonance the response falls as 1/ω²). At f₀/2: **1.333** (below resonance the boom follows the force, slightly amplified). **Rule:** keep the wheel's operating speeds, including its zero-crossings during manoeuvres and its harmonics, outside roughly 40–55 rpm, and either pass through the band quickly or add damping to the boom.

**Check.** Peak amplification ≈ Q for small ζ (50.0025 vs 50); decay time × f₀ × δ = ln 100 = 4.61 (73.3 × 0.0628 = 4.60).

**Common mistakes.** Treating the half-power band as the danger band; converting Hz to rpm by × 2π instead of × 60.

## 4. The wave solver ($s = 15$: T_s = 65 N)

(a) c = √(65/0.002) = **180.3 m/s**; f₁, f₂, f₃ = **138.7, 277.4, 416.0 Hz**.

(b)

| Courant | max error vs d'Alembert | energy spread |
|---|---|---|
| 1.0 | 7.5 × 10⁻¹⁵ | 5.0 × 10⁻¹⁵ |
| 0.8 | 2.84 × 10⁻² | 3.6 × 10⁻¹⁵ |
| 0.5 | 5.67 × 10⁻² | 5.2 × 10⁻¹⁵ |

At C = 1 the scheme is exact on the grid (each update is an identity for any sampled f(x − ct) + g(x + ct)), so only rounding remains. Below 1 it is **dispersive**: the short wavelengths in the 3 cm Gaussian (only about 9 cells across its half-width) travel slower than c and trail behind as ripples, and the error grows as C falls. The discrete energy is conserved to rounding at every C ≤ 1: dispersion moves energy between wavelengths, it does not create or destroy it.

(c) Mode 3 at C = 0.5, t = 1.3 periods:

| cells | 50 | 100 | 200 | 400 |
|---|---|---|---|---|
| max error | 8.66 × 10⁻³ | 2.16 × 10⁻³ | 5.40 × 10⁻⁴ | 1.35 × 10⁻⁴ |

**Observed order 2.00.** (Ending at exactly one period hides the phase error and gives a spurious higher order; that is why the run ends at 1.3 periods.)

(d) C = 1.02: max |u| = **3.5 × 10⁴⁷**; C = 1.05: **2.1 × 10⁸⁷**. The physical solution is bounded, but the scheme's amplification factor for the shortest grid wavelength exceeds 1 when C > 1; rounding seeds that mode at 10⁻¹⁶, and about 2 × 200/1.02 ≈ 392 steps per round trip multiply it beyond any bound. Instability is a property of the discrete map, not of the physics.

**Common mistakes.** Calling the C = 0.8 error a "bug" (it is dispersion, and it converges with the grid); measuring the order at a whole period.

## 5. SYNTHETIC pendulum dataset ($s = 15$, seed 15)

All numbers in this solution come from **SYNTHETIC** data generated by `make_synthetic_pendulum_data(seed=15)`.

(a) and (b) (u_T includes the frame term 13.6 ms/10):

| L (m) | T (s) | u_T (ms) | g (m/s²) | u_g | g corrected |
|---|---|---|---|---|---|
| 0.40 | 1.27133 | 1.59 | 9.770 | 0.055 | 9.798 |
| 0.60 | 1.55667 | 1.36 | 9.775 | 0.037 | 9.803 |
| 0.80 | 1.79667 | 1.72 | 9.784 | 0.031 | 9.812 |
| 1.00 | 2.00800 | 1.59 | 9.791 | 0.025 | 9.819 |
| 1.20 | 2.20067 | 1.52 | 9.782 | 0.021 | 9.810 |
| 1.50 | 2.46200 | 1.59 | 9.770 | 0.018 | 9.797 |

Correction factor (T/T₀)² at 0.15 rad = 1.00282. Weighted mean of the corrected values: **9.806 ± 0.011 m/s²**.

(c) p = **0.4998**; g from A = **9.779 m/s²** uncorrected, **9.807 m/s²** corrected.

(d) Corrected: E_n = |9.806 − 9.81|/0.011 = **0.36**, consistent. Uncorrected: weighted mean 9.779 ± 0.011, **E_n = 2.97**: you would have concluded, wrongly, that the data disagree with the local g by three standard uncertainties and gone looking for a timing or length error. The 0.28 % amplitude effect at 8.6° (0.027 m/s²) is about one u_g at L = 1 m, so no single length shows it, but the mean of six has u = 0.011 m/s² and it does: averaging sharpens the test until a small systematic effect dominates.

**Common mistakes.** Forgetting that the data are synthetic and writing "we measured"; omitting the frame term; applying the amplitude correction to T as a multiplication instead of a division (g then moves the wrong way).

## 6. Judgment

(a) **Air drag** is no longer small against the weight of a 10 g paper cone: the amplitude decays within a few swings, and drag lowers the frequency slightly and makes the decay non-exponential. Use a dense, compact bob; the model's "point mass in vacuum" needs a bob whose terminal speed is far above its swing speed.

(b) **The rigid support.** A curtain rail flexes and rocks, coupling the pendulum to the rail's own mode: the period shifts and the amplitude may beat (exchange energy with the rail). Test by touching the rail during a swing, or by comparing with a rigid support; fix by clamping to a door frame.

(c) **Viscous (linear) damping.** A decrement that grows with amplitude means the damping is nonlinear (friction in the hinge or joints, which dissipates relatively more at large amplitude). The resonance amplification at small amplitude (ζ ≈ 0.003, Q ≈ 160) is then much larger than a ζ measured at large amplitude suggests; design with the small-amplitude value.

(d) **The numerical model, not the physics.** At Courant 0.8 a sharp front disperses: ripples trail the shock and its arrival time smears. Run at C = 1 on a uniform cable (exact), or use a scheme designed for sharp fronts (`P03.12` finite volumes); check by halving Δx.

(e) **Constant structural damping.** Aerodynamic forces feed energy into the wing in a way that rises with speed; the total damping heads toward zero and, beyond it, negative: flutter. The trend in ζ against speed is the flight-test warning sign; extrapolating it to ζ = 0 gives the flutter speed, which must stay outside the flight envelope with margin.

## Check-yourself answers

1. T₀ = 2π√(L/g). It does not depend on the bob's mass or (in the small-angle limit) on the amplitude.
2. ζ = c/(2√(km)) = c/(2mω₀); Q = 1/(2ζ). At resonance the amplitude is Q times the static deflection.
3. δ = ln 1.37 = 0.315; ζ = δ/√(4π² + δ²) = 0.0500.
4. 1.74 % at 30°; 0.048 % at 5°.
5. c = √(T_s/μ); f_n = nc/(2L), n = 1, 2, 3, …
6. C = cΔt/Δx ≤ 1. At exactly 1 the scheme reproduces the exact solution on the grid (d'Alembert), apart from rounding.
7. A fixed end forces u = 0 there; the reflected wave must cancel the incident one at the boundary, which requires an inverted pulse (mathematically, the odd extension of the profile about the end).
8. u = 0.41/30 = 13.6 ms per timing, 1.36 ms per period over ten swings. No: repeating the run reduces only the scatter (Type A); the frame term is the same quantisation on every run.
