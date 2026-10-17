# P06.13  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for **s = 15**; numerical results are from `solution.py` and the probe script that calls it, given to the figures shown. Results from `synthetic_tap_test` are SYNTHETIC; all responses are simulations, not measurements. Material values are typical, not design allowables.

## 1. Your own bracket (s = 15: m_tip = 0.70 kg, L = 95 mm)

**Assumptions.** Uniform cantilever, rigid tip mass, first mode only, Rayleigh shape.

(a) I = 40 × 5³/12 = **416.7 mm⁴** (4.167 × 10⁻¹⁰ m⁴); k = 3 × 68.9 × 10⁹ × 4.167 × 10⁻¹⁰/0.095³ = **1.0045 × 10⁵ N/m**; bracket mass 2700 × 0.04 × 0.005 × 0.095 = **51.3 g**; effective mass 0.70 + 0.2357 × 0.0513 = **0.71209 kg**.

(b) **f_n = 59.78 Hz**; with the tip mass alone 60.29 Hz, **0.86 %** high.

(c) Static deflection under 1 g: g/ω_n² = **0.0695 mm**; root stress (9): **3.98 MPa**. At resonance with ζ = 0.02 the deflection is magnified 1/(2ζ) = 25 times, about 99.5 MPa per g of input, so yield (276 MPa) is reached at **2.77 g** of base input.

**Check.** The bracket's mass is 7 % of the tip mass and 33/140 of that is 1.7 % of the total: half of that in frequency, 0.86 %, as found.

**Common mistakes.** Using the width instead of the thickness cubed in I; forgetting that the resonant magnification multiplies the 1 g stress.

## 2. The pieces of the lesson

(a) With ẇ(x, t) = ẇ_tip φ(x), φ = (3x²L − x³)/(2L³), the kinetic energy is T = ½ ẇ_tip² ρA ∫₀ᴸ φ² dx. With u = x/L, ∫₀¹ ((3u² − u³)/2)² du = ¼(9/5 − 1 + 1/7) = ¼ × 33/35 = **33/140**, so T = ½ (33/140) ρAL ẇ_tip² and the beam adds **(33/140) m_beam** to the tip mass. The method is an upper-bound estimate of the frequency (the true shape minimises the Rayleigh quotient), which is why the continuous beam of `P06.15` gives a slightly lower f_n.

(b) With y = Y e^{iωt}, z = Z e^{iωt}: m(z̈ + ÿ) + cż + kz = 0 ⇒ (k − mω² + icω)Z = mω²Y ⇒ **Z/Y = r²/(1 − r² + 2iζr)**, modulus (7). Then X = Y + Z: X/Y = (1 + 2iζr)/(1 − r² + 2iζr), modulus (8). At r = √2 the denominator's modulus² is (1 − 2)² + 8ζ² = 1 + 8ζ², and the numerator's is 1 + 8ζ²: **|X/Y| = 1 for every ζ**.

(c) At r = 1, Z = −A/(2iζω_n²) = iA/(2ζω_n²), so the steady part is Im(Z e^{iω_n t}) = (A/(2ζω_n²)) cos ω_n t, with z_s(0) = A/(2ζω_n²) and ż_s(0) = 0. The homogeneous part starts from the opposite values: z_h = −(A/(2ζω_n²)) e^{−ζω_n t}(cos ω_d t + (ζω_n/ω_d) sin ω_d t). For small ζ, ω_d ≈ ω_n and the sine term is O(ζ), so **z ≈ (A/(2ζω_n²))(1 − e^{−ζω_n t}) cos ω_n t**. 99 % when e^{−ζω_n t} = 0.01: **t = ln 100/(ζω_n)** (0.468 s for the lesson's bracket).

(d) In a half cycle with ẋ < 0 (moving left), friction pushes right: m ẍ + kx = +F, so x = d + (x_s − d) cos ω_n t with d = F/k, x_s the starting turning point. After half a period (cos = −1) it reaches 2d − x_s: the amplitude has fallen by 2d, and the half period is π/ω_n, unchanged. Two half cycles: **4F/k per cycle**. Energy: friction dissipates F × (distance) = F × 4X per cycle; a viscous damper dissipates ∮cẋ dx = πcωX². Equating at ω = ω_n: c_eq = 4F/(πω_nX), and ζ_eq = c_eq/(2mω_n) = 4F/(2πmω_n²X) = **2(F/k)/(πX)**.

**Common mistakes.** In (b), writing the inertia term with ẍ − ÿ; in (c), forgetting that the homogeneous part must cancel *both* initial values.

## 3. Designing a vibration sensor

**Assumptions.** Linear SDOF sensor, rigidly attached case, 5 % flatness band.

(a) Below resonance (r ≪ 1) the denominator of Z = −A/(ω_n² − ω² + 2iζω_nω) is ≈ ω_n², so **z ≈ −ÿ/ω_n²**. Solving |1/√((1 − r²)² + (2ζr)²) − 1| = 0.05 for the largest r that keeps the gain inside ±5 % all the way down:

| ζ | largest r | f_n for flat to 2 kHz | phase lag at that r |
|---|---|---|---|
| 0.02 | 0.218 | 9.16 kHz | 0.5° |
| 0.6 | 0.844 | 2.37 kHz | 74° |
| 0.7 | 0.591 | 3.38 kHz | 52° |

(b) Above resonance (r ≫ 1), Z/Y = r²/(1 − r² + 2iζr) → −1, so **z ≈ −y**: the mass stands still and the case moves around it. |Z/Y| within 5 % of 1 for all r above **4.58** at ζ = 0.02 and **1.69** at ζ = 0.7.

(c) **ζ ≈ 0.6–0.7.** It widens the flat band (r up to 0.6–0.84 instead of 0.22) and, more importantly for a transient, makes the phase lag nearly proportional to frequency (a pure time delay), so a broadband pulse is delayed but not distorted. A lightly damped accelerometer would ring at its own resonance whenever the transient contained energy near it. (An undamped piezoelectric accelerometer avoids the problem differently: its f_n is so high, tens of kHz, that the band of interest sits at r ≪ 0.2.)

**Check.** At ζ = 1/√2 the gain is 1/√(1 + r⁴), maximally flat; 0.6 trades a small hump for more bandwidth.

**Common mistakes.** Reading the seismometer's requirement off the absolute transmissibility (8) instead of the relative one (7).

## 4. A SYNTHETIC tap test and a test prediction (s = 15)

(a) SYNTHETIC record, seed 15: ζ = **0.01698** (n = 1), **0.01880** (5), **0.01964** (10), **0.01913** (15), **0.01874** (20); damped frequency from 10 peak intervals **78.57 Hz** (exact 78.26 Hz; peak timing on a noisy record is the weak point, and a zero-crossing or FFT estimate would do better). Trust **n = 10**: experiment 1 showed the scatter at n = 1–5 and the noise bias beyond about 10 cycles, where the peaks fall below ten times the 1 μm noise.

(b) With ζ = 0.01964: camera **51.0 g**; deflection A/(2ζ√(1 − ζ²)ω_n²) = **2.065 mm**; root stress **128.1 MPa**.

(c) A 1 s resonant dwell from rest reaches max|z| = **2.0645 mm**, **99.97 %** of the steady value (ζω_n t = 9.7 time constants; the first-order envelope estimate gives 99.99 %).

(d) Stress ∝ 1/ζ: ζ = 0.01964 ± 0.0004 gives **125.5–130.7 MPa**.

**Common mistakes.** Using n = 1 because "it is the definition"; quoting the ζ that generated the record (0.02) instead of the one you estimated.

## 5. Friction and viscous damping together

(a) The integrator of `simulate_coulomb` with −2ζω_n ẋ added; the decay from 100 μm stops after about 14 cycles.

(b) Peaks (from the first full cycle on): 90.03, 80.67, 71.88, 63.63, 55.88, 48.60, 41.76, 35.34, 29.31, 23.65, 18.33, 13.34 μm, …. A straight line leaves residuals of **2.96 μm** (standard deviation) and an exponential **0.245** in ln x (25 %). **Neither is adequate**: the decrements shrink (viscous) but not in proportion to the amplitude (friction).

(c) In one half cycle the motion is a decaying oscillation about the friction offset ±d: the distance from the centre is multiplied by √q with √q = e^{−πζ/√(1 − ζ²)}. Starting at x₀ > 0 (centre +d), the next turning point is at distance √q(x₀ − d) beyond the centre on the other side: |x_{½}| = √q(x₀ − d) − d. The next half cycle (centre −d) gives x₁ = √q(|x_{½}| − d) − d = √q(√q(x₀ − d) − 2d) − d = q x₀ − d(q + 2√q + 1), so **x_{n+1} = q x_n − d(1 + √q)²**. Regression of x_{n+1} on x_n gives **q = 0.93910** (ζ = **0.010000**) and **D = 3.877 μm**; d(1 + √q)² = 1 × (1 + 0.96907)² = **3.877 μm**. The two mechanisms are separated cleanly.

(d) The motion stops at the first turning point inside |x| ≤ d = 1 μm, after about 14 cycles (final position 0.11 μm).

**Common mistakes.** Fitting only one of the two models and reading a damping value from it; forgetting that friction's offset changes sign each half cycle.

## 6. Before the test

(a) Response ∝ input, so 30 g at the camera needs **2 g × 30/50 = 1.2 g** at the base near 78 Hz. The risk: a notch placed by analysis alone relies on f_n and ζ from a model and a tap test; if f_n is 5 % off, the notch misses the peak (whose half-power band is only 3 Hz wide) and the camera sees the full 50 g. Standard practice is a low-level sweep first to locate the resonance and measure its Q, then notch on measured data.

(b) **Lower:** joint friction's equivalent damping falls with amplitude, (12). **Higher:** at larger amplitude more of the joint microslips, adding dissipation; viscoelastic or air damping can also grow. The tap value must be checked at test level.

(c) Stress concentrations at the bolt holes and the root fillet (factors of 2–3), and the bolt preload's bending and contact stresses; also the bracket's second mode and any rocking of the camera, which (9) does not include.

(d) Probably, to about 1 % in frequency here (Problem 1: the bracket mass changes f_n by under 1 %). To find out, solve the cantilever-with-tip-mass frequency equation of `P06.15` (or a beam FEM) and compare; the Rayleigh estimate is an upper bound, so the true f_n is a little lower.

## Check-yourself answers

1. Points along the bracket move less than the tip; weighting each piece by the square of its velocity in the deflection shape and integrating gives 33/140 of the bracket's mass at the tip velocity.
2. δ = (1/n) ln(x₀/x_n), ζ = δ/√(4π² + δ²). Late peaks are small and noisy: noise biases them upwards and ζ comes out low, so use only cycles with a good signal-to-noise ratio.
3. Receptance: 1/k, 1/(2ζk) at −90°, falling as 1/(mω²). Accelerance: rising as ω²/k from zero, peaking near ω_n, tending to 1/m (the mass line).
4. m z̈ + c ż + k z = −m ÿ with z = x − y: the base acceleration acts as an inertial force −mÿ on the mass.
5. r = √2. Above it, more damping transmits more motion (the damper becomes a second load path).
6. The relative transmissibility (deflection) sets the stress; the absolute transmissibility sets the camera's acceleration and the force into the base.
7. ln 100/(ζω_n), about 37 cycles at ζ = 0.02 (0.47 s at 78 Hz); a sweep that crosses the half-power band faster than that does not reach the steady peak, so a fast sweep under-tests (and an analysis that assumes steady state over-predicts).
8. Coulomb peaks fall by a constant difference (a straight line in linear axes) at the undamped frequency and stop dead; viscous peaks fall by a constant ratio (a straight line in log axes). Friction's equivalent ζ is 2(F/k)/(πX): inversely proportional to amplitude.
