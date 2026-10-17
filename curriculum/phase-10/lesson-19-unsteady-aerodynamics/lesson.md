# P10.19  Unsteady Aerodynamics

> Every change of lift leaves vorticity in the wake, and the wake answers back: lift lags, shrinks and, at high frequency, leads. Theodorsen in frequency, Wagner in time, and the Fourier transform between them.

| | |
|---|---|
| **Learning outcomes** | 1. Define the reduced frequency k = ωb/V and reduced time s = Vt/b, and decide from them whether a problem is quasi-steady. 2. Derive the structure of Theodorsen's lift and moment for a pitching and plunging thin airfoil (apparent-mass terms plus a circulatory term scaled by C(k)), and its quasi-steady limit 2πα₃/₄. 3. Derive, with the Laplace transform, the C(k) implied by an exponential indicial function, and use Duhamel's integral to compute the lift for any motion. 4. Implement and test Theodorsen's function from Hankel functions, Jones's Wagner approximation, Sears's Küssner approximation, a second-order Duhamel quadrature, lift loops and aerodynamic work, and check them against each other. 5. Judge where the attached-flow theory stops (dynamic stall, finite span, compressibility) and what it does and does not say about flutter. |
| **Prerequisites** | `P02.08` (Laplace transform, transfer functions, convolution), `P02.09` (Fourier series and the complex exponential), `P02.12` (complex variables), `P10.04` (thin-airfoil theory: lift at the quarter chord, the three-quarter-chord point), `P10.03` (circulation and the Kutta condition), `P03.05` (integrators and observed order). Skill check: you can write the response of a linear system to e^{iωt} as a complex amplitude and read off its gain and phase. |
| **Engineering question** | The reference trainer's wing (chord 1.2 m, cruise 55 m/s) twists in a 5 Hz torsional vibration about its quarter chord. How large is the lift oscillation compared with the quasi-steady 2πα, how far is it out of phase, and how many metres must the aircraft fly after a sudden change of angle of attack before the lift reaches 90 % of its new value? The answers: k = 0.343; the lift amplitude is 0.711 × 2πα̂ (|cl|/α̂ = 4.46 per radian), leading the pitch by 18.0° (10.0 ms); the wake alone cuts the circulatory lift to |C| = 0.669 and delays it by 15.0°; after a step, 90 % of the lift arrives after 12.7 semi-chords, 7.6 m or 139 ms. |
| **Assumptions and limits** | Two-dimensional thin airfoil (a flat plate), incompressible, inviscid, attached flow; small motions about a straight path; a flat wake that convects at the free-stream speed; the Kutta condition holds instantaneously. It fails with flow separation (dynamic stall), at Mach numbers where compressibility delays the pressure response, for low aspect ratio (the trailing vortices add their own lag), and for rotors whose wake returns under the blade. |
| **Analytical work** | Reduced frequency and time; the apparent-mass and circulatory lift and moment of Theodorsen's theory in coefficient form; the quasi-steady limit; the Laplace-transform route from φ(s) to C(k) and its limits; Duhamel's integral; the area of a lift loop and the work per cycle; the sign of the pitch damping about different pivots. |
| **Computational deliverable** | `starter/main.py`: `theodorsen`, `theodorsen_from_indicial`, `theodorsen_lift`, `theodorsen_moment`, `quasi_steady_lift`, `aero_work_per_cycle`, `wagner_jones`, `kussner_sears`, `duhamel`, `time_domain_lift`, `simulate_harmonic_pitch`, `rise_distance`, `gust_response`, `lift_loop`, `pitch_flutter_boundary` (`fit_harmonic`, `loop_area`, `observed_order`, `reduced_frequency` and `main` are given). Tests (17) prove: C(k) → 1 as k → 0 and → ½ as k → ∞, with F and G in range; the lift tends to 2πα₃/₄ as k → 0 for pitch and plunge and any pivot; linearity; the circulatory moment vanishes about the quarter chord; φ(0) = ½, ψ(0) = 0, both → 1; Jones's C(k) is within 0.02 of the exact one, and equals a direct numerical transform to 2 × 10⁻⁵; Duhamel converges at order 2.0 to a closed-form convolution; the time-domain pitch matches Theodorsen after the transient; a steady descent equals an angle of attack; lift-loop area equals its harmonic formula; plunge is always damped; pitch-only instability needs a forward pivot and low k; the gust response tends to its quasi-steady value. |
| **Practical activity** | **computational experiment**: C(k) and the accuracy of Jones's approximation; a time-domain simulation of harmonic pitch against Theodorsen across k; step and gust responses of the trainer; lift loops; the map of pitch damping against pivot position and k. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `unsteady_aero.png`, the derivation of C(k) from the indicial function, the five experiment tables, and a half-page note answering the engineering question with the assumption most likely to fail for the trainer's wing. |
| **Assessment** | `assessments/phase-10/grading-rubric.md`; competencies 1–4. Pass: the Laplace-transform route to C(k) is derived; all public tests pass; the time-domain and frequency-domain results are compared with the error attributed (discretisation against Jones's fit); the note states amplitude, phase and rise distance with their limits. |
| **References** | `@ntrs-naca-tr-496-theodorsen-flutter` Theodorsen, NACA TR-496, "Velocity Potentials, Forces, and Moments of the Noncirculatory Flow"; "Velocity Potentials, Forces, and Moments of the Circulatory Flow"; "The Magnitude of the Circulation"; "Value of the Function C(k)"; "Total Aerodynamic Forces and Moments" (supplementary-free; read the page images, the OCR is unreliable): the original derivation this lesson's formulas follow. `@ntrs-naca-tn-667-jones-operational-nonuniform-lift` R. T. Jones, NACA TN-667, Introduction (supplementary-free): lift growth after a sudden change by the method of operators, the source of the exponential approximation. `@leishman-intro-aerospace-flight-vehicles` "Unsteady Aerodynamics, Aeroelasticity, & Flutter" (https://eaglepubs.erau.edu/introductiontoaerospaceflightvehicles/chapter/unsteady-aerodynamics-flutter/): Unsteady Airfoil Behavior; Frequency Domain Theories; Time-Domain: The Indicial Response; Time- Vs. Frequency-Domain Functions; Dynamic Stall & Stall Flutter (required-free). `@ntrs-nasa-cr-3082-rotary-wing-aerodynamics` Chapter VI, 7. Unsteady Aerodynamic Effects (7.2 Physics of Unsteady Flow; 7.3 Analytical Approach to Oscillations of Thin, Flat Airfoils below Stall) (required-free): the rotor view. `@ntrs-naca-tr-681-jones-unsteady-lift-finite-wing` (supplementary-free): how finite span changes the indicial lift. `@scipy-userguide` "Special functions (scipy.special)" (required-free): `hankel2`. Optional textbooks for deeper study (paid, never required): Drela, *Flight Vehicle Aerodynamics*, for unsteady thin-airfoil theory in modern notation; Leishman, *Principles of Helicopter Aerodynamics*, for indicial methods, dynamic stall and rotor applications; library or used copies are fine. |
| **Extensions** | ★ Add the Sears function (lift on a section in a sinusoidal gust) and compare it with the Fourier transform of Küssner's approximation, as this lesson does for Wagner and Theodorsen. ★ Write Jones's approximation as two first-order ODEs (state-space form, `P02.08`) and integrate them with `solve_ivp`; show that it matches the Duhamel quadrature and costs O(N) instead of O(N log N). ★ Couple `theodorsen_lift` and `theodorsen_moment` to a two-degree-of-freedom typical section (plunge spring, pitch spring, mass offset) and find its flutter speed by sweeping V (the route `P13.18` takes). |
| **Time** | unmeasured |

## 1. The engineering problem

Steady aerodynamics assumes that the flow has forgotten its history: the lift at an angle of attack is 2πα whether the wing reached that angle a minute or a millisecond ago. That is true for a wing that moves slowly, and false in four places an aircraft designer cannot avoid. A gust changes the angle of attack in the time the aircraft flies a few chords. A wing in flutter or buffet vibrates at frequencies where the air cannot keep up. A helicopter blade changes its pitch once per revolution. A control surface deflects in a fraction of a second. In each case the lift is smaller than 2πα, arrives late, and sometimes arrives early.

The physics is circulation. When the lift on a section changes, its bound circulation changes, and Kelvin's theorem (`P10.03`) says the total circulation cannot: an equal and opposite vortex is shed from the trailing edge and carried away by the stream. While it is close it induces a downwash on the section that opposes the change. Only when it has drifted many chords away does the lift settle. That is Wagner's step response, and Theodorsen's function C(k) is its frequency response. Both belong to the same linear system, and part 3 checks that they agree.

The lesson answers the engineering question for the trainer's wing at 5 Hz, then uses the same tools on gusts and on the energy the air feeds into a pitching motion, which is where flutter starts. Success criteria, written first: Theodorsen's function from Hankel functions reaches its limits 1 and ½; the lift tends to quasi-steady theory as k → 0; the C(k) implied by Jones's Wagner approximation is within its known error (about 2 %) of the exact one; Duhamel's quadrature converges at second order; and a time-domain simulation of harmonic pitch, after its transient, reproduces Theodorsen's amplitude and phase to within that same 2 %.

## 2. Understand and derive

### 2.1 Reduced frequency and reduced time

Let the section have chord c = 2b (b is the semi-chord: Theodorsen's notation, used throughout), fly at speed V, and oscillate at angular frequency ω. The time a fluid particle takes to pass one semi-chord is b/V; the period is 2π/ω. Their ratio, made dimensionless,

**k = ωb/V**

is the reduced frequency: the phase angle (in radians) the motion advances while the air moves one semi-chord. For a transient, the matching variable is the reduced time **s = Vt/b**, the distance travelled in semi-chords. Rules of thumb: k < 0.05 is nearly quasi-steady, 0.05–0.2 is unsteady, above 0.2 is highly unsteady. The trainer at 5 Hz: k = 2π·5·0.6/55 = **0.343**.

### 2.2 Theodorsen's lift and moment

Take a flat plate pitching by α(t) (nose up) about a pivot a semi-chords aft of mid-chord (a = −½ is the quarter chord) and plunging by h(t), positive **down** (Theodorsen's sign: a downward velocity ḣ raises the effective angle of attack). Theodorsen (`@ntrs-naca-tr-496-theodorsen-flutter`) split the potential into two parts.

**Non-circulatory (apparent mass).** Accelerating the plate accelerates the air around it. This part acts at once, with no wake, like added mass:

L′_nc = πρb²(ḧ + Vα̇ − baα̈).

**Circulatory.** The bound circulation responds to the downwash at the three-quarter-chord point, the point where thin-airfoil theory (`P10.04`) says a flat plate "sees" its angle of attack:

α₃/₄ = α + ḣ/V + b(½ − a)α̇/V,  L′_c = 2πρVb C(k) V α₃/₄ = ½ρV²(2b) · 2π C(k) α₃/₄.

Every effect of the wake is in the single complex factor C(k). For harmonic motion α = Re(α̂ e^{iωt}), h = Re(ĥ e^{iωt}), d/dt becomes iω and, dividing by q c = ρV²b, the lift coefficient amplitude is

**ĉl = π[(ik)² ĥ/b + ik α̂ − a(ik)² α̂] + 2πC(k)[ik ĥ/b + α̂ + (½ − a) ik α̂].**   (1)

The moment about the pivot (nose-up, cm = M′/(qc²)) has the same structure: apparent-mass terms π[a(ik)²ĥ/b − (½ − a)ikα̂ − (⅛ + a²)(ik)²α̂]/2 plus the circulatory lift times its arm from the quarter chord to the pivot, (a + ½)/2 · 2πC(k)α̂₃/₄. About the quarter chord the circulatory part has no moment at all: the tests check it.

**Theodorsen's function** comes from the shed wake's integral, which Theodorsen evaluated in Bessel functions (his section "Value of the Function C(k)"). In Hankel functions of the second kind (`scipy.special.hankel2`, `@scipy-userguide`):

**C(k) = H₁⁽²⁾(k) / (H₁⁽²⁾(k) + i H₀⁽²⁾(k)) = F(k) + iG(k).**   (2)

| k | 0 | 0.01 | 0.05 | 0.1 | 0.2 | 0.343 | 0.5 | 1 | 2 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| F | 1 | 0.9824 | 0.9090 | 0.8319 | 0.7276 | 0.6458 | 0.5979 | 0.5394 | 0.5130 | 0.5006 |
| G | 0 | −0.0457 | −0.1306 | −0.1723 | −0.1886 | −0.1734 | −0.1507 | −0.1003 | −0.0577 | −0.0124 |
| phase | 0° | −2.7° | −8.2° | −11.7° | −14.5° | −15.0° | −14.1° | −10.5° | −6.4° | −1.4° |

F falls from 1 to ½: the wake removes up to half the circulatory lift. G is negative: the circulatory lift always lags α₃/₄, by at most about 15° (near k = 0.3). C approaches 1 slowly: |1 − C| is still 0.049 at k = 0.01.

**The quasi-steady limit.** As k → 0, C → 1 and the apparent-mass terms (proportional to k and k²) become small next to the circulatory term, so ĉl → 2πα̂₃/₄. That is the independent check of (1): the tests confirm, for pitch, plunge and three pivots, that the relative difference falls with k and is below 2 × 10⁻³ at k = 10⁻⁴.

**The trainer's wing.** Pitch about c/4 (a = −½), k = 0.343: α̂₃/₄ = α̂(1 + ik). Quasi-steady: |ĉl|/α̂ = 2π|1 + ik| = 6.64, leading by 18.9°. Theodorsen: circulatory part 4.441 at +3.9° (the wake has cut 2π|1 + ik| by |C| = 0.669 and rotated it by −15.0°), apparent-mass part 1.092 at +99.7°, total **4.465 = 0.711 × 2π, leading the pitch by 18.0°**, which at 5 Hz is 10.0 ms. The lift oscillation is 29 % smaller than steady theory says, and it leads, because the α̇ terms (three-quarter chord and apparent mass) outweigh the wake's lag at this k.

### 2.3 Wagner's step response and Küssner's gust response

Now change α suddenly from 0 to α₀ at s = 0. The apparent-mass lift is an impulse at the instant of the step; the circulatory lift jumps to half its final value and creeps up as the starting vortex drifts away. The ratio is **Wagner's function φ(s)**, with φ(0) = ½ and φ(∞) = 1. R. T. Jones (`@ntrs-naca-tn-667-jones-operational-nonuniform-lift`) gave a two-exponential approximation:

**φ(s) ≈ 1 − 0.165 e^{−0.0455 s} − 0.335 e^{−0.3 s}.**   (3)

It reaches 75 % at s = 3.70, **90 % at s = 12.70** (7.62 m, 139 ms for the trainer) and 95 % at s = 26.3. A gust is different: a sharp-edged vertical gust reaches the leading edge first and the trailing edge one chord later, so the lift starts from zero. **Küssner's function ψ(s)**, in Sears's approximation ψ ≈ 1 − 0.5e^{−0.13s} − 0.5e^{−s}, reaches 50 % at s = 1.65 and 90 % at s = 12.4.

### 2.4 Duhamel's integral and the link to C(k)

The system is linear, so any input is a sum of small steps. A step dα₃/₄ at σ contributes 2π φ(s − σ) dα₃/₄ at s. Adding them (Duhamel's integral, a convolution, `P02.08`):

**cl_c(s) = 2π [α₃/₄(0) φ(s) + ∫₀ˢ φ(s − σ) (dα₃/₄/dσ) dσ].**   (4)

Feed it α₃/₄ = e^{iks}, started long ago. With τ = s − σ,

cl_c = 2π e^{iks} · ik ∫₀^∞ φ(τ) e^{−ikτ} dτ = 2π α₃/₄ · ik Φ(ik),

where Φ(p) is the Laplace transform of φ. Comparing with (1): **C(k) = ik Φ(ik)**. Theodorsen's function is the frequency response of Wagner's system. For the exponential form (3), Φ(p) = 1/p − Σ Aᵢ/(p + βᵢ), so

**C_J(k) = 1 − 0.165 ik/(ik + 0.0455) − 0.335 ik/(ik + 0.3).**   (5)

Its limits are right: C_J(0) = 1 and C_J(∞) = 1 − 0.165 − 0.335 = ½. How close is it in between? Comparing (5) with (2) on 20 001 points of 0 < k ≤ 5, **the largest difference is 0.0145, at k = 0.41, which is 2.3 % of |C|**. That is the accuracy of Jones's fit, and it is the error to expect when a time-domain simulation built on (3) is compared with Theodorsen. Integrating by parts gives the form used for a direct numerical check: C(k) = φ(0) + ∫₀^∞ φ′(τ)e^{−ikτ}dτ, whose integrand decays.

### 2.5 Lift loops, work and the start of flutter

For α = ᾱ + α̂ sin θ, cl(α) traces an ellipse: a **lift loop**. Its area ∮cl dα = π Im(α̂* ĉl) is zero in steady flow and grows with k; its direction says whether lift leads or lags. For the trainer's c/4 pitch at α = 4° ± 2° the loop area is −0.00145 at k = 0.05 (lift lags) and +0.00527 at k = 0.343 (lift leads).

Lift loops are not energy. The work the air does on a motion is the generalized force times the displacement, ∮M dα for pitch and ∮L d(−h) for plunge (h is down), and for harmonic motion **W = π Im(x̂* f̂)**. Positive W means the air feeds energy into the motion: negative aerodynamic damping. Two results follow from (1) and the moment:

- **Pure plunge is always damped.** The work is −0.897, −3.94, −10.6 (ĥ/b = 1) at k = 0.05, 0.3, 1: negative, and the tests check 40 values of k.
- **Pure pitch is damped for pivots at or behind the quarter chord**, at every k from 10⁻⁴ to 10. For pivots ahead of it, energy is fed in below a small k: k < 0.0015 for a pivot at 0.20c, 0.0141 at 0.15c, 0.0268 at 0.10c, 0.0403 at the leading edge. The lagging circulatory lift acts at the quarter chord, behind such a pivot, and at low k its lag puts it in phase with the pitch velocity.

Neither single-degree-of-freedom result is classical flutter. That needs bending and torsion coupled by the structure's mass and stiffness, so that the lift from one does positive work on the other; the structural model comes in `P13.18`, and these functions are its aerodynamic input.

### 2.6 Where the model stops

- **Separation.** The theory needs attached flow. Near stall, oscillation produces dynamic stall: a leading-edge vortex, lift overshoot, a moment collapse and loops that can reverse and feed energy into pitch (stall flutter). Leishman's "Dynamic Stall & Stall Flutter" describes it.
- **Finite span.** The trailing vortices add their own lag; for low aspect ratio the lift builds more slowly and to a lower final value (`@ntrs-naca-tr-681-jones-unsteady-lift-finite-wing`).
- **Compressibility.** Pressure signals travel at the speed of sound; at high subsonic Mach numbers the indicial response changes shape.
- **Rotors.** The shed wake passes under the following blades; the flat-wake assumption fails.

## 3. Build the computational model

Open `starter/main.py`. `fit_harmonic`, `loop_area`, `observed_order`, `reduced_frequency` and `main` are given.

### Step 1: the frequency domain

`theodorsen(k)` is (2), returning exactly 1 at k = 0 (the Hankel functions are singular there). `theodorsen_from_indicial(k, A, beta)` is (5) for any number of exponentials. `theodorsen_lift(k, α̂, ĥ/b, a, C=None)` returns the circulatory and non-circulatory parts of (1) separately, so you can see each; passing `C` lets you swap in Jones's C_J. `theodorsen_moment` is the moment of §2.2, `quasi_steady_lift` is 2πα̂₃/₄, `aero_work_per_cycle(x̂, f̂)` is π Im(x̂* f̂).

### Step 2: the time domain

`wagner_jones` and `kussner_sears` are (3) and Sears's form. `duhamel(s, x, φ)` evaluates (4) on a uniform grid by the midpoint rule on each interval: y_n = x₀φ(s_n) + Σ_{j≤n} (x_j − x_{j−1}) φ(s_n − s_{j−½}). The sum is a discrete convolution, done with `scipy.signal.fftconvolve`:

```python
y = x[0] * indicial(s - s[0])
kernel = indicial((np.arange(n - 1) + 0.5) * ds)
y[1:] += fftconvolve(np.diff(x), kernel)[: n - 1]
```

`time_domain_lift(s, α, h/b, a)` forms α₃/₄ with second-order finite differences, applies `duhamel` with Wagner for the circulatory part, and adds the apparent mass π(h″/b + α′ − aα″). `simulate_harmonic_pitch` starts α = α̂ sin ks from rest, runs 12 cycles and fits the last two with `fit_harmonic`. `rise_distance` uses `brentq`. `gust_response` is Duhamel with Küssner on w(s)/V. `lift_loop` samples one cycle of (1). `pitch_flutter_boundary` finds where the pitch work changes sign.

### Verify

Running `main()` prints the table of §2.2 and

```
max |C_Jones - C| over 0 < k <= 5: 0.0145 at k = 0.410 (2.34 % of |C| at most)
trainer section, c = 1.2 m at 55.0 m/s, pitch about c/4 at 5.0 Hz: k = 0.3427
  total Theodorsen: |4.4647| = 0.711 x steady, phase +17.98 deg (time +9.99 ms)
Wagner (Jones): 90 % of steady lift at s = 12.70 semi-chords = 7.62 m = 138.5 ms at 55.0 m/s
Duhamel midpoint rule vs exact convolution (alpha = sin 0.3 s, 0 <= s <= 40): errors 3.84e-04, 9.61e-05, 2.40e-05, 6.00e-06; observed orders [2. 2. 2.]
```

The 17 tests in `tests/test_main.py` check the limits, linearity and checks listed in the metadata. The convergence check uses a closed form: for x = sin ks and Jones's φ, the convolution integrates term by term (each exponential gives k[β cos ks + k sin ks − β e^{−βs}]/(β² + k²)). The midpoint quadrature's error falls by 4.0 each time the step halves: **observed order 2.00**.

### Compare: Wagner against Theodorsen in the time domain

The decisive check runs the same physics two ways: a time-domain Duhamel simulation of harmonic pitch about c/4, started from rest, against Theodorsen's periodic solution.

| k | |cl|/α̂, time domain | with C_J | with exact C | phase TD | phase exact | TD vs C_J | TD vs exact |
|---|---|---|---|---|---|---|---|
| 0.05 | 5.7134 | 5.7134 | 5.7610 | −4.19° | −3.76° | 2.3 × 10⁻⁵ | 1.1 % |
| 0.1 | 5.3036 | 5.3036 | 5.3254 | −2.02° | −2.65° | 1.0 × 10⁻⁵ | 1.2 % |
| 0.2 | 4.8397 | 4.8397 | 4.7592 | +4.30° | +4.31° | 7.1 × 10⁻⁶ | 1.7 % |
| 0.343 | 4.4829 | 4.4829 | 4.4647 | +16.78° | +17.98° | 1.5 × 10⁻⁵ | 2.1 % |
| 0.5 | 4.5173 | 4.5174 | 4.5815 | +32.13° | +33.11° | 2.3 × 10⁻⁵ | 2.2 % |
| 1.0 | 6.2970 | 6.2970 | 6.3888 | +67.86° | +67.46° | 5.1 × 10⁻⁵ | 1.6 % |

Against Theodorsen built on the same Jones C_J the agreement is 10⁻⁵: the time stepping is right. Against the exact C(k) it is 1–2 % and up to 1.2° of phase, exactly the size of Jones's fit error from §2.4. Every difference has a named cause. (Lift amplitude rises again at k = 1 because the apparent-mass term grows as k².)

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction before running it.

1. **How good is Jones?** Compare C_J with the exact C(k) on 0 < k ≤ 5. Predict where the error peaks. *Finding:* 0.0096–0.0145 across 0.01 ≤ k ≤ 1, largest at k = 0.41; 0.0017 at k = 10.
2. **The trainer at 5 Hz.** Predict whether the unsteady lift is larger or smaller than 2πα̂ and whether it leads or lags. *Finding:* 0.711 × 2πα̂, leading by 18.0° (10.0 ms). Smaller and leading: the wake cuts the circulatory lift and delays it 15.0°, but the α̇ terms lead.
3. **Step and gust.** How far does the trainer fly before lift reaches 90 % after a step in α, and after entering a sharp-edged gust? *Finding:* 12.7 and 12.4 semi-chords (7.6 m and 7.4 m), but the gust starts from zero: 50 % takes 0 semi-chords after a step (φ(0) = ½) and 1.65 semi-chords in a gust. For a sharp-edged 5 m/s gust at 55 m/s (quasi-steady Δcl = 0.571), the section has 0.312 after one chord and 0.493 after five. For a "1 − cos" gust of gradient length H, the peak is 0.467, 0.642, 0.801, 0.903 and 0.966 of quasi-steady for H = 1, 2.5, 6.25, 12.5 and 25 chords: short gusts are softened most.
4. **Lift loops.** α = 4° ± 2° about c/4 at k = 0.05, 0.1, 0.343, 1. Predict when the loop direction reverses. *Finding:* areas −0.00145, −0.00094, +0.00527, +0.02259: lagging lift at low k, leading above about k = 0.2 (the phase table of §3 changes sign between 0.1 and 0.2).
5. **Pitch damping and pivot.** Work per cycle for pivots at 0.25c to 0.50c and ahead of 0.25c. *Finding:* the table of §2.5; damped for every pivot at or behind 0.25c, energy fed in at very low k for pivots ahead of 0.2c.

Raw data format: CSV tables `k,F,G,abs_C,phase_deg,dC_jones` (1); `k,abs_td,abs_jones,abs_exact,phase_td,phase_exact` (2 and §3); `level,s_wagner,s_kussner` (3); `H_chords,peak_ratio` (3); `k,loop_area` (4); `a,k,work` (5); plus `unsteady_aero.png`.

## 5. Validate and improve

| quantity | analysis | computation | difference | cause |
|---|---|---|---|---|
| C(k) as k → 0 | 1 | `theodorsen`: |1 − C| = 0.0489, 0.0072, 0.0009 at k = 10⁻², 10⁻³, 10⁻⁴ | → 0 | slow approach, as the theory says |
| lift as k → 0 | 2πα₃/₄ | relative difference < 2 × 10⁻³ at k = 10⁻⁴ (tests) | → 0 | apparent mass ~ k, wake ~ 1 − C |
| C(k) via Fourier transform of φ | exact (2) | C_J (5) | ≤ 0.0145 | Jones's two-term fit |
| C_J closed form vs numerical transform | (5) | trapezoid ∫φ′e^{−ikτ} | < 2 × 10⁻⁵ | quadrature |
| harmonic pitch, time vs frequency domain | Theodorsen with C_J | Duhamel simulation | 10⁻⁵ | discretisation, order 2.0 |
| same, exact C | Theodorsen | Duhamel simulation | 1.1–2.2 %, ≤ 1.2° | Jones's fit |

Every difference is explained. The engineering answer: **at 5 Hz the trainer's wing section produces 71 % of the quasi-steady lift amplitude per unit pitch, leading the pitch by 18° (10 ms), because the wake reduces and delays the circulatory lift (|C| = 0.669, −15.0°) while the pitch-rate terms lead; after a sudden change of angle of attack, 90 % of the lift arrives after 7.6 m (139 ms).**

What to improve: the wing is finite (AR 8.33), so the trailing vortices reduce the steady lift slope to about 4.7 per radian (`P10.20`) and add lag; the next step is a strip-theory or finite-wing indicial function (TR-681). What not to improve for this question: Jones's 2 % fit error, which is smaller than the finite-span effect.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of C(k) = ikΦ(ik) and of (5); `unsteady_aero.png`; the five experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) giving the amplitude ratio, phase and rise distance for the trainer, with the assumption that matters most.

### Adversarial self-checks

- **What if the frequency goes to zero?** C → 1, the apparent-mass terms vanish, and the lift is 2πα₃/₄ → 2πα: the tests confirm the limit for pitch, plunge and three pivots.
- **What if the frequency goes to infinity?** C → ½ and the apparent-mass terms (∝ k²) dominate: lift becomes an inertial reaction, in phase with acceleration, not with α.
- **Which assumption fails first for the trainer?** Finite span: a section analysis of an AR 8.33 wing overstates the lift slope by about a third, and the trailing vortices add lag. Then separation, if the vibration happens near stall.
- **What measurement could falsify the model?** Unsteady surface pressures on an oscillating wind-tunnel model below stall: the measured lift amplitude and phase against k should follow (1) with a finite-span correction; a phase that does not reverse sign between k = 0.1 and 0.2 for c/4 pitch would contradict it.
- **Does the result survive a change of method or step?** Frequency domain, time-domain Duhamel and a direct numerical transform agree to 10⁻⁵ with the same indicial function; halving the time step cuts the error by 4.
- **Which decision changes if the uncertainty doubles?** If Jones's error were 4 % instead of 2 %, it would still be smaller than the finite-span effect, so the choice of what to improve next would not change; a flutter margin computed with these coefficients would need a few per cent more speed margin.

## Check yourself

1. Define the reduced frequency. What is k for the trainer's wing (c = 1.2 m) at 55 m/s and 5 Hz, and is that quasi-steady?
2. What are the limits of C(k) as k → 0 and k → ∞, and what does each mean physically?
3. What is the difference between the non-circulatory and the circulatory lift?
4. What fraction of the final lift does a section have immediately after a step in angle of attack? After how many semi-chords does it reach 90 % (Jones)?
5. Why can Duhamel's integral give the lift for any motion, and when does it stop being valid?
6. What does the area of a lift loop measure, and why is it not the energy the air feeds into the motion?
7. Is pure plunge ever unstable in this theory? For which pivot positions can pure pitch be fed energy?
