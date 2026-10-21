# P02.07  Complex Numbers & Phasors

> A sinusoid is the shadow of a rotating arrow. Add the arrows, multiply the arrows, and the trigonometry does itself.

| | |
|---|---|
| **Learning outcomes** | 1. Derive Euler's formula from the power series, de Moivre's theorem and the n-th roots of a complex number, and show that multiplication by e^{iθ} is the rotation of P01.15. 2. Implement and test a phasor calculator: polar form, the exponential from its series, roots, the sum of same-frequency sinusoids, differentiation as multiplication by iω, the complex gain of a mass–spring–damper and the transmissibility of a mount, each checked against a time-domain integration with `solve_ivp`. 3. Measure the amplitude and phase of a simulated vibration record by least squares and compare them with the phasor prediction. 4. Decide what phase a propeller synchrophaser should hold and how soft an instrument mount must be to isolate a tray from 40 Hz cabin vibration. |
| **Prerequisites** | `P02.02` (damped oscillator), `P02.03` (forced response, gain and phase lag), `P01.05` (Taylor series of e^x, cos, sin), `P01.15` (rotation matrices), `P01.17` (least squares, normal equations), `P00.03`, `P00.07`. Skill check: you can add two vectors in the plane, use the cosine rule, and write cos(ωt + φ) in terms of cos ωt and sin ωt. |
| **Engineering question** | Two propellers at 2400 rpm shake the cabin floor of a twin at 40 Hz with 1.2 and 0.9 m/s²; what floor vibration results for each phase difference between them, what does an instrument tray on a 30 Hz mount see, and how soft must the mount be to pass only 20 % of it? |
| **Assumptions and limits** | Every signal is a single sinusoid at exactly 40 Hz: the phasor method is for one frequency at a time and for linear systems in steady state. The floor responds linearly; the tray is a single mass on a linear spring and viscous damper; the transient has died away (it decays as e^{−ζω_n t}, about 0.1 s here). The method says nothing about the transient, nothing about two sources at *different* frequencies (they beat, problem 3), and nothing about a nonlinear mount. |
| **Analytical work** | Euler's formula from the series; the exponential law e^{a+b} = e^a e^b and de Moivre; n-th roots and roots of unity; sinusoid as Re(Z e^{iωt}); phasor addition and the cosine rule; differentiation as multiplication by iω; the complex gain H = 1/(1 − r² + 2iζr) by substitution and why taking the real part is allowed; the peak of |H|; transmissibility and the crossover at r = √2. |
| **Computational deliverable** | `starter/main.py`: `to_polar`, `from_polar`, `exp_series`, `euler`, `de_moivre`, `nth_roots`, `rotate_2d`, `phasor`, `sinusoid`, `from_phasor`, `phasor_derivative`, `add_sinusoids`, `combined_amplitude`, `amplitude_range`, `fit_sinusoid`, `complex_gain`, `gain_and_phase_lag`, `steady_state_displacement`, `forced_oscillator_rhs`, `transmissibility`, `mount_frequency_for_transmissibility`. Tests prove: the series exponential matches `numpy.exp` to 10⁻¹³ and gives Euler's formula with unit modulus; e^{a+b} = e^a e^b; every n-th root raised to n returns z and the roots are equally spaced and sum to zero; complex multiplication equals the rotation matrix; phasor sums equal pointwise sums of the signals; the two-source amplitude obeys the cosine rule with the right extremes; the least-squares fit recovers amplitude and phase exactly from a clean record and within 10 σ from a noisy one; the complex gain has the static, resonant and mass-controlled limits and a monotone phase lag; the phasor steady state agrees with `solve_ivp` to 2 × 10⁻⁴ of the amplitude and satisfies the ODE; transmissibility equals 1 at r = √2 for every damping. |
| **Practical activity** | **computational experiment**: a sweep of combined floor vibration against propeller phase, of tray response against mount frequency, and a measurement of how fast the time-domain transient converges to the phasor prediction. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `phasor_calculator.png`, the hand derivation of Euler's formula and of the complex gain, the two sweep tables, and a short note recommending a synchrophaser setting and a mount frequency with their numbers. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: Euler's formula and the complex gain are derived with every step; all public tests pass; the hand phasor sum for one phase difference agrees with the code to three figures; the mount recommendation states the transmissibility it achieves and the frequency at which the new mount would resonate. |
| **References** | `@mit-ocw-18-03-differential-equations` L6 "Complex Numbers and Complex Exponentials", notes 5 "Complex numbers, roots of unity", 6 "Complex exponentials; sinusoidal functions", 7 "Linear system response to exponential and sinusoidal input; gain, phase lag" (required-free); `@openstax-precalculus-2e` §8.5 "Polar Form of Complex Numbers" (required-free, refresher); `@openstax-calculus-volume-2` §6.4 "Working with Taylor Series" (supplementary-free, the series behind Euler's formula). |
| **Extensions** | ★ Extend `fit_sinusoid` to fit several frequencies at once (a matrix with 2n columns) and separate the 40 Hz and 80 Hz content of a synthetic record. ★ Derive the impedance of a series RLC circuit and show it is the same complex number as k − mω² + icω with m ↔ L, c ↔ R, 1/k ↔ C. ★ Add a second mass (tray on a mount on a flexible floor) and show the transmissibility now has two peaks. |
| **Time** | unmeasured |

## 1. The engineering problem

In a twin-propeller aircraft each propeller shakes the structure once per revolution: at 2400 rpm that is 40 Hz. At the instrument tray the two excitations arrive with amplitudes a₁ = 1.2 m/s² and a₂ = 0.9 m/s², at the same frequency, with a phase difference that depends on where the two propellers are in their rotation relative to each other. A **synchrophaser** holds that relative angle fixed. The question for the cabin designer is simple to state and awkward to answer with trigonometry: what is the amplitude of 1.2 cos ωt + 0.9 cos(ωt + δ) for every δ, and which δ should the synchrophaser hold?

The second question concerns the instrument tray itself. It sits on a rubber mount, a spring and damper with natural frequency 30 Hz and damping ratio 0.08, and the mount was chosen for shock, not for 40 Hz vibration. Does it isolate the tray or make things worse? P02.03 answered such questions by undetermined coefficients, two unknowns and a page of algebra per case. This lesson replaces that page with one complex division. The trick is the **phasor**: a sinusoid A cos(ωt + φ) is the real part of the complex number A e^{iφ} e^{iωt}, a fixed arrow A e^{iφ} rotating at rate ω. Addition of same-frequency sinusoids is addition of arrows; differentiation is multiplication by iω; and the steady-state response of any linear system is the input arrow times a complex number, the **complex gain**.

Success criteria, written down first: the exponential computed from its series agrees with NumPy's to 10⁻¹³ and yields Euler's formula to 10⁻¹⁴; the phasor sum of the two propeller signals agrees with pointwise addition to 10⁻¹²; the steady-state tray response from the complex gain agrees with `solve_ivp` to 2 × 10⁻⁴ of the amplitude once the transient has died; and the mount frequency that gives a transmissibility of 0.2 at 40 Hz is found to 10⁻⁹.

## 2. Understand and derive

### 2.1 Complex numbers in two coordinates

A complex number z = x + iy, with i² = −1, is a point in the plane (`@openstax-precalculus-2e` §8.5). Its **polar form** is z = r(cos θ + i sin θ) with modulus r = |z| = √(x² + y²) and argument θ = atan2(y, x), which we keep in (−π, π]. Addition is vector addition. Multiplication is the surprise: (r₁, θ₁)(r₂, θ₂) = (r₁r₂, θ₁ + θ₂), moduli multiply and arguments add. Check with (cos θ₁ + i sin θ₁)(cos θ₂ + i sin θ₂) = (cos θ₁ cos θ₂ − sin θ₁ sin θ₂) + i(sin θ₁ cos θ₂ + cos θ₁ sin θ₂) = cos(θ₁ + θ₂) + i sin(θ₁ + θ₂): the addition formulae of trigonometry *are* complex multiplication.

In particular, multiplying by cos θ + i sin θ rotates a point anticlockwise by θ without changing its length. Written in components, (x + iy)(cos θ + i sin θ) = (x cos θ − y sin θ) + i(x sin θ + y cos θ), which is the rotation matrix of P01.15 applied to (x, y). `rotate_2d` does this with one complex multiplication and the test compares it with the matrix.

### 2.2 Euler's formula from the series

Define e^z for complex z by the same power series as for real z (P01.05; `@openstax-calculus-volume-2` §6.4), which converges for every z:

e^z = Σ_{k≥0} z^k/k! = 1 + z + z²/2 + z³/6 + …

Put z = iθ and use i² = −1, i³ = −i, i⁴ = 1:

e^{iθ} = (1 − θ²/2! + θ⁴/4! − …) + i(θ − θ³/3! + θ⁵/5! − …) = **cos θ + i sin θ**.   (1)

The real part is the series of cos θ and the imaginary part the series of sin θ, so **Euler's formula** is not an assumption but an identity between three power series. Two consequences: |e^{iθ}| = cos²θ + sin²θ = 1, every pure-imaginary exponent gives a point on the unit circle; and the polar form becomes z = r e^{iθ}. Because the series obeys e^{a+b} = e^a e^b (multiply the two series and collect by the binomial theorem; the test checks it numerically), the multiplication rule of 2.1 is now one line: r₁e^{iθ₁} r₂e^{iθ₂} = r₁r₂ e^{i(θ₁+θ₂)}. The code computes e^z from the series, term by term, and observes (1) rather than assuming it: `exp_series(1j * np.pi)` returns −1 + 3 × 10⁻¹⁶ i.

### 2.3 Powers and roots

**De Moivre:** (cos θ + i sin θ)ⁿ = (e^{iθ})ⁿ = e^{inθ} = cos nθ + i sin nθ. So z⁵ for z = 4.5 + 4i (r = 6.0208, θ = 41.63°) is 6.0208⁵ e^{i × 208.2°} = 7911.7 at −151.8°, i.e. −6975 − 3735i, with no binomial expansion.

**Roots:** w = z^{1/n} means wⁿ = z. Writing z = r e^{i(θ + 2πk)} for any integer k (the same point, k extra turns) and taking the n-th root, w_k = r^{1/n} e^{i(θ + 2πk)/n}, k = 0, …, n − 1: **n distinct roots**, equally spaced by 2π/n on the circle of radius r^{1/n}. For z = 1 they are the **roots of unity** ω_k = e^{2πik/n}, and their sum is zero for n > 1 (a geometric series with ratio e^{2πi/n} ≠ 1 summing to (1 − e^{2πi})/(1 − e^{2πi/n}) = 0). The cube roots of 8 are 2, −1 ± 1.732i. The test raises each root to the n and checks the spacing and the zero sum. Roots of unity are the phases of the blades of an N-bladed propeller, which is why a balanced propeller produces no once-per-revolution force (problem 6), and they are the heart of the DFT in P02.10.

### 2.4 Phasors: a sinusoid as a rotating arrow

Write a sinusoid as the real part of a complex exponential:

**A cos(ωt + φ) = Re(A e^{iφ} e^{iωt}) = Re(Z e^{iωt})**,   Z = A e^{iφ}.   (2)

Z is the **phasor** of the signal: the arrow at t = 0. As time runs the arrow Z e^{iωt} rotates at ω rad/s and the signal is its projection on the real axis (`@mit-ocw-18-03-differential-equations` notes 6).

![The complex plane with real axis Re and imaginary axis Im. The phasor Z = A e^{iφ} is the arrow at t = 0; as time runs the arrow Z e^{iωt} rotates counter-clockwise at ω around a circle of radius A, at angle φ + ωt from the real axis, and its projection onto the real axis is the signal x(t) = Re(Z e^{iωt}).](figures/rotating-phasor.svg)
*Figure 1. A sinusoid as a rotating arrow. The phasor Z fixes the amplitude and the phase at t = 0; the signal is the shadow of the rotating arrow on the real axis.*

Three operations become trivial.

**Addition.** For sinusoids at the same ω, Re(Z₁e^{iωt}) + Re(Z₂e^{iωt}) = Re((Z₁ + Z₂)e^{iωt}): add the arrows. For the two propellers, Z₁ = 1.2, Z₂ = 0.9 e^{iδ}, and the cosine rule gives

**|Z₁ + Z₂| = √(a₁² + a₂² + 2a₁a₂ cos δ)**,   (3)

from 2.1 m/s² in phase (δ = 0) to 0.3 m/s² in anti-phase (δ = 180°), with 1.50 m/s² at 90° (quadrature: Pythagoras). The synchrophaser is worth a factor 7 at this tray. The code's table:

| δ | combined amplitude [m/s²] | phase of the sum |
|---|---|---|
| 0° | 2.1000 | 0.00° |
| 45° | 1.9435 | 19.11° |
| 90° | 1.5000 | 36.87° |
| 135° | 0.8501 | 48.47° |
| 180° | 0.3000 | 0.00° |

**Differentiation.** d/dt Re(Z e^{iωt}) = Re(iωZ e^{iωt}): the derivative's phasor is iωZ, the arrow scaled by ω and turned by 90°. (Check: d/dt A cos(ωt + φ) = −Aω sin(ωt + φ) = Aω cos(ωt + φ + π/2). ✓) Integration is division by iω. This is why every constant-coefficient linear ODE becomes algebra for sinusoidal inputs.

**Taking the real part is allowed** because the operations are linear with real coefficients: for a real constant a, a Re(w) = Re(aw), and Re(w₁) + Re(w₂) = Re(w₁ + w₂). It is *not* allowed inside products: Re(w₁)Re(w₂) ≠ Re(w₁w₂), which is why the phasor method stops at linear systems and why power and energy calculations need the conjugate.

### 2.5 The complex gain

Apply this to the tray: m ẍ + c ẋ + k x = F₀ cos ωt, F₀ real. Seek the steady state as x = Re(X e^{iωt}). Substituting and using the differentiation rule, every term carries the factor e^{iωt}:

Re[(−mω² + icω + k) X e^{iωt}] = Re[F₀ e^{iωt}] for all t,   so   **X = F₀/(k − mω² + icω)**.   (4)

(Equality of the real parts for all t forces equality of the complex coefficients: evaluate at ωt = 0 and ωt = −π/2.) One division replaces the undetermined-coefficients algebra of P02.03. In dimensionless form, with ω_n = √(k/m), ζ = c/(2√(km)), r = ω/ω_n and the static deflection F₀/k,

**H(r) = X k/F₀ = 1/(1 − r² + 2iζr)**,   |H| = 1/√((1 − r²)² + (2ζr)²),   phase lag = atan2(2ζr, 1 − r²).   (5)

|H| is the dynamic magnification; its limits are 1 (static, r → 0), 1/(2ζ) at r = 1 (resonance, 6.25 for ζ = 0.08), and 1/r² (mass-controlled, r → ∞). The lag rises monotonically from 0 through 90° at resonance to 180°. The peak of |H| is not exactly at r = 1 but at r = √(1 − 2ζ²) = 0.9936 for ζ = 0.08 (set the derivative of (1 − r²)² + 4ζ²r² to zero); the test checks it to 2 × 10⁻⁵. For the tray at r = 40/30 = 1.3333: H = −1.1958 − 0.3280i, |H| = 1.2399, lag 164.66°. A 5 N force at 40 Hz gives 1.2399 × 5/88 826 m = 69.79 µm, 164.66° behind the force. `solve_ivp` from rest, fitted after 0.6 s: 69.795 µm at −164.66°, the difference 0.008 µm being the remaining transient (e^{−ζω_n t} = e^{−15.1 × 0.6} = 1.2 × 10⁻⁴).

### 2.6 Transmissibility of a mount

The tray is not driven by a force; its base moves. With base motion y = Y cos ωt the spring and damper act on the relative displacement: m ẍ + c(ẋ − ẏ) + k(x − y) = 0, i.e. m ẍ + c ẋ + k x = c ẏ + k y. Phasors: (k − mω² + icω) X = (k + icω) Y, so

**T(r) = |X/Y| = |1 + 2iζr| / |1 − r² + 2iζr|**.   (6)

T is the **transmissibility**: the fraction of the floor motion (equally, of the floor acceleration, since both are multiplied by ω²) that reaches the tray. Two facts decide every isolation design. T = 1 at r = 0 and *again* at r = √2 whatever ζ (set |1 + 2iζr|² = |1 − r² + 2iζr|²: 1 = (1 − r²)², r² = 2); between them the mount amplifies; above √2 it isolates, and there *more* damping transmits *more* (the damper is a stiff path at high frequency). For the tray at r = 1.333: T = 1.2678. The 30 Hz mount makes the 40 Hz vibration 27 % worse: 2.10 m/s² on the floor becomes 2.66 m/s² on the tray. To pass 20 % the mount must be softer: bisection on (6) gives f_n = 15.84 Hz (r = 2.525); 0.5 needs 22.87 Hz, 0.1 needs 11.33 Hz. The price of a soft mount is a large static sag (g/ω_n² = 1.0 mm at 15.8 Hz) and a resonance at 15.8 Hz that something else, an engine idling at 950 rpm say, may excite. That trade is the recommendation the note must make.

### 2.7 Where the model stops being valid

- **One frequency.** Two propellers at 2400 and 2415 rpm do not add as phasors: the sum is a 40 Hz signal whose amplitude sweeps between 0.3 and 2.1 every 4 s, a **beat** (problem 3). The synchrophaser's job is precisely to make the frequencies equal so that phasor addition applies.
- **Steady state.** The complex gain says nothing about the first few cycles after a change; the transient decays as e^{−ζω_n t} and is the subject of P02.08.
- **Linearity.** A rubber mount stiffens with amplitude and a floor has its own modes near 40 Hz; both make H depend on the amplitude or add peaks to (6).
- **Real coefficients.** The step "equality of real parts implies equality of coefficients" needs the coefficients in (4) to be constants; for a time-varying m or k the method fails.

## 3. Build the computational model

Open `starter/main.py`. `cabin_vibration_case()` holds the numbers; every other function takes its arguments explicitly.

### Step 1: polar form and the series exponential

```python
def to_polar(z):
    z = complex(z)
    return abs(z), float(np.angle(z))             # angle in (-pi, pi]

def exp_series(z, n_terms=40):
    term = total = complex(1.0, 0.0)
    for k in range(1, n_terms):
        term *= z / k                             # z^k / k! from the previous term
        total += term
    return total
```

Build each term from the previous one; computing z^k and k! separately overflows the factorial by k = 171 and loses accuracy long before. Forty terms give 10⁻¹⁵ for |z| ≤ 5; the test also checks that four terms are visibly *not* enough at z = 3, so that you see truncation error rather than trust the function. `euler(theta)` returns the closed form; `de_moivre` multiplies the unit phasor n times; `nth_roots` implements 2.3 with `np.exp(1j * (th + 2 * np.pi * k) / n)`.

### Step 2: phasors

```python
def phasor(amplitude, phase):              return from_polar(amplitude, phase)       # (2)
def from_phasor(Z, omega, t):              return np.real(Z * np.exp(1j * omega * t))
def phasor_derivative(Z, omega):           return 1j * omega * Z
def add_sinusoids(components):             return to_polar(sum(phasor(a, p) for a, p in components))
def combined_amplitude(a1, a2, delta):     return np.sqrt(a1**2 + a2**2 + 2 * a1 * a2 * np.cos(delta))   # (3)
```

`phasor` refuses a negative amplitude: a negative amplitude is a phase shift of π, and letting both representations exist makes comparisons ambiguous.

### Step 3: reading a phasor off a record

```python
def fit_sinusoid(t, x, omega):
    M = np.column_stack([np.cos(omega * t), np.sin(omega * t)])
    a, b = np.linalg.solve(M.T @ M, M.T @ x)          # normal equations (P01.17)
    return to_polar(complex(a, -b))                   # A cos(wt + phi) = a cos wt + b sin wt
```

Since A cos(ωt + φ) = A cos φ cos ωt − A sin φ sin ωt, the fitted a = A cos φ and b = −A sin φ, so Z = a − ib. This is the instrument that turns a time series, simulated here and measured in Phase 05, into an amplitude and a phase. With 5001 samples and noise of rms 0.1 the fit error in a and b is about 0.1√(2/5001) = 0.002; the test allows ten times that.

### Step 4: gain and transmissibility

```python
def complex_gain(r, zeta):          return 1.0 / (1.0 - r**2 + 2j * zeta * r)                      # (5)
def transmissibility(r, zeta):      return np.abs(1 + 2j * zeta * r) / np.abs(1 - r**2 + 2j * zeta * r)   # (6)
```

`gain_and_phase_lag` returns |H| and the lag in [0, π]; `steady_state_displacement` is (4) turned into a time signal; `forced_oscillator_rhs` is the same ODE for `solve_ivp`; `mount_frequency_for_transmissibility` bisects (6) on r ∈ [√2, 50], where T is monotone, and returns f_drive/r.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-07-complex-numbers-phasors`. The tests:

1. **Analytic cases.** `exp_series` equals `np.exp` to 10⁻¹³ at six points; e^{iθ} has real part cos θ, imaginary part sin θ and unit modulus at thirteen angles; e^{iπ} = −1; e^{a+b} = e^a e^b; de Moivre for n = 0, 5, 7; every n-th root raised to n returns z, the roots are equally spaced and sum to zero; `rotate_2d` equals the rotation matrix and preserves length.
2. **Phasor algebra.** `from_phasor` reproduces `sinusoid`; the derivative phasor has modulus Aω and leads by 90° and matches `np.gradient` to 10⁻⁴; the phasor sum of four random components equals the pointwise sum to 10⁻¹²; the two-source amplitude obeys (3), is even in δ, and lies within [|a₁ − a₂|, a₁ + a₂].
3. **The instrument.** `fit_sinusoid` recovers A and φ to 10⁻¹² from a clean record and to 0.02 and 0.03 rad from a seeded noisy one.
4. **Limiting cases.** H(0) = 1, |H(50)| < 10⁻³, |H(1)| = 1/(2ζ) with 90° lag, lag monotone from 0 to π, peak at √(1 − 2ζ²); T(0) = T(√2) = 1 for ζ = 0, 0.05, 0.3, 1; T(4, 0) = 1/15; more damping transmits more above √2 and less at resonance.
5. **Agreement with the time domain.** The phasor steady state and `solve_ivp` agree to 2 × 10⁻⁴ of the amplitude after 0.6 s; the fitted amplitude and phase agree with X to 10⁻³; the phasor solution satisfies the ODE under second-order finite differences.
6. **Design.** The mount frequency for T = 0.5, 0.25, 0.1 gives back the target transmissibility to 10⁻⁹, lies below f_drive/√2, and decreases with the target.

### Compare

`solve_ivp` integrates the same ODE from rest and the comparison is in `main()`: 69.794 µm at −164.66° from the phasor, 69.795 µm at −164.66° from the fitted integration, with a maximum pointwise difference of 0.008 µm over the last 0.4 s. The integrator needs about 10⁴ steps and a fit to deliver two numbers; the phasor needs one division. But the integrator delivers the *transient* as well, which the phasor cannot: at t = 0.05 s the tray is still 20 % away from its steady state. In practice you use the complex gain for every steady-state and frequency-response question (P02.08 generalises it to the transfer function) and the integrator whenever the start-up, a shock, or a nonlinear element matters. `scipy.signal.freqs` evaluates (5) for an arbitrary rational transfer function; we wrote the one-line version to see that nothing is hidden in it.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Three sweeps with `main.py` as the instrument; all records are synthetic and must be labelled as such in the report.

1. **Phase sweep.** Combined floor amplitude for δ = 0°, 30°, …, 360°. **Prediction to write down first:** the curve is (3), with a minimum of 0.30 at 180° and a maximum of 2.10 at 0° and 360°; check that the sum's phase returns to 0 at δ = 180°.
2. **Mount sweep.** Tray acceleration (T × 2.10 m/s²) for mount natural frequencies 10, 15, 20, 25, 30, 35, 40, 45 Hz at ζ = 0.08. Predict first where T = 1 (28.3 Hz) and where T peaks (near 40 Hz, at 6.25).
3. **Transient convergence.** For the 30 Hz mount driven by 5 N at 40 Hz from rest, fit `fit_sinusoid` over successive 0.05 s windows and record the fitted amplitude against window start time. Predict first the time at which it is within 1 % of 69.79 µm (e^{−ζω_n t} = 0.01 gives t = 0.31 s, modified by the transient's own oscillation).

Raw data format: three CSV tables, `delta_deg, amplitude`; `f_mount, T, tray_accel`; `t_window, A_fit`, each with a header naming the fixed parameters.

What you should find: the phase sweep follows (3) to 10⁻¹² (it is the same formula, which is the point: a phasor diagram and the cosine rule are one calculation); the mount sweep crosses T = 1 at 28.3 Hz and peaks at 39.5 Hz, so the 30 Hz mount sits on the amplifying side; the fitted amplitude oscillates about the steady value with period 1/|40 − 29.9| = 0.1 s (the beat between the forced and the free oscillation) and settles within 1 % after about 0.35 s.

## 5. Validate and improve

| quantity | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| e^{iπ} from the series | −1 | −1 + 3.3 × 10⁻¹⁶ i | 3 × 10⁻¹⁶ | rounding in 40 terms |
| combined amplitude at 90° | √(1.2² + 0.9²) = 1.5000 | 1.5000 | 0 | same formula |
| H at r = 4/3, ζ = 0.08 | 1/(1 − 1.7778 + 0.2133i) = −1.1958 − 0.3280i | same | 0 | — |
| tray amplitude, 5 N | 1.2399 × 5/88 826 = 69.79 µm | 69.794 (phasor), 69.795 (ODE fit) | 0.008 µm | transient e^{−15.1 t} at 0.6 s |
| T at r = 4/3 | 1.0225/0.8066 = 1.268 | 1.2678 | 0 | — |
| mount for T = 0.2 | r ≈ 2.5 (graph) | 15.84 Hz (r = 2.525) | reading accuracy | — |

Every discrepancy is rounding or the known transient. What the agreement does not show is whether a real cabin behaves this way: the floor acceleration at a seat is a sum over structural modes, the 40 Hz line sits on top of a 80 Hz blade-passing line and broadband noise, and the "phase between propellers" at one seat differs from that at another because the structural paths differ (problem 5). A synchrophaser setting found by measurement at one seat is a compromise for the cabin.

Improvements worth making: fit two frequencies at once to separate 40 and 80 Hz (extension); replace the single-mass tray by a two-mass model with the floor's own mode; use `fit_sinusoid` on a sliding window to track a slowly drifting phase (problem 3). Improvement worth not making: applying the phasor sum to two propellers that are *not* synchronised.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivation of (1) and (4); `phasor_calculator.png`; the three CSV tables; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a synchrophaser angle and a mount frequency, with the transmissibility achieved, the static sag, and the new resonance frequency.

### Adversarial self-checks

- **What if both amplitudes double?** Every entry in the phase table doubles; the phase of the sum is unchanged; T is unchanged (a ratio). What if the *frequency* doubles (80 Hz, the blade-passing line)? r = 2.667 and the 30 Hz mount now isolates: T = 0.177. Check both.
- **What if ζ → 0?** |H(1)| → ∞ and T(1) → ∞, but T(√2) = 1 still. What if ζ → ∞? T → 1 at every frequency: a rigid damper is a rigid link.
- **Which assumption fails first in the real world?** Single frequency: the two propellers drift unless synchronised, and the sum beats. Second: the linear mount.
- **What measurement could falsify the model?** Sweep the synchrophaser angle and record the tray amplitude: (3) predicts a cosine in δ with extremes 0.3 and 2.1. A measured curve with a different minimum, or a minimum not at 180°, means a third path (the fuselage structure) is contributing a phasor of its own.
- **Does the result survive a change in solver or sample rate?** The phasor answer has no solver; the `solve_ivp` check at rtol = 10⁻¹⁰ and 10⁻⁶ changes the fitted amplitude in the sixth figure. The fit is insensitive to sample rate as long as it covers whole cycles.
- **Which decision changes if the uncertainty doubles?** If a₂ were 0.9 ± 0.3 instead of ± 0.15, the anti-phase residual could be anywhere from 0 to 0.6 m/s², and the synchrophaser's benefit at this tray would be uncertain by a factor 2; the mount decision (T = 0.2 needs 15.8 Hz) does not change, because it depends on the frequency ratio and not on the amplitudes.

## Check yourself

1. Write 1 + i in polar form and compute (1 + i)⁸ by de Moivre.
2. State Euler's formula and name the three power series that prove it.
3. What are the three cube roots of unity, what is their sum, and why?
4. A signal is 3 cos(ωt − 30°). Give its phasor, the phasor of its derivative, and the phasor of its integral, for ω = 10 rad/s.
5. Two 40 Hz sinusoids of amplitudes 1.2 and 0.9 have phase difference 90°. What is the amplitude of their sum? At what phase difference is the sum 1.2?
6. Why may we take the real part of X e^{iωt} to get the solution of a linear ODE with real coefficients, and why does the trick fail for the kinetic energy ½mẋ²?
7. At r = 1 with ζ = 0.05, what are |H| and the phase lag? At r = 2 with ζ = 0.1, is the transmissibility above or below 1?
8. A mount isolates at 40 Hz with T = 0.2. The engine idles at 950 rpm. Why might the mount be a bad idea anyway?
