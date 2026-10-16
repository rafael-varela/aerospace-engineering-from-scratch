# P02.02  Second-Order Linear ODEs

> Guess e^{rt}, solve a quadratic, and the whole motion of a spring, a mass and a damper is two numbers: how fast it rings and how fast it dies.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the characteristic equation of m ẍ + c ẋ + k x = 0, classify the three damping regimes from its discriminant, and write the real closed-form solution in each, including the critical case by reduction of order. 2. Define natural frequency, damping ratio, damped frequency and logarithmic decrement, and extract the damping ratio from the peaks of a decaying record. 3. Implement and test a damped spring-mass module: closed forms in all regimes, energy, the Wronskian and Abel's identity, a reduction to a first-order system integrated by RK4, and peak-based damping estimation, with `scipy.integrate.solve_ivp` and `scipy.signal.find_peaks` as oracles. 4. Decide what damping ratio a landing-gear strut should have so that the wheel does not bounce and the airframe load is least, and state where the linear strut model fails. |
| **Prerequisites** | `P02.01` (Euler and RK4 steppers, first-order systems, `solve_ivp`), `P01.16` (eigenvalues of a 2×2 matrix; the damped oscillator as d**u**/dt = A**u**), `P01.05` (Euler's formula e^{iθ} = cos θ + i sin θ), `P00.03`, `P00.07`. Skill check: you can solve a quadratic with complex roots, differentiate e^{−σt} cos ωt, and read ln(x₁/x₂) off two numbers. |
| **Engineering question** | A light twin's main gear leg carries 1200 kg on a 150 kN/m strut and touches down at 3 m/s. How far does the strut compress, does the wheel leave the ground again, how hard is the airframe hit, and what damping ratio makes the best compromise? |
| **Assumptions and limits** | Single degree of freedom (the leg alone, airframe rigid, tyre rigid); linear spring and linear (viscous) damper with constant coefficients; weight fully on the gear from the instant of touchdown (no residual lift); motion about the static equilibrium. Real oleo struts have a gas spring (stiffening) and orifice damping proportional to v², the tyre adds a second, stiffer spring, and lift decays over the roll-out; the linear model gives the right regimes and orders of magnitude, not certification numbers. The model allows negative compression (tension), which a wheel on a runway cannot supply: a negative value is read as a bounce. |
| **Analytical work** | Characteristic equation and discriminant; the three regimes and their real solutions; Euler's formula to turn complex exponentials into damped sinusoids; reduction of order for the repeated root; superposition and the Wronskian with Abel's identity W = W₀e^{−ct/m}; energy balance dE/dt = −cẋ²; logarithmic decrement δ = 2πζ/√(1 − ζ²) and its inversion; the coordinate shift that absorbs weight into the equilibrium; the first-peak time of the touchdown response. |
| **Computational deliverable** | `starter/main.py`: `natural_frequency`, `damping_ratio`, `damper_for_ratio`, `damped_frequency`, `characteristic_roots`, `classify_damping`, `free_response`, `mechanical_energy`, `wronskian`, `to_first_order`, `rk4_step`, `rk4_free_response`, `find_peaks_simple`, `logarithmic_decrement`, `zeta_from_decrement`, `settling_time`, `synthetic_decay_record`, `estimate_damping_from_record`, `touchdown_response`, `max_stroke`, `strut_force`, `loses_contact`. Tests prove: the closed form satisfies the ODE and the initial conditions in every regime; the three regimes join continuously at ζ = 1; solutions superpose; undamped energy is conserved and damped energy falls at exactly −cẋ²; the Wronskian follows Abel's identity; RK4 converges at order 4 to the closed form and `solve_ivp` agrees to 10⁻⁸ m; peaks match SciPy's and the decrement recovers ζ to 2 %; the seeded synthetic record returns ζ within 5 %; the touchdown stroke is the true maximum, equals the energy formula when undamped, and falls with ζ; the strut force starts at c v₀ and ends at m g. |
| **Practical activity** | **computational experiment** (a damping-ratio sweep of stroke, settling time, bounce and peak airframe load) and **data-based investigation** of a **synthetic** accelerometer record of a decaying panel mode: how the noise level and the peak threshold bias the estimated damping ratio. All data in this lesson are simulated. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `landing_gear_response.png`, the hand derivation of the three regimes and of Abel's identity, the damping-ratio sweep table, the synthetic-record table with your explanation of the bias, and a short note recommending a damping ratio for the strut with its reasons and the model's limits. |
| **Assessment** | `assessments/phase-02/grading-rubric.md`; competencies 1–4. Pass: the derivation covers all three regimes and the repeated root by reduction of order; all public tests pass; hand values for the baseline strut (ω_n, ζ, stroke) agree with the code to three figures; the note gives a damping ratio with a numerical justification and names the assumption it trusts least. |
| **References** | `@openstax-calculus-volume-3` §7.1 "Second-Order Linear Equations", §7.3 "Applications" (required-free). `@mit-ocw-18-03-differential-equations` Lecture 9 "Solving Second-order Linear ODE's with Constant Coefficients", Lecture 10 "Continuation: Complex Characteristic Roots", Lecture 11 "Theory of General Second-order Linear Homogeneous ODEs" (required-free); notes 11 "Modes and the characteristic polynomial" and 12 "Good vibrations, damping conditions" (supplementary-free). |
| **Extensions** | ★ Replace the linear damper by orifice damping F = c₂ ẋ|ẋ| and integrate with RK4; show the peak airframe load at touchdown is lower for the same stroke. ★ Add the tyre as a second spring in series (two degrees of freedom) and find both natural frequencies with `P01.16`. ★ Let the lift decay as L = mg e^{−t/τ} during the roll-out and show how the static equilibrium itself moves. |
| **Time** | unmeasured |

## 1. The engineering problem

When an aircraft touches down, the kinetic energy of its descent, ½mv² for the mass on each gear leg, has to go somewhere in a fraction of a second. The **shock strut** turns it into spring energy (a column of compressed gas or a steel spring) and then into heat (oil forced through an orifice). If the strut is too soft and lightly damped, the aircraft bounces: the wheel leaves the runway and comes down again, which is how tyres burst and nose gears collapse. If it is too stiff or too heavily damped, the strut transmits a hammer blow to the wing attachments. Every landing-gear design is a compromise between stroke, load and bounce, and the simplest model that shows the compromise is a mass on a spring with a damper: **m ẍ + c ẋ + k x = 0** about the static equilibrium.

That equation is the most important second-order linear ODE in engineering. The same equation, with the same solutions, describes a vibrating panel, an accelerometer's proof mass, the short-period pitching of an aircraft, the voltage in an RLC circuit and the nutation of a spinning satellite. Learn it once, with its three regimes and its two dimensionless numbers, and all of these become the same calculation.

The concrete problem: a leg of a light twin carries m = 1200 kg on a strut of stiffness k = 150 kN/m and touches down with a sink rate of 3.0 m/s. For damping ratios from 0 to 2, find the maximum compression, whether the wheel bounces, how fast the motion settles, and the peak force into the airframe; recommend a damping ratio. Success criteria, written first: the closed-form responses satisfy the equation to 10⁻⁶ relative by finite differences; our RK4 integration agrees with them to 10⁻⁸ m with 1200 steps over 3 s and `solve_ivp` to 10⁻¹¹ m; the damping ratio recovered from the peaks of a decaying record is within 2 % of the one that generated it; and the undamped stroke equals the energy-balance value exactly.

## 2. Understand and derive

### 2.1 The equation and its two numbers

**Assumption 1 (linear spring, linear damper):** spring force −kx, damper force −cẋ. **Assumption 2 (one degree of freedom):** the leg is a point mass m moving vertically. Newton's second law for the displacement x from the position where the spring force balances the weight gives

**m ẍ + c ẋ + k x = 0**.   (1)

(Weight disappears from (1) because we measure x from the *static* equilibrium x_s = mg/k; 2.7 does the bookkeeping.) The equation is **second order** (highest derivative ẍ), **linear** (x and its derivatives appear to the first power, no products) and **homogeneous** (no forcing term). Divide by m and define

**ω_n = √(k/m)**  (the **natural frequency**, rad/s),   **ζ = c/(2√(km))**  (the **damping ratio**, dimensionless),   (2)

so that (1) reads ẍ + 2ζω_n ẋ + ω_n² x = 0. Every spring-mass-damper in the world is specified by these two numbers and a scale. For the strut: ω_n = √(150 000/1200) = 11.18 rad/s = 1.779 Hz, period 0.562 s; the static deflection is x_s = 1200 × 9.81/150 000 = 78.5 mm; and c = 2ζ√(km) = 2ζ × 13 416 N s/m, so ζ = 0.7 needs c = 18.8 kN s/m. Units of c: N/(m/s) = kg/s; √(km) = √(N/m × kg) = √(kg² s⁻²) = kg/s ✓.

### 2.2 The exponential guess and the characteristic equation

Because (1) has constant coefficients, try x = e^{rt} (`@openstax-calculus-volume-3` §7.1; 18.03 Lecture 9). Then ẋ = r e^{rt}, ẍ = r² e^{rt}, and (1) becomes (m r² + c r + k) e^{rt} = 0. The exponential is never zero, so

**m r² + c r + k = 0**,   r = [−c ± √(c² − 4km)]/(2m) = ω_n(−ζ ± √(ζ² − 1)).   (3)

This is the **characteristic equation**, and its **discriminant** c² − 4km, or equivalently the sign of ζ − 1, sorts the physics into three regimes. Note the product of the roots is k/m = ω_n² and their sum is −c/m = −2ζω_n, whatever the regime; the tests check both.

Linearity and homogeneity give **superposition**: if x₁ and x₂ solve (1), so does Ax₁ + Bx₂. Two independent solutions therefore give the general one, and the two constants are fixed by the two initial conditions x(0) = x₀, ẋ(0) = v₀. "Independent" is checked by the **Wronskian** W = x₁ẋ₂ − ẋ₁x₂ (Lecture 11): differentiate it, use (1) to replace ẍ, and the x₁x₂ terms cancel, leaving

**dW/dt = −(c/m) W  ⟹  W(t) = W(0) e^{−ct/m}**   (Abel's identity).   (4)

W is either always zero or never zero, so independence at t = 0 is independence forever; and the decay rate of W is exactly the damping c/m. For the fundamental pair x₁ (x₀ = 1, v₀ = 0) and x₂ (x₀ = 0, v₀ = 1), W(0) = 1, and the code's `wronskian` reproduces e^{−ct/m} to 10⁻⁸ in every regime.

### 2.3 Underdamped, ζ < 1: a decaying sinusoid

The roots are complex, r = −ζω_n ± iω_d with the **damped frequency** **ω_d = ω_n√(1 − ζ²)** (Lecture 10). The complex solutions e^{(−ζω_n ± iω_d)t} combine, through Euler's formula e^{iθ} = cos θ + i sin θ (`P01.05`), into the real pair e^{−ζω_n t} cos ω_d t and e^{−ζω_n t} sin ω_d t. Fitting x₀ and v₀:

**x(t) = e^{−ζω_n t} [x₀ cos ω_d t + ((v₀ + ζω_n x₀)/ω_d) sin ω_d t]**.   (5)

The motion oscillates at ω_d, slightly below ω_n, inside an envelope e^{−ζω_n t} with time constant 1/(ζω_n). For the strut at ζ = 0.2: roots −2.24 ± 10.95i s⁻¹, ω_d = 10.95 rad/s, period 0.574 s, time constant 0.447 s.

**Logarithmic decrement.** Successive positive peaks are one damped period T_d = 2π/ω_d apart, so their ratio is e^{ζω_n T_d}, and

**δ = ln(x_i/x_{i+1}) = 2πζ/√(1 − ζ²)**,   inverted exactly as   **ζ = δ/√(4π² + δ²)**.   (6)

This is how damping is *measured*: strike the structure, record the decay, read the peaks. For ζ = 0.2, δ = 1.283, and peaks shrink by a factor 3.6 per cycle; the code reads 5 positive peaks off the ζ = 0.2 touchdown record, finds δ = 1.2826, and returns ζ = 0.2000 and a period of 0.5735 s against 2π/ω_d = 0.5736 s.

### 2.4 Critically damped, ζ = 1: the repeated root

The discriminant is zero and both roots are r = −ω_n. One solution is e^{−ω_n t}; the exponential guess gives no second one. **Reduction of order** finds it: try x = u(t) e^{−ω_n t}, substitute into ẍ + 2ω_n ẋ + ω_n² x = 0, and everything cancels except ü e^{−ω_n t} = 0, so u = A + Bt. The second solution is t e^{−ω_n t}, and

**x(t) = e^{−ω_n t} [x₀ + (v₀ + ω_n x₀) t]**.   (7)

No oscillation, and the fastest return without overshoot for a *displaced* start with v₀ = 0 (a result the strut problem complicates, because v₀ ≠ 0). The three formulas must agree as ζ → 1 from either side; the test compares (5) and (8) at ζ = 1 ∓ 10⁻⁶ with (7) and finds them within 10⁻⁶ m, which is also the reason `free_response` treats |ζ − 1| < 10⁻⁹ as critical rather than letting ω_d → 0 divide by zero.

### 2.5 Overdamped, ζ > 1: two decaying exponentials

Both roots are real and negative, r₁ < r₂ < 0, and

**x(t) = A e^{r₁t} + B e^{r₂t}**,  A + B = x₀,  A r₁ + B r₂ = v₀.   (8)

The slow root r₂ = ω_n(−ζ + √(ζ² − 1)) dominates after a short time. For ζ = 2 the roots are −41.7 and −3.00 s⁻¹: the fast mode is gone in a tenth of a second and the strut then creeps toward equilibrium with a 0.33 s time constant, slower than the critically damped case. More damping is not faster.

### 2.6 Energy

Multiply (1) by ẋ: m ẋ ẍ + c ẋ² + k x ẋ = 0, i.e.

**d/dt [½mẋ² + ½kx²] = −c ẋ²**.   (9)

The **mechanical energy** E is conserved when c = 0 and otherwise decreases monotonically at exactly the rate the damper dissipates power. The test checks both statements numerically. At touchdown the strut holds E = ½ × 1200 × 3² + ½ × 150 000 × 0.0785² = 5400 + 462 = 5862 J (the second term is the spring energy of the *unloaded* strut relative to equilibrium); after 3 s at ζ = 0.2, 0.007 J remains.

Energy also gives the **undamped stroke** without solving anything: all of E ends up in the spring at maximum compression, so ½k y_max² = E, y_max = √(x_s² + mv₀²/k) = 0.2795 m about equilibrium, and the compression from the unloaded length is x_s + y_max = **358 mm**. The code's `max_stroke` at ζ = 0 returns 358.0 mm, which is the test's exact check of the root-finder.

### 2.7 The landing gear: shifting to the static equilibrium

Measure the compression X from the *unloaded* strut. Newton: m Ẍ = mg − kX − cẊ, a nonhomogeneous equation. Put X = x_s + x with x_s = mg/k: the constant terms cancel and x obeys (1). The initial conditions become x(0) = −x_s (the strut starts 78.5 mm *above* its equilibrium) and ẋ(0) = v₀ = 3.0 m/s. The wheel **bounces** if X ever goes negative, because the runway cannot pull on the tyre; the linear model happily reports negative X and we read it as a bounce flag. The **strut force** on the airframe is F = kX + cẊ; at t = 0 the spring is slack and F = c v₀, the damper alone: at ζ = 0.7 that is 56.3 kN before the strut has moved a millimetre, 4.8 times the weight. The maximum compression is where Ẋ = 0, which for the underdamped case happens at tan(ω_d t*) = (−σx₀ + bω_d)/(σb + x₀ω_d) with σ = ζω_n, b = (v₀ + σx₀)/ω_d (problem 6 derives it); the code brackets and solves it with `brentq`.

The results of `main()`:

```
 zeta  c [kN s/m]       regime                  roots [1/s]  max stroke [mm]  at t [s]  bounce?  settle 2% [s]  peak force [kN]
 0.00        0.00     undamped    -0.00-11.18j, 0.00+11.18j            358.0     0.166     True            inf             53.7
 0.20        5.37  underdamped   -2.24-10.95j, -2.24+10.95j            266.5     0.152     True          2.229             42.4
 0.50       13.42  underdamped     -5.59-9.68j, -5.59+9.68j            189.0     0.138    False          0.897             42.1
 0.70       18.78  underdamped     -7.83-7.98j, -7.83+7.98j            157.1     0.132    False          0.667             56.3
 1.00       26.83     critical               -11.18, -11.18            124.7     0.126    False          0.593             80.5
 2.00       53.67   overdamped                -41.73, -3.00             78.5       inf    False          0.504            161.0
```

Reading the table: the stroke falls steadily with damping (358 → 78.5 mm); the wheel bounces for ζ below about 0.25 (a sweep finds the threshold at ζ = 0.245); the settling time (to within 2 % of x_s) has a shallow minimum near ζ ≈ 0.78 at 0.45 s; and the peak airframe force has a *minimum* of 39.8 kN at ζ = 0.36, rising to 80 kN at critical damping and 161 kN at ζ = 2 because the damper force cv₀ at the instant of touchdown grows without bound. The engineering answer is ζ between 0.3 and 0.5: no bounce, peak load within 2 % of its minimum, stroke 190–235 mm. Critical damping, the textbook "ideal", is the wrong answer for a landing gear, because the gear is not released from rest: it arrives with velocity, and a viscous damper punishes velocity.

### 2.8 Where the model stops being valid

- **Linear damper.** Real oleo struts meter oil through an orifice: F ∝ ẋ|ẋ|, zero at zero speed and quadratic at high speed, which is exactly why they avoid the touchdown spike of the linear model. Extension 1.
- **Linear spring.** The gas spring stiffens as it compresses (polytropic), so the stroke is less and the force at full stroke more than (5) predicts.
- **Rigid tyre.** The tyre is a softer, undamped spring in series with the strut; it adds a second, higher-frequency mode (two degrees of freedom, `P01.16`) and does most of the bouncing in practice.
- **Full weight from t = 0.** The wing is still lifting at touchdown; the equilibrium x_s moves during the roll-out.
- **One leg.** The aircraft pitches as the mains and the nose gear load up at different times.

## 3. Build the computational model

Open `starter/main.py`. `landing_gear()` is given; everything else is yours.

### Step 1: parameters, roots, regimes

`natural_frequency`, `damping_ratio`, `damper_for_ratio` and `damped_frequency` are (2). `characteristic_roots` solves (3) in complex arithmetic (so a negative discriminant is not an error) and returns the pair sorted; `classify_damping` compares ζ with 1 to a tolerance of 10⁻⁹.

### Step 2: the closed form in every regime

```python
def free_response(t, m, c, k, x0, v0):
    t = np.asarray(t, dtype=float)
    wn, z = natural_frequency(m, k), damping_ratio(m, c, k)
    if abs(z - 1.0) < CRITICAL_TOL:                                    # (7)
        a = v0 + wn * x0
        e = np.exp(-wn * t)
        return e * (x0 + a * t), e * (a - wn * (x0 + a * t))
    if z < 1.0:                                                        # (5)
        wd, sig = wn * np.sqrt(1.0 - z * z), z * wn
        b = (v0 + sig * x0) / wd
        e, cs, sn = np.exp(-sig * t), np.cos(wd * t), np.sin(wd * t)
        return e * (x0 * cs + b * sn), e * ((-sig * x0 + b * wd) * cs + (-sig * b - x0 * wd) * sn)
    r1, r2 = characteristic_roots(m, c, k).real                        # (8)
    B = (v0 - r1 * x0) / (r2 - r1); A = x0 - B
    return A * np.exp(r1 * t) + B * np.exp(r2 * t), A * r1 * np.exp(r1 * t) + B * r2 * np.exp(r2 * t)
```

Return the velocity as well as the position: the energy, the Wronskian and the strut force all need it, and the test differentiates x numerically to check that the two agree. `mechanical_energy` is (9)'s bracket; `wronskian` builds the fundamental pair and returns x₁ẋ₂ − ẋ₁x₂.

### Step 3: first-order form and RK4

`to_first_order(m, c, k)` returns A = [[0, 1], [−k/m, −c/m]], so that d/dt (x, ẋ) = A (x, ẋ). Its eigenvalues are the characteristic roots (the test checks this against `np.linalg.eigvals`): a second-order equation *is* a 2×2 first-order system, which is how `solve_ivp` and every general-purpose integrator see it. `rk4_step` is the stepper from `P02.01`; `rk4_free_response` loops it on a uniform grid.

### Step 4: reading a decaying record

`find_peaks_simple` returns strict interior local maxima; with `min_separation > 0` it accepts candidates in order of height and drops any within that distance of an accepted peak (SciPy's `distance` rule). `logarithmic_decrement` uses the first and last peak, ln(p₀/p_N)/N, which averages over the record; `zeta_from_decrement` is the exact inversion in (6). `settling_time` returns the first time after which |x| stays inside a band. `synthetic_decay_record` makes a labelled synthetic accelerometer trace with seeded Gaussian noise, and `estimate_damping_from_record` turns it into f_d (median peak spacing) and two ζ estimates: first/last decrement and a least-squares slope of ln(peak) against time.

### Step 5: the landing gear

`touchdown_response` applies the shift of 2.7; `max_stroke` brackets the first zero of the velocity and refines it with `brentq`, returning (∞, mg/k) if the velocity never changes sign (heavily overdamped); `strut_force` is kX + cẊ; `loses_contact` flags any negative compression.

### Verify

Run `python -m pytest curriculum/phase-02/lesson-02-second-order-linear-odes`. The tests:

1. **Analytic cases.** For ζ = 0, 0.15, 0.7, 1.0, 1.6 the closed form satisfies m ẍ + c ẋ + k x = 0 to 10⁻⁶ relative (ẍ by central differences) and meets x(0) = x₀, ẋ(0) = v₀; the returned velocity is the derivative of the returned position; the roots match `np.roots` and obey Vieta; the real part is −ζω_n and the imaginary part ω_d.
2. **Limiting cases.** The three regimes agree to 10⁻⁶ m at ζ = 1 ± 10⁻⁶; c = 0 conserves energy to 10⁻¹⁰; ζ ≥ 1 gives ω_d = 0; the undamped stroke equals √(x_s² + mv₀²/k) + x_s; at ζ = 2 the stroke is the static deflection, reached only asymptotically.
3. **Conservation and structure.** dE/dt = −cẋ² along the response; W(t) = e^{−ct/m}; solutions superpose and scale; the eigenvalues of the state matrix are the roots.
4. **Convergence and oracle.** RK4 errors fall 16-fold per halving (order 4.0 ± 0.3) and reach 10⁻⁸ m at 1200 steps; `solve_ivp` at rtol 10⁻¹⁰ agrees with the closed form to 10⁻⁸ m.
5. **Measurement.** Peaks coincide with `scipy.signal.find_peaks`; the decrement returns ζ within 2 % for ζ = 0.05, 0.2, 0.4; the separated finder rejects noise ripples and lands on the true peaks; the seeded synthetic record returns ζ = 0.015 within 5 % and f_d within 2 %.
6. **Gear.** The response starts unloaded at the sink rate and relaxes to mg/k; the stroke is the true maximum and decreases with ζ; the force starts at cv₀, ends at mg, and satisfies m Ẍ = mg − F along the way; ζ = 0.1 bounces, ζ = 0.6 does not.

### Compare

`solve_ivp` (RK45, rtol 10⁻¹⁰) reproduces the ζ = 0.2 closed form to 8 × 10⁻¹² m; our fixed-step RK4 reaches 6.5 × 10⁻⁷ m with 300 steps and 1.6 × 10⁻¹⁰ m with 2400, ratios of 16.1, 16.1, 16.0. Both are far below any physical uncertainty, and the comparison teaches something else: for a *linear* constant-coefficient equation, stepping is the wrong tool. The closed form is exact, costs nothing, and is what `scipy.signal` uses internally for linear systems. We step anyway because the moment the damper becomes F ∝ ẋ|ẋ| (extension 1) or the spring stiffens, the closed form is gone and the stepper is all that is left; the closed form is then the test case that tells you the stepper is right before you trust it on the nonlinear problem. For the peak-finding, `scipy.signal.find_peaks` is the oracle, and the lesson's `min_separation` is its `distance` argument.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment** (sweep) and **data-based investigation** (synthetic record). Write each prediction down first.

1. **Damping sweep.** For ζ from 0.05 to 3 in 300 steps, compute the maximum stroke, the 2 % settling time, the bounce flag and the peak strut force; plot all three curves against ζ (`landing_gear_response.png`, right panel). Predictions to write first: where is the bounce threshold; does the settling time have a minimum and where; does the peak force have a minimum, and why is it not at ζ = 1.
2. **Sink-rate scaling.** Repeat the stroke at ζ = 0.4 for v₀ = 1, 2, 3, 4 m/s. Predict the dependence first (the undamped energy argument gives √(x_s² + mv₀²/k) + x_s; how far below it does damping put the real stroke?).
3. **Synthetic record.** `synthetic_decay_record(42.0, 0.015, 25.0, 2.0, 2000.0, noise, seed=202)` makes a **simulated** accelerometer trace of a 42 Hz panel mode with ζ = 0.015 (label it synthetic in every table). For noise 0.1, 0.5 and 1.0 m/s² and peak thresholds of 5σ and 20σ, run `estimate_damping_from_record` and tabulate the number of peaks, f_d and both ζ estimates. Predict first: does noise bias ζ up or down, and which threshold is better.
4. **Record length.** With noise 0.5 m/s² and a 20σ threshold, cut the record to 0.25, 0.5 and 1.0 s and repeat. Predict first whether the estimate changes (hint: where in the record do the peaks above 20σ lie?).

Raw data format: CSV tables (`zeta, stroke_m, settle_s, bounce, peak_force_N`; `v0, stroke_m`; `noise, threshold, n_peaks, f_d, zeta_decrement, zeta_fit`), each with a header naming the fixed parameters and the word SYNTHETIC on the record tables.

What you should find: a bounce threshold at ζ = 0.245; a settling-time minimum near ζ = 0.78; a peak-force minimum at ζ = 0.36 (39.8 kN) with the force rising linearly in ζ beyond about 0.6 because the touchdown spike cv₀ takes over; the stroke growing almost linearly with sink rate (116, 160, 210, 262 mm at ζ = 0.4 for 1–4 m/s), well below the undamped energy formula because the damper removes energy on the way in; and in the synthetic study, a *downward* bias in ζ that grows with noise and shrinks with the threshold (table in part 5), because noise lifts small peaks more than large ones and so flattens the apparent decay.

## 5. Validate and improve

| quantity (strut, ζ = 0.2 unless stated) | hand (§2) | code | difference | cause |
|---|---|---|---|---|
| ω_n (rad/s) | 11.18 | 11.180 | 0 | exact |
| static deflection (mm) | 78.5 | 78.5 | 0 | exact |
| roots (s⁻¹) | −2.24 ± 10.95i | −2.24 ± 10.95i | 0 | exact |
| undamped stroke (mm) | 358 (energy) | 358.0 (`brentq` on ẋ) | 0 | exact |
| δ from peaks | 1.283 (formula) | 1.2826 (5 peaks) | 0.03 % | peak on a 1 ms grid |
| damped period (s) | 0.5736 | 0.5735 | 0.02 % | same |
| W(1 s) | e^{−c/m} = 0.011423 | 0.011423 | 0 | exact |
| RK4 error, 1200 steps (m) | ∝ h⁴ | 2.5 × 10⁻⁹ | — | order 4.0 measured |
| `solve_ivp` error (m) | — | 8 × 10⁻¹² | — | rtol 10⁻¹⁰ |

The synthetic-record investigation (all rows simulated, seed 202, true ζ = 0.015, f_d = 41.995 Hz):

| noise (m/s²) | threshold | peaks | f_d (Hz) | ζ first/last | ζ fit |
|---|---|---|---|---|---|
| 0.1 | 5σ | 45 | 41.67 | 0.0137 | 0.0140 |
| 0.1 | 20σ | 27 | 42.11 | 0.0148 | 0.0148 |
| 0.5 | 5σ | 28 | 42.55 | 0.0129 | 0.0127 |
| 0.5 | 20σ | 10 | 42.55 | 0.0145 | 0.0144 |
| 1.0 | 5σ | 21 | 42.55 | 0.0090 | 0.0119 |
| 1.0 | 20σ | 4 | 44.44 | 0.0023 | 0.0037 |

Every discrepancy has a cause. The low-threshold rows are biased low by 7–40 %: a peak of the noisy signal is the maximum of signal plus noise over a neighbourhood, which exceeds the true peak by a few tenths of σ on average, and that fixed offset is a larger *fraction* of the small late peaks than of the large early ones, so the fitted decay is too slow. Raising the threshold to 20σ removes the bias (−1 % to −4 %) while there are still ten or more peaks, and fails when the noise is 4 % of the initial amplitude because only four peaks survive and the first/last estimate then rests on two numbers. The least-squares fit degrades more gracefully than the first/last decrement, as it should: it uses all the peaks. The rule you should extract is quantitative: keep peaks above about 20σ and need at least ten of them, which fixes the record length and the noise floor a damping test needs. On a real structure a further bias appears that no threshold removes: a second mode close in frequency beats against the first and modulates the peaks, which is why `P02.10` moves this measurement into the frequency domain.

For the gear, the numbers are self-consistent but validated against nothing physical. What would validate them is a **drop test**: the standard certification test drops a complete gear with its share of the aircraft mass from a height that gives the design sink rate onto a load platform and records stroke and platform force against time. Comparing that record with `touchdown_response` would show the linear model's two characteristic errors at once: the measured force has no spike at t = 0 (orifice damping is zero at zero velocity) and the measured stroke is shorter near full compression (gas spring stiffening). Both are improvements you can make in the extensions, and both are made by replacing a closed form with the RK4 stepper, which is exactly why the stepper was built and tested here.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand derivation of (3), (5), (7), (8) and Abel's identity (4); `landing_gear_response.png`; the sweep and synthetic-record tables with your predictions beside them; and a one-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a damping ratio, giving the stroke, bounce margin and peak load it implies, and naming the assumption of the strut model you trust least.

### Adversarial self-checks

- **What if the mass doubles (a heavier aircraft on the same gear)?** ω_n falls by √2, x_s doubles, and for the same c the damping ratio falls by √2: a strut tuned to ζ = 0.4 at 1200 kg is at ζ = 0.28 at 2400 kg, just above the bounce threshold. Check with the code; this is why gear is designed for the maximum landing weight.
- **What if the sink rate goes to zero?** The strut is released from 78.5 mm above equilibrium with no velocity; the critically damped case is then genuinely the fastest non-overshooting return and the peak force is kx_s = mg for every ζ. The whole trade-off of 2.7 came from v₀ ≠ 0.
- **Which assumption fails first in the real world?** The linear damper. The touchdown force spike cv₀ does not exist in an oleo strut; the real peak force comes near maximum stroke from the stiffening gas spring.
- **What measurement could falsify the model?** A drop test's force record: a linear model predicts F(0) = cv₀ > 0 and a maximum force at the first velocity zero; an orifice strut gives F(0) = 0 and a later maximum.
- **Does the result survive a change of solver?** Closed form, RK4 and `solve_ivp` agree to 10⁻⁸ m. Does the *recommendation* survive a change of criterion? Minimising peak force gives ζ = 0.36; minimising settling time gives ζ = 0.78 at twice the load. State which you would choose for a passenger aircraft and for a carrier aircraft.
- **Which decision changes if the uncertainty doubles?** If the sink rate at touchdown could be 4 m/s rather than 3 m/s (a hard landing), the stroke at ζ = 0.4 grows from 210 to 262 mm and the touchdown kinetic energy by 78 %; a strut with 250 mm of travel would bottom out, which the linear model cannot show but the energy balance can.

## Check yourself

1. Write m ẍ + c ẋ + k x = 0 in terms of ω_n and ζ, and give the units of c and of √(km).
2. What are the roots of the characteristic equation when ζ = 0.6, in terms of ω_n? What is the damped frequency as a fraction of ω_n?
3. Why does the critically damped case need the solution t e^{−ω_n t}, and how is it found?
4. State Abel's identity for the Wronskian of (1) and explain in one sentence why two solutions independent at t = 0 remain independent.
5. The peaks of a decaying record are 12.0, 8.5, 6.0, 4.3 mm. Find δ and ζ.
6. Multiply (1) by ẋ and derive the energy equation. What does it say when c = 0?
7. Why is the peak airframe force at touchdown *not* minimised by critical damping, and roughly where is it minimised for the lesson's strut?
8. A synthetic record with noise returned ζ = 0.0129 for a true 0.015 at a 5σ peak threshold. Which way did the noise bias the estimate, why, and what fixed it?
