# P09.01  Signals & Systems

> A system is a rule that turns one signal into another. Five yes-or-no questions about that rule (linear? time-invariant? causal? memoryless? stable?) decide which of the tools in the rest of this phase you are allowed to use.

| | |
|---|---|
| **Learning outcomes** | 1. Define and compute the energy and average power of a signal, split a signal into even and odd parts and prove they are orthogonal, and apply time shift, time scale and time reversal in the right order. 2. State the five system properties (linearity, time invariance, causality, memory, BIBO stability) precisely, and derive the verdict for a delay, a squarer, an accumulator, a time-varying gain and seven other blocks by hand. 3. Implement and test property testers that try to falsify each property on a sampled discrete-time operator, and show that they reproduce the hand verdicts. 4. Prove that an LTI system is BIBO stable if and only if Σ\|h[n]\| is finite, and explain, with numbers, how a numerical test can give a false verdict (amplitude, record length, untested inputs). |
| **Prerequisites** | `P02.08` (Laplace transform, impulse and step responses), `P02.10` (sampled signals, the DFT), `P01.13` (linear maps, superposition), `P00.07` (series, geometric sums), `P00.03` (NumPy arrays). Skill check: you can sum a geometric series and say what makes a map f linear. |
| **Engineering question** | The rate signal of a MEMS gyro in a flight computer passes through a delay, a bias correction, a converter that saturates, an integrator to angle, moving averages, a recursive low-pass, a centred smoother and a decimator. Which of these blocks are linear and time-invariant (so that impulse responses, transfer functions and Bode plots apply to them), which are causal (so they can run in real time), and which are BIBO stable? And how far can a numerical test of these properties be trusted? |
| **Assumptions and limits** | Discrete-time systems acting on finite records x[0..N−1] that start from rest (x[n] = 0 and all internal states zero for n < 0); continuous signals are represented by samples and their integrals by Riemann sums. The numerical property tests are *falsification* tests: one counterexample proves that a property fails, but passing every trial is evidence, not proof. Their verdicts depend on the test amplitude (a saturation looks linear below its limit), on the record length (a slow but stable filter looks unstable on a short record) and on which inputs are tried. |
| **Analytical work** | Energy of e^{−at}u(t) and power of a sinusoid; the even/odd split and its orthogonality; x(at − b) by two orders of operations; the five property definitions; hand verdicts for eleven blocks; superposition for the affine bias correction; the BIBO theorem for LTI systems (sufficiency by the triangle inequality, necessity by the worst-case input sign h[n − k]); the impulse response of y[n] = a y[n−1] + b x[n]. |
| **Computational deliverable** | `starter/main.py`: `unit_step`, `ramp`, `rect_pulse`, `unit_impulse`, `energy`, `average_power`, `even_part`, `odd_part`, `time_shift`, `time_scale`, `time_reverse`, `superposition_error`, `is_linear`, `is_time_invariant`, `is_causal`, `is_memoryless`, `is_bibo_stable`, `impulse_response`, `abs_sum`, `classify` (`sinusoid`, `shift_signal`, `bounded_test_inputs`, the gyro-path zoo and `main` are given). Tests prove: the ramp is the running integral of the step; the energy of e^{−at} converges to 1/(2a) at first order in dt; a sinusoid has power A²/2 whatever its frequency and phase; the even/odd parts are exact, symmetric and orthogonal; x(2t − 3) comes out the same by both orders of operations; the property testers reproduce the hand verdict for all eleven zoo systems; the saturation is linear inside its limits and not beyond; linear systems map zero to zero and the affine one does not; the recursion's impulse response is b aⁿ; every causal LTI zoo system equals convolution with its impulse response; the bound \|y\| ≤ Σ\|h\| max\|x\| holds and is reached by the worst-case input; the accumulator's Σ\|h\| grows without bound; the stability verdict is wrong on a record shorter than the time constant. |
| **Practical activity** | **computational experiment**: classify the gyro-path zoo and compare with the hand verdicts; sweep the test amplitude of the saturation; sweep the record length of the stability test for three time constants; integrate a gyro bias for an hour. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `signals_and_systems.png`, the hand verdict table for the zoo with one line of reasoning per entry, the BIBO proof, the two sweep tables, and a short note stating which blocks of the gyro path may be analysed with transfer functions and what each numerical verdict does and does not prove. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the BIBO theorem and three hand verdicts (accumulator, decimator, bias correction) are derived; all public tests pass; the amplitude and record-length failures of the numerical tests are explained with numbers. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 1 "Signals and Systems" (The Signals and Systems Abstraction; Example System: Leaky Tank; Analysis of the Leaky Tank), Lecture 2 "Discrete-Time (DT) Systems" (Difference Equations; Operator Notation; Operator Algebra; Example: Accumulator; Finite and Infinite Impulse Responses), Lecture 4 "Continuous-Time (CT) Systems" (Elementary Building-Block Signals; Unit-Impulse and Unit-Step Signals) and Lecture 8 "Convolution" (Superposition; Linearity; Time-Invariance) (required-free): a system as an operator on whole signals, the delay and accumulator as the canonical discrete-time operators, and linearity plus time invariance as the two properties that make superposition work. Causality, memory, BIBO stability, signal energy and power, the even/odd decomposition, the numerical property tests on sampled operators and the analytic verdicts (including the time-varying gain and the squarer) are original instruction in this lesson. |
| **Extensions** | ★ Add an "incremental linearity" tester (S(x₁) − S(x₂) linear in x₁ − x₂) and show the bias correction passes it while the squarer fails. ★ Extend the zoo with a median filter, a rate limiter and a deadband, and classify them by hand first. ★ Replace the fixed test inputs of `is_bibo_stable` by the worst-case input sign(h[n − k]) for LTI systems and measure how much shorter a record it needs. |
| **Time** | unmeasured |

## 1. The engineering problem

A MEMS rate gyro measures pitch rate as a voltage. Before the number reaches the control law it passes through a chain of processing blocks written by different people at different times: a sample delay from the converter pipeline, a bias correction from the last calibration, a converter that clips at full scale, an integrator that turns rate into angle, one or two moving averages, a recursive low-pass, perhaps a centred smoother copied from a post-processing script, and a decimator that drops every second sample to feed a slower loop. The rest of this phase builds powerful tools for such chains: the impulse response and convolution (`P09.02`), transfer functions and block-diagram algebra (`P09.03`), Bode plots (`P09.05`), the Nyquist criterion (`P09.10`). Every one of them assumes that the block is **linear and time-invariant** (LTI). Apply them to a block that is not, and the answer is wrong without warning.

So the first job of a control engineer reviewing a signal path is to classify each block. Is it **linear** (does the response to a sum equal the sum of the responses)? **Time-invariant** (does a delayed input give a delayed output)? **Causal** (does the output depend only on present and past input, so the block can run in real time)? **Memoryless** (does it depend on the present input only)? **BIBO stable** (does every bounded input give a bounded output)? The last matters physically: a gyro bias of 0.01 deg/s fed through the integrator is a bounded input, and after one hour it is a 36 deg attitude error.

Success criteria, written first: for every block in the gyro-path zoo the numerical property tests return the verdict derived by hand; the energy and power routines converge to their closed forms; the LTI blocks reproduce themselves as convolutions with their impulse responses to 10⁻¹⁰; and the lesson states, with numbers, at least two ways in which the numerical tests can be fooled.

## 2. Understand and derive

### 2.1 Signals, energy and power

A **continuous-time signal** x(t) is a function of time; a **discrete-time signal** x[n] is a sequence, usually samples x[n] = x(nT) taken every T seconds (`P02.10`). The elementary signals are the **unit step** u(t) (0 for t < 0, 1 for t ≥ 0), the **ramp** r(t) = t u(t), whose derivative is the step, the **rectangular pulse** u(t) − u(t − w), the **discrete impulse** δ[n] (1 at n = 0, 0 elsewhere), and the sinusoid A cos(2πft + φ). Products such as e^{−σt} cos(ωt) u(t) describe the decaying oscillations that every stable system in this phase produces.

The **energy** of a signal is E = ∫|x(t)|² dt (Σ|x[n]|² for a sequence) and its **average power** is P = lim (1/T)∫_{−T/2}^{T/2} |x|² dt. A signal with finite energy has zero power (an **energy signal**: a pulse, a decaying transient); a signal with finite non-zero power has infinite energy (a **power signal**: a sinusoid, a step, noise). The units are those of x² times seconds and of x², whatever x is: for a rate signal in deg/s, energy is in deg²/s.

Two closed forms are worth deriving once. For x(t) = e^{−at}u(t), a > 0: E = ∫₀^∞ e^{−2at} dt = **1/(2a)**. For x(t) = A cos(ωt + φ): x² = (A²/2)(1 + cos(2ωt + 2φ)), and the cosine averages to zero over whole periods, so **P = A²/2** for every ω and φ. A Riemann sum with step dt overestimates the first by a relative amount close to a·dt (the left end point of each strip is the highest), so with a = 2 and dt = 10⁻⁴ s the code gives 0.25005 instead of 0.25; the power of 2 cos(2π·5t) over 10 s comes out 2.00000.

### 2.2 Even and odd parts

Any signal splits uniquely into an **even part** x_e(t) = [x(t) + x(−t)]/2 and an **odd part** x_o(t) = [x(t) − x(−t)]/2, with x = x_e + x_o. The two are **orthogonal**: x_e x_o is odd, so ∫x_e x_o dt = 0 over a symmetric interval, and therefore E(x) = E(x_e) + E(x_o). For x = e^{−t}u(t), x_e = ½e^{−|t|} and x_o = ½ sgn(t) e^{−|t|}, each with energy ¼, together ½. A flight-test **elevator doublet** (+1 deg for 1 s, then −1 deg for 1 s, centred on t = 0) is odd: on a 40 001-point grid over [−2, 2] s the code gives energy 2.0000 deg²·s, of which the odd part carries 1.9998 and the even part 1.5 × 10⁻⁴ (the residue sits at the three grid points where the pulse jumps). Odd inputs matter in flight testing because they excite the aircraft without changing its mean attitude.

### 2.3 Operations on time

Three operations act on the time axis: the **shift** x(t − t₀) (a delay for t₀ > 0: the signal moves right), the **scale** x(at) (compressed for |a| > 1) and the **reversal** x(−t). They do not commute, and the order matters in practice when a recorded manoeuvre is replayed faster and later. To build y(t) = x(at − b) either shift first by b and then scale by a, giving x(a t − b), or scale first by a and then shift by **b/a**, giving x(a(t − b/a)) = x(at − b). Shifting by b after scaling gives x(a(t − b)) = x(at − ab), a different signal. In code these are higher-order functions: `time_shift(f, t0)` returns the new function t ↦ f(t − t₀), and the test checks that x(2t − 3) comes out the same by both correct orders and differently by the wrong one.

### 2.4 Systems and their five properties

A **system** S maps an input signal x to an output signal y = S(x) (the abstraction of 6.003 Lecture 1; the delay and the accumulator of Lecture 2 are the canonical discrete-time examples). In this lesson a system is a sampled operator: a function that takes the array x[0..N−1] and returns y[0..N−1], with the convention of **initial rest** (everything is zero before n = 0). The properties are statements about S for *all* inputs:

- **Linear**: S(a x₁ + b x₂) = a S(x₁) + b S(x₂) for all inputs and all constants (superposition = additivity + homogeneity). Taking a = b = 0 shows a linear system maps zero to zero.
- **Time-invariant**: if y = S(x), then delaying the input by k samples delays the output by k samples and changes nothing else.
- **Causal**: y[n] depends only on x[m] for m ≤ n. Only causal systems can run in real time.
- **Memoryless** (static): y[n] depends only on x[n].
- **BIBO stable** (bounded input, bounded output): every input with |x[n]| ≤ M for all n gives an output with |y[n]| ≤ M′ < ∞ for all n.

For an LTI system everything is carried by the **impulse response** h = S(δ). Because x = Σ_k x[k] δ[n − k], linearity and time invariance give (6.003 Lecture 8) y[n] = Σ_k x[k] h[n − k]: the convolution of `P09.02`. Causality is then h[n] = 0 for n < 0, memorylessness is h = h[0]δ, and stability has a clean test:

**Theorem (BIBO stability of LTI systems).** An LTI system is BIBO stable if and only if Σ_n |h[n]| < ∞.

*Sufficiency.* If |x| ≤ M then |y[n]| = |Σ_k h[k] x[n − k]| ≤ Σ_k |h[k]| M, a bound independent of n.
*Necessity.* Suppose Σ|h| = ∞. For any N choose the bounded input x[m] = sign(h[N − m]) (|x| ≤ 1). Then y[N] = Σ_k h[k] sign(h[k]) = Σ_{k≤N} |h[k]|, which exceeds any bound as N grows. ∎

The worst-case input of the necessity proof also shows that the bound is reached: for y[n] = −0.8 y[n−1] + x[n], Σ|h| = 1/(1 − 0.8) = 5, and the input sign(h[N − m]) drives y[N] to 5 exactly (a public test checks this).

### 2.5 The gyro-path zoo, by hand

The verdicts below are derived from the definitions, not from the code; they are the independent check that the property testers must reproduce.

| block | rule | linear | TI | causal | memoryless | BIBO | reason for the less obvious entries |
|---|---|---|---|---|---|---|---|
| delay 2 | x[n − 2] | yes | yes | yes | no | yes | h = δ[n − 2], Σ\|h\| = 1 |
| squarer | x[n]² | no | yes | yes | yes | yes | (x₁ + x₂)² ≠ x₁² + x₂²; \|x\| ≤ M gives \|y\| ≤ M² |
| accumulator | Σ_{m≤n} x[m] | yes | yes | yes | no | **no** | a step gives a ramp; h = u, Σ\|h\| = ∞ |
| gain n x[n] | n x[n] | yes | **no** | yes | yes | no | the gain depends on n; a constant input gives a ramp |
| bias correction | x[n] − 0.05 | **no** | yes | yes | yes | yes | affine: S(0) = −0.05 ≠ 0 |
| saturation ±1 | clip(x, −1, 1) | no | yes | yes | yes | yes | linear only while \|x\| ≤ 1 |
| moving average 4 | mean of x[n−3..n] | yes | yes | yes | no | yes | Σ\|h\| = 1 |
| low-pass a = 0.9 | 0.9 y[n−1] + 0.1 x[n] | yes | yes | yes | no | yes | h = 0.1·0.9ⁿ, Σ\|h\| = 1 |
| recursion a = 1.05 | 1.05 y[n−1] + x[n] | yes | yes | yes | no | **no** | h = 1.05ⁿ grows |
| centred smoother | (x[n−1] + x[n] + x[n+1])/3 | yes | yes | **no** | no | yes | uses one future sample |
| decimator | x[2n] | yes | **no** | **no** | no | yes | shifting x by 1 does not shift y by 1; y[n] reads x[2n], ahead of n |

Two entries surprise most people. The **bias correction** subtracts a constant, which feels harmless, but it is *affine*, not linear: S(x₁ + x₂) = x₁ + x₂ − b while S(x₁) + S(x₂) = x₁ + x₂ − 2b. (Its *increments* are linear, which is why the bias can be handled by superposition around an operating point, extension 1.) The **decimator** is linear but neither time-invariant nor causal: y[n] = x[2n] reads the input at index 2n, which for n ≥ 1 is in the future of n, and a one-sample delay of the input becomes a half-sample delay of the output, which does not exist. Multirate blocks therefore need their own theory (`P09.06`, `P09.17`).

### 2.6 Testing a property numerically: what a test can and cannot prove

Each property is a statement "for all inputs". A program can only try some inputs, so each tester is built to *falsify*:

- `is_linear` draws random Gaussian inputs of standard deviation `amplitude` and random weights, and measures the superposition error max|S(ax₁ + bx₂) − aS(x₁) − bS(x₂)| / (1 + max|aS(x₁) + bS(x₂)|). One trial above 10⁻⁹ proves non-linearity.
- `is_time_invariant` applies a random input that is zero on the second half of the record, delays it by k ≤ N/4 samples and compares S(shifted x)[k:] with S(x)[:N−k]. The zero tail keeps the record end from faking a failure for a block with a short look-ahead.
- `is_causal` changes the input from index k onward and checks that y[:k] did not move, for several k.
- `is_memoryless` changes one sample x[k] and checks that only y[k] moved.
- `is_bibo_stable` runs four bounded inputs (|x| ≤ 1: a step, ±1 alternating, random signs, a slow sine) for 4n samples and declares instability if the peak output over the whole record is more than 1.5 times the peak over the first n samples, or not finite.

Every one of these can be fooled, and the ways it can be fooled are the engineering content. A saturation tested with signals far below its limit is indistinguishable from a wire. A stable low-pass whose time constant is longer than the record is still rising at the end and looks unstable; an unstable mode that none of the test inputs excites is missed. A time-varying gain scheduled to change after 60 s passes any test on a 10 s record. Section 4 measures two of these failures.

### 2.7 Where the model stops being valid

Real signal paths run on hardware with fixed-point arithmetic, so even a "linear" filter is linear only to its quantisation step and until it overflows: every digital block is a saturation at some amplitude. Initial rest is a convention; a filter started with a non-zero state is affine until that state decays. A finite record cannot settle questions about infinite time, which is exactly what BIBO stability is: the numerical test is a screening tool, and for LTI blocks the decisive check is the analytic one (poles inside the unit circle in `P09.07` and `P09.17`, in the left half plane in `P09.04`).

## 3. Build the computational model

Open `starter/main.py`. The zoo of gyro-path blocks, `shift_signal`, `bounded_test_inputs` and `main` are given.

### Step 1: signals

`unit_step` is `np.where(t >= 0, 1.0, 0.0)`; `ramp`, `rect_pulse` (as `u(t) − u(t − w)`) and `unit_impulse` follow. `energy(x, dt)` is `np.sum(np.abs(x)**2) * dt`; `average_power(x)` is the mean of |x|². `even_part(x)` is `(x + x[::-1]) / 2` for a grid symmetric about zero. The time operations return functions:

```python
def time_shift(f, t0):
    return lambda t: f(np.asarray(t, dtype=float) - t0)

def time_scale(f, a):
    return lambda t: f(a * np.asarray(t, dtype=float))
```

### Step 2: the property testers

```python
def superposition_error(system, x1, x2, a, b):
    lhs = system(a * x1 + b * x2)
    rhs = a * system(x1) + b * system(x2)
    return np.max(np.abs(lhs - rhs)) / (1.0 + np.max(np.abs(rhs)))
```

`is_linear` loops over seeded trials and returns `False` at the first error above the tolerance. `is_time_invariant`, `is_causal` and `is_memoryless` follow the recipes of 2.6; use `numpy.random.default_rng(seed)` so a failure can be reproduced. `is_bibo_stable` compares the peak over the first n samples with the peak over 4n samples inside `np.errstate(over="ignore")`, because an unstable recursion overflows to `inf` (that is itself a verdict). `impulse_response(system, n)` is `system(unit_impulse(n))`, `abs_sum` is Σ|h|, and `classify` returns the five verdicts as a dictionary.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-01-signals-and-systems`. The tests:

1. **Signals**: the step takes only 0 and 1; the ramp equals the running sum of the step within 2 dt; the pulse area is its width; the impulse.
2. **Energy**: e^{−at} for a = 0.5, 2, 10 converges to 1/(2a) with relative error below 1.2 a·dt at dt = 10⁻² and 10⁻⁴.
3. **Power**: A²/2 for three amplitudes, frequencies and phases; the energy of a sinusoid doubles when the record doubles.
4. **Even and odd**: exact split, symmetry, antisymmetry, orthogonality, Pythagoras for energy; cos is even.
5. **Time operations**: x(2t − 3) by both correct orders, not by the wrong one; double reversal; a positive t₀ delays.
6. **The analytic verdicts**: `classify` reproduces the table of 2.5 for all eleven blocks (eleven parametrised tests).
7. **Saturation**: superposition error below 10⁻¹⁵ inside the limits and above 0.05 at five times the amplitude; `is_linear` says yes at amplitude 0.01 and no at 10.
8. **LTI checks**: zero maps to zero for linear blocks and to −0.05 for the bias correction; h = b aⁿ for three recursions with Σ|h| from the geometric sum; four causal LTI blocks equal `np.convolve(x, h)`; the BIBO bound and its worst-case input; the accumulator's Σ|h| = N and the 36 deg bias drift; the record-length false alarm.
9. **main** writes the figure.

### Compare

There is no library routine that decides "is this function linear?", and there cannot be one: the properties quantify over all inputs. The independent references are the hand derivations of 2.5 and, for the LTI blocks, the identity y = h * x checked against `numpy.convolve`. The numerical testers agree with the hand table on all eleven blocks with the default settings (amplitude 10, records of 256 samples, 2000 and 8000 for stability). The agreement is a property of those settings as much as of the code, which is what part 4 shows.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Classify the zoo.** Run `main()` and compare its property table with your hand table. Predict Σ|h| over the first 2000 samples for each LTI block (1 for the delay, moving average and low-pass, 2000 for the accumulator, astronomically large for a = 1.05). Note what Σ|h| reports for the centred smoother and explain it.
2. **Amplitude sweep.** Feed the ±1 saturation two seeded Gaussian signals scaled to amplitudes 0.1, 0.3, 0.5, 1 and 3 with weights 0.6 and 0.3, and record the superposition error. Predict the amplitude at which the error first becomes non-zero.
3. **Record length.** For the low-pass y = a y[n−1] + (1 − a)x[n] with a = 0.9, 0.99, 0.999 (time constants 10, 100, 1000 samples), run `is_bibo_stable` with n = 50, 500, 5000. Predict which verdicts are wrong.
4. **Bias drift.** Integrate a constant 0.01 deg/s at 100 Hz for one hour with the accumulator (dt = 0.01 s). Predict the attitude error.

Raw data format: CSV tables `system, linear, time_invariant, causal, memoryless, bibo_stable, sum_abs_h` (1); `amplitude, superposition_error` (2); `a, time_constant_samples, n, verdict` (3), each with a header line naming the seed and record length.

What you should find. (1) The code reproduces every hand verdict; Σ|h| is 1 for the delay, moving average and low-pass, 2000 for the accumulator, 4.8 × 10⁴³ for a = 1.05, and 0.6667 for the centred smoother: the impulse response measured from n = 0 misses the h[−1] = ⅓ that the smoother applies *before* the impulse arrives, a reminder that an impulse response recorded from rest cannot see a non-causal part. (2) The superposition error is exactly 0 at amplitude 0.1 (every sample of ax₁ + bx₂ stays inside ±1), 2.4 × 10⁻² at 0.3, 0.17 at 0.5, 0.43 at 1 and 0.68 at 3: the saturation is linear for small signals and the test only finds out when the signals reach the limit. (3) For a = 0.9 all three verdicts are right; for a = 0.99 the n = 50 verdict is wrong; for a = 0.999 the n = 50 and n = 500 verdicts are wrong. A stable filter is declared unstable whenever n is shorter than a few time constants, because its step response is still rising. (4) 36.00 deg after one hour: the gyro bias, a bounded input, produces an output that grows without bound, which is the accumulator's instability in physical form and the reason inertial systems need an external attitude reference.

## 5. Validate and improve

| quantity | hand (§2) | code | independent check | difference | cause |
|---|---|---|---|---|---|
| energy of e^{−2t}u(t) | 1/(2a) = 0.25 | 0.25005 (dt = 10⁻⁴) | — | 5 × 10⁻⁵ | left Riemann sum, ≈ a dt relative |
| power of 2 cos(2π·5t) | A²/2 = 2 | 2.00000 | — | < 10⁻⁵ | whole periods |
| doublet energy (deg²·s) | 2 | 2.0000 | odd part 1.9998 | 1.5 × 10⁻⁴ in the even part | grid points at the jumps |
| property verdicts, 11 blocks | table 2.5 | `classify` | hand table | none | — |
| Σ\|h\| of the low-pass | 1 | 1 | geometric sum | < 10⁻¹² | — |
| LTI output | y = h * x | the block itself | `np.convolve` | < 10⁻¹⁰ | rounding |
| bias drift, 1 h | 0.01 × 3600 = 36 deg | 36.00 deg | — | < 10⁻⁹ | — |
| stability verdict, a = 0.999, n = 50 | stable | "unstable" | Σ\|h\| = 1 | wrong verdict | record shorter than the time constant |

The numbers agree to the precision the theory claims, and the one disagreement in the table is the point of the lesson: the numerical stability test fails on a short record. The improvement for LTI blocks is to stop testing and compute: the impulse response of a recursion is a geometric sequence and Σ|h| is known in closed form, and from `P09.04` onwards the decision is made from the pole locations. For nonlinear blocks no such shortcut exists, and the improvement is better test design: amplitudes that span the whole operating range, records much longer than the slowest time constant, and inputs chosen to excite every mode (the random-sign input catches a marginally stable resonator that the step misses, Problem 6).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand verdict table with a reason per entry; the BIBO proof; `signals_and_systems.png`; the experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) that says which blocks of the gyro path may be analysed with transfer functions, which may run in real time, and what each numerical verdict proves.

### Adversarial self-checks

- **What if the input amplitude doubles?** For every linear block the output doubles; for the squarer it quadruples; the saturation clips. If your design assumes linearity, state the amplitude below which it holds (for the ±1 saturation, the superposition error is exactly zero at amplitude 0.1 and already 2.4 × 10⁻² at 0.3).
- **What are the limiting cases?** a → 0 turns the low-pass into a gain of 1 − a with no memory; a → 1 turns it into a scaled accumulator, and the record needed to see that it is stable grows without bound; |a| > 1 is unstable at any record length.
- **Which assumption fails first in the real world?** Initial rest and exact arithmetic: a fixed-point filter overflows (a saturation) and has a quantisation step (a nonlinearity), so every digital block is linear only over a finite amplitude range.
- **What measurement could falsify a verdict?** One input pair that breaks superposition, one shift that is not reproduced, one future sample that changes the past output: a single counterexample is decisive. No number of passing trials is.
- **Does the result survive a change of test signal or record length?** Not always, and part 4 measures it: the stability verdict flips with record length, the linearity verdict with amplitude. The hand verdicts do not change.
- **Which decision changes if the uncertainty doubles?** If the bias estimate is uncertain by twice as much, the attitude error after an hour without an external reference doubles too (it grows linearly with time and with the bias); the decision is how often the attitude must be reset by an accelerometer or a star tracker.

## Check yourself

1. Give the energy of e^{−3t}u(t) and the average power of 5 cos(100t + 1). Which is an energy signal and which a power signal?
2. Write the even and odd parts of x(t) = e^{−t}u(t). Why is E(x) = E(x_e) + E(x_o)?
3. To build x(3t − 6) from x(t), you scale by 3 first. By how much must you then shift, and in which direction?
4. Is y[n] = x[n] + 2 linear? Time-invariant? Give the one-line reason for each.
5. Why is the decimator y[n] = x[2n] not causal, although it never uses "future" samples of the recorded file?
6. State the BIBO stability condition for an LTI system in terms of h, and show that the accumulator fails it.
7. A numerical test declares a low-pass with time constant 1000 samples unstable on a 50-sample record. Is the filter unstable? What would you change in the test?
8. A gyro bias of 0.005 deg/s is integrated for 20 minutes. What attitude error results, and which system property explains why it keeps growing?
