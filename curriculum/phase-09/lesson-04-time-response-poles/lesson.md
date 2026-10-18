# P09.04  Time Response & Poles

> Where the poles sit decides how a system moves: the real part sets how fast it settles, the imaginary part how fast it rings, and their ratio how much it overshoots. Zeros and extra poles decide how far you can trust that sentence.

| | |
|---|---|
| **Learning outcomes** | 1. Measure rise time, peak time, overshoot, undershoot, settling time and steady-state value from a sampled step response, with interpolation between samples. 2. Derive the second-order step response for every damping ratio and the formulas for overshoot, peak time and the envelope settling time, and invert them to place poles from a specification. 3. Implement and test the step response of any transfer function with simple poles by partial fractions (residues), checked against an RK4 simulation and `scipy.signal.step`. 4. Decide, with numbers, when a higher-order loop is "second order enough": the dominance ratio, the residues, the effect of a third pole, and the effect of a left- or right-half-plane zero. |
| **Prerequisites** | `P09.03` (transfer functions, feedback, the pitch loop), `P02.08` (partial fractions, residues, overshoot and peak time of an actuator), `P02.02` (second-order ODEs, damping ratio), `P02.07` (complex numbers), `P03.05` (RK4). Skill check: you can find the residue of a simple pole by the cover-up rule and write e^{(σ + jω)t} in terms of cos and sin. |
| **Engineering question** | An elevator position servo, a DC motor and gearbox K/(s(τs + 1)) with τ = 50 ms under proportional feedback, must reach a commanded deflection with at most 10 % overshoot. What gain does that require, what are its rise, peak and settling times, and can any proportional gain give a 2 % settling time of 0.3 s or less? Then: when the same formulas are applied to the five-pole pitch-attitude loop of `P09.03`, to a servo with a feedback filter, or to a loop with a zero, how wrong are they? |
| **Assumptions and limits** | Linear, time-invariant models from rest with a unit step command; the servo is a DC motor whose electrical time constant is neglected, with no friction, backlash, saturation or rate limit; the partial-fraction code assumes simple (non-repeated) poles; measured metrics depend on sampling and on knowing the final value; the classic formulas are exact only for the standard second-order system with no zeros. |
| **Analytical work** | The second-order step response for ζ < 1, ζ = 1 and ζ > 1 by partial fractions; peak time π/ω_d and overshoot exp(−πζ/√(1 − ζ²)); the envelope bound −ln(0.02√(1 − ζ²))/(ζω_n) and the 4/(ζω_n) rule; ζ from a required overshoot; the closed-loop servo K/(τs² + s + K) with ω_n = √(K/τ), ζ = 1/(2√(Kτ)) and Re s = −1/(2τ); residues of G(s)/s; the step response with a zero, y = y₂ − (1/z) dy₂/dt, and the initial slope −ω_n²/z of a right-half-plane zero. |
| **Computational deliverable** | `starter/main.py`: `step_metrics`, `second_order_params`, `second_order_step`, `overshoot_formula`, `peak_time_formula`, `settling_time_envelope`, `rise_time_exact`, `zeta_from_overshoot`, `step_residues`, `residue_step_response`, `dominant_poles`, `zero_step_closed_form`, `servo_closed_loop`, `servo_gain_for_overshoot` (`tf`, `tf_to_ss`, `simulate_tf`, `simulate_step`, `second_order_tf`, `with_third_pole`, `with_zero`, `pitch_closed_loop` and `main` are given). Tests prove: measured overshoot, peak time and rise time match the formulas for ζ = 0.2, 0.45, 0.7, the settling time lies below the envelope bound and within 25 % of 4/(ζω_n); first-order rise and settling are τ ln 9 and τ ln 50; metrics scale as 1/ω_n and ignore amplitude and sign; the closed form matches `scipy.signal.step` for ζ = 0, 0.3, 1, 2.5 and is continuous through ζ = 1; ζ from overshoot inverts the formula; partial fractions match RK4 and SciPy on the pitch loop and on random systems; the residue at 0 is G(0) and the residues sum to y(0) = 0; repeated poles are refused; the servo gain gives the target overshoot with Re s = −1/(2τ) for every gain; a far third pole leaves the overshoot unchanged and a near one removes it; LHP zeros nearer the origin add overshoot, a RHP zero gives undershoot, a far zero changes nothing. |
| **Practical activity** | **computational experiment**: design the servo gain for 5, 10 and 20 % overshoot and compare measured and predicted metrics; sweep the gain for the shortest settling time; test the second-order formulas on the pitch loop of `P09.03`; sweep a third pole and a zero across the s-plane. |
| **Equipment** | None. A computer with the environment from `P00.01`. (The same servo plant is closed on real hardware in `P09.21`.) |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `time_response.png`, the derivations of the second-order formulas and of the servo's fixed real part, the design table (target, K, ζ, ω_n, measured and predicted metrics), the third-pole and zero tables, and a short note recommending a servo gain with its metrics and its sensitivity to τ. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the overshoot and peak-time formulas and the fixed real part are derived; all public tests pass; the measured and predicted metrics are compared in a table with the discrepancies explained; the note states why 4/(ζω_n) is misleading for this servo. |
| **References** | `@mit-ocw-16-06-principles-of-automatic-control` Lecture 5 "Dynamic response of closed-loop systems" (PDF "Dynamic Response"), Lecture 6 "Time-domain specifications" (the step response of H(s); desired pole locations) and Lecture 7 "Effect of zeros" (PDF "Effects of Zeros on Step Response", including "Non-Minimum Phase Zeros") (required-free): rise time, peak time, overshoot and settling time against ζ and ω_n, the desired-pole region, dominant poles, and why a left-half-plane zero adds overshoot and a right-half-plane zero gives undershoot. `@mit-ocw-6-003-signals-and-systems` Lecture 3 "Feedback, Poles, and Fundamental Modes" (Geometric Growth: Poles; Second-Order Systems; Partial Fractions; Complex Poles), Lecture 4 "Continuous-Time (CT) Systems" (CT Feedback; Convergent and Divergent Poles) and Lecture 6 "Laplace Transform" (Initial Value Theorem; Final Value Theorem) (required-free): partial fractions as a sum of modes, and the initial and final values read from G(s). |
| **Extensions** | ★ Extend `step_residues` to repeated poles with the derivative formula and check ζ = 1. ★ Add a velocity (tachometer) feedback K_v s to the servo and show that the real part −1/(2τ) can now be moved; find K and K_v for 10 % overshoot and 0.15 s settling. ★ Fit ζ and ω_n to a SYNTHETIC noisy step response by least squares and report their uncertainty (the method `P09.21` uses on measured data). |
| **Time** | unmeasured |

## 1. The engineering problem

A control surface is moved by a position servo: a small DC motor with a gearbox drives the surface, a potentiometer measures its angle, and an amplifier drives the motor in proportion to the error between commanded and measured angle. For small signals the motor and gearbox from voltage to angle behave as K_m/(s(τs + 1)), an integrator (the angle is the integral of the speed) behind a first-order lag (the speed responds to voltage with the mechanical time constant τ, here 50 ms). With a proportional gain the loop gain is K = K_p K_m and the closed loop is second order. The specification is written in the time domain: at most 10 % overshoot (a surface that overshoots loads the structure and excites the airframe) and a 2 % settling time as short as possible.

`P02.08` measured the overshoot and settling time of an actuator whose poles were given. This lesson turns that around: from a specification, place the poles; from the poles, predict the response; from a sampled response, measure the metrics; and compare all three. Then it asks the question every real loop raises: the formulas are exact for the standard second-order system with no zeros, but the pitch-attitude loop of `P09.03` has five poles and two zeros, a servo with a filtered potentiometer has three poles and a zero, and an aircraft's altitude response to the elevator has a zero in the right half plane. When are the formulas good enough?

Success criteria, written first: the measured overshoot and peak time of the designed servo agree with the formulas to the sampling resolution; the partial-fraction step response of the five-pole pitch loop matches an RK4 simulation to 10⁻⁸ and `scipy.signal.step` to 10⁻⁹; and the lesson states, with numbers, the conditions under which the second-order formulas mislead.

## 2. Understand and derive

### 2.1 Measuring a step response

Given samples y(t_k) of the response to a unit step applied at t = 0, with final value y_∞:

- **rise time** t_r: from 10 % to 90 % of y_∞, each crossing found by linear interpolation between samples;
- **peak** y_p and **peak time** t_p: the largest sample; **percentage overshoot** M_p = 100(y_p − y_∞)/y_∞ (zero if never exceeded);
- **undershoot**: how far y goes below zero, in percent of y_∞ (a signature of a right-half-plane zero, 2.6);
- **settling time** t_s: the time after which |y − y_∞| ≤ 0.02|y_∞| for good, found by interpolating the last exit from the band.

The final value must be known: from the model (G(0) by the final-value theorem) or from a long record. A biased estimate of y_∞ moves every metric.

### 2.2 The standard second-order system

G(s) = ω_n²/(s² + 2ζω_n s + ω_n²) has poles at **s = −ζω_n ± jω_n√(1 − ζ²) = −σ ± jω_d** for ζ < 1. Partial fractions of G(s)/s (`P02.08`) give the step response:

- ζ < 1: **y(t) = 1 − e^{−σt} sin(ω_d t + φ)/√(1 − ζ²)**, φ = arccos ζ;
- ζ = 1: y = 1 − e^{−ω_n t}(1 + ω_n t);
- ζ > 1: y = 1 + (p₂e^{p₁t} − p₁e^{p₂t})/(p₁ − p₂) with real poles p₁,₂ = −ζω_n ± ω_n√(ζ² − 1);
- ζ = 0: y = 1 − cos(ω_n t), oscillating between 0 and 2 for ever.

Setting ẏ = 0 (16.06 Lecture 6) gives ω_n e^{−σt} sin(ω_d t)/√(1 − ζ²) = 0, so the first peak is at **t_p = π/ω_d** and

**M_p = 100 exp(−πζ/√(1 − ζ²))**,   inverted:   **ζ = −ln M/√(π² + ln² M)**,  M = M_p/100.   (1)

Overshoot depends on ζ alone, which is the angle of the poles from the negative real axis (cos⁻¹ ζ). The decay envelope 1 ± e^{−σt}/√(1 − ζ²) enters the 2 % band at

**t_s ≤ −ln(0.02√(1 − ζ²))/(ζω_n) ≈ 4/(ζω_n)**,   (2)

an upper bound (the response itself can be inside the band before the envelope is). The rule 4/(ζω_n) drops the √(1 − ζ²) and rounds ln 50 = 3.91 to 4. Settling depends on the real part σ = ζω_n. The 10–90 % rise time has no closed form; `rise_time_exact` finds it by bisection on y(t).

### 2.3 The servo under proportional feedback

Unity feedback around K/(s(τs + 1)) gives

**T(s) = K/(τs² + s + K)**,  **ω_n = √(K/τ)**,  **ζ = 1/(2√(Kτ))**,  **ζω_n = 1/(2τ)**.   (3)

The last equality is the key fact of this lesson. For every gain above K = 1/(4τ) (where ζ < 1) the closed-loop poles are **−1/(2τ) ± jω_d**: raising K raises ω_d and the overshoot but **cannot move the real part**. With τ = 50 ms the real part is −10 s⁻¹ and the 4/(ζω_n) rule predicts 0.40 s for every gain. From (1), 10 % overshoot needs ζ = 0.5912, hence K = 1/(4ζ²τ) = **14.31 s⁻¹**, ω_n = 16.92 rad/s.

The measured 2 % settling time is 0.350 s, not 0.40 s: the envelope is an upper bound, and the true settling time depends on where the oscillation crosses the band. Sweeping K from 5 to 60 s⁻¹, the shortest settling is **0.282 s at K = 8.2 s⁻¹**, with 1.97 % overshoot: the response settles on its first peak because the peak only just exceeds the band. Settling ≤ 0.3 s is reached only for K between 7.9 and 8.2 s⁻¹, a window so narrow that a 20 % change in τ destroys it (Problem 5). A proportional servo cannot be made to settle reliably faster; velocity feedback (extension 2) or a lead compensator (`P09.12`) is needed to move the real part.

### 2.4 Any transfer function with simple poles: residues

For a proper G = N_G/D_G, the step response transform is Y(s) = G(s)/s = N(s)/D(s) with D(s) = s D_G(s). If the roots p_k of D are simple,

**y(t) = Σ_k r_k e^{p_k t},  r_k = N(p_k)/D′(p_k)**   (4)

(`P02.08`; 6.003 Lecture 3, "Partial Fractions"). The residue at p = 0 is G(0), the steady state. For a strictly proper G, y(0) = Σ r_k = 0. Each other residue measures **how strongly the step excites that mode**: a pole with a tiny residue hardly shows in the response wherever it sits. (4) is exact; the test compares it with RK4 and with `scipy.signal.step`.

### 2.5 Dominant poles and the third pole

The **dominant poles** are those closest to the imaginary axis, and the usual rule says a higher-order system behaves like its dominant pair if the next pole is at least 5 times further left (**dominance ratio** ≥ 5). Two things can break the rule.

*A third pole near the pair.* Take ζ = 0.5, ω_n = 1 (overshoot 16.3 %, 2 % settling 8.08 s) and add a real pole at −p₃ with unit DC gain. At p₃ = 10ζω_n the overshoot is 15.9 % and the settling 8.26 s: second order is enough. At p₃ = 2ζω_n: 8.2 % and 6.64 s. At p₃ = ζω_n the real pole is as slow as the pair: no overshoot at all, settling 9.14 s, and at 0.5ζω_n it dominates outright (16.5 s). A slow third pole acts like a filter on the command.

*Zeros and residues.* The pitch-attitude loop of `P09.03` has poles −100.9, −17.24, −2.00 ± 7.36j and −1.068, and zeros −1.5 and −100. By real part the "dominant" pole is the real one at −1.068, with a dominance ratio of only 1.87 to the pair. Its residue (0.32) is comparable with the pair's (0.37 each), because the zero at −1.5 sits next to it and partly cancels it. Neither the real pole alone nor the pair alone describes the response: the pair alone (ζ = 0.262, ω_n = 7.62) predicts 42.6 % overshoot and 0.427 s peak time, while the loop overshoots 7.8 % at 0.510 s and settles in 2.74 s on the slow tail of the −1.068 pole. When poles and zeros are interleaved, use (4), not the formulas.

### 2.6 The effect of a zero

Add a zero at s = z (16.06 Lecture 7) to the standard system with unit DC gain: G_z(s) = ω_n²(1 − s/z)/(s² + 2ζω_n s + ω_n²) = G(s) − (s/z)G(s). Multiplying by s differentiates, so

**y_z(t) = y(t) − (1/z) ẏ(t)**,   ẏ = (ω_n/√(1 − ζ²)) e^{−σt} sin(ω_d t),   (5)

the closed form `zero_step_closed_form`. A **left-half-plane zero** (z < 0) adds a positive multiple of ẏ: the response is faster and overshoots more, increasingly as the zero approaches the origin (ζ = 0.5, ω_n = 1: overshoot 16.7 % at z = −5, 19.1 % at −2, 29.8 % at −1, 69.9 % at −0.5). A **right-half-plane zero** (z > 0) subtracts ẏ: the initial slope is −ω_n²/z < 0 and the response **starts in the wrong direction** (undershoot 9.1 % at z = 2, 75 % at z = 0.5). Such a system is called non-minimum phase. An aircraft's altitude response to elevator is the classic example: pulling the elevator up first pushes the tail down, and the lift on the tailplane makes the aircraft sink briefly before it climbs.

### 2.7 Where the model stops being valid

The servo model neglects the motor's electrical time constant (a third pole, typically a few milliseconds), friction and backlash (which make the final value depend on the direction of approach), amplifier saturation (the large initial error of a big step demands more voltage than the supply has, and the response becomes rate-limited) and the potentiometer's filter (Problem 4). The metrics also depend on the measurement: a 2 % band on a 1 deg step is 0.02 deg, below the resolution of many hobby-servo potentiometers.

## 3. Build the computational model

Open `starter/main.py`. `tf`, `tf_to_ss`, `simulate_tf`, `simulate_step`, `second_order_tf`, `with_third_pole`, `with_zero`, `pitch_closed_loop` and `main` are given (the first four from `P09.03`).

### Step 1: the metrics

```python
def step_metrics(t, y, y_final=None, band=0.02):
    yf = y[-1] if y_final is None else y_final
    rise = crossing(t, y, 0.9 * yf) - crossing(t, y, 0.1 * yf)     # interpolated upward crossings
    i_pk = np.argmax(y)
    out = np.abs(y - yf) > band * abs(yf)
    j = np.nonzero(out)[0][-1]                                      # last sample outside the band
    settle = t[j] + (|y_j - yf| - tol) * (t[j+1] - t[j]) / (|y_j - yf| - |y_{j+1} - yf|)
    ...
```

Handle negative steps by working with sgn(y_∞)·y, return `nan` for the settling time if the response is outside the band at the end, and report the undershoot as −100 min(y)/y_∞ when positive.

### Step 2: formulas and the closed form

`second_order_params(den)` reads ζ and ω_n from a quadratic; `second_order_step` implements the four cases of 2.2; `overshoot_formula`, `peak_time_formula`, `settling_time_envelope` and `zeta_from_overshoot` are (1) and (2); `rise_time_exact` bisects y(t) = 0.1 and 0.9.

### Step 3: residues

`step_residues(G)` forms D(s) = s·den(s), finds its roots, refuses repeated roots, and returns N(p)/D′(p) with `np.polyval` and `np.polyder`; `residue_step_response` sums r_k e^{p_k t} and keeps the real part.

### Step 4: dominance, zeros, servo

`dominant_poles` sorts the poles by real part, returns the leading pole or conjugate pair and the ratio of the next real part to its real part. `zero_step_closed_form` is (5). `servo_closed_loop(K, τ)` is (3) and `servo_gain_for_overshoot` is K = 1/(4ζ²τ) with ζ from (1).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-04-time-response-poles`. The tests:

1. **Metrics against formulas** for ζ = 0.2, 0.45, 0.7 on a 0.5 ms grid: overshoot within 10⁻⁴ %, peak time within 0.6 ms, rise time against `rise_time_exact` within 1 μs, settling below the envelope bound and within 25 % of 4/(ζω_n), no undershoot.
2. **First order**: rise τ ln 9 and settling τ ln 50 to 10⁻⁵ relative.
3. **Scaling**: metrics of ω_n = 5 are those of ω_n = 1 divided by 5; a step of −3 has the same percentages.
4. **Closed form** against `scipy.signal.step` for ζ = 0, 0.3, 1, 2.5; continuity at ζ = 1 ± 10⁻⁶; the undamped case 1 − cos t; the formula inversions.
5. **Partial fractions** against RK4 (10⁻⁸) and `scipy.signal.step` (10⁻⁹) on the pitch loop and on three random stable systems; the residue at zero equals G(0) and the residues sum to zero; a double pole raises.
6. **Servo design** for 5, 10, 20 %: the simulated overshoot within 0.02 % of the target and Re s = −1/(2τ) exactly.
7. **Third pole**: at 20ζω_n the dominance ratio is 20 and the overshoot within 5 % of second order; at 0.5ζω_n the dominant pole is real at −0.25 and the overshoot vanishes.
8. **Zeros**: (5) equals the residue response; overshoot grows as an LHP zero approaches the origin; an RHP zero matches `scipy.signal.step` and undershoots by more than 5 %; a zero at −10⁶ changes nothing.
9. **main** writes the figure.

### Compare

`scipy.signal.step` simulates the system with the exact matrix exponential and is the oracle; `scipy.signal.residue` would give the same partial fractions as `step_residues`. Our RK4 and residue routes agree with it to 10⁻⁹ on every test case. What SciPy does not provide is the metric extraction: the definitions of 2.1 are conventions, and different tools (and textbooks) use 0–100 % instead of 10–90 % rise time, 5 % instead of 2 % settling, or the first rather than the last band crossing. State the convention with every number.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction first.

1. **Servo design.** For targets of 5, 10 and 20 % overshoot with τ = 50 ms, compute K, ζ, ω_n and ζω_n, and predict the overshoot, peak time, rise time, envelope settling and 4/(ζω_n). Then simulate (RK4, 0.5 ms grid) and measure.
2. **Settling sweep.** Sweep K from 5 to 60 s⁻¹ and record the measured 2 % settling time. Predict the shape of the curve from (2) before you look.
3. **Second order or not.** For the pitch loop of `P09.03`, list the poles with their residues; predict the overshoot from the complex pair alone and from the slow real pole alone; compare with the measured response.
4. **Third pole and zero.** For ζ = 0.5, ω_n = 1, add a third pole at p₃/(ζω_n) = 0.5, 1, 2, 5, 10 and, separately, a zero at z = −0.5, −1, −2, −5, +2, +0.5. Predict the trend of overshoot and settling in each case.

Raw data format: CSV tables `target_pct, K, zeta, wn, OS_meas, OS_pred, tp_meas, tp_pred, tr_meas, tr_pred, ts_meas, ts_env, ts_4` (1); `K, ts_meas, OS_pred` (2); `pole_re, pole_im, residue_abs` (3); `case, value, OS_pct, undershoot_pct, tp_s, ts_s` (4).

What you should find. (1)

| target | K (s⁻¹) | ζ | ω_n (rad/s) | overshoot meas./pred. | t_p meas./pred. (s) | t_r (s) | t_s meas. / envelope / 4/(ζω_n) (s) |
|---|---|---|---|---|---|---|---|
| 5 % | 10.50 | 0.6901 | 14.49 | 5.000 / 5.000 | 0.2995 / 0.2996 | 0.1447 | 0.414 / 0.424 / 0.400 |
| 10 % | 14.31 | 0.5912 | 16.92 | 10.000 / 10.000 | 0.2305 / 0.2303 | 0.1084 | 0.350 / 0.413 / 0.400 |
| 20 % | 24.05 | 0.4559 | 21.93 | 20.000 / 20.000 | 0.1610 / 0.1609 | 0.0710 | 0.380 / 0.403 / 0.400 |

ζω_n = 10.00 s⁻¹ in every row. The overshoot and peak time agree to the grid; the settling time does not follow 4/(ζω_n): it is 0.35–0.41 s and not monotonic in the target. (2) The sweep is a saw-tooth: settling falls as the overshoot grows toward the band, jumps up when an extra half-cycle must decay, and has its minimum, 0.282 s, at K = 8.2 s⁻¹ (1.97 % overshoot); only 7.9 ≤ K ≤ 8.2 s⁻¹ gives 0.3 s or less. (3) The residue table of 2.5: |r| = 0.32 for −1.068, 0.37 for each of −2.00 ± 7.36j, 0.14 for −17.24 and 7.7 × 10⁻⁶ for −100.9; the pair alone predicts 42.6 % overshoot against 7.8 % measured. (4) The third-pole and zero tables of 2.5 and 2.6: overshoot falls monotonically as the third pole approaches the pair, while settling first improves (6.64 s at ratio 2) and then degrades; LHP zeros near the origin multiply the overshoot, RHP zeros produce undershoot that grows as the zero approaches the origin.

## 5. Validate and improve

| quantity | formula (§2) | measured (code) | independent check | difference | cause |
|---|---|---|---|---|---|
| overshoot, 10 % design | 10.000 % | 10.000 % | RK4 vs closed form | < 10⁻³ % | grid |
| peak time, 10 % design | π/ω_d = 0.2303 s | 0.2305 s | — | 0.2 ms | 0.5 ms grid |
| rise time, 10 % design | exact (bisection) 0.1084 s | 0.1084 s | — | < 10⁻⁴ s | interpolation |
| settling, 10 % design | 4/(ζω_n) = 0.400 s | 0.350 s | envelope bound 0.413 s | −12 % | the rule bounds the envelope, not the response |
| Re of servo poles | −1/(2τ) = −10 s⁻¹ | −10 s⁻¹ for every K | `np.roots` | 0 | exact |
| pitch loop, partial fractions vs RK4 | equal | 1.1 × 10⁻¹⁰ | `scipy.signal.step` | — | RK4 error |
| pitch loop overshoot | pair alone: 42.6 % | 7.8 % | (4) | ×5.5 | zero at −1.5 and slow pole |

The formulas are exact where their assumptions hold (overshoot and peak time of the standard system) and misleading where they do not (settling time, any loop with zeros or interleaved poles). The improvement is not a better formula but the habit of computing (4) and measuring the response. For the servo design the improvement is structural: a proportional loop's real part is fixed at −1/(2τ), and a faster, robust settling time needs velocity feedback or a lead compensator, which move it.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), (3) and (5); `time_response.png`; the four experiment tables; and a half-page note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) recommending a proportional gain for the elevator servo, with its predicted and measured metrics, its sensitivity to a ±20 % error in τ, and what you would add to make it settle faster.

### Adversarial self-checks

- **What if the gain doubles?** From (3), ω_n grows by √2, ζ falls by √2 (10 % → 0.418, about 23.6 % overshoot), the real part stays at −10 s⁻¹. Peak time shortens, settling barely changes.
- **What are the limiting cases?** K → 0: two real poles, one at 0 and one at −1/τ; the response is slow and never overshoots. K → ∞: ζ → 0 and the servo rings at ever higher frequency with the same decay rate: no gain gives a fast, well-damped response.
- **Which assumption fails first in the real world?** Saturation of the amplifier: the loop commands a motor speed proportional to the error, K × 10 deg = 143 deg/s at the start of a 10 deg step with K = 14.31 s⁻¹, and the motor's speed limit turns the start of the response into a ramp. Then backlash, which leaves a dead zone around the target.
- **What measurement could falsify the model?** A step test on the servo (`P09.21`): if the overshoot matches but the period does not, τ is wrong; if the response has a slow tail, there is friction or a third pole; if small and large steps differ in shape, it is not linear.
- **Does the result survive a change of method?** The closed form, the residues, RK4 and SciPy agree to 10⁻⁹; the metric definitions do not survive a change of convention (5 % instead of 2 % settling changes every settling number), so state them.
- **Which decision changes if τ is uncertain by ±20 %?** With K = 14.31 s⁻¹, τ = 40 ms gives 6.3 % overshoot and τ = 60 ms gives 13.4 %, violating the specification. Design for the largest plausible τ, or measure it first.

## Check yourself

1. Write the 10–90 % rise time, peak time, percentage overshoot and 2 % settling time in words, and say which one depends only on ζ for the standard second-order system.
2. A second-order system has poles at −3 ± 4j. Give ζ, ω_n, the overshoot and the peak time.
3. What damping ratio gives 5 % overshoot?
4. Show that the servo K/(s(τs + 1)) under unity feedback has closed-loop poles with real part −1/(2τ) for every K > 1/(4τ). What does that imply for its settling time?
5. A system has a pole pair at −2 ± 5j and a real pole at −20. Is the pair dominant? What if the real pole were at −3?
6. Why can a pole with a large dominance ratio still be invisible in the step response, and how do you check?
7. What does a right-half-plane zero do to a step response, and give an aerospace example.
8. Why is 4/(ζω_n) an upper estimate of the 2 % settling time, and by how much was it wrong for the 10 % servo design?
