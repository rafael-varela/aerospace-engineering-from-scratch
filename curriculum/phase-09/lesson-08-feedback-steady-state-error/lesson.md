# P09.08  Feedback & Steady-State Error

> Feedback shrinks an error by the factor 1 + L. To make a constant error vanish you need an integrator in the loop; to follow a ramp you need two. The final-value theorem says how much error is left, but only if the loop is stable.

| | |
|---|---|
| **Learning outcomes** | 1. Write the closed-loop transfer function T = L/(1 + L), the sensitivity S = 1/(1 + L) and the disturbance path P/(1 + CP) of a unity-feedback loop from polynomial transfer functions, and show S + T = 1. 2. Determine the system type and the error constants K_p, K_v, K_a, and compute the steady-state error to a step, a ramp and a parabola. 3. Apply the final-value theorem correctly, including the stability condition under which it holds, and verify every prediction against a long from-scratch time simulation. 4. Design proportional and PI airspeed-hold controllers for a light aircraft, quantify their error to a speed command, a speed ramp and a drag increase, and find the gain at which the proportional loop goes unstable. |
| **Prerequisites** | `P09.03` (transfer functions and block diagrams, feedback reduction), `P09.04` (time response and poles), `P09.05` (frequency response), `P02.08` (Laplace transform, the final-value theorem), `P02.05` (the matrix exponential and the zero-order-hold step), `P02.06` (linearization). Skill check: you can reduce a single feedback loop to L/(1 + L) and evaluate lim_{s→0} of a rational function. |
| **Engineering question** | A light aircraft (1100 kg, 55 m/s, 1080 N of drag at trim) holds airspeed with an auto-throttle that commands thrust through a 0.3 s throttle servo and a 1.0 s engine lag. With a proportional gain of 400 N per m/s, how much airspeed error remains after a 1 m/s speed command and after a 200 N drag increase, how large can the gain be made before the loop oscillates, and what does adding an integrator (K_i = 40 N/m) change for a step, for a 0.5 m/s² speed ramp and for the drag increase? Answer: 8.94 % of the step and 0.455 m/s of speed after the drag increase (against 5.09 m/s with no feedback); the loop goes unstable above 4990 N per m/s; with PI the step and drag errors vanish and the ramp leaves 0.491 m/s. |
| **Assumptions and limits** | Level flight with altitude held by another loop (the phugoid exchange of speed and height is ignored); thrust independent of airspeed and drag linearised about trim, dD/dV = 2D₀/V₀ (constant drag coefficient); first-order throttle servo and engine lags; no throttle saturation, so large commands and ramps are only valid while the thrust stays within its range; airspeed measured perfectly. The numbers are representative of a light single-engine aircraft, not one type's data. The simulation is linear and exact for piecewise-constant inputs; it is not a flight test. |
| **Analytical work** | Linearised along-track equation m v̇ = ΔT − cv + d with c = 2D₀/V₀; the plant G(s) = 1/((τ_a s + 1)(τ_e s + 1)(ms + c)); T, S and the disturbance path by block-diagram algebra; E = S R; the final-value theorem and its condition; e_ss = A/(1 + K_p), A/K_v, A/K_a; the disturbance error d/(c + K) under P control and 0 under PI; the Routh condition for the cubic characteristic polynomial. |
| **Computational deliverable** | `starter/main.py`: `open_loop`, `closed_loop`, `sensitivity`, `disturbance_to_output`, `is_stable`, `system_type`, `error_constants`, `steady_state_error`, `final_value`, `critical_gain` (`poly_add`, `poly_trim`, `tf_eval`, `cruise_plant`, `tf_to_companion`, `expm_series`, `simulate_tf`, `pi_controller`, `p_controller`, `fvt_vs_simulation` and `main` are given). Tests prove: the given simulator matches `scipy.signal.lsim`; L, T, S and the disturbance path equal their definitions at four complex points and S + T = 1 on the imaginary axis; the type counts the integrators (0, 1, 2); the error constants and errors follow from the type; error falls monotonically with gain and is 1 at K → 0; the final-value theorem equals the 600 s simulation for steps and drag increases (exactly) and ramps (within the zero-order hold's A·dt/2); it refuses an unstable loop whose simulation diverges and an undamped oscillator; the critical gain equals Routh's a₂a₁/a₃ − a₀. |
| **Practical activity** | **computational experiment**: a gain sweep of the proportional loop (error, damping, stability limit), an integral-gain sweep of the PI loop (ramp error against damping), and a final-value-theorem against long-simulation comparison for five cases, including the step-size effect on the ramp. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `cruise_control.png`, the derivations of the plant, of E = SR, of the three error formulas and of the disturbance errors, the sweep tables, and a half-page recommendation of a controller with its step, ramp and disturbance errors and its stability margin in gain. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the plant and the error formulas are derived; all public tests pass; the final-value predictions are confirmed by simulation and the one refused case is explained; the recommendation states numbers and the model's limits. |
| **References** | `@mit-ocw-6-003-signals-and-systems` Lecture 16 "Fourier Transform" (supplementary-free: the s-domain view behind the final-value theorem). Feedback, system type, error constants and sensitivity: GAP in the mapped free sources; this lesson is original instruction. `@scipy-userguide` "Signal Processing (scipy.signal)" (supplementary-free: `lsim`, the simulation oracle). Map `references/phase-09.md` not yet available when written; only ids already in `sources.yaml` are cited, with lecture titles as recorded in `references/phase-02.md`. |
| **Extensions** | ★ Add throttle saturation (thrust between −800 N and +1500 N about trim) to the simulation and find the largest speed step the PI loop follows without overshoot growing (integrator windup, the subject of `P09.11`). ★ Add the phugoid (a second-order speed–height exchange) and see how much of the speed error the elevator loop must absorb. ★ Close the loop with a measured airspeed that includes a 0.5 s pitot-static lag in the feedback path, and show that the error is then no longer S R. |
| **Time** | unmeasured |

## 1. The engineering problem

An auto-throttle holds airspeed. The pilot (or the flight-management computer) sets a target; the controller compares it with the measured airspeed and moves the throttle. Two things disturb it: the target changes (a step to a new speed, or a scheduled acceleration, which is a ramp), and the aircraft changes (a gear or flap extension, icing, or a climb, all of which add drag or remove thrust). The question every such loop must answer is how much error remains once the transients have died out. That is the **steady-state error**. It decides whether the auto-throttle can hold an approach speed to a knot, whether it lags a schedule, and whether a drag change makes it sag.

The answer comes from the structure of the loop, not from its details. Proportional feedback divides a constant error by 1 + L(0), so a higher gain gives a smaller error but never zero, and too high a gain makes the loop oscillate. An integrator in the controller makes L(0) infinite and the constant error zero. That is why nearly every practical regulator, from an auto-throttle to a satellite's attitude controller, contains one. This lesson computes those errors three ways: by the error constants, by the final-value theorem, and by simulating the closed loop for ten minutes and reading the error at the end. It also shows where the final-value theorem silently lies.

Success criteria, written first: T, S and the disturbance transfer function agree with their definitions to 10⁻¹²; the final-value prediction agrees with a 600 s simulation to 10⁻⁶ m/s for steps and drag increases, and to the explainable zero-order-hold offset for ramps; the critical proportional gain agrees with the Routh condition to 10⁻⁸.

## 2. Understand and derive

### 2.1 The plant

Along the flight path, with altitude held, m dV/dt = T − D. Linearise about trim (V₀, T₀ = D₀) with v = V − V₀. With drag D = ½ρV²SC_D and C_D constant, dD/dV = 2D₀/V₀ = c, so

**m v̇ = ΔT − c v + d,**   (1)

where d is any extra force along the path (negative for a drag increase or for the weight component −mg sin γ in a climb). For the representative aircraft, c = 2 × 1080/55 = **39.27 N per m/s**, and the speed time constant m/c = **28.0 s**: left alone, the aircraft takes about half a minute to settle at a new speed. The throttle servo (τ_a = 0.3 s) and the engine and propeller (τ_e = 1.0 s) lag the commanded thrust u:

**G(s) = v/u = 1 / ((τ_a s + 1)(τ_e s + 1)(m s + c)),**   (2)

with poles at −3.333, −1.000 and −0.0357 s⁻¹. The DC gain is 1/c: one newton of thrust buys 1/39.27 m/s of speed.

### 2.2 The loop, the error and the sensitivity

With a controller C(s) acting on the error e = r − v and the disturbance d entering at the plant input, the loop gain is L = CG and

**V = T R + G S D,  E = R − V = S R − G S D,**   (3)

with **T = L/(1 + L)** (the closed loop), **S = 1/(1 + L)** (the sensitivity) and the disturbance path **G/(1 + CG) = GS**. By construction **S + T = 1** at every frequency. Where the loop gain is large, S ≈ 1/L is small: commands are followed and disturbances rejected. Where it is small, S ≈ 1: feedback does nothing. In polynomials, L = N/D gives T = N/(D + N) and S = D/(D + N): the closed-loop poles are the roots of D + N.

### 2.3 Steady-state error: system type and the error constants

For a reference R(s) = A/s^{q+1} (q = 0 step, 1 ramp, 2 parabola), E(s) = S(s)R(s). If the closed loop is stable, the **final-value theorem** (`P02.08`) gives

**e_ss = lim_{s→0} s E(s) = lim_{s→0} A / (s^q (1 + L(s))).**   (4)

Write L(s) = K₀ N₀(s)/(s^n D₀(s)) with N₀(0) = D₀(0) = 1. The integer n, the number of free integrators in the loop, is the **system type**. Define the **error constants**

**K_p = lim L(s),  K_v = lim s L(s),  K_a = lim s² L(s)  (s → 0).**   (5)

Then (4) gives

| type | step A/s | ramp A/s² | parabola A/s³ |
|---|---|---|---|
| 0 | A/(1 + K_p) | ∞ | ∞ |
| 1 | 0 | A/K_v | ∞ |
| 2 | 0 | 0 | A/K_a |

For the auto-throttle with C = K_p (proportional), the loop is type 0 with L(0) = K_p/c. With C = K_p + K_i/s (PI), it is type 1 with K_v = K_i/c. A pure integrator in the plant, as in a satellite's attitude θ = ∫∫τ/J (Problem 3), makes the loop type 2 with a plain PD controller.

### 2.4 Disturbances

From (3), a constant drag increase of size F (d = −F/s) leaves the speed error

**e_ss = lim_{s→0} F G S = F/(c + C(0)) for C(0) finite, 0 if C has an integrator.**   (6)

With no feedback, the aircraft loses F/c: **5.09 m/s** for F = 200 N. With K_p = 400 N per m/s it loses 200/(39.27 + 400) = **0.455 m/s**. With PI it loses nothing, because the integrator keeps adding thrust as long as any error remains and settles only when the extra thrust equals F.

### 2.5 When the final-value theorem lies

lim_{s→0} sE(s) is a number whenever the limit exists, but it equals lim_{t→∞} e(t) only if **every pole of sE(s) has a negative real part**. For a proportional gain above the stability limit, (4) still returns 1/(1 + K_p/c), a small positive error, while the real error grows without bound. With K_p = 6000 N per m/s the formula gives 0.65 % and the simulated error is 19.6 m/s after 60 s, and growing. `final_value` therefore checks the poles of sE(s) and raises `ValueError` if any lies in the closed right half plane. An undamped oscillator (poles on the imaginary axis) has no final value either.

### 2.6 How much gain? Routh for the cubic

With C = K, the characteristic polynomial is (τ_a s + 1)(τ_e s + 1)(ms + c) + K = a₃s³ + a₂s² + a₁s + (a₀ + K), with a₃ = 330, a₂ = 1441.78, a₁ = 1151.05, a₀ = 39.27. Routh's array for a cubic (`P09.09` builds it in general) has first column a₃, a₂, (a₂a₁ − a₃(a₀ + K))/a₂, a₀ + K. All positive means stable:

**K < K_crit = a₂a₁/a₃ − a₀ = 4989.7 N per m/s,**   (7)

at which a pair of poles sits on the imaginary axis at ω = √(a₁/a₃) = 1.868 rad/s. At K_crit the type-0 step error would be 0.78 %. It cannot be reached in practice, because the loop is oscillating long before that. The damping of the closed-loop poles falls steadily as K rises (2.7).

### 2.7 The numbers

Proportional control, from `main()`:

| K_p [N per m/s] | L(0) | step error | least damping ratio of the closed-loop poles |
|---|---|---|---|
| 100 | 2.55 | 28.20 % | 1.000 |
| 400 | 10.19 | 8.94 % | 0.722 |
| 1000 | 25.46 | 3.78 % | 0.387 |
| 2000 | 50.93 | 1.93 % | 0.194 |
| 4740 (0.95 K_crit) | 120.70 | 0.82 % | 0.009 |

To get the step error below 2 % the damping ratio must fall below 0.2: proportional control trades error against oscillation. PI control with K_p = 400 N per m/s and K_i = 40 N/m has closed-loop poles at −3.472, −0.385 ± 0.357j and −0.127 s⁻¹, K_v = 1.0185 s⁻¹, zero step and drag errors, and a ramp error of 0.5/1.0185 = **0.491 m/s** for a 0.5 m/s² speed schedule (a 10 m/s change in 20 s). The price is the slow real pole at −0.127 s⁻¹ (an 8 s time constant): a 1 m/s step overshoots by 18.0 % and takes 23.5 s to settle within 2 %, against 3.7 % overshoot and 9.9 s (to its own, wrong, final value of 0.911 m/s) for K_p = 400 alone.

The sensitivity of the PI loop peaks at \|S\| = 1.405 at 0.627 rad/s: disturbances near a 10 s period are *amplified* by 40 %, the price paid for rejecting slow ones (\|S(j0.01)\| = 0.0102). \|T(j5)\| = 0.0080: airspeed-sensor noise at 5 rad/s barely reaches the throttle loop's output.

## 3. Build the computational model

Open `starter/main.py`. The plant, the polynomial helpers, the simulator and `main()` are given.

### Step 1: loop algebra in polynomials

```python
def open_loop(num_c, den_c, num_p, den_p):
    return np.convolve(num_c, num_p), np.convolve(den_c, den_p)

def closed_loop(num_l, den_l):                       # T = N / (D + N)
    return poly_trim(num_l), poly_trim(poly_add(den_l, num_l))
```

`sensitivity` returns (D, D + N); `disturbance_to_output` returns (N_P D_C, D_P D_C + N_P N_C); `is_stable` checks `np.roots(den).real < 0`.

### Step 2: type, constants, errors

`system_type` counts trailing zero coefficients (roots at s = 0) of the denominator minus those of the numerator. `error_constants` strips the common powers of s and evaluates the ratio of the constant terms: that is lim sⁿL; the lower-index constants are ∞ and the higher ones 0. `steady_state_error` applies the table of 2.3 after checking that the closed loop is stable.

### Step 3: the final-value theorem and the critical gain

`final_value(num, den)` multiplies num by s, cancels common factors of s, then checks the roots of the remaining denominator before returning num_s[−1]/den[−1]. `critical_gain` bisects on K between a tiny gain and 10⁹, testing `is_stable(den + K num)`.

### Step 4: the independent check

`simulate_tf` (given) builds the companion state-space form and steps it with the exact zero-order-hold matrices from the augmented matrix exponential (`P02.05`). It equals `scipy.signal.lsim(..., interp=False)` to 10⁻¹⁰. `fvt_vs_simulation` runs a case for 600 s at dt = 0.05 s and returns both the theorem's answer and the simulated error at the end. Running `main()` prints the tables of 2.7 and

```
final-value theorem against a 600 s simulation (dt = 0.05 s):
  P 400, step 1 m/s       : FVT +0.08940 m/s, simulation at 600 s +0.08940 m/s
  P 400, 200 N drag step  : FVT +0.45530 m/s, simulation at 600 s +0.45530 m/s
  PI, step 1 m/s          : FVT +0.00000 m/s, simulation at 600 s +0.00000 m/s
  PI, ramp 0.5 m/s^2      : FVT +0.49091 m/s, simulation at 600 s +0.50341 m/s
  PI, 200 N drag step     : FVT -0.00000 m/s, simulation at 600 s +0.00000 m/s
```

and writes `cruise_control.png` (the step responses of K_p = 400, K_p = 1000 and PI; their response to a 200 N drag increase; \|S\| and \|T\| of the PI loop).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-08-feedback-steady-state-error`. The 19 tests check:

1. **Simulator** against `scipy.signal.lsim` with held input, to 10⁻¹⁰.
2. **Loop algebra**: L, T, S and the disturbance path at four complex points against the definitions; S + T = 1 on the jω axis; \|S\| at low frequency equals the step error; \|S\| → 1 at high frequency.
3. **Type and constants**: types 0, 1, 2 for P, PI and a double integrator; K_p = K/c, K_v = K_i/c, K_a = K/J for a satellite PD loop.
4. **Errors**: step, ramp and parabola errors from the constants; ValueError for an unstable loop; the error is 1 at K → 0, falls monotonically with K, and never reaches 0 for type 0.
5. **Final-value theorem against simulation** for five cases; the ramp case differs by exactly A·dt/2.
6. **The theorem's conditions**: an unstable loop (K_p = 6000) is refused and its simulation diverges; an undamped oscillator is refused.
7. **Routh**: `critical_gain` equals a₂a₁/a₃ − a₀; stable at 0.99 K_crit, unstable at 1.01 K_crit; a second-order plant is never destabilised by gain.
8. **Superposition** and zero input of the closed-loop simulation.

### Compare

`scipy.signal.lsim` is the oracle for the simulator. There is no SciPy routine for "steady-state error"; the independent check is the long simulation itself. One detail matters. For the ramp, the simulation's error exceeds the theorem's by 0.0125 m/s. The cause is not the theorem. The zero-order hold feeds the loop a staircase that lags the true ramp by half a step on average, A·dt/2 = 0.5 × 0.05/2 = 0.0125 m/s. Halving dt halves the difference: 0.0250, 0.0125 and 0.00625 m/s at dt = 0.1, 0.05 and 0.025 s.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down first.

1. **Proportional gain sweep.** For K_p = 100, 400, 1000, 2000 and 4000 N per m/s record L(0), the step error, the least damping ratio and whether the loop is stable. Predict the error from 1/(1 + K_p/c) and the stability limit from (7).
2. **Integral gain sweep.** With K_p = 400, try K_i = 10, 40, 100 and 200 N/m. Record the ramp error for 0.5 m/s², the least damping ratio and the slowest real pole. Predict the ramp error from c·A/K_i, then predict where the damping collapses.
3. **Theorem against simulation.** Run `fvt_vs_simulation` for the five cases of `main()`, then repeat the ramp with dt = 0.1, 0.05 and 0.025 s. Predict the offset A·dt/2.
4. **A climb.** Entering a 2° climb at constant thrust is a drag-like force mg sin γ = 376.5 N. Predict the speed lost with no feedback, with K_p = 400 and with PI, then simulate and record the peak loss and the time to come within 0.1 m/s of the final value.

Raw data format: CSV tables `Kp, L0, step_error_pct, zeta_min, stable`; `Ki, ramp_error_m_s, zeta_min, slowest_real_pole`; `case, fvt_m_s, sim_m_s, dt_s`; `controller, final_loss_m_s, peak_loss_m_s, t_settle_s`. Each table's header names the plant parameters and says that the results are from a linear simulation.

What you should find: errors of 28.20, 8.94, 3.78, 1.93 and 0.97 % with damping 1.000, 0.722, 0.387, 0.194 and 0.041, all stable (K_crit = 4990). Ramp errors of 1.964, 0.491, 0.196 and 0.098 m/s with damping 0.726, 0.734, 0.453 and 0.146. The slowest pole moves from −0.024 to −0.127, −0.182 and back to −0.070 s⁻¹, so the fastest settling is at an intermediate K_i. The theorem matches the simulation to 10⁻⁶ except for the ramp offset of exactly A·dt/2. In the climb, the aircraft with no feedback would lose 9.59 m/s with a 28 s time constant. K_p = 400 loses 0.857 m/s for good (peak 0.889 m/s at 7.7 s). PI peaks at 0.782 m/s at 5.9 s and recovers fully; it is within 0.1 m/s of trim after 22.7 s, against 4.5 s for the P loop to reach its (wrong) final value.

## 5. Validate and improve

| quantity | analysis (§2) | code | simulation (600 s) | difference |
|---|---|---|---|---|
| step error, K_p = 400 | 1/(1 + 400/39.27) = 0.08940 | 0.08940 | 0.08940 | < 10⁻⁶ |
| drag 200 N, K_p = 400 | 200/(39.27 + 400) = 0.4553 m/s | 0.45530 | 0.45530 | < 10⁻⁶ |
| step and drag error, PI | 0 | 0 | < 10⁻⁶ | — |
| ramp 0.5 m/s², PI | c A/K_i = 0.49091 m/s | 0.49091 | 0.50341 | A·dt/2 = 0.0125 |
| critical gain | (7): 4989.7 N per m/s | bisection 4989.7 | — | < 10⁻⁸ relative |
| S + T | 1 | 1 ± 4.6 × 10⁻¹⁶ | — | rounding |

Analysis, code and simulation agree, and the one difference is explained to the last digit by the input model of the simulation. All three are statements about the *linear model*. Two things in the real aircraft break it first. Throttle saturation: the 0.5 m/s² ramp needs 550 N of extra thrust for the acceleration plus up to 393 N for the extra drag at +10 m/s. A small engine may not have it, and an integrator that keeps integrating against a saturated throttle winds up (`P09.11`). And the phugoid: in a real aircraft a speed error also changes lift and height, and the elevator loop's response feeds back into the speed. The answer for this model is **PI with K_p = 400 N per m/s and K_i = 40 N/m: zero step and drag error, 0.49 m/s lag on a 0.5 m/s² schedule, a gain margin of 10.9 (the whole controller can be multiplied by 10.9 before the loop goes unstable, found by `critical_gain` applied to L), at the price of an 18 % overshoot and a 23 s settling time.**

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1)–(7); `cruise_control.png`; the four experiment tables; and a half-page recommendation (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) of an auto-throttle gain set, with its step, ramp and disturbance errors, its damping, and a statement of what the model leaves out.

### Adversarial self-checks

- **What if the aircraft is twice as heavy?** m doubles, c is unchanged (same drag at the same speed), so L(0) and every steady-state error are unchanged, but the slow pole halves and the transients take twice as long. Steady-state error depends only on the low-frequency gain.
- **What if K → 0?** The step error tends to 1 (no correction) and the drag error to F/c, the open-loop loss.
- **Which assumption fails first?** Throttle saturation on large commands and climbs; then thrust independence of airspeed (a fixed-pitch propeller's thrust falls with speed, which adds to c and *reduces* the steady-state error a little).
- **What measurement could falsify the model?** In a flight test, a step in commanded speed with the auto-throttle in P mode: the residual error should be 1/(1 + K_p/c), from which c can be identified, and the time constant of the response compared with (2).
- **Does the result survive a change of time step?** The step and disturbance results are independent of dt, because the zero-order hold is exact for a constant input. The ramp result shifts by A·dt/2, as shown.
- **Which decision changes if the drag derivative c is twice as large?** L(0) = K_p/c halves, so the proportional step error grows from 8.94 to 16.4 %, while the drag error F/(c + K_p) falls slightly, from 0.455 to 0.418 m/s. K_v = K_i/c halves, so the PI ramp error doubles to 0.98 m/s; K_i would have to double to keep it.

## Check yourself

1. Write T, S and the disturbance path for unity feedback with controller C and plant G. What is S + T?
2. What is the system type, and how does it determine which inputs leave zero, finite or infinite steady-state error?
3. A type-0 loop has L(0) = 9. What fraction of a step command remains as error?
4. A type-1 loop has K_v = 2 s⁻¹. What is the steady-state error to a ramp of 0.5 m/s per second?
5. Under what condition is the final-value theorem valid, and what does it return if you ignore the condition for an unstable loop?
6. Why does a PI auto-throttle lose no speed after a constant drag increase, while a proportional one does?
7. Why can the step error of a proportional loop not be made arbitrarily small by raising the gain in a real plant with lags?
