# P09.10  Nyquist & Stability Margins

> Count how many times the open-loop frequency response winds around −1, add the open-loop unstable poles, and you have the number of unstable closed-loop poles, without solving for a single one. How far the curve stays from −1 is how much gain, phase and delay the loop can lose.

| | |
|---|---|
| **Learning outcomes** | 1. State the Nyquist criterion Z = N + P from the argument principle, build the Nyquist D-contour with indentations around imaginary-axis poles, and count encirclements of −1 by a winding number. 2. Verify Z = N + P against the closed-loop poles for stable and unstable gains, for a loop with a double integrator, for an open-loop unstable airframe, and for loops with a pure time delay (where there is no characteristic polynomial, checked through a Padé approximant). 3. Define and compute the gain margin, the phase margin and the delay margin from the phase and gain crossovers, and verify them against the crossings read from `scipy.signal.bode` and against the Routh critical gain. 4. Judge a loop's robustness from its margins, and explain what the margins miss. |
| **Prerequisites** | `P09.05` (Bode plots, magnitude and phase), `P09.08` (closed loop, sensitivity), `P09.09` (Routh, the critical gain, root locus), `P02.07` (complex numbers, phasors), `P02.08` (Laplace transform, the delay e^{−sT}). Skill check: you can sketch L(jω) for K/(s(s + 1)) in the complex plane and read its phase at ω = 1 rad/s. |
| **Engineering question** | A control-surface servo with loop L = K/(s(0.05s + 1)(0.01s + 1)), a reaction-wheel satellite with L = 0.3162(10s + 1)/(10s²(s + 1)) and a statically unstable airframe with L = K(s + 3)/((s − 2)(s + 10)): which are stable, how much gain, phase and delay can each lose before it goes unstable, and what changes when the open loop itself is unstable? Answer: the servo at K = 20 s⁻¹ has a gain margin of 6.0 (15.6 dB), a phase margin of 43.2° and a delay margin of 48.4 ms; the satellite has an infinite gain margin, 54.9° of phase margin and can tolerate 3.03 s of delay; the airframe needs one counterclockwise encirclement of −1 and is stable only above K = 6.67, so at K = 20 its gain may fall to a third (−9.5 dB) before it diverges. |
| **Assumptions and limits** | Linear time-invariant loops, unity negative feedback, rational L(s) possibly times a pure delay e^{−sT}. The Nyquist contour uses a finite radius R = 10⁴ rad/s and indentations of radius 10⁻⁴ around imaginary-axis poles; that is enough because every L here is strictly proper (it vanishes on the large arc). The margins are single-loop, one-parameter-at-a-time measures; they say nothing about simultaneous gain and phase changes (the subject of the disk margin in `P09.18`), about nonlinearity, or about multiple loops. The servo, satellite and airframe numbers are representative, not one vehicle's data. |
| **Analytical work** | The argument principle for 1 + L on a clockwise contour; the D-contour and the indentation's image (a large arc of angle −kπ for a k-fold pole at s = 0); Z = N + P; gain crossover \|L(jω_gc)\| = 1 and phase margin 180° + ∠L(jω_gc); phase crossover ∠L(jω_pc) = −180° and gain margin 1/\|L(jω_pc)\|; the delay margin PM/ω_gc; the servo's ω_pc = 1/√(τ_mτ_f) and GM = (τ_m + τ_f)/(Kτ_mτ_f) in closed form; the Padé approximant of e^{−sT}. |
| **Computational deliverable** | `starter/main.py`: `imag_axis_poles`, `nyquist_contour`, `loop_on_contour`, `winding_number`, `rhp_poles`, `nyquist_count`, `phase_margin`, `phase_crossovers`, `gain_margin`, `delay_margin`, `margins` (`poly_add`, `tf_eval`, `servo_loop`, `satellite_loop`, `unstable_airframe_loop`, `pade_delay`, `closed_loop_rhp_roots` and `main` are given). Tests prove: the contour is closed, clockwise, and keeps the indented poles outside; the winding number counts known circles (±1, 2, 0); Z = N + P equals the closed-loop RHP root count for six servo gains, five airframe gains (N = −1 exactly when K > 20/3) and five satellite delays (checked by an order-8 Padé closed loop); the Padé approximant is all-pass and matches e^{−jωT} at low frequency; the gain and phase margins match the `scipy.signal.bode` crossings to 0.01 dB and 0.01°; GM equals the Routh critical gain over K; the satellite PM matches its closed form and its GM is infinite; the delay margin puts the curve through −1 to 10⁻⁹; half the delay margin costs half the phase margin; at small K the phase margin tends to 90°, and at the critical gain GM = 1 and PM = 0. |
| **Practical activity** | **computational experiment**: a gain sweep of the servo (margins, delay margin, peak sensitivity and Z against the closed-loop poles), a delay sweep of the satellite (Nyquist with an exact delay against a Padé closed loop), and a gain sweep of the open-loop unstable airframe. |
| **Equipment** | None. A computer with the environment from `P00.01`. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `nyquist_margins.png`, a hand sketch of the servo's Nyquist plot including the image of the indentation, the derivations of Z = N + P and of the servo's GM in closed form, the three sweep tables, and a half-page robustness note on one of the loops. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: the criterion and the servo's margins are derived; all public tests pass; the note states GM, PM and DM with numbers, says what each protects against, and names one thing the margins miss. |
| **References** | `@mit-ocw-16-06-principles-of-automatic-control` Lecture 20 "The Nyquist stability criterion" (PDF "Bode Plots With Complex Poles": neutral stability; The Argument Principle), Lecture 21 "The Nyquist stability criterion (cont.)" (P, the number of unstable poles; Examples), Lecture 22 "Nyquist with poles on imaginary axis" (also: Nyquist Plot of Open Loop Unstable System) and Lecture 23 "Stability margins, Bode gain-phase theorem" (Stability Margins; Relationship between PM and damping) (required-free): the argument principle and Z = N + P, indenting the contour around imaginary-axis poles, the open-loop unstable case, and gain and phase margins from the −1 point. `@mit-ocw-16-30-feedback-control-systems` Topic 3 "Frequency response methods": Frequency Stability Tests; Gain and Phase Margins; Nyquist Stability; FR: Summary (required-free): the criterion and the margins in compact form. `@mit-ocw-6-003-signals-and-systems` Lecture 10 "Feedback and Control": Analysis of wallFinder System: Adding Sensor Delay; Destabilizing Effect of Delay (required-free): how a delay in the loop destabilises it. None of these notes defines the delay margin (PM divided by crossover frequency) or counts encirclements numerically by a winding number: both are original instruction in this lesson. |
| **Extensions** | ★ Compute the peak sensitivity M_s = max 1/\|1 + L(jω)\| for the servo sweep and show that GM ≥ M_s/(M_s − 1) and PM ≥ 2 arcsin(1/(2M_s)) (a preview of `P09.18`). ★ Add a zero-order hold and a one-sample computation delay at 200 Hz to the servo and find the gain that keeps 45° of phase margin. ★ Draw the Nyquist plot of an inverted-pendulum launch stand under PD control and find the range of K for stability. |
| **Time** | unmeasured |

## 1. The engineering problem

A designer who has a model can compute the closed-loop poles and know whether the loop is stable. But a real loop contains things a polynomial cannot hold: a computation delay, a transport delay, a sensor with a frequency response measured on a bench rather than derived. What is available for every loop is the **open-loop frequency response** L(jω), the response measured with the loop broken. Nyquist's criterion answers the stability question from that curve alone. It does so for a loop with a delay, for a loop whose open-loop plant is itself unstable, and for a loop whose L(jω) is a table of measurements.

It also answers the question that matters more in practice: *how stable*. A loop that is stable on paper but sits close to instability will go unstable when the actuator ages, the payload changes or the computer gets busy. The distance of the curve L(jω) from the point −1 measures that. The **gain margin** says how much more loop gain is tolerated, the **phase margin** how much more phase lag, the **delay margin** how much more time delay. Typical aerospace requirements are a gain margin of at least 6 dB and a phase margin of at least 45° on each loop (flight-control specifications state them per loop and per flight condition). This lesson builds the Nyquist plot from scratch, counts encirclements, computes the three margins and checks every result two independent ways.

Success criteria, written first: Z = N + P agrees with the closed-loop poles in every one of 16 cases, including open-loop unstable and delayed loops; the gain and phase margins agree with the crossings of `scipy.signal.bode` to 0.01 dB and 0.01°; the gain margin agrees with the Routh critical gain to 10⁻⁸; the delay margin puts L(jω)e^{−jωT} through −1 to 10⁻⁹.

## 2. Understand and derive

### 2.1 The argument principle

Let F(s) be rational and C a closed curve, traversed clockwise, that passes through no zero or pole of F. As s goes once around C, the angle of F(s) changes by −2π(Z_C − P_C), where Z_C and P_C are the numbers of zeros and poles of F inside C. Each zero inside contributes one full clockwise turn of its factor (s − z) and each pole one counterclockwise turn. Factors whose roots lie outside C return to their starting angle. So **the image F(C) encircles the origin clockwise Z_C − P_C times.**

### 2.2 Nyquist's criterion

Apply it to F = 1 + L (16.06 Lectures 20 and 21). The zeros of 1 + L are the closed-loop poles; its poles are the open-loop poles. Choose C to enclose the whole right half plane: the **Nyquist D-contour**, up the imaginary axis from −jR to +jR and back along a semicircle of radius R → ∞. Encircling the origin by 1 + L is the same as encircling −1 by L. Hence

**Z = N + P,**   (1)

with Z the number of closed-loop poles in the right half plane, P the number of open-loop poles there, and N the number of **clockwise** encirclements of −1 by L(s) on the contour. For stability Z = 0, so N = −P: an open-loop stable loop must not encircle −1, and an open-loop unstable loop must encircle it counterclockwise once per unstable pole.

**Indentations.** The contour must not pass through a pole of L (16.06 Lecture 22). A servo's integrator (pole at 0) and a satellite's double integrator lie on the imaginary axis, so the contour detours around them on small semicircles s = jω_k + εe^{jφ}, φ from −90° to +90°, to the right, which keeps them outside (they are not counted in P). Near a k-fold pole at the origin, L ≈ c/s^k, and on the small semicircle L ≈ (c/ε^k)e^{−jkφ}: **a huge arc swinging clockwise through kπ**. For the servo (k = 1) that is a half circle at infinity from +90° through 0° to −90° (the image of φ = −90° … 90°); for the satellite (k = 2), a full circle. Forgetting the indentation is the most common error in Nyquist plots, and the winding count is wrong without it.

**On the big arc**, every strictly proper L tends to 0, so it maps to the origin and contributes nothing. With a delay e^{−sT}, \|e^{−sT}\| = e^{−T Re s} ≤ 1 on the right half of the arc, so the same holds.

### 2.3 Counting encirclements by a winding number

Sample the contour densely, evaluate z_k = L(s_k), and add up the angle increments of z_k + 1 between consecutive points, each wrapped to (−π, π]: arg((z_{k+1} + 1)/(z_k + 1)). The total divided by 2π is the counterclockwise winding number; N is minus that. The sum is exact as long as consecutive points never turn by more than π about −1, which a dense, log-spaced grid guarantees except when the curve passes almost through −1 (a loop at its stability limit, where the count is genuinely ambiguous).

### 2.4 The three margins

Write L(jω) = \|L\|e^{jφ} (16.06 Lecture 23; 16.30 Topic 3, "Gain and Phase Margins"). For a loop that is open-loop stable and stable as designed:

- **Gain crossover** ω_gc: \|L(jω_gc)\| = 1. **Phase margin PM = 180° + φ(ω_gc)**: the extra phase lag at ω_gc that would put the curve through −1.
- **Phase crossover** ω_pc: φ(ω_pc) = −180° (the curve crosses the negative real axis). **Gain margin GM = 1/\|L(jω_pc)\|**: the factor by which the gain can rise before the crossing reaches −1. In decibels, GM_dB = −20 log₁₀\|L(jω_pc)\|.
- **Delay margin DM = PM (in radians)/ω_gc**: a delay T multiplies L by e^{−jωT}, which leaves \|L\| unchanged and subtracts ωT of phase. At ω_gc the phase margin is used up when ω_gc T = PM. If the curve crosses \|L\| = 1 more than once, the delay margin is the smallest of PM_i/ω_gc,i.

**The servo in closed form.** L = K/(s(τ_m s + 1)(τ_f s + 1)) has phase −90° − arctan(τ_mω) − arctan(τ_fω). It reaches −180° where arctan(τ_mω) + arctan(τ_fω) = 90°, i.e. τ_mτ_fω² = 1: **ω_pc = 1/√(τ_mτ_f) = 44.72 rad/s** for τ_m = 0.05 s and τ_f = 0.01 s. There \|L\| = Kτ_mτ_f/(τ_m + τ_f), so

**GM = (τ_m + τ_f)/(Kτ_mτ_f) = K_crit/K,  K_crit = 120 s⁻¹,**   (2)

exactly the Routh limit of `P09.09` for the cubic τ_mτ_f s³ + (τ_m + τ_f)s² + s + K. Gain margin and Routh are two views of the same boundary.

### 2.5 Three loops

**Servo** (from `main()`):

| K [s⁻¹] | N | Z (roots) | GM | PM at ω_gc | DM |
|---|---|---|---|---|---|
| 5 | 0 | 0 (0) | 24.0 (27.60 dB) | 73.58° at 4.853 rad/s | 264.61 ms |
| 20 | 0 | 0 (0) | 6.00 (15.56 dB) | 43.21° at 15.587 | 48.38 ms |
| 60 | 0 | 0 (0) | 2.00 (6.02 dB) | 15.55° at 31.038 | 8.75 ms |
| 119 | 0 | 0 (0) | 1.0084 (0.07 dB) | 0.18° at 44.535 | 0.07 ms |
| 121 | 2 | 2 (2) | 0.9917 (−0.07 dB) | −0.18° at 44.907 | — |
| 200 | 2 | 2 (2) | 0.600 (−4.44 dB) | −10.53° at 57.247 | — |

The phase crossover is at 44.721 rad/s for every K (it depends only on the lags). Above K = 120 the curve wraps around −1 twice clockwise and two closed-loop poles are unstable. A negative margin on a loop already unstable says by how much it is unstable, not how much it can lose.

**Satellite**: L = K_p(T_d s + 1)/(J s²(τs + 1)) with J = 10 kg m², T_d = 10 s, τ = 1 s, K_p = 0.3162 N m/rad. The phase is −180° + arctan(T_dω) − arctan(τω), above −180° for every ω > 0, so the curve never crosses the negative real axis: **GM = ∞**. The double pole at the origin makes the indentation's image a full clockwise circle at infinity, which together with the rest of the curve leaves −1 unencircled: N = 0. The lead between the zero at 0.1 and the pole at 1 rad/s peaks at ω = 1/√(T_dτ) = 0.3162 rad/s, where K_p was chosen to put the crossover. **PM = 54.90° at 0.3162 rad/s, DM = 0.9582 rad/0.3162 rad/s = 3.030 s.** A gain margin of infinity does not mean the loop is safe: the satellite has a finite delay margin, and the delay is what a slow on-board computer adds.

| delay T | N (exact delay) | Z | Padé(8) closed loop RHP roots | PM with the delay |
|---|---|---|---|---|
| 0 | 0 | 0 | 0 | +54.90° |
| 0.500 s | 0 | 0 | 0 | +45.84° |
| 2.727 s (0.9 DM) | 0 | 0 | 0 | +5.49° |
| 3.333 s (1.1 DM) | 2 | 2 | 2 | −5.49° |
| 6.061 s (2 DM) | 2 | 2 | 2 | −54.90° |

At T = DM, \|1 + L(jω_gc)e^{−jω_gcT}\| = 3.5 × 10⁻¹⁶: the curve passes through −1.

**Unstable airframe**: L = K(s + 3)/((s − 2)(s + 10)), P = 1. The closed loop s² + (8 + K)s + (3K − 20) is stable for **K > 20/3 = 6.667**. Nyquist must show N = −1 there:

| K | L(0) | N | P | Z (roots) | PM | gain change to instability |
|---|---|---|---|---|---|---|
| 3 | −0.450 | 0 | 1 | 1 (1) | — | must rise × 2.222 |
| 6 | −0.900 | 0 | 1 | 1 (1) | — | must rise × 1.111 |
| 7 | −1.050 | −1 | 1 | 0 (0) | 34.26° at 0.858 rad/s | may fall to 0.952 (−0.42 dB) |
| 20 | −3.000 | −1 | 1 | 0 (0) | 103.49° at 17.505 rad/s | may fall to 0.333 (−9.54 dB) |
| 100 | −15.000 | −1 | 1 | 0 (0) | 92.86° at 99.524 rad/s | may fall to 0.067 (−23.52 dB) |

The curve starts at L(0) = −3K/20 on the negative real axis. Once that point lies to the left of −1, the curve encircles −1 counterclockwise once and the loop is stable. The relevant margin is a **lower** gain margin: how far the gain may *fall*. That is the Nyquist view of the root-locus fact in `P09.09` Problem 3.

### 2.6 Where the margins stop

A gain margin and a phase margin are two points on the curve. A curve can pass both tests and still come close to −1 somewhere else, which a simultaneous small change of gain and phase would exploit. The servo at K = 20 has GM = 6 and PM = 43°, yet its closest approach to −1 is 1/M_s with peak sensitivity M_s = 1.77 (Problem 4). Margins assume a single loop and linear behaviour; an actuator rate limit acts like an extra, amplitude-dependent delay and can erode them in large manoeuvres. And the margins are computed on the model: the 0.01 s filter lag is as uncertain as any number in it, and doubling it cuts the critical gain from 120 to 70 s⁻¹, roughly in half (Problem 6).

## 3. Build the computational model

Open `starter/main.py`. The three loops, the Padé approximant and `main()` are given.

### Step 1: the contour

`imag_axis_poles(den)` returns the distinct imaginary parts of the roots of den with \|Re\| < 10⁻⁸. `nyquist_contour(den, R, eps)` builds the axis points uniformly in asinh(ω) (logarithmic for large \|ω\|, linear near 0), adds log-spaced points near each indentation, removes the points within eps of each pole, and stitches the pieces together:

```python
for x in list(wk) + [None]:
    hi = R if x is None else x - eps
    pieces.append(1j * w[(w >= lo) & (w <= hi)])                 # up the axis
    if x is not None:
        pieces.append(1j * x + eps * np.exp(1j * phi))           # around the pole, to the right
        lo = x + eps
pieces.append(R * np.exp(1j * np.linspace(np.pi / 2, -np.pi / 2, n_arc)))   # back along the arc
```

### Step 2: counting

`loop_on_contour` evaluates L(s)e^{−sT}. `winding_number(z, point)` is `np.sum(np.angle(d[1:] / d[:-1])) / (2*np.pi)` on d = z − point with the first point appended to close the curve. `rhp_poles` counts roots with Re > 10⁻⁹. `nyquist_count` returns (N, P, N + P) with N = −round(winding).

### Step 3: margins

`phase_margin` finds the first sign change of log\|L(jω)\| on a log grid from 10⁻⁴ to 10⁴ rad/s, refines it by bisection in log ω, and wraps 180° + ∠L to (−180°, 180°]. `phase_crossovers` finds the sign changes of Im L where Re L < 0 and refines each by bisection; crossings with \|L\| < 10⁻³ (a margin above 60 dB) are skipped, because a delayed loop spirals through thousands of them on its way to the origin. `gain_margin` takes the crossing with the largest \|L\|. `delay_margin` is PM/ω_gc.

Running `main()` prints the three tables of 2.5 and writes `nyquist_margins.png` (the three Nyquist curves near −1 with the unit circle; the satellite with no delay and with T = DM).

### Verify

Run `python -m pytest curriculum/phase-09/lesson-10-nyquist-stability-margins`. The 28 tests check:

1. **Contour**: never within ε of the double pole, starts and ends at −jR, winds −1 around a right-half-plane point and 0 around left-half-plane points and the indented pole.
2. **Winding number** of circles about −1: +1, −1, +2, and 0 for a circle elsewhere; imaginary-axis poles found.
3. **Z = N + P against the closed-loop poles** (independent check 1): six servo gains (Z = 0 below 120, 2 above); five airframe gains (P = 1, N = −1 exactly when K > 20/3); five satellite delays from 0 to 2 DM against an order-8 Padé closed loop.
4. **Padé**: matches e^{−jωT} to 10⁻⁹ for ωT ≤ 1 and has \|·\| = 1 at all frequencies.
5. **Margins against `scipy.signal.bode`** (independent check 2) for three servo gains (ω_gc, PM, ω_pc, GM_dB) and for the satellite (PM, ω_gc).
6. **Gain margin against Routh**: GM = 120/K and ω_pc = 1/√(τ_mτ_f) to 10⁻⁸.
7. **Satellite**: GM infinite; PM equal to arctan(T_dω) − arctan(τω) at ω_gc.
8. **Delay margin**: \|1 + L e^{−jω_gc DM}\| < 10⁻⁹; half the delay costs exactly half the phase margin.
9. **Limits**: at K = 0.5 the crossover is at ω ≈ K and PM = 90° − arctan(τ_mω) − arctan(τ_fω) > 88°; at K = 120, GM = 1 and PM = 0.

### Compare

`scipy.signal.bode` returns magnitude in dB and unwrapped phase in degrees on a grid; reading its crossings by interpolation on 200 001 points agrees with our bisection to 0.01 dB and 0.01°. `bode` is a plotting tool, though: it does not count encirclements, does not handle delays, and its phase unwrapping can start at a different multiple of 360°. That is why the test compares phase margins modulo 360°. For the delayed loops, the Padé approximant turns e^{−sT} into a rational function so that `np.roots` can find closed-loop poles. That check is independent of the Nyquist count. It is approximate: order 8 has a phase error of 0.0001° at ωT = 6, 0.01° at ωT = 8 and 2.6° at ωT = 12, while the crossovers here have ωT ≤ 2.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Sketch each Nyquist plot by hand first, including the image of the indentation, and write each prediction down before running.

1. **Servo gain sweep.** For K = 5, 10, 20, 40, 60 and 100 s⁻¹, record GM (and dB), PM, ω_gc, DM, the peak sensitivity M_s = max 1/\|1 + L\| and Z against the roots. Predict GM from (2) and the trend of PM and M_s.
2. **Satellite delay sweep.** For T = 0, 0.5 s, 0.9 DM, 1.1 DM and 2 DM, record N, Z and the Padé closed-loop count. Predict the change of Z at T = DM and the phase margin with delay from PM − ω_gcT.
3. **Satellite gain sweep.** For K_p = 0.1, 0.3162, 1 and 3 N m/rad, record PM, ω_gc and DM. Predict which K_p maximises PM and which maximises DM.
4. **Airframe gain sweep.** For K = 3, 6, 7, 20 and 100, record L(0), N, Z and the roots. Predict the K at which N changes.

Raw data format: CSV tables `K, GM, GM_dB, PM_deg, w_gc, DM_ms, Ms, N, Z, Z_roots`; `delay_s, N, Z, Z_pade`; `Kp, PM_deg, w_gc, DM_s`; `K, L0, N, P, Z, Z_roots`, each with a header naming the loop.

What you should find: the servo table of 2.5 extended (Problem 4) by GM 12.0 (21.58 dB), PM 60.42°, DM 116.3 ms, M_s 1.390 at K = 10; GM 3.0, PM 25.39°, DM 18.06 ms, M_s 2.716 at K = 40; GM 1.2, PM 3.94°, DM 1.69 ms, M_s 15.74 at K = 100; and M_s 1.206, 1.769 and 4.207 at K = 5, 20 and 60. The delay margin falls much faster than the gain margin: from K = 5 to K = 60 GM falls by 12 while DM falls by 30. The satellite tables of 2.5, and for the gain sweep: PM 44.46°, 54.90°, 44.46° and 28.48° at ω_gc = 0.1265, 0.3162, 0.7907 and 1.596 rad/s, with DM 6.14, 3.03, 0.98 and 0.31 s. The phase margin peaks at K_p = 0.3162, where ω_gc sits at the lead's maximum, but the delay margin keeps growing as K_p falls, because ω_gc falls. The airframe table of 2.5: N jumps from 0 to −1 between K = 6 and 7, at L(0) = −1.

## 5. Validate and improve

| quantity | analysis (§2) | code | independent check | difference |
|---|---|---|---|---|
| servo ω_pc | 1/√(τ_mτ_f) = 44.721 rad/s | 44.721 | Bode crossing | < 10⁻⁴ relative |
| servo GM, K = 20 | (2): 120/20 = 6.000 (15.56 dB) | 6.0000 | Bode: 15.56 dB; Routh K_crit = 120 | < 0.01 dB |
| servo PM, K = 20 | — | 43.21° at 15.587 rad/s | Bode crossing | < 0.01° |
| servo Z, K = 121 | Routh: 2 | N = 2, Z = 2 | `np.roots`: 2 | 0 |
| satellite PM | arctan(10ω) − arctan(ω) at 0.3162 = 54.90° | 54.90° | Bode | < 10⁻⁶ ° |
| satellite DM | 0.9582/0.3162 = 3.030 s | 3.030 s | \|1 + Le^{−jωT}\| = 3.5 × 10⁻¹⁶ | — |
| satellite Z at 1.1 DM | 2 | Nyquist (exact delay) 2 | Padé(8) roots: 2 | 0 |
| airframe stability | K > 20/3 | N = −1 for K = 7, 20, 100 | roots | 0 |

Every Nyquist count agrees with the closed-loop poles, and every margin with the Bode crossings and with Routh. That establishes the method on these models; it does not make the models true. The engineering conclusions: **the servo at K = 20 s⁻¹ meets a 6 dB / 45° requirement on gain but not quite on phase (43.2°), and its 48 ms delay margin is the number to watch if the servo is driven from a slow digital bus. The satellite's infinite gain margin is irrelevant; its 3.0 s delay margin is ample for a 10 Hz attitude computer. The unstable airframe at K = 20 has a 9.5 dB lower gain margin, and that, not an upper margin, is its safety limit.**

What to improve: report M_s alongside GM and PM (extension 1), and include the digital implementation's delay (one sample plus half a sample for the zero-order hold) in L before computing the margins (extension 2 and `P09.17`).

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the hand sketch of the servo's Nyquist plot made before running the code; the derivations of (1) and (2); `nyquist_margins.png`; the four experiment tables; and a half-page robustness note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) on one loop, giving GM, PM, DM and M_s and the physical change each protects against.

### Adversarial self-checks

- **What if the servo's filter lag doubles to 0.02 s?** K_crit falls from 120 to 70 s⁻¹. At K = 20 the GM falls from 6.0 to 3.5 (10.9 dB), the PM from 43.2° to 35.8° and the DM from 48 to 41 ms. A "small" sensor lag decides the margins.
- **What if a 20 ms computation delay is added to the servo at K = 20?** The gain crossover stays at 15.59 rad/s and the PM falls by 0.02 × 15.59 rad = 17.9° to 25.3°. The phase crossover moves down to 23.6 rad/s and the GM falls to 1.87 (5.5 dB). Both now fail a 6 dB / 45° requirement.
- **Which assumption fails first?** Linearity: actuator rate and position limits act like extra delay in large commands. Then the single-loop view: in an aircraft the servo loop sits inside an attitude loop whose margins it affects.
- **What measurement could falsify the margins?** Inject a sine into the loop on a test rig, measure L(jω) by breaking it at the actuator command, and compare the measured crossings with the predicted ones; or raise the gain on the rig until it oscillates, which should happen at 120 s⁻¹ and 44.7 rad/s (7.1 Hz).
- **Does the result survive a change of contour?** Halving the indentation radius, doubling R, or both, changes none of the 16 counts of 2.5, the near-marginal servo gains 119 and 121 included.
- **Which decision changes if the margin requirement rises to 60°?** The servo gain must fall to about 10 s⁻¹ (PM 60.4°), halving the bandwidth; the satellite's 54.9° is the maximum its lead network can give, so its compensator, not its gain, must change.

## Check yourself

1. State the Nyquist criterion. What are Z, N and P, and in which direction is N counted?
2. Why does the Nyquist contour detour around a pole at s = 0, and what does the detour look like after mapping through L for a single integrator?
3. Define the gain margin and the phase margin in terms of the crossover frequencies.
4. A loop has PM = 45° at ω_gc = 10 rad/s. What is its delay margin?
5. The servo L = 20/(s(0.05s + 1)(0.01s + 1)): at what frequency is the phase −180°, and what is the gain margin in dB?
6. An open-loop unstable airframe has P = 1. How many times, and in which direction, must the Nyquist curve encircle −1 for stability?
7. Why can a loop with an infinite gain margin still go unstable from a delay?
8. Give one way a loop can satisfy GM ≥ 6 dB and PM ≥ 45° and still be close to instability.
