# P09.21  Home Lab: Servo Feedback Control 🔬

> A model earns trust by predicting something before you measure it. Identify the actuator from a step, write the closed-loop rise time, overshoot and settling time down with the date, then close the loop and find out which effect you forgot.

| | |
|---|---|
| **Learning outcomes** | 1. Identify a first-order-plus-dead-time model K e^{−Ls}/(τs + 1) and a first-order-plus-integrator model K/(s(τs + 1)) with an input dead zone from open-loop step tests, by the log-linear method and by least squares, and report each parameter with a standard error. 2. Design a P gain for an overshoot target and a PI gain for a rise-time target from the identified model, and predict the 10–90 % rise time, the percent overshoot and the 2 % settling time, by the second-order formulas and by simulation. 3. Close a real feedback loop on a hobby servo or a small gear motor with a microcontroller, inside a written safety envelope with a software duty (or command) limit and a stop button tested first, and with the prediction written and dated before the loop is closed. 4. Extract the metrics from sampled data, compare them with the prediction, and explain each discrepancy by a named effect (sample delay, dead zone, saturation, quantisation, friction, slew limit) that a refined model reproduces. 5. Redesign on the refined model and re-measure. |
| **Prerequisites** | `P09.04` (step-response metrics; ζ and ω_n), `P09.08` (closed-loop transfer function, system type, steady-state error), `P09.11` (PI control, actuator saturation, clamping anti-windup), `P09.17` (sampling, zero-order hold, delay in a digital loop), `P04.15` (the bench microcontroller, its ADC, quantisation LSB/√12; the p04 safety box), `P00.13` (safety box), `P00.14` (uncertainty, least-squares fits), `P03.14` (Gauss–Newton and damped least squares). Skill check: you can write the closed-loop transfer function of K/(s(τs + 1)) under a gain Kp and read ζ and ω_n from it. |
| **Engineering question** | From three or four open-loop step tests of a small actuator, can you predict, before closing the loop, the rise time, overshoot and settling time that a P or PI controller will give to within 15 % (and 3 percentage points of overshoot), and when the measurement disagrees, can you name and model the effect responsible? |
| **Assumptions and limits** | The actuator is linear about its operating point apart from a static input dead zone and a hard output limit; the potentiometer (or encoder) is a linear, noise-and-quantisation-limited angle sensor; the gearbox has no backlash; the battery voltage is constant during a session (K scales with it); the servo's internal loop behaves as a first-order lag with a frame latency (it also has a slew limit and a dead band, which the linear model omits); the continuous design ignores the loop's sampling. The models fail for steps large enough to saturate (the motor at 60°, the servo's slew limit at 40°), for very small steps inside the dead zone, and wherever stiction differs from the dead-zone picture. Battery or USB power only, motor supply ≤ 6 V, nothing on the shaft except a light pointer. |
| **Analytical work** | The step response of K/(τs + 1) with dead time and its log-linear form; the step response of K/(s(τs + 1)), its asymptote with lag τ + L and the exponential remainder; linearised least-squares standard errors; the dead zone from the slope–duty line; the P loop's ζ = 1/(2√(τKpK)) and ω_n = √(KpK/τ) and the gain for a given overshoot; PI pole-zero cancellation and the first-order closed loop; the exact second-order overshoot and peak time; why 4/(ζω_n) overestimates the settling time; the loop delay 1.5 Ts + L as a phase loss ω_c θ and as a first-order lag; the integrating loop with a lag as a second-order system with ζ = 1/(2√(θω_c)); the dead-zone steady-state error band; the saturated slew rate. |
| **Computational deliverable** | `starter/main.py`: `first_order_step_response`, `integrator_step_response`, `fit_first_order_loglinear`, `fit_integrator_loglinear`, `refine_least_squares`, `fit_dead_zone`, `step_metrics`, `second_order_step`, `second_order_metrics`, `zeta_from_overshoot`, `overshoot_from_zeta`, `design_p_position`, `design_pi_first_order`, `closed_loop_tf`, `predict_step`, `predict_closed_loop`, `simulate_closed_loop`, `simulate_open_loop`, `compare_metrics`, `design_p_sampled` (the firmware law `pi_update`, `tf_to_ss`, `rk4_map`, `quantize`, `dead_zone`, the SYNTHETIC bench, the CSV writer and reader, `analyse_dataset` and `main` are given). Tests (30) prove: the step models' limits and agreement with `scipy.signal.step`; identification recovering known K, τ, L to 10⁻⁶ from noiseless data, from a poor start by least squares, and within 4 standard errors from seeded noisy, quantised data; the dead zone recovered exactly; sampled metrics against the closed-form second-order overshoot and peak time and the exact rise and settling times, and the first-order τ ln 9 and τ ln 50; the P and PI designs meeting their targets; the closed loop against `scipy.signal.step` and against C G/(1 + C G) evaluated pointwise; a fast ideal sampled loop reproducing the prediction; the open-loop simulator against the closed form; the duty limit never exceeded; quantisation, seeding and zero gain; the one-sample delay arriving one sample late and costing damping; the sampled redesign hitting its target; the lesson's independent check (below); the SYNTHETIC files regenerated byte-identically from the seed and equal to the shipped copies. |
| **Practical activity** | **physical experiment**: wire a feedback servo (Option A) or a gear motor with a driver and a potentiometer (Option B) to the p04 microcontroller, run the four safety tests, record open-loop steps, identify the model, write and date the closed-loop prediction, close the loop, measure, explain, redesign and re-measure. See `projects/p09-servo-control/`. A learner without the kit completes the lesson with the SYNTHETIC dataset in `projects/p09-servo-control/data/synthetic/`, labelled synthetic everywhere. |
| **Equipment** | Level 2. Minimum (Option A): an analog-feedback micro servo, a 4 × AA NiMH pack with a 0.5 A resettable fuse, a push button, a capacitor, a card pointer and a printed protractor, plus the p04 board, breadboard, multimeter and USB cable. Recommended (Option B): an N20-class gear motor (stall ≤ 500 mA at 6 V), a low-voltage H-bridge module, the p04 10 kΩ potentiometer on the shaft. Enhanced: an encoder gear motor, a current-sense breakout, slow-motion phone video as an independent angle check. Tiers, substitutes and prices in `projects/p09-servo-control/bom.yaml` (reuses `p04-electronics-bench`). |
| **Safety** | See `projects/p09-servo-control/safety.md` (authoritative). Summary box at the start of part 4. Hard limits: battery or USB power only, never mains; motor supply ≤ 6 V; one motor, stall ≤ 500 mA, fused; nothing on the shaft except a light pointer arm (no propellers, no blades); fingers clear while the loop runs; a software duty or command limit and a stop button, tested before the first closed-loop run. |
| **Evidence** | `data/safety_tests.md` (the four tests, dated), `data/calibration.csv`, `data/prediction.md` (the closed-loop prediction dated before the first closed-loop run, and a second one before the redesign), `data/raw/*.csv`, the identification table with standard errors and repeat spread, the comparison tables, `servo_feedback_lab.png` from your data (or from the synthetic data, labelled as such), the ranked discrepancy causes with the refined-model result, and a report following `projects/p09-servo-control/report-template.md`. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; all four competencies. Threshold: the four safety tests logged before the first closed-loop run; a prediction of rise time, overshoot and settling time dated before the loop was closed; measured metrics with repeat spread; every discrepancy beyond that spread explained by a named effect that a refined model reproduces; a redesign with its own dated prediction and measurement; the completed safety box. A report built on the synthetic dataset must say so in its title and summary and must not call anything measured. |
| **References** | `@sparkfun-learn-tutorials` (required-free; CC BY-SA 4.0; ignore the purchase sections; the examples use Arduino, not MicroPython): "Hobby Servo Tutorial" (https://learn.sparkfun.com/tutorials/hobby-servo-tutorial/all): Servo Motor Background (Electrical Connection; Control signal; Powering Servos; Show Me The Guts) and A Handful of Distinctions (Range Constraints; Analog vs. Digital): the servo's own potentiometer-and-amplifier loop that Option A wraps (§2.1). "Pulse Width Modulation": Duty Cycle: the duty as the actuator command. "TB6612FNG Hookup Guide": Board Overview; Hardware Setup, and "Motors and Selecting the Right One": DC Brush Motors - The Classic (supplementary-free): the Option B driver and motor. "Analog to Digital Conversion": "Relating ADC Value to Voltage": counts to a physical quantity, used for the angle calibration. `@micropython-docs` (required-free; the lesson does not require that board): "Quick reference for the RP2" (https://docs.micropython.org/en/latest/rp2/quickref.html): Timers; PWM (pulse width modulation); ADC (analog to digital conversion); "class PWM – pulse width modulation": Methods; Limitations of PWM; "class Timer – control hardware timers": a fixed-rate control-loop tick; "class ADC – analog to digital conversion": Methods (`ADC.read_u16()`): the calls the firmware uses. `@mit-ocw-6-003-signals-and-systems` Lecture 12 "Continuous-Time (CT) Feedback and Control, Part 1": Motor Controller; Dominant Pole (supplementary-free): a feedback motor controller and the dominant-pole view of its closed loop. `@mit-ocw-16-06-principles-of-automatic-control` Lecture 6 "Time-domain specifications" (supplementary-free review): rise time, overshoot and settling time against ζ and ω_n. Partial GAP (original instruction): no verified free source gives a step-test identification procedure (fitting K, τ and L, or K/(s(τs + 1)), from measured data) or a home closed-loop servo experiment with predictions written before the run; the identification method, the P/PI design to a stated target, the procedure, the safety limits and the synthetic dataset are written here, building on `P09.04`, `P09.08`, `P09.11` and `P09.17` and following p04-electronics-bench practice. |
| **Extensions** | ★ Option B with an encoder (enhanced tier): no end stops and a resolution you can compute from the encoder's counts per revolution; repeat the redesign at 200 Hz and measure how much of the delay-induced overshoot disappears. ★ Add the current-sense breakout, log the motor current through the 60° step, and add a firmware current limit (for example 400 mA) beside the duty limit; show from the log which one binds. ★ Film the pointer at 240 frames per second against the protractor and compare the video rise time with the potentiometer's: an independent check of the sensor. ★ Identify the servo at three step sizes and fit the slew limit as a second parameter. |
| **Time** | unmeasured |

## 1. The engineering problem

Every control surface on a modern aircraft is moved by an actuator in a feedback loop: the flight computer commands an angle, a position sensor on the actuator reports the angle reached, and a controller drives the difference to zero. Before an aircraft flies, those loops are tested on the ground on an "iron bird" rig, and the measured step responses are compared with what the control-law designers predicted. A loop that overshoots more than predicted is a loop with less stability margin than anyone thought; a loop that rises more slowly is a loop that will lag the pilot. Either way the model the control law was designed on is wrong somewhere, and the disagreement is the most useful number on the test sheet, because it points to what the model left out.

This lab is the smallest honest version of that test. You will take a hobby servo with its feedback potentiometer brought out (Option A), or a small DC gear motor with a driver and a potentiometer on its shaft (Option B), and treat it as an unknown plant. You will step it open loop, fit a model, design a gain to a stated target, **write the predicted rise time, overshoot and settling time down and date them**, and only then close the loop. Then you will measure, and where the measurement and the prediction disagree, you will name the effect (the one-sample delay of the firmware, the dead zone, the duty limit, quantisation, the servo's own slew limit), add it to the model, and show that the refined model predicts what you measured. Finally you will redesign on the refined model and measure again.

Success criteria, written first: the four firmware safety tests logged before any closed loop; the identified parameters with standard errors and a repeat spread; a dated prediction; the measured metrics with their spread over three repeats; a refined model whose rise and settling times agree with the measurement within 15 % and whose overshoot agrees within 3 percentage points; every larger difference explained. A learner without the kit meets the same criteria on the SYNTHETIC dataset, says so on every page, and does not claim a measurement.

## 2. Understand and derive

### 2.1 Two plants

```
 Option A (feedback servo)              Option B (gear motor + driver + pot)
 command c ──► [servo's own loop] ──► θ  duty u ──► [dead zone] ──► K/(τs + 1) ──► ω ──► 1/s ──► θ
               K e^{−Ls}/(τs + 1)
 outer loop:  r ──►(+)──► PI ──► c       outer loop:  r ──►(+)──► P (+ dead-zone comp.) ──► u
                  ▲−                         ▲−
                  └── pot wiper ◄── θ         └── pot ◄── θ
```

In **Option A** the servo already contains a motor, a gearbox, a potentiometer and its own analog loop (SparkFun's "Hobby Servo Tutorial", "Show Me The Guts", opens one up). From the outside, a commanded angle c (a pulse width) produces an output angle θ that follows c with a lag. A first-order model with a dead time, θ/c = K e^{−Ls}/(τs + 1), captures most of it: K is the ratio of the angle reached to the angle commanded (not exactly 1: the servo's pulse-to-angle mapping and its internal dead band see to that), τ is the lag, and L is the latency between the board writing a pulse width and the servo acting on it (up to one 20 ms frame). The outer loop's job is accuracy: it removes the servo's gain error and dead band, which an integrator does. It is type 0, so it needs PI.

In **Option B** the duty u (the fraction of the battery voltage the H-bridge applies, signed for direction) sets the motor's speed through a first-order mechanical-electrical lag τ, and the angle is the integral of the speed: θ/u = K/(s(τs + 1)), K in degrees per second per unit duty. Friction and the driver's threshold mean that below some duty u_dz nothing moves; above it, the speed is K(u − u_dz). That **dead zone** is modelled as a static nonlinearity before the linear plant. The integrator makes the loop type 1: P control already gives zero steady-state error for a step, apart from what the dead zone leaves (2.5).

### 2.2 Identification from a step

**First order with dead time.** For a step of size u at t = 0 from rest,

**y(t) = K u (1 − e^{−(t − L)/τ}) for t ≥ L, and 0 before.**   (1)

The fraction still to go, z = 1 − y/y_ss, obeys **ln z = L/τ − t/τ**: a straight line in t. The **log-linear method**: take y₀ as the mean before the step and y_ss as the mean of the last fifth of the record, K = (y_ss − y₀)/u, fit a line to ln z for 0.1 < z < 0.9, and read τ = −1/slope and L = intercept × τ. The window avoids the start (where z ≈ 1 and the delay distorts the line) and the end (where noise dominates the small z).

**First order plus integrator.** For K/(s(τs + 1)) with dead time L,

**y(t) = K u [(t − L) − τ(1 − e^{−(t − L)/τ})] for t ≥ L**.   (2)

After a few τ the exponential has died and y is a ramp of slope m = Ku that lags the ideal ramp Kut by **T₀ = τ + L**. So fit a straight line to the tail: its slope is Ku and its intercept gives T₀. What is left, r(t) = y − m(t − T₀), is **Kuτ e^{−(t − L)/τ}**, a decaying exponential: a straight line through ln r gives τ, and then L = T₀ − τ.

**Least squares.** The log-linear fits weight the data unevenly (ln amplifies the noise on small values) and use only part of it. They are excellent starting points for a **Levenberg–Marquardt** fit of (1) or (2) to all the samples (`P03.14`): minimise Σ(y_i − f(t_i; p))², with a Gauss–Newton step damped by λ. At the minimum, the linearised covariance s²(JᵀJ)⁻¹, with s² the residual variance on n − 3 degrees of freedom and J the Jacobian, gives each parameter's **standard error**. It is a lower bound on the real uncertainty: it knows the noise, not the model error, which is why the procedure also asks for the spread over repeats.

**Dead zone.** Steps of several sizes give slopes m_i = K(u_i − u_dz): a straight line in u_i, whose slope is K and whose intercept on the u axis is u_dz.

On the SYNTHETIC dataset (Option B, duty steps of 0.30, 0.45 and 0.60 at 100 Hz) the least-squares fits give apparent gains of 360.3 ± 0.4, 420.2 ± 0.3 and 450.2 ± 0.3 deg/s per unit duty, τ = 44.7 ± 1.1, 44.9 ± 0.7 and 45.1 ± 0.6 ms, and L indistinguishable from zero (0.3 ± 0.8 ms at most), with residuals of 0.075–0.101° rms. The apparent gain rises with the step because the dead zone takes a fixed bite out of each duty; the line through the three slopes gives **K = 540.1 deg/s per unit duty and u_dz = 0.100**, and the mean **τ = 44.9 ms**. The log-linear first passes gave τ between 40.2 and 45.6 ms: good starts, visibly noisier.

### 2.3 Designing a P gain for an overshoot target

Under P control with dead-zone compensation (2.5) the loop sees the linear plant, and the closed loop from reference to angle is

**θ/r = KpK/(τs² + s + KpK)**, so **ω_n = √(KpK/τ)**, **ζ = 1/(2√(τKpK))**.   (3)

The overshoot of a standard second-order step is **M_p = exp(−πζ/√(1 − ζ²))**, peak time π/(ω_n√(1 − ζ²)); inverting, ζ = −ln M_p/√(π² + ln²M_p). For 10 %, ζ = 0.5912, and (3) gives **Kp = 1/(4τζ²K)**: with τ = 44.9 ms and K = 540.1, KpK = 15.94 1/s, **Kp = 0.02951 duty/deg**, ω_n = 18.85 rad/s. The firmware gets **Kp = 0.029**, rounded down (a lower gain can only lower the overshoot), for which ζ = 0.5963, ω_n = 18.68 rad/s, poles at −11.14 ± 15.00j, and the predicted step metrics are **overshoot 9.69 %, peak time 0.209 s, rise time (10–90 %) 0.0988 s, 2 % settling time 0.318 s**.

Two remarks. The 10–90 % rise time of a second-order system has no closed form; the code finds the crossings of the exact step response by bisection. And the textbook settling estimate 4/(ζω_n) is an envelope bound: here it gives 0.359 s against the exact 0.318 s, because the response crosses into the 2 % band before its envelope does. Predict with the exact value; quote the estimate only as a sanity check.

Before trusting the gain, check the duty it asks for: a 15° step starts at Kp × 15 + 0.09 = 0.525, inside the 0.6 limit; a 60° step asks for 1.83, so the duty saturates at 0.6 and the motor runs at K(0.6 − u_dz) = 270 deg/s until the error is small. The linear prediction cannot hold for that step.

### 2.4 Designing a PI gain for a rise-time target

For Option A, put the PI zero on the plant pole: C = Kp + Ki/s = Kp(s + Ki/Kp)/s with **Ki/Kp = 1/τ**. The loop becomes KKp/(τs), an integrator, and the closed loop is first order:

**θ/r = 1/(τ_cl s + 1), τ_cl = τ/(KKp)**.   (4)

Its rise time is τ_cl ln 9 = 2.197 τ_cl, its 2 % settling time τ_cl ln 50 = 3.912 τ_cl, and it never overshoots. On the SYNTHETIC servo (command step of 20° at 50 Hz) the least-squares fit gives **K = 0.9282 ± 0.0016, τ = 55.5 ± 0.8 ms, L = 12.6 ± 0.5 ms**. For a rise time of 0.4 s, τ_cl = 182.0 ms, **Kp = 0.3287, Ki = 5.918 1/s**, rounded to 0.33 and 5.9 for the firmware. With the rounded gains the zero sits at −17.88 rad/s and the plant pole at −18.0, so the closed loop has poles at −18.06 and −5.46 rad/s: almost exactly first order. Predicted: **rise 0.402 s, overshoot 0, settling 0.716 s**.

### 2.5 What the continuous model leaves out

**The loop delay.** The firmware writes the output computed at the previous sample, then reads the sensor (a fixed one-sample delay: deterministic, whatever the computation takes); the output is then held for a sample (the zero-order hold of `P09.17`, which behaves on average like half a sample of delay); and the plant may add its own dead time L. In total about **θ = 1.5 Ts + L**: 15.1 ms for Option B at 100 Hz, 42.6 ms for Option A at 50 Hz. A delay θ subtracts ω_c θ radians of phase at the crossover frequency ω_c. For the motor loop, ω_c = 13.4 rad/s and the phase margin falls from 58.9° to 47.3°: the overshoot will be larger than predicted. A delay is not a rational transfer function, but for ω θ ≪ 1, **e^{−θs} ≈ 1/(θs + 1)**, a lag the code can add to the model.

For Option A the same lag turns the integrating loop ω_c/s into ω_c/(s(θs + 1)), which is (3) again: **ζ = 1/(2√(θω_c))**, ω_n = √(ω_c/θ). With ω_c = KKp/τ = 5.51 rad/s and θ = 42.6 ms, θω_c = 0.235, so ζ = 1.03 and ω_n = 11.4 rad/s: a nearly **critically damped** second-order loop. Counter-intuitively its 10–90 % rise is *shorter* than the first-order loop's: the critically damped curve starts slowly but climbs steeply (its 10–90 % rise is 3.358/ω_n against ln 9/ω_c for the first order), and here that gives 0.309 s against 0.398 s. The delay makes this loop look faster and less smooth, not slower.

**The dead zone.** Without compensation, a P loop stops wherever Kp|e| < u_dz: a steady-state error band of ±u_dz/Kp, ±3.4° for Kp = 0.029. Adding a duty c sign(e) in the direction of the error cancels most of it; the firmware uses **c = 0.09** for a fitted 0.100, deliberately low, because over-compensation makes the output chatter around the target. The residual band is (u_dz − c)/Kp = 0.34°, and on the synthetic bench's negative side, where the true dead zone is 0.11, 0.69°.

**Quantisation and noise.** The 12-bit potentiometer reading resolves 0.0732° per count; the 2 % band of a 15° step is 0.3°, only 4 counts wide, and the wiper noise (0.08° rms on the synthetic bench) makes a sample leave the band now and then. Settling times on real data are therefore noisier than rise times.

**Saturation and slew.** The duty limit caps the motor's speed (270 deg/s here); the servo has its own internal slew limit (450 deg/s on the synthetic bench): a 40° open-loop command step asks it for 0.95 × 40/0.060 = 633 deg/s, so the 40° step identifies a slower, more delayed servo (τ = 58.9 ± 0.6 ms, L = 14.8 ± 0.3 ms) than the 20° step, which asks for 317 deg/s. Identify on small steps; check the large ones.

### 2.6 Where the model fails

In order of how soon they bite at home: the battery voltage (K is proportional to it; a sagging pack slows everything), stiction that is larger than the running friction (the motor sticks a little further from the target than the dead-zone model says), gearbox backlash (the pot on the output sees a small free play the model does not have, and a PI loop can hunt across it), a servo whose internal loop is not first order (digital servos especially), and the potentiometer's own nonlinearity near its ends.

## 3. Build the computational model

Open `starter/main.py`. The firmware control law `pi_update` (clamping anti-windup, dead-zone compensation; identical to the listing in `procedure.md`), the state-space helper `tf_to_ss`, the RK4 step map `rk4_map`, `quantize`, `dead_zone`, the SYNTHETIC bench, the CSV writer and reader, `analyse_dataset` and `main` are given.

### Step 1: the models and identification

`first_order_step_response` and `integrator_step_response` from (1) and (2). `fit_first_order_loglinear` (2.2: pre-step mean, tail mean, the ln z line in its window) and `fit_integrator_loglinear` (the tail line, then the ln r line before the tail). `refine_least_squares(model, t, y, p0, lower)`: Levenberg–Marquardt with a central-difference Jacobian, λ halved on success and multiplied by 10 on failure, optional lower bounds by clipping, and standard errors from s²(JᵀJ)⁻¹. `fit_dead_zone` by a straight line.

### Step 2: metrics

`step_metrics(t, y, ...)`: initial value from the pre-step mean, final value from the tail mean (or given), the 10 % and 90 % crossings interpolated linearly between samples, the overshoot from the maximum, the peak time, and the 2 % settling time as the last exit from the band (interpolated), nan if the record ends outside it. `second_order_step` (the closed form for every ζ), `second_order_metrics` (exact overshoot and peak time, exact rise and settling by bisection, and the 4/(ζω_n) estimate), `zeta_from_overshoot`, `overshoot_from_zeta`.

### Step 3: design and prediction

`design_p_position` from (3), `design_pi_first_order` from (4). `closed_loop_tf(num, den, kp, ki)` by polynomial algebra (C G/(1 + C G)); `predict_step` integrates its state-space form with the RK4 map at 0.1 ms; `predict_closed_loop` returns the response, the poles and the metrics, with the final value from the final-value theorem.

### Step 4: the sampled loop

`simulate_closed_loop(num, den, kp, ki, r, Ts, t_end, u_min, u_max, delay_samples, plant_delay, dz, dz_comp, quant, noise_sd, ...)`: every Ts, read the quantised noisy output, call `pi_update`, apply the value computed `delay_samples` earlier through a zero-order hold, delay it further by `plant_delay`, pass it through the dead zone, and integrate the plant with RK4 substeps. `simulate_open_loop` is the same without the controller. `design_p_sampled` finds by bisection the gain that gives an overshoot target on the sampled loop. `compare_metrics` puts two metric sets side by side.

Running `main()` writes the four SYNTHETIC CSV files and `servo_feedback_lab.png` and prints (excerpt):

```
SYNTHETIC DATASET (seeded simulation of the bench, not a measurement)
OPTION B: DC GEAR MOTOR, duty -> angle, model K e^{-Ls}/(s(tau s + 1)) with a dead zone
  step duty 0.30: log-linear K_app 358.8 deg/s, tau 40.2 ms, L 3.2 ms | LM K_app 360.3 +- 0.4, tau 44.7 +- 1.1 ms, L 0.3 +- 0.8 ms, rms 0.075 deg
  dead-zone fit: K = 540.1 deg/s per unit duty, dead zone 0.0999; mean tau 44.88 ms, mean L 0.09 ms
  P design for 10 % overshoot: zeta 0.5912, wn 18.85 rad/s, Kp 0.02951 duty/deg -> firmware Kp 0.029, dead-zone compensation 0.09
  redesign on the sampled model for 10 % overshoot: Kp 0.02265 -> firmware Kp 0.022 (run 3)
  step 15 deg, Kp 0.029: predicted (continuous linear model) | SYNTHETIC | difference
    rise time 10-90 % (s)         0.0988    0.0868   -0.0119
    overshoot (%)                 9.6919   19.8249  +10.1330
    settling 2 % (s)              0.3178    0.3959   +0.0782
OPTION A: HOBBY SERVO, command -> wiper angle, model K e^{-Ls}/(tau s + 1)
  step 20 deg: log-linear K 0.9228, tau 55.6 ms, L 12.7 ms | LM K 0.9282 +- 0.0016, tau 55.5 +- 0.8 ms, L 12.6 +- 0.5 ms, rms 0.130 deg
  PI design for rise time 0.4 s: tau_cl 182.0 ms, Kp 0.3287, Ki 5.918 1/s -> firmware Kp 0.33, Ki 5.9
```

### Verify

Run `python -m pytest curriculum/phase-09/lesson-21-home-lab-servo-feedback`. The 30 tests check the step models (nothing before the dead time, 63.2 % one τ later, the final value, the tail slope and its τ + L lag, and both against `scipy.signal.step`); identification recovering K, τ and L to 10⁻⁶ from noiseless records with an offset, least squares converging from a start 30–40 % off, and, on a seeded record with 0.08° noise and 12-bit quantisation, estimates within 4 standard errors of the truth with residuals at the expected √(σ² + q²/12); the dead zone exactly; the metrics of a densely sampled second-order step against exp(−πζ/√(1 − ζ²)) and π/ω_d and the bisection rise and settling times, and the first-order τ ln 9 and τ ln 50; the P design hitting 5 % and 20 % overshoot with poles at |p| = ω_n and −Re p/|p| = ζ; the PI design hitting its rise time with the zero on the pole; the closed loop against the oracle and against C G/(1 + C G) at three frequencies, and the type-0 steady-state error 1/(1 + KKp); the sampled loop at 0.5 ms with no delay reproducing the continuous prediction within 2–3 %; the duty limit never exceeded and actually reached; seeded noise reproducible, quantised output on the count grid, no output at zero gain; the one-sample delay arriving exactly one sample late and adding more than 5 points of overshoot; `design_p_sampled` hitting 10 % within 0.05 points with a gain below the continuous design; and the SYNTHETIC files.

### Compare

The oracle is `scipy.signal`: `lti` and `step` reproduce the predicted closed-loop step to 10⁻⁹ and the closed-form step models of the identified plants. The lesson's own independent check is the experiment: the identified model, run through the sampled-loop simulator with the firmware's Ts, delay, limit, dead zone and quantisation, must predict the closed-loop rise time, overshoot and settling time that were measured (or, without a kit, that the SYNTHETIC bench produced). That test is in `test_model_predicted_metrics_against_synthetic_closed_loop`; on your own data it is part 5.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box, the full procedure and the firmware listing are in `projects/p09-servo-control/` (`bom.yaml`, `safety.md`, `procedure.md`). Read the safety box first, every session. **Without the kit:** run the analysis on `projects/p09-servo-control/data/synthetic/` (SYNTHETIC data, written by `write_synthetic_csvs(folder, seed=2026)`), or generate your own with your seed; write your dated prediction from the synthetic open-loop files before you open the closed-loop ones; label every table, plot and report "synthetic", and do not describe any of it as measured.

> **Safety (summary; the authoritative box is `projects/p09-servo-control/safety.md`).** Hazards: a loop with the wrong sign or too much gain driving the output into a stall (heat in the motor, driver, wires and cells); pinched fingertips and caught hair at a moving shaft; anything attached to the shaft being flung; short-circuited packs and USB ports; a board browned out by a motor sharing its supply. Energy limits, carried over from `p04-electronics-bench` plus four more: **battery or USB power only, never mains; motor supply ≤ 6 V** (measured; 4 × NiMH or 3 × alkaline, never 4 fresh alkaline cells), board from USB only; one motor or micro servo, **stall current ≤ 500 mA**, with a 0.5 A resettable fuse in the motor pack's + lead; **nothing on the shaft except a light pointer arm (≤ 5 g, ≤ 80 mm): no propellers, no fans, no blades, no wheels**; **fingers, hair and sleeves clear while the loop is running**; **a software duty limit (|u| ≤ 0.6) or command limit (±60°), a travel limit and a latched stop button in the firmware, tested before the first closed-loop run** (stop test, limit test, sign test, travel-limit test, logged and dated). Every p04 limit still holds (no capacitors above 25 V or 1000 µF, no lithium cells or charging circuits, no soldering required). Controls: build with the pack off and USB unplugged; common ground; mount the actuator; first closed loop at half gain and a 10° step; one hand on the stop button; touch-test motor, driver and fuse after each run. Supervision: none for adults; children with an adult who has read the box. Stop criteria: sustained oscillation, an end stop or travel limit hit, a stall hum, any warm part or tripped fuse, a smell of burning, a pack above 6.0 V, a board reset during motion, a failed safety test, anything on the shaft but the pointer. Not permitted: mains or bench supplies, motors above 500 mA stall or 6 V, brushless motors and ESCs, any load on the shaft, removing or raising the limits, closing the loop before the four tests pass, an unattended powered loop.

The bench (Option B shown; Option A replaces driver, motor and pot by the servo and its feedback wire):

```
  USB ── board ── PWM ×2 ──► H-bridge ──► gear motor ══ shaft ══ tube ══ 10 kΩ pot ── ADC pin
          │  └── EN ───────►    ▲                │                          │
          │  └── button ── GND  │ VM             └ pointer arm (card) over a printed protractor
          │                     │
          └──── GND ────────────┴──── AA pack (≤ 6.0 V) ── 0.5 A PTC fuse ── (+)
```

Procedure in outline (full version in `procedure.md`):

1. **Build** with the power off; check continuity; measure the pack (≤ 6.0 V).
2. **Program the board** with the listed firmware (`KP = KI = 0`) and run the **four safety tests**; log them, dated.
3. **Calibrate** the angle at −45° and +45° on the protractor; record the noise at rest and the ADC's true resolution.
4. **Predict the open-loop step** roughly from the actuator's listing (dated), then record **open-loop steps** (B: duties 0.30, 0.45, 0.60 for 0.5 s; A: commands 20° and 40°), three repeats each.
5. **Identify** K, τ, L (and u_dz) with standard errors and the repeat spread.
6. **Design and predict, dated, before closing the loop**: the gain(s), and for each planned step (B: 15° and 60°; A: 20° and 45°) the rise time, overshoot, settling time, steady-state error and peak duty or command.
7. **Close the loop**: half gain and 10° first; then the designed gain, each step three times.
8. **Analyse, explain, refine, redesign, re-predict (dated), re-measure.**

Raw data go in `data/raw/` in the formats of `data/README.md`, never edited afterwards.

## 5. Validate and improve

**The first comparison.** On the SYNTHETIC Option B data, the dated continuous prediction for the 15° step was 9.69 % overshoot; the synthetic bench gave **19.82 %**, double. The rise time was 12 % shorter than predicted (0.0868 against 0.0988 s) and the settling time 25 % longer (0.396 against 0.318 s). Far outside any noise: the model is missing something, and the comparison table says what it looks like: less damping, as from a delay.

**Refine the model, step by step.** Three predictions of the same run, each including more of the loop:

| Option B, 15° step, Kp = 0.029 | rise 10–90 % (s) | overshoot (%) | settling 2 % (s) |
|---|---|---|---|
| continuous model (2.3), the dated prediction | 0.0988 | 9.69 | 0.318 |
| + loop delay as a lag, θ = 1.5 Ts + L = 15.1 ms | 0.0901 | 20.12 | 0.482 |
| sampled model: Ts, 1-sample delay, limit, dead zone, compensation, quantisation | 0.0850 | 21.89 | 0.360 |
| SYNTHETIC bench | 0.0868 | 19.82 | 0.396 |

The lag alone moves the overshoot from 9.7 % to 20.1 %: the delay is the effect. The sampled model then gets the rise time within 3 % and the settling time within 10 %. For the 60° step the continuous model is useless (rise 0.0988 s predicted, 0.1918 s on the bench: the duty sat at its 0.6 limit), and only the sampled model, which has the limit, gets it: rise 0.1921 s, overshoot 8.94 % against 9.19 %, settling 0.461 against 0.479 s. The logged peak duty was 0.600: the limit acted, and was never exceeded.

**Is the remaining difference real?** One run is one sample. Over the ten seeds 15–24, the synthetic 15° overshoot is 20.8 ± 0.3 % (standard deviation over runs) while the sampled model, identified afresh each time, predicts 22.0 ± 0.4 %. The 1.2-point difference is several times the scatter, so it is systematic. The candidate: the bench's dead zone is asymmetric (0.10 forward, 0.11 back), and the overshoot is corrected on the back side, where the extra friction damps it; the model has a symmetric 0.100. Falsification test: run the synthetic bench alone for ten 15° steps (seeds 15–24), first as it is (20.9 ± 0.5 % overshoot), then with a symmetric 0.10 dead zone: the overshoot rises to 21.8 ± 0.4 %, onto the model's prediction. On hardware the equivalent test is a step test in each direction.

**Redesign.** `design_p_sampled` finds the gain that gives 10 % on the sampled model: Kp = 0.02265, rounded to 0.022 for run 3, with a new dated prediction. The synthetic bench gave 6.27 % overshoot (target met), rise 0.1255 s against 0.1168 predicted, settling 0.451 against 0.408 s. The price is visible too: a steady-state error of −0.78°, because the uncompensated part of the dead zone now acts through a smaller gain (the band is 0.45° forward and 0.91° back at Kp = 0.022, against 0.34° and 0.69° at 0.029). A designer would now choose between a slightly higher gain, a better-matched compensation, or a small integral term with anti-windup (`P09.11`).

**Option A.** The dated PI prediction was a first-order loop with rise 0.402 s and settling 0.716 s. The synthetic bench gave **0.2645 s and 0.515 s** with 0.55 % overshoot: faster, not slower. The hand refinement of 2.5 (ζ = 1.03, ω_n = 11.4 rad/s) predicts 0.309 s and 0.544 s; the lag model in the code 0.312 s and 0.551 s; the sampled model with Ts = 20 ms, the one-sample delay and L = 12.6 ms, 0.279 s and 0.499 s, within 6 % and 4 %. The residual 0.55 % overshoot and the small steady-state error (0.046°) are the servo's internal dead band and slew limit, which no model here contains.

| check | hand (part 2) | code | independent |
|---|---|---|---|
| ζ for a given overshoot | 0.5912 for 10 % | `zeta_from_overshoot` 0.5912 | the predicted closed loop's overshoot equals the target to 10⁻³ points (test, 5 % and 20 %) |
| Kp for 10 % (Option B) | 1/(4τζ²K) = 0.02951 | `design_p_position` 0.02951 | closed-loop poles: |p| = ω_n, −Re p/|p| = ζ (test) |
| predicted step, Kp = 0.029 | M_p 9.69 %, t_p 0.209 s | `predict_closed_loop` 9.69 %, 0.2095 s | `scipy.signal.step` to 10⁻⁹ (test) |
| delay effect, Option B | PM 58.9° → 47.3° | lag model 20.1 %, sampled 21.9 % | SYNTHETIC bench 19.8 % (20.8 ± 0.3 % over 10 seeds) |
| saturation, 60° step | duty 1.83 → 0.6, 270 deg/s | sampled rise 0.192 s | SYNTHETIC 0.192 s; peak duty 0.600 |
| PI by cancellation (Option A) | τ_cl ln 9 = 0.400 s | 0.402 s with rounded gains | test: rise = target to 10⁻⁴ |
| delay effect, Option A | ζ = 1.03 → rise 0.309 s | sampled 0.279 s | SYNTHETIC 0.265 s |

## 6. Defend your work

Submit: the logged safety tests, `data/calibration.csv`, both dated predictions, the raw CSV files, your analysis script (using this lesson's functions), the identification table with standard errors and repeat spread, the comparison tables for the continuous, lag and sampled models, the redesign and its measurement, `servo_feedback_lab.png`, the ranked discrepancy causes with the test of the first one, the completed safety box, and a report following `projects/p09-servo-control/report-template.md`. If you used the synthetic dataset, the title, the summary and every figure say SYNTHETIC.

### Adversarial self-checks

- **What if the gain doubles?** On the continuous model ζ falls by √2 (0.596 → 0.42) and the overshoot rises from 9.7 % to 23.2 %; on the sampled model, with the one-sample delay, it rises to 42.6 %. The linear sampled loop loses stability at Kp = 0.131, 4.5 times the design (13 times without the one-sample delay): the delay, not the plant, sets the gain margin. The half-gain first run of the procedure exists for this reason.
- **What if the step goes to zero?** A step smaller than about (u_dz − c)/Kp produces no motion at all; the metrics are then undefined, and `step_metrics` refuses a zero change. Small-step behaviour is the dead zone's, not the linear model's.
- **Which assumption fails first in the real world?** The battery voltage: K is proportional to it, so a pack that has sagged 10 % since the step test lowers KpK by 10 % and changes ζ by 5 %. Re-measure the pack, or repeat one open-loop step, at the start of every closed-loop session.
- **What measurement could falsify the model?** Steps in both directions (an asymmetric dead zone or a gravity load shows as different overshoots); a step at twice the loop rate (a delay-dominated discrepancy must shrink); and the pointer filmed at 240 frames per second (a potentiometer that lags or is nonlinear would disagree with the video).
- **Does the result survive a change in method?** The log-linear and least-squares τ differ by up to 4.5 ms on single runs and by 1.5 ms on the mean of three; designing from the log-linear fits alone gives Kp = 0.0305 instead of 0.0295, 3.2 % higher, which the rounding down to 0.029 absorbs. The conclusion (the delay doubles the overshoot) is independent of the fitting method.
- **Which decision changes if the uncertainty doubles?** With a repeat spread of 0.7 points instead of 0.3, the 1.2-point residual on the 15° overshoot is no longer clearly systematic, and you could not claim to have found the asymmetric friction; the decision to redesign for the delay (a 10-point effect) stands.

## Check yourself

1. Write the step response of K e^{−Ls}/(τs + 1) and show which quantity is linear in t. What does the log-linear method read from the slope and the intercept?
2. For K/(s(τs + 1)), what are the slope and the time-axis intercept of the response's straight tail, and how do you then separate τ from L?
3. Give ζ and ω_n of the P loop on K/(s(τs + 1)), and the gain that gives 10 % overshoot when τ = 45 ms and K = 540 deg/s per unit duty.
4. Why does PI with Ki/Kp = 1/τ on K/(τs + 1) give a first-order closed loop, and what are its rise and settling times?
5. A loop runs at 100 Hz and writes the output computed one sample earlier. About how much delay does the loop see, and what does it do to the phase margin at a crossover of 13.4 rad/s?
6. A P loop on a gear motor with a dead zone of 0.10 and compensation of 0.09 uses Kp = 0.029 duty/deg. What steady-state error band do you expect, and why not compensate with 0.11?
7. Name the four firmware safety tests, the order they run in, and the hazard each one prevents.
8. Your SYNTHETIC analysis shows 19.8 % overshoot where the dated prediction said 9.7 %. What may you write in the report, and what may you not?
