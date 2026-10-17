# Phase 09: Signals, Systems & Control

A first course in signals and systems and a first and second course in automatic control, at the level of MIT 6.003, 16.06 and 16.30/16.31, taught through aerospace control loops and organised around one idea: *model, analyse, design, verify*. Every lesson writes the plant model with its assumptions, analyses it in the time domain and the frequency domain, designs something, and checks every main result by an independent route. The checks include convolution against FFT multiplication, a reduced transfer function against a simulation of the unreduced block diagram, a Bode plot against the steady state of a simulated sine test, and an aliasing prediction against the FFT of a sampled tone. Root-locus points must satisfy the angle condition, Nyquist encirclements must agree with the closed-loop poles, and a Riccati solution must agree with `scipy.linalg.solve_continuous_are` and with its own residual. A describing-function prediction is compared with a simulated limit cycle, unconstrained MPC with the Riccati recursion, and in lesson 21 an identified model against your own closed-loop measurement.

By the end you can:
- classify signals and systems, and convolve;
- write and reduce transfer functions, and read time and frequency responses;
- design anti-aliasing and digital filters;
- analyse feedback for steady-state error and disturbance rejection;
- sketch and compute root loci, apply the Nyquist criterion and read gain, phase and delay margins;
- design and tune a PID with anti-windup, and shape a loop with lead and lag;
- build state-space models, test controllability and observability, place poles and design observers;
- solve the LQR Riccati equation, discretise a controller and budget its sample rate and delay;
- quantify robustness with sensitivity peaks and the small-gain test, and analyse simple nonlinear loops;
- write a linear MPC with constraints;
- close a real feedback loop on a hobby servo or a small motor at home.

You get a tested Python module for each of these. The mathematics is not new: Phases 01–03 supplied the linear algebra, the Laplace and Fourier transforms and the numerical methods, and Phase 04 the circuits and the bench microcontroller. What is new is systems thinking: signals as objects, LTI systems as operators, feedback, stability and design. The aerospace loops recur across lessons: a pitch-attitude autopilot, a control-surface servo actuator, an altitude hold, an IMU's anti-aliasing and digital filtering, and a reaction-wheel satellite. Each is a representative model stated as typical, never a specific vehicle's data.

One lesson has a physical experiment (marked 🔬): the servo feedback home lab. A learner without the kit completes it with the synthetic dataset supplied, which is labelled synthetic wherever it appears; a simulation is never presented as a measurement. Every other lesson's part 4 is a computational experiment, such as a parameter sweep, a sensitivity study, a method comparison, or an analysis of a synthetic dataset clearly labelled synthetic (the IMU and gyro records of lessons 06 and 07). Lessons 18 and 20 are marked ★: they go beyond a first control course, and no other lesson lists them as a prerequisite.

| # | Lesson | Build it | Mode |
|---|---|---|---|
| 01 | [Signals & Systems](lesson-01-signals-and-systems/lesson.md) | Signal toolkit (elementary signals, energy and power, even/odd parts, time shift and scale; tests for linearity, time invariance, causality, memory and BIBO stability) | computational |
| 02 | [Convolution & Impulse Response](lesson-02-convolution-impulse-response/lesson.md) | Discrete convolution from scratch (direct sum against FFT; continuous convolution by Riemann sum; an RC IMU front end and a first-order servo) | computational |
| 03 | [Transfer Functions & Block Diagrams](lesson-03-transfer-functions-block-diagrams/lesson.md) | Block diagram reducer (series, parallel, feedback; poles, zeros, DC gain; a pitch-attitude loop by hand and by code) | computational |
| 04 | [Time Response & Poles](lesson-04-time-response-poles/lesson.md) | Step response analyzer (rise, peak, overshoot, settling; ζ–ω_n formulas; dominant poles; the effect of a zero; partial fractions) | computational |
| 05 | [Frequency Response & Bode Plots](lesson-05-frequency-response-bode/lesson.md) | Bode plotter from scratch (exact and asymptotic, unwrapped phase; checked against a simulated sine test) | computational |
| 06 | [Sampling & Aliasing](lesson-06-sampling-aliasing/lesson.md) | Anti-aliasing filter design (alias calculator, Butterworth prototype, order selection for an IMU vibration line, sinc reconstruction) | computational, on a synthetic IMU record |
| 07 | [Digital Filters](lesson-07-digital-filters/lesson.md) | FIR and IIR from scratch (windowed sinc, moving average, prewarped bilinear IIR, difference equations) | computational, on a synthetic gyro record |
| 08 | [Feedback & Steady-State Error](lesson-08-feedback-steady-state-error/lesson.md) | Cruise controller (airspeed hold: system type, error constants, disturbance rejection, S and T) | computational |
| 09 | [Root Locus](lesson-09-root-locus/lesson.md) | Root locus plotter (asymptotes, breakaway, jω crossings by Routh; gain for a damping target on an altitude hold) | computational |
| 10 | [Nyquist & Stability Margins](lesson-10-nyquist-stability-margins/lesson.md) | Nyquist plotter and margins (indented contour, winding number, Z = N + P; gain, phase and delay margins) | computational |
| 11 | [PID Design & Tuning](lesson-11-pid-design-tuning/lesson.md) | Pitch-attitude PID with anti-windup (derivative filter, saturation, clamping and back-calculation; Ziegler–Nichols and pole placement) | computational |
| 12 | [Lead-Lag Compensation](lesson-12-lead-lag-compensation/lesson.md) | Loop shaping a servo (lead for phase margin at a target crossover, lag for steady-state error) | computational |
| 13 | [State-Space Models](lesson-13-state-space-models/lesson.md) | State space from ODEs (companion and modal forms, similarity transforms, TF ↔ state space both ways) | computational |
| 14 | [Controllability & Observability](lesson-14-controllability-observability/lesson.md) | Rank tests (controllability and observability matrices, PBH, Gramians; a reaction-wheel satellite with a failed wheel or sensor) | computational |
| 15 | [Pole Placement & Observers](lesson-15-pole-placement-observers/lesson.md) | Full-state feedback and a Luenberger observer (Ackermann from scratch, duality, separation) | computational |
| 16 | [LQR](lesson-16-lqr/lesson.md) | Riccati solver and LQR design (Hamiltonian eigenvectors and Newton–Kleinman; guaranteed margins; Bryson's rule; satellite attitude) | computational |
| 17 | [Digital Control](lesson-17-digital-control/lesson.md) | Discretisation and the z-transform (ZOH by the augmented matrix exponential, Tustin, matched pole-zero; sample-rate and delay effects on a pitch loop) | computational |
| 18 | [Robustness & Uncertainty](lesson-18-robustness-uncertainty/lesson.md) ★ | Sensitivity functions and an H∞ introduction (Ms, Mt, disk margins, multiplicative uncertainty and the small-gain test, H∞ norm by sweep and by Hamiltonian bisection, the Bode sensitivity integral) | computational |
| 19 | [Nonlinear Control Basics](lesson-19-nonlinear-control-basics/lesson.md) | Lyapunov and feedback linearisation (Lyapunov equation, region-of-attraction estimate, describing functions of saturation and relay with limit-cycle prediction) | computational |
| 20 | [Model Predictive Control](lesson-20-model-predictive-control/lesson.md) ★ | Linear MPC from scratch (condensed QP, equivalence with finite- and infinite-horizon LQR, input and rate constraints by an active-set QP; receding horizon on an altitude hold) | computational |
| 21 | [Home Lab: Servo Feedback Control](lesson-21-home-lab-servo-feedback/lesson.md) 🔬 | Close the loop on a servo or small motor; identify a model from a step test; design a P or PI gain; predict, measure and compare rise time, overshoot and settling | physical experiment ([project](../../projects/p09-servo-control/)) |

Lesson ids are `P09.01` to `P09.21`. Each lesson's metadata table gives its exact prerequisites, and every lesson's metadata table reads:

| **Time** | unmeasured |
|---|---|

The Time row stays `unmeasured` until pilots measure it. No time estimate is given here for the same reason.

## Prerequisites

Phases 00–05 (no later phase is assumed), in particular:
- `P00.13` (safety for home experiments) and `P00.14` (measurement, calibration and uncertainty) for lesson 21;
- `P01.13` (linear systems), `P01.14` (rank) and `P01.16` (eigenvalues) for the state-space lessons;
- `P02.02` (second-order ODEs), `P02.03` (forced response and resonance), `P02.04` (systems of ODEs), `P02.05` (the matrix exponential) and `P02.06` (linearisation);
- `P02.07` (complex numbers and phasors), `P02.08` (the Laplace transform), `P02.09` (Fourier series) and `P02.10` (the Fourier transform and the FFT) for lessons 01–07 and 10, and the describing functions of lesson 19;
- `P02.19` (random processes) for the observer-noise study of lesson 15;
- `P03.05` (RK4) and `P03.06` (adaptive and stiff integrators) for the simulations, and `P03.14` and `P03.15` (unconstrained and constrained optimisation) for lessons 16, 20 and 21;
- `P04.10` (AC circuits and RC filters) for lessons 02 and 06, and `P04.15` (the electronics bench lab, its microcontroller and its safety box) for lesson 21.

Skill checks before starting: you can find the roots of a quadratic and the eigenvalues of a 2 × 2 matrix by hand, take the Laplace transform of e^{−at} and invert a partial-fraction expansion, compute the magnitude and angle of a complex number such as (1 + jω/a)⁻¹, write a second-order ODE as two first-order ones, and propagate an uncertainty through a product.

## The project kit

The kit lives in `projects/p09-servo-control/`. It holds:
- `bom.yaml`: three cost tiers, `minimum_cost`, `recommended` and `enhanced`, with dated price ranges, reusing the microcontroller, breadboard, potentiometer, multimeter and battery holder of `p04-electronics-bench`;
- `safety.md`: the full safety box;
- `procedure.md`, with the firmware listing and the four pre-run tests;
- `report-template.md`;
- `data/README.md`: the raw-data format, with the seeded synthetic dataset in `data/synthetic/`.

Read [`safety.md`](../../projects/p09-servo-control/safety.md) before buying anything and again before every session. The kit's `safety.md` is authoritative. The limits below summarise the hard rules that every version of this kit must keep, and nothing in the kit may relax them. The prediction, with its uncertainty, is written down and dated before the loop is first closed, raw data are never edited, and the report states which numbers were measured, which were computed, and which are synthetic.

- **`p09-servo-control`** ([folder](../../projects/p09-servo-control/), lesson 21). The p04 microcontroller closes a loop on either a hobby micro servo with its feedback wiper brought out to the ADC (Option A, PI control) or a small DC gear motor on a low-voltage motor-driver module with a potentiometer on its shaft (Option B, P control with dead-zone compensation). You run open-loop step tests, identify a first-order or integrator model with a stated fit method, design the gain, date the predicted rise time, overshoot and settling time, then close the loop, measure and compare. Hard limits (every limit of `p04-electronics-bench` still applies):
  - battery or USB power only, never mains;
  - motor (or servo) supply ≤ 6 V, measured before every session (four NiMH or three alkaline cells, never four fresh alkaline cells);
  - one motor or micro servo, stall current ≤ 500 mA, with a 0.5 A resettable fuse in the pack's + lead;
  - nothing on the shaft except a light pointer arm (≤ 5 g, ≤ 80 mm): no propellers, no fans, no blades, no wheels;
  - fingers, hair and loose sleeves clear of the shaft while the loop is running;
  - a software duty limit (|u| ≤ 0.6) or command limit (±60°), a travel limit and a stop button in the firmware, tested (stop, limit, sign and travel-limit tests, logged) before the first closed-loop run, and never removed or raised;
  - the first closed loop of each session at half the designed gain with a small step.

Nothing in this phase involves mains voltage, propellers, flight, lithium-cell charging, or anything attached to a moving shaft beyond the pointer. Every IMU, gyro and servo record outside lesson 21's own measurements is synthetic and seeded. **A report without a completed safety box and the logged pre-run tests is not graded.**

## Running the code

From the repo root, `pytest curriculum/phase-09/lesson-NN-<slug>` tests your `starter/main.py`. `AEFS_IMPL=solution pytest curriculum/phase-09/lesson-NN-<slug>` runs the same tests against the reference solution in `solutions/phase-09/lesson-NN-<slug>/solution.py`. Tests check properties rather than single numbers:
- each lesson's independent check (for example, FFT against direct convolution, simulation of the unreduced diagram, a simulated sine test, the angle condition, Z = N + P against the closed-loop poles, the ARE residual, step invariance at the samples, a describing-function prediction against a simulated limit cycle, unconstrained MPC against the Riccati recursion);
- agreement with the SciPy oracle (`scipy.signal` and `scipy.linalg` routines such as `lti`, `step`, `bode`, `freqz`, `butter`, `cont2discrete`, `place_poles`, `solve_continuous_are`, `solve_discrete_are`, `solve_continuous_lyapunov`, `expm`), which appears only in the tests;
- limiting cases (K → 0 and K → ∞, ζ = 0 and ζ = 1, ω → 0 and ω → ∞, a pure delay, a zero input);
- superposition, time invariance and residual identities (S + T = 1, Riccati and Lyapunov residuals);
- dimensional sanity, and convergence with step size where a numerical method is used.

Everything in the tested code is implemented from scratch in NumPy; `scipy.integrate.solve_ivp` is used for nonlinear simulation where a lesson says so, and the python-control package is not used anywhere. Anything random is seeded. Tests need no network and write files only to a temporary directory. `python scripts/check_lessons.py` checks every lesson's structure.

## The exam and the mini-project

**Phase exam and mini-project:** [assessments/phase-09/](../../assessments/phase-09/). The exam has four parts:
- **Part A, concepts:** signals and systems, transforms and sampling, feedback and robustness, and state space, estimation, optimal and nonlinear control.
- **Part B, analysis by hand:** every lesson from 01 to 20, on a gyro front end, a pitch-attitude autopilot, an elevator servo, a roll loop, an IMU's anti-aliasing filter, a satellite attitude hold, an altitude hold and an unstable airframe.
- **Part C, computation:** a tested unconstrained-MPC gain function checked against the Riccati recursion, and an end-to-end analysis of a servo design with your own lesson functions (margins, sensitivity peak, discretisation with delay, step metrics).
- **Part D, engineering interpretation:** its first question asks you to interpret your own `p09-servo-control` measurement; a seeded synthetic dataset is supplied for learners without a kit. Its second question reviews a synthetic flight-control design sheet for wrong margins, aliasing and disturbance errors, and tests the servo kit's safety limits.

The exam has not been piloted, so no time limit is suggested. Choose one in advance and record it. Results are reported per competency: (1) reasoning, (2) analysis, (3) computation, (4) interpretation. The worked solutions and pass thresholds are in [solutions/phase-09/exam/](../../solutions/phase-09/exam/worked-solutions.md). The mini-project is the servo feedback report from lesson 21, graded with [the rubric](../../assessments/phase-09/grading-rubric.md). Seeded unseen cases for the computational lessons are planned in [parameterized-cases/](../../assessments/phase-09/parameterized-cases/README.md). Outcomes are recorded as *practice completed*, *self-assessed mastery* or *externally reviewed mastery*. None of these is a credential, and nothing in this phase qualifies anyone to design, tune or approve the control software of a real aircraft or spacecraft.

**References:** the Phase 09 reference map (`references/phase-09.md`) has not been written yet, so the lessons cite only sources already in `references/sources.yaml` and write original instruction where no mapped free source covers a topic. Every required reading is free. The signals-and-systems lessons use the MIT 6.003 lecture notes (`@mit-ocw-6-003-signals-and-systems`) and the MIT 18.03 notes (`@mit-ocw-18-03-differential-equations`). The state-space lessons use MIT 18.06 (`@mit-ocw-18-06-linear-algebra`). LQR and robustness use the MIT 16.323 lecture notes (`@mit-ocw-16-323-principles-of-optimal-control`). The SciPy user guide (`@scipy-userguide`) documents the oracle routines. The home lab uses `@micropython-docs` and `@sparkfun-learn-tutorials` for the board's ADC. Each lesson's references row names the exact sections. The curriculum links to every source and writes its own explanations.
