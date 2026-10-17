# Author brief: Phase 09 lessons

This brief is what every lesson author in Phase 09 works from. Read it fully, then read `LESSON_TEMPLATE.md`, `CONTRIBUTING.md`, `references/phase-09.md`, `references/sources.yaml`, the Phase 09 table in `ROADMAP.md`, and the `curriculum/_template/lesson-template/`, `solutions/_template/` and `projects/_template/` folders. **Phases 00–05 are the quality bar:** before writing, read `curriculum/phase-02/lesson-05-matrix-exponential/` and `curriculum/phase-04/lesson-10-ac-circuits-electronics-basics/` (lesson, starter, tests, solution, hints, method, worked solutions) end to end and match their depth, tone and structure. For the home lab (lesson 21), also read `curriculum/phase-04/lesson-15-electronics-bench-lab/` and every file in `projects/p04-electronics-bench/`; that is the bar for the experiment, the bill of materials and the safety box.

## Where files go

```
curriculum/phase-09/lesson-NN-<slug>/lesson.md
curriculum/phase-09/lesson-NN-<slug>/problems.md
curriculum/phase-09/lesson-NN-<slug>/starter/main.py
curriculum/phase-09/lesson-NN-<slug>/tests/test_main.py
solutions/phase-09/lesson-NN-<slug>/hints.md
solutions/phase-09/lesson-NN-<slug>/method.md
solutions/phase-09/lesson-NN-<slug>/worked-solutions.md
solutions/phase-09/lesson-NN-<slug>/solution.py          <- complete, passing implementation of starter/main.py
projects/p09-servo-control/bom.yaml, safety.md, procedure.md, report-template.md, data/README.md   <- lesson 21 only
```

Slugs and builds (fixed, use exactly these):

| # | slug | title | build | mode of part 4 |
|---|---|---|---|---|
| 01 | signals-and-systems | Signals & Systems | Signal toolkit (elementary signals, energy and power, even/odd parts, time shift and scale; tests for linearity, time invariance, causality, memory and BIBO stability on sampled operators) | computational experiment |
| 02 | convolution-impulse-response | Convolution & Impulse Response | Discrete convolution from scratch (direct sum, checked against FFT multiplication; continuous convolution by Riemann sum against the closed form; impulse and step response of an RC IMU front end and a first-order servo) | computational experiment |
| 03 | transfer-functions-block-diagrams | Transfer Functions & Block Diagrams | Block diagram reducer (polynomial transfer-function algebra: series, parallel, feedback; poles, zeros, DC gain; a pitch-attitude loop reduced by hand and by code, checked by simulating the unreduced diagram) | computational experiment |
| 04 | time-response-poles | Time Response & Poles | Step response analyzer (rise time, peak time, overshoot, settling time, steady-state value from a sampled response; second-order ζ–ω_n formulas; dominant poles; effect of a zero; partial-fraction step response against simulation) | computational experiment |
| 05 | frequency-response-bode | Frequency Response & Bode Plots | Bode plotter from scratch (exact magnitude and phase from G(jω), asymptotic straight-line Bode, unwrapped phase; checked against `scipy.signal.bode` and against the steady-state amplitude and phase of a simulated sine response) | computational experiment |
| 06 | sampling-aliasing | Sampling & Aliasing | Anti-aliasing filter design (aliased frequency calculator, Butterworth analog prototype from scratch, order selection for an IMU vibration line, reconstruction by sinc interpolation) | computational experiment on a synthetic IMU record |
| 07 | digital-filters | Digital Filters | FIR and IIR from scratch (windowed-sinc FIR, moving average, bilinear-transform IIR with prewarping, direct-form difference equation; checked against `scipy.signal.freqz` and against the impulse response) | computational experiment on a synthetic gyro record |
| 08 | feedback-steady-state-error | Feedback & Steady-State Error | Cruise controller (airspeed hold: closed-loop TF, system type, error constants Kp, Kv, Ka, steady-state error to step and ramp, disturbance rejection, sensitivity S and T = 1 − S) | computational experiment |
| 09 | root-locus | Root Locus | Root locus plotter (closed-loop roots against gain by polynomial roots; asymptotes, centroid, breakaway points, jω crossings by Routh; every plotted point satisfies the angle condition; gain for a damping target on an altitude-hold loop) | computational experiment |
| 10 | nyquist-stability-margins | Nyquist & Stability Margins | Nyquist plotter, gain and phase margins (contour with indentation around imaginary-axis poles, encirclement count by winding number, Z = N + P agreeing with closed-loop poles; gain margin, phase margin, delay margin, checked against the Bode crossings) | computational experiment |
| 11 | pid-design-tuning | PID Design & Tuning | Pitch attitude PID with anti-windup (discrete PID with derivative filter, actuator saturation, back-calculation and clamping anti-windup; Ziegler–Nichols and pole-placement tuning; simulation of a short-period pitch model with an elevator servo) | computational experiment |
| 12 | lead-lag-compensation | Lead-Lag Compensation | Loop shaping a servo (lead design for phase margin at a target crossover, lag for steady-state error; max phase and frequency of a lead; design checked by the margins of the compensated loop and by its closed-loop step) | computational experiment |
| 13 | state-space-models | State-Space Models | State-space from ODEs (companion and modal forms, similarity transforms, TF ↔ state space both ways; checked: the same TF from C(sI − A)⁻¹B + D as from the ODE; simulation by matrix exponential against `scipy.signal.lti`) | computational experiment |
| 14 | controllability-observability | Controllability & Observability | Rank tests (controllability and observability matrices, PBH eigenvector test, Gramians by Lyapunov solves; checked: rank tests and PBH agree; minimum-energy steering from the Gramian; a reaction-wheel satellite with a failed wheel or sensor) | computational experiment |
| 15 | pole-placement-observers | Pole Placement & Observers | Full-state feedback + Luenberger observer (Ackermann's formula from scratch, checked against `scipy.signal.place_poles`; observer by duality; separation principle; observer error decays at the placed poles) | computational experiment |
| 16 | lqr | LQR | Riccati solver, LQR design (continuous algebraic Riccati equation solved from scratch by the Hamiltonian-matrix eigenvectors and by Newton–Kleinman, checked against `scipy.linalg.solve_continuous_are`; guaranteed margins; Bryson's rule; satellite attitude) | computational experiment |
| 17 | digital-control | Digital Control | Discretization, z-transform (ZOH by the augmented matrix exponential, Tustin, matched pole-zero; checked against `scipy.signal.cont2discrete` and against the matrix exponential; z-plane stability; sample-rate and delay effects on a pitch loop) | computational experiment |
| 18 | robustness-uncertainty | Robustness & Uncertainty ★ | Sensitivity functions, H∞ intro (S, T, KS, peaks Ms and Mt, disk margin from Ms, multiplicative uncertainty and the small-gain test ‖W T‖∞ < 1, H∞ norm by frequency sweep and by Hamiltonian bisection; Bode sensitivity integral checked numerically) | computational experiment |
| 19 | nonlinear-control-basics | Nonlinear Control Basics | Lyapunov, feedback linearization (Lyapunov equation for linear systems, region-of-attraction estimate for a pendulum, describing function of saturation and relay with limit-cycle prediction, feedback linearization of an inverted pendulum or a nonlinear actuator; checked against simulation) | computational experiment |
| 20 | model-predictive-control | Model Predictive Control ★ | Linear MPC from scratch (condensed QP from the prediction matrices, unconstrained solution checked equal to finite-horizon LQR and, with the Riccati terminal cost, to infinite-horizon LQR; box constraints on input and rate by a from-scratch projected-gradient or active-set QP; receding horizon on an altitude-hold model) | computational experiment |
| 21 | home-lab-servo-feedback | Home Lab: Servo Feedback Control | Close the loop on a servo or small motor; identify a model from a step test; design a P or PI gain; measure rise time, overshoot and settling and compare with the model | physical experiment, `projects/p09-servo-control` |

Lesson ids are `P09.NN`. Prerequisites: every lesson needs the previous Phase 09 lessons it uses and the earlier lessons it builds on; name them in the metadata row. Useful earlier lessons: `P01.13` linear systems, `P01.14` rank, `P01.16` eigenvalues, `P02.02` second-order ODEs, `P02.03` forced response and resonance, `P02.04` systems of ODEs, `P02.05` matrix exponential, `P02.06` linearization, `P02.07` complex numbers and phasors, `P02.08` Laplace transform, `P02.09` Fourier series, `P02.10` Fourier transform and the FFT, `P02.19` random processes, `P03.05` RK4, `P03.06` adaptive and stiff integrators, `P03.14` unconstrained optimization, `P03.15` constrained optimization, `P04.10` AC circuits and RC filters, `P04.15` electronics bench lab. Prerequisites in phases not yet written (for example Phase 06 dynamics) may be named where the roadmap names them; otherwise rely on Phases 00–05. The learner has finished Phases 00–05: Python, the scientific stack, calculus, linear algebra, ODEs, Laplace and Fourier transforms, numerical methods, introductory physics and mechanics are known. What is new is systems thinking: signals as objects, LTI systems as operators, feedback, stability, and design.

## What Phase 09 must achieve

Signals, systems and feedback control at the level of a first course in signals and systems and a first and second course in automatic control (MIT 6.003, 16.06 and 16.30/16.31 scope), taught through aerospace control loops. The organising idea of the phase is *model, analyse, design, verify*: every lesson writes the plant model with its assumptions, analyses it in the time domain and the frequency domain, designs something, and **checks every main result by an independent route**. Required independent checks (each lesson's tests must contain its row):

| lesson | main result | independent check |
|---|---|---|
| 01 | system property tests | analytic verdict on known systems (a delay, a squarer, an accumulator, a time-varying gain) |
| 02 | discrete convolution | FFT multiplication with zero padding; continuous closed form |
| 03 | reduced transfer function | time simulation of the unreduced block diagram; `scipy.signal.lti` |
| 04 | step-response metrics | ζ–ω_n closed-form formulas; partial-fraction inverse Laplace |
| 05 | Bode magnitude and phase | `scipy.signal.bode`; amplitude and phase of a simulated sine response at steady state |
| 06 | Butterworth prototype and aliasing | `scipy.signal.butter(analog=True)` as oracle; FFT of a sampled tone lands at the predicted alias |
| 07 | FIR and IIR responses | `scipy.signal.freqz`; the filter's impulse response through the difference equation |
| 08 | steady-state error | final-value theorem against a long simulation |
| 09 | root-locus points | every point satisfies the angle condition ∠KG = 180° and |KG| = 1; jω crossing agrees with Routh |
| 10 | Nyquist encirclements | Z = N + P agrees with the count of closed-loop poles in the right half plane; margins agree with Bode crossings |
| 11 | PID closed loop | linear closed loop without saturation against the transfer-function prediction; anti-windup reduces overshoot under saturation |
| 12 | compensated margins | margins measured on the Bode of the compensated loop equal the design targets; closed-loop step |
| 13 | state-space ↔ TF | C(sI − A)⁻¹B + D evaluated pointwise against the polynomial TF; simulation against `scipy.signal.lti` |
| 14 | rank tests | PBH test agrees with the rank test; Gramian from a Lyapunov solve against quadrature of the integral |
| 15 | pole placement and observer | eig(A − BK) equals the target; `scipy.signal.place_poles`; observer error decays at the placed observer poles (fitted decay rate) |
| 16 | Riccati solution | `scipy.linalg.solve_continuous_are`; residual of the ARE; Hamiltonian and Newton–Kleinman agree |
| 17 | discretization | `scipy.signal.cont2discrete`; the matrix exponential (augmented-matrix ZOH); step invariance at sample instants |
| 18 | H∞ norm and Ms | frequency sweep against Hamiltonian bisection; disk margin bound from Ms against measured gain and phase margins |
| 19 | Lyapunov and describing function | Lyapunov equation residual and `scipy.linalg.solve_continuous_lyapunov`; predicted limit-cycle amplitude and frequency against simulation |
| 20 | MPC | unconstrained MPC gain equals finite-horizon LQR (Riccati recursion); with terminal cost P from the ARE it equals infinite-horizon LQR; constraints never violated |
| 21 | identified model | model-predicted closed-loop rise time, overshoot and settling against measured (or synthetic) values |

A learner who finishes can classify signals and systems, convolve, write and reduce transfer functions, read time and frequency responses, design anti-aliasing and digital filters, analyse feedback for steady-state error, sketch and compute root loci, apply the Nyquist criterion and read margins, design and tune a PID with anti-windup, shape a loop with lead and lag, build state-space models, test controllability and observability, place poles, design observers, solve the LQR Riccati equation, discretize a controller, quantify robustness, analyse simple nonlinear loops, write a linear MPC, and close a real feedback loop on a hobby servo or small motor at home, with a tested Python module for each.

## Aerospace examples (use these; numbers must be realistic and stated as typical)

- **Pitch attitude.** A short-period approximation of a light aircraft or small UAV, for example Mα, Mq, Zα-based second-order pitch-rate dynamics with an integrator to pitch angle, and a first-order elevator servo. State the numbers as a representative model, not a specific aircraft's data.
- **Servo actuator.** A DC-motor servo, (K/(s(τs + 1))) position plant, used for control surfaces and the home lab.
- **Altitude hold.** Altitude as an integral of climb rate, with an inner pitch loop; an outer loop designed by root locus or MPC.
- **IMU anti-aliasing.** A MEMS gyro or accelerometer sampled at 100–1000 Hz with engine or rotor vibration lines above Nyquist; an RC or Butterworth front end and a digital FIR/IIR after it. All IMU records are **synthetic** (seeded) and labelled synthetic.
- **Satellite reaction-wheel loop.** Rigid-body single-axis attitude J θ̈ = τ with a wheel torque limit; momentum and saturation; LQR and observers.
- Others where natural: cruise (airspeed) hold, roll attitude, an inverted-pendulum launch stand, a thrust-vector gimbal.

## Test harness (mandatory)

Tests import the implementation through the repo's `conftest.py`:

```python
from conftest import load_impl
impl = load_impl(__file__)
```

`pytest curriculum/phase-09/lesson-NN-<slug>` runs against the learner's `starter/main.py`; `AEFS_IMPL=solution pytest curriculum/phase-09/lesson-NN-<slug>` runs against `solutions/phase-09/lesson-NN-<slug>/solution.py`. **Your solution.py must pass your tests**; run the second command from the repo root before you finish, and also run the first and confirm **every failure is `NotImplementedError`** (tests must not fail for any other reason against the starter; given helper functions in the starter must be complete). Run `python3 scripts/check_lessons.py`; it must report no ERROR for your lessons.

- The starter has the same function signatures and docstrings as the solution, with `raise NotImplementedError` bodies for the functions the learner writes, complete given helpers, and a working `main()` that produces the lesson's output when complete.
- Tests assert properties, not one magic number: the independent check in the table above, limiting cases (K → 0, K → ∞, ζ = 0 and ζ = 1, ω → 0 and ω → ∞, a pure delay, a zero input), superposition and time invariance, conservation or residual identities (ARE residual, Lyapunov residual, angle condition), agreement with the SciPy oracle, and dimensional sanity. **12–30 tests per lesson** (parametrized cases count individually). Each lesson's tests run in **under 10 s**.
- **Implement everything from scratch in NumPy.** `scipy.signal` and `scipy.linalg` appear **only as the oracle in tests**: `lti`, `lsim`, `step`, `bode`, `freqs`, `freqz`, `butter`, `cont2discrete`, `place_poles`, `solve_continuous_are`, `solve_discrete_are`, `solve_continuous_lyapunov`, `expm`. The solution may use `numpy.linalg` (eig, solve, roots via `np.roots`, matrix rank, SVD) and `scipy.integrate.solve_ivp` for nonlinear simulation where a lesson says so, but must not call `scipy.signal` or the `scipy.linalg` control routines. **Do not use the python-control package anywhere.**
- No network access in tests. No files written outside `tmp_path`. Matplotlib uses the `Agg` backend (`matplotlib.use("Agg")` before importing pyplot); every figure is closed (`plt.close(fig)`); `main()` takes an output directory argument and saves there, never calls `plt.show()`.
- Anything random is seeded (`numpy.random.default_rng(seed)`), with statistical tolerances justified in a comment.

## Physical experiment (lesson 21)

- Lesson 21 owns `projects/p09-servo-control/` with `bom.yaml` (three tiers: `minimum_cost`, `recommended`, `enhanced`, prices as ranges with the date checked, vendor-neutral search terms, substitutes, `reuse_from: [p04-electronics-bench]`), `safety.md` (the full safety box in the format of `projects/p04-electronics-bench/safety.md`), `procedure.md`, `report-template.md` and `data/README.md` with a seeded **synthetic** dataset in `data/synthetic/` (generated by the solution code; labelled synthetic in its file name, header and every place it appears).
- Apparatus: the bench microcontroller from p04 closes a loop on either (a) a hobby servo with its internal potentiometer wiper brought out (or an external potentiometer on the output), the servo's own loop bypassed or treated as part of the plant, or (b) a small DC gear motor driven by a small motor-driver module, with an encoder or a potentiometer on the output shaft. The learner runs an open-loop step test, identifies a first-order (velocity) or first-order-plus-integrator / second-order (position) model with a stated fit method, designs a P or PI gain to a stated rise-time or overshoot target, **writes the predicted rise time, overshoot and settling time down and dates them before closing the loop**, then measures and compares.
- Hard limits, stated in `safety.md` and repeated in the lesson (carried over from p04-electronics-bench plus four more): battery or USB power only, never mains; motor supply at most 6 V; nothing attached to the motor shaft except a light pointer arm — no propellers, no blades; fingers kept clear while the loop is running; a software current or duty-cycle limit and a stop button (or stop command) in the firmware, tested before the first closed-loop run.
- The learner without a kit completes the lesson with the supplied synthetic dataset, labelled synthetic everywhere it appears. Simulation is never presented as measurement.

## Content standards

- Follow the six-part arc and the metadata table from `LESSON_TEMPLATE.md` exactly, including `**Time** | unmeasured`.
- Audience: a motivated adult who has passed Phases 00–05. Plain language, short sentences, define each term on first use. SI units (feet and knots only where industry practice demands, and say so).
- Every lesson states a concrete *engineering question* in the first paragraph, derives in part 2, builds in part 3, runs an **experiment** in part 4 (physical in lesson 21; otherwise a computational experiment: a parameter sweep, a sensitivity study, a method comparison, or a synthetic dataset clearly labelled synthetic), validates in part 5 (a table of hand vs code vs independent check), and ends with adversarial self-checks and 5–8 check-yourself questions (answers in `worked-solutions.md`).
- **Every number in `lesson.md` and `worked-solutions.md` must be reproduced by `solution.py`** (its `main()` or a function it exposes) or by a probe script you run; do not write numbers you have not computed.
- References: cite only sources in `references/sources.yaml`, by `@id`, with exact section titles from `references/phase-09.md`. Where the map says GAP or partial GAP, write original instruction. Never cite or require a paid resource; paid books appear only under "Optional textbooks for deeper study" in the reference map. Never paste text from a source; write your own explanation and link.
- Problem sets: 5–7 problems, graded easy → hard, at least one derivation or proof, one unfamiliar scenario, one that uses the lesson's code, one judgment question about model limits. Prefer parameterized problems (a seed or the learner's own numbers).
- Solutions: `hints.md` (one hint per problem), `method.md` (steps, no final numbers), `worked-solutions.md` (assumptions, units, derivation, result, independent check, common mistakes, plus check-yourself answers).
- Code: Python 3.11+, NumPy/SciPy/Matplotlib/SymPy only; type hints; docstrings; `if __name__ == "__main__": main()`.
- Honesty rules: never claim the course grants a degree, credential or equivalence; never present simulation as physical measurement; label synthetic data as synthetic; say what a model assumes and where it fails.
- Length: `lesson.md` roughly 200–450 lines. Depth over padding.
- Order of work: **code and tests first** (solution, then starter, then tests passing), prose second.

## Do not

- Do not edit files outside your assigned lesson folders (and `projects/p09-servo-control/` for the lesson 21 author).
- Do not edit `ROADMAP.md`, `references/`, `conftest.py`, `pyproject.toml`, `requirements.txt`, `scripts/` or templates.
- Do not run git commands.
- Do not install packages beyond `requirements.txt`; never install or import python-control.
- Do not add time estimates, grades, or badges.
- Give helper scripts in the scratchpad unique names prefixed with your lesson numbers (for example `p09_l07_probe.py`).
