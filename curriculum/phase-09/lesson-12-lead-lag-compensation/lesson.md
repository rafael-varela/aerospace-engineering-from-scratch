# P09.12  Lead-Lag Compensation

> Lead buys phase where the loop crosses over; lag buys gain where it does not. Shape the loop in the frequency domain, then prove the shape by measuring it.

| | |
|---|---|
| **Learning outcomes** | 1. Derive the maximum phase of a lead compensator, sin φ_max = (1 − α)/(1 + α), its frequency 1/(T√α) and its gain 1/√α there, and the phase cost of a lag a decade above its zero. 2. Design a lead for a target crossover frequency and phase margin, and a lag for a velocity error constant, on a control-surface servo G = K/(s(τs + 1)), and implement both from scratch. 3. Implement and test a margin finder (gain crossover, phase margin, phase crossover, gain margin, delay margin) and verify that the margins measured on the compensated loop equal the design targets, independently of the design formulas and against `scipy.signal.bode`. 4. Relate phase margin to overshoot, and judge the costs of lead (high-frequency gain, noise) and lag (slow tail, phase loss). |
| **Prerequisites** | `P09.05` (Bode plots, phase of poles and zeros), `P09.08` (system type, Kv and ramp error), `P09.10` (gain, phase and delay margins), `P09.04` (step metrics, ζ and overshoot), `P09.03` (series and feedback connections). Skill check: you can sketch the asymptotic Bode plot of 20/(s(0.05s + 1)) and read its crossover. |
| **Engineering question** | An elevator servo (DC motor, gearbox, position feedback) is modelled as G = 20/(s(0.05s + 1)). It must follow commands with a crossover of 30 rad/s, a phase margin of 55° and a ramp-following error of no more than 5 mrad per rad/s of command rate (Kv ≥ 200 s⁻¹). Can a gain alone do it, and if not, what lead-lag compensator does? Answer: a gain that meets Kv crosses at 61.7 rad/s with only 18.0° of margin (60.5 % overshoot). A lead with α = 0.3896 centred at 30 rad/s plus a lag with β = 5.953 and its zero at 3 rad/s meets all three targets exactly: crossover 30.000 rad/s, PM 55.000°, Kv = 200.0 s⁻¹, with 16.5 % overshoot. |
| **Assumptions and limits** | Linear time-invariant servo: a DC motor whose electrical time constant is neglected (one mechanical lag τ), a stiff gearbox with no backlash, no friction, no current or voltage limit, a perfect position sensor, continuous-time control. The model is a representative servo, not a specific product. It fails for large steps (the amplifier saturates), for small ones (friction and backlash dominate), and above the first structural or electrical resonance, which the lead amplifies. |
| **Analytical work** | The phase of (Ts + 1)/(αTs + 1) and its maximum; the lead design equations (φ needed, α, T, Kc); the lag's DC gain and its phase at 10z; Kv of the lead-lag loop; the fixed-point iteration for β; the phase-margin to damping relation of the standard second-order loop; the delay margin PM/ω_c. |
| **Computational deliverable** | `starter/main.py`: `freq_response`, `unwrapped_phase_deg`, `margins`, `velocity_constant`, `lead_max_phase`, `lead_alpha`, `lead_compensator`, `lag_compensator`, `design_lead`, `design_lead_lag` (`servo_plant`, `servo_tf`, `series`, `unity_feedback`, `expm`, `tf_response`, `step_metrics`, `main` are given). Tests prove: the frequency response equals `scipy.signal.freqs`; the lead's numerical maximum phase equals asin((1 − α)/(1 + α)) at 1/(T√α) for three α; the lead and lag limits; the margin finder gives GM = 6 at √2 rad/s for 1/(s(s + 1)(s + 2)) and agrees with `scipy.signal.bode`; **the margins measured on every compensated loop equal the design targets** (to 10⁻⁷ in ω_c and 10⁻⁵ degrees in PM, and to 10⁻³ degrees on `scipy.signal.bode`); Kv meets its target; the closed-loop step equals `scipy.signal.step` and its overshoot falls as PM rises; the ramp error from `scipy.signal.lsim` is 1/Kv. |
| **Practical activity** | **computational experiment**: a phase-margin sweep (overshoot and settling against PM), a lag-placement sweep (lag phase at crossover, overshoot and settling tail against z = ω_c/ratio), and a gain-only versus lead versus lead-lag comparison. |
| **Equipment** | None. A computer with the environment from `P00.01`. The servo here is the same kind of plant the home lab `P09.21` identifies on real hardware. |
| **Safety** | None beyond normal computer use. |
| **Evidence** | Completed `main.py`, passing tests, `lead_lag.png`, the derivation of φ_max and of the design equations, the two sweep tables, and a half-page design note giving the compensator, its measured margins, Kv, overshoot and settling time. |
| **Assessment** | `assessments/phase-09/grading-rubric.md`; competencies 1–4. Pass: φ_max and the lead equations are derived; all public tests pass; the note shows the margins measured on the compensated loop, not just the design inputs, and explains the lag's slow tail with numbers. |
| **References** | `@mit-ocw-6-003-signals-and-systems` course lecture notes (required-free, course level: feedback and the Laplace-domain view of LTI systems). `@scipy-userguide` "Linear Algebra (scipy.linalg)" (supplementary-free: the oracle routines used in the tests). The Phase 09 reference map (`references/phase-09.md`) did not exist when this lesson was written, so no section of a control text is mapped yet: the PID, compensation, state-space and observer material here is original instruction (GAP), to be linked to exact sections when the map is published. |
| **Extensions** | ★ Add the motor's electrical time constant (a second pole near 200 rad/s) and a 1 ms computation delay, and redesign for the same targets; how much of the phase margin is left? ★ Replace the lag by a PI (lag with β → ∞) and compare the step tail. ★ Add measurement noise and compare the actuator activity of lead designs with PM 45°, 55° and 65°. |
| **Time** | unmeasured |

## 1. The engineering problem

A control surface is moved by an actuator: on small aircraft and UAVs usually an electric servo, a DC motor driving the surface through a gearbox, with a potentiometer or encoder on the output and a feedback loop around it. The outer loops (`P09.11` pitch, roll, autopilot) treat that servo as "a fast first-order lag". That is true only if the servo loop itself is well designed: fast enough that its lag does not eat the outer loop's phase margin, damped enough that it does not ring, and stiff enough to follow a ramping command without falling behind.

The position plant of such a servo is G(s) = K/(s(τs + 1)): an integrator (voltage sets speed, speed integrates to angle) and the motor's mechanical time constant. A proportional controller has one knob, and the three requirements pull it in different directions. A high gain gives a small ramp error (large Kv) but pushes the crossover up to where the plant's phase is near −180°, so the margin collapses. A low gain gives a good margin but a large ramp error. **Lead-lag compensation** separates the jobs. A lead network adds phase exactly at the chosen crossover. A lag network adds gain at low frequency only, a decade or more below the crossover, where it does not disturb the margin.

Success criteria, written first: with the compensator in the loop, the gain crossover measured on |L(jω)| is 30 rad/s, the phase margin measured there is 55°, both to the precision of the margin finder and confirmed by `scipy.signal.bode`; Kv = 200 s⁻¹; the closed-loop step overshoots less than 20 % and the ramp error is 1/Kv = 5 mrad per rad/s.

## 2. Understand and derive

### 2.1 Why a gain cannot do it

With C = Kc, L = Kc K/(s(τs + 1)) and Kv = lim s L = Kc K. The phase is −90° − arctan(ωτ), above −180° at every frequency, so the gain margin is infinite and the phase margin is 180° + arg L(jω_c) = 90° − arctan(ω_c τ). For τ = 0.05 s:

- Kv = 200 needs Kc = 10. The crossover is then where |L| = 1, at **61.7 rad/s**, and the margin is 90° − arctan(3.08) = **18.0°**, with 60.5 % overshoot in the closed-loop step.
- A crossover at 30 rad/s leaves 90° − arctan(1.5) = **33.7°** (the plant's phase there is −146.3°), and needs a much lower gain (Kv = 37 for the margin-55° lead design below).

No gain meets ω_c = 30 rad/s, PM = 55° and Kv = 200 at once.

### 2.2 The lead compensator

C_lead(s) = Kc (Ts + 1)/(αTs + 1), 0 < α < 1: a zero at 1/T below a pole at 1/(αT). Its phase is

φ(ω) = arctan(ωT) − arctan(αωT).

Setting dφ/dω = 0: T/(1 + ω²T²) = αT/(1 + α²ω²T²) ⇒ 1 + α²ω²T² = α + αω²T² ⇒ ω²T²α(1 − α) = 1 − α ⇒

**ω_max = 1/(T√α),** the geometric mean of the zero and the pole.   (1)

There, with x = ω_maxT = 1/√α, tan φ_max = (x − αx)/(1 + αx²) = (1 − α)/(2√α), so

**sin φ_max = (1 − α)/(1 + α),  α = (1 − sin φ)/(1 + sin φ),**   (2)

and |(jω_maxT + 1)/(jω_maxαT + 1)| = √((1 + 1/α)/(1 + α)) = **1/√α**. A lead with α = 0.1 gives 54.9° at 3.16 times its zero frequency and raises the high-frequency gain by 1/α = 10 (20 dB). One stage is practical up to about 60–65°; beyond that α is tiny and the high-frequency gain enormous, and two stages are used (Problem 6).

**Design for crossover ω_c and margin PM** (`design_lead`):

1. φ = PM − (180° + arg G(jω_c)): the phase the lead must add.
2. α from (2); T = 1/(ω_c√α), so the maximum sits at ω_c.
3. Kc = 1/|C_lead,1(jω_c) G(jω_c)|, so that |L(jω_c)| = 1.

For the servo, φ = 55° − 33.69° = **21.31°**, α = **0.4669**, zero 20.50 rad/s, pole 43.90 rad/s, Kc = **1.8478**. This loop has Kv = Kc K = **36.96 s⁻¹**, a ramp error of 27.1 mrad per rad/s. The lead has fixed the margin and left Kv more than five times too small.

### 2.3 The lag compensator

C_lag(s) = (s + z)/(s + z/β), β > 1: DC gain β, high-frequency gain 1, a pole below its zero. Placed a decade below the crossover (z = ω_c/10) it multiplies Kv by β, while at ω_c its gain is ≈ 1 and its phase is arctan(ω_c/z) − arctan(βω_c/z) = arctan 10 − arctan(10β), which lies between −5.7° and 0° for every β. For β = 5.95 it is **−4.75°**.

Because the lag changes the phase at ω_c slightly, the lead is designed on G C_lag, not on G. And because the lead's gain Kc then changes slightly, β must be iterated: start with β = 1, design the lead on G C_lag, compute Kv = Kc β K, update β ← β · Kv,target/Kv, and repeat until β settles (a few iterations; `design_lead_lag`). Result: **β = 5.9534** (lag zero 3.00 rad/s, pole 0.504 rad/s), and the lead needs α = **0.3896** (φ_max = 26.06°, the extra 4.75° covering the lag), Kc = **1.6797**. Kv = 1.6797 × 5.9534 × 20 = **200.0 s⁻¹**.

### 2.4 Phase margin, damping and delay

For the standard loop L = ω_n²/(s(s + 2ζω_n)), whose closed loop is the textbook second-order system (`P09.04`),

PM = arctan[2ζ / √(√(1 + 4ζ⁴) − 2ζ²)],

which is close to PM ≈ 100ζ (in degrees) for ζ ≤ 0.7: ζ = 0.5 gives 51.8°, ζ = 0.55 gives 55.7° and 12.6 % overshoot. The lead loop is not exactly this system (it has a zero and a third pole), but the rule predicts its 13.3 % overshoot well. The phase margin also sets the **delay margin**, the pure time delay that would destabilise the loop: PM (rad)/ω_c = 0.960/30 = **32.0 ms**. Sampling, computation and a slow sensor all eat into it (`P09.17`).

### 2.5 Where the model stops being valid

- **Motor electrical time constant and amplifier bandwidth** add poles typically at 100–1000 rad/s. The lead's pole (43.9 rad/s alone, 48.1 rad/s in the lead-lag) and its high-frequency gain Kc/α ≈ 4.3 bring them closer to mattering; extension 1.
- **Current and voltage limits.** For a large step, the linear loop asks the motor for more current than the amplifier supplies, the response becomes rate-limited, and overshoot and settling change.
- **Friction and backlash.** For small commands, static friction makes the servo stick and the lag's slow tail turns into a limit cycle (a small sustained oscillation, `P09.19`).
- **Structural modes.** The surface and its linkage have a first bending or torsion mode; above it the rigid model is wrong, and a lead that raises the gain there can excite it.

## 3. Build the computational model

Open `starter/main.py`. `servo_plant`, `servo_tf`, `series`, `unity_feedback`, `expm`, `tf_response` (exact response to a held input through the companion realisation) and `step_metrics` are given.

### Step 1: frequency response and margins

`freq_response` is `np.polyval(num, 1j*w) / np.polyval(den, 1j*w)`. `unwrapped_phase_deg` evaluates the phase on a fine log grid from a low starting frequency and unwraps it, so that −190° is reported as −190°, not +170°. `margins` sweeps 20 000 log-spaced frequencies, brackets every sign change of log|L| (gain crossovers) and of the phase relative to −180° (mod 360°), and refines each by bisection on a log scale:

```python
for i in np.where(np.sign(gmag[:-1]) != np.sign(gmag[1:]))[0]:
    wc = bisect(lambda v: np.log(abs(L(v))), w[i], w[i + 1])
    pm = 180.0 + phase(wc)          # wrapped into (-180, 180]
```

It returns the smallest PM and GM over the crossings, GM = ∞ when the phase never reaches −180°, and the delay margin.

### Step 2: lead and lag

`lead_max_phase`, `lead_alpha`, `lead_compensator`, `lag_compensator` are one line each from (1), (2) and §2.3. `design_lead` is the three steps of §2.2, raising `ValueError` when more than 65° are needed. `design_lead_lag` is the fixed-point iteration of §2.3; `velocity_constant` divides num(0) by the coefficient of s in the denominator, after checking there is exactly one pole at the origin.

### Verify

Run `python -m pytest curriculum/phase-09/lesson-12-lead-lag-compensation`. The tests:

1. **Building blocks**: `freq_response` against `scipy.signal.freqs`; for α = 0.1, 0.3 and 0.6 the numerically maximal phase on a 200 001-point grid equals (2) to 10⁻⁶ degrees at (1), with gain 1/√α, and `lead_alpha` inverts `lead_max_phase`; α → 1 and α → 0; the lag's DC and high-frequency gains and its phase arctan 10 − arctan 10β at 10z; Kv of a type-1 plant.
2. **Margins**: GM = 6 at √2 rad/s for 1/(s(s + 1)(s + 2)), |L| = 1 at the reported crossover, the delay margin; agreement with crossovers interpolated from `scipy.signal.bode`; infinite GM for the second-order servo.
3. **Independent check**: for three (ω_c, PM) targets the lead loop's measured crossover and margin equal the targets, by our finder and by `scipy.signal.bode`; for three lag placements the lead-lag loop also meets Kv = 200; a lax Kv needs no lag; 160° of lead raises.
4. **Closed loop**: the step response equals `scipy.signal.step` to 10⁻⁹ for PM = 40°, 50°, 60°, reaches 1, and its overshoot decreases with PM; the ramp error from `scipy.signal.lsim` is 1/Kv for the lead-lag (and five times larger for the lead alone); the gain that meets Kv alone has PM < 25°.

### Compare

`main()` prints the design and measures it:

```
lead: phi_max = 21.31 deg, alpha = 0.4669, zero 20.50 rad/s, pole 43.90 rad/s, Kc = 1.8478; measured wc = 30.0000, PM = 55.0000 deg, GM = inf, Kv = 36.96 1/s
lead-lag: beta = 5.9534, lag zero 3.00, lag pole 0.5039 rad/s; lead alpha = 0.3896, phi_max = 26.06 deg, Kc = 1.6797
          measured wc = 30.0000 rad/s, PM = 55.0000 deg, GM = inf, delay margin 32.00 ms, Kv = 200.00 1/s
  gain only (Kv met)  : overshoot  60.5 %, rise  18.5 ms, settle  366.0 ms, ramp error at 10 s 5.000 mrad per rad/s (1/Kv = 5.000)
  lead                : overshoot  13.3 %, rise  43.5 ms, settle  147.5 ms, ramp error at 10 s 27.060 mrad per rad/s (1/Kv = 27.060)
  lead-lag            : overshoot  16.5 %, rise  41.5 ms, settle  471.0 ms, ramp error at 10 s 5.000 mrad per rad/s (1/Kv = 5.000)
```

The margin finder and `scipy.signal.bode` are two different instruments reading the same curve. The design formulas set the targets, and the measurement confirms them to the finder's precision. The ramp errors come from simulating the closed loop with a ramp command and agree with 1/Kv to the printed digit. The ramp is held at the midpoint of each 0.5 ms step; held at the left end, the staircase would lag the true ramp by half a step and add 0.25 mrad.

## 4. Experiment at home, or investigate real data

**Mode: computational experiment.** Write each prediction down before running it.

1. **Phase-margin sweep.** With the lead alone at ω_c = 30 rad/s, design for PM = 40°, 50°, 60° and 70°; record α, overshoot and 2 % settling time. Predict overshoot from PM ≈ 100ζ and the second-order formula, and predict which PM settles fastest.
2. **Lag placement.** With the lead-lag at PM = 55° and Kv = 200, put the lag zero at ω_c/3, ω_c/5, ω_c/10 and ω_c/20; record the lag's phase at ω_c, the lead α needed to compensate, overshoot and settling time over a 4 s window. Predict the trade between phase stolen at crossover and the length of the slow tail.
3. **Gain only, lead, lead-lag.** Tabulate crossover, PM, Kv, overshoot, rise, settling and ramp error for the three designs.

Raw data format: CSV tables `pm_target_deg, alpha, overshoot_pct, settle_ms`; `lag_ratio, lag_phase_at_wc_deg, lead_alpha, overshoot_pct, settle_ms`; `design, wc_rad_s, pm_deg, kv, overshoot_pct, rise_ms, settle_ms, ramp_error_mrad`, each with a header naming the plant and stating that the values are from simulation of a representative model.

What you should find: overshoot falls from 29.3 % (PM 40°, α = 0.802) through 18.4 % (50°) and 8.5 % (60°) to 0.0 % (70°, α = 0.256), and settling is fastest at PM 60° (137 ms) rather than 70° (176.5 ms), because an overdamped loop creeps. Moving the lag zero closer to the crossover steals more phase (−16.0° at ω_c/3, −9.6° at /5, −4.75° at /10, −2.36° at /20). The lead must then supply more (α from 0.427 down to 0.245) and the overshoot grows from 15.1 % to 20.2 %. Settling times are 316, 383, 471 and 416 ms: not monotonic, because the lag's slow pole-zero pair leaves a tail whose amplitude and time constant both change with z, and the 2 % band catches it at different points. The lag's price is that tail. The settling time of the lead-lag (471 ms) is three times that of the lead alone (147.5 ms) for the same margin.

## 5. Validate and improve

| quantity | design target / hand | code (our finder) | independent check | difference |
|---|---|---|---|---|
| lead φ_max at α = 0.4669 | asin(0.5331/1.4669) = 21.31° | 21.31° | grid maximum (test) | < 10⁻⁶° |
| crossover, lead-lag | 30 rad/s | 30.0000 rad/s | `scipy.signal.bode` (test) | < 10⁻⁵ relative |
| phase margin, lead-lag | 55° | 55.0000° | `scipy.signal.bode` (test) | < 10⁻³° |
| Kv, lead-lag | 200 s⁻¹ | 200.00 s⁻¹ | ramp error 5.000 mrad (simulation), `lsim` (test) | — |
| delay margin | 0.960 rad/30 rad/s = 32.0 ms | 32.00 ms | — | — |
| overshoot at PM 55° (lead) | ζ ≈ 0.55 ⇒ 12.6 % | 13.3 % (simulation) | `scipy.signal.step` (test) | zero and third pole |

Everything the design formulas promised is confirmed by an independent measurement of the loop. The one approximate entry is the overshoot rule, which assumes a pure second-order loop; the 0.7-point difference is the lead's zero. These are checks of the analysis pipeline on a model: the hardware servo of `P09.21` is where the model meets friction, backlash and a current limit.

What to improve: include the motor's electrical pole and the sampling delay before trusting a 32 ms delay margin; consider a PI instead of the lag if the slow tail matters more than the ramp error; check the high-frequency gain Kc/α of the final design (1.6797/0.3896 = 4.31) against sensor noise. What not to improve: the lag's β beyond what Kv needs. More low-frequency gain lengthens the tail and adds nothing the specification asked for.

## 6. Defend your work

Submit: the completed `main.py` with passing tests; the derivations of (1), (2), the lead design steps and Kv of the lead-lag loop; `lead_lag.png`; the three experiment tables; and a half-page design note (template in `projects/_template/report-template.md`, sections 1, 2, 5, 6 and 7) for the servo compensator: the transfer function, the measured margins (not the design inputs), Kv, overshoot, settling time, the high-frequency gain, and the first unmodelled effect you would check on hardware.

### Adversarial self-checks

- **What if τ is 30 % larger (a heavier surface)?** The plant phase at 30 rad/s drops from −146.3° to −152.9° (τ = 0.065 s, Problem 3). The lead must supply 27.9° instead of 21.3°. A fixed compensator designed for τ = 0.05 would cross lower and with less margin. Redesign and re-measure.
- **What if the crossover target is raised to 60 rad/s?** The plant phase is −161.6° and the lead must add 36.6°. That is still one stage, but the high-frequency gain rises, and so does the noise reaching the motor. Where is the motor's electrical pole?
- **Which assumption fails first in the real world?** For small commands, friction and backlash; for large ones, the current limit. A real servo is linear only in between.
- **What measurement could falsify the model?** A frequency sweep of the closed servo with a small sine command: if the measured bandwidth or resonant peak differs from the model's closed loop, τ or K is wrong. The step test of `P09.21` is the cheaper version.
- **Does the result survive a change of tool?** The margins agree between our bisection finder and `scipy.signal.bode` to 10⁻³ degrees; the step agrees with `scipy.signal.step` to 10⁻⁹.
- **Which decision changes if the uncertainty doubles?** Scaling the loop gain of the lead-lag design by 0.8 and 1.2 moves the crossover to 25.0 and 34.5 rad/s and the margin to 58.7° and 51.8°. A ±40 % gain error (0.6 to 1.4) gives 19.5 to 38.7 rad/s and 62.5° to 49.0°. The margin stays acceptable; the bandwidth seen by the outer loop does not. If the outer pitch loop needs a servo of at least 25 rad/s, ±40 % uncertainty in K forces either a higher nominal crossover or an identified K (`P09.21`).

## Check yourself

1. Write the lead compensator and state where its phase is largest and by how much.
2. Why does a gain that meets Kv = 200 on this servo leave so little phase margin?
3. What does a lag compensator do to Kv, and why is its zero placed about a decade below the crossover?
4. A lead with α = 0.25 is centred at 40 rad/s. Give its zero, its pole, its maximum phase and its gain at 40 rad/s.
5. Why must the lead be designed on G·C_lag rather than on G when both are used?
6. State the rule of thumb linking phase margin and damping ratio, and the overshoot it predicts for PM = 55°.
7. What is the delay margin of a loop with PM = 55° at ω_c = 30 rad/s, and what in a digital servo consumes it?
8. Why is "measure the margins on the compensated loop" a stronger check than "the design formulas give PM = 55°"?
