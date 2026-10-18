# Procedure: Servo feedback control

**Engineering question.** From a few open-loop step tests, can you identify a model of a small servo or gear motor good enough to predict, before closing the loop, the rise time, overshoot and settling time that a P or PI controller will give, and when the measurement disagrees, can you name which effect (sample delay, dead zone, saturation, quantisation, friction, slew limit) is responsible and show it with a model that includes it?

The analysis code is `curriculum/phase-09/lesson-21-home-lab-servo-feedback` (your completed `starter/main.py`). Read `safety.md` first, every session. Build and rewire with the motor pack switched off and the USB cable unplugged. A learner without the kit skips steps 1–5 and 8–10 and analyses `data/synthetic/` (SYNTHETIC data) instead, labelling every result as synthetic; the dated prediction of step 7 is still written, from the synthetic step tests, before the synthetic closed-loop file is opened.

Choose one apparatus and write down which:

- **Option A (minimum_cost):** a micro servo with its feedback wiper brought out to an ADC pin. The servo's own loop is part of the plant: the board commands an angle, the wiper reports the angle reached. Model: K e^{−Ls}/(τs + 1). Controller: PI (removes the servo's own gain error and dead band). Loop rate 50 Hz, one update per servo frame.
- **Option B (recommended):** a small gear motor on a low-voltage H-bridge, the p04 potentiometer (or an encoder) on the output shaft. The board commands a duty; the angle is the integral of the speed. Model: K/(s(τs + 1)) with an input dead zone. Controller: P with dead-zone compensation. Loop rate 100 Hz.

## 1. Build (power off)

Option A:

```
 AA pack (<= 6.0 V measured) + ──[PTC 0.5 A]──┬──────────── servo red (V+)
                                              = 220 uF (+ to V+)
 AA pack − ───────────────────────────────────┴──┬───────── servo brown/black (GND)
 board GND ──────────────────────────────────────┘          (common ground)
 board PWM pin ──────────────────────────────────────────── servo signal (orange/yellow)
 board ADC pin ──[divider if needed, see below]──────────── servo feedback wire (wiper, 4th lead)
 board input pin ──[button]── GND                            (stop button, internal pull-up)
 board USB ── computer (logging)                             pointer arm (card, <= 5 g) on the horn
```

Option B:

```
 AA pack (<= 6.0 V) + ──[PTC 0.5 A]──┬── driver VM          driver OUT1/OUT2 ── motor (two leads)
                                     = 220 uF               driver IN1, IN2 ── two board PWM pins
 AA pack − ──────────────────────────┴── driver GND ──┬──── board GND (common ground)
 board 3.3 V ── driver logic V (if the module needs it)│     driver EN/SLEEP ── board output pin
 board 3.3 V ── pot end A,  pot end B ── GND ──────────┘     pot wiper ── board ADC pin
 board input pin ──[button]── GND                            pot shaft ── silicone tube ── motor shaft
```

Mount the servo or motor (and the potentiometer) on the block so nothing can move except the output and the pointer. Tape the printed protractor behind the pointer. Check every node with continuity and look for a path from the pack's + to ground. Measure the pack: it must read ≤ 6.0 V. Write the pack voltage in `data/calibration.csv`.

## 2. Program the board and run the four tests

Copy the firmware listing at the end of this file to the board (MicroPython path; the Arduino-style path is described after the listing). Set `OPTION`, the pin numbers for your board, and leave `KP = KI = DZ_COMP = 0`. First the two one-time checks below (the wiper check is Option A only); then, with the pointer clear and fingers away, and logging every result in `data/safety_tests.md` with the date:

**Check (a), Option A: wiper voltage before it reaches the ADC.** The servo's internal potentiometer may be fed from the motor pack (up to 6 V), so its wiper can exceed the board's supply, which p04's rule "no analog pin above the board's supply voltage" forbids. Leave the feedback wire **unconnected** from the board at first. Pack on, `arm()`, and command the servo to −60° and then +60° (`test_stop(u=-60, seconds=1)`, `test_stop(u=60, seconds=1)`); at each, measure the feedback wire against GND with the multimeter, then pack off. If either reading is above the board's supply (3.3 V on most boards), connect the wire through a divider: two equal p04 resistors (10 kΩ from the wiper to the ADC pin, 10 kΩ from the ADC pin to GND), which halves 6 V to 3.0 V; then measure the ADC-pin voltage at both ends again and confirm it stays below the board supply before connecting. Record both readings and whether the divider is fitted in `data/safety_tests.md`. The two-point calibration of step 3 absorbs the divider's factor. (Option B needs no check: its potentiometer is fed from the board's own 3.3 V.)

**Check (b): stall current (once, before the four tests; required if the listing does not state a stall current ≤ 500 mA at your supply).** The stall is **commanded against a fixed stop, never made by holding the shaft or horn by hand**, and lasts **≤ 1 s**. Pack off: remove the pointer arm; tape or clamp a rigid stop (a wooden block) to the base so the horn (A) or the shaft coupling (B) butts against it; put the multimeter, on a DC current range rated for at least 1 A (usually the 10 A input), in series in the pack's + lead after the fuse.
- Option A: place the stop about 20° from centre, then pack on, `arm()`, `test_stop(u=40, seconds=1)` (the servo drives into the stop for at most 1 s), read the meter, pack off.
- Option B: disconnect the motor leads from the driver and connect them straight to the pack (through the fuse and the meter); switch the pack on for **at most 1 s**, read the meter, switch it off. Reconnect the driver afterwards.

One reading only; do not repeat it to "get a better number". Record it in `data/safety_tests.md`. **If it is above 500 mA, or the fuse trips, reject the motor or servo**: it is not permitted in this project (see `safety.md`). Remove the stop and refit the pointer arm before going on.

The four tests:

1. **Stop test.** `arm()`, then `test_stop()`: the output runs at a small open-loop duty (B: 0.3) or command (A: 10°); press the button. The printed applied value must be 0 from the next sample and the pointer must stop. Re-arm only with `arm()`.
2. **Limit test.** Motor pack **off**: `test_limit()` asks for a duty of 1.0 (A: 90°) and prints what `output()` actually applied: it must be 0.6 (A: 60°). Then pack on, the same with `test_limit(apply=True)` for 0.2 s: the pointer moves, the log shows 0.6.
3. **Sign test.** `run(0, 0.25, closed=False, n_post=20)` (A: `run(0, 10, ...)`): the logged angle must increase. If it decreases, pack off, swap the motor leads (A: change the sign of `DEG_PER_U16`), repeat.
4. **Travel-limit test.** Motor pack off (and for B, the coupling loosened if the gearbox will not turn by hand): turn the pot past `TRAVEL_LIMIT_DEG` while `watch()` runs; it must print STOP.

Do not continue until all four pass. Repeat test 1 at the start of every session.

## 3. Calibrate the angle sensor

1. Pack off. Turn the output by hand (A: command the servo with the pack on to −45° and +45° marks if the horn cannot be turned) so that the pointer is on the −45° mark, then the +45° mark; at each, record the mean of 32 `adc.read_u16()` values with `calib()`.
2. DEG_PER_U16 = 90/(u16(+45) − u16(−45)); CENTER_U16 = the mean of the two. Put both in the firmware.
3. At rest, record 200 readings: their standard deviation is the sensor noise (degrees). Record the hardware resolution: one true ADC count is 2^(16−n) u16 units for an n-bit converter, so for a 12-bit ADC one count is 16 × DEG_PER_U16 degrees.
4. Write all of it in `data/calibration.csv`.

## 4. Predict the open-loop step (dated)

In `data/prediction.md`, before the first open-loop step: for Option B, the no-load output speed at your pack voltage from the motor's listing (rpm × 6 = deg/s), hence K ≈ that speed per unit duty, a guess of τ (tens of ms for a micro gear motor) and of the dead zone (a duty between 0.05 and 0.2); for Option A, K ≈ 1, a lag τ of a few tens of ms from the servo's speed rating, and the frame latency L ≤ 20 ms. These are rough; the point is to have written something you can be wrong about.

## 5. Open-loop step tests

- **Option B.** From rest, pointer at about −60° (so the run stays inside the travel limit): steps of duty 0.3, 0.45 and 0.6 for 0.5 s each with `run(i, duty, closed=False, n_post=50)`. Three repeats of each. Return the pointer to −60° between runs with a small duty or by hand with the pack off.
- **Option A.** From the centre: command steps of 20° and 40° (`run(i, 20, closed=False, n_post=30)`), three repeats each, returning to centre between runs.

Paste the printed rows into `data/raw/YYYY-MM-DD_<option>_open_loop.csv` with the header of `data/README.md`. Never edit a raw file afterwards.

## 6. Identify the model

With the lesson's functions (or `analyse_dataset(folder, files=MEASURED_FILES)` after copying the raw files to those names):

- B: `fit_integrator_loglinear` then `refine_least_squares` on each run (K_app, τ, L with standard errors); `fit_dead_zone` on the duties and the fitted slopes K_app × duty.
- A: `fit_first_order_loglinear` then `refine_least_squares` on the 20° runs (K, τ, L); fit the 40° runs too and compare (a larger τ or L means the servo's slew limit is active).
- Record the spread over the three repeats as well as the fit's standard errors.

## 7. Design and predict the closed loop (dated, before closing the loop)

- B: `design_p_position(K, τ, 10)` → Kp; round it for the firmware, down if anything (a lower gain lowers the overshoot); dead-zone compensation slightly below the fitted dead zone (for example 0.09 for a fitted 0.10: under-compensating leaves a small error, over-compensating makes the output chatter).
- A: `design_pi_first_order(K, τ, 0.4)` → Kp, Ki; round both.
- `predict_closed_loop` for each planned step (B: 15° and 60°; A: 20° and 45°): rise time (10–90 %), overshoot, 2 % settling time, steady-state error; the peak duty (B: Kp × step + compensation) or command; whether it reaches the limit.

Write the table in `data/prediction.md` with the **date and time**, the identified parameters and the gains. Commit it, or photograph the page, before the first closed-loop run. Do not revise it after measuring; add a dated correction below if you find a mistake.

## 8. Close the loop

1. Repeat the stop test. Fingers clear.
2. First run at **half** the designed gain and a 10° step (`run(i, 10, closed=True)`): it must settle. If it oscillates, stop and recheck the sign test.
3. The designed gain: each planned step three times (`n_post=100` for B at 100 Hz, `n_post=60` for A at 50 Hz). Return to the start between runs.
4. Paste the rows into `data/raw/YYYY-MM-DD_<option>_closed_loop.csv`.

## 9. Analyse and compare

- `step_metrics` on each closed-loop run (`ref` = the step); mean and standard deviation over the repeats.
- `compare_metrics` against the dated prediction.
- The peak applied duty or command from the log: was the limit reached?
- The figure `servo_feedback_lab.png` from your data.

## 10. Explain, improve, re-predict, re-measure

1. For every metric that differs from the prediction by more than its repeat spread, list candidate causes and rank them: the one-sample firmware delay and the hold (about 1.5 samples of lag), the plant dead time L, the duty or command limit, the dead zone and its compensation, quantisation and noise in the 2 % band, Coulomb friction asymmetry, the servo's slew limit or dead band, a battery voltage that has sagged since the step test (K scales with it).
2. Test the first cause with a model that includes it: `predict_closed_loop` with the delay as a lag 1/(θs + 1), θ = 1.5 Ts + L; then `simulate_closed_loop` with Ts, `delay_samples=1`, `plant_delay=L`, the limit, the dead zone and the quantisation.
3. B: redesign with `design_p_sampled` for the same 10 % target; A: choose new gains on the sampled model for your target. Write a **new dated prediction**, run it three times, compare.

"Experimental error" is not an explanation.

## 11. Tidy up

Stop button, motor pack off, USB unplugged. Discharge the capacitor through 1 kΩ if you remove it. Remove the cells from the pack. Touch-test the motor and the fuse: they should be at room temperature.

## 12. Report

Fill `report-template.md`. Submit: the logged safety tests, the calibration file, both dated predictions, the raw data, the analysis script, the plots, the comparison tables, the ranked causes with the improved-model result, the redesign and its measurement, the completed safety box.

---

## Firmware listing (MicroPython)

Pin numbers are examples for a Raspberry Pi Pico-style board; use your own board's PWM-capable, ADC-capable and input pins and check them in your port's quick reference. The control law in `pi_update` is the same as `pi_update` in the lesson code, so the simulation and the board run the same algorithm. Each loop **first writes the value computed in the previous loop, then reads the sensor and computes the next value**: the delay is exactly one sample and does not depend on how long the computation takes.

```python
# p09-servo-control firmware (MicroPython). OPTION "B": gear motor + H-bridge; "A": feedback servo.
from machine import Pin, PWM, ADC
import time

OPTION = "B"
TS_US = 10_000 if OPTION == "B" else 20_000   # 100 Hz motor loop; one update per 50 Hz servo frame
U_LIMIT = 0.6 if OPTION == "B" else 60.0      # software duty limit (B) / command limit in deg (A). Never raise.
TRAVEL_LIMIT_DEG = 110.0                      # stop beyond this angle from centre (inside the pot's end stops)
KP, KI, DZ_COMP = 0.0, 0.0, 0.0               # set from your design only after the four tests pass
DEG_PER_U16 = 300.0 / 65535                   # from your calibration (procedure step 3)
CENTER_U16 = 32768                            # from your calibration

stopped = True                                # latched; only arm() clears it


def _on_stop(pin):
    global stopped
    stopped = True


stop_pin = Pin(15, Pin.IN, Pin.PULL_UP)       # button from the pin to GND
stop_pin.irq(trigger=Pin.IRQ_FALLING, handler=_on_stop)
adc = ADC(Pin(26))
if OPTION == "B":
    in1, in2 = PWM(Pin(16)), PWM(Pin(17))
    for p in (in1, in2):
        p.freq(20_000)
        p.duty_u16(0)
    enable = Pin(18, Pin.OUT, value=0)        # driver disabled until armed
else:
    servo = PWM(Pin(16))
    servo.freq(50)
    servo.duty_u16(0)


def arm():
    global stopped
    stopped = False


def output(u):
    """Apply u through the software limit (zero if stopped); return the value actually applied."""
    u = max(-U_LIMIT, min(U_LIMIT, u))
    if stopped or stop_pin.value() == 0:
        u = 0.0
    if OPTION == "B":
        d = int(abs(u) * 65535)
        in1.duty_u16(d if u > 0 else 0)
        in2.duty_u16(d if u < 0 else 0)
        enable.value(0 if stopped else 1)
    elif stopped:
        servo.duty_u16(0)                     # no pulses (an analog servo goes limp; a digital one may hold)
    else:
        pulse_us = 1500 + u * 500 / 90        # check your servo's pulse range and direction
        servo.duty_u16(int(pulse_us * 65535 / 20_000))
    return u


def angle_deg():
    return (adc.read_u16() - CENTER_U16) * DEG_PER_U16


def pi_update(integ, e, kp, ki, ts, umin, umax, dzc):
    """P or PI with clamping anti-windup and dead-zone compensation (same law as the lesson code)."""
    i_new = integ + e * ts
    u = kp * e + ki * integ + (dzc if e > 0 else (-dzc if e < 0 else 0.0))
    if u > umax:
        if e > 0:
            i_new = integ
        u = umax
    elif u < umin:
        if e < 0:
            i_new = integ
        u = umin
    return u, i_new


def run(run_id, step, closed, n_pre=10, n_post=100):
    """One logged step from rest. Open loop: step = duty (B) or command deg (A) from t = 0.
    Closed loop: step = reference angle change in deg. Angles are relative to the starting angle.
    Rows are buffered and printed after the run so printing cannot disturb the timing."""
    global stopped
    ts = TS_US / 1e6
    th0 = sum(angle_deg() for _ in range(8)) / 8
    u_pending, integ, rows, overruns = 0.0, 0.0, [], 0
    t_next = time.ticks_us()
    for k in range(-n_pre, n_post + 1):
        u = output(u_pending if closed else (step if k >= 0 else 0.0))   # write first
        y = angle_deg()                                                  # then read
        if abs(y) > TRAVEL_LIMIT_DEG:
            stopped = True
            output(0.0)
        y -= th0
        r = step if (closed and k >= 0) else 0.0
        if closed:
            u_pending, integ = pi_update(integ, r - y, KP, KI, ts, -U_LIMIT, U_LIMIT, DZ_COMP)
            rows.append(f"{run_id},{k * ts:.3f},{r:.4f},{y:.4f},{u:.4f},{KP:.4f},{KI:.4f},measured")
        else:
            rows.append(f"{run_id},{k * ts:.3f},{u:.4f},{y:.4f},measured")
        if stopped:
            break
        t_next = time.ticks_add(t_next, TS_US)
        if time.ticks_diff(t_next, time.ticks_us()) < 0:
            overruns += 1
        while time.ticks_diff(t_next, time.ticks_us()) > 0:
            pass
    output(0.0)
    for row in rows:
        print(row)
    print("# overruns", overruns, "stopped" if stopped else "")


def test_stop(u=None, seconds=5):
    """Safety test 1: run a small open-loop output; press the button; the output must go to 0."""
    u = (0.3 if OPTION == "B" else 10.0) if u is None else u
    t0 = time.ticks_ms()
    while time.ticks_diff(time.ticks_ms(), t0) < seconds * 1000:
        applied = output(u)
        print(time.ticks_diff(time.ticks_ms(), t0), applied)
        if applied == 0.0:
            print("STOP OK")
            break
        time.sleep_ms(50)
    output(0.0)


def test_limit(apply=False):
    """Safety test 2: ask for more than the limit; print what is actually applied."""
    ask = 1.0 if OPTION == "B" else 90.0
    if not apply:
        print("asked", ask, "would apply", max(-U_LIMIT, min(U_LIMIT, ask)))
        return
    t0 = time.ticks_ms()
    while time.ticks_diff(time.ticks_ms(), t0) < 200:
        print("asked", ask, "applied", output(ask))
        time.sleep_ms(20)
    output(0.0)


def watch(seconds=20):
    """Safety test 4: print the angle; STOP when the travel limit is passed."""
    global stopped
    t0 = time.ticks_ms()
    while time.ticks_diff(time.ticks_ms(), t0) < seconds * 1000:
        y = angle_deg()
        if abs(y) > TRAVEL_LIMIT_DEG:
            stopped = True
            output(0.0)
            print("STOP at", y)
            return
        print(round(y, 2))
        time.sleep_ms(100)


def calib(n=32):
    """Mean of n raw readings, for the two-point calibration."""
    return sum(adc.read_u16() for _ in range(n)) / n
```

**Arduino-style path.** The same structure works in any Arduino-style sketch: a `loop()` paced by `micros()` that first writes the pending output (`analogWrite` to the two driver inputs, or a servo library's write-microseconds call), then reads the angle (`analogRead`, counts at the hardware resolution), then computes the next output with the same `pi_update`; a stop flag set by an interrupt on the button pin and latched until a deliberate re-arm command over the serial port; the limit enforced inside the one function that writes the output; rows buffered in an array and printed after the run. Write the four safety tests as serial commands and log their results the same way.
