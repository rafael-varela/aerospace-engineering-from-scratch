# P00.15  First Experiment: Predict, Measure, Compare 🔬

> Write the prediction down before you measure. Everything that follows is only honest if you did.

| | |
|---|---|
| **Learning outcomes** | 1. Derive t = √(2h/g) from constant acceleration and propagate the uncertainties of h and g into t. 2. Design a timing measurement, identify its dominant uncertainty (frame quantisation or reaction time), and reduce repeated timings to a mean with a Type A and a Type B term. 3. Compare a prediction with a measurement using the normalized error E_n and interpret values below and above 2. 4. Derive the first-order effect of air drag on fall time and decide, before dropping anything, whether your setup can detect it. 5. Rank candidate causes of a discrepancy and test one with a follow-up measurement. |
| **Prerequisites** | `P00.04`, `P00.13`, `P00.14`; skill check: can propagate uncertainty through a square root; has written a safety box and run the checker on it. |
| **Engineering question** | For a ball dropped from 0.6 m to 2.0 m, does the measured fall time agree with t = √(2h/g) within the combined uncertainty of prediction and measurement, and if not, which assumption failed? |
| **Assumptions and limits** | Constant g over 2 m (true to one part in 10⁶); release from rest; no air drag; a rigid floor and a point of impact you can identify to one frame. Drag is the first assumption to fail: for a tennis ball it lengthens a 1.8 m fall by about 0.6 %, below a 30 fps camera's resolution; for a ping-pong ball it is 4 % and detectable. A non-zero release velocity from a hand that pushes shifts every time by the same ~v₀/g. |
| **Analytical work** | Derivation of the fall-time law and its uncertainty; the fall with a non-zero release velocity; the quadratic-drag solution x(t) = (v_t²/g) ln cosh(gt/v_t), its inversion, and the small-drag expansion t ≈ t₀(1 + gh/(6v_t²)); the frame-quantisation uncertainty; the normalized error. |
| **Computational deliverable** | `starter/main.py`: `predict_fall_time`, `frame_timing`, `reduce_timings`, `normalized_error`, `compare`, `fit_g`, `fall_time_with_drag`, `plot_comparison`, plus `make_synthetic_timings` (simulation, for tests only). Tests check the formula and its √h scaling, the propagated uncertainty against a finite difference, the drag-free limit and monotonic drag effect, the distance law of the drag solution, 1/fps scaling of frame uncertainty, E_n properties, that a synthetic dataset consistent with free fall gives E_n < 2 at every height while a biased one gives E_n > 2, recovery of g from a fit, and that the plot file is written. |
| **Practical activity** | **physical experiment**: drop a ball from five heights, time each drop at least five times by counting video frames (or with a stopwatch as a comparison method), and compare with the prediction written beforehand. See `projects/p00-first-experiment/`. |
| **Equipment** | Level 1: a ball (tennis ball recommended; a ping-pong ball as a second object), a tape measure, a phone camera with a known frame rate, something to hold the phone still, a soft mat. Optional: a stopwatch, a 120 fps or 240 fps slow-motion mode. See `projects/p00-first-experiment/bom.yaml`. |
| **Safety** | See `projects/p00-first-experiment/safety.md`; the box you wrote in `P00.13` is the draft of it. |
| **Evidence** | `data/prediction.md` (dated before any drop), `data/raw/*.csv`, the comparison plot with k = 2 error bars, the E_n table, the ranked list of discrepancy causes with one follow-up test, the report. |
| **Assessment** | `assessments/phase-00/grading-rubric.md`; all four competencies. Threshold: a prediction with uncertainty recorded before measurement; E_n computed at every height; every E_n > 2 explained by a named, tested cause (never "experimental error"). |
| **References** | `@jcgm-100-2008` §5.1 "Uncorrelated input quantities" and §7 "Reporting uncertainty" (supplementary-free; the propagation and reporting rules used here were taught in `P00.14`). |
| **Extensions** | ★ Detect drag: drop a ping-pong ball and a tennis ball side by side at 120 fps and measure the arrival gap. ★ Solve the drag equation of motion with `scipy.integrate.solve_ivp` and compare with the closed form. ★ Estimate your own reaction time from stopwatch-versus-video differences and its scatter. |
| **Time** | unmeasured |

## 1. The engineering problem

Every validation in engineering has the same shape. Someone predicts a number with a model, someone measures it, and the two are compared with their uncertainties. A wind-tunnel drag polar is compared with a CFD prediction. A flight-test climb rate is compared with the performance model. A rocket motor's measured thrust curve is compared with the design curve. The comparison is only meaningful if the prediction was made *first*, with an uncertainty, and the measurement carries its own uncertainty. Otherwise a prediction can be adjusted until it matches, and nobody learns anything.

This lesson is that whole shape at the smallest scale that is still real. The model is the oldest in mechanics: an object released from rest falls a height h in time t = √(2h/g). You will predict t for five heights, with uncertainty, write the predictions down, then drop a ball and time it with your phone's camera, and then compare. The physics is simple on purpose. What you are learning is the method: how to budget a prediction, how to turn frame numbers into a time with an uncertainty, how to say in one number whether two values agree, and what to do when they do not.

Success criteria. Before any drop: a table of t_pred ± u_pred for each height, dated, in `data/prediction.md`. After: a table of t_meas ± u_meas at each height from at least five drops, and the normalized error E_n = |t_pred − t_meas| / √(u_pred² + u_meas²). Agreement means E_n < 2 at every height. Disagreement at any height is not failure; it is the interesting part, and it must end with a named cause and a follow-up measurement that tested it.

## 2. Understand and derive

### 2.1 Free fall from rest

Take a vertical axis x pointing down with origin at the release point. The only force on the ball is its weight mg (assumption 1: no drag; assumption 2: g constant over the drop, which holds to better than one part in a million over 2 m). Newton's second law gives dv/dt = g. Integrating from rest (assumption 3: v(0) = 0),

v(t) = g t,   x(t) = ½ g t².

Setting x = h, the height of the *bottom* of the ball above the floor at release (this matters: see §5), and solving:

t = √(2h / g).

Check the units: √(m / (m/s²)) = √(s²) = s. Check the limits: h → 0 gives t → 0; doubling h multiplies t by √2, not 2. For h = 1.800 m and g = 9.807 m/s², t = 0.6059 s.

### 2.2 The prediction's uncertainty

From `P00.14`, t ∝ h^½ g^−½, so the power rule gives

u_t / t = ½ √( (u_h/h)² + (u_g/g)² ).

Two inputs. You measure h with a tape: a careful reading to the bottom of a ball held against a mark is good to about 5 mm (u_h = 0.005 m, treated as a standard uncertainty from a ±9 mm rectangular band). You do not know your local g exactly: it ranges from 9.78 at the equator to 9.83 m/s² at the poles and falls with altitude. Unless you look it up for your city, use 9.807 with u_g = 0.03 m/s² (generous). At h = 1.8 m: u_h/h = 0.0028, u_g/g = 0.0031, so u_t/t = 0.0021 and u_t = 1.3 ms. At h = 0.6 m the height term triples to 0.0083 and u_t = 1.5 ms. The prediction is good to a couple of milliseconds. The measurement will not be, and that is the first thing a budget tells you: the comparison will be limited by the camera, not by the model.

### 2.3 A hand that pushes

If the ball leaves your fingers with a small downward speed v₀ instead of from rest, then h = v₀t + ½gt², and

t = ( −v₀ + √(v₀² + 2gh) ) / g.

Differentiating at v₀ = 0 gives dt/dv₀ = −1/g: every 0.1 m/s of release speed shortens the fall by 10 ms, *independent of height* to first order. At 30 fps a frame is 33 ms, so a slight push is a third of a frame, every time, in the same direction. This is a **systematic** effect: averaging more drops does not remove it. Its signature in the data is an offset that is the same at every height. Control it: open two fingers that were holding the ball from below, rather than letting go from above, or rest the ball on a ruler and tilt the ruler away.

### 2.4 Air drag as a limiting case

Drag on a ball at these speeds is roughly proportional to v²: F_D = ½ρ C_D A v² with ρ the air density (1.2 kg/m³), C_D about 0.5 for a sphere, and A = πd²/4 the frontal area. The equation of motion is m dv/dt = mg − ½ρC_D A v². The ball stops accelerating when the two terms balance, at the **terminal speed**

v_t = √( 2mg / (ρ C_D A) ).

Tennis ball (m = 58 g, d = 67 mm): v_t ≈ 23 m/s. Ping-pong ball (2.7 g, 40 mm): v_t ≈ 8.4 m/s. Writing the equation as dv/dt = g(1 − v²/v_t²) and integrating from rest gives v(t) = v_t tanh(gt/v_t), and integrating again,

x(t) = (v_t² / g) ln cosh(g t / v_t).

Setting x = h and inverting, t = (v_t/g) arccosh( exp(g h / v_t²) ). Two checks on this closed form. As v_t → ∞ (no drag), ln cosh(u) ≈ u²/2 for small u, so h ≈ ½gt² and t → √(2h/g). And expanding one order further, ln cosh u ≈ u²/2 − u⁴/12, and solving for t to first order in the correction gives

t ≈ t₀ ( 1 + g h / (6 v_t²) ),   t₀ = √(2h/g).

The drag correction grows with h and falls with v_t². Numbers at h = 1.8 m: tennis ball, gh/(6v_t²) = 17.7/(6 × 529) = 0.0056, so +3.4 ms on 606 ms; ping-pong ball, 17.7/(6 × 70) = 0.042, so +25 ms. At 30 fps (frame 33 ms, timing uncertainty about 15 ms after averaging) the tennis ball's drag is invisible and the ping-pong ball's is marginal; at 120 fps both the ping-pong ball's drag and the difference between the two balls are clearly measurable. That is the right way to think about an assumption: not "is drag negligible?" but "can *my instrument* see it?" You decide before the experiment, and you can design the experiment to make the answer yes (★ extension).

A numerical detail you will meet in the code: for weak drag, g h / v_t² is tiny and `acosh(exp(a))` loses almost all its digits to cancellation. Writing the same quantity as log1p(x + √(x(x+2))) with x = expm1(a) keeps full precision; a test checks the drag-free limit to one part in 10⁹.

### 2.5 Timing with a camera

A phone camera records frames at a fixed rate f (30, 60, 120 or 240 per second). Step through the video and find the first frame in which the ball has visibly left the fingers (release) and the first in which it touches the floor (impact). The fall time is t = N / f with N the frame difference. Each event happened somewhere inside its frame interval, anywhere with equal probability, so each end carries a rectangular uncertainty of half-width ½ frame: u = (0.5/f)/√3 per end, and the two ends combine in quadrature to

u_frame = √2 × 0.5 / (√3 f) = 0.41 / f.

At 30 fps that is 13.6 ms; at 120 fps, 3.4 ms. This is a Type B uncertainty that applies to every drop in the same way: it does **not** shrink when you average drops, because it comes from the instrument's quantisation, not from scatter between drops. (Strictly, different drops land at different phases within a frame, so part of it does average down; treating it as fully systematic is the conservative choice and keeps the budget honest.) The drop-to-drop scatter (your release, the ball's spin, the frame you chose) is the Type A part, s/√n, and it does shrink.

**Worked example: one height, five drops.** Camera at 30.0 fps, h = 1.200 m. Frame differences N = 15, 15, 14, 15, 15. Times t = N/30 = 0.5000, 0.5000, 0.4667, 0.5000, 0.5000 s. Mean 0.4933 s; sample standard deviation 0.0149 s; Type A u_A = 0.0149/√5 = 0.0067 s. Frame term u_frame = 0.41/30 = 0.0136 s. Combined u_meas = √(0.0067² + 0.0136²) = 0.0152 s. Prediction from §2.2: t_pred = 0.4946 s, u_pred = 0.0013 s. Then

E_n = |0.4946 − 0.4933| / √(0.0013² + 0.0152²) = 0.0013 / 0.0153 = 0.09.

Consistent, and the budget shows why the test is not sharper: the frame term is twice the drop-to-drop scatter and ten times the prediction's uncertainty. Averaging another twenty drops would bring u_A to 0.003 s and u_meas only to 0.014 s. Switching to 120 fps brings u_frame to 0.0034 s and u_meas to about 0.0075 s at once.

Two things to verify rather than trust. The frame rate: film a running stopwatch for ten seconds and count frames; "30 fps" is sometimes 29.97, and slow-motion modes sometimes deliver fewer unique frames than advertised. And the shutter: a fast-moving ball blurs into a streak; pick a consistent rule for which frame counts as "impact" (first frame with the ball in contact, or the frame of maximum squash) and apply it the same way every time.

### 2.6 Timing with a stopwatch, and why it is the comparison method

A hand-operated stopwatch adds your reaction time at start and at stop. Human visual reaction time is about 0.2 s with a scatter of 0.03 to 0.05 s. If you both release the ball and start the watch, the start is anticipated and nearly unbiased; the stop, reacting to the impact, is late by about 0.2 s. On a 0.6 s fall that is a 30 % systematic error. Repeating the measurement reduces the *scatter* (the Type A term) but not the *bias*: the mean of ten late stops is still late. The signature in the comparison: E_n large at every height, with the measured time always longer than predicted by about the same amount. A stopwatch is still worth using here, as the method that shows what a systematic error looks like next to a method that does not have it.

### 2.7 The comparison in one number

Prediction t_p ± u_p, measurement t_m ± u_m, independent. Their difference has uncertainty √(u_p² + u_m²), and the **normalized error** is

E_n = |t_p − t_m| / √(u_p² + u_m²).

If the budgets are honest and the model is right, |t_p − t_m| is below twice its uncertainty about 95 % of the time: E_n < 2 means "consistent at k = 2". E_n > 2 at one height could be chance (one in twenty); E_n > 2 at every height with the same sign is a systematic effect; E_n > 2 only at the tallest heights, growing with h, looks like drag; growing toward the *shortest* heights looks like a height reference error (§5). The number does not tell you what is wrong. It tells you, reliably, *that* something is, and its pattern across heights tells you where to look.

One honest caveat. E_n < 2 is weak praise when u_m is huge: a stopwatch with u_m = 0.1 s "agrees" with almost anything. Report the uncertainties next to E_n so a reader can see how sharp the test was.

## 3. Build the computational model

### Step 1: Predict with uncertainty

```python
def predict_fall_time(h, u_h, g=G_STANDARD, u_g=U_G_LOCAL):
    t = math.sqrt(2.0 * h / g)
    u_t = 0.5 * t * math.sqrt((u_h / h) ** 2 + (u_g / g) ** 2)
    return t, u_t
```

A test compares `u_t` with a central finite difference of √(2h/g) (the way `P00.14`'s `propagate` would do it); they must agree to one part in 10⁶.

### Step 2: Frame counts to times

```python
def frame_timing(n_frames, fps):
    return n_frames / fps, math.sqrt(2.0) * 0.5 / (math.sqrt(3.0) * fps)
```

### Step 3: Reduce repeated drops

```python
def reduce_timings(times, u_single=0.0):
    x = np.asarray(times, dtype=float)
    u_a = x.std(ddof=1) / math.sqrt(x.size)
    return float(x.mean()), math.hypot(u_a, u_single)
```

`u_single` is the frame term from step 2; `math.hypot` adds it in quadrature, and because it is added *after* the division by √n it does not average down. A test with fifty identical timings checks exactly that.

### Step 4: Compare

`normalized_error` is one line plus a guard for the zero-uncertainty case. `compare` loops over heights, calls steps 1 and 3, and returns a `Comparison` of arrays. `fit_g` fits h = (g/2) t² through the origin: slope b = Σ(t²h)/Σ(t²)², g = 2b, with an uncertainty from the residual scatter. It is a second look at the same data that does not assume a value of g at all.

### Step 5: Drag

Implement `fall_time_with_drag(h, v_terminal, g)` from §2.4 with the `expm1`/`log1p` form, returning √(2h/g) when `v_terminal` is infinite and using arccosh(eᵃ) ≈ a + ln 2 when a > 300 so `exp` cannot overflow.

### Step 6: Plot

`plot_comparison` draws the curve √(2h/g), the predictions and the measurements with k = 2 error bars against height, and a lower panel of E_n with a dashed line at 2. It saves to a path and closes the figure; `main()` never calls `plt.show()`.

### Step 7: A synthetic dataset for the tests

`make_synthetic_timings` produces frame-quantised timings that are, by construction, consistent with free fall: the true time plus an optional reaction-time error, floored to whole frames at a random phase. It is labelled as simulation in its docstring and exists so the tests can check the pipeline on data whose truth is known. It is not data, and nothing from it goes in your report.

### Verify

`tests/test_main.py` checks the formula and the √h scaling; zero input uncertainty → zero output uncertainty; u_t linear in u_h; the finite-difference agreement; the drag-free limit (both `inf` and v_t = 10⁶ m/s, to 10⁻⁹); drag lengthening the fall more for lower v_t; the tennis-ball correction under 5 % at 2 m; the drag solution satisfying x(t) = h to 10⁻⁹; frame uncertainty halving at double frame rate; Type A reduction; the Type B term not averaging down; E_n symmetry, zero and unity cases; four seeded synthetic datasets at 30 fps all giving E_n < 2 at every height with monotone times; a dataset with a 0.15 s stopwatch bias giving E_n > 2 everywhere; `fit_g` within 5 % of g; and the plot file written to `tmp_path`.

### Compare

Three independent cross-checks, each a few lines. (1) `P00.14`'s `propagate(lambda p: sqrt(2 p[0]/p[1]), [h, g], [u_h, u_g])` against `predict_fall_time`: same to five figures. (2) `np.linalg.lstsq(t2[:, None], h)` against `fit_g` for the slope through the origin. (3) `scipy.integrate.solve_ivp` on dv/dt = g(1 − v²/v_t²), dx/dt = v, with an event at x = h, against `fall_time_with_drag`: agreement to the solver's tolerance (set `rtol=1e-10`). If any of the three disagrees, one of the two sides has a bug, and you have just found it before it reached your report.

## 4. Experiment at home, or investigate real data

**Mode: physical experiment.** The bill of materials, the authoritative safety box and the full procedure are in `projects/p00-first-experiment/` (`bom.yaml`, `safety.md`, `procedure.md`). Read the safety box first, every session; it is the one you drafted in `P00.13`, finished.

> **Safety (summary; the authoritative box is `projects/p00-first-experiment/safety.md`).** Hazards: falling mass (foot); bounce (eye, phone); trip (tape measure). Energy limits: mass ≤ 0.2 kg, height ≤ 2.0 m (3.9 J), phone battery only. Controls: soft mat at the impact point; safety glasses; phone on a stand 1.5 m from impact; tape measure fixed flat and removed after marking; 1 m exclusion circle. Supervision: none needed for an adult; children only with an adult present. Location: indoors, hard floor with a mat, clear of furniture edges, door closed. Stop criteria: anyone in the circle; a bounce toward a person; any floor damage; experimenter tired or rushed. Not permitted: ladders, chairs or windows as release points; metal or glass objects; drops onto people or pets.

The setup, seen from the camera's side:

![Side view of the drop test: a wall or door frame with height marks at 0.6, 0.9, 1.2, 1.5 and 1.8 m; a ball held from below by two fingertips with its bottom at 1.8 m; the drop height h measured from the bottom of the ball to the floor; a mat on the floor with a 1 m exclusion circle around it; and a phone on a stand 1.5 m back with its lens level with mid-fall.](figures/drop-test-setup.svg)
*Figure 1. The drop-test setup seen from the camera’s side. The height h is measured to the bottom of the ball, and the lens is level with the middle of the fall.*

Raw data go in `data/raw/drops.csv`, one row per drop, exactly as read from the video and never edited:

```
height_m,drop,release_frame,impact_frame,fps_nominal,note
1.200,1,212,227,30,
1.200,2,418,433,30,
1.200,3,655,669,30,ball clipped mark on release
```

The `note` column is where a bad drop is recorded, not deleted: the analysis may exclude it, but the raw file keeps it.

Procedure in outline:

1. **Predict first.** For heights 0.60, 0.90, 1.20, 1.50, 1.80 m, compute t_pred ± u_pred with `predict_fall_time` and write them, with the date and the inputs (u_h, g, u_g), in `data/prediction.md`. This file is frozen from now on.
2. **Calibrate the camera.** Film a running stopwatch (or a phone's stopwatch) for 10 s and count frames; record the measured frame rate in `data/calibration.csv`. Decide and write down your frame rule for release and impact.
3. **Mark heights.** Fix the tape, mark each height on a wall or door frame, and remember: h is from the **bottom of the ball** to the floor.
4. **Drop.** Camera on its stand with the whole fall in frame and the phone perpendicular to the fall. For each height, at least five drops, releasing by opening two fingers from below. Record `height_m, drop, release_frame, impact_frame` in `data/raw/drops.csv`; never edit it afterwards. Randomise the order of heights if you can.
5. **Comparison method (optional).** Time five drops at one height with a stopwatch as well, started by the hand that releases and stopped by ear or eye.
6. **Second object (optional, ★).** Repeat at 1.8 m with a ping-pong ball, and if your phone has a 120 fps mode, drop both balls side by side and film the arrival gap.

## 5. Validate and improve

Run `compare` on your data and look at the table and the plot. Then work through the three comparisons.

**Prediction against measurement.** For each height, E_n. Suppose all five are below 2: the model survives at the sharpness your camera allowed, and you say so, with u_meas quoted so a reader can judge the sharpness. Suppose instead the pattern is one of these:

| Pattern across heights | Likely cause | Test |
|---|---|---|
| t_meas shorter than t_pred by a constant ~10 ms | release with a push (§2.3) | change the release method, repeat at one height |
| t_meas longer by a constant, large amount | stopwatch reaction bias; or the impact frame rule picks a late frame | compare stopwatch with video at the same drops; re-read the frames with the other rule |
| discrepancy grows with h, t_meas longer | drag (§2.4), or the floor is not where you think at the far marks | compute the drag correction; re-measure the tall heights |
| discrepancy largest at the *shortest* height, t_meas shorter | h measured to the ball's centre, not its bottom: a constant 34 mm error is 6 % of 0.6 m but 2 % of 1.8 m | re-measure one height to the bottom of the ball |
| one height off, the rest fine | one bad frame reading, or a mis-marked height | re-read that height's frames; re-measure the mark |

The last row is why you drop at five heights and not one: a pattern across heights is a diagnosis; a single number is not.

**Hand calculation against code.** Take one height and reproduce t_pred, u_pred, the mean, u_A, u_frame and E_n with a calculator from the raw frame counts. Everything must match to the printed digits.

**A second look at the same data.** `fit_g` returns g from the slope of h against t². With honest video data it lands within a few percent of 9.81 m/s², with an uncertainty that is itself a check on your budget: if u_g from the fit is much smaller than the scatter of your E_n values suggests, the fit's uncertainty (which assumes independent scatter) is missing a systematic term.

**Improve and re-run.** Choose the dominant uncertainty and attack it. If u_frame dominates (it will at 30 fps), the single most effective change is the slow-motion mode: 120 fps cuts it by four. If the Type A term dominates, your release or frame rule is inconsistent; fix the procedure before adding drops. Re-run `compare`, put the two tables side by side in the report, and say what changed. If a discrepancy remains, say what you *would* do next and why you did not (for example, "drag is the remaining candidate at 1.8 m; a 240 fps camera would resolve it; mine stops at 60").

## 6. Defend your work

Submit: `data/prediction.md` (dated before the first drop), `data/calibration.csv`, `data/raw/drops.csv`, your analysis script, the comparison plot, the E_n table with u_pred and u_meas, the ranked cause table with the follow-up test you ran and its result, and a report following `projects/_template/report-template.md`.

### Adversarial self-checks

- What if h doubles? t grows by √2; does u_pred grow, shrink, or stay? (Work it out: the height term shrinks, the g term is constant in relative terms.)
- What if the ball's mass doubles at the same diameter? Which quantity changes, by how much, and can your camera see it?
- Which assumption fails first in your real setup: release from rest, no drag, or "the floor is at h = 0"?
- What measurement could falsify the model outright? (A fall time that scales as h, not √h, would.)
- Does the conclusion survive the frame rule? Re-read every impact frame with the alternative rule and recompute E_n.
- Which decision changes if u_meas were twice as large? Would you still claim to have seen (or ruled out) drag?

## Check yourself

1. Write the fall-time law and state the three assumptions behind it.
2. By what factor does the fall time change when the height doubles?
3. Write u_t/t in terms of u_h/h and u_g/g.
4. What is the frame-quantisation standard uncertainty of a fall time at 60 fps?
5. Define E_n and state the threshold for "consistent at k = 2".
6. A stopwatch is stopped late by about 0.2 s on every drop. Does taking more drops fix this? What does E_n look like across heights?
7. Does air drag lengthen or shorten the fall, and which ball shows it more, a tennis ball or a ping-pong ball? Give the first-order formula.
8. For a 0.6 m drop filmed at 30 fps with u_h = 5 mm, which is larger: the prediction's uncertainty or the measurement's frame term? By roughly what factor?
