# Procedure: Time a falling object and test the free-fall prediction

**Engineering question.** Does the free-fall model t = √(2h/g) predict the fall time of a small dense object from 0.5–2.0 m to within the uncertainty of a phone-video measurement, and if not, which assumption failed?

This is the whole course in one afternoon: a model, a prediction with uncertainty written down first, a safe measurement with calibrated instruments, a comparison, and an explanation of the gap. Everything later is this, with harder physics.

## 0. Predict first

In `data/prediction.md`, before any drop:

1. Choose four drop heights between 0.5 m and 2.0 m. Measure them later; for now, write the targets.
2. For each height, compute the predicted fall time t = √(2h/g) with g = 9.81 m/s² (write down where your value of g comes from; local g varies by about ±0.03 m/s²), and propagate the uncertainty from your height uncertainty (tape measure, your reading, and how well you can define the release point: typically 5–10 mm) and from g.
3. Write down your measurement uncertainty *before* measuring: for video at frame rate f, a timing read from frames is quantized in steps of 1/f, so a reasonable Type B uncertainty for each of the start and end frames is (1/f)/√12, and for the difference √2 times that. For a stopwatch, reaction time of 0.1–0.2 s dominates.
4. Predict whether air drag matters: estimate the terminal speed of your object and compare it with the impact speed √(2gh). If impact speed is below about 10% of terminal speed, drag changes the fall time by well under 1%.

## 1. Set up

Fix the phone on its stand so the whole drop path is in frame, side-on, 2–3 m away, and level. Put the two reference tape marks (recommended tier) on the wall in frame and measure their separation. Set the camera to its highest frame rate and record the value in the notebook. Put the landing pad down. Clear the drop zone. Read the safety box once more.

## 2. Measure heights

Mark each release height on the wall with tape. Measure from the floor (top of the landing pad if the object stops there; be consistent and say which) to the bottom of the object at release. Three readings each; record all.

## 3. Drop

For each height, 8 drops (10 if time allows). Release from rest, without pushing: open fingers, or open the clip. Record the video file name against each drop in the notebook as you go.

Stopwatch variant: 10 timings per height. Record every timing, including the ones you think you fumbled, and mark them as such; do not delete data.

## 4. Extract times

Step through each video frame by frame (most phone players allow this, or use a free video player on a computer). Record the frame number at which the object first moves and the frame at which it first touches the pad. Fall time = (frames) / f. Record both frame numbers, not just the difference. Note any frame where the judgment was ambiguous.

If the reference marks are in frame, use them to check that the scale is what you measured, and to check for motion blur that could hide the exact release frame.

## 5. Analyze

Reduce the data with the lesson's code (`curriculum/phase-00/lesson-15-first-experiment`):

- For each height: mean fall time, standard deviation, Type A standard uncertainty of the mean, combined with the Type B frame-quantization term.
- Predicted time with its propagated uncertainty, from step 0, recomputed with the *measured* heights.
- Normalized error for each height: |t_pred − t_meas| / √(u_pred² + u_meas²). Values below about 2 mean agreement within uncertainty.
- Plot measured against predicted time with error bars, and the line of equality.
- Fit t² against h: the slope should be 2/g. Report your estimate of g with its uncertainty.

## 6. Compare and explain

Where normalized error exceeds 2, list candidate causes and rank them: release not from rest (a push or a drop from fingers gives an initial velocity), release frame later than the true release (motion too small to see for the first frame or two), pad compression changing the effective height, drag, frame rate different from the setting, height measured to the wrong reference. Test at least one with a follow-up: for example, re-run one height with the clip release instead of fingers, or at a higher frame rate. Say what the test showed.

If normalized error is comfortably below 1 everywhere, say what that implies about your uncertainty estimate: were you too pessimistic, and about which term?

## 7. Report

Fill `report-template.md`. Submit: prediction sheet, videos' frame tables as raw data, reduced data, the two plots, your estimate of g with uncertainty, ranked causes with the follow-up result, limitations, and what you would redesign. This report is the first entry in your portfolio.
