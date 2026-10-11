# Procedure: Calibrate a scale and a length measurement

**Engineering question.** What are the true reading error and the uncertainty of the cheapest instruments you own, and can you make a 1 g kitchen scale report mass to better than 1 g by calibrating it?

Every later experiment in this course quotes an uncertainty for every measured number. This project is where you learn to produce that number honestly, on the instruments you will actually use.

## 0. Predict first

Before any reading, write down in `data/prediction.md`:

- The resolution of each instrument (the smallest step it can show), and from it the Type B standard uncertainty you expect: resolution / √12 for a digital display (uniform distribution over one step), half a graduation / √3 for a ruler read by eye.
- Your predicted reading for each reference object (its official mass, or its known length).
- Whether you expect the scale to have an *offset* (reads non-zero with nothing on it), a *scale error* (reads 2% high across the range), or *nonlinearity* (error depends on load). Guess; you will find out.

## 1. Set up

Minimum tier: kitchen scale, 10 or more identical coins of one denomination (look up and record the official mass and its source), a ruler, a notebook. Record the room temperature if you can; note it if you cannot.

Put the scale on a flat, hard table, away from drafts. Switch it on and let it settle for one minute. Zero it.

## 2. Calibrate the scale

**Repeatability (Type A).** Place one coin on the scale, read, remove, wait for zero, repeat 10 times. Record all 10 readings, not the average. Do the same with a stack of 10 coins.

**Calibration curve.** Build reference masses from 1, 2, 5, 10, 20, (30, 40, 50) coins. For each, three readings, placing the stack near the center of the platform each time. Record stack size, nominal reference mass, and each reading. Then repeat the sequence going *down* from the largest stack, to detect hysteresis.

**Placement sensitivity.** With the 10-coin stack, read it at the center and at each of four positions near the platform edges. Record all five.

If you have a second scale (recommended tier), repeat the calibration curve on it, with the same stacks, so the two instruments can be compared.

## 3. Calibrate a length measurement

Measure the diameter of one coin 10 times with the ruler, re-positioning each time, reading to the nearest half-millimeter. Then measure a stack of 10 coins' total height 5 times. Look up the official diameter and thickness. If you have calipers, repeat both with the calipers.

## 4. Analyze

Reduce the data with the lesson's code (`curriculum/phase-00/lesson-14-measurement-calibration-uncertainty`):

- Type A standard uncertainty of the mean from the repeated readings; compare with the Type B resolution estimate from step 0. Which dominates?
- Fit `reading = a + b × reference` to the calibration data. Report `a` (offset), `b` (scale factor), their uncertainties, and the residuals. Plot residuals against reference mass. Is there a trend (nonlinearity) or a loop (hysteresis)?
- Combined standard uncertainty of a *corrected* mass measurement, and the expanded uncertainty with k = 2.
- For the length measurements: mean, Type A uncertainty, comparison with the official value, and with the calipers if you have them.

## 5. Compare and explain

Where does the scale disagree with the reference beyond the combined uncertainty? Rank the candidate causes (coin wear, placement, zero drift, temperature, nonlinearity, your reference value's own uncertainty) and test at least one with a follow-up measurement; for example, re-measure the 10-coin stack after 10 minutes to check drift, or weigh a single worn coin against a shiny one. "Experimental error" is not an explanation.

## 6. Report

Fill `report-template.md`. Submit: prediction sheet, raw readings, calibration fit with residual plot, uncertainty budget (a table of every contribution and its size), the corrected-measurement recipe you will use from now on, and what you would change. Record the final instrument uncertainties in `data/instruments.md`; later projects reuse them.
