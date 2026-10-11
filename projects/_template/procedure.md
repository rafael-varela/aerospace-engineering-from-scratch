# Procedure: <Project title>

**Engineering question.** …

## 0. Predict first

Before touching the apparatus, compute and write down your predicted result with its uncertainty, using the lesson's model. Record it in `data/prediction.md`. The comparison in step 5 is against *this* number; you may not revise it after measuring.

## 1. Build

Assembly steps by tier (minimum-cost / recommended / enhanced). Drawings or cutting patterns in `drawings/`. Photos are optional and only where they help.

## 2. Calibrate

How to calibrate every instrument against a known reference (a known mass, a known length, a known voltage). Record the calibration and its uncertainty in `data/calibration.csv`.

## 3. Measure

Baseline (unloaded) first. Then the test matrix: at least N levels, each repeated at least 3 times, randomized order where practical. Stay inside the energy limits in `safety.md`. Raw data to `data/raw/` in the format given in `data/README.md`; never edit raw files.

## 4. Analyze

Reduce the data with the lesson's code. Propagate uncertainty (see P00.14). Plot prediction and measurement on the same axes with error bars.

## 5. Compare and explain

Where do prediction and measurement disagree beyond their uncertainties? List candidate causes, rank them, and test at least one with a follow-up measurement. "Experimental error" is not an explanation.

## 6. Report

Fill `report-template.md`. Submit: assumptions sheet, calculations, code, raw data, calibration, plots, error analysis, interpretation, limitations, and what you would redesign.
