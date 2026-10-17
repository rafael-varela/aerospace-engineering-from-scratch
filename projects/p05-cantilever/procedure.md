# Procedure: Load-deflection and deflected shape of a clamped metal strip

**Engineering question.** Does a clamped metal strip deflect as δ = Pa³/(3EI), with a typical handbook modulus and the measured section, within the uncertainty of the readings? Can the measurement tell a compliant clamp from a low modulus?

This is the bench version of a structural ground test: load a structure in steps, measure its deflection, and compare with the analysis written down beforehand. Read `safety.md` first. If you have no kit, run every analysis step on the SYNTHETIC files in `data/synthetic/` (or on `make_synthetic_cantilever_data(seed=<your seed>)` from the lesson code) and label every result **SYNTHETIC**.

## 0. Predict first

In `data/prediction.md`, before any load, dated:

1. Measure the strip: width b at three points and thickness h at five points along it (caliper if you have one; write the instrument and its resolution). Use the mean and, as the standard uncertainty, the larger of the scatter's standard deviation and the instrument's resolution divided by √3.
2. Decide the clamp position and the load position a (measured from the **clamp edge**, not from the end of the strip under the clamp), and the free length L. For the lesson's strip: 300 mm free, a = 290 mm.
3. Compute I ± u_I, EI with E = 69 ± 2 GPa (typical for aluminium; for another metal, a typical handbook value with its source), the clamp stress 6Pa/(bh²) at each planned load, and the safe load for F_ty/2. Your top load is the smaller of 2 kg, the safe load, and the mass that makes δ/a = 0.10. Write it down.
4. Predict δ ± u at the load point for each planned load with `predict_deflection`, with the budget term by term (E, b, h, a). Say which term dominates. Write the predicted compliance (mm/N).
5. Choose five to seven marks along the strip for the profile (for example 50, 100, 150, 200, 250 mm, a, and the end), and predict the deflection at each under your top-but-one load.
6. Detectability of the clamp: run `fit_profile` on your noise-free predicted profile with u_w = √2 × (your reading uncertainty). Write the resulting u(θ₀) and say whether a 1 mrad clamp rotation would be detected.
7. Freeze the file. The comparison in step 5 is against these numbers.

## 1. Build

Minimum-cost tier: deburr the strip's ends. File a small notch across the top face at a. Mark the clamp edge position on the strip in pencil. Clamp the strip flat side up on the table with 10–30 mm overlap, its length perpendicular to the table edge, the clamp jaw's edge exactly at the pencil mark. Tape the reading ruler vertically to a heavy book (or hold it in a second clamp) beside the strip's end, zero downward or upward as you like, close to but not touching the strip. Put the cushion under the load point. Put the phone on its stand, lens at the height of the strip's top edge, 20–40 cm away, ruler graduations and strip edge both in focus.

Recommended tier: clamp the strip between the two steel plates, so the clamp bites on steel; measure b and h with the caliper. Enhanced tier: set the dial gauge on its stand with the plunger on the strip's top face at the load point (or at a profile mark), pre-loaded by about 2 mm so it can follow the strip down; the procedure is otherwise the same, with gauge readings in place of photos.

## 2. Calibrate

1. **Masses.** Weigh every load, including the bag, string and hook, on the kitchen scale; check the scale first against a known mass (p00-calibration). Record in `data/calibration.csv`.
2. **Reading.** Photograph the unloaded strip ten times, moving the phone away and back between photos, and read each. The standard deviation is your Type A reading uncertainty; if it is below the graduation divided by √12, use that instead.
3. **Lengths.** Measure a and L from the clamp edge three times each.

## 3. Measure

Envelope: total hanging mass ≤ 2 kg and ≤ the safe load; stop at the stop criteria. Check the clamp before each load.

1. **Load series at x = a.** For each load in increasing order: photograph (or read) unloaded; hang the load gently by hand; wait until it stops swinging (about 10 s); photograph loaded; lift the load off; photograph unloaded again and confirm it agrees with the first within 0.5 mm (permanent-set check). Then repeat the series in decreasing or shuffled order.
2. **Profile.** With the top-but-one load hanging, move the reading ruler (or the dial gauge) to each mark and read loaded; remove the load and read unloaded at each mark. Do it twice.
3. **Raw data.** Record every reading in `data/raw/YYYY-MM-DD_<tier>_loads.csv` and `..._profile.csv` (format in `data/README.md`), with photo file names; note any doubtful reading in the `note` column instead of deleting it.

## 4. Analyze

Reduce the data with the lesson's code (`curriculum/phase-05/lesson-18-home-lab-cantilever-deflection`):

- Differences loaded − unloaded; their uncertainty √2 × the reading uncertainty.
- `fit_line(P, d, u_d)`: slope ± u, intercept ± u, χ² and dof.
- `modulus_from_slope`: E ± u; the budget term by term.
- E_n of the compliance against the prediction and of E against the typical value (`en_number`).
- `fit_profile(x, w, P, a, u_w)`: EI, θ₀ and their uncertainties, χ².
- Plots: deflection against load with k = 2 error bars, the fit, and the prediction band; deflection against x with the predicted cubic and the fitted cubic plus θ₀x.

## 5. Compare and explain

Where E_n exceeds 2, or χ² is far from the dof, or the intercept differs from zero by more than 2u, list candidate causes and rank them: a measured from the wrong reference; the thickness taken from a label; clamp compliance; a hanger that slid; parallax in the photos; a swinging load; permanent set; the wrong alloy (last: aluminium alloys share a modulus within a few per cent). Use the pattern table in the lesson's part 5. Test at least one cause with a follow-up (re-measure a from the clamp edge; re-measure h at the clamp; clamp between steel plates and repeat; repeat the profile) and say what it showed.

For the clamp: did you detect θ₀, as predicted in step 0, or not? If you predicted no detection and got a "significant" θ₀, check the profile's χ² before believing it.

## 6. Report

Fill `report-template.md`. Submit: the frozen prediction, calibration file, raw data files, analysis script, both plots, the E and θ₀ table with E_n, ranked causes with the follow-up result, the completed safety box, limitations, and what you would redesign.
