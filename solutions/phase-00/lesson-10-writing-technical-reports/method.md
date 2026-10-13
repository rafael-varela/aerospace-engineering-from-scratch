# P00.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. For each uncertainty, keep one significant figure (two when the first digit is 1). 2. Round the value to the decimal place of the rounded uncertainty. 3. Write `value ± u unit`. 4. For the concise form, write the value to the same place and put the uncertainty's significant digits in parentheses immediately after.

## 2.
1. Write $\rho$ as a product of powers and read off the exponents. 2. State the law of propagation for independent inputs and small uncertainties (first-order Taylor expansion; cross terms vanish because covariances are zero). 3. Square each (exponent × relative uncertainty), add, take the square root. 4. Convert mm to m consistently or work in relative terms where units cancel. 5. Compute $\rho$, then $u_\rho$. 6. Compare the three squared terms to find the dominant one.

## 3.
1. Convert 20° to radians. 2. Evaluate $\theta_0^2/16$ and express as a percentage of $T_0$. 3. Compare with the observed 2 %. 4. Separately estimate whether drag on a dense bob over one swing can change the period by anything close to 2 % (it changes amplitude, and only through amplitude the period, by a second-order amount). 5. State the verdict and the section the attribution belongs in.

## 4.
1. Copy `projects/_template/report-template.md`. 2. Fill each section with one or two lines. 3. Introduce the three mistakes. 4. Run the linter, keep the output. 5. Fix and rerun. 6. Choose a fourth mistake from: a unit not in the list, a bare number in the Summary, a figure with a caption but no reference in the prose, a Prediction written after the measurement.

## 5.
1. Write the test for a referenced figure (expect `[]`) and for an unreferenced one (expect `[n]`). 2. Collect caption line numbers from `find_figures`. 3. Search all other lines for `Figure n` / `Fig. n` and collect the integers. 4. Return captioned numbers not in that set, in document order. 5. Run pytest.

## 6.
1. Keep all nine headings. 2. Make the Prediction and Result disagree by many combined uncertainties. 3. Write a Conclusion that claims agreement, or that claims something the data cannot address. 4. Explain which comparison a human makes that the code does not, and whether a numeric comparison between the Prediction and the Results sections could be automated (partly, if both are written in the same `value ± u unit` form).

## 7.
1. Go section by section. 2. For each, ask the reviewer's question from the lesson's table and answer it for a parameter sweep. 3. Note that Prediction still exists: what trend did you expect before running? 4. Note that Results, Discussion, Reproducibility are unchanged in meaning. 5. Fill "Apparatus" with environment, versions, seeds, solver settings; "Calibration" with a verification case; "Safety" with "none beyond computer use".
