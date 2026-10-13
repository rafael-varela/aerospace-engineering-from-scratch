# P00.11  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units and uncertainty on every numerical answer.

## 1. (Warm-up: decoding designations)

Decode NACA 2415, 0009, 6412, and 4418 into maximum camber, camber position, and maximum thickness as fractions of chord. Then, for a wing of chord 1.20 m, give the maximum camber and maximum thickness in millimetres. Answer form: a table with 4 rows and 5 numeric columns.

## 2. (Derivation: the rear camber arc)

Starting from a general parabola $y_c = a x^2 + b x + c$ on $p \le x \le 1$ and the three conditions $y_c(p) = m$, $y_c'(p) = 0$, $y_c(1) = 0$, derive equation (2) of the lesson. Then show that the slope of the front arc (equation (1)) and of the rear arc are both zero at $x = p$, and compute both slopes at $x = 0$ and $x = 1$ for the NACA 2412. Answer form: derivation plus two slope values (dimensionless).

## 3. (Unfamiliar scenario: a sailplane polar)

A sailplane section's drag polar is tabulated below (whole-aircraft values are different; these are section values).

| $C_l$ | 0.0 | 0.2 | 0.4 | 0.6 | 0.8 | 1.0 | 1.2 | 1.4 |
|---|---|---|---|---|---|---|---|---|
| $C_d$ | 0.0080 | 0.0068 | 0.0066 | 0.0072 | 0.0088 | 0.0115 | 0.0160 | 0.0240 |

Find $C_{d,\min}$ and the $C_l$ where it occurs, and the maximum $L/D$ and the $C_l$ where it occurs, first by the tangent construction on a sketch and then with `drag_polar_summary`. Explain in two sentences why the maximum L/D does not occur at minimum drag. Answer form: four numbers, a sketch, two sentences.

## 4. (Numerical, uses your code and the report)

Digitise about eight points from the NACA 2412 lift curve at Re ≈ 6 million (smooth) in Report 824, each with a reading uncertainty $u_{C_l}$. Fit with `lift_curve_fit`. Then refit twice with every $C_l$ shifted by $+u_{C_l}$ and by $-u_{C_l}$, and once with the first and last points removed. Report the slope and zero-lift angle with an uncertainty that covers all four fits, and compare the slope with $2\pi$ per radian. Answer form: two numbers with uncertainties and units, a one-sentence comparison, and the CSV.

## 5. (Numerical, uses your code: Reynolds-number effect)

From the same plot, read $C_{l,\max}$ at the lowest and highest Reynolds numbers shown, with reading uncertainties. Does $C_{l,\max}$ change by more than the combined uncertainty? Then read the lift-curve slope at both Reynolds numbers with `lift_curve_fit`. Which of the two quantities is sensitive to Reynolds number, and what does that tell you about where viscosity matters on the section? Answer form: four numbers with uncertainties and a three-sentence argument.

## 6. (Judgment: does the data apply?)

A hobbyist plans a model glider with a 0.30 m chord flying at 8 m/s at sea level (take $\rho = 1.225$ kg/m³, $\mu = 1.8 \times 10^{-5}$ Pa·s) and wants to use the NACA 2412 data from Report 824 for its stall speed. Compute the model's Reynolds number, compare it with the report's range, and decide whether the stall prediction can be trusted. State which quantity from the report (slope, $\alpha_{L0}$, $C_{l,\max}$, $C_{d,\min}$) you would trust most and least at the model's Reynolds number, with one reason each. Answer form: a Reynolds number, a verdict, and a ranked list.

## 7. (Argument: summary versus abstract)

Write a five-sentence summary of Report 824 in the form the `P00.10` report template demands (question, method, headline result with a quantitative statement, agreement with expectation, main limitation). Then write the two-sentence abstract that NTRS might show. Explain in one sentence what makes them different. Answer form: 5 + 2 + 1 sentences.
