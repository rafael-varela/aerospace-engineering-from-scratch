# P01.17  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units (coefficients are dimensionless; say so) and a standard error on every fitted number.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: a lift curve by hand)

Three tunnel points on the lift curve: α = 0°, 4°, 8° with C_L = 0.200 + 0.001s, 0.570 + 0.001s, 0.980 + 0.001s.

(a) Write the design matrix for C_L = c₀ + c₁α, form AᵀA and Aᵀ**b** by hand, and solve the 2×2 normal equations.
(b) Compute the three residuals and verify that they sum to zero and that Σ α_i r_i = 0.
(c) Check with `least_squares_qr` and with `np.linalg.lstsq`. Convert c₁ to per radian and compare with the thin-airfoil value 2π.

Answer form: two coefficients, three residuals with the two orthogonality checks, one comparison.

## 2. (Derivation: projections, normal equations, conditioning)

(a) Starting from |A**x** − **b**|², differentiate with respect to **x** and derive the normal equations. Then re-derive them in one line from the requirement that the residual be orthogonal to every column.
(b) Show that P = A(AᵀA)⁻¹Aᵀ satisfies Pᵀ = P, P² = P and PA = A, and that I − P projects onto the orthogonal complement.
(c) Show that if A = QR with QᵀQ = I and R invertible, the normal-equation solution equals R⁻¹Qᵀ**b**.
(d) Using the fact that the eigenvalues of AᵀA are the squares of the singular values of A (Lesson 18 proves it; take it as given here), show cond₂(AᵀA) = cond₂(A)². For the Läuchli matrix with ε = 10⁻⁸, compute AᵀA in exact arithmetic and in double precision and explain what the solver sees.

Answer form: four short derivations.

## 3. (Unfamiliar scenario: a cambered wing)

A cambered wing has its minimum drag at a non-zero lift coefficient: C_D = C_D,min + K(C_L − C_L0)². Simulate a sweep with `wind_tunnel_data`'s α grid, C_L = 0.095(α + 2°), C_D,min = 0.022, K = `TRUE_K`, C_L0 = s/100, noise 0.0008 in C_D and 0.005 in C_L (use `np.random.default_rng(300 + s)`).

(a) Fit the *parabolic* polar (two parameters) and the *general* polar c₀ + c₁C_L + c₂C_L² (three parameters) with your code. Report both residual rms values.
(b) Show algebraically that C_L0 = −c₁/(2c₂), C_D,min = c₀ − c₁²/(4c₂) and K = c₂, and recover the three physical parameters with standard errors (propagate from the covariance of c₁, c₂ at first order, or by Monte Carlo).
(c) Look at the signs of the parabolic-fit residuals in order of C_L. What do they tell you that the rms alone does not?
(d) The parabolic fit returns a value of K. Is it biased, and in which direction? Explain with the algebra of (b).

Answer form: two rms values, three parameters with uncertainties, a sign pattern with a sentence, one explanation.

## 4. (Numerical, uses your code: how many points?)

(a) For the default data, compute e and its standard error from the fit. Then regenerate the data with 34, 68 and 136 points over the same α range and tabulate the standard error of e against the point count. Fit the exponent of the power law.
(b) Run `monte_carlo_k_scatter` with 400 trials, seed s, for C_L noise 0 and 0.005. Report the ratio of the sample scatter to the formula's standard error in each case.
(c) Estimate the extra C_D-equivalent noise caused by the C_L noise at C_L = 0.8 and at the maximum C_L of the sweep, and show that it accounts for the residual rms of the default fit.
(d) How many points are needed for a standard error of 0.005 in e? Give the number and the assumption it rests on.

Answer form: a table and an exponent, two ratios, two noise estimates, one point count.

## 5. (Judgment: when the fit lies)

(a) The sweep is extended to α = 18°, past stall. Describe what happens to the fitted C_D0 and K and to the residual pattern, and say how you would detect it from the residuals alone.
(b) The balance's C_D uncertainty is 0.0008 at low α but 0.002 above α = 8° (unsteady flow). Which assumption of the standard-error formula breaks, and what replaces ordinary least squares?
(c) A colleague fits a tenth-degree polynomial in α (in degrees) with the normal equations and reports coefficients to eight figures. Using `condition_squaring` on `np.vander(alpha, 11, increasing=True)`, say how many of those figures can be trusted and what you would change.
(d) Name one quantity in this lesson that the simulated data cannot validate, because the truth was built in.

Answer form: four short paragraphs, each with a number.

## 6. (Hard, uses your code: repairing Gram–Schmidt)

"Twice is enough": classical Gram–Schmidt applied a second time to its own Q (CGS2) restores orthogonality.

(a) Implement `cgs2(A)` as Q₁, R₁ = CGS(A); Q₂, R₂ = CGS(Q₁); return Q₂ and R₂R₁. Show A = Q₂(R₂R₁).
(b) On the Läuchli matrix for ε = 10⁻⁵, 10⁻⁷, 10⁻⁸, tabulate `orthogonality_loss` for CGS, MGS and CGS2, alongside ε_mach cond(A).
(c) Explain why the MGS loss tracks ε_mach cond(A) and the CGS loss tracks ε_mach cond(A)².
(d) Solve the ε = 10⁻⁸ least-squares problem with CGS2 and an explicit Qᵀ**b**. Does the re-orthogonalisation rescue the explicit route? Explain the result with the argument of lesson 2.5, and say what that implies about *why* the single-pass MGS explicit route failed.

Answer form: code, a table, two explanations.
