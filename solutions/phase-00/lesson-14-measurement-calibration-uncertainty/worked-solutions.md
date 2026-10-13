# P00.14  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty, at least one independent check, and the common mistakes. Standard uncertainties unless marked k = 2.

## 1. One point of a budget

**Assumptions.** The eight readings are independent; resolution and tolerance are rectangular; the three sources are independent.

**Derivation.** Mean x̄ = 159.8 / 8 = 19.975 g. Deviations: −0.075, 0.125, 0.025, −0.175, 0.025, 0.125, −0.075, 0.025 g; squares sum to 0.0750 g²; s = √(0.0750 / 7) = 0.1035 g; u_A = s / √8 = 0.0366 g. u_B(res) = 0.1 / √12 = 0.0289 g. u_B(tol) = 0.05 / √3 = 0.0289 g.

| Source | Type | u (g) | u² (g²) |
|---|---|---|---|
| repeatability of the mean | A | 0.0366 | 0.00134 |
| display resolution | B | 0.0289 | 0.00083 |
| coin tolerance | B | 0.0289 | 0.00083 |
| **combined** | | **0.0548** | 0.00301 |

U = 2 × 0.0548 = 0.110 g.

**Result.** 19.98 ± 0.11 g (k = 2). The offset at 20 g is x̄ − 20.000 = −0.025 g, which is 0.46 u_c: consistent with zero.

**Check.** The repeatability term is the largest, so more readings would help first; 32 readings would halve it. Dimensional: all terms in grams before squaring.

**Common mistakes.** Using s (0.10 g) instead of s/√n as the Type A term. Dividing the tolerance by 2 as if it were a k = 2 bound rather than a hard limit. Reporting 19.975 ± 0.11 (digits beyond the uncertainty's place).

## 2. The rules from the law

**(a)** y = a ± b: ∂y/∂a = 1, ∂y/∂b = ±1. u_y² = 1²·u_a² + (±1)²·u_b² = u_a² + u_b². The sign disappears on squaring: a difference is as uncertain as a sum.

**(b)** y = ab: ∂y/∂a = b, ∂y/∂b = a. u_y² = b² u_a² + a² u_b². Divide by y² = a²b²: (u_y/y)² = (u_a/a)² + (u_b/b)². For y = a/b: ∂y/∂a = 1/b, ∂y/∂b = −a/b²; u_y² = u_a²/b² + a²u_b²/b⁴; divide by y² = a²/b²: the same result. Relative uncertainties add in quadrature whether you multiply or divide.

**(c)** y = aᵖ: ∂y/∂a = p aᵖ⁻¹. u_y = |p| aᵖ⁻¹ u_a; divide by y = aᵖ: u_y/y = |p| u_a/a.

**(d)** ρ = p/(RT) with R exact: ρ ∝ p¹ T⁻¹, so (u_ρ/ρ)² = (100/101325)² + (0.5/288.15)² = (9.87×10⁻⁴)² + (1.735×10⁻³)² = 9.74×10⁻⁷ + 3.01×10⁻⁶ = 3.98×10⁻⁶; u_ρ/ρ = 1.996×10⁻³. ρ = 101325 / (287.05 × 288.15) = 1.2250 kg/m³; u_ρ = 0.00245 kg/m³.

**Result.** ρ = 1.2250 ± 0.0049 kg/m³ (k = 2). Temperature dominates: its squared term is three times the pressure's. Half a kelvin matters more than a hectopascal.

**Check.** `propagate(lambda p: p[0]/(287.05*p[1]), [101325, 288.15], [100, 0.5])` gives u = 0.002445 kg/m³ with contributions 0.00121 (p) and 0.00213 (T). Agrees.

**Common mistakes.** Forgetting that T must be in kelvin before taking the relative uncertainty (0.5/15 °C would be absurd). Adding relative uncertainties linearly (gives 2.7×10⁻³, too large).

## 3. Wind-tunnel drag coefficient

**Assumptions.** All four inputs independent; the balance reads drag only; v is the free-stream speed at the model.

**Derivation.** C_D = D / (½ ρ v² S) = 0.120 / (0.5 × 1.20 × 144 × 0.0120) = 0.120 / 1.0368 = 0.1157. Relative uncertainties: D 0.0333; ρ 0.00833; v 0.025, entering with exponent 2 → 0.0500; S 0.0167.

| Input | relative u | exponent | term | term² |
|---|---|---|---|---|
| D | 0.0333 | 1 | 0.0333 | 0.00111 |
| ρ | 0.0083 | 1 | 0.0083 | 0.00007 |
| v | 0.0250 | 2 | 0.0500 | 0.00250 |
| S | 0.0167 | 1 | 0.0167 | 0.00028 |
| sum | | | | 0.00396 |

u_C/C = √0.00396 = 0.0629; u_C = 0.00728.

**Result.** C_D = 0.116 ± 0.015 (k = 2). The airspeed contributes 63 % of the variance (0.0025 of 0.00396); improve it first. Halving u_v: the v term becomes 0.025, term² 0.000625, sum 0.00209, u_C/C = 0.0457, so u_C falls to 0.00528: a 27 % reduction from fixing one input. Halving the drag uncertainty instead would reduce u_C by only about 11 %.

**Check.** `propagate` on the same function returns u = 0.007282 and relative contributions [0.0333, 0.0083, 0.0500, 0.0167]. Agrees.

**Common mistakes.** Forgetting the factor 2 on the velocity term (gives 0.047 instead of 0.063). Improving the balance because it is the instrument you built, when the airspeed is what limits you.

## 4. Numerical, with the toolkit

Numbers below are for seed 0; yours differ in the decimals but not in the conclusions.

**(a)** slope 0.98596 ± 0.00087, offset 0.293 ± 0.032 g, s_residual = 0.046 g. |slope − 1| / u = 16: the gain error is unmistakable. |offset| / u = 9.3: the zero error is real too. (Both are what the synthetic data were built with: 0.985 and 0.3.)

**(b)** u_c = √(0.0289² + 0.046²) = 0.055 g. The seven residuals scatter about zero within about ±0.07 g with no bow and no trend: the straight line is adequate.

**(c)** `corrected_value(23.4, 0.05, fit)`: mass 23.436 g, u = 0.063 g, so **23.44 ± 0.13 g (k = 2)**. Contributions: reading 0.051 g, slope 0.021 g, offset 0.032 g. The reading's own uncertainty dominates; the calibration adds about 25 % to the total.

**(d)** `np.polyfit(reference, reading, 1, cov=True)` returns the same slope and offset and the same uncertainties (0.00087, 0.0315) to all printed digits: `polyfit` scales the covariance by the residual variance with n − 2 degrees of freedom, exactly as `calibration_fit` does. With `cov="unscaled"` the uncertainties come out 0.019 and 0.68, a factor 1/s = 21.6 larger, because that option assumes each point has unit variance. Which one is right depends on whether you trust your per-point uncertainties (unscaled, with weights) or want the scatter to tell you (scaled). With no weights given, scaled is correct.

**Common mistakes.** Reading `polyfit`'s coefficients in the wrong order (it returns highest power first). Claiming the slope is "about 1, fine" without dividing the difference by its uncertainty.

## 5. When the straight line is wrong

| Physical cause | Sign in the data | Response |
|---|---|---|
| Non-linear spring or strain gauge, especially near full scale | residuals bow (sign pattern − + − or + − +) against reference | fit a quadratic, or restrict the calibrated range and say so |
| Hysteresis or creep in the mechanism or plastic case | readings on unloading differ from loading at the same reference; residual depends on the previous load | load in one direction only, pause a fixed time before reading, report the loop width as a Type B term |
| Drift: warm-up, battery voltage, room temperature | residuals trend with time or with the order of measurement, not with reference | warm up, re-zero at intervals, randomise order so drift does not masquerade as non-linearity |
| Auto-zero or auto-off firmware | sudden jump in residuals; small loads read exactly 0.0 | keep a tare mass on the pan, read quickly after loading, document the behaviour |
| Overload or off-centre loading | residuals large at the top end only, or dependent on where the mass sat | stay below 80 % of capacity, centre every load |

The measurable sign is always a *structure* in residuals that the error bars cannot explain; a line is "wrong" when the residuals are not random.

**Common mistakes.** Fitting a higher-order polynomial to what is actually drift. Calling any misfit "non-linearity".

## 6. Correlation, and a mistake that was not an uncertainty

**(a)** Readings r₁ = 31.2 g, r₂ = 18.7 g, each with u = 0.05 g; offset z with u = 0.3 g. The difference is y = (r₁ + z) − (r₂ + z) = 12.5 g.
(i) Treating the offset as independent for each reading: u_y² = 0.05² + 0.05² + 0.3² + 0.3² → u_y = 0.43 g.
(ii) Recognising that the same z enters both: the sensitivity coefficients on z are +1 and −1, so its contribution cancels exactly: u_y = √(0.05² + 0.05²) = 0.071 g.
`propagate(lambda p: (p[0]+p[2]) - (p[1]+p[3]), [31.2, 18.7, 0, 0], [0.05, 0.05, 0.3, 0.3], corr=C)` with C having ones at [2,3] and [3,2] returns 0.0707 g; without `corr` it returns 0.430 g. Six times smaller when the correlation is honoured. This is why "tare, then weigh" works, and why a difference measurement can be far better than either absolute reading.

**(b)** An uncertainty is a *spread* of plausible values around a believed value, with a known distribution; the Mars Climate Orbiter's thruster file had none of that. It had a value that was wrong by a definite factor of 4.45 (pound-force-seconds read as newton-seconds), with no spread at all. Every statistical analysis of the thruster data would have found it self-consistent: Type A scatter small, Type B from sensor specifications small, combined uncertainty small, and the result still wrong by a factor 4.45. Uncertainty analysis quantifies how well you know a number *given that the model is right*; it cannot detect that two pieces of software disagree about what the number means. The check that catches it is dimensional, not statistical: carry units through the calculation (pint in `P00.05`) so that lbf·s added to N·s raises a `DimensionalityError`, and verify interface specifications by test, not by reading. The lesson for every report you write: a budget is necessary and not sufficient; the independent sanity check (`P00.12`) is a separate step.

**Common mistakes.** Treating the two offsets as independent "to be safe"; it is not safe, it is wrong, and it hides a real advantage of differential measurement.

## 7. When k = 2 is not 95 %

**(a)** The quantity (x̄ − μ)/(s/√n) follows Student's t with ν = n − 1 degrees of freedom, because s is itself uncertain when n is small. For ν = 2 the 97.5 % quantile is 4.30: the 95 % interval is ±4.30 u_A, more than twice what k = 2 claims. Reporting k = 2 with n = 3 overstates your confidence by a factor of two.

**(b)** t₀.₉₇₅: ν = 3 → 3.18, 4 → 2.78, 5 → 2.57, 9 → 2.26, 19 → 2.09, 29 → 2.05. Within 10 % of 2 (that is, below 2.2) needs ν ≥ 11, so about 12 readings; within 5 % needs about 30. The course rule "at least 5, prefer 10" is a compromise; say in the report when a Type A term from few readings dominates.

**(c)** h = 0.30 ± 0.10 m, g = 9.81 m/s². Linear law: t₀ = 0.2473 s, u_t = ½ t₀ (0.10/0.30) = 0.0412 s. Monte Carlo (10⁵ normal samples; about 0.1 % fall at h ≤ 0 and must be discarded, which already says the normal model for h is poor here): standard deviation 0.0432 s, 5 % larger, and the mean of t is 0.2435 s, below t₀. The distribution of t is skewed to the left (skewness −0.58) because √h rises steeply for small h and flattens for large h: the low tail of h stretches out in t. The linear law is a tangent-line approximation and misses the curvature; at a 33 % relative input uncertainty the error is a few percent of the uncertainty, which is acceptable, but the asymmetry and the shift of the mean are things the linear law cannot report at all.

**Common mistakes.** Reading the Monte Carlo standard deviation as "the right answer" when 0.1 % of the samples were physically impossible: the input model was already stretched.

## Check-yourself answers

1. Type A: a standard uncertainty obtained by statistical analysis of repeated observations. Type B: one obtained by any other means (specifications, tolerances, resolution, physics, previous data). Type A is the "statistical" one.
2. 0.5 / √12 = 0.14 g.
3. 0.20 / √5 = 0.089 g; 20 readings (four times as many) halve it.
4. u_y = √(0.3² + 0.4²) = 0.5; the sign does not matter because the sensitivity coefficient is squared.
5. (u_q/q)² = (1 %)² + (2 × 2 %)² → 4.1 %.
6. u rounds to 0.10 (two significant figures), so the value rounds to two decimals: 3.14 ± 0.10.
7. The scale's response is not a straight line over the range used; either fit a curve (quadratic) and show the residuals become random, or restrict the calibration to the range where the line holds, and in both cases report the misfit as a Type B term.
8. No: the zero offset enters both readings identically and cancels in the difference (correlation coefficient +1). The exception is a scale that auto-zeros between the two readings, in which case the two offsets are no longer the same quantity.
