# P00.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result with uncertainty where applicable, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes.

## 1.

**Assumptions.** Standard four-digit decoding; chord 1.20 m exactly.
**Derivation.** $m = d_1/100$, $p = d_2/10$, $t = d_{34}/100$; millimetres $= $ fraction $\times 1200$ mm.
**Result.**

| Section | $m$ | $p$ | $t$ | max camber (mm) | max thickness (mm) |
|---|---|---|---|---|---|
| 2415 | 0.02 | 0.4 | 0.15 | 24 | 180 |
| 0009 | 0 | (undefined) | 0.09 | 0 | 108 |
| 6412 | 0.06 | 0.4 | 0.12 | 72 | 144 |
| 4418 | 0.04 | 0.4 | 0.18 | 48 | 216 |

**Check.** `parse_naca4("6412")` returns `(0.06, 0.4, 0.12)`. Thickness in mm is always larger than camber in mm for these sections, as it should be: a 6 % camber on a 12 % thick section is already a strongly arched shape.
**Common mistakes.** Reading the second digit as hundredths; forgetting that 0009's position digit means nothing; reporting thickness as the half-thickness.

## 2.

**Assumptions.** Parabolic arc; $0 < p < 1$.
**Derivation.** $y_c = ax^2 + bx + c$, $y_c' = 2ax + b$. From $y_c'(p) = 0$: $b = -2ap$. From $y_c(1) = 0$: $a + b + c = 0 \Rightarrow c = -a - b = -a + 2ap = a(2p - 1)$. From $y_c(p) = m$: $ap^2 - 2ap^2 + a(2p-1) = a(-p^2 + 2p - 1) = -a(1-p)^2 = m$, so $a = -m/(1-p)^2$. Then $b = 2mp/(1-p)^2$ and $c = m(1-2p)/(1-p)^2$, giving
$$y_c = \frac{m}{(1-p)^2}\left[(1-2p) + 2px - x^2\right],$$
which is equation (2). Slopes: front $y_c' = \frac{2m}{p^2}(p - x)$, rear $y_c' = \frac{2m}{(1-p)^2}(p - x)$; both are zero at $x = p$. For 2412 ($m = 0.02$, $p = 0.4$): at $x = 0$, $y_c' = 2(0.02)(0.4)/0.16 = 0.100$; at $x = 1$, $y_c' = 2(0.02)(0.4 - 1)/0.36 = -0.0667$.
**Result.** Leading-edge camber slope $+0.100$ (about 5.7°), trailing-edge slope $-0.067$ (about −3.8°). Dimensionless.
**Check.** The rear arc at $x = p$: $(1-2p) + 2p^2 - p^2 = (1-p)^2$, so $y_c = m$. Limiting case $p \to 1$ would make the rear arc degenerate, which is why the designation never uses a 9 or 10 in the second digit for the four-digit series. `camber_slope(np.array([0.0, 1.0]), 0.02, 0.4)` returns the two numbers above.
**Common mistakes.** Dropping the constant $c$ on the rear arc (it is nonzero); sign error at the trailing edge; evaluating the front-arc formula at $x = 1$.

## 3.

**Assumptions.** The table is the polar; linear interpolation between rows is not needed for the asked values.
**Derivation.** Minimum $C_d$ in the table is 0.0066 at $C_l = 0.4$. Ratios $C_l/C_d$: 0, 29.4, 60.6, 83.3, 90.9, 87.0, 75.0, 58.3. Maximum 90.9 at $C_l = 0.8$. On a sketch with $C_d$ horizontal and $C_l$ vertical, a line from the origin rotated down until it first touches the curve touches at (0.0088, 0.8); its slope is 90.9.
**Result.** $C_{d,\min} = 0.0066$ at $C_l = 0.4$; $(L/D)_{\max} = 91$ at $C_l = 0.8$ (reading these from a table, the true maximum lies somewhere between $C_l = 0.6$ and $1.0$, so quote $0.8 \pm 0.2$). `drag_polar_summary` returns the same four numbers.
**Check.** A parabola $C_d = 0.0066 + k(C_l - 0.4)^2$ fitted to the points near the minimum gives $k \approx 0.014$; the tangent point for an offset parabola is at $C_l^* = \sqrt{0.4^2 + 0.0066/0.014} = \sqrt{0.16 + 0.47} = 0.79$, consistent with 0.8.
**Common mistakes.** Reporting maximum L/D at the minimum-drag point; computing $C_d/C_l$; ignoring that tabulated data give the maximum only to the table's resolution. Two sentences: minimum drag is the lowest point of the polar, but moving up the polar from there, lift rises linearly while drag rises only quadratically at first, so the ratio keeps improving until the quadratic term catches up; the tangent from the origin marks where it does.

## 4.

**Assumptions.** Learner-specific readings; the numbers below are representative, not the answer key. Reading uncertainty $u_{C_l} = 0.02$, $u_\alpha = 0.2°$.
**Derivation.** Typical readings for NACA 2412, Re ≈ 6 million, smooth: α = −6, −4, −2, 0, 2, 4, 6, 8° with $C_l$ ≈ −0.42, −0.21, 0.00, 0.22, 0.43, 0.65, 0.86, 1.06. Fit: slope ≈ 0.106/deg, $\alpha_{L0}$ ≈ −2.0°. Shifting all $C_l$ by ±0.02 changes the slope by < 0.0005 but moves $\alpha_{L0}$ by ∓0.19°. Dropping the end points changes the slope by a few thousandths depending on the scatter.
**Result.** Slope $a_0 = 0.106 \pm 0.004$ per degree; $\alpha_{L0} = -2.0 \pm 0.3°$. Thin-airfoil: $2\pi$ rad⁻¹ $= 0.1097$ deg⁻¹; the measured slope is about 3 % lower, as expected for a real section of finite thickness with a boundary layer.
**Check.** The slope must lie between 0.09 and 0.115 per degree; the 0012 fitted the same way must give $\alpha_{L0} = 0.0 \pm 0.3°$.
**Common mistakes.** Including stall-region points; reading the "standard roughness" curve; forgetting to convert $2\pi$ to per degree; quoting the slope to five figures when the reading uncertainty allows three.

## 5.

**Assumptions.** Reading uncertainty 0.02 on $C_l$; the plot shows Re ≈ 3 and 9 million as the extremes.
**Derivation.** Typical: $C_{l,\max} \approx 1.5 \pm 0.02$ at 3 million, $\approx 1.7 \pm 0.02$ at 9 million; difference 0.2 against a combined uncertainty of 0.03, so the change is real. Slopes: both ≈ 0.105–0.107/deg, difference within ±0.004.
**Result.** $C_{l,\max}$ depends on Reynolds number (rises with it); the lift-curve slope does not, within reading uncertainty.
**Check.** This matches §2.6: the linear region is set by the pressure distribution (nearly inviscid), the stall by boundary-layer separation (viscous). At higher Reynolds number the boundary layer is thinner and transitions to turbulent earlier, which resists separation and delays stall.
**Common mistakes.** Comparing curves at different roughness instead of different Reynolds numbers; claiming a slope difference smaller than the reading uncertainty is significant.

## 6.

**Assumptions.** $\rho = 1.225$ kg/m³, $\mu = 1.8 \times 10^{-5}$ Pa·s, $V = 8$ m/s, $c = 0.30$ m.
**Derivation.** $Re = \rho V c / \mu = 1.225 \times 8 \times 0.30 / 1.8 \times 10^{-5} = 1.63 \times 10^5$. Dimensions: (kg/m³)(m/s)(m)/(Pa·s) = (kg/(m·s))/(kg/(m·s)) = 1. The report's lowest curve is 3 million, about 18 times higher.
**Result.** $Re \approx 1.6 \times 10^5$, far below the report's range. The stall speed prediction cannot be trusted: $C_{l,\max}$ was shown in problem 5 to depend on Reynolds number even within the report's range, and at $10^5$ the boundary layer may stay laminar and separate early, giving a lower $C_{l,\max}$ (and a stall speed higher than predicted, which is the dangerous direction). Trust ranking: most, the lift-curve slope (pressure-driven, weakly viscous) and $\alpha_{L0}$ (geometric); least, $C_{l,\max}$ and $C_{d,\min}$ (both boundary-layer-driven; drag at low Re is typically two to three times the report's value).
**Check.** Order of magnitude: a full-size light aircraft, $c = 1.5$ m at 50 m/s, gives $Re = 5 \times 10^6$, inside the range. The model is 5× smaller and 6× slower: 30× lower Re. Consistent.
**Common mistakes.** Forgetting that both chord and speed enter; using the report because "it is the same section"; trusting the stall value because it is the one they wanted.

## 7.

**Assumptions.** Written in the learner's own words; a model answer follows.
**Result.** Summary: "This report asks what lift, drag, and pitching-moment characteristics the NACA airfoil families have, and how each family's geometry is defined. Sections were tested as two-dimensional models in the Langley low-turbulence pressure tunnel at Reynolds numbers of about 3, 6, and 9 million, smooth and with standard leading-edge roughness. In the linear range all sections show lift-curve slopes near 0.10 to 0.11 per degree, with minimum drag coefficients of roughly 0.005 to 0.009 for smooth surfaces and about double with roughness. The slopes agree with thin-airfoil theory ($2\pi$ per radian) to within a few percent, as expected. The data are two-dimensional and low-speed; finite-wing, compressibility, and low-Reynolds-number effects are outside their scope." Abstract: "A compilation of the aerodynamic characteristics of NACA airfoil sections from wind-tunnel tests, with a description of the airfoil families and their designation systems. Section lift, drag, and moment data are presented in plotted form for a wide range of sections at several Reynolds numbers." Difference: the summary states what was found, with numbers; the abstract states what the document contains.
**Check.** Each sentence of the summary maps to one item of the `P00.10` template's Summary line.
**Common mistakes.** Writing the summary without a number; writing the abstract as a sales pitch; making either longer than asked.

## Check-yourself answers

1. $m = 0.06$, $p = 0.4$, $t = 0.12$.
2. A symmetric section has zero camber, so there is no maximum-camber position to locate.
3. $2\pi$ rad⁻¹ × (π/180 rad/deg) $= 2\pi^2/180 = 0.1097$ per degree.
4. Positive: at α = 0 the section is 2° above its zero-lift angle, so $C_l \approx 0.11 \times 2 = 0.22$.
5. $Re = \rho V c/\mu$, the ratio of inertial to viscous effects. Nine million means a full-scale chord at flight speed, achieved in the pressure tunnel by raising the air density rather than the model size.
6. The line through the origin tangent to the polar; its slope is $(L/D)_{\max}$ and the tangent point gives the $C_l$ at which it occurs.
7. Equation (3) evaluates to $0.0021 \times 5t$ at $x = 1$, a deliberate feature of the original polynomial; for a plot it is invisible, for a manufactured part the last coefficient is often adjusted to close the edge.
8. The notes. The validator checks that a field has text, not that the text is a result; "see figure" is a placeholder in meaning but not in form, and only a reader catches it.
