# P01.09  Worked solutions

Worked for seed $s = 15$ where a seed is used. Every number below is reproduced by `solution.py` or by the short scripts quoted inline.

## 1. Minimum drag at fixed weight

$W = 13\,500$ N, $S = 14$ m², $C_{D0} = 0.028$, $e = 0.8$, AR = 8.

(a) $K = 1/(\pi \times 0.8 \times 8) = 0.049736$. $V_{md} = \sqrt{2 \times 13\,500/(1.225 \times 14)}\,(0.049736/0.028)^{1/4} = 39.68 \times 1.1544 = 45.807$ m/s. $D_{\min} = 2 \times 13\,500\sqrt{0.049736 \times 0.028} = 1007.58$ N. $(L/D)_{\max} = 13\,500/1007.58 = 13.40$.

(b) $C_L = 2 \times 13\,500/(1.225 \times 45.807^2 \times 14) = 0.7503$; $\sqrt{C_{D0}/K} = \sqrt{0.028/0.049736} = 0.7503$ ✓.

(c) $V = 54.968$ m/s, $q = 1850.6$ Pa. Parasite $= 1850.6 \times 14 \times 0.028 = 725.45$ N (which is $1.44 \times 503.79$, half of $D_{\min}$ ✓); induced $= 0.049736 \times 13\,500^2/(1850.6 \times 14) = 349.85$ N ($= 503.79/1.44$ ✓). $D = 1075.31$ N. $\partial D/\partial V = 1.225 \times 54.968 \times 14 \times 0.028 - 4 \times 0.049736 \times 13\,500^2/(1.225 \times 54.968^3 \times 14) = 26.397 - 12.731 = 13.666$ N/(m/s). $\partial D/\partial\mathrm{AR} = K'W^2/(qS)$ with $K' = -K/8 = -0.0062170$: $-0.0062170 \times 13\,500^2/(1850.6 \times 14) = -43.73$ N per unit AR. So $\nabla D = (13.666, -43.732)$.

(d) $\hat{\mathbf u} = (1, 0.1)/1.00499 = (0.99504, 0.099504)$; $D_{\hat{\mathbf u}} = 13.666 \times 0.99504 - 43.732 \times 0.099504 = 9.247$ N per unit step. Largest available: $|\nabla D| = \sqrt{13.666^2 + 43.732^2} = 45.82$, along $(0.298, -0.954)$, i.e. mostly toward *lower* AR, which is the direction of fastest *increase*; the fastest decrease is the opposite way.

Code: `ac = Aircraft(S=14, cd0=0.028, e=0.8, W0=13500, c_w=0.0)`; `minimum_drag_speed(13500, 8, ac)` → 45.8066; `drag(54.9679, 8, ac)` → 1075.31; `drag_gradient(54.9679, 8, ac)` → `[13.6662, -43.7316]`; `directional_derivative(g, [1, 0.1])` → 9.2470.

## 2. The drag landscape in closed form

(a) With AR fixed, $K$ and $W$ are constants and $D = \tfrac12\rho SC_{D0}V^2 + (2KW^2/\rho S)V^{-2}$, so $\partial D/\partial V = \rho SC_{D0}V - 4KW^2/(\rho SV^3)$, which is (2). With $V$ fixed, $q$ is a constant and $D = qSC_{D0} + (KW^2)/(qS)$; $\partial D/\partial\mathrm{AR} = (KW^2)'/(qS) = (K'W^2 + 2KWW')/(qS)$ with $K' = -K/\mathrm{AR}$ and $W' = 2c_w\mathrm{AR}$, which is (3).

(b) $\partial D/\partial V = 0 \iff \rho SC_{D0}V = 4KW^2/(\rho SV^3) \iff \tfrac12\rho V^2SC_{D0} = 2KW^2/(\rho V^2S) \iff qSC_{D0} = KW^2/(qS)$: parasite = induced. Each equals $\sqrt{(qSC_{D0})(KW^2/qS)} = W\sqrt{KC_{D0}}$, so $D_{\min} = 2W\sqrt{KC_{D0}}$.

(c) $D_{\min}(\mathrm{AR}) = 2\sqrt{C_{D0}/\pi e}\;(W_0 + c_w\mathrm{AR}^2)\mathrm{AR}^{-1/2} = P\,(W_0\mathrm{AR}^{-1/2} + c_w\mathrm{AR}^{3/2})$ with $P = 2\sqrt{C_{D0}/\pi e}$. $D_{\min}' = P(-\tfrac12W_0\mathrm{AR}^{-3/2} + \tfrac32c_w\mathrm{AR}^{1/2}) = 0 \iff 3c_w\mathrm{AR}^2 = W_0 \iff \mathrm{AR}^* = \sqrt{W_0/3c_w}$. Wing penalty $c_w\mathrm{AR}^{*2} = W_0/3$ ✓.

(d) $D = AV^2 + B/V^2$, $A = \tfrac12\rho SC_{D0}$, $B = 2KW^2/(\rho S)$. $D_{VV} = 2A + 6B/V^4$; at $V_{md}$, $2AV = 2B/V^3 \Rightarrow B = AV^4$, so $D_{VV} = 8A = 4\rho SC_{D0} = 4 \times 1.225 \times 16 \times 0.025 = 1.960$ N/(m/s)². $D_{V,\mathrm{AR}} = \partial(2AV - 2B/V^3)/\partial\mathrm{AR} = -2B'/V^3$, and $B' = (2/\rho S)(KW^2)'$; but $(KW^2)' = 0$ at AR* is exactly $D_{\min}' = 0$ (since $D_{\min} \propto W\sqrt K \propto \sqrt{KW^2}$), so $D_{V,\mathrm{AR}} = 0$. $D_{\mathrm{AR,AR}}$: second derivative of $P(W_0\mathrm{AR}^{-1/2} + c_w\mathrm{AR}^{3/2})$ is $P(\tfrac34W_0\mathrm{AR}^{-5/2} + \tfrac34c_w\mathrm{AR}^{-1/2})$; with $W_0 = 3c_w\mathrm{AR}^{*2}$ this is $P \cdot 3c_w\mathrm{AR}^{*-1/2} = 6c_w\sqrt{C_{D0}/(\pi e\,\mathrm{AR}^*)} = 6c_w\sqrt{K^*C_{D0}}$. Numbers: $K^* = 1/(\pi \times 0.8 \times 10.4083) = 0.038225$; $6 \times 40 \times \sqrt{0.038225 \times 0.025} = 7.4194$ N. $H^* = \mathrm{diag}(1.960, 7.419)$; both eigenvalues positive; determinant $14.54 > 0$ and $D_{VV} > 0$: minimum. (The full second partial $\partial^2D/\partial\mathrm{AR}^2$ at fixed $V$ differs from $D_{\min}''$ in general, but at the optimum they coincide because $D_{V,\mathrm{AR}} = 0$ there; `hessian_fd` confirms 7.4194.)

(e) $D_V = \rho VS(C_{D0} + KC_L^2)$, $L_V = \rho VSC_L$; $D_{C_L} = 2qSKC_L$, $L_{C_L} = qS$. $\nabla D = \lambda\nabla L$: the $C_L$ equation gives $\lambda = 2KC_L$; the $V$ equation then reads $C_{D0} + KC_L^2 = 2KC_L^2$, so $C_L^* = \sqrt{C_{D0}/K}$ and $\lambda^* = 2K\sqrt{C_{D0}/K} = 2\sqrt{KC_{D0}}$. $(L/D)_{\max} = W/D_{\min} = 1/(2\sqrt{KC_{D0}})$, so $\lambda^* = 1/(L/D)_{\max}$. Directly: $\mathrm{d}(2W\sqrt{KC_{D0}})/\mathrm{d}W = 2\sqrt{KC_{D0}}$ ✓. For the lesson's aircraft, $\lambda^* = 0.061829$; `lagrange_min_drag(17333.3)` → `(46.767, 0.80869, 0.061829)`.

## 3. A glider hunting a thermal

$R^2 = 1.6 \times 10^5$ m², $r = \sqrt{300^2 + 200^2} = 360.6$ m.

(a) $w = 3\exp(-130\,000/160\,000) = 3 \times 0.44375 = 1.3312$ m/s. $\nabla w = -(2/R^2)w(x, y) = -1.6640 \times 10^{-5} \times (300, 200) = (-4.992 \times 10^{-3}, -3.328 \times 10^{-3})$ (m/s)/m $= $ 1/s. $|\nabla w| = 6.000 \times 10^{-3}$ 1/s. Direction: toward the core, bearing $\mathrm{atan2}(-300, -200) = 236.3°$ from north (south-west), at 6.0 mm/s of updraft per metre.

(b) $\hat{\mathbf u} = (1, 1)/\sqrt2$: $D_{\hat{\mathbf u}}w = (-4.992 - 3.328) \times 10^{-3}/\sqrt2 = -5.883 \times 10^{-3}$ 1/s. Flying north-east, the lift weakens by 5.9 mm/s per metre; the glider is leaving the thermal.

(c) Steepest ascent is bearing 236.3°; the glider is on 045°, so the turn is 191.3°, i.e. essentially reverse course (168.7° the other way is equivalent; either way a near-180° turn). For a radial field $w(r)$, $\nabla w = w'(r)\hat{\mathbf r}$ is always along the line to the core, so "toward the core" and "steepest ascent" coincide. For an elliptical thermal $w = w_0\exp[-(x^2/a^2 + y^2/b^2)]$ with $a \ne b$, $\nabla w \propto -(x/a^2, y/b^2)$, which is not parallel to $(x, y)$: the steepest-ascent direction leans toward the short axis.

(d) $w_{xx} = (-2/R^2 + 4x^2/R^4)w$, $w_{yy} = (-2/R^2 + 4y^2/R^4)w$, $w_{xy} = (4xy/R^4)w$. At the core: $H = -(2w_0/R^2)I = -3.75 \times 10^{-5}\,I$ (1/(m·s)), both eigenvalues negative: maximum ✓. At the glider: $w_{xx} = (-1.25 \times 10^{-5} + 1.4063 \times 10^{-5}) \times 1.3312 = 2.08 \times 10^{-6}$, $w_{yy} = (-1.25 \times 10^{-5} + 6.25 \times 10^{-6}) \times 1.3312 = -8.32 \times 10^{-6}$, $w_{xy} = 9.375 \times 10^{-6} \times 1.3312 = 1.248 \times 10^{-5}$; eigenvalues $-1.664 \times 10^{-5}$ and $+1.040 \times 10^{-5}$. Not a critical point ($\nabla w \ne 0$), so "saddle" does not apply; the mixed signs mean the surface curves down along the tangential direction (the ring $r = $ const is a crest in that direction) and *up* along the radial direction because $r = 361$ m is beyond the inflection radius $R/\sqrt2 = 283$ m, where the Gaussian's tail is convex. The Hessian describes curvature everywhere; it classifies only where the gradient vanishes.

(e) Circling at radius $\varrho$ about a centre $\mathbf c$, the variometer reads $w(\mathbf c + \varrho\hat{\mathbf u}(\theta)) \approx w(\mathbf c) + \varrho\,\nabla w \cdot \hat{\mathbf u}(\theta) = w(\mathbf c) + \varrho|\nabla w|\cos(\theta - \theta_g)$. The reading varies sinusoidally around the circle; its maximum occurs on the heading $\theta_g$ of the gradient and its peak-to-peak amplitude is $2\varrho|\nabla w|$. Each point on the circle is a directional-derivative sample along $\hat{\mathbf u}(\theta)$; one full circle samples every direction, and the pilot shifts the circle toward the heading where the vario peaked. That is gradient ascent with a scalar sensor and a 20 m probe.

## 4. Step size, scaling and conditioning

(a) Raw variables from $(35, 6)$, `tol=1e-6`:

| $s$ | $10^{-4}$ | $10^{-3}$ | $10^{-2}$ | 0.03 | 0.1 | 0.25 | 0.5 |
|---|---|---|---|---|---|---|---|
| iterations | > 50 000 (not converged) | 6943 | 752 | 259 | 74 | 93 | 93 |

(b) Scaled, `tol=1e-9`:

| $s$ | 0.05 | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|---|
| iterations | 385 | 197 | 97 | 64 | 323 | 46 |

(c) Raw Hessian at the optimum: diag(1.960, 7.419), eigenvalues 1.960 and 7.419, condition number 3.79, limit $2/\lambda_{\max} = 0.270$. Scaled: diag(4.900, 0.7419) (the $V$ entry is multiplied by $50^2/1000 = 2.5$, the AR entry by $10^2/1000 = 0.1$), eigenvalues 0.742 and 4.900, condition number 6.60, limit 0.408. For small $s$ each step shrinks the error along the slow eigen-direction by $(1 - s\lambda_{\min})$, so the count is about $\ln(e_0/\mathrm{tol})/(s\lambda_{\min}) \propto 1/s$: $10^{-3} \to 10^{-2}$ divides the count by 9.2 ✓. At $s = 0.25$, just under the raw limit, the stiff direction's factor $|1 - 0.25 \times 7.42| = 0.855$ is no longer small and the count rises. Above the limit a pure fixed-step iteration would diverge; the code's backtracking halves $s$ whenever $f$ does not decrease, so $s = 0.5$ behaves exactly like $s = 0.25$ (93 iterations) and scaled $s = 0.5$ behaves like 0.25 in the stiff direction while keeping the longer step where it helps (46). The scaled $s = 0.4$ is the worst case: just under the limit, the stiff direction barely contracts and nothing is halved.

(d) Scaling was worth doing because the gradient components became comparable and the step dimensionless: with raw variables you cannot know whether $s = 10^{-3}$ is small or large without computing a Hessian in mixed units, and in ten variables you cannot plot it. Conditioning is a separate property. The scaled Hessian is $\mathrm{diag}(1.96V_{\mathrm{ref}}^2, 7.419\mathrm{AR}_{\mathrm{ref}}^2)/D_{\mathrm{ref}}$; equal entries need $\mathrm{AR}_{\mathrm{ref}}/V_{\mathrm{ref}} = \sqrt{1.96/7.419} = 0.514$, e.g. $V_{\mathrm{ref}} = 50$, $\mathrm{AR}_{\mathrm{ref}} = 25.7$. With `scaled_drag(ac, V_ref=50, AR_ref=25.7)` the Hessian is $4.900\,I$ to three figures, condition number 1.00, and steepest descent with $s = 1/4.9 = 0.204$ converges in a handful of steps because the problem is locally a circular bowl. This is Newton's method in disguise: scaling by the Hessian's square root is what the Newton step does implicitly.

(e) Newton (scaled) from $(0.7, 0.6)$: $|\nabla f| = 2.69, 0.876, 0.211, 3.07 \times 10^{-2}, 1.11 \times 10^{-3}, 1.53 \times 10^{-6}, 2.9 \times 10^{-12}, 4 \times 10^{-16}$. Ratios $e_{k+1}/e_k^2$ over the last three usable pairs: $1.53 \times 10^{-6}/(1.11 \times 10^{-3})^2 = 1.24$, $2.9 \times 10^{-12}/(1.53 \times 10^{-6})^2 = 1.24$: constant, so the order is 2 (quadratic). The first three steps are far from the optimum and converge only linearly; quadratic convergence is a local property.

## 5. Is the optimum worth having?

(a) $\mathrm{AR}^* = \sqrt{W_0/3c_w} \propto c_w^{-1/2}$, independent of $e$: $\delta\mathrm{AR}^*/\mathrm{AR}^* = \tfrac12 \times 20\% = 10\%$, so $\mathrm{AR}^* = 10.41 \pm 1.0$ (code: $c_w = 32 \to 11.64$, $48 \to 9.50$). $D^* = 2W^*\sqrt{K^*C_{D0}}$ with $W^* = \tfrac43W_0$ fixed and $K^* = \sqrt{3c_w/W_0}/(\pi e)$, so $D^* \propto c_w^{1/4}e^{-1/2}$: $\delta D^*/D^* = \sqrt{(5\%)^2 + (5\%)^2} = 7.1\%$, $D^* = 1072 \pm 76$ N (code: $c_w = 32 \to 1013.5$, $48 \to 1121.7$; $e = 0.72 \to 1129.7$, $0.88 \to 1021.8$).

(b) $D_{\min}(8) = 1097.3$ N ($+25.7$), $D_{\min}(13) = 1093.2$ N ($+21.5$). A ±2.5 swing in AR costs 2.4 % of drag, a third of the ±7 % uncertainty in $D^*$ itself. There is no engineering content in the second decimal of $\mathrm{AR}^*$; there is barely any in the first. The useful statement is "AR between 9 and 12, drag within 1 % of the minimum".

(c) $C(\mathrm{AR}) = 15D_{\min}(\mathrm{AR}) + 4000\max(\mathrm{AR} - 9, 0)$ €. At AR = 9 the drag term's slope is $15 \times (-12.12) = -182$ €/AR, far smaller than the $+4000$ €/AR cost slope that switches on above 9. So the cost is decreasing up to 9 and increasing beyond: the optimum is at the kink, $\mathrm{AR} = 9.00$, $C = 16\,197$ €, against $21\,709$ € at $\mathrm{AR}^*$. `minimize_scalar` with bounds (9, 14) returns 9.0001. No Newton step is needed; the derivative does not vanish at the answer, which is the general rule for a kinked objective.

(d) *Recommend AR 9: the drag penalty against the unconstrained optimum is 12 N (1.1 %), while each unit of AR above 9 costs €4000 against a lifetime fuel saving of under €200. The two numbers that drive the choice are the shallowness of the drag valley ($D_{\mathrm{AR,AR}} = 7.4$ N per AR²) and the cost slope of the wing. The uncertainty that matters most is the wing-weight coefficient $c_w$ (±20 %), which moves $\mathrm{AR}^*$ by ±1 and $D^*$ by ±5 %; the span-efficiency $e$ matters as much for drag but not for the choice of AR. A weight breakdown of a built wing of the same construction, or a finite-element sizing of this one, would halve the $c_w$ uncertainty.*

## 6. A hangar door as a constraint

(a) $\mathrm{AR} = 6.25$: $W = 13\,000 + 40 \times 39.06 = 14\,562.5$ N, $K = 0.063662$, $V_{md} = \sqrt{2 \times 14\,562.5/(1.225 \times 16)}(0.063662/0.025)^{1/4} = 38.54 \times 1.2632 = 48.70$ m/s, $D_{\min} = 2 \times 14\,562.5\sqrt{0.063662 \times 0.025} = 1161.9$ N. Door penalty $1161.9 - 1071.7 = 90.2$ N (8.4 %).

(b) Central difference of `drag(minimum_drag_speed(weight(AR), AR), AR)` at 6.25 with $h = 10^{-4}$: $-53.06$ N per unit AR. Closed form: $P = 2\sqrt{0.025/(\pi \times 0.8)} = 0.19947$; $P(-\tfrac12 \times 13\,000 \times 6.25^{-1.5} + 1.5 \times 40 \times 6.25^{0.5}) = 0.19947 \times (-416.0 + 150.0) = -53.06$ ✓. $\mathrm{d}\mathrm{AR}/\mathrm{d}b = 2b/S = 1.25$ per metre, so $\mathrm{d}D_{\min}/\mathrm{d}b = -66.3$ N per metre of span.

(c) $\mathcal L = D_{\min}(\mathrm{AR}) + \mu\,g(\mathrm{AR})$, $g = \mathrm{AR} - 6.25$. Stationarity: $D_{\min}'(\mathrm{AR}) + \mu = 0 \Rightarrow \mu = -D_{\min}'(6.25) = +53.06$ N per unit AR. With the constraint active, $\mu > 0$ means the objective would decrease if the constraint were relaxed ($g$ allowed to go positive): the door is costing drag, 53 N per unit of AR or 66 N per metre. A negative $\mu$ would mean the constraint is not binding (the unconstrained optimum lies inside the feasible set) and should be dropped. The first-order condition for an inequality constraint is exactly this sign check (the Karush–Kuhn–Tucker condition, which Phase 03 uses in earnest).

(d) Linear estimate with the two multipliers at the constrained point: $\lambda_W = 2\sqrt{K(6.25)C_{D0}} = 2\sqrt{0.063662 \times 0.025} = 0.07979$ N/N, $\lambda_b = -66.32$ N/m. $\Delta D \approx 0.07979 \times 150 + (-66.32) \times 1.0 = 11.97 - 66.32 = -54.4$ N: the folding tips pay. Exact: `Aircraft(W0=13150)` at $\mathrm{AR} = 121/16 = 7.5625$: $D_{\min} = 1119.77$ N, change $-42.1$ N. The linear estimate overstates the saving by 12 N because $D_{\min}(\mathrm{AR})$ is convex ($D_{\mathrm{AR,AR}} > 0$): the slope flattens as AR grows, and a 1.3-unit step in AR is not small compared with the distance to the optimum. Both say yes; the multipliers gave the answer and its sign with two multiplications, the exact calculation gave the magnitude to keep.
