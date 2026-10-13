# P00.09  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; symbolic results are exact, numerical ones are given to four significant figures.

## 1. The package

**Assumptions.** Lesson model: constant $\rho$, $C_D$, $A$, $g$; released from rest; falls without tumbling.

**Derivation** ($s = 15$, $m = 2.0$ kg). (a) $k = \tfrac12 \times 1.225 \times 1.0 \times 0.010 = 0.006125$ kg/m. (b) $v_t = \sqrt{2.0 \times 9.80665/0.006125} = 56.59$ m/s. (c) $\tau = v_t/g = 5.770$ s. (d) $t_{99} = 2.6467\,\tau = 15.27$ s. (e) $x_{99} = v_t\tau\ln\cosh(2.6467) = 56.59 \times 5.770 \times 1.9587 = 639.5$ m. General: $v_t = 40.01\sqrt{m/\mathrm{kg}}$ m/s for this $k$.

**Result.** $k = 0.006125$ kg/m, $v_t = 56.59$ m/s, $\tau = 5.770$ s, $t_{99} = 15.27$ s, $x_{99} = 639.5$ m.

**Check.** Units: $\sqrt{\mathrm{kg \cdot m\,s^{-2}/(kg\,m^{-1})}} = \mathrm{m/s}$. Long-time asymptote $x \approx v_t t - v_t\tau\ln 2 = 864.1 - 226.3 = 637.8$ m, within 0.3 % of the exact 639.5 m at $t = 2.65\tau$, as it should be ($e^{-2u} = 0.5\,\%$ there). `lambdify_position()(15.27, 2.0, 9.80665, 0.006125)` returns 639.5.

**Common mistakes.** Forgetting the $\tfrac12$ in $k$ (makes $v_t$ too small by $\sqrt2$); using $\tanh^{-1}$ in degrees mode on a calculator; quoting $x_{99}$ as $v_t t_{99}$ (ignores the acceleration phase, 35 % too high).

## 2. Solving the equation of motion

**Derivation.** (a) Set $v_t^2 = mg/k$, so $g - (k/m)v^2 = (g/v_t^2)(v_t^2 - v^2)$. For $0 \le v < v_t$ the right side is positive: a body released from rest accelerates monotonically towards $v_t$ and never crosses it (at $v = v_t$ the acceleration is zero, so $v_t$ is an equilibrium). Separate variables, which requires $v_t^2 - v^2 \ne 0$, i.e. $v < v_t$:
$$\int_0^v \frac{\mathrm{d}v'}{v_t^2 - v'^2} = \frac{g}{v_t^2}\int_0^t \mathrm{d}t' \;\Rightarrow\; \frac{1}{v_t}\operatorname{artanh}\frac{v}{v_t} = \frac{gt}{v_t^2} \;\Rightarrow\; v = v_t\tanh\frac{gt}{v_t}.$$
(The antiderivative follows from partial fractions: $\frac{1}{v_t^2 - v^2} = \frac{1}{2v_t}\big[\frac{1}{v_t - v} + \frac{1}{v_t + v}\big]$, giving $\frac{1}{2v_t}\ln\frac{v_t + v}{v_t - v}$.)

(b) With $u = gt/v_t$, $\mathrm{d}t = (v_t/g)\,\mathrm{d}u$: $x = \int_0^t v\,\mathrm{d}t' = \frac{v_t^2}{g}\int_0^u \tanh u'\,\mathrm{d}u' = \frac{v_t^2}{g}\ln\cosh u$, since $(\ln\cosh)' = \tanh$ and $\ln\cosh 0 = 0$.

(c) $\cosh u = \tfrac12 e^{u}(1 + e^{-2u})$, so $\ln\cosh u = u - \ln 2 + \ln(1 + e^{-2u})$. For $u \gg 1$ the last term is $\approx e^{-2u} \to 0$ and $x \approx v_t t - (v_t^2/g)\ln 2$. The constant $(v_t^2/g)\ln 2 = 0.693\,v_t\tau$ is the distance "lost" relative to a body that had moved at $v_t$ from the start: the acceleration phase costs about 0.69 time constants' worth of travel.

(d) `sp.dsolve` returns an equivalent expression, typically written with $\sqrt g\sqrt m/\sqrt k$ factors and sometimes as a ratio of exponentials, $v_t\,(e^{2u} - 1)/(e^{2u} + 1)$. `sp.simplify(sol.rhs - velocity_expr())` gives 0; if a form resists `simplify`, `sp.simplify(sp.expand(..., trig=True))` or evaluating both at several numerical points with `sp.N` settles it. `sp.integrate(velocity_expr(), t)` returns the form $v_t\,[t - \tau\ln(\tanh u + 1)]$; it differs from $\ln\cosh$ by the constant $-v_t\tau\ln 2$ only, which `sp.simplify(sp.diff(...) - velocity_expr())` confirms (zero derivative difference) together with checking the value at $t = 0$.

**Check.** Differentiate $v$: $g(1 - \tanh^2 u) = g - (g/v_t^2)v^2$. ✓ Dimensions of $x$: $\mathrm{(m/s)^2/(m/s^2) = m}$. ✓

**Common mistakes.** Dropping the absolute value in $\ln|v_t - v|$ and then worrying about $v > v_t$; forgetting that `sp.integrate` fixes no integration constant; declaring two expressions unequal because `==` returned `False`.

## 3. Thrown upward

**Assumptions.** Upward positive; same constant $k$; $v_0 = 30$ m/s, $m = 1.2$ kg, $k = 0.006125$ kg/m, so $v_t = 43.83$ m/s and $\tau = 4.470$ s.

**Derivation.** (a) Moving up, gravity and drag both point down: $m\,\mathrm{d}v/\mathrm{d}t = -mg - kv^2$. Substituting $v \to -v$ in the lesson's falling equation would give $-mg + kv^2$, i.e. drag pointing up while the body moves up, which is wrong: drag always opposes velocity, so its sign depends on the direction of motion, not on the coordinate choice. (The general form is $-k\,v|v|$.)

(b) $\mathrm{d}v/\mathrm{d}t = -(g/v_t^2)(v_t^2 + v^2)$. Separate: $\int_{v_0}^{v}\mathrm{d}v'/(v_t^2 + v'^2) = -(g/v_t^2)t$, so $\frac{1}{v_t}[\arctan(v/v_t) - \arctan(v_0/v_t)] = -gt/v_t^2$ and $v(t) = v_t\tan\big(\arctan(v_0/v_t) - gt/v_t\big)$. Residual check in SymPy: `sp.simplify(sp.diff(cand, t) + g + k*cand**2/m)` returns 0 (it uses $1 + \tan^2 = \sec^2$).

(c) Apex at $v = 0$: $t_a = \tau\arctan(v_0/v_t) = 4.470 \times \arctan(0.6845) = 2.683$ s. Height: from $v\,\mathrm{d}v/\mathrm{d}x = -(g/v_t^2)(v_t^2 + v^2)$, $h_a = \frac{v_t^2}{2g}\ln\big(1 + v_0^2/v_t^2\big) = 97.95 \times \ln 1.4685 = 37.64$ m. Without drag: $v_0/g = 3.059$ s and $v_0^2/2g = 45.89$ m. Drag shortens the climb by 12 % in time and 18 % in height.

(d) At the apex the body is momentarily at rest and about to fall: that is exactly the initial condition $v(0) = 0$ of the lesson's solution, so the descent is $v_t\tanh(g(t - t_a)/v_t)$ (downward), with no further information needed.

**Check.** $k \to 0$: $v_t \to \infty$, $\arctan(v_0/v_t) \to v_0/v_t$ and $\tan(\epsilon) \approx \epsilon$, so $v \to v_0 - gt$. ✓ The descent from 37.64 m takes longer than the ascent and lands slower than 30 m/s (energy is lost both ways): a good extension check with `lambdify_position`.

**Common mistakes.** Using $\tanh$ for the upward leg; forgetting that $\arctan$ is in radians; computing $h_a$ as $v_t t_a/2$.

## 4. The error curve

**Result** (lesson parameters, $t_0 = 3$ s, $v_0 = 20$ m/s):

| $h$ | velocity check [m/s²] | $\mathrm{d}a/\mathrm{d}v$ check [1/s] | five-point [m/s²] |
|---|---|---|---|
| $10^{-1}$ | $3.14 \times 10^{-5}$ | $8 \times 10^{-16}$ | $5.9 \times 10^{-7}$ |
| $10^{-2}$ | $3.15 \times 10^{-7}$ | $3 \times 10^{-14}$ | $5.9 \times 10^{-11}$ |
| $10^{-3}$ | $3.15 \times 10^{-9}$ | $1 \times 10^{-13}$ | $2.7 \times 10^{-13}$ |
| $10^{-4}$ | $4.6 \times 10^{-11}$ | $3 \times 10^{-13}$ | $2.8 \times 10^{-11}$ |
| $10^{-5}$ | $2.4 \times 10^{-10}$ | | $4.2 \times 10^{-10}$ |
| $10^{-6}$ | $2.8 \times 10^{-9}$ | | |
| $10^{-7}$ | $2.8 \times 10^{-8}$ | | |
| $10^{-8}$ | $4.0 \times 10^{-7}$ | | |

(a) Left-hand slope ($h \ge 10^{-3}$): $2.00 \pm 0.01$. Right-hand slope ($h \le 10^{-5}$): $-1.05 \pm 0.1$ (rounding is noisy; expect scatter of a factor 2 to 3 between runs on different machines).

(b) $v(3) = 25.68$ m/s; $v'''(3) = 0.01893$ m/s⁴ from `sp.diff(velocity_expr(), t, 3)`. Predicted truncation error $h^2|v'''|/6 = 3.15 \times 10^{-3}\,h^2$: at $h = 0.1$ that is $3.15 \times 10^{-5}$, matching the table to three figures. $h_{\text{opt}} = (3 \times 2.2 \times 10^{-16} \times 25.68/0.01893)^{1/3} = 9.6 \times 10^{-5}$: the measured minimum is at $10^{-4}$. The textbook $\varepsilon^{1/3} = 6 \times 10^{-6}$ is a factor 16 too small because it hides the ratio $|v|/|v'''| \approx 1360$.

(c) The five-point stencil has left slope 4 (error $\propto h^4 f^{(5)}/30$), its minimum is $2.7 \times 10^{-13}$ at $h = 10^{-3}$, two orders of magnitude lower than the central difference and at a ten times larger step. The right-hand branch is unchanged: rounding does not care how clever the stencil is.

**Check.** The $\mathrm{d}a/\mathrm{d}v$ column shows no $h^2$ branch at all, because $a$ is quadratic in $v$ and the central difference of a quadratic is exact; what remains is pure rounding growing as $h$ shrinks, confirming the mechanism of the right-hand branch.

**Common mistakes.** Fitting one straight line through all eight points; estimating $v'''$ by finite differences (circular); reporting the right-hand slope to three figures when it is dominated by rounding noise.

## 5. When does drag matter?

**Derivation.** (a) From $v = gt - \frac{g^2k}{3m}t^3 + \dots$: $1 - v/(gt) \approx \frac{gk}{3m}t^2$.

(b) Set equal to 0.01: $t_1 = \sqrt{0.03\,m/(gk)} = \sqrt{0.03 \times 1.2/(9.80665 \times 0.006125)} = 0.7742$ s. Exact deviation at $t_1$: `1 - vel(0.7742, 1.2, 9.80665, 0.006125)/(9.80665*0.7742)` $= 0.988\,\%$. The series overestimates by 1.2 % of itself (the $t^5$ term is positive and reduces the deviation).

(c) Tennis ball: $k = 1.04 \times 10^{-3}$ kg/m, $v_t = 23.37$ m/s, $t_1 = 0.413$ s. Table-tennis ball: $k = 3.86 \times 10^{-4}$ kg/m, $v_t = 8.28$ m/s, $t_1 = 0.146$ s. The free-fall time from 2 m is $\sqrt{2x/g} = 0.6387$ s, longer than both $t_1$, so drag changes the *speed* by more than 1 % in both cases. But the experiment measures *time*: solving $x(t) = 2$ m with the drag model, $t = \tau\operatorname{arcosh}(e^{gx/v_t^2})$, gives 0.6425 s for the tennis ball (+3.8 ms, below the 10 ms resolution: **free fall adequate**) and 0.6695 s for the table-tennis ball (+31 ms, three times the resolution: **free fall not adequate**; the drag model is needed and the experiment can actually detect drag).

**Check.** The relative change in fall time is roughly half the relative change in speed ($x \propto t^2$ at leading order), consistent with $1 - v/gt \approx 2.3\,\%$ for the tennis ball at $t = 0.64$ s giving $\approx 0.6\,\%$ in time.

**Common mistakes.** Deciding from $t_1$ alone without asking what is measured and to what resolution; using the ball's diameter instead of its cross-section area.

## 6. Where the model stops

(a) **Constant density fails.** $\rho(30\ \mathrm{km})/\rho_0 = 0.0147$, so $k$ is 68 times smaller and the local terminal velocity $v_t \propto \rho^{-1/2}$ is $8.3\times$ larger: 362 m/s for the package, which is also supersonic at that altitude (the speed of sound there is about 300 m/s), so $C_D$ changes too. The body never reaches a single $v_t$; it decelerates as it descends into thicker air. With $k = k(x)$ the ODE is $v\,\mathrm{d}v/\mathrm{d}x = g - k(x)v^2/m$, linear in $v^2$ but with no closed form for the ISA's piecewise exponential density; integrate numerically.

(b) **Constant $C_D A$ fails.** $v_t \propto (C_D A)^{-1/2}$: with $C_D A$ between 0.006 and 0.014 m², $v_t$ ranges from 37.0 to 56.6 m/s, a ±20 % band around the nominal 43.8 m/s. A closed form exists for each fixed value, so the honest output is a range (or a time-averaged $C_D A$ if the tumbling is fast compared with $\tau$).

(c) **Quadratic drag fails.** Fog droplet, $r = 10\ \mathrm{\mu m}$: $m = \tfrac43\pi r^3 \times 1000 = 4.19 \times 10^{-12}$ kg, $A = 3.14 \times 10^{-10}$ m². The quadratic model with $C_D = 0.47$ gives $v_t = 0.67$ m/s and $Re = \rho v_t (2r)/\mu = 0.9$: at the edge of the Stokes regime, so the quadratic law is not self-consistent. Stokes' law gives $v_t = 2r^2\rho_w g/(9\mu) = 0.0121$ m/s with $Re = 0.016$, self-consistent. The quadratic model overpredicts by a factor of 56: fog floats. A 2 mm raindrop has $v_t \approx 6.7$ m/s and $Re \approx 900$: quadratic drag is fine. Linear drag has its own closed form (problem 7).

(d) **Constant $g$ fails** by 0.1 % when $(R/(R+h))^2 = 0.999$, i.e. $h = 3.19$ km. A 3 km drop therefore has $g$ varying by 0.1 % and $v_t \propto \sqrt g$ by 0.05 %; the density change over the same 3 km (26 %) is 500 times larger. Variable $g$ is the last assumption to worry about in the atmosphere.

**Common mistakes.** Treating $v_t$ at 30 km as a prediction of impact speed; computing $Re$ with the diameter in millimetres; forgetting that fixing one assumption (density) can break another (Mach number).

## 7. Linear versus quadratic drag

**Derivation.** (a) $m\,\mathrm{d}v/\mathrm{d}t = mg - cv$ is linear; with $v_t = mg/c$ and $\tau = m/c = v_t/g$: $\mathrm{d}v/\mathrm{d}t = (v_t - v)/\tau$, so $v = v_t(1 - e^{-t/\tau})$ for $v(0) = 0$; $x = v_t t - v_t\tau(1 - e^{-t/\tau})$.

(b) Matching $v_t = 43.83$ m/s: $c = mg/v_t = 0.2685$ kg/s. In units of $\tau = v_t/g$ (the same for both models when $v_t$ matches): quadratic $\tanh u = 0.9 \Rightarrow u = \operatorname{artanh}0.9 = 1.472$; linear $1 - e^{-u} = 0.9 \Rightarrow u = \ln 10 = 2.303$. The quadratic model reaches 90 % of terminal speed in 6.58 s, the linear one in 10.29 s. Physically: at the same terminal speed the two drag laws agree only at $v = v_t$; below it, $k v^2 = (mg)(v/v_t)^2 < (mg)(v/v_t) = cv$, so quadratic drag is *weaker* at every sub-terminal speed, the body accelerates harder early, and it closes on $v_t$ sooner.

(c) At $t = 2$ s: free fall 19.61 m/s, quadratic 18.40 m/s, linear 15.81 m/s. The quadratic model is faster by 2.59 m/s (16 %).

(d) The distinguishing signature is the *shape* of $v(t)$ early on: the quadratic model follows free fall closely (within 6 % at 2 s) and then bends sharply, the linear model departs immediately. A video record of a 10 m drop at 60 frames per second gives $v(t)$ from frame-to-frame displacement with a few percent uncertainty, enough to separate a 16 % difference at 2 s if the drop lasts that long (10 m takes about 1.5 s for this package; use a lighter body or a higher drop). A single fall time from 2 m is marginal: the two models differ by roughly 15 ms there, comparable to a 10 ms timing resolution, so several repeats would be needed.

**Check.** Both solutions give $v \to v_t$ as $t \to \infty$ and $v \approx gt$ as $t \to 0$ (expand: $v_t(1 - e^{-t/\tau}) \approx v_t t/\tau = gt$). ✓ Dimensions of $c$: kg/s, so $cv$ is $\mathrm{kg\,m\,s^{-2}} = \mathrm{N}$. ✓

**Common mistakes.** Matching $c$ to $k$ by setting $c = k$ (different units); concluding that linear drag is "weaker" because it has a lower power of $v$ without fixing what is held equal.

## Check-yourself answers

1. $m\,\mathrm{d}v/\mathrm{d}t = mg - kv^2$ (downward positive, $v \ge 0$) with $k = \tfrac12\rho C_D A$. $k$ is constant if the air density is constant (altitude change small compared with 8 km) and $C_D A$ is constant (no tumbling, low Mach, Reynolds number well above 1).
2. $[k] = \mathrm{kg\,m^{-3}\cdot m^2 = kg/m}$; $[v_t] = \sqrt{\mathrm{kg\cdot m\,s^{-2}/(kg\,m^{-1})}} = \mathrm{m/s}$.
3. `a == b` tests whether the two expression *trees* are identical, so $(x+1)^2$ and $x^2 + 2x + 1$ compare unequal although they are the same function. `sp.simplify(a - b) == 0` asks whether the difference reduces to the zero expression, which is a mathematical test (though `simplify` can fail to find a reduction; then `expand`, `trigsimp` or numerical spot checks help).
4. $k \to 0$: $v \to gt$, found with `sp.limit(v, k, 0)`. $t \to \infty$: $v \to v_t$, found with `sp.limit(v, t, sp.oo)`.
5. The *formula* is exact for a quadratic, but each of $f(x+h)$, $f(x-h)$, their difference and the division is rounded to 53 bits; the result differs from the symbolic $-2kv/m$ by a few ulps, growing as $h$ shrinks because the subtraction cancels leading digits.
6. $t_{99} = 2.647\,v_t/g = 2.647 \times 4.470 = 11.83$ s; $x_{99} = (v_t^2/g)\ln\cosh(2.647) = 195.9 \times 1.959 = 383.8$ m.
7. With all symbols positive, `sp.solve` returns only the root $\sqrt{gm/k}$ instead of $\pm\sqrt{gm/k}$, and `sp.simplify` can reduce $\sqrt{gm/k}^2$ to $gm/k$ and $\sqrt{g}\sqrt{m}/\sqrt{k}$ to a single radical without worrying about branch cuts.
8. Slope $+2$ on the large-$h$ side (truncation, $h^2 f'''/6$) and slope $-1$ on the small-$h$ side (rounding, $\varepsilon|f|/h$); the right-hand side is set by the machine epsilon $\varepsilon = 2^{-52} \approx 2.2 \times 10^{-16}$ of `float64`.
