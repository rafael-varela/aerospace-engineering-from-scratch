# P01.08  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results to four significant figures. Angles in radians inside every formula; $\rho = 1.225$ kg/m³.

## 1. A linear wing

**Assumptions.** Linear $C_L$ ($k_3 = 0$), incompressible, fixed density.

**Derivation** ($s = 15$: $S = 13.0$ m², $a = 4.8$ /rad). (a) $d = 7° = 0.12217$ rad; $C_L = 4.8 \times 0.12217 = 0.5864$; $q = \tfrac12 \times 1.225 \times 2500 = 1531.25$ Pa; $L = 1531.25 \times 13.0 \times 0.5864 = 11\,674$ N.
(b) $L_\alpha = qSa = 1531.25 \times 13 \times 4.8 = 95\,550$ N/rad $= 1667.7$ N/deg; $L_V = 2L/V = 466.9$ N/(m/s).
(c) $+1°$: $\mathrm{d}L = 1667.7$ N; exact $L(6°) - L(5°) = qSa \times 1° = 1667.7$ N (linear in $\alpha$, so exact). $+1$ m/s: $\mathrm{d}L = 466.9$ N; exact $L[(51/50)^2 - 1] = 11\,674 \times 0.0404 = 471.6$ N. Both: $\mathrm{d}L = 2134.6$ N; exact $11\,674 \times [(51/50)^2 \times 8/7 - 1] = 2206.7$ N.
(d) Errors: $0$ %, $-0.99$ %, $-3.26$ % (the differential underestimates because the missing terms $\tfrac12 L_{VV}\Delta V^2 = 4.67$ N and $L_{\alpha V}\Delta\alpha\Delta V = 2L_\alpha\Delta\alpha\Delta V/V = 66.7$ N are both positive).

**Result.** $C_L = 0.5864$, $L = 11\,674$ N, $L_\alpha = 95\,550$ N/rad $= 1667.7$ N/deg, $L_V = 466.9$ N/(m/s); increments 1667.7 / 471.6 / 2206.7 N exact against 1667.7 / 466.9 / 2134.6 N linear; errors 0 / $-1.0$ / $-3.3$ %.

**Check.** $L_V\Delta V/L = 2\Delta V/V = 4$ % and the exact relative change is $(1.02)^2 - 1 = 4.04$ %: the 0.04 % gap is $\tfrac12(2\Delta V/V)^2$ ✓. Units: N/rad × rad = N ✓.

**Common mistakes.** Using degrees in $C_L = a\,d$ (lift 57 times too big); forgetting that the combined increment includes the cross term, and reporting $-1$ % for case (c) instead of $-3.3$ %.

## 2. The calculus of the lift function

**Derivation.** (a) With $d = \alpha - \alpha_0$, $q = \tfrac12\rho V^2$: $L_\alpha = qS(a - 3k_3d^2)$ ($V$ fixed); $L_V = \rho VS(ad - k_3d^3) = 2L/V$ ($\alpha$ fixed); $L_{\alpha\alpha} = qS(-6k_3d)$; $L_{VV} = \rho S(ad - k_3d^3)$; $L_{\alpha V} = \partial_V[qS(a - 3k_3d^2)] = \rho VS(a - 3k_3d^2)$; $L_{V\alpha} = \partial_\alpha[\rho VS(ad - k_3d^3)] = \rho VS(a - 3k_3d^2)$. Equal: both are the derivative of the product $\tfrac12\rho V^2 \cdot SC_L(\alpha)$ once in each factor, and multiplication commutes.
(b) Implicit: $L_\alpha + L_VV' = 0 \Rightarrow V' = -qSC_L'/(\rho VSC_L) = -\tfrac12 V C_L'/C_L$. Explicit: $V = (2W/\rho S)^{1/2}C_L^{-1/2}$, so $V' = (2W/\rho S)^{1/2}(-\tfrac12)C_L^{-3/2}C_L' = -\tfrac12 V C_L'/C_L$. Same.
(c) $L(\alpha, V + \Delta V) = \tfrac12\rho(V^2 + 2V\Delta V + \Delta V^2)SC_L = L + \rho VSC_L\Delta V + \tfrac12\rho SC_L\Delta V^2$: the first-order term is $L_V\Delta V$ and the remainder is exactly $\tfrac12\rho SC_L\Delta V^2 = \tfrac12 L_{VV}\Delta V^2$. In $\alpha$: $C_L(d + \Delta) = a(d + \Delta) - k_3(d + \Delta)^3 = C_L + (a - 3k_3d^2)\Delta - 3k_3d\Delta^2 - k_3\Delta^3$; multiply by $qS$: $L_\alpha\Delta + \tfrac12(-6k_3dqS)\Delta^2 + \tfrac16(-6k_3qS)\Delta^3$, i.e. $\tfrac12 L_{\alpha\alpha}\Delta^2 + \tfrac16 L_{\alpha\alpha\alpha}\Delta^3$ with $L_{\alpha\alpha\alpha} = -6k_3qS$, and nothing beyond.
(d) $L_{\alpha V} = \rho VSC_L' = (2/V)(qSC_L') = 2L_\alpha/V$, so $L_{\alpha V}\Delta\alpha\Delta V = (L_\alpha\Delta\alpha)(2\Delta V/V)$: the extra lift from the angle change is itself scaled by the relative change in dynamic pressure, $2\Delta V/V$. A $10$ m/s gust at 55 m/s makes every degree of $\alpha$ worth 36 % more.

**Check.** (a) at $k_3 = 0$: $L_{\alpha\alpha} = 0$, $L_{\alpha V} = \rho VSa$ ✓. (c) matches the numbers of lesson §2.4 (851 N and 236 N) ✓.

**Common mistakes.** Differentiating $V^2$ as $V^2$ in $L_{\alpha V}$; claiming the Taylor series in $\alpha$ has an $O(\Delta^4)$ term (it is a cubic; there is none).

## 3. Airspeed from a Pitot tube, with uncertainty

**Assumptions.** Incompressible Bernoulli ($V < 100$ m/s); independent Gaussian uncertainties.

**Derivation.** (a) $V = \sqrt{2 \times 1500/1.225} = 49.49$ m/s. $V_q = V/(2q) = 0.016496$ (m/s)/Pa; $V_\rho = -V/(2\rho) = -20.20$ (m/s)/(kg/m³).
(b) $V_q\delta q = 0.2474$ m/s; $V_\rho\delta\rho = -0.4040$ m/s; $\delta V = \sqrt{0.2474^2 + 0.4040^2} = 0.4737$ m/s. Density dominates (73 % of the variance).
(c) $\ln V = \text{const} + \tfrac12\ln q - \tfrac12\ln\rho$, so the logarithmic sensitivities are $+\tfrac12$ and $-\tfrac12$. Relative uncertainties: $\tfrac12 \times 1.0\,\% = 0.50\,\%$ from $q$, $\tfrac12 \times 1.63\,\% = 0.82\,\%$ from $\rho$, combined $0.957\,\%$ of 49.49 m/s $= 0.474$ m/s ✓.
(d) Indicated airspeed assumes $\rho_0$: $V_{\text{IAS}} = \sqrt{2q/\rho_0}$. True airspeed $V = \sqrt{2q/\rho} = V_{\text{IAS}}\sqrt{\rho_0/\rho} = V_{\text{IAS}}/\sqrt{0.85} = 1.085\,V_{\text{IAS}}$: 8.5 % higher. Lift $= qSC_L$ depends only on the measured $q$, so at a given *indicated* airspeed and angle of attack the lift is unchanged: that is why pilots fly approach speeds by indicated airspeed regardless of altitude or temperature. (Ground speed, and therefore landing distance, is 8.5 % higher.)

**Result.** $V = 49.49 \pm 0.47$ m/s; $V_q = 0.0165$ (m/s)/Pa, $V_\rho = -20.2$ (m/s)/(kg/m³); density dominates; TAS 8.5 % above IAS; lift at fixed IAS unchanged.

**Check.** Units: $V/q$ is (m/s)/Pa ✓. Doubling $q$ raises $V$ by $\sqrt2$: the sensitivity $\tfrac12$ ✓.

**Common mistakes.** Adding the two uncertainties linearly (0.65 m/s instead of 0.47); applying $V_\rho$ with the wrong sign; concluding that lower density means less lift at the same *indicated* speed.

## 4. The step-size study

**Result** (reference point; errors in N/rad and N/(m/s)):

| $h$ | $\alpha$ error | $V$ error |
|---|---|---|
| $10^{-1}$ | $3.56 \times 10^3$ | $1.1 \times 10^{-11}$ |
| $10^{-2}$ | $3.56 \times 10^1$ | $2.2 \times 10^{-10}$ |
| $10^{-3}$ | $3.56 \times 10^{-1}$ | $1.2 \times 10^{-9}$ |
| $10^{-4}$ | $3.56 \times 10^{-3}$ | $1.7 \times 10^{-8}$ |
| $10^{-5}$ | $3.57 \times 10^{-5}$ | $1.9 \times 10^{-7}$ |
| $10^{-6}$ | $5.5 \times 10^{-7}$ | $6.3 \times 10^{-7}$ |
| $10^{-7}$ | $6.7 \times 10^{-6}$ | $2.8 \times 10^{-7}$ |
| $10^{-8}$ | $1.3 \times 10^{-4}$ | $1.5 \times 10^{-4}$ |

(b) Left-hand slope ($h \ge 10^{-4}$): $2.000$; the prediction $k_3qSh^2 = 12 \times 1852.8 \times 16 \times h^2 = 3.557 \times 10^5h^2$ gives $3.557 \times 10^3$ at $h = 0.1$, matching to four figures (the cubic has no higher derivatives, so the "error term" is exact). Right-hand slope ($h \le 10^{-7}$): about $-1.3$ (rounding is noisy; $-1$ is the expectation). Minimum near $h = 10^{-6}$, where $(3\varepsilon|L|/|L_{\alpha\alpha\alpha}|)^{1/3} = (3 \times 2.2 \times 10^{-16} \times 15\,114/2.13 \times 10^6)^{1/3} = 1.7 \times 10^{-6}$.
(c) $L$ is quadratic in $V$, so $L_{VVV} = 0$ and the central difference has no truncation error at all; the column is pure rounding, $\sim\varepsilon|L|/h$ with $\varepsilon|L| \approx 3 \times 10^{-12}$: $10^{-11}$ at $h = 0.1$, $10^{-7}$ at $10^{-5}$, visible as a monotone rise from the first row.
(d) The two orders agree to $10^{-13}$ relative at every $h$ down to $10^{-5}$ (both are the same arithmetic with the roles of the two steps swapped), so "agreement to $10^{-6}$" holds for all tested $h$ including $10^{-1}$; the *accuracy* against the analytic $4964.42$ is another matter: $4835$ at $h = 0.1$, $4963.1$ at $10^{-2}$, $4964.41$ at $10^{-3}$.

**Check.** The $\alpha$ error at $h = 10^{-2}$ (35.6) equals the lesson's §2.7 prediction ✓.

**Common mistakes.** Mixing up "the two orders agree" (a symmetry of the algorithm) with "the result is accurate" (a question of $h$); fitting a single line through both branches.

## 5. Where the linear sensitivity stops

Reference $L_0 = 15\,114$ N; threshold $2\,\% = 302$ N.

(a) $2°$ gust: linear 4766 N, exact 4614 N, error 152 N $= 1.0\,\%$ of $L_0$. **Adequate.** (Remainder: $\tfrac12 L_{\alpha\alpha}\Delta\alpha^2 = -136$ N plus cubic $-15$ N.)
(b) $6°$ and $12$ m/s: linear 20 892 N, exact 26 105 N, error $-5213$ N $= -34\,\%$ of $L_0$. **Not adequate**; the cross term alone is $4964 \times 0.1047 \times 12 = 6239$ N. Use the second-order estimate (8), which recovers all but the cubic terms, or evaluate $L$ directly.
(c) $4°$ gust at $\alpha = 14°$: there $C_L' = 5 - 36 \times 0.2793^2 = 2.19$ /rad, less than half the value at 4°, $L_\alpha = 64\,900$ N/rad; linear 4538 N, exact 2964 N, error 1574 N $= 4.7\,\%$ of $L_0$ (and 53 % of the increment). **Not adequate**, and the model itself is unreliable this close to its stall; use the full curve, preferably a measured one.
(d) $2°$ gust at 95 m/s: by scaling, the error is again $1.0\,\%$ of the local $L_0$ (both the linear estimate and the remainder scale with $V^2$). **The linearization is adequate**, but $M = 95/340.3 = 0.28$ and the Prandtl–Glauert factor $1/\sqrt{1 - M^2} = 1.041$ says the incompressible $C_L$ is already 4 % low, twice the linearization error. The number that decides is the Mach number, and the thing to fix is the model, not the differential.

**Common mistakes.** Evaluating the partials at 4° for case (c); treating "adequate" as a property of the gust size alone rather than of gust size *and* operating point.

## 6. Lift through altitude

**Assumptions.** ISA troposphere; $\alpha$ and $V$ as stated; $L_\rho = L/\rho$ exactly since $L$ is linear in $\rho$.

**Derivation.** (a) $\partial L/\partial\rho = \tfrac12 V^2SC_L = L/\rho$. $\rho'(h) = -\rho_0 \times 4.2559 \times 2.2558 \times 10^{-5}(1 - 2.2558 \times 10^{-5}h)^{3.2559}$. At 2000 m: $\rho = 1.0065$ kg/m³, $\rho' = -1.0119 \times 10^{-4}$ kg/m⁴ (central difference: $-1.0119 \times 10^{-4}$ ✓); $L = 12\,418$ N, $L_\rho = 12\,338$ N/(kg/m³).
(b) $\mathrm{d}L/\mathrm{d}t = L_\rho\rho'\dot h = 12\,338 \times (-1.0119 \times 10^{-4}) \times 5 = -6.242$ N/s; central difference in $t$ with $\Delta t = 10^{-3}$ s: $-6.242$ N/s ✓. The aircraft loses 0.05 % of its lift per second of climb.
(c) $L_\alpha$ at 2000 m $= 136\,522 \times 1.0065/1.225 = 112\,170$ N/rad; $\dot\alpha = 6.242/112\,170 = 5.565 \times 10^{-5}$ rad/s $= 0.00319°$/s: the pilot adds a degree every five minutes. Stall altitude: the lift to hold is the sea-level 15 114 N; $C_{L,\max} = 1.242$ needs $\rho = 15\,114/(\tfrac12 \times 1.242 \times 16 \times 55^2) = 0.5027$ kg/m³, i.e. $h = [1 - (0.5027/1.225)^{1/4.2559}]/2.2558 \times 10^{-5} = 8371$ m. (A real light aircraft has no engine power to get there at 55 m/s, which is a different limit.)
(d) $L_V\dot V = 549.6 \times (1.0065/1.225) \times (-0.1) = -45.2$ N/s against the density term $-6.2$ N/s: the speed term dominates by a factor of 7.2, and the required $\dot\alpha$ becomes $(6.24 + 45.2)/112\,170 = 4.58 \times 10^{-4}$ rad/s $= 0.0263°$/s. The chain rule adds the contributions of every route by which time reaches $L$ (through $\rho$, through $V$, through $\alpha$), each weighted by its own partial; a single-variable derivative could follow only one of them.

**Result.** $L_\rho = 12\,338$ N/(kg/m³), $\rho'(2000) = -1.012 \times 10^{-4}$ kg/m⁴; $\mathrm{d}L/\mathrm{d}t = -6.24$ N/s (FD agrees to $10^{-10}$); $\dot\alpha = 0.0032°$/s; stall altitude 8.4 km; with $\dot V = -0.1$ m/s² the speed term dominates 7:1 and $\dot\alpha = 0.026°$/s.

**Check.** Units: (N/(kg/m³))(kg/m⁴)(m/s) = N/s ✓. Sea-level limit: $\rho' \to -1.18 \times 10^{-4}$ kg/m⁴, 17 % steeper, as the exponent's base is larger ✓.

**Common mistakes.** Using $L_\alpha$ from sea level at 2000 m (18 % high); forgetting that holding lift means the *sum* of the chain-rule terms is zero, not each term; quoting the stall altitude as a performance ceiling.

## Check-yourself answers

1. $\partial L/\partial V = \rho VSC_L = 2L/V$ ($\alpha$ fixed); $C_L$ is a constant factor when $\alpha$ is held fixed, so only $V^2$ is differentiated.
2. $V$ is held fixed. Along the trim curve $V$ changes with $\alpha$, and $\mathrm{d}L/\mathrm{d}\alpha = L_\alpha + L_V\,\mathrm{d}V/\mathrm{d}\alpha = 0$ there.
3. $L_{\alpha V} = L_{V\alpha}$ whenever both second partials are continuous near the point (true for any polynomial, as here).
4. The cross term $L_{\alpha V}\Delta\alpha\Delta V$: the extra angle acting at the changed speed.
5. $\mathrm{d}V/\mathrm{d}\alpha = -L_\alpha/L_V$; negative below stall (more angle, less speed), zero at $C_{L,\max}$.
6. 2 % (since $L \propto V^2$); 1 % in the linear model ($L \propto \alpha - \alpha_0$).
7. $L$ is quadratic in $V$, so $L_{VVV} = 0$ and the $h^2$ truncation term vanishes; it is cubic in $\alpha$, so $L_{\alpha\alpha\alpha} = -6k_3qS \ne 0$.
8. $\mathrm{d}L/\mathrm{d}t = L_\alpha\dot\alpha + L_V\dot V$ [N/s] $= 136\,522 \times 0.0349 - 549.6 \times 1 = 4766 - 550 = 4216$ N/s; the $\alpha$ term dominates by nearly nine to one.
