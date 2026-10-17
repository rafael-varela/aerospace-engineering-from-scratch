# P01.06  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for $s = 15$; numerical results are given to four significant figures, with the Simpson panel count where one was used. $\mu = 3.986\,004\,418 \times 10^{14}$ m³/s², $R_E = 6378.137$ km.

## 1. A low Earth orbit from its apsides

**Assumptions.** Two-body conic, periapsis at $\theta = 0$, no perturbations.

**Derivation** ($s = 15$: altitudes 350 km and 2500 km). $r_p = 6728.137$ km, $r_a = 8878.137$ km.
(a) $e = (8878.137 - 6728.137)/(8878.137 + 6728.137) = 2150/15\,606.27 = 0.13777$; $p = 6728.137 \times 1.13777 = 7655.04$ km; $a = (r_p + r_a)/2 = 7803.14$ km (check: $p/(1 - e^2) = 7655.04/0.98102 = 7803.1$ km ✓).
(b) $r(90°) = p = 7655.04$ km.
(c) $T = 2\pi\sqrt{(7.80314 \times 10^6)^3/3.986 \times 10^{14}} = 6859.85$ s $= 114.33$ min.
(d) $h = \sqrt{\mu p} = \sqrt{3.986 \times 10^{14} \times 7.65504 \times 10^6} = 5.5239 \times 10^{10}$ m²/s; $v_p = h/r_p = 8210.1$ m/s; $v_a = h/r_a = 6221.9$ m/s.
(e) $\tan\gamma = e\sin\theta/(1 + e\cos\theta)$. Quotient rule: $\mathrm{d}(\tan\gamma)/\mathrm{d}\theta = e[\cos\theta(1 + e\cos\theta) + e\sin^2\theta]/(1 + e\cos\theta)^2 = e(\cos\theta + e)/(1 + e\cos\theta)^2$, zero when $\cos\theta = -e$. There $\tan\gamma = e\sqrt{1 - e^2}/(1 - e^2) = e/\sqrt{1 - e^2}$, i.e. $\sin\gamma_{\max} = e$. So $\gamma_{\max} = \arcsin 0.13777 = 7.919°$ at $\theta = \arccos(-0.13777) = 97.92°$.

**Result.** $e = 0.1378$, $p = 7655$ km, $a = 7803$ km, $r(90°) = 7655$ km, $T = 114.3$ min, $v_p = 8210$ m/s, $v_a = 6222$ m/s, $\gamma_{\max} = 7.919°$ at $\theta = 97.92°$. The code agrees to every figure shown; `flight_path_angle` on $10^5$ points peaks at $7.9185°$ at $\theta = 97.918°$.

**Check.** $r_p v_p = r_a v_a$: $6728.137 \times 8210.1 = 8878.137 \times 6221.9 = 5.524 \times 10^{10}$ m²/s $= h$ ✓. Units of $h$: m · m/s = m²/s ✓. For $e \to 0$, $\gamma_{\max} \to 0$ at $\theta \to 90°$, consistent with a circle.

**Common mistakes.** Using altitudes instead of radii in (3); quoting $T$ from $a$ in km with $\mu$ in m³/s²; taking $\gamma_{\max}$ at $\theta = 90°$ (it is past 90°, by $\arcsin e$).

## 2. The conic and its calculus

**Derivation.** (a) Point $P$ at $(r, \theta)$: distance to the focus $= r$; distance to the line $x = -d$ is $x + d = r\cos\theta + d$. The definition gives $r = e(d + r\cos\theta)$, so $r(1 - e\cos\theta) = ed$. (With the directrix on the *left* and $\theta$ measured so that periapsis is at $\theta = 0$ on the right, replace $\theta$ by $\pi - \theta$: $r(1 + e\cos\theta) = ed$.) Hence $r = p/(1 + e\cos\theta)$ with $p = ed$, and $r(\pm\pi/2) = p$ identifies $p$ as the semi-latus rectum.

(b) Partition $[\theta_0, \theta_1]$ into $N$ sub-intervals of width $\Delta\theta$. On the $k$-th, $r$ lies between $m_k = \min r$ and $M_k = \max r$ (it exists: $r$ is continuous). The swept region on that sub-interval contains the circular sector of radius $m_k$ and is contained in the one of radius $M_k$, whose areas are $\tfrac12 m_k^2\Delta\theta$ and $\tfrac12 M_k^2\Delta\theta$ (a sector of angle $\Delta\theta$ is the fraction $\Delta\theta/2\pi$ of a disc). Summing, $\tfrac12\sum m_k^2\Delta\theta \le A \le \tfrac12\sum M_k^2\Delta\theta$; both are Riemann sums of $\tfrac12 r^2$ and converge to $\tfrac12\int r^2\,\mathrm{d}\theta$ as $\Delta\theta \to 0$ because $r^2$ is continuous, so $A$ is squeezed to that value. "Almost circular" is used only in bounding each piece between two exact sectors; no approximation survives the limit. Now let $\theta = \theta(t)$ with $A(t) = \tfrac12\int_{\theta_0}^{\theta(t)} r^2\,\mathrm{d}\theta$; the fundamental theorem and the chain rule give $\dot A = \tfrac12 r(\theta)^2\dot\theta$. If $r^2\dot\theta = h$ is constant then $\dot A = h/2$ is constant: equal areas in equal times.

(c) $v^2 = \dot r^2 + (r\dot\theta)^2 = (h/p)^2[e^2\sin^2\theta + 1 + 2e\cos\theta + e^2\cos^2\theta] = (h/p)^2[1 + 2e\cos\theta + e^2]$. Then $rv = h\sqrt{1 + 2e\cos\theta + e^2}/(1 + e\cos\theta)$; $rv = h$ requires $1 + 2e\cos\theta + e^2 = 1 + 2e\cos\theta + e^2\cos^2\theta$, i.e. $e^2\sin^2\theta = 0$, so $\sin\theta = 0$: the apsides only (for $e > 0$).

(d) $x^2 + y^2 = a^2[\cos^2 E - 2e\cos E + e^2 + (1 - e^2)\sin^2 E] = a^2[1 - 2e\cos E + e^2\cos^2 E] = a^2(1 - e\cos E)^2$, so $r = a(1 - e\cos E)$ (positive since $e < 1$). Then $r + ex = a(1 - e\cos E) + ea(\cos E - e) = a(1 - e^2) = p$. Since $x = r\cos\theta$, this reads $r(1 + e\cos\theta) = p$. ✓

**Check.** (c) at $\theta = 0$: $v = (h/p)(1 + e) = h/r_p$ ✓. (d) at $E = 0$: $x = a(1 - e) = r_p$, $y = 0$ ✓.

**Common mistakes.** In (a), forgetting that the directrix side fixes the sign of $\cos\theta$; in (b), claiming the sector *is* circular rather than bounding it; in (d), expanding $(\cos E - e)^2$ and losing the $e^2$ term that completes the square.

## 3. A weaving ground track

**Assumptions.** Flat ground, constant $V$ along $x$, no wind; $g = 9.80665$ m/s². $\omega = 2\pi/20 = 0.31416$ rad/s, $A\omega = 15.708$ m/s, $A\omega^2 = 4.9348$ m/s².

**Derivation.** (a) $\mathrm{d}y/\mathrm{d}x = \dot y/\dot x = (A\omega/V)\cos\omega t$, largest in magnitude at the crossings ($\cos\omega t = \pm 1$): $\arctan(15.708/60) = 14.67°$. Speed $\sqrt{V^2 + A^2\omega^2\cos^2\omega t}$: maximum $\sqrt{60^2 + 15.708^2} = 62.02$ m/s at the crossings, minimum $V = 60.00$ m/s at the crests.
(b) $\ddot x = 0$, $\ddot y = -A\omega^2\sin\omega t$, so $\kappa = V A\omega^2|\sin\omega t|/(V^2 + A^2\omega^2\cos^2\omega t)^{3/2}$, maximal at the crests ($\sin\omega t = \pm 1$, $\cos\omega t = 0$): $\kappa_{\max} = A\omega^2/V^2 = 4.9348/3600 = 1.3708 \times 10^{-3}$ m⁻¹, minimum radius $729.5$ m.
(c) At a crest the velocity is $(V, 0)$ and the acceleration is $(0, \mp A\omega^2)$, purely lateral: $a_{\text{lat}} = 4.935$ m/s² $= 0.5032\,g$ (also $= V^2\kappa_{\max}$ ✓). Load factor $n = \sqrt{1 + 0.5032^2} = 1.119$.
(d) `simpson(lambda t: np.sqrt(60**2 + (15.708*np.cos(0.31416*t))**2), 0, 20, 200)` $= 1220.3$ m against $VT = 1200.0$ m: the weave adds 1.7 %.

**Result.** 14.67° (exact), 62.02 and 60.00 m/s (exact), $R_{\min} = 729.5$ m (exact), $a_{\text{lat}} = 4.935$ m/s², $n = 1.119$ (exact), path length 1220.3 m (Simpson, $n = 200$; changing to $n = 400$ alters it by less than $10^{-6}$ m).

**Check.** Units: $A\omega^2$ is m/s² ✓; $\kappa$ in 1/m ✓. For $A \to 0$ the path length $\to VT$ and $n \to 1$. ✓

**Common mistakes.** Using $\mathrm{d}y/\mathrm{d}t$ as the slope (it is $\mathrm{d}y/\mathrm{d}x$); computing curvature from $y(x)$ with the formula $|y''|/(1 + y'^2)^{3/2}$ but forgetting to convert $\omega t$ to $x$ ($x = Vt$ makes $y(x) = A\sin(\omega x/V)$, which gives the same answer when done right).

## 4. A hyperbolic flyby

**Assumptions.** Two-body, Earth point mass, $r_p = 6678.137$ km, $e = 1.5$.

**Derivation.** (a) $p = r_p(1 + e) = 16\,695.3$ km; $a = p/(1 - e^2) = 16\,695.3/(-1.25) = -13\,356.3$ km; $\theta_\infty = \arccos(-1/1.5) = 131.81°$; $\delta = 2\arcsin(1/1.5) = 83.62°$ (equivalently $2\theta_\infty - 180°$).
(b) $v_\infty = \sqrt{\mu/|a|} = \sqrt{3.986 \times 10^{14}/1.33563 \times 10^7} = 5462.9$ m/s; `orbital_speed(0, p, 1.5)` $= 12\,215.5$ m/s. Energy: $v_p^2 = 1.49218 \times 10^8$ m²/s²; $v_\infty^2 + 2\mu/r_p = 2.9844 \times 10^7 + 1.19374 \times 10^8 = 1.49218 \times 10^8$ ✓ (agreement to $10^{-15}$ relative in the code).
(c) Time from $-60°$ to $+60°$: 1474.61 s (24.6 min). Errors against $n = 4096$: $4.91 \times 10^{-2}$, $3.16 \times 10^{-3}$, $1.99 \times 10^{-4}$, $1.25 \times 10^{-5}$ s for $n = 16, 32, 64, 128$; observed orders 3.96, 3.99, 4.00.
(d) Centre at $x_c = -ae = 20\,034.4$ km on the $+x$ side; asymptotes through it at $\pm 131.81°$ from $+x$. Periapsis at $(6678.1, 0)$ km.

**Result.** $p = 16\,695$ km, $a = -13\,356$ km, $\theta_\infty = 131.8°$, $\delta = 83.62°$, $v_\infty = 5463$ m/s, $v_p = 12\,216$ m/s, $t_{\pm 60°} = 1474.6$ s (Simpson, order 4 observed).

**Check.** At $\theta = 60°$, $r = p/(1 + 0.75) = 9540.2$ km and `orbital_speed` gives 10 649 m/s; vis-viva with negative $a$: $\sqrt{\mu(2/r + 1/|a|)} = 10\,649$ m/s ✓. $\gamma(60°) = \arctan(1.5\sin 60°/1.75) = 36.59°$.

**Common mistakes.** Using $|a|$ in vis-viva with a plus sign already written; drawing the asymptotes through the focus instead of the centre; evaluating `conic_radius` at $\theta = \pm\theta_\infty$ exactly (division by zero, caught by the `ValueError`).

## 5. Where the conic model stops

(a) **No-drag assumption.** The density ratio between 150 km and 300 km is about 100. For a ballistic coefficient $m/(C_D A) = 100$ kg/m² and an effective periapsis-pass length of order 1500 km, the per-pass speed loss $\Delta v \approx \rho v L/(2\beta)$ is about 0.14 m/s at 150 km (about 8 km off the apoapsis per pass, via $\Delta a = 2a^2 v\,\Delta v/\mu$) against 0.002 m/s at 300 km (0.1 km per pass). Measurable sign: the period, which the conic says is constant, shortens every revolution; the apoapsis radius decays while the periapsis barely moves.
(b) **Point-mass central body.** At 800 km and $i = 90°$ the nodal rate $\propto\cos i$ vanishes but the argument of periapsis rotates at $\tfrac34 n J_2(R_E/p)^2(5\cos^2 i - 1) = -3.29°$ per day, a full turn in 109 days; for a sun-synchronous $i = 98°$ the node precesses $+0.92°$/day, which is the point of that orbit. Measurable sign: the direction of periapsis (or of the node) drifts secularly; a conic has none.
(c) **No third body.** At $r = 200\,000$ km the Moon's tidal acceleration $2\mu_{\text{Moon}}r/d^3 \approx 3.5 \times 10^{-5}$ m/s² is 0.35 % of Earth's $\mu/r^2 = 9.97 \times 10^{-3}$ m/s² (the Sun adds half as much again). Over one 10-day orbit the apoapsis moves by thousands of kilometres; the elements change orbit to orbit. Measurable sign: $p$ and $e$ fitted on successive revolutions are not constant.
(d) **Same conic, for planning.** $\Delta v = 2$ m/s on $v = 7657$ m/s at 420 km changes $a$ by $\Delta a = 2a^2 v\,\Delta v/\mu = 3.55$ km and the period by $\tfrac32 T\,\Delta a/a = 4.4$ s. The two conics are indistinguishable on a plot, but the along-track separation grows by $v\,\Delta T \approx 33$ km per orbit (9 km if the separation is radial, since then $\Delta a$ is second order), so after a day they are several hundred kilometres apart. Different conic for conjunction screening, same conic for everything else. Measurable sign: the range between the two objects grows linearly with time.

**Common mistakes.** Quoting the density ratio as the effect (the effect is the integral over the pass); forgetting that $J_2$ changes the *orientation* but not the shape at first order; treating (d) as "negligible" without saying for what decision.

## 6. Time from the parametric form

**Derivation.** (a) With $\cos E = (e + \cos\theta)/(1 + e\cos\theta)$, $\theta = 90°$ gives $\cos E = e$, $E_1 = \arccos e$, $\sin E_1 = \sqrt{1 - e^2}$. Time from periapsis to $\theta = 90°$ is $(T/2\pi)(E_1 - e\sin E_1)$; the periapsis half is symmetric, so $t_p = (T/\pi)(\arccos e - e\sqrt{1 - e^2})$ and $f_p = (\arccos e - e\sqrt{1 - e^2})/\pi$. For $e = 0.72654$: $\arccos e = 0.75745$, $e\sqrt{1 - e^2} = 0.72654 \times 0.68713 = 0.49923$, $f_p = 0.25822/\pi = 0.08222$. `time_of_flight(p, e, -np.pi/2, np.pi/2)/orbital_period(p, e)` returns 0.082219 (agreement to $10^{-9}$ at $n = 400$). The transfer orbit spends 8.2 % of its period in the near half and 91.8 % in the far half.

(b) $M = \pi/2$. Newton on $g(E) = E - e\sin E - \pi/2$, $g' = 1 - e\cos E$, from $E_0 = 1.570796$:

| iterate | $E$ |
|---|---|
| 1 | 2.297338 |
| 2 | 2.173593 |
| 3 | 2.170544 |
| 4 | 2.1705417267 |
| 5 | 2.1705417267 |

The number of correct digits roughly doubles each step (quadratic convergence). $\tan(\theta/2) = \sqrt{1.72654/0.27346}\tan(1.08527) = 2.5127 \times 1.8957 = 4.763$, so $\theta = 2\arctan(4.763) = 156.28°$. Check: `time_of_flight(p, e, 0, np.radians(156.28))/T` $= 0.2500$ ✓.

(c) Kepler's equation is a closed-form relation between $E$ and $t$ whose only cost is solving one transcendental equation, which Newton does to machine precision in four or five steps with no discretisation error. The area integral $\int r^2\,\mathrm{d}\theta$ has no elementary antiderivative for $0 < e$ ($\int\mathrm{d}\theta/(1 + e\cos\theta)^2$ is doable but leads back to the same $\arctan$ that defines $E$), so `time_of_flight` discretises and carries a panel-dependent error. Inside a trajectory code called a million times, use Kepler's equation.

**Check.** At $e = 0$, $f_p = (\pi/2)/\pi = 0.5$ ✓. As $e \to 1$, $\arccos e \to 0$ and $f_p \to 0$: a near-parabola spends all its time far away ✓.

**Common mistakes.** Using $\cos E = \cos\theta$ (the two anomalies coincide only at the apsides); converting $E$ to $\theta$ with $\arctan$ without the half-angle formula (loses the quadrant); starting Newton from $E_0 = 0$ for $M$ near $\pi$ (converges, but slower).

## Check-yourself answers

1. $r = p$, because $\cos 90° = 0$ kills the $e$ term: all conics with the same $p$ pass through the two points at $\pm 90°$.
2. $e = (3 - 1)/(3 + 1) = 0.5$.
3. $A = \tfrac12\int r^2\,\mathrm{d}\theta$; differentiating in time gives $\dot A = \tfrac12 r^2\dot\theta$, so if $r^2\dot\theta$ is constant the area is swept at a constant rate.
4. Periapsis and apoapsis ($\sin\theta = 0$); at $\theta = 90°$, $\tan\gamma = e$.
5. Simpson fits one parabola through each *pair* of panels (three points); with an odd count the last panel has no partner.
6. $rv$ is constant at the apsides, so $v_a = 10\,151 \times 6678/42\,164 \approx 10\,151/6.31 \approx 1600$ m/s (exact: 1608 m/s).
7. Where $1 + 2\cos\theta = 0$: $\theta = \pm 120°$.
8. $\theta$ (the true anomaly) is the angle at the focus, i.e. at the Earth's centre; $E$ is measured at the ellipse's centre, $ae$ away.
