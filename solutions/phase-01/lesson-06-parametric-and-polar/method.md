# P01.06  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. $r_p = R_E + h_p$, $r_a = R_E + h_a$ (convert km to m). 2. $e = (r_a - r_p)/(r_a + r_p)$, $p = r_p(1 + e)$, $a = (r_p + r_a)/2$ and check $a = p/(1 - e^2)$. 3. $r(90°) = p$. 4. $T = 2\pi\sqrt{a^3/\mu}$, divide by 60. 5. $h = \sqrt{\mu p}$; $v_p = h/r_p$, $v_a = h/r_a$. 6. Differentiate $\tan\gamma$: the derivative vanishes where $\cos\theta = -e$; substitute back to get $\sin\gamma_{\max} = e$. 7. Compare each with the code; for (e) evaluate `flight_path_angle` on $10^5$ points and take the maximum.

## 2.
(a) Write $r = e(d + r\cos\theta)$, collect $r$, identify $p = ed$. (b) Partition $[\theta_0, \theta_1]$; bound each sector between circular sectors of radii $\min r$ and $\max r$ on the sub-interval; both Riemann sums converge to the same integral as $\Delta\theta \to 0$ because $r$ is continuous. Then $\mathrm{d}A/\mathrm{d}t = \tfrac12 r^2\dot\theta$ by the chain rule; constant $r^2\dot\theta$ gives constant $\mathrm{d}A/\mathrm{d}t$. (c) $v^2 = \dot r^2 + (r\dot\theta)^2$; expand $(1 + e\cos\theta)^2$; combine $e^2\sin^2\theta + e^2\cos^2\theta$. For $rv = h$ note $r v = h\sqrt{1 + 2e\cos\theta + e^2}/(1 + e\cos\theta)$ and set the ratio to 1: this forces $\sin\theta = 0$. (d) $x^2 + y^2 = a^2[(\cos E - e)^2 + (1 - e^2)\sin^2 E] = a^2(1 - e\cos E)^2$, so $r = a(1 - e\cos E)$; then $r + e x = a(1 - e^2) = p$, which is $r(1 + e\cos\theta) = p$.

## 3.
(a) $\dot x = V$, $\dot y = A\omega\cos\omega t$; the slope is largest where $\cos\omega t = \pm 1$; speed is $\sqrt{\dot x^2 + \dot y^2}$, largest at the crossings and smallest ($= V$) at the crests. (b) $\ddot x = 0$, $\ddot y = -A\omega^2\sin\omega t$; $\kappa = V|\ddot y|/(V^2 + \dot y^2)^{3/2}$, maximal at the crests where $\dot y = 0$: $\kappa_{\max} = A\omega^2/V^2$. (c) $a_{\text{lat}} = A\omega^2$ at the crest; $n = \sqrt{1 + (a_{\text{lat}}/g)^2}$. (d) `simpson(lambda t: np.sqrt(V**2 + (A*w*np.cos(w*t))**2), 0, T, 200)`; compare with $VT$.

## 4.
(a) $p = r_p(1 + e)$, $a = p/(1 - e^2) < 0$, $\theta_\infty = \arccos(-1/e)$, $\delta = 2\arcsin(1/e)$ (the asymptotes make angle $\theta_\infty$ with the apse line; the velocity turns through $2\theta_\infty - \pi$). (b) $v_\infty = \sqrt{\mu/|a|}$; $v_p$ from `orbital_speed(0, p, e)`; check energy. (c) Loop over $n$; reference with $n = 4096$; order $= \log_2(\text{err}_n/\text{err}_{2n})$. (d) The hyperbola's centre is at $x_c = -ae = |a|e$, which is positive: it lies beyond periapsis on the $+x$ side, a distance $|a|$ past the periapsis point. The asymptotes pass through the centre at angles $\pm\theta_\infty$ to the $+x$ axis.

## 5.
(a) Drag: density ratio between 150 and 300 km, roughly 100; the sign is a periapsis-pass decay of the apoapsis (period decreasing). (b) $J_2$: nodal precession rate $\dot\Omega = -\tfrac32 n J_2 (R_E/p)^2\cos i$; at $i = 90°$ it is zero, but the argument of periapsis still rotates; sign: the ascending-node longitude, which the conic says is fixed, drifts. (c) Third body: compare lunar acceleration difference $2\mu_{\text{Moon}} r/d^3$ to Earth's $\mu/r^2$ at $r = 200\,000$ km. (d) Compare $\Delta v = 2$ m/s with the orbital speed of about 7.67 km/s and estimate $\Delta a$ from $v\,\mathrm{d}v = (\mu/2a^2)\,\mathrm{d}a$; the orbits differ by a few kilometres and the relative drift is a few kilometres per orbit.

## 6.
(a) At $\theta = 90°$, $\cos E = e$; time from periapsis to there is $(T/2\pi)(\arccos e - e\sqrt{1 - e^2})$; double it for the symmetric half; divide by $T$. (b) $g(E) = E - e\sin E - \pi/2$, $g'(E) = 1 - e\cos E$; iterate; convert with the half-angle formula; check with `time_of_flight(p, e, 0, theta)/orbital_period(p, e)`. (c) Kepler's equation is a closed-form relation between $E$ and $t$ (one transcendental equation, solved by Newton in a few iterations); the area integral has no elementary antiderivative in $\theta$ for $e \ne 0$, so it needs quadrature.
