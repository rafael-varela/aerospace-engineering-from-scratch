# P01.06  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and say which results are exact and which are numerical (with the panel count used).

Several problems use your personal seed **s**: the day of the month you were born (1–31). Use $\mu = 3.986\,004\,418 \times 10^{14}$ m³/s² and $R_E = 6378.137$ km.

## 1. (Warm-up: a low Earth orbit from its apsides)

A satellite's periapsis altitude is $200 + 10s$ km and its apoapsis altitude is $1000 + 100s$ km.

Find: (a) $e$, $p$ and $a$; (b) the radius at $\theta = 90°$; (c) the period in minutes; (d) the speed at periapsis and at apoapsis; (e) the largest flight-path angle on the orbit and the true anomaly at which it occurs (derive the condition $\mathrm{d}\gamma/\mathrm{d}\theta = 0$ from equation (6) of the lesson before you compute). Check (a)–(d) with `conic_from_apsides`, `orbital_period` and `orbital_speed`, and (e) with a fine grid of `flight_path_angle`. Answer form: seven numbers with units, four significant figures, plus the two-line derivation for (e).

## 2. (Derivation: the conic and its calculus)

(a) Starting from the focus–directrix definition (distance to focus $= e \times$ distance to directrix, directrix at $x = -d$), derive $r = p/(1 + e\cos\theta)$ and identify $p$.
(b) Derive the polar area formula $A = \tfrac12\int r^2\,\mathrm{d}\theta$ from a sum of thin sectors, then show that constant $r^2\dot\theta$ implies equal areas in equal times. State where you use that the sector is "almost circular" and why the error vanishes in the limit.
(c) From $\dot r = (h/p)e\sin\theta$ and $r\dot\theta = (h/p)(1 + e\cos\theta)$, derive the speed $v = (h/p)\sqrt{1 + 2e\cos\theta + e^2}$ and show that $rv = h$ exactly at the apsides and nowhere else on an ellipse with $e > 0$.
(d) Show that every point of $x = a(\cos E - e)$, $y = a\sqrt{1 - e^2}\sin E$ satisfies $r = p/(1 + e\cos\theta)$ with $r = \sqrt{x^2 + y^2}$ and $\cos\theta = x/r$. (Hint: compute $r$ first and show $r = a(1 - e\cos E)$.)

Answer form: four derivations with each assumption named where it is used.

## 3. (Unfamiliar scenario: a weaving ground track)

A survey aircraft flies a sinusoidal track over flat ground: $x(t) = Vt$, $y(t) = A\sin(\omega t)$ with $V = 60$ m/s, $A = 50$ m and a weave period of 20 s ($\omega = 2\pi/20$ rad/s).

(a) Find the maximum heading deviation from the $x$ axis (the maximum of $\arctan(\mathrm{d}y/\mathrm{d}x)$) and the maximum and minimum ground speed.
(b) The curvature of a parametric curve is $\kappa = |\dot x\ddot y - \dot y\ddot x|/(\dot x^2 + \dot y^2)^{3/2}$. Find the minimum turning radius $1/\kappa_{\max}$ and where on the weave it occurs.
(c) Find the lateral acceleration at that point and the load factor $n = \sqrt{1 + (a_{\text{lat}}/g)^2}$ of a coordinated turn there.
(d) Use your `simpson` to compute the ground distance flown in one weave period and compare it with the straight-line distance $VT$.

Answer form: six numbers with units; say which are exact and which are numerical.

## 4. (Numerical, uses your code: a hyperbolic flyby)

A spacecraft approaches Earth on a hyperbola with $e = 1.5$ and periapsis altitude 300 km.

(a) Compute $p$, $a$ (negative), the asymptote angle $\theta_\infty$ and the total turning angle of the velocity, $\delta = 2\arcsin(1/e)$.
(b) Compute the hyperbolic excess speed $v_\infty = \sqrt{\mu/(-a)}$ and the periapsis speed, and verify $v_p^2 = v_\infty^2 + 2\mu/r_p$ (energy conservation).
(c) With `time_of_flight`, compute the time from $\theta = -60°$ to $+60°$ for $n = 16, 32, 64, 128$ and report the observed convergence order.
(d) Plot the hyperbola with `orbit_points`, add the two asymptotes through the centre at $x = -ae$ (note $ae > 0$ for $e > 1$ and $a < 0$), and mark periapsis. Save as `flyby.png`.

Answer form: seven numbers with units, a convergence table with the observed order, one plot.

## 5. (Judgment: where the conic model stops)

For each of the following say which assumption of the two-body conic fails, estimate the size of the effect with one number, and name one measurable sign in tracking data that it has:

(a) a transfer orbit whose periapsis is at 150 km altitude instead of 300 km (air density at 150 km is roughly $2 \times 10^{-9}$ kg/m³, at 300 km roughly $2 \times 10^{-11}$ kg/m³);
(b) a polar orbit at 800 km altitude flown for a year (Earth's $J_2 = 1.083 \times 10^{-3}$ makes the orbit plane precess at a rate proportional to $J_2 (R_E/p)^2\cos i$);
(c) an apoapsis at 200 000 km;
(d) a CubeSat deployed from the International Space Station with a 2 m/s separation speed: is its orbit "a different conic" or "the same conic" for mission-planning purposes?

Answer form: four short paragraphs, each with a number and a measurable sign.

## 6. (Hard: time from the parametric form, Kepler's equation)

The two parameters of §2.7 are linked by $\cos E = (e + \cos\theta)/(1 + e\cos\theta)$, and the time since periapsis is $t = (T/2\pi)(E - e\sin E)$ (Kepler's equation; take it as given or derive it from $\mathrm{d}A/\mathrm{d}t = h/2$ with the parametric ellipse).

(a) Show that the fraction of the period spent in the periapsis half of the orbit ($-90° < \theta < 90°$) is $f_p = (\arccos e - e\sqrt{1 - e^2})/\pi$, and evaluate it for the lesson's transfer orbit ($e = 0.72654$). Compare with `time_of_flight`.
(b) For the transfer orbit, find the true anomaly at $t = T/4$ by solving Kepler's equation for $E$ with Newton's method (`P01.02`) starting from $E_0 = M = \pi/2$; list the iterates to show quadratic convergence; convert to $\theta$ with $\tan(\theta/2) = \sqrt{(1 + e)/(1 - e)}\tan(E/2)$. Check by computing `time_of_flight(p, e, 0, theta)/T`.
(c) Explain in two sentences why `time_of_flight` needs numerical integration while Kepler's equation does not, and which one you would put inside a trajectory code that calls it a million times.

Answer form: one derivation, one fraction, a Newton table, one angle in degrees, two sentences.
