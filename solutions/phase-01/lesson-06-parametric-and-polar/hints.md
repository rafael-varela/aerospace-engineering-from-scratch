# P01.06  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything follows from $r_p$ and $r_a$ through $e = (r_a - r_p)/(r_a + r_p)$ and $p = r_p(1 + e)$. For the largest flight-path angle, differentiate $\tan\gamma = e\sin\theta/(1 + e\cos\theta)$ with the quotient rule; the numerator of the derivative is $e(\cos\theta + e)$.
2. (a) The distance from $(r, \theta)$ to the line $x = -d$ is $d + r\cos\theta$. (b) A circular sector of angle $\Delta\theta$ has area $\tfrac12 r^2\Delta\theta$. (c) Square and add; $\sin^2 + \cos^2 = 1$. (d) Compute $x^2 + y^2$, use $b^2 = a^2(1 - e^2)$, and look for a perfect square.
3. $\mathrm{d}y/\mathrm{d}x = \dot y/\dot x = (A\omega/V)\cos\omega t$; curvature is largest where $\ddot y$ is largest and $\dot y = 0$, i.e. at the crests. Lateral acceleration there is just $|\ddot y|$.
4. $p = r_p(1 + e)$ still holds; $a = p/(1 - e^2)$ is negative. The asymptote angle is where $1 + e\cos\theta = 0$. For the convergence order, take $\log_2$ of successive error ratios against a run with $n = 4096$.
5. Ask for each case which of "no drag", "point-mass central body", "no third body" fails, and look for a quantity that should be constant under the conic model ($p$, $e$, the orbit plane, the period) but drifts.
6. At $\theta = 90°$, $\cos E = e$, so $E = \arccos e$ and $\sin E = \sqrt{1 - e^2}$; the periapsis half runs from $-E$ to $+E$. Newton's step for $g(E) = E - e\sin E - M$ is $E \leftarrow E - g/(1 - e\cos E)$.
