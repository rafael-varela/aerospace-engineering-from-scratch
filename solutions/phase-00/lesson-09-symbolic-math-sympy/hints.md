# P00.09  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. Everything follows from $k$ and $v_t$. The time to 99 % is $\tau\,\operatorname{artanh}(0.99)$, and the distance is $x(t) = (v_t^2/g)\ln\cosh(gt/v_t)$ evaluated at that time, so $\ln\cosh(\operatorname{artanh} 0.99)$ is a pure number you can compute once.
2. Factor $g - (k/m)v^2 = (g/v_t^2)(v_t^2 - v^2)$ and use partial fractions or the table integral $\int \mathrm{d}v/(a^2 - v^2) = \tfrac{1}{a}\operatorname{artanh}(v/a)$ for $|v| < a$. For (c), $\ln\cosh u = u - \ln 2 + \ln(1 + e^{-2u})$.
3. Upward, both gravity and drag point down, so the equation is $\mathrm{d}v/\mathrm{d}t = -g - (k/m)v^2$: the sign of the $v^2$ term does not flip with $v$. The integral $\int \mathrm{d}v/(a^2 + v^2)$ gives an arctangent instead of an inverse hyperbolic tangent.
4. The left slope is set by the truncation term $h^2 v'''/6$, the right by rounding $\varepsilon |v|/h$. Get $v'''$ from `sp.diff(velocity_expr(), t, 3)` and substitute numbers; do not estimate it by finite differences.
5. From the series, $1 - v/(gt) \approx gkt^2/(3m)$. For (c) the question is not whether drag is *present* but whether it moves the measured fall time by more than the timing resolution; compute the fall time with and without drag.
6. For each case ask: which of $g$, $\rho$, $C_D$, $A$ stops being constant, or which power of $v$ the drag follows? One number each: a density ratio, a $\sqrt{\cdot}$ spread in $v_t$, a Reynolds number, an altitude from $(R/(R+h))^2$.
7. Matching $v_t$ means $c = mg/v_t$. The linear model's solution is an exponential; compare $\tanh(u)$ and $1 - e^{-u}$ at the same $u$. Early on, which drag law bites harder at small $v$?
