# P01.08  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. $d = \alpha - \alpha_0$ in radians; $C_L = a\,d$. 2. $q = \tfrac12\rho V^2$; $L = qSC_L$. 3. $L_\alpha = qSa$ (N/rad; multiply by $\pi/180$ for N/deg); $L_V = \rho V S C_L = 2L/V$. 4. Increments: $\mathrm{d}L = L_\alpha\Delta\alpha + L_V\Delta V$; exact $\Delta L = L(\alpha + \Delta\alpha, V + \Delta V) - L$. 5. Percentage error $= (\mathrm{d}L - \Delta L)/\Delta L$. 6. Build `Wing(S=..., a=..., k3=0.0)` and compare.

## 2.
(a) Differentiate (1) treating $V$ as constant for $L_\alpha$ and $\alpha$ as constant for $L_V$; repeat for the second partials; $L_{\alpha V}$ from $L_\alpha$ in $V$, $L_{V\alpha}$ from $L_V$ in $\alpha$. (b) Implicit: $\mathrm{d}/\mathrm{d}\alpha$ of $L(\alpha, V(\alpha)) = W$ gives $L_\alpha + L_V V' = 0$; substitute (4), (5). Explicit: $V = (2W/\rho S)^{1/2}C_L^{-1/2}$, differentiate the power. (c) $L(\alpha, V + \Delta V) = \tfrac12\rho(V + \Delta V)^2 SC_L$; expand the square: the terms are $L + L_V\Delta V + \tfrac12\rho SC_L\Delta V^2$ and nothing more. For $\alpha$: expand $C_L(d + \Delta\alpha)$ as a cubic in $\Delta\alpha$. (d) $L_{\alpha V} = \rho VSC_L' = (2/V)\cdot\tfrac12\rho V^2SC_L' = 2L_\alpha/V$; multiply by $\Delta\alpha\Delta V$.

## 3.
(a) $V = (2q/\rho)^{1/2}$; $V_q = V/(2q)$, $V_\rho = -V/(2\rho)$. (b) Multiply each partial by its uncertainty, square, add, root. (c) $\ln V = \tfrac12\ln 2 + \tfrac12\ln q - \tfrac12\ln\rho$; relative uncertainty $= \sqrt{(\tfrac12\delta q/q)^2 + (\tfrac12\delta\rho/\rho)^2}$. (d) Indicated airspeed is $\sqrt{2q/\rho_0}$; true is $\sqrt{2q/\rho}$; ratio $\sqrt{\rho_0/\rho}$. Lift $= qSC_L$ depends on $q$, which the instrument measures directly.

## 4.
(a) Loop over `h`; `abs(partial_derivative(f, x0, i, h) - analytic)`; `ax.loglog`. (b) `np.polyfit(np.log10(h), np.log10(err), 1)` on each branch; prediction $\tfrac16 h^2 \cdot 6k_3 qS = k_3 qS h^2$. (c) $L_{VVV} = 0$; the only error is rounding, which grows as $h$ shrinks and becomes visible when $\varepsilon|L|/h$ exceeds the $10^{-10}$ level, i.e. for $h$ below about $10^{-5}$. (d) Loop `mixed_partial(f, x0, 0, 1, h)` and `(1, 0, h)`; relative difference.

## 5.
For each: compute $L_0$, the partials and second partials at the stated point, the linear estimate, the remainder terms of (8) plus the cubic term, and the exact $\Delta L$ from `lift`; compare $|\mathrm{d}L - \Delta L|$ with $0.02 L_0$. (d) also compute the Mach number $V/340.3$ and the Prandtl–Glauert factor $1/\sqrt{1 - M^2}$.

## 6.
(a) $L_\rho = \tfrac12 V^2SC_L = L/\rho$; $\rho'(h) = -\rho_0 \cdot 4.2559 \cdot 2.2558 \times 10^{-5}(1 - 2.2558 \times 10^{-5}h)^{3.2559}$. (b) $\mathrm{d}L/\mathrm{d}t = L_\rho\,\rho'(h)\,\dot h$; finite difference of `lift(alpha, V, Wing(rho=rho(h0 + h_dot*t)))` in $t$. (c) $\dot\alpha = -L_\rho\rho'\dot h/L_\alpha$; the stall altitude is where $\tfrac12\rho(h)V^2SC_{L,\max} = L_0$: solve for $\rho$, invert the ISA formula. (d) Add $L_V\dot V$; compare its size with $L_\rho\rho'\dot h$; new $\dot\alpha = -(L_\rho\rho'\dot h + L_V\dot V)/L_\alpha$.
