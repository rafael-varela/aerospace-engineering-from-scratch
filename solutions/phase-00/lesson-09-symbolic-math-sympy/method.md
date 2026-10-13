# P00.09  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. $k = \tfrac12\rho C_D A$. 2. $v_t = \sqrt{mg/k}$. 3. $\tau = v_t/g$. 4. $t_{99} = \tau\operatorname{artanh}(0.99)$. 5. $x_{99} = (v_t^2/g)\ln\cosh(\operatorname{artanh}0.99) = v_t\tau\,\ln\cosh(2.6467)$. 6. Check (b) by `float(terminal_velocity_expr().subs({m: ..., g: ..., k: ...}))` and (e) by `lambdify_position()(t99, m, g, k)`.

## 2.
(a) 1. Define $v_t^2 = mg/k$ and rewrite the right side as $(g/v_t^2)(v_t^2 - v^2)$. 2. Separate: $\mathrm{d}v/(v_t^2 - v^2) = (g/v_t^2)\,\mathrm{d}t$. 3. Integrate from $0$ to $v$ and $0$ to $t$, using the $\operatorname{artanh}$ antiderivative, valid for $v < v_t$ (state that a body released from rest never exceeds $v_t$ because $\mathrm{d}v/\mathrm{d}t > 0$ for $v < v_t$ and $= 0$ at $v_t$). 4. Invert. (b) Substitute $u = gt/v_t$ and integrate $\tanh u$ to $\ln\cosh u$; the constant is zero because $\cosh 0 = 1$. (c) Write $\cosh u = \tfrac12 e^u(1 + e^{-2u})$ and expand the logarithm. (d) `sp.dsolve(sp.Eq(sp.Derivative(V(t), t), g - k*V(t)**2/m), V(t), ics={V(0): 0})`; compare with `sp.simplify(sol.rhs - velocity_expr())` or by evaluating both at a few points.

## 3.
(a) Both forces point down when the body moves up: $m\,\mathrm{d}v/\mathrm{d}t = -mg - kv^2$. Note that replacing $v \to -v$ in the lesson's equation gives $-mg + kv^2$, not this. (b) Separate: $\mathrm{d}v/(v_t^2 + v^2) = -(g/v_t^2)\,\mathrm{d}t$; integrate to $\arctan$. In SymPy: build the candidate, compute `sp.simplify(sp.diff(cand, t) + g + k*cand**2/m)`. (c) Apex when $v = 0$: $t_a = \tau\arctan(v_0/v_t)$. Height: integrate $v$ or use $\mathrm{d}(v^2)/\mathrm{d}x = -2g(1 + v^2/v_t^2)$ to get $h_a = (v_t^2/2g)\ln(1 + v_0^2/v_t^2)$. Substitute numbers; compare with $v_0/g$ and $v_0^2/2g$. (d) After the apex the body is at rest and falling, which is exactly the lesson's initial condition.

## 4.
(a) Loop over `h`, store both checks, `ax.loglog`. Fit `np.polyfit(np.log10(h), np.log10(err), 1)` separately for $h \ge 10^{-3}$ and $h \le 10^{-5}$. (b) `v3 = sp.diff(velocity_expr(), t, 3).subs({...}).subs(t, 3.0)`; `h_opt = (3*eps*abs(v)/abs(v3))**(1/3)`. (c) Write a `five_point(f, x, h)` function; repeat; the truncation term is now $h^4 f^{(5)}/30$.

## 5.
(a) Divide the series by $gt$ and keep the first correction. (b) Set it equal to 0.01 and solve for $t_1$; evaluate `1 - vel(t1, m, g, k)/(g*t1)`. (c) For each ball compute $k$, $v_t$, $t_1$; then the drag-model fall time from $x(t) = 2$ m: $t = \tau\operatorname{arcosh}(e^{gx/v_t^2})$; compare with $\sqrt{2x/g}$ and with the 0.01 s timing resolution.

## 6.
(a) Density ratio and its effect on $v_t \propto \rho^{-1/2}$; constant-$k$ assumption gone; no closed form with $\rho(x)$. (b) $v_t \propto (C_D A)^{-1/2}$: compute the spread. (c) Quadratic-model $v_t$ for the droplet, then $Re$; if $Re < 1$ switch to Stokes and recompute; compare the two $v_t$. Repeat for the raindrop. (d) Solve $(R/(R+h))^2 = 0.999$ for $h$; compare with any realistic drop.

## 7.
(a) Linear ODE: integrating factor or separation; $v = v_t(1 - e^{-t/\tau})$ with $v_t = mg/c$, $\tau = m/c = v_t/g$. (b) Set $mg/c = \sqrt{mg/k}$ to fix $c$. Solve $\tanh(u) = 0.9$ and $1 - e^{-u} = 0.9$ for $u = gt/v_t$. (c) `sp.lambdify` both; evaluate at $t = 2$ s. (d) Think about which quantity differs most between the two curves at the times a 2 m or 10 m drop can reach, and what timing or video resolution resolves that difference.
