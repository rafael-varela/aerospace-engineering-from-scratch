# P01.08  Hints

One hint per problem. A hint points at the idea, not the answer. Read `method.md` next if the hint is not enough.

1. With $k_3 = 0$ the lift is *linear* in $\alpha$, so the total differential is exact for a pure $\alpha$ change; the only error comes from the $V^2$ and from the product of the two changes. Work in radians for $L_\alpha$ and convert to N/deg at the end.
2. Everything is a polynomial in $V$ (degree 2) and in $d = \alpha - \alpha_0$ (degree 3), so Taylor expansions terminate: write $L(\alpha + \Delta\alpha, V + \Delta V)$ out fully and collect powers. For (b), $\ln V = \tfrac12[\ln(2W/\rho S) - \ln C_L]$.
3. $V = \sqrt{2q/\rho}$ is a product of powers, so logarithmic derivatives are just the exponents: $\partial\ln V/\partial\ln q = \tfrac12$. Add relative uncertainties in quadrature with those weights. For (d), note that $L = qSC_L$ contains no $\rho$ at all.
4. The truncation term is $\tfrac16 h^2|f'''|$ and $L$ is cubic in $\alpha$ with $L_{\alpha\alpha\alpha} = -6k_3 qS$, a constant; $L$ is quadratic in $V$, so $L_{VVV} = 0$. The right-hand branch is rounding, $\varepsilon|L|/h$.
5. Use (8): compute the three second-order terms and the cubic $\tfrac16 L_{\alpha\alpha\alpha}\Delta\alpha^3$ and compare their sum with 2 % of $L_0$. For (c), evaluate the partials *at* 14°, not at 4°. For (d), the question is whether the *model* holds, not whether the linearization does.
6. $\partial L/\partial\rho = L/\rho$ (lift is linear in density). The chain rule for $L(\alpha(t), V(t), \rho(h(t)))$ has three terms; with $\alpha$ and $V$ constant only the density term survives. Holding $L$ constant means setting the sum of the terms to zero and solving for $\dot\alpha$.
