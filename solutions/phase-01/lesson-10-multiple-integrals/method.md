# P01.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. $w = (f_1 - f_0)c$, $h = \tau c$, $m = \rho whL$. 2. CG at $(\tfrac12(f_0 + f_1)c,\ L/2,\ 0)$. 3. Block about its CG: $I_{xx} = m(L^2 + h^2)/12$, $I_{yy} = m(w^2 + h^2)/12$, $I_{zz} = m(w^2 + L^2)/12$, products zero. 4. $I_O = I_{cg} + m(|\mathbf d|^2E - \mathbf d\mathbf d^{\mathsf T})$ with $\mathbf d = \mathbf r_{cg}$: diagonal entries gain $m(d_y^2 + d_z^2)$ etc., $I_{xy}$ gains $-m\,d_xd_y$. 5. `mass_properties_analytic(box)`; `inertia_about_cg(I, m, cg)` should return the block tensor; `parallel_axis(I_cg, m, cg)` should return $I$.

## 2.
(a) $c^2 = c_r^2 + 2c_r(c_t - c_r)y/L + (c_t - c_r)^2y^2/L^2$; integrate termwise; simplify to $L(c_r^2 + c_rc_t + c_t^2)/3$. Set $c_t = c_r$ and $c_t = 0$. (b) $\int yc^2 = L^2(c_r^2 + 2c_rc_t + 3c_t^2)/12$ by the same expansion; divide by (a). (c) Rows $(x, y, z)$, columns $(\xi, \eta, \zeta)$: $x_\xi = (f_1 - f_0)c$, $x_\eta = L\tan\Lambda + [f_0 + (f_1 - f_0)\xi](c_t - c_r)$, $x_\zeta = 0$; $y_\eta = L$ only; $z_\eta = (\zeta - \tfrac12)\tau(c_t - c_r)$, $z_\zeta = \tau c$. Expand the determinant along the second row. (d) $|\mathbf r|^2\delta_{ij} - r_ir_j$ with $\mathbf r = \mathbf d + \mathbf r'$; the terms linear in $\mathbf r'$ integrate to zero. Eigen-decomposition of $|\mathbf d|^2E - \mathbf d\mathbf d^{\mathsf T}$: it annihilates $\mathbf d$ and acts as $|\mathbf d|^2$ on anything perpendicular. (e) $I_{xx} + I_{yy} - I_{zz} = \int\rho[(y^2 + z^2) + (x^2 + z^2) - (x^2 + y^2)]\,\mathrm{d}V = 2\int\rho z^2\,\mathrm{d}V$.

## 3.
(a) $\mathrm{d}V = r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z$; $m = \rho\int_0^h\int_0^{2\pi}\int_0^{r_{\max}(z)}r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z$; $z_{cg}$ has an extra $z$; $I_z$ has an extra $r^2$. (b) Inner two integrals: $\pi r_{\max}^2$ and $\tfrac12\pi r_{\max}^4$ with $r_{\max}^2 = 2Rz - z^2$; then polynomial integrals in $z$. (c) Substitute. Full tank: $h = 2R$. (d) $I_{\text{roll}} = I_z + m(2\,\mathrm m)^2$ for each fill; subtract. (e) Prose: surface is normal to the net body acceleration, not to gravity; sloshing; propellant management devices.

## 4.
(a) Loop `mass_properties_midpoint(n)[2][0, 0] - I_a[0, 0]`; order $= \log_2|e_n/e_{2n}|$. Coefficients: $g = \eta^2c(\eta)^2$, $g' = 2\eta c^2 + 2\eta^2cc'$; $g'(0) = 0$, $g'(1) = 2c_t^2 + 2c_tc' = 2c_t(c_t + c_t - c_r)$, zero when $c_t = c_r/2$. For the $z^2$ part, $\zeta$-direction: $g = (\zeta - \tfrac12)^2$, $g'(1) - g'(0) = 2 > 0$, so the error is negative. (b) $(4I_{2n} - I_n)/3$ on both pairs; compare with the plain $I_{2n}$ error. (c) `WingBox(c_tip=0.6)`; redo (a). (d) `time.perf_counter()` around each call; `plt.loglog(times, errors)`.

## 5.
(a) Half-wing tensor about its CG from `inertia_about_cg`; shift by $\mathbf d = (0, y_{cg}, 0)$ with `parallel_axis`; double it and zero $I_{xy}$ (the mirrored half has the opposite sign); rod: $I_{xx} = 0$, $I_{yy} = I_{zz} = m_r\ell^2/12$. (b) Multiply the wing contribution to each entry by 0.15 (for $I_{xx}$) or 1.0 (for $I_{yy}$, i.e. doubling) and divide by the aircraft total. (c)–(d) Prose from (b).

## 6.
(a) $m = \int_0^L\rho(y)A(y)\,\mathrm{d}y$, $y_{cg} = \int y\rho A/m$, $I_{xx} = \int\rho A[y^2 + (\tau c)^2/12]\,\mathrm{d}y$ with $A = (f_1 - f_0)\tau c^2$; SymPy `integrate`. (b) `midpoint_3d(lambda x, y, z: rho(y), n)`, `midpoint_3d(lambda x, y, z: y * rho(y), n) / m`, `midpoint_3d(lambda x, y, z: (y**2 + z**2) * rho(y), n)`. (c) Ratios to the uniform values; $I_{xx}$ weights $y^2$, where the density has fallen most. (d) Give `mass_properties_midpoint` an optional `rho` callable defaulting to `lambda x, y, z: box.rho`; multiply every integrand by it. Tests that use the default still pass; `mass_closed_form` and `mass_properties_analytic` are unaffected because they do not take the callable.
