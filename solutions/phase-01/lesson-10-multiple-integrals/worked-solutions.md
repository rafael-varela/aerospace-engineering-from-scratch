# P01.10  Worked solutions

Worked for seed $s = 15$ where a seed is used. Every number below is reproduced by `solution.py` or by the short scripts quoted inline.

## 1. An untapered, unswept block

$L = 5.5$ m, $c = 1.2$ m.

(a) $w = 0.4 \times 1.2 = 0.48$ m, $h = 0.12 \times 1.2 = 0.144$ m, $m = 150 \times 0.48 \times 0.144 \times 5.5 = 57.02$ kg.

(b) $\mathbf r_{cg} = (0.4 \times 1.2,\ 2.75,\ 0) = (0.48, 2.75, 0)$ m.

(c) $I_{xx,cg} = 57.024(5.5^2 + 0.144^2)/12 = 57.024 \times 30.2707/12 = 143.85$; $I_{yy,cg} = 57.024(0.48^2 + 0.144^2)/12 = 1.193$; $I_{zz,cg} = 57.024(0.48^2 + 5.5^2)/12 = 144.84$ kg m²; products zero.

(d) $\mathbf d = (0.48, 2.75, 0)$, $|\mathbf d|^2 = 7.7929$. $I_{xx} = 143.85 + 57.024(2.75^2 + 0) = 143.85 + 431.24 = 575.09$; $I_{yy} = 1.193 + 57.024 \times 0.2304 = 14.33$; $I_{zz} = 144.84 + 57.024 \times 7.7929 = 589.23$; $I_{xy} = 0 - 57.024 \times 0.48 \times 2.75 = -75.27$ kg m²; $I_{xz} = I_{yz} = 0$.

$$I_O = \begin{pmatrix}575.09 & -75.27 & 0\\ -75.27 & 14.33 & 0\\ 0 & 0 & 589.23\end{pmatrix}\ \mathrm{kg\,m^2}.$$

Code: `mass_properties_analytic(WingBox(L=5.5, c_root=1.2, c_tip=1.2, f0=0.2, f1=0.6, tau=0.12, sweep=0.0, rho=150.0))` → $m = 57.024$, cg $= (0.48, 2.75, 0)$, $I_{xx} = 575.0905$, $I_{yy} = 14.3317$, $I_{zz} = 589.2252$, $I_{xy} = -75.2717$; `inertia_about_cg` → diag(143.8465, 1.1934, 144.8429) with zero products ✓.

## 2. The closed forms

(a) With $\delta = (c_t - c_r)/L$, $c = c_r + \delta y$ and $c^2 = c_r^2 + 2c_r\delta y + \delta^2y^2$. $\int_0^L c^2\,\mathrm{d}y = c_r^2L + c_r\delta L^2 + \delta^2L^3/3 = L[c_r^2 + c_r(c_t - c_r) + (c_t - c_r)^2/3] = L[c_r^2 + c_rc_t + c_t^2]/3$ after expanding $(c_t - c_r)^2 = c_t^2 - 2c_rc_t + c_r^2$ and collecting: $c_r^2 - c_r^2 + c_r^2/3 = c_r^2/3$, $c_rc_t - 2c_rc_t/3 = c_rc_t/3$, $c_t^2/3$ ✓. Multiply by $\rho(f_1 - f_0)\tau$ for (3). $c_t = c_r$: $m = \rho(f_1 - f_0)\tau Lc_r^2$, the block. $c_t = 0$: $m = \rho(f_1 - f_0)\tau Lc_r^2/3$, a third of the block, as for a pyramid-like solid whose cross-section shrinks quadratically.

(b) $\int_0^Lyc^2\,\mathrm{d}y = c_r^2L^2/2 + 2c_r\delta L^3/3 + \delta^2L^4/4 = L^2[c_r^2/2 + 2c_r(c_t - c_r)/3 + (c_t - c_r)^2/4] = L^2[c_r^2 + 2c_rc_t + 3c_t^2]/12$ (collect: $c_r^2(\tfrac12 - \tfrac23 + \tfrac14) = c_r^2/12$; $c_rc_t(\tfrac23 - \tfrac12) = c_rc_t/6$; $c_t^2/4$). Divide by (a): $y_{cg} = L(c_r^2 + 2c_rc_t + 3c_t^2)/[4(c_r^2 + c_rc_t + c_t^2)]$. $c_t = c_r$: $6/(4 \times 3) = 1/2$ ✓. $c_t = 0$: $1/4$ ✓.

(c) Partials, rows $(x, y, z)$, columns $(\xi, \eta, \zeta)$:

$$\frac{\partial(x, y, z)}{\partial(\xi, \eta, \zeta)} = \begin{pmatrix}(f_1 - f_0)c & L\tan\Lambda + [f_0 + (f_1 - f_0)\xi](c_t - c_r) & 0\\ 0 & L & 0\\ 0 & (\zeta - \tfrac12)\tau(c_t - c_r) & \tau c\end{pmatrix}.$$

Expand along the middle row: only the $L$ entry is non-zero, with cofactor $(f_1 - f_0)c \cdot \tau c - 0 = (f_1 - f_0)\tau c^2$. So $J = (f_1 - f_0)\tau Lc^2$. The sweep and taper terms sit in the $\eta$ column of the $x$ and $z$ rows, which multiply the zeros of the $y$ row's other entries: they describe how a slab at constant $\eta$ is displaced and sheared, and shearing preserves volume.

(d) $\mathbf r = \mathbf d + \mathbf r'$ with $\int\rho\mathbf r'\,\mathrm{d}V = \mathbf 0$. $|\mathbf r|^2 = |\mathbf d|^2 + 2\mathbf d\cdot\mathbf r' + |\mathbf r'|^2$ and $r_ir_j = d_id_j + d_ir_j' + r_i'd_j + r_i'r_j'$. Integrating with $\rho$, the terms linear in $\mathbf r'$ vanish, leaving $I_{ij} = \int\rho(|\mathbf r'|^2\delta_{ij} - r_i'r_j')\,\mathrm{d}V + m(|\mathbf d|^2\delta_{ij} - d_id_j)$, which is (10). The matrix $M = |\mathbf d|^2E - \mathbf d\mathbf d^{\mathsf T}$ gives $M\mathbf d = |\mathbf d|^2\mathbf d - \mathbf d|\mathbf d|^2 = \mathbf 0$, and for $\mathbf v \perp \mathbf d$, $M\mathbf v = |\mathbf d|^2\mathbf v$: eigenvalues $0, |\mathbf d|^2, |\mathbf d|^2$, times $m$. A point mass has no inertia about the axis through itself.

(e) $I_{xx} + I_{yy} - I_{zz} = \int\rho[(y^2 + z^2) + (x^2 + z^2) - (x^2 + y^2)]\,\mathrm{d}V = 2\int\rho z^2\,\mathrm{d}V \ge 0$, with equality only for a mass distribution entirely in the plane $z = 0$ (a lamina).

## 3. A partly filled spherical tank

(a) With $r_{\max}(z)^2 = R^2 - (z - R)^2 = 2Rz - z^2$ and $\mathrm{d}V = r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z$:

$$m = \rho\int_0^h\!\int_0^{2\pi}\!\int_0^{r_{\max}}r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z, \quad z_{cg} = \frac{\rho}{m}\int_0^h\!\int_0^{2\pi}\!\int_0^{r_{\max}}z\,r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z, \quad I_z = \rho\int_0^h\!\int_0^{2\pi}\!\int_0^{r_{\max}}r^2\cdot r\,\mathrm{d}r\,\mathrm{d}\theta\,\mathrm{d}z.$$

(b) Inner integrals: $\int_0^{2\pi}\int_0^{r_{\max}}r\,\mathrm{d}r\,\mathrm{d}\theta = \pi r_{\max}^2 = \pi(2Rz - z^2)$; $\int\int r^3 = \tfrac12\pi(2Rz - z^2)^2$. Then

$$V = \pi\int_0^h(2Rz - z^2)\,\mathrm{d}z = \frac{\pi h^2(3R - h)}{3}, \qquad z_{cg} = \frac{\pi\int_0^hz(2Rz - z^2)\,\mathrm{d}z}{V} = \frac{h(8R - 3h)}{4(3R - h)},$$

$$I_z = \frac{\rho\pi}{2}\int_0^h(2Rz - z^2)^2\,\mathrm{d}z = \frac{\rho\pi h^3(20R^2 - 15Rh + 3h^2)}{30}.$$

(c) $R = 0.6$ m, $\rho = 1141$ kg/m³:

| $h$ (m) | $V$ (m³) | $m$ (kg) | $z_{cg}$ (m) | $I_z$ (kg m²) |
|---|---|---|---|---|
| 0.3 | 0.14137 | 161.3 | 0.1950 | 15.39 |
| 0.6 | 0.45239 | 516.2 | 0.3750 | 74.33 |
| 0.9 | 0.76341 | 871.0 | 0.5250 | 133.27 |

Full tank $h = 2R$: $V = \pi 4R^2 R/3 = \tfrac43\pi R^3$ ✓; $z_{cg} = 2R(8R - 6R)/(4R) = R$ ✓; $I_z = \rho\pi 8R^3(20 - 30 + 12)R^2/30 = \tfrac{16}{30}\rho\pi R^5 = \tfrac25(\tfrac43\pi R^3\rho)R^2 = \tfrac25mR^2$ ✓. At $h = R$ the half-sphere has $z_{cg} = 5R/8 = 0.375$ m ✓ (the familiar centroid of a hemisphere measured from its flat face, here $R - 3R/8$).

(d) Tank axis parallel to the roll axis at $d = 2$ m: $I_{\text{roll}} = I_z + md^2$. At $h = 0.9$: $133.3 + 871.0 \times 4 = 133.3 + 3484.2 = 3617.5$; at $h = 0.3$: $15.4 + 645.2 = 660.6$ kg m². Change: $-2956.9$ kg m². The parallel-axis term dominates by $3484/133 = 26\times$ at full and $645/15 = 42\times$ at low fill; the tank's own inertia about its axis is almost irrelevant and what matters is *how much mass* sits 2 m off the axis. This is why vehicle inertia bookkeeping is mostly a mass-and-position exercise.

(e) The free surface is perpendicular to the net body-frame acceleration, not to gravity: during a burn with thrust acceleration $a$ along the axis and lateral acceleration $a_\perp$ the surface tilts by $\arctan(a_\perp/a)$, and when the acceleration changes the liquid sloshes with its own natural frequencies, which couple into the attitude control loop. In coasting flight there is no defined surface at all and the propellant is positioned by surface tension devices or a settling burn. The flat-surface, rigid-liquid model is a static approximation; a sloshing propellant is modelled as a pendulum or spring-mass of equivalent inertia (Phase 11 touches this).

## 4. Convergence of the inertia and Richardson extrapolation

(a) Signed error of midpoint $I_{xx}$ against SymPy 384.095232:

| $n$ | 8 | 16 | 32 | 64 | 128 |
|---|---|---|---|---|---|
| error (kg m²) | $+3.36 \times 10^{-3}$ | $-1.21 \times 10^{-4}$ | $-9.05 \times 10^{-5}$ | $-2.64 \times 10^{-5}$ | $-6.83 \times 10^{-6}$ |
| observed order | | 4.8 | 0.4 | 1.8 | 1.9 |

Split into parts: the $y^2$ part has errors $+5.1 \times 10^{-3}, +3.2 \times 10^{-4}, +2.0 \times 10^{-5}, +1.3 \times 10^{-6}$ (order 4, positive) and the $z^2$ part $-1.76 \times 10^{-3}, -4.42 \times 10^{-4}, -1.10 \times 10^{-4}, -2.76 \times 10^{-5}$ (order 2, negative). The one-dimensional midpoint error is $-\tfrac{h^2}{24}[g'(1) - g'(0)]$ per direction. For the $y^2$ part the mapped integrand's $\eta$-dependence is $g(\eta) = \eta^2c(\eta)^2$ with $c = c_r + (c_t - c_r)\eta$: $g' = 2\eta c^2 + 2\eta^2c(c_t - c_r)$, so $g'(0) = 0$ and $g'(1) = 2c_t^2 + 2c_t(c_t - c_r) = 2c_t(2c_t - c_r)$, which is zero exactly when $c_t = c_r/2$: our box. The $h^2$ term vanishes and the next term, $O(h^4)$, remains (positive here). For the $z^2$ part the $\zeta$-dependence $(\zeta - \tfrac12)^2$ has $g'(1) - g'(0) = 2$ and the $\eta$-dependence $c^4$ has $g'(1) - g'(0) = 4(c_t - c_r)(c_t^3 - c_r^3) > 0$ (the product of two negatives; numerically $7168/625 = 11.47$ with the lesson's chords), both giving negative errors of order 2. The total is positive while the $y^2$ term dominates ($n \le 8$) and negative once the slower $z^2$ term wins ($n \ge 16$); the crossing near $n = 16$ produces the meaningless "order 0.4". For comparison the mass integrand $c^2$ has $g'(1) - g'(0) = 32/25$ and the predicted error at $n = 32$, $-\tfrac{1}{24 \times 32^2}\times\rho(f_1 - f_0)\tau L \times 1.28 = -0.00234$ kg, is the observed one to three figures.

(b) $(4I_{32} - I_{16})/3 - I_{\text{exact}} = -8.0 \times 10^{-5}$ (worse than $I_{32}$ alone, $-9.0 \times 10^{-5}$, by nothing useful); $(4I_{64} - I_{32})/3 - I_{\text{exact}} = -5.0 \times 10^{-6}$ (against $-2.6 \times 10^{-5}$ for $I_{64}$ alone: a 5× gain). Richardson assumes the error is $Ch^2$ with a single $C$; between $n = 16$ and 32 the error is still a mixture of an $h^4$ term and an $h^2$ term of opposite sign and the extrapolation removes the wrong thing. By $n = 32$–$64$ the $h^2$ term dominates and the extrapolation works.

(c) `WingBox(c_tip=0.6)`: errors $+1.53, +0.358, +0.0880, +0.0219, +0.00547$ for $n = 4$–64, ratios 4.27, 4.07, 4.02, 4.00: a clean second order, positive throughout, because $g'(1) = 2 \times 0.6(1.2 - 1.6) < 0$ no longer vanishes. The sign change in (a) was a coincidence of the 2:1 taper, not a property of the rule.

(d) Timings on the reference machine (`mass_properties_midpoint`, ten integrals each): $n = 16$: 0.8 ms, $3.2 \times 10^{-7}$ relative on $I_{xx}$; 32: 4 ms, $2.4 \times 10^{-7}$; 64: 44 ms, $6.9 \times 10^{-8}$; 128: 0.57 s, $1.8 \times 10^{-8}$. `tplquad_mass_and_ixx` (two integrals): 0.34 s for $10^{-10}$ relative. Cost grows as $n^3$ for the grid and the error only as $n^{-2}$, so the grid's cost for a target error $\epsilon$ grows as $\epsilon^{-3/2}$: $10^{-10}$ would need $n \approx 10^4$ and hours. The adaptive routine is cheaper below about $10^{-8}$ relative; above that, the grid wins by a factor of 100 and vectorises over ten integrands at once. For mass properties, where the model error is percent-level, $n = 32$ is already far beyond what the inputs deserve.

## 5. How good does the tensor need to be?

(a) Half-wing about its own CG: $I_{xx} = 124.81$, $I_{yy} = 4.05$, $I_{zz} = 128.67$, $I_{xy} = -14.01$ kg m². Shift by $\mathbf d = (0, 1.9643, 0)$ to the aircraft CG: $I_{xx} \to 124.81 + 67.2 \times 3.8585 = 384.10$ (the same as about the root, since the shift is purely spanwise), $I_{yy} \to 4.05$, $I_{zz} \to 128.67 + 259.3 = 387.95$. Two mirrored halves: $I_{xx} = 768.2$, $I_{yy} = 8.1$, $I_{zz} = 775.9$, products cancel. Rod of 900 kg and 7 m along $x$: $I_{xx} = 0$, $I_{yy} = I_{zz} = 900 \times 49/12 = 3675.0$. Aircraft: $I_{xx} = 768.2$ (100 % wings), $I_{yy} = 3683.1$ (99.8 % rod), $I_{zz} = 4450.9$ (83 % rod, 17 % wings) kg m².

(b) A 15 % error in the wing's $I_{xx}$ is $\pm115$ kg m², 15 % of the aircraft's $I_{xx}$: the roll inertia *is* the wing model. Doubling the wing's $I_{yy}$ adds 8.1 to 3683.1: 0.2 %, invisible. (The 15 % also enters $I_{zz}$ through the wings' $I_{zz}$ share: 17 % × 15 % = 2.6 %.)

(c) Roll: worth a better wing model before Phase 09, because nothing else contributes and the smeared density misplaces mass along the span, which is what $I_{xx}$ weights. Pitch: no wing model will help; $I_{yy}$ is the fuselage, engine and tail, which were lumped as a uniform rod here, and only a swing test of the whole aircraft (or a proper mass breakdown of the fuselage) improves it. Yaw: in between; the wing's share is 17 % and its $I_{zz}$ depends on the same spanwise distribution as $I_{xx}$, so the better wing model fixes roll and most of yaw together.

(d) *The tensor handed over is diag(768, 3683, 4451) kg m² about the aircraft CG with body axes $x$ aft, $y$ starboard, $z$ up, products of inertia zero by symmetry. $I_{xx}$ comes entirely from the wing box model with a smeared density and is trusted to about 15 %; $I_{yy}$ and $I_{zz}$ come mostly from a uniform-rod fuselage and are trusted to no better than 30 %. The least trustworthy entry for its purpose is $I_{yy}$, because pitch dynamics are the first thing Phase 09 will tune and the rod is a placeholder, not a model. A component mass breakdown of the fuselage with positions, or a bifilar swing test of the assembled aircraft, would bring all three to 5 %.*

## 6. A density that varies along the span

(a) $A(y) = 0.05c(y)^2$ and $\rho(y) = 180(1 - y/10)$ for $L = 5$:

$$m = \int_0^5 180(1 - y/10) \times 0.05\,c^2\,\mathrm{d}y = 54.000\ \mathrm{kg}, \qquad y_{cg} = \frac{\int y\rho A}{m} = \frac{26}{15} = 1.7333\ \mathrm{m},$$
$$I_{xx} = \int_0^5\rho(y)A(y)\left[y^2 + \frac{(0.1c)^2}{12}\right]\mathrm{d}y = \frac{787\,752}{3125} = 252.081\ \mathrm{kg\,m^2}.$$

(SymPy: `integrate(180*(1 - y/10)*Rational(1,20)*c**2, (y, 0, 5))` with `c = Rational(8,5) - Rational(4,25)*y`.)

(b) `midpoint_3d(lambda x, y, z: 180*(1 - 0.5*y/5), n)` and the two weighted integrands:

| $n$ | $m$ | $y_{cg}$ | $I_{xx}$ |
|---|---|---|---|
| 8 | 53.9156 | 1.74473 | 252.322 |
| 16 | 53.9789 | 1.73618 | 252.139 |
| 32 | 53.9947 | 1.73405 | 252.095 |

Mass errors $-0.0844, -0.0211, -0.0053$: ratio 4.0 each time. $I_{xx}$ errors $+0.242, +0.059, +0.015$: ratio 4.1, 4.0. Second order, as the integrands are still polynomials in the cube coordinates.

(c) Uniform box: $m = 67.2$, $y_{cg} = 1.964$, $I_{xx} = 384.10$. Graded: the mass falls to 80.4 % (the mean of $1 - y/10$ weighted by $c^2$), $y_{cg}$ moves inboard by 0.23 m, and $I_{xx}$ falls to 65.6 %. $I_{xx}$ falls further because it weights the mass by $y^2$, and the density reduction is largest exactly where $y^2$ is largest: the outboard mass that the uniform model placed at large radius is the mass that is no longer there.

(d) Give `mass_properties_midpoint` an optional parameter `rho=None`; inside, `rho_fn = (lambda x, y, z: np.full_like(x, box.rho)) if rho is None else rho`, and multiply every integrand by `rho_fn(x, y, z)` instead of the scalar `box.rho`. Calls without the argument are unchanged, so every existing test passes; `mass_closed_form`, `spanwise_cg_closed_form` and `mass_properties_analytic` keep the uniform density and are not touched (their docstrings should say so). A new test would pass `rho=lambda x, y, z: 180*(1 - 0.1*y)` and check the $n = 32$ mass against 54.0 to $10^{-3}$. The change is six lines because the midpoint rule never assumed anything about the integrand beyond broadcasting, which is the practical reason to keep the quadrature routine general and the physics in the integrands.
