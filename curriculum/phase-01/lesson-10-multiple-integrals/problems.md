# P01.10  Problem set

Attempt every problem before opening `solutions/`. Record what you tried. State units on every numerical answer and the reference point and axes of every inertia value you write.

Several problems use your personal seed **s**: the day of the month you were born (1–31).

## 1. (Warm-up: an untapered, unswept block)

A rectangular wing box has half-span $L = 4 + 0.1s$ m, constant chord $c = 1.2$ m, walls at $0.2c$ and $0.6c$, depth $0.12c$, no sweep, and effective density 150 kg/m³. Model it as `WingBox(L=L, c_root=1.2, c_tip=1.2, f0=0.2, f1=0.6, tau=0.12, sweep=0.0, rho=150.0)`.

Find by hand: (a) the box's width $w$, height $h$ and mass; (b) the centre of mass; (c) the inertia tensor about the CG from the rectangular-block formulas $I_{xx,cg} = m(L^2 + h^2)/12$ etc.; (d) the tensor about the root leading edge by the parallel-axis theorem, including the product $I_{xy}$. Check (a)–(d) with `mass_properties_analytic`, `inertia_about_cg` and `parallel_axis`. Answer form: numbers with units, four significant figures; the full $3 \times 3$ tensor for (d).

## 2. (Derivation: the closed forms)

For the tapered box of the lesson with $c(y) = c_r + (c_t - c_r)y/L$:

(a) Derive $\int_0^Lc^2\,\mathrm{d}y = L(c_r^2 + c_rc_t + c_t^2)/3$ and hence the mass (3). Show it reduces to the untapered and triangular cases.
(b) Derive $y_{cg}$, equation (4), and show it equals $L/2$ for $c_t = c_r$ and $L/4$ for $c_t = 0$.
(c) Derive the Jacobian (7) of the map (6) by computing the full $3 \times 3$ matrix of partials, and explain in one sentence why the sweep and taper terms drop out.
(d) Starting from (2), derive the parallel-axis theorem (10) for the full tensor, and show that the added matrix has eigenvalues $0, m|\mathbf d|^2, m|\mathbf d|^2$ with eigenvector $\mathbf d$ for the zero.
(e) Show that for any mass distribution $I_{xx} + I_{yy} \ge I_{zz}$, and find the condition for equality.

Answer form: derivations.

## 3. (Unfamiliar scenario: a partly filled spherical tank)

A spherical oxidiser tank of radius $R = 0.6$ m contains liquid oxygen ($\rho = 1141$ kg/m³) to a depth $h$ above the bottom, under $1g$ so the surface is flat. Use cylindrical coordinates about the vertical axis with $z$ measured from the bottom of the tank, so that the liquid occupies $0 \le z \le h$, $0 \le r \le \sqrt{2Rz - z^2}$.

(a) Set up the triple integrals for the mass, the height $z_{cg}$ of the centre of mass and the moment of inertia $I_z$ about the vertical axis, with the $r$ from the cylindrical volume element.
(b) Evaluate all three in closed form as functions of $h$ (integrate over $r$ and $\theta$ first).
(c) Tabulate mass, $z_{cg}$ and $I_z$ for $h = 0.3, 0.6, 0.9$ m. Check the full-tank limit $h = 2R$ against $\tfrac25mR^2$ and $z_{cg} = R$.
(d) The tank sits 2 m from the vehicle's roll axis. By how much does $I_{\text{roll}}$ of the vehicle change as the tank empties from $h = 0.9$ to $h = 0.3$ m? Which term, the tank's own $I_z$ or the parallel-axis term, dominates, and by what factor?
(e) Explain why the flat-surface assumption fails during a burn and what replaces it (one paragraph).

Answer form: integrals in (a), closed forms in (b), a table for (c), numbers for (d).

## 4. (Uses your code: convergence of the inertia and Richardson extrapolation)

With the lesson's `DEFAULT_BOX`:

(a) Tabulate the signed error of the midpoint $I_{xx}$ against SymPy for $n = 8, 16, 32, 64, 128$ and the observed order between successive rows. Explain the sign change with the error coefficient of §2.5 (compute $g'(1) - g'(0)$ for $g(\eta) = \eta^2c(\eta)^2$ and for the $z^2$ part).
(b) Apply Richardson extrapolation $(4I_{2n} - I_n)/3$ to the $(16, 32)$ and $(32, 64)$ pairs. Which pair benefits, and why does the other not?
(c) Change the tip chord to 0.6 m and repeat (a). Is the order now a clean 2?
(d) Time `mass_properties_midpoint(n)` for $n = 16, 32, 64, 128$ and `tplquad_mass_and_ixx`. Plot error against time for the midpoint rule and mark `tplquad`. At what accuracy does the adaptive routine become the cheaper choice?

Answer form: two tables, one plot, a paragraph for (a) and (d).

## 5. (Judgment: how good does the tensor need to be?)

Phase 09 will need the inertia tensor of a 1200 kg aircraft whose two half-wings are the box of this lesson (mirrored), plus a fuselage, tail and engine that you may take as a uniform rod of 900 kg and length 7 m along $x$ through the CG.

(a) Assemble the aircraft's $I_{xx}$, $I_{yy}$ and $I_{zz}$ about its CG (take the aircraft CG at the wing root, $y = 0$, $x = x_{cg,\text{wing}}$, $z = 0$). Which component dominates each entry?
(b) The smeared-density wing model's $I_{xx}$ could be wrong by 15 % and its $I_{yy}$ by a factor of 2. Propagate each into the aircraft's tensor. Which aircraft entry is sensitive to the wing model and which is not?
(c) Roll response depends on $I_{xx}$, pitch on $I_{yy}$, yaw on $I_{zz}$. For which axis is it worth doing a better wing model before Phase 09, and for which would a swing test of the whole aircraft be the only useful improvement?
(d) Write five sentences for the flight-dynamics group: the tensor you are handing over, its reference point, the entry you trust least and why, and what would fix it.

Answer form: numbers with units for (a)–(b), prose for (c)–(d).

## 6. (Unfamiliar scenario: a density that varies along the span)

A real wing box is heavier at the root (thicker skins, fuel) than at the tip. Replace the uniform density by $\rho(y) = 180(1 - 0.5y/L)$ kg/m³ for the lesson's box.

(a) Write the mass, $y_{cg}$ and $I_{xx}$ as single integrals in $y$ using the cross-section area, and evaluate them in closed form (SymPy is fine; the integrands are polynomials).
(b) Evaluate the same three quantities with `midpoint_3d` and an integrand that includes $\rho(y)$, for $n = 8, 16, 32$, and confirm second-order convergence to (a).
(c) Compare with the uniform box: mass, $y_{cg}$ and $I_{xx}$ about the root. Why does $I_{xx}$ fall by a larger fraction than the mass?
(d) The lesson's closed forms (3)–(4) and `mass_closed_form` cannot be used here. What minimal change to `mass_properties_midpoint` would let it accept any `rho(x, y, z)`? Describe (or implement) it, and say which tests would still pass unchanged.

Answer form: closed forms and numbers for (a)–(c), a design note for (d).
