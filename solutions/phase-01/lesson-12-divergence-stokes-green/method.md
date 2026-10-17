# P01.12  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Apply (1) and (2) of the lesson term by term; note that the s terms cancel in the divergence of **v** for every s. Check: `divergence_sym([s*x, -y, (1-s)*z])`. (b) If ∇ · **F** = 0 then by the divergence theorem the flux through any closed surface is zero; no integral needed. (c) On x = 1 the flux is s × 1 × 1; on x = 0 it is 0; the y and z faces carry no flux because G_y = G_z = 0. Volume integral: ∫₀¹ 2sx dx = s. The midpoint rule is exact for the linear integrand 2sx, so the residual is rounding.

## 2.
(a) Q term: ∬ ∂Q/∂x dA = ∫_c^d [Q(b, y) − Q(a, y)] dy = ∮ Q dy (right side up, left side down, top and bottom have dy = 0). P term: ∬ ∂P/∂y dA = ∫_a^b [P(x, d) − P(x, c)] dx = −∮ P dx (bottom traversed with dx > 0, top with dx < 0). Subtract. (b) d/dt ∭ ρ dV = −∯ ρ**v** · **n** dA; divergence theorem on the right; d/dt inside the fixed volume; the integrand must vanish because V is arbitrary. (c) Control volume = stream tube segment; side walls have **v** · **n** = 0; steady ⇒ ∯ ρ**v** · **n** dA = 0 ⇒ outflow at section 2 equals inflow at section 1.

## 3.
(a) ∇ · **g** = −(4/3)πGρ × 3. Integrate over V: −4πGρV = −4πG M. (b) `box_flux(lambda X, Y, Z: uniform_sphere_gravity(X, Y, Z, 5515.0), ((3.0e6, 3.001e6), (0, 1e3), (0, 1e3)), 4)`; divide by −4πG; compare with ρ × 10⁹ m³. (c) g(R) = (4/3)πGρR ⇒ ρ = 3g/(4πGR). (d) g(r) = (4/3)πGρr, linear in r; the real Earth has a dense core, so more mass is enclosed at a given radius than the uniform model allows.

## 4.
(a) `res = divergence_theorem_convergence(f, d, box, [2, 4, 8, 16, 32, 64])`; slope from `np.polyfit`; `tplquad(lambda z, y, x: d(x, y, z), -0.2, 0.2, -0.2, 0.2, 1.0, 2.0)`. (b) For fixed (y_j, z_k), the volume sum of ∂F_x/∂x over the x column is Σ_i ∂_x F_x(x_i) Δx with x_i the cell centres; the face contribution is F_x(x₁) − F_x(x₀). If F_x is at most quadratic in x, ∂_x F_x is at most linear and the midpoint rule integrates it exactly, so the two agree. For x³ the derivative 3x² is quadratic and the midpoint rule has error Δx²/… per cell. (c) Loop over n; errors |circ − Γ| and |area − Γ|; slopes; the area error depends on how the grid cuts the circle r = a.

## 5.
(a) Loop over n calling `stokes_theorem_check(f, c, 1.3, n)`. (b) Difference cap − disc for each n; ratio between successive n; identity ∇ · (∇ × **F**) = 0 with the divergence theorem. (c) (∇ × **v**) · **n** = 2ωz/R on the hemisphere; ∬ z dA in spherical coordinates with z = R cos s, dA = R² sin s ds dt. (d) The normal reverses, so the cap flux changes sign; detect by comparing with the disc or the rim (they must agree), or by checking that **r** · **n** > 0 on the surface.

## 6.
For each case decide: ∂ρ/∂t zero or not; ρ uniform or not; surface impermeable or not. Volume flux balances only if ρ is uniform and steady; mass flux balances if steady and no leak.

## 7.
(a) ∇ · (**r**/3) = 1 ⇒ ∯ (**r**/3) · **n** dA = ∭ 1 dV = V. ∇ · (x²/2, 0, 0) = x ⇒ ∯ (x²/2) n_x dA = ∭ x dV = V x_c. (b) Define `ell(s, t)` and `third = lambda X, Y, Z: (X/3, Y/3, Z/3)`; `surface_flux(third, ell, (0, pi), (0, 2 pi), n, n)`. (c) Shift the parametrisation by +1 in x; flux of (X²/2, 0, 0); divide by V. (d) For a flat triangle **n** dA = ½(**p**₂ − **p**₁) × (**p**₃ − **p**₁) and **r** · **n** is constant on the plane, equal to **p**₁ · **n**; so the contribution is (1/6) **p**₁ · [(**p**₂ − **p**₁) × (**p**₃ − **p**₁)] = (1/6) **p**₁ · (**p**₂ × **p**₃). Reversing the vertex order flips the sign of that triangle's term.
