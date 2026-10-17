# P08.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) h = 1/32; U h Re; −2U/h; repeat U h Re for N = 17, 65, 129 and compare with 2. (b) ψ_yy at y = 0 is 2x²(1 − x)²; ω = −ψ_yy; build ψ on each grid, call `wall_vorticity(psi, U=0)` and read the node at x = 0.5; errors, then log₂ of successive ratios. (c) The wall error is O(h) but lives in a layer of one cell; its effect on the interior is a boundary-layer-sized perturbation whose integrated effect scales like h × h.

## 2.
(a) ∂/∂x of the v-momentum equation minus ∂/∂y of the u-momentum equation; collect u ω_x + v ω_y + ω(u_x + v_y); use continuity; the pressure's mixed derivatives cancel. (b) Substitute; mixed partials commute. (c) Taylor series of ψ normal to the wall to second order; ψ_w = 0, ψ_n = U_w, ψ_nn = −ω_w because the tangential second derivative of a constant is zero; the next term gives the order. (d) D_x = D⊗I, D_y = I⊗D; (D⊗I)(I⊗D) = D⊗D = (I⊗D)(D⊗I); the zero edge values make the wall velocities in `velocities` equal to the differenced ones. Flux: sum the trapezoid weights times centred differences; telescoping leaves two end terms; Taylor-expand ψ₁ and ψ_{N−2}.

## 3.
Four Newton solves; three tables (ψ_min and centre, u centre); `observed_order` on (33, 65, 129) and on (17, 33, 65); the GCI; the distance between centres on 65 and 129; the 65 value against the extrapolate as a percentage.

## 4.
(a) Run it; read the residual history. (b) Loop over Re with `initial=previous`. (c) Same loop on 33 and 129; compute U h Re at 1000 on 33. (d) Ratio of the two ψ_min values; the difference against the 100-case GCI.

## 5.
(a) Time the two functions; read `steps` and `dt`. (b) Evaluate the two limits on each grid. (c) Steps ∝ 1/dt at constant settling time; cost per step from the 65 run scaled by the Poisson solve size. (d) Max |ψ_march − ψ_Newton| against the marching stopping tolerance and the discretisation error.

## 6.
(a) Read the two entries of `omega` on each grid; ratio of successive corner values. (b) Compare the code's GCI on 129 with the 129 − Ghia difference and with the extrapolated − Ghia difference. (c) Separate what the code models (2-D, steady, laminar, a square) from the real well (3-D, unsteady, turbulent, an incoming boundary layer).
