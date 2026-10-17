# P05.11  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. V = 2000 + 100s. A = bh. τ_avg = V/A; τ_max = 1.5V/A; τ(h/4) from (2): 1.5V/A(1 − 4(h/4)²/h²) = (9/8)V/A. 2. Pin: A = πr², τ_avg = V/A, τ_max = 4V/(3A). 3. ∫_{−h/2}^{h/2} (3V/2bh)(1 − 4y²/h²) b dy = (3V/2h)(h − h/3) = V.

## 2.
(a) Net axial force on the part = ∫(dM y/I) dA = dM Q/I; balance with τt dz; dM/dz = −V gives |τ| = VQ/(It). Equilibrium (both steps) and the flexure stress (which itself used compatibility and the constitutive law). (b) Half-flange from the tip: x along the flange, y = ±h_c/2: q = −(S_y/I) t_f (h_c/2) s, magnitude linear to S_y t_f (b/2)(h_c/2)/I at the junction. Web from the top junction with q₀ = sum of the two half-flange flows and y = h_c/2 − s: q = q₀ − (S_y/I) t_w (h_c s/2 − s²/2). Integrate over the web length; use I = t_w h_c³/12 + 2bt_f(h_c/2)². Flanges are horizontal: their flow has no vertical component. (c) F = ½ b × S_y t_f b h/(2I); F h = S_y e; substitute I = t_w h³/12 + bt_f h²/2; limits.

## 3.
(a) All flows act along the legs; moment about the junction is zero ⇒ the shear centre is at the junction. Confirm with `shear_centre(segments)` with segments from free end to free end through the junction. (b) Segments ((a,0),(0,0),t), ((0,0),(0,a),t); `thin_wall_properties`; torque = P × (x_c − x_s) for a vertical load through the centroid. (c) I_xy ≠ 0 about axes parallel to the legs ⇒ unsymmetric bending; use the general formula of `P05.10`, or principal axes.

## 4.
1. V = 3.8 m g/2; I and the stresses from `i_section_shear_stress(V, h, b_f, t_f, t_w, y)` at y = 0 and y = h_w/2; web share = ∫τ t_w dy over the web / V. 2. `open_section_shear_flow(i_beam_segments(h − t_f, b_f, t_f, t_w), 0, V)`, peak |q| in the web segment / t_w. 3. Q_cap = b_f t_f (h − t_f)/2; q = VQ_cap/I; pitch = 2 × `fastener_pitch(q, 4000, 1.5)`. 4. MS = F_su/(1.5 τ_NA) − 1.

## 5.
(a) `channel_segments(h, 0.030, 0.0016, 0.0016)`, `shear_centre`, `channel_shear_centre_closed_form`, `thin_wall_properties`. (b) T = P × e. (c) Build the lipped list; repeat (a). (d) Hinge on the shear-centre line versus the aerodynamic and mass-balance requirements and the space available.

## 6.
(a) τ = V/(b t); MS = F_su/(1.5τ) − 1. (b) k_s with b/a = 0.8; τ_cr formula; ratio. (c) Diagonal tension: compressive loads on the stiffeners, inward pull and bending on the caps, rivet tension; list the other limits. (d) Choose t where τ/τ_cr at limit load is below 1 (no buckling at limit) or require a tension-field analysis for thinner webs.
