# P05.10  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. FBD: the arm, the tip load, the reactions (a force and a moment) at the support. M = P L. 2. I = bh³/12, c = h/2, S = I/c, σ = M/S. 3. MS_yield = F_ty/σ − 1; MS_ultimate = F_tu/(1.5σ) − 1. 4. Swap b and h: S′ = hb²/6; ratio σ′/σ = h/b.

## 2.
(a) ∫σ dA = −Eκ∫y dA = 0 ⇒ ∫y dA = 0 ⇒ origin at the centroid (equilibrium + constitutive). M = −∫σ y dA = EκI ⇒ κ = M/(EI). ε = −κy is compatibility (geometry). (b) Substitute σ = E(a + bx + cy) into the three conditions; with centroidal axes the first gives a = 0; the other two give E[[I_xy, I_xx], [I_yy, I_xy]]·(b, c) = −(M, M_h); solve by Cramer's rule; σ = E(bx + cy). Set σ = 0 for the neutral axis: y/x = −(M_h I_xx − M I_xy)/(M I_yy − M_h I_xy). (c) I = 2[A(h_c/2)² + At²/12]; σ_outer = M(h/2)/I; take the ratio of M/(A h_c) to σ_outer; simplify with h = h_c + t; let t → 0.

## 3.
(a) I = ∫(R sin φ)² Rt dφ over 0..2π = πR³t; compare with the exact annulus formula. (b) σ_N = N/A; σ_b = ±MR/I; total = σ_N ± σ_b; maximum compression on the side where bending compresses, maximum tension opposite, both in the plane of the bending moment. (c) Set N/A + MR/I = 0 with N < 0 and substitute the thin-wall A and I. (d) Thin shells buckle in compression at stresses far below yield.

## 4.
1. L = 3.8 m g / 2; M₀ = L · 4s/(3π) (or `elliptic_shear_moment`). 2. f(b_f) = 1.5 M₀ (h/2)/I(b_f) − F_tu with I from `section_properties(i_section(h, b_f, 0.010, 0.003))`; `brentq` on a bracket where f changes sign. A_cap = b_f × 0.010; two-flange: `cap_area_required(M0, h − 0.010, F_tu, 1.5)`. 3. σ_limit = F_tu/1.5; MS_yield = F_ty/σ_limit − 1. 4. `cantilever_fem(s, E I, 20, lambda y: elliptic_lift(y, L, s))`; compare E κ_mid[0] h/2 with M(x_mid[0]) (h/2)/I.

## 5.
(a) `section_properties(z_section(...))`, `neutral_axis_angle(M, 0, ...)`. (b) `bending_stress_general(M, 0, x, y, ...)` and `flexure_stress(M, y, I_xx)` at the four points. (c) With ε = cy: M = −∫Ecy·y dA = −EcI_xx ⇒ c = −M/(EI_xx) ⇒ σ = −My/I_xx; M_h = −∫Ecy·x dA = −EcI_xy = M I_xy/I_xx. Check with `bending_stress_general(M, M_h, ...)`. (d) Consider the stiffness of the skin and rivets in the lateral direction versus the stringer's lateral bending stiffness; bracket the answer by the two analyses.

## 6.
For each case compute the root moment or the net section and then σ = Mc/I and MS_ultimate. (a) M = Ls/2. (b) q = n m_w g/s, relief q s²/2. (c) Lift × 1.1. (d) Rectangles with the bottom cap reduced to 48 mm wide; `section_properties`; σ_bottom = M(c + y_c)/I. (e), (f): argue direction; (f) cannot be quantified with this lesson's tools; name the analysis (`P05.15`, plate buckling).
