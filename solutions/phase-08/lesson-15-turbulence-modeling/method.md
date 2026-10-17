# P08.15  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
1. Re_b = V(2h)/ν. 2. `duct_flow(4.0, gap)` and `duct_flow(4.0, gap, outer="nikuradse")`: Re_τ, u_τ, τ_w, dP/dx; Δp over 2 m = 2|dP/dx|. 3. Dean: c_f from `dean_cf`, τ_w = c_f ρV²/2; percentage differences. 4. y₁ = ν/u_τ, y₃₀ = 30ν/u_τ.

## 2.
(a) Integrate 0 = −dP/dx + dτ/dy from 0 to y; τ(0) = τ_w, τ(h) = 0 ⇒ dP/dx = −τ_w/h; nondimensionalise. CV: pressure force on the two ends (2h × Δp) against shear on two walls (2 × L × τ_w). (b) (1 + l⁺²g)g = τ⁺; quadratic; rationalise; limits l⁺ → 0 and l⁺ → ∞. (c) ν_t dU/dy = u_τ²; ν_t (dU/dy)² = ε; with dU/dy = u_τ/(κy) solve for ν_t, ε, then k from ν_t = C_μk²/ε; substitute in the ε equation; diffusion term with ε ∝ 1/y; set the sum to zero.

## 3.
(a) Integrate (1/κ)ln y⁺ + B from 0 to Re_τ, divide by Re_τ. Viscous layer: the real u⁺ is below the log law for y⁺ < 30 (lowers U_b⁺); wake: above it in the core (raises U_b⁺). (b) `brentq` on `log_law_friction(rt)[0] − Re_b`; then c_f; Dean. (c) Same with κ = 0.4327.

## 4.
(a) Loop over the cases with `mixing_length_profile`, `friction`, `fit_log_law`, `dean_cf`. (b) Rank by the range of c_f each produces. (c) Compare the spread with the correlation's scatter.

## 5.
(a) Force balance on a cylinder of radius r: πr²(−dP/dx) = 2πr τ(r); τ linear in r, hence in y = R − r. Flow rate ∫ U 2πr dr with r = R − y. (b) Profile from `mixing_length_profile(R⁺, outer=...)`; U_b⁺ with the area weight; Re_D, c_f; Blasius. (c) Compare errors; recall where Nikuradse's formula came from. (d) u⁺ at y⁺ = R⁺ over U_b⁺.

## 6.
(a) y₁ from Problem 1; n = ln(1 + (r − 1)h/y₁)/ln r for growth ratio r. (b) 30 ≤ y⁺ ≤ 0.1–0.2Re_τ; evaluate at Re_τ = 180. (c) Separation, secondary flow, curvature effects, relaminarisation. (d) The sensitivity band against the 2 % request.
