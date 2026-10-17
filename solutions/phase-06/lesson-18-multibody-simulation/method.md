# P06.18  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) `gas_force` and `gas_energy` with `Gear(n_poly=...)` at the three strokes; isothermal check with the logarithm. (b) `brentq` on gas_force(s) − Mg cos β over 0–0.3 m; stiffness by a central difference (or the formula n A_a² p/V). (c) fₙ = (1/2π)√(k/M); compare with 65 Hz and think about what the strut must absorb.

## 2.
(a) x₂ = −(ℓ₀ − s) sin β, y₂ = y₁ − (ℓ₀ − s) cos β; T; M_r; δW for δy₁ and δs. (b) Integrate p₀V₀ⁿ(V₀ − A_a s′)⁻ⁿ A_a ds′; take n → 1 with the limit (xᵉ − 1)/ε → ln x. (c) Partial derivatives of each Φ; for Φ₄ differentiate d·n twice with ṅ = θ̇₁n′, ṅ′ = −θ̇₁n, move the q̈ terms to the left. (d) Power of constraint forces = λᵀΦ_q q̇ = λᵀΦ̇ = 0; so the energy balance needs only applied forces.

## 3.
(a) E = ½mv² + m g_Mars × stroke; F_total = E/stroke; divide by 4. (b) a = F/m − g_Mars; divide by 9.80665; load factor F/(m g₀). (c) η = 1 for constant force; linear spring: ½k x² = E gives k x_max = 2E/x_max; ratio 2. Think about reuse, mass and the number of landings.

## 4.
(a) Loop over A_o with `dataclasses.replace(Gear(), A_o=...)`, `simulate_reduced(..., t_end=0.3)`, `drop_metrics`, max |`energy_terms`[...]['residual']|; take the minimum. (b) Plot n against A_o; identify the damping-dominated and gas-dominated sides. (c) Map a ±0.1 C_d band to an equivalent A_o band (C_d A_o is what matters) and look at the slope of the curve there.

## 5.
(a) f(v) = peak stroke(v) − 0.30 for lift 1 and lift 0; `brentq` between 1 and 6 m/s. (b) Energy to absorb with lift ≈ ½(M + m)v²; without lift add (M + m)g × (vertical travel). (c) Ratio of sink rates and of energies.

## 6.
(a) Three runs with `replace(Gear(), Cd=0.6)`, `replace(Gear(), k_t=...)`; peak of Ft and its time from a fine output grid. (b) Compare the times with 16 ms. (c) Think about the stroke and stroke-rate records and the tyre deflection. (d) Consider what a tuned parameter does to every other prediction.
