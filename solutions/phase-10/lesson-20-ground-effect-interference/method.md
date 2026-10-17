# P10.20  Method outlines

The steps, without final numbers. Read `worked-solutions.md` only after you have an answer of your own.

## 1.
(a) Real vortices at (±b_v/2, h), images at (±b_v/2, −h) with opposite sense. (b) Real pair at the centre: 2 × Γ/(2π b_v/2) down. Image pair: magnitude Γ/(2πr) each, vertical fraction (b_v/2)/r, upward; ratio 1 − (b_v/2)²/r². Evaluate `single_horseshoe_trefftz_ratio`. (c) h/b_v = 0.1/(π/4); `wieselsberger_ratio(0.1)`; compare with the lattice; explain with the concentrated tip vortices.

## 2.
(a) Velocity of a 2-D vortex, w = Γ(y − y₀)/(2π r²) for a vortex along +x; evaluate both at z = 0. For 3-D segments apply the reflection to the Biot–Savart integrand. (b) Image downwash opposes the real one in the Trefftz plane. (c) Write tangency: Σ w_ij Γ_j = −V sin α; images reduce each |w_ij|.

## 3.
1. `wing_in_ground_effect(60, 7, None, radians(5))` and with h = 5, 10, 20. 2. CL ratio; K ratio; `wieselsberger_ratio(h/60)`. 3. Argue (b) from K and CL_α; (c) from flaps, sweep, dihedral, fuselage.

## 4.
1. h from your seed. 2. Free and ground solutions at α = 5°. 3. CL = mg/(½ρV²S); Di_free = ½ρV²S CL² (CDi/CL²)_free; Di_ge = K Di_free. 4. α = CL/CL_α in both. 5. ΔDi/m. 6. Repeat K for ns = 20, 80.

## 5.
1. `tail_downwash(lat, gamma, 0.3 + 4.5, z_t, 2.83)/α` and `downwash_gradient_far_field`. 2. Ground cases with `ground=True` and z_t = h. 3. a_t from `wing_in_ground_effect(2.83, 0.707, None, α)`; V_H = 2 × 4.5/(12 × 1.2); Δε = (dε/dα_free − dε/dα_ge) α; ΔCm = −V_H a_t Δε (Δε reduction raises the tail angle of attack).

## 6.
Argue each from the assumption list: thickness and camber images, venturi suction under the wing; ground-board boundary layer and moving belts; time-dependent images; non-flat moving surface.
