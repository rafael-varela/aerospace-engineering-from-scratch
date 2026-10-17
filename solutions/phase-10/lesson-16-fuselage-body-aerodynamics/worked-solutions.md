# P10.16  Worked solutions

Each solution states: assumptions, units, the derivation or method, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15; numerical results come from `solution.py` and a probe script that calls it, given to the figures shown. The trainer's fuselage profile and the tail data are the lesson's stated assumptions.

## 1. An external fuel tank (s = 15: D = 0.55 m)

**Assumptions.** Potential flow, slender body, α = 4°, q = ½ × 1.225 × 55² = 1852.8 Pa.

(a) Vol = (π/6) × 3.0 × 0.55² = **0.47517 m³**; `body_volume` on a 4001-point cosine grid gives 0.47517 m³ (difference 5 × 10⁻⁸ m³). Wetted area **4.130 m²**.

(b) N = **0** (−7 × 10⁻¹⁴ N, round-off): the tank is closed. M = 2qα·Vol = **122.93 N·m** nose-up about the point 1.0 m behind the nose, and about every other point, because a pure couple has no reference point.

(c) L/D = 5.455: k₁ = 0.0520, k₂ = 0.9057, so the exact moment is **104.60 N·m**; slender-body theory over-predicts it by **17.5 %**.

**Independent check.** `slender_body_closed_form(0, Vol, 3.0, α, q, x_ref)` gives the same 122.93 N·m for any x_ref.

**Common mistakes.** Expecting a lift on a closed body; forgetting the radians.

## 2. Derivations

(a) The cross-flow w = V sin α ≈ Vα moves each section sideways relative to the air; the apparent mass of a circle is ρS, so a slab of air at station x carries vertical momentum ρS(x)w per unit length. The slab moves aft at V: in time dt it moves dx = V dt, and its momentum changes by ρw dS. The body exerts force per unit length dN/dx = ρwV dS/dx = **ρV²α dS/dx = 2qα dS/dx** (and receives the reaction).

(b) N = 2qα∫₀ᴸ S′ dx = 2qα(S(L) − S(0)) = **2qαS_base** (S(0) = 0 at a pointed nose). M = 2qα∫S′(x_ref − x)dx = 2qα{[S(x_ref − x)]₀ᴸ + ∫S dx} = **2qα(Vol − S_base(L − x_ref))**. With S_base = 0 the x_ref term disappears.

(c) S = πR²x/L, S′ = πR²/L: uniform loading 2qαπR²/L, N = 2qαπR², C_Nα = 2 on πR², centre of pressure at the centroid of a uniform load, **L/2**. Cone: S ∝ x², S′ ∝ x, loading triangular, centre of pressure at **2L/3** from the apex.

(d) For small e, artanh e = e + e³/3 + e⁵/5 + …: a₀ = (2(1 − e²)/e³)(e³/3 + …) → 2/3 and b₀ = 1/e² − ((1 − e²)/e³)(e + e³/3 + …) = 1/e² − 1/e² + 1 − 1/3 + O(e²) → 2/3, so **k₁ = k₂ = (2/3)/(4/3) = ½**. As L/D → ∞, e → 1: (1 − e²)artanh e → 0, so a₀ → 0 (k₁ → 0) and b₀ → 1 (k₂ → 1); then (6) gives q·Vol·sin 2α ≈ 2qα·Vol, which is (4).

**Common mistakes.** Taking the apparent mass of a circle as half the displaced mass (that is the sphere's 3-D value).

## 3. An airship hull

(a) D = 15 m, **Vol = 7068.6 m³**; **k₁ = 0.0816, k₂ = 0.8598**; q = 245 Pa, α = 5°: slender body **M = 302 257 N·m**, exact **M = 234 025 N·m**.

(b) dM/dα (exact, small α) = 2q·Vol(k₂ − k₁) = 2.695 × 10⁶ N·m/rad. Fin area for neutral stability: 2.695 × 10⁶/(245 × 3.0 × 25) = **146.7 m²**, about the projected area of a large house. Airships carry very large fins because the hull's Munk moment is proportional to its huge volume, while the fins' restoring moment grows only with their area and arm.

(c) Transverse apparent mass: exact k₂ρVol = **7445 kg**; slender body ρVol = **8659 kg**, equal to the mass of air displaced (8659 kg). An airship's own mass is about equal to the air it displaces (it floats), so its apparent mass is nearly as large as its real mass and changes its dynamics completely; an aeroplane is hundreds of times denser than the air it displaces.

**Common mistakes.** Using the slender-body moment for a fineness-4 body without the 23 % correction.

## 4. A longer fuselage (s = 15: L = 7.75 m)

(a) Volume **4.574 m³**, wetted area **20.12 m²**, fineness **7.05**.

(b) Spheroid-corrected slender body **Cm_α = 0.557 /rad**, neutral point **0.401 c**; cut loading **Cm_α = 0.216 /rad**, neutral point **0.459 c** (against 0.412 c and 0.457 c for the 7 m fuselage).

(c) **ΔC_D = 0.00507** (against 0.00490).

(d) The stretch is slightly destabilising in the potential model (more volume, so more Munk moment: the neutral point moves 0.011 c further forward) and costs 3.5 % more fuselage drag. In the cut model the change is negligible because the extra volume is mostly behind the maximum section. A real stretch usually also lengthens the tail arm, which raises the tail volume and more than compensates; here the tail position was kept, so the tail area would have to grow by about 5 % to recover the 0.012 c of margin (from (8): Δh_n ≈ (h_t − h_n) Δa_tail/Σa_i, so Δa_tail ≈ 0.012 × 5.32/(4.00 − 0.40) ≈ 0.018 /rad on a tail term of 0.36 /rad, about 5 %).

**Common mistakes.** Stretching the cabin but also scaling the diameter; moving the wing reference with the stretch.

## 5. How reliable is the optimum fineness?

(a) Optimum L/D = **8.79, 8.75, 8.71** at 30, 55, 90 m/s. The friction coefficient varies only logarithmically with Re, and the optimum is set by the balance between FF and S_wet, which does not depend on speed.

(b) With FF = 1 + 1.5/f^{1.5} + 7/f³: optimum **L/D = 6.15**; D/q = **0.0638, 0.0532, 0.0553 m²** at L/D = 3, 6.36, 10, against **0.1331, 0.0593, 0.0569 m²** with the lesson's correlation.

(c) Both correlations agree that bodies blunter than about L/D = 4 are costly, and that the minimum is flat between about 6 and 10 (a few per cent). They disagree on where in that range it lies, so the location of the optimum is not a reliable design output; the fineness is better chosen for volume, structure, tail arm and stability.

**Common mistakes.** Quoting the optimum to three figures as if it were physics.

## 6. Which neutral point would you sign off?

(a) Inside the bracket: the tail-cone boundary layer thickens and partly separates, so it carries some but not all of the potential negative loading. Left out of both: the wing's upwash ahead of the wing raises the local α of the forward fuselage (and its downwash lowers it behind), which increases the destabilising moment and moves the neutral point further forward.

(b) At 10° C_N is **42 % above** potential theory (0.482 against 0.341), at 20° **87 % above** (1.182 against 0.633); the viscous term exceeds the potential term above **22.7°**. At cruise (α of a few degrees) it adds about 8 % at 2°, small compared with the bracket of (a), but it matters for stall and spin.

(c) Circular cross-sections: a flat-sided section has a larger apparent mass in cross-flow than the circle of the same area, so Cm_α would be larger (more destabilising). A canopy adds local volume and area growth ahead of the wing, also destabilising.

(d) A wind-tunnel test of the fuselage alone and of the wing–body combination (and the complete model), measuring the pitching moment against α about the wing quarter chord; the slopes give Cm_α of the fuselage alone and with interference, and the complete-model slope gives the neutral point directly.

**Common mistakes.** Choosing the more favourable (aft) neutral point; for tail sizing the forward bound is the safe choice.

## Check-yourself answers

1. ρπr², the mass of fluid displaced per unit length.
2. dN/dx = 2qα dS/dx: positive where the body grows (the nose), negative where it shrinks (the tail cone).
3. The positive loading on the growing part and the negative loading on the shrinking part cancel because S returns to zero; it feels a pure nose-up couple 2qα·Vol.
4. A couple (two equal and opposite forces) has the same moment about every point; the net force is zero.
5. C_Nα = 2, based on the base (end) area.
6. A finite-fineness spheroid has k₂ < 1 and k₁ > 0, so its exact moment is q·Vol(k₂ − k₁) sin 2α < 2qα·Vol; at L/D = 6 the factor is about 0.87–0.88, i.e. slender-body theory is roughly 13–15 % high.
7. The drag of the cross-flow past the body (lee-side separation and vortices), proportional to sin²α, so it grows as α².
8. Forward (it is destabilising); the tail (area or arm) must be enlarged, or the centre-of-gravity range moved forward.
