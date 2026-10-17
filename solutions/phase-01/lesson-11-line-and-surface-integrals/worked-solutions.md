# P01.11  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check (limiting case, dimensional analysis, or alternative method), and the common mistakes. Seeded problems are worked for s = 15. Every number below is reproduced by `solution.py` with the stated n.

## 1. One vortex, three loops

**Assumptions.** Point vortex, counter-clockwise, Γ = 10 + 15 = 25 m²/s; potential flow.

**Derivation.** (a) Speed Γ/(2πr) = 25/(4π) = 1.989 m/s at r = 2 m. Direction (−y, x)/r: at (2, 0) it is (0, 1), at (0, −2) it is (1, 0). So **v**(2, 0) = (0, 1.989) m/s and **v**(0, −2) = (1.989, 0) m/s: the flow goes up on the right of the vortex and to the right below it, which is counter-clockwise. ✓

(b) The 3 m circle encloses the vortex: circulation 25 m²/s. The 1 m circle at (5, 0) does not: 0. The flux of a vortex is zero through both (no radial velocity anywhere).

(c) With n = 2000: circulation 25.0000206 (error 2.06 × 10⁻⁵, predicted Γπ²/(3n²) = 2.06 × 10⁻⁵ ✓); off-centre circle −1 × 10⁻¹⁶ (rounding); fluxes −9 × 10⁻¹⁸ and 0.

**Result.** 1.989 m/s along +y and +x respectively; 25, 0, 0, 0 m²/s.

**Check.** Dimensions: (m²/s)/(m) = m/s. The relative error 8 × 10⁻⁷ is exactly π²/(3 · 2000²) = 8.2 × 10⁻⁷.

**Common mistakes.** Forgetting the 2π (a factor of 6.3 in the speed); reporting the velocity at (0, −2) as (−1.989, 0) (the sign of −y/r² with y = −2 is positive); expecting a non-zero flux because the field is "strong".

## 2. Circulation of a gradient field

**Derivation.** (a) For x > 0, θ = arctan(y/x). ∂θ/∂x = [1/(1 + y²/x²)] · (−y/x²) = −y/(x² + y²) = −y/r². ∂θ/∂y = [1/(1 + y²/x²)] · (1/x) = x/r². On the other half-planes θ differs from arctan(y/x) by a constant (±π), which does not change the derivatives, and on the y axis the same formulas follow from θ = π/2 − arctan(x/y). Hence ∇(Γθ/2π) = (Γ/2π)(−y, x)/r² = **v** for all (x, y) ≠ (0, 0). ✓

(b) Let C be a closed loop not enclosing the origin. Because the origin is outside the region bounded by C, there is a ray from the origin that does not meet C (take any ray through a point of the exterior region's unbounded part, far from C). Remove that ray from the plane: on what remains, θ can be chosen continuous and single-valued, so φ = Γθ/2π is a genuine potential along all of C. The fundamental theorem gives ∮_C ∇φ · d**r** = φ(end) − φ(start) = 0 because the loop ends where it starts. The non-enclosing assumption is used exactly once: to find a ray that misses C, which is impossible when the origin is inside C (every ray from an interior point crosses the boundary).

(c) Follow θ along a counter-clockwise loop that encloses the origin once: θ increases continuously and returns to its starting direction having gained 2π, so φ has gained Γ, and ∮ ∇φ · d**r** = Γ. Alternatively, cut the annular region between C and a small circle C_ε centred on the origin along a ray; the boundary of the cut region (C counter-clockwise, C_ε clockwise, the two edges of the cut in opposite directions) is a loop not enclosing the origin, so by (b) its circulation is zero; the cut edges cancel; hence Γ_C = Γ_{C_ε} = Γ by the direct calculation of lesson §2.2.

(d) Figure-of-eight: left lobe counter-clockwise around the vortex contributes Γ, right lobe encloses nothing and contributes 0: total Γ. Loop winding twice: 2Γ. In general the circulation is Γ times the winding number of the loop about the vortex.

**Check.** The argument never used the circular shape: a square loop gives Γ too, which is what the code finds in problem 1 and the tests.

**Common mistakes.** Concluding "the field is a gradient, so every closed-loop integral is zero" without checking that the potential is single-valued on the loop; differentiating arctan(y/x) and forgetting the inner derivative −y/x².

## 3. A hurricane as a Rankine vortex

**Assumptions.** Axisymmetric, solid-body core of radius a = 30 km, 1/r wind outside; the real wind profile decays faster than 1/r far out and the storm is neither axisymmetric nor steady.

**Derivation.** (a) Γ = 2πaV_max = 2π × 30 000 × 50 = 9.42 × 10⁶ m²/s. Core vorticity 2V_max/a = 100/30 000 = 3.33 × 10⁻³ s⁻¹.
(b) At 15 km (inside): Γ(15/30)² = 2.36 × 10⁶ m²/s. At 60 km (outside): 9.42 × 10⁶ m²/s. Wind at 100 km: Γ/(2π × 10⁵) = 15.0 m/s.
(c) Γ ratio hurricane/wake = 9.42 × 10⁶/300 = 3.1 × 10⁴. Vorticity: wake 2 × 31.8/1.5 = 42.4 s⁻¹ against 3.33 × 10⁻³ s⁻¹: the wake core spins 1.3 × 10⁴ times faster. Circulation scales with speed × size, vorticity with speed/size: the hurricane is 20 000 times bigger and 1.6 times slower.
(d) Code (n = 2000): 2 356 191 and 9 424 786 m²/s, matching (b) to 10⁻⁶.

**Result.** Γ = 9.42 × 10⁶ m²/s; ω = 3.3 × 10⁻³ s⁻¹; 2.36 × 10⁶ and 9.42 × 10⁶ m²/s; 15 m/s.

**Check.** 15 m/s at 100 km is a strong breeze, consistent with the gale-force radius of a real hurricane being of order 100–200 km. Dimensions: 2πaV_max is m × m/s. ✓

**Common mistakes.** Using a in kilometres with V in m/s (Γ wrong by 10³); taking the 15 km circulation as Γ/2 (linear) instead of Γ/4 (area).

## 4. How many segments?

**Result.** (a) Point vortex, Γ = 1:

| n | circle error | closed form (9) | π²/(3n²) | square error |
|---|---|---|---|---|
| 64 | 8.040 × 10⁻⁴ | 8.040 × 10⁻⁴ | 8.032 × 10⁻⁴ | 4.145 × 10⁻⁴ |
| 128 | 2.008 × 10⁻⁴ | 2.008 × 10⁻⁴ | 2.008 × 10⁻⁴ | 1.036 × 10⁻⁴ |
| 256 | 5.020 × 10⁻⁵ | 5.020 × 10⁻⁵ | 5.020 × 10⁻⁵ | 2.590 × 10⁻⁵ |
| 512 | 1.255 × 10⁻⁵ | 1.255 × 10⁻⁵ | 1.255 × 10⁻⁵ | 6.476 × 10⁻⁶ |

Rates: circle 2.001, 2.000, 2.000; square 2.000, 2.000, 2.000. The closed form matches to all printed digits; its leading term is 0.1 % low at n = 64.

(b) Vertices at angles 2πk/n on radius R. The chord from vertex k to k+1 has midpoint at angle (2k+1)π/n and radius R cos(π/n), and length 2R sin(π/n); it is perpendicular to the radius through its midpoint, hence parallel to the velocity there, whose magnitude is Γ/(2πR cos(π/n)). Product: Γ tan(π/n)/π per segment; n segments give Γ(n/π) tan(π/n). Expand tan x = x + x³/3 + …: Γ[1 + π²/(3n²) + …].

(c) Cylinder on the 4 m square: relative error 4.14 × 10⁻⁴ at n = 64, 2.59 × 10⁻⁵ at n = 256, 1.62 × 10⁻⁶ at n = 1024. Scaling the n = 64 point as n⁻²: 10⁻⁴ needs n = 64 × √4.14 = 130, i.e. n = 132 (33 per side).

(d) On a circle centred on the cylinder the doublet velocity at two diametrically opposite midpoints is the same vector (the doublet field is even under **r** → −**r**), while the two chords are opposite vectors; the pair's doublet contributions cancel exactly, leaving only the vortex term, which is (9). With the centre moved to (0.3, 0.2) m the errors are 1.608 × 10⁻², 4.017 × 10⁻³, 1.004 × 10⁻³, 2.510 × 10⁻⁴ for n = 64 to 512 at Γ = 20: identical to the centred case to four figures. The pairwise cancellation no longer holds, but the doublet has a single-valued potential V∞a²x/r², and the midpoint rule on a polygon integrates a smooth gradient field around a closed loop with an error that is exponentially small in n (the integrand is periodic and analytic), so its contribution is invisible at this precision.

**Check.** The circle column is the derivation; the square column has the same slope with a constant about half as large, consistent with the square having no chord-versus-arc error.

**Common mistakes.** Fitting the slope including n = 8 (where n//4 = 2 points per side); quoting absolute errors for Γ = −20 as if relative.

## 5. Inlet mass flow with a boundary-layer profile

**Assumptions.** Axisymmetric profile, ρ uniform across the face, no swirl.

**Derivation.** (a) ṁ = ρV · 2π∫₀^R (1 − r/R)^{1/7} r dr. With u = 1 − r/R: r = R(1 − u), dr = −R du, limits r: 0 → R become u: 1 → 0, so ṁ = 2πρVR²∫₀¹ u^{1/7}(1 − u) du = 2πρVR² [7/8 − 7/15] = 2πρVR² × 49/120 = (49/60) ρVπR². Fraction 49/60 = 0.8167; ṁ = 0.8167 × 324.76 = 265.22 kg/s.

(b) `surface_flux` on the plain disc: 266.03, 265.58, 265.38, 265.29, 265.26 kg/s for n = 20, 40, 80, 160, 320. Successive differences 0.445, 0.198, 0.088, 0.040; ratios 2.25, 2.24, 2.23; order log₂ 2.24 = 1.16. The integrand ρV(1 − r/R)^{1/7} r has an unbounded first derivative at r = R, so the midpoint error on the last radial cell is O(h^{1 + 1/7}) = O(h^{8/7}), order 1.14; the measured 1.16 is that limit approached from above.

(c) Substitute r = R(1 − w⁷): the integrand becomes ρV w · R(1 − w⁷) · 7Rw⁶, a polynomial in w. With `disc_w(w, t) = (R(1 − w⁷) cos t, R(1 − w⁷) sin t, 0)` on w ∈ [0, 1] (the normal flips because w runs inward; take the absolute value) the errors are 3.25, 0.82, 0.21, 0.052, 0.013 kg/s: rates 1.98, 1.99, 2.00, 2.00. Second order is back, with a larger constant (the polynomial has degree 13 and a large second derivative), so the substitution only wins beyond n ≈ 300.

**Result.** 49/60; 265.2 kg/s; order 1.16 → 2.00 after substitution.

**Check.** The deficit, 18 %, equals 1 − 49/60 = 11/60: for a 1/7 profile the "displacement" of mass is concentrated in the outer 30 % of the radius, where most of the area is.

**Common mistakes.** Forgetting the Jacobian r in the polar area element (gives 7/8 instead of 49/60); forgetting that the new normal points the other way.

## 6. Where the models stop

(a) True in the model for loops outside the body and its wake; the field is a gradient with a single-valued potential except at the vortex. Real sign: PIV circulation that changes by more than a few percent between a loop hugging the airfoil and one twice as large, because the loop then includes part of the viscous wake, where vorticity of both signs is present; the inviscid assumption fails. Measurable: the difference between the two loop circulations (m²/s).

(b) False even in the model unless the profile is uniform; and the velocity at the fan face is not the flight speed (the inlet diffuses the flow to about Mach 0.5 at the fan). Real sign: a pitot rake across the face showing a 1/7-type profile; the uniform-profile assumption fails by 18 %.

(c) True for any conservative field, including Earth's with J₂. Fails only if a non-conservative force acts: atmospheric drag, whose work does depend on the path (longer path at low altitude, more drag). Sign: tracked orbital energy loss per revolution at 400 km of order 10⁻⁶ of the orbital energy per day, path dependent.

(d) False in reality: the velocity in the core is finite, peaking at about Γ/(2πa) ≈ 32 m/s for the lesson's wake at a = 1.5 m. Sign: a lidar or PIV profile that is linear in r through the centre, not 1/r; the point-vortex (inviscid, zero core) assumption fails.

**Common mistakes.** Saying "viscosity" for everything; (b) is a profile and compressibility issue, (c) is drag.

## 7. The wake vortex pair

**Assumptions.** Elliptic lift distribution (b′ = πb/4), fully rolled-up vortices, point-vortex model outside the cores, no ground, no decay.

**Derivation.** (a) Γ = 4 × 686 700/(π × 1.225 × 70 × 34) = 299.9 m²/s; b′ = π × 34/4 = 26.70 m.
(b) From equation (1), the right vortex (+Γ, counter-clockwise) induces at the left vortex, where dx = −b′ and dy = 0, the velocity v = Γ dx/(2πb′²) = −Γ/(2πb′): downward. The left vortex (−Γ) induces at the right vortex, where dx = +b′, the velocity v = (−Γ)(+b′)/(2πb′²) = −Γ/(2πb′): also downward. Each vortex is carried by the other's field (a vortex induces no velocity on itself), so the pair descends at Γ/(2πb′) = 299.9/(2π × 26.70) = 1.787 m/s.
(c) Code: circulation on a 5 m circle around the left core −299.89 m²/s, around the right +299.89, around a 60 m × 20 m rectangle enclosing both 0.000 m²/s. Circulation is additive: the loop enclosing both sees −Γ + Γ = 0. The total circulation of the wake is zero, as it must be: the aircraft started from rest in still air and Kelvin's theorem (P01.12) keeps the total at zero.
(d) At the midpoint each vortex is b′/2 away and both induce −y velocity: 2 × Γ/(2π · b′/2) = 2Γ/(πb′) = 7.149 m/s downward, exactly 4 × 1.787 m/s. The factor 4 is 2 (both vortices contribute) × 2 (half the distance).
(e) Swirl at 5 m: 299.9/(2π × 5) = 9.55 m/s, vertical and of opposite sign at the two wingtips if the aircraft is centred on the core: a velocity difference of 19.1 m/s across a 10 m span. Each tip sees an angle-of-attack change of atan(9.55/40) = 13°, one up and one down; a wing near its approach angle of attack stalls on one side and the rolling moment is beyond aileron authority. A 20 m/s swirl at 40 m/s is a 27° change: not survivable at approach height. Hence the several-minute, several-nautical-mile separations behind heavy aircraft.

**Result.** Γ = 300 m²/s; b′ = 26.7 m; descent 1.79 m/s; midpoint downwash 7.15 m/s; 9.5 m/s swirl at 5 m.

**Check.** Descent rate of a few m/s and swirl of tens of m/s match published wake-vortex measurements for medium transports. Units: (N)/((kg/m³)(m/s)(m)) = m²/s. ✓

**Common mistakes.** Using b instead of b′; adding the two vortices' induced velocities at the *vortex* locations (only the other vortex counts: a vortex does not move itself); forgetting the sign of Γ on the left loop.

## Check-yourself answers

1. Circulation m²/s; 2D source strength (flow per unit depth) m²/s; 3D source strength m³/s.
2. 2Γ: the winding number is 2.
3. The potential Γθ/2π is not single-valued on a loop around the origin; it gains 2π · Γ/2π = Γ per turn. The fundamental theorem needs a single-valued potential along the whole path.
4. (1/2)² = 1/4 of Γ.
5. π²/(3 × 100²) = 3.3 × 10⁻⁴, independent of the radius (equation (9) has no R in it).
6. The sphere's parametrisation has **r**_s × **r**_t pointing outward and **r** · **n** = R > 0 everywhere; swapping the order of the parameters (or of the cross product) would flip the sign.
7. About 265 kg/s, 49/60 of 325 kg/s: the profile goes to zero at the wall, and most of the disc's area is near the wall.
8. The same, −3.28 MJ/kg: gravity is conservative, so only the endpoints matter.
