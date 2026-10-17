# P10.02  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. One station by hand (s = 15: NACA 4416)

**Assumptions.** Standard (open) trailing edge; m = 0.04, p = 0.4, t = 0.16; chord 1.

(a) y_t = 5(0.16)(0.2969√0.25 − 0.1260(0.25) − 0.3516(0.0625) + 0.2843(0.015625) − 0.1015(0.00390625)) = 0.8 × 0.099021 = **0.079217**. On the forward parabola, y_c = (0.04/0.16)(2 × 0.4 × 0.25 − 0.0625) = **0.034375** and dy_c/dx = (0.08/0.16)(0.4 − 0.25) = **0.075**, so θ = **4.289°**. The surface points:

x_u = 0.25 − 0.079217 sin θ = **0.244075**, y_u = 0.034375 + 0.079217 cos θ = **0.113370**;
x_l = 0.25 + 0.079217 sin θ = **0.255925**, y_l = 0.034375 − 0.079217 cos θ = **−0.044620**.

The upper point sits 5.9 × 10⁻³ ahead of the station and the lower point the same distance behind: thickness laid perpendicular to a sloping camber line tilts the "rib".

(b) r = 1.1019 × 0.16² = **0.028208 c**; trailing-edge thickness 2y_t(1) = 0.021 × 0.16 = **0.00336 c**.

(c) `naca4_thickness(0.25, 0.16)` and `naca4_camber(0.25, 0.04, 0.4)` return the same values to round-off (the hand values above were copied from that call to six figures).

**Check.** The polygon area of the generated 4416 (n = 400) is 0.109995, against the thickness form's 0.68508 × 0.16 = 0.109613: camber adds 0.35 % here.

**Common mistakes.** Adding y_t vertically; using the aft parabola at x < p; forgetting that x_u ≠ x.

## 2. The constants behind the families

(a) ∫₀¹ 2y_t dx = 10t(2a₀/3 + a₁/2 + a₂/3 + a₃/4 + a₄/5) = 10t(0.197933 − 0.063 − 0.1172 + 0.071075 − 0.0203) = **0.68508 t**; with a₄ = −0.1036, **0.68088 t**. (`thickness_area_coefficient` returns these, and a test checks them against adaptive quadrature to 10⁻⁸.) Camber adds area only through the perpendicular construction, a second-order effect: the 2412's polygon area is 0.6857 t c², **0.09 %** above the 0012's.

(b) As x → 0 the √x term dominates: y_t ≈ 5t a₀ √x, so y_t² ≈ 25t²a₀² x = 2rx with **r = (5 × 0.2969 t)²/2 = 1.1019 t²**. A parabola y² = 2rx has radius of curvature r at its vertex.

(c) 0.2969 − 0.1260 − 0.3516 + 0.2843 − 0.1036 = **0.0000**, so y_t(1) = 0. With −0.1015 the sum is 0.0021 and the trailing edge is **0.021 t** thick (0.00252c for 12 %).

(d) dy_c/dx = (k₁/6)(3x² − 6rx + r²(3 − r)). Setting it to zero at x = p: 3p² − 6rp + 3r² − r³ = 0, i.e. **r³ − 3r² + 6pr − 3p² = 0**. For p = 0.15 the root between p and 1 is **r = 0.2027** (bisection; `naca5_r`). Every term of y_c carries the factor k₁, so the slope, A₁ and the ideal lift coefficient πA₁ are all proportional to k₁: k₁ = c_l,design/(πA₁ evaluated with k₁ = 1). With the quadrature split at x = r, this gives k₁ = 15.920 for the 230 line.

**Common mistakes.** Integrating √x as x^{1/2}/(1/2) (that is the derivative's reciprocal, not the integral); treating 0.685 as exact.

## 3. A battery pocket in a foam UAV wing

**Assumptions.** Rectangular wing, constant section; foam 30 kg/m³; vertical thickness for the battery (a flat slab).

(a)

| section | area (m²) | foam volume (L) | foam mass (g) |
|---|---|---|---|
| NACA 2412 | 0.0032913 | 4.937 | 148.1 |
| NACA 4412 | 0.0032999 | 4.950 | 148.5 |

(b) The 2412's vertical thickness is at least 20 mm from **x/c = 0.121 to 0.543** (24.0 mm, its maximum, at 30 %); the 4412's from 0.118 to 0.543. Camber barely changes it: thickness, not camber, sets the pocket.

(c) Between 20 % and 60 % chord the section holds 0.0017897 m² (2412), so the pocket region holds **2.68 L** of the wing's 4.94 L, **54 %** of the foam. A 20 mm slab fits only from 0.121c to 0.543c, so a full-depth pocket over 20–60 % would break through the skin aft of 54 %: shorten it or accept a shallower aft end.

**Check.** 0.0032913 m² / (0.12 × 0.2²) = 0.6857, the 2412 coefficient from main().

**Common mistakes.** Using the thickness perpendicular to the camber line for a flat slab; forgetting the c² scaling.

## 4. Does spacing matter?

(a) NACA 2412, error against n = 20 000:

| n | cosine | uniform |
|---|---|---|
| 10 | 1.333 × 10⁻³ | 2.460 × 10⁻³ |
| 20 | 3.348 × 10⁻⁴ | 8.590 × 10⁻⁴ |
| 40 | 8.380 × 10⁻⁵ | 3.009 × 10⁻⁴ |
| 80 | 2.095 × 10⁻⁵ | 1.057 × 10⁻⁴ |
| 160 | 5.235 × 10⁻⁶ | 3.719 × 10⁻⁵ |
| 320 | 1.305 × 10⁻⁶ | 1.310 × 10⁻⁵ |

Observed orders: cosine **1.99, 2.00, 2.00, 2.00, 2.00**; uniform **1.52, 1.51, 1.51, 1.51, 1.51**. (The reference itself changes by 5 × 10⁻⁹ between 5000 and 20 000 panels.)

(b) Away from the nose the polygon cuts each chord of a smooth curve and loses area of order h³ per cell, h² overall. In the first cell, y ≈ k√x is not smooth: the chord from (0, 0) to (h, k√h) misses ∫₀ʰ k√x dx − kh^{3/2}/2 = kh^{3/2}/6, an error of order h^{1.5}, which dominates. With x = (1 − cos β)/2 ≈ β²/4, √x ≈ β/2: in the variable β the surface is smooth, so uniform steps in β give second order everywhere.

(c) Uniform, n = 100: the first point is at x = 0.0100, **0.63 r** (r = 0.01587 for t = 0.12): the nose radius is spanned by one panel. Cosine: x = (1 − cos(π/100))/2 = 2.47 × 10⁻⁴, **0.016 r**: dozens of panels around the nose, where `P10.05`'s suction peak lives.

**Common mistakes.** Comparing against a reference no finer than the test grids; attributing order 1.5 to the camber kink (the 0012 shows the same 1.5).

## 5. The 5-digit constants

(a)

| line | derived r | derived k₁ | quoted (r, k₁) | ideal c_l of quoted pair |
|---|---|---|---|---|
| 210 | 0.0581 | 350.33 | (0.0580, 361.4) | 0.3084 |
| 220 | 0.1257 | 51.577 | (0.1260, 51.64) | 0.3019 |
| 230 | 0.2027 | 15.920 | (0.2025, 15.957) | 0.3000 |
| 240 | 0.2903 | 6.624 | (0.2900, 6.643) | 0.3001 |
| 250 | 0.3913 | 3.223 | (0.3910, 3.230) | 0.3000 |

(b) The **210** pair is furthest: its ideal c_l is **0.3084, 2.8 % high**. The discrepancy is in k₁: with the quoted r = 0.0580, the design condition needs k₁ = 351.6, not 361.4; the 0.0001 difference in r accounts for only 0.35 % of k₁. The 230 pair, the one in common use (23012, 23015), satisfies its design condition to 0.01 %. The lesson: verify a table by its defining condition rather than copying it.

(c) "430": **r = 0.2027** (unchanged: r depends only on where the maximum is, which does not depend on the scale factor) and **k₁ = 31.839**, twice the 230 value.

**Common mistakes.** Integrating the ideal c_l across the kink at x = r with one unsplit Gauss rule (a few × 10⁻⁵ error, larger than the 230 discrepancy you are trying to detect).

## 6. What does "12 % thick" mean?

(a) Vertical thickness and mean line of the generated surfaces (n = 600):

| section | max vertical t | at x | max mean line | at x |
|---|---|---|---|---|
| 0012 | 0.12003 | 0.300 | 0 | – |
| 2412 | 0.12007 | 0.300 | 0.02000 | 0.402 |
| 6412 | 0.12039 | 0.292 | 0.06000 | 0.402 |
| 9415 | 0.15110 | 0.284 | 0.09000 | 0.403 |

The designation gives thickness perpendicular to the camber line; the vertical thickness differs by up to 0.7 % for the 9415, and its position moves forward with camber.

(b) The 23012's mean line peaks at **x = 0.143**, not 0.150: the forward camber line is steep, and the perpendicular construction shifts the upper surface forward and the lower one aft by different amounts, so the vertical average of the two is not the camber line.

(c) The difference between the vertical and the perpendicular construction is largest where θ is largest, at the nose: up to 0.0018c on the 2412's upper surface and 0.0066c on the 6412's, against 0.0002c and 0.0011c over mid-chord. A panel method notices, because the suction peak and the leading-edge radius live there; the lift barely changes but the peak Cp, and with it the boundary layer and stall of `P10.06` and `P10.15`, do. For model-making, always ask which construction the coordinates used.

**Common mistakes.** Reading the mean line of tabulated coordinates as the camber line; assuming the maximum thickness of a cambered section is exactly t at exactly 30 %.

## Check-yourself answers

1. 2 % camber (m = 0.02) at 40 % chord (p = 0.4), 12 % thick (t = 0.12).
2. Half-cosine spacing x = (1 − cos β)/2 with uniform β: it clusters points at the leading edge, where the surface curves fastest and y ∝ √x, and at the trailing edge, where the Kutta condition acts.
3. From the trailing edge along the upper surface to the leading edge, then along the lower surface back to the trailing edge: counter-clockwise, so the shoelace area is positive.
4. Area = 0.68508 t c² (open trailing edge), from integrating the polynomial term by term; 0.68088 with the closed edge.
5. 1.1019 t² c, from y_t ≈ 5 × 0.2969 t √x near the nose: a parabola y² = 2rx.
6. The camber line is linear in k₁, so its ideal lift coefficient is too; k₁ is chosen to make that ideal c_l equal 0.15 times the first digit (0.3 for the 230 line). r is fixed first, by the maximum-camber position.
7. Polygon-area error falls at order 2 with cosine spacing and only 1.5 with uniform spacing, because of the √x nose.
