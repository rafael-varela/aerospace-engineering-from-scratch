# P10.04  Worked solutions

Each solution states: assumptions, units, the derivation, the result, at least one independent check, and the common mistakes. Seeded problems are worked for s = 15. Numbers come from `solution.py` and the probe scripts that call it, given to the figures shown.

## 1. A parabolic camber line by hand (s = 15: h = 0.025)

**Assumptions.** Thin airfoil, small angles, linearised boundary condition on the chord.

(a) dz/dx = 4h(1 − 2x), and with x = (1 − cos θ)/2, 1 − 2x = cos θ: **dz/dx = 4h cos θ**. Then A₀ = α − (1/π)∫₀^π 4h cos θ dθ = **α = 0.05236** (3°), A₁ = (2/π)∫₀^π 4h cos² θ dθ = **4h = 0.100**, and Aₙ = 0 for n ≥ 2 (orthogonality).

(b) α_L0 = −A₁/2 = −2h = −0.050 rad = **−2.865°**. c_l = π(2A₀ + A₁) = **0.6431**. c_m,c/4 = (π/4)(0 − 4h) = −πh = **−0.07854**. Ideal c_l = πA₁ = **0.3142**. x_cp = ¼(1 + π(A₁ − A₂)/c_l) = **0.3721c**.

(c) `thin_airfoil` returns A = (0.05236, 0.1000, 1 × 10⁻¹⁷, 8 × 10⁻¹⁶) and the same five numbers.

**Check.** c_l = 2π(α − α_L0) = 2π(0.05236 + 0.05) = 0.6431.

**Common mistakes.** Forgetting that A₀ contains α; quoting α_L0 with the wrong sign (positive camber gives a negative zero-lift angle).

## 2. From the vortex sheet to the Fourier coefficients

(a) With ξ = (c/2)(1 − cos φ), dξ = (c/2) sin φ dφ and x − ξ = (c/2)(cos φ − cos θ). The A₀ term: (1/2π)∫ 2VA₀(1 + cos φ) dφ/(cos φ − cos θ) = (VA₀/π)(0 + π) = VA₀ (the n = 0 Glauert integral is zero, n = 1 gives π). The Aₙ terms: sin nφ sin φ = ½[cos(n − 1)φ − cos(n + 1)φ], and the Glauert integral gives (V/π)Aₙ · π(sin(n − 1)θ − sin(n + 1)θ)/(2 sin θ) = −VAₙ cos nθ. So the left side is V(A₀ − Σ Aₙ cos nθ), and matching with V(α − dz/dx): dz/dx = (α − A₀) + Σ Aₙ cos nθ. Integrating over (0, π) gives **A₀ = α − (1/π)∫dz/dx dθ**; multiplying by cos nθ and integrating gives **Aₙ = (2/π)∫dz/dx cos nθ dθ**.

(b) Γ = ∫γ dξ = (c/2)∫γ sin θ dθ = cV(πA₀ + (π/2)A₁), so c_l = 2Γ/(Vc) = **π(2A₀ + A₁)**. The leading-edge moment −ρV∫γξ dξ gives c_m,le = −(π/2)(A₀ + A₁ − A₂/2) = −(c_l/4 + (π/4)(A₁ − A₂)). Transferring to c/4 (c_m,c/4 = c_m,le + c_l/4): **c_m,c/4 = (π/4)(A₂ − A₁)**. A₁ and A₂ depend only on the camber line, so c_m,c/4 is independent of α: the quarter chord is the aerodynamic centre.

(c) c_l = 0 ⇔ 2A₀ + A₁ = 0 ⇔ 2α − (2/π)∫dz/dx dθ + (2/π)∫dz/dx cos θ dθ = 0, so **α_L0 = −(1/π)∫dz/dx (cos θ − 1) dθ**.

(d) dz/dx = −δ for θ > θ_h. A₀ = α + δ(π − θ_h)/π; A₁ = −(2δ/π)∫cos θ dθ = (2δ/π) sin θ_h; A₂ = (δ/π) sin 2θ_h. Then c_l = 2πα + 2δ(π − θ_h + sin θ_h) = 2π(α + τδ) with **τ = 1 − (θ_h − sin θ_h)/π**, and c_m,c/4 = (δ/4)(sin 2θ_h − 2 sin θ_h) = **−½δ sin θ_h(1 − cos θ_h)**. For a small tab, cos θ_h = 1 − 2(1 − E) = −1 + 2E, so θ_h ≈ π − 2√E, sin θ_h ≈ 2√E, and **τ ≈ 4√E/π**: a 1 % tab still has τ = 0.13.

**Common mistakes.** Losing the 1/2 in sin nφ sin φ; transferring the moment the wrong way (c_m,c/4 = c_m,le + c_l/4, because the lift acts aft of the leading edge).

## 3. A reflexed section for a tailless glider

**Assumptions.** Thin-airfoil superposition of the 2412 camber line and a plain tab (a camber-line slope jump at x = 0.85).

(a) The 2412 line has c_m,c/4 = −0.05312; the tab's dc_m/dδ = −0.6070 /rad for E = 0.15. Zero moment needs **δ = −5.01°** (trailing edge up). The combined line's c_m,c/4 is 4 × 10⁻¹⁷.

(b) τ = 0.4805, so α_L0 rises from −2.077° to **+0.332°**. c_l = 0.5 needs **α = 4.89°** instead of 2.48°.

(c) The price is lift: 2.4° more incidence for the same c_l, or equivalently the loss of most of the camber's free lift, plus the drag of the deflected tab. Designing the reflex into a smooth camber line (an S-shaped mean line) gets the same zero moment with less local curvature at the hinge and so a better-behaved boundary layer; a tab is the quick, adjustable version (and an elevon, which is how tailless aircraft trim).

**Check.** Superposition: c_m = −0.05312 + (−0.6070)(−0.08751) = 0.0000.

**Common mistakes.** Deflecting the tab down (makes c_m more negative); forgetting that τ < 1 for a small flap when estimating the lift penalty.

## 4. Camber position

(a) α_L0 (deg) and c_m,c/4:

| m \ p | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 |
|---|---|---|---|---|---|
| 0.02 | −1.799, −0.03696 | −1.918, −0.04473 | −2.077, −0.05312 | −2.292, −0.06283 | −2.592, −0.07489 |
| 0.04 | −3.598, −0.07392 | −3.836, −0.08946 | −4.154, −0.10624 | −4.584, −0.12566 | −5.184, −0.14979 |
| 0.06 | −5.396, −0.11088 | −5.754, −0.13419 | −6.232, −0.15936 | −6.875, −0.18850 | −7.776, −0.22468 |

(b) Each column doubles and triples with m exactly: the coefficients are linear in the slope. Moving the camber aft makes both larger in magnitude: aft camber is more effective at producing lift at zero α (the trailing-edge region dominates the (cos θ − 1) weight) and produces a stronger nose-down moment.

(c) The 230 line gives α_L0 = **−1.094°** and c_m,c/4 = **−0.01284**, against −2.077° and −0.05312 for the 2412. Forward camber gives a quarter of the pitching moment, which means a smaller tail download to trim and less trim drag; it was also found to give a higher maximum lift. The cost is less lift at zero angle.

## 5. The discrete-vortex check

(a) Error of the lumped-vortex c_l, NACA 2412 line, α = 4°:

| n | uniform | cosine |
|---|---|---|
| 10 | 1.461 × 10⁻⁴ | 8.04 × 10⁻⁶ |
| 20 | 3.833 × 10⁻⁵ | 2.697 × 10⁻⁴ |
| 40 | 9.811 × 10⁻⁶ | 6.89 × 10⁻⁶ |
| 80 | 2.481 × 10⁻⁶ | 1.052 × 10⁻⁵ |
| 160 | 6.239 × 10⁻⁷ | 4.79 × 10⁻⁶ |
| 320 | 1.564 × 10⁻⁷ | 2.18 × 10⁻⁷ |

Uniform orders 1.93, 1.97, 1.98, 1.99, 2.00; cosine erratic (−5.1 to 5.3).

(b) The camber slope has a kink in its derivative at x = 0.4. On a uniform grid with n a multiple of 5 that point is always a panel edge, so every panel sees a smooth slope and the error expansion is clean. On cosine spacing the kink falls at a different place inside a panel for each n, and the local error depends on where, so the errors jump about while still trending down. The lesson: put a node at every kink.

(c) For z = 0.12x(1 − x) (h = 0.03): c_l is exact to 10⁻¹⁵ at n = 10, 40 and 160 on cosine spacing (and on uniform). c_m errors on uniform spacing: **9.42 × 10⁻⁴, 2.36 × 10⁻⁴, 5.89 × 10⁻⁵, 1.47 × 10⁻⁵**, order **2.0**. The quarter-chord vortex placement makes the total circulation exact for a parabola but not its distribution.

## 6. How far can you trust a thin-airfoil flap?

(a) A 30° plain flap breaks (1) **attached flow**: the boundary layer separates at the hinge, where the surface turns sharply, which is the first to fail and caps the lift increment far below 2.0; (2) **small slopes**: tan 30° = 0.58 is not 0.52, and the flap's projected chord shrinks; (3) the **planar wake and linearised boundary condition** on the chord. Theory says Δc_l is linear for ever; real plain flaps lose effectiveness as the deflection grows (`P10.08` uses NACA TR-938's data to quantify it).

(b) Thickness raises the lift slope (+11 % at 12 % thickness by the panel method), but a thick section's flap sits in a thicker boundary layer near the trailing edge, which reduces its effectiveness; the two effects partly cancel, and only the data say which wins. The nose-down moment also grows (the panel method's c_m,c/4 is −0.0612 against −0.0531 at 4° without a flap).

(c) No. 6.7° is the inviscid answer for a 2-D section; the finite wing needs more lift coefficient for the same aircraft lift and loses some to downwash (`P10.09`), and the real flap is less effective. Design the flap for a larger deflection range (20–30°) with a slotted flap that keeps the flow attached (`P10.08`), and check the approach attitude and stall margin, not just the c_l.

## Check-yourself answers

1. x = (c/2)(1 − cos θ); it maps the chord to (0, π), clusters near both edges, and turns the singular integral equation into Fourier series.
2. c_l = π(2A₀ + A₁) = 2π(α − α_L0); the lift slope is 2π per radian.
3. (π/4)(A₂ − A₁): it depends only on the camber line, so it does not change with α; the quarter chord is the aerodynamic centre.
4. α_L0 = −(1/π)∫₀^π (dz/dx)(cos θ − 1) dθ; −2.077° for the NACA 2412 camber line.
5. A plain flap is a camber line with a slope jump of −δ at the hinge; its effectiveness τ = 1 − (θ_h − sin θ_h)/π is 0.609 for a 25 % flap.
6. A discrete-vortex (lumped-vortex) model with vortices at the quarter points and collocation at the three-quarter points of each panel, converging at second order; and the vortex panel method of `P10.05` as the thickness goes to zero.
7. Because the camber slope (or the flap's slope jump) has a kink there: splitting the integral makes each piece smooth so Gauss–Legendre converges fast (to round-off with 8 points per piece, against 2.5 × 10⁻⁶ with 64 unsplit points).
